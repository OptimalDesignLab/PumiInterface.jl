# file for function used for node/dof numbering


# this can be generalized to 3D trivially
function populateDofNumbers(mesh::PumiMesh)
# populate the dofnums_Nptr with values calculated from
# nodenums_Nptr

#  resetAllIts2(mesh.m_ptr)
  # mesh iterator increment, retreval functions
  # TODO: checks if using mutable arrays is causing dynamic dispatch
#  iterators_inc = [incrementVertIt, incrementEdgeIt, incrementFaceIt, incrementElIt]
#  iterators_get = [getVert, getEdge, getFace, getEl]
  num_entities = mesh.numEntitiesPerType
  num_nodes_entity = mesh.numNodesPerType

  for etype = 1:(mesh.dim+1) # loop over entity types
    if (num_nodes_entity[etype] != 0)  # if no nodes on this type of entity, skip
      it = MeshIterator(mesh.m_ptr, etype - 1)
      for entity = 1:num_entities[etype]  # loop over all entities of this type
#	entity_ptr = iterators_get[etype]()  # get entity
        entity_ptr = iterate(mesh.m_ptr, it)

	for node = 1:num_nodes_entity[etype]
          nodenum = getNumberJ(mesh.nodenums_Nptr, entity_ptr, node-1, 0)
	  if nodenum != 0
	    for i=1:mesh.numDofPerNode
	      dofnum_i = (nodenum -1)*mesh.numDofPerNode + i
  	      numberJ(mesh.dofnums_Nptr, entity_ptr, node-1, i-1, dofnum_i)
	    end  # end loop over dofsPerNode
	  end   # end if nodenum != 0
	end  # end loop over node

#      iterators_inc[etype]()
      end  # end loop over entitiesa
      free(mesh.m_ptr, it)
    end  # end if 
  end  # end loop over entity types

#  resetAllIts2(mesh.m_ptr)

#  writeVtkFiles("dofs_numbered", mesh.m_ptr)
  return nothing
end


# this can be generalized trivially
function populateNodeStatus(mesh::PumiMesh)
# populate the nodestatus_Nptr with values
# currently we set all nodes to status 3 (free)

#  resetAllIts2(mesh.m_ptr)
  # mesh iterator increment, retreval functions
  # TODO: checks if using mutable arrays is causing dynamic dispatch
#  iterators_inc = [incrementVertIt, incrementEdgeIt, incrementFaceIt, incrementElIt]
#  iterators_get = [getVert, getEdge, getFace, getEl]
  num_entities = mesh.numEntitiesPerType
  num_nodes_entity = mesh.numNodesPerType

#  mesh.numNodesPerType = num_nodes_entity

  for etype = 1:(mesh.dim + 1) # loop over entity types
    if (num_nodes_entity[etype] != 0)  # if no nodes on this type of entity, skip
      it = MeshIterator(mesh.m_ptr, etype - 1)
      for entity = 1:num_entities[etype]  # loop over all entities of this type
#	entity_ptr = iterators_get[etype]()  # get entity
        entity_ptr = iterate(mesh.m_ptr, it)

	for node = 1:num_nodes_entity[etype]
	  numberJ(mesh.nodestatus_Nptr, entity_ptr, node-1, 0, 3)
	end  # end loop over node

#      iterators_inc[etype]()
      end  # end loop over entitiesa
      free(mesh.m_ptr, it)
    end  # end if 
  end  # end loop over entity types

#  resetAllIts2(mesh.m_ptr)
  return nothing
end




# this can be generalized to 3D trivially
function getDofNumbers(mesh::PumiMeshDG)
# populate array of dof numbers, in same shape as solution array u (or q)

#println("in getDofNumbers")
#println("numNodesPerElement = ", mesh.numNodesPerElement)

mesh.dofs = Array(Int, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.shared_element_offsets[end] - 1)

for i=1:mesh.numEl
  #TODO: use the non-allocating version
  dofnums = getGlobalNodeNumbers(mesh, i)

  for j=1:mesh.numNodesPerElement
    for k=1:mesh.numDofPerNode  # loop over dofs on the node
      mesh.dofs[k, j, i] = dofnums[k,j]
    end
  end
end
# get dof number of non-local elements
# post receives first
send_reqs = mesh.send_reqs
recv_reqs = mesh.recv_reqs
for i=1:mesh.npeers
  # get the pointer of the start location
  start_el = mesh.shared_element_offsets[i]
  numel = mesh.shared_element_offsets[i+1] - start_el
  ndata = mesh.numDofPerNode*mesh.numNodesPerElement*numel
  lin_idx = sub2ind(size(mesh.dofs), 1, 1, start_el)
  ptr_start_el = pointer(mesh.dofs, lin_idx)
  peer_i = mesh.peer_parts[i]
  recv_reqs[i] = MPI.Irecv!(ptr_start_el, ndata, peer_i, 1, mesh.comm)
end

# now send data
dof_sendbuf = Array(Array{Int, 3}, mesh.npeers)
for i=1:mesh.npeers
  elnums_i = mesh.local_element_lists[i]
  numel = length(elnums_i)
  dof_sendbuf[i] = Array(Int, mesh.numDofPerNode, mesh.numNodesPerElement, numel)
  sendbuf_i = dof_sendbuf[i]
  for j=1:length(elnums_i)
    elnum_j = elnums_i[j]
    for k=1:mesh.numNodesPerElement
      for p=1:mesh.numDofPerNode
        sendbuf_i[p, k, j] = mesh.dofs[p, k, elnum_j]
      end
    end
  end  # end loop over current list of elements

  # now do the send
  send_reqs[i] = MPI.Isend(sendbuf_i, mesh.peer_parts[i], 1, mesh.comm)
end

# figure out the local to global offset of the dof numbers
# compute number of local dofs
ndof = 0
for i=1:(mesh.dim + 1)
  ndof += mesh.numNodesPerType[i]*mesh.numEntitiesPerType[i]*mesh.numDofPerNode
end
println(mesh.f, "ndof = ", ndof)
# get all processes dof offsets
dof_offsets = MPI.Allgather(ndof, mesh.comm)
println(mesh.f, "dof_offsets = ", dof_offsets)
# figure out own dof_offset
dof_offset = 0
for i=1:mesh.myrank  # sum all dofs up to but *not* including self
  dof_offset += dof_offsets[i]
end
mesh.dof_offset = dof_offset

# figure out peer dof_offsets
peer_dof_offsets = Array(Int, mesh.npeers)
for i=1:mesh.npeers
  dof_offset_i = 0
  for j=1:mesh.peer_parts[i]
    dof_offset_i += dof_offsets[j]
  end
  peer_dof_offsets[i] = dof_offset_i
end

println(mesh.f, "peer_parts = ", mesh.peer_parts)
println(mesh.f, "peer_dof_offsets = ", peer_dof_offsets)
# adjust other elements dof such that dof + dof_offset of this process = global
# dof number
# store local dof number + (offset of owner - offset of host) in dofs array
dofs = mesh.dofs
for i=1:mesh.npeers
  idx, stat =  MPI.Waitany!(recv_reqs)
  # do subtraction first to avoid overflow
  offset = peer_dof_offsets[idx] - dof_offset
  start_el = mesh.shared_element_offsets[idx]
  end_el = mesh.shared_element_offsets[idx+1] - 1
  for j=start_el:end_el
    for k=1:mesh.numNodesPerElement
      for p=1:mesh.numDofPerNode
        dofs[p, k, j] += offset
      end
    end
  end

end  # end loop over peers

# wait for the sends to finish before returning
for i=1:mesh.npeers
  MPI.Wait!(send_reqs[i])
end

#println("mesh.dof = ", mesh.dofs)

return nothing

end



# this could be generalized with some clever array bookkeeping tricks
function numberNodes(mesh::PumiMesh, number_dofs=false)
# assign node numbers to entire mesh, or dof numbers if number_dofs=true,
# using the correct Numbering pointer
# assumes mesh elements have already been reordered
# this works for high order elements

  numberNodesElement(mesh, number_dofs=number_dofs, start_at_one=false)

  # calculate number of nodes, dofs
  numnodes = mesh.numNodes

  # initally number all dofs as numDof+1 to 2*numDof
  # this allows quick check to see if somthing is labelled or not

#  resetAllIts2(mesh.m_ptr)
  # mesh iterator increment, retreval functions
#  iterators_inc = [incrementVertIt, incrementEdgeIt, incrementFaceIt, incrementElIt]
#  iterators_get = [getVert, getEdge, getFace, getEl]
  num_entities = mesh.numEntitiesPerType
  num_nodes_entity = mesh.numNodesPerType

#  mesh.numNodesPerType = num_nodes_entity

  if number_dofs
    numbering_ptr = mesh.dofnums_Nptr
    curr_dof = mesh.numDof + 1
    dofpernode = mesh.numDofPerNode
    numDof = mesh.numDofPerNode*numnodes
  else  # do node numbering
    numbering_ptr = mesh.nodenums_Nptr
    curr_dof = mesh.numNodes + 1
    dofpernode = 1
    numDof = numnodes
  end


  verts_i = Array(Ptr{Void}, 12)
  edges_i = Array(Ptr{Void}, 12)
#  resetAllIts2(mesh.m_ptr)
  el_i_ptr = Ptr{Void}(0)  # hold current element
# TODO: move all if statements out one for loop (check only first dof on each node)
  curr_dof = 1
  it = MeshIterator(mesh.m_ptr, mesh.dim)
  for i=1:mesh.numEl
    el_i_ptr = iterate(mesh.m_ptr, it)
#    el_i_ptr = getFace()
#    incrementFaceIt()
    # get vertices, edges for this element
    numVert = getDownward(mesh.m_ptr, el_i_ptr, 0, verts_i)
#    println("verts_i = ", verts_i)
#    println("mesh.verts = ", mesh.verts)
    numEdge = getDownward(mesh.m_ptr, el_i_ptr, 1, edges_i)
    for j=1:3  # loop over vertices, edges
      vert_j = verts_i[j]
      edge_j = edges_i[j]
      for k=1:num_nodes_entity[1] # loop over vertex nodes
        for p=1:dofpernode
          dofnum_k = getNumberJ(numbering_ptr, vert_j, k-1, p-1)
          if dofnum_k > numDof  # still has initial number
            # give it new (final) number
            numberJ(numbering_ptr, vert_j, 0, k-1, curr_dof)
            curr_dof += 1
          end
        end
      end  # end loop over vertex dofs
      
      # loop over nodes on edge
      for k=1:num_nodes_entity[2]  # loop over nodes
	for p=1:dofpernode  # loop over dofs
	  dofnum_p = getNumberJ(numbering_ptr, edge_j, k-1, p-1)
	  if dofnum_p > numDof  # still has initial number
	    # give it new (final) number
	    numberJ(numbering_ptr, edge_j, k-1, p-1, curr_dof)
	    curr_dof += 1
	  end
	end
      end
    end  # end loop over vertices, edges
    # label face nodes
    for k=1:num_nodes_entity[3]  # loop over nodes on face
        # visit the nodes in the SBP ordering, so consider k to be the SBP
        # node number
        # this only works when there are no nodes on edges or vertices, because
        # we would have to do additional translation to take them into account

      for p=1:dofpernode  # loop over dofs
	dofnum_p = getNumberJ(numbering_ptr, el_i_ptr, k-1, p-1)
	if dofnum_p > numDof
	  numberJ(numbering_ptr, el_i_ptr, k-1, p-1, curr_dof)
	  curr_dof += 1
	end
      end
    end  # end loop over face nodes
  end  # end loop over elements

  free(mesh.m_ptr, it)
#  resetAllIts2(mesh.m_ptr)


  if (curr_dof -1) != numDof 
    println("Warning: number of dofs assigned is not equal to teh expected number")
    println("number of dofs assigned = ", curr_dof-1, " ,expected number = ", numDof)
  end

  return nothing

end

function numberNodesElement(mesh::PumiMesh; number_dofs=false, start_at_one=true)
# number all nodes (or dofs) by iterating over all entities, from lowest
# dimension to highest, and numbering the nodes in order

  # calculate number of nodes, dofs
  numnodes = mesh.numNodes

  # initally number all dofs as numDof+1 to 2*numDof
  # this allows quick check to see if somthing is labelled or not

#  resetAllIts2(mesh.m_ptr)
  # mesh iterator increment, retreval functions
#  iterators_inc = [incrementVertIt, incrementEdgeIt, incrementFaceIt, incrementElIt]
#  iterators_get = [getVert, getEdge, getFace, getEl]
  num_entities = mesh.numEntitiesPerType
  num_nodes_entity = mesh.numNodesPerType

#  mesh.numNodesPerType = num_nodes_entity

  if number_dofs
    numbering_ptr = mesh.dofnums_Nptr
    curr_dof = mesh.numDof + 1
    dofpernode = mesh.numDofPerNode
    numDof = mesh.numDofPerNode*numnodes
  else  # do node numbering
    numbering_ptr = mesh.nodenums_Nptr
    curr_dof = mesh.numNodes + 1
    dofpernode = 1
    numDof = numnodes
  end

  if (numDof > 2^30 || numDof < 0)
    println("Warning: too many dofs, renumbering will fail")
  end


#  println("num_entities = ", num_entities)
#  println("num_nodes_entity = ", num_nodes_entity)

  if start_at_one
    curr_dof = 1
  else
    curr_dof = numDof+1
  end

  if mesh.isDG
    nodemap =  mesh.nodemapSbpToPumi
  else
    nodemap = collect(eltype(mesh.nodemapSbpToPumi), 1:mesh.numNodesPerElement)
  end

  for etype = 1:(mesh.dim + 1) # loop over entity types
    if (num_nodes_entity[etype] != 0)  # if no nodes on this type of entity, skip
      it = MeshIterator(mesh.m_ptr, etype - 1)
      for entity = 1:num_entities[etype]  # loop over all entities of this type
#	entity_ptr = iterators_get[etype]()  # get entity
        entity_ptr = iterate(mesh.m_ptr, it)

	for node = 1:num_nodes_entity[etype]
          pumi_node = nodemap[node]
	  for dof = 1:dofpernode
	    numberJ(numbering_ptr, entity_ptr, pumi_node-1, dof-1, curr_dof)
	    curr_dof += 1
	  end  # end loop over dof
	end  # end loop over node

#      iterators_inc[etype]()
      end  # end loop over entitiesa

      free(mesh.m_ptr, it)
    end  # end if 
  end  # end loop over entity types

  return nothing
end


function getDofNumbers(mesh::PumiMesh2)
# populate array of dof numbers, in same shape as solution array u (or q)

mesh.dofs = Array(Int32, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)

for i=1:mesh.numEl
  dofnums = getGlobalNodeNumbers(mesh, i)

  for j=1:mesh.numNodesPerElement
    for k=1:mesh.numDofPerNode  # loop over dofs on the node
      mesh.dofs[k, j, i] = dofnums[k,j]
    end
  end
end

#println("mesh.dof = ", mesh.dofs)

return nothing

end


function numberNodesWindy(mesh::PumiMeshDG, start_coords, number_dofs=false)
# number nodes and elements starting from a particular location

  println("numbering nodes windy")
  println("startin from ", start_coords)

  # calculate number of nodes, dofs
  numDof = mesh.numDof
  if (numDof > 2^30 || numDof < 0)
    println("Warning: too many dofs, renumbering will fail")
  end


  if number_dofs
    numbering_ptr = mesh.dofnums_Nptr
    dofpernode = mesh.numDofPerNode
    numDof = mesh.numDofPerNode*numnodes
  else  # do node numbering
    numbering_ptr = mesh.nodenums_Nptr
    dofpernode = 1
    numDof = mesh.numNodes
  end
  dim = mesh.dim
  el_Nptr = mesh.el_Nptr
  adj_els = Array(Ptr{Void}, mesh.numFacesPerElement)
  nodemap = mesh.nodemapSbpToPumi
  
  # initially number all components in range (numEl+1):(2*numEl)
  curr_elnum = mesh.numEl+1
#  resetIt(dim)
  it = MeshIterator(mesh.m_ptr, dim)
  for i=1:mesh.numEl
    el_i = iterate(mesh.m_ptr, it)
#    el_i = getEntity(dim)

    numberJ(el_Nptr, el_i, 0, 0, curr_elnum)
    curr_elnum += 1
#    incrementIt(dim)
  end
  free(mesh.m_ptr, it)
  @assert (curr_elnum - 1) == 2*mesh.numEl

  # do the final numbering
  start_el = getStartEl(mesh, start_coords)
  que = FIFOQueue{Ptr{Void}}(size_hint = div(mesh.numEl, 2))
  push!(que, start_el)

  curr_elnum = 1
  curr_dof = 1
  # the dreaded while loop
  while (!isempty(que))
i#    print("\n")
    curr_el = pop!(que)

    # only unlabelled entities are added to the que, and they are added 
    # exactly once, so no need to check here
    numberJ(el_Nptr, curr_el, 0, 0, curr_elnum)
    curr_elnum += 1

    # now label dofs
    for i=1:mesh.numNodesPerElement
      pumi_node = nodemap[i]
      for j=1:dofpernode
        numberJ(numbering_ptr, curr_el, pumi_node-1, j-1, curr_dof)
        curr_dof += 1
      end
    end

    # add face adjacent neighbor elements to the que
    numadj = countBridgeAdjacent(mesh.m_ptr, curr_el, dim-1, dim)
    getBridgeAdjacent(adj_els)

    for i=1:numadj
      el_i = adj_els[i]
      elnum = getNumberJ(el_Nptr, el_i, 0, 0)

       # if not numbered and not in que
       # when adding to que, scale elnum so it is in range (2*numEl+1):3*numEl
       # as long as there are at least 3 nodes per element, node numbers will
       # overflow before this
      if (elnum > mesh.numEl) && (elnum <= 2*mesh.numEl)
        numberJ(el_Nptr, el_i, 0, 0, elnum + mesh.numEl)
        push!(que, el_i)
      end
    end  # end loop over adjacent elements

  end  # end while loop 

  if (curr_dof -1) != numDof 
    println("Warning: number of dofs assigned is not equal to the expected number")
    println("number of dofs assigned = ", curr_dof-1, " ,expected number = ", numDof)
    throw(ErrorException("dof numbering failed"))
  end

  if (curr_elnum - 1) != mesh.numEl
    throw(ErrorException("element numbering failed"))
  end

  # the element number should have been zero based (oops), so decrement
  # the element numbers
  for i=1:mesh.numEl
    el_i = mesh.elements[i]
    idx = getNumberJ(mesh.el_Nptr, el_i, 0, 0)
    numberJ(mesh.el_Nptr, el_i, 0, 0, idx - 1)
  end

  # need to update the element pointer array because element numbering has
  # changed
  #TODO: don't reallocate the arrays
  mesh.verts, mesh.edges, mesh.faces, mesh.elements = getEntityPointers(mesh)
  return nothing
end


function getStartEl(mesh::PumiMeshDG, start_coords)
# find the element with centroid closest to start_coords

  numVertsPerElement = mesh.numTypePerElement[1]
  coords = zeros(3, numVertsPerElement)
  centroid = zeros(3)
  dim = mesh.dim  # hoist the lookup

#  resetIt(dim)

  min_norm = typemin(Float64)
  min_el = Ptr{Void}(0)

  it = MeshIterator(mesh.m_ptr, dim)
  for i=1:mesh.numEl
    el_i = iterate(mesh.m_ptr, it)
#    el_i = getEntity(dim)
    getElementCoords(mesh, el_i, coords)

    # compute centroid
    for j=1:dim
      for k=1:numVertsPerElement
        centroid[j] += coords[j, k]
      end
      centroid[j] = centroid[j]/numVertsPerElement
    end
    
    # compute norm of distance from start_coords
    nrm = 0.0
    for j=1:dim
      nrm += (centroid[j] - start_coords[j])*(centroid[j] - start_coords[j])
    end
    # skip the square root because it doesn't affect the comparison

    if nrm > min_norm
      min_norm = nrm
      min_el = el_i
    end

    fill!(centroid, 0.0)
#    incrementIt(dim)
  end

  return min_el
end


