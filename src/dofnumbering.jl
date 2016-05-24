# file for function used for node/dof numbering



function populateDofNumbers(mesh::PumiMesh)
# populate the dofnums_Nptr with values calculated from
# nodenums_Nptr

  resetAllIts2()
  # mesh iterator increment, retreval functions
  # TODO: checks if using mutable arrays is causing dynamic dispatch
  iterators_inc = [incrementVertIt, incrementEdgeIt, incrementFaceIt]
  iterators_get = [getVert, getEdge, getFace]
  num_entities = [mesh.numVert, mesh.numEdge, mesh.numEl]
  num_nodes_entity = mesh.numNodesPerType

  for etype = 1:3 # loop over entity types
#    println("etype = ", etype)
    if (num_nodes_entity[etype] != 0)  # if no nodes on this type of entity, skip
      for entity = 1:num_entities[etype]  # loop over all entities of this type
#	println("  entity number: ", entity)
	entity_ptr = iterators_get[etype]()  # get entity

	for node = 1:num_nodes_entity[etype]
#	  println("    node : ", node)
          nodenum = getNumberJ(mesh.nodenums_Nptr, entity_ptr, node-1, 0)
	  if nodenum != 0
	    for i=1:mesh.numDofPerNode
	      dofnum_i = (nodenum -1)*mesh.numDofPerNode + i
  	      numberJ(mesh.dofnums_Nptr, entity_ptr, node-1, i-1, dofnum_i)
	    end  # end loop over dofsPerNode
	  end   # end if nodenum != 0
	end  # end loop over node

      iterators_inc[etype]()
      end  # end loop over entitiesa
    end  # end if 
  end  # end loop over entity types

  resetAllIts2()

#  writeVtkFiles("dofs_numbered", mesh.m_ptr)
  return nothing
end



function populateNodeStatus(mesh::PumiMesh)
# populate the nodestatus_Nptr with values
# currently we set all nodes to status 3 (free)

  resetAllIts2()
  # mesh iterator increment, retreval functions
  # TODO: checks if using mutable arrays is causing dynamic dispatch
  iterators_inc = [incrementVertIt, incrementEdgeIt, incrementFaceIt]
  iterators_get = [getVert, getEdge, getFace]
  num_entities = [mesh.numVert, mesh.numEdge, mesh.numEl]
  num_nodes_entity = mesh.numNodesPerType

#  mesh.numNodesPerType = num_nodes_entity

  println("num_entities = ", num_entities)
  println("num_nodes_entity = ", num_nodes_entity)

  for etype = 1:3 # loop over entity types
#    println("etype = ", etype)
    if (num_nodes_entity[etype] != 0)  # if no nodes on this type of entity, skip
      for entity = 1:num_entities[etype]  # loop over all entities of this type
#	println("  entity number: ", entity)
	entity_ptr = iterators_get[etype]()  # get entity

	for node = 1:num_nodes_entity[etype]
#	  println("    node : ", node)
	  numberJ(mesh.nodestatus_Nptr, entity_ptr, node-1, 0, 3)
	end  # end loop over node

      iterators_inc[etype]()
      end  # end loop over entitiesa
    end  # end if 
  end  # end loop over entity types

  resetAllIts2()
  return nothing
end





function getDofNumbers(mesh::PumiMeshDG2)
# populate array of dof numbers, in same shape as solution array u (or q)

#println("in getDofNumbers")
#println("numNodesPerElement = ", mesh.numNodesPerElement)

mesh.dofs = Array(Int, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.shared_element_offsets[end] - 1)

for i=1:mesh.numEl
#  println("element ", i)
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
  numel = mesh.shared_element_offsets[i+1] - mesh.shared_element_offsets[i]
  dof_sendbuf[i] = Array(Int, mesh.numDofPerNode, mesh.numNodesPerElement, numel)
  sendbuf_i = dof_sendbuf[i]
  elnums_i = mesh.local_element_lists[i]
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
for i=1:3
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




function numberNodes(mesh::PumiMesh, number_dofs=false)
# assign node numbers to entire mesh, or dof numbers if number_dofs=true,
# using the correct Numbering pointer
# assumes mesh elements have already been reordered
# this works for high order elements
  println("Entered numberNodes")

  # calculate number of nodes, dofs
  num_nodes_v = countNodesOn(mesh.mshape_ptr, 0) # on vert
  # number of nodes on a vertex
  num_nodes_e = countNodesOn(mesh.mshape_ptr, 1) # on edge
  num_nodes_f = countNodesOn(mesh.mshape_ptr, 2) # on face
  numnodes = num_nodes_v*mesh.numVert + num_nodes_e*mesh.numEdge + num_nodes_f*mesh.numEl



  # initally number all dofs as numDof+1 to 2*numDof
  # this allows quick check to see if somthing is labelled or not

  resetAllIts2()
  # mesh iterator increment, retreval functions
  iterators_inc = [incrementVertIt, incrementEdgeIt, incrementFaceIt]
  iterators_get = [getVert, getEdge, getFace]
  num_entities = [mesh.numVert, mesh.numEdge, mesh.numEl]
  num_nodes_entity = [num_nodes_v, num_nodes_e, num_nodes_f]

#  mesh.numNodesPerType = num_nodes_entity

  if number_dofs
    println("numbering degrees of freedom")
    numbering_ptr = mesh.dofnums_Nptr
    curr_dof = mesh.numDof + 1
    dofpernode = mesh.numDofPerNode
    numDof = mesh.numDofPerNode*numnodes
  else  # do node numbering
    println("numbering nodes")
    numbering_ptr = mesh.nodenums_Nptr
    curr_dof = mesh.numNodes + 1
    dofpernode = 1
    numDof = numnodes
  end

  println("expected number of dofs = ", numDof)
  if (numDof > 2^30 || numDof < 0)
    println("Warning: too many dofs, renumbering will fail")
  end


#  println("num_entities = ", num_entities)
#  println("num_nodes_entity = ", num_nodes_entity)

#  curr_dof = 1
  curr_dof = numDof+1
  for etype = 1:3 # loop over entity types
#    println("etype = ", etype)
    if (num_nodes_entity[etype] != 0)  # if no nodes on this type of entity, skip
      for entity = 1:num_entities[etype]  # loop over all entities of this type
#	println("  entity number: ", entity)
	entity_ptr = iterators_get[etype]()  # get entity

	for node = 1:num_nodes_entity[etype]
#	  println("    node : ", node)
	  for dof = 1:dofpernode
	    numberJ(numbering_ptr, entity_ptr, node-1, dof-1, curr_dof)
#	    println("      entity ", entity_ptr, " labelled ", curr_dof)
	    curr_dof += 1
	  end  # end loop over dof
	end  # end loop over node

      iterators_inc[etype]()
      end  # end loop over entitiesa
    end  # end if 
  end  # end loop over entity types

  println("Finished initial numbering of dofs") 


  println("Performing final numbering of dofs")

  verts_i = Array(Ptr{Void}, 12)
  edges_i = Array(Ptr{Void}, 12)
  resetAllIts2()
  el_i_ptr = Ptr{Void}(0)  # hold current element
# TODO: move all if statements out one for loop (check only first dof on each node)
  curr_dof = 1
  for i=1:mesh.numEl
#    println("element number: ", i)
    el_i_ptr = getFace()
#    println("element pointer = ", el_i_ptr)
    incrementFaceIt()
    # get vertices, edges for this element
    numVert = getDownward(mesh.m_ptr, el_i_ptr, 0, verts_i)
#    println("verts_i = ", verts_i)
#    println("mesh.verts = ", mesh.verts)
    numEdge = getDownward(mesh.m_ptr, el_i_ptr, 1, edges_i)
    for j=1:3  # loop over vertices, edges
#      println("  vertex and edge number: ", j)
      vert_j = verts_i[j]
      edge_j = edges_i[j]
#      println("  vert_j = ", vert_j)
#      println("  edge_j = ", edge_j)
      for k=1:num_nodes_entity[1] # loop over vertex nodes
        for p=1:dofpernode
  #	println("    dof number: ", k)
  #	println("    vert_j = ", vert_j)
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
        pumi_node = mesh.nodemapSbpToPumi[k]

      for p=1:dofpernode  # loop over dofs
	dofnum_p = getNumberJ(numbering_ptr, el_i_ptr, pumi_node-1, p-1)
	if dofnum_p > numDof
	  numberJ(numbering_ptr, el_i_ptr, pumi_node-1, p-1, curr_dof)
	  curr_dof += 1
	end
      end
    end  # end loop over face nodes
  end  # end loop over elements

  resetAllIts2()


  println("finished performing final dof numbering")

  println("number of dofs = ", curr_dof - 1)
  if (curr_dof -1) != numDof 
    println("Warning: number of dofs assigned is not equal to teh expected number")
    println("number of dofs assigned = ", curr_dof-1, " ,expected number = ", numDof)
  else
    println("Dof numbering is sane")
  end



  resetAllIts2()
  return nothing

end

#=
function numberNodes(mesh::PumiMesh2, number_dofs=false)
# assign node numbers to entire mesh, or dof numbers if number_dofs=true,
# using the correct Numbering pointer
# assumes mesh elements have already been reordered
# this works for high order elements
  println("Entered numberDofs")

  # calculate number of nodes, dofs
  num_nodes_v = countNodesOn(mesh.mshape_ptr, 0) # on vert
  # number of nodes on a vertex
  num_nodes_e = countNodesOn(mesh.mshape_ptr, 1) # on edge
  num_nodes_f = countNodesOn(mesh.mshape_ptr, 2) # on face
  println("num_nodes_v = ", num_nodes_v)
  println("num_nodes_e = ", num_nodes_e)
  println("num_nodes_f = ", num_nodes_f)
  numnodes = num_nodes_v*mesh.numVert + num_nodes_e*mesh.numEdge + num_nodes_f*mesh.numEl



  # initally number all dofs as numDof+1 to 2*numDof
  # this allows quick check to see if somthing is labelled or not

  resetAllIts2()
  # mesh iterator increment, retreval functions
  iterators_inc = [incrementVertIt, incrementEdgeIt, incrementFaceIt]
  iterators_get = [getVert, getEdge, getFace]
  num_entities = [mesh.numVert, mesh.numEdge, mesh.numEl]
  num_nodes_entity = [num_nodes_v, num_nodes_e, num_nodes_f]

#  mesh.numNodesPerType = num_nodes_entity

  if number_dofs
    println("numbering degrees of freedom")
    numbering_ptr = mesh.dofnums_Nptr
    curr_dof = mesh.numDof + 1
    dofpernode = mesh.numDofPerNode
    numDof = mesh.numDofPerNode*numnodes
  else  # do node numbering
    println("numbering nodes")
    numbering_ptr = mesh.nodenums_Nptr
    curr_dof = mesh.numNodes + 1
    dofpernode = 1
    numDof = numnodes
  end

  println("expected number of dofs = ", numDof)
  if (numDof > 2^30 || numDof < 0)
    println("Warning: too many dofs, renumbering will fail")
  end


  println("num_entities = ", num_entities)
  println("num_nodes_entity = ", num_nodes_entity)

#  curr_dof = 1
  curr_dof = numDof+1
  for etype = 1:3 # loop over entity types
#    println("etype = ", etype)
    if (num_nodes_entity[etype] != 0)  # if no nodes on this type of entity, skip
      for entity = 1:num_entities[etype]  # loop over all entities of this type
#	println("  entity number: ", entity)
	entity_ptr = iterators_get[etype]()  # get entity

	for node = 1:num_nodes_entity[etype]
#	  println("    node : ", node)
	  for dof = 1:dofpernode
	    numberJ(numbering_ptr, entity_ptr, node-1, dof-1, curr_dof)
#	    println("      entity ", entity_ptr, " labelled ", curr_dof)
	    curr_dof += 1
	  end  # end loop over dof
	end  # end loop over node

      iterators_inc[etype]()
      end  # end loop over entitiesa
    end  # end if 
  end  # end loop over entity types

  println("Finished initial numbering of dofs") 


  println("Performing final numbering of dofs")

  verts_i = Array(Ptr{Void}, 12)
  edges_i = Array(Ptr{Void}, 12)
  resetAllIts2()
  el_i_ptr = Ptr{Void}(0)  # hold current element
# TODO: move all if statements out one for loop (check only first dof on each node)
  curr_dof = 1
  for i=1:mesh.numEl
#    println("element number: ", i)
    el_i_ptr = getFace()
    incrementFaceIt()
    # get vertices, edges for this element
    numVert = getDownward(mesh.m_ptr, el_i_ptr, 0, verts_i)
#    println("verts_i = ", verts_i)
#    println("mesh.verts = ", mesh.verts)
    numEdge = getDownward(mesh.m_ptr, el_i_ptr, 1, edges_i)
    for j=1:3  # loop over vertices, edges
#      println("  vertex and edge number: ", j)
      vert_j = verts_i[j]
      edge_j = edges_i[j]
#      println("  vert_j = ", vert_j)
#      println("  edge_j = ", edge_j)
      for k=1:dofpernode # loop over vertex dofs
#	println("    dof number: ", k)
#	println("    vert_j = ", vert_j)
        dofnum_k = getNumberJ(numbering_ptr, vert_j, 0, k-1)
	if dofnum_k > numDof  # still has initial number
	  # give it new (final) number
	  numberJ(numbering_ptr, vert_j, 0, k-1, curr_dof)
	  curr_dof += 1
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
      for p=1:dofpernode  # loop over dofs
	dofnum_p = getNumberJ(numbering_ptr, el_i_ptr, k-1, p-1)
	if dofnum_p > numDof
	  numberJ(numbering_ptr, el_i_ptr, k-1, p-1, curr_dof)
	  curr_dof += 1
	end
      end
    end  # end loop over face nodes
  end  # end loop over elements

  resetAllIts2()


  println("finished performing final dof numbering")

  println("number of dofs = ", curr_dof - 1)
  if (curr_dof -1) != numDof 
    println("Warning: number of dofs assigned is not equal to teh expected number")
    println("number of dofs assigned = ", curr_dof-1, " ,expected number = ", numDof)
  else
    println("Dof numbering is sane")
  end



  resetAllIts2()
  return nothing

end

=#
function getDofNumbers(mesh::PumiMesh2)
# populate array of dof numbers, in same shape as solution array u (or q)

println("in getDofNumbers")
println("numNodesPerElement = ", mesh.numNodesPerElement)

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

