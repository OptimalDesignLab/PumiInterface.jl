# file for functions related to calculating the sparsity pattern of a matrix

# can be generalized with some topology constants
function getSparsityCounts(mesh::PumiMeshDG, sparse_bnds::AbstractArray{Int32, 2}; getdofs=true)
# count how many local, remote nodes/dofs each dof is related to

  numNodesPerElement = mesh.numNodesPerElement
  if getdofs || mesh.numDofPerNode == 1
    numDofPerNode = mesh.numDofPerNode
    div_factor = 1
    add_factor = 0
  else
    numDofPerNode = 1
    div_factor = mesh.numDofPerNode
    add_factor = 1
  end

  edges = Array{Ptr{Void}}(4)
  part_nums = Array{Cint}(4)
  shared_edges = Array{Ptr{Void}}(4)
  for i=1:mesh.numEl
    el_ptr = mesh.elements[i]
    apf.getDownward(mesh.m_ptr, el_ptr, mesh.dim-1, edges)
    nremotes = 0
    for j=1:mesh.numFacesPerElement
      edge_j = edges[j]
      ncopies = apf.countCopies(mesh.shr_ptr, edge_j)
      apf.getCopies(part_nums, shared_edges)
      for k=1:ncopies
        if part_nums[k] != mesh.myrank
          nremotes += 1
        end
      end
    end

    nlocal = mesh.numFacesPerElement + 1  # 3 neighbors + 1
    nlocal = nlocal - nremotes

    # apply this to all dofs on this element
    for j=1:mesh.numNodesPerElement
      for k =1:numDofPerNode
        dof_k = mesh.dofs[k, j, i]
        dof_k = div(dof_k, div_factor) + add_factor
        sparse_bnds[1, dof_k] = nlocal*numNodesPerElement*numDofPerNode
        sparse_bnds[2, dof_k] = nremotes*numNodesPerElement*numDofPerNode
      end
    end
  end

  return nothing
end

# WIP: tighter sparsity bounds for CG
#=
function getDofConnectivity(mesh::PumiMesh2)

  # find the maximum number of elements a dof is related to
  max_el = 0
  up_cnt = 400 # maximum number of upward adjacencies
  up = Array{Ptr{Void}}(up_cnt)  # equivalent of apf::Up
  verts = Array{Ptr{Void}}(3)
  nel = 0
  els_all = Array{Ptr{Void}}(3*up_cnt)
  els_tmp = Array{Ptr{Void}}(3*up_cnt)
  for i=1:mesh.numEl
    apf.getDownward(mesh.m_ptr, 0, verts)
    pos = 1
    for j=1:3
      nel[j] = apf.countAdjacent(mesh.m_ptr, verts[j], 2)
      apf.getAdjacent(up)
      # copy element pointer into a single array
      for k=1:nel[j]
        els_all[pos] = up[k]
        pos += 1
      end
    end
    numel = pos - 1
    
    # count number of non duplicates
    els_used = sview(els_all, 1:numel);
    sort!(els_used)
    cnt = 1
    for j=2:length(els_used)
      if els_used[j] != els_used[j-1]
        cnt += 1
      end
    end
    if cnt > max_el
      max_el = cnt
    end
  end

  dof_elements = zeros(Int32, max_el, mesh.numNodes)
  # now get the element numbers
  for i=1:mesh.numEl
    apf.getDownward(mesh.m_ptr, 0, verts)
    pos = 1
    for j=1:3
      nel[j] = apf.countAdjacent(mesh.m_ptr, verts[j], 2)
      apf.getAdjacent(up)
      # copy element pointer into a single array
      for k=1:nel[j]
        els_all[pos] = up[k]
        pos += 1
      end
    end
    numel = pos - 1

    # de-duplicate whle copying into node_elements
    els_used = sview(els_all, 1:numel);
    sort!(els_used)
    els_tmp[1] = els_used[1]
    pos = 2
    for j=1:length(els_used)
      if els_used[j] != els_used[j-1]
        els_tmp[pos] = els_used[j]
      end
    end

  end
end
=#
 





    # copy the element pointers into a single array


    
# or maybe replace this with a tighter sparsity bound function
# can be generalized with a few if statements
function getSparsityBounds(mesh::PumiMesh, sparse_bnds::AbstractArray{Int32, 2}; getdofs=true)
# sparse_bnds : 2 by numDof array of the minimum, maximum dofs connected to 
# a dof in the Jacobian
# or node numbers, if getdofs=false
# this works for high order elements despite not using mesh.elementNodeOffsets 
# because it does not reference what element a particular entity belongs to
# getDofBounds also does not need it because only the minimum and maximum
# are important, not the order of them

#resetAllIts2(mesh.m_ptr)
# mesh iterator increment, retreval functions
#iterators_inc = [incrementVertIt, incrementEdgeIt, incrementFaceIt]
#iterators_inc = (VertIterator(), EdgeIterator(), FaceIterator(), ElIterator())

#iterators_get = [getVert, getEdge, getFace]
#iterators_get = (VertGetter(), EdgeGetter(), FaceGetter(), ElGetter)
num_entities = mesh.numEntitiesPerType
#num_entities = [mesh.numVert, mesh.numEdge, mesh.numEl]
num_nodes_entity = mesh.numNodesPerType  # number of nodes on each type
                                         # of mesh entity

if getdofs
  numbering_Nptr = mesh.dofnums_Nptr
  numDofPerNode = mesh.numDofPerNode
  @assert size(sparse_bnds, 2) == mesh.numDof
else
  numbering_Nptr = mesh.nodenums_Nptr
  numDofPerNode = 1
  @assert size(sparse_bnds, 2) == mesh.numNodes
end

for etype=1:(mesh.dim+1)  # loop over mesh entity types
 
  if (num_nodes_entity[etype] != 0)  # there are nodes herea
    it = apf.MeshIterator(mesh.m_ptr, etype-1)
    for entity = 1:num_entities[etype]  # loop over all entities of this type
#      entity_ptr = iterators_get[etype]()  # get pointer to mesh entity
      entity_ptr = apf.iterate(it)

      min, max = getDofBounds(mesh, entity_ptr, getdofs=getdofs)
      for node = 1:num_nodes_entity[etype]  # loop over nodes on each entity
	# get the minimum, maximum related dof numbers
	# use the same min, max for all dofs on this node

	for dof = 1:numDofPerNode
	  dofnum = apf.getNumberJ(numbering_Nptr, entity_ptr, node - 1, dof-1)
          sparse_bnds[1, dofnum] = min
	  sparse_bnds[2, dofnum] = max

	end  # end loop over dofs
      end  # end loop over nodes on entity
#      iterators_inc[etype]()
    end  # end loops over entities of this type
    apf.free(it)
  end  # end if statement
end  # end loop over entity types

return nothing

end

# rewrite this because it is terrible
# this could be generalized 3d?
function getDofBounds(mesh::PumiMesh2D, entity_ptr::Ptr{Void}; getdofs=true) 
# gets the maximum, minimum dofs associated with the entity currently
# pointed to by the iterator specified by etype
# getDofs = true -> get dof numbers, false -> get node numbers
# this works for distance-0 and 1 colorings

#iterators_get = [getVert, getEdge, getFace]
#it = apf.MeshIterator(mesh.m_ptr, etype - 1)
#entity_ptr = iterators_get[etype]()

# get associated elements (distance-0 elements)
num_adj = apf.countAdjacent(mesh.m_ptr, entity_ptr, 2)
el_arr = apf.getAdjacent(num_adj)

# distance-1
#=
dofnums = zeros(Int32, mesh.numDofPerNode, mesh.numNodesPerElement, num_adj)

for i=1:num_adj
  el_i = apf.getNumberJ(mesh.el_Nptr, el_arr[i], 0, 0) + 1
  sub_arr = sub(dofnums, :, :, i)
  getGlobalNodeNumbers(mesh, el_i, sub_arr)
end

min, max = getMinandMax(dofnums)

return min, max
=#

  # get distance-2 elements
  # this is really distance-1 
  if mesh.coloringDistance >= 2
    edge_arr = Array{Ptr{Void}}(num_adj*3)  # enough space for all edges, including repeats
    for i=1:num_adj  # get the edges
      sub_arr = sview(edge_arr, (3*(i-1) + 1):(3*i))
      apf.getDownward(mesh.m_ptr, el_arr[i], 1, sub_arr)
    end

    # edge_arr now populated with all edges

    # print the edge numbers
  #=
    for i=1:num_adj*3
      edge_num_i = apf.getNumberJ(mesh.edge_Nptr, edge_arr[i], 0, 0)
      println("edge ", i, " has number ", edge_num_i)
    end
  =#
    # now get the elements the edges belong to
    # count the number of elements first, then get them
    num_els = zeros(Int, length(edge_arr) + 1)  # count number of elements each edge has

  #  println("length(num_els) = ", length(num_els))
  #  println("length(edge_arr) = ", length(edge_arr))
    for i=1:length(edge_arr)
      num_els[i] = apf.countAdjacent(mesh.m_ptr, edge_arr[i], 2)
    end

    # now get the elements
    num_adj = sum(num_els)
    el_arr = Array{Ptr{Void}}(num_adj)  # rebind el_arr to new array
    start_idx = 1
    end_idx = num_els[1]
    for i=1:length(edge_arr)
      edge_i = edge_arr[i]
      sub_arr = sview(el_arr, start_idx:end_idx)
      apf.countAdjacent(mesh.m_ptr, edge_i, 2)
      apf.getAdjacent(sub_arr)

      # update indices
      start_idx = end_idx + 1
      end_idx += num_els[i+1]
    end

    # now el_arr has all the elements, including repeats
    # num_adj = legnth(el_arr)

  end  # end if distance-2



    
  # get dofnums for the elements
  if getdofs
    numDofPerNode = mesh.numDofPerNode
  else  # get node numbers
    numDofPerNode = 1 
  end

  dofnums = zeros(Int32, numDofPerNode, mesh.numNodesPerElement, num_adj)

  for i=1:num_adj
    elnum_i = apf.getNumberJ(mesh.el_Nptr, el_arr[i], 0, 0) + 1
    getGlobalNodeNumbers(mesh, elnum_i, sview(dofnums, :, :, i), getdofs=getdofs)
  end

  min, max = getMinandMax(dofnums)

  return min, max

end

