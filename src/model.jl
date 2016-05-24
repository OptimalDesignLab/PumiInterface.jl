# file for functions related to the geometry model underlying the mesh

function getMeshEdgesFromModel{T}(mesh::PumiMesh, medges::AbstractArray{Int, 1}, offset::Integer, boundary_nums::AbstractArray{T, 2})
# get the global numbers of the mesh edges lying on the model edges in medges
# offset is the index in boundary_nums to start with
# this allows populating the entire array without making temporary copies
# offset should start as 1, not zero

#  bndry_edges = zeros(Int, mesh.numBoundaryEdges)
  index = 0  # relative index in bndry_edges
  faces = Array(Ptr{Void}, 2)  # an edge can have maximum 2 faces using it
  for i=1:mesh.numEdge
    edge_i = mesh.edges[i]
    me_i = toModel(mesh.m_ptr, edge_i)
    me_dim = getModelType(mesh.m_ptr, me_i)
    me_tag = getModelTag(mesh.m_ptr, me_i)
    if me_dim == 1  # edge
      onBoundary = findfirst(medges, me_tag)

      if onBoundary != 0  # if mesh edge is on model edge
        
	# get face number
        numFace = countAdjacent(mesh.m_ptr, edge_i, 2)  # should be count upward

	@assert( numFace == 1)

        getAdjacent(faces)
        facenum = getFaceNumber2(faces[1]) + 1
        edgenum = getEdgeNumber2(edge_i)  # unneeded?

	boundary_nums[offset + index, 1] = facenum
        boundary_nums[offset + index, 2] = i

	index += 1
      end
    end
  end

  return offset + index

end  # end function

function countBoundaryEdges(mesh::PumiMeshDG2, bndry_edges_all)
  # count boundary edges by checking if their model edge has a BC
  # count number of external edges by checking the number of upward adjacencies
  # store array of [element number, global edge number]
  resetEdgeIt()
  bnd_edges_cnt = 0
  external_edges_cnt = 0
  internal_edge_cnt = 0 # count the number of internal interfaces
  bnd_edges = Array(Int, mesh.numEdge, 2)
  faces = Array(Ptr{Void}, 2)  # edge has maximum 2 faces
  for i=1:mesh.numEdge
    edge_i = getEdge()

    # get  model edge info
    me_i = toModel(mesh.m_ptr, edge_i)
    me_dim = getModelType(mesh.m_ptr, me_i)
    me_tag = getModelTag(mesh.m_ptr, me_i)

    # get mesh face info
    numFace = countAdjacent(mesh.m_ptr, edge_i, 2)  # should be count upward
    nremotes = countRemotes(mesh.m_ptr, edge_i)

    if numFace == 1 && nremotes == 0  # external edges
      external_edges_cnt += 1
    elseif nremotes == 0  # internal interfaces (not including shared parallel edges)
      internal_edge_cnt += 1
    end

    if me_dim == 1  # if classified on model edge
      index = findfirst(bndry_edges_all, me_tag)

      if index != 0  # if model edge has a BC on i

	getAdjacent(faces)
	facenum = getFaceNumber2(faces[1]) + 1

	bnd_edges_cnt += 1
	bnd_edges[bnd_edges_cnt, 1] = facenum
	bnd_edges[bnd_edges_cnt, 2] = i
      end
    end  # end if me_dim == 1
    incrementEdgeIt()

  end  # end for loop


#  mesh.boundary_nums = bnd_edges[1:bnd_edges_cnt, :] # copy, bad but unavoidable
#  mesh.numBoundaryEdges = bnd_edges_cnt

return bnd_edges_cnt, external_edges_cnt, internal_edge_cnt

end  # end function


function countBoundaryEdges(mesh::PumiMesh2, bndry_edges_all)
  # count boundary edges by checking if their model edge has a BC
  # count number of external edges by checking the number of upward adjacencies
  # store array of [element number, global edge number]
  resetEdgeIt()
  bnd_edges_cnt = 0
  external_edges_cnt = 0
  bnd_edges = Array(Int, mesh.numEdge, 2)
  faces = Array(Ptr{Void}, 2)  # edge has maximum 2 faces
  for i=1:mesh.numEdge
    edge_i = getEdge()

    # get  model edge info
    me_i = toModel(mesh.m_ptr, edge_i)
    me_dim = getModelType(mesh.m_ptr, me_i)
    me_tag = getModelTag(mesh.m_ptr, me_i)

    # get mesh face info
    numFace = countAdjacent(mesh.m_ptr, edge_i, 2)  # should be count upward
    if numFace == 1  # external edges
      external_edges_cnt += 1
    end

    if me_dim == 1  # if not classified on model edge
      index = findfirst(bndry_edges_all, me_tag)



      if index != 0  # if model edge has a BC on i

	getAdjacent(faces)
	facenum = getFaceNumber2(faces[1]) + 1



	bnd_edges_cnt += 1
	bnd_edges[bnd_edges_cnt, 1] = facenum
	bnd_edges[bnd_edges_cnt, 2] = i
      end
    end  # end if me_dim == 1
    incrementEdgeIt()

  end  # end for loop


#  mesh.boundary_nums = bnd_edges[1:bnd_edges_cnt, :] # copy, bad but unavoidable
#  mesh.numBoundaryEdges = bnd_edges_cnt

return bnd_edges_cnt, external_edges_cnt

end  # end function



