# file for functions related to the geometry model underlying the mesh

"""
  This function gets an array describing the mesh edges (faces in 3D) on a set
  of geometric edges (faces in 3D).

  Inputs:
    mesh: a mesh object (2d or 3d)
    medges: an array of geometric mesh edges to get the mesh edges on
    offset: the current index in boundary_nums (TODO: rename this: its an
            index not an offset). This value should be 1 the first time this 
            function is called

  Inputs/Outputs:
    boundary_nums: a two dimensional array, 2 x mesh.numBoundaryFaces will be
                   updated with [element number, global face number] for each
                   mesh edge classified on medges

  Outputs:
    new_offset: the offset that should be passed in the next time this
                function is called
    print_warning: a Bool indicated whether an edge classified on medges is
                   periodic.  This function is typically used to get the mesh
                   edges to apply boundary conditions to.  Applying boundary
                   conditions to periodic faces is not allowed.  These
                   edges are not added to boundary_nums

  The setup with the offset and boundary_nums allows this function to be called
  repeatedly with different values in medges to get the mesh edges with
  different boundary conditions applied to them.
  TOOD: use sview to remove this machinery

"""
function getMeshEdgesFromModel{T}(mesh::PumiMesh, medges::AbstractArray{Int, 1}, offset::Integer, boundary_nums::AbstractArray{T, 2})
# get the global numbers of the mesh edges lying on the model edges in medges
# medges typically coresponds to all the model edges that have a particular
# boundary condition applied to them
# offset is the index in boundary_nums to start with
# this allows populating the entire array without making temporary copies
# offset should start as 1, not zero

#  bndry_edges = zeros(Int, mesh.numBoundaryFaces)
  index = 0  # relative index in bndry_edges
  faces = Array(Ptr{Void}, 2)  # an edge can have maximum 2 faces using it
  print_warning = false
  for i=1:mesh.numFace
    edge_i = mesh.faces[i]
    me_i = toModel(mesh.m_ptr, edge_i)
    me_dim = getModelType(mesh.m_ptr, me_i)
    me_tag = getModelTag(mesh.m_ptr, me_i)
    if me_dim == (mesh.dim-1)  # face (edge in 2d)
      onBoundary = findfirst(medges, me_tag)

      # check if face is periodic
      isPeriodic = countMatches(mesh.m_ptr, edge_i) > 0
      print_warning = print_warning || (isPeriodic && onBoundary != 0)

      if onBoundary != 0 && !isPeriodic # if mesh face is on any of the  model faces
        
	# get face number
        numFace = countAdjacent(mesh.m_ptr, edge_i, mesh.dim)  # should be count upward

	@assert( numFace == 1)

        getAdjacent(faces)
        facenum = getNumberJ(mesh.el_Nptr, faces[1], 0, 0) + 1
#        facenum = getFaceNumber2(faces[1]) + 1
#        edgenum = getEdgeNumber2(edge_i)  # unneeded?

	boundary_nums[offset + index, 1] = facenum
        boundary_nums[offset + index, 2] = i

	index += 1
      end
    end
  end

  return offset + index, print_warning

end  # end function

"""
  This function counts the number of mesh edges (faces in 3D) in several
  catagories

  Inputs:
    mesh: the Mesh object (2D or 3D)

  Outputs:
    bnd_edges_cnt:  The total number of mesh edges classified on 
                    geometric edges.
                    A set of matched edges (ie. periodic edges) are counted
                    only once (ie. treated as a single edge).
    internal_edge_cnt:  The number of edges than are shared by two element,
                        either because the edges are internal or they are 
                        periodic.  As with bnd_edges_cnt, a set of matched edges
                        is counted only once.
    periodic_edge_cnt: The number of periodic edges.  As with bnd_edges_cnt,
                       each pair of matched edges is counted only once.

    geo_edge_nums: numbers of the geometric edges that have non-periodic
                   mesh edges on them.
"""
function countBoundaryEdges(mesh::PumiMesh)
  # count boundary edges by checking if their model edge has a BC
  # count number of external edges by checking the number of upward adjacencies
  # store array of [element number, global edge number]a
#  resetAllIt()

  bnd_edges_cnt = 0  # number of edges with boundary conditions
  internal_edge_cnt = 0 # count the number of internal interfaces
  periodic_edge_cnt = 0  # number of *pairs* of periodic BCs
  geo_edge_nums = Array(Int, 0)
#  bnd_edges = Array(Int, mesh.numEdge, 2)
  elements = Array(Ptr{Void}, 2)  # edge has maximum 2 faces
  part_nums = Array(Cint, 1)
  matched_entities = Array(Ptr{Void}, 1)
  seen_matches = Set{Ptr{Void}}()
  for i=1:mesh.numFace
    edge_i = mesh.faces[i]

    # skip if the other half of this match has been seen
    if in(edge_i, seen_matches) 
      continue
    end
#    edge_i = getEdge()

    # get  model edge info
    me_i = toModel(mesh.m_ptr, edge_i)
    me_dim = getModelType(mesh.m_ptr, me_i)
    me_tag = getModelTag(mesh.m_ptr, me_i)

    # get mesh face info
    numEl = countAdjacent(mesh.m_ptr, edge_i, mesh.dim)  # should be count upward 
    getAdjacent(elements)
    elnum = getNumberJ(mesh.el_Nptr, elements[1], 0, 0) + 1


    nremotes = countRemotes(mesh.m_ptr, edge_i)
    nmatches = countMatches(mesh.m_ptr, edge_i)
    getMatches(part_nums, matched_entities)
    has_local_match = (nmatches > 0) && part_nums[1] == mesh.myrank

    println("edge ", i, " numEl = ", numEl)
    println("nmatches = ", nmatches)
    println("nremotes = ", nremotes)
    println("has_local_match = ", has_local_match)
    println("first element = ", elnum)

    if numEl == 2 || (has_local_match)  # internal interfaces (not including shared parallel edges)
      internal_edge_cnt += 1
    end

    if has_local_match
      periodic_edge_cnt += 1
      push!(seen_matches, matched_entities[1])
    end

    if me_dim == (mesh.dim-1) && nmatches == 0  # if classified on model edge
#      index = findfirst(bndry_edges_all, me_tag)
      # accumulate all geometric edges with non-matched mesh edges
      if !(me_tag in geo_edge_nums)
        println("adding edge ", me_tag, " to geo_edge_nums")
        println("nmatches = ", nmatches)
        println("nremotes = ", nremotes)
        println("elnum = ", elnum)
        push!(geo_edge_nums, me_tag)
      end

#      if index != 0  # if model edge has a BC on i
#	facenum = getFaceNumber2(faces[1]) + 1

	bnd_edges_cnt += 1

        # bnd_edges array is unused?
#	bnd_edges[bnd_edges_cnt, 1] = elnum
#	bnd_edges[bnd_edges_cnt, 2] = i
#      end
    end  # end if me_dim == 1
#    incrementEdgeIt()

  end  # end for loop


#  mesh.boundary_nums = bnd_edges[1:bnd_edges_cnt, :] # copy, bad but unavoidable
#  mesh.numBoundaryFaces = bnd_edges_cnt

return bnd_edges_cnt, internal_edge_cnt, periodic_edge_cnt, geo_edge_nums

end  # end function

#=
function countBoundaryEdges(mesh::PumiMesh3, bndry_edges_all)
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
#  mesh.numBoundaryFaces = bnd_edges_cnt

return bnd_edges_cnt, external_edges_cnt

end  # end function
=#
