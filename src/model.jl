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
    boundary_nums: a two dimensional array, mesh.numBoundaryFaces x 2 will be
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
function getMeshEdgesFromModel(mesh::PumiMesh, medges::AbstractArray{I, 1}, offset::Integer, boundary_nums::AbstractArray{T, 2}) where {T, I<:Integer}
# get the global numbers of the mesh edges lying on the model edges in medges
# medges typically coresponds to all the model edges that have a particular
# boundary condition applied to them
# offset is the index in boundary_nums to start with
# this allows populating the entire array without making temporary copies
# offset should start as 1, not zero

#  bndry_edges = zeros(Int, mesh.numBoundaryFaces)
  index = 0  # relative index in bndry_edges
  faces = Array{Ptr{Void}}(2)  # an edge can have maximum 2 faces using it
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
  geo_edge_nums = Array{Int}(0)

  # temporary arrays
  elements = Array{Ptr{Void}}(2)  # edge has maximum 2 faces
  part_nums = Array{Cint}(1)
  matched_entities = Array{Ptr{Void}}(1)
  seen_matches = Set{Ptr{Void}}()
  for i=1:mesh.numFace
    edge_i = mesh.faces[i]

    # skip if the other half of this match has been seen
    if in(edge_i, seen_matches) 
      continue
    end

    # get  model edge info
    me_i = toModel(mesh.m_ptr, edge_i)
    me_dim = getModelType(mesh.m_ptr, me_i)
    me_tag = getModelTag(mesh.m_ptr, me_i)

    # get mesh face info
    numEl = countAdjacent(mesh.m_ptr, edge_i, mesh.dim)  # should be count upward 
    nremotes = countRemotes(mesh.m_ptr, edge_i)
    nmatches = countMatches(mesh.m_ptr, edge_i)
    getMatches(part_nums, matched_entities)
    has_local_match = (nmatches > 0) && part_nums[1] == mesh.myrank

    # internal interfaces (not including shared parallel edges)
    if numEl == 2 || (has_local_match)  
      internal_edge_cnt += 1
    end

    if has_local_match
      periodic_edge_cnt += 1
      push!(seen_matches, matched_entities[1])
    end

    if me_dim == (mesh.dim-1) && nmatches == 0  # if classified on model edge
      getAdjacent(elements)
      elnum = getNumberJ(mesh.el_Nptr, elements[1], 0, 0) + 1
      bnd_edges_cnt += 1

      # accumulate all geometric edges with non-matched mesh edges
      if !(me_tag in geo_edge_nums)
        push!(geo_edge_nums, me_tag)
      end
    end  # end if me_dim == 1

  end  # end for loop

  return bnd_edges_cnt, internal_edge_cnt, periodic_edge_cnt, geo_edge_nums

end  # end function

"""
  Returns an array of Boundary objects for the faces on the specified
  geometric face (edges in 2D)

  **Inputs**
  
   * mesh: a PumiMesh object
   * geo_face_nums: an array of the geometric face numbers (edges in 2D)
   * geo_region_nums: an array of the geometric regions (faces in 2D)
                      corresponding to geo_face_nums.  This
                      is only required if the geometric edge is on the interior
                      of the domain, to uniquely identify the elements.

  **Outputs**

   * an array of Boundary objects
"""
function getBoundaries(mesh::PumiMesh, geo_face_nums::Array{Int, 1},
                       geo_region_nums::Array{Int, 1} = Array{Int}(0))

  bndries = Array{Boundary}(0)
  adjacent_els = Array{Ptr{Void}}(400)  # apf::Up

  #TODO: check input arrays for uniqueness

  for i=1:mesh.numFace
    edge_i = mesh.faces[i]

    me_i = toModel(mesh.m_ptr, edge_i)
    me_dim = getModelType(mesh.m_ptr, me_i)
    me_tag = getModelTag(mesh.m_ptr, me_i)

    if me_dim == mesh.dim-1
      idx = findfirst(geo_face_nums, me_tag)
      if idx != 0  # this mesh edge is on one of the geometric edges

        # get mesh face info
        numEl = countAdjacent(mesh.m_ptr, edge_i, mesh.dim)  # should be count upward 
        if length(geo_region_nums) == 0
          @assert numEl == 1
          getAdjacent(adjacent_els)
          el_ptr = adjacent_els[1]
          elnum = getNumberJ(mesh.el_Nptr, el_ptr, 0, 0) + 1
        else
          @assert numEl <= 400
          el_tag = geo_region_nums[idx]
          el_count = 0 # number of elements matching
          elnum = 0

          # find the element that is on the geo_region_nums
          for j=1:numEl
            el_j = adjacent_els[j]
            me_i = toModel(mesh.m_ptr, el_j)
            me_dim = getModelType(mesh.m_ptr, me_i)
            me_tag = getModelTag(mesh.m_ptr, me_i)

            if me_dim == mesh.dim && me_tag == el_tag
              elnum = getNumberJ(mesh.el_Nptr, el_j, 0, 0) + 1
              el_count += 1
            end
          end

          # verify exactly one element matches
          @assert el_count == 1
        end  # end if length == 0

        # get local number of the edge
        facenum_local = getFaceLocalNum(mesh, i, elnum)
        push!(bndries, Boundary(elnum, facenum_local))
      end  # end if idx != 0
    end  # end if me_dim
  end  # end loop over faces

  return bndries
end







