# interface functions 

"""
  Numbers all points on the specified surfaces.  All points on the
  surface of the domain not specified are numbered with the number of
  points on the specified surface + 1 (the numbering is 1-based)

  Currently this supports linear coordinate fields only

  **Inputs**

   * mesh
   * bc_nums: array of boundary condition numbers that define the surface(s)
              to number
   * isglobal: if true, do a globally consistent numbering of surface points,
               otherwise, do part-local numbering, default true
   * numbering_name: the name of the Pumi numbering to be created, any existing
                     numbering with this name will be deleted.
                     Default "warp_surf_nums"

  **Outputs**

   * numFaceNodes: number of nodes on the specified surface(s)
   * n_face: a Ptr{Void} (really an apf::Numbering*).  Note that this is
             not an apf::GlobalNumbering*, so be careful not to overflow
             a Cint.
   * face_verts: array of apf::MeshEntity* in the order they appear in the
                 surface numbering, length numFaceNodes
"""
function numberSurfacePoints(mesh::PumiMeshDG, bc_nums::AbstractVector{I}, isglobal::Bool=false, numbering_name::AbstractString="warp_surf_nums") where I<:Integer

  #TODO: should face_verts be only the MeshEntities owned by this part, or
  #      all the ones present on this part?

  @assert mesh.coord_order <= 2
  for i in bc_nums  # check that the BC numbers are valid
    @assert i <= length(mesh.bndry_offsets - 1)
    @assert i > 0
  end
  # number the (unique) coordinate nodes on the faces

  # we can't guarantee an existing numbering with the same name was
  # created using the same bc_nums, so delete any existing numbering
  n_old = findNumbering(mesh.m_ptr, numbering_name)
  if n_old != C_NULL
    destroyNumbering(n_old)
  end
  n_face = createNumberingJ(mesh.m_ptr, numbering_name, 
                             mesh.coordshape_ptr, 1)
  topo = mesh.topo
  num_i = 1
  verts = Array{Ptr{Void}}(4)
  edges = Array{Ptr{Void}}(12)  # all edges of the element
  face_verts = Array{Ptr{Void}}(0)  # TODO: add sizehint
#  coords = zeros(Float64, 3)
  for i in bc_nums
    start_idx = mesh.bndry_offsets[i]
    end_idx = mesh.bndry_offsets[i+1] - 1
    rng = start_idx:end_idx

    for j in rng
      bndry_j = mesh.bndryfaces[j]
      el_j = bndry_j.element
      el_ptr = mesh.elements[el_j]

      # get the vertices
      getDownward(mesh.m_ptr, el_ptr, 0, verts)

      for k=1:size(mesh.topo.face_verts, 1)
        v_k = verts[mesh.topo.face_verts[k, bndry_j.face]]

        if !isNumbered(n_face, v_k, 0, 0)
#          getPoint(mesh.m_ptr, v_k, 0, coords)
          numberJ(n_face, v_k, 0, 0, num_i)
          push!(face_verts, v_k)
          num_i += 1
        end
      end  # end loop k

      # get the edges nodes too
      if hasNodesIn(mesh.coordshape_ptr, 1)
        getDownward(mesh.m_ptr, el_ptr, 1, edges)
        for k=1:size(topo.face_edges, 1)
          edge_k = edges[topo.face_edges[k, bndry_j.face]]

          # up to 2nd order fields, we don't need to worry about edge orientation
          if !isNumbered(n_face, edge_k, 0, 0)

            numberJ(n_face, edge_k, 0, 0, num_i)
            push!(face_verts, edge_k)
            num_i += 1
          end  # end if
        end  # end loop k
      end  # end if

    end  # end loop j
  end  # end loop i

  # number all remaining points with the number of face verts + 1
  #TODO: this should be a loop 1:dim, iterative of meshentities of this dim
  for i=1:mesh.numVert
    v_i = mesh.verts[i]

    if !isNumbered(n_face, v_i, 0, 0)
      numberJ(n_face, v_i, 0, 0, num_i)
    end
  end

  if hasNodesIn(mesh.coordshape_ptr, 1)
    for i=1:mesh.numEdge
      edge_i = mesh.edges[i]
      if !isNumbered(n_face, edge_i, 0, 0)
        numberJ(n_face, edge_i, 0, 0, num_i)
      end
    end
  end

  numFacePts = num_i - 1

  return numFacePts, n_face, face_verts
end


#------------------------------------------------------------------------------
# Reduction operations

"""
  Assignment reduction (always returns the right hand value)
"""
struct AssignReduction{T} <: Reduction{T}
  neutral_element::T

  function AssignReduction{T}() where {T}
    return new(0)
  end
end

function (obj::AssignReduction)(valL::Number, valR::Number)
  return valR
end

"""
  Sum reduction
"""
struct SumReduction{T} <: Reduction{T}
  neutral_element::T

  function SumReduction{T}() where {T}
    return new(0)
  end
end

function (obj::SumReduction)(valL::Number, valR::Number)
  return valL + valR
end

"""
  Minimum reduction
"""
struct MinReduction{T} <: Reduction{T}
  neutral_element::T

  function MinReduction{T}() where {T}
    return new(typemax(T))
  end
end

function (obj::MinReduction)(valL::Number, valR::Number)
  return min(valL, valR)
end

"""
  Max reduction
"""
struct MaxReduction{T}  <: Reduction{T}
  neutral_element::T

  function MaxReduction{T}() where {T}
    return new(typemin(T))
  end
end

function (obj::MaxReduction)(valL::Number, valR::Number)
  return max(valL, valR)
end



"""
  This function applies a user-supplied operation to take an array the same
  shape as `mesh.vert_coords` and transform it into a vector of length
  `mesh.dim*mesh.coord_numNodes` (a flattened array of all the coordinates).
  The user supplied operation is applied pairwise to all the values that are
  mapped to the same location in the vector.  `mesh.coord_nodenums_Nptr` is
  used to define the mapping from the array to the vector.

  The user supplied operation must have the signature

  ```
    val = reduce_op(valL, valR)
  ```

  where `valL` is the value that is already in the vector, and `valR` is
  the new value from `coords_arr`.  The vector is initialized with 
  `reduce_op.neutral_element`

  **Inputs**
  
   * mesh: DG mesh object
   * coords_arr: array, same shape as `mesh.vert_coords`
   * reduce_op: a [`Reduction`](@ref) object, defaults to [`SumReduction`](@ref)


  **Inputs/Outputs**

   * coords_vec: vector to be populated

  **Keyword Arguments**

   * parallel: if true (default), do the reduction in parallel, sending data
               to the owning process.  The entry in the local vector for
               non-owned entities will be zero (in general,
               `reduce_op.neutral_element`)
"""
function coords3DTo1D(mesh::PumiMeshDG, coords_arr::AbstractArray{T, 3},
                     coords_vec::AbstractVector, reduce_op::Reduction{T}=SumReduction{T}(); parallel=true)  where {T}

  @assert mesh.coord_order <= 2
  @assert size(coords_arr, 3) == mesh.numEl
  @assert size(coords_arr, 2) == mesh.coord_numNodesPerElement
  @assert size(coords_arr, 1) == mesh.dim
  @assert length(coords_vec) == mesh.coord_numNodes*mesh.dim

  _parallel::Bool = parallel

  if _parallel
    sendParallelData(mesh.coordscatter, coords_arr, reduce_op)
  end

  fill!(coords_vec, reduce_op.neutral_element)
  node_entities = ElementNodeEntities(mesh.m_ptr, mesh.coordshape_ptr, mesh.dim)

 
  #TODO: not sure if this gives enough time for data to arrive, maybe combine
  #      with calcCoordinatesAndMetrics_rev?
  shr = mesh.normalshr_ptr
  for i=1:mesh.numEl
    getNodeEntities(node_entities, mesh.elements[i])
    for j=1:mesh.coord_numNodesPerElement
      entity = node_entities.entities[j]

      if _parallel && getOwner(shr, entity) == mesh.myrank
        for k=1:mesh.dim
          idx = getNumberJ(mesh.coord_nodenums_Nptr, entity, 0, k-1)

          coords_vec[idx] = reduce_op(coords_vec[idx], coords_arr[k, j, i])
        end
      end  # end if
    end  # end j
  end  # end i

  if _parallel
    # lambda function for receiving
    calc_func = (data::PeerData) -> receiveVecFunction(data, mesh, coords_vec,
                                                       reduce_op)

    receiveParallelData(mesh.coordscatter, calc_func)
  end

  return nothing
end

"""
  Like `coords1DTo3D`, but goes from the vector to the 3D array.  Unlike 
  `coords1DTo3D`, this is an entirely local operation (does not do parallel
  communications)

  TODO: make this do parallel communication if reduce_op is not assignment

  **Inputs**

   * mesh: DG mesh
   * coords_vec: vector of coordinate-like values
   * reduce_op: [`Reduction`](@ref) object, default [`AssignmnetReduction`](@ref)

  **Inputs/Outputs**

   * coords_arr: 3D array to be populated
"""
function coords1DTo3D(mesh::PumiMeshDG, coords_vec::AbstractVector,
                      coords_arr::AbstractArray{T, 3},
                      reduce_op::Reduction{T}=AssignReduction{T}()) where {T}

  @assert mesh.coord_order <= 2
  @assert size(coords_arr, 3) == mesh.numEl
  @assert size(coords_arr, 2) == mesh.coord_numNodesPerElement
  @assert size(coords_arr, 1) == mesh.dim
  @assert length(coords_vec) == mesh.coord_numNodes*mesh.dim


  fill!(coords_arr, reduce_op.neutral_element)
  node_entities = ElementNodeEntities(mesh.m_ptr, mesh.coordshape_ptr, mesh.dim)

  for i=1:mesh.numEl
    getNodeEntities(node_entities, mesh.elements[i])
    for j=1:mesh.coord_numNodesPerElement
      entity = node_entities.entities[j]
      for k=1:mesh.dim
        idx = getNumberJ(mesh.coord_nodenums_Nptr, entity, 0, k-1)

        coords_arr[k, j, i] = reduce_op(coords_arr[k, j, i], coords_vec[idx])
      end
    end
  end

  return nothing
end
