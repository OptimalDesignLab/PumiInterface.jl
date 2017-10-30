# functions to help with creating submeshes

"""
  Copy the balues from the second array to the first, using the mask to map
  the last dimension.  For example:

  A[:, :, idx] = B[:, :, list[i]] if i != 0

  This function always copies the elements of the arrays. and preserves the
  order of `list`

  **Inputs**

   * src: source array
   * list: list of elements to copy (used for last dimension of array)

  **Inputs/Outputs**

   * dest: destination array
"""
function copy_masked(dest::AbstractVector, src::AbstractVector,
                     list::AbstractVector)

  @assert length(list) == size(dest, 1)

  pos = 1
  for i in list
    dest[pos] = copy(src[i])
    pos += 1
  end
end 

function copy_masked(dest::AbstractMatrix, src::AbstractMatrix,
                     list::AbstractVector)

  @assert length(list) == size(dest, 2)
  @assert size(dest, 1) == size(src, 1)

  pos = 1
  for i in list
    for j=1:size(dest, 1)
      dest[j, pos] = copy(src[j, i])
    end
    pos += 1
  end
end 

function copy_masked(dest::Abstract3DArray, src::Abstract3DArray,
                     list::AbstractVector)

  @assert length(list) == size(dest, 3)
  @assert size(dest, 1) == size(src, 1)
  @assert size(dest, 2) == size(src, 2)

  pos = 1
  for i in list
    for k=1:size(dest, 2)
      for j=1:size(dest, 1)
        dest[j, k, pos] = copy(src[j, k, i])
      end
    end
    pos += 1
  end  # end if
end   # end function

function copy_masked(dest::Abstract4DArray, src::Abstract4DArray,
                     list::AbstractVector)

  @assert length(list) == size(dest, 4)
  @assert size(dest, 1) == size(src, 1)
  @assert size(dest, 2) == size(src, 2)
  @assert size(dest, 3) == size(src, 3)

  pos = 1
  for i in list
    for p=1:size(dest, 3)
      for k=1:size(dest, 2)
        for j=1:size(dest, 1)
          dest[j, k, p, pos] = copy(src[j, k, p, i])
        end
      end
    end
    pos += 1
  end  # end if
end   # end function

"""
  Copies the main data arrays from one mesh to a submesh

  TODO: list all arrays

  **Inputs**

   * oldmesh: the original mesh
   * el_list: list of elements to copy (order is preserved)

  **Inputs/Outputs**

   * mesh: the submesh
"""
function copy_data_arrays(mesh::PumiMesh, oldmesh::PumiMesh, el_list::AbstractVector)

  mesh.vert_coords = zeros(T1, mesh.dim, mesh.numVertPerElement, mesh.numEl)
  copy_masked(mesh.vert_coords, oldmesh.vert_coords, el_list)

  mesh.coords = zeros(T1, mesh.dim, mesh.numNodesPerElement, mesh.numEl)
  copy_masked(mesh.coords, oldmesh.coords, el_list)


  return nothing
end


"""
  This function counts the number of the different kinds of entities and
  gets their pointers on oldmesh.m_ptr.  This function reproduces the behavior
  of [`getEntityPointers`](@ref) for 2d-3d compatability

  **Fields Populated**

   * mesh.verts
   * mesh.edges
   * mesh.faces
   * mesh.elements
   * mesh.numVert
   * mesh.numEdge
   * mesh.numFace
   * mesh.numEl

  The entities in the arrays are in the same order as they are stored in Pumi.
  This is different than [`getEntityPointers`](@ref) which uses a Numbering
  object to assign them indices in the arrays.
"""
function count_entities(mesh::PumiMesh, oldmesh::PumiMesh, el_list::AbstractVector)
# get arrays of pointers to elements on the new mesh
# also count how many there are

  # make a set of el_list for fast test of containment
  el_set = Set(el_list)

  mesh.numEl = length(el_list)

  mesh.verts = verts = Array(Ptr{Void}, 0)
  mesh.edges = edges = Array(Ptr{Void}, 0)

  if mesh.dim == 2
    mesh.faced = edges
    entity_arrays = [verts, edges]
  else
    mesh.faces = faces = Array(Ptr{Void}, 0)
    entity_arrays = [verts, edges, faces]
  end

  # apf::Up
  up_adj = Array(Ptr{Void}, 400)

  # count vertices, edges, and faces with an element as a parent
  for dim=0:(mesh.dim - 1)
    obj_arr = entity_arrays[dim]
    
    itr = MeshIterator(mesh.m_ptr, dim)
    for i=1:oldmesh.numEntitiesPerType[dim+1]
      v_i = iterate(mesh.m_ptr, itr)

      # get elements
      nvals = countAdjacent(mesh.m_ptr, v_i, mesh.dim)
      @assert nvals <= 400
      getAdjacent(up_adj)

      for j=1:nvals
        if up_adj[j] in el_set
          push!(obj_arr, v_i)
          break
        end
      end  # end loop j

    end  # end loop i

    free(itr)
  end  # end loop dim

  # get the elements array
  els = Array(Ptr{Void}, mesh.numEl)
  copy_masked(els, mesh.elements, el_list)

  mesh.numVert = length(mesh.verts)
  mesh.numEdge = length(mesh.edges)
  mesh.numFace = length(mesh.faces)

  return nothing
end






