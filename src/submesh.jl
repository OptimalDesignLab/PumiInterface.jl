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

  mesh.verts = verts = Array{Ptr{Void}}(0)
  mesh.edges = edges = Array{Ptr{Void}}(0)

  if mesh.dim == 2
    mesh.faced = edges
    entity_arrays = [verts, edges]
  else
    mesh.faces = faces = Array{Ptr{Void}}(0)
    entity_arrays = [verts, edges, faces]
  end

  # apf::Up
  up_adj = Array{Ptr{Void}}(400)

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
  els = Array{Ptr{Void}}(mesh.numEl)
  copy_masked(els, mesh.elements, el_list)

  mesh.numVert = length(mesh.verts)
  mesh.numEdge = length(mesh.edges)
  mesh.numFace = length(mesh.faces)

  return nothing
end

"""
  original mesh -> submesh  for vectors (indexed using mesh.dofs)

  **Inputs**

   * oldmesh: the original mesh
   * oldq: the AbstractVector of values to be copied to the new mesh (for
           the elements that exist on the new mesh), length oldmesh.numdof
   * newmesh: the submesh

  **Inputs/Outputs**

   * newq: vector of length newmesh.numdof to put values into
"""
function injectionOperator(oldmesh::PumiMesh, oldq::AbstractVector,
                           newmesh::PumiMesh, newq::AbstractVector)

  @assert oldmesh.m_ptr == getOldMesh(newmesh.subdata)

  parent_N = getParentNumbering(newmesh.subdata)

  for i=1:newmesh.numEl
    el_i = newmesh.elements[i]
    oldel = getNumberJ(parent_N, el_i, 0, 0) + 1

    for j=1:newmesh.numNodesPerElement
      for k=1:newmesh.numDofPerNode
        olddof = oldmesh.dofs[k, j, oldel]
        newdof = newmesh.dofs[k, j, i]
        newq[newdof] = oldq[olddof]
      end
    end
  end

  return nothing
end

"""
  Injection operator for 3D arrays, see other method for details.
"""
function injectionOperator(oldmesh::PumiMesh, oldq::Abstract3DArray,
                           newmesh::PumiMesh, newq::Abstract3DArray)

  @assert oldmesh.m_ptr == getOldMesh(newmesh.subdata)

  parent_N = getParentNumbering(newmesh.subdata)

  for i=1:newmesh.numEl
    el_i = newmesh.elements[i]
    oldel = getNumberJ(parent_N, el_i, 0, 0) + 1

    for k=1:newmesh.numNodesPerElement
      for j=1:newmesh.numDofPerNode
        newq[j, k, i] = oldq[j, k, oldel]
      end
    end
  end

  return nothing
end

"""
  submesh -> original mesh for vectors, inverse function of
  [`injectionOperator`](@ref)

  **Inputs**

   * newmesh: the submesh
   * newq: vector of length newmesh.numdof
   * oldmesh: the original mesh

  **Inputs/Outputs**

   * oldq: vector of length oldmesh.numdof.  Entries that have corresponding
           dofs on the submesh are overwritten, other entries are unchanged.

"""
function rejectionOperator(newmesh::PumiMesh, newq::AbstractVector,
                           oldmesh::PumiMesh, oldq::AbstractVector)

  @assert oldmesh.m_ptr == getOldMesh(newmesh.subdata)

  parent_N = getParentNumbering(newmesh.subdata)

  for i=1:newmesh.numEl
    el_i = newmesh.elements[i]
    oldel = getNumberJ(parent_N, el_i, 0, 0) + 1

    for j=1:newmesh.numNodesPerElement
      for k=1:newmesh.numDofPerNode
        olddof = oldmesh.dofs[k, j, oldel]
        newdof = newmesh.dofs[k, j, i]
        oldq[olddof] = newq[newdof]
      end
    end
  end

  return nothing
end                  

"""
  rejectionOperator() method for 3D arrays, see other method for details
"""
function rejectionOperator(newmesh::PumiMesh, newq::Abstract3DArray,
                           oldmesh::PumiMesh, oldq::Abstract3DArray)

  @assert oldmesh.m_ptr == getOldMesh(newmesh.subdata)

  parent_N = getParentNumbering(newmesh.subdata)

  for i=1:newmesh.numEl
    el_i = newmesh.elements[i]
    oldel = getNumberJ(parent_N, el_i, 0, 0) + 1

    for k=1:newmesh.numNodesPerElement
      for j=1:newmesh.numDofPerNode
        oldq[j, k, oldel] = newq[j, k, i]
      end
    end
  end

  return nothing
end

"""
  This function returns an Array of Interface objects in the same order as
  newmesh.bndryfaces where the element that also exists on the submesh is
  elementL.  This can be useful for interpolating the solution for the new
  boundary condition needed for the submesh.
  This boundary condion must always be the last one.

  Note: 2D only (for now)

  **Inputs**

   * oldmesh: the original mesh
   * newmesh: the submesh

  **Outputs**

   * interface_arr: array of Interface objects
"""
function getBoundaryInterpArray(oldmesh::PumiMesh2D, newmesh::PumiMesh2D)

  rng = newmesh.bndry_offsets[end-1]:(newmesh.bndry_offsets[end]-1)
  interface_arr = Array{Interface}(length(rng))

  parent_N = getParentNumbering(newmesh.subdata)
  edges = Array{Ptr{Void}}(12)  # apf::Downward
  el_arr = Array{Ptr{Void}}(2)  # number of elements that share an edge
  pos = 1  # position in interface_arr

  for i in rng
    el_i = newmesh.bndryfaces[i].element
    face_i = newmesh.bndryfaces[i].face
    el_ptr = newmesh.elements[el_i]

    parent_el = getNumberJ(parent_N, el_ptr, 0, 0) + 1
    parent_elptr = oldmesh.elements[parent_el]

    getDownward(oldmesh.m_ptr, parent_elptr, oldmesh.dim-1, edges)
    edge_ptr = edges[face_i]

    nup = countAdjacent(oldmesh.m_ptr, edge_ptr, oldmesh.dim)
    @assert nup == 2
    getAdjacent(el_arr)

    # figure out which is left and which is right
    elL_ptr = parent_elptr
    if el_arr[1] == elL_ptr
      elR_ptr = el_arr[2]
    else
      elR_ptr = el_arr[1]
    end

   # get local face numbers
   elnumL = getNumberJ(oldmesh.el_Nptr, elL_ptr, 0, 0) + 1
   elnumR = getNumberJ(oldmesh.el_Nptr, elR_ptr, 0, 0) + 1
   edgenum = getNumberJ(oldmesh.face_Nptr, edge_ptr, 0, 0) + 1

   facelocalnumL = getFaceLocalNum(oldmesh, edgenum, elnumL)
   facelocalnumR = getFaceLocalNum(oldmesh, edgenum, elnumR)

   interface_arr[pos] = Interface(elnumL, elnumR, facelocalnumL, facelocalnumR, UInt8(1))  # orientation is always 1 for 2d
   pos += 1
 end

  return interface_arr
end
