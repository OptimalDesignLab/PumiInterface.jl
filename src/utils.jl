# file for utility functions
function flattenArray(A::AbstractArray{AbstractArray{T}, 1}) where T
# copies array-of-array A into a flattened version B

  n = length(A)
  num_entries = 0
  for i=1:n  # count number of entries total
    num_entries += length(A[i])
  end

  B = zeros(T, num_entries)

  B_index = 1
  for i=1:n
    for j = 1:length(A[i])
      B[B_index] = A[i][j]
      B_index += 1
    end
  end

  return B

end


function getElIndex(a::Array, b)

  for i = 1:length(a)
    if a[i] == b
      return i
    end
  end

  return 0
end

function getMinandMax(arr::AbstractArray{T}) where T
# need to check this function for type stability

min_entry = typemax(T)
max_entry = typemin(T)

@inbounds for i=1:length(arr)
  entry_i = arr[i]

  if entry_i < min_entry
    min_entry = entry_i
  end

  if entry_i > max_entry
    max_entry = entry_i
  end
end

return min_entry, max_entry

end

"""
  Count the number of instances of val in arr

  Inputs:
    val: the value to find
    arr: the array to search in

  Outputs:
    cnt: the number of instances of val in arr
"""
function countin(val, arr::AbstractArray)

  cnt = 0
  for i=1:length(arr)
    if arr[i] == val
      cnt += 1
    end
  end

  return cnt
end

function common(A::AbstractArray, B::AbstractArray)
# find common element in the two sets

  # get the larger array
  if length(A) > length(B)
    arr1 = A
    arr2 = B
  else
    arr1 = B
    arr2 = A
  end

  # the output array is no larger than the smaller array
  arr_ret = zeros(eltype(arr2), size(arr2))
  num_found = 0
  for i=1:length(arr1)  # loop over larger array
    # check every entry in arr2 for the current entry in arr1
    el_i = arr1[i]
    for j=1:length(arr2)
      if el_i == arr2[j]
	num_found += 1
	arr_ret[num_found] = el_i
	break
      end
    end
  end

    return arr_ret, num_found
end

"""
  Gets the first element in A that is also in B.  If none found, returns
  first element of A
"""
function first_common(A::AbstractArray, B::AbstractArray)

  for i=1:length(A)
    A_i = A[i]
    for j=1:length(B)
      if A_i == B[j]
        return A_i
      end
    end
  end

  throw(ErrorException("No common element found"))
  return A[1]
end


function calcNewNode(i, offset_pumi, offset_orient)
# this function calculates the new node index on the entity
# using offset_pumi, the offset
# that maps the SBP node index i to the pumi node index, and offset_orient,
# which maps the pumi node index to the pumi node index that takes into account
# the orientation of the mesh entity
# the returned node 
  tmp = abs(offset_pumi - i)
  tmp2 = abs(offset_orient - tmp) - 1
  return tmp2
end

"""
  Calculate the volume of the mesh as the sum of the volume of the elements

  Methods are available for 2D and 3D.
"""
function calcVolume(mesh::PumiMesh2DG)

  area = 0.0
  for i=1:mesh.numEl
    r1x = mesh.vert_coords[1, 2, i] - mesh.vert_coords[1, 1, i]
    r1y = mesh.vert_coords[2, 2, i] - mesh.vert_coords[2, 1, i]

    r2x = mesh.vert_coords[1, 3, i] - mesh.vert_coords[1, 1, i]
    r2y = mesh.vert_coords[2, 3, i] - mesh.vert_coords[2, 1, i]

    area_i = abs(0.5*(r1x*r2y - r1y*r2x))

    area +=  area_i
  end

  area = MPI.Allreduce(area, MPI.SUM, mesh.comm)

  return area
end

function calcVolume(mesh::PumiMesh3DG)

  volume = 0.0
  for i=1:mesh.numEl
    r1x = mesh.vert_coords[1, 2, i] - mesh.vert_coords[1, 1, i]
    r1y = mesh.vert_coords[2, 2, i] - mesh.vert_coords[2, 1, i]
    r1z = mesh.vert_coords[3, 2, i] - mesh.vert_coords[3, 1, i]

    r2x = mesh.vert_coords[1, 3, i] - mesh.vert_coords[1, 1, i]
    r2y = mesh.vert_coords[2, 3, i] - mesh.vert_coords[2, 1, i]
    r2z = mesh.vert_coords[3, 3, i] - mesh.vert_coords[3, 1, i]

    r3x = mesh.vert_coords[1, 4, i] - mesh.vert_coords[1, 1, i]
    r3y = mesh.vert_coords[2, 4, i] - mesh.vert_coords[2, 1, i]
    r3z = mesh.vert_coords[3, 4, i] - mesh.vert_coords[3, 1, i]

    # cross product
    c1 = r2y*r3z - r2z*r3y
    c2 = -(r2x*r3z - r2z*r3x)
    c3 = r2x*r3y - r2y*r3x

    volume_i = (r1x*c1 + r1y*c2 + r1z*c3)/6

    # dot product
    volume += volume_i
  end

  volume = MPI.Allreduce(volume, MPI.SUM, mesh.comm)
  println(mesh.f, "volume = ", volume)

  return volume
end

"""
  Calculate the volume of the mesh by integrating the 1/|dxidx|.
"""
function calcVolumeIntegral(mesh::PumiMesh, sbp)

  volume = 0.0
  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      volume += sbp.w[j]*real(1./mesh.jac[j, i])
    end
  end

  return volume
end


function accumulateAtVerts(mesh::PumiMeshDG, u_volume::Abstract3DArray, u_verts::AbstractMatrix)
# takes u_volume which is n x numVertsPerElement x numEl and does an addition
# reduction to u_verts, which is n x numVerts

  @assert size(u_volume, 1) == size(u_verts, 1)

  n = size(u_volume, 1)
  down_verts = Array{Ptr{Void}}(12)
  
  fill!(u_verts, 0.0)
  for i=1:mesh.numEl
    el_i = mesh.elements[i]
    nverts = apf.getDownward(mesh.m_ptr, el_i, 0, down_verts)

    for j=1:nverts
      vert_j = down_verts[j]
      vertnum_j = apf.getNumberJ(mesh.vert_Nptr, vert_j, 0, 0) + 1

      for k=1:n
        u_verts[k, vertnum_j] +=  u_volume[k, j, i]
      end
    end
  end

  return nothing
end

"""
  This function adds a new boundary condition for the specified geometric
  entity, and removes that geometric entity from any other boundary conditions.

  **Inputs**

   * new_geo: geometric entity number for the BC to be added
   * name: name of new boundary condition

  **Inputs/Outputs**

   * opts: options dictionary

   This function uses the following keys:
   
    * numBC
    * BC*
"""
function updateBCs(opts::Dict, new_geo::Integer, name::AbstractString)

  numBC = opts["numBC"]
  for i=1:numBC
    bndry_geo = opts["BC$i"]
    idx = findfirst(bndry_geo, new_geo)
    if idx != 0
      new_arr = zeros(eltype(bndry_geo), length(bndry_geo)-1)
      new_arr[1:(idx-1)] = bndry_geo[1:(idx-1)]
      new_arr[idx:end] = bndry_geo[(idx+1):end]
      opts["BC$i"] = new_arr
    end
  end

  numBC = numBC + 1
  opts["numBC"]
  opts["numBC"] = numBC
  opts["BC$numBC"] = [new_geo]
  opts["BC$(numBC)_name"] = name

  return nothing
end

