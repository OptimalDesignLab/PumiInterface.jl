# file for utility functions
function flattenArray{T}(A::AbstractArray{AbstractArray{T}, 1})
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

function getMinandMax{T}(arr::AbstractArray{T})
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

