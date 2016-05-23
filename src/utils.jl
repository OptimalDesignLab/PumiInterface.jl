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


