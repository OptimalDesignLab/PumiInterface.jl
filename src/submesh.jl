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


