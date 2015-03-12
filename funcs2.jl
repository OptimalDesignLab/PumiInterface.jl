#useful functions that do no directly call C code

function getAllVertCoords()
numV = num_entities[1]

vertCoords = zeros(3, numV)  # storage for all vertex coordinates
coords_tmp = zeros(3,1)  # temporary storage for vetex coordinates

for i=1:numV
  getVertCoords(coords_tmp, 3, 1)
  vertCoords[:, i] = coords_tmp
  incrementVertIt()
end

resetVertIt()

return vertCoords
end

function makeSBPMatrix()
# make matrix needed for summation by parts

numEl = num_entities[4]  # number of elements
numPerEl = downward_counts[4, 1]  # nodes per element 

coords_tmp = zeros(3, numPerEl)
sbpMatrix = zeros(3, numPerEl, numEl)

for i = 1:numEl 
    getElCoords( coords_tmp, 3,  numPerEl)
    sbpMatrix[:, :, i] = coords_tmp
    incrementElIt()
end

resetElIt()
return sbpMatrix
end

function incrementIt(dim::Int)
# increment an iterator of a given dimension
if (dim == 0)
  incrementVertIt()
elseif (dim==1)
  incrementEdgeIt()
elseif (dim == 2)
  incrementElIt()
elseif (dim == 3)
  incrementElIt()
end

end  # end function

