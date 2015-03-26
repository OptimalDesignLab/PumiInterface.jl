# functions that do not call the PUMI shared library, but are related to those
# function


export getAdjacentFull, resetAllIts2, countDownwards, countAllNodes, printEdgeVertNumbers, getValues2, getLocalGradients2, getJacobian2

function getAdjacentFull(m_ptr, entity, dimension::Integer)
# returns an array with the adjacencies of meshentity entity of specified dimension, and the number of entries in the array
# this encasuplates the calls to countAdjacent and getAdjacent

  n = countAdjacent(m_ptr, entity, dimension)
  adj = getAdjacent(n)

return adj, n

end

function resetAllIts2()
 # reset all the iterationr that exist for a 2d mesh
 resetVertIt()
 resetEdgeIt()
 resetFaceIt()

 return nothing
end

function countDownwards(m_ptr, entity)
# gets an array containing the number of downward adjacencies of each type for
# the given entity, in ascending order
# ie. array[1] = # of vertices, array[3] = # of faces

  entity_type = getType(m_ptr, entity)

  if entity_type == 0  # vertex
    return [0]  # return a non empty array

  elseif entity_type == 1  # edge
    return [2]

  elseif entity_type == 2 # triangle
    return [3 3]
  elseif entity_type == 3 # quadralaterial
    return [4, 4]
  elseif entity_type == 4 # tet
    return [4, 6, 4]
  elseif entity_type == 5 # hex (wedge?)
    return [6, 9, 5]
  elseif entity_type == 6 # prism (brick?)
    return [8, 12, 6]
  elseif entity_type == 7 #pyramid
    return [5, 8, 5]
  else
    println("unrecognized entity type")
    return [0]
  end

end

function countAllNodes(mshape_ptr, entity_type::Integer)
# gets the total number of nodes on an entity type, including downward dependencies

  eshape_ptr = getEntityShape(mshape_ptr, entity_type)
  i = countNodes(eshape_ptr)
  return i
end

function printEdgeVertNumbers(edgeN_ptr, vertN_ptr)
# prints the number of the verteces on the edge
# edgeN_ptr is a pointer to the numbering of edges
# vertN_ptr is a pointer to the number of vertices
# numEdges is the number of edges in the mesh
# this function uses iterators

  resetEdgeIt();
  m_ptr = getMesh(edgeN_ptr)  # get pointer to mesh
  m = countJ(m_ptr, 1)  # count number of edges on the mesh
  println("m = ", m)

  for i=1:m  # loop over edges
    edge_i = getEdge()
    edge_num = getNumberJ(edgeN_ptr, edge_i, 0, 0)

    (verts, nverts) = getDownward(m_ptr, edge_i, 0)
    n1 = getNumberJ(vertN_ptr, verts[1], 0, 0)
    n2 = getNumberJ(vertN_ptr, verts[2], 0, 0)
    println("edge ", edge_num, " has vertices ", n1, " , ", n2 )
    incrementEdgeIt()
  end

end


function getValues2(eshape_ptr, coords::Array{Float64,1})
# gets the values of the shape functions at coords
# output is a vector whose length is the number of nodes affecting the entity used to get eshape_ptr

  numN = countNodes(eshape_ptr)
  vals = getValues(eshape_ptr, coords, numN)
  return vals
end

function getLocalGradients2(eshape_ptr, coords)
# get the shape function derivative values at coords
# the output is a 3 x numN matrix, where numN is the number of nodes affecting the entity usued to get eshape_ptr

   numN = countNodes(eshape_ptr)
   vals = getLocalGradients(eshape_ptr, coords, numN)
   return vals
end

function getJacobian2(m_ptr, entity, coords)
# creates a mesh element, then gets its jacobian at the coordinates

   mel_ptr = createMeshElement(m_ptr, entity)
   jac = getJacobian(mel_ptr, coords)
   return jac
end
