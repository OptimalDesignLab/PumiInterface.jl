# functions that do not call the PUMI shared library, but are related to those
# function


export getAdjacentFull

function getAdjacentFull(m_ptr, entity, dimension::Integer)
# returns an array with the adjacencies of meshentity entity of specified dimension, and the number of entries in the array
# this encasuplates the calls to countAdjacent and getAdjacent

  n = countAdjacent(m_ptr, entity, dimension)
  adj = getAdjacent(dimension)

return adj, n

end
