# functions that do not call the PUMI shared library, but are related to those
# function


export getAdjacentFull, resetAllIts2, countDownwards, countAllNodes, printEdgeVertNumbers, printFaceVertNumbers,  getValues2, getLocalGradients2, getJacobian2, getNodeEntities, getEntityLocalNumber, printElementVertNumbers, isShared_sharing

export apfVERTEX, apfEDGE, apfTRIANGLE, apfQUAD, apfTET, apfHEX, apfPRISM, apfPYRAMIX, simplexTypes

# declare the enums
global const apfVERTEX=0
global const apfEDGE=1
global const apfTRIANGLE=2
global const apfQUAD=3
global const apfTET=4
global const apfHEX=5
global const apfPRISM=6
global const apfPYRAMID=7

global const simplexTypes = [apfVERTEX, apfEDGE, apfTRIANGLE, apfTET]

_tmp1, _tmp2, _tmp3 = getTopologyMaps()
global const tri_edge_verts = _tmp1 + 1
global const tet_edge_verts = _tmp2 + 1
global const tet_tri_verts = _tmp3 + 1

function getAdjacentFull(m_ptr, entity, dimension::Integer)
# returns an array with the adjacencies of meshentity entity of specified dimension, and the number of entries in the array
# this encasuplates the calls to countAdjacent and getAdjacent

  n = countAdjacent(m_ptr, entity, dimension)
  adj = getAdjacent(n)

return adj, n

end

function resetAllIts2()
 # reset all the iterationr that exist for a 2d mesh

 m_ptr = getMeshPtr()
 dim = getMeshDimension(m_ptr) 
 reset_funcs = [resetVertIt, resetEdgeIt, resetFaceIt, resetElIt]

 for i=1:(dim+1)
   reset_funcs[i]()
 end


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

function printEdgeVertNumbers(edgeN_ptr, vertN_ptr; fstream=STDOUT)
# prints the number of the verteces on the edge
# edgeN_ptr is a pointer to the numbering of edges
# vertN_ptr is a pointer to the number of vertices
# numEdges is the number of edges in the mesh
# fstream is the stream the data is written to, default STDOUT
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
    println(fstream, "edge ", edge_num, " has vertices ", n1, " , ", n2 )
    incrementEdgeIt()
  end

end

function printFaceVertNumbers(edges::AbstractArray{Ptr{Void}}, edgeN_ptr, vertN_ptr; fstream=STDOUT)

  m_ptr = getMesh(edgeN_ptr)
  m = length(edges)

  for i=1:m
    edge_i = edges[i]
    edge_num = getNumberJ(edgeN_ptr, edge_i, 0, 0)
    print(fstream, "face ", edge_num, " has vertices ")
    (verts, nverts) = getDownward(m_ptr, edge_i, 0)
    for j=1:nverts
      n = getNumberJ(vert_Nptr, verts[j], 0, 0)
      print(fstream, n1, ", ")
    end
    print(fstream, "\n")
  end

end
 
function printFaceVertNumbers(faceN_ptr, vertN_ptr; fstream=STDOUT)
#print the numbers of the vertices that compose a face
resetFaceIt();
m_ptr = getMesh(faceN_ptr)
m = countJ(m_ptr, 2)  # count number of faces on the mesh
println("numFaces = ", m)

 for i=1:m  # loop over faces
   face_i = getFace()
   face_num = getNumberJ(faceN_ptr, face_i, 0, 0)
   (verts, nverts) = getDownward(m_ptr, face_i, 0)
   vertnums = zeros(Int, nverts)
   for j=1:nverts
     vertnums[j] = getNumberJ(vertN_ptr, verts[j], 0, 0)
   end

   println(fstream, "face ", face_num, " has vertices $vertnums")
   incrementFaceIt()
 end

 return nothing
end

function printElementVertNumbers(el_Nptr, vert_Nptr; fstream=STDOUT)
  resetElIt()
  m_ptr = getMesh(el_Nptr)
  m = countJ(m_ptr, 3)  # count number of element

  for i=1:m
    el_i = getEl()
    elnum = getNumberJ(el_NPtr, el_i, 0, 0)
    (verts, nverts) = getDownward(m_ptr, el_i, 0)
    vertnums = zeros(Int, nverts)
    for j=1:nverts
      vertnums[j] = getNumberJ(vert_Nptr, verts[j], 0, 0)
    end
    println(fstream, "element ", elnum, " has vertices $vertnums")
    incrementElIt()
  end

  return nothing
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


function getNodeEntities(m_ptr, mshape_ptr, entity)

  entity_type = getType(m_ptr, entity)
  eshape_ptr = getEntityShape(mshape_ptr, entity_type)
  nnodes = countNodes(eshape_ptr)
  downward_entities = Array(Ptr{Void}, nnodes)  # holds mesh entities
  getNodeEntities(m_ptr, mshape_ptr, entity, downward_entities)

  return downward_entities
end

global const _getNodeEntities_retrieved_entities = Array(Ptr{Void}, 12)  # reusable storage
function getNodeEntities(m_ptr, mshape_ptr, entity, 
                         downward_entities::AbstractArray{Ptr{Void}})
# get the meshentities that have nodes on them 
# duplicates entities that have multiple nodes
# this only works for simplex elements
  
  entity_type = getType(m_ptr, entity)
  eshape_ptr = getEntityShape(mshape_ptr, entity_type)
  nnodes = countNodes(eshape_ptr)
#  downward_entities = Array(Ptr{Void}, nnodes)  # holds mesh entities
  num_vert_nodes = countNodesOn(mshape_ptr, 0)
  num_edge_nodes = countNodesOn(mshape_ptr, 1)
  num_face_nodes = countNodesOn(mshape_ptr, 2)
  num_region_nodes = countNodesOn(mshape_ptr, 4) # tetrahedron
  
  has_vert_nodes = (num_vert_nodes != 0)
  has_edge_nodes = (num_edge_nodes != 0)
  has_face_nodes = (num_face_nodes != 0)
  has_region_nodes = (num_region_nodes != 0)

  if entity_type == 2  # triangle
    num_type_per_entity = [3, 3, 1]
  elseif entity_type == 4  # tetrahedron
    num_type_per_entity = [4, 6, 4, 1]
    numRegions = num_type_per_entity[4]
  else
    throw(ErrorException("unsupported entity type in getNodeEntites, entity type = $entity_type"))
  end

  numVerts = num_type_per_entity[1]
  numEdges = num_type_per_entity[2]
  numFaces = num_type_per_entity[3]
#=

  if entity_type == 2  # triangle
    num_entities = 3*has_vert_nodes + 3*has_edge_nodes + has_face_nodes
  elseif entity_type == 4 # tetrahedron
    num_entities = 4*has_vert_nodes + 6*has_edge_nodes + 4*has_face_nodes + has_region_nodes
  else
    println("element type not supported by getNodeEntities")
  end
  
=#
  retrieved_entities = _getNodeEntities_retrieved_entities
  vert_offset = 0
  if has_vert_nodes  # vertices
    vert_offset = getDownward(m_ptr, entity, 0, downward_entities)
  end
    
	
  if has_edge_nodes # edges
    numEdges = getDownward(mshape_ptr, entity, 1, retrieved_entities)
        for i=1:numEdges
          insertN(downward_entities, retrieved_entities[i], vert_offset+num_edge_nodes*(i-1) + 1, num_edge_nodes)
        end
  end

  if has_face_nodes
     numFaces = getDownward(mshape_ptr, entity, 2, retrieved_entities)
    for i=1:numFaces
      insertN(downward_entities, retrieved_entities[i], vert_offset + num_edge_nodes*numEdges +num_edge_nodes*(i-1) + 1, num_face_nodes)
    end
  end

  if has_region_nodes
    for i=1:numRegions
      insertN(downward_entities, entity, vert_offset + num_edge_nodes*numEdges + num_face_nodes*numFaces + num_region_nodes*(i-1) + 1, num_region_nodes)
    end
  end

  return nothing
#  return downward_entities
end

function insertN{T}(vec::AbstractArray{T}, element::T,  index::Integer, n::Integer)
# insert element into the vector n times, starting at index

  for i=index:(index+n-1)
    vec[i] = element
  end

  return nothing
end

function getEntityLocalNumber(mesh::Ptr{Void}, entity::Ptr{Void}, parent::Ptr{Void},  entity_dim::Integer, parent_dim::Integer)
# get the local index of the entity among all the entities of the same dimension that
# belong to the same parent
# entity_dim is the dimension of the entity
# parent_dim is the dimension of the parent
# zero based indexing

  if parent_dim <= entity_dim
    println("parent dimension must be > entity dimension")
    return -1
  end

  (down_entities, numdown) = getDownward(mesh, parent, entity_dim)

  local_num = -1
  for i=1:numdown
    if (down_entities[i] == entity)
      local_num = i-1
    end
  end

  return local_num
end

# check if an entity is shared according to the supplied sharing
# this is a bit expensive because it has to actually get the 
# shared entities to find out how many there are
function isShared_sharing(shr::Ptr{Void},  entity::Ptr{Void})
  return  countCopies(shr, entity) != 0
end
