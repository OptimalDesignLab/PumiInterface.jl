# functions for gathering MeshEntity*s 

function getEdgeFaces(mesh::PumiMeshDG2, bndry_edges::AbstractArray{Int, 1})
# get the faces corresponding to the boundary edges

bndry_faces = zeros(bndry_edges)
faces = Array(Ptr{Void}, 400)  # equivilent to apf::Up
for i=1:bndry_edges
  edgenum_i = bndry_edges[i]
  edge_i = mesh.edges[edgenum_i]
  numFace = countAdjacent(mesh.m_ptr, edge_i, 2)  # should be count upward
  getAdjacent(faces)
  facenum = getFaceNumber2(faces[1]) + 1

  bndry_faces[i] = facenum
end

return bndry_faces

end

function getEdgeFaces(mesh::PumiMesh2, bndry_edges::AbstractArray{Int, 1})
# get the faces corresponding to the boundary edges

bndry_faces = zeros(bndry_edges)
faces = Array(Ptr{Void}, 400)  # equivilent to apf::Up
for i=1:bndry_edges
  edgenum_i = bndry_edges[i]
  edge_i = mesh.edges[edgenum_i]
  numFace = countAdjacent(mesh.m_ptr, edge_i, 2)  # should be count upward
  getAdjacent(faces)
  facenum = getFaceNumber2(faces[1]) + 1

  bndry_faces[i] = facenum
end

return bndry_faces

end


function getEntityPointers(mesh::PumiMesh)
# get the pointers to all the apf::MeshEntities and put them in arrays
# uses the Numberings to determine what the index in the array of each
# entity


  verts = Array(Ptr{Void}, mesh.numVert)
  edges = Array(Ptr{Void}, mesh.numEdge)
  elements = Array(Ptr{Void}, mesh.numEl)
  entity = Ptr{Void}(0)
  idx = 0
  # get pointers to all MeshEntities
  # also initilize the field to zero
  resetAllIts2()
#  comps = zeros(dofpernode)
  for i=1:mesh.numVert
    entity = getVert()
    idx = getNumberJ(mesh.vert_Nptr, entity, 0, 0) + 1
    verts[idx] = entity
    incrementVertIt()
  end

  for i=1:mesh.numEdge
    entity = getEdge()
    idx = getNumberJ(mesh.edge_Nptr, entity, 0, 0) + 1
    edges[idx] = entity
    incrementEdgeIt()
  end

  for i=1:mesh.numEl
    entity = getFace()
    idx = getNumberJ(mesh.el_Nptr, entity, 0, 0) + 1
    elements[idx] = entity
    incrementFaceIt()
  end

  resetAllIts2()

  return verts, edges, elements

end  # end getEntityPointers

function getEntityPointers(mesh::PumiMesh)
# get the pointers to all the apf::MeshEntities and put them in arrays
# uses the Numberings to determine what the index in the array of each
# entity


  verts = Array(Ptr{Void}, mesh.numVert)
  edges = Array(Ptr{Void}, mesh.numEdge)
  elements = Array(Ptr{Void}, mesh.numEl)
  entity = Ptr{Void}(0)
  idx = 0
  # get pointers to all MeshEntities
  # also initilize the field to zero
  resetAllIts2()
#  comps = zeros(dofpernode)
  for i=1:mesh.numVert
    entity = getVert()
    idx = getNumberJ(mesh.vert_Nptr, entity, 0, 0) + 1
    verts[idx] = entity
    incrementVertIt()
  end

  for i=1:mesh.numEdge
    entity = getEdge()
    idx = getNumberJ(mesh.edge_Nptr, entity, 0, 0) + 1
    edges[idx] = entity
    incrementEdgeIt()
  end

  for i=1:mesh.numEl
    entity = getFace()
    idx = getNumberJ(mesh.el_Nptr, entity, 0, 0) + 1
    elements[idx] = entity
    incrementFaceIt()
  end

  resetAllIts2()

  return verts, edges, elements

end  # end getEntityPointers


function getNodeIdx(e_type::Integer, e_idx::Integer, node_idx::Integer, typeOffsetsPerElement, numNodesPerType )
# calculate the local node number (in the list of all nodes on the element)
# from the entity the node is classified on and its index on the entity

# e_type is the type of the entity: 1 = vert, 2 = edge, 3 = face
# node_idx is the index of the node on the entity it is classified on
# e_idx is the index of the entity the node is classified on (ie 1st, 2nd or
# 3rd edge)

  type_offset = typeOffsetsPerElement[e_type]
  println("typeoffsetsPerElement = ", typeOffsetsPerElement)
  nodes_on_type = numNodesPerType[e_type]
  return type_offset -1 + (e_idx - 1)*nodes_on_type + node_idx
end

#TODO: stop using slice notation
function getCoordinates(mesh::PumiMeshDG2, sbp::AbstractSBP)
# populate the coords array of the mesh object

mesh.coords = Array(Float64, 2, sbp.numnodes, mesh.numEl)

#println("entered getCoordinates")

coords_i = zeros(3,3)
coords_it = zeros(3,2)
for i=1:mesh.numEl  # loop over elements
  
  el_i = mesh.elements[i]
  (sizex, sizey) = size(coords_i)
  getFaceCoords(el_i, coords_i, sizex, sizey)  # populate coords

#  println("coords_i = ", coords_i)

  coords_it[:,:] = coords_i[1:2, :].'
#  println("coords_it = ", coords_it)
  mesh.coords[:, :, i] = SummationByParts.SymCubatures.calcnodes(sbp.cub, coords_it)
#  println("mesh.coords[:,:,i] = ", mesh.coords[:,:,i])
end

return nothing

end


function getBndryCoordinates{Tmsh}(mesh::PumiMeshDG2{Tmsh}, 
                             bndryfaces::Array{Boundary}, 
                             coords_bndry::Array{Tmsh, 3})
# calculate the coordinates on the boundary for the specified faces
# and store in coords_bndry

#  println("----- Entered getBndryCoordinates -----")
  sbpface = mesh.sbpface

  coords_i = zeros(3, 3)
  coords_it = zeros(3, 2)
  coords_edge = zeros(2, 2)
#  facemap = [1 2 1; 2 3 3]

  #TODO: undo once SBP node ordering convention is decided
  facemap = [1 2 3; 2 3 1]
  for i=1:length(bndryfaces)
    bndry_i = bndryfaces[i]

    el = bndry_i.element
    el_ptr = mesh.elements[el]
    face = bndry_i.face

    sizex, sizey = size(coords_i)
    getFaceCoords(el_ptr, coords_i, sizex, sizey)

    coords_it[:, :] = coords_i[1:2, :].'

    # extract the needed vertex coords
#    v1 = facemap[1, face]
#    v2 = facemap[2, face]
    v1 = face
    v2 = mod(face,3) + 1
    coords_edge[1, 1] = coords_it[v1, 1]
    coords_edge[1, 2] = coords_it[v1, 2]
    coords_edge[2, 1] = coords_it[v2, 1]
    coords_edge[2, 2] = coords_it[v2 ,2]

    coords_bndry[:, :, i] = SummationByParts.SymCubatures.calcnodes(sbpface.cub, coords_edge)

  end

end


function getElementVertCoords(mesh::PumiMeshDG2, elnum::Integer, coords::AbstractArray{Float64,2})
# get the coordinates of the vertices of an element
# elnum is the number of the element
# coords is the array to be populated with the coordinates
# each column of coords contains the coordinates for a vertex
# coords must be 3x3

  el_j = mesh.elements[elnum]
  (sizex, sizey) = size(coords)
  getFaceCoords(el_j, coords, sizex, sizey)  # populate coords

  return nothing

end # end function

function getElementVertCoords(mesh::PumiMeshDG2,  elnums::Array{Int,1})
# elnums = vector of element numbbers
# return array of size 3x3xn, where each column (first index) contains the coordinates of a vertex
# n is the number of elements in elnums


# get the number of verts 
n = length(elnums)

coords = zeros(3, 3, n)  # first dimension = x,y, or z, second dimension = vertex number

for j = 1:n  # loop over elements in elnums
#  println("j = ", j)
  elnum_j = elnums[j]
#  println("elnum_j = ", elnum_j)
  el_j = mesh.elements[elnum_j]
  sub_array = sub(coords, :, :, j)
  getFaceCoords(el_j, sub_array, 3, 3)  # populate coords
end

return coords

end

function getGlobalNodeNumbers(mesh::PumiMeshDG2, elnum::Integer; getdofs=true)

  if getdofs
    numDofPerNode = mesh.numDofPerNode
  else
    numDofPerNode = 1
  end
  dofnums = zeros(Int32, numDofPerNode, mesh.numNodesPerElement)
  getGlobalNodeNumbers(mesh, elnum, dofnums)

  return dofnums
end

function getGlobalNodeNumbers(mesh::PumiMeshDG2, elnum::Integer, dofnums::AbstractArray{Int32}; getdofs=true)
# gets global node numbers of all dof on all nodes of the element
# output formap is array [numdofpernode, nnodes]  (each column contains dof numbers for a node)
# if getdofs=false, then dofnums need only have 1 column, or can be a vector
# dofnums must be passable to C
# 2D only
# getdofs specifies whether to get dof numbers or node numbers
# 

el_i = mesh.elements[elnum]
type_i = getType(mesh.m_ptr, el_i)  # what is this used for?

#println("elnum = ", elnum)
# calculate total number of nodes
#nnodes = 3 + 3*mesh.numNodesPerType[2] + mesh.numNodesPerType[3]
nnodes = mesh.numNodesPerElement

if getdofs
  numdof = nnodes*mesh.numDofPerNode  # what is this used for?
  numbering_ptr = mesh.dofnums_Nptr
else
  numdof = nnodes
  numbering_ptr = mesh.nodenums_Nptr
end

# get entities in the Pumi order
node_entities = getNodeEntities(mesh.m_ptr, mesh.mshape_ptr, el_i)

# get node offsets in the SBP order
node_offsets = view(mesh.elementNodeOffsets[:, elnum])
#println("node_offsets = ", [Int(i) for i in node_offsets])
#println("size(node_offsets) = ", size(node_offsets))
#println("node_entities = ", node_entities)
#println("size(node_entities) = ", size(node_entities))
#dofnums = zeros(Int,mesh.numDofPerNode, nnodes)  # to be populated with global dof numbers
#println("nnodes = ", nnodes)
#println("dofpernode = ", mesh.numDofPerNode)



num_entities = [3, 3, 1] # number of vertices, edges, faces

#=
col = 1 # current node
for i=1:3  # loop over verts, edges, faces
  for j=1:num_entities[i]  # loop over all entities of this type
    for k=1:mesh.numNodesPerType[i]  # loop over nodes on this entity
      entity = node_entities[col]  # current entity
      for p=1:mesh.numDofPerNode  # loop over all dofs
	dofnums[p, col] = getNumberJ(mesh.dofnums_Nptr, entity, k-1, p-1)
      end
      col += 1
    end  # end loop over nodes on curren entity
  end  # end loop over entities of current type
end  # end loop over entity types
=#

PumiInterface.getDofNumbers(numbering_ptr, node_entities, node_offsets, mesh.nodemapPumiToSbp, el_i, dofnums)  # C implimentation


#=

for i=1:nnodes
  for j=1:mesh.numDofPerNode
    dofnums[j, i] = getNumberJ(mesh.dofnums_Nptr, node_entities[i], 0, j-1)
  end
end
=#

return nothing


end

function getAdjacentEntityNums(mesh::PumiMesh, entity_index::Integer, input_dimension::Integer, output_dimension::Integer)
# gets the numbers of the adjacent entities (upward or downward) of the given entity
# entity_index specifies the index of the input entity
# input_dimension specifies the dimension of the input entity (0 =vert, 1 = edge ...)
# output_dimension is the dimension of the adjacent  entities to be retrieved
# returns an array containing the indicies of the adjacent entities, and the number of them
# the length of the array might be greater than the number of adjacencies, so use the second returned value for the number of adjacencies


if (input_dimension == output_dimension)
  println("cannot get same level adjacencies")
end


if typeof(mesh) <: PumiMeshDG2

  array1 = [mesh.vert_Nptr, mesh.edge_Nptr, mesh.el_Nptr]
  #array2 = [mesh.verts; mesh.edges; mesh.elements]  # array of arrays
  array2 = Array(Array{Ptr{Void},1}, 3)
  array2[1] = mesh.verts
  array2[2] = mesh.edges
  array2[3] = mesh.elements
else # 3d mesh
   array1 = [mesh.vert_Nptr, mesh.edge_Nptr, mesh.face_Nptr,  mesh.el_Nptr]
  #array2 = [mesh.verts; mesh.edges; mesh.elements]  # array of arrays
  array2 = Array(Array{Ptr{Void},1}, 4)
  array2[1] = mesh.verts
  array2[2] = mesh.edges
  array2[3] = mesh.faces
  array2[4] = mesh.elements
end

# choose which numbering to use
numbering_ptr = array1[output_dimension + 1]
#println("numbering_ptr = ", numbering_ptr)


# get the the entity
#println("array2 = ", array2)
entity_array = array2[input_dimension + 1]
entity = entity_array[entity_index]
#entity = array2[input_dimension+1][entity_index]
#println("entity = ", entity)

#=
entity_type = getType(mesh.m_ptr, entity)
println("entity_type = ", entity_type)
entity_num = getNumberJ(mesh.edge_Nptr, entity, 0, 0)
println("edge number = ", entity_num)
=#

if input_dimension > output_dimension # downward adjacencies
#  println("getting downward adjacencies")
  adjacent_entities, num_adjacent = getDownward(mesh.m_ptr, entity, output_dimension)

else  # upward adjacencies
#  println("getting upward adjacencies")
  num_adjacent = countAdjacent(mesh.m_ptr, entity, output_dimension)
  adjacent_entities = getAdjacent(num_adjacent)
end

# get their numbers
adjacent_nums = zeros(Int, num_adjacent)
for i=1:num_adjacent # loop over adjacent entities
  
  # get their numbers here
  adjacent_nums[i] = getNumberJ(numbering_ptr, adjacent_entities[i], 0, 0) + 1
end

return adjacent_nums, num_adjacent


end



function getBoundaryEdgeLocalNum(mesh::PumiMeshDG2, edge_num::Integer)
# gets the local edge number of a specified edge that is on the boundary
# of the mesh
# edge_num is an edge num from the output of getBoundaryEdgeNums() (ie. the global edge number)
# the local edge number is the edge number within the element (1st,s 2nd, 3rd...)
  edge_i = mesh.edges[edge_num]

  # get mesh face associated with edge
  countAdjacent(mesh.m_ptr, edge_i, 2)
  face = getAdjacent(1)[1]  # get the single face (not an array)

  facenum_i = getFaceNumber2(face) + 1  # convert to 1 based indexing
  down_edges = Array(Ptr{Void}, 3)
  numedges = getDownward(mesh.m_ptr, face, 1, down_edges)
  edgenum_local = 0
  for j = 1:3  # find which local edge is edge_i
    if down_edges[j] == edge_i
      edgenum_local = j
    end
  end

  return edgenum_local

end

function getEdgeLocalNum(mesh::PumiMeshDG2, edge_num::Integer, element_num::Integer)
# find the local edge number of a specified edge on a specified element

  edge = mesh.edges[edge_num]
  element = mesh.elements[element_num]

  down_edges = Array(Ptr{Void}, 3)
   numedges = getDownward(mesh.m_ptr, element, 1, down_edges)

  edgenum_local = 0
   for j = 1:3  # find which local edge is edge_i
    if down_edges[j] == edge
      edgenum_local = j
    end
  end

  return edgenum_local

end




#TODO: stop using slice notation
function getCoordinates(mesh::PumiMesh2, sbp::AbstractSBP)
# populate the coords array of the mesh object

mesh.coords = Array(Float64, 2, sbp.numnodes, mesh.numEl)

#println("entered getCoordinates")

coords_i = zeros(3,3)
coords_it = zeros(3,2)
for i=1:mesh.numEl  # loop over elements
  
  el_i = mesh.elements[i]
  (sizex, sizey) = size(coords_i)
  getFaceCoords(el_i, coords_i, sizex, sizey)  # populate coords

#  println("coords_i = ", coords_i)

  coords_it[:,:] = coords_i[1:2, :].'
#  println("coords_it = ", coords_it)
  mesh.coords[:, :, i] = calcnodes(sbp, coords_it)
#  println("mesh.coords[:,:,i] = ", mesh.coords[:,:,i])
end

return nothing

end




function getElementVertCoords(mesh::PumiMesh2, elnum::Integer, coords::AbstractArray{Float64,2})
# get the coordinates of the vertices of an element
# elnum is the number of the element
# coords is the array to be populated with the coordinates
# each column of coords contains the coordinates for a vertex
# coords must be 3x3

  el_j = mesh.elements[elnum]
  (sizex, sizey) = size(coords)
  getFaceCoords(el_j, coords, sizex, sizey)  # populate coords

  return nothing

end # end function

function getElementVertCoords(mesh::PumiMesh2,  elnums::Array{Int,1})
# elnums = vector of element numbbers
# return array of size 3x3xn, where each column (first index) contains the coordinates of a vertex
# n is the number of elements in elnums


# get the number of verts 
n = length(elnums)

coords = zeros(3, 3, n)  # first dimension = x,y, or z, second dimension = vertex number

for j = 1:n  # loop over elements in elnums
#  println("j = ", j)
  elnum_j = elnums[j]
#  println("elnum_j = ", elnum_j)
  el_j = mesh.elements[elnum_j]
  sub_array = sub(coords, :, :, j)
  getFaceCoords(el_j, sub_array, 3, 3)  # populate coords
end

return coords

end

function getGlobalNodeNumbers(mesh::PumiMesh2, elnum::Integer; getdofs=true)

  if getdofs
    numDofPerNode = mesh.numDofPerNode
  else
    numDofPerNode = 1
  end
  dofnums = zeros(Int32, numDofPerNode, mesh.numNodesPerElement)
  getGlobalNodeNumbers(mesh, elnum, dofnums)

  return dofnums
end

function getGlobalNodeNumbers(mesh::PumiMesh2, elnum::Integer, dofnums::AbstractArray{Int32}; getdofs=true)
# gets global node numbers of all dof on all nodes of the element
# output formap is array [numdofpernode, nnodes]  (each column contains dof numbers for a node)
# if getdofs=false, then dofnums need only have 1 column, or can be a vector
# dofnums must be passable to C
# 2D only
# getdofs specifies whether to get dof numbers or node numbers
# 

el_i = mesh.elements[elnum]
type_i = getType(mesh.m_ptr, el_i)  # what is this used for?

#println("elnum = ", elnum)
# calculate total number of nodes
#nnodes = 3 + 3*mesh.numNodesPerType[2] + mesh.numNodesPerType[3]
nnodes = mesh.numNodesPerElement

if getdofs
  numdof = nnodes*mesh.numDofPerNode  # what is this used for?
  numbering_ptr = mesh.dofnums_Nptr
else
  numdof = nnodes
  numbering_ptr = mesh.nodenums_Nptr
end

# get entities in the Pumi order
node_entities = getNodeEntities(mesh.m_ptr, mesh.mshape_ptr, el_i)
# get node offsets in the SBP order
node_offsets = view(mesh.elementNodeOffsets[:, elnum])
#println("node_offsets = ", [Int(i) for i in node_offsets])
#println("size(node_offsets) = ", size(node_offsets))
#println("node_entities = ", node_entities)
#println("size(node_entities) = ", size(node_entities))
#dofnums = zeros(Int,mesh.numDofPerNode, nnodes)  # to be populated with global dof numbers
#println("nnodes = ", nnodes)
#println("dofpernode = ", mesh.numDofPerNode)



num_entities = [3, 3, 1] # number of vertices, edges, faces

#=
col = 1 # current node
for i=1:3  # loop over verts, edges, faces
  for j=1:num_entities[i]  # loop over all entities of this type
    for k=1:mesh.numNodesPerType[i]  # loop over nodes on this entity
      entity = node_entities[col]  # current entity
      for p=1:mesh.numDofPerNode  # loop over all dofs
	dofnums[p, col] = getNumberJ(mesh.dofnums_Nptr, entity, k-1, p-1)
      end
      col += 1
    end  # end loop over nodes on curren entity
  end  # end loop over entities of current type
end  # end loop over entity types
=#

PumiInterface.getDofNumbers(numbering_ptr, node_entities, node_offsets, mesh.nodemapPumiToSbp, el_i, dofnums)  # C implimentation


#=

for i=1:nnodes
  for j=1:mesh.numDofPerNode
    dofnums[j, i] = getNumberJ(mesh.dofnums_Nptr, node_entities[i], 0, j-1)
  end
end
=#

return nothing


end

function getAdjacentEntityNums(mesh::PumiMesh, entity_index::Integer, input_dimension::Integer, output_dimension::Integer)
# gets the numbers of the adjacent entities (upward or downward) of the given entity
# entity_index specifies the index of the input entity
# input_dimension specifies the dimension of the input entity (0 =vert, 1 = edge ...)
# output_dimension is the dimension of the adjacent  entities to be retrieved
# returns an array containing the indicies of the adjacent entities, and the number of them
# the length of the array might be greater than the number of adjacencies, so use the second returned value for the number of adjacencies


if (input_dimension == output_dimension)
  println("cannot get same level adjacencies")
end


if mesh.dim == 2

  array1 = [mesh.vert_Nptr, mesh.edge_Nptr, mesh.el_Nptr]
  #array2 = [mesh.verts; mesh.edges; mesh.elements]  # array of arrays
  array2 = Array(Array{Ptr{Void},1}, 3)
  array2[1] = mesh.verts
  array2[2] = mesh.edges
  array2[3] = mesh.elements
else # 3d mesh
   array1 = [mesh.vert_Nptr, mesh.edge_Nptr, mesh.face_Nptr,  mesh.el_Nptr]
  #array2 = [mesh.verts; mesh.edges; mesh.elements]  # array of arrays
  array2 = Array(Array{Ptr{Void},1}, 4)
  array2[1] = mesh.verts
  array2[2] = mesh.edges
  array2[3] = mesh.faces
  array2[4] = mesh.elements
end

# choose which numbering to use
numbering_ptr = array1[output_dimension + 1]
#println("numbering_ptr = ", numbering_ptr)


# get the the entity
#println("array2 = ", array2)
entity_array = array2[input_dimension + 1]
entity = entity_array[entity_index]
#entity = array2[input_dimension+1][entity_index]
#println("entity = ", entity)

#=
entity_type = getType(mesh.m_ptr, entity)
println("entity_type = ", entity_type)
entity_num = getNumberJ(mesh.edge_Nptr, entity, 0, 0)
println("edge number = ", entity_num)
=#

if input_dimension > output_dimension # downward adjacencies
#  println("getting downward adjacencies")
  adjacent_entities, num_adjacent = getDownward(mesh.m_ptr, entity, output_dimension)

else  # upward adjacencies
#  println("getting upward adjacencies")
  num_adjacent = countAdjacent(mesh.m_ptr, entity, output_dimension)
  adjacent_entities = getAdjacent(num_adjacent)
end

# get their numbers
adjacent_nums = zeros(Int, num_adjacent)
for i=1:num_adjacent # loop over adjacent entities
  
  # get their numbers here
  adjacent_nums[i] = getNumberJ(numbering_ptr, adjacent_entities[i], 0, 0) + 1
end

return adjacent_nums, num_adjacent


end



function getBoundaryEdgeLocalNum(mesh::PumiMesh2, edge_num::Integer)
# gets the local edge number of a specified edge that is on the boundary
# of the mesh
# edge_num is an edge num from the output of getBoundaryEdgeNums() (ie. the global edge number)
# the local edge number is the edge number within the element (1st,s 2nd, 3rd...)
  edge_i = mesh.edges[edge_num]

  # get mesh face associated with edge
  countAdjacent(mesh.m_ptr, edge_i, 2)
  face = getAdjacent(1)[1]  # get the single face (not an array)

  facenum_i = getFaceNumber2(face) + 1  # convert to 1 based indexing
  down_edges = Array(Ptr{Void}, 3)
  numedges = getDownward(mesh.m_ptr, face, 1, down_edges)
  edgenum_local = 0
  for j = 1:3  # find which local edge is edge_i
    if down_edges[j] == edge_i
      edgenum_local = j
    end
  end

  return edgenum_local

end

function getEdgeLocalNum(mesh::PumiMesh2, edge_num::Integer, element_num::Integer)
# find the local edge number of a specified edge on a specified element

  edge = mesh.edges[edge_num]
  element = mesh.elements[element_num]

  down_edges = Array(Ptr{Void}, 3)
   numedges = getDownward(mesh.m_ptr, element, 1, down_edges)

  edgenum_local = 0
   for j = 1:3  # find which local edge is edge_i
    if down_edges[j] == edge
      edgenum_local = j
    end
  end

  return edgenum_local

end


