# functions for gathering MeshEntity*s 
"""
  This function calculates the numNodesPerType and typeOffsetsPerElement arrays
  and returns them

  Inputs:
    fshape: a FieldShape*
    dim: the dimensionality of the mesh (2 or 3)
    numTypePerElement: number of each dimension entity per element (see
                       interfaces.md)

    Outputs:
      numNodesPerType: array of length 3 or 4 (2d or 3d), containing the number
                       of nodes on each vert, edge, face (or region).

      typeOffsetsPerElement: array of length dim + 2 containing the index of 
                             the first node of each type.  The last element
                             is one more than the number of nodes

"""
function getNodeInfo(fshape::Ptr{Void}, dim::Integer, numTypePerElement::Array{I, 1}) where I <: Integer

  num_nodes_v = apf.countNodesOn(fshape, 0)  # number of nodes on a vertex
  num_nodes_e = apf.countNodesOn(fshape, 1) # on edge
  num_nodes_f = apf.countNodesOn(fshape, 2) # on face

  if dim == 3
    num_nodes_r = apf.countNodesOn(fshape, apf.TET)  # on region
    numNodesPerType = [num_nodes_v, num_nodes_e, num_nodes_f, num_nodes_r]
  else
    numNodesPerType = [num_nodes_v, num_nodes_e, num_nodes_f]
  end
  # count numbers of different things per other thing
  # use for bookkeeping
  typeOffsetsPerElement = zeros(Int, dim+2)
  pos = 1
  typeOffsetsPerElement[1] = pos
  for i=2:(dim + 2)
    pos += numTypePerElement[i-1]*numNodesPerType[i-1]
    typeOffsetsPerElement[i] = pos
  end

  return numNodesPerType, typeOffsetsPerElement
end


# can be generalized trivially
function getBoundaryElements(mesh::PumiMeshDG2, bndry_edges::AbstractArray{Int, 1})
# get the faces corresponding to the boundary edges

bndry_faces = zeros(bndry_edges)
faces = Array{Ptr{Void}}(400)  # equivilent to apf::Up
numbering = mesh.entity_Nptrs[mesh.dim + 1]
for i=1:bndry_edges
  edgenum_i = bndry_edges[i]
  face_i = mesh.faces[edgenum_i]
  numFace = apf.countAdjacent(mesh.m_ptr, face_i, mesh.dim)  # should be count upward
  apf.getAdjacent(faces)
  facenum = apf.getNumberJ(numbering, face_i, 0, 0) + 1
#  facenum = getFaceNumber2(faces[1]) + 1

  bndry_faces[i] = facenum
end

return bndry_faces

end

# add if statements to generalized this to 3D
function getEntityPointers(mesh::PumiMesh)
# get the pointers to all the apf::MeshEntities and put them in arrays
# uses the Numberings to determine what the index in the array of each
# entity


  verts = Array{Ptr{Void}}(mesh.numVert)
  edges = Array{Ptr{Void}}(mesh.numEdge)
  elements = Array{Ptr{Void}}(mesh.numEl)
  entity = Ptr{Void}(0)
  idx = 0
  # get pointers to all MeshEntities
  # also initilize the field to zero
#  resetAllIts2(mesh.m_ptr)
#  comps = zeros(dofpernode)
  it = apf.MeshIterator(mesh.m_ptr, 0)
  for i=1:mesh.numVert
#    entity = getVert()
    entity = apf.iterate(mesh.m_ptr, it)
    idx = apf.getNumberJ(mesh.vert_Nptr, entity, 0, 0) + 1
    verts[idx] = entity
#    incrementVertIt()
  end
  apf.free(mesh.m_ptr, it)

  it = apf.MeshIterator(mesh.m_ptr, 1)
  for i=1:mesh.numEdge
    entity = apf.iterate(mesh.m_ptr, it)
#    entity = getEdge()
    idx = apf.getNumberJ(mesh.edge_Nptr, entity, 0, 0) + 1
    edges[idx] = entity
#    incrementEdgeIt()
  end
  apf.free(mesh.m_ptr, it)

  if mesh.dim == 3
    faces = Array{Ptr{Void}}(mesh.numFace)
    it = apf.MeshIterator(mesh.m_ptr, 2)
    for i=1:mesh.numFace
      entity = apf.iterate(mesh.m_ptr, it)
#      entity = getFace()
      idx = apf.getNumberJ(mesh.face_Nptr, entity, 0, 0) + 1
      faces[idx] = entity
#      incrementFaceIt()
    end
    apf.free(mesh.m_ptr, it)

    it = apf.MeshIterator(mesh.m_ptr, 3)
    for i=1:mesh.numEl
      entity = apf.iterate(mesh.m_ptr, it)
#      entity = getEl()
      idx = apf.getNumberJ(mesh.el_Nptr, entity, 0, 0) + 1
      elements[idx] = entity
#      incrementElIt()
    end
    apf.free(mesh.m_ptr, it)


  else
    faces = edges
    it = apf.MeshIterator(mesh.m_ptr, 2)
    for i=1:mesh.numEl
      entity = apf.iterate(mesh.m_ptr, it)
#      entity = getFace()
      idx = apf.getNumberJ(mesh.el_Nptr, entity, 0, 0) + 1
      elements[idx] = entity
#      incrementFaceIt()
    end
    apf.free(mesh.m_ptr, it)


  end

#  resetAllIts2(mesh.m_ptr)

  return verts, edges, faces, elements

end  # end getEntityPointers

function getNodeIdx(e_type::Integer, e_idx::Integer, node_idx::Integer, typeOffsetsPerElement, numNodesPerType )
# calculate the local node number (in the list of all nodes on the element)
# from the entity the node is classified on and its index on the entity

# e_type is the type of the entity: 1 = vert, 2 = edge, 3 = face
# node_idx is the index of the node on the entity it is classified on
# e_idx is the index of the entity the node is classified on (ie 1st, 2nd or
# 3rd edge)

  type_offset = typeOffsetsPerElement[e_type]
  nodes_on_type = numNodesPerType[e_type]
  return type_offset -1 + (e_idx - 1)*nodes_on_type + node_idx
end

#------------------------------------------------------------------------------
# misc. helper functions
function getGlobalNodeNumbers(mesh::PumiMesh, elnum::Integer; getdofs=true)

  if getdofs
    numDofPerNode = mesh.numDofPerNode
  else
    numDofPerNode = 1
  end
  dofnums = zeros(Int32, numDofPerNode, mesh.numNodesPerElement)
  getGlobalNodeNumbers(mesh, elnum, dofnums)

  return dofnums
end

# should be generalizable with entity counts
function getGlobalNodeNumbers(mesh::PumiMesh, elnum::Integer, dofnums::AbstractArray{Int32}; getdofs=true)
# gets global node numbers of all dof on all nodes of the element
# output formap is array [numdofpernode, nnodes]  (each column contains dof numbers for a node)
# if getdofs=false, then dofnums need only have 1 column, or can be a vector
# dofnums must be passable to C
# 2D only
# getdofs specifies whether to get dof numbers or node numbers
# 

  el_i = mesh.elements[elnum]
  type_i = apf.getType(mesh.m_ptr, el_i)  # what is this used for?

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
  node_entities = apf.getNodeEntities(mesh.m_ptr, mesh.mshape_ptr, el_i)

  # get node offsets in the SBP order
  node_offsets = sview(mesh.elementNodeOffsets, :, elnum)
  #node_offsets = sview(mesh.elementNodeOffsets[:, elnum])

  PumiInterface.apf.getDofNumbers(numbering_ptr, node_entities, node_offsets, mesh.nodemapPumiToSbp, el_i, dofnums)  # C implimentation

  return nothing
end


# is this really needed?
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


if typeof(mesh) <: Union{PumiMeshDG2, PumiMesh2}

  array1 = [mesh.vert_Nptr, mesh.edge_Nptr, mesh.el_Nptr]
  #array2 = [mesh.verts; mesh.edges; mesh.elements]  # array of arrays
  array2 = Array{Array{Ptr{Void},1}}(3)
  array2[1] = mesh.verts
  array2[2] = mesh.edges
  array2[3] = mesh.elements
else # 3d mesh
   array1 = [mesh.vert_Nptr, mesh.edge_Nptr, mesh.face_Nptr,  mesh.el_Nptr]
  #array2 = [mesh.verts; mesh.edges; mesh.elements]  # array of arrays
  array2 = Array{Array{Ptr{Void},1}}(4)
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
entity_type = apf.getType(mesh.m_ptr, entity)
println("entity_type = ", entity_type)
entity_num = apf.getNumberJ(mesh.edge_Nptr, entity, 0, 0)
println("edge number = ", entity_num)
=#

if input_dimension > output_dimension # downward adjacencies
#  println("getting downward adjacencies")
  adjacent_entities, num_adjacent = apf.getDownward(mesh.m_ptr, entity, output_dimension)

else  # upward adjacencies
#  println("getting upward adjacencies")
  num_adjacent = apf.countAdjacent(mesh.m_ptr, entity, output_dimension)
  adjacent_entities = apf.getAdjacent(num_adjacent)
end

# get their numbers
adjacent_nums = zeros(Int, num_adjacent)
for i=1:num_adjacent # loop over adjacent entities
  
  # get their numbers here
  adjacent_nums[i] = apf.getNumberJ(numbering_ptr, adjacent_entities[i], 0, 0) + 1
end

return adjacent_nums, num_adjacent


end

# this can be generalized with mesh.dim
function getBoundaryFaceLocalNum(mesh::PumiMesh, edge_num::Integer)
# gets the local edge number of a specified edge that is on the boundary
# of the mesh
# edge_num is an edge num from the output of getBoundaryEdgeNums() (ie. the global edge number)
# the local edge number is the edge number within the element (1st,s 2nd, 3rd...)
#TODO: usage of this function might be a problem because it returns the local face number
#      of the *Pumi* reference topology, not the SBP topology.
  face_i = mesh.faces[edge_num]

  # get mesh face associated with edge
  apf.countAdjacent(mesh.m_ptr, face_i, mesh.dim)
  el = apf.getAdjacent(1)[1]  # get the single face (not an array)

#  elnum_i = apf.getNumberJ(mesh.el_Nptr, face, 0, 0) + 1
#  facenum_i = getFaceNumber2(face) + 1  # convert to 1 based indexing
  down_faces = Array{Ptr{Void}}(12)
  numedges = apf.getDownward(mesh.m_ptr, el, mesh.dim-1, down_faces)
  facenum_local = 0
  for j = 1:mesh.numFacesPerElement  # find which local edge is edge_i
    if down_faces[j] == face_i
      facenum_local = j
    end
  end

  return facenum_local

end

# this can be generalized with some topology and dimension information
function getFaceLocalNum(mesh::PumiMesh, edge_num::Integer, element_num::Integer)
# find the local edge number of a specified edge on a specified element

  edge = mesh.faces[edge_num]
  element = mesh.elements[element_num]

  down_edges = Array{Ptr{Void}}(12)
   numedges = apf.getDownward(mesh.m_ptr, element, mesh.dim-1, down_edges)

  edgenum_local = 0
   for j = 1:mesh.numFacesPerElement  # find which local edge is edge_i
    if down_edges[j] == edge
      edgenum_local = j
    end
  end

  return edgenum_local

end

function getElementVertMap(mesh::PumiMesh)
# get array that maps from elements to their vertex numbers
# the array is numVertPerElement x numEl
  numVertPerElement = mesh.numTypePerElement[1]

  elvertmap = zeros(Int32, numVertPerElement, mesh.numEl)
  verts = Array{Ptr{Void}}(12)  # apf::Downward
  vert_dim = 0

  for i=1:mesh.numEl
    el_i = mesh.elements[i]
    apf.getDownward(mesh.m_ptr, el_i,  vert_dim, verts)

    for j=1:numVertPerElement
      vert_j = verts[j]
      elvertmap[j, i] = apf.getNumberJ(mesh.vert_Nptr, vert_j, 0, 0) + 1
    end
  end

  return elvertmap
end


