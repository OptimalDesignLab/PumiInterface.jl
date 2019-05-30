# functions that do not call the PUMI shared library, but are related to those
# function


export getAdjacentFull, resetAllIts2, countDownwards, countAllNodes, printEdgeVertNumbers, printFaceVertNumbers,  getValues2, getLocalGradients2, getJacobian2, getNodeEntities, getEntityLocalNumber, printElementVertNumbers, isShared_sharing

export VERTEX, EDGE, TRIANGLE, QUAD, TET, HEX, PRISM, PYRAMIX, simplexTypes, getTypeDimension, getDimension

export ElementNodeEntities

# declare the enums
global const VERTEX=0
global const EDGE=1
global const TRIANGLE=2
global const QUAD=3
global const TET=4
global const HEX=5
global const PRISM=6
global const PYRAMID=7

global const simplexTypes = [VERTEX, EDGE, TRIANGLE, TET]

_tmp1, _tmp2, _tmp3 = getTopologyMaps()
global const tri_edge_verts = _tmp1 + 1
global const tet_edge_verts = _tmp2 + 1
global const tet_tri_verts = _tmp3 + 1

# this is not exported, use the getter function
global const typeDimension = Array{Int}(8)
typeDimension[VERTEX + 1]   = 0
typeDimension[EDGE + 1]     = 1
typeDimension[TRIANGLE + 1] = 2
typeDimension[QUAD + 1]     = 2
typeDimension[TET + 1]      = 3
typeDimension[HEX + 1]      = 3
typeDimension[PRISM + 1]    = 3
typeDimension[PYRAMID + 1]  = 3

"""
  Returns the dimension of a given type (one of the apf::Type enums)
"""
function getTypeDimension(typ::Integer)

  return typeDimension[typ + 1]
end

function getDimension(m_ptr::Ptr{Void}, entity::Ptr{Void})

  typ = getType(m_ptr, entity)
  return getTypeDimension(typ)
end

function getAdjacentFull(m_ptr, entity, dimension::Integer)
# returns an array with the adjacencies of meshentity entity of specified dimension, and the number of entries in the array
# this encasuplates the calls to countAdjacent and getAdjacent

  n = countAdjacent(m_ptr, entity, dimension)
  adj = getAdjacent(n)

  return adj, n

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

  m_ptr = getNumberingMesh(edgeN_ptr)  # get pointer to mesh
  it = MeshIterator(m_ptr, 1)
  m = countJ(m_ptr, 1)  # count number of edges on the mesh
  println("m = ", m)

  for i=1:m  # loop over edges
    edge_i = iterate(it)
    edge_num = getNumberJ(edgeN_ptr, edge_i, 0, 0)

    (verts, nverts) = getDownward(m_ptr, edge_i, 0)
    n1 = getNumberJ(vertN_ptr, verts[1], 0, 0)
    n2 = getNumberJ(vertN_ptr, verts[2], 0, 0)
    println(fstream, n1, " ", n2)
  end

  free(it)

end

function printFaceVertNumbers(edges::AbstractArray{Ptr{Void}}, edgeN_ptr, vertN_ptr; fstream=STDOUT)

  m_ptr = getNumberingMesh(edgeN_ptr)
  m = length(edges)

  for i=1:m
    edge_i = edges[i]
    edge_num = getNumberJ(edgeN_ptr, edge_i, 0, 0)
#    print(fstream, "face ", edge_num, " has vertices ")
    (verts, nverts) = getDownward(m_ptr, edge_i, 0)
    for j=1:nverts
      n = getNumberJ(vert_Nptr, verts[j], 0, 0)
      print(fstream, n1, " ")
    end
    print(fstream, "\n")
  end

end
 
function printFaceVertNumbers(faceN_ptr, vertN_ptr; fstream=STDOUT)
#print the numbers of the vertices that compose a face
  m_ptr = getNumberingMesh(faceN_ptr)
  it = MeshIterator(m_ptr, 2)
  m = countJ(m_ptr, 2)  # count number of faces on the mesh

 for i=1:m  # loop over faces
   face_i = iterate(it)
   face_num = getNumberJ(faceN_ptr, face_i, 0, 0)
   (verts, nverts) = getDownward(m_ptr, face_i, 0)
   vertnums = zeros(Int, nverts)
   for j=1:nverts
     vertnums[j] = getNumberJ(vertN_ptr, verts[j], 0, 0)
     print(fstream, vertnums[j], " ")
   end
   print(fstream, "\n")
 end

 free(it)

 return nothing
end

function printElementVertNumbers(el_Nptr, vert_Nptr; fstream=STDOUT)
  m_ptr = getNumberingMesh(el_Nptr)
  it = MeshIterator(m_ptr, 3)
  m = countJ(m_ptr, 3)  # count number of element

  for i=1:m
    el_i = iterate(it)
    elnum = getNumberJ(el_NPtr, el_i, 0, 0)
    (verts, nverts) = getDownward(m_ptr, el_i, 0)
    vertnums = zeros(Int, nverts)
    for j=1:nverts
      vertnums[j] = getNumberJ(vert_Nptr, verts[j], 0, 0)
      print(fstream, vertnums[j], " ")
    end
    print(fstream, "\n")
  end

  free(it)

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

"""
  Struct that contains all the MeshEntities that have nodes on them, as
  well as the number of nodes on each one

  It might be helpful to generalize this so it works on any MeshEntity,
  not just an element, but that would require more topology information

  **Public Fields**

   * entities: array of MeshEntities that have nodes on them, length `nentities`
   * nodecounts: array containing the number of nodes on each entity,
                 length `nentities`
   * nentities: number of entities
   * nnodes: total number of nodes on all the entities
"""
struct ElementNodeEntities
  # public fields 
  entities::Vector{Ptr{Void}}
  nodecounts::Vector{Int}  # number of nodes on each entity
  dimensions::Vector{Int}  # dimension of each entity
  nentities::Int  # = length of above arrays
  nnodes::Int # total number of nodes

  # private fields
  m_ptr::Ptr{Void}
  has_nodes_in_dim::Vector{Bool}
  dim::Int
  
  # temporary arrays
  down_entities::Array{Ptr{Void}, 1}
end

function ElementNodeEntities(m_ptr::Ptr{Void}, mshape_ptr::Ptr{Void},
                             dim::Integer)

  if dim == 2  # triangle
    num_type_per_entity = [3, 3, 1]
    dim = 2
  elseif dim == 3  # tetrahedron
    num_type_per_entity = [4, 6, 4, 1]
    numRegions = num_type_per_entity[4]
    dim = 3
  else
    throw(ErrorException("unsupported dimension in getNodeEntites, dim = $dim"))
  end

  num_nodes_per_entity = zeros(Int, 4)
  for i=1:3
    num_nodes_per_entity[i] = countNodesOn(mshape_ptr, i-1)
  end
  num_nodes_per_entity[4] = countNodesOn(mshape_ptr, 4)  # tetrahedron

  has_nodes_in_dim = Array{Bool}(4)
  for i=1:4
    has_nodes_in_dim[i] = num_nodes_per_entity[i] != 0
  end

  # compute number of entities
  nentities = 0
  nnodes = 0
  for i=0:dim
    if num_nodes_per_entity[i+1] > 0
      nentities += num_type_per_entity[i+1]
      nnodes += num_type_per_entity[i+1]*num_nodes_per_entity[i+1]
    end
  end

  entities = Array{Ptr{Void}}(nentities)
  nodecounts = Array{Int}(nentities)
  dimensions = Array{Int}(nentities)
  down_entities = Array{Ptr{Void}}(12)  # apf::Downward

  # get the nodecounts
  idx = 1
  for i=0:dim
    if num_nodes_per_entity[i+1] > 0
      for j=1:num_type_per_entity[i+1]
        nodecounts[idx] = num_nodes_per_entity[i+1]
        dimensions[idx] = i
        idx += 1
      end
    end
  end

  @assert idx == nentities + 1

  return ElementNodeEntities(entities, nodecounts, dimensions,  nentities,
                             nnodes, m_ptr,
                             has_nodes_in_dim, dim, down_entities)
end


"""
  Gets the apf::MeshEntities of the entities that have nodes on them.
  The entities are retrieve in order: vertices, edges, faces, and
  regions.  Entities of each dimension are ordered consistently with
  Pumi's `getDownward` function.

  **Inputs**

   * obj: an [`ElementNodeEntities`](@ref) object.  The `entities` field
          will be overwritten.
   * entity: an apf::MeshEntity* for the element.  The MeshEntities of
             the downward adjacent entities that have nodes on them will
             be retrieve
"""
function getNodeEntities(obj::ElementNodeEntities, entity::Ptr{Void})

  # entity must be of the same type passed to the constructor of obj

  # get entities in order of increasing dimension
  idx = 1
  for i=1:(obj.dim+1)
    if obj.has_nodes_in_dim[i]
      n = getDownward(obj.m_ptr, entity, i-1, obj.down_entities)
      for j=1:n
        obj.entities[idx] = obj.down_entities[j]
        idx += 1
      end
    end
  end

  return nothing
end



function getNodeEntities(m_ptr, mshape_ptr, entity)

  entity_type = getType(m_ptr, entity)
  eshape_ptr = getEntityShape(mshape_ptr, entity_type)
  nnodes = countNodes(eshape_ptr)
  downward_entities = Array{Ptr{Void}}(nnodes)  # holds mesh entities
  getNodeEntities(m_ptr, mshape_ptr, entity, downward_entities)

  return downward_entities
end

global const _getNodeEntities_retrieved_entities = Array{Ptr{Void}}(12)  # reusable storage
function getNodeEntities(m_ptr, mshape_ptr, entity, 
                         downward_entities::AbstractArray{Ptr{Void}})
# get the meshentities that have nodes on them 
# duplicates entities that have multiple nodes
# this only works for simplex elements
  
  entity_type = getType(m_ptr, entity)
  eshape_ptr = getEntityShape(mshape_ptr, entity_type)
  nnodes = countNodes(eshape_ptr)
#  downward_entities = Array{Ptr{Void}}(nnodes)  # holds mesh entities
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
    retrieved_entities[1] = Ptr{Void}(42)
    numEdges = getDownward(m_ptr, entity, 1, retrieved_entities)
        for i=1:numEdges
          insertN(downward_entities, retrieved_entities[i], vert_offset+num_edge_nodes*(i-1) + 1, num_edge_nodes)
        end
  end

  if has_face_nodes
     numFaces = getDownward(m_ptr, entity, 2, retrieved_entities)
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

function insertN(vec::AbstractArray{T}, element::T,  index::Integer, n::Integer) where T
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


#------------------------------------------------------------------------------
# Iteration over entities with nodes on them (according to a FieldShape)

import Base: start, next, done, iteratorsize, iteratoreltype, eltype, length, size
using Base: HasLength, HasEltype

"""
  This iterator iterates over MeshEntities that have nodes on them in a given
  FieldShape.  The iterator returns the MeshEntity* and its dimension.  The
  iterator does not need to be freed unless the loop is exited early using a
  `break` statement.  In this case, call `apf.free` on the iterator.

  **Constructor**

   * m_ptr: apf::Mesh*
   * fshape_ptr: apf::FieldShape*

  **Usage**

   ```
     for (e, e_dim) in apf.FieldEntityIt(m_ptr, fshape_ptr)
       # do stuff
     end
   ```
"""
mutable struct FieldEntityIt
  m_ptr::Ptr{Void}
  fshape_ptr::Ptr{Void}
  dim::Int  # current dimension
  lastdim::Int  # greatest dimension entity that has nodes on it
  mesh_dim::Int  # mesh dimension
  its::Vector{MeshIterator}
  has_nodes_in::Vector{Bool}
  len::Int  # number of entities
end

mutable struct FieldEntityItState
  e::Ptr{Void}
  dim::Int
end


using ODLCommonTools

function FieldEntityIt(m_ptr::Ptr{Void}, fshape_ptr::Ptr{Void})

  dim = 0
  lastdim = 0
  mesh_dim = getDimension(m_ptr)
  its = Vector{MeshIterator}(mesh_dim + 1)
  has_nodes_in = Vector{Bool}(mesh_dim + 1)
  for i=0:mesh_dim
    its[i+1] = MeshIterator(m_ptr, i)
    has_nodes_in[i+1] = hasNodesIn(fshape_ptr, i)
  end

  # find first dimension with nodes
  for i=1:length(has_nodes_in)
    if has_nodes_in[i]
      dim = i-1
      break
    end
  end

  for i=length(has_nodes_in):-1:1
    if has_nodes_in[i]
      lastdim = i - 1
      break
    end
  end

  len = 0
  for i=1:length(has_nodes_in)
    if has_nodes_in[i]
      len += apf.countJ(m_ptr, i-1)
    end
  end

  return FieldEntityIt(m_ptr, fshape_ptr, dim, lastdim, mesh_dim, its,
                       has_nodes_in, len)
end


function free(iter::FieldEntityIt)

  for it in iter.its
    free(it)
  end

  return nothing
end



function start(iter::FieldEntityIt)

  dim = iter.dim
  e = iterate(iter.its[dim+1])

  return FieldEntityItState(e, dim)
end


function next(iter::FieldEntityIt, state::FieldEntityItState)

  # the state always leads the current element by 1, so we can do the done
  # check
  new_val = (state.e, state.dim)
  state.e = iterate(iter.its[iter.dim+1])
  state.dim = iter.dim
  return new_val, state
end


function done(iter::FieldEntityIt, state)
  
  e = state.e
  isdone = e == C_NULL && iter.dim == (iter.lastdim)
  if e == C_NULL && !isdone
    # start iterating over next dimension entities
    iter.dim += 1
    state.e = iterate(iter.its[iter.dim+1])
    state.dim = iter.dim
  end

  if isdone
    free(iter)
  end

  return isdone
end

function iteratorsize(::Type{FieldEntityIt})
  return HasLength()
end

function length(iter::FieldEntityIt)
  return iter.len
end

function iteratoreltype(::Type{FieldEntityIt})
  return HasEltype()
end

function eltype(::Type{FieldEntityIt})
  return Tuple{Ptr{Void}, Int}
end

