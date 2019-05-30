__precompile__(false)
# functions to test the julia/PUMI interface
module apf

using MPI
using PumiConfig
using apf_types
using gmi_types

# import these names, so they can be accessed as apf.foo
import apf_types: IsoFuncJ, SolutionTransfers, MAInput, ModelEntity, MeshIterator, SubMeshData

# no names should exported because there should be higher level functions
# wrapping these
# but they are going to be exported anyways

#=
  declareNames() declares global constant variables that are used to construct
  the first argument to the ccall function (it must be a constant expression)

  names declared in C++ code as 'extern "C" ' will not be mangled, but some
  libraries do not do this, so these variables can be used to declare a human
  readable variable name for the mangled name
=#

#function declareNames()
# declare variables that hold the (possible mangled) names of c++ library functions
global const init_name = "initABC"
global const init2_name = "initABC2"
global const pushMeshRef_name = "pushMeshRef"
global const popMeshRef_name = "popMeshRef"

global const getConstantShapePtr_name = "getConstantShapePtr"
global const getMeshShapePtr_name = "getMeshShapePtr"

global const count_name = "count"
global const writeVtkFiles_name = "writeVtkFiles"

global const getModel_name = "getModel"
global const toModel_name = "toModel"
global const getModelType_name = "getModelType"
global const getModelTag_name = "getModelTag"

global const getMeshDimension_name = "getMeshDimension"
global const getType_name = "getType"
global const getDownward_name = "getDownward"
global const countAdjacent_name = "countAdjacent"
global const getAdjacent_name = "getAdjacent"
global const countBridgeAdjacent_name = "countBridgeAdjacent"
global const getBridgeAdjacent_name = "getBridgeAdjacent"
global const getAlignment_name = "getAlignment"

global const hasNodesIn_name = "hasNodesIn"
global const countNodesOn_name = "countNodesOn"
global const getEntityShape_name = "getEntityShape"
global const getOrder_name = "getOrder"

global const createMeshElement_name = "createMeshElement"
global const countIntPoints_name = "countIntPoints"
global const getIntPoint_name = "getIntPoint"
global const getIntWeight_name = "getIntWeight"
global const getJacobian_name = "getJacobian"

global const countNodes_name = "countNodes"
global const getValues_name = "getValues"
global const getLocalGradients_name = "getLocalGradients"
global const alignSharedNodes_name = "alignSharedNodes"

global const checkNums_name = "checkNums"
global const getVertCoords_name = "getVertCoords"
global const getVertCoords2_name = "getVertCoords2"
global const getEdgeCoords_name = "getEdgeCoords"
global const getEdgeCoords2_name = "getEdgeCoords2"
global const getFaceCoords_name = "getFaceCoords"
global const getFaceCoords2_name = "getFaceCoords2"
global const getElCoords_name = "getElCoords"
global const getElCoords2_name = "getElCoords2"
global const getAllEntityCoords_name = "getAllEntityCoords"
global const createNumberingJ_name = "createNumberingJ"
global const destroyNumbering_name = "destroyNumbering"
global const findNumbering_name = "findNumbering"
global const destroyNumberings_name = "destroyNumberings"
global const getNumberingShape_name = "getNumberingShape"
global const numberJ_name = "numberJ"
global const getNumberJ_name = "getNumberJ"
global const isNumbered_name = "isNumbered"
global const getDofNumbers_name = "getDofNumbers"
global const setNumberingOffset_name = "setNumberingOffset"
global const getElementNumbers_name = "getElementNumbers"
global const getMesh_name = "getNumberingMesh"
global const printNumberingName_name = "printNumberingName"

global const createDoubleTag_name = "createDoubleTag"
global const setDoubleTag_name = "setDoubleTag"
global const getDoubleTag_name = "getDoubleTag"


global const reorder_name = "reorder"
global const reorderXi_name = "reorderXi"

global const createIsoFunc_name = "createIsoFunc"
global const deleteIsoFunc_name = "deleteIsoFunc"
global const createSolutionTransfers_name = "createSolutionTransfers"
global const deleteSolutionTransfers_name = "deleteSolutionTransfers"
global const addSolutionTransfer_name = "addSolutionTransfer"
global const configureMAInput_name = "configureMAInput"
global const runMA_name = "runMA"
global const getAvgElementSize_name = "getAvgElementSize"

#global const createAnisoFunc_name = "createAnisoFunc"

#:field related functions
global const createPackedField_name = "createPackedField"
global const setComponents_name = "setComponents"
global const getComponents_name = "getComponents"
global const zeroField_name = "zeroField"
global const reduceField_name = "reduceField"
global const getCoordinateField_name = "getCoordinateField"
global const findField_name = "findField"
global const destroyField_name = "destroyField"
global const destroyFields_name = "destroyFields"


global const createSubMesh_name = "createSubMesh"
global const transferField_name = "transferField"
global const createSubMeshDG_name = "createSubMeshDG"
global const transferFieldDG_name = "transferFieldDG"

global const getFieldShape_name = "getFieldShape"

global const countPeers_name = "countPeers"
global const getPeers_name = "getPeers"
global const isShared_name = "isShared"
global const countRemotes_name = "countRemotes"
global const getRemotes_name = "getRemotes"

global const setPoint_name = "setPoint"
global const acceptChanges_name = "acceptChanges"
global const Verify_name = "Verify"
global const getPoint_name = "getPoint"

global const hasMatching_name = "hasMatching"
global const getSharing_name = "getSharing"
global const getNormalSharing_name = "getNormalSharing"
global const freeSharing_name = "freeSharing"
global const isOwned_name = "isOwned"
global const countCopies_name = "countCopies"
global const getCopies_name = "getCopies"
global const getOwner_name = "getOwner"
global const isSharedShr_name = "isSharedShr"

global const countMatches_name = "countMatches"
global const getMatches_name = "getMatches"
#end


# export low level interface functions
export init, loadMesh, initMesh, initGeometry, snapEdgeNodes, pushMeshRef,
       popMeshRef,
       getConstantShapePtr, getMeshShapePtr, countJ, writeVtkFiles,
       getMeshDimension, getType, getDownward, countAdjacent, getAdjacent,
       getAlignment, hasNodesIn, countNodesOn, getEntityShape, getOrder,
       getFieldShapeName,
       createMeshElement, countIntPoints, getIntPoint, getIntWeight,
       getJacobian, countNodes, getValues, getLocalGradients, alignSharedNodes,
       getVertCoords, getEdgeCoords, getFaceCoords, getElCoords,
       getAllEntityCoords, createNumberingJ, destroyNumbering, findNumbering,
       destroyNumberings,
       getNumberingShape, countNumberings, getNumbering,
       numberJ, getNumberJ, isNumbered, getDofNumbers,
       getElementNumbers, getNumberingMesh, getNumberingName,
       printNumberingName, createDoubleTag,
       setDoubleTag, getDoubleTag, reorder, reorderXi,
       createIsoFunc, createAnisoFunc,
       deleteIsoFunc, createSolutionTransfers, deleteSolutionTransfers,
       addSolutionTransfer, configureMAInput, runMA, getAvgElementSize, IsoFuncJ,
       SolutionTransfers, MAInput,
       createPackedField, setComponents,
       getComponents, zeroField, getFieldMesh,
       reduceField, getCoordinateField, findField,
       countFields, getField, destroyField, destroyFields,
       countBridgeAdjacent,
       getBridgeAdjacent, setNumberingOffset, createSubMesh, transferField

# iterator functors
export MeshIterator, iterate, iteraten, free, deref, getDimension, count

export createSubMeshDG, transferFieldDG, getFieldShape

export getModel, toModel, getModelType, getModelTag
export countPeers, getPeers
export countPeers, getPeers, countRemotes, getRemotes, isShared
export getEntity, incrementIt, resetIt

export setPoint, setParam,  acceptChanges, Verify, getPoint, getParam
export hasMatching, getSharing, getNormalSharing, freeSharing, isOwned,
       countCopies, getCopies, countMatches, getMatches, getOwner, isSharedShr

# SubMesh creation
export createSubMesh, getNewMesh, getOldMesh, writeNewMesh, getParentNumbering, 
       getNewMeshData, getGeoTag, SubMeshData




@doc """
  initilize the state of the interface library

  Naming conventions: mesh entity type 0 = vertex, 1 = edge, 2 = face, 3 = element

  Output variables:
    * downward_counts: returns the number of downward adjacencies of the first
                       element (useful if all elements are the same)
		       downward_counts[i,j] gives the number of downward
		       adjacenies entity type i has of type j
		       upward adjacencies ( i >= j) are not likely to be the
		       same for different elements, so they are set to zero
    * num_Entities : returns column vector containing number of entities of
                       each type in the mesh
    * m_ptr :        pointer to :Mesh2 that was loaded from the smb file
    * mshape_ptr :   pointer to the :FieldShape of the mesh
"""

function init(dmg_name::AbstractString, smb_name::AbstractString, order::Integer; load_mesh=true, shape_type=0 )
# initilize mesh interface
# initilize pointers to some value
# this is hack-ish -- there should be a better way to do this
#m_ptr = Ptr{Void}
#mshape_ptr = Ptr{Void}
num_Entities = zeros(Int32, 4)

m_ptr_array = Array{Ptr{Void}}(1)
mshape_ptr_array = Array{Ptr{Void}}(1)
n_arr = Array{Ptr{Void}}(4)

i = ccall( (init_name, pumi_libname), Int32, (Ptr{UInt8}, Ptr{UInt8},Ptr{Int32}, Ptr{Void}, Ptr{Void}, Ptr{Ptr{Void}}, Int32, Int32, Int32), dmg_name, smb_name, num_Entities, m_ptr_array, mshape_ptr_array, n_arr, order, load_mesh, shape_type )  # call init in interface library



#i = ccall( (init_name, pumi_libname), Int32, (Ptr{UInt8}, Ptr{UInt8},Ptr{Int32}, Ptr{Void}, Ptr{Void}), dmg_name, smb_name, num_Entities, m_ptr_array, mshape_ptr_array )  # call init in interface library

if ( i != 0)
  error("init failed, exiting ...")
end


return num_Entities, m_ptr_array[1], mshape_ptr_array[1], n_arr
end

#=
# 2d initilization
function init2(dmg_name::AbstractString, smb_name::AbstractString, order::Integer; load_mesh=true, shape_type=0)
# initilize mesh interface
# initilize pointers to some value
# order = order of shape functions to use, currently only necessary because
#      pumi does not properly support saving shape functions to files
# load_mesh : load mesh from file, or perform initilization functions on existing mesh (for re-initilizing after mesh adaptation)
# shape_type: 0 = lagrange, 1 = SBP

# this is hack-ish -- there should be a better way to do this
#m_ptr = Ptr{Void}
#mshape_ptr = Ptr{Void}
#downward_counts = zeros(Int32, 3,3);
num_Entities = zeros(Int32, 4, 1)

m_ptr_array = Array{Ptr{Void}}(1)
mshape_ptr_array = Array{Ptr{Void}}(1)
dim = Ref{Cint}()
n_arr = Array{Ptr{Void}}(4)

i = ccall( (init2_name, pumi_libname), Int32, (Ptr{UInt8}, Ptr{UInt8},Ptr{Int32}, Ptr{Void}, Ptr{Void}, Ptr{Cint}, Ptr{Ptr{Void}}, Int32, Int32, Int32), dmg_name, smb_name, num_Entities, m_ptr_array, mshape_ptr_array, dim, n_arr, order, load_mesh, shape_type )  # call init in interface library

if ( i != 0)
  throw(ErrorException("Init failed"))
end


return num_Entities, m_ptr_array[1], mshape_ptr_array[1], dim[], n_arr
end
=#
"""
  This function loads a mesh from a file and ensure the coordinate FieldShape
  is correct.  It does not perform any additional initialization, see
  [`initMesh`](@ref) for that.

  **Inputs**

   * dmg_name: name of geometry file, can be .null for null geometric model
   * smb_name: name of mesh file, not including file number (ie. abc.smb not
               abc0.smb for process 0)
   * order: order of coordinate field
   * shape_type: keyword argument, FieldShape identifier,
                 see getFieldShape in funcs1.cc, default 0

  **Outputs**

   * m_ptr: :Mesh2*
   * dim: dimension of the loaded mesh

  [`popMeshRef`](@ref) should be called on this pointer to free it.
"""
function loadMesh(dmg_name::AbstractString, smb_name::AbstractString,
                  order::Integer; shape_type::Integer=0)

  dim_ret = Array{Cint}(1)
  m_ptr = ccall( (:loadMesh, pumi_libname), Ptr{Void},
                 (Cstring, Cstring, Cint, Cint, Ptr{Cint}),
                 dmg_name, smb_name, shape_type, order, dim_ret)

  return m_ptr, dim_ret[1]
end

"""
  This function performs some initial numbering and counting of various
  mesh entities

  **Inputs**

   * m_ptr: a :Mesh* of an already loaded mesh
 
  **Outputs**

   * mshape_ptr: a :FieldShape* for the mesh coordinate field
   * num_entities: array of length 4 containing the number of entities of
                   each dimension, low to high, on this partition of the mesh.
                   All entities, including non-owned ones, are counted
   * n_array: array of length 4 containing a local numbering of the entities
              of each dimension, low to high.  If the mesh is 2 dimensional,
              accessing the 4th element of the array is undefined
"""
function initMesh(m_ptr::Ptr{Void})

  mshape_ptr_array = Array{Ptr{Void}}(1)
  num_entities = Array{Cint}(4)
  n_arr = Array{Ptr{Void}}(4)

  ccall( (:initMesh, pumi_libname), Void,
    (Ptr{Void}, Ptr{Cint}, Ptr{Ptr{Void}}, Ptr{Ptr{Void}}),
    m_ptr, num_entities, mshape_ptr_array, n_arr)

  return mshape_ptr_array[1], num_entities, n_arr
end


"""
  This function initializes the field used for storing the CAD parametric
  coordinates of the mid-edge nodes.

  **Inputs**

   * m_ptr: the apf::Mesh2 *

  **Outputs**

   * f_ptr: pointer to the apf::Field holding the CAD parametric coordinates
            of the mid-edge nodes.  If geometric model does not have parametric
            coordinates or if the mesh is linear, the null pointer is returned.
            The field has one node per edge with 2 values at each node.
"""
function initGeometry(m_ptr::Ptr{Void})

  f_ptr = ccall( (:initGeometry, pumi_libname), Ptr{Void}, (Ptr{Void},), m_ptr)

  return f_ptr
end


"""
  This function snaps the mid-edge nodes to the geometry and updates the
  apf::Field with the new parametric values.  Only call this function if the
  geometry model is capable of snapping.

  **Inputs**

   * m_ptr: the mesh pointer
   * f_ptr: the Field pointer
"""
function snapEdgeNodes(m_ptr::Ptr{Void}, f_ptr::Ptr{Void})

  ccall( (:snapEdgeNodes, pumi_libname), Void, (Ptr{Void}, Ptr{Void}), m_ptr, f_ptr)

  return nothing
end


function pushMeshRef(m_ptr::Ptr{Void})

  ccall( (pushMeshRef_name, pumi_libname), Void, (Ptr{Void},), m_ptr)
end

function popMeshRef(m_ptr::Ptr{Void})
  ccall( (popMeshRef_name, pumi_libname), Void, (Ptr{Void},), m_ptr)
end

function countMeshRefs(m_ptr::Ptr{Void})
  n = ccall( (:countMeshRefs, pumi_libname), Cint, (Ptr{Void},), m_ptr)

  return n
end



# no longer needed
function getMeshShapePtr(m_ptr::Ptr{Void})

  mshape_ptr = ccall( (getMeshShapePtr_name, pumi_libname), Ptr{Void}, (Ptr{Void},), m_ptr )
  return mshape_ptr
end


function getConstantShapePtr(dimension::Integer)
  mshape_ptr = ccall( ( getConstantShapePtr_name, pumi_libname), Ptr{Void}, (Int32,), dimension)

  return mshape_ptr

end

#------------------------------------------------------------------------------
# Iterators

import apf_types.MeshIterator


function MeshIterator(m_ptr::Ptr{Void}, dim::Integer)

  it = ccall( (:begin, pumi_libname), _MeshIterator, (Ptr{Void}, Cint,), m_ptr, dim)

  len = countJ(m_ptr, dim)
  return MeshIterator(it, m_ptr, len)
end

function free(it::MeshIterator)

  ccall((:end, pumi_libname), Void, (Ptr{Void}, _MeshIterator), it.m_ptr, it.p)

end

function iterate(it::MeshIterator)

  me = ccall((:iterate, pumi_libname), Ptr{Void}, (Ptr{Void}, _MeshIterator), it.m_ptr, it.p)

  return me
end

function iteraten(it::MeshIterator, n::Integer)

  me = ccall((:iteraten, pumi_libname), Ptr{Void}, (Ptr{Void}, _MeshIterator, Cint), it.m_ptr, it.p, n)

  return me
end


function deref(it::MeshIterator)

  me = ccall((:deref, pumi_libname), Ptr{Void}, (Ptr{Void}, _MeshIterator), it.m_ptr, it.p)

  return me
end

import Base: start, next, done, iteratorsize, iteratoreltype, eltype, length, size
using Base: HasLength, HasEltype


function start(it::MeshIterator)

  return iterate(it)
end


function next(iter::MeshIterator, state)

  # the state always leads the current element by 1, so we can do the done
  # check
  e_next = iterate(iter)
  new_state = e_next
  return state, new_state
end


function done(iter::MeshIterator, state)

  isdone = state == C_NULL
  if isdone
    free(iter)
  end

  return isdone
end


function iteratorsize(::Type{MeshIterator})
  return HasLength()
end

function length(iter::MeshIterator)
  return iter.len
end

function iteratoreltype(::Type{MeshIterator})
  return HasEltype()
end

function eltype(::Type{MeshIterator})
  return Ptr{Void}
end


#------------------------------------------------------------------------------

function getDimension(m_ptr::Ptr{Void})

  val = ccall( (:getDimension, pumi_libname), Cint, (Ptr{Void},), m_ptr)

  return val
end



function countJ(m_ptr, dimension::Integer)
# returns the number of entities of a particular dimension

  i = ccall( (count_name, pumi_libname), Cint, (Ptr{Void}, Cint), m_ptr, dimension)
  return i
end


"""
  Writes VTK files.  The keyword argument `writeall` determines which Numberings
  and Fields are written to the file.  By default only those with a FieldShape
  compatible with the the coordinate field are written.

  **Inputs**

   * name: file name
   * m_ptr: mesh pointer
   * f_ptrs: Vector of apf::Field*, defaults to array of zero length
   * n_ptrs: Vector of apf::Numbering*, defaults to array of zero length
   * gn_ptrs: Vector of apf::GlobalNumbering*, defaults to array of zero length

  **Keyword Arguments**

   * writeall: if true, write fields that are not representable cleanly in
               vtk file, default false

  For `f_ptrs`, `n_ptrs`, and `gn_ptrs`, only the objects in the vector will
  be considered for writing to the vtk file (according to `writeall`).  If the
  vector has zero length, all objects will be considered.
"""
function writeVtkFiles(name::AbstractString, m_ptr,
                       f_ptrs=Vector{Ptr{Void}}(0),
                       n_ptrs=Vector{Ptr{Void}}(0),
                       gn_ptrs=Vector{Ptr{Void}}(0);
                       writeall::Bool=false)
# write vtk files to be read by paraview

ccall( (writeVtkFiles_name, pumi_libname), Void,
      (Ptr{UInt8}, Ptr{Void}, CppBool, Ptr{Ptr{Void}}, Cint, Ptr{Ptr{Void}},
       Cint, Ptr{Ptr{Void}}, Cint), name, m_ptr, writeall,
      f_ptrs, length(f_ptrs), n_ptrs, length(n_ptrs), gn_ptrs, length(gn_ptrs))
  return nothing
end

function getModel(m_ptr::Ptr{Void})

  g = ccall( (getModel_name, pumi_libname), Ptr{Void}, (Ptr{Void},), m_ptr)
  return Model(g)
end


function toModel(m_ptr, entity_ptr)
# get the model entity a mesh entity is classified on

  model_entity_ptr = ccall( (toModel_name, pumi_libname), Ptr{Void}, (Ptr{Void}, Ptr{Void}), m_ptr, entity_ptr)

  return ModelEntity(model_entity_ptr)
end

function getModelType(m_ptr, model_entity::ModelEntity)
# get the model entity *dimension*

  dim = ccall( (getModelType_name, pumi_libname), Int32, (Ptr{Void}, ModelEntity), m_ptr, model_entity)

  return dim
end

function getModelTag(m_ptr, model_entity::ModelEntity)
# get the dimension unique identifier of the model entity

  tag = ccall( (getModelTag_name, pumi_libname), Int32, (Ptr{Void}, ModelEntity), m_ptr, model_entity)

  return tag
end





function getMeshDimension(m_ptr)
# get mesh dimension

  i = ccall( (getMeshDimension_name, pumi_libname), Int32, (Ptr{Void},), m_ptr)
  return i
end


function getType(m_ptr, entity)
# get mesh dimension

  i = ccall( (getType_name, pumi_libname), Int32, (Ptr{Void}, Ptr{Void}), m_ptr, entity)
  return i
end

function getDownward(m_ptr, entity, dimension::Integer)
  
  # create array that can fit max number of downward adjacencies
  downwards = Array{Ptr{Void}}(12)


  i = ccall( (getDownward_name, pumi_libname), Int32, (Ptr{Void}, Ptr{Void}, Int32, Ptr{Void}), m_ptr, entity, dimension, downwards)

  return downwards[1:i], i

end

function getDownward(m_ptr, entity, dimension::Integer, arr::AbstractArray{Ptr{Void}})
# populate arr with the downward entiteis of the specified dimension
# arr is not checked for size

  i = ccall( (getDownward_name, pumi_libname), Int32, (Ptr{Void}, Ptr{Void}, Int32, Ptr{Ptr{Void}}), m_ptr, entity, dimension, arr)

  return  i

end



function countAdjacent(m_ptr, entity, dimension::Integer)
# counts the number of upward adjacencies of meshentity with given dimension
# use with getAdjacent to fetch them

i = ccall( (countAdjacent_name, pumi_libname), Int32, (Ptr{Void}, Ptr{Void}, Int32), m_ptr, entity, dimension)

return i

end

function getAdjacent(num_adjacent::Integer)
# returns an array of MeshEntity* that are the upward adjacencies fetched by countAdjacent
# this could be made  more efficient by taking in an array and resizing it if it is too small
  adjacencies_ret = Array{Ptr{Void}}(num_adjacent)  # create the array

#  ccall( (getAdjacent_name, pumi_libname), Void, (Ptr{Void},), adjacencies_ret)
   getAdjacent(adjacencies_ret)

  return adjacencies_ret
end

function getAdjacent(arr::AbstractArray{Ptr{Void}})


  ccall( (getAdjacent_name, pumi_libname), Void, (Ptr{Ptr{Void}},), arr)

  return nothing
end


function countBridgeAdjacent(m_ptr, entity, bridge_dimension::Integer, target_dimension::Integer)
# counts the number of entities to be retrieved
# this is a second order adjacency

  i = ccall( (countBridgeAdjacent_name, pumi_libname), Int32, (Ptr{Void}, Ptr{Void}, Int32, Int32), m_ptr, entity, bridge_dimension, target_dimension)

  return i
end

function getBridgeAdjacent(arr::AbstractArray{Ptr{Void}})

  ccall( (getBridgeAdjacent_name, pumi_libname), Void, (Ptr{Void},), arr)

  return nothing
end

function getAlignment(m_ptr, elem, elem_boundary)
# gets which, flip, and rotate for elem_boundary, the mesh entity shared
# between elem and another element

# create variables
#=
which = convert(Int32, 42)
flip = convert(UInt8, 42)
rotate = convert(Int32, 42)
=#
which = Array{Int32}(1)
flip = Array{CppBool}(1)
rotate = Array{Int32}(1)

  ccall( (getAlignment_name, pumi_libname), Void, (Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Int32}, Ptr{CppBool}, Ptr{Int32}), m_ptr, elem, elem_boundary, which, flip, rotate)

  return which[1], convert(Bool, flip[1]), rotate[1]
end


function hasNodesIn(mshape_ptr, dimension::Integer)
# check whether this FieldShape* has nodes on entities of a given dimension

  i = ccall( (hasNodesIn_name, pumi_libname), CppBool, (Ptr{Void}, Int32), mshape_ptr, dimension)

  i_bool = convert(Bool, i)
#  println("hasNodesIn returned", i_bool)

  return i_bool

end

function countNodesOn(mshape_ptr, entity_type::Integer)
# count the number of nodes on an entity of the specified type (:Mesh::Type)
  i = ccall( (countNodesOn_name, pumi_libname), Int32, (Ptr{Void}, Int32), mshape_ptr, entity_type)

  return i

end

function getEntityShape(mshape_ptr, entity_type::Integer)
# get the EntityShape* (object describing shape functions) ofa given entity type

  eshape_ptr = ccall( (getEntityShape_name, pumi_libname), Ptr{Void}, (Ptr{Void}, Int32), mshape_ptr, entity_type)

  return eshape_ptr
end

function getOrder(mshape_ptr::Ptr{Void})

  order = ccall( (getOrder_name, pumi_libname), Cint, (Ptr{Void},), mshape_ptr)

  return Int(order)  # make order an Int for convenience
end


function getFieldShapeName(mshape_ptr::Ptr{Void})

  str= ccall( (:getFieldShapeName, pumi_libname), Cstring, (Ptr{Void},), mshape_ptr)

  return unsafe_string(str)
end



function createMeshElement(m_ptr, entity)
# creates a MeshElement from a MeshEntity

  mel_ptr = ccall( ( createMeshElement_name, pumi_libname), Ptr{Void}, (Ptr{Void}, Ptr{Void}), m_ptr, entity)
  return mel_ptr
end

function countIntPoints(mel_ptr, order::Integer)
# count the number of integration points needed to achieve specified order of
# accuracy.  MeshElement can be edges as well as 2D and 3D regions.
# not sure what happens if it is a vertex.

  i = ccall( (countIntPoints_name, pumi_libname), Int32, (Ptr{Void}, Int32), mel_ptr, order)
  return i
end

function getIntPoint(mel_ptr, order::Integer, point::Integer)
# returns a vector containing the parent coordinates of the integration point
# order = order of accuracy of integration
# point = which integration point (1 thru number omef points )

   coords = zeros(3)
   ccall( (getIntPoint_name, pumi_libname), Void, (Ptr{Void}, Int32, Int32, Ptr{Float64}), mel_ptr, order, point-1, coords)
   return coords
end

function getIntWeight(mel_ptr, order::Integer, point::Integer)
# gets the weight corresponding to the integration point (see getIntPoint)

  i = ccall( (getIntWeight_name, pumi_libname), Float64, (Ptr{Void}, Int32, Int32), mel_ptr, order, point-1)
  return i
end


function getJacobian(mel_ptr, coords::Array{Float64,1})
# gets the jacobian of a mesh element at a particular location in parametric space
# coords should be a vector of length 3

  jac = zeros(3,3)

  ccall( (getJacobian_name, pumi_libname), Void, (Ptr{Void}, Ptr{Float64}, Ptr{Float64}), mel_ptr, coords, jac)

  return jac.'  # return the transpose because of row major ordering
end

function countNodes(eshape_ptr)
# get the total number of nodes related to an entity (including downward adjacencies)

  i = ccall( (countNodes_name, pumi_libname), Int32, (Ptr{Void},), eshape_ptr)
  return i
end

function getValues(m_ptr::Ptr{Void},  eshape_ptr, coords::Array{Float64,1}, numN::Integer)
#  gets an array of shape function values at the specified coordinates
# numN is the number of points affecting the entity that was used to get eshape_ptr
# coords must be a vector of length 3
# the output is a vector of length numN

  vals = zeros(numN)
  ccall( (getValues_name, pumi_libname), Void, (Ptr{Void}, Ptr{Void}, Ptr{Float64}, Ptr{Float64}), m_ptr, eshape_ptr, coords, vals)

  return vals
end


function getLocalGradients(m_ptr::Ptr{Void}, eshape_ptr, coords::Array{Float64,1}, numN::Integer)
#  gets an array of shape function derivatives at the specified coordinates
# numN is the number of points affecting the entity that was used to get eshape_ptr
# coords must be a vector of length 3
# the output is a matrix of dimension 3 x numN

  vals = zeros(3, numN)
  ccall( (getLocalGradients_name, pumi_libname), Void, (Ptr{Void}, Ptr{Void}, Ptr{Float64}, Ptr{Float64}), m_ptr, eshape_ptr, coords, vals)

  return vals
end


function alignSharedNodes(eshape_ptr, m_ptr, elem, shared, order::Array{Int32, 1})
# get the array that transofrms teh local element order to the canonical order
# the length of order must be the number of nodes classified on 
# the shared entity
# elem is a pointer to the MeshEntity that is the element
# shared is a poniter to the MeshEntity that is shared between elements

  # bounds check
 mshape_ptr = getMeshShapePtr(m_ptr)
 shared_type = getType(m_ptr, shared)
 numnodes = countNodesOn(mshape_ptr, shared_type)

 if length(order) < numnodes
   println("alignSharedNodes order array too short")
   return nothing
 end

 ccall( (alignSharedNodes_name, pumi_libname), Void, (Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Int32}), esahpe_ptr, m_ptr, elem, shared, order)

 return nothing
end


function getVertCoords(m_ptr::Ptr{Void}, entity, coords::Array{Float64, 2}, m::Integer, n::Integer)
# coords is array to put coordsinates in, must be 3 by 1,
# m, n are number of rows, columns in coords, respectively

#coords = Array{Float64}(3, 2)   # pass an array 3 by n (3 coordinates each for n points)
#(m,n) = size(coords)
# pass reversed m,n because C arrays are row-major
ccall( (getVertCoords2_name, pumi_libname), Void, (Ptr{Void}, Ptr{Void}, Ptr{Float64}, Int32, Int32), m_ptr, entity, coords, n, m) 

#println("\n from julia, coords = ", coords)

return coords
end



function getEdgeCoords(m_ptr::Ptr{Void}, entity, coords::Array{Float64, 2}, m::Integer, n::Integer)
# coords is array to put coordinates in, must be 3 by 2
# m,n = number of rows, columns in coords, respectively

i = ccall( (getEdgeCoords2_name, pumi_libname), Int, (Ptr{Void}, Ptr{Void}, Ptr{Float64}, Int32, Int32), m_ptr, entity, coords, n, m);

if ( i != 0)
  throw(ErrorException("Error in getEdgeCoords"))
end

#println("in julia, coords = ", coords)

end



function getFaceCoords(m_ptr::Ptr{Void}, entity::Ptr{Void}, coords::AbstractMatrix, m::Integer, n::Integer)
# get coordinates of points on a face, in order
# coords is array to populate with coordinates, must by 3 by number of points
# on a face
# m,n = number of rows, columns in coords, respectively

# reverse m and n because C is row major
i = ccall( (getFaceCoords2_name, pumi_libname), Int32, (Ptr{Void}, Ptr{Void}, Ptr{Float64}, Int32, Int32), m_ptr, entity, coords, n, m);

if ( i != 0)
  throw(ErrorException("Error in getFaceCoords"))
end

#println("in julia, coords = ", coords)

end


function getElCoords(m_ptr::Ptr{Void}, entity, coords::Array{Float64, 2}, m::Integer, n::Integer)
# get coordinates of points in an element, in order
# coords is array to populate with coordinates, must by 3 by number of points
# on a element
# m,n = number of rows, columns in coords, respectively

# reverse m and n because C is row major
i = ccall( (getElCoords2_name, pumi_libname), Int, (Ptr{Void}, Ptr{Void}, Ptr{Float64}, Int32, Int32), m_ptr, entity, coords, n, m);

if ( i != 0)
  throw(ErrorException("Error in getElCoords"))
end

#println("in julia, coords = ", coords)

end

# get the coordinates of all nodes downward adjacent to an element
# the array is populated contiguously in memory, and should be
# dim x numNodesPerElement in the FieldShape of the coordinates
function getAllEntityCoords(m_ptr::Ptr{Void}, entity::Ptr{Void}, 
                            coords::AbstractArray{Float64})

  dims = length(size(coords))
  @assert stride(coords, 1) == 1
  for i=2:dims
    @assert stride(coords, i) == size(coords, i-1)
  end

  ccall( (getAllEntityCoords_name, pumi_libname), Void, (Ptr{Void}, Ptr{Void}, Ptr{Float64}), m_ptr, entity, coords);

  return nothing
end

function createNumberingJ(m_ptr::Ptr{Void}, name::AbstractString, field::Ptr{Void}, components::Integer)
# create a generally defined numbering, get a pointer to it
# this just passes through to :createNumbering
# field is an :FieldShape*

numbering_ptr = ccall( (createNumberingJ_name, pumi_libname), Ptr{Void}, (Ptr{Void}, Ptr{UInt8}, Ptr{Void}, Int32), m_ptr, name, field, components)

return numbering_ptr
end

function destroyNumbering(n_ptr::Ptr{Void})

  ccall( (destroyNumbering_name, pumi_libname), Void, (Ptr{Void},), n_ptr)

  return nothing
end

function findNumbering(m_ptr::Ptr{Void}, name::String)

  n_ptr = ccall( (findNumbering_name, pumi_libname), Ptr{Void}, (Ptr{Void}, Cstring), m_ptr, name)

  return n_ptr
end

"""
  Destroys all :Numberings associated with a given mesh, except
  those specified

  **Inputs**

   * m_ptr: :Mesh*
   * n_save: array of :Numbering* objects to not destroy (optional)
"""
function destroyNumberings(m_ptr::Ptr{Void}, n_save=Array{Ptr{Void}}(0))

  ccall( (destroyNumberings_name, pumi_libname), Void, (Ptr{Void}, Ptr{Ptr{Void}}, Cint), m_ptr, n_save, length(n_save))

  return nothing
end

function getNumberingShape(n_ptr::Ptr{Void})

  fshape = ccall( (getNumberingShape_name, pumi_libname), Ptr{Void}, (Ptr{Void},), n_ptr)
end

function countNumberings(m_ptr::Ptr{Void})

  n = ccall( (:countNumberings, pumi_libname), Cint, (Ptr{Void},), m_ptr)

  return n
end

function getNumbering(m_ptr::Ptr{Void}, i::Integer)

  n_ptr = ccall( (:getNumbering, pumi_libname), Ptr{Void}, (Ptr{Void}, Cint), m_ptr, i)

  return n_ptr
end

function numberJ(numbering_ptr, entity, node::Integer, component::Integer, number::Integer)
#  assign a number to a node defined on a mesh entity

  i = ccall( (numberJ_name, pumi_libname), Int32, (Ptr{Void}, Ptr{Void}, Int32, Int32, Int32), numbering_ptr, entity, node, component, number)

# if i == number
#   println("numbering successful")
# else
#   println("numbering failure, number sent = $number, number returnee = $i")
# end

return nothing

end


function getNumberJ(n_ptr, entity, node::Integer, component::Integer)
#  retrieve the number  of a node defined on a mesh entity for a particular numbering

i = ccall( (getNumberJ_name, pumi_libname), Int32, (Ptr{Void}, Ptr{Void}, Int32, Int32), n_ptr, entity, node, component)
return i

end

function isNumbered(n_ptr, entity, node::Integer, component::Integer)
# check whether a node is numbered or not, returns a Julia Bool

i = ccall( (isNumbered_name, pumi_libname), CppBool, (Ptr{Void}, Ptr{Void}, Int32, Int32), n_ptr, entity, node, component)
return i != 0

end



function getDofNumbers(n_ptr, entities::AbstractArray{Ptr{Void}}, node_offsets::AbstractArray{UInt8}, nodemap::AbstractArray{UInt8}, element::Ptr{Void}, dofnums::AbstractArray{Int32})

  ccall( (getDofNumbers_name, pumi_libname), Int32, (Ptr{Void}, Ptr{Void}, Ptr{UInt8}, Ptr{UInt8}, Ptr{Void}, Ptr{Int32}), n_ptr, entities, node_offsets, nodemap, element, dofnums)

  return nothing
end


function setNumberingOffset(n_ptr, offset::Integer)

  ccall( (setNumberingOffset_name, pumi_libname), Void, (Ptr{Void}, Int32), n_ptr, offset)

  return nothing
end



function getElementNumbers(n_ptr, entity, num_dof::Integer, nums::Array{Int32, 1})
# get the numbers of all dofs affecting an element in the canonical order
# n_ptr is the numbering to use
# entity is the entity to use
# num_dof is the number of degrees of freedom (nnodes*ndofpernode)
# nums is an array of Int32s long enough to hold all the numbers
# the dofs numbers are returned in a column vector, all dofs of a node
# are sequential

  if length(nums) < num_dof
    println("Warning: getElementNumbers nums array too short for number of nodes")
    return nothing
  end

  ccall( (getElementNumbers_name, pumi_libname), Void, (Ptr{Void}, Ptr{Void}, Int32, Ptr{Int32}), n_ptr, entity, num_dof, nums)

  return nothing

end

function getNumberingMesh(n_ptr)
# get the mesh a numbering is defined on

  m_ptr = ccall( (getMesh_name, pumi_libname), Ptr{Void}, (Ptr{Void},), n_ptr)
  return m_ptr
end


function getNumberingName(n_ptr::Ptr{Void})

  ptr = ccall( (:getNumberingName, pumi_libname), Cstring, (Ptr{Void},), n_ptr)
  return unsafe_string(ptr)
end

function printNumberingName(numbering)
# print the name of a numbering

  ccall( (printNumberingName_name, pumi_libname), Void, (Ptr{Void},), numbering)
  return nothing
end

function createDoubleTag( m_ptr, name::AbstractString, components::Integer)

  tag_ptr = ccall( (createDoubleTag_name, pumi_libname), Ptr{Void}, (Ptr{Void}, Ptr{UInt8}, Int32), m_ptr, name, components)

  return tag_ptr
end

function setDoubleTag( m_ptr, entity, tag_ptr, data)
  # data is the array of data to be stored in the tag
  ccall( (setDoubleTag_name, pumi_libname), Void, (Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Float64}), m_ptr, entity, tag_ptr, data)

  return nothing
end

function getDoubleTag( m_ptr, entity, tag_ptr, data)

ccall( (getDoubleTag_name, pumi_libname), Void, (Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Float64}), m_ptr, entity, tag_ptr, data)


  return nothing
end



# these functions are not to be used meshes that are either large or 3d
@doc """
### PumiInterface.reorder

  This function performs node reordering according to the algorithm in 
  "Adjacency-based Data Reordering Algoirthm for Acceleration of Finite 
   Element Computations" by Zhou, Sahni, Shephard, Carothers, Jansen

   Arguments:
   m_ptr: mesh pointer
   ndof: total number of degrees of freedom to be labelled (exlucding any
         fixed dofs not to be given a number)
   comp: number of dofs per node
   node_statusNumbering: :Numbering* that tells what the status of a dof
                         is.  If the value for a dof >= 2, then it gets 
			 labelled.  Must have ncomp components per node.
                         If C_NULL, all dofs will be numbered.
  nodeNums:  an already created :Numbering* that the dof numbers are
             written to. Its :Fieldshape must be the same as the mesh.
	     The numbering is 1-based.  Any dof with status < 2 will be given
	     the number zero.

  elNums: an already created :Numbering* that the element numbers are 
          written to.  If C_NULL, this numbering will not be used.
  start_coords: coordinates of a point. The mesh vertex classified on a model
         vertex closest to this point is used as the starting entity for the 
	 reordering.  Note that the reordering algorithm assigns dof numbers
	 in reverse order (ie. high to low rather than low to high).  Array
         of length m->getDimension()

"""->
function reorder(m_ptr, ndof::Integer, ncomp::Integer, node_statusN_ptr,
                 nodeNums, elNums, start_coords::Vector{Cdouble})

  @assert length(start_coords) == 3
  ccall( (reorder_name, pumi_libname), Void, (Ptr{Void}, Int32, Int32, Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Cdouble}),  m_ptr, ndof,  ncomp, node_statusN_ptr, nodeNums, elNums, start_coords)

  return nothing

end


"""
  This function populates the `xiNums` apf::Numbering with the numbering of
  the geometric degrees of freedom. The number of degrees of freedom a 
  MeshEntity has is exactly equal to the dimension of the geometry entity it
  is classified on.  Therefore, a mesh vertex classified on a geometric vertex
  has zero degrees of freedom, while a vertex classified on a geometric edge
  has one.  This effectively constrains the MeshEntity to remain on the
  geometric entity it is classified on.  Any unused entries in the `xiNums`
  have value `numXiDof + 1`.

  **Inputs**

   * m_ptr: apf::Mesh*
   * xiNums: apf::Numbering* to be written to
   * start_coords: the MeshEntity closest to this these coordinates will be
                   used to start the numbering process.  Note that reverse
                   numbering is performed, so this mesh entity will have the
                   highest dof number.  This array must be of length 3, even in
                   2D.
  **Outputs**

   * numXiDof: the number of geometric degrees of freedom
"""
function reorderXi(m_ptr::Ptr{Void}, xiNums::Ptr{Void},
                   start_coords::Vector{Cdouble})

  @assert length(start_coords) == 3

  val = ccall( (reorderXi_name, pumi_libname), Cint, (Ptr{Void}, Ptr{Void},
          Ptr{Cdouble}), m_ptr, xiNums, start_coords)

  return val
end

#------------------------------------------------------------------------------
# mesh adapation functions

"""
  Create an [`IsoFuncJ`](@ref)

  **Inputs**

   * m_ptr: an :Mesh2*
   * f: an :Field* specifying the desired mesh size at each coordinate
        node of mesh (note: coordinate node, not solution node)
"""
function createIsoFunc(m_ptr::Ptr{Void}, f::Ptr{Void})

  ptr = ccall( (createIsoFunc_name, pumi_libname), Ptr{Void}, (Ptr{Void}, Ptr{Void}), m_ptr, f)

  return IsoFuncJ(ptr)
end

"""
  Free an [`IsoFuncj`](@ref)

  **Inputs**

   * func: an IsoFuncJ

  **Outputs**

    none
"""
function deleteIsoFunc(func::IsoFuncJ)

  ccall( (deleteIsoFunc_name, pumi_libname), Void, (IsoFuncJ,), func)

  return nothing
end

"""
  Create an object that manages the transfer of fields from the original
  mesh to the adapted mesh.  Use [`addSolutionTransfer`](@ref) to add
  fields to be transfered one by one.

  **Inputs**

    * none

  **Outputs**

   * a [`SolutionTransfers`](@ref)
"""
function createSolutionTransfers()

  ptr = ccall( (createSolutionTransfers_name, pumi_libname), Ptr{Void}, () )

  return SolutionTransfers(ptr)
end

"""
  Frees a [`SolutionTransfers`](@ref) object.

  **Inputs**

   * soltrans: the object to be deleted

  **Outputs**

   * none

  **Implementation Notes**

   This actually deletes the `SolutionTransfers` object as well as all the
   `SolutionTransfer` object inside of it.  Currently the individual
   `SolutionTransfer` objects are not exposed to Julia, so this is a
   non-issue.
"""
function deleteSolutionTransfers(soltrans::SolutionTransfers)

  ccall( (deleteSolutionTransfers_name, pumi_libname), Void, (SolutionTransfers,), soltrans)

  return nothing
end

"""
  Add an :Field* to be transferred to the new mesh

  **Inputs**

   * soltrans: the [`SolutionTransfers`](@ref) object
   * f: the :Field*

  **Outputs**

   none
"""
function addSolutionTransfer(soltrans::SolutionTransfers, f::Ptr{Void})

  ccall( (addSolutionTransfer_name, pumi_libname), Void, (SolutionTransfers, Ptr{Void}), soltrans, f)

  return nothing
end

"""
  Create the input configuration object for MeshAdapt

  **Inputs**

   * m_ptr: an :Mesh2*
   * isofunc: an [`IsoFuncJ`](@ref)
   * soltrans: a [`SolutionTransfers`](@ref)

  **Outputs**

   * an [`MAInput`](@ref)
"""
function configureMAInput(m_ptr::Ptr{Void}, isofunc::IsoFuncJ, soltrans::SolutionTransfers)

  ptr = ccall( (configureMAInput_name, pumi_libname), Ptr{Void},
               (Ptr{Void}, IsoFuncJ, SolutionTransfers), m_ptr, isofunc, soltrans)

  return MAInput(ptr)
end

"""
  Run MeshAdapt

  **Inputs**

   * input: an [`MAInput](@ref)

  **Outputs**

   none
"""
function runMA(input::MAInput)

  ccall( (runMA_name, pumi_libname), Void, (MAInput,), input)

  return nothing
end

"""
  Gets the average size of each element (ie. the average edge length)

  **Inputs**

   * m_ptr: :Mesh*
   * el_N: an :Numbering* for the elements

  **Inputs/Outputs**

   * el_sizes: vector to be overwritten with the average size of each element
"""
function getAvgElementSize(m_ptr::Ptr{Void}, el_N::Ptr{Void}, el_sizes::AbstractVector{Cdouble})

  ccall( (getAvgElementSize_name, pumi_libname), Void, (Ptr{Void}, Ptr{Void}, Ptr{Cdouble}), m_ptr, el_N, el_sizes)

  return nothing
end

#------------------------------------------------------------------------------
# :Field related function
# needed for automagical solution transfer

function createPackedField(m_ptr, fieldname::AbstractString, numcomponents::Integer, fshape=C_NULL)
# create a field with the specified number of componenets per node
# returns a pointer to the field

  f_ptr = ccall((createPackedField_name, pumi_libname), Ptr{Void}, (Ptr{Void}, Ptr{UInt8}, Int32, Ptr{Void}), m_ptr, fieldname, numcomponents, fshape)

  return f_ptr
end

function setComponents(f_ptr, entity_ptr, node::Integer, components::AbstractVector)
  ccall((setComponents_name, pumi_libname), Void, (Ptr{Void}, Ptr{Void}, Int32, Ptr{Float64}), f_ptr, entity_ptr, node, components)

  return nothing
end

function getComponents(f_ptr, entity_ptr, node::Integer, components::AbstractVector)

  ccall( (getComponents_name, pumi_libname), Void, (Ptr{Void}, Ptr{Void}, Int32, Ptr{Float64}), f_ptr, entity_ptr, node, components)

  return nothing
end

function zeroField(f_ptr)
  ccall( (zeroField_name, pumi_libname), Void, (Ptr{Void},), f_ptr)
end

function getFieldMesh(f_ptr::Ptr{Void})
# get the mesh a numbering is defined on

  m_ptr = ccall( (:getFieldMesh, pumi_libname), Ptr{Void}, (Ptr{Void},), f_ptr)
  return m_ptr
end



"""
  Apply a reduction operation along the partition boundaries of a Field
  to make the field globally consistent.  The field must have values
  written to it by all processes before this function is applied.

  **Inputs**

   * f_ptr: :Field*
   * shr_ptr:: :Sharing*
   * reduce_op: 0 = sum, 1 = min, 2 = max
"""
function reduceField(f_ptr::Ptr{Void}, shr_ptr::Ptr{Void}, reduce_op::Integer)

  ccall( (reduceField_name, pumi_libname), Void, (Ptr{Void}, Ptr{Void}, Cint), f_ptr, shr_ptr, reduce_op)

end

function getCoordinateField(m_ptr::Ptr{Void})

  ccall( (getCoordinateField_name, pumi_libname), Ptr{Void}, (Ptr{Void},), m_ptr)
end

function findField(m_ptr::Ptr{Void}, fieldname::String)

  f_ptr = ccall( (findField_name, pumi_libname), Ptr{Void}, (Ptr{Void}, Cstring), m_ptr, fieldname)

  return f_ptr
end

function countFields(m_ptr::Ptr{Void})

  n = ccall( (:countFields, pumi_libname), Cint, (Ptr{Void},), m_ptr)

  return n
end

function getField(m_ptr::Ptr{Void}, i::Integer)

  f_ptr = ccall( (:getField, pumi_libname), Ptr{Void}, (Ptr{Void}, Cint), m_ptr, i)

  return f_ptr
end



function destroyField(f_ptr::Ptr{Void})

  ccall( (destroyField_name, pumi_libname), Void, (Ptr{Void},), f_ptr)

  return nothing
end

"""
  Destroys all :Fields associated with a given mesh, except
  those specified

  **Inputs**

   * m_ptr: :Mesh*
   * n_save: array of :Field* objects to not destroy (optional)
"""
function destroyFields(m_ptr::Ptr{Void}, n_save=Array{Ptr{Void}}(0))

  ccall( (destroyFields_name, pumi_libname), Void, (Ptr{Void}, Ptr{Ptr{Void}}, Cint), m_ptr, n_save, length(n_save))

  return nothing
end


@doc """
###PumiInterface.createSubMesh

  This function takes in a (presumable high order) mesh, subtriangulates 
  each element into a set of linear triangles, writes a vtk file with the
  name newmesh_linear, and returns a :Mesh2* pointer as a Void pointer.
  Note that Pumi is fully capabable of having multiple meshes loaded at one
  time (although only one of them can reference geometry).  The Julia 
  interface to Pumi does not always handle multiple meshes correctly.  In
  general, any function that takes in a mesh pointer will work correctly,
  but some of the ones that do not take in a mesh pointer will not.  This
  deficiency needs to be corrected eventually.

  It is useful to note that vertices of the original mesh that are also
  vertices of the new mesh will have the same number assigned to them by
  Pumi.  This makes it easier to compare ParaView images of the meshes
  side-by-side.

  Argument:
    m_ptr : pointer to existing mesh
    triangulation : n x 3 array of Int32s that tell which nodes of an existing
                    element should be used as the vertices of each subtriangle,
		    where n is the number of subtriangles.  These should be
		    1-based indices, with values in the range 1 to number of
		    nodes per element
    elementNodeOffsets:  array of UInt8s, dimensions number of nodes per element 
                         by number of elements in the mesh, used to figure out
			 how each element should access the nodes of a mesh 
			 entity.  See PdePumiInterface for details.
    typeOffsetsPerElement: array of Int32s, of length 4, where entry i
                           tells where the nodes on entities of dimension
			   i start in an array of all the nodes of an element.
			   The last element should be the number of nodes on 
			   an element + 1
    numberings:  array of length 3 containings numbers of all the verts, edges
                 and faces in the mesh m_ptr.

  Outputs:
    mnew_ptr : pointer to the new mesh
"""->
function createSubMesh(m_ptr, triangulation::AbstractArray{Int32, 2}, elementNodeOffsets::AbstractArray{UInt8, 2}, typeOffsetsPerElement::AbstractArray{Int32, 1}, numberings::AbstractArray{Ptr{Void}, 1})

  # check the the triangulation array is oriented correctly
  @assert size(triangulation, 1) == 3
  
 mnew_ptr =  ccall( (createSubMesh_name, pumi_libname), Ptr{Void}, (Ptr{Void}, Int32, Ptr{Int32}, Ptr{UInt8}, Ptr{Int32}, Ptr{Ptr{Void}}), m_ptr, size(triangulation, 2), triangulation, elementNodeOffsets, typeOffsetsPerElement, numberings)

 return mnew_ptr
end


# DG version of the above functions

function createSubMeshDG(m_ptr, mshape_ptr, triangulation::AbstractArray{Int32, 2}, elementNodeOffsets::AbstractArray{UInt8, 2}, typeOffsetsPerElement::AbstractArray{Int32, 1}, nodemapPumiToSbp::Array{UInt8, 1}, numberings::AbstractArray{Ptr{Void}, 1}, coords::AbstractArray{Float64, 3})

  # check the the triangulation array is oriented correctly
  @assert size(triangulation, 1) == 3
  
 mnew_ptr =  ccall( (createSubMeshDG_name, pumi_libname), Ptr{Void}, (Ptr{Void}, Ptr{Void}, Int32, Ptr{Int32}, Ptr{UInt8}, Ptr{Int32}, Ptr{UInt8}, Ptr{Ptr{Void}}, Ptr{Float64}), m_ptr, mshape_ptr, size(triangulation, 2), triangulation, elementNodeOffsets, typeOffsetsPerElement, nodemapPumiToSbp, numberings, coords)

 return mnew_ptr
end



function transferFieldDG(m_ptr, mnew_ptr, triangulation::AbstractArray{Int32, 2}, elementNodeOffsets::AbstractArray{UInt8, 2}, typeOffsetsPerElement::AbstractArray{Int32, 1}, numberings::AbstractArray{Ptr{Void}, 1}, field_old, interp_op::AbstractArray{Float64, 2}, field_new)

  # check the the triangulation array is oriented correctly
  @assert size(triangulation, 1) == 3
  @assert size(interp_op, 2) == 3  # the interpolation operator must be 
                                   # 3 x numnodesperlement in row major land
                                   # so, make sure it is numnodesperelement x 3
  
 ccall( (transferFieldDG_name, pumi_libname), Void, (Ptr{Void}, Ptr{Void}, Int32, Ptr{Int32}, Ptr{UInt8}, Ptr{Int32}, Ptr{Ptr{Void}}, Ptr{Cdouble}, Ptr{Void}, Ptr{Void}), m_ptr, mnew_ptr, size(triangulation, 2), triangulation, elementNodeOffsets, typeOffsetsPerElement, numberings, interp_op, field_old, field_new)

 return nothing
end


@doc """
### PumiInterface.transferField

  Transfers the specified field from the hold mesh to the new mesh.

  See createSubMesh for the meanings of the arguments
"""->
function transferField(m_ptr, mnew_ptr, triangulation::AbstractArray{Int32, 2}, elementNodeOffsets::AbstractArray{UInt8, 2}, typeOffsetsPerElement::AbstractArray{Int32, 1}, numberings::AbstractArray{Ptr{Void}, 1}, field_old, field_new)

  # check the the triangulation array is oriented correctly
  @assert size(triangulation, 1) == 3
#  @assert size(interp_op, 2) == 3  # the interpolation operator must be 
                                   # 3 x numnodesperlement in row major land
                                   # so, make sure it is numnodesperelement x 3
  
 ccall( (transferField_name, pumi_libname), Void, (Ptr{Void}, Ptr{Void}, Int32, Ptr{Int32}, Ptr{UInt8}, Ptr{Int32}, Ptr{Ptr{Void}}, Ptr{Void}, Ptr{Void}), m_ptr, mnew_ptr, size(triangulation, 2), triangulation, elementNodeOffsets, typeOffsetsPerElement, numberings, field_old, field_new)

 return nothing
end



function getFieldShape(shape_type::Integer, order::Integer, dim::Integer)
# shape_type: 0 = Lagrange, 1 = SBP, 2 = SBPDG1

  change_shape = Ref{Bool}()
  ccall( (getFieldShape_name, pumi_libname), Ptr{Void}, (Cint, Cint, Cint, Ref{Bool}), shape_type, order, dim, change_shape)

end


function countPeers(m_ptr, dim)
# get the number of other parts that have the entities of the specified
# dimension in common with this part

 npeers = ccall( (countPeers_name, pumi_libname), Csize_t, (Ptr{Void}, Cint), m_ptr, dim)

 return npeers
end

function getPeers(m_ptr, part_nums::Array{Cint})

  ccall( (getPeers_name, pumi_libname), Void, (Ptr{Void}, Ptr{Cint}), m_ptr, part_nums)

  return nothing
end

function isShared(m_ptr, e_ptr)

  val = ccall( (isShared_name, pumi_libname), Cint, (Ptr{Void}, Ptr{Void}), m_ptr, e_ptr)
  return val != 0
end

function countRemotes(m_ptr, entity)

  nremotes = ccall( (countRemotes_name, pumi_libname), Csize_t, (Ptr{Void}, Ptr{Void}), m_ptr, entity)

  return nremotes
end

function getRemotes(partnums::Array{Cint}, entities::Array{Ptr{Void}})
  
  ccall( (getRemotes_name, pumi_libname), Void, (Ptr{Cint}, Ptr{Ptr{Void}}), partnums, entities)

  return nothing
end

function setPoint(m_ptr, entity, node, coords::AbstractArray{Float64})
# coords must be of length 3, even in 2D
  @assert length(coords) == 3

  ccall((setPoint_name, pumi_libname), Void, (Ptr{Void}, Ptr{Void}, Cint, Ptr{Float64}), m_ptr, entity, node, coords)

  return nothing
end

function setParam(m_ptr, entity, coords::AbstractArray{Float64})
  @assert length(coords) >= 2

  ccall((:setParam, pumi_libname), Void, (Ptr{Void}, Ptr{Void}, Ptr{Float64}), m_ptr, entity, coords)

  return nothing
end



function acceptChanges(m_ptr)

  ccall( (acceptChanges_name, pumi_libname), Void, (Ptr{Void},), m_ptr)
  return nothing
end

function Verify(m_ptr)

  ccall( (Verify_name, pumi_libname), Void, (Ptr{Void},), m_ptr)
  return nothing
end

function getPoint(m_ptr, entity, node, coords::AbstractArray{Float64})
# coords must be of length 3, even in 2D
  @assert length(coords) == 3

  ccall((getPoint_name, pumi_libname), Void, (Ptr{Void}, Ptr{Void}, Cint, Ptr{Float64}), m_ptr, entity, node, coords)

  return nothing
end

function getParam(m_ptr, entity, coords::AbstractArray{Float64})
  @assert length(coords) >= 2

  ccall((:getParam, pumi_libname), Void, (Ptr{Void}, Ptr{Void}, Ptr{Float64}), m_ptr, entity, coords)

  return nothing
end



function hasMatching(m_ptr::Ptr{Void})

  val = ccall( (hasMatching_name, pumi_libname), CppBool, (Ptr{Void},), m_ptr)

  return val != 0
end

function getSharing(m_ptr)

  shr = ccall( (getSharing_name, pumi_libname), Ptr{Void}, (Ptr{Void},), m_ptr)

  return shr
end


function getNormalSharing(m_ptr)

  shr = ccall( (getNormalSharing_name, pumi_libname), Ptr{Void}, (Ptr{Void},), m_ptr)

  return shr
end


function freeSharing(shr_ptr::Ptr{Void})

  val = ccall( (freeSharing_name, pumi_libname), Void, (Ptr{Void},), shr_ptr)

end

function isOwned(shr_ptr, entity)

  val = ccall( (isOwned_name, pumi_libname), CppBool, (Ptr{Void}, Ptr{Void}), shr_ptr, entity)

  return val != 0
end

function countCopies(shr_ptr, entity)

  val = ccall( (countCopies_name, pumi_libname), Csize_t, (Ptr{Void}, Ptr{Void}), shr_ptr, entity)

  return val
end

function getCopies(part_nums::AbstractArray{Cint}, entities::AbstractArray{Ptr{Void}})
# as usual, the arrays must be the right length, as defined by countCopies
# the entities returned correspond to the entity passed in from the most recent
# call to countCopies

  ccall( (getCopies_name, pumi_libname), Void, (Ptr{Cint}, Ptr{Ptr{Void}}), part_nums, entities)

end


function getOwner(shr_ptr, entity)

  val = ccall( (getOwner_name, pumi_libname), Cint, (Ptr{Void}, Ptr{Void}), shr_ptr, entity)

  return val
end


function isSharedShr(shr_ptr, entity)

  val = ccall( (isSharedShr_name, pumi_libname), CppBool, (Ptr{Void}, Ptr{Void}), shr_ptr, entity)

  return val != 0
end



function countMatches(m_ptr, entity)

  val = ccall( (countMatches_name, pumi_libname), Csize_t, (Ptr{Void}, Ptr{Void}), m_ptr, entity)

  return val
end

function getMatches(part_nums::AbstractArray{Cint}, entities::AbstractArray{Ptr{Void}})
# as usual, the arrays must be the right length, as defined by countMatches
# the entities returned correspond to the entity passed in from the most recent
# call to countMatches

  ccall( (getMatches_name, pumi_libname), Void, (Ptr{Cint}, Ptr{Ptr{Void}}), part_nums, entities)

end


"""
  Get information about the Pumi reference element
"""
function getTopologyMaps()

  tri_edge_verts = Array{Cint}(3, 2)
  tet_edge_verts = Array{Cint}(6, 2)
  tet_tri_verts = Array{Cint}(4, 3)

  ccall( (:getTopologyMaps, pumi_libname), Void, (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), tri_edge_verts, tet_edge_verts, tet_tri_verts)

  return tri_edge_verts, tet_edge_verts, tet_tri_verts
end


"""
  Create a submesh, returning an object that contains information about the
  relation between the elements of the two meshes

  **Inputs**

   * m_ptr: :Mesh* for existing mesh
   * numberings: array of 0-based :Numberings, of length mesh.dim + 1, numbering
     the entities of dimension 0 to dim
   * el_list: array of Cint element numbers (in the numberings[end] numbering)
              that will exist on the new mesh

  **Output**
  
   * SubMeshData: type containing information relating the old mesh and the
                  new mesh
"""
function createSubMesh(m_ptr::Ptr{Void}, numberings::AbstractArray{Ptr{Void}},
                       el_list::AbstractArray{Cint})

  sdata = ccall( (:createSubMesh2, pumi_libname), Ptr{Void}, (Ptr{Void}, Ptr{Ptr{Void}}, Ptr{Cint}, Cint), m_ptr, numberings, el_list, length(el_list))

  return SubMeshData(sdata)
end

"""
  Writes the new mesh (created by [`createSubMesh`](@ref) to a file

  **Inputs**

   * sdata: SubMeshData
   * fname: file name, including .smb extension

"""
function writeNewMesh(sdata::SubMeshData, fname::AbstractString)

  ccall( (:writeNewMesh, pumi_libname), Void, (SubMeshData, Cstring), sdata, fname)
end

"""
  Get the :Mesh* of the new mesh

  **Inputs**

   * sdata: SubMeshData

  **Outputs**

   * a :Mesh*
"""
function getNewMesh(sdata::SubMeshData)

  m_ptr = ccall( (:getNewMesh, pumi_libname), Ptr{Void}, (SubMeshData,), sdata)
  return m_ptr
end

"""
  Get the original mesh pointer.
"""
function getOldMesh(sdata::SubMeshData)

  m_ptr = ccall( (:getOldMesh, pumi_libname), Ptr{Void}, (SubMeshData,), sdata)
  return m_ptr
end


"""
  Returns an :Numbering* for a Numbering on the submesh containing the
  element number on the original mesh that each element came from.
"""
function getParentNumbering(sdata::SubMeshData)

  n_ptr = ccall( (:getParentNumbering, pumi_libname), Ptr{Void}, (SubMeshData,), sdata)
  return n_ptr
end

function getGeoTag(sdata::SubMeshData)

  new_geo = ccall( (:getGeoTag, pumi_libname), Cint, (SubMeshData,), sdata)
  return new_geo
end


function free(sdata::SubMeshData)

  n_ptr = ccall( (:freeSubMesh2, pumi_libname), Void, (SubMeshData,), sdata)
  return n_ptr
end

function printTags(m_ptr::Ptr{Void})

  ccall( (:printTags, pumi_libname), Void, (Ptr{Void},), m_ptr)
end


include("apf2.jl")  # higher level functions
end  # end of module
