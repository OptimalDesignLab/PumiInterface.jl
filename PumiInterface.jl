# functions to test the julia/PUMI interface
module PumiInterface
include("PumiInterface2.jl")  # higher level functions

# no names should exported because there should be higher level functions
# wrapping these
# but they are going to be exported anyways
@doc """
  declareNames() declares global constant variables that are used to construct
  the first argument to the ccall function (it must be a constant expression)

  names declared in C++ code as 'extern "C" ' will not be mangled, but some
  libraries do not do this, so these variables can be used to declare a human
  readable variable name for the mangled name
"""

function declareNames()
# declare variables that hold the (possible mangled) names of c++ library functions
global const pumi_libname = "libfuncs1"
global const init_name = "initABC"
global const init2_name = "initABC2"
global const getMeshPtr_name = "getMeshPtr"
global const getConstantShapePtr_name = "getConstantShapePtr"
global const getMeshShapePtr_name = "getMeshShapePtr"
global const getVertNumbering_name = "getVertNumbering"
global const getEdgeNumbering_name = "getEdgeNumbering"
global const getFaceNumbering_name = "getFaceNumbering"
global const getElNumbering_name = "getElNumbering"

global const resetVertIt_name = "resetVertIt"
global const resetEdgeIt_name = "resetEdgeIt"
global const resetFaceIt_name = "resetFaceIt"
global const resetElIt_name  = "resetElIt"
global const incrementVertIt_name = "incrementVertIt"
global const incrementVertItn_name = "incrementVertItn"

global const incrementEdgeIt_name = "incrementEdgeIt"
global const incrementEdgeItn_name = "incrementEdgeItn"

global const incrementFaceIt_name = "incrementFaceIt"
global const incrementFaceItn_name = "incrementFaceItn"

global const incrementElIt_name = "incrementElIt"
global const incrementElItn_name = "incrementElItn"

global const count_name = "count"
global const writeVtkFiles_name = "writeVtkFiles"

global const setGlobalVertNumber_name = "setGlobalVertNumber"
global const getGlobalVertNumber_name = "getGlobalVertNumber"
global const getVertNumber_name = "getVertNumber"
global const getEdgeNumber_name = "getEdgeNumber"
global const getFaceNumber_name = "getFaceNumber"
global const getElNumber_name = "getElNumber"
global const getVert_name = "getVert"
global const getEdge_name = "getEdge"
global const getFace_name = "getFace"
global const getEl_name = "getEl"
global const getVertNumber2_name = "getVertNumber2"
global const getEdgeNumber2_name = "getEdgeNumber2"
global const getFaceNumber2_name = "getFaceNumber2"
global const getElNumber2_name = "getElNumber2"
global const getMeshDimension_name = "getMeshDimension"
global const getType_name = "getType"
global const getDownward_name = "getDownward"
global const countAdjacent_name = "countAdjacent"
global const getAdjacent_name = "getAdjacent"

global const hasNodesIn_name = "hasNodesIn"
global const countNodesOn_name = "countNodesOn"
global const getEntityShape_name = "getEntityShape"

global const createMeshElement_name = "createMeshElement"
global const countIntPoints_name = "countIntPoints"
global const getIntPoint_name = "getIntPoint"
global const getIntWeight_name = "getIntWeight"
global const getJacobian_name = "getJacobian"

global const countNodes_name = "countNodes"
global const getValues_name = "getValues"
global const getLocalGradients_name = "getLocalGradients"

global const checkVars_name = "checkVars"
global const checkNums_name = "checkNums"
global const getVertCoords_name = "getVertCoords"
global const getVertCoords2_name = "getVertCoords2"
global const getEdgeCoords_name = "getEdgeCoords"
global const getEdgeCoords2_name = "getEdgeCoords2"
global const getFaceCoords_name = "getFaceCoords"
global const getFaceCoords2_name = "getFaceCoords2"
global const getElCoords_name = "getElCoords"
global const createNumberingJ_name = "createNumberingJ"
global const numberJ_name = "numberJ"
global const getNumberJ_name = "getNumberJ"
global const getMesh_name = "getMesh"
global const printNumberingName_name = "printNumberingName"

global const createDoubleTag_name = "createDoubleTag"
global const setDoubleTag_name = "setDoubleTag"
global const getDoubleTag_name = "getDoubleTag"


global const reorder_name = "reorder"
end


# export low level interface functions
export declareNames, init, init2, getMeshPtr, getConstantShapePtr, getMeshShapePtr, getVertNumbering, getEdgeNumbering, getFaceNumbering, getElNumbering, resetVertIt, resetEdgeIt, resetFaceIt, resetElIt, incrementVertIt, incrementVertItn, incrementEdgeIt, incrementEdgeItn, incrementFaceIt, incrementFaceItn, incrementElIt, incrementElItn, countJ, writeVtkFiles, getVertNumber, getEdgeNumber, getFaceNumber, getElNumber, getVert, getEdge, getFace, getEl, getVertNumber2, getEdgeNumber2, getElNumber2, getMeshDimension, getType, getDownward, countAdjacent, getAdjacent, hasNodesIn, countNodesOn, getEntityShape, createMeshElement, countIntPoints, getIntPoint, getIntWeight, getJacobian, countNodes, getValues, getLocalGradients, checkVars, checkNums, getVertCoords, getEdgeCoords, getFaceCoords, getElCoords, createNumberingJ, numberJ, getNumberJ, getMesh, printNumberingName, createDoubleTag, setDoubleTag, getDoubleTag, reorder

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
    * m_ptr :        pointer to apf::Mesh2 that was loaded from the smb file
    * mshape_ptr :   pointer to the apf::FieldShape of the mesh
"""

function init(dmg_name::AbstractString, smb_name::AbstractString)
# initilize mesh interface
# initilize pointers to some value
# this is hack-ish -- there should be a better way to do this
#m_ptr = Ptr{Void}
#mshape_ptr = Ptr{Void}
downward_counts = zeros(Int32, 4,4);
num_Entities = zeros(Int32, 4, 1)

m_ptr_array = Array(Ptr{Void}, 1)
mshape_ptr_array = Array(Ptr{Void}, 1)

i = ccall( (init_name, pumi_libname), Int32, (Ptr{UInt8}, Ptr{UInt8}, Ptr{Int32},Ptr{Int32}, Ptr{Void}, Ptr{Void}), dmg_name, smb_name, downward_counts, num_Entities, m_ptr_array, mshape_ptr_array )  # call init in interface library

if ( i != 0)
  println("init failed, exiting ...")
  exit()
end


return downward_counts, num_Entities, m_ptr_array[1], mshape_ptr_array[1]
end


# 2d initilization
function init2(dmg_name::AbstractString, smb_name::AbstractString, order::Integer)
# initilize mesh interface
# initilize pointers to some value
# order = order of shape functions to use, currently only necessary because
#      pumi does not properly support saving shape functions to files

# this is hack-ish -- there should be a better way to do this
#m_ptr = Ptr{Void}
#mshape_ptr = Ptr{Void}
downward_counts = zeros(Int32, 3,3);
num_Entities = zeros(Int32, 4, 1)

m_ptr_array = Array(Ptr{Void}, 1)
mshape_ptr_array = Array(Ptr{Void}, 1)

i = ccall( (init2_name, pumi_libname), Int32, (Ptr{UInt8}, Ptr{UInt8}, Ptr{Int32},Ptr{Int32}, Ptr{Void}, Ptr{Void}, Int32), dmg_name, smb_name, downward_counts, num_Entities, m_ptr_array, mshape_ptr_array, order )  # call init in interface library

if ( i != 0)
  println("init failed, exiting ...")
  exit()
end


return downward_counts, num_Entities, m_ptr_array[1], mshape_ptr_array[1]
end





# no longer needed
function getMeshPtr()

  m_ptr = ccall( (getMeshPtr_name, pumi_libname), Ptr{Void}, () )
  return m_ptr
end

# no longer needed
function getMeshShapePtr()

  m_ptr = ccall( (getMeshShapePtr_name, pumi_libname), Ptr{Void}, () )
  return m_ptr
end

function getConstantShapePtr(dimension::Integer)
  mshape_ptr = ccall( ( getConstantShapePtr_name, pumi_libname), Ptr{Void}, (Int32,), dimension)

  return mshape_ptr

end


function getVertNumbering()

  numbering_ptr = ccall( (getVertNumbering_name, pumi_libname), Ptr{Void}, () )
  return numbering_ptr
end

function getEdgeNumbering()

  numbering_ptr = ccall( (getEdgeNumbering_name, pumi_libname), Ptr{Void}, () )
  return numbering_ptr
end

function getFaceNumbering()

  numbering_ptr = ccall( (getFaceNumbering_name, pumi_libname), Ptr{Void}, () )
  return numbering_ptr
end

function getElNumbering()

  numbering_ptr = ccall( (getElNumbering_name, pumi_libname), Ptr{Void}, () )
  return numbering_ptr
end

function resetVertIt()
# reset the vertex iterator to the beginning

ccall( (resetVertIt_name, pumi_libname), Void, () );
return nothing

end


function resetEdgeIt()
# reset the edge iterator to the beginning

ccall( (resetEdgeIt_name, pumi_libname), Void, () );
return nothing

end


function resetFaceIt()
# reset the face iterator to the beginning

ccall( (resetFaceIt_name, pumi_libname), Void, () );
return nothing

end


function resetElIt()
# reset the element iterator to the beginning

ccall( (resetElIt_name, pumi_libname), Void, () );
return nothing

end

function incrementVertIt()
# increment it vertex iterator

ccall( (incrementVertIt_name, pumi_libname), Void, () );
return nothing

end


function incrementVertItn(n::Integer)
# increment it vertex iterator

ccall( (incrementVertItn_name, pumi_libname), Void, (Int32,), n );
return nothing

end




function incrementEdgeIt()
# increment it edge iterator

ccall( (incrementEdgeIt_name, pumi_libname), Void, () );
return nothing

end

function incrementEdgeItn(n::Integer)
# increment it vertex iterator

ccall( (incrementEdgeItn_name, pumi_libname), Void, (Int32,), n );
return nothing

end




function incrementFaceIt()
# increment it Face iterator

ccall( (incrementFaceIt_name, pumi_libname), Void, () );
return nothing

end

function incrementFaceItn(n::Integer)
# increment it face iterator n times

ccall( (incrementFaceItn_name, pumi_libname), Void, (Int32,), n );
return nothing

end




function incrementElIt()
# increment it element iterator

ccall( (incrementElIt_name, pumi_libname), Void, () );
return nothing

end

function incrementElItn(n::Integer)
# increment it element iterator

ccall( (incrementElItn_name, pumi_libname), Void, (Int32,), n );
return nothing

end


function countJ(m_ptr, dimension::Integer)
# returns the number of entities of a particular dimension

  i = ccall( (count_name, pumi_libname), Int32, (Ptr{Void},Int32), m_ptr, dimension)
  return i
end

function writeVtkFiles(name::AbstractString, m_ptr)
# write vtk files to be read by paraview

  ccall( (writeVtkFiles_name, pumi_libname), Void, (Ptr{UInt8}, Ptr{Void}), name, m_ptr)
  return nothing
end

function setGlobalVertNumber(val::Integer)
# set global Vertex number of the current node

ccall( (setGlobalVertNumber_name, pumi_libname), Void, (Int32,), val)
return nothing

end

function getGlobalVertNumber()
# get global vertex number of the current

i = ccall( (getGlobalVertNumber_name, pumi_libname), Int32, () )
return i

end


function getVertNumber()
# get the number of the current vertex
  i = ccall( (getVertNumber_name, pumi_libname), Int32, () )
  return i
end


function getEdgeNumber()
# get the number of the current edge
  i = ccall( (getEdgeNumber_name, pumi_libname), Int32, () )
  return i
end


function getFaceNumber()
# get the number of the current face
  i = ccall( (getFaceNumber_name, pumi_libname), Int32, () )
  return i
end

function getElNumber()
# get the number of the current element
  i = ccall( (getElNumber_name, pumi_libname), Int32, () )
  return i
end


function getVert()
# get pointer to current vertex
  entity = ccall( (getVert_name, pumi_libname), Ptr{Void}, () )

  return entity
end

function getEdge()
# get pointer to current Edge
  entity = ccall( (getEdge_name, pumi_libname), Ptr{Void}, () )

  return entity
end

function getFace()
# get pointer to current face
  entity = ccall( (getFace_name, pumi_libname), Ptr{Void}, () )

  return entity
end

function getEl()
# get pointer to current element
  entity = ccall( (getEl_name, pumi_libname), Ptr{Void}, () )

  return entity
end




function getVertNumber2(entity)

  i = ccall( (getVertNumber2_name, pumi_libname), Int32, (Ptr{Void},), entity)
  return i
end


function getEdgeNumber2(entity)

  i = ccall( (getEdgeNumber2_name, pumi_libname), Int32, (Ptr{Void},), entity)
  return i
end


function getFaceNumber2(entity)

  i = ccall( (getFaceNumber2_name, pumi_libname), Int32, (Ptr{Void},), entity)
  return i
end


function getElNumber2(entity)

  i = ccall( (getElNumber2_name, pumi_libname), Int32, (Ptr{Void},), entity)
  return i
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
  downwards = Array(Ptr{Void}, 12)


  i = ccall ( (getDownward_name, pumi_libname), Int32, (Ptr{Void}, Ptr{Void}, Int32, Ptr{Void}), m_ptr, entity, dimension, downwards)

  return downwards[1:i], i

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
  adjacencies_ret = Array(Ptr{Void}, num_adjacent)  # create the array

  ccall( (getAdjacent_name, pumi_libname), Void, (Ptr{Void},), adjacencies_ret)

  return adjacencies_ret
end


function hasNodesIn(mshape_ptr, dimension::Integer)
# check whether this FieldShape* has nodes on entities of a given dimension

  i = ccall ( (hasNodesIn_name, pumi_libname), Cuchar, (Ptr{Void}, Int32), mshape_ptr, dimension)

  i_bool = convert(Bool, i)
#  println("hasNodesIn returned", i_bool)

  return i_bool

end

function countNodesOn(mshape_ptr, entity_type::Integer)
# count the number of nodes on an entity of the specified type (apf::Mesh::Type)
  i = ccall ( (countNodesOn_name, pumi_libname), Int32, (Ptr{Void}, Int32), mshape_ptr, entity_type)

  return i

end

function getEntityShape(mshape_ptr, entity_type::Integer)
# get the EntityShape* (object describing shape functions) ofa given entity type

  eshape_ptr = ccall ( (getEntityShape_name, pumi_libname), Ptr{Void}, (Ptr{Void}, Int32), mshape_ptr, entity_type)

  return eshape_ptr
end


function createMeshElement(m_ptr, entity)
# creates a MeshElement from a MeshEntity

  mel_ptr = ccall ( ( createMeshElement_name, pumi_libname), Ptr{Void}, (Ptr{Void}, Ptr{Void}), m_ptr, entity)
  return mel_ptr
end

function countIntPoints(mel_ptr, order::Integer)
# count the number of integration points needed to achieve specified order of
# accuracy.  MeshElement can be edges as well as 2D and 3D regions.
# not sure what happens if it is a vertex.

  i = ccall ( (countIntPoints_name, pumi_libname), Int32, (Ptr{Void}, Int32), mel_ptr, order)
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

  i = ccall ( (countNodes_name, pumi_libname), Int32, (Ptr{Void},), eshape_ptr)
  return i
end

function getValues( eshape_ptr, coords::Array{Float64,1}, numN::Integer)
#  gets an array of shape function values at the specified coordinates
# numN is the number of points affecting the entity that was used to get eshape_ptr
# coords must be a vector of length 3
# the output is a vector of length numN

  vals = zeros(numN)
  ccall ( (getValues_name, pumi_libname), Void, (Ptr{Void}, Ptr{Float64}, Ptr{Float64}), eshape_ptr, coords, vals)

  return vals
end


function getLocalGradients( eshape_ptr, coords::Array{Float64,1}, numN::Integer)
#  gets an array of shape function derivatives at the specified coordinates
# numN is the number of points affecting the entity that was used to get eshape_ptr
# coords must be a vector of length 3
# the output is a matrix of dimension 3 x numN

  vals = zeros(3, numN)
  ccall ( (getLocalGradients_name, pumi_libname), Void, (Ptr{Void}, Ptr{Float64}, Ptr{Float64}), eshape_ptr, coords, vals)

  return vals
end

function checkVars()
ccall ( (checkVars_name, pumi_libname), Void, () );

return nothing
end






function checkVars()
ccall ( (checkVars_name, pumi_libname), Void, () );

return nothing
end


function checkNums()
ccall( (checkNums_name, pumi_libname), Void, () )

end

function getVertCoords(coords::Array{Float64, 2}, m::Integer, n::Integer)
# coords is array to put coordsinates in, must be 3 by 1,
# m, n are number of rows, columns in coords, respectively

#coords = Array(Float64, 3, 2)   # pass an array 3 by n (3 coordinates each for n points)
#(m,n) = size(coords)
# pass reversed m,n because C arrays are row-major
ccall( (getVertCoords_name, pumi_libname), Void, (Ptr{Float64}, Int32, Int32), coords, n, m) 

#println("\n from julia, coords = ", coords)

return coords
end


function getVertCoords(entity, coords::Array{Float64, 2}, m::Integer, n::Integer)
# coords is array to put coordsinates in, must be 3 by 1,
# m, n are number of rows, columns in coords, respectively

#coords = Array(Float64, 3, 2)   # pass an array 3 by n (3 coordinates each for n points)
#(m,n) = size(coords)
# pass reversed m,n because C arrays are row-major
ccall( (getVertCoords2_name, pumi_libname), Void, (Ptr{Void}, Ptr{Float64}, Int32, Int32), entity, coords, n, m) 

#println("\n from julia, coords = ", coords)

return coords
end



function getEdgeCoords(coords::Array{Float64, 2}, m::Integer, n::Integer)
# coords is array to put coordinates in, must be 3 by 2
# m,n = number of rows, columns in coords, respectively

i = ccall( (getEdgeCoords_name, pumi_libname), Int, (Ptr{Float64}, Int32, Int32), coords, n, m);

if ( i != 0)
  println("Error in getEdgeCoords... exiting")
  exit()
end

#println("in julia, coords = ", coords)

end

function getEdgeCoords(entity, coords::Array{Float64, 2}, m::Integer, n::Integer)
# coords is array to put coordinates in, must be 3 by 2
# m,n = number of rows, columns in coords, respectively

i = ccall( (getEdgeCoords2_name, pumi_libname), Int, (Ptr{Void}, Ptr{Float64}, Int32, Int32), entity, coords, n, m);

if ( i != 0)
  println("Error in getEdgeCoords... exiting")
  exit()
end

#println("in julia, coords = ", coords)

end





function getFaceCoords(coords::Array{Float64, 2}, m::Integer, n::Integer)
# get coordinates of points on a face, in order
# coords is array to populate with coordinates, must by 3 by number of points
# on a face
# m,n = number of rows, columns in coords, respectively

# reverse m and n because C is row major
i = ccall( (getFaceCoords_name, pumi_libname), Int32, (Ptr{Float64}, Int32, Int32), coords, n, m);

if ( i != 0)
  println("Error in getEdgeCoords... exiting")
  exit()
end

#println("in julia, coords = ", coords)

end


function getFaceCoords(entity, coords::Array{Float64, 2}, m::Integer, n::Integer)
# get coordinates of points on a face, in order
# coords is array to populate with coordinates, must by 3 by number of points
# on a face
# m,n = number of rows, columns in coords, respectively

# reverse m and n because C is row major
i = ccall( (getFaceCoords2_name, pumi_libname), Int32, (Ptr{Void}, Ptr{Float64}, Int32, Int32), entity, coords, n, m);

if ( i != 0)
  println("Error in getEdgeCoords... exiting")
  exit()
end

#println("in julia, coords = ", coords)

end


function getElCoords(coords::Array{Float64, 2}, m::Integer, n::Integer)
# get coordinates of points in an element, in order
# coords is array to populate with coordinates, must by 3 by number of points
# on a element
# m,n = number of rows, columns in coords, respectively

# reverse m and n because C is row major
i = ccall( (getElCoords_name, pumi_libname), Int, (Ptr{Float64}, Int32, Int32), coords, n, m);

if ( i != 0)
  println("Error in getEdgeCoords... exiting")
  exit()
end

#println("in julia, coords = ", coords)

end



function createNumberingJ(m_ptr, name::AbstractString, field, components::Integer)
# create a generally defined numbering, get a pointer to it
# this just passes through to apf::createNumbering

numbering_ptr = ccall ( (createNumberingJ_name, pumi_libname), Ptr{Void}, (Ptr{Void}, Ptr{UInt8}, Ptr{Void}, Int32), m_ptr, name, field, components)

return numbering_ptr
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


function getMesh(n_ptr)
# get the mesh a numbering is defined on

  m_ptr = ccall( (getMesh_name, pumi_libname), Ptr{Void}, (Ptr{Void},), n_ptr)
  return m_ptr
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
function reorder(m_ptr, ndof::Integer, nnodes::Integer, dof_statusN_ptr, nodeNums, elNums, els_reordered)

  nodeN_array = [nodeNums] # pass by reference, holds numbering of dofs
  elN_array = [elNums] # pass by reference, hold numbering of elements

  ccall( (reorder_name, pumi_libname), Void, (Ptr{Void}, Int32, Int32, Int32, Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Void}),  m_ptr, ndof, nnodes,  2, dof_statusN_ptr, nodeNums, elNums, els_reordered)

  flush_cstdio()
  return nothing

end



declareNames()  # will this execute when module is compiled?

end  # end of module
