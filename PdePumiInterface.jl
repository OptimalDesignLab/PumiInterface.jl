module PdePumiInterface
push!(LOAD_PATH, "/users/creanj/julialib_fork/PUMI.jl")

using PumiInterface
using SummationByParts

export AbstractMesh,PumiMesh2, reinitPumiMesh2, getElementVertCoords, getShapeFunctionOrder, getGlobalNodeNumber, getGlobalNodeNumbers, getNumEl, getNumEdges, getNumVerts, getNumNodes, getNumDofPerNode, getAdjacentEntityNums, getBoundaryEdgeNums, getBoundaryFaceNums, getBoundaryEdgeLocalNum, getEdgeLocalNum, getBoundaryArray, saveSolutionToMesh, retrieveSolutionFromMesh, retrieveNodeSolution, getAdjacentEntityNums, getNumBoundaryElements, getInterfaceArray

abstract AbstractMesh

type PumiMesh2 <: AbstractMesh   # 2d pumi mesh, triangle only
  m_ptr::Ptr{Void}  # pointer to mesh
  mshape_ptr::Ptr{Void} # pointer to mesh's FieldShape
  f_ptr::Ptr{Void} # pointer to apf::field for storing solution during mesh adaptation
  vert_Nptr::Ptr{Void}  # numbering of vertices
  edge_Nptr::Ptr{Void}  # numbering of edges
  el_Nptr::Ptr{Void}  # numbering of elements (faces)

  numVert::Int  # number of vertices in the mesh
  numEdge::Int # number of edges in the mesh
  numEl::Int  # number of elements (faces)
  order::Int # order of shape functions
  numDof::Int # number of degrees of freedom
  numNodes::Int  # number of nodes
  numDofPerNode::Int  # number of dofs per node
  numBoundaryEdges::Int # number of edges on the exterior boundary

  # hold pointers to mesh entities
  verts::Array{Ptr{Void},1}  # holds pointers to mesh vertices
  edges::Array{Ptr{Void},1}  # pointers to mesh edges
  elements::Array{Ptr{Void},1}  # pointers to faces

  dofnums_Nptr::Ptr{Void}  # pointer to Numbering of dofs (result of reordering)
  boundary_nums::Array{Int, 2}  # array of [element number, edgenumber] for each edge on the boundary

  
  
end

function PumiMesh2(dmg_name::AbstractString, smb_name::AbstractString, order; dofpernode=1)
  # construct pumi mesh by loading the files named
  # dmg_name = name of .dmg (geometry) file to load (use .null to load no file)
  # smb_name = name of .smb (mesh) file to load
  # order = order of shape functions
  # dofpernode = number of dof per node, default = 1


  tmp, num_Entities, m_ptr, mshape_ptr = init2(dmg_name, smb_name, order)
  f_ptr = createPackedField(m_ptr, "solution_field", dofpernode)

  numVert = convert(Int, num_Entities[1])
  numEdge =convert(Int,  num_Entities[2])
  numEl = convert(Int, num_Entities[3])

  # get pointers to mesh entity numberings
  vert_Nptr = getVertNumbering()
  edge_Nptr = getEdgeNumbering()
  el_Nptr = getFaceNumbering()

  verts = Array(Ptr{Void}, numVert)
  edges = Array(Ptr{Void}, numEdge)
  elements = Array(Ptr{Void}, numEl)
  dofnums_Nptr = createNumberingJ(m_ptr, "reordered dof numbers", mshape_ptr, dofpernode)  # 1 dof per node

  # get pointers to all MeshEntities
  # also initilize the field to zero
  resetAllIts2()
#  comps = zeros(dofpernode)
  comps = [1.0, 2, 3, 4]
  for i=1:numVert
    verts[i] = getVert()
    setComponents(f_ptr, verts[i], 0, comps)  # initilize the field
    incrementVertIt()
  end

  for i=1:numEdge
    edges[i] = getEdge()
    incrementEdgeIt()
  end

  for i=1:numEl
    elements[i] = getFace()
    incrementFaceIt()
  end

  resetAllIts2()
  println("performing initial numbering of dofs")
  # calculate number of nodes, dofs (works for first and second order)
  numnodes = order*numVert 
  numdof = numnodes*dofpernode
  # number dofs
  ctr= 1
  for i=1:numVert
    for j=1:dofpernode
      numberJ(dofnums_Nptr, verts[i], 0, j-1, ctr)
      println("vertex ", i,  " numbered ", ctr)
      ctr += 1
    end
  end

  if order >= 2
    for i=1:numEdges
      for j=1:dofpernode
        numberJ(dofnums_Nptr, edges[i], 0, j-1, ctr)
        ctr += 1
      end
    end
  end

 

  # count boundary edges
  bnd_edges_cnt = 0
  bnd_edges = Array(Int, numEdge, 2)
  for i=1:numEdge
    edge_i = getEdge()
    numFace = countAdjacent(m_ptr, edge_i, 2)  # should be count upward

    if numFace == 1  # if an exterior edge
      faces = getAdjacent(numFace)
      facenum = getFaceNumber2(faces[1]) + 1

      bnd_edges_cnt += 1
      bnd_edges[bnd_edges_cnt, 1] = facenum
      bnd_edges[bnd_edges_cnt, 2] = i
    end
    incrementEdgeIt()
  end

  bnd_edges_small = bnd_edges[1:bnd_edges_cnt, :]
#=
  println("typeof m_ptr = ", typeof(m_ptr))
  println("typeof mshape_ptr = ", typeof(mshape_ptr))
  println("typeof numVerg = ", typeof(numVert))
  println("typeof numEdge = ", typeof(numEdge))
  println("typeof numEl = ", typeof(numEl))
  println("typeof order = ", typeof(order))
  println("typeof numdof = ", typeof(numdof))
  println("typeof bnd_edges_cnt = ", typeof(bnd_edges_cnt))
  println("typeof verts = ", typeof(verts))
  println("typeof edges = ", typeof(edges))
  println("typeof element = ", typeof(elements))
  println("typeof dofnums_Nptr = ", typeof(dofnums_Nptr))
  println("typeof bnd_edges_small = ", typeof(bnd_edges_small))
=#
  println("numVert = ", numVert)
  println("numEdge = ", numEdge)
  println("numEl = ", numEl)
  println("numDof = ", numdof)
  println("numNodes = ", numnodes)



  writeVtkFiles("mesh_complete", m_ptr)
  return PumiMesh2(m_ptr, mshape_ptr, f_ptr, vert_Nptr, edge_Nptr, el_Nptr, numVert, numEdge, numEl, order, numdof, numnodes, dofpernode, bnd_edges_cnt, verts, edges, elements, dofnums_Nptr, bnd_edges_small)
end


# for reinitilizeing after mesh adaptation
function reinitPumiMesh2(mesh::PumiMesh2)
  # construct pumi mesh by loading the files named
  # dmg_name = name of .dmg (geometry) file to load (use .null to load no file)
  # smb_name = name of .smb (mesh) file to load
  # order = order of shape functions
  # dofpernode = number of dof per node, default = 1

  println("Reinitilizng PumiMesh2")

  # create random filenames because they are not used
  smb_name = "a"
  dmg_name = "b"
  order = mesh.order
  dofpernode = mesh.numDofPerNode
  tmp, num_Entities, m_ptr, mshape_ptr = init2(dmg_name, smb_name, order, false) # do not load new mesh
  f_ptr = mesh.f_ptr  # use existing solution field

  numVert = convert(Int, num_Entities[1])
  numEdge =convert(Int,  num_Entities[2])
  numEl = convert(Int, num_Entities[3])

  mesh.vert_Nptr = getVertNumbering()
  mesh.edge_Nptr = getEdgeNumbering()
  mesh.el_Nptr = getFaceNumbering()




  verts = Array(Ptr{Void}, numVert)
  edges = Array(Ptr{Void}, numEdge)
  elements = Array(Ptr{Void}, numEl)
  dofnums_Nptr = mesh.dofnums_Nptr  # use existing dof pointers

  # get pointers to all MeshEntities
  # also initilize the field to zero
  resetAllIts2()
#  comps = zeros(dofpernode)
  comps = [1.0, 2, 3, 4]
  for i=1:numVert
    verts[i] = getVert()
    incrementVertIt()
  end

  for i=1:numEdge
    edges[i] = getEdge()
    incrementEdgeIt()
  end

  for i=1:numEl
    elements[i] = getFace()
    incrementFaceIt()
  end

  resetAllIts2()
  println("performing initial numbering of dofs")
  # calculate number of nodes, dofs (works for first and second order)
  numnodes = order*numVert 
  numdof = numnodes*dofpernode
  # number dofs
  ctr= 1
  for i=1:numVert
    for j=1:dofpernode
      numberJ(dofnums_Nptr, verts[i], 0, j-1, ctr)
#      println("vertex ", i,  " numbered ", ctr)
      ctr += 1
    end
  end

  if order >= 2
    for i=1:numEdges
      for j=1:dofpernode
        numberJ(dofnums_Nptr, edges[i], 0, j-1, ctr)
        ctr += 1
      end
    end
  end

 

  # count boundary edges
  bnd_edges_cnt = 0
  bnd_edges = Array(Int, numEdge, 2)
  for i=1:numEdge
    edge_i = getEdge()
    numFace = countAdjacent(m_ptr, edge_i, 2)  # should be count upward

    if numFace == 1  # if an exterior edge
      faces = getAdjacent(numFace)
      facenum = getFaceNumber2(faces[1]) + 1

      bnd_edges_cnt += 1
      bnd_edges[bnd_edges_cnt, 1] = facenum
      bnd_edges[bnd_edges_cnt, 2] = i
    end
    incrementEdgeIt()
  end

  bnd_edges_small = bnd_edges[1:bnd_edges_cnt, :]

  # replace exising fields with new values
  mesh.numVert = numVert
  mesh.numEdge = numEdge
  mesh.numEl = numEl
  mesh.numDof = numdof
  mesh.numNodes= numnodes
  mesh.numBoundaryEdges = bnd_edges_cnt
  mesh.verts = verts  # does this need to be a deep copy?
  mesh.edges = edges
  mesh.elements = elements
  mesh.boundary_nums = bnd_edges_small

#=
  println("typeof m_ptr = ", typeof(m_ptr))
  println("typeof mshape_ptr = ", typeof(mshape_ptr))
  println("typeof numVerg = ", typeof(numVert))
  println("typeof numEdge = ", typeof(numEdge))
  println("typeof numEl = ", typeof(numEl))
  println("typeof order = ", typeof(order))
  println("typeof numdof = ", typeof(numdof))
  println("typeof bnd_edges_cnt = ", typeof(bnd_edges_cnt))
  println("typeof verts = ", typeof(verts))
  println("typeof edges = ", typeof(edges))
  println("typeof element = ", typeof(elements))
  println("typeof dofnums_Nptr = ", typeof(dofnums_Nptr))
  println("typeof bnd_edges_small = ", typeof(bnd_edges_small))
=#
  println("numVert = ", numVert)
  println("numEdge = ", numEdge)
  println("numEl = ", numEl)
  println("numDof = ", numdof)
  println("numNodes = ", numnodes)

  writeVtkFiles("mesh_complete", m_ptr)
end






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

function getShapeFunctionOrder(mesh::PumiMesh2)

return mesh.order
end


function getGlobalNodeNumber(mesh::PumiMesh2, el_num::Integer, local_node_num::Integer)
# el_num = element number
# local_node_num = vector of dof numbers on the specified node within element
# this only works for first and second order
# returns a vector with dof numbers of all dofs on the node


#println("el_num = ", el_num, " local_node_num = ", local_node_num)
el = mesh.elements[el_num]  # get element
node_entities = getNodeEntities(mesh.m_ptr, mesh.mshape_ptr, el)
#println("node_entities = ", node_entities)
node_entity = node_entities[local_node_num]
#println("node_entity = ", node_entity)
dofnums = zeros(Int, mesh.numDofPerNode)
for i=1:mesh.numDofPerNode
  dofnums[i] = getNumberJ(mesh.dofnums_Nptr, node_entities[local_node_num], 0, i-1)

#println("global dof number = ", number_i)
end

return dofnums

end

function getGlobalNodeNumbers(mesh::PumiMesh2, elnum::Integer)
# gets global node numbers of all dof on all nodes of the element
# this only works for first and second order
# output formap is array [numdofpernode, nnodes]  (each column contains dof numbers for a node)
el_i = mesh.elements[elnum]
type_i = getType(mesh.m_ptr, el_i)


nnodes = countAllNodes(mesh.mshape_ptr, type_i)
numdof = nnodes*mesh.numDofPerNode
node_entities = getNodeEntities(mesh.m_ptr, mesh.mshape_ptr, el_i)
dofnums = zeros(Int,mesh.numDofPerNode, nnodes)  # to be populated with global dof numbers
#println("nnodes = ", nnodes)
#println("dofpernode = ", mesh.numDofPerNode)
for i=1:nnodes
  for j=1:mesh.numDofPerNode
    dofnums[j, i] = getNumberJ(mesh.dofnums_Nptr, node_entities[i], 0, j-1)
  end
end

return dofnums


end

function getEdgeNumber(mesh::PumiMesh2, edge_num::Integer)
# get the number of an edge

  i = getNumberJ(mesh.edge_Nptr, mesh.edges[edge_num], 0, 0)

  return i

end


function getNumEl(mesh::PumiMesh2)
# returns the number of elements in the mesh

return mesh.numEl
end

function getNumEdges(mesh::PumiMesh2)
# retrns the number of edges in a mesh

return mesh.numEdge

end

function getNumVerts(mesh::PumiMesh2)
# returns number of vertices in a mesh

return mesh.numVert
end

function getNumNodes(mesh::PumiMesh2)
# returns total number of nodes in the mesh

return mesh.numDof
end

function getNumDofPerNode(mesh::PumiMesh2)
# get the number of dof on a node

return mesh.numDofPerNode

end

function getNumBoundaryEdges(mesh::PumiMesh2)
# return the number of edges on the boundary

  return mesh.numBoundaryEdges
end

function getNumBoundaryElements(mesh::PumiMesh2)
# count the number of elements on the boundary

   return length(unique(mesh.boundary_nums[:,1]))
end

function getBoundaryEdgeNums(mesh::PumiMesh2)
# get vector of edge numbers that are on boundary

edge_nums = mesh.boundary_nums
return edge_nums
end



# this doesn't exist for a 2D mesh
function getBoundaryFaceNums(mesh::PumiMesh2)
# get array of face numbers that are on boundary
# this only works in 2d
# this function doesn't work yet
#face_nums = zeros(Int, numEdges_on_boundary)
face_nums = zeros(Int, 0)
return face_nums
end


function getAdjacentEntityNums(mesh::PumiMesh2, entity_index::Integer, input_dimension::Integer, output_dimension::Integer)
# gets the numbers of the adjacent entities (upward or downward) of the given entity
# entity_index specifies the index of the input entity
# input_dimension specifies the dimension of the input entity (0 =vert, 1 = edge ...)
# output_dimension is the dimension of the adjacent  entities to be retrieved
# returns an array containing the indicies of the adjacent entities, and the number of them
# the length of the array might be greater than the number of adjacencies, so use the second returned value for the number of adjacencies


if (input_dimension == output_dimension)
  println("cannot get same level adjacencies")
end


array1 = [mesh.vert_Nptr, mesh.edge_Nptr, mesh.el_Nptr]
#array2 = [mesh.verts; mesh.edges; mesh.elements]  # array of arrays
array2 = Array(Array{Ptr{Void},1}, 3)
array2[1] = mesh.verts
array2[2] = mesh.edges
array2[3] = mesh.elements

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

  (down_edges, numedges) = getDownward(mesh.m_ptr, face, 1)
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

  (down_edges, numedges) = getDownward(mesh.m_ptr, element, 1)

  edgenum_local = 0
   for j = 1:3  # find which local edge is edge_i
    if down_edges[j] == edge
      edgenum_local = j
    end
  end

  return edgenum_local

end

function getBoundaryArray(mesh::PumiMesh2, bnd_array::Array{Boundary, 1})
# get an array of type Boundary for SBP
# creating an an array of a user defined type seems like a waste of memory operations
# bnd_array is a vector of type Boundary with length equal to the number of edges on the boundary of the mesh

#  bnd_array = Array(Boundary, mesh.numBoundaryEdges)

  for i=1:mesh.numBoundaryEdges
    facenum = mesh.boundary_nums[i,1]
    edgenum_global = mesh.boundary_nums[i,2]
    edgenum_local = getBoundaryEdgeLocalNum(mesh, edgenum_global)
    bnd_array[i] = Boundary(facenum, edgenum_local)
  end

  return nothing
end


function getInterfaceArray(mesh::PumiMesh2, interfaces::AbstractArray{Interface,1})
# get array of [elementL, elementR, edgeL, edgeR] for each internal edge,
# where elementL and R are the elements that use the edge, edgeL R are the
# local number of the edge within the element
# interfaces is the array to be populated with the data, it must be 
# number of internal edges by 1

  # only need internal boundaries (not external)
#  num_ext_edges = size(getBoundaryEdgeNums(mesh))[1]  # bad memory efficiency
#  num_int_edges = getNumEdges(mesh) - num_ext_edges

#  new_bndry = Boundary(2, 3)
#   println("new_bndry = ", new_bndry)

#  new_interface = Interface(1, 2, 3, 4)
#   println("new_interface = ", new_interface)

#   println("num_int_edges = ", num_int_edges)

#  interfaces = Array(typeof(new_interface), num_int_edges)

  pos = 1 # current position in interfaces
  for i=1:getNumEdges(mesh)
#     println("i = ", i)
#     println("pos = ", pos)
    # get number of elements using the edge
    adjacent_nums, num_adjacent = getAdjacentEntityNums(mesh, i, 1, 2)
#     println("num_adjacent = ", num_adjacent)
#     println("adjacent_nums = ", adjacent_nums)
    if num_adjacent > 1  # internal edge
#       println("this is an internal edge")
      element1 = adjacent_nums[1]
      element2 = adjacent_nums[2]

#      coords_1 = x[:, :, element1]
#      coords_2 = x[:, :, element2]

      coords_1 = zeros(3,3)
      coords_2 = zeros(3,3)
      getFaceCoords(mesh.elements[element1], coords_1, 3, 3)
      getFaceCoords(mesh.elements[element2], coords_2, 3, 3)

      # calculate centroid
      centroid1 = sum(coords_1, 2)
      centroid2 = sum(coords_2, 2)

      if abs(centroid1[1] - centroid2[2]) < 1e-10  # if big enough difference
        if centroid1[1] < centroid2[1]
    elementL = element1
    elementR = element2
        else
    elementL = element2
    elementR = element1
        end
      else  # use y coordinate to decide
        if centroid1[2] < centroid2[2]
    elementL = element1
    elementR = element2
        else
    elementL = element2
    elementR = element1
        end
      end

      edgeL = getEdgeLocalNum(mesh, i, elementL)
      edgeR = getEdgeLocalNum(mesh, i, elementR)

      interfaces[pos] = Interface(elementL, elementR, edgeL, edgeR)
#       println("updating pos")
      pos += 1

  #    print("\n")

    end  # end if internal edge

#     print("\n")
  end  # end loop over edges

  return nothing

end  # end function





#=
function getEdgeInterfaceData(i::Integer)

#     println("i = ", i)
#     println("pos = ", pos)
    # get number of elements using the edge
    adjacent_nums, num_adjacent = getAdjacentEntityNums(mesh, i, 1, 2)
#     println("num_adjacent = ", num_adjacent)
#     println("adjacent_nums = ", adjacent_nums)
    if num_adjacent > 1  # internal edge
#       println("this is an internal edge")
      element1 = adjacent_nums[1]
      element2 = adjacent_nums[2]

      coords_1 = x[:, :, element1]
      coords_2 = x[:, :, element2]

      # calculate centroid
      centroid1 = sum(coords_1, 2)
      centroid2 = sum(coords_2, 2)

      if abs(centroid1[1] - centroid2[2]) < 1e-10  # if big enough difference
        if centroid1[1] < centroid2[1]
    elementL = element1
    elementR = element2
        else
    elementL = element2
    elementR = element1
        end
      else  # use y coordinate to decide
        if centroid1[2] < centroid2[2]
    elementL = element1
    elementR = element2
        else
    elementL = element2
    elementR = element1
        end
      end

      edgeL = getEdgeLocalNum(mesh, i, elementL)
      edgeR = getEdgeLocalNum(mesh, i, elementR)

      return elementL, elementR, edgeL, edgeR
    end

=#



function saveSolutionToMesh(mesh::PumiMesh2, u::AbstractVector)
# saves the solution in the vector u to the mesh (in preparation for mesh adaptation
# linear meshes only

  dofnums = zeros(Int, mesh.numDofPerNode)
  u_vals = zeros(mesh.numDofPerNode)
  for i=1:mesh.numVert
    for j = 1:mesh.numDofPerNode  # get dof numbers of this vertex
      dofnums[j] = getNumberJ(mesh.dofnums_Nptr, mesh.verts[i], 0, j-1)
    end
#    println("dofnums = ", dofnums)
    u_vals[:] = u[dofnums]
    setComponents(mesh.f_ptr, mesh.verts[i], 0, u_vals)
  end

  return nothing
end  # end function saveSolutionToMesh


function retrieveSolutionFromMesh(mesh::PumiMesh2, u::AbstractVector)
# retrieve solution from mesh (after mesh adaptation)
# mesh needs to have been reinitilized after mesh adaptation, u needs to be the right size for the new mesh
  # linear meshes only

  dofnums = zeros(Int, mesh.numDofPerNode)
  u_vals = zeros(mesh.numDofPerNode)
  for i=1:mesh.numVert
    for j = 1:mesh.numDofPerNode  # get dof numbers of this vertex
      dofnums[j] = getNumberJ(mesh.dofnums_Nptr, mesh.verts[i], 0, j-1)
    end
    getComponents(mesh.f_ptr, mesh.verts[i], 0, u_vals)
    u[dofnums] = u_vals
  end

  return nothing
end  # end function saveSolutionToMesh


function retrieveNodeSolution(f_ptr, entity, u_node::AbstractVector)
# retrieve solution on a particular entity, stores it in u_node
# u_node must be a vector of length mesh.numDofPerNode
# used during mesh adaptation
# this is a low level function because t takes in f_ptr, entity rather than the mesh object and an index

  getComponents(f_ptr, entity, 0, u_node)

end

end  # end of module
