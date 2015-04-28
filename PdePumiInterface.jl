module PdePumiInterface
push!(LOAD_PATH, "/users/creanj/julialib_fork/PUMI.jl")

using PumiInterface
using SummationByParts

export AbstractMesh,PumiMesh2, getElementVertCoords, getShapeFunctionOrder, getGlobalNodeNumber, getGlobalNodeNumbers, getNumEl, getNumEdges, getNumVerts, getNumNodes, getNumDofPerNode, getBoundaryEdgeNums, getBoundaryFaceNums, getBoundaryEdgeLocalNum, getBoundaryArray

abstract AbstractMesh

type PumiMesh2 <: AbstractMesh   # 2d pumi mesh, triangle only
  m_ptr::Ptr{Void}  # pointer to mesh
  mshape_ptr::Ptr{Void} # pointer to mesh's FieldShape

  numVert::Int  # number of vertices in the mesh
  numEdge::Int # number of edges in the mesh
  numEl::Int  # number of elements (faces)
  order::Int # order of shape functions
  numDof::Int # number of degrees of freedom
  numNodes::Int  # number of nodes
  numDofPerNode::Int  # number of dofs per node
  numBoundaryEdges::Int # number of edges on the exterior boundary

  verts::Array{Ptr{Void},1}  # holds pointers to mesh entities
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

  numVert = convert(Int, num_Entities[1])
  numEdge =convert(Int,  num_Entities[2])
  numEl = convert(Int, num_Entities[3])

  verts = Array(Ptr{Void}, numVert)
  edges = Array(Ptr{Void}, numEdge)
  elements = Array(Ptr{Void}, numEl)
  dofnums_Nptr = createNumberingJ(m_ptr, "reordered dof numbers", mshape_ptr, dofpernode)  # 1 dof per node

  # get pointers to all MeshEntities
  resetAllIts2()

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
    faces = getAdjacent(numFace)

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

  writeVtkFiles("mesh_complete", m_ptr)
  return PumiMesh2(m_ptr, mshape_ptr, numVert, numEdge, numEl, order, numdof, numnodes, dofpernode, bnd_edges_cnt, verts, edges, elements, dofnums_Nptr, bnd_edges_small)
end








function getElementVertCoords(mesh::PumiMesh2,  elnums::Array{Int,1})
# elnums = vector of element numbbers
# return array of size 3x3, where each column contains the coordinates of a vertex


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

function getBoundaryEdgeNums(mesh::PumiMesh2)
# get vector of edge numbers that are on boundary

edge_nums = mesh.boundary_edge_nums
return edge_nums
end



function getBoundaryFaceNums(mesh::PumiMesh2)
# get array of face numbers that are on boundary
# this only works in 2d
# this function doesn't work yet
#face_nums = zeros(Int, numEdges_on_boundary)
face_nums = zeros(Int, 0)
return face_nums
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

function getBoundaryArray(mesh::PumiMesh2)
# get an array of type Boundary for SBP
# creating an an array of a user defined type seems like a waste of memory operations

  bnd_array = Array(Boundary, mesh.numBoundaryEdges)

  for i=1:mesh.numBoundaryEdges
    facenum = mesh.boundary_nums[i,1]
    edgenum_global = mesh.boundary_nums[i,2]
    edgenum_local = getBoundaryEdgeLocalNum(mesh, edgenum_global)
    bnd_array[i] = Boundary(facenum, edgenum_local)
  end

  return bnd_array
end



 

end  # end of module
