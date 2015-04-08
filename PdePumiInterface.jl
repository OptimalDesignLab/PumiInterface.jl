module PdePumiInterface
push!(LOAD_PATH, "/users/creanj/julialib_fork/PUMI.jl")

using PumiInterface


export PumiMesh2, getElementVertCoords, getShapeFunctionOrder, getGlobalNodeNumber, getNumEl, getNumEdges, getNumVerts, getNumNodes, getBoundaryEdgeNums, getBoundaryFaceNums

abstract AbstractMesh

type PumiMesh2   # 2d pumi mesh, triangle only
  m_ptr::Ptr{Void}  # pointer to mesh
  mshape_ptr::Ptr{Void} # pointer to mesh's FieldShape

  numVert::Int  # number of vertices in the mesh
  numEdge::Int # number of edges in the mesh
  numEl::Int  # number of elements (faces)
  order::Int # order of shape functions
  numDof::Int # number of degrees of freedom
  numBoundaryEdges::Int # number of edges on the exterior boundary

  verts::Array{Ptr{Void},1}  # holds pointers to mesh entities
  edges::Array{Ptr{Void},1}  # pointers to mesh edges
  elements::Array{Ptr{Void},1}  # pointers to faces

  dofnums_Nptr::Ptr{Void}  # pointer to Numbering of dofs (result of reordering)
  boundary_edge_nums::Array{Int, 1}
  
  
end

function PumiMesh2(dmg_name::AbstractString, smb_name::AbstractString, order)
  # construct pumi mesh by loading the files named


  tmp, num_Entities, m_ptr, mshape_ptr = init2(dmg_name, smb_name, order)

  numVert = convert(Int, num_Entities[1])
  numEdge =convert(Int,  num_Entities[2])
  numEl = convert(Int, num_Entities[3])

  verts = Array(Ptr{Void}, numVert)
  edges = Array(Ptr{Void}, numEdge)
  elements = Array(Ptr{Void}, numEl)
  dofnums_Nptr = createNumberingJ(m_ptr, "reordered dof numbers", mshape_ptr, 1)  # 1 dof per node

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
  numdof = order*numVert  # works for first and second order
  # number dofs
  ctr= 1
  for i=1:numVert
      numberJ(dofnums_Nptr, verts[i], 0, 0, ctr)
      println("vertex ", i,  " numbered ", ctr)
      ctr += 1
    end

  if order >= 2
    for i=1:numEdges
      numberJ(dofnums_Nptr, edges[i], 0, 0, ctr)
      ctr += 1
    end
  end

 

  # count boundary edges
  bnd_edges_cnt = 0
  bnd_edges = Array(Int, numEdge)
  for i=1:numEdge
    edge_i = getEdge()
    numFace = countAdjacent(m_ptr, edge_i, 2)  # should be count upward
    if numFace == 1  # if an exterior edge
      bnd_edges_cnt += 1
      bnd_edges[bnd_edges_cnt] = i
    end
    incrementEdgeIt()
  end

  bnd_edges_small = bnd_edges[1:bnd_edges_cnt]

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

  return PumiMesh2(m_ptr, mshape_ptr, numVert, numEdge, numEl, order, numdof, bnd_edges_cnt, verts, edges, elements, dofnums_Nptr, bnd_edges_small)
end








function getElementVertCoords(mesh::PumiMesh2,  elnums::Array{Int,1})
# elnums = array of element numbbers
# return array 


# get the number of verts 
n = length(elnums)

coords = zeros(3, 3, n)  # first dimension = x,y, or z, second dimension = vertex number

for i=1:n
  println("i = ", i)
  elnum_i = elnums[i]
  el_i = mesh.elements[elnum_i]
  sub_array = sub(coords, :, :, i)
  getFaceCoords(el_i, sub_array, 3, 3)  # populate coords
end

return coords

end

function getShapeFunctionOrder(mesh::PumiMesh2)

return mesh.order
end


function getGlobalNodeNumber(mesh::PumiMesh2, el_num::Integer, local_node_num::Integer)
# el_num = element number
# local_node_num = node number within element
# this only works for first and second order


println("el_num = ", el_num, " local_node_num = ", local_node_num)
el = mesh.elements[el_num]  # get element
node_entities = getNodeEntities(mesh.m_ptr, mesh.mshape_ptr, el)
println("node_entities = ", node_entities)
node_entity = node_entities[local_node_num]
println("node_entity = ", node_entity)
number_i = getNumberJ(mesh.dofnums_Nptr, node_entities[local_node_num], 0, 0)
println("global dof number = ", number_i)

return number_i

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

return mesh.numdof
end

function getBoundaryEdgeNums(mesh::PumiMesh2)
# get array of edge numbers that are on boundary

edge_nums = mesh.boundary_edge_nums
return edge_nums
end

function getBoundaryFaceNums(mesh::PumiMesh2)
# get array of face numbers that are on boundary
# this only works in 2d
#face_nums = zeros(Int, numEdges_on_boundary)
face_nums = zeros(Int, 0)
return face_nums
end

end
