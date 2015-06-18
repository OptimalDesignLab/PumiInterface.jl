export PumiMesh3

type PumiMesh3{T1} <: AbstractMesh{T1}   # 3d pumi mesh, tetrahedron only
  m_ptr::Ptr{Void}  # pointer to mesh
  mshape_ptr::Ptr{Void} # pointer to mesh's FieldShape
  f_ptr::Ptr{Void} # pointer to apf::field for storing solution during mesh adaptation
  vert_Nptr::Ptr{Void}  # numbering of vertices
  edge_Nptr::Ptr{Void}  # numbering of edges
  face_Nptr::Ptr{Void} # numbering of faces
  el_Nptr::Ptr{Void}  # numbering of elements (faces)

  numVert::Int  # number of vertices in the mesh
  numEdge::Int # number of edges in the mesh
  numFace::Int # number of faces in the mesh
  numEl::Int  # number of elements (regions)
  order::Int # order of shape functions
  numDof::Int # number of degrees of freedom
  numNodes::Int  # number of nodes
  numDofPerNode::Int  # number of dofs per node
  numBoundaryFaces::Int # number of faces on the exterior boundary
  numInterfaces::Int # number of internal interfaces
  numNodesPerElement::Int  # number of nodes per element
  numNodesPerType::Array{Int, 1}  # number of nodes classified on each vertex, edge, face

  # hold pointers to mesh entities
  verts::Array{Ptr{Void},1}  # holds pointers to mesh vertices
  edges::Array{Ptr{Void},1}  # pointers to mesh edges
  faces::Array{Ptr{Void},1}  # pointers to mesh faces
  elements::Array{Ptr{Void},1}  # pointers to regions

  dofnums_Nptr::Ptr{Void}  # pointer to Numbering of dofs (result of reordering)
  boundary_nums::Array{Int, 2}  # array of [element number, edgenumber] for each edge on the boundary

  bndryfaces::Array{Boundary, 1}  # store data on external boundary of mesh
  interfaces::Array{Interface, 1}  # store data on internal edges

  coords::Array{T1, 3}  # store coordinates of all nodes
  dxidx::Array{T1, 4}  # store scaled mapping jacobian
  jac::Array{T1,2}  # store mapping jacobian output

  dofs::Array{Int, 3}  # store dof numbers of solution array to speed assembly



 function PumiMesh3(dmg_name::AbstractString, smb_name::AbstractString, order, sbp::SBPOperator; dofpernode=1)
  # construct pumi mesh by loading the files named
  # dmg_name = name of .dmg (geometry) file to load (use .null to load no file)
  # smb_name = name of .smb (mesh) file to load
  # order = order of shape functions
  # dofpernode = number of dof per node, default = 1

  mesh = new()
  mesh.numDofPerNode = dofpernode
  mesh.order = order
  mesh.numNodesPerElement = getNumNodes(order)
  num_Entities, mesh.m_ptr, mesh.mshape_ptr = init(dmg_name, smb_name, order)
  mesh.f_ptr = createPackedField(mesh.m_ptr, "solution_field", dofpernode)

  mesh.numVert = convert(Int, num_Entities[1])
  mesh.numEdge = convert(Int,  num_Entities[2])
  mesh.numFace = convert(Int, num_Entities[3])
  mesh.numEl = convert(Int, num_Entities[4])

  # get pointers to mesh entity numberings
  mesh.vert_Nptr = getVertNumbering()
  mesh.edge_Nptr = getEdgeNumbering()
  mesh.face_Nptr = getFaceNumbering()
  mesh.el_Nptr = getElNumbering()

  mesh.verts = Array(Ptr{Void}, mesh.numVert)
  mesh.edges = Array(Ptr{Void}, mesh.numEdge)
  mesh.faces = Array(Ptr{Void}, mesh.numFace)
  mesh.elements = Array(Ptr{Void}, mesh.numEl)
  mesh.dofnums_Nptr = createNumberingJ(mesh.m_ptr, "reordered dof numbers", mesh.mshape_ptr, dofpernode)  # 1 dof per node

  # get pointers to all MeshEntities
  # also initilize the field to zero
  resetAllIts2()
#  comps = zeros(dofpernode)
  for i=1:mesh.numVert
    mesh.verts[i] = getVert()
    incrementVertIt()
  end

  for i=1:mesh.numEdge
    mesh.edges[i] = getEdge()
    incrementEdgeIt()
  end

  println("getting Face pointers")
  for i=1:mesh.numFace
    mesh.faces[i] = getFace()
    incrementFaceIt()
  end

  println("getting region pointers")
  for i=1:mesh.numEl
    mesh.elements[i] = getEl()
    incrementElIt()
  end

  # use partially constructed mesh object to populate arrays

  
#  numberDofs(mesh)

  countBoundaryFaces(mesh)
#=
  mesh.bndryfaces = Array(Boundary, mesh.numBoundaryEdges)
  getBoundaryArray(mesh)

  mesh.numInterfaces = mesh.numEdge - mesh.numBoundaryEdges
  mesh.interfaces = Array(Interface, mesh.numInterfaces)
  getInterfaceArray(mesh)

  getCoordinates(mesh, sbp)  # store coordinates of all nodes into array
  getDofNumbers(mesh)  # store dof numbers

  mesh.dxidx = Array(T1, 2, 2, sbp.numnodes, mesh.numEl)
  mesh.jac = Array(T1, sbp.numnodes, mesh.numEl)
  mappingjacobian!(sbp, mesh.coords, mesh.dxidx, mesh.jac)
=#
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
  println("numVert = ", mesh.numVert)
  println("numEdge = ", mesh.numEdge)
  println("numEl = ", mesh.numEl)
  println("numDof = ", mesh.numDof)
  println("numNodes = ", mesh.numNodes)



  writeVtkFiles("mesh_complete", mesh.m_ptr)
  return mesh
  # could use incomplete initilization to avoid copying arrays
#  return PumiMesh2(m_ptr, mshape_ptr, f_ptr, vert_Nptr, edge_Nptr, el_Nptr, numVert, numEdge, numEl, order, numdof, numnodes, dofpernode, bnd_edges_cnt, verts, edges, elements, dofnums_Nptr, bnd_edges_small)
  end

 
  
end

function getNumNodes(order::Integer)
# get the number of nodes on an element
# assumes SBP elements

  nnodes = 0
  if order == 1
    nnodes = 4
  elseif order == 2
    nnodes = 11
  elseif order == 3
    nnodes = 24
  elseif order == 4
    nnodes = 45
  else
    println("Warning, unsupported order elements are requested")
  end

  return nnodes
end




function countBoundaryFaces(mesh::PumiMesh3)
  # count boundary face
  # store array of [element number, global edge number]
  resetFaceIt()
  bnd_faces_cnt = 0
  bnd_faces = Array(Int, mesh.numFace, 2)
  for i=1:mesh.numFace
    face_i = getFace()
    numRegion = countAdjacent(mesh.m_ptr, face_i, 3)  # should be count upward

    if numRegion == 1  # if an exterior edge
      elements = getAdjacent(numRegion)
      elnum = getElNumber2(elements[1]) + 1

      println("face number ", i, " is part of only one element")

      bnd_faces_cnt += 1
      bnd_faces[bnd_faces_cnt, 1] = elnum
      bnd_faces[bnd_faces_cnt, 2] = i
    end
    incrementFaceIt()
  end


  mesh.boundary_nums = bnd_faces[1:bnd_faces_cnt, :] # copy, bad but unavoidable
  mesh.numBoundaryFaces = bnd_faces_cnt

return nothing

end  # end function



