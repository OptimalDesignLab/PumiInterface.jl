# PdePumiInterfaceDG.jl: implements AbstractMesh for a 2D Discontinuous
#                        Galerkin mesh

export AbstractMesh,PumiMeshDG2, PumiMeshDG2Preconditioning, reinitPumiMeshDG2, getElementVertCoords, getShapeFunctionOrder, getGlobalNodeNumber, getGlobalNodeNumbers, getNumEl, getNumEdges, getNumVerts, getNumNodes, getNumDofPerNode, getAdjacentEntityNums, getBoundaryEdgeNums, getBoundaryFaceNums, getBoundaryEdgeLocalNum, getEdgeLocalNum, getBoundaryArray, saveSolutionToMesh, retrieveSolutionFromMesh, retrieveNodeSolution, getAdjacentEntityNums, getNumBoundaryElements, getInterfaceArray, printBoundaryEdgeNums, printdxidx, getdiffelementarea, writeVisFiles

# Element = an entire element (verts + edges + interior face)
# Type = a vertex or edge or interior face
# 

# the vert, edge, face Nptrs are numberings that map a MeshEntity to its 
# index the the arrays verts, edges, faces, except that they are zero based indices

# internally, all operations and loops are done by looping over nodes 
# in the Pumi order (ie. vertex nodes counter clockwise, edge nodes
# counterclockwise, the face nodes counterclockwise starting in with the
# node nearest the first vertex).
# The external interface remaps the nodes into an order specified
# by the nodemaps nodemapSbptoPumi and nodemapPumitoSbp.
# The external interface is compreised of the arrays:
# dofs, coords, dxidx, jac, sparsity_bnds 
# Care must be taken with the visualization files to map back to the Pumi
# order

# For DG meshes, two different FieldShapes are used: the coordinate FieldShape,
# which coveres vertices (and might cover mid edge nodes in the future for 
# curved meshes), and the solution Fieldshape, which describes the nodes
# Pumi *only* stores the coordinates of nodes in the coordinate FieldShape,
# so the interface between what Pumi knows about elements and what Julia
# knows about element is the follow:
#  Pumi keeps track of vertices (and possibly mid edge nodes)
#  SBP calculates the location of all nodes from the vertex coordinates
#  SBP calculates the mapping jacobian from the node coordinates
#
# In order to store data on the mesh in the solution FieldShape, Pumi
# must know at least the number and classification of the nodes.
# In practice, we define a Pumi node ordering as follows:
#  vertex nodes, from the bottom left, going counter-clockwise
#  edge nodes: edge 1, edge 2, edge 3
#  face nodes: starting at the node closest to vertex 1, proceeding 
#              counter-clockwise, outermost to innermost
# This requirement is stronger than what is actually necessary, but 
# having an explicitly defined node ordering convention and a general
# mechanism to convert from Pumi to SBP ordering is a better interface
#         
export PumiMeshDG2, PumiMeshDG
#abstract AbstractMesh
@doc """
### PumiInterface.PumiMeshDG2

  This is an implementation of AbstractMesh for a 2 dimensional equation.  
  The constructor for this type extracts all the needed information from Pumi,
  so the solver never needs access to Pumi.

  Fields:
    m_ptr: a Ptr{Void} to the Pumi mesh
    mnew_ptr: a Ptr{Void} to the subtriangulated mesh used for high
              order visualization, null pointer if not needed
    mshape_ptr: a Ptr{Void} do the FieldShape of the solution field (not the
                coordinate field)
    coordshape_ptr: a pointer to the FieldShape of the coordinate field

    f_ptr:  a Ptr{Void} to the Pumi apf::Field holding the solution
    fnew_ptr: a pointer to the solution Field on the subtriangulated mesh,
              null if not needed
    shape_type: an enum describing what shape functions to use
    min_node_dist: the minimum distance between nodes on a reference element
    vert_NPtr: a pointer to the apf::Numbering object that numbers the 
               vertices (0-based numbering)
    edge_Nptr: a pointer to the Pumi apf::Numbering that numbers the mesh 
               edges (0-based numbering)
    el_Nptr:  a pointer to the Pumi apf::Numbering that numbers the mesh
              elements (0-based numbering)
    coloring_Nptr: a pointer to the Pumi Numbering containing the element
                   coloring

    numVert: number of vertices in the mesh
    numEdge: number of edges in the mesh
    numEl: number of elements in the mesh
    order: degree of the shape functions
    numDof: total number of degrees of freedom
    numNodes: number of nodes in the mesh
    numDofPerNode: number of degrees of freedom on each node
    numBoundaryFaces: number of edges on the boundary of the domain
    numInterfaces: number of internal edges (edges not on boundary)
    numNodesPerElements: number of nodes on an element
    numNodesPerType: array of length 3 that tells how many nodes are on a mesh
                     vertex, edge, and element
    numNodesPerFace: number of nodes on each face of the element (edges in 2D)
    numEntitiesPertype: [numVert, numEdges, numEl]
    numTypesPerElement: array of length 3 telling how many vertices, edges
                         and faces in each element 
    
    typeOffsetsPerElement: array of length 3 telling the starting index of the 
                           vertices, edges, and faces in an array where
                           all quantities for the vertices are stored, then
                           the edges, then the faces

    typeOffsetsPerElement_: Int32 version of the above
    nodemapSBPtoPumi: array of UInt8s that maps the SBP ordering of nodes to 
                      the Pumi ordering of nodes 
    nodemapPumiToSBP: array of UInt8s that maps the Pumi node ordering to the
                      SBP one.

    coloringDistance: the distance k of the distance-k graph coloring used to
                      color the elements (graph vertices are elements and graph
                      edges exist where elements share an edge)

    dofs: a 3D array holding the degree of freedom numbers for the elements,
          both locally owned and the ghost elements.  The first mesh.numEl
          elements are the locally owned ones, the rest are for the ghost 
          element.  Their offsets are specified by mesh.shared_element_offsets.
          The global dof number is the number stored in this array + 
          dof_offset  (even for the non-local elements)
"""->
type PumiMeshDG2{T1} <: PumiMesh2DG{T1}   # 2d pumi mesh, triangle only
  m_ptr::Ptr{Void}  # pointer to mesh
  mnew_ptr::Ptr{Void}  # pointer to subtriangulated mesh (high order only)
  mshape_ptr::Ptr{Void} # pointer to the FieldShape of the node field
  coordshape_ptr::Ptr{Void}  # pointer to FieldShape of the coordinate 
                             # field
  f_ptr::Ptr{Void} # pointer to apf::field for storing solution during mesh adaptation
  fnew_ptr::Ptr{Void}  # pointer to field on mnew_ptr
  shape_type::Int  #  type of shape functions
  min_node_dist::Float64  # minimum distance between nodes
  min_el_size::Float64 # size of the smallest element (units of length)
  f::IOStream
  vert_Nptr::Ptr{Void}  # numbering of vertices (zero based)
  edge_Nptr::Ptr{Void}  # numbering of edges (zero based)
  el_Nptr::Ptr{Void}  # numbering of elements (faces)  (zero based)
  coloring_Nptr::Ptr{Void}  # coloring of mesh for sparse jacobian calculation
  entity_Nptrs::Array{Ptr{Void}, 1}  # [vert_Nptr, edge_Nptr, el_Nptr], 0-based
  numVert::Int  # number of vertices in the mesh
  numEdge::Int # number of edges in the mesh
  numFace::Int # number of faces (=edge) in the mesh
  numEl::Int  # number of elements (faces)
  numGlobalEl::Int  # number of local elements + number of non-local elements
  numSharedEl::Int  # number of non-local elements
  order::Int # order of shape functions
  numDof::Int # number of degrees of freedom
  numNodes::Int  # number of nodes
  numDofPerNode::Int  # number of dofs per node
  numBoundaryFaces::Int # number of edges on the exterior boundary
  numInterfaces::Int # number of internal interfaces
  numNodesPerElement::Int  # number of nodes per element
  numNodesPerType::Array{Int, 1}  # number of nodes classified on each vertex, edge, face
  numNodesPerFace::Int  # number of face nodes
  numEntitiesPerType::Array{Int, 1} # [numVert, numEdge, numEl]
  numTypePerElement::Array{Int, 1}  # number of verts, edges, faces per element
  typeOffsetsPerElement::Array{Int, 1} # the starting index of the vert, edge, and face nodes in an element 
  typeOffsetsPerElement_::Array{Int32, 1}  # Int32 version of above
  nodemapSbpToPumi::Array{UInt8, 1}  # maps nodes of SBP to Pumi order
  nodemapPumiToSbp::Array{UInt8, 1}  # maps nodes of Pumi to SBP order

  # constants needed by Pumi
  el_type::Int  # apf::Type for the elements of the mesh
  face_type::Int # apf::Type for the faces of the mesh

  # parallel bookkeeping info
  comm::MPI.Comm  # MPI Communicator
  myrank::Int  # MPI rank, zero based
  commsize::Int # MPI comm size
  peer_parts::Array{Int, 1}  # array of part numbers that share entities with
                             # the current part
  npeers::Int  # length of the above array
  peer_face_counts::Array{Int, 1}  # number of edges on each part
  send_waited::Array{Bool, 1}  # indicates whether the send_reqs have been
                               # waited on already
  send_reqs::Array{MPI.Request, 1}  # Request objects for sends
  send_stats::Array{MPI.Status, 1} # Status objects for sends
  recv_waited::Array{Bool, 1}  # similar to send_waited
  recv_reqs::Array{MPI.Request, 1}  # Request objects for sends
  recv_stats::Array{MPI.Status, 1}  # status objects for sends

  ref_verts::Array{Float64, 2}  # 2 x 3 array giving the coordinates
                                # of the vertices of the reference
                                # element in parametric space
  dim::Int  # dimension of mesh (2 or 3D)
  isDG::Bool  # is this a DG mesh (always true)
  isInterpolated::Bool # does the solution need to be interpolated to
                       # the boundaries
  coloringDistance::Int  # distance between elements of the same color, measured in number of edges
  numColors::Int  # number of colors
  maxColors::Int  # maximum number of colors on any processor
  numBC::Int  # number of boundary conditions

  # hold pointers to mesh entities
  verts::Array{Ptr{Void},1}  # holds pointers to mesh vertices
  edges::Array{Ptr{Void},1}  # pointers to mesh edges
  faces::Array{Ptr{Void},1}  # alias for edges
  elements::Array{Ptr{Void},1}  # pointers to faces

  # used for high order elements to determine the orientations of edges and 
  # faces (and vertices, too, for consistency)
  # stored in Pumi node ordering, not SBP
  elementNodeOffsets::Array{UInt8, 2}
  # truth values if the entity is oriented consistently with the element
  # Array of length 3, whose elements are numTypePerelement x numEl
  typeNodeFlags::Array{BitArray{2}, 1}
  triangulation::Array{Int32, 2}  # sub triangulation of an element

  nodestatus_Nptr::Ptr{Void}  # node status pointer, status >= 2 -> give it a
                              # node number
  nodenums_Nptr::Ptr{Void}    # node numbering itself
  dofnums_Nptr::Ptr{Void}  # pointer to Numbering of dofs (result of reordering)
#  boundary_nums::Array{Int, 2}  # array of [element number, edgenumber] for each edge on the boundary

  bndry_funcs::Array{BCType, 1}  # array of boundary functors (Abstract type)
  bndry_normals::Array{T1, 3}  # array of normals to each face on the boundary
#  bndry_facenums::Array{Array{Int, 1}, 1}  # hold array of faces corresponding to each boundary condition
  bndry_offsets::Array{Int, 1}  # location in bndryfaces where new type of BC starts
                                # and one past the end of the last BC type
				# array has length numBC + 1

  bndryfaces::Array{Boundary, 1}  # store data on external boundary of mesh
  interfaces::Array{Interface, 1}  # store data on internal edges
  interface_normals::Array{T1, 4}  # array of interior face normals, 
                                   # indexing: [norm_component, Left or right, left face node, left face element]

  coords::Array{T1, 3}  # store coordinates of all nodes
  coords_bndry::Array{T1, 3}  # store coordinates of nodes on boundary,
                              # 2 x numFaceNodes x numBoundaryFaces
  coords_sharedface::Array{Array{T1, 3}, 1}  # coordinates of shared interface nodes
  dxidx::Array{T1, 4}  # store scaled mapping jacobian
  dxidx_face::Array{T1, 4} # store scaled mapping jacobian at face nodes
                           # 2 x 2 x numfacenodes x numInterfaces
  dxidx_sharedface::Array{Array{T1, 4}, 1}  # array of arrays for dxidx
                                            # on shared edges
  dxidx_bndry::Array{T1, 4} # store scaled mapping jacobian at boundary nodes,
                            # similar to dxidx_face
  jac::Array{T1,2}  # store mapping jacobian output
  jac_face::Array{T1,2}  # store jacobian determanent at face nodes
                         # numfacenodes x numInterfaces
  jac_sharedface::Array{Array{T1, 2}, 1}  # array of arrays for shared
                                          # edge jacobian determinent

  jac_bndry::Array{T1, 2} # store jacobian determinant at boundry nodes
                          # similar to jac_bndry

  dof_offset::Int  # local to global offset for dofs
  dofs::Array{Int, 3}  # store dof numbers of solution array to speed assembly
  sparsity_bnds::Array{Int32, 2}  # store max, min dofs for each dof
  sparsity_nodebnds::Array{Int32, 2}  # store min, max nodes for each node

  sparsity_counts::Array{Int32, 2}  # store number of local, remote dofs 
                                    # for each dof
  sparsity_counts_node::Array{Int32, 2}  # store number of local, remote dofs
                                         # for each node
  color_masks::Array{BitArray{1}, 1}  # array of bitarray masks used to control element perturbations when forming jacobian, number of arrays = number of colors
  neighbor_colors::Array{UInt8, 2}  # 4 by numEl array, holds colors of edge-neighbor elements + own color
  neighbor_nums::Array{Int32, 2}  # 4 by numEl array, holds element numbers of neighbors + own number, in same order as neighbor_colors
  pertNeighborEls::Array{Int32, 2}  # numEl by numcolors array, for each color,
                                    # stores the element number of the element
				    # whose perturbation is affecting the 
				    # current element
  pertNeighborEls_edge::Array{Int32, 2}  # numEl by 3 (number of edges per element)
                                         # stores the element number of the 
					 # neighboring element that shares 
					 #  edge, where edge is the local edge
					 # index (ie. edge 1,2, or 3)
  color_cnt::Array{Int32, 1}  # number of elements in each color

  interp_op::Array{Float64, 2}  # 3 x numNodesPerEl matrix that interpolates
                                # solution values to the vertices

  ###### parallel data #####

  # array of arrays of Boundary objects, one array for each peer part
  # the Boundary describes the local part of a shared interface
  bndries_local::Array{Array{Boundary, 1}, 1} 
  # array of arrays of Boundary objects from the *remote* process describing its
  # part of a shared interface, on array for each peer part
  # ie. the element number is the element number of the remote part.  Use with
  # caution.
  bndries_remote::Array{Array{Boundary, 1}, 1} # hold a Boundary object from

  # array of arrays of Interface objects giving a complete description of the
  # interfaces that span part boundaries.  elementL is always the local element
  # and elementR is always the element that lives on the remote process.
  # elementR has been given an element number > numEl
  shared_interfaces::Array{Array{Interface, 1}, 1}

  # array of length npeers + 1, holding the first element number for 
  # each peer part, with the last element being 1 more than the highest element
  # number for the last peer
  shared_element_offsets::Array{Int, 1}
  local_element_counts::Array{Int, 1}  # number of element on the local
                                       # side of each peer boundary
  remote_element_counts::Array{Int, 1}

  # element numbers of the elements on each peer boundary
  local_element_lists::Array{Array{Int32, 1}, 1}

  # color masks for the non-local elements, contains npeers x numColor
  # BitArrays, where each BitArray has length = the number of ghost elements
  # on the current peer boundary
  shared_element_colormasks::Array{Array{BitArray{1}, 1}, 1}                               
  sbpface::TriFace{Float64}  # SBP object needed to do interpolation


 function PumiMeshDG2(dmg_name::AbstractString, smb_name::AbstractString, order, sbp::AbstractSBP, opts, interp_op, sbpface; dofpernode=1, shape_type=2, coloring_distance=2, comm=MPI.COMM_WORLD)
  # construct pumi mesh by loading the files named
  # dmg_name = name of .dmg (geometry) file to load (use .null to load no file)
  # smb_name = name of .smb (mesh) file to load
  # order = order of shape functions
  # opts: dictionary of options
  # dofpernode = number of dof per node, default = 1
  # shape_type = type of shape functions, 0 = lagrange, 1 = SBP, 2 = SBP DG1
  # coloring_distance : distance between elements of the same color, where distance is the minimum number of edges that connect the elements, default = 2

  println("\nConstructing PumiMeshDG2 Object")
  println("  sbp_name = ", smb_name)
  println("  dmg_name = ", dmg_name)
  mesh = new()
  mesh.isDG = true
  mesh.dim = 2
  mesh.numDofPerNode = dofpernode
  mesh.order = order
  mesh.shape_type = shape_type
  mesh.coloringDistance = coloring_distance
  mesh.interp_op = interp_op
  mesh.ref_verts = [0.0 1 0; 0 0 1]  # ???
  mesh.numNodesPerFace = sbpface.numnodes
  mesh.comm = comm

  if !MPI.Initialized()
    MPI.Init()
  end

  mesh.myrank = MPI.Comm_rank(mesh.comm)
  mesh.commsize = MPI.Comm_size(mesh.comm)
  myrank = mesh.myrank
  mesh.f = open("meshlog_$myrank.dat", "w")

  if sbp.numfacenodes == 0
    mesh.sbpface = sbpface
    mesh.isInterpolated = true
  else
    mesh.isInterpolated = false
    # leave mesh.sbpface undefined - bad practice
  end

  # figure out coordinate FieldShape, node FieldShape
  coord_shape_type = 0 # integer to indicate the FieldShape of the coordinates
  field_shape_type = 0 # integer to indicate the FieldShape of the nodes
  mesh_order = order  # order of the coordinate field
  if shape_type == 2
    coord_shape_type = 0  # lagrange
    field_shape_type = shape_type
    mesh_order = 1
  else  # same coordinate, field shape
    coord_shape_type = shape_type
    field_shape_type = shape_type
    mesh_order = order
  end

  num_Entities, mesh.m_ptr, mesh.coordshape_ptr = init2(dmg_name, smb_name, mesh_order, shape_type=coord_shape_type)

  # create the solution field
  mesh.mshape_ptr = getSBPShapes(field_shape_type, order)
  mesh.f_ptr = createPackedField(mesh.m_ptr, "solution_field", dofpernode, mesh.mshape_ptr)
  mesh.min_node_dist = minNodeDist(sbp, mesh.isDG)

  # count the number of all the different mesh attributes
  mesh.numVert = convert(Int, num_Entities[1])
  mesh.numEdge =convert(Int,  num_Entities[2])
  mesh.numFace = mesh.numEdge
  mesh.numEl = convert(Int, num_Entities[3])
  mesh.numEntitiesPerType = [mesh.numVert, mesh.numEdge, mesh.numEl]
  mesh.numTypePerElement = [3, 3, 1]
  mesh.el_type = apfTRIANGLE
  mesh.face_type = apfEDGE

  num_nodes_v = countNodesOn(mesh.mshape_ptr, 0)  # number of nodes on a vertex
  num_nodes_e = countNodesOn(mesh.mshape_ptr, 1) # on edge
  num_nodes_f = countNodesOn(mesh.mshape_ptr, 2) # on face
  num_nodes_entity = [num_nodes_v, num_nodes_e, num_nodes_f]
  # count numbers of different things per other thing
  # use for bookkeeping
  mesh.numNodesPerType = num_nodes_entity
  mesh.typeOffsetsPerElement = zeros(Int, 4)
  pos = 1
  mesh.typeOffsetsPerElement[1] = pos
  for i=2:4
    pos += mesh.numTypePerElement[i-1]*mesh.numNodesPerType[i-1]
    mesh.typeOffsetsPerElement[i] = pos
  end
  mesh.typeOffsetsPerElement_ = [Int32(i) for i in mesh.typeOffsetsPerElement]

  mesh. numNodesPerElement = mesh.typeOffsetsPerElement[end] - 1
  numnodes = mesh.numNodesPerElement*mesh.numEl
  println("numNodesPerType = ", mesh.numNodesPerType)
  println("numEntitesPerType = ", mesh.numEntitiesPerType)
  mesh.numNodes = numnodes      # we assume there are no non-free nodes/dofs
  mesh.numDof = numnodes*dofpernode

  # get nodemaps
  mesh.nodemapSbpToPumi, mesh.nodemapPumiToSbp = getNodeMaps(order, shape_type, mesh.numNodesPerElement)

 
  # get pointers to mesh entity numberings
  mesh.vert_Nptr = getVertNumbering()
  mesh.edge_Nptr = getEdgeNumbering()
  mesh.el_Nptr = getFaceNumbering()
  mesh.entity_Nptrs = [mesh.vert_Nptr, mesh.edge_Nptr, mesh.el_Nptr]

  # create the coloring_Nptr
  el_mshape = getConstantShapePtr(2)
  mesh.coloring_Nptr = createNumberingJ(mesh.m_ptr, "coloring", el_mshape, 1)

  # create node status numbering (node, not dof)
  mesh.nodestatus_Nptr = createNumberingJ(mesh.m_ptr, "dof status", 
                         mesh.mshape_ptr, 1)   
  # create node numbering
  mesh.nodenums_Nptr = createNumberingJ(mesh.m_ptr, "reordered node numbers",
                       mesh.mshape_ptr, 1)

  # create dof numbering
  mesh.dofnums_Nptr = createNumberingJ(mesh.m_ptr, "reordered dof numbers", 
                      mesh.mshape_ptr, dofpernode)

  # populate node status numbering
  populateNodeStatus(mesh)


  # do node reordering
#=
 if opts["reordering_algorithm"] == "adjacency"
    start_coords = opts["reordering_start_coords"]
    # tell the algorithm there is only 1 dof per node because we only
    # want to label nodes
    reorder(mesh.m_ptr, mesh.numNodes, 1, 
            mesh.nodestatus_Nptr, mesh.nodenums_Nptr, mesh.el_Nptr, 
	    start_coords[1], start_coords[2])
=#
# elseif opts["reordering_algorithm"] == "default"
#    println("about to number nodes")
  numberNodes(mesh)

#  else
#    println(STDERR, "Error: invalid dof reordering algorithm requested")
#  end

#  println("finished numbering nodes")

#  println("about to number dofs")
  # do dof numbering
  populateDofNumbers(mesh)
#  println("finished numbering dofs")
 

  # get entity pointers
#  println("about to get entity pointers")
  mesh.verts, mesh.edges, mesh.elements = getEntityPointers(mesh)
#  println("finished getting entity pointers")

#  println("about to get boundary edge list")
  mesh.numBC = opts["numBC"]

  # create array of all model edges that have a boundary condition
  bndry_edges_all = Array(Int, 0)
  for i=1:mesh.numBC
    key_i = string("BC", i)
    bndry_edges_all = [ bndry_edges_all; opts[key_i]]  # ugly but easy
  end

#  println("finished getting boundary edge list")

#  println("about to count boundary edges")
 mesh.numBoundaryFaces, num_ext_edges, mesh.numInterfaces =  countBoundaryEdges(mesh, bndry_edges_all)
#  println("finished counting boundary edges")

  # populate mesh.bndry_faces from options dictionary
#  mesh.bndry_faces = Array(Array{Int, 1}, mesh.numBC)
#  println("about to get boudnary offets")
  mesh.bndry_offsets = Array(Int, mesh.numBC + 1)
  mesh.bndry_funcs = Array(BCType, mesh.numBC)
  boundary_nums = Array(Int, mesh.numBoundaryFaces, 2)

  offset = 1
  for i=1:mesh.numBC
    key_i = string("BC", i)
    println("opts[key_i] = ", opts[key_i])
#    println("typeof(opts[key_i]) = ", typeof(opts[key_i]))
    mesh.bndry_offsets[i] = offset
    offset = getMeshEdgesFromModel(mesh, opts[key_i], offset, boundary_nums)  # get the mesh edges on the model edge
    # offset is incremented by getMeshEdgesFromModel
  end


  mesh.bndry_offsets[mesh.numBC + 1] = offset # = num boundary edges
#  println("finished getting boundary offsets")

  # get array of all boundary mesh edges in the same order as in mesh.bndry_faces
#  boundary_nums = flattenArray(mesh.bndry_faces[i])
#  boundary_edge_faces = getEdgeFaces(mesh, mesh.bndry_faces)
  # use partially constructed mesh object to populate arrays

#  println("about to get entity orientations")
  mesh.elementNodeOffsets, mesh.typeNodeFlags = getEntityOrientations(mesh)
#  println("finished getting entity orientations")


  # start parallel initializiation
  colordata = getParallelInfo(mesh)

#  println("about to get degree of freedom numbers")
  getDofNumbers(mesh)  # store dof numbers
#  println("finished getting degree of freedom numbers")


  MPI.Barrier(mesh.comm)
  if coloring_distance == 2
    numc = colorMesh2(mesh, colordata)
    mesh.numColors = numc
    mesh.maxColors = MPI.Allreduce(numc, MPI.MAX, mesh.comm)
    println(mesh.f, "max colors = ", mesh.maxColors)
    mesh.color_masks = Array(BitArray{1}, numc)  # one array for every color
    mesh.neighbor_colors = zeros(UInt8, 4, mesh.numEl)
    mesh.neighbor_nums = zeros(Int32, 4, mesh.numEl)
    println("about to get distance 1 coloring")
    cnt, mesh.shared_element_colormasks = getColors1(mesh, colordata, mesh.color_masks, 
                                mesh.neighbor_colors, mesh.neighbor_nums; verify=opts["verify_coloring"] )
    mesh.pertNeighborEls = getPertNeighbors1(mesh)

  elseif coloring_distance == 0  # do a distance-0 coloring
    numc = colorMesh0(mesh)
    @assert numc == 1
    mesh.numColors = numc
    mesh.maxColors = numc
    mesh.color_masks = Array(BitArray{1}, numc)
    mesh.neighbor_colors = zeros(UInt8, 0, 0)  # unneeded array for distance-0
    mesh.neighbor_nums = zeros(Int32, 0, 0)  # unneeded for distance-0
    getColors0(mesh, mesh.color_masks)
    mesh.pertNeighborEls = getPertNeighbors0(mesh)

  else
    println(STDERR, "Error: unsupported coloring distance requested")
  end

  # get sparsity information
  # this takes into account the coloring distance
  #TODO: make getting sparsity bounds faster
  if opts["run_type"] != 1  # no a rk4 run
    println("getting sparsity bounds")
    mesh.sparsity_bnds = zeros(Int32, 0, 0)
#    mesh.sparsity_bnds = zeros(Int32, 2, mesh.numDof)
#    @time getSparsityBounds(mesh, mesh.sparsity_bnds)
#    mesh.sparsity_nodebnds = zeros(Int32, 2, mesh.numNodes)
#    @time getSparsityBounds(mesh, mesh.sparsity_nodebnds, getdofs=false)
#    println("finished getting sparsity bounds")

    mesh.sparsity_counts = zeros(Int32, 2, mesh.numDof)
    mesh.sparsity_counts_node = zeros(Int32, 2, mesh.numNodes)
    @time getSparsityCounts(mesh, mesh.sparsity_counts)
    @time getSparsityCounts(mesh, mesh.sparsity_counts_node, getdofs=false)
    println("finished getting sparsity counts")


  end


  # TODO: make this a separate option from use_edge_res, make decision
  #       in read_input.jl
  if opts["use_edge_res"]  # if using edge based residual
    mesh.pertNeighborEls_edge = getPertEdgeNeighbors(mesh)
  end

  # get boundary information for entire mesh
#  println("getting boundary info")
  mesh.bndryfaces = Array(Boundary, mesh.numBoundaryFaces)
  getBoundaryArray(mesh, boundary_nums)

  # need to count the number of internal interfaces - do this during boundary edge counting
#  println("getting interface info")
  println("num_ext_edges = ", num_ext_edges)
  println("numEdges = ", mesh.numEdge)
  println("numInterfaces = ", mesh.numInterfaces)
  mesh.interfaces = Array(Interface, mesh.numInterfaces)
  getInterfaceArray(mesh)
  sort!(mesh.interfaces)

  getCoordinates(mesh, sbp)  # store coordinates of all nodes into array

  mesh.dxidx = Array(T1, 2, 2, sbp.numnodes, mesh.numEl)
  mesh.jac = Array(T1, sbp.numnodes, mesh.numEl)
  mappingjacobian!(sbp, mesh.coords, mesh.dxidx, mesh.jac)

  mesh.min_el_size = getMinElementSize(mesh)

  # get face normals
  mesh.bndry_normals = Array(T1, 2, sbp.numfacenodes, mesh.numBoundaryFaces)
  getBoundaryFaceNormals(mesh, sbp, mesh.bndryfaces, mesh.bndry_normals)

  mesh.interface_normals = Array(T1, 2, 2, sbp.numfacenodes, mesh.numInterfaces)
  getInternalFaceNormals(mesh, sbp, mesh.interfaces, mesh.interface_normals)

  if mesh.isInterpolated
    mesh.dxidx_face, mesh.jac_face, mesh.dxidx_sharedface, mesh.jac_sharedface, mesh.dxidx_bndry, mesh.jac_bndry = interpolateMapping(mesh)

    mesh.coords_bndry = zeros(T1, 2, sbpface.numnodes, mesh.numBoundaryFaces)
    getBndryCoordinates(mesh, mesh.bndryfaces, mesh.coords_bndry)
    mesh.coords_sharedface = Array(Array{T1, 3}, mesh.npeers)
    for i=1:mesh.npeers
      mesh.coords_sharedface[i] = zeros(T1, 2, sbpface.numnodes, mesh.peer_face_counts[i])
      getBndryCoordinates(mesh, mesh.bndries_local[i], mesh.coords_sharedface[i])
    end
  end

  @time createSubtriangulatedMesh(mesh)
  println("finished creating sub mesh\n")

  println("printin main mesh statistics")

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


  # write data if requested
  myrank = mesh.myrank
  if opts["write_edge_vertnums"]
    rmfile("edge_vertnums_$myrank.dat")
    f = open("edge_vertnums_$myrank.dat", "a+")
    printEdgeVertNumbers(mesh.edge_Nptr, mesh.vert_Nptr, fstream=f)
    close(f)
  end

  if opts["write_face_vertnums"]
    rmfile("face_vertnums_$myrank.dat")
    f = open("face_vertnums_$myrank.dat", "a+")
    printFaceVertNumbers(mesh.el_Nptr, mesh.vert_Nptr, fstream=f)
    close(f)
  end

  if opts["write_boundarynums"]
    rmfile("boundary_nums_$myrank.dat")
    f = open("boundary_nums_$myrank.dat", "a+")
    println(f, boundary_nums)
    close(f)
  end

  if opts["write_dxidx"]
    rmfile("dxidx_$myrank.dat")
    printdxidx("dxidx_$myrank.dat", mesh.dxidx)
  end


  if opts["write_coords"]
    rmfile("coords_$myrank.dat")
    println("size(coords) = ", size(mesh.coords))
    writedlm("coords_$myrank.dat", mesh.coords)
#    printcoords("coords.dat", mesh.coords)
  end

  if opts["write_sparsity"]
    rmfile("sparsity_bnds_$myrank.dat")
    writedlm("sparsity_bnds_$myrank.dat", mesh.pertNeighborEls)
  end
#=
  if opts["write_sparsity_nodebnds"]
    println("writing sparsiy node bounds")
    rmfile("sparsity_nodebnds_$myrank.dat")
    writedlm("sparsity_nodebnds_$myrank.dat", mesh.sparsity_nodebnds)
  end
=#
  if opts["write_offsets"]
    rmfile("entity_offsets_$myrank.dat")
    writedlm("entity_offsets_$myrank.dat", mesh.elementNodeOffsets)
  end

  if opts["write_dofs"]
    rmfile("dofs_$myrank.dat")
    println("size(mesh.dofs) = ", size(mesh.dofs))
    writedlm("dofs_$myrank.dat", mesh.dofs)
  end

  if opts["write_counts"]
    writeCounts(mesh)
  end

  if opts["write_interfaces"]
    f = open("interface_$myrank.dat", "a+")
    ifaces = mesh.interfaces
    for i=1:length(ifaces)
      face_i = ifaces[i]
      println(f, (face_i.elementL), " ", Int(face_i.elementR), " ", Int(face_i.faceL), " ", Int(face_i.faceR), " ", Int(face_i.orient))
    end
    close(f)
  end

  if opts["write_boundaries"]
    f = open("boundary_$myrank.dat", "a+")
    bndries = mesh.bndryfaces
    for i=1:length(mesh.bndryfaces)
      bndry_i = mesh.bndryfaces[i]
      println(f, Int(bndry_i.element), " ", Int(bndry_i.face))
    end
    close(f)
  end

  if opts["write_sharedboundaries"]
    for i=1:mesh.npeers
      fname = string("sharedboundaries_", i, "_", myrank, ".dat")
      f = open(fname, "a+")
      bndries_i = mesh.bndries_local[i]
      for j=1:length(bndries_i)
        bndry_j = bndries_i[j]
        println(f, Int(bndry_j.element), " ", Int(bndry_j.face))
      end
      close(f)
    end
  end

  writeVisFiles(mesh, "mesh_complete")

  myrank = mesh.myrank
  f = open("load_balance_$myrank.dat", "a+")
  println(f, mesh.numVert)
  println(f, mesh.numEdge)
  println(f, mesh.numEl)
  println(f, mesh.numInterfaces)
  println(f, mesh.numBoundaryFaces)
  println(f, sum(mesh.peer_face_counts))
  close(f)

  close(mesh.f)
  return mesh
  # could use incomplete initilization to avoid copying arrays
#  return PumiMeshDG2(m_ptr, mshape_ptr, f_ptr, vert_Nptr, edge_Nptr, el_Nptr, numVert, numEdge, numEl, order, numdof, numnodes, dofpernode, bnd_edges_cnt, verts, edges, elements, dofnums_Nptr, bnd_edges_small)
end

 
  
end



function PumiMeshDG2Preconditioning(mesh_old::PumiMeshDG2, sbp::AbstractSBP, opts; 
                                  coloring_distance=0)
# construct pumi mesh for preconditioner residual evaluations
# this operates by copying an existing mesh object (hence not loading a new
#  pumi mesh), and updating a few of its fields
# once there big picture of preconditioning is more clear, there should be 
# options to control these steps.

  println("Creating Pumi mesh for preconditioning")
  mesh = deepcopy(mesh_old)  # create new PumiMesh object, sharing all
                        # mutable fields (like arrays)

  mesh.coloringDistance = coloring_distance
  # create the coloring_Nptr
  el_mshape = getConstantShapePtr(2)
  mesh.coloring_Nptr = createNumberingJ(mesh.m_ptr, "preconditioning coloring", el_mshape, 1)

  if coloring_distance == 2
    numc = colorMesh2(mesh)
    mesh.numColors = numc

    mesh.color_masks = Array(BitArray{1}, numc)  # one array for every color
    mesh.neighbor_colors = zeros(UInt8, 4, mesh.numEl)
    mesh.neighbor_nums = zeros(Int32, 4, mesh.numEl)
    getColors1(mesh, colordata, mesh.color_masks, mesh.neighbor_colors, mesh.neighbor_nums; verify=opts["verify_coloring"] )
    mesh.pertNeighborEls = getPertNeighbors1(mesh)

  elseif coloring_distance == 0  # do a distance-0 coloring
    numc = colorMesh0(mesh)
    @assert numc == 1
    mesh.numColors = numc
    mesh.color_masks = Array(BitArray{1}, numc)
    mesh.neighbor_colors = zeros(UInt8, 0, 0)  # unneeded array for distance-0
    mesh.neighbor_nums = zeros(Int32, 0, 0)  # unneeded for distance-0
    getColors0(mesh, mesh.color_masks)
    mesh.pertNeighborEls = getPertNeighbors0(mesh)

  else
    println(STDERR, "Error: unsupported coloring distance requested")
  end

  # get sparsity information
  # this takes into account the coloring distanc
  mesh.sparsity_bnds = zeros(Int32, 2, mesh.numDof)
  getSparsityBounds(mesh, mesh.sparsity_bnds)
  mesh.sparsity_nodebnds = zeros(Int32, 2, mesh.numNodes)
  getSparsityBounds(mesh, mesh.sparsity_nodebnds, getdofs=false)

  mesh.sparsity_counts = zeros(Int32, 2, mesh.numDof)
  mesh.sparsity_counts_node = zeros(Int32, 2, mesh.numNodes)
  getSparsityCounts(mesh, mesh.sparsity_counts)
  getSparsityCounts(mesh, mesh.sparsity_counts_node, getdofs=false)

  if opts["write_counts"]
    writeCounts(mesh, fname="countsp.txt")
  end


  return mesh

end


# for reinitilizeing after mesh adaptation
function reinitPumiMeshDG2(mesh::PumiMeshDG2)
  # construct pumi mesh by loading the files named
  # dmg_name = name of .dmg (geometry) file to load (use .null to load no file)
  # smb_name = name of .smb (mesh) file to load
  # order = order of shape functions
  # dofpernode = number of dof per node, default = 1

  println("Reinitilizng PumiMeshDG2")

  # create random filenames because they are not used
  smb_name = "a"
  dmg_name = "b"
  order = mesh.order
  dofpernode = mesh.numDofPerNode
  tmp, num_Entities, m_ptr, coordshape_ptr = init2(dmg_name, smb_name, order, load_mesh=false, shape_type=mesh.shape_type) # do not load new mesh
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
  mesh.numBoundaryFaces = bnd_edges_cnt
  mesh.verts = verts  # does this need to be a deep copy?
  mesh.edges = edges
  mesh.elements = elements
#  mesh.boundary_nums = bnd_edges_small

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


