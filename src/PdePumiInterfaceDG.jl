# PdePumiInterfaceDG.jl: implements AbstractMesh for a 2D Discontinuous
#                        Galerkin mesh


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

# for visualization:
#   2D SBPOmega meshes are subtriangulated and the solution field interpolated
#   onto the new mesh
#   SBPGamma and 3D SBP Omega meshes create a new, linear field on the 
#   existing mesh and interpolate teh solution to it
#
# The fields fnew_ptr, mnew_ptr, and fnewshape_ptr should always refer to the
# field/mesh which is used for visualization, even if it is the same as 
# m_ptr, f_ptr, mshape_ptr or coordshape_ptr

export PumiMeshDG2
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
mutable struct PumiMeshDG2{T1, Tface <: AbstractFace{Float64}} <: PumiMesh2DG{T1}   # 2d pumi mesh, triangle only
  m_ptr::Ptr{Void}  # pointer to mesh
  mnew_ptr::Ptr{Void}  # pointer to mesh used for visualization, which might be
                       # m_ptr or a subtriangulated mesh (high order only)
  mshape_ptr::Ptr{Void} # pointer to the FieldShape of the node field
  coordshape_ptr::Ptr{Void}  # pointer to FieldShape of the coordinate 
                             # field
  f_ptr::Ptr{Void} # pointer to apf::field for storing solution
  fnew_ptr::Ptr{Void}  # pointer to field on mnew_ptr
  fnewshape_ptr::Ptr{Void}  # fieldshape of fnew_ptr
  mexact_ptr::Ptr{Void}  # pointer to subtriangulated mesh used for
                         # exact visualization)
  fexact_ptr::Ptr{Void}  # pointer to field on mexact_ptr
  fexactshape_ptr::Ptr{Void}  # apf::FieldShape of fexact_ptr
  shr_ptr::Ptr{Void}  # pointer to the apf::Sharing object
  normalshr_ptr::Ptr{Void}  # pointer to the NormalSharing object
  shape_type::Int  #  type of shape functions
  subdata::apf.SubMeshData  # if this is a submesh, the submeshdata object,
                        # otherwise NULL
  min_node_dist::Float64  # minimum distance between nodes
  min_el_size::Float64 # size of the smallest element (units of length)
  volume::T1  # volume of entire mesh
  f::IOStream
  vert_Nptr::Ptr{Void}  # numbering of vertices (zero based)
  edge_Nptr::Ptr{Void}  # numbering of edges (zero based)
  face_Nptr::Ptr{Void}  # number of face s(edges), (zero based)
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
  numInterfaces::Int # number of internal interfaces (including periodic)
  numPeriodicInterfaces::Int  # number of periodic interfaces
  numNodesPerElement::Int  # number of nodes per element
  numFacesPerElement::Int  # number of faces (edges) on an element
  numNodesPerType::Array{Int, 1}  # number of nodes classified on each vertex, edge, face
  numNodesPerFace::Int  # number of face nodes
  numEntitiesPerType::Array{Int, 1} # [numVert, numEdge, numEl]
  numTypePerElement::Array{Int, 1}  # number of verts, edges, faces per element
  typeOffsetsPerElement::Array{Int, 1} # the starting index of the vert, edge, and face nodes in an element 
  typeOffsetsPerElement_::Array{Int32, 1}  # Int32 version of above
  nodemapSbpToPumi::Array{UInt8, 1}  # maps nodes of SBP to Pumi order
  nodemapPumiToSbp::Array{UInt8, 1}  # maps nodes of Pumi to SBP order

  # minimal bookkeeping info for coordinate field
  coord_order::Int
  coord_numNodesPerElement::Int
  coord_numNodesPerType::Array{Int, 1} 
  coord_typeOffsetsPerElement::Array{Int, 1}
  coord_numNodes::Int  # total number of coordinate nodes (on this part)
  coord_numNodesPerFace::Int  
  coord_xi::Array{Float64, 2}  # xi coordinates of nodes on reference element
                                # in the Pumi order
                                # dim x coordsnumNodesPerElement
  coord_facexi::Array{Float64, 2}  # like coord_xi, but for the face of an
                                     # element
  coord_nodenums_Nptr::Ptr{Void}  # numbering for nodes of coordinate field,
                                  # number of components = mesh.dim
  geoNums::GeometricDofs  # mapping between coordinate dofs and geometric
                          # dofs

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
  bndry_funcs_revm::Array{BCType_revm, 1}  # reverse mode functors
  bndry_funcs_revq::Array{BCType_revq, 1}
  bndry_normals::Array{T1, 3}  # array of normals to each face on the boundary
#  bndry_facenums::Array{Array{Int, 1}, 1}  # hold array of faces corresponding to each boundary condition
  bndry_offsets::Array{Int, 1}  # location in bndryfaces where new type of BC 
                                # starts
                                # and one past the end of the last BC type
				# array has length numBC + 1

  bndry_geo_nums::Array{Array{Int, 1}, 1}  # array of arrays, where the 
                                           # outer array is of length numBC
                                           # and the inner arrays contains
                                           # the geometric edge numbers of this
                                           # BC

  bndryfaces::Array{Boundary, 1}  # store data on external boundary of mesh
  interfaces::Array{Interface, 1}  # store data on internal edges

  vert_coords::Array{T1, 3}  # dim x coords_numNodesPerElement x numEl array of 
                             # coordinates of vertices of each element,
                             # ordered vertices, then edge nodes, then face
                             # nodes
  vert_coords_bar::Array{T1, 3}  # adjoint part
  coords::Array{T1, 3}  # store coordinates of all nodes

  coords_bndry::Array{T1, 3}  # store coordinates of nodes on boundary,
                              # 2 x numFaceNodes x numBoundaryFaces
  coords_bndry_bar::Array{T1, 3}  # adjoint part

  coords_interface::Array{T1, 3} # store coordinates of nodes on interfaces
                                 # 2 x numFaceNodes x numInterfaces

  coords_sharedface::Array{Array{T1, 3}, 1}  # coordinates of shared interface nodes
  dxidx::Array{T1, 4}  # store scaled mapping jacobian
  dxidx_bar::Array{T1, 4}
  dxidx_face::Array{T1, 4} # store scaled mapping jacobian at face nodes
  dxidx_face_bar::Array{T1, 4}
                           # 2 x 2 x numfacenodes x numInterfaces
  dxidx_sharedface::Array{Array{T1, 4}, 1}  # array of arrays for dxidx
  dxidx_sharedface_bar::Array{Array{T1, 4}, 1}
                                            # on shared edges
  dxidx_bndry::Array{T1, 4} # store scaled mapping jacobian at boundary nodes,
                            # similar to dxidx_face
  dxidx_bndry_bar::Array{T1, 4}

  jac::Array{T1,2}  # store mapping jacobian output
  jac_bar::Array{T1, 2}
  jac_face::Array{T1,2}  # store jacobian determanent at face nodes
                         # numfacenodes x numInterfaces
  jac_face_bar::Array{T1, 2}
  jac_sharedface::Array{Array{T1, 2}, 1}  # array of arrays for shared
                                          # edge jacobian determinent
  jac_sharedface_bar::Array{Array{T1, 2}, 1}

  jac_bndry::Array{T1, 2} # store jacobian determinant at boundry nodes
                          # similar to jac_bndry
  jac_bndry_bar::Array{T1, 2}

  nrm_bndry::Array{T1, 3}  # dim x numfacenodes x numBoundaryFaces array holding
                           # normal vector in x-y to each face node  on the boundary
  nrm_bndry_bar::Array{T1, 3}

  nrm_face::Array{T1, 3}  # dim x numfacenodes x numInterfaces array holding
                          # normal vector to each face node on elementL of an
                          # interface
  nrm_face_bar::Array{T1, 3}

  nrm_sharedface::Array{Array{T1, 3}, 1}  # array of arrays.  Outer array is of
                          # lench npeers.  Inner arrays are of size dim x
                          # numfacenodes x number of faces shared with current
                          # peer
  nrm_sharedface_bar::Array{Array{T1, 3}, 1}

  dof_offset::Int  # local to global offset for dofs
  dofs::Array{Int, 3}  # store dof numbers of solution array to speed assembly
  element_vertnums::Array{Int32, 2}  # map from elements to their vertex numbers
                                  # numVertPerElement x numEl
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

  interp_op::Array{Float64, 2}  # coord_numNodesPerElement x numNodesPerEl
                                # matrix that interpolates
                                # solution values to the vertices
  interp_op2::Array{Float64, 2}  # numNodesPerEl x coord_numNodesPerElement
                                 # array that interpolates from the
                                 # coordinate field to the solution field
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

  remote_metrics::Array{RemoteMetrics{T1}, 1}  # metric information for
                                                 # remote elements
  remote_metrics_bar::Array{RemoteMetrics{T1}, 1}

  sbpface::Tface  # SBP object needed to do interpolation
  topo::ElementTopology{2}  # topology of the SBP element
  topo_pumi::ElementTopology{2}  # topology of the pumi element
  # add field that maps back and forth?

  vert_sharing::VertSharing
  coordscatter::ScatterData{T1}  # abstract type, used for coords3DTo1D

  # temporarily allow nested meshes for the staggered grid work
  mesh2::AbstractMesh
  sbp2::AbstractSBP

  # interpolation operators
  # these are usually stored in the solution mesh but not in the flux mesh
  I_S2F::Matrix{Float64}  # interpolation operator from solution to flux grid,
                          # numNodesPerElement_f x numNodesPerElement_s
  I_S2FT::Matrix{Float64}  # transpose of above
  I_F2S::Matrix{Float64}  # interpolation operator from flux grid to solution
                          # grid, numNodesPerElement_s x numNodesPerElement_f
  I_F2ST::Matrix{Float64} # transpose of above

  fields::AttachedData
  numberings::AttachedData

  """
    This inner constructor loads a Pumi mesh from files and sets a few
    essential fields that must be consistent with how the mesh was loaded

    **Inputs**

     * dmg_name: Pumi geometry file name, including extension
     * smb_name: Pumi mesh file name, including extension, not including
                 the file number (ex. abc.smb not abc0.smb on process 0)

  **Keywords**
   
   * shape_type: the integer identifying the coordinate field shape,
                 default=-1, meaning keep existing coordinate field
   * order: the coordinate field order (only used if the coordinate field
            is changed from what is stored in the mesh, default 1
   * comm: MPI communicator the mesh should be defined on, must be the same
           number of processes as the mesh file is partitioned into.

  **Outputs**

   * mesh

  """
  function PumiMeshDG2{T1, Tface}(dmg_name::AbstractString,
                       smb_name::AbstractString;
                       order::Integer=1, shape_type::Integer=-1,
                       comm=MPI.COMM_WORLD) where {T1, Tface <: AbstractFace{Float64}}
    # this constructor does all the low-level Pumi initialization
    # shape_type should be the *coordinate* shape type
    # order should be the shape function order

    if !MPI.Initialized()
      MPI.Init()
    end


    mesh = new{T1, Tface}()
    mesh.isDG = true
    mesh.dim = 2

    mesh.comm = comm
    mesh.topo_pumi = ElementTopology{2}(apf.tri_edge_verts.')
    mesh.myrank = MPI.Comm_rank(mesh.comm)
    mesh.commsize = MPI.Comm_size(mesh.comm)
    myrank = mesh.myrank
    mesh.subdata = apf.SubMeshData(C_NULL)
    mesh.fields = AttachedData()
    mesh.numberings = AttachedData()

    if myrank == 0
      println("\nConstructing PumiMeshDG2 Object")
      println("  smb_name = ", smb_name)
      println("  dmg_name = ", dmg_name)
    end

    mesh.m_ptr, dim = apf.loadMesh(dmg_name, smb_name, order, 
                               shape_type=shape_type)
    apf.pushMeshRef(mesh.m_ptr)
    recordAllFields(mesh)
    recordAllNumberings(mesh)

    if dim != mesh.dim
      throw(ErrorException("loaded mesh is not 2 dimensions"))
    end


    return mesh
  end # end inner constructor

  """
    Returns uninitialized PumiMeshDG2 object.
  """
  function PumiMeshDG2{T1, Tface}() where {T1, Tface}
    return new{T1, Tface}()
  end

end  # end PumiMeshDG2 declaration

  
"""
  This outer constructor populates the essential fields of the mesh
  that are populated when a mesh is loaded from a file, but instead
  of loading a mesh from a file, it copies the fields from the old mesh.

  The mesh.m_ptr field is not populated because this function is used
  for creating submeshes.
"""
function PumiMeshDG2(old_mesh::PumiMeshDG2{T},
                     sbpface::Tface=old_mesh.sbpface) where {T, Tface}

  mesh = PumiMeshDG2{T, Tface}()  # get uninitialized object

  # set essential fields from old_mesh
  mesh.isDG = true
  mesh.dim = 2
  mesh.comm = old_mesh.comm
  mesh.topo_pumi = old_mesh.topo_pumi
  mesh.sbpface = sbpface
  mesh.myrank = old_mesh.myrank
  mesh.commsize = old_mesh.commsize
  mesh.fields = AttachedData()
  mesh.numberings = AttachedData()

  return mesh
end


"""
   This outer constructor loads a mesh from a file, according to the options
   specified in the options dictionary
"""
function PumiMeshDG2(::Type{T}, sbp::AbstractSBP, opts, 
                     sbpface::Tface; dofpernode=1, shape_type=2,
                     comm=MPI.COMM_WORLD) where {T, Tface}

  set_defaults(opts)  # get default arguments

  # distinguish between the coordinate field and the solution field
  coord_shape_type = -1  # keep coordinate field already present
  coord_order = 1  # this is only used if the shape type is changed

  # unpack options keys
  # TODO: doc this
  smb_name = opts["smb_name"]
  dmg_name = opts["dmg_name"]

  mesh = PumiMeshDG2{T, Tface}(dmg_name, smb_name, shape_type=coord_shape_type,
                             order=coord_order)
  mesh.sbpface = sbpface
  finishMeshInit(mesh, sbp, opts, dofpernode=dofpernode,
                 shape_type=shape_type)

  return mesh 
end  # end outer constructor


"""
  This outer constructor makes a submesh from an existing mesh. Mesh edges that
  were previously on the interior will be classified on a new geometric edge,
  and a new options dictionary will be created that has a new boundary condition.
  Serial meshes only (for now)
  No periodic BCs (for now)

  **Inputs**

   * old_mesh: the existing mesh
   * sbp: the SBP operator used to create old_mesh
   * opts_old: the options dictionary used to create old_mesh
   * newbc_name: the name of the new boundary condition that will be created
   * el_list: vector of elements to include on the new mesh (preferably a
              vector of Cints)

  **Outputs**

   * mesh: the new mesh object, fully initialized
   * opts: the new options dictionary
"""
function PumiMeshDG2(old_mesh::PumiMeshDG2{T, Tface}, sbp, opts_old,
                    newbc_name::AbstractString, el_list::AbstractVector) where {T, Tface}

  @assert length(el_list) > 0

  if old_mesh.commsize != 1
    throw(ErrorException("Submesh not supported in parallel"))
  end

  if opts_old["reordering_algorithm"] == "adjacency"
    throw(ErrorException("numberNodesWindy not supported for submesh"))
  end

  # get new mesh with essential fields populated

  mesh = PumiMeshDG2(old_mesh)

  # construct new mesh
  if eltype(el_list) != Cint
    _el_list = Array{Cint}(length(el_list))
    copy!(_el_list, el_list)
  else
    _el_list = el_list
  end

  mesh.subdata = apf.createSubMesh(old_mesh.m_ptr, old_mesh.entity_Nptrs, _el_list)
  mesh.m_ptr = apf.getNewMesh(mesh.subdata)
  apf.pushMeshRef(mesh.m_ptr)
  new_geo = apf.getGeoTag(mesh.subdata)

  # create new options dictionary, updating BCs
  opts = deepcopy(opts_old)
  updateBCs(opts, new_geo, newbc_name) 

  finishMeshInit(mesh, sbp, opts, dofpernode=old_mesh.numDofPerNode,
                 shape_type=old_mesh.shape_type)


  return mesh, opts
end

"""
  This constructor creates a new mesh object after mesh adaptation has
  run.

  **Inputs**

   * old_mesh: old mesh object, mesh.m_ptr points to the adapted mesh (
               which means that mesh.m_ptr and the mesh object are inconsistent
               while this constructor runs
   * sbp: SBP operator
   * opts: options dictionary
   * sbpface: an `AbstractFace` object, defaults to `old_mesh.sbpface`

  **Outputs**

   * mesh: new mesh object
"""
function PumiMeshDG2(old_mesh::PumiMeshDG2{T}, sbp, opts, sbpface::Tface=old_mesh.sbpface) where {T, Tface}

  mesh = PumiMeshDG2(old_mesh, sbpface)
  mesh.m_ptr = old_mesh.m_ptr
  apf.pushMeshRef(mesh.m_ptr)
  mesh.subdata = apf.SubMeshData(C_NULL)
  attachOrigFields(mesh, old_mesh.fields.orig)
  attachOrigNumberings(mesh, old_mesh.numberings.orig)

  finishMeshInit(mesh, sbp, opts, dofpernode=old_mesh.numDofPerNode,
                 shape_type=old_mesh.shape_type)

  return mesh
end

"""
  This function finishes initializing the mesh object.  This does the
  bulk of the work, and is used by most of the constructors.  The following
  fields must already be populated:

   * isDG
   * dim
   * comm
   * topo_pumi
   * myrank
   * commsize
   * m_ptr

  **Inputs**

   * mesh: mesh with fields initialized.  All fields will be initialized on
           exit
   * sbp: the SBP operator to be used with the mesh (must match shape_type)
   * opts: options dictionary, the keys "order" and "coloring_distance"
           are used
  
  **Keywords**

   * dofpernode: number of dofs on each node, default 1
   * shape_type: integer describing solution field

"""
function finishMeshInit(mesh::PumiMeshDG2{T1},  sbp::AbstractSBP, opts; dofpernode=1, shape_type=2) where T1
  # construct pumi mesh by loading the files named
  # dmg_name = name of .dmg (geometry) file to load (use .null to load no file)
  # smb_name = name of .smb (mesh) file to load
  # order = order of shape functions
  # opts: dictionary of options
  # dofpernode = number of dof per node, default = 1
  # shape_type = type of shape functions, 0 = lagrange, 1 = SBP, 2 = SBP DG1
  #              3 = DG2
  # coloring_distance : distance between elements of the same color, where distance is the minimum number of edges that connect the elements, default = 2


  # call inner constructor here

  # unpack options keys
  # TODO: doc this
  field_shape_type = shape_type
  order = opts["order"]  # SBP degree
  coloring_distance = opts["coloring_distance"]
  myrank = mesh.myrank
  sbpface = mesh.sbpface

  # for now, leave shape_type has an argument, but eventually this function
  # should take in the SBP operator name and return a shape_type

  mesh.coordshape_ptr, num_Entities, n_arr = apf.initMesh(mesh.m_ptr)
  for n in n_arr[1:(mesh.dim + 1)]
    attachUserNumbering(mesh, n)
  end
  #  num_Entities, mesh.m_ptr, mesh.coordshape_ptr, dim, n_arr = apf.init2(dmg_name, smb_name, mesh_order, shape_type=coord_shape_type)

  mesh.coloringDistance = coloring_distance
  mesh.shape_type = shape_type
#  mesh.interp_op = interp_op
  mesh.ref_verts = [0.0 1 0; 0 0 1]  # ???
  mesh.numNodesPerFace = sbpface.numnodes
  mesh.topo = ElementTopology2() # create default topology because it isn't
                                 # important for 2D
  mesh.numDofPerNode = dofpernode
  mesh.order = order
  
  mesh.mshape_ptr = apf.getFieldShape(field_shape_type, order, mesh.dim)

  #TODO: is mesh.f_ptr used for anything? visualization uses mesh.f_new,
  # which has the same fieldshape?
  mesh.f_ptr = apf.createPackedField(mesh, "solution_field", dofpernode, mesh.mshape_ptr)

  mesh.shr_ptr = apf.getSharing(mesh.m_ptr)
  mesh.normalshr_ptr = apf.getNormalSharing(mesh.m_ptr)
  # count the number of all the different mesh attributes
  mesh.numVert = convert(Int, num_Entities[1])
  mesh.numEdge =convert(Int,  num_Entities[2])
  mesh.numFace = mesh.numEdge
  mesh.numEl = convert(Int, num_Entities[3])
  mesh.numEntitiesPerType = [mesh.numVert, mesh.numEdge, mesh.numEl]
  mesh.numTypePerElement = [3, 3, 1]
  mesh.el_type = apf.TRIANGLE
  mesh.face_type = apf.EDGE
  mesh.numFacesPerElement = mesh.numTypePerElement[end-1]


  num_nodes_v = apf.countNodesOn(mesh.mshape_ptr, 0)  # number of nodes on a vertex
  num_nodes_e = apf.countNodesOn(mesh.mshape_ptr, 1) # on edge
  num_nodes_f = apf.countNodesOn(mesh.mshape_ptr, 2) # on face
  num_nodes_entity = [num_nodes_v, num_nodes_e, num_nodes_f]
  mesh.numNodesPerType = num_nodes_entity


  checkTopologyConsistency(mesh.topo, mesh.topo_pumi)

  mesh.f = open("meshlog_$myrank.dat", "w")

  # temporary testing of SBP-Gamma DG
#  if sbp.numfacenodes == 0
  mesh.isInterpolated = true
#  else
#    mesh.isInterpolated = false
#    # leave mesh.sbpface undefined - bad practice
#  end

 
  # create the adjoint part of the coordinate field
  # 3 components always, even in 2D for consistency with Pumi's coordinate field

  # create the solution field
  mesh.min_node_dist = minNodeDist(sbp, mesh.isDG)

  # count numbers of different things per other thing
  # use for bookkeeping
  mesh.typeOffsetsPerElement = zeros(Int, 4)
  pos = 1
  mesh.typeOffsetsPerElement[1] = pos
  for i=2:4
    pos += mesh.numTypePerElement[i-1]*mesh.numNodesPerType[i-1]
    mesh.typeOffsetsPerElement[i] = pos
  end
  mesh.numNodesPerType, mesh.typeOffsetsPerElement = getNodeInfo(mesh.mshape_ptr, mesh.dim, mesh.numTypePerElement)

  # get coordinate field info
  mesh.coord_numNodesPerType, mesh.coord_typeOffsetsPerElement = getNodeInfo(mesh.coordshape_ptr, mesh.dim, mesh.numTypePerElement)
  mesh.coord_numNodesPerElement = mesh.coord_typeOffsetsPerElement[end] - 1
  mesh.coord_order = apf.getOrder(mesh.coordshape_ptr)
  mesh.coord_xi = getXiCoords(mesh.coord_order, mesh.dim)
  mesh.coord_facexi = getXiCoords(mesh.coord_order, mesh.dim-1)
  mesh.coord_numNodesPerFace = size(mesh.coord_facexi, 2)

  # build interpolation operator
  ref_coords = baryToXY(mesh.coord_xi, sbp.vtx)
  mesh.interp_op = SummationByParts.buildinterpolation(sbp, ref_coords)
  mesh.interp_op2 = SummationByParts.buildinterpolation(ref_coords, calcnodes(sbp), mesh.coord_order)
  
  mesh.typeOffsetsPerElement_ = [Int32(i) for i in mesh.typeOffsetsPerElement]

  mesh.numNodesPerElement = mesh.typeOffsetsPerElement[end] - 1
  numnodes = mesh.numNodesPerElement*mesh.numEl
  mesh.numNodes = numnodes      # we assume there are no non-free nodes/dofs
  mesh.numDof = numnodes*dofpernode

  # get nodemaps
  mesh.nodemapSbpToPumi, mesh.nodemapPumiToSbp = getNodeMaps(order, shape_type, mesh.numNodesPerElement, mesh.dim, mesh.isDG)

 
  # get pointers to mesh entity numberings
  mesh.vert_Nptr = n_arr[1] #getVertNumbering()
  mesh.edge_Nptr = n_arr[2] #getEdgeNumbering()
  mesh.face_Nptr = n_arr[2] #mesh.edge_Nptr
  mesh.el_Nptr = n_arr[3] #getFaceNumbering()
  mesh.entity_Nptrs = [mesh.vert_Nptr, mesh.edge_Nptr, mesh.el_Nptr]

  # create the coloring_Nptr
  el_mshape = apf.getConstantShapePtr(2)
  mesh.coloring_Nptr = apf.createNumberingJ(mesh, "coloring", 1, el_mshape)

  # create node status numbering (node, not dof)
  mesh.nodestatus_Nptr = apf.createNumberingJ(mesh, "dof status", 
                                              1, mesh.mshape_ptr)   
  # create node numbering
  mesh.nodenums_Nptr = apf.createNumberingJ(mesh, "reordered node numbers",
                                            1, mesh.mshape_ptr)

  # create dof numbering
  mesh.dofnums_Nptr = apf.createNumberingJ(mesh, "reordered dof numbers", 
                                           dofpernode, mesh.mshape_ptr)

  # coordinate node numbering
  mesh.coord_nodenums_Nptr = apf.createNumberingJ(mesh, "coord node numbers",
                                                  mesh.dim, mesh.coordshape_ptr)

  xiNums_Nptr = apf.createNumberingJ(mesh, "geometric dof numbers",
                                     mesh.dim, mesh.coordshape_ptr)

  # get entity pointers
#  println("about to get entity pointers")
  mesh.verts, mesh.edges, mesh.faces, mesh.elements = getEntityPointers(mesh)
#  println("finished getting entity pointers")
  checkConnectivity(mesh)


  # populate node status numbering
  populateNodeStatus(mesh)


  # do node reordering
  coord_numnodes = 0
  for i=1:(mesh.dim+1)
    coord_numnodes += mesh.coord_numNodesPerType[i]*mesh.numEntitiesPerType[i]
  end
  mesh.coord_numNodes = coord_numnodes

 if opts["reordering_algorithm"] == "adjacency"
    start_coords = opts["reordering_start_coords"]
    numberNodesWindy(mesh, start_coords)
    # tell the algorithm there is only 1 dof per node because we only
    # want to label nodes

    # coordinate node numbering
    apf.reorder(mesh.m_ptr, mesh.dim*mesh.coord_numNodes, mesh.dim, 
            C_NULL, mesh.coord_nodenums_Nptr, C_NULL, 
	    start_coords)

    numXiDof = apf.reorderXi(mesh.m_ptr, xiNums_Nptr, start_coords)
 elseif opts["reordering_algorithm"] == "default"
#    println("about to number nodes")
    numberNodesElement(mesh)

    # coordinate node numbering
    # have to use getPoint here because mesh.geoNums hasn't been assigned yet
    start_coords = zeros(3)
    apf.getPoint(mesh.m_ptr, mesh.verts[1], 0, start_coords)  
    apf.reorder(mesh.m_ptr, mesh.dim*mesh.coord_numNodes, mesh.dim, 
            C_NULL, mesh.coord_nodenums_Nptr, C_NULL, 
	    start_coords)

    numXiDof = apf.reorderXi(mesh.m_ptr, xiNums_Nptr, start_coords)
  else
    throw(ErrorException("invalid dof reordering algorithm requested"))
  end

  mesh.geoNums = GeometricDofs(mesh.m_ptr, mesh.coord_nodenums_Nptr,
                               xiNums_Nptr, mesh.dim*mesh.coord_numNodes,
                               numXiDof)



#  println("finished numbering nodes")

#  println("about to number dofs")
  # do dof numbering
  populateDofNumbers(mesh)
#  println("finished numbering dofs")
 

  mesh.element_vertnums = getElementVertMap(mesh)
  boundary_nums = getAllFaceData(mesh, opts)


  # get array of all boundary mesh edges in the same order as in mesh.bndry_faces
  # use partially constructed mesh object to populate arrays
#  println("about to get entity orientations")
  mesh.elementNodeOffsets, mesh.typeNodeFlags = getEntityOrientations(mesh)
#  println("finished getting entity orientations")


  # start parallel initializiation
  colordata = getParallelInfo(mesh)
  mesh.vert_sharing = getVertexParallelInfo(mesh)

#  println("about to get degree of freedom numbers")
  getDofNumbers(mesh)  # store dof numbers
#  println("finished getting degree of freedom numbers")
  mesh.coordscatter = initSendToOwner(mesh, mesh.coordshape_ptr, (mesh.dim,))



  MPI.Barrier(mesh.comm)
  if coloring_distance == 2
    numc = colorMesh2(mesh, colordata)
    mesh.numColors = numc
    mesh.maxColors = MPI.Allreduce(numc, MPI.MAX, mesh.comm)
    mesh.color_masks = Array{BitArray{1}}(numc)  # one array for every color
    mesh.neighbor_colors = zeros(UInt8, 4, mesh.numEl)
    mesh.neighbor_nums = zeros(Int32, 4, mesh.numEl)
    cnt, mesh.shared_element_colormasks = getColors1(mesh, colordata, mesh.color_masks, 
                                mesh.neighbor_colors, mesh.neighbor_nums; verify=opts["verify_coloring"] )
    mesh.pertNeighborEls = getPertNeighbors1(mesh)

  elseif coloring_distance == 0  # do a distance-0 coloring
    numc = colorMesh0(mesh)
    @assert numc == 1
    mesh.numColors = numc
    mesh.maxColors = numc
    mesh.color_masks = Array{BitArray{1}}(numc)
    mesh.neighbor_colors = zeros(UInt8, 0, 0)  # unneeded array for distance-0
    mesh.neighbor_nums = zeros(Int32, 0, 0)  # unneeded for distance-0
    getColors0(mesh, mesh.color_masks)
    mesh.pertNeighborEls = getPertNeighbors0(mesh)

  else
    error("Error: unsupported coloring distance requested")
  end

  # get sparsity information
  # this takes into account the coloring distance
  #TODO: make getting sparsity bounds faster
  if opts["run_type"] != 1  # no a rk4 run
    mesh.sparsity_bnds = zeros(Int32, 0, 0)
#    mesh.sparsity_bnds = zeros(Int32, 2, mesh.numDof)
#    @time getSparsityBounds(mesh, mesh.sparsity_bnds)
#    mesh.sparsity_nodebnds = zeros(Int32, 2, mesh.numNodes)
#    @time getSparsityBounds(mesh, mesh.sparsity_nodebnds, getdofs=false)
#    println("finished getting sparsity bounds")

    mesh.sparsity_counts = zeros(Int32, 2, mesh.numDof)
    mesh.sparsity_counts_node = zeros(Int32, 2, mesh.numNodes)
    getSparsityCounts(mesh, mesh.sparsity_counts)
    getSparsityCounts(mesh, mesh.sparsity_counts_node, getdofs=false)


  end


  # TODO: make this a separate option from use_edge_res, make decision
  #       in read_input.jl
  if opts["use_edge_res"]  # if using edge based residual
    mesh.pertNeighborEls_edge = getPertEdgeNeighbors(mesh)
  end

  getAllCoordinatesAndMetrics(mesh, sbp, opts)

  createSubtriangulatedMesh(mesh, opts)

  printStats(mesh)

  # write data if requested
  myrank = mesh.myrank
  if opts["write_edge_vertnums"]
    rmfile("edge_vertnums_$myrank.dat")
    f = open("edge_vertnums_$myrank.dat", "a+")
    apf.printEdgeVertNumbers(mesh.edge_Nptr, mesh.vert_Nptr, fstream=f)
    close(f)
  end

  if opts["write_face_vertnums"]
    rmfile("face_vertnums_$myrank.dat")
    f = open("face_vertnums_$myrank.dat", "a+")
    apf.printFaceVertNumbers(mesh.el_Nptr, mesh.vert_Nptr, fstream=f)
    close(f)
  end

  if opts["write_boundarynums"]
    writedlm("boundary_nums_$myrank.dat", boundary_nums)
  end

  if opts["write_dxidx"]
    rmfile("dxidx_$myrank.dat")
    writedlm("dxidx_$myrank.dat", mesh.dxidx)
  end

  if opts["write_jac"]
    rmfile("jac_$myrank.dat")
    writedlm("jac_$myrank.dat", mesh.jac)
  end


  if opts["write_coords"]
    rmfile("coords_$myrank.dat")
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

  checkFinalMesh(mesh)

  myrank = mesh.myrank
  f = open("load_balance_$myrank.dat", "a+")
  println(f, mesh.numVert)
  println(f, mesh.numEdge)
  println(f, mesh.numEl)
  println(f, mesh.numInterfaces)
  println(f, mesh.numBoundaryFaces)
  println(f, sum(mesh.peer_face_counts))
  close(f)

  registerFinalizer(mesh)

  #close(mesh.f)
  return mesh
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
  el_mshape = apf.getConstantShapePtr(2)
  mesh.coloring_Nptr = apf.createNumberingJ(mesh, "preconditioning coloring",
                                            1, el_mshape)

  if coloring_distance == 2
    numc = colorMesh2(mesh)
    mesh.numColors = numc

    mesh.color_masks = Array{BitArray{1}}(numc)  # one array for every color
    mesh.neighbor_colors = zeros(UInt8, 4, mesh.numEl)
    mesh.neighbor_nums = zeros(Int32, 4, mesh.numEl)
    getColors1(mesh, colordata, mesh.color_masks, mesh.neighbor_colors, mesh.neighbor_nums; verify=opts["verify_coloring"] )
    mesh.pertNeighborEls = getPertNeighbors1(mesh)

  elseif coloring_distance == 0  # do a distance-0 coloring
    numc = colorMesh0(mesh)
    @assert numc == 1
    mesh.numColors = numc
    mesh.color_masks = Array{BitArray{1}}(numc)
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

