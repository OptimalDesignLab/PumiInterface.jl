module PdePumiInterface
#push!(LOAD_PATH, "/users/creanj/julialib_fork/PUMI.jl")
#push!(LOAD_PATH, "/users/creanj/.julia/v0.4/PDESolver/src/common")
using PumiInterface
using SummationByParts
using ODLCommonTools
import ODLCommonTools.sview
using ArrayViews
using MPI

import ODLCommonTools.sview

include("nodecalc.jl")
#include("bary.jl")
include("options.jl")
#include(joinpath(Pkg.dir("PDESolver"), "src/tools/misc.jl"))

export AbstractMesh,PumiMesh2, PumiMesh2Preconditioning, reinitPumiMesh2, getShapeFunctionOrder, getGlobalNodeNumber, getGlobalNodeNumbers, getNumEl, getNumEdges, getNumVerts, getNumNodes, getNumDofPerNode, getAdjacentEntityNums, getBoundaryEdgeNums, getBoundaryFaceNums, getBoundaryFaceLocalNum, getFaceLocalNum, getBoundaryArray, saveSolutionToMesh, retrieveSolutionFromMesh, retrieveNodeSolution, getAdjacentEntityNums, getNumBoundaryElements, getInterfaceArray, printBoundaryEdgeNums, printdxidx, getdiffelementarea, writeVisFiles, update_coords, commit_coords

export zeroBarArrays, getAllCoordinatesAndMetrics_rev
# Element = an entire element (verts + edges + interior face)
# Type = a vertex or edge or interior face
# 

# the vert, edge, face Nptrs are numberings that map a MeshEntity to its 
# index the the arrays verts, edges, faces, except that they are zero based indices

# internally, all operations and loops are done by looping over nodes 
# in the Pumi order (ie. vertex nodes counter clockwise, edge nodes
# counterclockwise, the face nodes counterclockwise starting in with the
# node nearest the first vertex).
# The external interface uses remaps the nodes into an order specified
# by the nodemaps nodemapSbptoPumi and nodemapPumitoSbp.
# for SBP elements, this amounts to changing the order of the face nodes.
# The external interface is compreised of the arrays:
# dofs, coords, dxidx, jac, sparsity_bnds 
# Care must be taken with the visualization files to map back to the Pumi
# order


export PumiMesh
#abstract AbstractMesh

@doc """
### PdePumiInterface.PumiMeshCG

  This is the abstract supertype for all pumi 2D continuous galerkin type meshes

  The static parameter T1 is the datatype of the mesh variables (coords, dxidx,
  jac).
"""->
abstract PumiMesh2CG{T1} <: AbstractCGMesh{T1}

@doc """
### PdePumiInterface.PumiMeshDG

  This abstract type is the supertype of all 2D Pumi mesh objects for
  discontinuous Galerkin type meshes

  The static parameter T1 is the datatype of the mesh variables (coords, dxidx,
  jac).
"""->
abstract PumiMesh2DG{T1} <: AbstractDGMesh{T1}

"""
### PdePumiInterface.PumiMesh3CG

  This abstract type is the supertype for all 3D Pumi Mesh object for 
  continuous Galrerkin type meshes
"""
abstract PumiMesh3CG{T1} <: AbstractCGMesh{T1}

"""
### PdePumiInterface.PumiMesh3CG

  This abstract type is the supertype of all 3D Pumi mesh object for discontinuous 
  Galerkin type meshes
"""
abstract PumiMesh3DG{T1} <: AbstractDGMesh{T1}

@doc """
### PumiInterface.PumiMesh

  This type is the union of all Pumi mesh types

"""->
typealias PumiMesh{T1} Union{PumiMesh2CG{T1}, PumiMesh2DG{T1}, PumiMesh3CG{T1}, PumiMesh3DG{T1}}

"""
### PumiInterface.PumiMesh2D

  This type is the union of all 2D Pumi mesh types
"""
typealias PumiMesh2D{T1} Union{PumiMesh2CG{T1}, PumiMesh2DG{T1}}

"""
### PumiInterface.PumiMesh3D

  This type is the union of all 3D Pumi mesh types
"""
typealias PumiMesh3D{T1} Union{PumiMesh3CG{T1}, PumiMesh3DG{T1}}

"""
### PumiInterface.PumiMeshCG

  This type is the union of all CG Pumi meshes
"""
typealias PumiMeshCG{T1} Union{PumiMesh2CG{T1}, PumiMesh3CG{T1}}

"""
### PumiInterface.PumiMeshDG

  This type is the union of all DG Pumi meshes
"""
typealias PumiMeshDG{T1} Union{PumiMesh2DG{T1}, PumiMesh3DG{T1}}

"""
  Holds data describing vertices shared between parts

  Fields:
    npeers: number of other processes the current process shared vertices with
    peer_nums: the MPI ranks of the peer processes
    counts: the number of vertices shared with each peer process
            (note: a vertex can be shared with several peer processes)
    vert_nums: Array of arrays.  For each peer process, the inner array contains
               the 1-based numbers of the vertices shared with each process.
               The order of the vertices in each inner array is the same on
               the local and remote process.
    rev_mapping: a dictionary mapping from vert_nums to a
                 Pair{Vector{Cint}, Vector{Int}).
                 The fist vector contains the MPI ranks of all peer processes
                 this vertex is shared with and the second vector is the index
                 of vertex in vert_nums
"""
type VertSharing
  npeers::Int  # number of peers this process shares verts with
  peer_nums::Array{Cint, 1}  # the Part numbers of the peer processes
  counts::Array{Int, 1}  # the number of vertices shared with each peer
  # an array of arrays.  The outer array is of length npeers, the inner arrays
  # are of length counts[i] and hold the 1-based number of each vertex shared
  # with the current peer, in the mutually agreed to ordering
  vert_nums::Array{Array{Int, 1}, 1}

  # mapping from the 1-based number of a vertex to a Pair object containing
  # two arrays, the parts numbers and the corresponding indices of
  # the vertices in the mutually agreed to ordering
  rev_mapping::Dict{Int, Pair{Array{Cint, 1}, Array{Int, 1}}}

end


"""
  This function copies the data fields of one mesh object to another
  mesh object.  Pumi fields are not modified.  The two meshes should
  generally represent the same Pumi mesh.  Any fields not present on both
  meshes are skipped.

  Data fields means: vert_coords, coords*, dxidx*, jac*, nrm*
"""
function copy_data!(dest::PumiMesh, src::PumiMesh)

  # make sure these mesh objects have compatable numbers of things
  @assert src.shape_type == dest.shape_type
  @assert src.order == dest.order
  @assert src.numEl == dest.numEl
  @assert src.numDofPerNode == dest.numDofPerNode

  src_fields = fieldnames(src)
  dest_fields = fieldnames(dest)

  fieldnames_serial = [:vert_coords, :coords_bndry, :coords_interface, :dxidx, 
                :jac, :nrm_face, :nrm_bndry]
  fieldnames_shared = [:coords_sharedface, :dxidx_sharedface, :jac_sharedface,
                       :nrm_sharedface]

  if src.coord_order == 1
    push!(fieldnames_serial, :dxidx_face, :dxidx_bndry, :jac_face, :jac_bndry)
    push!(fieldnames_shared, :dxidx_sharedface, :jac_sharedface)
  end

  for i in fieldnames_serial
    if i in src_fields && i in dest_fields
      f_src = getfield(src, i)
      f_dest = getfield(dest, i)
      copy!(f_dest, f_src)
    end
  end

   for i in fieldnames_shared
     if i in src_fields && i in dest_fields
       f_src = getfield(src, i)
       f_dest = getfield(dest, i)

       for peer=1:src.npeers
        copy!(f_dest[peer], f_src[peer])
      end
    end
  end


  return nothing
end

"""
  Zero out all fields than end with _bar.  Useful for resetting everything
  (including intermediate arrays) after doing a reverse mode evaluation.
"""
function zeroBarArrays(mesh::PumiMesh)

  fnames = fieldnames(mesh)

  for fname in fnames
    fname_str = string(fname)
    if endswith(fname_str, "_bar")
      field = getfield(mesh, fname)

      # treat parallel (Arrays-of-Arrays) separately
      if contains(fname_str, "sharedface")
        for i=1:mesh.npeers
          fill!(field[i], 0.0)
        end
      else
        fill!(field, 0.0)
      end
    end
  end

  return nothing
end

"""
  Registers a finalizer that cleans up the pumi mesh underlying the julia
  mesh object, including all the other meshes (subtriangulated etc.).
  This uses reference counting, so the mesh is only freed if there are no 
  more references to it.
"""
function registerFinalizer(mesh::PumiMesh)

  finalizer(mesh, finalizeMesh)

end

"""
  The finalizer function itself used by registerFinalizer().  If users
  want to free the Pumi datastructure early, they can call this function.
"""
function finalizeMesh(mesh::PumiMesh)

  fnames = fieldnames(mesh)

  # only do the main mesh for now
  if mesh.m_ptr != C_NULL
    popMeshRef(mesh.m_ptr)

    # figure out which other mesh pointers need to be zeroed
    zero_mnew = ( :mnew_ptr in fnames) && mesh.mnew_ptr == mesh.m_ptr
    zero_mexact = ( :mexact_ptr in fnames) && mesh.mexact_ptr == mesh.m_ptr

    mesh.m_ptr = C_NULL
    if zero_mnew
      mesh.mnew_ptr = C_NULL
    end
    if zero_mexact
      mesh.mexact_ptr = C_NULL
    end
  end

  if ( :mnew_ptr in fnames) && mesh.mnew_ptr != C_NULL
    popMeshRef(mesh.mnew_ptr)

    zero_mexact = ( :mexact_ptr in fnames) && mesh.mexact_ptr == mesh.mnew_ptr

    mesh.mnew_ptr = C_NULL
    if zero_mexact
      mesh.mexact_ptr == C_NULL
    end
  end

  if ( :mexect_ptr in fnames) && mesh.mexact_ptr != C_NULL
    popMeshRef(mesh.mexact_ptr)
    mesh.mexact_ptr = C_NULL
  end



  return nothing
end


include("elements.jl")
include("./PdePumiInterface3.jl")
include("PdePumiInterfaceDG.jl")
include("PdePumiInterface3DG.jl")

@doc """
### PumiInterface.PumiMesh2


  This is an implementation of AbstractMesh for a 2 dimensional equation.  
  The constructor for this type extracts all the needed information from Pumi,
  so the solver never needs access to Pumi.

  CG meshes do not support periodic boundary conditions or reverse mode, despite
  the presece of some fields that indicate otherwise.

  Fields:
    m_ptr: a Ptr{Void} to the Pumi mesh
    mnew_ptr: a Ptr{Void} to the subtriangulated mesh used for high
              order visualization, null pointer if not needed

    f_ptr:  a Ptr{Void} to the Pumi apf::Field holding the solution
    vert_NPtr: a pointer to the apf::Numbering object that numbers the 
               vertices (0-based numbering)
    edge_Nptr: a pointer to the Pumi apf::Numbering that numbers the mesh 
               edges (0-based numbering)
    el_Nptr:  a pointer to the Pumi apf::Numbering that numbers the mesh
              elements

    numVert: number of vertices in the mesh
    numEdge: number of edges in the mesh
    numEl: number of elements in the mesh
    order: degree of the shape functions
    numDof: total number of degrees of freedom
    numNodes: number of nodes in the mesh
    numDofPerNode: number of degrees of freedom on each node
    numBoundaryFaces: number of edges on the boundary of the domain
    numInterfaces: number of internal edges (edges not on boundary)
    numPeriodicInterfaces: number of matched interfaces
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

    coloringDistance: the distance k of the distance-k graph coloring used to
                      color the elements (graph vertices are elements and graph
                      edges exist where elements share an edge)
"""->
type PumiMesh2{T1, Tface} <: PumiMesh2CG{T1}   # 2d pumi mesh, triangle only
  m_ptr::Ptr{Void}  # pointer to mesh
  mnew_ptr::Ptr{Void}  # pointer to subtriangulated mesh (high order only)
  mshape_ptr::Ptr{Void} # pointer to mesh's FieldShape
  coordshape_ptr::Ptr{Void}  # pointer to FieldShape of coordinate field
  f_ptr::Ptr{Void} # pointer to apf::field for storing solution during mesh adaptation  (FieldShape = mshape_ptr)
  fnew_ptr::Ptr{Void}  # pointer to field on mnew_ptr
  fnewshape_ptr::Ptr{Void}  # fieldshape of fnew_ptr
  shape_type::Int  #  type of shape functions
  min_node_dist::Float64  # minimum distance between nodes
  min_el_size::Float64 # the minimum size of an element

  vert_Nptr::Ptr{Void}  # numbering of vertices (zero based)
  edge_Nptr::Ptr{Void}  # numbering of edges (zero based)
  face_Nptr::Ptr{Void}  # number of faces (edges), (zero based)
  el_Nptr::Ptr{Void}  # numbering of elements (faces)  (zero based)
  coloring_Nptr::Ptr{Void}  # coloring of mesh for sparse jacobian calculation
  entity_Nptrs::Array{Ptr{Void}, 1}  # [vert_Nptr, edge_Nptr, el_Nptr]
  numVert::Int  # number of vertices in the mesh
  numEdge::Int # number of edges in the mesh
  numFace::Int # alias for numEdge
  numEl::Int  # number of elements (faces)
  order::Int # order of shape functions
  numDof::Int # number of degrees of freedom
  numNodes::Int  # number of nodes
  numDofPerNode::Int  # number of dofs per node
  numBoundaryFaces::Int # number of edges on the exterior boundary
  numInterfaces::Int # number of internal interfaces (including periodic)
  numPeriodicInterfaces::Int  # number of periodic interfaces
  numNodesPerElement::Int  # number of nodes per element
  numFacesPerElement::Int  # number of faces per element
  numNodesPerType::Array{Int, 1}  # number of nodes classified on each vertex, edge, face
  numNodesPerFace::Int  # number of nodes per face
  numEntitiesPerType::Array{Int, 1} # [numVert, numEdge, numEl]
  numTypePerElement::Array{Int, 1}  # number of verts, edges, faces per element
  typeOffsetsPerElement::Array{Int, 1} # the starting index of the vert, edge, and face nodes in an element 
  typeOffsetsPerElement_::Array{Int32, 1}  # Int32 version of above
  nodemapSbpToPumi::Array{UInt8, 1}  # maps nodes of SBP to Pumi order
  nodemapPumiToSbp::Array{UInt8, 1}  # maps nodes of Pumi to SBP order

  # constants needed by Pumi
  el_type::Int  # apf::Type for the elements of the mesh
  face_type::Int # apf::Type for the faces of the mesh

  # parallel info
  comm::MPI.Comm  # MPI Communicator
  myrank::Int  # MPI rank, zero based
  commsize::Int # MPI comm size
  peer_parts::Array{Int, 1}  # array of part numbers that share entities with
                             # the current part
  npeers::Int  # length of the above array

  ref_verts::Array{Float64, 2}  # 2 x 3 array holding the coordinates of 
                                # the vertices of the reference element
                                # in reference coordinates
  dim::Int  # dimension of mesh (2D or 3D)
  isDG::Bool  # is this a DG mesh (always false)
  isInterpolated::Bool # is the field interpolated to the faces
  coloringDistance::Int  # distance between elements of the same color, measured in number of edges
  numColors::Int  # number of colors
  maxColors::Int  # maximum number of colors on any process
  numBC::Int  # number of boundary conditions

  # hold pointers to mesh entities
  verts::Array{Ptr{Void},1}  # holds pointers to mesh vertices
  edges::Array{Ptr{Void},1}  # pointers to mesh edges
  faces::Array{Ptr{Void}, 1}  # alias for edges
  elements::Array{Ptr{Void},1}  # pointers to faces

  # used for high order elements to determine the orientations of edges and 
  # faces (and vertices, too, for consistency)
  elementNodeOffsets::Array{UInt8, 2}
  # truth values if the entity is oriented consistently with the element
  typeNodeFlags::Array{BitArray{2}, 1}
  triangulation::Array{Int32, 2}  # sub triangulation of an element

  nodestatus_Nptr::Ptr{Void}  # node status pointer, status >= 2 -> give it a
                              # node number
  nodenums_Nptr::Ptr{Void}    # node numbering itself
  dofnums_Nptr::Ptr{Void}  # pointer to Numbering of dofs (result of reordering)
#  boundary_nums::Array{Int, 2}  # array of [element number, edgenumber] for each edge on the boundary

  bndry_funcs::Array{BCType, 1}  # array of boundary functors (Abstract type)
  bndry_funcs_revm::Array{BCType_revm, 1}  # reverse mode functors
  bndry_normals::Array{T1, 3}  # array of normals to each face on the boundary
#  bndry_facenums::Array{Array{Int, 1}, 1}  # hold array of faces corresponding to each boundary condition
  bndry_offsets::Array{Int, 1}  # location in bndryfaces where new type of BC starts
                                # and one past the end of the last BC type
				# array has length numBC + 1
  bndry_geo_nums::Array{Array{Int, 1}, 1}  # array of arrays, where the 
                                           # outer array is of length numBC
                                           # and the inner arrays contains
                                           # the geometric edge numbers of this
                                           # BC


  bndryfaces::Array{Boundary, 1}  # store data on external boundary of mesh
  interfaces::Array{Interface, 1}  # store data on internal edges
  interface_normals::Array{T1, 4}  # array of interior face normals, 
                                   # indexing: [norm_component, Left or right, left face node, left face element]

  coords::Array{T1, 3}  # store coordinates of all nodes
  dxidx::Array{T1, 4}  # store scaled mapping jacobian
  dxidx_bar::Array{T1, 4}  # needed for compatability with DG
  jac::Array{T1,2}  # store mapping jacobian output
  jac_bar::Array{T1, 2}  # needed for compatability with DG

  dofs::Array{Int32, 3}  # store dof numbers of solution array to speed assembly
  element_vertnums::Array{Int32, 2}  # map from elements to their vertex numbers
                                     # numVertPerElement x numEl

  sparsity_bnds::Array{Int32, 2}  # store max, min dofs for each dof
  sparsity_nodebnds::Array{Int32, 2}  # store min, max nodes for each node
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

  dof_offset::Int  # local to global dof offset
  sbpface::Tface
  topo::ElementTopology{2}

  facenodes::Array{Int, 2}  # array of numNodesPerFace x numFacesPerElement 
                            # giving the index of the volume node that
                            # corresponds to each face node
                            # this is a temporary hack to keep the PDESolver
                            # test running
 function PumiMesh2(dmg_name::AbstractString, smb_name::AbstractString, order, sbp::AbstractSBP, opts, sbpface; dofpernode=1, shape_type=1, coloring_distance=2)
  # construct pumi mesh by loading the files named
  # dmg_name = name of .dmg (geometry) file to load (use .null to load no file)
  # smb_name = name of .smb (mesh) file to load
  # order = order of shape functions
  # opts: dictionary of options
  # dofpernode = number of dof per node, default = 1
  # shape_type = type of shape functions, 0 = lagrange, 1 = SBP
  # coloring_distance : distance between elements of the same color, where distance is the minimum number of edges that connect the elements, default = 2

  println("\nConstructing PumiMesh2 Object")
  println("  smb_name = ", smb_name)
  println("  dmg_name = ", dmg_name)
  if !MPI.Initialized()
    MPI.Init()
  end
  @assert MPI.Comm_size(MPI.COMM_WORLD) == 1

  mesh = new()
  mesh.isDG = false
  mesh.isInterpolated = false
  mesh.dim = 2
  mesh.npeers = 0
  mesh.numDofPerNode = dofpernode
  mesh.order = order
  mesh.shape_type = shape_type
  mesh.coloringDistance = coloring_distance
  mesh.dof_offset = 0
  mesh.numPeriodicInterfaces = 0
  mesh.topo = ElementTopology2()  # get default topology because it isn't
                                  # important for 2d
  num_Entities, mesh.m_ptr, mesh.mshape_ptr, dim, n_arr = init2(dmg_name, smb_name, order, shape_type=shape_type)
  mesh.coordshape_ptr = mesh.mshape_ptr  # coordinate shape is same as mesh
                                         # field shape for CG
  if dim != mesh.dim
    throw(ErrorException("loaded mesh is not 2 dimensions"))
  end
  mesh.sbpface = sbpface
  if !haskey(ENV, "PDEPUMIINTERFACE_TESTING")
    @assert mesh.order == 1
  end
  mesh.facenodes = Int[1 2 3; 2 3 1]

  mesh.f_ptr = createPackedField(mesh.m_ptr, "solution_field", dofpernode)
  mesh.min_node_dist = minNodeDist(sbp, mesh.isDG)
  mesh.comm = MPI.COMM_WORLD
  mesh.commsize = MPI.Comm_size(MPI.COMM_WORLD)
  mesh.myrank = MPI.Comm_rank(MPI.COMM_WORLD)
  mesh.ref_verts = [0.0 1 0; 0 0 1]

  # count the number of all the different mesh attributes
  mesh.numVert = convert(Int, num_Entities[1])
  mesh.numEdge =convert(Int,  num_Entities[2])
  mesh.numFace = mesh.numEdge
  mesh.numEl = convert(Int, num_Entities[3])
  mesh.numEntitiesPerType = [mesh.numVert, mesh.numEdge, mesh.numEl]
  mesh.numTypePerElement = [3, 3, 1]
  mesh.numFacesPerElement = mesh.numTypePerElement[end-1]

  num_nodes_v = countNodesOn(mesh.mshape_ptr, 0) # on vert
  num_nodes_e = countNodesOn(mesh.mshape_ptr, 1) # on edge
  num_nodes_f = countNodesOn(mesh.mshape_ptr, 2) # on face
  num_nodes_entity = [num_nodes_v, num_nodes_e, num_nodes_f]
  # count numbers of different things per other thing
  # use for bookkeeping
  mesh.numNodesPerType = num_nodes_entity
  mesh.numNodesPerFace = 2*num_nodes_v + num_nodes_e
  mesh.typeOffsetsPerElement = zeros(Int, 4)
  pos = 1
  mesh.typeOffsetsPerElement[1] = pos
  for i=2:4
    pos += mesh.numTypePerElement[i-1]*mesh.numNodesPerType[i-1]
    mesh.typeOffsetsPerElement[i] = pos
  end
  mesh.typeOffsetsPerElement_ = [Int32(i) for i in mesh.typeOffsetsPerElement]
  mesh.numNodesPerElement = mesh.typeOffsetsPerElement[end]-1
  numnodes = 0
  for i=1:3
    numnodes += mesh.numNodesPerType[i]*mesh.numEntitiesPerType[i]
  end
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

  println("about to get entity pointers")
  mesh.verts, mesh.edges, mesh.faces, mesh.elements = getEntityPointers(mesh)
  println("finished getting entity pointers")
  # TODO: check for dangling elements
  checkConnectivity(mesh)


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
    numberNodes(mesh)

#  else
#    println(STDERR, "Error: invalid dof reordering algorithm requested")
#  end

  # do dof numbering
  populateDofNumbers(mesh)
 

  # get entity pointers
  mesh.element_vertnums = getElementVertMap(mesh)

  # get face data
  boundary_nums = getAllFaceData(mesh, opts)

  # use partially constructed mesh object to populate arrays
  mesh.elementNodeOffsets, mesh.typeNodeFlags = getEntityOrientations(mesh)

  getDofNumbers(mesh)  # store dof numbers

  if coloring_distance == 2
    numc = colorMesh2(mesh)
    mesh.numColors = numc
    mesh.maxColors = numc

    mesh.color_masks = Array(BitArray{1}, numc)  # one array for every color
    mesh.neighbor_colors = zeros(UInt8, 4, mesh.numEl)
    mesh.neighbor_nums = zeros(Int32, 4, mesh.numEl)
    getColors1(mesh, mesh.color_masks, mesh.neighbor_colors, mesh.neighbor_nums; verify=opts["verify_coloring"] )
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
    mesh.sparsity_bnds = zeros(Int32, 2, mesh.numDof)
    getSparsityBounds(mesh, mesh.sparsity_bnds)
    mesh.sparsity_nodebnds = zeros(Int32, 2, mesh.numNodes)
    getSparsityBounds(mesh, mesh.sparsity_nodebnds, getdofs=false)
  end


  # TODO: make this a separate option from use_edge_res, make decision
  #       in read_input.jl
  if opts["use_edge_res"]  # if using edge based residual
    mesh.pertNeighborEls_edge = getPertEdgeNeighbors(mesh)
  end

  getCoordinatesAndMetrics(mesh, sbp)  # store coordinates of all nodes into array

  mesh.min_el_size = getMinElementSize(mesh)
  # get face normals
  mesh.bndry_normals = Array(T1, 0, 0, 0)
#  mesh.bndry_normals = Array(T1, 2, sbp.numfacenodes, mesh.numBoundaryFaces)
#  getBoundaryFaceNormals(mesh, sbp, mesh.bndryfaces, mesh.bndry_normals)

  mesh.interface_normals = Array(T1, 0, 0, 0, 0)
#  mesh.interface_normals = Array(T1, 2, 2, sbp.numfacenodes, mesh.numInterfaces)
#  getInternalFaceNormals(mesh, sbp, mesh.interfaces, mesh.interface_normals)

  # create subtriangulated mesh
  createSubtriangulatedMesh(mesh, opts)

  checkFinalMesh(mesh)

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

  if opts["write_edge_vertnums"]
    rmfile("edge_vertnums.dat")
    f = open("edge_vertnums.dat", "a+")
    printEdgeVertNumbers(mesh.edge_Nptr, mesh.vert_Nptr, fstream=f)
    close(f)
  end

  if opts["write_face_vertnums"]
    rmfile("face_vertnums.dat")
    f = open("face_vertnums.dat", "a+")
    printFaceVertNumbers(mesh.el_Nptr, mesh.vert_Nptr, fstream=f)
    close(f)
  end


  if opts["write_boundarynums"]
    rmfile("boundary_nums.dat")
    f = open("boundary_nums.dat", "a+")
    println(f, boundary_nums)
    close(f)
  end


  if opts["write_dxidx"]
    rmfile("dxidx.dat")
    printdxidx("dxidx.dat", mesh.dxidx)
  end


  if opts["write_coords"]
    rmfile("coords.dat")
    writedlm("coords.dat", mesh.coords)
#    printcoords("coords.dat", mesh.coords)
  end

  if opts["write_sparsity"]
    rmfile("sparsity_bnds.dat")
    writedlm("sparsity_bnds.dat", mesh.sparsity_bnds.')
  end

  if opts["write_sparsity_nodebnds"]
    rmfile("sparsity_nodebnds.dat")
    writedlm("sparsity_nodebnds.dat", mesh.sparsity_nodebnds)
  end

  if opts["write_offsets"]
    rmfile("entity_offsets.dat")
    writedlm("entity_offsets.dat", mesh.elementNodeOffsets)
  end

  if opts["write_dofs"]
    rmfile("dofs.dat")
    writedlm("dofs.dat", mesh.dofs)
  end

  if opts["write_counts"]
    writeCounts(mesh)
  end

  writeVisFiles(mesh, "mesh_complete")
  return mesh
  # could use incomplete initilization to avoid copying arrays
#  return PumiMesh2(m_ptr, mshape_ptr, f_ptr, vert_Nptr, edge_Nptr, el_Nptr, numVert, numEdge, numEl, order, numdof, numnodes, dofpernode, bnd_edges_cnt, verts, edges, elements, dofnums_Nptr, bnd_edges_small)
end

 
  
end

include("adapt.jl")
include("bary.jl")
include("coloring.jl")
include("dofnumbering.jl")
include("elements.jl")
include("entities.jl")
include("faces.jl")
include("interpolation.jl")
include("metrics.jl")
include("model.jl")
include("orientation.jl")
include("output.jl")
include("parallel.jl")
include("sparsity.jl")
include("utils.jl")
include("visualization.jl")
include("warp.jl")
include("verify.jl")
include("submesh.jl")

function PumiMesh2Preconditioning(mesh_old::PumiMesh2, sbp::AbstractSBP, opts; 
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
    mesh.maxColors = numc

    mesh.color_masks = Array(BitArray{1}, numc)  # one array for every color
    mesh.neighbor_colors = zeros(UInt8, 4, mesh.numEl)
    mesh.neighbor_nums = zeros(Int32, 4, mesh.numEl)
    getColors1(mesh, mesh.color_masks, mesh.neighbor_colors, mesh.neighbor_nums; verify=opts["verify_coloring"] )
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
    throw(ErrorException("unsupported coloring distance requested"))
  end

  # get sparsity information
  # this takes into account the coloring distanc
  mesh.sparsity_bnds = zeros(Int32, 2, mesh.numDof)
  getSparsityBounds(mesh, mesh.sparsity_bnds)
  mesh.sparsity_nodebnds = zeros(Int32, 2, mesh.numNodes)
  getSparsityBounds(mesh, mesh.sparsity_nodebnds, getdofs=false)

  if opts["write_counts"]
    writeCounts(mesh, fname="countsp.txt")
  end


  return mesh

end



function getMinElementSize(mesh::AbstractMesh)
  # TODO; write a for loop, don't do real(mesh.jac)
  local_min = real(((1./maximum(real(mesh.jac)) )^(1/mesh.dim) )*mesh.min_node_dist)
  global_min = MPI.Allreduce(real(local_min), MPI.MIN, mesh.comm)  #TODO: is the complex part needed in general
  return global_min
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
  tmp, num_Entities, m_ptr, mshape_ptr, dim, n_arr = init2(dmg_name, smb_name, order, load_mesh=false, shape_type=mesh.shape_type) # do not load new mesh
  f_ptr = mesh.f_ptr  # use existing solution field

  numVert = convert(Int, num_Entities[1])
  numEdge =convert(Int,  num_Entities[2])
  numEl = convert(Int, num_Entities[3])

  mesh.vert_Nptr = n_arr[1] #getVertNumbering()
  mesh.edge_Nptr = n_arr[2] #getEdgeNumbering()
  mesh.el_Nptr = n_arr[3] #getFaceNumbering()




  verts = Array(Ptr{Void}, numVert)
  edges = Array(Ptr{Void}, numEdge)
  elements = Array(Ptr{Void}, numEl)
  dofnums_Nptr = mesh.dofnums_Nptr  # use existing dof pointers

  # get pointers to all MeshEntities
  # also initilize the field to zero
#  comps = zeros(dofpernode)
  comps = [1.0, 2, 3, 4]
  it = MeshIterator(mesh.m_ptr, 0)
  for i=1:numVert
    verts[i] = iterate(mesh.m_ptr, it)
  end
  free(mesh.m_ptr, it)

  it = MeshIterator(mesh.m_ptr, 1)
  for i=1:numEdge
    edges[i] = iterate(mesh.m_ptr, it)
  end
  free(mesh.m_ptr, it)

  it = MeshIterator(mesh.m_ptr, it)
  for i=1:numEl
    elements[i] = iterate(mesh.m_ptr, it)
  end
  free(mesh.m_ptr, it)

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
  it = MeshIterator(mesh.m_ptr, it)
  for i=1:numEdge
    edge_i = iterate(mesh.m_ptr, it)
    numFace = countAdjacent(m_ptr, edge_i, 2)  # should be count upward

    if numFace == 1  # if an exterior edge
      faces = getAdjacent(numFace)
      facenum = getNumberJ(mesh.el_Nptr, faces[1], 0, 0) + 1

      bnd_edges_cnt += 1
      bnd_edges[bnd_edges_cnt, 1] = facenum
      bnd_edges[bnd_edges_cnt, 2] = i
    end
  end
  free(mesh.m_ptr, it)

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

  registerFinalizer(mesh)

  writeVtkFiles("mesh_complete", m_ptr)
end


function getShapeFunctionOrder(mesh::PumiMesh2)

return mesh.order
end


end  # end of module



