module PdePumiInterface
push!(LOAD_PATH, "/users/creanj/julialib_fork/PUMI.jl")
push!(LOAD_PATH, "/users/creanj/.julia/v0.4/PDESolver/src/common")
using PumiInterface
using SummationByParts
using PDESolverCommon
using ArrayViews
#include(joinpath(Pkg.dir("PDESolver"), "src/tools/misc.jl"))

export AbstractMesh,PumiMesh2, reinitPumiMesh2, getElementVertCoords, getShapeFunctionOrder, getGlobalNodeNumber, getGlobalNodeNumbers, getNumEl, getNumEdges, getNumVerts, getNumNodes, getNumDofPerNode, getAdjacentEntityNums, getBoundaryEdgeNums, getBoundaryFaceNums, getBoundaryEdgeLocalNum, getEdgeLocalNum, getBoundaryArray, saveSolutionToMesh, retrieveSolutionFromMesh, retrieveNodeSolution, getAdjacentEntityNums, getNumBoundaryElements, getInterfaceArray, printBoundaryEdgeNums, printdxidx, getdiffelementarea, writeVisFiles

# Element = an entire element (verts + edges + interior face)
# Type = a vertex or edge or interior face
# 

# the vert, edge, face Nptrs are numberings that map a MeshEntity to its 
# index the the arrays verts, edges, faces, except that they are zero based indices
export PumiMesh
#abstract AbstractMesh
abstract PumiMesh{T1} <: AbstractMesh{T1}
include("./PdePumiInterface3.jl")
type PumiMesh2{T1} <: PumiMesh{T1}   # 2d pumi mesh, triangle only
  m_ptr::Ptr{Void}  # pointer to mesh
  mshape_ptr::Ptr{Void} # pointer to mesh's FieldShape
  f_ptr::Ptr{Void} # pointer to apf::field for storing solution during mesh adaptation
  shape_type::Int  #  type of shape functions

  vert_Nptr::Ptr{Void}  # numbering of vertices (zero based)
  edge_Nptr::Ptr{Void}  # numbering of edges (zero based)
  el_Nptr::Ptr{Void}  # numbering of elements (faces)  (zero based)
  coloring_Nptr::Ptr{Void}  # coloring of mesh for sparse jacobian calculation
  entity_Nptrs::Array{Ptr{Void}, 1}  # [vert_Nptr, edge_Nptr, el_Nptr]
  numVert::Int  # number of vertices in the mesh
  numEdge::Int # number of edges in the mesh
  numEl::Int  # number of elements (faces)
  order::Int # order of shape functions
  numDof::Int # number of degrees of freedom
  numNodes::Int  # number of nodes
  numDofPerNode::Int  # number of dofs per node
  numBoundaryEdges::Int # number of edges on the exterior boundary
  numInterfaces::Int # number of internal interfaces
  numNodesPerElement::Int  # number of nodes per element
  numNodesPerType::Array{Int, 1}  # number of nodes classified on each vertex, edge, face
  numEntitiesPerType::Array{Int, 1} # [numVert, numEdge, numEl]
  numTypePerElement::Array{Int, 1}  # number of verts, edges, faces per element
  typeOffsetsPerElement::Array{Int, 1} # the starting index of the vert, edge, and face nodes in an element 
  coloringDistance::Int  # distance between elements of the same color, measured in number of edges
  numColors::Int  # number of colors
  numBC::Int  # number of boundary conditions

  # hold pointers to mesh entities
  verts::Array{Ptr{Void},1}  # holds pointers to mesh vertices
  edges::Array{Ptr{Void},1}  # pointers to mesh edges
  elements::Array{Ptr{Void},1}  # pointers to faces

  # used for high order elements to determine the orientations of edges and 
  # faces (and vertices, too, for consistency)
  elementNodeOffsets::Array{Uint8, 2}
  # truth values if the entity is oriented consistently with the element
  typeNodeFlags::Array{BitArray{2}, 1}

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
  dxidx::Array{T1, 4}  # store scaled mapping jacobian
  jac::Array{T1,2}  # store mapping jacobian output

  dofs::Array{Int32, 3}  # store dof numbers of solution array to speed assembly
  sparsity_bnds::Array{Int32, 2}  # store max, min dofs for each dof
  color_masks::Array{BitArray{1}, 1}  # array of bitarray masks used to control element perturbations when forming jacobian, number of arrays = number of colors
  neighbor_colors::Array{UInt8, 2}  # 4 by numEl array, holds colors of edge-neighbor elements + own color
  neighbor_nums::Array{Int32, 2}  # 4 by numEl array, holds element numbers of neighbors + own number, in same order as neighbor_colors
  color_cnt::Array{Int32, 1}  # number of elements in each color



 function PumiMesh2(dmg_name::AbstractString, smb_name::AbstractString, order, sbp::SBPOperator, opts; dofpernode=1, shape_type=1, coloring_distance=2)
  # construct pumi mesh by loading the files named
  # dmg_name = name of .dmg (geometry) file to load (use .null to load no file)
  # smb_name = name of .smb (mesh) file to load
  # order = order of shape functions
  # opts: dictionary of options
  # dofpernode = number of dof per node, default = 1
  # shape_type = type of shape functions, 0 = lagrange, 1 = SBP
  # coloring_distance : distance between elements of the same color, where distance is the minimum number of edges that connect the elements, default = 2

  println("\nConstructing PumiMesh2 Object")
  println("  sbp_name = ", smb_name)
  println("  dmg_name = ", dmg_name)
  mesh = new()
  mesh.numDofPerNode = dofpernode
  mesh.order = order
  mesh.numNodesPerElement = getNumNodes(order)
  mesh.shape_type = shape_type
  mesh.coloringDistance = coloring_distance
  num_Entities, mesh.m_ptr, mesh.mshape_ptr = init2(dmg_name, smb_name, order, shape_type=shape_type)
  mesh.f_ptr = createPackedField(mesh.m_ptr, "solution_field", dofpernode)

  # count the number of all the different mesh attributes
  mesh.numVert = convert(Int, num_Entities[1])
  mesh.numEdge =convert(Int,  num_Entities[2])
  mesh.numEl = convert(Int, num_Entities[3])
  mesh.numEntitiesPerType = [mesh.numVert, mesh.numEdge, mesh.numEl]
  mesh.numTypePerElement = [3, 3, 1]

  num_nodes_v = 1  # number of nodes on a vertex
  num_nodes_e = countNodesOn(mesh.mshape_ptr, 1) # on edge
  num_nodes_f = countNodesOn(mesh.mshape_ptr, 2) # on face
  num_nodes_entity = [num_nodes_v, num_nodes_e, num_nodes_f]
  mesh.numNodesPerType = num_nodes_entity
  mesh.typeOffsetsPerElement = zeros(Int, 4)
  pos = 1
  mesh.typeOffsetsPerElement[1] = pos
  for i=2:4
    pos += mesh.numTypePerElement[i-1]*mesh.numNodesPerType[i-1]
    mesh.typeOffsetsPerElement[i] = pos
  end
  println("mesh.typeOffsetsPerElement = ", mesh.typeOffsetsPerElement)



 
  # get pointers to mesh entity numberings
  mesh.vert_Nptr = getVertNumbering()
  mesh.edge_Nptr = getEdgeNumbering()
  mesh.el_Nptr = getFaceNumbering()
  mesh.entity_Nptrs = [mesh.vert_Nptr, mesh.edge_Nptr, mesh.el_Nptr]

  # create the coloring_Nptr
  el_mshape = getConstantShapePtr(2)
  mesh.coloring_Nptr = createNumberingJ(mesh.m_ptr, "coloring", el_mshape, 1)


  mesh.verts = Array(Ptr{Void}, mesh.numVert)
  mesh.edges = Array(Ptr{Void}, mesh.numEdge)
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

  for i=1:mesh.numEl
    mesh.elements[i] = getFace()
    incrementFaceIt()
  end

  mesh.numBC = opts["numBC"]

  # create array of all model edges that have a boundary condition
  bndry_edges_all = Array(Int, 0)
  for i=1:mesh.numBC
    key_i = string("BC", i)
    bndry_edges_all = [ bndry_edges_all; opts[key_i]]  # ugly but easy
  end

 mesh.numBoundaryEdges, num_ext_edges =  countBoundaryEdges(mesh, bndry_edges_all)

  # populate mesh.bndry_faces from options dictionary
#  mesh.bndry_faces = Array(Array{Int, 1}, mesh.numBC)
  mesh.bndry_offsets = Array(Int, mesh.numBC + 1)
  mesh.bndry_funcs = Array(BCType, mesh.numBC)
  boundary_nums = Array(Int, mesh.numBoundaryEdges, 2)

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

  # get array of all boundary mesh edges in the same order as in mesh.bndry_faces
#  boundary_nums = flattenArray(mesh.bndry_faces[i])
#  boundary_edge_faces = getEdgeFaces(mesh, mesh.bndry_faces)
  # use partially constructed mesh object to populate arrays

  mesh.elementNodeOffsets, mesh.typeNodeFlags = getEntityOrientations(mesh)
  numberDofs(mesh)
  getDofNumbers(mesh)  # store dof numbers

  # get sparsing information
  mesh.sparsity_bnds = zeros(Int32, 2, mesh.numDof)
  getSparsityBounds(mesh, mesh.sparsity_bnds)

  mesh.color_masks = Array(BitArray{1}, 4)  # one array for every color
  #colorMesh1(mesh, mesh.color_masks)
  numc = colorMesh2(mesh, mesh.color_masks)
  mesh.numColors = numc

  mesh.color_masks = Array(BitArray{1}, numc)  # one array for every color
  mesh.neighbor_colors = zeros(UInt8, 4, mesh.numEl)
  mesh.neighbor_nums = zeros(Int32, 4, mesh.numEl)
  getColors1(mesh, mesh.color_masks, mesh.neighbor_colors, mesh.neighbor_nums; verify=opts["verify_coloring"] )

  # get boundary information for entire mesh
  mesh.bndryfaces = Array(Boundary, mesh.numBoundaryEdges)
  getBoundaryArray(mesh, boundary_nums)

  # need to count the number of internal interfaces - do this during boundary edge counting
  mesh.numInterfaces = mesh.numEdge - num_ext_edges
  mesh.interfaces = Array(Interface, mesh.numInterfaces)
  getInterfaceArray(mesh)

  getCoordinates(mesh, sbp)  # store coordinates of all nodes into array

  mesh.dxidx = Array(T1, 2, 2, sbp.numnodes, mesh.numEl)
  mesh.jac = Array(T1, sbp.numnodes, mesh.numEl)
  mappingjacobian!(sbp, mesh.coords, mesh.dxidx, mesh.jac)

  # get face normals
  mesh.bndry_normals = Array(T1, 2, sbp.numfacenodes, mesh.numBoundaryEdges)
  getBoundaryFaceNormals(mesh, sbp, mesh.bndryfaces, mesh.bndry_normals)

  mesh.interface_normals = Array(T1, 2, 2, sbp.numfacenodes, mesh.numInterfaces)
  getInternalFaceNormals(mesh, sbp, mesh.interfaces, mesh.interface_normals)
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
    printcoords("coords.dat", mesh.coords)
  end

  if opts["write_sparsity"]
    rmfile("sparsity_bnds.dat")
    writedlm("sparsity_bnds.dat", mesh.sparsity_bnds.')
  end

  if opts["write_counts"]
    writeCounts(mesh)
  end

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
    nnodes = 3
  elseif order == 2
    nnodes = 7
  elseif order == 3
    nnodes = 12
  elseif order == 4
    nnodes = 18
  else
    println("Warning, unsupported order elements are requested")
  end

  return nnodes
end

function getEdgeFaces(mesh::PumiMesh2, bndry_edges::AbstractArray{Int, 1})
# get the faces corresponding to the boundary edges

bndry_faces = zeros(bndry_edges)
faces = Array(Ptr{Void}, 400)  # equivilent to apf::Up
for i=1:bndry_edges
  edgenum_i = bndry_edges[i]
  edge_i = mesh.edges[edgenum_i]
  numFace = countAdjacent(mesh.m_ptr, edge_i, 2)  # should be count upward
  getAdjacent(faces)
  facenum = getFaceNumber2(faces[1]) + 1

  bndry_faces[i] = facenum
end

return bndry_faces

end



  


function flattenArray{T}(A::AbstractArray{AbstractArray{T}, 1})
# copies array-of-array A into a flattened version B

  n = length(A)
  num_entries = 0
  for i=1:n  # count number of entries total
    num_entries += length(A[i])
  end

  B = zeros(T, num_entries)

  B_index = 1
  for i=1:n
    for j = 1:length(A[i])
      B[B_index] = A[i][j]
      B_index += 1
    end
  end

  return B

end



function getSparsityBounds(mesh::PumiMesh2, sparse_bnds::AbstractArray{Int32, 2})
# sparse_bnds : 2 by numDof array of the minimum, maximum dofs connected to 
# a dof in the Jacobian
# this works for high order elements despite not using mesh.elementNodeOffsets 
# because it does not reference what element a particular entity belongs to
# getDofBounds also does not need it because only the minimum and maximum
# are important, not the order of them

resetAllIts2()
# mesh iterator increment, retreval functions
iterators_inc = [incrementVertIt, incrementEdgeIt, incrementFaceIt]
iterators_get = [getVert, getEdge, getFace]
num_entities = [mesh.numVert, mesh.numEdge, mesh.numEl]
num_nodes_entity = mesh.numNodesPerType  # number of nodes on each type
                                         # of mesh entity

for etype=1:3  # loop over mesh entity types

  
  if (num_nodes_entity[etype] != 0)  # there are nodes here
    for entity = 1:num_entities[etype]  # loop over all entities of this type
      entity_ptr = iterators_get[etype]()  # get pointer to mesh entity

      for node = 1:num_nodes_entity[etype]  # loop over nodes on each entity
	# get the minimum, maximum related dof numbers
        min, max = getDofBounds(mesh, etype)
	# use the same min, max for all dofs on this node
	for dof = 1:mesh.numDofPerNode
	  dofnum = getNumberJ(mesh.dofnums_Nptr, entity_ptr, node - 1, dof-1)
          sparse_bnds[1, dofnum] = min
	  sparse_bnds[2, dofnum] = max

	end  # end loop over dofs
      end  # end loop over nodes on entity
      iterators_inc[etype]()
    end  # end loops over entities of this type
  end  # end if statement
end  # end loop over entity types

return nothing

end

# this could be generalized 3d?
function getDofBounds(mesh::PumiMesh2, etype::Integer) 
# gets the maximum, minimum dofs associated with the entity currently
# pointed to by the iterator specified by etype


iterators_get = [getVert, getEdge, getFace]
entity_ptr = iterators_get[etype]()

# get associated elements (distance-1 elements)
num_adj = countAdjacent(mesh.m_ptr, entity_ptr, 2)
el_arr = getAdjacent(num_adj)

# distance-1
#=
dofnums = zeros(Int32, mesh.numDofPerNode, mesh.numNodesPerElement, num_adj)

for i=1:num_adj
  el_i = getNumberJ(mesh.el_Nptr, el_arr[i], 0, 0) + 1
  sub_arr = sub(dofnums, :, :, i)
  getGlobalNodeNumbers(mesh, el_i, sub_arr)
end

min, max = getMinandMax(dofnums)

return min, max
=#

# get distance-2 elements
if mesh.coloringDistance >= 2
  edge_arr = Array(Ptr{Void}, num_adj*3)  # enough space for all edges, including repeats
  for i=1:num_adj  # get the edges
    sub_arr = sub(edge_arr, (3*(i-1) + 1):(3*i))
    getDownward(mesh.m_ptr, el_arr[i], 1, sub_arr)
  end

  # edge_arr now populated with all edges

  # print the edge numbers
#=
  for i=1:num_adj*3
    edge_num_i = getNumberJ(mesh.edge_Nptr, edge_arr[i], 0, 0)
    println("edge ", i, " has number ", edge_num_i)
  end
=#
  # now get the elements the edges belong to
  # count the number of elements first, then get them
  num_els = zeros(Int, length(edge_arr) + 1)  # count number of elements each edge has

#  println("length(num_els) = ", length(num_els))
#  println("length(edge_arr) = ", length(edge_arr))
  for i=1:length(edge_arr)
    num_els[i] = countAdjacent(mesh.m_ptr, edge_arr[i], 2)
  end

#  println("finished counting distance-2 elements")

  # now get them elements
  num_adj = sum(num_els)
  el_arr = Array(Ptr{Void}, num_adj)  # rebind el_arr to new array
  start_idx = 1
  end_idx = num_els[1]
  for i=1:length(edge_arr)
    edge_i = edge_arr[i]
    sub_arr = sub(el_arr, start_idx:end_idx)
    countAdjacent(mesh.m_ptr, edge_i, 2)
    getAdjacent(sub_arr)

    # update indices
    start_idx = end_idx + 1
    end_idx += num_els[i+1]
  end

  # now el_arr has all the elements, including repeats
  # num_adj = legnth(el_arr)

end  # end if distance-2



    
# get dofnums for the elements
dofnums = zeros(Int32, mesh.numDofPerNode, mesh.numNodesPerElement, num_adj)

for i=1:num_adj
  elnum_i = getNumberJ(mesh.el_Nptr, el_arr[i], 0, 0) + 1
  getGlobalNodeNumbers(mesh, elnum_i, sub(dofnums, :, :, i))
end

min, max = getMinandMax(dofnums)

return min, max

end


function getMinandMax{T}(arr::AbstractArray{T})
# need to check this function for type stability

min_entry = typemax(T)
max_entry = typemin(T)

@inbounds for i=1:length(arr)
  entry_i = arr[i]

  if entry_i < min_entry
    min_entry = entry_i
  end

  if entry_i > max_entry
    max_entry = entry_i
  end
end

return min_entry, max_entry

end

# perform distance-1 coloring of mesh 
function colorMesh1(mesh::PumiMesh2, masks::Array{BitArray{1}})
# each element must have a different color than its neighbors with which it 
# shares and edge

# figure out the number of colors
# this is a lot of random memory access
#=
num_neigh_max = 0
for i=1:mesh.numEl  # loop over elements
  el_i = mesh.elements[i]
  num_neigh_i = countBridgeAdjacent(mesh.m_ptr, element, 1, 2)

  if num_neigh_i > num_neigh_max
    num_neigh_max = num_neigh_i
  end
end
=#


# now perform the coloring
# visit each element, get colors of its neighbors
# give current element the lowest color possible
# also construct BitArray masks
#setNumberingOffset(mesh.coloring_Nptr, 1)  # set all values to -1 + 1 = 0

#=
# initialize masks to zero
for i=1:length(masks)
  masks[i] = falses(mesh.numEl)
  println("masks[i] = ", masks[i])
end
=#

for i=1:mesh.numEl
  numberJ(mesh.coloring_Nptr, mesh.elements[i], 0, 0, 0)
end

adj_size = 6  # guess number of neighboring faces
numc = 4  # guess number of colors
adj = Array(Ptr{Void}, adj_size)  # hold neighboring faces
adj_color =zeros(Int32, adj_size)  # colors of neighboring faces
cnt_colors = zeros(Int32, numc)  # count how many of each color

for i=1:mesh.numEl
  el_i = mesh.elements[i]
  # get faces that share a vert with el_i
  # this is actually more restrictive than we want because
  # edge stabilization only causes interaction between 
  # elements that share an edge
  # However, this is much easier to write an algorithm for

  # This is more restrictive than we want, but it ensure
  # that we void assigning the same color to two elements
  # bordering the same element
  # ie. we must avoid vertex neighbor non uniqueness,
  # which is not the same thing as having uniqueness because
  # there are unassigned values remaining
  # thus, this is a logical system where a double negative
  # is not a positive
 
  num_adj = countBridgeAdjacent(mesh.m_ptr, el_i, 0, 2)

  if num_adj > adj_size
    println("resizing adj")
    println("element number = ", i)
    resize!(adj, num_adj)
    resize!(adj_color, num_adj)
    adj_size = num_adj
  end


  # need to verify this works in parallel (proper ghosting)
  getBridgeAdjacent(adj)

  for j=1:num_adj
    adj_color[j] = getNumberJ(mesh.coloring_Nptr, adj[j], 0, 0)
  end

  min_color = getMinColor2(adj_color, numc)

  if min_color > numc
    resize!(cnt_colors, min_color)
    cnt_colors[min_color] = 0  # initialize new value to zero
    numc = min_color
  end

  numberJ(mesh.coloring_Nptr, el_i, 0, 0, min_color)

  cnt_colors[min_color] += 1  # update counts
#  masks[min_color][i] = true  # update mask

  fill!(adj_color, 0)

end
println("maximum number of neighbor faces = ", adj_size)
println("number of colors = ", numc)
println("number of each color = ", cnt_colors)
mesh.color_cnt = cnt_colors

return nothing

end



# perform distance-1 coloring of mesh 
function colorMesh2(mesh::PumiMesh2, masks::Array{BitArray{1}})
# each element must have a different color than its neighbors with which it 
# shares and edge

# figure out the number of colors
# this is a lot of random memory access
#=
num_neigh_max = 0
for i=1:mesh.numEl  # loop over elements
  el_i = mesh.elements[i]
  num_neigh_i = countBridgeAdjacent(mesh.m_ptr, element, 1, 2)

  if num_neigh_i > num_neigh_max
    num_neigh_max = num_neigh_i
  end
end
=#


# now perform the coloring
# visit each element, get colors of its neighbors
# give current element the lowest color possible
# also construct BitArray masks
#setNumberingOffset(mesh.coloring_Nptr, 1)  # set all values to -1 + 1 = 0

#=
# initialize masks to zero
for i=1:length(masks)
  masks[i] = falses(mesh.numEl)
  println("masks[i] = ", masks[i])
end
=#

for i=1:mesh.numEl
  numberJ(mesh.coloring_Nptr, mesh.elements[i], 0, 0, 0)
end

#adj_size = 6  # guess number of neighboring faces
numc = 4  # guess number of colors
#adj = Array(Ptr{Void}, adj_size)  # hold neighboring faces
#adj_color =zeros(Int32, adj_size)  # colors of neighboring faces

adj = Array(Ptr{Void}, 3)  # distance-1 edge neighbors
adj2 = Array(Array{Ptr{Void}, 1}, 3)  # distance-2 edge neighbors + distance-1 

for i=1:3
  adj2[i] = Array(Ptr{Void}, 4)
end

colors = zeros(Int32, 3*4)

cnt_colors = zeros(Int32, numc)  # count how many of each color
for i=1:mesh.numEl
#  println("Processing element ", i)
  el_i = mesh.elements[i]
  # get faces that share an edge with el_i, and all elements
  # those elements share an edge with
  # this is actually more restrictive than we want because
  # a unique distance-1 coloring, not distance-2
  # However, this is much easier to write an algorithm for

  # This is more restrictive than we want, but it ensure
  # that we avoid assigning the same color to two elements
  # bordering the same element
  # ie. we must avoid distance-2 neighbor non uniqueness,
  # which is not the same thing as having uniqueness because
  # there are unassigned values remaining
  # thus, this is a logical system where a double negative
  # is not a positive
 
  getDistance2Colors(mesh, i, adj, adj2, colors)

#  println("colors = \n", colors)
#  num_adj = countBridgeAdjacent(mesh.m_ptr, el_i, 0, 2)

#=
  if num_adj > adj_size
    println("resizing adj")
    println("element number = ", i)
    resize!(adj, num_adj)
    resize!(adj_color, num_adj)
    adj_size = num_adj
  end


  # need to verify this works in parallel (proper ghosting)
  getBridgeAdjacent(adj)

  for j=1:num_adj
    adj_color[j] = getNumberJ(mesh.coloring_Nptr, adj[j], 0, 0)
  end
=#
  min_color = getMinColor2(colors, numc)

  if min_color > numc
    resize!(cnt_colors, min_color)
    cnt_colors[min_color] = 0  # initialize new value to zero
    numc = min_color
  end

  numberJ(mesh.coloring_Nptr, el_i, 0, 0, min_color)

  cnt_colors[min_color] += 1  # update counts
#  masks[min_color][i] = true  # update mask

  fill!(colors, 0)

end
println("number of colors = ", numc)
println("number of each color = ", cnt_colors)
mesh.color_cnt = cnt_colors

return numc

end


function getDistance2Colors(mesh, elnum::Integer, adj, adj2, colors)
# get the distance-2 neighbors of a given element
# repeats included
# adj : array to be populated with distance-1 edge neighbors
# adj2 : array to be populated with sum of distance-1 and distance-2 edge
# neighbors
# colors : array to be populated with colors of elements in adj2

  el_i = mesh.elements[elnum]

  num_adj = countBridgeAdjacent(mesh.m_ptr, el_i, 1, 2)
#  adj = Array(Ptr{Void}, 3)
  getBridgeAdjacent(adj)

#  adj2 = Array(Array{Ptr{Void}, 1}, num_adj)  # hold distance 2 neighbors

  adj_cnt = zeros(Int, num_adj)
  for j=1:num_adj
    num_adj_j = countBridgeAdjacent(mesh.m_ptr, adj[j], 1, 2)
    adj_cnt[j] = num_adj_j + 1
 #   adj2[j] = Array(Ptr{Void}, num_adj_j + 1)
    getBridgeAdjacent(adj2[j])
    adj2[j][num_adj_j + 1] = adj[j]  # include the distance-1 neighbor
  end

  num_adj_total = sum(adj_cnt)  # total number of neighbors, including repeats

#  colors = zeros(Int32, num_adj_total)

  # get the colors of all the elements
  # also flatten adj2 array of arrays into a single array
  index = 1
  for i=1:num_adj
    for j=1:adj_cnt[i]
      colors[index] = getNumberJ(mesh.coloring_Nptr, adj2[i][j], 0, 0)
      index += 1
    end
  end


  return nothing
end




function getColors1{T, T2}(mesh, masks::AbstractArray{BitArray{1}, 1}, neighbor_colors::AbstractArray{T, 2}, neighbor_nums::AbstractArray{T2, 2}; verify=true )
# verify edge neighbor faces have different colors (ie. distance-1 coloring
# of graph where elements are the vertices connected by edges if 
# there is a mesh edge connecting them
# also construct the BitArray masks that describe which elements to perturb
# for which colors, and the list of neighbor + self colors
# masks is an array of bitarrays, number of arrays = number of colors
# neighbor_colors is 4 by numEl array of integers holding the colors
# of the nieghbors + own color
# neighbor_nums is 4 by numEl array of integers holding the element numbers
# of the neighbors + self

adj = Array(Ptr{Void}, 4)   # pointers to element i + 3 neighbors
adj_color = Array(Int32, 4)  # element colors

# initialize masks to 0 (false)
for i=1:length(masks)
  masks[i] = falses(mesh.numEl)
end

if verify
  println("verifying distance-1 coloring was sucessful")
end


cnt = 0
for i=1:mesh.numEl
  el_i = mesh.elements[i]

  # check edge neighbors only
 num_adj = countBridgeAdjacent(mesh.m_ptr, el_i, 1, 2)
  @assert num_adj <= 3
  getBridgeAdjacent(adj)

  adj[num_adj + 1] = el_i  # insert the current element

  # get color, element numbers
  for j=1:(num_adj + 1)
    adj_color[j] = getNumberJ(mesh.coloring_Nptr, adj[j], 0, 0)
    neighbor_colors[j, i] = adj_color[j]
    neighbor_nums[j, i] = getNumberJ(mesh.el_Nptr, adj[j], 0, 0) + 1
  end

  color_i = adj_color[num_adj + 1]  # color of current element
  masks[color_i][i] = true  # indicate element i gets perturbed by color_i

  if verify
    sort!(adj_color)

    # remove leading zeros
    nnz = countnz(adj_color)
    nz_arr = zeros(eltype(adj_color), nnz)
    start_idx = length(adj_color) - nnz + 1
    for i=1:nnz
      nz_arr[i] = adj_color[start_idx]
      start_idx += 1
    end


    if nz_arr != unique(nz_arr)
      println("element ", i, " has non unique colors")
      println("adj_color = ", nz_arr)
      cnt += 1
    end
  end   # end if verify
  fill!(adj_color, 0)
end

if verify
  println("color-1 verification finished")
  println("number of element with non unique coloring = ", cnt)
end

return cnt

end





function getMinColor{T}(adj::AbstractArray{T})
# adj contains colors of adjacent elements
  min_color = 1
  sort!(adj)  # adj must be in increasing order for this to work
  for i=1:length(adj)
    if min_color < adj[i]
      continue
    elseif min_color == adj[i]
	min_color = adj[i] + 1
    end  # if min_color > adj[i] do nothing
  end  # end for loop

  return min_color
end



function getMinColor2{T}(adj::AbstractArray{T}, numc::Integer)
# ensure uniqueness of neighboring colors
# adj is array of colors of adjacent faces
# numc is the current number of colors

mask = zeros(Bool, numc)

min_color = 0
for i=1:length(adj) # identify already used colors
  if adj[i] != 0
    mask[adj[i]] = true
  end
end

mask_sum = sum(mask)

if mask_sum == numc  # all existing colors used, so add another
  println("adding color ", numc + 1)
  min_color = numc + 1
elseif mask_sum == (numc - 1)  # there is exactly 1 color remaining
#  println("exactly 1 color remaining")
  # find out which color is missing and use it
  for i=1:numc
    if !mask[i]  # if mask is false
      min_color = i
    end
  end
else  # some colors are missing
#  println("getting minimum color")
  min_color = getMinColor(adj)  # get the minimum
end

@assert min_color != 0

return min_color

end



function getEntityOrientations(mesh::PumiMesh2)
# get the offset for node access of each mesh entity of each element
# to read/write from a Pumi field, accessing the data stored on an edge looks like:
# getComponents(f, e, abs(offset - node) - 1, vals)
# where f is the field, e is a mesh entity, node is the node index (one based),
# and vals is an array of values to be written to the field
# offset = 0 corresponds to the original ordering
# offset = number of nodes on the entity + 1 corresponds to reversing the 
#  ordering
# other values of offset should corresponds to rotations in 3d
# for 2d, the offsets for all nodes on an entity would typically be the same,
# corresponding to either flipping or not flipping the edge, but in 3d they
# might be different for a rotation
  # create array of arrays
  # hold an offset for every node of every entity of every element
  offsets = zeros(Uint8, mesh.numNodesPerElement, mesh.numEl)
  # hold a bit telling if the entity is in the proper orientation for this element
  # ie. it is "owned" by this element, so the nodes can be accessed in order
  flags = Array(BitArray{2}, 3)
  # create the arrays
  for i=1:3
    flags[i] = trues(mesh.numTypePerElement[i], mesh.numEl)
  end

  # now do edges

  edge_flags = flags[2]
  edges_i = Array(Ptr{Void}, 12)
  edgenode_range = mesh.typeOffsetsPerElement[2]:(mesh.typeOffsetsPerElement[3]-1)
  for i=1:mesh.numEl
    getDownward(mesh.m_ptr, mesh.elements[i], 1, edges_i)
 
    for j=1:mesh.numTypePerElement[2]  # loop over edges on element i
      # get orientation of edge j on element i
      # get global number of edge j
      edgenum_global = getNumberJ(mesh.edge_Nptr, edges_i[j], 0, 0) + 1

      orient = getEdgeOrientation(mesh, i, edgenum_global)

      # save value to flag array
      edge_flags[j, i] = div(orient + 1, 2)
      # write n = mesh.numNodesPerType[2] + 1 of orientation = -1, or n=0 if orientation=1
      if orient == 1
	offsets[edgenode_range, i] = 0
      else
	offsets[edgenode_range, i] = mesh.numNodesPerType[2] + 1
      end

    end  # end loop over edges
  end # end loop over elements

  return offsets, flags
end

function getEdgeOrientation(mesh::PumiMesh2, elnum::Integer, edgenum::Integer)
# figure out what the orientation of the specified edge is relative ot the element
# if the edge is in the same orientation as the element, return 1, otherwise -1

#  println("\nEntered getEdgeOrientation")
  # get all the vertices
  el_verts, tmp = getDownward(mesh.m_ptr, mesh.elements[elnum], 0)
  edge_verts, tmp = getDownward(mesh.m_ptr, mesh.edges[edgenum], 0)

  pos1 = 0  # position of edge_verts[1] in el_verts
  pos2 = 0  # position of edge_verts[2] in el_verts


  for i=1:3  # loop over el_verts
    if el_verts[i] == edge_verts[1]
      pos1 = i
    elseif el_verts[i] == edge_verts[2]
      pos2 = i
    end
  end

  @assert pos1 != 0
  @assert pos2 != 0

  if pos2 - pos1 > 0  # positively oriented
    return 1
  elseif pos2 - pos1 < 0  # negatively oriented
    return -1
  else
    println(STDERR, "Warning, bad orientation determination in PdePumiInterface getEdgeOrientation")
    return 0
  end

end  # end function



function getMeshEdgesFromModel{T}(mesh::PumiMesh2, medges::AbstractArray{Int, 1}, offset::Integer, boundary_nums::AbstractArray{T, 2})
# get the global numbers of the mesh edges lying on the model edges in medges
# offset is the index in boundary_nums to start with
# this allows populating the entire array without making temporary copies
# offset should start as 1, not zero

#  bndry_edges = zeros(Int, mesh.numBoundaryEdges)
  index = 0  # relative index in bndry_edges
  faces = Array(Ptr{Void}, 2)  # an edge can have maximum 2 faces using it
  for i=1:mesh.numEdge
    edge_i = mesh.edges[i]
    me_i = toModel(mesh.m_ptr, edge_i)
    me_dim = getModelType(mesh.m_ptr, me_i)
    me_tag = getModelTag(mesh.m_ptr, me_i)
    if me_dim == 1  # edge
      onBoundary = findfirst(medges, me_tag)

      if onBoundary != 0  # if mesh edge is on model edge
        
	# get face number
        numFace = countAdjacent(mesh.m_ptr, edge_i, 2)  # should be count upward

	@assert( numFace == 1)

        getAdjacent(faces)
        facenum = getFaceNumber2(faces[1]) + 1
        edgenum = getEdgeNumber2(edge_i)  # unneeded?

	boundary_nums[offset + index, 1] = facenum
        boundary_nums[offset + index, 2] = i

	index += 1
      end
    end
  end

  return offset + index

end  # end function


function numberDofs(mesh::PumiMesh2)
# assign dof numbers to entire mesh
# calculates numNodes, numDof
# assumes mesh elements have already been reordered
# this works for high order elements
  println("Entered numberDofs")

  # calculate number of nodes, dofs
  num_nodes_v = 1  # number of nodes on a vertex
  num_nodes_e = countNodesOn(mesh.mshape_ptr, 1) # on edge
  num_nodes_f = countNodesOn(mesh.mshape_ptr, 2) # on face
  println("num_nodes_v = ", num_nodes_v)
  println("num_nodes_e = ", num_nodes_e)
  println("num_nodes_f = ", num_nodes_f)
  numnodes = num_nodes_v*mesh.numVert + num_nodes_e*mesh.numEdge + num_nodes_f*mesh.numEl
  numDof = mesh.numDofPerNode*numnodes
  println("expected number of dofs = ", numDof)


  if (numDof > 2^30)
    println("Warning: too many dofs, renumbering will fail")
  end

  # save numbers to mesh
  mesh.numNodes = numnodes
  mesh.numDof = numDof

  # initally number all dofs as numDof+1 to 2*numDof
  # this allows quick check to see if somthing is labelled or not

  resetAllIts2()
  # mesh iterator increment, retreval functions
  iterators_inc = [incrementVertIt, incrementEdgeIt, incrementFaceIt]
  iterators_get = [getVert, getEdge, getFace]
  num_entities = [mesh.numVert, mesh.numEdge, mesh.numEl]
  num_nodes_entity = [num_nodes_v, num_nodes_e, num_nodes_f]

#  mesh.numNodesPerType = num_nodes_entity

  println("num_entities = ", num_entities)
  println("num_nodes_entity = ", num_nodes_entity)

#  curr_dof = 1
  curr_dof = numDof+1
  for etype = 1:3 # loop over entity types
#    println("etype = ", etype)
    if (num_nodes_entity[etype] != 0)  # if no nodes on this type of entity, skip
      for entity = 1:num_entities[etype]  # loop over all entities of this type
#	println("  entity number: ", entity)
	entity_ptr = iterators_get[etype]()  # get entity

	for node = 1:num_nodes_entity[etype]
#	  println("    node : ", node)
	  for dof = 1:mesh.numDofPerNode
	    numberJ(mesh.dofnums_Nptr, entity_ptr, node-1, dof-1, curr_dof)
#	    println("      entity ", entity_ptr, " labelled ", curr_dof)
	    curr_dof += 1
	  end  # end loop over dof
	end  # end loop over node

      iterators_inc[etype]()
      end  # end loop over entitiesa
    end  # end if 
  end  # end loop over entity types

  println("Finished initial numbering of dofs") 


  println("Performing final numbering of dofs")

  verts_i = Array(Ptr{Void}, 12)
  edges_i = Array(Ptr{Void}, 12)
# TODO: move all if statements out one for loop (check only first dof on each node)
  curr_dof = 1
  for i=1:mesh.numEl
#    println("element number: ", i)
    # get vertices, edges for this element
    numVert = getDownward(mesh.m_ptr, mesh.elements[i], 0, verts_i)
#    println("verts_i = ", verts_i)
#    println("mesh.verts = ", mesh.verts)
    numEdge = getDownward(mesh.m_ptr, mesh.elements[i], 1, edges_i)
    el_i = mesh.elements[i]
    for j=1:3  # loop over vertices, edges
#      println("  vertex and edge number: ", j)
      vert_j = verts_i[j]
      edge_j = edges_i[j]
#      println("  vert_j = ", vert_j)
#      println("  edge_j = ", edge_j)
      for k=1:mesh.numDofPerNode # loop over vertex dofs
#	println("    dof number: ", k)
#	println("    vert_j = ", vert_j)
        dofnum_k = getNumberJ(mesh.dofnums_Nptr, vert_j, 0, k-1)
	if dofnum_k > numDof  # still has initial number
	  # give it new (final) number
	  numberJ(mesh.dofnums_Nptr, vert_j, 0, k-1, curr_dof)
	  curr_dof += 1
	end
      end  # end loop over vertex dofs
      
      # loop over nodes on edge
      for k=1:num_nodes_entity[2]  # loop over nodes
	for p=1:mesh.numDofPerNode  # loop over dofs
	  dofnum_p = getNumberJ(mesh.dofnums_Nptr, edge_j, k-1, p-1)
	  if dofnum_p > numDof  # still has initial number
	    # give it new (final) number
	    numberJ(mesh.dofnums_Nptr, edge_j, k-1, p-1, curr_dof)
	    curr_dof += 1
	  end
	end
      end
    end  # end loop over vertices, edges
    # label face nodes
    for k=1:num_nodes_entity[3]  # loop over nodes on face
      for p=1:mesh.numDofPerNode  # loop over dofs
	dofnum_p = getNumberJ(mesh.dofnums_Nptr, el_i, k-1, p-1)
	if dofnum_p > numDof
	  numberJ(mesh.dofnums_Nptr, el_i, k-1, p-1, curr_dof)
	  curr_dof += 1
	end
      end
    end  # end loop over face nodes
  end  # end loop over elements




  println("finished performing final dof numbering")

  println("number of dofs = ", curr_dof - 1)
  if (curr_dof -1) != numDof 
    println("Warning: number of dofs assigned is not equal to teh expected number")
    println("number of dofs assigned = ", curr_dof-1, " ,expected number = ", numDof)
  else
    println("Dof numbering is sane")
  end



  resetAllIts2()
  return nothing

end

#=
function numberDofs(mesh::PumiMesh2)
# number the degrees of freedom of the mesh, using the apf::Numbering* stroed
# in the mesh

  # move this into a function
  resetAllIts2()
  println("performing initial numbering of dofs")
  # calculate number of nodes, dofs (works for first and second order)
  numnodes = mesh.order*mesh.numVert 
  numdof = numnodes*mesh.numDofPerNode
  # number dofs
  ctr= 1
  for i=1:mesh.numVert
    for j=1:mesh.numDofPerNode
      numberJ(mesh.dofnums_Nptr, mesh.verts[i], 0, j-1, ctr)
      println("vertex ", i,  " numbered ", ctr)
      ctr += 1
    end
  end

  if mesh.order >= 2
    for i=1:mesh.numEdges
      for j=1:mesh.numDofPerNode
        numberJ(mesh.dofnums_Nptr, mesh.edges[i], 0, j-1, ctr)
        ctr += 1
      end
    end
  end

  mesh.numNodes = numnodes
  mesh.numDof = numdof

return nothing

end  # end function
=#


function countBoundaryEdges(mesh::PumiMesh2, bndry_edges_all)
  # count boundary edges by checking if their model edge has a BC
  # count number of external edges by checking the number of upward adjacencies
  # store array of [element number, global edge number]
  resetEdgeIt()
  bnd_edges_cnt = 0
  external_edges_cnt = 0
  bnd_edges = Array(Int, mesh.numEdge, 2)
  faces = Array(Ptr{Void}, 2)  # edge has maximum 2 faces
  for i=1:mesh.numEdge
    edge_i = getEdge()

    # get  model edge info
    me_i = toModel(mesh.m_ptr, edge_i)
    me_dim = getModelType(mesh.m_ptr, me_i)
    me_tag = getModelTag(mesh.m_ptr, me_i)

    # get mesh face info
    numFace = countAdjacent(mesh.m_ptr, edge_i, 2)  # should be count upward
    if numFace == 1  # external edges
      external_edges_cnt += 1
    end

    if me_dim == 1  # if not classified on model edge
      index = findfirst(bndry_edges_all, me_tag)



      if index != 0  # if model edge has a BC on i

	getAdjacent(faces)
	facenum = getFaceNumber2(faces[1]) + 1



	bnd_edges_cnt += 1
	bnd_edges[bnd_edges_cnt, 1] = facenum
	bnd_edges[bnd_edges_cnt, 2] = i
      end
    end  # end if me_dim == 1
    incrementEdgeIt()

  end  # end for loop


#  mesh.boundary_nums = bnd_edges[1:bnd_edges_cnt, :] # copy, bad but unavoidable
#  mesh.numBoundaryEdges = bnd_edges_cnt

return bnd_edges_cnt, external_edges_cnt

end  # end function




function getDofNumbers(mesh::PumiMesh2)
# populate array of dof numbers, in same shape as solution array u (or q)

println("in getDofNumbers")
println("numNodesPerElement = ", mesh.numNodesPerElement)

mesh.dofs = Array(Int32, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)

for i=1:mesh.numEl
  dofnums = getGlobalNodeNumbers(mesh, i)

  for j=1:mesh.numNodesPerElement
    for k=1:mesh.numDofPerNode  # loop over dofs on the node
      mesh.dofs[k, j, i] = dofnums[k,j]
    end
  end
end

#println("mesh.dof = ", mesh.dofs)

return nothing

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
  tmp, num_Entities, m_ptr, mshape_ptr = init2(dmg_name, smb_name, order, load_mesh=false, shape_type=mesh.shape_type) # do not load new mesh
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


function getCoordinates(mesh::PumiMesh2, sbp::SBPOperator)
# populate the coords array of the mesh object

mesh.coords = Array(Float64, 2, sbp.numnodes, mesh.numEl)

#println("entered getCoordinates")

coords_i = zeros(3,3)
coords_it = zeros(3,2)
for i=1:mesh.numEl  # loop over elements
  
  el_i = mesh.elements[i]
  (sizex, sizey) = size(coords_i)
  getFaceCoords(el_i, coords_i, sizex, sizey)  # populate coords

#  println("coords_i = ", coords_i)

  coords_it[:,:] = coords_i[1:2, :].'
#  println("coords_it = ", coords_it)
  mesh.coords[:, :, i] = calcnodes(sbp, coords_it)
#  println("mesh.coords[:,:,i] = ", mesh.coords[:,:,i])
end

return nothing

end




function getElementVertCoords(mesh::PumiMesh2, elnum::Integer, coords::AbstractArray{Float64,2})
# get the coordinates of the vertices of an element
# elnum is the number of the element
# coords is the array to be populated with the coordinates
# each column of coords contains the coordinates for a vertex
# coords must be 3x3

  el_j = mesh.elements[elnum]
  (sizex, sizey) = size(coords)
  getFaceCoords(el_j, coords, sizex, sizey)  # populate coords

  return nothing

end # end function

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

#=
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
=#


function getGlobalNodeNumbers(mesh::PumiMesh2, elnum::Integer)

  dofnums = zeros(Int32, mesh.numDofPerNode, mesh.numNodesPerElement)
  getGlobalNodeNumbers(mesh, elnum, dofnums)

  return dofnums
end

function getGlobalNodeNumbers(mesh::PumiMesh2, elnum::Integer, dofnums::AbstractArray{Int32, 2})
# gets global node numbers of all dof on all nodes of the element
# output formap is array [numdofpernode, nnodes]  (each column contains dof numbers for a node)
# dofnums must be passable to C
# 2D only
el_i = mesh.elements[elnum]
type_i = getType(mesh.m_ptr, el_i)

# calculate total number of nodes
#nnodes = 3 + 3*mesh.numNodesPerType[2] + mesh.numNodesPerType[3]
nnodes = mesh.numNodesPerElement
numdof = nnodes*mesh.numDofPerNode

node_entities = getNodeEntities(mesh.m_ptr, mesh.mshape_ptr, el_i)
node_offsets = view(mesh.elementNodeOffsets[:, elnum])
#println("node_offsets = ", [Int(i) for i in node_offsets])
#println("size(node_offsets) = ", size(node_offsets))
#println("node_entities = ", node_entities)
#println("size(node_entities) = ", size(node_entities))
#dofnums = zeros(Int,mesh.numDofPerNode, nnodes)  # to be populated with global dof numbers
#println("nnodes = ", nnodes)
#println("dofpernode = ", mesh.numDofPerNode)



num_entities = [3, 3, 1] # number of vertices, edges, faces

#=
col = 1 # current node
for i=1:3  # loop over verts, edges, faces
  for j=1:num_entities[i]  # loop over all entities of this type
    for k=1:mesh.numNodesPerType[i]  # loop over nodes on this entity
      entity = node_entities[col]  # current entity
      for p=1:mesh.numDofPerNode  # loop over all dofs
	dofnums[p, col] = getNumberJ(mesh.dofnums_Nptr, entity, k-1, p-1)
      end
      col += 1
    end  # end loop over nodes on curren entity
  end  # end loop over entities of current type
end  # end loop over entity types
=#

PumiInterface.getDofNumbers(mesh.dofnums_Nptr, node_entities, node_offsets, el_i, dofnums)  # C implimentation


#=

for i=1:nnodes
  for j=1:mesh.numDofPerNode
    dofnums[j, i] = getNumberJ(mesh.dofnums_Nptr, node_entities[i], 0, j-1)
  end
end
=#

return nothing


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

#=
function getBoundaryEdgeNums(mesh::PumiMesh2)
# get vector of edge numbers that are on boundary

edge_nums = mesh.boundary_nums
return edge_nums
end
=#


# this doesn't exist for a 2D mesh
function getBoundaryFaceNums(mesh::PumiMesh2)
# get array of face numbers that are on boundary
# this only works in 2d
# this function doesn't work yet
#face_nums = zeros(Int, numEdges_on_boundary)
face_nums = zeros(Int, 0)
return face_nums
end


function getAdjacentEntityNums(mesh::PumiMesh, entity_index::Integer, input_dimension::Integer, output_dimension::Integer)
# gets the numbers of the adjacent entities (upward or downward) of the given entity
# entity_index specifies the index of the input entity
# input_dimension specifies the dimension of the input entity (0 =vert, 1 = edge ...)
# output_dimension is the dimension of the adjacent  entities to be retrieved
# returns an array containing the indicies of the adjacent entities, and the number of them
# the length of the array might be greater than the number of adjacencies, so use the second returned value for the number of adjacencies


if (input_dimension == output_dimension)
  println("cannot get same level adjacencies")
end


if typeof(mesh) <: PumiMesh2

  array1 = [mesh.vert_Nptr, mesh.edge_Nptr, mesh.el_Nptr]
  #array2 = [mesh.verts; mesh.edges; mesh.elements]  # array of arrays
  array2 = Array(Array{Ptr{Void},1}, 3)
  array2[1] = mesh.verts
  array2[2] = mesh.edges
  array2[3] = mesh.elements
else # 3d mesh
   array1 = [mesh.vert_Nptr, mesh.edge_Nptr, mesh.face_Nptr,  mesh.el_Nptr]
  #array2 = [mesh.verts; mesh.edges; mesh.elements]  # array of arrays
  array2 = Array(Array{Ptr{Void},1}, 4)
  array2[1] = mesh.verts
  array2[2] = mesh.edges
  array2[3] = mesh.faces
  array2[4] = mesh.elements
end

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
  down_edges = Array(Ptr{Void}, 3)
  numedges = getDownward(mesh.m_ptr, face, 1, down_edges)
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

  down_edges = Array(Ptr{Void}, 3)
   numedges = getDownward(mesh.m_ptr, element, 1, down_edges)

  edgenum_local = 0
   for j = 1:3  # find which local edge is edge_i
    if down_edges[j] == edge
      edgenum_local = j
    end
  end

  return edgenum_local

end

function getBoundaryArray(mesh::PumiMesh2, boundary_nums::AbstractArray{Int, 2})
# get an array of type Boundary for SBP
# creating an an array of a user defined type seems like a waste of memory operations
# bnd_array is a vector of type Boundary with length equal to the number of edges on the boundary of the mesh

#  mesh.bndryfaces = Array(Boundary, mesh.numBoundaryEdges)

  for i=1:mesh.numBoundaryEdges
    facenum = boundary_nums[i, 1]
    edgenum_global = boundary_nums[i, 2]
#    println("edgenum_global = ", edgenum_global)
    edgenum_local = getBoundaryEdgeLocalNum(mesh, edgenum_global)
    mesh.bndryfaces[i] = Boundary(facenum, edgenum_local)
  end

  return nothing
end


function getInterfaceArray(mesh::PumiMesh2)
# get array of [elementL, elementR, edgeL, edgeR] for each internal edge,
# where elementL and R are the elements that use the edge, edgeL R are the
# local number of the edge within the element
# interfaces is the array to be populated with the data, it must be 
# number of internal edges by 1
# the nodemap is used to determine the node ordering of elementR relative to
# elementL
# the edge is "owned" by elementL, ie the first node on the edge is the first
# node from the perspective of elementL

  # only need internal boundaries (not external)
#  num_ext_edges = size(getBoundaryEdgeNums(mesh))[1]  # bad memory efficiency
#  num_int_edges = getNumEdges(mesh) - num_ext_edges

#  new_bndry = Boundary(2, 3)
#   println("new_bndry = ", new_bndry)

#  new_interface = Interface(1, 2, 3, 4)
#   println("new_interface = ", new_interface)

#   println("num_int_edges = ", num_int_edges)

#  interfaces = Array(typeof(new_interface), num_int_edges)
  # get number of nodes affecting an edge
  num_edge_nodes = countAllNodes(mesh.mshape_ptr, 1)

  nodemap = Array(num_edge_nodes:(-1):1)

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

      if abs(centroid1[1] - centroid2[1]) > 1e-10  # if big enough difference
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

      mesh.interfaces[pos] = Interface(elementL, elementR, edgeL, edgeR, nodemap)
#       println("updating pos")
      pos += 1

  #    print("\n")

    end  # end if internal edge

#     print("\n")
  end  # end loop over edges

  return nothing

end  # end function



function getBoundaryFaceNormals{Tmsh}(mesh::PumiMesh2, sbp::SBPOperator, bndry_faces::AbstractArray{Boundary, 1}, face_normals::Array{Tmsh, 3})

  nfaces = length(bndry_faces)

  alpha = zeros(Tmsh, 2,2)
  dxidx = mesh.dxidx
  jac = mesh.jac
  for i=1:nfaces
    element_i = bndry_faces[i].element
    face_i = bndry_faces[i].face
    for j=1:sbp.numfacenodes
      node_index = sbp.facenodes[j, face_i]

      # calculate the mysterious parameter alpha
      for di1=1:2
	for di2=1:2
	  alpha[di1, di2] = (dxidx[di1,1, node_index, element_i].*dxidx[di2,1, node_index, element_i] + dxidx[di1,2, node_index, element_i].*dxidx[di2,2, node_index, element_i])*jac[node_index,element_i]
	end
      end

    # call SBP function
    smallmatvec!(alpha, view(sbp.facenormal, :, face_i), view(face_normals, :, j, i))

    end  # end loop over face nodes
  end  # end loop over faces

  return nothing
end


function getInternalFaceNormals{Tmsh}(mesh::PumiMesh2, sbp::SBPOperator, internal_faces::AbstractArray{Interface, 1}, face_normals::Array{Tmsh, 4})

  nfaces = length(internal_faces)

  alpha = zeros(Tmsh, 2,2)
  dxidx = mesh.dxidx
  jac = mesh.jac
  for i=1:nfaces
    element_iL = internal_faces[i].elementL
    face_iL = internal_faces[i].faceL
    element_iR = internal_faces[i].elementR
    face_iR = internal_faces[i].faceR
    for j=1:sbp.numfacenodes

      # calculate left face normal
      node_index = sbp.facenodes[j, face_iL]

      # calculate the mysterious parameter alpha
      for di1=1:2
	for di2=1:2
	  alpha[di1, di2] = (dxidx[di1,1, node_index, element_iL].*dxidx[di2,1, node_index, element_iL] + dxidx[di1,2, node_index, element_iL].*dxidx[di2,2, node_index, element_iL])*jac[node_index,element_iL]
	end
      end

    # call SBP function
    smallmatvec!(alpha, view(sbp.facenormal, :, face_iL), view(face_normals, :, 1, j, i))
      # calculate right fae normal
      node_index = sbp.facenodes[j, face_iR]

      # calculate the mysterious parameter alpha
      for di1=1:2
	for di2=1:2
	  alpha[di1, di2] = (dxidx[di1,1, node_index, element_iR].*dxidx[di2,1, node_index, element_iR] + dxidx[di1,2, node_index, element_iR].*dxidx[di2,2, node_index, element_iR])*jac[node_index,element_iR]
	end
      end

    # call SBP function
    smallmatvec!(alpha, view(sbp.facenormal, :, face_iR), view(face_normals, :, 2, j, i))




    end  # end loop over face nodes
  end  # end loop over faces

  return nothing
end







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
# it uses mesh.elementNodeOffsets to access the pumi field values in the 
# right order of the given element
# not that this performs a loop over elements, so some values get
# written several times.  This is why is is necessary to write to the
# field in the right order
 # dofnums = zeros(Int, mesh.numDofPerNode)
#  u_vals = zeros(mesh.numDofPerNode)


  num_entities = [3, 3, 1] # number of vertices, edges, faces

  q_vals = zeros(4)

  for el=1:mesh.numEl
    el_i = mesh.elements[el]
    node_entities = getNodeEntities(mesh.m_ptr, mesh.mshape_ptr, el_i)
    col = 1 # current node of the element
    for i=1:3  # loop over verts, edges, faces
      for j=1:num_entities[i]  # loop over all entities of this type
	for k=1:mesh.numNodesPerType[i]  # loop over nodes on this entity
	  entity = node_entities[col]  # current entity
	  offset_k = mesh.elementNodeOffsets[col, el] # offset for current node

          # get solution values
	  for p=1:mesh.numDofPerNode  # loop over all dofs
	    new_node = abs(offset_k - k) - 1
	    dofnum_p = getNumberJ(mesh.dofnums_Nptr, entity, new_node, p-1)
	    q_vals[p] = u[dofnum_p]
	  end
          # save to mesh
          setComponents(mesh.f_ptr, entity, abs(offset_k - k) - 1, q_vals)

	  col += 1
	end  # end loop over nodes on curren entity
      end  # end loop over entities of current type
    end  # end loop over entity types
  end  # end loop over elements


  return nothing
end  # end function saveSolutionToMesh


function retrieveSolutionFromMesh(mesh::PumiMesh2, u::AbstractVector)
# retrieve solution from mesh (after mesh adaptation)
# mesh needs to have been reinitilized after mesh adaptation, u needs to be the right size for the new mesh




  num_entities = [3, 3, 1] # number of vertices, edges, faces

  q_vals = zeros(4)

  for el=1:mesh.numEl
    el_i = mesh.elements[el]
    node_entities = getNodeEntities(mesh.m_ptr, mesh.mshape_ptr, el_i)
    col = 1 # current node of the element
    for i=1:3  # loop over verts, edges, faces
      for j=1:num_entities[i]  # loop over all entities of this type
	for k=1:mesh.numNodesPerType[i]  # loop over nodes on this entity
	  entity = node_entities[col]  # current entity
	  offset_k = mesh.elementNodeOffsets[col, el]
         
	  # get values from mesh
          getComponents(mesh.f_ptr, entity, abs(offset_k - k) - 1, q_vals)

          # put solution values into vector u
	  for p=1:mesh.numDofPerNode  # loop over all dofs
	    dofnum_p = getNumberJ(mesh.dofnums_Nptr, entity, abs(offset_k - k) -1, p-1)
	    u[dofnum_p] = q_vals[p]
	  end
	  col += 1
	end  # end loop over nodes on curren entity
      end  # end loop over entities of current type
    end  # end loop over entity types
  end  # end loop over elements

  

  return nothing
end  # end function saveSolutionToMesh


function retrieveNodeSolution(f_ptr, entity, u_node::AbstractVector)
# retrieve solution on a particular entity, stores it in u_node
# u_node must be a vector of length mesh.numDofPerNode
# used during mesh adaptation
# this is a low level function because t takes in f_ptr, entity rather than the mesh object and an index

  getComponents(f_ptr, entity, 0, u_node)

end


function printBoundaryEdgeNums(mesh::PumiMesh)

  n = mesh.numBC

  bndry = 1
  edges = Array(Ptr{Void}, 12)
  for i=1:n
    fname = string("boundary_edge_verts", i, ".dat")
    println("printing ", fname)
    f = open(fname, "a+")

    start_index = mesh.bndry_offsets[i]
    end_index = mesh.bndry_offsets[i+1] - 1
    num_edge = end_index - start_index + 1
    arr = Array(Ptr{Void}, num_edge)


    for i=1:num_edge  # get the mesh edge pointers
      el = mesh.bndryfaces[bndry].element
      local_face = mesh.bndryfaces[bndry].face
      el_ptr = mesh.elements[el]
      getDownward(mesh.m_ptr, el_ptr, 1, edges)
      arr[i] = edges[local_face]
      bndry += 1
    end
    
    printEdgeVertNumbers(arr, mesh.edge_Nptr, mesh.vert_Nptr; fstream=f)
    close(f)
  end

  return nothing
end  # end function


function printdxidx(name::AbstractString, mat)
# print the 4d matrix formatted like dxidx

(p, p, m, n) = size(mat)

f = open(name, "a+")

for i=1:n  # loop over elements
  for j=1:m  # loop over nodes on an element
    println(f, "el ", i, " node ", j, " xi_vec = ", mat[1, :, j, i], " eta_vec = ", mat[2, :, j, i])
  end
end

close(f)

return nothing
end


function printcoords(name::AbstractString, coords)

  f = open(name, "a+")

  (d, numnodes, numel) = size(coords)

  for i=1:numel  # loop over elements
    for j=1:numnodes # loop over nodes
      print(f, i, " ", j)

      for k=1:d
	print(f, " ",  coords[k, j, i])
      end

      print(f, "\n")
    end
  end

  return nothing
end



function getdiffelementarea{T, T2, T3}(nrm::AbstractArray{T,1}, dxidx::AbstractArray{T2,2},
                               workvec::AbstractArray{T3,1})
  fill!(workvec, zero(T3))
  for di1 = 1:size(nrm,1)
    for di2 = 1:size(nrm,1)
      workvec[di2] += nrm[di1]*dxidx[di1,di2]
    end
  end
  return norm(workvec)
end


function writeVisFiles(mesh::PumiMesh, fname::AbstractString)
  # writes vtk files 

  writeVtkFiles(fname, mesh.m_ptr)

  return nothing
end



function writeCounts(mesh::PumiMesh2)
# write values needed for memory usage estimate
vals = Array(Int, 9)
vals[1] = mesh.numVert
vals[2] = mesh.numEdge
vals[3] = mesh.numEl
vals[4] = mesh.numBoundaryEdges
vals[5] = mesh.numNodesPerType[1]
vals[6] = mesh.numNodesPerType[2]
vals[7] = mesh.numNodesPerType[3]
vals[8] = mesh.numDofPerNode


# estimate jacobian storage size
acc = 0
for i=1:mesh.numDof
  acc += mesh.sparsity_bnds[2, i] - mesh.sparsity_bnds[1, i]
end

size_nz = 64*acc
size_rowval = 32*acc
size_colptr = 32*mesh.numDof

vals[9] = size_nz + size_rowval + size_colptr

writedlm("counts.txt", vals)
			  
return nothing
end



end  # end of module



