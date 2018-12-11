export PumiMesh3

mutable struct PumiMesh3{T1} <: PumiMesh3CG{T1}   # 3d pumi mesh, tetrahedron only
  m_ptr::Ptr{Void}  # pointer to mesh
  mshape_ptr::Ptr{Void} # pointer to mesh's FieldShape
  f_ptr::Ptr{Void} # pointer to apf::field for storing solution during mesh adaptation

  shape_type::Int  # type of shape functions

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
  numNodesPerType::Array{Int, 1}  # number of nodes classified on each vertex, edge, face, length 4
  numEntities::Array{Int, 1}  # number of entities of each type in the mesh, length 4
  numEntitiesPerElement::Array{Int, 1}  # number of vertices, edges, and regions in each element

  # parallel info
  comm::MPI.Comm  # MPI Communicator
  myrank::Int  # MPI rank, zero based
  commsize::Int # MPI comm size

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



 function PumiMesh3{T1}(dmg_name::AbstractString, smb_name::AbstractString, order, sbp::AbstractSBP; dofpernode=1, shape_type=0)  where {T1}
  # construct pumi mesh by loading the files named
  # dmg_name = name of .dmg (geometry) file to load (use .null to load no file)
  # smb_name = name of .smb (mesh) file to load
  # order = order of shape functions
  # dofpernode = number of dof per node, default = 1
  # shape_type = type of shape functions to use, 0 = lagrange, 1 = SBP

  mesh = new{T1}()
  mesh.comm = MPI.COMM_WORLD
  mesh.commsize = MPI.Comm_size(MPI.COMM_WORLD)
  mesh.myrank = MPI.Comm_rank(MPI.COMM_WORLD)

  @assert mesh.commsize == 1  # not paralllezed yet

  mesh.numDofPerNode = dofpernode
  mesh.order = order
  mesh.numNodesPerElement = getNumNodes(order)
  mesh.shape_type = shape_type
  @time mesh.numEntities, mesh.m_ptr, mesh.mshape_ptr, n_arr = apf.init(dmg_name, smb_name, order, shape_type=shape_type)
  mesh.f_ptr = apf.findField(mesh.m_ptr, "solution_field")
  if mesh.f_ptr == C_NULL
    @time mesh.f_ptr = apf.createPackedField(mesh.m_ptr, "solution_field", dofpernode)
  end

  mesh.numVert = convert(Int, mesh.numEntities[1])
  mesh.numEdge = convert(Int,  mesh.numEntities[2])
  mesh.numFace = convert(Int, mesh.numEntities[3])
  mesh.numEl = convert(Int, mesh.numEntities[4])

  # get number of nodes on each type
  mesh.numNodesPerType = zeros(Int, 4)
  mesh.numNodesPerType[1] = apf.countNodesOn(mesh.mshape_ptr, 0)  # number of nodes on a vertex
  mesh.numNodesPerType[2]= apf.countNodesOn(mesh.mshape_ptr, 1) # on edge
  mesh.numNodesPerType[3] = apf.countNodesOn(mesh.mshape_ptr, 2) # on face
  mesh.numNodesPerType[4]= apf.countNodesOn(mesh.mshape_ptr, 4) # on element interior 
  mesh.numEntitiesPerElement = [4, 6, 4, 1]

  mesh.numNodesPerElement = 0
  for i=1:4
    mesh.numNodesPerElement += mesh.numNodesPerType[i]*mesh.numEntitiesPerElement[i]
  end


  # get pointers to mesh entity numberings
  mesh.vert_Nptr = n_arr[1] #getVertNumbering()
  mesh.edge_Nptr = n_arr[2] #getEdgeNumbering()
  mesh.face_Nptr = n_arr[3] #getFaceNumbering()
  mesh.el_Nptr = n_arr[4] #getElNumbering()

   mesh.verts = Array{Ptr{Void}}(mesh.numVert)
  mesh.edges = Array{Ptr{Void}}(mesh.numEdge)
  mesh.faces = Array{Ptr{Void}}(mesh.numFace)
  mesh.elements = Array{Ptr{Void}}(mesh.numEl)
  mesh.dofnums_Nptr = apf.createNumberingJ(mesh.m_ptr, "reordered dof numbers", mesh.mshape_ptr, dofpernode)  # 1 dof per node




  # get pointers to all MeshEntities
  # also initilize the field to zero
  resetAllIts2(mesh.m_ptr)
#  comps = zeros(dofpernode)
  it = apf.MeshIterator(mesh.m_ptr, 0)
  for i=1:mesh.numVert
    mesh.verts[i] = apf.iterate(mesh.m_ptr, it)
  end
  apf.free(mesh.m_ptr, it)

  it = apf.MeshIterator(mesh.m_ptr, 1)
  for i=1:mesh.numEdge
    mesh.edges[i] = apf.iterate(mesh.m_ptr, it)
  end
  apf.free(mesh.m_ptr, it)

  it = apf.MeshIterator(mesh.m_ptr, 2)
  for i=1:mesh.numFace
    mesh.faces[i] = apf.iterate(mesh.m_ptr, it)
  end
  apf.free(mesh.m_ptr, it)

  it = apf.MeshIterator(mesh.m_ptr, 3)
  for i=1:mesh.numEl
    mesh.elements[i] = apf.iterate(mesh.m_ptr, it)
  end
  apf.free(mesh.m_ptr, it)

  # use partially constructed mesh object to populate arrays
  checkConnectivity(mesh)

#  numberDofs(mesh)

  countBoundaryFaces(mesh)
  mesh.bndryfaces = Array{Boundary}(mesh.numBoundaryFaces)
  getBoundaryArray(mesh)

  numberDofs(mesh)
  getDofNumbers(mesh)

  countBoundaryFaces(mesh)
  mesh.bndryfaces = Array{Boundary}(mesh.numBoundaryFaces)
  getBoundaryArray(mesh)

  mesh.numInterfaces = mesh.numFace - mesh.numBoundaryFaces
  mesh.interfaces = Array{Interface}(mesh.numInterfaces)
  getInterfaceArray(mesh)


  getCoordinatesAndMetrics(mesh, sbp)
#

#=
  mesh.bndryfaces = Array{Boundary}(mesh.numBoundaryFaces)
  getBoundaryArray(mesh)

  mesh.numInterfaces = mesh.numEdge - mesh.numBoundaryFaces
  mesh.interfaces = Array{Interface}(mesh.numInterfaces)
  getInterfaceArray(mesh)

  getCoordinates(mesh, sbp)  # store coordinates of all nodes into array
  getDofNumbers(mesh)  # store dof numbers

  mesh.dxidx = Array{T1}(2, 2, sbp.numnodes, mesh.numEl)
  mesh.jac = Array{T1}(sbp.numnodes, mesh.numEl)
  mappingjacobian!(sbp, mesh.coords, mesh.dxidx, mesh.jac)
=#

  printStats(mesh)

  apf.writeVtkFiles("mesh_complete", mesh.m_ptr)
  return mesh
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
#  resetFaceIt()
  bnd_faces_cnt = 0
  bnd_faces = Array{Int}(mesh.numFace, 2)
  it = apf.MeshIterator(mesh.m_ptr, 2)
  for i=1:mesh.numFace
    face_i = apf.iterate(mesh.m_ptr, it)
    numRegion = apf.countAdjacent(mesh.m_ptr, face_i, 3)  # should be count upward

    if numRegion == 1  # if an exterior edge
      elements = apf.getAdjacent(numRegion)
      elnum = apf.getNumberJ(mesh.el_Nptr, elements[1], 0, 0) + 1

#      println("face number ", i, " is part of only one element")

      bnd_faces_cnt += 1
      bnd_faces[bnd_faces_cnt, 1] = elnum
      bnd_faces[bnd_faces_cnt, 2] = i
    end
  end
  apf.free(mesh.m_ptr, it)


  mesh.boundary_nums = bnd_faces[1:bnd_faces_cnt, :] # copy, bad but unavoidable
  mesh.numBoundaryFaces = bnd_faces_cnt

return nothing

end  # end function


function getBoundaryArray(mesh::PumiMesh3)
# get an array of type Boundary for SBP
# creating an an array of a user defined type seems like a waste of memory operations
# bnd_array is a vector of type Boundary with length equal to the number of edges on the boundary of the mesh

#  mesh.bndryfaces = Array{Boundary}(mesh.numBoundaryFaces)

  for i=1:mesh.numBoundaryFaces
    elnum = mesh.boundary_nums[i,1]
    facenum_global = mesh.boundary_nums[i,2]
    facenum_local = getBoundaryEntityLocalNum(mesh, facenum_global, 2)
    mesh.bndryfaces[i] = Boundary(elnum, facenum_local)
  end

  return nothing
end

function getBoundaryEntityLocalNum(mesh::PumiMesh3, entity_num::Integer, entity_dim::Integer)
# gets the local edge number of a specified edge that is on the boundary
# of the mesh
# entity_num is an the  num from the output of getBoundaryEdgeNums() (ie. the global entity number)
# entity_dim is the dimension of the entity (0 = vertex, 1 = edge, 2 = face)
# the local edge number is the edge number within the element (1st,s 2nd, 3rd...)
  if entity_dim == 1
    entity_i = mesh.edges[entity_num]
  elseif entity_dim == 2
    entity_i = mesh.faces[entity_num]
  else
    println("unsupported entity dimension requested for getBoundaryEntityLocalNum")
    println("requested entity dimension = ", entity_dim)
  end

  # get mesh face or region associated with edge
  apf.countAdjacent(mesh.m_ptr, entity_i, entity_dim+1)
  parent = apf.getAdjacent(1)[1]  # get the single face (not an array)

  local_num = apf.getEntityLocalNumber(mesh.m_ptr, entity_i, parent, entity_dim, entity_dim+1)

  return local_num

end

function numberDofs(mesh::PumiMesh3)
# assign dof numbers to entire mesh
# calculates numNodes, numDof
# assumes mesh elements have already been reordered
# this method minimizes the assembly/disassembly time
# because blocks of entries in res go into res_vec
 
  # calculate total number of nodes
  # Continuous Galerkin only
  numnodes = 0
  for i=1:4
    numnodes += mesh.numNodesPerType[i]*mesh.numEntities[i]
  end

  numDof = mesh.numDofPerNode*numnodes
 
  if (numDof > 2^30)
    println("Warning: too many dofs, renumbering will fail")
  end

  # save numbers to mesh
  mesh.numNodes = numnodes
  mesh.numDof = numDof

  # initally number all dofs as numDof+1 to 2*numDof
  # this allows quick check to see if somthing is labelled or not

#  resetAllIts2(mesh.m_ptr)
  # mesh iterator increment, retreval functions
#  iterators_inc = [incrementVertIt, incrementEdgeIt, incrementFaceIt, incrementElIt]
#  iterators_get = [getVert, getEdge, getFace, getEl]
  num_entities = mesh.numEntities
  num_nodes_entity = mesh.numNodesPerType

  #  curr_dof = 1
  curr_dof = numDof+1
  for etype = 1:length(num_nodes_entity) # loop over entity types
    if (num_nodes_entity[etype] != 0)  # if no nodes on this type of entity, skip
      it = apf.MeshIterator(mesh.m_ptr, etype - 1)
      for entity = 1:num_entities[etype]  # loop over all entities of this type
#	println("  entity number: ", entity)
        entity_ptr = apf.iterate(mesh.m_ptr, it)
#	entity_ptr = iterators_get[etype]()  # get entity

	for node = 1:num_nodes_entity[etype]
#	  println("    node : ", node)
	  for dof = 1:mesh.numDofPerNode
	    apf.numberJ(mesh.dofnums_Nptr, entity_ptr, node-1, dof-1, curr_dof)
#	    println("      entity ", entity_ptr, " labelled ", curr_dof)
	    curr_dof += 1
	  end  # end loop over dof
	end  # end loop over node

#      iterators_inc[etype]()
      end  # end loop over entitiesaa
      apf.free(mesh.m_ptr, it)
    end  # end if 
  end  # end loop over entity types

  # final dof numbering
# move all if statements out one for loop (check only first dof on each node)
  curr_dof = 1
  for i=1:mesh.numEl
#    println("element number: ", i)

    el_i = mesh.elements[i]
    # get vertices, edges for this element
    (verts_i, numVert) = apf.getDownward(mesh.m_ptr, mesh.elements[i], 0)
#    println("verts_i = ", verts_i)
#    println("mesh.verts = ", mesh.verts)
    (edges_i, numEdge) = apf.getDownward(mesh.m_ptr, mesh.elements[i], 1)

    if mesh.numEntitiesPerElement[3] > 1  # don't call get downward if el_i is a face
      (face_i, numFace) = apf.getDownward(mesh.m_ptr, mesh.elements[i], 2)
      regions_i = [mesh.elements[i]]  # make this an array for consistency
      numRegions = 1
    else  # this is a 2d mesh
      face_i = [mesh.elements[i]]
      numFace = 1
      numRegion = 0
    end

    # loop over vertices
    for j=1:numVert  # loop over vertices, edges
#      println("  vertex number: ", j)
      vert_j = verts_i[j]
#      println("  vert_j = ", vert_j)
#      println("  edge_j = ", edge_j)
      for k=1:mesh.numDofPerNode # loop over vertex dofs
#	println("    dof number: ", k)
#	println("    vert_j = ", vert_j)
        dofnum_k = apf.getNumberJ(mesh.dofnums_Nptr, vert_j, 0, k-1)
	if dofnum_k > numDof  # still has initial number
	  # give it new (final) number
	  apf.numberJ(mesh.dofnums_Nptr, vert_j, 0, k-1, curr_dof)
#	  println("numbering a vertex ", curr_dof)
	  curr_dof += 1
	end
      end  # end loop over vertex dofs
    end  # end loop over vertices

    # now loop over edges
    for j=1:numEdge
#      println("edge number: ", j)
      edge_j = edges_i[j]
      # loop over nodes on edge
      for k=1:num_nodes_entity[2]  # loop over all nodes on the edge
	for p=1:mesh.numDofPerNode  # loop over dofs
	  dofnum_p = apf.getNumberJ(mesh.dofnums_Nptr, edge_j, k-1, p-1)
	  if dofnum_p > numDof  # still has initial number
	    # give it new (final) number
	    apf.numberJ(mesh.dofnums_Nptr, edge_j, k-1, p-1, curr_dof)
	    curr_dof += 1
	  end
	end
      end  # en loop over nodes on edge
    end  # end loop over edges

    # now loop over faces
    for j=1:numFace
#      println("face number: ", j)
      face_j = face_i[j]
      for k=1:num_nodes_entity[3]  # loop over nodes on faces
	for p=1:mesh.numDofPerNode  # loop over dofs
	  dofnum_p = apf.getNumberJ(mesh.dofnums_Nptr, face_j, k-1, p-1)
	  if dofnum_p > numDof
	    apf.numberJ(mesh.dofnums_Nptr, face_j, k-1, p-1, curr_dof)
	    curr_dof += 1
	  end
	end
      end  # end loop over face nodes
    end  # end loop over faces

    # now loop over region
    for j=1:numRegions
#      println("region number: ", j)
      region_j = regions_i[j]
      for k=1:num_nodes_entity[4]
	for p=1:mesh.numDofPerNode
	  dofnum_p = apf.getNumberJ(mesh.dofnums_Nptr, region_j, k-1, p-1)
	  if dofnum_p > numDof
	    apf.numberJ(mesh.dofnums_Nptr, region_j, k-1, p-1, curr_dof)
	    curr_dof += 1
	  end
	end
      end  # end loop over nodes on the region
    end  # end loop over region

  end  # end loop over elements


  if (curr_dof -1) != numDof 
    println("Warning: number of dofs assigned is not equal to teh expected number")
    println("number of dofs assigned = ", curr_dof-1, " ,expected number = ", numDof)
  end


#  resetAllIts2(mesh.m_ptr)
  return nothing

end



function getDofNumbers(mesh::PumiMesh3)
# populate array of dof numbers, in same shape as solution array u (or q)

mesh.dofs = Array{Int}(mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)

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


function getGlobalNodeNumbers(mesh::PumiMesh3, elnum::Integer)
# gets global node numbers of all dof on all nodes of the element
# output formap is array [numdofpernode, nnodes]  (each column contains dof numbers for a node)
el_i = mesh.elements[elnum]
#type_i = apf.getType(mesh.m_ptr, el_i)


# use apf::getElementNumbers to ensure everything is in the proper orientation

  nnodes = mesh.numNodesPerElement
  numdof = nnodes*mesh.numDofPerNode
  nums = zeros(Int32, numdof)

apf.getElementNumbers(mesh.dofnums_Nptr, el_i, numdof, nums::Array{Int32, 1})

# reshape into output format
# underlying memory is shared
dofnums_reshaped = reshape(nums, mesh.numDofPerNode, nnodes)

return dofnums_reshaped

end

function getCoordinates(mesh::PumiMesh3, sbp::AbstractSBP)
# populate the coords array of the mesh object

mesh.coords = Array{Float64}(3, sbp.numnodes, mesh.numEl)

println("entered getCoordinates")

coords_i = zeros(3, 4)
coords_it = zeros(4, 3)
for i=1:mesh.numEl  # loop over elements
  
  el_i = mesh.elements[i]
  (sizex, sizey) = size(coords_i)
  apf.getElCoords(mesh.m_ptr, el_i, coords_i, sizex, sizey)  # populate coords

#  println("coords_i = ", coords_i)

  coords_it[:,:] = coords_i.'
 # println("coords_it = ", coords_it)
  mesh.coords[:, :, i] = calcnodes(sbp, coords_it)
#  println("mesh.coords[:,:,i] = ", mesh.coords[:,:,i])
end

return nothing

end

function getInterfaceArray(mesh::PumiMesh3)
# get array of [elementL, elementR, faceL, faceR] for each internal face,
# where elementL and R are the elements that use the face, faceL R are the
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

#  interfaces = Array{typeof(new_interface}(num_int_edges)
  # get number of nodes affecting an edge
  num_face_nodes = apf.countAllNodes(mesh.mshape_ptr, 2)

  # store whether the node map has been calculated for each case
  # index = relative rotation of faces + 1
  nodemap_mask = zeros(Bool, 3)
  nodemaps = Array{Array{Uint8,1}}(3)
  

#  nodemap = collect(num_edge_nodes:(-1):1)

  pos = 1 # current position in interfaces
  for i=1:mesh.numFace
#     println("i = ", i)
#     println("pos = ", pos)
    # get number of elements using the edge
    adjacent_nums, num_adjacent = getAdjacentEntityNums(mesh, i, 2, 3)
#     println("num_adjacent = ", num_adjacent)
#     println("adjacent_nums = ", adjacent_nums)
    if num_adjacent > 1  # internal edge
#       println("this is an internal edge")
      element1 = adjacent_nums[1]
      element2 = adjacent_nums[2]

#      coords_1 = x[:, :, element1]
#      coords_2 = x[:, :, element2]

      coords_1 = zeros(3,4)
      coords_2 = zeros(3,4)
      apf.getElCoords(mesh.m_ptr, mesh.elements[element1], coords_1, 3, 4)
      apf.getElCoords(mesh.m_ptr, mesh.elements[element2], coords_2, 3, 4)

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

      faceL = apf.getEntityLocalNumber(mesh.m_ptr, mesh.faces[i], mesh.elements[elementL], 2, 3)
      faceR = apf.getEntityLocalNumber(mesh.m_ptr, mesh.faces[i], mesh.elements[elementR], 2, 3)
#      edgeL = getFaceLocalNum(mesh, i, elementL)
#      edgeR = getFaceLocalNum(mesh, i, elementR)

      # construct nodemap if needed
      rel_rotate = calcRelRotate(mesh, elementL, elementR, i)
      if !nodemap_mask[rel_rotate + 1]  # if nodemap not yet constructed
	nodemaps[rel_rotate + 1] = constructNodemap(mesh, rel_rotate)
	nodemap_mask[rel_rotate + 1] = true
      end


      mesh.interfaces[pos] = Interface(elementL, elementR, faceL, faceR, nodemaps[rel_rotate + 1])
#       println("updating pos")
      pos += 1

  #    print("\n")

    end  # end if internal edge

#     print("\n")
  end  # end loop over edges

  return nothing

end  # end function



function calcRelRotate(mesh::PumiMesh3, elementL::Integer, elementR::Integer, facenum::Integer)
# calculate the rotation of faceR relative to faceL, where faceL and faceR
# are shared between the two element
# this uses getAligmnet, which describes how to transform the face from
# its current orientation to the canonical orientation as part of its parent
# element
# the difference in the required rotations is the relative rotation
# elementL : global element index
# elementR : global element index
# faceL : local face number of shared face on elementL
# faceR : local face number of shared face on elementR
# face : global face number of the shared face

eL = mesh.elements[elementL]
eR = mesh.elements[elementR]
face = mesh.faces[facenum]

# get rotations to bring faces into canonical orientation
# for each element
# we can ignore flip because one of the faces will always be flipped
which, flipL, rotateL = apf.getAlignment(mesh.m_ptr, eL, face)
which, flipR, rotateR = apf.getAlignment(mesh.m_ptr, eR, face)

if flipL && flipR
  println("Warning both faces are flipped, nodemap may  be incorrect")
end

# sum rotations because they each rotate ccw in their parent
# element's orientation, so they rotate in opposite directions
rel_rotate = rotateR + rotateL
#println("rel_rotate before wrapping = ", rel_rotate)
rel_rotate = wrapNumber(rel_rotate, 0, 2)
#println("rel_rotate after wrapping = ", rel_rotate)

return rel_rotate 

end  # end function


function constructNodemap(mesh::PumiMesh3, rel_rotate::Integer)
# calculate node map from elementL in canonical order to the canonical order
# on elementR
# both orders are canonical because the dofs array ensure everything
# is in canonical order for the element
# element must be a tetrahedron
# make nodemap UInt8?
# calculate total number of nodes on a face
num_nodes_face = 0

numEntitiesPerFace = [3, 3, 1]

for i=1:3
  num_nodes_face += numEntitiesPerFace[i]*mesh.numNodesPerType[i]
end

# get number of nodes classified on edges and face
nnodes_edge = mesh.numNodesPerType[2]
nnodes_face = mesh.numNodesPerType[3]

nodemap = zeros(Uint8, num_nodes_face)

for i=1:3  # loop over verts
  index_i = wrapNumber(i + rel_rotate, 1, 3)  # correct for relative rotation
#  println("after first wrap, index_i = ", index_i)
  index_i = 5 - index_i  # permuate because one triangle is flipped
#  println("pre wrapping, vertex index_i = ", index_i)
  index_i = wrapNumber(index_i, 1, 3)
#  println("post wrapping, vert index_i = ", index_i)
  nodemap[i] = index_i
end

for i=1:3  # loop over edges
  edgeR = 3 - i + 1 - rel_rotate

#  println("before wrap, edgeR = ", edgeR)
  edgeR = wrapNumber(edgeR, 1, 3)
#  println("after wrap, edgeR = ", edgeR)
  for j=1:nnodes_edge
    nodemap_index = 3 + nnodes_edge*(i-1) + j
    nodemap_val = 3 + nnodes_edge*(edgeR-1) + (nnodes_edge -j + 1)
#    println("nodemap_index = ", nodemap_index)
#    println("nodemap_val = ", nodemap_val)
    nodemap[nodemap_index] = nodemap_val
  end
end

#=
# threshold = location of wraparound
threshold = rel_rotate + 1
sum1 = 2 + (mesh.order - 2)*rel_rotate  # index sum below threshold
sum2 = 2 + (mesh.order - 2)*rel_rotate + nnodes_face  # index sum above threshold
=#

if ((rel_rotate % 2) == 0)
  threshold = rel_rotate + 1
  sum1 = nnodes_face - 1
  sum2 = 2*nnodes_face - 1
else  # rel_rotate is odd
  threshold = 0  # always above threshold
  sum1 = 0  # irrelevent
  sum2 = nnodes_face + 1
end


if mesh.order == 3
  threshold = 0
  sum1 = 2  # irrelevent
  sum2 = 2
end

if mesh.order == 4
  if rel_rotate == 0
    threshold = 1
    sum1 = 2
    sum2 = 5
  elseif rel_rotate == 1
    threshold == 0
    sum1 = 4  # irrelevent
    sum2 = 4
  else  # rel_rotate = 2
    threshold = 2
    sum1 = 3
    sum2 = 5
  end
end



nodemap_index = 3 + 3*nnodes_edge
for j=1:nnodes_face  # loop over interior face points
  if j <= threshold  # choose which sum to use
    sum = sum1
  else
    sum = sum2
  end
#  println("face loop j = ", j, " sum = ", sum)
  nodemap[nodemap_index + j] = nodemap_index + sum - j
end

return nodemap

end

