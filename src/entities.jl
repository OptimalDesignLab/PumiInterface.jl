# functions for gathering MeshEntity*s 
"""
  This function calculates the numNodesPerType and typeOffsetsPerElement arrays
  and returns them

  Inputs:
    fshape: a FieldShape*
    dim: the dimensionality of the mesh (2 or 3)
    numTypePerElement: number of each dimension entity per element (see
                       interfaces.md)

    Outputs:
      numNodesPerType: array of length 3 or 4 (2d or 3d), containing the number
                       of nodes on each vert, edge, face (or region).

      typeOffsetsPerElement: array of length dim + 2 containing the index of 
                             the first node of each type.  The last element
                             is one more than the number of nodes

"""
function getNodeInfo{I <: Integer}(fshape::Ptr{Void}, dim::Integer, numTypePerElement::Array{I, 1})

  num_nodes_v = countNodesOn(fshape, 0)  # number of nodes on a vertex
  num_nodes_e = countNodesOn(fshape, 1) # on edge
  num_nodes_f = countNodesOn(fshape, 2) # on face

  if dim == 3
    num_nodes_r = countNodesOn(fshape, apfTET)  # on region
    numNodesPerType = [num_nodes_v, num_nodes_e, num_nodes_f, num_nodes_r]
  else
    numNodesPerType = [num_nodes_v, num_nodes_e, num_nodes_f]
  end
  # count numbers of different things per other thing
  # use for bookkeeping
  typeOffsetsPerElement = zeros(Int, dim+2)
  pos = 1
  typeOffsetsPerElement[1] = pos
  for i=2:(dim + 2)
    pos += numTypePerElement[i-1]*numNodesPerType[i-1]
    typeOffsetsPerElement[i] = pos
  end

  return numNodesPerType, typeOffsetsPerElement
end


# can be generalized trivially
function getBoundaryElements(mesh::PumiMeshDG2, bndry_edges::AbstractArray{Int, 1})
# get the faces corresponding to the boundary edges

bndry_faces = zeros(bndry_edges)
faces = Array(Ptr{Void}, 400)  # equivilent to apf::Up
numbering = mesh.entity_Nptrs[mesh.dim + 1]
for i=1:bndry_edges
  edgenum_i = bndry_edges[i]
  face_i = mesh.faces[edgenum_i]
  numFace = countAdjacent(mesh.m_ptr, face_i, mesh.dim)  # should be count upward
  getAdjacent(faces)
  facenum = getNumberJ(numbering, face_i, 0, 0) + 1
#  facenum = getFaceNumber2(faces[1]) + 1

  bndry_faces[i] = facenum
end

return bndry_faces

end

# add if statements to generalized this to 3D
function getEntityPointers(mesh::PumiMesh)
# get the pointers to all the apf::MeshEntities and put them in arrays
# uses the Numberings to determine what the index in the array of each
# entity


  verts = Array(Ptr{Void}, mesh.numVert)
  edges = Array(Ptr{Void}, mesh.numEdge)
  elements = Array(Ptr{Void}, mesh.numEl)
  entity = Ptr{Void}(0)
  idx = 0
  # get pointers to all MeshEntities
  # also initilize the field to zero
  resetAllIts2()
#  comps = zeros(dofpernode)
  for i=1:mesh.numVert
    entity = getVert()
    idx = getNumberJ(mesh.vert_Nptr, entity, 0, 0) + 1
    verts[idx] = entity
    incrementVertIt()
  end

  for i=1:mesh.numEdge
    entity = getEdge()
    idx = getNumberJ(mesh.edge_Nptr, entity, 0, 0) + 1
    edges[idx] = entity
    incrementEdgeIt()
  end

  if mesh.dim == 3
    faces = Array(Ptr{Void}, mesh.numFace)
    for i=1:mesh.numFace
      entity = getFace()
      idx = getNumberJ(mesh.face_Nptr, entity, 0, 0) + 1
      faces[idx] = entity
      incrementFaceIt()
    end
    for i=1:mesh.numEl
      entity = getEl()
      idx = getNumberJ(mesh.el_Nptr, entity, 0, 0) + 1
      elements[idx] = entity
      incrementElIt()
    end


  else
    faces = edges
    for i=1:mesh.numEl
      entity = getFace()
      idx = getNumberJ(mesh.el_Nptr, entity, 0, 0) + 1
      elements[idx] = entity
      incrementFaceIt()
    end


  end

  resetAllIts2()

  return verts, edges, faces, elements

end  # end getEntityPointers

function getNodeIdx(e_type::Integer, e_idx::Integer, node_idx::Integer, typeOffsetsPerElement, numNodesPerType )
# calculate the local node number (in the list of all nodes on the element)
# from the entity the node is classified on and its index on the entity

# e_type is the type of the entity: 1 = vert, 2 = edge, 3 = face
# node_idx is the index of the node on the entity it is classified on
# e_idx is the index of the entity the node is classified on (ie 1st, 2nd or
# 3rd edge)

  type_offset = typeOffsetsPerElement[e_type]
  nodes_on_type = numNodesPerType[e_type]
  return type_offset -1 + (e_idx - 1)*nodes_on_type + node_idx
end

#------------------------------------------------------------------------------
# allocators for straight sides meshes
"""
  This function allocates the arrays that store the mesh node coordinates,
  the mesh vertex coordinates and assigns them to the corresponding field of
  the mesh object.
"""
function allocateCoordinateArrays{Tmsh}(mesh::PumiMeshDG{Tmsh}, 
                                                 sbp::AbstractSBP)

  num_coord_nodes = mesh.coord_numNodesPerElement
  mesh.coords = Array(Float64, mesh.dim, sbp.numnodes, mesh.numEl)
  mesh.vert_coords = Array(Float64, mesh.dim, num_coord_nodes, mesh.numEl)

  return nothing
end

"""
  This function allocates the arrays that store the dxidx and jac and
  assigns them to the corresponding field of the mesh object.
"""
function allocateMetricsArrays{Tmsh}(mesh::PumiMeshDG{Tmsh}, sbp::AbstractSBP)

  mesh.dxidx = Array(Tmsh, mesh.dim, mesh.dim, sbp.numnodes, mesh.numEl)
  mesh.jac = Array(Tmsh, sbp.numnodes, mesh.numEl)

  return nothing
end

#------------------------------------------------------------------------------
# allocators for curvilinear meshes

function allocateMeshCoordinateArray{Tmsh}(mesh::PumiMeshDG{Tmsh}, sbp::AbstractSBP)

  num_coord_nodes = mesh.coord_numNodesPerElement
  mesh.vert_coords = Array(Float64, mesh.dim, num_coord_nodes, mesh.numEl)

  return nothing
end

function allocateNodeCoordinateArray{Tmsh}(mesh::PumiMeshDG{Tmsh}, sbp::AbstractSBP)


  mesh.coords = Array(Float64, mesh.dim, sbp.numnodes, mesh.numEl)

  return nothing
end


function allocateCurvilinearCoordinateAndMetricArrays{Tmsh}(mesh::PumiMeshDG{Tmsh},                                                         sbp::AbstractSBP)

  dim = mesh.dim
  sbpface = mesh.sbpface
  
  mesh.coords = zeros(Tmsh, mesh.dim, mesh.numNodesPerElement, mesh.numEl)
  mesh.dxidx = zeros(Tmsh, mesh.dim, mesh.dim, mesh.numNodesPerElement, mesh.numEl)
  mesh.jac = zeros(Tmsh, mesh.numNodesPerElement, mesh.numEl)

  # these arrays are not used for 3D curvilinear meshes

  # interior faces
  mesh.dxidx_face = zeros(Tmsh, 0, 0, 0, 0)
  mesh.jac_face = zeros(Tmsh, 0, 0)

  # boundary faces
  mesh.dxidx_bndry = zeros(Tmsh, 0, 0, 0, 0)
  mesh.jac_bndry = zeros(Tmsh, 0, 0)

  # parallel shared faces
  mesh.dxidx_sharedface = Array(Array{Tmsh, 4}, 0)
  mesh.jac_sharedface = Array(Array{Tmsh, 2}, 0)


  return nothing
end
#------------------------------------------------------------------------------
# straight-sided coordinate and metric functions

"""
  This function gets the coordinates of both the mesh coordinate field and the
  node field and stores them into the member fields of the mesh object.  
  It uses a smart allocator to allocate
  these arrays needed, or avoid doing so if they have already been allocated.
  See allocateCoordinateArrays for which arrays are populated by
  this function.  Non curvilinear meshes only.
"""
#TODO: stop using slice notation
# can be generalized with numVertsPerElement
function getCoordinates(mesh::PumiMeshDG, sbp::AbstractSBP)
# populate the coords array of the mesh object

#  mesh.coords = Array(Float64, mesh.dim, sbp.numnodes, mesh.numEl)
#  mesh.vert_coords = Array(Float64, mesh.dim, nvert_per_el, mesh.numEl)

  if !isFieldDefined(mesh, :coords, :vert_coords)
    allocateCoordinateArrays(mesh, sbp)
  end

  nvert_per_el = mesh.numTypePerElement[1]
  @assert size(mesh.coords) == (mesh.dim, sbp.numnodes, mesh.numEl)
  @assert size(mesh.vert_coords) == (mesh.dim, nvert_per_el, mesh.numEl)
  @assert mesh.coord_order == 1

  #println("entered getCoordinates")
  numVertsPerElement = mesh.numTypePerElement[1]
  coords_i = zeros(mesh.dim ,numVertsPerElement)
  coords_it = zeros(numVertsPerElement, mesh.dim)
  for i=1:mesh.numEl  # loop over elements
    
    el_i = mesh.elements[i]
    getAllEntityCoords(mesh.m_ptr, el_i, coords_i)
    mesh.vert_coords[:, :, i] = coords_i[:, :]
    coords_it[:,:] = coords_i[1:mesh.dim, :].'
    mesh.coords[:, :, i] = SummationByParts.SymCubatures.calcnodes(sbp.cub, coords_it)
  end

  return nothing

end

"""
  This function calculates the dxidx and jac values for the entire mesh.  It
  uses a smart allocator to allocate the arrays if needed.  Non-curvilinear
  meshes only.
"""
function getMetrics(mesh::PumiMeshDG, sbp::AbstractSBP)

  if !isFieldDefined(mesh, :jac, :dxidx)
    allocateMetricsArrays(mesh, sbp)
  end

  mappingjacobian!(sbp, mesh.coords, mesh.dxidx, mesh.jac)

  return nothing
end


    

function getElementCoords(mesh::PumiMesh2D, entity::Ptr{Void}, coords::AbstractMatrix)
  # coords must be 3 x numVertsPerElement
  sx, sy = size(coords)
  getFaceCoords(entity, coords, sx, sy)
end


function getElementCoords(mesh::PumiMesh3D, entity::Ptr{Void}, coords::AbstractMatrix)

  sx, sy = size(coords)
  getElCoords(entity, coords, sx, sy)
end



#TODO: stop using slice notation
function getCoordinatesAndMetrics{Tmsh}(mesh::PumiMeshCG{Tmsh}, sbp::AbstractSBP)
# populate the coords array of the mesh object
mesh.coords = Array(Float64, mesh.dim, sbp.numnodes, mesh.numEl)

#println("entered getCoordinates")

numVertsPerElement = mesh.numTypePerElement[1]
coords_i = zeros(3,numVertsPerElement)
coords_it = zeros(numVertsPerElement, mesh.dim)
for i=1:mesh.numEl  # loop over elements
  el_i = mesh.elements[i]
  getElementCoords(mesh, el_i, coords_i)

  coords_it[:,:] = coords_i[1:mesh.dim, :].'
  mesh.coords[:, :, i] = calcnodes(sbp, coords_it)
end

  mesh.dxidx = Array(Tmsh, mesh.dim, mesh.dim, sbp.numnodes, mesh.numEl)
  mesh.jac = Array(Tmsh, sbp.numnodes, mesh.numEl)
  mappingjacobian!(sbp, mesh.coords, mesh.dxidx, mesh.jac)


return nothing

end


"""
  Populates the input array with the coordinates of the nodes on list of
  boundary faces.  Non-curvilinear meshes only.

  Inputs:
    mesh: a mesh object
    bndryfaces: an array of Boundary objects to calculate the coordinates of

  Inputs/Outputs:
    coords_bndry: an array dim x numfacenodes x length(bndryfaces) to populate
                  with the coordinates of the nodes on the faces
"""
function getBndryCoordinates{Tmsh}(mesh::PumiMeshDG2{Tmsh}, 
                             bndryfaces::Array{Boundary}, 
                             coords_bndry::Array{Tmsh, 3})
# calculate the coordinates on the boundary for the specified faces
# and store in coords_bndry

#  println("----- Entered getBndryCoordinates -----")
  sbpface = mesh.sbpface

  coords_i = zeros(3, 3)
  coords_it = zeros(3, 2)
  coords_edge = zeros(2, 2)

  for i=1:length(bndryfaces)
    bndry_i = bndryfaces[i]

    el = bndry_i.element
    el_ptr = mesh.elements[el]
    face = bndry_i.face

    sizex, sizey = size(coords_i)
    getFaceCoords(el_ptr, coords_i, sizex, sizey)

    coords_it[:, :] = coords_i[1:2, :].'

    # extract the needed vertex coords
#    v1 = facemap[1, face]
#    v2 = facemap[2, face]
    v1 = face
    v2 = mod(face,3) + 1
    coords_edge[1, 1] = coords_it[v1, 1]
    coords_edge[1, 2] = coords_it[v1, 2]
    coords_edge[2, 1] = coords_it[v2, 1]
    coords_edge[2, 2] = coords_it[v2 ,2]

    coords_bndry[:, :, i] = SummationByParts.SymCubatures.calcnodes(sbpface.cub, coords_edge)

  end

end

function getBndryCoordinates{Tmsh}(mesh::PumiMesh3DG{Tmsh}, bndryfaces::Array{Boundary}, coords_bndry::Array{Tmsh, 3})

#  println("----- entered getBndryCoordinates -----")
  sbpface = mesh.sbpface
  coords_i = zeros(3, mesh.numEntitiesPerType[1])
#  coords_i = zeros(3, mesh.numTypePerElement[1])
  coords_face = zeros(3, 3)
  vertmap = mesh.topo.face_verts
  for i=1:length(bndryfaces)
    bndry_i = bndryfaces[i]
    el = bndry_i.element
    el_ptr = mesh.elements[el]
    face = bndry_i.face

    getElementCoords(mesh, el_ptr, coords_i)

    # extract the vertices of the face
    for v=1:3
      vidx = vertmap[v, face]
      for dim=1:3
        coords_face[v, dim] = coords_i[dim, vidx]
      end
    end

    coords_bndry[:, :, i] = SummationByParts.SymCubatures.calcnodes(sbpface.cub, coords_face)
  end

  return nothing
end

function getInterfaceCoordinates{Tmsh}(mesh::PumiMeshDG2{Tmsh}, 
                             bndryfaces::Array{Interface}, 
                             coords_bndry::Array{Tmsh, 3})
# calculate the coordinates on the boundary for the specified faces
# and store in coords_bndry

#  println("----- Entered getBndryCoordinates -----")
  sbpface = mesh.sbpface

  coords_i = zeros(3, 3)
  coords_it = zeros(3, 2)
  coords_edge = zeros(2, 2)  # one point per row

  for i=1:length(bndryfaces)
    bndry_i = bndryfaces[i]

    el = bndry_i.elementL
    el_ptr = mesh.elements[el]
    face = bndry_i.faceL

    sizex, sizey = size(coords_i)
    getFaceCoords(el_ptr, coords_i, sizex, sizey)

    coords_it[:, :] = coords_i[1:2, :].'

    # extract the needed vertex coords
#    v1 = facemap[1, face]
#    v2 = facemap[2, face]
    v1 = face
    v2 = mod(face,3) + 1
    coords_edge[1, 1] = coords_it[v1, 1]
    coords_edge[1, 2] = coords_it[v1, 2]
    coords_edge[2, 1] = coords_it[v2, 1]
    coords_edge[2, 2] = coords_it[v2 ,2]

    coords_bndry[:, :, i] = SummationByParts.SymCubatures.calcnodes(sbpface.cub, coords_edge)

  end

end



function getInterfaceCoordinates{Tmsh}(mesh::PumiMesh3DG{Tmsh}, bndryfaces::Array{Interface}, coords_bndry::Array{Tmsh, 3})

#  println("----- entered getInterfaceCoordinates -----")
  sbpface = mesh.sbpface
  coords_i = zeros(3, mesh.numTypePerElement[1])
  coords_face = zeros(3, 3)
  vertmap = mesh.topo.face_verts
  for i=1:length(bndryfaces)
    bndry_i = bndryfaces[i]
    el = bndry_i.elementL
    el_ptr = mesh.elements[el]
    face = bndry_i.faceL

    getElementCoords(mesh, el_ptr, coords_i)

    # extract the vertices of the face
    for v=1:3
      vidx = vertmap[v, face]
      for dim=1:3
        coords_face[v, dim] = coords_i[dim, vidx]
      end
    end


    coords_bndry[:, :, i] = SummationByParts.SymCubatures.calcnodes(sbpface.cub, coords_face)
  end

  return nothing
end



#=
function getElementVertCoords(mesh::PumiMesh, elnum::Integer, coords::AbstractArray{Float64,2})
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

# generalizes with numVertsPerElement (or just get rid of this?)
function getElementVertCoords(mesh::PumiMesh,  elnums::Array{Int,1})
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
=#


#------------------------------------------------------------------------------
# curvilinear functions

"""
  This function gets the coordinates that define the mesh 
  (not the node coordinates) and puts them the mesh.vert_coords field
  of the mesh object.  For each element, the ordering of the coordinates
  is verts, then edges, then faces, then regions.  This function uses
  smart allocators to allocate the array if needed.

  Input:
    mesh: a DG mesh
    sbp: an SBP operators
"""
function getMeshCoordinates{Tmsh}(mesh::PumiMeshDG{Tmsh}, sbp::AbstractSBP)

  if !isFieldDefined(mesh, :vert_coords)
    allocateMeshCoordinateArray(mesh, sbp)
  else
    fill!(mesh.vert_coords, 0.0)
  end

  for i=1:mesh.numEl
    el_i = mesh.elements[i]
    coords_i = sview(mesh.vert_coords, :, :, i)
    getAllEntityCoords(mesh.m_ptr, el_i, coords_i)
  end

  return nothing
end

"""
  This function calculates the fields of the mesh that hold coordinates of the
  face nodes for boundaries, interfaces, and sharedfaces.  This function uses
  smart allocators to allocate the arrays if needed
"""
function getFaceCoordinatesAndNormals{Tmsh}(mesh::PumiMeshDG{Tmsh}, sbp::AbstractSBP)
  if !isFieldDefined(mesh, :coords_bndry, :coords_sharedface)
    allocateFaceCoordinates(mesh)
  else  # zero out arrays that might need it
    fill!(mesh.coords_bndry, 0.0)
    fill!(mesh.coords_interface, 0.0)
    for i=1:mesh.npeers
      fill!(mesh.coords_sharedface[i], 0.0)
    end
  end

  if !isFieldDefined(mesh, :nrm_bndry, :nrm_face, :nrm_sharedface)
    allocateNormals(mesh, sbp)
  else
    fill!(mesh.nrm_bndry, 0.0)
    fill!(mesh.nrm_face, 0.0)
    for i=1:mesh.npeers
      fill!(mesh.coords_sharedface[i], 0.0)
    end
  end

  calcFaceCoordinatesAndNormals(mesh, sbp, mesh.bndryfaces, mesh.coords_bndry, 
                               mesh.nrm_bndry)
  calcFaceCoordinatesAndNormals(mesh, sbp, mesh.interfaces, 
                               mesh.coords_interface, mesh.nrm_face)
  for i=1:mesh.npeers
    calcFaceCoordinatesAndNormals(mesh, sbp, mesh.bndries_local[i], 
                               mesh.coords_sharedface[i], mesh.nrm_sharedface[i])
  end

  return nothing
end

"""
  This is the inner function used by getFaceCoordinatesandNormals.  It
  gets the coordinates of the face nodes and their normal vectors for
  a given array of Interfaces/Boundaries

  Inputs:
    mesh: a DG mesh object
    faces: an array of Boundaries or Interfaces

  Inputs/Outputs
    coords_face: array Tdim x numNodesPerFace x length(faces)
    nrm_face: array, same size as coords_face

  Aliasing restrictions: none, although it would be weird if coords_face and
                         nrm_face aliased
"""
function calcFaceCoordinatesAndNormals{Tmsh, I <: Union{Boundary, Interface}}(
                    mesh::PumiMeshDG{Tmsh}, sbp::AbstractSBP,
                    faces::AbstractArray{I, 1}, 
                    coords_face::AbstractArray{Tmsh, 3}, 
                    nrm_face::AbstractArray{Tmsh, 3})

  #TODO: do this in a block format to avoid a temporary array of size O(numEl)
  nfaces = length(faces)
  numNodesPerElement = mesh.coord_numNodesPerElement
  numNodesPerFace = mesh.coord_numNodesPerType[mesh.dim]

  # some temporary arrays
  down_faces = Array(Ptr{Void}, 12)
  coords_lag_face = Array(Float64, mesh.dim, mesh.coord_numNodesPerFace, nfaces)

  # get the parametic coordinates of the face nodes
  face_xi = mesh.coord_facexi
  ref_verts = baryToXY(face_xi, mesh.sbpface.vtx)
  for i=1:nfaces
    el_i = getElementL(faces[i])
    face_i = getFaceL(faces[i])
    el_ptr = mesh.elements[el_i]

    coords_i = sview(coords_lag_face, :, :, i)
    getMeshFaceCoordinates(mesh, el_i, face_i, coords_i)
    #=
    # get the MeshEntity* for the face
    getDownward(mesh.m_ptr, el_ptr, mesh.dim-1, down_faces)
    face_ptr = down_faces[face_i]

    # get the lagrangian node coordinates
    # this might not be in the orientation specified by mesh.topo.vertmap
    getAllEntityCoords(mesh.m_ptr, face_ptr, coords_i)
    =#
  end

  # call SBP
  calcFaceNormals!(mesh.sbpface, mesh.coord_order, ref_verts, coords_lag_face, 
                   coords_face, nrm_face)

  fixOutwardNormal(mesh, faces, nrm_face)

  return nothing
end

"""
  This function check to make sure each face normal vector is oriented
  outwards, and flips it if needed.
"""
function fixOutwardNormal{I <: Union{Boundary, Interface}, Tmsh}(mesh, 
                          faces::AbstractArray{I, 1},
                          nrm_face::AbstractArray{Tmsh, 3})


  tmp = zeros(3)  # temporary vector to hold coordinates
  topo = mesh.topo
  numVertPerElement = mesh.numTypePerElement[1]
  numVertPerFace = numVertPerElement - 1

  # temporary arrays
  el_verts = Array(Ptr{Void}, numVertPerElement)
  other_vert_coords = zeros(mesh.dim)
  face_verts = Array(Ptr{Void}, numVertPerElement - 1)
  face_vert_coords = zeros(mesh.dim, numVertPerFace)

  nfaces = length(faces)
  for i=1:nfaces
    iface_i = faces[i]
    elnum = getElementL(iface_i)
    facenum_local = getFaceL(iface_i)

    el_i = mesh.elements[elnum]
    getDownward(mesh.m_ptr, el_i, 0, el_verts)

    for j=1:numVertPerFace
      face_verts[j] = el_verts[topo.face_verts[j, facenum_local]]
      getPoint(mesh.m_ptr, face_verts[j], 0, tmp)

      for p=1:mesh.dim
        face_vert_coords[p, j] = tmp[p]
      end
    end

    # get the vert not on the face
    other_vert = Ptr{Void}(0)
    for j=1:numVertPerElement
      if !(el_verts[j] in face_verts)
        other_vert = el_verts[j]
      end
    end

    getPoint(mesh.m_ptr, other_vert, 0, tmp)
    for p=1:mesh.dim
      other_vert_coords[p] = tmp[p]
    end

    # check that the face normal is in the opposite direction as the
    # vectors from a vertex on the face to the vertex not on the face
    for j=1:mesh.numNodesPerFace
      outward_count = 0  # count number of calculations that showed outward
      for k=1:numVertPerFace
        val = zero(Float64)
        for p=1:mesh.dim
          r1_p = other_vert_coords[p] - face_vert_coords[p, k]
          val += nrm_face[p, j, i]*r1_p  # accumulate dot product
        end
        
        if val < 0
          outward_count += 1
        end  # end if
      end  # end k

      # reverse face if needed, throw exception is unclear
      if outward_count == 0  # no calculation found the normal is outward
        for p=1:mesh.dim
          nrm_face[p, j, i] = -nrm_face[p, j, i]
        end
      # some, but not all, calculations found outward
      elseif outward_count < numVertPerFace          
        throw(ErrorException("face $iface_i, node $j has indeterminate orientation"))
      end  # else the face is oriented outwards, do nothing

    end  # end j

  end  # end loop i

   return nothing
end



"""
  This function calculates the node coordinates, scaled mapping jacobian, and
  mapping jacobian determinant for a curvilinear mesh and stores them to
  the fields of the mesh object.  This function uses smart allocators to
  allocate the arrays if needed

  Inputs:
    mesh: a DG mesh object
    sbp: an SBP operator
"""
function getCurvilinearCoordinatesAndMetrics{Tmsh}(mesh::PumiMeshDG{Tmsh}, 
                                                  sbp::AbstractSBP)

  if !isFieldDefined(mesh, :coords, :dxidx, :jac)
    allocateCurvilinearCoordinateAndMetricArrays(mesh, sbp)
  else
    fill!(mesh.coords, 0.0)
    fill!(mesh.dxidx, 0.0)
    fill!(mesh.jac, 0.0)
  end

  ref_vtx = baryToXY(mesh.coord_xi, sbp.vtx)
  calcMappingJacobian!(sbp, mesh.coord_order, ref_vtx, mesh.vert_coords, mesh.coords, mesh.dxidx, mesh.jac)
end
#------------------------------------------------------------------------------
# misc. helper functions

function getGlobalNodeNumbers(mesh::PumiMesh, elnum::Integer; getdofs=true)

  if getdofs
    numDofPerNode = mesh.numDofPerNode
  else
    numDofPerNode = 1
  end
  dofnums = zeros(Int32, numDofPerNode, mesh.numNodesPerElement)
  getGlobalNodeNumbers(mesh, elnum, dofnums)

  return dofnums
end

# should be generalizable with entity counts
function getGlobalNodeNumbers(mesh::PumiMesh, elnum::Integer, dofnums::AbstractArray{Int32}; getdofs=true)
# gets global node numbers of all dof on all nodes of the element
# output formap is array [numdofpernode, nnodes]  (each column contains dof numbers for a node)
# if getdofs=false, then dofnums need only have 1 column, or can be a vector
# dofnums must be passable to C
# 2D only
# getdofs specifies whether to get dof numbers or node numbers
# 

el_i = mesh.elements[elnum]
type_i = getType(mesh.m_ptr, el_i)  # what is this used for?

#println("elnum = ", elnum)
# calculate total number of nodes
#nnodes = 3 + 3*mesh.numNodesPerType[2] + mesh.numNodesPerType[3]
nnodes = mesh.numNodesPerElement

if getdofs
  numdof = nnodes*mesh.numDofPerNode  # what is this used for?
  numbering_ptr = mesh.dofnums_Nptr
else
  numdof = nnodes
  numbering_ptr = mesh.nodenums_Nptr
end

# get entities in the Pumi order
node_entities = getNodeEntities(mesh.m_ptr, mesh.mshape_ptr, el_i)

# get node offsets in the SBP order
node_offsets = sview(mesh.elementNodeOffsets, :, elnum)
#node_offsets = sview(mesh.elementNodeOffsets[:, elnum])

PumiInterface.getDofNumbers(numbering_ptr, node_entities, node_offsets, mesh.nodemapPumiToSbp, el_i, dofnums)  # C implimentation

return nothing


end


# is this really needed?
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


if typeof(mesh) <: Union{PumiMeshDG2, PumiMesh2}

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

# this can be generalized with mesh.dim
function getBoundaryFaceLocalNum(mesh::PumiMesh, edge_num::Integer)
# gets the local edge number of a specified edge that is on the boundary
# of the mesh
# edge_num is an edge num from the output of getBoundaryEdgeNums() (ie. the global edge number)
# the local edge number is the edge number within the element (1st,s 2nd, 3rd...)
  face_i = mesh.faces[edge_num]

  # get mesh face associated with edge
  countAdjacent(mesh.m_ptr, face_i, mesh.dim)
  el = getAdjacent(1)[1]  # get the single face (not an array)

#  elnum_i = getNumberJ(mesh.el_Nptr, face, 0, 0) + 1
#  facenum_i = getFaceNumber2(face) + 1  # convert to 1 based indexing
  down_faces = Array(Ptr{Void}, 12)
  numedges = getDownward(mesh.m_ptr, el, mesh.dim-1, down_faces)
  facenum_local = 0
  for j = 1:mesh.numFacesPerElement  # find which local edge is edge_i
    if down_faces[j] == face_i
      facenum_local = j
    end
  end

  return facenum_local

end

# this can be generalized with some topology and dimension information
function getFaceLocalNum(mesh::PumiMesh, edge_num::Integer, element_num::Integer)
# find the local edge number of a specified edge on a specified element

  edge = mesh.faces[edge_num]
  element = mesh.elements[element_num]

  down_edges = Array(Ptr{Void}, 12)
   numedges = getDownward(mesh.m_ptr, element, mesh.dim-1, down_edges)

  edgenum_local = 0
   for j = 1:mesh.numFacesPerElement  # find which local edge is edge_i
    if down_edges[j] == edge
      edgenum_local = j
    end
  end

  return edgenum_local

end

function getElementVertMap(mesh::PumiMesh)
# get array that maps from elements to their vertex numbers
# the array is numVertPerElement x numEl
  numVertPerElement = mesh.numTypePerElement[1]

  elvertmap = zeros(Int32, numVertPerElement, mesh.numEl)
  verts = Array(Ptr{Void}, 12)  # apf::Downward
  vert_dim = 0

  for i=1:mesh.numEl
    el_i = mesh.elements[i]
    getDownward(mesh.m_ptr, el_i,  vert_dim, verts)

    for j=1:numVertPerElement
      vert_j = verts[j]
      elvertmap[j, i] = getNumberJ(mesh.vert_Nptr, vert_j, 0, 0) + 1
    end
  end

  return elvertmap
end

# get the coordinates of the nodes on a given face (in the coordinate field, not
# the solution field)
# this only works up to second order, and uses the ElementTopology to figure
# out the orientation of the edge
# elnum is the global element number
# facenum is the local face number
"""
  This function gets the coordinates of the coordinate field nodes of the
  face of an element, in the orientation specified by mesh.topo.  This only
  works up to second order.  In 3D it makes the assumption that edge 1 of
  the triangle is defined from v1 -> v2, edge 2 is v2 -> v3, and edge 3 is
  v3 -> v1.

  Inputs:
    mesh: a 2D mesh
    elnum: the global element number
    facenum: the local face number

  Inputs/Outputs:
    coords: array to be populated with the coordinates, Tdim x 
            number of coordinate nodes on a face (2 for linear 2D, 
            3 for quadratic 2D)
"""
function getMeshFaceCoordinates(mesh::PumiMesh2DG, elnum, facenum, coords::AbstractMatrix)

# coords = a 2 x number of nodes on the face array to be populated with the
# coordinates

  vertmap = mesh.topo.face_verts
  el = mesh.elements[elnum]
  el_verts = Array(Ptr{Void}, 12)
  coords_tmp = Array(Float64, 3)
  
  getDownward(mesh.m_ptr, el, 0, el_verts)

  v1 = el_verts[vertmap[1, facenum]]
  v2 = el_verts[vertmap[2, facenum]]

  getPoint(mesh.m_ptr, v1, 0, coords_tmp)

  coords[1, 1] = coords_tmp[1]
  coords[2, 1] = coords_tmp[2]

  getPoint(mesh.m_ptr, v2, 0, coords_tmp)
  coords[1, 2] = coords_tmp[1]
  coords[2, 2] = coords_tmp[2]

  if hasNodesIn(mesh.coordshape_ptr, 1)
    getDownward(mesh.m_ptr, el, 1, el_verts)
    edge = el_verts[facenum]
    getPoint(mesh.m_ptr, edge, 0, coords_tmp)

    coords[1, 3] = coords_tmp[1]
    coords[2, 3] = coords_tmp[2]
  end

  return nothing
end

#TODO: 3D version
