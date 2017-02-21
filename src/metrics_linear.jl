# functions for calculating the node coordinates and metrics of straight
# sided elements (ie. meshes with linear coordinate fields)

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


"""
  This function dispatches on the type of the mesh to get the coordinates of
  a mesh element (a face in 2d or a region in 3d)
"""
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
"""
  This function calculates the node coordinates, dxidx, jac, for a 
  CG mesh
"""
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
  boundary faces.  Non-curvilinear meshes only.  In 3D it uses mesh.topo
  to get the vertices of the correct order.

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

"""
  This is basically the same as getBndryCoordinates, except it operates on
  interfaces.  It calculates all quantities from the perspective of elementL

  I think this can be combined with getBndryCoordinates
"""
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

"""
  This function calculates the face normal vectors for all interfaces, 
  boundary faces, and shared faces.  A smart allocator is used to allocate
  the arrays if needed.  Non curvilinear meshes only.

  The interpolated dxidx arrays must be populated before this function is
  called.

  Inputs:
    mesh: a PumiMeshdG
    sbp: an SBP operator
"""
function getFaceNormals(mesh::PumiMeshDG, sbp)

  if !isFieldDefined(mesh, :nrm_bndry, :nrm_face, :nrm_sharedface)
    allocateNormals(mesh, sbp)
  end

  calcFaceNormal(mesh, sbp, mesh.bndryfaces, mesh.dxidx_bndry, mesh.nrm_bndry)
  calcFaceNormal(mesh, sbp, mesh.interfaces, mesh.dxidx_face, mesh.nrm_face)
  for i=1:mesh.npeers
    calcFaceNormal(mesh, sbp, mesh.bndries_local[i], mesh.dxidx_sharedface[i],
                   mesh.nrm_sharedface[i])
  end

  return nothing
end



"""
  This function calculates the face normal for a set of faces given the
  dxidx values for each face node and the facenormal in parametric space

  Inputs:
    mesh: a DG mesh
    sbp: an SBP operator containing the face normals in parametric space
    faces: list of Interfaces or Boundaries, should be dim x numFacesPerElement
    dxidx: d(xi)/d(x) at the face nodes, should be dim x dim x numfaceNodes x 
           length(faces)

  Inputs/Outputs:
    nrm: array, dim x numfacenodes x length(faces) containing the face normal
         vector in x-y space at each face node of each face

  Aliasing restrictions: none
"""
function calcFaceNormal{Tmsh, Tbndry <: Union{Boundary, Interface}}(mesh::PumiMeshDG, sbp, faces::AbstractArray{Tbndry, 1}, dxidx::AbstractArray{Tmsh, 4}, nrm::AbstractArray{Tmsh, 3})

  nfaces = length(faces)
  Tdim = mesh.dim
  for i=1:nfaces
    faceL = getFaceL(faces[i])

    for j=1:mesh.numNodesPerFace
      for dim=1:Tdim
        nrm_d = zero(Tmsh)
        for d=1:Tdim
          nrm_d += mesh.sbpface.normal[d, faceL]*dxidx[d, dim, j, i]
        end  # end loop d
        nrm[dim, j, i] = nrm_d
      end  # end loop dim
    end  # end loop j

  end  # end loop i

  return nothing
end
