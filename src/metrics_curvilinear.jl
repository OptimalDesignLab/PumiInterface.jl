# functions for calculating node coordinates and metrics for curvilinear
# meshes

#------------------------------------------------------------------------------
# allocators for curvilinear meshes

"""
  Allocates mesh.vert_coords
"""
function allocateMeshCoordinateArray{Tmsh}(mesh::PumiMeshDG{Tmsh}, sbp::AbstractSBP)

  num_coord_nodes = mesh.coord_numNodesPerElement

  if !isFieldDefined(mesh, :vert_coords)
    mesh.vert_coords = Array(Float64, mesh.dim, num_coord_nodes, mesh.numEl)
  else
    fill!(mesh.vert_coords, 0.0)
  end

  return nothing
end

"""
  Allocates the following fields of mesh:
    coords
    dxidx
    jac

  Also allocates the following fields with size zero (because they are not used):
   dxidx_face
   jac_face
   dxidx_bndry
   jac_bndry
   dxidx_sharedface
   jac_sharedface
"""
function allocateCurvilinearCoordinateAndMetricArrays{Tmsh}(mesh::PumiMeshDG{Tmsh},                                                         sbp::AbstractSBP)

  dim = mesh.dim
  sbpface = mesh.sbpface
  
  if !isFieldDefined(mesh, :coords, :dxidx, :jac)
    mesh.coords = zeros(Tmsh, mesh.dim, mesh.numNodesPerElement, mesh.numEl)
    mesh.dxidx = zeros(Tmsh, mesh.dim, mesh.dim, mesh.numNodesPerElement, mesh.numEl)
    mesh.jac = zeros(Tmsh, mesh.numNodesPerElement, mesh.numEl)

    # these arrays are not used for curvilinear meshes

    # interior faces
    mesh.dxidx_face = zeros(Tmsh, 0, 0, 0, 0)
    mesh.jac_face = zeros(Tmsh, 0, 0)

    # boundary faces
    mesh.dxidx_bndry = zeros(Tmsh, 0, 0, 0, 0)
    mesh.jac_bndry = zeros(Tmsh, 0, 0)

    # parallel shared faces
    mesh.dxidx_sharedface = Array(Array{Tmsh, 4}, 0)
    mesh.jac_sharedface = Array(Array{Tmsh, 2}, 0)

  else
    fill!(mesh.coords, 0.0)
    fill!(mesh.dxidx, 0.0)
    fill!(mesh.jac, 0.0)
  end


  return nothing
end

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

  allocateMeshCoordinateArray(mesh, sbp)

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

  allocateFaceCoordinates(mesh)

  allocateNormals(mesh, sbp)

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

  allocateCurvilinearCoordinateAndMetricArrays(mesh, sbp)

  ref_vtx = baryToXY(mesh.coord_xi, sbp.vtx)
  calcMappingJacobian!(sbp, mesh.coord_order, ref_vtx, mesh.vert_coords, mesh.coords, mesh.dxidx, mesh.jac)
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
