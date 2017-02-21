# entry point for calculating node coordinates and metrics for either
# linear or curvilinear meshes
# also functions used by both linear and curvilienar calulcation

# files containing the code that do the calculations
include("metrics_linear.jl")
include("metrics_curvilinear.jl")

#------------------------------------------------------------------------------
# entry point for both linear and curvilinear metric calculation
"""
  This function gets all the node coordinates, face coordinates, normal
  vectors, and metrics terms for both linear and curvilinear elements

  For linear meshes, dxidx and jac are calculated at the volume nodes and
  interpolated to the face.  For curvilinear meshes the face normals are
  calculated directly, and dxidx and jac are not calculated (the arrays have
  dimensions of all zeros).  This is because it is not clear how to define
  jac and dxidx at the face such that they satisfy the metric invarients.
  All code should eventually use mesh.nrm_face rather than dxidx_face and
  jac_face (and the corresponding fields for boundaryfaces and shared faces),
  but until that is done, linear meshes will still calculate the face
  quantities but curvilinear meshes will not

  All functions called by this function use smart allocators, so it is
  efficient to call this function multiple times (for example, after
  updating the mesh coordinates)

  Fields populated by this function:
    Linear meshes:
      vert_coords
      coords
      coords_bndry
      coords_interface
      coords_sharedface
      dxidx
      dxidx_face
      dxidx_sharedface
      dxidx_bndry
      jac
      jac_face
      jac_sharedface
      jac_bndry
      nrm_bndry
      nrm_face
      nrm_sharedface
      min_el_size
      volume

    Curvilinear meshes:
      vert_coords
      coords
      coords_bndry
      coords_interface
      coords_sharedface
      dxidx
      jac
      nrm_bndry
      nrm_face
      nrm_sharedface
      min_el_size
      volume


"""
function getAllCoordinatesAndMetrics(mesh, sbp)

  if mesh.coord_order == 1  # DEBUGGING: disable linear calculation
    getCoordinates(mesh, sbp)  # store coordinates of all nodes into array
    getMetrics(mesh, sbp)

    if mesh.isInterpolated
      interpolateCoordinatesAndMetrics(mesh)
    end

    getFaceNormals(mesh, sbp)
  else  # curvilinear

    # do other things
    # getMeshCoordinates
    getMeshCoordinates(mesh, sbp)
    getFaceCoordinatesAndNormals(mesh, sbp)
    getCurvilinearCoordinatesAndMetrics(mesh, sbp)

  end

  mesh.min_el_size = getMinElementSize(mesh)
  mesh.volume = calcVolumeIntegral(mesh, sbp)

  return nothing
end

#------------------------------------------------------------------------------
# functions common to linear and curvilinear
"""
  This function allocates the arrays of normal vector for interfaces, boundary
  faces, and shared faces.

  For interfaces, the normal is calculated for elementL
"""
function allocateNormals{Tmsh}(mesh::PumiMeshDG{Tmsh}, sbp)

  dim = mesh.dim
  numfacenodes = mesh.numNodesPerFace
  mesh.nrm_bndry = Array(Tmsh, dim, numfacenodes, mesh.numBoundaryFaces )


  mesh.nrm_face = Array(Tmsh, mesh.dim, numfacenodes, mesh.numInterfaces)
  mesh.nrm_sharedface = Array(Array{Tmsh, 3}, mesh.npeers)
  for i=1:mesh.npeers
    mesh.nrm_sharedface[i] = Array(Tmsh, mesh.dim, numfacenodes, mesh.peer_face_counts[i])
  end

  return nothing
end


