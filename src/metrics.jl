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

  This function is collective mesh.comm.

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
      remote_metrics
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
      remote_metrics
      volume


"""
function getAllCoordinatesAndMetrics(mesh, sbp, opts)

  if opts["use_linear_metrics"]
#   if mesh.coord_order == 1
    @assert mesh.coord_order == 1
    getCoordinates(mesh, sbp)  # store coordinates of all nodes into array
    getMetrics(mesh, sbp)

    if mesh.isInterpolated
      interpolateCoordinatesAndMetrics(mesh)
    end

    getFaceNormals(mesh, sbp)
  else  # curvilinear
    # 3rd order and above are not supported yet
    @assert mesh.coord_order <= 2

    getMeshCoordinates(mesh, sbp)
    getFaceCoordinatesAndNormals(mesh, sbp)
    getCurvilinearCoordinatesAndMetrics(mesh, sbp)

  end

  # make sure the mapping jacobian is > 0
  checkMapping(mesh)

  # calculate things that depend on the above
  mesh.min_el_size = getMinElementSize(mesh)
  mesh.volume = calcVolumeIntegral(mesh, sbp)

  # send metric information in parallel
  exchangeMetricInfo(mesh, sbp)

  return nothing
end

"""
  This function recalculates the volume node coordinates and metrics from
  the vertex coordinates in mesh.vert_coords.  This function populates all
  the same fields as [`getAllCoordinatesAndMetics`](@ref) except for 
  `mesh.vert_coords`.

  This function does not support opts["use_linear_metrics"].

  This function is collective on mesh.comm

  **Inputs**

   * mesh: the mesh
   * sbp: the SBP operator
   * opts: options dictionary
"""
function recalcCoordinatesAndMetrics(mesh, sbp, opts)

  @assert !opts["use_linear_metrics"]
  @assert mesh.coord_order <= 2

  getFaceCoordinatesAndNormals(mesh, sbp)
  getCurvilinearCoordinatesAndMetrics(mesh, sbp)

  # make sure the mapping jacobian is > 0
  checkMapping(mesh)

  # calculate things that depend on the above
  mesh.min_el_size = getMinElementSize(mesh)
  mesh.volume = calcVolumeIntegral(mesh, sbp)

  # send metric information in parallel
  exchangeMetricInfo(mesh, sbp)

  return nothing
end


"""
  Reverse mode of getAllCoordinateAndMetrics.  Back propigates
  mesh.nrm_*_bar, mesh.dxidx_bar, mesh.jac_bar to mesh.vert_coords_bar.

  Users should generally call zeroBarArrays() before calling this function
  a second time, to zero out all intermediate arrays.

  This function also recalculates the face normal vectors and coordinates,
  overwriting the relevent arrays in the mesh object.

  This function does not yet support the reverse mode of mesh.remote_metrics.
"""
function getAllCoordinatesAndMetrics_rev(mesh, sbp, opts)

  @assert mesh.commsize == 1
  if opts["use_linear_metrics"]
#   if mesh.coord_order == 1
    @assert mesh.coord_order == 1
    interpolateMapping_rev(mesh)
    getVertCoords_rev(mesh, sbp)
  else
    @assert mesh.coord_order <= 2

    getCurvilinearCoordinatesAndMetrics_rev(mesh, sbp)
    getFaceCoordinateAndNormals_rev(mesh, sbp)  
  end

  return nothing
end

#------------------------------------------------------------------------------
# functions common to linear and curvilinear
"""
  This function allocates the arrays of normal vector for interfaces, boundary
  faces, and shared faces.

  For interfaces, the normal is calculated for elementL
"""
function allocateNormals(mesh::PumiMeshDG{Tmsh}, sbp) where Tmsh

  dim = mesh.dim
  numfacenodes = mesh.numNodesPerFace

  if !isFieldDefined(mesh, :nrm_bndry, :nrm_face, :nrm_sharedface)
    mesh.nrm_bndry = Array(Tmsh, dim, numfacenodes, mesh.numBoundaryFaces )
    mesh.nrm_face = Array(Tmsh, mesh.dim, numfacenodes, mesh.numInterfaces)
    mesh.nrm_sharedface = Array(Array{Tmsh, 3}, mesh.npeers)

    for i=1:mesh.npeers
      mesh.nrm_sharedface[i] = Array(Tmsh, mesh.dim, numfacenodes, mesh.peer_face_counts[i])
    end

    # adjoint parts
    mesh.nrm_bndry_bar = zeros(mesh.nrm_bndry)
    mesh.nrm_face_bar = zeros(mesh.nrm_face)
    mesh.nrm_sharedface_bar = Array(Array{Tmsh, 3}, mesh.npeers)
    for i=1:mesh.npeers
      mesh.nrm_sharedface_bar[i] = zeros(mesh.nrm_sharedface[i])
    end

  else
    fill!(mesh.nrm_bndry, 0.0)
    fill!(mesh.nrm_face, 0.0)
    for i=1:mesh.npeers
      fill!(mesh.nrm_sharedface[i], 0.0)
    end
#=
    fill!(mesh.nrm_bndry_bar, 0.0)
    fill!(mesh.nrm_face_bar, 0.0)
    for i=1:mesh.npeers
      fill!(mesh.nrm_sharedface[i], 0.0)
    end
=#
  end

  return nothing
end

"""
  Send and recieve the metric information for elements on parallel boundaries.
  This function is collective

  **Inputs**

   * mesh: a DG mesh. mesh.remote_metrics is populated by this function
"""
function exchangeMetricInfo(mesh::PumiMeshDG{Tmsh}, sbp) where Tmsh

  remote_metrics = Array(RemoteMetrics{Tmsh}, mesh.npeers)  # receive buffers
  local_metrics = Array(RemoteMetrics{Tmsh}, mesh.npeers)  # send buffers
  send_reqs = Array(MPI.Request, 4, mesh.npeers)
  recv_reqs = Array(MPI.Request, 4, mesh.npeers)


  # allocate the arrays and post the MPI receives
  for i=1:mesh.npeers
    remote_metrics[i] = RemoteMetrics(mesh, i, islocal=false)
    local_metrics[i] = RemoteMetrics(mesh, i, islocal=true)

    obj = remote_metrics[i]
    peernum = mesh.peer_parts[i]
    recv_reqs[1, i] = MPI.Irecv!(obj.vert_coords, peernum, 1, mesh.comm)
    recv_reqs[2, i] = MPI.Irecv!(obj.coords, peernum, 2, mesh.comm)
    recv_reqs[3, i] = MPI.Irecv!(obj.dxidx, peernum, 3, mesh.comm)
    recv_reqs[4, i] = MPI.Irecv!(obj.jac, peernum, 4, mesh.comm)
  end


  # extract the local values and send them
  for i=1:mesh.npeers
    obj = local_metrics[i]
    peernum = mesh.peer_parts[i]
    getLocalMetrics(mesh, obj)

    send_reqs[1, i] = MPI.Isend(obj.vert_coords, peernum, 1, mesh.comm)
    send_reqs[2, i] = MPI.Isend(obj.coords, peernum, 2, mesh.comm)
    send_reqs[3, i] = MPI.Isend(obj.dxidx, peernum, 3, mesh.comm)
    send_reqs[4, i] = MPI.Isend(obj.jac, peernum, 4, mesh.comm)
  end

  # wait for all communication to finish
  # this is not strictly necessary, the real requirement is that the
  # communication finishes before the mesh constructor exits
  for i=1:length(send_reqs)
    MPI.Wait!(send_reqs[i])
    MPI.Wait!(recv_reqs[i])
  end

  mesh.remote_metrics = remote_metrics

  return nothing
end

"""
  Copy the metric values from the main mesh arrays into the RemoteMetrics
  arrays.

  **Inputs**

   * mesh: a mesh object
  
  **Inputs/Outputs**

   * obj: a RemoteMetrics object
"""
function getLocalMetrics(mesh, obj::RemoteMetrics)

  @assert obj.islocal

  numEl = size(obj.vert_coords, 3)
  elnums = mesh.local_element_lists[obj.peer_idx]
  for i=1:numEl
    elnum = elnums[i]

    # vert_coords
    for j=1:size(obj.vert_coords, 2)
      for k=1:size(obj.vert_coords, 1)
        obj.vert_coords[k, j, i] = mesh.vert_coords[k, j, elnum]
      end
    end

    # coords
    for j=1:size(obj.coords, 2)
      for k=1:size(obj.coords, 1)
        obj.coords[k, j, i] = mesh.coords[k, j, elnum]
      end
    end

    # dxidx 
    for j=1:mesh.numNodesPerElement
      for d2=1:mesh.dim
        for d1=1:mesh.dim
          obj.dxidx[d1, d2, j, i] = mesh.dxidx[d1, d2, j, elnum]
        end
      end
    end

    # jac
    for j=1:mesh.numNodesPerElement
      obj.jac[j, i] = mesh.jac[j, elnum]
    end

  end  # end loop i

  return nothing
end



