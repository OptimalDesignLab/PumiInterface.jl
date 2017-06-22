# file for functions related to interpolating from solution points to 
# faces

"""
  Holds all the temporary arrays needed to interpolate the metrics
"""
type Interpolation{Tmsh, Tdim}
  dxdxi_el::Array{Tmsh, 3}
  dxdxi_el_bar::Array{Tmsh, 3}

  dxdxi_elface::Array{Tmsh, 3}
  dxdxi_elface_bar::Array{Tmsh, 3}

  dxidx_node::Array{Tmsh, 2}
  dxidx_node_bar::Array{Tmsh, 2}

  dxdxi_node::Array{Tmsh, 2}
  dxdxi_node_bar::Array{Tmsh, 2}

  bndry_arr::Array{Boundary, 1}
end

function Interpolation{Tmsh}(mesh::PumiMesh2D{Tmsh})
  dxdxi_el = zeros(Tmsh, 4, mesh.numNodesPerElement, 1)
  dxdxi_el_bar = zeros(dxdxi_el)

  dxdxi_elface = zeros(Tmsh, 4, mesh.sbpface.numnodes, 1)
  dxdxi_elface_bar = zeros(dxdxi_elface)

  dxidx_node = zeros(Tmsh, 2, 2)
  dxidx_node_bar = zeros(dxidx_node)

  dxdxi_node = zeros(Tmsh, 2, 2)
  dxdxi_node_bar = zeros(dxdxi_node)

  bndry_arr = Array(Boundary, 1)
  return Interpolation{Tmsh, 2}(dxdxi_el, dxdxi_el_bar,
                                dxdxi_elface, dxdxi_elface_bar,
                                dxidx_node, dxidx_node_bar,
                                dxdxi_node, dxdxi_node_bar,
                                bndry_arr)
end

function Interpolation{Tmsh}(mesh::PumiMesh2D, ::Type{Tmsh})
  dxdxi_el = zeros(Tmsh, 4, mesh.numNodesPerElement, 1)
  dxdxi_el_bar = zeros(dxdxi_el)

  dxdxi_elface = zeros(Tmsh, 4, mesh.sbpface.numnodes, 1)
  dxdxi_elface_bar = zeros(dxdxi_elface)

  dxidx_node = zeros(Tmsh, 2, 2)
  dxidx_node_bar = zeros(dxidx_node)

  dxdxi_node = zeros(Tmsh, 2, 2)
  dxdxi_node_bar = zeros(dxdxi_node)

  bndry_arr = Array(Boundary, 1)
  return Interpolation{Tmsh, 2}(dxdxi_el, dxdxi_el_bar,
                                dxdxi_elface, dxdxi_elface_bar,
                                dxidx_node, dxidx_node_bar,
                                dxdxi_node, dxdxi_node_bar,
                                bndry_arr)
end


function Interpolation{Tmsh}(mesh::PumiMesh3D{Tmsh})
  dxdxi_el = zeros(Tmsh, 9, mesh.numNodesPerElement, 1)
  dxdxi_el_bar = zeros(dxdxi_el)

  dxdxi_elface = zeros(Tmsh, 9, mesh.sbpface.numnodes, 1)
  dxdxi_elface_bar = zeros(dxdxi_elface)

  dxidx_node = zeros(Tmsh, 3, 3)
  dxidx_node_bar = zeros(dxidx_node)

  dxdxi_node = zeros(Tmsh, 3, 3)
  dxdxi_node_bar = zeros(dxdxi_node)

  bndry_arr = Array(Boundary, 1)
  return Interpolation{Tmsh, 3}(dxdxi_el, dxdxi_el_bar,
                                dxdxi_elface, dxdxi_elface_bar,
                                dxidx_node, dxidx_node_bar,
                                dxdxi_node, dxdxi_node_bar,
                                bndry_arr)
end

function Interpolation{Tmsh}(mesh::PumiMesh3D, ::Type{Tmsh})
  dxdxi_el = zeros(Tmsh, 9, mesh.numNodesPerElement, 1)
  dxdxi_el_bar = zeros(dxdxi_el)

  dxdxi_elface = zeros(Tmsh, 9, mesh.sbpface.numnodes, 1)
  dxdxi_elface_bar = zeros(dxdxi_elface)

  dxidx_node = zeros(Tmsh, 3, 3)
  dxidx_node_bar = zeros(dxidx_node)

  dxdxi_node = zeros(Tmsh, 3, 3)
  dxdxi_node_bar = zeros(dxdxi_node)

  bndry_arr = Array(Boundary, 1)
  return Interpolation{Tmsh, 3}(dxdxi_el, dxdxi_el_bar,
                                dxdxi_elface, dxdxi_elface_bar,
                                dxidx_node, dxidx_node_bar,
                                dxdxi_node, dxdxi_node_bar,
                                bndry_arr)
end


"""
  Interpolates dxidx and jac to the face nodes, faces shared in parallel and
  boundary faces, as well as gets the coordinates of the nodes for
  the shared faces and boundaryfaces.  All the functions use smart allocators,
  so the fields of the mesh object will be allocated if needed, but not if 
  they were previously allocated.

  This is the main entry point for everything in this file.
"""
function interpolateCoordinatesAndMetrics(mesh::PumiMeshDG)

  # interpolate metric terms
  interpolateMapping(mesh)

  #  allocate coordinate arrays if needed
  allocateFaceCoordinates(mesh)

  # get the coordinates
  getBndryCoordinates(mesh, mesh.bndryfaces, mesh.coords_bndry)
  getInterfaceCoordinates(mesh, mesh.interfaces, mesh.coords_interface)
  for i=1:mesh.npeers
    getBndryCoordinates(mesh, mesh.bndries_local[i], 
                        mesh.coords_sharedface[i])
  end



  return nothing
end

"""
  Allocates the dxidx and jac arrays for the face nodes/ shared face
  data.  Non-curvilinear meshes only.  Specifically:
    mesh.dxidx_face, mesh.jac_face
    mesh.dxidx_bndry, mesh.jac_bndry
    mesh.dxidx_sharedface, mesh.jac_shareface (and their internal arrays)

    Also the corresponding bar variables
      mesh.dxidx_face_bar, mesh.jac_face_bar
      mesh.dxidx_bndry_bar, mesh.jac_bndry_bar
      mesh.dxidx_sharedface_bar, mesh.jac_sharedface_bar
"""
function allocateInterpolatedMetrics{Tmsh}(mesh::PumiMeshDG{Tmsh})

  dim = mesh.dim
  sbpface = mesh.sbpface
   if !isFieldDefined(mesh, :dxidx_face, :jac_face, :dxidx_sharedface, 
                           :jac_sharedface, :dxidx_bndry, :jac_bndry)
    # interior faces
    mesh.dxidx_face = zeros(Tmsh, dim, dim, sbpface.numnodes, mesh.numInterfaces)
    mesh.jac_face = zeros(Tmsh, sbpface.numnodes, mesh.numInterfaces)

    mesh.dxidx_face_bar = zeros(mesh.dxidx_face)
    mesh.jac_face_bar = zeros(mesh.jac_face)

    # boundary faces
    mesh.dxidx_bndry = zeros(Tmsh, dim, dim, sbpface.numnodes, 
                                   mesh.numBoundaryFaces)
    mesh.jac_bndry = zeros(Tmsh, sbpface.numnodes, mesh.numBoundaryFaces)
    mesh.dxidx_bndry_bar = zeros(mesh.dxidx_bndry)
    mesh.jac_bndry_bar = zeros(mesh.jac_bndry)

    # parallel shared faces
    mesh.dxidx_sharedface = Array(Array{Tmsh, 4}, mesh.npeers)
    mesh.jac_sharedface = Array(Array{Tmsh, 2}, mesh.npeers)
    mesh.dxidx_sharedface_bar = Array(Array{Tmsh, 4}, mesh.npeers)
    mesh.jac_sharedface_bar = Array(Array{Tmsh, 2}, mesh.npeers)
    for i=1:mesh.npeers
      mesh.dxidx_sharedface[i] = zeros(Tmsh, dim, dim, sbpface.numnodes, 
                                       mesh.peer_face_counts[i])

      mesh.jac_sharedface[i] = Array(Tmsh, sbpface.numnodes, mesh.peer_face_counts[i])
      mesh.dxidx_sharedface_bar[i] = zeros(mesh.dxidx_sharedface[i])
      mesh.jac_sharedface_bar[i] = zeros(mesh.jac_sharedface[i])

    end
  else
    # zero out metrics?
    fill!(mesh.dxidx_face, 0.0)
    fill!(mesh.jac_face, 0.0)
    fill!(mesh.dxidx_bndry, 0.0)
    fill!(mesh.jac_bndry, 0.0)
    for i=1:mesh.npeers
      fill!(mesh.dxidx_sharedface[i], 0.0)
      fill!(mesh.jac_sharedface[i], 0.0)
    end
  end


  return nothing
end

"""
  Allocates coordinate arays for the boundary and shared face nodes.
  Specifically:
    mesh.coords_bndry
    mesh.coords_sharedface (and its internal arrays)
"""
function allocateFaceCoordinates{Tmsh}(mesh::PumiMeshDG{Tmsh})

  sbpface = mesh.sbpface

  if !isFieldDefined(mesh, :coords_bndry, :coords_sharedface)
    mesh.coords_bndry = zeros(Tmsh, mesh.dim, sbpface.numnodes, 
                                  mesh.numBoundaryFaces)
    mesh.coords_interface = zeros(Tmsh, mesh.dim, sbpface.numnodes, mesh.numInterfaces)
    mesh.coords_sharedface = Array(Array{Tmsh, 3}, mesh.npeers)
    for i=1:mesh.npeers
      mesh.coords_sharedface[i] = zeros(Tmsh, mesh.dim, sbpface.numnodes, 
                                             mesh.peer_face_counts[i])
    end
  else
    fill!(mesh.coords_bndry, 0.0)
    fill!(mesh.coords_interface, 0.0)
    for i=1:mesh.npeers
      fill!(mesh.coords_sharedface[i], 0.0)
    end

  end


  return nothing
end

"""
  Interpolates dxidx, jac to the face nodes of local interfaces, shared 
  faces, and boundary faces.  This function uses a smart allocator to
  allocate the relevent fields of the mesh object if they have been been
  allocated previously and to not allocate them if they have been previously.

  By convention, the metrics of elementL of an interface are interpolated
  (recall elementL is always the local elements for shared faces).

  See allocateInterpolatedMetrics for the fields that are populated by this
  function.
"""
function interpolateMapping{Tmsh}(mesh::PumiMeshDG{Tmsh})

  allocateInterpolatedMetrics(mesh)

  sbpface = mesh.sbpface
  dim = mesh.dim

  # pull out the needed arrays
  dxidx_face = mesh.dxidx_face
  dxidx_sharedface = mesh.dxidx_sharedface
  jac_face = mesh.jac_face

  dxidx_bndry = mesh.dxidx_bndry
  jac_bndry = mesh.jac_bndry
  jac_sharedface = mesh.jac_sharedface

  # temporary storage
  interp_data = Interpolation(mesh)

  for i=1:mesh.numInterfaces
    dxidx_i = sview(dxidx_face, :, :, :, i)
    jac_i = sview(jac_face, :, i)

    interface_i = mesh.interfaces[i]
    el = interface_i.elementL
    face = interface_i.faceL
    interp_data.bndry_arr[1] =  Boundary(1, face)

    dxidx_in = sview(mesh.dxidx, :, :, :, el)
    jac_in = sview(mesh.jac, :, el)

    interpolateFace(interp_data, mesh.sbpface, dxidx_in, jac_in, dxidx_i, jac_i)

  end
  # now do shared edges
  for peer = 1:mesh.npeers
    dxidx_p = dxidx_sharedface[peer]
    jac_p = jac_sharedface[peer]
    interfaces_p = mesh.shared_interfaces[peer]
    for i=1:mesh.peer_face_counts[peer]
      dxidx_i = sview(dxidx_p, :, :, :, i)
      jac_i = sview(jac_p, :, i)

      interface_i = interfaces_p[i]
      el = interface_i.elementL
      face = interface_i.faceL

      interp_data.bndry_arr[1] = Boundary(1, face)


      dxidx_in = sview(mesh.dxidx, :, :, :, el)
      jac_in = sview(mesh.jac, :, el)

      interpolateFace(interp_data, mesh.sbpface, dxidx_in, jac_in, dxidx_i, jac_i)
    end
  end

  # now do boundary
  for i=1:mesh.numBoundaryFaces
    bndry = mesh.bndryfaces[i]

    dxidx_i = sview(dxidx_bndry, :, :, :, i)
    jac_i = sview(jac_bndry, :, i)

    el = bndry.element
    dxidx_in = sview(mesh.dxidx, :, :, :, el)
    jac_in = sview(mesh.jac, :, el)
    interp_data.bndry_arr[1] = bndry
    interpolateFace(interp_data, mesh.sbpface, dxidx_in, jac_in, dxidx_i, jac_i)
  end

  return nothing

end  # end function

"""
  Reverse mode of interpolateMapping.  This function propigates the adjoint
  parts of the metrics at the interfaces, boundary faces, and shared faces
  back to the volume nodes.

  This function back propigates both dxidx and jac, although for 2D it turns
  out jac has no effect on the computation.

"""
function interpolateMapping_rev{Tmsh}(mesh::PumiMeshDG{Tmsh})

  sbpface = mesh.sbpface
  dim = mesh.dim

  # pull out the needed arrays
  dxidx_face = mesh.dxidx_face
  dxidx_face_bar = mesh.dxidx_face_bar

  dxidx_sharedface = mesh.dxidx_sharedface
  dxidx_sharedface_bar = mesh.dxidx_sharedface_bar

  dxidx_bndry = mesh.dxidx_bndry
  dxidx_bndry_bar = mesh.dxidx_bndry_bar

  jac_face = mesh.jac_face
  jac_face_bar = mesh.jac_face_bar


  jac_bndry = mesh.jac_bndry
  jac_bndry_bar = mesh.jac_bndry_bar

  jac_sharedface = mesh.jac_sharedface
  jac_sharedface_bar = mesh.jac_sharedface_bar

  # temporary storage
  interp_data = Interpolation(mesh)

  for i=1:mesh.numInterfaces
    dxidx_i = sview(dxidx_face, :, :, :, i)
    dxidx_bar_i = sview(dxidx_face_bar, :, :, :, i)
    jac_i = sview(jac_face, :, i)
    jac_bar_i = sview(jac_face_bar, :, i)

    interface_i = mesh.interfaces[i]
    el = interface_i.elementL
    face = interface_i.faceL
    interp_data.bndry_arr[1] =  Boundary(1, face)

    dxidx_in = sview(mesh.dxidx, :, :, :, el)
    dxidx_bar_in = sview(mesh.dxidx_bar, :, :, :, el)
    jac_in = sview(mesh.jac, :, el)
    jac_bar_in = sview(mesh.jac_bar, :, el)

    interpolateFace_rev(interp_data, mesh.sbpface, 
                         dxidx_in, dxidx_bar_in, jac_in, jac_bar_in,
                         dxidx_i, dxidx_bar_i, jac_i, jac_bar_i)

  end
  # now do shared edges
  for peer = 1:mesh.npeers
    dxidx_p = dxidx_sharedface[peer]
    dxidx_bar_p = dxidx_sharedface_bar[peer]
    jac_p = jac_sharedface[peer]
    jac_bar_p = jac_sharedface_bar[perr]

    interfaces_p = mesh.shared_interfaces[peer]
    for i=1:mesh.peer_face_counts[peer]
      dxidx_i = sview(dxidx_p, :, :, :, i)
      dxidx_bar_i = sview(dxidx_bar_p, :, :, :, i)
      jac_i = sview(jac_p, :, i)
      jac_bar_i = sview(jac_bar_p, :, i)

      interface_i = interfaces_p[i]
      el = interface_i.elementL
      face = interface_i.faceL

      interp_data.bndry_arr[1] = Boundary(1, face)


      dxidx_in = sview(mesh.dxidx, :, :, :, el)
      dxidx_bar_in = sview(mesh.dxidx_bar, :, :, :, el)
      jac_in = sview(mesh.jac, :, el)
      jac_bar_in = sview(mesh.jac_bar, :, el)

      interpolateFace_rev(interp_data, mesh.sbpface,
                      dxidx_in, dxidx_bar_in, jac_in, jac_bar_in,
                      dxidx_i, dxidx_bar_i, jac_i, jac_bar_i)
    end
  end

  # now do boundary
  for i=1:mesh.numBoundaryFaces
    bndry = mesh.bndryfaces[i]

    dxidx_i = sview(dxidx_bndry, :, :, :, i)
    dxidx_bar_i = sview(dxidx_bndry_bar, :, :, :, i)
    jac_i = sview(jac_bndry, :, i)
    jac_bar_i = sview(jac_bndry_bar, :, i)

    el = bndry.element

    dxidx_in = sview(mesh.dxidx, :, :, :, el)
    dxidx_bar_in = sview(mesh.dxidx_bar, :, :, :, el)
    jac_in = sview(mesh.jac, :, el)
    jac_bar_in = sview(mesh.jac_bar, :, el)

    interp_data.bndry_arr[1] = bndry
    interpolateFace_rev(interp_data, mesh.sbpface, 
                    dxidx_in, dxidx_bar_in, jac_in, jac_bar_in,
                    dxidx_i, dxidx_bar_i, jac_i, jac_bar_i)
  end


  return nothing
end


# 2D version
function interpolateFace{Tmsh}(interp_data::Interpolation{Tmsh, 2}, sbpface, 
                               dxidx_hat_in::AbstractArray{Tmsh, 3}, 
                               jac_in::AbstractVector{Tmsh}, 
                               dxidx_hat_out::AbstractArray{Tmsh, 3}, 
                               jac_out::AbstractVector{Tmsh})

  # unpack argumetns
  dxdxi_el = interp_data.dxdxi_el
  dxdxi_elface = interp_data.dxdxi_elface
  dxdxi_node = interp_data.dxdxi_node
  dxidx_node = interp_data.dxidx_node

  dim = 2
  numNodesPerElement = size(dxidx_hat_in, 3)
  # get the data
  for j=1:numNodesPerElement
    dxidx_hat = sview(dxidx_hat_in, :, :, j)
    detJ = jac_in[j]

    adjugate2(dxidx_hat, dxidx_node)
    for k=1:dim
      for p=1:dim
        pos = p + dim*(k-1)
        dxdxi_el[pos,j,1] = dxidx_node[k, p]
      end
    end
  end


  # interpolate to the face
  face = interp_data.bndry_arr[1].face
  interp_data.bndry_arr[1] = Boundary(1, face)
  boundaryinterpolate!(sbpface, interp_data.bndry_arr, dxdxi_el, dxdxi_elface)


  # now store dxidx, |J| at the boundary nodesa
  for j=1:sbpface.numnodes
    for k=1:dim
      for p=1:dim
        pos = p + dim*(k-1)
        dxdxi_node[k,p] = dxdxi_elface[pos, j, 1]
      end
    end

    # inv(A) = adj(A)/|A|
    det_dxdxi = det2(dxdxi_node)
    adjugate2(dxdxi_node, dxidx_node)
    scale!(dxidx_node, 1./det_dxdxi)

    detJ = det2(dxidx_node)
    fac = 1./detJ
    for k=1:dim
      for p=1:dim
        dxidx_hat_out[p, k, j] = dxidx_node[p, k]*fac
      end
    end
    jac_out[j] = detJ
  end

  return nothing

end

"""
  Reverse mode of interpolateFace.

  Propigates the adjoint parts of the inputs to the adjoint parts of the
  outputs.  This function also verifies the primal parts have not been modified
  (ie. it computes values dxidx_hat_out and verifies they are the same as
  the ones already present in the array).

  Inputs (see interpolateFace):
    interp_data
    sbpface
    dxidx_hat_in
    jac_in
    dxidx_hat_out
    dxidx_hat_out_bar: the adjoint part of dxidx_hat_out
    jac_out
    jac_out_bar: the adjoint part of jac_out

  Outputs:
    dxidx_hat_in_bar: the adjoint part of dxidx_hat_in
    jac_in_bar: the adjoint part of jac_in

  Aliasing restrictions: please don't
"""
function interpolateFace_rev{Tmsh}(interp_data::Interpolation{Tmsh, 2}, 
                               sbpface, 
                               dxidx_hat_in::AbstractArray{Tmsh, 3}, 
                               dxidx_hat_in_bar::AbstractArray{Tmsh, 3},
                               jac_in::AbstractVector{Tmsh}, 
                               jac_in_bar::AbstractVector{Tmsh},
                               dxidx_hat_out::AbstractArray{Tmsh, 3}, 
                               dxidx_hat_out_bar::AbstractArray{Tmsh, 3},
                               jac_out::AbstractVector{Tmsh}, 
                               jac_out_bar::AbstractVector{Tmsh})

  #----------------------------------------------------------------------------
  # forward sweep
  # unpack arguments
  dxdxi_el = interp_data.dxdxi_el
  dxdxi_elface = interp_data.dxdxi_elface
  dxdxi_node = interp_data.dxdxi_node
  dxidx_node = interp_data.dxidx_node

  fill!(dxdxi_el, 0.0)
  fill!(dxdxi_elface, 0.0)

  dim = 2
  numNodesPerElement = size(dxidx_hat_in, 3)

  # get the data
  for j=1:numNodesPerElement
    dxidx_hat = sview(dxidx_hat_in, :, :, j)
#    detJ = jac_in[j]

    adjugate2(dxidx_hat, dxidx_node)
    for k=1:dim
      for p=1:dim
        pos = p + dim*(k-1)
        dxdxi_el[pos,j,1] = dxidx_node[k, p]
      end
    end
  end


  # interpolate to the face
  face = interp_data.bndry_arr[1].face
  interp_data.bndry_arr[1] = Boundary(1, face)
  boundaryinterpolate!(sbpface, interp_data.bndry_arr, dxdxi_el, dxdxi_elface)


# everything in here is local to the loop, no need to compute it here
  for j=1:sbpface.numnodes
    for k=1:dim
      for p=1:dim
        pos = p + dim*(k-1)
        dxdxi_node[k,p] = dxdxi_elface[pos, j, 1]
      end
    end

    # inv(A) = adj(A)/|A|
    det_dxdxi = det2(dxdxi_node)
    adjugate2(dxdxi_node, dxidx_node)
    for i=1:length(dxidx_node)
       dxidx_node[i] = dxidx_node[i] / det_dxdxi
    end

    detJ = det2(dxidx_node)
    fac = 1./detJ

    for k=1:dim
      for p=1:dim
        # verify this forward sweep is consistent with the original function
        @assert (dxidx_hat_out[p, k, j] - dxidx_node[p, k]*fac) < 1e-13
#        dxidx_hat_out[p, k, j] = dxidx_node[p, k]*fac
      end
    end
    jac_out[j] = detJ
  end


  #----------------------------------------------------------------------------
  # reverse sweep

  # unpack arguments
  dxdxi_el_bar = interp_data.dxdxi_el_bar
  dxdxi_elface_bar = interp_data.dxdxi_elface_bar
  dxdxi_node_bar = interp_data.dxdxi_node_bar
  dxidx_node_bar = interp_data.dxidx_node_bar

  # zero them out
  fill!(dxdxi_el_bar, 0.0)
  fill!(dxdxi_elface_bar, 0.0)

  # each interation of this loop is independent - combine forward and reverse
  for j=sbpface.numnodes:-1:1

    #--------------------------------------------------------------------------
    # forward sweep
    fill!(dxdxi_node_bar, 0.0)
    fill!(dxidx_node_bar, 0.0)

    for k=1:dim
      for p=1:dim
        pos = p + dim*(k-1)
        dxdxi_node[k,p] = dxdxi_elface[pos, j, 1]
      end
    end

    # inv(A) = adj(A)/|A|
    det_dxdxi = det2(dxdxi_node)
    adjugate2(dxdxi_node, dxidx_node)
    for i=1:length(dxidx_node)
       dxidx_node[i] = dxidx_node[i] / det_dxdxi
    end

    detJ = det2(dxidx_node)
    fac = 1./detJ
#=
    for k=1:dim
      for p=1:dim
        @assert (dxidx_hat_out[p, k, j] - dxidx_node[p, k]*fac) < 1e-13
        dxidx_hat_out[p, k, j] = dxidx_node[p, k]*fac
      end
    end
    jac_out[j] = detJ
=#
    #--------------------------------------------------------------------------
    # reverse sweep

    det_dxdxi_bar = zero(det_dxdxi)
    fac_bar = zero(eltype(dxidx_hat_out))

    detJ_bar = jac_out_bar[j]*1

    for k=dim:-1:1
      for p=dim:-1:1
        dxidx_node_bar[p, k] += dxidx_hat_out_bar[p, k, j]*fac

        fac_bar += dxidx_hat_out_bar[p, k, j]*dxidx_node[p, k]
      end
    end

    detJ_bar -= fac_bar*1/(detJ*detJ)
    det2_rev(dxidx_node, dxidx_node_bar, detJ_bar)


    # uncomment this to introdue the bug
    for i=1:length(dxidx_node)
       dxidx_node[i] *= det_dxdxi  # need to reverse the primal variables too
       det_dxdxi_bar += dxidx_node_bar[i]*dxidx_node[i]*(-1/(det_dxdxi*det_dxdxi))
       dxidx_node_bar[i] = dxidx_node_bar[i]/det_dxdxi
    end


    adjugate2_rev(dxdxi_node, dxdxi_node_bar, dxidx_node, dxidx_node_bar)
    det2_rev(dxdxi_node, dxdxi_node_bar, det_dxdxi_bar)


    for k=dim:-1:1
      for p=dim:-1:1
        pos = p + dim*(k-1)
        dxdxi_elface_bar[pos, j, 1] += dxdxi_node_bar[k, p]
      end
    end

  end  # end loop j
  


  # interpolate to the face
  face = interp_data.bndry_arr[1].face
  interp_data.bndry_arr[1] = Boundary(1, face)
  boundaryinterpolate_rev!(sbpface, interp_data.bndry_arr, dxdxi_el_bar, dxdxi_elface_bar)

  # propagate back to the inputs
  for j=numNodesPerElement:-1:1
    fill!(dxidx_node_bar, 0.0)
    fill!(dxdxi_node_bar, 0.0)

    dxidx_hat = sview(dxidx_hat_in, :, :, j)
    dxidx_hat_bar = sview(dxidx_hat_in_bar, :, :, j)

    for k=dim:-1:1
      for p=dim:-1:1
        pos = p + dim*(k-1)
        dxidx_node_bar[k, p] += dxdxi_el_bar[pos, j, 1]
      end
    end

    adjugate2_rev(dxidx_hat, dxidx_hat_bar, dxidx_node, dxidx_node_bar)

    # it turns out jac_in isn't used for anything
#    jac_in_bar[j] = detJ_bar
  end  # end loop j
    
  return nothing

end


# 3d version
"""
  Interpolates the metric terms from the volume nodes of an element to 
  the nodes of a particular face.

  Inputs:
    interp_data:  Interpolation object containing the temporary arrays.
                  Which face is interpolated to depends 
                  interp_data.bndry_arr.face.  The element field does
                  not matter.  This array is overwritten.

    sbpface: an SBPFace
    dxidx_hat_in: 3D array holding dxidx/|J| at the volume nodes,
                   2 x 2 x numNodesPerElement in 2D or 
                   3 x 3 x numNodesPerElemente in 3D.
    jac_hat_in: determanent of the dxidx at the volume nodes, vector of
                length numNodesPerElement

  Outputs:
    dxidx_i: array to be populated with dxidx/|J| at the face nodes,
             2 x 2 x numNodesPerFace in 2 or
             3 x 3 x numNodesPerFace in 3D

    jac_i: array to be populated with |J|, vector of length numNodesPerFAce

  Aliasing restrictions: do not alias
"""
function interpolateFace{Tmsh}(interp_data::Interpolation{Tmsh, 3}, sbpface, 
                               dxidx_hat_in::AbstractArray{Tmsh, 3}, 
                               jac_in::AbstractVector{Tmsh}, 
                               dxidx_hat_out::AbstractArray{Tmsh, 3}, 
                               jac_out::AbstractVector{Tmsh})
# herein, dxidx_hat is dxidx/|J|, and dxidx is the true dxidx

  dxdxi_el = interp_data.dxdxi_el
  dxdxi_elface = interp_data.dxdxi_elface
  dxdxi_node = interp_data.dxdxi_node
  dxidx_node = interp_data.dxidx_node

  dim = 3
  numNodesPerElement = size(dxidx_hat_in, 3)
  # get the data
  for j=1:numNodesPerElement
    dxidx_hat = sview(dxidx_hat_in, :, :, j)
    detJ = jac_in[j]

    adjugate3(dxidx_hat, dxidx_node)
    for k=1:dim
      for p=1:dim
        dxidx_node[k, p] *= detJ
        pos = p + dim*(k-1)
        dxdxi_el[pos,j,1] = dxidx_node[k, p]
      end
    end
  end


  # interpolate to the face
#  face = bndry.face
#  bndry_arr = [Boundary(1, face)]
  face = interp_data.bndry_arr[1].face
  interp_data.bndry_arr[1] = Boundary(1, face)

  boundaryinterpolate!(sbpface, interp_data.bndry_arr, dxdxi_el, dxdxi_elface)


  # now store dxidx, |J| at the boundary nodesa
  for j=1:sbpface.numnodes
    for k=1:dim
      for p=1:dim
        pos = p + dim*(k-1)
        dxdxi_node[k,p] = dxdxi_elface[pos, j, 1]
      end
    end

    # inv(A) = adj(A)/|A|
    det_dxdxi = det3(dxdxi_node)
    adjugate3(dxdxi_node, dxidx_node)
    scale!(dxidx_node, 1./det_dxdxi)

    detJ = det3(dxidx_node)
    fac = 1./detJ
    for k=1:dim
      for p=1:dim
        dxidx_hat_out[p, k, j] = dxidx_node[p, k]*fac
      end
    end
    jac_out[j] = detJ
  end

  return nothing

end

function interpolateFace_rev{Tmsh}(interp_data::Interpolation{Tmsh, 3}, 
                               sbpface, 
                               dxidx_hat_in::AbstractArray{Tmsh, 3}, 
                               dxidx_hat_in_bar::AbstractArray{Tmsh, 3},
                               jac_in::AbstractVector{Tmsh}, 
                               jac_in_bar::AbstractVector{Tmsh},
                               dxidx_hat_out::AbstractArray{Tmsh, 3}, 
                               dxidx_hat_out_bar::AbstractArray{Tmsh, 3},
                               jac_out::AbstractVector{Tmsh}, 
                               jac_out_bar::AbstractVector{Tmsh})

  #----------------------------------------------------------------------------
  # forward sweep
  # unpack arguments
  dxdxi_el = interp_data.dxdxi_el
  dxdxi_elface = interp_data.dxdxi_elface
  dxdxi_node = interp_data.dxdxi_node
  dxidx_node = interp_data.dxidx_node

  fill!(dxdxi_el, 0.0)
  fill!(dxdxi_elface, 0.0)

  dim = 3
  numNodesPerElement = size(dxidx_hat_in, 3)

  # get the data
  for j=1:numNodesPerElement
    dxidx_hat = sview(dxidx_hat_in, :, :, j)
    detJ = jac_in[j]

    adjugate3(dxidx_hat, dxidx_node)
    for k=1:dim
      for p=1:dim
        dxidx_node[k, p] *= detJ
        pos = p + dim*(k-1)
        dxdxi_el[pos,j,1] = dxidx_node[k, p]
      end
    end
  end


  # interpolate to the face
  face = interp_data.bndry_arr[1].face
  interp_data.bndry_arr[1] = Boundary(1, face)
  boundaryinterpolate!(sbpface, interp_data.bndry_arr, dxdxi_el, dxdxi_elface)


  # everything in here is local to the loop, no need to compute it here
  for j=1:sbpface.numnodes
    for k=1:dim
      for p=1:dim
        pos = p + dim*(k-1)
        dxdxi_node[k,p] = dxdxi_elface[pos, j, 1]
      end
    end

    # inv(A) = adj(A)/|A|
    det_dxdxi = det3(dxdxi_node)
    adjugate3(dxdxi_node, dxidx_node)
    for i=1:length(dxidx_node)
       dxidx_node[i] = dxidx_node[i] / det_dxdxi
    end

    detJ = det3(dxidx_node)
    fac = 1./detJ

    for k=1:dim
      for p=1:dim
        # verify this forward sweep is consistent with the original function
        @assert (dxidx_hat_out[p, k, j] - dxidx_node[p, k]*fac) < 1e-13
#        dxidx_hat_out[p, k, j] = dxidx_node[p, k]*fac
      end
    end
#    jac_out[j] = detJ
  end

  #----------------------------------------------------------------------------
  # reverse sweep

  # unpack arguments
  dxdxi_el_bar = interp_data.dxdxi_el_bar
  dxdxi_elface_bar = interp_data.dxdxi_elface_bar
  dxdxi_node_bar = interp_data.dxdxi_node_bar
  dxidx_node_bar = interp_data.dxidx_node_bar

  # zero them out
  fill!(dxdxi_el_bar, 0.0)
  fill!(dxdxi_elface_bar, 0.0)

  # each interation of this loop is independent - combine forward and reverse
  for j=sbpface.numnodes:-1:1

    #--------------------------------------------------------------------------
    # forward sweep
    fill!(dxdxi_node_bar, 0.0)
    fill!(dxidx_node_bar, 0.0)

    for k=1:dim
      for p=1:dim
        pos = p + dim*(k-1)
        dxdxi_node[k,p] = dxdxi_elface[pos, j, 1]
      end
    end

    # inv(A) = adj(A)/|A|
    det_dxdxi = det3(dxdxi_node)
    adjugate3(dxdxi_node, dxidx_node)
    for i=1:length(dxidx_node)
       dxidx_node[i] = dxidx_node[i] / det_dxdxi
    end

    detJ = det3(dxidx_node)
    fac = 1./detJ
#=
    for k=1:dim
      for p=1:dim
        @assert (dxidx_hat_out[p, k, j] - dxidx_node[p, k]*fac) < 1e-13
        dxidx_hat_out[p, k, j] = dxidx_node[p, k]*fac
      end
    end
    jac_out[j] = detJ
=#
    #--------------------------------------------------------------------------
    # reverse sweep
    det_dxdxi_bar = zero(det_dxdxi)
    fac_bar = zero(eltype(dxidx_hat_out))

    detJ_bar = jac_out_bar[j]*1

    for k=dim:-1:1
      for p=dim:-1:1
        dxidx_node_bar[p, k] += dxidx_hat_out_bar[p, k, j]*fac

        fac_bar += dxidx_hat_out_bar[p, k, j]*dxidx_node[p, k]
      end
    end

    detJ_bar -= fac_bar*1/(detJ*detJ)
    det3_rev(dxidx_node, dxidx_node_bar, detJ_bar)


    # uncomment this to introdue the bug
    for i=1:length(dxidx_node)
       dxidx_node[i] *= det_dxdxi  # need to reverse the primal variables too
       det_dxdxi_bar += dxidx_node_bar[i]*dxidx_node[i]*(-1/(det_dxdxi*det_dxdxi))
       dxidx_node_bar[i] = dxidx_node_bar[i]/det_dxdxi
    end


    adjugate3_rev(dxdxi_node, dxdxi_node_bar, dxidx_node, dxidx_node_bar)
    det3_rev(dxdxi_node, dxdxi_node_bar, det_dxdxi_bar)


    for k=dim:-1:1
      for p=dim:-1:1
        pos = p + dim*(k-1)
        dxdxi_elface_bar[pos, j, 1] += dxdxi_node_bar[k, p]
      end
    end

  end  # end loop j
  

  # interpolate to the face
  face = interp_data.bndry_arr[1].face
  interp_data.bndry_arr[1] = Boundary(1, face)
  boundaryinterpolate_rev!(sbpface, interp_data.bndry_arr, dxdxi_el_bar, dxdxi_elface_bar)

  # propagate back to the inputs
  for j=numNodesPerElement:-1:1
    fill!(dxidx_node_bar, 0.0)
    fill!(dxdxi_node_bar, 0.0)
    detJ = jac_in[j]
    detJ_bar = zero(detJ)

    dxidx_hat = sview(dxidx_hat_in, :, :, j)
    dxidx_hat_bar = sview(dxidx_hat_in_bar, :, :, j)

    for k=dim:-1:1
      for p=dim:-1:1
        pos = p + dim*(k-1)
        dxidx_node_bar[k, p] += dxdxi_el_bar[pos, j, 1]
        
        dxidx_node[k, p] /= detJ
        detJ_bar += dxidx_node_bar[k, p]*dxidx_node[k, p]
        dxidx_node_bar[k, p] = dxidx_node_bar[k, p]*detJ
      end
    end

    adjugate3_rev(dxidx_hat, dxidx_hat_bar, dxidx_node, dxidx_node_bar)

    jac_in_bar[j] = detJ_bar
  end  # end loop j
    
  return nothing

end


#TODO: reorder these operations for column major arrays
function det2(A::AbstractMatrix)
# determinet of 2x2 matrix

  return A[1,1]*A[2,2] - A[1,2]*A[2,1]
end

"""
  Computes the reverse mode of det2, returning only the derivative and not
  the value of the determinant.
  
  Inputs:
    A: The matrix
    y_bar: the adjoint part of the determinant (the primal part is not required
           for the reverse calculation)

  Inputs/Outputs:
    A_bar: the adjoint part of A, which is incremented, not overwritten

"""
function det2_rev(A::AbstractMatrix, A_bar::AbstractMatrix, y_bar::Number)

  A_bar[1,1] += y_bar*A[2,2]
  A_bar[1,2] += -y_bar*A[2,1]
  A_bar[2,1] += -y_bar*A[1,2]
  A_bar[2,2] += y_bar*A[1,1]
  
  return nothing
end


function det3(A::AbstractMatrix)
# determinent of 3x3 matrix
# cofactor expansion about first column
  val1 = A[1,1]*(A[2,2]*A[3,3] - A[2,3]*A[3,2])
  val2 = -A[2,1]*(A[1,2]*A[3,3] - A[1,3]*A[3,2])
  val3 = A[3,1]*(A[1,2]*A[2,3] - A[1,3]*A[2,2])
  return val1 + val2 + val3
end

function det3_rev(A::AbstractMatrix, A_bar::AbstractMatrix, y_bar::Number)
  
  # initialize intermediate adjoints to zero
  val1_bar = y_bar
  val2_bar = y_bar
  val3_bar = y_bar

  # val3 expression
  # for true alias compatability, need to apply reverse mode to each line 
  # right to left, rather than left to right
  A_bar[1,1] += val3_bar*(A[2,2]*A[3,3] - A[2,3]*A[3,2])
  A_bar[2,2] += val3_bar*(A[1,1]*A[3,3])
  A_bar[3,3] += val3_bar*(A[1,1]*A[2,2])

  A_bar[2,3] += -val3_bar*A[1,1]*A[3,2]
  A_bar[3,2] += -val3_bar*A[1,1]*A[2,3]

  # val2 expression
  A_bar[2,1] += -val2_bar*(A[1,2]*A[3,3] - A[1,3]*A[3,2])
  A_bar[1,2] += -val2_bar*A[2,1]*A[3,3]
  A_bar[3,3] += -val2_bar*A[2,1]*A[1,2]

  A_bar[1,3] += val2_bar*A[2,1]*A[3,2]
  A_bar[3,2] += val2_bar*A[2,1]*A[1,3]

  # val1 expression
  A_bar[3,1] += val1_bar*(A[1,2]*A[2,3] - A[1,3]*A[2,2])
  A_bar[1,2] += val1_bar*A[3,1]*A[2,3]
  A_bar[2,3] += val1_bar*A[3,1]*A[1,2]

  A_bar[1,3] += -val1_bar*A[3,1]*A[2,2]
  A_bar[2,2] += -val1_bar*A[3,1]*A[1,3]

  return nothing
end



function adjugate2(A::AbstractMatrix, B::AbstractMatrix)
  B[1,1] = A[2,2]
  B[1,2] = -A[1,2]
  B[2,1] = -A[2,1]
  B[2,2] = A[1,1]
end

"""
  Computes the reverse mode of adjugate2, without modifying A or B

  Inputs:
    A: The original matrix
    B: The adjuage matrix (presumably computed by adjugate2)
    B_bar: the adjoint part of B

  Inputs/Outputs:
    A_bar: The adjoint part of A, which is incremented, not overwritten

  Aliasing restrictions: I'm not sure what happens if things alias, so don't
                         do it.
"""
function adjugate2_rev(A::AbstractMatrix, A_bar::AbstractMatrix, 
                       B::AbstractMatrix, B_bar::AbstractMatrix)

  A_bar[2,2] += B_bar[1,1]
  A_bar[1,2] += -B_bar[1,2]
  A_bar[2,1] += -B_bar[2,1]
  A_bar[1,1] += B_bar[2,2]


end

function adjugate3(A::AbstractMatrix, B::AbstractMatrix)

  B[1,1] = A[2,2]*A[3,3] - A[2,3]*A[3,2]
  B[1,2] = -(A[1,2]*A[3,3] - A[1,3]*A[3,2])
  B[1,3] = A[1,2]*A[2,3] - A[1,3]*A[2,2]

  B[2,1] = -(A[2,1]*A[3,3] - A[2,3]*A[3,1])
  B[2,2] = A[1,1]*A[3,3] - A[1,3]*A[3,1]
  B[2,3] = -(A[1,1]*A[2,3] - A[1,3]*A[2,1])

  B[3,1] = A[2,1]*A[3,2] - A[2,2]*A[3,1]
  B[3,2] = -(A[1,1]*A[3,2] - A[1,2]*A[3,1])
  B[3,3] = A[1,1]*A[2,2] - A[1,2]*A[2,1]
end


"""
  Computes the reverse mode of adjugate2, without modifying A or B
  
  Inputs:
    A: the original matrix
    B: the adjugate matrix (presumably computed by adjugate3)
    B_bar: the adjoint part of B

  Inputs/Outputs
    A_bar: the adjoint part of A, which is incremented, not overwritten

  Aliasing restrictions: Please don't
"""
function adjugate3_rev(A::AbstractMatrix, A_bar::AbstractMatrix,
                       B::AbstractMatrix, B_bar::AbstractMatrix)

  # B11
  B11_bar = B_bar[1,1]
  A_bar[2,2] += B11_bar*A[3,3]
  A_bar[3,3] += B11_bar*A[2,2]
  A_bar[2,3] += -B11_bar*A[3,2]
  A_bar[3,2] += -B11_bar*A[2,3]

  bar = B_bar[1,2]
  A_bar[1,2] += -bar*A[3,3]
  A_bar[3,3] += -bar*A[1,2]
  A_bar[1,3] += bar*A[3,2]
  A_bar[3,2] += bar*A[1,3]

  bar = B_bar[1,3]
  A_bar[1,2] += bar*A[2,3]
  A_bar[2,3] += bar*A[1,2]
  A_bar[1,3] += -bar*A[2,2]
  A_bar[2,2] += -bar*A[1,3]

  bar = B_bar[2,1]
  A_bar[2,1] += -bar*A[3,3]
  A_bar[3,3] += -bar*A[2,1]
  A_bar[2,3] += bar*A[3,1]
  A_bar[3,1] += bar*A[2,3]

  bar = B_bar[2,2]
  A_bar[1,1] += bar*A[3,3]
  A_bar[3,3] += bar*A[1,1]
  A_bar[1,3] += -bar*A[3,1]
  A_bar[3,1] += -bar*A[1,3]

  bar = B_bar[2,3]
  A_bar[1,1] += -bar*A[2,3]
  A_bar[2,3] += -bar*A[1,1]
  A_bar[1,3] += bar*A[2,1]
  A_bar[2,1] += bar*A[1,3]

  bar = B_bar[3,1]
  A_bar[2,1] += bar*A[3,2]
  A_bar[3,2] += bar*A[2,1]
  A_bar[2,2] += -bar*A[3,1]
  A_bar[3,1] += -bar*A[2,2]

  bar = B_bar[3,2]
  A_bar[1,1] += -bar*A[3,2]
  A_bar[3,2] += -bar*A[1,1]
  A_bar[1,2] += bar*A[3,1]
  A_bar[3,1] += bar*A[1,2]

  bar = B_bar[3,3]
  A_bar[1,1] += bar*A[2,2]
  A_bar[2,2] += bar*A[1,1]
  A_bar[1,2] += -bar*A[2,1]
  A_bar[2,1] += -bar*A[1,2]

  return nothing
end

