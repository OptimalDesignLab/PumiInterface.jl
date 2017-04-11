# file for functions related to interpolating from solution points to 
# faces

type Interpolation{Tmsh, Tdim}
  dxdxi_el::Array{Tmsh, 3}
  dxdxi_elface::Array{Tmsh, 3}
  dxidx_node::Array{Tmsh, 2}
  dxdxi_node::Array{Tmsh, 2}
  bndry_arr::Array{Boundary, 1}
end

function Interpolation{Tmsh}(mesh::PumiMesh2D{Tmsh})
  dxdxi_el = zeros(Tmsh, 4, mesh.numNodesPerElement, 1)
  dxdxi_elface = zeros(Tmsh, 4, mesh.sbpface.numnodes, 1)
  dxidx_node = zeros(Tmsh, 2, 2)
  dxdxi_node = zeros(Tmsh, 2, 2)
  bndry_arr = Array(Boundary, 1)
  return Interpolation{Tmsh, 2}(dxdxi_el, dxdxi_elface, dxidx_node, dxdxi_node, bndry_arr)
end

function Interpolation{Tmsh}(mesh::PumiMesh3D{Tmsh})
  dxdxi_el = zeros(Tmsh, 9, mesh.numNodesPerElement, 1)
  dxdxi_elface = zeros(Tmsh, 9, mesh.sbpface.numnodes, 1)
  dxidx_node = zeros(Tmsh, 3, 3)
  dxdxi_node = zeros(Tmsh, 3, 3)
  bndry_arr = Array(Boundary, 1)
  return Interpolation{Tmsh, 3}(dxdxi_el, dxdxi_elface, dxidx_node, dxdxi_node, bndry_arr)
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
  if !isFieldDefined(mesh, :coords_bndry, :coords_sharedface)
    allocateFaceCoordinates(mesh)
  else  # zero out arrays that might need it
    fill!(mesh.coords_bndry, 0.0)
    for i=1:mesh.npeers
      fill!(mesh.coords_sharedface[i], 0.0)
    end
  end

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
  data. Specifically:
    mesh.dxidx_face, mesh.jac_face
    mesh.dxidx_bndry, mesh.jac_bndry
    mesh.dxidx_sharedface, mesh.jac_shareface (and their internal arrays)

"""
function allocateInterpolatedMetrics{Tmsh}(mesh::PumiMeshDG{Tmsh})

  dim = mesh.dim
  sbpface = mesh.sbpface
  
  # interior faces
  mesh.dxidx_face = zeros(Tmsh, dim, dim, sbpface.numnodes, mesh.numInterfaces)
  mesh.jac_face = zeros(Tmsh, sbpface.numnodes, mesh.numInterfaces)
  mesh.dxidx_face_bar = zeros(mesh.dxidx_face)

  # boundary faces
  mesh.dxidx_bndry = zeros(Tmsh, dim, dim, sbpface.numnodes, 
                                 mesh.numBoundaryFaces)
  mesh.jac_bndry = zeros(Tmsh, sbpface.numnodes, mesh.numBoundaryFaces)
  mesh.dxidx_bndry = zeros(mesh.dxidx_bndry)

  # parallel shared faces
  mesh.dxidx_sharedface = Array(Array{Tmsh, 4}, mesh.npeers)
  mesh.jac_sharedface = Array(Array{Tmsh, 2}, mesh.npeers)
  mesh.dxidx_sharedface_bar = Array(Array{Tmsh, 4}, mesh.npeers)
  for i=1:mesh.npeers
    mesh.dxidx_sharedface[i] = zeros(Tmsh, dim, dim, sbpface.numnodes, 
                                     mesh.peer_face_counts[i])

    mesh.jac_sharedface[i] = Array(Tmsh, sbpface.numnodes, mesh.peer_face_counts[i])
    mesh.dxidx_sharedface_bar[i] = zeros(mesh.dxidx_sharedface[i])
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
  mesh.coords_bndry = zeros(Tmsh, mesh.dim, sbpface.numnodes, 
                                mesh.numBoundaryFaces)
  mesh.coords_interface = zeros(Tmsh, mesh.dim, sbpface.numnodes, mesh.numInterfaces)
  mesh.coords_sharedface = Array(Array{Tmsh, 3}, mesh.npeers)
  for i=1:mesh.npeers
    mesh.coords_sharedface[i] = zeros(Tmsh, mesh.dim, sbpface.numnodes, 
                                           mesh.peer_face_counts[i])
  end


  return nothing
end

"""
  Interpolates dxidx, jac to the face nodes of local interfaces, shared 
  faces, and boundary faces.  This function uses a smart allocator to
  allocate the relevent fields of the mesh object if they have been been
  allocated previously and to not allocate them if they have been previously.

  See allocateInterpolatedMetrics for the fields that are populated by this
  function.
"""
# interpolates dxidx, jac to the face nodes
# only do this for elementL of each interface
function interpolateMapping{Tmsh}(mesh::PumiMeshDG{Tmsh})

#  if (!isdefined(mesh, :dxidx_face) || !isdefined(mesh, :jac_face) ||
#      !isdefined(mesh, :dxidx_sharedface) || !isdefined(mesh, :jac_sharedface
#      || !isdefined(mesh, :dxidx_bndry) || !isdefined(mesh, :jac_bndry)

  if !isFieldDefined(mesh, :dxidx_face, :jac_face, :dxidx_sharedface, 
                           :jac_sharedface, :dxidx_bndry, :jac_bndry)
    allocateInterpolatedMetrics(mesh)
  else  # zero out any arrays that might need it
    fill!(mesh.dxidx_face, 0.0)
    fill!(mesh.jac_face, 0.0)
    fill!(mesh.dxidx_bndry, 0.0)
    fill!(mesh.jac_bndry, 0.0)
    for i=1:mesh.npeers
      fill!(mesh.dxidx_sharedface[i], 0.0)
    end
  end
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

# 2D version
function interpolateFace{Tmsh}(interp_data::Interpolation{Tmsh, 2}, sbpface, dxidx_hat_in, jac_in, dxidx_i, jac_i)

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
        dxidx_i[p, k, j] = dxidx_node[p, k]*fac
      end
    end
    jac_i[j] = detJ
  end

  return nothing

end

# 3d version
function interpolateFace{Tmsh}(interp_data::Interpolation{Tmsh, 3}, sbpface, dxidx_hat_in, jac_in, dxidx_i, jac_i)

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
        dxidx_i[p, k, j] = dxidx_node[p, k]*fac
      end
    end
    jac_i[j] = detJ
  end

  return nothing

end


function det2(A::AbstractMatrix)
# determinet of 2x2 matrix

  return A[1,1]*A[2,2] - A[1,2]*A[2,1]
end

function det3(A::AbstractMatrix)
# determinent of 3x3 matrix
# cofactor expansion about first column
  val1 = A[1,1]*(A[2,2]*A[3,3] - A[2,3]*A[3,2])
  val2 = -A[2,1]*(A[1,2]*A[3,3] - A[1,3]*A[3,2])
  val3 = A[3,1]*(A[1,2]*A[2,3] - A[1,3]*A[2,2])
  return val1 + val2 + val3
end

function adjugate2(A::AbstractMatrix, B::AbstractMatrix)
  B[1,1] = A[2,2]
  B[1,2] = -A[1,2]
  B[2,1] = -A[2,1]
  B[2,2] = A[1,1]
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


