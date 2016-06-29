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




# interpolates dxidx, jac to the face nodes
# only do this for elementL of each interface
function interpolateMapping{Tmsh}(mesh::PumiMeshDG{Tmsh})
  sbpface = mesh.sbpface
  dim = mesh.dim
  dxidx_face = zeros(Tmsh, dim, dim, sbpface.numnodes, mesh.numInterfaces)
  dxidx_sharedface = Array(Array{Tmsh, 4}, mesh.npeers)
  for i=1:mesh.npeers
    dxidx_sharedface[i] = zeros(Tmsh, dim, dim, sbpface.numnodes, mesh.peer_face_counts[i])
  end
  jac_face = zeros(Tmsh, sbpface.numnodes, mesh.numInterfaces)

  dxidx_bndry = zeros(Tmsh, dim, dim, sbpface.numnodes, mesh.numBoundaryFaces)
  jac_bndry = zeros(Tmsh, sbpface.numnodes, mesh.numBoundaryFaces)
  jac_sharedface = Array(Array{Tmsh, 2}, mesh.npeers)
  for i=1:mesh.npeers
    jac_sharedface[i] = Array(Tmsh, sbpface.numnodes, mesh.peer_face_counts[i])
  end

  # temporary storage
  interp_data = Interpolation(mesh)
  myrank = mesh.myrank
  f = open("interpolation$myrank.dat", "w")

  for i=1:mesh.numInterfaces
    dxidx_i = view(dxidx_face, :, :, :, i)
    jac_i = view(jac_face, :, i)

    interface_i = mesh.interfaces[i]
    el = interface_i.elementL
    face = interface_i.faceL
    println(f, "interface ", i, ": element ", el, ", face ", face)
    interp_data.bndry_arr[1] =  Boundary(1, face)

    dxidx_in = view(mesh.dxidx, :, :, :, el)
    jac_in = view(mesh.jac, :, el)

    interpolateFace(interp_data, mesh.sbpface, dxidx_in, jac_in, dxidx_i, jac_i)
    println(f, "dxidx_in = ", dxidx_in)
    println(f, "dxidx_interp = ", dxidx_i)
    println(f, "jac_in = ", jac_in)
    println(f, "jac_interp = ", jac_i)

  end
  # now do shared edges
  for peer = 1:mesh.npeers
    dxidx_p = dxidx_sharedface[peer]
    jac_p = jac_sharedface[peer]
    interfaces_p = mesh.shared_interfaces[peer]
    for i=1:mesh.peer_face_counts[peer]
      println(f, "peer ", peer, ", face ", i)
      dxidx_i = view(dxidx_p, :, :, :, i)
      jac_i = view(jac_p, :, i)

      interface_i = interfaces_p[i]
      el = interface_i.elementL
      face = interface_i.faceL

      println(f, "element ", el, ", face ", face)
      interp_data.bndry_arr[1] = Boundary(1, face)


      dxidx_in = view(mesh.dxidx, :, :, :, el)
      jac_in = view(mesh.jac, :, el)

      interpolateFace(interp_data, mesh.sbpface, dxidx_in, jac_in, dxidx_i, jac_i)
      println(f, "dxidx_in = ", dxidx_in)
      println(f, "dxidx_interp = ", dxidx_i)
      println(f, "jac_in = ", jac_in)
      println(f, "jac_interp = ", jac_i)
    end
  end
  close(f)

  # now do boundary
  for i=1:mesh.numBoundaryFaces
    bndry = mesh.bndryfaces[i]

    dxidx_i = view(dxidx_bndry, :, :, :, i)
    jac_i = view(jac_bndry, :, i)

    el = bndry.element
    dxidx_in = view(mesh.dxidx, :, :, :, el)
    jac_in = view(mesh.jac, :, el)
    interp_data.bndry_arr[1] = bndry
    interpolateFace(interp_data, mesh.sbpface, dxidx_in, jac_in, dxidx_i, jac_i)
  end

  return dxidx_face, jac_face, dxidx_sharedface, jac_sharedface, dxidx_bndry, jac_bndry

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
    dxidx_hat = view(dxidx_hat_in, :, :, j)
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
    dxidx_hat = view(dxidx_hat_in, :, :, j)
    detJ = jac_in[j]

    adjugate3(dxidx_hat, dxidx_node)
    for k=1:dim
      for p=1:dim
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


