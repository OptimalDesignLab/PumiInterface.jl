# file for functions related to interpolating from solution points to 
# faces

# interpolates dxidx, jac to the face nodes
# only do this for elementL of each interface
function interpolateMapping{Tmsh}(mesh::PumiMeshDG2{Tmsh})
  sbpface = mesh.sbpface

  dxidx_face = zeros(Tmsh, 2, 2, sbpface.numnodes, mesh.numInterfaces)
  dxidx_sharedface = Array(Array{Tmsh, 4}, mesh.npeers)
  for i=1:mesh.npeers
    dxidx_sharedface[i] = zeros(Tmsh, 2, 2, sbpface.numnodes, mesh.peer_face_counts[i])
  end
  jac_face = zeros(Tmsh, sbpface.numnodes, mesh.numInterfaces)

  dxidx_bndry = zeros(Tmsh, 2, 2, sbpface.numnodes, mesh.numBoundaryEdges)
  jac_bndry = zeros(Tmsh, sbpface.numnodes, mesh.numBoundaryEdges)
  jac_sharedface = Array(Array{Tmsh, 2}, mesh.npeers)
  for i=1:mesh.npeers
    jac_sharedface[i] = Array(Tmsh, sbpface.numnodes, mesh.peer_face_counts[i])
  end

  dxdxi_el = zeros(Tmsh, 4, mesh.numNodesPerElement, 1)
  dxdxi_elface = zeros(Tmsh, 4, sbpface.numnodes, 1)
  dxidx_node = zeros(Tmsh, 2, 2)
  dxdxi_node = zeros(Tmsh, 2, 2)
  bndry_arr = Array(Boundary, 1)
  for i=1:mesh.numInterfaces
    dxidx_i = view(dxidx_face, :, :, :, i)
    jac_i = view(jac_face, :, i)

    interface_i = mesh.interfaces[i]
    el = interface_i.elementL
    face = interface_i.faceL
    bndry = Boundary(1, face)

    dxidx_in = view(mesh.dxidx, :, :, :, el)
    jac_in = view(mesh.jac, :, el)

    interpolateFace(bndry, mesh.sbpface, dxidx_in, jac_in, dxdxi_el, dxdxi_elface, dxdxi_node, dxidx_node, dxidx_i, jac_i)

  end

  # now do shared edges
  for peer = 1:mesh.npeers
    dxidx_p = dxidx_sharedface[peer]
    jac_p = jac_sharedface[peer]
    interfaces_p = mesh.shared_interfaces[peer]
    for i=1:mesh.peer_face_counts[peer]
      dxidx_i = view(dxidx_p, :, :, :, i)
      jac_i = view(jac_p, :, i)

      interface_i = interfaces_p[i]
      el = interface_i.elementL
      face = interface_i.faceL
      bndry = Boundary(1, face)

      dxidx_in = view(mesh.dxidx, :, :, :, el)
      jac_in = view(mesh.jac, :, el)

      interpolateFace(bndry, mesh.sbpface, dxidx_in, jac_in, dxdxi_el, dxdxi_elface, dxdxi_node, dxidx_node, dxidx_i, jac_i)
    end
  end

  # now do boundary
  for i=1:mesh.numBoundaryEdges
    bndry = mesh.bndryfaces[i]

    dxidx_i = view(dxidx_bndry, :, :, :, i)
    jac_i = view(jac_bndry, :, i)

    el = bndry.element
    dxidx_in = view(mesh.dxidx, :, :, :, el)
    jac_in = view(mesh.jac, :, el)

    interpolateFace(bndry, mesh.sbpface, dxidx_in, jac_in, dxdxi_el, dxdxi_elface, dxdxi_node, dxidx_node, dxidx_i, jac_i)
  end

  return dxidx_face, jac_face, dxidx_sharedface, jac_sharedface, dxidx_bndry, jac_bndry

end  # end function

function interpolateFace(bndry::Boundary, sbpface, dxidx_hat_in, jac_in, dxdxi_el, dxdxi_elface, dxdxi_node, dxidx_node, dxidx_i, jac_i)

  numNodesPerElement = size(dxidx_hat_in, 3)
  # get the data
  for j=1:numNodesPerElement
    dxidx_hat = view(dxidx_hat_in, :, :, j)
    detJ = jac_in[j]

    # dxdxi = inv(dxidx) = adj(dxidx_hat)
    dxdxi_node[1,1] = dxidx_hat[2,2]
    dxdxi_node[1,2] = -dxidx_hat[1,2]
    dxdxi_node[2,1] = -dxidx_hat[2,1]
    dxdxi_node[2,2] = dxidx_hat[1,1]

    dxdxi_el[1,j,1] = dxdxi_node[1,1]
    dxdxi_el[2,j,1] = dxdxi_node[1,2]
    dxdxi_el[3,j,1] = dxdxi_node[2,1]
    dxdxi_el[4,j,1] = dxdxi_node[2,2]

  end


  # interpolate to the face
  face = bndry.face
  bndry_arr = [Boundary(1, face)]
#    bndry_arr[1] = Boundary(1, face)
  boundaryinterpolate!(sbpface, bndry_arr, dxdxi_el, dxdxi_elface)


  # now store dxidx, |J| at the boundary nodes
  for j=1:sbpface.numnodes
    dxdxi_node[1,1] = dxdxi_elface[1, j, 1]
    dxdxi_node[1,2] = dxdxi_elface[2, j, 1]
    dxdxi_node[2,1] = dxdxi_elface[3, j, 1]
    dxdxi_node[2,2] = dxdxi_elface[4, j, 1]


    # inv(A) = adj(A)/|A|
    det_dxdxi = dxdxi_node[1,1]*dxdxi_node[2,2] - dxdxi_node[1,2]*dxdxi_node[2,1]
    dxidx_node[1,1] = dxdxi_node[2,2]/det_dxdxi
    dxidx_node[1,2] = -dxdxi_node[1,2]/det_dxdxi
    dxidx_node[2,1] = -dxdxi_node[2,1] /det_dxdxi
    dxidx_node[2,2] = dxdxi_node[1,1]/det_dxdxi

    detJ = dxidx_node[1,1]*dxidx_node[2,2] - dxidx_node[1,2]*dxidx_node[2,1]
    dxidx_i[1,1,j] = dxidx_node[1,1]/detJ
    dxidx_i[1,2,j] = dxidx_node[1,2]/detJ
    dxidx_i[2,1,j] = dxidx_node[2,1]/detJ
    dxidx_i[2,2,j] = dxidx_node[2,2]/detJ
    jac_i[j] = detJ
  end

  return nothing

end


