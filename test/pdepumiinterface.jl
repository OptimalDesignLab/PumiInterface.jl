function filter_noise(A::AbstractArray)
  for i=1:length(A)
    if abs(A[i]) < 1e-8
      A[i] = 0
    end
  end

  return A
end

"""
  Test reverse mode of metric interpolation function
"""
function test_interp_rev(mesh)
  facts("\n ----- Testing interp rev -----") do

    dim = mesh.dim

    interp_data = PdePumiInterface.Interpolation(mesh)
    interp_datac = PdePumiInterface.Interpolation(mesh, Complex128)
    interp_data.bndry_arr[1] = Boundary(2, 1)
    interp_datac.bndry_arr[1] = Boundary(2, 1)

    if mesh.dim == 2
      dxidx_in = zeros(dim, dim, mesh.numNodesPerElement)
      jac_in = zeros(mesh.numNodesPerElement)

      for i=1:length(dxidx_in)
        dxidx_in[i] = i
      end
      for i=1:length(jac_in)
        jac_in[i] = i
      end
    else
      dxidx_in = mesh.dxidx[ :, :, :, 1]
      jac_in = mesh.jac[ :, 1]
    end


    dxidx_out = zeros(dim, dim, mesh.numNodesPerFace)
    jac_out = zeros(mesh.numNodesPerFace)

    PdePumiInterface.interpolateFace(interp_data, mesh.sbpface, dxidx_in, jac_in, dxidx_out, jac_out)

    # compute forward jacobian
    jac_size = dim*dim*mesh.numNodesPerFace
    jac_forward = zeros(jac_size, jac_size)
    
    dxidx_in2 = zeros(Complex128, size(dxidx_in)...)
    copy!(dxidx_in2, dxidx_in)
    jac_in2 = zeros(Complex128, size(jac_in)...)
    copy!(jac_in2, jac_in)

    dxidx_out2 = zeros(Complex128, dim, dim, mesh.numNodesPerFace)
    jac_out2 = zeros(Complex128, mesh.numNodesPerFace)


    h =1e-20
    pert = Complex128(0, h)

    for i=1:jac_size
      dxidx_in2[i] += pert
      PdePumiInterface.interpolateFace(interp_datac, mesh.sbpface, dxidx_in2, jac_in2, dxidx_out2, jac_out2)
      
      for j=1:jac_size
        jac_forward[j, i] = imag(dxidx_out2[j])/h
      end

      dxidx_in2[i] -= pert
    end

    # compute reverse mode
    jac_rev = zeros(jac_size, jac_size)
    dxidx_in_bar = zeros(dxidx_in)
    jac_in_bar = zeros(jac_in)
    dxidx_out_bar = zeros(dxidx_out)
    jac_out_bar = zeros(jac_out)

    for i=1:jac_size
      dxidx_out_bar[i] = 1
      fill!(dxidx_in_bar, 0.0)
      fill!(jac_in_bar, 0.0)
      fill!(jac_out_bar, 0.0)

      PdePumiInterface.interpolateFace_rev(interp_data, mesh.sbpface, dxidx_in, 
                                       dxidx_in_bar, jac_in, jac_in_bar, 
                                       dxidx_out, dxidx_out_bar, jac_out, 
                                       jac_out_bar)

      for j=1:jac_size
        jac_rev[i, j] = dxidx_in_bar[j]
      end

      dxidx_out_bar[i] -= 1
    end


    @fact norm(jac_forward - jac_rev) --> roughly(0.0, atol=1e-11)

    # now test the full routine
    fill!(mesh.dxidx_bar, 0.0)
#    mesh.dxidx_face_bar[1] = 1
    fill!(mesh.dxidx_face_bar, 1.0)

    # check regular interfaces
    PdePumiInterface.interpolateMapping_rev(mesh)
    for i=1:mesh.numInterfaces
      el = mesh.interfaces[i].elementL
      for j=1:mesh.numNodesPerFace
        @fact norm(mesh.dxidx_bar[:, :, j, el]) --> greater_than(0.0)
      end
    end
    
    # check boundary faces
    fill!(mesh.dxidx_bar, 0.0)
    fill!(mesh.dxidx_bndry_bar, 1.0)
    PdePumiInterface.interpolateMapping_rev(mesh)
    for i=1:mesh.numBoundaryFaces
      el = mesh.bndryfaces[i].element
      for j=1:mesh.numNodesPerFace
        @fact norm(mesh.dxidx_bar[:, :, j, el]) --> greater_than(0.0)
      end
    end
 
    fill!(mesh.dxidx_bndry_bar, 1.0)
    
    # check sharedfaces
    for peer = 1:mesh.npeers
      fill!(mesh.dxidx_bar, 0.0)
      fill!(mesh.dxidx_sharedface[peer], 1.0)
      PdePumiInterface.interpolateMapping_rev(mesh)
      for i=1:mesh.numBoundaryFaces
        el = mesh.shared_interfacesfaces[p][i].elementL
        for j=1:mesh.numNodesPerFace
          @fact norm(mesh.dxidx_bar[:, :, j, el]) --> greater_than(0.0)
        end
      end

      fill!(mesh.dxidx_bar, 0.0)
      fill!(mesh.dxidx_sharedface[peer]. 0.0)
    end
   
  end

end  # end function

"""
  Test the reverse mode of the curvilinear metric calculation
"""
function test_metric_rev(mesh, sbp)

  facts("----- testing metric reverse mode -----") do
    sbpface = mesh.sbpface
    nrm = rand(Complex128, mesh.dim, mesh.numNodesPerFace)
    nrm[:, :] = real(nrm)
    Rone = ones(mesh.numNodesPerFace)
    tmp = zeros(Complex128, mesh.numNodesPerFace)

    # test calcEoneElement
    Eone_el_c = zeros(Complex128, sbpface.stencilsize, mesh.dim)
    
    nin = length(nrm)
    nout = length(Eone_el_c)
    jac = zeros(nin, nout)
    jac2 = zeros(jac)

    h = 1e-20
    pert = Complex128(0, h)

    println("computing forward mode")
    for i=1:nin
      nrm[i] += pert
      fill!(Eone_el_c, 0.0)
      PdePumiInterface.calcEoneElement(sbpface, nrm, Rone, tmp, Eone_el_c)

      for j=1:nout
        jac[i, j] = imag(Eone_el_c[j])/h
      end

      nrm[i] -= pert
    end

    nrm_bar = zeros(mesh.dim, mesh.numNodesPerFace)
    tmp_bar = zeros(Float64, size(tmp)...)
    Eone_el_bar = zeros(Float64, size(Eone_el_c)...)

    println("computing reverse mode")
    for j=1:nout
      Eone_el_bar[j] = 1
      fill!(nrm_bar, 0.0)
      PdePumiInterface.calcEoneElement_rev(sbpface, nrm_bar, Rone, tmp_bar, Eone_el_bar)

      for i=1:nin
        jac2[i, j] = nrm_bar[i]
      end

      Eone_el_bar[j] = 0
    end

    @fact norm(jac - jac2)/length(jac2) --> roughly(0.0, atol=1e-12)


    # test assembleEone
    Eone_c = zeros(Complex128, mesh.numNodesPerElement, mesh.dim, 1)
    Eone_el_c[:, :] = real(Eone_el_c)
    elnum = 1
    facenum_local = 1
    
    nin = length(Eone_el_c)
    nout = length(Eone_c)

    jac = zeros(nin, nout)
    jac2 = zeros(jac)

    println("testing forward mode")
    for i=1:nin
      Eone_el_c[i] += pert
      fill!(Eone_c, 0.0)
      PdePumiInterface.assembleEone(sbpface, elnum, facenum_local, Eone_el_c, Eone_c)
      for j=1:nout
        jac[i, j] = imag(Eone_c[j])/h
      end

      Eone_el_c[i] -= pert
    end

    println("testing reverse mode")
    Eone_bar = zeros(Float64, mesh.numNodesPerElement, mesh.dim, 1)


    for j=1:nout
      Eone_bar[j] = 1

      fill!(Eone_el_bar, 0.0)
      PdePumiInterface.assembleEone_rev(sbpface, elnum, facenum_local, Eone_el_bar, Eone_bar)

      for i=1:nin
        jac2[i, j] = Eone_el_bar[i]
      end


      Eone_bar[j] = 0
    end

    @fact norm(jac - jac2) --> roughly(0.0, atol=1e-12)


    test_metric2_rev(mesh, sbp)
    test_metrics3_rev(mesh, sbp)
    test_metrics4_rev(mesh, sbp)

  end  # end facts block

  return nothing
end

"""
  Test calcEone_rev
"""
function test_metric2_rev(mesh, sbp)

  facts("----- Testing second part of Metric reverse mode -----") do
    # test calcEone_rev
    # use finite differences because the perturbation is applied to the fields
    # of the mesh

    pert = 1e-6
    nin = length(mesh.nrm_face) + length(mesh.nrm_bndry)
    for i=1:mesh.npeers
      nin += length(mesh.nrm_sharedface[i])
    end
    nout = mesh.numNodesPerElement*mesh.dim*mesh.numEl

    Eone = zeros(mesh.numNodesPerElement, mesh.dim, mesh.numEl)
    Eone_orig = zeros(Eone)

    jac = zeros(nin, nout)
    jac2 = zeros(jac)
    element_range = 1:mesh.numEl

    PdePumiInterface.calcEone(mesh, sbp, element_range, Eone_orig)

    in_idx = 1
    for i=1:length(mesh.nrm_face)
      mesh.nrm_face[i] += pert
      fill!(Eone, 0.0)
      PdePumiInterface.calcEone(mesh, sbp, element_range, Eone)

      for j=1:nout
        jac[in_idx, j] = (Eone[j] - Eone_orig[j])/pert
      end

      in_idx += 1
      mesh.nrm_face[i] -= pert
    end

    for i=1:length(mesh.nrm_bndry)
      mesh.nrm_bndry[i] += pert
      fill!(Eone, 0.0)
      PdePumiInterface.calcEone(mesh, sbp, element_range, Eone)
      for j=1:nout
        jac[in_idx, j] = (Eone[j] - Eone_orig[j])/pert
      end

      in_idx += 1
      mesh.nrm_bndry[i] -= pert
    end

    for peer=1:mesh.npeers
      nrm_peer = mesh.nrm_sharedface[peer]
      for i=1:length(nrm_peer)
        nrm_peer[i] += pert
        fill!(Eone, 0.0)
        PdePumiInterface.calcEone(mesh, sbp, element_range, Eone)
        for j=1:nout
          jac[in_idx, j] = (Eone[j] - Eone_orig[j])/pert
        end

        in_idx += 1
        nrm_peer[i] -= pert
      end
    end

    # reverse mode
    Eone_bar = zeros(Eone)
    for j=1:nout
      Eone_bar[j] = 1

      fill!(mesh.nrm_face_bar, 0.0)
      fill!(mesh.nrm_bndry_bar, 0.0)
      for i=1:mesh.npeers
        fill!(mesh.nrm_sharedface[i], 0.0)
      end
      PdePumiInterface.calcEone_rev(mesh, sbp, element_range, Eone_bar)

      in_idx = 1
      for i=1:length(mesh.nrm_face_bar)
        jac2[in_idx, j] = mesh.nrm_face_bar[i]
        in_idx += 1
      end

      for i=1:length(mesh.nrm_bndry_bar)
        jac2[in_idx, j] = mesh.nrm_bndry_bar[i]
        in_idx += 1
      end

      for peer=1:mesh.npeers
        nrm_bar_peer = mesh.nrm_sharedface_bar[peer]
        for i=1:length(nrm_bar_peer)
          jac2[in_idx, j] = nrm_bar_peer[i]
          in_idx += 1
        end
      end

      Eone_bar[j] = 0

    end  # end loop j

    @fact norm(jac - jac2)/length(jac) --> roughly(0.0, atol=1e-5)



  end  # end facts block

  return nothing
end

"""
  Test getCurvilinearCoordinatesAndMetrics_rev (dxidx -> vert_coord, nrm)
"""
function test_metrics3_rev(mesh, sbp)

  facts("----- testing metrics reverse mode 3 -----") do
#    nout = mesh.dim*mesh.numNodesPerElement*2
    nout = length(mesh.dxidx) + length(mesh.jac)
    nin = length(mesh.vert_coords) + length(mesh.nrm_face) + length(mesh.nrm_bndry)
    for i=1:mesh.npeers
      nin += length(mesh.nrm_sharedface[i])
    end
#    nin = length(mesh.vert_coords) + mesh.dim*mesh.numNodesPerElement*2

    jac = zeros(nin, nout)
    jac2 = zeros(jac)

    pert = 1e-6

    dxidx_orig = copy(mesh.dxidx)
    jac_orig = copy(mesh.jac)
#    for i=1:mesh.npeers
#      nrm_sharedface[i] = copy(mesh.nrm_sharedface[i])
#    end

    # forward mode
    in_idx = 1
    for i=1:length(mesh.vert_coords)
      mesh.vert_coords[i] += pert

      PdePumiInterface.getCurvilinearCoordinatesAndMetrics(mesh, sbp)

      out_idx = 1
      for j=1:length(mesh.dxidx)
        jac[in_idx, out_idx] = (mesh.dxidx[j] - dxidx_orig[j])/pert
        out_idx += 1
      end

      # TODO: uncomment when this is fixed
      for j=1:length(mesh.jac)
#        jac[in_idx, out_idx] = (mesh.jac[j] - jac_orig[j])/pert
        out_idx += 1
      end

      in_idx += 1
      mesh.vert_coords[i] -= pert
    end

    for i=1:length(mesh.nrm_face)
      mesh.nrm_face[i] += pert
      PdePumiInterface.getCurvilinearCoordinatesAndMetrics(mesh, sbp)

      out_idx = 1
      for j=1:length(mesh.dxidx)
        jac[in_idx, out_idx] = (mesh.dxidx[j] - dxidx_orig[j])/pert
        out_idx += 1
      end

      for j=1:length(mesh.jac)
        jac[in_idx, out_idx] = (mesh.jac[j] - jac_orig[j])/pert
        out_idx += 1
      end

      in_idx += 1
      mesh.nrm_face[i] -= pert
    end

    for i=1:length(mesh.nrm_bndry)
      mesh.nrm_bndry[i] += pert
      PdePumiInterface.getCurvilinearCoordinatesAndMetrics(mesh, sbp)

      out_idx = 1
      for j=1:length(mesh.dxidx)
        jac[in_idx, out_idx] = (mesh.dxidx[j] - dxidx_orig[j])/pert
        out_idx += 1
      end

      # TODO: uncomment when this is fixed
      for j=1:length(mesh.jac)
#        jac[in_idx, out_idx] = (mesh.jac[j] - jac_orig[j])/pert
        out_idx += 1
      end

      in_idx += 1
      mesh.nrm_bndry[i] -= pert
    end

    for peer=1:mesh.npeers
      nrm_peer = mesh.nrm_sharedface[peer]
      for i=1:length(nrm_peer)
        nrm_peer[i] += pert
        PdePumiInterface.getCurvilinearCoordinatesAndMetrics(mesh, sbp)

        for j=1:length(mesh.dxidx)
          jac[in_idx, out_idx] = (mesh.dxidx[j] - dxidx_orig[j])/pert
          out_idx += 1
        end

        # uncomment when this is fixed
        for j=1:length(mesh.jac)
#          jac[in_idx, out_idx] = (mesh.jac[j] - jac_orig[j])/pert
          out_idx += 1
        end

        in_idx += 1
        nrm_peer[i] -= pert
      end
    end

    # reverse mode
    println("running reverse mode")
    fill!(mesh.dxidx_bar, 0.0)
    fill!(mesh.jac_bar, 0.0)
    out_idx = 1
    for j=1:length(mesh.dxidx)
      mesh.dxidx_bar[j] = 1
#      Eone_bar[j] = 1

      fill!(mesh.vert_coords_bar, 0.0)
      fill!(mesh.nrm_face_bar, 0.0)
      fill!(mesh.nrm_bndry_bar, 0.0)
      for i=1:mesh.npeers
        fill!(mesh.nrm_sharedface_bar[i], 0.0)
      end
      fill!(mesh.jac_bar, 0.0)
      PdePumiInterface.getCurvilinearCoordinatesAndMetrics_rev(mesh, sbp)

      in_idx = 1
      for i=1:length(mesh.vert_coords)
        jac2[in_idx, out_idx] = mesh.vert_coords_bar[i]
        in_idx += 1
      end

      for i=1:length(mesh.nrm_face)
        jac2[in_idx, out_idx] = mesh.nrm_face_bar[i]
        in_idx += 1
      end

      
      for i=1:length(mesh.nrm_bndry)
        jac2[in_idx, out_idx] = mesh.nrm_bndry_bar[i]
        in_idx += 1
      end

      for peer=1:mesh.npeers
        nrm_peer_bar = mesh.nrm_sharedface_bar[peer]
        for i=1:length(nrm_peer)
          jac2[in_idx, out_idx] = nrm_peer_bar[i]
          in_idx += 1
        end
      end
      
      out_idx += 1
      mesh.dxidx_bar[j] = 0
    end  # end loop j

    fill!(mesh.jac_bar, 0.0)
    fill!(mesh.dxidx_bar, 0.0)
    for j=1:length(mesh.jac)
      #TODO: uncomment when this is fixed
#      mesh.jac_bar[j] = 1

      fill!(mesh.vert_coords_bar, 0.0)
      fill!(mesh.nrm_face_bar, 0.0)
      fill!(mesh.nrm_bndry_bar, 0.0)
      for i=1:mesh.npeers
        fill!(mesh.nrm_sharedface_bar[i], 0.0)
      end
      PdePumiInterface.getCurvilinearCoordinatesAndMetrics_rev(mesh, sbp)

      in_idx = 1
      for i=1:length(mesh.vert_coords)
        jac2[in_idx, out_idx] = mesh.vert_coords_bar[i]
        in_idx += 1
      end

      for i=1:length(mesh.nrm_face)
        jac2[in_idx, out_idx] = mesh.nrm_face_bar[i]
        in_idx += 1
      end

      
      for i=1:length(mesh.nrm_bndry)
        jac2[in_idx, out_idx] = mesh.nrm_bndry_bar[i]
        in_idx += 1
      end

      for peer=1:mesh.npeers
        nrm_peer_bar = mesh.nrm_sharedface_bar[peer]
        for i=1:length(nrm_peer)
          jac2[in_idx, out_idx] = nrm_peer_bar[i]
          in_idx += 1
        end
      end
      
      out_idx += 1
      mesh.jac_bar[j] = 0
    end  # end loop j

#    println("jac = \n", jac)
#    println("jac2 = \n", jac2)
#    println("diff = \n", jac - jac2)

    @fact norm(jac - jac2)/length(jac) --> roughly(0.0, atol=1e-5)

  end

  return nothing
end


function test_metrics4_rev(mesh, sbp)

  # zero out all bar variables, just in case
  fill!(mesh.dxidx_bar, 0.0)
  fill!(mesh.jac_bar, 0.0)
  fill!(mesh.vert_coords_bar, 0.0)
  fill!(mesh.nrm_face_bar, 0.0)
  fill!(mesh.nrm_bndry_bar, 0.0)
  for i=1:mesh.npeers
    fill!(mesh.nrm_sharedface_bar[i], 0.0)
  end

  facts("----- Testing metrics4_rev -----") do
    nin = length(mesh.vert_coords)
    nout = length(mesh.nrm_bndry) + length(mesh.nrm_face)
    for i=1:mesh.npeers
      nout += length(mesh.nrm_sharedface[i])
    end

    jac = zeros(nin, nout)
    jac2 = zeros(jac)

    PdePumiInterface.getFaceCoordinatesAndNormals(mesh, sbp)
    nrm_bndry_orig = copy(mesh.nrm_bndry)
    nrm_face_orig = copy(mesh.nrm_face)
    nrm_sharedface_orig = copy(mesh.nrm_sharedface)
    for i=1:mesh.npeers
      nrm_sharedface_orig[i] = copy(mesh.nrm_sharedface[i])
    end

    pert = 1e-6

    # forward mode
    in_idx = 1
    for i=1:length(mesh.vert_coords)
      mesh.vert_coords[i] += pert

      PdePumiInterface.getFaceCoordinatesAndNormals(mesh, sbp)

      out_idx = 1

      for j=1:length(mesh.nrm_bndry)
        jac[in_idx, out_idx] = (mesh.nrm_bndry[j] - nrm_bndry_orig[j])/pert
        out_idx += 1
      end

      for j=1:length(mesh.nrm_face)
        jac[in_idx, out_idx] = (mesh.nrm_face[j] - nrm_face_orig[j])/pert
        out_idx += 1
      end

      for peer=1:mesh.npeers
        for j=1:length(mesh.nrm_sharedface[peer])
          jac[in_idx, out_idx] = (mesh.nrm_sharedface[peer][j] - nrm_sharedface_orig[peer][j])/pert
          out_idx += 1
        end
      end


      in_idx += 1
      mesh.vert_coords[i] -= pert
    end

    # reverse mode
    println("reverse mode")
    out_idx = 1
    for j=1:length(mesh.nrm_bndry)
      mesh.nrm_bndry_bar[j] = 1
      fill!(mesh.vert_coords_bar, 0.0)

      PdePumiInterface.getFaceCoordinatesAndNormals_rev(mesh, sbp)

      in_idx = 1
      for i=1:length(mesh.vert_coords_bar)
        jac2[in_idx, out_idx] = mesh.vert_coords_bar[i]
        in_idx += 1
      end

      out_idx += 1
      mesh.nrm_bndry_bar[j] = 0
    end


    for j=1:length(mesh.nrm_face)
      mesh.nrm_face_bar[j] = 1
      fill!(mesh.vert_coords_bar, 0.0)

      PdePumiInterface.getFaceCoordinatesAndNormals_rev(mesh, sbp)

      in_idx = 1
      for i=1:length(mesh.vert_coords_bar)
        jac2[in_idx, out_idx] = mesh.vert_coords_bar[i]
        in_idx += 1
      end

      out_idx += 1
      mesh.nrm_face_bar[j] = 0
    end

    for peer=1:mesh.npeers
      nrm_sharedface_bar = mesh.nrm_sharedface_bar[peer]
      for j=1:length(nrm_sharedface_bar)
        nrm_sharedface_bar[j] = 1
        fill!(mesh.vert_coords_bar, 0.0)

        PdePumiInterface.getFaceCoordinatesAndNormals_rev(mesh, sbp)

        in_idx = 1
        for i=1:length(mesh.vert_coords_bar)
          jac2[in_idx, out_idx] = mesh.vert_coords_bar[i]
          in_idx += 1
        end

        out_idx += 1
        nrm_sharedface_bar[j] = 0
      end
    end

    diffnorm = norm(jac - jac2)/length(jac)
    @fact diffnorm --> roughly(0.0, atol=1e-5)
  end

  error("stop here")
  return nothing
end


  

function test_coords_rev(mesh, sbp)

#  sbp_complex = TriSBP{Complex128}(degree=sbp.degree, reorder=sbp.reorder, internal=sbp.internal, vertices=sbp.vertices)

#  function coords_forward(vert_coords)
#    coords_it = vert_coords.'
#    coords_volume = SummationByParts.SymCubatures.calcnodes(sbp.cub, coords_it)
#    return coords_volume
#  end
  numVertsPerElement = mesh.numTypePerElement[1]
  # define indexing for input and output dimensions
  idx1(dim, vert, el) = dim + mesh.dim*(vert-1) + (el-1)*numVertsPerElement*mesh.dim
  idx2(dim, node, el) = dim + mesh.dim*(node-1) + (el-1)*mesh.numNodesPerElement*mesh.dim

  facts("----- Testing coordiate reverse mode -----") do
  # compute the jacobian with forward mode
  nin = mesh.dim*numVertsPerElement*mesh.numEl
  nout = mesh.dim*mesh.numNodesPerElement*mesh.numEl
  jac = zeros(nin, nout)

  vertcoords_in = zeros(Complex128, size(mesh.vert_coords)...)
#  coords_out = zeros(Complex128, size(mesh.coords)...)

  copy!(vertcoords_in, mesh.vert_coords)
  h = 1e-20
  pert = Complex128(0, h)
  for el=1:mesh.numEl
    for vert=1:numVertsPerElement
      for d=1:mesh.dim
        idx_in = idx1(d, vert, el)

        vertcoords_in[d, vert, el] += pert

        vertcoords_el = sview(vertcoords_in, :, :, el)
        coords_el = PdePumiInterface.vertToVolumeCoords(mesh, sbp, vertcoords_el)

        vertcoords_in[d, vert, el] -= pert

        # take advantage of the fact that the operation is element local
        for node_out=1:mesh.numNodesPerElement
          for d2=1:mesh.dim
            idx_out = idx1(d2, node_out, el)
            jac[idx_out, idx_in] = imag(coords_el[d2, node_out])/h
          end
        end

      end
    end
  end

  jac2 = zeros(jac)
  coords_bar = zeros(mesh.dim, mesh.numNodesPerElement, mesh.numEl)
  vertcoords_bar = zeros(size(mesh.vert_coords))
  # construct reverse mode jacobian
  for el=1:mesh.numEl
    for node=1:mesh.numNodesPerElement
      for d=1:mesh.dim
        coords_bar[d, node, el] = 1
        PdePumiInterface.volumeToVertCoords_rev(mesh, sbp, vertcoords_bar, coords_bar)

        coords_bar[d, node, el] = 0
        idx_out = idx1(d, node, el)
        for vert_in = 1:numVertsPerElement
          for d2=1:mesh.dim
            idx_in = idx1(d2, vert_in, el)
            jac2[idx_out, idx_in] = vertcoords_bar[d2, vert_in, el]
          end
        end

        fill!(vertcoords_bar, 0.0)

      end
    end
  end

  for i=1:nout
    for j=1:nin
      @fact jac[j, i] --> roughly(jac2[j, i], atol=1e-13)
    end
  end


  end

  return nothing
end



facts("--- Testing PdePumiInterface --- ") do

  opts = Dict{Any, Any}(
    "numBC" => 1,
    "BC1" =>  [0],
    "run_type" => 4,
    "verify_coloring" => true,
    "use_edge_res" => false,
    "write_edge_vertnums" => true,
    "write_face_vertnums" => true,
    "write_boundarynums" => true,
    "write_dxidx" => true,
    "write_coords" => true,
    "write_sparsity" => true,
    "write_offsets" => true,
    "write_counts" => true,
    "write_sparsity_nodebnds" => true,
    "write_dofs" => false,
    "write_interfaces" => false,
    "write_boundaries" => false,
    "write_sharedboundaries" => false,
    )

    # names of output files
    fnames = ["boundary_nums", "face_vertnums", "edge_vertnums", "sparsity_bnds", "sparsity_bnds", "sparsity_nodebnds"]

    smb_name = "tri2l.smb"
    dmg_name = ".null"
    for order = 1:4
      println("testing order ", order, " CG mesh")
      sbp = TriSBP{Float64}(degree=order)
      ref_verts = [-1. 1 -1; -1 -1 1]
      sbpface = TriFace{Float64}(order, sbp.cub, ref_verts.')


    mesh =  PumiMesh2{Float64}(dmg_name, smb_name, order, sbp, opts, sbpface, coloring_distance=2, dofpernode=4)
    @fact mesh.numVert --> 4
    @fact mesh.numEdge --> 5
    @fact mesh.numFace --> mesh.numEdge
    @fact mesh.numEl --> 2
    @fact mesh.numEntitiesPerType --> [4, 5, 2]
    @fact mesh.numTypePerElement --> [3, 3, 1]
    @fact mesh.numDofPerNode --> 4
    @fact mesh.numBoundaryFaces --> 4
    @fact mesh.numInterfaces --> 1

    @fact mesh.bndryfaces[1].element --> 1
    @fact mesh.bndryfaces[1].face --> 3
    @fact mesh.bndryfaces[2].element --> 2
    @fact mesh.bndryfaces[2].face --> 1
    @fact mesh.bndryfaces[3].element --> 1
    @fact mesh.bndryfaces[3].face --> 2
    @fact mesh.bndryfaces[4].element --> 2
    @fact mesh.bndryfaces[4].face --> 2

  #  println("mesh.interfaces = ",  mesh.interfaces)
    @fact mesh.interfaces[1].elementL --> 1
    @fact mesh.interfaces[1].elementR --> 2
    @fact mesh.interfaces[1].faceL --> 1
    @fact mesh.interfaces[1].faceR --> 3

    @fact mesh.order --> order
    @fact length(mesh.bndry_funcs) --> 1
    @fact mesh.bndry_offsets --> [1, 5]
    @fact mesh.bndry_geo_nums[1] --> opts["BC1"]

#=
    for i=1:mesh.numBoundaryFaces
      for j=1:(sum(mesh.numNodesPerType[1:2]))
        if i == 1
          @fact mesh.bndry_normals[:, j, i] --> roughly([-1.0, 1.0], atol=1e-13)
        elseif i == 2
          @fact mesh.bndry_normals[:, j, i] --> roughly([1.0, -1.0], atol=1e-13)
        elseif i == 3
          @fact mesh.bndry_normals[:, j, i] --> roughly([0.0, 1.0], atol=1e-13)
        elseif i == 4
          @fact mesh.bndry_normals[:, j, i] --> roughly([1.0, 0.0], atol=1e-13)
        end
      end
    end
=#
#=
    for i=1:mesh.numInterfaces
      for j=1:sbp.numfacenodes
        @fact mesh.interface_normals[:, 1, j, i] --> roughly([1.0, -2.0], atol=1e-13)
        @fact mesh.interface_normals[:, 2, j, i] --> roughly([-2.0, 1.0], atol=1e-13)
      end
    end
=#
    # verify that dofs on a node are numbered consecutively
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        start_dof = mesh.dofs[1, j, i]
        for k=1:mesh.numDofPerNode
          @fact mesh.dofs[k, j, i] --> start_dof + k - 1 
        end
      end
    end

    @fact mesh.color_masks[1][1] --> 1
    @fact mesh.color_masks[1][2] --> 0
    @fact mesh.color_masks[2][1] --> 0
    @fact mesh.color_masks[2][2] --> 1

    @fact mesh.neighbor_colors[1,1] --> 2
    @fact mesh.neighbor_colors[2,1] --> 1
    @fact mesh.neighbor_colors[1,2] --> 1
    @fact mesh.neighbor_colors[2,2] --> 2

    @fact mesh.neighbor_nums[1,1] --> 2
    @fact mesh.neighbor_nums[2,1] --> 1
    @fact mesh.neighbor_nums[2,1] --> 1
    @fact mesh.neighbor_nums[2,2] --> 2

    @fact mesh.pertNeighborEls[1, 1] --> 1
    @fact mesh.pertNeighborEls[2, 1] --> 1
    @fact mesh.pertNeighborEls[1, 2] --> 2
    @fact mesh.pertNeighborEls[2, 2] --> 2

    @fact mesh.color_cnt --> [1,1,0, 0]
    if order == 1
      @fact mesh.numDof --> 16
      @fact mesh.numNodes --> 4
      @fact mesh.numNodesPerElement --> 3
      @fact mesh.numNodesPerType --> [1, 0 , 0]
      @fact mesh.typeOffsetsPerElement --> [1, 4, 4, 4]
    elseif order == 2
      @fact mesh.numNodes --> 11
      @fact mesh.numDof --> 44
      @fact mesh.numNodesPerElement --> 7
      @fact mesh.numNodesPerType --> [1, 1, 1]
      @fact mesh.typeOffsetsPerElement --> [1, 4, 7, 8]
    elseif order == 3
      @fact mesh.numNodes --> 20
      @fact mesh.numDof --> 80
      @fact mesh.numNodesPerElement --> 12
      @fact mesh.numNodesPerType --> [1, 2, 3]
      @fact mesh.typeOffsetsPerElement --> [1, 4, 10, 13]
      # check orientation of 3rd edge
#      println("mesh.dofs = ", mesh.dofs)
      @fact reshape(mesh.dofs[1, 4:5, 1], 2) --> reshape(mesh.dofs[1, 8:9, 2], 2)
    elseif order == 4
      @fact mesh.numNodes --> 31
      @fact mesh.numDof --> 124
      @fact mesh.numNodesPerElement --> 18
      @fact mesh.numNodesPerType --> [1, 3, 6]
      @fact mesh.typeOffsetsPerElement --> [1, 4, 13, 19]
#      println("mesh.dofs = ", mesh.dofs)
      # check orientation of 3rd edge
      @fact reshape(mesh.dofs[1, 10:12, 2], 3) -->reshape(mesh.dofs[1, 4:6, 1], 3)
    end

    @fact mesh.typeOffsetsPerElement_ --> mesh.typeOffsetsPerElement

    @fact length(mesh.verts) --> mesh.numVert
    @fact length(mesh.edges) --> mesh.numEdge
    @fact length(mesh.elements) --> mesh.numEl

#    println("mesh.coords = ", mesh.coords)
    
    @fact mesh.jac --> roughly(ones(mesh.numNodesPerElement ,2))


    fnames = ["boundary_nums", "face_vertnums", "edge_vertnums"]
    for name in fnames
      name_code = string(name, ".dat")
      name_ref = string(name, "_p", order, "true.dat")
      println("checking file ", name_code)
      data_code = readdlm(name_code)
      data_ref = readdlm(name_ref)


      for i=1:length(data_code)
        data_i = data_code[i]
#        println("data_i = ", data_i)
#        println("data_ref[i] = ", data_ref[i])
        if typeof(data_i) <: Number
          @fact data_i --> roughly(data_ref[i], atol=1e-13)
        else
          @fact data_i --> data_ref[i]
        end
      end
    end

  # test vertmap
  println("testing vertmap")

  nedges_interior = 0
  for i=1:mesh.numEl
    vertnums_i = mesh.element_vertnums[:, i]
    for j=1:3
      @fact vertnums_i[j] --> greater_than(0)
      @fact vertnums_i[j] --> less_than(mesh.numVert + 1)
    end

    # test there are exactly 3 elements the current element shares 2 vertices
    # with
    nedges = 0
    for j=1:mesh.numEl
      vertnums_j = mesh.element_vertnums[:, j]
      if length(intersect(vertnums_i, vertnums_j)) == 2
        nedges += 1
        nedges_interior += 1
      end
    end

    @fact nedges --> less_than(4)
    @fact nedges --> greater_than(0)
  end  # end loop i

  @fact nedges_interior --> 2*(mesh.numEdge - mesh.numBoundaryFaces)


  # test getNodeXi
  if order <= 2
    fshape = getFieldShape(0, order, 2)
    eshape = getEntityShape(fshape, 2)  # triangle
    nodexi = PdePumiInterface.getXiCoords(order, 2)
    numnodes = size(nodexi, 2)
    for i=1:numnodes
      vals = getValues(eshape, nodexi[:, i], numnodes)
      for j=1:numnodes
        if i == j
          @fact abs(vals[j] - 1) --> less_than(1e-12)
        else
          @fact abs(vals[j]) --> less_than(1e-12)
        end  # end if else
      end   # end loop j
    end  # end loop i
  end  # end if p < 2




  end  # end loop over p=1:4

  # verify adj and det are correct
  for i=1:10  # 10 random matrices
    A2 = rand(2,2)
    d = det(A2)
    @fact PdePumiInterface.det2(A2) --> roughly(d, atol=1e-12)
    B2 = zeros(A2)
    PdePumiInterface.adjugate2(A2, B2)
    @fact B2./d --> roughly(inv(A2), atol=1e-12)

    A3 = rand(3,3)
    d = det(A3)
    @fact PdePumiInterface.det3(A3) --> roughly(d, atol=1e-12)
    B3 = zero(A3)
    PdePumiInterface.adjugate3(A3, B3)
    @fact B3/d --> roughly(inv(A3), atol=1e-12)
  end

  # now test on a fully unstructured mesh
  smb_name = "vortex.smb"
  dmg_name = "vortex.dmg"

  opts = Dict{Any, Any}(
    "numBC" => 2,
    "BC1" =>  [4, 10],
    "BC2" =>  [7, 13],
    "run_type" => 4,
    "verify_coloring" => true,
    "use_edge_res" => false,
    "write_edge_vertnums" => true,
    "write_face_vertnums" => true,
    "write_boundarynums" => true,
    "write_dxidx" => true,
    "write_coords" => true,
    "write_sparsity" => true,
    "write_offsets" => true,
    "write_counts" => true,
    "write_sparsity_nodebnds" => true,
    "write_dofs" => false,
    "write_interfaces" => false,
    "write_boundaries" => false,
    "write_sharedboundaries" => false,
    )



  for order = 1:4
    println("testing order ", order, " CG mesh against files")
    sbp = TriSBP{Float64}(degree=order)
    ref_verts = [-1. 1 -1; -1 -1 1]
    sbpface = TriFace{Float64}(order, sbp.cub, ref_verts.')


    mesh =  PumiMesh2{Float64}(dmg_name, smb_name, order, sbp, opts, sbpface, coloring_distance=2, dofpernode=4)

    
    for name in fnames
      name_code = string(name, ".dat")
      name_ref = string(name, "vortex", "_p", order, "true.dat")
      println("checking file ", name_code)
      println("against reference file ", name_ref)
      data_code = readdlm(name_code)
      data_ref = readdlm(name_ref)


      for i=1:length(data_code)
        data_i = data_code[i]
        if typeof(data_i) <: Number
          @fact data_i --> roughly(data_ref[i], atol=1e-13)
        else
          @fact data_i --> data_ref[i]
        end
      end
    end
  
  end   # end loop order=1:4


end

facts("----- Testing PdePumiInterfaceDG -----") do

  opts = Dict{Any, Any}(
    "numBC" => 1,
    "BC1" =>  [0],
    "run_type" => 4,
    "verify_coloring" => true,
    "use_edge_res" => false,
    "write_edge_vertnums" => true,
    "write_face_vertnums" => true,
    "write_boundarynums" => true,
    "write_dxidx" => true,
    "write_coords" => true,
    "write_sparsity" => true,
    "write_offsets" => true,
    "write_counts" => true,
    "write_sparsity_nodebnds" => true,
    "write_dofs" => false,
    "write_interfaces" => false,
    "write_boundaries" => false,
    "write_sharedboundaries" => false,
    "exact_visualization" => true,
    )
    order = 1
    smb_name = "tri2l.smb"
    dmg_name = ".null"
#    interp_op = [0.5 0 0; 0 0.5 0; 0 0 0.5]

    sbp = TriSBP{Float64}(degree=order, internal=true)
    vtx = sbp.vtx
#    interp_op = SummationByParts.buildinterpolation(sbp, vtx.')

    sbpface = TriFace{Float64}(order, sbp.cub, vtx)
    println("sbpface.numnodes = ", sbpface.numnodes)
    mesh =  PumiMeshDG2{Float64}(dmg_name, smb_name, order, sbp, opts, sbpface, coloring_distance=2, dofpernode=4)

   @fact mesh.m_ptr --> not(C_NULL)
   @fact mesh.mnew_ptr --> not(C_NULL)
   @fact mesh.numVert --> 4
   @fact mesh.numEdge --> 5
   @fact mesh.numFace --> mesh.numEdge
   @fact mesh.numEl --> 2
   @fact mesh.numDof --> 24
   @fact mesh.numNodes --> 6
   @fact mesh.numDofPerNode --> 4
   @fact mesh.numBoundaryFaces --> 4
   @fact mesh.numInterfaces --> 1
   @fact mesh.numNodesPerElement --> 3
   @fact mesh.numNodesPerType --> [0, 0, 3]
   @fact mesh.numEntitiesPerType --> [4, 5, 2]
   @fact mesh.numTypePerElement --> [3, 3, 1]
   @fact mesh.typeOffsetsPerElement --> [1, 1, 1, 4]
   @fact mesh.typeOffsetsPerElement_ --> mesh.typeOffsetsPerElement
   @fact mesh.volume --> roughly(4.0, atol=1e-12)
   volume2 = PdePumiInterface.calcVolume(mesh)
   @fact volume2 --> roughly(mesh.volume, atol=1e-12)
   @fact mesh.dim --> 2
   @fact mesh.isDG --> true
   @fact mesh.coloringDistance --> 2
   @fact mesh.numColors --> 4
   @fact mesh.numBC --> 1
   @fact mesh.bndry_geo_nums[1] --> opts["BC1"]
   @fact mesh.elementNodeOffsets --> zeros(mesh.numNodesPerElement, mesh.numEl)
   @fact mesh.typeNodeFlags[1] --> trues(3, mesh.numEl)
   tmp = trues(3, 2); tmp[1, 1] = false
   @fact mesh.typeNodeFlags[2] --> tmp
   @fact mesh.bndry_offsets --> [1, 5]
   iface = mesh.interfaces[1]
   @fact iface.elementL --> 1
   @fact iface.elementR --> 2
   @fact iface.faceL --> 1
   @fact iface.faceR --> 3
   @fact mesh.coords[:, :, 1] --> roughly([-2/3 -2/3 1/3; 2/3 -1/3 2/3], atol=1e-13)
   @fact mesh.coords[:, :, 2] --> roughly([2/3 -1/3 2/3; 1/3 -2/3 -2/3], atol=1e-13)
   @fact sort(unique(mesh.dofs)) --> collect(1:mesh.numDof)

   @fact mesh.color_masks[1][1] --> true
   @fact mesh.color_masks[1][2] --> false
   @fact mesh.color_masks[2][1] --> false
   @fact mesh.color_masks[2][2] --> true
   @fact mesh.neighbor_colors[1, 1] --> 2
   @fact mesh.neighbor_colors[2, 1] --> 1
   @fact mesh.neighbor_colors[1, 2] --> 1
   @fact mesh.neighbor_colors[2, 2] --> 2

   @fact mesh.neighbor_nums[1, 1] --> 2
   @fact mesh.neighbor_nums[2, 1] --> 1
   @fact mesh.neighbor_nums[1, 2] --> 1
   @fact mesh.neighbor_nums[2, 2] --> 2

   @fact mesh.pertNeighborEls[1, 1] --> 1
   @fact mesh.pertNeighborEls[2, 1] --> 1
   @fact mesh.pertNeighborEls[1, 2] --> 2
   @fact mesh.pertNeighborEls[2, 2] --> 2

   @fact mesh.color_cnt[1] --> 1
   @fact mesh.color_cnt[2] --> 1

   @fact mesh.coord_order --> 1
   @fact mesh.coord_numNodesPerElement --> 3
   @fact length(mesh.typeOffsetsPerElement) --> 4
   @fact size(mesh.coord_xi, 1) --> 2
   @fact size(mesh.coord_xi, 2) --> 3
     
   # check that adjoint variables are right size
  @fact size(mesh.dxidx) --> size(mesh.dxidx_bar)
  @fact size(mesh.dxidx_bndry) --> size(mesh.dxidx_bndry_bar)
  @fact size(mesh.dxidx_face) --> size(mesh.dxidx_face_bar)
  @fact size(mesh.dxidx_sharedface_bar) --> size(mesh.dxidx_sharedface_bar)


  @fact mesh.jac --> roughly(ones(mesh.numNodesPerElement ,2))

  # check if dxidx is consistent with the old way of calculating it
  dxidx2 = zeros(mesh.dxidx)
  jac2 = zeros(mesh.jac)

  SummationByParts.mappingjacobian!(sbp, mesh.coords, dxidx2, jac2)

  @fact norm(mesh.jac - jac2)/length(mesh.jac) --> roughly(0.0, atol=1e-13)

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      @fact norm(mesh.dxidx[:, :, j, i]  - dxidx2[:, :, j, i]) --> roughly(0.0, atol=1e-13)
    end
  end

  u = ones(mesh.numDofPerNode, mesh.numTypePerElement[1], mesh.numEl)
  PdePumiInterface.saveNodalSolution(mesh, u)
  writeVisFiles(mesh, "nodal_solution")

  # test accumulation at vertices
  u_volume = ones(6, mesh.numTypePerElement[1], mesh.numEl)
  u_verts = ones(6, mesh.numVert)

  PdePumiInterface.accumulateAtVerts(mesh, u_volume, u_verts)

  # check that value at each vert is the number of elements using that vert
  for i=1:mesh.numVert
    vert_i = mesh.verts[i]
    nel = countAdjacent(mesh.m_ptr, vert_i, mesh.dim)
    for j=1:6
      @fact u_verts[j, i] --> nel
    end
  end


  # check reverse mode
  # SBP testing the correctness, these tests only verify values get to the right place

  fill!(mesh.dxidx_bar, 1.0)
  getVertCoords_rev(mesh, sbp)


  for i=1:mesh.numEl
    @fact norm(mesh.vert_coords_bar[:, :, i]) --> greater_than(0.0)
  end

  test_metric_rev(mesh, sbp)

   function test_interp{Tmsh}(mesh::AbstractMesh{Tmsh})
     sbpface = mesh.sbpface
     dxdxi_element = zeros(2, 2, mesh.numNodesPerElement, 1)
     dxdxi_face = zeros(4, sbpface.numnodes, 1)
     dxidx_face = zeros(2,2, sbpface.numnodes)
     jac_face = zeros(sbpface.numnodes)

     for i=1:mesh.numInterfaces
       el = mesh.interfaces[i].elementL
       face = mesh.interfaces[i].elementR

       for j=1:mesh.numNodesPerElement
         dxidx_hat = mesh.dxidx[:, :, j, el]
         dxidx = dxidx_hat*mesh.jac[j, el]
         dxdxi = inv(dxidx)
         dxdxi_element[:, :, j, el] = dxdxi
       end

       bndry_arr = [ Boundary(1, face)]
       dxdxi_element_rshape = reshape(dxdxi_element, 4, mesh.numNodesPerElement, 1)
       boundaryinterpolate!(sbpface, bndry_arr, dxdxi_element_rshape, dxdxi_face)
       dxdxi_face_rshape = reshape(dxdxi_face, 2, 2, sbpface.numnodes)
    for j=1:sbpface.numnodes
      dxidx = inv(dxdxi_face_rshape[:, :, j])
      jac_face[j] = det(dxidx)
      dxidx_face[:, :, j] = dxidx/det(dxidx)

      @fact jac_face[j] --> roughly(mesh.jac_face[j, i], atol=1e-13)
      @fact dxidx_face[:, :, j] --> roughly(mesh.dxidx_face[:, :, j, i], atol=1e-13)
    end

    end  # end loop over interfaces

  end  # end function


  println("testing dxidx_face")
  if size(mesh.dxidx_face, 1) > 0
    println("running the test")
    test_interp(mesh)
  end

   # check dxidx_bndry
   # for straight sided elements, dxidx is constant
   println("testing dxidx_bndry")
   if size(mesh.jac_bndry, 1) > 0
     println("running the test")
     for i=1:mesh.numBoundaryFaces
       el = mesh.bndryfaces[i].element
       dxidx_test = mesh.dxidx[:, :, 1, el]
       jac_test = mesh.jac[1, el]
       for j=1:mesh.sbpface.numnodes
         dxidx_code = mesh.dxidx_bndry[:, :, j, i]
         jac_code = mesh.jac_bndry[j, i]
         @fact dxidx_code --> roughly(dxidx_test, atol=1e-13)
         @fact jac_code --> roughly(jac_test, atol=1e-13)
       end
     end
   end


  println("testing coords bndry")

   # check coords_bndry
   vert_coords = [-1.0 -1.0; 1  1; -1 1]
   
   function getBndry(arr, el, face)
     for i=1:length(arr)
       bndry_i = arr[i]
       if bndry_i.element == el && bndry_i.face == face
         return bndry_i, i
       end
     end
   end

   bndry, idx = getBndry(mesh.bndryfaces, 1, 2)
   # get vertex coordinates
   v1 = vert_coords[2, :]
   v2 = vert_coords[3, :]
   diff = v2 - v1
   slope = diff[2]/diff[1]
   b = v1[2] - slope*v1[1]

   # vertify the face node coordinates are on the line
   for i=1:mesh.sbpface.numnodes
     x_i = mesh.coords_bndry[1, i, idx]
     y_i = mesh.coords_bndry[2, i, idx]
    @fact y_i --> roughly(slope*x_i + b, atol=1e-13)
   end

   function test_normal_orientation{I <: Union{Boundary, Interface}, T}(mesh, ifaces::AbstractArray{I}, nrm::AbstractArray{T, 3})
     println("-----entered test_normal_orientation-----")
     topo = mesh.topo
     numVertPerElement = mesh.numTypePerElement[1]
     numVertPerFace = numVertPerElement - 1 
     println("numVertPerFace = ", numVertPerFace)
     el_verts = Array(Ptr{Void}, numVertPerElement)
     other_vert_coords = zeros(mesh.dim)
     face_verts = Array(Ptr{Void}, numVertPerElement - 1)
     face_vert_coords = zeros(mesh.dim, numVertPerFace)

     println("face_verts = ", topo.face_verts)
     nfaces = length(ifaces)
     for i=1:nfaces
       println("i = ", i)
       iface_i = ifaces[i]
       elnum = getElementL(iface_i)
       facenum_local = getFaceL(iface_i)

       el_i = mesh.elements[elnum]
       nverts = getDownward(mesh.m_ptr, el_i, 0, el_verts)
       println("nverts = ", nverts)

       for j=1:numVertPerFace
         face_verts[j] = el_verts[topo.face_verts[j, facenum_local]]
         tmp = zeros(3)
         getPoint(mesh.m_ptr, face_verts[j], 0, tmp)
         face_vert_coords[1:mesh.dim, j] = tmp[1:mesh.dim]
       end

       # get the vert not on the face
       other_vert = Ptr{Void}(0)
       for j=1:numVertPerElement
         if !(el_verts[j] in face_verts)
           other_vert = el_verts[j]
         end
       end

       tmp = zeros(3)
       getPoint(mesh.m_ptr, other_vert, 0, tmp)
       other_vert_coords[1:mesh.dim] = tmp[1:mesh.dim]

       # check that the face normal is in the opposite direction as the
       # vectors from a vertex on the face to the vertex not on the face
       for j=1:mesh.numNodesPerFace
         for k=1:numVertPerFace
           r1 = other_vert_coords - face_vert_coords[:, k]

           val = dot(nrm[:, j, i], r1)
           @fact val --> less_than(0.0)
         end
       end

     end  # end loop i

     return nothing


   end  # end function test_normal_orientation

   test_normal_orientation(mesh, mesh.interfaces, mesh.nrm_face)
   test_normal_orientation(mesh, mesh.bndryfaces, mesh.nrm_bndry)


   println("testing number nodes windy")
   # check adjacency reordering algorithm doesn't error out
   PdePumiInterface.numberNodesWindy(mesh, [0.0, 0.0, 0.0])

    smb_name = "tri8l.smb"
    mesh =  PumiMeshDG2{Float64}(dmg_name, smb_name, order, sbp, opts, sbpface, coloring_distance=2, dofpernode=4)

  # check mapping interpolation
  # should be constant within an element for straight-sided elementsa
  if size(mesh.jac_face, 1) > 0
    for i=1:mesh.numInterfaces
      iface_i = mesh.interfaces[i]
      el_i = iface_i.elementL
      dxidx_el = mesh.dxidx[:, :, 1, el_i]
      jac_el = mesh.jac[:, el_i]
      jac_face = mesh.jac_face[:, i]

      for j=1:mesh.numNodesPerFace
        dxidx_face = mesh.dxidx_face[:, :, j, i]
        for k=1:2
          for p=1:2
            @fact dxidx_face[p, k] --> roughly(dxidx_el[p, k], atol=1e-13)
          end
        end
        @fact jac_face[j] --> roughly(jac_el[j], atol=1e-13)
      end
    end  # end loop over interfaces
  end

  # test reverse mode interpolation

  test_interp_rev(mesh)
  test_coords_rev(mesh, sbp)

  # test update_coords
  println("testing update_coords")
  coords_orig = copy(mesh.vert_coords)
  
  # double all coordinates
  coords_i = zeros(Float64, 2, 3)
  for i=1:mesh.numEl
    for j=1:3
      coords_i[1, j] = 2*coords_orig[1, j, i]
      coords_i[2, j] = 2*coords_orig[2, j, i]
    end
    update_coords(mesh, i, coords_i)
  end

  commit_coords(mesh, sbp)

#  PdePumiInterface.getCoordinatesAndMetrics(mesh, sbp)
  for i = 1:length(mesh.vert_coords)
    @fact abs(2*coords_orig[i] - mesh.vert_coords[i]) --> less_than(1e-10)
  end

  # test vertmap
  println("testing vertmap")
  nedges_interior = 0
  for i=1:mesh.numEl
    vertnums_i = mesh.element_vertnums[:, i]
    for j=1:3
      @fact vertnums_i[j] --> greater_than(0)
      @fact vertnums_i[j] --> less_than(mesh.numVert + 1)
    end

    # test there are exactly 3 elements the current element shares 2 vertices
    # with
    nedges = 0
    for j=1:mesh.numEl
      vertnums_j = mesh.element_vertnums[:, j]
      if length(intersect(vertnums_i, vertnums_j)) == 2
        nedges += 1
        nedges_interior += 1
      end
    end

    @fact nedges --> less_than(4)
    @fact nedges --> greater_than(0)
  end  # end loop i

  @fact nedges_interior --> 2*(mesh.numEdge - mesh.numBoundaryFaces)

  # test saveSolutionToMesh interpolation
  println("testing saveSolutionToMesh interpolation")
  u_vals = zeros(mesh.numDof)
  for i=1:mesh.numDof

    # get dof, node, element
    idx = findfirst(mesh.dofs, i)
    dofidx, node, el = ind2sub(mesh.dofs, idx)

    x = mesh.coords[1, node, el]
    y = mesh.coords[2, node, el]
    order = mesh.order

    u_vals[i] = x^order + y^order + 1
  end

  saveSolutionToMesh(mesh, u_vals)
  writeVisFiles(mesh, "dg_vis_test")

  # check that the solution is interpolated exactly
  down_verts = Array(Ptr{Void}, mesh.numTypePerElement[1])
  coords_vert = zeros(Float64, 3)
  interp_vals = zeros(mesh.numDofPerNode)
  for i=1:mesh.numEl
    order = mesh.order
    el_ptr = mesh.elements[i]
    getDownward(mesh.m_ptr, el_ptr, 0, down_verts)

    for j=1:mesh.numTypePerElement[1]  # loop over vertices
      vert_j = down_verts[j]
      getPoint(mesh.m_ptr, vert_j, 0, coords_vert)

      x = coords_vert[1]
      y = coords_vert[2]

      # note: this assumes mnew = m
      getComponents(mesh.fnew_ptr, vert_j, 0, interp_vals)

      val_expected = x^order + y^order + 1
      for k=1:mesh.numDofPerNode
        @fact abs(interp_vals[k] - val_expected) --> less_than(1e-12)
      end
    end
  end


  # test periodic
  println("testing periodic")
  @fact mesh.numPeriodicInterfaces --> 0

  smb_name = "tri3_px.smb"
  opts["BC1"] = [0, 2]
  mesh = PumiMeshDG2{Float64}(dmg_name, smb_name, order, sbp, opts, sbpface, coloring_distance=2, dofpernode=4)

  @fact mesh.numPeriodicInterfaces --> 3
  @fact length(mesh.interfaces) --> 24
  @fact mesh.numInterfaces --> 24
  @fact mesh.numBoundaryFaces --> 6

  for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]
    @fact iface_i.elementL --> greater_than(0)
    @fact iface_i.elementL --> less_than(mesh.numEl + 1)
    @fact iface_i.elementR --> greater_than(0)
    @fact iface_i.elementR --> less_than(mesh.numEl + 1)
    @fact iface_i.faceL --> less_than(4)
    @fact iface_i.faceR --> greater_than(0)
    @fact iface_i.faceR --> less_than(4)
  end

  # test curvilinear
  println("testing curvilinear")
  # a 0 - 5 square that used a sin wave to remap the nondimensionalized
  # coordinates
  smb_name = "square_05_curve.smb"
  mesh =  PumiMeshDG2{Float64}(dmg_name, smb_name, order, sbp, opts, sbpface, coloring_distance=2, dofpernode=4)

  function test_volume_curvilinear(mesh, sbp)

    volume = 0.0
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        dxidx_scaled = mesh.dxidx[:, :, j, i]
        dxidx = dxidx_scaled*mesh.jac[j, i]

        dxdxi = inv(dxidx)
        jac = det(dxdxi)
        volume += sbp.w[j]*jac
      end
    end

    println("volume = ", volume)
    @fact volume --> roughly(25.0, atol=1e-12)
  end  # end function

  test_volume_curvilinear(mesh, sbp)


  # translate the mesh and verify all quantities are the same
  dxidx_orig = copy(mesh.dxidx)
  jac_orig = copy(mesh.jac)
  nrm_bndry_orig = copy(mesh.nrm_bndry)
  nrm_face_orig = copy(mesh.nrm_face)

  for i=1:mesh.numEl
    coords_i = mesh.vert_coords[:, :, i]
    coords_i += 1
    update_coords(mesh, i, coords_i)
  end

  commit_coords(mesh, sbp)

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      for k=1:mesh.dim
        for p=1:mesh.dim
          @fact mesh.dxidx[p, k, j, i] --> roughly(dxidx_orig[p, k, j, i], atol=1e-12)
        end
      end

      @fact mesh.jac[j, i] --> roughly(jac_orig[j, i], atol=1e-12)
    end
  end

  for i=1:mesh.numBoundaryFaces
    for j=1:mesh.numNodesPerFace
      for k=1:mesh.dim
        @fact mesh.nrm_bndry[k, j, i] --> roughly(nrm_bndry_orig[k, j, i], atol=1e-12)
      end
    end
  end

  for i=1:mesh.numInterfaces
    for j=1:mesh.numNodesPerFace
      for k=1:mesh.dim
        @fact mesh.nrm_face[k, j, i] --> roughly(nrm_face_orig[k, j, i], atol=1e-12)
      end
    end
  end
      

  println("finished")

end
