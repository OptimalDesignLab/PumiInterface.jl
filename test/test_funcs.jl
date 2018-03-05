# test functions, potentially used by 2D and 3D tests


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
  Test the reverse mode of the curvilinear metric calculation.

  This function tests the reverse mode of calcEoneElement and then
  calls other functions to test the other parts of the calculation.
"""
function test_metric_rev(mesh, mesh_c, sbp)

  facts("----- testing metric reverse mode -----") do
    sbpface = mesh.sbpface
    nrm = rand(Complex128, mesh.dim, mesh.numNodesPerFace)
    nrm[:, :] = real(nrm)
    Rone = ones(mesh.numNodesPerFace)
    tmp = zeros(Complex128, mesh.numNodesPerFace)

    # test calcEoneElement
    
    Eone_el_c = zeros(Complex128, size(sbpface.perm, 1), mesh.dim)
    
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


    test_metric2_rev(mesh, mesh_c, sbp)
    test_metrics3_rev(mesh, mesh_c, sbp)
    test_metrics4_rev(mesh, mesh_c, sbp)
    test_metrics5_rev(mesh, mesh_c, sbp)

  end  # end facts block

  return nothing
end

"""
  Test calcEone_rev
"""
function test_metric2_rev(mesh, mesh_c, sbp)

  facts("----- Testing second part of Metric reverse mode -----") do

    # test calcEone_rev
    # use finite differences because the perturbation is applied to the fields
    # of the mesh

    h = 1e-20
    pert = Complex128(0, h)

    nin = length(mesh.nrm_face) + length(mesh.nrm_bndry)
    for i=1:mesh.npeers
      nin += length(mesh.nrm_sharedface[i])
    end
    Eone = zeros(Complex128, mesh.numNodesPerElement, mesh.dim, mesh.numEl)
    nout = length(Eone)
#    nout = mesh.numNodesPerElement*mesh.dim*mesh.numEl

#    Eone_orig = zeros(Eone)

    jac = zeros(nin, nout)
    jac2 = zeros(jac)
    element_range = 1:mesh.numEl
 #   PdePumiInterface.calcEone(mesh, sbp, element_range, Eone_orig)


    println("forward mode")
    in_idx = 1
    for i=1:length(mesh_c.nrm_face)
      mesh_c.nrm_face[i] += pert
      fill!(Eone, 0.0)
      PdePumiInterface.calcEone(mesh_c, sbp, element_range, Eone)

      for j=1:nout
        jac[in_idx, j] = imag(Eone[j])/h  #(Eone[j] - Eone_orig[j])/pert
      end

      in_idx += 1
      mesh_c.nrm_face[i] -= pert
    end

    for i=1:length(mesh_c.nrm_bndry)
      mesh_c.nrm_bndry[i] += pert
      fill!(Eone, 0.0)
      PdePumiInterface.calcEone(mesh_c, sbp, element_range, Eone)
      for j=1:nout
        jac[in_idx, j] = imag(Eone[j])/h #(Eone[j] - Eone_orig[j])/pert
      end

      in_idx += 1
      mesh_c.nrm_bndry[i] -= pert
    end

    for peer=1:mesh_c.npeers
      nrm_peer = mesh_c.nrm_sharedface[peer]
      for i=1:length(nrm_peer)
        nrm_peer[i] += pert
        fill!(Eone, 0.0)
        PdePumiInterface.calcEone(mesh_c, sbp, element_range, Eone)
        for j=1:nout
          jac[in_idx, j] = imag(Eone[j])/h #(Eone[j] - Eone_orig[j])/pert
        end

        in_idx += 1
        nrm_peer[i] -= pert
      end
    end
    

    # reverse mode
    println("reverse mode")
    Eone_bar = zeros(mesh.numNodesPerElement, mesh.dim, mesh.numEl)
    for j=1:nout
      Eone_bar[j] = 1

      fill!(mesh.nrm_face_bar, 0.0)
      fill!(mesh.nrm_bndry_bar, 0.0)
      for i=1:mesh.npeers
        fill!(mesh.nrm_sharedface_bar[i], 0.0)
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

    @fact maximum(abs(jac - jac2)) --> roughly(0.0, atol=1e-12)
  end  # end facts block

  return nothing
end

"""
  Test getCurvilinearCoordinatesAndMetrics_rev (dxidx -> vert_coord, nrm)
"""
function test_metrics3_rev(mesh, mesh_c, sbp)

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

    h = 1e-20
    pert = Complex128(0, h)

#    dxidx_orig = copy(mesh.dxidx)
#    jac_orig = copy(mesh.jac)
#    for i=1:mesh.npeers
#      nrm_sharedface[i] = copy(mesh.nrm_sharedface[i])
#    end

    # forward mode
    println("forward mode verts")
    in_idx = 1
    for i=1:length(mesh_c.vert_coords)
      mesh_c.vert_coords[i] += pert

      PdePumiInterface.getCurvilinearCoordinatesAndMetrics(mesh_c, sbp)

      out_idx = 1
      for j=1:length(mesh_c.dxidx)
        jac[in_idx, out_idx] = imag(mesh_c.dxidx[j])/h  #(mesh_c.dxidx[j] - dxidx_orig[j])/pert
        out_idx += 1
      end

      # TODO: uncomment when this is fixed
      for j=1:length(mesh_c.jac)
#        jac[in_idx, out_idx] = imag(mesh_c.jac[j])/h (mesh_c.jac[j] - jac_orig[j])/pert
        out_idx += 1
      end

      in_idx += 1
      mesh_c.vert_coords[i] -= pert
    end

    for i=1:length(mesh_c.nrm_face)
      mesh_c.nrm_face[i] += pert
      PdePumiInterface.getCurvilinearCoordinatesAndMetrics(mesh_c, sbp)

      out_idx = 1
      for j=1:length(mesh_c.dxidx)
        jac[in_idx, out_idx] = imag(mesh_c.dxidx[j])/h #(mesh_c.dxidx[j] - dxidx_orig[j])/pert
        out_idx += 1
      end

      for j=1:length(mesh_c.jac)
        jac[in_idx, out_idx] = imag(mesh_c.jac[j])/h  #(mesh_c.jac[j] - jac_orig[j])/pert
        out_idx += 1
      end

      in_idx += 1
      mesh_c.nrm_face[i] -= pert
    end

    for i=1:length(mesh_c.nrm_bndry)
      mesh_c.nrm_bndry[i] += pert
      PdePumiInterface.getCurvilinearCoordinatesAndMetrics(mesh_c, sbp)

      out_idx = 1
      for j=1:length(mesh_c.dxidx)
        jac[in_idx, out_idx] = imag(mesh_c.dxidx[j])/h  #(mesh_c.dxidx[j] - dxidx_orig[j])/pert
        out_idx += 1
      end

      # TODO: uncomment when this is fixed
      for j=1:length(mesh_c.jac)
#        jac[in_idx, out_idx] = imag(mesh_c.jac[j])/h #(mesh_c.jac[j] - jac_orig[j])/pert
        out_idx += 1
      end

      in_idx += 1
      mesh_c.nrm_bndry[i] -= pert
    end

    for peer=1:mesh_c.npeers
      println("peer ", peer)
      nrm_peer = mesh_c.nrm_sharedface[peer]
      for i=1:length(nrm_peer)
        nrm_peer[i] += pert
        PdePumiInterface.getCurvilinearCoordinatesAndMetrics(mesh_c, sbp)

        out_idx = 1
        for j=1:length(mesh_c.dxidx)
          jac[in_idx, out_idx] = imag(mesh_c.dxidx[j])/h #(mesh_c.dxidx[j] - dxidx_orig[j])/pert
          out_idx += 1
        end

        # uncomment when this is fixed
        for j=1:length(mesh_c.jac)
#          jac[in_idx, out_idx] = imag(mesh_c.dxidx[j])/h  #(mesh_c.jac[j] - jac_orig[j])/pert
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
        for i=1:length(nrm_peer_bar)
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
        for i=1:length(nrm_peer_bar)
          jac2[in_idx, out_idx] = nrm_peer_bar[i]
          in_idx += 1
        end
      end
      
      out_idx += 1
      mesh.jac_bar[j] = 0
    end  # end loop j


    @fact maximum(abs(jac - jac2)) --> roughly(0.0, atol=1e-5)
  end

  return nothing
end


function test_metrics4_rev(mesh, mesh_c, sbp)

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

    h = 1e-20
    pert = Complex128(0, h)

    # forward mode
    in_idx = 1
    for i=1:length(mesh_c.vert_coords)
      mesh_c.vert_coords[i] += pert

      PdePumiInterface.getFaceCoordinatesAndNormals(mesh_c, sbp)

      out_idx = 1

      for j=1:length(mesh_c.nrm_bndry)
        jac[in_idx, out_idx] = imag(mesh_c.nrm_bndry[j])/h  #(mesh_c.nrm_bndry[j] - nrm_bndry_orig[j])/pert
        out_idx += 1
      end

      for j=1:length(mesh_c.nrm_face)
        jac[in_idx, out_idx] = imag(mesh_c.nrm_face[j])/h #(mesh_c.nrm_face[j] - nrm_face_orig[j])/pert
        out_idx += 1
      end

      for peer=1:mesh_c.npeers
        for j=1:length(mesh_c.nrm_sharedface[peer])
          jac[in_idx, out_idx] = imag(mesh_c.nrm_sharedface[peer][j])/h #(mesh_c.nrm_sharedface[peer][j] - nrm_sharedface_orig[peer][j])/pert
          out_idx += 1
        end
      end


      in_idx += 1
      mesh_c.vert_coords[i] -= pert
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

     
    diffnorm = maximum(abs(jac - jac2))
    @fact diffnorm --> roughly(0.0, atol=1e-5)
  end

  return nothing
end


"""
  Test everything together
"""
function test_metrics5_rev(mesh, mesh_c, sbp)
  
  facts("----- testing metrics5_rev -----") do

    nin = length(mesh.vert_coords)

    nout = length(mesh.dxidx) + length(mesh.jac) + length(mesh.nrm_face) + length(mesh.nrm_bndry)
    for i=1:mesh.npeers
      nout += length(mesh.nrm_sharedface[i])
    end


    h = 1e-20
    pert = Complex128(0, h)

    jac = zeros(nin, nout)
    jac2 = zeros(jac)


    println("forward mode")
    for i=1:length(mesh_c.vert_coords)
      mesh_c.vert_coords[i] += pert
      PdePumiInterface.getFaceCoordinatesAndNormals(mesh_c, sbp)
      PdePumiInterface.getCurvilinearCoordinatesAndMetrics(mesh_c, sbp)

      out_idx = 1
      for j=1:length(mesh_c.dxidx)
        jac[i, out_idx] = imag(mesh_c.dxidx[j])/h
        out_idx += 1
      end

      
      for j=1:length(mesh_c.jac)
  #      jac[i, out_idx] = imag(mesh_c.jac[j])/h
        out_idx += 1
      end

      for j=1:length(mesh_c.nrm_face)
        jac[i, out_idx] = imag(mesh_c.nrm_face[j])/h
        out_idx += 1
      end

      for j=1:length(mesh_c.nrm_bndry)
        jac[i, out_idx] = imag(mesh_c.nrm_bndry[j])/h
        out_idx += 1
      end

      for peer=1:mesh_c.npeers
        nrm_peer = mesh_c.nrm_sharedface[peer]
        for j=1:length(nrm_peer)
          jac[i, out_idx] = imag(nrm_peer[j])/h
          out_idx += 1
        end
      end

      mesh_c.vert_coords[i] -= pert
    end


    println("reverse mode")
    zeroBarArrays(mesh)
    out_idx = 1
    for j=1:length(mesh.dxidx)
      mesh.dxidx_bar[j] = 1

      PdePumiInterface.getCurvilinearCoordinatesAndMetrics_rev(mesh, sbp)
      PdePumiInterface.getFaceCoordinatesAndNormals_rev(mesh, sbp)

      for i=1:length(mesh.vert_coords)
        jac2[i, out_idx] = mesh.vert_coords_bar[i]
      end

      out_idx += 1
      zeroBarArrays(mesh)
    end

    for j=1:length(mesh.jac)
      mesh.jac_bar[j] = 0  # 1

#      getAllCoordinatesAndMetrics_rev(mesh, sbp)
      PdePumiInterface.getCurvilinearCoordinatesAndMetrics_rev(mesh, sbp)
      PdePumiInterface.getFaceCoordinatesAndNormals_rev(mesh, sbp)


      for i=1:length(mesh.vert_coords)
        jac2[i, out_idx] = mesh.vert_coords_bar[i]
      end

      out_idx += 1
      zeroBarArrays(mesh)
    end

    for j=1:length(mesh.nrm_face)
      mesh.nrm_face_bar[j] = 1

#      getAllCoordinatesAndMetrics_rev(mesh, sbp)
      PdePumiInterface.getCurvilinearCoordinatesAndMetrics_rev(mesh, sbp)
      PdePumiInterface.getFaceCoordinatesAndNormals_rev(mesh, sbp)


      for i=1:length(mesh.vert_coords)
        jac2[i, out_idx] = mesh.vert_coords_bar[i]
      end

      out_idx += 1
      zeroBarArrays(mesh)
    end

    for j=1:length(mesh.nrm_bndry)
      mesh.nrm_bndry_bar[j] = 1

#      getAllCoordinatesAndMetrics_rev(mesh, sbp)
      PdePumiInterface.getCurvilinearCoordinatesAndMetrics_rev(mesh, sbp)
      PdePumiInterface.getFaceCoordinatesAndNormals_rev(mesh, sbp)


      for i=1:length(mesh.vert_coords)
        jac2[i, out_idx] = mesh.vert_coords_bar[i]
      end

      out_idx += 1
      zeroBarArrays(mesh)
    end


    for peer=1:mesh.npeers
      nrm_bar_peer = mesh.nrm_sharedface_bar[peer]
      for j=1:length(nrm_bar_peer)
        nrm_bar_peer[j] = 1

#        getAllCoordinatesAndMetrics_rev(mesh, sbp)
        PdePumiInterface.getCurvilinearCoordinatesAndMetrics_rev(mesh, sbp)
        PdePumiInterface.getFaceCoordinatesAndNormals_rev(mesh, sbp)

        for i=1:length(mesh.vert_coords)
          jac2[i, out_idx] = mesh.vert_coords_bar[i]
        end

        out_idx += 1
        zeroBarArrays(mesh)
      end
    end

    @fact maximum(abs(jac - jac2)) --> roughly(0.0, atol=1e-12)
  end

  return nothing
end




function test_coords_rev(mesh, sbp)

#  sbp_complex = TriSBP{Complex128}(degree=sbp.degree, reorder=sbp.reorder, internal=sbp.internal, vertices=sbp.vertices)

#  function coords_forward(vert_coords)
#    coords_it = vert_coords.'
#    coords_volume = SummationByParts.SymCubatures.calcnodes(sbp.cub, coords_it)
#    return coords_volume
# : end
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

function test_submesh()

  facts("----- Testing SubMesh -----") do
    order = 1
    sbp = getTriSBPOmega(degree=order)
    vtx = sbp.vtx
    sbpface = TriFace{Float64}(order, sbp.cub, vtx)

    opts = PdePumiInterface.get_defaults()
    opts["order"] = order
    opts["coloring_distance"] = 2
    opts["dmg_name"] = ".null"
    opts["smb_name"] = "tri8l.smb"
    opts["numBC"] = 1
    opts["BC1"] = [0, 1, 2, 3]
    opts["BC1_name"] = "testBC"

    mesh = PumiMeshDG2(Float64, sbp, opts, sbpface, dofpernode=4)

    el_list = Cint[1, 2, 3, 4]

    submesh, subopts = PumiMeshDG2(mesh, sbp, opts, "resolveBC", el_list)

    @fact submesh.order --> mesh.order
    @fact submesh.numEl --> length(el_list)
    @fact subopts["numBC"] --> 2
    @fact subopts["BC2_name"] --> "resolveBC"
    @fact (subopts["BC2"][1] in subopts["BC1"]) --> false

    # do a a dangling mesh
    el_list = Cint[1, 2, 5, 6]

    submesh, subopts = PumiMeshDG2(mesh, sbp, opts, "resolveBC", el_list)

    @fact submesh.order --> mesh.order
    @fact submesh.numEl --> length(el_list)
    @fact subopts["numBC"] --> 2
    @fact subopts["BC2_name"] --> "resolveBC"
    @fact (subopts["BC2"][1] in subopts["BC1"]) --> false


    # test injection and rejection
    qold = zeros(Float64, mesh.numDof)
    qold2 = zeros(Float64, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        for k=1:mesh.numDofPerNode
          dof = mesh.dofs[k, j, i]
          qold[dof] = i + j + k
          qold2[k, j, i] = i + j + k
        end
      end
    end
    qnew = zeros(Float64, submesh.numDof)
    qnew2 = zeros(Float64, submesh.numDofPerNode, submesh.numNodesPerElement, submesh.numEl)

    injectionOperator(mesh, qold, submesh, qnew)
    injectionOperator(mesh, qold2, submesh, qnew2)

    for i=1:submesh.numDof
      qnew[i] += 1
      qnew2[i] += 1
    end

    qold3 = zeros(qold)
    qold4 = zeros(qold2)
    rejectionOperator(submesh, qnew, mesh, qold3)
    rejectionOperator(submesh, qnew2, mesh, qold4)

    # only the submesh elements got incremented
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        for k=1:mesh.numDofPerNode
          dof = mesh.dofs[k, j, i]
          if i in el_list
            @fact qold3[dof] --> i + j + k + 1
            @fact qold4[k, j, i] --> i + j + k + 1
          else
            @fact qold3[dof] --> 0.0
            @fact qold4[k, j, i] --> 0.0
          end
        end
      end
    end

    # test boundary interpolation array
    interp_arr = getBoundaryInterpArray(mesh, submesh)

    for iface in interp_arr
      @fact iface.elementL in el_list --> true
      @fact iface.elementR in el_list --> false
      @fact iface.elementL == iface.elementR --> false
    end


  end  # end facts


end


function test_parallel_metrics(mesh, sbp, opts)

  # put known values into mesh arrays, do data exchange, check the values
  # came out correctly

  # need to skip any element that is on multiple boundaries
  skip_els = Set{Int}()
  for i=1:mesh.npeers
    el_list_i = mesh.local_element_lists[i]

    for j=1:(mesh.npeers-1)
      el_list_j = mesh.local_element_lists[i]
      common_els = intersect(el_list_i, el_list_j)
      for k=1:length(common_els)
        push!(skip_els, common_els[k])
      end
    end
  end
  # do one peer at a time in case a single element is shared with multiple peers
  facts("----- Testing exchangeMetricInfo -----") do 
    for peer=1:mesh.npeers
      # vert_coords
      numel = mesh.local_element_counts[peer]
      rank = mesh.myrank

      for i=1:numel
        elnum = mesh.local_element_lists[peer]

        if elnum in skip_els
          skip_i = true
        else
          skip_i = false
        end

        val = 1 + rank
        for j=1:mesh.coord_numNodesPerElement
          for k=1:mesh.dim
            if skip_i
              tmp = -1
            else
              tmp = val
            end
            mesh.vert_coords[k, j, elnum] = tmp
            val += 1
          end
        end

        val = 2 + rank
        for j=1:mesh.numNodesPerElement
          for k=1:mesh.dim
            if skip_i
              tmp = -1
            else
              tmp = val
            end
            mesh.coords[k, j, elnum] = tmp
            val += 1
          end
        end

        val = 3 + rank
        for j=1:mesh.numNodesPerElement
          for d2=1:mesh.dim
            for d1=1:mesh.dim
              if skip_i
                tmp = -1
              else
                tmp = val
              end
              mesh.dxidx[d1, d2, j, elnum] = tmp
              val += 1
            end
          end
        end

        val = 4 + rank
        for j=1:mesh.numNodesPerElement
          if skip_i
            tmp = -1
          else
            tmp = val
          end
          mesh.jac[j, elnum] = tmp
          val += 1
        end
      end  # end loop i

    end  # end loop peer

    # do the data exchange
    PdePumiInterface.exchangeMetricInfo(mesh, sbp)

    # check that is came out right on the other end
    for peer=1:mesh.npeers
      numel = mesh.remote_element_counts[peer]
      obj = mesh.remote_metrics[peer]
      rank = mesh.peer_parts[peer]

      for i=1:numel

        val = 1 + rank
        for j=1:mesh.coord_numNodesPerElement
          for k=1:mesh.dim
            tmp = obj.vert_coords[k, j, i]
            if tmp >= 0 
              @fact tmp --> val
            end
            val += 1
          end
        end

        val = 2 + rank
        for j=1:mesh.numNodesPerElement
          for k=1:mesh.dim
            tmp = obj.coords[k, j, i]
            if tmp >= 0
              @fact tmp --> val
            end
            val += 1
          end
        end

        val = 3 + rank
        for j=1:mesh.numNodesPerElement
          for d2=1:mesh.dim
            for d1=1:mesh.dim
              tmp = obj.dxidx[d1, d2, j, i]
              if tmp >= 0
                @fact tmp --> val
              end
              val += 1
            end
          end
        end

        val = 4 + rank
        for j=1:mesh.numNodesPerElement
          tmp = obj.jac[j, i]
          if tmp >= 0
            @fact tmp --> val
          end
          val += 1
        end

      end  # end loop i

    end  # end loop peer

  end  # end facts block

  return nothing
end

function testSurfaceNumbering(mesh, sbp, opts)

  facts("----- Testing Surface Numbering -----") do
    nBCs = opts["numBC"]
    bc_nums =  [nBCs]  # only do the last BC
    numFacePts, n_face, face_verts = numberSurfacePoints(mesh, 1:nBCs)

    @fact length(face_verts) --> numFacePts

    geo_tags = opts["BC$nBCs"]
    
    # check that every face verts is on the specified geometry
    for vert in face_verts
      me =  toModel(mesh.m_ptr, vert)
      me_type = getModelType(mesh.m_ptr, me)  # dimension of model entity
      me_tag = getModelTag(mesh.m_ptr, me)

      @fact me_type --> mesh.dim-1  # is on a face somewhere
      @fact me_tag in geo_tags --> true
    end

    # check that all verts with number numFacePts + 1 are not on the geometry
    for vert in mesh.verts
      n_v = getNumberJ(n_face, vert, 0, 0)
      me =  toModel(mesh.m_ptr, vert)
      me_type = getModelType(mesh.m_ptr, me)
      me_tag = getModelTag(mesh.m_ptr, me)

      if n_v > numFacePts
        if me_type == mesh.dim - 1
          @fact !(me_tag in geo_tags) --> true
        end
      else
        if me_type == mesh.dim - 1
          @fact me_tag in geo_tags --> true
        end
      end
    end
  end  # end facts block

  return nothing
end
