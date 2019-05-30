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
  @testset "\n ----- Testing interp rev -----" begin

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


    @test isapprox( norm(jac_forward - jac_rev), 0.0) atol=1e-11

    # now test the full routine
    fill!(mesh.dxidx_bar, 0.0)
#    mesh.dxidx_face_bar[1] = 1
    fill!(mesh.dxidx_face_bar, 1.0)

    # check regular interfaces
    PdePumiInterface.interpolateMapping_rev(mesh)
    for i=1:mesh.numInterfaces
      el = mesh.interfaces[i].elementL
      for j=1:mesh.numNodesPerFace
        @test  norm(mesh.dxidx_bar[:, :, j, el])  > 0.0
      end
    end
    
    # check boundary faces
    fill!(mesh.dxidx_bar, 0.0)
    fill!(mesh.dxidx_bndry_bar, 1.0)
    PdePumiInterface.interpolateMapping_rev(mesh)
    for i=1:mesh.numBoundaryFaces
      el = mesh.bndryfaces[i].element
      for j=1:mesh.numNodesPerFace
        @test  norm(mesh.dxidx_bar[:, :, j, el])  > 0.0
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
          @test  norm(mesh.dxidx_bar[:, :, j, el])  > 0.0
        end
      end

      fill!(mesh.dxidx_bar, 0.0)
      fill!(mesh.dxidx_sharedface[peer], 0.0)
    end
   
  end

end  # end function

"""
  Test the reverse mode of the curvilinear metric calculation.

  This function tests the reverse mode of calcEoneElement and then
  calls other functions to test the other parts of the calculation.
"""
function test_metrics_rev(mesh, mesh_c, sbp, opts)

  @testset "----- testing metric reverse mode -----" begin
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

    @test isapprox( norm(jac - jac2)/length(jac2), 0.0) atol=1e-12


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

    @test isapprox( norm(jac - jac2), 0.0) atol=1e-12


    test_metric2_rev(mesh, mesh_c, sbp, opts)
    test_metrics3_rev(mesh, mesh_c, sbp, opts)
    test_metrics4_rev(mesh, mesh_c, sbp, opts)
    test_metrics5_rev(mesh, mesh_c, sbp, opts)


  end  # end facts block

  return nothing
end

"""
  Test calcEone_rev
"""
function test_metric2_rev(mesh, mesh_c, sbp, opts)

  @testset "----- Testing second part of Metric reverse mode -----" begin

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

    @test isapprox( maximum(abs.(jac - jac2)), 0.0) atol=1e-12
  end  # end facts block

  return nothing
end

"""
  Test getCurvilinearCoordinatesAndMetrics_rev (dxidx -> vert_coord, nrm)
"""
function test_metrics3_rev(mesh, mesh_c, sbp, opts)

  @testset "----- testing metrics reverse mode 3 -----" begin
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
        #jac[in_idx, out_idx] = imag(mesh_c.jac[j])/h #(mesh_c.jac[j] - jac_orig[j])/pert
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
        #jac[in_idx, out_idx] = imag(mesh_c.jac[j])/h  #(mesh_c.jac[j] - jac_orig[j])/pert
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
        #jac[in_idx, out_idx] = imag(mesh_c.jac[j])/h #(mesh_c.jac[j] - jac_orig[j])/pert
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
          #jac[in_idx, out_idx] = imag(mesh_c.jac[j])/h  #(mesh_c.jac[j] - jac_orig[j])/pert
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
      #mesh.jac_bar[j] = 1

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


    println("max diff = ", maximum(abs.(jac - jac2)))
    @test isapprox( maximum(abs.(jac - jac2)), 0.0) atol=1e-13
  end

  return nothing
end


function test_metrics4_rev(mesh, mesh_c, sbp, opts)

  # zero out all bar variables, just in case
  fill!(mesh.dxidx_bar, 0.0)
  fill!(mesh.jac_bar, 0.0)
  fill!(mesh.vert_coords_bar, 0.0)
  fill!(mesh.nrm_face_bar, 0.0)
  fill!(mesh.nrm_bndry_bar, 0.0)
  for i=1:mesh.npeers
    fill!(mesh.nrm_sharedface_bar[i], 0.0)
  end

  @testset "----- Testing metrics4_rev -----" begin
    nin = length(mesh.vert_coords)
    nout = length(mesh.nrm_bndry) + length(mesh.nrm_face) + length(mesh.coords_bndry)
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

      for j=1:length(mesh_c.coords_bndry)
        jac[in_idx, out_idx] = imag(mesh_c.coords_bndry[j])/h
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

    for j=1:length(mesh.coords_bndry)
      mesh.coords_bndry_bar[j] = 1
      fill!(mesh.vert_coords_bar, 0.0)

      PdePumiInterface.getFaceCoordinatesAndNormals_rev(mesh, sbp)
      in_idx = 1
      for i=1:length(mesh.vert_coords_bar)
        jac2[in_idx, out_idx] = mesh.vert_coords_bar[i]
        in_idx += 1
      end

      out_idx += 1 
      mesh.coords_bndry_bar[j] = 0
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

     
    diffnorm = maximum(abs.(jac - jac2))
    println("diffnorm = ", diffnorm)
    @test isapprox( diffnorm, 0.0) atol=1e-13
  end

  return nothing
end


"""
  Test everything together
"""
function test_metrics5_rev(mesh, mesh_c, sbp, opts)
  
  @testset "----- testing metrics5_rev -----" begin

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
        jac[i, out_idx] = imag(mesh_c.jac[j])/h
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
      mesh.jac_bar[j] = 1  #TODO: uncomment when this is fixed

      #TODO: should these be in the other order?
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

    @test isapprox( maximum(abs.(jac - jac2)), 0.0) atol=1e-12
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

  @testset "----- Testing coordiate reverse mode -----" begin
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
      @test isapprox( jac[j, i], jac2[j, i]) atol=1e-13
    end
  end


  end

  return nothing
end

function test_submesh()

  @testset "----- Testing SubMesh -----" begin
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

    @test ( submesh.order )== mesh.order
    @test ( submesh.numEl )== length(el_list)
    @test ( subopts["numBC"] )== 2
    @test ( subopts["BC2_name"] )== "resolveBC"
    @test ( (subopts["BC2"][1] in subopts["BC1"]) )== false

    # do a a dangling mesh
    el_list = Cint[1, 2, 5, 6]

    submesh, subopts = PumiMeshDG2(mesh, sbp, opts, "resolveBC", el_list)

    @test ( submesh.order )== mesh.order
    @test ( submesh.numEl )== length(el_list)
    @test ( subopts["numBC"] )== 2
    @test ( subopts["BC2_name"] )== "resolveBC"
    @test ( (subopts["BC2"][1] in subopts["BC1"]) )== false


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
            @test ( qold3[dof] )== i + j + k + 1
            @test ( qold4[k, j, i] )== i + j + k + 1
          else
            @test ( qold3[dof] )== 0.0
            @test ( qold4[k, j, i] )== 0.0
          end
        end
      end
    end

    # test boundary interpolation array
    interp_arr = getBoundaryInterpArray(mesh, submesh)

    for iface in interp_arr
      @test ( iface.elementL in el_list )== true
      @test ( iface.elementR in el_list )== false
      @test ( iface.elementL == iface.elementR )== false
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
  @testset "----- Testing exchangeMetricInfo -----" begin
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
              @test ( tmp )== val
            end
            val += 1
          end
        end

        val = 2 + rank
        for j=1:mesh.numNodesPerElement
          for k=1:mesh.dim
            tmp = obj.coords[k, j, i]
            if tmp >= 0
              @test ( tmp )== val
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
                @test ( tmp )== val
              end
              val += 1
            end
          end
        end

        val = 4 + rank
        for j=1:mesh.numNodesPerElement
          tmp = obj.jac[j, i]
          if tmp >= 0
            @test ( tmp )== val
          end
          val += 1
        end

      end  # end loop i

    end  # end loop peer

  end  # end facts block

  return nothing
end

function testSurfaceNumbering(mesh, sbp, opts)

  @testset "----- Testing Surface Numbering -----" begin
    nBCs = opts["numBC"]
    bc_nums =  [nBCs]  # only do the last BC  #TODO: re-enable this
#    bc_nums = 1:nBCs
    numFacePts, n_face, face_verts = numberSurfacePoints(mesh, bc_nums)

    @test ( length(face_verts) )== numFacePts

    geo_tags = Array{Int}(0)
    for bc in bc_nums
      geo_tags_bc = opts["BC$bc"]
      for j=1:length(geo_tags_bc)
        push!(geo_tags, geo_tags_bc[j])
      end
    end

    
    # check that if the MeshEntity is classified in an edge/face, then it
    # is part of the geometry
    # this is a bit too loose because some MeshEntities are classified on
    # the geometric entities that are on the closure of the specified geometry
    for vert in face_verts
      me =  apf.toModel(mesh.m_ptr, vert)
      me_type = apf.getModelType(mesh.m_ptr, me)  # dimension of model entity
      me_tag = apf.getModelTag(mesh.m_ptr, me)

      if me_type == mesh.dim - 1
        @test ( me_tag in geo_tags )== true
      end

    end

    # check that all entities with number numFacePts + 1 are not on the geometry

    seen_nums = zeros(Int, numFacePts)
    entities = Array{Vector{Ptr{Void}}}(2)
    entities[1] = mesh.verts
    entities[2] = mesh.edges
    for edim = 1:length(entities)
      if apf.hasNodesIn(mesh.coordshape_ptr, edim - 1)
        for vert in entities[edim]
          n_v = apf.getNumberJ(n_face, vert, 0, 0)
          me =  apf.toModel(mesh.m_ptr, vert)
          me_type = apf.getModelType(mesh.m_ptr, me)
          me_tag = apf.getModelTag(mesh.m_ptr, me)

          if n_v > numFacePts
            if me_type == mesh.dim - 1
              @test ( !(me_tag in geo_tags) )== true
            end
          else 
            seen_nums[n_v] += 1  # for uniqueness check
            if me_type == mesh.dim - 1
              @test ( me_tag in geo_tags )== true
            end
          end
        end  # end loop vert
      end  # end if
    end  # end loop edim

    @test ( maximum(seen_nums) )== 1
    @test ( minimum(seen_nums) )== 1
  end  # end facts block

  return nothing
end

"""
  Do sanity checks on the coordinate field number.  Currently tests

   * all coordinate field dofs have numbers assigned
   * all the numbers are in the range 1 to mesh.dim * coord_numNodesPerElement
   * all numbers from 1 to coord_numNodesPerElement exist exactly once
"""
function test_coordNumbering(mesh)

  @testset "Testing coordinate field Numbering" begin

    nums = zeros(mesh.dim*mesh.coord_numNodes)
    # check verts
    for entity in mesh.verts
      for i=1:mesh.dim
        val = apf.getNumberJ(mesh.coord_nodenums_Nptr, entity, 0, i-1)
        @test val >= 1
        @test val <= mesh.dim*mesh.coord_numNodes

        nums[val] += 1
      end
    end

    if mesh.coord_order == 2
      for entity in mesh.edges
        for i=1:mesh.dim
          val = apf.getNumberJ(mesh.coord_nodenums_Nptr, entity, 0, i-1)
          @test val >= 1
          @test val <= mesh.dim*mesh.coord_numNodes

          nums[val] += 1
        end
      end
    end

    # every number exists once
    for val in nums
      @test val == 1
    end

  end  # end testset

  return nothing
end


"""
  Test coordinate field transformation functions
"""
function test_coord_field(mesh::PumiMeshDG{T}) where {T}
  
  @testset "Testing coordinate field transformation" begin

    # test round tripping the vert_coords array
    coord_arr = copy(mesh.vert_coords)
    coord_vec = zeros(T, mesh.dim*mesh.coord_numNodes)

    coords3DTo1D(mesh, coord_arr, coord_vec, PdePumiInterface.AssignReduction{T}())
    fill!(coord_arr, 0.0)
    coords1DTo3D(mesh, coord_vec, coord_arr, PdePumiInterface.AssignReduction{T}())

    @test vecnorm(coord_arr - mesh.vert_coords) < 1e-13

  end

  return nothing
end


"""
  Returns an array of the specified size with random values for the real part
  and zeros for the imaginary part
"""
function rand_realpart(dims...)

  a = rand(Complex128, dims...)
  for i=1:length(a)
    a[i] = real(a[i])
  end

  return a
end


"""
  Tests back propigating the metrics to the 1D vector form
"""
function test_metrics_rev_1d(mesh::PumiMesh{T}, sbp, opts) where {T}

  println("testing metrics_rev_1d")
  h = 1e-20
  pert = Complex128(0, h)

  xvec = zeros(Complex128, mesh.dim*mesh.coord_numNodes)
  xvec_dot = rand_realpart(size(xvec))
  xvec_bar = zeros(Complex128, length(xvec))

  dxidx_bar       = rand_realpart(size(mesh.dxidx))
  jac_bar         = rand_realpart(size(mesh.jac))
  nrm_bndry_bar   = rand_realpart(size(mesh.nrm_bndry))
  nrm_face_bar    = rand_realpart(size(mesh.nrm_face_bar))
  coords_bndry_bar    = rand_realpart(size(mesh.coords_bndry_bar))
  nrm_sharedface_bar = Array{Array{Complex128, 3}}(mesh.npeers)
  for i=1:mesh.npeers
    nrm_sharedface_bar[i] = rand_realpart(size(mesh.nrm_sharedface[i]))
  end


  zeroBarArrays(mesh)
  vert_coords_orig = copy(mesh.vert_coords)

  # get unique-ified xvec_dot
  coords3DTo1D(mesh, mesh.vert_coords, xvec, PdePumiInterface.AssignReduction{T}(), parallel=false)
  xvec .+= pert*xvec_dot
  coords1DTo3D(mesh, xvec, mesh.vert_coords, parallel=true)
  xvec_dot = imag(xvec)/h  # unique-ified version

  diff = maximum(abs.(vert_coords_orig - real(mesh.vert_coords)))
  @assert diff < 1e-15


  # test that coords3DTo1D is the reverse mode of coords1DTo3D
  vert_coords2 = zeros(mesh.vert_coords)
  vert_coords_bar = rand_realpart(size(mesh.vert_coords))
  xvec_bar = zeros(Complex128, length(xvec))

  # forward mode
  coords1DTo3D(mesh, xvec, vert_coords2, parallel=false)
  #xvec .-= pert*xvec_dot
  for i=1:length(xvec)
    xvec[i] = real(xvec[i])
  end
  val1 = sum(imag(vert_coords2)/h .* vert_coords_bar)

  # reverse mode
  coords3DTo1D(mesh, vert_coords_bar, xvec_bar, parallel=true)
  val2 = sum(xvec_bar .* xvec_dot)

  val1 = MPI.Allreduce(val1, MPI.SUM, mesh.comm)
  val2 = MPI.Allreduce(val2, MPI.SUM, mesh.comm)
  @test abs(val1 - val2) < max(abs(val1)*1e-13, 1e-12)

  #----------------------------------------------------------------------------
  # test back-propigation of the metrics (locally)
  zeroBarArrays(mesh)
  fill!(xvec_bar, 0)
  vert_coords_dot = rand_realpart(size(mesh.vert_coords))

  for i=1:length(mesh.vert_coords)
    mesh.vert_coords[i] = real(mesh.vert_coords[i])
  end

  mesh.vert_coords .+= pert*vert_coords_dot
  PdePumiInterface.recalcCoordinatesAndMetrics(mesh, sbp, opts)
  mesh.vert_coords .-= pert*vert_coords_dot

  val1 = sum(imag(mesh.dxidx)/h .* dxidx_bar)               +
         sum(imag(mesh.jac)/h .* jac_bar)                  +
         sum(imag(mesh.nrm_bndry)/h .* nrm_bndry_bar)      +
         sum(imag(mesh.nrm_face)/h .* nrm_face_bar)          +
         sum(imag(mesh.coords_bndry)/h .* coords_bndry_bar)
  for i=1:mesh.npeers
    val1 += sum(imag(mesh.nrm_sharedface[i])/h .* nrm_sharedface_bar[i])
  end

  for i=1:length(mesh.vert_coords)
    mesh.vert_coords[i] = real(mesh.vert_coords[i])
  end


  copy!(mesh.dxidx_bar, dxidx_bar)
  copy!(mesh.jac_bar, jac_bar)
  copy!(mesh.nrm_bndry_bar, nrm_bndry_bar)
  copy!(mesh.nrm_face_bar, nrm_face_bar)
  copy!(mesh.coords_bndry_bar, coords_bndry_bar)
  for i=1:mesh.npeers  # DEBUGGING
    copy!(mesh.nrm_sharedface_bar[i], nrm_sharedface_bar[i])
  end

  getAllCoordinatesAndMetrics_rev(mesh, sbp, opts)
#  getAllCoordinatesAndMetrics_rev(mesh, sbp, opts, xvec_bar, parallel=true)
  val2 = sum(mesh.vert_coords_bar .* vert_coords_dot)
  #val2 = sum(xvec_bar .* xvec_dot)

  @test abs(val1 - val2) < max(abs(val1)*1e-12, 1e-12)

  #----------------------------------------------------------------------------
  # test back-propigation of metrics (parallel)
  zeroBarArrays(mesh)
  fill!(xvec_bar, 0)

  # forward mode
  xvec .+= pert*xvec_dot
  coords1DTo3D(mesh, xvec, mesh.vert_coords, parallel=true)
  PdePumiInterface.recalcCoordinatesAndMetrics(mesh, sbp, opts)
  xvec .-= pert*xvec_dot

  val1 = sum(imag(mesh.dxidx)/h .* dxidx_bar)               +
         sum(imag(mesh.jac)/h .* jac_bar)                  +
         sum(imag(mesh.nrm_bndry)/h .* nrm_bndry_bar)      +
         sum(imag(mesh.nrm_face)/h .* nrm_face_bar)          +
         sum(imag(mesh.coords_bndry)/h .* coords_bndry_bar)
  for i=1:mesh.npeers
    val1 += sum(imag(mesh.nrm_sharedface[i])/h .* nrm_sharedface_bar[i])
  end

  # remove complex parts
  for i=1:length(mesh.vert_coords)
    mesh.vert_coords[i] = real(mesh.vert_coords[i])
  end
  PdePumiInterface.recalcCoordinatesAndMetrics(mesh, sbp, opts)

  # reverse mode
  copy!(mesh.dxidx_bar, dxidx_bar)
  copy!(mesh.jac_bar, jac_bar)
  copy!(mesh.nrm_bndry_bar, nrm_bndry_bar)
  copy!(mesh.nrm_face_bar, nrm_face_bar)
  copy!(mesh.coords_bndry_bar, coords_bndry_bar)
  for i=1:mesh.npeers
    copy!(mesh.nrm_sharedface_bar[i], nrm_sharedface_bar[i])
  end

  getAllCoordinatesAndMetrics_rev(mesh, sbp, opts, xvec_bar)
  val2 = sum(xvec_bar .* xvec_dot)
  val1 = MPI.Allreduce(val1, MPI.SUM, mesh.comm)
  val2 = MPI.Allreduce(val2, MPI.SUM, mesh.comm)

  @test abs(val1 - val2) < max(abs(val1)*1e-12, 1e-12)

  return nothing
end


function test_geoNums(mesh)


  @testset "testing geometric numbering" begin
    geonums = mesh.geoNums

    for dim=0:mesh.dim
      for e in apf.MeshIterator(mesh.m_ptr, dim)
        me = apf.toModel(mesh.m_ptr, e)
        me_dim = apf.getModelType(mesh.m_ptr, me)

        for j=1:mesh.coord_numNodesPerType[dim+1]

          for k=0:(me_dim-1)
            @assert apf.isNumbered(geonums.xiNums, e, j-1, k)
            n_k = apf.getNumberJ(geonums.xiNums, e, j-1, k)
            @test n_k <= geonums.numXiDofs
          end

          for k=me_dim:(mesh.dim-1)
            @assert apf.isNumbered(geonums.xiNums, e, j-1, k)
            n_k = apf.getNumberJ(geonums.xiNums, e, j-1, k)
            @test n_k == geonums.numXiDofs + 1
          end
        end   # end j
      end # end iterator
    end  # end dim

  end # end testset

  return nothing
end


function test_geoMapping(mesh)
# mesh must have loaded a CAD-based geometric model


  @testset "testing geoMapping" begin
    println("testing geometric mapping")
    xivec = zeros(mesh.geoNums.numXiDofs)
    coords_arr = copy(mesh.vert_coords)

    coords_xyzToXi(mesh, coords_arr, xivec)
    fill!(coords_arr, 0)
    coords_XiToXYZ(mesh, xivec, coords_arr)

    println("maxdiff = ", maximum(abs.(coords_arr - mesh.vert_coords)))
    @test maximum(abs.(coords_arr - mesh.vert_coords)) < 1e-8 # CAD is not that accurate

  end

  return nothing
end


function test_geoWrapping(mesh)
# mesh must have loaded a CAD-based geometric model

  println("testing geometric wrapping")
  @testset "testing geoWrapping" begin
    # test the kernel function
    rng = [0.0, 1.0]
    @test abs(PdePumiInterface.wrapPeriodicXi(rng, 1.5) - 0.5) < 1e-12
    @test abs(PdePumiInterface.wrapPeriodicXi(rng, 1.4) - 0.4) < 1e-12
    @test abs(PdePumiInterface.wrapPeriodicXi(rng, 2.4) - 0.4) < 1e-12
    @test abs(PdePumiInterface.wrapPeriodicXi(rng, 3.4) - 0.4) < 1e-12

    rng = [1.5, 3.5]
    @test abs(PdePumiInterface.wrapPeriodicXi(rng, 4.0) - 2.0) < 1e-12

    val = Complex128(4.0, 1.5)
    val2 = Complex128(2.0, 1.5)
    @test abs(PdePumiInterface.wrapPeriodicXi(rng, val) - val2) < 1e-12


    # test the overall function
    
    # get periodic model entities
    g = apf.getModel(mesh.m_ptr)
    periodic_mes = Array{gmi.ModelEntity}(0)
    for me in gmi.Gmi_iter(g, mesh.dim-1)
      if gmi.periodic(g, me, 0)  # only check first dimension
        push!(periodic_mes, me)
      end
    end

    xvec = zeros(Float64, mesh.geoNums.numCoordDofs)
    op = PdePumiInterface.AssignReduction{Float64}()
    coords3DTo1D(mesh, mesh.vert_coords, xvec, op, parallel=false)

    xivec = zeros(mesh.geoNums.numXiDofs)
    coords_xyzToXi(mesh, xvec, xivec)
    xivec2 = copy(xivec)

    rng = zeros(Float64, 2)
    geonums = mesh.geoNums
    for (e, dim) in apf.FieldEntityIt(mesh.m_ptr, mesh.coordshape_ptr)
      me = apf.toModel(mesh.m_ptr, e)
      me_dim = apf.getModelType(mesh.m_ptr, me)
      if me in periodic_mes
        # add a multiple of the period to the xivec2
        gmi.range(g, me, 0, rng)
        delta = rng[2] - rng[1]

        idx = apf.getNumberJ(geonums.xiNums, e, 0, 0)
        xivec2[idx] += 2*delta
      end
    end

    # call wrapXiVals, check xivec2 == xivec
    PdePumiInterface.wrapXiCoords(mesh, xivec2)
    @test maximum(abs.(xivec - xivec2)) < 1e-12


  end  # end testset

  return nothing
end



function test_geoDerivative(mesh)
# mesh must have loaded a CAD-based geometric model.
# test against finite difference

  function f_x(i, xvec::AbstractVector)
    # compute f(x)

      coeff1 = i % 10
      coeff2 = i % 5

      x_i = xvec[i]
#      val = coeff1*x_i*x_i + coeff2*x_i + i
      val = coeff1*x_i 

    return val
  end  # end function

  function f_x_deriv(i, xvec::AbstractVector, df_dx::AbstractVector)
    # compute f(x)

    #for i=1:length(xvec)
#    i = 19
      coeff1 = i % 10
      coeff2 = i % 5

      x_i = xvec[i]
      #df_dx[i] = 2*coeff1*x_i + coeff2
      df_dx[i] = coeff1
    #end
  end  # end function

  function f2_x(xvec::AbstractVector)

    val = 0.0
    for i=1:length(xvec)
      val += xvec[i]*xvec[i]
    end

    return val
  end

  function f2_x_deriv(xvec::AbstractVector, df_dx::AbstractVector)
    for i=1:length(xvec)
      df_dx[i] = 2*xvec[i]
    end
  end



  @testset "Geometric Derivative" begin
    println("\ntesting geometric derivative")
    xvec = getXCoords(mesh)
    xvec_pert = copy(xvec)
    xivec = getXiCoords(mesh)
    #=
    xivec = zeros(mesh.geoNums.numXiDofs)
    coords_xyzToXi(mesh, xvec, xivec)
    coords_XiToXYZ(mesh, xivec, xvec)  # hopefully this reduces numerical
                                       # error when compared to coords_XiToXYZ
    =#
    df_dx = zeros(mesh.geoNums.numCoordDofs)
    df_dxi = zeros(mesh.geoNums.numXiDofs)

    xi_indices = constructGeoMapping(mesh)
    h = 1e-6

#=
    # test each component individually
    # this test is rather slow, so don't run it.
    maxdiff = 0.0
    for i=1:length(xvec)
      # function at initial x value
      val1 = f_x(i, xvec)
      fill!(df_dx, 0)
      f_x_deriv(i, xvec, df_dx)
      coords_dXTodXi(mesh, df_dx, df_dxi)

      # only do j values related to i
      for j in xi_indices[i]
        if length(xi_indices[i]) != 1
          continue
        end
        # compute finite difference
        xivec[j] += h
        coords_XiToXYZ(mesh, xivec, xvec_pert)
        val2 = f_x(i, xvec_pert)
        xivec[j] -= h
        
        dfdxi_fd = (val2 - val1)/h
        #println("df_dxi = ", df_dxi[j])
        #println("df_dxi_fd = ", dfdxi_fd)

        diff = abs(df_dxi[j] - dfdxi_fd)
        #@test abs(df_dxi[j] - dfdxi_fd) < 1e-5

        if diff > maxdiff
          maxdiff = diff
        end
      end  # end j
    end  # end i
    println("maxdiff = ", maxdiff)
=#
    h = 1e-6
    pert = rand(length(xivec))

    # finite difference
    J1 = f2_x(xvec)
    for i=1:length(pert)
      xivec[i] += h*pert[i]
    end
    coords_XiToXYZ(mesh, xivec, xvec_pert)
    J2 = f2_x(xvec_pert)
    for i=1:length(pert)
      xivec[i] -= h*pert[i]
    end
    val1 = (J2 - J1)/h


    # using CAD derivative
    dJdx = zeros(xvec)
    dJdxi = zeros(xivec)
    f2_x_deriv(xvec, dJdx)
    coords_dXTodXi(mesh, dJdx, dJdxi)

    val2 = dot(dJdxi, pert)

    @test abs(val1 - val2) < 1e-2

  end  # end testset

  return nothing
end


function constructGeoMapping(mesh)
# get array of all xi indices that correspond to a given x index


  geonums = mesh.geoNums
  indices = Array{Array{Int, 1}}(geonums.numCoordDofs)
  for i=1:mesh.geoNums.numCoordDofs
    indices[i] = Array{Int}(0)
  end

  for (e, dim) in apf.FieldEntityIt(mesh.m_ptr, mesh.coordshape_ptr)
    for j=1:mesh.coord_numNodesPerType[dim+1]
      for k=1:mesh.dim
        idx_in = apf.getNumberJ(geonums.coordNums, e, j-1, k-1)
        for d=1:mesh.dim
          idx_out = apf.getNumberJ(geonums.xiNums, e, j-1, d-1)
          if idx_out != geonums.numXiDofs + 1
            push!(indices[idx_in], idx_out)
          end  # end if
        end  # end d
      end   # end k
    end  # end j
  end  # end iterator

  return indices
end


function test_setPoint(mesh, opts)

  coords_orig = zeros(Float64, 3)
  coords = zeros(Float64, 3)
  coords2 = zeros(Float64, 3)
  coords3 = zeros(Float64, 3)
  xi_orig = zeros(Float64, 3)
  xi = zeros(Float64, 3)
  xi3 = zeros(Float64, 3)

  geoname = opts["dmg_name"]
  if endswith(geoname, ".x_t") || endswith(geoname, ".smd")
    println("\nTesting setPoint")
    @test mesh.geoNums.can_eval

    println("testing setCoordsXi")
    for (e, dim) in apf.FieldEntityIt(mesh.m_ptr, mesh.coordshape_ptr)
      for j=1:mesh.coord_numNodesPerType[dim+1]
        g = apf.getModel(mesh.m_ptr)
        me = apf.toModel(mesh.m_ptr, e)
        me_dim = apf.getModelType(mesh.m_ptr, me)

        # skip things that don't have parametric coords/can't change them
        if me_dim == 0 || me_dim == 3
          continue
        end

        getCoordsXi(mesh, e, j-1, xi_orig, coords_orig)
        getCoordsXi(mesh, e, j-1, xi, coords)

        for i=1:me_dim
          xi[i] += 0.1
        end
      
        setCoordsXi(mesh, e, j-1, xi, coords2)
        getCoordsXi(mesh, e, j-1, xi3, coords3)

        @test maximum(abs.(coords2 - coords3)) < 1e-8
        @test maximum(abs.(xi - xi3)) < 1e-13
      end
    end

  end  # end if

  # test setCoords
  println("testing setCoords")
  for (e, dim) in apf.FieldEntityIt(mesh.m_ptr, mesh.coordshape_ptr)
    for j=1:mesh.coord_numNodesPerType[dim+1]
      g = apf.getModel(mesh.m_ptr)
      me = apf.toModel(mesh.m_ptr, e)
      me_dim = apf.getModelType(mesh.m_ptr, me)

      getCoords(mesh, e, j-1, coords_orig)

      if mesh.geoNums.can_eval
        if me_dim > 0 && me_dim < 3
          getCoordsXi(mesh, e, j-1, xi_orig, coords)
        else
          getCoords(mesh, e, j-1, coords)
        end

        # test snapping
        for i=1:mesh.dim
          coords[i] += 0.1
        end

        setCoords(mesh, e, j-1, coords, true)

        if me_dim == 0
          getCoords(mesh, e, j-1, coords3)
          @test maximum(abs.(coords3 - coords_orig)) < 1e-8
        elseif me_dim == 3
          getCoords(mesh, e, j-1, coords3)
          @test maximum(abs.(coords3 - coords)) < 1e-8
        else
          getCoordsXi(mesh, e, j-1, xi3, coords3)
          gmi.geval(g, me, xi3, coords2)
          @test maximum(abs.(coords2 - coords)) < 1e-8
          @test maximum(abs.(coords3 - coords)) < 1e-8
        end

        # test non snapping
        fill!(coords2, 0); fill!(coords3, 0)
        fill!(xi, 0); fill!(xi3, 0)

        setCoords(mesh, e, j-1, coords_orig, true)
        for i=1:3
          coords[i] = coords_orig[i] + 0.1
        end
        setCoords(mesh, e, j-1, coords, false)

        if me_dim == 0
          getCoords(mesh, e, j-1, coords3)
          @test maximum(abs.(coords3 - coords)) < 1e-8
        elseif me_dim == 3
          getCoords(mesh, e, j-1, coords3)
          @test maximum(abs.(coords3 - coords)) < 1e-8
        else
          # check that the xi coordinates changed, even though they are
          # not consistent with the xyz coordinates
          getCoordsXi(mesh, e, j-1, xi3, coords3)
          gmi.geval(g, me, xi3, coords2)
          coords_tmp = zeros(Float64, 3)
          gmi.geval(g, me, xi_orig, coords_tmp)
          @test maximum(abs.(coords2 - coords)) > 1e-8
          @test maximum(abs.(coords3 - coords)) < 1e-8
        end
      else  # cannot eval

        for i=1:mesh.dim
          coords[i] += 0.1
        end

        setCoords(mesh, e, j-1, coords, true)
        getCoords(mesh, e, j-1, coords2)

        @test maximum(abs.(coords2 - coords)) < 1e-13
      end
    end
  end

  return nothing
end


"""
  Tests the single-element version of update_coords with snapping
"""
function test_geoWarp(mesh, sbp, opts)

  # test single-element update_coords
  coords_orig = copy(mesh.vert_coords)
  
  # double all coordinates
  coords_i = zeros(Float64, 2, 3)
  for i=1:mesh.numEl
    for j=1:3
      coords_i[1, j] = coords_orig[1, j, i] + 0.0001
      coords_i[2, j] = coords_orig[2, j, i] + 0.0001
    end
    update_coords(mesh, i, coords_i)
  end

  commit_coords(mesh, sbp, opts)

  # check xyz and xi are consistent
  xyz_j = zeros(Float64, 3)
  xi_j = zeros(Float64, 2)
  xyz2_j = zeros(Float64, 3)
  g = apf.getModel(mesh.m_ptr)
  for e in apf.MeshIterator(mesh.m_ptr, 0)
    apf.getPoint(mesh.m_ptr, e, 0, xyz_j)
    apf.getParam(mesh.m_ptr, e, xi_j)
    me = apf.toModel(mesh.m_ptr, e)
    gmi.geval(g, me, xi_j, xyz2_j)

    for k=1:mesh.dim
      @test abs(xyz2_j[k] - xyz_j[k]) < 1e-5  # closestPoint is not that accurate
    end
  end

  return nothing
end


function test_update_coords(mesh, sbp, opts)

  @assert mesh.geoNums.can_eval

  # test updateCoords works
  xivec = getXiCoords(mesh)
  xvec = getXCoords(mesh)
  println("initially, xvec = \n", xvec)

  xivec .+= 1e-3*rand(length(xivec))
  update_coordsXi(mesh, sbp, opts, xivec)
  update_coords(mesh, sbp, opts, xvec)
  xvec2 = getXCoords(mesh)

  println("xvec = \n", xvec)
  println("xvec2 = \n", xvec2)
  println("diff = \n", xvec2 - xvec)
  @test maximum(abs.(xvec2 - xvec)) < 1e-12

  return nothing
end


function test_coords_file(opts, sbp, sbpface, mesh_ctor=PumiMeshDG2)


  mesh = mesh_ctor(Float64, sbp, opts, sbpface)

  # perturb coordinates
  if mesh.geoNums.can_eval
    xivec = getXiCoords(mesh)
    pert = 0.001*rand(size(xivec))
    xivec .+= pert
    update_coordsXi(mesh, sbp, opts, xivec)
  else
    xvec = getXCoords(mesh)
    pert = 0.01*rand(size(xvec))
    xvec .+= pert
    update_coords(mesh, sbp, opts, xvec)
  end
    
  writeCoordsBinary(mesh, "coordsb.dat")

  # load new mesh and check coordinates
  mesh2 = mesh_ctor(Float64, sbp, opts, sbpface)
  readCoordsBinary(mesh2, sbp, opts, "coordsb.dat")

  println("coords maxdiff = ", maximum(abs.(mesh.coords - mesh2.coords)))
  @test maximum(abs.(mesh.coords - mesh2.coords)) < 1e-15

  if mesh.geoNums.can_eval
    xivec2 = getXiCoords(mesh2)
    @test maximum(abs.(xivec2 - xivec)) < 1e-15
  end

  # test loading when mesh does not support xi coordinates
  opts2 = copy(opts)
  opts2["dmg_name"] = ".null"
  mesh3 = mesh_ctor(Float64, sbp, opts, sbpface)
  readCoordsBinary(mesh3, sbp, opts, "coordsb.dat")
  @test maximum(abs.(mesh.coords - mesh3.coords)) < 1e-15

  if mesh.geoNums.can_eval
    # save from mesh that does not support xi coordinates and load on mesh
    # that does
    writeCoordsBinary(mesh3, "coordsb2.dat")

    mesh4 = mesh_ctor(Float64, sbp, opts, sbpface)
    readCoordsBinary(mesh4, sbp, opts, "coordsb2.dat")

    @test maximum(abs.(mesh.coords - mesh4.coords)) < 1e-15
  
    xivec2 = getXiCoords(mesh4)
    @test maximum(abs.(xivec2 - xivec)) < 1e-7
  end

  #TODO: test Complex128
  return nothing
end


function test_refcounting()

  order = 1
  smb_name = "vortex.smb"
  dmg_name = "vortex.dmg"

  opts = Dict{Any, Any}(
    "dmg_name" => dmg_name,
    "smb_name" => smb_name,
    "order" => order,
    "coloring_distance" => 2,
    "numBC" => 1,
    "BC1" =>  [4, 7, 10, 13],
    "run_type" => 4,
    )
#    interp_op = [0.5 0 0; 0 0.5 0; 0 0 0.5]

  sbp = getTriSBPOmega(degree=order)
  vtx = sbp.vtx
  sbpface = TriFace{Float64}(order, sbp.cub, vtx)

  mesh = PumiMeshDG2(Float64, sbp, opts, sbpface, dofpernode=4)

  #fshape_ptr = apf.getConstantShapePtr(0)
  fshape_ptr = mesh.coordshape_ptr
  f_ptr = apf.createPackedField(mesh, "testfield1", 1, fshape_ptr)

  @test PdePumiInterface.countFieldRef(f_ptr) == 1

  f_ptr2 = apf.createPackedField(mesh, "testfield1", 1, fshape_ptr)
  @test f_ptr2 != f_ptr
  @test PdePumiInterface.countFieldRef(f_ptr2) == 1

  # load the same mesh from disk again, check the fields are not shared
  mesh2 = PumiMeshDG2(Float64, sbp, opts, sbpface, dofpernode=4)

  f_ptr3 = apf.createPackedField(mesh2, "testfield1", 1, fshape_ptr)

  @test mesh.m_ptr != mesh2.m_ptr
  for f in mesh.fields.orig
    @test !(f in mesh2.fields.orig)
  end
  for f in mesh.fields.user
    @test !(f in mesh2.fields.user)
  end

  m1_ptr = mesh.m_ptr
  finalize(mesh)
  @test apf.countMeshRefs(m1_ptr) == 0
  @test PdePumiInterface.countFieldRef(f_ptr) == 0

  # make sure finalizer is freeing fields
  finalize(mesh2)
  @test PdePumiInterface.countFieldRef(f_ptr3) == 0


  # create 2 Julia meshes sharing the same apf::Mesh
  mesh1 = PumiMeshDG2(Float64, sbp, opts, sbpface, dofpernode=4)

  #fshape_ptr = apf.getConstantShapePtr(0)
  fshape_ptr = mesh.coordshape_ptr
  f_orig_mesh1 = apf.createPackedField(mesh1.m_ptr, "origfield_mesh1",
                                       1, fshape_ptr)
  PdePumiInterface.attachOrigField(mesh1, f_orig_mesh1)

  n_orig_mesh1 = apf.createNumberingJ(mesh1.m_ptr, "orignum_mesh1", fshape_ptr,
                                      1)
  PdePumiInterface.attachOrigNumbering(mesh1, n_orig_mesh1)


  f_user_mesh1 = apf.createPackedField(mesh1.m_ptr, "userfield_mesh1",
                                       1, fshape_ptr)
  PdePumiInterface.attachUserField(mesh1, f_user_mesh1)

  n_user_mesh1 = apf.createNumberingJ(mesh1.m_ptr, "usernum_mesh1", fshape_ptr,
                                      1)
  PdePumiInterface.attachUserNumbering(mesh1, n_user_mesh1)


  mesh2 = PumiMeshDG2(mesh1, sbp, opts)

  n_user_mesh2 = apf.createNumberingJ(mesh2, "usernum_mesh2", 1, fshape_ptr)
  f_user_mesh2 = apf.createPackedField(mesh2, "userfield_mesh2", 1, fshape_ptr)

  # test orig fields are attached to both, but user fields are not
  @test f_orig_mesh1 in mesh2.fields.orig
  @test !(f_user_mesh1 in mesh2.fields.user)
  @test f_user_mesh1 in mesh1.fields.user

  @test n_orig_mesh1 in mesh2.numberings.orig
  @test !(n_user_mesh1 in mesh2.numberings.user)
  @test n_user_mesh1 in mesh1.numberings.user

  @test PdePumiInterface.countFieldRef(f_orig_mesh1) == 2
  @test PdePumiInterface.countFieldRef(f_user_mesh1) == 1

  @test PdePumiInterface.countNumberingRef(n_orig_mesh1) == 2
  @test PdePumiInterface.countNumberingRef(n_user_mesh1) == 1

  # write values to field
  vals = Float64[1.0]
  for (e, dim) in apf.FieldEntityIt(mesh1.m_ptr, fshape_ptr)
    apf.numberJ(n_orig_mesh1, e, 0, 0, 1)
    apf.numberJ(n_user_mesh1, e, 0, 0, 1)
    apf.numberJ(n_user_mesh2, e, 0, 0, 1)

    apf.setComponents(f_orig_mesh1, e, 0, vals)
    apf.setComponents(f_user_mesh1, e, 0, vals)
    apf.setComponents(f_user_mesh2, e, 0, vals)
  end


  writeVisFiles(mesh1, "mesh1_vis")
  writeVisFiles(mesh2, "mesh2_vis")

  mesh1_has_f_orig1 = false
  mesh1_has_n_orig1 = false
  mesh1_has_f_user1 = false
  mesh1_has_n_user1 = false

  mesh2_has_f_orig1 = false
  mesh2_has_n_orig1 = false
  mesh2_has_f_user1 = false
  mesh2_has_n_user1 = false

  mesh1_has_f_user2 = false
  mesh1_has_n_user2 = false
  mesh2_has_f_user2 = false
  mesh2_has_n_user2 = false

  for line in eachline("mesh1_vis/mesh1_vis.pvtu")
    if contains(line, "origfield_mesh1")
      mesh1_has_f_orig1 = true
    end
    if contains(line, "orignum_mesh1")
      mesh1_has_n_orig1 = true
    end
    if contains(line, "userfield_mesh1")
      mesh1_has_f_user1 = true
    end
    if contains(line, "usernum_mesh1")
      mesh1_has_n_user1 = true
    end
    if contains(line, "userfield_mesh2")
      mesh1_has_f_user2 = true
    end
    if contains(line, "usernum_mesh2")
      mesh1_has_n_user2 = true
    end
  end

  @test mesh1_has_f_orig1
  @test mesh1_has_n_orig1
  @test mesh1_has_f_user1
  @test mesh1_has_n_user1
  @test !mesh1_has_f_user2
  @test !mesh1_has_n_user2

  for line in eachline("mesh2_vis/mesh2_vis.pvtu")
    if contains(line, "origfield_mesh1")
      mesh2_has_f_orig1 = true
    end
    if contains(line, "orignum_mesh1")
      mesh2_has_n_orig1 = true
    end
    if contains(line, "userfield_mesh1")
      mesh2_has_f_user1 = true
    end
    if contains(line, "usernum_mesh1")
      mesh2_has_n_user1 = true
    end
    if contains(line, "userfield_mesh2")
      mesh2_has_f_user2 = true
    end
    if contains(line, "usernum_mesh2")
      mesh2_has_n_user2 = true
    end
  end

  @test mesh2_has_f_orig1
  @test mesh2_has_n_orig1
  @test !mesh2_has_f_user1
  @test !mesh2_has_n_user1
  @test mesh2_has_f_user2
  @test mesh2_has_n_user2

  finalize(mesh1)
  @test PdePumiInterface.countFieldRef(f_orig_mesh1) == 1
  @test PdePumiInterface.countFieldRef(f_user_mesh1) == 0

  @test PdePumiInterface.countNumberingRef(n_orig_mesh1) == 1
  @test PdePumiInterface.countNumberingRef(n_user_mesh1) == 0

  finalize(mesh2)

  return nothing
end
