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

     #= 
      for block_el = 1:size(Eone_orig, 3)
        println("block_el = ", block_el, ", i = ", i)
        for d=1:mesh.dim
          println("d = ", d)
          val = sum(Eone[:, d, block_el]) 
          println("val = ", val)
          @assert val < 1e-14
        end
      end
    =#
      

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
    println("forward mode verts")
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

    println("forward mode nrm_face")
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

    println("forward mode nrm_bndry")
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


