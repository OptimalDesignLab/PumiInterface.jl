"""
  Compare the fields of two meshes where the first has Tmsh = Float64 and
  the second one has Tmsh = Complex128
"""
function compare_meshes(mesh, mesh_c)

  # decide which fields to check
  fieldnames = [:vert_coords, :coords_bndry, :coords_interface, :dxidx, 
                :jac, :nrm_face, :nrm_bndry]

  if mesh.coord_order == 1
    push!(fieldnames, :dxidx_face, :dxidx_bndry, :jac_face, :jac_bndry)
  end

  fieldnames_shared = [:coords_sharedface, :dxidx_sharedface, :jac_sharedface,
                       :nrm_sharedface]

  if mesh.coord_order == 1
    push!(fieldnames_shared, :dxidx_sharedface, :jac_sharedface)
  end

  # check that the real parts of the field are the same and the imaginary
  # parts are zero
  @testset "----- Comparing meshes -----" begin
    for i in fieldnames

      f_real = getfield(mesh, i)
      f_imag = getfield(mesh_c, i)

      if length(f_real) != 0
        @test isapprox( vecnorm(f_real - real(f_imag))/sqrt(length(f_real)), 0.0) atol=1e-14
        @test isapprox( vecnorm(imag(f_imag)), 0.0) atol=1e-14
      end
    end

    for i in fieldnames_shared

      f_real = getfield(mesh, i)
      f_imag = getfield(mesh_c, i)

      for peer=1:mesh.npeers
        f_real_peer = f_real[peer]
        f_imag_peer = f_imag[peer]
        
        @test isapprox( vecnorm(f_real_peer - real(f_imag_peer))/sqrt(length(f_real_peer)), 0.0) atol=1e-14
        @test isapprox( vecnorm(imag(f_imag_peer)), 0.0) atol=1e-14
      end
    end
  end

  return nothing
end


