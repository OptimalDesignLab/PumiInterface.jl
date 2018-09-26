# test mesh adaptation

function test_adapt_2d(;parallel=false)

  @testset "Testing 2D mesh adaptation" begin

    if parallel
      smb_name = "parallel2.smb"
    else
      smb_name = "tri8l.smb"
    end
    order = 1
    opts = Dict{Any, Any}(
      "dmg_name" => ".null",
      "smb_name" => smb_name,
      "order" => order,
      "coloring_distance" => 2,
      "numBC" => 1,
      "BC1" =>  [0],
      "run_type" => 4)

    PdePumiInterface.set_defaults(opts)

    sbp = getTriSBPOmega(degree=order)
    vtx = sbp.vtx
    sbpface = TriFace{Float64}(order, sbp.cub, vtx)

    mesh = PumiMeshDG2(Float64, sbp, opts, sbpface, dofpernode=1)


    # write linear field to mesh
    u = zeros(mesh.numDof)
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        x = mesh.coords[1, j, i]; y = mesh.coords[2, j, i]
        u[mesh.dofs[1, j, i]] = x + y + 1
      end
    end

    el_sizes = getElementSizes(mesh)
    scale!(el_sizes, 0.5)

    numEl_initial = mesh.numEl

    newmesh, unew = adaptMesh(mesh, sbp, opts, el_sizes, u)

    @test newmesh.numEl == 4*numEl_initial

    for i=1:newmesh.numEl
      for j=1:newmesh.numNodesPerElement
        x = newmesh.coords[1, j, i]; y = newmesh.coords[2, j, i]
        @test isapprox(unew[newmesh.dofs[1, j, i]], x + y + 1) atol=1e-13
      end
    end

  end  # end testset

  return nothing
end


function test_adapt_3d(; parallel=false)

  @testset "Testing 3D mesh adaptation" begin
    degree = 1
    dmg_name = ".null"
    if parallel
      smb_name = "pcube10.smb"
    else
      smb_name = "unitcube.smb"
    end

    opts = PdePumiInterface.get_defaults()
    opts["smb_name"] = smb_name
    opts["dmg_name"] = dmg_name
    opts["order"] = degree
    opts["coloring_distance"] = 2
    opts["use_linear_metrics"] = true
    opts["numBC"] = 1
    opts["BC1"] = [0,1,2,3,4,5]


    sbp = getTetSBPOmega(degree=degree, Tsbp=Float64)
    ref_verts = sbp.vtx
    face_verts = SummationByParts.SymCubatures.getfacevertexindices(sbp.cub)
    topo = ElementTopology{3}(face_verts)
    sbpface = TetFace{Float64}(degree, sbp.cub, ref_verts)

    mesh = PumiMeshDG3(Float64, sbp, opts, sbpface, topo)

    # write linear field to mesh
    println("writing field")
    u = zeros(mesh.numDof)
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        x = mesh.coords[1, j, i]; y = mesh.coords[2, j, i]; z = mesh.coords[3, j, i]
        u[mesh.dofs[1, j, i]] = x + 2*y + 3*z + 1
      end
    end

    el_sizes = getElementSizes(mesh)
    scale!(el_sizes, 0.5)

    numEl_initial = mesh.numEl

    println("adapting mesh")
    newmesh, unew = adaptMesh(mesh, sbp, opts, el_sizes, u)

    @test newmesh.numEl == 8*numEl_initial

    # check that the linear field was interpolated exactly
    for i=1:newmesh.numEl
      for j=1:newmesh.numNodesPerElement
        x = newmesh.coords[1, j, i]; y = newmesh.coords[2, j, i]; z = newmesh.coords[3, j, i]
        @test isapprox(unew[newmesh.dofs[1, j, i]], x + 2*y + 3*z + 1) atol=1e-13
      end
    end

  end  # end testset

  return nothing
end
