push!(LOAD_PATH, "../src")
using Base.Test
using PumiInterface
using ODLCommonTools
using SummationByParts
using PdePumiInterface

function not_isapprox(args...; kwargs...)
  return !isapprox(args...; kwargs...)
end

include("defs.jl")
include("common_functions.jl")
include("test_funcs.jl")
include("test_adapt.jl")

@testset "----- Testing 4 process PDEPumiInterface3DG -----" begin
  degree = 1
  Tsbp = Float64
  sbp = getTetSBPOmega(degree=degree)
#  sbp = TetSBP{Tsbp}(degree=degree, internal=true)
  ref_verts = sbp.vtx
  interp_op = SummationByParts.buildinterpolation(sbp, ref_verts.')
  face_verts = SummationByParts.SymCubatures.getfacevertexindices(sbp.cub)
  topo = ElementTopology{3}(face_verts)
  sbpface = TetFace{Tsbp}(degree, sbp.cub, ref_verts)

  dmg_name = ".null"
  smb_name = "pcube10.smb"

  opts = PdePumiInterface.get_defaults()
  opts["dmg_name"] = dmg_name
  opts["smb_name"] = smb_name
  opts["use_linear_metrics"] = true
  opts["coloring_distance"] = 2
  opts["order"] = degree
  opts["numBC"] = 1
  opts["BC1"] = [0,1,2,3,4,5]

#  interp_op = eye(4)
  mesh = PumiMeshDG3(Float64, sbp, opts, sbpface, topo)
#  mesh = PumiMeshDG3{Float64, typeof(sbpface)}(dmg_name, smb_name, degree, sbp, opts, sbpface, topo)

  # check coloring
  @test  mesh.maxColors  < 18
  @test  mesh.numColors  < mesh.maxColors+1
  @test  mesh.numColors  > 6

  for i=1:mesh.numEl
    el_i = mesh.elements[i]
    cnt = countBridgeAdjacent(mesh.m_ptr, el_i, mesh.dim-1, mesh.dim)
    @test  countnz(mesh.pertNeighborEls[i, :])  > cnt
  end

  # check mapping interpolation
  for i=1:mesh.npeers

    dxidx_peer = mesh.dxidx_sharedface[i]
    jac_peer = mesh.jac_sharedface[i]
    for j=1:mesh.peer_face_counts[i]
      iface_i = mesh.shared_interfaces[i][j]
      el_j = iface_i.elementL
      dxidx_el = mesh.dxidx[:, :, 1, el_j]
      jac_el = mesh.jac[:, el_j]
      jac_face = jac_peer[:, j]

      for k=1:mesh.numNodesPerFace
        dxidx_face = dxidx_peer[:, :, k, j]
        for p=1:3
          for n=1:3
            @test isapprox( dxidx_face[n, p], dxidx_el[n, p]) atol=1e-13
          end
        end
        @test isapprox( jac_face[k], jac_el[k]) atol=1e-13
      end
    end   # end loop over peer faces
  end  # end loop over peers


  mesh_c = PumiMeshDG3(Complex128, sbp, opts, sbpface, topo)
#  mesh_c = PumiMeshDG3{Complex128, typeof(sbpface)}(dmg_name, smb_name, degree, sbp, opts, sbpface, topo)
  mesh = PumiMeshDG3(Float64, sbp, opts, sbpface, topo)
#  mesh = PumiMeshDG3{Float64, typeof(sbpface)}(dmg_name, smb_name, degree, sbp, opts, sbpface, topo)
  compare_meshes(mesh, mesh_c)


  # note: this trashes the mesh arrays
  test_parallel_metrics(mesh, sbp, opts)

  # this allocates a ton of memory - don't run
#  test_metric_rev(mesh, mesh_c, sbp)
end

MPI.Barrier(MPI.COMM_WORLD)
if MPI.Initialized()
  MPI.Finalize()
end
