push!(LOAD_PATH, "../src")
using FactCheck
using PumiInterface
using ODLCommonTools
using SummationByParts
using PdePumiInterface
include("defs.jl")


facts("----- Testing 4 process PDEPumiInterface3DG -----") do
  degree = 1
  Tsbp = Float64
  sbp = TetSBP{Tsbp}(degree=degree, reorder=false, internal=true)
  ref_verts = sbp.vtx
  interp_op = SummationByParts.buildinterpolation(sbp, ref_verts)
  face_verts = SummationByParts.SymCubatures.getfacevertexindices(sbp.cub)
  topo = ElementTopology{3}(face_verts)
  sbpface = TetFace{Tsbp}(degree, sbp.cub, ref_verts)

  dmg_name = ".null"
  smb_name = "pcube10.smb"

  opts = PdePumiInterface.get_defaults()
  opts["numBC"] = 1
  opts["BC1"] = [0,1,2,3,4,5]

  interp_op = eye(4)
  mesh = PumiMeshDG3{Float64}(dmg_name, smb_name, degree, sbp, opts, interp_op, sbpface, topo)

  # check coloring
  @fact mesh.maxColors --> less_than(18)
  @fact mesh.numColors --> less_than(mesh.maxColors+1)
  @fact mesh.numColors --> greater_than(6)

  for i=1:mesh.numEl
    el_i = mesh.elements[i]
    cnt = countBridgeAdjacent(mesh.m_ptr, el_i, mesh.dim-1, mesh.dim)
    @fact countnz(mesh.pertNeighborEls[i, :]) --> greater_than(cnt)
  end



end

