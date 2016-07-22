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
  interp_op = SummationByParts.buildinterpolation(sbp, ref_verts.')
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
            @fact dxidx_face[n, p] --> roughly(dxidx_el[n, p], atol=1e-13)
          end
        end
        @fact jac_face[k] --> roughly(jac_el[k], atol=1e-13)
      end
    end   # end loop over peer faces
  end  # end loop over peers



end
