push!(LOAD_PATH, "../src")
using FactCheck
using PumiInterface
using ODLCommonTools
using SummationByParts
using PdePumiInterface



  dmg_name = "parallel2.dmg"
  smb_name = "./parallel2.smb"
  order = 1
facts("Testing PUMIInterface.jl") do
  num_Entities, m_ptr, mshape_ptr = init2(dmg_name, smb_name, order)

  myrank = MPI.Comm_rank(MPI.COMM_WORLD)
  @fact num_Entities[1] --> 6
  @fact num_Entities[2] --> 9
  @fact num_Entities[3] --> 4

  @fact countPeers(m_ptr, apfVERTEX) --> 1
  peers = Array(Int32, 1)
  getPeers(m_ptr, peers)


  @fact peers[1] --> 1-myrank
  resetVertIt()

  resetAllIts2()
  verts = Array(Ptr{Void}, num_Entities[1])
  edges = Array(Ptr{Void}, num_Entities[2])
  faces = Array(Ptr{Void}, num_Entities[3])

  for i=1:length(verts)
    verts[i] = getVert()
    incrementVertIt()
  end

  for i=1:length(edges)
    edges[i] = getEdge()
    incrementEdgeIt()
  end

  for i=1:length(faces)
    faces[i] = getFace()
    incrementFaceIt()
  end

  resetAllIts2()

  # check the countRemotes function
  # the only thing that can really be tested is the number of total remotes
  remote_cnt = 0 # count the total number of remotes
  remote_sum = 0 # sum of the remote values
  part_nums = Array(Cint, 1)
  entities = Array(Ptr{Void}, 1)

  sleep(5*myrank)
  for i=1:length(verts)
    nremotes = countRemotes(m_ptr, verts[i])
    remote_cnt += Int(nremotes)
    if nremotes != 0
      getRemotes(part_nums, entities)
      remote_sum += Int(part_nums[1])
    end
  end

  @fact remote_cnt --> 3
  @fact remote_sum --> 3*(1 - myrank)

  remote_cnt = 0
  remote_sum = 0
  for i=1:length(edges)
    nremotes =  countRemotes(m_ptr, edges[i])
    remote_cnt += Int(nremotes)
    if nremotes != 0
      getRemotes(part_nums, entities)
      remote_sum += Int(part_nums[1])
    end
  end

  @fact remote_cnt --> 2
  @fact remote_sum --> 2*(1-myrank)
end


facts("----- Testing PdePumiInterfaceDG -----") do
  opts = Dict{Any, Any}(
    "numBC" => 1,
    "BC1" =>  [0,1,2,3],
    "run_type" => 1,
    "verify_coloring" => true,
    "use_edge_res" => false,
    "write_edge_vertnums" => true,
    "write_face_vertnums" => true,
    "write_boundarynums" => true,
    "write_dxidx" => false,
    "write_coords" => false,
    "write_sparsity" => false,
    "write_offsets" => false,
    "write_counts" => false,
    "write_sparsity_nodebnds" => false,
    "write_dofs" => false
    )

  sbp = TriSBP{Float64}(degree=order, reorder=false, internal=true)

  vtx = [-1. 1 -1; -1 -1 1]
  interp_op = SummationByParts.buildinterpolation(sbp, vtx)
  sbpface = TriFace{Float64}(order, sbp.cub, vtx.')

  mesh = PumiMeshDG2{Float64}(dmg_name, smb_name, order, sbp, opts, interp_op, sbpface)

  @fact mesh.numVert --> 6
  @fact mesh.numEdge --> 9
  @fact mesh.numEl --> 4
  @fact mesh.numInterfaces --> 3
  @fact mesh.myrank --> MPI.Comm_rank(mesh.comm)
  @fact mesh.commsize --> 2  # these tests only work for 2 ranks
  @fact mesh.peer_parts[1] --> 1-(mesh.myrank)
  @fact mesh.peer_face_counts[1] --> 2
  @fact length(mesh.bndries_local[1]) --> 2
  @fact length(mesh.bndries_remote[1]) --> 2
  @fact length(mesh.shared_interfaces[1]) --> 2

  # check the interfaces
  for i=1:length(mesh.shared_interfaces[1])
    interface_i = mesh.shared_interfaces[1][i]
    @fact interface_i.elementL --> greater_than(0)
    @fact interface_i.elementL --> less_than(mesh.numEl + 1)
    @fact interface_i.elementR --> greater_than(mesh.numEl)
    @fact interface_i.elementR --> less_than(mesh.numEl + 3)
  end

  @fact mesh.shared_element_offsets[1] --> (mesh.numEl + 1)
  @fact mesh.shared_element_offsets[2] --> (mesh.numEl + 3)

  # check the interpolation
  # use the fact that dxidx, jac are constant across the element
  for p=1:mesh.npeers
    interfaces_p = mesh.shared_interfaces[p]
    dxidx_p = mesh.dxidx_sharedface[p]
    jac_p = mesh.jac_sharedface[p]
    for i=1:mesh.peer_face_counts[p]
      interface_i = interfaces_p[i]
      dxidx_ref = mesh.dxidx[:, :, 1, interface_i.elementL]
      jac_ref= mesh.jac[1, interface_i.elementL]

      for j=1:mesh.sbpface.numnodes
        for k=1:2
          for m=1:2
            @fact dxidx_p[k, m, j, i] --> roughly(dxidx_p[k,m], atol=1e-12)
          end
        end
        @fact jac_p[j, i] --> roughly(jac_ref, atol=1e-12)
      end
    end
  end

end
MPI.Barrier( MPI.COMM_WORLD)
