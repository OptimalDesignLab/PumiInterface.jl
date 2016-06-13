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
    "write_dofs" => false,
    "write_interfaces" => false,
    "write_boundaries" => false,
    "write_sharedboundaries" => false,
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
    @fact interface_i.orient --> 1
  end

  @fact mesh.shared_element_offsets[1] --> (mesh.numEl + 1)
  @fact mesh.shared_element_offsets[2] --> (mesh.numEl + 3)

  # check the gathering of the elements on the boundaries
  for i=1:mesh.npeers
    peernum = mesh.peer_parts[i]
    local_els = mesh.local_element_lists[i]
    for j = 1:length(local_els)
      el = local_els[j]
      el_ptr = mesh.elements[el]
      down, tmp = PumiInterface.getDownward(mesh.m_ptr, el_ptr, apfEDGE)
      cnt = 0 # count number of shared edges
      for k=1:3
        if isShared(mesh.m_ptr, down[k])
          nremotes = PumiInterface.countRemotes(mesh.m_ptr, down[k])
          partnums = Array(Cint, nremotes)
          entities = Array(Ptr{Void}, nremotes)
          PumiInterface.getRemotes(partnums, entities)
          for p=1:nremotes
            if partnums[p] == peernum
              cnt += 1
            end  # end if
          end  # end loop p
        end  # end if isShared
      end  # end loop k

      @fact cnt --> greater_than(0)
    end  # end loop j

    # check element numbers are unique
    local_elnums = copy(local_els)
    sort!(local_elnums)
    @fact local_elnums --> unique(local_elnums)
  end  # end loop i


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

  # check the adjacency dictonary
  adj_dict, revadj = PdePumiInterface.getLocalAdjacency(mesh)
  @fact length(adj_dict) --> mesh.numSharedEl
  # now check that no element appears more than twice
  count_dict = Dict{Int, Int}()
  for (key, val) in adj_dict
    el1 = val[1]
    el2 = val[2]
    for elnums = 1:2
      val_i = val[elnums]
      if !haskey(count_dict, val_i)
        count_dict[val_i] = 1
      else
        count_dict[val_i] = count_dict[val_i] + 1
      end
    end
  end

  for (key, val) in count_dict
    @fact val  --> less_than(3)
    @fact val[1] --> not(0)
  end

  # check the revadj
  count_dict2 = Dict{Int, Int}()
  for i=1:mesh.numSharedEl
    val1 = revadj[i, 1]
    val2 = revadj[i, 2]
    @fact val1 --> not(0)
    if !haskey(count_dict2, val1)
      count_dict2[val1] = 1
    else
      count_dict2[val1] += 1
    end

    if val2 != 0
      if !haskey(count_dict2, val2)
        count_dict2[val1] = 1
      else
        count_dict2[val1] += 1
      end
    end

  end

  # check no value appears more than twice
  for (key, val) in count_dict2
    @fact val --> less_than(3)
    @fact val --> not(0)
  end




  # check neighbor_colors for uniqueness
  myrank = mesh.myrank
  for i=1:mesh.numEl
    colors_i = mesh.neighbor_colors[:, i]
    sort!(colors_i)
    nzs = countnz(colors_i)
    nz_start = length(colors_i) - nzs + 1
    colors_i2 = colors_i[nz_start:end]
    @fact colors_i2 --> unique(colors_i2)
  end

  for i=1:mesh.npeers
    numel_i = mesh.shared_element_offsets[i+1] - mesh.shared_element_offsets[i]
    masks = mesh.shared_element_colormasks[i]
    for j = 1:numel_i
      val = 0
      for k = 1:mesh.numColors
        val += masks[k][j]
      end
      @fact val --> 1  # each element should have exactly 1 color
    end
  end

  # check dof numbers for non-local elements
  npeer_els = mesh.shared_element_offsets[end] = mesh.shared_element_offsets[1]
  dof_min = mesh.numDofPerNode*mesh.numNodesPerElement*mesh.numEl
  dof_max = mesh.numDofPerNode*mesh.numNodesPerElement*(npeer_els + mesh.numEl)  + 1
  for i=(mesh.numEl + 1):size(mesh.dofs, 3)
    for j=1:mesh.numNodesPerElement
      for k=1:mesh.numDofPerNode
        if mesh.myrank == 0
          @fact (mesh.dofs[k, j, i] + mesh.dof_offset) --> greater_than(dof_min)
          @fact (mesh.dofs[k, j, i] + mesh.dof_offset) --> less_than(dof_max)
        else
          @fact mesh.dofs[k, j, i] --> less_than(1)
        end
      end
    end
  end

end
MPI.Barrier( MPI.COMM_WORLD)
if MPI.Initialized()
  MPI.Finalize()
end
