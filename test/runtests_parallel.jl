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

  dmg_name = "parallel2.dmg"
  smb_name = "./parallel2.smb"
  order = 1
@testset "Testing PUMIInterface.jl" begin

  m_ptr, dim = loadMesh(dmg_name, smb_name, order)
  mshape_ptr, num_Entities, n_arr = initMesh(m_ptr)
#  num_Entities, m_ptr, mshape_ptr = init2(dmg_name, smb_name, order)

  myrank = MPI.Comm_rank(MPI.COMM_WORLD)
  @test ( num_Entities[1] )== 6
  @test ( num_Entities[2] )== 9
  @test ( num_Entities[3] )== 4

  @test ( countPeers(m_ptr, apfVERTEX) )== 1
  peers = Array{Int32}(1)
  getPeers(m_ptr, peers)


  @test ( peers[1] )== 1-myrank

  verts = Array{Ptr{Void}}(num_Entities[1])
  edges = Array{Ptr{Void}}(num_Entities[2])
  faces = Array{Ptr{Void}}(num_Entities[3])

  it = MeshIterator(m_ptr, 0)
  for i=1:length(verts)
    verts[i] = iterate(m_ptr, it)
  end
  free(m_ptr, it)

  it = MeshIterator(m_ptr, 1)
  for i=1:length(edges)
    edges[i] = iterate(m_ptr, it)
  end
  free(m_ptr, it)

  it = MeshIterator(m_ptr, 2)
  for i=1:length(faces)
    faces[i] = iterate(m_ptr, it)
  end
  free(m_ptr, it)

  # check the countRemotes function
  # the only thing that can really be tested is the number of total remotes
  remote_cnt = 0 # count the total number of remotes
  remote_sum = 0 # sum of the remote values
  part_nums = Array{Cint}(1)
  entities = Array{Ptr{Void}}(1)

  sleep(5*myrank)
  for i=1:length(verts)
    nremotes = countRemotes(m_ptr, verts[i])
    remote_cnt += Int(nremotes)
    if nremotes != 0
      getRemotes(part_nums, entities)
      remote_sum += Int(part_nums[1])
    end
  end

  @test ( remote_cnt )== 3
  @test ( remote_sum )== 3*(1 - myrank)

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

  @test ( remote_cnt )== 2
  @test ( remote_sum )== 2*(1-myrank)
end


@testset "----- Testing PdePumiInterfaceDG -----" begin
  opts = Dict{Any, Any}(
    "dmg_name" => dmg_name,
    "smb_name" => smb_name,
    "coloring_distance" => 2,
    "order" => order,
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
    "use_linear_metrics" => true,
    )

  sbp = getTriSBPOmega(degree=order)
#  sbp = TriSBP{Float64}(degree=order, internal=true)

  vtx = [-1. 1 -1; -1 -1 1]
#  interp_op = SummationByParts.buildinterpolation(sbp, vtx)
  sbpface = TriFace{Float64}(order, sbp.cub, vtx.')
  mesh_c = PumiMeshDG2(Complex128, sbp, opts, sbpface)
#  mesh_c = PumiMeshDG2{Complex128, typeof(sbpface)}(dmg_name, smb_name, order, sbp, opts, sbpface)
  mesh = PumiMeshDG2(Float64, sbp, opts, sbpface)
#  mesh = PumiMeshDG2{Float64, typeof(sbpface)}(dmg_name, smb_name, order, sbp, opts, sbpface)

  @test ( mesh.numVert )== 6
  @test ( mesh.numEdge )== 9
  @test ( mesh.numEl )== 4
  @test ( mesh.numInterfaces )== 3
  @test ( mesh.myrank )== MPI.Comm_rank(mesh.comm)
  @test ( mesh.commsize )== 2  # these tests only work for 2 ranks
  @test ( mesh.peer_parts[1] )== 1-(mesh.myrank)
  @test ( mesh.peer_face_counts[1] )== 2
  @test ( length(mesh.bndries_local[1]) )== 2
  @test ( length(mesh.bndries_remote[1]) )== 2
  @test ( length(mesh.shared_interfaces[1]) )== 2

  # check the interfaces
  for i=1:length(mesh.shared_interfaces[1])
    interface_i = mesh.shared_interfaces[1][i]
    @test  interface_i.elementL  > 0
    @test  interface_i.elementL  < mesh.numEl + 1
    @test  interface_i.elementR  > mesh.numEl
    @test  interface_i.elementR  < mesh.numEl + 3
    @test ( interface_i.orient )== 1
  end

  @test ( mesh.shared_element_offsets[1] )== (mesh.numEl + 1)
  @test ( mesh.shared_element_offsets[2] )== (mesh.numEl + 3)

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
          partnums = Array{Cint}(nremotes)
          entities = Array{Ptr{Void}}(nremotes)
          PumiInterface.getRemotes(partnums, entities)
          for p=1:nremotes
            if partnums[p] == peernum
              cnt += 1
            end  # end if
          end  # end loop p
        end  # end if isShared
      end  # end loop k

      @test  cnt  > 0
    end  # end loop j

    # check element numbers are unique
    local_elnums = copy(local_els)
    sort!(local_elnums)
    @test ( local_elnums )== unique(local_elnums)
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
            @test isapprox( dxidx_p[k, m, j, i], dxidx_p[k,m]) atol=1e-12
          end
        end
        @test isapprox( jac_p[j, i], jac_ref) atol=1e-12
      end
    end
  end

  # check the adjacency dictonary
  adj_dict, revadj = PdePumiInterface.getLocalAdjacency(mesh)
  @test ( length(adj_dict) )== mesh.numSharedEl
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
    @test  val   < 3
    @test ( val[1] )!=0
  end

  # check the revadj
  count_dict2 = Dict{Int, Int}()
  for i=1:mesh.numSharedEl
    val1 = revadj[i, 1]
    val2 = revadj[i, 2]
    @test ( val1 )!=0
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
    @test  val  < 3
    @test ( val )!=0
  end




  # check neighbor_colors for uniqueness
  myrank = mesh.myrank
  for i=1:mesh.numEl
    colors_i = mesh.neighbor_colors[:, i]
    sort!(colors_i)
    nzs = countnz(colors_i)
    nz_start = length(colors_i) - nzs + 1
    colors_i2 = colors_i[nz_start:end]
    @test ( colors_i2 )== unique(colors_i2)
  end

  for i=1:mesh.npeers
    numel_i = mesh.shared_element_offsets[i+1] - mesh.shared_element_offsets[i]
    masks = mesh.shared_element_colormasks[i]
    for j = 1:numel_i
      val = 0
      for k = 1:mesh.numColors
        val += masks[k][j]
      end
      @test ( val )== 1  # each element should have exactly 1 color
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
          @test  (mesh.dofs[k, j, i] + mesh.dof_offset)  > dof_min
          @test  (mesh.dofs[k, j, i] + mesh.dof_offset)  < dof_max
        else
          @test  mesh.dofs[k, j, i]  < 1
        end
      end
    end
  end


  vshare = mesh.vert_sharing
  @test ( vshare.npeers )== 1
  @test ( vshare.peer_nums )== [1 - mesh.myrank]
  @test ( vshare.counts[1] )== 3
  @test ( length(vshare.vert_nums[1]) )== 3
  @test ( length(keys(vshare.rev_mapping)) )== 3

  # check vert_nums and rev_mapping in more detail
  for i=1:3
    @test  vshare.vert_nums[1][i]  > 0
    @test  vshare.vert_nums[1][i]  < mesh.numVert + 1
  end

  for (key, val) in vshare.rev_mapping
    @test ( length(val.first) )== 1
    @test ( length(val.second) )== 1
    @test ( val.first[1] )== 1-mesh.myrank  # part number
    @test  val.second[1]  > 0
    @test  val.second[1]  < 4# number of shared vertices + 1
  end

  compare_meshes(mesh, mesh_c)
  test_metric_rev(mesh, mesh_c, sbp, opts)

  MPI.Barrier(mesh.comm)

end


@testset "----- Testing PdePumiInterface3DG -----" begin

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
  smb_name = "pcube2.smb"

  opts = PdePumiInterface.get_defaults()
  opts["dmg_name"] = dmg_name
  opts["smb_name"] = smb_name
  opts["use_linear_metrics"] = true
  opts["coloring_distance"] = 2
  opts["order"] = degree
  opts["numBC"] = 1
  opts["BC1"] = [0,1,2,3,4,5]

#  interp_op = eye(4)
  mesh_c = PumiMeshDG3(Complex128, sbp, opts, sbpface, topo)
#  mesh_c = PumiMeshDG3{Complex128, typeof(sbpface)}(dmg_name, smb_name, degree, sbp, opts, sbpface, topo)
  mesh = PumiMeshDG3(Float64, sbp, opts, sbpface, topo)
#  mesh = PumiMeshDG3{Float64, typeof(sbpface)}(dmg_name, smb_name, degree, sbp, opts, sbpface, topo)

  @test ( mesh.numEl )== 12
  @test ( mesh.numGlobalEl )== 16
#  @test ( mesh.numBoundaryFaces )==#  @fact mesh.numBoundaryFaces --> 4
  @test ( mesh.numInterfaces )== 14
  @test ( mesh.numBoundaryFaces )== 16
  @test ( mesh.myrank )== MPI.Comm_rank(mesh.comm)
  @test ( mesh.commsize )== MPI.Comm_size(mesh.comm)
  
  if mesh.myrank == 0
    @test ( mesh.peer_parts )== [1]
  else
    @test ( mesh.peer_parts )== [0]
  end

  @test ( mesh.peer_face_counts[1] )== 4


  # TODO: test dxidx, jac sharedface
  # check interface array
  for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]
    @test  iface_i.elementL  > 0
    @test  iface_i.elementL  < mesh.numEl+1
    @test  iface_i.elementR  > 0
    @test  iface_i.elementR  < mesh.numEl+1
    @test  iface_i.faceL  > 0
    @test  iface_i.faceL  < 5
    @test  iface_i.faceR  > 0
    @test  iface_i.faceR  < 5
    @test  iface_i.orient  > 0
    @test  iface_i.orient  < 4
  end

  # check boundary array
  for i=1:mesh.numBoundaryFaces
    bndry_i = mesh.bndryfaces[i]
    @test  bndry_i.element  > 0
    @test  bndry_i.element  < mesh.numEl + 1
    @test  bndry_i.face  > 0
    @test  bndry_i.face  < 5
  end

  # check mapping interpolation
  # should be constant within an element for straight-sided elements
  for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]
    el_i = iface_i.elementL
    dxidx_el = mesh.dxidx[:, :, 1, el_i]
    jac_el = mesh.jac[:, el_i]
    jac_face = mesh.jac_face[:, i]

    for j=1:mesh.numNodesPerFace
      dxidx_face = mesh.dxidx_face[:, :, j, i]
      for k=1:3
        for p=1:3
          @test isapprox( dxidx_face[p, k], dxidx_el[p, k]) atol=1e-13
        end
      end
      @test isapprox( jac_face[j], jac_el[j]) atol=1e-13
    end
  end  # end loop over interfaces

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

  # test curvilinear metrics reverse mode
  compare_meshes(mesh, mesh_c)
  test_metric_rev(mesh, mesh_c, sbp, opts)



  # test periodic things
  smb_name = "tet2_pxz_p2_.smb"

  opts = PdePumiInterface.get_defaults()
  opts["smb_name"] = smb_name
  opts["dmg_name"] = dmg_name
  opts["use_linear_metrics"] = true
  opts["coloring_distance"] = 2
  opts["order"] = degree
  opts["numBC"] = 1
  opts["BC1"] = [0,2,4,5]

  mesh = PumiMeshDG3(Float64, sbp, opts, sbpface, topo)
#  mesh = PumiMeshDG3{Float64, typeof(sbpface)}(dmg_name, smb_name, degree, sbp, opts, sbpface, topo)

  @test ( mesh.numPeriodicInterfaces )== 4
  @test ( mesh.numInterfaces )== (mesh.numFace - 8 - 3*4 - 8)
  @test ( mesh.peer_face_counts[1] )== 8

  smb_name = "tet2_pxy_p2_.smb"
  opts["smb_name"] = smb_name
  opts["BC1"] = [1,2,3,4]
  mesh = PumiMeshDG3(Float64, sbp, opts, sbpface, topo)
#  mesh = PumiMeshDG3{Float64, typeof(sbpface)}(dmg_name, smb_name, degree, sbp, opts, sbpface, topo)

  @test ( mesh.numPeriodicInterfaces )== 0
  @test ( mesh.numInterfaces )== (mesh.numFace - 2*8 - 4*4)
  @test ( mesh.peer_face_counts[1] )== 16

  vshare = mesh.vert_sharing
  @test ( vshare.npeers )== 1
  @test ( vshare.peer_nums )== [1 - mesh.myrank]
  @test ( vshare.counts[1] )== 18
  @test ( length(vshare.vert_nums[1]) )== 18
  @test ( length(keys(vshare.rev_mapping)) )== 18

  # check vert_nums and rev_mapping in more detail
  for i=1:18
    @test  vshare.vert_nums[1][i]  > 0
    @test  vshare.vert_nums[1][i]  < mesh.numVert + 1
  end

  for (key, val) in vshare.rev_mapping
    @test ( length(val.first) )== 1
    @test ( length(val.second) )== 1
    @test ( val.first[1] )== 1-mesh.myrank  # part number
    @test  val.second[1]  > 0
    @test  val.second[1]  < 19# number of shared vertices + 1
  end


  MPI.Barrier(mesh.comm)

  # check assertion for ambiguous parallel things
  dmg_name = ".null"
  smb_name = "tet111_.smb"

  opts = PdePumiInterface.get_defaults()
  opts["dmg_name"] = dmg_name
  opts["smb_name"] = smb_name
  opts["order"] = degree
  opts["numBC"] = 0
  opts["coloring_distance"] = 2
 
  @test_throws Exception  PumiMeshDG3(Float64, sbp, opts, sbpface, topo)
#  @test_throws Exception #  @fact_throws PumiMeshDG3{Float64, typeof(sbpface)}(dmg_name, smb_name, degree, sbp, opts, sbpface, topo)

  # check 2 x 2 x 2 mesh works
  dmg_name = ".null"
  smb_name = "tet222_.smb"

  opts["dmg_name"] = dmg_name
  opts["smb_name"] = smb_name

  PumiMeshDG3(Float64, sbp, opts, sbpface, topo)
#  PumiMeshDG3{Float64, typeof(sbpface)}(dmg_name, smb_name, degree, sbp, opts, sbpface, topo)

  # this doesn't work because each partition is a 1 x 1 x 1 cube with periodic faces
  dmg_name = ".null"
  smb_name = "tet211_.smb"
  opts["dmg_name"] = dmg_name
  opts["smb_name"] = smb_name

  @test_throws Exception  PumiMeshDG3(Float64, sbp, opts, sbpface, topo)
#  @test_throws Exception #  @fact_throws PumiMeshDG3{Float64, typeof(sbpface)}(dmg_name, smb_name, degree, sbp, opts, sbpface, topo)



end

MPI.Barrier(MPI.COMM_WORLD)
if MPI.Initialized()
  MPI.Finalize()
end

