# Unit tests

push!(LOAD_PATH, "../src")
ENV["PDEPUMIINTERFACE_TESTING"] = true

using Base.Test
using PumiInterface
using ODLCommonTools
import ODLCommonTools.sview
using SummationByParts
using PdePumiInterface

import Base.isapprox

# generic fallback
function isapprox(a, b)
  return a == b
end

function not_isapprox(args...; kwargs...)
  return !isapprox(args...; kwargs...)
end





@testset "Testing PUMIInterface.jl" begin

  # do tests here
  @testset "Testing apf.initilization" begin
  # create mesh
  dmg_name = ".null"
  #smb_name = "../../mesh_files/quarter_vortex3l.smb"
  # smb_name = "../../mesh_files/quarter_vortex8l.smb"
  smb_name = "./tri8l.smb"


  order = 1  # linear elements
  m_ptr, dim = apf.loadMesh(dmg_name, smb_name, order)
  mshape_ptr, num_Entities, n_arr = apf.initMesh(m_ptr)
#  num_Entities, m_ptr, mshape_ptr, dim, n_arr = apf.init2(dmg_name, smb_name, order)
  println("length(n_arr) = ", length(n_arr))

  @test ( num_Entities[1] )== 9
  @test ( num_Entities[2] )== 16
  @test ( num_Entities[3] )== 8
  @test ( m_ptr )!=C_NULL
  @test ( mshape_ptr )!=C_NULL

  # check the apf.countJ() function
  @test ( apf.countJ(m_ptr, 0) )== 9
  @test ( apf.countJ(m_ptr, 1) )== 16
  @test ( apf.countJ(m_ptr, 2) )== 8
  @test ( apf.countJ(m_ptr, 3) )== 0

  # check numberings
  vertN_ptr = n_arr[1] #getVertNumbering()
  edgeN_ptr = n_arr[2] #getEdgeNumbering()
  faceN_ptr = n_arr[3] #getFaceNumbering()

  numberings = [vertN_ptr, edgeN_ptr, faceN_ptr]


#  apf.printEdgeVertNumbers(edgeN_ptr, vertN_ptr)

#  resetVertIt()
#  resetEdgeIt()
#  resetFaceIt()

  vertit = apf.MeshIterator(m_ptr, 0)
  edgeit = apf.MeshIterator(m_ptr, 1)
  faceit = apf.MeshIterator(m_ptr, 2)

  vert = apf.deref(m_ptr, vertit)
  edge = apf.deref(m_ptr, edgeit)
  face = apf.deref(m_ptr, faceit)
#  vert = getVert()
#  edge = getEdge()
#  face = getFace()

  vertnum = apf.getNumberJ(vertN_ptr, vert, 0, 0) 

  edgenum = apf.getNumberJ(edgeN_ptr, edge, 0, 0) 

  facenum = apf.getNumberJ(faceN_ptr, face, 0, 0) 

  @test ( apf.getNumberJ(vertN_ptr, vert, 0, 0) )== 0
  @test ( apf.getNumberJ(edgeN_ptr, edge, 0, 0) )== 0
  @test ( apf.getNumberJ(faceN_ptr, face, 0, 0) )== 0


  for i=1:num_Entities[1]  # loop over verts
    entity = apf.iterate(m_ptr, vertit)
    @test ( apf.getNumberJ(vertN_ptr, entity, 0, 0) )== i-1
#    incrementVertIt()
  end
  apf.free(m_ptr, vertit)
#  resetVertIt()

  for i=1:num_Entities[2]  # loop over edges
    entity = apf.iterate(m_ptr, edgeit)
    @test ( apf.getNumberJ(edgeN_ptr, entity, 0, 0) )== i-1
#    incrementEdgeIt()
  end
  apf.free(m_ptr, edgeit)
#  resetEdgeIt()


  for i=1:num_Entities[3]  # loop over faces
    entity = apf.iterate(m_ptr, faceit)
    @test ( apf.getNumberJ(faceN_ptr, entity, 0, 0) )== i-1
#    incrementFaceIt()
  end
  apf.free(m_ptr, faceit)
#  resetFaceIt()




  @test ( apf.getMeshDimension(m_ptr) )== 2
  @test ( apf.getType(m_ptr, vert) )== 0
  @test ( apf.getType(m_ptr, edge) )== 1
  @test ( apf.getType(m_ptr, face) )== 2


  apf.printEdgeVertNumbers(edgeN_ptr, vertN_ptr)
  apf.printFaceVertNumbers( faceN_ptr, vertN_ptr)
  # need to recreate these test post apf.reordering

  # check downward adjacency function
  down_entities, num_down = apf.getDownward(m_ptr, face, 1) # face -> edges
  down_nums = zeros(Int, num_down)
  for i=1:num_down
    down_nums[i] = apf.getNumberJ(edgeN_ptr, down_entities[i], 0, 0)
  end
#  println("down_nums = ", down_nums)
   @test ( down_nums )== [1, 5, 0]

  down_entities, num_down = apf.getDownward(m_ptr, face, 0) # face -> verts
  down_nums = zeros(Int, num_down)
  for i=1:num_down
    down_nums[i] = apf.getNumberJ(vertN_ptr, down_entities[i], 0, 0)
  end
   @test ( down_nums )== [0, 2, 1]

  down_entities, num_down = apf.getDownward(m_ptr, edge, 0) # edge -> verts
  down_nums = zeros(Int, num_down)
  for i=1:num_down
    down_nums[i] = apf.getNumberJ(vertN_ptr, down_entities[i], 0, 0)
  end
   @test ( down_nums )== [1, 0]


  # check upward adjacency function
  println("testing vert -> edge")
  num_up = apf.countAdjacent(m_ptr, vert, 1)  # vert -> edge
  up_entities = apf.getAdjacent( num_up)
  up_nums = zeros(Int, num_up)
  for i=1:num_up
    up_nums[i] = apf.getNumberJ(edgeN_ptr, up_entities[i], 0, 0)
  end

#  println("up_nums = ", up_nums)
  @test ( up_nums )== [2, 1, 0]

 num_up = apf.countAdjacent(m_ptr, vert, 2)  # vert -> face
  up_entities = apf.getAdjacent( num_up)
  up_nums = zeros(Int, num_up)
  for i=1:num_up
    up_nums[i] = apf.getNumberJ(faceN_ptr, up_entities[i], 0, 0)
  end

#  println("up_nums = ", up_nums)
  @test ( up_nums )== [1, 0]


  num_up = apf.countAdjacent(m_ptr, edge, 2)  # edge -> face
  up_entities = apf.getAdjacent( num_up)
  up_nums = zeros(Int, num_up)
  for i=1:num_up
    up_nums[i] = apf.getNumberJ(faceN_ptr, up_entities[i], 0, 0)
  end

#  println("up_nums = ", up_nums)
  @test ( up_nums )== [0]


  # check second order adjacencies


  # check FieldShape and EntityShape functions
  @test ( apf.hasNodesIn(mshape_ptr, 0) )== true
  @test ( apf.hasNodesIn(mshape_ptr, 1) )== false
  @test ( apf.hasNodesIn(mshape_ptr, 2) )== false
  
  @test ( apf.countNodesOn(mshape_ptr, 0) )== 1
  @test ( apf.countNodesOn(mshape_ptr, 1) )== 0
  @test ( apf.countNodesOn(mshape_ptr, 2) )== 0

  # check shape functions, jacobians, integration points here?

  # check second order adjacencies

  # check vertices connected through edges
  in_dim = 1
  out_dim = 0
  cnt = apf.countBridgeAdjacent(m_ptr, vert, in_dim, out_dim) 
  adj = Array{Ptr{Void}}(cnt)
  adj_num = Array{Int}(cnt)
  apf.getBridgeAdjacent(adj)
  for i=1:cnt
    adj_num[i] = apf.getNumberJ(numberings[out_dim + 1], adj[i], 0, 0)
  end
  sort!(adj_num)
  @test ( adj_num )== [1,2,3]

  # target_dimension = bridge dimension not defined
  in_dim = 1
  out_dim = 2
  cnt = apf.countBridgeAdjacent(m_ptr, vert, in_dim, out_dim) 
  adj = Array{Ptr{Void}}(cnt)
  adj_num = Array{Int}(cnt)
  apf.getBridgeAdjacent(adj)
  for i=1:cnt
    adj_num[i] = apf.getNumberJ(numberings[out_dim + 1], adj[i], 0, 0)
  end
  sort!(adj_num)
  @test ( adj_num )== [0,1]

  # edge connected to vertices by edges
  in_dim = 0
  out_dim = 1
  cnt = apf.countBridgeAdjacent(m_ptr, edge, in_dim, out_dim) 
  adj = Array{Ptr{Void}}(cnt)
  adj_num = Array{Int}(cnt)
  apf.getBridgeAdjacent(adj)
  for i=1:cnt
    adj_num[i] = apf.getNumberJ(numberings[out_dim + 1], adj[i], 0, 0)
  end
  sort!(adj_num)
  @test ( adj_num )== sort!([1,2,3,4,5])

  # edge connected to edges by faces
  in_dim = 2
  out_dim = 1
  cnt = apf.countBridgeAdjacent(m_ptr, edge, in_dim, out_dim) 
  adj = Array{Ptr{Void}}(cnt)
  adj_num = Array{Int}(cnt)
  apf.getBridgeAdjacent(adj)
  for i=1:cnt
    adj_num[i] = apf.getNumberJ(numberings[out_dim + 1], adj[i], 0, 0)
  end
  sort!(adj_num)
  @test ( adj_num )== [1, 5]

  # face connected to edges by vertex
  in_dim = 0
  out_dim = 1
  cnt = apf.countBridgeAdjacent(m_ptr, face, in_dim, out_dim) 
  adj = Array{Ptr{Void}}(cnt)
  adj_num = Array{Int}(cnt)
  apf.getBridgeAdjacent(adj)
  for i=1:cnt
    adj_num[i] = apf.getNumberJ(numberings[out_dim + 1], adj[i], 0, 0)
  end
  sort!(adj_num)
  ans = [0,1,2,3,4,5,6,7,8,9]

  @test ( adj_num )== ans

  # face connected to faces by vert
  in_dim = 0
  out_dim = 2
  cnt = apf.countBridgeAdjacent(m_ptr, face, in_dim, out_dim) 
  adj = Array{Ptr{Void}}(cnt)
  adj_num = Array{Int}(cnt)
  apf.getBridgeAdjacent(adj)
  for i=1:cnt
    adj_num[i] = apf.getNumberJ(numberings[out_dim + 1], adj[i], 0, 0)
  end
  sort!(adj_num)
  @test ( adj_num )== sort!([1, 2,3,4,5,6])

  # face connected to faces by edge
  in_dim = 1
  out_dim = 2
  cnt = apf.countBridgeAdjacent(m_ptr, face, in_dim, out_dim) 
  adj = Array{Ptr{Void}}(cnt)
  adj_num = Array{Int}(cnt)
  apf.getBridgeAdjacent(adj)
  for i=1:cnt
    adj_num[i] = apf.getNumberJ(numberings[out_dim + 1], adj[i], 0, 0)
  end
  sort!(adj_num)
  @test ( adj_num )== [1, 3]


  mshape2 = apf.getFieldShape(1, 4, 2)
  node_entities = apf.getNodeEntities(m_ptr, mshape2, face)
  verts = node_entities[1:3]
  edges = node_entities[4:12]
  faces = node_entities[13:18]
  for i in verts
    @test ( apf.getType(m_ptr, i) )== apf.VERTEX
  end
  for i in edges
    @test ( apf.getType(m_ptr, i) )== apf.EDGE
  end
  for i in faces
    @test ( apf.getType(m_ptr, i) )== apf.TRIANGLE
  end

  @test ( apf.countPeers(m_ptr, apf.VERTEX) )== 0
  @test ( apf.countPeers(m_ptr, apf.EDGE) )== 0

  for i in verts
    @test ( apf.countRemotes(m_ptr, i) )== 0
  end
  for i in edges
    @test ( apf.countRemotes(m_ptr, i) )== 0
  end












  # check coordinates
#=
  coords_vert = zeros(3,1)
  apf.getVertCoords(vert, coords_vert, 3, 1)
  coords = coords_vert[:,1]
  @test isapprox( coords, [-1.0, -1.0, 0]) 
  apf.getVertCoords(coords_vert, 3, 1)
  coords = coords_vert[:,1]
  @test isapprox( coords, [-1.0, -1.0, 0]) 


  coords_edge = zeros(3,2)
  apf.getEdgeCoords(edge, coords_edge, 3, 2)
  @test isapprox( coords_edge, [-1.0 0; -1 -1; 0 0]) 
  apf.getEdgeCoords(coords_edge, 3, 2)
  @test isapprox( coords_edge, [-1.0 0; -1 -1; 0 0]) 



  coords_el = zeros(3,3)
  apf.getFaceCoords(face, coords_el, 3, 3)
  @test isapprox( coords_el, [-1.0 0 0; -1 -1 0; 0 0 0]) 
  apf.getFaceCoords(coords_el, 3, 3)
 @test isapprox( coords_el, [-1.0 0 0; -1 -1 0; 0 0 0]) 

=#

 # check creating numberings
 n_ptr = apf.createNumberingJ(m_ptr, "testnumbering", mshape_ptr, 1)

 apf.numberJ(n_ptr, vert, 0, 0, 1)
 @test ( apf.getNumberJ(n_ptr, vert, 0, 0) )== 1


 # test mesh warping functions

 # find vertex at -1, -1
# resetVertIt()
 it = apf.MeshIterator(m_ptr, 0)
 coords = zeros(3)
 entity = C_NULL
 for i=1:num_Entities[1]
   vert_i = apf.iterate(m_ptr, it)
   apf.getPoint(m_ptr, vert_i, 0, coords)
   if coords[1] < -0.5 && coords[2] < -0.5
     entity = vert_i
     break
   end
 end
  apf.free(m_ptr, it)
 @test ( entity )!=C_NULL

 coords[1] *= 2
 coords[2] *= 2
 apf.setPoint(m_ptr, entity, 0, coords)
 apf.acceptChanges(m_ptr)
 apf.Verify(m_ptr)

 coords2 = copy(coords)
 apf.getPoint(m_ptr, entity, 0, coords2)
 for i=1:2
   @test  abs(coords[i] - coords2[i])  < 1e-12
 end

 apf.writeVtkFiles("warp_test", m_ptr)


 # test periodic mesh
  @test ( apf.hasMatching(m_ptr) )== false

  smb_name = "tri3_px.smb"
  m_ptr, dim = apf.loadMesh(dmg_name, smb_name, order)
  mshape_ptr, num_Entities, n_arr = apf.initMesh(m_ptr)
#  num_Entities, m_ptr, mshape_ptr = apf.init2(dmg_name, smb_name, order)
  shr_ptr = apf.getSharing(m_ptr)
  numVert = num_Entities[1]
  numEdge = num_Entities[2]
  numFace = num_Entities[3]
  # verify all meshentities are owned
#  resetAllIts2(m_ptr)
  nowned = 0

  @test ( apf.hasMatching(m_ptr) )== true

  it = apf.MeshIterator(m_ptr, 0)
  for i =1:numVert
    entity = apf.iterate(m_ptr, it)
    if apf.isOwned(shr_ptr, entity)
      nowned += 1
    end

#    incrementVertIt()
  end
  apf.free(m_ptr, it)
  @test ( nowned )== numVert - 4

  nowned = 0
  it = apf.MeshIterator(m_ptr, 1)
  for i=1:numEdge
    entity = apf.iterate(m_ptr, it)
    if apf.isOwned(shr_ptr, entity)
      nowned += 1
    end
#    incrementEdgeIt()
  end
  apf.free(m_ptr, it)
  @test ( nowned )== numEdge - 3

  ncopies = zeros(Int, 2)
  nmatches = zeros(Int, 2)
  part_nums = Array{Cint}(1)
  copies = Array{Ptr{Void}}(1)
  matches = Array{Ptr{Void}}(1)
#  resetAllIts2(m_ptr)
  it = apf.MeshIterator(m_ptr, 1)
  for i=1:numEdge
    entity = apf.iterate(m_ptr, it)
    n = apf.countCopies(shr_ptr, entity)
    n2 = apf.countMatches(m_ptr, entity)
    ncopies[n+1] += 1  # either 0 or 1 copy
    nmatches[n+1] += 1 # either 0 or 1 match
    if n > 0
      apf.getCopies(part_nums, copies)
      @test ( part_nums[1] )== 0
      e_type = apf.getType(m_ptr, copies[1])
      @test ( e_type )== apf.EDGE
    end
    if n2 > 0
      apf.getMatches(part_nums, matches)
      @test ( part_nums[1] )== 0
      e_type = apf.getType(m_ptr, copies[1])
      @test ( e_type )== apf.EDGE
    end
#    incrementEdgeIt()
  end
  apf.free(m_ptr, it)

  @test ( ncopies[2] )== 6  # 3 x 3 element mesh has 6 shared edges
  @test ( nmatches[2] )== 6

end  # end context
end

#=
@testset "Testing PdePumiInterface3.jl" begin

  # do tests here
  @testset "Testing apf.initilization" begin
  dmg_name =  ".null"
  smb_name = "cube8l.smb"
  order = 1
  sbp = TetSBP{Float64}(degree=order)

  mesh = PumiMesh3{Float64}(dmg_name, smb_name, order, sbp, dofpernode=4)

  @test ( PdePumiInterface.constructNodemap(mesh, 0) )== [1, 3, 2]
  @test ( PdePumiInterface.constructNodemap(mesh, 1) )== [3, 2, 1]
  @test ( PdePumiInterface.constructNodemap(mesh, 2) )== [2, 1, 3]

  order = 2
  mesh = PumiMesh3{Float64}(dmg_name, smb_name, order, sbp, dofpernode=4)

  @test ( PdePumiInterface.constructNodemap(mesh, 0) )== [1, 3, 2, 6, 5, 4]
  @test ( PdePumiInterface.constructNodemap(mesh, 1) )== [3, 2, 1, 5, 4, 6]
  @test ( PdePumiInterface.constructNodemap(mesh, 2) )== [2, 1, 3, 4, 6, 5]

  order = 3
  mesh = PumiMesh3{Float64}(dmg_name, smb_name, order, sbp, dofpernode=4)

  @test ( PdePumiInterface.constructNodemap(mesh, 0) )== [1, 3, 2, 9, 8, 7, 6, 5, 4, 10]
  @test ( PdePumiInterface.constructNodemap(mesh, 1) )== [3, 2, 1, 7, 6, 5, 4, 9, 8, 10]
  @test ( PdePumiInterface.constructNodemap(mesh, 2) )== [2, 1, 3, 5, 4, 9, 8, 7, 6, 10]

  order = 4
  mesh = PumiMesh3{Float64}(dmg_name, smb_name, order, sbp, dofpernode=4)

  @test ( PdePumiInterface.constructNodemap(mesh, 0) )== [1, 3, 2, 12, 11, 10, 9, 8, 7, 6, 5, 4, 13, 15, 14]
  @test ( PdePumiInterface.constructNodemap(mesh, 1) )== [3, 2, 1, 9, 8, 7, 6, 5, 4, 12, 11, 10, 15, 14, 13]
  @test ( PdePumiInterface.constructNodemap(mesh, 2) )== [2, 1, 3, 6, 5, 4, 12, 11, 10, 9, 8, 7, 14, 13, 15]


  end

end
=#


println("about to test pdepumiinterface")
include("test_math.jl")
include("test_funcs.jl")  # include functions used by both 2d and 3D
include("common_functions.jl")
include("test_adapt.jl") # mesh adaptation
include("test_parallel_types.jl")
include("pdepumiinterface.jl")
include("pdepumiinterface3.jl")

MPI.Barrier(MPI.COMM_WORLD)
if MPI.Initialized()
  MPI.Finalize()
end
