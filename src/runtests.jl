# testing using FactCheck
# because PUMI is not available on Travis, all tests are run locally


using FactCheck
using PumiInterface
using PDESolverCommon
using SummationByParts
using PdePumiInterface



# create mesh
dmg_name = ".null"
#smb_name = "../../mesh_files/quarter_vortex3l.smb"
# smb_name = "../../mesh_files/quarter_vortex8l.smb"
smb_name = "./tri8l.smb"


facts("Testing PUMIInterface.jl") do

  # do tests here
  context("Testing initilization") do

  order = 1  # linear elements
  num_Entities, m_ptr, mshape_ptr = init2(dmg_name, smb_name, order)
  @fact num_Entities[1] => 9
  @fact num_Entities[2] => 16
  @fact num_Entities[3] => 8
  @fact num_Entities[4] => 0
  @fact m_ptr => truthy
  @fact mshape_ptr => truthy 

  # check the countJ() function
  @fact countJ(m_ptr, 0) => 9
  @fact countJ(m_ptr, 1) => 16
  @fact countJ(m_ptr, 2) => 8
  @fact countJ(m_ptr, 3) => 0

  # check numberings
  vertN_ptr = getVertNumbering()
  edgeN_ptr = getEdgeNumbering()
  faceN_ptr = getFaceNumbering()


#  printEdgeVertNumbers(edgeN_ptr, vertN_ptr)

  resetVertIt()
  resetEdgeIt()
  resetFaceIt()

  vert = getVert()
  edge = getEdge()
  face = getFace()

  @fact getNumberJ(vertN_ptr, vert, 0, 0) => 0
  @fact getNumberJ(edgeN_ptr, edge, 0, 0) => 0
  @fact getNumberJ(faceN_ptr, face, 0, 0) => 0


  for i=1:num_Entities[1]  # loop over verts
    entity = getVert()
    @fact getNumberJ(vertN_ptr, entity, 0, 0) => i-1
    @fact getVertNumber() => i-1
    @fact getVertNumber2(entity) => i-1
    incrementVertIt()
  end
  resetVertIt()

  for i=1:num_Entities[2]  # loop over edges
    entity = getEdge()
    @fact getNumberJ(edgeN_ptr, entity, 0, 0) => i-1
    @fact getEdgeNumber() => i-1
    @fact getEdgeNumber2(entity) => i-1
    incrementEdgeIt()
  end
  resetEdgeIt()


  for i=1:num_Entities[3]  # loop over faces
    entity = getFace()
    @fact getNumberJ(faceN_ptr, entity, 0, 0) => i-1
    @fact getFaceNumber() => i-1
    @fact getFaceNumber2(entity) => i-1
    incrementFaceIt()
  end
  resetFaceIt()




  @fact getMeshDimension(m_ptr) => 2
  @fact getType(m_ptr, vert) => 0
  @fact getType(m_ptr, edge) => 1
  @fact getType(m_ptr, face) => 2

  # need to recreate these test post reordering
#=
  # check downward adjacency function
  down_entities, num_down = getDownward(m_ptr, face, 1) # face -> edges
  down_nums = zeros(Int, num_down)
  for i=1:num_down
    down_nums[i] = getNumberJ(edgeN_ptr, down_entities[i], 0, 0)
  end
#  println("down_nums = ", down_nums)
   @fact down_nums => [0, 1, 2]

  down_entities, num_down = getDownward(m_ptr, face, 0) # face -> verts
  down_nums = zeros(Int, num_down)
  for i=1:num_down
    down_nums[i] = getNumberJ(vertN_ptr, down_entities[i], 0, 0)
  end
   @fact down_nums => [0, 1, 4]

  down_entities, num_down = getDownward(m_ptr, edge, 0) # edge -> verts
  down_nums = zeros(Int, num_down)
  for i=1:num_down
    down_nums[i] = getNumberJ(vertN_ptr, down_entities[i], 0, 0)
  end
   @fact down_nums => [0, 1]

=#
  # check upward adjacency function
  num_up = countAdjacent(m_ptr, vert, 1)  # vert -> edge
  up_entities = getAdjacent( num_up)
  up_nums = zeros(Int, num_up)
  for i=1:num_up
    up_nums[i] = getNumberJ(edgeN_ptr, up_entities[i], 0, 0)
  end

#  println("up_nums = ", up_nums)
  @fact up_nums => [4, 2, 0]

 num_up = countAdjacent(m_ptr, vert, 2)  # vert -> face
  up_entities = getAdjacent( num_up)
  up_nums = zeros(Int, num_up)
  for i=1:num_up
    up_nums[i] = getNumberJ(faceN_ptr, up_entities[i], 0, 0)
  end

#  println("up_nums = ", up_nums)
  @fact up_nums => [1, 0]


  num_up = countAdjacent(m_ptr, edge, 2)  # edge -> face
  up_entities = getAdjacent( num_up)
  up_nums = zeros(Int, num_up)
  for i=1:num_up
    up_nums[i] = getNumberJ(faceN_ptr, up_entities[i], 0, 0)
  end

#  println("up_nums = ", up_nums)
  @fact up_nums => [0]



  # check FieldShape and EntityShape functions
  @fact hasNodesIn(mshape_ptr, 0) => true
  @fact hasNodesIn(mshape_ptr, 1) => false
  @fact hasNodesIn(mshape_ptr, 2) => false
  
  @fact countNodesOn(mshape_ptr, 0) => 1
  @fact countNodesOn(mshape_ptr, 1) => 0
  @fact countNodesOn(mshape_ptr, 2) => 0

  # check shape functions, jacobians, integration points here?

  # check coordinates
#=
  coords_vert = zeros(3,1)
  getVertCoords(vert, coords_vert, 3, 1)
  coords = coords_vert[:,1]
  @fact coords => roughly([-1.0, -1.0, 0])
  getVertCoords(coords_vert, 3, 1)
  coords = coords_vert[:,1]
  @fact coords => roughly([-1.0, -1.0, 0])


  coords_edge = zeros(3,2)
  getEdgeCoords(edge, coords_edge, 3, 2)
  @fact coords_edge => roughly([-1.0 0; -1 -1; 0 0])
  getEdgeCoords(coords_edge, 3, 2)
  @fact coords_edge => roughly([-1.0 0; -1 -1; 0 0])



  coords_el = zeros(3,3)
  getFaceCoords(face, coords_el, 3, 3)
  @fact coords_el => roughly([-1.0 0 0; -1 -1 0; 0 0 0])
  getFaceCoords(coords_el, 3, 3)
 @fact coords_el => roughly([-1.0 0 0; -1 -1 0; 0 0 0])

=#

 # check creating numberings
 n_ptr = createNumberingJ(m_ptr, "testnumbering", mshape_ptr, 1)

 numberJ(n_ptr, vert, 0, 0, 1)
 @fact getNumberJ(n_ptr, vert, 0, 0) => 1


  end

end


facts("Testing PdePumiInterface3.jl") do

  # do tests here
  context("Testing initilization") do
  dmg_name =  ".null"
  smb_name = "cube8l.smb"
  order = 1
  sbp = TetSBP{Float64}(degree=order)

  mesh = PumiMesh3{Float64}(dmg_name, smb_name, order, sbp, dofpernode=4)

  @fact PdePumiInterface.constructNodemap(mesh, 0) => [1, 3, 2]
  @fact PdePumiInterface.constructNodemap(mesh, 1) => [3, 2, 1]
  @fact PdePumiInterface.constructNodemap(mesh, 2) => [2, 1, 3]

  order = 2
  mesh = PumiMesh3{Float64}(dmg_name, smb_name, order, sbp, dofpernode=4)

  @fact PdePumiInterface.constructNodemap(mesh, 0) => [1, 3, 2, 6, 5, 4]
  @fact PdePumiInterface.constructNodemap(mesh, 1) => [3, 2, 1, 5, 4, 6]
  @fact PdePumiInterface.constructNodemap(mesh, 2) => [2, 1, 3, 4, 6, 5]

  order = 3
  mesh = PumiMesh3{Float64}(dmg_name, smb_name, order, sbp, dofpernode=4)

  @fact PdePumiInterface.constructNodemap(mesh, 0) => [1, 3, 2, 9, 8, 7, 6, 5, 4, 10]
  @fact PdePumiInterface.constructNodemap(mesh, 1) => [3, 2, 1, 7, 6, 5, 4, 9, 8, 10]
  @fact PdePumiInterface.constructNodemap(mesh, 2) => [2, 1, 3, 5, 4, 9, 8, 7, 6, 10]

  order = 4
  mesh = PumiMesh3{Float64}(dmg_name, smb_name, order, sbp, dofpernode=4)

  @fact PdePumiInterface.constructNodemap(mesh, 0) => [1, 3, 2, 12, 11, 10, 9, 8, 7, 6, 5, 4, 13, 15, 14]
  @fact PdePumiInterface.constructNodemap(mesh, 1) => [3, 2, 1, 9, 8, 7, 6, 5, 4, 12, 11, 10, 15, 14, 13]
  @fact PdePumiInterface.constructNodemap(mesh, 2) => [2, 1, 3, 6, 5, 4, 12, 11, 10, 9, 8, 7, 14, 13, 15]


  end

end


