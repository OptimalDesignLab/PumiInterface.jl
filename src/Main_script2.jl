# julia script to test PUMI interface



using PumiInterface
include("funcs2.jl")
apf.declareNames();  # declare global variable names
# apf.initilize mesh
dmg_name = "cube.dmg"
smb_name = "tet-mesh-1.smb"
#dmg_name = "apf.reorder_a.dmg"
#smb_name = "apf.reorder_a.smb"
#dmg_name = ".null"
#smb_name = ".smb"
downward_counts_tmp, num_entities_tmp, m2_ptr, mshape_ptr = apf.init(dmg_name, smb_name)
m_ptr = apf.getMeshPtr()
mshape_ptr = apf.getMeshShapePtr()

vertN_ptr = getVertNumbering()
edgeN_ptr = getEdgeNumbering()
faceN_ptr = getFaceNumbering()
elN_ptr = getElNumbering()

global const downward_counts = downward_counts_tmp
global const num_entities = num_entities_tmp

println("downward_counts = ", downward_counts)
println("num_entities = ", num_entities)

checkVars();

resetVertIt();
resetEdgeIt();
resetFaceIt();
resetElIt();


checkNums()

resetVertIt();
resetEdgeIt();
resetFaceIt();
resetElIt();

coords = Array{Float64}(3, 1)   # pass an array 3 by n (3 coordinates each for n points)
(m,n) = size(coords)
apf.getVertCoords(coords, m, n)

coords = zeros(3,2)
(m,n) = size(coords)
apf.getEdgeCoords(coords, m, n)

# get face coordinates
numP = downward_counts[3, 1]
coords = zeros(3,numP)
(m,n) = size(coords)
apf.getFaceCoords(coords, m, n)
println("from julia, face coords = ", coords)
print("\n")


# get element coordinates
numP = downward_counts[4,1]
coords = zeros(3, numP)
(m,n) = size(coords)
apf.getElCoords(coords, m,n)
println("from julia, element coords = ", coords)

resetElIt()

print("\n")

# get all vertex coordinates
vertCoords = getAllVertCoords()
println("verts coords = ", vertCoords)

resetVertIt()
resetEdgeIt()
resetFaceIt()
resetElIt()

i = getVertNumber()
println(" first vertex number = " , i)
incrementVertItn(3);
i = getVertNumber();
println(" 4th vertex number = ", i )

i = getEdgeNumber()
println(" first edge number = " , i)
incrementEdgeItn(3);
i = getEdgeNumber();
println(" 4th Edge number = ", i )

i = getFaceNumber()
println(" first face number = " , i)
incrementFaceItn(3);
i = getFaceNumber();
println(" 4th Face number = ", i )

i = getElNumber()
println(" first element number = " , i)
incrementElItn(3);
i = getElNumber();
println(" 4th element number = ", i )

entity = getVert()
println("typeof(entity) in main = ", typeof(entity) )
i = getVertNumber2(entity)
println("vertex number = ", i)

i = apf.getMeshDimension(m_ptr);
println("mesh dimension from julia = ", i )

i = apf.countNodesOn(mshape_ptr, 2);
println("number of nodes on type 2 = ", i)

numbering_ptr = apf.createNumberingJ(m_ptr, "numberingJ1", mshape_ptr, 1)
println("created numbering")

apf.numberJ(numbering_ptr, entity, 0, 0, 5)
i = apf.getNumberJ(numbering_ptr, entity, 0, 0)
println("number of entity = ", i)

i = apf.getNumberJ(vertN_ptr, entity, 0, 0)
println("global number of vertex entity = ", i)

entity = getEdge()
i = apf.getNumberJ(edgeN_ptr, entity, 0, 0)
println("global number of edge entity = ", i)

entity = getFace()
i = apf.getNumberJ(faceN_ptr, entity, 0, 0)
println("global number of face entity = ", i)

entity = getEl()
i = apf.getNumberJ(elN_ptr, entity, 0, 0)
println("global number of element entity = ", i)

i = apf.getMeshDimension(m2_ptr)
println("mesh2 dimension = ", i)

downwards, npoints = apf.getDownward(m_ptr, entity, 0)
println("npoints = ", npoints)
j = 1
for entity2 in downwards
  num = apf.getNumberJ(vertN_ptr, entity2, 0, 0);
  println("entity $j on current element has number", num)
  j += 1
end

entity_type = apf.getType(m_ptr, entity)
println("type of current entity = ", entity_type)
i = apf.countNodesOn(mshape_ptr, entity_type)
println(" number of nodes on current entity = ", i)

apf.printNumberingName(faceN_ptr)

resetVertIt()
incrementIt(3)
entity = getEl()
i = apf.getNumberJ(elN_ptr, entity, 0, 0)
println("global number of vertex entity = ", i)


tag_ptr = apf.createDoubleTag( m_ptr, "tag1", 2)
apf.setDoubleTag( m_ptr, entity, tag_ptr, [3.5, 2.5])
data_ret = zeros(1,2)
data_ret[1] = 2.0
apf.getDoubleTag( m_ptr, entity, tag_ptr, data_ret)
println("data_ret = ", data_ret)


entity = getVert()
n = apf.countAdjacent(m_ptr, entity, 3)
#println("vertex has ", i, " 3d regions ")

adj = apf.getAdjacent(n)

(adj, n) = apf.getAdjacentFull(m_ptr, entity, 3)
println("vertex has ", n, " 3d regions ")
println("size(adj) = ", size(adj))
for k in adj
  n = apf.getNumberJ(elN_ptr, k, 0, 0)
  m = apf.countDownwards(m_ptr, k)
  println("element number = ", n)
  println("downward adjacency counts = ", m)
end

apf.hasNodesIn(mshape_ptr, 0)
apf.hasNodesIn(mshape_ptr, 1)

i = apf.countNodesOn(mshape_ptr, apf.getType(m_ptr, entity))
println("number of nodes on vertex = ", i)

for i = 0:4
  j = apf.countNodesOn(mshape_ptr, i)
  println("entity of type ", i, "has ", j, " nodes")
end


eshape_ptr = apf.getEntityShape(mshape_ptr, 3)
i = apf.countNodes(eshape_ptr)
println("number of nodes on entity of type 3 = ", i)

i = apf.countAllNodes(mshape_ptr, 3)
println("number of nodes on entity of type 3 = ", i)

apf.printEdgeVertNumbers(edgeN_ptr, vertN_ptr)

apf.writeVtkFiles("vtkTest", m_ptr)

entity = getFace()
etype = apf.getType(m_ptr, entity)
eshape_ptr = apf.getEntityShape(mshape_ptr, etype)


mel_ptr = apf.createMeshElement(m_ptr, entity)
i = apf.countIntPoints(mel_ptr, 5)
println("face requires ", i, " points for 5th order accurate integration")
numN = apf.countNodes(eshape_ptr)

for k=1:i
  coords = apf.getIntPoint(mel_ptr, 5, k)
  weight = apf.getIntWeight(mel_ptr, 5,k)
  sh_vals = apf.getValues2(eshape_ptr, coords)
  shd_vals = apf.getLocalGradients2(eshape_ptr, coords)
  jac = apf.getJacobian(mel_ptr, coords)
  println("point ", k, " has coordinates ", coords, " and weight ", weight)
  println("shape function values = ", sh_vals)
  println("shape function derivative values = ", shd_vals)
  println("jacobian = \n", jac)
  print("\n")
end

resetVertIt()
entity = getVert()
coords = zeros(3,1)
apf.getVertCoords(entity, coords, 3, 1)
println("coords = ", coords)

resetEdgeIt()
entity = getEdge()
coords = zeros(3,2)
apf.getEdgeCoords(entity, coords, 3, 2)
println("coords = ", coords)

resetFaceIt()
entity = getFace()
coords = zeros(3,3)
apf.getFaceCoords(entity, coords, 3, 4)

