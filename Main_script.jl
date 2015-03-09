# julia script to test PUMI interface



#open all the needed libraries
#dlopen("libapf_mpas", RTLD_GLOBAL)
#dlopen("libapf_sim", RTLD_GLOBAL)
#dlopen("libapf", RTLD_GLOBAL)
#dlopen("libapf_zoltan", RTLD_GLOBAL)
#dlopen("libgmi_sim", RTLD_GLOBAL)
#dlopen("libgmi", RTLD_GLOBAL)
#dlopen("libma", RTLD_GLOBAL)
#dlopen("libmds", RTLD_GLOBAL)
#dlopen("libparma", RTLD_GLOBAL)
#dlopen("libpcu", RTLD_GLOBAL)
#dlopen("libph", RTLD_GLOBAL)
#dlopen("libspr", RTLD_GLOBAL)

include("funcs1.jl")
include("funcs2.jl")
declareNames();  # declare global variable names
# initilize mesh
downward_counts_tmp, num_entities_tmp = init()

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

coords = Array(Float64, 3, 1)   # pass an array 3 by n (3 coordinates each for n points)
(m,n) = size(coords)
getVertCoords(coords, m, n)

coords = zeros(3,2)
(m,n) = size(coords)
getEdgeCoords(coords, m, n)

# get face coordinates
numP = downward_counts[3, 1]
coords = zeros(3,numP)
(m,n) = size(coords)
getFaceCoords(coords, m, n)
println("from julia, face coords = ", coords)
print("\n")


# get element coordinates
numP = downward_counts[4,1]
coords = zeros(3, numP)
(m,n) = size(coords)
getElCoords(coords, m,n)
println("from julia, element coords = ", coords)

resetElIt()

print("\n")

# get all vertex coordinates
vertCoords = getAllVertCoords()
println("verts coords = ", vertCoords)

sbpMatrix = makeSBPMatrix()
println("sbpMatrix = ", sbpMatrix)


resetVertIt()
i = 3
setGlobalVertNumber(i)
j = getGlobalVertNumber()
println("j = ", j)


i += 1
setGlobalVertNumber(i)
j = getGlobalVertNumber()
println("j = ", j)


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
