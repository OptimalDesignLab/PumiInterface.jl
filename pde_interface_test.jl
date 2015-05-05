push!(LOAD_PATH, "/users/creanj/julialib_for/PUMI.jl")
#include("pde_interface.jl")
using PdePumiInterface
using PumiInterface
dmg_name = ".null"
smb_name = "tri8l.smb"

mesh = PumiMesh2(dmg_name, smb_name, 1)

el_coords = getElementVertCoords(mesh, [1])
println("el_coords = ", el_coords)

i = getShapeFunctionOrder(mesh)
println("shape function order = ", i)

i = getGlobalNodeNumber(mesh, 1, 2)
println("global node number = ", i)


bnd_array = getBoundaryArray(mesh::PumiMesh2)

for i=1:mesh.numBoundaryEdges
  bndry = bnd_array[i]
  println("element ", bndry.element, " has local edge number ", bndry.face, " on the external boundary")
end

printEdgeVertNumbers(mesh.edge_Nptr, mesh.vert_Nptr)

adjacent_nums, num_adjacent = getAdjacentEntityNums(mesh, 1, 1, 2)
println("adjacent_nums = ", adjacent_nums, " num_adjacent = ", num_adjacent)

edgenum_local = getEdgeLocalNum(mesh, 1, adjacent_nums[1])
println("edgenum_local = ", edgenum_local)

