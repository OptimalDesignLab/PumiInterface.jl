

using PDESolverCommon
using SummationByParts
using PumiInterface
using PdePumiInterface

dmg_name =  ".null"
smb_name = "cube8f.smb"
order = 4
shape_type = 1

sbp = TetSBP{Float64}(degree=order)

mesh = PumiMesh3{Float64}(dmg_name, smb_name, order, sbp, dofpernode=4, shape_type=shape_type)

#=
for i=1:mesh.numBoundaryFaces
  println("Element ", mesh.boundary_nums[i,1] -1, " face ", mesh.boundary_nums[i,2])
end

el = mesh.elements[1]
downwards, numdown = getDownward(mesh.m_ptr, el, 2)
face = downwards[4]


#=
for i=1:numdown
  numup = countAdjacent(mesh.m_ptr, downwards[i], 3)
  println("i = ", i, " numup = ", numup)
end
=#

which, flip, rotate  = getAlignment(mesh.m_ptr, el, face) 
=#
