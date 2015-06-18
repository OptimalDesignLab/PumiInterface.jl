

using PDESolverCommon
using SummationByParts
using PdePumiInterface

dmg_name =  ".null"
smb_name = "cube8l.smb"
order = 1
sbp = TetSBP{Float64}(degree=order)

mesh = PumiMesh3{Float64}(dmg_name, smb_name, order, sbp)

for i=1:mesh.numBoundaryFaces
  println("Element ", mesh.boundary_nums[i,1] -1, " face ", mesh.boundary_nums[i,2])
end

