using PDESolverCommon
using SummationByParts
using PumiInterface
using PdePumiInterface

dmg_name =  ".null"
smb_name = "tri8q.smb"
order = 4

sbp = TriSBP{Float64}(degree=order)

mesh = PumiMesh2{Float64}(dmg_name, smb_name, order, sbp, dofpernode=4, shape_type=1)


