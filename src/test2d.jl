using ODLCommonTools
using SummationByParts
using PumiInterface
using PdePumiInterface

dmg_name =  ".null"
smb_name = "tri8q.smb"
order = 1

sbp = TriSBP{Float64}(degree=order)

#mesh = PumiMesh2{Float64}(dmg_name, smb_name, order, sbp, dofpernode=4, shape_type=1)

arr = [3 4 5; 7 6 1]
min, max = PdePumiInterface.getMinandMax(arr)

println("min = ", min, ", max = ", max)
