using SummationByParts

vtx = [0.0 0; 1 0; 0 1]
r1 = vtx[1, :]
r2 = vtx[2, :]
r3 = vtx[3, :]

T = zeros(2,2)
T[:, 1] = r2 - r1
T[:, 2] = r3 - r1

# create operator
sbp1 = TriSBP{Float64}(degree=4)

coords = calcnodes(sbp1, vtx)

xi = zeros(coords)

for i=1:size(coords,2)
  xi[:, i] = T\(coords[:, i] - r1.')
end

xi

