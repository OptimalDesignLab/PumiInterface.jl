using SummationByParts


function nodecalc(sbp::TetSBP, isDG::Bool)
  vtx = [0.0 0 0; 1 0 0; 0 1 0; 0 0 1]
  r1 = vtx[1, :]
  r2 = vtx[2, :]
  r3 = vtx[3, :]
  r4 = vtx[4, :]

  T = zeros(3,3)
  T[:, 1] = r2 - r1
  T[:, 2] = r3 - r1
  T[:, 3] = r4 - r1

  # create operator
  if isDG
    println("getting DG mesh coordinates")
    coords = SummationByParts.SymCubatures.calcnodes(sbp.cub, vtx)
  else
    coords = calcnodes(sbp, vtx)
  end

  xi = zeros(coords)

  for i=1:size(coords,2)
    xi[:, i] = T\(coords[:, i] - r1.')
  end

  return xi, coords
end

#=
(tmp, nnodes) = size(xi)


for i=1:nnodes
  println("node ", i, " coords: ", coords[:,i])
end

print("\n")

for i=1:nnodes
  xi4 = 1 - sum(xi[:,i])
  println("node ", i, " coords: ", xi[:,i], " , xi4 = ", xi4)
end
=#
