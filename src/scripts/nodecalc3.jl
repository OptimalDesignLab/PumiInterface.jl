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

function printCaseStatement(xi)
  nnodes = size(xi, 2)
  
  f = open("caseout.txt", "w")
  for i=1:nnodes
    xi1 = xi[1, i]
    xi2 = xi[2, i]
    if size(xi, 1) == 2
      xi3 = 0.0
    else
      xi3 = xi[3, i]
    end
    println(f, "case $(i-1):")
    println(f, "  {xi = Vector3($xi1, $xi2, $xi3); break; }")
  end

  println(f, "default:")
  println(f, "  {xi = Vector3(0, 0, 0); break; }")
  close(f)
end

function printBaryCoords(xi, coords)

  (tmp, nnodes) = size(xi)

  for i=1:nnodes
    println("node ", i, " coords: ", coords[:,i])
  end

  print("\n")

  for i=1:nnodes
    xi4 = 1 - sum(xi[:,i])
    println("node ", i, " coords: ", xi[:,i], " , xi4 = ", xi4)
  end
end


#=
sbp = TetSBP{Float64}(degree=4, reorder=false, internal=false)
xi, coords = nodecalc(sbp, true)
printBaryCoords(xi, coords)
#printCaseStatement(xi)
=#
