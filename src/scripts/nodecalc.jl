# thsi script calculates the barycentric coordinates for 2D and 3D
# SBP operator nodes, and prints the case statement needed by 
# apf::FieldShape

using SummationByParts
include("nodecalc3.jl")

function nodecalc(sbp::TriSBP, isDG::Bool)
  vtx = [0.0 0; 1 0; 0 1]
  r1 = vtx[1, :]
  r2 = vtx[2, :]
  r3 = vtx[3, :]

  T = zeros(2,2)
  T[:, 1] = r2 - r1
  T[:, 2] = r3 - r1


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

function nodecalc(sbp::TetSBP, isDG::Bool)

  vtx = [-1.0 -1.0 -1.0
        1.0 -1.0 -1.0
       -1.0 1.0 -1.0
       -1.0 -1.0 1.0]

  r1 = vtx[1, :]
  r2 = vtx[2, :]
  r3 = vtx[3, :]
  r4 = vtx[4, :]
  T = zeros(3,3)
  T[:, 1] = r2 - r1
  T[:, 2] = r3 - r1
  T[:, 3] = r4 - r1


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


function minNodeDist(sbp, isDG::Bool)
# get the minimum distance between nodes on a reference element of degree p

    xi, coords = nodecalc(sbp, isDG)
    min_dist = typemax(Float64)
#    println("coords = ", coords)
    for i=1:size(coords, 2)
      for j=(i+1):size(coords, 2)
	# calculate distance between node i and node j
	dist_j = norm(coords[:, i] - coords[:, j])

	if dist_j < min_dist
	  min_dist = dist_j
	end

      end  # end loop j
    end  # end loop i

#    println("for p=$p elements, min node distance = ", min_dist)
  
  return min_dist
end


#sbp = TriSBP{Float64}(degree=1, reorder=false, internal=false)
sbp = TetSBP{Float64}(degree=4, reorder=false, internal=true)
xi, coords = nodecalc(sbp, true)
printCaseStatement(xi)
printBaryCoords(xi, coords)
#writedlm("coordsout.dat", coords, ' ')
#minNodeDist(2)

