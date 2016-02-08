using SummationByParts

function nodecalc(p)
  vtx = [0.0 0; 1 0; 0 1]
  r1 = vtx[1, :]
  r2 = vtx[2, :]
  r3 = vtx[3, :]

  T = zeros(2,2)
  T[:, 1] = r2 - r1
  T[:, 2] = r3 - r1

  # create operator
  sbp1 = TriSBP{Float64}(degree=p)

  coords = calcnodes(sbp1, vtx)

  xi = zeros(coords)

  for i=1:size(coords,2)
    xi[:, i] = T\(coords[:, i] - r1.')
  end

  return xi, coords
end

function minNodeDist(p)
# get the minimum distance between nodes on a reference element of degree p

    xi, coords = nodecalc(p)
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

#minNodeDist(2)
