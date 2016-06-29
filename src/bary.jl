# function for converting back and forth between xyz and barycentric coordinates

#=
function calc3dnodes(vtx::Array{Float64, 2})
# the columns of vtx are the x, y, and z coordinates of the nodes

  v1 = vtx[:, 2] - vtx[:, 1]
  v2 = vtx[:, 3] - vtx[:, 1]
  v3 = vtx[:, 4] - vtx[:, 1]
  xi = [1/8 5/8 1/8 1/8;
        1/8 1/8 5/8 1/8;
        1/8 1/8 1/8 5/8]
  T = zeros(3,3)
  T[:, 1] = v1
  T[:, 2] = v2
  T[:, 3] = v3

  coords = zeros(3, 4)
  r1 = vtx[:, 1]
  for i=1:4
    coords[:, i] = T*xi[:, i] + r1
  end

  return coords
end
#=
vtx = [1. 2 1 1;
       1  1 2 1;
       1  1 1 2]
coords = calcnodes(vtx)
println("vtx = \n", vtx)
println("coords = \n", coords)
=#
=#
