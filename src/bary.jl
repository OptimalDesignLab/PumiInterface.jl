# function for converting back and forth between xyz and barycentric coordinates

"""
  This function calculates xyz coordinates given the barycentric corodinates
  and the vertices of the triangle xyz coordinates

  Inputs:
    bary: a dim x numpoints array of barycentric coordinates
    vtx: a nvertex x dim array of xyz coordinates of the vertices
         of the triangle (or tetrahedron

  Outputs:
    coords: a dim x numpoints array of xyz coordinates corresponding to 
            the barycentric coordinates
"""
function baryToXY{T2}(bary::AbstractMatrix{T2}, vtx::AbstractMatrix)

  @assert (size(bary, 1) >= 1 && size(bary, 1) <= 3)
  @assert (size(vtx, 2) >= 1 && size(vtx, 2) <= 3)
  @assert (size(vtx, 1) >= 2 && size(vtx, 1) <= 4)

  npoints = size(bary, 2)
  dim = size(bary, 1)
  coords = zeros(bary)

  T = zeros(dim, dim)
  r1 = vtx[1, :]

  for i=1:dim
    T[i, :] = vtx[i+1, :] - r1
  end

  # compute xyz coordinates
  for i=1:npoints
    coords[:, i] = T*bary[:, i] + r1
  end

  return coords
end


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
