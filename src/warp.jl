# functions related to mesh warping/coordinate updating

function update_coords(mesh::PumiMesh, elnum::Integer,  coords_new::AbstractArray{Float64, 2})

  @assert size(coords_new, 1) == mesh.dim
  @assert size(coords_new, 2) == mesh.dim + 1 # number of verts

  verts = Array(Ptr{Void}, 12)
  coords_j = zeros(Float64, 3)
  
  el_i = mesh.elements[elnum]
  getDownward(mesh.m_ptr, el_i, 0, verts)

  for j=1:(mesh.dim + 1)  # loop over verts
    for k=1:mesh.dim
      coords_j[k] = coords_new[k, j]
    end
    setPoint(mesh.m_ptr, verts[j], 0, coords_j)
  end

  return nothing
end


function commit_coords(mesh::PumiMesh)
# users must call this function when they have finished updating the coordinates

  acceptChanges(mesh.m_ptr)
  Verify(mesh.m_ptr)  # TODO: only do this for small meshes

  # TODO: update mesh.coords, dxidx, jac etc.

end
