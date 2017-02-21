# functions related to mesh warping/coordinate updating

"""
  This function updates the coordinates of the vertices of a given element.
  In parallel, it is the responsibility of the caller to update the coordinates
  consistenly for shared elements.  Note that the user must not change 
  the topology of the mesh (ie. the adjacency of mesh entities), only 
  their positions.

  elnum is the element number to update
  coords_new is the mesh.dim x numVertsPerElement array of new coordinates.
    Each column contains the coordinates for a vertex

  After completeing all calles the update_coords, users *must* call
  commit_coords()
"""
function update_coords(mesh::PumiMesh, elnum::Integer,  coords_new::AbstractArray{Float64, 2})

#  @assert size(coords_new, 1) == mesh.dim
#  @assert (size(coords_new, 2) == mesh.dim + 1 || size(coords_new, 2) == mesh.dim+1 + mesh.numTypePerElement[2]) # number of verts or number of verts + edges

  verts = Array(Ptr{Void}, 12)
  coords_j = zeros(Float64, 3)
  
  el_i = mesh.elements[elnum]
  getDownward(mesh.m_ptr, el_i, 0, verts)

  if mesh.coord_order == 1
    @assert (size(coords_new, 2) >= mesh.numTypePerElement[1])
  elseif mesh.coord_order == 2
    @assert (size(coords_new, 2) == (mesh.numTypePerElement[1] + mesh.numTypePerElement[2]))
  else
    throw(ErrorException("Unsupported corodinate order"))
  end

  #TODO: make this a proper master loop construct
  for j=1:(mesh.dim + 1)  # loop over verts
    for k=1:mesh.dim
      coords_j[k] = coords_new[k, j]
    end
    setPoint(mesh.m_ptr, verts[j], 0, coords_j)
  end

  if hasNodesIn(mesh.coordshape_ptr, 1)
    offset = mesh.dim + 1
    nedges = getDownward(mesh.m_ptr, el_i, 1, verts)
    for j=1:nedges
      for k=1:mesh.dim
        coords_j[k] = coords_new[k, j + offset]
      end
      setPoint(mesh.m_ptr, verts[j], 0, coords_j)
    end
  end


  return nothing
end

"""
  This must be called after all calls to update_coords are complete. It
  update all fields of the mesh that are derived from the coordinates
  (dxidx, jac, node coordinate etc.), given than update_coords does not
  change the mesh topology.

  The verify keyword argument determines whether or not to run the Pumi
  verifier on the updated mesh.  It can be expensive, so it is recommended
  to only run it on small meshes.
"""
function commit_coords(mesh::PumiMesh, sbp; verify=true)
# users must call this function when they have finished updating the coordinates
# the verify kwarg determines if the Pumi verifier on the new mesh
# that should check for negative volumes

  acceptChanges(mesh.m_ptr)
  if verify
    Verify(mesh.m_ptr)
  end

  getAllCoordinatesAndMetrics(mesh, sbp)

end
