# functions related to mesh warping/coordinate updating

"""
  This function updates the coordinates of the vertices of a given element.
  In parallel, it is the responsibility of the caller to update the coordinates
  consistenly for shared elements.  Note that the user must not change 
  the topology of the mesh (ie. the adjacency of mesh entities), only 
  their positions.

  elnum is the element number to update
  coords_new is the mesh.dim x coord_numNodesPerElement array of new coordinates.
    Each column contains the coordinate of a coordinate node, starting
    with the vertices, followed by the edges (if quadratic).
    The entities must be ordered the same as the downward adjacencies of
    the elements.  The mesh.topo field combined with mesh.element_vertnums
    can be used to figure this out.

  Note that even if `coords_new` is complex valued, this function only uses
  the real part.  If you need to propigate complex number through the
  metric calculation, see [`recalcCoordinatesAndMetrics`](@ref).

  After completeing all calles the update_coords, users *must* call
  commit_coords()
"""
function update_coords(mesh::PumiMesh, elnum::Integer,  coords_new::AbstractMatrix)

#  @assert size(coords_new, 1) == mesh.dim
#  @assert (size(coords_new, 2) == mesh.dim + 1 || size(coords_new, 2) == mesh.dim+1 + mesh.numTypePerElement[2]) # number of verts or number of verts + edges

  verts = Array{Ptr{Void}}(12)
  coords_j = zeros(Float64, 3)
  
  el_i = mesh.elements[elnum]
  apf.getDownward(mesh.m_ptr, el_i, 0, verts)

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
      coords_j[k] = real(coords_new[k, j])
    end
    apf.setPoint(mesh.m_ptr, verts[j], 0, coords_j)
  end

  if apf.hasNodesIn(mesh.coordshape_ptr, 1)
    offset = mesh.dim + 1
    nedges = apf.getDownward(mesh.m_ptr, el_i, 1, verts)
    for j=1:nedges
      for k=1:mesh.dim
        coords_j[k] = real(coords_new[k, j + offset])
      end
      apf.setPoint(mesh.m_ptr, verts[j], 0, coords_j)
    end
  end


  return nothing
end

"""
  This must be called after all calls to update_coords are complete. It
  update all fields of the mesh that are derived from the coordinates
  (dxidx, jac, node coordinate etc.), given that update_coords does not
  change the mesh topology.

  The verify keyword argument determines whether or not to run the Pumi
  verifier on the updated mesh.  It can be expensive, so it is recommended
  to only run it on small meshes.

  The `write_vis` keyword controls whether a vtk file, named `pre_commit` is
  written just before the metrics are recalculated (if an element is turned 
  inside out, the metric recalculation is where the error gets thrown).
"""
function commit_coords(mesh::PumiMesh, sbp, opts; verify=true, write_vis=false)
# users must call this function when they have finished updating the coordinates
# the verify kwarg determines if the Pumi verifier on the new mesh
# that should check for negative volumes

  apf.acceptChanges(mesh.m_ptr)
  if verify
    apf.Verify(mesh.m_ptr)
    # TODO; call our verification code too
  end

  if write_vis
    writeVisFiles(mesh, "pre_commit")
  end

  getAllCoordinatesAndMetrics(mesh, sbp, opts)

end
