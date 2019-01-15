# functions related to mesh warping/coordinate updating

"""
  This function updates the coordinates of the vertices of a given element.
  In parallel, it is the responsibility of the caller to update the coordinates
  consistenly for shared elements.  Note that the user must not change 
  the topology of the mesh (ie. the adjacency of mesh entities), only 
  their positions.

  **Inputs**

   * elnum: the element number to update
   * coords_new: the `mesh.dim` x `coord_numNodesPerElement` array of new
            coordinates.
   * snap: if true, snap the coordinates to the geometry (if the geometric
           model supports snapping).  Setting this to false is strongly
           discouraged and can cause inaccurate coordinate and derivative
           evaluations in other function.  If true, coords_new will be
           overwritten with the snapped coordinates


  For `coords_new`, each column contains the coordinate of a coordinate node,
  starting
  with the vertices, followed by the edges (if quadratic).
  The entities must be ordered the same as the downward adjacencies of
  the elements.  The mesh.topo field combined with mesh.element_vertnums
  can be used to figure this out.

  Note that even if `coords_new` is complex valued, this function only uses
  the real part.  If you need to propigate complex number through the
  metric calculation, see [`recalcCoordinatesAndMetrics`](@ref).

  After completeing all calls to update_coords, users *must* call
  commit_coords().

  This function is supported for backwards compatability, see the other
  methods of `update_coords` which work on the entire mesh at once, which may
  be more efficient.
"""
function update_coords(mesh::PumiMesh, elnum::Integer,  coords_new::AbstractMatrix, snap::Bool=true)

#  @assert size(coords_new, 1) == mesh.dim
#  @assert (size(coords_new, 2) == mesh.dim + 1 || size(coords_new, 2) == mesh.dim+1 + mesh.numTypePerElement[2]) # number of verts or number of verts + edges

  verts = Array{Ptr{Void}}(12)
  coords_j = zeros(Float64, 3)
  newx_j = zeros(Float64, 3)
  xi_j = zeros(Float64, 2)
  g = apf.getModel(mesh.m_ptr)
  can_eval = gmi.can_eval(g)
  snap = snap && can_eval

  el_i = mesh.elements[elnum]
  apf.getDownward(mesh.m_ptr, el_i, 0, verts)

  if mesh.coord_order == 1
    @assert (size(coords_new, 2) >= mesh.numTypePerElement[1])
  elseif mesh.coord_order == 2
    @assert (size(coords_new, 2) == (mesh.numTypePerElement[1] + mesh.numTypePerElement[2]))
  else
    throw(ErrorException("Unsupported coordinate order"))
  end

  #TODO: make this a proper master loop construct
  for j=1:(mesh.dim + 1)  # loop over verts
    for k=1:mesh.dim
      coords_j[k] = real(coords_new[k, j])
    end

    # get parametric coordinates if possible
    if can_eval
      has_xi = getSnappedCoords(mesh, verts[j], snap, coords_j, newx_j, xi_j)

      # set the values
      if has_xi
        apf.setParam(mesh.m_ptr, verts[j], xi_j)
      end
      apf.setPoint(mesh.m_ptr, verts[j], 0, newx_j)
    else
      apf.setPoint(mesh.m_ptr, verts[j], 0, coords_j)
    end  # if can_eval

  end  # end j

  if apf.hasNodesIn(mesh.coordshape_ptr, 1)
    offset = mesh.dim + 1
    nedges = apf.getDownward(mesh.m_ptr, el_i, 1, verts)
    for j=1:nedges
      for k=1:mesh.dim
        coords_j[k] = real(coords_new[k, j + offset])
      end

      if can_eval
        has_xi = getSnappedCoords(mesh, verts[j], snap, coords, newx, xi)

        # set the values
        #TODO: I think this is broken in Pumi for edges, re-enable when fixed
        #if has_xi
        #  apf.setParam(mesh.m_ptr, verts[j], xi_j)
        #end
        apf.setPoint(mesh.m_ptr, verts[j], 0, newx_j)
      else
        apf.setPoint(mesh.m_ptr, verts[j], 0, coords_j)
      end  # if can eval
    end  # end j
  end


  return nothing
end


"""
  Internal function to snap coordinates to the geometry.  Only call this
  function if the geometric model supports evaluating coordinates.

  **Inputs**

   * mesh
   * entity: the MeshEntity
   * snap: if the xyz coordinates should be snapped
   * coords: vector of length 3 containing xyz coordinates

  **Inputs/Outputs**
  
   * newx: new xyz coordinates, snapped to the geometry if `snap` is true,
           same as `coords` if false
   * xi: the CAD parametric coordinates

  **Outputs**

   * has_xi: Bool, if true, the model entity that `entity` is classified on
             has parametric coordinates.  If false, the value of `xi` on exit
             is undefined.
"""
function getSnappedCoords(mesh::PumiMesh, entity::Ptr{Void}, snap::Bool,
                          coords::AbstractVector, newx::AbstractVector,
                          xi::AbstractVector)

  g = apf.getModel(mesh.m_ptr)
  me = apf.toModel(mesh.m_ptr, entity)
  me_dim = apf.getModelType(mesh.m_ptr, me)
  has_xi = me_dim != 3  # only snap if the parametric coordinates exist

  if has_xi
    if me_dim == 0  # Simmetrix doesn't have closest_point for vertices
      fill!(xi, 0)
      gmi.geval(g, me, xi, newx)
    else
      gmi.closest_point(g, me, coords, newx, xi)
    end
  end

  if !snap || !has_xi
    for k=1:mesh.dim
      newx[k] = coords[k]
    end
  end

  return has_xi
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
