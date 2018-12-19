# interface functions for mapping between xyz and CAD parametric coordinates

"""
  This function takes a vector of the mesh xyz coordinate field and turns it
  into a parametric representation using the CAD system.  The transformation
  has the following rules:

   1. MeshEntities that do not have any degrees of freedom in the CAD
      parametric system do not appear in the output (verticies classified on
      geometric vertices)
   2. MeshEntities that are on the interior of the domain retain their xyz
      coordinates (not all CAD system have a parametric representation of 
      3D regions)
   3. All other entities (those that lie on a curve/surface that bounds a
      higher-dimension geometric entity) have their xyz coordinates transformed
      to the CAD parametric coordinates (which are fewer number than the xyz
      coordinates).

  It is assumed that the input coordinates have not been modified in such a way
  that MeshEntities no longer lie on the geometric entity they are classified
  on.

  This function cannot be complex-stepped.  The complex part of the input
  vector will be ignored, and the output vector will have zero complex part.

  **Inputs**

   * mesh: DG mesh
   * xvec: a vector of mesh coordinates, perhaps obtained by
           [`coords3DTo1D`](@ref), length `mesh.dim*mesh.coord_numNodes`

  **Inputs/Outputs*

   * xivec: vector to be overwritten with parametric coordinates, length
            `mesh.geoNums.numXiDofs`
"""
function coords_xyzToXi(mesh::PumiMeshDG, xvec::AbstractVector,
                        xivec::AbstractVector)

  g = apf.getModel(mesh.m_ptr)
  if !gmi.can_eval(g)
    error("geometric model does not support evaluating parametric coordinates")
  end

  geonums = mesh.geoNums

  fill!(xivec, 0)
  xyz_j = zeros(Float64, 3)  # Pumi uses Float64 only
  xi_j = zeros(Float64, 3)

  newxyz_j = zeros(Float64, 3)


  for dim=0:mesh.dim

    it = apf.MeshIterator(mesh.m_ptr, dim)
    for i=1:mesh.numEntitiesPerType[dim+1]
      e = apf.iterate(mesh.m_ptr, it)
      me = apf.toModel(mesh.m_ptr, e)
      me_dim = apf.getModelType(mesh.m_ptr, me)

      for j=1:mesh.coord_numNodesPerType[dim+1]
        firstidx = apf.getNumberJ(geonums.xiNums, e, j-1, 0)

        # get xyz coordinates
        for k=1:mesh.dim
          idx = apf.getNumberJ(geonums.coordNums, e, j-1, k-1)
          xyz_j[k] = real(xvec[idx])
        end

        # decide parametric representation
        if me_dim == mesh.dim  # keep xyz coordinates
          for k=1:mesh.dim
            xi_j[k] = xyz_j[k]
          end
          
        elseif firstidx == geonums.numXiDofs + 1  # no geometric dofs
          for k=1:mesh.dim
            xi_j[k] = 0
          end
        else  # this is a constrained entity, find parametric representation
          gmi.closest_point(g, me, xyz_j, newxyz_j, sview(xi_j, 1:2))
#          @assert norm(newxyz_j - xyz_j) < 1e-8  # xyz_j should lie on this 
#                                                 # geometric entity
        end

        # put in output vector
        for k=1:me_dim
          idx = apf.getNumberJ(geonums.xiNums, e, j-1, k-1)
          xivec[idx] = xi_j[k]
        end

      end  # end j
    end  # end i

    apf.free(mesh.m_ptr, it)
  end  # end dim

  return nothing
end


"""
  Methods that takes a `mesh.vert_coords` shaped array as input, rather
  than a vector.

  **Inputs**

   * mesh
   * xarr: 3D array of xyz coordinates

  **Inputs/Outputs**

   * xivec: vector, length `mesh.geoNums,numXiDofs`
"""
function coords_xyzToXi(mesh::PumiMeshDG, xarr::Abstract3DArray{T},
                        xivec::AbstractVector) where {T}

  xvec = zeros(T, mesh.geoNums.numCoordDofs)
  op = AssignReduction{T}()
  coords3DTo1D(mesh, xarr, xvec, op, parallel=false)

  coords_xyzToXi(mesh, xvec, xivec)
end



"""
  This function is the inverse of [`coords_XiToXYZ`](@ref).  It recomputes
  the xyz coordinates from the parametric coordinates.  For MeshEntities which
  do not have any geometric degrees of freedom (and therefore no entires in
  `xivec`, the coordinates are retrieved from Pumi.

  **Inputs**

   * mesh
   * xivec: vector of CAD parametric coordinates, length `mesh.geoNums.numXiDofs`.
   
  **Inputs/Outputs**

   * xvec: vector to be overwritten with xyz coordinates, length `mesh.dim*mesh.coord_numNodes`.
"""
function coords_XiToXYZ(mesh::PumiMeshDG, xivec::AbstractVector,
                        xvec::AbstractVector)
  g = apf.getModel(mesh.m_ptr)
  if !gmi.can_eval(g)
    error("geometric model does not support evaluating parametric coordinates")
  end

  geonums = mesh.geoNums

  xyz_j = zeros(Float64, 3)
  xi_j = zeros(Float64, 2)
  xi_idx = zeros(Cint, mesh.dim)

  for dim=0:mesh.dim
    it = apf.MeshIterator(mesh.m_ptr, dim)
    for i=1:mesh.numEntitiesPerType[dim+1]
      e = apf.iterate(mesh.m_ptr, it)
      me = apf.toModel(mesh.m_ptr, e)
      me_dim = apf.getModelType(mesh.m_ptr, me)

      for j=1:mesh.coord_numNodesPerType[dim+1]

        for k=1:mesh.dim
          xi_idx[k] = apf.getNumberJ(geonums.xiNums, e, j-1, k-1)
        end

        # put the new coordinates in xyz_j
        if me_dim == mesh.dim  # xivec has xyz coordinates
          for k=1:mesh.dim
            xyz_j[k] = xivec[xi_idx[k]]
          end
        elseif xi_idx[1] == geonums.numXiDofs + 1  # no geometric dofs
          apf.getPoint(mesh.m_ptr, e, j-1, xyz_j)
        else  # constrained entity, convert parametric to xyz
          for k=1:me_dim
            xi_j[k] = xivec[xi_idx[k]]
          end
          gmi.geval(g, me, xi_j, xyz_j)
        end

        # write xyz_j to xvec
        for k=1:mesh.dim
          idx = apf.getNumberJ(geonums.coordNums, e, j-1, k-1)
          xvec[idx] = xyz_j[k]
        end

      end   # end j
    end  # end i

    apf.free(mesh.m_ptr, it)

  end  # end dim

  return nothing
end


"""
  Method that takes a 3D array of xyz coordinates as the output.

  Note that this method does not update the Pumi coordinate database nor
  does it recalculate the metrics.
  Users should call [`update_coords`](@ref) if this is required.

  **Inputs**

   * mesh
   * xivec: vector if CAD parametric coordinates, length `mesh.geoNums.numXiDofs`.


  **Inputs/Outputs**

   * xarr: array, same shape as `mesh.vert_coords`, to be overwritten with
           new coordinates.
"""
function coords_XiToXYZ(mesh::PumiMeshDG, xivec::AbstractVector,
                        xarr::Abstract3DArray{T}) where {T}

  xvec = zeros(T, mesh.geoNums.numCoordDofs)
  coords_XiToXYZ(mesh, xivec, xvec)

  coords1DTo3D(mesh, xvec, xarr)
end

