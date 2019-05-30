# interface functions for mapping between xyz and CAD parametric coordinates

"""
  This function gets the vector of CAD parametric coordinates from the
  mesh database.  See [`coords_xyzToXi`](@ref) for a description of the
  parametric coordinate system.

  **Inputs**

   * mesh

  **Inputs/Outputs**

   * xivec: vector to be overwritten with parametric coordinates.  Length
            `mesh.geoNums.numXiDofs`.
"""
function getXiCoords(mesh::PumiMeshDG, xivec::AbstractVector)
# retrieve Xi coordinates

  g = apf.getModel(mesh.m_ptr)
  if !gmi.can_eval(g)
    error("geometric model does not support evaluating parametric coordinates")
  end

  xi_j = zeros(Float64, 3)
  xyz_j = zeros(Float64, 3)

  for (e, dim) in apf.FieldEntityIt(mesh.m_ptr, mesh.coordshape_ptr)
    me = apf.toModel(mesh.m_ptr, e)
    me_dim = apf.getModelType(mesh.m_ptr, me)

    for j=1:mesh.coord_numNodesPerType[dim+1]
      firstidx = apf.getNumberJ(mesh.geoNums.xiNums, e, j-1, 0)


      # put the new coordinates in xyz_j
      if me_dim == mesh.dim  # xyz coordinates == xi coordinates
        getCoords(mesh, e, j-1, xi_j)
        ndof = mesh.dim
      elseif firstidx == mesh.geoNums.numXiDofs + 1  # no geometric dofs
        fill!(xi_j, 0)
        ndof = 0
      else  # constrained entity, convert parametric to xyz
        ndof = me_dim
        getCoordsXi(mesh, e, j-1, xi_j)
      end

      # write into output vector
      for k=1:ndof
        idx = apf.getNumberJ(mesh.geoNums.xiNums, e, j-1, k-1)
        xivec[idx] = xi_j[k]
      end
    end  # end j
  end  # end iterator

  return nothing
end

"""
  Similar to the other method, but allocates a new vector and returns it.

  **Inputs**

   * mesh

  **Outputs**

   * vector containing xi coordinates.  The vector has element type Float64
"""
function getXiCoords(mesh::PumiMeshDG)

  xivec = zeros(Float64, mesh.geoNums.numXiDofs)
  getXiCoords(mesh, xivec)

  return xivec
end


"""
  Gets the vector of xyz coordinates for the entire mesh.

  **Inputs**

   * mesh

  **Inputs/Outputs**

   * xvec: vector to be overwritten, length `mesh.dim*mesh.coord_numNodes`.
"""
function getXCoords(mesh::PumiMeshDG, xvec::AbstractVector)

  xyz_j = zeros(Float64, 3)
  for (e, dim) in apf.FieldEntityIt(mesh.m_ptr, mesh.coordshape_ptr)
    for j=1:mesh.coord_numNodesPerType[dim+1]

      getCoords(mesh, e, j-1, xyz_j)
      for k=1:mesh.dim
        idx = apf.getNumberJ(mesh.geoNums.coordNums, e, j-1, k-1)
        xvec[idx] = xyz_j[k]
      end
    end  # end j
  end # end iterator

  return nothing
end

"""
  Similar to the other method, but allocates a new vector and returns it.

  **Inputs**

   * mesh

  **Outputs**

   * vector containing the xyz coordinates.  Element type Float64.
"""
function getXCoords(mesh::PumiMeshDG)

  xvec = zeros(Float64, mesh.dim*mesh.coord_numNodes)
  getXCoords(mesh, xvec)

  return xvec
end

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

  This function has limited accuracy for entities classified on boundaries.
  TODO: what function to use instead.

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

  for (e, dim) in apf.FieldEntityIt(mesh.m_ptr, mesh.coordshape_ptr)
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
        getSnappedCoords(mesh, e, true, xyz_j, newxyz_j, xi_j)
#        gmi.closest_point(g, me, xyz_j, newxyz_j, sview(xi_j, 1:2))
#          @assert norm(newxyz_j - xyz_j) < 1e-8  # xyz_j should lie on this 
#                                                 # geometric entity
      end

      # put in output vector
      for k=1:me_dim
        idx = apf.getNumberJ(geonums.xiNums, e, j-1, k-1)
        xivec[idx] = xi_j[k]
      end

    end  # end j
  end  # end iterator

  return nothing
end


"""
  Methods that takes a `mesh.vert_coords` shaped array as input, rather
  than a vector.

  **Inputs**

   * mesh
   * xarr: 3D array of xyz coordinates

  **Inputs/Outputs**

   * xivec: vector, length `mesh.geoNums.numXiDofs`
"""
function coords_xyzToXi(mesh::PumiMeshDG, xarr::Abstract3DArray{T},
                        xivec::AbstractVector) where {T}

  xvec = zeros(T, mesh.geoNums.numCoordDofs)
  op = AssignReduction{T}()
  coords3DTo1D(mesh, xarr, xvec, op, parallel=false)

  coords_xyzToXi(mesh, xvec, xivec)
end


"""
  For non-periodic geometric entities, this function asserts the xi coordinates
  are in the allowed range.  If periodic geometric entities, this function wraps
  the coordinates to ensure they are in the correct range.  This function
  correctly preserves complex values if `xivec` is complex.

  This function is called by the other geometric interface functions when
  needed, users generally do not need to call it.

  **Inputs**

   * mesh

  **Inputs/Outputs**

   * xivec: vector of xi coordinates, length `mesh.geoNums.numXiDofs`.
"""
function wrapXiCoords(mesh::PumiMeshDG, xivec::AbstractVector)

  g = apf.getModel(mesh.m_ptr)
  if !gmi.can_eval(g)
    error("geometric model does not support evaluating parametric coordinates")
  end

  xi_idx = zeros(Int, 2)
  rng = zeros(Float64, 2)
  geonums = mesh.geoNums

  for (e, dim) in apf.FieldEntityIt(mesh.m_ptr, mesh.coordshape_ptr)
    me = apf.toModel(mesh.m_ptr, e)
    me_dim = apf.getModelType(mesh.m_ptr, me)

    for j=1:mesh.coord_numNodesPerType[dim+1]

      for k=1:2
        xi_idx[k] = apf.getNumberJ(geonums.xiNums, e, j-1, k-1)
      end

      # if not in geometric interior and has geometric degrees of freedom
      if !(me_dim == mesh.dim) && !(xi_idx[1] == geonums.numXiDofs + 1)
        for k=1:me_dim
          idx = xi_idx[k]
          gmi.range(g, me, k-1, rng)

          if gmi.periodic(g, me, k-1)
            xivec[idx] = wrapPeriodicXi(rng, xivec[idx])
          else
            @assert real(xivec[idx]) >= rng[1]
            @assert real(xivec[idx]) <= rng[2]
          end
        end  # end for
      end  # end if
    end  # end j
  end  # end iterator

  return nothing
end


"""
  Wraps a periodic xi value such that it represents the same location, but has
  value within the given range.  Methods are available for real and complex
  numbers.

  **Inputs**

   * rng: vector of length 2, containing the range of parametric values
   * xival: the xi value

  **Outputs**

   * xival: a new xi value
"""
function wrapPeriodicXi(rng::AbstractVector, xival::Real)

  # the algorithm is to shift to the range (0, rng[2] - rng[1]),
  # compute the modulus with the length of the parameter value range, then
  # shift back to the original range.
  delta = rng[2] - rng[1]
  t1 = (xival - rng[1]) % delta

  t1 += rng[1]

  if t1 < rng[1]
    t1 += delta
  end

  return t1
end


function wrapPeriodicXi(rng::AbstractVector, xival::T) where {T <: Complex}

  xival_imag = imag(xival)
  xival_real = wrapPeriodicXi(rng, real(xival))

  # because the geometric entity is periodic, wrapping by a multiple of the
  # period shouldn't change the derivative
  return T(xival_real, xival_imag)
end


"""
  This function is the inverse of [`coords_XiToXYZ`](@ref).  It recomputes
  the xyz coordinates from the parametric coordinates.  For MeshEntities which
  do not have any geometric degrees of freedom (and therefore no entries in
  `xivec`, the coordinates are retrieved from Pumi.

  **Inputs**

   * mesh
   * xivec: vector of CAD parametric coordinates, length
            `mesh.geoNums.numXiDofs`.  This vector may be modified if any
            periodic parametric values are out of range.  see
            [`wrapXiCoords`](@ref)
   
  **Inputs/Outputs**

   * xvec: vector to be overwritten with xyz coordinates, length `mesh.dim*mesh.coord_numNodes`.
"""
function coords_XiToXYZ(mesh::PumiMeshDG, xivec::AbstractVector,
                        xvec::AbstractVector)
  g = apf.getModel(mesh.m_ptr)
  if !gmi.can_eval(g)
    error("geometric model does not support evaluating parametric coordinates")
  end

  # make sure all values are valid
  wrapXiCoords(mesh, xivec)

  geonums = mesh.geoNums

  xyz_j = zeros(Float64, 3)
  xi_j = zeros(Float64, 2)
  xi_idx = zeros(Cint, mesh.dim)

  for (e, dim) in apf.FieldEntityIt(mesh.m_ptr, mesh.coordshape_ptr)
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
        getCoords(mesh, e, j-1, xyz_j)
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
  end  # end iterator

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


"""
  This function takes a vector that contains derivatives wrt to the xyz
  coordinates and transforms it to derivative wrt geometric coordinates xi.
  The derivative dx/dxi is evaluated at the x coordinates stored in the Pumi
  database.

  **Inputs**

   * mesh
   * xvec: vector containing dJ/dx

  **Inputs/Outputs**

   * xivec: vector to be overwritten with dJ/dxi
"""
function coords_dXTodXi(mesh::PumiMeshDG, xvec::AbstractVector,
                        xivec::AbstractVector)

  g = apf.getModel(mesh.m_ptr)
  if !gmi.can_eval(g)
    error("geometric model does not support evaluating parametric coordinates")
  end

  geonums = mesh.geoNums

  fill!(xivec, 0)
  dx_j = zeros(Float64, 3)
  dxi_j = zeros(Float64, 3)
  xidx = zeros(Cint, 3)

#  x_j = zeros(Float64, 3)
#  xnew_j = zeros(Float64, 3)
  xi_j = zeros(Float64, 2)
  dx_dxi = zeros(Float64, 3, 2)
  dx_dxi1 = sview(dx_dxi, :, 1)
  dx_dxi2 = sview(dx_dxi, :, 2)

  for (e, dim) in apf.FieldEntityIt(mesh.m_ptr, mesh.coordshape_ptr)
    me = apf.toModel(mesh.m_ptr, e)
    me_dim = apf.getModelType(mesh.m_ptr, me)

    for j=1:mesh.coord_numNodesPerType[dim+1]
      firstidx = apf.getNumberJ(geonums.xiNums, e, j-1, 0)

       # get xyz derivative
      for k=1:mesh.dim
        idx = apf.getNumberJ(geonums.coordNums, e, j-1, k-1)
        dx_j[k] = real(xvec[idx])
      end
        
      # decide parametric representation
      if me_dim == mesh.dim  # keep xyz coordinates
        for k=1:mesh.dim
          dxi_j[k] = dx_j[k]
        end
      elseif firstidx == geonums.numXiDofs + 1  # no geometric dofs
        for k=1:mesh.dim
          dxi_j[k] = 0
        end
      else  # this is a constrained entity, compute dx/dxi
        getCoordsXi(mesh, e, j-1, xi_j)
        gmi.first_derivative(g, me, xi_j, dx_dxi1, dx_dxi2)

        # apply chain rule to compute dJ/dxi
        fill!(dxi_j, 0)
        for k=1:me_dim
          for d=1:mesh.dim
            dxi_j[k] += dx_j[d]*dx_dxi[d, k]
          end
        end

      end  # end if
      
      # put in output vector
      for k=1:me_dim
        idx = apf.getNumberJ(geonums.xiNums, e, j-1, k-1)
        xivec[idx] = dxi_j[k]
      end

    end  # end j
  end  # end iterator

  return nothing
end


"""
  Takes a vector of the xyz coordinates for the entire mesh and updates the
  mesh with them.  Prefer the per-element method of this function when possible,
  this method is only preferable when the vector of all coordinates is required
  for some other reason.

  Unlike the per-element method, this function calls [`commit_coords`](@ref).

  Prefer the method that takes both `xvec` and `xivec` as arguments.  This
  method may lead to somewhat less accurate derivatives computed by
  [`coords_dXTodXi`](@ref).  This method is useful only when the caller
  has no interest in the parametric coordinates.

  **Inputs**

   * mesh
   * sbp
   * opts

  **Inputs/Outputs**

   * xvec: vector of new coordinates, length `mesh.dim*mesh.coord_numNodes`

  **Keyword Arguments**

   * snap: if true, the coordinates in `xvec` will be snapped to the
           nearest point on the geometry.  The updates coordinates will
           be saved to the mesh database and written back into `xvec`.
           If false, the coordinates in `xvec` will be entered into the
           mesh database without snapping.  This can lead to inconsistent
           parametric and Cartesian coordinates, as well as inaccurate
           derivatives.  Use only if the caller absolutely requires the
           coordinates be entered exactly as they appear in `xvec`
   * verify: run Pumi's verifier on the new mesh, default true
   * write_vis: write visulzation files after updating the coordinates but
                before recalculating the metrics (which errors if an element
                is turned inside out), default false.
"""
function update_coords(mesh::PumiMesh, sbp, opts, xvec::AbstractVector;
                       snap::Bool=true, verify=true, write_vis=false)

  @assert length(xvec) == mesh.dim*mesh.coord_numNodes

  xyz_j = zeros(Float64, 3)

  g = apf.getModel(mesh.m_ptr)
  _snap::Bool = snap && gmi.can_eval(g)   # make type-stable
  can_eval = gmi.can_eval(g)

  for (e, dim) in apf.FieldEntityIt(mesh.m_ptr, mesh.coordshape_ptr)
    me = apf.toModel(mesh.m_ptr, e)
    for j=1:mesh.coord_numNodesPerType[dim+1]

      for k=1:mesh.dim
        idx = apf.getNumberJ(mesh.coord_nodenums_Nptr, e, j-1, k-1)
        xyz_j[k] = xvec[idx]
      end

      setCoords(mesh, e, j-1, xyz_j, _snap)
      # write back into xvec
      for k=1:mesh.dim
        idx = apf.getNumberJ(mesh.coord_nodenums_Nptr, e, j-1, k-1)
        xvec[idx] = xyz_j[k]
      end
    end  # end j
  end  # end iterator


  commit_coords(mesh, sbp, opts, verify=verify, write_vis=write_vis)

  return nothing
end

#TODO: single elemet version of update_coords

# It looks like the MeshAdapt Snapper updates the parametric coords for
# vertices, but not mid edge nodes (I also recall that MeshAdapt doesn't
# support quadratic triangles either).

"""
  Convenience method that computes the xyz coordinates from the parametric
  coordinates and updates the mesh database.  This method is more accurate than
  the function that updates the mesh database using the xyz coordinates

  **Inputs**
  
   * mesh
   * sbp
   * opts
   * xivec: vector of xi coordinates, length `mesh.geoNums.numXiDofs`

  **Keyword Arguments**

   * verify: run Pumi's verifier on the new mesh, default true
   * write_vis: write visulzation files after updating the coordinates but
                before recalculating the metrics (which errors if an element
                is turned inside out), default false.
"""
function update_coordsXi(mesh::PumiMesh, sbp, opts, xivec::AbstractVector;
                       verify=true, write_vis=false)
  @assert length(xivec) == mesh.geoNums.numXiDofs

  xi_j = zeros(Float64, 3)
  geonums = mesh.geoNums
  for (e, dim) in apf.FieldEntityIt(mesh.m_ptr, mesh.coordshape_ptr)

    me = apf.toModel(mesh.m_ptr, e)
    me_dim = apf.getModelType(mesh.m_ptr, me)

    for j=1:mesh.coord_numNodesPerType[dim+1]

      # decide parametric representation
      firstidx = apf.getNumberJ(geonums.xiNums, e, j-1, 0)
      if me_dim == mesh.dim
        # parametric coords == xyz coordinates 
        fill!(xi_j, 0)
        for k=1:me_dim
          xi_j[k] = xivec[apf.getNumberJ(geonums.xiNums, e, j-1, k-1)]
        end
        setCoords(mesh, e, j-1, xi_j)
      elseif firstidx == geonums.numXiDofs + 1  # no geometric dofs
        # this must be classified on a vertex, which always has parametric
        # coordinates 0
        fill!(xi_j, 0)
      else  # this is a constrained entity, get the values from xivec
        fill!(xi_j, 0)
        for k=1:me_dim
          xi_j[k] = xivec[apf.getNumberJ(geonums.xiNums, e, j-1, k-1)]
        end
        setCoordsXi(mesh, e, j-1, xi_j)
      end  # end if
    end  # end j
  end  # end iterator

  commit_coords(mesh, sbp, opts, verify=verify, write_vis=write_vis)

  return nothing
end
