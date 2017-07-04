# functions for calculating node coordinates and metrics for curvilinear
# meshes

#------------------------------------------------------------------------------
# allocators for curvilinear meshes

"""
  Allocates mesh.vert_coords
"""
function allocateMeshCoordinateArray{Tmsh}(mesh::PumiMeshDG{Tmsh}, sbp::AbstractSBP)

  num_coord_nodes = mesh.coord_numNodesPerElement

  if !isFieldDefined(mesh, :vert_coords)
    mesh.vert_coords = Array(Float64, mesh.dim, num_coord_nodes, mesh.numEl)
  else
    fill!(mesh.vert_coords, 0.0)
  end

  return nothing
end

"""
  Allocates the following fields of mesh:
    coords
    dxidx
    jac

  Also allocates the following fields with size zero (because they are not used):
   dxidx_face
   jac_face
   dxidx_bndry
   jac_bndry
   dxidx_sharedface
   jac_sharedface
"""
function allocateCurvilinearCoordinateAndMetricArrays{Tmsh}(mesh::PumiMeshDG{Tmsh},                                                         sbp::AbstractSBP)

  dim = mesh.dim
  sbpface = mesh.sbpface
  
  if !isFieldDefined(mesh, :coords, :dxidx, :jac)
    mesh.coords = zeros(Tmsh, mesh.dim, mesh.numNodesPerElement, mesh.numEl)
    mesh.dxidx = zeros(Tmsh, mesh.dim, mesh.dim, mesh.numNodesPerElement, mesh.numEl)
    mesh.jac = zeros(Tmsh, mesh.numNodesPerElement, mesh.numEl)

    # these arrays are not used for curvilinear meshes

    # interior faces
    mesh.dxidx_face = zeros(Tmsh, 0, 0, 0, 0)
    mesh.jac_face = zeros(Tmsh, 0, 0)

    # boundary faces
    mesh.dxidx_bndry = zeros(Tmsh, 0, 0, 0, 0)
    mesh.jac_bndry = zeros(Tmsh, 0, 0)

    # parallel shared faces
    mesh.dxidx_sharedface = Array(Array{Tmsh, 4}, 0)
    mesh.jac_sharedface = Array(Array{Tmsh, 2}, 0)

  else
    fill!(mesh.coords, 0.0)
    fill!(mesh.dxidx, 0.0)
    fill!(mesh.jac, 0.0)
  end


  return nothing
end

#------------------------------------------------------------------------------
# curvilinear functions

"""
  This function gets the coordinates that define the mesh 
  (not the solution node coordinates) and puts them the mesh.vert_coords field
  of the mesh object.  For each element, the ordering of the coordinates
  is verts, then edges, then faces, then regions.  This function uses
  smart allocators to allocate the array if needed.

  Input:
    mesh: a DG mesh
    sbp: an SBP operators
"""
function getMeshCoordinates{Tmsh}(mesh::PumiMeshDG{Tmsh}, sbp::AbstractSBP)

  allocateMeshCoordinateArray(mesh, sbp)

  for i=1:mesh.numEl
    el_i = mesh.elements[i]
    coords_i = sview(mesh.vert_coords, :, :, i)
    getAllEntityCoords(mesh.m_ptr, el_i, coords_i)
  end

  return nothing
end

"""
  This function calculates the fields of the mesh that hold coordinates of the
  face nodes for boundaries, interfaces, and sharedfaces.  This function uses
  smart allocators to allocate the arrays if needed
"""
function getFaceCoordinatesAndNormals{Tmsh}(mesh::PumiMeshDG{Tmsh}, sbp::AbstractSBP)

  allocateFaceCoordinates(mesh)

  allocateNormals(mesh, sbp)

  if length(mesh.bndryfaces) > 0  # debugging: don't call if unneeded
    calcFaceCoordinatesAndNormals(mesh, sbp, mesh.bndryfaces, 
                                  mesh.coords_bndry, mesh.nrm_bndry)
  end

  calcFaceCoordinatesAndNormals(mesh, sbp, mesh.interfaces, 
                               mesh.coords_interface, mesh.nrm_face)
  for i=1:mesh.npeers
    calcFaceCoordinatesAndNormals(mesh, sbp, mesh.bndries_local[i], 
                               mesh.coords_sharedface[i], mesh.nrm_sharedface[i])
  end

  return nothing
end

"""
  This function back propigates mesh.nrm_*_bar and mesh.coords_*_bar to
  mesh.vert_coords_bar.  mesh.vert_coords_bar is updated (not overwritten)
  with the results.
"""
function getFaceCoordinatesAndNormals_rev{Tmsh}(mesh::PumiMeshDG{Tmsh},
                                                sbp::AbstractSBP)

  if length(mesh.bndryfaces) > 0  # debugging: don't call if unneeded
    #TODO: don't allocate this every time
    coords_bndry_bar = zeros(mesh.dim, mesh.numNodesPerFace, mesh.numBoundaryFaces)
    calcFaceCoordinatesAndNormals_rev(mesh, sbp, mesh.bndryfaces, 
                                      coords_bndry_bar,
                                      mesh.nrm_bndry,
                                      mesh.nrm_bndry_bar)
  end


  coords_face_bar = zeros(mesh.dim, mesh.numNodesPerFace, mesh.numInterfaces)
  calcFaceCoordinatesAndNormals_rev(mesh, sbp, mesh.interfaces, 
                                coords_face_bar, mesh.nrm_face,
                                mesh.nrm_face_bar)
  for i=1:mesh.npeers
    coords_sharedface_bar = zeros(mesh.dim, mesh.numNodesPerFace, length(mesh.bndries_local[i]))
    calcFaceCoordinatesAndNormals_rev(mesh, sbp, mesh.bndries_local[i], 
                                  coords_sharedface_bar[i], 
                                  mesh.nrm_sharedface[i],
                                  mesh.nrm_sharedface_bar[i])
  end

  return nothing
end

"""
  This is the inner function used by getFaceCoordinatesandNormals.  It
  gets the coordinates of the face nodes and their normal vectors for
  a given array of Interfaces/Boundaries

  Inputs:
    mesh: a DG mesh object
    faces: an array of Boundaries or Interfaces

  Inputs/Outputs
    coords_face: array Tdim x numNodesPerFace x length(faces)
    nrm_face: array, same size as coords_face

  Aliasing restrictions: none, although it would be weird if coords_face and
                         nrm_face aliased
"""
function calcFaceCoordinatesAndNormals{Tmsh, I <: Union{Boundary, Interface}}(
                    mesh::PumiMeshDG{Tmsh}, sbp::AbstractSBP,
                    faces::AbstractArray{I, 1}, 
                    coords_face::AbstractArray{Tmsh, 3}, 
                    nrm_face::AbstractArray{Tmsh, 3})

  blocksize = 1000  # magic parameter: number of faces to do in a group
  nfaces = length(faces)

  # calculate number of blocks
  nblocks_full = div(nfaces, blocksize)
  nrem = nfaces % blocksize

  numNodesPerElement = mesh.coord_numNodesPerElement
  numNodesPerFace = mesh.coord_numNodesPerType[mesh.dim]

  # some temporary arrays
  down_faces = Array(Ptr{Void}, 12)
  coords_lag_face = Array(Float64, mesh.dim, mesh.coord_numNodesPerFace, blocksize)

  # get the parametic coordinates of the face nodes
  face_xi = mesh.coord_facexi
  ref_verts = baryToXY(face_xi, mesh.sbpface.vtx)

  face_idx = 1  # index in faces array
  for block=1:nblocks_full
    start_idx = (block - 1)*blocksize + 1
    end_idx = block*blocksize
    coords_face_block = sview(coords_face, :, :, start_idx:end_idx)
    nrm_face_block = sview(nrm_face, :, :, start_idx:end_idx)

    # load up data for current block
    for i=1:blocksize
      el_i = getElementL(faces[face_idx])
      face_i = getFaceL(faces[face_idx])
      el_ptr = mesh.elements[el_i]

      coords_i = sview(coords_lag_face, :, :, i)
      getMeshFaceCoordinates(mesh, el_i, face_i, coords_i)
      face_idx += 1
    end

    # populate output array for current block
    calcFaceNormals!(mesh.sbpface, mesh.coord_order, ref_verts, coords_lag_face,
                     coords_face_block, nrm_face_block)
    fill!(coords_lag_face, 0.0)

  end  # end loop over full blocks

  # do remainder loop
  start_idx = nblocks_full*blocksize + 1
  end_idx = nfaces
  @assert end_idx - start_idx + 1 <= blocksize
  @assert end_idx - start_idx + 1 == nrem

  coords_face_block = sview(coords_face, :, :, start_idx:end_idx)
  nrm_face_block = sview(nrm_face, :, :, start_idx:end_idx)
  coords_lag_face_block = sview(coords_lag_face, :, :, 1:nrem)

  for i=1:nrem  # loop over remaining faces
    el_i = getElementL(faces[face_idx])
    face_i = getFaceL(faces[face_idx])
    el_ptr = mesh.elements[el_i]

    coords_i = sview(coords_lag_face, :, :, i)
    getMeshFaceCoordinates(mesh, el_i, face_i, coords_i)
    face_idx += 1


  end

  calcFaceNormals!(mesh.sbpface, mesh.coord_order, ref_verts, 
                   coords_lag_face_block, coords_face_block, nrm_face_block)


  # make sure the normal vectors point outwards
  fixOutwardNormal(mesh, faces, nrm_face)

  return nothing
end

"""
  Reverse mode of calcFaceCoordinatesAndNormals_rev, back propigates
  coords_face_bar and nrm_face_bar to mesh.coords_bar_ptr (which is
  incremented, not overwritten).
"""
function calcFaceCoordinatesAndNormals_rev{Tmsh, I <: Union{Boundary, Interface}}(
                    mesh::PumiMeshDG{Tmsh}, sbp::AbstractSBP,
                    faces::AbstractArray{I, 1}, 
                    coords_face_bar::AbstractArray{Tmsh, 3},
                    nrm_face::AbstractArray{Tmsh, 3},
                    nrm_face_bar::AbstractArray{Tmsh, 3})

  blocksize = 1000  # magic parameter: number of faces to do in a group
  nfaces = length(faces)

  # calculate number of blocks
  nblocks_full = div(nfaces, blocksize)
  nrem = nfaces % blocksize

  numNodesPerElement = mesh.coord_numNodesPerElement
  numNodesPerFace = mesh.coord_numNodesPerType[mesh.dim]

  # some temporary arrays
  down_faces = Array(Ptr{Void}, 12)
  coords_lag_face = Array(Float64, mesh.dim, mesh.coord_numNodesPerFace, blocksize)
  coords_lag_face_bar = zeros(coords_lag_face)

  # get the parametic coordinates of the face nodes
  face_xi = mesh.coord_facexi
  ref_verts = baryToXY(face_xi, mesh.sbpface.vtx)

  for block=1:nblocks_full
    start_idx = (block - 1)*blocksize + 1
    end_idx = block*blocksize
    faces_block = sview(faces, start_idx:end_idx)
    coords_face_bar_block = sview(coords_face_bar, :, :, start_idx:end_idx)
    nrm_face_bar_block = sview(nrm_face_bar, :, :, start_idx:end_idx)

    # load up data for current block
    for i=1:blocksize
      el_i = getElementL(faces_block[i])
      face_i = getFaceL(faces_block[i])
      el_ptr = mesh.elements[el_i]

      coords_i = sview(coords_lag_face, :, :, i)
      getMeshFaceCoordinates(mesh, el_i, face_i, coords_i)
    end

    # populate output array for current block
    fill!(coords_lag_face_bar, 0.0)
    calcFaceNormals_rev!(mesh.sbpface, mesh.coord_order, ref_verts,
                         coords_lag_face, coords_lag_face_bar,
                         coords_face_bar_block, nrm_face_bar_block)

    fill!(coords_lag_face, 0.0)

    # write coords_lag_face_bar to Pumi field (to be retrieved later)
    for i=1:blocksize
      el_i = getElementL(faces_block[i])
      face_i = getFaceL(faces_block[i])
      el_ptr = mesh.elements[el_i]

      coords_bar_i = sview(coords_lag_face_bar, :, :, i)
      getMeshFaceCoordinates_rev(mesh, el_i, face_i, coords_bar_i)
    end

  end  # end loop over full blocks

  # do remainder loop
  start_idx = nblocks_full*blocksize + 1
  end_idx = nfaces
  @assert end_idx - start_idx + 1 <= blocksize
  @assert end_idx - start_idx + 1 == nrem

  faces_block = sview(faces, start_idx:end_idx)
  coords_face_bar_block = sview(coords_face_bar, :, :, start_idx:end_idx)
  nrm_face_bar_block = sview(nrm_face_bar, :, :, start_idx:end_idx)
  coords_lag_face_block = sview(coords_lag_face, :, :, 1:nrem)
  coords_lag_face_bar_block = sview(coords_lag_face_bar, :, :, 1:nrem)

  for i=1:nrem  # loop over remaining faces
    el_i = getElementL(faces_block[i])
    face_i = getFaceL(faces_block[i])
    el_ptr = mesh.elements[el_i]

    coords_i = sview(coords_lag_face, :, :, i)
    getMeshFaceCoordinates(mesh, el_i, face_i, coords_i)
  end

  calcFaceNormals_rev!(mesh.sbpface, mesh.coord_order, ref_verts, 
                   coords_lag_face_block, coords_lag_face_bar_block,
                   coords_face_bar_block, nrm_face_bar_block)


  # write coords_lag_face_bar to Pumi field (to be retrieved later)
  for i=1:nrem
    el_i = getElementL(faces_block[i])
    face_i = getFaceL(faces_block[i])
    el_ptr = mesh.elements[el_i]

    coords_bar_i = sview(coords_lag_face_bar_block, :, :, i)
    getMeshFaceCoordinates_rev(mesh, el_i, face_i, coords_bar_i)
  end


  fixOutwardNormal_rev(mesh, faces, nrm_face, nrm_face_bar)

  return nothing
end



"""
  This function check to make sure each face normal vector is oriented
  outwards, and flips it if needed.

  Inputs:
    mesh
    faces: array of Boundary of Interfaces to check for normal orientation

  Inputs/Outputs:
    nrm_face: array of size mesh.dim x mesh.numNodesPerFace x length(faces)
              containing the normal vector at each node of each face.
              Updated in place to make normal vector point outward
"""
function fixOutwardNormal{I <: Union{Boundary, Interface}, Tmsh}(mesh, 
                          faces::AbstractArray{I, 1},
                          nrm_face::AbstractArray{Tmsh, 3})


  nfaces = length(faces)
  should_flip_node = Array(Bool, mesh.numNodesPerFace)
  for i=1:nfaces
    is_inward_normal(mesh, faces[i], sview(nrm_face, :, :, i), should_flip_node)

    for j=1:mesh.numNodesPerFace
      if should_flip_node[j]
        for p=1:mesh.dim
          nrm_face[p, j, i] = -nrm_face[p, j, i]
        end
      end
    end  # end j

  end  # end loop i

   return nothing
end

"""
  Reverse mode of fixOutwardNormal.

  Inputs:
    mesh
    faces
    nrm_face

  Inputs/Outputs:
    nrm_face_bar: adjoint part of nrm_face
"""
function fixOutwardNormal_rev{Tmsh, I <: Union{Boundary, Interface}}(mesh,
                          faces::AbstractArray{I, 1},
                          nrm_face::AbstractArray{Tmsh, 3},
                          nrm_face_bar::AbstractArray{Tmsh, 3})

  nfaces = length(faces)
  should_flip_node = Array(Bool, mesh.numNodesPerFace)
  for i=1:nfaces
    is_inward_normal(mesh, faces[i], sview(nrm_face, :, :, i), should_flip_node)

    for j=1:mesh.numNodesPerFace
      if should_flip_node[j]
        for p=1:mesh.dim
          nrm_face_bar[p, j, i] = nrm_face_bar[p, j, i]
        end
      end
    end  # end j

  end  # end loop i

  return nothing
end



"""
  Check if the normal vector on each node of a face is pointing inward

  The algorithm determines orientation by comparing the normal vector against
  the vector along an edge of the simplex (from the vertex not on the face to
  the vertex on the face).  For highly skewed elements, using different
  vertices on the face to compute the vector can give different results, so
  only the vector with the largest magnitude dot product with the normal vector
  is considered.

  Inputs:
    mesh: the mesh
    iface: either a Boundary or an Interface
    nrm_face: mesh.dim x mesh.numNodesPerFace array containing the normal
              vectors for each face node

  Inputs/Outputs:
    should_flip_node: an Bool array of length numNodesPerFace specifying whether
                      the normal vector at each node is pointing inward
"""
function is_inward_normal(mesh, iface::Union{Boundary, Interface},
                          nrm_face::AbstractMatrix,
                          should_flip_node::AbstractVector{Bool})

  tmp = zeros(3)  # temporary vector to hold coordinates
  topo = mesh.topo
  numVertPerElement = mesh.numTypePerElement[1]
  numVertPerFace = numVertPerElement - 1

  # temporary arrays
  el_verts = Array(Ptr{Void}, numVertPerElement)
  other_vert_coords = zeros(mesh.dim)
  face_verts = Array(Ptr{Void}, numVertPerElement - 1)
  face_vert_coords = zeros(mesh.dim, numVertPerFace)

  elnum = getElementL(iface)
  facenum_local = getFaceL(iface)

  el_i = mesh.elements[elnum]
  getDownward(mesh.m_ptr, el_i, 0, el_verts)

  for j=1:numVertPerFace
    face_verts[j] = el_verts[topo.face_verts[j, facenum_local]]
    getPoint(mesh.m_ptr, face_verts[j], 0, tmp)

    for p=1:mesh.dim
      face_vert_coords[p, j] = tmp[p]
    end
  end

  # get the vert not on the face
  other_vert = Ptr{Void}(0)
  for j=1:numVertPerElement
    if !(el_verts[j] in face_verts)
      other_vert = el_verts[j]
    end
  end

  getPoint(mesh.m_ptr, other_vert, 0, tmp)
  for p=1:mesh.dim
    other_vert_coords[p] = tmp[p]
  end

  # check that the face normal is in the opposite direction as the
  # vectors from a vertex on the face to the vertex not on the face

  # in some cases the normal vector can be nearly orthoogonal to vertex
  # vectors, so use the one with the greatest magnitude dot product
  should_flip = false
  max_mag = 0.0  # maximum dot product
  for j=1:mesh.numNodesPerFace
    outward_count = 0  # count number of calculations that showed outward
    for k=1:numVertPerFace
      val = zero(Float64)
      for p=1:mesh.dim
        r1_p = other_vert_coords[p] - face_vert_coords[p, k]
        val += nrm_face[p, j]*r1_p  # accumulate dot product
      end
      
      if abs(val) > max_mag
        max_mag = abs(val)
        # flip if the value is greater than 0
        should_flip = val > 0
        #=
        if val < 0
          should_flip = false
        else 
          should_flip = true
        end  # end if
        =#

      end  
    end  # end k

    should_flip_node[j] = should_flip
  end  # end loop j

  return nothing
end




"""
  This function calculates the node coordinates, scaled mapping jacobian, and
  mapping jacobian determinant for a curvilinear mesh and stores them to
  the fields of the mesh object.  This function uses smart allocators to
  allocate the arrays if needed

  Inputs:
    mesh: a DG mesh object
    sbp: an SBP operator
"""
function getCurvilinearCoordinatesAndMetrics{Tmsh}(mesh::PumiMeshDG{Tmsh}, 
                                                  sbp::AbstractSBP)

  allocateCurvilinearCoordinateAndMetricArrays(mesh, sbp)

  ref_vtx = baryToXY(mesh.coord_xi, sbp.vtx)
#  if mesh.dim == 2
#    calcMappingJacobian!(sbp, mesh.coord_order, ref_vtx, mesh.vert_coords, 
#                         mesh.coords, mesh.dxidx, mesh.jac)
#  else  # need to calculate Eone

    # block format
    blocksize = 1000  # number of elements per block
    nblocks_full = div(mesh.numEl, blocksize)
    nrem = mesh.numEl % blocksize

    Eone = zeros(Tmsh, mesh.numNodesPerElement, mesh.dim, blocksize)

    for block=1:nblocks_full
      start_idx = (block - 1)*blocksize + 1
      end_idx = block*blocksize
      element_range = start_idx:end_idx
      
      getCurvilinearMetricsAndCoordinates_inner(mesh, sbp, element_range, Eone)
     
    end  # end loop over blocks

    if nrem != 0
      # remainder block
      start_idx = nblocks_full*blocksize + 1
      end_idx = mesh.numEl
      @assert end_idx - start_idx + 1 <= blocksize
      @assert end_idx - start_idx + 1 == nrem

      element_range = start_idx:end_idx

      # make sure dimensions of Eone is correct (number of elements might be
      # less than before)
      Eone_rem = sview(Eone, :, :, 1:nrem)

      getCurvilinearMetricsAndCoordinates_inner(mesh, sbp, element_range, Eone_rem)
    end
#  end  # end if dim == 2

  return nothing
end

"""
  Back propigate dxidx to mesh.vert_coords and mesh.nrm.  See the primal method
  for details.

  Inputs:
    mesh: vert_coords_bar, nrm_face_bar, nrm_bndry_bar, nrm_sharedface are
          updated with the results
    sbp:
"""
function getCurvilinearCoordinatesAndMetrics_rev{Tmsh}(mesh::PumiMeshDG{Tmsh},
                                                       sbp::AbstractSBP)

  # we have to compute E1_bar for both 2d and 3D
  blocksize = 1000  # number of elements per block
  nblocks_full = div(mesh.numEl, blocksize)
  nrem = mesh.numEl % blocksize

  Eone_bar = zeros(Tmsh, mesh.numNodesPerElement, mesh.dim, blocksize)

  for block=1:nblocks_full
    start_idx = (block - 1)*blocksize + 1
    end_idx = block*blocksize
    element_range = start_idx:end_idx
    
    getCurvilinearMetricsAndCoordinates_inner_rev(mesh, sbp, element_range,
                                                  Eone_bar)
  end

 if nrem != 0
    # remainder block
    start_idx = nblocks_full*blocksize + 1
    end_idx = mesh.numEl
    @assert end_idx - start_idx + 1 <= blocksize
    @assert end_idx - start_idx + 1 == nrem

    element_range = start_idx:end_idx

    # make sure dimensions of Eone is correct (number of elements might be
    # less than before)
    Eone_bar_rem = sview(Eone_bar, :, :, 1:nrem)

    getCurvilinearMetricsAndCoordinates_inner_rev(mesh, sbp, element_range,
                                                  Eone_bar_rem)
  end

  return Eone_bar
end


"""
  This function calculates the metrics and coordinates for one block of 
  elements.  Used by getCurvilinearMetricsAndCoordinates
"""
function getCurvilinearMetricsAndCoordinates_inner{T}(mesh, sbp, 
                             element_range::UnitRange, Eone::AbstractArray{T, 3})

  # calculate Eone for current range

  vert_coords_block = sview(mesh.vert_coords, :, :, element_range)
  coords_block = sview(mesh.coords, :, :, element_range)
  dxidx_block = sview(mesh.dxidx, :, :, :, element_range)
  jac_block = sview(mesh.jac, :, element_range)
  ref_vtx = baryToXY(mesh.coord_xi, sbp.vtx)

  calcEone(mesh, sbp, element_range, Eone)

  println("checking Eone")
  for el=1:size(Eone, 3)
    println("el = ", el)
    println("Eone = \n", Eone[:, :, el])
    for d=1:mesh.dim
      println("  d = ", d)
      val = abs(sum(Eone[:, d, el]))
      println("  val = ", val)
      @assert val < 1e-14
    end
  end
  calcMappingJacobian!(sbp, mesh.coord_order, ref_vtx, vert_coords_block, 
                       coords_block, dxidx_block, jac_block, Eone)

  fill!(Eone, 0.0)
  return nothing
end

"""
  This function calculates the metrics and coordinates for one block of 
  elements.  Used by getCurvilinearMetricsAndCoordinates
"""
function getCurvilinearMetricsAndCoordinates_inner_rev{T}(mesh, sbp, 
                             element_range::UnitRange, Eone_bar::AbstractArray{T, 3})


  vert_coords_block = sview(mesh.vert_coords, :, :, element_range)
  # we don't allow perturbing mesh coordinate directly (yet)
  coords_bar_block = zeros(T, mesh.dim, mesh.numNodesPerElement, length(element_range))
  dxidx_block = sview(mesh.dxidx, :, :, :, element_range)
  dxidx_bar_block = sview(mesh.dxidx_bar, :, :, :, element_range)
#  jac_block = sview(mesh.jac, :, element_range)
  jac_bar_block = sview(mesh.jac_bar, :, element_range)
  @assert vecnorm(jac_bar_block) < 1e-13  # this is currently broken in SBP

  # outputs
  vert_coords_bar_block = sview(mesh.vert_coords_bar, :, :, element_range)
  fill!(Eone_bar, 0.0)

  ref_vtx = baryToXY(mesh.coord_xi, sbp.vtx)

  # back propigate dxidx to vert_coords, E1
  calcMappingJacobian_rev!(sbp, mesh.coord_order, ref_vtx, vert_coords_block, 
                           vert_coords_bar_block, coords_bar_block,
                           dxidx_bar_block, jac_bar_block, Eone_bar)


  # back propigate E1 to the face normals
  calcEone_rev(mesh, sbp, element_range, Eone_bar)

  return nothing
end

function calcEone{Tmsh}(mesh::PumiMeshDG{Tmsh}, sbp, element_range, 
                        Eone::AbstractArray{Tmsh, 3})

  # search interfaces and bndryfaces

  first_el = element_range[1]
  last_el = element_range[end]

  offset = first_el - 1

  # R times vector of ones
  # the permutation doesn't matter because it is being multiplied by a constant
  # vector.
  # also R1 = 1 by definition
  Rone = vec(sum(mesh.sbpface.interp.', 2))
  sbpface = mesh.sbpface
  tmp = zeros(Rone)
  nrmL = zeros(Tmsh, mesh.dim, sbpface.numnodes)

  # accumulate E1 for a given element
  Eone_el = zeros(Tmsh, sbpface.stencilsize, mesh.dim)

  for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]

    if iface_i.elementL >= first_el && iface_i.elementL <= last_el
      elnum = iface_i.elementL
      facenum_local = iface_i.faceL
      nrm = sview(mesh.nrm_face, :, :, i)

      # call inner functions
      calcEoneElement(sbpface, nrm, Rone, tmp, Eone_el)
      assembleEone(sbpface, elnum - offset, facenum_local, Eone_el, Eone)

    end

    if iface_i.elementR >= first_el && iface_i.elementR <= last_el
      # do same as above, negating and permuting nrm
      elnum = iface_i.elementR
      facenum_local = iface_i.faceR
      orient = iface_i.orient
      for j=1:sbpface.numnodes
        for d=1:mesh.dim
          nrmL[d, j] = -mesh.nrm_face[d, sbpface.nbrperm[j, orient], i]
        end
      end
 
      calcEoneElement(sbpface, nrmL, Rone, tmp, Eone_el)
      assembleEone(sbpface, elnum - offset, facenum_local, Eone_el, Eone)
    end  # end if/else

  end  # end loop over interfaces

  for i=1:mesh.numBoundaryFaces
    bface_i = mesh.bndryfaces[i]

    if bface_i.element >= first_el && bface_i.element <= last_el
      elnum = bface_i.element
      facenum_local = bface_i.face
      nrm = sview(mesh.nrm_bndry, :, :, i)

      # call inner functions
      calcEoneElement(sbpface, nrm, Rone, tmp, Eone_el)
      assembleEone(sbpface, elnum - offset, facenum_local, Eone_el, Eone)
    end
  end  # end loop over boundary faces

  # check shared faces
  for peer=1:mesh.npeers
    bndryfaces_peer = mesh.bndries_local[peer]
    nrm_peer = mesh.nrm_sharedface[peer]
    for i=1:mesh.peer_face_counts[peer]
      bface_i = bndryfaces_peer[i]

      if bface_i.element >= first_el && bface_i.element <= last_el
        elnum = bface_i.element
        facenum_local = bface_i.face
        nrm = sview(nrm_peer, :, :, i)

        calcEoneElement(sbpface, nrm, Rone, tmp, Eone_el)
        assembleEone(sbpface, elnum - offset, facenum_local, Eone_el, Eone)
      end
    end
  end



  return nothing
end

"""
  Back propigates Eone_bar to the various normal vectors stored in the
  mesh objects (bndry, face, sharedface)

  Inputs:
    mesh: mesh.nrm_face_bar, mesh.nrm_bndry_bar, mesh.nrm_sharedface_bar are
          updated
    sbp
    element range: range of elements to compute
    Eone_bar: the adjoint part of Eone, see the primal method for details
"""
function calcEone_rev{Tmsh}(mesh::PumiMeshDG{Tmsh}, sbp, element_range, 
                            Eone_bar::AbstractArray{Tmsh, 3})

  # search interfaces and bndryfaces

  first_el = element_range[1]
  last_el = element_range[end]

  offset = first_el - 1

  # R times vector of ones
  # the permutation doesn't matter because it is being multiplied by a constant
  # vector.
  # also R1 = 1 by definition
  Rone = vec(sum(mesh.sbpface.interp.', 2))
  sbpface = mesh.sbpface
  tmp = zeros(Rone)
  nrmL_bar = zeros(Tmsh, mesh.dim, sbpface.numnodes)

  # accumulate E1 for a given element
  Eone_el_bar = zeros(Tmsh, sbpface.stencilsize, mesh.dim)

  for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]

    if iface_i.elementL >= first_el && iface_i.elementL <= last_el
      elnum = iface_i.elementL
      facenum_local = iface_i.faceL
      nrm_bar = sview(mesh.nrm_face_bar, :, :, i)

      # call inner functions
      fill!(Eone_el_bar, 0.0)
      assembleEone_rev(sbpface, elnum - offset, facenum_local, Eone_el_bar,
                       Eone_bar)
      calcEoneElement_rev(sbpface, nrm_bar, Rone, tmp, Eone_el_bar)

    end

    if iface_i.elementR >= first_el && iface_i.elementR <= last_el
      # do same as above, negating and permuting nrm
      elnum = iface_i.elementR
      facenum_local = iface_i.faceR
      orient = iface_i.orient

      fill!(Eone_el_bar, 0.0)
      assembleEone_rev(sbpface, elnum - offset, facenum_local, Eone_el_bar,
                       Eone_bar)
      calcEoneElement_rev(sbpface, nrmL_bar, Rone, tmp, Eone_el_bar)

      for j=1:sbpface.numnodes
        for d=1:mesh.dim
          mesh.nrm_face_bar[d, sbpface.nbrperm[j, orient], i] -= nrmL_bar[d, j]
        end
      end
 
    end  # end if/else

  end  # end loop over interfaces

  for i=1:mesh.numBoundaryFaces
    bface_i = mesh.bndryfaces[i]

    if bface_i.element >= first_el && bface_i.element <= last_el
      elnum = bface_i.element
      facenum_local = bface_i.face
      nrm_bar = sview(mesh.nrm_bndry_bar, :, :, i)

      # call inner functions
      fill!(Eone_el_bar, 0.0)
      assembleEone_rev(sbpface, elnum - offset, facenum_local, Eone_el_bar, 
                       Eone_bar)
      calcEoneElement_rev(sbpface, nrm_bar, Rone, tmp, Eone_el_bar)
    end
  end  # end loop over boundary faces

  # check shared faces
  for peer=1:mesh.npeers
    bndryfaces_peer = mesh.bndries_local[peer]
    nrm_peer_bar = mesh.nrm_sharedface_bar[peer]
    for i=1:mesh.peer_face_counts[peer]
      bface_i = bndryfaces_peer[i]

      if bface_i.element >= first_el && bface_i.element <= last_el
        elnum = bface_i.element
        facenum_local = bface_i.face
        nrm_bar = sview(nrm_peer_bar, :, :, i)

        fill!(Eone_el_bar, 0.0)
        assembleEone_rev(sbpface, elnum - offset, facenum_local, Eone_el_bar,
                     Eone_bar)
        calcEoneElement_rev(sbpface, nrm_bar, Rone, tmp, Eone_el_bar)
      end
    end
  end

  return nothing
end



"""
  Calculates E1 for a given interface.  Used by calcEone.

  Actually what it computes is:
      (R^T)*N*B*R*P*1

  Note that there should be a factor of P^T on the left.  That factor is
  applied in assembleEone().

  Inputs:
    sbpface: an AbstractFace
    nrm: the normal vectors for each node of the face, dim x numFaceNodes
    Rone: the interpolation operator R times the vector of ones
    tmp: a temporary vector of length numFaceNodes, overwritten

  Inputs/Outputs:
    Eone_el: sbp.stencilsize x dim matrix to be populated with the E
             contribution for this face.  Overwritten.

"""
function calcEoneElement(sbpface::AbstractFace, nrm::AbstractMatrix, 
                         Rone::AbstractVector, tmp::AbstractVector, 
                         Eone_el::AbstractMatrix)

  dim = size(Eone_el, 2)
  numFaceNodes = length(Rone)
  for d=1:dim
    for i=1:numFaceNodes
      tmp[i] = Rone[i]*nrm[d, i]*sbpface.wface[i]
    end

    Eone_dim = sview(Eone_el, :, d)
    smallmatvec!(sbpface.interp, tmp, Eone_dim)
  end

  return nothing
end

# back propigate Eone_el_bar to nrm_bar
"""
  This function uses reverse mode to back-propigate Eone_el_bar to nrm_bar

  Input:
    sbpface: an AbstractFace
    Rone: the interpolation operator R times a vector of ones
    tmp: a temporary vector, overwritten
    Eone_el_bar: the adjoint part of E1

  Inputs/Outputs:
    nrm_bar: the adjoint part of the scaled face normal vector in x-y space
             updated with the contribution from this function

"""
function calcEoneElement_rev(sbpface::AbstractFace,
                         nrm_bar::AbstractMatrix,
                         Rone::AbstractVector, tmp_bar::AbstractVector, 
                         Eone_el_bar::AbstractMatrix)

  dim = size(Eone_el_bar, 2)
  numFaceNodes = length(Rone)

  for d=1:dim
    # this function is linear, so no need for a forward sweep

    # reverse sweep
    Eone_dim_bar = sview(Eone_el_bar, :, d)
    # only back propigate to tmp_bar (sbpface.interp is unimportant)
    fill!(tmp_bar, 0.0)
    smallmatvec_revv!(sbpface.interp, tmp_bar, Eone_dim_bar)

    for i=1:numFaceNodes
      nrm_bar[d, i] += tmp_bar[i]*Rone[i]*sbpface.wface[i]
    end
  end

  return nothing
end





"""
  Takes an Eone_el matrix from calcEoneElement and assemble it into the big
  Eone.

  Inputs:
    sbpface: an SBP face
    elnum: the element index of Eone to put the values into.  Note that if
           Eone is for the entire mesh, then this is the element number.  If
           Eone is for a range of elements, then elnum is the index of the
           element within the range.
    Eone_el: the E1 contribution of an element

  Inputs/Output
    Eone: a numNodesPerElement x dim x blocksize array to store E1 for a
          range of elements, updated with the new contribution
"""
function assembleEone{Tmsh}(sbpface::AbstractFace, elnum::Integer, 
                      facenum_local::Integer, Eone_el::AbstractMatrix{Tmsh}, 
                      Eone::AbstractArray{Tmsh, 3})

  dim = size(Eone, 2)
  for d=1:dim
    for i=1:sbpface.stencilsize
      p_i = sbpface.perm[i, facenum_local]
      Eone[p_i, d, elnum] += Eone_el[i, d]
    end
  end

  return nothing
end

"""
  Back propigates Eone_bar to Eone_el_bar.

  Inputs:
    sbpface: an AbstractFace
    elnum: the element number
    facenum_local: the local number of the face that Eone_el is calculated for
    Eone_bar: the adjoint part

  Inputs/Outputs:
    Eone_el_bar: the adjoint part
"""
function assembleEone_rev{Tmsh}(sbpface::AbstractFace, elnum::Integer, 
                      facenum_local::Integer,
                      Eone_el_bar::AbstractMatrix{Tmsh},
                      Eone_bar::AbstractArray{Tmsh, 3})

  dim = size(Eone_bar, 2)
  for d=1:dim  # TODO: switch loops: turn an indexed store into a strided store
    for i=1:sbpface.stencilsize
      p_i = sbpface.perm[i, facenum_local]
      # reverse mode step
      Eone_el_bar[i, d] += Eone_bar[p_i, d, elnum]
    end
  end

  return nothing
end


# get the coordinates of the nodes on a given face (in the coordinate field, not
# the solution field)
# this only works up to second order, and uses the ElementTopology to figure
# out the orientation of the edge
# elnum is the global element number
# facenum is the local face number
"""
  This function gets the coordinates of the coordinate field nodes of the
  face of an element, in the orientation specified by mesh.topo.  This only
  works up to second order.  It uses the SBP element topology to determine
  the ordering of the nodes.  I guess this makes sense because calcFaceNormals
  is an SBP function and requires its definition of the vertices that define the
  face.  Because ref_vtx is passed into calcFaceNormals, it should be ok to 
  use the Pumi definition of the coordinates of the reference nodes (precisely
  because we use the SBP definition of which vertices define the face, the
  face node coordinates are relative to those vertices).
  
  Inputs:
    mesh: a 2D mesh
    elnum: the global element number
    facenum: the local face number

  Inputs/Outputs:
    coords: array to be populated with the coordinates, Tdim x 
            number of coordinate nodes on a face (2 for linear 2D, 
            3 for quadratic 2D)
"""
function getMeshFaceCoordinates(mesh::PumiMesh2DG, elnum::Integer,
                                facenum::Integer, 
                                coords::AbstractMatrix)

# coords = a 2 x number of nodes on the face array to be populated with the
# coordinates

  topo = mesh.topo

  # equivalent to topo.face_verts[:, facenum]
  v1_idx = topo.face_verts[1, facenum]
  v2_idx = topo.face_verts[2, facenum]
#  v1_idx = facenum
#  v2_idx = mod(facenum, 3) + 1

  coords[1, 1] = mesh.vert_coords[1, v1_idx, elnum]
  coords[2, 1] = mesh.vert_coords[2, v1_idx, elnum]

  coords[1, 2] = mesh.vert_coords[1, v2_idx, elnum]
  coords[2, 2] = mesh.vert_coords[2, v2_idx, elnum]

  if hasNodesIn(mesh.coordshape_ptr, 1)
#    edge_idx = facenum + 3  # for simplex elements
    edge_idx = 3 + topo.face_edges[1, facenum]
    coords[1, 3] = mesh.vert_coords[1, edge_idx, elnum]
    coords[2, 3] = mesh.vert_coords[2, edge_idx, elnum]
  end

  return nothing
end

"""
  This function stores the adjoint part of the vert_coords to a Pumi Field

  This stores to the adjoint part of the coordinate field, mesh.coords_bar_ptr.
  2nd order coordinate fields only.

  This function *accumulates* into mesh.coords_bar_ptr.  Users should generally
  call zeroField() on that pointer before using this function the first time.

  Inputs
    mesh::PumiMesh2DG
    elnum: the element number
    facenum: the local face number
    storage :FaceCoordinateStorage for some temporary arrays
    coords_bar:  the adjoint part of the coordinate field for the given face,
                 2 x coord_numNodesPerFace.
                 The data should be ordered vertices, then mid edge nodes.
"""
function getMeshFaceCoordinates_rev(mesh::PumiMesh2DG, elnum::Integer,
                                facenum::Integer, 
                                coords_bar::AbstractMatrix)

#  vertmap = mesh.topo_pumi.face_verts
  topo = mesh.topo_pumi

  # equivalent to topo.face_verts[:, facenum]
  v1_idx = topo.face_verts[1, facenum]
  v2_idx = topo.face_verts[2, facenum]

  # equivalent to topo.face_verts[:, facenum]
#  v1_idx = facenum
#  v2_idx = mod(facenum, 3) + 1

  mesh.vert_coords_bar[1, v1_idx, elnum] += coords_bar[1, 1]
  mesh.vert_coords_bar[2, v1_idx, elnum] += coords_bar[2, 1]

  mesh.vert_coords_bar[1, v2_idx, elnum] += coords_bar[1, 2]
  mesh.vert_coords_bar[2, v2_idx, elnum] += coords_bar[2, 2]

  if hasNodesIn(mesh.coordshape_ptr, 1)
    edge_idx = facenum + 3  # for simplex elements
    mesh.vert_coords_bar[1, edge_idx, elnum] += coords_bar[1, 3]
    mesh.vert_coords_bar[2, edge_idx, elnum] += coords_bar[2, 3]
  end


  return nothing
end


function getMeshFaceCoordinates(mesh::PumiMesh3DG, elnum::Integer,
                                facenum::Integer, coords::AbstractMatrix)

#  println("\ngetting coordinate for element ", elnum, " face ", facenum)
#  vertmap = mesh.topo.face_verts
  topo = mesh.topo
#  println("topo_pumi.face_verts = \n", topo.face_verts)
  # get the face vertices
  face_vert_idx = sview(topo.face_verts, :, facenum)
#  println("face ", facenum, " has verts \n", face_vert_idx)

  for i=1:3  # 3 vertices per face
    vert_i = face_vert_idx[i]
#    println("vertex ", vert_i)
    for j=1:3  # 3 coordinates at each vertex
      coords[j, i] = mesh.vert_coords[j, face_vert_idx[i], elnum]
    end
  end

#  println("coords = \n", coords)
#=
  # check against pumi
  el_i = mesh.elements[elnum]
  down_verts = Array(Ptr{Void}, 12)
  getDownward(mesh.m_ptr, el_i, 0, down_verts)

  # get all vertex coordinates
  coords_pumi = zeros(3, 4)
  for i=1:4
    getPoint(mesh.m_ptr, down_verts[i], 0, sview(coords_pumi, :, i))
  end


  println("all vertex coords (array) = \n", mesh.vert_coords[:, :, elnum])
  println("all vertex coords (pumi) = \n", coords_pumi)

  coords_tmp = zeros(3, 3)
  failflag = false
  vertmap = mesh.topo.face_verts
  println("topo_sbp.face_verts = \n", vertmap)
  for i=1:3  # 3 vertices per face
    v_i = down_verts[ vertmap[i, facenum] ]
    getPoint(mesh.m_ptr, v_i, 0, sview(coords_tmp, :, i))
#    println("vertex ", vertmap[i, facenum], " has coords = \n", sview(coords_tmp, :, i))
    for j=1:3
      if abs(coords_tmp[j, i] - coords[j, i]) > 1e-13
        #=
        println("setting failflag")
        println("i = ", i)
        println("coords_tmp[:, i] = \n", coords_tmp[:, i])
        println("coords[:, i] = \n", coords[:, i])
        =#
        failflag = failflag || true
      end
    end
  end

  println("array vertices ", face_vert_idx, " have coordinates = \n", coords)
  println("pumi vertices ", vertmap[:, facenum], " have coordinates = \n", coords_tmp)
#  @assert !failflag
  =#
  if hasNodesIn(mesh.coordshape_ptr, 1)
#    println("getting 2nd order node")

    for i=1:3  # 3 edges per face
      edge_idx = topo.face_edges[i, facenum] + 4 # 4 = numVerts
      for j=1:3  # 3 coordinates per face
        coords[j, 3 + i] = mesh.vert_coords[j, edge_idx]
      end
    end

  end

  return nothing
end

"""
  Similar to the 2d method.

  coords_bar should be ordered vertices, edge nodes, then face nodes
"""
function getMeshFaceCoordinates_rev(mesh::PumiMesh3DG, elnum::Integer,
                                facenum::Integer, 
                                coords_bar::AbstractMatrix)

#  vertmap = mesh.topo.face_verts
  topo = mesh.topo
  
  # get the face vertices
  face_vert_idx = sview(topo.face_verts, :, facenum)

  for i=1:3  # 3 vertices per face
    for j=1:3  # 3 coordinates at each vertex
      mesh.vert_coords_bar[j, face_vert_idx[i]] += coords_bar[j, i]
    end
  end

   if hasNodesIn(mesh.coordshape_ptr, 1)
    offset = 3 + (facenum - 1)*3  # 3 vertices + 1 edge node per face

    for i=1:3  # 3 edges per face
      edge_idx = topo.face_edges[i, facenum] + 4  # 4 = numVerts
      for j=1:3  # 3 coordinates per face
        mesh.vert_coords_bar[j, edge_idx] += coords_bar[j, 3 + i]  # 3 verts
      end
    end

  end


  return nothing
end
