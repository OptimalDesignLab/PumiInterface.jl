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
  (not the node coordinates) and puts them the mesh.vert_coords field
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

  calcFaceCoordinatesAndNormals(mesh, sbp, mesh.bndryfaces, mesh.coords_bndry, 
                               mesh.nrm_bndry)
  calcFaceCoordinatesAndNormals(mesh, sbp, mesh.interfaces, 
                               mesh.coords_interface, mesh.nrm_face)
  for i=1:mesh.npeers
    calcFaceCoordinatesAndNormals(mesh, sbp, mesh.bndries_local[i], 
                               mesh.coords_sharedface[i], mesh.nrm_sharedface[i])
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
  storage = FaceCoordinateStorage()

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
      getMeshFaceCoordinates(mesh, el_i, face_i, storage, coords_i)
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
    getMeshFaceCoordinates(mesh, el_i, face_i, storage, coords_i)
    face_idx += 1
  end

  calcFaceNormals!(mesh.sbpface, mesh.coord_order, ref_verts, 
                   coords_lag_face_block, coords_face_block, nrm_face_block)


  # make sure the normal vectors point outwards
  fixOutwardNormal(mesh, faces, nrm_face)

  return nothing
end

"""
  This function check to make sure each face normal vector is oriented
  outwards, and flips it if needed.
"""
function fixOutwardNormal{I <: Union{Boundary, Interface}, Tmsh}(mesh, 
                          faces::AbstractArray{I, 1},
                          nrm_face::AbstractArray{Tmsh, 3})


  tmp = zeros(3)  # temporary vector to hold coordinates
  topo = mesh.topo
  numVertPerElement = mesh.numTypePerElement[1]
  numVertPerFace = numVertPerElement - 1

  # temporary arrays
  el_verts = Array(Ptr{Void}, numVertPerElement)
  other_vert_coords = zeros(mesh.dim)
  face_verts = Array(Ptr{Void}, numVertPerElement - 1)
  face_vert_coords = zeros(mesh.dim, numVertPerFace)

  nfaces = length(faces)
  for i=1:nfaces
    iface_i = faces[i]
    elnum = getElementL(iface_i)
    facenum_local = getFaceL(iface_i)

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
    for j=1:mesh.numNodesPerFace
      outward_count = 0  # count number of calculations that showed outward
      for k=1:numVertPerFace
        val = zero(Float64)
        for p=1:mesh.dim
          r1_p = other_vert_coords[p] - face_vert_coords[p, k]
          val += nrm_face[p, j, i]*r1_p  # accumulate dot product
        end
        
        if val < 0
          outward_count += 1
        end  # end if
      end  # end k

      # reverse face if needed, throw exception is unclear
      if outward_count == 0  # no calculation found the normal is outward
        for p=1:mesh.dim
          nrm_face[p, j, i] = -nrm_face[p, j, i]
        end
      # some, but not all, calculations found outward
      elseif outward_count < numVertPerFace          
        throw(ErrorException("face $iface_i, node $j has indeterminate orientation"))
      end  # else the face is oriented outwards, do nothing

    end  # end j

  end  # end loop i

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
  if mesh.dim == 2
    calcMappingJacobian!(sbp, mesh.coord_order, ref_vtx, mesh.vert_coords, 
                         mesh.coords, mesh.dxidx, mesh.jac)
  else  # need to calculate Eone
    println("doing 3d curvilinear coordinate and metric calcullation")

    # block format
    blocksize = 1000  # number of elements per block
    nblocks_full = div(mesh.numEl, blocksize)
    nrem = mesh.numEl % blocksize

    Eone = zeros(mesh.numNodesPerElement, mesh.dim, blocksize)

    for block=1:nblocks_full
      println("full block ", block)
      start_idx = (block - 1)*blocksize + 1
      end_idx = block*blocksize
      element_range = start_idx:end_idx
      
      getCurvilinearMetricsAndCoordinates_inner(mesh, sbp, element_range, Eone)
     
    end  # end loop over blocks

    # remainder block
    println("doing remainder block")
    start_idx = nblocks_full*blocksize + 1
    end_idx = mesh.numEl
    @assert end_idx - start_idx + 1 <= blocksize
    @assert end_idx - start_idx + 1 == nrem

    element_range = start_idx:end_idx

    getCurvilinearMetricsAndCoordinates_inner(mesh, sbp, element_range, Eone)
  end  # end if dim == 2

  return nothing
end

global call_cnt = 1
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

  # for the remainder loop, make sure the last dimensions of Eone matches the
  # other arrays
  Eone_block = sview(Eone, :, :, element_range)

  global call_cnt
  println("call_cnt = ", call_cnt)
  calcEone(mesh, sbp, element_range, Eone_block)
  
  writedlm("Eone_$call_cnt.dat", Eone_block)
  call_cnt += 1
  calcMappingJacobian!(sbp, mesh.coord_order, ref_vtx, vert_coords_block, 
                       coords_block, dxidx_block, jac_block, Eone_block)

  fill!(Eone, 0.0)
  return nothing
end

function calcEone{Tmsh}(mesh::PumiMeshDG{Tmsh}, sbp, element_range, 
                        Eone::AbstractArray{Tmsh, 3})

  # search interfaces and bndryfaces

  first_el = element_range[1]
  last_el = element_range[end]

  # R times vector of ones
  println("R = \n", mesh.sbpface.interp.')
  Rone = vec(sum(mesh.sbpface.interp.', 2))
  println("Rone = \n", Rone)
  sbpface = mesh.sbpface
  tmp = zeros(Rone)
  nrmL = zeros(Tmsh, mesh.dim, sbpface.numnodes)

  # accumulate E1 for a given element
  Eone_el = zeros(Tmsh, sbpface.stencilsize, mesh.dim)

  for i=1:mesh.numInterfaces
    println("interface ", i)
    iface_i = mesh.interfaces[i]

    if iface_i.elementL >= first_el && iface_i.elementL <= last_el
      println("first element in range")
      elnum = iface_i.elementL
      facenum_local = iface_i.faceL
      nrm = sview(mesh.nrm_face, :, :, i)

      # call inner functions
      calcEoneElement(sbpface, nrm, Rone, tmp, Eone_el)
      println("iface = ", iface_i)
      println("nrm = \n", nrm)
      println("Eone_el = \n", Eone_el)
      assembleEone(sbpface, elnum, facenum_local, Eone_el, Eone)
      println("Eone[:, :, $elnum] = \n", Eone[:, :, elnum])

    end

    if iface_i.elementR >= first_el && iface_i.elementR <= last_el
      println("second element in ranage")
      # do same as above, negating and permuting nrm
      elnum = iface_i.elementR
      facenum_local = iface_i.faceR
      orient = iface_i.orient
      for j=1:sbpface.numnodes
        for d=1:mesh.dim
          #TODO: is this the right use of nbrperm?  or should it the
          # receiving array be permuted
          nrmL[d, j] = -mesh.nrm_face[d, sbpface.nbrperm[j, orient], i]
        end
      end
 
      calcEoneElement(sbpface, nrmL, Rone, tmp, Eone_el)
      println("iface = ", iface_i)
      println("nrm = \n", nrm)
      println("Eone_el = \n", Eone_el)

      assembleEone(sbpface, elnum, facenum_local, Eone_el, Eone)
      println("Eone[:, :, $elnum] = \n", Eone[:, :, elnum])
    end  # end if/else

  end  # end loop over interfaces

  for i=1:mesh.numBoundaryFaces
    println("boundary ", i)
    bface_i = mesh.bndryfaces[i]

    if bface_i.element >= first_el && bface_i.element <= last_el
      println("boundary face in range")
      elnum = bface_i.element
      facenum_local = bface_i.face
      nrm = sview(mesh.nrm_bndry, :, :, i)

      # call inner functions
      calcEoneElement(sbpface, nrm, Rone, tmp, Eone_el)
      println("iface = ", bface_i)
      println("nrm = \n", nrm)
      println("Eone_el = \n", Eone_el)

      assembleEone(sbpface, elnum, facenum_local, Eone_el, Eone)

      println("Eone[:, :, $elnum] = \n", Eone[:, :, elnum])
    end
  end  # end loop over boundary faces

  return nothing
end

"""
  Calculates E1 for a given interface.  Used by calcEone.

  Inputs:
    sbpface: an AbstractFace
    nrm: the normal vectors for each node of the face, dim x numFaceNodes
    Rone: the interpolation operator R times the vector of ones
    tmp: a temporary vector of length numFaceNodes, overwritten

  Inputs/Outputs:
    Eone_el: sbp.stencilsize x dim matrix to be populated with the E
             contribution for this face.  Overwritten.

"""
function calcEoneElement{Tmsh}(sbpface::AbstractFace, nrm::AbstractMatrix, 
                         Rone::AbstractVector, tmp::AbstractVector, 
                         Eone_el::AbstractMatrix{Tmsh})

  dim = size(Eone_el, 2)
  println("wface = \n", sbpface.wface)
  numFaceNodes = length(Rone)
  for d=1:dim
    println("d = ", d)
    for i=1:numFaceNodes
      tmp[i] = Rone[i]*nrm[d, i]*sbpface.wface[i]
    end

    println("tmp = \n", tmp)

    Eone_dim = sview(Eone_el, :, d)
    smallmatvec!(sbpface.interp, tmp, Eone_dim)
    println("Eone_dim = \n", Eone_dim)
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

  println("assembling into element ", elnum)
#  println("Eone_el = \n", Eone_el)
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
  Type to store the temporary arrays needed by getMeshFaceCoordinates
"""
immutable FaceCoordinateStorage
  el_verts::Array{Ptr{Void}, 1}
  coords_tmp::Array{Float64, 1}
  v1_edges::Array{Ptr{Void}, 1}
  v2_edges::Array{Ptr{Void}, 1}

  function FaceCoordinateStorage()
    el_verts = Array(Ptr{Void}, 12)
    coords_tmp = Array(Float64, 3)
    v1_edges = Array(Ptr{Void}, 400)
    v2_edges = Array(Ptr{Void}, 400)

    return new(el_verts, coords_tmp, v1_edges, v2_edges)
  end
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
  works up to second order.  In 3D it makes the assumption that edge 1 of
  the triangle is defined from v1 -> v2, edge 2 is v2 -> v3, and edge 3 is
  v3 -> v1.

  Inputs:
    mesh: a 2D mesh
    elnum: the global element number
    facenum: the local face number

  Inputs/Outputs:
    coords: array to be populated with the coordinates, Tdim x 
            number of coordinate nodes on a face (2 for linear 2D, 
            3 for quadratic 2D)
"""
function getMeshFaceCoordinates(mesh::PumiMesh2DG, elnum, facenum, storage::FaceCoordinateStorage, coords::AbstractMatrix)

# coords = a 2 x number of nodes on the face array to be populated with the
# coordinates

  vertmap = mesh.topo.face_verts
  el = mesh.elements[elnum]
  el_verts = storage.el_verts
  coords_tmp = storage.coords_tmp
#  el_verts = Array(Ptr{Void}, 12)
#  coords_tmp = Array(Float64, 3)
  
  getDownward(mesh.m_ptr, el, 0, el_verts)

  v1 = el_verts[vertmap[1, facenum]]
  v2 = el_verts[vertmap[2, facenum]]

  getPoint(mesh.m_ptr, v1, 0, coords_tmp)

  coords[1, 1] = coords_tmp[1]
  coords[2, 1] = coords_tmp[2]

  getPoint(mesh.m_ptr, v2, 0, coords_tmp)
  coords[1, 2] = coords_tmp[1]
  coords[2, 2] = coords_tmp[2]

  if hasNodesIn(mesh.coordshape_ptr, 1)
    getDownward(mesh.m_ptr, el, 1, el_verts)
    edge = el_verts[facenum]
    getPoint(mesh.m_ptr, edge, 0, coords_tmp)

    coords[1, 3] = coords_tmp[1]
    coords[2, 3] = coords_tmp[2]
  end

  return nothing
end


function getMeshFaceCoordinates(mesh::PumiMesh3DG, elnum, facenum, storage::FaceCoordinateStorage, coords::AbstractMatrix)

  vertmap = mesh.topo.face_verts
  el = mesh.elements[elnum]
  el_verts = storage.el_verts
#  el_verts = Array(Ptr{Void}, 12)
  
  getDownward(mesh.m_ptr, el, 0, el_verts)

  for i=1:3  # 3 vertices per face
    v_i = el_verts[ vertmap[i, facenum] ]
    coords_i = sview(coords, :, i)
    getPoint(mesh.m_ptr, v_i, 0, coords_i)
  end

  if hasNodesIn(mesh.coordshape_ptr, 1)
    offset = 3
    # edge nodes
    v1_edges = storage.v1_edges
    v2_edges = storage.v2_edges
#    v1_edges = Array(Ptr{Void}, 400)  # apf::Up
#    v2_edges = Array(Ptr{Void}, 400)

    for i=1:3  # 3 edges per face
      i2 = mod(i, 3) + 1
      v1 = el_verts[i]
      v2 = el_verts[i2]
      countAdjacent(mesh.m_ptr, v1, 1)
      getAdjacent(v1_edges)
      countAdjacent(mesh.m_ptr, v2, 1)
      getAdjacent(v2_edges)

      # simple linear search for common edge
      # in the average case this is fine because there are only 20 or so
      # edges upward ajacent to a vertex, but in the worst case there could
      # be 400
      common_edge = first_common(v1_edges, v2_edges)

      coords_i = sview(coords, :, i + offset)
      getPoint(mesh.m_ptr, common_edge, 0, coords_i)
    end
  end

  return nothing
end


    








