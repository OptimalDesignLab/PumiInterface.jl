# functiosn related to data about faces (edges in 2D)

function getBoundaryArray(mesh::PumiMesh, boundary_nums::AbstractArray{Int, 2})
# get an array of type Boundary for SBP
# creating an an array of a user defined type seems like a waste of memory operations
# bnd_array is a vector of type Boundary with length equal to the number of edges on the boundary of the mesh

#  mesh.bndryfaces = Array(Boundary, mesh.numBoundaryFaces)

  for i=1:mesh.numBoundaryFaces
    facenum = boundary_nums[i, 1]
    edgenum_global = boundary_nums[i, 2]
    edgenum_local = getBoundaryFaceLocalNum(mesh, edgenum_global)
    mesh.bndryfaces[i] = Boundary(facenum, edgenum_local)
  end

  return nothing
end

"""
  This function allocates the arrays of normal vector for interfaces, boundary
  faces, and shared faces.

  For interfaces, the normal is calculated for elementL
"""
function allocateNormals{Tmsh}(mesh::PumiMeshDG{Tmsh}, sbp)

  dim = mesh.dim
  numfacenodes = sbp.numfacenodes
  mesh.nrm_bndry = Array(Tmsh, dim, numfacenodes, mesh.numBoundaryFaces )


  mesh.nrm_face = Array(Tmsh, mesh.dim, numfacenodes, mesh.numInterfaces)
  mesh.nrm_sharedface = Array(Array{Tmsh, 3}, mesh.npeers)
  for i=1:mesh.npeers
    mesh.nrm_sharedface[i] = Array(Tmsh, mesh.dim, numfacenodes.mesh.peer_face_counts[i])
  end

  return nothing
end


"""
  This function calculates the face normal vectors for all interfaces, 
  boundary faces, and shared faces.  A smart allocator is used to allocate
  the arrays if needed.

  The interpolated dxidx arrays must be populated before this function is
  called.

  Inputs:
    mesh: a PumiMeshdG
    sbp: an SBP operator
"""
function getFaceNormals(mesh::PumiMeshDG, sbp)

  if !isFieldDefined(mesh, :nrm_bndry, :nrm_face, :nrm_sharedface)
    allocateNormals(mesh, sbp)
  end

  calcFaceNormal(mesh, sbp, mesh.bndryfaces, mesh.dxidx_bndry, mesh.nrm_bndry)
  calcFaceNormal(mesh, sbp, mesh.interfaces, mesh.dxidx_face, mesh.nrm_face)
  for i=1:mesh.npeers
    calcFaceNormal(mesh, sbp, mesh.bndries_local[i], mesh.dxidx_sharedface[i],
                   mesh.nrm_sharedface[i])
  end

  return nothing
end



"""
  This function calculates the face normal for a set of faces given the
  dxidx values for each face node and the facenormal in parametric space

  Inputs:
    mesh: a DG mesh
    sbp: an SBP operator containing the face normals in parametric space
    faces: list of Interfaces or Boundaries, should be dim x numFacesPerElement
    dxidx: d(xi)/d(x) at the face nodes, should be dim x dim x numfaceNodes x 
           length(faces)

  Inputs/Outputs:
    nrm: array, dim x numfacenodes x length(faces) containing the face normal
         vector in x-y space at each face node of each face

  Aliasing restrictions: none
"""
function calcFaceNormal{Tmsh, Tbndry <: Union{Boundary, Interface}}(mesh::PumiMeshDG, sbp, faces::AbstractArray{Tbndry, 1}, dxidx::AbstractArray{Tmsh, 4}, nrm::AbstractArray{Tmsh, 3})

  nfaces = length(faces)
  dim = mesh.dim
  for i=1:nfaces
    faceL = getFaceL(faces[i])
#    faceL = faces[i].faceL

    for j=1:sbp.numfacenodes
      nrm_j = sview(nrm, :, j, i)
      for dim=1:dim
        nrm_d = zero(Tmsh)
        for d=1:dim
          nrm_d += sbpface.normal[d, faceL]*dxidx[d, dim, j, i]
        end  # end loop d
        nrm_j[dim] = nrm_d
      end  # end loop dim
    end  # end loop j

  end  # end loop i

  return nothing
end

# deprecated
function getBoundaryFaceNormals{Tmsh}(mesh::PumiMesh2D, sbp::AbstractSBP, bndry_faces::AbstractArray{Boundary, 1}, face_normals::Array{Tmsh, 3})

  nfaces = length(bndry_faces)

  alpha = zeros(Tmsh, 2,2)
  dxidx = mesh.dxidx
  jac = mesh.jac
  for i=1:nfaces
    element_i = bndry_faces[i].element
    face_i = bndry_faces[i].face
    for j=1:sbp.numfacenodes
      node_index = sbp.facenodes[j, face_i]

      # calculate the mysterious parameter alpha
      for di1=1:2
	for di2=1:2
	  alpha[di1, di2] = (dxidx[di1,1, node_index, element_i].*dxidx[di2,1, node_index, element_i] + dxidx[di1,2, node_index, element_i].*dxidx[di2,2, node_index, element_i])*jac[node_index,element_i]
	end
      end

    # call SBP function
    smallmatvec!(alpha, sview(sbp.facenormal, :, face_i), sview(face_normals, :, j, i))

    end  # end loop over face nodes
  end  # end loop over faces

  return nothing
end

# deprecated
function getInternalFaceNormals{Tmsh}(mesh::PumiMesh2D, sbp::AbstractSBP, internal_faces::AbstractArray{Interface, 1}, face_normals::Array{Tmsh, 4})

  nfaces = length(internal_faces)

  alpha = zeros(Tmsh, 2,2)
  dxidx = mesh.dxidx
  jac = mesh.jac
  for i=1:nfaces
    element_iL = internal_faces[i].elementL
    face_iL = internal_faces[i].faceL
    element_iR = internal_faces[i].elementR
    face_iR = internal_faces[i].faceR
    for j=1:sbp.numfacenodes

      # calculate left face normal
      node_index = sbp.facenodes[j, face_iL]

      # calculate the mysterious parameter alpha
      for di1=1:2
	for di2=1:2
	  alpha[di1, di2] = (dxidx[di1,1, node_index, element_iL].*dxidx[di2,1, node_index, element_iL] + dxidx[di1,2, node_index, element_iL].*dxidx[di2,2, node_index, element_iL])*jac[node_index,element_iL]
	end
      end

    # call SBP function
    smallmatvec!(alpha, sview(sbp.facenormal, :, face_iL), sview(face_normals, :, 1, j, i))
      # calculate right fae normal
      node_index = sbp.facenodes[j, face_iR]

      # calculate the mysterious parameter alpha
      for di1=1:2
	for di2=1:2
	  alpha[di1, di2] = (dxidx[di1,1, node_index, element_iR].*dxidx[di2,1, node_index, element_iR] + dxidx[di1,2, node_index, element_iR].*dxidx[di2,2, node_index, element_iR])*jac[node_index,element_iR]
	end
      end

    # call SBP function
    smallmatvec!(alpha, sview(sbp.facenormal, :, face_iR), sview(face_normals, :, 2, j, i))




    end  # end loop over face nodes
  end  # end loop over faces

  return nothing
end


function getInterfaceArray(mesh::PumiMesh2D)
# get array of [elementL, elementR, edgeL, edgeR] for each internal edge,
# where elementL and R are the elements that use the edge, edgeL R are the
# local number of the edge within the element
# interfaces is the array to be populated with the data, it must be 
# number of internal edges by 1
# the nodemap is used to determine the node ordering of elementR relative to
# elementL (ie. the relative orientation is determined between the nodes of 
# the 2 element in SBP order)
# the edge is "owned" by elementL, ie the first node on the edge is the first
# node from the perspective of elementL

  # only need internal boundaries (not external)
  # get number of nodes affecting an edge
  num_edge_nodes = countAllNodes(mesh.mshape_ptr, 1)

  # unused variable?
  nodemap = Array(num_edge_nodes:(-1):1)
  part_nums = Array(Cint, 1)
  matched_entities = Array(Ptr{Void}, 1)
  seen_entities = Set{Ptr{Void}}()
  sizehint!(seen_entities, mesh.numPeriodicInterfaces)

  pos = 1 # current position in interfaces
  for i=1:mesh.numEdge
    edge_i = mesh.edges[i]
    if in(edge_i, seen_entities)
      continue
    end
    # get number of elements using the edge
    adjacent_nums, num_adjacent = getAdjacentEntityNums(mesh, i, 1, 2)
    n = countMatches(mesh.m_ptr, mesh.edges[i])
    getMatches(part_nums, matched_entities)
    has_local_match = part_nums[1] == mesh.myrank && n > 0

    if num_adjacent > 1 || has_local_match # internal edge
      if !has_local_match
        element1 = adjacent_nums[1]
        element2 = adjacent_nums[2]
        edge1 = i
        edge2 = i
      else  # does have local match
        edge2_ptr = matched_entities[1]
        push!(seen_entities, edge2_ptr)
        element1 = adjacent_nums[1]
        num_adjacent = countAdjacent(mesh.m_ptr, edge2_ptr, mesh.dim)
        adjacent_entities = getAdjacent(num_adjacent)
        element2 = getNumberJ(mesh.el_Nptr, adjacent_entities[1], 0, 0) + 1
        edge1 = i
        edge2 = getNumberJ(mesh.edge_Nptr, edge2_ptr, 0, 0) + 1
      end

      coords_1 = zeros(3,3)
      coords_2 = zeros(3,3)
      getFaceCoords(mesh.elements[element1], coords_1, 3, 3)
      getFaceCoords(mesh.elements[element2], coords_2, 3, 3)

      # calculate centroid
      centroid1 = sum(coords_1, 2)
      centroid2 = sum(coords_2, 2)

      if abs(centroid1[1] - centroid2[1]) > 1e-10  # if big enough difference
        if centroid1[1] < centroid2[1]
    elementL = element1
    elementR = element2
    edgeL = edge1
    edgeR = edge2
        else
    elementL = element2
    elementR = element1
    edgeL = edge2
    edgeR = edge1
        end
      else  # use y coordinate to decide
        if centroid1[2] < centroid2[2]
    elementL = element1
    elementR = element2
    edgeL = edge1
    edgeR = edge2
        else
    elementL = element2
    elementR = element1
    edgeL = edge2
    edgeR = edge1
        end
      end

      edgeL = getFaceLocalNum(mesh, edgeL, elementL)
      edgeR = getFaceLocalNum(mesh, edgeR, elementR)

      # here we set the orientation flag to 1 for all edges, because in 2D
      # the edges always have opposite orientation
      mesh.interfaces[pos] = Interface(elementL, elementR, edgeL, edgeR, UInt8(1))
      pos += 1

  #    print("\n")

    end  # end if internal edge

#     print("\n")
  end  # end loop over edges

  return nothing

end  # end function


function getInterfaceArray(mesh::PumiMesh3D)
  adj_elements = Array(Ptr{Void}, 2)
  coords1 = zeros(3, 4)
  coords2 = zeros(3, 4)
  centroid1 = zeros(3)
  centroid2 = zeros(3)
  vertsL = Array(Ptr{Void}, 4)
  vertsR = Array(Ptr{Void}, 4)
  facevertsL = Array(Ptr{Void}, 3)
  facevertsR = Array(Ptr{Void}, 3)
  # there can be up to 400 elements using a vertex + 8 
  # corners of the domain, so in the worst case a parallel partition could
  # generate 400 + 8 matches
  part_nums = Array(Cint, 400 + 8)
  matched_entities = Array(Ptr{Void}, 400 + 8)  # there can be up to 400 
                                           
  pos = 1 # position in mesh.interfaces
  seen_entities = Set{Ptr{Void}}()
  sizehint!(seen_entities, mesh.numPeriodicInterfaces)
 
  for i=1:mesh.numFace
    face_i = mesh.faces[i]
    if face_i in seen_entities  # avoid counting matched entities twice
      continue
    end

    num_adjacent = countAdjacent(mesh.m_ptr, face_i, mesh.dim)

    n = countMatches(mesh.m_ptr, face_i)
    getMatches(part_nums, matched_entities)
    has_local_match = part_nums[1] == mesh.myrank && n > 0

    if num_adjacent == 2 || has_local_match  # this is a shared interface

      getAdjacent(adj_elements)
      if has_local_match
        # get the parent element
        el1 = adj_elements[1]
        edgenum1 = i
        elnum1 = getNumberJ(mesh.el_Nptr, el1 , 0, 0) + 1

        # get the matched face
        other_face = matched_entities[1]
        # get the parent element of the matched face
        countAdjacent(mesh.m_ptr, other_face, mesh.dim)
        getAdjacent(adj_elements)
        el2 = adj_elements[1]
        edgenum2 = getNumberJ(mesh.face_Nptr, other_face, 0, 0) + 1
        elnum2 = getNumberJ(mesh.el_Nptr, el2, 0, 0) + 1
        push!(seen_entities, other_face)
      else
        # get both parent elements
        el1 = adj_elements[1]
        el2 = adj_elements[2]
        edgenum1 = i
        edgenum2 = i
        elnum1 = getNumberJ(mesh.el_Nptr, el1, 0, 0) + 1
        elnum2 = getNumberJ(mesh.el_Nptr, el2, 0, 0) + 1
      end

      # decide which one is elementL
      getElementCoords(mesh, el1, coords1)
      getElementCoords(mesh, el2, coords2)

      for j=1:4
        for k=1:3
          centroid1[k] += coords1[k, j]
          centroid2[k] += coords2[k, j]
        end
      end

      flag = getLR(centroid1, centroid2)
      if flag
        coordsL = coords2
        coordsR = coords1
        elnumL = elnum2
        elnumR = elnum1
        edgenumL = edgenum2
        edgenumR = edgenum1
        elL = el2
        elR = el1
      else
        coordsL = coords1
        coordsR = coords2
        elnumL = elnum1
        elnumR = elnum2
        edgenumL = edgenum1
        edgenumR = edgenum2
        elL = el1
        elR = el2
      end

      localfacenumL = getFaceLocalNum(mesh, edgenumL, elnumL)
      localfacenumR = getFaceLocalNum(mesh, edgenumR, elnumR)
      fdata = FaceData(elnumL, elL, elnumR, elR, localfacenumL, localfacenumR,
                       vertsL, vertsR, facevertsL, facevertsR, part_nums, 
                       matched_entities, has_local_match )

      rel_rotate = getRelativeOrientation(fdata, mesh)
      mesh.interfaces[pos] = Interface(elnumL, elnumR, localfacenumL, localfacenumR, UInt8(rel_rotate))
      pos += 1

      fill!(centroid1, 0.0)
      fill!(centroid2, 0.0)

    end  # end if 
  end

  return nothing

end

function getLR(coordsL, coordsR)
# return true if coords should be reversed
  tol = 1e-10

  if abs(coordsL[1] - coordsR[1]) > tol
    if coordsL[1] > coordsR[1]
      return false
    else
      return true
    end
  elseif abs(coordsL[2] - coordsR[2]) > tol
    if coordsL[2] > coordsR[2]
      return false
    else
      return true
    end
  else  # if the first two coordinates are not decisive, use the last one
    if coordsL[3] > coordsR[3]
      return false
    else
      return true
    end
  end

  return false  # this should never be reached
end
