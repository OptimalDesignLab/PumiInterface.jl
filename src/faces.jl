# functiosn related to data about faces (edges in 2D)


"""
  This function gets the information about boundary and interior faces
  and stores it to the mesh object.  In particular, the following fields of
  the mesh are allocated and populatedby this function:
    numBoundaryFaces
    numInterfaces
    numPeriodicInterfaces
    numBC
    bndry_geo_nums
    bndryfaces
    interfaces

  The following fields are allocated but not populated:
    bndry_funcs
    bndry_funcs_revm

  Because the curvilinear metric calculation requires knowing about *all*
  the mesh faces (to calculate Eone), this function creates a new boundary
  condition and assigns all faces not used by one of the existing BCs to it.
  The BC_name for this new BC is defaultBC, which allows the solver to do the
  face integral for faces that don't have a regular BC.  The new boundary
  condition is only created if there are faces not previously assigned to a BC.

  Inputs/Outputs:
    mesh: a mesh object
    opts: the options dictionary, which contains BC info and is updated with the
           new BC info if needed

  Outputs:
    bndry_nums: array of [element number, global face number] for each
                face on a geometric face.  This is returned to make testing
                easier.
"""
function getAllFaceData(mesh::PumiMesh, opts)

  # TODO: change countBoundaryEdges to only take 1 argument
  nbndryfaces, ninterfaces, npinterfaces, geo_edge_nums = countBoundaryEdges(mesh)
  mesh.numBoundaryFaces = nbndryfaces
  mesh.numInterfaces = ninterfaces
  mesh.numPeriodicInterfaces = npinterfaces

  if typeof(mesh) <: PumiMeshCG
    @assert mesh.numPeriodicInterfaces == 0
  end

  unused_geo_edge_nums = popBCEdges(geo_edge_nums, opts)
  # create an additional BC if needed to make sure the face integral gets
  # done for the unused edges
  add_bc = length(unused_geo_edge_nums) != 0
  if add_bc
    numBC = opts["numBC"] + 1

    opts["numBC"] = numBC
    opts[string("BC", numBC)] = unused_geo_edge_nums
    opts[string("BC", numBC, "_name")] = "defaultBC"
  end

  # populate mesh.bndry_faces from options dictionary
  mesh.numBC = opts["numBC"]
  allocateBoundaryAndInterfaceArrays(mesh, opts)
  boundary_nums = Array(Int, mesh.numBoundaryFaces, 2)

  # get the element number, global face number of the mesh edges that have BCs
  offset = 1
  for i=1:mesh.numBC
    key_i = string("BC", i)
    model_edges = opts[key_i]
    # record geometric edges
    ngeo = length(model_edges)
    mesh.bndry_geo_nums[i] = Array(Int, ngeo)
    copy!(mesh.bndry_geo_nums[i], model_edges)
    mesh.bndry_offsets[i] = offset
    offset, print_warning = getMeshEdgesFromModel(mesh, model_edges, offset, boundary_nums)  # get the mesh edges on the model edge
    if !(add_bc && i == mesh.numBC)  # don't error for the added BC
      if print_warning
        throw(ErrorException("Cannot apply boundary conditions to periodic boundary, model entity $model_edges"))
      end
    end
    # offset is incremented by getMeshEdgesFromModel
  end

  mesh.bndry_offsets[mesh.numBC + 1] = offset # = num boundary edges

  # get boundary information for entire mesh
  getBoundaryArray(mesh, boundary_nums)

  # get interface info for the entire mesh
  getInterfaceArray(mesh)
  sort!(mesh.interfaces)  #TODO: see if this is actually a good idea

  return boundary_nums
end

"""
  Smart allocator for the following fields of the mesh:

    bndry_offsets
    bndry_funcs
    bndry_funcs_revm
    bndry_geo_nums
    bndryfaces
    interfaces

  Because there is no reaonsable default value for bndry_funcs and
  bndry_funcs_revm, these fields are not modified if this function is called
  repeatedly.
"""
function allocateBoundaryAndInterfaceArrays(mesh, opts)

  if !isFieldDefined(mesh, :bndry_offsets, :bndry_funcs, :bndry_funcs_revm,
                     :bndry_geo_nums, :bndryfaces, :interfaces)

    mesh.bndry_offsets = Array(Int, mesh.numBC + 1)
    mesh.bndry_funcs = Array(BCType, mesh.numBC)
    mesh.bndry_funcs_revm = Array(BCType_revm, mesh.numBC)
    mesh.bndry_geo_nums = Array(Array{Int, 1}, mesh.numBC)

    mesh.bndryfaces = Array(Boundary, mesh.numBoundaryFaces)
    mesh.interfaces = Array(Interface, mesh.numInterfaces)
  else
    fill!(mesh.bndry_offsets, 0)
    fill!(mesh.bndry_geo_nums, 0)
    fill!(mesh.bndryfaces, Boundary(0, 0))
    fill!(mesh.interfaces, Interface(0, 0, 0, 0, 0))
  end

  return nothing
end


"""
  This function returns an array containing only the geometric edge numbers
  that do not have boundary conditions applied to them

  Inputs:
    geo_edge_nums: vector containing the numbers of all the geometric edges
    opts: the options dictionary

  Outputs:
    unused_geo_edge_nums: numbers of the geometric edges that do not have
                          boundary conditions applied to them

"""
function popBCEdges{I <: Integer}(geo_edge_nums::AbstractArray{I, 1}, opts)

  numBC = opts["numBC"]

  # create array of all model edges that have a boundary condition
  bndry_edges_BC = Array(Int, 0)
  for i=1:numBC
    key_i = string("BC", i)
    bndry_edges_i = opts[key_i]
    for j=1:length(bndry_edges_i)
      push!(bndry_edges_BC, bndry_edges_i[j])
    end
  end

  # find the (first) duplicate BC, throw error
  for i=1:length(bndry_edges_BC)
    for j=(i+1):length(bndry_edges_BC)
      if bndry_edges_BC[i] == bndry_edges_BC[j]
        dup_edge = bndry_edges_BC[i]
        #TODO; find which BCs it is assigned to for an even better error 
        #       message
        throw(ErrorException("Only 1 boundary condition may be applied to each geometric edge/face, edge/face $dup_edge appears twice"))
      end
    end
  end

  unused_geo_edge_nums = Array(Int, 0)

  for i=1:length(geo_edge_nums)
    edge_i = geo_edge_nums[i]
    if !(edge_i in bndry_edges_BC)
      push!(unused_geo_edge_nums, edge_i)
      continue
    end
  end

  @assert unused_geo_edge_nums == unique(unused_geo_edge_nums)

  return unused_geo_edge_nums
end

      


"""
  This function populates the mesh.bndryfaces field.

  Inputs:
    mesh: a mesh object
    boundary_nums: the 2 x mesh.numBoundaryFaces array containing the
                   [element number, global face number] for each
                   boundary faces

  Notes:
    the local face number is found using the *Pumi* topology.  It really should
    use the SBP topology, but they are the same for this purpose, so it
    hasn't caused a problem (yet)
"""
function getBoundaryArray(mesh::PumiMesh, boundary_nums::AbstractArray{Int, 2})
# get an array of type Boundary for SBP
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

        # get the element of the other edge
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
    centroidL = centroid1
    centroidR = centroid2
        else
    elementL = element2
    elementR = element1
    edgeL = edge2
    edgeR = edge1
    centroidL = centroid2
    centroidR = centroid1
        end
      else  # use y coordinate to decide
        if centroid1[2] < centroid2[2]
    elementL = element1
    elementR = element2
    edgeL = edge1
    edgeR = edge2
    centroidL = centroid1
    centroidR = centroid2
        else
    elementL = element2
    elementR = element1
    edgeL = edge2
    edgeR = edge1
    centroidL = centroid2
    centroidR = centroid1
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

"""
  This function populates the mesh.interfaces array

  Inputs/Outputs:
    mesh: a 3D DG mesh object.  The mesh.interfaces field should already have 
          been allocated to be of length mesh.numInterfaces.  This function
          populates it with one Interface object for each face shared by
          two elements, including periodic faces.  The face numbers are the
          local face numbers, and the orientation is determined using the
          SBP topology.

    Notes:
      The local face numbers use the Pumi topology.  This face numbering
      (but not orientation) is consistent between Pumi and the SummationByParts
      package, but this could be a problem if that ever changes.
      
"""
function getInterfaceArray(mesh::PumiMesh3D)
# Interface.orient is calculated using the *SBP* topology of the faces

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

      # Note: this must be consistent with the topo definition of face numbering
      #       (I think)
      localfacenumL = getFaceLocalNum(mesh, edgenumL, elnumL)
      localfacenumR = getFaceLocalNum(mesh, edgenumR, elnumR)
      fdata = FaceData(elnumL, elL, elnumR, elR, localfacenumL, localfacenumR,
                       vertsL, vertsR, facevertsL, facevertsR, part_nums, 
                       matched_entities, has_local_match)

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
