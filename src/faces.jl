# functiosn related to data about faces (edges in 2D)

function getBoundaryArray(mesh::PumiMesh, boundary_nums::AbstractArray{Int, 2})
# get an array of type Boundary for SBP
# creating an an array of a user defined type seems like a waste of memory operations
# bnd_array is a vector of type Boundary with length equal to the number of edges on the boundary of the mesh

#  mesh.bndryfaces = Array(Boundary, mesh.numBoundaryFaces)

  for i=1:mesh.numBoundaryFaces
    facenum = boundary_nums[i, 1]
    edgenum_global = boundary_nums[i, 2]
#    println("edgenum_global = ", edgenum_global)
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
    smallmatvec!(alpha, view(sbp.facenormal, :, face_i), view(face_normals, :, j, i))

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
    smallmatvec!(alpha, view(sbp.facenormal, :, face_iL), view(face_normals, :, 1, j, i))
      # calculate right fae normal
      node_index = sbp.facenodes[j, face_iR]

      # calculate the mysterious parameter alpha
      for di1=1:2
	for di2=1:2
	  alpha[di1, di2] = (dxidx[di1,1, node_index, element_iR].*dxidx[di2,1, node_index, element_iR] + dxidx[di1,2, node_index, element_iR].*dxidx[di2,2, node_index, element_iR])*jac[node_index,element_iR]
	end
      end

    # call SBP function
    smallmatvec!(alpha, view(sbp.facenormal, :, face_iR), view(face_normals, :, 2, j, i))




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
#  num_ext_edges = size(getBoundaryEdgeNums(mesh))[1]  # bad memory efficiency
#  num_int_edges = getNumEdges(mesh) - num_ext_edges

#  new_bndry = Boundary(2, 3)
#   println("new_bndry = ", new_bndry)

#  new_interface = Interface(1, 2, 3, 4)
#   println("new_interface = ", new_interface)

#   println("num_int_edges = ", num_int_edges)

#  interfaces = Array(typeof(new_interface), num_int_edges)
  # get number of nodes affecting an edge
  num_edge_nodes = countAllNodes(mesh.mshape_ptr, 1)

  # unused variable?
  nodemap = Array(num_edge_nodes:(-1):1)

  pos = 1 # current position in interfaces
  for i=1:mesh.numEdge
    # get number of elements using the edge
    adjacent_nums, num_adjacent = getAdjacentEntityNums(mesh, i, 1, 2)
    if num_adjacent > 1  # internal edge
#       println("this is an internal edge")
      element1 = adjacent_nums[1]
      element2 = adjacent_nums[2]

#      coords_1 = x[:, :, element1]
#      coords_2 = x[:, :, element2]

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
        else
    elementL = element2
    elementR = element1
        end
      else  # use y coordinate to decide
        if centroid1[2] < centroid2[2]
    elementL = element1
    elementR = element2
        else
    elementL = element2
    elementR = element1
        end
      end

      edgeL = getFaceLocalNum(mesh, i, elementL)
      edgeR = getFaceLocalNum(mesh, i, elementR)

      # here we set the orientation flag to 1 for all edges, because in 2D
      # the edges always have opposite orientation
      mesh.interfaces[pos] = Interface(elementL, elementR, edgeL, edgeR, UInt8(1))
#       println("updating pos")
      pos += 1

  #    print("\n")

    end  # end if internal edge

#     print("\n")
  end  # end loop over edges

  return nothing

end  # end function


function getInterfaceArray(mesh::PumiMesh3D)
  println("----- entered getInterfaceArray -----")

  adj_elements = Array(Ptr{Void}, 2)
  coords1 = zeros(3, 4)
  coords2 = zeros(3, 4)
  centroid1 = zeros(3)
  centroid2 = zeros(3)
  pos = 1 # position in mesh.interfaces
  for i=1:mesh.numFace
    println("face ", i)
    face_i = mesh.faces[i]
    num_adjacent = countAdjacent(mesh.m_ptr, face_i, mesh.dim)

    if num_adjacent == 2  # this is a shared interface
      getAdjacent(adj_elements)
      el1 = adj_elements[1]
      el2 = adj_elements[2]
      elnum1 = getNumberJ(mesh.el_Nptr, el1, 0, 0) + 1
      elnum2 = getNumberJ(mesh.el_Nptr, el2, 0, 0) + 1

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
      else
        coordsL = coords1
        coordsR = coords2
        elnumL = elnum1
        elnumR = elnum2
      end

      whichL, flipL, rotateL = getAlignment(mesh.m_ptr, el1, face_i)
      r1 = EntityOrientation(whichL, flipL, rotateL)
      whichR, flipR, rotateR = getAlignment(mesh.m_ptr, el2, face_i)
      r2 = EntityOrientation(whichR, flipR, rotateR)

      rel_rotate = calcRelRotation(mesh, r1, r2)

      mesh.interfaces[pos] = Inteface(el1, el2, whichL, whichR, UInt8(rel_rotate))
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
