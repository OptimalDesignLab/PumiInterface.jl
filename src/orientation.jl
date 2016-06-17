# file for determining MeshEntity orientations
# not generalizable (but may not be necessary for DG)
function getEntityOrientations(mesh::PumiMesh)
# get the offset for node access of each mesh entity of each element
# to read/write from a Pumi field, accessing the data stored on an edge looks like:
# getComponents(f, e, abs(offset - node) - 1, vals)
# where f is the field, e is a mesh entity, node is the node index (one based),
# and vals is an array of values to be written to the field
# offset = 0 corresponds to the original ordering
# offset = number of nodes on the entity + 1 corresponds to reversing the 
#  ordering
# other values of offset should corresponds to rotations in 3d
# for 2d, the offsets for all nodes on an entity would typically be the same,
# corresponding to either flipping or not flipping the edge, but in 3d they
# might be different for a rotation
  # create array of arrays
  # hold an offset for every node of every entity of every element
  offsets = zeros(UInt8, mesh.numNodesPerElement, mesh.numEl)
  # hold a bit telling if the entity is in the proper orientation for this element
  # ie. it is "owned" by this element, so the nodes can be accessed in order

  # the elementNodeOffsets array is internal, therefore should be 
  # stored in Pumi node order, not SBP order
  flags = Array(BitArray{2}, 3)
  # create the arrays
  for i=1:3
    flags[i] = trues(mesh.numTypePerElement[i], mesh.numEl)
  end

  # now do edges

  edge_flags = flags[2]
  edges_i = Array(Ptr{Void}, 12)
  edges_start = mesh.typeOffsetsPerElement[2]
  nnodes_per_edge = mesh.numNodesPerType[2]
  for i=1:mesh.numEl
    getDownward(mesh.m_ptr, mesh.elements[i], 1, edges_i)
 
    for j=1:mesh.numTypePerElement[2]  # loop over edges on element i
      # get orientation of edge j on element i
      # get global number of edge j
      edgenum_global = getNumberJ(mesh.edge_Nptr, edges_i[j], 0, 0) + 1

      orient, edge_idx = getEdgeOrientation(mesh, i, edgenum_global)
      
      # figure out offset value
      # write n = mesh.numNodesPerType[2] + 1 if orientation = -1, or n=0 if orientation=1
      if orient == 1
	offset_val = 0
      else
        offset_val = mesh.numNodesPerType[2] + 1
      end

      # starting index for nodes on this edge in Pumi ordering
      start_idx = edges_start + (edge_idx - 1)*nnodes_per_edge

      # take into account Pumi to SBP mapping
      for k=1:mesh.numNodesPerType[2]  # loop over nodes on the edge
	# figure out the SBP index of this node
	node_idx = start_idx + k - 1  # pumi node_idx
        offsets[node_idx, i] = offset_val
      end

      # save value to flag array
      edge_flags[j, i] = div(orient + 1, 2)

    end  # end loop over edges
  end # end loop over elements

  # now do faces
#  getFaceOffsets(mesh, offsets, flags)

  return offsets, flags
end


function getEdgeOrientation(mesh::PumiMesh, elnum::Integer, edgenum::Integer)
# figure out what the orientation of the specified edge is relative ot the 
# element by looking at the order of the vertices that define the edge
# if the edge is in the same orientation as the element, return 1, otherwise -1
# because we are dealing with simplex elements, the third edge has to be 
# handled specially

#  println("\nEntered getEdgeOrientation")
  # get all the vertices
  el_verts, tmp = getDownward(mesh.m_ptr, mesh.elements[elnum], 0)
  edge_verts, tmp = getDownward(mesh.m_ptr, mesh.edges[edgenum], 0)
  el_edges, tmp = getDownward(mesh.m_ptr, mesh.elements[elnum], 1)

  # find out which edge of the face it is
  edge_idx = 0
  for i=1:3
    edgenum_ret = getNumberJ(mesh.edge_Nptr, el_edges[i], 0, 0) + 1
    if edgenum_ret == edgenum
      edge_idx = i
      break
    end
  end

  # an edge is defined by its two vertices
  # therefore, if the vertices are in the same order on the edge as on
  # the element, then the edge is in the same orientation as the element,
  # otherwise opposite orientation.
  # For simpleces, the 3rd edge has to be special cased because it uses
  # vertices 3 and 1
  pos1 = 0  # position of edge_verts[1] in el_verts
  pos2 = 0  # position of edge_verts[2] in el_verts

  # find the positions of the edge vertices on the element
  for i=1:3  # loop over el_verts
    if el_verts[i] == edge_verts[1]
      pos1 = i
    elseif el_verts[i] == edge_verts[2]
      pos2 = i
    end
  end

  @assert pos1 != 0
  @assert pos2 != 0

  # use the vertex positions to determine oridentation
  if edge_idx <= 2
    if pos2 - pos1 > 0  # positively oriented
      return 1, edge_idx
    elseif pos2 - pos1 < 0  # negatively oriented
      return -1, edge_idx
    else
      println(STDERR, "Warning, bad orientation determination in PdePumiInterface getEdgeOrientation")
      return 0, edge_idx
    end
  elseif edge_idx == 3  # ordering is reversed for 3rd edge
    if -(pos2 - pos1) > 0  # positively oriented
      return 1, edge_idx
    elseif -(pos2 - pos1) < 0  # negatively oriented
      return -1, edge_idx
    else
      println(STDERR, "Warning, bad orientation determination in PdePumiInterface getEdgeOrientation")
      return 0, edge_idx
    end
  else
    println(STDERR, "Warning, bad edge index determination inf PdePumiInterface getEdgeOrientation")
  end
 
end  # end function

immutable EntityOrientation
  which::Cint
  flip::Bool
  rotate::Cint
end


function getEntityOrientations(mesh::PumiMesh3DG)
# for DG (without subtriangulation), these shouldn't ever be needed

  offsets = zeros(UInt8, mesh.numNodesPerElement, mesh.numEl)
  flags = Array(BitArray{2}, 4)  # on array for Verts, edges, faces, regiosn
  for i=1:4
    flags[i] = falses(mesh.numTypePerElement[i], mesh.numEl)
  end

  return offsets, flags
end

function getRelRotate(mesh::PumiMesh3, elementL::Integer, elementR::Integer, facenum::Integer)
# calculate the rotation of faceR relative to faceL, where faceL and faceR
# are shared between the two element
# this uses getAligmnet, which describes how to transform the face from
# its current orientation to the canonical orientation as part of its parent
# element
# the difference in the required rotations is the relative rotation
# elementL : global element index
# elementR : global element index
# faceL : local face number of shared face on elementL
# faceR : local face number of shared face on elementR
# face : global face number of the shared face

eL = mesh.elements[elementL]
eR = mesh.elements[elementR]
face = mesh.faces[facenum]

# get rotations to bring faces into canonical orientation
# for each element
# we can ignore flip because one of the faces will always be flipped
whichL, flipL, rotateL = getAlignment(mesh.m_ptr, eL, face)
whichR, flipR, rotateR = getAlignment(mesh.m_ptr, eR, face)

r1 = EntityOrientation(whichL, flipL, rotateL)
r2 = EntityOrientation(whichR, flipR, rotateR)

rel_rotate = calcRelRotation(mesh, r1, r2)
return rel_rotate

end  # end function

# calculate the relative rotation of 2 faces
function calcRelRotation(mesh::PumiMesh, r1::EntityOrientation, 
                         r2::EntityOrientation, localfaces=true)

  if r1.flip == r2.flip && localfaces
    throw(ErrorException("Both faces cannot be flipped"))
  end

  # sum rotations because they each rotate ccw in their parent
  # element's orientation, so they rotate in opposite directions
  rel_rotate = r1.rotate + r2.rotate
  #println("rel_rotate before wrapping = ", rel_rotate)
  rel_rotate = wrapNumber(rel_rotate, 0, 2)
  #println("rel_rotate after wrapping = ", rel_rotate)

  # add 1 so output is in range [1,3]
  return rel_rotate + 1
end
 
function wrapNumber(num::Integer, lower::Integer, upper::Integer)
# make num perioid, where upper and lower are the max and min values allowable
# ie. num can only be in the range [1 3], if num == 0, then it gets mapped to 3
# similarly, 4 gets mapped to 1



range = upper - lower + 1

if num < lower
  diff = lower - num
  return upper - (diff % range) + 1
  
elseif  num > upper
  diff = num - upper
  return lower + (diff % range) - 1
else
  return num
end

end
"""
  Type to hold data about 2 elements that share a face, along with some 
  arrays needed by getRelativeOrientation
"""
immutable FaceData
  elnumL::Int
  elL::Ptr{Void}
  elnumR::Int
  elR::Ptr{Void}
  faceL::Int  # local face numbers
  faceR::Int
  vertsL::Array{Ptr{Void}, 1}
  vertsR::Array{Ptr{Void}, 1}
  facevertsL::Array{Ptr{Void}, 1}
  facevertsR::Array{Ptr{Void}, 1}
end

"""
  Maps the position where vertex 1 of faceL is found in faceR to the relative 
  rotation of the face.
"""
global const pos_to_rotation = [1, 2, 3]


"""
  This function takes the data from element that share a face and determines 
  the relative orientation of the shared face (ie. how many rotations cw 
  faceR is from faceL when looking from the perspective of faceL).  The 
  user supplied topology information is used to determine the ordering of the 
  vertices on both faces, therefore the relative orientation is consistent with
  the users definition of the faces
"""
function getRelativeOrientation(fdata::FaceData, mesh::PumiMesh3DG)

  faceL = fdata.faceL
  faceR = fdata.faceR
  vertmap = mesh.topo.face_verts

  # get all vertices of the tet
  getDownward(mesh.m_ptr, fdata.elL, 0, fdata.vertsL)
  getDownward(mesh.m_ptr, fdata.elR, 0, fdata.vertsR)

  # extract the face vertices 
  # uses the user suppied topology information to order the face verts
  for i=1:3
    fdata.facevertsL[i] = fdata.vertsL[vertmap[i, faceL]]
    fdata.facevertsR[i] = fdata.vertsR[vertmap[i, faceR]]
  end

  return calcRelativeOrientation(fdata.facevertsL, fdata.facevertsR)

end

function calcRelativeOrientation(facevertsL::AbstractArray, facevertsR::AbstractArray)

  # find vert1 of faceL in faceR
  idx = 0
  v1 = facevertsL[1]
  for i=1:3
    if facevertsL[i] == v1
      idx = i
      break
    end
  end

  return pos_to_rotation[idx]
end


function getRelativeOrientation(fdata::FaceData, mesh::PumiMesh2DG)

  return 1
end




