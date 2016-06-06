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


