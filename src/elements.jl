function getTriangulation(order::Integer, shape_type::Integer)
# get the sub-triangulation for an element of the specified order element
# the first 3 values indicate the verticies, later values refer to the nodes
# triangulation must be a 3 x n array of Int32s, so when it gets transposed
# by passing it to C, it becomes a n x 3 array of ints

  if shape_type == 1
    if order == 1  # no need to subtriangulate this case
      triangulation = zeros(Int32, 0, 0)
    elseif order == 2
      triangulation = Int32[1 1 4 2 5 6; 4 7 2 5 3 7; 7 6 7 7 7 3]
    elseif order == 3
      triangulation = Int32[1 1 4 5 5 2 6 7 10 12 12 10 9; 4 10 5 11 2 6 7 12 11 7 3 12 10; 10 9 10 10 11 11 11 11 12 3 8 8 8]
    elseif order == 4
      triangulation = Int32[1 1 4 5 5 6 6 2 7 7 15 14 14 8 8 9 17 11 18 13 13 18; 4 13 5 14 6 15 2 7 8 16 16 16 18 17 9 3 3 17 17 18 11 16; 13 12 13 13 14 14 15 15 16 15 14 18 13 16 17 17 10 10 11 11 12 17]
    else
      println(STDERR, "Warning, unsupported triangulation requested")
      triangulation = zeros(Int32, 0, 0)
    end
     
  elseif shape_type == 2
    if order == 1
      triangulation = Int32[1 2 2 6 1 5 5; 2 6 3 3 5 4 6; 5 5 6 4 3 3 4]
    elseif order == 2
      triangulation = Int32[1 8 4 1 1 9 2 5 6 1 8 6 1; 4 4 9 9 2 2 3 3 7 6 5 5 8; 8 5 5 4 9 5 5 7 3 3 6 7 6]
    elseif order == 3
      triangulation = Int32[7 1 8 8 8 7 13 7 5 12 8 5 13 13 12 11 12 10 3 11 6;12 5 1 3 4 3 8 4 10 6 13 1 12 9 13 13 11 1 12 1 11;3 9 9 1 3 4 9 8 9 2 7 10 7 10 11 10 6 11 2 2 2]
    #  triangulation = Int32[1 1 1 1 2 6 10 10 5 5 1 2 2 2 3 3 3 9 9 13; 10 11 6 2 12 12 11 12 10 13 5 7 4 3 8 9 1 7 7 13 12; 5 10 11 6 6 11 12 13 13 9 9 12 7 4 4 8 9 4 8 7 7]
    elseif order == 4
      triangulation = Int32[6 11 18 3 12 5 15 12 1 1 4 1 13 3 9 3 16 11 14 2 14 10 18 7 12 11 15 12 7 8 6;
      2 9 7 14 7 15 5 17 16 5 14 9 14 13 10 9 15 10 13 8 10 11 6 11 8 7 1 18 1 2 18;
      17 16 2 8 18 16 1 8 9 16 3 3 4 4 13 13 11 9 10 17 8 8 17 15 11 12 7 17 2 3 2]
    else
      println(STDERR, "Warning, unsupported triangulation requested")
      triangulation = zeros(Int32, 0, 0)
    end
  end  # end if shape_type

  return triangulation
end

function getNodeMaps(order::Integer, shape_type::Integer, numNodesPerElement)
# get the mappings between the SBP and Pumi node orderings
# having to do the mapping at all is inelegent to say the least
# store mappings in both directions in case they are needed
# use UInt8s to save cache space during loops

  if shape_type == 1
    if order == 1
      sbpToPumi = UInt8[1,2,3]
      pumiToSbp = UInt8[1,2,3]
    elseif order == 2
      sbpToPumi = UInt8[1,2,3,4,5,6,7]
      pumiToSbp = UInt8[1,2,3,4,5,6,7]
    elseif order == 3
      sbpToPumi = UInt8[1,2,3,4,5,6,7,9,8,12,10,11]
      pumiToSbp= UInt8[1,2,3,4,5,6,7,9,8,11,12,10]
    elseif order == 4 
      sbpToPumi = UInt8[1,2,3,4,5,6,7,8,9,12,11,10,17,13,15,14,16,18]
      pumiToSbp = UInt8[1,2,3,4,5,6,7,8,9,12,11,10,14,16,15,17,13,18]
    else
      println(STDERR, "Warning: Unsupported element order requestion in getNodeMaps")

      # default to 1:1 mapping
      sbpToPumi = UInt8[1:numNodesPerElement]
      pumiToSbp = UInt8[1:numNodesPerElement]
    end

  elseif shape_type == 2
    if order <= 4
      sbpToPumi = collect(UInt8, 1:numNodesPerElement)
      pumiToSbp = collect(UInt8, 1:numNodesPerElement)
    else

      println(STDERR, "Warning: Unsupported element order requestion in getFaceOffsets")
      # default to 1:1 mapping
      sbpToPumi = UInt8[1:numNodesPerElement;]
      pumiToSbp = UInt8[1:numNodesPerElement;]
    end
  end  # end if shape_type

  return sbpToPumi, pumiToSbp
end  # end getNodeMaps

function createSubtriangulatedMesh(mesh::AbstractMesh)
  order = mesh.order
  shape_type = mesh.shape_type
  dofpernode = mesh.numDofPerNode

  if mesh.shape_type == 1
    if order >= 3

      mesh.triangulation = getTriangulation(order, shape_type)
      flush(STDOUT)
      flush(STDERR)
      mesh.mnew_ptr = createSubMesh(mesh.m_ptr, mesh.triangulation, mesh.elementNodeOffsets, mesh.typeOffsetsPerElement_, mesh.entity_Nptrs)

      println("creating solution field on new mesh")
      mesh.fnew_ptr = createPackedField(mesh.mnew_ptr, "solution_field", dofpernode)
    else
      mesh.triangulation = zeros(Int32, 0, 0)
      mesh.mnew_ptr = C_NULL
      mesh.fnew_ptr = C_NULL
    end

  elseif mesh.shape_type == 2
    if order >= 1

      mesh.triangulation = getTriangulation(order, shape_type)
      flush(STDOUT)
      flush(STDERR)
      println("size(mesh.triangulation) = ", size(mesh.triangulation))
      mesh.mnew_ptr = createSubMeshDG(mesh.m_ptr, mesh.mshape_ptr, mesh.triangulation, mesh.elementNodeOffsets, mesh.typeOffsetsPerElement_, mesh.nodemapPumiToSbp, mesh.entity_Nptrs, mesh.coords)

      println("creating solution field on new mesh")
      mesh.fnew_ptr = createPackedField(mesh.mnew_ptr, "solution_field", dofpernode)
    else
      mesh.triangulation = zeros(Int32, 0, 0)
      mesh.mnew_ptr = C_NULL
      mesh.fnew_ptr = C_NULL
    end

  end  # end if mesh.shape_type 

  return nothing
end

function transferFieldToSubmesh(mesh::AbstractMesh)

  if mesh.isInterpolated
    println("transferring field to submesh")
    transferFieldDG(mesh.m_ptr, mesh.mnew_ptr, mesh.triangulation, mesh.elementNodeOffsets, mesh.typeOffsetsPerElement_, mesh.entity_Nptrs, mesh.f_ptr, mesh.interp_op.', mesh.fnew_ptr)
  else
    if mesh.order >= 3
      transferField(mesh.m_ptr, mesh.mnew_ptr, mesh.triangulation, mesh.elementNodeOffsets, mesh.typeOffsetsPerElement_, mesh.entity_Nptrs, mesh.f_ptr, mesh.fnew_ptr)
    end
  end

  return nothing
end
