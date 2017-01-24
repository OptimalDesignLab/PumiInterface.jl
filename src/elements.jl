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
      throw(ErrorException("unsupported triangulation requested"))
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
      throw(ErrorException("unsupported triangulation requested"))
    end
  elseif shape_type == 3
#    if order  == 1
#      triangulation = Int32 [ 3 1 2;].'
#    else
       triangulation = zeros(Int32, 0, 0)
#      throw(ErrorException("unsupported triangulation requested"))
#    end
  else
      throw(ErrorException("unsupported triangulation requested"))
  end  # end if shape_type

  return triangulation
end

function getNodeMaps(order::Integer, shape_type::Integer, numNodesPerElement, dim=2)
# get the mappings between the SBP and Pumi node orderings
# having to do the mapping at all is inelegent to say the least
# store mappings in both directions in case they are needed
# use UInt8s to save cache space during loops

  if dim == 2
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

    elseif shape_type == 2 || shape_type == 3
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
  else  # dim == 3
    sbpToPumi = collect(UInt8, 1:numNodesPerElement)
    pumiToSbp = collect(UInt8, 1:numNodesPerElement)
  end

  return sbpToPumi, pumiToSbp
end  # end getNodeMaps

function createSubtriangulatedMesh(mesh::AbstractMesh, opts)
# this function is used to figure out how to do subtriangulation/visualization
# for 2D meshes only
# 3D meshes don't subtriangulate

  @assert mesh.dim == 2

  order = mesh.order
  shape_type = mesh.shape_type
  dofpernode = mesh.numDofPerNode

  if mesh.shape_type == 1
    if order >= 3
      mesh.mnew_ptr, mesh.fnew_ptr, mesh.fnewshape_ptr = createCGSubMesh(mesh)

    else
      mesh.triangulation = zeros(Int32, 0, 0)
      # maybe makes these equal f_ptr and m_ptr for consistency?
      mesh.mnew_ptr = C_NULL
      mesh.fnew_ptr = C_NULL
    end
#=
  elseif mesh.shape_type == 2
    if order >= 1
      mesh.mnew_ptr, mesh,fnew_ptr, mesh.fnewshape_ptr = createDGSubMesh(mesh)

    else
      throw(ErrorException("Congratulations: you have reached and unreachable case"))
    end
=#
  elseif mesh.shape_type == 3 || mesh.shape_type == 2  # DG SBP Omega or Gamma
    # I don't think subtriangulation will work in this case, so create a 
    # linear field on the existing mesh, and interpolate the solution onto it
    fshape_new = getFieldShape(0, 1, mesh.dim)

    mesh.mnew_ptr = mesh.m_ptr
    mesh.fnew_ptr = createPackedField(mesh.mnew_ptr, "solution_field_interp", dofpernode)
    mesh.fnewshape_ptr = fshape_new

    if mesh.shape_type == 2 && opts["exact_visualization"]
      mesh.mexact_ptr, mesh.fexact_ptr, mesh.fexactshape_ptr = createDGSubMesh(mesh)
    else
      mesh.mexact_ptr = C_NULL
      mesh.fexact_ptr = C_NULL
      mesh.fexactshape_ptr = C_NULL
    end

    if mesh.shape_type != 2 && opts["exact_visualization"]
      throw(ErrorException("exact visualization only supported for SBP Omega"))
    end


  else
    throw(ErrorException("Unsupported shape_type"))
  end  # end if mesh.shape_type 



  return nothing
end


function createCGSubMesh(mesh::PumiMeshCG)
  order = mesh.order
  shape_type = mesh.shape_type
  dofpernode = mesh.numDofPerNode
  
  mesh.triangulation = getTriangulation(order, shape_type)
  mnew_ptr = createSubMesh(mesh.m_ptr, mesh.triangulation, mesh.elementNodeOffsets, mesh.typeOffsetsPerElement_, mesh.entity_Nptrs)

  fnew_ptr = createPackedField(mnew_ptr, "solution_field", dofpernode)

  fnewshape_ptr = getMeshShapePtr(mnew_ptr)

  return mnew_ptr, fnew_ptr, fnewshape_ptr
end


function createDGSubMesh(mesh::PumiMeshDG)
  order = mesh.order
  shape_type = mesh.shape_type
  dofpernode = mesh.numDofPerNode

  mesh.triangulation = getTriangulation(order, shape_type)
  mnew_ptr = createSubMeshDG(mesh.m_ptr, mesh.mshape_ptr, mesh.triangulation, mesh.elementNodeOffsets, mesh.typeOffsetsPerElement_, mesh.nodemapPumiToSbp, mesh.entity_Nptrs, mesh.coords)

  fnew_ptr = createPackedField(mnew_ptr, "solution_field", dofpernode)
  fnewshape_ptr = getMeshShapePtr(mnew_ptr)

  return mnew_ptr, fnew_ptr, fnewshape_ptr
end

function transferFieldToSubmesh(mesh::AbstractMesh, u, mnew_ptr=mesh.mnew_ptr, 
                                fnew_ptr=mesh.fnew_ptr)

  if mesh.isInterpolated
    if mesh.shape_type != 3
      transferFieldDG(mesh.m_ptr, mnew_ptr, mesh.triangulation, 
                      mesh.elementNodeOffsets, mesh.typeOffsetsPerElement_, 
                      mesh.entity_Nptrs, mesh.f_ptr, mesh.interp_op.', fnew_ptr)
    else
      throw(ErrorException("Submesh not supported for SBPGamma"))
    end

      # this is a bit wasteful, because the un-interpolated solution was
      # already saved to mesh.f_ptr
      # interpolateToMesh(mesh, u)
    #end
  else  # CG mesh
    if mesh.order >= 3
      transferField(mesh.m_ptr, mesh.mnew_ptr, mesh.triangulation, mesh.elementNodeOffsets, mesh.typeOffsetsPerElement_, mesh.entity_Nptrs, mesh.f_ptr, mesh.fnew_ptr)
    end
  end

  return nothing
end

"""
  This function gets the xi coordinates of the nodes of a lagrangian element
  in order verts, edges, faces, regions

  Inputs:
    order: order of the element
    dim: dimensionality of the mesh (2d or 3d)

  Outputs:
    node_xi: a dim x numnodes array containing the xi coordinates of each node

  Note that this follows the Pumi convention that the reference element is 
  a right triangle with the right angle in the lower left corner.  Barycentric
  coordinates are used and span from 0 to 1.  The right angle corner has
  coordinate xi3, the next corner counter-clowckwise is xi1, and the next is
  xi2.  The array contains xi1 and xi2 along the first dimension.  However,
  the ordering of the vertices is counter-clockwise from the bottom left, 
  ie. the 1st vertex has coordinate xi3.  This is rather confusing and the
  source of many problems.

  The easiest way to see the Pumi convention is to look at the shape functions
  for the lagrangian elements in apfShape.cc
"""
function getXiCoords(order::Integer, dim::Integer)

  # TODO: write test
  # idea: verify that shape function i is == 1 at node i, all others zero
  if dim == 2
    if order == 1
      xicoords = [0. 1 0;
                  0 0 1]
    elseif order == 2
      xicoords = [0. 1 0 0.5 0.5 0.0;
                  0 0 1 0.0 0.5 0.5]
    else
      throw(ErrorException("unsupported order $order for dimension $dim"))
    end
  elseif dim == 3
    if order == 1
      xicoords = [0. 1 0 0;
                  0 0 1 0;
                  0  0 0 1]
    elseif order == 2
      xicoords = [0. 1 0 0 0.5 0.5 0.0 0.0 0.5 0.0;
                  0 0 1 0 0.0 0.5 0.5 0.0 0.0 0.5;
                  0 0 0 1 0.0 0.0 0.0 0.5 0.5 0.5]
    else
      throw(ErrorException("unsupported order $order for dimension $dim"))
    end
  else
    throw(ErrorException("unsupported dimension $dim"))
  end


end
