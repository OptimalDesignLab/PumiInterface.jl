# test 2D PdePumiInterface

using PumiConfig

@testset "--- Testing PdePumiInterface --- " begin

  test_math()

  # test masked copy
  mask = [1, 3]
  a = rand(4)
  a2 = zeros(2)
  PdePumiInterface.copy_masked(a2, a, mask)
  @test ( a2 )== a[mask]

  a = rand(3, 4)
  a2 = rand(3, 2)
  PdePumiInterface.copy_masked(a2, a, mask)
  @test ( a2 )== a[:, mask]

  a = rand(3, 2, 4)
  a2 = rand(3, 2,  2)
  PdePumiInterface.copy_masked(a2, a, mask)
  @test ( a2 )== a[:, :, mask]
  
  a = rand(3,2, 5, 4)
  a2 = rand(3, 2, 5, 2)
  PdePumiInterface.copy_masked(a2, a, mask)
  @test ( a2 )== a[:, :, :, mask]

  a = rand(2, 3)
  a2 = rand(1, 3)
  @test_throws Exception  PdePumiInterface.copy_masked(a2, a, mask)

  test_refcounting()




  opts = Dict{Any, Any}(
    "numBC" => 1,
    "BC1" =>  [0],
    "run_type" => 4,
    "verify_coloring" => true,
    "use_edge_res" => false,
    "write_edge_vertnums" => true,
    "write_face_vertnums" => true,
    "write_boundarynums" => true,
    "write_dxidx" => true,
    "write_coords" => true,
    "write_sparsity" => true,
    "write_offsets" => true,
    "write_counts" => true,
    "write_sparsity_nodebnds" => true,
    "write_dofs" => false,
    "write_interfaces" => false,
    "write_boundaries" => false,
    "write_sharedboundaries" => false,
    "use_linear_metrics" => true,
    )

    # names of output files
    fnames = ["boundary_nums", "face_vertnums", "edge_vertnums", "sparsity_bnds", "sparsity_bnds", "sparsity_nodebnds"]

    smb_name = "tri2l.smb"
    dmg_name = ".null"
    for order = 1:4
      println("testing order ", order, " CG mesh")
      sbp = getTriSBPGamma(degree=order)
#      sbp = TriSBP{Float64}(degree=order)
      ref_verts = [-1. 1 -1; -1 -1 1]
      sbpface = TriFace{Float64}(order, sbp.cub, ref_verts.')


    mesh =  PumiMesh2{Float64, typeof(sbpface)}(dmg_name, smb_name, order, sbp, opts, sbpface, coloring_distance=2, dofpernode=4)
    @test ( mesh.numVert )== 4
    @test ( mesh.numEdge )== 5
    @test ( mesh.numFace )== mesh.numEdge
    @test ( mesh.numEl )== 2
    @test ( mesh.numEntitiesPerType )== [4, 5, 2]
    @test ( mesh.numTypePerElement )== [3, 3, 1]
    @test ( mesh.numDofPerNode )== 4
    @test ( mesh.numBoundaryFaces )== 4
    @test ( mesh.numInterfaces )== 1

    @test ( mesh.bndryfaces[1].element )== 1
    @test ( mesh.bndryfaces[1].face )== 3
    @test ( mesh.bndryfaces[2].element )== 2
    @test ( mesh.bndryfaces[2].face )== 1
    @test ( mesh.bndryfaces[3].element )== 1
    @test ( mesh.bndryfaces[3].face )== 2
    @test ( mesh.bndryfaces[4].element )== 2
    @test ( mesh.bndryfaces[4].face )== 2

  #  println("mesh.interfaces = ",  mesh.interfaces)
    @test ( mesh.interfaces[1].elementL )== 1
    @test ( mesh.interfaces[1].elementR )== 2
    @test ( mesh.interfaces[1].faceL )== 1
    @test ( mesh.interfaces[1].faceR )== 3

    @test ( mesh.order )== order
    @test ( length(mesh.bndry_funcs) )== 1
    @test ( mesh.bndry_offsets )== [1, 5]
    @test ( mesh.bndry_geo_nums[1] )== opts["BC1"]

#=
    for i=1:mesh.numBoundaryFaces
      for j=1:(sum(mesh.numNodesPerType[1:2]))
        if i == 1
          @test isapprox( mesh.bndry_normals[:, j, i], [-1.0, 1.0]) atol=1e-13
        elseif i == 2
          @test isapprox( mesh.bndry_normals[:, j, i], [1.0, -1.0]) atol=1e-13
        elseif i == 3
          @test isapprox( mesh.bndry_normals[:, j, i], [0.0, 1.0]) atol=1e-13
        elseif i == 4
          @test isapprox( mesh.bndry_normals[:, j, i], [1.0, 0.0]) atol=1e-13
        end
      end
    end
=#
#=
    for i=1:mesh.numInterfaces
      for j=1:sbp.numfacenodes
        @test isapprox( mesh.interface_normals[:, 1, j, i], [1.0, -2.0]) atol=1e-13
        @test isapprox( mesh.interface_normals[:, 2, j, i], [-2.0, 1.0]) atol=1e-13
      end
    end
=#
    # verify that dofs on a node are numbered consecutively
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        start_dof = mesh.dofs[1, j, i]
        for k=1:mesh.numDofPerNode
          @test ( mesh.dofs[k, j, i] )== start_dof + k - 1 
        end
      end
    end

    @test ( mesh.color_masks[1][1] )== 1
    @test ( mesh.color_masks[1][2] )== 0
    @test ( mesh.color_masks[2][1] )== 0
    @test ( mesh.color_masks[2][2] )== 1

    @test ( mesh.neighbor_colors[1,1] )== 2
    @test ( mesh.neighbor_colors[2,1] )== 1
    @test ( mesh.neighbor_colors[1,2] )== 1
    @test ( mesh.neighbor_colors[2,2] )== 2

    @test ( mesh.neighbor_nums[1,1] )== 2
    @test ( mesh.neighbor_nums[2,1] )== 1
    @test ( mesh.neighbor_nums[2,1] )== 1
    @test ( mesh.neighbor_nums[2,2] )== 2

    @test ( mesh.pertNeighborEls[1, 1] )== 1
    @test ( mesh.pertNeighborEls[2, 1] )== 1
    @test ( mesh.pertNeighborEls[1, 2] )== 2
    @test ( mesh.pertNeighborEls[2, 2] )== 2

    @test ( mesh.color_cnt )== [1,1,0, 0]
    if order == 1
      @test ( mesh.numDof )== 16
      @test ( mesh.numNodes )== 4
      @test ( mesh.numNodesPerElement )== 3
      @test ( mesh.numNodesPerType )== [1, 0 , 0]
      @test ( mesh.typeOffsetsPerElement )== [1, 4, 4, 4]
    elseif order == 2
      @test ( mesh.numNodes )== 11
      @test ( mesh.numDof )== 44
      @test ( mesh.numNodesPerElement )== 7
      @test ( mesh.numNodesPerType )== [1, 1, 1]
      @test ( mesh.typeOffsetsPerElement )== [1, 4, 7, 8]
    elseif order == 3
      @test ( mesh.numNodes )== 20
      @test ( mesh.numDof )== 80
      @test ( mesh.numNodesPerElement )== 12
      @test ( mesh.numNodesPerType )== [1, 2, 3]
      @test ( mesh.typeOffsetsPerElement )== [1, 4, 10, 13]
      # check orientation of 3rd edge
#      println("mesh.dofs = ", mesh.dofs)
      @test ( reshape(mesh.dofs[1, 4:5, 1], 2) )== reshape(mesh.dofs[1, 8:9, 2], 2)
    elseif order == 4
      @test ( mesh.numNodes )== 31
      @test ( mesh.numDof )== 124
      @test ( mesh.numNodesPerElement )== 18
      @test ( mesh.numNodesPerType )== [1, 3, 6]
      @test ( mesh.typeOffsetsPerElement )== [1, 4, 13, 19]
#      println("mesh.dofs = ", mesh.dofs)
      # check orientation of 3rd edge
      @test ( reshape(mesh.dofs[1, 10:12, 2], 3) )==reshape(mesh.dofs[1, 4:6, 1], 3)
    end

    @test ( mesh.typeOffsetsPerElement_ )== mesh.typeOffsetsPerElement

    @test ( length(mesh.verts) )== mesh.numVert
    @test ( length(mesh.edges) )== mesh.numEdge
    @test ( length(mesh.elements) )== mesh.numEl

#    println("mesh.coords = ", mesh.coords)
    
    @test isapprox( mesh.jac, ones(mesh.numNodesPerElement ,2)) 

    fnames = ["boundary_nums", "face_vertnums", "edge_vertnums"]
    for name in fnames
      name_code = string(name, ".dat")
      name_ref = string(name, "_p", order, "true.dat")
      println("checking file ", name_code)
      println("against ", name_ref)
      data_code = readdlm(name_code)
      data_ref = readdlm(name_ref)


      for i=1:length(data_code)
        data_i = data_code[i]
        if typeof(data_i) <: Number
          @test isapprox( data_i, data_ref[i]) atol=1e-13
        else
          @test ( data_i )== data_ref[i]
        end
      end
    end

  # test vertmap
  println("testing vertmap")

  nedges_interior = 0
  for i=1:mesh.numEl
    vertnums_i = mesh.element_vertnums[:, i]
    for j=1:3
      @test  vertnums_i[j]  > 0
      @test  vertnums_i[j]  < mesh.numVert + 1
    end

    # test there are exactly 3 elements the current element shares 2 vertices
    # with
    nedges = 0
    for j=1:mesh.numEl
      vertnums_j = mesh.element_vertnums[:, j]
      if length(intersect(vertnums_i, vertnums_j)) == 2
        nedges += 1
        nedges_interior += 1
      end
    end

    @test  nedges  < 4
    @test  nedges  > 0
  end  # end loop i

  @test ( nedges_interior )== 2*(mesh.numEdge - mesh.numBoundaryFaces)


  # test getNodeXi
  if order <= 2
    fshape = apf.getFieldShape(0, order, 2)
    eshape = apf.getEntityShape(fshape, 2)  # triangle
    nodexi = PdePumiInterface.getXiCoords(order, 2)
    numnodes = size(nodexi, 2)
    for i=1:numnodes
      vals = apf.getValues(mesh.m_ptr, eshape, nodexi[:, i], numnodes)
      for j=1:numnodes
        if i == j
          @test  abs(vals[j] - 1)  < 1e-12
        else
          @test  abs(vals[j])  < 1e-12
        end  # end if else
      end   # end loop j
    end  # end loop i
  end  # end if p < 2


  # test getBoundary
  geo_bndries = opts["BC1"]
  bndries = getBoundaries(mesh, geo_bndries)
  bndries_range = mesh.bndry_offsets[1]:(mesh.bndry_offsets[2]-1)
  @test ( length(bndries) )== length(bndries_range)
  for i=1:length(bndries)
    bndry1 = bndries[i]
    bndry2 = mesh.bndryfaces[bndries_range[i]]
    @test ( bndry1.face )== bndry2.face
    @test ( bndry1.element )== bndry2.element
  end

  end  # end loop over p=1:4

  # verify adj and det are correct
  for i=1:10  # 10 random matrices
    A2 = rand(2,2)
    d = det(A2)
    @test isapprox( PdePumiInterface.det2(A2), d) atol=1e-12
    B2 = zeros(A2)
    PdePumiInterface.adjugate2(A2, B2)
    @test isapprox( B2./d, inv(A2)) atol=1e-12

    A3 = rand(3,3)
    d = det(A3)
    @test isapprox( PdePumiInterface.det3(A3), d) atol=1e-12
    B3 = zero(A3)
    PdePumiInterface.adjugate3(A3, B3)
    @test isapprox( B3/d, inv(A3)) atol=1e-12
  end

  # now test on a fully unstructured mesh
  smb_name = "vortex.smb"
  dmg_name = "vortex.dmg"

  opts = Dict{Any, Any}(
    "numBC" => 2,
    "BC1" =>  [4, 10],
    "BC2" =>  [7, 13],
    "run_type" => 4,
    "verify_coloring" => true,
    "use_edge_res" => false,
    "write_edge_vertnums" => true,
    "write_face_vertnums" => true,
    "write_boundarynums" => true,
    "write_dxidx" => true,
    "write_coords" => true,
    "write_sparsity" => true,
    "write_offsets" => true,
    "write_counts" => true,
    "write_sparsity_nodebnds" => true,
    "write_dofs" => false,
    "write_interfaces" => false,
    "write_boundaries" => false,
    "write_sharedboundaries" => false,
    "use_linear_metrics" => true,
    )



  for order = 1:4
    println("testing order ", order, " CG mesh against files")
    sbp = getTriSBPGamma(degree=order)
#    sbp = TriSBP{Float64}(degree=order)
    ref_verts = [-1. 1 -1; -1 -1 1]
    sbpface = TriFace{Float64}(order, sbp.cub, ref_verts.')


    mesh =  PumiMesh2{Float64, typeof(sbpface)}(dmg_name, smb_name, order, sbp, opts, sbpface, coloring_distance=2, dofpernode=4)

    
    for name in fnames
      name_code = string(name, ".dat")
      name_ref = string(name, "vortex", "_p", order, "true.dat")
      println("checking file ", name_code)
      println("against reference file ", name_ref)
      data_code = readdlm(name_code)
      data_ref = readdlm(name_ref)


      for i=1:length(data_code)
        data_i = data_code[i]
        if typeof(data_i) <: Number
          @test isapprox( data_i, data_ref[i]) atol=1e-13
        else
          @test ( data_i )== data_ref[i]
        end
      end
    end
  
  end   # end loop order=1:4


end

@testset "----- Testing PdePumiInterfaceDG -----" begin

  order = 1
  smb_name = "tri2l.smb"
  dmg_name = ".null"

  opts = Dict{Any, Any}(
    "dmg_name" => dmg_name,
    "smb_name" => smb_name,
    "order" => order,
    "coloring_distance" => 2,
    "numBC" => 1,
    "BC1" =>  [0],
    "run_type" => 4,
    "verify_coloring" => true,
    "use_edge_res" => false,
    "write_edge_vertnums" => true,
    "write_face_vertnums" => true,
    "write_boundarynums" => true,
    "write_dxidx" => true,
    "write_coords" => true,
    "write_sparsity" => true,
    "write_offsets" => true,
    "write_counts" => true,
    "write_sparsity_nodebnds" => true,
    "write_dofs" => false,
    "write_interfaces" => false,
    "write_boundaries" => false,
    "write_sharedboundaries" => false,
    "exact_visualization" => true,
    "use_linear_metrics" => true,
    )
#    interp_op = [0.5 0 0; 0 0.5 0; 0 0 0.5]

    sbp = getTriSBPOmega(degree=order)
#    sbp = TriSBP{Float64}(degree=order, internal=true)
    vtx = sbp.vtx
#    interp_op = SummationByParts.buildinterpolation(sbp, vtx.')

    sbpface = TriFace{Float64}(order, sbp.cub, vtx)
    println("sbpface.numnodes = ", sbpface.numnodes)

    # make a complex mesh because it is needed later
    # create it first so the underlying Pumi mesh gets destroyed and
    # replaced by the Float64 version
    mesh_c = PumiMeshDG2(Complex128, sbp, opts, sbpface, dofpernode=4)
#    mesh_c = PumiMeshDG2{Complex128, typeof(sbpface)}(dmg_name, smb_name, order, sbp, opts, sbpface, coloring_distance=2, dofpernode=4)

    mesh = PumiMeshDG2(Float64, sbp, opts, sbpface, dofpernode=4)
#    mesh =  PumiMeshDG2{Float64, typeof(sbpface)}(dmg_name, smb_name, order, sbp, opts, sbpface, coloring_distance=2, dofpernode=4)

   @test ( mesh.m_ptr )!=C_NULL
   @test ( mesh.mnew_ptr )!=C_NULL
   @test ( mesh.numVert )== 4
   @test ( mesh.numEdge )== 5
   @test ( mesh.numFace )== mesh.numEdge
   @test ( mesh.numEl )== 2
   @test ( mesh.numDof )== 24
   @test ( mesh.numNodes )== 6
   @test ( mesh.numDofPerNode )== 4
   @test ( mesh.numBoundaryFaces )== 4
   @test ( mesh.numInterfaces )== 1
   @test ( mesh.numNodesPerElement )== 3
   @test ( mesh.numNodesPerType )== [0, 0, 3]
   @test ( mesh.numEntitiesPerType )== [4, 5, 2]
   @test ( mesh.numTypePerElement )== [3, 3, 1]
   @test ( mesh.typeOffsetsPerElement )== [1, 1, 1, 4]
   @test ( mesh.typeOffsetsPerElement_ )== mesh.typeOffsetsPerElement
   @test isapprox( mesh.volume, 4.0) atol=1e-12
   volume2 = PdePumiInterface.calcVolume(mesh)
   @test isapprox( volume2, mesh.volume) atol=1e-12
   @test ( mesh.dim )== 2
   @test ( mesh.isDG )== true
   @test ( mesh.coloringDistance )== 2
   @test ( mesh.numColors )== 4
   @test ( mesh.numBC )== 1
   @test ( mesh.bndry_geo_nums[1] )== opts["BC1"]
   @test ( mesh.elementNodeOffsets )== zeros(mesh.numNodesPerElement, mesh.numEl)
   @test ( mesh.typeNodeFlags[1] )== trues(3, mesh.numEl)
   tmp = trues(3, 2); tmp[1, 1] = false
   @test ( mesh.typeNodeFlags[2] )== tmp
   @test ( mesh.bndry_offsets )== [1, 5]
   iface = mesh.interfaces[1]
   @test ( iface.elementL )== 1
   @test ( iface.elementR )== 2
   @test ( iface.faceL )== 1
   @test ( iface.faceR )== 3
   @test isapprox( mesh.coords[:, :, 1], [-2/3 -2/3 1/3; 2/3 -1/3 2/3]) atol=1e-13
   @test isapprox( mesh.coords[:, :, 2], [2/3 -1/3 2/3; 1/3 -2/3 -2/3]) atol=1e-13
   @test ( sort(unique(mesh.dofs)) )== collect(1:mesh.numDof)

   @test ( mesh.color_masks[1][1] )== true
   @test ( mesh.color_masks[1][2] )== false
   @test ( mesh.color_masks[2][1] )== false
   @test ( mesh.color_masks[2][2] )== true
   @test ( mesh.neighbor_colors[1, 1] )== 2
   @test ( mesh.neighbor_colors[2, 1] )== 1
   @test ( mesh.neighbor_colors[1, 2] )== 1
   @test ( mesh.neighbor_colors[2, 2] )== 2

   @test ( mesh.neighbor_nums[1, 1] )== 2
   @test ( mesh.neighbor_nums[2, 1] )== 1
   @test ( mesh.neighbor_nums[1, 2] )== 1
   @test ( mesh.neighbor_nums[2, 2] )== 2

   @test ( mesh.pertNeighborEls[1, 1] )== 1
   @test ( mesh.pertNeighborEls[2, 1] )== 1
   @test ( mesh.pertNeighborEls[1, 2] )== 2
   @test ( mesh.pertNeighborEls[2, 2] )== 2

   @test ( mesh.color_cnt[1] )== 1
   @test ( mesh.color_cnt[2] )== 1

   @test ( mesh.coord_order )== 1
   @test ( mesh.coord_numNodesPerElement )== 3
   @test ( length(mesh.typeOffsetsPerElement) )== 4
   @test ( size(mesh.coord_xi, 1) )== 2
   @test ( size(mesh.coord_xi, 2) )== 3
     
   # check that adjoint variables are right size
  @test ( size(mesh.dxidx) )== size(mesh.dxidx_bar)
  @test ( size(mesh.dxidx_bndry) )== size(mesh.dxidx_bndry_bar)
  @test ( size(mesh.dxidx_face) )== size(mesh.dxidx_face_bar)
  @test ( size(mesh.dxidx_sharedface_bar) )== size(mesh.dxidx_sharedface_bar)


  @test isapprox( mesh.jac, ones(mesh.numNodesPerElement ,2))
  test_coordNumbering(mesh)

  # check if dxidx is consistent with the old way of calculating it
  dxidx2 = zeros(mesh.dxidx)
  jac2 = zeros(mesh.jac)

  SummationByParts.mappingjacobian!(sbp, mesh.coords, dxidx2, jac2)

  @test isapprox( norm(mesh.jac - jac2)/length(mesh.jac), 0.0) atol=1e-13

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      @test isapprox( norm(mesh.dxidx[:, :, j, i]  - dxidx2[:, :, j, i]), 0.0) atol=1e-13
    end
  end

  u = ones(mesh.numDofPerNode, mesh.numTypePerElement[1], mesh.numEl)
  PdePumiInterface.saveNodalSolution(mesh, u)
  writeVisFiles(mesh, "nodal_solution")

  # test accumulation at vertices
  u_volume = ones(6, mesh.numTypePerElement[1], mesh.numEl)
  u_verts = ones(6, mesh.numVert)

  PdePumiInterface.accumulateAtVerts(mesh, u_volume, u_verts)

  # check that value at each vert is the number of elements using that vert
  for i=1:mesh.numVert
    vert_i = mesh.verts[i]
    nel = apf.countAdjacent(mesh.m_ptr, vert_i, mesh.dim)
    for j=1:6
      @test ( u_verts[j, i] )== nel
    end
  end

  # test surface numbering
  testSurfaceNumbering(mesh, sbp, opts)

  # check reverse mode
  # SBP testing the correctness, these tests only verify values get to the right place

  fill!(mesh.dxidx_bar, 1.0)
  PdePumiInterface.getVertCoords_rev(mesh, sbp)


  for i=1:mesh.numEl
    @test  norm(mesh.vert_coords_bar[:, :, i])  > 0.0
  end

  # make sure the data fields are the same
  fill!(mesh_c.vert_coords, -1)
  PdePumiInterface.copy_data!(mesh_c, mesh)
  compare_meshes(mesh, mesh_c)

  # test zero_bar_arrays
  fill!(mesh.vert_coords_bar, 1.0)
  zeroBarArrays(mesh)
  @test ( vecnorm(mesh.vert_coords_bar) )== 0.0



   function test_interp(mesh::AbstractMesh{Tmsh}) where Tmsh
     sbpface = mesh.sbpface
     dxdxi_element = zeros(2, 2, mesh.numNodesPerElement, 1)
     dxdxi_face = zeros(4, sbpface.numnodes, 1)
     dxidx_face = zeros(2,2, sbpface.numnodes)
     jac_face = zeros(sbpface.numnodes)

     for i=1:mesh.numInterfaces
       el = mesh.interfaces[i].elementL
       face = mesh.interfaces[i].elementR

       for j=1:mesh.numNodesPerElement
         dxidx_hat = mesh.dxidx[:, :, j, el]
         dxidx = dxidx_hat*mesh.jac[j, el]
         dxdxi = inv(dxidx)
         dxdxi_element[:, :, j, el] = dxdxi
       end

       bndry_arr = [ Boundary(1, face)]
       dxdxi_element_rshape = reshape(dxdxi_element, 4, mesh.numNodesPerElement, 1)
       boundaryinterpolate!(sbpface, bndry_arr, dxdxi_element_rshape, dxdxi_face)
       dxdxi_face_rshape = reshape(dxdxi_face, 2, 2, sbpface.numnodes)
    for j=1:sbpface.numnodes
      dxidx = inv(dxdxi_face_rshape[:, :, j])
      jac_face[j] = det(dxidx)
      dxidx_face[:, :, j] = dxidx/det(dxidx)

      @test isapprox( jac_face[j], mesh.jac_face[j, i]) atol=1e-13
      @test isapprox( dxidx_face[:, :, j], mesh.dxidx_face[:, :, j, i]) atol=1e-13
    end

    end  # end loop over interfaces

  end  # end function


  println("testing dxidx_face")
  if size(mesh.dxidx_face, 1) > 0
    println("running the test")
    test_interp(mesh)
  end

   # check dxidx_bndry
   # for straight sided elements, dxidx is constant
   println("testing dxidx_bndry")
   if size(mesh.jac_bndry, 1) > 0
     println("running the test")
     for i=1:mesh.numBoundaryFaces
       el = mesh.bndryfaces[i].element
       dxidx_test = mesh.dxidx[:, :, 1, el]
       jac_test = mesh.jac[1, el]
       for j=1:mesh.sbpface.numnodes
         dxidx_code = mesh.dxidx_bndry[:, :, j, i]
         jac_code = mesh.jac_bndry[j, i]
         @test isapprox( dxidx_code, dxidx_test) atol=1e-13
         @test isapprox( jac_code, jac_test) atol=1e-13
       end
     end
   end


  println("testing coords bndry")

   # check coords_bndry
   vert_coords = [-1.0 -1.0; 1  1; -1 1]
   
   function getBndry(arr, el, face)
     for i=1:length(arr)
       bndry_i = arr[i]
       if bndry_i.element == el && bndry_i.face == face
         return bndry_i, i
       end
     end
   end

   bndry, idx = getBndry(mesh.bndryfaces, 1, 2)
   # get vertex coordinates
   v1 = vert_coords[2, :]
   v2 = vert_coords[3, :]
   diff = v2 - v1
   slope = diff[2]/diff[1]
   b = v1[2] - slope*v1[1]

   # vertify the face node coordinates are on the line
   for i=1:mesh.sbpface.numnodes
     x_i = mesh.coords_bndry[1, i, idx]
     y_i = mesh.coords_bndry[2, i, idx]
    @test isapprox( y_i, slope*x_i + b) atol=1e-13
   end

   function test_normal_orientation(mesh, ifaces::AbstractArray{I}, nrm::AbstractArray{T, 3}) where {I <: Union{Boundary, Interface}, T}
     println("-----entered test_normal_orientation-----")
     topo = mesh.topo
     numVertPerElement = mesh.numTypePerElement[1]
     numVertPerFace = numVertPerElement - 1 
     println("numVertPerFace = ", numVertPerFace)
     el_verts = Array{Ptr{Void}}(numVertPerElement)
     other_vert_coords = zeros(mesh.dim)
     face_verts = Array{Ptr{Void}}(numVertPerElement - 1)
     face_vert_coords = zeros(mesh.dim, numVertPerFace)

     println("face_verts = ", topo.face_verts)
     nfaces = length(ifaces)
     for i=1:nfaces
       println("i = ", i)
       iface_i = ifaces[i]
       elnum = getElementL(iface_i)
       facenum_local = getFaceL(iface_i)

       el_i = mesh.elements[elnum]
       nverts = apf.getDownward(mesh.m_ptr, el_i, 0, el_verts)
       println("nverts = ", nverts)

       for j=1:numVertPerFace
         face_verts[j] = el_verts[topo.face_verts[j, facenum_local]]
         tmp = zeros(3)
         apf.getPoint(mesh.m_ptr, face_verts[j], 0, tmp)
         face_vert_coords[1:mesh.dim, j] = tmp[1:mesh.dim]
       end

       # get the vert not on the face
       other_vert = Ptr{Void}(0)
       for j=1:numVertPerElement
         if !(el_verts[j] in face_verts)
           other_vert = el_verts[j]
         end
       end

       tmp = zeros(3)
       apf.getPoint(mesh.m_ptr, other_vert, 0, tmp)
       other_vert_coords[1:mesh.dim] = tmp[1:mesh.dim]

       # check that the face normal is in the opposite direction as the
       # vectors from a vertex on the face to the vertex not on the face
       for j=1:mesh.numNodesPerFace
         for k=1:numVertPerFace
           r1 = other_vert_coords - face_vert_coords[:, k]

           val = dot(nrm[:, j, i], r1)
           @test  val  < 0.0
         end
       end

     end  # end loop i

     return nothing


   end  # end function test_normal_orientation

   test_normal_orientation(mesh, mesh.interfaces, mesh.nrm_face)
   test_normal_orientation(mesh, mesh.bndryfaces, mesh.nrm_bndry)


   println("testing number nodes windy")
   # check adjacency apf.reordering algorithm doesn't error out
   PdePumiInterface.numberNodesWindy(mesh, [0.0, 0.0, 0.0])

    opts["smb_name"] = "tri8l.smb"
    mesh = PumiMeshDG2(Float64, sbp, opts, sbpface, dofpernode=4)
#    mesh =  PumiMeshDG2{Float64, typeof(sbpface)}(dmg_name, smb_name, order, sbp, opts, sbpface, coloring_distance=2, dofpernode=4)

  # check mapping interpolation
  # should be constant within an element for straight-sided elementsa
  if size(mesh.jac_face, 1) > 0
    for i=1:mesh.numInterfaces
      iface_i = mesh.interfaces[i]
      el_i = iface_i.elementL
      dxidx_el = mesh.dxidx[:, :, 1, el_i]
      jac_el = mesh.jac[:, el_i]
      jac_face = mesh.jac_face[:, i]

      for j=1:mesh.numNodesPerFace
        dxidx_face = mesh.dxidx_face[:, :, j, i]
        for k=1:2
          for p=1:2
            @test isapprox( dxidx_face[p, k], dxidx_el[p, k]) atol=1e-13
          end
        end
        @test isapprox( jac_face[j], jac_el[j]) atol=1e-13
      end
    end  # end loop over interfaces
  end

  # test reverse mode interpolation

  test_interp_rev(mesh)
  test_coords_rev(mesh, sbp)

  # test update_coords
  println("testing update_coords")
  coords_orig = copy(mesh.vert_coords)
  
  # double all coordinates
  coords_i = zeros(Float64, 2, 3)
  for i=1:mesh.numEl
    for j=1:3
      coords_i[1, j] = 2*coords_orig[1, j, i]
      coords_i[2, j] = 2*coords_orig[2, j, i]
    end
    update_coords(mesh, i, coords_i)
  end

  commit_coords(mesh, sbp, opts)

#  PdePumiInterface.getCoordinatesAndMetrics(mesh, sbp)
  for i = 1:length(mesh.vert_coords)
    @test  abs(2*coords_orig[i] - mesh.vert_coords[i])  < 1e-10
  end

  # test vertmap
  println("testing vertmap")
  nedges_interior = 0
  for i=1:mesh.numEl
    vertnums_i = mesh.element_vertnums[:, i]
    for j=1:3
      @test  vertnums_i[j]  > 0
      @test  vertnums_i[j]  < mesh.numVert + 1
    end

    # test there are exactly 3 elements the current element shares 2 vertices
    # with
    nedges = 0
    for j=1:mesh.numEl
      vertnums_j = mesh.element_vertnums[:, j]
      if length(intersect(vertnums_i, vertnums_j)) == 2
        nedges += 1
        nedges_interior += 1
      end
    end

    @test  nedges  < 4
    @test  nedges  > 0
  end  # end loop i

  @test ( nedges_interior )== 2*(mesh.numEdge - mesh.numBoundaryFaces)

  # test saveSolutionToMesh interpolation
  println("testing saveSolutionToMesh interpolation")
  u_vals = zeros(mesh.numDof)
  for i=1:mesh.numDof

    # get dof, node, element
    idx = findfirst(mesh.dofs, i)
    dofidx, node, el = ind2sub(mesh.dofs, idx)

    x = mesh.coords[1, node, el]
    y = mesh.coords[2, node, el]
    order = mesh.order

    u_vals[i] = x^order + y^order + 1
  end

  saveSolutionToMesh(mesh, u_vals)
  writeVisFiles(mesh, "dg_vis_test")

  # check that the solution is interpolated exactly
  down_verts = Array{Ptr{Void}}(mesh.numTypePerElement[1])
  coords_vert = zeros(Float64, 3)
  interp_vals = zeros(mesh.numDofPerNode)
  for i=1:mesh.numEl
    order = mesh.order
    el_ptr = mesh.elements[i]
    apf.getDownward(mesh.m_ptr, el_ptr, 0, down_verts)

    for j=1:mesh.numTypePerElement[1]  # loop over vertices
      vert_j = down_verts[j]
      apf.getPoint(mesh.m_ptr, vert_j, 0, coords_vert)

      x = coords_vert[1]
      y = coords_vert[2]

      # note: this assumes mnew = m
      apf.getComponents(mesh.fnew_ptr, vert_j, 0, interp_vals)

      val_expected = x^order + y^order + 1
      for k=1:mesh.numDofPerNode
        @test  abs(interp_vals[k] - val_expected)  < 1e-12
      end
    end
  end

  # test retrieveing solution from mesh
  u_vals2 = zeros(u_vals)
  PdePumiInterface.retrieveSolutionFromMesh_interp(mesh, u_vals2)

  for i=1:length(u_vals2)
    @test abs(u_vals2[i] - u_vals[i]) < 1e-12
  end


  # just for good measure, create a new mesh
  mesh = PumiMeshDG2(Float64, sbp, opts, sbpface, dofpernode=4)
#  mesh =  PumiMeshDG2{Float64, typeof(sbpface)}(dmg_name, smb_name, order, sbp, opts, sbpface, coloring_distance=2, dofpernode=4)
  mesh_c = PumiMeshDG2(Complex128, sbp, opts, sbpface, dofpernode=4)
#  mesh_c =  PumiMeshDG2{Complex128, typeof(sbpface)}(dmg_name, smb_name, order, sbp, opts, sbpface, coloring_distance=2, dofpernode=4)

  compare_meshes(mesh, mesh_c)

  # test periodic
  println("testing periodic")
  @test ( mesh.numPeriodicInterfaces )== 0

  opts["smb_name"] = "tri3_px.smb"
  opts["BC1"] = [0, 2]
  mesh = PumiMeshDG2(Float64, sbp, opts, sbpface, dofpernode=4)
  mesh_c = PumiMeshDG2(Complex128, sbp, opts, sbpface, dofpernode=4)
#  mesh = PumiMeshDG2{Float64, typeof(sbpface)}(dmg_name, smb_name, order, sbp, opts, sbpface, coloring_distance=2, dofpernode=4)

  @test ( mesh.numPeriodicInterfaces )== 3
  @test ( length(mesh.interfaces) )== 24
  @test ( mesh.numInterfaces )== 24
  @test ( mesh.numBoundaryFaces )== 6
  @test ( opts["numBC"] )== 1  # the default BC should not be created for periodic

  for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]
    @test  iface_i.elementL  > 0
    @test  iface_i.elementL  < mesh.numEl + 1
    @test  iface_i.elementR  > 0
    @test  iface_i.elementR  < mesh.numEl + 1
    @test  iface_i.faceL  < 4
    @test  iface_i.faceR  > 0
    @test  iface_i.faceR  < 4
  end

  test_coordNumbering(mesh)

  
  # test curvilinear
  println("testing curvilinear")
  # check curvilinear routines reproduce linear results
  coords_face_orig = copy(mesh.coords_interface)
  nrm_face_orig = copy(mesh.nrm_face)

  PdePumiInterface.getCoordinates(mesh, sbp)
  PdePumiInterface.getFaceCoordinatesAndNormals(mesh, sbp)
  PdePumiInterface.getCurvilinearCoordinatesAndMetrics(mesh, sbp)

  for i=1:mesh.numInterfaces
    for j=1:mesh.numNodesPerFace
      @test isapprox( mesh.coords_interface[:, j, i], coords_face_orig[:, j, i]) atol=1e-12
      @test isapprox( mesh.nrm_face[:, j, i], nrm_face_orig[:, j, i]) atol=1e-12
    end
  end

  # test metrics reverse
  test_metrics_rev(mesh, mesh_c, sbp, opts)



  # a 0 - 5 square that used a sin wave to remap the nondimensionalized
  # coordinates
  opts["smb_name"] = "square_05_curve.smb"
  opts["use_linear_metrics"] = false
  opts["BC1"] = [0, 1, 2, 3]
  mesh = PumiMeshDG2(Float64, sbp, opts, sbpface, dofpernode=4)
  mesh_c = PumiMeshDG2(Complex128, sbp, opts, sbpface, dofpernode=4)
#  mesh =  PumiMeshDG2{Float64, typeof(sbpface)}(dmg_name, smb_name, order, sbp, opts, sbpface, coloring_distance=2, dofpernode=4)

  testSurfaceNumbering(mesh, sbp, opts)
  test_coordNumbering(mesh)
  test_coord_field(mesh)
  # TODO: re-enable when Pumi's ENABLE_ZOLTAN flag is available
  #test_split(mesh)
  test_submesh_transfer(mesh)

  #TODO: check sizes of arrays

  function test_volume_curvilinear(mesh, sbp)

    volume = 0.0
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        dxidx_scaled = mesh.dxidx[:, :, j, i]
        dxidx = dxidx_scaled*mesh.jac[j, i]

        dxdxi = inv(dxidx)
        jac = det(dxdxi)
        volume += sbp.w[j]*jac
      end
    end

    println("volume = ", volume)
    @test isapprox( volume, 25.0) atol=1e-12
  end  # end function

  test_volume_curvilinear(mesh, sbp)


  # translate the mesh and verify all quantities are the same
  dxidx_orig = copy(mesh.dxidx)
  jac_orig = copy(mesh.jac)
  nrm_bndry_orig = copy(mesh.nrm_bndry)
  nrm_face_orig = copy(mesh.nrm_face)

  for i=1:mesh.numEl
    coords_i = mesh.vert_coords[:, :, i]
    coords_i += 1
    update_coords(mesh, i, coords_i)
  end

  commit_coords(mesh, sbp, opts)

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      for k=1:mesh.dim
        for p=1:mesh.dim
          @test isapprox( mesh.dxidx[p, k, j, i], dxidx_orig[p, k, j, i]) atol=1e-12
        end
      end

      @test isapprox( mesh.jac[j, i], jac_orig[j, i]) atol=1e-12
    end
  end

  for i=1:mesh.numBoundaryFaces
    for j=1:mesh.numNodesPerFace
      for k=1:mesh.dim
        @test isapprox( mesh.nrm_bndry[k, j, i], nrm_bndry_orig[k, j, i]) atol=1e-12
      end
    end
  end

  for i=1:mesh.numInterfaces
    for j=1:mesh.numNodesPerFace
      for k=1:mesh.dim
        @test isapprox( mesh.nrm_face[k, j, i], nrm_face_orig[k, j, i]) atol=1e-12
      end
    end
  end

  # recalculate and verify all quantities are the same
  recalcCoordinatesAndMetrics(mesh, sbp, opts)

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      for k=1:mesh.dim
        for p=1:mesh.dim
          @test isapprox( mesh.dxidx[p, k, j, i], dxidx_orig[p, k, j, i]) atol=1e-12
        end
      end

      @test isapprox( mesh.jac[j, i], jac_orig[j, i]) atol=1e-12
    end
  end

  for i=1:mesh.numBoundaryFaces
    for j=1:mesh.numNodesPerFace
      for k=1:mesh.dim
        @test isapprox( mesh.nrm_bndry[k, j, i], nrm_bndry_orig[k, j, i]) atol=1e-12
      end
    end
  end

  for i=1:mesh.numInterfaces
    for j=1:mesh.numNodesPerFace
      for k=1:mesh.dim
        @test isapprox( mesh.nrm_face[k, j, i], nrm_face_orig[k, j, i]) atol=1e-12
      end
    end
  end

  # test metrics reverse
  test_metrics_rev(mesh, mesh_c, sbp, opts)
  test_metrics_rev_1d(mesh_c, sbp, opts)

      

  mesh = PumiMeshDG2(Float64, sbp, opts, sbpface, dofpernode=4)
#  mesh =  PumiMeshDG2{Float64, typeof(sbpface)}(dmg_name, smb_name, order, sbp, opts, sbpface, coloring_distance=2, dofpernode=4)
  mesh_c = PumiMeshDG2(Complex128, sbp, opts, sbpface, dofpernode=4)
#  mesh_c =  PumiMeshDG2{Complex128, typeof(sbpface)}(dmg_name, smb_name, order, sbp, opts, sbpface, coloring_distance=2, dofpernode=4)

  compare_meshes(mesh, mesh_c)

  test_submesh()

  test_adapt_2d()

  test_ScatterData(mesh)

  # do this last since it rewrites the mesh coordinate field
  test_setPoint(mesh, opts)

  # test high order DG interpolation to Pumi fields
  test_ho_interpolation()

  opts["smb_name"] = "meshes/tri4x4_.smb"
  opts["dmg_name"] = "meshes/tri4x4_.dmg"
  mesh = PumiMeshDG2(Float64, sbp, opts, sbpface, dofpernode=4)
  test_geoNums(mesh)

  println("\nTesting coords file")
  test_coords_file(opts, sbp, sbpface)

  if HAVE_SIMMETRIX

    # load airfoil CAD mesh
    opts = Dict{Any, Any}(
    "dimensions" => 2,
    "run_type" => 5,
    "jac_type" => 2,
    "order" => 1,
    "use_DG" => true,
    "coloring_distance" => 2,
    "numBC" => 2,
    "BC1" => [40],
    "BC1_name" => "FreeStreamBC",
    "BC2" => [5],
    "BC2_name" => "noPenetrationBC",
    "smb_name" => "meshes/airfoil2_.smb",
    "dmg_name" => "meshes/airfoil2.smd",
    )

    mesh = PumiMeshDG2(Float64, sbp, opts, sbpface, dofpernode=4)
    test_geoMapping(mesh)
    test_geoWrapping(mesh)

    # the airfoil has large curvature at the trailing edge, so its hard to find
    # a finite difference step size that works.  Use a smoother shape instead.

    # load airfoil CAD mesh
    opts = Dict{Any, Any}(
    "dimensions" => 2,
    "run_type" => 5,
    "jac_type" => 2,
    "order" => 1,
    "use_DG" => true,
    "coloring_distance" => 2,
    "numBC" => 1,
    "BC1" => [4, 7, 10, 13],
    "BC1_name" => "FreeStreamBC",
    "smb_name" => "meshes/UnitSquare/UnitSquare.smb",
    "dmg_name" => "meshes/UnitSquare/UnitSquare.x_t",
    )
    
    mesh = PumiMeshDG2(Float64, sbp, opts, sbpface, dofpernode=4)
    @time test_geoDerivative(mesh)
    println("test_geoDerivative @time printed above")

    println("\n\nTesting geometric derivative on curved geometry")
    # load smooth CAD mesh
    opts = Dict{Any, Any}(
    "dimensions" => 2,
    "run_type" => 5,
    "jac_type" => 2,
    "order" => 1,
    "use_DG" => true,
    "coloring_distance" => 2,
    "numBC" => 1,
    "BC1" => [4, 7, 10, 13],
    "BC1_name" => "FreeStreamBC",
    "smb_name" => "meshes/UnitSquareCurve/UnitSquareCurve_linear.smb",
    "dmg_name" => "meshes/UnitSquareCurve/UnitSquareCurve.x_t",
    )
 
    mesh = PumiMeshDG2(Float64, sbp, opts, sbpface, dofpernode=4)
    test_geoDerivative(mesh)
  
    test_coords_file(opts, sbp, sbpface)

    xivec = getXiCoords(mesh)
    xivec_pert = rand(length(xivec))
    for i=1:length(xivec)
      xivec[i] += 1e-3*xivec[i]
    end
    update_coordsXi(mesh, sbp, opts, xivec)
    test_geoDerivative(mesh)
    test_geoWarp(mesh, sbp, opts)
    test_update_coords(mesh, sbp, opts)

    # this following test messes up the coordinate field, so load a new mesh
    # for this test
    mesh = PumiMeshDG2(Float64, sbp, opts, sbpface, dofpernode=4)
    test_setPoint(mesh, opts)


    #TODO: test quadratic mesh once that works
 
    opts = Dict{Any, Any}(
    "dimensions" => 2,
    "run_type" => 5,
    "jac_type" => 2,
    "order" => 1,
    "use_DG" => true,
    "coloring_distance" => 2,
    "numBC" => 1,
    "BC1" => [4, 7, 10, 13],
    "BC1_name" => "FreeStreamBC",
    "smb_name" => "meshes/UnitSquareCurve/UnitSquareCurve_quadratic.smb",
    "dmg_name" => "meshes/UnitSquareCurve/UnitSquareCurve.x_t",
    )

    mesh = PumiMeshDG2(Float64, sbp, opts, sbpface, dofpernode=4)
    println("testing quadratic mesh")
    test_setPoint(mesh, opts)
    test_coords_file(opts, sbp, sbpface)

    # load new mesh because test_setPoint messes up the coordinate field
    mesh = PumiMeshDG2(Float64, sbp, opts, sbpface, dofpernode=4)
    test_geoDerivative(mesh)

   
    # test mesh with geometry that has an internal edge
    println("\n\nTesting internal edge")
    opts = Dict{Any, Any}(
    "dimensions" => 2,
    "run_type" => 5,
    "jac_type" => 2,
    "order" => 1,
    "use_DG" => true,
    "coloring_distance" => 2,
    "numBC" => 1,
    "BC1" => [7, 56, 60, 13, 4, 10, 16],
    "BC1_name" => "FreeStreamBC",
    "smb_name" => "meshes/corner20_split/corner20_split_mesh0_.smb",
    "dmg_name" => "meshes/corner20_split/corner20_split.smd",
    )

    mesh = PumiMeshDG2(Float64, sbp, opts, sbpface, dofpernode=4)
    test_interioredge(mesh)
  end 

end
