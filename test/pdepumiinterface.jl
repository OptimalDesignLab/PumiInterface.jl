
facts("--- Testing PdePumiInterface --- ") do

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
    "write_dofs" => false
    )

    # names of output files
    fnames = ["boundary_nums", "face_vertnums", "edge_vertnums", "sparsity_bnds", "sparsity_bnds", "sparsity_nodebnds"]

    smb_name = "tri2l.smb"
    dmg_name = ".null"
    for order = 1:4
      println("testing order ", order, " mesh")
      sbp = TriSBP{Float64}(degree=order)

    mesh =  PumiMesh2{Float64}(dmg_name, smb_name, order, sbp, opts, coloring_distance=2, dofpernode=4)
    @fact mesh.numVert --> 4
    @fact mesh.numEdge --> 5
    @fact mesh.numEl --> 2
    @fact mesh.numEntitiesPerType --> [4, 5, 2]
    @fact mesh.numTypePerElement --> [3, 3, 1]
    @fact mesh.numDofPerNode --> 4
    @fact mesh.numBoundaryEdges --> 4
    @fact mesh.numInterfaces --> 1

    @fact mesh.bndryfaces[1].element --> 1
    @fact mesh.bndryfaces[1].face --> 3
    @fact mesh.bndryfaces[2].element --> 2
    @fact mesh.bndryfaces[2].face --> 1
    @fact mesh.bndryfaces[3].element --> 1
    @fact mesh.bndryfaces[3].face --> 2
    @fact mesh.bndryfaces[4].element --> 2
    @fact mesh.bndryfaces[4].face --> 2

  #  println("mesh.interfaces = ",  mesh.interfaces)
    @fact mesh.interfaces[1].elementL --> 1
    @fact mesh.interfaces[1].elementR --> 2
    @fact mesh.interfaces[1].faceL --> 1
    @fact mesh.interfaces[1].faceR --> 3

    @fact mesh.order --> order
    @fact length(mesh.bndry_funcs) --> 1
    @fact mesh.bndry_offsets --> [1, 5]

    for i=1:mesh.numBoundaryEdges
      for j=1:(sum(mesh.numNodesPerType[1:2]))
        if i == 1
          @fact mesh.bndry_normals[:, j, i] --> roughly([-1.0, 1.0], atol=1e-13)
        elseif i == 2
          @fact mesh.bndry_normals[:, j, i] --> roughly([1.0, -1.0], atol=1e-13)
        elseif i == 3
          @fact mesh.bndry_normals[:, j, i] --> roughly([0.0, 1.0], atol=1e-13)
        elseif i == 4
          @fact mesh.bndry_normals[:, j, i] --> roughly([1.0, 0.0], atol=1e-13)
        end
      end
    end


    for i=1:mesh.numInterfaces
      for j=1:sbp.numfacenodes
        @fact mesh.interface_normals[:, 1, j, i] --> roughly([1.0, -2.0], atol=1e-13)
        @fact mesh.interface_normals[:, 2, j, i] --> roughly([-2.0, 1.0], atol=1e-13)
      end
    end

    # verify that dofs on a node are numbered consecutively
    for i=1:mesh.numEl
      for j=1:mesh.numNodesPerElement
        start_dof = mesh.dofs[1, j, i]
        for k=1:mesh.numDofPerNode
          @fact mesh.dofs[k, j, i] --> start_dof + k - 1 
        end
      end
    end

    @fact mesh.color_masks[1][1] --> 1
    @fact mesh.color_masks[1][2] --> 0
    @fact mesh.color_masks[2][1] --> 0
    @fact mesh.color_masks[2][2] --> 1

    @fact mesh.neighbor_colors[1,1] --> 2
    @fact mesh.neighbor_colors[2,1] --> 1
    @fact mesh.neighbor_colors[1,2] --> 1
    @fact mesh.neighbor_colors[2,2] --> 2

    @fact mesh.neighbor_nums[1,1] --> 2
    @fact mesh.neighbor_nums[2,1] --> 1
    @fact mesh.neighbor_nums[2,1] --> 1
    @fact mesh.neighbor_nums[2,2] --> 2

    println("mesh.pertNeighborEls = ", mesh.pertNeighborEls)
    @fact mesh.pertNeighborEls[1, 1] --> 1
    @fact mesh.pertNeighborEls[2, 1] --> 1
    @fact mesh.pertNeighborEls[1, 2] --> 2
    @fact mesh.pertNeighborEls[2, 2] --> 2

    @fact mesh.color_cnt --> [1,1,0, 0]
    if order == 1
      @fact mesh.numDof --> 16
      @fact mesh.numNodes --> 4
      @fact mesh.numNodesPerElement --> 3
      @fact mesh.numNodesPerType --> [1, 0 , 0]
      @fact mesh.typeOffsetsPerElement --> [1, 4, 4, 4]
    elseif order == 2
      @fact mesh.numNodes --> 11
      @fact mesh.numDof --> 44
      @fact mesh.numNodesPerElement --> 7
      @fact mesh.numNodesPerType --> [1, 1, 1]
      @fact mesh.typeOffsetsPerElement --> [1, 4, 7, 8]
    elseif order == 3
      @fact mesh.numNodes --> 20
      @fact mesh.numDof --> 80
      @fact mesh.numNodesPerElement --> 12
      @fact mesh.numNodesPerType --> [1, 2, 3]
      @fact mesh.typeOffsetsPerElement --> [1, 4, 10, 13]
      # check orientation of 3rd edge
      println("mesh.dofs = ", mesh.dofs)
      @fact reshape(mesh.dofs[1, 4:5, 1], 2) --> reshape(mesh.dofs[1, 8:9, 2], 2)
    elseif order == 4
      @fact mesh.numNodes --> 31
      @fact mesh.numDof --> 124
      @fact mesh.numNodesPerElement --> 18
      @fact mesh.numNodesPerType --> [1, 3, 6]
      @fact mesh.typeOffsetsPerElement --> [1, 4, 13, 19]
      println("mesh.dofs = ", mesh.dofs)
      # check orientation of 3rd edge
      @fact reshape(mesh.dofs[1, 10:12, 2], 3) -->reshape(mesh.dofs[1, 4:6, 1], 3)
    end

    @fact mesh.typeOffsetsPerElement_ --> mesh.typeOffsetsPerElement

    @fact length(mesh.verts) --> mesh.numVert
    @fact length(mesh.edges) --> mesh.numEdge
    @fact length(mesh.elements) --> mesh.numEl

    
    @fact mesh.jac --> roughly(ones(mesh.numNodesPerElement ,2))


    for name in fnames
      name_code = string(name, ".dat")
      name_ref = string(name, "_p", order, "true.dat")
      println("checking file ", name_code)
      data_code = readdlm(name_code)
      data_ref = readdlm(name_ref)


      for i=1:length(data_code)
        data_i = data_code[i]
        println("data_i = ", data_i)
        println("data_ref[i] = ", data_ref[i])
        if typeof(data_i) <: Number
          @fact data_i --> roughly(data_ref[i], atol=1e-13)
        else
          @fact data_i --> data_ref[i]
        end
      end
    end


  end  # end loop over p=1:4


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
    "write_dofs" => false
    )



  for order = 1:4
    println("testing order ", order, " mesh")
    sbp = TriSBP{Float64}(degree=order)

    mesh =  PumiMesh2{Float64}(dmg_name, smb_name, order, sbp, opts, coloring_distance=2, dofpernode=4)

    
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
          @fact data_i --> roughly(data_ref[i], atol=1e-13)
        else
          @fact data_i --> data_ref[i]
        end
      end
    end
  
  end   # end loop order=1:4


end

facts("----- Testing PdePumiInterfaceDG -----") do

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
    "write_dofs" => false
    )
    order = 1
    smb_name = "tri2l.smb"
    dmg_name = ".null"
    interp_op = [0.5 0 0; 0 0.5 0; 0 0 0.5]
    sbp = TriSBP{Float64}(degree=order, reorder=false, internal=true)

    vtx = [0. 0; 1 0; 0 1]
    sbpface = TriFace{Float64}(order, sbp.cub, vtx)
    mesh =  PumiMeshDG2{Float64}(dmg_name, smb_name, order, sbp, opts, interp_op, sbpface, coloring_distance=2, dofpernode=4)

   @fact mesh.m_ptr --> not(C_NULL)
   @fact mesh.mnew_ptr --> not(C_NULL)
   @fact mesh.numVert --> 4
   @fact mesh.numEdge --> 5
   @fact mesh.numEl --> 2
   @fact mesh.numDof --> 24
   @fact mesh.numNodes --> 6
   @fact mesh.numDofPerNode --> 4
   @fact mesh.numBoundaryEdges --> 4
   @fact mesh.numInterfaces --> 1
   @fact mesh.numNodesPerElement --> 3
   @fact mesh.numNodesPerType --> [0, 0, 3]
   @fact mesh.numEntitiesPerType --> [4, 5, 2]
   @fact mesh.numTypePerElement --> [3, 3, 1]
   @fact mesh.typeOffsetsPerElement --> [1, 1, 1, 4]
   @fact mesh.typeOffsetsPerElement_ --> mesh.typeOffsetsPerElement
   @fact mesh.dim --> 2
   @fact mesh.isDG --> true
   @fact mesh.coloringDistance --> 2
   @fact mesh.numColors --> 4
   @fact mesh.numBC --> 1
   @fact mesh.elementNodeOffsets --> zeros(mesh.numNodesPerElement, mesh.numEl)
   @fact mesh.typeNodeFlags[1] --> trues(3, mesh.numEl)
   tmp = trues(3, 2); tmp[1, 1] = false
   @fact mesh.typeNodeFlags[2] --> tmp
   @fact mesh.bndry_offsets --> [1, 5]
   iface = mesh.interfaces[1]
   @fact iface.elementL --> 1
   @fact iface.elementR --> 2
   @fact iface.faceL --> 1
   @fact iface.faceR --> 3
   @fact mesh.coords[:, :, 1] --> roughly([-2/3 -2/3 1/3; 2/3 -1/3 2/3], atol=1e-13)
   @fact mesh.coords[:, :, 2] --> roughly([2/3 -1/3 2/3; 1/3 -2/3 -2/3], atol=1e-13)
   @fact sort(unique(mesh.dofs)) --> collect(1:mesh.numDof)
   for i=1:size(mesh.sparsity_bnds, 2)
     @fact mesh.sparsity_bnds[1,i] --> 1
     @fact mesh.sparsity_bnds[2,i] --> 24
   end

   for i=1:size((mesh.sparsity_nodebnds), 2)
     @fact mesh.sparsity_nodebnds[1,i] --> 1
     @fact mesh.sparsity_nodebnds[2,i] --> 6
   end

   @fact mesh.color_masks[1][1] --> true
   @fact mesh.color_masks[1][2] --> false
   @fact mesh.color_masks[2][1] --> false
   @fact mesh.color_masks[2][2] --> true
   @fact mesh.neighbor_colors[1, 1] --> 2
   @fact mesh.neighbor_colors[2, 1] --> 1
   @fact mesh.neighbor_colors[1, 2] --> 1
   @fact mesh.neighbor_colors[2, 2] --> 2

   @fact mesh.neighbor_nums[1, 1] --> 2
   @fact mesh.neighbor_nums[2, 1] --> 1
   @fact mesh.neighbor_nums[1, 2] --> 1
   @fact mesh.neighbor_nums[2, 2] --> 2

   @fact mesh.pertNeighborEls[1, 1] --> 1
   @fact mesh.pertNeighborEls[2, 1] --> 1
   @fact mesh.pertNeighborEls[1, 2] --> 2
   @fact mesh.pertNeighborEls[2, 2] --> 2

   @fact mesh.color_cnt[1] --> 1
   @fact mesh.color_cnt[2] --> 1

 
   function test_interp{Tmsh}(mesh::AbstractMesh{Tmsh})
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
#       println("bndry = ", bndry_arr[1])
       dxdxi_element_rshape = reshape(dxdxi_element, 4, mesh.numNodesPerElement, 1)
#       println("dxdxi_element_rshape = \n", dxdxi_element_rshape)
#       println("dxdxi_face = \n", dxdxi_face)
#       println("sbpface.numnodes = ", sbpface.numnodes)
#       println("sbpface.stencilsize = ", sbpface.stencilsize)
#       println("size(dxdxi_element_rshape) = ", size(dxdxi_element_rshape))
#       println("size(dxdxi_face) = ", size(dxdxi_face))
#       println("sbpface = ", sbpface)
#       println("about to call boundaryinterpolate")
       boundaryinterpolate!(sbpface, bndry_arr, dxdxi_element_rshape, dxdxi_face)
#       println("finished calling boundary interpolate")

#       println("dxdxi_face = \n", dxdxi_face)
       dxdxi_face_rshape = reshape(dxdxi_face, 2, 2, sbpface.numnodes)
#       println("dxidx_face_rshape = \n", dxdxi_face_rshape)
    for j=1:sbpface.numnodes
      dxidx = inv(dxdxi_face_rshape[:, :, j])
      jac_face[j] = det(dxidx)
      dxidx_face[:, :, j] = dxidx/det(dxidx)

      @fact jac_face[j] --> roughly(mesh.jac_face[j, i], atol=1e-13)
      @fact dxidx_face[:, :, j] --> roughly(mesh.dxidx_face[:, :, j, i], atol=1e-13)
    end

    end  # end loop over interfaces

  end  # end function


   test_interp(mesh)


end
