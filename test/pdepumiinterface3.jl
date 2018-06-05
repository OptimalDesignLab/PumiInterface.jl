push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
using Base.Test
using PdePumiInterface
using ODLCommonTools
using SummationByParts
using PumiInterface
include("defs.jl")

@testset "----- Testing PdePumiInterface3DG -----" begin
  degree = 1
  Tsbp = Float64
  sbp = getTetSBPOmega(degree=degree, Tsbp=Tsbp)
#  sbp = TetSBP{Tsbp}(degree=degree, internal=true)
  ref_verts = sbp.vtx
  interp_op = SummationByParts.buildinterpolation(sbp, ref_verts.')
  face_verts = SummationByParts.SymCubatures.getfacevertexindices(sbp.cub)
  topo = ElementTopology{3}(face_verts)
  sbpface = TetFace{Tsbp}(degree, sbp.cub, ref_verts)

  dmg_name = ".null"
#  smb_name = "tet1.smb"
  smb_name = "unitcube.smb"

  opts = PdePumiInterface.get_defaults()
  opts["smb_name"] = smb_name
  opts["dmg_name"] = dmg_name
  opts["order"] = degree
  opts["coloring_distance"] = 2
  opts["use_linear_metrics"] = true
  opts["numBC"] = 1
  opts["BC1"] = [0,1,2,3,4,5]

#  interp_op = eye(4)
  mesh = PumiMeshDG3(Float64, sbp, opts, sbpface, topo)
#  mesh = PumiMeshDG3{Float64, typeof(sbpface)}(dmg_name, smb_name, degree, sbp, opts, sbpface, topo)


  @test ( mesh.m_ptr )!=C_NULL
  @test ( mesh.mshape_ptr )!=C_NULL
  @test ( mesh.coordshape_ptr )!=mesh.mshape_ptr
  @test ( mesh.vert_Nptr )!=C_NULL
  @test ( mesh.edge_Nptr )!=C_NULL
  @test ( mesh.face_Nptr )!=C_NULL
  @test ( mesh.el_Nptr )!=C_NULL
  @test ( mesh.numBC )== 1
  @test ( mesh.bndry_geo_nums[1] )== opts["BC1"]

  function checkNumbering(nshape_ptr, dim, cnt)
    for i=1:4
      nodes = countNodesOn(nshape_ptr, simplexTypes[i])
      if i == dim
        @test ( nodes )== cnt
      else
        @test ( nodes )== 0
      end
    end
  end

  for i=1:4
    numbering_shape = getNumberingShape(mesh.entity_Nptrs[i])
    checkNumbering(numbering_shape, i, 1)
  end

  @test ( mesh.numVert )== 8
  @test ( mesh.numEdge )== 19
  @test ( mesh.numFace )== 18
  @test ( mesh.numEl )== 6
  @test ( mesh.numBoundaryFaces )== 12
  @test ( mesh.numInterfaces )== 6
  @test ( mesh.numFacesPerElement )== 4
  @test isapprox( mesh.volume, 1.0) atol=1e-12
  @test ( mesh.numEntitiesPerType )== [mesh.numVert, mesh.numEdge, mesh.numFace, mesh.numEl]
  @test ( mesh.numTypePerElement )== [4, 6, 4, 1]
  @test ( mesh.el_type )== apfTET
  @test ( mesh.face_type )== apfTRIANGLE
  @test ( mesh.isDG )== true
  @test ( mesh.dim )== 3
  
  
  for i=1:min(length(mesh.faces), length(mesh.edges))
    @test ( mesh.edges[i] )!=mesh.faces[i]
  end

  nodeshape = getNumberingShape(mesh.nodenums_Nptr)
  checkNumbering(nodeshape, 4, 4)
  dofshape = getNumberingShape(mesh.dofnums_Nptr)
  checkNumbering(dofshape, 4, 4)

  @test ( mesh.bndry_offsets[1] )== 1
  @test ( mesh.bndry_offsets[2] )== 13

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      for k=1:mesh.dim
        @test  mesh.coords[k, j, i]  > 0.99
        @test  mesh.coords[k, j, i]  < 3.01
      end
    end
  end
  dofs_sorted = reshape(mesh.dofs, mesh.numDof)
  sort!(dofs_sorted)
  @test ( unique(dofs_sorted) )== dofs_sorted

  @test ( maximum(mesh.dofs) )== mesh.numEl*mesh.numNodesPerElement*mesh.numDofPerNode
  @test ( minimum(mesh.dofs) )== 1

  for i=1:mesh.numEl
    @test  mesh.sparsity_counts[1, i]  > 0
    @test  mesh.sparsity_counts[2, i]  < 5
  end

  @test ( size(mesh.neighbor_colors, 1) )== 5
  @test ( size(mesh.neighbor_nums, 1) )== 5


  # check interface array
  for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]
    @test  iface_i.elementL  > 0
    @test  iface_i.elementL  < mesh.numEl+1
    @test  iface_i.elementR  > 0
    @test  iface_i.elementR  < mesh.numEl+1
    @test  iface_i.faceL  > 0
    @test  iface_i.faceL  < 5
    @test  iface_i.faceR  > 0
    @test  iface_i.faceR  < 5
    @test  iface_i.orient  > 0
    @test  iface_i.orient  < 4
  end

  # check boundary array
  for i=1:mesh.numBoundaryFaces
    bndry_i = mesh.bndryfaces[i]
    @test  bndry_i.element  > 0
    @test  bndry_i.element  < mesh.numEl + 1
    @test  bndry_i.face  > 0
    @test  bndry_i.face  < 5
  end

  testSurfaceNumbering(mesh, sbp, opts)

  # check mapping interpolation
  # should be constant within an element for straight-sided elements
  for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]
    el_i = iface_i.elementL
    dxidx_el = mesh.dxidx[:, :, 1, el_i]
    jac_el = mesh.jac[:, el_i]
    jac_face = mesh.jac_face[:, i]

    for j=1:mesh.numNodesPerFace
      dxidx_face = mesh.dxidx_face[:, :, j, i]
      for k=1:3
        for p=1:3
          @test isapprox( dxidx_face[p, k], dxidx_el[p, k]) atol=1e-13
        end
      end
      @test isapprox( jac_face[j], jac_el[j]) atol=1e-13
    end
  end  # end loop over interfaces

  for i=1:mesh.numBoundaryFaces
    bndry_i = mesh.bndryfaces[i]
    el_i = bndry_i.element
    dxidx_el = mesh.dxidx[:, :, 1, el_i]
    jac_el = mesh.jac[:, el_i]
    jac_face = mesh.jac_bndry[:, i]

    for j=1:mesh.numNodesPerFace
      dxidx_face = mesh.dxidx_bndry[:, :, j, i]
      for k=1:3
        for p=1:3
          @test isapprox( dxidx_face[p, k], dxidx_el[p, k]) atol=1e-13
        end
      end
      @test isapprox( jac_face[j], jac_el[j]) atol=1e-13
    end
  end  # end loop over interfaces


  # the 2D tests should work for this too
  test_interp_rev(mesh)

  # test update_coords
  coords_orig = copy(mesh.vert_coords)
  
  # double all coordinates
  coords_i = zeros(Float64, 3, 4)
  for i=1:mesh.numEl
    for j=1:4
      coords_i[1, j] = 2*coords_orig[1, j, i]
      coords_i[2, j] = 2*coords_orig[2, j, i]
      coords_i[3, j] = 2*coords_orig[3, j, i]
    end
    update_coords(mesh, i, coords_i)
  end

  commit_coords(mesh, sbp, opts)

#  PdePumiInterface.getCoordinatesAndMetrics(mesh, sbp)
  for i = 1:length(mesh.vert_coords)
    @test  abs(2*coords_orig[i] - mesh.vert_coords[i])  < 1e-10
  end

  # test vertmap

  nfaces_interior = 0
  for i=1:mesh.numEl
    vertnums_i = mesh.element_vertnums[:, i]
    for j=1:4
      @test  vertnums_i[j]  > 0
      @test  vertnums_i[j]  < mesh.numVert + 1
    end

    # test there are exactly 3 elements the current element shares 2 vertices
    # with
    nfaces = 0
    for j=1:mesh.numEl
      vertnums_j = mesh.element_vertnums[:, j]
      if length(intersect(vertnums_i, vertnums_j)) == 3
        nfaces += 1
        nfaces_interior += 1
      end
    end

    @test  nfaces  < 4
    @test  nfaces  > 0
  end  # end loop i

  @test ( nfaces_interior )== 2*(mesh.numFace - mesh.numBoundaryFaces)


  # test periodic 
  @test ( mesh.numPeriodicInterfaces )== 0

  opts["smb_name"] = "tet3_pxz.smb"
  opts["BC1"] = [0,2,4,5,6]
  mesh = PumiMeshDG3(Float64, sbp, opts, sbpface, topo)
#  mesh = PumiMeshDG3{Float64, typeof(sbpface)}(dmg_name, smb_name, degree, sbp, opts, sbpface, topo)

  @test ( mesh.numPeriodicInterfaces )== 18
  for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]
    @test  iface_i.elementL  > 0
    @test  iface_i.elementL  < mesh.numEl + 1
    @test  iface_i.elementR  > 0
    @test  iface_i.elementR  < mesh.numEl + 1
    @test  iface_i.faceL  > 0
    @test  iface_i.faceL  < 5
    @test  iface_i.faceR  > 0
    @test  iface_i.faceR  < 5
  end

  # this mesh test that matched vertices are used for orientation 
  # determiniation

  opts["smb_name"] = "tet5_periodic.smb"
  opts["numBC"] = 0
  delete!(opts, "BC1")

  mesh = PumiMeshDG3(Float64, sbp, opts, sbpface, topo)
#  mesh = PumiMeshDG3{Float64, typeof(sbpface)}(dmg_name, smb_name, degree, sbp, opts, sbpface, topo)

  @test ( mesh.numPeriodicInterfaces )== 3*(5*5*2)
  for i=1:mesh.numInterfaces
    @test ( mesh.interfaces[i].orient )!=0
    @test  mesh.interfaces[i].orient  < 4
  end


  for order=1:2
    fshape = getFieldShape(0, order, 3)
    eshape = getEntityShape(fshape, apfTET)  # triangle
    nodexi = PdePumiInterface.getXiCoords(order, 3)
    numnodes = size(nodexi, 2)
    for i=1:numnodes
      vals = getValues(mesh.m_ptr, eshape, nodexi[:, i], numnodes)
      for j=1:numnodes
        if i == j
          @test  abs(vals[j] - 1)  < 1e-12
        else
          @test  abs(vals[j])  < 1e-12
        end  # end if else
      end   # end loop j
    end  # end loop i

  end  # end loop order

  mesh = PumiMeshDG3(Float64, sbp, opts, sbpface, topo)
#  mesh = PumiMeshDG3{Float64, typeof(sbpface)}(dmg_name, smb_name, degree, sbp, opts, sbpface, topo)


  mesh_c = PumiMeshDG3(Complex128, sbp, opts, sbpface, topo)
#  mesh_c = PumiMeshDG3{Complex128, typeof(sbpface)}(dmg_name, smb_name, degree, sbp, opts, sbpface, topo)

  compare_meshes(mesh, mesh_c)

  # test curvilinear is same as linear for linear mesh
  println("testing curvilinear")
  degree = 2
  Tsbp = Float64
  sbp = getTetSBPOmega(degree=degree, Tsbp=Tsbp)
#  sbp = TetSBP{Tsbp}(degree=degree, internal=true)
#  ref_verts = sbp.vtx
  interp_op = SummationByParts.buildinterpolation(sbp, ref_verts.')
  face_verts = SummationByParts.SymCubatures.getfacevertexindices(sbp.cub)
  topo = ElementTopology{3}(face_verts)
  sbpface = TetFace{Tsbp}(degree, sbp.cub, ref_verts)


  opts = PdePumiInterface.get_defaults()
  opts["dmg_name"] = ".null"
#  smb_name = "tet1.smb"
  opts["smb_name"] = "unitcube.smb"
  opts["coloring_distance"] = 2


  opts["use_linear_metrics"] = true
  opts["order"] = degree
  opts["numBC"] = 1
  opts["BC1"] = [0,1,2,3,4,5]

#  interp_op = eye(4)
  mesh_c = PumiMeshDG3(Complex128, sbp, opts, sbpface, topo)
#  mesh_c = PumiMeshDG3{Complex128, typeof(sbpface)}(dmg_name, smb_name, degree, sbp, opts, sbpface, topo)
  mesh = PumiMeshDG3(Float64, sbp, opts, sbpface, topo)
#  mesh = PumiMeshDG3{Float64, typeof(sbpface)}(dmg_name, smb_name, degree, sbp, opts, sbpface, topo)

  
  dxidx_orig = copy(mesh.dxidx)
  jac_orig = copy(mesh.jac)
  nrm_bndry_orig = copy(mesh.nrm_bndry)
  nrm_face_orig = copy(mesh.nrm_face)
  coords_orig = copy(mesh.coords)
  coords_bndry_orig = copy(mesh.coords_bndry)
  coords_interface_orig = copy(mesh.coords_interface)
  
  # use curvilinear routines to do calculation
  PdePumiInterface.getMeshCoordinates(mesh, sbp)
  PdePumiInterface.getFaceCoordinatesAndNormals(mesh, sbp)
 # PdePumiInterface.getCurvilinearCoordinatesAndMetrics(mesh, sbp)
  

  # verify metrics are constant for a uniform linear mesh

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      @test isapprox( norm(mesh.dxidx[:, :, j, i] - mesh.dxidx[:, :, 1, i]), 0.0) atol=1e-12
    end
  end

  for i=1:mesh.numBoundaryFaces
    for j=1:mesh.numNodesPerFace
      @test isapprox( norm(mesh.nrm_bndry[:, j, i] - mesh.nrm_bndry[:, 1, i]), 0.0) atol=1e-12
    end
  end

  for i=1:mesh.numInterfaces
    for j=1:mesh.numNodesPerFace
      @test isapprox( norm(mesh.nrm_face[:, j, i] - mesh.nrm_face[:, 1, i]), 0.0) atol=1e-12
    end
  end

  

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
        @test isapprox( mesh.coords_bndry[k, j, i], coords_bndry_orig[k, j, i]) atol=1e-12
      end
    end
  end

  for i=1:mesh.numInterfaces
    for j=1:mesh.numNodesPerFace
      for k=1:mesh.dim
        @test isapprox( mesh.nrm_face[k, j, i], nrm_face_orig[k, j, i]) atol=1e-12
        @test isapprox( mesh.coords_interface[k, j, i], coords_interface_orig[k, j, i]) atol=1e-12
      end
    end
  end

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      for d=1:mesh.dim
        @test isapprox( mesh.coords[d, j, i], coords_orig[d, j, i]) atol=1e-13
      end
    end
  end


  testSurfaceNumbering(mesh, sbp, opts)

  # test reverse mode
  PdePumiInterface.copy_data!(mesh_c, mesh)
  test_metric_rev(mesh, mesh_c, sbp, opts)


  # test with a diagonal E operator
  degree = 2
  Tsbp = Float64
  sbp = getTetSBPDiagE(degree=degree, Tsbp=Tsbp)
#  sbp = TetSBP{Tsbp}(degree=degree, internal=true)
#  ref_verts = sbp.vtx
  interp_op = SummationByParts.buildinterpolation(sbp, ref_verts.')
  face_verts = SummationByParts.SymCubatures.getfacevertexindices(sbp.cub)
  topo = ElementTopology{3}(face_verts)
  sbpface = getTetFaceForDiagE(degree, sbp.cub, ref_verts)
#  sbpface = TetFace{Tsbp}(degree, sbp.cub, ref_verts)

  opts["dmg_name"] = ".null"
#  smb_name = "tet1.smb"
  opts["smb_name"] = "unitcube.smb"

  opts["use_linear_metrics"] = false

#  interp_op = eye(4)

  mesh2_c = PumiMeshDG3(Complex128, sbp, opts, sbpface, topo; shape_type=4)
#  mesh2_c = PumiMeshDG3{Complex128, typeof(sbpface)}(dmg_name, smb_name, degree, sbp, opts, sbpface, topo, shape_type=4)
  
  mesh2 = PumiMeshDG3(Float64, sbp, opts, sbpface, topo, shape_type=4)
#  mesh2 = PumiMeshDG3{Float64, typeof(sbpface)}(dmg_name, smb_name, degree, sbp, opts, sbpface, topo, shape_type=4)

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      @test isapprox( norm(mesh.dxidx[:, :, j, i] - mesh2.dxidx[:, :, 1, i]), 0.0) atol=1e-12
    end
  end

  for i=1:mesh.numBoundaryFaces
    for j=1:mesh.numNodesPerFace
      @test isapprox( norm(mesh.nrm_bndry[:, j, i] - mesh2.nrm_bndry[:, 1, i]), 0.0) atol=1e-12
    end
  end

  for i=1:mesh.numInterfaces
    for j=1:mesh.numNodesPerFace
      @test isapprox( norm(mesh.nrm_face[:, j, i] - mesh2.nrm_face[:, 1, i]), 0.0) atol=1e-12
    end
  end

  

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      for k=1:mesh.dim
        for p=1:mesh.dim
          @test isapprox( mesh.dxidx[p, k, j, i], mesh2.dxidx[p, k, 1, i]) atol=1e-12
        end
      end

      @test isapprox( mesh.jac[j, i], mesh2.jac[j, i]) atol=1e-12
    end
  end

  for i=1:mesh.numBoundaryFaces
    for j=1:mesh.numNodesPerFace
      for k=1:mesh.dim
        @test isapprox( mesh.nrm_bndry[k, j, i], nrm_bndry_orig[k, 1, i]) atol=1e-12
      end
    end
  end

  for i=1:mesh.numInterfaces
    for j=1:mesh.numNodesPerFace
      for k=1:mesh.dim
        @test isapprox( mesh.nrm_face[k, j, i], nrm_face_orig[k, 1, i]) atol=1e-12
      end
    end
  end

  # test reverse mode
  PdePumiInterface.copy_data!(mesh2_c, mesh2)
  test_metric_rev(mesh2, mesh2_c, sbp, opts)


 

end



