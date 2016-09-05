push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
using FactCheck
using PdePumiInterface
using ODLCommonTools
using SummationByParts
using PumiInterface
include("defs.jl")

facts("----- Testing PdePumiInterface3DG -----") do
  degree = 1
  Tsbp = Float64
  sbp = TetSBP{Tsbp}(degree=degree, reorder=false, internal=true)
  ref_verts = sbp.vtx
  interp_op = SummationByParts.buildinterpolation(sbp, ref_verts.')
  face_verts = SummationByParts.SymCubatures.getfacevertexindices(sbp.cub)
  topo = ElementTopology{3}(face_verts)
  sbpface = TetFace{Tsbp}(degree, sbp.cub, ref_verts)

  dmg_name = ".null"
#  smb_name = "tet1.smb"
  smb_name = "unitcube.smb"

  opts = PdePumiInterface.get_defaults()
  opts["numBC"] = 1
  opts["BC1"] = [0,1,2,3,4,5]

  interp_op = eye(4)
  mesh = PumiMeshDG3{Float64}(dmg_name, smb_name, degree, sbp, opts, interp_op, sbpface, topo)

  println("finished")

  @fact mesh.m_ptr --> not(C_NULL)
  @fact mesh.mshape_ptr --> not(C_NULL)
  @fact mesh.coordshape_ptr --> not(mesh.mshape_ptr)
  @fact mesh.vert_Nptr --> not(C_NULL)
  @fact mesh.edge_Nptr --> not(C_NULL)
  @fact mesh.face_Nptr --> not(C_NULL)
  @fact mesh.el_Nptr --> not(C_NULL)

  function checkNumbering(nshape_ptr, dim, cnt)
    for i=1:4
      nodes = countNodesOn(nshape_ptr, simplexTypes[i])
      if i == dim
        @fact nodes --> cnt
      else
        @fact nodes --> 0
      end
    end
  end

  for i=1:4
    numbering_shape = getNumberingShape(mesh.entity_Nptrs[i])
    checkNumbering(numbering_shape, i, 1)
  end

  @fact mesh.numVert --> 8
  @fact mesh.numEdge --> 19
  @fact mesh.numFace --> 18
  @fact mesh.numEl --> 6
  @fact mesh.numBoundaryFaces --> 12
  @fact mesh.numInterfaces --> 6
  @fact mesh.numFacesPerElement --> 4
  @fact mesh.numEntitiesPerType --> [mesh.numVert, mesh.numEdge, mesh.numFace, mesh.numEl]
  @fact mesh.numTypePerElement --> [4, 6, 4, 1]
  @fact mesh.el_type --> apfTET
  @fact mesh.face_type --> apfTRIANGLE
  @fact mesh.isDG --> true
  @fact mesh.dim --> 3
  
  
  for i=1:min(length(mesh.faces), length(mesh.edges))
    @fact mesh.edges[i] --> not(mesh.faces[i])
  end

  nodeshape = getNumberingShape(mesh.nodenums_Nptr)
  checkNumbering(nodeshape, 4, 4)
  dofshape = getNumberingShape(mesh.dofnums_Nptr)
  checkNumbering(dofshape, 4, 4)

  @fact mesh.bndry_offsets[1] --> 1
  @fact mesh.bndry_offsets[2] --> 13

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      for k=1:mesh.dim
        @fact mesh.coords[k, j, i] --> greater_than(0.99)
        @fact mesh.coords[k, j, i] --> less_than(3.01)
      end
    end
  end
  dofs_sorted = reshape(mesh.dofs, mesh.numDof)
  sort!(dofs_sorted)
  @fact unique(dofs_sorted) --> dofs_sorted

  @fact maximum(mesh.dofs) --> mesh.numEl*mesh.numNodesPerElement*mesh.numDofPerNode
  @fact minimum(mesh.dofs) --> 1

  for i=1:mesh.numEl
    @fact mesh.sparsity_counts[1, i] --> greater_than(0)
    @fact mesh.sparsity_counts[2, i] --> less_than(5)
  end

  @fact size(mesh.neighbor_colors, 1) --> 5
  @fact size(mesh.neighbor_nums, 1) --> 5


  # check interface array
  for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]
    @fact iface_i.elementL --> greater_than(0)
    @fact iface_i.elementL --> less_than(mesh.numEl+1)
    @fact iface_i.elementR --> greater_than(0)
    @fact iface_i.elementR --> less_than(mesh.numEl+1)
    @fact iface_i.faceL --> greater_than(0)
    @fact iface_i.faceL --> less_than(5)
    @fact iface_i.faceR --> greater_than(0)
    @fact iface_i.faceR --> less_than(5)
    @fact iface_i.orient --> greater_than(0)
    @fact iface_i.orient --> less_than(4)
  end

  # check boundary array
  for i=1:mesh.numBoundaryFaces
    bndry_i = mesh.bndryfaces[i]
    @fact bndry_i.element --> greater_than(0)
    @fact bndry_i.element --> less_than(mesh.numEl + 1)
    @fact bndry_i.face --> greater_than(0)
    @fact bndry_i.face --> less_than(5)
  end

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
          @fact dxidx_face[p, k] --> roughly(dxidx_el[p, k], atol=1e-13)
        end
      end
      @fact jac_face[j] --> roughly(jac_el[j], atol=1e-13)
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
          @fact dxidx_face[p, k] --> roughly(dxidx_el[p, k], atol=1e-13)
        end
      end
      @fact jac_face[j] --> roughly(jac_el[j], atol=1e-13)
    end
  end  # end loop over interfaces

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

  commit_coords(mesh)

  PdePumiInterface.getCoordinates(mesh, sbp)
  for i = 1:length(mesh.vert_coords)
    @fact abs(2*coords_orig[i] - mesh.vert_coords[i]) --> less_than(1e-10)
  end

  # test periodic 
  @fact mesh.numPeriodicInterfaces --> 0

  smb_name = "tet3_pxz.smb"
  opts["BC1"] = [0,2,4,5,6]
  mesh = PumiMeshDG3{Float64}(dmg_name, smb_name, degree, sbp, opts, interp_op, sbpface, topo)

  @fact mesh.numPeriodicInterfaces --> 18
  for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]
    @fact iface_i.elementL --> greater_than(0)
    @fact iface_i.elementL --> less_than(mesh.numEl + 1)
    @fact iface_i.elementR --> greater_than(0)
    @fact iface_i.elementR --> less_than(mesh.numEl + 1)
    @fact iface_i.faceL --> greater_than(0)
    @fact iface_i.faceL --> less_than(5)
    @fact iface_i.faceR --> greater_than(0)
    @fact iface_i.faceR --> less_than(5)
  end

  # this mesh test that matched vertices are used for orientation 
  # determiniation
  smb_name = "tet5_periodic.smb"
  opts["numBC"] = 0
  delete!(opts, "BC1")

  mesh = PumiMeshDG3{Float64}(dmg_name, smb_name, degree, sbp, opts, interp_op, sbpface, topo)

  @fact mesh.numPeriodicInterfaces --> 3*(5*5*2)
  for i=1:mesh.numInterfaces
    @fact mesh.interfaces[i].orient --> not(0)
    @fact mesh.interfaces[i].orient --> less_than(4)
  end






end



