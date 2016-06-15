push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
using FactCheck
using PdePumiInterface
using ODLCommonTools
using SummationByParts
using PumiInterface

immutable MySBP{Tsbp} <: AbstractSBP{Tsbp}
  degree::Int
  numnodes::Int
  numfacenodes::Int
end

immutable MyFace{Tsbp} <: SummationByParts.AbstractFace{Tsbp}
  degree::Int
  numnodes::Int
  stencilsize::Int
end

facts("----- Testing PdePumiInterface3DG -----") do
  degree = 1
  numnodes = 4
  Tsbp = Float64
  sbp = MySBP{Tsbp}(degree, numnodes, 0)
  sbpface = MyFace{Tsbp}(degree, 3, numnodes)
  topo = ElementTopology3()

  dmg_name = ".null"
  smb_name = "tet1.smb"

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

  function checkNumbering(nshape_ptr, dim)
    for i=1:4
      nodes = countNodesOn(nshape_ptr, simplexTypes[i])
      if i == dim
        @fact nodes --> 1
      else
        @fact nodes --> 0
      end
    end
  end

  for i=1:4
    numbering_shape = getNumberingShape(mesh.entity_Nptrs[i])
    checkNumbering(numbering_shape, i)
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
  @fact mesh.face_type--> apfTRIANGLE
  @fact mesh.dim --> 3
  
  
  for i=1:min(length(mesh.faces), length(mesh.edges))
    @fact mesh.edges[i] --> not(mesh.faces[i])
  end



end



