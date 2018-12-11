# types for apf interface

module apf_types

export IsoFuncJ, SolutionTransfers, MAInput, ModelEntity, MeshIterator, SubMeshData

# struct declarations
# these have same memory layout as their contents, so they can be passed in
# in place of a Ptr{Void}

"""
  Isotropic mesh size function used for mesh adapation

  In reality, it is just a (typed) container for a pointer to a C++ 
  IsotropicFunctionJ.
"""
struct IsoFuncJ
  p::Ptr{Void}
end

"""
  A ma::SolutionTransfers*
"""
struct SolutionTransfers
  p::Ptr{Void}
end

"""
  A ma::Input*
"""
struct MAInput
  p::Ptr{Void}
end

"""
  An :ModelEntity* aka. gmi_ent*
"""
struct ModelEntity
  p::Ptr{Void}
end

"""
  Immutable wrapper for a MeshIterator of any dimension entity.

  This grants some type safety to the interface
"""
struct MeshIterator
  p::Ptr{Void}
end

"""
  Type to encapsulate a pointer to a SubMeshData class
"""
struct SubMeshData
  pobj::Ptr{Void}
end



end # end module
