# types for apf interface

module apf_types

using PumiConfig

export pumi_libname, IsoFuncJ, SolutionTransfers, MAInput, ModelEntity,
       _MeshIterator, MeshIterator, SubMeshData

global const pumi_libname = joinpath(CONFIG_PATHS["PUMIINTERFACE_LIBDIR"], "libpumiInterface")

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
  apf::MeshIterator* object itself.  Prefer [`MeshIterator`](@ref).
"""
struct _MeshIterator
  p::Ptr{Void}
end


"""
  Immutable wrapper for a MeshIterator of any dimension entity.

  This type also supports Julia's iteration protocol, ie.
  ```
    for e in MeshIterator(m_ptr, dim)
      # do stuff
    end
  ```
"""
struct MeshIterator
  p::_MeshIterator
  m_ptr::Ptr{Void}
  len::Int
end

"""
  Type to encapsulate a pointer to a SubMeshData class
"""
struct SubMeshData
  pobj::Ptr{Void}
end



end # end module
