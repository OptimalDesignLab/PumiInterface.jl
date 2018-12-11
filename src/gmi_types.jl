module gmi_types

using apf_types
import apf_types.ModelEntity  # == gmi_ent*

export ModelEntity, NullModelEntity, Model, Gmi_set, _Gmi_iter, Gmi_iter

const NullModelEntity = ModelEntity(C_NULL)

"""
  gmi_model*
"""
struct Model
  p::Ptr{Void}
end

"""
  Gmi_set*.  Has an extra field `v` that contains the ModelEntities themselves.

  Users can call `finalize` on this object to free it if needed.
"""
struct Gmi_set
  p::Ptr{Void}
  v::Vector{ModelEntity}
end
#TODO implement AbstractAray interface?


"""
  gmi_iter*
"""
struct _Gmi_iter
  p::Ptr{Void}
end

"""
  A Julia iterator that behaves like `_Gmi_iter` does in C

  This type has to be mutable because of the finalizer
"""
mutable struct Gmi_iter
  it::_Gmi_iter
  g::Model
end


end  # end module
