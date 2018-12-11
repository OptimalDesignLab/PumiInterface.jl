module gmi_types

using apf_types
import apf_types.ModelEntity  # == gmi_ent*

export ModelEntity, Model, Gmi_set

"""
  gmi_model*
"""
struct Model
  p::Ptr{Void}
end

"""
  Gmi_set*.  Has an extra field `v` that contains the ModelEntities themselves.
"""
struct Gmi_set
  p::Ptr{Void}
  v::Vector{ModelEntity}
end
#TODO: add finalizer?
#TODO implement AbstractAray interface?

end  # end module
