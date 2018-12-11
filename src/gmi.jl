# Interface to Pumi's Geometric Modelling Interface

module gmi

using PumiConfig
using gmi_types

# allow types to be accessed as gmi.foo
import gmi_types: ModelEntity, Model, Gmi_set, Gmi_iter


# iteration
import Base: start, next, done, eltype, iteratorsize, eltype


const gmi_lib = joinpath(CONFIG_PATHS["PUMIINTERFACE_LIBDIR"], "libgmi")

#------------------------------------------------------------------------------
# gmi_set stuff

"""
  ccall wrapper.  Users should not call this, call the constructor for
  Gmi_set instead.
"""
function _make_set(n::Integer)

  p = ccall( (:gmi_make_set, gmi_lib), Ptr{Void}, (Cint,), n)
  return p
end

function free_set(s::Ptr{Void})

  ccall( (:gmi_free_set, gmi_lib), Void, (Ptr{Void},), s)
end

function free_set(s::Gmi_set)

  free_set(s.p)
end

function Gmi_set(n::Integer)

  p = _make_set(n)
  v = Vector{ModelEntity}(n)

  obj = Gmi_set(p, v)
  finalizer(obj, free_set)

  return obj
end



#------------------------------------------------------------------------------
# gmi_iter

function _gmi_begin(g::Model, dim::Integer)
  p = ccall( (:gmi_begin, gmi_lib), Ptr{Void}, (Model, Cint), g, dim)
  return _Gmi_iter(p)
end


function _gmi_next(g::Model, it::_Gmi_iter)
  ge = ccall( (:gmi_next, gmi_lib), Ptr{Void}, (Model, _Gmi_iter), g, it)
  return ModelEntity(ge)
end


function gmi_end(g::Model, it::_Gmi_iter)
  ccall( (:gmi_end, gmi_lib), Void, (Model, _Gmi_iter), g, it)
end

function gmi_end(obj::Gmi_iter)
  gmi_end(obj.it)
end


function Gmi_iter(g::Model, dim::Integer)
  it = _gmi_begin(g, dim)
  obj = new(it, g)

  # if not for the finalizer, this struct could be immutable
  finalizer(obj, gmi_end)
  return obj
end


function start(obj::Gmi_iter) return  _gmi_next(obj.g, obj.it) end


function next(obj::Gmi_iter, state::ModelEntity)

  # the state always leads the current element by 1, so we can do
  # the done check
  next_el = _gmi_next(obj.g, obj.it)
  return (state, next_el)
end


function done(obj::Gmi_iter, state::ModelEntity)

  isdone = state == NullModelEntity
  if isdone
    gmi_end(obj)
  end

  return state == NullModelEntity
end


function eltype(::Gmi_iter) return ModelEntity end



#------------------------------------------------------------------------------
# more methods

function dim(g::Model, ge::ModelEntity)
  return ccall((:gmi_dim, gmi_lib), Cint, (Model, ModelEntity), g, ge)
end


function tag(g::Model, ge::ModelEntity)
  return ccall((:gmi_tag, gmi_lib), Cint, (Model, ModelEntity), g, ge)
end


function find(g::Model, dim::Integer, tag::Integer)
  ge = ccall( (:gmi_find, gmi_lib), Ptr{Void}, (Model, Cint, Cint), g, dim, tag)
  return ModelEntity(ge)
end


function adjacent(g::Model, ge::ModelEntity, dim::Integer)

  # stuff these functions into libPumiInterface
  n = ccall( (:gmi_adjacent_count, pumi_libname), Cint, (Model, ModelEntity, Cint), g, ge, dim)

  v = Vector{ModelEntity}(n)
  gset = ccall( (:gmi_adjacent_get, pumi_libname), Ptr{Void}, (Ptr{ModelEntity},), v)

  obj = Gmi_set(gset, v)
  finalizer(obj, free_set)

  return obj
end


"""
  Note the prefix `g`, to avoid conflict with Julia's eval
"""
function geval(g::Model, ge::ModelEntity, p::AbstractVector{Cdouble},
               x::AbstractVector{Cdouble})

  @assert length(p) == 2
  @assert length(x) == 3
  ccall( (:gmi_eval, gmi_lib), Void, (Model, ModelEntity, Ptr{Cdouble}, Ptr{Cdouble}), g, ge, p, x)
end


function reparam(g::Model, ge_from::ModelEntity,
                 from_p::AbstractVector{Cdouble},
                 ge_to::ModelEntity, to_p::AbstractVector{Cdouble})

  ccall( (:gmi_reparam, gmi_lib), Void, (Model, ModelEntity, Ptr{Cdouble}, ModelEntity, Ptr{Cdouble}), g, ge_from,  from_p, ge_to, to_p)
end


function periodic(g::Model, ge::ModelEntity, dim::Integer)

  val = ccall( (:gmi_periodic, gmi_lib), Cint, (Model, ModelEntity, Cint), g, ge, dim)

  return val != 0
end


function range(g::Model, ge::ModelEntity, dim::Integer, r::AbstractVector{Cdouble})

  @assert length(r) == 2

  ccall( (:gmi_range, gmi_lib), Void, (Model, ModelEntity, Cint, Ptr{Cdouble}), g, ge, dim, r)
end


function closest_point(g::Model, ge::ModelEntity,
                           from::AbstractVector{Cdouble},
                           to::AbstractVector{Cdouble},
                           to_p::AbstractVector{Cdouble})

  @assert length(from) == 3
  @assert length(to) == 3
  @assert length(to_p) == 2

  ccall( (:gmi_closest_point, gmi_lib), Void, (Model, ModelEntity, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), g, ge, from, to, to_p)
end


"""
  For GeomSim models, this only works for geometric faces (not edges, even
  if the model is 2D).
"""
function normal(g::Model, ge::ModelEntity, p::AbstractVector{Cdouble},
                n::AbstractVector{Cdouble})

  ccall( (:gmi_normal, gmi_lib), Void, (Model, ModelEntity, Ptr{Cdouble}, Ptr{Cdouble}), g, ge, p, n)
end


"""
  If `ge` is an edge, only `t0` is populated with dX/dp[1].  If `get` is
  a face`, then `t1` is also populated with dX/dp[2]
"""
function first_derivative(g::Model, ge::ModelEntity, p::AbstractVector{Cdouble},
                      t0::AbstractVector{Cdouble}, t1::AbstractVector{Cdouble})

  ccall( (:gmi_first_derivative, gmi_lib), Void, (Model, ModelEntity, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), g, ge, p, t0, t1)

end


function is_point_in_region(g::Model, ge::ModelEntity,
                            point::AbstractVector{Cdouble})

  val = ccall( (:gmi_is_point_in_region, gmi_lib), Cint, (Model, ModelEntity, Ptr{Cdouble}), g, ge, point)

  return val != 0
end


function is_in_closure_of(g::Model, ge1::ModelEntity, ge2::ModelEntity)

  val = vvall( (:gmi_is_in_closure_of, gmi_lib), Cint, (Model, ModelEntity, ModelEntity), g, ge1, g2)

  return val != 0
end


end # end module
