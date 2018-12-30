# Interface to Pumi's Geometric Modelling Interface

module gmi

using PumiConfig
using apf_types
using gmi_types

# allow types to be accessed as gmi.foo
import gmi_types: ModelEntity, Model, Gmi_set, Gmi_iter


# iteration
import Base: start, next, done, eltype, iteratorsize, iteratoreltype, eltype
using Base: SizeUnknown, HasEltype


const gmi_lib = joinpath(CONFIG_PATHS["PUMI_LIBDIR"], "libgmi")

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
  gmi_end(obj.g, obj.it)
end


function Gmi_iter(g::Model, dim::Integer)
  it = _gmi_begin(g, dim)
  obj = Gmi_iter(it, g)

  return obj
end


function start(obj::Gmi_iter)
  return  _gmi_next(obj.g, obj.it)
end


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


function iteratorsize(::Type{Gmi_iter})
  return SizeUnknown()
end

function iteratoreltype(::Type{Gmi_iter})
  return HasEltype()
end

function eltype(::Type{Gmi_iter})
  return ModelEntity
end



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


function can_eval(g::Model)

  val = ccall( (:gmi_can_eval, gmi_lib), Cint, (Model,), g)
  return val != 0
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
  If `ge` is an edge, only `t0` is populated with dX/dp[1].  If `ge` is
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

  val = ccall( (:gmi_is_in_closure_of, gmi_lib), Cint, (Model, ModelEntity, ModelEntity), g, ge1, ge2)

  return val != 0
end


"""
  Note that users must manually call `destroy` on the resulting model (ie.
  no finalizer is attached).
"""
function load(fname::String)

  g = ccall( (:gmi_load, gmi_lib), Ptr{Void}, (Cstring,), fname)
  return Model(g)
end


function destroy(g::Model)

  ccall( (:gmi_destroy, gmi_lib), Void, (Model,), g)
end

#------------------------------------------------------------------------------
# gmi_sim


function sim_start()
  ccall( (:gmi_sim_startJ, pumi_libname), Void, (), )
end

function sim_stop()
  ccall( (:gmi_sim_stopJ, pumi_libname), Void, (), )
end


function register_sim()
  ccall( (:gmi_register_simJ, pumi_libname), Void, (), )
end


# the functions gmi_sim_startJ/stopJ check if Simmetrix is supported, so call
# these unconditionally
sim_start()
atexit(sim_stop)


end # end module
