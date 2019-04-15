# wrappers around apf::Field and apf::Numbering creation that deal with field
# naming and recording of which fields are attached to which Julia mesh
# objects (even though they are all attached to the same apf::Mesh object)

#------------------------------------------------------------------------------
# Internal reference counting functions
mutable struct ApfRef
  m::Ptr{Void}  # the mesh
  count::Int
end

global const FieldRefs = Dict{Ptr{Void}, ApfRef}()
global const NumberingRefs = Dict{Ptr{Void}, ApfRef}()

function pushRef(f_ptr::Ptr{Void}, m_ptr::Ptr{Void}, Refs::Dict)

  if haskey(Refs, f_ptr)
    @assert Refs[f_ptr].m == m_ptr
    Refs[f_ptr].count += 1
  else
    Refs[f_ptr] = ApfRef(m_ptr, 1)
  end

  apf.pushMeshRef(m_ptr)  # need to keep the mesh alive so the mesh doesn't try
                          # to free the field

  return nothing
end


function popRef(f_ptr::Ptr{Void}, Refs::Dict, destructor::Function)

  @assert haskey(Refs, f_ptr)

  obj = Refs[f_ptr]
  obj.count -= 1

  if obj.count == 0
    destructor(f_ptr)
    delete!(Refs, f_ptr)
  end

  apf.popMeshRef(obj.m)

  return nothing
end


#------------------------------------------------------------------------------
# Front ends for Fields and Numberings

function pushFieldRef(f_ptr::Ptr{Void}, m_ptr::Ptr{Void})

  pushRef(f_ptr, m_ptr, FieldRefs)

  return nothing
end

function popFieldRef(f_ptr::Ptr{Void})

  popRef(f_ptr, FieldRefs, apf.destroyField)

  return nothing
end

function countFieldRef(f_ptr::Ptr{Void})

  if haskey(FieldRefs, f_ptr)
    return FieldRefs[f_ptr].count
  else
    return 0
  end
end


function pushNumberingRef(f_ptr::Ptr{Void}, m_ptr::Ptr{Void})

  pushRef(f_ptr, m_ptr, NumberingRefs)

  return nothing
end

function popNumberingRef(f_ptr::Ptr{Void})

  popRef(f_ptr, NumberingRefs, apf.destroyNumbering)

  return nothing
end

function countNumberingRef(f_ptr::Ptr{Void})

  if haskey(NumberingRefs, f_ptr)
    return NumberingRefs[f_ptr].count
  else
    return 0
  end
end



#------------------------------------------------------------------------------
# User-visible API
# Extend the apf API replacting the mesh pointer argumetn with a PumiMesh
import apf.createPackedField

"""
  Wrapper around apf.createPackedFields that does the required bookkeeping
  (reference counting and mesh association)

  **Inputs**

   * mesh: PumiMesh object
   * fieldname: the requested name for the field.  If a field with this name
                already exists, a suffix of "_i`, where i is an integer
                will be added to generate a unique field name.
   * ncomponents: number of components on each node
   * fshape: an apf::FieldShape*

  **Outputs**
  
   * f_ptr: the newly created apf::Field*
"""
function createPackedField(mesh::PumiMesh, fieldname::AbstractString,
                           ncomponents::Integer, fshape::Ptr{Void}=C_NULL)


  # createPackedField does weird things if the name already exists, so find
  # a name that doesn't
  f_ptr = apf.findField(mesh.m_ptr, fieldname)
  if f_ptr != C_NULL
    idx = 1
    while true  # this scares me
      fname2 = string(fieldname, "_", idx)
      f_ptr = apf.findField(mesh.m_ptr, fname2)
      if f_ptr == C_NULL
        fieldname = fname2
        break
      end
      idx += 1
    end
  end

  f_ptr = apf.createPackedField(mesh.m_ptr, fieldname, ncomponents, fshape)

  # add to reference count and list of associated fields
  pushFieldRef(f_ptr, mesh.m_ptr)

  push!(mesh.fields.user, f_ptr)

  return f_ptr
end

"""
  Takes a field that exists on one Julia mesh object and attach it to
  another Julia mesh object as an original field).
  The field will then be attached to both mesh objects at the same time.

  **Inputs**

   * mesh: mesh to attach to
   * f_ptr: the apf::Field*
"""
function attachOrigField(mesh::PumiMesh, f_ptr::Ptr{Void})

  if f_ptr in mesh.fields.orig
    error("cannot attach a field that is already associated with the mesh")
  end

  @assert mesh.m_ptr == apf.getFieldMesh(f_ptr)

  pushFieldRef(f_ptr, mesh.m_ptr)
  push!(mesh.fields.orig, f_ptr)
end


"""
  Like `attachOrigField`, but attaches the field as a user field.
"""
function attachUserField(mesh::PumiMesh, f_ptr::Ptr{Void})

  if f_ptr in mesh.fields.user
    error("cannot attach a field that is already associated with the mesh")
  end

  @assert mesh.m_ptr == apf.getFieldMesh(f_ptr)

  pushFieldRef(f_ptr, mesh.m_ptr)
  push!(mesh.fields.user, f_ptr)
end


"""
  Attach a vector of fields as original fields to a mesh

  **Inputs**

   * mesh: the mesh to attach the fields to
   * fields: AbstractVector of fields.
"""
function attachOrigFields(mesh::PumiMesh, fields::AbstractVector{Ptr{Void}})

  for f in fields
    attachOrigField(mesh, f)
  end

  return nothing
end


import apf.destroyField

function destroyField(mesh::PumiMesh, f_ptr::Ptr{Void})

  if f_ptr == C_NULL
    return nothing
  end

  idx = findfirst(mesh.fields.user, f_ptr)
  arr = mesh.fields.user
  if idx == 0
    idx = findfirst(mesh.fields.orig, f_ptr)
    arr = mesh.fields.orig
  end
  @assert idx != 0

  arr[idx] = C_NULL


  popFieldRef(f_ptr)

  if length(arr) > 50
    shrinkVector(arr, C_NULL)
  end
  
  return nothing
end




#------------------------------------------------------------------------------
# Numberings

import apf: createNumberingJ, destroyNumbering

function createNumberingJ(mesh::PumiMesh, fieldname::AbstractString,
                           ncomponents::Integer, fshape::Ptr{Void})

  f_ptr = apf.findNumbering(mesh.m_ptr, fieldname)
  if f_ptr != C_NULL
    idx = 1
    while true  # this scares me
      fname2 = string(fieldname, "_", idx)
      f_ptr = apf.findNumbering(mesh.m_ptr, fname2)
      if f_ptr == C_NULL
        fieldname = fname2
        break
      end
      idx += 1
    end
  end

  f_ptr = apf.createNumberingJ(mesh.m_ptr, fieldname, fshape, ncomponents)

  # add to reference count and list of associated fields
  pushNumberingRef(f_ptr, mesh.m_ptr)

  push!(mesh.numberings.user, f_ptr)

  return f_ptr
end
 

"""
"""
function attachOrigNumbering(mesh::PumiMesh, f_ptr::Ptr{Void})

  if f_ptr in mesh.numberings.orig
    error("cannot attach a field that is already associated with the mesh")
  end

  @assert mesh.m_ptr == apf.getNumberingMesh(f_ptr)

  pushNumberingRef(f_ptr, mesh.m_ptr)
  push!(mesh.numberings.orig, f_ptr)
end


"""
"""
function attachUserNumbering(mesh::PumiMesh, f_ptr::Ptr{Void})

  if f_ptr in mesh.numberings.user
    error("cannot attach a field that is already associated with the mesh")
  end

  @assert mesh.m_ptr == apf.getNumberingMesh(f_ptr)

  pushNumberingRef(f_ptr, mesh.m_ptr)
  push!(mesh.numberings.user, f_ptr)
end


"""
"""
function attachOrigNumberings(mesh::PumiMesh, numberings::AbstractVector{Ptr{Void}})

  for f in numberings
    attachOrigNumbering(mesh, f)
  end

  return nothing
end


import apf.destroyNumbering

function destroyNumbering(mesh::PumiMesh, f_ptr::Ptr{Void})

  if f_ptr == C_NULL
    return nothing
  end

  idx = findfirst(mesh.numberings.user, f_ptr)
  arr = mesh.numberings.user
  if idx == 0
    idx = findfirst(mesh.numberings.orig, f_ptr)
    arr = mesh.numberings.orig
  end
  @assert idx != 0

  arr[idx] = C_NULL

  popNumberingRef(f_ptr)

  if length(arr) > 50
    shrinkVector(arr, C_NULL)
  end
  
  return nothing
end



#------------------------------------------------------------------------------
# Utilities

"""
  Shrinks a vector by removing all entries that have the given value

  **Inputs**

   * A: the vector
   * val: the value to remove
"""
function shrinkVector(A::AbstractVector, val)

  # indices to keep
  idxs = Vector{Int}(0)
  for i=1:length(A)
    if A[i] != val
      push!(idxs, i)
    end
  end

  for i=1:length(idxs)
    A[i] = A[idx]
  end

  resize!(A, length(idxs))

  return nothing
end

"""
  Adds reference counts for all fields currently attached to the mesh as
  `orig` fields.

  **Inputs**

   * mesh: a PumiMesh object
"""
function recordAllFields(mesh::PumiMesh)

  m_ptr = mesh.m_ptr
  nfields = apf.countFields(m_ptr)
  for i=0:(nfields-1)
    f_ptr = apf.getField(m_ptr, i)
    pushFieldRef(f_ptr, m_ptr)
    push!(mesh.fields.orig, f_ptr)
  end

  return nothing
end

function recordAllNumberings(mesh::PumiMesh)

  m_ptr = mesh.m_ptr
  nfields = apf.countNumberings(m_ptr)
  for i=0:(nfields-1)
    f_ptr = apf.getNumbering(m_ptr, i)
    attachOrigNumbering(mesh, f_ptr)
    #pushNumberingRef(f_ptr, m_ptr)
    #push!(mesh.numberings.orig, f_ptr)
  end

  return nothing
end


import Base.collect

function collect(obj::AttachedData)

  f_arr = Vector{Ptr{Void}}(0)
  for val in obj.orig
    if val != C_NULL
      push!(f_arr, val)
    end
  end

  for val in obj.user
    if val !=  C_NULL
      push!(f_arr, val)
    end
  end

  return f_arr
end
