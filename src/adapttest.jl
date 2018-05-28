push!(LOAD_PATH, "/users/creanj/julialib_fork/PUMI.jl")
using PumiInterface

dmg_name = ".null"
smb_name = "tri8l.smb"

tmp, num_Entities, m_ptr, mshape_ptr = init2(dmg_name, smb_name, 1)

writeVtkFiles("output_pre", m_ptr)

#=
println("about to create field")
f_ptr = createPackedField(m_ptr, "testfield", 4)
println("created field")
comps = [1.0, 2, 3, 4]

resetVertIt()
vert_ptr = getVert()  # get a vertex
setComponents(f_ptr, vert_ptr, 0, comps)
println("set components")
comps_retrieved = [0.0, 0, 0, 0]
getComponents(f_ptr, vert_ptr, 0, comps_retrieved)
println("comps_retrieved = ", comps_retrieved)
=#

function func2(entity_ptr, m_ptr, u_::Ptr{T}) where T
# determine the mesh size at a vertex
# this function can have arbitrary arguments, but can only access data that is passed in as an argument (cannot violate trival closure condition)

  u = pointer_to_array(u_, 9)  # pointer_to_array(p, dims), where p is the pointer, dims can either be an integer that is the number of elements for a 1D array, or an array (tuple?) that tells the number of elements in each dimension of the array.  This enables the C code that is behind pointer_to_array to make sense of the data laid out in memory
  coords = zeros(3,1)
  getVertCoords(entity_ptr, coords, 3, 1)
  x_coords = coords[1]
  frac = (x_coords + 1.5)/2.0
  h_value = 1.25*frac
  println("u = \n", u)

  return h_value  # a Cdouble
end

function func3(entity_ptr, r_::Ptr{T}, h_::Ptr{T}, m_ptr, u_::Ptr{T}) where T
# an anisotropic function
# populates h with the desired mesh size in all three dimension

  println("entered func3")
  # load arrays from C pointers
  r = pointer_to_array(r_, (3,3))  # REMEMBER: r is transposed when it is passed back to C
  println("r = \n", r)
  h = pointer_to_array(h_, 3)
  println("h = \n", h)
  u = pointer_to_array(u_, 9)
  println("u = \n", u)
 
  r[1,1] = 1.0
  r[2,2] = 1.0
  r[3,3] = 1.0

  println("in julia, r = ", r)

  h[1] = 0.5
  h[2] = 1.0
  h[3] = 1.0


return nothing
end

u = ones(Float64, 9)

#=
cfunc2 = cfunction(func2, Cdouble, (Ptr{Void}, Ptr{Void}, Ptr{Float64}))
println("typeof(cfunc2 = ", typeof(cfunc2))
createIsoFunc(m_ptr, cfunc2, u)
runIsoAdapt(m_ptr)
=#
cfunc3 = cfunction(func3, Void, (Ptr{Void}, Ptr{Float64}, Ptr{Float64}, Ptr{Void}, Ptr{Float64}))

createAnisoFunc(m_ptr, cfunc3, u)
runAnisoAdapt(m_ptr)


writeVtkFiles("output_post", m_ptr)
tmp, num_Entities, m_ptr, mshape_ptr = init2(dmg_name, smb_name, 1, 0)  # re-initilize the same mesh

