# utilities for managing Pumi lbraries


# make MeshEntity* passable by MPI
if sizeof(Ptr{Void}) == sizeof(Int32)
  T = Int32
elseif  sizeof(Ptr{Void}) == sizeof(Int64)
  T = Int64
else
  error("cannot find compatible size of Ptr{Void}")
end

if !haskey(MPI.mpitype_dict, Ptr{Void})
  MPI.mpitype_dict[Ptr{Void}] = MPI.mpitype_dict[T]
end


# make PdePumiInterface findable
push!(LOAD_PATH, dirname(@__FILE__))

"""
  Typealias for C++ Bool.  This is implementation defined, so it might need
  to be manually changed for different systems
"""
const CppBool = UInt8  # this is implementation defined


"""
  Holds paths that were determined at configure time
"""
global const CONFIG_PATHS = Dict{String, String}()


"""
  Reads the PUMI_RPATH file and updates [`CONFIG_PATHS`](@ref)
"""
function readConfigFile()

  # this is basically like RPATH
  fname = joinpath(dirname(@__FILE__), "PUMI_RPATH")

  for line in eachline(fname)
    substr = split(chomp(line))
    @assert length(substr) == 2
    CONFIG_PATHS[substr[1]] = substr[2]
  end

  # String manipulation ins CMake is not easy, so CMake writes whatever
  # paths are convenient and this function updates them to what we really want
  updatePaths()

  return nothing
end

function updatePaths()

  # PUMI_CMAKE_PATH -> PUMI_LIBDIR
  str = CONFIG_PATHS["PUMI_CMAKE_PATH"]
  suffix = "cmake/SCOREC/SCORECConfig.cmake"
  CONFIG_PATHS["PUMI_LIBDIR"] =  replace(str, suffix, "")

  # PUMIINTERFACE_PREFIX -> PUMIINTERFACE_LIBDIR
  str = CONFIG_PATHS["PUMIINTERFACE_PREFIX"]
  CONFIG_PATHS["PUMIINTERFACE_LIBDIR"] = joinpath(str, "lib")

  return nothing
end


readConfigFile()


