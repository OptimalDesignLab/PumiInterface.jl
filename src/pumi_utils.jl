# utilities for managing Pumi lbraries

"""
  Returns the path to the /lib directory of the Pumi installation
"""
function getPumiLibDir()
  f = open(joinpath(dirname(@__FILE__), "PUMI_CMAKE_PATH"), "r")
  str = chomp(read(f, String))
  close(f)

  suffix = "cmake/SCOREC/SCORECConfig.cmake"
  str = replace(str, suffix, "")

  return str
end


