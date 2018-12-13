# a very simple MPI tag manager

"""
  Type that manages assigning MPI tags
"""
struct MPITagManager
  used_tags::Vector{Cint}
  start_tag::Int
end


"""
  Constructor

  **Inputs**

   * start_tag: tag to start counting from, default 5000
"""
function MPITagManager(start_tag=5000)
  used_tags = Array{Cint}(0)

  return MPITagManager(used_tags, start_tag)
end


"""
  Returns the next free tag
"""
function getNextTag(mgr::MPITagManager)

  if length(mgr.used_tags) == 0
    lasttag = mgr.start_tag
  else
    lasttag = mgr.used_tags[end]
  end

  newtag = lasttag + 1
  push!(mgr.used_tags, newtag)

  return newtag
end


#------------------------------------------------------------------------------
# default manager

global const TagManager = MPITagManager()
