# set default values for most options
# these values should never really be used because PDESolver should set them

function set_defaults(opts)

  get!(opts, "run_type", 4)
  get!(opts, "verify_coloring", true)
  get!(opts, "use_edge_res", false)
  get!(opts, "write_edge_vertnums", false)
  get!(opts, "write_boundarynums", false)
  get!(opts, "write_dxidx", false)
  get!(opts, "write_sparsity", false)
  get!(opts, "write_offsets", false)
  get!(opts, "write_counts", false)
  get!(opts, "write_sparsity_nodebnds", false)
  get!(opts, "write_dofs", false)
  get!(opts, "write_interfaces", false)
  get!(opts, "write_boundaries", false)
  get!(opts, "write_sharedboundaries", false)
  get!(opts, "write_face_vertnums", false)
  get!(opts, "write_coords", false)
  get!(opts, "reordering_algorithm", "default")

  return nothing
end

function get_defaults()
  opts = Dict{Any, Any}()
  set_defaults(opts)
  return opts
end
