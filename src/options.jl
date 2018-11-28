# set default values for most options
# these values should never really be used because PDESolver should set them

"""
  This function supplies default values for the options dictionary.
  Options already defined are not overwritten.

  The keys it does not supply default values for are:

   * smb_name
   * dmg_name
   * order
   * coloring_distance
   * numBC
   * other boundary condition related keys

  **Inputs**

   * opts: a dictionary
"""
function set_defaults(opts)

  get!(opts, "run_type", 4)
  get!(opts, "verify_coloring", true)
  get!(opts, "use_edge_res", false)
  get!(opts, "write_edge_vertnums", false)
  get!(opts, "write_boundarynums", false)
  get!(opts, "write_dxidx", false)
  get!(opts, "write_jac", false)
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
  get!(opts, "write_element_vertnums", false)
  get!(opts, "exact_visualization", false)

  # force the use of the metric calculation method for linear elements
  # typically the curvilinear method is used even for linear (straight-sided)
  # elements.
  get!(opts, "use_linear_metrics", false)
  get!(opts, "reordering_start_coords", [0.0, 0.0, 0.0])


  return nothing
end

"""
  Create an options dictionary with the default options

  **Outputs**

   * opts: the dictionary
"""
function get_defaults()
  opts = Dict{Any, Any}()
  set_defaults(opts)
  return opts
end
