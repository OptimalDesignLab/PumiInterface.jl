# functions for creating a new mesh object using the same underlying
# apf::Mesh as an existing mesh object.


"""
  This function returns a new mesh object that uses the same underlying
  apf::Mesh as an existing mesh object

  **Inputs**

   * mesh: the existing mesh
   * sbp: an SBP operator
   * opts: options dictonary
   * sbpface: an `AbstractFace` that is compatible with `sbp`. If `sbp` is
              the same `sbp` used to create `mesh`, then this argument is not
              required

  **Outputs**

   * newmesh: the new mesh object

  Note that `sbp` can be a different SBP operator than was used to create `mesh`,
  either a different type of a different degree.

"""
function copy_mesh(mesh::PumiMesh, sbp::AbstractOperator, opts,
                   sbpface::AbstractFace=mesh.sbpface)

  error("generic fallback reached")
end


function copy_mesh(mesh::PumiMeshDG2, sbp::AbstractOperator, opts,
                   sbpface::AbstractFace=mesh.sbpface)

  return PumiMeshDG2(mesh, sbp, opts, sbpface)
end


function copy_mesh(mesh::PumiMeshDG3, sbp::AbstractOperator, opts,
                   sbpface::AbstractFace=mesh.sbpface)

  return PumiMeshDG3(mesh, sbp, opts, sbpface)

end
