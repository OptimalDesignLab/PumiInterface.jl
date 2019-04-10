# mesh adaptation functions

"""
  This function takes a vector containing the desired size for each element
  and turns in into an apf::Field for mesh adaptation.

  This is required because MeshAdapt runs several passes, so for all passes
  after the first, some of the mesh entities it will query will not exist
  on the original mesh.  Using the apf::Field allows the size field to be
  interpolated to the new mesh after each pass.

  The apf::Field has the same apf::FieldShape as the *coordinate* field
  of the mesh (not the solution field).

  **Inputs**

   * mesh: a mesh object
   * el_sizes: desired size for each element

  **Outputs*

   * an apf::Field object
"""
function getSizeField(mesh::PumiMesh, el_sizes::AbstractVector)

  # create apf::Field
  fsize = apf.createPackedField(mesh, "size_field", 1, mesh.coordshape_ptr)

  apf.zeroField(fsize)

  # figure out which entities have nodes on them
  has_nodes = Array{Bool}(mesh.dim+1)
  for d=0:mesh.dim
    has_nodes[d+1] = apf.hasNodesIn(mesh.coordshape_ptr, d)
  end

  vals = zeros(Cdouble, 1)
  down_entities = Array{Ptr{Void}}(12)  # apf::Downward
  # assign size to field
  # for each entity, we specify the *minimum* size of any element it
  # bounds
  for i=1:mesh.numEl
    el_i = mesh.elements[i]
    for d=0:mesh.dim
      if has_nodes[d+1]
        ndown = apf.getDownward(mesh.m_ptr, el_i, d, down_entities)
        for j=1:ndown
          apf.getComponents(fsize, down_entities[j], 0, vals)
          if vals[1] == 0
            vals[1] = el_sizes[i]
          else
            vals[1] = min(vals[1], el_sizes[i])
          end

          apf.setComponents(fsize, down_entities[j], 0, vals)
        end  # end loop j
      end  # end if
    end  # end loop d
  end  # end loop i

  # take minimum along partition boundaries
  apf.reduceField(fsize, mesh.shr_ptr, 1);

  # return field
  return fsize
end


#TODO: expose more of MA::Inputs options
"""
  This is the main entry point for runing mesh adaptation (h refinement).
  
  This function takes in a mesh and a vector describing the desired size
  for each element and returns a new mesh object.  The original mesh object
  is finalized.

  **Inputs**

   * mesh: a mesh object
   * sbp: SBP operator (must be the same operator used to construct
          `mesh`
   * opts: options dictionary
   * el_sizes: desired size for each element (vector of Float64)
   * u_vec: solution field to interpolate to the adapted mesh (optional)

  **Outputs**

   * mesh_new: a new mesh object
   * u_vec_new: `u_vec` interpolated to the new mesh.  A vector of length
                 zero if `u_vec` was not supplied

  **Keyword Arguments**

   * free_mesh: if true, free the old mesh object.  Note that the old mesh
                object cannot correctly call Pumi functions after adaptation.
                Default true,  Only set to false if you are knowledgeable
                about PdePumiInterface's internals.

  **Implementation Notes**

    The original and new mesh objects have the same apf::Mesh*.

"""
function adaptMesh(oldmesh::PumiMeshDG2, sbp, opts, el_sizes::AbstractVector, u_vec::AbstractVector=zeros(0); free_mesh::Bool=true)

  # run the mesh adaptation
  _adaptMesh(oldmesh, el_sizes, u_vec; free_mesh=free_mesh)

  # now construct new mesh object
  # To avoid the mesh reference count reaching zero and destroying the
  # mesh, we have to create the new mesh object before finalizing the
  # old one
  newmesh = PumiMeshDG2(oldmesh, sbp, opts)

  if length(u_vec) > 0
    u_vec_new = zeros(Float64, newmesh.numDof)
    retrieveSolutionFromMesh_interp(newmesh, u_vec_new)
  else
    u_vec_new = u_vec
  end

  if free_mesh
    finalize(oldmesh)
  end


  return newmesh, u_vec_new
end

# method for 3D DG meshes
function adaptMesh(oldmesh::PumiMeshDG3, sbp, opts, el_sizes::AbstractVector, u_vec::AbstractVector=zeros(0); free_mesh::Bool=true)

  # run the mesh adaptation
  _adaptMesh(oldmesh, el_sizes, u_vec, free_mesh=free_mesh)

  newmesh = PumiMeshDG3(oldmesh, sbp, opts)

  if length(u_vec) > 0
    u_vec_new = zeros(Float64, newmesh.numDof)
    retrieveSolutionFromMesh_interp(newmesh, u_vec_new)
  else
    u_vec_new = u_vec
  end

  if free_mesh
    finalize(oldmesh)
  end


  return newmesh, u_vec_new
end



"""
  This does most of the work involved in mesh adaptation.  Specifically,
  it does all of the work that does not require knowing exactly which
  type of mesh object is being adapted (so, basically constructing the
  new mesh object)

  This function deletes most of the existing Numberings and Fields on the
  mesh in preparation for mesh reinitialization.

  **Inputs**

   * mesh: a mesh object
   * el_sizes: desired size for each element (vector)
   * u_vec: solution field to interpolate to the adapted mesh (optional)

  **Keyword Arguments**

   * free_mesh: if true, destroys old data on the mesh data structure after
                adaptation
"""
function _adaptMesh(mesh::PumiMesh, el_sizes::AbstractVector, u_vec::AbstractVector; free_mesh::Bool=true)

  # get size function
  size_f = getSizeField(mesh, el_sizes)
  isofunc = apf.createIsoFunc(mesh.m_ptr, size_f)

  # interpolate size and possibly the solution to new mesh
  soltrans = apf.createSolutionTransfers()
  apf.addSolutionTransfer(soltrans, size_f)
  if length(u_vec) > 0
    saveSolutionToMesh(mesh, u_vec)
    apf.addSolutionTransfer(soltrans, mesh.fnew_ptr)  # transfer field to new mesh
  end

  # configure input
  ma_input = apf.configureMAInput(mesh.m_ptr, isofunc, soltrans)

  # run mesh adaptation
  apf.runMA(ma_input)  # this function deletes ma_input

  # cleanup intermediate data
  apf.deleteSolutionTransfers(soltrans)
  apf.deleteIsoFunc(isofunc)
  apf.destroyField(mesh, size_f) 


  # because we are going to reinitialize the mesh, and Pumi will return
  # an existing Numbering (and possibly Field?) if a new one is created
  # with the same name, we have to delete all existing Numberings/Fields
  if free_mesh
    apf.destroyNumberings(mesh.m_ptr)

    for f in copy(mesh.fields.user)
      if f == mesh.fnew_ptr
        continue
      end
      destroyField(mesh, f)
    end
  end

  return nothing
end


"""
  Returns a vector containing the size of each element.  Useful as
  the starting point for the size vector required by [`adaptMesh`](@ref)

  **Inputs**

   * mesh: mesh object

  **Outputs**

   * el_sizes: vector of length mesh.numEl, containing the size of each
               element (really the average edge length)
"""
function getElementSizes(mesh::PumiMesh)

  el_sizes = zeros(Cdouble, mesh.numEl)
  apf.getAvgElementSize(mesh.m_ptr, mesh.el_Nptr, el_sizes)

  return el_sizes
end
