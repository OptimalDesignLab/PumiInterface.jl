# interface functions 

"""
  Numbers all points on the specified surfaces.  All points on the
  surface of the domain not specified are numbered with the number of
  points on the specified surface + 1 (the numbering is 1-based)

  Currently this supports linear coordinate fields only

  **Inputs**

   * mesh
   * bc_nums: array of boundary condition numbers that define the surface(s)
              to number

  **Outputs**

   * numFaceNodes: number of nodes on the specified surface(s)
   * n_face: a Ptr{Void} (really an apf::Numbering*) 
   * face_verts: array of apf::MeshEntity* in the order they appear in the
                 surface numbering, length numFaceNodes
"""
function numberSurfacePoints{I<:Integer}(mesh::PumiMeshDG, bc_nums::AbstractVector{I})

  @assert mesh.coord_order == 1
  for i in bc_nums  # check that the BC numbers are valid
    @assert i <= length(mesh.bndry_offsets - 1)
    @assert i > 0
  end
  # number the (unique) coordinate nodes on the faces

  # we can't guarantee an existing numbering with the same name was
  # created using the same bc_nums, so delete any existing numbering
  const numbering_name = "warp_surf_nums"
  n_old = findNumbering(mesh.m_ptr, numbering_name)
  if n_old != C_NULL
    destroyNumbering(n_old)
  end
  n_face = createNumberingJ(mesh.m_ptr, numbering_name, 
                             mesh.coordshape_ptr, 1)
  num_i = 1
  verts = Array(Ptr{Void}, 4)
  face_verts = Array(Ptr{Void}, 0)  # TODO: add sizehint
  coords = zeros(Float64, 3)
  for i in bc_nums
    start_idx = mesh.bndry_offsets[i]
    end_idx = mesh.bndry_offsets[i+1] - 1
    rng = start_idx:end_idx

    for j in rng
      bndry_j = mesh.bndryfaces[j]
      el_j = bndry_j.element
      el_ptr = mesh.elements[el_j]

      # get the vertices
      getDownward(mesh.m_ptr, el_ptr, 0, verts)

      for k=1:size(mesh.topo.face_verts, 1)
        v_k = verts[mesh.topo.face_verts[k, bndry_j.face]]

        if !isNumbered(n_face, v_k, 0, 0)
          getPoint(mesh.m_ptr, v_k, 0, coords)
          numberJ(n_face, v_k, 0, 0, num_i)
          push!(face_verts, v_k)
          num_i += 1
        end
      end  # end loop k
    end  # end loop j
  end  # end loop i

  # number all remaining vertices with the number of face verts + 1
  for i=1:mesh.numVert
    v_i = mesh.verts[i]

    if !isNumbered(n_face, v_i, 0, 0)
      numberJ(n_face, v_i, 0, 0, num_i)
    end
  end

  numFacePts = num_i - 1

  return numFacePts, n_face, face_verts
end


