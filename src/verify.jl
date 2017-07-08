# run time verification of mesh topology

"""
  This function performs some checks on the final mesh object to make sure
  it is sane
"""
function checkFinalMesh(mesh::PumiMesh)

  checkFaceCount(mesh)

  # add more checks here

  return nothing
end


function checkConnectivity(mesh::PumiMesh)

  checkVertConnectivity(mesh)
  checkEdgeConnectivity(mesh)
  if mesh.dim > 2
    checkFaceConnectivity(mesh)
  end

  return nothing
end

"""
  This function verifies that at least one vertex of every element is shared
  with a different element (ie. no elements dangling by a vertex)m, otherwise
  an exception is thrown

  If this doesn't hold, there is a chance Pumi will call abort(), but 
  either way the user will be notified something is wrong with the mesh.

  Inputs:
    mesh: a Pumi mesh, must have the mesh.elements field populated
"""
function checkVertConnectivity(mesh::PumiMesh)

  down_verts = Array(Ptr{Void}, 12)

  for i=1:mesh.numEl
    el_i = mesh.elements[i]

    nverts = getDownward(mesh.m_ptr, el_i, 0, down_verts)

    nel_sum = 0
    for j=1:nverts
      nel = countAdjacent(mesh.m_ptr, down_verts[j], mesh.dim)
      nel_sum += nel
    end

    if nel_sum <= nverts
      throw(ErrorException("element $i is dangling by a vertex"))
    end
  end

  return nothing
end

"""
  This function verifies no elements are dangling by an edge (see 
  checkVertConnectivity for the idea).  This test may not be distinct from
  checkVertConnectivity

  Inputs:
    mesh: a PumiMesh
"""
function checkEdgeConnectivity(mesh::PumiMesh)

  down_edges = Array(Ptr{Void}, 12)

  for i=1:mesh.numEl
    el_i = mesh.elements[i]
    nedges = getDownward(mesh.m_ptr, el_i, 1, down_edges)

    nel_sum = 0
    for j=1:nedges
      nel = countAdjacent(mesh.m_ptr, down_edges[j], mesh.dim)
      nel_sum += nel
    end

    if nel_sum <= nedges
      throw(ErrorException("element $i is dangling by an edge"))
    end

  end

  return nothing
end

"""
  This function verifies no elements are dangling by a face.  This works
  for both 2D and 3D meshes because msh.faces is an alias for mesh.edges.
  So it will check the same thing as checkEdgeConnectivity.

  Inputs:
    mesh: a PumiMesh
"""
function checkFaceConnectivity(mesh::PumiMesh)

  down_faces = Array(Ptr{Void}, 12)

  for i=1:mesh.numEl
    el_i = mesh.elements[i]
    nfaces = getDownward(mesh.m_ptr, el_i, 1, down_faces)

    nel_sum = 0
    for j=1:nfaces
      nel = countAdjacent(mesh.m_ptr, down_faces[j], mesh.dim)
      nel_sum += nel
    end

    if nel_sum <= nfaces
      throw(ErrorException("element $i is dangling by a face"))
    end

  end

  return nothing
end


"""
  This function verifies than all the faces of each element appear
  in mesh.interfaces, mesh.bndryfaces, or mesh.bndries_local exactly once
"""
function checkFaceCount(mesh::PumiMesh)


  seen_faces = falses(mesh.numFacesPerElement, mesh.numEl)

  for i=1:mesh.numInterfaces
    iface_i = mesh.interfaces[i]
    j = iface_i.faceL
    k = iface_i.elementL
    if !seen_faces[j, k]
      seen_faces[j, k] = true
    else
      throw(ErrorException("face $j of element $k appears twice"))
    end

    j = iface_i.faceR
    k = iface_i.elementR
    if !seen_faces[j, k]
      seen_faces[j, k] = true
    else
      throw(ErrorException("face $j of element $k appears twice"))
    end
  end

  for i=1:mesh.numBoundaryFaces
    bndry_i = mesh.bndryfaces[i]
    j = bndry_i.face
    k = bndry_i.element

    if !seen_faces[j, k]
      seen_faces[j, k] = true
    else
      throw(ErrorException("face $j of element $k appears twice"))
    end
  end

  for peer=1:mesh.npeers
    bndries_peer = mesh.bndries_local[peer]
    for i=1:length(bndries_peer)
      bndry_i = bndries_peer[i]
      j = bndry_i.face
      k = bndry_i.element
      if !seen_faces[j, k]
        seen_faces[j, k] = true
      else
        throw(ErrorException("face $j of element $k appears twice"))
      end
    end
  end

  for i=1:mesh.numEl
    for j=1:mesh.numFacesPerElement
      if !seen_faces[j, i]
        println("element $i vert coords = \n", mesh.vert_coords[:, :, i])
        println("face_verts = \n", mesh.topo.face_verts)
        println("mesh.interfaces = \n")
        for k=1:mesh.numInterfaces
          println("interface ", k, " = ", mesh.interfaces[k])
        end
        println("mesh.bndryfaces = \n")
        for k=1:mesh.numBoundaryFaces
          println("boundary ", k, " = ", mesh.bndryfaces[k])
        end
        throw(ErrorException("face $j of element $i is missing"))
      end
    end
  end

  return nothing
end

