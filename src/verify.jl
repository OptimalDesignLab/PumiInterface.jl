# run time verification of mesh topology

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

