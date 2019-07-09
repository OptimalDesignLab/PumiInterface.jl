# run time verification of mesh topology

"""
  This function performs some checks on the final mesh object to make sure
  it is sane
"""
function checkFinalMesh(mesh::PumiMesh)

  checkFaceCount(mesh)
  if mesh.isDG
    checkElementNumbering(mesh, mesh.el_Nptr)
    checkDofNumbering(mesh, mesh.dofnums_Nptr)
  end

  # add more checks here

  return nothing
end

"""
  This function runs all the test of mesh connectivity
"""
function checkConnectivity(mesh::PumiMesh)
#=
  checkVertConnectivity(mesh)
  checkEdgeConnectivity(mesh)
  if mesh.dim > 2
    checkFaceConnectivity(mesh)
  end
=#
#  checkContiguity(mesh)

  return nothing
end

"""
  This function verifies that at least one vertex of every element is shared
  with a different element (ie. no elements dangling by a vertex), otherwise
  an exception is thrown

  If this doesn't hold, there is a chance Pumi will call abort(), but 
  either way the user will be notified something is wrong with the mesh.

  Inputs:
    mesh: a Pumi mesh, must have the mesh.elements field populated
"""
function checkVertConnectivity(mesh::PumiMesh)

  down_verts = Array{Ptr{Void}}(12)

  for i=1:mesh.numEl
    el_i = mesh.elements[i]

    nverts = apf.getDownward(mesh.m_ptr, el_i, 0, down_verts)

    nel_sum = 0
    for j=1:nverts
      nel = apf.countAdjacent(mesh.m_ptr, down_verts[j], mesh.dim)
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

  down_edges = Array{Ptr{Void}}(12)

  for i=1:mesh.numEl
    el_i = mesh.elements[i]
    nedges = apf.getDownward(mesh.m_ptr, el_i, 1, down_edges)

    nel_sum = 0
    for j=1:nedges
      nel = apf.countAdjacent(mesh.m_ptr, down_edges[j], mesh.dim)
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

  down_faces = Array{Ptr{Void}}(12)

  for i=1:mesh.numEl
    el_i = mesh.elements[i]
    nfaces = apf.getDownward(mesh.m_ptr, el_i, 1, down_faces)

    nel_sum = 0
    for j=1:nfaces
      nel = apf.countAdjacent(mesh.m_ptr, down_faces[j], mesh.dim)
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
        throw(ErrorException("face $j of element $i is missing"))
      end
    end
  end

  return nothing
end

"""
  This function checks that every mesh element is reachable via face
  adjacencies.  Throws an exception otherwise

  The main case this guards against is having two disjoint sets of elements.

  This algorithm walks the mesh using face adjacencies and records if
  each element was seen.  If after finishing the walk some elements are not
  seen yet, throw an exception
"""
function checkContiguity(mesh::PumiMesh)

  curr_el = mesh.elements[1]
  el_Nptr = mesh.el_Nptr
  curr_elnum = apf.getNumberJ(el_Nptr, curr_el, 0, 0) + 1

  # temporary array for apf.getBridgeAdjacent
  adjacent_els = Array{Ptr{Void}}(mesh.numFacesPerElement)

  seen_els = falses(mesh.numEl)  # record whether or not an element either is
                                 # or has previously been in the queue

  # the queue of elements to process
  queue = FIFOQueue{Ptr{Void}}(size_hint = div(mesh.numEl, 4))

  # put first element in the queue
  push!(queue, curr_el)
  seen_els[curr_elnum] = true


  while length(queue) > 0
    curr_el = pop!(queue)

    # get face adjacent elements
    n = apf.countBridgeAdjacent(mesh.m_ptr, curr_el, mesh.dim - 1, mesh.dim)
    @assert n <= mesh.numFacesPerElement
    apf.getBridgeAdjacent(adjacent_els)

    # if not already seen, add to queue
    for i=1:n
      el_i = adjacent_els[i]
      new_elnum = apf.getNumberJ(el_Nptr, el_i, 0, 0) + 1

      # use seen_els to avoid growing the size of the queue unnecessarily
      if !seen_els[new_elnum]
        seen_els[new_elnum] = true
        push!(queue, el_i)
      end
    end

  end  # end while loop
        
  # make sure all elements were seen
  unseen_count = 0
  for i=1:mesh.numEl
    if !seen_els[i] 
      unseen_count += 1
    end
  end

  if unseen_count > 0
    throw(ErrorException("Error: Mesh is not contiguous: $(unseen_count) elements were not seen while walking the mesh"))
  end

  return nothing
end

"""
  This function verifies that the element topologies are sufficiently consistent
  to support second order coordinate fields.  In particular, it means that the
  element faces are defined by the same vertices (but not necessarily in the
  same order

  Inputs
    topo1: one topology
    topo2: the other topology

  Outputs:
    none
"""
function checkTopologyConsistency(topo1::ElementTopology, topo2::ElementTopology)

  nfaces = size(topo1.face_verts, 2)
  for i=1:nfaces
    @assert sort!(topo1.face_verts[:, i]) == sort!(topo2.face_verts[:, i])
  end

  return nothing
end

"""
  Verify the determinent of the mapping jacobian is positive
"""
function checkMapping(mesh::PumiMesh)

  for i=1:mesh.numEl
    for j=1:mesh.numNodesPerElement
      if real(mesh.jac[j, i]) < 0
        println("coords_i = \n", real(mesh.coords[:, :, i]))
        println("vert_coords_i = \n", real(mesh.vert_coords[:, :, i]))
        println("dxidx_i = \n", real(mesh.dxidx[:, :, :, i]))
        println("jac_i = \n", real(mesh.jac[:, i]))
        throw(ErrorException("element $i node $j has negative mapping jacobian determinant"))
      end
    end
  end

  return nothing
end
