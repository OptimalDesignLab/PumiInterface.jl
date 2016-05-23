# file for functions related to parallel bookkeeping

function getParallelInfo(mesh::PumiMeshDG2)
# get information on number of shared entities
#TODO: update fields of mesh

  println("----- Entered getParallelInfo -----")

  myrank = mesh.myrank
  npeers = countPeers(mesh.m_ptr, 1)  # get edge peers
  peer_nums = zeros(Cint, npeers)
  getPeers(mesh.m_ptr, peer_nums)
  mesh.peer_parts = peer_nums
  mesh.npeers =  npeers
  mesh.send_reqs = Array(MPI.Request, npeers)  # array of Requests for sends
  mesh.send_waited = Array(Bool, npeers)
  mesh.recv_reqs = Array(MPI.Request, npeers)  # array of Requests for receives
  mesh.recv_waited = Array(Bool, npeers)
  mesh.send_stats = Array(MPI.Status, npeers)
  mesh.recv_stats = Array(MPI.Status, npeers)
  # count the number of edges shared with each peer
  partnums = zeros(Cint, 1)
  remotes = Array(Ptr{Void}, 1)
  counts = zeros(Int, npeers)  # hold the counts of the number of edges shared
                               # with each peer
  for i=1:mesh.numEdge
    edge = mesh.edges[i]
    if isShared(mesh.m_ptr, edge)
      nremotes = countRemotes(mesh.m_ptr, edge)
      @assert nremotes == 1
      getRemotes(partnums, remotes)

      part_idx = getElIndex(peer_nums, partnums[1])
      counts[part_idx] += 1
    end
  end
  mesh.peer_face_counts = counts

  # get the edges
  edges_local = Array(Array{Ptr{Void}, 1}, npeers)
  edges_remote = Array(Array{Ptr{Void}, 1}, npeers)  # is this still used?
  mesh.bndries_local = bndries_local = Array(Array{Boundary,1}, npeers)
  mesh.bndries_remote = bndries_remote = Array(Array{Boundary,1}, npeers)
  mesh.shared_interfaces = shared_interfaces = Array(Array{Interface,1}, npeers)
  for i=1:npeers
    edges_local[i] = Array(Ptr{Void}, counts[i])
    edges_remote[i] = Array(Ptr{Void}, counts[i])
    bndries_local[i] = Array(Boundary, counts[i])
    bndries_remote[i] = Array(Boundary, counts[i])
  end

  # get all the (local pointers to) edges in a single pass, 
  # even though they might not be needed
  curr_pos = ones(Int, npeers)  # hold the current position in each edge array
  for i=1:mesh.numEdge
    edge_i = mesh.edges[i]
    if isShared(mesh.m_ptr, edge_i)
      nremotes = countRemotes(mesh.m_ptr, edge_i)
      getRemotes(partnums, remotes)
      peer_i = partnums[1]

      idx = getElIndex(peer_nums, peer_i)  # get the part boundary index
      edges_local[idx][curr_pos[idx]] = edge_i
      curr_pos[idx] += 1
    end
  end

  # get boundary info for the edges
  for i=1:npeers
    if mesh.myrank > peer_nums[i]
      getEdgeBoundaries(mesh, edges_local[i], bndries_local[i])
      sort!(bndries_local[i])
    end
  end

  # get the remote edge pointers if I am the higher number process
  # tell MPI that the pointers are really Ints
  dtype = MPI.mpitype(Int)
  MPI.mpitype_dict[Ptr{Void}] = dtype

  down = Array(Ptr{Void}, 12)  # hold downward edges
  nrecvs = 0
  for i=1:npeers
    edges_i = edges_remote[i]  # array to hold the remote edge pointers
    if mesh.myrank > peer_nums[i]
      bndries_i = bndries_local[i]
      for j = 1:length(bndries_i)
        bndry_j = bndries_i[j]
        # get the edge
        el_j = mesh.elements[bndry_j.element]
        getDownward(mesh.m_ptr, el_j, apfEDGE, down)
        edge_j = down[bndry_j.face]

        # get the remote edge pointer
        nremotes = countRemotes(mesh.m_ptr, edge_j)
        @assert nremotes == 1
        getRemotes(partnums, remotes)
        edges_i[j] = remotes[1]
      end

      # now send the edges to the lower rank process
      mesh.send_reqs[i] = MPI.Isend(edges_i, peer_nums[i], 1, mesh.comm)
    else  # I am the receiving process
      mesh.recv_reqs[i] = MPI.Irecv!(edges_i, peer_nums[i], 1, mesh.comm)
      nrecvs += 1
    end
  end  # end loop over peers

  # get arrays of Requests
  recv_reqs_reduced = Array(MPI.Request, nrecvs)
  recv_peers = Array(Int, nrecvs)  # which process they came from
  pos = 1
  for i=1:npeers
    if mesh.myrank < peer_nums[i]
      recv_reqs_reduced[pos] = mesh.recv_reqs[i]
      recv_peers[pos] = i
      pos += 1
    end
  end


  # wait for all communication to finish
  for i=1:npeers
    if mesh.myrank > peer_nums[i]
      mesh.send_stats[i] = MPI.Wait!(mesh.send_reqs[i])
    else
      j, stat = MPI.Waitany!(recv_reqs_reduced)
      mesh.recv_stats[j] = stat
      peernum = recv_peers[j]
      getEdgeBoundaries(mesh, edges_remote[peernum], bndries_local[peernum])
    end
  end

  # now send Boundary info
  MPI.type_create(Boundary)
  for i=1:npeers
    mesh.send_reqs[i] = MPI.Isend(bndries_local[i], peer_nums[i], 1, mesh.comm)
    mesh.recv_reqs[i] = MPI.Irecv!(bndries_remote[i], peer_nums[i], 1, mesh.comm)
  end

  for i=1:npeers
    mesh.recv_stats[i] = MPI.Wait!(mesh.recv_reqs[i])
  end

  # now create Interfaces from the two Boundary arrays
  peer_offsets = Array(Int, npeers+1)  # record the starting element number
                                       # of the elements belonging to each peer
  mesh.local_element_lists = Array(Array{Int32, 1}, mesh.npeers)
  curr_elnum = mesh.numEl + 1
  for i=1:npeers
    peer_offsets[i] = curr_elnum
    curr_elnum, shared_interfaces[i] = numberBoundaryEls(curr_elnum, 
                                 bndries_local[i], bndries_remote[i])
    curr_elnum += 1

    mesh.local_element_lists[i] = getBoundaryElList(bndries_local[i])
  end
  peer_offsets[npeers+1] = curr_elnum
  mesh.shared_element_offsets = peer_offsets
  mesh.numGlobalEl = curr_elnum - 1
  mesh.numSharedEl = curr_elnum -1 - mesh.numEl  # number of non-local shared
                                                 # elnums
  mesh.local_element_counts = Array(Int, mesh.npeers)
  mesh.remote_element_counts = Array(Int, mesh.npeers)
  for i=1:mesh.npeers
    mesh.remote_element_counts[i] = peer_offsets[i+1] - peer_offsets[i]
    mesh.local_element_counts[i] = length(mesh.local_element_lists[i])
  end


  # create the dictonary that maps from locally owned element's numbers
  # to their adjacent non-local neighbors
  adj_dict, revadj = getLocalAdjacency(mesh)
  colordata = ColoringData(adj_dict, revadj)

  # delete the dangerous Ptr{Void} -> Int MPI type
  delete!(MPI.mpitype_dict, Ptr{Void})

  # wait for the sends to finish before exiting
  for i=1:npeers
    mesh.send_stats[i] = MPI.Wait!(mesh.send_reqs[i])
  end

  return colordata
end


function getEdgeBoundaries(mesh::PumiMeshDG2, edges::Array{Ptr{Void}}, 
                           bndries::Array{Boundary})
# get the array of Boundaryies for an array of edges
# edges is the array of the edge MeshEnities
# bndries is the array to be populated with the Boundaryies

  faces = Array(Ptr{Void}, 1)
  for i=1:length(edges)
    edge_i = edges[i]

    numFace = countAdjacent(mesh.m_ptr, edge_i, 2)  # should be count upward

    @assert( numFace == 1)

    getAdjacent(faces)
    facenum = getFaceNumber2(faces[1]) + 1
    edgenum = getEdgeNumber2(edge_i) + 1  # unneeded?
    edgenum_local = getEdgeLocalNum(mesh, edgenum, facenum)

    bndries[i] = Boundary(facenum, edgenum_local)
  end

  return nothing
end

function numberBoundaryEls(startnum, bndries_local::Array{Boundary}, bndries_remote::Array{Boundary})
# create Interfaces out of the local + remote Boundary arrays
# also numbers the remote elements with numbers > numEl, storing them in
# the elementR field of the Interface

  ninterfaces = length(bndries_local)
  interfaces = Array(Interface, ninterfaces)
  curr_elnum = startnum  # counter for 
  new_elR = 0 # number of new element to create
  for i=1:ninterfaces
    bndry_l = bndries_local[i]
    bndry_r = bndries_remote[i]
    old_iface_idx = isRepeated(bndries_remote, i)
    if old_iface_idx == 0
      new_elR = curr_elnum
      curr_elnum += 1
    else
      new_elR = interfaces[old_iface_idx].elementR
    end

    interfaces[i] = Interface(bndry_l.element, new_elR,  bndry_l.face, bndry_r.face, UInt8(1))
  end

  last_elnum = curr_elnum - 1
  return last_elnum, interfaces
end

function getBoundaryElList(bndries_local::Array{Boundary})
# get the list of elemements on the boundary

  nfaces = length(bndries_local)
  elnums = Array(Int32, nfaces)  # nfaces is the upper bound on 
                                 # the number of elements
  pos = 1 # current position in elnums
  for i=1:nfaces
    old_iface_idx = isRepeated(bndries_local, i)
    if old_iface_idx == 0
      elnums[pos] = bndries_local[i].element
      pos += 1
    end
  end

  return elnums[1:(pos-1)]
end

function isRepeated(bndries::Array{Boundary}, idx)
# figures out if a specified element has been seen before, returning
# the index of the boundary containing it if so, and returning zero if not
  el = bndries[idx].element
  for i=1:(idx-1)
    if bndries[i].element == el
      return i
    end
  end

  return 0
end

function getLocalAdjacency(mesh::PumiMeshDG2)

  # map from an element number of all the element numbers of the 
  # non-local elements
  # an element can have a maximum of 2 non-local neighbor elements
  adj_dict = Dict{Int, NTuple{2, Int}}()
  sizehint!(adj_dict, mesh.numSharedEl)

  # map from non-local element to local neighbors
  revadj = zeros(Int, mesh.numSharedEl, 2)

  for i=1:mesh.npeers
    ifaces_i = mesh.shared_interfaces[i]
    for j=1:length(ifaces_i)
      iface_j = ifaces_i[j]
      el_local = iface_j.elementL
      el_nonlocal = iface_j.elementR

      if !haskey(adj_dict, el_local)
        adj_dict[el_local] = (Int(el_nonlocal), 0)
      else
        old_tuple = adj_dict[el_local]
        @assert old_tuple[2] == 0
        adj_dict[el_local] = (old_tuple[1], Int(el_nonlocal))
      end  # end if-else

      nonlocal_idx = el_nonlocal - mesh.numEl
      first_neighbor = revadj[nonlocal_idx, 1]
      second_neighbor = revadj[nonlocal_idx, 2]
      if first_neighbor == 0
        revadj[nonlocal_idx, 1] = el_local
      elseif second_neighbor == 0
        revadj[nonlocal_idx, 2] = el_local
      else
        throw(ErrorException("Too many adjacent elements for peer $i interface $j"))
      end

    end # end loop over interfaces
  end  # end loop over peers

  return adj_dict, revadj
end



