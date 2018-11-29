# data structures for parallel communication

"""
  Abstraction around an MPI send or receive buffer (either one, not both at
  the same time).  It enables transforming from an N dimensional array
  representation of the data to an N2 represtation of the data, where the
  trailing dimensions of the array are flattened.  Most commonly, this
  enables transforming between the 3D array form and the 1D array form.
  As a further generalization, each index can have an arbitrary rank array
  of data.  This data will be packed into a contiguous array to be send via
  MPI communication.

  **Static Parameters**

   * T: element type of the buffer
   * N: dimensionality of the buffer
   * N2: number of dimensions of the indicies in `local_indices`

  The dimension of the arrays are:

   * The data at each index is a N - 1 dimension array
   * the MPI buffer has dimension N (N - 1 for the data, final dimension for
     the apf::MeshEntity* the data lives on)
   * The original array has dimension N + N2 - 1.  When transforming from the
     higher dimensional array to the MPI buffer, the final N2 dimensions are
     flattened to a single dimension according to the apf::MeshEntity the
     data lives on.

  Overall, this data structure enables sending an array of arbitrary rank
  arrays over MPI.

  **Fields**

   * peernum: MPI rank of peer to send/receive from
   * entities: the apf::MeshEntity* whose associated data will be send/received
   * local_indices: (for send buffer only) array of `CartesianIndex`es of
                      
"""
mutable struct PeerData{T, N, N2}
  peernum::Int
  entities::Vector{Ptr{Void}}  # remote apf::MeshEntity*
  _entities_local::Vector{Ptr{Void}}  # local apf::MeshEntity* (unneeded)
  local_indices::Vector{CartesianIndex{N2}}
  colptr::Vector{Int}  # colptr[i]:colptr[i+1]-1 describes the range of
                       # indices associated with entities[i]
  vals::Array{T, N}      # currently we only support sending one scalar per
                         # entity.  Eventually we should support sending an
                         # arbitrary dimension tensor (packed into an array)
  req::MPI.Request
  req_waited::Bool
  tag::Int
  comm::MPI.Comm
  myrank::Int
  commsize::Int
end

"""
  Constructor for [`PeerData`](@ref) objects that will be used for sending data
"""
function PeerData(::Type{T}, peernum::Integer, tag::Integer, comm::MPI.Comm) where {T}

  N = 2
  N2 = 2
  entities = Array{Ptr{Void}}(0)
  _entities_local = Array{Ptr{Void}}(0)
  local_indices = Array{CartesianIndex{N2}}(0)
  colptr = Array{Int}(1); colptr[1] = 1
  vals = Array{T, N}(ntuple(i -> 0, Val{N}))
  req = MPI.REQUEST_NULL
  req_waited = true
  myrank = MPI.Comm_rank(comm)
  commsize = MPI.Comm_size(comm)

  return PeerData{T, N, N2}(peernum, entities, _entities_local, local_indices, colptr,
                            vals, req, req_waited, tag, comm, myrank, commsize)
end


"""
  Allocates the send or receive buffer

  **Inputs**

   * data: the `PeerData` object
   * dims: a tuple of length N.  The first N-1 entries will be used as the
           leading dimensions of the buffer.  The final dimension is the length
           of the `data.entities` array.
"""
function allocateBuffer(data::PeerData{T, N, N2}, dims::NTuple) where {T, N, N2}

  final_dim = length(data.entities)
  _dims = (dims..., final_dim)
  @assert length(_dims) == N
  data.vals = zeros(T, _dims)

  return nothing
end



"""
  Constructor for [`PeerData`](@ref) that will be used for receiving data

  **Inputs**

   * T: element type of buffer
   * peernum: MPI rank of process to receive data from
   * nval: number of entities to receive
   * comm: MPI communicator
"""
function PeerData(::Type{T}, peernum::Integer, dims::NTuple, nval::Integer,
                  tag::Integer, comm::MPI.Comm) where {T}

  N = 2
  N2 = 2
  entities = Array{Ptr{Void}}(nval)
  _entities_local = Array{Ptr{Void}}(0)  #TODO: document
  local_indices = Array{CartesianIndex{N2}}(0)  #TODO: document
  colptr = Array{Int}(0)  # TODO: document
  vals = Array{T, N}(dims..., nval)  # this might be type-unstable
                                              # is the overall constructor
                                              # still type-unstable?
  req = MPI.REQUEST_NULL
  req_waited = true
  myrank = MPI.Comm_rank(comm)
  commsize = MPI.Comm_size(comm)

  return PeerData{T, N, N2}(peernum, entities, _entities_local, local_indices,
                            colptr, vals, req, req_waited, tag, comm, myrank,
                            commsize)
end


import MPI: Isend, Irecv!, Wait!

function Isend(data::PeerData)

  @assert data.req_waited

  data.req = MPI.Isend(data.vals, data.peernum, data.tag, data.comm)
  data.req_waited = false

  return nothing
end


function Irecv!(data::PeerData)

  @assert data.req_waited

  data.req = MPI.Irecv!(data.vals, data.peernum, data.tag, data.comm)
  data.req_waited = false

  return nothing
end


function Wait!(data::PeerData)

  if !data.req_waited
    MPI.Wait!(data.req)
    data.req_waited = true
  end

  return nothing
end



"""
  Object to manage communication with all other MPI processes.

  Using this type has two stages.  The first is a setup stage where the
  user supplies the (MPI rank, remote MeshEntity*) where data will be sent
  to.  This is the purpose of the [`pushKey!`](@ref) function.  Once this
  is done, [`exchangeKeys`](@ref) must be called.  This completes the setup
  phase.

  The second phase is the repeated communication phase.  In this phase,
  [`pushValue`](@ref) is used to add values to the MPI send buffer.
  The data must be added to the buffer in the same order as the entities
  were added during the setup phase.  This is necessary for efficiency.
  Use `ScatterData.send[i].local_indices`to help with this.
  After all values have been added, call [`sendParallelData`](@ref) to
  start communication and [`receiveParallelData`](@ref) to finish it.

  See [`initSendToOwner`](@ref) for an example of adding entities during
  the setup phase.

  **Static Parameters**

   * T: element type of the buffer
   * N: dimensionality of the buffer
   * N2: number of dimensions of the indicies in `local_indices`
   * N3: N - 1 (dimension of data at each index in the buffer)

  **Fields**
  
   * send: vector of [`PeerData`](@ref) objects used for sending data
   * recv: vector of [`PeerData`](@ref) objects used for receiving data
   * peernums_send: array of Cints containing the MPI ranks of the processes
                    to send to
   * peernums_recv: array of Cints containing the MPI ranks of the processes to
                    receive from
   * comm: MPI communicator
   * tag: MPI tag used for communication
   * dims: NTuple containing the leading dimensions of the MPI buffer
"""
struct ScatterData{T, N, N2, N3}
  send::Vector{PeerData{T, N, N2}}
  recv::Vector{PeerData{T, N, N2}}
  peernums_send::Vector{Cint}  # MPI ranks of processes to send data to
  peernums_recv::Vector{Cint}  # MPI ranks of processes to receive data from
  comm::MPI.Comm
  tag::Int
  dims::NTuple{N3, Int}  # leading dimensions of send and receive buffers
  curridx::Vector{Int}  # used for pushing into pre-sized array, contains the
                        # current index
end


function ScatterData(::Type{T}, dims::NTuple, comm::MPI.Comm) where {T}

  N = 2
  N2 = 2
  N3 = N - 1

  send = Array{PeerData{T, N, N2}}(0)
  recv = Array{PeerData{T, N, N2}}(0)
  peernums_send = Array{Cint}(0)
  peernums_recv = Array{Cint}(0)
  tag = getNextTag(TagManager)
  curridx = Array{Int}(0)

  return ScatterData{T, N, N2, N3}(send, recv, peernums_send, peernums_recv,
                                   comm, tag, dims, curridx)
end



#------------------------------------------------------------------------------
# ScatterData API

"""
  Get the peer index from the peer MPI rank (using `data.peernums_send`)

  **Inputs**

   * data: a [`ScatterData`](@ref)
   * peernum: MPI rank

  **Outputs**

   * peeridx: index of peer in `data.peernums`. 0 if not found
"""
function getPeerIdx(data::ScatterData, peernum::Integer)

  idx = 0
  for i=1:length(data.peernums_send)
    if data.peernums_send[i] == peernum
      idx = i
    end
  end

  return idx
end

"""
  Like [`getPeerIdx`](@ref), but searches the `data.peernums_recv` array.
"""
function getRemotePeerIdx(data::ScatterData, peernum::Integer)

  idx = 0
  for i=1:length(data.peernums_recv)
    if data.peernums_recv[i] == peernum
      idx = i
    end
  end

  return idx
end



"""
  Adds a MeshEntity to the list of entities that will send data the given
  process.

  **Inputs**

   * data: a [`ScatterData`](@ref) object
   * e_local: the local apf::MeshEntity*
   * e_remote: the remote apf::MeshEntity*
   * peernum: the MPI rank of the process that `e_remote` lives on
   * localidx: either a single `CartesianIndex` or a vector of
               `CartesianIndex`es`describing the index(es) in the
               multi-dimensional array corresponding to `e_local`
"""
function pushKey!(data::ScatterData{T, N, N2}, e_local::Ptr{Void},
                  e_remote::Ptr{Void}, peernum::Integer,
                  localidx::Union{T2, AbstractVector{T2}}) where {T, N, N2, T2 <: CartesianIndex}

  peeridx = getPeerIdx(data, peernum)

  if peeridx == 0
    push!(data.send, PeerData(T, peernum, data.tag, data.comm))
    push!(data.peernums_send, peernum)
    peeridx = length(data.send)
  end

  data_i = data.send[peeridx]
  push!(data_i.entities, e_remote)
  push!(data_i._entities_local, e_local)
  pushorappend!(data_i.local_indices, localidx)
  push!(data_i.colptr, length(data_i.local_indices) + 1)
  push!(data.curridx, 1)

  return nothing
end


"""
  Helper function for appending to vectors.

  **Inputs**

   * arr: vector to append to
   * val: either a single value of an array of values.  All the values will
          be pushed onto the end of `arr`
"""
function pushorappend!(arr::AbstractArray, val::AbstractArray)
  return append!(arr, val)
end

function pushorappend!(arr::AbstractArray, val)
  return push!(arr, val)
end


"""
  Adds a new MPI process to the list of processes that data will be received
  from.  It is allowed to call this function with the same `peernum` more than
  once

  **Inputs**

   * data: [`ScatterData`](@ref) object
   * peernum: MPI rank of new process, or a vector of MPI ranks, all of which
              will be added

"""
function pushReceivePart!(data::ScatterData, peernum::Integer)

  peeridx = getRemotePeerIdx(data, peernum)
  if peeridx == 0
    push!(data.peernums_recv, peernum)
  end

  return nothing
end


function pushReceivePart!(data::ScatterData, peernum::AbstractArray)

  for p in peernum
    pushReceivePart!(data, p)
  end

  return nothing
end



"""
  Performs parallel communication to send the entities supplied to
  [`pushKey!`](@ref) to the remote processes.  The order the keys were
  supplied to `pushKey!` defines the order in which values must be to
  [`pushValue!`](@ref).

  **Inputs**

   * data: [`ScatterData`](@ref) object
"""
function exchangeKeys(data::ScatterData{T}) where {T}

  # in principle, the user could figure out how many values will be received
  # while calling pushKey!(), but we do some extra parallel communication here
  # to make things easier for them

  # first communication: figure out how many values will be received
  # second communication: send the MeshEntity* themselves
 
  if !haskey(MPI.mpitype_dict, Ptr{Void})
    @assert sizeof(Ptr{Void}) == sizeof(Int)
    MPI.mpitype_dict[Ptr{Void}] = MPI.mpitype_dict[Int]
  end

  # receive first communication
  tag = data.tag + 1
  tag2 = data.tag
  nrecvs = length(data.peernums_recv)
  recv_vals = Array{Array{Int, 1}}(nrecvs)
  recv_reqs = Array{MPI.Request}(nrecvs)
  for i=1:nrecvs
    recv_vals[i] = Array{Int}(1)
    recv_reqs[i] = MPI.Irecv!(recv_vals[i], data.peernums_recv[i], tag, data.comm)
  end

  # start first communication: because these are small messages, send them
  #       all before starting the second communication
  nsends = length(data.peernums_send)
  send_reqs = Array{MPI.Request}(nsends)
  for i=1:nsends
    nvals_i = length(data.send[i].entities)
    send_reqs[i] = MPI.Isend(nvals_i, data.peernums_send[i], tag, data.comm)
  end

  # start second communication: do this after all first communications have been
  # sent, because those will be done eagerly
  send_reqs2 = Array{MPI.Request}(nsends)
  for i=1:nsends
    data_i = data.send[i]
    send_reqs2[i] = MPI.Isend(data_i.entities, data.peernums_send[i], tag2, data.comm)
  end

  # wait for first communication and post receives for the second
  recv_reqs2 = Array{MPI.Request}(nrecvs)
  resize!(data.recv, nrecvs)
  for i=1:nrecvs
    idx, stat = MPI.Waitany!(recv_reqs)
    recv_reqs[idx] = MPI.REQUEST_NULL
    data.recv[idx] = PeerData(T, data.peernums_recv[idx], data.dims, recv_vals[idx][1], data.tag, data.comm)


    # post receive for second communication
    recv_reqs2[idx] = MPI.Irecv!(data.recv[idx].entities, data.peernums_recv[idx], tag2, data.comm)

  end

  # allocate value buffers.  Not needed in this function, but gives us
  # something to do while communication progresses
  for i=1:nsends
    allocateBuffer(data.send[i], data.dims)
  end



  MPI.Waitall!(recv_reqs2)
  MPI.Waitall!(send_reqs)
  MPI.Waitall!(send_reqs2)

  return nothing
end


"""
  Add a value to the send buffer.  Values must be added in the same order
  is the indices in `data.sends[peeridx].local_indices`

  **Inputs**

   * data: the [`ScatterData`](@ref) object
   * peeridx: the index (not MPI rank) of process ot send to
   * val: value to copy into receive buffer
"""
function pushValue!(data::ScatterData{T, 2}, peeridx::Integer, val::Number) where {T}

  idx = data.curridx[peeridx]
  data.send[peeridx].vals[idx] = val
  data.curridx[peeridx] = idx + 1

  return nothing
end


function pushValue!(data::ScatterData{T, N}, peeridx::Integer,
                    val::AbstractArray{T2, N2}) where {T, N, T2, N2}

  idx = data.curridx[peeridx]
  tpl = getDataColons(data)
  dest_view = sview(data.send[peeridx].vals, tpl..., idx)
  @assert size(dest_view) == size(val)
  for i=1:length(val)
    dest_view[i] = val[i]
  end
  data.curridx[peeridx] = idx + 1

  return nothing
end


"""
  Helper function for [`pushValue!`](@ref).  Returns a (type-stable) typle
  of `N` colons, where `N` is the dimensionality of the array

  **Inputs**

   * AbstractArray{T, N}: an abstract array.

  **Outputs**

   * tpl: a tuple of `N` colons
"""
function getDataColons(data::ScatterData{T, N, N2, N3})where {T, N, N2, N3}
  return ntuple(x -> :, Val{N3})
end


"""
  Resets the data structure to prepare for another communication cycle

  **Inputs**

   * data: [`ScatterData`](@ref)
"""
function resetBuffers(data::ScatterData)

  # make sure previous send has finished so it is safe to overwrite the send
  # buffer

  for i=1:length(data.peernums_send)
    data_i = data.send[i]
    Wait!(data_i)
    data.curridx[i] = 1
  end

  for i=1:length(data.peernums_recv)
    data_i = data.recv[i]
    Wait!(data_i)
  end

  return nothing
end


"""
  Wrapper around MPI.Waitany on a vector of `PeerData` objects

  **Inputs**

  * data: vector of [`PeerData`](@ref) objects

  **Outputs**

   * idx: index of `data.recv` object that was waited on
"""
function waitAny(data::Vector{P}) where {P <: PeerData}

  npeers = length(data)
  recv_reqs = Array{MPI.Request}(npeers)
  for i=1:npeers
    recv_reqs[i] = data[i].req
  end

  idx, stat = MPI.Waitany!(recv_reqs)

  data[idx].req = MPI.REQUEST_NULL
  data[idx].req_waited = true

  return idx
end


"""
  Extracts the data from the array and starts MPI communication

  **Inputs**

   * data: [`ScatterData`](@ref)
   * arr: array to extract data from.  Only the values that will be sent in
          parallel will be accessed
   * reduce_op: a [`Reduction`](@ref) to apply to values corresponding to
                the same mesh entity
"""
function sendParallelData(data::ScatterData{T, N, N2}, arr::AbstractArray{T2, N3}, reduce_op::Reduction{T}=SumReduction{T}()) where {T, T2, N, N2, N3}

  @assert N3 == (N + N2 - 1)

  resetBuffers(data)


  # post receives
  for data_i in data.recv
    Irecv!(data_i)
  end

  if length(data.send) == 0
    return nothing
  end

  tpl = getDataColons(data)

  # this is a trick to get a temporary array of the right size, because I 
  # can't figure out a type-stable way of getting the dimensions without
  # constructing an intermediate view
  idx = data.send[1].local_indices[1]
  _tmp_view = sview(arr, tpl..., idx.I...)
  vals_arr = Array{T}(size(_tmp_view))

  for peeridx=1:length(data.peernums_send)
    data_i = data.send[peeridx]
    for i=1:length(data_i.entities)
      # loop over all values associated with entity i
      fill!(vals_arr, reduce_op.neutral_element)
      for j=data_i.colptr[i]:(data_i.colptr[i+1]-1)
        idx = data_i.local_indices[j]

        vals = sview(arr, tpl..., idx.I...) # splat the colons, and the integer
                                            # indices inside the Cartesian
                                            # index
        # apply reduction
        for k=1:length(vals)
          vals_arr[k] = reduce_op(vals_arr[k], vals[k])
        end

      end  # end j
      # save result to buffer
      pushValue!(data, peeridx, vals_arr)

    end # end i

    Isend(data_i)
  end # end peeridx

  return nothing
end


"""
  Receives the parallel communication started by [`sendParallelData`](@ref).

  **Inputs**

   * data: [`ScatterData`](@ref) object
   * calc_func: function to be called with each set of received data.
                This function must have the signature:
                `calc_func(data::PeerData)`, where `data` is the
                `PeerData` object for the receive.
"""
function receiveParallelData(data::ScatterData, calc_func::Function)

  for i=1:length(data.recv)
    idx = waitAny(data.recv)
    calc_func(data.recv[idx])
  end

  return nothing
end


"""
  Function for receiving data shaped liked the coordinate field and putting
  it into a vector.  Note that this does *not* meet the interface defined by
  [`receiveParallelData`](@ref), so it will have to be called from inside a
  lambda function

  **Inputs**

   * data: [`ScatterData`](@ref) object
   * mesh: mesh object
   * reduce_op: a [`Reduction`](@ref) object to apply
   * vec: vector to put the result in
"""
function receiveVecFunction(data::PeerData{T, 2},
                            mesh::PumiMesh, vec::AbstractVector,
                            reduce_op::Reduction=SumReduction{T}()) where {T}

  pos = 1
  for i=1:length(data.entities)
    entity = data.entities[i]
    typ = getType(mesh.m_ptr, entity)
    dim = getTypeDimension(typ)
    for j=1:mesh.coord_numNodesPerType[dim+1]
      for k=1:mesh.dim
        idx = getNumberJ(mesh.coord_nodenums_Nptr, entity, j-1, k-1)
        vec[idx] = reduce_op(data.vals[k, pos], vec[idx])
      end
      pos += 1
    end
  end

  return nothing
end




#------------------------------------------------------------------------------
# Sample functions for 3D arrays


"""
  Initializes a [`ScatterData`] object for sending data in a 3D array format
  to the owning MeshEntity, which receives it in vector format.

  The 3D array *must* have the data in the second dimension ordered according to
  [`getNodeEntities`](@ref)

  **Inputs**

   * mesh
   * fshape: the apf::FieldShape of the data
   * dims: the tuple of dimension describing size of the data at each node of the
           field.

  **Outputs**

   * data: the `ScatterData` object
"""
function initSendToOwner(mesh::PumiMesh{T}, fshape::Ptr{Void}, dims::NTuple) where {T}

  data = ScatterData(T, dims, mesh.comm)
  shr = mesh.normalshr_ptr
#  shr = getNormalSharing(mesh.m_ptr)

  part_nums = Array{Cint}(400)  # apf::Up
  remote_entities = Array{Ptr{Void}}(400)  # apf::Up
  _adj_els = Array{Ptr{Void}}(400)  # apf::Up
  node_entities = ElementNodeEntities(mesh.m_ptr, fshape, mesh.dim)

  coords = Array{Float64}(3)  # DEBUGGING
  local_indices = Array{CartesianIndex{2}}(400)
  for dim=0:(mesh.dim-1)  # elements cant be shared in parallel
    it = MeshIterator(mesh.m_ptr, dim)
    nnodes = countNodesOn(fshape, dim)  # dim == apf::Type up to faces

    if nnodes == 0
      continue
    end

    for i=1:mesh.numEntitiesPerType[dim+1]

      entity = iterate(mesh.m_ptr, it)

      # bail out early if entity is not shared
      if !isSharedShr(shr, entity)
        continue
      end

      ncopy = countCopies(shr, entity)
      @assert ncopy <= 400
      getCopies(part_nums, remote_entities)

      owner = getOwner(shr, entity)
      if owner == mesh.myrank
        pushReceivePart!(data, sview(part_nums, 1:Int(ncopy)))
      else # figure out all entries in the 3D array that correspond to this
           # entity

        # get remote_entity
        owner_idx = findfirst(sview(part_nums, 1:Int(ncopy)), owner)
        remote_entity = remote_entities[owner_idx]
        

        nel = countAdjacent(mesh.m_ptr, entity, mesh.dim)
        @assert nel <= 400
        getAdjacent(_adj_els); adj_els = sview(_adj_els, 1:Int(nel))

        addEntityKeys(mesh, data, nnodes, node_entities, adj_els,
                      entity, remote_entity, owner, local_indices)

      end  # end if
    end  # end loop i

    free(mesh.m_ptr, it)
  end  #end loop dim


  freeSharing(shr)

  return data
end


"""
  Calls [`pushKey`](@ref) for every node on the given element (and computes
  the `local_indices` required)

  **Inputs**

   * mesh
   * data: [`ScatterData`](@ref) object
   * nnodes: number of nodes on the entity
   * node_entities: an [`ElementNodeEntities`](@ref) object (reusable storage)
   * adj_els: vector of elements that are upward adjacent of `entity`
   * entity: the local MeshEntity*
   * remote_entity: the MeshEntity* on the owner process
   * owner: the MPI rank of the owning process
   * local_indices: array to be overwritten with the `CartesianIndex`es of
                    of each (node, element) of the 3D array
"""
function addEntityKeys(mesh::PumiMesh, data::ScatterData, nnodes::Integer,
                       node_entities::ElementNodeEntities,
                       adj_els::AbstractVector, entity::Ptr{Void},
                       remote_entity::Ptr{Void}, owner::Integer,
                       local_indices::AbstractVector)

  #TODO: this may need to be updated for higher order meshes, to verify
  #      the shared entity has the same orientation on both processes
  for k=1:nnodes  # add the current entity nnodes times
    local_idx = 1  # reset index into local_indices
    for j=1:length(adj_els)
      getNodeEntities(node_entities, adj_els[j])
      e_idx = findfirst(node_entities.entities, entity)

      # compute index in second dimension of 3D array
      idx = k
      for p=1:(e_idx-1)
        idx += node_entities.nodecounts[p]
      end

      elnum = getNumberJ(mesh.el_Nptr, adj_els[j], 0, 0) + 1
      local_indices[local_idx] = CartesianIndex(idx, elnum)
      local_idx += 1
    end  # end loop j

    pushKey!(data, entity, remote_entity, owner,
             sview(local_indices, 1:(local_idx-1)))
  end  # end loop /

  return nothing
end


# this function is used by coords1DTo3D
"""
  This function is used to scatter the data by [`coords1DTo3D`](@ref).
  This function reverses the roles of `data.send` and `data.recv`

  **Inputs**

   * mesh
   * data: [`ScatterData`](@ref)
   * xvec: vector of coordinate data (indexed by `mesh.coord_nodenums_Nptr`)
           that will be scattered to all users of the data.  This is somewhat
           wasteful, because we only need to send it if this process is the
           owner.
"""
function sendParallelData_rev(mesh::PumiMesh, data::ScatterData, xvec::AbstractVector)

  resetBuffers(data)

  # post receives
  for data_i in data.send
    Irecv!(data_i)
  end

  for data_i in data.recv
    idx_dest = 1
    for entity in data_i.entities
      dim = getDimension(mesh.m_ptr, entity)
      for j=1:mesh.coord_numNodesPerType[dim+1]
        for k=1:mesh.dim
          idx_src = getNumberJ(mesh.coord_nodenums_Nptr, entity, j-1, k-1)
          data_i.vals[k, idx_dest] = xvec[idx_src]
        end
        idx_dest += 1
      end  # end j
    end  # end entity

    Isend(data_i)
  end

  return nothing
end


"""
  Function to receive data send by [`sendParallelData_rev`](@ref).
"""
function receiveParallelData_rev(data::ScatterData, calc_func::Function)

  for i=1:length(data.peernums_send)
    idx = waitAny(data.send)
    calc_func(data.send[idx])
  end

  return nothing
end


"""
  Function to help [`coords1DTo3D`](@ref) with reveiving data.  It only assigns
  the data to the vector if the remote process owns it

  **Inputs**

   * data: a [`PeerData`](@ref) object
   * mesh
  
  **Inputs/Outputs**

   * xvec: coordinate vector, indexed by `mesh.coord_nodenums_Nptr`.  Only
           entries corresponding to MeshEntities not owned by this process
           will be overwritten.  They will be overwritten by the data from the
           owning process.
"""
function receiveFromOwner(data::PeerData{T, 2}, mesh::PumiMesh, xvec::AbstractVector) where {T}

  idx_src = 1
  for entity in data._entities_local
    if data.peernum == getOwner(mesh.normalshr_ptr, entity)
      dim = getDimension(mesh.m_ptr, entity)
      for j=1:mesh.coord_numNodesPerType[dim+1]
        for k=1:mesh.dim
          idx_dest = getNumberJ(mesh.coord_nodenums_Nptr, entity, j-1, k-1)
          xvec[idx_dest] = data.vals[k, idx_src]
        end
        idx_src += 1
      end  # end j
    end  # end if
  end  # end entity

  return nothing
end



