# test parallel communication types

# import the API to make things easier
import PdePumiInterface: PeerData, ScatterData, allocateBuffer, pushKey!,
                         getPeerIdx, getRemotePeerIdx,
                         pushReceivePart!, exchangeKeys, pushValue!,
                         getDataColons, resetBuffers,
                         sendParallelData, receiveParallelData

"""
  Test the API in serial
"""
function test_ScatterData(mesh::PumiMesh{T}) where {T}

  @testset "ScatterData" begin
    data = ScatterData(T, (mesh.dim,), mesh.comm)

    # test the pushKey! API
    @test getPeerIdx(data, 1) == 0
    @test getRemotePeerIdx(data, 1) == 0

    idx = CartesianIndex(2, 3)
    peernum = 2
    e_local = C_NULL
    e_remote = C_NULL
    pushKey!(data, e_local, e_remote, peernum, idx)

    @test getPeerIdx(data, 1) == 0
    @test getRemotePeerIdx(data, 2) == 0
    @test getPeerIdx(data, 2) == 1
    @test getRemotePeerIdx(data, 2) == 0

    pushReceivePart!(data, 3)
    @test getPeerIdx(data, 1) == 0
    @test getRemotePeerIdx(data, 3) == 1
    @test getPeerIdx(data, 2) == 1
    @test getRemotePeerIdx(data, 2) == 0

    peernum = 4
    pushKey!(data, e_local, e_remote, peernum, idx)
    pushReceivePart!(data, 5)

    @test getPeerIdx(data, 1) == 0
    @test getRemotePeerIdx(data, 3) == 1
    @test getPeerIdx(data, 2) == 1
    @test getRemotePeerIdx(data, 2) == 0

    @test getPeerIdx(data, 4) == 2
    @test getRemotePeerIdx(data, 5) == 2

    pushKey!(data, e_local, e_remote, 4, [CartesianIndex(2, 4), CartesianIndex(3, 5)])
    @test length(data.send[2].local_indices) == 3

    
    for data_i in data.send
      allocateBuffer(data_i, data.dims)
    end


    # test pushValue!
    pushValue!(data, 2, [2.0, 3.0])
    pushValue!(data, 2, [4.0, 5.0])

    @test data.curridx[2] == 3
    @test maximum(abs.(data.send[2].vals[:, 1] - [2.0, 3.0])) < 1e-13
    @test maximum(abs.(data.send[2].vals[:, 2] - [4.0, 5.0])) < 1e-13

  end

  return nothing
end


"""
  Test `initSendToOwner` is correct by testing against the `mesh.vert_coords`
  array (parallel only)
"""
function test_initSendToOwner(mesh::PumiMesh{T}) where {T}

  data_recv = zeros(T, mesh.coord_numNodes*mesh.dim)
  # create lambda function for receiving data
  function calc_func(data::PdePumiInterface.PeerData)
    PdePumiInterface.receiveVecFunction(data, mesh, data_recv)
  end

  @testset "testing initSendToOwner" begin
    data = PdePumiInterface.initSendToOwner(mesh, mesh.coordshape_ptr, (mesh.dim,))

    shr = apf.getNormalSharing(mesh.m_ptr)
    # test that all entities owned by another process are present
    for dim=0:(mesh.dim-1)
      it = apf.MeshIterator(mesh.m_ptr, dim)
      nnodes = apf.countNodesOn(mesh.coordshape_ptr, dim) != 0


      for j=1:mesh.numEntitiesPerType[dim+1]
        entity = apf.iterate(mesh.m_ptr, it)
        owner = apf.getOwner(shr, entity)
        if apf.isSharedShr(shr, entity) && owner != mesh.myrank && apf.countNodesOn(mesh.coordshape_ptr, dim) != 0
          owner_idx = getPeerIdx(data, owner)
          @test entity in data.send[owner_idx]._entities_local
        end
      end
      apf.free(mesh.m_ptr, it)
    end

    # test that the local indices were computed correctly
    coords = Array{Float64}(3)
    for data_i in data.send
      for j=1:length(data_i._entities_local)
        apf.getPoint(mesh.m_ptr, data_i._entities_local[j], 0, coords)

        for j=data_i.colptr[j]:(data_i.colptr[j+1]-1)
          idx = data_i.local_indices[j]
          @test maximum(abs.(coords - mesh.vert_coords[:, idx])) < 1e-13
        end
      end
    end

    # test that sending the data works correctly
    exchangeKeys(data)

    # check that all entities received are owned by this process
    idx = 1
    for data_i in data.recv
      for entity in data_i.entities
        @test apf.isSharedShr(shr, entity)
        @test apf.getOwner(shr, entity) == mesh.myrank
      end
    end

    # send the coordinates using a SumReduction, so the receiver will get the
    # n * coordinates, where n is the number of remotes + self
    reduce_op = PdePumiInterface.AssignReduction{T}()
    sendParallelData(data, mesh.vert_coords, reduce_op)

    # also test coords3DTo1D at same time
    coords_vec = zeros(T, mesh.dim*mesh.coord_numNodes)
    coords3DTo1D(mesh, mesh.vert_coords, coords_vec, reduce_op)

    receiveParallelData(data, calc_func)
    coords = Array{Float64}(3)
    for dim=0:mesh.dim

      if !apf.hasNodesIn(mesh.coordshape_ptr, dim)
        continue
      end

      it = apf.MeshIterator(mesh.m_ptr, dim)
      for i=1:mesh.numEntitiesPerType[dim+1]
        entity = apf.iterate(mesh.m_ptr, it)
        owner = apf.getOwner(shr, entity)
        #typ = apf.getType(mesh.m_ptr, entity)

        if apf.isSharedShr(shr, entity)
          ncopies = apf.countCopies(shr, entity)
          nlocals = apf.countAdjacent(mesh.m_ptr, entity, mesh.dim)
          for j=1:mesh.coord_numNodesPerType[dim+1]
            apf.getPoint(mesh.m_ptr, entity, j-1, coords)
            for k=1:mesh.dim
              idx = apf.getNumberJ(mesh.coord_nodenums_Nptr, entity, j-1, k-1)
              if owner == mesh.myrank
                @test abs(data_recv[idx] - ncopies*coords[k]) < 1e-13
                @test abs(coords_vec[idx] - coords[k]) < 1e-13
              elseif owner != mesh.myrank
                idx = apf.getNumberJ(mesh.coord_nodenums_Nptr, entity, j-1, k-1)
                @test coords_vec[idx] == 0
              end  # end if
            end  # end k
          end  # end j
        end  # end if
      end  # end i
      apf.free(mesh.m_ptr, it)
    end  # end dim


    # test scattering data
    coords = zeros(Float64, 3)
    coords_vec = zeros(T, mesh.dim*mesh.coord_numNodes)
    coords3DTo1D(mesh, mesh.vert_coords, coords_vec, reduce_op)
    for i=1:length(coords_vec)  
      coords_vec[i] += mesh.myrank  # add rank-specific offset
    end
    coords_vec2 = copy(coords_vec)

    PdePumiInterface.sendParallelData_rev(mesh, data, coords_vec)
    calc_func2 = (data::PeerData) -> PdePumiInterface.receiveFromOwner(data, mesh, coords_vec2)
    PdePumiInterface.receiveParallelData_rev(data, calc_func2)
    # check the rank-specific offset
    for dim=0:mesh.dim
      it = apf.MeshIterator(mesh.m_ptr, dim)
      for i=1:mesh.numEntitiesPerType[dim+1]
        entity = apf.iterate(mesh.m_ptr, it)
        for j=1:mesh.coord_numNodesPerType[dim+1]
          apf.getPoint(mesh.m_ptr, entity, j-1, coords)
          for k=1:mesh.dim
            idx = apf.getNumberJ(mesh.coord_nodenums_Nptr, entity, j-1, k-1)

            if apf.isOwned(mesh.normalshr_ptr, entity)
              @test abs(coords_vec2[idx] - (coords[k] + mesh.myrank)) < 1e-13
            else
              owner = apf.getOwner(mesh.normalshr_ptr, entity)
              @test abs(coords_vec2[idx] - (coords[k] + owner)) < 1e-13
            end
          end  # end k
        end  # end j
      end   # end i

      apf.free(mesh.m_ptr, it)
    end  # end dim


  end  # end testset
#  error("stop here")

  return nothing
end

