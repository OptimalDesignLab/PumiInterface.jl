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

    println("data.peernums_send = ", data.peernums_send)
    println("data.peernums_recv = ", data.peernums_recv)

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
    println("data.send[2].vals = ", data.send[2].vals)
    pushValue!(data, 2, [4.0, 5.0])
    println("data.send[2].vals = ", data.send[2].vals)

    @test data.curridx[2] == 3
    println("data.send[2].vals = ", data.send[2].vals)
    @test maximum(abs.(data.send[2].vals[:, 1] - [2.0, 3.0])) < 1e-13
    @test maximum(abs.(data.send[2].vals[:, 2] - [4.0, 5.0])) < 1e-13

  end

  return nothing
end



