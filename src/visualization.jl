# visualization related functions

# can be generalized with a few constants
function _saveSolutionToMesh(mesh::PumiMesh, u::AbstractVector, 
                             mnew_ptr=mesh.mnew_ptr, fnew_ptr=mesh.fnew_ptr)
# saves the solution in the vector u to the mesh.f_ptr, and then
# callcs transferFieldToSubmesh to transfer the field fo 
# fnew_ptr on mesh mnew_ptr, if needed
# it uses mesh.elementNodeOffsets to access the pumi field values in the 
# right order of the given element
# note that this performs a loop over elements, so some values get
# written several times.  This is why is is necessary to write to the
# field in the right order
# because this function used getNumberJ to get dof numbers, and set 
# values, it doesn't really need to use the nodemap because values
# are set/get *consistently*, even if not in the same order for every
# element depending on the orientation
 # dofnums = zeros(Int, mesh.numDofPerNode)
#  u_vals = zeros(mesh.numDofPerNode)


  num_entities = [3, 3, 1] # number of vertices, edges, faces
  q_vals = zeros(mesh.numDofPerNode)

  for el=1:mesh.numEl
    el_i = mesh.elements[el]
    node_entities = getNodeEntities(mesh.m_ptr, mesh.mshape_ptr, el_i)
    col = 1 # current node of the element
    for i=1:3  # loop over verts, edges, faces
      for j=1:num_entities[i]  # loop over all entities of this type
	for k=1:mesh.numNodesPerType[i]  # loop over nodes on this entity
	  entity = node_entities[col]  # current entity
	  offset_k = mesh.elementNodeOffsets[col, el] # offset for current node
	  new_node = abs(offset_k - k) - 1

          # get solution values
	  for p=1:mesh.numDofPerNode  # loop over all dofs
	    dofnum_p = getNumberJ(mesh.dofnums_Nptr, entity, new_node, p-1)
	    q_vals[p] = u[dofnum_p]
	  end
          # save to mesh
          setComponents(mesh.f_ptr, entity, abs(offset_k - k) - 1, q_vals)

	  col += 1
	end  # end loop over nodes on curren entity
      end  # end loop over entities of current type
    end  # end loop over entity types
  end  # end loop over elements

  transferFieldToSubmesh(mesh, u, mnew_ptr, fnew_ptr)

  return nothing
end  # end function saveSolutionToMesh

function saveSolutionToMesh(mesh::PumiMesh, u::AbstractVector)
  if mesh.isDG
    # all DG meshes interpolate directly
    interpolateToMesh(mesh, u)
  else
    _saveSolutionToMesh(mesh, u)
  end

  if mesh.isDG && mesh.mexact_ptr != C_NULL
    _saveSolutionToMesh(mesh, u, mesh.mexact_ptr, mesh.fexact_ptr)
  end
end

function interpolateToMesh{T}(mesh::PumiMesh{T}, u::AbstractVector)
# interpolate solution stored in u to the field mesh.fnew_ptr which resides on
# mesh.mnew_ptr 
# u is the vector containing the solution, even though for DG the vector and 
# the 3D array forms of the solution have identical memory layout
# change this in the future with ReshapedArrays?

# the interpolation operator *must* interpolate the solution in a particular
# order.  The order must be the order of the nodes returned by 
# getNodeEntities(mnew_ptr)

  @assert mesh.m_ptr == mesh.mnew_ptr  # needed for averaging step

  interp = mesh.interp_op
  u_el = zeros(Float64, mesh.numNodesPerElement, mesh.numDofPerNode)
  u_verts = zeros(Float64, size(interp, 1), mesh.numDofPerNode)
  jac_verts = zeros(T, size(interp, 1))
  u_node = zeros(Float64, mesh.numDofPerNode)  # hold new node values
  u_node2 = zeros(u_node)  # hold existing node values
  node_entities = Array(Ptr{Void}, mesh.numNodesPerElement)
  dofs = mesh.dofs
  numTypePerElement = mesh.numTypePerElement
  numEntitiesPerType = mesh.numEntitiesPerType
  fshape_ptr = mesh.fnewshape_ptr
  vshare = mesh.vert_sharing

  # count nodes on solution field 
  numNodesPerType = Array(Int, mesh.dim + 1)
  numNodesPerType[1] = countNodesOn(fshape_ptr, 0)
  numNodesPerType[2] = countNodesOn(fshape_ptr, 1)
  numNodesPerType[3] = countNodesOn(fshape_ptr, 2)
  if mesh.dim == 3
    numNodesPerType[4] = countNodesOn(fshape_ptr, 4) # tetrahedron
  end

  zeroField(mesh.fnew_ptr)

  # data for MPI
  # values solution values in first n indices of inner array, followed
  # by the weighting factor
  peer_vals_send = Array(Array{Float64, 2}, vshare.npeers)
  peer_vals_recv = Array(Array{Float64, 2}, vshare.npeers)
  for i=1:vshare.npeers
    peer_vals_send[i] = zeros(Float64, mesh.numDofPerNode + 1, vshare.counts[i])
    peer_vals_recv[i] = zeros(peer_vals_send[i])

  end



  for el=1:mesh.numEl
    el_i = mesh.elements[el]

    # get the solution values out of u
    for j=1:mesh.numNodesPerElement
      for k=1:mesh.numDofPerNode
        dof_k = dofs[k, j, el]
        u_el[j, k] = real(u[dof_k])
      end
    end

    # interpolate solution
    smallmatmat!(interp, u_el, u_verts)

    # interpolate jacobian
    jac_el = sview(mesh.jac, :, el)
    smallmatvec!(interp, jac_el, jac_verts)


    # save values to mesh
    # assumes the solution fieldshape is the same as the coordinate fieldshape
    getNodeEntities(mesh.mnew_ptr, fshape_ptr, el_i, node_entities)
    col = 1  # current node of the element

    for i=1:(mesh.dim+1)
      for j=1:numTypePerElement[i]
        for k=1:numNodesPerType[i]
          entity = node_entities[col]
          @assert i == 1
          vertnum = getNumberJ(mesh.vert_Nptr, entity, 0, 0) + 1

          # pack array to send to other processes
          if haskey(vshare.rev_mapping, vertnum)
            pair = vshare.rev_mapping[vertnum]
            for peer=1:length(pair.first)
              shared_vert_idx = pair.second[peer]
              arr_peer = peer_vals_send[peer]
              # accumulate weighted solution values, weighting factor

              for p=1:mesh.numDofPerNode
                arr_peer[p, shared_vert_idx] += real(jac_verts[col])*u_verts[col, p]
              end
              arr_peer[mesh.numDofPerNode + 1, shared_vert_idx] += real(jac_verts[col])
            end  # end loop over peers
          end  # end if haskey


          # skip elementNodeOffsets - maximum of 1 node per entity
          getComponents(mesh.fnew_ptr, entity, 0, u_node2)
          for p=1:mesh.numDofPerNode
            # multiply by element volume as a weighting factor
            u_node[p] = real(jac_verts[col])*u_verts[col, p] + u_node2[p]
          end

          setComponents(mesh.fnew_ptr, entity, 0, u_node)
          col += 1
        end
      end
    end  # end loop over entity dimensions

  end  # end loop over elements

  # send the data
  send_reqs = Array(MPI.Request, vshare.npeers)
  recv_reqs = Array(MPI.Request, vshare.npeers)
  for i=1:vshare.npeers
    # use the tag 20 in case any other parts of the program have active
    # communications
    send_reqs[i] = MPI.Isend(peer_vals_send[i], vshare.peer_nums[i], 20, mesh.comm)
    recv_reqs[i] = MPI.Irecv!(peer_vals_recv[i], vshare.peer_nums[i], 20, mesh.comm)
  end

  # wait for all recieves to finish
  for i=1:vshare.npeers
    MPI.Wait!(recv_reqs[i])
  end

  # add in the contributions from each peer to the weighted sum
  #TODO: use MPI.WaitAny! to combine this with the above loop
  for peer=1:vshare.npeers
    peer_vals_p = peer_vals_recv[peer]
    vertnums_p = vshare.vert_nums[peer]

    for i=1:length(vertnums_p)
      vert_i = mesh.verts[vertnums_p[i]]
      getComponents(mesh.fnew_ptr, vert_i, 0, u_node)
      weight_i = peer_vals_p[mesh.numDofPerNode + 1, i]
      for p=1:mesh.numDofPerNode
        u_node[p] += weight_i*peer_vals_p[p, i]
      end
      setComponents(mesh.fnew_ptr, vert_i, 0, u_node)
    end
  end





  up_els = Array(Ptr{Void}, 400)  # equivalent of apf::Up
  # divide by the total volume of elements that contributed to each node
  # so the result is the average value
  for dim=1:(mesh.dim + 1)
    if numNodesPerType[dim] > 0
      @assert dim == 1
      resetIt(dim - 1) 
      for i=1:numEntitiesPerType[dim]
        entity_i = getEntity(dim - 1)
        entity_num = getNumberJ(mesh.vert_Nptr, entity_i, 0, 0) + 1
        nel = countAdjacent(mesh.mnew_ptr, entity_i, mesh.dim)
        getAdjacent(up_els)
        # compute the sum of the volumes of the elements
        jac_sum = 0.0
        for j=1:nel
          el_j = up_els[j]
          elnum_j = getNumberJ(mesh.el_Nptr, el_j, 0, 0) + 1

          jac_el = sview(mesh.jac, :, elnum_j)
          # searching getNodeEntities for the index guarantees we can
          # identify the right interpolated point for all elements, 
          # independent of topology
          getNodeEntities(mesh.m_ptr, fshape_ptr, el_j, node_entities)
          node_idx = findfirst(node_entities, entity_i)

          # interpolate the jacobian to this point
          jac_entity = 0.0
          for p=1:size(interp, 2)
            jac_entity += interp[node_idx, p]*real(jac_el[p])
          end

          jac_sum += jac_entity
        end  # end loop j

        # add in the parallel contributions
        if haskey(vshare.rev_mapping, entity_num)
          pair = vshare.rev_mapping[entity_num]
          for i=1:length(pair.first)  # loop over peers that share this vert
            peer_idx = pair.first[i]
            vert_idx = pair.second[i]
            jac_sum += peer_vals_recv[peer_idx][mesh.numDofPerNode + 1, vert_idx]
          end
        end  # end if haskey

        getComponents(mesh.fnew_ptr, entity_i, 0, u_node)
        fac = 1/jac_sum
        for p=1:mesh.numDofPerNode
          u_node[p] = fac*u_node[p]
        end
        setComponents(mesh.fnew_ptr, entity_i, 0, u_node)
        incrementIt(dim - 1)
      end  # end loop i
    end  # end if
  end  # end loop dim

  # wait for sends to finish before exiting
  for i=1:vshare.npeers
    MPI.Wait!(send_reqs)
  end


  return nothing
end  # end function


function writeVisFiles(mesh::PumiMeshDG, fname::AbstractString)
  # writes vtk files 

  writeVtkFiles(fname, mesh.mnew_ptr)

  if mesh.mexact_ptr != C_NULL
    fname_exact = fname*"_exact"
    writeVtkFiles(fname_exact, mesh.mexact_ptr)
  end

  return nothing
end


function writeVisFiles(mesh::PumiMesh2CG, fname::AbstractString)
  # writes vtk files 

  if mesh.order <= 2
    println("writing original mesh vtk file")
    writeVtkFiles(fname, mesh.m_ptr)
  else
    println("writing subtriangulated mesh vtk file")
    writeVtkFiles(fname, mesh.mnew_ptr)
  end

  return nothing
end
