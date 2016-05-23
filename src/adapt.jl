function retrieveSolutionFromMesh(mesh::PumiMeshDG2, u::AbstractVector)
# retrieve solution from mesh (after mesh adaptation)
# mesh needs to have been reinitilized after mesh adaptation, u needs to be the right size for the new mesh




  num_entities = [3, 3, 1] # number of vertices, edges, faces

  q_vals = zeros(4)

  for el=1:mesh.numEl
    el_i = mesh.elements[el]
    node_entities = getNodeEntities(mesh.m_ptr, mesh.mshape_ptr, el_i)
    col = 1 # current node of the element
    for i=1:3  # loop over verts, edges, faces
      for j=1:num_entities[i]  # loop over all entities of this type
	for k=1:mesh.numNodesPerType[i]  # loop over nodes on this entity
	  entity = node_entities[col]  # current entity
	  offset_k = mesh.elementNodeOffsets[col, el]
         
	  # get values from mesh
          getComponents(mesh.f_ptr, entity, abs(offset_k - k) - 1, q_vals)

          # put solution values into vector u
	  for p=1:mesh.numDofPerNode  # loop over all dofs
	    dofnum_p = getNumberJ(mesh.dofnums_Nptr, entity, abs(offset_k - k) -1, p-1)
	    u[dofnum_p] = q_vals[p]
	  end
	  col += 1
	end  # end loop over nodes on curren entity
      end  # end loop over entities of current type
    end  # end loop over entity types
  end  # end loop over elements

  

  return nothing
end  # end function saveSolutionToMesh

#=
function retrieveNodeSolution(f_ptr, entity, u_node::AbstractVector)
# retrieve solution on a particular entity, stores it in u_node
# u_node must be a vector of length mesh.numDofPerNode
# used during mesh adaptation
# this is a low level function because t takes in f_ptr, entity rather than the mesh object and an index

  getComponents(f_ptr, entity, 0, u_node)

end
=#


function retrieveSolutionFromMesh(mesh::PumiMesh2, u::AbstractVector)
# retrieve solution from mesh (after mesh adaptation)
# mesh needs to have been reinitilized after mesh adaptation, u needs to be the right size for the new mesh




  num_entities = [3, 3, 1] # number of vertices, edges, faces

  q_vals = zeros(4)

  for el=1:mesh.numEl
    el_i = mesh.elements[el]
    node_entities = getNodeEntities(mesh.m_ptr, mesh.mshape_ptr, el_i)
    col = 1 # current node of the element
    for i=1:3  # loop over verts, edges, faces
      for j=1:num_entities[i]  # loop over all entities of this type
	for k=1:mesh.numNodesPerType[i]  # loop over nodes on this entity
	  entity = node_entities[col]  # current entity
	  offset_k = mesh.elementNodeOffsets[col, el]
         
	  # get values from mesh
          getComponents(mesh.f_ptr, entity, abs(offset_k - k) - 1, q_vals)

          # put solution values into vector u
	  for p=1:mesh.numDofPerNode  # loop over all dofs
	    dofnum_p = getNumberJ(mesh.dofnums_Nptr, entity, abs(offset_k - k) -1, p-1)
	    u[dofnum_p] = q_vals[p]
	  end
	  col += 1
	end  # end loop over nodes on curren entity
      end  # end loop over entities of current type
    end  # end loop over entity types
  end  # end loop over elements

  

  return nothing
end  # end function saveSolutionToMesh


function retrieveNodeSolution(f_ptr, entity, u_node::AbstractVector)
# retrieve solution on a particular entity, stores it in u_node
# u_node must be a vector of length mesh.numDofPerNode
# used during mesh adaptation
# this is a low level function because t takes in f_ptr, entity rather than the mesh object and an index

  getComponents(f_ptr, entity, 0, u_node)

end


