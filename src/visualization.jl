# visualization related functions

# can be generalized with a few constants
function saveSolutionToMesh(mesh::PumiMesh, u::AbstractVector)
# saves the solution in the vector u to the mesh (in preparation for mesh adaptation
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

  if !(mesh.isDG && mesh.shape_type == 3)
    transferFieldToSubmesh(mesh)
  end

  return nothing
end  # end function saveSolutionToMesh

function saveSolutionToMesh(mesh::PumiMesh3DG, u::AbstractVector)
# interpolate solution to the vertices 
# u is the vector containing the solution, even though for DG the vector and the 
# 3D array forms of the solution have identical memory layout
# change this in the future with ReshapedArrays?

  interp = mesh.interp_op
  println("size(interp, 1) = ", size(interp, 1))
  u_el = zeros(Float64, mesh.numNodesPerElement, mesh.numDofPerNode)
  u_verts = zeros(Float64, size(interp, 1), mesh.numDofPerNode)
  u_node = zeros(Float64, mesh.numDofPerNode)
  dofs = mesh.dofs
  numEntitiesPerType = mesh.numTypePerElement
  fshape_ptr = mesh.coordshape_ptr

  println("size(u_verts) = ", size(u_verts))
  println("size(u_node) = ", size(u_node))

  # count nodes on solution field 
  numNodesPerType = Array(Int, 4)
  numNodesPerType[1] = countNodesOn(fshape_ptr, 0)
  numNodesPerType[2] = countNodesOn(fshape_ptr, 1)
  numNodesPerType[3] = countNodesOn(fshape_ptr, 2)
  numNodesPerType[4] = countNodesOn(fshape_ptr, 4) # tetrahedron
 
  for el=1:mesh.numEl
    el_i = mesh.elements[el]

    # get the solution values out of u
    for j=1:mesh.numNodesPerElement
      for k=1:mesh.numDofPerNode
        dof_k = dofs[k, j, el]
        u_el[j, k] = real(u[dof_k])
      end
    end

    # interpolate
    smallmatmat!(interp, u_el, u_verts)

    # save values to mesh
    # assumes the solution fieldshape is the same as the coordinate fieldshape
    node_entities = getNodeEntities(mesh.m_ptr, fshape_ptr, el_i)
    col = 1  # current node of the element

    for i=1:(mesh.dim+1)
      for j=1:numEntitiesPerType[i]
        for k=1:numNodesPerType[i]
          entity = node_entities[col]
          # skip elementNodeOffsets - maximum of 1 node per entity
          for p=1:mesh.numDofPerNode
            u_node[p] = u_verts[col, p]
          end

          setComponents(mesh.f_ptr, entity, 0, u_node)
          col += 1
        end
      end
    end  # end loop over entity dimensions

  end  # end loop over elements


  return nothing
end  # end function


function writeVisFiles(mesh::PumiMesh2DG, fname::AbstractString)
  # writes vtk files 

  if mesh.shape_type != 3
    writeVtkFiles(fname, mesh.mnew_ptr)
  else
    println(STDERR, "Warning: not printing visualization file for 2D SBP-Gamma")
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

function writeVisFiles(mesh::PumiMesh3DG, fname::AbstractString)
  writeVtkFiles(fname, mesh.m_ptr)
end



