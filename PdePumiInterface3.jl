export PumiMesh3

type PumiMesh3{T1} <: AbstractMesh{T1}   # 3d pumi mesh, tetrahedron only
  m_ptr::Ptr{Void}  # pointer to mesh
  mshape_ptr::Ptr{Void} # pointer to mesh's FieldShape
  f_ptr::Ptr{Void} # pointer to apf::field for storing solution during mesh adaptation
  vert_Nptr::Ptr{Void}  # numbering of vertices
  edge_Nptr::Ptr{Void}  # numbering of edges
  face_Nptr::Ptr{Void} # numbering of faces
  el_Nptr::Ptr{Void}  # numbering of elements (faces)

  numVert::Int  # number of vertices in the mesh
  numEdge::Int # number of edges in the mesh
  numFace::Int # number of faces in the mesh
  numEl::Int  # number of elements (regions)
  order::Int # order of shape functions
  numDof::Int # number of degrees of freedom
  numNodes::Int  # number of nodes
  numDofPerNode::Int  # number of dofs per node
  numBoundaryFaces::Int # number of faces on the exterior boundary
  numInterfaces::Int # number of internal interfaces
  numNodesPerElement::Int  # number of nodes per element
  numNodesPerType::Array{Int, 1}  # number of nodes classified on each vertex, edge, face, length 4
  numEntities::Array{Int, 1}  # number of entities of each type in the mesh, length 4
  numEntitiesPerElement::Array{Int, 1}  # number of vertices, edges, and regions in each element

  # hold pointers to mesh entities
  verts::Array{Ptr{Void},1}  # holds pointers to mesh vertices
  edges::Array{Ptr{Void},1}  # pointers to mesh edges
  faces::Array{Ptr{Void},1}  # pointers to mesh faces
  elements::Array{Ptr{Void},1}  # pointers to regions

  dofnums_Nptr::Ptr{Void}  # pointer to Numbering of dofs (result of reordering)
  boundary_nums::Array{Int, 2}  # array of [element number, edgenumber] for each edge on the boundary

  bndryfaces::Array{Boundary, 1}  # store data on external boundary of mesh
  interfaces::Array{Interface, 1}  # store data on internal edges

  coords::Array{T1, 3}  # store coordinates of all nodes
  dxidx::Array{T1, 4}  # store scaled mapping jacobian
  jac::Array{T1,2}  # store mapping jacobian output

  dofs::Array{Int, 3}  # store dof numbers of solution array to speed assembly



 function PumiMesh3(dmg_name::AbstractString, smb_name::AbstractString, order, sbp::SBPOperator; dofpernode=1)
  # construct pumi mesh by loading the files named
  # dmg_name = name of .dmg (geometry) file to load (use .null to load no file)
  # smb_name = name of .smb (mesh) file to load
  # order = order of shape functions
  # dofpernode = number of dof per node, default = 1

  mesh = new()
  mesh.numDofPerNode = dofpernode
  mesh.order = order
  mesh.numNodesPerElement = getNumNodes(order)
  mesh.numEntities, mesh.m_ptr, mesh.mshape_ptr = init(dmg_name, smb_name, order)
  mesh.f_ptr = createPackedField(mesh.m_ptr, "solution_field", dofpernode)

  mesh.numVert = convert(Int, mesh.numEntities[1])
  mesh.numEdge = convert(Int,  mesh.numEntities[2])
  mesh.numFace = convert(Int, mesh.numEntities[3])
  mesh.numEl = convert(Int, mesh.numEntities[4])

  # get number of nodes on each type
  mesh.numNodesPerType = zeros(Int, 4)
  mesh.numNodesPerType[1] = countNodesOn(mesh.mshape_ptr, 0)  # number of nodes on a vertex
  mesh.numNodesPerType[2]= countNodesOn(mesh.mshape_ptr, 1) # on edge
  mesh.numNodesPerType[3] = countNodesOn(mesh.mshape_ptr, 2) # on face
  mesh.numNodesPerType[4]= countNodesOn(mesh.mshape_ptr, 4) # on element interior 
  mesh.numEntitiesPerElement = [4, 6, 4, 1]

  mesh.numNodesPerElement = 0
  for i=1:4
    mesh.numNodesPerElement += mesh.numNodesPerType[i]*mesh.numEntitiesPerElement[i]
  end


  println("numEntities = ", mesh.numEntities)
  println("numNodesPerType = ", mesh.numNodesPerType)




  # get pointers to mesh entity numberings
  mesh.vert_Nptr = getVertNumbering()
  mesh.edge_Nptr = getEdgeNumbering()
  mesh.face_Nptr = getFaceNumbering()
  mesh.el_Nptr = getElNumbering()

  mesh.verts = Array(Ptr{Void}, mesh.numVert)
  mesh.edges = Array(Ptr{Void}, mesh.numEdge)
  mesh.faces = Array(Ptr{Void}, mesh.numFace)
  mesh.elements = Array(Ptr{Void}, mesh.numEl)
  mesh.dofnums_Nptr = createNumberingJ(mesh.m_ptr, "reordered dof numbers", mesh.mshape_ptr, dofpernode)  # 1 dof per node

  # get pointers to all MeshEntities
  # also initilize the field to zero
  resetAllIts2()
#  comps = zeros(dofpernode)
  for i=1:mesh.numVert
    mesh.verts[i] = getVert()
    incrementVertIt()
  end

  for i=1:mesh.numEdge
    mesh.edges[i] = getEdge()
    incrementEdgeIt()
  end

  println("getting Face pointers")
  for i=1:mesh.numFace
    mesh.faces[i] = getFace()
    incrementFaceIt()
  end

  println("getting region pointers")
  for i=1:mesh.numEl
    mesh.elements[i] = getEl()
    incrementElIt()
  end

  # use partially constructed mesh object to populate arrays

  
#  numberDofs(mesh)

  countBoundaryFaces(mesh)
   mesh.bndryfaces = Array(Boundary, mesh.numBoundaryFaces)
  getBoundaryArray(mesh)

  numberDofs(mesh)
#=
  mesh.bndryfaces = Array(Boundary, mesh.numBoundaryEdges)
  getBoundaryArray(mesh)

  mesh.numInterfaces = mesh.numEdge - mesh.numBoundaryEdges
  mesh.interfaces = Array(Interface, mesh.numInterfaces)
  getInterfaceArray(mesh)

  getCoordinates(mesh, sbp)  # store coordinates of all nodes into array
  getDofNumbers(mesh)  # store dof numbers

  mesh.dxidx = Array(T1, 2, 2, sbp.numnodes, mesh.numEl)
  mesh.jac = Array(T1, sbp.numnodes, mesh.numEl)
  mappingjacobian!(sbp, mesh.coords, mesh.dxidx, mesh.jac)
=#
#=
  println("typeof m_ptr = ", typeof(m_ptr))
  println("typeof mshape_ptr = ", typeof(mshape_ptr))
  println("typeof numVerg = ", typeof(numVert))
  println("typeof numEdge = ", typeof(numEdge))
  println("typeof numEl = ", typeof(numEl))
  println("typeof order = ", typeof(order))
  println("typeof numdof = ", typeof(numdof))
  println("typeof bnd_edges_cnt = ", typeof(bnd_edges_cnt))
  println("typeof verts = ", typeof(verts))
  println("typeof edges = ", typeof(edges))
  println("typeof element = ", typeof(elements))
  println("typeof dofnums_Nptr = ", typeof(dofnums_Nptr))
  println("typeof bnd_edges_small = ", typeof(bnd_edges_small))
=#
  println("numVert = ", mesh.numVert)
  println("numEdge = ", mesh.numEdge)
  println("numEl = ", mesh.numEl)
  println("numDof = ", mesh.numDof)
  println("numNodes = ", mesh.numNodes)



  writeVtkFiles("mesh_complete", mesh.m_ptr)
  return mesh
  # could use incomplete initilization to avoid copying arrays
#  return PumiMesh2(m_ptr, mshape_ptr, f_ptr, vert_Nptr, edge_Nptr, el_Nptr, numVert, numEdge, numEl, order, numdof, numnodes, dofpernode, bnd_edges_cnt, verts, edges, elements, dofnums_Nptr, bnd_edges_small)
  end

 
  
end

function getNumNodes(order::Integer)
# get the number of nodes on an element
# assumes SBP elements

  nnodes = 0
  if order == 1
    nnodes = 4
  elseif order == 2
    nnodes = 11
  elseif order == 3
    nnodes = 24
  elseif order == 4
    nnodes = 45
  else
    println("Warning, unsupported order elements are requested")
  end

  return nnodes
end




function countBoundaryFaces(mesh::PumiMesh3)
  # count boundary face
  # store array of [element number, global edge number]
  resetFaceIt()
  bnd_faces_cnt = 0
  bnd_faces = Array(Int, mesh.numFace, 2)
  for i=1:mesh.numFace
    face_i = getFace()
    numRegion = countAdjacent(mesh.m_ptr, face_i, 3)  # should be count upward

    if numRegion == 1  # if an exterior edge
      elements = getAdjacent(numRegion)
      elnum = getElNumber2(elements[1]) + 1

      println("face number ", i, " is part of only one element")

      bnd_faces_cnt += 1
      bnd_faces[bnd_faces_cnt, 1] = elnum
      bnd_faces[bnd_faces_cnt, 2] = i
    end
    incrementFaceIt()
  end


  mesh.boundary_nums = bnd_faces[1:bnd_faces_cnt, :] # copy, bad but unavoidable
  mesh.numBoundaryFaces = bnd_faces_cnt

return nothing

end  # end function


function getBoundaryArray(mesh::PumiMesh3)
# get an array of type Boundary for SBP
# creating an an array of a user defined type seems like a waste of memory operations
# bnd_array is a vector of type Boundary with length equal to the number of edges on the boundary of the mesh

#  mesh.bndryfaces = Array(Boundary, mesh.numBoundaryEdges)

  for i=1:mesh.numBoundaryFaces
    elnum = mesh.boundary_nums[i,1]
    facenum_global = mesh.boundary_nums[i,2]
    facenum_local = getBoundaryEntityLocalNum(mesh, facenum_global, 2)
    mesh.bndryfaces[i] = Boundary(elnum, facenum_local)
  end

  return nothing
end

function getBoundaryEntityLocalNum(mesh::PumiMesh3, entity_num::Integer, entity_dim::Integer)
# gets the local edge number of a specified edge that is on the boundary
# of the mesh
# entity_num is an the  num from the output of getBoundaryEdgeNums() (ie. the global entity number)
# entity_dim is the dimension of the entity (0 = vertex, 1 = edge, 2 = face)
# the local edge number is the edge number within the element (1st,s 2nd, 3rd...)
  if entity_dim == 1
    entity_i = mesh.edges[entity_num]
  elseif entity_dim == 2
    entity_i = mesh.faces[entity_num]
  else
    println("unsupported entity dimension requested for getBoundaryEntityLocalNum")
    println("requested entity dimension = ", entity_dim)
  end

  # get mesh face or region associated with edge
  countAdjacent(mesh.m_ptr, entity_i, entity_dim+1)
  parent = getAdjacent(1)[1]  # get the single face (not an array)

  local_num = getEntityLocalNumber(mesh.m_ptr, entity_i, parent, entity_dim, entity_dim+1)

  return local_num

end

function numberDofs(mesh::PumiMesh3)
# assign dof numbers to entire mesh
# calculates numNodes, numDof
# assumes mesh elements have already been reordered
# this method minimizes the assembly/disassembly time
# because blocks of entries in res go into res_vec
  println("Entered numberDofs")

 
  # calculate total number of nodes
  # Continuous Galerkin only
  numnodes = 0
  for i=1:4
    numnodes += mesh.numNodesPerType[i]*mesh.numEntities[i]
  end

  numDof = mesh.numDofPerNode*numnodes
 
  println("expected number of dofs = ", numDof)

  if (numDof > 2^30)
    println("Warning: too many dofs, renumbering will fail")
  end

  # save numbers to mesh
  mesh.numNodes = numnodes
  mesh.numDof = numDof

  # initally number all dofs as numDof+1 to 2*numDof
  # this allows quick check to see if somthing is labelled or not

  resetAllIts2()
  # mesh iterator increment, retreval functions
  iterators_inc = [incrementVertIt, incrementEdgeIt, incrementFaceIt, incrementElIt]
  iterators_get = [getVert, getEdge, getFace, getEl]
  num_entities = mesh.numEntities
  num_nodes_entity = mesh.numNodesPerType

  println("num_nodes_entity = ", num_nodes_entity)
  println("num_entities = ", num_entities)
  #  curr_dof = 1
  curr_dof = numDof+1
  for etype = 1:length(num_nodes_entity) # loop over entity types
    println("etype = ", etype)
    if (num_nodes_entity[etype] != 0)  # if no nodes on this type of entity, skip
      for entity = 1:num_entities[etype]  # loop over all entities of this type
	println("  entity number: ", entity)
	entity_ptr = iterators_get[etype]()  # get entity

	for node = 1:num_nodes_entity[etype]
	  println("    node : ", node)
	  for dof = 1:mesh.numDofPerNode
	    numberJ(mesh.dofnums_Nptr, entity_ptr, node-1, dof-1, curr_dof)
#	    println("      entity ", entity_ptr, " labelled ", curr_dof)
	    curr_dof += 1
	  end  # end loop over dof
	end  # end loop over node

      iterators_inc[etype]()
      end  # end loop over entitiesa
    end  # end if 
  end  # end loop over entity types

  println("Finished initial numbering of dofs") 


  println("Performing final numbering of dofs")
# move all if statements out one for loop (check only first dof on each node)
  curr_dof = 1
  for i=1:mesh.numEl
    println("element number: ", i)

    el_i = mesh.elements[i]
    # get vertices, edges for this element
    (verts_i, numVert) = getDownward(mesh.m_ptr, mesh.elements[i], 0)
#    println("verts_i = ", verts_i)
#    println("mesh.verts = ", mesh.verts)
    (edges_i, numEdge) = getDownward(mesh.m_ptr, mesh.elements[i], 1)

    if mesh.numEntitiesPerElement[3] > 1  # don't call get downward if el_i is a face
      (face_i, numFace) = getDownward(mesh.m_ptr, mesh.elements[i], 2)
      regions_i = [mesh.elements[i]]  # make this an array for consistency
      numRegions = 1
    else  # this is a 2d mesh
      face_i = [mesh.elements[i]]
      numFace = 1
      numRegion = 0
    end

    # loop over vertices
    for j=1:numVert  # loop over vertices, edges
      println("  vertex number: ", j)
      vert_j = verts_i[j]
#      println("  vert_j = ", vert_j)
#      println("  edge_j = ", edge_j)
      for k=1:mesh.numDofPerNode # loop over vertex dofs
#	println("    dof number: ", k)
#	println("    vert_j = ", vert_j)
        dofnum_k = getNumberJ(mesh.dofnums_Nptr, vert_j, 0, k-1)
	if dofnum_k > numDof  # still has initial number
	  # give it new (final) number
	  numberJ(mesh.dofnums_Nptr, vert_j, 0, k-1, curr_dof)
	  curr_dof += 1
	end
      end  # end loop over vertex dofs
    end  # end loop over vertices

    # now loop over edges
    for j=1:numEdge
      println("edge number: ", j)
      edge_j = edges_i[j]
      # loop over nodes on edge
      for k=1:num_nodes_entity[2]  # loop over all nodes on the edge
	for p=1:mesh.numDofPerNode  # loop over dofs
	  dofnum_p = getNumberJ(mesh.dofnums_Nptr, edge_j, k-1, p-1)
	  if dofnum_p > numDof  # still has initial number
	    # give it new (final) number
	    numberJ(mesh.dofnums_Nptr, edge_j, k-1, p-1, curr_dof)
	    curr_dof += 1
	  end
	end
      end  # en loop over nodes on edge
    end  # end loop over edges

    # now loop over faces
    for j=1:numFace
      println("face number: ", j)
      face_j = face_i[j]
      for k=1:num_nodes_entity[3]  # loop over nodes on faces
	for p=1:mesh.numDofPerNode  # loop over dofs
	  dofnum_p = getNumberJ(mesh.dofnums_Nptr, region_j, k-1, p-1)
	  if dofnum_p > numDof
	    numberJ(mesh.dofnums_Nptr, face_j, k-1, p-1, curr_dof)
	    curr_dof += 1
	  end
	end
      end  # end loop over face nodes
    end  # end loop over faces

    # now loop over region
    for j=1:numRegions
      println("region number: ", j)
      region_j = regions_i[j]
      for k=1:num_nodes_entity[4]
	for p=1:mesh.numDofPerNode
	  dofnum_p = getNumberJ(mesh.dofnums_Nptr, region_j, k-1, p-1)
	  if dofnum_p > numDof
	    numberJ(mesh.dofnums_Nptr, region_j, k-1, p-1, curr_dof)
	    curr_dof += 1
	  end
	end
      end  # end loop over nodes on the region
    end  # end loop over region

  end  # end loop over elements




  println("finished performing final dof numbering")

  println("number of dofs = ", curr_dof - 1)
  if (curr_dof -1) != numDof 
    println("Warning: number of dofs assigned is not equal to teh expected number")
    println("number of dofs assigned = ", curr_dof-1, " ,expected number = ", numDof)
  else
    println("Dof numbering is sane")
  end



  resetAllIts2()
  return nothing

end

function getDofNumbers(mesh::PumiMesh3)
# populate array of dof numbers, in same shape as solution array u (or q)

println("in getDofNumbers")
println("numNodesPerElement = ", mesh.numNodesPerElement)

mesh.dofs = Array(Int, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)

for i=1:mesh.numEl
  dofnums = getGlobalNodeNumbers(mesh, i)

  for j=1:mesh.numNodesPerElement
    for k=1:mesh.numDofPerNode  # loop over dofs on the node
      mesh.dofs[k, j, i] = dofnums[k,j]
    end
  end
end

#println("mesh.dof = ", mesh.dofs)

return nothing

end
function getGlobalNodeNumbers(mesh::PumiMesh3, elnum::Integer)
# gets global node numbers of all dof on all nodes of the element
# output formap is array [numdofpernode, nnodes]  (each column contains dof numbers for a node)
# 2D only
el_i = mesh.elements[elnum]
type_i = getType(mesh.m_ptr, el_i)


# use apf::getElementNumbers to ensure everything is in the proper orientation
# calculate total number of nodes
nnodes = 0
for i=1:4
  nnodes += mesh.numNodesPerType[i]*mesh.numEntitiesPerElement[i]
end

nnodes = 3 + 3*mesh.numNodesPerType[2] + mesh.numNodesPerType[3]

numdof = nnodes*mesh.numDofPerNode
node_entities = getNodeEntities(mesh.m_ptr, mesh.mshape_ptr, el_i)
dofnums = zeros(Int,mesh.numDofPerNode, nnodes)  # to be populated with global dof numbers
#println("nnodes = ", nnodes)
#println("dofpernode = ", mesh.numDofPerNode)



#num_entities = [3, 3, 1] # number of vertices, edges, faces
num_entities = mesh.num_entities
col = 1 # current node
for i=1:3  # loop over verts, edges, faces
  for j=1:num_entities[i]  # loop over all entities of this type
    for k=1:mesh.numNodesPerType[i]  # loop over nodes on this entity
      entity = node_entities[col]  # current entity
      for p=1:mesh.numDofPerNode  # loop over all dofs
	dofnums[p, col] = getNumberJ(mesh.dofnums_Nptr, entity, k-1, p-1)
      end
      col += 1
    end  # end loop over nodes on curren entity
  end  # end loop over entities of current type
end  # end loop over entity types



#=

for i=1:nnodes
  for j=1:mesh.numDofPerNode
    dofnums[j, i] = getNumberJ(mesh.dofnums_Nptr, node_entities[i], 0, j-1)
  end
end
=#

return dofnums


end


