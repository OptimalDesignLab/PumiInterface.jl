module PdePumiInterface
push!(LOAD_PATH, "/users/creanj/julialib_fork/PUMI.jl")
push!(LOAD_PATH, "/users/creanj/.julia/v0.4/PDESolver/src/common")
using PumiInterface
using SummationByParts
using PDESolverCommon
using ArrayViews

export AbstractMesh,PumiMesh2, reinitPumiMesh2, getElementVertCoords, getShapeFunctionOrder, getGlobalNodeNumber, getGlobalNodeNumbers, getNumEl, getNumEdges, getNumVerts, getNumNodes, getNumDofPerNode, getAdjacentEntityNums, getBoundaryEdgeNums, getBoundaryFaceNums, getBoundaryEdgeLocalNum, getEdgeLocalNum, getBoundaryArray, saveSolutionToMesh, retrieveSolutionFromMesh, retrieveNodeSolution, getAdjacentEntityNums, getNumBoundaryElements, getInterfaceArray

export PumiMesh
#abstract AbstractMesh
abstract PumiMesh{T1} <: AbstractMesh{T1}
include("./PdePumiInterface3.jl")
type PumiMesh2{T1} <: PumiMesh{T1}   # 2d pumi mesh, triangle only
  m_ptr::Ptr{Void}  # pointer to mesh
  mshape_ptr::Ptr{Void} # pointer to mesh's FieldShape
  f_ptr::Ptr{Void} # pointer to apf::field for storing solution during mesh adaptation
  shape_type::Int  #  type of shape functions

  vert_Nptr::Ptr{Void}  # numbering of vertices
  edge_Nptr::Ptr{Void}  # numbering of edges
  el_Nptr::Ptr{Void}  # numbering of elements (faces)

  numVert::Int  # number of vertices in the mesh
  numEdge::Int # number of edges in the mesh
  numEl::Int  # number of elements (faces)
  order::Int # order of shape functions
  numDof::Int # number of degrees of freedom
  numNodes::Int  # number of nodes
  numDofPerNode::Int  # number of dofs per node
  numBoundaryEdges::Int # number of edges on the exterior boundary
  numInterfaces::Int # number of internal interfaces
  numNodesPerElement::Int  # number of nodes per element
  numNodesPerType::Array{Int, 1}  # number of nodes classified on each vertex, edge, face
  numBC::Int  # number of boundary conditions

  # hold pointers to mesh entities
  verts::Array{Ptr{Void},1}  # holds pointers to mesh vertices
  edges::Array{Ptr{Void},1}  # pointers to mesh edges
  elements::Array{Ptr{Void},1}  # pointers to faces

  dofnums_Nptr::Ptr{Void}  # pointer to Numbering of dofs (result of reordering)
#  boundary_nums::Array{Int, 2}  # array of [element number, edgenumber] for each edge on the boundary

  bndry_funcs::Array{BCType, 1}  # array of boundary functors (Abstract type)
  bndry_normals::Array{T1, 3}  # array of normals to each face on the boundary
#  bndry_facenums::Array{Array{Int, 1}, 1}  # hold array of faces corresponding to each boundary condition
  bndry_offsets::Array{Int, 1}  # location in bndryfaces where new type of BC starts
                                # and one past the end of the last BC type
				# array has length numBC + 1

  bndryfaces::Array{Boundary, 1}  # store data on external boundary of mesh
  interfaces::Array{Interface, 1}  # store data on internal edges
  interface_normals::Array{T1, 4}  # array of interior face normals, 
                                   # indexing: [norm_component, Left or right, left face node, left face element]

  coords::Array{T1, 3}  # store coordinates of all nodes
  dxidx::Array{T1, 4}  # store scaled mapping jacobian
  jac::Array{T1,2}  # store mapping jacobian output

  dofs::Array{Int, 3}  # store dof numbers of solution array to speed assembly



 function PumiMesh2(dmg_name::AbstractString, smb_name::AbstractString, order, sbp::SBPOperator, opts; dofpernode=1, shape_type=1)
  # construct pumi mesh by loading the files named
  # dmg_name = name of .dmg (geometry) file to load (use .null to load no file)
  # smb_name = name of .smb (mesh) file to load
  # order = order of shape functions
  # opts: dictionary of options
  # dofpernode = number of dof per node, default = 1
  # shape_type = type of shape functions, 0 = lagrange, 1 = SBP

  mesh = new()
  mesh.numDofPerNode = dofpernode
  mesh.order = order
  mesh.numNodesPerElement = getNumNodes(order)
  mesh.shape_type = shape_type
  num_Entities, mesh.m_ptr, mesh.mshape_ptr = init2(dmg_name, smb_name, order, shape_type=shape_type)
  mesh.f_ptr = createPackedField(mesh.m_ptr, "solution_field", dofpernode)

  mesh.numVert = convert(Int, num_Entities[1])
  mesh.numEdge =convert(Int,  num_Entities[2])
  mesh.numEl = convert(Int, num_Entities[3])

  # get pointers to mesh entity numberings
  mesh.vert_Nptr = getVertNumbering()
  mesh.edge_Nptr = getEdgeNumbering()
  mesh.el_Nptr = getFaceNumbering()

  mesh.verts = Array(Ptr{Void}, mesh.numVert)
  mesh.edges = Array(Ptr{Void}, mesh.numEdge)
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

  for i=1:mesh.numEl
    mesh.elements[i] = getFace()
    incrementFaceIt()
  end

  mesh.numBC = opts["numBC"]

  # create array of all model edges that have a boundary condition
  bndry_edges_all = Array(Int, 0)
  for i=1:mesh.numBC
    key_i = string("BC", i)
    bndry_edges_all = [ bndry_edges_all; opts[key_i]]  # ugly but easy
  end

 mesh.numBoundaryEdges, num_ext_edges =  countBoundaryEdges(mesh, bndry_edges_all)

  # populate mesh.bndry_faces from options dictionary
#  mesh.bndry_faces = Array(Array{Int, 1}, mesh.numBC)
  mesh.bndry_offsets = Array(Int, mesh.numBC + 1)
  mesh.bndry_funcs = Array(BCType, mesh.numBC)
  boundary_nums = Array(Int, mesh.numBoundaryEdges, 2)

  offset = 1
  for i=1:mesh.numBC
    key_i = string("BC", i)
    println("opts[key_i] = ", opts[key_i])
#    println("typeof(opts[key_i]) = ", typeof(opts[key_i]))
    mesh.bndry_offsets[i] = offset
    offset = getMeshEdgesFromModel(mesh, opts[key_i], offset, boundary_nums)  # get the mesh edges on the model edge
    # offset is incremented by getMeshEdgesFromModel
  end


  mesh.bndry_offsets[mesh.numBC + 1] = offset # = num boundary edges

  # get array of all boundary mesh edges in the same order as in mesh.bndry_faces
#  boundary_nums = flattenArray(mesh.bndry_faces[i])
#  boundary_edge_faces = getEdgeFaces(mesh, mesh.bndry_faces)
  # use partially constructed mesh object to populate arrays

  numberDofs(mesh)

  # get boundary information for entire mesh
  mesh.bndryfaces = Array(Boundary, mesh.numBoundaryEdges)
  getBoundaryArray(mesh, boundary_nums)

  # need to count the number of internal interfaces - do this during boundary edge counting
  mesh.numInterfaces = mesh.numEdge - num_ext_edges
  mesh.interfaces = Array(Interface, mesh.numInterfaces)
  getInterfaceArray(mesh)

  getCoordinates(mesh, sbp)  # store coordinates of all nodes into array
  getDofNumbers(mesh)  # store dof numbers

  mesh.dxidx = Array(T1, 2, 2, sbp.numnodes, mesh.numEl)
  mesh.jac = Array(T1, sbp.numnodes, mesh.numEl)
  mappingjacobian!(sbp, mesh.coords, mesh.dxidx, mesh.jac)

  # get face normals
  mesh.bndry_normals = Array(T1, 2, sbp.numfacenodes, mesh.numBoundaryEdges)
  getBoundaryFaceNormals(mesh, sbp, mesh.bndryfaces, mesh.bndry_normals)

  mesh.interface_normals = Array(T1, 2, 2, sbp.numfacenodes, mesh.numInterfaces)
  getInternalFaceNormals(mesh, sbp, mesh.interfaces, mesh.interface_normals)
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

  # write edge and face vertex correspondence to files

  # delete existing files
  if isfile("edge_vertdofs.txt")
    rm("edge_vertdofs.txt")
  end

  if isfile("face_vertdofs.txt")
    rm("face_vertdofs.txt")
  end

  f = open("edge_vertdofs.txt", "a+")
  printEdgeVertNumbers(mesh.edge_Nptr, mesh.dofnums_Nptr, fstream=f)
  close(f)

  f = open("face_vertdofs.txt", "a+")
  printFaceVertNumbers(mesh.el_Nptr, mesh.dofnums_Nptr, fstream=f)
  close(f)

  if isfile("boundary_nums.txt")
    rm("boundary_nums.txt")
  end

  f = open("boundary_nums.txt", "a+")
  println(f, boundary_nums)
  close(f)



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
    nnodes = 3
  elseif order == 2
    nnodes = 7
  elseif order == 3
    nnodes = 12
  elseif order == 4
    nnodes = 18
  else
    println("Warning, unsupported order elements are requested")
  end

  return nnodes
end

function getEdgeFaces(mesh::PumiMesh2, bndry_edges::AbstractArray{Int, 1})
# get the faces corresponding to the boundary edges

bndry_faces = zeros(bndry_edges)

for i=1:bndry_edges
  edgenum_i = bndry_edges[i]
  edge_i = mesh.edges[edgenum_i]
  numFace = countAdjacent(mesh.m_ptr, edge_i, 2)  # should be count upward
  faces = getAdjacent(numFace)
  facenum = getFaceNumber2(faces[1]) + 1

  bndry_faces[i] = facenum
end

return bndry_faces

end



  


function flattenArray{T}(A::AbstractArray{AbstractArray{T}, 1})
# copies array-of-array A into a flattened version B

  n = length(A)
  num_entries = 0
  for i=1:n  # count number of entries total
    num_entries += length(A[i])
  end

  B = zeros(T, num_entries)

  B_index = 1
  for i=1:n
    for j = 1:length(A[i])
      B[B_index] = A[i][j]
      B_index += 1
    end
  end

  return B

end


function getMeshEdgesFromModel{T}(mesh::PumiMesh2, medges::AbstractArray{Int, 1}, offset::Integer, boundary_nums::AbstractArray{T, 2})
# get the global numbers of the mesh edges lying on the model edges in medges
# offset is the index in boundary_nums to start with
# this allows populating the entire array without making temporary copies
# offset should start as 1, not zero

#  bndry_edges = zeros(Int, mesh.numBoundaryEdges)
  index = 0  # relative index in bndry_edges
  for i=1:mesh.numEdge
    edge_i = mesh.edges[i]
    me_i = toModel(mesh.m_ptr, edge_i)
    me_dim = getModelType(mesh.m_ptr, me_i)
    me_tag = getModelTag(mesh.m_ptr, me_i)
    if me_dim == 1  # edge
      onBoundary = findfirst(medges, me_tag)

      if onBoundary != 0  # if mesh edge is on model edge
        
	# get face number
        numFace = countAdjacent(mesh.m_ptr, edge_i, 2)  # should be count upward

        faces = getAdjacent(numFace)
        facenum = getFaceNumber2(faces[1]) + 1
        edgenum = getEdgeNumber2(edge_i)

	boundary_nums[offset + index, 1] = facenum
        boundary_nums[offset + index, 2] = i

	index += 1
      end
    end
  end

  return offset + index

end  # end function


function numberDofs(mesh::PumiMesh2)
# assign dof numbers to entire mesh
# calculates numNodes, numDof
# assumes mesh elements have already been reordered

  println("Entered numberDofs")

  # calculate number of nodes, dofs
  num_nodes_v = 1  # number of nodes on a vertex
  num_nodes_e = countNodesOn(mesh.mshape_ptr, 1) # on edge
  num_nodes_f = countNodesOn(mesh.mshape_ptr, 2) # on face
  println("num_nodes_v = ", num_nodes_v)
  println("num_nodes_e = ", num_nodes_e)
  println("num_nodes_f = ", num_nodes_f)
  numnodes = num_nodes_v*mesh.numVert + num_nodes_e*mesh.numEdge + num_nodes_f*mesh.numEl
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
  iterators_inc = [incrementVertIt, incrementEdgeIt, incrementFaceIt]
  iterators_get = [getVert, getEdge, getFace]
  num_entities = [mesh.numVert, mesh.numEdge, mesh.numEl]
  num_nodes_entity = [num_nodes_v, num_nodes_e, num_nodes_f]

  mesh.numNodesPerType = num_nodes_entity

  println("num_entities = ", num_entities)
  println("num_nodes_entity = ", num_nodes_entity)

#  curr_dof = 1
  curr_dof = numDof+1
  for etype = 1:3 # loop over entity types
#    println("etype = ", etype)
    if (num_nodes_entity[etype] != 0)  # if no nodes on this type of entity, skip
      for entity = 1:num_entities[etype]  # loop over all entities of this type
#	println("  entity number: ", entity)
	entity_ptr = iterators_get[etype]()  # get entity

	for node = 1:num_nodes_entity[etype]
#	  println("    node : ", node)
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
#    println("element number: ", i)
    # get vertices, edges for this element
    (verts_i, numVert) = getDownward(mesh.m_ptr, mesh.elements[i], 0)
#    println("verts_i = ", verts_i)
#    println("mesh.verts = ", mesh.verts)
    (edges_i, numEdge) = getDownward(mesh.m_ptr, mesh.elements[i], 1)
    el_i = mesh.elements[i]
    for j=1:3  # loop over vertices, edges
#      println("  vertex and edge number: ", j)
      vert_j = verts_i[j]
      edge_j = edges_i[j]
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
      
      # loop over nodes on edge
      for k=1:num_nodes_entity[2]  # loop over nodes
	for p=1:mesh.numDofPerNode  # loop over dofs
	  dofnum_p = getNumberJ(mesh.dofnums_Nptr, edge_j, k-1, p-1)
	  if dofnum_p > numDof  # still has initial number
	    # give it new (final) number
	    numberJ(mesh.dofnums_Nptr, edge_j, k-1, p-1, curr_dof)
	    curr_dof += 1
	  end
	end
      end
    end  # end loop over vertices, edges
    # label face nodes
    for k=1:num_nodes_entity[3]  # loop over nodes on face
      for p=1:mesh.numDofPerNode  # loop over dofs
	dofnum_p = getNumberJ(mesh.dofnums_Nptr, el_i, k-1, p-1)
	if dofnum_p > numDof
	  numberJ(mesh.dofnums_Nptr, el_i, k-1, p-1, curr_dof)
	  curr_dof += 1
	end
      end
    end  # end loop over face nodes
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

#=
function numberDofs(mesh::PumiMesh2)
# number the degrees of freedom of the mesh, using the apf::Numbering* stroed
# in the mesh

  # move this into a function
  resetAllIts2()
  println("performing initial numbering of dofs")
  # calculate number of nodes, dofs (works for first and second order)
  numnodes = mesh.order*mesh.numVert 
  numdof = numnodes*mesh.numDofPerNode
  # number dofs
  ctr= 1
  for i=1:mesh.numVert
    for j=1:mesh.numDofPerNode
      numberJ(mesh.dofnums_Nptr, mesh.verts[i], 0, j-1, ctr)
      println("vertex ", i,  " numbered ", ctr)
      ctr += 1
    end
  end

  if mesh.order >= 2
    for i=1:mesh.numEdges
      for j=1:mesh.numDofPerNode
        numberJ(mesh.dofnums_Nptr, mesh.edges[i], 0, j-1, ctr)
        ctr += 1
      end
    end
  end

  mesh.numNodes = numnodes
  mesh.numDof = numdof

return nothing

end  # end function
=#


function countBoundaryEdges(mesh::PumiMesh2, bndry_edges_all)
  # count boundary edges by checking if their model edge has a BC
  # count number of external edges by checking the number of upward adjacencies
  # store array of [element number, global edge number]
  resetEdgeIt()
  bnd_edges_cnt = 0
  external_edges_cnt = 0
  bnd_edges = Array(Int, mesh.numEdge, 2)
  for i=1:mesh.numEdge
    edge_i = getEdge()

    # get  model edge info
    me_i = toModel(mesh.m_ptr, edge_i)
    me_dim = getModelType(mesh.m_ptr, me_i)
    me_tag = getModelTag(mesh.m_ptr, me_i)

    # get mesh face info
    numFace = countAdjacent(mesh.m_ptr, edge_i, 2)  # should be count upward
    if numFace == 1  # external edges
      external_edges_cnt += 1
    end

    if me_dim == 1  # if not classified on model edge
      index = findfirst(bndry_edges_all, me_tag)



      if index != 0  # if model edge has a BC on i

	faces = getAdjacent(numFace)
	facenum = getFaceNumber2(faces[1]) + 1



	bnd_edges_cnt += 1
	bnd_edges[bnd_edges_cnt, 1] = facenum
	bnd_edges[bnd_edges_cnt, 2] = i
      end
    end  # end if me_dim == 1
    incrementEdgeIt()

  end  # end for loop


#  mesh.boundary_nums = bnd_edges[1:bnd_edges_cnt, :] # copy, bad but unavoidable
#  mesh.numBoundaryEdges = bnd_edges_cnt

return bnd_edges_cnt, external_edges_cnt

end  # end function




function getDofNumbers(mesh::PumiMesh2)
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
# for reinitilizeing after mesh adaptation
function reinitPumiMesh2(mesh::PumiMesh2)
  # construct pumi mesh by loading the files named
  # dmg_name = name of .dmg (geometry) file to load (use .null to load no file)
  # smb_name = name of .smb (mesh) file to load
  # order = order of shape functions
  # dofpernode = number of dof per node, default = 1

  println("Reinitilizng PumiMesh2")

  # create random filenames because they are not used
  smb_name = "a"
  dmg_name = "b"
  order = mesh.order
  dofpernode = mesh.numDofPerNode
  tmp, num_Entities, m_ptr, mshape_ptr = init2(dmg_name, smb_name, order, load_mesh=false, shape_type=mesh.shape_type) # do not load new mesh
  f_ptr = mesh.f_ptr  # use existing solution field

  numVert = convert(Int, num_Entities[1])
  numEdge =convert(Int,  num_Entities[2])
  numEl = convert(Int, num_Entities[3])

  mesh.vert_Nptr = getVertNumbering()
  mesh.edge_Nptr = getEdgeNumbering()
  mesh.el_Nptr = getFaceNumbering()




  verts = Array(Ptr{Void}, numVert)
  edges = Array(Ptr{Void}, numEdge)
  elements = Array(Ptr{Void}, numEl)
  dofnums_Nptr = mesh.dofnums_Nptr  # use existing dof pointers

  # get pointers to all MeshEntities
  # also initilize the field to zero
  resetAllIts2()
#  comps = zeros(dofpernode)
  comps = [1.0, 2, 3, 4]
  for i=1:numVert
    verts[i] = getVert()
    incrementVertIt()
  end

  for i=1:numEdge
    edges[i] = getEdge()
    incrementEdgeIt()
  end

  for i=1:numEl
    elements[i] = getFace()
    incrementFaceIt()
  end

  resetAllIts2()
  println("performing initial numbering of dofs")
  # calculate number of nodes, dofs (works for first and second order)
  numnodes = order*numVert 
  numdof = numnodes*dofpernode
  # number dofs
  ctr= 1
  for i=1:numVert
    for j=1:dofpernode
      numberJ(dofnums_Nptr, verts[i], 0, j-1, ctr)
#      println("vertex ", i,  " numbered ", ctr)
      ctr += 1
    end
  end

  if order >= 2
    for i=1:numEdges
      for j=1:dofpernode
        numberJ(dofnums_Nptr, edges[i], 0, j-1, ctr)
        ctr += 1
      end
    end
  end

 

  # count boundary edges
  bnd_edges_cnt = 0
  bnd_edges = Array(Int, numEdge, 2)
  for i=1:numEdge
    edge_i = getEdge()
    numFace = countAdjacent(m_ptr, edge_i, 2)  # should be count upward

    if numFace == 1  # if an exterior edge
      faces = getAdjacent(numFace)
      facenum = getFaceNumber2(faces[1]) + 1

      bnd_edges_cnt += 1
      bnd_edges[bnd_edges_cnt, 1] = facenum
      bnd_edges[bnd_edges_cnt, 2] = i
    end
    incrementEdgeIt()
  end

  bnd_edges_small = bnd_edges[1:bnd_edges_cnt, :]

  # replace exising fields with new values
  mesh.numVert = numVert
  mesh.numEdge = numEdge
  mesh.numEl = numEl
  mesh.numDof = numdof
  mesh.numNodes= numnodes
  mesh.numBoundaryEdges = bnd_edges_cnt
  mesh.verts = verts  # does this need to be a deep copy?
  mesh.edges = edges
  mesh.elements = elements
#  mesh.boundary_nums = bnd_edges_small

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
  println("numVert = ", numVert)
  println("numEdge = ", numEdge)
  println("numEl = ", numEl)
  println("numDof = ", numdof)
  println("numNodes = ", numnodes)

  writeVtkFiles("mesh_complete", m_ptr)
end


function getCoordinates(mesh::PumiMesh2, sbp::SBPOperator)
# populate the coords array of the mesh object

mesh.coords = Array(Float64, 2, sbp.numnodes, mesh.numEl)

println("entered getCoordinates")

coords_i = zeros(3,3)
coords_it = zeros(3,2)
for i=1:mesh.numEl  # loop over elements
  
  el_i = mesh.elements[i]
  (sizex, sizey) = size(coords_i)
  getFaceCoords(el_i, coords_i, sizex, sizey)  # populate coords

#  println("coords_i = ", coords_i)

  coords_it[:,:] = coords_i[1:2, :].'
#  println("coords_it = ", coords_it)
  mesh.coords[:, :, i] = calcnodes(sbp, coords_it)
#  println("mesh.coords[:,:,i] = ", mesh.coords[:,:,i])
end

return nothing

end




function getElementVertCoords(mesh::PumiMesh2, elnum::Integer, coords::AbstractArray{Float64,2})
# get the coordinates of the vertices of an element
# elnum is the number of the element
# coords is the array to be populated with the coordinates
# each column of coords contains the coordinates for a vertex
# coords must be 3x3

  el_j = mesh.elements[elnum]
  (sizex, sizey) = size(coords)
  getFaceCoords(el_j, coords, sizex, sizey)  # populate coords

  return nothing

end # end function

function getElementVertCoords(mesh::PumiMesh2,  elnums::Array{Int,1})
# elnums = vector of element numbbers
# return array of size 3x3xn, where each column (first index) contains the coordinates of a vertex
# n is the number of elements in elnums


# get the number of verts 
n = length(elnums)

coords = zeros(3, 3, n)  # first dimension = x,y, or z, second dimension = vertex number

for j = 1:n  # loop over elements in elnums
#  println("j = ", j)
  elnum_j = elnums[j]
#  println("elnum_j = ", elnum_j)
  el_j = mesh.elements[elnum_j]
  sub_array = sub(coords, :, :, j)
  getFaceCoords(el_j, sub_array, 3, 3)  # populate coords
end

return coords

end

function getShapeFunctionOrder(mesh::PumiMesh2)

return mesh.order
end


function getGlobalNodeNumber(mesh::PumiMesh2, el_num::Integer, local_node_num::Integer)
# el_num = element number
# local_node_num = vector of dof numbers on the specified node within element
# this only works for first and second order
# returns a vector with dof numbers of all dofs on the node


#println("el_num = ", el_num, " local_node_num = ", local_node_num)
el = mesh.elements[el_num]  # get element
node_entities = getNodeEntities(mesh.m_ptr, mesh.mshape_ptr, el)
#println("node_entities = ", node_entities)
node_entity = node_entities[local_node_num]
#println("node_entity = ", node_entity)
dofnums = zeros(Int, mesh.numDofPerNode)
for i=1:mesh.numDofPerNode
  dofnums[i] = getNumberJ(mesh.dofnums_Nptr, node_entities[local_node_num], 0, i-1)

#println("global dof number = ", number_i)
end

return dofnums

end

function getGlobalNodeNumbers(mesh::PumiMesh2, elnum::Integer)
# gets global node numbers of all dof on all nodes of the element
# output formap is array [numdofpernode, nnodes]  (each column contains dof numbers for a node)
# 2D only
el_i = mesh.elements[elnum]
type_i = getType(mesh.m_ptr, el_i)

# calculate total number of nodes
nnodes = 3 + 3*mesh.numNodesPerType[2] + mesh.numNodesPerType[3]

numdof = nnodes*mesh.numDofPerNode
node_entities = getNodeEntities(mesh.m_ptr, mesh.mshape_ptr, el_i)
dofnums = zeros(Int,mesh.numDofPerNode, nnodes)  # to be populated with global dof numbers
#println("nnodes = ", nnodes)
#println("dofpernode = ", mesh.numDofPerNode)



num_entities = [3, 3, 1] # number of vertices, edges, faces

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

function getEdgeNumber(mesh::PumiMesh2, edge_num::Integer)
# get the number of an edge

  i = getNumberJ(mesh.edge_Nptr, mesh.edges[edge_num], 0, 0)

  return i

end


function getNumEl(mesh::PumiMesh2)
# returns the number of elements in the mesh

return mesh.numEl
end

function getNumEdges(mesh::PumiMesh2)
# retrns the number of edges in a mesh

return mesh.numEdge

end

function getNumVerts(mesh::PumiMesh2)
# returns number of vertices in a mesh

return mesh.numVert
end

function getNumNodes(mesh::PumiMesh2)
# returns total number of nodes in the mesh

return mesh.numDof
end

function getNumDofPerNode(mesh::PumiMesh2)
# get the number of dof on a node

return mesh.numDofPerNode

end

function getNumBoundaryEdges(mesh::PumiMesh2)
# return the number of edges on the boundary

  return mesh.numBoundaryEdges
end

function getNumBoundaryElements(mesh::PumiMesh2)
# count the number of elements on the boundary

   return length(unique(mesh.boundary_nums[:,1]))
end

#=
function getBoundaryEdgeNums(mesh::PumiMesh2)
# get vector of edge numbers that are on boundary

edge_nums = mesh.boundary_nums
return edge_nums
end
=#


# this doesn't exist for a 2D mesh
function getBoundaryFaceNums(mesh::PumiMesh2)
# get array of face numbers that are on boundary
# this only works in 2d
# this function doesn't work yet
#face_nums = zeros(Int, numEdges_on_boundary)
face_nums = zeros(Int, 0)
return face_nums
end


function getAdjacentEntityNums(mesh::PumiMesh, entity_index::Integer, input_dimension::Integer, output_dimension::Integer)
# gets the numbers of the adjacent entities (upward or downward) of the given entity
# entity_index specifies the index of the input entity
# input_dimension specifies the dimension of the input entity (0 =vert, 1 = edge ...)
# output_dimension is the dimension of the adjacent  entities to be retrieved
# returns an array containing the indicies of the adjacent entities, and the number of them
# the length of the array might be greater than the number of adjacencies, so use the second returned value for the number of adjacencies


if (input_dimension == output_dimension)
  println("cannot get same level adjacencies")
end


if typeof(mesh) <: PumiMesh2

  array1 = [mesh.vert_Nptr, mesh.edge_Nptr, mesh.el_Nptr]
  #array2 = [mesh.verts; mesh.edges; mesh.elements]  # array of arrays
  array2 = Array(Array{Ptr{Void},1}, 3)
  array2[1] = mesh.verts
  array2[2] = mesh.edges
  array2[3] = mesh.elements
else # 3d mesh
   array1 = [mesh.vert_Nptr, mesh.edge_Nptr, mesh.face_Nptr,  mesh.el_Nptr]
  #array2 = [mesh.verts; mesh.edges; mesh.elements]  # array of arrays
  array2 = Array(Array{Ptr{Void},1}, 4)
  array2[1] = mesh.verts
  array2[2] = mesh.edges
  array2[3] = mesh.faces
  array2[4] = mesh.elements
end

# choose which numbering to use
numbering_ptr = array1[output_dimension + 1]
#println("numbering_ptr = ", numbering_ptr)


# get the the entity
#println("array2 = ", array2)
entity_array = array2[input_dimension + 1]
entity = entity_array[entity_index]
#entity = array2[input_dimension+1][entity_index]
#println("entity = ", entity)

#=
entity_type = getType(mesh.m_ptr, entity)
println("entity_type = ", entity_type)
entity_num = getNumberJ(mesh.edge_Nptr, entity, 0, 0)
println("edge number = ", entity_num)
=#

if input_dimension > output_dimension # downward adjacencies
#  println("getting downward adjacencies")
  adjacent_entities, num_adjacent = getDownward(mesh.m_ptr, entity, output_dimension)

else  # upward adjacencies
#  println("getting upward adjacencies")
  num_adjacent = countAdjacent(mesh.m_ptr, entity, output_dimension)
  adjacent_entities = getAdjacent(num_adjacent)
end

# get their numbers
adjacent_nums = zeros(Int, num_adjacent)
for i=1:num_adjacent # loop over adjacent entities
  
  # get their numbers here
  adjacent_nums[i] = getNumberJ(numbering_ptr, adjacent_entities[i], 0, 0) + 1
end

return adjacent_nums, num_adjacent


end



function getBoundaryEdgeLocalNum(mesh::PumiMesh2, edge_num::Integer)
# gets the local edge number of a specified edge that is on the boundary
# of the mesh
# edge_num is an edge num from the output of getBoundaryEdgeNums() (ie. the global edge number)
# the local edge number is the edge number within the element (1st,s 2nd, 3rd...)
  edge_i = mesh.edges[edge_num]

  # get mesh face associated with edge
  countAdjacent(mesh.m_ptr, edge_i, 2)
  face = getAdjacent(1)[1]  # get the single face (not an array)

  facenum_i = getFaceNumber2(face) + 1  # convert to 1 based indexing

  (down_edges, numedges) = getDownward(mesh.m_ptr, face, 1)
  edgenum_local = 0
  for j = 1:3  # find which local edge is edge_i
    if down_edges[j] == edge_i
      edgenum_local = j
    end
  end

  return edgenum_local

end

function getEdgeLocalNum(mesh::PumiMesh2, edge_num::Integer, element_num::Integer)
# find the local edge number of a specified edge on a specified element

  edge = mesh.edges[edge_num]
  element = mesh.elements[element_num]

  (down_edges, numedges) = getDownward(mesh.m_ptr, element, 1)

  edgenum_local = 0
   for j = 1:3  # find which local edge is edge_i
    if down_edges[j] == edge
      edgenum_local = j
    end
  end

  return edgenum_local

end

function getBoundaryArray(mesh::PumiMesh2, boundary_nums::AbstractArray{Int, 2})
# get an array of type Boundary for SBP
# creating an an array of a user defined type seems like a waste of memory operations
# bnd_array is a vector of type Boundary with length equal to the number of edges on the boundary of the mesh

#  mesh.bndryfaces = Array(Boundary, mesh.numBoundaryEdges)

  for i=1:mesh.numBoundaryEdges
    facenum = boundary_nums[i, 1]
    edgenum_global = boundary_nums[i, 2]
#    println("edgenum_global = ", edgenum_global)
    edgenum_local = getBoundaryEdgeLocalNum(mesh, edgenum_global)
    mesh.bndryfaces[i] = Boundary(facenum, edgenum_local)
  end

  return nothing
end


function getInterfaceArray(mesh::PumiMesh2)
# get array of [elementL, elementR, edgeL, edgeR] for each internal edge,
# where elementL and R are the elements that use the edge, edgeL R are the
# local number of the edge within the element
# interfaces is the array to be populated with the data, it must be 
# number of internal edges by 1

  # only need internal boundaries (not external)
#  num_ext_edges = size(getBoundaryEdgeNums(mesh))[1]  # bad memory efficiency
#  num_int_edges = getNumEdges(mesh) - num_ext_edges

#  new_bndry = Boundary(2, 3)
#   println("new_bndry = ", new_bndry)

#  new_interface = Interface(1, 2, 3, 4)
#   println("new_interface = ", new_interface)

#   println("num_int_edges = ", num_int_edges)

#  interfaces = Array(typeof(new_interface), num_int_edges)
  # get number of nodes affecting an edge
  num_edge_nodes = countAllNodes(mesh.mshape_ptr, 1)

  nodemap = Array(num_edge_nodes:(-1):1)

  pos = 1 # current position in interfaces
  for i=1:getNumEdges(mesh)
#     println("i = ", i)
#     println("pos = ", pos)
    # get number of elements using the edge
    adjacent_nums, num_adjacent = getAdjacentEntityNums(mesh, i, 1, 2)
#     println("num_adjacent = ", num_adjacent)
#     println("adjacent_nums = ", adjacent_nums)
    if num_adjacent > 1  # internal edge
#       println("this is an internal edge")
      element1 = adjacent_nums[1]
      element2 = adjacent_nums[2]

#      coords_1 = x[:, :, element1]
#      coords_2 = x[:, :, element2]

      coords_1 = zeros(3,3)
      coords_2 = zeros(3,3)
      getFaceCoords(mesh.elements[element1], coords_1, 3, 3)
      getFaceCoords(mesh.elements[element2], coords_2, 3, 3)

      # calculate centroid
      centroid1 = sum(coords_1, 2)
      centroid2 = sum(coords_2, 2)

      if abs(centroid1[1] - centroid2[2]) < 1e-10  # if big enough difference
        if centroid1[1] < centroid2[1]
    elementL = element1
    elementR = element2
        else
    elementL = element2
    elementR = element1
        end
      else  # use y coordinate to decide
        if centroid1[2] < centroid2[2]
    elementL = element1
    elementR = element2
        else
    elementL = element2
    elementR = element1
        end
      end

      edgeL = getEdgeLocalNum(mesh, i, elementL)
      edgeR = getEdgeLocalNum(mesh, i, elementR)

      mesh.interfaces[pos] = Interface(elementL, elementR, edgeL, edgeR, nodemap)
#       println("updating pos")
      pos += 1

  #    print("\n")

    end  # end if internal edge

#     print("\n")
  end  # end loop over edges

  return nothing

end  # end function



function getBoundaryFaceNormals{Tmsh}(mesh::PumiMesh2, sbp::SBPOperator, bndry_faces::AbstractArray{Boundary, 1}, face_normals::Array{Tmsh, 3})

  nfaces = length(bndry_faces)

  alpha = zeros(Tmsh, 2,2)
  dxidx = mesh.dxidx
  jac = mesh.jac
  for i=1:nfaces
    element_i = bndry_faces[i].element
    face_i = bndry_faces[i].face
    for j=1:sbp.numfacenodes
      node_index = sbp.facenodes[j, face_i]

      # calculate the mysterious parameter alpha
      for di1=1:2
	for di2=1:2
	  alpha[di1, di2] = (dxidx[di1,1, node_index, element_i].*dxidx[di2,1, node_index, element_i] + dxidx[di1,2, node_index, element_i].*dxidx[di2,2, node_index, element_i])*jac[node_index,element_i]
	end
      end

    # call SBP function
    getdir!(alpha, view(sbp.facenormal, :, face_i), view(face_normals, :, j, i))

    end  # end loop over face nodes
  end  # end loop over faces

  return nothing
end


function getInternalFaceNormals{Tmsh}(mesh::PumiMesh2, sbp::SBPOperator, internal_faces::AbstractArray{Interface, 1}, face_normals::Array{Tmsh, 4})

  nfaces = length(internal_faces)

  alpha = zeros(Tmsh, 2,2)
  dxidx = mesh.dxidx
  jac = mesh.jac
  for i=1:nfaces
    element_iL = internal_faces[i].elementL
    face_iL = internal_faces[i].faceL
    element_iR = internal_faces[i].elementR
    face_iR = internal_faces[i].faceR
    for j=1:sbp.numfacenodes

      # calculate left face normal
      node_index = sbp.facenodes[j, face_iL]

      # calculate the mysterious parameter alpha
      for di1=1:2
	for di2=1:2
	  alpha[di1, di2] = (dxidx[di1,1, node_index, element_iL].*dxidx[di2,1, node_index, element_iL] + dxidx[di1,2, node_index, element_iL].*dxidx[di2,2, node_index, element_iL])*jac[node_index,element_iL]
	end
      end

    # call SBP function
    getdir!(alpha, view(sbp.facenormal, :, face_iL), view(face_normals, :, 1, j, i))
      # calculate right fae normal
      node_index = sbp.facenodes[j, face_iR]

      # calculate the mysterious parameter alpha
      for di1=1:2
	for di2=1:2
	  alpha[di1, di2] = (dxidx[di1,1, node_index, element_iR].*dxidx[di2,1, node_index, element_iR] + dxidx[di1,2, node_index, element_iR].*dxidx[di2,2, node_index, element_iR])*jac[node_index,element_iR]
	end
      end

    # call SBP function
    getdir!(alpha, view(sbp.facenormal, :, face_iR), view(face_normals, :, 2, j, i))




    end  # end loop over face nodes
  end  # end loop over faces

  return nothing
end







#=
function getEdgeInterfaceData(i::Integer)

#     println("i = ", i)
#     println("pos = ", pos)
    # get number of elements using the edge
    adjacent_nums, num_adjacent = getAdjacentEntityNums(mesh, i, 1, 2)
#     println("num_adjacent = ", num_adjacent)
#     println("adjacent_nums = ", adjacent_nums)
    if num_adjacent > 1  # internal edge
#       println("this is an internal edge")
      element1 = adjacent_nums[1]
      element2 = adjacent_nums[2]

      coords_1 = x[:, :, element1]
      coords_2 = x[:, :, element2]

      # calculate centroid
      centroid1 = sum(coords_1, 2)
      centroid2 = sum(coords_2, 2)

      if abs(centroid1[1] - centroid2[2]) < 1e-10  # if big enough difference
        if centroid1[1] < centroid2[1]
    elementL = element1
    elementR = element2
        else
    elementL = element2
    elementR = element1
        end
      else  # use y coordinate to decide
        if centroid1[2] < centroid2[2]
    elementL = element1
    elementR = element2
        else
    elementL = element2
    elementR = element1
        end
      end

      edgeL = getEdgeLocalNum(mesh, i, elementL)
      edgeR = getEdgeLocalNum(mesh, i, elementR)

      return elementL, elementR, edgeL, edgeR
    end

=#



function saveSolutionToMesh(mesh::PumiMesh2, u::AbstractVector)
# saves the solution in the vector u to the mesh (in preparation for mesh adaptation

 # dofnums = zeros(Int, mesh.numDofPerNode)
#  u_vals = zeros(mesh.numDofPerNode)


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
          
          # get solution values
	  for p=1:mesh.numDofPerNode  # loop over all dofs
	    dofnum_p = getNumberJ(mesh.dofnums_Nptr, entity, k-1, p-1)
	    q_vals[p] = u[dofnum_p]
	  end
          # save to mesh
          setComponents(mesh.f_ptr, entity, k-1, q_vals)

	  col += 1
	end  # end loop over nodes on curren entity
      end  # end loop over entities of current type
    end  # end loop over entity types
  end  # end loop over elements


  return nothing
end  # end function saveSolutionToMesh


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
         
	  # get values from mesh
          getComponents(mesh.f_ptr, entity, k-1, q_vals)

          # put solution values into vector u
	  for p=1:mesh.numDofPerNode  # loop over all dofs
	    dofnum_p = getNumberJ(mesh.dofnums_Nptr, entity, k-1, p-1)
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

end  # end of module
