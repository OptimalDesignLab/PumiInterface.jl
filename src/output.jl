# functions for writing output

function writeCounts(mesh::PumiMeshDG2; fname="counts")
# write values needed for memory usage estimate
vals = Array(Int, 8)
vals[1] = mesh.numVert
vals[2] = mesh.numEdge
vals[3] = mesh.numEl
vals[4] = mesh.numBoundaryEdges
vals[5] = mesh.numNodesPerType[1]
vals[6] = mesh.numNodesPerType[2]
vals[7] = mesh.numNodesPerType[3]
vals[8] = mesh.numDofPerNode

#=
# estimate jacobian storage size
acc = 0
for i=1:mesh.numDof
  acc += mesh.sparsity_bnds[2, i] - mesh.sparsity_bnds[1, i]
end

size_nz = 64*acc
size_rowval = 64*acc
size_colptr = 64*mesh.numDof

vals[9] = size_nz + size_rowval + size_colptr
=#
# append mpi_rank, file extension
myrank = mesh.myrank
fname2 = string(fname, "_", myrank, ".txt")
writedlm(fname2, vals)
			  
return nothing
end


function printBoundaryEdgeNums(mesh::PumiMesh)

  n = mesh.numBC

  bndry = 1
  edges = Array(Ptr{Void}, 12)
  for i=1:n
    fname = string("boundary_edge_verts", i, ".dat")
    println("printing ", fname)
    f = open(fname, "a+")

    start_index = mesh.bndry_offsets[i]
    end_index = mesh.bndry_offsets[i+1] - 1
    num_edge = end_index - start_index + 1
    arr = Array(Ptr{Void}, num_edge)


    for i=1:num_edge  # get the mesh edge pointers
      el = mesh.bndryfaces[bndry].element
      local_face = mesh.bndryfaces[bndry].face
      el_ptr = mesh.elements[el]
      getDownward(mesh.m_ptr, el_ptr, 1, edges)
      arr[i] = edges[local_face]
      bndry += 1
    end
    
    printEdgeVertNumbers(arr, mesh.edge_Nptr, mesh.vert_Nptr; fstream=f)
    close(f)
  end

  return nothing
end  # end function


function printdxidx(name::AbstractString, mat)
# print the 4d matrix formatted like dxidx

(p, p, m, n) = size(mat)

f = open(name, "a+")

for i=1:n  # loop over elements
  for j=1:m  # loop over nodes on an element
    println(f, "el ", i, " node ", j, " xi_vec = ", mat[1, :, j, i], " eta_vec = ", mat[2, :, j, i])
  end
end

close(f)

return nothing
end


function printcoords(name::AbstractString, coords)

  f = open(name, "a+")

  (d, numnodes, numel) = size(coords)

  for i=1:numel  # loop over elements
    for j=1:numnodes # loop over nodes
      print(f, i, " ", j)

      for k=1:d
	print(f, " ",  coords[k, j, i])
      end

      print(f, "\n")
    end
  end

  return nothing
end

function writeCounts(mesh::PumiMesh2; fname="counts_0.txt")
# write values needed for memory usage estimate
vals = Array(Int, 9)
vals[1] = mesh.numVert
vals[2] = mesh.numEdge
vals[3] = mesh.numEl
vals[4] = mesh.numBoundaryEdges
vals[5] = mesh.numNodesPerType[1]
vals[6] = mesh.numNodesPerType[2]
vals[7] = mesh.numNodesPerType[3]
vals[8] = mesh.numDofPerNode


# estimate jacobian storage size
acc = 0
for i=1:mesh.numDof
  acc += mesh.sparsity_bnds[2, i] - mesh.sparsity_bnds[1, i]
end

size_nz = 64*acc
size_rowval = 64*acc
size_colptr = 64*mesh.numDof

vals[9] = size_nz + size_rowval + size_colptr

writedlm(fname, vals)
			  
return nothing
end



