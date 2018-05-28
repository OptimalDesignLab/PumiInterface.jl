# a type to hold data related to coloring the mesh
#TODO: reverse dimensions of revadj
mutable struct ColoringData
  adj_dict::Dict{Int, Array{Int, 1}}  # map local element to non-locals
  revadj::Array{Int32, 2}  # map non-local elements to local elements

  nonlocal_colors::Array{Int32, 1}
  # put other fields here

  function ColoringData(dict, revadj)
    obj = new()
    obj.adj_dict = dict
    obj.revadj = revadj
    return obj
  end
end

#------------------------------------------------------------------------------
# Distance-0 coloring
#------------------------------------------------------------------------------

# perform a distance-0 coloring of the mesh (ie. all elements same color)
function colorMesh0(mesh::PumiMesh)

  # make every element color 1
  for i=1:mesh.numEl
    numberJ(mesh.coloring_Nptr, mesh.elements[i], 0, 0, 1)
  end

  return 1  # return number of colors
end


#------------------------------------------------------------------------------
# Distance- coloring
#------------------------------------------------------------------------------

# perform distance-1 coloring of mesh 
# not sure if this works correctly
#=
function colorMesh1(mesh::PumiMesh, masks::Array{BitArray{1}})
# each element must have a different color than its neighbors with which it 
# shares and edge

for i=1:mesh.numEl
  numberJ(mesh.coloring_Nptr, mesh.elements[i], 0, 0, 0)
end

adj_size = 6  # guess number of neighboring faces
numc = 4  # guess number of colors
adj = Array{Ptr{Void}}( adj_size)  # hold neighboring faces
adj_color =zeros(Int32, adj_size)  # colors of neighboring faces
cnt_colors = zeros(Int32, numc)  # count how many of each color

for i=1:mesh.numEl
  el_i = mesh.elements[i]
  # get faces that share a vert with el_i
  # this is actually more restrictive than we want because
  # edge stabilization only causes interaction between 
  # elements that share an edge
  # However, this is much easier to write an algorithm for

  # This is more restrictive than we want, but it ensure
  # that we void assigning the same color to two elements
  # bordering the same element
  # ie. we must avoid vertex neighbor non uniqueness,
  # which is not the same thing as having uniqueness because
  # there are unassigned values remaining
  # thus, this is a logical system where a double negative
  # is not a positive
 
  num_adj = countBridgeAdjacent(mesh.m_ptr, el_i, 0, 2)

  if num_adj > adj_size
    println("resizing adj")
    println("element number = ", i)
    resize!(adj, num_adj)
    resize!(adj_color, num_adj)
    adj_size = num_adj
  end


  # need to verify this works in parallel (proper ghosting)
  getBridgeAdjacent(adj)

  for j=1:num_adj
    adj_color[j] = getNumberJ(mesh.coloring_Nptr, adj[j], 0, 0)
  end

  min_color = getMinColor2(adj_color, numc)

  if min_color > numc
    resize!(cnt_colors, min_color)
    cnt_colors[min_color] = 0  # initialize new value to zero
    numc = min_color
  end

  numberJ(mesh.coloring_Nptr, el_i, 0, 0, min_color)

  cnt_colors[min_color] += 1  # update counts
#  masks[min_color][i] = true  # update mask

  fill!(adj_color, 0)

end
println("maximum number of neighbor faces = ", adj_size)
println("number of colors = ", numc)
println("number of each color = ", cnt_colors)
mesh.color_cnt = cnt_colors

return nothing

end
=#




#------------------------------------------------------------------------------
# Distance-2 coloring
#------------------------------------------------------------------------------

# perform distance-1 coloring of mesh
# this can be generalized to 3D via mesh.dim field
function colorMesh2(mesh::PumiMeshDG, colordata::ColoringData)
# each element must have a different color than its neighbors with which it 
# shares and edge

println(mesh.f, "----- Entered colorMesh2 -----")

for i=1:mesh.numEl
  numberJ(mesh.coloring_Nptr, mesh.elements[i], 0, 0, 0)
end

nfaces = mesh.numFacesPerElement
numc = nfaces + 1  # guess number of colors

adj = Array{Ptr{Void}}(nfaces)  # distance-1 edge neighbors
adj2 = Array{Array{Ptr{Void}, 1}}(nfaces)  # distance-2 edge neighbors + distance-1 

for i=1:nfaces
  adj2[i] = Array{Ptr{Void}}(nfaces+1)
end

el_faces = Array{Ptr{Void}}(mesh.numFacesPerElement)
part_nums = Array{Cint}(1)
matched_entities = Array{Ptr{Void}}(1)
adj_entities = Array{Ptr{Void}}(1)
matchdata = LocalNeighborMatches(el_faces, part_nums, matched_entities, adj_entities)

colors = zeros(Int32, mesh.numFacesPerElement*(mesh.numFacesPerElement+1))

cnt_colors = zeros(Int32, numc)  # count how many of each color
adj_dict = colordata.adj_dict

for i=1:mesh.numEl
#  println("Processing element ", i)
  if haskey(adj_dict, i)  # skip any element on a parallel boundary
    continue
  end

  el_i = mesh.elements[i]

  # get faces that share an edge with el_i, and all elements
  # those elements share an edge with
  # this is actually more restrictive than we want because
  # a unique distance-1 coloring, not distance-2
  # However, this is much easier to write an algorithm for

  # This is more restrictive than we want, but it ensure
  # that we avoid assigning the same color to two elements
  # bordering the same element
  # ie. we must avoid distance-2 neighbor non uniqueness,
  # which is not the same thing as having uniqueness because
  # there are unassigned values remaining
  # thus, this is a logical system where a double negative
  # is not a positive
 
  getDistance2Colors(mesh, i, adj, adj2, colors, matchdata)

  min_color = getMinColor2(colors, numc)

  if min_color > numc
    println(mesh.f, "adding color ", min_color)
    resize!(cnt_colors, min_color)
    cnt_colors[min_color] = 0  # initialize new value to zero
    numc = min_color
  end

  numberJ(mesh.coloring_Nptr, el_i, 0, 0, min_color)

  cnt_colors[min_color] += 1  # update counts

  fill!(colors, 0)

end

# now color the parallel boundary elements and their neighbors
numc = colorMeshBoundary2(mesh, colordata, numc, cnt_colors)

println(mesh.f, "number of colors = ", numc)
println(mesh.f, "number of each color = ", cnt_colors)
mesh.color_cnt = cnt_colors

return numc

end

# perform distance-2 coloring of mesh 
function colorMesh2(mesh::PumiMesh2)
# each element must have a different color than its neighbors with which it 
# shares and edge

# figure out the number of colors
# this is a lot of random memory access

for i=1:mesh.numEl
  numberJ(mesh.coloring_Nptr, mesh.elements[i], 0, 0, 0)
end

nfaces = mesh.numFacesPerElement

#adj_size = 6  # guess number of neighboring faces
numc = nfaces + 1  # guess number of colors
#adj = Array{Ptr{Void}}(adj_size)  # hold neighboring faces
#adj_color =zeros(Int32, adj_size)  # colors of neighboring faces

adj = Array{Ptr{Void}}(nfaces)  # distance-1 edge neighbors
adj2 = Array{Array{Ptr{Void}, 1}}(nfaces)  # distance-2 edge neighbors + distance-1 

for i=1:3
  adj2[i] = Array{Ptr{Void}}(nfaces+1)
end

el_faces = Array{Ptr{Void}}(mesh.numFacesPerElement)
part_nums = Array{Cint}(1)
matched_entities = Array{Ptr{Void}}(1)
adj_entities = Array{Ptr{Void}}(1)
matchdata = LocalNeighborMatches(el_faces, part_nums, matched_entities, adj_entities)


# space for distance-2 neighbors + self
colors = zeros(Int32, nfaces*(nfaces+1))

cnt_colors = zeros(Int32, numc)  # count how many of each color
for i=1:mesh.numEl
#  println("Processing element ", i)
  el_i = mesh.elements[i]
  # get faces that share an edge with el_i, and all elements
  # those elements share an edge with
  # this is actually more restrictive than we want because
  # a unique distance-1 coloring, not distance-2
  # However, this is much easier to write an algorithm for

  # This is more restrictive than we want, but it ensure
  # that we avoid assigning the same color to two elements
  # bordering the same element
  # ie. we must avoid distance-2 neighbor non uniqueness,
  # which is not the same thing as having uniqueness because
  # there are unassigned values remaining
  # thus, this is a logical system where a double negative
  # is not a positive
 
  getDistance2Colors(mesh, i, adj, adj2, colors, matchdata)

  min_color = getMinColor2(colors, numc)

  if min_color > numc
    resize!(cnt_colors, min_color)
    cnt_colors[min_color] = 0  # initialize new value to zero
    numc = min_color
  end

  numberJ(mesh.coloring_Nptr, el_i, 0, 0, min_color)

  cnt_colors[min_color] += 1  # update counts

  fill!(colors, 0)

end
mesh.color_cnt = cnt_colors

return numc

end

#TODO: investigate whether duplicating any non-local elements that share a face
#      with more than one local elements produces better coloring
function colorMeshBoundary2(mesh::PumiMeshDG, colordata::ColoringData, numc, cnt_colors)

  println(mesh.f, "----- Entered colorMeshBoundary2 -----")
  nfaces = mesh.numFacesPerElement
  colordata.nonlocal_colors = zeros(Int32, mesh.numSharedEl)
  adj = Array{Ptr{Void}}(nfaces)  # distance-1 edge neighbors
  adj2 = Array{Array{Ptr{Void}, 1}}(nfaces)  # distance-2 edge neighbors + distance-1
  self = Array{Ptr{Void}}(1)  # the current element pointer
  for i=1:nfaces
    adj2[i] = Array{Ptr{Void}}(nfaces+1)
  end

  el_faces = Array{Ptr{Void}}(mesh.numFacesPerElement)
  part_nums = Array{Cint}(1)
  matched_entities = Array{Ptr{Void}}(1)
  adj_entities = Array{Ptr{Void}}(1)
  matchdata = LocalNeighborMatches(el_faces, part_nums, matched_entities, adj_entities)


  const local_start = 1
  const nonlocal_start = nfaces*(nfaces+1) + local_start
  const nonlocal_d2_start = (nfaces-1)*nfaces + nonlocal_start
  const final_start = (nfaces-2)*(nfaces-1) + nonlocal_d2_start

  colors = zeros(Int32, final_start-1)  

  # these calculations are for a 2D mesh (triangles), but are similar for
  # 3d
  # the first 3 sets of four elements are for the colors of the distance-2 
  # neighbors + their common distance-2 neighbor
  # the last 2*3 entries are for the non-local neighbors of the distance-1 
  # neighbors
  # the last 2 entries are for the local distance-2 neighbors connected by a
  # non-local element
  # (for the case where a single element has 2 non local neighbors)
  local_colors = sview(colors, local_start:(nonlocal_start-1))
  nonlocal_d1neighborcolors = sview(colors, nonlocal_start:(nonlocal_d2_start-1))
  nonlocal_neighborcolors = sview(colors, nonlocal_d2_start:(final_start-1))


  for i in keys(colordata.adj_dict)
    el_i = mesh.elements[i]
    self[1] = el_i

    num_adj = getDistance2Colors(mesh, i, adj, adj2, local_colors, matchdata)
    getNonLocalColors(mesh, sview(adj, 1:num_adj), colordata, nonlocal_d1neighborcolors)
    getNonlocalDistance2Colors(mesh, i, colordata, nonlocal_neighborcolors)
    # the non-local elements have not been colored yet, so no need to get their colors
#    getNonLocalColors(mesh, self, colordata, nonlocal_neighborcolors)

    min_color = getMinColor2(colors, numc)
    if min_color > numc
      resize!(cnt_colors, min_color)
      cnt_colors[min_color] = 0  # initialize new value to zero
      numc = min_color
    end
    # record the color
    numberJ(mesh.coloring_Nptr, el_i, 0, 0, min_color)
    cnt_colors[min_color] += 1  # update counts

    fill!(colors, 0)
  end

  # now do non-local elements
  # it is preferable to do these last because they have strictly fewer 
  # adjacent elements than other elemenets, so the odds are better they can 
  # get a lower color by doing coloring all the boundary elements first

  # store local d1 adjacencies in adj (there are nfaces-1 maximum)
  # store colors of their local d1 neighbors + self color in next 
  # (nfaces-1)(nfaces) entries of colors
  # store their nonlocal d1 neighbors in next (nfaces-1)(nfaces-1) entries of
  # colors
  revadj = colordata.revadj
  colors = zeros(Int32, (nfaces-1)*nfaces + (nfaces-1)*(nfaces))

  if mesh.dim == 2
    @assert length(colors) == nfaces*(nfaces+1)
  end

  d1_neighbors = Array{Int32}(nfaces-1)
  local_d2_neighbors = Array{ContiguousView{Int32,1, Array{Int32, 1}}}(nfaces-1)
  nonlocal_d2_neighbors = Array{ContiguousView{Int32,1, Array{Int32, 1}}}(nfaces-1)
  pos = 1
  for i=1:nfaces-1
    # local d2 neighbor + d1 neighbor
    section_start = nfaces + pos
    local_d2_neighbors[i] = view(colors, pos:(section_start-1))

    # nonloca d2 neighbors
    pos = section_start
    section_start = nfaces - 1 + pos
    nonlocal_d2_neighbors[i] = view(colors, pos:(section_start-1))
    pos = section_start
  end

  if mesh.dim == 2
    @assert pos-1 == 10
  end

  d1_ptr = Array{Ptr{Void}}(1)
  for i=1:mesh.numSharedEl

    for j=1:(nfaces-1)
      neighbor = revadj[i, j]
      if neighbor != 0
        d2_local = local_d2_neighbors[j]
        d2_nonlocal = nonlocal_d2_neighbors[j]

        getDistance1Colors(mesh, neighbor, adj, d2_local)
        d1_ptr[1] = mesh.elements[neighbor]
        d2_local[end] = getNumberJ(mesh.coloring_Nptr, d1_ptr[1], 0, 0)
        getNonLocalColors(mesh, d1_ptr, colordata, d2_nonlocal)
      end
    end

    min_color = getMinColor2(colors, numc)

    if min_color > numc
      println(mesh.f, "adding color ", min_color)
      resize!(cnt_colors, min_color)
      cnt_colors[min_color] = 0  # initialize new value to zero
      numc = min_color
    end
    colordata.nonlocal_colors[i] = min_color
    cnt_colors[min_color] += 1
  end

  fill!(colors, 0)

  return numc
end

struct LocalNeighborMatches
  faces::Array{Ptr{Void}, 1}
  part_nums::Array{Cint, 1}
  matched_entities::Array{Ptr{Void}, 1}
  adj_entities::Array{Ptr{Void}, 1}
end

function getNeighborMatches(mesh::PumiMesh, el::Ptr{Void}, adj::AbstractArray{Ptr{Void}}, num_used::Integer, thunk::LocalNeighborMatches)
# check to see if any faces of the element are matched, and get their parent
# elements
# adj is the partially populated array of neighbor elements
# num_used is the number of entries in adj already used
  faces = Array{Ptr{Void}}(mesh.numFacesPerElement)
  n = getDownward(mesh.m_ptr, el, mesh.dim-1, faces)

  part_nums = Array{Cint}(1)
  matched_entities = Array{Ptr{Void}}(1)
  adj_entities = Array{Ptr{Void}}(1)
  for i=1:n
    nmatches = countMatches(mesh.m_ptr, el)
    @assert nmatches <= 1
    if nmatches == 1 && part_nums[1] == mesh.myrank
      other_entity = matched_entities[1]
      nel = countAdjacent(mesh.m_ptr, other_entity, mesh.dim)
      @assert nel == 1
      getAdjacent(adj_entities)

      # put it in the next spot in the array
      num_used += 1
      adj[num_used] = adj_entities[1]
    end
  end

  return num_used
end


function getDistance2Colors(mesh::PumiMesh, elnum::Integer, adj::AbstractArray{Ptr{Void}}, adj2::AbstractArray{Array{Ptr{Void}, 1}}, colors, matchdata::LocalNeighborMatches)
# get the distance-2 neighbors of a given element
# repeats included
# adj : array to be populated with distance-1 edge neighbors
# adj2 : array of arrays to be populated with sum of distance-1 and distance-2 edge
# neighbors
# colors : array to be populated with colors of elements in adj2
# getting the distance-2 colors (instead of the distance-1 colors) is necessary
# to avoid a problem where 2 elements that have a common neighbor element
# will be visited without visiting the common neighbor first, coloring
# the 2 elements the same

  el_i = mesh.elements[elnum]

  num_adj = countBridgeAdjacent(mesh.m_ptr, el_i, mesh.dim-1, mesh.dim)
#  adj = Array{Ptr{Void}}(3)
  getBridgeAdjacent(adj)

  # check for matches
  if num_adj < mesh.numFacesPerElement
    num_adj = getNeighborMatches(mesh, el_i, adj, num_adj, matchdata)
  end

#  adj2 = Array{Array{Ptr{Void}, 1}}(num_adj)  # hold distance 2 neighbors

  adj_cnt = zeros(Int, num_adj)
  for j=1:num_adj
    num_adj_j = countBridgeAdjacent(mesh.m_ptr, adj[j], mesh.dim-1, mesh.dim)
 #   adj2[j] = Array{Ptr{Void}}(num_adj_j + 1)
    getBridgeAdjacent(adj2[j])
    if num_adj_j < mesh.numFacesPerElement
      num_adj_j = getNeighborMatches(mesh, el_i, adj2[j], num_adj_j, matchdata)
    end
    adj2[j][num_adj_j + 1] = adj[j]  # include the distance-1 neighbor
    adj_cnt[j] = num_adj_j + 1
  end

  num_adj_total = sum(adj_cnt)  # total number of neighbors, including repeats

#  colors = zeros(Int32, num_adj_total)

  # get the colors of all the elements
  # also flatten adj2 array of arrays into a single array
  index = 1
  for i=1:num_adj
    for j=1:adj_cnt[i]
      colors[index] = getNumberJ(mesh.coloring_Nptr, adj2[i][j], 0, 0)
      index += 1
    end
  end


  return num_adj
end



function getNonlocalDistance2Colors(mesh::PumiMeshDG, elnum::Integer, colordata::ColoringData, colors)
# get the distance 2 neighbors of a local element that are connected via a non-local element, not including the current element

  el = mesh.elements[elnum]
  nfaces = mesh.numFacesPerElement
  pos = 1
  if haskey(colordata.adj_dict, elnum)
    vals = colordata.adj_dict[elnum]
    for i=1:length(vals)
      if (vals[i] != 0)
        for j=1:(nfaces-1)
          neighbor = colordata.revadj[vals[i] - mesh.numEl, j]
          if neighbor != 0 && neighbor != elnum
            neighbor_ptr = mesh.elements[neighbor]
            colors[pos] = getNumberJ(mesh.coloring_Nptr, neighbor_ptr, 0, 0)
            pos += 1
          end
        end
      end
    end
  end

  return nothing
end
            

function getDistance1Colors(mesh::PumiMeshDG, elnum::Integer, adj, colors)
# get the distance1 colors of a specified element

  el = mesh.elements[elnum]
  num_adj = countBridgeAdjacent(mesh.m_ptr, el, mesh.dim-1, mesh.dim)
  getBridgeAdjacent(adj)

  for i=1:num_adj
    colors[i] = getNumberJ(mesh.coloring_Nptr, adj[i], 0, 0)
  end

  return num_adj
end

function getNonLocalColors(mesh, adj::AbstractArray{Ptr{Void}}, colordata::ColoringData, colors::AbstractArray)
# gets the colors of the non-local neighbors of the specified elements
# adj holds the pointers to the elements

  adj_dict = colordata.adj_dict
  pos = 1
  for i=1:length(adj)
    elnum = getNumberJ(mesh.el_Nptr, adj[i], 0, 0) + 1
    if haskey(adj_dict, elnum)
      vals = adj_dict[elnum]
      for j=1:length(vals)
        if vals[j] != 0
          val_i = vals[j]
          colors[pos] = colordata.nonlocal_colors[val_i - mesh.numEl]
          pos += 1
        end
      end
    end
  end

  return nothing
end


function getMinColor(adj::AbstractArray{T}) where T
# adj contains colors of adjacent elements
  min_color = 1
  sort!(adj)  # adj must be in increasing order for this to work
  for i=1:length(adj)
    if min_color < adj[i]
      continue
    elseif min_color == adj[i]
	min_color = adj[i] + 1
    end  # if min_color > adj[i] do nothing
  end  # end for loop

  return min_color
end

#TODO: see if output can be removed
function getMinColor2(adj::AbstractArray{T}, numc::Integer) where T
# ensure uniqueness of neighboring colors
# adj is array of colors of adjacent faces
# numc is the current number of colors

  mask = zeros(Bool, numc)
  sort!(adj)
  min_color = 0
  for i=1:length(adj) # identify already used colors
    if adj[i] != 0
      mask[adj[i]] = true
    end
  end

  mask_sum = sum(mask)

  if mask_sum == numc  # all existing colors used, so add another
    min_color = numc + 1
  elseif mask_sum == (numc - 1)  # there is exactly 1 color remaining
    # find out which color is missing and use it
    for i=1:numc
      if !mask[i]  # if mask is false
        min_color = i
      end
    end
  else  # some colors are missing
    min_color = getMinColor(adj)  # get the minimum
  end

  @assert min_color != 0

  return min_color

end

#------------------------------------------------------------------------------
# functions that get the bookkeeping information after coloring is done
#------------------------------------------------------------------------------

function getColors0(mesh, masks::AbstractArray{BitArray{1}, 1})
# populate the masks that indicate which elements are perturbed by which
# colors
# for a distance-0 coloring, every element is perturbed and there is only
# 1 color

  masks[1] = trues(mesh.numEl)
  return nothing
end


function getColors1(mesh, masks::AbstractArray{BitArray{1}, 1}, neighbor_colors::AbstractArray{T, 2}, neighbor_nums::AbstractArray{T2, 2}; verify=true ) where {T, T2}
# verify edge neighbor faces have different colors (ie. distance-1 coloring
# of graph where elements are the vertices connected by edges if 
# there is a mesh edge connecting them
# also construct the BitArray masks that describe which elements to perturb
# for which colors, and the list of neighbor + self colors
# masks is an array of bitarrays, number of arrays = number of colors
# neighbor_colors is 4 by numEl array of integers holding the colors
# of the nieghbors + own color
# neighbor_nums is 4 by numEl array of integers holding the element numbers
# of the neighbors + self

nfaces = mesh.numFacesPerElement
adj = Array{Ptr{Void}}(nfaces + 1) # pointers to element i + 3 neighbors
adj_color = zeros(Int32, nfaces + 1)  # element colors

el_faces = Array{Ptr{Void}}(mesh.numFacesPerElement)
part_nums = Array{Cint}(1)
matched_entities = Array{Ptr{Void}}(1)
adj_entities = Array{Ptr{Void}}(1)
matchdata = LocalNeighborMatches(el_faces, part_nums, matched_entities, adj_entities)



# initialize masks to 0 (false)
for i=1:length(masks)
  masks[i] = falses(mesh.numEl)
end

#if verify
#  println("verifying distance-1 coloring")
#end


cnt = 0
for i=1:mesh.numEl
  el_i = mesh.elements[i]

  # check edge neighbors only
  num_adj = countBridgeAdjacent(mesh.m_ptr, el_i, mesh.dim-1, mesh.dim)
  @assert num_adj <= nfaces
  getBridgeAdjacent(adj)

  if num_adj < mesh.numFacesPerElement
    num_adj = getNeighborMatches(mesh, el_i, adj, num_adj, matchdata)
  end

  adj[num_adj + 1] = el_i  # insert the current element
  # get color, element numbers
  for j=1:(num_adj + 1)
    adj_color[j] = getNumberJ(mesh.coloring_Nptr, adj[j], 0, 0)
    neighbor_colors[j, i] = adj_color[j]
    neighbor_nums[j, i] = getNumberJ(mesh.el_Nptr, adj[j], 0, 0) + 1
  end

  color_i = adj_color[num_adj + 1]  # color of current element
  masks[color_i][i] = true  # indicate element i gets perturbed by color_i

  if verify
    sort!(adj_color)

    # remove leading zeros
    nnz = countnz(adj_color)
    nz_arr = zeros(eltype(adj_color), nnz)
    start_idx = length(adj_color) - nnz + 1
    for k=1:nnz
      nz_arr[k] = adj_color[start_idx]
      start_idx += 1
    end

    #TODO: do this by verifying adjacent values are not equal, rather than
    #      calling unique which allocates memory
    if nz_arr != unique(nz_arr)
      println("element ", i, " has non unique colors")
      println("adj_color = ", nz_arr)
      cnt += 1
    end
  end   # end if verify
  fill!(adj_color, 0)
end

if verify
#  println("color-1 verification finished")
  if cnt != 0
    println(STDERR, "number of element with non unique coloring = ", cnt)
    throw(ErrorException("non unique element coloring"))
  end
end

return cnt

end

function getColors1(mesh, colordata::ColoringData, masks::AbstractArray{BitArray{1}, 1}, neighbor_colors::AbstractArray{T, 2}, neighbor_nums::AbstractArray{T2, 2}; verify=true ) where {T, T2}
# verify edge neighbor faces have different colors (ie. distance-1 coloring
# of graph where elements are the vertices connected by edges if 
# there is a mesh edge connecting them
# also construct the BitArray masks that describe which elements to perturb
# for which colors, and the list of neighbor + self colors
# masks is an array of bitarrays, number of arrays = number of colors
# neighbor_colors is 4 by numEl array of integers holding the colors
# of the nieghbors + own color
# neighbor_nums is 4 by numEl array of integers holding the element numbers
# of the neighbors + self

println(mesh.f, "entered getColors")

nfaces = mesh.numFacesPerElement
adj = Array{Ptr{Void}}(nfaces + 1)   # pointers to element i + 3 neighbors
adj_color = zeros(Int32, nfaces + 1)  # element colors
adj_elnum = zeros(Int32, nfaces + 1)  # element numbers
adj_dict = colordata.adj_dict
nonlocal_colors = colordata.nonlocal_colors

# initialize masks to 0 (false)
for i=1:length(masks)
  masks[i] = falses(mesh.numEl)
end

if verify
  println(mesh.f, "verifying distance-1 coloring")
end


cnt = 0
for i=1:mesh.numEl
  el_i = mesh.elements[i]
  elnum = getNumberJ(mesh.el_Nptr, el_i, 0, 0) + 1

  # check edge neighbors only
  num_adj = countBridgeAdjacent(mesh.m_ptr, el_i, mesh.dim-1, mesh.dim)
  @assert num_adj <= nfaces
  getBridgeAdjacent(adj)  
  adj[num_adj + 1] = el_i  # insert the current elementa

  # get color, element numbers for non-local elements
  pos = 1  #current index in the adj_color, adj_elnum arrays
  if haskey(adj_dict, elnum)
  nonlocal_els = adj_dict[elnum]
    for j=1:(nfaces-1)
      nonlocal_elnum = nonlocal_els[j]
      if nonlocal_elnum != 0  # if this is a real ghost neighbor
        elcolor = nonlocal_colors[nonlocal_elnum - mesh.numEl]
        adj_color[pos] = elcolor
        adj_elnum[pos] = nonlocal_elnum
        pos += 1
      end
    end
  end

  # get color, element numbers for local elements + self
  for j=1:(num_adj + 1)
    adj_color[pos] = getNumberJ(mesh.coloring_Nptr, adj[j], 0, 0)
    adj_elnum[pos] = getNumberJ(mesh.el_Nptr, adj[j], 0, 0) + 1
    pos += 1
  end

  # copy them into the global arrays
  for j=1:(nfaces + 1)
    neighbor_colors[j, i] = adj_color[j]
    neighbor_nums[j, i] = adj_elnum[j]
  end

  color_i = adj_color[pos - 1]  # color of current element
  masks[color_i][i] = true  # indicate element i gets perturbed by color_i

  if verify
    sort!(adj_color)

    # remove leading zeros
    nnz = countnz(adj_color)
    nz_arr = zeros(eltype(adj_color), nnz)
    start_idx = length(adj_color) - nnz + 1
    for k=1:nnz
      nz_arr[k] = adj_color[start_idx]
      start_idx += 1
    end

    #TODO: do this by comparing adjacent entries, without allocating memory
    if nz_arr != unique(nz_arr)
      println(mesh.f, "element ", i, " has non unique colors")
      println(mesh.f, "adj_color = ", nz_arr)
      cnt += 1
    end
  end   # end if verify
  fill!(adj_color, 0)
  fill!(adj_elnum, 0)
end

if verify
  println(mesh.f, "color-1 verification finished"); flush(mesh.f)
  if cnt != 0
    throw(ErrorException("Non unique distance-1 coloring: $cnt elements have bad colors"))
  end
  println(mesh.f, "number of element with non unique coloring = ", cnt)
  println(mesh.f, "abc")
end

# now create color masks for non-local elements
println(mesh.f, "getting non-local color masks")
nonlocal_masks = Array{Array{BitArray{1}, 1}}(mesh.npeers)
start_elnum = 1
for i=1:mesh.npeers
  nonlocal_masks[i] = Array{BitArray{1}}(mesh.numColors)
  masks_i = nonlocal_masks[i]  # starting element number for this peer
  numel_i = mesh.shared_element_offsets[i+1] - mesh.shared_element_offsets[i]

  # create the masks BitArrays
  for j=1:mesh.numColors
    masks_i[j] = falses(numel_i)
  end

  # set flag to true if the element has the color
  for k=1:numel_i
    global_elnum = start_elnum + k - 1
    color_k = nonlocal_colors[global_elnum]
    masks_i[color_k][k] = true
  end

  start_elnum += numel_i
end


return cnt, nonlocal_masks

end




function getPertNeighbors0(mesh::PumiMesh)
# get the element that is perturbed for each element for each color
# of a distance-0 coloring
# for a distance-0 coloring, each element is only perturbed by 
# other dofs on the element

  pertNeighborEls = zeros(Int32, mesh.numEl, mesh.numColors)

  for i=1:mesh.numEl
    pertNeighborEls[i, 1] = i
  end

  return pertNeighborEls
end


function getPertNeighbors1(mesh::PumiMesh)
# populate the array with the element that is perturbed for each element 
# for each color for a distance-1 coloring
# element number == 0 if no perturbation

#  println("getting neighbor list")
  pertNeighborEls = zeros(Int32, mesh.numEl, mesh.numColors)

for color = 1:mesh.numColors
  num_neigh = size(mesh.neighbor_colors, 1)
#  fill!(arr, 0)
  for i=1:mesh.numEl
    # find out if current element or its neighbors have the current color
    pos = 0
    for j=1:num_neigh
      if color == mesh.neighbor_colors[j, i]
	pos = j
	break
      end
    end

    if pos != 0
      pertNeighborEls[i, color] = mesh.neighbor_nums[pos, i]
    else
       pertNeighborEls[i, color] = 0
     end

   end  # end loop over elements

end  # end loop over colors

  return pertNeighborEls
end


# perterubed neighbors for edge-based residual
# TODO: make this 3d
function getPertEdgeNeighbors(mesh::PumiMesh)

  neighbor_nums = zeros(Int32, mesh.numEl, 3)

  edges = Array{Ptr{Void}}(3)  # hold element edges
  adj = Array{Ptr{Void}}(2)  # hold adjacent elements

  for i=1:mesh.numEl
    el_i = mesh.elements[i]

    getDownward(mesh.m_ptr, el_i, 1, edges)


    #=
    # check edge neighbors only
    num_adj = countBridgeAdjacent(mesh.m_ptr, el_i, 1, 2)
    @assert num_adj <= 3
    getBridgeAdjacent(adj)
    =#

    # get color, element numbers
    for j=1:3 # loop over edges
      # in use, this array is traveres by numEl first, so we have to 
      # populate it by rows here
      num_adj = countAdjacent(mesh.m_ptr, edges[j], 2)
      @assert num_adj <= 2
      if num_adj == 2  # if there is another adjacnet element
	getAdjacent(adj)

	# figure out which adjacent element is the *other* one
	if adj[1] == el_i
	  other_el = adj[2]
	else
	  other_el = adj[1]
	end

        neighbor_nums[i, j] = getNumberJ(mesh.el_Nptr, other_el, 0, 0) + 1
      end  # end if num_adj == 2
    end  # end loop j=1:3
  end  # end i=1:mesh.numEl

  return neighbor_nums

end  # end getPertEdgeNeighbors


