# perform a distance-0 coloring of the mesh (ie. all elements same color)
function colorMesh0(mesh::PumiMeshDG2)

  # make every element color 1
  for i=1:mesh.numEl
    numberJ(mesh.coloring_Nptr, mesh.elements[i], 0, 0, 1)
  end

  return 1  # return number of colors
end



# perform distance-1 coloring of mesh 
# not sure if this works correctly
function colorMesh1(mesh::PumiMeshDG2, masks::Array{BitArray{1}})
# each element must have a different color than its neighbors with which it 
# shares and edge

# figure out the number of colors
# this is a lot of random memory access
#=
num_neigh_max = 0
for i=1:mesh.numEl  # loop over elements
  el_i = mesh.elements[i]
  num_neigh_i = countBridgeAdjacent(mesh.m_ptr, element, 1, 2)

  if num_neigh_i > num_neigh_max
    num_neigh_max = num_neigh_i
  end
end
=#


# now perform the coloring
# visit each element, get colors of its neighbors
# give current element the lowest color possible
# also construct BitArray masks
#setNumberingOffset(mesh.coloring_Nptr, 1)  # set all values to -1 + 1 = 0

#=
# initialize masks to zero
for i=1:length(masks)
  masks[i] = falses(mesh.numEl)
  println("masks[i] = ", masks[i])
end
=#

for i=1:mesh.numEl
  numberJ(mesh.coloring_Nptr, mesh.elements[i], 0, 0, 0)
end

adj_size = 6  # guess number of neighboring faces
numc = 4  # guess number of colors
adj = Array(Ptr{Void}, adj_size)  # hold neighboring faces
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





# perform distance-1 coloring of mesh 
function colorMesh2(mesh::PumiMeshDG2, colordata::ColoringData)
# each element must have a different color than its neighbors with which it 
# shares and edge

println(mesh.f, "----- Entered colorMesh2 -----")

for i=1:mesh.numEl
  numberJ(mesh.coloring_Nptr, mesh.elements[i], 0, 0, 0)
end

numc = 4  # guess number of colors

adj = Array(Ptr{Void}, 3)  # distance-1 edge neighbors
adj2 = Array(Array{Ptr{Void}, 1}, 3)  # distance-2 edge neighbors + distance-1 

for i=1:3
  adj2[i] = Array(Ptr{Void}, 4)
end

colors = zeros(Int32, 3*4)

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
 
  getDistance2Colors(mesh, i, adj, adj2, colors)

#  println("colors = \n", colors)
#  num_adj = countBridgeAdjacent(mesh.m_ptr, el_i, 0, 2)

#=
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
=#
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

function colorMeshBoundary2(mesh::PumiMeshDG2, colordata::ColoringData, numc, cnt_colors)

  println(mesh.f, "----- Entered colorMeshBoundary2 -----")
  colordata.nonlocal_colors = zeros(Int32, mesh.numSharedEl)
  adj = Array(Ptr{Void}, 3)  # distance-1 edge neighbors
  adj2 = Array(Array{Ptr{Void}, 1}, 3)  # distance-2 edge neighbors + distance-1
  self = Array(Ptr{Void}, 1)  # the current element pointer
  for i=1:3
    adj2[i] = Array(Ptr{Void}, 4)
  end

  const local_start = 1
  const nonlocal_start = 12

  colors = zeros(Int32, 5*4 + 2)  
  # the first 3 sets of four elements are for the colors of the distance-2 
  # neighbors + their common distance-2 neighbor
  # the last 2*3 entries are for the non-local neighbors of the distance-1 
  # neighbors
  # the last two entries are for the local distance-2 neighbors connected by a
  # non-local element
  # (for the case where a single element has 2 non local neighbors)
  local_colors = view(colors, 1:12)
  nonlocal_d1neighborcolors = view(colors, 13:18)
  nonlocal_neighborcolors = view(colors, 20:22)


  for i in keys(colordata.adj_dict)
    el_i = mesh.elements[i]
    self[1] = el_i

    num_adj = getDistance2Colors(mesh, i, adj, adj2, local_colors)
    getNonLocalColors(mesh, view(adj, 1:num_adj), colordata, nonlocal_d1neighborcolors)
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
  # it is preferable to do these last because they have strictly fewer adjacent elements than
  # other elemenets, so the odds are better they can get a lower color by doing coloring all
  # the boundary elements first

  revadj = colordata.revadj
  colors = zeros(Int32, 3*4)
  first_neighborcolors = view(colors, 1:4)
  second_neighborcolors = view(colors, 5:8)
  first_nonlocal_neighborcolors = view(colors, 9:10)
  second_nonlocal_neighborcolors = view(colors, 11:12)

  d1_ptr = Array(Ptr{Void}, 1)
  for i=1:mesh.numSharedEl

    # get the neighbors
    neighbor1 = revadj[i, 1]
    neighbor2 = revadj[i, 2]
    getDistance1Colors(mesh, neighbor1, adj, first_neighborcolors)
    d1_ptr[1] = mesh.elements[neighbor1]
    first_neighborcolors[4] = getNumberJ(mesh.coloring_Nptr, d1_ptr[1], 0, 0)
    getNonLocalColors(mesh, d1_ptr, colordata, first_nonlocal_neighborcolors)

    if neighbor2 != 0
      getDistance1Colors(mesh, neighbor2, adj, second_neighborcolors)
      d1_ptr[1] = mesh.elemenets[neighbor2]
      second_neighborcolors[4] = getNumberJ(mesh.coloring_Nptr, d1_ptr[1], 0, 0)
      getNonLocalColors(mesh, d1_ptr, colordata, second_nonloca_neighborcolors)
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


function getDistance2Colors(mesh::PumiMeshDG2, elnum::Integer, adj, adj2, colors)
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

  num_adj = countBridgeAdjacent(mesh.m_ptr, el_i, 1, 2)
#  adj = Array(Ptr{Void}, 3)
  getBridgeAdjacent(adj)

#  adj2 = Array(Array{Ptr{Void}, 1}, num_adj)  # hold distance 2 neighbors

  adj_cnt = zeros(Int, num_adj)
  for j=1:num_adj
    num_adj_j = countBridgeAdjacent(mesh.m_ptr, adj[j], 1, 2)
    adj_cnt[j] = num_adj_j + 1
 #   adj2[j] = Array(Ptr{Void}, num_adj_j + 1)
    getBridgeAdjacent(adj2[j])
    adj2[j][num_adj_j + 1] = adj[j]  # include the distance-1 neighbor
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

function getNonlocalDistance2Colors(mesh::PumiMeshDG2, elnum::Integer, colordata::ColoringData, colors)
# get the distance 2 neighbors of a local element that are connected via a non-local element

  el = mesh.elements[elnum]
  pos = 1
  if haskey(colordata.adj_dict, elnum)
    vals = colordata.adj_dict[elnum]
    for i=1:length(vals)
      if (vals[i] != 0)
        for j=1:2
          neighbor = colordata.revadj[vals[i] - mesh.numEl, j]
          if neighbor != 0
            neighbor_ptr = mesh.elements[neighbor]
            colors[pos] = getNumberJ(mesh.coloring_Nptr, neighbor_ptr, 0, 0)
          end
        end
      end
    end
  end

  return nothing
end
            

function getDistance1Colors(mesh::PumiMeshDG2, elnum::Integer, adj, colors)
# get the distance1 colors of a specified element

  el = mesh.elements[elnum]
  num_adj = countBridgeAdjacent(mesh.m_ptr, el, 1, 2)
  getBridgeAdjacent(adj)

  for i=1:num_adj
    colors[i] = getNumberJ(mesh.coloring_Nptr, adj[i], 0, 0)
  end

  return nothing
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



function getPertEdgeNeighbors(mesh::PumiMeshDG2)

  neighbor_nums = zeros(Int32, mesh.numEl, 3)

  edges = Array(Ptr{Void}, 3)  # hold element edges
  adj = Array(Ptr{Void}, 2)  # hold adjacent elements

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


function getColors0(mesh, masks::AbstractArray{BitArray{1}, 1})
# populate the masks that indicate which elements are perturbed by which
# colors
# for a distance-0 coloring, every element is perturbed and there is only
# 1 color

  masks[1] = trues(mesh.numEl)
  return nothing
end


function getColors1{T, T2}(mesh, colordata::ColoringData, masks::AbstractArray{BitArray{1}, 1}, neighbor_colors::AbstractArray{T, 2}, neighbor_nums::AbstractArray{T2, 2}; verify=true )
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

adj = Array(Ptr{Void}, 4)   # pointers to element i + 3 neighbors
adj_color = zeros(Int32, 4)  # element colors
adj_elnum = zeros(Int32, 4)  # element numbers
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
  num_adj = countBridgeAdjacent(mesh.m_ptr, el_i, 1, 2)
  @assert num_adj <= 3
  getBridgeAdjacent(adj)  
  adj[num_adj + 1] = el_i  # insert the current elementa

  # get color, element numbers for non-local elements
  pos = 1  #current index in the adj_color, adj_elnum arrays
  if haskey(adj_dict, elnum)
  nonlocal_els = adj_dict[elnum]
    for j=1:2
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
  for j=1:4
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
nonlocal_masks = Array(Array{BitArray{1}, 1}, mesh.npeers)
start_elnum = 1
for i=1:mesh.npeers
  nonlocal_masks[i] = Array(BitArray{1}, mesh.numColors)
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


function getMinColor{T}(adj::AbstractArray{T})
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

function getPertNeighbors0(mesh::PumiMeshDG2)
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


function getPertNeighbors1(mesh::PumiMeshDG2)
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





# perform a distance-0 coloring of the mesh (ie. all elements same color)
function colorMesh0(mesh::PumiMesh2)

  # make every element color 1
  for i=1:mesh.numEl
    numberJ(mesh.coloring_Nptr, mesh.elements[i], 0, 0, 1)
  end

  return 1  # return number of colors
end



# perform distance-1 coloring of mesh 
# not sure if this works correctly
function colorMesh1(mesh::PumiMesh2, masks::Array{BitArray{1}})
# each element must have a different color than its neighbors with which it 
# shares and edge

# figure out the number of colors
# this is a lot of random memory access
#=
num_neigh_max = 0
for i=1:mesh.numEl  # loop over elements
  el_i = mesh.elements[i]
  num_neigh_i = countBridgeAdjacent(mesh.m_ptr, element, 1, 2)

  if num_neigh_i > num_neigh_max
    num_neigh_max = num_neigh_i
  end
end
=#


# now perform the coloring
# visit each element, get colors of its neighbors
# give current element the lowest color possible
# also construct BitArray masks
#setNumberingOffset(mesh.coloring_Nptr, 1)  # set all values to -1 + 1 = 0

#=
# initialize masks to zero
for i=1:length(masks)
  masks[i] = falses(mesh.numEl)
  println("masks[i] = ", masks[i])
end
=#

for i=1:mesh.numEl
  numberJ(mesh.coloring_Nptr, mesh.elements[i], 0, 0, 0)
end

adj_size = 6  # guess number of neighboring faces
numc = 4  # guess number of colors
adj = Array(Ptr{Void}, adj_size)  # hold neighboring faces
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





# perform distance-1 coloring of mesh 
function colorMesh2(mesh::PumiMesh2)
# each element must have a different color than its neighbors with which it 
# shares and edge

# figure out the number of colors
# this is a lot of random memory access
#=
num_neigh_max = 0
for i=1:mesh.numEl  # loop over elements
  el_i = mesh.elements[i]
  num_neigh_i = countBridgeAdjacent(mesh.m_ptr, element, 1, 2)

  if num_neigh_i > num_neigh_max
    num_neigh_max = num_neigh_i
  end
end
=#


# now perform the coloring
# visit each element, get colors of its neighbors
# give current element the lowest color possible
# also construct BitArray masks
#setNumberingOffset(mesh.coloring_Nptr, 1)  # set all values to -1 + 1 = 0

#=
# initialize masks to zero
for i=1:length(masks)
  masks[i] = falses(mesh.numEl)
  println("masks[i] = ", masks[i])
end
=#

for i=1:mesh.numEl
  numberJ(mesh.coloring_Nptr, mesh.elements[i], 0, 0, 0)
end

#adj_size = 6  # guess number of neighboring faces
numc = 4  # guess number of colors
#adj = Array(Ptr{Void}, adj_size)  # hold neighboring faces
#adj_color =zeros(Int32, adj_size)  # colors of neighboring faces

adj = Array(Ptr{Void}, 3)  # distance-1 edge neighbors
adj2 = Array(Array{Ptr{Void}, 1}, 3)  # distance-2 edge neighbors + distance-1 

for i=1:3
  adj2[i] = Array(Ptr{Void}, 4)
end

colors = zeros(Int32, 3*4)

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
 
  getDistance2Colors(mesh, i, adj, adj2, colors)

#  println("colors = \n", colors)
#  num_adj = countBridgeAdjacent(mesh.m_ptr, el_i, 0, 2)

#=
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
=#
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
println("number of colors = ", numc)
println("number of each color = ", cnt_colors)
mesh.color_cnt = cnt_colors

return numc

end


function getDistance2Colors(mesh, elnum::Integer, adj, adj2, colors)
# get the distance-2 neighbors of a given element
# repeats included
# adj : array to be populated with distance-1 edge neighbors
# adj2 : array to be populated with sum of distance-1 and distance-2 edge
# neighbors
# colors : array to be populated with colors of elements in adj2

  el_i = mesh.elements[elnum]

  num_adj = countBridgeAdjacent(mesh.m_ptr, el_i, 1, 2)
#  adj = Array(Ptr{Void}, 3)
  getBridgeAdjacent(adj)

#  adj2 = Array(Array{Ptr{Void}, 1}, num_adj)  # hold distance 2 neighbors

  adj_cnt = zeros(Int, num_adj)
  for j=1:num_adj
    num_adj_j = countBridgeAdjacent(mesh.m_ptr, adj[j], 1, 2)
    adj_cnt[j] = num_adj_j + 1
 #   adj2[j] = Array(Ptr{Void}, num_adj_j + 1)
    getBridgeAdjacent(adj2[j])
    adj2[j][num_adj_j + 1] = adj[j]  # include the distance-1 neighbor
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


  return nothing
end

function getPertEdgeNeighbors(mesh::PumiMesh2)

  neighbor_nums = zeros(Int32, mesh.numEl, 3)

  edges = Array(Ptr{Void}, 3)  # hold element edges
  adj = Array(Ptr{Void}, 2)  # hold adjacent elements

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



function getColors0(mesh, masks::AbstractArray{BitArray{1}, 1})
# populate the masks that indicate which elements are perturbed by which
# colors
# for a distance-0 coloring, every element is perturbed and there is only
# 1 color

  masks[1] = trues(mesh.numEl)
  return nothing
end


function getColors1{T, T2}(mesh, masks::AbstractArray{BitArray{1}, 1}, neighbor_colors::AbstractArray{T, 2}, neighbor_nums::AbstractArray{T2, 2}; verify=true )
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

adj = Array(Ptr{Void}, 4)   # pointers to element i + 3 neighbors
adj_color = zeros(Int32, 4)  # element colors

# initialize masks to 0 (false)
for i=1:length(masks)
  masks[i] = falses(mesh.numEl)
end

if verify
  println("verifying distance-1 coloring")
end


cnt = 0
for i=1:mesh.numEl
  el_i = mesh.elements[i]

  # check edge neighbors only
  num_adj = countBridgeAdjacent(mesh.m_ptr, el_i, 1, 2)
  @assert num_adj <= 3
  getBridgeAdjacent(adj)

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


    if nz_arr != unique(nz_arr)
      println("element ", i, " has non unique colors")
      println("adj_color = ", nz_arr)
      cnt += 1
    end
  end   # end if verify
  fill!(adj_color, 0)
end

if verify
  println("color-1 verification finished")
  if cnt != 0
    println(STDERR, "Warning: non unique element coloring")
  end
  println("number of element with non unique coloring = ", cnt)
end

return cnt

end





function getMinColor{T}(adj::AbstractArray{T})
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



function getMinColor2{T}(adj::AbstractArray{T}, numc::Integer)
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
  println("adding color ", numc + 1)
  min_color = numc + 1
elseif mask_sum == (numc - 1)  # there is exactly 1 color remaining
#  println("exactly 1 color remaining")
  # find out which color is missing and use it
  for i=1:numc
    if !mask[i]  # if mask is false
      min_color = i
    end
  end
else  # some colors are missing
#  println("getting minimum color")
  min_color = getMinColor(adj)  # get the minimum
end

@assert min_color != 0

return min_color

end


function getPertNeighbors0(mesh::PumiMesh2)
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


function getPertNeighbors1(mesh::PumiMesh2)
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


