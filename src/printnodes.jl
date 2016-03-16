include("nodecalc.jl")

order = parse(ARGS[1])
sbp = TriSBP{Float64}(degree=order, reorder=false, internal=true)
ref_verts = [0. 1 0; 0 0 1]

xi, coords = nodecalc(sbp, true)

nnodes = sbp.numnodes
for i=1:nnodes
  println("node $i xi coordinates = ", xi[:, i])
end
