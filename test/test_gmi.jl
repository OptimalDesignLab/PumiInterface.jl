# test for gmi

using gmi

function test_gmi()

  @testset "GMI" begin
    fname = "meshes/airfoil.smd"
    gmi.sim_start()
    gmi.register_sim()

    g = gmi.load(fname)
    
    it = gmi._gmi_begin(g, 1)
    ge = gmi._gmi_next(g, it)
    edge_tags = Array{Cint}(0)
    while ( ge != gmi.NullModelEntity)
      @test gmi.dim(g, ge) == 1
      push!(edge_tags, gmi.tag(g, ge))

      ge = gmi._gmi_next(g, it)
    end

    gmi.gmi_end(g, it)
    sort!(edge_tags)
    @test edge_tags == Cint[5, 8]

    # test iteration on Gmi_iter
    edge_tags2 = Array{Cint}(0)
    edges = Array{gmi.ModelEntity}(0)
    for ge in gmi.Gmi_iter(g, 1)
      @test gmi.dim(g, ge) == 1
      push!(edge_tags2, gmi.tag(g, ge))
      push!(edges, ge)
    end

    _edge_tags = sort(edge_tags2)
    @test _edge_tags == Cint[5, 8]

    @test eltype(gmi.Gmi_iter) == gmi.ModelEntity

    for i=1:2
      ge = gmi.find(g, 1, edge_tags2[1])
      @test ge == edges[1]
    end

    # edge 8 is the free stream (circle)
    circle = gmi.find(g, 1, 8)
    r = zeros(Float64, 2)
    gmi.range(g, circle, 0, r)
    @test gmi.periodic(g, circle, 0)
    @test gmi.can_eval(g)
    # the circle is centered at the origin, so traversing half the parameter
    # range should rotate the point 180 degrees about the origin
    p1 = [0.25, 0.0]
    p2 = p1 + [ (r[2]-r[1])/2, 0]
    x1 = zeros(Float64, 3)
    x2 = zeros(Float64, 3)
    gmi.geval(g, circle, p1, x1)
    gmi.geval(g, circle, p2, x2)
    println("x1 = ", x1)
    println("x2 = ", x2)
    @test maximum(abs.(x1 + x2)) < 1e-13
    @test abs(norm(x1) - 10) < 1e-12
    @test abs(norm(x2) - 10) < 1e-12

    # lets see what happens when we go beyond the parametric range for a 
    # periodic entity
    # It appears for p values greater than r[2], it returns the x values at
    # r[2].  It also appears the parameter is not linear
    # p = 0 -> (10.0, 0)
    # p = 0.25 -> (0, 10.0)
    # p = 0.125 -> (0.8, 0.6) ???
    #=
    p3 = copy(p1)
    p3[1] += 1
    x3 = zeros(Float64, 3)
    gmi.geval(g, circle, p3, x3)
    println("p1 = ", p1)
    println("p3 = ", p3)
    println("x1 = ", x1)
    println("x3 = ", x3)

    p0 = zeros(2)
    x0 = zeros(3)
    gmi.geval(g, circle, p0, x0)
    println("p0 = ", p0)
    println("x0 = ", x0)
    =#

    # test reparam
    p0 = [0.125, 0.0]
    x0 = zeros(3)
    gmi.geval(g, circle, p0, x0)
    faces = gmi.adjacent(g, circle, 2)
    @test length(faces.v) == 1
    face = faces.v[1]
    
    p1 = zeros(2)
    x1 = zeros(3)
    gmi.reparam(g, circle, p0, face, p1)
    gmi.geval(g, face, p1, x1)

    @test maximum(abs.(x1 - x0)) < 1e-10  # face evaluations are not so accurate

    # test closest point, normal
    x0 = zeros(3)
    x1 = zeros(3)
    p1 = zeros(2)
    gmi.closest_point(g, circle, x0, x1, p1)
    @test (norm(x1) - 10) < 1e-10

    # normal only works for geometric faces
    p = zeros(2)
    n = zeros(3)
    gmi.normal(g, face, p, n)
    println("n = ", n)
    @test abs(norm(n) - 1) < 1e-13
    @test abs(n[3] -1) < 1e-13

    # testing first_derivative is tricky, because the parameter is not
    # theta, do finite difference instead
    p = zeros(2)
    t0 = zeros(3)
    t1 = zeros(3)
    gmi.first_derivative(g, circle, p, t0, t1)

    h = 1e-6
    p1 = [h, 0.0]
    x0 = zeros(3)
    x1 = zeros(3)
    gmi.geval(g, circle, p, x0)
    gmi.geval(g, circle, p1, x1)
    t0_fd = (x1 - x0)/h
    @test maximum(abs.(t0_fd - t0)) < 1e-5*maximum(t0_fd)

    # is_in_closure_of
    @test gmi.is_in_closure_of(g, circle, face)












    gmi.destroy(g)
  end

end
