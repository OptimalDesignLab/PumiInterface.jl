# test some auxiliary math functions

function test_math()

  @testset "Auxiliary Math" begin
    h = 1e-20
    pert = Complex128(0, h)

    A = rand_realpart((2, 2))
    test_det(A)
    test_adjugate(A)

    A = rand_realpart((3, 3))
    test_det(A)
    test_adjugate(A)

  end

  return nothing
end

function test_det(A::AbstractMatrix)

  m = size(A, 1)
  n = size(A, 2)
  h = 1e-20
  pert = Complex128(0, h)


  A_dot = zeros(Complex128, m, n)
  A_dot2 = zeros(Complex128, m, n)
  # test det
  for i=1:m
    for j=1:n
      A[i, j] += pert
      A_dot[i, j] = imag(PdePumiInterface.det2(A))/h
      A[i, j] -= pert
    end
  end

  PdePumiInterface.det2_rev(A, A_dot2, 1.0)

  @test maximum(abs.(A_dot - A_dot2)) < 1e-13

  return nothing
end

function test_adjugate(A::AbstractMatrix)

  m = size(A, 1)
  n = size(A, 2)
  h = 1e-20
  pert = Complex128(0, h)



  A_dot = zeros(Complex128, m, n, m, n)
  A_dot2 = zeros(Complex128, m, n, m, n)
  tmp = zeros(Complex128, m, n)
  tmp_bar = zeros(Complex128, m, n)
  for i=1:m
    for j=1:n
      A[i, j] += pert
      if m == 2
        PdePumiInterface.adjugate2(A, tmp)
      else
        PdePumiInterface.adjugate3(A, tmp)
      end
      A_dot[:, :, i, j] = imag(tmp)/h
      A[i, j] -= pert
    end
  end

  fill!(tmp, 0)
  A_bar = zeros(Complex128, m, n)
  for i=1:m
    for j=1:n
      tmp_bar[i, j] = 1
      fill!(A_bar, 0)
      if m == 2
        PdePumiInterface.adjugate2(A, tmp)
        PdePumiInterface.adjugate2_rev(A, A_bar, tmp, tmp_bar)
      else
        PdePumiInterface.adjugate3(A, tmp)
        PdePumiInterface.adjugate3_rev(A, A_bar, tmp, tmp_bar)
      end

      A_dot2[i, j, :, :] = A_bar
      tmp_bar[i, j] = 0
    end
  end

  println("A_dot = \n", real(A_dot))
  println("A_dot2 = \n", real(A_dot2))
  @test maximum(abs.(A_dot - A_dot2)) < 1e-13

  return nothing
end
