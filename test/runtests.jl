using Romberg, Test

@testset "Romberg.jl" begin
    # Test interfaces
    x = range(0, π, length=2^8+1)
    y = sin.(x)

    @test romberg(x, y) == romberg(x, y, maxeval=9)

    # Test ispow2(length(x) - 1) == false
    x = range(0, π, length=2^8)
    y = sin.(x)
    @test romberg(x, y)[1] ≈ 2   rtol=1e-12

    x = range(0, π, length=2^8-1)
    y = sin.(x)
    @test romberg(x, y)[1] ≈ 2   rtol=1e-9

    # Test length(x) == 1
    x = 0:0
    y = sin.(x)
    @test romberg(x, y)[1] == 0

    # Test length(x) == 2
    x = range(0, 1, length=2)
    y = x.^2
    @test romberg(x,y)[1] == 0.5   # integration is inaccurate

    # Test length(x) == 3
    x = range(0, 1, length=3)
    y = x.^2
    @test romberg(x, y)[1] ≈ 1/3   rtol=3e-16  # should be exact up to roundoff error

    # Integrate different functions
    x = range(0, π, length=2^8+1)
    y = sin.(x)
    @test romberg(x, y)[1] ≈ 2    rtol=1e-15

    x = range(0, 1, length=2^8+1)
    y = x.^3
    @test romberg(x, y)[1] ≈ 1/4   rtol=3e-16  # should be exact up to roundoff error

    x = range(0, π/2, length=2^5+1)
    y = sin.(x).^2
    @test romberg(x, y)[1] ≈ π/4    rtol=1e-15

    m = 3
    n = 4
    x = range(0, π, length=2^8+1)
    y = sin.(m*x).*cos.(n*x)
    v = romberg(x, y)
    @test v[1] ≈ 2*m/(m^2 - n^2)  rtol=1e-13
    @test v[2] < 1e-10

    # the following are singular integrands where the underlying
    # theory of Romberg integration starts to break down:

    # this one requires lots of samples...
    a = 15
    x = range(0, a, length=2^16+1)
    y = sqrt.(a^2 .- x.^2)
    @test romberg(x, y)[1] ≈ π*a^2/4   rtol=1e-8
    # we can do much better if we put in the correct convergence
    # rate for this singular integrand, but this requires some
    # understanding of the theory that most people won't have…
    x = range(0, a, length=2^7+1)
    y = sqrt.(a^2 .- x.^2)
    @test romberg(x, y, power=0.5)[1] ≈ π*a^2/4   rtol=1e-8

    # tricky to integrate
    x = range(1e-15, 1, length=2^16+1)
    y = log.(x)./(1 .+ x)
    @test romberg(x, y)[1] ≈ -π^2/12  atol=1e-4

    # make sure it works for abstractly-typed y and integer Δx
    @test romberg(1, Any[3,3,3,3,3,3,3])[1] == 6*3

    @test @inferred(romberg(1, [1//2, 1//4])) === (0.375, 0.375)
    @test @inferred(romberg(1, [[1//2,1//1], [1//4,1//1]])) == ([0.375,1.0], hypot(0.375,1.0))
end
