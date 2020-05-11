using Romberg
using Test
using Trapz

@testset "Romberg.jl" begin
    # Test interfaces
    x = range(0, π, length=2^8+1)
    y = sin.(x)

    @test romberg(x, y) == romberg(x, y, 8)
    @test_throws AssertionError romberg(x, y, 9)

    max_steps = 7
    R = zeros(max_steps+1, max_steps+1)
    @test romberg!(R, x, y)[end] == romberg(x, y, max_steps)

    x = range(0, π, length=2^8)
    y = sin.(x)
    @test_throws AssertionError romberg(x, y)


    # Integrate different functions
    x = range(0, π, length=2^8+1)
    y = sin.(x)
    @test romberg(x, y) ≈ 2

    x = range(0, 1, length=2^8+1)
    y = x.^3
    @test romberg(x, y) ≈ 0.25

    x = range(0, π/2, length=2^5+1)
    y = sin.(x).^2
    @test romberg(x, y) ≈ π/4

    # this one requires lots of samples...
    a = 15
    x = range(0, a, length=2^16+1)
    y = sqrt.(a^2 .- x.^2)
    @test romberg(x, y) ≈ π*a^2/4

    # tricky to integrate
    x = range(1e-15, 1, length=2^16+1)
    y = log.(x)./(1 .+ x)
    @test romberg(x, y) ≈ -π^2/12 atol=1e-4
end
