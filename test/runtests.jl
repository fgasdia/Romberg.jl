using Romberg
using Test
using trapz

@testset "Romberg.jl" begin
    # `max_steps = 0` means no Richardson extrapolation (`trapz` only)
    @test romberg(vx, vy, 0) == trapz(vx, vy)
end
