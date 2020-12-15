"""
The `Romberg` module provides the function `romberg` that
performs a Romberg integration of discrete 1-dimensional data `y` sampled at
equally spaced points over `x`.

Romberg integration combines the trapezoid method of integration with Richardson
extrapolation. Significantly lower error estimates of definite integrals can be
obtained using Romberg integration over the trapezoid method alone, at only
slightly higher computational cost.
"""
module Romberg

export romberg

import Primes, Richardson

@views function _romberg(Δx, y, endsum, factors, numfactors; kws...)
    b, e = firstindex(y), lastindex(y)
    v = (Δx * (sum(y[b+1:e-1]) + endsum), Δx)
    vals = Vector{typeof(v)}(undef, numfactors+1)
    i = numfactors+1
    vals[i] = v
    step = 1
    for (f,K) in factors
        for k = 1:K
            step *= f
            sΔx = step*Δx
            if i == 2 # last iteration (empty sum)
                vals[1] = (sΔx * endsum, sΔx)
            else
                vals[i -= 1] = (sΔx * (sum(y[b+step:step:e-step]) + endsum), sΔx)
            end
        end
    end
    return Richardson.extrapolate!(vals, power=2; kws...)
end

function romberg(Δx::Real, y::AbstractVector; kws...)
    n = length(y)
    endsum = n ≥ 2 ? (y[begin]+y[end])/2 : (n == 1 ? zero(y[1])/1 : zero(eltype(y))/1)
    n <= 2 && return endsum * Δx
    m = n - 1
    if ispow2(m)
        # fast path: no need to allocate factors array
        k = 0
        while m > 1
            k += 1
            m >>= 1
        end
        return _romberg(float(Δx), y, endsum, (2=>k,), k; kws...)
    else
        factors = Primes.factor(m)
        return _romberg(float(Δx), y, endsum, factors, sum(values(factors)); kws...)
    end
end

romberg(x::AbstractRange, y::AbstractVector; kws...) = romberg(step(x), y; kws...)

"""
    romberg(Δx::Real, y::AbstractVector)
    romberg(x::AbstractRange, y::AbstractVector)

Integrate `y` sampled at points with equal spacing `Δx`, or equivalently
a range `x` (`Δx=step(x)`) by Romberg integration, which corresponds
to Richardson extrapolation of a trapezoidal rule to smaller and smaller `Δx`.

The return value is a tuple `(I, E)` of the estimated integral `I` and
an estimated upper bound `E` on the error in `I`.

The algorithm is most effective if `length(x) - 1` has many small factors,
ideally being a power of two, but it can handle any length.  If `length(x) - 1`
is a prime number, it is nearly equivalent to the trapezoidal rule without extrapolation.

# Examples

```jldoctest
julia> x = range(0, π, length=2^8+1);

julia> romberg(x, sin.(x))
(2.0000000000000018, 1.9984014443252818e-15)
```
"""
romberg

end # module
