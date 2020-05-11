"""
The `Romberg` module provides the functions `romberg` and `romberg!` that
perform a Romberg integration of discrete 1-dimensional data `y` sampled at
equally spaced points over `x`.

Romberg integration combines the trapezoid method of integration with Richardson
extrapolation. Significantly lower error estimates of definite integrals can be
obtained using Romberg integration over the trapezoid method alone, at only
slightly higher computational cost.
"""
module Romberg

using Trapz

export romberg, romberg!

maxsteps(N) = trunc(Int, log2(prevpow(2, N)))

"""
    romberg(x::AbstractRange, y::AbstractVector)

Integrate `y` sampled at `x` by Romberg integration applying the maximum possible
number of Richardson extrapolation steps.

Romberg integration requires equally spaced points. This is enforced by
requiring `x` to be an `AbstractRange`, rather than a dense vector.
An additional requirement is that `length(x) == 2ⁿ + 1` for any positive integer
`n`.

# Examples

```jldoctest
julia> x = range(0, π, length=2^8+1);
julia> romberg(x, sin.(x))
1.9999999999999996
```
"""
function romberg(x::AbstractRange, y::AbstractVector)
    N = length(x)

    @boundscheck begin
        ispow2(N-1) || throw(DomainError(length(x), "`length(x) - 1` must be a power of 2"))
        N == length(y) || throw(DimensionMismatch("length of `y` not equal to length of `x`"))

        # NOTE: by requiring x::AbstractRange, a fixed step size is guaranteed
    end

    # Automatically select max_steps
    max_steps = maxsteps(N)
    @assert max_steps <= log2(prevpow(2, N)) "`max_steps` cannot exceed `log2(prevpow(2, length(x)))`"

    return integrate(x, y, max_steps)
end

"""
    romberg(x::AbstractRange, y::AbstractVector, max_steps::Integer)

Integrate `y` sampled at `x` by Romberg integration applying `max_steps` of
Richardson extrapolation.

`max_steps` is the number of Richardson extrapolation steps applied (also equal
to the number of trapezoid integrations - 1). In general, the larger it is, the
lower error in the integration. The largest `max_steps` can be is
`log2(prevpow(2, length(x)))`.

# Examples

```jldoctest
julia> x = range(0, π, length=2^8+1);
julia> romberg(x, sin.(x), 8)
1.9999999999999996
```

```jldoctest
julia> romberg(x, sin.(x), 9)
ERROR: DomainError with 9:
`max_steps` cannot exceed `log2(prevpow(2, length(x)))` = 8
[...]
```
"""
function romberg(x::AbstractRange, y::AbstractVector, max_steps::Integer)
    N = length(x)

    @boundscheck begin
        max_steps <= log2(prevpow(2, N)) || throw(DomainError(max_steps, "`max_steps` cannot exceed `log2(prevpow(2, length(x)))` = $(Int(log2(prevpow(2, length(x)))))"))
        ispow2(N-1) || throw(DomainError(length(x), "`length(x) - 1` must be a power of 2"))
        N == length(y) || throw(DimensionMismatch("length of `y` not equal to length of `x`"))

        # NOTE: by requiring x::AbstractRange, a fixed step size is guaranteed
    end

    return integrate(x, y, max_steps)
end

"""
    romberg!(R::AbstractMatrix, x::AbstractRange, y::AbstractVector)

Integrate `y` sampled at `x` by Romberg integration, filling in the square
matrix `R` in place.

The size of `R` determines the number of Richardson extrapolation steps. If `R`
has size 6×6, 5 extrapolation steps will be performed. In general, for a matrix
with dimensions `L×L`, `max_steps` corresponds to `L - 1`.

This function is useful to see how the Romberg integration progressed.

# Examples

```jldoctest
julia> x = range(0, π, length=2^8+1);
julia> R = zeros(4, 4);
julia> romberg!(R, x, sin.(x))
4×4 Array{Float64,2}:
 1.92367e-16  0.0      0.0      0.0
 1.5708       2.0944   0.0      0.0
 1.89612      2.00456  1.99857  0.0
 1.97423      2.00027  1.99998  2.00001
```
"""
function romberg!(R::AbstractMatrix, x::AbstractRange, y::AbstractVector)
    N = length(x)

    @boundscheck begin
        # Assume `size(R, 1)` is L (`max_steps + 1`)
        size(R,1) <= log2(prevpow(2, N)) || throw(DomainError(size(R,1), "dimenensions of `R` cannot be greater than `log2(prevpow(2, length(x)))`"))
        isequal(size(R)...) || throw(DimensionMismatch("`R` must be square"))
        N == length(y) || throw(DimensionMismatch("length of `y` not equal to length of `x`"))

        # NOTE: by requiring x::AbstractRange, a fixed step size is guaranteed
    end

    integrate!(R, x, y)
end

@inline function integrate(x::AbstractRange{Tx}, y::AbstractVector{Ty}, max_steps::Integer) where {Tx,Ty}
    # `max_steps` are extrapolation steps, size of `R` is `max_steps + 1`
    L = max_steps + 1

    # `Rp` is "Rprevious", `Rc` is "Rcurrent"
    Rp = Array{promote_type(Tx,Ty)}(undef, L)
    Rc = similar(Rp)

    trapezoid!(Rp, x, y)
    extrapolate!(Rp, Rc)

    # best estimate
    return Rc[end]
end

@inline function integrate!(R::AbstractMatrix, x::AbstractRange, y::AbstractVector)
    L = size(R, 1)

    trapezoid!(view(R,:,1), x, y)
    extrapolate!(R)

    return R
end

@inline function trapezoid!(R, x, y)
    N = length(x)
    @inbounds for i in eachindex(R)
        idxs = 1:div(N-1, 2^(i-1)):N

        R[i] = trapz(x[idxs], view(y,idxs))
    end
end

@inline function extrapolate!(Rp, Rc)
    L = length(Rp)
    @inbounds for j = 2:L
        # Precompute reused values
        pow_prev = 4^(j-1)
        den = inv(pow_prev - 1)

        @inbounds for i = j:L
            Rc[i] = (pow_prev*Rp[i] - Rp[i-1]) * den
        end

        # Only copying the relevant values! (not necessary on last loop)
        if j < L
            copyto!(Rp, j, Rc, j, L-j+1)
        end
    end
end

@inline function extrapolate!(R)
    L = size(R,1)
    @inbounds for j = 2:L
        # Precompute reused values
        pow_prev = 4^(j-1)
        den = inv(pow_prev - 1)

        @inbounds for i = j:L
            R[i,j] = (pow_prev*R[i,j-1] - R[i-1,j-1]) * den
        end
    end
end

function error_test(R, max_steps)
    # See http://www.sci.utah.edu/~beiwang/teaching/cs6210-fall-2016/lecture24.pdf

    # The ratio of the difference between successive entries in column `j`
    # should be approximately 4^j.

    ratio = zeros(max_steps-2, max_steps-2)
    for j = 1:max_steps-2
        ratio[j:end,j] .= (R[1+j:end-1,j] - R[j:end-2,j]) ./ (R[2+j:end,j] - R[1+j:end-1,j])
    end

    return ratio
end

end # module
