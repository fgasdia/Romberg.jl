module Romberg

using Trapz

export romberg, romberg!

maxsteps(N) = trunc(Int, log2(prevpow(2, N)))

function romberg(x::AbstractRange, y::AbstractVector, max_steps::Integer)
    # If this function is called with `@inbounds`, these assertions will be skipped
    N = length(x)

    @boundscheck begin
        @assert max_steps <= log2(prevpow(2, N)) "`max_steps` cannot exceed `log2(nextpow(2, length(x)))-1`"
        @assert ispow2(N-1) "`length(x) - 1` must be a power of 2"
        @assert N == length(y) "Integration over `y` is incompatible with `x`. Make sure their lengths match!"

        # NOTE: by requiring x::AbstractRange, a fixed step size is guaranteed
    end

    return integrate(x, y, max_steps)
end

function romberg(x::AbstractRange, y::AbstractVector)
    N = length(x)

    @boundscheck begin
        @assert ispow2(N-1) "`length(x) - 1` must be a power of 2"
        @assert N == length(y) "Integration over `y` is incompatible with `x`. Make sure their lengths match!"

        # NOTE: by requiring x::AbstractRange, a fixed step size is guaranteed
    end

    # Automatically select max_steps
    max_steps = maxsteps(N)

    return integrate(x, y, max_steps)
end

"""
R = zeros(L,L)
"""
function romberg!(R::AbstractMatrix, x::AbstractRange, y::AbstractVector)
    N = length(x)

    @boundscheck begin
        # Assume `size(R, 1)` is L (`max_steps + 1`)
        @assert size(R,1) <= log2(prevpow(2, N))
        @assert isequal(size(R)...) "`R` must be square"
        @assert N == length(y) "Integration over `y` is incompatible with `x`. Make sure their lengths match!"

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
    N = length(x)
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
