module Romberg

using StaticArrays
using Trapz

function romberg(vx, vy, max_steps)
     integrate()
end


function integrate(x::AbstractRange{Tx}, y::AbstractVector{Ty}, max_steps::Integer) where {Tx,Ty}
    N = length(x)

    # `max_steps` are extrapolation steps, size of `R` is `max_steps + 1`
    L = max_steps + 1

    # If this function is called with `@inbounds`, these assertions will be skipped
    @boundscheck begin
        @assert max_steps <= log2(nextpow(2, N)) "`max_steps` cannot exceed `log2(nextpow(2, length(x)))`"
        @assert max_steps <= log2(prevpow(2, typemax(max_steps))) "`max_steps` cannot exceed `log2(prevpow(2, typemax(max_steps)))`"

        @assert ispow2(N-1) "`length(x) - 1` must be a power of 2"
        @assert N == length(y) "Integration over `y` is incompatible with `x`. Make sure their lengths match!"

        # NOTE: by requiring x::AbstractRange, a fixed step size is guaranteed
    end

    # `Rp` is "Rprevious", `Rc` is "Rcurrent"
    Rp = Array{promote_type(Tx,Ty)}(undef, L)
    Rc = similar(Rp)

    # `trapz` integration for `L` divisions of `x`
    @inbounds for i = 1:L
        # `div` is safe because we assert ispow2
        idxs = 1:div(N-1, 2^(i-1)):N

        Rp[i] = trapz(x[idxs], view(y,idxs))
    end

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

    # best estimate
    return Rc[end]
end


@inline function integratea(x::AbstractRange, y::AbstractVector, max_steps::Integer)
    # http://www.sci.utah.edu/~beiwang/teaching/cs6210-fall-2016/lecture24.pdf
    # see also: http://www.sci.utah.edu/~beiwang/teaching/cs6210-fall-2016/lecture23.pdf

    N = length(x)
    L = max_steps + 1

    R = zeros(L,L)
    for i = 1:L
        # `div` is safe because we assert ispow2
        idxs = 1:div(N-1, 2^(i-1)):N

        R[i,1] = trapz(x[idxs], view(y, idxs))
    end

    for j = 2:L
        # Precompute reused values
        jm1 = j - 1
        pow4_jm1 = 4^jm1
        den = 1/(pow4_jm1 - 1)

        for i = j:L
            R[i,j] = (pow4_jm1*R[i,jm1] - R[i-1,jm1]) * den
        end
    end

    return LowerTriangular(R)
end

function integrate(x::AbstractRange, y)
    # Automatically select max_steps
    max_steps = log2(nextpow(2, vxlength)) - 1
    return integrate(vx, vy, max_steps)
end

# TODO: allow mutable R argument, which then determines max_steps
function integrate!(R, vx, vy)

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
