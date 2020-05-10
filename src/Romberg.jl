module Romberg

using StaticArrays
using Trapz

function romberg() = integrate()

# TODO: Don't create a full R matrix - we only need to store an Rcurrent and Rprevious column
# TODO: allow mutable R argument
# TODO: by default, autoselect max_steps
function integrate(vx::AbstractRange, vy, max_steps::Integer)
    # http://www.sci.utah.edu/~beiwang/teaching/cs6210-fall-2016/lecture24.pdf
    # see also: http://www.sci.utah.edu/~beiwang/teaching/cs6210-fall-2016/lecture23.pdf

    vxlength = length(vx)

    @assert max_steps <= log2(prevpow(2, typemax(max_steps))) "`max_steps` cannot exceed `log2(prevpow(2, typemax(max_steps)))`"
    @assert max_steps <= log2(nextpow(2, vxlength)) "`max_steps` cannot exceed `log2(nextpow(2, vxlength))`"

    @assert ispow2(vxlength-1) "`length(vx) - 1` must be a power of 2"
    @assert vxlength == length(vy) "`length(vx)` must equal `length(vy)`"
    # NOTE: by requiring vx::AbstractRange, we can guarantee fixed step size

    R = zeros(max_steps,max_steps)
    for i = 1:max_steps
        # `div` is safe because we assert ispow2
        s = 1:div(vxlength - 1, 2^(i-1)):vxlength

        # R[i,1] = trapz(vx[s], vy[s])
        R[i,1] = trapz(vx[s], view(vy, s))
    end

    for j = 2:max_steps
        # Precompute reused values
        jm1 = j - 1
        pow4_jm1 = 4^jm1
        den = 1/(pow4_jm1 - 1)

        for i = j:max_steps
            R[i,j] = (pow4_jm1*R[i,jm1] - R[i-1,jm1]) * den
        end
    end

    return LowerTriangular(R)
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
