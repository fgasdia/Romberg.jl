# Romberg.jl

[![Build Status](https://travis-ci.com/fgasdia/Romberg.jl.svg?branch=master)](https://travis-ci.com/fgasdia/Romberg.jl)

A simple Julia package to perform Romberg integration over discrete 1-dimensional
data.

[Romberg integration](https://en.wikipedia.org/wiki/Romberg's_method) combines trapezoidal integration with Richardson extrapolation for improved accuracy. This package is
meant to be used for integrating discrete data sampled with equal spacing. If
your integrand is in functional form, then other methods are probably a better
choice, e.g. [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl), [QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl), or
[HCubature.jl](https://github.com/JuliaMath/HCubature.jl), among others. Similarly,
if you have discrete data sampled at _unequally_ spaced points, you will be limited
to applying pure trapezoidal integration with [Trapz.jl](https://github.com/francescoalemanno/Trapz.jl).

## Usage

First, install with

```jl
]add https://github.com/fgasdia/Romberg.jl
```

Then, the interface is somewhat similar to `trapz`, except the primary function
to call is `romberg`:
```jl
x = range(0, pi, length=2^8+1)
y = sin.(x)

romberg(x, y)
```

By default, `romberg(x, y)` will use the maximum number of extrapolation steps
possible given the length of `x`.

If lower accuracy is acceptable or shorter compute time required, the number of
extrapolation steps can be specified as an additional argument, up to
`log2(nextpow(2, length(x)))`:
```jl
romberg(x, y, 6)
```

Finally, an in place form is available for filling in the entire triangular
array formed by the Romberg integration. This can be used to analyze the
convergence of the integration. The size of the mutable _square_ matrix argument
`R` determines the number of extrapolation steps. If `size(R)` is `(L, L)`, then
`L - 1` extrapolation steps are taken:
```jl
R = zeros(6, 6)
romberg!(R, x, y)
```
gives the output
```jl
6×6 Array{Float64,2}:
 1.92367e-16  0.0      0.0      0.0      0.0  0.0
 1.5708       2.0944   0.0      0.0      0.0  0.0
 1.89612      2.00456  1.99857  0.0      0.0  0.0
 1.97423      2.00027  1.99998  2.00001  0.0  0.0
 1.99357      2.00002  2.0      2.0      2.0  0.0
 1.99839      2.0      2.0      2.0      2.0  2.0
```

The best estimate is the lower right corner of the matrix, accessible through
linear indexing as
```jl
julia> R[end]
2.0000000000013207
```

Note that the forms of `romberg` that do _not_ have the in place argument do not
internally preallocate a full dense square array `R`. However, the same total number
of mutations (allocations) must occur.

## Limitations

***READ THIS SECTION!***

### Equally spaced `x`

Unlike [Trapz.jl](https://github.com/francescoalemanno/Trapz.jl), Romberg
integration is a [Newton-Cotes](https://en.wikipedia.org/wiki/Newton%E2%80%93Cotes_formulas)
formula which requires each element of `x` be equally spaced. This is indirectly
enforced in `Romberg` by requiring `x::AbstractRange`.

### `length(x) - 1` must be a power of 2

The Romberg integration algorithm is coded assuming `ispow2(length(x) - 1)`.
This is obviously a major downside of the implementation. It is _not_ generally
recommended to simply expand `x` and pad `y` with zeros in order to meet this
criteria because the extrapolation step assumes the integrand has at least a few
continuous derivatives. Unless the function being integrated goes to zero at the
lower or upper limit, then padding either side of `y` with zeros will likely
cause a discontinuity and affect the accuracy of the result.

### 1-dimensional

Currently `Romberg` only allows integration over a single dimension, so
`y::AbstractVector`.

## Comparison to `trapz`

Given the limitations of `Romberg`, why use it over `Trapz`? For discrete
samples of an underlying smooth function, Romberg integration can obtain
_significantly higher accuracy_ estimates at relatively _little additional
computational cost_ over trapezoidal integration.

Here are some examples:

<p align="center"><img src="/tex/bc7e8af11561e520460810c128417ea8.svg?invert_in_darkmode&sanitize=true" align=middle width=121.86696059999998pt height=38.242408049999995pt/></p>

```jl
using BenchmarkTools
using Trapz

x = range(0, π, length=2^6+1)
y = sin.(x)
exact_answer = 2

tans = trapz(x, y)
rans = romberg(x, y)
```

```jl
julia> exact_answer - tans
0.0004016113599627502

julia> exact_answer - rans
1.3322676295501878e-15
```

```jl
julia> b = @benchmarkable trapz(<img src="/tex/db413d26a2e0e9608ff4980da96a053f.svg?invert_in_darkmode&sanitize=true" align=middle width=13.96121264999999pt height=14.15524440000002pt/>y);

julia> run(b)
BenchmarkTools.Trial:
  memory estimate:  96 bytes
  allocs estimate:  1
  --------------
  minimum time:     1.118 μs (0.00% GC)
  median time:      1.397 μs (0.00% GC)
  mean time:        1.366 μs (0.00% GC)
  maximum time:     23.480 μs (0.00% GC)
  --------------
  samples:          10000
  evals/sample:     1
```

```jl
julia> b = @benchmarkable romberg(<img src="/tex/db413d26a2e0e9608ff4980da96a053f.svg?invert_in_darkmode&sanitize=true" align=middle width=13.96121264999999pt height=14.15524440000002pt/>y);
BenchmarkTools.Trial:
  memory estimate:  1.38 KiB
  allocs estimate:  16
  --------------
  minimum time:     2.235 μs (0.00% GC)
  median time:      2.514 μs (0.00% GC)
  mean time:        2.544 μs (0.00% GC)
  maximum time:     36.324 μs (0.00% GC)
  --------------
  samples:          10000
  evals/sample:     1
```

```jl
R = zeros(7, 7);
b = @benchmarkable romberg!(r, <img src="/tex/db413d26a2e0e9608ff4980da96a053f.svg?invert_in_darkmode&sanitize=true" align=middle width=13.96121264999999pt height=14.15524440000002pt/>y) setup=(r=copy(<img src="/tex/055334325c1dc16b45893849c9ad163d.svg?invert_in_darkmode&sanitize=true" align=middle width=700.7312251499999pt height=203.6529759pt/><img src="/tex/2260d4577cdbb05f0161be1ac778faa6.svg?invert_in_darkmode&sanitize=true" align=middle width=110.25105794999999pt height=33.187449900000026pt/><img src="/tex/11fb91139daaab5729827cd894ae1587.svg?invert_in_darkmode&sanitize=true" align=middle width=700.2746222999999pt height=236.7123297pt/><img src="/tex/5b8876bf3b39873c7c21812baaeb73aa.svg?invert_in_darkmode&sanitize=true" align=middle width=224.26618335pt height=28.26507089999998pt/>$

```jl
m = 3
n = 4
x = range(0, π, length=2^6+1)
y = sin.(m*x).*cos.(n*x)
exact_answer = 2*m/(m^2 - n^2)

tans = trapz(x, y)
rans = romberg(x, y)
```

```jl
julia> exact_answer - tans
0.0012075513578178043

julia> exact_answer - rans
6.515672331675049e-7
```

`romberg` is ~4 digits better accuracy than `trapz` and ~2× the run time.
