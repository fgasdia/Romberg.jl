# Romberg.jl

[![Build Status](https://travis-ci.com/fgasdia/Romberg.jl.svg?branch=master)](https://travis-ci.com/fgasdia/Romberg.jl)

**Note** [NumericalIntegration.jl](https://github.com/dextorious/NumericalIntegration.jl) offers a faster Romberg integration routine as `integrate(x, y, RombergEven())`. 

A simple Julia package to perform Romberg integration over discrete 1-dimensional
data.

[Romberg integration](https://en.wikipedia.org/wiki/Romberg's_method) combines trapezoidal integration with Richardson extrapolation for improved accuracy. This package is
meant to be used for integrating discrete data sampled with equal spacing. If
your integrand is in functional form, then other methods are probably a better
choice, e.g. [QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl) or
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
`log2(prevpow(2, length(x)))`:
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

### Equally spaced `x`

Unlike [Trapz.jl](https://github.com/francescoalemanno/Trapz.jl), Romberg
integration is a [Newton-Cotes](https://en.wikipedia.org/wiki/Newton%E2%80%93Cotes_formulas)
formula which requires each element of `x` be equally spaced. This is indirectly
enforced in `Romberg` by requiring `x::AbstractRange`.

### Most efficient if `length(x) - 1` is a power of 2

The Romberg integration algorithm is coded assuming `ispow2(length(x) - 1) == true`.
Therefore, `Romberg` is most efficient if these criteria can be met. However,
`romberg(x,y)` can be used for `x` of any length at reduced accuracy and increased
runtime. This works by applying Romberg integration piecewise across sections of
`x` that do meet the criteria.

### 1-dimensional

Currently `Romberg` only allows integration over a single dimension, so
`y::AbstractVector`.

## Comparison to `trapz`

Given the limitations of `Romberg`, why use it over `Trapz`? For discrete
samples of an underlying smooth function, Romberg integration can obtain
_significantly higher accuracy_ estimates at relatively _low additional
computational cost_ over trapezoidal integration for a given number of samples.

Here are some examples:

### 1)

![\displaystyle \int_0^\pi \sin(x) \,\mathrm{d}x = 2](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle%20%5Cint_0%5E%5Cpi%20%5Csin(x)%20%5C%2C%5Cmathrm%7Bd%7Dx%20%3D%202)

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
julia> b = @benchmarkable trapz($x, $y);

julia> run(b)
BenchmarkTools.Trial:
  memory estimate:  96 bytes
  allocs estimate:  1
  --------------
  minimum time:     276.000 ns (0.00% GC)
  median time:      290.000 ns (0.00% GC)
  mean time:        301.079 ns (0.00% GC)
  maximum time:     20.295 μs (0.00% GC)
  --------------
  samples:          10000
  evals/sample:     1
```

```jl
b = @benchmarkable romberg($x, $y);

julia> run(b)
BenchmarkTools.Trial:
  memory estimate:  1.38 KiB
  allocs estimate:  16
  --------------
  minimum time:     1.354 μs (0.00% GC)
  median time:      1.415 μs (0.00% GC)
  mean time:        1.587 μs (0.00% GC)
  maximum time:     42.352 μs (0.00% GC)
  --------------
  samples:          10000
  evals/sample:     1
```

```jl
R = zeros(7, 7);
b = @benchmarkable romberg!(r, $x, $y) setup=(r=copy($R))

julia> run(b)
BenchmarkTools.Trial:
  memory estimate:  1.09 KiB
  allocs estimate:  14
  --------------
  minimum time:     1.298 μs (0.00% GC)
  median time:      1.348 μs (0.00% GC)
  mean time:        1.439 μs (0.00% GC)
  maximum time:     33.786 μs (0.00% GC)
  --------------
  samples:          10000
  evals/sample:     1
```

So `romberg` is ~5× slower than `trapz`, but nearly at machine precision accuracy,
~10 digits more accurate than `trapz`. Even if 5 times as many samples of the
function are used in `trapz`, it's still ~8 digits accuracy behind the `romberg`.

### 2)

![\displaystyle \int_0^1 x^3 \,\mathrm{d}x = 0.25](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle%20%5Cint_0%5E1%20x%5E3%20%5C%2C%5Cmathrm%7Bd%7Dx%20%3D%200.25)

```jl
x = range(0, 1, length=2^4+1)
y = x.^3
exact_answer = 0.25

tans = trapz(x, y)
rans = romberg(x, y)
```

```jl
julia> exact_answer - tans
-0.0009765625

julia> exact_answer - rans
0.0
```

`romberg` was able to obtain the exact answer, compared to ~3 digits of accuracy
for `trapz`, at the cost of ~7× the run time.

### 3)

![\displaystyle \int_0^\pi \sin(mx)\cos(nx) \,\mathrm{d}x = \frac{2m}{m^2 - n^2}, \quad mn \; \text{odd}](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle%20%5Cint_0%5E%5Cpi%20%5Csin(mx)%5Ccos(nx)%20%5C%2C%5Cmathrm%7Bd%7Dx%20%3D%20%5Cfrac%7B2m%7D%7Bm%5E2%20-%20n%5E2%7D%2C%20%5Cquad%20mn%20%5C%3B%20%5Ctext%7Bodd%7D)

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

`romberg` is ~4 digits better accuracy than `trapz` and ~5× the run time.
