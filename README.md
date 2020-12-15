# Romberg.jl

[![Build Status](https://travis-ci.com/fgasdia/Romberg.jl.svg?branch=master)](https://travis-ci.com/fgasdia/Romberg.jl)

A simple Julia package to perform Romberg integration over discrete 1-dimensional
data.

[Romberg integration](https://en.wikipedia.org/wiki/Romberg's_method) combines trapezoidal integration with Richardson extrapolation for improved accuracy. This package is
meant to be used for integrating discrete data sampled with equal spacing. If
your integrand can be evaluated at arbitrary points, then other methods are probably a better
choice, e.g. [QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl). Similarly,
if you have discrete data sampled at generic _unequally_ spaced points, you probably
need to use a low-order method like [Trapz.jl](https://github.com/francescoalemanno/Trapz.jl).

## Usage

First, install with

```jl
] add https://github.com/fgasdia/Romberg.jl
```

The Romberg module exports a single function, `romberg(x,y)`, or alternatively `romberg(Δx,y)`,
that returns a tuple `(I,E)` of the estimated integral `I` and a rough upper bound `E` on
the error.
```jl
using Romberg

x = range(0, pi, length=2^8+1)
y = sin.(x)

romberg(x, y)
```

## Limitations

### Equally spaced `x`

Unlike [Trapz.jl](https://github.com/francescoalemanno/Trapz.jl), Romberg
integration is a [Newton-Cotes](https://en.wikipedia.org/wiki/Newton%E2%80%93Cotes_formulas)
formula which requires each element of `x` be equally spaced. This is indirectly
enforced in `Romberg` by requiring `x::AbstractRange`, or directly by passing the
spacing `Δx` between points.

### Most effective if `length(x) - 1` a power of 2

Romberg integration works by recursively breaking the integral down into
trapezoidal-rule evaluations using larger and larger spacings `Δx` and then
extrapolating back towards `Δx → 0`.   This works by factorizing `length(x) - 1`,
and therefore works best when `length(x) - 1` has **many small factors**, ideally
being a power of two.

(In the even that `length(x) - 1` is prime, the `romberg` function is nearly
equivalent to the trapezoidal rule, since it extrapolates only from 2 points to
`length(x)` points.)

### 1-dimensional

Currently `Romberg` only allows integration over a single dimension, so
`y::AbstractVector`.

## Comparison to `trapz`

Given the limitations of `Romberg`, why use it over `Trapz`? For discrete
samples of an underlying smooth function, Romberg integration can obtain
_significantly more accurate_ estimates at relatively _low additional
computational cost_ over trapezoidal integration for a given number of samples.
Moreover, unlike the trapezoidal rule, the `romberg` function also *returns an error estimate*.

Here are some examples:

### 1)

![\displaystyle \int_0^\pi \sin(x) \,\mathrm{d}x = 2](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle%20%5Cint_0%5E%5Cpi%20%5Csin(x)%20%5C%2C%5Cmathrm%7Bd%7Dx%20%3D%202)

```jl
using BenchmarkTools
using Trapz, Romberg

x = range(0, π, length=2^6+1)
y = sin.(x)
exact_answer = 2

tans = trapz(x, y)
rans, _ = romberg(x, y)
```

```jl
julia> exact_answer - tans
0.0004016113599627502

julia> exact_answer - rans
4.440892098500626e-16
```

```jl
julia> @btime trapz($x, $y);
  340.834 ns (1 allocation: 96 bytes)

julia> @btime romberg($x, $y);
  515.078 ns (1 allocation: 192 bytes)
```

So `romberg` is ~50% slower than `trapz`, but achieves nearly machine-precision accuracy,
~12 digits more accurate than `trapz`. Even if 500 times as many samples of the
function were to be used in `trapz`, it would still be ~7 digits less accurate than `romberg`.

### 2)

![\displaystyle \int_0^1 x^3 \,\mathrm{d}x = 0.25](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle%20%5Cint_0%5E1%20x%5E3%20%5C%2C%5Cmathrm%7Bd%7Dx%20%3D%200.25)

```jl
x = range(0, 1, length=2^4+1)
y = x.^3
exact_answer = 0.25

tans = trapz(x, y)
rans, _ = romberg(x, y)
```

```jl
julia> exact_answer - tans
-0.0009765625

julia> exact_answer - rans
0.0
```

`romberg` is able to obtain the exact answer (and in general is exact for polynomials
of sufficiently low degree), compared to ~3 digits of accuracy
for `trapz`, at the cost of ~2× the run time.

### 3)

![\displaystyle \int_0^\pi \sin(mx)\cos(nx) \,\mathrm{d}x = \frac{2m}{m^2 - n^2}, \quad mn \; \text{odd}](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle%20%5Cint_0%5E%5Cpi%20%5Csin(mx)%5Ccos(nx)%20%5C%2C%5Cmathrm%7Bd%7Dx%20%3D%20%5Cfrac%7B2m%7D%7Bm%5E2%20-%20n%5E2%7D%2C%20%5Cquad%20mn%20%5C%3B%20%5Ctext%7Bodd%7D)

```jl
m = 3
n = 4
x = range(0, π, length=2^6+1)
y = sin.(m*x).*cos.(n*x)
exact_answer = 2*m/(m^2 - n^2)

tans = trapz(x, y)
rans, _ = romberg(x, y)
```

```jl
julia> exact_answer - tans
0.0012075513578178043

julia> exact_answer - rans
-1.2385595582475872e-7
```

`romberg` is ~4 digits better accuracy than `trapz` and ~50% greater run time.
