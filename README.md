# Romberg

[![Build Status](https://travis-ci.com/fgasdia/Romberg.jl.svg?branch=master)](https://travis-ci.com/fgasdia/Romberg.jl)


# Romberg integration

# usage

generally similar to trapz

# Limitations

- 2^n + 1
- equally spaced x
  - ranges only
- 1 dimensional

# Comparison to trapz

- timing
- accuracy


# todo

- relax 2n+1 -> temp fix: expand x to nextpow(2) and fill y with 0s
- automatically stop when given tolerance is achieved?
