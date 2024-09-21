# AstroCoords

[![Build Status](https://github.com/jmurphy6895/AstroCoords.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/jmurphy6895/AstroCoords.jl/actions/workflows/CI.yml?query=branch%3Amaster)


AstroCoords.jl
================================

This package contains various propagators of satellite trajectories for the **HAMMERHEAD.jl** ecosystem. Currently this package implements:
- [x] Cartesian
- [x] Keplerian
- [x] Delaunay
- [x] Modified Equinoctial
- [x] Spherical
- [x] Cylindrical
- [x] Unified State Model
- [x] Milankovich
- [] EDROMO
- [] Kustaanheimo-Stiefel
- [] Stiefel-Scheifel

This package may eventually support Attitude Coordinates as well.

## Installation

```julia
julia> using Pkg
julia> Pkg.add("AstroCoords")
```