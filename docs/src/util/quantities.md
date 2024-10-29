# Quantities

This package implments a number of interfaces to get common astrodynamics quantities regardless of coordinate definition used. All of these have a common interface where the function simply takes the AstroCoords struct and a gravitational parameter.

The quantities currently supported are:
* meanMotion
* orbitalPeriod
* orbitalNRG
* angularMomentumVector
* angularMomentumQuantity

A sample function call is provided below

```julia
state = [
    -1076.225324679696
    -6765.896364327722
    -332.3087833503755
    9.356857417032581
    -3.3123476319597557
    -1.1880157328553503
]

μ = 3.986004415e5

cart_state = Cartesian(state)

NRG = orbitalNRG(cart_state, μ)
```