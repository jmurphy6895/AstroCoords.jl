# Usage

The main way to convert between the elements is by using the AstroCoords struct. After one has been instantiated simply pass it and a gravitational parameter to a constructor of the desired element set and the package will handle the rest.

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
kep_state = Keplerian(cart_state, μ)
```

While not explicitly export if the user desired to avoid the structs, simply find the appropriate conversion inside of the coordinate_changes.jl file. Note, it make take multiple conversions to get to the desired set when using this approach.

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

kep_state = AstroCoords.cart2koe(state, μ)
```