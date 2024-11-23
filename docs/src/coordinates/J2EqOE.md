# J2 Perturbed Modified Equinoctial

The J2 equinoctial orbital elements (J2EqOE) are a specialized set of six non-singular orbital elements designed to account for the dominant J2 perturbation caused by Earth's oblateness. This set is particularly suited for uncertainty propagation in satellite dynamics, as it absorbs significant nonlinearities in the equations of motion.

![image](../assets/equinoctial-frame.png)
*Modified Equinoctial Orbital Element Frame [2]*

## Components

This set is very similar to the Modified Equinoctial Orbital element set except this set uses mean motion instead of the semi-major axis. This set accounts for the J2=perturbation via the Brouwer-Lyddane theory.

* **Orbital Shape and Size**
    * mean motion (n): The mean motion of the orbit.
    * first component of eccentricity (h): First of two parameters that replace the scalar eccentricity. f=e*cos(ω + Ω) These components prevent singularities when the eccentricity is zero
    * second component of eccentricity (k): Second of two parameters that replace the scalar eccentricity. g=e*sin(ω + Ω)

* **Orbital Orientation**
    * first component of inclination (p): First of two parameters that replace the scalar inclination. h=tan(i/2)*cos(Ω) These components ensures the elements remain well-defined for equatorial orbits
    * second component of inclination (q): Second of two parameters that replace the scalar inclination. h=tan(i/2)*sin(Ω)

* **Position within the Orbit**
    * longitude of periapsis (L): This replaces true anomaly and combines information about the satellite’s position within the orbit. L=ω + Ω + f

## References
[1]: https://link.springer.com/article/10.1007/s10569-021-10004-0
[2]: https://degenerateconic.com/modified-equinoctial-elements.html