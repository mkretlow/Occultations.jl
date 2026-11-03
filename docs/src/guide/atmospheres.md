# Atmospheres

```julia
atm = xv93_atmosphere(:CH4)       # 47 K, +2.7 K/km up to 20 km, free molecular flow above 25 km
atm = pluto_atmosphere()          # or Atmosphere(:N2; T0=40, dTdz=3, z_inv=22, z_fmf=Inf, width=2)
prof = atmosphere_profile(atm, 124.0, 244.0, body_gm(body), dist_km; λ=0.6, star_diam=0.09)
plot(prof); half_light_radius(prof)

# joint fit of the light curves: radius, centre and surface pressure
r0 = fit_limb(ev; model=CircleLimb(), photometry=true)                      # airless first
ra = fit_limb(ev; model=CircleLimb(), photometry=true, atmosphere=atm,
              init=Dict(:radius => r0[:radius][1], :center_f => r0[:center_f][1],
                        :center_g => r0[:center_g][1], :pressure => 50.0), starts=6)
ra[:pressure]                                                               # nbar ± σ
r0.chi2 - ra.chi2                                                           # significance of the atmosphere

# test for an atmosphere: profile likelihood in the pressure → detection or upper limit
t = atmosphere_test(ev; atmosphere=atm, init=Dict(:pressure => 50.0, ...), airless_init=Dict(:radius => r0[:radius][1]))
t.significance, t.interval, t.upper_limit                                   # σ, (lower, upper) nbar at 95 %
atmosphere_test(ev; atmosphere=atm, confidence=0.99865)                     # 3σ limits
```

The bending angle is computed by integrating the refractivity gradient along the ray
(`ω(r) = 4r ∫ ν′(r+u²)/√(2r+u²) du`); the free-molecular-flow tail `n ∝ r⁻²`, which
dominates the bending of tenuous extended atmospheres, is integrated analytically to
infinity. The top of the inversion layer and the onset of the free molecular flow are
smooth transitions (`width`, default 2 km): a kink in the density gradient produces a
spurious jump of the refracted flux (φ ∝ 1/(1 + D dω/dr)). The diffraction by the surface
edge multiplies the refracted flux (an approximation valid when the refractive scale is
much larger than the Fresnel scale). Start atmospheric fits from an airless photometric
fit: the χ² surface has local minima far from the solution.

## Reading the plots

### Light curves with an atmosphere

Each panel shows the data with
the atmospheric (red) and airless (dashed) models.
* The airless model drops abruptly at the surface, with Fresnel fringes.
* Refraction makes the drop gradual (bending of the star light), starting tens of km
  above the surface.
* A near miss, such as Fukushima for 2002 XV93, can show a dip without a disappearance.

### Atmosphere test

`atmosphere_test(...)` gives a profile-likelihood curve:

* **Axes:** x is the surface pressure (logarithmic; p = 0 is plotted at the left edge). y is
  Δχ² = χ²(p) − χ²_min. For each pressure, the limb (radius, centre) is refitted with the
  pressure fixed.
* **Minimum:** the best-fitting pressure.
* **Horizontal lines:** Δχ² = 2.71 and 9. The pressures where the curve is below a line are
  allowed by the data. Where it crosses the 2.71 line gives one-sided 95 % limits; where it
  crosses the 9 line gives 3σ limits.
* **Detection:** compare the left end (an airless body) with the minimum. The significance
  is √Δχ²(0) σ (e.g. 2002 XV93: Δχ² = 64.9, 8.1σ).
* **No atmosphere:**
  * The minimum is at or near p = 0, or the curve is flat on the left. Below a certain
    pressure, the refraction is too weak to change the light curve.
  * Where the curve crosses Δχ² = 2.71 on the right is the **95 % upper limit** on the
    surface pressure (`upper_limit`).
  * Steep walls show that higher pressures would bend the light visibly before the surface
    contact.
* **Assumptions:**
  * The result assumes the atmosphere model (gas, temperature profile) and the error bars of
    the light curve. Check χ²/dof ≈ 1 of the best fit.
  * Δχ² scales with 1/σ²: underestimated errors inflate both the significance and the
    exclusion.

### Atmosphere profile

`plot(atmosphere_profile(...))` shows the normalised flux against
distance from the shadow centre. The flux falls gradually where refraction defocuses the star
light, then to zero at the shadow of the surface. Comparing gases or pressures shows how
distinguishable they are.

![2002 XV93 light curves](../assets/figures/2002xv93_lightcurves.png)

![2002 XV93 atmosphere test](../assets/figures/2002xv93_atmosphere_test.png)
