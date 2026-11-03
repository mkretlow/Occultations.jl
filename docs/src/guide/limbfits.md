# Limb fits

## Measurement model

* **Contacts**: by default the residual of a contact is the distance *along the chord*
  between the observed and the predicted contact (timing residual × speed), with
  uncertainty `√((σ_t |v|)² + σ_limb²)` (`σ_limb`: an allowance for limb topography and
  unmodelled timing errors; `fit_limb(...; σ_limb=:auto)` sets it so that χ²/dof = 1). This is
  the measurement model: a grazing chord, whose timing hardly constrains the radius, gets
  its proper weight, and the residual stays finite at tangency; a chord that misses the
  limb is penalised by the distance by which it misses. `residual=:radial` gives SORA's
  radial residuals; the difference is explained in
  [Residuals along the chords or radial residuals?](../theory/residuals.md).
* **Negative chords**: `negative=:hard` (default) forbids a limb crossing the observed
  track (zero posterior probability); `:soft` penalises the depth of the crossing, in units
  of `sensitivity × speed`; `:none` ignores negative chords.
* **Light curves** (`photometry=true`): the flux of every light curve is modelled directly
  from the limb (predicted contacts and velocities normal to the limb, diffraction, stellar
  disc, passband, exposure) or from the `atmosphere`; the flux levels are profiled. This
  needs no separate timing step and handles grazes, partial and negative light curves
  uniformly.
* **Clock offsets**: `fit_dt=["station", …]` fits `dt_station` (true = recorded + dt).
* **Priors**: `priors=Dict(:oblateness => (0.1, 0.02))`, `fix=[...]`, `init=Dict(...)`,
  `bounds=Dict(...)`, or `setparam!` on a problem.

## Limb models

| Model | Parameters | Use |
|---|---|---|
| `EllipseLimb()` (default) | `equatorial_radius` *a*, `oblateness` ε = 1 − b/a, `position_angle` (of the minor axis = apparent pole, N→E) | one event; or several events with `per_event=(:position_angle,)` etc. |
| `CircleLimb()` | `radius` | few chords, atmospheres |
| `EllipsoidLimb()` | `a`, `b`, `c` (body x, y, z), `pole_ra`, `pole_dec`, `W0` (prime meridian at the reference epoch), `period` (fixed by default) | several events at different rotation phases/aspects; `RotationalAmplitude` and `RotationalPhase` constraints; `principal=true` for rotation about the shortest axis |
| `MeshLimb(mesh)` | `scale`, `rotation` (sky), optionally `pole_ra`, `pole_dec`, `W0` (`fit_orientation=true`) | shape models (`read_obj`), `ellipsoid_mesh(a, b, c)` |
| `BinaryLimb(primary, secondary)` | those of both models (`sec_` prefix for the secondary), `sep_f`, `sep_g` | a body with a satellite, a binary system ([Binaries and contact binaries](binaries.md)) |
| `ContactBinaryLimb()` | `a1`, `b1`, `a2`, `b2`, `axis_pa`, `neck` | a contact binary in one event |
| `ContactBinaryEllipsoids()` | `a1`, `b1`, `c1`, `a2`, `b2`, `c2`, `neck`, pole, `W0`, `period` | a contact binary over several events |
| `RingedLimb(body; rings=n)` | those of the body, `ring_pole_ra`, `ring_pole_dec`, `ringk_radius`, `ringk_width`, `ringk_tau` | a body with rings ([Ring systems](rings.md)) |

Every event has its own centre (`center_f_k`, `center_g_k`: offset of the body from its
ephemeris in the fundamental plane), so multi-event fits are not biased by ephemeris
errors. A new model only has to declare its parameters and return an `Outline` per event.

```julia
# three events of a rotating triaxial body plus a light-curve amplitude
amp = RotationalAmplitude("2023-06-01", 0.31, 0.02)          # Δm ± σ
phs = RotationalPhase(ev1.tref, 0.63, 0.05)                   # rotations since maximum light ± σ
r3 = fit(LimbProblem([ev1, ev2, ev3]; model=EllipsoidLimb(principal=true), amplitudes=[amp],
                     phases=[phs], priors=Dict(:pole_dec => isotropic_pole())))
r3[:a], r3[:b], r3[:c], r3[:volume_equivalent_radius], r3[:apparent_a_2], r3[:aspect_2]

# a shape model
mesh = read_obj("shape.obj"; scale=1.0)                       # km
rm = fit_limb([ev1, ev2]; model=MeshLimb(mesh; fit_orientation=true), fix=[:pole_ra, :pole_dec])
```

The classical formulas of [Dunham2019](@citet) and [Magnusson1986](@citet) for the apparent axes,
and the rotational amplitude of a triaxial ellipsoid, are available
(`projected_axes_dunham`, `projected_axes_magnusson`, `rotational_amplitude`); the package
uses the exact projection (`project`), which agrees with both to machine precision.

## Derivatives

Fits use **automatic differentiation** (ForwardDiff) for their Jacobians:
* where it applies: limb fits from contacts and negative chords (ellipse, circle, 3-D
  ellipsoid) and timing fits;
* it gives exact derivatives, and is faster (≈ 3× for an ellipse fit);
* joint light-curve, atmosphere and mesh fits fall back to central finite differences
  automatically;
* `fit(problem; ad=false)` forces finite differences.

## Limb fit or posterior?

Both use the same likelihood: contacts, negative chords, light
curves and constraints. They differ in how they describe the uncertainty.

* **Limb fit** (`fit_limb`, least squares):
  * it finds the best-fitting parameters and approximates the uncertainty by a Gaussian
    around them, from the covariance (the curvature of χ² at the minimum);
  * this is exact for a linear model with Gaussian errors, and fast;
  * it fails where χ² is not a paraboloid: few chords, a limb constrained on one side only,
    parameters at a bound, hard constraints (negative chords, a ≥ b ≥ c), or degeneracies
    such as a mirror solution or the 3-D valley;
  * the 68.3 % region of the fit plot is drawn from this Gaussian, with the draws that violate a hard
    constraint rejected (a truncated Gaussian).
* **Posterior** (`sample`, ensemble MCMC):
  * it explores the actual probability distribution of the parameters, including the
    priors (e.g. an isotropic pole, clock-offset priors);
  * it is slower (thousands of model evaluations);
  * it is the reliable description whenever the two disagree: asymmetric errors, banana-
    shaped correlations, multiple solutions.
* The plots differ in their black curve: the best fit, or the maximum-posterior sample.
  With flat priors these coincide up to sampling noise.

If the regions and the errors of both plots agree, the Gaussian approximation is adequate.

## Negative chords

Negative chords are hard constraints by default (`negative = :hard`): a limb that crosses
the track of a negative chord is excluded from the fit, the 68.3 % region and the posterior. The
option `negative = :soft` accepts a crossing with a penalty, scaled by
`sensitivity × speed`. It suits a chord whose detection limit is uncertain, where a short
event within one long exposure or in noisy data could have been missed. The region can still
touch a negative chord where many admissible limbs pass close to it.
