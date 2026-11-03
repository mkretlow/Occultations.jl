# Ring systems

[`RingedLimb`](@ref)`(body; rings=n)` fits a body together with `n` circular rings. The rings
lie in the plane of the ring pole and are centred on the body.

| Kind | Parameters |
|---|---|
| body | any [`LimbModel`](@ref) (`EllipseLimb`, `CircleLimb`, `EllipsoidLimb`, …), its parameters unchanged |
| ring plane | the ring pole `ring_pole_ra`, `ring_pole_dec` [deg] |
| ring `k` | mid radius `ringk_radius` [km], radial width `ringk_width` [km], normal optical depth `ringk_tau` |

## Geometry

On the sky a circular ring is an **elliptical annulus**, a [`ProjectedRing`](@ref). Two angles
describe it in each event:

* the **opening angle** B: the angle between the line of sight and the ring plane;
* the **position angle** P of the projected ring pole.

Its edges are ellipses with semi-axes `r ± W/2` and `(r ± W/2)|sin B|`, with the minor axis
along P. The ring is partly transparent. Its **apparent opacity** is

```
p′ = 1 − exp(−τ_N / |sin B|),
```

the fraction of the stellar flux that it blocks.

[`ring_geometry`](@ref)`(pole_ra, pole_dec, u)` gives B and P for a line of sight `u`.
[`project_ring`](@ref) gives the projected ring.

**Ambiguities of the ring plane:**

* A pole and its antipode describe the same ring plane. The fitted pole is the one in the
  northern hemisphere.
* **One event** fixes only the ellipse: the radius, |B| and P. The plane and its mirror image
  (opening angle −B, the pole tilted towards or away from the observer) give the same ellipse.
* **Several events** at different geometries decide between the two planes, but only if the
  line of sight changes enough. For distant bodies it changes by about a degree per year. See
  the [three-event example](../generated/sim_ring_pole.md).

## Observations

A ring dip of a chord is a pair of contacts labelled with its ring: the passage of the
star through the annulus.

```julia
Chord("Station 3", obs; contacts=[
    Contact(:immersion, t1, 0.05; body=:ring1), Contact(:emersion, t2, 0.05; body=:ring1),   # ring, ingress
    Contact(:immersion, t3, 0.03), Contact(:emersion, t4, 0.03),                             # the body
    Contact(:immersion, t5, 0.05; body=:ring1), Contact(:emersion, t6, 0.05; body=:ring1)])  # ring, egress
```

* **Ring contacts must be labelled** (`:ring1` … `:ring9`; `:ring` is `:ring1`). Each is compared
  with the nearest crossing of the same kind of that ring's annulus, along the chord.
* **Unlabelled contacts belong to the body.**
* Negative chords and the visible segments of the positive chords constrain the **body only**:
  a ring dip may escape detection.

**Ring dips in light curves.** A ring dip is short and partly transparent.

* **Finding them:** the box search `detect_occultations(lc; snr_limit=6, maxdur=3)` outside
  the body drop.
* **Identifying them:** by their symmetry about the body and by comparison with the
  neighbouring chords. A faint ring is seen only in the better light curves, so a blind
  assignment by order is error-prone.
* **Timing them:** with the opacity free, `fit_timing(lc; tmin, tmax, fit_opacity=true)`, in a
  window between the neighbouring drops.

The [two-ring example](../generated/sim_rings.md) shows the whole chain, from simulated light
curves to the fits.

**Which data constrain what:**

* The **timings** constrain the radii and widths of the rings.
* The **optical depths** come from the apparent opacities. They are fitted in the joint fit
  of the light curves (`photometry = true`). There the rings are partly transparent strips
  with diffraction: the complex amplitudes of the strips add, each weighted by its
  transmission, as for the gaps of a contact binary. The stellar disc and the exposure are
  included too.
* Fix the optical depths (`fix=[:ring1_tau, …]`) in a fit of the timings alone.
* **Where the exposures are comparable to the duration of a ring crossing**, the width and
  the optical depth of a ring are correlated. Their product, the equivalent depth W τ_N, is
  then the robust quantity.
* **A faint ring is found only in the best light curves.** A ring without any detected dip has
  no data; leave it out of the model.

**Simulation.** [`simulate_lightcurve`](@ref)`(t, drops; opacity=[…])` simulates light curves
with ring dips, and `predicted_contacts(event, observer, ring)` predicts ring crossings.

## Derived quantities

| Name | Meaning |
|---|---|
| `ring_opening`, `ring_pa` | B and P (suffix `_k` per event in multi-event fits) |
| `ringk_equivalent_depth` | W τ_N [km] |
| `ringk_opacity` | the apparent opacity p′, per event |
| `ringk_equivalent_width` | W p′ [km], per event |

## Plots

`plot(fit)` and `plot(posterior)` draw:

* the edges of the rings, in orange;
* the body;
* the chords, with each ring dip and each body drop as a segment.

The plot box covers the rings.

## Not yet included

* eccentric or inclined rings, and rings with their own centre;
* radial profiles within a ring;
* maps of the ring shadow paths.

## Examples

* [Synthetic: a body with two narrow rings](../generated/sim_rings.md), Chariklo-like:
  * light curves with ring dips;
  * detection and timing of the dips;
  * the fit of the timings, then the joint fit of the light curves for the optical depths.
* [Synthetic: a distant ring seen in three events](../generated/sim_ring_pole.md), Quaoar-like:
  * the ring ellipse of each event and its mirror ambiguity;
  * the ring pole from the three events together.
