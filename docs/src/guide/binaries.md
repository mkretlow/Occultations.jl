# Binaries and contact binaries

Three limb models describe bodies that are not a single convex outline:

| Model | Shape | Typical use |
|---|---|---|
| [`BinaryLimb`](@ref)`(primary, secondary)` | two separate bodies, each with its own limb model | a body and its satellite; a binary system of similar bodies |
| [`ContactBinaryLimb`](@ref)`()` | two overlapping ellipses along a common axis (2-D) | a bilobed body seen in one event |
| [`ContactBinaryEllipsoids`](@ref)`()` | two triaxial ellipsoids on the body x-axis, rotating (3-D) | several events at different rotational phases, light-curve amplitudes |

Their apparent limb is a [`CompositeOutline`](@ref): the union of the outlines of the
components, each with its own centre. The union is in general not star-shaped. A chord can
cross it several times, and the boundary is not a radius per direction from a centre.
Everything else (negative chords, clock offsets, priors, posteriors, multi-event fits) works
as for single bodies.

## Observations: chords with several drops

A chord that crosses two bodies, or both lobes of a contact binary with a gap between them,
has several drops. Give its contacts as a vector:

```julia
Chord("Station 4", obs; contacts=[Contact(:immersion, t1, 0.05), Contact(:emersion, t2, 0.05),
                                  Contact(:immersion, t3, 0.05), Contact(:emersion, t4, 0.05)])
```

Immersions and emersions must alternate. [`occulted_intervals`](@ref)`(chord)` lists the drops,
and `plot(event)` draws each of them.

**Contact labels.** A contact can be assigned to a component:

* `Contact(...; body=:primary)` or `body=:secondary` for a [`BinaryLimb`](@ref);
* `body=:lobe1` or `body=:lobe2` for the contact-binary models.

A labelled contact is compared with the crossing of that component alone. An unlabelled
contact (`body=:auto`, the default) is matched to the nearest crossing of the same kind
(immersion or emersion) of the whole outline. Labels make the fit robust when a satellite is
crossed by only one or two chords, or when the assignment is otherwise ambiguous.

**Light curves with several drops.** [`detect_occultations`](@ref)`(lc; method=:threshold)` finds the
drops as runs of points below half the baseline; the default box search would find two deep
drops with a short gap between them as one. Time each one in its own window with
[`fit_timing`](@ref)`(lc; tmin, tmax, immersion, emersion, …)`, between the neighbouring drops,
as in the [contact-binary example](../generated/sim_contact_binary.md).

[`simulate_lightcurve`](@ref) and [`occultation_flux`](@ref) take a list of drops
`[(t_imm, t_emm), …]`:

* drops more than about ten Fresnel scales apart are independent;
* closer drops are computed together. The complex amplitudes of the opaque strips add
  (Babinet's principle), so the light diffracted into a narrow gap between two lobes is part
  of the model.

The joint fit of the light curves (`photometry=true`) uses the same model. It is exact where
the timing of each drop on its own is not: the timing treats the gap as fully lit.

## The visible segments of the chords

Where a chord saw the star (before the first drop, between drops, after the last drop), no
part of the body can have been. For the composite models these **visible segments** are
constraints, just like negative chords. Each segment is shortened at every contact by 3σ of
its timing (at least the chord's `sensitivity`).

| Value | Behaviour |
|---|---|
| `visible=:soft` | penalty for the depth by which the limb enters a segment, in units of `sensitivity × speed` |
| `:hard` | no limb may enter a segment |
| `:none` | the segments are ignored |
| `:auto` (the default) | `:soft` for the composite models; `:none` for single convex limbs, which cannot produce a drop where none was seen, so the published results are unchanged |

The visible segments decide between solutions that the contacts alone allow:

* a chord with two drops proves the gap between the lobes;
* a satellite cannot grow across a chord that saw the star.

## A body with a satellite, a binary system

```julia
r = fit_limb(ev; model=BinaryLimb())                              # ellipse + circle (default)
r = fit_limb(ev; model=BinaryLimb(EllipseLimb(), EllipseLimb()))  # an ellipse for the satellite too
```

**Parameters.**

* The primary's parameters keep their names; `center_f` and `center_g` are the centre of the
  primary.
* The secondary's parameters take the prefix `sec_`: `sec_radius`, or
  `sec_equatorial_radius`, `sec_oblateness` and `sec_position_angle`.
* `sep_f` and `sep_g` [km] are the position of the secondary relative to the primary in each
  event (suffix `_k` in multi-event fits).

**Starting values.** The contacts of each event are split into two groups (by their labels,
otherwise by clustering), and each group gives the starting values of one body.

**Circle or ellipse for the satellite.** The circle (three parameters) is the default: few
chords usually cross a satellite. An ellipse has five parameters. It needs at least three
chords on the satellite, better four or more, and the visible segments to exclude elongated
solutions. With two chords it is undetermined: the fit reports a singular covariance. Compare
the two with the Bayesian information criterion, χ² + k ln n.

**Derived quantities.**

* those of each component;
* `size_ratio`: the ratio of the area-equivalent radii;
* `mass_ratio`: `size_ratio³`, for equal densities;
* `system_equivalent_radius`;
* per event: `separation` [km], `separation_mas` and `separation_pa` (position angle of the
  secondary from the primary).

**Astrometry.** [`astrometry`](@ref)`(fit; reference=…)` gives one of three positions:

| `reference` | Position |
|---|---|
| `:center` (default) | the centre of the primary |
| `:barycenter` | equal densities |
| `:photocenter` | equal albedos |

For a satellite of a quarter of the primary's size, the three differ by several
milliarcseconds at TNO distances. For an orbit fit of the system, use the barycentre.

The separation in each event is free. A Keplerian orbit of the secondary (and with it the
mass of the system) over several events is a natural extension.

## Contact binaries

**In 2-D**, [`ContactBinaryLimb`](@ref) has these parameters:

* `a1`, `b1`, `a2`, `b2`: the lobes, with `a` along the common axis;
* `axis_pa`: the position angle of the axis, from lobe 1 to lobe 2;
* `neck`: the centres are `d = (a1 + a2)(1 − neck)` apart.

| `neck` | Lobes |
|---|---|
| 0 | the lobes touch |
| > 0 | the lobes overlap, with a waist |
| < 0 | there is a gap (down to −0.5) |

The centre is the area-weighted centre of the lobes. Lobe 1 is always the larger one; the
fit swaps the lobes and turns the axis by 180° if needed. The derived quantities are:

* `equivalent_radius` (of the area of the union);
* `lobe1_equivalent_radius` and `lobe2_equivalent_radius`;
* `size_ratio` and `mass_ratio`;
* `separation`;
* `length` (along the axis);
* `neck_width` (the width of the waist).

The starting values come from the ellipse fitted to all contacts, whose major axis gives the
axis of the lobes.

**In 3-D**, [`ContactBinaryEllipsoids`](@ref) has two triaxial ellipsoids (`a1`, `b1`, `c1`,
`a2`, `b2`, `c2`) with a `neck` as above. They rotate with `pole_ra`, `pole_dec`, `W0` and
`period`, as in [`EllipsoidLimb`](@ref); `center_f` and `center_g` are the centre of mass. The
shape is a [`BilobedEllipsoid`](@ref).

**Exact silhouettes.** The silhouette of two ellipsoids is exactly the union of their two
projected ellipses, so each event constrains the shape without approximation. Rotational
amplitudes and phases come from the area of the union over a rotation.

**Volume.** The volume (and `volume_equivalent_radius`) counts the overlap of the lobes once.

**Output.** `ellipsoid_mesh(shape)` writes the fitted body as an OBJ file, and
[`occview`](@ref) shows it in 3-D.

## Fitting strategy

Composite fits can have several minima. Which contacts belong to which body, mirror solutions
and the orientation of the axis are all ambiguous. These help:

* labels on the contacts where the assignment is clear;
* the default multistarts (`starts=8`); for 3-D models, a grid of starting poles (see the
  [3-D example](../generated/sim_contact_binary_3d.md));
* hard negative chords and the visible segments, which exclude many minima;
* distinct solutions are best compared side by side by their χ². The ensemble sampler does
  not move between separated modes.

## Plots and maps

`plot(fit)` and `plot(posterior)` draw:

* the union of the outlines;
* a cross at the reference centre, and small crosses at the centres of the components;
* each drop of the chords.

The plot box covers the whole system.

**Uncertainty band.** The bands of single bodies are built in polar coordinates around one
centre, which fails for a non-star-shaped limb. For composite outlines the band is the
**occupancy** of the drawn limbs: the region covered by between 15.9 % and 84.1 % of the limbs
drawn from the covariance or the posterior (for `nsigma=1`). It is labelled "occupancy" in the
legend.

**Maps.** On maps ([`occ_map`](@ref), [`geomap`](@ref)) the shadow path of a composite body
spans its whole extent across the track. For a separated satellite this includes the gap
between the bodies.

## Examples

* [Huya and its satellite, 2023 June 24](../generated/huya_2023_satellite.md): the published
  timings of [Rommel2025](@citet); Huya and the satellite (one chord) fitted together, the two
  solutions on either side of the satellite chord, the lower limit of the satellite's size,
  the astrometry of the satellite and the barycentre.
* [Synthetic: a satellite and a binary system](../generated/sim_binary.md): the satellite as a
  circle and as an ellipse, two chords on the satellite, astrometry of the barycentre, a binary
  of two similar bodies with two drops per chord, and the coverage of the uncertainties over
  100 simulations.
* [Synthetic: a contact binary in one event](../generated/sim_contact_binary.md): light curves
  with two drops, their timing, `ContactBinaryLimb`, its posterior, and a joint fit of the
  light curves with the diffraction of the narrow gaps.
* [Synthetic: a contact binary in 3-D](../generated/sim_contact_binary_3d.md): four events and
  a light-curve amplitude, `ContactBinaryEllipsoids`, a known and an unknown pole.
