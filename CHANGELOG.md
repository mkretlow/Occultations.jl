# Changelog

All notable changes to Occultations.jl are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the package uses
[semantic versioning](https://semver.org/).

## [0.1.0] – unreleased

First version.

### Added

- **Light curves**: readers (text, PyMovie, Tangra, AstroImageJ), normalisation, event
  detection; a physical model with Fresnel diffraction over the passband, the stellar disc,
  the exposure time, finite opacity and several drops per light curve; fitted timings.
- **Atmospheres**: hydrostatic, spherically symmetric atmospheres (N₂, CH₄, CO) ray-traced to
  the refracted flux; the surface pressure as a fit parameter.
- **Limb fits** for one or several events: circle, ellipse, triaxial ellipsoid with pole and
  rotation, triangle meshes; binaries and satellites (`BinaryLimb`), contact binaries
  (`ContactBinaryLimb`, `ContactBinaryEllipsoids`), ring systems (`RingedLimb`, ring
  detections as contacts or mid-times); negative chords as hard or soft constraints, clock
  offsets, rotational amplitudes and phases, priors, limb and ring allowances (`σ_limb`,
  `σ_ring`); joint fits of the light curves with the limb.
- **Uncertainties**: covariance, propagation to derived quantities, χ² profiles and
  confidence intervals, Bayesian posterior sampling.
- **Astrometry** of the body from the fitted centre, including the star's position error.
- **Geometry and predictions**: ERFA-based reduction, ground and space observers, JPL
  Horizons, SPICE, Gaia, SBDB, Keplerian orbits, closest approach and shadow paths.
- **Plots and maps**: Plots.jl recipes (light curves, fits, chords, posterior limbs, corner
  plots), 3-D views (`meshview`, `occview`, `occanimate`), maps and globes with GeoMakie.
  Limb plots draw the occulted parts of the positive chords; their tracks beyond the limb
  only with `extend=true`.
- **Examples**: re-reductions of 22 published studies (19 objects, 43 events), synthetic
  examples (binaries, contact binaries, rings; stations at random offsets, on land and in the
  dark, with their own equipment) and a tutorial, with their data and outputs;
  the comparison with the publications and with SORA in the manual.
- **Comparison with SORA on the same data** (`examples/sora/`, outputs in
  `examples/output/sora/`): nine published events reduced with both packages (published –
  Occultations.jl – SORA), and the computing time of both.
