# Changelog

All notable changes to Occultations.jl are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the package uses
[semantic versioning](https://semver.org/).

## [0.1.0] – unreleased

First version.

### Fixed

- Direct light-curve fits: an event shorter than the sampling, lying between two exposure
  midpoints, was missed (model flux 1); a light curve wholly inside the occultation gave NaN.
- `profile` of the only free parameter returned NaN, and `confidence_interval` then reported a
  vanishing interval; profiled values outside the parameter bounds got a finite χ². Now:
  `Inf` outside the bounds, `NaN` for a failed refit, intervals limited by the bounds, `±Inf`
  where Δχ² does not reach the threshold.
- Covariance from the singular value decomposition of the Jacobian instead of the
  pseudoinverse of JᵀJ: directions the data do not constrain get σ = `Inf` instead of 0.

### Added

- **Observation reports** in the CORA format: `read_cora_report`, `cora_chord`.
- **Light curves**: readers (text, PyMovie, Tangra, AstroImageJ), normalisation, event
  detection; a physical model with Fresnel diffraction over the passband, the stellar disc,
  the exposure time, finite opacity and several drops per light curve; fitted timings.
- **Double stars**: `double_star=true` in joint light-curve fits: two point sources with a free
  flux ratio and separation, each occulted by the limb; separation, position angle and
  magnitude difference derived. Example: 2007 OC10, 2022 (Gómez-Limón et al. 2025).
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
  only with `extend=true`. The fundamental-plane plots have a compass (E, N; east to the left, as
  on the sky); the dimensions of apparent ellipses are labelled a′, b′.
- **Examples**: re-reductions of 23 published studies (20 objects, 44 events), synthetic
  examples (binaries, contact binaries, rings; stations at random offsets, on land and in the
  dark, with their own equipment) and a tutorial, with their data and outputs;
  the comparison with the publications and with SORA in the manual. The published events are
  grouped by the kind of analysis (reproduced, with differences), as in the paper.
- **Quick-start Jupyter notebook** (`examples/quickstart.ipynb`): prediction, chords, limb fit,
  astrometry, posterior and maps for the Varda occultation of 2018, offline.
- **Comparison with SORA on the same data** (`examples/sora/`, outputs in
  `examples/output/sora/`): nine published events reduced with both packages (published –
  Occultations.jl – SORA), and the computing time of both.
