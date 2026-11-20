# Changelog

All notable changes to Occultations.jl are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the package uses
[semantic versioning](https://semver.org/).

## [Unreleased]

### Added

- **Atmospheres**: any temperature profile `T(r)` (`Atmosphere(...; temperature=f)`), the
  template of Dias-Oliveira et al. (2015) (`DO15Temperature`), the pressure at a reference
  radius (`rref`, e.g. Triton's p₁₄₀₀), Triton's 2017 model of Marques Oliveira et al. (2022)
  (`triton_atmosphere`), the surface pressure of a profile (`surface_pressure`); the stellar
  disc is averaged in two dimensions near the central flash; `background=:free` in joint
  light-curve fits profiles the zero level of each light curve (the flux of the occulting body,
  e.g. Triton). Tests against the Baum & Code (1953) solution and the half-light pressure of
  Sicardy (2022).
- Manual page comparing the atmosphere model with Baum & Code, Sicardy and Arimatsu et al.
- Data of the Triton occultation of 2017 October 5 (stations, star, ephemeris; light curves
  pending).
- The 3-D guide page points to the interactive demo `examples/occview_demo.jl`.
- `contact_gaps(timing_fit)`: whether each fitted contact falls in a gap between exposures, and
  the gap bounds (the honest interval there; the covariance error is not meaningful).

### Changed

- `fit`: for problems with constraints (negative chords and visible segments of limb fits,
  not light-curve fits) every start is also solved by a staged search (`relax=true`, the
  default): the constraints as weak penalties first, tightened in steps, then the problem itself
  from the nearest admissible point; the better solution is kept. Since negative chords are
  checked exactly, hard constraints can split the admissible parameters of composite limbs into
  pieces that a local search from a random start cannot leave (the synthetic contact binary:
  χ²/dof 3316 → 1.8 in its first pass). `relax=false` restores the previous search.
- Example `sim_contact_binary`: soft negative chords, a two-step contact-binary fit (without,
  then with the visible segments), the second pass started from the first.
- 2002 GZ32 moved to the reproduced examples (one event); its comparison updated (radial
  residuals and residuals along the chords now give the same limb).
- `fit_timing` / `TimingProblem`: the starting values are also searched on a joint (t₁, t₂) grid
  (±4 exposure times), solutions whose event dims no exposure are excluded from the start, and a
  warning is issued if the fit ends in such a degenerate solution. Sparse light curves (one point
  in the event, long dead times) could end with both edges in one gap between exposures.
- `LimbProblem` with `photometry`: if an event has fewer than three timed contacts and `init`
  gives no centre, light curves without timed contacts give starting values from their main drop
  (`detect_occultations`, SNR ≥ 5; not with an atmosphere or a double star), also for `BinaryLimb`
  and `ContactBinaryLimb`; before, such fits started at (0, 0) with the body radius and could
  stop with every chord missing the limb. Fits with contacts or a given centre are unchanged.
- `LimbProblem`: a warning when a light-curve chord without contacts that is not fitted
  photometrically (hence a negative chord) shows a clear drop. `check_negative_chords(fit)`
  skips the light curves of the fit; `check_negative_chords(event, outline; exclude=…)`.
- Docstring of `LimbProblem`: `dt_<chord>` is the correction added to the recorded times (minus
  the clock error).

### Fixed

- `TimingProblem`: the bottom flux level is profiled only while the model has at least three deep
  points (as in `LimbProblem`), not merely because the start had them; a solution that barely
  dimmed any exposure could fit a drop with an absurd bottom level (e.g. −900) and a good χ².
- `confidence_interval`: profile values whose refit fails (an inadmissible start far from the
  best fit, e.g. an emersion fixed before the immersion) are treated as beyond the domain and the
  search backs off towards the best fit, instead of returning `NaN`; the result could depend on
  the state of the random number generator through the replacement of inadmissible starts.
- Atmospheres without a free-molecular-flow tail: the altitude grid could be far too coarse
  for thin layers above an extended atmosphere (2002 XV93-like bodies: flux at the surface
  edge 0.63 instead of ≈ 0.44); the step is now at most 1/50 of the smallest density scale
  height and the rays are spaced evenly in density e-folds. The bending integral resolves the
  deepest scale height also for a high model top.
- Dispersion of the N₂ refractivity: K = 1.091×10⁻²⁹ + 6.282×10⁻³²/λ² m³ (Washburn 1930, as
  used by Dias-Oliveira et al. 2015); the dispersion term was 9 % too large (0.1 % in K).
- Joint light-curve fits with an atmosphere: one ray tracing per event and passband, rescaled
  exactly to each observer's distance (previously one per chord, and the cache was cleared on
  every call with more than 16 light curves).

## [0.1.0] – 2026-11-07

First version.

### Added

- Synthetic diffraction timing example with passband, stellar-disc and exposure smoothing
  comparisons, fitted contacts and residuals; figures exported as PNG and vector PDF.

### Fixed

- Direct light-curve fits: an event shorter than the sampling, lying between two exposure
  midpoints, was missed (model flux 1); a light curve wholly inside the occultation gave NaN.
- `profile` of the only free parameter returned NaN, and `confidence_interval` then reported a
  vanishing interval; profiled values outside the parameter bounds got a finite χ². Now:
  `Inf` outside the bounds, `NaN` for a failed refit, intervals limited by the bounds, `±Inf`
  where Δχ² does not reach the threshold.
- Covariance from the singular value decomposition of the Jacobian instead of the
  pseudoinverse of JᵀJ: directions the data do not constrain get σ = `Inf` instead of 0.
- Negative chords (and the visible segments of positive chords) are checked exactly, not only
  at their 400 samples: a limb could cross a chord between two samples (up to 40 km apart).
  The intervals where a crossing is possible follow from a bound on the rate of change of the
  radial distance and are searched by golden section; also in `check_negative_chords`.

### Added

- `sample(…; init)`: the starting positions of the walkers, as a matrix or a function that
  draws one start (e.g. from broad ranges, to search for other modes).
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
