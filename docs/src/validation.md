# Validation

`julia --project -e 'using Pkg; Pkg.test()'` runs the test suite (373 tests, ≈ 8 min). Among the tests:

* Tracks reproduce the full reduction to < 2 m; the analytic and numerical limb crossings
  agree to 1 ms.
* The exact projection of an ellipsoid equals the [Dunham2019](@citet) formulas, and those of
  [Magnusson1986](@citet) with the phase shifted by 90°, to 10⁻⁹; the Haumea example (1161 × 852 ×
  513 km, aspect 76.3°) gives 852 × 569 km; mesh silhouettes agree to 0.5 % in area; the
  closed-form rotational amplitude equals the numerical one.
* Synthetic events: timings, ellipse parameters, clock offsets (0.7 s), 3-D axes (340, 290,
  250 km) and rotation phase from three events, mesh scale and phase, and the surface
  pressure of a 2002 XV93-like atmosphere (124 nbar) are recovered within their errors;
  the posterior widths agree with the covariance and with χ² profiles.
* Ray tracing: the bending angle of an isothermal atmosphere agrees with
  `ω ≈ −ν √(2πr/H)` to 2 % (the expected O(H/r) correction); the free-molecular-flow limit
  `ω = −πν` is exact.
* Real events: the Huya occultation of 2019 March 18 reproduces the published ellipse
  [SantosSanz2022](@citep) within its uncertainties; the five 2002 KX14 events
  [Rizos2025](@citep) reproduce the published geometry, semi-major axis and diameter within 1–2σ; the
  2003 UY117 occultation [Kretlow2024](@citep) reproduces the published limb, χ²/dof and
  astrometry. For 2002 XV93 [Arimatsu2026](@citep), the ray-traced atmosphere reproduces the
  authors' model light curves to 0.01 rms, and the joint light-curve fit recovers their
  surface pressures for CH₄, N₂ and CO within about half their errors.
* Generic engine: χ² profiles and confidence intervals at parameter bounds and for flat
  profiles; infinite σ for directions the data do not constrain; light curves with an event
  between two exposures or wholly inside it; negative chords crossed between their samples;
  sampler starts given explicitly.

The live services (MPC codes, Horizons, SBDB; `using HTTP`) were checked by hand and are not
part of the automated tests. The real-event atmosphere results come from the 2002 XV93
example; the atmosphere tests use synthetic data.

Synthetic-data validation is no substitute for real events: compare with published
reductions (e.g. SORA, Occult) before relying on new results.

## Not implemented (yet)

Eccentric or inclined rings and radial ring profiles; central flashes beyond the
spherical, transparent atmosphere model (haze, oblate atmospheres); limb topography models;
light-curve readers for other programs (LiMovie, Occular); prediction services
(occultation catalogue queries, Occult Watcher export); automatic differentiation of the
atmosphere fits (these use finite differences).
