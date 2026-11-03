# Package design and conventions

## Architecture

The package is organised in layers; each depends only on the ones above it.

| Layer | Contents |
|---|---|
| `core/` | `Epoch` (UTC, two-part JD via ERFA), angles, constants |
| `astro/` | `Star`, ephemerides (`EphemTable`, `EphemKeplerian`, `EphemKernel`, `EphemPlane`), `Body`, observers (`Observer`, `Spacecraft`), fundamental-plane geometry, shadow paths, predictions |
| `body/` | `BodyFrame` (pole, rotation), 3-D shapes (`Ellipsoid`, `TriMesh`), their projections, **`Outline`s** (`EllipseOutline`, `PolarOutline`): the apparent limb with crossings, normals, areas |
| `fit/` | the **generic inference engine**: `Param`s and priors, `fit` (multistart Levenberg–Marquardt, covariance, error propagation), `profile`/`confidence_interval`, `sample` (ensemble MCMC) |
| `photometry/` | `LightCurve` (immutable), diffraction model, `Atmosphere` ray tracing, `TimingProblem` |
| `event/` | `Contact`, `Chord`, `Event`, and **`Track`s**: each observer's path (f, g)(t) from the full reduction, stored as Chebyshev polynomials (≪ 1 m), so that fits never touch ephemerides |
| `limb/` | `LimbModel`s and the `LimbProblem` (contacts, negative chords, light curves, amplitudes, clock offsets) and their results |
| `plots/` | Plots recipes; maps with GeoMakie (`geomap`, module `GeoMaps`); Makie mesh viewer (module `MeshViews`). `ext/` (package extensions, loaded with the optional package): HTTP services, SPICE |

Every fit is an `AbstractProblem` (parameters + normalised residuals), so all of them share
the same solver, uncertainties and Bayesian sampler: a single light-curve timing fit, a 2-D
limb fit, a joint light-curve/atmosphere fit and a multi-event 3-D fit are used identically.
Results are returned, never stored in the inputs.

## Conventions

* `Occultation` is an alias of `Event` (one occultation: star, body, reference epoch, chords).
* Times: `Epoch` (UTC); light curves in seconds from their `tref`; TDB internally.
* Fundamental plane: (f, g) is the position of the star relative to the body (ephemeris)
  centre seen by the observer, in km at the body, f towards east, g towards north. A fitted
  centre (f₀, g₀) is the offset of the body from its ephemeris; `astrometry` converts it to
  `Δα cos δ = f₀/D`, `Δδ = g₀/D`.
* Position angles from north towards east; an ellipse's `position_angle` is that of its
  minor axis (the apparent pole of an oblate body).
* Body frame: IAU convention [Archinal2018](@citep) (pole α, δ; prime meridian W = W0 + 360°·24/P·d), right-handed.
