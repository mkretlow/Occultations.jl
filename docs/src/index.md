# Occultations.jl

Reduction and analysis of stellar occultations by solar-system bodies in Julia:

* **light curves** – readers (text, PyMovie, Tangra, AstroImageJ), normalisation, event detection, and a
  physical model with Fresnel diffraction (integrated over the passband), the projected
  stellar disc, the exposure time, finite opacity and the diffraction of short (grazing)
  events; fitted timings with profiled flux levels, covariance, χ² profiles or posterior;
* **atmospheres** – hydrostatic, spherically symmetric atmospheres (N₂, CH₄, CO; isothermal,
  Pluto-like inversion layer, free molecular flow) ray-traced to the refracted flux in the
  shadow (near and far limb, limb-curvature focusing, surface edge), e.g. Pluto or
  (612533) 2002 XV93 [Arimatsu2026](@citep); the surface pressure is a fit parameter;
* **limb fits** – one or several events; a 2-D ellipse (default), circle, triaxial
  **3-D ellipsoid** with pole and rotation phase (the apparent ellipse of every event follows
  from the body orientation at that epoch), or the silhouette of a **triangle mesh** (OBJ
  files, or a generated ellipsoid mesh); timed contacts with residuals measured along the
  chords, **negative chords** as hard (default) or soft constraints, partial chords, clock offsets,
  **rotational light-curve amplitudes** as constraints, Gaussian priors; or a **joint fit of
  the light curves themselves** with the limb (diffraction or atmosphere model);
* **uncertainties everywhere** – covariance matrix, linear propagation to derived
  quantities (equivalent radius, axes, apparent ellipses, astrometric offsets), χ² profiles
  and confidence intervals, and **Bayesian posterior sampling** (affine-invariant ensemble
  sampler) with samples of all derived quantities;
* **astrometry** – the body's ICRS position from the fitted centre, including the star's
  position uncertainty;
* **geometry and predictions** – ERFA-based reduction (UTC/TT/TDB, IAU 2006/2000A Earth
  orientation, light time per observer, stellar proper motion and parallax), ground and
  space observers, MPC codes, JPL Horizons, SPICE kernels, Gaia, SBDB, Keplerian orbits,
  Besselian/CORA elements, closest approach and shadow paths on the rotating Earth;
* **plots and maps** – Plots.jl recipes (light curves, fits, chords with 1σ limb draws,
  posterior limbs, corner plots with credible regions and reference values, atmosphere
  profiles), 3-D views of shape models and of 3-D fits in the geometry of the occultation
  (`meshview`, `occview`, `occanimate`; Makie) and publication-quality maps
  with **GeoMakie** (any PROJ projection, Natural Earth coastlines, borders and capitals –
  the Julia counterpart of cartopy).

The package reduces stellar occultations from the light curves to the size, shape,
astrometry and atmosphere of the occulting body, with uncertainties throughout. Every fit
uses the same inference engine, so light-curve timings, limb fits, 3-D shape fits and
atmosphere fits all come with covariances, likelihood profiles and Bayesian posteriors.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/mkretlow/Occultations.jl")
# or, for development, a local clone:  Pkg.develop(path="/path/to/Occultations.jl")
```

GeoMakie (maps) and Makie are dependencies. Optional: `Plots` (recipes), a Makie backend to
display or save maps (`CairoMakie` for files, `GLMakie` for interactive windows), `HTTP`
(Horizons, Gaia, SBDB, MPC), `SPICE` (kernels).

## Contents of this manual

* **User guide:** a [quick start](guide/quickstart.md), then
  [light curves](guide/lightcurves.md), [limb fits](guide/limbfits.md),
  [binaries and contact binaries](guide/binaries.md), [ring systems](guide/rings.md),
  [atmospheres](guide/atmospheres.md), [plots and maps](guide/plots.md) and
  [the occultation in 3-D](guide/occview.md); the [package design and conventions](guide/design.md).
* **Methods:** [the mathematics of the limb fit](theory/limbfit.md),
  [residuals along the chords or radial](theory/residuals.md),
  [the limb allowance σ_limb](theory/sigma_limb.md) and
  [the uncertainty of the limb in the plots](theory/uncertainty_bands.md).
* **[Examples](examples.md):** ten examples, nine of them published occultations reduced
  again, and their [comparison with the publications](comparison/publications.md) and
  [with SORA](comparison/sora.md).
* **[Validation](validation.md)**, the [API reference](api.md) and the
  [reference papers](references.md).

## Licence and citation

Occultations.jl is released under the [MIT licence](https://github.com/mkretlow/Occultations.jl/blob/main/LICENSE). If you use it,
please cite it as given in [`CITATION.cff`](https://github.com/mkretlow/Occultations.jl/blob/main/CITATION.cff). Changes are listed in
the [changelog](https://github.com/mkretlow/Occultations.jl/blob/main/CHANGELOG.md).

The maps use Natural Earth data (public domain) through GeoMakie.
