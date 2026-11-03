# Occultations.jl and SORA: the fitting methods compared

SORA (Stellar Occultation Reduction and Analysis; [GomesJunior2022](@citet), MNRAS 511, 1167,
[doi:10.1093/mnras/stac032](https://doi.org/10.1093/mnras/stac032);
[arXiv:2201.01799](https://arxiv.org/abs/2201.01799)) is the Python package used for most of
the published reductions reproduced in [Comparison with published reductions](publications.md).

This document compares the fitting methods of the two packages: limb fits, shape-model fits,
negative chords, light-curve timing, and the uncertainties and plots derived from them. The
description of SORA is based on its source code: the master branch of
[github.com/riogroup/SORA](https://github.com/riogroup/SORA), read on 2026-11-02. The
relevant modules are `occultation/fitting.py`, `occultation/utils.py`, `extra/chisquare.py`
and `lightcurve/core.py`. Earlier versions, used by older papers, may differ in details.

## Corresponding functions

| Task | SORA | Occultations.jl |
|---|---|---|
| Ellipse (or circle) fit to the chords | `fit_ellipse(occ...; equatorial_radius, dequatorial_radius, center_f, dcenter_f, …, method, sigma_result, ellipse_error)` | `fit_limb(events; model=EllipseLimb() / CircleLimb(), residual, σ_limb, priors, …)` |
| Several occultations in one fit | `fit_ellipse(occ1, occ2, …)` | `fit_limb([ev1, ev2, …])`: shared shape, one centre per event; per-event parameters and clock offsets optional |
| Negative chords | `filter_negative_chord(chord, chisquare; step, sigma)` (post-filter) | part of the likelihood: `negative = :hard` (default) or `:soft` |
| Shape model (mesh) | `fit_shape(occ; center, scale, …)`, `fit_to_limb(limb, fg, error; …)` | `MeshLimb(mesh; fit_orientation)`: centre, scale, pole, W0, in-plane rotation, over several events |
| Triaxial ellipsoid over several events | – | `EllipsoidLimb()` with `RotationalAmplitude`, `RotationalPhase` |
| Light-curve timing | `LightCurve.occ_lcfit(; loop, method, sigma_result, …)` | `fit_timing(lc; …)` |
| Joint fit of light curves and limb | – | `fit_limb(…; photometry=true)` |
| Atmosphere | – (not in the modules examined) | `Atmosphere`, `atmosphere_profile`, `photometry=true, atmosphere=…`, `atmosphere_test` |
| Posterior (MCMC) | – | `sample(fit)`, `cornerplot`, posterior plots |

## Ellipse fit

**SORA – `fit_ellipse`**

* **Data:** the contact points of the positive chords. The uncertainty of a point is its
  timing error × |v|, plus `ellipse_error` (km) in quadrature.
* **Residual:** radial. The model point lies on the ray from the ellipse centre through the
  data point, and χ² is the sum of the squared distances between the data points and these
  model points, divided by (σᵢ² + `ellipse_error`²).
* **Parameters:** centre (f, g), equatorial radius a, oblateness ε, position angle. Each is
  searched within a window `value ± dvalue`; a zero width fixes the parameter.
* **Search** (`method`):
  * `'least_squares'` (the default): trust-region least squares, then a Monte Carlo
    refinement (`fastchi`) around the solution;
  * `'differential_evolution'`: the same, with differential evolution as the first stage;
  * `'chisqr'` / `'fastchi'`: a pure Monte Carlo search with uniform random draws in the
    windows, until `number_chi` samples are collected.
* **Result** (`ChiSquare.get_nsigma(sigma=sigma_result)`):
  * all samples with χ² < χ²_min + σ² are selected;
  * for each parameter, the result is the **midpoint** of the range of the selected samples,
    and the uncertainty is **half that range**.

**Occultations.jl – `fit_limb`**

* **Residual:** along the chord by default (`residual=:chord`). This is the measurement model:
  the timing residual × |v|, with σ = √((σₜ|v|)² + σ_limb²). A chord that misses the limb is
  penalised by the distance of the miss. `residual=:radial` gives SORA's definition.
* **Search:** Levenberg–Marquardt with multiple starts, a Jacobian by automatic
  differentiation, and the parameters reduced to a canonical form.
* **Result:**
  * the minimum of χ², with errors from the covariance (JᵀJ)⁻¹;
  * likelihood profiles (`profile`);
  * the MCMC posterior (`sample`): the median and the 15.87–84.13 % quantiles of the
    marginal distribution.
* **Derived quantities** (b, R_eq, D_eq, astrometry, …) are propagated with the covariance
  or taken from the posterior.

**Differences that matter when comparing numbers**

1. **The definition of the uncertainty.** SORA's "1σ" is half the range of the parameter
   among all sampled solutions with Δχ² < 1.
   * This is a likelihood-profile interval for one parameter at a time (Δχ² = 1). It is
     correct for a single parameter of a linear, Gaussian problem.
   * It is estimated from the extremes of a finite random sample. Few samples near the
     edges make it too small, and it depends on the windows `d*`, `number_chi` and the
     random draws.
   * The reported value is the midpoint of the interval, not the χ² minimum. For an
     asymmetric χ² surface, this differs from the best fit.
   * Errors of quantities derived afterwards (e.g. ε from a and b, R_eq) are not always
     consistent with those of a and b. Examples are 1995 SM55 and Ixion in
     [Comparison with published reductions](publications.md).
2. **Radial residuals** understate the misfit of oblique and grazing contacts by cos θ
   (see [Residuals along the chords or radial residuals?](../theory/residuals.md)). χ²/dof is therefore
   smaller with radial residuals (2003 UY117: 0.41 radial vs 2.10 along the chords).
3. **No scaling with χ².** With χ²/dof ≫ 1, SORA's Δχ² < 1 intervals are taken at face
   value. The 2019 refit of Huya in [Rommel2025](@citet) has a′ ± 0.11 km at χ²/dof ≈ 211.
   Here, `σ_limb = :auto`, `scale_errors = true` or clock-offset priors make the excess
   scatter explicit.

## Negative chords

**SORA – `filter_negative_chord(chord, chisquare, step=1, sigma=0)`.** After the fit, each
accepted solution is tested against one negative chord:
* the chord's full path is sampled every `step` seconds (`step='exposure'` samples the
  exposures of its light curve);
* points farther than 1.5 a from the centre are ignored;
* the solution is kept only if every remaining point lies outside the ellipse, by more than
  −`sigma` km.

The filter returns a reduced `ChiSquare`, from which the results are recomputed. It has to be
applied to each negative chord explicitly, and the minimisation itself does not know about
negative chords.

**Occultations.jl.** Negative chords are part of the likelihood:
* `:hard`: a limb crossing the observed track is inadmissible (χ² = ∞, zero posterior);
* `:soft`: a penalty scaled by the exposure × speed.

The best fit, the covariance, the profiles and the posterior all respect them.
`check_negative_chords(fit)` lists the distance of each negative chord from the limb.

The examples show why this matters. [Rommel2025](@citet) filtered the 2023 June 24 Huya
solutions with the Belesta negative chord only; their free solution crosses the Guirguillano
negative chord by 41 km (see [Comparison with published reductions](publications.md)).

## Shape models

**SORA – `fit_shape(occ; center_f, dcenter_f, center_g, dcenter_g, scale, dscale, loop, sigma_result, model_error)`**
* The orientation is fixed: the limb of the body's shape model at the time of closest
  approach (`occ.body.shape.get_limb`).
* A Monte Carlo search (`loop` draws) fits only the centre and the scale, with radial
  residuals (`fit_to_limb`).
* One occultation at a time.

**Occultations.jl – `MeshLimb(mesh; fit_orientation=true)`**
* Fits the centre, the scale and the orientation: pole, W0, and optionally an in-plane
  rotation per event.
* Uses all events at once, each with its own light-time corrected rotation phase.
* The triaxial-ellipsoid model (`EllipsoidLimb`) adds rotational light-curve amplitudes and
  phases as constraints. `occview` shows the result in 3-D, with the chords and the
  rotational light curve.

## Light-curve timing

**SORA – `occ_lcfit`**
* **Model:** the Fresnel diffraction of an opaque (or partly transparent) edge, averaged
  over the two edges of the passband, the stellar disc and the exposure time.
* **Flux levels:** fixed from the event detection (`occ_detect`), not fitted.
* **Search:** a Monte Carlo search of immersion and emersion times by default (`'chisqr'`),
  within ±`delta_t`. Least squares or differential evolution followed by Monte Carlo are
  options.
* **Errors:** half the range of the solutions with Δχ² < σ², as for the ellipse.

**Occultations.jl – `fit_timing`**
* **Model:** diffraction integrated over the passband (sampled), the stellar disc and the
  exposure.
* **Flux levels:** profiled, i.e. the linear least-squares solution at each step, so they are
  fitted rather than fixed.
* **Errors:** Levenberg–Marquardt with automatic differentiation, a covariance, and an
  optional posterior (`sample`) for asymmetric or non-Gaussian timing errors.
* **Joint fits:** with `fit_limb(…; photometry=true)` the light curves are fitted together
  with the limb, and with an atmosphere if given. No separate timing step is needed, and
  grazes, partial and negative light curves are treated in the same way.

## Uncertainty regions in the figures

In SORA-based papers, the gray "1σ uncertainty region" around the fitted ellipse is drawn
from the solutions with χ² < χ²_min + 1 ("3σ" with + 9). It is the envelope of the ellipses
compatible with the data for one parameter at a time.

The corresponding display in Occultations.jl is `uncertainty=:pointwise` (with `nsigma=3`
for "3σ"); the examples write it as `*_limb_fit_paper.png`. Its default display
(`uncertainty=:shape`) and the 68.3 % region of whole limbs (`uncertainty=:region`) are
explained in [Uncertainty of the limb in the plots](../theory/uncertainty_bands.md). For 2003 UY117, a whole
limb stays inside the pointwise band only 7 % of the time. The 68.3 % region is about as wide
as SORA's "3σ" envelope.

## Summary

| Aspect | SORA | Occultations.jl |
|---|---|---|
| Limb residual | radial | along the chord (radial optional) |
| Minimisation | least squares / differential evolution + Monte Carlo, or Monte Carlo only | Levenberg–Marquardt, multiple starts, automatic differentiation |
| Reported value | midpoint of the Δχ² < σ² range | χ² minimum (fit), or posterior median |
| Uncertainty | half range of sampled solutions with Δχ² < σ² | covariance, likelihood profiles, MCMC posterior quantiles |
| Excess scatter | `ellipse_error` (fixed by the user) | `σ_limb` (fixed or `:auto`), clock offsets with priors, `scale_errors` |
| Negative chords | post-filter of the accepted solutions | in the likelihood (hard or soft) |
| Shape models | centre and scale, fixed orientation, one event | centre, scale, pole, rotation; several events; ellipsoid with amplitudes and phases |
| Light curves | timing fit with fixed flux levels; Monte Carlo errors | timing fit with profiled flux levels; joint light-curve + limb (+ atmosphere) fits |
| Bayesian analysis | – | posteriors, priors, corner plots, credible regions |

Both packages reduce the same geometry: tracks, the fundamental plane, Fresnel diffraction
and astrometry. Where the published reductions are reproduced, the best-fitting limbs agree
(see [Comparison with published reductions](publications.md)). The differences
lie in the uncertainties, the treatment of negative chords and the joint 3-D and photometric
models.
