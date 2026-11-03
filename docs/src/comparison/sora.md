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
SORA 0.3.3 was also run on nine of the published events and on a speed benchmark, with the
same input as Occultations.jl (sections [The same events reduced with both
packages](#The-same-events-reduced-with-both-packages) and [Computing time](#Computing-time)).

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
    refinement (`fastchi`) around the solution; the least-squares stage minimises Σ r⁴
    instead of χ² (see [below](#The-same-events-reduced-with-both-packages));
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

## The same events reduced with both packages

Nine limb fits of the examples were repeated with SORA 0.3.3 on identical input
([`examples/sora/`](https://github.com/mkretlow/Occultations.jl/blob/main/examples/sora/README.md), outputs in
[`examples/output/sora/`](https://github.com/mkretlow/Occultations.jl/blob/main/examples/output/sora/table.md)):

* the stored JPL Horizons table of each event (in SORA as `EphemPlanete`), the same star
  (position, proper motion, parallax), the positive chords and contact times of the example;
* radial residuals (SORA's definition), the same limb allowance (`σ_limb` = SORA's
  `ellipse_error`), no negative chords (`fit_ellipse` does not use them);
* Occultations.jl: `fit_limb(ev; residual=:radial, σ_limb, negative=:none)`, values ± covariance
  errors;
* SORA: `fit_ellipse(…, method='least_squares', dchi_min=1, number_chi=2000)` in a window of
  ±4σ around the solution, values and errors as SORA reports them (midpoint ± half range of
  the samples with Δχ² < 1).

Each cell gives published / Occultations.jl / SORA. PA is the position angle of the minor axis
(north towards east, modulo 180°). Published 3σ errors (2002 MS4) are divided by 3, the full
axes of 2007 UK126 halved. The published fits differ in details (chords, offsets, σ_limb; see
[Comparison with published reductions](publications.md)); the 1995 SM55 paper used no allowance.

| Event | σ_limb [km] | N | a′ [km] | ε′ | PA [°] |
|---|---|---|---|---|---|
| Varda 2018 | 0.0 | 10 | 383.0 ± 3.0 / 383.4 ± 3.6 / 383.4 ± 3.1 | 0.066 ± 0.047 / 0.071 ± 0.054 / 0.069 ± 0.048 | 67.0 ± 8.0 / 72.4 ± 4.7 / 68.8 ± 7.1 |
| 1995 SM55 2024 | 0.0 | 14 | 104.3 ± 0.4 / 102.8 ± 1.4 / 100.9 ± 0.4 | 0.200 ± 0.002 / 0.201 ± 0.006 / 0.230 ± 0.000 | 44.1 ± 0.4 / 41.0 ± 1.8 / 36.1 ± 0.4 |
| 1995 SM55 2024, σ_limb | 13.5 | 14 | 104.3 ± 0.4 / 106.4 ± 27.7 / 112.6 ± 21.3 | 0.200 ± 0.002 / 0.176 ± 0.232 / 0.175 ± 0.175 | 44.1 ± 0.4 / 95.4 ± 40.3 / 117.9 ± 88.3 |
| Huya 2023 | 0.0 | 22 | 222.5 ± 9.1 / 224.3 ± 12.1 / 226.3 ± 12.0 | 0.107 ± 0.075 / 0.116 ± 0.064 / 0.122 ± 0.063 | 31.5 ± 15.5 / 45.4 ± 27.6 / 43.4 ± 22.5 |
| 2002 TC302 2018 | 10.0 | 24 | 271.6 ± 9.0 / 270.8 ± 12.7 / 274.1 ± 12.7 | 0.154 ± 0.035 / 0.147 ± 0.033 / 0.152 ± 0.032 | 3.0 ± 1.0 / 7.6 ± 13.8 / 7.3 ± 12.4 |
| 2007 UK126 2014 | 0.0 | 6 | 322.9 ± 2.8 / 322.9 ± 4.1 / 323.9 ± 3.8 | 0.074 ± 0.021 / 0.078 ± 0.071 / 0.072 ± 0.053 | 158.8 ± 5.7 / 159.5 ± 19.2 / 151.9 ± 26.9 |
| 2002 GZ32 2017 | 0.0 | 10 | 146.0 ± 12.0 / 140.2 ± 10.0 / 141.4 ± 7.3 | 0.500 ± 0.045 / 0.510 ± 0.028 / 0.510 ± 0.025 | 138.0 ± 3.0 / 143.3 ± 6.2 / 142.9 ± 4.4 |
| Makemake 2011 (circle) | 0.0 | 16 | 715.0 ± 4.5 / 714.9 ± 4.8 / 714.9 ± 4.7 | – | – |
| 2002 MS4 2020 | 7.0 | 26 | 412.0 ± 3.3 / 411.3 ± 3.3 / 411.3 ± 3.3 | 0.066 ± 0.011 / 0.066 ± 0.012 / 0.066 ± 0.011 | 121.0 ± 5.3 / 121.2 ± 5.3 / 121.2 ± 4.8 |

**The two packages agree.** a′, ε′ and the position angle agree within a fraction of their
errors, and both agree with the publications as discussed in
[Comparison with published reductions](publications.md). The exception is 1995 SM55 with the
timing errors alone (χ²/dof ≈ 40), where SORA's Monte Carlo stage does not reach the minimum
(below).

χ² of the same data, radial residuals:

| Event | Occultations.jl minimum | SORA's formula there | SORA least squares | scipy, Σ r² | scipy, Σ r⁴ | SORA Monte Carlo minimum | samples with Δχ² < 1 |
|---|---|---|---|---|---|---|---|
| Varda 2018 | 7.46 | 7.46 | 8.76 | 7.46 | 8.76 | 7.61 | 232 of 2002 |
| 1995 SM55 2024 | 353.27 | 360.97 | 457.19 | 352.71 | 457.19 | 375.67 | 2 of 2000 |
| 1995 SM55 2024, σ_limb | 9.00 | 9.00 | 10.53 | 9.01 | 10.53 | 9.08 | 297 of 2004 |
| Huya 2023 | 12.75 | 12.75 | 13.43 | 12.75 | 13.43 | 12.86 | 524 of 2001 |
| 2002 TC302 2018 | 19.00 | 19.00 | 19.64 | 19.00 | 19.64 | 19.23 | 199 of 2000 |
| 2007 UK126 2014 | 0.03 | 0.03 | 0.03 | 0.03 | 0.03 | 0.03 | 119 of 2006 |
| 2002 GZ32 2017 | 43.27 | 43.34 | 52.85 | 43.28 | 52.85 | 43.51 | 66 of 2000 |
| Makemake 2011 (circle) | 7.36 | 7.38 | 7.36 | 7.36 | 7.36 | 7.36 | 1619 of 2013 |
| 2002 MS4 2020 | 19.12 | 19.13 | 19.98 | 19.13 | 19.98 | 19.23 | 540 of 2002 |

* **The same χ².** SORA's formula evaluated at the Occultations.jl solution gives the same χ²,
  to < 0.1 % where the contact errors are ≳ 1 km. For 1995 SM55 (contact errors ≈ 0.15 km) it
  differs by 2 %: the contact points of the two packages differ by up to ≈ 0.1 km, which is
  negligible elsewhere but not at this precision.
* **SORA's least squares minimises Σ r⁴.** Its residual function for the ellipse
  (`ellipseError` in `occultation/fitting.py`) returns the squared normalised residuals,
  ((Δf)² + (Δg)²)/σ², and `scipy.optimize.least_squares` squares them again. Minimising Σ r⁴
  with scipy from the same start reproduces SORA's least-squares χ² in every case (columns
  "SORA least squares" and "scipy, Σ r⁴"), while Σ r² gives the Occultations.jl minimum. The
  least-squares solution therefore has a higher χ² than the minimum (Varda 8.76 vs 7.46,
  2002 GZ32 52.9 vs 43.3, 1995 SM55 457 vs 353); only for the circle of Makemake are the two
  nearly the same. The residual of the light-curve fit (`_occ_model_fitError` in
  `lightcurve/core.py`, `occ_lcfit(method='least_squares')`) has the same form. The master
  branch of 2026-08-03 has the same code.
* **The Monte Carlo stage** samples around the least-squares solution and mostly repairs its
  bias, but does not reach the minimum where χ² is steep: for 1995 SM55 only 2 of 2000 samples
  lie within Δχ² < 1, so SORA reports ε′ = 0.230 ± 0.000 and a′ ± 0.4 km. Errors of this size
  are published for this event (ε′ = 0.200 ± 0.002); with χ²/dof ≈ 40 they are not meaningful.
  With a limb allowance (σ_limb = 13.5 km, χ²/dof = 1) both packages give errors fifty times
  larger.
* SORA's random draws are not seeded; its Monte Carlo results change slightly from run to run.

Setting up SORA for these runs needed care: with `EphemPlanete` (a local table) SORA's default
chord projection ignores the observer and puts every station on the geocentre, so the
scripts use the older "geocenter" projection (`occ.chords._method = "geocenter"`); and
`Star(…, local=True)` still queries VizieR unless `nomad=False`.

## Computing time

The Varda event of 2018 September 10 (five positive chords) and one synthetic light curve
(0.5 s exposures, noise 0.05), with the same input in both packages
([`examples/sora/bench_julia.jl`](https://github.com/mkretlow/Occultations.jl/blob/main/examples/sora/bench_julia.jl),
[`bench_sora.py`](https://github.com/mkretlow/Occultations.jl/blob/main/examples/sora/bench_sora.py); outputs in
[`examples/output/sora/benchmark/`](https://github.com/mkretlow/Occultations.jl/blob/main/examples/output/sora/benchmark/julia_output.txt)).
Intel i7-1165G7 (2.8 GHz), one thread each; Julia 1.13.2; SORA 0.3.3 with Python 3.14.8 and
NumPy 2.5.4. Medians of repeated calls; Julia compiles each function on its first call, which
takes 0.5–12 s once per session (in parentheses).

| Task | Occultations.jl | SORA | Ratio |
|---|---|---|---|
| Event set-up (star, body, chords) | 0.97 ms (1.6 s) | 133 ms | 140 |
| (f, g) of the 10 contacts | 0.01 ms (0.5 s) | 77 ms | ≈ 7000 |
| Ellipse: best fit | 0.96 ms (12 s), LM + covariance | 190 ms, least squares (of Σ r⁴) | 200 |
| Ellipse: uncertainty region | 141 ms (1.5 s), posterior 20 × 2000; 57 ms, χ² profiles of 5 parameters | 1.6 s (10⁴ samples, 16 within Δχ² < 1); 2.7 s (10³ samples within Δχ² < 1 of the running minimum, 36 of the final one) | 10–50 |
| Ellipse with 6 negative chords (hard constraints) | 40 ms | – (post-filter only) | |
| Light-curve timing | 35 ms (4.3 s) | 84 s (Monte Carlo), 86 s (least squares + Monte Carlo) | ≈ 2400 |

The results agree: ellipse a′ = 383.4 ± 3.6 km (Occultations.jl) and 382.8 ± 2.4 km (SORA,
Monte Carlo); contact times −31.741 ± 0.033 and +31.718 ± 0.035 s (Occultations.jl) and
−31.733 ± 0.083 and +31.714 ± 0.053 s (SORA), the truth being ∓31.700 s. SORA's timing errors
are about twice the covariance errors here; a coverage test would decide which are right.

Speed matters where fits are repeated: posteriors, χ² profiles, joint light-curve fits,
multi-event 3-D fits, and coverage tests such as those of the synthetic examples, which
refit hundreds of simulated events.

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
| Minimisation | least squares (of Σ r⁴, see above) / differential evolution + Monte Carlo, or Monte Carlo only | Levenberg–Marquardt, multiple starts, automatic differentiation |
| Reported value | midpoint of the Δχ² < σ² range | χ² minimum (fit), or posterior median |
| Uncertainty | half range of sampled solutions with Δχ² < σ² | covariance, likelihood profiles, MCMC posterior quantiles |
| Excess scatter | `ellipse_error` (fixed by the user) | `σ_limb` (fixed or `:auto`), clock offsets with priors, `scale_errors` |
| Negative chords | post-filter of the accepted solutions | in the likelihood (hard or soft) |
| Shape models | centre and scale, fixed orientation, one event | centre, scale, pole, rotation; several events; ellipsoid with amplitudes and phases |
| Light curves | timing fit with fixed flux levels; Monte Carlo errors | timing fit with profiled flux levels; joint light-curve + limb (+ atmosphere) fits |
| Bayesian analysis | – | posteriors, priors, corner plots, credible regions |

Both packages reduce the same geometry: tracks, the fundamental plane, Fresnel diffraction
and astrometry. On the same data the best-fitting limbs agree (see
[The same events reduced with both packages](#The-same-events-reduced-with-both-packages)), and
Occultations.jl is two to three orders of magnitude faster (see [Computing time](#Computing-time)). The differences
lie in the uncertainties, the treatment of negative chords and the joint 3-D and photometric
models.
