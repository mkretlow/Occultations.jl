# The limb allowance σ_limb

The timing error σₜ of a contact describes only the measurement of the time. Real contacts
also deviate from the fitted limb for other reasons:
* **topography:** mountains and depressions of a few km, which an ellipse cannot follow;
* **shape:** a limb that is not exactly an ellipse;
* **timing:** errors not captured by σₜ, such as clock offsets, latencies, or a wrong
  exposure time.

`σ_limb` [km] is an allowance for all of these. It is added in quadrature to the uncertainty
of every contact:

    σ = √((σₜ |v|)² + σ_limb²)

* **Units:** km at the body, i.e. along the chord for `residual=:chord` and radially for
  `residual=:radial`.
* **Scope:** it applies to the contacts only. Negative chords keep their own tolerance
  (hard, or `sensitivity × speed` for `:soft`).
* **Values:**
  * `σ_limb = 0` (the default) uses the timing errors alone;
  * a number fixes it, e.g. from known topography ([Kilic2026](@citet) use 4 km for Ixion);
  * `fit_limb(...; σ_limb=:auto)` finds the value for which χ²/dof = 1 (by bisection)
    and fits with it.

**What it does to the fit.**
* **Weights:** they become more equal. Without it, a contact timed to 0.01 s dominates one
  timed to 1 s by a factor of 10⁴ in weight. With σ_limb = 5 km, neither can be trusted
  better than 5 km, which is what topography or a non-elliptical limb implies.
* **Errors:** a, b, the position angle, the centre and everything derived from them get
  larger, consistent with the actual scatter of the contacts about the limb.
* **Scaling:** `scale_errors = true` multiplies all uncertainties by √(χ²/dof) and keeps the
  relative weights. σ_limb instead adds the same amount everywhere, which suits scatter that
  does not depend on the timing precision.
* **Systematic offsets:** a clock offset of one station is better modelled by its own
  parameter (`fit_dt`, with a prior). σ_limb is for random deviations.
* **SORA:** the corresponding option of `fit_ellipse` is `ellipse_error`, fixed by the user.

**How to read it.** σ_limb from `:auto` estimates the rms deviation of the contacts from the
fitted limb, beyond their timing errors. It is a property of the data together with the
model:

| Example | σ_limb | Meaning |
|---|---|---|
| Huya 2019 | 7–10 km | timing problems of some stations, and possibly topography |
| 1995 SM55 | 13.5 km | driven by two inconsistent Wise timings |

With few degrees of freedom the estimate itself is uncertain: with one degree of freedom
(2007 UK126) it means little.

**In the plots.**
* The red contact segments show ±1σ of the **timing** only, i.e. σₜ |v|, without σ_limb, so
  that the timing precision of each station stays visible.
* The bands (shape band, pointwise band, 68.3 % region) and the centre ellipse come from the
  fit's covariance or posterior, which include σ_limb. They are wider than with the timing
  errors alone.
* A limb may therefore lie outside a station's red segment and still be consistent: the
  difference is absorbed by σ_limb.
* When σ_limb > 0, its value is printed with the fitted dimensions below the legend
  ("σ_limb = … km").
