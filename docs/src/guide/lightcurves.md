# Light curves and timing

## Reading light curves

```julia
lc = read_pymovie("pymovie.csv")                        # date from the header (or date="YYYY-MM-DD")
lc = read_tangra("tangra.csv"; date="2024-01-10")        # SignalMinusBackground (1)
lc = read_astroimagej("measurements.tbl"; target=1)      # JD_UTC, rel_flux_T1 ± rel_flux_err_T1
lc = read_lightcurve("lc.txt"; timescale=:jd, cols=(1, 2, 3))   # generic text file
```

`read_lightcurve` reads generic text files; the readers for PyMovie, Tangra and AstroImageJ
take the columns and the date from the files of these programs. `normalize`, `select`,
`detect_occultation` and `simulate_lightcurve` prepare, cut, search and simulate light
curves (see the [API reference](../api.md#Light-curves-and-photometry)).

## Fitting the contact times

`fit_timing(lc; vel, dist, star_diam)` fits the immersion and emersion times with a physical
model of the light curve: Fresnel diffraction integrated over the passband, the projected
stellar disc and the exposure time, with the flux levels profiled (fitted at each step by
linear least squares). The result has the covariance of the times; `sample` gives their
posterior. Short (grazing) events and partial transparency are part of the model.

## Reading the timing-fit plot

`plot(fit_timing(lc; ...))`:
* The light curve (points) with the model, which includes diffraction, the stellar diameter
  and the exposure.
* Dotted vertical lines mark the fitted immersion and emersion times.
* Residuals are shown below, offset; they should scatter without structure. Systematic
  wiggles at the contacts point to a wrong velocity, stellar diameter or exposure time.

![Timing fit](../assets/figures/timing_fit.png)
