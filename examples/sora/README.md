# Occultations.jl and SORA on the same data

Scripts for the comparison in the manual (Comparison → Occultations.jl and SORA): published
occultations reduced with both packages, and their computing time. The outputs are in
[`examples/output/sora/`](../output/sora/).

SORA is a Python package and not a dependency of Occultations.jl. Install it in its own
virtual environment:

    python3 -m venv <venv>
    <venv>/bin/pip install sora-astro==0.3.3

## Published events (paper – Occultations.jl – SORA)

    julia --project=. examples/sora/export.jl                  # events, Occultations.jl fits
    <venv>/bin/python -I examples/sora/compare_sora.py         # SORA fits
    <venv>/bin/python -I examples/sora/table.py > examples/output/sora/table.md

* `export.jl` writes, for each case, `examples/output/sora/cases/<case>/`: `ephem.txt` (the
  stored Horizons table as JD UTC, RA, Dec, distance in au, for SORA's `EphemPlanete`),
  `event.txt` (the star with its proper motion and parallax, the used contact times of the
  positive chords) and `config.txt` (σ_limb and the Occultations.jl fit: radial residuals, no
  negative chords).
* `compare_sora.py` builds the SORA `Occultation` from these files and writes `sora.txt`:
  SORA's least squares from a rough start and from the Occultations.jl solution; the
  minimum of Σ r² and of Σ r⁴ with scipy from the same start and bounds; SORA's Monte Carlo
  region (`dchi_min=1`, 2000 samples, ±4σ around the better solution); SORA's χ² at the
  Occultations.jl solution. SORA's random draws are not seeded, so its Monte Carlo results
  vary slightly between runs.
* `table.py` prints the tables of the manual.

## Computing time (Varda, 2018 September 10)

    julia --project=. examples/sora/bench_julia.jl             # first: writes the shared inputs
    <venv>/bin/python -I examples/sora/bench_sora.py

Event set-up, the (f, g) of the contacts, the ellipse fit with its uncertainties and the
timing of one synthetic light curve; medians of repeated calls, Julia's first call (with
compilation) reported separately. Outputs in `examples/output/sora/benchmark/`.

## Notes on SORA (0.3.3, and the master branch of 2026-08-03)

* `EphemPlanete` returns geocentric positions for any observer. With SORA's default chord
  projection (topocentric) every station then falls on the geocentre; the scripts set
  `occ.chords._method = "geocenter"`, the projection meant for local tables.
* `Star(..., local=True)` queries VizieR for NOMAD magnitudes unless `nomad=False`.
* From Python 3.14 the default start method of `multiprocessing` is forkserver, which
  re-imports the calling script; the scripts use `fork` (`fit_ellipse` always opens a process
  pool).
* `fit_ellipse`'s least-squares stage and `LightCurve.occ_lcfit(method='least_squares')`
  minimise Σ r⁴ instead of χ² = Σ r²: their residual functions (`ellipseError`,
  `_occ_model_fitError`) return squared residuals, which `scipy.optimize.least_squares` squares
  again.
