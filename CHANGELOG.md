# Changelog

All notable changes to Occultations.jl are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the package uses
[semantic versioning](https://semver.org/).

## [Unreleased]

### Added

- Real-data examples: 2002 TC302 (Ortiz et al. 2020), 2003 AZ84 (Dias-Oliveira et al. 2017;
  2012 and 2014, and the single-chord events of 2011 and 2013 with their stars identified in
  Gaia DR3), 2003 VS2 (Vara-Lubiano et al. 2022), Haumea and its ring (Ortiz et al. 2017),
  2002 GZ32 (Santos-Sanz et al. 2021), Makemake (Ortiz et al. 2012), Bienor (Fernández-Valenzuela
  et al. 2023; Rizos et al. 2024), Chariklo and its rings (Morgado et al. 2021), 2002 MS4
  (Rommel et al. 2023), Quaoar and its rings (Braga-Ribas et al. 2013; Pereira et al. 2023) and
  Chiron (Braga-Ribas et al. 2023), with their data (CDS tables where published there),
  stored ephemerides and manual pages.
- `σ_ring`: an uncertainty allowance for ring contacts (default `σ_limb`), for rings that are
  not exactly circular.
- Ring detections given by their mid-times: `Contact(:mid, t, σ; body=:ring1)`; a chord that
  saw only ring dips is a negative chord for the body.
- **Ring systems**: `RingedLimb(body; rings=n)` – circular rings in the plane of the ring pole
  around any body model (ring pole, radius, width and normal optical depth of each ring);
  `ProjectedRing`, `project_ring`, `ring_geometry`; ring dips as contacts labelled
  `body=:ring1` …; partly transparent strips with diffraction in the light-curve model
  (`occultation_flux`/`simulate_lightcurve` with an opacity per drop) and in joint light-curve
  fits; derived opening and position angles, apparent opacities, equivalent widths and depths;
  rings in the limb plots; synthetic examples `sim_rings.jl` (two narrow rings, light curves)
  and `sim_ring_pole.jl` (a distant ring in three events); user-guide page and tests.
- **Binaries and contact binaries**:
  - `CompositeOutline`: the union of outlines with their own centres (signed distance,
    occulted `intervals` of a track, union area, boundary as polylines);
  - `BinaryLimb(primary, secondary)`: a body with a satellite or a binary system, any limb
    model per component (default: ellipse and circle; `BinaryLimb(EllipseLimb(), EllipseLimb())`
    fits an ellipse to the satellite); separation, size and mass ratio; astrometry of the
    barycentre or photocentre (`astrometry(fit; reference=:barycenter)`);
  - `ContactBinaryLimb`: two overlapping ellipses along a common axis (neck, waist width);
  - `ContactBinaryEllipsoids` and the shape `BilobedEllipsoid`: two triaxial ellipsoids,
    rotating, with exact silhouettes, union volume, light-curve amplitudes and phases, OBJ
    export and `occview`;
  - chords with several drops (`Chord(...; contacts=[…])`, `occulted_intervals`); contacts
    assigned to a component (`Contact(...; body=:secondary)`) or matched to the nearest
    crossing; the visible segments of positive chords as constraints (`visible`, `:soft` by
    default for the composite models, `:none` otherwise);
  - light curves with several drops (`occultation_flux`, `simulate_lightcurve` with a list of
    drops; close drops by adding the complex amplitudes of the strips) and joint light-curve
    fits of composite limbs;
  - limb plots of composite fits (occupancy band, component centres);
  - synthetic examples `sim_binary.jl`, `sim_contact_binary.jl`, `sim_contact_binary_3d.jl`
    and the real example `huya_2023_satellite.jl` (Huya and its satellite, Rommel et al.
    2025), with manual pages, a user-guide page, and tests;
  - `detect_occultations(lc; method=:threshold)`: drops as runs of points below a level,
    for drops separated by short gaps.
- **Manual** built with Documenter.jl (`docs/make.jl`, sources in `docs/src/`): user guide,
  methods, examples, comparisons, validation, API reference and reference papers. The
  figures are copied from `examples/output/`. Doctests in the docstrings of `Epoch`,
  `parse_ra`, `ellipse_outline`, `projected_axes_dunham`, `rotational_amplitude` and
  `AbstractProblem` (a complete fit of a user-defined problem).
- **Example pages** in the manual, one per example, made from the scripts with Literate.jl:
  the comments of a script are the text of its page, each `section("…")` is a heading,
  and the output stored in `examples/output/<example>/output.txt` follows each section.
  The code is not run when the manual is built. The scripts' code is unchanged.
- **Bibliography** with DocumenterCitations.jl: `docs/src/refs.bib` (19 entries from their
  DOIs) and author-year citations in the manual and the example pages.
- Docstrings for `parameters`, `residuals`, `derived`, `errors`, `freenames`,
  `reduced_chi2`, `ellipse_outline` and `enable`; every exported name is now documented.
- Package infrastructure: `LICENSE` (MIT), `CITATION.cff` and this changelog.
- Real-data examples:
  - 1995 SM55 (Ortiz et al. 2025);
  - Huya 2023 with its satellite (Rommel et al. 2025);
  - Ixion, six events (Kilic et al. 2026);
  - 2007 UK126 (Schindler et al. 2017);
  - Varda (Souami et al. 2020).

  Each example comes with its data file, stored JPL Horizons ephemerides and a section in
  `docs/src/comparison/publications.md`.
- Example 2002 XV93 (Arimatsu et al. 2026, Nature Astronomy), based on the authors' light curves
  (Zenodo, CC BY 4.0). It covers:
  - timings;
  - the circular shadow;
  - the ray-traced CH₄/N₂/CO models compared with the authors' profiles;
  - a joint fit of the light curves with the limb and the pressure;
  - the test for an atmosphere.
- `atmosphere_test`: profile likelihood in the surface pressure. It gives the detection
  significance, one-sided limits and upper limits at a chosen confidence, and accepts
  `airless_init` as the start of the p = 0 refit.
- `RotationalPhase` and `rotational_phase`: the rotational phase at the occultation time as a
  constraint on 3-D ellipsoid fits. The UY117 example uses it.
- `LogPrior` and `isotropic_pole` priors; `summarize(posterior)`.
- `meshview`: a shape model in its body frame, or as seen from the Earth.
- `occview` and `occanimate`: 3-D limb fits and posteriors in the geometry of the occultation.
  - The body at the event, lit by the Sun, with the pole, the sub-observer and sub-solar
    points and the chords in front of it.
  - The silhouette with the chords and posterior silhouettes.
  - The rotational light curve from the projected area, with the observed amplitudes.
  - Sliders for the rotational phase, the event and, with `explore`, the pole, W0 and the
    size, with the live χ².
  - Animations through one rotation.
- `cornerplot` gains:
  - credible-region contours and the maximum-posterior sample;
  - reference values (`truths`, `truth_label`) and `axislabels`;
  - a legend with the sampler statistics.
- `geomap` (GeoMakie extension) gains:
  - Natural Earth coastlines at 110, 50 and 10 m;
  - country borders and capitals;
  - the shaded shadow band;
  - numbered stations with a station legend (`observer_labels`).
- Limb plots gain:
  - the station number at each chord (`numbers`);
  - the contact-time uncertainty at the chord ends;
  - the uncertainty of the limb, in both the fit and the posterior plots:
    - `uncertainty=:shape` (the default): the 1σ band of the shape, consistent with σ_a and
      σ_b, plus the 1σ error ellipse of the centre;
    - `:pointwise`: like SORA's 1σ envelope;
    - `:region`: the simultaneous 68.3 % region of whole limbs;
    - `:draws`: the individual limbs;
    - `nsigma` sets the level of each.
  - each example also writes a plot in the style of its paper (`*_paper.png`).
- The fitted 3-D ellipsoids are written as OBJ files to the example output (`write_obj`).
  `ellipsoid_mesh` gains a `name` keyword.
- Limb-fit and posterior plots print the dimensions with 1σ errors in the right margin below the legend:
  a, b and D_eq for an ellipse, or R and D for a circle (`dimensions=true`, the default),
  and the limb allowance σ_limb if it is used.
- Automatic differentiation (ForwardDiff) for the Jacobians of the limb fits (ellipse,
  circle, 3-D ellipsoid) and the timing fits, with an automatic fallback to finite
  differences (joint light-curve, atmosphere, mesh fits); `fit(...; ad=false)`.
- `read_astroimagej` for AstroImageJ measurement tables; `read_lightcurve(format=:astroimagej)`.
- Prediction maps: `occ_map(event; error_mas)` takes the cross-track ephemeris uncertainty in mas.
  `geomap(prediction, fit)` draws the predicted and the post-fit shadow paths together.
  `geomap(...; overlays, overlay_labels, overlay_colors, path_label)` overlays further paths.

- `Occultation`, an alias of `Event`.
- `docs/src/comparison/sora.md`: SORA's fitting methods compared with this package (based on its
  source code).
- `README.md`: a "Reading the plots" section and an explanation of limb fits versus
  posteriors.
- `CITATION.cff` with the ORCID, the affiliation and the repository URL. The installation
  instructions use the GitHub URL.

### Changed

- Maps: grid spacing adapted to the map extent (1–30°, at most about six lines), "0°" instead
  of GeoMakie's "1.9e-6°"; the shadow band follows the Earth's limb where only part of the
  shadow falls on the Earth (it was cut off between the limit lines); globe legend text 10 pt;
  night-side shading 0.38.
- Light curves: data black, model red, residuals grey (as SORA). occview: the sunlit side of
  the body light grey. Corner plots: a title (the event), the reference values (`truths`) as
  an orange line in the histograms and a diamond in the 2-D panels, both in the legend.
- Quick start: the Varda event (a 2-D ellipse close to the published one) instead of 2002 XV93;
  no mesh viewer.
- Maps (`geomap`): a larger title (`titlesize`, default 20) with more space above the frame;
  on the globe no grid labels (`ticklabels`) and a compact legend; time marks on whole UTC
  minutes, the first labelled in full, the others with the minute only; an arrow for the
  direction of motion; a red star where the star stands highest on the centre line (at the
  zenith, if reached; `OccMap.zenith`); the shadow band clipped to the visible hemisphere.
- Examples, their data, outputs and figures are named with the full provisional designation:
  `2002kx14_2020_2023.jl`, `2002xv93_2024.jl`, `2003uy117_2020.jl`, `1995sm55_2024.jl`,
  `2007uk126_2014.jl`, `examples/data/2002xv93_zenodo/`, `2002kx14_limb_fit.png`, … (named
  objects keep their names: Huya, Ixion, Varda).
- `fit` keeps the refit from the canonical form of the parameters (e.g. ε ≥ 0) also when its
  χ² equals the best one within rounding (it could keep the non-canonical form before).
- The README is an overview (features, installation, quick start, examples); the
  detailed text moved to the manual. `docs/comparison_with_publications.md` and
  `docs/sora_comparison.md` moved to `docs/src/comparison/`.
- Repository layout: the reference papers moved from `docs/` to `references/` (with
  `references/README.md`, formerly `docs/references.md`) and are named
  `<object>_<year>_<first author>.pdf` (e.g. `2003UY117_2024_Kretlow.pdf`). The example
  outputs moved from `output/` to `examples/output/`. `docs/` is kept for the documentation.
- Negative chords are hard constraints by default (`negative = :hard`; previously `:soft`).
  The covariance band of the limb-fit plot rejects draws that violate a hard constraint.
- XV93 example: the joint light-curve fits search a deterministic grid of starting radii and
  centres. The χ² surface has several minima, and the random multistarts found different
  ones depending on the seed.
- Examples pass their seeded random-number generator to `fit_timing`, which makes them
  reproducible.
- Every example writes its figures and its printed output to its own folder,
  `examples/output/<example>/`.
- **GeoMakie and Makie are now dependencies.** `geomap` and `meshview` are part of the
  package (modules `GeoMaps`, `MeshViews`), no longer extensions.
- The GSHHG shorelines (`data/`), `coastlines` and the Plots map recipe (`plot(occ_map(...))`)
  are removed; maps are drawn with `geomap`.
- Chord plots (`plot(event)`) zoom on the contacts.
- All plots have a closed frame: the Plots recipes use `framestyle = :box` (axes on all four
  sides), and maps have a frame along the map area.
- Globe maps (orthographic):
  - the Earth's limb is drawn;
  - land is filled and clipped at the horizon on the sphere;
  - the night side is shaded transparently, with the twilight as a gradient;
  - the terminator is drawn (`night`, `terminator`). Other
  projections draw neither by default.
- Every example writes the same standard figures (`examples/standard_figures.jl`):
  - the limb fit (default, 68.3 % region, in the paper's style) and the posterior;
  - the predicted, the post-fit and the combined shadow path;
  - the globe.

- UY117 3-D posterior: the aspect angle is folded to min(ψ, 180° − ψ). The pole flip leaves
  the data unchanged, so only the folded angle is determined.

### Fixed

- Multi-event 3-D fits with rotational `phases` or `amplitudes`: each constraint now uses the
  ephemeris of the event closest in time (previously that of the first event, extrapolated
  far outside its table, which could fail).
- Chord plots: a positive chord with ring mid-times only no longer fails.
- Maps of a limb fit (`occ_map(fit)`, `geomap(fit)`): the shadow limits are the extent of the
  fitted limb across the motion, previously ±R_eq. The ±1σ lines come from the fit covariance
  (centre and shape, separately for each edge), previously only the error of the centre.

- `read_tangra` took the raw `Signal (n)` column even when `SignalMinusBackground (n)` was
  present.

- `atmosphere_test`: the airless (p = 0) refit could stay near the atmospheric solution,
  which overstated the detection. It now uses multistarts and two sweeps of the pressure
  grid.
- The chord-number labels of distant negative chords were drawn outside the plot.

## [0.1.0] – 2026-11-01

### Added

- Initial version with a new architecture, not modelled on SORA. Its layers:
  - core and astro: time, star, ephemerides, observers, geometry;
  - body: outlines, shapes, rotation;
  - a generic fitting engine: least squares with multistarts, covariance, likelihood
    profiles, ensemble MCMC;
  - photometry: diffraction, timing, ray-traced atmospheres;
  - events: chords and shadow tracks;
  - limb models: ellipse, circle, multi-event 3-D ellipsoid with rotational amplitude, mesh;
  - plots.
- Real-data examples:
  - Huya (Santos-Sanz et al. 2022);
  - 2002 KX14, five events with the corrigendum (Rizos et al. 2025, 2026);
  - 2003 UY117 (Kretlow et al. 2024).
- `docs/comparison_with_publications.md`, comparing the examples with the published results.
