# Plots and maps

## Limb fit

`plot(fit; event=1, uncertainty=:region, numbers=true)`. The fundamental
plane as seen from the Earth: f towards east (drawn to the *left*, as on the sky), g towards
north, in km at the body, relative to the ephemeris position. The black curve is the fitted
limb and the cross its centre (the offset from the ephemeris, i.e. the astrometry). Coloured
segments are the occulted parts of the positive chords (`extend=true` adds their tracks
beyond the limb as dotted gray lines); dashed gray lines are negative chords. The number at each chord is its index in the legend (colours
repeat). The red segments at the chord ends ("contact-time uncertainty") are ±1σ of the
contact times, converted to km with the shadow velocity: a limb passing within them fits
that contact.

The gray area shows the uncertainty of the limb. There are three kinds; the keywords are
described in [Uncertainty of the limb in the plots](../theory/uncertainty_bands.md). By default (`uncertainty=:shape`), the gray band is the
**1σ shape band**, and the small black ellipse at the centre cross is the **1σ error ellipse of
the centre**.

The limb can be compared with the chords in three ways:
* the limb should pass through the red segments;
* no negative chord should cross it;
* where no chord constrains the limb, the region widens and can be asymmetric. For example,
  if only the southern limb is observed, the size and centre trade off against each other.

The fitted dimensions with 1σ errors (a, b and D_eq of an ellipse; R and D of a circle;
median −/+ for a posterior) are printed in the right margin below the legend
(`dimensions=false` omits them).

## Posterior limbs

`plot(posterior)`. The same view and the same uncertainty displays,
computed from limbs drawn from the MCMC posterior instead of the fit covariance. The black
curve is the maximum-posterior limb.

![2003 UY117 limb fit](../assets/figures/2003uy117_limb_fit.png)

## Corner plot

`cornerplot(posterior; vars, axislabels, truths)`:

* **Diagonal:** the marginal posterior of each parameter. The title gives the median with
  the 16 % and 84 % quantiles as −/+ errors (dashed lines).
* **Off-diagonal:** the joint posterior of a pair: samples, with contours enclosing 39 % and
  86 % of the probability. These are the 1σ and 2σ regions of a 2-D Gaussian.
* The star marks the maximum-posterior sample. Reference values (`truths`, e.g. a
  publication) are drawn in orange: a vertical line in the histograms, a diamond in the 2-D
  panels; the legend names them (`truth_label`).
* The title is the name of the event (`plot_title` sets another).
* The legend panel gives the number of samples and walkers, the acceptance fraction and the
  autocorrelation time τ. Sampling is reliable if the chain length is much greater than 50τ.
* Elongated or curved contours show correlations and degeneracies (e.g. size against
  aspect). A reference value outside the 86 % contour disagrees with the fit.

![Huya 2019 corner plot](../assets/figures/huya_corner.png)

## 3-D posterior histograms

Marginal posteriors of a triaxial
ellipsoid (a ≥ b ≥ c, rotating about c), fitted to a single occultation with the pole
unknown (isotropic prior). The two histograms compare two sets of constraints:

* **occultation + amplitude** (blue):
  * the chords fix the apparent ellipse at one instant;
  * the light-curve amplitude Δm fixes the elongation the body shows over a rotation,
    √[(a²cos²ψ + c²sin²ψ)/(b²cos²ψ + c²sin²ψ)].

  Together they do not determine a, b, c and the pole. One apparent ellipse can be cut from
  many ellipsoids, seen at different aspects ψ and rotational orientations. The least-squares
  fits form a valley of equally good solutions, and the posterior spreads along it. The
  valley includes long, thin bodies seen nearly end-on, which give the tail to large
  diameters.
* **+ rotational phase** (orange): adds the phase of the light curve at the occultation
  time, i.e. where the body was in its rotation (0.63 rotations after maximum light for
  2003 UY117). The phase says how far the long axis was from the sky plane during the
  occultation, and removes the valley solutions in which the chords cut the body in a
  different orientation.

How to read the panels:

* **Left, volume-equivalent diameter D_eq,V** (the diameter of a sphere of the same volume,
  2(abc)^⅓): the phase cuts the tail of large diameters. For 2003 UY117 the median moves
  from 243 −24 +33 km to 241 −21 +20 km. With hard negative chords, most of the very
  elongated solutions are already excluded by the occultation itself. The paper's
  235 ± 25 km (dashed line) agrees.
* **Right, aspect angle**: the angle between the spin axis and the line of sight (90° =
  equator-on, the full light-curve amplitude; 0° = pole-on, no amplitude).
  * A pole flip (ψ → 180° − ψ) leaves the occultation and the amplitude unchanged, so only
    min(ψ, 180° − ψ) is determined, and that is what is plotted.
  * The amplitude alone sets a lower bound on the aspect. Even a very elongated body must
    be seen far enough from pole-on to produce Δm = 0.36 mag.
  * The phase does not narrow the aspect (2003 UY117: 71 −19 +14° without, 69 −16 +13° with
    the phase; paper 70°). Its effect is on the shape and size.
  * The sampler visits the two mirror poles unequally: 66 % and 59 % of the samples are in
    one hemisphere, instead of 50 %. Folding the aspect removes this sampling artefact.
* The posterior of a single event is wide in any case: only further occultations at
  other rotational phases, or a pole from photometry, determine a, b and c individually.
  A narrow peak would indicate a constraint (or a prior) that is stronger than the data.

![2003 UY117 3-D posterior](../assets/figures/2003uy117_3d_posterior.png)

## Maps

`geomap(occ_map(...))` (GeoMakie):

* **Shadow path:** the black lines are the shadow limits and the blue line is the centre
  line. Red dashed lines are offsets outside the shadow (e.g. 1σ of the ephemeris). The blue
  dots are time marks along the centre line, every `marks` s.
* **Night side:** on a globe (orthographic projection, the default of `occ_map`) the
  Earth's limb is drawn as a circle, the night side is shaded transparently (with the
  twilight, Sun 0–18° below the horizon, as a gradient) and the terminator is drawn
  dash-dotted. `night` and `terminator` switch them on or off for any projection.
* **Stations:** positive stations are green circles and negative ones red crosses, numbered
  as in the legend (`observer_labels=:name` writes the names on the map).
* `geomap` adds the shaded shadow band, country borders and capitals.
* **Prediction and post-fit path:**
  * `occ_map(event; error_mas=20)` is the prediction from the ephemeris alone, with the
    nominal radius and the cross-track ephemeris uncertainty (mas) as ±1σ lines;
  * `occ_map(fit)` is the reconstruction after the event: the fitted centre, the limits as
    the limb's extent across the motion, and their 1σ from the fit;
  * `geomap(prediction, fit)` draws both, so the shift of the shadow relative to the
    prediction (the astrometric offset) can be seen at a glance. `overlays`,
    `overlay_labels`, `overlay_colors` and `path_label` add further paths;
  * the post-fit 1σ lines are much closer to the limits than those of a prediction: after
    the event the chords fix the shadow to kilometres.

![Huya 2019 prediction and post-fit path](../assets/figures/huya_geomap_comparison.png)

![Huya 2019 on the globe](../assets/figures/huya_globe.png)

## Mesh views

`meshview(shape)` shows a shape model in its body frame with the x, y and z
axes (z = spin axis). `meshview(body, t)` shows the body as seen from the Earth at time t,
with its projected limb and the position angle of the pole.
