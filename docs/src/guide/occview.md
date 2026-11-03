# The occultation in 3-D

`occview(fit)` and `occview(posterior)` show a 3-D limb fit (`EllipsoidLimb`, `MeshLimb`) in
the geometry of the occultation. Everything comes from the fitted problem: shape, pole,
rotation (light-time corrected), chords and χ². The figure has four parts:

1. **The body at the event, as seen from the Earth** (top left; with GLMakie it can be
   rotated with the mouse):
   * the shape oriented for the event, lit by the Sun, with the night side dark;
   * the sunlit side (light) and the night side (dark): near opposition, i.e. at a small
     phase angle, the whole visible side is lit;
   * the pole (+z of the body, blue), the Sun direction (orange, an arrow from the centre;
     near opposition it points almost at the viewer), the sub-observer point (cyan) and the
     sub-solar point (yellow); the title gives the phase angle;
   * the positive chords (red) and the negative chords (gray, dashed) on the sky plane in
     front of the body.

   This shows where the chords cut the shape, which negative chord constrains which part
   of the limb, and where topography or a different orientation would be needed.
2. **The sky plane** (top right):
   * the silhouette (the fitted limb) with the chords;
   * for a posterior, `ndraw` silhouettes drawn from the samples (thin blue), which show the
     uncertainty of the 3-D solution as seen at the event.
3. **The rotational light curve** (bottom right): the magnitude from the projected area
   through one rotation, Δm = −2.5 log(A/Ā), for a geometric scatterer.
   * It is plotted against the **rotational phase** in [0, 1). Phase 0 is the last maximum
     of the light curve before the event, and the phase increases with time.
   * The event's phase is then the fraction of a rotation since that maximum, as in
     `RotationalPhase`. An ellipsoid has two equal maxima per rotation, so a published phase
     φ corresponds to φ mod 0.5 (2003 UY117: 0.63 → 0.13).
   * Markers show the event's phase (dashed) and the phase displayed (red); both values are
     also printed next to the sliders.
   * The title compares the model amplitude with the observed ones (`RotationalAmplitude`).
   * The curve is the model's prediction for the light curve: its shape and the phase of
     its maxima at the event (`RotationalPhase`).
4. **Sliders** (bottom; interactive with GLMakie):
   * *rotational phase* (0 = the maximum before the event; it starts at the event's phase):
     turns the body about its spin axis. All three panels follow: the 3-D view, the
     silhouette against the chords, and the marker on the light curve;
   * *event* (multi-event fits): steps through the events, each with its own orientation,
     chords and silhouette;
   * with `explore = true`, also the pole (RA, Dec), *W0* and the size (semi-axes *a*, *b*,
     *c*, or the *scale* of a mesh). W0 is the prime-meridian angle at the reference epoch
     (IAU convention), measured along the body's equator from its ascending node on the ICRS
     equator. With the rotation period P it gives the orientation at every time,
     W(t) = W0 + 360° (t − t_ref)/P. The live χ² of the problem (all events, all constraints) is shown
     next to them. This is meant for finding starting values of a mesh fit, or for seeing
     how well determined the pole is.
   * the *reset* button sets all sliders, and with them the parameters, back to their start
     values: the fitted solution, the event's rotational phase and the first event shown.

![2003 UY117 in 3-D](../assets/figures/2003uy117_occview.png)

`occanimate(fit_or_posterior, "rotation.gif"; frames=120, framerate=24)` records the figure
through one rotation (GIF, MP4, WebM).

The views are static with CairoMakie (`CairoMakie.save("view.png", occview(fit))`) and
interactive with GLMakie (`using GLMakie; display(occview(fit; explore=true))`).
