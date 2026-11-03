# Residuals along the chords or radial residuals?

Both fit the same limb model to the same contact times. They differ in how a contact's misfit
is measured and weighted.

**Along the chords** (`residual=:chord`, the default) uses the timing measurement itself. For
each contact the model predicts the time at which the station's chord crosses the fitted
limb.
* The residual is the distance along the chord between the observed and the predicted
  contact: the timing residual × the shadow speed.
* Its uncertainty is σ_t·v, plus `σ_limb` if set.

This is the actual measurement model: a timing error moves the contact point along the
chord, and nothing else. So χ² is the correct likelihood of the timings, and the errors and
posteriors follow from it. Special cases are handled explicitly:
* a chord that misses the fitted limb has no predicted contact and is penalised by the
  distance by which it misses;
* a grazing chord, whose timing hardly constrains the size, is weighted correctly, and its
  residual stays finite at tangency.

**Radial residuals** (`residual=:radial`, as in SORA and most papers) convert the observed time
into a point in the sky plane.
* The residual is the distance of that point from the limb, measured along the line from the
  ellipse centre.
* It is divided by the same σ_t·v.

The error, however, lies along the chord, not along the radius. If the chord meets the limb at
an angle θ to the limb normal, a timing error ε moves the point radially by only about
v·ε·cos θ. Each normalised residual is therefore understated by the factor cos θ:
* central chords (cos θ ≈ 1): the two definitions agree;
* grazing contacts (cos θ → 0): their misfit almost vanishes, however wrong the timing, and
  they get much less weight.

Consequences:

* **χ² is smaller with radial residuals**, so a good radial χ² can hide timing misfits:

  | Example | radial | along the chords | published |
  |---|---|---|---|
  | 2003 UY117 | 0.41 | 2.10 | 0.43 (radial) |
  | Varda | 1.49 | 1.69 | 1.53 (radial) |
  | 1995 SM55 | 39 | 46 | 1.55 |

* **The parameters usually agree within their errors.** They differ where grazing or oblique
  contacts matter, because those carry different weights (1995 SM55: a′ = 102.9 vs 105.4 km).
* **Errors:** a radial fit has no correct likelihood behind it, so its covariance and posterior
  are only approximate.

Which to use:
* `:chord` for results and uncertainties;
* `:radial` to reproduce or compare with published fits; the examples print both.
