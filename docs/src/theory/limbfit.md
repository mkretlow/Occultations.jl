# The mathematics of the limb fit

All fits minimise χ² = Σ rᵢ² over normalised residuals rᵢ. The sections below define these
residuals and how the minimum, the uncertainties and the posterior are obtained.

## 1. Geometry: the fundamental plane and the tracks

For an observer and a time t, the
fundamental-plane coordinates (f, g) [km] give the position of the star relative to the body's
ephemeris centre. They are measured perpendicular to the line of sight at the body's
distance, f towards east and g towards north. The full reduction (topocentric observer,
light time, aberration, proper motion and parallax of the star) is done once per chord. The
resulting track (f, g)(t) is stored as Chebyshev polynomials accurate to < 1 m. A fitted
clock offset Δtᵢ of chord i shifts its times, t → t + Δtᵢ.

## 2. The limb

In event k the body's centre is offset from the ephemeris by
cₖ = (f₀, g₀)ₖ, the astrometric result. The limb is a closed outline around cₖ, given by its
radius r(θ) in each direction:

* **Ellipse** (`EllipseLimb`): semi-major axis a, oblateness ε = 1 − b/a, and the position
  angle φ of the minor axis (from north towards east). With (u, w) the coordinates of a
  direction along the major and minor axes, r(θ) = [(u/a)² + (w/b)²]^(−1/2), and the area-equivalent radius is √(ab).
* **Circle** (`CircleLimb`): r = R.
* **Triaxial ellipsoid** (`EllipsoidLimb`): semi-axes a, b, c in the body frame, oriented by
  the matrix M (ICRS → body) of the pole (α, δ) and the prime-meridian angle W.
  * W(t) = W0 + 360°·(t_emit − t_ref)/P, at the light-time corrected emission time of each
    event.
  * The silhouette is exact. P is a sky basis (east, north, line of sight). The quadratic
    form Q = Pᵀ Mᵀ diag(a⁻², b⁻², c⁻²) M P is reduced to the sky plane by the Schur
    complement S = Q₁₂,₁₂ − q qᵀ/Q₃₃, which eliminates the line-of-sight coordinate. The
    apparent ellipse is xᵀ S x = 1, with semi-axes 1/√λ for the eigenvalues λ of S.
* **Mesh** (`MeshLimb`): the projected vertices are rotated to the event, and each projected
  edge raises the radius r(θ) in the directions it spans. The result is a polar outline
  sampled in n directions (the mesh must be star-shaped as seen from its centre); it is
  scaled by `scale` and rotated by `rotation`.

The radial signed distance of a point x (relative to cₖ) from the limb is
h(x) = |x| − r(arg x): negative inside, positive outside.

## 3. Contact residuals

A contact of chord i at the observed time t (error σₜ) lies at
xᵢ = (f, g)ᵢ(t + Δtᵢ) − cₖ, where the track moves with velocity v. With a limb allowance
σ_limb, its uncertainty is σ = √((σₜ |v|)² + σ_limb²).

* **Along the chord** (`residual=:chord`): s is the solution of x + s·v on the limb, the
  entry root for an immersion and the exit root for an emersion.
  * Ellipse: a quadratic in s, A s² + B s + C = 0, solved analytically.
  * Other outlines: bracketing of h along the line, then bisection.
  * The residual is r = −s|v|/σ, i.e. the timing residual × |v|/σ. If the line misses the
    limb, s is taken at the closest approach and a second residual (distance of the miss)/σ
    pulls the limb towards the chord.
* **Radial** (`residual=:radial`): r = h(x)/σ. The difference from the residual along the
  chord is explained in [Residuals along the chords or radial residuals?](residuals.md).

## 4. Negative chords

Along the observed track of a negative chord, d = min h(x) is the
closest approach of the track to the limb; d < 0 means the limb crosses the track.
* `:hard` (the default): parameters with d < 0 are inadmissible (χ² = ∞, zero posterior).
* `:soft`: residual max(−d, 0)/(sensitivity·|v|), where sensitivity is the exposure plus
  dead time of the chord.

## 5. Further constraints

Each is a residual (model − observation)/σ:
* rotational amplitudes: Δm = 2.5 log(a/b) − 1.25 log[(a² cos²ψ + c² sin²ψ)/(b² cos²ψ + c² sin²ψ)],
  for the aspect ψ of the pole at the epoch of the light curve;
* rotational phases: the phase since the last maximum of the projected area at the event,
  wrapped to (−½, ½] of the light-curve period (half a rotation for an ellipsoid);
* Gaussian priors: (x − μ)/σ. Other priors (`LogPrior`, e.g. `isotropic_pole`) enter only
  the posterior.

## 6. The minimum: Levenberg–Marquardt with multiple starts


* **Step:** δ from (JᵀJ + λ D) δ = −Jᵀr, with D = diag(JᵀJ). λ starts at 10⁻³ and is
  divided by 5 after a successful step, multiplied by 10 after a failed one.
* **Convergence:** relative χ² change < 10⁻⁹, or every |δⱼ| < 10⁻⁷ × the parameter's scale.
* **Jacobian:** J = ∂r/∂x by forward-mode automatic differentiation; otherwise central finite
  differences with steps of 10⁻⁶ × scale (joint light-curve, atmosphere and mesh fits).
* **Starts:** the first start is the initial value; the others are drawn uniformly within
  ±2 × scale around it. An inadmissible start is replaced by a random search nearby. The
  best χ² wins.
* **Canonical form:** the solution is brought to a canonical form, then wrapped:

  | Model | Canonical form |
  |---|---|
  | Ellipse | ε ≥ 0 (otherwise the axes are swapped and φ + 90°), φ mod 180° |
  | Ellipsoid | a ≥ b (otherwise swapped and W + 90°), W0 mod 180° |
  | Periodic parameters | wrapped |

  It is then polished by a final Levenberg–Marquardt run.

## 7. Uncertainties of the fit


* **Covariance:** C = (JᵀJ)⁻¹ at the minimum (pseudo-inverse if singular).
  `scale_errors = true` multiplies it by χ²/dof.
* **Derived quantities** (b, R_eq, D_eq, albedo inputs, astrometric offsets, apparent axes
  of an ellipsoid, …) are propagated linearly: C_d = G C Gᵀ, with G = ∂d/∂x by central
  differences.
* **Limb allowance:** `σ_limb = :auto` chooses it by bisection so that χ²/dof = 1.
* **Likelihood profile:** `profile(fit, name)` refits all other parameters on a grid of fixed
  values. Δχ² = 1 gives the 1σ interval, which can be asymmetric.

## 8. The posterior


* **Density:** log p(x) = −χ²(x)/2 + Σ log prior(xⱼ), and −∞ for inadmissible x (hard
  constraints, bounds).
* **Sampler:** the affine-invariant ensemble sampler of [GoodmanWeare2010](@citet), as in
  emcee [ForemanMackey2013](@citep):
  * 4n walkers (at least 16; n free parameters) start in a small Gaussian ball around the
    least-squares solution;
  * each walker moves by the stretch move y = xⱼ + z (x − xⱼ) towards a random other walker
    xⱼ, with z drawn from g(z) ∝ 1/√z on [1/a, a] (a = 2);
  * the move is accepted with probability min(1, zⁿ⁻¹ p(y)/p(x)).
* **Output:** the first third of the steps is burn-in. The integrated autocorrelation time
  τ (the automatic window of [Sokal1997](@citet), c = 5) shows how many independent samples the chain holds. Summaries
  are the median and the 15.87–84.13 % quantiles.

## 9. Several events


* The shape parameters are shared; each event has its own centre cₖ. `per_event` makes
  individual shape parameters event-specific, and per-chord clock offsets Δtᵢ can be added.
* For an ellipsoid or a mesh, each event has its own orientation M(W(tₖ)) and line of sight.
  Events at different rotational phases and aspects therefore constrain the 3-D shape
  together.

## 10. Light curves instead of contacts

With `photometry = true`, the flux of each light curve is modelled directly:
* **Model:** F(t) = f_b + (f_t − f_b) I(t), where I(t) is the normalised intensity at the
  observer's position on the track.
* **Airless body:** I comes from Fresnel diffraction at the limb, at the predicted contacts
  and velocities normal to the limb. It is integrated over the passband, the stellar disc
  and the exposure time.
* **Atmosphere:** I is the ray-traced refractive flux (see [Atmospheres](../guide/atmospheres.md)), multiplied by the
  diffraction of the surface edge.
* **Flux levels:** f_t and f_b are the linear least-squares solution for given I, so they are
  profiled out.
* **Residuals:** (F_obs − F)/σ_F, one per data point.
