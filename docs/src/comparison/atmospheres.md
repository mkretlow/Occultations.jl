# Atmosphere models: Baum & Code, Sicardy, Arimatsu and Occultations.jl

This page compares four descriptions of a stellar occultation by an atmosphere:

* the analytic model of [BaumCode1953](@citet), derived for the occultation of σ Arietis by
  Jupiter;
* the formalism and the ray-tracing code of B. Sicardy's group, summarised in the review of
  [Sicardy2022](@citet) and used for Pluto ([DiasOliveira2015](@citet), "DO15") and Triton
  ([MarquesOliveira2022](@citet));
* the model of [Arimatsu2026](@citet) for the tenuous atmosphere of (612533) 2002 XV93;
* the implementation in Occultations.jl ([`Atmosphere`](@ref), [`atmosphere_profile`](@ref),
  [`atmosphere_flux`](@ref) and the joint light-curve fit of [`LimbProblem`](@ref)).

The descriptions of the other models follow the papers; equation numbers refer to them.
Occultations.jl is described at the state of this manual. Section [Validation](#Validation)
lists the numerical comparisons. All of them can be rerun from the tests.

## Common picture

All four models describe the same physics. A transparent atmosphere bends the starlight by a
small angle ω (negative: towards the body). Deeper rays are bent more, so a beam of rays
spreads out behind the body and the star fades gradually instead of disappearing abruptly.
The quantities are:

* the refractivity ν = n − 1 = K n_g (molecular refractivity K, number density n_g);
* the bending angle ω(r) of a ray whose closest approach to the centre is r;
* the position of the ray in the shadow plane at distance D, z = r + Dω (Sicardy's Eq. 8);
* the normalised flux φ. Energy conservation gives φ = dr/dz (Eq. 9) for a straight limb.

The models differ in how the density profile is built, how ω is computed, and which
corrections are applied to φ: limb curvature, the far-limb image, the stellar disc,
diffraction, extinction.

## Baum & Code (1953)

**Assumptions.** The scale height H is constant and small compared with the radius, so the
refractivity decreases exponentially, ν(r) = ν₀ e^(−(r−r₀)/H). Each ray is a straight line.
The limb is straight: there is no focusing by its curvature and no far-limb image. The star is
a point.

**Steps** (as derived in [Sicardy2022](@citet), Sect. 5):

1. Bending angle of an exponential atmosphere (Eq. 18): ω(r) ≈ −ν(r) √(2πr/H). Most of the
   bending happens along a path of length √(2πrH) around the closest approach (Eq. 19).
2. Since dω/dr ≈ −ω/H, the flux is 1/φ = 1 + D dω/dr = 1 − Dω/H (Eqs. 20–21).
3. The flux has dropped to one half where the ray is bent by one scale height:
   ω₁/₂ = −H/D, so z₁/₂ = r₁/₂ − H (Eq. 22). This happens at the refractivity
   ν₁/₂ = √(H³/(2π r₁/₂ D²)) (Eq. 23).
4. Eliminating r between 1/φ − 1 = e^(−(r−r₁/₂)/H) and z = r + Dω gives the **Baum & Code
   equation** (Eq. 27):

   (1/φ − 2) + ln(1/φ − 1) = −Δz/H,  Δz = z − z₁/₂.

**Consequences.**

* Above the half-light level the flux returns to 1 exponentially, φ ≈ 1 − e^(−Δz/H)
  (Eq. 29). With 1 % photometry the drop is detectable only up to about 4.6 H above the
  half-light level.
* Deep in the shadow the flux falls slowly, φ ≈ H/|Δz| (Eq. 30).
* Differentiating Eq. 27 at φ = 1/2 gives dφ/dz = 1/(8H). The slope of the light curve at
  half light, divided by the shadow velocity, gives the scale height. With an assumed
  temperature it then gives the molecular mass, which is how Baum & Code estimated the mean
  molecular mass of Jupiter's upper atmosphere.
* The pressure at half light follows from Eq. 23 and the ideal gas law (Eq. 31):

  P₁/₂ ≈ (GMμ/(KD)) √(H⁵/(2π r₁/₂⁵)).

  For bodies from Jupiter to Pluto this gives a fraction of a pascal (a few µbar) in
  Sicardy's Table 1. Ground-based occultations always probe roughly the same pressure range.

**Limits.** The model needs a constant H with H ≪ r, and it has no surface. It ignores the limb
curvature (r/|z| ≈ 1 near half light) and the far-limb image, which matter only near the
shadow centre. It cannot describe a central flash, an inversion layer, or an atmosphere so
tenuous that the surface cuts the light curve before half light (2002 XV93).

## Sicardy (2022) and the ray-tracing code of his group

[Sicardy2022](@citet) derives the general formalism and shows that Baum & Code is its limit for
an exponential atmosphere.

**Forward model** (the direct approach of DO15 and [MarquesOliveira2022](@citet)):

1. **Temperature.** A prescribed profile T(r): Eq. 13 needs one boundary value of T, which the
   occultation does not provide. DO15 use a parametric template, a conic from the surface
   through a temperature maximum to an inflection point, a cubic above it, then an isothermal
   branch. For Triton 2017 the upper branch has a constant gradient (0.105 K/km), and the
   lower part was tuned to fit the central flashes (Table B.1).
2. **Density.** Hydrostatic equilibrium, dp/dr = −μ n_g g(r), with p = n_g k_B T (Eqs. 12–13).
   The boundary condition is the pressure at a reference radius: p₁₂₇₅ for Pluto in DO15,
   p₁₄₀₀ for Triton.
3. **Bending angle.** Exactly (Bouguer's rule, Eq. 6), or in the straight-line approximation
   ω(r) = ∫ (r/ρ) dν/dρ dl (Eq. 16). The latter is adequate because ω ≲ 10⁻⁵ rad for
   ground-based observers: the ray is displaced by less than 0.2 m inside the atmosphere (DO15,
   App. A.2).
4. **Mapping** z = r + Dω (Eq. 8).
5. **Flux of one image** (Eq. 32):

   φ = [1/(1 + D dω/dr)] · (r/|z|) · e^(−τ).

   The factor r/|z| is the focusing by the limb curvature, and e^(−τ) the extinction by haze.
   The fluxes of the primary (near-limb) and secondary (far-limb) images are summed.
6. **Stellar disc.** If r⋆ ≪ H, Eq. 32 holds for a point source. Otherwise, and always near the
   shadow centre, the flux is the area of the refracted stellar image divided by πr⋆²
   (Clausius' theorem), computed numerically in two dimensions.
7. **Central flash.** Near z = 0 the two images merge into an Einstein–Chwolson ring of width
   2r⋆φ_c (Eq. 33). The flash height is φ_cf = 4 (r_cf/r⋆) φ_c ≈ 4H/r⋆ (Eqs. 34–35), where r_cf
   is the radius of the flash layer and φ_c the flux there without the focusing term. A
   non-spherical atmosphere gives up to four images and caustics (Sect. 6.3).

**Inverse approach.** An Abel inversion of the light curve (Eqs. 7, 10) gives ν(r), and from it
n(r), p(r) and T(r). It assumes a single (near-limb) image, which fails below about 20 km on
Triton, where the secondary image contributes ([MarquesOliveira2022](@citet), App. B.3).

**Fitting practice for Triton 2017** ([MarquesOliveira2022](@citet), Sect. 4.2):

* 52 light curves are fitted simultaneously with the fixed template.
* Free parameters: p₁₄₀₀, the cross-track offset Δρ of the shadow centre, and Triton's
  fraction φ₀ of the total flux for every light curve, 54 in all, with 68 446 data points.
* χ² is mapped over a grid in (Δρ, p₁₄₀₀). The along-track offset is decoupled.
* The central flashes are first excluded, then added to test the shape of the deep atmosphere.
* Result: p₁₄₀₀ = 1.18 ± 0.03 µbar, χ²/dof = 0.85.

**Modified Baum & Code equation** (App., Eqs. 36–41). For a source at a finite distance D′
behind the body (a transit), a geometric factor g = D′/(D + D′) appears:
g(1/φ − 1) + ln(1/φ − 1) = −(y_i − r₁/₂)/H. It also describes "deep images", which probe only
a few scale heights below r₁/₂ because of the logarithm in Eq. 41.

## Arimatsu et al. (2026): 2002 XV93

[Arimatsu2026](@citet) (Methods, "Atmosphere refraction model fit") use the same forward model
for a body of about 250 km:

1. **Temperature.** A Pluto-like profile: 47 K at the surface (thermal model), rising by
   2.7 K/km up to about 20 km.
2. **Density.** Hydrostatic equation integrated upward from the surface (their Eq. 3), with a
   bulk density of 1500 kg/m³. The Jeans parameter is λ ≈ 1, so above about 25 km the density
   follows the free-molecular-flow law n ∝ (R + z)⁻². They note that the upper profile hardly
   matters if it is connected smoothly.
3. **Refractivity** ν = Kn, for pure CH₄ or N₂.
4. **Bending angle.** The straight-line integral (their Eq. 4, the same as Sicardy's Eq. 16).
5. **Mapping** y = r + ωD with D = 37.0 au.
6. **Flux.** φ = (dr/dy)(r/y) (their Eq. 5): Sicardy's Eq. 32 with τ = 0. The second factor is
   the limb-curvature focusing.
7. **Fit.** R, the surface pressure p and the centre are fitted to the Kiso light curve (χ²,
   their Eq. 6). For Fukushima, which saw no disappearance, only the distance to the shadow
   centre is fitted, with the Kiso atmosphere fixed. Result: CH₄ p = 124 (+15/−12) nbar.

Baum & Code does not apply here. Above the inversion the hydrostatic scale height (about 600 km)
is larger than the radius, and the free-molecular-flow tail, not a hydrostatic layer, does
most of the bending. The surface cuts the light curve at a flux of about 0.46, before the
half-light level would be reached.

## Occultations.jl

The implementation follows the forward model of Sicardy and Arimatsu, step by step:

| Step | Implementation |
|---|---|
| Temperature | `Atmosphere(gas; T0, dTdz, z_inv, z_fmf, width)`: T0 at the surface, a gradient up to `z_inv`, isothermal above; corners smoothed over `width`. Or any `temperature = T(r)`, e.g. a [`DO15Temperature`](@ref) template (conic, cubic, straight line, as DO15 Appendix A.1). [`triton_atmosphere`](@ref) is Table B.1 of Marques Oliveira et al. |
| Density | d ln n/dz = −μg/(k_B T) − (dT/dz)/T (Sicardy Eq. 13), trapezoidal rule. The step is at most 1/50 of the smallest density scale height. The boundary value is the pressure at the surface, or at `rref` (e.g. p₁₄₀₀). |
| Free molecular flow | Above `z_fmf`, d ln n/dz blends smoothly into −2/(R + z) (Arimatsu). The tail ν ∝ ρ⁻² is integrated analytically to infinity, and rays above the table have ω = −πA/r² exactly. |
| Refractivity | ν = K n; for N₂, K = 1.091×10⁻²⁹ + 6.282×10⁻³²/λ² m³ (Washburn, as DO15), at the central wavelength of the light curve. |
| Bending angle | Straight-line integral (Sicardy Eq. 16) written as ω(r) = 4r ∫₀^∞ ν′(r+u²)/√(2r+u²) du. A sinh mapping of u resolves the deepest scale height. |
| Mapping | y = r + ωD, one ray tracing per event and passband, rescaled exactly to each observer's distance. |
| Flux | φ = (r/\|y\|)/\|dy/dr\|, summed over the near (+ρ) and far (−ρ) limb: Sicardy Eq. 32 with τ = 0, Arimatsu Eq. 5 with the far-limb image. |
| Stellar disc | Uniform disc. Within 20 stellar radii of the shadow centre, the exact two-dimensional mean of the circularly symmetric flux. For φ ∝ 1/ρ this gives φ(0) = 4 r_cf φ_c / r⋆, which is Sicardy's Eq. 34. Farther out, strips across the limb. |
| Surface | The ray grazing the surface defines the shadow edge. The refracted flux is multiplied by the Fresnel diffraction pattern of that edge, integrated over the passband and the stellar disc (not in the other models). |
| Exposure | The model is averaged over each exposure, finely near the surface edge and the central flash. |
| Fit | All light curves jointly with the limb ([`LimbProblem`](@ref) with `photometry=true, atmosphere=…`): the centre of each event, the radius (or fixed), the pressure. Per light curve the unocculted level is profiled, and with `background=:free` also the zero level (Triton's φ₀). Covariance, χ² profiles ([`atmosphere_test`](@ref)) and posterior sampling. |

## Side by side

| | Baum & Code | Sicardy / DO15 | Arimatsu et al. | Occultations.jl |
|---|---|---|---|---|
| Temperature | constant H | template T(r) | linear rise, isothermal | linear rise or any T(r), DO15 template |
| Upper atmosphere | exponential | hydrostatic | free molecular flow | hydrostatic or free molecular flow |
| Pressure boundary | – (P₁/₂ estimate) | p at a reference radius | surface | surface or reference radius |
| Bending | analytic, exponential | exact or straight line | straight line | straight line, analytic tail |
| Limb curvature (r/\|z\|) | no | yes | yes | yes |
| Far-limb image | no | yes | – (not stated) | yes |
| Stellar disc | point | 2-D image area | – (not stated) | 1-D strips; exact 2-D near the centre |
| Central flash | no | yes, also non-spherical | not relevant | spherical only |
| Haze e^(−τ) | no | yes | no | no |
| Surface edge | no surface | rays lost (geometric) | rays lost | Fresnel diffraction × refraction |
| Fitted | H from the slope | p_ref, Δρ, φ₀ per light curve | R, p, centre (Kiso) | centre, R, p, levels (and φ₀); covariance, profiles, posterior |

## Validation

Model tests, all in `test/runtests.jl` ("atmosphere: temperature templates, Baum & Code,
Triton"):

* **Baum & Code.** For an atmosphere with constant pressure scale height (T ∝ r⁻² on a
  Neptune-sized body, near-limb flux without focusing), the ray-traced flux agrees with the
  solution of Eq. 27 for Δz from −10H to +5H:

  | H | largest \|φ_ray − φ_BC\| |
  |---|---|
  | 20 km | 1.2 × 10⁻³ |
  | 50 km | 2.1 × 10⁻⁴ |
  | 100 km | 4.9 × 10⁻⁴ |

  The residual differences are of order H/r.
* **Half-light pressure** (Sicardy Eq. 31). The model pressure at the ray-traced half-light
  radius agrees with Eq. 31 within 0.7 % (H = 20, 50, 100 km: 0.0471, 0.465, 2.62 Pa against
  0.0474, 0.466, 2.64 Pa).
* **DO15 template.** With the Pluto parameters of DO15 (Table 3), the conic implies
  dT/dr(r₃) = −0.2057 K/km (DO15: −0.206). The Triton template passes through
  T(1453 km) = 54.55 K (Marques Oliveira et al.: 54.5 K).
* **Triton 2017** (p₁₄₀₀ = 1.18 µbar, the template of Table B.1):

  | quantity | Occultations.jl | Marques Oliveira et al. 2022 |
  |---|---|---|
  | p_surf / p₁₄₀₀ | 12.00 | 12.0 |
  | La Palma at 685 km from the centre: primary flux | 0.0614 | 0.061 |
  | La Palma: secondary flux | 0.0069 | 0.0072 |
  | deepest radius probed at La Palma | 1373.8 km | ≈ 1375 km |
  | altitude of the flash layer | 7.7 km | ≈ 8 km |
  | half-light distance from the shadow centre | 1419 km | "radius ≈ 1415 km" |

  The paper's half-light "radius" agrees with our distance in the shadow plane, and its
  0.55 µbar is the pressure at r ≈ 1415 km. The half-light ray itself passes about one scale
  height higher (r ≈ 1441 km, 0.17 µbar), as z₁/₂ = r₁/₂ − H implies (our reading).
* **Central flash.** For the 0.65 km stellar diameter of 2017, φ(0) = 184.4, equal to Eq. 34
  evaluated with the model's flash layer.
* **2002 XV93** (CH₄, 124 nbar, R = 244 km, against the published profile of Arimatsu et al.
  on Zenodo): shadow edge of the surface 229.9 km (230.1 km); flux at the edge 0.459 (0.465);
  rms difference 0.01 over 200–320 km. The joint fit of the published light curves gives
  p = 129 ± 12 nbar for CH₄ (124 +15/−12). See the example
  [2002 XV93, 2024-01-10](../generated/2002xv93_2024.md).

The Triton light curves themselves (CDS, J/A+A/659/A136) have not been fitted yet.

## Differences and limitations of Occultations.jl

* **No extinction term** e^(−τ): haze or dust cannot be fitted (Sicardy Eq. 32).
* **No diffraction of the central flash.** The flash is finite only through the stellar disc.
  [SicardyDettwiller2026](@citet) treat diffraction and interference at the flash; this matters
  when the stellar radius is smaller than the Fresnel scale (Triton 2017: r⋆ = 0.33 km,
  F ≈ 1.2 km).
* **Spherical atmospheres only:** no oblate flash layers, four images or caustics.
* **No Abel inversion:** the temperature profile is prescribed, not retrieved.
* **Template fixed in radius:** use [`DO15Temperature`](@ref) with a fixed body radius.
* **Surface edge:** diffraction × refraction is an approximation, valid when the refracted
  flux varies slowly over a Fresnel scale.
* **Refractivity** is evaluated at the central wavelength. For N₂ the band average of K over a
  broad visible band (0.45–0.85 µm) differs from that value by less than 0.1 %.
* **Speed:** about 0.05 s per profile (Triton template) and one profile per passband and
  pressure. The Jacobian with respect to the pressure uses finite differences.
