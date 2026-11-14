# ---------------------------------------------------------------------------
# Refraction by a thin, spherically symmetric atmosphere (Pluto, 2002 XV93, ...):
# hydrostatic density profile, ray tracing of the bending angle, flux in the shadow
# (near and far limb, limb-curvature focusing), as e.g. in Sicardy et al. and
# Arimatsu et al. (2026, arXiv:2605.02243)
# ---------------------------------------------------------------------------

# molecular mass [amu] and molecular refractivity K [m³] (ν = K n) in the visible
const _GASES = Dict(
    :N2  => (μ=28.0134, K=1.091e-29, disp=6.282e-32 / 1.091e-29),   # Washburn (1930): K + 6.282e-32 m³/λ²[µm]
    :CH4 => (μ=16.0425, K=1.637e-29, disp=0.0),
    :CO  => (μ=28.0101, K=1.258e-29, disp=0.0),
)

"""
    DO15Temperature(r1, T1, dTdr1, r2, T2, r3, T3, r4, T4, dTdr4=0)

Temperature template `T(r)` [K] versus radius `r` [km] of Dias-Oliveira et al. (2015,
Appendix A.1), as used for Pluto and Triton by the Sicardy group:

* `r1 ≤ r ≤ r3`: a conic `c₁r + c₂T + c₃rT + c₄r² + c₅T² = 1` through `(r1, T1)` with gradient
  `dTdr1`, through the maximum `(r2, T2)` (`dT/dr = 0`) and through `(r3, T3)`;
* `r3 ≤ r ≤ r4`: a cubic, continuous in `T` and `dT/dr` at `r3` and reaching `T4` with gradient
  `dTdr4` at `r4`;
* `r ≥ r4`: a straight line `T4 + dTdr4 (r − r4)` (isothermal for `dTdr4 = 0`, as in DO15;
  Marques Oliveira et al. 2022 use 0.105 K/km for Triton).

Below `r1` the temperature is held at `T1`. Call it with a radius: `tp(1400.0)`.
"""
struct DO15Temperature
    r1::Float64; T1::Float64; dTdr1::Float64
    r2::Float64; T2::Float64
    r3::Float64; T3::Float64
    r4::Float64; T4::Float64; dTdr4::Float64
    c::NTuple{5,Float64}        # conic coefficients in the scaled variables x = r/r1, y = T/T1
    sgn::Float64                # branch of the quadratic in y
    dTdr3::Float64
    cubic::NTuple{4,Float64}    # T = a + b t + c t² + d t³, t = r − r3
end

function DO15Temperature(r1::Real, T1::Real, dTdr1::Real, r2::Real, T2::Real, r3::Real, T3::Real, r4::Real, T4::Real,
                         dTdr4::Real=0.0)
    r1 < r2 < r3 < r4 || throw(ArgumentError("the radii must increase: r1 < r2 < r3 < r4"))
    xs, ys = r1, T1                                    # scaling keeps the form of the conic
    pts = ((r1 / xs, T1 / ys), (r2 / xs, T2 / ys), (r3 / xs, T3 / ys))
    A = zeros(5, 5); b = zeros(5)
    for (i, (x, y)) in enumerate(pts)
        A[i, :] = [x, y, x * y, x^2, y^2]; b[i] = 1
    end
    # dF/dx + dF/dy y' = 0 at point 1 (gradient dTdr1) and point 2 (maximum)
    for (i, (x, y), yp) in ((4, pts[1], dTdr1 * xs / ys), (5, pts[2], 0.0))
        A[i, :] = [1, yp, y + x * yp, 2x, 2y * yp]
    end
    c = Tuple(A \ b)
    # the root of the quadratic in y that passes through point 1
    sgn = abs(_conic_y(c, pts[1][1], 1.0) - pts[1][2]) <= abs(_conic_y(c, pts[1][1], -1.0) - pts[1][2]) ? 1.0 : -1.0
    for (x, y) in pts
        isapprox(_conic_y(c, x, sgn), y; rtol=1e-8) ||
            throw(ArgumentError("no single conic branch passes through points 1-3; check the template parameters"))
    end
    x3, y3 = pts[3]
    dTdr3 = -(c[1] + c[3] * y3 + 2c[4] * x3) / (c[2] + c[3] * x3 + 2c[5] * y3) * ys / xs
    # cubic Hermite from (r3, T3, dTdr3) to (r4, T4, dTdr4)
    h = r4 - r3
    cc = (3(T4 - T3) / h - 2dTdr3 - dTdr4) / h
    dd = (dTdr3 + dTdr4 - 2(T4 - T3) / h) / h^2
    return DO15Temperature(r1, T1, dTdr1, r2, T2, r3, T3, r4, T4, dTdr4, c, sgn, dTdr3, (T3, dTdr3, cc, dd))
end

function _conic_y(c, x, sgn)
    qa, qb, qc = c[5], c[2] + c[3] * x, c[1] * x + c[4] * x^2 - 1
    abs(qa) < 1e-14 * abs(qb) && return -qc / qb
    return (-qb + sgn * sqrt(max(qb^2 - 4qa * qc, zero(qb)))) / (2qa)
end

function (tp::DO15Temperature)(r)
    r <= tp.r1 && return tp.T1 + zero(r)
    r <= tp.r3 && return tp.T1 * _conic_y(tp.c, r / tp.r1, tp.sgn)
    if r <= tp.r4
        t = r - tp.r3
        a, b, c, d = tp.cubic
        return a + t * (b + t * (c + t * d))
    end
    return tp.T4 + tp.dTdr4 * (r - tp.r4)
end

Base.show(io::IO, tp::DO15Temperature) =
    @printf(io, "DO15Temperature(%.1f km %.1f K, max %.1f K at %.1f km, %.1f K at %.1f km, %.1f K at %.1f km, %.3f K/km above)",
            tp.r1, tp.T1, tp.T2, tp.r2, tp.T3, tp.r3, tp.T4, tp.r4, tp.dTdr4)

"""
    Atmosphere(gas=:N2; T0=40, dTdz=0, z_inv=0, z_fmf=Inf, μ, K, ztop=NaN, width=2,
               temperature=nothing, rref=NaN)

Spherically symmetric, transparent atmosphere in hydrostatic equilibrium.

* `gas` – `:N2`, `:CH4` or `:CO` (molecular mass `μ` [amu] and refractivity `K` [m³];
  both can be given explicitly)
* temperature: `T0` [K] at the surface, rising (or falling) by `dTdz` [K/km] up to the
  altitude `z_inv` [km], isothermal above (a Pluto-like inversion layer; `dTdz = 0` gives
  an isothermal atmosphere)
* or any temperature profile `temperature(r)` [K] versus the radius `r` [km] (a function or
  callable object, e.g. a [`DO15Temperature`](@ref) template); `T0`, `dTdz` and `z_inv` are
  then ignored. The profile is fixed in radius, so it suits fits with a fixed body radius.
* `z_fmf` [km]: above this altitude the density follows the free-molecular-flow profile
  `n ∝ (R + z)⁻²` (low Jeans parameter, Arimatsu et al. 2026)
* `ztop` [km]: top of the model (default: 30 scale heights above the inversion; with a
  `temperature` profile, the altitude at which the density has dropped by e⁻³⁰)
* `width` [km]: the top of the inversion layer and the onset of the free molecular flow are
  smooth transitions of this width (a kink in the density gradient would produce a
  spurious jump of the refracted flux, since φ ∝ 1/(1 + D dω/dr))
* `rref` [km]: the pressure refers to this radius (e.g. `p₁₄₀₀` for Triton) instead of the
  surface

The free parameter of a fit is the pressure (`pressure`, nbar) at the surface or at `rref`;
the radius comes from the limb model and the gravity from the body (`GM` or `density`). See
[`atmosphere_profile`](@ref), [`xv93_atmosphere`](@ref), [`pluto_atmosphere`](@ref),
[`triton_atmosphere`](@ref).
"""
struct Atmosphere
    gas::Symbol
    μ::Float64
    K::Float64
    T0::Float64
    dTdz::Float64
    z_inv::Float64
    z_fmf::Float64
    ztop::Float64
    width::Float64
    temperature::Any
    rref::Float64
end

function Atmosphere(gas::Symbol=:N2; T0::Real=40.0, dTdz::Real=0.0, z_inv::Real=0.0, z_fmf::Real=Inf,
                    μ=nothing, K=nothing, ztop::Real=NaN, width::Real=2.0, temperature=nothing, rref::Real=NaN)
    g = get(_GASES, gas, nothing)
    g === nothing && (μ === nothing || K === nothing) &&
        throw(ArgumentError("unknown gas $gas: give μ [amu] and K [m³]"))
    width > 0 || throw(ArgumentError("width must be positive"))
    return Atmosphere(gas, something(μ, g === nothing ? NaN : g.μ), something(K, g === nothing ? NaN : g.K),
                      T0, dTdz, z_inv, z_fmf, ztop, width, temperature, rref)
end

function Base.show(io::IO, a::Atmosphere)
    ref = isfinite(a.rref) ? @sprintf(", pressure at r = %.1f km", a.rref) : ""
    fmf = isfinite(a.z_fmf) ? @sprintf(", FMF above %.1f km", a.z_fmf) : ""
    if a.temperature === nothing
        @printf(io, "Atmosphere(%s, T0 = %.1f K, dT/dz = %.2f K/km up to %.1f km%s%s)", a.gas, a.T0, a.dTdz, a.z_inv, fmf, ref)
    else
        print(io, "Atmosphere(", a.gas, ", T(r) = ", a.temperature, fmf, ref, ")")
    end
end

"""
    xv93_atmosphere(gas=:CH4)

The model of Arimatsu et al. (2026) for (612533) 2002 XV93: surface temperature 47 K,
inversion layer with 2.7 K/km up to 20 km, free molecular flow above 25 km.
"""
xv93_atmosphere(gas::Symbol=:CH4) = Atmosphere(gas; T0=47.0, dTdz=2.7, z_inv=20.0, z_fmf=25.0)

"""
    pluto_atmosphere()

A simple Pluto-like N₂ atmosphere: 40 K at the surface, a strong inversion of 3 K/km up to
22 km (≈ 106 K), isothermal above.
"""
pluto_atmosphere() = Atmosphere(:N2; T0=40.0, dTdz=3.0, z_inv=22.0)

"""
    triton_atmosphere()

Triton's N₂ atmosphere of Marques Oliveira et al. (2022, Table B.1): a
[`DO15Temperature`](@ref) template from 38 K at the surface (r = 1353 km, 5 K/km) to a
maximum of 49.8 K at 1363 km, 48.2 K at 1375 km and 49.3 K at 1403 km, rising by 0.105 K/km
above. The pressure refers to r = 1400 km (`p₁₄₀₀`; they find 1.18 ± 0.03 µbar = 1180 nbar
in October 2017). Use it with a radius of 1353 km and GM = 1428 km³/s².
"""
triton_atmosphere() = Atmosphere(:N2; temperature=DO15Temperature(1353.0, 38.0, 5.0, 1363.0, 49.8, 1375.0, 48.2,
                                                                   1403.0, 49.3, 0.105), rref=1400.0)

_softplus(x) = x > 30 ? x : log1p(exp(x))
_logistic(x) = 1 / (1 + exp(-x))
# smooth min(z, z_inv): the temperature gradient switches off over `width`
_temperature(a::Atmosphere, z) = a.T0 + a.dTdz * (a.z_inv - a.width * _softplus((a.z_inv - z) / a.width) +
                                                  a.width * _softplus(a.z_inv / a.width) - a.z_inv)
_dTdz(a::Atmosphere, z) = a.dTdz * _logistic((a.z_inv - z) / a.width)
# temperature and its gradient at altitude z above a body of radius R
_temperature(a::Atmosphere, R, z) = a.temperature === nothing ? _temperature(a, z) : a.temperature(R + z)
_dTdz(a::Atmosphere, R, z) = a.temperature === nothing ? _dTdz(a, z) : ForwardDiff.derivative(a.temperature, R + z)

function _refractivity_K(a::Atmosphere, λ)
    g = get(_GASES, a.gas, nothing)
    return g === nothing || g.disp == 0 ? a.K : a.K * (1 + g.disp / λ^2)
end

"""
    AtmosphereProfile

Result of [`atmosphere_profile`](@ref): bending angle `ω(r)` [rad] versus closest-approach
radius `r` [km] of the rays and their position `y = r + ω D` in the shadow plane; `pressure`
is the given pressure (at the surface or at `rref`), `surface_pressure` [nbar] the pressure at
the surface and `scale_height` [km] the pressure scale height there.
[`atmosphere_flux`](@ref) evaluates the flux at any distance from the shadow centre.
"""
struct AtmosphereProfile
    R::Float64
    D::Float64
    pressure::Float64
    r::Vector{Float64}
    ω::Vector{Float64}
    y::Vector{Float64}
    dydr::Vector{Float64}
    y_surface::Float64
    star_radius::Float64
    edge::EdgeProfile
    scale_height::Float64
    surface_pressure::Float64
end

"""
    atmosphere_profile(atm, pressure, R, GM, D; λ=0.7, Δλ=0.3, star_diam=0, nz=3000, nr=1500, nu=300)
        -> AtmosphereProfile

Ray tracing through the atmosphere `atm` with surface `pressure` [nbar] on a body of radius
`R` [km] and `GM` [km³/s²], seen from distance `D` [km] (`pressure` refers to the surface,
or to the radius `atm.rref`):

1. number density from the hydrostatic equation `dp/dz = -μ n g(z)`, `p = n k T(z)`,
   `g = GM/(R+z)²`, blended smoothly into `n ∝ (R+z)⁻²` above the free-molecular-flow
   altitude (this tail is integrated analytically to infinity: it dominates the bending of
   tenuous atmospheres); the altitude step is at most 1/50 of the smallest density scale
   height (`nz` is the minimum number of steps);
2. refractivity `ν = K n`; bending angle of the ray with closest approach `r`,
   `ω(r) = ∫ (r/ρ) dν/dρ dx = 4r ∫₀^∞ ν′(r+u²)/√(2r+u²) du`;
3. the ray reaches the shadow plane at `y = r + ω(r) D`; flux `φ = (dr/dy)(r/|y|)` summed over
   the near (`y = ρ`) and far (`y = -ρ`) limb (the factor `r/|y|` is the focusing by the
   limb curvature); rays with `r < R` hit the surface.
"""
function atmosphere_profile(atm::Atmosphere, pressure::Real, R::Real, GM::Real, D::Real; λ::Real=0.7, Δλ::Real=0.3,
                            star_diam::Real=0.0, nz::Integer=3000, nr::Integer=1500, nu::Integer=300)
    μ = atm.μ * AMU
    K = _refractivity_K(atm, λ)
    gs = GM * 1e9 / (R * 1e3)^2                                          # m/s²
    T0 = _temperature(atm, R, 0.0)
    H0 = K_BOLTZMANN * T0 / (μ * gs) / 1e3
    fmf = isfinite(atm.z_fmf)
    # d ln n / dz [1/km]: hydrostatic, blended smoothly into the free-molecular-flow law -2/(R+z)
    dlog = _dlog_function(atm, atm.temperature, Float64(R), Float64(GM), μ)
    ztab = if fmf
        atm.z_fmf + 10atm.width
    elseif !isnan(atm.ztop)
        atm.ztop
    elseif atm.temperature === nothing
        Hinv = K_BOLTZMANN * _temperature(atm, atm.z_inv) / (μ * GM * 1e9 / ((R + atm.z_inv) * 1e3)^2) / 1e3
        atm.z_inv + 30Hinv
    else                                       # where the density has dropped by e⁻³⁰
        _efold_altitude(dlog, 30.0, 100R)
    end
    isfinite(atm.rref) && !(R <= atm.rref <= R + ztab) &&
        throw(ArgumentError("reference radius rref = $(atm.rref) km lies outside the model (R = $R km, top $(R + ztab) km)"))
    # step: a fraction of the smallest density scale height (and of the transition width)
    Hmin = minimum(1 / abs(dlog(zz)) for zz in range(0.0, ztab; length=2001))
    dzmax = min(Hmin / 50, atm.temperature === nothing ? atm.width / 10 : Inf)
    ng = clamp(ceil(Int, ztab / dzmax) + 1, nz, 2_000_000)
    z = collect(range(0.0, ztab; length=ng))
    dz = z[2] - z[1]
    lnn = _integrate_lnn(dlog, z)
    # boundary condition: the pressure [nbar → Pa] at the surface or at the reference radius
    zref = isfinite(atm.rref) ? atm.rref - R : 0.0
    kref = clamp(floor(Int, zref / dz) + 1, 1, ng - 1)
    wref = (zref - z[kref]) / dz
    nref = pressure * 1e-4 / (K_BOLTZMANN * _temperature(atm, R, zref))   # m⁻³
    lnn .+= log(max(nref, 1e-300)) - ((1 - wref) * lnn[kref] + wref * lnn[kref+1])
    ν = K .* exp.(lnn)
    dν = [(ν[min(i + 1, ng)] - ν[max(i - 1, 1)]) / (z[min(i + 1, ng)] - z[max(i - 1, 1)]) for i in 1:ng]
    ρf = R + ztab
    A = fmf ? ν[end] * ρf^2 : 0.0                                         # ν = A/ρ² above ρf
    νp_tab(ρ) = (x = (ρ - R) / dz; x >= ng - 1 ? dν[ng] : (k = max(floor(Int, x), 0); w = x - k; (1 - w) * dν[k+1] + w * dν[k+2]))
    # largest radius worth tracing: where |D dω/dr| < 1e-5
    # outer radius of the traced rays (independent of the pressure, so that the grid does not
    # move with it): far enough for the free-molecular-flow tail of any plausible atmosphere
    rmax = fmf ? max(R + 15R, ρf + 10H0) : ρf
    cu = sqrt(Hmin) / 2                                                   # √km
    rs = if fmf
        [R + (rmax - R) * x^2 for x in range(0.0, 1.0; length=nr)]
    else                                       # equal steps in the number of density e-folds (independent of the pressure)
        cum = cumsum([0.0; abs.(diff(lnn))])
        [R + _interp_lin(cum, z, s) for s in range(0.0, cum[end]; length=nr)]
    end
    ω = zeros(nr)
    Threads.@threads for j in 1:nr
        r = rs[j]
        if fmf && r >= ρf
            ω[j] = -π * A / r^2                                           # exact for ν ∝ ρ⁻²
            continue
        end
        uf = sqrt(max(ρf - r, 0.0))
        acc = 0.0
        if uf > 0                              # u = c sinh(s): steps ≪ √H near closest approach, also for a high top
            sf = asinh(uf / cu)
            h = sf / nu
            for i in 0:nu
                u = cu * sinh(i * h)
                acc += ((i == 0 || i == nu) ? 0.5 : 1.0) * νp_tab(r + u^2) / sqrt(2r + u^2) * cu * cosh(i * h) * h
            end
        end
        if fmf                                                             # analytic tail beyond ρf
            nt = 200
            for i in 0:nt-1
                t = (i + 0.5) / nt
                u = uf + t / (1 - t)
                ρ = r + u^2
                acc += (-2A / ρ^3) / sqrt(2r + u^2) / (1 - t)^2 / nt
            end
        end
        ω[j] = 4r * acc
    end
    y = rs .+ ω .* D
    dydr = [(y[min(j + 1, nr)] - y[max(j - 1, 1)]) / (rs[min(j + 1, nr)] - rs[max(j - 1, 1)]) for j in 1:nr]
    edge = _cached_profile(D, star_diam, λ, Δλ)
    psurf = exp(lnn[1]) * K_BOLTZMANN * T0 * 1e4                        # m⁻³ → Pa → nbar
    return AtmosphereProfile(R, D, pressure, rs, ω, y, dydr, y[1], star_diam / 2, edge, H0, psurf)
end

# d ln n / dz [1/km] at altitude z: hydrostatic, blended smoothly into the free-molecular-flow
# law -2/(R+z) (a closure specialised on the type of the temperature profile)
function _dlog_function(atm::Atmosphere, tf, R::Float64, GM::Float64, μ::Float64)
    fmf, zf, wd = isfinite(atm.z_fmf), atm.z_fmf, atm.width
    T(z) = tf === nothing ? _temperature(atm, z) : tf(R + z)
    dT(z) = tf === nothing ? _dTdz(atm, z) : ForwardDiff.derivative(tf, R + z)
    return function (zz)
        Tz = T(zz)
        g = GM * 1e9 / ((R + zz) * 1e3)^2
        hyd = -(μ * g / (K_BOLTZMANN * Tz)) * 1e3 - dT(zz) / Tz
        fmf || return hyd
        w = _logistic((zz - zf) / wd)
        return (1 - w) * hyd + w * (-2 / (R + zz))
    end
end

# altitude at which the density has dropped by e⁻ⁿ from the surface (at most zmax)
function _efold_altitude(dlog, n, zmax)
    z, acc = 0.0, 0.0
    while acc > -n && z < zmax
        h = 0.05 / abs(dlog(z))
        acc += dlog(z + h / 2) * h
        z += h
    end
    return z
end

function _integrate_lnn(dlog, z)
    nz = length(z)
    dl = Vector{Float64}(undef, nz)
    Threads.@threads for i in 1:nz
        dl[i] = dlog(z[i])
    end
    lnn = zeros(nz)
    for i in 2:nz
        lnn[i] = lnn[i-1] + (dl[i-1] + dl[i]) / 2 * (z[i] - z[i-1])
    end
    return lnn
end

# linear interpolation of y(x) on increasing x
function _interp_lin(x, y, xi)
    j = clamp(searchsortedlast(x, xi), 1, length(x) - 1)
    while j < length(x) - 1 && x[j+1] == x[j]
        j += 1
    end
    w = x[j+1] == x[j] ? 0.0 : (xi - x[j]) / (x[j+1] - x[j])
    return (1 - w) * y[j] + w * y[j+1]
end

# the same atmosphere seen from the distance D (y = r + ωD; exact) with the edge diffraction `edge`
# computed for that distance
function _rescale(p::AtmosphereProfile, D::Real, edge::EdgeProfile)
    D == p.D && edge === p.edge && return p
    y = p.r .+ p.ω .* D
    dydr = 1 .+ (p.dydr .- 1) .* (D / p.D)
    return AtmosphereProfile(p.R, D, p.pressure, p.r, p.ω, y, dydr, y[1], p.star_radius, edge, p.scale_height,
                             p.surface_pressure)
end

# refracted flux at distance ρ ≥ 0 from the shadow centre (near + far limb)
function _refracted_flux(p::AtmosphereProfile, ρ::Real)
    ρ = max(abs(ρ), 1e-6)
    ρ >= p.y[end] && return 1.0
    acc = 0.0
    n = length(p.y)
    for target in (ρ, -ρ)
        (target < p.y[1] || target > p.y[end]) && continue
        j = clamp(searchsortedlast(p.y, target), 1, n - 1)
        w = (target - p.y[j]) / (p.y[j+1] - p.y[j])
        r = (1 - w) * p.r[j] + w * p.r[j+1]
        dd = (1 - w) * p.dydr[j] + w * p.dydr[j+1]
        acc += (r / ρ) / abs(dd)
    end
    return acc
end

"""
    atmosphere_flux(profile, ρ) -> flux

Normalised stellar flux at distance `ρ` [km] from the shadow centre: the refracted flux
(averaged over the stellar disc: in two dimensions within 20 stellar radii of the shadow
centre, where the central flash forms, as strips across the limb farther out) times the
diffraction pattern of the opaque surface edge at `profile.y_surface` in the shadow plane.
"""
function atmosphere_flux(p::AtmosphereProfile, ρ::Real)
    # inside the geometric shadow of the surface the refracted flux is continued with its
    # value at the shadow edge; the Fresnel pattern of the edge then cuts it off smoothly
    ys = max(p.y_surface, 0.0)
    φref(x) = _refracted_flux(p, max(abs(x), ys))
    rs = p.star_radius
    φ = if rs <= 0
        φref(ρ)
    elseif abs(ρ) < 20rs
        _disc_average(φref, abs(ρ), rs)
    else                                       # far from the centre: strips across the limb
        ξ, w = _disc_nodes(rs, 0.05)
        sum(w[j] * φref(ρ + rs * ξ[j]) for j in eachindex(ξ))
    end
    return φ * _P(p.edge.x, p.edge.Pb, abs(ρ) - p.y_surface)
end

# mean of a circularly symmetric flux f(ρ′) over a uniform stellar disc of radius a centred at
# distance ρ from the shadow centre: (1/πa²) ∫ f(ρ′) L(ρ′) dρ′, L = length of the circle of
# radius ρ′ inside the disc (exact near the central flash, where the 1-D strip average is not)
function _disc_average(f, ρ, a; n::Integer=200)
    lo, hi = max(ρ - a, 0.0), ρ + a
    c, h = (lo + hi) / 2, (hi - lo) / 2
    acc = 0.0
    for i in 1:n
        t = (i - 0.5) * π / n
        ρ′ = c + h * cos(t)
        L = ρ′ <= a - ρ ? 2π * ρ′ : 2ρ′ * acos(clamp((ρ′^2 + ρ^2 - a^2) / (2ρ * ρ′), -1.0, 1.0))
        acc += f(ρ′) * L * h * sin(t)
    end
    return acc * (π / n) / (π * a^2)
end

"""
    half_light_radius(profile) -> km

Distance from the shadow centre at which the refracted flux is 1/2 (NaN if the flux at the
surface ray exceeds 1/2).
"""
function half_light_radius(p::AtmosphereProfile)
    f(ρ) = _refracted_flux(p, ρ) - 0.5
    # outermost crossing of 1/2, scanning inwards from the top of the atmosphere
    ρs = range(p.y[end], max(p.y_surface, 0.0); length=4000)
    prev = ρs[1]
    for ρ in ρs[2:end]
        if f(ρ) < 0
            lo, hi = ρ, prev
            for _ in 1:60
                m = (lo + hi) / 2
                f(m) < 0 ? (lo = m) : (hi = m)
            end
            return (lo + hi) / 2
        end
        prev = ρ
    end
    return NaN
end
