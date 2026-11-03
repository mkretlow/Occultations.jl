# ---------------------------------------------------------------------------
# Refraction by a thin, spherically symmetric atmosphere (Pluto, 2002 XV93, ...):
# hydrostatic density profile, ray tracing of the bending angle, flux in the shadow
# (near and far limb, limb-curvature focusing), as e.g. in Sicardy et al. and
# Arimatsu et al. (2026, arXiv:2605.02243)
# ---------------------------------------------------------------------------

# molecular mass [amu] and molecular refractivity K [m³] (ν = K n) in the visible
const _GASES = Dict(
    :N2  => (μ=28.0134, K=1.091e-29, disp=6.282e-3),   # Washburn (1948): K (1 + 6.282e-3/λ²[µm])
    :CH4 => (μ=16.0425, K=1.637e-29, disp=0.0),
    :CO  => (μ=28.0101, K=1.258e-29, disp=0.0),
)

"""
    Atmosphere(gas=:N2; T0=40, dTdz=0, z_inv=0, z_fmf=Inf, μ, K, ztop=NaN)

Spherically symmetric, transparent atmosphere in hydrostatic equilibrium.

* `gas` – `:N2`, `:CH4` or `:CO` (molecular mass `μ` [amu] and refractivity `K` [m³];
  both can be given explicitly)
* temperature: `T0` [K] at the surface, rising (or falling) by `dTdz` [K/km] up to the
  altitude `z_inv` [km], isothermal above (a Pluto-like inversion layer; `dTdz = 0` gives
  an isothermal atmosphere)
* `z_fmf` [km]: above this altitude the density follows the free-molecular-flow profile
  `n ∝ (R + z)⁻²` (low Jeans parameter, Arimatsu et al. 2026)
* `ztop` [km]: top of the model (default: 30 scale heights above the inversion)
* `width` [km]: the top of the inversion layer and the onset of the free molecular flow are
  smooth transitions of this width (a kink in the density gradient would produce a
  spurious jump of the refracted flux, since φ ∝ 1/(1 + D dω/dr))

The free parameter of a fit is the surface pressure (`pressure`, nbar); the radius comes
from the limb model and the gravity from the body (`GM` or `density`). See
[`atmosphere_profile`](@ref), [`xv93_atmosphere`](@ref), [`pluto_atmosphere`](@ref).
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
end

function Atmosphere(gas::Symbol=:N2; T0::Real=40.0, dTdz::Real=0.0, z_inv::Real=0.0, z_fmf::Real=Inf,
                    μ=nothing, K=nothing, ztop::Real=NaN, width::Real=2.0)
    g = get(_GASES, gas, nothing)
    g === nothing && (μ === nothing || K === nothing) &&
        throw(ArgumentError("unknown gas $gas: give μ [amu] and K [m³]"))
    width > 0 || throw(ArgumentError("width must be positive"))
    return Atmosphere(gas, something(μ, g === nothing ? NaN : g.μ), something(K, g === nothing ? NaN : g.K),
                      T0, dTdz, z_inv, z_fmf, ztop, width)
end

Base.show(io::IO, a::Atmosphere) = @printf(io, "Atmosphere(%s, T0 = %.1f K, dT/dz = %.2f K/km up to %.1f km%s)",
                                           a.gas, a.T0, a.dTdz, a.z_inv, isfinite(a.z_fmf) ? @sprintf(", FMF above %.1f km", a.z_fmf) : "")

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

_softplus(x) = x > 30 ? x : log1p(exp(x))
_logistic(x) = 1 / (1 + exp(-x))
# smooth min(z, z_inv): the temperature gradient switches off over `width`
_temperature(a::Atmosphere, z) = a.T0 + a.dTdz * (a.z_inv - a.width * _softplus((a.z_inv - z) / a.width) +
                                                  a.width * _softplus(a.z_inv / a.width) - a.z_inv)
_dTdz(a::Atmosphere, z) = a.dTdz * _logistic((a.z_inv - z) / a.width)

function _refractivity_K(a::Atmosphere, λ)
    g = get(_GASES, a.gas, nothing)
    return g === nothing || g.disp == 0 ? a.K : a.K * (1 + g.disp / λ^2)
end

"""
    AtmosphereProfile

Result of [`atmosphere_profile`](@ref): bending angle `ω(r)` [rad] versus closest-approach
radius `r` [km] of the rays and their position `y = r + ω D` in the shadow plane.
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
end

"""
    atmosphere_profile(atm, pressure, R, GM, D; λ=0.7, Δλ=0.3, star_diam=0, nz=3000, nr=1500, nu=300)
        -> AtmosphereProfile

Ray tracing through the atmosphere `atm` with surface `pressure` [nbar] on a body of radius
`R` [km] and `GM` [km³/s²], seen from distance `D` [km]:

1. number density from the hydrostatic equation `dp/dz = -μ n g(z)`, `p = n k T(z)`,
   `g = GM/(R+z)²`, blended smoothly into `n ∝ (R+z)⁻²` above the free-molecular-flow
   altitude (this tail is integrated analytically to infinity: it dominates the bending of
   tenuous atmospheres);
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
    H0 = K_BOLTZMANN * atm.T0 / (μ * gs) / 1e3
    Hinv = K_BOLTZMANN * _temperature(atm, atm.z_inv) / (μ * GM * 1e9 / ((R + atm.z_inv) * 1e3)^2) / 1e3
    fmf = isfinite(atm.z_fmf)
    ztab = fmf ? atm.z_fmf + 10atm.width : (isnan(atm.ztop) ? atm.z_inv + 30Hinv : atm.ztop)
    z = collect(range(0.0, ztab; length=nz))
    dz = z[2] - z[1]
    lnn = zeros(nz)
    n0 = pressure * 1e-4 / (K_BOLTZMANN * atm.T0)                         # nbar → Pa → m⁻³
    lnn[1] = log(max(n0, 1e-300))
    # d ln n / dz [1/km]: hydrostatic, blended smoothly into the free-molecular-flow law -2/(R+z)
    function dlog(zz)
        T = _temperature(atm, zz)
        g = GM * 1e9 / ((R + zz) * 1e3)^2
        hyd = -(μ * g / (K_BOLTZMANN * T)) * 1e3 - _dTdz(atm, zz) / T
        fmf || return hyd
        w = _logistic((zz - atm.z_fmf) / atm.width)
        return (1 - w) * hyd + w * (-2 / (R + zz))
    end
    for i in 2:nz
        lnn[i] = lnn[i-1] + (dlog(z[i-1]) + dlog(z[i])) / 2 * dz
    end
    ν = K .* exp.(lnn)
    dν = [(ν[min(i + 1, nz)] - ν[max(i - 1, 1)]) / (z[min(i + 1, nz)] - z[max(i - 1, 1)]) for i in 1:nz]
    ρf = R + ztab
    A = fmf ? ν[end] * ρf^2 : 0.0                                         # ν = A/ρ² above ρf
    νp_tab(ρ) = (x = (ρ - R) / dz; x >= nz - 1 ? dν[nz] : (k = max(floor(Int, x), 0); w = x - k; (1 - w) * dν[k+1] + w * dν[k+2]))
    # largest radius worth tracing: where |D dω/dr| < 1e-5
    # outer radius of the traced rays (independent of the pressure, so that the grid does not
    # move with it): far enough for the free-molecular-flow tail of any plausible atmosphere
    rmax = fmf ? max(R + 15R, ρf + 10H0) : ρf
    xs = range(0.0, 1.0; length=nr)
    rs = [R + (rmax - R) * x^2 for x in xs]
    ω = zeros(nr)
    Threads.@threads for j in 1:nr
        r = rs[j]
        if fmf && r >= ρf
            ω[j] = -π * A / r^2                                           # exact for ν ∝ ρ⁻²
            continue
        end
        uf = sqrt(max(ρf - r, 0.0))
        acc = 0.0
        if uf > 0
            h = uf / nu
            for i in 0:nu
                u = i * h
                acc += ((i == 0 || i == nu) ? 0.5 : 1.0) * νp_tab(r + u^2) / sqrt(2r + u^2) * h
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
    return AtmosphereProfile(R, D, pressure, rs, ω, y, dydr, y[1], star_diam / 2, edge, H0)
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
(averaged over the stellar disc) times the diffraction pattern of the opaque surface edge
at `profile.y_surface` in the shadow plane.
"""
function atmosphere_flux(p::AtmosphereProfile, ρ::Real)
    # inside the geometric shadow of the surface the refracted flux is continued with its
    # value at the shadow edge; the Fresnel pattern of the edge then cuts it off smoothly
    ys = max(p.y_surface, 0.0)
    φref(x) = _refracted_flux(p, max(abs(x), ys))
    φ = if p.star_radius > 0
        ξ, w = _disc_nodes(p.star_radius, 0.05)
        sum(w[j] * φref(ρ + p.star_radius * ξ[j]) for j in eachindex(ξ))
    else
        φref(ρ)
    end
    return φ * _P(p.edge.x, p.edge.Pb, abs(ρ) - p.y_surface)
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
