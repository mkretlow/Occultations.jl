# ---------------------------------------------------------------------------
# Diffraction model of an occultation by an airless body: Fresnel diffraction at the
# limb (integrated over the passband), the finite stellar disc and the exposure time
# ---------------------------------------------------------------------------

"""
    fresnel_scale(λ, dist) -> km

Fresnel scale `sqrt(λ D / 2)` for wavelength `λ` [µm] and distance `dist` [km].
"""
fresnel_scale(λ::Real, dist::Real) = sqrt(λ * 1e-9 * dist / 2)

# complex Fresnel integral C(z) + i S(z)
@inline _fresnelE(z::Real) = (1 + im) / 2 * SpecialFunctions.erf((sqrt(π) / 2) * (1 - im) * z)

# knife-edge amplitude through the open half-plane; x > 0 on the illuminated side.
# For a body of amplitude transmission τ the intensity is |τ + (1-τ)A|²
#   = τ² + 2τ(1-τ) Re(A) + (1-τ)² |A|²
@inline function _edge_amplitude(x, F)
    E = _fresnelE(-x / F)
    return (1 - im) / 2 * ((0.5 - real(E)) + im * (0.5 - imag(E)))
end

# opaque strip [0, L] along track; x is the along-track coordinate
@inline function _strip_intensity(x, L, F, τ)
    A = (1 - im) / 2 * (_fresnelE((L - x) / F) - _fresnelE(-x / F))
    return abs2(1 - (1 - τ) * A)
end

# Gauss–Chebyshev (2nd kind) nodes/weights for a uniform stellar disc: ∫ g(ξ) √(1-ξ²) dξ
function _disc_nodes(R, F)
    R <= 0 && return [0.0], [1.0]
    n = clamp(ceil(Int, 6R / F), 5, 101)
    ξ = [cos(k * π / (n + 1)) for k in 1:n]
    w = [sin(k * π / (n + 1))^2 for k in 1:n]
    return ξ, w ./ sum(w)
end

function _band_nodes(λ, Δλ)
    Δλ <= 0 && return [λ]
    n = 7
    return collect(range(λ - Δλ / 2, λ + Δλ / 2; length=n))
end

# Tabulated knife-edge profile, averaged over band and stellar disc: mean Re(A) and
# mean |A|² (independent of the opacity) and their cumulative integrals
struct EdgeProfile
    x::Vector{Float64}
    Pa::Vector{Float64}
    Pb::Vector{Float64}
    Qa::Vector{Float64}
    Qb::Vector{Float64}
end

function _edge_profile(dist, star_diam, λ, Δλ)
    lams = _band_nodes(λ, Δλ)
    Fs = fresnel_scale.(lams, dist)
    Fmin, Fmax = extrema(Fs)
    R = star_diam / 2
    ξ, w = _disc_nodes(R, Fmin)
    X = 100Fmax + 2R
    dx = min(Fmin / 20, R > 0 ? R / 10 : Inf)
    n = min(ceil(Int, 2X / dx) + 1, 40_001)
    xs = collect(range(-X, X; length=n))
    Pa = zeros(n); Pb = zeros(n)
    Threads.@threads for i in 1:n
        sa = 0.0; sb = 0.0
        for F in Fs, (ξj, wj) in zip(ξ, w)
            A = _edge_amplitude(xs[i] + R * ξj, F)
            sa += wj * real(A); sb += wj * abs2(A)
        end
        Pa[i] = sa / length(Fs); Pb[i] = sb / length(Fs)
    end
    cum(P) = (Q = zeros(n); for i in 2:n; Q[i] = Q[i-1] + (P[i] + P[i-1]) / 2 * (xs[i] - xs[i-1]); end; Q)
    return EdgeProfile(xs, Pa, Pb, cum(Pa), cum(Pb))
end

function _interp(xs, ys, x)
    i = clamp(searchsortedlast(xs, x), 1, length(xs) - 1)
    w = (x - xs[i]) / (xs[i+1] - xs[i])
    return (1 - w) * ys[i] + w * ys[i+1]
end

# ∫ profile dx from x[1]; beyond the table the amplitudes are 0 (shadow) or 1 (light)
# Beyond the table the profile is continued continuously towards its limits (0 deep in the
# shadow, 1 in the light) with the asymptotic decay of the fringes (∝ 1/x² resp. 1/x); a
# hard cut would make the model jump when a sample crosses the end of the table.
function _Q(xs, Q, x)
    x1, xe = xs[1], xs[end]
    x <= x1 && return _P1(Q, xs) * x1^2 * (1 / x1 - 1 / x)
    x >= xe && return Q[end] + (x - xe) + (_Pe(Q, xs) - 1) * xe * log(x / xe)
    return _interp(xs, Q, x)
end
# end values of the profile, from the slope of its cumulative integral
_P1(Q, xs) = (Q[2] - Q[1]) / (xs[2] - xs[1])
_Pe(Q, xs) = (Q[end] - Q[end-1]) / (xs[end] - xs[end-1])

function _P(xs, P, x)
    x <= xs[1] && return P[1] * (xs[1] / x)^2
    x >= xs[end] && return 1 + (P[end] - 1) * xs[end] / x
    return _interp(xs, P, x)
end

# mean intensity over [xa, xb] for amplitude transmission τ
function _Pmean(p::EdgeProfile, τ, xa, xb)
    xa > xb && ((xa, xb) = (xb, xa))
    if xb - xa < 1e-9
        x = (xa + xb) / 2
        a, b = _P(p.x, p.Pa, x), _P(p.x, p.Pb, x)
    else
        a = (_Q(p.x, p.Qa, xb) - _Q(p.x, p.Qa, xa)) / (xb - xa)
        b = (_Q(p.x, p.Qb, xb) - _Q(p.x, p.Qb, xa)) / (xb - xa)
    end
    return τ^2 + 2τ * (1 - τ) * a + (1 - τ)^2 * b
end

# short-chord model: direct strip computation with exposure sub-sampling
function _strip_flux(t, texp, tim, tem, v, dist, star_diam, λ, Δλ, τ)
    lams = _band_nodes(λ, Δλ)
    Fs = fresnel_scale.(lams, dist)
    R = star_diam / 2
    ξ, w = _disc_nodes(R, minimum(Fs))
    L = v * (tem - tim)
    ns = clamp(ceil(Int, v * texp / (minimum(Fs) / 4)), 1, 200)
    out = similar(t, Float64)
    for k in eachindex(t)
        acc = 0.0
        for s in 1:ns
            ts = t[k] + texp * ((s - 0.5) / ns - 0.5)
            x = v * (ts - tim)
            a = 0.0
            for F in Fs, (ξj, wj) in zip(ξ, w)
                a += wj * _strip_intensity(x + R * ξj, L, F, τ)
            end
            acc += a / length(Fs)
        end
        out[k] = acc / ns
    end
    return out
end

# amplitude transmission √(1 − opacity); exactly 0 for an opaque body (no infinite derivative
# of the square root at 0 under automatic differentiation)
_transmission(opacity) = (x = clamp(1 - opacity, 0, 1); x > 0 ? sqrt(x) : zero(x))

function _model_intensity(t, tim, tem, vi, ve, texp, prof, dist, star_diam, λ, Δλ, opacity)
    F = fresnel_scale(λ + Δλ / 2, dist)
    vm = (vi + ve) / 2
    L = vm * (tem - tim)
    if isfinite(L) && L < 10F + 2star_diam && L > 0
        # short event: diffraction by an opaque strip, crossed at the mean normal velocity
        return _strip_flux(t, texp, tim, tem, vm, dist, star_diam, λ, Δλ, _transmission(opacity))
    end
    tmid = isfinite(tim) && isfinite(tem) ? (tim + tem) / 2 : (isfinite(tem) ? -Inf : Inf)
    out = similar(t, promote_type(Float64, typeof(tim), typeof(tem), typeof(opacity)))   # generic for AD
    h = texp / 2
    τ = _transmission(opacity)
    @inbounds for k in eachindex(t)
        tk = t[k]
        if tem <= tim
            out[k] = 1.0
        elseif tk < tmid
            out[k] = _Pmean(prof, τ, vi * (tim - tk - h), vi * (tim - tk + h))
        else
            out[k] = _Pmean(prof, τ, ve * (tk - h - tem), ve * (tk + h - tem))
        end
    end
    return out
end


# several drops: intervals [tims[k], tems[k]] (sorted) crossed at the normal velocities
# vims[k], vems[k]. Drops separated by more than ~10 Fresnel scales and the stellar diameter
# are independent: the product of their single-drop models (each ≈ 1 far from its drop). Closer
# drops (e.g. the narrow gap between the lobes of a contact binary, partly filled by
# diffraction) are computed together: the complex amplitudes of the opaque strips add
# (Babinet), U = 1 − (1 − τ) Σₖ Aₖ.
function _multi_intensity(t, tims, tems, vims, vems, texp, prof, dist, star_diam, λ, Δλ, opacity)
    n = length(tims)
    op(k) = opacity isa AbstractVector ? opacity[k] : opacity
    n == 0 && return ones(promote_type(Float64, eltype(t)), length(t))
    n == 1 && return _model_intensity(t, tims[1], tems[1], vims[1], vems[1], texp, prof, dist, star_diam, λ, Δλ, op(1))
    F = fresnel_scale(λ + Δλ / 2, dist)
    sep = all(k -> (tims[k+1] - tems[k]) * min(vems[k], vims[k+1]) > 10F + 2star_diam + vems[k] * texp, 1:n-1)
    if sep
        I = _model_intensity(t, tims[1], tems[1], vims[1], vems[1], texp, prof, dist, star_diam, λ, Δλ, op(1))
        for k in 2:n
            I = I .* _model_intensity(t, tims[k], tems[k], vims[k], vems[k], texp, prof, dist, star_diam, λ, Δλ, op(k))
        end
        return I
    end
    return _multistrip_flux(t, texp, tims, tems, (vims .+ vems) ./ 2, dist, star_diam, λ, Δλ,
                            [_transmission(op(k)) for k in 1:n])
end

# several strips (amplitudes added), each with its amplitude transmission τ[k] (0: opaque), with
# exposure sub-sampling, passband and stellar disc
function _multistrip_flux(t, texp, tims, tems, vs, dist, star_diam, λ, Δλ, τ)
    lams = _band_nodes(λ, Δλ)
    Fs = fresnel_scale.(lams, dist)
    R = star_diam / 2
    ξ, w = _disc_nodes(R, minimum(Fs))
    ns = clamp(ceil(Int, maximum(vs) * texp / (minimum(Fs) / 4)), 1, 200)
    out = similar(t, Float64)
    for k in eachindex(t)
        acc = 0.0
        for s in 1:ns
            ts = t[k] + texp * ((s - 0.5) / ns - 0.5)
            a = 0.0
            for F in Fs, (ξj, wj) in zip(ξ, w)
                A = zero(ComplexF64)
                for m in eachindex(tims)
                    x = vs[m] * (ts - tims[m]) + R * ξj
                    L = vs[m] * (tems[m] - tims[m])
                    A += (1 - τ[m]) * (1 - im) / 2 * (_fresnelE((L - x) / F) - _fresnelE(-x / F))
                end
                a += wj * abs2(1 - A)
            end
            acc += a / length(Fs)
        end
        out[k] = acc / ns
    end
    return out
end

"""
    Optics(; λ=0.7, Δλ=0.3, exptime=0)

Observation parameters entering the light-curve model: central wavelength `λ` and
bandwidth `Δλ` [µm] (the band is sampled uniformly) and exposure (integration) time
`exptime` [s].
"""
Base.@kwdef struct Optics
    λ::Float64 = 0.7
    Δλ::Float64 = 0.3
    exptime::Float64 = 0.0
end

"""
    occultation_flux(t, t_imm, t_emm; vel, dist, star_diam=0, λ=0.7, Δλ=0.3, exptime=0,
                     vel_emersion=vel, opacity=1) -> normalised flux

Model intensity (1 = unocculted) at times `t` [s, mid-exposure] for immersion and emersion
instants `t_imm`, `t_emm` [s]: Fresnel diffraction at each limb crossing (band `λ ± Δλ/2`
[µm]), smoothing by the projected stellar disc of diameter `star_diam` [km] and integration
over the exposure `exptime` [s]. `vel`, `vel_emersion` are the shadow velocities normal to
the limb [km/s], `dist` the distance of the body [km], `opacity` the fraction of the
stellar flux blocked by the body. Events shorter than about ten Fresnel scales are
computed as the diffraction pattern of an opaque strip. Use `t_imm = -Inf` or
`t_emm = Inf` for a single edge.

`t_imm` and `t_emm` are geometric contacts of the stellar centre with the projected
limb, before stellar-disc and exposure averaging. At an isolated opaque straight edge,
a monochromatic point source with zero exposure averaging has intensity 0.25 at contact;
this is not a general threshold for the averaged or partly transparent light curve.

    occultation_flux(t, drops; vel, kwargs...) -> normalised flux

Several drops (a binary system, the lobes of a contact binary): `drops = [(t_imm, t_emm), …]`
[s]; `vel` (and `vel_emersion`) and `opacity` one value or one per drop (e.g. the partly
transparent dips of rings next to the opaque drop of the body). Drops closer than about ten
Fresnel scales are computed together (the amplitudes of the opaque strips add, so the light
diffracted into a narrow gap is included).
"""
function occultation_flux(t::AbstractVector{<:Real}, drops::AbstractVector{<:Tuple}; vel, dist::Real,
                          star_diam::Real=0.0, λ::Real=0.7, Δλ::Real=0.3, exptime::Real=0.0,
                          vel_emersion=vel, opacity=1.0, profile=nothing)
    o = sortperm(collect(drops); by=first)
    d = collect(drops)[o]
    n = length(d)
    vi = vel isa Real ? fill(abs(vel), n) : abs.(collect(vel))[o]
    ve = vel_emersion isa Real ? fill(abs(vel_emersion), n) : abs.(collect(vel_emersion))[o]
    op = opacity isa Real ? fill(Float64(opacity), n) : collect(Float64, opacity)[o]
    prof = profile === nothing ? _edge_profile(dist, star_diam, λ, Δλ) : profile
    return _multi_intensity(t, first.(d), last.(d), vi, ve, exptime, prof, dist, star_diam, λ, Δλ, op)
end

function occultation_flux(t::AbstractVector{<:Real}, t_imm::Real, t_emm::Real; vel::Real, dist::Real,
                          star_diam::Real=0.0, λ::Real=0.7, Δλ::Real=0.3, exptime::Real=0.0,
                          vel_emersion::Real=vel, opacity::Real=1.0, profile=nothing)
    prof = profile === nothing ? _edge_profile(dist, star_diam, λ, Δλ) : profile
    return _model_intensity(t, t_imm, t_emm, abs(vel), abs(vel_emersion), exptime, prof,
                            dist, star_diam, λ, Δλ, opacity)
end

# cached edge profiles (they depend only on the distance, stellar diameter and band)
const _PROFILE_CACHE = Dict{NTuple{4,Float64},EdgeProfile}()
const _PROFILE_LOCK = ReentrantLock()
function _cached_profile(dist, star_diam, λ, Δλ)
    key = (round(dist; sigdigits=6), round(star_diam; sigdigits=6), Float64(λ), Float64(Δλ))
    lock(_PROFILE_LOCK) do
        length(_PROFILE_CACHE) > 200 && empty!(_PROFILE_CACHE)
        get!(() -> _edge_profile(key...), _PROFILE_CACHE, key)
    end
end
