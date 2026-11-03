# ---------------------------------------------------------------------------
# Fundamental-plane geometry
#
# (f, g) is the orthographic projection of the star relative to the occulting
# body, in km at the distance of the body, as seen by the observer:
# f positive towards local East (increasing RA), g towards North. The body
# (ephemeris position) is at the origin.
# ---------------------------------------------------------------------------

# f, g of a star with unit direction `us`, for a body at topocentric position `Bt` [km]
@inline function _fg(us::AbstractVector, Bt::AbstractVector)
    D = norm(Bt)
    ub = Bt ./ D
    e, n = sky_basis(ub)
    return D * dot(us, e), D * dot(us, n), D
end

"Barycentric velocity of the Earth [km/s] (ICRS axes)."
function earth_velocity(t::Epoch)
    t1, t2 = tdb_jd(t)
    _, pvb = ERFA.epv00(t1, t2)
    return pvb[2] .* (AU_KM / DAY_S)
end

# Approximate barycentric velocity of the body [km/s]: d/dt of the astrometric geocentric
# vector plus the barycentric velocity of the Earth.
function _body_velocity(eph::AbstractEphemeris, t::Epoch; h::Real=30.0)
    return (geocentric_vector(eph, t + h) .- geocentric_vector(eph, t - h)) ./ (2h) .+ earth_velocity(t)
end

# Topocentric astrometric vector of the body from the astrometric geocentric vector `B`,
# the barycentric body velocity `V` and the observer position `ρ`: the body is taken at
# the emission time of the light reaching the observer (first order in the light-time
# difference, which is ≤ 21 ms for ground-based observers).
@inline function _topocentric(B::AbstractVector, V::AbstractVector, ρ::AbstractVector)
    Bt = B .- ρ
    δτ = (norm(Bt) - norm(B)) / C_KMS
    return Bt .- V .* δτ
end

function _topocentric_vector(eph::AbstractEphemeris, obs::AbstractObserver, t::Epoch; dut1::Real=0.0)
    B = geocentric_vector(eph, t)
    _isgeocenter(obs) && return B
    return _topocentric(B, _body_velocity(eph, t), gcrs_position(obs, t; dut1=dut1))
end
_topocentric_vector(b::Body, obs::AbstractObserver, t::Epoch; kw...) = _topocentric_vector(_ephem(b), obs, t; kw...)
"""
    topocentric_vector(body, observer, t; dut1=0) -> [x, y, z] km

Astrometric position of the body relative to `observer` (ICRS axes), with the
body taken at the light emission time for that observer.
"""
topocentric_vector(body, obs::AbstractObserver, t; kw...) = _topocentric_vector(_ephem(body), obs, _epoch(t); kw...)

# observers farther than this from the geocentre [km] get the stellar parallax applied
const _PARALLAX_DISTANCE = 1e5

function _fg_at(star::Star, eph::AbstractEphemeris, obs::AbstractObserver, t::Epoch; dut1=0.0)
    ρ = gcrs_position(obs, t; dut1=dut1)
    ra, dec = norm(ρ) > _PARALLAX_DISTANCE ? geocentric(star, t; offset=ρ) : geocentric(star, t)
    us = radec2vec(ra, dec)
    B = geocentric_vector(eph, t)
    Bt = _isgeocenter(obs) ? B : _topocentric(B, _body_velocity(eph, t), ρ)
    return _fg(us, Bt)
end

"""
    positionv(star, body, observer, t; dt=1.0, dut1=0.0)
        -> (f, g, vf, vg, dist)

Position (f, g) [km] of the star relative to the body in the fundamental plane
as seen by `observer` at `t`, its velocity (vf, vg) [km/s] (central difference
with step `dt` seconds) and the topocentric distance of the body [km].
`body` may be a [`Body`](@ref) or any ephemeris. `t` may be a vector.
"""
function positionv(star::Star, body, obs::AbstractObserver, t; dt::Real=1.0, dut1::Real=0.0)
    eph = _ephem(body)
    t = _epoch(t)
    f, g, D = _fg_at(star, eph, obs, t; dut1=dut1)
    f1, g1, _ = _fg_at(star, eph, obs, t - dt / 2; dut1=dut1)
    f2, g2, _ = _fg_at(star, eph, obs, t + dt / 2; dut1=dut1)
    return (f=f, g=g, vf=(f2 - f1) / dt, vg=(g2 - g1) / dt, dist=D)
end
positionv(star::Star, body, obs::AbstractObserver, t::AbstractVector; kw...) = [positionv(star, body, obs, ti; kw...) for ti in t]
positionv(star::Star, body, t; kw...) = positionv(star, body, GEOCENTER, t; kw...)

"""
    shadow_velocity(star, body, t; observer=GEOCENTER) -> km/s

Magnitude of the velocity of the star relative to the body in the fundamental plane.
"""
function shadow_velocity(star::Star, body, t; observer::AbstractObserver=GEOCENTER)
    p = positionv(star, body, observer, t)
    return hypot(p.vf, p.vg)
end

"""
    closest_approach(star, body, tstart, tend; observer=GEOCENTER, step=60)
        -> (epoch, ca_km, ca_arcsec, pa, vel, dist_km)

Instant of minimum (f, g) distance between `tstart` and `tend`. `pa` is the
position angle of the body relative to the star [deg, from north through east].
"""
function closest_approach(star::Star, body, tstart, tend; observer::AbstractObserver=GEOCENTER, step::Real=60.0)
    eph = _ephem(body)
    t0, t1 = _epoch(tstart), _epoch(tend)
    span = t1 - t0
    n = max(2, ceil(Int, span / step) + 1)
    ts = range(0, span; length=n)
    dfun(s) = (r = _fg_at(star, eph, observer, t0 + s); hypot(r[1], r[2]))
    ds = dfun.(ts)
    i = argmin(ds)
    a, b = ts[max(i - 1, 1)], ts[min(i + 1, n)]
    s, _ = _golden_min(dfun, a, b; tol=1e-4)
    tca = t0 + s
    p = positionv(star, eph, observer, tca)
    d = hypot(p.f, p.g)
    return (epoch=tca, ca_km=d, ca_arcsec=d / p.dist * RAD2ARCSEC,
            pa=mod(atan(-p.f, -p.g) * RAD2DEG, 360.0), vel=hypot(p.vf, p.vg), dist_km=p.dist)
end

# ---------------------------------------------------------------------------
# Shadow paths on the Earth
# ---------------------------------------------------------------------------

"""
    ShadowPath

Result of [`shadow_path`](@ref): `lon[i, k]`, `lat[i, k]` [deg] of the point
of the Earth reached by the line `k` (cross-track offset `offsets[k]` km)
at `times[i]`. Points off the Earth are `NaN`.
"""
struct ShadowPath
    times::Vector{Epoch}
    offsets::Vector{Float64}
    lon::Matrix{Float64}
    lat::Matrix{Float64}
    sun_alt::Matrix{Float64}
    body::String
end

function Base.show(io::IO, ::MIME"text/plain", p::ShadowPath)
    println(io, "ShadowPath of ", p.body, ": ", length(p.times), " epochs × ", length(p.offsets),
            " lines (offsets ", p.offsets, " km)")
    k = argmin(abs.(p.offsets))
    for i in eachindex(p.times)
        isnan(p.lon[i, k]) && continue
        @printf(io, "  %s  lon=%9.4f  lat=%8.4f  sun_alt=%6.1f\n", isot(p.times[i]; digits=1),
                p.lon[i, k], p.lat[i, k], p.sun_alt[i, k])
    end
end

# point of the ellipsoid where the observer sees the star at (ftarget, gtarget)
function _ground_point(us, e, n, M, B, V, ftarget, gtarget, f0, g0)
    a = R_EARTH_KM
    b = a * (1 - F_EARTH)
    Dg = [1.0, 1.0, a / b]
    ξ, η = ftarget - f0, gtarget - g0
    d = Dg .* (M * us)
    ρ = fill(NaN, 3)
    for _ in 1:4
        w = Dg .* (M * (ξ .* e .+ η .* n))
        A = dot(d, d); Bq = dot(w, d); C = dot(w, w) - a^2
        disc = Bq^2 - A * C
        disc < 0 && return fill(NaN, 3)
        ζ = (-Bq + sqrt(disc)) / A
        ρ = ξ .* e .+ η .* n .+ ζ .* us
        f, g, _ = _fg(us, _topocentric(B, V, ρ))
        ξ -= f - ftarget
        η -= g - gtarget
    end
    return ρ
end

"""
    shadow_path(star, body, tstart, tend; step=10, radius=body.radius,
                offsets=[0, -radius, radius], dut1=0) -> ShadowPath

Compute the ground track of the shadow centre and of parallel lines at
cross-track `offsets` [km] (positive to the left of the motion in the
(f, g) plane), between `tstart` and `tend` every `step` seconds.
"""
function shadow_path(star::Star, body, tstart, tend; step::Real=10.0,
                     radius::Real=_radius(body), offsets=nothing, dut1::Real=0.0)
    eph = _ephem(body)
    t0, t1 = _epoch(tstart), _epoch(tend)
    offs = offsets === nothing ? (radius > 0 ? [0.0, -radius, radius] : [0.0]) : Float64.(collect(offsets))
    n = max(2, floor(Int, (t1 - t0) / step) + 1)
    times = [t0 + (i - 1) * step for i in 1:n]
    lon = fill(NaN, n, length(offs)); lat = similar(lon); salt = similar(lon)
    for (i, t) in enumerate(times)
        ra, dec = geocentric(star, t)
        us = radec2vec(ra, dec)
        e, nn = sky_basis(us)
        B = geocentric_vector(eph, t)
        V = _body_velocity(eph, t)
        f0, g0, _ = _fg(us, B)
        p = positionv(star, eph, GEOCENTER, t)
        v = hypot(p.vf, p.vg)
        pf, pg = -p.vg / v, p.vf / v          # left-hand normal to the motion
        M = celestial_to_terrestrial(t; dut1=dut1)
        for (k, δ) in enumerate(offs)
            ρ = _ground_point(us, e, nn, M, B, V, δ * pf, δ * pg, f0, g0)
            isnan(ρ[1]) && continue
            λ, φ, _ = ERFA.gc2gd(ERFA.WGS84, (M * ρ) .* 1000)
            lon[i, k] = λ * RAD2DEG
            lat[i, k] = φ * RAD2DEG
            salt[i, k] = sun_altitude(lon[i, k], lat[i, k], t)
        end
    end
    name = body isa Body ? body.name : eph.name
    return ShadowPath(times, offs, lon, lat, salt, name)
end

"""
    sun_altitude(lon, lat, t) -> degrees

Geometric altitude of the Sun for a site (geodetic lon/lat in degrees).
"""
function sun_altitude(lon::Real, lat::Real, t)
    t = _epoch(t)
    s = -earth_heliocentric(t)
    st = celestial_to_terrestrial(t) * (s ./ norm(s))
    λ, φ = lon * DEG2RAD, lat * DEG2RAD
    up = [cos(φ) * cos(λ), cos(φ) * sin(λ), sin(φ)]
    return asin(clamp(dot(up, st), -1, 1)) * RAD2DEG
end
sun_altitude(o::Observer, t) = sun_altitude(o.lon, o.lat, t)
