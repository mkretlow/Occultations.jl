"""
    AbstractEphemeris

Supertype of all ephemerides. A concrete ephemeris `E` must implement

    _geocentric_vector(eph::E, t::Epoch) -> Vector{Float64}

returning the astrometric (light-time corrected) geocentric position of the
body in km, ICRS axes. User offsets (see [`add_offset!`](@ref)) are applied by
[`geocentric_vector`](@ref).

Every ephemeris has fields `offset_ra`, `offset_dec` [mas].
"""
abstract type AbstractEphemeris end

"""
    geocentric_vector(eph, t) -> [x, y, z] km

Astrometric geocentric position of the body (ICRS axes), including any offset.
"""
function geocentric_vector(eph::AbstractEphemeris, t)
    t = _epoch(t)
    v = _geocentric_vector(eph, t)
    if eph.offset_ra != 0 || eph.offset_dec != 0
        ra, dec, r = vec2radec(v)
        ra, dec = offset_radec(ra, dec, eph.offset_ra, eph.offset_dec)
        v = r .* radec2vec(ra, dec)
    end
    return v
end

"""
    ephem_radec(eph, t) -> (ra, dec, dist_km)

Astrometric geocentric RA/Dec [deg] and distance [km].
"""
ephem_radec(eph::AbstractEphemeris, t) = vec2radec(geocentric_vector(eph, t))

"""
    add_offset!(eph, dα, dδ)

Add an offset to the ephemeris, `dα = Δα cos δ` and `dδ` in mas (cumulative).
"""
function add_offset!(eph::AbstractEphemeris, dα::Real, dδ::Real)
    eph.offset_ra += dα
    eph.offset_dec += dδ
    return eph
end

# ---------------------------------------------------------------------------
# Tabulated ephemeris
# ---------------------------------------------------------------------------

"""
    EphemTable(times, ra, dec, dist; name="", order=8)

Tabulated astrometric geocentric ephemeris (e.g. from JPL Horizons, quantity 1
with "astrometric RA/Dec" and observer range, or from any other source).

* `times` – vector of [`Epoch`](@ref)s (or strings/DateTimes), UTC
* `ra`, `dec` – astrometric ICRS coordinates [deg]
* `dist`  – geocentric distance [km]  (use `dist_unit=:au` for au)

Positions are interpolated in cartesian coordinates with a Lagrange
polynomial of `order` points.
"""
mutable struct EphemTable <: AbstractEphemeris
    name::String
    t0::Epoch
    times::Vector{Float64}        # seconds from t0
    xyz::Matrix{Float64}          # 3 × n, km
    order::Int
    offset_ra::Float64
    offset_dec::Float64
end

function EphemTable(times::AbstractVector, ra::AbstractVector, dec::AbstractVector,
                    dist::AbstractVector; name::AbstractString="", order::Integer=8,
                    dist_unit::Symbol=:km)
    n = length(times)
    n == length(ra) == length(dec) == length(dist) || throw(DimensionMismatch("inconsistent table lengths"))
    n >= 2 || throw(ArgumentError("at least two ephemeris points are required"))
    ep = _epoch.(times)
    p = sortperm(ep)
    t0 = ep[p[1]]
    scale = dist_unit === :au ? AU_KM : 1.0
    ts = [ep[i] - t0 for i in p]
    xyz = reduce(hcat, [Float64(dist[i]) * scale .* radec2vec(parse_ra(ra[i]), parse_dec(dec[i])) for i in p])
    return EphemTable(String(name), t0, ts, xyz, min(order, n), 0.0, 0.0)
end

"""
    EphemTable(times, xyz::AbstractMatrix; name="", order=8)

Tabulated geocentric position vectors (`3 × n`, km, ICRS axes).
"""
function EphemTable(times::AbstractVector, xyz::AbstractMatrix; name::AbstractString="", order::Integer=8)
    size(xyz, 1) == 3 || throw(DimensionMismatch("xyz must be 3 × n"))
    n = length(times)
    n == size(xyz, 2) || throw(DimensionMismatch("inconsistent table lengths"))
    n >= 2 || throw(ArgumentError("at least two ephemeris points are required"))
    ep = _epoch.(times)
    p = sortperm(ep)
    t0 = ep[p[1]]
    return EphemTable(String(name), t0, [ep[i] - t0 for i in p], Matrix{Float64}(xyz[:, p]), min(order, n), 0.0, 0.0)
end

function _lagrange(ts::AbstractVector, ys::AbstractMatrix, x::Real)
    n = length(ts)
    out = zeros(size(ys, 1))
    for j in 1:n
        w = 1.0
        for k in 1:n
            k == j && continue
            w *= (x - ts[k]) / (ts[j] - ts[k])
        end
        out .+= w .* view(ys, :, j)
    end
    return out
end

function _geocentric_vector(eph::EphemTable, t::Epoch)
    x = t - eph.t0
    ts = eph.times
    n = length(ts)
    tol = 0.5 * (ts[end] - ts[1]) / max(n - 1, 1)
    (x < ts[1] - tol || x > ts[end] + tol) &&
        @warn "EphemTable $(eph.name): $(isot(t)) is outside the tabulated interval; extrapolating" maxlog = 3
    i = searchsortedlast(ts, x)
    m = eph.order
    lo = clamp(i - m ÷ 2 + 1, 1, n - m + 1)
    idx = lo:lo+m-1
    return _lagrange(view(ts, idx), view(eph.xyz, :, idx), x)
end

Base.show(io::IO, e::EphemTable) = print(io, "EphemTable(", e.name, ", ", length(e.times),
    " points, ", isot(e.t0), " + ", round(e.times[end] / 3600; digits=2), " h)")

# ---------------------------------------------------------------------------
# Two-body ephemeris from osculating elements
# ---------------------------------------------------------------------------

const GM_SUN_AU3_D2 = 2.9591220828559115e-4    # k² [au³/day²]

"""
    EphemKeplerian(; a, e, i, Ω, ω, M, epoch, name="")

Heliocentric two-body ephemeris from osculating elements referred to the
J2000 ecliptic: semimajor axis `a` [au], eccentricity `e`, inclination `i`,
longitude of node `Ω`, argument of perihelion `ω`, mean anomaly `M` [deg] at
`epoch` (an [`Epoch`](@ref), or a TDB Julian date). Earth positions come from
ERFA `epv00`; light time is iterated. Adequate for tests, planning and
simulations — use Horizons or SPICE for real reductions.
"""
mutable struct EphemKeplerian <: AbstractEphemeris
    name::String
    a::Float64
    e::Float64
    i::Float64
    Ω::Float64
    ω::Float64
    M::Float64
    epoch_tdb::Float64
    offset_ra::Float64
    offset_dec::Float64
end

function EphemKeplerian(a::Real, e::Real, i::Real, Ω::Real, ω::Real, M::Real, epoch; name="")
    jd = epoch isa Real ? Float64(epoch) : sum(tdb_jd(_epoch(epoch)))
    return EphemKeplerian(String(name), a, e, i, Ω, ω, M, jd, 0.0, 0.0)
end
EphemKeplerian(; a, e, i, Ω, ω, M, epoch, name="") = EphemKeplerian(a, e, i, Ω, ω, M, epoch; name=name)

function _kepler(M, e)
    E = e < 0.8 ? M : π
    for _ in 1:50
        dE = (E - e * sin(E) - M) / (1 - e * cos(E))
        E -= dE
        abs(dE) < 1e-14 && break
    end
    return E
end

"Heliocentric ICRS position [au] at TDB Julian date `jd`."
function heliocentric_position(k::EphemKeplerian, jd::Real)
    n = sqrt(GM_SUN_AU3_D2 / k.a^3)                      # rad/day
    M = mod2pi(k.M * DEG2RAD + n * (jd - k.epoch_tdb))
    E = _kepler(M, k.e)
    xo = k.a * (cos(E) - k.e)
    yo = k.a * sqrt(1 - k.e^2) * sin(E)
    Ω, ω, i = k.Ω * DEG2RAD, k.ω * DEG2RAD, k.i * DEG2RAD
    cΩ, sΩ, cω, sω, ci, si = cos(Ω), sin(Ω), cos(ω), sin(ω), cos(i), sin(i)
    x = (cΩ * cω - sΩ * sω * ci) * xo + (-cΩ * sω - sΩ * cω * ci) * yo
    y = (sΩ * cω + cΩ * sω * ci) * xo + (-sΩ * sω + cΩ * cω * ci) * yo
    z = (sω * si) * xo + (cω * si) * yo
    return ecl2equ([x, y, z])
end

function _geocentric_vector(k::EphemKeplerian, t::Epoch)
    t1, t2 = tdb_jd(t)
    jd = t1 + t2
    pvh, _ = ERFA.epv00(t1, t2)
    rE = pvh[1]
    τ = 0.0
    v = heliocentric_position(k, jd) .- rE
    for _ in 1:3
        τ = norm(v) * AU_KM / C_KMS / DAY_S
        v = heliocentric_position(k, jd - τ) .- rE
    end
    return v .* AU_KM
end

Base.show(io::IO, k::EphemKeplerian) = print(io, "EphemKeplerian(", k.name, ", a=", k.a, " au, e=", k.e, ")")

# ---------------------------------------------------------------------------
# SPICE kernels (methods in the SPICE extension)
# ---------------------------------------------------------------------------

"""
    EphemKernel(target, kernels; name=target)

Ephemeris computed with SPICE from the given kernel files (a planetary kernel
such as `de440.bsp` plus the body's SPK). `target` is the NAIF id or name.
`abcorr` is the SPICE aberration correction: `"CN"` (astrometric, for occulting
bodies) or `"NONE"` (geometric, for [`Spacecraft`](@ref) observers).
Requires `using SPICE`.
"""
mutable struct EphemKernel <: AbstractEphemeris
    name::String
    target::String
    kernels::Vector{String}
    abcorr::String
    offset_ra::Float64
    offset_dec::Float64
    loaded::Bool
end

function EphemKernel(target, kernels::AbstractVector; name=string(target), abcorr::AbstractString="CN")
    return EphemKernel(String(name), string(target), String.(kernels), String(abcorr), 0.0, 0.0, false)
end
EphemKernel(target, kernel::AbstractString; kw...) = EphemKernel(target, [kernel]; kw...)

# implemented by the SPICE extension
function _spice_position end

function _geocentric_vector(k::EphemKernel, t::Epoch)
    isempty(methods(_spice_position)) && error("EphemKernel requires the SPICE package: run `using SPICE` first")
    return _spice_position(k, t)
end

# ---------------------------------------------------------------------------
# Online services (methods in the HTTP extension)
# ---------------------------------------------------------------------------

"""
    horizons_ephem(target, tstart, tend; step="1m", name=target, id_type=nothing)

Query JPL Horizons for the astrometric geocentric ephemeris of `target`
(e.g. `"136199"`, `"Chariklo"`, `"2060;"`, `"599"`) and return an [`EphemTable`](@ref).
Requires `using HTTP`.
"""
function horizons_ephem end

"""
    query_sbdb(target) -> NamedTuple

Physical parameters of a small body from the JPL Small-Body Database
(`radius`, `H`, `G`, `rotation_period`, `pole`, `density` [kg/m³], `GM`; missing values
are left out), usable as keyword arguments of
[`Body`](@ref): `Body("Chariklo", eph; query_sbdb("Chariklo")...)`. Requires `using HTTP`.
"""
function query_sbdb end


"""
    query_gaia(source_id; release="dr3") -> Star
    query_gaia(ra, dec; radius_arcsec=2, release="dr3") -> Vector{Star}

Retrieve stars from the Gaia archive (ESA TAP service). Requires `using HTTP`.
"""
function query_gaia end

# ---------------------------------------------------------------------------
# Body
# ---------------------------------------------------------------------------

"""
    Body(name, ephem; radius=NaN, diameter=NaN, H=NaN, G=0.15, pole=nothing,
         rotation_period=NaN, W0=0, frame=nothing, shape=nothing, GM=NaN, density=NaN)

An occulting body: an ephemeris plus physical parameters.

* `radius` – reference (equivalent) radius [km] (or give `diameter`), used for predictions
  and as the default starting value of limb fits
* `H`, `G` – absolute magnitude and slope ([`apparent_magnitude`](@ref))
* `pole = (ra, dec)` [deg], `rotation_period` [h] (negative: retrograde) and `W0` (prime
  meridian at J2000 TDB) define the body frame ([`BodyFrame`](@ref)) unless `frame` is given
* `shape` – a 3-D shape ([`Ellipsoid`](@ref) or [`TriMesh`](@ref)) in the body frame
* `GM` [km³/s²] or `density` [kg/m³] – needed for atmosphere models
"""
Base.@kwdef mutable struct Body
    name::String
    ephem::AbstractEphemeris
    radius::Float64 = NaN
    H::Float64 = NaN
    G::Float64 = 0.15
    frame::Union{Nothing,BodyFrame} = nothing
    shape::Union{Nothing,AbstractShape} = nothing
    GM::Float64 = NaN
    density::Float64 = NaN
end

function Body(name::AbstractString, ephem::AbstractEphemeris; diameter::Real=NaN, pole=nothing,
              rotation_period::Real=NaN, W0::Real=0.0, frame=nothing, kw...)
    kw = Dict{Symbol,Any}(kw)
    isnan(diameter) || haskey(kw, :radius) || (kw[:radius] = diameter / 2)
    for k in (:radius, :H, :G, :GM, :density)
        haskey(kw, k) && (kw[k] = Float64(kw[k]))
    end
    if frame === nothing && pole !== nothing
        frame = BodyFrame(parse_ra(pole[1]), parse_dec(pole[2]), rotation_period; W0=W0)
    end
    return Body(; name=String(name), ephem=ephem, frame=frame, kw...)
end

geocentric_vector(b::Body, t) = geocentric_vector(b.ephem, t)
ephem_radec(b::Body, t) = ephem_radec(b.ephem, t)

"""
    body_gm(body) -> km³/s²

`GM` of the body, from `GM` or from `density` and `radius`.
"""
function body_gm(b::Body)
    isfinite(b.GM) && return b.GM
    (isfinite(b.density) && isfinite(b.radius)) ||
        throw(ArgumentError("body $(b.name): give `GM` or `density` (and `radius`)"))
    return G_NEWTON * b.density * 4π / 3 * (b.radius * 1e3)^3 * 1e-9
end

Base.show(io::IO, b::Body) = print(io, "Body(", b.name, isnan(b.radius) ? "" : ", R=$(b.radius) km", ")")

function Base.show(io::IO, ::MIME"text/plain", b::Body)
    println(io, "Body ", b.name)
    println(io, "  ephemeris: ", b.ephem)
    isnan(b.radius) || @printf(io, "  radius = %.3f km\n", b.radius)
    isnan(b.H) || @printf(io, "  H = %.3f, G = %.3f\n", b.H, b.G)
    b.frame === nothing || println(io, "  frame: ", b.frame)
    b.shape === nothing || println(io, "  shape: ", b.shape)
    isnan(b.GM) || @printf(io, "  GM = %.6g km³/s²\n", b.GM)
    isnan(b.density) || @printf(io, "  density = %.1f kg/m³\n", b.density)
end

"""
    hg_magnitude(H, G, r, Δ, α) -> V

Apparent magnitude in the IAU H, G system (Bowell et al. 1989) for heliocentric
and observer distances `r`, `Δ` [au] and phase angle `α` [deg].
"""
function hg_magnitude(H::Real, G::Real, r::Real, Δ::Real, α::Real)
    t = tand(abs(α) / 2)
    Φ1 = exp(-3.33 * t^0.63)
    Φ2 = exp(-1.87 * t^1.22)
    return H + 5 * log10(r * Δ) - 2.5 * log10((1 - G) * Φ1 + G * Φ2)
end

"""
    apparent_magnitude(body, t; observer=GEOCENTER) -> V

Apparent V magnitude of `body` from its `H` and `G` (H, G system); the distances and
phase angle come from the ephemeris.
"""
function apparent_magnitude(b::Body, t; observer::AbstractObserver=GEOCENTER)
    isnan(b.H) && throw(ArgumentError("body $(b.name): H is not defined"))
    t = _epoch(t)
    Bt = _topocentric_vector(b.ephem, observer, t)
    Bh = geocentric_vector(b.ephem, t) .+ earth_heliocentric(t) .* AU_KM
    α = acosd(clamp(dot(Bh, Bt) / (norm(Bh) * norm(Bt)), -1, 1))
    return hg_magnitude(b.H, b.G, norm(Bh) / AU_KM, norm(Bt) / AU_KM, α)
end

_ephem(b::Body) = b.ephem
_ephem(e::AbstractEphemeris) = e
_radius(b::Body) = isnan(b.radius) ? 0.0 : b.radius
_radius(::AbstractEphemeris) = 0.0

# ---------------------------------------------------------------------------
# Fundamental-plane (Besselian) ephemeris, e.g. from a published prediction
# ---------------------------------------------------------------------------

"""
    EphemPlane(star, epoch, x, y; dist, units=:km, name="")
    EphemPlane(star, epoch, x0, y0, dx, dy; dist, units=:km, name="")

Ephemeris of an occulting body given only by its shadow on the fundamental plane, as in
the Besselian elements of an occultation prediction: the shadow axis (body centre) is at

    x(τ) = x[1] + x[2] τ + x[3] τ² + …,   y(τ) likewise,   τ = t − epoch,

relative to the geocentre, in the plane through the geocentre perpendicular to the
direction of `star`, `x` towards east and `y` towards north (so the geocentric (f, g) of
[`positionv`](@ref) are `(−x, −y)`). With `units=:km` the coefficients are in km, km/s,
km/s², …; with `units=:earth_radii` in Earth radii, Earth radii/hour, … (CORA, Occult).
`dist` is the geocentric distance of the body [km] (needed for the light time, the
Fresnel scale and the angular conversion of the astrometry).

The body is placed so that the occultation of `star` (the same [`Star`](@ref) must be used
for the [`Occultation`](@ref)) has exactly this geometry; topocentric observers, maps, light
curves and limb fits then work as with any other ephemeris. The predicted path is
reproduced, not the body's orbit: use it only near `epoch`.
See [`parse_cora`](@ref) for CORA predictions.
"""
mutable struct EphemPlane <: AbstractEphemeris
    name::String
    star::Star
    epoch::Epoch
    x::Vector{Float64}          # km, km/s, km/s², …
    y::Vector{Float64}
    dist::Float64
    offset_ra::Float64
    offset_dec::Float64
end

function EphemPlane(star::Star, epoch, x::AbstractVector{<:Real}, y::AbstractVector{<:Real}; dist::Real,
                    units::Symbol=:km, name::AbstractString="")
    units in (:km, :earth_radii) || throw(ArgumentError("units must be :km or :earth_radii"))
    dist > 0 || throw(ArgumentError("dist must be positive [km]"))
    conv(c) = units === :km ? Float64.(c) : [R_EARTH_KM * c[k] / 3600.0^(k - 1) for k in eachindex(c)]
    return EphemPlane(String(name), star, _epoch(epoch), conv(x), conv(y), Float64(dist), 0.0, 0.0)
end
EphemPlane(star::Star, epoch, x0::Real, y0::Real, dx::Real, dy::Real; kw...) =
    EphemPlane(star, epoch, [x0, dx], [y0, dy]; kw...)

function Base.show(io::IO, e::EphemPlane)
    @printf(io, "EphemPlane(%s at %s: x = %.1f km, y = %.1f km, v = (%.4f, %.4f) km/s)", e.name,
            isot(e.epoch; digits=1), e.x[1], e.y[1], get(e.x, 2, 0.0), get(e.y, 2, 0.0))
end

_poly(c, τ) = foldr((ck, acc) -> ck + τ * acc, c; init=0.0)

function _geocentric_vector(e::EphemPlane, t::Epoch)
    τ = t - e.epoch
    ft, gt = -_poly(e.x, τ), -_poly(e.y, τ)
    us = radec2vec(geocentric(e.star, t)...)
    eb, nb = sky_basis(us)
    B = e.dist .* us .- ft .* eb .- gt .* nb
    for _ in 1:3                                     # make (f, g) exact (second-order terms)
        f, g, _ = _fg(us, B)
        B .+= (f - ft) .* eb .+ (g - gt) .* nb
    end
    return B
end

"""
    parse_cora(row; scale=:utc) -> NamedTuple

Read one event of a CORA prediction (a `|`-separated row): the epoch of the elements
(field 24, a Julian date in time scale `scale`), the shadow position `x`, `y` [Earth radii]
and velocity `dx`, `dy` [Earth radii/hour] (fields 26–29), and the relative velocity
`vel` [km/s] (field 19), as read by Occultations1.jl's `occ_from_cora_csv`. Make an
ephemeris from it with

    c = parse_cora(row)
    eph = EphemPlane(star, c.epoch, c.x, c.y, c.dx, c.dy; units=:earth_radii, dist=…)

and compare `shadow_velocity(star, eph, c.epoch)` with `c.vel` as a check.
"""
function parse_cora(row::AbstractString; scale::Symbol=:utc)
    h = strip.(split(row, '|'))
    length(h) >= 29 || throw(ArgumentError("a CORA row needs at least 29 `|`-separated fields"))
    x, y, dx, dy = parse.(Float64, h[26:29])
    return (epoch=Epoch(parse(Float64, h[24]); scale=scale), x=x, y=y, dx=dx, dy=dy, vel=parse(Float64, h[19]))
end
