"""
    Observer(name, lon, lat, height=0.0)

Ground-based observing site. `lon` (east positive) and `lat` (geodetic, WGS84)
in degrees or sexagesimal strings (`"-45 30 12.5"`), `height` above the
ellipsoid in metres. See also [`parse_mpc_obscode`](@ref) and [`observer_from_mpc`](@ref).
"""
struct Observer <: AbstractObserver
    name::String
    lon::Float64
    lat::Float64
    height::Float64
    geocenter::Bool
end

Observer(name::AbstractString, lon, lat, height::Real=0.0) =
    Observer(String(name), parse_angle(lon), parse_angle(lat), Float64(height), false)

"The geocentre as a (virtual) observer."
const GEOCENTER = Observer("Geocenter", 0.0, 0.0, 0.0, true)

Base.show(io::IO, o::Observer) = o.geocenter ? print(io, "Observer(Geocenter)") :
    @printf(io, "Observer(%s, lon=%.6f°, lat=%.6f°, h=%.1f m)", o.name, o.lon, o.lat, o.height)

_isgeocenter(o::AbstractObserver) = false
_isgeocenter(o::Observer) = o.geocenter

"Terrestrial (ITRS) position of the observer [km]."
function itrs_position(o::Observer)
    o.geocenter && return zeros(3)
    return ERFA.gd2gc(ERFA.WGS84, o.lon * DEG2RAD, o.lat * DEG2RAD, o.height) ./ 1000
end

"""
    gcrs_position(observer, t; dut1=0.0) -> km

Geocentric position of the observer in the GCRS (ICRS-aligned axes),
IAU 2006/2000A precession-nutation; `dut1 = UT1-UTC` in seconds.
"""
function gcrs_position(o::Observer, t; dut1::Real=0.0)
    o.geocenter && return zeros(3)
    M = celestial_to_terrestrial(_epoch(t); dut1=dut1)
    return M' * itrs_position(o)
end

"""
    parse_mpc_obscode(line; name=nothing) -> Observer

Observer from one line of the MPC list of observatory codes
(`ObsCodes.html`: code, longitude [deg E], ρ cos φ′, ρ sin φ′ [Earth radii], name).
"""
function parse_mpc_obscode(line::AbstractString; name=nothing)
    line = rstrip(line)
    # fixed columns (code 1-3, longitude 4-13, ρ cos φ′ 14-21, ρ sin φ′ 22-30, name 31-);
    # whitespace-separated values are accepted as well
    v = length(line) >= 30 ? tryparse.(Float64, strip.((line[4:13], line[14:21], line[22:30]))) : (nothing,)
    if all(!isnothing, v)
        code, (lon, ρc, ρs), site = strip(line[1:3]), v, strip(line[31:end])
    else
        m = match(r"^(.{3})\s*(\d+\.\d*)\s+(\d\.\d*)\s*([+-]?\d\.\d*)\s*(.*)$", line)
        m === nothing && throw(ArgumentError("not an MPC observatory line with coordinates: \"$line\""))
        code, site = strip(m.captures[1]), strip(m.captures[5])
        lon, ρc, ρs = parse.(Float64, m.captures[2:4])
    end
    λ = lon * DEG2RAD
    x = R_EARTH_KM * 1000 .* [ρc * cos(λ), ρc * sin(λ), ρs]
    elon, φ, h = ERFA.gc2gd(ERFA.WGS84, x)
    nm = name === nothing ? (isempty(site) ? String(code) : "$(site) ($(code))") : String(name)
    return Observer(nm, mod(elon * RAD2DEG + 180, 360) - 180, φ * RAD2DEG, h)
end

"""
    sidereal_time(observer, t; mode=:local, dut1=0) -> hours

Apparent sidereal time (IAU 2006/2000A), local (`mode=:local`) or at Greenwich (`:greenwich`).
"""
function sidereal_time(o::Observer, t; mode::Symbol=:local, dut1::Real=0.0)
    e = _epoch(t)
    u1, u2 = ut1_jd(e; dut1=dut1)
    t1, t2 = tt_jd(e)
    gst = ERFA.gst06a(u1, u2, t1, t2) * RAD2DEG / 15
    mode === :greenwich && return gst
    mode === :local || throw(ArgumentError("mode must be :local or :greenwich"))
    return mod(gst + o.lon / 15, 24)
end

"""
    altaz(observer, t, ra, dec; dut1=0) -> (alt, az)
    altaz(observer, t, star::Star)

Geometric altitude and azimuth (from north through east) [deg] of the direction
`ra`, `dec` [deg, ICRS] for a ground-based observer (no refraction).
"""
function altaz(o::Observer, t, ra, dec; dut1::Real=0.0)
    o.geocenter && throw(ArgumentError("altaz is undefined for the geocentre"))
    M = celestial_to_terrestrial(_epoch(t); dut1=dut1)
    v = M * radec2vec(parse_ra(ra), parse_dec(dec))
    λ, φ = o.lon * DEG2RAD, o.lat * DEG2RAD
    up = [cos(φ) * cos(λ), cos(φ) * sin(λ), sin(φ)]
    east = [-sin(λ), cos(λ), 0.0]
    north = cross(up, east)
    return asind(clamp(dot(v, up), -1, 1)), mod(atand(dot(v, east), dot(v, north)), 360)
end
altaz(o::Observer, t, s::Star; kw...) = altaz(o, t, geocentric(s, t)...; kw...)

"""
    Spacecraft(name, ephem)
    Spacecraft(name, times, xyz)

An observer in space. `ephem` is any [`AbstractEphemeris`](@ref) that returns the
**geometric** geocentric position of the spacecraft (ICRS axes, km): a table of
state vectors (`xyz` is `3 × n` in km at UTC `times`), an [`EphemKernel`](@ref)
with `abcorr="NONE"`, or [`spacecraft_horizons`](@ref) (HTTP extension).
"""
struct Spacecraft{E<:AbstractEphemeris} <: AbstractObserver
    name::String
    ephem::E
end

Spacecraft(name::AbstractString, times::AbstractVector, xyz::AbstractMatrix; order::Integer=8) =
    Spacecraft(String(name), EphemTable(times, xyz; name=name, order=order))

Base.show(io::IO, s::Spacecraft) = print(io, "Spacecraft(", s.name, ")")

gcrs_position(s::Spacecraft, t; dut1::Real=0.0) = geocentric_vector(s.ephem, _epoch(t))

"""
    observer_from_mpc(code; name=nothing) -> Observer

Look up an observatory code in the MPC list (downloaded once per session) and
return the site as an [`Observer`](@ref). Requires `using HTTP`.
"""
function observer_from_mpc end

"""
    spacecraft_horizons(id, tstart, tend; step="10m", name=id) -> Spacecraft

[`Spacecraft`](@ref) observer from geometric geocentric state vectors computed by
JPL Horizons (e.g. `id = "-98"` for New Horizons). Requires `using HTTP`.
"""
function spacecraft_horizons end
