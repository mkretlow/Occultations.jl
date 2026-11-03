# ---------------------------------------------------------------------------
# Angles and spherical/cartesian helpers. All public angles are in degrees.
# ---------------------------------------------------------------------------

function _parse_sexagesimal(s::AbstractString)
    t = strip(replace(s, r"[hdmsHDM:'\"°]" => " "))
    sgn = startswith(t, "-") ? -1.0 : 1.0
    parts = parse.(Float64, split(lstrip(t, ['+', '-'])))
    isempty(parts) && throw(ArgumentError("cannot parse angle \"$s\""))
    v = 0.0
    for (i, p) in enumerate(parts)
        v += p / 60.0^(i - 1)
    end
    return sgn * v
end

"""
    parse_ra(x) -> degrees

Right ascension from a number (degrees) or a sexagesimal string in hours
(`"12 34 56.789"`, `"12:34:56.789"`, `"12h34m56.789s"`).

# Examples
```jldoctest
julia> ra = parse_ra("16 04 37.609"); round(ra; digits=8)
241.15670417

julia> format_ra(ra)
"16 04 37.6090"

julia> format_dec(parse_dec("-16 51 19.78"))
"-16 51 19.780"
```
"""
parse_ra(x::Real) = Float64(x)
parse_ra(s::AbstractString) = 15.0 * _parse_sexagesimal(s)

"""
    parse_dec(x) -> degrees

Declination from a number (degrees) or a sexagesimal string (`"-12 34 56.78"`).
"""
parse_dec(x::Real) = Float64(x)
parse_dec(s::AbstractString) = _parse_sexagesimal(s)

"Longitude/latitude parser (degrees or sexagesimal degrees)."
parse_angle(x::Real) = Float64(x)
parse_angle(s::AbstractString) = _parse_sexagesimal(s)

function _sexa(x::Real, digits::Integer)
    s = sign(x) < 0 ? "-" : "+"
    x = abs(x)
    d = floor(x)
    m = floor((x - d) * 60)
    sec = ((x - d) * 60 - m) * 60
    sec = round(sec; digits=digits)
    if sec >= 60
        sec -= 60; m += 1
    end
    if m >= 60
        m -= 60; d += 1
    end
    w = digits > 0 ? digits + 3 : 2
    return s, Int(d), Int(m), lpad(@sprintf("%.*f", digits, sec), w, '0')
end

"Format a right ascension in degrees as `hh mm ss.ssss`."
function format_ra(ra::Real; digits::Integer=4)
    _, h, m, s = _sexa(mod(ra, 360) / 15, digits)
    return @sprintf("%02d %02d %s", h, m, s)
end

"Format a declination in degrees as `±dd mm ss.sss`."
function format_dec(dec::Real; digits::Integer=3)
    sg, d, m, s = _sexa(dec, digits)
    return @sprintf("%s%02d %02d %s", sg, d, m, s)
end

"Unit vector from RA/Dec in degrees."
function radec2vec(ra::Real, dec::Real)
    α, δ = ra * DEG2RAD, dec * DEG2RAD
    cδ = cos(δ)
    return [cδ * cos(α), cδ * sin(α), sin(δ)]
end

"RA [0,360), Dec in degrees and norm of a cartesian vector."
function vec2radec(v::AbstractVector)
    r = norm(v)
    ra = mod(atan(v[2], v[1]) * RAD2DEG, 360.0)
    dec = asin(clamp(v[3] / r, -1.0, 1.0)) * RAD2DEG
    return ra, dec, r
end

"East and north unit vectors of the local tangent plane at direction `u`."
function sky_basis(u::AbstractVector)
    ρ = hypot(u[1], u[2])
    e = ρ > 0 ? [-u[2] / ρ, u[1] / ρ, 0.0] : [0.0, 1.0, 0.0]
    n = cross(u / norm(u), e)
    return e, n
end

"""
    angular_separation(ra1, dec1, ra2, dec2) -> degrees
"""
function angular_separation(ra1, dec1, ra2, dec2)
    u1, u2 = radec2vec(ra1, dec1), radec2vec(ra2, dec2)
    return atan(norm(cross(u1, u2)), dot(u1, u2)) * RAD2DEG
end

"""
    position_angle(ra1, dec1, ra2, dec2) -> degrees

Position angle of point 2 relative to point 1, from north towards east, in [0,360).
"""
function position_angle(ra1, dec1, ra2, dec2)
    e, n = sky_basis(radec2vec(ra1, dec1))
    u2 = radec2vec(ra2, dec2)
    return mod(atan(dot(u2, e), dot(u2, n)) * RAD2DEG, 360.0)
end

"""
    offset_radec(ra, dec, dα, dδ) -> (ra, dec)

Apply a tangent-plane offset (`dα = Δα cos δ`, `dδ`) in milliarcseconds.
"""
function offset_radec(ra::Real, dec::Real, dα::Real, dδ::Real)
    (dα == 0 && dδ == 0) && return Float64(ra), Float64(dec)
    u = radec2vec(ra, dec)
    e, n = sky_basis(u)
    v = u .+ (dα * MAS2RAD) .* e .+ (dδ * MAS2RAD) .* n
    r, d, _ = vec2radec(v)
    return r, d
end

"Rotation matrix ITRS -> GCRS transpose (i.e. `GCRS = M' * ITRS` with `M = c2t`)."
function celestial_to_terrestrial(t::Epoch; dut1::Real=0.0, xp::Real=0.0, yp::Real=0.0)
    tt1, tt2 = tt_jd(t)
    u1, u2 = ut1_jd(t; dut1=dut1)
    return ERFA.c2t06a(tt1, tt2, u1, u2, xp * ARCSEC2RAD, yp * ARCSEC2RAD)
end

"Mean obliquity of the ecliptic (IAU 2006) at J2000 [rad]."
const OBLIQUITY_J2000 = 84381.406 * ARCSEC2RAD

"Rotation ecliptic (J2000) -> equatorial (ICRS)."
function ecl2equ(v::AbstractVector)
    ε = OBLIQUITY_J2000
    c, s = cos(ε), sin(ε)
    return [v[1], c * v[2] - s * v[3], s * v[2] + c * v[3]]
end
