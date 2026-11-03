"""
    Star(; ra, dec, pmra=0, pmdec=0, parallax=0, rv=0, epoch=2016.0,
           ra_err=0, dec_err=0, pmra_err=0, pmdec_err=0, parallax_err=0,
           mag=Dict(), diameter=NaN, name="", source_id="")

A catalogue star.

* `ra`, `dec`  – ICRS coordinates at `epoch` [deg, or sexagesimal strings]
* `pmra`       – proper motion in RA, μα* = μα cos δ [mas/yr]
* `pmdec`      – proper motion in Dec [mas/yr]
* `parallax`   – [mas]; `rv` – radial velocity [km/s]
* `epoch`      – catalogue epoch as a Julian year (Gaia DR3: 2016.0)
* `*_err`      – uncertainties [mas, mas/yr]
* `mag`        – magnitudes, e.g. `Dict("G"=>12.3, "V"=>12.8, "K"=>10.1)`
* `diameter`   – angular diameter [mas]; estimated from magnitudes if NaN
                 and V and K are available (see [`star_diameter_kervella`](@ref)).
"""
Base.@kwdef mutable struct Star
    name::String = ""
    source_id::String = ""
    ra::Float64
    dec::Float64
    pmra::Float64 = 0.0
    pmdec::Float64 = 0.0
    parallax::Float64 = 0.0
    rv::Float64 = 0.0
    epoch::Float64 = 2016.0
    ra_err::Float64 = 0.0
    dec_err::Float64 = 0.0
    pmra_err::Float64 = 0.0
    pmdec_err::Float64 = 0.0
    parallax_err::Float64 = 0.0
    mag::Dict{String,Float64} = Dict{String,Float64}()
    diameter::Float64 = NaN
end

function Star(ra, dec; kwargs...)
    kw = Dict{Symbol,Any}(kwargs)
    for k in (:pmra, :pmdec, :parallax, :rv, :epoch, :ra_err, :dec_err, :pmra_err,
              :pmdec_err, :parallax_err, :diameter)
        haskey(kw, k) && (kw[k] = Float64(kw[k]))
    end
    haskey(kw, :source_id) && (kw[:source_id] = string(kw[:source_id]))
    if haskey(kw, :mag)
        kw[:mag] = Dict{String,Float64}(string(k) => Float64(v) for (k, v) in kw[:mag])
    end
    s = Star(; ra=parse_ra(ra), dec=parse_dec(dec), kw...)
    if isnan(s.diameter)
        d = _auto_diameter(s)
        isnan(d) || (s.diameter = d)
    end
    return s
end

function _auto_diameter(s::Star)
    m = s.mag
    if haskey(m, "V") && haskey(m, "K")
        return star_diameter_kervella(V=m["V"], K=m["K"])
    elseif haskey(m, "B") && haskey(m, "K")
        return star_diameter_kervella(B=m["B"], K=m["K"])
    end
    return NaN
end

_label(s::Star) = !isempty(s.name) ? s.name : !isempty(s.source_id) ? s.source_id : "unnamed"

function Base.show(io::IO, ::MIME"text/plain", s::Star)
    println(io, "Star ", _label(s))
    @printf(io, "  ICRS (epoch J%.3f): RA = %s  Dec = %s\n", s.epoch, format_ra(s.ra), format_dec(s.dec))
    @printf(io, "  pmRA* = %.3f ± %.3f mas/yr, pmDec = %.3f ± %.3f mas/yr\n", s.pmra, s.pmra_err, s.pmdec, s.pmdec_err)
    @printf(io, "  parallax = %.4f ± %.4f mas, RV = %.2f km/s\n", s.parallax, s.parallax_err, s.rv)
    isempty(s.mag) || println(io, "  magnitudes: ", join(["$k=$(round(v; digits=3))" for (k, v) in sort(collect(s.mag))], ", "))
    isnan(s.diameter) || @printf(io, "  angular diameter: %.4f mas\n", s.diameter)
end
Base.show(io::IO, s::Star) = print(io, "Star(", _label(s),
                                   ", RA=", format_ra(s.ra), ", Dec=", format_dec(s.dec), ")")

_jd_of_julian_year(y) = (JD_J2000, (y - 2000.0) * JULIAN_YEAR)

"""
    barycentric(star, t) -> (ra, dec, parallax)

Barycentric ICRS position [deg] and parallax [mas] of the star at epoch `t`,
using rigorous space motion (ERFA `pmsafe`).
"""
function barycentric(s::Star, t)
    t = _epoch(t)
    e1a, e1b = _jd_of_julian_year(s.epoch)
    e2a, e2b = tdb_jd(t)
    δ = s.dec * DEG2RAD
    pmr = s.pmra * MAS2RAD / cos(δ)      # dα/dt [rad/yr]
    pmd = s.pmdec * MAS2RAD
    px = max(s.parallax, 0.0) / 1000      # arcsec
    if pmr == 0 && pmd == 0 && s.rv == 0
        return s.ra, s.dec, s.parallax
    end
    ra2, dec2, _, _, px2, _ = ERFA.pmsafe(s.ra * DEG2RAD, δ, pmr, pmd, px, s.rv, e1a, e1b, e2a, e2b)
    return mod(ra2 * RAD2DEG, 360.0), dec2 * RAD2DEG, s.parallax > 0 ? px2 * 1000 : s.parallax
end

"Barycentric position of the Earth [au] (ICRS axes) at TDB epoch."
function earth_barycentric(t::Epoch)
    t1, t2 = tdb_jd(t)
    _, pvb = ERFA.epv00(t1, t2)
    return pvb[1]
end

"Heliocentric position of the Earth [au] (ICRS axes)."
function earth_heliocentric(t::Epoch)
    t1, t2 = tdb_jd(t)
    pvh, _ = ERFA.epv00(t1, t2)
    return pvh[1]
end

"""
    geocentric(star, t; offset=nothing) -> (ra, dec)

Astrometric geocentric ICRS position [deg] of the star at `t`: space motion
plus annual parallax (no aberration or light deflection, which are common to
the star and the occulting body). `offset` is an optional observer position
relative to the geocentre [km] (e.g. a spacecraft), included in the parallax.
"""
function geocentric(s::Star, t; offset=nothing)
    t = _epoch(t)
    ra, dec, px = barycentric(s, t)
    px <= 0 && return ra, dec
    d = PC_AU / (px / 1000)               # distance in au
    v = d .* radec2vec(ra, dec) .- earth_barycentric(t)
    offset === nothing || (v = v .- offset ./ AU_KM)
    r, de, _ = vec2radec(v)
    return r, de
end

"""
    position_error(star, t) -> (σα*, σδ) [mas]

Catalogue position uncertainty propagated to epoch `t` with the proper-motion errors.
"""
function position_error(s::Star, t)
    Δt = julian_year(_epoch(t)) - s.epoch
    return hypot(s.ra_err, Δt * s.pmra_err), hypot(s.dec_err, Δt * s.pmdec_err)
end

"""
    star_diameter_van_belle(; V, K, class=:ms) -> angular diameter [mas]

Angular diameter from the V–K colour following van Belle (1999, PASP 111, 1515).
`class` is `:ms` (main sequence), `:sg` (giants/supergiants) or `:vs` (variable stars).
"""
function star_diameter_van_belle(; V::Real, K::Real, class::Symbol=:ms)
    a, b = class === :sg ? (0.669, 0.223) :
           class === :ms ? (0.500, 0.264) :
           class === :vs ? (0.789, 0.218) :
           throw(ArgumentError("class must be :ms, :sg or :vs"))
    return 10.0^(a + b * (V - K) - 0.2 * V)
end

"""
    star_diameter_kervella(; K, V=nothing, B=nothing) -> limb-darkened angular diameter [mas]

Surface-brightness relations of Kervella et al. (2004, A&A 426, 297) for dwarfs and
sub-giants, using V–K (preferred) or B–K.
"""
function star_diameter_kervella(; K::Real, V=nothing, B=nothing)
    if V !== nothing
        return 10.0^(0.0755 * (V - K) + 0.5170 - 0.2 * K)
    elseif B !== nothing
        return 10.0^(0.0535 * (B - K) + 0.5159 - 0.2 * K)
    end
    throw(ArgumentError("V or B magnitude required"))
end

"""
    apparent_diameter(star, distance_km) -> km

Projected stellar diameter at the distance of the occulting body.
"""
apparent_diameter(s::Star, dist_km::Real) = isnan(s.diameter) ? 0.0 : s.diameter * MAS2RAD * dist_km
