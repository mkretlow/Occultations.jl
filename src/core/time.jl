"""
    Epoch

An instant of time in UTC, stored as an ERFA two-part quasi Julian Date
(`jd1 + jd2`), which keeps sub-microsecond resolution and handles leap seconds.

Constructors

    Epoch("2020-06-03T12:30:45.1234")   # ISO string (UTC); 'T' or space separator
    Epoch(DateTime(2020, 6, 3, 12))     # Dates.DateTime (UTC)
    Epoch(2459004.02; scale=:utc)       # Julian date in :utc, :tt or :tdb
    Epoch(2020, 6, 3, 12, 30, 45.1234)  # calendar components (UTC)

Arithmetic: `e + 1.5` (seconds, SI), `e2 - e1` (seconds, SI), comparisons.

# Examples
```jldoctest
julia> t = Epoch("2024-01-10T15:00:00");

julia> t + 1.5
Epoch("2024-01-10T15:00:01.5000")

julia> round(Epoch("2024-01-10 15:00:30") - t; digits=6)
30.0

julia> round(Epoch("2017-01-01T00:00:00") - Epoch("2016-12-31T23:59:59"); digits=6)   # leap second
2.0
```
"""
struct Epoch
    jd1::Float64
    jd2::Float64
end

function Epoch(y::Integer, m::Integer, d::Integer, h::Integer=0, mi::Integer=0, s::Real=0.0)
    d1, d2 = ERFA.dtf2d("UTC", y, m, d, h, mi, Float64(s))
    return Epoch(d1, d2)
end

Epoch(dt::DateTime) = Epoch(year(dt), month(dt), day(dt), hour(dt), minute(dt),
                            second(dt) + millisecond(dt) / 1000)

const _ISO_RE = r"^\s*([+-]?\d{4})-(\d{1,2})-(\d{1,2})(?:[T\s]+(\d{1,2}):(\d{1,2})(?::(\d{1,2}(?:\.\d*)?))?)?\s*(?:Z|UTC)?\s*$"

function Epoch(s::AbstractString)
    m = match(_ISO_RE, s)
    m === nothing && throw(ArgumentError("cannot parse \"$s\" as an ISO time (YYYY-MM-DD[Thh:mm[:ss.sss]])"))
    y, mo, d = parse.(Int, m.captures[1:3])
    h  = m.captures[4] === nothing ? 0 : parse(Int, m.captures[4])
    mi = m.captures[5] === nothing ? 0 : parse(Int, m.captures[5])
    sec = m.captures[6] === nothing ? 0.0 : parse(Float64, m.captures[6])
    return Epoch(y, mo, d, h, mi, sec)
end

function Epoch(jd::Real; scale::Symbol=:utc)
    j1 = floor(jd) + 0.5 > jd ? floor(jd) - 0.5 : floor(jd) + 0.5
    j2 = jd - j1
    return _from_scale(j1, Float64(j2), scale)
end

function _from_scale(j1, j2, scale::Symbol)
    if scale === :utc
        return Epoch(j1, j2)
    elseif scale === :tai
        return Epoch(ERFA.taiutc(j1, j2)...)
    elseif scale === :tt
        return Epoch(ERFA.taiutc(ERFA.tttai(j1, j2)...)...)
    elseif scale === :tdb
        # TT ≈ TDB - dtdb; one iteration is ample (|dtdb| < 2 ms)
        dt = ERFA.dtdb(j1, j2, mod(j2 + 0.5, 1.0), 0.0, 0.0, 0.0) / DAY_S
        return Epoch(ERFA.taiutc(ERFA.tttai(j1, j2 - dt)...)...)
    else
        throw(ArgumentError("unknown time scale $scale (use :utc, :tai, :tt or :tdb)"))
    end
end

Epoch(e::Epoch) = e

"Two-part TAI Julian date."
tai_jd(e::Epoch) = ERFA.utctai(e.jd1, e.jd2)
"Two-part TT Julian date."
tt_jd(e::Epoch) = ERFA.taitt(tai_jd(e)...)
"Two-part UTC (quasi) Julian date."
utc_jd(e::Epoch) = (e.jd1, e.jd2)
"Two-part UT1 Julian date for a given `dut1 = UT1 - UTC` [s]."
ut1_jd(e::Epoch; dut1::Real=0.0) = ERFA.utcut1(e.jd1, e.jd2, Float64(dut1))

"Two-part TDB Julian date (geocentric TDB-TT from ERFA.dtdb)."
function tdb_jd(e::Epoch)
    t1, t2 = tt_jd(e)
    dt = ERFA.dtdb(t1, t2, mod(e.jd1 - 0.5, 1.0) + e.jd2, 0.0, 0.0, 0.0)
    return t1, t2 + dt / DAY_S
end

"Ephemeris time (TDB seconds past J2000), as used by SPICE."
function et(e::Epoch)
    t1, t2 = tdb_jd(e)
    return ((t1 - JD_J2000) + t2) * DAY_S
end

"Julian epoch (year) in TDB."
function julian_year(e::Epoch)
    t1, t2 = tdb_jd(e)
    return 2000.0 + ((t1 - JD_J2000) + t2) / JULIAN_YEAR
end

function Base.:+(e::Epoch, s::Real)
    t1, t2 = tai_jd(e)
    t2 += s / DAY_S
    # renormalise so that the second part stays small
    k = floor(t2)
    return Epoch(ERFA.taiutc(t1 + k, t2 - k)...)
end
Base.:-(e::Epoch, s::Real) = e + (-s)
Base.:+(s::Real, e::Epoch) = e + s

function Base.:-(a::Epoch, b::Epoch)
    a1, a2 = tai_jd(a)
    b1, b2 = tai_jd(b)
    return ((a1 - b1) + (a2 - b2)) * DAY_S
end

Base.isless(a::Epoch, b::Epoch) = (a - b) < 0
Base.:(==)(a::Epoch, b::Epoch) = abs(a - b) < 1e-9
Base.hash(e::Epoch, h::UInt) = hash(round((e.jd1 - 2.4e6 + e.jd2) * DAY_S * 1e6), h)
Base.broadcastable(e::Epoch) = Ref(e)

"""
    isot(e::Epoch; digits=4)

ISO-8601 representation of `e` in UTC.
"""
function isot(e::Epoch; digits::Integer=4)
    iy, im, id, h, mi, s, f = ERFA.d2dtf("UTC", digits, e.jd1, e.jd2)
    fs = digits > 0 ? "." * lpad(string(f), digits, '0') : ""
    return @sprintf("%04d-%02d-%02dT%02d:%02d:%02d", iy, im, id, h, mi, s) * fs
end

Base.show(io::IO, e::Epoch) = print(io, "Epoch(\"", isot(e), "\")")
Base.print(io::IO, e::Epoch) = print(io, isot(e))

"Convert to `Dates.DateTime` (millisecond resolution)."
function Dates.DateTime(e::Epoch)
    iy, im, id, h, mi, s, f = ERFA.d2dtf("UTC", 3, e.jd1, e.jd2)
    return DateTime(iy, im, id, h, mi, min(s, 59), f)
end

_epoch(t::Epoch) = t
_epoch(t::AbstractString) = Epoch(t)
_epoch(t::DateTime) = Epoch(t)
