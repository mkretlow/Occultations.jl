module OccultationsHTTPExt

using Occultations
using Occultations: Epoch, EphemTable, Star, Spacecraft, isot, _epoch, AU_KM, parse_mpc_obscode
using HTTP
using Dates
using Printf

const HORIZONS_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"

# request timeout keyword: `readtimeout` (HTTP 1.x) was replaced by `request_timeout` in HTTP 2
_timeout(s) = pkgversion(HTTP) >= v"2" ? (request_timeout=s,) : (readtimeout=s,)
const GAIA_TAP_URL = "https://gea.esac.esa.int/tap-server/tap/sync"

const _MONTHS = Dict("Jan" => 1, "Feb" => 2, "Mar" => 3, "Apr" => 4, "May" => 5, "Jun" => 6,
                     "Jul" => 7, "Aug" => 8, "Sep" => 9, "Oct" => 10, "Nov" => 11, "Dec" => 12)

function _horizons_time(s::AbstractString)
    m = match(r"(\d{4})-(\w{3})-(\d{2})\s+(\d{2}):(\d{2}):?(\d{2}(?:\.\d*)?)?", strip(s))
    m === nothing && error("cannot parse Horizons time \"$s\"")
    sec = m.captures[6] === nothing ? 0.0 : parse(Float64, m.captures[6])
    return Epoch(parse(Int, m.captures[1]), _MONTHS[m.captures[2]], parse(Int, m.captures[3]),
                 parse(Int, m.captures[4]), parse(Int, m.captures[5]), sec)
end

_hq(s) = "'" * s * "'"

function Occultations.horizons_ephem(target, tstart, tend; step::AbstractString="1m",
                                     name::AbstractString=string(target), center::AbstractString="500@399")
    t1, t2 = _epoch(tstart), _epoch(tend)
    fmt(t) = replace(isot(t; digits=3), "T" => " ")
    query = Dict(
        "format" => "text", "COMMAND" => _hq(string(target)), "OBJ_DATA" => "'NO'",
        "MAKE_EPHEM" => "'YES'", "EPHEM_TYPE" => "'OBSERVER'", "CENTER" => _hq(center),
        "START_TIME" => _hq(fmt(t1)), "STOP_TIME" => _hq(fmt(t2)), "STEP_SIZE" => _hq(step),
        "QUANTITIES" => "'1,20'", "ANG_FORMAT" => "'DEG'", "EXTRA_PREC" => "'YES'",
        "CSV_FORMAT" => "'YES'", "TIME_DIGITS" => "'FRACSEC'", "TIME_TYPE" => "'UT'",
        "REF_SYSTEM" => "'ICRF'", "APPARENT" => "'AIRLESS'")
    r = HTTP.get(HORIZONS_URL; query=query, _timeout(120)...)
    text = String(r.body)
    lines = split(text, '\n')
    i0 = findfirst(l -> startswith(strip(l), "\$\$SOE"), lines)
    i1 = findfirst(l -> startswith(strip(l), "\$\$EOE"), lines)
    if i0 === nothing || i1 === nothing
        error("Horizons did not return an ephemeris for \"$target\":\n" * text[1:min(end, 2000)])
    end
    hdr_i = findlast(l -> occursin("Date__", l), lines[1:i0])
    hdr = strip.(split(lines[hdr_i], ','))
    ira = findfirst(h -> startswith(h, "R.A."), hdr)
    idec = findfirst(h -> startswith(h, "DEC"), hdr)
    idel = findfirst(==("delta"), hdr)
    (ira === nothing || idec === nothing || idel === nothing) && error("unexpected Horizons header: $(lines[hdr_i])")
    times = Epoch[]; ra = Float64[]; dec = Float64[]; dist = Float64[]
    for l in lines[i0+1:i1-1]
        c = strip.(split(l, ','))
        length(c) < idel && continue
        push!(times, _horizons_time(c[1]))
        push!(ra, parse(Float64, c[ira]))
        push!(dec, parse(Float64, c[idec]))
        push!(dist, parse(Float64, c[idel]) * AU_KM)
    end
    return EphemTable(times, ra, dec, dist; name=name)
end

# --------------------------------------------------------------------------- Gaia

const _GAIA_COLS = "source_id, ra, dec, ra_error, dec_error, pmra, pmra_error, pmdec, pmdec_error, " *
                   "parallax, parallax_error, radial_velocity, phot_g_mean_mag, phot_bp_mean_mag, " *
                   "phot_rp_mean_mag, ref_epoch"

function _tap_csv(adql::AbstractString)
    body = HTTP.URIs.escapeuri(Dict("REQUEST" => "doQuery", "LANG" => "ADQL", "FORMAT" => "csv", "QUERY" => adql))
    r = HTTP.post(GAIA_TAP_URL, ["Content-Type" => "application/x-www-form-urlencoded"], body; _timeout(120)...)
    lines = filter(!isempty, split(String(r.body), r"\r?\n"))
    hdr = String.(split(lines[1], ','))
    rows = [Dict(zip(hdr, split(l, ','))) for l in lines[2:end]]
    return rows
end

_num(d, k) = (v = get(d, k, ""); isempty(v) ? NaN : parse(Float64, v))
_num0(d, k) = (v = _num(d, k); isnan(v) ? 0.0 : v)

function _gaia_star(d)
    mag = Dict{String,Float64}()
    for (k, b) in (("phot_g_mean_mag", "G"), ("phot_bp_mean_mag", "BP"), ("phot_rp_mean_mag", "RP"))
        v = _num(d, k)
        isnan(v) || (mag[b] = v)
    end
    return Star(_num(d, "ra"), _num(d, "dec"); source_id=d["source_id"], name="Gaia " * d["source_id"],
                pmra=_num0(d, "pmra"), pmdec=_num0(d, "pmdec"), parallax=_num0(d, "parallax"),
                rv=_num0(d, "radial_velocity"), epoch=_num(d, "ref_epoch"),
                ra_err=_num0(d, "ra_error"), dec_err=_num0(d, "dec_error"),
                pmra_err=_num0(d, "pmra_error"), pmdec_err=_num0(d, "pmdec_error"),
                parallax_err=_num0(d, "parallax_error"), mag=mag)
end

_gaia_table(release) = release == "dr3" ? "gaiadr3.gaia_source" :
                       release == "dr2" ? "gaiadr2.gaia_source" :
                       throw(ArgumentError("release must be \"dr3\" or \"dr2\""))

function Occultations.query_gaia(source_id::Union{Integer,AbstractString}; release::AbstractString="dr3")
    rows = _tap_csv("SELECT $(_GAIA_COLS) FROM $(_gaia_table(release)) WHERE source_id = $(source_id)")
    isempty(rows) && error("Gaia $(release) source $(source_id) not found")
    return _gaia_star(rows[1])
end

function Occultations.query_gaia(ra::Real, dec::Real; radius_arcsec::Real=2.0, release::AbstractString="dr3",
                                 maxmag::Real=Inf)
    q = "SELECT $(_GAIA_COLS) FROM $(_gaia_table(release)) WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), " *
        "CIRCLE('ICRS', $(ra), $(dec), $(radius_arcsec / 3600)))"
    isfinite(maxmag) && (q *= " AND phot_g_mean_mag < $(maxmag)")
    return [_gaia_star(d) for d in _tap_csv(q)]
end

# --------------------------------------------------------------------------- minimal JSON reader

# Parses the subset of JSON returned by the JPL APIs into Dict/Vector/String/Float64/Bool/nothing.
function _json(s::AbstractString)
    v, i = _jvalue(s, _jskip(s, firstindex(s)))
    return v
end

function _jskip(s, i)
    while i <= ncodeunits(s) && isspace(s[i])
        i = nextind(s, i)
    end
    return i
end

function _jvalue(s, i)
    c = s[i]
    if c == '{'
        d = Dict{String,Any}()
        i = _jskip(s, i + 1)
        s[i] == '}' && return d, i + 1
        while true
            k, i = _jstring(s, _jskip(s, i))
            i = _jskip(s, i)
            s[i] == ':' || error("JSON: ':' expected at $i")
            v, i = _jvalue(s, _jskip(s, i + 1))
            d[k] = v
            i = _jskip(s, i)
            s[i] == '}' && return d, i + 1
            s[i] == ',' || error("JSON: ',' expected at $i")
            i += 1
        end
    elseif c == '['
        a = Any[]
        i = _jskip(s, i + 1)
        s[i] == ']' && return a, i + 1
        while true
            v, i = _jvalue(s, _jskip(s, i))
            push!(a, v)
            i = _jskip(s, i)
            s[i] == ']' && return a, i + 1
            s[i] == ',' || error("JSON: ',' expected at $i")
            i += 1
        end
    elseif c == '"'
        return _jstring(s, i)
    elseif startswith(SubString(s, i), "true")
        return true, i + 4
    elseif startswith(SubString(s, i), "false")
        return false, i + 5
    elseif startswith(SubString(s, i), "null")
        return nothing, i + 4
    else
        m = match(r"^-?\d+(\.\d+)?([eE][+-]?\d+)?", SubString(s, i))
        m === nothing && error("JSON: unexpected character '$c' at $i")
        return parse(Float64, m.match), i + ncodeunits(m.match)
    end
end

function _jstring(s, i)
    s[i] == '"' || error("JSON: string expected at $i")
    io = IOBuffer()
    i = nextind(s, i)
    while s[i] != '"'
        if s[i] == '\\'
            i = nextind(s, i)
            e = s[i]
            if e == 'u'
                write(io, Char(parse(UInt16, s[i+1:i+4]; base=16)))
                i += 4
            else
                write(io, get(Dict('n' => '\n', 't' => '\t', 'r' => '\r', 'b' => '\b', 'f' => '\f'), e, e))
            end
        else
            write(io, s[i])
        end
        i = nextind(s, i)
    end
    return String(take!(io)), nextind(s, i)
end

# --------------------------------------------------------------------------- SBDB

const SBDB_URL = "https://ssd-api.jpl.nasa.gov/sbdb.api"

_sbdb_num(x) = x === nothing ? NaN : x isa Number ? Float64(x) : something(tryparse(Float64, strip(string(x))), NaN)

function Occultations.query_sbdb(target::AbstractString)
    r = HTTP.get(SBDB_URL; query=Dict("sstr" => target, "phys-par" => "1"), _timeout(60)..., status_exception=false)
    d = _json(String(r.body))
    haskey(d, "object") || error("SBDB: \"$target\" not found" * (haskey(d, "message") ? ": " * string(d["message"]) : ""))
    obj = d["object"]
    phys = Dict{String,Any}()
    for p in get(d, "phys_par", Any[])
        phys[p["name"]] = p["value"]
    end
    out = Dict{Symbol,Any}()
    v = _sbdb_num(get(phys, "diameter", nothing)); isnan(v) || (out[:radius] = v / 2)
    for (k, n, f) in ((:H, "H", 1.0), (:G, "G", 1.0), (:rotation_period, "rot_per", 1.0),
                      (:density, "density", 1000.0), (:GM, "GM", 1.0))
        v = _sbdb_num(get(phys, n, nothing))
        isnan(v) || (out[k] = f * v)
    end
    if haskey(phys, "pole") && phys["pole"] !== nothing
        w = tryparse.(Float64, split(replace(string(phys["pole"]), r"[^0-9.+\-]+" => " ")))
        length(w) >= 2 && all(!isnothing, w[1:2]) && (out[:pole] = (w[1], w[2]))
    end
    return NamedTuple(out)
end

# --------------------------------------------------------------------------- MPC observatory codes

const MPC_OBSCODES_URL = "https://minorplanetcenter.net/iau/lists/ObsCodes.html"
const _MPC_CODES = Ref{Vector{String}}(String[])

function Occultations.observer_from_mpc(code::AbstractString; name=nothing)
    if isempty(_MPC_CODES[])
        r = HTTP.get(MPC_OBSCODES_URL; _timeout(60)...)
        _MPC_CODES[] = String.(split(String(r.body), '\n'))
    end
    c = uppercase(strip(code))
    i = findfirst(l -> length(l) >= 4 && l[1:3] == c, _MPC_CODES[])
    i === nothing && error("MPC observatory code \"$code\" not found")
    line = _MPC_CODES[][i]
    length(strip(line[4:min(end, 30)])) == 0 && error("MPC code $code has no geocentric coordinates (space or roving observatory)")
    return parse_mpc_obscode(line; name=name)
end

# --------------------------------------------------------------------------- Horizons state vectors

function Occultations.spacecraft_horizons(id, tstart, tend; step::AbstractString="10m", name::AbstractString=string(id))
    t1, t2 = _epoch(tstart), _epoch(tend)
    # vector tables are in TDB: request the interval as TDB Julian dates
    fmt(t) = @sprintf("JD %.9f", sum(Occultations.tdb_jd(t)))
    query = Dict(
        "format" => "text", "COMMAND" => _hq(string(id)), "OBJ_DATA" => "'NO'", "MAKE_EPHEM" => "'YES'",
        "EPHEM_TYPE" => "'VECTORS'", "CENTER" => "'500@399'", "START_TIME" => _hq(fmt(t1)),
        "STOP_TIME" => _hq(fmt(t2)), "STEP_SIZE" => _hq(step), "REF_SYSTEM" => "'ICRF'",
        "REF_PLANE" => "'FRAME'", "VEC_CORR" => "'NONE'", "OUT_UNITS" => "'KM-S'", "VEC_TABLE" => "'1'",
        "CSV_FORMAT" => "'YES'")
    r = HTTP.get(HORIZONS_URL; query=query, _timeout(120)...)
    lines = split(String(r.body), '\n')
    i0 = findfirst(l -> startswith(strip(l), "\$\$SOE"), lines)
    i1 = findfirst(l -> startswith(strip(l), "\$\$EOE"), lines)
    (i0 === nothing || i1 === nothing) && error("Horizons did not return vectors for \"$id\"")
    times = Epoch[]; xyz = Float64[]
    for l in lines[i0+1:i1-1]
        c = strip.(split(l, ','))
        length(c) < 5 && continue
        push!(times, Epoch(parse(Float64, c[1]); scale=:tdb))
        append!(xyz, parse.(Float64, c[3:5]))
    end
    return Spacecraft(name, times, reshape(xyz, 3, :))
end

end # module
