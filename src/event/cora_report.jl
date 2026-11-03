# ---------------------------------------------------------------------------
# Observation reports in the CORA format ("Key = value   {comment}" lines), as collected by
# the CORA portal and the OccultPortal: one observer at one site, with the result (O+, O−,
# O?, ND), the observed interval and, for a positive observation, the contact times.
# ---------------------------------------------------------------------------

"""
    CoraReport

One observation report in the CORA format, as read by [`read_cora_report`](@ref). Fields:
`date`, `object`, `star`, `result` (`:positive`, `:negative`, `:uncertain` or `:nodata`),
`observer` (names), `site` (an [`Observer`](@ref) named after the site), `optics`, `method`,
`timesource`, `start`, `stop` (the recorded interval), `immersion`, `emersion` (`nothing` if
not given) and their 1σ errors `immersion_err`, `emersion_err` [s], `source`, `comment`, and
`fields` (all key–value pairs of the file).
"""
struct CoraReport
    date::String
    object::String
    star::String
    result::Symbol
    observer::String
    site::Observer
    optics::String
    method::String
    timesource::String
    start::Union{Nothing,Epoch}
    stop::Union{Nothing,Epoch}
    immersion::Union{Nothing,Epoch}
    emersion::Union{Nothing,Epoch}
    immersion_err::Float64
    emersion_err::Float64
    source::String
    comment::String
    fields::Dict{String,String}
end

function Base.show(io::IO, r::CoraReport)
    print(io, "CoraReport(", r.date, ", ", r.object, ", ", r.site.name, ", ", r.result)
    r.immersion === nothing || print(io, ", D ", isot(r.immersion; digits=2), " ± ", r.immersion_err, " s")
    r.emersion === nothing || print(io, ", R ", isot(r.emersion; digits=2), " ± ", r.emersion_err, " s")
    print(io, ")")
end

const _CORA_RESULT = Dict("O+" => :positive, "O-" => :negative, "O−" => :negative, "O?" => :uncertain, "ND" => :nodata)

# "N 37:03:50.8896", "W 03 23 04.92", "-3.3847" → degrees (north and east positive)
function _cora_angle(s::AbstractString, pos::Char, neg::Char)
    t = strip(s)
    isempty(t) && throw(ArgumentError("empty coordinate"))
    sgn = 1.0
    if uppercase(t[1]) in (pos, neg)
        sgn = uppercase(t[1]) == neg ? -1.0 : 1.0
        t = strip(t[nextind(t, 1):end])
    end
    parts = split(t, r"[:\s]+")
    v = parse(Float64, parts[1])
    sgn *= startswith(parts[1], "-") ? -1.0 : 1.0
    v = abs(v)
    length(parts) > 1 && (v += parse(Float64, parts[2]) / 60)
    length(parts) > 2 && (v += parse(Float64, parts[3]) / 3600)
    return sgn * v
end

# the time of day `s` on `date`, the next day if it is earlier than `after` (midnight passed)
function _cora_time(date::AbstractString, s::AbstractString, after)
    t = strip(s)
    isempty(t) && return nothing
    e = Epoch(date * "T" * t)
    after !== nothing && e < after && (e = e + 86400.0)
    return e
end

"""
    read_cora_report(file) -> CoraReport

Read an observation report in the CORA format: lines `Key = value   {comment}` with the keys
`Date`, `Object`, `Star`, `Occultation` (O+, O−, O?, ND), `Observer`, `Sitename`, `Latitude`,
`Longitude` (`N`/`S`, `E`/`W` and `dd:mm:ss.ss` or `dd mm ss.ss`), `Altitude` [m], `Optics`,
`Method`, `Timesource`, `UT1`/`UT4` (start and end of the recording), `UT2`/`UT3` (disappearance
and reappearance) with `UT2ERR`/`UT3ERR` [s], `Source` and `Comment`. Times after midnight are
taken on the next day. The contact times of a report that is not positive are ignored, with a
warning. [`cora_chord`](@ref) converts the report into a [`Chord`](@ref).

# Examples
```jldoctest
julia> file = tempname();

julia> write(file, \"\"\"
       Date         = 2023-02-02                {yyyy-mm-dd}
       Object       = (447178) 2005 RO43        {object}
       Occultation  = O+                        {result}
       Sitename     = Test site                 {site}
       Latitude     = N 37:03:50.89             {latitude}
       Longitude    = W 03:23:04.92             {longitude}
       Altitude     = 2930                      {m}
       UT1          = 23:55:00                  {start}
       UT4          = 00:20:00                  {end}
       UT2          = 00:03:17.24               {D}
       UT3          = 00:03:54.32               {R}
       UT2ERR       = 0.2                       {s}
       UT3ERR       = 0.2                       {s}
       \"\"\");

julia> r = read_cora_report(file)
CoraReport(2023-02-02, (447178) 2005 RO43, Test site, positive, D 2023-02-03T00:03:17.24 ± 0.2 s, R 2023-02-03T00:03:54.32 ± 0.2 s)

julia> round(r.site.lon; digits=4), round(r.site.lat; digits=4)
(-3.3847, 37.0641)
```
"""
function read_cora_report(file::AbstractString)
    f = Dict{String,String}()
    for line in eachline(file)
        i = findfirst('=', line)
        i === nothing && continue
        key = strip(line[1:prevind(line, i)])
        val = line[nextind(line, i):end]
        j = findfirst('{', val)                            # the comment in braces
        j === nothing || (val = val[1:prevind(val, j)])
        f[key] = strip(val)
    end
    get(k) = Base.get(f, k, "")
    date = get("Date")
    isempty(date) && throw(ArgumentError("$file: no Date"))
    res = Base.get(_CORA_RESULT, uppercase(replace(get("Occultation"), " " => "")), nothing)
    res === nothing && throw(ArgumentError("$file: unknown Occultation result \"$(get("Occultation"))\""))
    name = isempty(get("Sitename")) ? basename(file) : get("Sitename")
    alt = isempty(get("Altitude")) ? 0.0 : parse(Float64, get("Altitude"))
    site = Observer(name, _cora_angle(get("Longitude"), 'E', 'W'), _cora_angle(get("Latitude"), 'N', 'S'), alt)
    start = _cora_time(date, get("UT1"), nothing)
    stop = _cora_time(date, get("UT4"), start)
    imm = _cora_time(date, get("UT2"), start)
    emm = _cora_time(date, get("UT3"), something(imm, start, Some(nothing)))
    err(k) = isempty(get(k)) ? 0.0 : parse(Float64, get(k))
    if res !== :positive && (imm !== nothing || emm !== nothing)
        @warn "CORA report $(basename(file)) ($name): result $(get("Occultation")), the contact times UT2/UT3 are ignored"
        imm = emm = nothing
    end
    return CoraReport(date, get("Object"), get("Star"), res, get("Observer"), site, get("Optics"), get("Method"),
                      get("Timesource"), start, stop, imm, emm, err("UT2ERR"), err("UT3ERR"), get("Source"),
                      get("Comment"), f)
end

"""
    cora_chord(report; sensitivity=nothing) -> Union{Chord,Nothing}

The [`Chord`](@ref) of a [`CoraReport`](@ref): a positive report gives the timed contacts
(UT2, UT3 with their errors) and the recorded interval (UT1–UT4) as the span; a negative
report a negative chord over UT1–UT4. A report with an uncertain result (O?) or without data
(ND) gives `nothing`. `sensitivity` [s] (e.g. the exposure time) is passed to the chord.

```julia
ev = Event(star, body, "2023-02-02T00:43")
for f in readdir("reports"; join=true)
    c = cora_chord(read_cora_report(f))
    c === nothing || push!(ev, c)
end
```
"""
function cora_chord(r::CoraReport; sensitivity=nothing)
    span = r.start === nothing || r.stop === nothing ? nothing : (r.start, r.stop)
    if r.result === :positive
        (r.immersion === nothing && r.emersion === nothing) &&
            throw(ArgumentError("CORA report $(r.site.name): positive, but no contact times"))
        return Chord(r.site.name, r.site; immersion=r.immersion, emersion=r.emersion, immersion_err=r.immersion_err,
                     emersion_err=r.emersion_err, span=span, sensitivity=sensitivity)
    elseif r.result === :negative
        span === nothing && throw(ArgumentError("CORA report $(r.site.name): negative, but no UT1/UT4 interval"))
        return Chord(r.site.name, r.site; span=span, sensitivity=sensitivity)
    end
    return nothing
end
