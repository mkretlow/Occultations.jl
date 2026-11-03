# ---------------------------------------------------------------------------
# Observations of an event: contacts and chords
# ---------------------------------------------------------------------------

"""
    Contact(kind, time, σ; use=true, body=:auto)

A timed contact: `kind` is `:immersion` (disappearance) or `:emersion`, `time` an
[`Epoch`](@ref) (UTC), `σ` its 1σ uncertainty [s]. `use=false` keeps it out of fits. A ring
detection reported only by its mid-time is a contact of kind `:mid` with a ring label
(`Contact(:mid, t, σ; body=:ring1)`).
`body` assigns the contact to a component of a binary or contact-binary model:
`:primary`/`:secondary` ([`BinaryLimb`](@ref)), `:lobe1`/`:lobe2` (contact binaries), a ring
`:ring1` … `:ring9` (`:ring` = `:ring1`; [`RingedLimb`](@ref)), or `:auto` (default: matched to
the nearest crossing of the same kind of the whole outline).
"""
struct Contact
    kind::Symbol
    time::Epoch
    σ::Float64
    use::Bool
    body::Symbol
    function Contact(kind::Symbol, time, σ::Real; use::Bool=true, body::Symbol=:auto)
        kind in (:immersion, :emersion, :mid) || throw(ArgumentError("contact kind must be :immersion, :emersion or :mid"))
        kind === :mid && !_isring(body) && throw(ArgumentError("a :mid contact is a ring detection: give it a ring label"))
        σ >= 0 || throw(ArgumentError("contact uncertainty must be ≥ 0"))
        body in (:auto, :primary, :secondary, :lobe1, :lobe2) || _isring(body) ||
            throw(ArgumentError("contact body must be :auto, :primary, :secondary, :lobe1, :lobe2 or :ring1 … :ring9"))
        return new(kind, _epoch(time), Float64(σ), use, body)
    end
end

# whether a contact label is that of a ring (:ring, :ring1 … :ring9)
_isring(b::Symbol) = (s = string(b); s == "ring" || (length(s) == 5 && startswith(s, "ring") && s[5] in '1':'9'))

# a copy of the contact with some fields changed
_contact(x::Contact; kind=x.kind, time=x.time, σ=x.σ, use=x.use, body=x.body) = Contact(kind, time, σ; use=use, body=body)

Base.show(io::IO, c::Contact) = print(io, "Contact(", c.kind, ", ", isot(c.time), " ± ", c.σ, " s",
                                      c.body === :auto ? "" : ", $(c.body)", c.use ? "" : ", unused", ")")

"""
    Chord(name, observer; immersion=nothing, emersion=nothing, immersion_err=0, emersion_err=0,
          contacts=nothing, span=nothing, sensitivity=0.1, lightcurve=nothing)
    Chord(name, observer, timing::FitResult{TimingProblem}; span=timespan(lc))
    Chord(name, observer, lc::LightCurve)

One station's observation of an event.

* **positive**: one or two timed contacts (`immersion`, `emersion` as Epochs or ISO strings
  with their 1σ errors [s]); with only one contact it is a *partial* chord (recording started
  or ended during the event);
* **several drops** (a binary system, the two lobes of a contact binary): a vector of
  [`Contact`](@ref)s, e.g. `contacts = [Contact(:immersion, t1, σ), Contact(:emersion, t2, σ),
  Contact(:immersion, t3, σ), Contact(:emersion, t4, σ)]`, optionally assigned to a body
  (`Contact(...; body=:secondary)`); the star was seen between the drops;
* **negative**: no contact; `span = (start, stop)` is the observed interval;
* from a timing fit of a light curve ([`fit_timing`](@ref)), or from a light curve alone for
  joint photometric fits (`LimbProblem(...; photometry=true)`).

`span` is the time interval covered by the observation (default: the contacts, or the light
curve). `sensitivity` [s] is the timing resolution used for the soft constraint of a negative
chord (default: the exposure time of the light curve, else 0.1 s).
"""
struct Chord
    name::String
    observer::AbstractObserver
    contacts::Vector{Contact}
    span::Tuple{Epoch,Epoch}
    lightcurve::Union{Nothing,LightCurve}
    sensitivity::Float64
end

function Chord(name::AbstractString, obs::AbstractObserver; immersion=nothing, emersion=nothing,
               immersion_err::Real=0.0, emersion_err::Real=0.0, contacts=nothing, span=nothing, sensitivity=nothing,
               lightcurve::Union{Nothing,LightCurve}=nothing)
    cs = Contact[]
    immersion === nothing || push!(cs, Contact(:immersion, immersion, immersion_err))
    emersion === nothing || push!(cs, Contact(:emersion, emersion, emersion_err))
    if length(cs) == 2 && cs[2].time <= cs[1].time
        throw(ArgumentError("chord $name: emersion must follow immersion"))
    end
    if contacts !== nothing
        isempty(cs) || throw(ArgumentError("chord $name: give either `contacts` or `immersion`/`emersion`"))
        cs = sort!(collect(Contact, contacts); by=x -> x.time)
        ie = [x for x in cs if x.kind !== :mid]            # mid-times of ring dips are points, not edges
        for i in 2:length(ie)
            ie[i].kind === ie[i-1].kind &&
                throw(ArgumentError("chord $name: immersions and emersions must alternate in time"))
        end
    end
    sp = if span !== nothing
        (_epoch(span[1]), _epoch(span[2]))
    elseif lightcurve !== nothing
        timespan(lightcurve)
    elseif !isempty(cs)
        (minimum(c.time for c in cs), maximum(c.time for c in cs))
    else
        throw(ArgumentError("chord $name: a negative chord needs its observed `span`"))
    end
    sp[2] >= sp[1] || throw(ArgumentError("chord $name: span end before start"))
    sens = sensitivity === nothing ? (lightcurve === nothing ? 0.1 : max(exptime(lightcurve), 1e-3)) : Float64(sensitivity)
    return Chord(String(name), obs, cs, sp, lightcurve, sens)
end

function Chord(name::AbstractString, obs::AbstractObserver, r::FitResult{TimingProblem}; kw...)
    c = contact_times(r)
    lc = r.problem.lc
    return Chord(name, obs; immersion=c.immersion, emersion=c.emersion,
                 immersion_err=isnan(c.immersion_err) ? 0.0 : c.immersion_err,
                 emersion_err=isnan(c.emersion_err) ? 0.0 : c.emersion_err, lightcurve=lc, kw...)
end

Chord(name::AbstractString, obs::AbstractObserver, lc::LightCurve; kw...) = Chord(name, obs; lightcurve=lc, kw...)

"`true` if the chord has at least one contact."
ispositive(c::Chord) = !isempty(c.contacts)
"`true` if every immersion of the chord has its emersion (one or several drops, none partial)."
iscomplete(c::Chord) = (e = [x for x in c.contacts if x.kind !== :mid]; n = count(x -> x.kind === :immersion, e);
                        n >= 1 && n == count(x -> x.kind === :emersion, e) && first(e).kind === :immersion)

"""
    occulted_intervals(chord) -> Vector{Tuple{Epoch,Epoch}}

The intervals in which the star was occulted: one per drop, from its immersion to its
emersion (from the start or to the end of the span for a partial chord). Empty for a negative
chord.
"""
function occulted_intervals(c::Chord)
    out = Tuple{Epoch,Epoch}[]
    start = nothing
    for (i, x) in enumerate(c.contacts)
        x.kind === :mid && continue
        if x.kind === :immersion
            start = x.time
        else
            push!(out, (something(start, c.span[1]), x.time))
            start = nothing
        end
    end
    start === nothing || push!(out, (start, c.span[2]))
    return out
end
isnegative(c::Chord) = isempty(c.contacts)

"""
    disable(chord, kind=:both) -> Chord

Copy of the chord with its immersion and/or emersion (`:immersion`, `:emersion`, `:both`)
excluded from fits. The contacts are still plotted. See [`enable`](@ref).
"""
disable(c::Chord, kind::Symbol=:both) = _setuse(c, kind, false)

"""
    enable(chord, kind=:both) -> Chord

Copy of the chord with its immersion and/or emersion (`:immersion`, `:emersion`, `:both`)
included in fits again; the reverse of [`disable`](@ref).
"""
enable(c::Chord, kind::Symbol=:both) = _setuse(c, kind, true)
function _setuse(c::Chord, kind::Symbol, use::Bool)
    kind in (:both, :immersion, :emersion) || throw(ArgumentError("kind must be :both, :immersion or :emersion"))
    cs = [kind in (:both, x.kind) ? _contact(x; use=use) : x for x in c.contacts]
    return Chord(c.name, c.observer, cs, c.span, c.lightcurve, c.sensitivity)
end

"""
    shift(chord, dt) -> Chord

Copy of the chord with all its instants (contacts, span, light curve) shifted by `dt` [s],
e.g. to apply a clock correction.
"""
function shift(c::Chord, dt::Real)
    cs = [_contact(x; time=x.time + dt) for x in c.contacts]
    lc = c.lightcurve === nothing ? nothing :
         LightCurve(c.lightcurve.name, c.lightcurve.tref + dt, c.lightcurve.time, c.lightcurve.flux,
                    c.lightcurve.σ, c.lightcurve.optics)
    return Chord(c.name, c.observer, cs, (c.span[1] + dt, c.span[2] + dt), lc, c.sensitivity)
end

function Base.show(io::IO, c::Chord)
    nd = count(x -> x.kind === :immersion, c.contacts)
    kind = isnegative(c) ? "negative" : !iscomplete(c) ? "partial" : nd > 1 ? "positive, $nd drops" : "positive"
    print(io, "Chord(", c.name, ", ", kind, c.lightcurve === nothing ? "" : ", with light curve", ")")
end

function Base.show(io::IO, ::MIME"text/plain", c::Chord)
    println(io, c)
    println(io, "  observer: ", c.observer)
    println(io, "  span: ", isot(c.span[1]; digits=2), " – ", isot(c.span[2]; digits=2))
    for x in c.contacts
        println(io, "  ", x)
    end
end
