# ---------------------------------------------------------------------------
# An occultation event: star, body, reference epoch and the chords observed
# ---------------------------------------------------------------------------

"""
    Event(star, body, t; chords=Chord[], name="", window=3600)

One occultation of `star` by `body`. The reference epoch `tref` is the geocentric closest
approach searched within `t ± window` seconds (pass `window=0` to use `t` itself). Add
observations with `push!(event, chord)`. Per-chord [`Track`](@ref)s are computed when
first needed and cached.
"""
struct Event
    name::String
    star::Star
    body::Body
    tref::Epoch
    ca::NamedTuple
    chords::Vector{Chord}
    tracks::Dict{String,Track}
end

function Event(star::Star, body::Body, t; chords=Chord[], name::AbstractString="", window::Real=3600.0)
    t = _epoch(t)
    ca = window > 0 ? closest_approach(star, body, t - window, t + window; step=min(60.0, window / 10)) :
         (p = positionv(star, body, GEOCENTER, t); d = hypot(p.f, p.g);
          (epoch=t, ca_km=d, ca_arcsec=d / p.dist * RAD2ARCSEC, pa=mod(atand(-p.f, -p.g), 360.0),
           vel=hypot(p.vf, p.vg), dist_km=p.dist))
    nm = isempty(name) ? "$(body.name), $(isot(ca.epoch; digits=0))" : String(name)
    ev = Event(nm, star, body, ca.epoch, ca, Chord[], Dict{String,Track}())
    foreach(c -> push!(ev, c), chords)
    return ev
end

function Base.push!(ev::Event, c::Chord)
    any(x -> x.name == c.name, ev.chords) && throw(ArgumentError("event $(ev.name): chord $(c.name) already exists"))
    push!(ev.chords, c)
    return ev
end

function Base.getindex(ev::Event, name::AbstractString)
    i = findfirst(c -> c.name == name, ev.chords)
    i === nothing && throw(KeyError(name))
    return ev.chords[i]
end

"Replace the chord of the same name (e.g. by `disable(chord)` or `shift(chord, dt)`)."
function replace_chord!(ev::Event, c::Chord)
    i = findfirst(x -> x.name == c.name, ev.chords)
    i === nothing && throw(KeyError(c.name))
    ev.chords[i] = c
    delete!(ev.tracks, c.name)
    return ev
end

"Remove a chord."
function Base.delete!(ev::Event, name::AbstractString)
    ev[name]
    filter!(c -> c.name != name, ev.chords)
    delete!(ev.tracks, name)
    return ev
end

"""
    track(event, chord; margin=120) -> Track

Track of the chord's observer covering its observed span plus `margin` seconds on each
side (room for fitted clock offsets); cached in the event.
"""
function track(ev::Event, c::Chord; margin::Real=120.0)
    tr = get(ev.tracks, c.name, nothing)
    a, b = c.span[1] - ev.tref - margin, c.span[2] - ev.tref + margin
    if tr === nothing || tr.τmin > a || tr.τmax < b || tr.observer !== c.observer
        tr = Track(ev.star, ev.body, c.observer, ev.tref, a, b)
        ev.tracks[c.name] = tr
    end
    return tr
end
track(ev::Event, name::AbstractString; kw...) = track(ev, ev[name]; kw...)

"""
    contact_points(event; all=false) -> Vector{NamedTuple}

Positions of the timed contacts in the fundamental plane: `(chord, kind, epoch, f, g, vf,
vg, σ_t, σ, body)` with `σ = σ_t |v|` [km] and the `body` label of the contact. Contacts with `use=false` are left out unless `all`.
"""
function contact_points(ev::Event; all::Bool=false)
    out = NamedTuple[]
    for c in ev.chords, x in c.contacts
        (x.use || all) || continue
        tr = track(ev, c)
        τ = x.time - ev.tref
        f, g = tr(τ)
        vf, vg = velocity(tr, τ)
        push!(out, (chord=c.name, kind=x.kind, epoch=x.time, f=f, g=g, vf=vf, vg=vg, σ_t=x.σ, σ=x.σ * hypot(vf, vg),
                    body=x.body))
    end
    return out
end

"""
    chord_path(event, chord; n=200, segment=:observed) -> (f, g)

(f, g) path of a chord over its observed span (`segment=:observed`) or between its contacts
(`:occulted`; for a partial chord from/to the end of the span; with several drops, one
piece per drop, separated by `NaN`).
"""
function chord_path(ev::Event, c::Chord; n::Integer=200, segment::Symbol=:observed, dt::Real=0.0)
    tr = track(ev, c)
    pieces = segment === :occulted ? occulted_intervals(c) : [c.span]
    fs, gs = Float64[], Float64[]
    for (k, (a, b)) in enumerate(pieces)
        k > 1 && (push!(fs, NaN); push!(gs, NaN))
        p = tr.(range(a - ev.tref + dt, b - ev.tref + dt; length=n))
        append!(fs, first.(p)); append!(gs, last.(p))
    end
    return fs, gs
end
chord_path(ev::Event, name::AbstractString; kw...) = chord_path(ev, ev[name]; kw...)

function Base.show(io::IO, ev::Event)
    print(io, "Event(", ev.name, ", ", length(ev.chords), " chords)")
end

function Base.show(io::IO, ::MIME"text/plain", ev::Event)
    println(io, "Event ", ev.name, ": ", ev.body.name, " occults ", ev.star)
    @printf(io, "  closest approach %s: %.4f\" (%.1f km), PA %.2f°, %.3f km/s, distance %.4f au\n",
            isot(ev.tref; digits=3), ev.ca.ca_arcsec, ev.ca.ca_km, ev.ca.pa, ev.ca.vel, ev.ca.dist_km / AU_KM)
    for c in ev.chords
        println(io, "  ", c)
    end
end

"""
    Occultation

Alias of [`Event`](@ref).
"""
const Occultation = Event

# ---------------------------------------------------------------------------
# geometry of an event used by the limb models
# ---------------------------------------------------------------------------

"""
    EventContext

Geometry of an event at its reference epoch: the geocentric line of sight `u` (unit vector
observer → body, ICRS), the distance [km], the light time [s], and the TDB Julian date of
light emission (for body rotation models).
"""
struct EventContext
    tref::Epoch
    u::Vector{Float64}
    dist::Float64
    light_time::Float64
    jd_emit::Float64
end

function EventContext(ev::Event)
    B = geocentric_vector(ev.body.ephem, ev.tref)
    D = norm(B)
    lt = D / C_KMS
    t1, t2 = tdb_jd(ev.tref - lt)
    return EventContext(ev.tref, B ./ D, D, lt, t1 + t2)
end
