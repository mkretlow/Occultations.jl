# ---------------------------------------------------------------------------
# The limb-fitting problem: one or several events, any limb model, timed contacts,
# negative chords, light curves (joint photometric fit, with or without an atmosphere),
# clock offsets, rotational-amplitude constraints and priors
# ---------------------------------------------------------------------------

"""
    RotationalAmplitude(epoch, amplitude, σ)

Peak-to-peak amplitude [mag] (± σ) of the rotational light curve observed at `epoch`, a
constraint on the shape and pole of 3-D limb models (geometric scattering: brightness ∝
projected area). For an ellipsoid this is the classical
`Δm = 2.5 log(a/b) - 1.25 log[(a² cos²θ + c² sin²θ)/(b² cos²θ + c² sin²θ)]`
with the aspect angle θ following from the pole and the viewing direction at `epoch`.
"""
struct RotationalAmplitude
    epoch::Epoch
    amplitude::Float64
    σ::Float64
end
RotationalAmplitude(epoch, amplitude::Real, σ::Real) = RotationalAmplitude(_epoch(epoch), Float64(amplitude), Float64(σ))

"""
    RotationalPhase(epoch, phase, σ)

Rotational phase at `epoch` measured from the preceding **maximum of light**, in rotations
(e.g. 0.63 = 0.63 rotations after the maximum), ± σ, from a rotational light curve; a
constraint on the rotation phase (`W0`) of 3-D limb models. Brightness is taken ∝ projected
area (geometric scattering): the model phase is the time since the last maximum of the
projected area. For an ellipsoid the two maxima per rotation are equal, so the phase is
compared modulo ½ rotation.
"""
struct RotationalPhase
    epoch::Epoch
    phase::Float64
    σ::Float64
end
RotationalPhase(epoch, phase::Real, σ::Real) = RotationalPhase(_epoch(epoch), Float64(phase), Float64(σ))

struct _PhaseData
    u::Vector{Float64}
    jd_emit::Float64
    a::RotationalPhase
end

struct _ContactData
    k::Int
    tr::Track
    τ::Float64
    entry::Bool
    σt::Float64
    dt::Int
    chord::String
    kind::Symbol
    body::Int                      # component of a composite outline (0: nearest crossing)
end

struct _NegativeData
    k::Int
    tr::Track
    τ::Vector{Float64}
    xy::Matrix{Float64}            # precomputed track positions (2 × n) when no clock offset is fitted
    σ::Float64
    dt::Int
    chord::String
end

struct _PhotData
    k::Int
    tr::Track
    chord::String
    lc::LightCurve
    off::Float64                   # τ = off + t_lc (+ dt)
    dt::Int
    star_diam::Float64
    profile::EdgeProfile
end

struct _AmpData
    u::Vector{Float64}
    a::RotationalAmplitude
end

"""
    LimbProblem(events; model=EllipseLimb(), residual=:chord, negative=:hard, visible=:auto,
                σ_limb=0, σ_ring=σ_limb, fit_dt=(), photometry=false, atmosphere=nothing, amplitudes=(), phases=(),
                fix=(), init=Dict(), priors=Dict(), bounds=Dict())

Fit of the apparent limb to one event or a vector of events (multi-event fit).

**Parameters**: the centre of each event (`center_f`, `center_g` [km] in the fundamental
plane relative to the ephemeris; suffix `_k` for event `k` in multi-event fits), the
parameters of the `model` ([`LimbModel`](@ref)), a clock offset `dt_<chord>` [s] for every
chord in `fit_dt`, and `pressure` [nbar] with an `atmosphere`. Change their settings with
`fix`, `init`, `priors` (`name => Gaussian(μ, σ)` or `(μ, σ)`), `bounds` (`name => (lo, hi)`)
or [`setparam!`](@ref).

**Data**:
* timed contacts of the chords (`residual = :chord`: distance along the chord between the
  observed and the predicted contact, i.e. timing residual × speed, which is the
  measurement model and well behaved for grazing chords; a chord that misses the limb is
  penalised by the distance by which it misses. `:radial`: radial distance of the contact
  point from the limb, as in SORA). Uncertainty per contact: `√((σ_t |v|)² + σ_limb²)`,
  `σ_limb` [km] being an allowance for topography; ring contacts take `σ_ring` instead (default:
  `σ_limb`), an allowance for the ring's eccentricity or width variations;
* negative chords (`negative = :hard`, the default: limbs crossing the track of a negative
  chord are excluded; `:soft`: a penalty proportional to the depth by which the limb crosses
  the track, in units of `sensitivity × speed`, for chords whose detection limit is
  uncertain, e.g. long exposures or low signal-to-noise; `:none`);
* the **visible segments** of the positive chords: where the star was seen before, between
  and after the drops (shortened at each contact by 3σ of its timing, at least the chord's
  `sensitivity`). They are treated as negative chords: `visible = :soft` (penalty as for
  negative chords), `:hard` or `:none`. The default `:auto` is `:soft` for binary and
  contact-binary models ([`BinaryLimb`](@ref), [`ContactBinaryLimb`](@ref),
  [`ContactBinaryEllipsoids`](@ref)), whose outline can produce a drop where none was seen,
  and `:none` otherwise (a convex limb cannot);
* light curves (`photometry = true`, or a list of chord names): the light curve itself is
  fitted with the model flux predicted by the limb — diffraction, stellar diameter,
  passband, exposure (contacts of these chords are not used separately), or the
  refractive flux of the `atmosphere` ([`Atmosphere`](@ref)) whose surface pressure is
  then fitted; the flux levels are profiled per light curve;
* rotational light-curve `amplitudes` ([`RotationalAmplitude`](@ref)) and `phases`
  ([`RotationalPhase`](@ref)) for 3-D models.

Solve with [`fit`](@ref), sample the posterior with [`sample`](@ref).
"""
struct LimbProblem{M<:LimbModel} <: AbstractProblem
    setup::_Setup
    model::M
    ps::Vector{Param}
    contacts::Vector{_ContactData}
    negatives::Vector{_NegativeData}
    visibles::Vector{_NegativeData}
    photometry::Vector{_PhotData}
    amplitudes::Vector{_AmpData}
    phases::Vector{_PhaseData}
    residual::Symbol
    negative::Symbol
    visible::Symbol
    σ_limb::Float64
    σ_ring::Float64
    atmosphere::Union{Nothing,Atmosphere}
    centers::Vector{Tuple{Int,Int}}
    pressure::Int
    atm_cache::Dict{NTuple{5,Float64},AtmosphereProfile}
end

parameters(p::LimbProblem) = p.ps
events(p::LimbProblem) = p.setup.events

LimbProblem(ev::Event; kw...) = LimbProblem([ev]; kw...)

function LimbProblem(evs::AbstractVector{Event}; model::LimbModel=EllipseLimb(), residual::Symbol=:chord,
                     negative::Symbol=:hard, visible::Symbol=:auto, σ_limb::Real=0.0, σ_ring::Real=NaN, fit_dt=(), photometry=false,
                     atmosphere::Union{Nothing,Atmosphere}=nothing, amplitudes=(), phases=(), fix=(), init=Dict(),
                     priors=Dict(), bounds=Dict(), epoch=nothing)
    residual in (:chord, :radial) || throw(ArgumentError("residual must be :chord or :radial"))
    negative in (:soft, :hard, :none) || throw(ArgumentError("negative must be :soft, :hard or :none"))
    visible in (:auto, :soft, :hard, :none) || throw(ArgumentError("visible must be :auto, :soft, :hard or :none"))
    visible === :auto && (visible = _composite(model) ? :soft : :none)
    evs = collect(Event, evs)
    nev = length(evs)
    nev >= 1 || throw(ArgumentError("no events"))
    ctxs = [EventContext(e) for e in evs]
    photchords = photometry === true ? Set(c.name for e in evs for c in e.chords if c.lightcurve !== nothing) :
                 photometry === false ? Set{String}() : Set(String.(collect(photometry)))
    # starting estimates of each event from its contact points
    ests = map(evs) do ev
        pts = filter(q -> !_isring(q.body), contact_points(ev))      # the body only, not its rings
        f = [p.f for p in pts]; g = [p.g for p in pts]
        R0 = isfinite(ev.body.radius) ? ev.body.radius : 100.0
        if length(f) >= 3
            c = algebraic_circle(f, g)
            ae = algebraic_ellipse(f, g)
            ae !== nothing && !(0 <= ae.oblateness < 0.9) && (ae = nothing)
            (center_f=c.center_f, center_g=c.center_g, radius=c.radius, ellipse=ae)
        elseif !isempty(f)
            (center_f=mean(f), center_g=mean(g), radius=R0, ellipse=nothing)
        else
            (center_f=0.0, center_g=0.0, radius=R0, ellipse=nothing)
        end
    end
    ref = epoch !== nothing ? _epoch(epoch) :
          (hasproperty(model, :epoch) && model.epoch !== nothing) ? model.epoch : evs[1].tref
    setup = _Setup(evs, ctxs, ests, Dict{Symbol,Int}(), ref)
    _prepare!(model, setup)
    ests = setup.estimates                  # possibly updated by the model (e.g. a binary: the primary)
    ps = Param[]
    for k in 1:nev
        sc = max(ests[k].radius / 10, 1.0)
        push!(ps, Param(_pn(:center_f, k, nev), ests[k].center_f; scale=sc),
              Param(_pn(:center_g, k, nev), ests[k].center_g; scale=sc))
    end
    append!(ps, _model_params(model, setup))
    # clock offsets
    dtnames = String.(collect(fit_dt))
    dtidx = Dict{Tuple{Int,String},Symbol}()
    for n in dtnames
        found = false
        for (k, e) in enumerate(evs), c in e.chords
            c.name == n || continue
            pname = Symbol("dt_", n, nev > 1 && count(x -> any(c -> c.name == n, x.chords), evs) > 1 ? "_$k" : "")
            push!(ps, Param(pname, 0.0; scale=0.2))
            dtidx[(k, n)] = pname
            found = true
        end
        found || throw(ArgumentError("fit_dt: no chord named $n"))
    end
    if atmosphere !== nothing
        _composite(model) && throw(ArgumentError("atmospheres are not available for binary or contact-binary models"))
        isempty(photchords) && throw(ArgumentError("an atmosphere needs light curves: set `photometry`"))
        push!(ps, Param(:pressure, 100.0; lower=0.0, scale=20.0))
    end
    _configure!(ps; fix=fix, init=init, priors=priors, bounds=bounds)
    for (i, p) in enumerate(ps)
        setup.index[p.name] = i
    end
    _finalize!(model, setup)
    # data
    contacts = _ContactData[]; negs = _NegativeData[]; vis = _NegativeData[]; phots = _PhotData[]
    for (k, ev) in enumerate(evs), c in ev.chords
        tr = track(ev, c)
        di = haskey(dtidx, (k, c.name)) ? setup.index[dtidx[(k, c.name)]] : 0
        if c.name in photchords
            c.lightcurve === nothing && throw(ArgumentError("chord $(c.name) has no light curve"))
            lc = c.lightcurve
            off = lc.tref - ev.tref
            tr = track(ev, c; margin=120.0)
            sd = apparent_diameter(ev.star, tr.dist)
            push!(phots, _PhotData(k, tr, c.name, lc, off, di, sd, _cached_profile(tr.dist, sd, lc.optics.λ, lc.optics.Δλ)))
            continue
        end
        # a chord that saw only ring dips missed the body: a negative chord for the body
        if !isnegative(c) && negative !== :none && all(x -> _isring(x.body), c.contacts)
            τs = collect(range(c.span[1] - ev.tref, c.span[2] - ev.tref; length=400))
            xy = reduce(hcat, [collect(tr(τ)) for τ in τs])
            v = hypot(velocity(tr, τs[200])...)
            push!(negs, _NegativeData(k, tr, τs, xy, c.sensitivity * v, di, c.name))
        end
        if isnegative(c)
            negative === :none && continue
            τs = collect(range(c.span[1] - ev.tref, c.span[2] - ev.tref; length=400))
            xy = reduce(hcat, [collect(tr(τ)) for τ in τs])
            v = hypot(velocity(tr, τs[200])...)
            push!(negs, _NegativeData(k, tr, τs, xy, c.sensitivity * v, di, c.name))
        else
            for x in c.contacts
                x.use || continue
                push!(contacts, _ContactData(k, tr, x.time - ev.tref, x.kind === :immersion, x.σ, di, c.name, x.kind,
                                             _body_index(model, x.body)))
            end
            visible === :none || append!(vis, _visible_segments(ev, c, tr, k, di))
        end
    end
    # the ephemeris of the event closest in time (each event may carry its own short table)
    eph(t) = evs[argmin([abs(t - e.tref) for e in evs])].body.ephem
    amps = _AmpData[_AmpData(normalize(geocentric_vector(eph(a.epoch), a.epoch)), a) for a in amplitudes]
    phs = map(collect(phases)) do a
        B = geocentric_vector(eph(a.epoch), a.epoch)
        t1, t2 = tdb_jd(a.epoch - norm(B) / C_KMS)
        _PhaseData(normalize(B), t1 + t2, a)
    end
    phs = convert(Vector{_PhaseData}, phs)
    (isempty(amps) && isempty(phs)) || hasmethod(_shape, Tuple{typeof(model),Vector{Float64},_Setup}) ||
        throw(ArgumentError("rotational amplitudes and phases need a 3-D limb model (EllipsoidLimb or MeshLimb)"))
    centers = [(setup.index[_pn(:center_f, k, nev)], setup.index[_pn(:center_g, k, nev)]) for k in 1:nev]
    isempty(contacts) && isempty(phots) && throw(ArgumentError("no timed contacts or light curves to fit"))
    return LimbProblem(setup, model, ps, contacts, negs, vis, phots, amps, phs, residual, negative, visible, Float64(σ_limb),
                       Float64(isnan(σ_ring) ? σ_limb : σ_ring),
                       atmosphere, centers, atmosphere === nothing ? 0 : setup.index[:pressure],
                       Dict{NTuple{5,Float64},AtmosphereProfile}())
end

# outlines of all events for the values v
_outlines(p::LimbProblem, v) = [_outline(p.model, v, k, p.setup) for k in 1:length(p.setup.events)]

function residuals(p::LimbProblem, v)
    _admissible(p.model, v, p.setup) || return nothing
    O = _outlines(p, v)
    R = [_rings(p.model, v, k, p.setup) for k in 1:length(p.setup.events)]
    r = eltype(v)[]
    sizehint!(r, 2length(p.contacts) + length(p.negatives) + sum(d -> length(d.lc), p.photometry; init=0))
    # contacts
    for d in p.contacts
        O_ = O[d.k]
        cf, cg = v[p.centers[d.k][1]], v[p.centers[d.k][2]]
        τ = d.τ + (d.dt == 0 ? 0.0 : v[d.dt])
        f, g = d.tr(τ)
        vf, vg = velocity(d.tr, τ)
        x, y = f - cf, g - cg
        sp = hypot(vf, vg)
        σ = sqrt((d.σt * sp)^2 + (d.body < 0 ? p.σ_ring : p.σ_limb)^2)
        σ > 0 || throw(ArgumentError("contact of chord $(d.chord) has zero uncertainty: set its σ or σ_limb"))
        if d.body < 0                          # a ring: along the chord, nearest crossing of its annulus
            s, miss = d.kind === :mid ? _ring_mid(R[d.k][-d.body], x, y, vf, vg) :
                      _ring_crossing(R[d.k][-d.body], x, y, vf, vg, d.entry)
            push!(r, -s * sp / σ, miss / σ)
        elseif p.residual === :radial
            push!(r, _part_distance(O_, x, y, d.body) / σ)
        else
            s, miss = _contact_crossing(O_, x, y, vf, vg, d.entry, d.body)
            push!(r, -s * sp / σ, miss / σ)
        end
    end
    # negative chords, and the visible segments of the positive chords
    for (mode, ds) in ((p.negative, p.negatives), (p.visible, p.visibles))
        for d in ds
            dmin = _negative_depth(p, d, O, v)
            if mode === :hard
                dmin < 0 && return nothing
            else
                push!(r, max(-dmin, 0.0) / max(d.σ, 1e-6))
            end
        end
    end
    # light curves
    for d in p.photometry
        m = _phot_model(p, d, O[d.k], v, R[d.k])
        m === nothing && return nothing
        I, w = m
        lc = d.lc
        ft, fb = _flux_levels(I, lc.flux, w, 1.0, 0.0, true, p.atmosphere === nothing && count(<(0.1), I) >= 3)
        append!(r, @. (lc.flux - fb - (ft - fb) * I) / lc.σ)
    end
    # rotational amplitudes and phases
    if !isempty(p.amplitudes) || !isempty(p.phases)
        sh = _shape(p.model, v, p.setup)
        α, δ = v[p.setup.index[:pole_ra]], v[p.setup.index[:pole_dec]]
        for d in p.amplitudes
            push!(r, (_amplitude(sh, α, δ, d.u) - d.a.amplitude) / d.a.σ)
        end
        for d in p.phases
            per = sh isa Union{Ellipsoid,BilobedEllipsoid} ? 0.5 : 1.0
            f = _phase_since_maximum(sh, α, δ, _meridian_at(v, p.setup, d.jd_emit), v[p.setup.index[:period]], d.u)
            Δ = mod(f - d.a.phase + per / 2, per) - per / 2
            push!(r, Δ / d.a.σ)
        end
    end
    return r
end

# smallest signed distance of a negative chord's track from the limb (< 0: the limb crosses it)
function _negative_depth(p::LimbProblem, d, O, v)
    O_ = O[d.k]
    cf, cg = v[p.centers[d.k][1]], v[p.centers[d.k][2]]
    dmin = Inf
    if d.dt == 0
        @inbounds for j in axes(d.xy, 2)
            dmin = min(dmin, signed_distance(O_, d.xy[1, j] - cf, d.xy[2, j] - cg))
        end
    else
        for τ in d.τ
            f, g = d.tr(τ + v[d.dt])
            dmin = min(dmin, signed_distance(O_, f - cf, g - cg))
        end
    end
    return dmin
end

# the visible segments of a positive chord (before, between and after its drops), shortened
# at every contact by 3σ of its timing (at least the chord's sensitivity), as negative-chord
# data; σ of the soft penalty: sensitivity × speed
function _visible_segments(ev::Event, c::Chord, tr::Track, k::Integer, di::Integer; n::Integer=200)
    out = _NegativeData[]
    cs = [x for x in c.contacts if x.kind !== :mid]
    isempty(cs) && return out
    margin(x) = max(3x.σ, c.sensitivity)
    pieces = Tuple{Float64,Float64}[]
    t0 = c.span[1] - ev.tref
    first(cs).kind === :immersion && push!(pieces, (t0, cs[1].time - ev.tref - margin(cs[1])))
    for i in 1:length(cs)-1
        cs[i].kind === :emersion && cs[i+1].kind === :immersion || continue
        push!(pieces, (cs[i].time - ev.tref + margin(cs[i]), cs[i+1].time - ev.tref - margin(cs[i+1])))
    end
    last(cs).kind === :emersion && push!(pieces, (cs[end].time - ev.tref + margin(cs[end]), c.span[2] - ev.tref))
    for (j, (a, b)) in enumerate(pieces)
        b > a || continue
        τs = collect(range(a, b; length=n))
        xy = reduce(hcat, [collect(tr(τ)) for τ in τs])
        v = hypot(velocity(tr, (a + b) / 2)...)
        push!(out, _NegativeData(k, tr, τs, xy, c.sensitivity * v, di, "$(c.name) (visible $j)"))
    end
    return out
end

# the parts of the interval [a, b] outside the intervals `cut` (sorted, disjoint)
function _subtract(a, b, cut)
    out = Tuple{typeof(a),typeof(a)}[]
    lo = a
    for (c, d) in cut
        d <= lo && continue
        c >= b && break
        c > lo && push!(out, (lo, min(c, b)))
        lo = max(lo, d)
    end
    lo < b && push!(out, (lo, b))
    return out
end

# whether the parameters `v` satisfy the hard constraints (model admissibility and, with
# `negative = :hard`, no limb crossing a negative chord) – cheap, without the light curves
function _satisfies_constraints(p::LimbProblem, v)
    _admissible(p.model, v, p.setup) || return false
    (p.negative === :hard || p.visible === :hard) || return true
    O = _outlines(p, v)
    return (p.negative !== :hard || all(d -> _negative_depth(p, d, O, v) >= 0, p.negatives)) &&
           (p.visible !== :hard || all(d -> _negative_depth(p, d, O, v) >= 0, p.visibles))
end

# prime meridian [deg] at the TDB Julian date of light emission `jd`
function _meridian_at(v, s::_Setup, jd)
    t1, t2 = tdb_jd(s.epoch)
    return _v(v, s, :W0) + 360 * 24 / _v(v, s, :period) * (jd - (t1 + t2))
end

"""
    rotational_phase(shape, pole_ra, pole_dec, W, period, u; n=120) -> rotations

Time since the last maximum of the projected area (geometric-scattering light curve) of a
shape with pole (`pole_ra`, `pole_dec`) and prime meridian `W` [deg] at the epoch, seen
along `u`, in rotations (`period` [h] only sets the sense of rotation).
"""
rotational_phase(sh::AbstractShape, α, δ, W, P, u::AbstractVector; n::Integer=120) =
    _phase_since_maximum(sh, α, δ, W, P, u; n=n)

function _phase_since_maximum(sh, α, δ, W, P, u; n::Integer=120)
    sgn = P < 0 ? -1.0 : 1.0
    A(x) = area(project(sh, _body_matrix(α, δ, W - sgn * 360 * x), u; n=180))     # x rotations ago
    xs = [(k - 1) / n for k in 1:n]
    as = A.(xs)
    k = argmax(as)
    a0, am, ap = as[k], as[mod1(k - 1, n)], as[mod1(k + 1, n)]
    den = am - 2a0 + ap
    δx = den < 0 ? 0.5 * (am - ap) / den : 0.0
    return mod(xs[k] + δx / n, 1.0)
end

function _amplitude(e::Ellipsoid, α, δ, u)
    a, b = max(e.a, e.b), min(e.a, e.b)
    return rotational_amplitude(a, b, e.c, aspect_angle(α, δ, u))
end
_amplitude(s::AbstractShape, α, δ, u) = rotational_amplitude(s, α, δ, u)

nobs(p::LimbProblem) = length(p.contacts) + sum(d -> length(d.lc) - 2, p.photometry; init=0) + length(p.amplitudes) +
                       length(p.phases)

# model intensity of a light curve for the outline O; returns (I, weights)
function _phot_model(p::LimbProblem, d::_PhotData, O::Outline, v, rings=_rings(p.model, v, d.k, p.setup))
    lc = d.lc
    cf, cg = v[p.centers[d.k][1]], v[p.centers[d.k][2]]
    dt = d.dt == 0 ? 0.0 : v[d.dt]
    w = 1 ./ lc.σ .^ 2
    if p.atmosphere !== nothing
        R = equivalent_radius(O)
        v[p.pressure] >= 0 || return nothing
        key = (Float64(v[p.pressure]), Float64(R), d.tr.dist, lc.optics.λ, lc.optics.Δλ)
        prof = get(p.atm_cache, key, nothing)
        if prof === nothing
            length(p.atm_cache) > 16 && empty!(p.atm_cache)
            prof = atmosphere_profile(p.atmosphere, key[1], R, body_gm(p.setup.events[d.k].body), d.tr.dist;
                                      λ=lc.optics.λ, Δλ=lc.optics.Δλ, star_diam=d.star_diam)
            p.atm_cache[key] = prof
        end
        texp = exptime(lc)
        # sub-sample the exposure finely enough to resolve the diffraction at the surface edge
        speed = hypot(velocity(d.tr, d.off + lc.time[end ÷ 2 + 1])...)
        F = fresnel_scale(lc.optics.λ, d.tr.dist)
        ns = texp > 0 ? clamp(ceil(Int, speed * texp / (F / 3)), 3, 128) : 1
        I = map(lc.time) do t
            acc = 0.0
            for j in 1:ns
                τ = d.off + t + dt + (ns == 1 ? 0.0 : texp * ((j - 0.5) / ns - 0.5))
                f, g = d.tr(τ)
                x, y = f - cf, g - cg
                ρ = hypot(x, y) * R / limb_radius(O, x, y)        # elliptical distance scaled to R
                acc += atmosphere_flux(prof, ρ)
            end
            acc / ns
        end
        return I, w
    end
    # several bodies or lobes, or rings: the occulted intervals along the (locally straight) track,
    # each an opaque (body) or partly transparent (ring) strip
    if O isa CompositeOutline || !isempty(rings)
        tm = d.off + lc.time[(end + 1) ÷ 2] + dt
        f0, g0 = d.tr(tm)
        vf, vg = velocity(d.tr, tm)
        x0, y0 = f0 - cf, g0 - cg
        vnorm(o, s) = ((nx, ny) = limb_normal(o, x0 + s * vf, y0 + s * vg); max(abs(vf * nx + vg * ny), 1e-3))
        T = promote_type(eltype(v), Float64)
        strips = Tuple{T,T,T,T,T}[]                        # start, end, normal velocities, opacity
        body = intervals(O, x0, y0, vf, vg)
        for (a, b) in body
            push!(strips, (a, b, vnorm(O, a), vnorm(O, b), one(T)))
        end
        for rg in rings
            edge, _ = _ring_edges(rg)
            for (a, b) in intervals(rg, x0, y0, vf, vg), (a2, b2) in _subtract(a, b, body)
                push!(strips, (a2, b2, vnorm(edge, a2), vnorm(edge, b2), rg.opacity))
            end
        end
        isempty(strips) && return ones(length(lc.time)), w
        sort!(strips; by=first)
        t0 = tm - d.off - dt
        I = _multi_intensity(lc.time, [t0 + q[1] for q in strips], [t0 + q[2] for q in strips], [q[3] for q in strips],
                             [q[4] for q in strips], exptime(lc), d.profile, d.tr.dist, d.star_diam, lc.optics.λ,
                             lc.optics.Δλ, [q[5] for q in strips])
        return I, w
    end
    # contacts predicted by the limb along the light curve
    τs = d.off .+ lc.time .+ dt
    h(τ) = ((f, g) = d.tr(τ); signed_distance(O, f - cf, g - cg))
    hs = h.(τs)
    n = length(τs)
    ie = findfirst(i -> hs[i] >= 0 && hs[i+1] < 0, 1:n-1)
    ix = findlast(i -> hs[i] < 0 && hs[i+1] >= 0, 1:n-1)
    function root(i)
        lo, hi = τs[i], τs[i+1]
        hlo = hs[i]
        for _ in 1:40
            m = (lo + hi) / 2
            hm = h(m)
            (hm >= 0) == (hlo >= 0) ? ((lo, hlo) = (m, hm)) : (hi = m)
        end
        return (lo + hi) / 2
    end
    function vn(τ)
        f, g = d.tr(τ)
        nx, ny = limb_normal(O, f - cf, g - cg)
        vf, vg = velocity(d.tr, τ)
        return max(abs(vf * nx + vg * ny), 1e-3)
    end
    τi = ie === nothing ? (hs[1] < 0 ? -Inf : NaN) : root(ie)
    τe = ix === nothing ? (hs[end] < 0 ? Inf : NaN) : root(ix)
    if isnan(τi) || isnan(τe)                  # no event along this light curve
        return ones(n), w
    end
    vi = isfinite(τi) ? vn(τi) : vn(τe)
    ve = isfinite(τe) ? vn(τe) : vi
    ti, te = τi - d.off - dt, τe - d.off - dt
    I = _model_intensity(lc.time, ti, te, vi, ve, exptime(lc), d.profile, d.tr.dist, d.star_diam,
                         lc.optics.λ, lc.optics.Δλ, 1.0)
    return I, w
end

function derived(p::LimbProblem, v)
    s = p.setup
    nev = length(s.events)
    d = Pair{Symbol,Float64}[]
    for k in 1:nev
        D = s.contexts[k].dist
        push!(d, _pn(:offset_ra, k, nev) => v[p.centers[k][1]] / D * RAD2MAS,
              _pn(:offset_dec, k, nev) => v[p.centers[k][2]] / D * RAD2MAS)
    end
    return merge(NamedTuple(d), _model_derived(p.model, v, s))
end

canonical(p::LimbProblem, v) = _canonical!(p.model, copy(v), p.setup)
