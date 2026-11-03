# ---------------------------------------------------------------------------
# Results of limb fits: astrometry, outlines, residuals, negative-chord checks,
# predicted contacts
# ---------------------------------------------------------------------------

const LimbFit = FitResult{<:LimbProblem}
const LimbPosterior = Posterior{<:LimbProblem}

function _event_index(p::LimbProblem, event)
    event isa Integer && return event
    i = findfirst(e -> e === event || e.name == event, p.setup.events)
    i === nothing && throw(ArgumentError("unknown event $event"))
    return i
end

"""
    outlines(fit) -> Vector{(center, outline)}
    outlines(problem, values)

Fitted limb of each event: centre `(f, g)` [km] and [`Outline`](@ref).
"""
outlines(p::LimbProblem, v) = [((v[p.centers[k][1]], v[p.centers[k][2]]), _outline(p.model, v, k, p.setup))
                               for k in 1:length(p.setup.events)]
outlines(r::LimbFit) = outlines(r.problem, r.values)

# offsets of the astrometric reference point (the fitted centre, or for a binary the barycentre or
# the photocentre)
_astrometry_names(reference::Symbol) =
    reference in (:center, :primary) ? (:offset_ra, :offset_dec) :
    reference in (:barycenter, :photocenter) ? (Symbol(reference, "_offset_ra"), Symbol(reference, "_offset_dec")) :
    throw(ArgumentError("reference must be :center, :primary, :barycenter or :photocenter"))

"""
    astrometry(fit; event=1, reference=:center) -> NamedTuple
    astrometry(posterior; event=1, reference=:center)

Astrometric geocentric ICRS position of the body at the reference epoch of the event,
from the fitted centre (offset from the ephemeris):
`(epoch, ra, dec, ra_err, dec_err, offset_ra, offset_dec, offset_ra_err, offset_dec_err,
ra_str, dec_str)`. Offsets are Δα cos δ and Δδ [mas]; the errors [mas] combine the fit
uncertainty with the propagated position uncertainty of the star (catalogue errors and
proper motion). For a posterior, the fit uncertainty is the posterior standard deviation.

`reference` is the point whose position is given: the fitted centre (`:center`; for a
[`BinaryLimb`](@ref) the centre of the primary), or for a binary the `:barycenter` (equal
densities) or the `:photocenter` (equal albedos).
"""
function astrometry(r::LimbFit; event=1, reference::Symbol=:center)
    p = r.problem
    k = _event_index(p, event)
    nev = length(p.setup.events)
    na, nd = _astrometry_names(reference)
    haskey(r.derived, _pn(na, k, nev)) || throw(ArgumentError("reference $reference needs a BinaryLimb model"))
    dα, σα = r[_pn(na, k, nev)]
    dδ, σδ = r[_pn(nd, k, nev)]
    return _astrometry(p, k, dα, dδ, σα, σδ)
end

function astrometry(post::LimbPosterior; event=1, reference::Symbol=:center)
    p = post.problem
    k = _event_index(p, event)
    nev = length(p.setup.events)
    na, nd = _astrometry_names(reference)
    d = derived_samples(post)
    haskey(d, _pn(na, k, nev)) || throw(ArgumentError("reference $reference needs a BinaryLimb model"))
    a, b = d[_pn(na, k, nev)], d[_pn(nd, k, nev)]
    return _astrometry(p, k, median(a), median(b), std(a), std(b))
end

function _astrometry(p::LimbProblem, k, dα, dδ, σα, σδ)
    ev = p.setup.events[k]
    ra, dec, _ = ephem_radec(ev.body, ev.tref)
    r2, d2 = offset_radec(ra, dec, dα, dδ)
    sα, sδ = position_error(ev.star, ev.tref)
    return (epoch=ev.tref, ra=r2, dec=d2, ra_err=hypot(σα, sα), dec_err=hypot(σδ, sδ),
            offset_ra=dα, offset_dec=dδ, offset_ra_err=σα, offset_dec_err=σδ,
            ra_str=format_ra(r2; digits=5), dec_str=format_dec(d2; digits=4))
end

"""
    contact_residuals(fit) -> Vector{NamedTuple}

Observed − computed residual of every timed contact: `(event, chord, kind, dt [s], dr [km])`
(time difference between the observed and the predicted contact, and the corresponding
distance along the chord; for a chord that misses the limb, relative to its closest
approach and `miss` [km]).
"""
function contact_residuals(r::LimbFit)
    p, v = r.problem, r.values
    O = _outlines(p, v)
    R = [_rings(p.model, v, k, p.setup) for k in 1:length(p.setup.events)]
    map(p.contacts) do d
        cf, cg = v[p.centers[d.k][1]], v[p.centers[d.k][2]]
        τ = d.τ + (d.dt == 0 ? 0.0 : v[d.dt])
        f, g = d.tr(τ)
        vf, vg = velocity(d.tr, τ)
        s, miss = d.body < 0 ? (d.kind === :mid ? _ring_mid(R[d.k][-d.body], f - cf, g - cg, vf, vg) :
                                _ring_crossing(R[d.k][-d.body], f - cf, g - cg, vf, vg, d.entry)) :
                  _contact_crossing(O[d.k], f - cf, g - cg, vf, vg, d.entry, d.body)
        (event=d.k, chord=d.chord, kind=d.kind, dt=-s, dr=-s * hypot(vf, vg), miss=miss)
    end
end

"""
    clock_offsets(fit) -> Dict{String,Tuple{Float64,Float64}}

Fitted clock offsets `chord => (dt, σ)` [s]; apply them with `shift(chord, dt)`.
"""
function clock_offsets(r::LimbFit)
    out = Dict{String,Tuple{Float64,Float64}}()
    for n in r.names
        s = string(n)
        startswith(s, "dt_") && (out[s[4:end]] = r[n])
    end
    return out
end

"""
    check_negative_chords(event, outline; center=(0, 0), n=400) -> Vector{NamedTuple}
    check_negative_chords(fit; event=1)

Check the negative chords of an event against a limb (any [`Outline`](@ref): a fitted one,
an ellipse from the literature via [`ellipse_outline`](@ref), or the projection of a 3-D
shape via [`project`](@ref)), placed at `center` in the fundamental plane. For each
negative chord: `(chord, conflict, distance, epoch, immersion, emersion, duration)` —
whether the chord enters the limb, its smallest signed distance from it [km] (negative:
inside) at `epoch`, and for a conflict the predicted, unobserved contacts and duration [s].
Sorted by distance (worst first).
"""
function check_negative_chords(ev::Event, O::Outline; center=(0.0, 0.0), n::Integer=400)
    out = NamedTuple[]
    for c in ev.chords
        isnegative(c) || continue
        tr = track(ev, c)
        a, b = c.span[1] - ev.tref, c.span[2] - ev.tref
        h(τ) = ((f, g) = tr(τ); signed_distance(O, f - center[1], g - center[2]))
        τs = range(a, b; length=max(n, 2))
        hs = h.(τs)
        i = argmin(hs)
        τm, hm = _golden_min(h, τs[max(i - 1, 1)], τs[min(i + 1, end)]; tol=1e-4)
        hm > hs[i] && ((τm, hm) = (τs[i], hs[i]))
        function root(lo, hi)
            hlo = h(lo)
            for _ in 1:50
                m = (lo + hi) / 2
                hmid = h(m)
                (hmid >= 0) == (hlo >= 0) ? ((lo, hlo) = (m, hmid)) : (hi = m)
            end
            return (lo + hi) / 2
        end
        conflict = hm < 0
        ti = te = nothing
        if conflict
            ti = h(a) < 0 ? nothing : ev.tref + root(a, τm)
            te = h(b) < 0 ? nothing : ev.tref + root(τm, b)
        end
        dur = conflict ? (something(te, ev.tref + b) - something(ti, ev.tref + a)) : 0.0
        push!(out, (chord=c.name, conflict=conflict, distance=hm, epoch=ev.tref + τm, immersion=ti,
                    emersion=te, duration=dur))
    end
    return sort!(out; by=x -> x.distance)
end

function check_negative_chords(r::LimbFit; event=1, n::Integer=400)
    k = _event_index(r.problem, event)
    c, O = outlines(r)[k]
    return check_negative_chords(r.problem.setup.events[k], O; center=c, n=n)
end

"""
    predicted_contacts(event, observer, outline; center=(0, 0), window=600, step=0.5)
        -> Vector{Tuple{Epoch,Epoch}}

Immersion/emersion instants predicted for an observer by a limb (e.g. a fit, to plan or
check observations), searched within `event.tref ± window` seconds.
"""
function predicted_contacts(ev::Event, obs::AbstractObserver, O::Outline; center=(0.0, 0.0),
                            window::Real=600.0, step::Real=0.5)
    tr = Track(ev.star, ev.body, obs, ev.tref, -window, window)
    h(τ) = ((f, g) = tr(τ); signed_distance(O, f - center[1], g - center[2]))
    τs = collect(-window:step:window)
    hs = h.(τs)
    function root(lo, hi)
        hlo = h(lo)
        for _ in 1:50
            m = (lo + hi) / 2
            hm = h(m)
            (hm >= 0) == (hlo >= 0) ? ((lo, hlo) = (m, hm)) : (hi = m)
        end
        return (lo + hi) / 2
    end
    out = Tuple{Epoch,Epoch}[]
    tin = nothing
    for i in 1:length(τs)-1
        if hs[i] >= 0 && hs[i+1] < 0
            tin = root(τs[i], τs[i+1])
        elseif hs[i] < 0 && hs[i+1] >= 0 && tin !== nothing
            push!(out, (ev.tref + tin, ev.tref + root(τs[i], τs[i+1])))
            tin = nothing
        end
    end
    return out
end

function predicted_contacts(r::LimbFit, obs::AbstractObserver; event=1, kw...)
    k = _event_index(r.problem, event)
    c, O = outlines(r)[k]
    return predicted_contacts(r.problem.setup.events[k], obs, O; center=c, kw...)
end

"""
    model_lightcurve(fit, chord_name) -> (time, flux)

Model flux of a light curve in a joint photometric fit (at the light-curve times, with the
profiled flux levels).
"""
function model_lightcurve(r::LimbFit, name::AbstractString)
    p, v = r.problem, r.values
    i = findfirst(d -> d.chord == name, p.photometry)
    i === nothing && throw(ArgumentError("chord $name is not part of the photometric fit"))
    d = p.photometry[i]
    I, w = _phot_model(p, d, _outlines(p, v)[d.k], v)
    ft, fb = _flux_levels(I, d.lc.flux, w, 1.0, 0.0, true, p.atmosphere === nothing && count(<(0.1), I) >= 3)
    return d.lc.time, fb .+ (ft - fb) .* I
end

"""
    fit_limb(events; starts=8, rng, σ_limb=0, kwargs...) -> FitResult

Shortcut for `fit(LimbProblem(events; kwargs...); starts, rng)`.

`σ_limb = :auto` determines the limb uncertainty allowance [km] (added in quadrature to every
contact, see [`LimbProblem`](@ref)) for which χ²/dof = 1, as SORA's `ellipse_error`: the
excess scatter of the contacts (limb topography, unmodelled timing errors) is then reflected
in the uncertainties. The value used is `fit.problem.σ_limb`.
"""
function fit_limb(evs; starts::Integer=8, rng=Random.default_rng(), scale_errors::Bool=false, σ_limb=0.0, kw...)
    σ_limb === :auto || return fit(LimbProblem(evs; σ_limb=σ_limb, kw...); starts=starts, rng=rng, scale_errors=scale_errors)
    r = fit(LimbProblem(evs; σ_limb=0.0, kw...); starts=starts, rng=rng)
    reduced_chi2(r) <= 1 && return r
    init = Dict(n => v for (n, v) in zip(r.names, r.values))
    f(σ) = fit(LimbProblem(evs; σ_limb=σ, kw..., init=init); starts=1, rng=rng)
    lo, hi = 0.0, 1.0
    while reduced_chi2(f(hi)) > 1
        lo, hi = hi, 2hi
        hi > 1e5 && error("σ_limb = :auto: no solution below 10⁵ km")
    end
    for _ in 1:30
        m = (lo + hi) / 2
        reduced_chi2(f(m)) > 1 ? (lo = m) : (hi = m)
        hi - lo < 1e-3 * hi && break
    end
    return fit(LimbProblem(evs; σ_limb=hi, kw..., init=init); starts=starts, rng=rng, scale_errors=scale_errors)
end

"""
    atmosphere_test(events; atmosphere, model=CircleLimb(), pressures=nothing, confidence=0.95,
                    init=Dict(), airless_init=Dict(), starts=3, rng, kwargs...) -> NamedTuple

Test light curves for an atmosphere: a profile likelihood in the surface pressure. The chords
must carry light curves (joint photometric fit, `photometry = true`); the body needs `GM` or
`density`. The limb (and centre) is refitted for every pressure of the grid `pressures` [nbar]
(default 0 and 0.3–3000 nbar, logarithmic), which gives

* `pressure` – the best-fitting surface pressure and `chi2_min`;
* `chi2_airless` – χ² for p = 0 (an airless body with the same limb model; refitted with
  `starts` multistarts from the free fit and from `init` updated by `airless_init`, e.g. the
  values of an airless fit, since its limb can lie far from the atmospheric one) and the
  detection `significance` √(χ²(0) − χ²_min) [σ] (likelihood ratio, one parameter);
* `upper_limit` – the one-sided upper limit at `confidence` (Δχ² = z², z the normal quantile:
  2.71 for 95 %, 5.41 for 99 %, 9 for 3σ) and `interval = (lower, upper)`, where `lower > 0`
  only if an atmosphere is detected at that confidence;
* `pressures`, `chi2` – the profile, and `fit` – the free fit.

Other keywords go to [`LimbProblem`](@ref) (e.g. `fit_dt`, `σ_limb`, `priors`).
"""
function atmosphere_test(evs; atmosphere::Atmosphere, model::LimbModel=CircleLimb(), pressures=nothing,
                         confidence::Real=0.95, init=Dict(), airless_init=Dict(), starts::Integer=3,
                         rng=Random.default_rng(), kw...)
    prob = LimbProblem(evs; model=model, photometry=true, atmosphere=atmosphere, init=init, kw...)
    ps = parameters(prob)
    ip = findfirst(q -> q.name === :pressure, ps)
    saved = [(q.value, q.free) for q in ps]
    saved_values = [Float64(get(airless_init, q.name, q.value)) for q in ps]   # start for p = 0
    restore!() = for (j, q) in enumerate(ps); q.value, q.free = saved[j]; end
    # χ² with the pressure fixed at p, starting from the parameter values v
    function prof(p, v; n=1)
        for (j, q) in enumerate(ps)
            q.value = v[j]
        end
        ps[ip].value, ps[ip].free = p, false
        r = try
            fit(prob; starts=n, rng=rng)
        catch
            nothing
        end
        ps[ip].free = true
        return r === nothing ? (Inf, v) : (r.chi2, r.values)
    end
    best = fit(prob; starts=starts, rng=rng)
    grid = pressures === nothing ? [0.0; exp10.(range(-0.5, 3.5; length=17))] : sort(collect(Float64, pressures))
    χ = similar(grid)
    sol = Vector{Vector{Float64}}(undef, length(grid))
    try
        v = copy(best.values)
        for (k, p) in enumerate(grid)
            # the airless solution may lie far from the atmospheric one: search it with multistarts
            χ[k], sol[k] = prof(p, v; n=p == 0 ? starts : 1)
            if p == 0
                c, s0 = prof(p, saved_values; n=starts)     # also from the initial values
                c < χ[k] && ((χ[k], sol[k]) = (c, s0))
            end
            isfinite(χ[k]) && (v = sol[k])                 # continue from the last solution
        end
        # and downwards from the free fit: keep the lower envelope of both sweeps
        v = copy(best.values)
        kb = searchsortedfirst(grid, best[:pressure][1])
        for k in reverse(1:min(kb - 1, length(grid)))
            c, s1 = prof(grid[k], v)
            c < χ[k] && ((χ[k], sol[k]) = (c, s1))
            isfinite(χ[k]) && (v = sol[k])
        end
        v = copy(best.values)
        for k in kb:length(grid)
            c, s1 = prof(grid[k], v)
            c < χ[k] && ((χ[k], sol[k]) = (c, s1))
            isfinite(χ[k]) && (v = sol[k])
        end
        # the free fit may have stopped in a local minimum: restart from the best grid point
        k0 = argmin(χ)
        if χ[k0] < best.chi2 - 1e-6
            for (j, q) in enumerate(ps)
                q.value = sol[k0][j]
            end
            ps[ip].free = true
            r = fit(prob; starts=1, rng=rng)
            r.chi2 < best.chi2 && (best = r)
        end
    finally
        restore!()
    end
    χmin = min(best.chi2, minimum(filter(isfinite, χ); init=Inf))
    z = sqrt(2) * SpecialFunctions.erfinv(2confidence - 1)
    thr = χmin + z^2
    χ0 = grid[1] == 0 ? χ[1] : NaN
    # limits of the region χ²(p) ≤ χ²_min + z² around the best pressure p̂: on each side, the
    # nearest grid point outside the region, then bisection (with refits) between it and p̂ –
    # the limits bracket p̂ by construction, even where the profile is not smooth
    p̂ = best[:pressure][1]
    vbest = copy(best.values)
    outside(k) = !(isfinite(χ[k]) && χ[k] <= thr)
    function crossing(pin, pout)                  # pin inside the region, pout outside
        a, b = pin, pout
        v = vbest
        uselog = a > 0 && b > 0
        for _ in 1:16
            m = uselog ? sqrt(a * b) : (a + b) / 2
            c, vm = prof(m, v)
            c <= thr ? (a = m; v = vm) : (b = m)
        end
        return uselog ? sqrt(a * b) : (a + b) / 2
    end
    up, lo = NaN, 0.0
    try
        if best.chi2 <= thr
            ku = findfirst(k -> grid[k] > p̂ && outside(k), eachindex(grid))
            up = ku === nothing ? Inf : crossing(p̂, grid[ku])
            kl = findlast(k -> grid[k] < p̂ && outside(k), eachindex(grid))
            lo = kl === nothing ? 0.0 : crossing(p̂, grid[kl])
        end
    finally
        restore!()
    end
    sig = isfinite(χ0) ? sqrt(max(χ0 - χmin, 0.0)) : NaN
    return (pressure=best[:pressure][1], chi2_min=χmin, chi2_airless=χ0, significance=sig, upper_limit=up,
            interval=(lo, up), confidence=Float64(confidence), pressures=grid, chi2=χ, fit=best)
end
