# ---------------------------------------------------------------------------
# Ring systems: circular rings in the equatorial plane of the ring pole, centred on the body.
# On the sky a ring is an elliptical annulus (opening angle B, position angle P of the
# projected pole), partly transparent: apparent opacity p′ = 1 − exp(−τ_N / |sin B|).
# ---------------------------------------------------------------------------

"""
    ProjectedRing(radius, width, opacity, pa, opening)
    project_ring(radius, width, tau, pole_ra, pole_dec, u) -> ProjectedRing

A ring on the sky plane: mid radius `radius` and radial width `width` [km] of a circular ring
seen at the opening angle `opening` B [deg] (the angle between the line of sight and the ring
plane; B > 0: the ring pole tilted towards the observer) with the projected ring
pole at the position angle `pa` [deg, north towards east]; `opacity` is the apparent opacity
p′ (the fraction of the stellar flux blocked). Its edges are the ellipses with semi-axes
`r ± W/2` and `(r ± W/2)|sin B|`, the minor axis along `pa`.

[`project_ring`](@ref) gives it from the ring pole and the normal optical depth.

[`intervals`](@ref)`(ring, x, y, vx, vy)` gives the crossings of a track (two per crossing of
the annulus); `outline_points(ring)` both edges (for plotting).
"""
struct ProjectedRing{T<:Real}
    radius::T
    width::T
    opacity::T
    pa::T
    opening::T
end
ProjectedRing(r, w, p, pa, B) = ProjectedRing(promote(r, w, p, pa, B)...)

Base.show(io::IO, r::ProjectedRing) = @printf(io, "ProjectedRing(r = %.2f km, W = %.2f km, p′ = %.3f, P = %.2f°, B = %.2f°)",
                                              r.radius, r.width, r.opacity, r.pa, r.opening)

# the two edges as ellipses (outer, inner); the inner one is `nothing` for W ≥ 2r
function _ring_edges(r::ProjectedRing)
    s = abs(sind(r.opening))
    ro, ri = r.radius + r.width / 2, r.radius - r.width / 2
    return EllipseOutline(ro, ro * s, r.pa), ri > 0 ? EllipseOutline(ri, ri * s, r.pa) : nothing
end

function intervals(r::ProjectedRing, x, y, vx, vy)
    T = promote_type(typeof(x), typeof(y), typeof(vx), typeof(vy), typeof(r.radius), Float64)
    out, inn = _ring_edges(r)
    a = _interval(out, x, y, vx, vy)
    a === nothing && return Tuple{T,T}[]
    b = inn === nothing ? nothing : _interval(inn, x, y, vx, vy)
    b === nothing && return Tuple{T,T}[(a[1], a[2])]
    return Tuple{T,T}[(a[1], b[1]), (b[2], a[2])]
end

function outline_points(r::ProjectedRing; center=(0.0, 0.0), n::Integer=361)
    out, inn = _ring_edges(r)
    x, y = outline_points(out; center=center, n=n)
    inn === nothing && return x, y
    xi, yi = outline_points(inn; center=center, n=n)
    return [x; NaN; xi], [y; NaN; yi]
end

"""
    ring_geometry(pole_ra, pole_dec, u) -> (opening, pa)

Opening angle B [deg] of a ring plane with the pole (`pole_ra`, `pole_dec`) seen along `u`
(observer → body), and the position angle P [deg] of the projected pole.
"""
function ring_geometry(α, δ, u::AbstractVector)
    p = radec2vec(α, δ)
    û = u ./ norm(u)
    e, n = sky_basis(û)
    B = asind(clamp(-dot(p, û), -1, 1))
    return B, mod(atand(dot(p, e), dot(p, n)), 360)
end

"""
    project_ring(radius, width, tau, pole_ra, pole_dec, u) -> ProjectedRing

The [`ProjectedRing`](@ref) of a circular ring of mid radius `radius` and radial width `width`
[km] with the normal optical depth `tau`, in the plane of the pole (`pole_ra`, `pole_dec`
[deg]), seen along `u` (observer → body): the opening angle and position angle from
[`ring_geometry`](@ref), the apparent opacity `1 − exp(−tau/|sin B|)`.
"""
function project_ring(radius, width, tau, α, δ, u::AbstractVector)
    B, P = ring_geometry(α, δ, u)
    s = max(abs(sind(B)), 1e-6)
    return ProjectedRing(radius, width, 1 - exp(-tau / s), P, B)
end

# the ring pole (ra, dec) [deg] for an opening angle B and a projected-pole position angle P
# seen along u
function _ring_pole(B, P, u::AbstractVector)
    û = u ./ norm(u)
    e, n = sky_basis(û)
    p = -sind(B) .* û .+ cosd(B) .* (sind(P) .* e .+ cosd(P) .* n)
    return mod(atand(p[2], p[1]), 360), asind(clamp(p[3], -1, 1))
end

"""
    predicted_contacts(event, observer, ring::ProjectedRing; center=(0, 0), window=600)
        -> Vector{Tuple{Epoch,Epoch}}

The crossings of a ring (start and end of each passage through the annulus) predicted for an
observer, e.g. to simulate observations.
"""
function predicted_contacts(ev::Event, obs::AbstractObserver, r::ProjectedRing; center=(0.0, 0.0), window::Real=600.0)
    tr = Track(ev.star, ev.body, obs, ev.tref, -window, window)
    # the track is straight over the ring passage: intervals at the closest approach to the centre
    τs = range(-window, window; length=2001)
    d = [hypot((tr(τ) .- center)...) for τ in τs]
    τ0 = τs[argmin(d)]
    f, g = tr(τ0)
    vf, vg = velocity(tr, τ0)
    # each crossing refined on the actual (curved) track: the nearest crossing of the same kind
    # from the track position and velocity at the current estimate
    function refine(τ, entry)
        for _ in 1:4
            x, y = tr(τ)
            ux, uy = velocity(tr, τ)
            s, _ = _ring_crossing(r, x - center[1], y - center[2], ux, uy, entry)
            τ += s
            abs(s) < 1e-6 && break
        end
        return τ
    end
    return [(ev.tref + refine(τ0 + a, true), ev.tref + refine(τ0 + b, false))
            for (a, b) in intervals(r, f - center[1], g - center[2], vf, vg)]
end

# ---------------------------------------------------------------------------
# the limb model: a body with rings
# ---------------------------------------------------------------------------

"""
    RingedLimb(body=EllipseLimb(); rings=1, epoch=nothing)

A body with `rings` circular rings in the plane of the ring pole, centred on the body. The
body is any [`LimbModel`](@ref) (its parameters unchanged). The rings add the ring pole
`ring_pole_ra`, `ring_pole_dec` [deg] and for each ring `k` its mid radius `ringk_radius`,
radial width `ringk_width` [km] and normal optical depth `ringk_tau`.

**Data**: the ring dips of the chords are contacts labelled with their ring,
`Contact(:immersion, t, σ; body=:ring1)` and its emersion (the first ring dip of a chord is
its passage through the annulus, the edges are the contacts): each is compared with the
nearest crossing of the same kind of that ring's annulus along the chord (residuals along the
chord). A dip known only by its mid-time is a `Contact(:mid, t, σ; body=:ring1)`, compared
with the mid-point of the nearest passage through the annulus (it constrains the ring's
geometry, not its width). A chord that saw only ring dips is a negative chord for the body. Unlabelled contacts belong to the body. Negative chords and the visible segments of
the positive chords constrain the body only (a ring may escape detection). In a joint fit of
the light curves (`photometry = true`) the rings are partly transparent strips (apparent
opacity p′), with diffraction.

A pole and its antipode describe the same ring plane; the fitted pole is the one in the
northern hemisphere (`ring_pole_dec` ≥ 0). A single event leaves a twofold ambiguity: the
mirror plane (opening angle −B, the same projected ellipse); several events at different
geometries, or other information, decide between them.

**Derived**: per event the opening angle `ring_opening` B and the position angle `ring_pa` P
of the projected pole; per ring the `ringk_equivalent_depth` W τ_N [km] and per event its
apparent opacity `ringk_opacity` and equivalent width `ringk_equivalent_width` W p′ [km].

The starting values come from the ring contacts (an ellipse through the mid-points of the
ring dips gives the radius, the opening angle and the position angle) unless given with `init`.
"""
struct RingedLimb{M<:LimbModel} <: LimbModel
    body::M
    nrings::Int
    epoch::Union{Nothing,Epoch}
end
function RingedLimb(body::LimbModel=EllipseLimb(); rings::Integer=1, epoch=nothing)
    1 <= rings <= 9 || throw(ArgumentError("1 to 9 rings"))
    return RingedLimb(body, Int(rings), epoch === nothing ? nothing : _epoch(epoch))
end

_composite(::RingedLimb) = true
function _body_index(m::RingedLimb, b::Symbol)
    s = string(b)
    if startswith(s, "ring")
        k = s == "ring" ? 1 : parse(Int, s[5:end])
        1 <= k <= m.nrings || throw(ArgumentError("contact label $b: the model has $(m.nrings) ring(s)"))
        return -k
    end
    return _body_index(m.body, b)
end
_prepare!(m::RingedLimb, s::_Setup) = _prepare!(m.body, s)
_finalize!(m::RingedLimb, s::_Setup) = _finalize!(m.body, s)
_outline(m::RingedLimb, v, k::Integer, s::_Setup) = _outline(m.body, v, k, s)
_admissible(m::RingedLimb, v, s::_Setup) = _admissible(m.body, v, s)
_shape(m::RingedLimb, v, s::_Setup) = _shape(m.body, v, s)

_rn(k, name) = Symbol("ring", k, "_", name)

# mid-points of the ring dips (by ring) of event k in the fundamental plane
function _ring_points(ev::Event, ring::Integer)
    out = Tuple{Float64,Float64,Float64}[]                  # f, g, length of the dip [km]
    lab = (Symbol("ring", ring), ring == 1 ? :ring : Symbol("ring", ring))
    for c in ev.chords
        cs = [x for x in c.contacts if x.body in lab]
        tr = track(ev, c)
        for x in cs
            x.kind === :mid || continue
            a = tr(x.time - ev.tref)
            push!(out, (a[1], a[2], 0.0))
        end
        for i in 1:length(cs)-1
            (cs[i].kind === :immersion && cs[i+1].kind === :emersion) || continue
            a, b = tr(cs[i].time - ev.tref), tr(cs[i+1].time - ev.tref)
            push!(out, ((a[1] + b[1]) / 2, (a[2] + b[2]) / 2, hypot(b[1] - a[1], b[2] - a[2])))
        end
    end
    return out
end

function _model_params(m::RingedLimb, s::_Setup)
    ps = _model_params(m.body, s)
    est = s.estimates[1]
    # ring 1: an ellipse through the mid-points of its dips (relative to the body centre)
    pts = _ring_points(s.events[1], 1)
    B0, P0, r1 = 30.0, 0.0, 3 * est.radius
    if length(pts) >= 5
        e = algebraic_ellipse([p[1] for p in pts], [p[2] for p in pts])
        if e !== nothing
            r1 = e.equatorial_radius
            B0 = asind(clamp(1 - e.oblateness, 0.05, 1.0))
            P0 = e.position_angle
        end
    elseif !isempty(pts)
        r1 = mean(hypot(p[1] - est.center_f, p[2] - est.center_g) for p in pts)
    end
    α0, δ0 = _ring_pole(B0, P0, s.contexts[1].u)
    push!(ps, Param(:ring_pole_ra, α0; scale=5.0, period=360.0),
          Param(:ring_pole_dec, δ0; lower=-90.0, upper=90.0, scale=5.0))
    for k in 1:m.nrings
        q = _ring_points(s.events[1], k)
        # ring 1 from its ellipse; the others scaled by the mean distance of their dips
        dist(ps_) = mean(hypot(p[1] - est.center_f, p[2] - est.center_g) for p in ps_)
        rk = k == 1 ? r1 : isempty(q) ? r1 * (1 + 0.05(k - 1)) : (isempty(pts) ? dist(q) : r1 * dist(q) / dist(pts))
        w0 = isempty(q) ? 5.0 : max(median([p[3] for p in q]) * abs(sind(B0)), 0.5)
        push!(ps, Param(_rn(k, :radius), rk; lower=1e-3, scale=max(rk / 100, 1.0)),
              Param(_rn(k, :width), w0; lower=1e-3, scale=max(w0 / 5, 0.1)),
              Param(_rn(k, :tau), 0.3; lower=0.0, upper=20.0, scale=0.05))
    end
    return ps
end

# a pole and its antipode describe the same ring plane: canonical pole in the northern hemisphere
function _canonical!(m::RingedLimb, v, s::_Setup)
    _canonical!(m.body, v, s)
    ia, id = s.index[:ring_pole_ra], s.index[:ring_pole_dec]
    if v[id] < 0
        v[ia] += 180
        v[id] = -v[id]
    end
    v[ia] = mod(v[ia], 360)
    return v
end

# the rings of event k: Vector{ProjectedRing}
function _rings(m::RingedLimb, v, k::Integer, s::_Setup)
    α, δ = _v(v, s, :ring_pole_ra), _v(v, s, :ring_pole_dec)
    return [project_ring(_v(v, s, _rn(j, :radius)), _v(v, s, _rn(j, :width)), _v(v, s, _rn(j, :tau)), α, δ,
                         s.contexts[k].u) for j in 1:m.nrings]
end
_rings(::LimbModel, v, k::Integer, s::_Setup) = ProjectedRing[]

function _model_derived(m::RingedLimb, v, s::_Setup)
    nev = _nev(s)
    d = Pair{Symbol,Float64}[pairs(_model_derived(m.body, v, s))...]
    for j in 1:m.nrings
        push!(d, _rn(j, :equivalent_depth) => _v(v, s, _rn(j, :width)) * _v(v, s, _rn(j, :tau)))
    end
    for k in 1:nev
        R = _rings(m, v, k, s)
        push!(d, _pn(:ring_opening, k, nev) => R[1].opening, _pn(:ring_pa, k, nev) => R[1].pa)
        for j in 1:m.nrings
            push!(d, _pn(_rn(j, :opacity), k, nev) => R[j].opacity,
                  _pn(_rn(j, :equivalent_width), k, nev) => R[j].width * R[j].opacity)
        end
    end
    return NamedTuple(d)
end

# a ring dip reported by its mid-time: the nearest mid-point of a passage through the annulus
function _ring_mid(r::ProjectedRing, x, y, vx, vy)
    iv = intervals(r, x, y, vx, vy)
    if isempty(iv)
        out, _ = _ring_edges(r)
        return crossing(out, x, y, vx, vy, true)
    end
    ss = [(a + b) / 2 for (a, b) in iv]
    s = ss[argmin(abs.(ss))]
    return s, zero(s)
end

# crossing of a ring contact: the nearest crossing of the same kind of the annulus of ring j
function _ring_crossing(r::ProjectedRing, x, y, vx, vy, entry::Bool)
    iv = intervals(r, x, y, vx, vy)
    if isempty(iv)
        out, _ = _ring_edges(r)
        return crossing(out, x, y, vx, vy, entry)
    end
    ss = [entry ? a : b for (a, b) in iv]
    s = ss[argmin(abs.(ss))]
    return s, zero(s)
end
