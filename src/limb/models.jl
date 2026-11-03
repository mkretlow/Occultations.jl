# ---------------------------------------------------------------------------
# Limb models: parameterised apparent limbs. A model declares its parameters and, for
# parameter values, returns the outline of the body in each event.
# ---------------------------------------------------------------------------

"""
    LimbModel

Parameterised model of the apparent limb, used by [`LimbProblem`](@ref):

* [`EllipseLimb`](@ref) – 2-D ellipse (default), optionally with event-specific parameters
* [`CircleLimb`](@ref) – circle
* [`EllipsoidLimb`](@ref) – triaxial ellipsoid with pole and rotation phase: the apparent
  ellipse of each event follows from the body orientation at that epoch (multi-event fits,
  rotational-amplitude constraints)
* [`MeshLimb`](@ref) – silhouette of a triangle mesh (e.g. an OBJ shape model), scaled
  and optionally rotated, with fixed or fitted orientation
* [`BinaryLimb`](@ref) – two bodies (a binary system, a body and its satellite), each with
  its own model
* [`ContactBinaryLimb`](@ref), [`ContactBinaryEllipsoids`](@ref) – a contact binary: two
  overlapping ellipses (2-D) or two rotating ellipsoids (3-D)

A new model implements `_model_params(model, setup)` and `_outline(model, values, k, setup)`.
"""
abstract type LimbModel end

# per-event parameter name
_pn(base::Symbol, k::Integer, nev::Integer) = nev == 1 ? base : Symbol(base, "_", k)

# information available to the models when the problem is built and evaluated
struct _Setup
    events::Vector{Event}
    contexts::Vector{EventContext}
    estimates::Vector{NamedTuple}       # per event: center_f, center_g, radius, ellipse (or nothing)
    index::Dict{Symbol,Int}             # parameter name -> position (filled after construction)
    epoch::Epoch                        # reference epoch of rotation phases
    aux::Dict{Symbol,Any}               # model-specific data (e.g. the setups of the components of a binary)
end
_Setup(events, contexts, estimates, index, epoch) = _Setup(events, contexts, estimates, index, epoch, Dict{Symbol,Any}())
_nev(s::_Setup) = length(s.events)
@inline _v(v, s::_Setup, name::Symbol) = v[s.index[name]]
@inline _v(v, s::_Setup, base::Symbol, k::Integer, per::Bool) = v[s.index[per ? _pn(base, k, _nev(s)) : base]]

_model_derived(::LimbModel, v, s::_Setup) = NamedTuple()
# hooks of composite models: before the parameters are declared (e.g. to split the starting
# estimates between the components) and once the parameter index is known
_prepare!(::LimbModel, s::_Setup) = nothing
_finalize!(::LimbModel, s::_Setup) = nothing
# whether the model's outlines are composite (several bodies or lobes)
_composite(::LimbModel) = false
# index of the component of a contact's `body` label (0: matched automatically)
_body_index(::LimbModel, ::Symbol) = 0
_admissible(::LimbModel, v, s::_Setup) = true
_canonical!(::LimbModel, v, s::_Setup) = v

# ---------------------------------------------------------------------------

"""
    EllipseLimb(; per_event=())

Apparent ellipse with parameters `equatorial_radius` (semi-major axis, km), `oblateness`
(`1 - b/a`) and `position_angle` (of the semi-minor axis = apparent pole, from north
towards east, deg). In a multi-event fit the parameters listed in `per_event` get one
value per event (`position_angle_1`, `position_angle_2`, …), the others are shared, e.g.
`EllipseLimb(per_event=(:position_angle,))` for a body seen at different rotation phases
with the same apparent axes.
"""
struct EllipseLimb <: LimbModel
    per_event::Vector{Symbol}
end
function EllipseLimb(; per_event=())
    pe = Symbol.(collect(per_event))
    all(in((:equatorial_radius, :oblateness, :position_angle)), pe) ||
        throw(ArgumentError("per_event may contain :equatorial_radius, :oblateness, :position_angle"))
    return EllipseLimb(pe)
end

function _model_params(m::EllipseLimb, s::_Setup)
    nev = _nev(s)
    est = s.estimates
    ps = Param[]
    # ε < 0 (the same ellipse with the axes swapped) removes the hard wall at ε = 0;
    # allowed when a and ε share the same per-event structure (canonicalised after the fit)
    swap = (:equatorial_radius in m.per_event) == (:oblateness in m.per_event)
    for (base, init, kw) in ((:equatorial_radius, k -> est[k].ellipse === nothing ? est[k].radius : est[k].ellipse.equatorial_radius,
                              (lower=1e-3,)),
                             (:oblateness, k -> est[k].ellipse === nothing ? 0.0 : est[k].ellipse.oblateness,
                              (lower=swap ? -0.95 : 0.0, upper=0.95, scale=0.05)),
                             (:position_angle, k -> est[k].ellipse === nothing ? 0.0 : est[k].ellipse.position_angle,
                              (scale=10.0, period=180.0)))
        ks = base in m.per_event ? (1:nev) : (0:0)
        for k in ks
            name = k == 0 ? base : _pn(base, k, nev)
            v0 = k == 0 ? mean(init(j) for j in 1:nev) : init(k)
            sc = base === :equatorial_radius ? max(v0 / 10, 1.0) : get(kw, :scale, 1.0)
            push!(ps, Param(name, v0; lower=get(kw, :lower, -Inf), upper=get(kw, :upper, Inf), scale=sc,
                            period=get(kw, :period, 0.0)))
        end
    end
    return ps
end

@inline function _ellipse_from(a, ε, pa)
    ε >= 0 && return EllipseOutline(a, a * (1 - ε), pa)
    return EllipseOutline(a * (1 - ε), a, pa + 90)          # axes swapped
end

function _outline(m::EllipseLimb, v, k::Integer, s::_Setup)
    a = _v(v, s, :equatorial_radius, k, :equatorial_radius in m.per_event)
    ε = _v(v, s, :oblateness, k, :oblateness in m.per_event)
    pa = _v(v, s, :position_angle, k, :position_angle in m.per_event)
    return _ellipse_from(a, ε, pa)
end

function _canonical!(m::EllipseLimb, v, s::_Setup)
    nev = _nev(s)
    pera, pere, perp = (x in m.per_event for x in (:equatorial_radius, :oblateness, :position_angle))
    pera == pere || return v
    for k in (pere ? (1:nev) : (0:0))
        nε = k == 0 ? :oblateness : _pn(:oblateness, k, nev)
        na = k == 0 ? :equatorial_radius : _pn(:equatorial_radius, k, nev)
        ε = v[s.index[nε]]
        ε < 0 || continue
        a = v[s.index[na]]
        v[s.index[na]] = a * (1 - ε)
        v[s.index[nε]] = 1 - 1 / (1 - ε)
        pks = perp ? (k == 0 ? (1:nev) : (k:k)) : (0:0)
        for j in pks
            np = j == 0 ? :position_angle : _pn(:position_angle, j, nev)
            v[s.index[np]] = mod(v[s.index[np]] + 90, 180)
        end
    end
    return v
end

function _model_derived(m::EllipseLimb, v, s::_Setup)
    nev = _nev(s)
    per = any(in(m.per_event), (:equatorial_radius, :oblateness))
    d = Pair{Symbol,Float64}[]
    for k in (per ? (1:nev) : (1:1))
        o = _outline(m, v, k, s)
        sfx(n) = per ? _pn(n, k, nev) : n
        push!(d, sfx(:equivalent_radius) => sqrt(o.a * o.b), sfx(:semi_minor_axis) => o.b)
    end
    return NamedTuple(d)
end

# ---------------------------------------------------------------------------

"""
    CircleLimb()

Circular limb with parameter `radius` [km] (shared by all events).
"""
struct CircleLimb <: LimbModel end

function _model_params(::CircleLimb, s::_Setup)
    r0 = mean(e.radius for e in s.estimates)
    return [Param(:radius, r0; lower=1e-3, scale=max(r0 / 10, 1.0))]
end
_outline(::CircleLimb, v, k::Integer, s::_Setup) = (r = _v(v, s, :radius); EllipseOutline(r, r, zero(r)))

# ---------------------------------------------------------------------------
# rotating bodies
# ---------------------------------------------------------------------------

# ICRS → body matrix for pole (α, δ) and prime meridian W [deg]
_body_matrix(α, δ, W) = _R3(W) * _R1(90 - δ) * _R3(90 + α)

# prime meridian of event k: W = W0 + 360°·24/P · (days since the reference epoch), at the
# light-emission time; P [h] (negative: retrograde)
function _meridian(v, s::_Setup, k::Integer)
    W0, P = _v(v, s, :W0), _v(v, s, :period)
    t1, t2 = tdb_jd(s.epoch)
    days = s.contexts[k].jd_emit - (t1 + t2)
    return W0 + 360 * 24 / P * days
end

function _rotation_params(b::Body, s::_Setup; free::Bool=true)
    fr = b.frame
    if fr === nothing
        α0, δ0, W0, P = 0.0, 90.0, 0.0, 10.0
    else
        α0, δ0, W = orientation_at(fr, s.epoch)
        W0, P = W, fr.rotation_rate == 0 ? Inf : 360 * 24 / fr.rotation_rate
    end
    return [Param(:pole_ra, α0; free=free, scale=10.0, period=360.0),
            Param(:pole_dec, δ0; free=free, lower=-90.0, upper=90.0, scale=10.0),
            Param(:W0, W0; free=free, scale=20.0, period=360.0),
            Param(:period, P; free=false, scale=max(abs(P) * 1e-4, 1e-6))]
end

"""
    EllipsoidLimb(; epoch=nothing)

Triaxial ellipsoid with semi-axes `a`, `b`, `c` [km] (body x, y, z; z = rotation axis),
pole `pole_ra`, `pole_dec` [deg], prime-meridian angle `W0` [deg] at the reference `epoch`
(default: the first event) and sidereal `period` [h] (fixed by default; from the body's
frame). For every event the apparent ellipse is the exact projection of the ellipsoid
oriented for that epoch (light-time corrected), so several events at different rotation
phases and aspects constrain the 3-D shape. Rotational light-curve amplitudes can be
added as constraints ([`RotationalAmplitude`](@ref)). Starting values come from the
body's `shape` (an `Ellipsoid`) and `frame`. With `principal=true` only rotation about the
shortest axis (c ≤ b ≤ a) is admitted. For an unknown pole use a flat prior on `pole_ra`
and [`isotropic_pole`](@ref) on `pole_dec` when sampling the posterior.
"""
struct EllipsoidLimb <: LimbModel
    epoch::Union{Nothing,Epoch}
    principal::Bool
end
EllipsoidLimb(; epoch=nothing, principal::Bool=false) = EllipsoidLimb(epoch === nothing ? nothing : _epoch(epoch), principal)

# rotation about the shortest axis: c ≤ b ≤ a (b ≤ a by the canonical form)
_admissible(m::EllipsoidLimb, v, s::_Setup) = !m.principal || _v(v, s, :c) <= min(_v(v, s, :a), _v(v, s, :b))

function _model_params(m::EllipsoidLimb, s::_Setup)
    b = s.events[1].body
    r0 = mean(e.radius for e in s.estimates)
    a0, b0, c0 = b.shape isa Ellipsoid ? (b.shape.a, b.shape.b, b.shape.c) : (1.1r0, r0, 0.9r0)
    sc = max(r0 / 10, 1.0)
    rp = _rotation_params(b, s)
    rp[3].period = 180.0                    # the ellipsoid is symmetric under W → W + 180°
    rp[3].value = mod(rp[3].value, 180)
    return vcat([Param(:a, a0; lower=1e-3, scale=sc), Param(:b, b0; lower=1e-3, scale=sc),
                 Param(:c, c0; lower=1e-3, scale=sc)], rp)
end

# the ellipsoid is unchanged by W → W + 180°, and (a, b, W) ≡ (b, a, W + 90°): canonical a ≥ b, W0 ∈ [0, 180)
function _canonical!(m::EllipsoidLimb, v, s::_Setup)
    ia, ib, iw = s.index[:a], s.index[:b], s.index[:W0]
    if v[ib] > v[ia]
        v[ia], v[ib] = v[ib], v[ia]
        v[iw] += 90
    end
    v[iw] = mod(v[iw], 180)
    return v
end

function _orientation(v, s::_Setup, k::Integer)
    return _body_matrix(_v(v, s, :pole_ra), _v(v, s, :pole_dec), _meridian(v, s, k))
end

_outline(m::EllipsoidLimb, v, k::Integer, s::_Setup) =
    project(Ellipsoid(_v(v, s, :a), _v(v, s, :b), _v(v, s, :c)), _orientation(v, s, k), s.contexts[k].u)

_shape(::EllipsoidLimb, v, s::_Setup) = Ellipsoid(_v(v, s, :a), _v(v, s, :b), _v(v, s, :c))

function _model_derived(m::EllipsoidLimb, v, s::_Setup)
    a, b, c = _v(v, s, :a), _v(v, s, :b), _v(v, s, :c)
    nev = _nev(s)
    d = Pair{Symbol,Float64}[:volume_equivalent_radius => cbrt(a * b * c), :b_over_a => b / a, :c_over_b => c / b]
    pole = radec2vec(_v(v, s, :pole_ra), _v(v, s, :pole_dec))
    for k in 1:nev
        o = _outline(m, v, k, s)
        push!(d, _pn(:apparent_a, k, nev) => o.a, _pn(:apparent_b, k, nev) => o.b,
              _pn(:apparent_pa, k, nev) => o.pa,
              _pn(:aspect, k, nev) => acosd(clamp(dot(pole, s.contexts[k].u), -1, 1)),
              _pn(:rotation_phase, k, nev) => mod(_meridian(v, s, k), 360))
    end
    return NamedTuple(d)
end

"""
    MeshLimb(mesh; matrix=nothing, fit_orientation=false, per_event=(), n=720, epoch=nothing)

Silhouette of a [`TriMesh`](@ref) (e.g. from [`read_obj`](@ref) or
[`ellipsoid_mesh`](@ref)). Parameters: `scale` (size factor of the mesh), `rotation` (in
position angle on the sky, deg; fixed at 0 unless freed, per event if in `per_event`) and
the orientation `pole_ra`, `pole_dec`, `W0`, `period` as in [`EllipsoidLimb`](@ref), free
only with `fit_orientation=true`. The orientation comes from the body's frame unless a fixed
ICRS → body `matrix` is given. The silhouette is sampled in `n` directions (the mesh must
be star-shaped as seen from its origin).
"""
struct MeshLimb <: LimbModel
    mesh::TriMesh
    matrix::Union{Nothing,Matrix{Float64}}
    fit_orientation::Bool
    per_event::Vector{Symbol}
    n::Int
    epoch::Union{Nothing,Epoch}
    cache::Dict{Tuple{Int,Float64,Float64,Float64},PolarOutline{Float64}}
end
MeshLimb(mesh::TriMesh; matrix=nothing, fit_orientation::Bool=false, per_event=(), n::Integer=720, epoch=nothing) =
    MeshLimb(mesh, matrix === nothing ? nothing : Matrix{Float64}(matrix), fit_orientation, Symbol.(collect(per_event)),
             n, epoch === nothing ? nothing : _epoch(epoch), Dict{Tuple{Int,Float64,Float64,Float64},PolarOutline{Float64}}())

function _model_params(m::MeshLimb, s::_Setup)
    nev = _nev(s)
    r0 = mean(e.radius for e in s.estimates)
    s0 = r0 / volume_equivalent_radius(m.mesh)
    ps = [Param(:scale, s0; lower=1e-6, scale=0.05 * s0)]
    if :rotation in m.per_event
        append!(ps, [Param(_pn(:rotation, k, nev), 0.0; free=false, scale=5.0) for k in 1:nev])
    else
        push!(ps, Param(:rotation, 0.0; free=false, scale=5.0))
    end
    m.matrix === nothing && append!(ps, _rotation_params(s.events[1].body, s; free=m.fit_orientation))
    return ps
end

function _outline(m::MeshLimb, v, k::Integer, s::_Setup)
    key = m.matrix === nothing ? (k, Float64(_v(v, s, :pole_ra)), Float64(_v(v, s, :pole_dec)), Float64(_meridian(v, s, k))) :
          (k, 0.0, 0.0, 0.0)
    base = get(m.cache, key, nothing)
    if base === nothing
        length(m.cache) > 2000 && empty!(m.cache)
        M = m.matrix === nothing ? _body_matrix(key[2], key[3], key[4]) : m.matrix
        base = project(m.mesh, M, s.contexts[k].u; n=m.n)
        m.cache[key] = base
    end
    return transform(base; scale=_v(v, s, :scale), rotation=_v(v, s, :rotation, k, :rotation in m.per_event))
end

_shape(m::MeshLimb, v, s::_Setup) = scale(m.mesh, _v(v, s, :scale))

function _model_derived(m::MeshLimb, v, s::_Setup)
    nev = _nev(s)
    d = Pair{Symbol,Float64}[:volume_equivalent_radius => _v(v, s, :scale) * volume_equivalent_radius(m.mesh)]
    for k in 1:nev
        push!(d, _pn(:area_equivalent_radius, k, nev) => equivalent_radius(_outline(m, v, k, s)))
    end
    return NamedTuple(d)
end
