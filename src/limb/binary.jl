# ---------------------------------------------------------------------------
# Composite limb models: a binary system or a body with a satellite (BinaryLimb), and contact
# binaries, in 2-D (ContactBinaryLimb) and in 3-D (ContactBinaryEllipsoids). Their outlines
# are CompositeOutlines; contacts can be assigned to a component (Contact(...; body=…)) or are
# matched to the nearest crossing; the visible segments of the positive chords are constraints.
# ---------------------------------------------------------------------------

# starting estimates from contact points (f, g): centre, radius, ellipse (or nothing)
function _point_estimate(f, g, R0)
    if length(f) >= 3
        c = algebraic_circle(f, g)
        ae = length(f) >= 5 ? algebraic_ellipse(f, g) : nothing
        ae !== nothing && !(0 <= ae.oblateness < 0.9) && (ae = nothing)
        return (center_f=c.center_f, center_g=c.center_g, radius=c.radius, ellipse=ae)
    elseif !isempty(f)
        r = maximum(hypot(f[i] - f[j], g[i] - g[j]) for i in eachindex(f), j in eachindex(f)) / 2
        return (center_f=mean(f), center_g=mean(g), radius=max(r, 0.6R0), ellipse=nothing)
    end
    return nothing
end

# the contact points of an event split into two groups: by their body labels (1 or 2), the
# others by 2-means clustering in the fundamental plane, seeded by the labelled points or by
# the two points farthest apart; group 1 is the group with the larger extent
function _split_points(pts)
    n = length(pts)
    f = [p.f for p in pts]; g = [p.g for p in pts]
    lab = [p.body in (:primary, :lobe1) ? 1 : p.body in (:secondary, :lobe2) ? 2 : 0 for p in pts]
    grp = copy(lab)
    n == 0 && return f, g, grp
    c = Vector{Tuple{Float64,Float64}}(undef, 2)
    for j in 1:2
        idx = findall(==(j), lab)
        c[j] = isempty(idx) ? (NaN, NaN) : (mean(f[idx]), mean(g[idx]))
    end
    if any(isnan, c[1]) || any(isnan, c[2])
        n == 1 && return f, g, max.(grp, 1)
        i, j = Tuple(argmax([hypot(f[i] - f[j], g[i] - g[j]) for i in 1:n, j in 1:n]))
        any(isnan, c[1]) && (c[1] = isnan(c[2][1]) || hypot(f[i] - c[2][1], g[i] - c[2][2]) > hypot(f[j] - c[2][1], g[j] - c[2][2]) ?
                             (f[i], g[i]) : (f[j], g[j]))
        any(isnan, c[2]) && (c[2] = c[1] == (f[i], g[i]) ? (f[j], g[j]) : (f[i], g[i]))
    end
    for _ in 1:20
        for i in 1:n
            lab[i] == 0 || continue
            grp[i] = hypot(f[i] - c[1][1], g[i] - c[1][2]) <= hypot(f[i] - c[2][1], g[i] - c[2][2]) ? 1 : 2
        end
        for j in 1:2
            idx = findall(==(j), grp)
            isempty(idx) || (c[j] = (mean(f[idx]), mean(g[idx])))
        end
    end
    # without labels, the group of larger extent is the primary
    if all(==(0), lab)
        ext(j) = (idx = findall(==(j), grp); length(idx) < 2 ? 0.0 :
                  maximum(hypot(f[a] - f[b], g[a] - g[b]) for a in idx, b in idx))
        ext(2) > ext(1) && (grp = 3 .- grp)
    end
    return f, g, grp
end

# ---------------------------------------------------------------------------
# a binary system / a body with a satellite
# ---------------------------------------------------------------------------

"""
    BinaryLimb(primary=EllipseLimb(), secondary=CircleLimb())

Two separate bodies: a binary system, or a body and its satellite. Each component is a
[`LimbModel`](@ref): by default an ellipse for the primary and a circle for the satellite,
which few chords usually cross; `BinaryLimb(EllipseLimb(), EllipseLimb())` fits an ellipse to
the satellite too (it needs at least three, better four or more chords on it). 3-D models
(`EllipsoidLimb`) can be used in multi-event fits.

**Parameters**: those of the primary (unchanged names; `center_f`, `center_g` is the centre
of the primary), those of the secondary with the prefix `sec_` (`sec_radius`, or
`sec_equatorial_radius`, `sec_oblateness`, `sec_position_angle`), and the position of the
secondary relative to the primary in each event, `sep_f`, `sep_g` [km] (suffix `_k` in
multi-event fits).

**Contacts** are assigned to a body by their label (`Contact(...; body=:primary)` or
`:secondary`, robust when the satellite is crossed by one or two chords) or, unlabelled,
matched to the nearest crossing of the same kind of either body. The starting values come
from the contacts split into two groups (labels, else clustering). The visible segments of
the positive chords are constraints by default (`visible=:soft`, see [`LimbProblem`](@ref)).

**Derived**: the derived quantities of each component (`sec_` prefix for the secondary);
`size_ratio` (ratio of the area-equivalent radii, secondary/primary), `mass_ratio`
(`size_ratio³`, equal densities), `system_equivalent_radius` (from the total area); per
event `separation` [km], `separation_mas`, `separation_pa` (position angle of the secondary
from the primary, deg) and the offsets from the ephemeris of the barycentre
(`barycenter_offset_ra`, `barycenter_offset_dec`, equal densities) and of the photocentre
(`photocenter_offset_ra`, …, equal albedos; mas), used by
[`astrometry`](@ref)`(fit; reference=:barycenter)`.
"""
struct BinaryLimb{P<:LimbModel,S<:LimbModel} <: LimbModel
    primary::P
    secondary::S
end
BinaryLimb() = BinaryLimb(EllipseLimb(), CircleLimb())

_composite(::BinaryLimb) = true
_body_index(::BinaryLimb, b::Symbol) = b in (:primary, :lobe1) ? 1 : b in (:secondary, :lobe2) ? 2 : 0

const _SEC = "sec_"

function _prepare!(m::BinaryLimb, s::_Setup)
    est2 = similar(s.estimates)
    for (k, ev) in enumerate(s.events)
        R0 = isfinite(ev.body.radius) ? ev.body.radius : 100.0
        f, g, grp = _split_points(filter(q -> !_isring(q.body), contact_points(ev)))
        e1 = _point_estimate(f[grp.==1], g[grp.==1], R0)
        e1 === nothing && (e1 = s.estimates[k])
        e2 = _point_estimate(f[grp.==2], g[grp.==2], e1.radius / 4)
        e2 === nothing && (e2 = (center_f=e1.center_f + 4 * e1.radius, center_g=e1.center_g, radius=e1.radius / 4,
                                 ellipse=nothing))
        s.estimates[k] = e1
        est2[k] = e2
    end
    s.aux[:secondary] = _Setup(s.events, s.contexts, est2, Dict{Symbol,Int}(), s.epoch)
    return nothing
end

function _finalize!(m::BinaryLimb, s::_Setup)
    sub = s.aux[:secondary]
    empty!(sub.index)
    for (n, i) in s.index
        str = string(n)
        startswith(str, _SEC) && (sub.index[Symbol(str[length(_SEC)+1:end])] = i)
    end
    # the secondary's model may use the event centres (not the case for the 2-D models)
    return nothing
end

function _model_params(m::BinaryLimb, s::_Setup)
    ps = _model_params(m.primary, s)
    sub = s.aux[:secondary]
    for p in _model_params(m.secondary, sub)
        p.name = Symbol(_SEC, p.name)
        push!(ps, p)
    end
    nev = _nev(s)
    for k in 1:nev
        e1, e2 = s.estimates[k], sub.estimates[k]
        sc = max(e2.radius / 5, 1.0)
        push!(ps, Param(_pn(:sep_f, k, nev), e2.center_f - e1.center_f; scale=sc),
              Param(_pn(:sep_g, k, nev), e2.center_g - e1.center_g; scale=sc))
    end
    return ps
end

function _outline(m::BinaryLimb, v, k::Integer, s::_Setup)
    nev = _nev(s)
    O1 = _outline(m.primary, v, k, s)
    O2 = _outline(m.secondary, v, k, s.aux[:secondary])
    df, dg = _v(v, s, _pn(:sep_f, k, nev)), _v(v, s, _pn(:sep_g, k, nev))
    return CompositeOutline([O1, O2], [(zero(df), zero(dg)), (df, dg)])
end

_admissible(m::BinaryLimb, v, s::_Setup) = _admissible(m.primary, v, s) && _admissible(m.secondary, v, s.aux[:secondary])

function _canonical!(m::BinaryLimb, v, s::_Setup)
    _canonical!(m.primary, v, s)
    _canonical!(m.secondary, v, s.aux[:secondary])
    return v
end

function _model_derived(m::BinaryLimb, v, s::_Setup)
    nev = _nev(s)
    d = Pair{Symbol,Float64}[pairs(_model_derived(m.primary, v, s))...]
    for (n, x) in pairs(_model_derived(m.secondary, v, s.aux[:secondary]))
        push!(d, Symbol(_SEC, n) => x)
    end
    O = _outline(m, v, 1, s)
    A1, A2 = area(O.parts[1]), area(O.parts[2])
    q = sqrt(A2 / A1)
    push!(d, :size_ratio => q, :mass_ratio => q^3, :system_equivalent_radius => sqrt((A1 + A2) / π))
    for k in 1:nev
        Ok = _outline(m, v, k, s)
        a1, a2 = area(Ok.parts[1]), area(Ok.parts[2])
        df, dg = _v(v, s, _pn(:sep_f, k, nev)), _v(v, s, _pn(:sep_g, k, nev))
        cf, cg = _v(v, s, _pn(:center_f, k, nev)), _v(v, s, _pn(:center_g, k, nev))
        D = s.contexts[k].dist
        μ = q^3 / (1 + q^3)                       # barycentre: fraction of the separation
        ϕ = a2 / (a1 + a2)                        # photocentre (equal albedos)
        push!(d, _pn(:separation, k, nev) => hypot(df, dg),
              _pn(:separation_mas, k, nev) => hypot(df, dg) / D * RAD2MAS,
              _pn(:separation_pa, k, nev) => mod(atand(df, dg), 360),
              _pn(:barycenter_offset_ra, k, nev) => (cf + μ * df) / D * RAD2MAS,
              _pn(:barycenter_offset_dec, k, nev) => (cg + μ * dg) / D * RAD2MAS,
              _pn(:photocenter_offset_ra, k, nev) => (cf + ϕ * df) / D * RAD2MAS,
              _pn(:photocenter_offset_dec, k, nev) => (cg + ϕ * dg) / D * RAD2MAS)
    end
    return NamedTuple(d)
end

# ---------------------------------------------------------------------------
# a contact binary seen in projection: two ellipses along a common axis
# ---------------------------------------------------------------------------

"""
    ContactBinaryLimb(; per_event=())

Apparent limb of a contact binary: two ellipses (the lobes) whose axes `a1`, `a2` [km] lie
along a common axis at the position angle `axis_pa` [deg, north towards east, from lobe 1
to lobe 2], with the semi-axes `b1`, `b2` across it. The distance between their centres is
`d = (a1 + a2)(1 − neck)`: `neck = 0` means the lobes just touch, `neck > 0` that they
overlap (a waist), `neck < 0` (down to −0.5) a gap between them. `center_f`, `center_g` is
the area-weighted centre of the two lobes (the centre of light for equal albedos, overlap
ignored). Lobe 1 is the larger one (canonical form). With `per_event=(:axis_pa,)` the axis
has its own position angle in every event (a rotating body seen at the same aspect); use
[`ContactBinaryEllipsoids`](@ref) for the 3-D shape from several events.

**Derived**: `equivalent_radius` (area of the union), `lobe1_equivalent_radius`,
`lobe2_equivalent_radius`, `size_ratio` (lobe 2/lobe 1), `mass_ratio` (`size_ratio³`),
`separation` `d`, `length` (total extent along the axis) and `neck_width` (width of the
waist, 0 for detached lobes).
"""
struct ContactBinaryLimb <: LimbModel
    per_event::Vector{Symbol}
end
function ContactBinaryLimb(; per_event=())
    pe = Symbol.(collect(per_event))
    all(==(:axis_pa), pe) || throw(ArgumentError("per_event may contain :axis_pa"))
    return ContactBinaryLimb(pe)
end

_composite(::ContactBinaryLimb) = true
_body_index(::ContactBinaryLimb, b::Symbol) = b in (:primary, :lobe1) ? 1 : b in (:secondary, :lobe2) ? 2 : 0

function _model_params(m::ContactBinaryLimb, s::_Setup)
    nev = _nev(s)
    est = s.estimates
    # the overall ellipse of the contacts: total length along its major axis
    function axis0(k)
        e = est[k].ellipse
        e !== nothing && return e.equatorial_radius, e.equatorial_radius * (1 - e.oblateness), mod(e.position_angle + 90, 360)
        f, g, grp = _split_points(contact_points(s.events[k]))
        (count(==(1), grp) > 0 && count(==(2), grp) > 0) || return est[k].radius, est[k].radius, 0.0
        c1 = (mean(f[grp.==1]), mean(g[grp.==1])); c2 = (mean(f[grp.==2]), mean(g[grp.==2]))
        return est[k].radius, 0.7est[k].radius, mod(atand(c2[1] - c1[1], c2[2] - c1[2]), 360)
    end
    ax = [axis0(k) for k in 1:nev]
    A = mean(first.(ax)); B = mean(getindex.(ax, 2))
    a0 = A / 1.9                                  # neck 0.1: d + a1 + a2 = 2A
    sc = max(a0 / 10, 1.0)
    ps = [Param(:a1, a0; lower=1e-3, scale=sc), Param(:b1, B; lower=1e-3, scale=sc),
          Param(:a2, a0; lower=1e-3, scale=sc), Param(:b2, B; lower=1e-3, scale=sc)]
    if :axis_pa in m.per_event
        append!(ps, [Param(_pn(:axis_pa, k, nev), ax[k][3]; scale=10.0, period=360.0) for k in 1:nev])
    else
        push!(ps, Param(:axis_pa, ax[1][3]; scale=10.0, period=360.0))
    end
    push!(ps, Param(:neck, 0.1; lower=-0.5, upper=0.95, scale=0.05))
    return ps
end

# lobe centres along the axis (area-weighted centre at 0) and the separation
function _lobes_2d(a1, b1, a2, b2, η)
    d = (a1 + a2) * (1 - η)
    A1, A2 = a1 * b1, a2 * b2
    return -d * A2 / (A1 + A2), d * A1 / (A1 + A2), d
end

function _outline(m::ContactBinaryLimb, v, k::Integer, s::_Setup)
    a1, b1, a2, b2, η = _v(v, s, :a1), _v(v, s, :b1), _v(v, s, :a2), _v(v, s, :b2), _v(v, s, :neck)
    pa = _v(v, s, :axis_pa, k, :axis_pa in m.per_event)
    x1, x2, _ = _lobes_2d(a1, b1, a2, b2, η)
    sn, cs = sincosd(pa)                          # the axis: (f, g) = (sin, cos) PA
    return CompositeOutline([EllipseOutline(a1, b1, pa + 90), EllipseOutline(a2, b2, pa + 90)],
                            [(x1 * sn, x1 * cs), (x2 * sn, x2 * cs)])
end

# lobe 1 is the larger: otherwise swap the lobes and turn the axis by 180°
function _canonical!(m::ContactBinaryLimb, v, s::_Setup)
    i = s.index
    if v[i[:a2]] * v[i[:b2]] > v[i[:a1]] * v[i[:b1]]
        v[i[:a1]], v[i[:a2]] = v[i[:a2]], v[i[:a1]]
        v[i[:b1]], v[i[:b2]] = v[i[:b2]], v[i[:b1]]
        for n in keys(i)
            startswith(string(n), "axis_pa") && (v[i[n]] += 180)
        end
    end
    for n in keys(i)
        startswith(string(n), "axis_pa") && (v[i[n]] = mod(v[i[n]], 360))
    end
    return v
end

# half-width of the waist: where the two lobe boundaries cross (lobe 1 at x1, lobe 2 at x2 on the axis)
function _waist(a1, b1, x1, a2, b2, x2)
    h1(x) = b1 * sqrt(max(1 - ((x - x1) / a1)^2, 0.0))
    h2(x) = b2 * sqrt(max(1 - ((x - x2) / a2)^2, 0.0))
    lo, hi = max(x1 - a1, x2 - a2), min(x1 + a1, x2 + a2)
    lo < hi || return 0.0
    # h1 - h2 decreases from lo to hi (lobe 1 left of lobe 2)
    f(x) = h1(x) - h2(x)
    f(lo) * f(hi) > 0 && return min(h1((lo + hi) / 2), h2((lo + hi) / 2))
    for _ in 1:60
        mid = (lo + hi) / 2
        (f(mid) > 0) == (f(lo) > 0) ? (lo = mid) : (hi = mid)
    end
    return h1((lo + hi) / 2)
end

function _model_derived(m::ContactBinaryLimb, v, s::_Setup)
    a1, b1, a2, b2, η = _v(v, s, :a1), _v(v, s, :b1), _v(v, s, :a2), _v(v, s, :b2), _v(v, s, :neck)
    x1, x2, d = _lobes_2d(a1, b1, a2, b2, η)
    q = sqrt(a2 * b2 / (a1 * b1))
    O = _outline(m, v, 1, s)
    return (equivalent_radius=sqrt(area(O) / π), lobe1_equivalent_radius=sqrt(a1 * b1),
            lobe2_equivalent_radius=sqrt(a2 * b2), size_ratio=q, mass_ratio=q^3, separation=d,
            length=max(x1 + a1, x2 + a2) - min(x1 - a1, x2 - a2), neck_width=2_waist(a1, b1, x1, a2, b2, x2))
end

# ---------------------------------------------------------------------------
# a contact binary in 3-D: two ellipsoids, rotating
# ---------------------------------------------------------------------------

"""
    ContactBinaryEllipsoids(; epoch=nothing)

3-D contact binary: two triaxial ellipsoids (the lobes, semi-axes `a1`, `b1`, `c1` and `a2`,
`b2`, `c2` [km]) on the body x-axis with parallel axes, at the distance
`d = (a1 + a2)(1 − neck)` (see [`BilobedEllipsoid`](@ref); `neck` > 0: overlapping lobes),
rotating about the body z-axis with the pole `pole_ra`, `pole_dec`, the prime meridian `W0`
at the reference `epoch` and the `period` [h] as [`EllipsoidLimb`](@ref). The centre
`center_f`, `center_g` is the centre of mass (equal densities). Each event sees the exact
silhouette, the union of the two projected ellipses; several events at different rotational
phases, rotational amplitudes ([`RotationalAmplitude`](@ref)) and phases
([`RotationalPhase`](@ref)) constrain the shape. Starting values come from the body's
`shape` if it is a [`BilobedEllipsoid`](@ref), else from the contacts. Lobe 1 is the larger
one (canonical form).

**Derived**: `volume_equivalent_radius` (of the union), `lobe1_volume_equivalent_radius`,
`lobe2_volume_equivalent_radius`, `mass_ratio` (V₂/V₁), `separation`, `length` (along x);
per event the area-equivalent radius `area_equivalent_radius`, the `aspect` and the
`rotation_phase`.
"""
struct ContactBinaryEllipsoids <: LimbModel
    epoch::Union{Nothing,Epoch}
end
ContactBinaryEllipsoids(; epoch=nothing) = ContactBinaryEllipsoids(epoch === nothing ? nothing : _epoch(epoch))

_composite(::ContactBinaryEllipsoids) = true
_body_index(::ContactBinaryEllipsoids, b::Symbol) = b in (:primary, :lobe1) ? 1 : b in (:secondary, :lobe2) ? 2 : 0

function _model_params(m::ContactBinaryEllipsoids, s::_Setup)
    b = s.events[1].body
    r0 = mean(e.radius for e in s.estimates)
    if b.shape isa BilobedEllipsoid
        sh = b.shape
        l1, l2 = sh.lobe1, sh.lobe2
        η0 = 1 - sh.d / (l1.a + l2.a)
        init = (l1.a, l1.b, l1.c, l2.a, l2.b, l2.c)
    else
        η0 = 0.1
        init = (0.6r0, 0.5r0, 0.45r0, 0.5r0, 0.45r0, 0.4r0)
    end
    sc = max(r0 / 10, 1.0)
    ps = [Param(n, x; lower=1e-3, scale=sc) for (n, x) in zip((:a1, :b1, :c1, :a2, :b2, :c2), init)]
    push!(ps, Param(:neck, η0; lower=-0.5, upper=0.95, scale=0.05))
    return vcat(ps, _rotation_params(b, s))
end

_shape(::ContactBinaryEllipsoids, v, s::_Setup) =
    BilobedEllipsoid(Ellipsoid(_v(v, s, :a1), _v(v, s, :b1), _v(v, s, :c1)),
                     Ellipsoid(_v(v, s, :a2), _v(v, s, :b2), _v(v, s, :c2)),
                     (_v(v, s, :a1) + _v(v, s, :a2)) * (1 - _v(v, s, :neck)))

_outline(m::ContactBinaryEllipsoids, v, k::Integer, s::_Setup) = project(_shape(m, v, s), _orientation(v, s, k), s.contexts[k].u)

# lobe 1 is the larger: otherwise swap the lobes and turn by 180° about the pole
function _canonical!(m::ContactBinaryEllipsoids, v, s::_Setup)
    i = s.index
    if v[i[:a2]] * v[i[:b2]] * v[i[:c2]] > v[i[:a1]] * v[i[:b1]] * v[i[:c1]]
        for (p, q) in ((:a1, :a2), (:b1, :b2), (:c1, :c2))
            v[i[p]], v[i[q]] = v[i[q]], v[i[p]]
        end
        v[i[:W0]] += 180
    end
    v[i[:W0]] = mod(v[i[:W0]], 360)
    return v
end

function _model_derived(m::ContactBinaryEllipsoids, v, s::_Setup)
    sh = _shape(m, v, s)
    nev = _nev(s)
    V1, V2 = volume(sh.lobe1), volume(sh.lobe2)
    x1, x2 = _lobe_x(sh)
    d = Pair{Symbol,Float64}[:volume_equivalent_radius => volume_equivalent_radius(sh),
                             :lobe1_volume_equivalent_radius => cbrt(3V1 / 4π),
                             :lobe2_volume_equivalent_radius => cbrt(3V2 / 4π), :mass_ratio => V2 / V1,
                             :separation => sh.d,
                             :length => max(x1 + sh.lobe1.a, x2 + sh.lobe2.a) - min(x1 - sh.lobe1.a, x2 - sh.lobe2.a)]
    pole = radec2vec(_v(v, s, :pole_ra), _v(v, s, :pole_dec))
    for k in 1:nev
        push!(d, _pn(:area_equivalent_radius, k, nev) => equivalent_radius(_outline(m, v, k, s)),
              _pn(:aspect, k, nev) => acosd(clamp(dot(pole, s.contexts[k].u), -1, 1)),
              _pn(:rotation_phase, k, nev) => mod(_meridian(v, s, k), 360))
    end
    return NamedTuple(d)
end
