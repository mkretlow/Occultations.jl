# ---------------------------------------------------------------------------
# Composite outlines: the union of several outlines with their own centres, e.g. a body and
# its satellite, or the two lobes of a contact binary. The union is in general not
# star-shaped: crossings of a track are a set of intervals, and the boundary is a set of
# polylines.
# ---------------------------------------------------------------------------

"""
    CompositeOutline(parts, offsets)

Union of the outlines `parts`, the `k`-th centred at `offsets[k] = (f, g)` [km] relative to
the reference centre of the composite (for [`BinaryLimb`](@ref) the centre of the primary,
for the contact-binary models the centre of mass). Produced by the binary and
contact-binary limb models.

* [`signed_distance`](@ref) is the smallest signed distance from the parts (negative inside
  any of them), so negative chords work unchanged;
* [`intervals`](@ref) gives the occulted intervals of a straight track (exact for ellipses);
  [`crossing`](@ref) the first entry or the last exit;
* [`area`](@ref) is the area of the union (overlaps counted once);
* [`outline_points`](@ref) the boundary of the union: the parts of each component outline
  outside the others, separated by `NaN`;
* `limb_radius(o, θ)` is the outermost limb point in the direction θ (the union is not
  star-shaped in general: use it only where this is meant).
"""
struct CompositeOutline{P<:Tuple,O<:Tuple} <: Outline
    parts::P
    offsets::O
    function CompositeOutline(parts, offsets)
        length(parts) == length(offsets) || throw(DimensionMismatch("one offset per part"))
        isempty(parts) && throw(ArgumentError("a composite outline needs at least one part"))
        all(p -> p isa Outline, parts) || throw(ArgumentError("the parts must be outlines"))
        ps = Tuple(parts)
        os = Tuple((o[1], o[2]) for o in offsets)
        return new{typeof(ps),typeof(os)}(ps, os)
    end
end

function Base.show(io::IO, c::CompositeOutline)
    print(io, "CompositeOutline(")
    for (k, (p, o)) in enumerate(zip(c.parts, c.offsets))
        k > 1 && print(io, ", ")
        print(io, p)
        @printf(io, " at (%.2f, %.2f) km", o...)
    end
    print(io, ")")
end

signed_distance(c::CompositeOutline, x::Real, y::Real) =
    minimum(signed_distance(p, x - o[1], y - o[2]) for (p, o) in zip(c.parts, c.offsets))

# index of the part nearest to (x, y) (smallest signed distance)
_nearest_part(c::CompositeOutline, x, y) = argmin([signed_distance(p, x - o[1], y - o[2]) for (p, o) in zip(c.parts, c.offsets)])

function limb_normal(c::CompositeOutline, x::Real, y::Real)
    k = _nearest_part(c, x, y)
    o = c.offsets[k]
    return limb_normal(c.parts[k], x - o[1], y - o[2])
end

# outermost limb point along the ray from the reference centre
function limb_radius(c::CompositeOutline, θ::Real)
    s, w = sincos(θ)
    r = zero(promote_type(typeof(θ), typeof(c.offsets[1][1])))
    for (p, o) in zip(c.parts, c.offsets)
        iv = _interval(p, -o[1], -o[2], w, s)
        iv === nothing || (r = max(r, iv[2]))
    end
    return r
end

_rmax(c::CompositeOutline) = maximum(hypot(o...) + _rmax(p) for (p, o) in zip(c.parts, c.offsets))

# ---------------------------------------------------------------------------
# intervals of a straight track inside an outline
# ---------------------------------------------------------------------------

# (s_in, s_out) of the track (x, y) + s (vx, vy) inside the outline, or nothing
function _interval(e::EllipseOutline, x, y, vx, vy)
    u0, w0 = _ellipse_uv(e, x, y)
    du, dw = _ellipse_uv(e, vx, vy)
    A = (du / e.a)^2 + (dw / e.b)^2
    B = 2 * (u0 * du / e.a^2 + w0 * dw / e.b^2)
    C = (u0 / e.a)^2 + (w0 / e.b)^2 - 1
    disc = B^2 - 4A * C
    disc > 0 || return nothing
    q = sqrt(disc)
    return ((-B - q) / (2A), (-B + q) / (2A))
end

function _interval(o::Outline, x, y, vx, vy)
    s1, m1 = crossing(o, x, y, vx, vy, true)
    m1 > 0 && return nothing
    s2, _ = crossing(o, x, y, vx, vy, false)
    return s2 > s1 ? (s1, s2) : nothing
end

function _interval(c::CompositeOutline, x, y, vx, vy)
    iv = intervals(c, x, y, vx, vy)
    return isempty(iv) ? nothing : (first(iv)[1], last(iv)[2])
end

"""
    intervals(outline, x, y, vx, vy) -> Vector{Tuple}

Intervals `(s_in, s_out)` of the straight track `(x, y) + s (vx, vy)` inside the outline,
sorted and merged (for a [`CompositeOutline`](@ref), the union of the intervals of its parts:
a chord crossing two separated bodies has two intervals, one crossing two overlapping lobes
one). Empty if the track misses the outline. With `(vx, vy)` the shadow velocity [km/s],
`s` is a time [s].

# Examples
```jldoctest
julia> c = CompositeOutline([EllipseOutline(100.0, 100.0, 0.0), EllipseOutline(50.0, 50.0, 0.0)],
                            [(0.0, 0.0), (300.0, 0.0)]);

julia> [round.(iv; digits=6) for iv in intervals(c, -500.0, 0.0, 10.0, 0.0)]   # 10 km/s along f: two drops
2-element Vector{Tuple{Float64, Float64}}:
 (40.0, 60.0)
 (75.0, 85.0)
```
"""
function intervals(o::Outline, x, y, vx, vy)
    iv = _interval(o, x, y, vx, vy)
    T = promote_type(typeof(x), typeof(y), typeof(vx), typeof(vy), Float64)
    return iv === nothing ? Tuple{T,T}[] : [convert(Tuple{T,T}, iv)]
end

function intervals(c::CompositeOutline, x, y, vx, vy)
    raw = [iv for (p, o) in zip(c.parts, c.offsets) for iv in (_interval(p, x - o[1], y - o[2], vx, vy),)
           if iv !== nothing]
    T = promote_type(typeof(x), typeof(y), typeof(vx), typeof(vy), Float64,
                     (typeof(v) for iv in raw for v in iv)...)
    out = Tuple{T,T}[]
    for iv in sort(raw; by=first)
        if !isempty(out) && iv[1] <= out[end][2]
            out[end] = (out[end][1], max(out[end][2], iv[2]))
        else
            push!(out, (iv[1], iv[2]))
        end
    end
    return out
end

function crossing(c::CompositeOutline, x, y, vx, vy, entry::Bool)
    iv = intervals(c, x, y, vx, vy)
    if !isempty(iv)
        s = entry ? first(iv)[1] : last(iv)[2]
        return s, zero(s)
    end
    # missed: the closest approach to the nearest part
    best = nothing
    for (p, o) in zip(c.parts, c.offsets)
        r = crossing(p, x - o[1], y - o[2], vx, vy, entry)
        (best === nothing || r[2] < best[2]) && (best = r)
    end
    return best
end

# crossing of a contact with a composite outline: with a part index `body` > 0 the crossing of
# that part; otherwise the crossing of the union of the same kind (entry/exit) nearest to s = 0
function _contact_crossing(O::Outline, x, y, vx, vy, entry::Bool, body::Integer)
    O isa CompositeOutline || return crossing(O, x, y, vx, vy, entry)
    if body > 0
        o = O.offsets[body]
        return crossing(O.parts[body], x - o[1], y - o[2], vx, vy, entry)
    end
    iv = intervals(O, x, y, vx, vy)
    isempty(iv) && return crossing(O, x, y, vx, vy, entry)
    ss = [entry ? a : b for (a, b) in iv]
    s = ss[argmin(abs.(ss))]
    return s, zero(s)
end

# signed distance from one part of a composite (body > 0) or from the outline
_part_distance(O::Outline, x, y, body::Integer) =
    O isa CompositeOutline && body > 0 ? signed_distance(O.parts[body], x - O.offsets[body][1], y - O.offsets[body][2]) :
    signed_distance(O, x, y)

# ---------------------------------------------------------------------------
# area and boundary of the union
# ---------------------------------------------------------------------------

# area of the union: the sum of the areas less the overlaps, the latter integrated along g over
# the lines where the parts overlap (n lines)
function area(c::CompositeOutline; n::Integer=600)
    A = sum(area(p) for p in c.parts)
    length(c.parts) == 1 && return A
    # g range in which at least two parts overlap
    ext = [(o[2] - _rmax(p), o[2] + _rmax(p)) for (p, o) in zip(c.parts, c.offsets)]
    lo, hi = Inf, -Inf
    for i in eachindex(ext), j in i+1:length(ext)
        a, b = max(ext[i][1], ext[j][1]), min(ext[i][2], ext[j][2])
        a < b && (lo = min(lo, a); hi = max(hi, b))
    end
    lo < hi || return A
    xmin = minimum(o[1] - _rmax(p) for (p, o) in zip(c.parts, c.offsets)) - 1
    h = (hi - lo) / n
    over = zero(A)
    length(c.parts) == 2 && return A - _overlap2(c.parts[1], c.parts[2], c.offsets[1], c.offsets[2], xmin, lo, h, n, A)
    for i in 1:n
        y = lo + (i - 0.5) * h                         # midpoint rule (no weight at the tangent points)
        tot = zero(A)
        for (p, o) in zip(c.parts, c.offsets)
            iv = _interval(p, xmin - o[1], y - o[2], one(A), zero(A))
            iv === nothing || (tot += iv[2] - iv[1])
        end
        tot > 0 || continue
        uni = sum(b - a for (a, b) in intervals(c, xmin, y, one(A), zero(A)); init=zero(A))
        over += tot - uni
    end
    return A - over * h
end

# two parts: the summed overlap of their intervals on n lines (a function barrier on the concrete
# outline types, no allocations)
function _overlap2(p1, p2, o1, o2, xmin, lo, h, n, A)
    over = zero(A)
    for i in 1:n
        y = lo + (i - 0.5) * h
        a = _interval(p1, xmin - o1[1], y - o1[2], one(A), zero(A))
        a === nothing && continue
        b = _interval(p2, xmin - o2[1], y - o2[2], one(A), zero(A))
        b === nothing && continue
        over += max(min(a[2], b[2]) - max(a[1], b[1]), zero(A))
    end
    return over * h
end

function outline_points(c::CompositeOutline; center=(0.0, 0.0), n::Integer=361)
    fs, gs = Float64[], Float64[]
    for (k, (p, o)) in enumerate(zip(c.parts, c.offsets))
        pt(θ) = (r = limb_radius(p, θ); (o[1] + r * cos(θ), o[2] + r * sin(θ)))
        vis(xy) = all(j == k || signed_distance(c.parts[j], xy[1] - c.offsets[j][1], xy[2] - c.offsets[j][2]) >= 0
                      for j in eachindex(c.parts))
        θs = collect(range(-π, π; length=n))[1:end-1]
        out = [vis(pt(θ)) for θ in θs]
        # the boundary point between a visible and a hidden direction (bisection)
        function edge(θa, θb)                          # θa visible, θb hidden
            for _ in 1:40
                θm = (θa + θb) / 2
                vis(pt(θm)) ? (θa = θm) : (θb = θm)
            end
            return pt(θa)
        end
        # walk once round from a hidden direction, so that every visible arc is one piece
        m = length(θs)
        i0 = something(findfirst(!, out), 1)
        h = 2π / m
        isempty(fs) || (push!(fs, NaN); push!(gs, NaN))
        piece = false
        for q in 1:m+1
            j = mod1(i0 + q - 1, m)
            θ = θs[i0] + (q - 1) * h
            if out[j]
                if !piece
                    (isempty(fs) || isnan(fs[end])) || (push!(fs, NaN); push!(gs, NaN))
                    if !out[mod1(j - 1, m)]
                        e = edge(θ, θ - h)
                        push!(fs, center[1] + e[1]); push!(gs, center[2] + e[2])
                    end
                    piece = true
                end
                x, y = pt(θ)
                push!(fs, center[1] + x); push!(gs, center[2] + y)
            elseif piece
                e = edge(θ - h, θ)
                push!(fs, center[1] + e[1]); push!(gs, center[2] + e[2])
                piece = false
            end
        end
    end
    return fs, gs
end
