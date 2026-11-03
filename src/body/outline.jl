# ---------------------------------------------------------------------------
# Outlines: the apparent limb of a body on the sky plane, in (f, g) km relative to its
# centre (f towards east, g towards north). Every limb model produces outlines; all
# fitting code works only through this interface.
# ---------------------------------------------------------------------------

"""
    Outline

Apparent limb on the sky plane, relative to the body centre, in km. Implementations:
[`EllipseOutline`](@ref) (analytic) and [`PolarOutline`](@ref) (tabulated, e.g. the
silhouette of a mesh). Interface: [`limb_radius`](@ref), [`inside`](@ref),
[`signed_distance`](@ref), [`outline_points`](@ref), [`area`](@ref).
"""
abstract type Outline end

"""
    EllipseOutline(a, b, pa)
    ellipse_outline(; equatorial_radius, oblateness=0, position_angle=0)

Ellipse with semi-major axis `a`, semi-minor axis `b` [km]; `pa` is the position angle of
the **semi-minor** axis (the apparent pole for an oblate body), from north towards east
[deg]. `oblateness = 1 - b/a`.
"""
struct EllipseOutline{T<:Real} <: Outline
    a::T
    b::T
    pa::T
end
EllipseOutline(a, b, pa) = EllipseOutline(promote(a, b, pa)...)

"""
    ellipse_outline(; equatorial_radius, oblateness=0.0, position_angle=0.0) -> EllipseOutline

An [`EllipseOutline`](@ref) from the parameters of [`EllipseLimb`](@ref): semi-major axis
`equatorial_radius` [km], `oblateness` ε = 1 − b/a and the `position_angle` of the minor
axis [deg, north towards east]. Useful to test a published ellipse against the chords, e.g.
with [`check_negative_chords`](@ref).

# Examples
```jldoctest
julia> e = ellipse_outline(equatorial_radius=250.0, oblateness=0.2)
EllipseOutline(a = 250.000 km, b = 200.000 km, PA(minor) = 0.00°)

julia> limb_radius(e, 0.0), limb_radius(e, π/2)          # along f (east) and g (north)
(250.0, 200.0)

julia> equivalent_radius(e) ≈ sqrt(250 * 200)
true

julia> crossing(e, -1000.0, 0.0, 1.0, 0.0, true)         # entry of a track along f: s = 750 km
(750.0, 0.0)

julia> s, miss = crossing(e, -1000.0, 300.0, 1.0, 0.0, true); round(miss; digits=6)   # misses by 100 km
100.0
```
"""
ellipse_outline(; equatorial_radius, oblateness=0.0, position_angle=0.0) =
    EllipseOutline(equatorial_radius, equatorial_radius * (1 - oblateness), position_angle)

Base.show(io::IO, e::EllipseOutline) =
    @printf(io, "EllipseOutline(a = %.3f km, b = %.3f km, PA(minor) = %.2f°)", e.a, e.b, e.pa)

# (major, minor)-axis coordinates of (x, y) = (f, g)
@inline function _ellipse_uv(e::EllipseOutline, x, y)
    s, c = sincosd(e.pa)
    return x * c - y * s, x * s + y * c
end

"""
    PolarOutline(θ, r; scale=1, rotation=0)

Star-shaped outline given by its radius `r[k]` [km] in the directions `θ[k] = atan(g, f)`
(uniform grid from -π), multiplied by `scale` and rotated by `rotation` [deg] in position
angle (from north towards east). Produced by [`project`](@ref) for meshes.
"""
struct PolarOutline{T<:Real} <: Outline
    θ::Vector{Float64}
    r::Vector{Float64}
    scale::T
    rotation::T
end
PolarOutline(θ, r; scale=1.0, rotation=0.0) = PolarOutline(collect(Float64, θ), collect(Float64, r), promote(scale, rotation)...)
transform(o::PolarOutline; scale=1.0, rotation=0.0) = PolarOutline(o.θ, o.r, promote(o.scale * scale, o.rotation + rotation)...)

Base.show(io::IO, o::PolarOutline) = @printf(io, "PolarOutline(%d directions, r = %.2f – %.2f km)",
                                             length(o.θ), (o.scale .* extrema(o.r))...)

"""
    limb_radius(outline, θ) -> km
    limb_radius(outline, x, y)

Radius of the outline in the direction `θ = atan(g, f)` [rad] (or of the point (x, y)).
"""
function limb_radius(e::EllipseOutline, θ::Real)
    u, w = _ellipse_uv(e, cos(θ), sin(θ))
    return 1 / sqrt((u / e.a)^2 + (w / e.b)^2)
end

function limb_radius(o::PolarOutline, θ::Real)
    n = length(o.θ)
    # a rotation by +ρ in position angle moves the direction θ to θ - ρ
    x = (θ + o.rotation * DEG2RAD - o.θ[1]) / (2π / n)
    k = floor(Int, _value(x))
    w = x - k
    return o.scale * ((1 - w) * o.r[mod1(k + 1, n)] + w * o.r[mod1(k + 2, n)])
end
limb_radius(o::Outline, x::Real, y::Real) = limb_radius(o, atan(y, x))

_value(x::Real) = x                          # extended for dual numbers if needed

"`true` if (x, y) lies inside the outline."
inside(o::Outline, x::Real, y::Real) = signed_distance(o, x, y) < 0

"""
    signed_distance(outline, x, y) -> km

Radial distance of (x, y) from the outline: `hypot(x, y) - limb_radius(outline, x, y)`
(negative inside).
"""
signed_distance(o::Outline, x::Real, y::Real) = hypot(x, y) - limb_radius(o, x, y)

# elliptical "level" (< 1 inside) of an ellipse
@inline function _level(e::EllipseOutline, x, y)
    u, w = _ellipse_uv(e, x, y)
    return (u / e.a)^2 + (w / e.b)^2
end

"""
    limb_normal(outline, x, y) -> (nx, ny)

Outward unit normal of the outline at the limb point in the direction of (x, y).
"""
function limb_normal(e::EllipseOutline, x::Real, y::Real)
    u, w = _ellipse_uv(e, x, y)
    s, c = sincosd(e.pa)
    du, dw = u / e.a^2, w / e.b^2               # gradient in (u, w)
    nx, ny = du * c + dw * s, -du * s + dw * c
    h = hypot(nx, ny)
    return nx / h, ny / h
end

function limb_normal(o::Outline, x::Real, y::Real)
    θ = atan(y, x)
    h = 1e-4
    r = limb_radius(o, θ)
    dr = (limb_radius(o, θ + h) - limb_radius(o, θ - h)) / (2h)
    # limb point P(θ) = r(θ) (cos θ, sin θ); tangent T = dP/dθ; outward normal ⟂ T
    tx, ty = dr * cos(θ) - r * sin(θ), dr * sin(θ) + r * cos(θ)
    h2 = hypot(tx, ty)
    return ty / h2, -tx / h2
end

"""
    area(outline) -> km²
"""
area(e::EllipseOutline) = π * e.a * e.b
area(o::PolarOutline) = o.scale^2 * sum(abs2, o.r) * π / length(o.r)

"Radius of the circle with the same area as the outline."
equivalent_radius(o::Outline) = sqrt(area(o) / π)

"""
    outline_points(outline; center=(0, 0), n=361) -> (f, g)

Closed polyline of the outline (for plotting).
"""
function outline_points(o::Outline; center=(0.0, 0.0), n::Integer=361)
    θ = range(-π, π; length=n)
    r = [limb_radius(o, t) for t in θ]
    return center[1] .+ r .* cos.(θ), center[2] .+ r .* sin.(θ)
end

# ---------------------------------------------------------------------------
# crossings of a straight track with an outline
# ---------------------------------------------------------------------------

"""
    crossing(outline, x, y, vx, vy, entry::Bool) -> (s, miss)

Crossing of the straight track `(x, y) + s (vx, vy)` with the outline: the offset `s` of
the entry (`entry = true`, immersion) or exit (emersion) point. If the track misses the
outline, `s` is the offset of the closest approach and `miss` (> 0) the radial distance
by which it misses [km]; otherwise `miss = 0`.
"""
function crossing(e::EllipseOutline, x, y, vx, vy, entry::Bool)
    u0, w0 = _ellipse_uv(e, x, y)
    du, dw = _ellipse_uv(e, vx, vy)
    A = (du / e.a)^2 + (dw / e.b)^2
    B = 2 * (u0 * du / e.a^2 + w0 * dw / e.b^2)
    C = (u0 / e.a)^2 + (w0 / e.b)^2 - 1
    disc = B^2 - 4A * C
    if disc >= 0
        q = sqrt(disc)
        return (entry ? (-B - q) : (-B + q)) / (2A), zero(disc)
    end
    s = -B / (2A)
    ρ = sqrt(C - B^2 / (4A) + 1)                    # elliptical radius at closest approach
    return s, hypot(x + s * vx, y + s * vy) * (1 - 1 / ρ)
end

function crossing(o::Outline, x, y, vx, vy, entry::Bool; n::Integer=81)
    v = hypot(vx, vy)
    S = (hypot(x, y) + 2 * _rmax(o)) / v
    h(s) = signed_distance(o, x + s * vx, y + s * vy)
    ss = range(-S, S; length=n)
    hs = h.(ss)
    k = entry ? findfirst(i -> hs[i] >= 0 && hs[i+1] < 0, 1:n-1) :
                findlast(i -> hs[i] < 0 && hs[i+1] >= 0, 1:n-1)
    if k !== nothing
        lo, hi = ss[k], ss[k+1]
        hlo = hs[k]
        for _ in 1:45
            m = (lo + hi) / 2
            hm = h(m)
            (hm >= 0) == (hlo >= 0) ? ((lo, hlo) = (m, hm)) : (hi = m)
        end
        return (lo + hi) / 2, 0.0
    end
    i = argmin(hs)
    hs[i] < 0 && return ss[i], 0.0
    s, hmin = _golden_min(h, ss[max(i - 1, 1)], ss[min(i + 1, n)]; tol=1e-7 * S)
    return s, max(hmin, 0.0)
end

_rmax(e::EllipseOutline) = max(e.a, e.b)
_rmax(o::PolarOutline) = o.scale * maximum(o.r)

function _golden_min(fun, a, b; tol=1e-6, maxiter=200)
    φ = (sqrt(5) - 1) / 2
    c = b - φ * (b - a); d = a + φ * (b - a)
    fc, fd = fun(c), fun(d)
    for _ in 1:maxiter
        abs(b - a) < tol && break
        if fc < fd
            b, d, fd = d, c, fc
            c = b - φ * (b - a); fc = fun(c)
        else
            a, c, fc = c, d, fd
            d = a + φ * (b - a); fd = fun(d)
        end
    end
    x = (a + b) / 2
    return x, fun(x)
end

"""
    algebraic_ellipse(f, g) -> NamedTuple or nothing

Direct least-squares conic through at least 5 points (data centred and scaled; smallest
singular vector of the design matrix). Returns `(center_f, center_g, equatorial_radius,
oblateness, position_angle)` or `nothing` if the conic is not an ellipse. A quick look
without uncertainties; used as a starting point of the limb fits.
"""
function algebraic_ellipse(f::AbstractVector{<:Real}, g::AbstractVector{<:Real})
    length(f) == length(g) || throw(DimensionMismatch("f and g lengths differ"))
    length(f) < 5 && return nothing
    f0, g0 = mean(f), mean(g)
    s = sqrt(mean(abs2, f .- f0) + mean(abs2, g .- g0))
    s > 0 || return nothing
    x, y = (f .- f0) ./ s, (g .- g0) ./ s
    D = hcat(x .^ 2, x .* y, y .^ 2, x, y, ones(length(x)))
    rank(D) >= 5 || return nothing
    v = svd(D; full=true).V[:, end]
    Q = [v[1] v[2]/2; v[2]/2 v[3]]
    det(Q) > 1e-12 * sum(abs2, Q) || return nothing
    c = -(Q \ v[4:5]) ./ 2
    k = dot(c, Q * c) - v[6]
    E = eigen(Symmetric(Q))
    r2 = k ./ E.values
    all(z -> isfinite(z) && z > 0, r2) || return nothing
    i = argmin(r2)
    a, b = sqrt(maximum(r2)) * s, sqrt(r2[i]) * s
    w = E.vectors[:, i]
    return (center_f=f0 + s * c[1], center_g=g0 + s * c[2], equatorial_radius=a, oblateness=1 - b / a,
            position_angle=mod(atand(w[1], w[2]), 180.0))
end

"Algebraic (Kåsa) circle through the points: `(center_f, center_g, radius)`."
function algebraic_circle(f::AbstractVector{<:Real}, g::AbstractVector{<:Real})
    A = [f g ones(length(f))]
    c = A \ (f .^ 2 .+ g .^ 2)
    f0, g0 = c[1] / 2, c[2] / 2
    return (center_f=f0, center_g=g0, radius=sqrt(max(c[3] + f0^2 + g0^2, 0.0)))
end
