# ---------------------------------------------------------------------------
# Tracks: the path of one observer in the fundamental plane, (f, g)(t), computed once
# (with the full astrometric reduction) and represented by Chebyshev polynomials, so
# that all fits work on cheap, smooth functions.
# ---------------------------------------------------------------------------

struct _ChebSeg
    lo::Float64
    hi::Float64
    cf::Vector{Float64}
    cg::Vector{Float64}
end

"""
    Track

Position (f, g) [km] of the star relative to the body (ephemeris) centre seen by one
observer, as a function of time `τ` [s] relative to the epoch `t0`, valid on
`[τmin, τmax]`. Evaluate with `track(τ) -> (f, g)` and [`velocity`](@ref)`(track, τ)`.
`dist` is the observer–body distance [km] at the middle of the interval.
"""
struct Track
    t0::Epoch
    τmin::Float64
    τmax::Float64
    segs::Vector{_ChebSeg}
    dist::Float64
    observer::AbstractObserver
end

function _cheb_fit(fun, lo, hi, deg)
    n = deg + 1
    θ = [(k - 0.5) * π / n for k in 1:n]
    x = cos.(θ)
    vals = [fun(lo + (xi + 1) * (hi - lo) / 2) for xi in x]
    F = first.(vals); G = last.(vals)
    cf = [2 / n * sum(F[k] * cos(j * θ[k]) for k in 1:n) for j in 0:deg]
    cg = [2 / n * sum(G[k] * cos(j * θ[k]) for k in 1:n) for j in 0:deg]
    cf[1] /= 2; cg[1] /= 2
    return cf, cg
end

@inline function _clenshaw(c, x)
    b1 = zero(x); b2 = zero(x)
    for j in length(c):-1:2
        b1, b2 = 2x * b1 - b2 + c[j], b1
    end
    return x * b1 - b2 + c[1]
end

"""
    Track(star, body, observer, t0, τmin, τmax; segment=600, degree=12)

Track of `observer` for `τ ∈ [τmin, τmax]` seconds from `t0`, from the full reduction
(star propagation and parallax, light-time-corrected topocentric body position, Earth
orientation), fitted by Chebyshev polynomials of `degree` on segments of at most
`segment` seconds (accuracy ≪ 1 m).
"""
function Track(star::Star, body, obs::AbstractObserver, t0::Epoch, τmin::Real, τmax::Real;
               segment::Real=600.0, degree::Integer=12)
    eph = _ephem(body)
    τmax > τmin || throw(ArgumentError("empty track interval"))
    nseg = max(1, ceil(Int, (τmax - τmin) / segment))
    edges = range(τmin, τmax; length=nseg + 1)
    fun(τ) = (r = _fg_at(star, eph, obs, t0 + τ); (r[1], r[2]))
    segs = [_ChebSeg(edges[i], edges[i+1], _cheb_fit(fun, edges[i], edges[i+1], degree)...) for i in 1:nseg]
    D = _fg_at(star, eph, obs, t0 + (τmin + τmax) / 2)[3]
    return Track(t0, Float64(τmin), Float64(τmax), segs, D, obs)
end

@inline function _seg(tr::Track, τ)
    n = length(tr.segs)
    n == 1 && return tr.segs[1]
    i = clamp(floor(Int, (τ - tr.τmin) / (tr.τmax - tr.τmin) * n) + 1, 1, n)
    return tr.segs[i]
end

function (tr::Track)(τ::Real)
    s = _seg(tr, τ)
    x = (2τ - s.lo - s.hi) / (s.hi - s.lo)
    return _clenshaw(s.cf, x), _clenshaw(s.cg, x)
end

"Velocity (df/dt, dg/dt) [km/s] of a track at `τ`."
function velocity(tr::Track, τ::Real; h::Real=1e-2)
    f1, g1 = tr(τ - h)
    f2, g2 = tr(τ + h)
    return (f2 - f1) / 2h, (g2 - g1) / 2h
end

Base.show(io::IO, t::Track) = print(io, "Track(", t.observer isa Observer ? t.observer.name : t.observer,
                                    @sprintf(", [%.1f, %.1f] s from ", t.τmin, t.τmax), isot(t.t0; digits=1), ")")
