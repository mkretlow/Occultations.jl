# ---------------------------------------------------------------------------
# Timing of a single light curve (the classical first step of a reduction):
# immersion/emersion instants from a diffraction model fit
# ---------------------------------------------------------------------------

"""
    TimingProblem(lc; vel, dist, star_diam=0, vel_emersion=vel, edges=:both,
                  immersion=nothing, emersion=nothing, fit_opacity=false, opacity=1,
                  fit_flux=true, flux_top=nothing, flux_bottom=nothing,
                  window=nothing, tmin=-Inf, tmax=Inf, sigma_model=0)

Fit of the contact times of one light curve with the diffraction model
([`occultation_flux`](@ref)): Fresnel diffraction, stellar disc (`star_diam`, km at the
body), passband and exposure of `lc`, shadow velocities normal to the limb `vel`,
`vel_emersion` [km/s] and distance `dist` [km].

Parameters (seconds from `lc.tref`): `immersion`, `emersion` (`edges = :both`,
`:immersion` or `:emersion` for a single edge) and optionally `opacity`. The flux levels
outside/inside the event are solved for linearly at every evaluation (profiled) when
`fit_flux` (levels given explicitly are fixed); the bottom level is held at
`flux_bottom` (default 0) for events too short to reach it or when the opacity is fitted.
Initial times come from [`detect_occultation`](@ref) unless given. Solve with
[`fit_timing`](@ref) / [`fit`](@ref), sample with [`sample`](@ref).
"""
struct TimingProblem <: AbstractProblem
    lc::LightCurve
    ps::Vector{Param}
    t::Vector{Float64}
    f::Vector{Float64}
    w::Vector{Float64}
    σ::Vector{Float64}
    vel::Float64
    vel_emersion::Float64
    dist::Float64
    star_diam::Float64
    edges::Symbol
    fit_top::Bool
    fit_bot::Bool
    flux_top::Float64
    flux_bottom::Float64
    profile::EdgeProfile
end

parameters(p::TimingProblem) = p.ps

function TimingProblem(lc::LightCurve; vel::Real, dist::Real, star_diam::Real=0.0, vel_emersion::Real=vel,
                       edges::Symbol=:both, immersion=nothing, emersion=nothing, fit_opacity::Bool=false,
                       opacity::Real=1.0, fit_flux::Bool=true, flux_top=nothing, flux_bottom=nothing,
                       window=nothing, tmin::Real=-Inf, tmax::Real=Inf, sigma_model::Real=0.0)
    edges in (:both, :immersion, :emersion) || throw(ArgumentError("edges must be :both, :immersion or :emersion"))
    if immersion === nothing || emersion === nothing
        d = detect_occultation(lc; tmin=tmin, tmax=tmax)
        immersion = something(immersion, d.immersion)
        emersion = something(emersion, d.emersion)
    end
    immersion isa Epoch && (immersion = immersion - lc.tref)
    emersion isa Epoch && (emersion = emersion - lc.tref)
    sel = (lc.time .>= tmin) .& (lc.time .<= tmax)
    if window !== nothing
        w = falses(length(lc))
        edges != :emersion && (w .|= abs.(lc.time .- immersion) .<= window)
        edges != :immersion && (w .|= abs.(lc.time .- emersion) .<= window)
        sel .&= w
    end
    t, f = lc.time[sel], lc.flux[sel]
    σ = sqrt.(lc.σ[sel] .^ 2 .+ sigma_model^2)
    sc = max(exptime(lc), 1e-3)
    ps = Param[]
    edges != :emersion && push!(ps, Param(:immersion, immersion; scale=sc, lower=t[1] - 10sc, upper=t[end] + 10sc))
    edges != :immersion && push!(ps, Param(:emersion, emersion; scale=sc, lower=t[1] - 10sc, upper=t[end] + 10sc))
    push!(ps, Param(:opacity, opacity; free=fit_opacity, lower=0.0, upper=1.0, scale=0.05))
    prof = _cached_profile(dist, star_diam, lc.optics.λ, lc.optics.Δλ)
    fit_top = fit_flux && flux_top === nothing
    fit_bot = fit_flux && flux_bottom === nothing && !fit_opacity
    ti = edges == :emersion ? -Inf : immersion
    te = edges == :immersion ? Inf : emersion
    if fit_bot
        I0 = _model_intensity(t, ti, te, vel, vel_emersion, exptime(lc), prof, dist, star_diam,
                              lc.optics.λ, lc.optics.Δλ, opacity)
        count(<(0.1), I0) >= 3 || (fit_bot = false)
    end
    p = TimingProblem(lc, ps, t, f, 1 ./ σ .^ 2, σ, abs(vel), abs(vel_emersion), dist, star_diam, edges,
                      fit_top, fit_bot, Float64(something(flux_top, 1.0)), Float64(something(flux_bottom, 0.0)), prof)
    _scan_edges!(p)
    return p
end

function _timing_model(p::TimingProblem, v)
    k = 1
    ti = p.edges == :emersion ? -Inf : v[(k += 1) - 1]
    te = p.edges == :immersion ? Inf : v[(k += 1) - 1]
    op = v[k]
    te > ti || return nothing
    I = _model_intensity(p.t, ti, te, p.vel, p.vel_emersion, exptime(p.lc), p.profile, p.dist, p.star_diam,
                         p.lc.optics.λ, p.lc.optics.Δλ, op)
    ft, fb = _flux_levels(I, p.f, p.w, p.flux_top, p.flux_bottom, p.fit_top, p.fit_bot)
    return I, ft, fb
end

function residuals(p::TimingProblem, v)
    m = _timing_model(p, v)
    m === nothing && return nothing
    I, ft, fb = m
    return @. (p.f - fb - (ft - fb) * I) / p.σ
end
nobs(p::TimingProblem) = length(p.t) - p.fit_top - p.fit_bot

function derived(p::TimingProblem, v)
    m = _timing_model(p, v)
    _, ft, fb = m === nothing ? (nothing, NaN, NaN) : m
    p.edges == :both || return (flux_top=ft, flux_bottom=fb)
    return (duration=v[2] - v[1], flux_top=ft, flux_bottom=fb)
end

"Model flux of a timing problem for the parameter values `v` (at the light-curve times)."
function model_flux(p::TimingProblem, v)
    I, ft, fb = _timing_model(p, v)
    return fb .+ (ft - fb) .* I
end

# coarse independent scans of each edge for robust starting values
function _scan_edges!(p::TimingProblem)
    sc = max(exptime(p.lc), 1e-3)
    v = [q.value for q in p.ps]
    χ(v) = (r = residuals(p, v); r === nothing ? Inf : sum(abs2, r))
    for _ in 1:2, k in 1:(p.edges == :both ? 2 : 1)
        grid = v[k] .+ sc .* range(-4, 4; length=161)
        c = [χ((u = copy(v); u[k] = g; u)) for g in grid]
        v[k] = grid[argmin(c)]
    end
    for (k, q) in enumerate(p.ps)
        q.value = v[k]
    end
    return p
end

"""
    fit_timing(lc; kwargs...) -> FitResult

Fit the contact times of a light curve ([`TimingProblem`](@ref); keywords as there plus
those of [`fit`](@ref)). Use [`contact_times`](@ref) for the instants as Epochs,
[`sample`](@ref) for the posterior (asymmetric, non-Gaussian timing errors),
[`profile`](@ref) / [`confidence_interval`](@ref) for χ² profiles.
"""
function fit_timing(lc::LightCurve; starts::Integer=3, rng=Random.default_rng(), scale_errors::Bool=false, kw...)
    return fit(TimingProblem(lc; kw...); starts=starts, rng=rng, scale_errors=scale_errors)
end

"""
    contact_times(timing_fit) -> NamedTuple

`(immersion, immersion_err, emersion, emersion_err)` as Epochs and seconds (`nothing` for
a missing edge).
"""
function contact_times(r::FitResult{TimingProblem})
    p = r.problem
    get(n) = n in r.names ? (p.lc.tref + r[n][1], r[n][2]) : (nothing, NaN)
    (ti, σi), (te, σe) = get(:immersion), get(:emersion)
    return (immersion=ti, immersion_err=σi, emersion=te, emersion_err=σe)
end

"""
    contact_times(posterior) -> NamedTuple

Contact instants (posterior medians) with asymmetric errors `(-σ, +σ)`.
"""
function contact_times(post::Posterior{TimingProblem})
    s = summarize(post)
    p = post.problem
    get(n) = haskey(s, n) ? (p.lc.tref + s[n].median, (s[n].minus, s[n].plus)) : (nothing, (NaN, NaN))
    (ti, σi), (te, σe) = get(:immersion), get(:emersion)
    return (immersion=ti, immersion_err=σi, emersion=te, emersion_err=σe)
end
