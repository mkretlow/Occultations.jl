# ---------------------------------------------------------------------------
# Parameters and priors shared by all fitting problems
# ---------------------------------------------------------------------------

"""
    Prior

Prior distribution of a parameter: [`Gaussian`](@ref)`(μ, σ)` or [`Flat`](@ref)`()` (uniform
within the parameter bounds). Gaussian priors enter the least-squares fits as additional
residuals `(p - μ)/σ` (maximum a posteriori) and the posterior sampling exactly.
"""
abstract type Prior end

"Gaussian prior with mean `μ` and standard deviation `σ`."
struct Gaussian <: Prior
    μ::Float64
    σ::Float64
end

"Uniform prior within the parameter bounds."
struct Flat <: Prior end

"""
    LogPrior(logpdf)

Prior with an arbitrary (unnormalised) log density `logpdf(x)`. It is used exactly by
[`sample`](@ref); the least-squares [`fit`](@ref) treats it as flat. See
[`isotropic_pole`](@ref).
"""
struct LogPrior{F} <: Prior
    logpdf::F
end

"""
    isotropic_pole()

Prior for the declination of a pole [deg] that makes the pole direction uniform on the sky
(density ∝ cos δ), for use with a flat prior on the right ascension.
"""
isotropic_pole() = LogPrior(δ -> log(max(cosd(δ), 1e-300)))

_prior_residual(::Flat, x) = nothing
_prior_residual(p::Gaussian, x) = (x - p.μ) / p.σ
_prior_residual(::LogPrior, x) = nothing
_prior_logpdf(::Prior, x) = 0.0
_prior_logpdf(p::LogPrior, x) = p.logpdf(x)

"""
    Param(name, value; free=true, lower=-Inf, upper=Inf, scale=1, prior=Flat(), period=0)

One model parameter: its initial (or fixed) `value`, whether it is `free`, its bounds,
a typical `scale` (step size for the starts and the numerical derivatives), its `prior`,
and a `period` (e.g. 180 for a position angle, 360 for a rotation phase; 0 = none).
"""
mutable struct Param
    name::Symbol
    value::Float64
    free::Bool
    lower::Float64
    upper::Float64
    scale::Float64
    prior::Prior
    period::Float64
end

Param(name::Symbol, value::Real; free::Bool=true, lower::Real=-Inf, upper::Real=Inf, scale::Real=1.0,
      prior::Prior=Flat(), period::Real=0.0) =
    Param(name, Float64(value), free, Float64(lower), Float64(upper), Float64(scale), prior, Float64(period))

Base.show(io::IO, p::Param) = @printf(io, "Param(%s = %.6g%s, [%.4g, %.4g]%s)", p.name, p.value,
                                      p.free ? "" : " fixed", p.lower, p.upper,
                                      p.prior isa Gaussian ? @sprintf(", prior N(%.4g, %.4g)", p.prior.μ, p.prior.σ) : "")

_inbounds(p::Param, x) = p.lower <= x <= p.upper

"""
    AbstractProblem

A fitting problem. Implementations provide

    parameters(problem) -> Vector{Param}
    residuals(problem, values::Vector{Float64}) -> Vector or nothing

(`values` holds all parameters in the order of `parameters`; the residuals are normalised
by their uncertainties so that χ² = Σ r²; `nothing` marks an inadmissible point, e.g. a
limb crossing a negative chord with hard constraints) and optionally

    nobs(problem)                 # number of data points (default: length of the residuals)
    derived(problem, values)      # NamedTuple of derived quantities (default: none)

They can then be solved with [`fit`](@ref), sampled with [`sample`](@ref) and profiled
with [`profile`](@ref).

# Examples
A straight line y = a + b x as a problem, with the abscissa x₀ = −a/b where it crosses zero
as a derived quantity:
```jldoctest line
julia> struct Line <: AbstractProblem
           ps::Vector{Param}
           x::Vector{Float64}
           y::Vector{Float64}
           σ::Float64
       end

julia> Occultations.parameters(p::Line) = p.ps

julia> Occultations.residuals(p::Line, v) = (v[1] .+ v[2] .* p.x .- p.y) ./ p.σ

julia> Occultations.derived(p::Line, v) = (x0 = -v[1] / v[2],)

julia> p = Line([Param(:a, 0.0), Param(:b, 1.0)], [0.0, 1, 2, 3], [1.1, 2.9, 5.2, 6.8], 0.1);

julia> r = fit(p);

julia> round.(r[:b]; digits=4)            # (value, σ)
(1.94, 0.0447)

julia> round.(r[:x0]; digits=4)           # derived, with propagated σ
(-0.5619, 0.0541)

julia> round(reduced_chi2(r); digits=3), freenames(r)
(4.1, [:a, :b])
```
"""
abstract type AbstractProblem end

"""
    parameters(problem) -> Vector{Param}

The parameters of a problem, free and fixed, in the order in which `values` are passed to
[`residuals`](@ref) and [`derived`](@ref). The `Param`s are those of the problem itself:
changing them (or using [`setparam!`](@ref)) changes the problem.
"""
function parameters end

"""
    residuals(problem, values) -> Vector or nothing

The normalised residuals (model − observation)/σ of a problem for the parameter values
`values` (all parameters, in the order of [`parameters`](@ref)), so that χ² = Σ r².
Returns `nothing` for an inadmissible point, e.g. a limb crossing a negative chord with
`negative = :hard`; such a point has χ² = ∞ and zero posterior probability.
"""
function residuals end

"""
    derived(problem, values) -> NamedTuple

Quantities derived from the parameter values `values` (in the order of
[`parameters`](@ref)), e.g. the minor axis, the equivalent radius or the astrometric offsets
of a limb fit. [`fit`](@ref) propagates their uncertainties linearly through the covariance;
[`sample`](@ref) evaluates them for every posterior sample. The default is an empty
NamedTuple.
"""
derived(::AbstractProblem, values) = NamedTuple()
# map equivalent parameter values to a canonical representation (e.g. swapped ellipse axes)
canonical(::AbstractProblem, values) = values
nobs(prob::AbstractProblem) = length(something(residuals(prob, [p.value for p in parameters(prob)]), Float64[]))

"""
    setparam!(problem, name; value, free, lower, upper, prior, scale)

Change the settings of a parameter of a problem (e.g. fix it, give a prior).
"""
function setparam!(prob::AbstractProblem, name::Symbol; value=nothing, free=nothing, lower=nothing,
                   upper=nothing, prior=nothing, scale=nothing)
    i = findfirst(p -> p.name == name, parameters(prob))
    i === nothing && throw(ArgumentError("no parameter $name; parameters are $(paramnames(prob))"))
    p = parameters(prob)[i]
    value === nothing || (p.value = value)
    free === nothing || (p.free = free)
    lower === nothing || (p.lower = lower)
    upper === nothing || (p.upper = upper)
    prior === nothing || (p.prior = prior)
    scale === nothing || (p.scale = scale)
    return prob
end

"Names of the parameters of a problem."
paramnames(prob::AbstractProblem) = [p.name for p in parameters(prob)]

# apply `fix`, `init`, `priors`, `bounds` keyword conventions to a parameter list
function _configure!(ps::Vector{Param}; fix=Symbol[], init=Dict(), priors=Dict(), bounds=Dict())
    names = [p.name for p in ps]
    get(i) = (k = findfirst(==(Symbol(i)), names); k === nothing &&
              throw(ArgumentError("unknown parameter $(i); parameters are $(names)")); ps[k])
    for n in fix
        get(n).free = false
    end
    for (n, v) in pairs(init)
        get(n).value = v
    end
    for (n, pr) in pairs(priors)
        p = get(n)
        p.prior = pr isa Prior ? pr : Gaussian(pr[1], pr[2])
    end
    for (n, (lo, hi)) in pairs(bounds)
        p = get(n)
        p.lower, p.upper = lo, hi
    end
    return ps
end
