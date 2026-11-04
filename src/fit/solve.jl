# ---------------------------------------------------------------------------
# Least-squares solution of a problem: multistart Levenberg–Marquardt, covariance,
# propagation of the uncertainties to derived quantities, profile likelihoods
# ---------------------------------------------------------------------------

# free-parameter vector <-> full value vector
struct _Free
    ps::Vector{Param}
    idx::Vector{Int}
    base::Vector{Float64}
end
_Free(ps::Vector{Param}) = _Free(ps, findall(p -> p.free, ps), [p.value for p in ps])
function _full(F::_Free, x)
    v = Vector{promote_type(eltype(F.base), eltype(x))}(F.base)     # generic: dual numbers for AD
    v[F.idx] .= x
    return v
end

# total residual vector (data + Gaussian priors) for the free parameters x; nothing if inadmissible
function _total_residuals(prob, F::_Free, x)
    for (k, i) in enumerate(F.idx)
        _inbounds(F.ps[i], x[k]) || return nothing
    end
    v = _full(F, x)
    r = residuals(prob, v)
    r === nothing && return nothing
    any(!isfinite, r) && return nothing
    pr = eltype(r)[]
    for (i, p) in enumerate(F.ps)
        q = _prior_residual(p.prior, v[i])
        q === nothing || push!(pr, q)
    end
    return isempty(pr) ? r : vcat(r, pr)
end

# Jacobian of the residuals: forward-mode automatic differentiation (ForwardDiff) when
# `ad[]` is true and the residual code accepts dual numbers; otherwise, and from then on
# (`ad[] = false`), central finite differences with steps 10⁻⁶ × scale
function _jacobian(res, x, r0, scale; ad::Ref{Bool}=Ref(false))
    if ad[]
        J = try
            ForwardDiff.jacobian(y -> (r = res(y); r === nothing ? throw(DomainError(y, "inadmissible")) : r), x)
        catch
            nothing
        end
        J !== nothing && size(J) == (length(r0), length(x)) && all(isfinite, J) && return J
        ad[] = false
    end
    J = zeros(length(r0), length(x))
    for k in eachindex(x)
        h = 1e-6 * scale[k]
        xp = copy(x); xp[k] += h
        xm = copy(x); xm[k] -= h
        rp, rm = res(xp), res(xm)
        J[:, k] = rp !== nothing && rm !== nothing && length(rp) == length(rm) == length(r0) ? (rp .- rm) ./ (2h) :
                  rp !== nothing && length(rp) == length(r0) ? (rp .- r0) ./ h :
                  rm !== nothing && length(rm) == length(r0) ? (r0 .- rm) ./ h : zeros(length(r0))
    end
    return J
end

"""
    levenberg_marquardt(res, x0, scale; maxiter=200, tol=1e-10)
        -> (x, χ², J, converged, iterations)

Minimise `sum(abs2, res(x))` from `x0` (`res` returns `nothing` for inadmissible points);
`scale` sets the finite-difference steps.
"""
function levenberg_marquardt(res, x0::AbstractVector, scale::AbstractVector; maxiter::Integer=200, tol::Real=1e-9,
                             ad::Ref{Bool}=Ref(true))
    x = collect(Float64, x0)
    r = res(x)
    r === nothing && return x, Inf, zeros(0, length(x)), false, 0
    c = sum(abs2, r)
    λ = 1e-3
    converged = false
    it = 0
    J = _jacobian(res, x, r, scale; ad=ad)
    for outer it in 1:maxiter
        g = J' * r
        H = J' * J
        D = Diagonal(max.(diag(H), 1e-12 * max(maximum(diag(H); init=0.0), 1e-300)))
        dx = -((H + λ * D) \ g)
        xn = x .+ dx
        rn = all(isfinite, xn) ? res(xn) : nothing
        cn = rn === nothing || length(rn) != length(r) ? Inf : sum(abs2, rn)
        if cn < c
            small = c - cn <= tol * max(c, 1.0) || all(abs.(dx) .<= 1e-7 .* scale)
            x, r, c = xn, rn, cn
            λ = max(λ / 5, 1e-12)
            J = _jacobian(res, x, r, scale; ad=ad)
            small && (converged = true; break)
        else
            λ *= 10
            λ > 1e12 && (converged = true; break)     # no further decrease possible
        end
    end
    return x, c, J, converged, it
end

"""
    FitResult

Result of [`fit`](@ref): best-fit parameters and their covariance, χ² and the derived
quantities with uncertainties propagated linearly through the covariance.

* `fit[name]` → `(value, σ)` for a parameter or a derived quantity
* `values(fit)`, `errors(fit)` → NamedTuples of all parameters (σ = 0 for fixed ones)
* `fit.covariance` (free parameters, order `freenames(fit)`), `fit.chi2`, `fit.dof`,
  `fit.converged`, `fit.derived`, `fit.derived_errors`
"""
struct FitResult{P<:AbstractProblem}
    problem::P
    names::Vector{Symbol}
    values::Vector{Float64}
    free::Vector{Bool}
    covariance::Matrix{Float64}
    chi2::Float64
    dof::Int
    converged::Bool
    residuals::Vector{Float64}
    derived::NamedTuple
    derived_errors::NamedTuple
    derived_covariance::Matrix{Float64}
end

"""
    freenames(fit) -> Vector{Symbol}

Names of the free parameters of a fit, in the order of the rows and columns of
`fit.covariance`.
"""
freenames(r::FitResult) = r.names[r.free]
Base.values(r::FitResult) = NamedTuple{Tuple(r.names)}(Tuple(r.values))

"""
    errors(fit) -> NamedTuple

The 1σ uncertainties of all parameters of a fit, √(diagonal of the covariance), as a
NamedTuple in the order of the parameters; fixed parameters have σ = 0. Together with
`values(fit)`; `fit[name]` returns `(value, σ)` of one parameter or derived quantity.
"""
function errors(r::FitResult)
    σ = zeros(length(r.names))
    σ[r.free] .= sqrt.(max.(diag(r.covariance), 0.0))
    return NamedTuple{Tuple(r.names)}(Tuple(σ))
end
function Base.getindex(r::FitResult, name::Symbol)
    i = findfirst(==(name), r.names)
    i !== nothing && return (r.values[i], errors(r)[i])
    haskey(r.derived, name) && return (r.derived[name], r.derived_errors[name])
    throw(KeyError(name))
end
"""
    reduced_chi2(fit) -> Float64

χ²/dof of a fit, with dof = number of residuals − number of free parameters (at least 1).
Values well above 1 mean that the scatter exceeds the uncertainties (see `σ_limb` and
`scale_errors` of [`fit_limb`](@ref)).
"""
reduced_chi2(r::FitResult) = r.chi2 / max(r.dof, 1)

"""
    correlation(fit) -> Matrix

Correlation matrix of the free parameters.
"""
function correlation(r::FitResult)
    σ = sqrt.(max.(diag(r.covariance), 1e-300))
    return r.covariance ./ (σ * σ')
end

function Base.show(io::IO, ::MIME"text/plain", r::FitResult)
    println(io, nameof(typeof(r.problem)), " fit: χ² = ", @sprintf("%.3f", r.chi2), ", dof = ", r.dof,
            @sprintf(", χ²/dof = %.3f", reduced_chi2(r)), r.converged ? "" : "  (NOT converged)")
    e = errors(r)
    for (i, n) in enumerate(r.names)
        @printf(io, "  %-22s = %14.6f %s\n", n, r.values[i], r.free[i] ? @sprintf("± %.6f", e[i]) : "(fixed)")
    end
    if !isempty(r.derived)
        println(io, "  derived:")
        for k in keys(r.derived)
            @printf(io, "  %-22s = %14.6f ± %.6f\n", k, r.derived[k], r.derived_errors[k])
        end
    end
end

_starts_rng(rng) = rng === nothing ? Random.default_rng() : rng

"""
    fit(problem; starts=8, maxiter=200, rng=Random.default_rng(), scale_errors=false, ad=true) -> FitResult

Least-squares (maximum a posteriori with Gaussian priors) solution: Levenberg–Marquardt
from the initial values and from `starts - 1` random starting points around them (within
±2 parameter `scale`s); the best solution is refined once more. The covariance is
`(JᵀJ)⁻¹` of the normalised residuals at the minimum (multiplied by χ²/dof if
`scale_errors`). If the data do not constrain some direction of the parameters (singular
`J`), the fit is not `converged` and the parameters and derived quantities along that
direction get σ = `Inf` (their covariances `NaN`). The uncertainties of the derived quantities follow by linear
propagation. See [`sample`](@ref) for the full posterior.

The Jacobian is computed by forward-mode automatic differentiation (ForwardDiff) if the
problem supports it, otherwise (or with `ad = false`) by central finite differences.
"""
function fit(prob::AbstractProblem; starts::Integer=8, maxiter::Integer=200, rng=Random.default_rng(),
             scale_errors::Bool=false, ad::Bool=true)
    ps = parameters(prob)
    F = _Free(ps)
    nfree = length(F.idx)
    nfree == 0 && throw(ArgumentError("no free parameters"))
    scale = [ps[i].scale for i in F.idx]
    res = x -> _total_residuals(prob, F, x)
    x0 = F.base[F.idx]
    best = (x=x0, c=Inf, J=zeros(0, nfree), conv=false)
    useAD = Ref(ad)                         # switched off for good if the problem does not support AD
    for s in 1:max(starts, 1)
        xs = s == 1 ? x0 : x0 .+ 2 .* scale .* (2 .* rand(rng, nfree) .- 1)
        xs = [clamp(xs[k], ps[i].lower, ps[i].upper) for (k, i) in enumerate(F.idx)]
        r0 = res(xs)
        if r0 === nothing                        # inadmissible start: random search nearby
            for k in 1:200
                xt = xs .+ scale .* (1 + k / 50) .* randn(rng, nfree)
                res(xt) === nothing || (xs = xt; break)
            end
            res(xs) === nothing && continue
        end
        x, c, J, conv, _ = levenberg_marquardt(res, xs, scale; maxiter=maxiter, ad=useAD)
        c < best.c && (best = (x=x, c=c, J=J, conv=conv))
    end
    isfinite(best.c) || error("fit: no admissible solution found (check the initial values and constraints)")
    xb = canonical(prob, _full(F, best.x))[F.idx]
    res(xb) === nothing && (xb = best.x)
    x, c, J, conv, _ = levenberg_marquardt(res, xb, scale; maxiter=maxiter, ad=useAD)
    # keep the refit from the canonical form unless it is worse (beyond rounding)
    c <= best.c + 1e-8 * max(best.c, 1.0) || ((x, c, J, conv) = (best.x, best.c, best.J, best.conv))
    # wrap periodic parameters
    for (k, i) in enumerate(F.idx)
        P = ps[i].period
        P > 0 && (x[k] = mod(x[k], P))
    end
    r = res(x)
    J = _jacobian(res, x, r, scale; ad=useAD)
    C, N = _covariance(J)
    ok = size(N, 2) == 0
    n = nobs(prob)
    dof = n - nfree
    scale_errors && dof > 0 && (C .*= c / dof)
    v = _full(F, x)
    d, de, dC = _propagate(prob, F, x, C, scale; null=N)
    _unconstrained!(C, N)
    return FitResult(prob, [p.name for p in ps], v, [p.free for p in ps], C, c, dof, conv && ok,
                     residuals(prob, v), d, de, dC)
end

# covariance (JᵀJ)⁻¹ from the singular value decomposition of J (not of JᵀJ, which squares its
# condition number), and an orthonormal basis N of the directions the data do not constrain
# (singular values ≤ √(n eps) of the largest; empty if J has full rank)
function _covariance(J)
    m = size(J, 2)
    size(J, 1) == 0 && return zeros(m, m), Matrix{Float64}(I, m, m)
    S = svd(J; full=true)
    k = count(>(sqrt(m * eps()) * S.S[1]), S.S)
    Vk = S.V[:, 1:k]
    return Vk * Diagonal(1 ./ S.S[1:k] .^ 2) * Vk', S.V[:, k+1:end]
end

# unconstrained parameters (a component along the null directions N): infinite variance,
# undefined covariances
function _unconstrained!(C, N)
    size(N, 2) == 0 && return C
    for k in axes(C, 1)
        norm(N[k, :]) > 1e-6 || continue
        C[k, :] .= NaN
        C[:, k] .= NaN
        C[k, k] = Inf
    end
    return C
end

# derived quantities and their linearly propagated uncertainties (infinite for a quantity
# that changes along a null direction of the fit)
function _propagate(prob, F::_Free, x, C, scale; null=zeros(length(x), 0))
    d = derived(prob, _full(F, x))
    isempty(d) && return d, d, zeros(0, 0)
    f(x) = collect(Float64, values(derived(prob, _full(F, x))))
    d0 = f(x)
    G = zeros(length(d0), length(x))
    for k in eachindex(x)
        h = 1e-5 * scale[k]
        xp = copy(x); xp[k] += h
        xm = copy(x); xm[k] -= h
        G[:, k] = (f(xp) .- f(xm)) ./ (2h)
    end
    DC = G * C * G'
    GN = G * null
    for j in axes(G, 1)
        norm(GN[j, :]) > 1e-6 * norm(G[j, :]) || continue
        DC[j, :] .= NaN
        DC[:, j] .= NaN
        DC[j, j] = Inf
    end
    σ = sqrt.(max.(diag(DC), 0.0))
    return d, NamedTuple{keys(d)}(Tuple(σ)), DC
end

"""
    profile(fit, name; values=nothing, n=21, width=4) -> (values, chi2)

Profile χ² of parameter `name`: for each value (default `n` values within ±`width` σ of the
best fit) the parameter is fixed and all others are refitted. The interval where
`Δχ² ≤ 1` is the 1σ confidence interval (asymmetric errors are revealed this way).
Values outside the parameter's bounds, and inadmissible values of a problem without other
free parameters, give `Inf`; a refit that fails gives `NaN`.
"""
function profile(r::FitResult, name::Symbol; values=nothing, n::Integer=21, width::Real=4.0)
    prob = r.problem
    ps = parameters(prob)
    i = findfirst(p -> p.name == name, ps)
    i === nothing && throw(ArgumentError("no parameter $name"))
    v0, σ = r[name]
    vals = values === nothing ? collect(range(v0 - width * σ, v0 + width * σ; length=n)) : collect(Float64, values)
    saved = [(p.value, p.free) for p in ps]
    out = similar(vals)
    try
        for (k, v) in enumerate(vals)
            for (j, p) in enumerate(ps)
                p.value = r.values[j]
            end
            ps[i].value, ps[i].free = v, false
            out[k] = if !_inbounds(ps[i], v)
                Inf                                         # excluded by the bounds
            elseif !any(p -> p.free, ps)                    # nothing to refit: χ² at this value
                rr = _total_residuals(prob, _Free(ps), Float64[])
                rr === nothing ? Inf : sum(abs2, rr)
            else
                try
                    fit(prob; starts=1).chi2
                catch
                    NaN
                end
            end
        end
    finally
        for (j, p) in enumerate(ps)
            p.value, p.free = saved[j]
        end
    end
    return vals, out
end

"""
    confidence_interval(fit, name; level=1) -> (lower, upper)

`level`-σ confidence interval of a parameter from its χ² profile (`Δχ² ≤ level²`),
found by bisection on each side (`(value, value)` for a fixed parameter). A side is `±Inf` if Δχ² stays below `level²` out to
`maxwidth × level` σ (unbounded there), `NaN` if a refit of the profile fails, and at most
the parameter's bound. The search needs a finite σ of the fit.
"""
function confidence_interval(r::FitResult, name::Symbol; level::Real=1, maxwidth::Real=10)
    v0, σ = r[name]
    i = findfirst(==(name), r.names)
    i !== nothing && !r.free[i] && return (v0, v0)                # a fixed parameter
    isfinite(σ) && σ > 0 || throw(ArgumentError("confidence_interval: σ($name) = $σ; the fit does not constrain it"))
    target = r.chi2 + level^2
    χ(v) = profile(r, name; values=[v])[2][1]
    function edge(dir)
        lo, hi = 0.0, level * σ
        c = χ(v0 + dir * hi)
        while c < target && hi < maxwidth * level * σ
            lo, hi = hi, 2hi
            c = χ(v0 + dir * hi)
        end
        isnan(c) && return NaN
        c < target && return dir * Inf                 # not bracketed: unbounded within maxwidth
        for _ in 1:25
            m = (lo + hi) / 2
            cm = χ(v0 + dir * m)
            isnan(cm) && return NaN
            cm < target ? (lo = m) : (hi = m)
        end
        return v0 + dir * (lo + hi) / 2
    end
    return edge(-1), edge(1)
end
