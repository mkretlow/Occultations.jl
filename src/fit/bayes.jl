# ---------------------------------------------------------------------------
# Bayesian inference: affine-invariant ensemble sampler (Goodman & Weare 2010, as in
# emcee) on the posterior  log p = -χ²/2 + log prior
# ---------------------------------------------------------------------------

"""
    Posterior

Result of [`sample`](@ref): Markov chains of the free parameters.

* `chain` – `nsteps × nwalkers × nfree` array, `logp` – `nsteps × nwalkers`
* `names` – free-parameter names, `burn` – discarded initial steps, `acceptance`
* `samples(post)` – flattened samples after burn-in (`nsamples × nfree`)
* `post[name]` – samples of a parameter; `summarize(post)` – median and 68.3 % intervals
* [`derived_samples`](@ref) – samples of the derived quantities (radius, axes,
  astrometry, ...); [`autocorrelation_time`](@ref) – convergence diagnostics
"""
struct Posterior{P<:AbstractProblem}
    problem::P
    names::Vector{Symbol}
    chain::Array{Float64,3}
    logp::Matrix{Float64}
    burn::Int
    acceptance::Float64
    base::Vector{Float64}
    idx::Vector{Int}
end

"Flattened samples after burn-in, `nsamples × nfree` (thinned by `thin`)."
function samples(p::Posterior; thin::Integer=1)
    c = p.chain[p.burn+1:thin:end, :, :]
    return reshape(c, :, size(c, 3))
end

function Base.getindex(p::Posterior, name::Symbol)
    k = findfirst(==(name), p.names)
    k === nothing && throw(KeyError(name))
    return samples(p)[:, k]
end

"""
    summarize(post; quantiles=(0.16, 0.5, 0.84)) -> NamedTuple of (median, minus, plus)

Median and asymmetric 68.3 % credible interval of each free parameter.
"""
function summarize(p::Posterior; quantiles=(0.15865525393145707, 0.5, 0.8413447460685429))
    S = samples(p)
    vals = map(eachindex(p.names)) do k
        x = S[:, k]
        _periodic_center!(x, p, k)
        lo, md, hi = Statistics.quantile(x, collect(quantiles))
        (median=md, minus=md - lo, plus=hi - md)
    end
    return NamedTuple{Tuple(p.names)}(Tuple(vals))
end

# centre the samples of a periodic parameter on their circular mean
function _periodic_center!(x, p::Posterior, k)
    P = parameters(p.problem)[p.idx[k]].period
    P > 0 || return x
    m = atan(mean(sin.(2π .* x ./ P)), mean(cos.(2π .* x ./ P))) * P / 2π
    x .= mod.(x .- m .+ P / 2, P) .- P / 2 .+ mod(m, P)
    return x
end

function Base.show(io::IO, ::MIME"text/plain", p::Posterior)
    ns, nw, _ = size(p.chain)
    println(io, "Posterior of ", nameof(typeof(p.problem)), ": ", nw, " walkers × ", ns, " steps (burn-in ",
            p.burn, "), acceptance ", @sprintf("%.2f", p.acceptance))
    s = summarize(p)
    for n in p.names
        v = s[n]
        @printf(io, "  %-22s = %14.6f  -%.6f +%.6f\n", n, v.median, v.minus, v.plus)
    end
end

function _logpost(prob, F::_Free, x)
    r = _total_residuals(prob, F, x)
    r === nothing && return -Inf
    lp = -sum(abs2, r) / 2
    for (k, i) in enumerate(F.idx)
        lp += _prior_logpdf(F.ps[i].prior, x[k])
    end
    return lp
end

"""
    sample(problem_or_fit; nwalkers=4nfree (≥ 16), nsteps=2000, burn=nsteps÷3, a=2,
           init=nothing, rng=Random.default_rng(), progress=false) -> Posterior

Sample the posterior distribution `p(θ | data) ∝ exp(-χ²/2) · prior` of the free
parameters with the affine-invariant ensemble ("stretch move") sampler. Flat priors
within the parameter bounds unless Gaussian priors are set ([`setparam!`](@ref),
`priors=` of the problem constructors); hard constraints (e.g. negative chords with
`negative=:hard`) have zero probability. By default (`init = nothing`) the walkers start in a
small Gaussian ball around the least-squares solution (from `fit`, computed if a problem is
given). A ball explores one mode only; to search for others, `init` gives the starts: a
`nwalkers × nfree` matrix (columns in the order of `freenames(fit)`), or a function
`rng -> x` that draws one start (drawn again where inadmissible), e.g. from broad ranges.
"""
function sample(r::FitResult; nwalkers::Union{Nothing,Integer}=nothing, nsteps::Integer=2000,
                burn::Integer=nsteps ÷ 3, a::Real=2.0, init=nothing, rng=Random.default_rng(),
                progress::Bool=false)
    prob = r.problem
    ps = parameters(prob)
    saved = [p.value for p in ps]
    for (i, p) in enumerate(ps)
        p.value = r.values[i]
    end
    F = _Free(ps)
    for (i, p) in enumerate(ps)
        p.value = saved[i]
    end
    nfree = length(F.idx)
    nw = init isa AbstractMatrix ? size(init, 1) : nwalkers === nothing ? max(16, 4nfree) : nwalkers
    init isa AbstractMatrix && (size(init, 2) == nfree && iseven(nw) ||
        throw(ArgumentError("init: need an even number of rows and $nfree columns (one per free parameter)")))
    iseven(nw) || (nw += 1)
    nw >= 2nfree || throw(ArgumentError("need at least 2 × (number of free parameters) walkers"))
    x0 = r.values[F.idx]
    σ = sqrt.(max.(diag(r.covariance), 0.0))
    scale = [ps[i].scale for i in F.idx]
    σ = [isfinite(s) && s > 0 ? s : 1e-3 * scale[k] for (k, s) in enumerate(σ)]
    X = zeros(nw, nfree)
    L = fill(-Inf, nw)
    for w in 1:nw
        for k in 1:(init isa AbstractMatrix ? 1 : 1000)
            x = init === nothing ? x0 .+ 0.1 .* σ .* randn(rng, nfree) :
                init isa AbstractMatrix ? collect(Float64, init[w, :]) : collect(Float64, init(rng))
            l = _logpost(prob, F, x)
            if isfinite(l)
                X[w, :] .= x; L[w] = l
                break
            end
        end
        isfinite(L[w]) || error("sample: cannot initialise walker $w (no admissible start",
                                init === nothing ? " near the fit)" : ")")
    end
    chain = zeros(nsteps, nw, nfree)
    logp = zeros(nsteps, nw)
    acc = 0
    half = nw ÷ 2
    for s in 1:nsteps
        for (S, C) in ((1:half, half+1:nw), (half+1:nw, 1:half))
            for w in S
                j = rand(rng, C)
                z = ((a - 1) * rand(rng) + 1)^2 / a
                y = X[j, :] .+ z .* (X[w, :] .- X[j, :])
                l = _logpost(prob, F, y)
                if log(rand(rng)) < (nfree - 1) * log(z) + l - L[w]
                    X[w, :] .= y; L[w] = l
                    acc += 1
                end
            end
        end
        chain[s, :, :] .= X
        logp[s, :] .= L
        progress && s % max(1, nsteps ÷ 10) == 0 && @info "sample: step $s / $nsteps"
    end
    # canonical representation of every sample (e.g. swapped ellipse axes)
    for s in 1:nsteps, w in 1:nw
        chain[s, w, :] .= canonical(prob, _full(F, chain[s, w, :]))[F.idx]
    end
    return Posterior(prob, [ps[i].name for i in F.idx], chain, logp, min(burn, nsteps - 1),
                     acc / (nsteps * nw), F.base, F.idx)
end

sample(prob::AbstractProblem; rng=Random.default_rng(), starts::Integer=8, kw...) =
    sample(fit(prob; rng=rng, starts=starts); rng=rng, kw...)

"""
    derived_samples(post; n=2000) -> NamedTuple of vectors

Samples of the derived quantities of the problem (e.g. equivalent radius, axes,
astrometric position) computed for `n` posterior samples.
"""
function derived_samples(p::Posterior; n::Integer=2000, rng=Random.default_rng())
    S = samples(p)
    rows = size(S, 1) <= n ? (1:size(S, 1)) : rand(rng, 1:size(S, 1), n)
    vals = map(rows) do i
        v = copy(p.base)
        v[p.idx] .= S[i, :]
        derived(p.problem, v)
    end
    (isempty(vals) || isempty(vals[1])) && return NamedTuple()
    ks = keys(vals[1])
    return NamedTuple{ks}(Tuple([Float64(v[k]) for v in vals] for k in ks))
end

"""
    derived_summary(post; n=2000) -> NamedTuple of (median, minus, plus)
"""
function derived_summary(p::Posterior; n::Integer=2000, rng=Random.default_rng())
    d = derived_samples(p; n=n, rng=rng)
    q = (0.15865525393145707, 0.5, 0.8413447460685429)
    return map(x -> (lo = Statistics.quantile(x, q[1]); md = Statistics.quantile(x, q[2]);
                     hi = Statistics.quantile(x, q[3]); (median=md, minus=md - lo, plus=hi - md)), d)
end

"""
    autocorrelation_time(post) -> Vector

Integrated autocorrelation time [steps] of each parameter (averaged over the walkers,
Sokal's automatic window). The chain should be ≳ 50 times longer.
"""
function autocorrelation_time(p::Posterior; c::Real=5)
    c_ = p.chain[p.burn+1:end, :, :]
    n = size(c_, 1)
    return map(1:size(c_, 3)) do k
        ρ = zeros(n)
        for w in 1:size(c_, 2)
            x = c_[:, w, k] .- mean(c_[:, w, k])
            v = sum(abs2, x)
            v > 0 || continue
            for lag in 0:n-1
                ρ[lag+1] += sum(x[1:n-lag] .* x[1+lag:n]) / v
            end
        end
        ρ ./= max(ρ[1], eps())
        τ = 1.0
        for m in 1:n-1
            τ += 2ρ[m+1]
            m >= c * τ && break
        end
        τ
    end
end
