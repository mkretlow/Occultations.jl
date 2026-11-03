# ---------------------------------------------------------------------------
# Light curves: an immutable photometric time series with its observation parameters
# ---------------------------------------------------------------------------

"""
    LightCurve(time, flux; tref, σ=nothing, exptime=median(diff(time)), λ=0.7, Δλ=0.3, name="")
    LightCurve(epochs::Vector{Epoch}, flux; tref=epochs[1], kwargs...)

Photometry of one station: `time` [s] relative to `tref` (an [`Epoch`](@ref), UTC) at
mid-exposure, `flux`, per-point uncertainties `σ` (default: a robust noise estimate from
the point-to-point scatter), the exposure time [s] and the passband (`λ`, `Δλ`, µm).
Light curves are immutable: [`normalize`](@ref), [`select`](@ref) return new ones.
"""
struct LightCurve
    name::String
    tref::Epoch
    time::Vector{Float64}
    flux::Vector{Float64}
    σ::Vector{Float64}
    optics::Optics
end

function LightCurve(time::AbstractVector{<:Real}, flux::AbstractVector{<:Real}; tref, σ=nothing,
                    exptime::Real=NaN, λ::Real=0.7, Δλ::Real=0.3, name::AbstractString="")
    length(time) == length(flux) || throw(DimensionMismatch("time and flux lengths differ"))
    p = sortperm(time)
    t = Float64.(time[p]); f = Float64.(flux[p])
    s = σ === nothing ? fill(length(f) > 3 ? _robust_noise(f) : 1.0, length(f)) :
        σ isa Real ? fill(Float64(σ), length(f)) : Float64.(σ[p])
    length(s) == length(f) || throw(DimensionMismatch("σ and flux lengths differ"))
    exp = isnan(exptime) ? (length(t) > 1 ? median(diff(t)) : 0.0) : Float64(exptime)
    length(t) > 2 && exp > 1.05 * median(diff(t)) &&
        @warn "LightCurve $name: exposure time $exp s exceeds the cycle time $(median(diff(t))) s"
    return LightCurve(String(name), _epoch(tref), t, f, s, Optics(λ, Δλ, exp))
end

function LightCurve(epochs::AbstractVector{Epoch}, flux::AbstractVector{<:Real}; tref=nothing, kw...)
    t0 = tref === nothing ? minimum(epochs) : _epoch(tref)
    return LightCurve([e - t0 for e in epochs], flux; tref=t0, kw...)
end

"Robust per-point noise from the median absolute deviation of first differences."
function _robust_noise(f::AbstractVector)
    d = diff(f)
    s = 1.4826 * median(abs.(d .- median(d))) / sqrt(2)
    return s > 0 ? s : max(std(f), eps())
end

Base.length(lc::LightCurve) = length(lc.time)
exptime(lc::LightCurve) = lc.optics.exptime
"UTC epochs of the light-curve points."
epochs(lc::LightCurve) = [lc.tref + t for t in lc.time]
"Observed time span `(first, last)` as Epochs (including half an exposure)."
timespan(lc::LightCurve) = (lc.tref + (lc.time[1] - exptime(lc) / 2), lc.tref + (lc.time[end] + exptime(lc) / 2))

Base.show(io::IO, lc::LightCurve) = print(io, "LightCurve(", lc.name, ", ", length(lc), " points)")
function Base.show(io::IO, ::MIME"text/plain", lc::LightCurve)
    println(io, "LightCurve ", lc.name)
    if length(lc) > 0
        @printf(io, "  %d points, %s + [%.3f, %.3f] s, exposure %.4f s, median σ = %.4g\n", length(lc),
                isot(lc.tref), lc.time[1], lc.time[end], exptime(lc), median(lc.σ))
    end
    @printf(io, "  band λ = %.3f ± %.3f µm", lc.optics.λ, lc.optics.Δλ / 2)
end

"""
    select(lc, tmin, tmax) -> LightCurve

Points with `tmin ≤ time ≤ tmax` (seconds from `tref`, or Epochs).
"""
function select(lc::LightCurve, tmin, tmax)
    a = tmin isa Epoch ? tmin - lc.tref : Float64(tmin)
    b = tmax isa Epoch ? tmax - lc.tref : Float64(tmax)
    k = (lc.time .>= a) .& (lc.time .<= b)
    return LightCurve(lc.name, lc.tref, lc.time[k], lc.flux[k], lc.σ[k], lc.optics)
end

"""
    with_optics(lc; λ, Δλ, exptime) -> LightCurve

Copy with changed passband or exposure time.
"""
with_optics(lc::LightCurve; λ=lc.optics.λ, Δλ=lc.optics.Δλ, exptime=lc.optics.exptime) =
    LightCurve(lc.name, lc.tref, lc.time, lc.flux, lc.σ, Optics(λ, Δλ, exptime))

"""
    normalize(lc; exclude=nothing, order=:auto, margin=2exptime) -> LightCurve

Divide by a polynomial fitted to the baseline (outside `exclude = (t1, t2)` [s from
`tref`]; by default the event found by [`detect_occultation`](@ref)). `order` is the
degree, or `:auto` for the degree 0–6 with the lowest BIC. The uncertainties are
rescaled to the scatter of the normalised baseline.
"""
function normalize(lc::LightCurve; exclude=nothing, order=:auto, margin::Real=2 * exptime(lc))
    t, f = lc.time, lc.flux
    if exclude === nothing
        d = detect_occultation(lc)
        exclude = d.snr > 5 ? (d.immersion, d.emersion) : (Inf, -Inf)
    end
    a, b = exclude
    sel = (t .< a - margin) .| (t .> b + margin)
    nsel = count(sel)
    tm = mean(t[sel]); ts = max(std(t[sel]), eps())
    x = (t .- tm) ./ ts
    poly(k) = (A = reduce(hcat, [x .^ j for j in 0:k]); A * (A[sel, :] \ f[sel]))
    if order === :auto
        best = (Inf, 0)
        for k in 0:min(6, nsel - 2)
            rss = sum(abs2, f[sel] .- poly(k)[sel])
            bic = nsel * log(max(rss / nsel, eps())) + (k + 1) * log(nsel)
            bic < best[1] && (best = (bic, k))
        end
        order = best[2]
    end
    nsel > order + 1 || throw(ArgumentError("not enough baseline points to normalise"))
    base = poly(order)
    fn = f ./ base
    s = lc.σ ./ base
    r = std(fn[sel]) / median(s[sel])
    return LightCurve(lc.name, lc.tref, t, fn, s .* r, lc.optics)
end

# ---------------------------------------------------------------------------
# readers and writer
# ---------------------------------------------------------------------------

"""
    read_lightcurve(file; tref=nothing, cols=(1, 2), timescale=:seconds, format=:text, kwargs...)

Read a light curve from a text file (`#` comments). `cols` = time, flux and optionally
flux-error columns; `timescale` of the time column: `:seconds` (from `tref`, required),
`:jd` (UTC Julian date) or `:iso`. `format=:pymovie`, `:tangra` or `:astroimagej` read the
files of these programs ([`read_pymovie`](@ref), [`read_tangra`](@ref),
[`read_astroimagej`](@ref)). Other keywords
(`exptime`, `λ`, `Δλ`, `name`) are passed to [`LightCurve`](@ref).
"""
function read_lightcurve(file::AbstractString; tref=nothing, cols=(1, 2), timescale::Symbol=:seconds,
                         skipstart::Integer=0, delim=nothing, format::Symbol=:text, name=basename(file), kw...)
    format === :pymovie && return read_pymovie(file; tref=tref, name=name, kw...)
    format === :tangra && return read_tangra(file; tref=tref, name=name, kw...)
    format === :astroimagej && return read_astroimagej(file; tref=tref, name=name, kw...)
    format === :text || throw(ArgumentError("unknown format $format (use :text, :pymovie, :tangra or :astroimagej)"))
    lines = [l for l in readlines(file)[skipstart+1:end] if !isempty(strip(l)) && !startswith(strip(l), "#")]
    rows = [delim === nothing ? split(strip(l)) : split(strip(l), delim) for l in lines]
    tcol = [r[cols[1]] for r in rows]
    flux = [parse(Float64, r[cols[2]]) for r in rows]
    σ = length(cols) >= 3 ? [parse(Float64, r[cols[3]]) for r in rows] : nothing
    if timescale === :seconds
        tref === nothing && throw(ArgumentError("`tref` is required with timescale=:seconds"))
        return LightCurve(parse.(Float64, tcol), flux; tref=tref, σ=σ, name=name, kw...)
    elseif timescale === :jd
        return LightCurve(Epoch.(parse.(Float64, tcol)), flux; tref=tref, σ=σ, name=name, kw...)
    elseif timescale === :iso
        return LightCurve(Epoch.(tcol), flux; tref=tref, σ=σ, name=name, kw...)
    end
    throw(ArgumentError("unknown timescale $timescale"))
end

function _hms_seconds(s::AbstractString)
    p = split(strip(s, ['[', ']', ' ', '"']), ':')
    length(p) == 3 || throw(ArgumentError("cannot parse time \"$s\""))
    return parse(Int, p[1]) * 3600 + parse(Int, p[2]) * 60 + parse(Float64, p[3])
end

function _unwrap_day(sod::AbstractVector{<:Real})
    out = Float64.(sod)
    for i in 2:length(out)
        while out[i] < out[i-1] - 43200
            out[i:end] .+= 86400
        end
    end
    return out
end

_csv_split(l) = strip.(split(l, ','))

"""
    read_pymovie(file; date=nothing, column=nothing, tref=nothing, kwargs...)

PyMovie CSV: the `timeInfo` column (`[hh:mm:ss.ssss]`) combined with `date`
(`"YYYY-MM-DD"`, or read from a `#` header line); crossings of 0h UT are handled.
`column` selects the flux column by name (default: the first `signal-*`).
"""
function read_pymovie(file::AbstractString; date=nothing, column=nothing, tref=nothing,
                      name::AbstractString=basename(file), kw...)
    lines = readlines(file)
    if date === nothing
        for l in lines
            startswith(strip(l), "#") || continue
            m = match(r"date.*?(\d{4}-\d{2}-\d{2})", l)
            m === nothing || (date = m.captures[1]; break)
        end
        date === nothing && throw(ArgumentError("no date found in $file: pass `date=\"YYYY-MM-DD\"`"))
    end
    ih = findfirst(l -> !startswith(strip(l), "#") && occursin("timeInfo", l), lines)
    ih === nothing && throw(ArgumentError("$file: no PyMovie header line with `timeInfo`"))
    hdr = _csv_split(lines[ih])
    it = findfirst(==("timeInfo"), hdr)
    ic = column === nothing ? findfirst(h -> startswith(h, "signal"), hdr) : findfirst(==(String(column)), hdr)
    ic === nothing && throw(ArgumentError("$file: flux column $(something(column, "signal-*")) not found"))
    rows = [_csv_split(l) for l in lines[ih+1:end] if !isempty(strip(l)) && !startswith(strip(l), "#")]
    sod = _unwrap_day([_hms_seconds(r[it]) for r in rows])
    flux = [parse(Float64, r[ic]) for r in rows]
    d0 = Epoch(string(date))
    t0 = tref === nothing ? d0 + floor(sod[1]) : _epoch(tref)
    return LightCurve(sod .- (t0 - d0), flux; tref=t0, name=name, kw...)
end

"""
    read_tangra(file; date, object=1, tref=nothing, kwargs...)

Tangra CSV export: `Time (UT)` (`hh:mm:ss.sss`) or `Time (Days UT)` column plus `date`;
`object` selects `SignalMinusBackground (object)`, or `Signal (object)` if there is no
background-subtracted column.
"""
function read_tangra(file::AbstractString; date, object::Integer=1, tref=nothing,
                     name::AbstractString=basename(file), kw...)
    lines = readlines(file)
    ih = findfirst(l -> startswith(strip(l), "FrameNo"), lines)
    ih === nothing && throw(ArgumentError("$file: no Tangra header line starting with `FrameNo`"))
    hdr = _csv_split(lines[ih])
    it = findfirst(h -> startswith(h, "Time"), hdr)
    ic = something(findfirst(h -> startswith(h, "SignalMinusBackground ($object)"), hdr),
                   findfirst(h -> startswith(h, "Signal ($object)"), hdr), Some(nothing))
    (it === nothing || ic === nothing) && throw(ArgumentError("$file: time or object-$object column not found"))
    days = occursin("Days", hdr[it])
    rows = [r for r in (_csv_split(l) for l in lines[ih+1:end] if !isempty(strip(l))) if length(r) >= max(it, ic)]
    sod = _unwrap_day([days ? 86400 * parse(Float64, r[it]) : _hms_seconds(r[it]) for r in rows])
    flux = [parse(Float64, r[ic]) for r in rows]
    d0 = Epoch(string(date))
    t0 = tref === nothing ? d0 + floor(sod[1]) : _epoch(tref)
    return LightCurve(sod .- (t0 - d0), flux; tref=t0, name=name, kw...)
end

"""
    read_astroimagej(file; target=1, column=nothing, error_column=nothing, tref=nothing, kwargs...)

AstroImageJ measurement table (`.tbl`, `.xls`, `.csv`: tab- or comma-separated, header in the
first line). Time from `JD_UTC`, or `J.D.-2400000` (UTC); the barycentric `BJD_TDB` is not
used, since contact times must be geocentric UTC. Flux and error from
`rel_flux_T<target>` / `rel_flux_err_T<target>` (relative to the comparison stars), or from
`column` / `error_column` (e.g. `"Source-Sky_T1"`, `"Source_Error_T1"`).
"""
function read_astroimagej(file::AbstractString; target::Integer=1, column=nothing, error_column=nothing,
                          tref=nothing, name::AbstractString=basename(file), kw...)
    lines = [l for l in readlines(file) if !isempty(strip(l))]
    isempty(lines) && throw(ArgumentError("$file is empty"))
    sep = occursin('\t', lines[1]) ? '\t' : ','
    hdr = strip.(split(lines[1], sep))
    col(n) = findfirst(==(n), hdr)
    jd, off = if col("JD_UTC") !== nothing
        col("JD_UTC"), 0.0
    elseif col("J.D.-2400000") !== nothing
        col("J.D.-2400000"), 2400000.0
    else
        throw(ArgumentError("$file: no `JD_UTC` or `J.D.-2400000` column (AstroImageJ measurement table)"))
    end
    fname = column === nothing ? "rel_flux_T$target" : String(column)
    ename = error_column === nothing ? (column === nothing ? "rel_flux_err_T$target" : nothing) : String(error_column)
    ic = col(fname)
    ic === nothing && throw(ArgumentError("$file: flux column `$fname` not found"))
    ie = ename === nothing ? nothing : col(ename)
    rows = [strip.(split(l, sep)) for l in lines[2:end]]
    rows = [r for r in rows if length(r) >= max(jd, ic) && tryparse(Float64, r[jd]) !== nothing]
    t = [Epoch(parse(Float64, r[jd]) + off) for r in rows]
    flux = [parse(Float64, r[ic]) for r in rows]
    σ = ie === nothing ? nothing : [parse(Float64, r[ie]) for r in rows]
    return LightCurve(t, flux; tref=tref, σ=σ, name=name, kw...)
end

"""
    write_lightcurve(file, lc; model=nothing)

Write time [s from `tref`], flux, σ (and a model column).
"""
function write_lightcurve(file::AbstractString, lc::LightCurve; model=nothing)
    open(file, "w") do io
        println(io, "# ", lc.name, "  tref = ", isot(lc.tref; digits=6), " UTC  exptime = ", exptime(lc),
                " s  λ = ", lc.optics.λ, " ± ", lc.optics.Δλ / 2, " µm")
        println(io, "# time[s] flux sigma", model === nothing ? "" : " model")
        for i in eachindex(lc.time)
            @printf(io, "%.6f %.8g %.4g", lc.time[i], lc.flux[i], lc.σ[i])
            model === nothing || @printf(io, " %.8g", model[i])
            println(io)
        end
    end
    return file
end

# ---------------------------------------------------------------------------
# detection and simulation
# ---------------------------------------------------------------------------

function _box_search(t, f, w, exptime; maxdur=Inf, mindur=0.0)
    N = length(f)
    S = [0.0; cumsum(f .* w)]
    W = [0.0; cumsum(w)]
    best = (snr=-Inf, i=1, j=1, depth=0.0, baseline=NaN)
    for i in 1:N
        for j in i:N-1
            (t[j] - t[i] > maxdur) && break
            win, wout = W[j+1] - W[i], W[end] - (W[j+1] - W[i])
            (j - i + 1 > N - 2) && break
            (t[j] - t[i] + exptime < mindur) && continue
            min_ = (S[j+1] - S[i]) / win
            base = (S[end] - (S[j+1] - S[i])) / wout
            depth = base - min_
            snr = depth / sqrt(1 / win + 1 / wout)
            snr > best.snr && (best = (snr=snr, i=i, j=j, depth=depth, baseline=base))
        end
    end
    return best
end

"""
    detect_occultations(lc; n=1, snr_limit=nothing, maxdur=Inf, mindur=0, tmin=-Inf, tmax=Inf,
                        method=:box, level=0.5, minpoints=2)
        -> Vector{NamedTuple}

Search for flux drops. Each result: `(immersion, emersion, time_err, duration, depth,
baseline, snr)` [s from `tref`].

* `method = :box` (default): box search (weighted) for up to `n` drops (e.g. a body and its
  rings); after each detection its points are replaced by the baseline. With `snr_limit`
  the search goes on while the signal-to-noise ratio exceeds it. Sensitive to faint drops,
  but two deep drops with a short gap between them are found as one box.
* `method = :threshold`: every run of at least `minpoints` points below `level` × the
  baseline (the median flux; the occulted flux is taken as 0) is one drop, e.g. the drops of
  a chord crossing two bodies or the two lobes of a contact binary. For drops deeper than a
  few times the noise (`n` and `snr_limit` limit the result to the strongest drops).
"""
function detect_occultations(lc::LightCurve; n=nothing, snr_limit=nothing, maxdur::Real=Inf, mindur::Real=0.0,
                             tmin::Real=-Inf, tmax::Real=Inf, method::Symbol=:box, level::Real=0.5,
                             minpoints::Integer=2)
    method in (:box, :threshold) || throw(ArgumentError("method must be :box or :threshold"))
    method === :threshold && return _threshold_drops(lc; n=n, snr_limit=snr_limit, tmin=tmin, tmax=tmax, level=level,
                                                     minpoints=minpoints)
    sel = (lc.time .>= tmin) .& (lc.time .<= tmax)
    t, f, w = lc.time[sel], copy(lc.flux[sel]), 1 ./ lc.σ[sel] .^ 2
    N = length(f)
    N < 4 && return NamedTuple[]
    nmax = n === nothing ? (snr_limit === nothing ? 1 : typemax(Int)) : n
    dtc = median(diff(t))
    out = NamedTuple[]
    while length(out) < nmax
        b = _box_search(t, f, w, exptime(lc); maxdur=maxdur, mindur=mindur)
        isfinite(b.snr) || break
        snr_limit !== nothing && b.snr < snr_limit && break
        i, j = b.i, b.j
        ti = i > 1 ? (t[i-1] + t[i]) / 2 : t[i] - exptime(lc) / 2
        te = j < N ? (t[j] + t[j+1]) / 2 : t[j] + exptime(lc) / 2
        push!(out, (immersion=ti, emersion=te, time_err=dtc / 2, duration=te - ti, depth=b.depth,
                    baseline=b.baseline, snr=b.snr))
        f[i:j] .= b.baseline
    end
    return sort!(out; by=d -> d.immersion)
end

function _threshold_drops(lc::LightCurve; n=nothing, snr_limit=nothing, tmin=-Inf, tmax=Inf, level=0.5, minpoints=2)
    sel = (lc.time .>= tmin) .& (lc.time .<= tmax)
    t, f, σ = lc.time[sel], lc.flux[sel], lc.σ[sel]
    N = length(f)
    N < 4 && return NamedTuple[]
    base = median(f)
    low = f .< level * base
    dtc = median(diff(t))
    out = NamedTuple[]
    i = 1
    while i <= N
        if !low[i]
            i += 1
            continue
        end
        j = i
        while j < N && low[j+1]
            j += 1
        end
        if j - i + 1 >= minpoints
            ti = i > 1 ? (t[i-1] + t[i]) / 2 : t[i] - exptime(lc) / 2
            te = j < N ? (t[j] + t[j+1]) / 2 : t[j] + exptime(lc) / 2
            depth = base - mean(f[i:j])
            snr = depth / (median(σ) / sqrt(j - i + 1))
            push!(out, (immersion=ti, emersion=te, time_err=dtc / 2, duration=te - ti, depth=depth, baseline=base, snr=snr))
        end
        i = j + 1
    end
    snr_limit === nothing || filter!(d -> d.snr >= snr_limit, out)
    if n !== nothing && length(out) > n
        out = sort(out; by=d -> -d.snr)[1:n]
    end
    return sort!(out; by=d -> d.immersion)
end

"""
    detect_occultation(lc; kwargs...) -> NamedTuple

The most significant flux drop (see [`detect_occultations`](@ref)).
"""
function detect_occultation(lc::LightCurve; kw...)
    r = detect_occultations(lc; n=1, kw...)
    isempty(r) && throw(ArgumentError("light curve too short"))
    return r[1]
end

"""
    simulate_lightcurve(time, t_imm, t_emm; tref, vel, dist, star_diam=0, λ=0.7, Δλ=0.3,
                        exptime, σ=0, flux_top=1, flux_bottom=0, vel_emersion=vel,
                        opacity=1, rng, name="sim") -> LightCurve

    simulate_lightcurve(time, drops; vel, kwargs...)

Synthetic light curve from [`occultation_flux`](@ref) with Gaussian noise `σ`; with
`drops = [(t_imm, t_emm), …]` several drops (binary systems, contact binaries).
"""
function simulate_lightcurve(time::AbstractVector{<:Real}, t_imm::Real, t_emm::Real; tref, vel::Real, dist::Real,
                             star_diam::Real=0.0, λ::Real=0.7, Δλ::Real=0.3,
                             exptime::Real=(length(time) > 1 ? median(diff(time)) : 0.0), σ::Real=0.0,
                             flux_top::Real=1.0, flux_bottom::Real=0.0, vel_emersion::Real=vel,
                             opacity::Real=1.0, rng=Random.default_rng(), name::AbstractString="sim")
    I = occultation_flux(time, t_imm, t_emm; vel=vel, dist=dist, star_diam=star_diam, λ=λ, Δλ=Δλ,
                         exptime=exptime, vel_emersion=vel_emersion, opacity=opacity)
    f = flux_bottom .+ (flux_top - flux_bottom) .* I .+ σ .* randn(rng, length(I))
    return LightCurve(time, f; tref=tref, σ=σ > 0 ? σ : nothing, exptime=exptime, λ=λ, Δλ=Δλ, name=name)
end

function simulate_lightcurve(time::AbstractVector{<:Real}, drops::AbstractVector{<:Tuple}; tref, vel, dist::Real,
                             star_diam::Real=0.0, λ::Real=0.7, Δλ::Real=0.3,
                             exptime::Real=(length(time) > 1 ? median(diff(time)) : 0.0), σ::Real=0.0,
                             flux_top::Real=1.0, flux_bottom::Real=0.0, vel_emersion=vel,
                             opacity=1.0, rng=Random.default_rng(), name::AbstractString="sim")
    I = occultation_flux(time, drops; vel=vel, dist=dist, star_diam=star_diam, λ=λ, Δλ=Δλ,
                         exptime=exptime, vel_emersion=vel_emersion, opacity=opacity)
    f = flux_bottom .+ (flux_top - flux_bottom) .* I .+ σ .* randn(rng, length(I))
    return LightCurve(time, f; tref=tref, σ=σ > 0 ? σ : nothing, exptime=exptime, λ=λ, Δλ=Δλ, name=name)
end

"""
    magnitude_drop(mag_star, mag_body) -> (Δm, bottom_flux)

Magnitude drop of an occultation and the normalised residual flux of the body.
"""
function magnitude_drop(mag_star::Real, mag_body::Real)
    contrast = 1 / (1 + 10^((mag_star - mag_body) * 0.4))
    mag_combined = mag_star - 2.5 * log10(1 / contrast)
    return mag_body - mag_combined, 10^((mag_combined - mag_body) * 0.4)
end

# weighted least-squares flux levels of the model fb + (ft - fb) I (free ones only), ft ≥ fb
function _flux_levels(I, f, w, ft, fb, fit_top::Bool, fit_bot::Bool)
    (fit_top || fit_bot) || return ft, fb
    J = 1 .- I
    sII, sJJ, sIJ = sum(w .* I .^ 2), sum(w .* J .^ 2), sum(w .* I .* J)
    sIf, sJf = sum(w .* I .* f), sum(w .* J .* f)
    if fit_top && fit_bot
        det = sII * sJJ - sIJ^2
        if det > 1e-10 * sII * sJJ
            ft, fb = (sIf * sJJ - sJf * sIJ) / det, (sJf * sII - sIf * sIJ) / det
        elseif sII > 0
            ft = (sIf - fb * sIJ) / sII
        end
    elseif fit_top
        sII > 0 && (ft = (sIf - fb * sIJ) / sII)
    else
        sJJ > 0 && (fb = (sJf - ft * sIJ) / sJJ)
    end
    ft < fb && (ft = fb = fit_top && fit_bot ? sum(w .* f) / sum(w) : fit_top ? fb : ft)
    return ft, fb
end
