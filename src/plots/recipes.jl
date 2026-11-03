# ---------------------------------------------------------------------------
# Plot recipes for Plots.jl (RecipesBase): `using Plots; plot(x)`
# ---------------------------------------------------------------------------

RecipesBase.@recipe function f(lc::LightCurve; model=nothing)
    left_margin --> (5, :mm)                      # room for the axis labels (an outer legend narrows the plot)
    bottom_margin --> (3, :mm)
    framestyle --> :box
    xguide --> "time − $(isot(lc.tref; digits=2)) UTC [s]"
    yguide --> "flux"
    title --> lc.name
    legend --> :bottomleft
    @RecipesBase.series begin
        seriestype := :scatter
        markersize --> 2
        markerstrokewidth --> 0
        markercolor --> :black
        label --> "data"
        lc.time, lc.flux
    end
    if model !== nothing
        @RecipesBase.series begin
            seriestype := :path
            linewidth --> 1.5
            linecolor --> :red
            label --> "model"
            lc.time, model
        end
    end
end

RecipesBase.@recipe function f(r::FitResult{TimingProblem}; residuals=true)
    left_margin --> (5, :mm)                      # room for the axis labels (an outer legend narrows the plot)
    bottom_margin --> (3, :mm)
    framestyle --> :box
    p = r.problem
    lc = p.lc
    m = model_flux(p, r.values)
    xguide --> "time − $(isot(lc.tref; digits=2)) UTC [s]"
    yguide --> "flux"
    title --> lc.name
    legend --> :bottomleft
    @RecipesBase.series begin
        seriestype := :scatter
        markersize --> 2
        markerstrokewidth --> 0
        markercolor --> :black
        label --> "data"
        p.t, p.f
    end
    @RecipesBase.series begin
        seriestype := :path
        linewidth --> 1.5
        linecolor --> :red
        label --> "model"
        p.t, m
    end
    if residuals
        off = minimum(p.f) - 4 * median(p.σ)
        @RecipesBase.series begin
            seriestype := :scatter
            markersize --> 1.5
            markerstrokewidth --> 0
            markercolor --> :gray50
            label --> "residuals (offset)"
            p.t, p.f .- m .+ off
        end
    end
    ts = [r[n][1] for n in (:immersion, :emersion) if n in r.names]
    @RecipesBase.series begin
        seriestype := :vline
        linestyle --> :dot
        label --> ""
        ts
    end
end

# chords of an event in the fundamental plane
# the plotted area: the axis limits in effect (given by the user or set by the recipe), or `box`
function _plot_box(attr, box)
    box === nothing && return nothing
    xl = get(attr, :xlims, nothing); yl = get(attr, :ylims, nothing)
    x = xl isa Tuple && length(xl) == 2 && all(isfinite, xl) ? xl : (box[1], box[2])
    y = yl isa Tuple && length(yl) == 2 && all(isfinite, yl) ? yl : (box[3], box[4])
    return (min(x...), max(x...), min(y...), max(y...))
end

# a compass in the lower right corner of the plotted area `box` (f grows to the left, as on the
# sky): arrows towards east (+f) and north (+g), labelled E and N
function _compass_series!(series, box)
    w, h = box[2] - box[1], box[4] - box[3]
    L = 0.09 * min(w, h)
    f0, g0 = box[1] + 0.06w, box[3] + 0.06h
    # shafts and heads as plain lines (GR's arrow heads ignore the flipped f axis)
    k, b = 0.28L, 0.14L
    x = [f0, f0 + L, NaN, f0 + L - k, f0 + L, f0 + L - k, NaN, f0, f0, NaN, f0 - b, f0, f0 + b]
    y = [g0, g0, NaN, g0 - b, g0, g0 + b, NaN, g0, g0 + L, NaN, g0 + L - k, g0 + L, g0 + L - k]
    push!(series, (x, y, (linecolor=:black, linewidth=1.2, label="")))
    push!(series, ([f0 + 1.35L, f0], [g0, g0 + 1.35L],
                   (seriestype=:scatter, markersize=0, markeralpha=0, label="",
                    series_annotations=[("E", 9, :center, :black), ("N", 9, :center, :black)])))
    return series
end

function _chord_series!(series, ev::Event; dts=Dict{String,Float64}(), numbers::Bool=true, ref=(0.0, 0.0),
                        box=nothing, extend::Bool=false, compass=nothing)
    nx = Float64[]; ny = Float64[]; nt = String[]
    for (i, c) in enumerate(ev.chords)
        dt = get(dts, c.name, 0.0)
        tag = numbers ? "$i  " : ""
        if isnegative(c)
            f, g = chord_path(ev, c; dt=dt)
            push!(series, (f, g, (linestyle=:dash, linecolor=:gray, label=tag * c.name * " (−)")))
            k = argmin(hypot.(f .- ref[1], g .- ref[2]))          # label next to the closest point
            push!(nx, f[k]); push!(ny, g[k]); push!(nt, string(i))
        else
            f0, g0 = chord_path(ev, c; dt=dt)
            extend && push!(series, (f0, g0, (linestyle=:dot, linecolor=:gray, label="")))   # the track beyond the limb
            f, g = chord_path(ev, c; segment=:occulted, dt=dt)
            if isempty(f)                                        # only ring mid-times: no occulted segment
                push!(series, (f0, g0, (linestyle=:dot, linecolor=:gray, label=tag * c.name)))
                k = argmin(hypot.(f0 .- ref[1], g0 .- ref[2]))
                push!(nx, f0[k]); push!(ny, g0[k]); push!(nt, string(i))
                continue
            end
            push!(series, (f, g, (linewidth=2, label=tag * c.name)))
            # label just outside the first contact, along the chord
            L = hypot(f[end] - f[1], g[end] - g[1])
            u = L > 0 ? ((f[1] - f[end]) / L, (g[1] - g[end]) / L) : (0.0, 0.0)
            δ = max(0.08L, 8.0)
            push!(nx, f[1] + δ * u[1]); push!(ny, g[1] + δ * u[2]); push!(nt, string(i))
        end
    end
    if box !== nothing                                   # no labels outside the plotted area
        k = findall(j -> box[1] <= nx[j] <= box[2] && box[3] <= ny[j] <= box[4], eachindex(nx))
        nx, ny, nt = nx[k], ny[k], nt[k]
    end
    if numbers && !isempty(nt)
        push!(series, (nx, ny, (seriestype=:scatter, markersize=0, markeralpha=0, label="",
                                series_annotations=[(t, 9, :center, :black) for t in nt])))
    end
    # ±1σ timing uncertainty of every contact, as a red segment along the chord
    ex = Float64[]; ey = Float64[]
    for c in ev.chords, x in c.contacts
        x.σ > 0 || continue
        tr = track(ev, c)
        τ = x.time - ev.tref + get(dts, c.name, 0.0)
        a, b = tr(τ - x.σ), tr(τ + x.σ)
        append!(ex, (a[1], b[1], NaN)); append!(ey, (a[2], b[2], NaN))
    end
    isempty(ex) || push!(series, (ex, ey, (linecolor=:red, linewidth=3, label="contact-time uncertainty")))
    compass === nothing || _compass_series!(series, compass)
    return series
end

RecipesBase.@recipe function f(ev::Event; limb=nothing, center=(0.0, 0.0), numbers=true, extend=false, compass=true,
                              dimensions=true)
    left_margin --> (5, :mm)                      # room for the axis labels (an outer legend narrows the plot)
    bottom_margin --> (3, :mm)
    aspect_ratio --> :equal
    framestyle --> :box
    xguide --> "f (km)"
    yguide --> "g (km)"
    xflip --> true
    title --> ev.name
    legend --> :outertopright
    # zoom on the contacts (and the limb): negative chords can extend far beyond the body
    pts = contact_points(ev; all=true)
    box = nothing; ref = center
    if !isempty(pts)
        fs = [p.f for p in pts]; gs = [p.g for p in pts]
        if limb !== nothing
            x, y = outline_points(limb; center=center)
            append!(fs, filter(isfinite, x)); append!(gs, filter(isfinite, y))
        end
        fc, gc = (minimum(fs) + maximum(fs)) / 2, (minimum(gs) + maximum(gs)) / 2
        h = 0.65 * max(maximum(fs) - minimum(fs), maximum(gs) - minimum(gs), 10.0)
        box = (fc - h, fc + h, gc - h, gc + h)
        ref = limb === nothing ? (fc, gc) : center
        xlims --> (box[1], box[2])
        ylims --> (box[3], box[4])
    end
    for (x, y, kw) in _chord_series!([], ev; numbers=numbers, ref=ref, box=box, extend=extend,
                                     compass=compass ? _plot_box(plotattributes, box) : nothing)
        @RecipesBase.series begin
            seriestype := :path
            for (k, v) in pairs(kw)
                plotattributes[k] = v
            end
            x, y
        end
    end
    if !isempty(pts)
        @RecipesBase.series begin
            seriestype := :scatter
            markersize --> 3
            markercolor --> :black
            label := ""
            [p.f for p in pts], [p.g for p in pts]
        end
    end
    if limb !== nothing
        @RecipesBase.series begin
            seriestype := :path
            linecolor --> :red
            label --> "limb"
            outline_points(limb; center=center)
        end
    end
    # the dimensions of the given limb, in the right margin below the legend
    if limb !== nothing && dimensions && box !== nothing
        dl = _dims_lines(limb)
        @RecipesBase.series begin
            seriestype := :scatter
            markersize := 0
            markeralpha := 0
            markerstrokewidth := 0
            label := ""
            series_annotations := [(l, 9, :left, :black) for l in dl]
            _dims_xy(box, length(dl))
        end
    end
end

# 68.3 % band of a set of limbs (vectors x, y) around c0, as a closed polygon (outer, inner):
# * pointwise: in every direction, the 16–84 % quantiles of the distances of the limbs;
# * simultaneous: the pointwise band widened by the factor k for which 68.3 % of the limbs lie
#   entirely inside (k: 68.3 % quantile over the limbs of the largest deviation from the median,
#   in units of the pointwise half-widths)
function _limb_band(draws, c0; n=361, simultaneous::Bool=false, nsigma::Real=1)
    θ = range(0, 2π; length=n)
    ρ = reduce(hcat, [_polar_radius(x .- c0[1], y .- c0[2], θ) for (x, y) in draws])
    pl = (1 - SpecialFunctions.erf(nsigma / sqrt(2))) / 2           # 0.1587 for 1σ, 0.00135 for 3σ
    md = [quantile(ρ[i, :], 0.5) for i in eachindex(θ)]
    lo = [quantile(ρ[i, :], simultaneous ? 0.1587 : pl) for i in eachindex(θ)]
    hi = [quantile(ρ[i, :], simultaneous ? 0.8413 : 1 - pl) for i in eachindex(θ)]
    if simultaneous
        dl, du = max.(md .- lo, 1e-9), max.(hi .- md, 1e-9)
        T = [maximum(i -> ρ[i, j] >= md[i] ? (ρ[i, j] - md[i]) / du[i] : (md[i] - ρ[i, j]) / dl[i], eachindex(θ))
             for j in axes(ρ, 2)]
        k = quantile(T, 1 - 2pl)
        # never beyond the drawn limbs (e.g. across a hard negative chord); coverage unchanged
        ρmin, ρmax = vec(minimum(ρ; dims=2)), vec(maximum(ρ; dims=2))
        lo, hi = max.(md .- k .* dl, ρmin), min.(md .+ k .* du, ρmax)
    end
    return [c0[1] .+ hi .* sin.(θ); c0[1] .+ reverse(lo .* sin.(θ))],
           [c0[2] .+ hi .* cos.(θ); c0[2] .+ reverse(lo .* cos.(θ))]
end

# square plot box (f1, f2, g1, g2) around a limb at centre c: ±1.6 R, or for a composite outline
# (several bodies or lobes) 1.6 × the half-extent of the union around its middle
function _limb_box(c, O)
    if O isa CompositeOutline
        x, y = outline_points(O; center=c, n=181)
        x, y = filter(isfinite, x), filter(isfinite, y)
        fc, gc = (minimum(x) + maximum(x)) / 2, (minimum(y) + maximum(y)) / 2
        h = 0.8 * max(maximum(x) - minimum(x), maximum(y) - minimum(y))
        return (fc - h, fc + h, gc - h, gc + h)
    end
    R = _rmax(O)
    return (c[1] - 1.6R, c[1] + 1.6R, c[2] - 1.6R, c[2] + 1.6R)
end

# the plot box enlarged to the rings (outer edges) around c
function _ring_box(box, c, rings)
    isempty(rings) && return box
    xs, ys = Float64[], Float64[]
    for rg in rings
        x, y = outline_points(_ring_edges(rg)[1]; center=c, n=181)
        append!(xs, x); append!(ys, y)
    end
    f1, f2 = min(box[1], minimum(xs)), max(box[2], maximum(xs))
    g1, g2 = min(box[3], minimum(ys)), max(box[4], maximum(ys))
    fc, gc, h = (f1 + f2) / 2, (g1 + g2) / 2, 0.55 * max(f2 - f1, g2 - g1)
    return (fc - h, fc + h, gc - h, gc + h)
end

# positions of the text lines of the dimensions, in the right margin of the box (x is flipped)
function _dims_xy(box, n)
    w = box[2] - box[1]
    return fill(box[1] - 0.0375w, n), [box[3] + 0.047w + 0.0406w * (n - k) for k in 1:n]
end

# band of limbs that are not star-shaped (composite outlines): the occupancy of a grid, i.e. the
# fraction of the drawn limbs (centre, outline) covering each point; the band is where it lies
# between the nσ quantiles (e.g. 15.9 % and 84.1 %). Returns the band as polygons (x, y separated
# by NaN): one rectangle per run of band cells in a row. (A heatmap would be simpler, but GR draws
# heatmaps without the flip of the f axis.)
function _occupancy_band(limbs, box; n::Integer=400, nsigma::Real=1)
    fs = range(box[1], box[2]; length=n)
    gs = range(box[3], box[4]; length=n)
    cnt = zeros(n, n)
    for (c, O) in limbs, (i, f) in enumerate(fs), (j, g) in enumerate(gs)
        signed_distance(O, f - c[1], g - c[2]) < 0 && (cnt[j, i] += 1)
    end
    pl = (1 - SpecialFunctions.erf(nsigma / sqrt(2))) / 2
    inband = [pl <= q / length(limbs) <= 1 - pl for q in cnt]
    df, dg = step(fs) / 2, step(gs) / 2
    xs, ys = Float64[], Float64[]
    for j in 1:n
        i = 1
        while i <= n
            if inband[j, i]
                k = i
                while k < n && inband[j, k+1]
                    k += 1
                end
                f1, f2, g1, g2 = fs[i] - df, fs[k] + df, gs[j] - dg, gs[j] + dg
                append!(xs, (f1, f2, f2, f1, f1, NaN)); append!(ys, (g1, g1, g2, g2, g1, NaN))
                i = k + 1
            else
                i += 1
            end
        end
    end
    return xs, ys
end

# `uncertainty` of the limb plots:
# * :shape (default) – the 1σ band of the shape alone: the drawn limbs placed on the fitted
#   centre, pointwise 16–84 %; its half-width is σ_a along the major and σ_b along the minor
#   axis, as the quoted errors; the 1σ ellipse of the centre is drawn separately
# * :region (synonyms :simultaneous, :band) – the simultaneous 68.3 % region of whole limbs
#   (centre and shape): 68.3 % of the limbs lie entirely inside
# * :pointwise – pointwise 16–84 % of whole limbs (centre and shape); :draws; :none / false
_simultaneous(u) = u in (:region, :simultaneous, :band)
function _band_label(u, nsigma=1; composite::Bool=false)
    pc = string(round(100 * SpecialFunctions.erf(nsigma / sqrt(2)); digits=1))      # 68.3, 95.4, 99.7
    sg = isinteger(nsigma) ? string(Int(nsigma)) : string(nsigma)
    composite && return u === :shape ? "$(sg)σ shape (occupancy)" : "$(sg)σ band (occupancy)"
    return _simultaneous(u) ? "$pc % ($(sg)σ) region" : u === :shape ? "$(sg)σ shape (a, b)" : "$pc % band (pointwise, $(sg)σ)"
end

# nσ error ellipse of the centre: the spread (covariance) of the drawn centres cs, drawn
# around the plotted centre μ (the cross) – not around their mean, which hard constraints
# (rejected draws) or an asymmetric posterior shift away from the fitted centre
function _centre_ellipse(cs, μ; nsigma=1)
    Σ = cov(hcat(first.(cs), last.(cs)))
    L = cholesky(Symmetric(Σ + 1e-12I); check=false)
    issuccess(L) || return nothing
    t = range(0, 2π; length=121)
    P = nsigma .* (L.L * vcat(cos.(t)', sin.(t)'))
    return μ[1] .+ P[1, :], μ[2] .+ P[2, :]
end

# distance from the origin of a closed curve (x, y) in the directions θ (from +g towards +f)
function _polar_radius(x, y, θ)
    φ = atan.(x, y); r = hypot.(x, y)
    o = sortperm(φ); φ, r = φ[o], r[o]
    φ = [φ[end] - 2π; φ; φ[1] + 2π]; r = [r[end]; r; r[1]]
    return map(θ) do t
        t = mod(t + π, 2π) - π
        j = clamp(searchsortedlast(φ, t), 1, length(φ) - 1)
        w = φ[j+1] == φ[j] ? 0.0 : (t - φ[j]) / (φ[j+1] - φ[j])
        (1 - w) * r[j] + w * r[j+1]
    end
end

# dimensions of a fitted ellipse or circle with 1σ errors, as lines of text; `get(name)`
# returns (value, σ) or (value, −σ, +σ), or nothing
function _dims_lines(get)
    fmt(x, f=1.0) = length(x) == 2 ? @sprintf("%.1f ± %.1f", f * x[1], f * x[2]) :
                                      @sprintf("%.1f −%.1f +%.1f", f * x[1], f * x[2], f * x[3])
    a, b, req, R = get(:equatorial_radius), get(:semi_minor_axis), get(:equivalent_radius), get(:radius)
    sep = get(:separation)
    l1, l2 = get(:lobe1_equivalent_radius), get(:lobe2_equivalent_radius)
    v1, v2, rv = get(:lobe1_volume_equivalent_radius), get(:lobe2_volume_equivalent_radius), get(:volume_equivalent_radius)
    if l1 !== nothing && l2 !== nothing && req !== nothing            # contact binary, 2-D
        return ["D_eq = $(fmt(req, 2.0)) km", "lobe 1: D_eq = $(fmt(l1, 2.0)) km", "lobe 2: D_eq = $(fmt(l2, 2.0)) km",
                "length = $(fmt(get(:length))) km"]
    elseif v1 !== nothing && v2 !== nothing && rv !== nothing          # contact binary, 3-D
        return ["D_V = $(fmt(rv, 2.0)) km", "lobe 1: D_V = $(fmt(v1, 2.0)) km", "lobe 2: D_V = $(fmt(v2, 2.0)) km",
                "length = $(fmt(get(:length))) km"]
    end
    out = String[]
    if a !== nothing && b !== nothing && req !== nothing
        out = ["a′ = $(fmt(a)) km", "b′ = $(fmt(b)) km", "D_eq = $(fmt(req, 2.0)) km"]
    elseif R !== nothing
        out = ["R = $(fmt(R)) km", "D = $(fmt(R, 2.0)) km"]
    end
    if sep !== nothing                                                # binary: the secondary
        sa, sr, sq = get(:sec_equatorial_radius), get(:sec_radius), get(:sec_equivalent_radius)
        sa !== nothing && sq !== nothing && push!(out, "sec. a′ = $(fmt(sa)) km", "sec. D_eq = $(fmt(sq, 2.0)) km")
        sr !== nothing && push!(out, "sec. D = $(fmt(sr, 2.0)) km")
        push!(out, "separation = $(fmt(sep)) km")
    end
    return out
end
# the dimensions of a given outline (e.g. the simulated truth), without errors; of a composite
# outline: two lobes if the parts overlap, otherwise a primary and its satellite
_outer_radius(O::EllipseOutline) = O.a
_outer_radius(O::Outline) = maximum(limb_radius(O, θ) for θ in range(0, 2π; length=181))
function _dims_lines(O::Outline)
    fm(x) = @sprintf("%.1f", x)
    if O isa CompositeOutline && length(O.parts) == 2
        (p1, p2), (o1, o2) = O.parts, O.offsets
        d = hypot((o2 .- o1)...)
        if d < _outer_radius(p1) + _outer_radius(p2)
            return ["D_eq = $(fm(2equivalent_radius(O))) km", "lobe 1: D_eq = $(fm(2equivalent_radius(p1))) km",
                    "lobe 2: D_eq = $(fm(2equivalent_radius(p2))) km"]
        end
        return [_dims_lines(p1); "sec. D_eq = $(fm(2equivalent_radius(p2))) km"; "separation = $(fm(d)) km"]
    elseif O isa EllipseOutline
        return ["a′ = $(fm(O.a)) km", "b′ = $(fm(O.b)) km", "D_eq = $(fm(2equivalent_radius(O))) km"]
    end
    return ["D_eq = $(fm(2equivalent_radius(O))) km"]
end

# the limb allowance of a limb problem, as an extra line (if > 0)
_slimb_line(p) = p isa LimbProblem && p.σ_limb > 0 ? [@sprintf("σ_limb = %.1f km", p.σ_limb)] : String[]
_dims_lines(r::FitResult) = (l = _dims_lines(n -> try r[n] catch; nothing end); isempty(l) ? l : [l; _slimb_line(r.problem)])
function _dims_lines(post::Posterior; rng=Random.default_rng())
    S = samples(post); ds = derived_samples(post; rng=rng)
    q(x) = (m = quantile(x, 0.5); (m, m - quantile(x, 0.15865525393145707), quantile(x, 0.8413447460685429) - m))
    l = _dims_lines(n -> (k = findfirst(==(n), post.names); k !== nothing ? q(S[:, k]) :
                                haskey(ds, n) ? q(ds[n]) : nothing))
    return isempty(l) ? l : [l; _slimb_line(post.problem)]
end


RecipesBase.@recipe function f(r::LimbFit; event=1, uncertainty=:shape, nsigma=1, ndraw=100, numbers=true, dimensions=true,
                              extend=false, compass=true,
                              rng=Random.default_rng())
    left_margin --> (5, :mm)                      # room for the axis labels (an outer legend narrows the plot)
    bottom_margin --> (3, :mm)
    p = r.problem
    k = _event_index(p, event)
    ev = p.setup.events[k]
    dts = Dict(c => v[1] for (c, v) in clock_offsets(r))
    aspect_ratio --> :equal
    framestyle --> :box
    xguide --> "f (km)"
    yguide --> "g (km)"
    xflip --> true
    title --> ev.name
    legend --> :outertopright
    (c0, O0) = outlines(r)[k]
    rings = _rings(p.model, r.values, k, p.setup)
    box = _ring_box(_limb_box(c0, O0), c0, rings)
    comp = O0 isa CompositeOutline
    xlims --> (box[1], box[2])
    ylims --> (box[3], box[4])
    if uncertainty !== false && uncertainty !== :none && all(isfinite, r.covariance)
        F = _Free(parameters(p))
        C = Symmetric(r.covariance + 1e-12I)
        L = cholesky(C; check=false)
        if issuccess(L)
            draws = Tuple{Vector{Float64},Vector{Float64}}[]
            limbs = Tuple{Tuple{Float64,Float64},Outline}[]
            centres = Tuple{Float64,Float64}[]
            for _ in 1:(uncertainty === :draws ? ndraw : max(ndraw, nsigma > 2 ? 4000 : 400))
                v = copy(r.values)
                v[F.idx] .+= L.L * randn(rng, length(F.idx))
                # the Gaussian of the covariance ignores the constraints: keep admissible draws
                p isa LimbProblem && !_satisfies_constraints(p, v) && continue
                try
                    (c, O) = outlines(p, v)[k]
                    cc = uncertainty === :shape ? c0 : c
                    push!(draws, outline_points(O; center=cc, n=361))
                    push!(limbs, ((cc[1], cc[2]), O))
                    push!(centres, (c[1], c[2]))
                catch
                end
            end
            if uncertainty === :draws
                xs = reduce(vcat, [[x; NaN] for (x, _) in draws]; init=Float64[])
                ys = reduce(vcat, [[y; NaN] for (_, y) in draws]; init=Float64[])
                @RecipesBase.series begin
                    seriestype := :path
                    linecolor --> :gray50
                    linealpha --> 0.15
                    label --> "limbs drawn from the covariance"
                    xs, ys
                end
            elseif length(draws) >= 20 && comp
                BX, BY = _occupancy_band(limbs, box; nsigma=nsigma)
                @RecipesBase.series begin
                    seriestype := :shape
                    fillcolor --> :gray45
                    fillalpha --> 0.4
                    linewidth --> 0
                    linealpha --> 0
                    label := ""
                    BX, BY
                end
                @RecipesBase.series begin                      # its legend entry
                    seriestype := :shape
                    fillcolor --> :gray45
                    fillalpha --> 0.4
                    linewidth --> 0
                    linealpha --> 0
                    label --> _band_label(uncertainty, nsigma; composite=true)
                    [NaN], [NaN]
                end
            elseif length(draws) >= 20
                @RecipesBase.series begin
                    seriestype := :shape
                    fillcolor --> :gray45
                    fillalpha --> 0.4
                    linewidth --> 0
                    linealpha --> 0
                    label --> _band_label(uncertainty, nsigma)
                    _limb_band(draws, c0; simultaneous=_simultaneous(uncertainty), nsigma=nsigma)
                end
                ce = uncertainty === :shape ? _centre_ellipse(centres, c0; nsigma=nsigma) : nothing
                if ce !== nothing
                    @RecipesBase.series begin
                        seriestype := :path
                        linecolor --> :black
                        linewidth --> 1
                        label --> (nsigma == 1 ? "1σ centre" : "$(nsigma)σ centre")
                        ce
                    end
                end
            end
        end
    end
    for (x, y, kw) in _chord_series!([], ev; dts=dts, numbers=numbers, ref=c0, box=box, extend=extend,
                                     compass=compass ? _plot_box(plotattributes, box) : nothing)
        @RecipesBase.series begin
            seriestype := :path
            for (kk, v) in pairs(kw)
                plotattributes[kk] = v
            end
            x, y
        end
    end
    @RecipesBase.series begin
        seriestype := :path
        linecolor --> :black
        linewidth --> 2
        label --> "fitted limb"
        outline_points(O0; center=c0)
    end
    for (j, rg) in enumerate(rings)
        @RecipesBase.series begin
            seriestype := :path
            linecolor --> :darkorange
            linewidth --> 1
            label --> (j == 1 ? "fitted rings" : "")
            outline_points(rg; center=c0)
        end
    end
    @RecipesBase.series begin
        seriestype := :scatter
        markershape --> :cross
        markercolor --> :black
        label := ""
        [c0[1]], [c0[2]]
    end
    if comp                                     # the centres of the components
        @RecipesBase.series begin
            seriestype := :scatter
            markershape --> :xcross
            markercolor --> :black
            markersize --> 3
            label := ""
            [c0[1] + o[1] for o in O0.offsets], [c0[2] + o[2] for o in O0.offsets]
        end
    end
    # dimensions as text in the right margin, below the legend (annotations are not clipped)
    dl = dimensions ? _dims_lines(r) : String[]
    if !isempty(dl)
        @RecipesBase.series begin
            seriestype := :scatter
            markersize := 0
            markeralpha := 0
            markerstrokewidth := 0
            label := ""
            series_annotations := [(l, 9, :left, :black) for l in dl]
            _dims_xy(box, length(dl))
        end
    end
end

RecipesBase.@recipe function f(post::LimbPosterior; event=1, uncertainty=:shape, nsigma=1, ndraw=400, numbers=true, dimensions=true,
                                    extend=false, compass=true,
                              rng=Random.default_rng())
    left_margin --> (5, :mm)                      # room for the axis labels (an outer legend narrows the plot)
    bottom_margin --> (3, :mm)
    p = post.problem
    k = _event_index(p, event)
    ev = p.setup.events[k]
    S = samples(post)
    aspect_ratio --> :equal
    framestyle --> :box
    xguide --> "f (km)"
    yguide --> "g (km)"
    xflip --> true
    title --> ev.name * " (posterior)"
    legend --> :outertopright
    vmed = copy(post.base); vmed[post.idx] .= vec(median(S; dims=1))
    (c0, O0) = outlines(p, vmed)[k]
    box = _ring_box(_limb_box(c0, O0), c0, _rings(p.model, vmed, k, p.setup))
    comp = O0 isa CompositeOutline
    xlims --> (box[1], box[2])
    ylims --> (box[3], box[4])
    # the maximum-posterior limb and centre (the cross): the shape draws and the centre
    # ellipse are placed on it
    Lp = vec(post.logp[post.burn+1:end, :])
    vbest = copy(post.base); vbest[post.idx] .= S[argmax(Lp), :]
    (cb, Ob) = outlines(p, vbest)[k]
    centres = Tuple{Float64,Float64}[]
    limbs = Tuple{Tuple{Float64,Float64},Outline}[]
    draws = map(rand(rng, 1:size(S, 1), uncertainty === :draws ? ndraw : max(ndraw, nsigma > 2 ? 4000 : 400))) do i
        v = copy(post.base); v[post.idx] .= S[i, :]
        (c, O) = outlines(p, v)[k]
        push!(centres, (c[1], c[2]))
        cc = uncertainty === :shape ? cb : c
        push!(limbs, ((cc[1], cc[2]), O))
        outline_points(O; center=cc, n=361)
    end
    if uncertainty === :draws
        @RecipesBase.series begin
            seriestype := :path
            linecolor --> :gray50
            linealpha --> 0.05
            label --> "limbs drawn from the posterior"
            reduce(vcat, [[x; NaN] for (x, _) in draws]), reduce(vcat, [[y; NaN] for (_, y) in draws])
        end
    elseif uncertainty !== false && uncertainty !== :none && comp
        BX, BY = _occupancy_band(limbs, box; nsigma=nsigma)
        @RecipesBase.series begin
            seriestype := :shape
            fillcolor --> :gray45
            fillalpha --> 0.4
            linewidth --> 0
            linealpha --> 0
            label := ""
            BX, BY
        end
        @RecipesBase.series begin                      # its legend entry
            seriestype := :shape
            fillcolor --> :gray45
            fillalpha --> 0.4
            linewidth --> 0
            linealpha --> 0
            label --> _band_label(uncertainty, nsigma; composite=true)
            [NaN], [NaN]
        end
    elseif uncertainty !== false && uncertainty !== :none
        @RecipesBase.series begin
            seriestype := :shape
            fillcolor --> :gray45
            fillalpha --> 0.4
            linewidth --> 0
            linealpha --> 0
            label --> _band_label(uncertainty, nsigma)
            _limb_band(draws, uncertainty === :shape ? cb : c0; simultaneous=_simultaneous(uncertainty), nsigma=nsigma)
        end
        ce = uncertainty === :shape ? _centre_ellipse(centres, cb; nsigma=nsigma) : nothing
        if ce !== nothing
            @RecipesBase.series begin
                seriestype := :path
                linecolor --> :black
                linewidth --> 1
                label --> (nsigma == 1 ? "1σ centre" : "$(nsigma)σ centre")
                ce
            end
        end
    end
    for (x, y, kw) in _chord_series!([], ev; numbers=numbers, ref=c0, box=box, extend=extend,
                                     compass=compass ? _plot_box(plotattributes, box) : nothing)
        @RecipesBase.series begin
            seriestype := :path
            for (kk, v) in pairs(kw)
                plotattributes[kk] = v
            end
            x, y
        end
    end
    # the maximum-posterior limb (computed above)
    @RecipesBase.series begin
        seriestype := :path
        linecolor --> :black
        linewidth --> 2
        label --> "maximum-posterior limb"
        outline_points(Ob; center=cb)
    end
    for (j, rg) in enumerate(_rings(p.model, vbest, k, p.setup))
        @RecipesBase.series begin
            seriestype := :path
            linecolor --> :darkorange
            linewidth --> 1
            label --> (j == 1 ? "maximum-posterior rings" : "")
            outline_points(rg; center=cb)
        end
    end
    @RecipesBase.series begin
        seriestype := :scatter
        markershape --> :cross
        markercolor --> :black
        label := ""
        [cb[1]], [cb[2]]
    end
    if comp                                     # the centres of the components
        @RecipesBase.series begin
            seriestype := :scatter
            markershape --> :xcross
            markercolor --> :black
            markersize --> 3
            label := ""
            [cb[1] + o[1] for o in Ob.offsets], [cb[2] + o[2] for o in Ob.offsets]
        end
    end
    # dimensions as text in the right margin, below the legend (annotations are not clipped)
    dl = dimensions ? _dims_lines(post; rng=rng) : String[]
    if !isempty(dl)
        @RecipesBase.series begin
            seriestype := :scatter
            markersize := 0
            markeralpha := 0
            markerstrokewidth := 0
            label := ""
            series_annotations := [(l, 9, :left, :black) for l in dl]
            _dims_xy(box, length(dl))
        end
    end
end

"""
    cornerplot(posterior; vars=posterior.names, axislabels=Dict(), truths=Dict(), truth_label="reference",
               bins=40, maxpoints=4000)                                    (with Plots loaded)

Corner plot of the posterior samples: histograms with the median and the 16 % / 84 %
quantiles (titled `median −σ +σ`), pairwise scatter plots with the 1σ and 2σ credible
regions (39 % and 86 % of the samples, 2-D), the maximum-posterior sample, optional reference
values `truths` (`name => value`, e.g. published results; an orange vertical line in the
histograms and an orange diamond in the 2-D panels, labelled `truth_label` in the legend)
and, in the top-right corner, a legend with the number of samples, the acceptance rate and
the autocorrelation time. `axislabels` maps parameter names to axis labels (e.g.
`:a => "a [km]"`). The title is the event name(s) (`plot_title` overrides it).
"""
RecipesBase.@userplot CornerPlot

function _hpd_levels(H, fracs)
    v = sort(vec(H); rev=true)
    c = cumsum(v) ./ max(sum(v), eps())
    return [v[clamp(searchsortedfirst(c, f), 1, length(v))] for f in fracs]
end

function _hist2d(x, y, bins)
    xe = range(minimum(x), maximum(x); length=bins + 1)
    ye = range(minimum(y), maximum(y); length=bins + 1)
    H = zeros(bins, bins)
    for (u, w) in zip(x, y)
        i = clamp(searchsortedlast(xe, u), 1, bins)
        j = clamp(searchsortedlast(ye, w), 1, bins)
        H[i, j] += 1
    end
    return (xe[1:end-1] .+ xe[2:end]) ./ 2, (ye[1:end-1] .+ ye[2:end]) ./ 2, H
end

RecipesBase.@recipe function f(cp::CornerPlot; vars=nothing, axislabels=Dict(), truths=Dict(), truth_label="reference",
                               bins=40, maxpoints=4000)
    left_margin --> (5, :mm)                      # room for the axis labels (an outer legend narrows the plot)
    bottom_margin --> (3, :mm)
    post = cp.args[1]
    names = vars === nothing ? post.names : collect(vars)
    S = samples(post)
    for k in eachindex(post.names)
        _periodic_center!(view(S, :, k), post, k)
    end
    idx = [findfirst(==(n), post.names) for n in names]
    any(isnothing, idx) && throw(ArgumentError("unknown parameter in vars; parameters are $(post.names)"))
    n = length(idx)
    L = vec(post.logp[post.burn+1:end, :])
    imap = argmax(L)
    rows = size(S, 1) > maxpoints ? round.(Int, range(1, size(S, 1); length=maxpoints)) : (1:size(S, 1))
    lab(k) = string(get(axislabels, names[k], names[k]))
    tr(k) = get(truths, names[k], nothing)
    layout := (n, n)
    size --> (280n, 280n)
    evs = post.problem isa LimbProblem ? events(post.problem) : Event[]
    plot_title --> (isempty(evs) ? "posterior" : join(unique(e.name for e in evs), "; ") * ": posterior")
    plot_titlefontsize --> 12
    q = (0.15865525393145707, 0.5, 0.8413447460685429)
    # subplot attributes have to be given with every series of the subplot
    function sub(sp; kw...)
        d = Dict{Symbol,Any}(:subplot => sp, :tickfontsize => 6, :guidefontsize => 8, :xrotation => 45,
                             :legend => false, :titlefontsize => 8)
        for (k, v) in kw
            d[k] = v
        end
        return d
    end
    for i in 1:n, j in 1:n
        sp = (i - 1) * n + j
        xi, xj = S[:, idx[i]], S[:, idx[j]]
        if j > i
            if i == 1 && j == n && n > 1             # legend and chain statistics
                τ = maximum(autocorrelation_time(post))
                info = @sprintf("%d samples\n%d walkers × %d steps\nacceptance %.2f\nτ ≤ %.0f steps",
                                size(S, 1), size(post.chain, 2), size(post.chain, 1) - post.burn, post.acceptance, τ)
                A = sub(sp; framestyle=:none, legend=:top, legendfontsize=8, xlims=(0, 1), ylims=(0, 1),
                        annotations=[(0.08, 0.12, (info, 8, :left))])
                entries = Any[("samples", (seriestype=:scatter, markersize=3, markercolor=:gray40, markerstrokewidth=0)),
                              ("1σ, 2σ credible regions", (seriestype=:path, linecolor=:black)),
                              ("median", (seriestype=:path, linecolor=:steelblue)),
                              ("16 %, 84 % quantiles", (seriestype=:path, linecolor=:steelblue, linestyle=:dash)),
                              ("maximum posterior", (seriestype=:scatter, markershape=:star5, markercolor=:red, markersize=7))]
                ## the reference values (e.g. published): a vertical line in the histograms, a diamond
                ## in the 2-D panels
                isempty(truths) || push!(entries, (truth_label, (seriestype=:path, linecolor=:darkorange, linewidth=2,
                                                                 markershape=:diamond, markercolor=:darkorange,
                                                                 markersize=6)))
                for (lbl, kw) in entries
                    @RecipesBase.series begin
                        for (kk, vv) in A
                            plotattributes[kk] = vv
                        end
                        label := lbl
                        for (kk, vv) in pairs(kw)
                            plotattributes[kk] = vv
                        end
                        [NaN], [NaN]
                    end
                end
            else
                @RecipesBase.series begin
                    for (kk, vv) in sub(sp; framestyle=:none)
                        plotattributes[kk] = vv
                    end
                    seriestype := :scatter
                    markersize := 0
                    label := ""
                    [NaN], [NaN]
                end
            end
        elseif i == j
            lo, md, hi = Statistics.quantile(xi, collect(q))
            A = sub(sp; title=@sprintf("%s = %.4g −%.2g +%.2g", lab(i), md, md - lo, hi - md),
                    xguide=(i == n ? lab(j) : ""), yticks=false)
            @RecipesBase.series begin
                for (kk, vv) in A
                    plotattributes[kk] = vv
                end
                seriestype := :stephist
                bins := bins
                normalize := :pdf
                linecolor := :black
                label := ""
                xi
            end
            for (x, ls) in ((lo, :dash), (md, :solid), (hi, :dash))
                @RecipesBase.series begin
                    for (kk, vv) in A
                        plotattributes[kk] = vv
                    end
                    seriestype := :vline
                    linecolor := :steelblue
                    linestyle := ls
                    label := ""
                    [x]
                end
            end
            if tr(i) !== nothing
                @RecipesBase.series begin
                    for (kk, vv) in A
                        plotattributes[kk] = vv
                    end
                    seriestype := :vline
                    linecolor := :darkorange
                    linewidth := 2
                    label := ""
                    [tr(i)]
                end
            end
        else
            A = sub(sp; xguide=(i == n ? lab(j) : ""), yguide=(j == 1 ? lab(i) : ""))
            @RecipesBase.series begin
                for (kk, vv) in A
                    plotattributes[kk] = vv
                end
                seriestype := :scatter
                markersize := 1
                markerstrokewidth := 0
                markeralpha := 0.25
                markercolor := :gray40
                label := ""
                xj[rows], xi[rows]
            end
            xc, yc, H = _hist2d(xj, xi, bins)
            l2, l1 = _hpd_levels(H, (0.8647, 0.3935))
            if l1 > l2 > 0
                @RecipesBase.series begin
                    for (kk, vv) in A
                        plotattributes[kk] = vv
                    end
                    seriestype := :contour
                    levels := [l2, l1]
                    linecolor := :black
                    linewidth := 1.2
                    colorbar := false
                    label := ""
                    xc, yc, permutedims(H)
                end
            end
            @RecipesBase.series begin
                for (kk, vv) in A
                    plotattributes[kk] = vv
                end
                seriestype := :scatter
                markershape := :star5
                markercolor := :red
                markersize := 6
                label := ""
                [xj[imap]], [xi[imap]]
            end
            if tr(i) !== nothing && tr(j) !== nothing
                @RecipesBase.series begin
                    for (kk, vv) in A
                        plotattributes[kk] = vv
                    end
                    seriestype := :scatter
                    markershape := :diamond
                    markercolor := :darkorange
                    markersize := 5
                    label := ""
                    [tr(j)], [tr(i)]
                end
            end
        end
    end
end

RecipesBase.@recipe function f(o::Outline; center=(0.0, 0.0))
    left_margin --> (5, :mm)                      # room for the axis labels (an outer legend narrows the plot)
    bottom_margin --> (3, :mm)
    framestyle --> :box
    aspect_ratio --> :equal
    xflip --> true
    xguide --> "f (km)"
    yguide --> "g (km)"
    label --> ""
    outline_points(o; center=center)
end

RecipesBase.@recipe function f(r::ProjectedRing; center=(0.0, 0.0))
    aspect_ratio --> :equal
    seriestype := :path
    label --> "ring"
    outline_points(r; center=center)
end

RecipesBase.@recipe function f(p::AtmosphereProfile; ρmax=nothing)
    left_margin --> (5, :mm)                      # room for the axis labels (an outer legend narrows the plot)
    bottom_margin --> (3, :mm)
    framestyle --> :box
    xguide --> "distance from shadow centre [km]"
    yguide --> "flux"
    hi = ρmax === nothing ? max(p.y_surface, 0.0) + 2 * (p.R - min(p.y_surface, p.R)) + 10p.scale_height : ρmax
    ρs = range(0.0, min(hi, p.y[end] * 1.2); length=8000)
    @RecipesBase.series begin
        label --> "refraction"
        ρs, [Occultations._refracted_flux(p, ρ) for ρ in ρs]
    end
    @RecipesBase.series begin
        label --> "with surface edge and stellar disc"
        ρs, [atmosphere_flux(p, ρ) for ρ in ρs]
    end
end

RecipesBase.@recipe function f(p::ShadowPath; observers=Observer[])
    left_margin --> (5, :mm)                      # room for the axis labels (an outer legend narrows the plot)
    bottom_margin --> (3, :mm)
    framestyle --> :box
    xguide --> "longitude (deg)"
    yguide --> "latitude (deg)"
    title --> "$(p.body) – $(isot(p.times[1]; digits=0))"
    aspect_ratio --> :equal
    for (k, δ) in enumerate(p.offsets)
        @RecipesBase.series begin
            seriestype := :path
            linestyle --> (δ == 0 ? :solid : :dash)
            linecolor --> (δ == 0 ? :blue : :black)
            label --> (δ == 0 ? "centre" : "offset $(round(δ; digits=1)) km")
            p.lon[:, k], p.lat[:, k]
        end
    end
    if !isempty(observers)
        annotations --> [(o.lon, o.lat, ("  " * o.name, 8, :left)) for o in observers]
        @RecipesBase.series begin
            seriestype := :scatter
            label --> "observers"
            [o.lon for o in observers], [o.lat for o in observers]
        end
    end
end
