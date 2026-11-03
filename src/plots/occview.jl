module OccViews

# The occultation in 3-D: a 3-D limb fit (ellipsoid or mesh) or its posterior shown as the
# body oriented at the event, with the chords, the Sun, the rotational light curve, sliders
# for the phase, the event and the parameters (live χ²), and animations.

using ..Occultations
using ..Occultations: LimbFit, LimbPosterior, LimbProblem, EllipsoidLimb, MeshLimb, ContactBinaryEllipsoids, Ellipsoid,
                      BilobedEllipsoid, TriMesh, ellipsoid_mesh,
                      chord_path, isnegative, ispositive, outline_points, project, area, earth_heliocentric, samples
using LinearAlgebra
using Printf
using Random
using Makie
import ..Occultations: occview, occanimate

const O = Occultations
const GB = Makie.GeometryBasics
const TriangleFace = GB.TriangleFace

# ------------------------------------------------------------------------------------------
# geometry of event k for the parameter values v, with an extra rotation dW [deg]
# ------------------------------------------------------------------------------------------

function _check(p::LimbProblem)
    p.model isa Union{EllipsoidLimb,MeshLimb,ContactBinaryEllipsoids} ||
        throw(ArgumentError("occview needs a 3-D limb model (EllipsoidLimb, MeshLimb or ContactBinaryEllipsoids), not $(typeof(p.model))"))
end

# ICRS → body matrix of event k (rotation dW added to the prime meridian)
function _matrix(p::LimbProblem, v, k; dW=0.0)
    s, m = p.setup, p.model
    m isa MeshLimb && m.matrix !== nothing && return m.matrix * O._R3(-dW)'
    return O._body_matrix(O._v(v, s, :pole_ra), O._v(v, s, :pole_dec), O._meridian(v, s, k) + dW)
end

_mesh(p::LimbProblem, v) = (sh = O._shape(p.model, v, p.setup);
                             sh isa Union{Ellipsoid,BilobedEllipsoid} ? ellipsoid_mesh(sh; subdivisions=4) : sh)

# in-plane rotation of a mesh limb (parameter `rotation`, deg) – sign as in `transform`
function _inplane(p::LimbProblem, v, k)
    m = p.model
    m isa MeshLimb || return 0.0
    haskey(p.setup.index, :rotation) || any(n -> startswith(string(n), "rotation"), keys(p.setup.index)) || return 0.0
    return O._v(v, p.setup, :rotation, k, :rotation in m.per_event)
end

# sky coordinates (f east, g north, w away from the observer) [km] of the mesh vertices, relative
# to the body centre; vertex normals; Sun direction and pole in sky coordinates
function _geometry(p::LimbProblem, v, k; dW=0.0)
    s = p.setup
    ctx = s.contexts[k]
    M = _matrix(p, v, k; dW=dW)
    S = O._sky_matrix(ctx.u)
    mesh = _mesh(p, v)
    X = S * (M' * mesh.vertices)
    ψ = _inplane(p, v, k)
    if ψ != 0
        c, sn = cosd(ψ), sind(ψ)
        X[1:2, :] = [c sn; -sn c] * X[1:2, :]
    end
    # vertex normals (sky frame) from the faces
    N = zeros(3, size(X, 2))
    for j in axes(mesh.faces, 2)
        a, b, c = mesh.faces[:, j]
        n = cross(X[:, b] - X[:, a], X[:, c] - X[:, a])
        N[:, a] .+= n; N[:, b] .+= n; N[:, c] .+= n
    end
    N ./= max.(sqrt.(sum(abs2, N; dims=1)), 1e-300)
    # make the normals point outwards
    dot(N[:, 1], X[:, 1]) < 0 && (N .*= -1)
    # Sun as seen from the body, in sky coordinates
    au = 1.495978707e8
    body = earth_heliocentric(ctx.tref) .+ ctx.u .* (ctx.dist / au)
    sun = S * (-body ./ norm(body))
    pole = S * (M' * [0.0, 0.0, 1.0])
    return (X=X, F=mesh.faces, N=N, sun=sun, pole=pole, M=M, S=S)
end

# silhouette (outline points relative to the centre) of event k
_silhouette(p::LimbProblem, v, k; dW=0.0) = dW == 0 ? outline_points(O.outlines(p, v)[k][2]; center=(0.0, 0.0), n=361) :
    _outline_dW(p, v, k, dW)

function _outline_dW(p::LimbProblem, v, k, dW)
    s = p.setup
    sh = O._shape(p.model, v, s)
    o = project(sh, _matrix(p, v, k; dW=dW), s.contexts[k].u)
    ψ = _inplane(p, v, k)
    ψ != 0 && (o = O.transform(o; rotation=ψ))
    return outline_points(o; center=(0.0, 0.0), n=361)
end

# sense of rotation: +1 if W increases with time (prograde period), −1 otherwise
_sense(p::LimbProblem, v) = haskey(p.setup.index, :period) ? sign(O._v(v, p.setup, :period)) : 1.0

# rotational light curve through one rotation, as magnitudes from the projected area, versus the
# time offset x from the event in rotations (x ∈ [−1, 0]); and the phase of the event: the
# rotations elapsed since the last maximum of the light curve (maximum projected area) before it
function _lightcurve(p::LimbProblem, v, k; n=721)
    s = p.setup
    sh = O._shape(p.model, v, s)
    sg = _sense(p, v)
    x = range(-1, 0; length=n)
    A = [area(project(sh, _matrix(p, v, k; dW=360 * sg * xi), s.contexts[k].u)) for xi in x]
    Ā = sum(A) / length(A)
    # local maxima (with a tolerance for the two equal maxima of an ellipsoid); the last before the event
    Amax = maximum(A)
    peaks = [i for i in 2:n-1 if A[i] >= A[i-1] && A[i] >= A[i+1] && A[i] >= Amax - 1e-3 * (Amax - minimum(A))]
    i0 = isempty(peaks) ? argmax(A) : last(peaks)
    φev = -x[i0]                                   # rotations since that maximum
    # the light curve over phase φ ∈ [0, 1): φ = (x + φev) mod 1
    φ = mod.(collect(x) .+ φev, 1.0)
    o = sortperm(φ)
    return φ[o], (-2.5 .* log10.(A ./ Ā))[o], φev
end

# the chords of event k relative to the fitted centre: (f, g, kind) per chord, clipped to
# |f|, |g| ≤ box (NaN outside)
function _chords(p::LimbProblem, v, k; box=Inf)
    ev = p.setup.events[k]
    cf, cg = v[p.centers[k][1]], v[p.centers[k][2]]
    out = []
    for c in ev.chords
        if isnegative(c)
            f, g = chord_path(ev, c; n=2000)
            push!(out, (_clip(f .- cf, g .- cg, box)..., :negative, c.name))
        else
            f, g = chord_path(ev, c; segment=:occulted)
            isempty(f) && continue                       # only ring mid-times
            push!(out, (_clip(f .- cf, g .- cg, box)..., :positive, c.name))
        end
    end
    return out
end

function _clip(f, g, box)
    ok = (abs.(f) .<= box) .& (abs.(g) .<= box)
    return ifelse.(ok, f, NaN), ifelse.(ok, g, NaN)
end

χ2(p::LimbProblem, v) = (r = O.residuals(p, v); r === nothing ? Inf : sum(abs2, r))

# ------------------------------------------------------------------------------------------
# the figure
# ------------------------------------------------------------------------------------------

function _occview(p::LimbProblem, v0::Vector{Float64}; event::Integer=1, explore::Bool=false, ensemble=Vector{Float64}[],
                  title::AbstractString="", size=(1500, 950), body_color=:gray80, night_color=:gray25)
    _check(p)
    s = p.setup
    nev = length(s.events)
    fig = Figure(; size=size)
    Label(fig[0, 1:2], isempty(title) ? s.events[event].name : title; fontsize=18, font=:bold)
    ax3 = Axis3(fig[1, 1]; aspect=:data, xreversed=true, elevation=π / 2 - 0.001, azimuth=-π / 2,
                xlabel="f (km)", ylabel="g (km)", zlabel="towards observer (km)",
                title="as seen from the Earth (drag to rotate with GLMakie)")
    right = fig[1, 2] = GridLayout()
    ax2 = Axis(right[1, 1]; aspect=DataAspect(), xreversed=true, xlabel="f (km)", ylabel="g (km)",
               title="sky plane")
    Occultations.MeshViews._sky_compass!(ax2)
    axl = Axis(right[2, 1]; xlabel="rotational phase (0 = maximum before the event)", ylabel="Δmag (projected area; brighter ↑)",
               yreversed=true,
               title="rotational light curve")
    rowsize!(right, 2, Relative(0.35))
    colsize!(fig.layout, 1, Relative(0.55))

    # controls
    controls = fig[2, 1:2] = GridLayout()
    sl = Any[(label="rotational phase (0 = maximum before the event)", range=0:0.002:1, startvalue=0.0)]
    nev > 1 && push!(sl, (label="event", range=1:nev, startvalue=event))
    pnames = Symbol[]
    if explore
        for (n, lo, hi) in ((:pole_ra, 0.0, 360.0), (:pole_dec, -90.0, 90.0), (:W0, 0.0, 360.0))
            haskey(s.index, n) || continue
            push!(pnames, n)
            lab = n === :pole_ra ? "pole RA [°]" : n === :pole_dec ? "pole Dec [°]" :
                  "W0 [°] (prime-meridian angle at the reference epoch)"
            push!(sl, (label=lab, range=range(lo, hi; length=361), startvalue=v0[s.index[n]]))
        end
        for n in (:a, :b, :c, :scale)
            haskey(s.index, n) || continue
            x = v0[s.index[n]]
            push!(pnames, n)
            push!(sl, (label=n === :scale ? "scale" : "$n [km] (semi-axis)", range=range(0.5x, 1.5x; length=201),
                       startvalue=x))
        end
    end
    sg = SliderGrid(controls[1, 1], sl...; tellheight=true)
    kk = Observable(event)
    nev > 1 && on(x -> (kk[] = Int(x)), sg.sliders[2].value)
    v = Observable(copy(v0))
    # rotational phase: 0 at the last light-curve maximum before the event; the phase of the event
    lc = lift((v, k) -> _lightcurve(p, v, k), v, kk)
    φev = lift(x -> x[3], lc)
    φ = Observable(φev[])                                     # displayed rotational phase
    set_close_to!(sg.sliders[1], φev[])
    on(x -> (φ[] = Float64(x)), sg.sliders[1].value)
    # rotation applied to the event's orientation [deg]: (φ − φ_event) rotations, in the sense of rotation
    dW = lift((f, fe, v) -> 360 * _sense(p, v) * (f - fe), φ, φev, v)
    if explore
        off = nev > 1 ? 2 : 1
        for (j, n) in enumerate(pnames)
            on(sg.sliders[off+j].value) do x
                w = copy(v[]); w[s.index[n]] = x; v[] = w
            end
        end
    end
    info = lift(v, kk, φ, φev) do v, k, f, fe
        c = χ2(p, v)
        @sprintf("rotational phase %.3f (event: %.3f)\nχ² = %s (all events)   event %d of %d", f, fe,
                 isfinite(c) ? @sprintf("%.2f", c) : "– (inadmissible)", k, nev)
    end
    Label(controls[1, 2], info; tellwidth=true, fontsize=14)
    # reset: all sliders (and with them the parameters) back to their start values
    reset = Button(controls[2, 2]; label="reset", tellwidth=false)
    on(reset.clicks) do _
        v[] = copy(v0); kk[] = event
        for (sld, spec) in zip(sg.sliders, sl)
            set_close_to!(sld, spec.startvalue)
        end
        set_close_to!(sg.sliders[1], φev[]); φ[] = φev[]
    end

    # 3-D: the body, lit/night
    geo = lift((v, k, d) -> _geometry(p, v, k; dW=d), v, kk, dW)
    pts = lift(g -> [Point3f(g.X[1, i], g.X[2, i], -g.X[3, i]) for i in axes(g.X, 2)], geo)
    fcs = lift(g -> [TriangleFace(g.F[1, j], g.F[2, j], g.F[3, j]) for j in axes(g.F, 2)], geo)
    cols = lift(g -> [dot(g.N[:, i], g.sun) > 0 ? Makie.to_color(body_color) : Makie.to_color(night_color)
                      for i in axes(g.X, 2)], geo)
    mesh!(ax3, lift((P, F) -> Makie.GeometryBasics.Mesh(P, F), pts, fcs); color=cols, shading=Makie.NoShading)
    R = lift(g -> maximum(sqrt.(sum(abs2, g.X; dims=1))), geo)
    # pole, Sun, sub-observer and sub-solar points
    arr(dir, len) = lift((g, r) -> (d = dir(g); [Point3f(0, 0, 0), Point3f(len * r * d[1], len * r * d[2], -len * r * d[3])]), geo, R)
    lines!(ax3, arr(g -> g.pole, 1.5); color=:blue, linewidth=3, label="pole (+z of the body)")
    lines!(ax3, arr(g -> g.sun, 1.8); color=:orange, linewidth=3, label="Sun direction (from the centre)")
    sub(f) = lift(g -> (i = argmax([f(g, i) for i in axes(g.X, 2)]); [Point3f(g.X[1, i], g.X[2, i], -g.X[3, i])]), geo)
    scatter!(ax3, sub((g, i) -> -g.N[3, i]); color=:cyan, markersize=14, label="sub-observer point")
    scatter!(ax3, sub((g, i) -> dot(g.N[:, i], g.sun)); color=:yellow, strokecolor=:black, strokewidth=1, markersize=14,
             label="sub-solar point")
    # chords in front of the body
    chords = lift((v, k, r) -> _chords(p, v, k; box=2.2r), v, kk, R)
    front = lift(r -> 1.3r, R)
    for kind in (:positive, :negative)
        seg = lift(chords, front) do C, z
            P = Point3f[]
            for (f, g, kd, _) in C
                kd === kind || continue
                append!(P, [Point3f(f[i], g[i], z) for i in eachindex(f)]); push!(P, Point3f(NaN, NaN, NaN))
            end
            isempty(P) ? [Point3f(NaN, NaN, NaN)] : P
        end
        lines!(ax3, seg; color=kind === :positive ? :red : :gray40, linestyle=kind === :positive ? :solid : :dash,
               linewidth=kind === :positive ? 2.5 : 1)
    end
    leg = [LineElement(color=:blue, linewidth=3),
           MarkerElement(marker=:circle, color=:cyan, markersize=12),
           MarkerElement(marker=:circle, color=:yellow, strokecolor=:black, strokewidth=1, markersize=12),
           PolyElement(color=body_color), PolyElement(color=night_color),
           LineElement(color=:red, linewidth=2.5), LineElement(color=:gray40, linestyle=:dash)]
    labs = ["pole (+z of the body)", "sub-observer point", "sub-solar point",
            "sunlit side", "night side", "chords", "negative chords"]
    axislegend(ax3, leg, labs; position=:lb, labelsize=10)
    # phase angle (Sun–body–observer): the observer is at −w
    on(geo; update=true) do g
        ph = acosd(clamp(-g.sun[3], -1, 1))
        ax3.title[] = @sprintf("as seen from the Earth; phase angle %.2f° (Sun arrow ≈ towards the viewer when small)", ph)
    end

    # sky plane: silhouette(s) and chords
    for w in ensemble
        lines!(ax2, lift((k, d) -> (x = _silhouette(p, w, k; dW=d); Point2f.(x[1], x[2])), kk, dW); color=(:steelblue, 0.25),
               linewidth=0.8)
    end
    lines!(ax2, lift((v, k, d) -> (x = _silhouette(p, v, k; dW=d); Point2f.(x[1], x[2])), v, kk, dW); color=:black,
           linewidth=2, label=isempty(ensemble) ? "silhouette" : "maximum posterior")
    for kind in (:positive, :negative)
        seg = lift(chords) do C
            P = Point2f[]
            for (f, g, kd, _) in C
                kd === kind || continue
                append!(P, Point2f.(f, g)); push!(P, Point2f(NaN, NaN))
            end
            isempty(P) ? [Point2f(NaN, NaN)] : P
        end
        lines!(ax2, seg; color=kind === :positive ? :red : :gray40, linestyle=kind === :positive ? :solid : :dash,
               linewidth=kind === :positive ? 2 : 1, label=kind === :positive ? "chords" : "negative chords")
    end
    on(R; update=true) do r
        limits!(ax2, -2.2r, 2.2r, -2.2r, 2.2r)
        limits!(ax3, -2.2r, 2.2r, -2.2r, 2.2r, -2.2r, 2.2r)
    end
    axislegend(ax2; position=:lb, labelsize=10)

    # rotational light curve over the rotational phase
    lines!(axl, lift(x -> Point2f.(x[1], x[2]), lc); color=:black)
    vlines!(axl, lift(f -> [Float32(f)], φ); color=:red, label="displayed phase")
    vlines!(axl, lift(f -> [Float32(f)], φev); color=:gray50, linestyle=:dash, label="event")
    amp = lift(x -> maximum(x[2]) - minimum(x[2]), lc)
    obs = isempty(p.amplitudes) ? "" : " | observed Δm = " * join((@sprintf("%.3f ± %.3f", d.a.amplitude, d.a.σ) for d in p.amplitudes), ", ")
    on(a -> (axl.title[] = @sprintf("rotational light curve: model Δm = %.3f mag", a) * obs), amp; update=true)
    xlims!(axl, 0, 1)
    axislegend(axl; position=:rb, labelsize=10)
    return fig, φ, kk, v
end

_values(r::LimbFit) = copy(r.values)
function _values(post::LimbPosterior)
    S = samples(post)
    Lp = vec(post.logp[post.burn+1:end, :])
    v = copy(post.base); v[post.idx] .= S[argmax(Lp), :]
    return v
end
function _ensemble(post::LimbPosterior, ndraw, rng)
    S = samples(post)
    return [(w = copy(post.base); w[post.idx] .= S[i, :]; w) for i in rand(rng, 1:size(S, 1), ndraw)]
end

function occview(r::LimbFit; event::Integer=1, explore::Bool=false, kw...)
    return _occview(r.problem, _values(r); event=event, explore=explore, kw...)[1]
end
function occview(post::LimbPosterior; event::Integer=1, explore::Bool=false, ndraw::Integer=24,
                 rng=MersenneTwister(1), kw...)
    return _occview(post.problem, _values(post); event=event, explore=explore, ensemble=_ensemble(post, ndraw, rng),
                    kw...)[1]
end

function occanimate(x::Union{LimbFit,LimbPosterior}, file::AbstractString; event::Integer=1, frames::Integer=120,
                    framerate::Integer=24, ndraw::Integer=12, rng=MersenneTwister(1), kw...)
    p = x.problem
    ens = x isa LimbPosterior ? _ensemble(x, ndraw, rng) : Vector{Float64}[]
    fig, φ, _, _ = _occview(p, _values(x); event=event, ensemble=ens, kw...)
    record(fig, file, range(0, 1; length=frames); framerate=framerate) do f
        φ[] = f
    end
    return file
end

end
