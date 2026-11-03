module GeoMaps

# Cartopy-like maps of occultation paths with GeoMakie: PROJ projections, Natural Earth
# land, coastlines, country borders and capitals, the shadow band, time marks, observers

using ..Occultations
using ..Occultations: OccMap, isot, LimbFit
using GeoMakie
using Printf
import ..Occultations: geomap

const GB = GeoMakie.GeometryBasics
const NE = GeoMakie.NaturalEarth
const _CACHE = Dict{Tuple{String,Int},Any}()

_default_projection(m::OccMap) = m.projection === :ortho ?
    "+proj=ortho +lon_0=$(round(m.center[1]; digits=4)) +lat_0=$(round(m.center[2]; digits=4))" :
    "+proj=eqearth +lon_0=$(round(m.center[1]; digits=4))"

# Natural Earth layers (scale 110 is bundled with GeoMakie; 50 and 10 are downloaded once and
# cached by NaturalEarth.jl); `nothing` if unavailable
function _layer(name::String, scale::Int)
    get!(_CACHE, (name, scale)) do
        try
            name == "coastline" ? (scale == 110 ? GeoMakie.coastlines() : GeoMakie.coastlines(scale)) :
            name == "land" ? (scale == 110 ? GeoMakie.land() : GeoMakie.land(scale)) :
            NE.naturalearth(name, scale)
        catch e
            @warn "geomap: Natural Earth layer $name at 1:$(scale)m is not available ($(sprint(showerror, e)))" maxlog = 1
            nothing
        end
    end
end

# labels of the time marks: the first in full (hh:mm:ss UTC), the others only the minute (ᵐ)
# or, for marks less than a minute apart, the second (ˢ)
function _mark_labels(ts)
    isempty(ts) && return String[]
    full = [isot(t; digits=0)[12:19] for t in ts]
    sub = length(ts) > 1 && any(i -> full[i][7:8] != full[1][7:8], eachindex(full))
    return [i == 1 ? " " * full[i] * " UTC" : sub ? " " * full[i][7:8] * "ˢ" : " " * full[i][4:5] * "ᵐ"
            for i in eachindex(full)]
end

# lon/lat polylines (NaN-separated) of any line geometry
function _xy!(xs, ys, g)
    if g isa GB.LineString
        for p in GB.coordinates(g)
            push!(xs, p[1]); push!(ys, p[2])
        end
        push!(xs, NaN); push!(ys, NaN)
    elseif g isa GB.MultiLineString
        foreach(l -> _xy!(xs, ys, l), g.linestrings)
    elseif g isa AbstractVector
        foreach(l -> _xy!(xs, ys, l), g)
    end
    return xs, ys
end
_xy(g) = _xy!(Float64[], Float64[], g)

# break a polyline where it jumps across the antimeridian or has gaps
function _segments(lon, lat)
    xs = Float64[]; ys = Float64[]
    prev = NaN
    for (λ, φ) in zip(lon, lat)
        if !isnan(prev) && !isnan(λ) && abs(λ - prev) > 180
            push!(xs, NaN); push!(ys, NaN)
        end
        push!(xs, λ); push!(ys, φ)
        prev = λ
    end
    return xs, ys
end

# Sutherland–Hodgman clipping of a polygon ring to a lon/lat box
function _clip_ring(pts::Vector{NTuple{2,Float64}}, x1, x2, y1, y2)
    for (inside, cut) in ((p -> p[1] >= x1, (a, b) -> (x1, a[2] + (b[2] - a[2]) * (x1 - a[1]) / (b[1] - a[1]))),
                          (p -> p[1] <= x2, (a, b) -> (x2, a[2] + (b[2] - a[2]) * (x2 - a[1]) / (b[1] - a[1]))),
                          (p -> p[2] >= y1, (a, b) -> (a[1] + (b[1] - a[1]) * (y1 - a[2]) / (b[2] - a[2]), y1)),
                          (p -> p[2] <= y2, (a, b) -> (a[1] + (b[1] - a[1]) * (y2 - a[2]) / (b[2] - a[2]), y2)))
        isempty(pts) && return pts
        out = NTuple{2,Float64}[]
        prev = pts[end]
        for cur in pts
            if inside(cur)
                inside(prev) || push!(out, cut(prev, cur))
                push!(out, cur)
            elseif inside(prev)
                push!(out, cut(prev, cur))
            end
            prev = cur
        end
        pts = out
    end
    return pts
end

_polygons(g::GB.Polygon) = [g]
_polygons(g::GB.MultiPolygon) = collect(g.polygons)
_polygons(g::AbstractVector) = reduce(vcat, [_polygons(x) for x in g]; init=GB.Polygon[])

# land clipped to a box well beyond the map limits (polygons that reach around the globe,
# such as Antarctica, are otherwise mangled by non-cylindrical projections)
function _clipped_land(land, limits)
    x1, x2, y1, y2 = Float64.(limits)
    w, h = x2 - x1, y2 - y1
    x1, x2, y1, y2 = x1 - w, x2 + w, max(y1 - h, -90.0), min(y2 + h, 90.0)
    out = GB.Polygon{2,Float64}[]
    for p in _polygons(land)
        ring = [(Float64(q[1]), Float64(q[2])) for q in GB.coordinates(p.exterior)]
        c = _clip_ring(ring, x1, x2, y1, y2)
        length(c) >= 3 && push!(out, GB.Polygon([GB.Point2(x...) for x in c]))
    end
    return out
end

# land polygons clipped to the visible hemisphere of a globe centred on (λ0, φ0): clipped on
# the sphere against the plane of a small circle `margin` degrees inside the horizon, with the
# horizon arc inserted where a polygon is cut; longitudes unwrapped around λ0
_unit(λ, φ) = [cosd(φ) * cosd(λ), cosd(φ) * sind(λ), sind(φ)]
_proj_centre(dest) = ((x = match(r"\+lon_0=(-?[\d.]+)", dest); x === nothing ? 0.0 : parse(Float64, x[1])),
                      (x = match(r"\+lat_0=(-?[\d.]+)", dest); x === nothing ? 0.0 : parse(Float64, x[1])))
function _hemisphere_land(land, λ0, φ0; margin=0.3)
    c = _unit(λ0, φ0)
    h = sind(margin)                                          # the clip plane: x ⋅ c = h
    function onbound(x)
        xp = x .- (x' * c) .* c
        return h .* c .+ sqrt(1 - h^2) .* xp ./ sqrt(sum(abs2, xp))
    end
    function arc(p, q)                                        # along the clip circle from p to q
        u = p .- (p' * c) .* c; u ./= sqrt(sum(abs2, u))
        w = q .- (q' * c) .* c; w ./= sqrt(sum(abs2, w))
        v = [c[2] * u[3] - c[3] * u[2], c[3] * u[1] - c[1] * u[3], c[1] * u[2] - c[2] * u[1]]
        α = atan(w' * v, w' * u)                              # signed angle u → w about c
        m = max(2, ceil(Int, abs(α) / deg2rad(2)))
        return [h .* c .+ sqrt(1 - h^2) .* (cos(t) .* u .+ sin(t) .* v) for t in range(0, α; length=m)]
    end
    lonlat(x) = (λ = atand(x[2], x[1]); (λ0 + mod(λ - λ0 + 180, 360) - 180, asind(clamp(x[3], -1, 1))))
    out = GB.Polygon{2,Float64}[]
    for p in _polygons(land)
        ring = [_unit(Float64(q[1]), Float64(q[2])) for q in GB.coordinates(p.exterior)]
        length(ring) >= 3 && ring[end] ≈ ring[1] && pop!(ring)
        n = length(ring)
        n < 3 && continue
        i0 = findfirst(x -> x' * c >= h, ring)
        i0 === nothing && continue                            # entirely on the far side
        ring = circshift(ring, 1 - i0)                        # start at a visible vertex
        pts = Vector{Vector{Float64}}()
        exitpt = nothing
        for i in 1:n
            a, b = ring[i], ring[mod1(i + 1, n)]
            da, db = a' * c - h, b' * c - h
            da >= 0 && push!(pts, a)
            if (da >= 0) != (db >= 0)
                t = da / (da - db)
                x = onbound((1 - t) .* a .+ t .* b)
                if da >= 0                                    # leaving the visible side
                    push!(pts, x); exitpt = x
                else                                          # re-entering: follow the horizon
                    exitpt === nothing || append!(pts, arc(exitpt, x)[2:end-1])
                    push!(pts, x); exitpt = nothing
                end
            end
        end
        length(pts) >= 3 && push!(out, GB.Polygon([GB.Point2(lonlat(x)...) for x in pts]))
    end
    return out
end

_inside(lon, lat, limits) = limits === nothing ||
    (limits[1] <= lon <= limits[2] && limits[3] <= lat <= limits[4])

# country capitals (Natural Earth populated places) in the map region
function _capitals(scale, limits)
    fc = _layer("populated_places_simple", scale == 110 ? 110 : 50)
    fc === nothing && return (Float64[], Float64[], String[])
    lon = Float64[]; lat = Float64[]; name = String[]
    for f in fc
        p = GeoMakie.GeoJSON.properties(f)
        get(p, :featurecla, "") == "Admin-0 capital" || continue
        λ, φ = Float64(p[:longitude]), Float64(p[:latitude])
        _inside(λ, φ, limits) || continue
        push!(lon, λ); push!(lat, φ); push!(name, string(p[:name]))
    end
    return lon, lat, name
end

# observers at (nearly) the same place share one marker: groups of indices
function _observer_groups(obs; tol=0.02)
    groups = Vector{Vector{Int}}()
    for (i, o) in enumerate(obs)
        g = findfirst(g -> abs(obs[g[1]].lon - o.lon) < tol && abs(obs[g[1]].lat - o.lat) < tol, groups)
        g === nothing ? push!(groups, [i]) : push!(groups[g], i)
    end
    return groups
end

# the shadow band between the body limits as a polygon, from the lines across the shadow
# (`m.band`): at each epoch its outermost points on the ground, so that the band follows the
# Earth's limb where only part of the shadow falls on the Earth; on the globe (`centre` = the
# projection centre) only its part on the visible hemisphere, where the projection is defined
function _band(m; centre=nothing)
    q = m.band
    q === nothing && return nothing
    vis(λ, φ) = isfinite(λ) && isfinite(φ) && (centre === nothing ||
                sind(centre[2]) * sind(φ) + cosd(centre[2]) * cosd(φ) * cosd(λ - centre[1]) > cosd(89.9))
    left = GB.Point2{Float64}[]; right = GB.Point2{Float64}[]
    for i in eachindex(q.times)
        ks = [k for k in eachindex(q.offsets) if vis(q.lon[i, k], q.lat[i, k])]
        isempty(ks) && continue
        push!(left, GB.Point2(q.lon[i, ks[1]], q.lat[i, ks[1]]))
        push!(right, GB.Point2(q.lon[i, ks[end]], q.lat[i, ks[end]]))
    end
    length(left) < 2 && return nothing
    (maximum(abs, diff([p[1] for p in left])) > 180 || maximum(abs, diff([p[1] for p in right])) > 180) && return nothing
    return GB.Polygon([left; reverse(right)])
end

# GeoMakie labels the grid lines with 3 significant digits of the inverse-projected position,
# which turns 0° into e.g. "1.91e-6°": such labels are rewritten as "0°"
function _round_ticklabels!(ga)
    tiny = r"^-?[0-9.]+e-[0-9]+°$"
    for p in ga.blockscene.plots
        p isa Makie.Text || continue
        Makie.on(p.text; update=true) do t
            t isa AbstractVector{<:AbstractString} && any(x -> occursin(tiny, x), t) || return
            p.text[] = [occursin(tiny, x) ? "0°" : x for x in t]
        end
    end
    return ga
end

# grid spacing [deg] for a map extent of `span` degrees: at most about six lines
_grid_step(span) = (for st in (1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0)
                        span / st <= 6 && return st
                    end; 30.0)

"""
    geomap(map::OccMap; projection=nothing, size=(1000, 850), resolution=nothing, land=nothing,
           coastlines=true, borders=true, capitals=nothing, band=true, limits=nothing,
           observer_labels=:number, terminator=false, title=map.title, titlesize=20,
           ticklabels=nothing) -> Figure
    geomap(event_or_fit_or_prediction; <same keywords>, kwargs...)

Map of an occultation path with GeoMakie.

* `projection` – any PROJ string (default: orthographic view centred on the sub-stellar
  point for `projection=:ortho` maps of [`occ_map`](@ref), Equal Earth otherwise), e.g.
  `"+proj=laea +lon_0=10 +lat_0=45"`; `limits = (lonmin, lonmax, latmin, latmax)` zooms in.
* Natural Earth background at `resolution` 110, 50 or 10 (1:`resolution` million; default
  50 for regional maps with `limits`, else 110): `land` (default off for orthographic views,
  where PROJ does not clip polygons at the horizon), `coastlines`, country `borders` and
  country `capitals` (default on for regional maps). The 1:110m coastlines and land are
  bundled with GeoMakie; other layers are downloaded once (Natural Earth, public domain).
* `band` – shade the shadow band between the body limits.
* `observer_labels` – `:number` (default): stations are numbered on the map and listed with
  their names in a legend (stations closer than ~2 km share a marker), which keeps crowded
  maps readable; `:name`: names written next to the markers; `:none`. Positive stations are
  green circles, negative ones red crosses.

* `title` in `titlesize` points, above the frame; `ticklabels` – the longitude and latitude
  labels of the grid (default: on regional maps, off on the globe, where they crowd the limb).

* `overlays`: further maps drawn as lines in `overlay_colors`, labelled `overlay_labels`
  (e.g. a prediction; see `geomap(prediction, fit)`); `path_label` prefixes the legend
  entries of the main path; `terminator = true` draws the day–night terminator.

Other keywords are passed to [`occ_map`](@ref) (`window`, `step`, `error`, `marks`, ...).
Save with `save("map.png", fig)`.
"""
function geomap(m::OccMap; projection=nothing, size=(1000, 850), resolution=nothing, land=nothing,
                coastlines::Bool=true, borders::Bool=true, capitals=nothing, band::Bool=true, limits=nothing,
                observer_labels::Symbol=:number, terminator=nothing, night=nothing, title::AbstractString=m.title,
                titlesize::Real=20, ticklabels=nothing,
                path_color=:blue,
                limit_color=:black, error_color=:red, band_color=(:steelblue, 0.18), path_label::AbstractString="",
                overlays=OccMap[], overlay_labels=["prediction"], overlay_colors=[:darkorange, :purple, :teal])
    observer_labels in (:number, :name, :none) || throw(ArgumentError("observer_labels must be :number, :name or :none"))
    dest = projection === nothing ? _default_projection(m) : String(projection)
    scale = something(resolution, limits === nothing ? 110 : 50)
    scale in (110, 50, 10) || throw(ArgumentError("resolution must be 110, 50 or 10"))
    fig = Figure(; size=size)
    globe = occursin("ortho", dest)
    # grid labels on regional maps only (on the globe they crowd the limb); the title well above the frame
    tl = something(ticklabels, !globe)
    # grid lines: a spacing of 1, 2, 5, 10, 15, 20 or 30° with at most about six lines across a
    # regional map; integer labels
    gx, gy = limits === nothing ? (30.0, 15.0) : (_grid_step(limits[2] - limits[1]), _grid_step(limits[4] - limits[3]))
    xt = limits === nothing ? (-180.0:30.0:180.0) : (floor(limits[1] / gx) * gx):gx:(ceil(limits[2] / gx) * gx)
    yt = limits === nothing ? (-90.0:15.0:90.0) : (floor(limits[3] / gy) * gy):gy:(ceil(limits[4] / gy) * gy)
    ga = GeoAxis(fig[1, 1]; dest=dest, title=title, titlesize=titlesize, titlegap=18,
                 xticklabelsvisible=tl, yticklabelsvisible=tl, xticks=collect(xt), yticks=collect(yt))
    _round_ticklabels!(ga)
    # a frame along the drawing area of the map (GeoAxis has no spines)
    frame = lift(ga.scene.viewport) do r
        x0, y0 = r.origin; w, h = r.widths
        Point2f[(x0, y0), (x0 + w, y0), (x0 + w, y0 + h), (x0, y0 + h), (x0, y0)]
    end
    lines!(ga.blockscene, frame; color=:black, linewidth=1.5, space=:pixel)
    # background
    if something(land, true)
        L = _layer("land", scale)
        L === nothing && scale != 110 && (L = _layer("land", 110))
        if L !== nothing
            polys = occursin("ortho", dest) ? _hemisphere_land(L, _proj_centre(dest)...) :
                    limits === nothing ? L : _clipped_land(L, limits)
            isempty(polys) || poly!(ga, polys; color=(:tan, 0.35), strokecolor=:transparent)
        end
    end
    if borders
        B = _layer("admin_0_boundary_lines_land", scale)
        B === nothing || lines!(ga, _xy(GeoMakie.geo2basic(B))...; color=:gray45, linewidth=0.6, linestyle=:dash)
    end
    if coastlines
        C = _layer("coastline", scale)
        C === nothing && scale != 110 && (C = _layer("coastline", 110))
        C === nothing || lines!(ga, _xy(C)...; color=:gray25, linewidth=0.6)
    end
    if something(capitals, limits !== nothing)
        lon, lat, name = _capitals(scale, limits)
        if !isempty(lon)
            scatter!(ga, lon, lat; marker=:rect, markersize=5, color=:gray20)
            text!(ga, lon, lat; text=[" " * n for n in name], fontsize=8, color=:gray30, align=(:left, :center))
        end
    end
    # the night side (Sun below the horizon), shaded transparently with the twilight as a
    # gradient (the Sun 0–18° below the horizon); by default on the globe
    if something(night, globe)
        t = m.path.times[cld(length(m.path.times), 2)]
        sv = -Occultations.earth_heliocentric(t)
        st = Occultations.celestial_to_terrestrial(t) * (sv ./ sqrt(sum(abs2, sv)))
        lons = range(-180, 180; length=361); lats = range(-90, 90; length=181)
        α = [clamp(-asind(clamp(cosd(φ) * cosd(λ) * st[1] + cosd(φ) * sind(λ) * st[2] + sind(φ) * st[3], -1, 1)) / 18, 0, 1)
             for λ in lons, φ in lats]
        meshimage!(ga, -180 .. 180, -90 .. 90, [RGBAf(0.02, 0.02, 0.15, 0.38a) for a in α]; npoints=400)
    end
    something(terminator, globe) && lines!(ga, _segments(m.terminator...)...; color=:navy, linestyle=:dashdot,
                                           linewidth=1, label="terminator")
    # the limb of the Earth on the globe: the horizon circle of the projection centre
    if globe
        λ0 = (x = match(r"\+lon_0=(-?[\d.]+)", dest); x === nothing ? 0.0 : parse(Float64, x[1]))
        φ0 = (x = match(r"\+lat_0=(-?[\d.]+)", dest); x === nothing ? 0.0 : parse(Float64, x[1]))
        d = 89.9
        az = range(0, 360; length=721)
        φs = [asind(sind(φ0) * cosd(d) + cosd(φ0) * sind(d) * cosd(a)) for a in az]
        λs = [λ0 + atand(sind(a) * sind(d) * cosd(φ0), cosd(d) - sind(φ0) * sind(φ)) for (a, φ) in zip(az, φs)]
        lines!(ga, λs, φs; color=:black, linewidth=1.8)
    end
    # shadow
    p = m.path
    if band
        poly_ = _band(m; centre=globe ? _proj_centre(dest) : nothing)
        poly_ === nothing || poly!(ga, poly_; color=band_color, strokecolor=:transparent, label="shadow")
    end
    kinds = Occultations._offset_kinds(p.offsets)
    seen = Set{Symbol}()
    for (k, δ) in enumerate(p.offsets)
        col, ls, lab = kinds[k] === :centre ? (path_color, :solid, "centre line") :
                       kinds[k] === :limit ? (limit_color, :solid, "body limits") : (error_color, :dash, "1σ uncertainty")
        lines!(ga, _segments(p.lon[:, k], p.lat[:, k])...; color=col, linestyle=ls, linewidth=δ == 0 ? 2 : 1.2,
               label=(kinds[k] in seen ? nothing : (isempty(path_label) ? lab : path_label * ": " * lab)))
        push!(seen, kinds[k])
    end
    # further paths (e.g. the prediction), as lines in their own colour
    for (j, om) in enumerate(overlays)
        q = om.path
        kq = Occultations._offset_kinds(q.offsets)
        col = overlay_colors[mod1(j, length(overlay_colors))]
        pre = j <= length(overlay_labels) ? overlay_labels[j] : "path $j"
        seenq = Set{Symbol}()
        for (k, δ) in enumerate(q.offsets)
            ls, lw, lab = kq[k] === :centre ? (:solid, 2.0, "centre line") :
                          kq[k] === :limit ? (:solid, 1.2, "body limits") : (:dot, 1.4, "1σ uncertainty")
            lines!(ga, _segments(q.lon[:, k], q.lat[:, k])...; color=col, linestyle=ls, linewidth=lw,
                   label=(kq[k] in seenq ? nothing : pre * ": " * lab))
            push!(seenq, kq[k])
        end
    end
    k0 = argmin(abs.(p.offsets))
    if !isempty(m.marks)
        ok = [i for i in m.marks if !isnan(p.lon[i, k0]) && _inside(p.lon[i, k0], p.lat[i, k0], limits)]
        scatter!(ga, p.lon[ok, k0], p.lat[ok, k0]; color=path_color, markersize=6)
        text!(ga, p.lon[ok, k0], p.lat[ok, k0]; text=_mark_labels(p.times[ok]), fontsize=9, align=(:left, :center),
              color=path_color)
    end
    # the direction of motion: an arrow along the centre line
    vis = [i for i in eachindex(p.times) if !isnan(p.lon[i, k0]) && _inside(p.lon[i, k0], p.lat[i, k0], limits)]
    if length(vis) >= 8
        n = length(vis)
        i1, i2 = vis[clamp(round(Int, 0.62n), 1, n)], vis[clamp(round(Int, 0.62n) + max(2, n ÷ 20), 1, n)]
        f = ga.transform_func[]
        P1 = Makie.apply_transform(f, Point2f(p.lon[i1, k0], p.lat[i1, k0]))
        P2 = Makie.apply_transform(f, Point2f(p.lon[i2, k0], p.lat[i2, k0]))
        if all(isfinite, P1) && all(isfinite, P2) && P1 != P2
            lines!(ga, p.lon[i1:i2, k0], p.lat[i1:i2, k0]; color=path_color, linewidth=2.5)
            scatter!(ga, [p.lon[i2, k0]], [p.lat[i2, k0]]; marker=:rtriangle, markersize=12, color=path_color,
                     rotation=atan(P2[2] - P1[2], P2[1] - P1[1]),
                     label="direction of motion" => (; markersize=8))
        end
    end
    # the point of the centre line where the star stands highest (the zenith, if reached)
    zl, zb, za = m.zenith
    if isfinite(zl) && _inside(zl, zb, limits)
        scatter!(ga, [zl], [zb]; marker=:star5, markersize=16, color=:gold, strokecolor=:black, strokewidth=0.8,
                 label=(za > 89.5 ? "star at the zenith" : @sprintf("star highest (%.0f°)", za)) => (; markersize=10))
    end
    # observers
    colour(s) = s === :positive ? :forestgreen : s === :negative ? :firebrick : :darkorange
    mk(s) = s === :negative ? :xcross : :circle
    obs_entries = Any[]; obs_labels = String[]
    if !isempty(m.observers)
        for g in _observer_groups(m.observers)
            o = m.observers[g[1]]
            st = m.observer_status[g[1]]
            scatter!(ga, [o.lon], [o.lat]; marker=mk(st), color=colour(st), markersize=st === :negative ? 11 : 9,
                     strokecolor=:white, strokewidth=0.5)
            if observer_labels === :number
                text!(ga, [o.lon], [o.lat]; text=" " * join(g, ","), fontsize=10, font=:bold, color=colour(st),
                      align=(:left, :bottom))
            elseif observer_labels === :name
                text!(ga, [o.lon], [o.lat]; text=" " * join((m.observers[i].name for i in g), ", "), fontsize=9,
                      color=colour(st), align=(:left, :bottom))
            end
        end
        if observer_labels === :number
            for (i, o) in enumerate(m.observers)
                st = m.observer_status[i]
                push!(obs_entries, MarkerElement(marker=mk(st), color=colour(st), markersize=9))
                push!(obs_labels, "$i  $(o.name)")
            end
        end
    end
    if limits !== nothing
        GeoMakie.xlims!(ga, limits[1], limits[2]); GeoMakie.ylims!(ga, limits[3], limits[4])
    end
    # a compact legend (on the globe it sits in the lower-left corner, outside most of the disc)
    axislegend(ga; position=:lb, merge=true, framevisible=true, labelsize=globe ? 10 : 10,
               patchsize=globe ? (16, 9) : (20, 10), rowgap=globe ? 1 : 2, padding=globe ? (5, 5, 4, 4) : (6, 6, 6, 6),
               margin=globe ? (4, 4, 4, 4) : (10, 10, 10, 10), backgroundcolor=(:white, 0.85))
    if !isempty(obs_entries)
        Legend(fig[1, 2], obs_entries, obs_labels, "Stations"; labelsize=10, titlesize=11, framevisible=false,
               rowgap=0, patchsize=(10, 10))
    end
    return fig
end

"""
    geomap(prediction::OccMap, fit::OccMap; kwargs...) -> Figure

The post-fit shadow path (`fit`, e.g. `occ_map(limbfit)`) with the predicted one
(`prediction`, e.g. `occ_map(event; error_mas=20)`) overlaid in another colour.
"""
geomap(pred::OccMap, fit::OccMap; kw...) =
    geomap(fit; overlays=[pred], overlay_labels=["prediction"], path_label="post-fit", kw...)

function geomap(x::Union{Occultations.Event,LimbFit}; projection=nothing, size=(1000, 850), resolution=nothing,
                land=nothing, coastlines::Bool=true, borders::Bool=true, capitals=nothing, band::Bool=true,
                limits=nothing, observer_labels::Symbol=:number, kw...)
    return geomap(occ_map(x; kw...); projection=projection, size=size, resolution=resolution, land=land,
                  coastlines=coastlines, borders=borders, capitals=capitals, band=band, limits=limits,
                  observer_labels=observer_labels)
end

function geomap(p::Occultations.OccultationPrediction, body; projection=nothing, size=(1000, 850), resolution=nothing,
                land=nothing, coastlines::Bool=true, borders::Bool=true, capitals=nothing, band::Bool=true,
                limits=nothing, observer_labels::Symbol=:number, kw...)
    return geomap(occ_map(p, body; kw...); projection=projection, size=size, resolution=resolution,
                  land=land, coastlines=coastlines, borders=borders, capitals=capitals, band=band, limits=limits,
                  observer_labels=observer_labels)
end

end # module
