# ---------------------------------------------------------------------------
# Maps of occultation paths (SORA `plot_occ_map`): the data; drawn by `geomap` (GeoMakie)
# ---------------------------------------------------------------------------

"""
    OccMap

Data for a map of an occultation path, produced by [`occ_map`](@ref) and drawn
with [`geomap`](@ref).
"""
struct OccMap
    path::ShadowPath
    center::Tuple{Float64,Float64}
    projection::Symbol
    zoom::Float64
    observers::Vector{Observer}
    observer_status::Vector{Symbol}       # :positive, :negative or :unknown
    marks::Vector{Int}
    terminator::Tuple{Vector{Float64},Vector{Float64}}
    grid::Bool
    title::String
    zenith::Tuple{Float64,Float64,Float64}   # (lon, lat, altitude) where the star is highest on the centre line
    band::Union{Nothing,ShadowPath}          # lines across the shadow between the limits (for the band)
end

Base.show(io::IO, m::OccMap) = print(io, "OccMap(", m.path.body, ", ", length(m.path.times), " epochs, ",
                                     m.projection, " centred at ", round.(m.center; digits=2), ")")

"""
    substellar_point(star, t; dut1=0) -> (lon, lat) [deg]

Geographic point where `star` is at the zenith at `t`.
"""
function substellar_point(star::Star, t; dut1::Real=0.0)
    t = _epoch(t)
    v = celestial_to_terrestrial(t; dut1=dut1) * radec2vec(geocentric(star, t)...)
    return atand(v[2], v[1]), asind(clamp(v[3], -1, 1))
end

# UTC seconds of the day of an epoch
_utc_sod(e::Epoch) = (s = isot(e; digits=4); parse(Int, s[12:13]) * 3600 + parse(Int, s[15:16]) * 60 + parse(Float64, s[18:end]))

# great circle of the terminator (Sun altitude `alt`) at t, as lon/lat [deg]
function _terminator(t::Epoch; alt::Real=0.0, n::Integer=361)
    s = -earth_heliocentric(t)
    st = celestial_to_terrestrial(t) * (s ./ norm(s))
    e, nn = sky_basis(st)
    lon = Float64[]; lat = Float64[]
    for θ in range(0, 2π; length=n)
        v = sind(alt) .* st .+ cosd(alt) .* (cos(θ) .* e .+ sin(θ) .* nn)
        push!(lon, atand(v[2], v[1])); push!(lat, asind(clamp(v[3], -1, 1)))
    end
    return lon, lat
end

"""
    occ_map(star, body, t; window=600, step=5, radius=body.radius, offsets=nothing,
            error=0, observers=Observer[], observer_status=nothing, center=nothing, projection=:ortho, zoom=1,
            marks=60, title="") -> OccMap
    occ_map(prediction::OccultationPrediction, body; kwargs...)
    occ_map(event::Event; error_mas=nothing, kwargs...)
    occ_map(fit; event=1, ndraw=400, kwargs...)

Map of the shadow path between `t - window` and `t + window` (computed every `step`
seconds with [`shadow_path`](@ref), i.e. on the rotating Earth): central line, body
limits (`±radius`) and optionally error lines (`±(radius + error)` km), time marks every
`marks` seconds of UTC (whole minutes for `marks = 60`) along the central line (for an
event, `error_mas` gives the cross-track
ephemeris uncertainty of a prediction in mas), observers and the day/night terminator for
`t`. Draw it with [`geomap`](@ref).

`projection` (`:ortho`, orthographic and centred by default on the sub-stellar point at `t`,
like SORA, or `:lonlat`) and `center` (`(lon, lat)` in degrees) set the default projection
of `geomap`. For an [`Event`](@ref) the chord observers are shown; for a limb
fit the path is that of the fitted centre (ephemeris + fitted offset); the body limits are
the extent of the fitted limb across the motion, and the error lines their ±1σ from the fit
covariance (centre and shape; `ndraw` draws).
"""
function occ_map(star::Star, body, t; window::Real=600.0, step::Real=5.0, radius::Real=_radius(body),
                 offsets=nothing, error::Real=0.0, observers=Observer[], observer_status=nothing, center=nothing,
                 projection::Symbol=:ortho, zoom::Real=1.0, marks::Real=60.0, grid::Bool=true,
                 title::AbstractString="", dut1::Real=0.0)
    projection in (:ortho, :lonlat) || throw(ArgumentError("projection must be :ortho or :lonlat"))
    t = _epoch(t)
    offs = if offsets !== nothing
        Float64.(collect(offsets))
    else
        o = radius > 0 ? [0.0, -radius, radius] : [0.0]
        error > 0 && append!(o, [-(radius + error), radius + error])
        o
    end
    # the path sampled on whole multiples of `step` in UTC, so that the time marks fall on them
    t0 = t - window
    t0 = t0 - mod(_utc_sod(t0), step)
    path = shadow_path(star, body, t0, t + window; step=step, offsets=offs, dut1=dut1)
    c = center === nothing ? substellar_point(star, t; dut1=dut1) : (Float64(center[1]), Float64(center[2]))
    sod = _utc_sod.(path.times)
    mk = marks > 0 ? [i for i in eachindex(sod) if abs(mod(sod[i] + marks / 2, marks) - marks / 2) < step / 2] : Int[]
    # where the star stands highest on the centre line (at the zenith if the line reaches it)
    k0 = argmin(abs.(path.offsets))
    zen = (NaN, NaN, -90.0)
    for (i, s) in enumerate(path.times)
        λ, φ = path.lon[i, k0], path.lat[i, k0]
        isnan(λ) && continue
        ra, dec = geocentric(star, s)
        alt = altaz(Observer("", λ, φ, 0.0), s, ra, dec; dut1=dut1)[1]
        alt > zen[3] && (zen = (λ, φ, alt))
    end
    ttl = isempty(title) ? "$(path.body) occults $(_label(star)) – $(isot(t; digits=0)) UTC" : String(title)
    obs = collect(Observer, observers)
    st = observer_status === nothing ? fill(:unknown, length(obs)) : collect(Symbol, observer_status)
    length(st) == length(obs) || throw(DimensionMismatch("observer_status and observers differ in length"))
    # the shadow sampled across, between the body limits: where only part of the shadow falls on
    # the Earth, the band follows the Earth's limb
    lim = [path.offsets[k] for k in findall(==(:limit), _offset_kinds(path.offsets))]
    band = length(lim) == 2 ? shadow_path(star, body, t0, t + window; step=step,
                                          offsets=range(minimum(lim), maximum(lim); length=17), dut1=dut1) : nothing
    return OccMap(path, c, projection, Float64(zoom), obs, st, mk, _terminator(t), grid, ttl, zen, band)
end

occ_map(p::OccultationPrediction, body; kw...) = occ_map(p.star, body, p.epoch; kw...)

_map_chords(ev::Event) = [c for c in ev.chords if c.observer isa Observer]
_status(c::Chord) = isnegative(c) ? :negative : :positive

function occ_map(ev::Event; observers=nothing, error_mas=nothing, kw...)
    cs = _map_chords(ev)
    obs, st = observers === nothing ? ([c.observer for c in cs], _status.(cs)) : (observers, nothing)
    # a cross-track ephemeris uncertainty in mas → km at the body
    err = error_mas === nothing ? () : (error=error_mas / RAD2MAS * ev.ca.dist_km,)
    return occ_map(ev.star, ev.body, ev.tref; observers=obs, observer_status=st, err..., kw...)
end

function occ_map(r::LimbFit; event=1, ndraw::Integer=400, kw...)
    p = r.problem
    k = _event_index(p, event)
    ev = p.setup.events[k]
    a = astrometry(r; event=k)
    b = deepcopy(ev.body)
    add_offset!(b.ephem, a.offset_ra, a.offset_dec)
    (c0, O0) = outlines(r)[k]
    # the shadow limits: the extent of the limb across the motion (left-hand normal n, as in
    # `shadow_path`), measured from the fitted centre; their 1σ from draws of the fit
    # covariance (centre and shape together)
    pv = positionv(ev.star, _ephem(ev.body), GEOCENTER, ev.tref)
    n = (-pv.vg, pv.vf) ./ hypot(pv.vf, pv.vg)
    edges(c, O) = (x = outline_points(O; center=c, n=361); s = filter(isfinite, x[1] .* n[1] .+ x[2] .* n[2]);
                   (minimum(s), maximum(s)))
    s0 = c0[1] * n[1] + c0[2] * n[2]
    lo0, hi0 = edges(c0, O0) .- s0
    offs = [0.0, lo0, hi0]
    F = _Free(parameters(p))
    L = cholesky(Symmetric(r.covariance + 1e-12I); check=false)
    if all(isfinite, r.covariance) && issuccess(L)
        rng = MersenneTwister(1)
        E = Tuple{Float64,Float64}[]
        for _ in 1:ndraw
            v = copy(r.values)
            v[F.idx] .+= L.L * randn(rng, length(F.idx))
            p isa LimbProblem && !_satisfies_constraints(p, v) && continue
            try
                (c, O) = outlines(p, v)[k]
                push!(E, edges(c, O) .- s0)
            catch
            end
        end
        if length(E) >= 20
            σlo, σhi = std(first.(E)), std(last.(E))
            append!(offs, [lo0 - σlo, hi0 + σhi])
        end
    end
    cs = _map_chords(ev)
    return occ_map(ev.star, b, ev.tref; offsets=get(Dict(kw), :offsets, offs), observers=[c.observer for c in cs],
                   observer_status=_status.(cs), (x for x in kw if first(x) !== :offsets)...)
end

# kind of each cross-track offset line: :centre (0), :limit (the nearest line on each side of
# the centre line), :error (farther lines, e.g. ±1σ of the limits)
function _offset_kinds(offs)
    kinds = fill(:error, length(offs))
    for (k, δ) in enumerate(offs)
        δ == 0 && (kinds[k] = :centre)
    end
    for side in (-1, 1)
        ks = [k for (k, δ) in enumerate(offs) if sign(δ) == side]
        isempty(ks) || (kinds[ks[argmin(abs.(offs[ks]))]] = :limit)
    end
    return kinds
end

