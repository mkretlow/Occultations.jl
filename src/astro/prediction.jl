"""
    OccultationPrediction

One predicted event (geocentric closest approach).

Fields: `epoch`, `star`, `body`, `ca` [arcsec], `ca_km`, `pa` [deg, body relative
to star], `vel` [km/s], `dist` [au], `mag` (G or first available),
`lon`, `lat` [deg] of the shadow centre at closest approach (NaN when it
misses the Earth), `sun_alt` there.
"""
struct OccultationPrediction
    epoch::Epoch
    star::Star
    body::String
    ca::Float64
    ca_km::Float64
    pa::Float64
    vel::Float64
    dist::Float64
    mag::Float64
    lon::Float64
    lat::Float64
    sun_alt::Float64
end

_star_mag(s::Star) = haskey(s.mag, "G") ? s.mag["G"] : isempty(s.mag) ? NaN : first(values(s.mag))

function Base.show(io::IO, p::OccultationPrediction)
    @printf(io, "%s  %s %s  CA=%7.4f\"  PA=%6.2f  vel=%6.2f km/s  D=%7.3f au  mag=%5.2f  lon=%8.3f lat=%7.3f",
            isot(p.epoch; digits=3), format_ra(p.star.ra; digits=4), format_dec(p.star.dec; digits=3),
            p.ca, p.pa, p.vel, p.dist, p.mag, p.lon, p.lat)
end

function Base.show(io::IO, ::MIME"text/plain", ps::AbstractVector{OccultationPrediction})
    println(io, length(ps), " predicted occultation(s)")
    println(io, " epoch (UTC)                star RA (ICRS)  Dec               C/A [\"]    P/A    vel        dist      mag   shadow centre")
    for p in ps
        println(io, " ", p)
    end
end

"""
    predict(star_or_stars, body, tstart, tend; step=60, radius=body.radius,
            margin=0.0, observer=GEOCENTER) -> Vector{OccultationPrediction}

Search for occultations of one or several stars by `body` between `tstart`
and `tend`. An event is reported when the geocentric closest approach in the
fundamental plane is smaller than `R_earth + radius + margin` [km].
"""
function predict(stars::AbstractVector{Star}, body, tstart, tend; step::Real=60.0,
                 radius::Real=_radius(body), margin::Real=0.0)
    eph = _ephem(body)
    name = body isa Body ? body.name : eph.name
    t0, t1 = _epoch(tstart), _epoch(tend)
    span = t1 - t0
    n = max(3, ceil(Int, span / step) + 1)
    h = span / (n - 1)
    times = [t0 + (i - 1) * h for i in 1:n]
    B = [geocentric_vector(eph, t) for t in times]
    D = norm.(B)
    U = [b ./ d for (b, d) in zip(B, D)]
    limit = R_EARTH_KM + radius + margin
    out = OccultationPrediction[]
    # star directions sampled every `stride` points and linearly interpolated
    stride = max(1, round(Int, 3600 / h))
    knots = unique([1:stride:n; n])
    for star in stars
        uk = [radec2vec(geocentric(star, times[k])...) for k in knots]
        d = Vector{Float64}(undef, n)
        for i in 1:n
            j = min((i - 1) ÷ stride + 1, length(knots) - 1)
            w = (i - knots[j]) / (knots[j+1] - knots[j])
            us = (1 - w) .* uk[j] .+ w .* uk[j+1]
            c = dot(us ./ norm(us), U[i])
            d[i] = c > 0 ? D[i] * sqrt(max(1 - c^2, 0.0)) : Inf
        end
        for i in 2:n-1
            # 30 km/s: generous bound on the shadow speed between grid points
            (d[i] <= d[i-1] && d[i] < d[i+1] && d[i] < limit + 30 * h) || continue
            ca = closest_approach(star, eph, times[i-1], times[i+1]; step=h / 4)
            ca.ca_km < limit || continue
            push!(out, _make_prediction(star, eph, name, ca))
        end
    end
    sort!(out; by=p -> p.epoch)
    return out
end
predict(star::Star, body, tstart, tend; kw...) = predict([star], body, tstart, tend; kw...)

function _make_prediction(star, eph, name, ca)
    path = shadow_path(star, eph, ca.epoch, ca.epoch; offsets=[0.0])
    return OccultationPrediction(ca.epoch, star, name, ca.ca_arcsec, ca.ca_km, ca.pa, ca.vel,
                                 ca.dist_km / AU_KM, _star_mag(star), path.lon[1, 1], path.lat[1, 1],
                                 path.sun_alt[1, 1])
end
