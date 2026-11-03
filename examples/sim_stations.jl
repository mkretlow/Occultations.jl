# Stations of the synthetic examples (included by the sim_*.jl scripts): placed on land, at
# night, with the star above 20°, spread along the track, each with its own equipment.
#
#   stations(star, body, tref, offsets; rng)   -> Vector of (observer, τ, equipment)
#   random_offsets(rng, n, lo, hi; mingap)     -> n random cross-track offsets [km]
#
# For every cross-track offset (km, as in `shadow_path`) the station sits on the ground track of
# that offset at a random instant τ (s from `tref`) within the longest stretch of the path over
# land and in the dark (within `spread` seconds along the track around its middle; else the
# nearest usable point along the path; with `skip=true` an offset without any is left out).
# Its equipment sets the timing accuracy: `σ` [s] of a contact time for timings, `exptime` [s]
# and `noise` for simulated light curves.

using Occultations, Random
import GeoMakie

const LAND = [[(Float64(p[1]), Float64(p[2])) for p in GeoMakie.GeometryBasics.coordinates(poly.exterior)]
              for poly in GeoMakie.land()]

# point in polygon (ray casting)
function _inside(ring, x, y)
    c = false
    j = length(ring)
    for i in eachindex(ring)
        (xi, yi), (xj, yj) = ring[i], ring[j]
        ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi) && (c = !c)
        j = i
    end
    return c
end
onland(lon, lat) = any(r -> _inside(r, mod(lon + 180, 360) - 180, lat), LAND)

# Cross-track offsets of a campaign: random within [lo, hi], sorted, at least `mingap` km apart.
# Real stations are not evenly spaced: some cluster, others leave gaps.
function random_offsets(rng, n, lo, hi; mingap=5.0)
    for _ in 1:10_000
        x = sort(lo .+ (hi - lo) .* rand(rng, n))
        all(diff(x) .>= mingap) && return x
    end
    error("no $n offsets at least $mingap km apart in [$lo, $hi]")
end

# typical equipment: name, timing accuracy of a contact [s], exposure [s], relative noise
const EQUIPMENT = [("video, GPS-timed", 0.02, 0.04, 0.10), ("CMOS camera, GPS", 0.03, 0.1, 0.05),
                   ("CMOS camera, NTP", 0.10, 0.2, 0.04), ("CCD, long cycle", 0.25, 0.5, 0.03),
                   ("visual / poorly timed", 0.40, 1.0, 0.03)]

# ground point of the cross-track offset δ at tref + τ, or nothing
function _ground(star, body, tref, τ, δ)
    p = shadow_path(star, body, tref + τ, tref + τ; offsets=[δ])
    isnan(p.lon[1, 1]) && return nothing
    return p.lon[1, 1], p.lat[1, 1]
end

# usable for an observation: on land, Sun below −12°, star above 20°
function _usable(star, body, tref, τ, δ)
    g = _ground(star, body, tref, τ, δ)
    g === nothing && return false
    lon, lat = g
    onland(lon, lat) || return false
    sun_altitude(lon, lat, tref + τ) < -12 || return false
    ra, dec = geocentric(star, tref + τ)
    return altaz(Observer("x", lon, lat, 0.0), tref + τ, ra, dec)[1] > 20
end

function stations(star, body, tref, offsets; rng, window=600.0, spread=60.0, prefix="S", equipment=EQUIPMENT,
                  skip=false)
    # the longest stretch of the centre line over land and in the dark
    τs = collect(-window:5.0:window)
    ok = [_usable(star, body, tref, τ, 0.0) for τ in τs]
    best, cur = (0, 0), (0, 0)
    for i in eachindex(ok)
        cur = ok[i] ? (cur[1] == 0 ? (i, i) : (cur[1], i)) : (0, 0)
        cur[1] > 0 && cur[2] - cur[1] > best[2] - best[1] && (best = cur)
    end
    best[1] == 0 && error("the path crosses no land in the dark within ±$window s")
    a, b = τs[best[1]], τs[best[2]]
    # a regional campaign: `spread` seconds along the track around the middle of the stretch
    m = (a + b) / 2
    a, b = max(a, m - spread / 2), min(b, m + spread / 2)
    out = []
    for (j, δ) in enumerate(offsets)
        τ = nothing
        for _ in 1:500
            t = a + (b - a) * rand(rng)
            _usable(star, body, tref, t, δ) && (τ = t; break)
        end
        # otherwise anywhere along the path (an outlying station), or none
        if τ === nothing
            ok = [t for t in -window:5.0:window if _usable(star, body, tref, t, δ)]
            isempty(ok) || (τ = ok[argmin(abs.(ok .- (a + b) / 2))])
        end
        if τ === nothing
            skip && continue
            error("no usable site for the offset $δ km")
        end
        lon, lat = _ground(star, body, tref, τ, δ)
        push!(out, (observer=Observer("$prefix$j", lon, lat, 0.0), τ=τ, equipment=equipment[rand(rng, 1:length(equipment))]))
    end
    return out
end
