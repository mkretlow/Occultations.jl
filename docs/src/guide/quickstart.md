# Quick start

A complete reduction in outline, with a real event whose elliptical limb agrees closely with
the publication: star and body, published timings and a negative chord, the limb fit with its
posterior and astrometry, and maps. Each step is described in detail in the following
chapters.

```julia
using Occultations, Plots

# the occultation of Gaia DR2 4367203805493754752 by (174567) Varda on 2018 September 10
# (Souami et al. 2020); the star at the event, the body with its ephemeris from JPL Horizons
# (`using HTTP`; or a stored table, an orbit or a SPICE kernel)
using HTTP
star = Star("17 18 25.12481", "-02 05 14.4282"; epoch=2018.69, ra_err=0.2, dec_err=0.1, mag=Dict("G" => 14.7))
body = Body("Varda", horizons_ephem("174567", "2018-09-10T03:00", "2018-09-10T04:00"); radius=370.0)
ev   = Event(star, body, "2018-09-10T03:34")              # closest approach searched around this time

# the published timings of the positive chords (UTC, 1σ) and a negative chord
day = "2018-09-10T"
for (name, lon, lat, alt, ti, σi, te, σe) in [
        ("TWF", "-114 28 13.1", "42 35 01.9", 1133, "03:32:57.10", 0.60, "03:33:58.83", 0.60),
        ("MHV", "-114 35 48.9", "35 01 54.1",  184, "03:33:56.99", 0.76, "03:35:05.27", 0.30),
        ("YMA", "-114 26 10.5", "32 39 34.0",   97, "03:34:17.10", 0.60, "03:35:26.39", 0.59),
        ("CAR", "-111 57 07.9", "33 48 42.9",  654, "03:34:12.39", 0.47, "03:35:19.52", 0.47),
        ("FLO", "-111 21 00.6", "33 00 54.3",  484, "03:34:23.30", 0.13, "03:35:26.75", 0.15)]
    push!(ev, Chord(name, Observer(name, lon, lat, alt); immersion=day * ti, immersion_err=σi,
                    emersion=day * te, emersion_err=σe))
end
push!(ev, Chord("GAR", Observer("GAR", "-119 40 20.3", "38 53 23.5", 1502);
                span=(day * "03:28", day * "03:42"), sensitivity=1.068))

# your own light curve instead: timings from the diffraction model
# lc = read_lightcurve("station.csv"; format=:pymovie, exptime=0.5, λ=0.6, Δλ=0.3)
# tf = fit_timing(normalize(lc); vel=ev.ca.vel, dist=ev.ca.dist_km,
#                 star_diam=apparent_diameter(star, ev.ca.dist_km))
# plot(tf); push!(ev, Chord("station", Observer("station", lon, lat, alt), tf))

# the elliptical limb (negative chords as hard constraints)
r = fit_limb(ev)
r[:equatorial_radius], r[:oblateness], r[:equivalent_radius]   # (value, σ): ≈ 383 km, 0.07, 370 km
reduced_chi2(r)
astrometry(r)                                             # RA/Dec of the body, offsets, errors incl. the star
plot(r)                                                   # chords, fitted limb, 1σ limb draws

# the posterior, with the published values for comparison
post = sample(r; nsteps=3000)
summarize(post); derived_summary(post)
cornerplot(post; vars=[:equatorial_radius, :oblateness], axislabels=Dict(:equatorial_radius => "a′ [km]"),
           truths=Dict(:equatorial_radius => 383.0, :oblateness => 0.066), truth_label="Souami et al. 2020")
confidence_interval(r, :oblateness)                       # from the χ² profile

# maps
import CairoMakie                                         # backend for saving Makie figures
CairoMakie.save("varda_map.png", geomap(r; window=600, projection="+proj=laea +lon_0=-110 +lat_0=38",
                                        limits=(-130, -85, 25, 50)))
CairoMakie.save("varda_globe.png", geomap(r; projection="+proj=ortho +lon_0=-110 +lat_0=30"))
```

The complete reduction of this event, with its comparison to the paper, is the
[Varda example](../generated/varda_2018.md); [`examples/tutorial.jl`](https://github.com/mkretlow/Occultations.jl/blob/main/examples/tutorial.jl)
runs a reduction from synthetic light curves offline.
