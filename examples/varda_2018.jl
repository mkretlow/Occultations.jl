# # Varda, 2018 September 10
#
# The stellar occultation by (174567) Varda on 2018 September 10, reduced from the published
# timings of [Souami2020](@citet), A&A 643, A125,
# "A multi-chord stellar occultation by the large trans-Neptunian object (174567) Varda": five
# positive chords and six negative chords.
#
# The example covers the chord lengths, the limb fit (Table 5), the posterior, the geometric
# albedo and the astrometry. The limb and χ² are reproduced, and the astrometry to 0.01 mas.
#
# The data are transcribed in [`examples/data/varda_20180910.jl`](https://github.com/mkretlow/Occultations.jl/blob/main/examples/data/varda_20180910.jl),
# with the stored JPL Horizons ephemeris, so the script runs offline:
#
#     julia --project=. examples/varda_2018.jl [output_dir]     (default: examples/output/varda_2018)

using Occultations, Random, Printf, Statistics
using Plots
import CairoMakie                     # saves the GeoMakie maps (Makie figures)
include(joinpath(@__DIR__, "standard_figures.jl"))
ENV["GKSwstype"] = "100"   #src
include(joinpath(@__DIR__, "data", "varda_20180910.jl"))

out = length(ARGS) >= 1 ? ARGS[1] : mkpath(joinpath(@__DIR__, "output", "varda_2018"))
rng = MersenneTwister(2018)
section(s) = println("\n", "="^78, "\n", s, "\n", "="^78)   #src
P = VARDA_PAPER
q(x) = quantile(x, [0.15865525393145707, 0.5, 0.8413447460685429])
pm(x; d=1) = (v = q(x); Printf.format(Printf.Format("%.$(d)f -%.$(d)f +%.$(d)f"), v[2], v[2] - v[1], v[3] - v[2]))

section("Event, geometry and chords (paper values in parentheses)")
ev = varda_event()
D = ev.ca.dist_km
@printf("%d positive, %d negative chords; distance %.4f au (%.4f), shadow velocity %.2f km/s (%.2f), stellar diameter %.2f km (%.2f)\n",
        count(ispositive, ev.chords), count(isnegative, ev.chords), D / 1.495978707e8, P.distance_au, ev.ca.vel, P.velocity,
        apparent_diameter(ev.star, D), P.star_diameter_km)
for c in ev.chords
    ispositive(c) || continue
    tr = track(ev, c)
    ti, te = c.contacts[1].time, c.contacts[2].time
    f1, g1 = tr(ti - ev.tref); f2, g2 = tr(te - ev.tref)
    @printf("  %-4s duration %.2f s (%.2f), chord %.1f km (%.1f ± %.1f)\n", c.name, te - ti, P.durations[c.name][1],
            hypot(f2 - f1, g2 - g1), P.chords[c.name]...)
end

section("Limb fit (Table 5; P′ of the minor axis from north towards east)")
row(label, r) = @printf("%-30s a′ = %5.1f ± %4.1f  ε′ = %.3f ± %.3f  P′ = %5.1f ± %4.1f°  R′_eq = %5.1f ± %4.1f km  χ²/dof = %.2f\n",
                        label, r[:equatorial_radius]..., r[:oblateness]..., mod(r[:position_angle][1] + 90, 180) - 90,
                        r[:position_angle][2], r[:equivalent_radius]..., reduced_chi2(r))
@printf("%-30s a′ = %5.1f ± %4.1f  ε′ = %.3f ± %.3f  P′ = %5.1f ± %4.1f°  R′_eq = %5.1f ± %4.1f km  χ²/dof = %.2f\n",
        "published (N = 10)", P.a..., P.oblateness..., P.position_angle..., P.R_equiv..., P.chi2_pdf)
rr = fit_limb(ev; residual=:radial, rng=rng)
rc = fit_limb(ev; residual=:chord, rng=rng)
row("radial residuals (as the paper)", rr)
row("residuals along the chords", rc)
println("closest negative chords:")
for x in first(check_negative_chords(rr), 3)
    @printf("  %-6s %s, %6.1f km from the limb\n", x.chord, x.conflict ? "CROSSES THE LIMB" : "consistent", x.distance)
end

section("Posterior and albedo (H_V = 3.81 ± 0.01)")
post = sample(rr; nsteps=5000, rng=rng)
S = samples(post)
col(n) = S[:, findfirst(==(n), post.names)]
A, ε = col(:equatorial_radius), col(:oblateness)
Req = A .* sqrt.(1 .- ε)
@printf("a′ = %s km, ε′ = %s, P′ = %s°, R′_eq = %s km\n", pm(A), pm(ε; d=3), pm(mod.(col(:position_angle) .+ 90, 180) .- 90),
        pm(Req))
## p_V = (au / R_equiv)² 10^{0.4 (H_sun − H)}, H_sun = −26.74 (footnote 8 of the paper)
pV = [(149597870.7 / R)^2 * 10^(0.4 * (-26.74 - (P.H_V[1] + P.H_V[2] * randn(rng)))) for R in Req]
@printf("geometric albedo %s (paper %.3f ± %.3f; R′_eq ± 7 km alone gives ± %.3f)\n", pm(pV; d=3), P.albedo...,
        2 * P.albedo[1] * P.R_equiv[2] / P.R_equiv[1])

section("Astrometry (Table 5)")
t = Epoch(P.astrometry.epoch)
a = astrometry(rr)
ra_, dec_, _ = ephem_radec(ev.body, t)
r_, d_ = Occultations.offset_radec(ra_, dec_, a.offset_ra, a.offset_dec)
@printf("  here      %s UTC: RA %s  Dec %s  (σ %.2f, %.2f mas)\n", P.astrometry.epoch, format_ra(r_; digits=5),
        format_dec(d_; digits=4), a.ra_err, a.dec_err)
@printf("  published                          RA %s  Dec %s  (σ %.3f, %.3f mas)\n", P.astrometry.ra, P.astrometry.dec,
        P.astrometry.ra_err, P.astrometry.dec_err)

section("Figures")
savefig(plot(ev; title="Varda, 2018-09-10: chords", size=(1000, 750)), joinpath(out, "varda_chords.png"))
savefig(cornerplot(post; vars=[:equatorial_radius, :oblateness, :position_angle],
                   axislabels=Dict(:equatorial_radius => "a′ [km]", :oblateness => "ε′", :position_angle => "P′ [°]"),
                   truths=Dict(:equatorial_radius => P.a[1], :oblateness => P.oblateness[1],
                               :position_angle => P.position_angle[1]), truth_label="Souami et al. 2020"),
        joinpath(out, "varda_corner.png"))
## the standard figures: limb fits (default, 68.3 % region, as the paper), posterior, maps (prediction,
## post-fit, comparison) and the globe
standard_figures("varda", rr; out=out, posterior=post, error_mas=43.1, error_note="NIMA v6 across-track", 
                 projection="+proj=laea +lon_0=-108 +lat_0=38", limits=(-128, -84, 28, 50), window=600, marks=60)
println("figures written to ", out)

# The limb fit (radial residuals, as the paper), and the prediction and post-fit shadow path;
# the prediction uses the across-track uncertainty of the NIMA v6 ephemeris (43.1 mas):
#
#md # ![limb fit](../assets/figures/varda_limb_fit.png)
#md # ![maps](../assets/figures/varda_geomap_comparison.png)
#
# The results are discussed in [Comparison with published reductions](@ref "(174567) Varda, 2018 September 10").
