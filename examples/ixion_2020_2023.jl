# # Ixion, six events 2020–2023
#
# Six multi-chord stellar occultations by (28978) Ixion, 2020–2023, reduced from the published
# timings of [Kilic2026](@citet), A&A 707, A70,
# "Constraining the size, shape, and albedo of the large trans-Neptunian object (28978) Ixion
# with multi-chord stellar occultations".
#
# The example covers:
# * the geometry of each event and ellipse fits per event (Table 4);
# * the global fit of one ellipse with a centre per event and a limb allowance of 4 km
#   (Eq. 4), reproduced (a = 362.6 vs 363.4 km, R_eq = 349.8 vs 348.4 km);
# * the posterior, the geometric albedo and the astrometry (Table 5).
#
# The paper mixes two position-angle conventions, and four of its six astrometric positions
# differ from the reduction here by 30–50 mas.
#
# The data are transcribed in [`examples/data/ixion_2020_2023.jl`](https://github.com/mkretlow/Occultations.jl/blob/main/examples/data/ixion_2020_2023.jl),
# with the stored JPL Horizons ephemerides, so the script runs offline:
#
#     julia --project=. examples/ixion_2020_2023.jl [output_dir]   (default: examples/output/ixion_2020_2023)

using Occultations, Random, Printf, Statistics
using Plots
import CairoMakie                     # saves the GeoMakie maps (Makie figures)
include(joinpath(@__DIR__, "standard_figures.jl"))
ENV["GKSwstype"] = "100"   #src
include(joinpath(@__DIR__, "data", "ixion_2020_2023.jl"))

out = length(ARGS) >= 1 ? ARGS[1] : mkpath(joinpath(@__DIR__, "output", "ixion_2020_2023"))
rng = MersenneTwister(2026)
section(s) = println("\n", "="^78, "\n", s, "\n", "="^78)   #src
P = IXION_PAPER
q(x) = quantile(x, [0.15865525393145707, 0.5, 0.8413447460685429])
pm(x; d=1) = (v = q(x); Printf.format(Printf.Format("%.$(d)f -%.$(d)f +%.$(d)f"), v[2], v[2] - v[1], v[3] - v[2]))
dates = [e[1] for e in IXION_EVENTS]

section("Events and geometry (Table 3 in parentheses)")
evs = ixion_events()
for (d, ev) in zip(dates, evs)
    g = P.geometry[d]
    @printf("%s: %d positive, %d negative chords; distance %.2f au (%.2f), velocity %.2f km/s (%.2f)\n", d,
            count(ispositive, ev.chords), count(isnegative, ev.chords), ev.ca.dist_km / 1.495978707e8, g[1], ev.ca.vel, g[2])
end

section("Ellipse fits per event (Table 4; radial residuals as SORA)")
# ϕ of the paper (SORA's position angle) is that of the minor axis (pole) from north towards
# east, as the position angle here; it is reduced to (−90°, 90°]:
φ(r) = mod(r[:position_angle][1] + 90, 180) - 90
for (d, ev) in zip(dates, evs)
    p = P.events[d]
    r = try
        fit_limb(ev; residual=:radial, σ_limb=P.sigma_model, rng=rng)    # σ_model = 4 km as in the global fit
    catch e
        nothing
    end
    @printf("%s paper: a = %5.1f ± %4.1f  b = %5.1f ± %4.1f  ϕ = %6.1f ± %5.1f  R_eq = %5.1f ± %4.1f  χ²pdf = %.3f\n",
            d, p.a..., p.b..., p.phi..., p.R_equiv..., p.chi2_pdf)
    r === nothing && (println("           (no solution)"); continue)
    @printf("%s here:  a = %5.1f ± %4.1f  b = %5.1f ± %4.1f  ϕ = %6.1f ± %5.1f  R_eq = %5.1f ± %4.1f  χ²/dof = %s\n",
            d, r[:equatorial_radius]..., r[:semi_minor_axis]..., φ(r), r[:position_angle][2], r[:equivalent_radius]...,
            r.dof > 0 ? @sprintf("%.3f", reduced_chi2(r)) : "– (no dof)")
end

section("Global fit: one ellipse, a centre per event, σ_model = 4 km (Eq. 4)")
rg = fit_limb(evs; residual=:radial, σ_limb=P.sigma_model, rng=rng, starts=8)
n = sum(length(c.contacts) for ev in evs for c in ev.chords if ispositive(c))
@printf("paper: a = %.2f +%.2f −%.2f  b = %.2f +%.2f −%.2f  ϕ = %.2f +%.2f −%.2f  R_eq = %.2f +%.2f −%.2f km  χ²ν = %.2f (N = %d, ν = %d)\n",
        P.a..., P.b..., mod(P.phi_major[1] + 90, 180) - 90, P.phi_major[2:3]..., P.R_equiv..., P.chi2_nu, P.ndata, P.dof)
println("(the paper's ϕ = 110.97° is given as −69.03° here, the same axis)")
@printf("here:  a = %.2f ± %.2f  b = %.2f ± %.2f  ϕ = %.2f ± %.2f  R_eq = %.2f ± %.2f km  χ²ν = %.2f (N = %d contact times, ν = %d)\n",
        rg[:equatorial_radius]..., rg[:semi_minor_axis]..., φ(rg), rg[:position_angle][2], rg[:equivalent_radius]...,
        reduced_chi2(rg), n, rg.dof)
println("(the paper's N = 56 and ν = 41 imply six events, M = 3 + 2·6 = 15; Table A.1 lists 55")
println(" contact times for the six multi-chord events)")
println("negative chords closest to the limb:")
for (k, ev) in enumerate(evs)
    for x in first(check_negative_chords(rg; event=k), 1)
        @printf("  %s %-18s %s, %6.1f km from the limb\n", dates[k], x.chord, x.conflict ? "CROSSES THE LIMB" : "consistent",
                x.distance)
    end
end

section("Posterior, diameter and albedo (H_V = 3.845 ± 0.006)")
post = sample(rg; nsteps=4000, rng=rng)
dp = derived_samples(post; n=4000, rng=rng)
S = samples(post)
col(nm) = S[:, findfirst(==(nm), post.names)]
Req = dp.equivalent_radius
pV = [(1329.0 / (2R) * 10^(-(P.H_V[1] + P.H_V[2] * randn(rng)) / 5))^2 for R in Req]
@printf("a = %s km, b = %s km, R_eq = %s km (paper %.2f +%.2f −%.2f), D_eq = %s km\n", pm(col(:equatorial_radius)),
        pm(dp.semi_minor_axis), pm(Req), P.R_equiv..., pm(2 .* Req))
@printf("geometric albedo %s (paper %.3f ± %.3f)\n", pm(pV; d=3), P.albedo...)

section("Astrometry (Table 5)")
for (k, (d, ev)) in enumerate(zip(dates, evs))
    a = astrometry(rg; event=k)
    p = P.events[d].astrometry
    t = Epoch(d * "T" * p[1])
    ra_, dec_, _ = ephem_radec(ev.body, t)
    r_, d_ = Occultations.offset_radec(ra_, dec_, a.offset_ra, a.offset_dec)
    dra = (r_ - parse_ra(p[2])) * cosd(dec_) * 3.6e6
    ddec = (d_ - parse_dec(p[3])) * 3.6e6
    @printf("%s %s: RA %s Dec %s (σ %.2f, %.2f mas)  paper − here: %+6.2f, %+6.2f mas\n", d, p[1],
            format_ra(r_; digits=5), format_dec(d_; digits=4), a.ra_err, a.dec_err, -dra, -ddec)
end

section("Figures")
for (k, d) in enumerate(dates)
    savefig(plot(rg; event=k, title="Ixion, " * d * " (global fit)", size=(1000, 750)), joinpath(out, "ixion_$(d).png"))
    ## as in the paper (Fig. E.2 / Fig. 2: 1σ envelope = pointwise 1σ)
    savefig(plot(rg; event=k, uncertainty=:pointwise, title="Ixion, " * d * " (as the paper: pointwise 1σ region)",
                 size=(1000, 750)), joinpath(out, "ixion_$(d)_paper.png"))
end
savefig(cornerplot(post; vars=[:equatorial_radius, :oblateness, :position_angle],
                   axislabels=Dict(:equatorial_radius => "a [km]", :oblateness => "ε′", :position_angle => "ϕ [°]"),
                   truths=Dict(:equatorial_radius => P.a[1], :oblateness => P.oblateness[1],
                               :position_angle => P.phi_major[1]), truth_label="Kilic et al. 2026"),
        joinpath(out, "ixion_corner.png"))
## the standard figures: limb fits (default, 68.3 % region, as the paper), posterior, maps (prediction,
## post-fit, comparison) and the globe
standard_figures("ixion", rg; out=out, posterior=post, event=2, 
                 projection="+proj=laea +lon_0=-105 +lat_0=35", limits=(-125, -85, 25, 47), window=600, marks=30)
println("figures written to ", out)

# The global fit for the event of 2020 October 13, and the posterior of the ellipse with the
# published values:
#
#md # ![limb fit](../assets/figures/ixion_limb_fit.png)
#md # ![corner plot](../assets/figures/ixion_corner.png)
#
# The results are discussed in [Comparison with published reductions](@ref "(28978) Ixion, six events 2020–2023").
