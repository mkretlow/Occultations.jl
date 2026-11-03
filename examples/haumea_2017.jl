# # Haumea and its ring, 2017 January 21
#
# The stellar occultation by the dwarf planet (136108) Haumea on 2017 January 21, reduced from
# the published timings of [Ortiz2017](@citet), Nature 550, 219, "The size, shape, density and
# ring of the dwarf planet Haumea from a stellar occultation": twelve chords of the body from
# ten sites in central Europe, and brief ring events before and after the main occultation at
# eight sites (one of which, Mount Agliale, saw only the ring).
#
# The example covers:
# * the limb of the body (an ellipse);
# * the ring from the mid-times of the ring events, with [`RingedLimb`](@ref): its radius, opening
#   angle and position angle, and its pole – a single event leaves two mirror solutions, as in
#   the paper; the ring's width and opacity (70 km, 0.5) are those of the paper (mid-times do
#   not constrain them);
# * the 3-D ellipsoid, with the ring plane as Haumea's equator (rotation pole = ring pole), the
#   rotation period, the rotational phase at the event (absolute brightness minimum) and the
#   2009 light-curve amplitude of 0.32 mag;
# * the astrometry.
#
# The data are transcribed in [`examples/data/haumea_20170121.jl`](https://github.com/mkretlow/Occultations.jl/blob/main/examples/data/haumea_20170121.jl),
# with the stored JPL Horizons ephemeris, so the script runs offline:
#
#     julia --project=. examples/haumea_2017.jl [output_dir]     (default: examples/output/haumea_2017)

using Occultations, Random, Printf, Statistics, LinearAlgebra
using Plots
import CairoMakie                     # saves the GeoMakie maps (Makie figures)
include(joinpath(@__DIR__, "standard_figures.jl"))
ENV["GKSwstype"] = "100"   #src
include(joinpath(@__DIR__, "data", "haumea_20170121.jl"))

out = length(ARGS) >= 1 ? ARGS[1] : mkpath(joinpath(@__DIR__, "output", "haumea_2017"))
rng = MersenneTwister(2017)
section(s) = println("\n", "="^78, "\n", s, "\n", "="^78)   #src
P = HAUMEA_PAPER

section("Event and geometry (paper values in parentheses)")
ev = haumea_event()
D = ev.ca.dist_km
@printf("%d chords (%d with ring events); distance %.2f au, shadow velocity %.2f km/s (%.1f)\n", length(ev.chords),
        count(c -> any(x -> x.body === :ring1, c.contacts), ev.chords), D / 1.495978707e8, ev.ca.vel, P.velocity)
# The geocentric shadow velocity from the Horizons ephemeris is about 10 % above the paper's
# 13.1 km/s; Earth's rotation at the European sites changes it by only ~0.3 km/s. The fits use
# the ephemeris velocity throughout, and the sizes below agree with the paper's: the paper's
# chord lengths must rest on a velocity close to this one, and its 13.1 km/s is unexplained here.

section("The limb of the body")
# The twelve chords; Mount Agliale, which saw only the ring, is a negative chord for the body.
rb = fit_limb(haumea_event(rings=false); rng=rng)
@printf("a′ = %.1f ± %.1f km (%.0f ± %.0f), b′ = %.1f ± %.1f km (%.0f ± %.0f), PA = %.1f ± %.1f° (%.1f ± %.1f), χ²/dof = %.2f\n",
        rb[:equatorial_radius]..., P.a..., rb[:semi_minor_axis]..., P.b..., mod(rb[:position_angle][1] + 90, 180) - 90,
        rb[:position_angle][2], P.position_angle..., reduced_chi2(rb))
rbs = fit_limb(haumea_event(rings=false); σ_limb=:auto, rng=rng)
@printf("with σ_limb = %.1f km: a′ = %.1f ± %.1f, b′ = %.1f ± %.1f km, PA = %.1f ± %.1f°\n", rbs.problem.σ_limb,
        rbs[:equatorial_radius]..., rbs[:semi_minor_axis]..., mod(rbs[:position_angle][1] + 90, 180) - 90,
        rbs[:position_angle][2])

section("The ring")
# The body and a circular ring in the plane of the ring pole, concentric with it. The ring
# events are mid-times: they fix the projected ring ellipse (radius, |B| and P), not its
# width; width and optical depth are fixed at the paper's 70 km and apparent opacity 0.5
# (τ_N = −|sin B| ln(1 − p′)). The plane and its mirror image give the same ellipse, and the
# same χ²: the fit starts from the paper's preferred pole, then from the mirror plane of its
# result (and back, if the mirror fits better – the χ² surface has local minima).
τN = -sind(P.ring_B[1]) * log(1 - P.ring_opacity)
model = RingedLimb(EllipseLimb(); rings=1)
ringfit(pole) = fit_limb(ev; model=model, fix=[:ring1_width, :ring1_tau],
                         init=Dict(:ring1_width => P.ring_width, :ring1_tau => τN, :ring_pole_ra => pole[1],
                                   :ring_pole_dec => pole[2], :ring1_radius => P.ring_a[1]), starts=4, rng=rng)
r1 = ringfit((P.ring_pole1[1][1], P.ring_pole1[2][1]))
u = EventContext(ev).u
## the mirror plane: all other parameters from the given solution, the pole mirrored (opening −B)
function mirrorfit(r)
    init = Dict(n => v for (n, v) in zip(r.names, r.values))
    init[:ring_pole_ra], init[:ring_pole_dec] = Occultations._ring_pole(-r[:ring_opening][1], r[:ring_pa][1], u)
    return fit_limb(ev; model=model, fix=[:ring1_width, :ring1_tau], init=init, starts=4, rng=rng)
end
r2 = mirrorfit(r1)
## the fit is multimodal: if the mirror found a lower χ², solution A restarts from its mirror
r2.chi2 < r1.chi2 - 1e-3 && (r1 = mirrorfit(r2))
antipode(α, δ) = (mod(α + 180, 360), -δ)
for (label, r) in (("solution A", r1), ("solution B (mirror)", r2))
    α, δ = r[:ring_pole_ra][1], r[:ring_pole_dec][1]
    @printf("%-20s χ² = %6.2f: ring radius %.0f ± %.0f km, B = %+.2f ± %.2f°, P = %.2f°; pole (%.1f, %+.1f)°, antipode (%.1f, %+.1f)°\n",
            label, r.chi2, r[:ring1_radius]..., r[:ring_opening]..., r[:ring_pa][1], α, δ, antipode(α, δ)...)
end
R, B, Pr = r1[:ring1_radius][1], r1[:ring_opening][1], r1[:ring_pa][1]
@printf("apparent ring ellipse: a′ = %.0f km (%.0f −%.0f +%.0f), b′ = %.0f km (%.0f ± %.0f), PA of the minor axis %.1f° (%.1f ± %.1f), |B| = %.1f° (%.1f ± %.1f)\n",
        R, P.ring_a..., R * abs(sind(B)), P.ring_b..., mod(Pr + 90, 180) - 90, P.ring_pa..., abs(B), P.ring_B...)
@printf("published poles: solution 1 (%.1f, %.1f)°, solution 2 (%.1f, %.1f)°; Hi'iaka's orbit pole (%.1f, %.1f)°\n",
        P.ring_pole1[1][1], P.ring_pole1[2][1], P.ring_pole2[1][1], P.ring_pole2[2][1], P.hiiaka_pole...)
println("(a pole and its antipode describe the same plane; the fit reports the northern one)")
for x in check_negative_chords(rbs)                  # Mount Agliale, which saw only the ring
    @printf("  %s (ring only): %s, %.0f km from the limb of the body\n", x.chord, x.conflict ? "CROSSES" : "consistent",
            x.distance)
end

section("3-D shape: the ring plane as the equator")
# Haumea's equator in the ring plane: the rotation pole is the ring pole (each mirror
# solution in turn, fixed). The rotation period is known (3.915341 h); at the event Haumea was
# at its absolute brightness minimum (rotational phase 0.25 after a maximum, modulo ½ for an
# ellipsoid); the HST light curve of 2009 gives an amplitude of 0.32 mag (mid-2009 assumed).
t = ev.tref
amp = RotationalAmplitude(Epoch("2009-07-01T00:05:00"), P.amplitude_2009, 0.02)
phase = RotationalPhase(t, 0.25, 0.02)
for (label, r) in (("solution A", r1), ("solution B", r2))
    α, δ = r[:ring_pole_ra][1], r[:ring_pole_dec][1]
    prob = LimbProblem(haumea_event(rings=false); model=EllipsoidLimb(principal=true), amplitudes=[amp], phases=[phase],
                       σ_limb=rbs.problem.σ_limb, fix=[:pole_ra, :pole_dec],
                       init=Dict(:a => 1100.0, :b => 850.0, :c => 550.0, :pole_ra => α, :pole_dec => δ, :W0 => 0.0,
                                 :period => P.period_h, :center_f => rbs[:center_f][1], :center_g => rbs[:center_g][1]),
                       bounds=Dict(:a => (300.0, 2500.0), :b => (300.0, 2500.0), :c => (200.0, 2500.0)))
    sols = [fit(prob; rng=MersenneTwister(k), starts=12) for k in 1:3]
    r3 = sols[argmin([s.chi2 for s in sols])]
    @printf("%s: a = %.0f ± %.0f, b = %.0f ± %.0f, c = %.0f ± %.0f km, D_eq,V = %.0f km, aspect %.1f°, χ² = %.2f\n", label,
            r3[:a]..., r3[:b]..., r3[:c]..., 2r3[:volume_equivalent_radius][1], r3[:aspect][1], r3.chi2)
    label == "solution A" && (global r3A = r3)
end
@printf("published: a = %.0f ± %.0f, b = %.0f ± %.0f, c = %.0f ± %.0f km, D_eq,V = %.0f ± %.0f km\n", P.axes_3d[1]...,
        P.axes_3d[2]..., P.axes_3d[3]..., P.D_eq_volume...)
# Solution A reproduces the paper's ellipsoid. The mirror pole, which fits the ring equally
# well, gives a nearly equator-on view (aspect 76°) in which the axes a and c are degenerate:
# the occultation alone cannot choose between the two, the paper prefers solution 1 (A) because
# the ring is then coplanar with Hi'iaka's orbit.

section("Astrometry")
a = astrometry(rbs)
@printf("RA %s  Dec %s at %s UTC (σ %.2f, %.2f mas); offset from JPL Horizons Δα cos δ = %+.2f, Δδ = %+.2f mas\n",
        a.ra_str, a.dec_str, isot(a.epoch; digits=1), a.ra_err, a.dec_err, a.offset_ra, a.offset_dec)
println("(the position refers to the star position of the paper; Hi'iaka shifts Haumea from the system's")
println(" barycentre by up to ≈ 1232/(1 + r) mas, r the mass ratio)")

section("Figures")
savefig(plot(r1; title="Haumea, 2017-01-21: body and ring", size=(1000, 900)), joinpath(out, "haumea_ring.png"))
savefig(plot(rbs; title="Haumea, 2017-01-21: limb", size=(1000, 750)), joinpath(out, "haumea_limb.png"))
CairoMakie.save(joinpath(out, "haumea_occview.png"), occview(r3A))
## the standard figures: limb fits (default, 68.3 % region, as the paper), maps (prediction, post-fit,
## comparison) and the globe
standard_figures("haumea", rbs; out=out, paper=(fit=rb, uncertainty=:pointwise, nsigma=1),
                 projection="+proj=laea +lon_0=14 +lat_0=47", limits=(2, 28, 38, 56), window=600, marks=60)
println("figures written to ", out)

# The body and the ring with the chords and the ring events:
#
#md # ![ring](../assets/figures/haumea_ring.png)
#
# The limb of the body, and the 3-D ellipsoid with the ring pole as its rotation axis:
#
#md # ![limb](../assets/figures/haumea_limb.png)
#md # ![occview](../assets/figures/haumea_occview.png)
