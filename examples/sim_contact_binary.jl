# # Synthetic: a contact binary in one event
#
# A simulated occultation by a bilobed trans-Neptunian object (two lobes of 260 and 200 km,
# 430 km long; Arrokoth-like, scaled up), reduced with [`ContactBinaryLimb`](@ref):
#
# 1. nine stations at random offsets, on land, with light curves (exposures of 0.04–1 s,
#    diffraction, stellar diameter) and timing errors of 0.05–0.1 s; the chords near the waist see two drops;
# 2. the drops are detected and timed in each light curve; chords with several drops; a
#    second timing with the velocities normal to the fitted limb;
# 3. a single ellipse cannot reproduce two drops on one chord;
# 4. the contact-binary model: the lobes, the neck and the waist, with the posterior;
# 5. a joint fit of the light curves: the diffraction of the narrow gaps between the drops
#    (complex amplitudes of the two opaque strips added) is modelled exactly, which a timing
#    of each drop on its own is not.
#
# The truth is known and printed next to the fitted values. Run it with:
#
#     julia --project=. examples/sim_contact_binary.jl [output_dir]      (default: examples/output/sim_contact_binary)

using Occultations, Random, Printf, Statistics
using Plots
ENV["GKSwstype"] = "100"   #src

out = length(ARGS) >= 1 ? ARGS[1] : mkpath(joinpath(@__DIR__, "output", "sim_contact_binary"))
rng = MersenneTwister(11)
section(s) = println("\n", "="^78, "\n", s, "\n", "="^78)   #src

section("1. the simulated event")
eph = EphemKeplerian(a=44.0, e=0.03, i=3.0, Ω=120.0, ω=40.0, M=300.0, epoch=Epoch("2024-01-01"), name="SimCB")
body = Body("SimCB", eph; radius=150.0)
t0 = Epoch("2024-12-30T20:00:00")
ra, dec, _ = ephem_radec(body, t0)
star = Star(ra + 0.001 / 3600 / cosd(dec), dec; name="SimStar", mag=Dict("V" => 13.0, "K" => 12.0), ra_err=0.1, dec_err=0.1)
ev = Event(star, body, t0; name="SimCB")                # simulated: no date in the plots
D = ev.ca.dist_km
sd = apparent_diameter(star, D)
@printf("distance %.2f au, shadow velocity %.2f km/s, stellar diameter %.2f km, Fresnel scale %.2f km at 0.6 µm\n",
        D / 1.495978707e8, ev.ca.vel, sd, fresnel_scale(0.6, D))
# The truth: lobe 1 with semi-axes 130 × 95 km, lobe 2 with 100 × 80 km, along an axis at
# position angle 255° (19° from the direction of motion, so that chords near the waist see two
# drops); the centres are `d = (a₁ + a₂)(1 − neck)` apart, `neck = 0.08` (the
# lobes overlap by 8 % of their length); the system is 8 km west and 14 km north of the
# ephemeris position.
a1, b1, a2, b2, neck, axis_pa = 130.0, 95.0, 100.0, 80.0, 0.08, 255.0
center = (-8.0, 14.0)
d = (a1 + a2) * (1 - neck)
x1, x2 = -d * a2 * b2 / (a1 * b1 + a2 * b2), d * a1 * b1 / (a1 * b1 + a2 * b2)   # lobe centres, area-weighted centre at 0
truth = CompositeOutline([EllipseOutline(a1, b1, axis_pa + 90), EllipseOutline(a2, b2, axis_pa + 90)],
                         [(x1 * sind(axis_pa), x1 * cosd(axis_pa)), (x2 * sind(axis_pa), x2 * cosd(axis_pa))])
println("true limb: ", truth)
@printf("area-equivalent diameter %.1f km\n", 2 * equivalent_radius(truth))

section("2. light curves → drops → timings")
# Nine stations at random offsets across the path, on land and in the dark, spread along the
# track ([`examples/sim_stations.jl`](https://github.com/mkretlow/Occultations.jl/blob/main/examples/sim_stations.jl)),
# each with its own equipment, from GPS-timed video to visual timing: exposures of 0.04 to 1 s
# and different noise. Each station has a systematic timing error with σ of 0.05–0.1 s (GPS-
# synchronised timing), unknown to the reduction; as in real reports, this accuracy is added in
# quadrature to the error of each contact time. Each light
# curve is simulated with all its drops. The drops are found as runs of points below half the
# baseline (a box search would find two deep drops with a short gap as one box), and each is
# timed in its own window (between the neighbouring drops) with the diffraction model.
include(joinpath(@__DIR__, "sim_stations.jl"))
offsets = random_offsets(MersenneTwister(45), 9, -220.0, 210.0; mingap=12.0)
pv = positionv(star, eph, GEOCENTER, ev.tref)
n̂ = (-pv.vg, pv.vf) ./ hypot(pv.vf, pv.vg)
sites = stations(star, body, ev.tref, offsets .+ (center[1] * n̂[1] + center[2] * n̂[2]); rng=rng)
## the velocity of the shadow normal to a limb O (centre c) at the instant τ of a track
function vnormal(tr, O, c, τ)
    f, g = tr(τ)
    nx, ny = limb_normal(O, f - c[1], g - c[2])
    vf, vg = velocity(tr, τ)
    return abs(vf * nx + vg * ny)
end
## the contacts of the drops `found` in a light curve, each drop timed in its own window with the
## normal velocities vel(τ) at its edges
function time_drops(lc, found, vel; σt=0.0)
    cs = Contact[]
    for (k, x) in enumerate(found)
        lo = k == 1 ? -Inf : (found[k-1].emersion + x.immersion) / 2
        hi = k == length(found) ? Inf : (x.emersion + found[k+1].immersion) / 2
        tf = fit_timing(lc; vel=vel(x.immersion), vel_emersion=vel(x.emersion), dist=D, star_diam=sd, tmin=lo, tmax=hi,
                        immersion=x.immersion, emersion=x.emersion, rng=rng)
        c = contact_times(tf)
        push!(cs, Contact(:immersion, c.immersion, hypot(c.immersion_err, σt)),
              Contact(:emersion, c.emersion, hypot(c.emersion_err, σt)))
    end
    return cs
end
detections = Dict{String,Any}()
σclks = Dict{String,Float64}()                            # the timing accuracy of each station [s]
dclocks = Dict{String,Float64}()                          # and its actual timing error
for (j, st) in enumerate(sites)
    obs = st.observer
    _, _, texp, noise = st.equipment
    σclock = 0.05 + 0.05rand(rng)                         # the timing accuracy of this station
    dclock = σclock * randn(rng)                          # and its actual timing error
    σclks[obs.name] = σclock
    dclocks[obs.name] = dclock
    drops = [(a - ev.tref + dclock, b - ev.tref + dclock) for (a, b) in predicted_contacts(ev, obs, truth; center=center,
                                                                          window=abs(st.τ) + 200, step=0.02)]
    ## the velocity of the shadow normal to the limb at each contact sets the time scale of its diffraction fringes
    tr = Track(star, body, obs, ev.tref, st.τ - 60, st.τ + 60)
    t = collect(range(round(st.τ) - 30, round(st.τ) + 30; step=texp))
    lc = simulate_lightcurve(t, drops; tref=ev.tref, vel=[vnormal(tr, truth, center, a) for (a, _) in drops],
                             vel_emersion=[vnormal(tr, truth, center, b) for (_, b) in drops], dist=D, star_diam=sd,
                             λ=0.6, σ=noise, exptime=texp, rng=rng, name=obs.name)
    ## a drop must span at least three exposures: one or two low points are not a detection
    found = filter(x -> x.emersion - x.immersion >= 3texp, detect_occultations(lc; method=:threshold))
    detections[obs.name] = found
    ## first pass: the normal velocities are not known yet; the shadow velocity is used
    cs = time_drops(lc, found, τ -> ev.ca.vel; σt=σclock)
    push!(ev, isempty(cs) ? Chord(obs.name, obs, lc) : Chord(obs.name, obs; contacts=cs, lightcurve=lc))
    @printf("%-4s %7.2f° %6.2f°  %-21s %.2f s, σ %.2f, timing error %+.2f s  cross-track %+5.0f km: %d drop(s) simulated, %d detected  %s\n",
            obs.name, obs.lon, obs.lat, st.equipment[1], texp, noise, dclock, offsets[j], length(drops), length(found),
            join([@sprintf("[%.2f, %.2f] s", x.immersion, x.emersion) for x in found], " "))
end
display(ev)

section("3. a single ellipse")
# A convex limb gives at most one drop per chord: the second drop of a chord cannot be
# reproduced, whatever the ellipse. With the visible segments as constraints (the star seen
# between the drops) the fit is worse still.
r1 = fit_limb(ev; rng=rng)
r1v = fit_limb(ev; visible=:soft, rng=rng)
@printf("ellipse: χ²/dof = %.0f; with the visible segments χ²/dof = %.0f\n", reduced_chi2(r1), reduced_chi2(r1v))

section("4. the contact-binary model")
# Two ellipses along a common axis; the starting values come from the ellipse of all the
# contacts (its major axis gives the axis of the lobes). Unlabelled contacts are matched to
# the nearest crossing of the same kind of the two-lobe outline, and the visible segments of
# the chords are constraints by default.
r2a = fit_limb(ev; model=ContactBinaryLimb(), rng=rng)
@printf("first pass (timings with the shadow velocity): χ²/dof = %.2f, b₁ = %.2f ± %.2f km (truth %.0f)\n",
        reduced_chi2(r2a), r2a[:b1]..., b1)
# The timings of the first pass assume the shadow velocity at every edge, while the time scale
# of the diffraction fringes is set by the velocity normal to the limb, slower at oblique
# edges: the contacts are slightly biased. Second pass: each drop is timed again with the
# normal velocities at its edges predicted by the fitted limb, and the limb is refitted.
(cf, Of) = outlines(r2a)[1]
for c in filter(ispositive, collect(ev.chords))
    tr = track(ev, c)
    off = c.lightcurve.tref - ev.tref                      # light-curve time → time from ev.tref
    cs = time_drops(c.lightcurve, detections[c.name], τ -> vnormal(tr, Of, cf, τ + off); σt=σclks[c.name])
    replace_chord!(ev, Chord(c.name, c.observer; contacts=cs, lightcurve=c.lightcurve))
end
# With these errors the chords still scatter more than expected: the short second drop of a
# chord through the waist is timed less well than its formal error says (a few exposures, a
# local minimum of the timing fit), and a gap narrower than an exposure merges two drops into
# one. The errors of the fit are therefore scaled by √(χ²/dof) (`scale_errors=true`). (A limb
# allowance `σ_limb=:auto` fails here: with so few chords a large allowance lets the two-lobe
# limb wander to a different shape.)
r2 = fit_limb(ev; model=ContactBinaryLimb(), scale_errors=true, rng=rng)
@printf("χ²/dof = %.1f: the errors are scaled by √(χ²/dof)\n", reduced_chi2(r2))
display(r2)
println("fitted vs true:")
row(n, tv) = @printf("  %-22s %9.3f ± %6.3f   (truth %8.3f)\n", n, r2[n]..., tv)
for (n, tv) in ((:a1, a1), (:b1, b1), (:a2, a2), (:b2, b2), (:axis_pa, axis_pa), (:neck, neck),
                (:center_f, center[1]), (:center_g, center[2]))
    row(n, tv)
end
w = 2 * Occultations._waist(a1, b1, x1, a2, b2, x2)          # the true width of the waist
for (n, tv) in ((:equivalent_radius, equivalent_radius(truth)), (:length, max(x1 + a1, x2 + a2) - min(x1 - a1, x2 - a2)),
                (:neck_width, w), (:size_ratio, sqrt(a2 * b2 / (a1 * b1))))
    row(n, tv)
end
# The posterior samples the likelihood with the timing errors as given, not the scaled ones:
# it is narrower than the scaled covariance errors above and does not contain the excess scatter.
post = sample(r2; nsteps=3000, rng=rng)
display(post)
ds = derived_summary(post)
@printf("posterior: D_eq = %.1f −%.1f +%.1f km, neck width %.1f −%.1f +%.1f km (truth %.1f), size ratio %.3f −%.3f +%.3f\n",
        2ds.equivalent_radius.median, 2ds.equivalent_radius.minus, 2ds.equivalent_radius.plus, ds.neck_width.median,
        ds.neck_width.minus, ds.neck_width.plus, w, ds.size_ratio.median, ds.size_ratio.minus, ds.size_ratio.plus)
savefig(cornerplot(post; vars=[:a1, :a2, :axis_pa, :neck],
                   axislabels=Dict(:a1 => "a₁ [km]", :a2 => "a₂ [km]", :axis_pa => "axis PA [°]", :neck => "neck"),
                   truths=Dict(:a1 => a1, :a2 => a2, :axis_pa => axis_pa, :neck => neck), truth_label="simulated truth"),
        joinpath(out, "contact_binary_corner.png"))

section("5. joint fit of the light curves")
# The light curves themselves against the two-lobe limb (diffraction, stellar disc, exposure),
# starting from the fit of the timings. Drops less than ~10 Fresnel scales apart are computed
# together: the light diffracted into a narrow gap is part of the model.
# Each light curve carries the timing error of its station (0.05–0.1 s), which the photometric
# model would otherwise take as exact time stamps: every positive station gets a time offset with
# a Gaussian prior of its timing accuracy.
init = Dict(n => v for (n, v) in zip(r2.names, r2.values))
pos = [c.name for c in ev.chords if ispositive(c)]
priors = Dict(Symbol("dt_", n) => (0.0, σclks[n]) for n in pos)
r3 = fit_limb(ev; model=ContactBinaryLimb(), photometry=true, init=init, fit_dt=pos, priors=priors, starts=1, rng=rng)
@printf("χ²/dof = %.3f for %d light-curve points\n", reduced_chi2(r3), r3.dof + length(freenames(r3)))
for (n, tv) in ((:a1, a1), (:b1, b1), (:a2, a2), (:b2, b2), (:axis_pa, axis_pa), (:neck, neck))
    @printf("  %-8s timings %8.3f ± %5.3f   light curves %8.3f ± %5.3f   (truth %7.3f)\n", n, r2[n]..., r3[n]..., tv)
end
println("time offsets of the joint fit (correction; the timing error has the opposite sign):")
for n in pos
    @printf("  %-4s %+6.3f ± %5.3f s   (timing error %+6.3f s)\n", n, r3[Symbol("dt_", n)]..., dclocks[n])
end

# **The joint fit is overconfident here.** It reaches χ²/dof ≈ 1, but with nine stations, of
# which five see the body, its formal errors (a few 0.01 km) are far smaller than its deviations
# from the truth (up to 10 km in b₁), and the fitted time offsets do not recover the timing
# errors. The fit of the timings with scaled errors is the reliable result of this campaign; a
# joint fit needs more, and better sampled, light curves.

section("Figures")
savefig(plot(ev; limb=truth, center=center, title="SimCB: chords and the true limb", size=(1000, 750)),
        joinpath(out, "contact_binary_chords.png"))
savefig(plot(r2; title="SimCB: ContactBinaryLimb fit", size=(1000, 750)), joinpath(out, "contact_binary_fit.png"))
savefig(plot(post; title="SimCB: posterior limbs", size=(1000, 750)), joinpath(out, "contact_binary_posterior.png"))
savefig(plot(r1; title="SimCB: a single ellipse", size=(1000, 750)), joinpath(out, "contact_binary_ellipse.png"))
## a light curve with two drops, and the model of the joint fit
two = [c.name for c in ev.chords if count(x -> x.kind === :immersion, c.contacts) == 2]
for name in two[1:min(2, end)]
    lc = ev[name].lightcurve
    dd = detections[name]
    savefig(plot(lc; model=model_lightcurve(r3, name)[2], title="$name: two drops", xguide="time [s]",
                 xlims=(dd[1].immersion - 5, dd[end].emersion + 5)),
            joinpath(out, "contact_binary_lc_$name.png"))
end
println("figures written to ", out)

# The chords with the true limb, and the fit of the two lobes with its 1σ shape band:
#
#md # ![chords](../assets/figures/contact_binary_chords.png)
#md # ![fit](../assets/figures/contact_binary_fit.png)
#
# The best single ellipse, and the posterior of the lobes:
#
#md # ![ellipse](../assets/figures/contact_binary_ellipse.png)
#md # ![corner plot](../assets/figures/contact_binary_corner.png)
