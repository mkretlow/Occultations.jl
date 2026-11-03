# # Tutorial: an offline tour with synthetic data
#
# A tour of the package on simulated data. It needs no network and no data files, and prints
# the expected results next to the fitted ones:
#
# 1. light curves → timings (diffraction model) → chords;
# 2. a limb fit (2-D ellipse) with negative chords and a clock offset; uncertainties by
#    covariance, χ² profile and Bayesian posterior; astrometry; a joint fit of the light curves;
# 3. a published ellipse checked against the negative chords;
# 4. three events of a rotating triaxial body: a 3-D ellipsoid fit with a rotational-amplitude
#    constraint, and the same with a triangle mesh (OBJ);
# 5. an atmosphere (2002 XV93-like): a joint fit of the light curves for the surface pressure;
# 6. maps (GeoMakie).
#
# Run it with:
#
#     julia --project=. examples/tutorial.jl [output_dir]      (default: examples/output/tutorial)

using Occultations, Random, Printf, LinearAlgebra
using Plots
import CairoMakie                     # saves the GeoMakie maps (Makie figures)
include(joinpath(@__DIR__, "standard_figures.jl"))
include(joinpath(@__DIR__, "sim_stations.jl"))       # stations on land, in the dark, with their equipment
ENV["GKSwstype"] = "100"   #src

out = length(ARGS) >= 1 ? ARGS[1] : mkpath(joinpath(@__DIR__, "output", "tutorial"))
rng = MersenneTwister(42)
section(s) = println("\n", "="^78, "\n", s, "\n", "="^78)   #src

section("1. light curves → timings → chords")
eph = EphemKeplerian(a=43.0, e=0.05, i=12.0, Ω=80.0, ω=150.0, M=10.0, epoch=Epoch("2023-01-01"), name="SimTNO")
body = Body("SimTNO", eph; radius=300.0, pole=(30.0, 10.0), rotation_period=7.0, W0=40.0,
            shape=Ellipsoid(340.0, 290.0, 250.0))
t0 = Epoch("2023-05-01T05:00:00")
ra, dec, _ = ephem_radec(body, t0)
star = Star(ra + 0.002 / 3600 / cosd(dec), dec; name="SimStar", mag=Dict("V" => 14.0, "K" => 12.5),
            ra_err=0.15, dec_err=0.12)
ev = Event(star, body, t0; name="SimTNO")              # simulated: no date in the plots
display(ev)
path = shadow_path(star, body, ev.tref, ev.tref; offsets=[0.0])
lon0, lat0 = path.lon[1, 1], path.lat[1, 1]
truth = project(body, ev.tref)                      # the true apparent limb
center = (25.0, -12.0)                              # ephemeris offset [km]
println("true limb: ", truth, ", centre ", center)

# Seven stations at random offsets across the path (km from the body's centre), in central South
# America, spread along the track by up to ±15 s, each with its own camera
# (`examples/sim_stations.jl`). The second station that sees the occultation has a clock 0.9 s late:
pv = positionv(star, eph, GEOCENTER, ev.tref)
n̂ = (-pv.vg, pv.vf) ./ hypot(pv.vf, pv.vg)          # cross-track unit vector
offsets = random_offsets(MersenneTwister(7), 7, -400.0, 400.0; mingap=20.0)
sites = stations(star, body, ev.tref, offsets .+ (center[1] * n̂[1] + center[2] * n̂[2]); rng=rng, spread=30.0,
                 prefix="Station ", equipment=EQUIPMENT[1:3])
npos, clocked = 0, ""
for (j, st) in enumerate(sites)
    obs = st.observer
    _, _, texp, noise = st.equipment
    t = collect(round(st.τ) - 45:texp:round(st.τ) + 45)
    pc = predicted_contacts(ev, obs, truth; center=center)
    if isempty(pc)
        lc = LightCurve(t, 1 .+ noise .* randn(rng, length(t)); tref=ev.tref, σ=noise, exptime=texp, name=obs.name)
        push!(ev, Chord(obs.name, obs, lc))
        continue
    end
    ti, te = pc[1][1] - ev.tref, pc[1][2] - ev.tref
    global npos += 1
    if npos == 2                                     # this station's clock is 0.9 s late
        ti, te = ti + 0.9, te + 0.9
        global clocked = obs.name
    end
    v = hypot(velocity(Track(star, body, obs, ev.tref, -60, 60), ti)...)
    lc = simulate_lightcurve(t, ti, te; tref=ev.tref, vel=v, dist=ev.ca.dist_km,
                             star_diam=apparent_diameter(star, ev.ca.dist_km), σ=noise, exptime=texp, rng=rng,
                             name=obs.name)
    tf = fit_timing(lc; vel=v, dist=ev.ca.dist_km, star_diam=apparent_diameter(star, ev.ca.dist_km), rng=rng)
    c = contact_times(tf)
    @printf("%-10s immersion %s ± %.3f s   emersion %s ± %.3f s\n", obs.name, isot(c.immersion; digits=3),
            c.immersion_err, isot(c.emersion; digits=3), c.emersion_err)
    npos == 1 && savefig(plot(tf; xguide="time [s]"), joinpath(out, "timing_fit.png"))
    push!(ev, Chord(obs.name, obs, tf))
end

section("2. limb fit: ellipse, negative chords, clock offset, uncertainties")
println("clock 0.9 s late: ", clocked)
fitE = fit_limb(ev; fit_dt=[clocked], rng=rng)
display(fitE)
println("\ncontact residuals [s]: ", [round(r.dt; digits=4) for r in contact_residuals(fitE)])
println("clock offsets: ", clock_offsets(fitE))
println("χ² profile 1σ interval of the equatorial radius: ", confidence_interval(fitE, :equatorial_radius))
a = astrometry(fitE)
@printf("astrometry: RA %s  Dec %s  (σ %.3f, %.3f mas)  offsets %+.3f, %+.3f mas\n", a.ra_str, a.dec_str,
        a.ra_err, a.dec_err, a.offset_ra, a.offset_dec)
post = sample(fitE; nsteps=2000, rng=rng)
display(post)
println("derived (posterior): ", derived_summary(post))
println("autocorrelation times [steps]: ", round.(autocorrelation_time(post); digits=1))
savefig(cornerplot(post; vars=[:equatorial_radius, :oblateness, :position_angle, Symbol("dt_", clocked)],
                   axislabels=Dict(:equatorial_radius => "a [km]", :oblateness => "ε", :position_angle => "PA [°]",
                                   Symbol("dt_", clocked) => "clock offset [s]"),
                   truths=Dict(:equatorial_radius => truth.a, :oblateness => 1 - truth.b / truth.a,
                               :position_angle => truth.pa, Symbol("dt_", clocked) => -0.9), truth_label="simulated truth"),
        joinpath(out, "corner.png"))

# A joint photometric fit: the light curves themselves against the limb.
fitP = fit_limb(ev; photometry=true, fit_dt=[clocked], rng=rng, starts=3)
println("\njoint photometric fit: χ²/dof = ", round(reduced_chi2(fitP); digits=3),
        ", a = ", fitP[:equatorial_radius], ", ε = ", fitP[:oblateness])

section("3. a published ellipse against the negative chords")
lit = ellipse_outline(; equatorial_radius=385.0, oblateness=0.15, position_angle=67.0)   # the truth, 16 % larger
for r in check_negative_chords(ev, lit; center=(fitE[:center_f][1], fitE[:center_g][1]))
    @printf("%-10s %s  min distance %7.2f km%s\n", r.chord, r.conflict ? "CONFLICT" : "ok      ", r.distance,
            r.conflict ? @sprintf("  (would have seen %.2f s)", r.duration) : "")
end

section("4. three events of a rotating triaxial body: 3-D ellipsoid and mesh")
evs = Event[]
## the same hour on three nights (P = 7 h): rotational phases 0, 0.43 and 0.86
for (k, dt) in enumerate([0.0, 1.0, 2.0] .* 86400)
    tk = t0 + dt
    rak, deck, _ = ephem_radec(body, tk)
    st = Star(rak + 0.001 / 3600 / cosd(deck), deck; name="star $k")
    e = Event(st, body, tk; name="event $k")
    Ok = project(body, e.tref)
    ck = (5.0k, -3.0k)
    pk = positionv(st, eph, GEOCENTER, e.tref)
    nk = (-pk.vg, pk.vf) ./ hypot(pk.vf, pk.vg)
    offs = random_offsets(MersenneTwister(10 + k), 5, -280.0, 280.0; mingap=20.0)
    for s in stations(st, body, e.tref, offs .+ (ck[1] * nk[1] + ck[2] * nk[2]); rng=rng, prefix="E$k-S")
        o = s.observer
        pc = predicted_contacts(e, o, Ok; center=ck)
        isempty(pc) && continue
        push!(e, Chord(o.name, o; immersion=pc[1][1] + 0.03randn(rng), emersion=pc[1][2] + 0.03randn(rng),
                       immersion_err=0.03, emersion_err=0.03))
    end
    push!(evs, e)
end
u = normalize(geocentric_vector(eph, t0 + 30 * 86400))
amp = RotationalAmplitude(t0 + 30 * 86400, rotational_amplitude(340.0, 290.0, 250.0, aspect_angle(30.0, 10.0, u)), 0.02)
body.shape = Ellipsoid(310.0, 300.0, 280.0)            # deliberately poor starting shape
fit3 = fit(LimbProblem(evs; model=EllipsoidLimb(), amplitudes=[amp], fix=[:pole_ra, :pole_dec]); rng=rng)
display(fit3)
println("(true axes 340, 290, 250 km)")
write_obj(joinpath(out, "simtno_fitted_ellipsoid.obj"),
          ellipsoid_mesh(fit3[:a][1], fit3[:b][1], fit3[:c][1]; name="SimTNO fitted ellipsoid"))
for k in 1:3
    savefig(plot(fit3; event=k), joinpath(out, "ellipsoid_event$k.png"))
end
# The 3-D view of each event. With GLMakie, the event, the rotational phase and (with
# `explore=true`) the pole, W0 and the axes are sliders, with the live χ²:
for k in 1:3
    CairoMakie.save(joinpath(out, "ellipsoid_occview_event$k.png"), occview(fit3; event=k))
end
# The same body as a triangle mesh, written to and read from an OBJ file:
objfile = joinpath(out, "simtno.obj")
write_obj(objfile, ellipsoid_mesh(340.0, 290.0, 250.0; subdivisions=4))
mesh = read_obj(objfile)
fitM = fit_limb(evs; model=MeshLimb(mesh; fit_orientation=true), fix=[:pole_ra, :pole_dec], rng=rng, starts=4)
println("\nmesh fit: scale = ", fitM[:scale], ", W0 = ", fitM[:W0], ", volume-equivalent radius = ",
        fitM[:volume_equivalent_radius])

section("5. atmosphere (2002 XV93-like, Arimatsu et al. 2026)")
ephA = EphemKeplerian(a=37.5, e=0.02, i=12.0, Ω=80.0, ω=150.0, M=10.0, epoch=Epoch("2024-01-01"), name="XV93-like")
bodyA = Body("XV93-like", ephA; radius=244.0, density=1500.0)
tA = Epoch("2024-05-20T04:00:00")
raA, decA, _ = ephem_radec(bodyA, tA)
starA = Star(raA + 0.001 / 3600 / cosd(decA), decA; name="StarA", diameter=0.0032)
evA = Event(starA, bodyA, tA; name="XV93-like")
## a central chord, one halfway to the limb and a near miss through the atmosphere
pvA = positionv(starA, ephA, GEOCENTER, evA.tref)
nA = (-pvA.vg, pvA.vf) ./ hypot(pvA.vf, pvA.vg)
for (j, s) in enumerate(stations(starA, bodyA, evA.tref, [20.0, 110.0, 250.0] .+ (15.0 * nA[1] - 20.0 * nA[2]); rng=rng,
                                 prefix="A"))
    tt = collect(round(s.τ) - 30:0.5:round(s.τ) + 30)
    push!(evA, Chord("A$j", s.observer,
                     LightCurve(tt, ones(length(tt)); tref=evA.tref, σ=0.03, exptime=0.5, λ=0.6, name="A$j")))
end
atm = xv93_atmosphere(:CH4)
sim = LimbProblem(evA; model=CircleLimb(), photometry=true, atmosphere=atm,
                  init=Dict(:center_f => 15.0, :center_g => -20.0, :radius => 244.0, :pressure => 124.0))
v = [q.value for q in parameters(sim)]
rngA = MersenneTwister(3)
for d in sim.photometry                                  # synthetic light curves from the model
    I, _ = Occultations._phot_model(sim, d, Occultations._outlines(sim, v)[1], v)
    d.lc.flux .= I .+ 0.03 .* randn(rngA, length(I))
end
airless = fit(LimbProblem(evA; model=CircleLimb(), photometry=true,
                          init=Dict(:center_f => 0.0, :center_g => 0.0, :radius => 230.0)); starts=3, rng=rng)
fitA = fit(LimbProblem(evA; model=CircleLimb(), photometry=true, atmosphere=atm,
                       init=Dict(:center_f => airless[:center_f][1], :center_g => airless[:center_g][1],
                                 :radius => airless[:radius][1], :pressure => 50.0)); starts=6, rng=rng)
display(fitA)
χtruth = sum(abs2, residuals(fitA.problem, [15.0, -20.0, 244.0, 124.0]))
@printf("airless model: χ² = %.1f; with atmosphere: χ² = %.1f (at the simulated truth, 124 nbar and R = 244 km: %.1f)\n",
        airless.chi2, fitA.chi2, χtruth)
prof = atmosphere_profile(atm, fitA[:pressure][1], fitA[:radius][1], body_gm(bodyA), evA.ca.dist_km;
                          λ=0.6, star_diam=apparent_diameter(starA, evA.ca.dist_km))
savefig(plot(prof; ρmax=400), joinpath(out, "atmosphere_profile.png"))
lc1 = evA["A1"].lightcurve
savefig(plot(lc1; model=model_lightcurve(fitA, "A1")[2], xguide="time [s]"), joinpath(out, "atmosphere_lightcurve.png"))

section("6. maps")
## the standard figures: limb fits (default, 68.3 % region, as the paper), posterior, maps (prediction,
## post-fit, comparison) and the globe
standard_figures("simtno", fitE; out=out, posterior=post, 
                 projection="+proj=laea +lon_0=$(round(lon0)) +lat_0=$(round(lat0))", limits=(lon0 - 30, lon0 + 30, lat0 - 20, lat0 + 20), window=600, marks=0)
println("\nfigures written to ", out)

# The timing fit of station 3, and the posterior of the ellipse and the clock offset with the
# simulated truth:
#
#md # ![timing fit](../assets/figures/timing_fit.png)
#md # ![corner plot](../assets/figures/corner.png)
#
# The fitted 3-D ellipsoid in the geometry of the first event:
#
#md # ![occview](../assets/figures/ellipsoid_occview_event1.png)
#
# The refracted flux of the fitted atmosphere, and a light curve with the model:
#
#md # ![atmosphere profile](../assets/figures/atmosphere_profile.png)
#md # ![atmosphere light curve](../assets/figures/atmosphere_lightcurve.png)
