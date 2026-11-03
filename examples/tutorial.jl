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
ENV["GKSwstype"] = "100"   #src

out = length(ARGS) >= 1 ? ARGS[1] : mkpath(joinpath(@__DIR__, "output", "tutorial"))
rng = MersenneTwister(42)
section(s) = println("\n", "="^78, "\n", s, "\n", "="^78)   #src

section("1. light curves → timings → chords")
eph = EphemKeplerian(a=43.0, e=0.05, i=12.0, Ω=80.0, ω=150.0, M=10.0, epoch=Epoch("2023-01-01"), name="SimTNO")
body = Body("SimTNO", eph; radius=300.0, pole=(30.0, 10.0), rotation_period=7.0, W0=40.0,
            shape=Ellipsoid(340.0, 290.0, 250.0))
t0 = Epoch("2023-05-01T03:00:00")
ra, dec, _ = ephem_radec(body, t0)
star = Star(ra + 0.002 / 3600 / cosd(dec), dec; name="SimStar", mag=Dict("V" => 14.0, "K" => 12.5),
            ra_err=0.15, dec_err=0.12)
ev = Event(star, body, t0)
display(ev)
path = shadow_path(star, body, ev.tref, ev.tref; offsets=[0.0])
lon0, lat0 = path.lon[1, 1], path.lat[1, 1]
truth = project(body, ev.tref)                      # the true apparent limb
center = (25.0, -12.0)                              # ephemeris offset [km]
println("true limb: ", truth, ", centre ", center)

t = collect(-45:0.25:45.0)
for (j, dlat) in enumerate([-3.2, -2.3, -1.0, 0.4, 1.7, 2.45, 3.6])
    obs = Observer("Station $j", lon0 + 0.3j, lat0 + dlat, 200.0)
    pc = predicted_contacts(ev, obs, truth; center=center)
    if isempty(pc)
        lc = LightCurve(t, 1 .+ 0.05 .* randn(rng, length(t)); tref=ev.tref, σ=0.05, exptime=0.25, name=obs.name)
        push!(ev, Chord(obs.name, obs, lc))
        continue
    end
    ti, te = pc[1][1] - ev.tref, pc[1][2] - ev.tref
    j == 2 && ((ti, te) = (ti + 0.9, te + 0.9))      # station 2: clock 0.9 s late
    v = hypot(velocity(Track(star, body, obs, ev.tref, -60, 60), ti)...)
    lc = simulate_lightcurve(t, ti, te; tref=ev.tref, vel=v, dist=ev.ca.dist_km,
                             star_diam=apparent_diameter(star, ev.ca.dist_km), σ=0.05, exptime=0.25, rng=rng,
                             name=obs.name)
    tf = fit_timing(lc; vel=v, dist=ev.ca.dist_km, star_diam=apparent_diameter(star, ev.ca.dist_km), rng=rng)
    c = contact_times(tf)
    @printf("%-10s immersion %s ± %.3f s   emersion %s ± %.3f s\n", obs.name, isot(c.immersion; digits=3),
            c.immersion_err, isot(c.emersion; digits=3), c.emersion_err)
    j == 3 && savefig(plot(tf), joinpath(out, "timing_fit.png"))
    push!(ev, Chord(obs.name, obs, tf))
end

section("2. limb fit: ellipse, negative chords, clock offset, uncertainties")
fitE = fit_limb(ev; fit_dt=["Station 2"], rng=rng)
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
savefig(cornerplot(post; vars=[:equatorial_radius, :oblateness, :position_angle, Symbol("dt_Station 2")],
                   axislabels=Dict(:equatorial_radius => "a [km]", :oblateness => "ε", :position_angle => "PA [°]",
                                   Symbol("dt_Station 2") => "clock offset [s]"),
                   truths=Dict(:equatorial_radius => truth.a, :oblateness => 1 - truth.b / truth.a,
                               :position_angle => truth.pa, Symbol("dt_Station 2") => -0.9), truth_label="simulated truth"),
        joinpath(out, "corner.png"))

# A joint photometric fit: the light curves themselves against the limb.
fitP = fit_limb(ev; photometry=true, fit_dt=["Station 2"], rng=rng, starts=3)
println("\njoint photometric fit: χ²/dof = ", round(reduced_chi2(fitP); digits=3),
        ", a = ", fitP[:equatorial_radius], ", ε = ", fitP[:oblateness])

section("3. a published ellipse against the negative chords")
lit = ellipse_outline(; equatorial_radius=360.0, oblateness=0.1, position_angle=120.0)
for r in check_negative_chords(ev, lit; center=(fitE[:center_f][1], fitE[:center_g][1]))
    @printf("%-10s %s  min distance %7.2f km%s\n", r.chord, r.conflict ? "CONFLICT" : "ok      ", r.distance,
            r.conflict ? @sprintf("  (would have seen %.2f s)", r.duration) : "")
end

section("4. three events of a rotating triaxial body: 3-D ellipsoid and mesh")
evs = Event[]
for (k, dt) in enumerate([0.0, 2.2 * 3600, 4.9 * 3600])
    tk = t0 + dt
    rak, deck, _ = ephem_radec(body, tk)
    st = Star(rak + 0.001 / 3600 / cosd(deck), deck; name="star $k")
    e = Event(st, body, tk; name="event $k")
    p = shadow_path(st, body, e.tref, e.tref; offsets=[0.0])
    Ok = project(body, e.tref)
    for (j, dlat) in enumerate([-2.0, -0.7, 0.6, 1.9])
        o = Observer("E$k-S$j", p.lon[1, 1], p.lat[1, 1] + dlat, 0.0)
        pc = predicted_contacts(e, o, Ok; center=(5.0k, -3.0k))
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
tA = Epoch("2024-01-10T15:00:00")
raA, decA, _ = ephem_radec(bodyA, tA)
starA = Star(raA + 0.001 / 3600 / cosd(decA), decA; name="StarA", diameter=0.0032)
evA = Event(starA, bodyA, tA)
pA = shadow_path(starA, bodyA, evA.tref, evA.tref; offsets=[0.0])
for (j, dlat) in enumerate([0.2, 1.0, 2.3])
    tt = collect(-30:0.5:30.0)
    push!(evA, Chord("A$j", Observer("A$j", pA.lon[1, 1], pA.lat[1, 1] + dlat, 0.0),
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
savefig(plot(lc1; model=model_lightcurve(fitA, "A1")[2]), joinpath(out, "atmosphere_lightcurve.png"))

section("6. maps")
## the standard figures: limb fits (default, 68.3 % region, as the paper), posterior, maps (prediction,
## post-fit, comparison) and the globe
standard_figures("simtno", fitE; out=out, posterior=post, 
                 projection="+proj=laea +lon_0=$(round(lon0)) +lat_0=$(round(lat0))", limits=(lon0 - 30, lon0 + 30, lat0 - 20, lat0 + 20), window=600, marks=60)
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
