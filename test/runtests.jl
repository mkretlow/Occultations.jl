using Occultations
using Occultations: _fg_at, K_BOLTZMANN, AMU, AU_KM, _refracted_flux, _refractivity_K
using Test, Random, LinearAlgebra, Statistics
ENV["GKSwstype"] = "100"
using Plots

# --- common synthetic scene -------------------------------------------------------------
const EPH = EphemKeplerian(a=43.0, e=0.05, i=12.0, Ω=80.0, ω=150.0, M=10.0, epoch=Epoch("2023-01-01"), name="SimTNO")
const T0 = Epoch("2023-05-01T03:00:00")

function scene(; radius=300.0, kw...)
    body = Body("SimTNO", EPH; radius=radius, kw...)
    ra, dec, _ = ephem_radec(body, T0)
    star = Star(ra + 0.002 / 3600 / cosd(dec), dec; name="SimStar", ra_err=0.2, dec_err=0.2)
    ev = Event(star, body, T0)
    p = shadow_path(star, body, ev.tref, ev.tref; offsets=[0.0])
    return body, star, ev, p.lon[1, 1], p.lat[1, 1]
end

# chords from a known outline, with Gaussian timing noise σt
function add_chords!(ev, O, center, lon0, lat0, dlats; σt=0.05, rng, shifts=Dict{Int,Float64}(), prefix="S")
    for (j, dlat) in enumerate(dlats)
        o = Observer("$prefix$j", lon0, lat0 + dlat, 0.0)
        pc = predicted_contacts(ev, o, O; center=center)
        if isempty(pc)
            push!(ev, Chord("$prefix$j", o; span=(ev.tref - 150, ev.tref + 150)))
        else
            s = get(shifts, j, 0.0)
            push!(ev, Chord("$prefix$j", o; immersion=pc[1][1] + σt * randn(rng) + s, emersion=pc[1][2] + σt * randn(rng) + s,
                            immersion_err=σt, emersion_err=σt))
        end
    end
    return ev
end

@testset "Occultations.jl" begin

@testset "time, angles, star" begin
    e = Epoch("2020-06-03T12:30:45.123456")
    @test isot(e; digits=6) == "2020-06-03T12:30:45.123456"
    @test (e + 0.25) - e ≈ 0.25 atol = 1e-6
    @test Epoch("2017-01-01T00:00:00") - Epoch("2016-12-31T23:59:59") ≈ 2.0 atol = 1e-6
    @test parse_ra("12 30 00") ≈ 187.5
    @test parse_dec("-12 30 00") ≈ -12.5
    s = Star(10.0, 20.0; pmra=100.0, pmdec=-50.0, parallax=10.0, epoch=2016.0)
    ra, dec, _ = barycentric(s, Epoch("2026-01-01"))
    @test (ra - 10) * cosd(20) * 3.6e6 ≈ 1000 rtol = 0.01
    @test (dec - 20) * 3.6e6 ≈ -500 rtol = 0.01
end

@testset "outlines and 3-D projection" begin
    O = ellipse_outline(; equatorial_radius=300.0, oblateness=0.2, position_angle=35.0)
    @test O.b ≈ 240 && area(O) ≈ π * 300 * 240
    # minor axis along PA 35°: point (sin35, cos35)·240 is on the limb
    @test limb_radius(O, 240sind(35), 240cosd(35)) ≈ 240 atol = 1e-9
    # analytic and numerical crossings agree
    P = project(ellipsoid_mesh(300.0, 300.0, 240.0; subdivisions=5), Matrix(1.0I, 3, 3), [1.0, 0, 0])
    for entry in (true, false)
        s1, m1 = crossing(O, -500.0, 30.0, 20.0, 1.0, entry)
        s2, m2 = crossing(Occultations.PolarOutline(collect(range(-π, π; length=2001)[1:2000]),
                                                     [limb_radius(O, t) for t in range(-π, π; length=2001)[1:2000]]),
                          -500.0, 30.0, 20.0, 1.0, entry)
        @test s1 ≈ s2 atol = 1e-3
        @test m1 == 0 && m2 == 0
    end
    s, miss = crossing(O, -500.0, 400.0, 20.0, 0.0, true)
    @test miss > 0
    θ = range(0, 2π; length=12)[1:end-1]
    u = 300 .* cos.(θ); w = 240 .* sin.(θ)
    ae = algebraic_ellipse(10 .+ u .* cosd(35) .+ w .* sind(35), -5 .- u .* sind(35) .+ w .* cosd(35))
    @test collect(values(ae)) ≈ [10.0, -5.0, 300.0, 0.2, 35.0] atol = 1e-8
    @test algebraic_ellipse([0.0, 1, 2, 3, 4], [0.0, 1, 2, 3, 4]) === nothing
    # exact projection of an ellipsoid vs Dunham et al. (2019) and Magnusson (1986)
    a, b, c = 1161.0, 852.0, 513.0                      # Haumea
    @test collect(projected_axes_dunham(a, b, c, 0.0, 76.3)) ≈ [852.0, 569.0] atol = 1.0
    for (phase, aspect) in ((0.0, 76.3), (30.0, 60.0), (75.0, 40.0), (120.0, 85.0))
        # line of sight in the body frame: aspect from the pole, phase from the long axis (minimum light)
        ub = [sind(aspect) * cosd(phase), sind(aspect) * sind(phase), cosd(aspect)]
        Ob = project(Ellipsoid(a, b, c), Matrix(1.0I, 3, 3), ub)
        @test collect(projected_axes_dunham(a, b, c, phase, aspect)) ≈ [Ob.a, Ob.b] rtol = 1e-9
        @test collect(projected_axes_magnusson(a, b, c, phase + 90, aspect)) ≈ [Ob.a, Ob.b] rtol = 1e-9
        Om = project(ellipsoid_mesh(a, b, c; subdivisions=5), Matrix(1.0I, 3, 3), ub)
        @test area(Om) ≈ area(Ob) rtol = 0.005
    end
    # rotational amplitude: closed formula vs numerical projected areas
    u = normalize([0.3, -0.5, 0.2])
    for (α, δ) in ((30.0, 10.0), (100.0, -40.0))
        A1 = rotational_amplitude(340.0, 290.0, 250.0, aspect_angle(α, δ, u))
        A2 = rotational_amplitude(Ellipsoid(340.0, 290.0, 250.0), α, δ, u; nphase=360)
        @test A1 ≈ A2 atol = 1e-3
    end
end

@testset "meshes and OBJ files" begin
    m = ellipsoid_mesh(30.0, 20.0, 10.0; subdivisions=4)
    @test volume(m) ≈ 4π / 3 * 6000 rtol = 0.01
    f = tempname() * ".obj"
    write_obj(f, m)
    m2 = read_obj(f)
    @test size(m2.faces) == size(m.faces)
    @test maximum(abs.(m2.vertices .- m.vertices)) < 1e-5
    @test volume_equivalent_radius(Occultations.scale(m, 2.0)) ≈ 2volume_equivalent_radius(m)
    @test ellipsoid_mesh(1.0; subdivisions=0, name="unit").name == "unit"
end

@testset "automatic differentiation" begin
    include(joinpath(@__DIR__, "..", "examples", "data", "2003uy117_20201023.jl"))
    ev = Base.invokelatest(uy117_event)
    for prob in (LimbProblem(ev), LimbProblem(ev; residual=:radial), LimbProblem(ev; model=CircleLimb()))
        r = fit(prob; rng=MersenneTwister(1), ad=false)
        ps = parameters(prob); F = Occultations._Free(ps)
        x = r.values[F.idx]; sc = [ps[i].scale for i in F.idx]
        res = y -> Occultations._total_residuals(prob, F, y)
        ad = Ref(true)
        Ja = Occultations._jacobian(res, x, res(x), sc; ad=ad)
        Jf = Occultations._jacobian(res, x, res(x), sc; ad=Ref(false))
        @test ad[]                                   # AD was used, not the fallback
        @test maximum(abs.(Ja .- Jf)) < 1e-6 * maximum(abs.(Jf))
        @test fit(prob; rng=MersenneTwister(1)).chi2 ≈ r.chi2 rtol = 1e-6
    end
end

@testset "limb bands" begin
    rng = MersenneTwister(7)
    θ = range(0, 2π; length=361)
    draws = [(c = (2randn(rng), 2randn(rng)); r = 100 + 3randn(rng);
              (c[1] .+ r .* sin.(θ), c[2] .+ r .* cos.(θ))) for _ in 1:1000]
    for (sim, lo_, hi_) in ((false, 0.0, 0.5), (true, 0.64, 0.72))
        x, y = Occultations._limb_band(draws, (0.0, 0.0); simultaneous=sim)
        hi = hypot.(x[1:361], y[1:361]); lo = reverse(hypot.(x[362:end], y[362:end]))
        f = count(d -> (ρ = Occultations._polar_radius(d[1], d[2], θ); all(lo .- 1e-6 .<= ρ .<= hi .+ 1e-6)), draws) / 1000
        @test lo_ <= f <= hi_
    end
end

@testset "light-curve readers" begin
    d = mktempdir()
    # PyMovie: date in the header, crossing of 0h UT
    f = joinpath(d, "pymovie.csv")
    write(f, """# PyMovie Version 3.8
# date at frame 0: 2024-01-10
FrameNum,timeInfo,signal-target,appsum-target
0,[23:59:59.9000],100.0,120.0
1,[23:59:59.9400],101.0,121.0
2,[00:00:00.0200],50.0,70.0
""")
    lc = read_pymovie(f)
    @test length(lc) == 3 && lc.flux[3] == 50.0
    @test (lc.tref + lc.time[3]) - Epoch("2024-01-11T00:00:00.02") ≈ 0 atol = 1e-6
    @test read_lightcurve(f; format=:pymovie).flux == lc.flux
    # Tangra
    f = joinpath(d, "tangra.csv")
    write(f, """Tangra v3.7
FrameNo,Time (UT),Signal (1),Background (1),SignalMinusBackground (1)
10,[13:12:50.100],210.0,10.0,200.0
11,[13:12:50.140],110.0,10.0,100.0
""")
    lc = read_tangra(f; date="2024-01-10")
    @test lc.flux == [200.0, 100.0]
    @test (lc.tref + lc.time[2]) - Epoch("2024-01-10T13:12:50.14") ≈ 0 atol = 1e-6
    # AstroImageJ, tab-separated
    f = joinpath(d, "aij.tbl")
    write(f, " \tLabel\tJD_UTC\trel_flux_T1\trel_flux_err_T1\tSource-Sky_T1\n" *
             "1\tim1\t2460320.050000000\t1.01\t0.02\t5000\n" *
             "2\tim2\t2460320.050011574\t0.20\t0.02\t1000\n")
    lc = read_astroimagej(f)
    @test lc.flux == [1.01, 0.20] && lc.σ == [0.02, 0.02]
    @test lc.time[2] - lc.time[1] ≈ 1.0 atol = 1e-3
    @test read_astroimagej(f; column="Source-Sky_T1").flux == [5000.0, 1000.0]
    @test read_lightcurve(f; format=:astroimagej).flux == [1.01, 0.20]
end

@testset "Occultation alias" begin
    @test Occultation === Event
end

@testset "tracks" begin
    body, star, ev, lon0, lat0 = scene()
    obs = Observer("T", lon0, lat0 + 1, 300)
    tr = Track(star, body, obs, ev.tref, -300, 900)
    err = maximum(begin
                      r = _fg_at(star, EPH, obs, ev.tref + τ)
                      f, g = tr(τ)
                      hypot(r[1] - f, r[2] - g)
                  end for τ in -300:7.3:900)
    @test err < 0.002                          # km
    p = positionv(star, body, obs, ev.tref + 100)
    @test collect(velocity(tr, 100.0)) ≈ [p.vf, p.vg] atol = 1e-4
end

@testset "light curves and timing" begin
    rng = MersenneTwister(1)
    lc = simulate_lightcurve(collect(-30:0.2:30.0), -5.0, 7.0; tref=T0, σ=0.03, vel=20.0, dist=6e9, rng=rng,
                             flux_top=1.07, flux_bottom=0.15)
    d = detect_occultation(lc)
    @test d.immersion ≈ -5 atol = 0.3
    @test d.emersion ≈ 7 atol = 0.3
    n = normalize(lc)
    @test median(n.flux[n.time .< -10]) ≈ 1 atol = 0.01
    r = fit_timing(lc; vel=20.0, dist=6e9)
    @test r[:immersion][1] ≈ -5 atol = 3r[:immersion][2] + 1e-3
    @test r[:emersion][1] ≈ 7 atol = 3r[:emersion][2] + 1e-3
    @test r[:flux_top][1] ≈ 1.07 atol = 0.01
    @test r[:flux_bottom][1] ≈ 0.15 atol = 0.02
    @test r[:duration][1] ≈ 12 atol = 0.05
    @test 0.5 < reduced_chi2(r) < 1.5
    ct = contact_times(r)
    @test ct.immersion - T0 ≈ r[:immersion][1]
    # single edge
    r1 = fit_timing(select(lc, -30, 2); vel=20.0, dist=6e9, edges=:immersion)
    @test r1[:immersion][1] ≈ -5 atol = 0.05
    @test contact_times(r1).emersion === nothing
    # posterior of the timing
    post = sample(r; nsteps=400, rng=rng)
    s = summarize(post)
    @test s.immersion.median ≈ r[:immersion][1] atol = 2r[:immersion][2]
    @test 0.5 < (s.immersion.minus + s.immersion.plus) / 2 / r[:immersion][2] < 2
    # file round trip
    f = tempname()
    write_lightcurve(f, lc)
    lc2 = read_lightcurve(f; tref=T0, cols=(1, 2, 3))
    @test lc2.flux ≈ lc.flux rtol = 1e-6
    @test lc2.σ ≈ lc.σ rtol = 1e-3
    # magnitude drop
    Δm, fb = magnitude_drop(12.0, 12.0)
    @test fb ≈ 0.5 && Δm ≈ 2.5log10(2)
end

@testset "ellipse limb: clock offset, negatives, uncertainties" begin
    rng = MersenneTwister(11)
    body, star, ev, lon0, lat0 = scene()
    O = ellipse_outline(; equatorial_radius=310.0, oblateness=0.15, position_angle=30.0)
    add_chords!(ev, O, (20.0, -15.0), lon0, lat0, [-3.0, -2.2, -1.0, 0.3, 1.6, 2.35, 2.6, 4.0]; rng=rng,
                shifts=Dict(3 => 0.7))
    @test count(isnegative, ev.chords) == 3 && count(iscomplete, ev.chords) == 5
    r = fit_limb(ev; fit_dt=["S3"], rng=rng)
    truth = (center_f=20.0, center_g=-15.0, equatorial_radius=310.0, oblateness=0.15, position_angle=30.0)
    for k in keys(truth)
        @test r[k][1] ≈ truth[k] atol = 3r[k][2]
    end
    @test r[:dt_S3][1] ≈ -0.7 atol = 3r[:dt_S3][2]
    @test r.converged && r.dof == 10 - 6
    @test clock_offsets(r)["S3"][1] == r[:dt_S3][1]
    @test all(x -> abs(x.dt) < 4 * 0.05, contact_residuals(r))
    @test all(x -> !x.conflict, check_negative_chords(r))
    # astrometry: offsets in mas, star error included
    a = astrometry(r)
    @test a.offset_ra ≈ r[:center_f][1] / ev.ca.dist_km * 206264806.247 rtol = 1e-3
    @test a.ra_err ≥ 0.2
    # χ² profile interval ≈ ±1σ
    lo, hi = confidence_interval(r, :equatorial_radius)
    @test (hi - lo) / 2 ≈ r[:equatorial_radius][2] rtol = 0.25
    # hard negative constraint, radial residuals, circle
    @test fit_limb(ev; fit_dt=["S3"], negative=:hard, rng=rng).chi2 ≈ r.chi2 rtol = 1e-3
    @test fit_limb(ev; fit_dt=["S3"], residual=:radial, rng=rng)[:equatorial_radius][1] ≈ 310 atol = 3
    @test reduced_chi2(fit_limb(ev; model=CircleLimb(), fit_dt=["S3"], rng=rng)) > 10
    # a literature ellipse that is too large crosses negative chords; contacts as predicted
    big = ellipse_outline(; equatorial_radius=420.0)
    chk = check_negative_chords(ev, big; center=(20.0, -15.0))
    bad = filter(x -> x.conflict, chk)
    @test !isempty(bad) && chk[1].distance < 0
    for x in bad
        pc = predicted_contacts(ev, ev[x.chord].observer, big; center=(20.0, -15.0))
        @test pc[1][1] - x.immersion ≈ 0 atol = 1e-3
        @test pc[1][2] - x.emersion ≈ 0 atol = 1e-3
    end
    # disabled contacts, partial chords, priors
    ev2 = deepcopy(ev)
    replace_chord!(ev2, disable(ev2["S4"], :emersion))
    @test length(contact_points(ev2)) == length(contact_points(ev)) - 1
    r2 = fit_limb(ev2; fit_dt=["S3"], priors=Dict(:oblateness => (0.15, 0.001)), rng=rng)
    @test r2[:oblateness][2] < r[:oblateness][2]
    # Bayesian posterior agrees with the covariance
    post = sample(r; nsteps=800, rng=rng)
    s = summarize(post)
    @test s.equatorial_radius.median ≈ r[:equatorial_radius][1] atol = r[:equatorial_radius][2]
    @test (s.equatorial_radius.minus + s.equatorial_radius.plus) / 2 ≈ r[:equatorial_radius][2] rtol = 0.4
    ds = derived_summary(post)
    @test ds.equivalent_radius.median ≈ r[:equivalent_radius][1] atol = 2r[:equivalent_radius][2]
    @test length(autocorrelation_time(post)) == 6
    @test astrometry(post).offset_ra ≈ a.offset_ra atol = 3a.offset_ra_err
    # plots
    @test plot(r) isa Plots.Plot
    @test plot(post) isa Plots.Plot
    @test cornerplot(post; vars=[:equatorial_radius, :oblateness]) isa Plots.Plot
    @test plot(ev; limb=O, center=(20.0, -15.0)) isa Plots.Plot
    m = occ_map(r; window=300)
    @test m.path.offsets[1] == 0 && length(m.path.offsets) == 5
    @test Occultations._offset_kinds(m.path.offsets) == [:centre, :limit, :limit, :error, :error]
    @test geomap(m) isa Occultations.GeoMaps.Figure
    @test geomap(r; window=300) isa Occultations.GeoMaps.Figure
    @test geomap(ev; window=120) isa Occultations.GeoMaps.Figure
    @test geomap(occ_map(ev; window=120, error_mas=10.0), m) isa Occultations.GeoMaps.Figure
end

@testset "timing fits → chords, and joint photometric fit" begin
    rng = MersenneTwister(5)
    body, star, ev, lon0, lat0 = scene()
    O = ellipse_outline(; equatorial_radius=333.0, oblateness=0.15, position_angle=125.0)
    t = collect(-40:0.2:40.0)
    for (j, dlat) in enumerate([-2.5, -1.0, 0.3, 1.6, 2.6, 4.0])
        o = Observer("S$j", lon0, lat0 + dlat, 0.0)
        pc = predicted_contacts(ev, o, O; center=(10.0, -5.0))
        if isempty(pc)
            push!(ev, Chord("S$j", o, LightCurve(t, 1 .+ 0.04 .* randn(rng, length(t)); tref=ev.tref, σ=0.04, exptime=0.2)))
        else
            ti, te = pc[1][1] - ev.tref, pc[1][2] - ev.tref
            tr = Track(star, body, o, ev.tref, -60, 60)
            v = hypot(velocity(tr, ti)...)
            lc = simulate_lightcurve(t, ti, te; tref=ev.tref, vel=v, dist=ev.ca.dist_km, σ=0.04, exptime=0.2, rng=rng)
            push!(ev, Chord("S$j", o, fit_timing(lc; vel=v, dist=ev.ca.dist_km)))
        end
    end
    r = fit_limb(ev; rng=rng)
    @test r[:equatorial_radius][1] ≈ 333 atol = 3r[:equatorial_radius][2] + 0.05
    @test r[:position_angle][1] ≈ 125 atol = 3r[:position_angle][2] + 0.05
    rp = fit_limb(ev; photometry=true, rng=rng, starts=2)
    @test 0.8 < reduced_chi2(rp) < 1.2
    @test rp[:equatorial_radius][1] ≈ 333 atol = 3rp[:equatorial_radius][2] + 0.05
    @test rp[:center_f][1] ≈ 10 atol = 3rp[:center_f][2] + 0.05
    tt, mf = model_lightcurve(rp, "S3")
    @test length(mf) == length(t) && minimum(mf) < 0.1
end

@testset "3-D ellipsoid over several events, rotational amplitude, meshes" begin
    rng = MersenneTwister(5)
    truthE = Ellipsoid(340.0, 290.0, 250.0)
    body = Body("SimTNO", EPH; radius=300.0, pole=(30.0, 10.0), rotation_period=7.0, W0=40.0, shape=truthE)
    evs = Event[]
    for (k, dt) in enumerate([0.0, 2.2 * 3600, 4.9 * 3600])
        tk = T0 + dt
        ra, dec, _ = ephem_radec(body, tk)
        st = Star(ra + 0.001 / 3600 / cosd(dec), dec; name="star$k")
        e = Event(st, body, tk; name="E$k")
        p = shadow_path(st, body, e.tref, e.tref; offsets=[0.0])
        add_chords!(e, project(body, e.tref), (5.0k, -3.0k), p.lon[1, 1], p.lat[1, 1], [-2.0, -0.7, 0.6, 1.9];
                    σt=0.03, rng=rng, prefix="E$(k)S")
        push!(evs, e)
    end
    u = normalize(geocentric_vector(EPH, T0 + 30 * 86400))
    amp = RotationalAmplitude(T0 + 30 * 86400, rotational_amplitude(340.0, 290.0, 250.0, aspect_angle(30.0, 10.0, u)), 0.02)
    body.shape = Ellipsoid(310.0, 300.0, 280.0)                   # poor start
    r = fit(LimbProblem(evs; model=EllipsoidLimb(), amplitudes=[amp], fix=[:pole_ra, :pole_dec]); rng=rng)
    @test r[:a][1] ≥ r[:b][1]                                     # canonical form
    @test r[:a][1] ≈ 340 atol = 3r[:a][2] + 0.5
    @test r[:b][1] ≈ 290 atol = 3r[:b][2] + 0.5
    @test r[:c][1] ≈ 250 atol = 3r[:c][2] + 0.5
    Wtrue = mod(orientation_at(body.frame, evs[1].tref)[3], 180)     # W0 refers to the first event's epoch
    @test mod(r[:W0][1] - Wtrue + 90, 180) - 90 ≈ 0 atol = 3r[:W0][2] + 0.5
    for k in 1:3
        Ok = project(Body("x", EPH; pole=(30.0, 10.0), rotation_period=7.0, W0=40.0, shape=truthE), evs[k].tref)
        @test r[Symbol("apparent_a_$k")][1] ≈ Ok.a atol = 3r[Symbol("apparent_a_$k")][2] + 0.5
    end
    # the 3-D view: the projected mesh of every event agrees with the fitted silhouette
    for k in 1:3
        g = Occultations.OccViews._geometry(r.problem, r.values, k)
        x, y = Occultations.OccViews._silhouette(r.problem, r.values, k)
        @test maximum(g.X[1, :]) ≈ maximum(x) rtol = 0.02
        @test maximum(g.X[2, :]) ≈ maximum(y) rtol = 0.02
        @test minimum(g.X[2, :]) ≈ minimum(y) rtol = 0.02
    end
    @test occview(r; event=2, explore=true) isa Occultations.OccViews.Figure
    # rotational phase since maximum light: model function vs brute force with the full
    # orientation code (light time included), also for a retrograde rotator
    for (Pr, pole, W0) in ((12.376, (337.0, 62.0), 40.0), (-9.3, (120.0, -35.0), 77.0))
        bb = Body("x", EPH; pole=pole, rotation_period=Pr, W0=W0, shape=Ellipsoid(166.0, 108.0, 90.0))
        τs = range(0, abs(Pr) * 3600; length=2001)[1:end-1]
        f_bf = mod(τs[argmax([area(project(bb, T0 - τ)) for τ in τs])] / (abs(Pr) * 3600), 0.5)
        B = geocentric_vector(EPH, T0)
        W = orientation_at(bb.frame, T0 - norm(B) / 299792.458)[3]
        @test mod(rotational_phase(bb.shape, pole..., W, Pr, normalize(B)), 0.5) ≈ f_bf atol = 1e-3
    end
    # the phase constraint selects the rotation phase: refit event 1 alone with W0 free
    ph = RotationalPhase(evs[1].tref, mod(rotational_phase(truthE, 30.0, 10.0, orientation_at(body.frame,
                         evs[1].tref - norm(geocentric_vector(EPH, evs[1].tref)) / 299792.458)[3], 7.0,
                         normalize(geocentric_vector(EPH, evs[1].tref))), 0.5), 0.02)
    body.shape = truthE
    rp = fit(LimbProblem(evs[1]; model=EllipsoidLimb(), phases=[ph], fix=[:pole_ra, :pole_dec, :a, :b, :c]); rng=rng)
    @test abs(residuals(rp.problem, rp.values)[end]) < 0.5
    # per-event 2-D ellipses sharing the axes but not the orientation
    r2 = fit_limb(evs; model=EllipseLimb(per_event=(:position_angle,)), rng=rng)
    @test :position_angle_3 in r2.names && :equatorial_radius in r2.names
    # mesh: the ellipsoid as a triangle mesh, scale and rotation phase
    body.shape = truthE
    mesh = ellipsoid_mesh(truthE; subdivisions=4)
    rm = fit_limb(evs[1]; model=MeshLimb(mesh), init=Dict(:scale => 0.9), rng=rng)
    @test rm[:scale][1] ≈ 1 atol = 0.005
    rm3 = fit_limb(evs; model=MeshLimb(mesh; fit_orientation=true), fix=[:pole_ra, :pole_dec],
                   init=Dict(:W0 => mod(Wtrue + 20, 360)), rng=rng, starts=4)
    @test mod(rm3[:W0][1] - Wtrue + 90, 180) - 90 ≈ 0 atol = 1.5
    @test rm3[:scale][1] ≈ 1 atol = 0.005
end

@testset "atmosphere" begin
    # isothermal: ω(r) ≈ -ν(r) √(2π r / H) for r ≫ H
    R, GM, D = 1200.0, 869.6, 33.0 * AU_KM
    atm = Atmosphere(:N2; T0=40.0)
    p = atmosphere_profile(atm, 10_000.0, R, GM, D)
    μ = 28.0134 * AMU
    for z in (10.0, 30.0, 60.0)
        j = argmin(abs.(p.r .- (R + z)))
        r = p.r[j]
        n = 1.0 / (K_BOLTZMANN * 40) * exp(-(μ * GM * 1e9 / (K_BOLTZMANN * 40)) * (1 / (R * 1e3) - 1 / (r * 1e3)))
        H = K_BOLTZMANN * 40 / (μ * GM * 1e9 / (r * 1e3)^2) / 1e3
        @test p.ω[j] / (-_refractivity_K(atm, 0.7) * n * sqrt(2π * r / H)) ≈ 1 atol = 0.02
    end
    @test 1250 < half_light_radius(p) < 1290
    # tenuous, extended atmosphere (2002 XV93-like): free-molecular-flow tail ω = -πν
    R = 244.0
    GM = 6.6743e-11 * 1500 * 4π / 3 * (R * 1e3)^3 * 1e-9
    q = atmosphere_profile(xv93_atmosphere(:CH4), 124.0, R, GM, 37.0 * AU_KM)
    @test q.y_surface < R
    @test atmosphere_flux(q, q.y_surface - 5) < 0.01       # Fresnel fringes ~ (F/πx)²/2 inside the shadow
    @test 0.3 < _refracted_flux(q, q.y_surface + 1) < 0.7
    @test atmosphere_flux(q, R + 300) > 0.99
    # joint fit of synthetic light curves recovers the surface pressure
    rng = MersenneTwister(2)
    eph = EphemKeplerian(a=37.5, e=0.02, i=12.0, Ω=80.0, ω=150.0, M=10.0, epoch=Epoch("2024-01-01"))
    body = Body("XV93-like", eph; radius=244.0, density=1500.0)
    t0 = Epoch("2024-01-10T15:00:00")
    ra, dec, _ = ephem_radec(body, t0)
    star = Star(ra + 0.001 / 3600 / cosd(dec), dec; diameter=0.0032)
    ev = Event(star, body, t0)
    pth = shadow_path(star, body, ev.tref, ev.tref; offsets=[0.0])
    for (j, dlat) in enumerate([0.2, 1.0, 2.3])
        t = collect(-30:0.5:30.0)
        push!(ev, Chord("A$j", Observer("A$j", pth.lon[1, 1], pth.lat[1, 1] + dlat, 0.0),
                        LightCurve(t, ones(length(t)); tref=ev.tref, σ=0.03, exptime=0.5, λ=0.6)))
    end
    truth = Dict(:center_f => 15.0, :center_g => -20.0, :radius => 244.0, :pressure => 124.0)
    pr = LimbProblem(ev; model=CircleLimb(), photometry=true, atmosphere=xv93_atmosphere(:CH4), init=truth)
    v = [x.value for x in parameters(pr)]
    for d in pr.photometry
        I, _ = Occultations._phot_model(pr, d, Occultations._outlines(pr, v)[1], v)
        d.lc.flux .= I .+ 0.03 .* randn(rng, length(I))
    end
    # two steps: the airless photometric fit gives the geometry, then the atmosphere is added
    r0 = fit(LimbProblem(ev; model=CircleLimb(), photometry=true, init=Dict(:center_f => 0.0, :center_g => 0.0, :radius => 230.0));
             starts=3, rng=rng)
    r = fit(LimbProblem(ev; model=CircleLimb(), photometry=true, atmosphere=xv93_atmosphere(:CH4),
                        init=Dict(:center_f => r0[:center_f][1], :center_g => r0[:center_g][1], :radius => r0[:radius][1],
                                  :pressure => 50.0)); starts=6, rng=rng)
    vt = [truth[n] for n in (:center_f, :center_g, :radius, :pressure)]
    @test r.chi2 <= sum(abs2, residuals(r.problem, vt)) + 1e-6       # at least as good as the truth
    @test 0.8 < reduced_chi2(r) < 1.2
    @test r[:pressure][1] ≈ 124 atol = 4r[:pressure][2]
    @test r[:radius][1] ≈ 244 atol = 4r[:radius][2]
    @test 1 < r[:pressure][2] < 10                   # realistic uncertainty (smooth model)
    @test r0.chi2 - r.chi2 > 100                     # the airless model is rejected
    @test plot(q) isa Plots.Plot
end

@testset "fundamental-plane ephemeris (CORA / Besselian elements)" begin
    body, star, ev, lon0, lat0 = scene()
    p = positionv(star, body, GEOCENTER, T0; dt=2.0)
    p2 = positionv(star, body, GEOCENTER, T0 + 1; dt=2.0)
    ep = EphemPlane(star, T0, [-p.f, -p.vf, -(p2.vf - p.vf) / 2], [-p.g, -p.vg, -(p2.vg - p.vg) / 2]; dist=p.dist)
    obs = Observer("X", lon0, lat0 + 1.0, 500.0)
    for s in (-120.0, 0.0, 60.0)
        q1 = positionv(star, body, obs, T0 + s); q2 = positionv(star, ep, obs, T0 + s)
        @test hypot(q2.f - q1.f, q2.g - q1.g) < 0.02
    end
    row = join([i == 19 ? "22.5" : i == 24 ? "2460000.5" : i == 26 ? "0.1" : i == 27 ? "-0.2" :
                i == 28 ? "3.0" : i == 29 ? "0.5" : "x" for i in 1:30], "|")
    c = parse_cora(row)
    @test c.epoch == Epoch(2460000.5) && (c.x, c.y, c.dx, c.dy, c.vel) == (0.1, -0.2, 3.0, 0.5, 22.5)
end

@testset "real event: (38628) Huya, 2019-03-18 (Santos-Sanz et al. 2022)" begin
    include(joinpath(@__DIR__, "..", "examples", "data", "huya_20190318.jl"))
    ev = huya_event()
    P = HUYA_PAPER
    @test ev.ca.dist_km / AU_KM ≈ P.distance_au atol = 0.01
    @test apparent_diameter(ev.star, ev.ca.dist_km) ≈ P.star_diameter_km atol = 0.02
    @test fresnel_scale(0.6, ev.ca.dist_km) ≈ P.fresnel_km atol = 0.01
    @test count(iscomplete, ev.chords) == 18 && count(isnegative, ev.chords) == 2
    # published ellipse (Table A.4): radial residuals, limb allowance for χ²/dof = 1
    r = fit_limb(ev; residual=:radial, σ_limb=:auto, rng=MersenneTwister(1))
    @test reduced_chi2(r) ≈ 1 atol = 0.02
    @test r[:equatorial_radius][1] ≈ P.a[1] atol = P.a[2]
    @test r[:semi_minor_axis][1] ≈ P.b[1] atol = P.b[2]
    @test r[:position_angle][1] ≈ P.position_angle[1] atol = P.position_angle[2]
    @test 2r[:equivalent_radius][1] ≈ P.D_eq[1] atol = P.D_eq[2]
    @test 0.5 < 2r[:equivalent_radius][2] / P.D_eq[2] < 1.5
    @test all(x -> !x.conflict, check_negative_chords(r))
end

@testset "real events: (119951) 2002 KX14, 2020–2023 (Rizos et al. 2025)" begin
    include(joinpath(@__DIR__, "..", "examples", "data", "2002kx14_2020_2023.jl"))
    evs = kx14_events()
    P = KX14_PAPER
    for (ev, s) in zip(evs, KX14_STARS)
        @test ev.ca.dist_km / AU_KM ≈ P.distance_au[s[1]] atol = 0.01
        @test ev.ca.vel ≈ P.velocity[s[1]] atol = 0.02
    end
    @test sum(c -> count(ispositive, c.chords), evs) == 15
    rng = MersenneTwister(3)
    # event D (two chords) has a mirror solution; the global minimum puts it on the far side
    r1 = fit_limb(evs; rng=rng)
    @test r1[:center_g_4][1] < 0 < r1[:center_g_3][1]
    # a prior that D's ephemeris offset is close to C's (two days apart) selects the other one
    c3 = (r1[:center_f_3][1], r1[:center_g_3][1])
    rJ = fit_limb(evs; priors=Dict(:center_f_4 => (c3[1], 30.0), :center_g_4 => (c3[2], 30.0)),
                  init=Dict(:center_f_4 => c3[1], :center_g_4 => c3[2]), scale_errors=true, rng=rng)
    @test abs(rJ[:center_g_4][1] - rJ[:center_g_3][1]) < 100
    @test rJ[:equatorial_radius][1] ≈ P.a[1] atol = 2P.a[2]
    @test 2rJ[:equivalent_radius][1] ≈ P.D_eq[1] atol = 3P.D_eq[2]
    @test mod(rJ[:position_angle][1] + 90, 180) ≈ P.A.pa[1] atol = 3P.A.pa[2]   # paper: angle of the major axis
    @test all(x -> !x.conflict, check_negative_chords(rJ; event=1))
    @test abs(rJ[:center_f_5][1]) > 1000                       # event E: offset star (Gaia RUWE 1.5)
end

@testset "real event: (143707) 2003 UY117, 2020-10-23 (Kretlow et al. 2024)" begin
    include(joinpath(@__DIR__, "..", "examples", "data", "2003uy117_20201023.jl"))
    ev = uy117_event()
    P = UY117_PAPER
    @test ev.ca.dist_km / AU_KM ≈ P.distance_au atol = 0.01
    @test ev.ca.vel ≈ P.velocity atol = 0.02
    @test count(ispositive, ev.chords) == 4
    @test count(isnegative, ev.chords) == 13
    r = fit_limb(ev; residual=:radial, rng=MersenneTwister(1))
    @test r[:equatorial_radius][1] ≈ P.a[1] atol = P.a[2]
    @test r[:semi_minor_axis][1] ≈ P.b[1] atol = P.b[2]
    @test mod(r[:position_angle][1] - P.position_angle[1] + 90, 180) - 90 ≈ 0 atol = P.position_angle[2]
    @test 2r[:equivalent_radius][1] ≈ P.D_eq[1] atol = P.D_eq[2]
    @test reduced_chi2(r) ≈ P.chi2_pdf atol = 0.1
    @test all(x -> !x.conflict, check_negative_chords(r))
    # 3-D ellipsoid with the light-curve amplitude: a valley of equally good solutions
    amp = RotationalAmplitude(ev.tref, P.amplitude...)
    prob = LimbProblem(ev; model=EllipsoidLimb(principal=true), amplitudes=[amp],
                       init=Dict(:a => 160.0, :b => 110.0, :c => 90.0, :pole_ra => 337.0, :pole_dec => 62.0,
                                 :W0 => 30.0, :period => P.period_h, :center_f => r[:center_f][1], :center_g => r[:center_g][1]),
                       priors=Dict(:pole_dec => isotropic_pole()),
                       bounds=Dict(:a => (20.0, 400.0), :b => (20.0, 400.0), :c => (20.0, 400.0)))
    sols = [fit(prob; rng=MersenneTwister(k), starts=8) for k in 1:3]
    @test all(s -> isapprox(s.chi2, sols[1].chi2; rtol=1e-3), sols)
    for s in sols
        @test s[:c][1] <= s[:b][1] <= s[:a][1]                         # rotation about the shortest axis
        @test abs(residuals(prob, s.values)[end]) < 0.01                 # the amplitude is matched
        @test 2s[:volume_equivalent_radius][1] ≈ P.D_eq_volume[1] atol = 2P.D_eq_volume[2]
    end
end

@testset "plots of light curves and timing fits" begin
    lc = simulate_lightcurve(collect(0.1:0.2:60.0), 20.0, 30.0; tref=T0, vel=20.0, dist=6e9, σ=0.02)
    @test plot(lc) isa Plots.Plot
    @test plot(fit_timing(lc; vel=20.0, dist=6e9)) isa Plots.Plot
    body, star, ev, lon0, lat0 = scene()
    @test plot(shadow_path(star, body, T0 - 30, T0 + 30; step=10)) isa Plots.Plot
end


@testset "binary systems and contact binaries" begin
    # --- geometry of composite outlines
    c = CompositeOutline([EllipseOutline(100.0, 100.0, 0.0), EllipseOutline(50.0, 50.0, 0.0)], [(0.0, 0.0), (300.0, 0.0)])
    iv = intervals(c, -500.0, 0.0, 10.0, 0.0)
    @test length(iv) == 2 && all(isapprox.(iv[1], (40.0, 60.0); atol=1e-9)) && all(isapprox.(iv[2], (75.0, 85.0); atol=1e-9))
    @test isempty(intervals(c, -500.0, 200.0, 10.0, 0.0))
    @test area(c) ≈ π * (100^2 + 50^2)
    @test signed_distance(c, 300.0, 0.0) ≈ -50 && signed_distance(c, 200.0, 0.0) ≈ 50
    lens = CompositeOutline([EllipseOutline(1.0, 1.0, 0.0), EllipseOutline(1.0, 1.0, 0.0)], [(0.0, 0.0), (1.0, 0.0)])
    @test area(lens) ≈ 2π - (2acos(0.5) - sqrt(3) / 2) rtol = 1e-5
    @test length(intervals(lens, -5.0, 0.0, 1.0, 0.0)) == 1                        # overlapping parts: merged
    @test crossing(lens, -5.0, 0.0, 1.0, 0.0, true) == (4.0, 0.0)
    s, miss = crossing(c, -500.0, 200.0, 10.0, 0.0, true)
    @test miss ≈ 100.0
    x, y = outline_points(lens)
    @test count(isnan, x) == 1 && all(i -> isnan(x[i]) || abs(signed_distance(lens, x[i], y[i])) < 1e-6, eachindex(x))
    # two equal spheres: the lens volume π(4R + d)(2R − d)²/12
    R, d = 10.0, 15.0
    b = BilobedEllipsoid(Ellipsoid(R), Ellipsoid(R), d)
    @test volume(b) ≈ 2 * 4π / 3 * R^3 - π * (4R + d) * (2R - d)^2 / 12 rtol = 1e-4        # quadrature of the overlap
    far = BilobedEllipsoid(Ellipsoid(20.0, 15.0, 10.0), Ellipsoid(12.0, 10.0, 8.0), 200.0)
    M = Occultations._body_matrix(30.0, 20.0, 40.0)
    o = project(far, M, [0.6, 0.0, 0.8])
    @test o isa CompositeOutline && area(o) ≈ area(o.parts[1]) + area(o.parts[2])
    @test volume(ellipsoid_mesh(far; subdivisions=4)) ≈ volume(far) rtol = 0.01
    # --- chords with several drops
    obs = Observer("X", 10.0, 45.0, 0.0)
    ch = Chord("X", obs; contacts=[Contact(:immersion, T0, 0.1), Contact(:emersion, T0 + 5, 0.1),
                                   Contact(:immersion, T0 + 8, 0.1; body=:secondary), Contact(:emersion, T0 + 9, 0.1)])
    @test iscomplete(ch) && length(occulted_intervals(ch)) == 2 && ch.contacts[3].body === :secondary
    @test_throws ArgumentError Chord("Y", obs; contacts=[Contact(:immersion, T0, 0.1), Contact(:immersion, T0 + 1, 0.1)])
    @test_throws ArgumentError Contact(:immersion, T0, 0.1; body=:third)
    @test shift(ch, 1.0).contacts[3].body === :secondary
    # --- light curves with several drops: two adjacent strips are one strip (Babinet)
    t = collect(-2:0.02:4.0)
    I2 = occultation_flux(t, [(0.0, 1.0), (1.0, 2.0)]; vel=1.0, dist=6e9, λ=0.6, Δλ=0.0)
    I1 = occultation_flux(t, 0.0, 2.0; vel=1.0, dist=6e9, λ=0.6, Δλ=0.0)
    @test maximum(abs.(I2 .- I1)) < 1e-8
    I3 = occultation_flux(t, [(0.0, 1.0), (60.0, 61.0)]; vel=20.0, dist=6e9)
    @test I3 ≈ occultation_flux(t, 0.0, 1.0; vel=20.0, dist=6e9) .* occultation_flux(t, 60.0, 61.0; vel=20.0, dist=6e9)

    # --- a body with a satellite
    rng = MersenneTwister(21)
    body, star, ev, lon0, lat0 = scene(; radius=300.0)
    p = positionv(star, EPH, GEOCENTER, ev.tref)
    al = (p.vf, p.vg) ./ hypot(p.vf, p.vg); ac = (-al[2], al[1])
    sep = 900.0 .* al .+ 450.0 .* ac
    truth = CompositeOutline([EllipseOutline(320.0, 290.0, 40.0), EllipseOutline(80.0, 80.0, 0.0)], [(0.0, 0.0), sep])
    center = (10.0, -5.0)
    offs = [-360.0, -250, -120, 0, 120, 250, 330, 390, 430, 470, 510, 560] .+ (center[1] * ac[1] + center[2] * ac[2])
    path = shadow_path(star, body, ev.tref, ev.tref; offsets=offs)
    for j in eachindex(offs)
        o = Observer("B$j", path.lon[1, j], path.lat[1, j], 0.0)
        pc = predicted_contacts(ev, o, truth; center=center, window=300, step=0.05)
        cs = reduce(vcat, [[Contact(:immersion, a + 0.05randn(rng), 0.05), Contact(:emersion, b + 0.05randn(rng), 0.05)]
                           for (a, b) in pc]; init=Contact[])
        push!(ev, isempty(cs) ? Chord("B$j", o; span=(ev.tref - 200, ev.tref + 200)) :
                  Chord("B$j", o; contacts=cs, span=(ev.tref - 200, ev.tref + 200)))
    end
    @test count(c -> ispositive(c), ev.chords) >= 9
    r = fit_limb(ev; model=BinaryLimb(), rng=rng)
    @test r.problem.visible === :soft && !isempty(r.problem.visibles)
    @test reduced_chi2(r) < 3
    for (n, tv) in ((:equatorial_radius, 320.0), (:semi_minor_axis, 290.0), (:sec_radius, 80.0), (:sep_f, sep[1]),
                    (:sep_g, sep[2]), (:center_f, center[1]), (:center_g, center[2]))
        @test abs(r[n][1] - tv) < 4r[n][2] + 1e-6
    end
    @test r[:separation][1] ≈ hypot(sep...) atol = 4r[:separation][2]
    @test r[:mass_ratio][1] ≈ r[:size_ratio][1]^3
    a0, ab, ap = astrometry(r), astrometry(r; reference=:barycenter), astrometry(r; reference=:photocenter)
    @test hypot(ab.offset_ra - a0.offset_ra, ab.offset_dec - a0.offset_dec) <
          hypot(ap.offset_ra - a0.offset_ra, ap.offset_dec - a0.offset_dec)              # q³ < area ratio
    @test_throws ArgumentError astrometry(fit_limb(ev; starts=1, rng=rng); reference=:barycenter)   # not a binary model
    @test all(x -> abs(x.dt) < 0.3, contact_residuals(r))
    re = fit_limb(ev; model=BinaryLimb(EllipseLimb(), EllipseLimb()), rng=rng)
    @test :sec_oblateness in re.names
    @test re[:sec_equivalent_radius][1] ≈ 80.0 atol = 3
    @test plot(r) isa Plots.Plot

    # --- a contact binary in one event (2-D)
    rng = MersenneTwister(5)
    body, star, ev, lon0, lat0 = scene(; radius=150.0)
    p = positionv(star, EPH, GEOCENTER, ev.tref)
    ac = (-p.vg, p.vf) ./ hypot(p.vf, p.vg)
    a1, b1, a2, b2, η, pa = 120.0, 90.0, 95.0, 75.0, 0.1, 80.0
    dd = (a1 + a2) * (1 - η)
    x1, x2 = -dd * a2 * b2 / (a1 * b1 + a2 * b2), dd * a1 * b1 / (a1 * b1 + a2 * b2)
    truth = CompositeOutline([EllipseOutline(a1, b1, pa + 90), EllipseOutline(a2, b2, pa + 90)],
                             [(x1 * sind(pa), x1 * cosd(pa)), (x2 * sind(pa), x2 * cosd(pa))])
    offs = collect(-150.0:20.0:150.0)
    path = shadow_path(star, body, ev.tref, ev.tref; offsets=offs)
    for j in eachindex(offs)
        o = Observer("C$j", path.lon[1, j], path.lat[1, j], 0.0)
        pc = predicted_contacts(ev, o, truth; window=300, step=0.02)
        cs = reduce(vcat, [[Contact(:immersion, a + 0.02randn(rng), 0.02), Contact(:emersion, b + 0.02randn(rng), 0.02)]
                           for (a, b) in pc]; init=Contact[])
        push!(ev, isempty(cs) ? Chord("C$j", o; span=(ev.tref - 100, ev.tref + 100)) :
                  Chord("C$j", o; contacts=cs, span=(ev.tref - 100, ev.tref + 100)))
    end
    rc = fit_limb(ev; model=ContactBinaryLimb(), rng=rng)
    @test reduced_chi2(rc) < 3
    for (n, tv) in ((:a1, a1), (:b1, b1), (:a2, a2), (:b2, b2), (:axis_pa, pa), (:neck, η))
        @test abs(rc[n][1] - tv) < 4rc[n][2] + 0.5
    end
    @test rc[:equivalent_radius][1] ≈ equivalent_radius(truth) rtol = 0.01
    @test rc[:neck_width][1] > 0
    @test reduced_chi2(fit_limb(ev; rng=rng)) > 10reduced_chi2(rc)                   # a single ellipse cannot fit
    @test plot(rc) isa Plots.Plot

    # --- a contact binary in 3-D: outlines, derived quantities, amplitude
    body, star, ev, lon0, lat0 = scene(; radius=150.0, pole=(35.0, 40.0), rotation_period=9.5, W0=20.0,
                                       shape=BilobedEllipsoid(Ellipsoid(130.0, 100.0, 85.0), Ellipsoid(100.0, 85.0, 75.0), 207.0))
    O = project(body, ev.tref)
    add_chords!(ev, O, (0.0, 0.0), lon0, lat0, [-1.0, -0.5, 0.0, 0.5, 1.0]; rng=MersenneTwister(1))
    amp = RotationalAmplitude(ev.tref, rotational_amplitude(body.shape, 35.0, 40.0, normalize(geocentric_vector(EPH, ev.tref))), 0.02)
    pr = LimbProblem(ev; model=ContactBinaryEllipsoids(), amplitudes=[amp], fix=[:pole_ra, :pole_dec])
    v = [q.value for q in parameters(pr)]
    @test outlines(pr, v)[1][2] isa CompositeOutline
    @test area(outlines(pr, v)[1][2]) ≈ area(O) rtol = 1e-6                          # the start is the truth here
    @test abs(residuals(pr, v)[end]) < 1e-6                                          # the amplitude is matched
    dv = derived(pr, v)
    @test dv.volume_equivalent_radius ≈ volume_equivalent_radius(body.shape)
    @test dv.mass_ratio < 1
end


@testset "ring systems" begin
    # --- geometry
    r = ProjectedRing(400.0, 10.0, 0.3, 0.0, 90.0)                      # seen pole-on: circles
    iv = intervals(r, -1000.0, 0.0, 10.0, 0.0)                            # a track through the centre
    @test length(iv) == 2
    @test all(isapprox.(iv[1], (59.5, 60.5); atol=1e-9)) && all(isapprox.(iv[2], (139.5, 140.5); atol=1e-9))
    @test length(intervals(r, -1000.0, 403.0, 10.0, 0.0)) == 1               # through the annulus only
    @test isempty(intervals(r, -1000.0, 410.0, 10.0, 0.0))
    u = [0.3, -0.5, 0.8]; u = u ./ norm(u)
    B, P = ring_geometry(151.3, 41.5, u)
    α, δ = Occultations._ring_pole(B, P, u)
    @test α ≈ 151.3 atol = 1e-8
    @test δ ≈ 41.5 atol = 1e-8
    pr = project_ring(391.0, 6.2, 0.4, 151.3, 41.5, u)
    @test pr.opacity ≈ 1 - exp(-0.4 / abs(sind(B))) && pr.opening ≈ B && pr.pa ≈ P
    am, bm = Occultations._ring_pole(-B, P, u)                         # the mirror plane: the same ellipse
    @test project_ring(391.0, 6.2, 0.4, am, bm, u).pa ≈ P atol = 1e-8
    @test Contact(:immersion, T0, 0.1; body=:ring2).body === :ring2
    @test_throws ArgumentError Contact(:immersion, T0, 0.1; body=:ringX)
    # --- light curves: a wide, partly transparent strip lets 1 − p′ through
    t = collect(-5:0.05:15.0)
    I = occultation_flux(t, [(0.0, 10.0)]; vel=10.0, dist=6e9, opacity=[0.4])
    @test median(I[(t .> 3) .& (t .< 7)]) ≈ 0.6 atol = 0.02
    I2 = occultation_flux(t, [(0.0, 1.0), (1.0, 2.0)]; vel=1.0, dist=6e9, λ=0.6, Δλ=0.0, opacity=[0.3, 0.3])
    @test maximum(abs.(I2 .- occultation_flux(t, 0.0, 2.0; vel=1.0, dist=6e9, λ=0.6, Δλ=0.0, opacity=0.3))) < 1e-8

    # --- a body with a ring, from timings
    rng = MersenneTwister(31)
    body, star, ev, lon0, lat0 = scene(; radius=200.0)
    U = EventContext(ev).u
    ring = project_ring(900.0, 8.0, 0.3, 120.0, 35.0, U)
    pv = positionv(star, EPH, GEOCENTER, ev.tref)
    offs = [-980.0, -930, -850, -600, -300, -150, 0, 120, 250, 500, 820, 900, 960]
    path = shadow_path(star, body, ev.tref, ev.tref; offsets=offs)
    for j in eachindex(offs)
        o = Observer("R$j", path.lon[1, j], path.lat[1, j], 0.0)
        bod = predicted_contacts(ev, o, EllipseOutline(200.0, 200.0, 0.0); window=400, step=0.05)
        cs = reduce(vcat, [[Contact(:immersion, a + 0.02randn(rng), 0.02), Contact(:emersion, b + 0.02randn(rng), 0.02)]
                           for (a, b) in bod]; init=Contact[])
        for (a, b) in predicted_contacts(ev, o, ring; window=600)
            any(x -> a < x[2] && b > x[1], bod) && continue
            ti, te = minmax(a + 0.02randn(rng), b + 0.02randn(rng))
            push!(cs, Contact(:immersion, ti, 0.02; body=:ring1), Contact(:emersion, te, 0.02; body=:ring1))
        end
        sort!(cs; by=x -> x.time)
        push!(ev, isempty(cs) ? Chord("R$j", o; span=(ev.tref - 300, ev.tref + 300)) :
                  Chord("R$j", o; contacts=cs, span=(ev.tref - 300, ev.tref + 300)))
    end
    @test count(c -> any(x -> x.body === :ring1, c.contacts), ev.chords) >= 6
    rr = fit_limb(ev; model=RingedLimb(CircleLimb(); rings=1), fix=[:ring1_tau], rng=rng)
    @test reduced_chi2(rr) < 3
    @test rr[:radius][1] ≈ 200 atol = 4rr[:radius][2] + 0.1
    @test rr[:ring1_radius][1] ≈ 900 atol = 4rr[:ring1_radius][2] + 0.1
    @test rr[:ring1_width][1] ≈ 8 atol = 4rr[:ring1_width][2] + 0.1
    @test rr[:ring_pa][1] ≈ ring.pa atol = 0.5
    @test abs(rr[:ring_opening][1]) ≈ abs(ring.opening) atol = 0.5          # the plane or its mirror
    @test rr[:ring_pole_dec][1] >= 0                                       # canonical pole
    @test all(x -> abs(x.dt) < 0.2, contact_residuals(rr))
    @test haskey(rr.derived, :ring1_equivalent_depth) && haskey(rr.derived, :ring1_opacity)
    @test plot(rr) isa Plots.Plot
end


@testset "real events: 2002 TC302, 2003 AZ84, 2003 VS2, 2002 GZ32, Makemake, Haumea" begin
    D = joinpath(@__DIR__, "..", "examples", "data")
    # 2002 TC302 (Ortiz et al. 2020): the ellipse with the limb allowance
    include(joinpath(D, "2002tc302_20180128.jl"))
    ev = tc302_event()
    @test ev.ca.vel ≈ TC302_PAPER.velocity atol = 0.02
    r = fit_limb(ev; σ_limb=:auto, rng=MersenneTwister(1))
    @test r[:equatorial_radius][1] ≈ TC302_PAPER.a[1] atol = r[:equatorial_radius][2]
    @test r[:semi_minor_axis][1] ≈ TC302_PAPER.b[1] atol = r[:semi_minor_axis][2]
    @test 2r[:equivalent_radius][1] ≈ TC302_PAPER.D_eq atol = 2r[:equivalent_radius][2]
    @test all(x -> !x.conflict, check_negative_chords(r))
    # 2003 AZ84 (Dias-Oliveira et al. 2017), 2014: the same star position as the paper, the same body position
    include(joinpath(D, "2003az84_2012_2014.jl"))
    e14 = az84_event(2014)
    @test e14.ca.vel ≈ 9.4 atol = 0.05
    r14 = fit_limb(e14; rng=MersenneTwister(1))
    @test r14[:equatorial_radius][1] ≈ AZ84_PAPER.e2014.a[1] atol = 3
    @test r14[:equivalent_radius][1] ≈ AZ84_PAPER.e2014.R_equiv[1] atol = 3
    t = Epoch(AZ84_PAPER.e2014.position.epoch)
    ra_, dec_, _ = ephem_radec(e14.body, t)
    a = astrometry(r14)
    r_, d_ = Occultations.offset_radec(ra_, dec_, a.offset_ra, a.offset_dec)
    @test abs(r_ - parse_ra(AZ84_PAPER.e2014.position.ra)) * cosd(d_) * 3.6e6 < 5
    @test abs(d_ - parse_dec(AZ84_PAPER.e2014.position.dec)) * 3.6e6 < 5
    # Bienor (Rizos et al. 2024), Occ. B: the ellipse; a 3-D fit over three events with phases
    include(joinpath(D, "bienor_2019_2023.jl"))
    rB = fit_limb(bienor_event("B"); rng=MersenneTwister(1))
    @test rB[:equatorial_radius][1] ≈ BIENOR_PAPER.B.u[1] atol = 3BIENOR_PAPER.B.u[2]
    @test rB[:semi_minor_axis][1] ≈ BIENOR_PAPER.B.v[1] atol = 2BIENOR_PAPER.B.v[2]
    ebs = [bienor_event(n) for n in ("A", "B", "C")]
    pb = LimbProblem(ebs; model=EllipsoidLimb(principal=true),
                     phases=[RotationalPhase(e.tref, BIENOR_PAPER.phases[n]...) for (e, n) in zip(ebs, ("A", "B", "C"))],
                     fix=[:pole_ra, :pole_dec], init=Dict(:a => 120.0, :b => 60.0, :c => 45.0, :pole_ra => 3.58,
                                                          :pole_dec => 58.23, :period => 9.1736))
    @test length(pb.phases) == 3                                  # each phase with its own event's ephemeris
    # Chariklo (Morgado et al. 2021), 2017-07-23: rings from the mid-times (CDS table)
    include(joinpath(D, "chariklo_2017_2020.jl"))
    ec = chariklo_event("2017-07-23")
    rc = fit_limb(ec; model=RingedLimb(EllipseLimb(); rings=2),
                  fix=[:ring1_width, :ring1_tau, :ring2_width, :ring2_tau],
                  init=Dict(:ring_pole_ra => 151.0, :ring_pole_dec => 41.8, :ring1_radius => 386.0, :ring1_width => 7.0,
                            :ring1_tau => 0.4, :ring2_radius => 400.0, :ring2_width => 3.0, :ring2_tau => 0.05),
                  starts=2, rng=MersenneTwister(1))
    @test rc[:ring1_radius][1] ≈ 385.9 atol = 2
    @test rc[:ring_pole_ra][1] ≈ 151.03 atol = 0.5
    @test rc[:ring_pole_dec][1] ≈ 41.81 atol = 0.5
    # σ_ring: ring contacts get their own allowance, body contacts keep σ_limb
    pr = LimbProblem(ec; model=RingedLimb(EllipseLimb(); rings=2), σ_limb=2.0, σ_ring=20.0)
    @test pr.σ_limb == 2.0 && pr.σ_ring == 20.0
    @test LimbProblem(ec; model=RingedLimb(EllipseLimb(); rings=2), σ_limb=3.0).σ_ring == 3.0
    # Quaoar 2022 (Pereira et al. 2023): the limb with an allowance
    include(joinpath(D, "quaoar_2011_2022.jl"))
    rq = fit_limb(quaoar_event(2022; rings=false); σ_limb=6.6, rng=MersenneTwister(1))
    @test rq[:equivalent_radius][1] ≈ QUAOAR_PAPER.e2022.R_equiv[1] atol = 3
    # Chiron 2019 (Braga-Ribas et al. 2023): the Paris clock offset
    include(joinpath(D, "chiron_2018_2019.jl"))
    rch = fit_limb(chiron_event(2019); fit_dt=["Paris"], rng=MersenneTwister(1))
    @test clock_offsets(rch)["Paris"][1] ≈ CHIRON_PAPER.paris_shift atol = 0.6
    # 2002 MS4 (Rommel et al. 2023): the global ellipse of the 13 selected chords
    include(joinpath(D, "2002ms4_20200808.jl"))
    rm4 = fit_limb(ms4_event(selected=true); σ_limb=7.0, rng=MersenneTwister(1))
    @test rm4[:equivalent_radius][1] ≈ MS4_PAPER.ellipse.R_equiv[1] atol = MS4_PAPER.ellipse.R_equiv[2]
    # the single-chord events with the Gaia DR3 stars: the paper's chord lengths
    for (y, (L, σ)) in zip((2011, 2013), AZ84_PAPER.single_chords_km)
        e1 = az84_event(y)
        f1, g1 = chord_path(e1, e1.chords[1]; segment=:occulted)
        @test hypot(f1[end] - f1[1], g1[end] - g1[1]) ≈ L atol = σ
    end
    # 2003 VS2 (Vara-Lubiano et al. 2022): the clock offset of chord 8 (the paper: −14 s)
    include(joinpath(D, "2003vs2_20191022.jl"))
    ev = vs2_event()
    @test ev.ca.vel ≈ VS2_PAPER.velocity atol = 0.02
    r = fit_limb(ev; fit_dt=["Romania 8"], rng=MersenneTwister(1))
    @test r[Symbol("dt_Romania 8")][1] ≈ -14.0 atol = 1.0
    @test r[:equatorial_radius][1] ≈ VS2_PAPER.original.a[1] atol = 2VS2_PAPER.original.a[2]
    @test 2r[:equivalent_radius][1] ≈ VS2_PAPER.original.D_eq[1] atol = VS2_PAPER.original.D_eq[2]
    # 2002 GZ32 (Santos-Sanz et al. 2021): an elongated limb
    include(joinpath(D, "2002gz32_20170520.jl"))
    ev = gz32_event()
    r = fit_limb(ev; rng=MersenneTwister(1))
    @test r[:semi_minor_axis][1] ≈ GZ32_PAPER.shifted.b[1] atol = 2GZ32_PAPER.shifted.b[2]
    @test r[:equatorial_radius][1] / r[:semi_minor_axis][1] > 1.6
    # Makemake (Ortiz et al. 2012): the circle
    include(joinpath(D, "makemake_20110423.jl"))
    r = fit_limb(makemake_event(); model=CircleLimb(), rng=MersenneTwister(1))
    @test 2r[:radius][1] ≈ MAKEMAKE_PAPER.circle_D[1] atol = 2MAKEMAKE_PAPER.circle_D[2]
    # Haumea (Ortiz et al. 2017): the limb, and the ring from the mid-times of its events
    include(joinpath(D, "haumea_20170121.jl"))
    rb = fit_limb(haumea_event(rings=false); rng=MersenneTwister(1))
    @test rb[:equatorial_radius][1] ≈ HAUMEA_PAPER.a[1] atol = 3HAUMEA_PAPER.a[2]
    @test rb[:semi_minor_axis][1] ≈ HAUMEA_PAPER.b[1] atol = 2HAUMEA_PAPER.b[2]
    ev = haumea_event()
    @test count(c -> any(x -> x.kind === :mid, c.contacts), ev.chords) == 8      # incl. Mount Agliale (ring only)
    rr = fit_limb(ev; model=RingedLimb(EllipseLimb(); rings=1), fix=[:ring1_width, :ring1_tau],
                  init=Dict(:ring1_width => 70.0, :ring1_tau => 0.165, :ring_pole_ra => 312.3, :ring_pole_dec => -18.6,
                            :ring1_radius => 2287.0), starts=2, rng=MersenneTwister(2))
    @test abs(rr[:ring_opening][1]) ≈ HAUMEA_PAPER.ring_B[1] atol = 1.5
    @test rr[:ring1_radius][1] ≈ HAUMEA_PAPER.ring_a[1] atol = 100
    @test mod(rr[:ring_pa][1] - (HAUMEA_PAPER.ring_pa[1] + 360) + 90, 180) - 90 ≈ 0 atol = 3
end

@testset "CORA observation reports" begin
    dir = mktempdir()
    report(result, extra) = """
        Date         = 2023-02-02                {yyyy-mm-dd}
        Object       = (447178) 2005 RO43        {object}
        Star         = GDR2 2951529567657635072  {star}
        Occultation  = $result                   {result}
        Sitename     = Site $result              {site}
        Latitude     = S 31 16 27.5              {latitude}
        Longitude    = E 149:03:42.0             {longitude}
        Altitude     = 1150                      {m}
        UT1          = 23:50:00                  {start}
        UT4          = 00:10:00                  {end}
        $extra
        Comment      =                           {comment}
        """
    fp = joinpath(dir, "pos.txt")
    write(fp, report("O+", "UT2 = 23:59:58.50 {D}\nUT3 = 00:00:21.75 {R}\nUT2ERR = 0.05 {s}\nUT3ERR = 0.08 {s}"))
    r = read_cora_report(fp)
    @test r.result === :positive
    @test r.site.lat ≈ -(31 + 16 / 60 + 27.5 / 3600)
    @test r.site.lon ≈ 149 + 3 / 60 + 42 / 3600
    @test r.site.height == 1150.0
    @test r.emersion - r.immersion ≈ 23.25 atol = 1e-6            # across midnight
    @test r.stop - r.start ≈ 1200.0 atol = 1e-6
    c = cora_chord(r; sensitivity=0.5)
    @test ispositive(c) && length(c.contacts) == 2 && c.contacts[2].σ == 0.08 && c.sensitivity == 0.5
    fn = joinpath(dir, "neg.txt")
    write(fn, report("O-", "UT2 = 00:00:01 {D}"))
    rn = @test_logs (:warn,) read_cora_report(fn)
    @test rn.result === :negative && rn.immersion === nothing
    @test isnegative(cora_chord(rn))
    fd = joinpath(dir, "nd.txt")
    write(fd, report("ND", ""))
    @test cora_chord(read_cora_report(fd)) === nothing
end

@testset "double star (photometry)" begin
    eph = EphemKeplerian(a=40.0, e=0.05, i=10.0, Ω=80.0, ω=150.0, M=10.0, epoch=Epoch("2022-01-01"), name="SimDS")
    body = Body("SimDS", eph; radius=150.0)
    t0 = Epoch("2022-06-20T02:00:00")
    ra, dec, _ = ephem_radec(body, t0)
    star = Star(ra + 0.001 / 3600 / cosd(dec), dec; name="double")
    ev = Event(star, body, t0; window=900)
    path = shadow_path(star, body, ev.tref, ev.tref; offsets=[0.0])
    lon0, lat0 = path.lon[1, 1], path.lat[1, 1]
    tt = collect(-40:0.5:40.0)
    for (j, dlat) in enumerate(-6.0:1.0:6.0)                  # stations across the shadows of both stars
        push!(ev, Chord("S$j", Observer("S$j", lon0, lat0 + dlat, 0.0),
                        LightCurve(tt, ones(length(tt)); tref=ev.tref, σ=0.03, exptime=0.5, name="S$j")))
    end
    truth = Dict(:center_f => 20.0, :center_g => -10.0, :equatorial_radius => 160.0, :oblateness => 0.1,
                 :position_angle => 40.0, :star_f1 => 0.75, :star_x1 => -120.0, :star_y1 => -150.0)
    prob = LimbProblem(ev; photometry=true, double_star=true, init=truth)
    v = [q.value for q in parameters(prob)]
    rng = MersenneTwister(5)
    for d in prob.photometry                                  # simulated light curves from the model
        I, _ = Occultations._phot_model(prob, d, Occultations._outlines(prob, v)[1], v)
        d.lc.flux .= I .+ 0.03 .* randn(rng, length(I))
    end
    @test count(d -> minimum(d.lc.flux) < 0.5, prob.photometry) >= 1     # the primary is occulted somewhere
    @test count(d -> 0.6 < minimum(d.lc.flux) < 0.9 && maximum(d.lc.flux) > 0.9, prob.photometry) >= 1   # the secondary
    ## a start near the solution, as from the timings (joint light-curve fits have local minima)
    start = merge(truth, Dict(:star_f1 => 0.72, :star_x1 => -110.0, :star_y1 => -140.0, :equatorial_radius => 155.0))
    r = fit(LimbProblem(ev; photometry=true, double_star=true, init=start); starts=4, rng=rng)
    @test reduced_chi2(r) < 1.3
    @test r[:star_f1][1] ≈ 0.75 atol = 4r[:star_f1][2] + 0.01
    @test r[:star_x1][1] ≈ -120 atol = 4r[:star_x1][2] + 5
    @test r[:star_y1][1] ≈ -150 atol = 4r[:star_y1][2] + 5
    @test r[:equatorial_radius][1] ≈ 160 atol = 4r[:equatorial_radius][2] + 5
    @test r[:star_dmag][1] ≈ 2.5log10(3) atol = 0.1
    @test_throws ArgumentError LimbProblem(ev; double_star=true, photometry=false)
end

end
