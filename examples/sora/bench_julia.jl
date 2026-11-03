# Speed of Occultations.jl on the Varda occultation of 2018 September 10 (Souami et al. 2020),
# for the comparison with SORA (bench_sora.py, same inputs). Writes the synthetic light curve
# used by both and prints the timings.
#
#   julia --project=. examples/sora/bench_julia.jl

using Occultations, Random, Printf, Statistics, LinearAlgebra, Logging
include(joinpath(@__DIR__, "..", "data", "varda_20180910.jl"))

const DIR = mkpath(joinpath(@__DIR__, "..", "output", "sora", "benchmark"))
# the event window (±1200 s) reaches past the stored ephemeris table; the extrapolation warnings
# would dominate the timings
Logging.disable_logging(Logging.Warn)

# first call (includes compilation) and the median of n further calls [s]
function bench(f, n)
    t1 = @elapsed r = f()
    ts = [@elapsed(f()) for _ in 1:n]
    return r, t1, median(ts), n
end
report(name, b) = @printf("%-44s first %9.4f s   median %9.5f s  (n = %d)\n", name, b[2], b[3], b[4])

# ---- common input: the ephemeris table for SORA (JD UTC, RA, Dec [deg], distance [au])
open(joinpath(DIR, "varda_ephem_sora.txt"), "w") do io
    for l in eachline(joinpath(VARDA_DIR, "varda_horizons_20180910.txt"))
        startswith(l, "#") && continue
        t, ra, dec, d = split(l)
        @printf(io, "%.10f %.10f %.10f %.12f\n", sum(utc_jd(Epoch(t))), parse(Float64, ra), parse(Float64, dec),
                parse(Float64, d) / 1.495978707e8)
    end
end

println("Occultations.jl, Julia ", VERSION, ", ", Sys.cpu_info()[1].model, ", threads ", Threads.nthreads())

# ---- 1. the event: star, body, stations, chords (closest approach included)
b1 = bench(() -> varda_event(negatives=:none), 20)
report("1. event set-up (5 positive chords)", b1)
ev = b1[1]

# ---- 2. (f, g) of the ten contact times
function contacts_fg(ev)
    out = NTuple{2,Float64}[]
    for c in ev.chords, k in c.contacts
        push!(out, track(ev, c)(k.time - ev.tref))
    end
    return out
end
b2 = bench(() -> contacts_fg(ev), 20)
report("2. (f, g) of the 10 contacts", b2)

# ---- 3. ellipse fit, radial residuals (as SORA), covariance errors
rng = MersenneTwister(1)
b3 = bench(() -> fit_limb(ev; residual=:radial, rng=rng), 20)
report("3. ellipse fit (radial, LM + covariance)", b3)
r = b3[1]
@printf("   a = %.1f ± %.1f km, ε = %.3f ± %.3f, PA = %.1f ± %.1f°, χ² = %.2f\n", r[:equatorial_radius]...,
        r[:oblateness]..., r[:position_angle]..., r.chi2)

# ---- 4. posterior sampling of the same fit (the counterpart of SORA's Monte-Carlo region)
b4 = bench(() -> sample(r; nsteps=2000, rng=rng), 5)
report("4. posterior (20 walkers × 2000 steps)", b4)
p = b4[1]
show(stdout, MIME"text/plain"(), p); println()

# ---- 4b. χ² profiles (Δχ² = 1 intervals) of the five parameters, SORA's kind of interval
pars = (:center_f, :center_g, :equatorial_radius, :oblateness, :position_angle)
b4b = bench(() -> [confidence_interval(r, q) for q in pars], 5)
report("4b. χ² profiles of the 5 parameters", b4b)
println("   ", join([@sprintf("%s %.3f–%.3f", q, ci...) for (q, ci) in zip(pars, b4b[1])], ", "))

# ---- 5. with the six negative chords as hard constraints
evn = varda_event()
b5 = bench(() -> fit_limb(evn; residual=:radial, rng=rng), 20)
report("5. ellipse fit with 6 negative chords", b5)

# ---- 6. light-curve timing: a synthetic FLO-like light curve (written for SORA)
D = VARDA_PAPER.distance_au * 1.495978707e8
vel, sd, texp = VARDA_PAPER.velocity, VARDA_PAPER.star_diameter_km, 0.5
t = collect(-60:texp:60.0)
lc = simulate_lightcurve(t, -31.7, 31.7; tref=Epoch("2018-09-10T03:34:55"), vel=vel, dist=D, star_diam=sd,
                         σ=0.05, exptime=texp, λ=0.7, Δλ=0.3, rng=MersenneTwister(7), name="FLO-like")
open(joinpath(DIR, "lightcurve.txt"), "w") do io
    println(io, "# t [s from 2018-09-10T03:34:55 UTC], flux; exposure 0.5 s, σ 0.05, true contacts −31.7 and +31.7 s")
    for (x, y) in zip(lc.time, lc.flux)
        @printf(io, "%.3f %.6f\n", x, y)
    end
end
b6 = bench(() -> fit_timing(lc; vel=vel, dist=D, star_diam=sd, rng=rng), 20)
report("6. light-curve timing (diffraction model)", b6)
c = contact_times(b6[1])
@printf("   immersion %+.3f ± %.3f s, emersion %+.3f ± %.3f s (truth −31.700, +31.700)\n",
        c.immersion - lc.tref, c.immersion_err, c.emersion - lc.tref, c.emersion_err)
