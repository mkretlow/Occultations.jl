# Published reductions re-done with SORA: export of the events and the Occultations.jl fits.
#
# For every event: <dir>/ephem.txt (JD UTC, RA, Dec [deg], distance [au]; for SORA's
# EphemPlanete), <dir>/event.txt (star, the used contacts of the positive chords) and
# <dir>/config.txt (the fit set-up and the Occultations.jl result). The fits here use the
# radial residuals (SORA's definition), the timing errors and, where given, a limb allowance
# σ_limb (SORA's `ellipse_error`); negative chords are not used (SORA's fit_ellipse ignores
# them).
#
#   julia --project=. examples/sora/export.jl

using Occultations, Random, Printf, Logging
Logging.disable_logging(Logging.Warn)
const ROOT = joinpath(@__DIR__, "..", "..")
const DATA = joinpath(ROOT, "examples", "data")
const OUT = joinpath(ROOT, "examples", "output", "sora", "cases")
for f in ("varda_20180910", "1995sm55_20240225", "huya_20230624", "2002tc302_20180128", "2007uk126_20141115",
          "2002gz32_20170520", "makemake_20110423", "2002ms4_20200808")
    include(joinpath(DATA, f * ".jl"))
end

# name, event, Horizons table, model (:ellipse or :circle), σ_limb [km]
const CASES = [
    ("varda_2018", () -> varda_event(negatives=:none), "varda_horizons_20180910.txt", :ellipse, 0.0),
    ("1995sm55_2024", () -> sm55_event(negatives=:none), "1995sm55_horizons_20240225.txt", :ellipse, 0.0),
    ("1995sm55_2024_slimb", () -> sm55_event(negatives=:none), "1995sm55_horizons_20240225.txt", :ellipse, :auto),
    ("huya_2023", () -> huya2023_event(negatives=:none), "huya_horizons_20230624.txt", :ellipse, 0.0),
    ("2002tc302_2018", () -> tc302_event(negatives=:none), "2002tc302_horizons_20180128.txt", :ellipse, :auto),
    ("2007uk126_2014", () -> uk126_event(), "2007uk126_horizons_20141115.txt", :ellipse, 0.0),
    ("2002gz32_2017", () -> gz32_event(negatives=:none), "2002gz32_horizons_20170520.txt", :ellipse, 0.0),
    ("makemake_2011", () -> makemake_event(), "makemake_horizons_20110423.txt", :circle, 0.0),
    ("2002ms4_2020", () -> ms4_event(negatives=:none, selected=true), "2002ms4_horizons_20200808.txt", :ellipse, 7.0),
]

function write_ephem(src, dst)
    open(dst, "w") do io
        for l in eachline(joinpath(DATA, src))
            startswith(l, "#") && continue
            t, ra, dec, d = split(l)
            @printf(io, "%.10f %.10f %.10f %.12f\n", sum(utc_jd(Epoch(t))), parse(Float64, ra), parse(Float64, dec),
                    parse(Float64, d) / 1.495978707e8)
        end
    end
end

function write_event(ev, dst)
    s = ev.star
    open(dst, "w") do io
        println(io, "# star: ra dec [deg] pmra pmdec [mas/yr] parallax [mas] epoch [Julian year]")
        @printf(io, "star %.12f %.12f %.6f %.6f %.6f %.8f\n", s.ra, s.dec, s.pmra, s.pmdec, s.parallax, s.epoch)
        println(io, "tref ", isot(ev.tref; digits=3))
        println(io, "# chord name | lon lat [deg] height [m] | kind time σ [s] ...")
        for c in ev.chords
            isnegative(c) && continue
            cs = [k for k in c.contacts if k.use && k.kind in (:immersion, :emersion)]
            isempty(cs) && continue
            o = c.observer
            @printf(io, "chord|%s|%.8f %.8f %.2f", c.name, o.lon, o.lat, o.height)
            for k in cs
                @printf(io, "|%s %s %.4f", k.kind, isot(k.time; digits=4), k.σ)
            end
            println(io)
        end
    end
end

rng = MersenneTwister(1)
for (name, mk, eph, model, σl) in CASES
    dir = mkpath(joinpath(OUT, name))
    ev = mk()
    write_ephem(eph, joinpath(dir, "ephem.txt"))
    write_event(ev, joinpath(dir, "event.txt"))
    m = model === :circle ? CircleLimb() : EllipseLimb()
    t = @elapsed r = fit_limb(ev; model=m, residual=:radial, σ_limb=σl, negative=:none, rng=rng)
    σlimb = r.problem.σ_limb
    function v(q)
        model === :circle && q === :equatorial_radius && return r[:radius]
        model === :circle && q in (:oblateness, :position_angle) && return (0.0, 0.0)
        return r[q]
    end
    open(joinpath(dir, "config.txt"), "w") do io
        println(io, "model ", model)
        @printf(io, "ellipse_error %.4f\n", σlimb)
        for q in (:center_f, :center_g, :equatorial_radius, :oblateness, :position_angle)
            @printf(io, "oj_%s %.6f %.6f\n", q, v(q)...)
        end
        @printf(io, "oj_chi2 %.6f\noj_dof %d\noj_time %.6f\n", r.chi2, r.dof, t)
    end
    @printf("%-22s σ_limb %5.1f  a %8.2f ± %5.2f  ε %6.3f ± %5.3f  PA %7.2f ± %5.2f  χ² %8.2f\n", name, σlimb,
            v(:equatorial_radius)..., v(:oblateness)..., v(:position_angle)..., r.chi2)
end
