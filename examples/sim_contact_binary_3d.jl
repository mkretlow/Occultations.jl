# # Synthetic: a contact binary in 3-D
#
# Four simulated occultations by a rotating bilobed body, reduced with
# [`ContactBinaryEllipsoids`](@ref): two triaxial ellipsoids on the body x-axis, rotating about
# z. Three events in one night see the body at different rotational phases, a fourth one a
# month later; a rotational light-curve amplitude is added. Each event sees the exact
# silhouette of the two ellipsoids (the union of their projected ellipses), so the events
# together constrain the 3-D shape of both lobes:
#
# 1. the simulation: the shape, the pole, the rotation, the chords of each event (eleven
#    stations per event on land and in the dark, with timing accuracies of 0.02–0.4 s);
# 2. the fit with the pole known, from a deliberately poor starting shape;
# 3. the fit with the pole unknown, from a grid of starting poles;
# 4. the shape as an OBJ file and the 3-D view of each event.
#
# Run it with:
#
#     julia --project=. examples/sim_contact_binary_3d.jl [output_dir]      (default: examples/output/sim_contact_binary_3d)

using Occultations, Random, Printf, LinearAlgebra
using Plots
import CairoMakie                     # saves the 3-D views (Makie figures)
ENV["GKSwstype"] = "100"   #src
include(joinpath(@__DIR__, "sim_stations.jl"))       # stations on land, in the dark, with their equipment

out = length(ARGS) >= 1 ? ARGS[1] : mkpath(joinpath(@__DIR__, "output", "sim_contact_binary_3d"))
rng = MersenneTwister(5)
section(s) = println("\n", "="^78, "\n", s, "\n", "="^78)   #src

section("1. the simulation")
# Lobe 1 with semi-axes 130 × 100 × 85 km, lobe 2 with 100 × 85 × 75 km, their centres
# (a₁ + a₂)(1 − 0.1) = 207 km apart; the pole at RA 35°, Dec 40°, a period of 9.5 h.
truth = BilobedEllipsoid(Ellipsoid(130.0, 100.0, 85.0), Ellipsoid(100.0, 85.0, 75.0), (130.0 + 100.0) * (1 - 0.1))
println(truth)
@printf("volume-equivalent diameter %.1f km, mass ratio %.3f (equal densities)\n", 2volume_equivalent_radius(truth),
        volume(truth.lobe2) / volume(truth.lobe1))
eph = EphemKeplerian(a=44.0, e=0.03, i=3.0, Ω=120.0, ω=40.0, M=300.0, epoch=Epoch("2024-01-01"), name="SimCB3")
body = Body("SimCB3", eph; radius=140.0, pole=(35.0, 40.0), rotation_period=9.5, W0=20.0, shape=truth)
t0 = Epoch("2024-12-30T19:00:00")
evs = Event[]
for (k, dt) in enumerate([0.0, 2.1 * 3600, 3.6 * 3600, 30 * 86400 + 1.3 * 3600])
    tk = t0 + dt
    rak, deck, _ = ephem_radec(body, tk)
    st = Star(rak + 0.001 / 3600 / cosd(deck), deck; name="star $k")
    e = Event(st, body, tk; name="event $k")
    Ok = project(body, e.tref)                           # the true silhouette at this event
    pv = positionv(st, eph, GEOCENTER, e.tref)
    n̂ = (-pv.vg, pv.vf) ./ hypot(pv.vf, pv.vg)
    c = (6.0k, -4.0k)                                    # ephemeris offset of the event [km]
    offs = random_offsets(MersenneTwister(20 + k), 11, -200.0, 200.0; mingap=10.0) .+ (c[1] * n̂[1] + c[2] * n̂[2])
    for site in stations(st, body, e.tref, offs; rng=rng, prefix="E$k-S")
        o, σ = site.observer, site.equipment[2]
        pc = predicted_contacts(e, o, Ok; center=c, window=abs(site.τ) + 200, step=0.05)
        span = (e.tref + site.τ - 100, e.tref + site.τ + 100)
        cs = reduce(vcat, [[Contact(:immersion, a + σ * randn(rng), σ), Contact(:emersion, b + σ * randn(rng), σ)]
                           for (a, b) in pc]; init=Contact[])
        push!(e, isempty(cs) ? Chord(o.name, o; span=span) : Chord(o.name, o; contacts=cs, span=span))
    end
    @printf("%s: %d positive chords (%d with two drops), area-equivalent diameter %.1f km\n", e.name,
            count(ispositive, e.chords), count(c -> count(x -> x.kind === :immersion, c.contacts) == 2, e.chords),
            2equivalent_radius(Ok))
    push!(evs, e)
end
# A rotational light curve observed ten days after the first event, amplitude ± 0.02 mag:
u = normalize(geocentric_vector(eph, t0 + 10 * 86400))
amp = RotationalAmplitude(t0 + 10 * 86400, rotational_amplitude(truth, 35.0, 40.0, u), 0.02)
@printf("rotational amplitude %.3f mag\n", amp.amplitude)

section("2. the pole known")
# The starting shape: two equal lobes, far from the truth; the period and the pole are fixed.
body.shape = BilobedEllipsoid(Ellipsoid(110.0, 95.0, 85.0), Ellipsoid(110.0, 95.0, 85.0), 200.0)
r1 = fit(LimbProblem(evs; model=ContactBinaryEllipsoids(), amplitudes=[amp], fix=[:pole_ra, :pole_dec]); rng=rng)
display(r1)
println("fitted vs true:")
for (n, tv) in ((:a1, 130.0), (:b1, 100.0), (:c1, 85.0), (:a2, 100.0), (:b2, 85.0), (:c2, 75.0), (:neck, 0.1),
                (:volume_equivalent_radius, volume_equivalent_radius(truth)),
                (:mass_ratio, volume(truth.lobe2) / volume(truth.lobe1)))
    @printf("  %-26s %8.3f ± %6.3f   (truth %8.3f)\n", n, r1[n]..., tv)
end

section("3. the pole unknown")
# The pole free, from a grid of starting poles (the fit is non-linear in the pole); the best
# solution is kept.
sols = map([(ra, dec) for ra in (0.0, 90.0, 180.0, 270.0) for dec in (-30.0, 30.0)]) do (pra, pdec)
    p = LimbProblem(evs; model=ContactBinaryEllipsoids(), amplitudes=[amp],
                    init=merge(Dict(zip(r1.names, r1.values)), Dict(:pole_ra => pra, :pole_dec => pdec)),
                    bounds=Dict(:pole_dec => (-90.0, 90.0)))
    r = fit(p; starts=2, rng=rng)
    @printf("  start (%5.1f°, %+5.1f°): pole (%6.1f°, %+5.1f°), χ² = %8.2f\n", pra, pdec, r[:pole_ra][1], r[:pole_dec][1], r.chi2)
    r
end
# Most starts end far away (χ² of thousands); one other start finds a second, worse pole
# (Δχ² ≈ 14): with four events the pole is determined, but a grid of starts is needed to find it.
r2 = sols[argmin([r.chi2 for r in sols])]
@printf("best: pole (%.1f ± %.1f°, %.1f ± %.1f°), truth (35°, 40°); χ² = %.2f for %d degrees of freedom\n",
        r2[:pole_ra]..., r2[:pole_dec]..., r2.chi2, r2.dof)
for (n, tv) in ((:a1, 130.0), (:b1, 100.0), (:c1, 85.0), (:a2, 100.0), (:b2, 85.0), (:c2, 75.0), (:neck, 0.1))
    @printf("  %-6s %8.3f ± %6.3f   (truth %7.3f)\n", n, r2[n]..., tv)
end

section("4. shape model and 3-D views")
sh = Occultations._shape(r1.problem.model, r1.values, r1.problem.setup)
write_obj(joinpath(out, "simcb3_fitted.obj"), ellipsoid_mesh(sh; name="SimCB3 fitted contact binary"))
for k in 1:4
    savefig(plot(r1; event=k, title="SimCB3 event $k", size=(1000, 750)), joinpath(out, "contact_binary_3d_event$k.png"))
    CairoMakie.save(joinpath(out, "contact_binary_3d_occview_event$k.png"), occview(r1; event=k))
end
println("figures written to ", out)

# The four events with the fitted silhouettes (the 1σ shape band from the covariance):
#
#md # ![event 1](../assets/figures/contact_binary_3d_event1.png)
#md # ![event 2](../assets/figures/contact_binary_3d_event2.png)
#
# The fitted body in the geometry of the first event, with its rotational light curve:
#
#md # ![occview](../assets/figures/contact_binary_3d_occview_event1.png)
