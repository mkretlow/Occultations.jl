# Interactive 3-D view of the 2003 UY117 occultation (needs GLMakie and a display)
#
# GLMakie is not a dependency of the package: install it in your default environment
# (`julia -e 'using Pkg; Pkg.add("GLMakie")'`), which is stacked on the project environment.
#
# The triaxial ellipsoid fitted to the occultation, the light-curve amplitude and the
# rotational phase (as in examples/2003uy117_2020.jl) is shown in the geometry of the event:
# * drag the 3-D panel to rotate the view, scroll to zoom;
# * the slider "rotational phase" (0–1, 0 = light-curve maximum before the event) turns the
#   body about its spin axis – the silhouette against the chords and the marker on the
#   rotational light curve follow; the phase of the event is shown next to it;
# * the sliders pole RA, pole Dec, W0 (prime-meridian angle at the reference epoch), a, b, c
#   change the solution; the χ² of the problem is shown live; "reset" restores the fit.
# A second window shows the posterior (sampled silhouettes in the sky plane).
#
#   julia --project=. examples/occview_demo.jl

using Occultations, Random
using GLMakie
include(joinpath(@__DIR__, "data", "2003uy117_20201023.jl"))

ev = uy117_event()
P = UY117_PAPER
r2 = fit_limb(ev; residual=:radial, rng=MersenneTwister(1))
amp = RotationalAmplitude(ev.tref, P.amplitude...)
phase = RotationalPhase("2020-10-23T22:21:46", P.rotational_phase, 0.05)
prob = LimbProblem(ev; model=EllipsoidLimb(principal=true), amplitudes=[amp], phases=[phase],
                   init=Dict(:a => 160.0, :b => 110.0, :c => 90.0, :pole_ra => 337.0, :pole_dec => 62.0, :W0 => 30.0,
                             :period => P.period_h, :center_f => r2[:center_f][1], :center_g => r2[:center_g][1]),
                   priors=Dict(:pole_dec => isotropic_pole()),
                   bounds=Dict(:a => (20.0, 400.0), :b => (20.0, 400.0), :c => (20.0, 400.0)))
r = fit(prob; rng=MersenneTwister(3), starts=16)
println("least squares: a = ", round(r[:a][1]; digits=1), ", b = ", round(r[:b][1]; digits=1), ", c = ",
        round(r[:c][1]; digits=1), " km; pole (", round(r[:pole_ra][1]), "°, ", round(r[:pole_dec][1]), "°)")

GLMakie.activate!(title="2003 UY117 – occultation in 3-D (explore)")
s1 = display(GLMakie.Screen(), occview(r; explore=true, title="2003 UY117, 2020-10-23 – least-squares ellipsoid (explore)"))

println("sampling the posterior for the second window …")
post = sample(r; nsteps=3000, nwalkers=32, burn=1000, rng=MersenneTwister(4))
s2 = display(GLMakie.Screen(), occview(post; ndraw=30, title="2003 UY117, 2020-10-23 – posterior"))
println("close both windows to end")
wait(s1); wait(s2)
