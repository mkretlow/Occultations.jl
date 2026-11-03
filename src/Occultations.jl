"""
    Occultations

Reduction and analysis of stellar occultations by solar-system bodies: light-curve
timing with diffraction and atmosphere models, limb fits of 2-D ellipses, 3-D
ellipsoids and shape meshes, binary systems, contact binaries and ring systems, over one or
several events (with negative chords, clock
offsets and rotational-amplitude constraints), uncertainties by covariance, χ²
profiles and Bayesian posterior sampling, astrometry, predictions and maps.

Layers:

* astrometry (`Epoch`, `Star`, ephemerides, observers, fundamental-plane geometry)
* bodies (frames, 3-D shapes, projected `Outline`s)
* photometry (`LightCurve`, diffraction, `Atmosphere`, timing fits)
* events (`Chord`, `Event`, cached `Track`s)
* limb models (`EllipseLimb`, `CircleLimb`, `EllipsoidLimb`, `MeshLimb`, `BinaryLimb`,
  `ContactBinaryLimb`, `ContactBinaryEllipsoids`, `RingedLimb`) and `LimbProblem`
* a generic inference engine (`fit`, `sample`, `profile`) used by all problems
"""
module Occultations

using Dates
using LinearAlgebra
import LinearAlgebra: normalize
using Printf
using Random
using Statistics
using ERFA
using SpecialFunctions
using RecipesBase
import ForwardDiff

# Julia ≥ 1.13 has `Base.ispositive`/`Base.isnegative` for numbers: extend them
@static if isdefined(Base, :ispositive)
    import Base: ispositive
end
@static if isdefined(Base, :isnegative)
    import Base: isnegative
end

# core
include("core/constants.jl")
include("core/time.jl")
include("core/angles.jl")

# astrometry and bodies
"An observer: a ground-based site ([`Observer`](@ref)) or a [`Spacecraft`](@ref)."
abstract type AbstractObserver end
"A 3-D body shape in the body frame (z = rotation axis, x = prime meridian), km."
abstract type AbstractShape end
include("body/rotation.jl")
include("body/outline.jl")
include("body/composite.jl")
include("astro/star.jl")
include("astro/ephemeris.jl")
include("astro/observer.jl")
include("astro/geometry.jl")
include("astro/prediction.jl")
include("body/shape.jl")
include("body/bilobed.jl")

# inference engine
include("fit/params.jl")
include("fit/solve.jl")
include("fit/bayes.jl")

# photometry
include("photometry/diffraction.jl")
include("photometry/lightcurve.jl")
include("photometry/atmosphere.jl")
include("photometry/timing.jl")

# events and limb fits
include("event/chord.jl")
include("event/track.jl")
include("event/event.jl")
include("event/cora_report.jl")
include("limb/models.jl")
include("limb/binary.jl")
include("limb/rings.jl")
include("limb/problem.jl")
include("limb/results.jl")

# plots and maps
include("plots/maps.jl")
include("plots/recipes.jl")

"""
    geomap(x; kwargs...) -> Makie Figure

Publication-quality map of an occultation path with GeoMakie (map projections via PROJ,
Natural Earth coastlines and land, as cartopy in Python). `x` is an [`Event`](@ref), a limb
fit or an [`OccMap`](@ref). Display or save the figure with a Makie backend, e.g.
`import CairoMakie; CairoMakie.save("map.png", geomap(fit))`. See `geomap(::OccMap)` for the
keywords.
"""
function geomap end

"""
    meshview(shape; kwargs...) -> Makie Figure
    meshview(body, t; observer=GEOCENTER, kwargs...)

3-D view of a shape model ([`TriMesh`](@ref), e.g. from [`read_obj`](@ref), or an
[`Ellipsoid`](@ref)), in the body frame with its axes, or – for a body with a frame and a
shape at epoch `t` – as seen by the observer, together with its projected limb. Requires a
Makie backend: `using GLMakie` (interactive: rotate and zoom with the mouse) or
`using CairoMakie` (static images). See the documentation of its methods for the keywords.
"""
function meshview end

"""
    occview(fit_or_posterior; event=1, explore=false, ndraw=24, kwargs...) -> Makie Figure

The occultation in 3-D for a limb fit with a 3-D model (`EllipsoidLimb`, `MeshLimb`) or
its posterior. The figure has four parts:

* the body as oriented at the event, seen from the Earth, lit by the Sun (night side
  dark), with its pole, the Sun direction, the sub-observer and sub-solar points, and the
  chords in front of it;
* the silhouette on the sky plane with the chords; for a posterior, `ndraw` sampled
  silhouettes and shapes as well;
* the rotational light curve: the magnitude from the projected area through one rotation,
  with the phase of the event and the observed amplitudes;
* sliders for the rotational phase and the event. With `explore = true` there are also
  sliders for the pole, W0 and the size, with the live χ² of the problem.

Static with CairoMakie, interactive with GLMakie. See also [`occanimate`](@ref).
"""
function occview end

"""
    occanimate(fit_or_posterior, file; event=1, frames=120, framerate=24, kwargs...) -> file

Record [`occview`](@ref) through one full rotation as a video or animated GIF (extension
of `file`: `.mp4`, `.gif`, `.webm`, …).
"""
function occanimate end

include("plots/geomap.jl")         # module GeoMaps: geomap (GeoMakie)
include("plots/meshview.jl")       # module MeshViews: meshview (Makie)
include("plots/occview.jl")        # module OccViews: occview, occanimate (Makie)

# time and angles
export Epoch, isot, tdb_jd, tt_jd, utc_jd, julian_year, parse_ra, parse_dec, format_ra, format_dec,
       angular_separation, position_angle
# stars, ephemerides, bodies, observers
export Star, geocentric, barycentric, position_error, apparent_diameter, star_diameter_kervella,
       star_diameter_van_belle
export AbstractEphemeris, EphemTable, EphemKeplerian, EphemKernel, EphemPlane, parse_cora, Body, body_gm,
       geocentric_vector, topocentric_vector, ephem_radec, add_offset!, horizons_ephem, query_gaia,
       query_sbdb, hg_magnitude, apparent_magnitude
export AbstractObserver, Observer, Spacecraft, GEOCENTER, gcrs_position, observer_from_mpc,
       parse_mpc_obscode, spacecraft_horizons, altaz, sidereal_time
export positionv, shadow_velocity, closest_approach, shadow_path, ShadowPath, sun_altitude,
       OccultationPrediction, predict
# bodies
export BodyFrame, orientation, orientation_at, body_matrix, pole_position_angle
export AbstractShape, Ellipsoid, TriMesh, read_obj, write_obj, ellipsoid_mesh, volume,
       volume_equivalent_radius, project, projected_axes_dunham, projected_axes_magnusson,
       rotational_amplitude, aspect_angle, BilobedEllipsoid
export Outline, EllipseOutline, PolarOutline, CompositeOutline, intervals, ellipse_outline, limb_radius, signed_distance, inside,
       limb_normal, area, equivalent_radius, outline_points, crossing, algebraic_ellipse, algebraic_circle
# photometry
export LightCurve, Optics, read_lightcurve, read_pymovie, read_tangra, read_astroimagej, write_lightcurve, normalize, select,
       with_optics, epochs, timespan, detect_occultation, detect_occultations, simulate_lightcurve,
       magnitude_drop, fresnel_scale, occultation_flux
export Atmosphere, AtmosphereProfile, atmosphere_profile, atmosphere_flux, half_light_radius,
       xv93_atmosphere, pluto_atmosphere
export TimingProblem, fit_timing, contact_times, model_flux
# events
export CoraReport, read_cora_report, cora_chord
export Contact, Chord, Event, Occultation, Track, EventContext, ispositive, iscomplete, isnegative, occulted_intervals,
       enable, disable, shift,
       replace_chord!, track, velocity, contact_points, chord_path
# limb fits
export LimbModel, EllipseLimb, CircleLimb, EllipsoidLimb, MeshLimb, BinaryLimb, ContactBinaryLimb, ContactBinaryEllipsoids,
       RingedLimb, ProjectedRing, project_ring, ring_geometry,
       LimbProblem, RotationalAmplitude, RotationalPhase, rotational_phase,
       fit_limb, atmosphere_test, outlines, astrometry, contact_residuals, clock_offsets, check_negative_chords,
       predicted_contacts, model_lightcurve
# inference
export AbstractProblem, Param, Prior, Gaussian, Flat, LogPrior, isotropic_pole, parameters, residuals, derived, setparam!, paramnames,
       fit, FitResult, errors, freenames, reduced_chi2, correlation, profile, confidence_interval,
       sample, Posterior, samples, summarize, derived_samples, derived_summary, autocorrelation_time, cornerplot
# maps
export OccMap, occ_map, substellar_point, geomap, meshview, occview, occanimate

end # module
