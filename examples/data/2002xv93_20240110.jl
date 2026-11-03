# Stellar occultation by (612533) 2002 XV93 on 2024 January 10 – an atmosphere
#
# K. Arimatsu et al. (2026), "Detection of an atmosphere on a trans-Neptunian object beyond
# Pluto", Nature Astronomy, doi:10.1038/s41550-026-02846-1 (preprint arXiv:2605.02243):
#   stations: Table 1; chord durations, closest approaches and fits: Methods.
# Light curves and best-fit model profiles: the authors' dataset, Zenodo
#   doi:10.5281/zenodo.18976431 (CC BY 4.0), in 2002xv93_zenodo/ (see its README.md).
# Star: Gaia DR3 940732910252854272 (from the Gaia archive); ephemeris: JPL Horizons (default
# solution), stored in xv93_horizons_20240110.txt; the paper used JPL#18.
#
# `xv93_star()`, `xv93_ephemeris()`, `xv93_lightcurve(station)`, `xv93_model_profile(gas, ρ)`,
# `XV93_STATIONS`, `XV93_PAPER`.

using Occultations

const XV93_DIR = @__DIR__

"Published results of Arimatsu et al. (2026, Nature Astronomy)."
const XV93_PAPER = (
    distance_au=37.0, velocity=25.2, star_diameter_mas=0.0032, density=1500.0,
    durations=Dict("Kyoto" => (18.0, 0.9), "Kiso" => (16.53, 0.08)),               # s
    chord_lengths=Dict("Kyoto" => (458.0, 22.0), "Kiso" => (422.0, 2.0)),          # km
    closest_approach=Dict("Kyoto" => (41.0, 45.0, 41.0), "Kiso" => (103.0, 45.0, 43.0),
                          "Fukushima" => (255.0, 45.0, 43.0)),                     # km (value, +, −)
    fukushima_ca_time="2024-01-10T13:12:50.818",
    circle_radius=(235.0, 22.0, 15.0), circle_chi2=1.40,                           # km (value, +, −)
    radius_atm=(244.0, 1.0),                                                       # km (Kiso fit)
    pressure=Dict(:CH4 => (124.0, 15.0, 12.0), :N2 => (177.0, 18.0, 19.0), :CO => (159.0, 15.0, 17.0)),  # nbar, ρ = 1500
    pressure_rho1000=Dict(:CH4 => (131.0, 11.0, 17.0), :N2 => (181.0, 20.0, 18.0), :CO => (163.0, 17.0, 17.0)),
    chi2_atm=Dict(:CH4 => 28.7, :N2 => 28.5), dof_atm=29, chi2_airless=89.2, dof_airless=30,
    surface_T=47.0, dTdz=2.7, z_inv=20.0, z_fmf=25.0,
)

# station, longitude (E), latitude, altitude [m], exposure [s], status
const XV93_STATIONS = [
    ("Kyoto", "135 46 59.9", "35 01 51.2", 62.0, 1.0, :positive),
    ("Kiso", "137 37 31.5", "35 47 50.0", 1132.0, 0.5, :positive),
    ("Fukushima", "140 26 04.2", "37 25 36.7", 274.0, 0.49, :negative),
]
xv93_station(name) = (s = XV93_STATIONS[findfirst(x -> x[1] == name, XV93_STATIONS)];
                      (observer=Observer(s[1], s[2], s[3], s[4]), exptime=s[5], status=s[6]))

"Gaia DR3 940732910252854272 (Gaia archive), with the angular diameter of the paper."
xv93_star() = Star(104.4072176861, 36.3662884068; name="Gaia DR3 940732910252854272", source_id="940732910252854272",
                   epoch=2016.0, pmra=-0.1855, pmdec=-1.3677, parallax=0.2505, ra_err=0.0362, dec_err=0.0335,
                   pmra_err=0.042, pmdec_err=0.035, mag=Dict("G" => 15.784), diameter=XV93_PAPER.star_diameter_mas)

"Ephemeris of 2002 XV93 from the stored Horizons table."
function xv93_ephemeris()
    rows = [split(l) for l in eachline(joinpath(XV93_DIR, "2002xv93_horizons_20240110.txt")) if !startswith(l, "#")]
    return EphemTable([r[1] for r in rows], parse.(Float64, getindex.(rows, 2)), parse.(Float64, getindex.(rows, 3)),
                      parse.(Float64, getindex.(rows, 4)); name="(612533) 2002 XV93")
end

"Body with the paper's bulk density (for the gravity of the atmosphere model)."
xv93_body(; radius=XV93_PAPER.radius_atm[1], density=XV93_PAPER.density) =
    Body("(612533) 2002 XV93", xv93_ephemeris(); radius=radius, density=density)

"""
    xv93_lightcurve(station; λ=0.6, Δλ=0.3) -> LightCurve

The reduced, normalised light curve of a station (Zenodo dataset), with its exposure time.
"""
function xv93_lightcurve(name::AbstractString; λ::Real=0.6, Δλ::Real=0.3)
    rows = [split(l, ',') for l in eachline(joinpath(XV93_DIR, "2002xv93_zenodo", "lightcurve_$name.csv"))][2:end]
    t = [Epoch(parse(Float64, r[1]) + 2400000.5) for r in rows]
    return LightCurve(t, parse.(Float64, getindex.(rows, 2)); σ=parse.(Float64, getindex.(rows, 3)),
                      exptime=xv93_station(name).exptime, λ=λ, Δλ=Δλ, name=name)
end

"""
    xv93_model_profile(gas, ρ=1500) -> (distance, intensity)

Best-fit refraction profile of the paper (`gas` = `:ch4`, `:n2`, `:co`; bulk density 1500 or
1000 kg/m³): normalised intensity versus distance from the body centre [km].
"""
function xv93_model_profile(gas::Symbol, ρ::Integer=1500)
    rows = [split(l, ',') for l in eachline(joinpath(XV93_DIR, "2002xv93_zenodo", "bestfitmodel_$(lowercase(string(gas)))_rho$ρ.csv"))][2:end]
    return parse.(Float64, strip.(getindex.(rows, 1))), parse.(Float64, strip.(getindex.(rows, 2)))
end
