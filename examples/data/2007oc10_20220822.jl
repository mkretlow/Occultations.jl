# Stellar occultation by (470316) 2007 OC10 on 2022 August 22: a double star
#
# Data from J. M. Gómez-Limón et al. (2025), "Size and shape of the trans-Neptunian object
# (470316) 2007 OC10: Comparison with thermal data", A&A 697, A157
# (doi:10.1051/0004-6361/202453540): the star (Table 1), the stations with their exposure and
# cycle times (Table C.1), the results (Tables 3 and 5) and the light curves of 13 stations at
# the CDS (J/A+A/697/A157), stored in 2007oc10_cds/ (time in s after 2022-08-22T00:00 UTC at
# mid-exposure, normalised flux and its error).
# Ephemeris: JPL Horizons (default solution), stored in 2007oc10_horizons_20220822.txt; the paper
# refers its fitted centre to an ephemeris updated with astrometry from Calar Alto.
#
# Notes on the transcription:
# * Stations without data (Cancelada, Forcarei, Sabadell: clouds) and with unusable data (OAO:
#   300 s exposures; LCO Tenerife: drift scan, low S/N) are omitted, as in the paper.
# * Botorrita has a clock drift of about 3 s and Çanakkale time stamps truncated to whole
#   seconds; both are negative, and the paper uses them as they are.
# * The light curves are cut to ±`OC10_WINDOW` s around the predicted geocentric closest approach.
#
# `oc10_event()` returns the `Event` with one light-curve chord per station; `OC10_PAPER` the
# published results.

using Occultations

const OC10_DIR = @__DIR__
const OC10_WINDOW = 300.0                   # s around the predicted closest approach (the drops lie 100–170 s before it)
const OC10_CA = "2022-08-22T01:32:37"       # predicted closest approach (JPL#12, CORA)

"Published results of Gómez-Limón et al. (2025), sampling 1 (uniform priors), Tables 3 and 5."
const OC10_PAPER = (
    # the ellipse: semi-major axis, axis ratio, tilt ϕ from east towards north to the major axis;
    # posterior medians and 68 % credible intervals (−, +)
    a=(215.0, 7.0, 10.0), b_over_a=(0.58, 0.16, 0.16), tilt=(16.0, 11.0, 12.0), D_eq=(330.0, 55.0, 56.0),
    # the centre relative to the paper's reference ephemeris (CAHA astrometry) [km]
    xc=(-653.0, 33.0, 66.0), yc=(98.0, 180.0, 58.0),
    # the double star: flux fraction of the primary, its offset from the photocentre [km],
    # separation [mas], position angle of the secondary from the primary [deg], Δmag
    f1=(0.748, 0.012, 0.011), x1=(-305.0, 38.0, 31.0), y1=(-209.0, 146.0, 22.0),
    sep_mas=(57.0, 11.0, 4.0), pa=(56.0, 17.0, 3.0), dmag=(1.18, 0.07, 0.07),
    # Table 5: astrometric positions at the reference epoch
    astrometry=(epoch="2022-08-22T01:32:37", ra="22 08 33.0652", dec="+12 05 55.7402",
                ra1="22 08 33.0672", dec1="+12 05 55.6714", ra2="22 08 33.0703", dec2="+12 05 55.7035"),
    # geometry (Sect. 3.1)
    distance_au=35.532, velocity=23.6, fresnel_km=1.3,
)

# station, CDS file, latitude [deg], longitude (E) [deg], altitude [m], exposure [s], cycle [s]
# (Table C.1)
const OC10_STATIONS = [
    ("Albox",          "albox",    37.405575,  -2.151800,  538, 3.0,   3.193),
    ("Botorrita",      "botorrit", 41.497375,  -1.020867,  456, 0.5,   0.5007),
    ("El Guijo",       "el_guijo", 40.607537,  -4.031010,  968, 1.0,   1.002),
    ("Guirguillano",   "guirguil", 42.711764,  -1.865110,  645, 10.24, 10.24),
    ("Çanakkale",      "canakkle", 40.098990,  26.474490,  437, 3.0,   3.726),
    ("Liverpool",      "liverpoo", 28.762574, -17.879200, 2333, 1.183, 1.223),
    ("Javalambre",     "javalamb", 40.042464,  -1.016482, 1988, 4.0,   4.444),
    ("Oukaimeden",     "oukaimed", 31.206302,  -7.866467, 2779, 2.0,   5.0),
    ("Wise",           "wise",     30.596673,  34.762208,  883, 1.0,   1.003),
    ("Artemis",        "artemis",  28.299985, -16.511551, 2437, 0.3,   1.822),
    ("GAL Hassin",     "gal_hass", 37.939377,  14.020610,  648, 9.0,   10.339),
    ("SLN",            "sln",      37.692467,  14.973149, 1729, 1.0,   1.244),
    ("TAR2",           "tar2",     28.298441, -16.510282, 2416, 1.0,   1.02),
]

"Gaia DR3 2727866328215869952 (Table 1; catalogue epoch 2016.0, with proper motion and parallax)."
oc10_star() = Star("22 08 33.066148", "+12 05 55.694318"; name="Gaia DR3 2727866328215869952",
                   source_id="2727866328215869952", epoch=2016.0, pmra=4.00, pmdec=-2.25, parallax=0.43,
                   ra_err=0.05, dec_err=0.05, pmra_err=0.06, pmdec_err=0.06, parallax_err=0.06,
                   mag=Dict("V" => 15.23))

"Ephemeris of 2007 OC10 from the stored Horizons table."
function oc10_ephemeris()
    rows = [split(l) for l in eachline(joinpath(OC10_DIR, "2007oc10_horizons_20220822.txt")) if !startswith(l, "#")]
    return EphemTable([r[1] for r in rows], parse.(Float64, getindex.(rows, 2)), parse.(Float64, getindex.(rows, 3)),
                      parse.(Float64, getindex.(rows, 4)); name="(470316) 2007 OC10")
end

"The light curve of a station from the CDS file, cut to ±`window` s around the closest approach."
function oc10_lightcurve(file, name, texp; window=OC10_WINDOW)
    rows = [split(l, '\t') for l in eachline(joinpath(OC10_DIR, "2007oc10_cds", file * ".tsv"))
            if !startswith(l, "#") && occursin(r"^\s*\d", l)]
    t = parse.(Float64, strip.(getindex.(rows, 2)))         # s after 2022-08-22T00:00 UTC, mid-exposure
    f = parse.(Float64, strip.(getindex.(rows, 3)))
    e = parse.(Float64, strip.(getindex.(rows, 4)))
    tref = Epoch("2022-08-22T00:00:00")
    tca = Epoch(OC10_CA) - tref
    k = findall(x -> abs(x - tca) <= window, t)
    return LightCurve(t[k], f[k]; tref=tref, σ=e[k], exptime=texp, λ=0.6, Δλ=0.3, name=name)
end

"""
    oc10_event(; window=OC10_WINDOW) -> Event

The occultation of 2022 August 22 with one light-curve chord for each of the 13 stations.
"""
function oc10_event(; window=OC10_WINDOW)
    body = Body("(470316) 2007 OC10", oc10_ephemeris(); radius=165.0)
    ev = Event(oc10_star(), body, OC10_CA; window=1200, name="2007 OC10, 2022-08-22")
    for (name, file, lat, lon, alt, texp, _) in OC10_STATIONS
        push!(ev, Chord(name, Observer(name, lon, lat, alt), oc10_lightcurve(file, name, texp; window=window)))
    end
    return ev
end
