# Stellar occultation by the dwarf planet (136472) Makemake on 2011 April 23
#
# Data from J. L. Ortiz et al. (2012), "Albedo and atmospheric constraints of dwarf planet
# Makemake from a stellar occultation", Nature 491, 566–569 (doi:10.1038/nature11597), and its
# Supplementary Information (SI):
#   stations and timings: Table 1 (coordinates as in SI Table S1); the circular and elliptical
#   fits: Fig. 2; the star and the shadow velocity: SI Sects. 1 and 3; albedo, atmosphere,
#   density: main text.
# Ephemeris: JPL Horizons (default solution), stored in makemake_horizons_20110423.txt.
#
# Notes on the transcription:
# * The star (NOMAD 1181-0235723, mR = 18.22): the SI gives its J2000 position from La Hita
#   astrometry (dispersion 0.13″ in RA, 0.06″ in Dec). It is Gaia DR3 4010602734012333056
#   (G = 18.22), used here with its proper motion (retrieved 2026-11-02); the SI position is
#   0.10″ from it.
# * The NTT observed in two channels (r′ and g′) of ULTRACAM, San Pedro de Atacama with two
#   telescopes: eight chords from five sites.
# * Negative chords: the eight negative stations of SI Table S1 have no data (clouds, technical
#   or instrument problems, acquisition started too late): none is included.
# * The sensitivity of a chord is its integration plus dead time (Table 1).
#
# `makemake_event()` returns the `Event`; `MAKEMAKE_PAPER` the published results.

using Occultations

const MAKEMAKE_DIR = @__DIR__

"Published results of Ortiz et al. (2012)."
const MAKEMAKE_PAPER = (
    # Fig. 2: circle of diameter 38.28 ± 0.22 mas = 1430 ± 9 km (χ²/dof 1.032); ellipse with a minor
    # axis of 1428 ± 17 km, axial ratio 1.15 ± 0.17, long axis tilted 9 ± 24° from north (χ²/dof
    # 1.027); preferred projected axes 1430 ± 9 km and 1502 ± 45 km (abstract)
    circle_D=(1430.0, 9.0), circle_chi2_pdf=1.032,
    ellipse_minor=(1428.0, 17.0), axial_ratio=(1.15, 0.17), major_axis_pa=(9.0, 24.0), ellipse_chi2_pdf=1.027,
    preferred_axes=((1430.0, 9.0), (1502.0, 45.0)),
    albedo=(0.77, 0.03), density=(1.7, 0.3), pressure_limit_nbar=(4.0, 12.0),
    distance_au=51.5,
    # SI: star position (J2000, La Hita astrometry), velocity of Makemake relative to the star
    star_si=("12 36 11.4018", "+28 11 10.448"), velocity=22.12,
)

# station, longitude (E), latitude, altitude [m], integration + dead time [s],
# immersion, σ, emersion, σ (UTC, 2011-04-23)
const MAKEMAKE_CHORDS = [
    ("La Silla NTT ULTRACAM r'", "289 15 58.5", "-29 15 31.8", 2345.4, 0.2756, "01:35:44.59", 0.07, "01:36:43.51", 0.08),
    ("La Silla NTT ULTRACAM g'", "289 15 58.5", "-29 15 31.8", 2345.4, 0.2756, "01:35:44.64", 0.04, "01:36:43.66", 0.07),
    ("La Silla TRAPPIST", "289 15 38.2", "-29 15 16.6", 2317.7, 6.435, "01:35:46.82", 1.6, "01:36:45.47", 1.6),
    ("Paranal VLT ISAAC", "289 35 50.1", "-24 37 30.3", 2635.0, 1.521, "01:35:46.00", 0.35, "01:36:49.60", 0.35),
    ("Armazones 0.84m", "289 48 13.6", "-24 35 51.9", 2705.7, 13.5, "01:35:46.30", 1.1, "01:36:48.52", 3.0),
    ("San Pedro de Atacama Harlingten", "291 49 13.0", "-22 57 12.2", 2305.0, 6.048, "01:35:37.86", 2.7, "01:36:43.56", 3.1),
    ("San Pedro de Atacama ASH2", "291 49 13.0", "-22 57 12.2", 2305.0, 20.966, "01:35:38.66", 4.0, "01:36:41.16", 2.0),
    ("Pico dos Dias", "314 25 02.5", "-22 32 07.8", 1810.0, 8.851, "01:33:57.27", 1.6, "01:35:01.08", 2.2),
]

"Gaia DR3 4010602734012333056 (= NOMAD 1181-0235723; catalogue epoch 2016.0, with proper motion)."
makemake_star() = Star(189.0474701418674, 28.186223580978275; name="Gaia DR3 4010602734012333056",
                       source_id="4010602734012333056", epoch=2016.0, pmra=-4.560498208477412,
                       pmdec=-10.270692813110907, parallax=0.18970045142629222, ra_err=0.1376, dec_err=0.0778,
                       pmra_err=0.1671, pmdec_err=0.1376, parallax_err=0.1366, mag=Dict("G" => 18.217, "R" => 18.22))

"Ephemeris of Makemake from the stored Horizons table."
function makemake_ephemeris()
    rows = [split(l) for l in eachline(joinpath(MAKEMAKE_DIR, "makemake_horizons_20110423.txt")) if !startswith(l, "#")]
    return EphemTable([r[1] for r in rows], parse.(Float64, getindex.(rows, 2)), parse.(Float64, getindex.(rows, 3)),
                      parse.(Float64, getindex.(rows, 4)); name="(136472) Makemake")
end

"""
    makemake_event() -> Event

The occultation of 2011 April 23 with the eight chords of Table 1 (five sites).
"""
function makemake_event()
    body = Body("(136472) Makemake", makemake_ephemeris(); radius=715.0)
    ev = Event(makemake_star(), body, "2011-04-23T01:36:00"; window=1200, name="Makemake, 2011-04-23")
    day = "2011-04-23T"
    for (name, lon, lat, alt, sens, ti, σi, te, σe) in MAKEMAKE_CHORDS
        push!(ev, Chord(name, Observer(name, lon, lat, alt); immersion=day * ti, immersion_err=σi, emersion=day * te,
                        emersion_err=σe, sensitivity=sens))
    end
    return ev
end
