# Stellar occultation by (229762) 2007 UK126 on 2014 November 15
#
# Data from K. Schindler et al. (2017), "Results from a triple chord stellar occultation and
# far-infrared photometry of the trans-Neptunian object (229762) 2007 UK126", A&A 600, A12
# (doi:10.1051/0004-6361/201628620):
#   stations and timings: Table 2; ellipse: Table 5; albedo: Sect. 4.3; Maclaurin: Table 6.
# Star: UCAC4 448-006503 in the paper (J2000 04 29 30.6, −00 28 20.9); the Gaia DR3
# counterpart 3206533286529807488 from the Gaia archive is used here. Ephemeris: JPL Horizons
# (default solution), stored in uk126_horizons_20141115.txt.
#
# Notes on the transcription:
# * The paper's "axes" a_El and b_El are full axes, not semi-axes: the central chord alone is
#   641 km long, and the albedo of Sect. 4.3 (15.0 %) follows from a_El · b_El as full axes.
# * θ_El is the angle of the major axis from east towards north.
# * The sensitivity of a chord is its integration time. 16 stations reported; only the three
#   successful ones are tabulated, so no negative chords are included.
#
# `uk126_event()` returns the `Event`; `UK126_PAPER` the published results.

using Occultations

const UK126_DIR = @__DIR__

"Published results of Schindler et al. (2017)."
const UK126_PAPER = (
    # Table 5 (full axes; centre relative to the authors' ephemeris; θ of the major axis from E to N)
    axis_major=(645.80, 5.68), axis_minor=(597.81, 12.74), center_x=(6.57, 2.06), center_y=(-10.82, 3.13),
    theta=(21.25, 5.65), axis_ratio=(1.080, 0.025),
    # Table 2: durations [s] and chord lengths [km]
    durations=Dict("Schindler & Wolf" => (25.893, 0.225), "Olsen" => (21.355, 1.459), "Bardecker" => (19.219, 0.639)),
    chords=Dict("Schindler & Wolf" => (641.0, 3.9), "Olsen" => (525.9, 25.4), "Bardecker" => (475.7, 11.2)),
    distance_au=42.572, H_V=(3.69, 0.044), albedo=(0.150, 0.016),
    # Table 6: Maclaurin spheroid, a = b (full axis) and c [km]
    maclaurin=(a=(640.0, 651.0), c=(373.0, 611.0)),
)

# station, longitude (E), latitude, altitude [m], integration time [s], D, σ, R, σ (UTC)
const UK126_CHORDS = [
    ("Schindler & Wolf", "-119 24 45.00", "37 04 13.50", 1405.0, 2.000, "10:19:24.356", 0.159, "10:19:50.249", 0.159),
    ("Olsen", "-88 11 46.30", "40 05 12.40", 224.0, 4.271, "10:18:05.447", 1.032, "10:18:26.802", 1.032),
    ("Bardecker", "-119 40 20.32", "38 53 23.53", 1534.0, 2.135, "10:19:29.151", 0.452, "10:19:48.370", 0.452),
]

"Gaia DR3 3206533286529807488 (= UCAC4 448-006503), from the Gaia archive."
uk126_star() = Star(67.37756526478476, -0.4724628577726251; name="Gaia DR3 3206533286529807488",
                    source_id="3206533286529807488", epoch=2016.0, pmra=-0.7572, pmdec=-12.3297, parallax=1.0220,
                    ra_err=0.0299, dec_err=0.0217, pmra_err=0.0441, pmdec_err=0.0338,
                    mag=Dict("G" => 15.715, "V" => 15.86), diameter=0.01)

"Ephemeris of 2007 UK126 from the stored Horizons table."
function uk126_ephemeris()
    rows = [split(l) for l in eachline(joinpath(UK126_DIR, "2007uk126_horizons_20141115.txt")) if !startswith(l, "#")]
    return EphemTable([r[1] for r in rows], parse.(Float64, getindex.(rows, 2)), parse.(Float64, getindex.(rows, 3)),
                      parse.(Float64, getindex.(rows, 4)); name="(229762) 2007 UK126")
end

"The occultation of 2014 November 15 with the three positive chords of Table 2."
function uk126_event()
    body = Body("(229762) 2007 UK126", uk126_ephemeris(); radius=310.0)
    ev = Event(uk126_star(), body, "2014-11-15T10:19:00"; window=1200, name="2007 UK126, 2014-11-15")
    day = "2014-11-15T"
    for (name, lon, lat, alt, sens, ti, σi, te, σe) in UK126_CHORDS
        push!(ev, Chord(name, Observer(name, lon, lat, alt); immersion=day * ti, immersion_err=σi,
                        emersion=day * te, emersion_err=σe, sensitivity=sens))
    end
    return ev
end
