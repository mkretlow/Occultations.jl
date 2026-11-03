# Stellar occultation by 2002 TC302 on 2018 January 28
#
# Data from J. L. Ortiz et al. (2020), "The large trans-Neptunian object 2002 TC302 from
# combined stellar occultation, photometry, and astrometry data", A&A 639, A134
# (doi:10.1051/0004-6361/202038046):
#   star: Table 1; stations: Table 2; timings: Table 3; limb fit: Sect. 5; rotation: Sect. 6.
# Ephemeris: JPL Horizons (default solution), stored in 2002tc302_horizons_20180128.txt; the
# paper does not give the fitted centre.
#
# Notes on the transcription:
# * The star: Table 1 gives Gaia DR2 130957817758443648 (the abstract 130957813463146112, which
#   is not in Gaia DR2) at "J2000 coordinates": its position is the Gaia DR2 position
#   propagated to the epoch J2000.0 (the differences, −0.55″ in RA and +0.11″ in Dec, are 15.5
#   years of proper motion). The star is therefore taken from Gaia DR2 itself (epoch 2015.5,
#   with proper motion and parallax, retrieved 2026-11-02); at the event it was 0.64″
#   (≈ 20 000 km at the body) from the Table 1 position.
# * The altitudes of the stations are not given: those of Crni Vrh, Asiago, San Marcello
#   Pistoiese, Mount Agliale and Sabadell are taken from the other examples (same sites),
#   0 m elsewhere (a few hundred metres change the chords by less than 0.5 km).
# * The observed intervals of the near misses are not tabulated; the observations started
#   ~15 min before and ended ~15 min after the event (Sect. 3): 21:38–22:08 UTC is used.
# * The sensitivity of a chord is its cycle time (Table 2).
# * Tavolaia had no time stamps in the images (file times of the operating system, Sect. 4):
#   ± 2.5 s in Table 3.
#
# `tc302_event()` returns the `Event`; `TC302_PAPER` the published results.

using Occultations

const TC302_DIR = @__DIR__

"Published results of Ortiz et al. (2020)."
const TC302_PAPER = (
    # Sect. 5: the ellipse, axes 543.2 ± 18 km and 459.5 ± 11 km (full axes), position angle of
    # the minor axis 3 ± 1°; equivalent diameter 499.6 km (abstract: 543 × 460 km, √ = 499.8 km)
    a=(543.2 / 2, 18.0 / 2), b=(459.5 / 2, 11.0 / 2), position_angle=(3.0, 1.0), D_eq=499.6,
    # albedo (abstract), absolute magnitude, rotation (Sects. 6–8)
    albedo=(0.147, 0.005), period_h=56.1, amplitude=(0.06, 0.01),
    thermal_diameter=(584.1, 105.0, 88.0),
    velocity=4.77, star_diameter_mas=0.009,
)

# station, longitude (E), latitude, altitude [m], cycle time [s],
# immersion, σ, emersion, σ (UTC, 2018-01-28; `nothing` for a near miss)
const TC302_CHORDS = [
    ("Crni Vrh", 14.071083, 45.945833, 713.0, 4.995, "21:52:51.687", 1.150, "21:53:35.510", 0.350),
    ("Asiago", 11.568806, 45.849444, 1370.0, 8.301, "21:52:43.400", 1.400, "21:54:08.875", 0.200),
    ("S. Marcello Pistoiese", 10.803889, 44.064167, 965.0, 3.933, "21:53:22.837", 0.475, "21:54:54.051", 0.163),
    ("Tavolaia", 10.673306, 43.736500, 0.0, 7.23, "21:53:36.3", 2.5, "21:55:08.7", 2.5),
    ("Mount Agliale", 10.514944, 43.995278, 760.0, 3.718, "21:53:23.150", 0.300, "21:54:57.450", 0.100),
    ("La Spezia", 9.853528, 44.126278, 0.0, 7.830, "21:53:18.850", 0.500, "21:54:55.550", 2.400),
    ("Gnosca", 9.024028, 46.231444, 0.0, 2.568, "21:52:37.170", 1.250, "21:54:21.850", 1.250),
    ("Varese, Schiaparelli Obs.", 8.770278, 45.868056, 0.0, 5.716, "21:52:42.585", 0.100, "21:54:24.775", 0.100),
    ("OCA, Nice", 7.299833, 43.725806, 0.0, 2.2005, "21:53:39.661", 0.280, "21:55:23.301", 0.720),
    ("Aosta Valley", 7.478333, 45.789444, 0.0, 2.695, "21:52:49.240", 0.440, "21:54:31.565", 0.420),
    ("Biot, Nice", 7.077778, 43.617222, 0.0, 26.463, "21:53:46.800", 5.500, "21:55:24.750", 3.500),
    ("Vinon sur Verdon", 5.796111, 43.737778, 0.0, 30.467, "21:53:52.08", 6.70, "21:55:21.13", 6.00),
    ("Salvatore Di Giacomo", 14.564056, 40.623944, 0.0, 3.4, nothing, NaN, nothing, NaN),
    ("Agerola", 14.571556, 40.626083, 0.0, 2.5, nothing, NaN, nothing, NaN),
    ("Sabadell", 2.090167, 41.550000, 224.0, 2.56, nothing, NaN, nothing, NaN),
    ("San Esteve", 1.872528, 41.493750, 0.0, 2.0, nothing, NaN, nothing, NaN),
]

"Gaia DR2 130957817758443648 (catalogue epoch 2015.5, with proper motion and parallax)."
tc302_star() = Star(35.45594509618496, 28.403702753684115; name="Gaia DR2 130957817758443648",
                    source_id="130957817758443648", epoch=2015.5, pmra=35.33722155018179,
                    pmdec=-7.059985509111287, parallax=1.4208908869002408, ra_err=0.0578, dec_err=0.0493,
                    pmra_err=0.0956, pmdec_err=0.0894, parallax_err=0.0787, mag=Dict("G" => 15.73),
                    diameter=TC302_PAPER.star_diameter_mas)

"Ephemeris of 2002 TC302 from the stored Horizons table."
function tc302_ephemeris()
    rows = [split(l) for l in eachline(joinpath(TC302_DIR, "2002tc302_horizons_20180128.txt")) if !startswith(l, "#")]
    return EphemTable([r[1] for r in rows], parse.(Float64, getindex.(rows, 2)), parse.(Float64, getindex.(rows, 3)),
                      parse.(Float64, getindex.(rows, 4)); name="2002 TC302")
end

"""
    tc302_event(; negatives=:all) -> Event

The occultation of 2018 January 28 with the twelve positive chords (Table 3) and the four
near misses of Table 2 (`negatives = :all` or `:none`).
"""
function tc302_event(; negatives::Symbol=:all)
    negatives in (:all, :none) || throw(ArgumentError("negatives must be :all or :none"))
    body = Body("2002 TC302", tc302_ephemeris(); radius=250.0)
    ev = Event(tc302_star(), body, "2018-01-28T21:54:00"; window=1200, name="2002 TC302, 2018-01-28")
    day = "2018-01-28T"
    for (name, lon, lat, alt, sens, ti, σi, te, σe) in TC302_CHORDS
        obs = Observer(name, lon, lat, alt)
        if ti === nothing
            negatives === :all || continue
            push!(ev, Chord(name, obs; span=(day * "21:38:00", day * "22:08:00"), sensitivity=sens))
        else
            push!(ev, Chord(name, obs; immersion=day * ti, immersion_err=σi, emersion=day * te, emersion_err=σe,
                            sensitivity=sens))
        end
    end
    return ev
end
