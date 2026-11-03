# Stellar occultations by the centaur (2060) Chiron on 2018 November 28 and 2019 September 8
#
# Data from F. Braga-Ribas et al. (2023), "Constraints on (2060) Chiron's size, shape, and surrounding
# material from the November 2018 and September 2019 stellar occultations", A&A 676, A72
# (doi:10.1051/0004-6361/202346749):
#   2019: star: Sect. 2; stations: Table 1; timings: Table 2; limb: Table 3; 3-D shape: Sect. 2.2;
#   2018: star: Sect. 3; stations: Table 4; timings: Table 5; astrometry: Table 7.
# J. L. Ortiz et al. (2023), "Changing material around (2060) Chiron revealed by an occultation on
# December 15, 2022" (A&A 676, L12, doi:10.1051/0004-6361/202347025), observed the
# event of 2022 December 15 (one chord of the body at Wise Observatory, extinction features at
# Kottamia); its timings are not tabulated, so that event is not included.
# Ephemeris: JPL Horizons (default solution), stored in chiron_horizons_<yyyymmdd>.txt.
#
# Notes on the transcription:
# * The stars are the paper's Gaia DR3 positions at the events.
# * Paris (2019) had no reliable time source; the paper shifts its chord by +0.9 s.
# * Observed intervals of the negative chords are not tabulated: ±10 min around the event.
#   ATOM (2018) started after the local closest approach and is left out.
# * The sensitivity of a chord is its cycle time.
#
# `chiron_event(year)` returns the `Event` of 2018 or 2019; `CHIRON_PAPER` the published results.

using Occultations

const CHIRON_DIR = @__DIR__

"Published results of Braga-Ribas et al. (2023)."
const CHIRON_PAPER = (
    # Table 3 (R_equiv from Lellouch et al. 2017, imposed: 98–111 km); the 3-D Jacobi ellipsoid
    limb=(R_equiv=(105.0, 7.0, 6.0), a=(126.0, 22.0), position_angle=(11.0, 33.0), oblateness=(0.27, 0.18)),
    axes_3d=((126.0, 22.0), (109.0, 19.0), (68.0, 13.0)), paris_shift=0.9, velocity_2019=22.8, velocity_2018=5.77,
    chord_2018=(180.5, 0.9), distance_2019_au=17.85, distance_2018_au=18.3538,
    # Table 7: astrometry
    positions=Dict(2018 => ("2018-11-28T20:50:43.14", "23 46 04.35419", "+02 13 05.25822"),
                   2019 => ("2019-09-08T23:04:26.20", "00 10 12.73280", "+04 37 05.26207")),
)

# station, longitude (E), latitude, altitude [m], cycle [s], immersion, σ, emersion, σ (UTC)
const CHIRON_CHORDS = Dict(
    2019 => [
        ("Dourbes", "04 34 56.0", "50 05 25.9", 195.0, 1.28, "23:04:12.10", 0.42, "23:04:17.79", 0.70),
        ("Nandrin", "05 26 29.5", "50 31 24.8", 261.0, 2.56, "23:04:08.00", 0.43, "23:04:13.40", 0.92),
        ("Ludiver", "-01 43 42.1", "49 37 46.9", 179.0, 0.50, "23:04:27.96", 0.22, "23:04:33.94", 0.42),
        ("Paris", "02 20 26.5", "48 49 26.9", 79.0, 1.50, "23:04:20.39", 0.56, "23:04:25.47", 0.79),
        ("Forcarei", "-08 22 15.2", "42 36 38.4", 670.0, 2.56, nothing, NaN, nothing, NaN),
    ],
    2018 => [
        ("SAAO", "20 48 37.8", "-32 22 41.9", 1760.0, 1.5, "20:50:12.51", 0.06, "20:50:43.80", 0.10),
        ("Boyden", "26 24 17.0", "-29 02 19.8", 1372.0, 10.0, nothing, NaN, nothing, NaN),
    ],
)

"The occulted star (Gaia DR3 at the event, as given in the paper)."
chiron_star(year::Integer) = year == 2019 ?
    Star("00 10 12.74303", "+04 37 04.9141"; name="Gaia DR3 2741744913738549760", source_id="2741744913738549760",
         epoch=julian_year(Epoch("2019-09-08T23:04")), ra_err=0.16, dec_err=0.28, mag=Dict("G" => 16.52)) :
    Star("23 46 04.33575", "+02 13 05.5270"; name="Gaia DR3 2646598228351156352", source_id="2646598228351156352",
         epoch=julian_year(Epoch("2018-11-28T20:50")), ra_err=0.25, dec_err=0.34, mag=Dict("G" => 17.28))

"Ephemeris of Chiron for 2018 or 2019 from the stored Horizons tables."
function chiron_ephemeris(year::Integer)
    file = year == 2019 ? "chiron_horizons_20190908.txt" : "chiron_horizons_20181128.txt"
    rows = [split(l) for l in eachline(joinpath(CHIRON_DIR, file)) if !startswith(l, "#")]
    return EphemTable([r[1] for r in rows], parse.(Float64, getindex.(rows, 2)), parse.(Float64, getindex.(rows, 3)),
                      parse.(Float64, getindex.(rows, 4)); name="(2060) Chiron")
end

"""
    chiron_event(year) -> Event

The occultation of 2019 September 8 (four positive chords, Forcarei negative) or of
2018 November 28 (SAAO positive, Boyden negative).
"""
function chiron_event(year::Integer)
    year in (2018, 2019) || throw(ArgumentError("the events are 2018 and 2019"))
    body = Body("(2060) Chiron", chiron_ephemeris(year); radius=105.0)
    date, tc = year == 2019 ? ("2019-09-08", "23:05:00") : ("2018-11-28", "20:51:00")
    ev = Event(chiron_star(year), body, date * "T" * tc; window=1200, name="Chiron, " * date)
    day = date * "T"
    span = year == 2019 ? ("22:55:00", "23:15:00") : ("20:40:00", "21:00:00")
    for (n, lon, lat, alt, sens, ti, σi, te, σe) in CHIRON_CHORDS[year]
        obs = Observer(n, lon, lat, alt)
        push!(ev, ti === nothing ? Chord(n, obs; span=(day * span[1], day * span[2]), sensitivity=sens) :
                  Chord(n, obs; immersion=day * ti, immersion_err=σi, emersion=day * te, emersion_err=σe,
                        sensitivity=sens))
    end
    return ev
end
