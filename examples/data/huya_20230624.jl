# Stellar occultation by (38628) Huya and its satellite on 2023 June 24
#
# Data from F. L. Rommel et al. (2025), "Stellar Occultation Observations of (38628) Huya and
# Its Satellite: A Detailed Look into the System", Planet. Sci. J. 6, 48
# (doi:10.3847/PSJ/adabc1):
#   star: Table 1 (Gaia DR3 position at the closest approach); timings: Table 2;
#   stations: Tables 7 and 8; limb fits: Table 3; astrometry: Tables 4 and 5.
# Ephemeris: JPL Horizons (default solution), stored in huya_horizons_20230624.txt; the paper
# refers the fitted centre to NIMA v10.
#
# Notes on the transcription:
# * The satellite was detected from La Palma only (00:58:26.42–00:58:34.43); for Huya this
#   station is a negative chord.
# * Stations with inconclusive results (Starhopper, Sant Esteve Sesrovires) or without data are
#   omitted. The observed intervals of the negatives are not tabulated; 00:52–01:02 UTC,
#   which brackets the event (00:56:20–00:57:12 at the positive sites), is used.
# * The sensitivity of a chord is its cycle time (exposure + readout, Tables 7, 8).
# * The time offsets of Table 7 (Sabadell −2.54 s, PixelSkies +2.5 s) are assumed to be
#   applied in the published instants.
# * Satellite astrometry (Table 4, 00:58:10.9 UTC): the southern and northern solutions of a
#   circle of radius 106.5 km on the La Palma chord.
#
# `huya2023_event()` returns the `Event`; `HUYA2023_PAPER` the published results.

using Occultations

const HUYA2023_DIR = @__DIR__

"Published results of Rommel et al. (2025)."
const HUYA2023_PAPER = (
    # Table 3 (centre relative to NIMA v10; PA of the minor axis from north)
    free=(center_f=(22.3, 6.7), center_g=(-24.9, 7.9), a=(222.5, 9.1), b=(198.7, 15.2), oblateness=(0.107, 0.058),
          position_angle=(31.5, 15.5), R_equiv=(210.3, 11.0), chi2_pdf=0.682),
    restrict=(center_f=(22.2, 6.6), center_g=(-25.7, 5.0), a=(218.7, 8.1), b=(200.3, 14.9), oblateness=(0.084, 0.059),
              position_angle=(53.7, 2.2), R_equiv=(209.3, 10.3), chi2_pdf=0.745),
    # the 2019 data of Santos-Sanz et al. (2022) refitted with PA = 51.6 ± 2.2° (Appendix B)
    restrict2019=(center_f=(49.9, 0.15), center_g=(26.61, 0.08), a=(218.05, 0.11), b=(195.59, 0.24),
                  oblateness=(0.103, 0.001), position_angle=(51.6, 2.2), R_equiv=(206.5, 0.16), chi2_pdf=210.9),
    npoints=22,
    # Table 4 astrometry (2023 June 24, 00:58:10.9 UTC)
    astrometry=(epoch="2023-06-24T00:58:10.9", ra="17 16 43.653342", dec="-07 00 20.11836", ra_err=0.51, dec_err=0.72),
    # satellite: Table 2 (La Palma), Table 5 (offset from Huya), chord length (Sect. 3.2)
    satellite=(immersion=("00:58:26.42", 0.29), emersion=("00:58:34.43", 0.30), chord_km=(179.0, 14.0),
               offset_arcsec=((0.06546, 0.00300), (-0.04440, 0.00500)), separation_km=1603.0, radius=106.5,
               astrometry=(southern=(ra="17 16 43.657737", dec="-07 00 20.16006", ra_err=0.48, dec_err=0.82),
                           northern=(ra="17 16 43.657738", dec="-07 00 20.16537", ra_err=0.48, dec_err=0.82))),
    velocity=22.43, star_diameter_km=0.18, fresnel_km=1.20,
)

# station, longitude (E), latitude, altitude [m], cycle time [s],
# immersion, σ, emersion, σ (UTC, 2023-06-24; `nothing` for a negative chord)
const HUYA2023_CHORDS = [
    ("Montsec", 0.72965, 42.051655, 1564.582, 14.69, "00:56:25.30", 0.54, "00:56:40.90", 0.55),
    ("Botorrita", -1.020867, 41.497375, 403.0, 4.32, "00:56:29.5", 1.7, "00:56:47.8", 1.1),
    ("Sabadell", 2.09013, 41.550043, 224.0, 5.18, "00:56:21.2", 2.1, "00:56:37.5", 1.3),
    ("Alto de la Vega", -1.9106369708525, 40.4172959742339, 1533.0, 10.37, "00:56:32.8", 1.2, "00:56:51.2", 1.4),
    ("Javalambre", -1.0163, 40.0418, 1957.0, 7.78, "00:56:26.4", 1.0, "00:56:44.7", 0.85),
    ("Linhaceira", -8.3838, 39.522688, 90.0, 12.96, "00:56:55.8", 2.2, "00:57:11.9", 2.4),
    ("La Hita", -3.1833, 39.568, 770.0, 11.21, "00:56:33.4", 2.4, "00:56:54.6", 2.4),
    ("Cala d'Hort", 1.2408, 38.891102, 160.0, 3.42, "00:56:20.3", 1.2, "00:56:35.4", 7.8),
    ("Arroyo", -1.675721238961636, 38.0968434847573, 468.7, 13.0, "00:56:28.9", 5.0, "00:56:37.0", 5.0),
    ("La Sagra", -2.564, 37.981, 1530.0, 5.0, "00:56:36.5", 0.6, "00:56:45.6", 0.9),
    ("PixelSkies", -2.643889, 37.739722, 805.0, 7.2, "00:56:39.9", 1.7, "00:56:44.7", 1.5),
    ("La Palma", -17.8792, 28.762516, 2387.63, 1.2355, nothing, NaN, nothing, NaN),
    ("IOTA Scorpii", 9.856022392325654, 44.12703306956611, 52.0, 18.68, nothing, NaN, nothing, NaN),
    ("Belesta", 1.8175, 43.445408, 247.397, 3.0, nothing, NaN, nothing, NaN),
    ("Guirguillano", -1.865, 42.712053, 594.0, 10.008, nothing, NaN, nothing, NaN),
    ("Otivar", -3.6802975274933303, 36.81611111880951, 314.37, 5.45, nothing, NaN, nothing, NaN),
    ("Crni Vrh", 14.071284495949174, 45.94585244301794, 713.03, 5.9, nothing, NaN, nothing, NaN),
    ("Albox", -2.1518, 37.405564, 491.0, 8.03, nothing, NaN, nothing, NaN),
    ("Estelia", -5.4449518, 43.20358286, 630.0, 15.008, nothing, NaN, nothing, NaN),
    ("BOOTES-1", -6.734117424827392, 37.10408250638086, 61.0, 21.7, nothing, NaN, nothing, NaN),
    ("BOOTES-2", -4.04097, 36.759241, 70.0, 30.007, nothing, NaN, nothing, NaN),
    ("Calar Alto", -2.5461, 37.22361, 2168.0, 0.864, nothing, NaN, nothing, NaN),
    ("Sierra Nevada", -3.3847, 37.064136, 2930.527, 2.0, nothing, NaN, nothing, NaN),
    ("San Marcello Pistoiese", 10.8042, 44.063036, 965.411, 7.0, nothing, NaN, nothing, NaN),
    ("Asiago", 11.568824, 45.849444, 1370.0, 12.0, nothing, NaN, nothing, NaN),
]

"Gaia DR3 4360429542557512064 at the closest approach (Rommel et al. 2025, Table 1)."
huya2023_star() = Star("17 16 43.65331", "-07 00 20.3449"; name="Gaia DR3 4360429542557512064",
                       source_id="4360429542557512064", epoch=julian_year(Epoch("2023-06-24T00:58")),
                       mag=Dict("V" => 17.60), diameter=0.18 / 4.194e9 * 206264806.247)

"Ephemeris of Huya from the stored Horizons table."
function huya2023_ephemeris()
    rows = [split(l) for l in eachline(joinpath(HUYA2023_DIR, "huya_horizons_20230624.txt")) if !startswith(l, "#")]
    return EphemTable([r[1] for r in rows], parse.(Float64, getindex.(rows, 2)), parse.(Float64, getindex.(rows, 3)),
                      parse.(Float64, getindex.(rows, 4)); name="(38628) Huya")
end

"""
    huya2023_event(; negatives=:all, satellite=false) -> Event

The occultation of 2023 June 24 with the eleven positive Huya chords (Table 2) and the
negative chords (Table 8 and La Palma; `negatives = :all` or `:none`). With `satellite = true`
the La Palma chord is positive: the occultation by the satellite (Table 2), its contacts
assigned to the secondary body of a `BinaryLimb` fit.
"""
function huya2023_event(; negatives::Symbol=:all, satellite::Bool=false)
    negatives in (:all, :none) || throw(ArgumentError("negatives must be :all or :none"))
    body = Body("(38628) Huya", huya2023_ephemeris(); radius=205.0)
    ev = Event(huya2023_star(), body, "2023-06-24T00:58:00"; window=900, name="Huya, 2023-06-24")
    day = "2023-06-24T"
    for (name, lon, lat, alt, sens, ti, σi, te, σe) in HUYA2023_CHORDS
        obs = Observer(name, lon, lat, alt)
        if satellite && name == "La Palma"
            S = HUYA2023_PAPER.satellite
            push!(ev, Chord(name, obs; contacts=[Contact(:immersion, day * S.immersion[1], S.immersion[2]; body=:secondary),
                                                 Contact(:emersion, day * S.emersion[1], S.emersion[2]; body=:secondary)],
                            span=(day * "00:52:00", day * "01:02:00"), sensitivity=sens))
        elseif ti === nothing
            negatives === :all || continue
            push!(ev, Chord(name, obs; span=(day * "00:52:00", day * "01:02:00"), sensitivity=sens))
        else
            push!(ev, Chord(name, obs; immersion=day * ti, immersion_err=σi, emersion=day * te, emersion_err=σe,
                            sensitivity=sens))
        end
    end
    return ev
end
