# Stellar occultation by (174567) Varda on 2018 September 10
#
# Data from D. Souami et al. (2020), "A multi-chord stellar occultation by the large
# trans-Neptunian object (174567) Varda", A&A 643, A125 (doi:10.1051/0004-6361/202038526):
#   circumstances and star: Table 2; stations: Table 3; timings: Table 4; limb fit and
#   astrometry: Table 5; MacLaurin solutions: Table 6.
# Ephemeris: JPL Horizons (default solution), stored in varda_horizons_20180910.txt; the paper
# refers the fitted centre to NIMA v6.
#
# Notes on the transcription:
# * The two stations without data (IDY, CLP) are omitted. The observed intervals of the
#   negative chords are not tabulated; 03:28–03:42 UTC is used. The sensitivity of a chord is
#   its exposure time, equal to the cycle time (Table 3).
# * Table 5: p_V = 0.099 ± 0.002, although R′_equiv = 370 ± 7 km alone gives ± 0.004 (and the
#   footnote adds ± 0.03 for the rotational amplitude).
#
# `varda_event()` returns the `Event`; `VARDA_PAPER` the published results.

using Occultations

const VARDA_DIR = @__DIR__

"Published results of Souami et al. (2020)."
const VARDA_PAPER = (
    # Table 5 (P′: position angle of the apparent semi-minor axis, from north towards east)
    a=(383.0, 3.0), oblateness=(0.066, 0.047), R_equiv=(370.0, 7.0), position_angle=(67.0, 8.0), chi2_pdf=1.53,
    npoints=10, H_V=(3.81, 0.01), albedo=(0.099, 0.002),
    astrometry=(epoch="2018-09-10T03:34:00", ra="17 18 25.12972", dec="-02 05 14.3055", ra_err=0.197, dec_err=0.143),
    # Table 4: durations [s] and chord lengths [km]
    durations=Dict("TWF" => (61.73, 0.85), "MHV" => (68.28, 0.82), "YMA" => (69.29, 0.84), "CAR" => (67.13, 0.66),
                   "FLO" => (63.45, 0.20)),
    chords=Dict("TWF" => (670.4, 9.2), "MHV" => (741.5, 8.9), "YMA" => (752.5, 9.1), "CAR" => (729.0, 7.2),
                "FLO" => (689.1, 2.2)),
    # Table 2
    distance_au=46.6335, velocity=10.86, star_diameter_km=0.82,
)

# station, longitude (E), latitude, altitude [m], exposure = cycle time [s], ingress, σ, egress, σ
const VARDA_CHORDS = [
    ("TWF", "-114 28 13.1", "42 35 01.9", 1133.0, 4.271, "03:32:57.10", 0.60, "03:33:58.83", 0.60),
    ("MHV", "-114 35 48.9", "35 01 54.1", 184.0, 2.135, "03:33:56.99", 0.76, "03:35:05.27", 0.30),
    ("YMA", "-114 26 10.5", "32 39 34.0", 97.0, 2.135, "03:34:17.10", 0.60, "03:35:26.39", 0.59),
    ("CAR", "-111 57 07.9", "33 48 42.9", 654.0, 1.068, "03:34:12.39", 0.47, "03:35:19.52", 0.47),
    ("FLO", "-111 21 00.6", "33 00 54.3", 484.0, 1.068, "03:34:23.30", 0.13, "03:35:26.75", 0.15),
    ("CPSLO", "-120 39 36.0", "35 18 01.8", 119.0, 2.135, nothing, NaN, nothing, NaN),
    ("GAR", "-119 40 20.3", "38 53 23.5", 1502.0, 1.068, nothing, NaN, nothing, NaN),
    ("UMA", "-119 17 53.9", "45 55 19.9", 140.0, 1.068, nothing, NaN, nothing, NaN),
    ("BDR", "-105 09 46.8", "40 15 09.0", 1587.0, 0.5, nothing, NaN, nothing, NaN),
    ("DOD", "-90 08 31.1", "42 57 36.9", 390.0, 1.068, nothing, NaN, nothing, NaN),
    ("NAP", "-88 06 59.7", "41 45 32.5", 230.0, 1.0, nothing, NaN, nothing, NaN),
]

"Gaia DR2 4367203805493754752 at the occultation epoch (Souami et al. 2020, Table 2)."
varda_star() = Star("17 18 25.12481", "-02 05 14.4282"; name="Gaia DR2 4367203805493754752",
                    source_id="4367203805493754752", epoch=julian_year(Epoch("2018-09-10T03:34")),
                    ra_err=0.2, dec_err=0.1, mag=Dict("G" => 14.7, "RP" => 13.8),
                    diameter=VARDA_PAPER.star_diameter_km / (VARDA_PAPER.distance_au * 1.495978707e8) * 206264806.247)

"Ephemeris of Varda from the stored Horizons table."
function varda_ephemeris()
    rows = [split(l) for l in eachline(joinpath(VARDA_DIR, "varda_horizons_20180910.txt")) if !startswith(l, "#")]
    return EphemTable([r[1] for r in rows], parse.(Float64, getindex.(rows, 2)), parse.(Float64, getindex.(rows, 3)),
                      parse.(Float64, getindex.(rows, 4)); name="(174567) Varda")
end

"""
    varda_event(; negatives=:all) -> Event

The occultation of 2018 September 10: five positive chords (Table 4) and the negative chords
of Table 3 (`negatives = :all` or `:none`).
"""
function varda_event(; negatives::Symbol=:all)
    negatives in (:all, :none) || throw(ArgumentError("negatives must be :all or :none"))
    body = Body("(174567) Varda", varda_ephemeris(); radius=370.0)
    ev = Event(varda_star(), body, "2018-09-10T03:34:00"; window=1200, name="Varda, 2018-09-10")
    day = "2018-09-10T"
    for (name, lon, lat, alt, sens, ti, σi, te, σe) in VARDA_CHORDS
        obs = Observer(name, lon, lat, alt)
        if ti === nothing
            negatives === :all || continue
            push!(ev, Chord(name, obs; span=(day * "03:28:00", day * "03:42:00"), sensitivity=sens))
        else
            push!(ev, Chord(name, obs; immersion=day * ti, immersion_err=σi, emersion=day * te, emersion_err=σe,
                            sensitivity=sens))
        end
    end
    return ev
end
