# Stellar occultation by (143707) 2003 UY117 on 2020 October 23
#
# Data from M. Kretlow et al. (2024), "Physical properties of trans-Neptunian object (143707)
# 2003 UY117 derived from stellar occultation and photometric observations", A&A 691, A31
# (doi:10.1051/0004-6361/202451329):
#   star: Table 1 (Gaia DR3 position at the occultation epoch); stations: Table A.1;
#   timings: Table 2; limb fit: Table 3; rotational light curve and 3-D shape: Sects. 3–4.
# Ephemeris: JPL Horizons (default solution), stored in uy117_horizons_20201023.txt; the paper
# refers the fitted centre to JPL#20.
#
# Notes on the transcription:
# * The observed intervals of the negative chords are not tabulated; 22:16–22:28 UTC, which
#   brackets the event (≈ 22:21:40–22:22:00 at the positive sites), is used.
# * The sensitivity of a negative chord (for the soft constraint) is its exposure plus dead
#   time (Table A.1).
#
# `uy117_event()` returns the `Event`; `UY117_PAPER` the published results.

using Occultations

const UY117_DIR = @__DIR__

"Published results of Kretlow et al. (2024)."
const UY117_PAPER = (
    # Table 3: limb fit (centre relative to JPL#20); φ′ from north, clockwise positive in the
    # (f, g) plane with f to the right, i.e. towards +f (east) as the position angles here
    center_f=(119.0, 7.0), center_g=(-10.0, 8.0), a=(141.0, 9.0), b=(94.0, 16.0), position_angle=(-31.0, 11.0),
    oblateness=(0.33, 0.11), D_eq=(228.0, 21.0), chi2_pdf=0.43,
    # rotational light curve (Sect. 3.2) and 3-D solution (Sect. 4.1); aspect from the pole
    period_h=12.376, amplitude=(0.36, 0.13), rotational_phase=0.63,
    axes_3d=((166.0, 12.0), (108.0, 12.0), (90.0, 13.0)), aspect=70.0, D_eq_volume=(235.0, 25.0),
    pole=((337.0, 10.0), (62.0, 5.0)),
    # albedo and astrometry (Sects. 4.2, 4.3)
    H_V=(5.97, 0.07), albedo=(0.139, 0.027),
    astrometry=(epoch="2020-10-23T22:19:21.8", ra="03 23 26.21404", dec="+22 47 19.3184", ra_err=0.3, dec_err=0.4),
    distance_au=33.47, velocity=21.80, star_diameter_km=0.4, fresnel_km_700nm=1.32,
)

# station, longitude (E), latitude, altitude [m], exposure + dead time [s],
# immersion, σ, emersion, σ (UTC, 2020-10-23; `nothing` for a negative chord)
const UY117_CHORDS = [
    ("Arenosillo", "-06 44 06.0", "37 05 53.0", 40, 7.567, "22:21:49.83", 0.69, "22:22:00.44", 0.54),
    ("Granada", "-03 38 22.1", "37 06 41.0", 730, 2.061, "22:21:42.36", 0.14, "22:21:50.48", 0.08),
    ("OSN90", "-03 23 09.35", "37 03 46.40", 2841, 5.181, "22:21:41.12", 1.25, "22:21:48.45", 0.23),
    ("OSN150", "-03 23 09.67", "37 03 46.52", 2841, 4.848, "22:21:42.18", 1.03, "22:21:49.13", 0.26),
    ("Leeds", "-01 36 28.0", "53 50 15.4", 114, 0.64, nothing, NaN, nothing, NaN),
    ("Stevenage", "-00 03 51.7", "51 57 04.4", 120, 0.34, nothing, NaN, nothing, NaN),
    ("Oxford", "-01 18 46.3", "51 55 41.2", 199, 1.0, nothing, NaN, nothing, NaN),
    ("Abingdon", "-01 16 55.2", "51 37 53.1", 59, 0.64, nothing, NaN, nothing, NaN),
    ("Dourbes", "04 34 56.0", "50 05 25.9", 195, 0.32, nothing, NaN, nothing, NaN),
    ("Cannet", "-00 02 40.7", "43 37 15.2", 180, 1.28, nothing, NaN, nothing, NaN),
    ("IST40 Istanbul", "28 57 56.34", "41 00 42.48", 60, 0.9, nothing, NaN, nothing, NaN),
    ("Ibiza", "01 14 26.9", "38 53 28.0", 166, 6.0, nothing, NaN, nothing, NaN),
    ("Albox", "-02 09 06.5", "37 24 20.0", 493, 3.3, nothing, NaN, nothing, NaN),
    ("CAHA220", "-02 32 46.30", "37 13 23.30", 2168, 9.0, nothing, NaN, nothing, NaN),
    ("Las Negras", "-02 00 54.0", "36 52 50.0", 50, 12.5, nothing, NaN, nothing, NaN),
    ("Teide MUSCAT2", "-16 18 29.0", "28 28 28.6", 2380, 23.0, nothing, NaN, nothing, NaN),
    ("Teide TAR", "-16 18 29.0", "28 28 28.6", 2380, 0.5, nothing, NaN, nothing, NaN),
]

"Gaia DR3 62553763421986304 at the occultation epoch (Kretlow et al. 2024, Table 1)."
uy117_star() = Star("03 23 26.2155", "+22 47 19.194"; name="Gaia DR3 62553763421986304",
                    source_id="62553763421986304", epoch=julian_year(Epoch("2020-10-23T22:20")),
                    ra_err=0.24, dec_err=0.14, mag=Dict("G" => 14.54, "V" => 14.54, "K" => 12.22), diameter=0.0178)

"Ephemeris of 2003 UY117 from the stored Horizons table."
function uy117_ephemeris()
    rows = [split(l) for l in eachline(joinpath(UY117_DIR, "2003uy117_horizons_20201023.txt")) if !startswith(l, "#")]
    return EphemTable([r[1] for r in rows], parse.(Float64, getindex.(rows, 2)), parse.(Float64, getindex.(rows, 3)),
                      parse.(Float64, getindex.(rows, 4)); name="(143707) 2003 UY117")
end

"""
    uy117_event(; negatives=:all) -> Event

The occultation of 2020 October 23 with the four positive chords (Table 2) and the negative
chords of Table A.1 (`negatives = :all`, `:nearest` for CAHA220 and Cannet as in the
paper's Fig. 3, or `:none`). The body carries the rotation period of the light curve.
"""
function uy117_event(; negatives::Symbol=:all)
    negatives in (:all, :nearest, :none) || throw(ArgumentError("negatives must be :all, :nearest or :none"))
    body = Body("(143707) 2003 UY117", uy117_ephemeris(); radius=115.0)
    ev = Event(uy117_star(), body, "2020-10-23T22:20:00"; window=900, name="2003 UY117, 2020-10-23")
    day = "2020-10-23T"
    for (name, lon, lat, alt, sens, ti, σi, te, σe) in UY117_CHORDS
        obs = Observer(name, lon, lat, alt)
        if ti === nothing
            (negatives === :all || (negatives === :nearest && name in ("CAHA220", "Cannet"))) || continue
            push!(ev, Chord(name, obs; span=(day * "22:16:00", day * "22:28:00"), sensitivity=sens))
        else
            push!(ev, Chord(name, obs; immersion=day * ti, immersion_err=σi, emersion=day * te, emersion_err=σe,
                            sensitivity=sens))
        end
    end
    return ev
end
