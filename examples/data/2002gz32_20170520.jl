# Stellar occultation by the centaur (95626) 2002 GZ32 on 2017 May 20
#
# Data from P. Santos-Sanz et al. (2021), "The 2017 May 20 stellar occultation by the elongated
# centaur (95626) 2002 GZ32", MNRAS 501, 6062–6075 (doi:10.1093/mnras/staa3881):
#   stations: Table 3; star: Table 4; timings and the shifts of the aligned chords: Table 5;
#   limb fits: Table 6; albedo, rotation and 3-D shape: Sect. 4.
# Ephemeris: JPL Horizons (default solution), stored in 2002gz32_horizons_20170520.txt; the paper
# refers the fitted centre to JPL#23.
#
# Notes on the transcription:
# * The star: Gaia DR2 4332853000658702848 at the event (Table 4, coordinates b).
# * La Hita observed with two telescopes; the 40 cm video chord has an unknown time offset
#   (lost frames, computer problem) and was not used by the authors: kept here but disabled.
# * Observed intervals of the negative chords are not tabulated; 01:20–01:40 UTC is used.
# * The sensitivity of a chord is its cycle time (Table 3).
#
# `gz32_event()` returns the `Event`; `GZ32_PAPER` the published results.

using Occultations

const GZ32_DIR = @__DIR__

"Published results of Santos-Sanz et al. (2021)."
const GZ32_PAPER = (
    # Table 6 (centre relative to JPL#23; P′ = position angle of the minor axis, N → E)
    original=(a=(146.0, 12.0), b=(73.0, 4.6), f=(-500.0, 9.0), g=(-112.0, 6.0), position_angle=(318.0, 3.0), chi2=29.5,
              D_eq=(206.0, 15.0), albedo=(0.044, 0.009)),
    shifted=(a=(152.6, 8.4), b=(73.0, 4.2), f=(-501.0, 6.0), g=(-107.0, 5.0), position_angle=(322.0, 3.0), chi2=4.0,
             D_eq=(211.0, 12.0), albedo=(0.043, 0.007)),
    # Table 5: time shifts of the aligned chords [s]
    shifts=Dict("Allariz" => 0.013, "Javalambre OAJ" => 0.006, "Aras de los Olmos" => -0.771, "La Hita 77cm" => -0.258,
                "Univ. of Athens" => 2.048),
    # absolute magnitude, rotation (light-curve amplitude 0.13 mag, period 5.80 h, the object close
    # to its brightness minimum at the occultation), 3-D shape (Sect. 4.2)
    H_V=(7.39, 0.06), amplitude=0.13, period_h=5.80, axes_3d=(366.0, 306.0, 120.0),
    velocity=21.19, distance_au=17.1033, star_diameter_mas=0.0472, fresnel_km=0.88,
)

# station, longitude (E), latitude, altitude [m], cycle time [s],
# immersion, σ, emersion, σ (UTC, 2017-05-20; `nothing` for a negative chord)
const GZ32_CHORDS = [
    ("Allariz", "-07 46 13.0", "42 11 57.0", 514.0, 2.48510, "01:30:51.50", 0.50, "01:30:55.75", 0.15),
    ("Javalambre OAJ", "-01 00 58.6", "40 02 30.6", 1957.0, 2.68829, "01:30:19.00", 0.10, "01:30:27.15", 0.08),
    ("Aras de los Olmos", "-01 06 05.4", "39 56 42.0", 1280.0, 1.49875, "01:30:19.50", 1.70, "01:30:28.25", 0.14),
    ("La Hita 77cm", "-03 11 09.8", "39 34 07.0", 674.0, 3.61480, "01:30:26.16", 0.23, "01:30:34.19", 1.70),
    ("La Hita 40cm", "-03 11 09.8", "39 34 07.0", 674.0, 0.301, "01:30:25.86", 0.90, "01:30:34.32", 0.69),
    ("Univ. of Athens", "23 47 00.1", "37 58 06.8", 250.0, 3.17352, "01:28:37.02", 2.66, "01:28:43.31", 2.26),
    ("Sant Esteve", "01 52 21.1", "41 29 37.5", 180.0, 0.16, nothing, NaN, nothing, NaN),
    ("La Sagra", "-02 33 55.0", "37 58 56.6", 1530.0, 0.301, nothing, NaN, nothing, NaN),
]

"Gaia DR2 4332853000658702848 at the occultation epoch (Santos-Sanz et al. 2021, Table 4)."
gz32_star() = Star("16 50 08.5378", "-13 05 12.7381"; name="Gaia DR2 4332853000658702848",
                   source_id="4332853000658702848", epoch=julian_year(Epoch("2017-05-20T01:30")),
                   mag=Dict("B" => 14.900, "V" => 14.020, "R" => 13.810, "K" => 10.290),
                   diameter=GZ32_PAPER.star_diameter_mas)

"Ephemeris of 2002 GZ32 from the stored Horizons table."
function gz32_ephemeris()
    rows = [split(l) for l in eachline(joinpath(GZ32_DIR, "2002gz32_horizons_20170520.txt")) if !startswith(l, "#")]
    return EphemTable([r[1] for r in rows], parse.(Float64, getindex.(rows, 2)), parse.(Float64, getindex.(rows, 3)),
                      parse.(Float64, getindex.(rows, 4)); name="(95626) 2002 GZ32")
end

"""
    gz32_event(; negatives=:all) -> Event

The occultation of 2017 May 20 with the six positive chords (Table 5; the La Hita 40 cm chord
disabled, as in the paper) and the two closest negative chords (`negatives = :all` or `:none`).
"""
function gz32_event(; negatives::Symbol=:all)
    negatives in (:all, :none) || throw(ArgumentError("negatives must be :all or :none"))
    body = Body("(95626) 2002 GZ32", gz32_ephemeris(); radius=105.0)
    ev = Event(gz32_star(), body, "2017-05-20T01:30:00"; window=1200, name="2002 GZ32, 2017-05-20")
    day = "2017-05-20T"
    for (name, lon, lat, alt, sens, ti, σi, te, σe) in GZ32_CHORDS
        obs = Observer(name, lon, lat, alt)
        if ti === nothing
            negatives === :all || continue
            push!(ev, Chord(name, obs; span=(day * "01:20:00", day * "01:40:00"), sensitivity=sens))
        else
            c = Chord(name, obs; immersion=day * ti, immersion_err=σi, emersion=day * te, emersion_err=σe, sensitivity=sens)
            push!(ev, name == "La Hita 40cm" ? disable(c) : c)
        end
    end
    return ev
end
