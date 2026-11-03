# Stellar occultation by (84922) 2003 VS2 on 2019 October 22
#
# Data from M. Vara-Lubiano et al. (2022), "The multichord stellar occultation on 2019 October 22
# by the trans-Neptunian object (84922) 2003 VS2", A&A 663, A121
# (doi:10.1051/0004-6361/202141842):
#   star: Table 2; stations: Tables 3 (positive) and 4 (negative); timings and the shifts of the
#   aligned chords: Table 5; limb fits: Table 6; rotation: Sect. 3.3; 3-D shape: Table 7.
# Ephemeris: JPL Horizons (default solution), stored in 2003vs2_horizons_20191022.txt; the paper
# refers the fitted centre to JPL#33.
#
# Notes on the transcription:
# * Chords 3 and 4 (two telescopes at one site, exposures interleaved) were merged by the
#   authors into chord 3′; its timings are used with the coordinates of chord 3.
# * Chord 8: the acquisition clock was off by more than a minute (Sect. 3.1); only its length
#   is meaningful. The paper shifts it by −14 s onto the line of the other chord centres.
# * Chord 12: times only to the second (Sect. 3.1), ± 2 s added by the authors.
# * The Romanian sites without a name are labelled by their chord number.
# * The negative chord of H. Mikuž (Vratnik, Croatia), the closest one in the north, has in
#   Table 4 the same coordinates as PDlink Observatory (Cadca, Slovakia) – a copy error – and is
#   omitted: its true position is not known.
# * Observed intervals of the negative chords are not tabulated ("several minutes before and
#   after"); 20:30–20:55 UTC is used. Their sensitivity is not given; 1 s is assumed.
# * Cycle times (exposure + dead time, Table 3) are the sensitivities of the positive chords.
#
# `vs2_event()` returns the `Event`; `VS2_PAPER` the published results.

using Occultations

const VS2_DIR = @__DIR__

"Published results of Vara-Lubiano et al. (2022)."
const VS2_PAPER = (
    # Table 6: limb fits; the tilt angle (PA of the minor axis) is counted from north towards
    # the WEST in the paper: −11° there is +11° (north through east) here
    original=(a=(292.0, 3.0), b=(231.0, 6.0), oblateness=(0.21, 0.03), position_angle=(11.0, 2.0),
              f=(-1010.0, 2.0), g=(770.0, 5.0), D_eq=(519.0, 10.0)),
    aligned=(a=(293.0, 3.0), b=(230.0, 6.0), oblateness=(0.21, 0.03), position_angle=(6.0, 2.0),
             f=(-1007.0, 2.0), g=(769.0, 5.0), D_eq=(519.0, 10.0)),
    # Table 5: shifts of the aligned chords [s] (positive: towards the east)
    shifts=Dict("ROASTERR-1 (1)" => 0.7, "ROASTERR-1 (2)" => 0.6, "Romania 3'" => 0.03, "Romania 5" => -0.06,
                "Berthelot" => 1.6, "Stardust" => 0.07, "Romania 8" => -14.0, "Stardreams" => 1.1,
                "Vidojevica" => 0.3, "Belogradchik" => -0.7, "Romania 12" => -2.0),
    D_mean=(545.0, 13.0), albedo=(0.134, 0.010), H_V=(4.14, 0.07),
    # rotation (Sect. 3.3): P = 7.41753 h, Δm = 0.264 ± 0.017 mag, phase 0.32 after the maximum
    period_h=7.41753, amplitude=(0.264, 0.017), rotational_phase=0.32,
    # Table 7 (observed Δm, original chords): a, b, c [km], aspect [deg]
    axes_3d=((339.0, 5.0), (235.0, 6.0), (226.0, 8.0)), aspect=(59.0, 2.0), D_eq_volume=(524.0, 7.0),
    pole_ecliptic=(228.0, 39.0),
    velocity=13.01, fresnel_km=1.33, star_diameter_km=(0.86, 0.91),
)

# station, longitude (E), latitude, altitude [m], cycle time [s],
# immersion, σ, emersion, σ (UTC, 2019-10-22; `nothing` for a negative chord)
const VS2_CHORDS = [
    ("ROASTERR-1 (1)", "23 35 47", "46 49 15.6", 391.0, 5.18, "20:41:34.800", 1.155, "20:42:14.995", 0.225),
    ("ROASTERR-1 (2)", "23 35 47", "46 49 15.6", 391.0, 3.01, "20:41:35.54", 0.46, "20:42:14.47", 0.55),
    ("Romania 3'", "23 35 35.647", "46 42 37.557", 783.405, 3.0, "20:41:35.609", 0.109, "20:42:16.585", 0.077),
    ("Romania 5", "21 26 12.78", "45 42 11.94", 92.0, 5.14, "20:41:44.375", 0.525, "20:42:24.95", 0.40),
    ("Berthelot", "22 53 19.7", "45 36 59.4", 400.0, 8.64, "20:41:37.325", 1.100, "20:42:24.062", 1.287),
    ("Stardust", "25 37 19", "45 38 30", 597.0, 10.28, "20:41:35.420", 0.425, "20:42:17.411", 0.382),
    ("Romania 8", "23 47 19.0", "44 19 24.5", 109.0, 3.34, "20:42:02.05", 0.73, "20:42:37.76", 0.76),
    ("Stardreams", "26 02 44.2", "45 12 13.3", 382.0, 5.52, "20:41:37.491", 0.590, "20:42:14.897", 0.160),
    ("Vidojevica", "21 33 20.4", "43 08 24.6", 1150.0, 1.55, "20:41:57.914", 0.110, "20:42:31.944", 0.136),
    ("Belogradchik", "22 40 30", "43 37 22", 650.0, 6.86, "20:41:54.60", 0.15, "20:42:28.76", 0.12),
    ("Romania 12", "25 59 8.69", "44 19 19.14", 70.0, 2.50, "20:41:53.551", 2.186, "20:42:12.949", 2.387),
    ("Bremthal", "08 21 50.4", "50 08 17.4", 256.0, 1.0, nothing, NaN, nothing, NaN),
    ("PDlink, Cadca", "18 42 09.47", "49 24 15.20", 680.0, 1.0, nothing, NaN, nothing, NaN),
    ("Mt. Suhora", "20 04 03", "49 34 09", 1009.0, 1.0, nothing, NaN, nothing, NaN),
    ("Skalnate Pleso", "20 14 02.1", "49 11 21.8", 1786.0, 1.0, nothing, NaN, nothing, NaN),
    ("Massa", "10 07 56.65", "44 01 16.96", 30.0, 1.0, nothing, NaN, nothing, NaN),
    ("Gothard, Szombathely", "16 36 15.67", "47 15 29.83", 220.0, 1.0, nothing, NaN, nothing, NaN),
    ("Crni Vrh", "14 04 16.4", "45 56 45.1", 726.0, 1.0, nothing, NaN, nothing, NaN),
    ("Monte Agliale", "10 30 53.7", "43 59 43.0", 756.0, 1.0, nothing, NaN, nothing, NaN),
    ("Basovizza", "13 52 23", "45 38 33", 400.0, 1.0, nothing, NaN, nothing, NaN),
    ("Piszkesteto", "19 53 41.7", "47 55 6.0", 924.0, 1.0, nothing, NaN, nothing, NaN),
    ("Konkoly, Normafa", "18 57 51.12", "47 29 59.28", 468.0, 1.0, nothing, NaN, nothing, NaN),
    ("Romania (Dumitrescu)", "26 13 04.3", "44 07 56.8", 88.0, 1.0, nothing, NaN, nothing, NaN),
    ("Mount Etna", "14 58 04", "37 41 05", 1735.0, 1.0, nothing, NaN, nothing, NaN),
]

"Gaia DR2 3449076721168026496 at the occultation epoch (Vara-Lubiano et al. 2022, Table 2)."
vs2_star() = Star("05 30 38.0442", "+33 07 01.748"; name="Gaia DR2 3449076721168026496",
                  source_id="3449076721168026496", epoch=julian_year(Epoch("2019-10-22T20:36:31")),
                  pmra=-0.34, pmdec=-1.76, parallax=0.2084, ra_err=0.0401, dec_err=0.0330,
                  mag=Dict("G" => 14.1625, "B" => 15.030, "V" => 14.150, "K" => 10.968),
                  diameter=0.88 / 5.40e9 * 206264806.247)

"Ephemeris of 2003 VS2 from the stored Horizons table."
function vs2_ephemeris()
    rows = [split(l) for l in eachline(joinpath(VS2_DIR, "2003vs2_horizons_20191022.txt")) if !startswith(l, "#")]
    return EphemTable([r[1] for r in rows], parse.(Float64, getindex.(rows, 2)), parse.(Float64, getindex.(rows, 3)),
                      parse.(Float64, getindex.(rows, 4)); name="(84922) 2003 VS2")
end

"""
    vs2_event(; negatives=:all) -> Event

The occultation of 2019 October 22 with the eleven positive chords (Table 5, chords 3 and 4
merged) and the negative chords of Table 4 (`negatives = :all` or `:none`).
"""
function vs2_event(; negatives::Symbol=:all)
    negatives in (:all, :none) || throw(ArgumentError("negatives must be :all or :none"))
    body = Body("(84922) 2003 VS2", vs2_ephemeris(); radius=260.0)
    ev = Event(vs2_star(), body, "2019-10-22T20:42:00"; window=1200, name="2003 VS2, 2019-10-22")
    day = "2019-10-22T"
    for (name, lon, lat, alt, sens, ti, σi, te, σe) in VS2_CHORDS
        obs = Observer(name, lon, lat, alt)
        if ti === nothing
            negatives === :all || continue
            push!(ev, Chord(name, obs; span=(day * "20:30:00", day * "20:55:00"), sensitivity=sens))
        else
            push!(ev, Chord(name, obs; immersion=day * ti, immersion_err=σi, emersion=day * te, emersion_err=σe,
                            sensitivity=sens))
        end
    end
    return ev
end
