# Stellar occultations by the centaur (54598) Bienor, 2019–2023
#
# Data from
# * E. Fernández-Valenzuela et al. (2023), "The multichord stellar occultation by the centaur
#   Bienor on January 11, 2019", A&A 669, A112 (doi:10.1051/0004-6361/202243214):
#   star: Table 1; stations: Table 2; timings and the chord shifts: Table 3; ellipses
#   (scenarios S1, S2 and the binary S3): Table 4; rotation: Sect. 5;
# * J. L. Rizos et al. (2024), "A study of centaur (54598) Bienor from multiple stellar
#   occultations and rotational light curves", A&A 689, A82 (doi:10.1051/0004-6361/202450833):
#   stars: Table 1; ellipse of Occ. B: Table 2; aspect and pole position angle: Table 3;
#   rotational phases: Table 5; stations: Table B.1; timings: Table B.2; 3-D shape: Sect. 3.4.
# Ephemeris: JPL Horizons (default solution), stored in bienor_horizons_<yyyymmdd>.txt.
#
# Notes on the transcription:
# * 2019: the star is Gaia DR2 236278279941456256, used with its proper motion (the paper gives
#   the DR2 position at the event). Table 2 lists longitudes as "W"; Sabadell (357° 54′ 57.4″)
#   and the TJO (0° 43′ 46.74″) are east longitudes of the real sites, as used below.
#   Observed intervals of the negative chords are not tabulated: all recordings started about
#   15 min before and ended 15 min after the predicted time; 00:50–01:17 UTC is used.
# * 2022–2023: the stars are Gaia DR3 (Table 1), used with their proper motions (retrieved
#   2026-11-02). Rizos et al. (2024) give no negative chords.
# * Occ. C took place on 2023 February 14 (Tables B.1, B.2; "2024" in Table 5 is a misprint).
# * The rotation pole: Table 3 gives the aspect angle ψ (from the direction to the observer) and
#   the position angle θp of the pole for each event; ψ and θp of Occ. B give the prograde pole
#   (α, δ) = (3.58°, +58.23°), which reproduces those of Occ. A and C within 0.1°.
# * The sensitivity of a chord is its cycle time.
#
# `bienor_event(name)` returns the `Event` of `"2019"`, `"A"`, `"B"` or `"C"`; `BIENOR_PAPER`
# the published results.

using Occultations

const BIENOR_DIR = @__DIR__

"Published results of Fernández-Valenzuela et al. (2023) and Rizos et al. (2024)."
const BIENOR_PAPER = (
    # 2019 (FV23, Table 4): ap, bp [km], tilt θ [deg] (PA of the minor axis = 180° − θ: S1 52°,
    # S2 61°), centre (x, y) [km], instantaneous and rotationally averaged area-equivalent
    # diameters, albedo [%]; S3: two ellipses E1 (east) and E2 (west)
    S1=(a=(94.0, 27.0), b=(66.0, 5.0), theta=(128.0, 13.0), D_eq=(160.0, 20.0)),
    S2=(a=(97.0, 24.0), b=(64.0, 4.0), theta=(119.0, 7.0), D_eq=(160.0, 20.0)),
    S3=(E1=(a=(84.0, 6.0), b=(57.0, 4.0), theta=(8.0, 10.0), x=(158.0, 8.0), y=(112.0, 3.0)),
        E2=(a=(67.0, 8.0), b=(60.0, 6.0), theta=(7.0, 10.0), x=(83.0, 10.0), y=(210.0, 9.0)),
        D_eq=(188.0, 9.0), D_eq_mean=(183.0, 8.0)),
    D_eq_mean_2019=(150.0, 20.0), albedo_2019=(0.08, 0.02),
    # Table 3: chord shifts [s] of S1 and S2 (positive: towards the south-west)
    shifts_S1=Dict("Manzanares" => -4.06),
    shifts_S2=Dict("Javalambre" => 2.6, "La Hita T77" => 1.01, "La Hita T40" => 1.01, "Manzanares" => -0.5,
                   "Constancia" => 3.0),
    velocity_2019=13.47, distance_2019_au=14.16, amplitude_2019=(0.172, 0.003), amplitude_at_event_2019=0.058,
    H_V=(7.47, 0.04), period_2019_h=9.1719,
    # 2022 Occ. B (Rizos et al. 2024, Table 2): semi-axes u, v, position angle, centre
    B=(u=(107.0, 3.7), v=(51.9, 2.1), position_angle=(19.3, 3.4), center=((15.5, 2.7), (12.9, 3.8)),
       D_eq=(149.0, 4.0), albedo=(0.065, 0.005), H_V_inst=(7.72, 0.04)),
    # aspect angle ψ and pole position angle θp (prograde solution), rotational phases
    # (0 = maximum of light), sidereal period [h]
    aspect=Dict("A" => 133.63, "B" => 127.14, "C" => 127.98), pole_pa=Dict("A" => 315.58, "B" => 318.69, "C" => 318.18),
    pole=(3.58, 58.23), phases=Dict("A" => (0.888, 0.006), "B" => (0.15, 0.02), "C" => (0.93, 0.02)),
    period_h=(9.1736, 0.0002), amplitude_2021_2023=(0.23, 0.30),
    # 3-D shape (Sect. 3.4); Van Vleck shifted by +4.0 s in Fig. C.1
    axes_3d=((127.0, 5.0), (55.0, 4.0), (45.0, 4.0)), van_vleck_shift=4.0,
)

# station, longitude (E), latitude, altitude [m], cycle time [s],
# immersion, σ, emersion, σ (UTC; `nothing` for a negative chord)
const BIENOR_CHORDS = Dict(
    "2019" => [
        ("Javalambre", "-01 00 58.58", "40 02 30.58", 1957.0, 6.928, "01:03:27.8", 1.5, "01:03:34.7", 1.5),
        ("La Hita T77", "-03 11 09.7", "39 34 06.8", 695.0, 2.957, "01:03:33.2", 0.7, "01:03:43.2", 0.2),
        ("La Hita T40", "-03 11 09.7", "39 34 06.8", 695.0, 4.0, "01:03:33.5", 0.5, "01:03:42.8", 0.6),
        ("Manzanares", "-03 27 23.3", "39 04 55.2", 132.0, 1.495, "01:03:33.9", 0.1, "01:03:41.7", 0.2),
        ("Constancia", "-08 19 25.2", "39 29 41.6", 147.0, 2.56, "01:03:58.2", 0.7, "01:04:04.9", 0.5),
        ("Sabadell", "02 05 02.6", "41 33 00.2", 224.0, 2.56, nothing, NaN, nothing, NaN),
        ("La Sagra", "-02 33 52", "37 59 02", 1530.0, 3.0, nothing, NaN, nothing, NaN),
        ("Calar Alto", "-02 32 46", "37 13 25", 2168.0, 1.246, nothing, NaN, nothing, NaN),
        ("TJO", "00 43 46.74", "42 03 05.95", 1570.0, 6.0, nothing, NaN, nothing, NaN),
    ],
    "A" => [
        ("Van Vleck", "-72 39 33.12", "41 33 18.6192", 68.9, 2.32, "08:44:44.48", 1.56, "08:44:51.67", 1.56),
        ("Westport", "-73 19 39.2999", "41 10 15.8999", 88.0, 0.25, "08:44:46.7731", 0.0461, "08:44:50.7746", 0.0461),
    ],
    "B" => [
        ("Salvia", "-00 24 26.3196", "47 58 57.5184", 99.0, 0.5, "20:16:18.0034", 0.1078, "20:16:21.0038", 0.1078),
        ("Schreurs", "05 33 27.4386", "50 38 35.0919", 107.0, 0.32, "20:15:56.6541", 0.0883, "20:16:00.8141", 0.0883),
        ("Dourbes", "04 34 55.992", "50 05 25.872", 192.0, 0.16, "20:15:59.9112", 0.0174, "20:16:04.3342", 0.0174),
        ("Tsu Mie", "136 27 28.5012", "34 39 58.1004", 5.4, 0.49, "20:09:39.6022", 0.1217, "20:09:43.5234", 0.1217),
        ("Pressigny", "00 47 17.88", "46 55 30.606", 58.0, 0.3, "20:16:16.3726", 0.0356, "20:16:19.9788", 0.0356),
        ("Uda Nara", "136 00 40.7016", "34 34 23.9016", 407.0, 0.338, "20:09:39.7147", 0.2112, "20:09:42.4180", 0.2112),
        ("Owase Mie", "136 11 34.6992", "34 04 20.6264", 37.0, 0.499, "20:09:37.7470", 0.4980, "20:09:38.7440", 0.4980),
    ],
    "C" => [
        ("Cottered", "-00 03 51.6996", "51 57 04.3992", 120.0, 1.603, "01:37:18.7117", 0.8173, "01:37:26.7260", 0.8173),
        ("Nandrin", "05 26 29.3604", "50 31 24.9312", 261.0, 2.56, "01:36:59.9112", 0.5737, "01:37:07.5912", 0.5737),
    ],
)

# date, reference time, ephemeris file of each event
const BIENOR_EVENTS = Dict(
    "2019" => ("2019-01-11", "01:04:00", "bienor_horizons_20190111.txt"),
    "A" => ("2022-02-06", "08:45:00", "bienor_horizons_20220206.txt"),
    "B" => ("2022-12-26", "20:13:00", "bienor_horizons_20221226.txt"),
    "C" => ("2023-02-14", "01:37:00", "bienor_horizons_20230214.txt"),
)

"""
    bienor_star(name)

The occulted star of event `name` (`"2019"`, `"A"`, `"B"`, `"C"`): Gaia DR2 (2019) or Gaia DR3
with proper motion.
"""
function bienor_star(name::AbstractString)
    name == "2019" &&
        return Star(49.990729644389255, 40.68145402669284; name="Gaia DR2 236278279941456256",
                    source_id="236278279941456256", epoch=2015.5, pmra=1.4159659135641327, pmdec=-3.6929475155421807,
                    parallax=0.47373831615420353, ra_err=0.0442, dec_err=0.0343, pmra_err=0.1075, pmdec_err=0.0927,
                    parallax_err=0.0558, mag=Dict("G" => 16.002, "V" => 16.329, "K" => 14.417))
    name == "A" &&
        return Star(78.15663333992994, 43.83321389573166; name="Gaia DR3 202357903847210880",
                    source_id="202357903847210880", epoch=2016.0, pmra=-0.9193745493858563, pmdec=-5.101328284551705,
                    parallax=0.8348808914164769, ra_err=0.0110, dec_err=0.0070, pmra_err=0.0151, pmdec_err=0.0115,
                    parallax_err=0.0122, mag=Dict("G" => 12.68, "V" => 12.93, "K" => 10.36))
    name == "B" &&
        return Star(91.56952782017574, 43.990014581484004; name="Gaia DR3 961233167113907328",
                    source_id="961233167113907328", epoch=2016.0, pmra=-1.9120392847271543, pmdec=-3.436637935988146,
                    parallax=0.4964806271737549, ra_err=0.0194, dec_err=0.0178, pmra_err=0.0251, pmdec_err=0.0188,
                    parallax_err=0.0230, mag=Dict("G" => 14.43, "V" => 13.98, "K" => 12.87))
    name == "C" &&
        return Star(88.62126849540857, 43.22482953922541; name="Gaia DR3 193308850132660224",
                    source_id="193308850132660224", epoch=2016.0, pmra=-0.2641131935021791, pmdec=-1.444301323430254,
                    parallax=0.4154285761151138, ra_err=0.0642, dec_err=0.0571, pmra_err=0.0798, pmdec_err=0.0580,
                    parallax_err=0.0690, mag=Dict("G" => 16.86, "V" => 16.45, "K" => 15.26))
    throw(ArgumentError("unknown event $name"))
end

"Ephemeris of Bienor for event `name` from the stored Horizons tables."
function bienor_ephemeris(name::AbstractString)
    rows = [split(l) for l in eachline(joinpath(BIENOR_DIR, BIENOR_EVENTS[name][3])) if !startswith(l, "#")]
    return EphemTable([r[1] for r in rows], parse.(Float64, getindex.(rows, 2)), parse.(Float64, getindex.(rows, 3)),
                      parse.(Float64, getindex.(rows, 4)); name="(54598) Bienor")
end

"""
    bienor_event(name; negatives=:all) -> Event

The occultation `name`: `"2019"` (2019 January 11, five positive and four negative chords),
`"A"` (2022 February 6), `"B"` (2022 December 26) or `"C"` (2023 February 14).
"""
function bienor_event(name::AbstractString; negatives::Symbol=:all)
    negatives in (:all, :none) || throw(ArgumentError("negatives must be :all or :none"))
    date, tc, _ = BIENOR_EVENTS[name]
    body = Body("(54598) Bienor", bienor_ephemeris(name); radius=75.0)
    ev = Event(bienor_star(name), body, date * "T" * tc; window=1200, name="Bienor, " * date)
    day = date * "T"
    for (n, lon, lat, alt, sens, ti, σi, te, σe) in BIENOR_CHORDS[name]
        obs = Observer(n, lon, lat, alt)
        if ti === nothing
            negatives === :all || continue
            push!(ev, Chord(n, obs; span=(day * "00:50:00", day * "01:17:00"), sensitivity=sens))
        else
            push!(ev, Chord(n, obs; immersion=day * ti, immersion_err=σi, emersion=day * te, emersion_err=σe,
                            sensitivity=sens))
        end
    end
    return ev
end
