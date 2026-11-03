# Stellar occultation by (38628) Huya on 2019 March 18
#
# Data from P. Santos-Sanz et al. (2022), "Physical properties of the trans-Neptunian
# object (38628) Huya from a multi-chord stellar occultation", A&A 664, A130
# (doi:10.1051/0004-6361/202141546, arXiv:2205.12882):
#   star: Table 3; stations: Table A.1; timings: Table A.3; published ellipse: Table A.4.
# Ephemeris: JPL Horizons (default solution), stored in huya_horizons_20190318.txt; the
# paper used JPL#28 + DE431, so the fitted centre (an offset from the ephemeris) differs
# from the published one while the shape does not.
#
# `huya_event()` returns the `Event` with all chords; `HUYA_PAPER` the published results.

using Occultations

const HUYA_DIR = @__DIR__

"Published best-fitting ellipse (Santos-Sanz et al. 2022, Table A.4) and other results."
const HUYA_PAPER = (a=(217.6, 3.5), b=(194.1, 6.1), position_angle=(55.2, 9.1), center_f=(2984.7, 3.2),
                    center_g=(-1850.9, 1.7), D_eq=(411.0, 7.3), velocity=8.07, distance_au=28.35,
                    fresnel_km=1.13, star_diameter_km=0.70, rotational_amplitude=(0.031, 0.006))

# chord number (Table A.3), station, latitude, longitude (E), altitude [m], exposure [s],
# immersion, σ, emersion, σ (UTC on 2019-03-18; `nothing` = not observed)
const HUYA_CHORDS = [
    (1, "TUBITAK-1 (Turkey)", "36 49 17.1", "30 20 08.0", 2470, 2.0, "00:52:34.387", 2.592, "00:53:00.703", 0.142),
    (2, "TUBITAK-2 (Turkey)", "36 49 29.1", "30 20 08.2", 2455, 15.0, "00:52:30.407", 5.500, "00:52:52.439", 5.500),
    (3, "Cluj-Feleacu (Romania)", "46 42 37.6", "23 35 35.7", 783, 4.0, "00:53:51.515", 0.115, "00:54:37.293", 0.089),
    (4, "ROASTERR-1 (Romania)", "46 49 15.6", "23 35 47.0", 390, 0.2, nothing, NaN, "00:54:38.748", 0.030),
    (5, "W-FAST Wise (Israel)", "30 35 48.0", "34 45 48.0", 900, 1.0, "00:51:21.416", 0.016, "00:52:07.918", 0.023),
    (6, "C18 Wise (Israel)", "30 35 48.0", "34 45 48.0", 900, 3.0, "00:51:21.110", 0.040, "00:52:07.851", 0.040),
    (7, "AIRA-1 Bucharest (Romania)", "44 24 43.2", "26 05 38.2", 83, 2.0, "00:53:30.014", 0.998, "00:54:18.863", 0.098),
    (8, "AIRA-2 Bucharest (Romania)", "44 24 48.0", "26 05 48.0", 81, 1.0, "00:53:31.111", 0.041, "00:54:18.814", 0.041),
    (9, "ISTEK Belde (Turkey)", "41 01 49.3", "29 02 33.6", 150, 0.8, "00:52:59.908", 0.062, "00:53:50.451", 0.083),
    (10, "Amateur Obs.-1 (Romania)", "44 55 05.0", "25 58 12.1", 167, 4.0, "00:53:33.562", 0.444, "00:54:23.630", 0.174),
    (11, "St. George (Romania)", "45 00 25.2", "25 58 42.2", 242, 3.0, "00:53:34.446", 0.083, "00:54:24.868", 0.077),
    (12, "Amateur Obs.-2 (Romania)", "44 26 36.5", "26 31 12.5", 65, 0.2, "00:53:30.484", 0.020, "00:54:21.217", 0.200),
    (13, "Stardust (Romania)", "45 38 30.0", "25 37 19.0", 597, 8.0, "00:53:38.100", 0.150, "00:54:29.800", 1.850),
    (14, "Stardreams (Romania)", "45 12 13.3", "26 02 44.2", 379, 1.0, "00:53:37.111", 0.098, "00:54:29.122", 0.099),
    (15, "Kraar/Weizmann (Israel)", "31 54 29.0", "34 48 45.0", 107, 1.0, "00:51:30.046", 0.100, "00:52:23.578", 1.148),
    (16, "Amateur Obs.-4 (Romania)", "44 45 41.5", "27 20 25.1", 41, 0.2, "00:53:33.840", 0.027, "00:54:25.060", 0.027),
    (17, "Bacau (Romania)", "46 33 56.3", "26 54 15.0", 170, 2.0, "00:53:46.520", 3.400, "00:54:33.919", 0.083),
    (18, "Galati (Romania)", "45 25 07.9", "28 01 57.0", 31, 20.0, "00:53:40.500", 0.440, "00:54:25.420", 0.440),
    (19, "Barlad (Romania)", "46 13 54.1", "27 40 10.2", 70, 0.5, "00:53:51.192", 0.047, "00:54:31.753", 0.073),
]

# negative chords (Table A.1); the observed interval is not given in the paper: 00:50–00:56 UT
# brackets the event at all stations
const HUYA_NEGATIVES = [
    ("QOS (Ukraine)", "48 50 54.0", "26 43 12.0", 352, 0.2),
    ("Cukurova (Turkey)", "37 03 23.0", "35 21 04.0", 130, 90.0),
]

"Gaia DR2 4352760586390566400, at the occultation epoch (Santos-Sanz et al. 2022, Table 3)."
huya_star() = Star("16 41 06.4260", "-06 43 34.5756"; name="Gaia DR2 4352760586390566400",
                   source_id="4352760586390566400", epoch=2019.211, pmra=-4.506, pmdec=-14.340, parallax=1.5717,
                   pmra_err=0.090, pmdec_err=0.055, parallax_err=0.0478, ra_err=0.1474, dec_err=0.1143,
                   mag=Dict("B" => 12.296, "V" => 11.751, "R" => 11.380, "J" => 10.541, "H" => 10.384,
                            "K" => 10.294, "G" => 11.538), diameter=0.0342)

"Ephemeris of Huya from the stored Horizons table."
function huya_ephemeris()
    rows = [split(l) for l in eachline(joinpath(HUYA_DIR, "huya_horizons_20190318.txt")) if !startswith(l, "#")]
    return EphemTable([r[1] for r in rows], parse.(Float64, getindex.(rows, 2)), parse.(Float64, getindex.(rows, 3)),
                      parse.(Float64, getindex.(rows, 4)); name="(38628) Huya")
end

"""
    huya_event(; negatives=true) -> Event

The 2019 March 18 occultation by Huya with the 19 positive chords of Table A.3 (chord 4 only
emersion; the unused chord "?" is left out) and the two closest negative chords.
"""
function huya_event(; negatives::Bool=true)
    body = Body("(38628) Huya", huya_ephemeris(); radius=203.0)
    ev = Event(huya_star(), body, "2019-03-18T00:53:00"; window=1200, name="Huya, 2019-03-18")
    day = "2019-03-18T"
    for (n, name, lat, lon, alt, texp, ti, σi, te, σe) in HUYA_CHORDS
        obs = Observer(name, lon, lat, alt)
        push!(ev, Chord("#$n $name", obs; immersion=ti === nothing ? nothing : day * ti, immersion_err=isnan(σi) ? 0.0 : σi,
                        emersion=day * te, emersion_err=σe, sensitivity=texp))
    end
    if negatives
        for (name, lat, lon, alt, texp) in HUYA_NEGATIVES
            push!(ev, Chord(name, Observer(name, lon, lat, alt); span=(day * "00:50:00", day * "00:56:00"), sensitivity=texp))
        end
    end
    return ev
end
