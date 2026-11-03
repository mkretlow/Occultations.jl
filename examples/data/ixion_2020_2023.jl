# Stellar occultations by (28978) Ixion, 2020–2023
#
# Data from Y. Kilic et al. (2026), "Constraining the size, shape, and albedo of the large
# trans-Neptunian object (28978) Ixion with multi-chord stellar occultations", A&A 707, A70
# (doi:10.1051/0004-6361/202557970):
#   stars: Table 2 (Gaia DR3 positions at the occultation epochs); events: Table 3;
#   timings: Table A.1; stations: Table D.1; limb fits: Table 4, Eq. (4); astrometry: Table 5.
# Ephemerides: JPL Horizons (default solution), stored in ixion_horizons_2020_2023.txt; the
# paper refers the fitted centres to NIMA v8–v12.
#
# The six multi-chord events of the paper's global fit are included (the single-chord events
# of 2021-08-17 and 2022-06-02 only give astrometry with the global shape).
#
# Notes on the transcription:
# * The observed intervals of the negative chords are not tabulated; ±5 min around the event
#   is used. The sensitivity of a chord is its exposure plus dead time (Table D.1).
# * "CAO" of Table A.1 is the Campocatino Austral Observatory of Table D.1; Lowell
#   Observatory (2020-10-13) recorded only the emersion.
# * Table 4, 28 April 2021: R_equiv = 386.3 km and ε′ = 0.231 do not follow from a = 391.9 and
#   b = 322.2 km (√(ab) = 355.3 km, 1 − b/a = 0.178); the other columns are consistent.
#
# `ixion_events()` returns the six `Event`s; `IXION_PAPER` the published results.

using Occultations

const IXION_DIR = @__DIR__

"Published results of Kilic et al. (2026)."
const IXION_PAPER = (
    # global fit, Eq. (4)–(7): ϕ (SORA's position angle: the minor axis from north towards
    # east, as the per-event values of Table 4 show); σ_model = 4 km
    a=(363.42, 3.53, 3.85), b=(333.98, 7.07, 4.96), phi_major=(110.97, 7.00, 6.54), R_equiv=(348.39, 5.37, 4.43),
    oblateness=(0.081, 0.004, 0.010), chi2_nu=1.17, ndata=56, dof=41, sigma_model=4.0,
    # Table 4: per-event fits (a, b, ϕ, R_equiv, χ²pdf); Table 5: astrometry
    events=Dict(
        "2020-08-17" => (a=(375.1, 11.7), b=(330.2, 20.1), phi=(14.86, 8.10), R_equiv=(351.9, 12.0), chi2_pdf=0.816,
                         astrometry=("01:37:11.820", "17 49 53.5853", "-29 48 12.3881", 0.284, 0.282)),
        "2020-10-13" => (a=(366.8, 6.2), b=(346.1, 7.4), phi=(34.11, 15.70), R_equiv=(356.3, 4.9), chi2_pdf=0.926,
                         astrometry=("01:57:46.320", "17 50 20.3585", "-29 43 32.9060", 0.141, 0.274)),
        "2021-04-28" => (a=(391.9, 9.4), b=(322.2, 17.2), phi=(35.25, 17.50), R_equiv=(386.3, 59.8), chi2_pdf=0.393,
                         astrometry=("07:36:45.400", "18 05 34.6798", "-30 09 09.7432", 0.567, 0.240)),
        "2021-05-20" => (a=(375.9, 20.8), b=(338.2, 41.6), phi=(-4.23, 63.65), R_equiv=(356.6, 24.1), chi2_pdf=0.303,
                         astrometry=("06:43:21.420", "18 04 05.2324", "-30 13 44.4493", 1.200, 3.868)),
        "2022-06-30" => (a=(389.5, 16.8), b=(325.6, 31.5), phi=(31.35, 13.51), R_equiv=(356.1, 18.9), chi2_pdf=0.098,
                         astrometry=("00:22:49.920", "18 06 53.6124", "-30 46 56.5116", 0.642, 0.351)),
        "2023-07-28" => (a=(371.1, 33.5), b=(334.0, 47.8), phi=(-0.04, 89.81), R_equiv=(352.0, 29.8), chi2_pdf=0.938,
                         astrometry=("06:13:19.080", "18 11 03.2660", "-31 15 22.5483", 0.399, 2.046)),
    ),
    # Table 3: distance [au], velocity [km/s]
    geometry=Dict("2020-08-17" => (38.37, 11.61), "2020-10-13" => (39.26, 16.12), "2021-04-28" => (38.15, 13.54),
                  "2021-05-20" => (37.89, 20.66), "2022-06-30" => (37.46, 24.66), "2023-07-28" => (37.34, 20.21)),
    H_V=(3.845, 0.006), albedo=(0.106, 0.003),
)

# date => (star: Gaia DR3 id, RA, Dec, G, diameter [mas]), reference time, chords:
# station, longitude, latitude, altitude [m], exposure + dead time [s], immersion, σ, emersion, σ
const IXION_EVENTS = [
    ("2020-08-17", ("4056438835548480896", "17 49 53.5830", "-29 48 12.589", 14.52, 0.0530), "01:38:00", [
        ("Chester", "-83 12 07.2", "32 22 15.0996", 104.0, 0.3, "01:38:02.48", 0.10, "01:38:31.71", 0.09),
        ("Naperville", "-88 07 00.0084", "41 45 32.4", 230.0, 5.0, "01:38:21.50", 1.10, "01:39:22.24", 1.06),
        ("Lake Placid", "-74 01 59.4192", "44 15 18.3996", 608.0, 3.00088, "01:36:59.90", 0.90, "01:37:38.54", 0.92)]),
    ("2020-10-13", ("4056440205544338944", "17 50 20.3595", "-29 43 33.087", 10.31, 0.6750), "02:01:00", [
        ("Lowell (TiMo)", "-111 40 01.4628", "35 12 10.3644", 2198.0, 0.05, nothing, NaN, "02:01:30.22", 0.02),
        ("Scottsdale", "-111 52 07.2984", "33 49 00.0984", 843.0, 0.033, "02:00:47.33", 0.02, "02:01:31.79", 0.02),
        ("Tempe", "-111 57 20.574", "33 23 24.3348", 358.0, 0.134, "02:00:48.17", 0.03, "02:01:30.69", 0.02),
        ("Glendale", "-112 12 46.0008", "33 41 39.0012", 400.0, 0.134, "02:00:45.99", 0.02, "02:01:29.34", 0.76),
        ("Tucson", "-110 57 46.62", "32 18 28.314", 750.0, 0.067, "02:00:57.68", 0.12, "02:01:34.16", 0.09),
        ("Running Springs", "-117 07 52.3992", "34 13 07.8996", 1874.0, 0.033, "02:00:21.52", 0.01, "02:01:05.06", 0.39),
        ("Carefree", "-111 57 07.974", "33 48 42.858", 654.0, 0.033, "02:00:47.00", 0.04, "02:01:30.59", 0.01),
        ("Fountain Hills", "-111 43 35.0004", "33 37 21", 515.0, 0.033, "02:00:48.66", 0.01, "02:01:31.61", 0.01),
        ("Gardnerville", "-119 40 20.3016", "38 53 23.4996", 1524.0, 0.334, nothing, NaN, nothing, NaN),
        ("Clear Lake Shores", "-95 02 06.57996", "29 32 55.5101", 15.8, 0.334, nothing, NaN, nothing, NaN)]),
    ("2021-04-28", ("4050000576486883840", "18 05 34.6814", "-30 09 09.703", 14.84, 0.0471), "07:37:00", [
        ("OPD", "-45 34 57.54", "-22 32 07.7532", 1810.71, 0.15, "07:34:38.43", 0.005, "07:34:56.47", 0.022),
        ("CTIO One-Meter", "-70 48 21.69", "-30 10 07.7916", 2201.0, 1.0, "07:37:49.83", 0.83, "07:38:02.49", 0.32),
        ("La Canelilla", "-70 47 45", "-30 32 03", 1548.0, 0.50112, "07:37:44.56", 0.19, "07:38:09.72", 0.19),
        ("CASLEO ASH", "-69 18 23.868", "-31 47 13.2", 2591.0, 12.0, "07:37:28.12", 1.12, "07:38:15.29", 0.57),
        ("PROMPT-6", "-70 48 18.8208", "-30 10 03.72", 2166.0, 1.0, "07:37:50.05", 0.35, "07:38:02.51", 0.01),
        ("Reconquista", "-59 38 36.6086", "-29 08 25.1826", 50.0, 5.0, "07:36:21.60", 0.78, "07:37:02.93", 0.81),
        ("TRAPPIST-South", "-70 44 21.84", "-29 15 16.56", 2315.0, 3.7, nothing, NaN, nothing, NaN),
        ("A.A.A.A.", "-58 26 04.6", "-34 36 18.6984", 39.4, 5.0, nothing, NaN, nothing, NaN),
        ("El Catalejo", "-64 19 27.48", "-36 38 15.99", 182.0, 2.0, nothing, NaN, nothing, NaN),
        ("Los Cabezones", "-64 17 16.8299", "-36 38 08.28996", 180.0, 5.0, nothing, NaN, nothing, NaN),
        ("Cruz del Sur", "-58 34 24.276", "-34 40 10.56", 39.0, 6.0, nothing, NaN, nothing, NaN)]),
    ("2021-05-20", ("4050180170543871104", "18 04 05.2362", "-30 13 44.615", 14.18, 0.0320), "06:45:00", [
        ("Fountain Hills", "-111 43 39.5879", "33 37 27.8868", 518.0, 2.0, "06:45:14.85", 0.39, "06:45:48.62", 0.40),
        ("Scottsdale", "-111 52 07.2984", "33 49 00.0984", 843.0, 0.53386, "06:45:14.84", 0.03, "06:45:48.97", 0.03),
        ("Clear Lake Shores", "-95 02 06.57996", "29 32 55.5101", 15.8, 1.068, "06:44:19.36", 0.27, "06:44:52.59", 0.21),
        ("Bartlett Lake", "-111 49 58.08", "33 50 54.8412", 996.0, 0.0333, "06:45:14.26", 0.26, "06:45:48.96", 0.13),
        ("Westport", "-73 19 39.2999", "41 10 15.8999", 88.0, 1.0, nothing, NaN, nothing, NaN)]),
    ("2022-06-30", ("4049905056379149824", "18 06 53.6148", "-30 46 56.476", 15.00, 0.0335), "00:25:00", [
        ("UEPG", "-50 05 56.7901", "-25 05 22.671", 923.0, 1.0, "00:25:01.67", 0.23, "00:25:22.52", 0.24),
        ("Orion", "-47 12 40.5144", "-23 34 10.2036", 884.0, 6.0, "00:24:48.17", 2.46, "00:25:15.84", 2.63),
        ("YPO", "-45 31 05.99988", "-22 41 49", 1640.0, 3.2, "00:24:44.11", 1.37, "00:25:10.19", 1.47),
        ("UTFPR", "-49 18 57.1841", "-25 30 31.7373", 937.67, 3.987, "00:25:00.41", 0.84, "00:25:19.50", 1.78),
        ("SONEAR", "-44 47 05.82", "-20 42 54.27", 1113.0, 3.0, nothing, NaN, nothing, NaN)]),
    ("2023-07-28", ("4049162851669604096", "18 11 03.2654", "-31 15 22.477", 14.53, 0.0500), "06:09:30", [
        ("SOAR", "-70 44 21.12", "-30 14 16.89", 2693.95, 0.2, "06:09:16.73", 0.13, "06:09:45.77", 0.14),
        ("CAO", "-70 45 53.9015", "-30 28 15.2905", 1525.0, 1.0, "06:09:17.98", 0.03, "06:09:44.95", 0.05),
        ("Danish/ESO", "-70 44 01.464", "-29 15 31.7592", 2345.44, 0.1, "06:09:12.23", 0.04, "06:09:45.28", 0.02),
        ("OPD 60cm", "-45 26 43.8041", "-22 32 49.5658", 1112.15, 0.3007, nothing, NaN, nothing, NaN),
        ("CASLEO HSH", "-69 18 24.12", "-31 47 13.2", 2591.0, 9.0, nothing, NaN, nothing, NaN),
        ("UC Observatory", "-70 32 04.1352", "-33 16 08.9364", 1475.2, 5.0, nothing, NaN, nothing, NaN),
        ("Astroquinta", "-71 14 55.7412", "-32 53 00.2796", 154.14, 12.7, nothing, NaN, nothing, NaN)]),
]

"Ephemeris of Ixion around the event of `date` (\"YYYY-MM-DD\") from the stored Horizons table."
function ixion_ephemeris(date::AbstractString)
    rows = [split(l) for l in eachline(joinpath(IXION_DIR, "ixion_horizons_2020_2023.txt"))
            if !startswith(l, "#") && startswith(l, date)]
    return EphemTable([r[1] for r in rows], parse.(Float64, getindex.(rows, 2)), parse.(Float64, getindex.(rows, 3)),
                      parse.(Float64, getindex.(rows, 4)); name="(28978) Ixion")
end

"""
    ixion_events() -> Vector{Event}

The six multi-chord occultations of 2020–2023 (Table A.1, positive and negative chords of
Table D.1).
"""
function ixion_events()
    map(IXION_EVENTS) do (date, st, tmid, chords)
        star = Star(st[2], st[3]; name="Gaia DR3 " * st[1], source_id=st[1], epoch=julian_year(Epoch(date * "T" * tmid)),
                    mag=Dict("G" => st[4]), diameter=st[5])
        body = Body("(28978) Ixion", ixion_ephemeris(date); radius=350.0)
        ev = Event(star, body, date * "T" * tmid; window=1200, name="Ixion, " * date)
        t0 = Epoch(date * "T" * tmid)
        span = (t0 - 300.0, t0 + 300.0)
        for (name, lon, lat, alt, sens, ti, σi, te, σe) in chords
            obs = Observer(name, lon, lat, alt)
            if ti === nothing && te === nothing
                push!(ev, Chord(name, obs; span=span, sensitivity=sens))
            elseif ti === nothing                          # emersion only
                push!(ev, Chord(name, obs; emersion=date * "T" * te, emersion_err=σe, sensitivity=sens))
            else
                push!(ev, Chord(name, obs; immersion=date * "T" * ti, immersion_err=σi, emersion=date * "T" * te,
                                emersion_err=σe, sensitivity=sens))
            end
        end
        ev
    end
end
