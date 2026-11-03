# Stellar occultation by (24835) 1995 SM55 on 2024 February 25
#
# Data from J. L. Ortiz et al. (2025), "A high geometric albedo and small size for the Haumea
# cluster member (24835) 1995 SM55 determined from a stellar occultation and photometric
# observations", A&A 703, A147 (doi:10.1051/0004-6361/202556498):
#   star: Table 1 (Gaia DR3 position at the occultation epoch); timings: Table 2;
#   stations: Table A.1; limb fit: Table 3; light curve, 3-D shape, albedo: Sects. 3–4.
# Ephemeris: JPL Horizons (default solution), stored in sm55_horizons_20240225.txt; the paper
# refers the fitted centre to NIMA v10.
#
# Notes on the transcription:
# * Table 2 and Table A.1 number the stations differently; the chords are matched by name
#   (Table 2 "Türkiye National Observatories-150cm" = RTT150, "-100cm" = T100; Wise H80 = the
#   80 cm telescope, Wise T45.7 = the 45.7 cm telescope).
# * The observed intervals of the negative chords are not tabulated; 18:08–18:24 UTC, which
#   brackets the event at all sites (18:13–18:18 at the positive ones), is used.
# * The sensitivity of a chord is its exposure plus dead time (Table A.1; the dead time only
#   where no exposure is given).
# * The text (Sect. 2.2) states that all positive detections used exposures ≥ 2 s; Table A.1
#   lists 0.1 s (RTT150, Chalin) and 0.3 s (T100).
#
# `sm55_event()` returns the `Event`; `SM55_PAPER` the published results.

using Occultations

const SM55_DIR = @__DIR__

"Published results of Ortiz et al. (2025)."
const SM55_PAPER = (
    # Table 3: limb fit (centre relative to NIMA v10; φ′ of the minor axis from north)
    center_f=(303.8, 0.17), center_g=(-787.8, 0.4), a=(104.3, 0.4), b=(83.5, 0.5), position_angle=(44.1, 0.4),
    oblateness=(0.200, 0.002), D_eq=(186.7, 1.8), chi2_pdf=1.55, npoints=14,
    # rotational light curve (Sect. 3.2): P = 52.52 h (or 26.26 h), Δm = 0.05 ± 0.02 mag,
    # rotational phase since maximum light at the occultation 0.09
    period_h=52.52, amplitude=(0.05, 0.02), rotational_phase=0.09,
    # 3-D solution of the grid search (Sect. 4.1); aspect from the pole
    axes_3d=((106.0, 2.0, 2.0), (92.0, 12.0, 12.0), (76.0, 10.0, 10.0)), aspect=(52.0, 38.0, 28.0),
    D_eq_volume=(181.0, 12.0),
    # absolute magnitude, albedo, astrometry (Sects. 3.1, 4.2, 4.3)
    H_V=(4.55, 0.03), albedo=(0.80, 0.04),
    astrometry=(epoch="2024-02-25T18:13:53.460", ra="03 40 30.5344080", dec="+39 08 06.708383", ra_err=0.1, dec_err=0.1),
    velocity=10.42, star_diameter_mas=0.0228, fresnel_km_700nm=1.32,
)

# station, longitude (E), latitude, altitude [m], exposure + dead time [s],
# immersion, σ, emersion, σ (UTC, 2024-02-25; `nothing` for a negative chord)
const SM55_CHORDS = [
    ("Cluj-Napoca", "23 35 35.6479", "46 42 37.557", 783.40, 2.6, "18:14:40.87", 0.82, "18:14:57.66", 0.40),
    ("Wise H80", "34 45 44.1366", "30 35 48.5862", 862.27, 8.0, "18:17:41.091", 0.014, "18:18:00.3", 2.4),
    ("Wise T45.7", "34 45 44.1366", "30 35 48.5862", 862.27, 7.1, "18:17:40.344", 0.035, "18:17:57.8", 2.1),
    ("TUG RTT150", "30 20 07.25715", "36 49 32.1948", 2458.59, 3.88, "18:16:32.59", 0.95, "18:16:48.3", 1.6),
    ("TUG T100", "30 20 07.98456", "36 49 17.0672", 2538.72, 2.66, "18:16:36.9", 1.2, "18:16:52.9", 1.2),
    ("Berthelot", "22 53 19.6894", "45 36 59.2515", 386.36, 5.7, "18:14:48.534", 0.024, "18:15:02.8", 2.3),
    ("Chalin", "16 02 33.4403", "52 36 14.027", 65.00, 0.127, "18:13:26.209", 0.004, "18:13:29.618", 0.006),
    ("Wendelstein 2.1m", "12 00 47.88", "47 42 14.5368", 1931.0, 1.4336, nothing, NaN, nothing, NaN),
    ("KAO", "31 49 37.2", "29 56 02.4", 476.0, 12.3, nothing, NaN, nothing, NaN),
    ("ELTE GAO MKK", "16 36 12.1097", "47 15 28.5471", 209.0, 0.4001, nothing, NaN, nothing, NaN),
    ("Bavarian Public Obs.", "11 36 25.9992", "48 07 18.9984", 500.0, 0.04, nothing, NaN, nothing, NaN),
    ("Ondrejov", "14 47 01.104", "49 54 38.016", 528.0, 4.2, nothing, NaN, nothing, NaN),
    ("Teplice", "13 50 48.3", "50 38 17.9556", 277.94, 0.3, nothing, NaN, nothing, NaN),
    ("San Marcello Pistoiese", "10 48 15.12", "44 03 46.9296", 965.41, 2.0, nothing, NaN, nothing, NaN),
    ("Crni Vrh", "14 04 16.6242", "45 56 45.0688", 713.03, 1.0, nothing, NaN, nothing, NaN),
    ("Montarrenti", "11 11 02.04", "43 13 57", 347.0, 2.3, nothing, NaN, nothing, NaN),
    ("Monte Agliale", "10 30 53.496", "43 59 43.008", 760.0, 5.0, nothing, NaN, nothing, NaN),
    ("Harpoint", "13 21 06.12", "47 54 33.12", 700.0, 0.5, nothing, NaN, nothing, NaN),
    ("Radebeul", "13 37 19.92", "51 06 58.6728", 185.0, 1.2, nothing, NaN, nothing, NaN),
    ("Wendelstein 0.43m", "12 00 47.88", "47 42 14.5368", 1931.0, 0.4, nothing, NaN, nothing, NaN),
    ("Starhopper", "25 46 08.0004", "45 51 56.0016", 588.0, 5.5, nothing, NaN, nothing, NaN),
    ("IST40", "28 57 56.5848", "41 00 42.2964", 60.0, 2.0, nothing, NaN, nothing, NaN),
    ("Karrenkneul", "09 21 38.6691", "51 23 36.8938", 266.49, 0.3, nothing, NaN, nothing, NaN),
    ("DROT", "18 42 09.33781", "49 24 15.2341", 635.0, 0.06, nothing, NaN, nothing, NaN),
    ("Ondrejov Cosmic Lab", "14 46 47.7012", "49 54 36.2016", 524.0, 0.1001, nothing, NaN, nothing, NaN),
    ("Kaufering", "10 50 57.9958", "48 05 22.7793", 596.0, 0.101, nothing, NaN, nothing, NaN),
    ("Strasice", "13 44 56.3093", "49 44 35.7118", 536.61, 0.15, nothing, NaN, nothing, NaN),
    ("Siena", "11 20 12.12", "43 18 45", 290.0, 4.0, nothing, NaN, nothing, NaN),
    ("Plzen-Valcha", "13 19 55.4668", "49 42 26.4403", 325.38, 0.1, nothing, NaN, nothing, NaN),
    ("Valasske Mezirici", "17 58 25.2984", "49 27 49.2984", 338.0, 0.6, nothing, NaN, nothing, NaN),
    ("Volkssternwarte Muenchen", "11 36 25.5974", "48 07 19.0733", 565.0, 0.08, nothing, NaN, nothing, NaN),
    ("Frosinone", "13 20 15.7991", "41 39 11.6387", 198.98, 0.394, nothing, NaN, nothing, NaN),
    ("Brtonigla", "13 37 47", "45 23 04.99999", 125.0, 1.0, nothing, NaN, nothing, NaN),
    ("Foligno", "12 42 17.928", "42 57 36.972", 220.0, 1.1, nothing, NaN, nothing, NaN),
    ("ADYU60", "38 13 31.188", "37 45 06.05999", 690.0, 1.05, nothing, NaN, nothing, NaN),
    ("TURKSAT", "32 48 14.9652", "39 38 11.8752", 950.0, 2.0, nothing, NaN, nothing, NaN),
    ("UZAYMER", "35 20 52.71", "37 03 20.088", 126.0, 5.0, nothing, NaN, nothing, NaN),
    ("Eskisehir", "30 27 38.4804", "39 53 07.6992", 1005.0, 10.0, nothing, NaN, nothing, NaN),
]

"Gaia DR3 224489389287197184 at the occultation epoch (Ortiz et al. 2025, Table 1)."
sm55_star() = Star("03 40 30.5277", "+39 08 06.647"; name="Gaia DR3 224489389287197184",
                   source_id="224489389287197184", epoch=julian_year(Epoch("2024-02-25T18:14")),
                   ra_err=0.1, dec_err=0.1, mag=Dict("G" => 12.4), diameter=SM55_PAPER.star_diameter_mas)

"Ephemeris of 1995 SM55 from the stored Horizons table."
function sm55_ephemeris()
    rows = [split(l) for l in eachline(joinpath(SM55_DIR, "1995sm55_horizons_20240225.txt")) if !startswith(l, "#")]
    return EphemTable([r[1] for r in rows], parse.(Float64, getindex.(rows, 2)), parse.(Float64, getindex.(rows, 3)),
                      parse.(Float64, getindex.(rows, 4)); name="(24835) 1995 SM55")
end

"""
    sm55_event(; negatives=:all) -> Event

The occultation of 2024 February 25 with the seven positive chords (Table 2) and the
negative chords of Table A.1 (`negatives = :all` or `:none`).
"""
function sm55_event(; negatives::Symbol=:all)
    negatives in (:all, :none) || throw(ArgumentError("negatives must be :all or :none"))
    body = Body("(24835) 1995 SM55", sm55_ephemeris(); radius=95.0)
    ev = Event(sm55_star(), body, "2024-02-25T18:15:00"; window=900, name="1995 SM55, 2024-02-25")
    day = "2024-02-25T"
    for (name, lon, lat, alt, sens, ti, σi, te, σe) in SM55_CHORDS
        obs = Observer(name, lon, lat, alt)
        if ti === nothing
            negatives === :all || continue
            push!(ev, Chord(name, obs; span=(day * "18:08:00", day * "18:24:00"), sensitivity=sens))
        else
            push!(ev, Chord(name, obs; immersion=day * ti, immersion_err=σi, emersion=day * te, emersion_err=σe,
                            sensitivity=sens))
        end
    end
    return ev
end
