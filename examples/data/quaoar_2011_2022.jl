# Stellar occultations by (50000) Quaoar on 2011 May 4 and 2022 August 9
#
# Data from
# * F. Braga-Ribas et al. (2013), "The size, shape, albedo, density, and atmospheric limit of
#   transneptunian object (50000) Quaoar from multi-chord stellar occultations", ApJ 773, 26
#   (doi:10.1088/0004-637X/773/1/26): star: Eq. (1); stations and timings: Table 1; limb fits:
#   Tables 2 and 3;
# * C. L. Pereira et al. (2023), "The two rings of (50000) Quaoar", A&A 673, L4
#   (doi:10.1051/0004-6361/202346365): star: Table 1; body and ring geometry: Table 2;
#   timings of the body: Table A.1; stations: Table B.1; ring detections: Table E.1;
# * K. Arimatsu et al. (2019), "New constraint on the atmosphere of (50000) Quaoar from a
#   stellar occultation", AJ 158, 236 (doi:10.3847/1538-3881/ab5058): a single chord on
#   2019 June 28 from Kiso Observatory, duration 45.50 ± 0.05 s; absolute times are not given,
#   so the event is not included.
# Ephemeris: JPL Horizons (default solution), stored in quaoar_horizons_<yyyymmdd>.txt.
#
# Notes on the transcription:
# * 2011: the star is the paper's astrometry (UCAC2 26257103, ± 7 and 14 mas). The longitude of
#   Pico dos Dias (45° 43′ 58″ W in Table 1) is taken as the site's 45° 34′ 57.5″ W. Observed
#   intervals of the negative chords: as given (Pico dos Dias, Salto), else 02:30–02:55 UTC.
#   The sensitivity of a chord is its exposure plus readout time.
# * 2022: the two bands of Gemini North (r′, z′) are separate chords, as in the paper. Ring
#   detections are used by their mid-times (Table E.1), labelled `:ring1` (Q1R, the outer ring)
#   and `:ring2` (Q2R). Woodside (inconclusive) is left out; Summerland is the negative chord.
#   Its observed interval is not given: 06:25–06:45 UTC.
#
# `quaoar_event(year)` returns the `Event` of 2011 or 2022; `QUAOAR_PAPER` the published results.

using Occultations

const QUAOAR_DIR = @__DIR__

"Published results of Braga-Ribas et al. (2013) and Pereira et al. (2023)."
const QUAOAR_PAPER = (
    # 2011 (Table 2): crater, mountain (apparent ellipses) and Maclaurin solutions
    e2011=(crater=(a=(656.0, 12.0), oblateness=(0.29, 0.04), position_angle=(21.0, 3.0), R_equiv=(552.0, 14.0)),
           mountain=(a=(568.0, 12.0), oblateness=(0.14, 0.05), position_angle=(7.5, 11.0), R_equiv=(524.0, 20.0)),
           maclaurin=(a=(569.0, 24.0), R_equiv=(555.0, 2.5))),
    # 2022 (Table 2): body limb (P: position angle of the minor axis), Q1R and Q2R
    e2022=(a=(579.5, 4.0), oblateness=(0.12, 0.01), position_angle=(345.2, 1.2), R_equiv=(543.0, 2.0), chi2_pdf=0.986,
           position=(epoch="2022-08-09T06:34:03.560", ra="18 21 42.8677703", dec="-15 12 45.829691"),
           Q1R=(radius=(4057.0, 6.0), B=(-20.0, 0.3), P=(350.2, 0.2), pole=((259.82, 0.23), (53.45, 0.30)),
                mirror_pole=(150.0, 79.35)),
           Q2R=(radius=(2520.0, 20.0),), velocity=17.57, distance_au=41.983157),
    arimatsu_2019=(duration=(45.50, 0.05), max_duration=45.1),
)

# station, longitude (E), latitude, altitude [m], sensitivity [s],
# immersion, σ, emersion, σ (UTC; `nothing`: negative chord, with its observed interval)
const QUAOAR_CHORDS_2011 = [
    ("Ponta Grossa UEPG", "-50 05 56", "-25 05 22", 909.0, 18.25, "02:41:47.2", 2.7, "02:42:42.1", 3.1),
    ("Florianopolis UFSC", "-48 31 21", "-27 36 12", 20.0, 12.75, "02:41:28.7", 3.6, "02:42:31.9", 3.8),
    ("San Pedro de Atacama Harlingten", "-68 10 48", "-22 57 12", 2397.0, 4.1, "02:42:33.1", 0.13, "02:43:31.4", 0.15),
    ("San Pedro de Atacama ASH2", "-68 10 48", "-22 57 12", 2397.0, 16.1, "02:42:32.5", 1.0, "02:43:30.3", 4.3),
    ("Armazones", "-70 11 46", "-24 35 52", 2706.0, 9.4, "02:42:41.5", 1.2, "02:43:26.3", 1.8),
    ("Rivera", "-55 36 00", "-30 51 48", 200.0, 6.3, "02:41:58.2", 0.4, "02:42:32.5", 0.4),
    ("CEAMIG-REA", "-43 59 51", "-19 49 49", 825.0, 9.0, nothing, ("02:30:00", "02:55:00"), nothing, NaN),
    ("Pico dos Dias", "-45 34 57.5", "-22 32 08", 1864.0, 2.0, nothing, ("02:26:10", "02:58:20"), nothing, NaN),
    ("Salto OLASU", "-57 58 43", "-31 23 33", 37.0, 6.7, nothing, ("02:20:29", "02:54:04"), nothing, NaN),
    ("TRAPPIST", "-70 44 22", "-29 15 17", 2315.0, 7.55, nothing, ("02:30:00", "02:55:00"), nothing, NaN),
    ("Los Molinos", "-56 11 25", "-34 45 19", 80.0, 5.0, nothing, ("02:30:00", "02:55:00"), nothing, NaN),
    ("CASLEO", "-69 17 45", "-31 47 56", 2492.0, 10.5, nothing, ("02:30:00", "02:55:00"), nothing, NaN),
    ("Santa Martina", "-70 32 04", "-33 16 09", 1450.0, 1.0, nothing, ("02:30:00", "02:55:00"), nothing, NaN),
]

# 2022: station, longitude (E), latitude, altitude [m], cycle [s], immersion, σ, emersion, σ;
# ring detections (mid-time, σ, ring)
const QUAOAR_CHORDS_2022 = [
    ("CFHT", "-155 28 07.93", "19 49 31.01", 4206.1, 0.113, "06:34:23.34", 0.02, "06:35:11.42", 0.01,
     [("06:31:15.8", 0.2, :ring1), ("06:38:32.581", 0.003, :ring1), ("06:32:43.9", 0.2, :ring2), ("06:37:00.20", 0.01, :ring2)]),
    ("Gemini North z'", "-155 28 08.6", "19 49 26", 4213.0, 0.101, "06:34:23.29", 0.01, "06:35:11.29", 0.01,
     [("06:31:15.72", 0.04, :ring1), ("06:38:32.444", 0.003, :ring1), ("06:32:43.61", 0.06, :ring2), ("06:37:00.09", 0.02, :ring2)]),
    ("Gemini North r'", "-155 28 08.6", "19 49 26", 4213.0, 0.101, "06:34:23.27", 0.03, "06:35:11.28", 0.02,
     [("06:31:15.7", 0.1, :ring1), ("06:38:32.429", 0.005, :ring1), ("06:36:59.93", 0.02, :ring2)]),
    ("TUHO", "-156 15 27", "20 42 25.14", 3051.8, 1.634, "06:34:23.5", 0.2, "06:35:18.1", 0.2,
     [("06:38:37.4", 0.7, :ring1)]),
    ("TAOS II", "-115 27 47.88", "31 02 36.26", 2824.8, 0.205, "06:30:38.3", 0.2, "06:31:35.1", 0.1,
     [("06:34:55.8", 0.5, :ring1)]),
    ("Dunrhomin", "-105 09 46.57", "40 15 08.46", 1592.8, 2.0, "06:30:03.3", 0.4, "06:31:02.8", 0.9, []),
    ("Sommers-Bausch", "-105 15 47.09", "40 00 13.32", 1645.0, 0.9, "06:30:02.6", 0.2, "06:31:03.8", 0.5, []),
    ("UCSC", "-122 04 46.20", "37 01 03.47", 344.0, 0.534, "06:31:20.3", 0.1, "06:32:00.9", 0.1, []),
    ("Nederland", "-105 26 44.05", "39 59 13.95", 2492.6, 2.0, "06:30:03.4", 0.5, "06:31:04.4", 0.8, []),
    ("Bonny Doon", "-122 08 18.62", "37 03 03.12", 488.4, 0.26, "06:31:20.7", 0.4, "06:32:01.3", 0.4, []),
    ("Summerland", "-119 40 12", "49 36 00", 485.0, 3.0, nothing, NaN, nothing, NaN, []),
]

"The occulted star of 2011 (the paper's astrometry) or 2022 (Gaia DR3 at the event)."
quaoar_star(year::Integer) = year == 2011 ?
    Star("17 28 50.801", "-15 27 42.770"; name="UCAC2 26257103", epoch=julian_year(Epoch("2011-05-04T02:42")),
         ra_err=7.0, dec_err=14.0, mag=Dict("R" => 15.8)) :
    Star("18 21 42.86965", "-15 12 45.9639"; name="Gaia DR3 4098214367441486592", source_id="4098214367441486592",
         epoch=julian_year(Epoch("2022-08-09T06:34:04")), ra_err=0.2404, dec_err=0.2191,
         mag=Dict("G" => 15.3, "J" => 12.0, "H" => 11.0, "K" => 10.7), diameter=0.44)

"Ephemeris of Quaoar for 2011 or 2022 from the stored Horizons tables."
function quaoar_ephemeris(year::Integer)
    file = year == 2011 ? "quaoar_horizons_20110504.txt" : "quaoar_horizons_20220809.txt"
    rows = [split(l) for l in eachline(joinpath(QUAOAR_DIR, file)) if !startswith(l, "#")]
    return EphemTable([r[1] for r in rows], parse.(Float64, getindex.(rows, 2)), parse.(Float64, getindex.(rows, 3)),
                      parse.(Float64, getindex.(rows, 4)); name="(50000) Quaoar")
end

"""
    quaoar_event(year; rings=true) -> Event

The occultation of 2011 May 4 (`year = 2011`, Braga-Ribas et al. 2013) or of 2022 August 9
(`2022`, Pereira et al. 2023; with `rings = true` the mid-times of the ring detections).
"""
function quaoar_event(year::Integer; rings::Bool=true)
    year in (2011, 2022) || throw(ArgumentError("the events are 2011 and 2022"))
    body = Body("(50000) Quaoar", quaoar_ephemeris(year); radius=550.0)
    if year == 2011
        ev = Event(quaoar_star(2011), body, "2011-05-04T02:42:00"; window=1200, name="Quaoar, 2011-05-04")
        day = "2011-05-04T"
        for (n, lon, lat, alt, sens, ti, σi, te, σe) in QUAOAR_CHORDS_2011
            obs = Observer(n, lon, lat, alt)
            push!(ev, ti === nothing ? Chord(n, obs; span=(day * σi[1], day * σi[2]), sensitivity=sens) :
                      Chord(n, obs; immersion=day * ti, immersion_err=σi, emersion=day * te, emersion_err=σe,
                            sensitivity=sens))
        end
        return ev
    end
    ev = Event(quaoar_star(2022), body, "2022-08-09T06:35:00"; window=1200, name="Quaoar, 2022-08-09")
    day = "2022-08-09T"
    for (n, lon, lat, alt, sens, ti, σi, te, σe, rs) in QUAOAR_CHORDS_2022
        obs = Observer(n, lon, lat, alt)
        if ti === nothing
            push!(ev, Chord(n, obs; span=(day * "06:25:00", day * "06:45:00"), sensitivity=sens))
            continue
        end
        cs = [Contact(:immersion, day * ti, σi), Contact(:emersion, day * te, σe)]
        rings && for (t, σ, ring) in rs
            push!(cs, Contact(:mid, day * t, σ; body=ring))
        end
        ts = sort([c.time for c in cs])
        push!(ev, Chord(n, obs; contacts=cs, span=(ts[1] - 60.0, ts[end] + 60.0), sensitivity=sens))
    end
    return ev
end
