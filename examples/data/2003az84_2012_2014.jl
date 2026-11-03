# Stellar occultations by (208996) 2003 AZ84 on 2012 February 3 and 2014 November 15, and the
# single-chord events of 2011 January 8 and 2013 December 2
#
# Data from A. Dias-Oliveira et al. (2017), "Study of the Plutino Object (208996) 2003 AZ84 from
# Stellar Occultations: Size, Shape, and Topographic Features", AJ 154, 22
# (doi:10.3847/1538-3881/aa74e9):
#   star positions: Eqs. (1) and (2); stations: Tables 2 and 3 and Sect. 3.1.1 (single-chord
#   events); timings: Table 4; limb fits and the Jacobi solution: Table 5; corrected body
#   positions: notes to Table 5.
# Ephemeris: JPL Horizons (default solution), stored in 2003az84_horizons_<yyyymmdd>.txt for the
# four events; the paper refers the fitted centres to NIMA V1.
#
# Notes on the transcription:
# * The star positions are the paper's own (pre-Gaia) astrometry at the event, ± 22 and 23 mas.
# * Tables 2 and 3 give "Longitude (W)" but mix the conventions: the sites in India, Israel,
#   Thailand, China and Japan are listed with negative (i.e. east) values, the European sites
#   with their east longitudes, and the Liverpool Telescope (La Palma, 17° 52′ 45.2″ W) with a
#   positive value. The longitudes below are east longitudes of the real sites: Liverpool
#   −17° 52′ 45.2″, Maidenhead −0° 48′ 58.2″, Alicante "−359° 33′ 18.0″" taken as 0° 26′ 42″ W.
# * IUCAA Girawali (2012): a single exposure with a partial drop between two 14.5 s gaps; the
#   paper represents the ingress (between 19:45:09.88 and 19:45:24.43) by its mid-time ± 7.3 s.
# * Moriyama (2014): the CCD and the video timings agree; the paper fits only the CCD chord,
#   the video chord is kept here but disabled.
# * Negative chords: observed intervals are not tabulated; ± 10 min around the event is used.
#   Sites that stopped before the event (Inabe, Tarui, Maibara) or observed the wrong field
#   (Hitachi) and the clouded Wise Observatory are omitted.
# * The single-chord events of 2011 January 8 (San Pedro de Atacama) and 2013 December 2
#   (Rockhampton): the paper does not identify the stars. Each is the only Gaia DR3 source
#   within 5″ of the body's position at the event (Horizons; retrieved 2026-11-02), 0.1″ from
#   it: 3151679712277335168 (G = 18.4) and 3147373651082806272 (G = 14.9); used with their
#   proper motions. The paper's chord lengths are 573 ± 21 km and 657 ± 35 km. The exposure
#   times of 2011 are not given.
#
# `az84_events()` returns the two `Event`s; `AZ84_PAPER` the published results.

using Occultations

const AZ84_DIR = @__DIR__

"Published results of Dias-Oliveira et al. (2017)."
const AZ84_PAPER = (
    # Table 5 (asymmetric errors as (value, minus, plus); centres relative to NIMA V1; PL = position
    # angle of the minor axis, N → E)
    e2012=(a=(426.0, 16.0, 26.0), R_equiv=(346.0, 6.0, 6.0), oblateness=(0.340, 0.086, 0.097),
           f=(6800.0, 17.0, 23.0), g=(2850.0, 11.0, 14.0), position_angle=(41.0, 5.0, 8.0), chi2_pdf=0.73,
           # corrected geocentric position at 2012-02-03 19:45 UTC (notes to Table 5)
           position=(epoch="2012-02-03T19:45:00", ra="07 45 54.80510", dec="+11 12 43.1310")),
    e2014=(a=(393.0, 1.0, 1.0), R_equiv=(382.0, 3.0, 3.0), oblateness=(0.054, 0.003, 0.003),
           f=(2862.0, 1.8, 2.1), g=(336.0, 2.7, 4.0), position_angle=(36.0, 4.0, 4.0), chi2_pdf=0.98,
           position=(epoch="2014-11-15T17:55:00", ra="08 03 51.29112", dec="+09 57 18.8060")),
    # Jacobi solution (Table 5), rotation (Sect. 4.2), single-chord lower limits (Sect. 3.1.1)
    jacobi=(a=(470.0, 20.0), b_over_a=(0.82, 0.05), c_over_a=(0.52, 0.02), density=(0.87, 0.01)),
    albedo=(0.097, 0.009), H_V=(3.74, 0.09), period_h=6.75, amplitude=0.07,
    single_chords_km=((573.0, 21.0), (657.0, 35.0)),
    velocity=(25.0, 9.4), star_diameter_km=(0.24, 0.20),
)

# station, longitude (E), latitude, altitude [m], cycle time [s],
# immersion, σ, emersion, σ (UTC; `nothing` for a negative chord)
const AZ84_CHORDS_2012 = [
    ("Mt. Abu", "72 46 47.18", "24 39 10.34", 1680.0, 4.98, "19:45:13.89", 0.70, "19:45:35.95", 2.5),
    ("IUCAA Girawali", "73 40 00.0", "19 05 00.0", 1005.0, 16.5, "19:45:17.16", 7.3, "19:45:25.43", 0.10),
    ("Kraar", "34 48 45.76", "31 54 29.1", 107.0, 5.0, "19:47:39.20", 0.70, "19:48:06.15", 0.68),
    ("Liverpool Tel.", "-17 52 45.2", "28 45 48.8", 2457.0, 1.5, nothing, NaN, nothing, NaN),
    ("Nikaya", "77 43 33.0", "12 36 23.0", 1005.0, 6.2, nothing, NaN, nothing, NaN),
    ("Maidenhead", "-00 48 58.2", "51 30 23.8", 75.0, 2.56, nothing, NaN, nothing, NaN),
    ("Paris Obs. mobile", "01 50 19.9", "46 31 08.2", 305.0, 2.56, nothing, NaN, nothing, NaN),
    ("Cabrils", "02 23 07.3", "41 31 11.3", 114.0, 13.0, nothing, NaN, nothing, NaN),
    ("Alicante", "-00 26 42.0", "38 28 33.0", 187.3, 5.0, nothing, NaN, nothing, NaN),
    ("Tourrette Levens", "07 15 47.2", "43 47 22.2", 385.0, 5.68, nothing, NaN, nothing, NaN),
    ("Bellinzona", "09 01 26.5", "46 13 53.2", 260.0, 7.0, nothing, NaN, nothing, NaN),
    ("Berlin (Wortmann)", "13 28 30.8", "52 29 12.5", 41.0, 2.56, nothing, NaN, nothing, NaN),
    ("Berlin (Andersson)", "13 39 41.5", "52 24 35.4", 45.0, 1.28, nothing, NaN, nothing, NaN),
    ("Nonndorf", "15 14 08.7", "48 47 13.5", 593.0, 2.56, nothing, NaN, nothing, NaN),
]

# the single-chord events (Sect. 3.1.1, Table 4): station, longitude (E), latitude, altitude [m],
# cycle time [s] (0: not given), date, reference time, immersion, σ, emersion, σ (UTC)
const AZ84_SINGLE = Dict(
    2011 => ("San Pedro de Atacama (SPACE)", "-68 10 48.70", "-22 57 09.80", 2400.0, 0.0, "2011-01-08T", "06:31:00",
             "06:29:48.30", 0.80, "06:30:10.00", 0.8),
    2013 => ("Rockhampton", "150 30 00.80", "-23 16 09.60", 50.0, 2.56, "2013-12-02T", "14:54:00",
             "14:52:48.44", 1.40, "14:53:28.17", 1.0),
)

const AZ84_CHORDS_2014 = [
    ("Thai National Telescope", "98 25 56.06", "18 34 25.41", 2457.0, 0.801169, "17:58:09.57", 0.02, "17:59:13.80", 0.035),
    ("Yunnan", "100 01 51.0", "26 42 32.0", 3193.0, 4.3164, "17:57:21.65", 0.13, "17:57:36.95", 0.135),
    ("Moriyama (CCD)", "135 59 23.8", "35 02 59.1", 105.0, 11.65, "17:52:15.32", 0.18, "17:52:59.90", 0.20),
    ("Moriyama (video)", "135 59 23.8", "35 02 59.1", 105.0, 8.53, "17:52:15.25", 3.6, "17:52:55.69", 4.2),
    ("Hamamatsu", "137 44 23.0", "34 43 07.0", 17.0, 8.53, nothing, NaN, nothing, NaN),
]

"""
    az84_star(year)

The star of an event: the paper's astrometry at the event (Eqs. 1 and 2) for 2012 and 2014,
Gaia DR3 for the single-chord events of 2011 and 2013.
"""
az84_star(year::Integer) = year == 2011 ?
    Star(115.9242759644801, 11.506543527957033; name="Gaia DR3 3151679712277335168", source_id="3151679712277335168",
         epoch=2016.0, pmra=-1.7997793993138973, pmdec=-1.5984502818738122, parallax=0.2305421561939682,
         ra_err=0.1506, dec_err=0.1202, pmra_err=0.1956, pmdec_err=0.1401, parallax_err=0.1898,
         mag=Dict("G" => 18.427, "BP" => 18.826, "RP" => 17.842)) :
                           year == 2013 ?
    Star(119.73636844201489, 10.318054123276468; name="Gaia DR3 3147373651082806272", source_id="3147373651082806272",
         epoch=2016.0, pmra=-2.4652673834774324, pmdec=6.844528908724911, parallax=0.5058291450257313,
         ra_err=0.2654, dec_err=0.1831, pmra_err=0.3298, pmdec_err=0.2011, parallax_err=0.2941,
         mag=Dict("G" => 14.883, "BP" => 15.164, "RP" => 14.335)) :
                           year == 2012 ?
    Star("07 45 54.7696", "+11 12 43.093"; name="2003 AZ84 star 2012", epoch=julian_year(Epoch("2012-02-03T19:45")),
         ra_err=22.0, dec_err=23.0, mag=Dict("B" => 15.8, "V" => 15.35, "K" => 14.1),
         diameter=0.24 / 6.6209e9 * 206264806.247) :
    Star("08 03 51.2980", "+09 57 18.729"; name="2003 AZ84 star 2014", epoch=julian_year(Epoch("2014-11-15T17:55")),
         ra_err=22.0, dec_err=23.0, mag=Dict("B" => 15.5, "V" => 15.5, "K" => 14.0),
         diameter=0.20 / 6.6634e9 * 206264806.247)

"Ephemeris of 2003 AZ84 for an event (2011, 2012, 2013 or 2014) from the stored Horizons tables."
function az84_ephemeris(year::Integer)
    file = "2003az84_horizons_" * Dict(2011 => "20110108", 2012 => "20120203", 2013 => "20131202",
                                       2014 => "20141115")[year] * ".txt"
    rows = [split(l) for l in eachline(joinpath(AZ84_DIR, file)) if !startswith(l, "#")]
    return EphemTable([r[1] for r in rows], parse.(Float64, getindex.(rows, 2)), parse.(Float64, getindex.(rows, 3)),
                      parse.(Float64, getindex.(rows, 4)); name="(208996) 2003 AZ84")
end

"""
    az84_event(year; negatives=:all) -> Event

The occultation of 2012 February 3 (`year = 2012`) or of 2014 November 15 (`2014`) with its
positive and negative chords (Tables 2–4; `negatives = :all` or `:none`). The Moriyama video
chord (2014) is included but disabled, as in the paper. `year = 2011` or `2013` gives the
single-chord event of 2011 January 8 or 2013 December 2.
"""
function az84_event(year::Integer; negatives::Symbol=:all)
    negatives in (:all, :none) || throw(ArgumentError("negatives must be :all or :none"))
    year in (2011, 2012, 2013, 2014) || throw(ArgumentError("the events are 2011, 2012, 2013 and 2014"))
    body = Body("(208996) 2003 AZ84", az84_ephemeris(year); radius=380.0)
    if haskey(AZ84_SINGLE, year)
        name, lon, lat, alt, sens, day, tc, ti, σi, te, σe = AZ84_SINGLE[year]
        ev = Event(az84_star(year), body, day * tc; window=1200, name="2003 AZ84, " * day[1:10])
        push!(ev, Chord(name, Observer(name, lon, lat, alt); immersion=day * ti, immersion_err=σi, emersion=day * te,
                        emersion_err=σe, sensitivity=sens))
        return ev
    end
    t, day, span = year == 2012 ? ("2012-02-03T19:46:00", "2012-02-03T", ("19:36:00", "19:58:00")) :
                                  ("2014-11-15T17:55:00", "2014-11-15T", ("17:45:00", "18:09:00"))
    ev = Event(az84_star(year), body, t; window=1200, name="2003 AZ84, " * day[1:10])
    for (name, lon, lat, alt, sens, ti, σi, te, σe) in (year == 2012 ? AZ84_CHORDS_2012 : AZ84_CHORDS_2014)
        obs = Observer(name, lon, lat, alt)
        if ti === nothing
            negatives === :all || continue
            push!(ev, Chord(name, obs; span=(day * span[1], day * span[2]), sensitivity=sens))
        else
            c = Chord(name, obs; immersion=day * ti, immersion_err=σi, emersion=day * te, emersion_err=σe,
                      sensitivity=sens)
            push!(ev, name == "Moriyama (video)" ? disable(c) : c)
        end
    end
    return ev
end
