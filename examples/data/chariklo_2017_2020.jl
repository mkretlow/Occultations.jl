# Stellar occultations by the centaur (10199) Chariklo and its rings, 2017–2020
#
# Data from B. E. Morgado et al. (2021), "Refined physical parameters for Chariklo's body and
# rings from stellar occultations observed between 2013 and 2020", A&A 652, A141
# (doi:10.1051/0004-6361/202141543):
#   stars: Table 1; ring pole per event: Table 2; astrometry: Table 5; summary: Table 6;
#   stations: Table A.1; timings of the main body and the ring detections: Table C.1 (CDS,
#   J/A+A/652/A141, stored in chariklo_cds/tablec1.tsv; retrieved 2026-11-02 from the VizieR
#   mirror at IDIA).
# Ephemeris: JPL Horizons (default solution), stored in chariklo_horizons_<yyyymmdd>.txt.
#
# Notes on the transcription:
# * The stars are the Gaia EDR3 positions propagated to the event (Table 1).
# * Table C.1 gives, for each detection of the main body (MB), the immersion and emersion
#   times; for each ring detection (C1R, C2R; "1st" before and "2nd" after the body) the start
#   and end of the dip. A ring detection is used here by its mid-time (σ: half the quadratic
#   sum of the two errors); the dips of C1R are 0.1–1.5 s long.
# * The two grazing ring detections of 2017 July 23 (SARA, PROMPT; one long dip each) are left
#   out, as in the paper's fit (Table 2, note b).
# * The two SSO telescopes at Paranal (Io: first, Europa: second detection) are one chord here.
# * Stations without detection (negative chords, overcast) are not included: their observed
#   intervals are not given.
# * The 2019 September 4 event (Mauna Kea, Mauna Loa: instrumental issues) has no entry in
#   Table C.1.
#
# `chariklo_event(date)` returns the `Event` of a date ("2017-04-09", …); `CHARIKLO_PAPER`
# the published results.

using Occultations

const CHARIKLO_DIR = @__DIR__

"Published results of Morgado et al. (2021)."
const CHARIKLO_PAPER = (
    # Table 2: ring C1R radius [km], position angle P and opening angle B of the ring [deg],
    # pole (α, δ) [deg] from the events with multi-chord ring detections
    ring_fits=Dict("2017-06-22" => (r=(386.6, 12.3), P=(-62.13, 3.46), B=(52.33, 1.97), pole=((149.46, 2.60), (40.98, 2.12))),
                   "2017-07-23" => (r=(385.9, 0.4), P=(-60.85, 0.11), B=(49.91, 0.11), pole=((151.03, 0.14), (41.81, 0.07)))),
    # Table 6: rings and the 3-D shape of the body
    C1R=(radius=(385.9, 0.4), pole=((151.03, 0.14), (41.81, 0.07)), width=6.86),
    C2R=(radius=(399.8, 0.6), pole=((150.91, 0.22), (41.60, 0.12)), equivalent_width=(0.117, 0.080)),
    gap=13.9, axes_3d=((143.8, 1.5), (135.2, 2.8), (99.1, 2.7)), R_vol=(124.8, 2.7),
    # Table 5: astrometric positions (with the stars of Table 1)
    positions=Dict("2020-06-19" => ("2020-06-19T15:51:00", "20 02 32.0002322", "-22 20 41.488487"),
                   "2019-08-08" => ("2019-08-08T21:41:00", "19 31 49.7910264", "-25 34 26.785014"),
                   "2019-08-02" => ("2019-08-02T10:02:00", "19 33 07.4812764", "-25 34 42.519146"),
                   "2017-08-24" => ("2017-08-24T02:59:00", "18 42 35.2371826", "-31 09 50.561462"),
                   "2017-07-23" => ("2017-07-23T05:58:00", "18 48 09.2288885", "-31 26 32.437598"),
                   "2017-06-22" => ("2017-06-22T21:18:00", "18 55 15.6602082", "-31 31 21.621802"),
                   "2017-04-09" => ("2017-04-09T02:24:00", "19 04 03.6255610", "-31 17 15.257638")),
)

# Table 1: Gaia EDR3 source and position at the event
const CHARIKLO_STARS = Dict(
    "2020-06-19" => ("6852214948874014720", "20 02 31.99612", "-22 20 41.4784", 15.569),
    "2019-08-08" => ("6766365630734160896", "19 31 49.78464", "-25 34 26.7724", 15.158),
    "2019-08-02" => ("6766356933423996544", "19 33 07.47470", "-25 34 42.3668", 17.766),
    "2017-08-24" => ("6737213175138988672", "18 42 35.22729", "-31 09 50.6837", 17.468),
    "2017-07-23" => ("6737020112089260672", "18 48 09.22138", "-31 26 32.4591", 14.005),
    "2017-06-22" => ("6760223758801661440", "18 55 15.65251", "-31 31 21.6706", 14.223),
    "2017-04-09" => ("6757456459825426560", "19 04 03.62801", "-31 17 15.2304", 14.094),
)

# Table A.1: station → longitude (E), latitude, altitude [m] (names as in Table C.1)
const CHARIKLO_SITES = Dict(
    "Glenlee" => ("150 30 01.6", "-23 16 10.1", 50.0),
    "Ste. Marie - Mondon" => ("55 34 00.2", "-20 53 48.5", 50.0),
    "Ste. Marie - Hesler" => ("55 31 08.6", "-20 54 38.9", 85.0),
    "Maido" => ("55 23 14.5", "-21 04 15.5", 2195.0),
    "Les Makes" => ("55 24 36.3", "-21 11 56.1", 990.0),
    "Langevin" => ("55 38 47.7", "-21 23 12.0", 20.0),
    "Yass" => ("148 58 35.1", "-34 51 50.9", 536.0),
    "Murrumbateman" => ("148 59 54.8", "-34 57 31.5", 594.0),
    "OPD" => ("-45 34 57.5", "-22 32 07.8", 1864.0),
    "S. P. A. - Fabrega" => ("-68 10 48.0", "-22 57 09.0", 2396.0),
    "VLT" => ("-70 24 16.5", "-24 37 39.3", 2635.0),
    "SSO" => ("-70 24 16.5", "-24 37 39.3", 2635.0),
    "Tolar Grande" => ("-67 23 43.6", "-24 35 23.5", 3524.0),
    "Cono de Arita" => ("-67 46 00.2", "-25 01 59.8", 3532.0),
    "El Rodeo" => ("-66 08 35.4", "-24 52 37.0", 3022.0),
    "Cachi Adentro" => ("-66 14 16.7", "-25 05 33.9", 2684.0),
    "Augusto de Lima" => ("-44 04 08.4", "-18 01 20.7", 620.0),
    "El Salvador" => ("-69 45 03.8", "-26 18 48.3", 1607.0),
    "Sao J. Rio Preto" => ("-49 17 52.2", "-20 42 43.4", 482.0),
    "Inca de Oro" => ("-69 54 45.6", "-26 44 45.3", 1592.0),
    "Bilac" => ("-50 30 20.9", "-21 23 38.6", 426.0),
    "Danish Visual" => ("-70 44 20.2", "-29 15 21.3", 2337.0),
    "Danish Red" => ("-70 44 20.2", "-29 15 21.3", 2337.0),
    "1m La Silla" => ("-70 44 20.2", "-29 15 21.3", 2337.0),
    "TRAPPIST-South" => ("-70 44 20.2", "-29 15 21.3", 2337.0),
    "Outeniqua Lodge" => ("16 49 17.7", "-21 17 58.2", 1416.0),
    "Onduruquea G. F." => ("15 59 33.7", "-21 36 26.0", 1220.0),
    "Windhoek - Backes" => ("17 06 31.9", "-22 41 55.2", 1902.0),
    "Windhoek - Meza" => ("17 06 31.9", "-22 41 55.2", 1902.0),
    "Tivoli Lodge" => ("18 01 01.2", "-23 27 40.2", 1344.0),
    "Hakos 2017-06-22" => ("16 21 41.3", "-23 14 11.0", 1843.0),
    "Wabi Lodge" => ("17 32 00.6", "-20 20 38.9", 1382.0),
    "Weaver R. L." => ("16 50 29.9", "-20 41 56.3", 1647.0),
    "Hakos 2017-04-09" => ("16 21 42.0", "-23 14 10.0", 1881.0),
)

"The rows of Table C.1: (date, station, detection, immersion, σ, emersion, σ)."
function chariklo_table()
    rows = NamedTuple[]
    for l in eachline(joinpath(CHARIKLO_DIR, "chariklo_cds", "tablec1.tsv"))
        c = split(l, '\t')
        (length(c) < 8 || startswith(l, "#") || !occursin(r"^\s*\d+$", c[1])) && continue
        push!(rows, (date=strip(c[3]), station=strip(c[2]), det=strip(c[4]), imm=strip(c[5]),
                     e_imm=parse(Float64, c[6]), em=strip(c[7]), e_em=parse(Float64, c[8])))
    end
    return rows
end

"The Gaia EDR3 star of the event of `date`, at the event (Table 1)."
function chariklo_star(date::AbstractString)
    id, ra, dec, G = CHARIKLO_STARS[date]
    return Star(ra, dec; name="Gaia EDR3 " * id, source_id=id, epoch=julian_year(Epoch(CHARIKLO_PAPER.positions[date][1])),
                mag=Dict("G" => G))
end

"Ephemeris of Chariklo for the event of `date` from the stored Horizons table."
function chariklo_ephemeris(date::AbstractString)
    file = joinpath(CHARIKLO_DIR, "chariklo_horizons_" * replace(date, "-" => "") * ".txt")
    rows = [split(l) for l in eachline(file) if !startswith(l, "#")]
    return EphemTable([r[1] for r in rows], parse.(Float64, getindex.(rows, 2)), parse.(Float64, getindex.(rows, 3)),
                      parse.(Float64, getindex.(rows, 4)); name="(10199) Chariklo")
end

"""
    chariklo_event(date; rings=true) -> Event

The occultation of `date` ("2017-04-09", "2017-06-22", "2017-07-23", "2017-08-24",
"2019-08-02", "2019-08-08", "2020-06-19") with the chords of Table C.1: immersion and emersion
of the main body, and (`rings = true`) the mid-times of the ring detections, labelled `:ring1`
(C1R) and `:ring2` (C2R).
"""
function chariklo_event(date::AbstractString; rings::Bool=true)
    rows = filter(r -> r.date == date, chariklo_table())
    isempty(rows) && throw(ArgumentError("no data for $date"))
    body = Body("(10199) Chariklo", chariklo_ephemeris(date); radius=125.0)
    ev = Event(chariklo_star(date), body, date * "T" * Dict("2020-06-19" => "15:52:00", "2019-08-08" => "21:39:00",
                                                             "2019-08-02" => "10:05:00", "2017-08-24" => "02:57:00",
                                                             "2017-07-23" => "05:56:00", "2017-06-22" => "21:21:00",
                                                             "2017-04-09" => "02:17:00")[date];
               window=1200, name="Chariklo, " * date)
    day = date * "T"
    for st in unique(r.station for r in rows)
        rs = filter(r -> r.station == st, rows)
        any(r -> r.det == "C1R", rs) && continue                 # the grazing detections (SARA, PROMPT)
        key = haskey(CHARIKLO_SITES, st) ? st : st * " " * date
        lon, lat, alt = CHARIKLO_SITES[key]
        cs = Contact[]
        for r in rs
            if r.det == "MB"
                push!(cs, Contact(:immersion, day * r.imm, r.e_imm), Contact(:emersion, day * r.em, r.e_em))
            elseif rings
                tm = (Epoch(day * r.imm) + (Epoch(day * r.em) - Epoch(day * r.imm)) / 2)
                push!(cs, Contact(:mid, tm, hypot(r.e_imm, r.e_em) / 2; body=startswith(r.det, "C2R") ? :ring2 : :ring1))
            end
        end
        isempty(cs) && continue
        # observed interval: the detections ± 30 s (a chord with ring detections only is a
        # negative chord for the body)
        ts = sort([c.time for c in cs])
        push!(ev, Chord(st, Observer(st, lon, lat, alt); contacts=cs, span=(ts[1] - 30.0, ts[end] + 30.0)))
    end
    return ev
end
