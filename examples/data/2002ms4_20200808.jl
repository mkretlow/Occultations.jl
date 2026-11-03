# Stellar occultation by (307261) 2002 MS4 on 2020 August 8
#
# Data from F. L. Rommel et al. (2023), "A large topographic feature on the surface of the
# trans-Neptunian object (307261) 2002 MS4 measured from stellar occultations", A&A 678, A167
# (doi:10.1051/0004-6361/202346892):
#   star: Table 2; timings (corrected for the time shifts of Table B.1): Table 3; the global
#   ellipse: Table 4; astrometry: Table 5; topography: Sect. 3.2; stations: Tables B.1
#   (positive) and B.2 (negative) at the CDS (J/A+A/678/A167, stored in 2002ms4_cds/; retrieved
#   2026-11-02 from the VizieR mirror at IDIA).
# Ephemeris: JPL Horizons (default solution), stored in 2002ms4_horizons_20200808.txt.
#
# Notes on the transcription:
# * The rows of Table B.1 follow the order of Table 3 (59 positive chords); the station names
#   are those of Table 3, the coordinates those of Table B.1.
# * The sensitivity of a chord is its cycle time (Table B.1, B.2).
# * Observed intervals of the negative chords are not tabulated; 20:38–20:51 UTC is used.
# * The paper's global ellipse uses 13 selected chords (MS4_SELECTED) with a tolerance of
#   7 km for topography.
#
# `ms4_event(; negatives=:all, selected=false)` returns the `Event`; `MS4_PAPER` the published
# results.

using Occultations

const MS4_DIR = @__DIR__

"Published results of Rommel et al. (2023)."
const MS4_PAPER = (
    # Table 4 (3σ errors): centre (f, g) [km] relative to NIMA, a′, apparent oblateness, PA of
    # the minor axis, area-equivalent radius; topographic tolerance 7 km
    ellipse=(f=(43.0, 6.0), g=(7.0, 9.0), a=(412.0, 10.0), oblateness=(0.066, 0.034), position_angle=(121.0, 16.0),
             R_equiv=(398.0, 12.0)), tolerance=7.0,
    # Sect. 3.2: depression of 45.1 ± 1.5 km depth over 322 ± 39 km, centred at PA 40 ± 20°
    # (the prior range); an elevation of 25 (+4 −5) km; a depression of ≈11 km at Varages' egress
    depression=(depth=(45.1, 1.5), extent=(322.0, 39.0)), elevation=25.0,
    # Table 5: geocentric astrometric position at the closest approach
    position=(epoch="2020-08-08T20:44:27.26", ra="18 47 29.961308", dec="-06 16 31.34442"),
    star_diameter_km=1.19, fresnel_km=1.54, distance_au=45.70,
)

"The 13 chords of the paper's global ellipse (Sect. 3.1)."
const MS4_SELECTED = ["Méo station", "Valbonne", "Mátraszentistván", "Catalonia", "Massa 02", "Rome", "Hvar", "Sassari",
                      "Odesa", "Agerola", "Algiers", "La Palma", "Çanakkale"]

# Table 3 (2020 August 8; the order of the rows of Table B.1): station, immersion, σ, emersion, σ
const MS4_TIMINGS = [
    ("Varages", "20:43:44.858", 0.058, "20:43:53.359", 0.028),
    ("TAROT North", "20:43:39.5", 2.0, "20:43:51.5", 2.0),
    ("Méo station", "20:43:39.488", 0.021, "20:43:53.960", 0.027),
    ("Caussols", "20:43:39.27", 0.13, "20:43:54.18", 0.56),
    ("Lleida", "20:44:01.4", 6.6, "20:44:25.13", 0.12),
    ("Cannes", "20:43:37.23", 0.21, "20:43:55.64", 0.39),
    ("Nice", "20:43:36.444", 0.029, "20:43:52.834", 0.071),
    ("Valbonne", "20:43:37.376", 0.055, "20:43:54.07", 0.13),
    ("Črni Vrh", "20:43:09.71", 0.94, "20:43:26.58", 0.35),
    ("Mátraszentistván", "20:42:46.866", 0.015, "20:43:04.739", 0.045),
    ("Ljubljana", "20:43:06.70", 0.20, "20:43:25.20", 0.33),
    ("Budapest", "20:42:49.59", 0.13, "20:43:08.61", 0.20),
    ("Trieste", "20:43:07.7", 3.3, "20:43:27.41", 0.85),
    ("Catalonia", "20:43:57.15", 0.21, "20:44:17.84", 0.14),
    ("Bologna", "20:43:17.40", 0.26, "20:43:38.77", 0.67),
    ("Massa 02", "20:43:20.52", 0.39, "20:43:43.83", 0.32),
    ("Castelvecchio Pascoli", "20:43:18.742", 0.088, "20:43:44.1", 1.4),
    ("Borgo a Mozzano", "20:43:19.224", 0.078, "20:43:42.80", 0.12),
    ("San Marcello Pistoiese", "20:43:16.200", 0.072, "20:43:42.6", 1.2),
    ("Sta. Maria a Monte", "20:43:17.68", 0.12, "20:43:42.80", 0.11),
    ("Signa", "20:43:15.61", 0.57, "20:43:40.38", 0.67),
    ("Lastra a Signa", "20:43:13.815", 0.084, "20:43:41.9", 1.2),
    ("Khmelnytskyi", "20:42:01.81", 0.82, "20:42:49.7", 15.7),
    ("Sevilla", "20:44:29.03", 0.13, "20:44:59.40", 0.19),
    ("Huelva", "20:44:32.74", 0.37, "20:45:02.72", 0.28),
    ("Cluj-Napoca", "20:42:25.0", 1.2, "20:42:56.61", 0.30),
    ("Fuensanta de Martos", "20:44:22.0", 3.0, "20:44:55.0", 3.1),
    ("Fiastra", "20:43:05.1", 1.8, "20:43:37.89", 0.52),
    ("Dragsina", "20:42:32.938", 0.087, "20:43:05.82", 0.32),
    ("Ibiza", "20:43:56.6", 1.8, "20:44:28.55", 0.24),
    ("Alhendín", "20:44:17.89", 0.57, "20:44:54.83", 0.18),
    ("Granada (150 cm)", "20:44:17.571", 0.008, "20:44:52.08", 0.51),
    ("Granada (90 cm)", "20:44:16.53", 0.51, "20:44:51.755", 0.026),
    ("Estepona", "20:44:24.37", 0.45, "20:45:00.4", 1.9),
    ("Marbella", "20:44:26.1", 4.3, "20:44:59.64", 0.22),
    ("Rome", "20:43:08.37", 0.37, "20:43:43.79", 0.37),
    ("Hvar", "20:42:51.69", 0.21, "20:43:27.96", 0.35),
    ("Bacau", "20:42:12.60", 0.30, "20:42:50.49", 0.72),
    ("Sassari", "20:43:19.35", 0.35, "20:43:56.60", 0.50),
    ("Guarcino", "20:43:03.86", 0.15, "20:43:41.13", 0.15),
    ("Kharkiv (70 cm)", "20:41:40.00", 0.13, "20:42:17.16", 0.12),
    ("Kharkiv (36 cm)", "20:41:39.6", 1.9, "20:42:18.3", 1.2),
    ("Ceccano", "20:43:02.39", 0.41, "20:43:40.76", 0.37),
    ("Brasov", "20:42:14.2", 2.7, "20:42:51.5", 1.6),
    ("Bârlad", "20:42:09.23", 0.55, "20:42:47.00", 0.49),
    ("Valenii de Munte", "20:42:15.14", 0.24, "20:42:52.61", 0.50),
    ("Ploiesti", "20:42:14.25", 2.1, "20:42:53.2", 1.5),
    ("Odesa", "20:42:00.00", 0.14, "20:42:38.49", 0.11),
    ("Sorrento", "20:43:00.66", 0.30, "20:43:37.7", 1.1),
    ("Agerola", "20:42:59.603", 0.050, "20:43:37.842", 0.093),
    ("Algiers", "20:43:50.844", 0.021, "20:44:29.115", 0.075),
    ("La Palma", "20:45:35.032", 0.014, "20:46:08.313", 0.018),
    ("Tijarafe", "20:45:35.68", 0.19, "20:46:08.08", 0.11),
    ("Ariana", "20:43:23.87", 0.29, "20:43:56.80", 0.25),
    ("Artemis", "20:45:31.938", 0.022, "20:46:01.297", 0.024),
    ("TAR 1", "20:45:32.360", 0.057, "20:46:01.595", 0.063),
    ("Catania", "20:43:04.59", 0.92, "20:43:34.564", 0.085),
    ("Kuban", "20:41:36.2", 1.2, "20:42:00.63", 0.31),
    ("Çanakkale", "20:42:23.39", 0.56, "20:42:38.830", 0.072),
]

"Rows of a CDS table of stations: (name, latitude, longitude, altitude, cycle time)."
function _ms4_stations(table)
    rows = NamedTuple[]
    for l in eachline(joinpath(MS4_DIR, "2002ms4_cds", table * ".tsv"))
        c = split(l, '\t')
        (startswith(l, "#") || length(c) < 16 || !occursin(r"^\s*\d+$", c[1])) && continue
        cyc = tryparse(Float64, strip(c[16])); exp_ = tryparse(Float64, strip(c[13]))
        push!(rows, (obs=strip(c[2]), city=strip(c[3]), lat=parse(Float64, c[5]), lon=parse(Float64, c[6]),
                     alt=parse(Float64, c[7]), cycle=cyc !== nothing ? cyc : exp_ !== nothing ? exp_ : 0.0))
    end
    return rows
end

"Gaia DR3 4253248324549054464 at the closest approach (Table 2)."
ms4_star() = Star("18 47 29.96384", "-06 16 31.4727"; name="Gaia DR3 4253248324549054464",
                  source_id="4253248324549054464", epoch=julian_year(Epoch("2020-08-08T20:44:27")),
                  ra_err=0.12, dec_err=0.10, mag=Dict("G" => 14.62, "V" => 14.62, "K" => 11.13),
                  diameter=MS4_PAPER.star_diameter_km / (MS4_PAPER.distance_au * 1.495978707e8) * 206264806.247)

"Ephemeris of 2002 MS4 from the stored Horizons table."
function ms4_ephemeris()
    rows = [split(l) for l in eachline(joinpath(MS4_DIR, "2002ms4_horizons_20200808.txt")) if !startswith(l, "#")]
    return EphemTable([r[1] for r in rows], parse.(Float64, getindex.(rows, 2)), parse.(Float64, getindex.(rows, 3)),
                      parse.(Float64, getindex.(rows, 4)); name="(307261) 2002 MS4")
end

"""
    ms4_event(; negatives=:all, selected=false) -> Event

The occultation of 2020 August 8: the 59 positive chords of Table 3 (or, with `selected=true`,
the 13 of the paper's global ellipse) and the 40 negative chords of Table B.2
(`negatives = :all` or `:none`).
"""
function ms4_event(; negatives::Symbol=:all, selected::Bool=false)
    negatives in (:all, :none) || throw(ArgumentError("negatives must be :all or :none"))
    body = Body("(307261) 2002 MS4", ms4_ephemeris(); radius=400.0)
    ev = Event(ms4_star(), body, "2020-08-08T20:44:00"; window=1200, name="2002 MS4, 2020-08-08")
    day = "2020-08-08T"
    pos = _ms4_stations("tableb1")
    length(pos) == length(MS4_TIMINGS) || error("Table B.1 and Table 3 differ in length")
    for (s, (name, ti, σi, te, σe)) in zip(pos, MS4_TIMINGS)
        selected && name ∉ MS4_SELECTED && continue
        push!(ev, Chord(name, Observer(name, s.lon, s.lat, s.alt); immersion=day * ti, immersion_err=σi,
                        emersion=day * te, emersion_err=σe, sensitivity=s.cycle))
    end
    if negatives === :all
        for s in _ms4_stations("tableb2")
            name = s.obs == "-" ? s.city : s.obs * " (" * s.city * ")"
            any(c -> c.name == name, ev.chords) && (name *= " 2")
            push!(ev, Chord(name, Observer(name, s.lon, s.lat, s.alt); span=(day * "20:38:00", day * "20:51:00"),
                            sensitivity=s.cycle))
        end
    end
    return ev
end
