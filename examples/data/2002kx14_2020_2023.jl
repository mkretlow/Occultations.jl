# Five stellar occultations by (119951) 2002 KX14, 2020–2023
#
# Data from J. L. Rizos et al. (2025), "The trans-Neptunian object (119951) 2002 KX14
# revealed via multiple stellar occultations", A&A 697, A62
# (doi:10.1051/0004-6361/202554154, arXiv:2503.19179), with the corrigendum A&A 708, C1
# (2026; doi:10.1051/0004-6361/202659687e): the Occ. D timings of GOA and Sonear3 in Table B.2
# were typographical errors; the corrected values are used (also in arXiv:2503.19179v3). The
# corrigendum states that the results of the paper are not affected.
#   stars: Table 1 (positions at the occultation epochs) and Sect. 2 (angular diameters);
#   stations: Table B.1; timings: Table B.2 (with the corrections of −2.5 s and −0.66 s for
#   Millen-Bruch and Nandrin already applied by the authors); ellipses: Tables 2, B.3, Fig. 1.
# Ephemerides: JPL Horizons (default solution), stored in 2002kx14_horizons.txt; the paper used
# JPL#18, so the fitted centres (offsets from the ephemeris) differ from the published ones.
#
# Notes on the transcription:
# * Table B.1 gives Apache Point Observatory at 32°46'48.3" S; the observatory is at
#   32°46'48" N, which is used here.
# * For Occ. E the paper quotes a stellar diameter of 0.04 mas and 0.65 km at 37.85 au, which
#   are inconsistent (0.04 mas is 1.10 km there); the angular diameter is used.
# * The observed intervals of the negative chords are not tabulated; the paper states that
#   acquisition covered at least ±5 min around the predicted time, which is used.
#
# `kx14_events()` returns the five `Event`s; `KX14_PAPER` the published results.

using Occultations

const KX14_DIR = @__DIR__

"Published results of Rizos et al. (2025)."
const KX14_PAPER = (
    # combined fit of all events (Fig. 1; u, v semi-major and semi-minor axes)
    a=(241.0, 7.2), b=(157.1, 5.2), D_eq=(389.2, 8.7),
    # Occ. A alone (Table 2); position angle in the paper's convention (see the example)
    A=(a=(241.5, 11.0), b=(159.0, 8.1), pa=(105.8, 5.8)),
    # geometric albedo with H_V = 4.978 ± 0.017 (Alvarez-Candal et al. 2016)
    H_V=(4.978, 0.017), albedo=(0.119, 0.007),
    distance_au=Dict("A" => 38.01, "B" => 37.88, "C" => 37.97, "D" => 37.96, "E" => 37.85),
    velocity=Dict("A" => 23.41, "B" => 24.33, "C" => 20.38, "D" => 20.83, "E" => 23.74),
    star_diameter_km=Dict("A" => 1.26, "B" => 0.58, "C" => 1.25, "D" => 0.89, "E" => 0.65),
)

# event, predicted time, Gaia DR3 source, RA, Dec (at the epoch), G, V, K, angular diameter [mas]
const KX14_STARS = [
    ("A", "2020-05-26T22:29:12", "4111560308371475840", "17 24 20.6730", "-23 21 50.203", 14.54, 14.75, 10.52, 0.046),
    ("B", "2022-06-23T00:51:33", "4116601049409385472", "17 35 07.9754", "-23 28 26.010", 16.82, 17.78, 12.31, 0.019),
    ("C", "2023-05-17T08:51:46", "4068746489444765440", "17 45 03.2702", "-23 32 34.722", 15.55, 16.47, 10.71, 0.046),
    ("D", "2023-05-19T03:08:49", "4068793394867887872", "17 44 54.9895", "-23 32 31.595", 16.28, 16.39, 11.36, 0.032),
    ("E", "2023-07-01T07:09:22", "4116448256084153472", "17 41 01.9529", "-23 30 54.310", 16.15, 15.90, 11.90, 0.040),
]

# event, station, longitude (E), latitude, altitude [m], cycle time [s], immersion, σ, emersion, σ
# (times UTC on the date of the event; immersion = nothing for a negative chord)
const KX14_CHORDS = [
    ("A", "Borowiec (Poland)", "17 04 30.8", "52 16 36.8", 83.0, 5.0, "22:30:07.863", 0.697, "22:30:22.307", 0.313),
    ("A", "Odesa-Mayaki (Ukraine)", "30 16 18.1", "46 23 50.2", 23.28, 1.9, "22:29:33.218", 0.736, "22:29:50.632", 0.752),
    ("A", "Konkoly (Hungary)", "19 53 41.7", "47 55 06.0", 477.0, 1.0, "22:30:04.601", 0.258, "22:30:23.218", 0.257),
    ("A", "Gothard (Hungary)", "16 36 15.7", "47 15 29.8", 200.0, 0.7, "22:30:16.151", 0.044, "22:30:33.211", 0.045),
    ("A", "Millen-Bruch (Germany)", "05 53 47.2", "51 01 38.8", 56.0, 0.7, "22:30:42.094", 0.954, "22:30:58.064", 1.110),
    ("A", "Nandrin (Belgium)", "05 26 29.5", "50 31 24.8", 260.9, 1.3, "22:30:45.845", 0.191, "22:30:57.746", 0.185),
    ("A", "Buelach (Switzerland)", "08 34 14.3", "47 31 10.4", 550.0, 1.3, nothing, NaN, nothing, NaN),
    ("A", "Cukurova (Turkey)", "35 21 12.6", "37 03 15.6", 260.9, 10.0, nothing, NaN, nothing, NaN),
    ("B", "Javalambre (Spain)", "-01 00 58.7", "40 02 30.5", 1957.0, 5.4, "00:50:34.107", 1.064, "00:50:45.536", 1.419),
    ("B", "La Hita (Spain)", "-03 10 59.9", "39 34 04.8", 770.0, 12.0, "00:50:36.114", 2.276, "00:50:55.567", 2.930),
    ("B", "Calar Alto (Spain)", "-02 32 44.9", "37 13 24.7", 2173.2, 1.5, "00:50:33.812", 0.391, "00:50:52.265", 0.336),
    ("C", "Shiawassee (USA)", "-84 01 59.2", "43 19 30.8", 178.0, 2.0, "08:50:42.397", 0.344, "08:51:03.586", 0.542),
    ("C", "Comer Field (USA)", "-105 00 35.6", "40 05 29.2", 1567.0, 2.5, "08:52:08.143", 0.652, "08:52:24.456", 1.355),
    ("C", "Nederland (USA)", "-105 26 44.0", "39 59 13.9", 2492.6, 2.5, "08:52:07.448", 0.778, "08:52:25.670", 0.689),
    # Occ. D: corrected timings of the corrigendum (A&A 708, C1, 2026);
    # originally 03:10:38.001 ± 1.623 / 03:10:48.765 ± 0.840 (GOA) and
    # 03:10:49.196 ± 2.328 / 03:11:09.773 ± 1.098 (Sonear3)
    ("D", "GOA (Brazil)", "-40 19 02.0", "-20 18 01.9", 24.0, 6.8, "03:10:38.132", 0.850, "03:10:50.035", 1.345),
    ("D", "Sonear3 (Brazil)", "-43 41 24.0", "-19 49 27.2", 1498.9, 8.0, "03:10:52.356", 0.594, "03:11:08.267", 0.596),
    ("E", "Apache Point (USA)", "-105 49 13.1", "32 46 48.3", 2788.0, 0.5, "07:09:08.865", 0.045, "07:09:22.961", 0.047),
]

"Ephemeris of 2002 KX14 around event `k` (\"A\"…\"E\") from the stored Horizons tables."
function kx14_ephemeris(k::AbstractString)
    rows = [split(l) for l in eachline(joinpath(KX14_DIR, "2002kx14_horizons.txt")) if !startswith(l, "#") && startswith(l, k * " ")]
    return EphemTable([r[2] for r in rows], parse.(Float64, getindex.(rows, 3)), parse.(Float64, getindex.(rows, 4)),
                      parse.(Float64, getindex.(rows, 5)); name="(119951) 2002 KX14")
end

"""
    kx14_events(; negatives=true) -> Vector{Event}

The five occultations A–E (2020-05-26, 2022-06-23, 2023-05-17, 2023-05-19, 2023-07-01) with
the 15 positive chords of Table B.2 and, for event A, the two negative chords.
"""
function kx14_events(; negatives::Bool=true)
    evs = Event[]
    for (k, t, src, ra, dec, G, V, K, θ) in KX14_STARS
        star = Star(ra, dec; name="Gaia DR3 $src", source_id=src, epoch=julian_year(Epoch(t)),
                    mag=Dict("G" => G, "V" => V, "K" => K), diameter=θ)
        body = Body("(119951) 2002 KX14", kx14_ephemeris(k); radius=195.0)
        ev = Event(star, body, t; window=900, name="2002 KX14, Occ. $k, $(t[1:10])")
        day = t[1:11]
        for (e, name, lon, lat, alt, cycle, ti, σi, te, σe) in KX14_CHORDS
            e == k || continue
            obs = Observer(name, lon, lat, alt)
            if ti === nothing
                negatives || continue
                push!(ev, Chord(name, obs; span=(Epoch(t) - 300, Epoch(t) + 300), sensitivity=cycle))
            else
                push!(ev, Chord(name, obs; immersion=day * ti, immersion_err=σi, emersion=day * te, emersion_err=σe,
                                sensitivity=cycle))
            end
        end
        push!(evs, ev)
    end
    return evs
end
