# Stellar occultation by the dwarf planet (136108) Haumea and its ring on 2017 January 21
#
# Data from J. L. Ortiz et al. (2017), "The size, shape, density and ring of the dwarf planet
# Haumea from a stellar occultation", Nature 550, 219–223 (doi:10.1038/nature24051):
#   star: Extended Data Table 1; stations and the timings of the body: Table 1; the ring
#   events: Extended Data Table 2; limb and ring fits, 3-D shape: main text and Methods.
# Ephemeris: JPL Horizons (default solution), stored in haumea_horizons_20170121.txt (the event,
# and ten minutes of 2009 July 1 for the direction of the 2009 light-curve amplitude; the date of
# the HST observations within 2009 is not given, mid-2009 is assumed).
#
# Notes on the transcription:
# * The star position (URAT1 533-182543, Gaia DR1 1233009038221203584) is the paper's own at the
#   epoch 2017.0575 (the event).
# * Konkoly and Wendelstein observed with two telescopes each: 12 chords from 10 sites.
# * Crni Vrh: 300 s drifted exposures, ± 8 and ± 10 s.
# * Mount Agliale did not see the body, only the ring: a chord with ring detections only.
# * The ring events are given by their mid-times (Extended Data Table 2). The egress at
#   Skalnate Pleso fell in a readout gap and "is not a real detection of the ring": it is kept
#   but disabled.
# * Observed intervals are not tabulated: the light curves (Fig. 1) cover 03:06–03:12:40 UTC;
#   03:05–03:13 is used.
# * The sensitivity of a chord is its cycle time.
#
# `haumea_event()` returns the `Event`; `HAUMEA_PAPER` the published results.

using Occultations

const HAUMEA_DIR = @__DIR__

"Published results of Ortiz et al. (2017)."
const HAUMEA_PAPER = (
    # the limb: semi-axes a′ = 852 ± 2 km, b′ = 569 ± 13 km, PA of the minor axis −76.3 ± 1.2°
    a=(852.0, 2.0), b=(569.0, 13.0), position_angle=(-76.3, 1.2),
    # the ring: apparent semi-axes 2287 +75 −45 km and 541 ± 15 km, PA of the minor axis
    # −74.3 ± 1.3°, opening angle 13.8 ± 0.5° (13.7° in Fig. 2), two poles (solution 1
    # preferred, close to the orbit pole of Hi'iaka (283.0°, −10.6°)); width 70 km, apparent
    # opacity 0.5; radial width at Konkoly 70 km, equivalent widths 41 and 25 km
    ring_a=(2287.0, 45.0, 75.0), ring_b=(541.0, 15.0), ring_pa=(-74.3, 1.3), ring_B=(13.8, 0.5),
    ring_pole1=((285.1, 0.5), (-10.6, 1.2)), ring_pole2=((312.3, 0.3), (-18.6, 1.2)),
    ring_width=70.0, ring_opacity=0.5, hiiaka_pole=(283.0, -10.6),
    # 3-D shape (with the HST amplitude of 0.32 mag in 2009 and the rotational phase: Haumea at
    # its absolute brightness minimum), density, albedo
    axes_3d=((1161.0, 30.0), (852.0, 4.0), (513.0, 16.0)), D_eq_volume=(1595.0, 11.0),
    density=(1885.0, 80.0), albedo=(0.51, 0.02), period_h=3.915341, amplitude_2009=0.32,
    velocity=13.1, star_diameter_mas=0.007,
)

# station, longitude (E), latitude, altitude [m], cycle time [s],
# immersion, σ, emersion, σ of the body (UTC, 2017-01-21; `nothing`: not seen), ring events
# (kind, mid-time, σ, used)
const HAUMEA_CHORDS = [
    ("Skalnate Pleso 1.3m", "20 14 02.1", "49 11 21.8", 1826.0, 15.5, "03:08:26.79", 0.96, "03:10:24.56", 0.8,
     [("03:07:04.1", 5.6, true), ("03:11:50", 2.4, false)]),
    ("Konkoly 1.0m", "19 53 41.5", "47 55 01.6", 935.0, 1.007, "03:08:20.3", 0.2, "03:10:17.39", 0.07,
     [("03:06:58.1", 0.3, true), ("03:11:40.7", 0.2, true)]),
    ("Konkoly 0.6m", "19 53 41.5", "47 55 01.6", 935.0, 2.944, "03:08:19.5", 0.8, "03:10:16.4", 1.3, []),
    ("Ondrejov", "14 46 53.3", "49 54 32.6", 526.0, 9.721, "03:08:29.2", 0.8, "03:10:12.2", 0.8,
     [("03:06:48.4", 4.2, true), ("03:10:51.6", 3.0, true)]),
    ("Crni Vrh", "14 04 15.9", "45 56 45.0", 713.0, 315.0, "03:07:54", 8.0, "03:09:57", 10.0, []),
    ("Wendelstein 2.0m", "12 00 44.0", "47 42 13.6", 1838.0, 14.536, "03:08:27.9", 2.8, "03:09:34.1", 0.5,
     [("03:06:39.5", 3.9, true)]),
    ("Wendelstein 0.4m", "12 00 44.0", "47 42 13.6", 1838.0, 53.096, "03:08:18.8", 6.0, "03:09:38.6", 6.0, []),
    ("Bavarian Public Obs.", "11 36 25.2", "48 07 19.2", 538.0, 20.304, "03:08:30.0", 3.3, "03:09:30.0", 4.9, []),
    ("Asiago", "11 34 08.4", "45 50 54.9", 1376.0, 5.026, "03:08:20.17", 0.08, "03:09:13.27", 1.5,
     [("03:06:35.4", 0.3, true)]),
    ("S. Marcello Pistoiese", "10 48 14.0", "44 03 51.0", 965.0, 11.877, "03:08:22.9", 0.9, "03:08:42.8", 0.9,
     [("03:06:38.8", 1.5, true)]),
    ("Lajatico", "10 43 01.2", "43 25 44.7", 433.0, 16.254, "03:08:19.9", 1.4, "03:08:34.3", 1.4,
     [("03:06:37.9", 5.1, true)]),
    ("Mount Agliale", "10 30 53.8", "43 59 43.1", 758.0, 16.724, nothing, NaN, nothing, NaN,
     [("03:06:56.7", 1.4, true), ("03:08:15.3", 4.2, true)]),
]

"URAT1 533-182543 (Gaia DR1 1233009038221203584) at the epoch of the event (Extended Data Table 1)."
haumea_star() = Star("14 12 03.2034", "+16 33 58.642"; name="URAT1 533-182543", source_id="1233009038221203584",
                     epoch=2017.0575, mag=Dict("B" => 20.88, "V" => 17.97, "R" => 18.0, "K" => 14.532),
                     diameter=HAUMEA_PAPER.star_diameter_mas)

"Ephemeris of Haumea from the stored Horizons table."
function haumea_ephemeris()
    rows = [split(l) for l in eachline(joinpath(HAUMEA_DIR, "haumea_horizons_20170121.txt")) if !startswith(l, "#")]
    return EphemTable([r[1] for r in rows], parse.(Float64, getindex.(rows, 2)), parse.(Float64, getindex.(rows, 3)),
                      parse.(Float64, getindex.(rows, 4)); name="(136108) Haumea")
end

"""
    haumea_event(; rings=true) -> Event

The occultation of 2017 January 21: twelve chords of the body (Table 1) and, with
`rings = true`, the mid-times of the ring events (Extended Data Table 2, labelled `:ring1`;
the egress at Skalnate Pleso disabled), including Mount Agliale, which saw only the ring.
"""
function haumea_event(; rings::Bool=true)
    body = Body("(136108) Haumea", haumea_ephemeris(); radius=800.0)
    ev = Event(haumea_star(), body, "2017-01-21T03:09:00"; window=1200, name="Haumea, 2017-01-21")
    day = "2017-01-21T"
    for (name, lon, lat, alt, sens, ti, σi, te, σe, ring) in HAUMEA_CHORDS
        obs = Observer(name, lon, lat, alt)
        cs = Contact[]
        ti === nothing || push!(cs, Contact(:immersion, day * ti, σi), Contact(:emersion, day * te, σe))
        rings && for (t, σ, use) in ring
            push!(cs, Contact(:mid, day * t, σ; body=:ring1, use=use))
        end
        push!(ev, isempty(cs) ? Chord(name, obs; span=(day * "03:05:00", day * "03:13:00"), sensitivity=sens) :
                  Chord(name, obs; contacts=cs, span=(day * "03:05:00", day * "03:13:00"), sensitivity=sens))
    end
    return ev
end
