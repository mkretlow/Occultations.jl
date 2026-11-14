# Stellar occultation by Triton on 2017 October 5: a central flash
#
# Data from J. Marques Oliveira et al. (2022), "Constraints on the structure and seasonal
# variations of Triton's atmosphere from the 5 October 2017 stellar occultation and previous
# observations", A&A 659, A136 (doi:10.1051/0004-6361/202141443): the star and the geometry
# (Tables 1 and 2), the stations with their exposure and cycle times (Table A.1), the temperature
# template (Table B.1, `triton_atmosphere()`) and the result (Table 4). The 90 light curves are at
# the CDS (J/A+A/659/A136, subdirectory lc/); they are to be stored in triton_cds/ (not yet
# downloaded: the CDS file server was unreachable on 2026-11-14).
# Ephemeris: JPL Horizons (Triton satellite solution nep098), stored in
# triton_horizons_20171005.txt; the paper refers its offsets to DE435/NEP081.
#
# Notes on the transcription:
# * Coordinates as in Table A.1 (WGS84 assumed), altitudes in m. Calar Alto (VIS, SWIR), Kryoneri
#   (R, I), Javalambre (Tx40, Excalibur), Saint Caprais (0.20 m, 0.94 m), La Hita (0.77 m, 0.40 m)
#   and Pic du Midi (0.60 m, 1.05 m) are separate light curves with the same coordinates.
# * La Sagra (four telescopes, one light curve): exposure and cycle of the light curve unknown
#   (tetra1/2: 1.00 s every 3.79/3.56 s; tetra3/4: 0.20 s); Massa: 5.12 s in the first half,
#   2.56 s in the second half of the light curve.
# * The effective wavelength λ and bandwidth Δλ [µm] are estimated from the detector and filter
#   (not given in the paper): unfiltered CCD/CMOS/video 0.65/0.4, R 0.65/0.15, V 0.55/0.1, V+R
#   0.6/0.25, I 0.8/0.15, I+z 0.85/0.25, visible EMCCD 0.6/0.4, SWIR 1.35/0.7, InGaAs behind
#   RG1000 1.3/0.6 (Table A.1: effective λ = 1.3 µm).
#
# `TRITON_STATIONS` lists all 90 light curves; `used = true` marks the 52 of the global fit.
# `TRITON_PAPER` holds the published results.

using Occultations

const TRITON_DIR = @__DIR__
const TRITON_CA = "2017-10-05T23:51:28.92"   # retrieved geocentric closest approach (Table 1)

"Published results of Marques Oliveira et al. (2022)."
const TRITON_PAPER = (
    p1400=(1.18, 0.03),            # µbar, pressure at r = 1400 km (Table 4), χ²/dof = 0.85
    p_surf_ratio=12.0,             # p_surf = 12.0 p1400 (Sect. 4.2)
    Δρ=(-359.3, 1.0),              # km, cross-track offset from DE435/NEP081 (Sect. 4.2)
    ca=(23 * 3600 + 51 * 60 + 28.92, 0.02),     # s after 0h UTC, geocentric closest approach (Table 1)
    ca_mas=(196.58, 0.05),         # geocentric closest approach distance Triton–star [mas]
    radius=1353.0, GM=1428.0,      # km, km³/s² (Table 2)
    distance=4.3506e9,             # km
    pole_pa=305.7,                 # deg, position angle of the north pole (Table 2)
    star_diam_km=0.65,             # stellar diameter projected at Triton (Table 1)
    φ0_lapalma=0.35885,            # Triton's contribution at La Palma adopted for the template (App. B.2)
    lapalma_ca_km=685.0,           # closest approach of La Palma to the shadow centre (App. B.3)
    lapalma_flux=(0.061, 0.0072),  # primary and secondary stellar flux there (App. B.3)
    n_points=68446, n_params=54,   # data points and fitted parameters of the global fit
)

# CDS file, station, latitude, longitude, altitude [m], aperture [m], instrument/filter,
# exposure [s], cycle [s], λ [µm], Δλ [µm], used in the global fit (Table A.1)
const TRITON_STATIONS = [
    ("agerolat28cm",           "Agerola 0.28 m",          "40 37 26.2",  "14 33 50.6",   660.0, 0.28, "CCD/clear",          0.86,  0.86,  0.65, 0.4,  true),
    ("albox",                  "Albox",                   "37 24 20.0",  "-02 09 07.0",  485.0, 0.41, "CCD/clear",          2.00,  3.34,  0.65, 0.4,  true),
    ("alcantarilha",           "Alcantarilha",            "37 07 58.8",  "-08 21 54.0",   65.0, 0.36, "CCD/red",            0.33,  0.64,  0.65, 0.15, true),
    ("algiers",                "Algiers",                 "36 47 52.1",  "03 01 56.2",   331.0, 0.81, "video/clear",        0.04,  0.04,  0.65, 0.4,  true),
    ("ariana",                 "Ariana",                  "36 53 03.0",  "10 11 42.0",     5.0, 0.20, "CCD/clear",          0.50,  0.60,  0.65, 0.4,  true),
    ("athens",                 "Athens",                  "37 58 06.8",  "23 47 00.1",   250.0, 0.40, "CCD/clear",          1.00,  3.50,  0.65, 0.4,  true),
    ("belesta",                "Belesta",                 "43 26 39.4",  "01 48 58.6",   235.0, 0.30, "video/clear",        0.64,  0.64,  0.65, 0.4,  true),
    ("boissettes",             "Boissettes",              "48 31 41.0",  "02 37 28.0",    75.0, 0.30, "CCD/clear",          0.50,  0.50,  0.65, 0.4,  true),
    ("calaraltovis",           "Calar Alto VIS",          "37 13 24.7",  "-02 32 44.9", 2160.0, 1.23, "sCMOS/0.4-1.0 µm",   0.063, 0.063, 0.7,  0.6,  true),
    ("calaraltoswir",          "Calar Alto SWIR",         "37 13 24.7",  "-02 32 44.9", 2160.0, 1.23, "SWIR/1.0-1.7 µm",    0.11,  0.11,  1.35, 0.7,  true),
    ("calern",                 "Calern",                  "43 45 13.2",  "06 55 22.4",  1268.0, 1.04, "CCD/I'",             0.10,  0.11,  0.8,  0.15, true),
    ("cataniat80cm",           "Catania 0.80 m",          "37 41 35.8",  "14 58 29.1",  1725.0, 0.80, "CCD/clear",          0.80,  1.70,  0.65, 0.4,  true),
    ("clanfield",              "Clanfield",               "50 56 19.2",  "-01 01 10.6",  155.0, 0.61, "video/clear",        0.16,  0.16,  0.65, 0.4,  true),
    ("constancia",             "Constância",              "39 29 41.6",  "-08 19 25.2",  147.0, 0.51, "video/clear",        0.64,  0.64,  0.65, 0.4,  true),
    ("cotesdemeuse",           "Côtes de Meuse",          "49 00 07.2",  "05 41 06.6",   284.0, 0.83, "EMCCD/clear",        0.06,  0.06,  0.6,  0.4,  true),
    ("cuqlesvielmur",          "Cuq les Vielmur",         "44 29 05.0",  "01 55 05.0",   300.0, 0.20, "CCD/clear",          2.00,  4.00,  0.65, 0.4,  true),
    ("darkskyobs",             "Dark Sky Obs.",           "36 15 09.6",  "-81 25 01.2",  960.0, 0.36, "CCD/clear",          1.50,  2.50,  0.65, 0.4,  true),
    ("elvas",                  "Elvas",                   "38 50 47.6",  "-07 12 27.4",  285.0, 0.27, "CCD/clear",          0.30,  0.30,  0.65, 0.4,  true),
    ("felsinaobs",             "Felsina Obs.",            "44 21 22.0",  "11 09 09.0",   652.0, 0.40, "video/clear",        0.64,  0.64,  0.65, 0.4,  true),
    ("forcarei",               "Forcarei",                "42 36 38.4",  "-08 22 15.1",  673.0, 0.51, "video/clear",        0.64,  0.64,  0.65, 0.4,  true),
    ("hartleywintney",         "Hartley Wintney",         "51 18 23.2",  "-00 54 44.9",   65.0, 0.36, "CCD/Johnson V",      1.00,  3.31,  0.55, 0.1,  true),
    ("helmos",                 "Helmos",                  "37 59 08.1",  "22 11 54.6",  2323.0, 2.28, "CCD/V+R",            0.60,  0.674, 0.6,  0.25, true),
    ("hornchurch",             "Hornchurch",              "51 33 06.4",  "00 11 38.8",    14.0, 0.24, "video/clear",        0.64,  0.64,  0.65, 0.4,  true),
    ("ithaca",                 "Ithaca",                  "42 27 29.4",  "-76 23 05.5",  530.0, 0.60, "CCD/clear",          0.20,  0.36,  0.65, 0.4,  true),
    ("javalambreobstx40",      "Javalambre Tx40",         "40 02 30.6",  "-01 00 57.6", 1955.0, 0.40, "CCD/clear",          2.00,  2.50,  0.65, 0.4,  true),
    ("javalambreobsexcalibur", "Javalambre Excalibur",    "40 02 30.6",  "-01 00 57.6", 1955.0, 0.28, "CCD/clear",          2.00,  3.00,  0.65, 0.4,  true),
    ("kryonerir",              "Kryoneri R",              "37 58 19.0",  "22 37 07.0",   930.0, 1.20, "sCMOS/R",            0.19,  0.20,  0.65, 0.15, true),
    ("kryonerii",              "Kryoneri I",              "37 58 19.0",  "22 37 07.0",   930.0, 1.20, "sCMOS/I",            0.19,  0.20,  0.8,  0.15, true),
    ("lapalma",                "La Palma",                "28 45 43.2",  "-17 52 39.3", 2340.0, 2.00, "CCD/I+z",            0.60,  0.635, 0.85, 0.25, true),
    ("lacarolina",             "La Carolina",             "38 16 27.0",  "-03 36 55.0",  595.0, 0.64, "CCD/clear",          0.15,  0.15,  0.65, 0.4,  true),
    ("larocheyon",             "La Roche-sur-Yon",        "46 32 48.4",  "-01 19 37.8",   48.0, 0.40, "CMOS/clear",         0.20,  0.20,  0.65, 0.4,  true),
    ("lebeausset",             "Le Beausset",             "43 11 38.1",  "05 48 05.9",   197.0, 0.25, "video/clear",        0.64,  0.64,  0.65, 0.4,  true),
    ("marrakech",              "Marrakech",               "31 35 16.2",  "-08 00 46.9",  494.0, 0.60, "CCD/clear",          2.50,  5.50,  0.65, 0.4,  true),
    ("nancy",                  "Nancy",                   "48 39 51.6",  "06 09 28.7",   284.0, 0.25, "CCD/luminance",      1.50,  2.10,  0.55, 0.3,  true),
    ("newark",                 "Newark",                  "43 00 23.7",  "-77 07 06.5",  165.0, 0.25, "video/clear",        0.53,  0.53,  0.65, 0.4,  true),
    ("oukaimeden",             "Oukaïmeden",              "31 12 23.0",  "-07 51 59.0", 2727.0, 0.50, "CCD/clear",          2.00,  6.80,  0.65, 0.4,  true),
    ("paris",                  "Paris",                   "48 52 17.0",  "02 23 07.0",    92.0, 0.15, "CCD/clear",          0.25,  0.25,  0.65, 0.4,  true),
    ("picdumidit60cm",         "Pic du Midi 0.60 m",      "42 56 12.0",  "00 08 31.9",  2862.0, 0.60, "EMCCD/visible",      0.05,  0.05,  0.6,  0.4,  true),
    ("puimichel",              "Puimichel",               "43 58 53.1",  "06 02 10.0",   714.0, 0.60, "EMCCD/visible",      0.10,  0.10,  0.6,  0.4,  true),
    ("reading",                "Reading",                 "51 30 23.8",  "-00 48 58.1",   80.0, 0.30, "video/red",          0.64,  0.64,  0.65, 0.15, true),
    ("sabadell",               "Sabadell",                "41 33 00.2",  "02 05 24.8",   224.0, 0.50, "video/clear",        0.16,  0.16,  0.65, 0.4,  true),
    ("stcapraist20cm",         "Saint Caprais 0.20 m",    "43 52 25.9",  "01 43 05.8",   185.0, 0.20, "CCD/Johnson V",      1.00,  1.58,  0.55, 0.1,  true),
    ("stcapraist94cm",         "Saint Caprais 0.94 m",    "43 52 25.9",  "01 43 05.8",   185.0, 0.94, "SWIR/clear",         0.05,  0.08,  1.35, 0.7,  true),
    ("seyssest31cm",           "Seysses 0.31 m",          "43 30 06.3",  "01 17 20.3",   183.0, 0.31, "video/UV-IR block",  0.32,  0.32,  0.6,  0.3,  true),
    ("snevadat150cm",          "Sierra Nevada 1.50 m",    "37 03 51.0",  "-03 23 05.0", 2925.8, 1.50, "video/clear",        0.04,  0.04,  0.65, 0.4,  true),
    ("southampton",            "Southampton",             "50 55 18.9",  "-01 22 28.1",   16.0, 0.40, "CMOS/clear",         0.25,  0.25,  0.65, 0.4,  true),
    ("stuttgart",              "Sternwarte Stuttgart",    "48 46 56.7",  "09 11 47.4",   346.0, 0.41, "video/clear",        0.64,  0.64,  0.65, 0.4,  true),
    ("steyning",               "Steyning",                "50 53 16.0",  "-00 19 57.0",   20.0, 0.28, "CCD/luminance",      2.00,  4.00,  0.55, 0.3,  true),
    ("tournefeuille",          "Tournefeuille",           "43 34 58.8",  "01 19 35.4",   163.0, 0.32, "CCD/clear",          0.25,  0.25,  0.65, 0.4,  true),
    ("uclobservatory",         "UCL Observatory",         "51 36 47.6",  "-00 14 32.3",   82.0, 0.35, "CCD/luminance",      1.00,  2.38,  0.55, 0.3,  true),
    ("valencia",               "Valencia",                "39 56 42.0",  "-01 06 05.4", 1280.0, 0.50, "CCD/clear",          1.00,  2.50,  0.65, 0.4,  true),
    ("varages",                "Varages",                 "43 36 44.6",  "05 57 49.1",   300.0, 0.50, "InGaAs/RG1000",      2.30,  2.30,  1.3,  0.6,  true),
    # light curves not used in the global fit (low S/N, contamination by Neptune)
    ("abingdon",               "Abingdon",                "51 37 53.1",  "-01 16 55.2",   59.0, 0.30, "video/clear",        0.16,  0.16,  0.65, 0.4,  false),
    ("agerolat50cm",           "Agerola 0.50 m",          "40 37 26.2",  "14 33 50.6",   660.0, 0.50, "CMOS/clear",         0.50,  0.50,  0.65, 0.4,  false),
    ("baronnies",              "Obs. des Baronnies Provençales", "44 24 29.7", "05 30 54.4", 820.0, 0.80, "CCD/clear",     2.00,  2.00,  0.65, 0.4,  false),
    ("calaraltot36cm",         "Calar Alto 0.36 m",       "37 13 24.7",  "-02 32 44.9", 2160.0, 0.36, "CCD/clear",          1.00,  1.00,  0.65, 0.4,  false),
    ("caserta",                "Caserta",                 "41 16 11.3",  "13 56 28.9",   407.0, 0.24, "CCD/clear",          2.00,  2.46,  0.65, 0.4,  false),
    ("cataniat28cm",           "Catania 0.28 m",          "37 30 24.0",  "15 05 35.5",    50.0, 0.28, "video/clear",        0.64,  0.64,  0.65, 0.4,  false),
    ("fehrenbachobs",          "Charles Fehrenbach Obs.", "50 05 01.5",  "03 46 34.4",   191.0, 0.36, "CCD/clear",          0.28,  0.30,  0.65, 0.4,  false),
    ("glyfadaathens",          "Glyfada-Athens",          "37 52 33.6",  "23 46 12.0",   100.0, 0.35, "CCD/clear",          1.00,  1.15,  0.65, 0.4,  false),
    ("guirguillano",           "Guirguillano",            "42 42 42.2",  "-01 51 54.3",  595.0, 0.31, "video/clear",        2.56,  2.56,  0.65, 0.4,  false),
    ("hamseygreen",            "Hamsey Green",            "51 19 09.4",  "-00 04 01.4",  170.0, 0.28, "video/clear",        1.00,  1.00,  0.65, 0.4,  false),
    ("houten",                 "Houten",                  "52 01 59.3",  "05 09 44.2",     2.0, 0.36, "CCD/685 nm long pass", 0.41, 0.41, 0.8,  0.2,  false),
    ("lahitat77cm",            "La Hita 0.77 m",          "39 34 06.8",  "-03 11 09.5",  695.0, 0.77, "CCD/luminance",      1.00,  3.55,  0.55, 0.3,  false),
    ("lahitat40cm",            "La Hita 0.40 m",          "39 34 06.8",  "-03 11 09.5",  695.0, 0.40, "video/clear",        0.30,  0.30,  0.65, 0.4,  false),
    ("lamandine",              "Lamandine",               "44 12 35.0",  "01 42 01.7",   297.0, 0.30, "CCD/clear",          0.06,  0.06,  0.65, 0.4,  false),
    ("lasagra",                "La Sagra",                "37 58 58.1",  "-02 34 05.1", 1530.0, 0.36, "CCD or video/clear", NaN,   NaN,   0.65, 0.4,  false),
    ("leeds",                  "Leeds",                   "53 50 15.6",  "-01 36 27.8",  113.0, 0.20, "video/IR-UV block",  1.28,  1.28,  0.6,  0.3,  false),
    ("liasurania",             "Lias",                    "43 33 33.3",  "01 06 34.4",   300.0, 0.40, "CCD/clear",          0.13,  0.13,  0.65, 0.4,  false),
    ("louargat",               "Louargat",                "48 32 16.7",  "-03 21 30.5",  196.0, 0.30, "CMOS/L",             0.20,  0.20,  0.55, 0.3,  false),
    ("massa",                  "Massa",                   "44 01 17.2",  "10 07 56.7",    30.0, 0.20, "video/clear",        5.12,  5.12,  0.65, 0.4,  false),
    ("montignybret",           "Montigny le Bretonneux",  "48 45 54.0",  "02 00 52.0",   168.0, 0.28, "CCD/clear",          0.50,  0.51,  0.65, 0.4,  false),
    ("mtagliale",              "Mount Agliale Obs.",      "43 59 43.1",  "10 30 53.8",   750.0, 0.50, "CCD/clear",          1.50,  3.23,  0.65, 0.4,  false),
    ("nerpio",                 "Nerpio",                  "38 09 56.0",  "-02 19 35.0", 1650.0, 0.32, "CCD/clear",          1.00, 16.11,  0.65, 0.4,  false),
    ("neutraubling",           "Neutraubling",            "48 59 23.1",  "12 12 57.3",   333.0, 0.28, "video/clear",        0.16,  0.16,  0.65, 0.4,  false),
    ("nice",                   "Nice",                    "43 43 32.9",  "07 17 59.4",   350.0, 0.13, "video/clear",        2.00,  2.00,  0.65, 0.4,  false),
    ("visker",                 "Obs. de Biscarmiau (Visker)", "43 08 40.4", "00 03 31.8", 488.0, 0.31, "video/clear",      0.64,  0.64,  0.65, 0.4,  false),
    ("monplaisirarles",        "Obs. Monplaisir",         "43 40 58.3",  "04 38 32.6",     5.0, 0.28, "CCD/clear",          1.03,  1.49,  0.65, 0.4,  false),
    ("orcemont",               "Orcemont",                "48 35 28.0",  "01 48 45.0",   165.0, 0.10, "CCD/IR cut",        10.00, 10.15,  0.55, 0.3,  false),
    ("overath",                "Overath",                 "50 57 11.4",  "07 14 53.1",   200.0, 0.20, "video/clear",        0.32,  0.32,  0.65, 0.4,  false),
    ("pfunz",                  "Pfünz",                   "48 53 26.0",  "11 16 23.0",   488.0, 0.28, "CCD/clear",          1.00,  1.20,  0.65, 0.4,  false),
    ("picdumidit1m",           "Pic du Midi 1.05 m",      "42 56 12.0",  "00 08 31.9",  2862.0, 1.05, "sCMOS/clear",        0.10,  0.10,  0.65, 0.4,  false),
    ("rosenfeldbritt",         "Rosenfeld-Brittheim",     "48 17 17.8",  "08 40 38.9",   700.0, 0.80, "CMOS/UV-IR block",   0.40,  0.40,  0.6,  0.3,  false),
    ("saintcheron",            "Saint Chéron",            "48 32 16.6",  "02 07 51.7",   160.0, 0.18, "video/clear",        0.08,  0.08,  0.65, 0.4,  false),
    ("stmauricenav",           "Saint-Maurice-Navacelles", "43 50 39.3", "03 33 47.3",   583.0, 0.21, "CCD/clear",          0.64,  0.64,  0.65, 0.4,  false),
    ("stmichel",               "Saint Michel",            "43 54 23.3",  "05 43 34.7",   564.0, 0.60, "EMCCD/clear",        0.05,  0.05,  0.6,  0.4,  false),
    ("selztal",                "Selztal Obs.",            "49 50 00.0",  "08 15 00.0",   110.0, 0.51, "CMOS/clear",         0.80,  0.80,  0.65, 0.4,  false),
    ("seyssest23cm",           "Seysses 0.23 m",          "43 30 06.3",  "01 17 20.3",   183.0, 0.23, "video/UV-IR block",  0.32,  0.32,  0.6,  0.3,  false),
    ("snevadat90cm",           "Sierra Nevada 0.90 m",    "37 03 51.0",  "-03 23 05.0", 2925.8, 0.90, "CCD/clear",          1.00,  3.33,  0.65, 0.4,  false),
    ("valledaosta",            "Valle d'Aosta",           "45 47 23.1",  "07 28 42.3",  1675.0, 0.81, "EMCCD/clear",        0.20,  0.20,  0.6,  0.4,  false),
]

"""
    triton_star()

UCAC4 410143659 = Gaia DR2 2610107911326516992: geocentric ICRF position at the epoch of the
event (Table 1: proper motion and parallax already applied), V = 12.7, G = 12.5; apparent
diameter 0.65 km at Triton's distance (0.031 mas).
"""
triton_star() = Star("22 54 18.4364", "-08 00 08.318"; name="Gaia DR2 2610107911326516992",
                     source_id="2610107911326516992", epoch=2017.7638, ra_err=0.2, dec_err=0.2,
                     diameter=0.65 / TRITON_PAPER.distance * 206264806.247, mag=Dict("V" => 12.655, "K" => 11.080))

"Ephemeris of Triton from the stored Horizons table."
function triton_ephemeris()
    rows = [split(l) for l in eachline(joinpath(TRITON_DIR, "triton_horizons_20171005.txt")) if !startswith(l, "#")]
    return EphemTable([r[1] for r in rows], parse.(Float64, getindex.(rows, 2)), parse.(Float64, getindex.(rows, 3)),
                      parse.(Float64, getindex.(rows, 4)); name="Triton")
end

"Triton as a `Body`: radius 1353 km, GM = 1428 km³/s² (Table 2)."
triton_body() = Body("Triton", triton_ephemeris(); radius=TRITON_PAPER.radius,
                     density=TRITON_PAPER.GM * 1e9 / (6.6743e-11 * 4π / 3 * (TRITON_PAPER.radius * 1e3)^3))

# The light curves (triton_cds/<file>.dat) are read once they have been downloaded and their
# format checked against the CDS ReadMe; until then `triton_event` builds the event without them.

"""
    triton_event(; used_only=true) -> Event

The occultation of 2017 October 5 with the stations of Table A.1 (without light curves until the
CDS files are available).
"""
function triton_event(; used_only::Bool=true)
    ev = Event(triton_star(), triton_body(), TRITON_CA; window=1200, name="Triton, 2017-10-05")
    for s in TRITON_STATIONS
        used_only && !s[12] && continue
        file, name, lat, lon, alt = s[1], s[2], s[3], s[4], s[5]
        push!(ev, Chord(name, Observer(name, lon, lat, alt); span=(Epoch(TRITON_CA) - 300, Epoch(TRITON_CA) + 300)))
    end
    return ev
end
