# Speed of SORA on the Varda occultation of 2018 September 10 (Souami et al. 2020), with the
# same inputs as bench_julia.jl (run that first: it writes varda_ephem_sora.txt and
# lightcurve.txt).
#
#   <venv>/bin/python -I examples/sora/bench_sora.py

import os, sys, time, platform, warnings, multiprocessing
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import Angle

warnings.filterwarnings('ignore')
# SORA's fit_ellipse opens a process pool; Python 3.14 defaults to forkserver, which re-imports this
# script. 'fork' was the Linux default up to Python 3.13.
multiprocessing.set_start_method('fork')
import sora
from sora import Star, Body, Observer, LightCurve, Occultation
from sora.ephem import EphemPlanete
from sora.occultation import fit_ellipse
warnings.simplefilter('ignore')      # SORA resets the filters on import

DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'output', 'sora', 'benchmark')


def bench(f, n):
    """First call and the median of n further calls [s]."""
    t0 = time.perf_counter()
    r = f()
    t1 = time.perf_counter() - t0
    ts = []
    for _ in range(n):
        t0 = time.perf_counter()
        f()
        ts.append(time.perf_counter() - t0)
    return r, t1, float(np.median(ts)), n


def report(name, b):
    print(f"{name:<44s} first {b[1]:9.4f} s   median {b[2]:9.5f} s  (n = {b[3]})")


# station, longitude (E), latitude, altitude [m], ingress, σ, egress, σ (Souami et al. 2020, Tables 3, 4)
CHORDS = [
    ("TWF", "-114 28 13.1", "42 35 01.9", 1133.0, "03:32:57.10", 0.60, "03:33:58.83", 0.60),
    ("MHV", "-114 35 48.9", "35 01 54.1", 184.0, "03:33:56.99", 0.76, "03:35:05.27", 0.30),
    ("YMA", "-114 26 10.5", "32 39 34.0", 97.0, "03:34:17.10", 0.60, "03:35:26.39", 0.59),
    ("CAR", "-111 57 07.9", "33 48 42.9", 654.0, "03:34:12.39", 0.47, "03:35:19.52", 0.47),
    ("FLO", "-111 21 00.6", "33 00 54.3", 484.0, "03:34:23.30", 0.13, "03:35:26.75", 0.15),
]
DAY = "2018-09-10 "


def dms(s):
    d, m, x = s.split()
    sign = -1 if d.startswith("-") else 1
    return sign * (abs(float(d)) + float(m) / 60 + float(x) / 3600)


def make_occ():
    star = Star(coord="17 18 25.12481 -02 05 14.4282", local=True, nomad=False, epoch=Time("2018-09-10 03:34:00"), verbose=False)
    eph = EphemPlanete(os.path.join(DIR, "varda_ephem_sora.txt"), name="Varda")
    body = Body(name="Varda", spkid="2174567", database=None, ephem=eph)
    occ = Occultation(star=star, body=body, time="2018-09-10 03:34:00")
    # EphemPlanete (a local table) gives geocentric positions only; SORA's default topocentric
    # projection would put every station at the geocentre. The older "geocenter" projection
    # (station offsets in ksi, eta) is the one meant for local tables.
    occ.chords._method = "geocenter"
    for name, lon, lat, alt, ti, si, te, se in CHORDS:
        obs = Observer(name=name, lon=dms(lon) * u.deg, lat=dms(lat) * u.deg, height=alt * u.m)
        lc = LightCurve(name=name, immersion=DAY + ti, immersion_err=si, emersion=DAY + te, emersion_err=se)
        occ.chords.add_chord(name=name, observer=obs, lightcurve=lc)
    return occ


print(f"SORA {sora.__version__}, Python {platform.python_version()}, NumPy {np.__version__}")

# ---- 1. the event
b1 = bench(make_occ, 20)
report("1. event set-up (5 positive chords)", b1)
occ = b1[0]


# ---- 2. (f, g) of the ten contact times
def contacts_fg():
    out = []
    for name, chord in occ.chords.items():
        for k in ("immersion", "emersion"):
            out.append(chord.get_fg(time=k))
    return out


b2 = bench(contacts_fg, 20)
report("2. (f, g) of the 10 contacts", b2)
fg = np.array(b2[0], dtype=float)
cf, cg = fg[:, 0].mean(), fg[:, 1].mean()


# ---- 3. ellipse fit. The χ² is the same as Occultations.jl's radial χ² (checked: 7.458 at its
# solution). Search box a few σ around the solution, as a SORA user would set after a first run.
BOX = dict(center_f=186.5, center_g=-302.7, dcenter_f=20, dcenter_g=20, equatorial_radius=381, dequatorial_radius=15,
           oblateness=0.05, doblateness=0.15, position_angle=65, dposition_angle=40, verbose=False,
           method="least_squares")


def summary(chi):
    r = chi.get_nsigma(sigma=1)
    return (f"   χ²_min {chi.data['chi2'].min():.3f}, {len(chi.data['chi2'])} samples, {r['n_points']} within Δχ² < 1: "
            f"a = {r['equatorial_radius'][0]:.1f} ± {r['equatorial_radius'][1]:.1f}, "
            f"ε = {r['oblateness'][0]:.3f} ± {r['oblateness'][1]:.3f}, "
            f"PA = {r['position_angle'][0]:.1f} ± {r['position_angle'][1]:.1f}")


for label, kw, n in (("3a. least squares only (trf)", dict(loop=10, number_chi=1), 20),
                     ("3b. LS + Monte Carlo, dchi_min 10, 10^4", dict(loop=10000, number_chi=10000, dchi_min=10), 5),
                     ("3c. LS + Monte Carlo, 1000 in Δχ² < 1", dict(loop=10000, number_chi=1000, dchi_min=1), 3)):
    b3 = bench(lambda: fit_ellipse(occ, **BOX, **kw), n)
    report(label, b3)
    print(summary(b3[0]))

# ---- 6. light-curve timing on the same synthetic light curve
t, flux = np.loadtxt(os.path.join(DIR, "lightcurve.txt"), unpack=True)
D = 46.6335


def lcfit(method):
    lc = LightCurve(name="FLO-like", time=t, flux=flux, tref="2018-09-10 03:34:55", exptime=0.5, vel=10.86,
                    dist=D, d_star=0.82)
    return lc.occ_lcfit(tmin=-60, tmax=60, immersion_time=-31.0, emersion_time=31.0, delta_t=2.0, loop=10000,
                        method=method)


for method, n in (("chisqr", 3), ("least_squares", 5)):
    b6 = bench(lambda: lcfit(method), n)
    report(f"6. light-curve timing ({method})", b6)
    r = b6[0].get_nsigma(sigma=1)
    print("   immersion", r['immersion'], " emersion", r['emersion'], "(truth −31.700, +31.700)")
