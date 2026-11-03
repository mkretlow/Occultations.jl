# Published reductions re-done with SORA (fit_ellipse), from the events exported by export.jl.
#
#   <venv>/bin/python -I examples/sora/compare_sora.py [case ...]      (after export.jl)
#
# Per case: (1) SORA's least squares from a rough start (centre at the mean of the chord ends,
# a from the longest chord, ε = 0.1, PA = 0) and (2) from the Occultations.jl solution;
# (3) SORA's Monte Carlo region (Δχ² < 1, midpoint ± half range) in a box of ±4σ (Occultations.jl
# covariance) around the better least-squares solution; (4) SORA's χ² at the Occultations.jl
# solution (it must equal Occultations.jl's radial χ²).

import os, sys, time, warnings, multiprocessing
import numpy as np
from scipy.optimize import least_squares
import astropy.units as u
from astropy.time import Time

warnings.filterwarnings('ignore')
multiprocessing.set_start_method('fork')     # fit_ellipse's process pool (Python 3.14)
from sora import Star, Body, Observer, LightCurve, Occultation
from sora.ephem import EphemPlanete
from sora.occultation import fit_ellipse
warnings.simplefilter('ignore')

DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'output', 'sora', 'cases')
PARS = ['center_f', 'center_g', 'equatorial_radius', 'oblateness', 'position_angle']


def read_case(name):
    d = os.path.join(DIR, name)
    ev = {'chords': []}
    for l in open(os.path.join(d, 'event.txt')):
        if l.startswith('#'):
            continue
        if l.startswith('star '):
            ev['star'] = [float(x) for x in l.split()[1:]]
        elif l.startswith('tref '):
            ev['tref'] = l.split()[1]
        elif l.startswith('chord|'):
            parts = l.rstrip('\n').split('|')
            lon, lat, h = (float(x) for x in parts[2].split())
            cs = {}
            for p in parts[3:]:
                k, t, s = p.split()
                cs[k] = (t, float(s))
            ev['chords'].append((parts[1], lon, lat, h, cs))
    cfg = {}
    for l in open(os.path.join(d, 'config.txt')):
        k, *v = l.split()
        cfg[k] = v[0] if k == 'model' else [float(x) for x in v]
    return d, ev, cfg


def make_occ(d, ev):
    ra, dec, pmra, pmdec, plx, ep = ev['star']
    star = Star(ra=ra * u.deg, dec=dec * u.deg, pmra=pmra, pmdec=pmdec, parallax=plx,
                epoch=Time(ep, format='jyear'), local=True, nomad=False, verbose=False)
    body = Body(name='target', spkid='1', database=None, ephem=EphemPlanete(os.path.join(d, 'ephem.txt'), name='target'))
    occ = Occultation(star=star, body=body, time=ev['tref'].replace('T', ' '))
    occ.chords._method = 'geocenter'       # EphemPlanete is geocentric: SORA's ksi/eta projection
    for name, lon, lat, h, cs in ev['chords']:
        kw = {}
        for k, (t, s) in cs.items():
            kw[k] = t.replace('T', ' ')
            kw[k + '_err'] = s
        if len(cs) == 1:                    # a single contact: give the light curve's span
            t0 = Time(list(cs.values())[0][0])
            kw['initial_time'], kw['end_time'] = t0 - 300 * u.s, t0 + 300 * u.s
        obs = Observer(name=name, lon=lon * u.deg, lat=lat * u.deg, height=h * u.m)
        occ.chords.add_chord(name=name, observer=obs, lightcurve=LightCurve(name=name, **kw))
    return occ


def points(occ):
    out = []
    for _, c in occ.chords.items():
        for k in ('immersion', 'emersion'):
            if c.is_able[k]:
                f, g, vf, vg = c.get_fg(time=k, vel=True)
                out.append([f, g, np.hypot(vf, vg) * getattr(c.lightcurve, k + '_err')])
    return np.array(out, dtype=float)


def chi2(P, f0, g0, a, obla, phi_deg, err=0.0):
    """SORA's χ² (fit_ellipse) for the points P = (f, g, σ)."""
    fi, gi, si = P.T
    b = a - a * obla
    theta = np.arctan2(gi - g0, fi - f0)
    ang = theta + np.radians(phi_deg)
    r = a * b / np.sqrt((a * np.sin(ang)) ** 2 + (b * np.cos(ang)) ** 2)
    return np.sum(((fi - f0 - r * np.cos(theta)) ** 2 + (gi - g0 - r * np.sin(theta)) ** 2) / (si ** 2 + err ** 2))


def minimise(P, start, err, circle, power):
    """scipy least squares (trf) of Σ r² (power=1) or Σ r⁴ (power=2, SORA's objective) from SORA's
    start and within SORA's bounds; returns χ² = Σ r² and the parameters."""
    x0 = [start[k] for k in PARS]
    d = [start['d' + k] for k in PARS]
    lo = [x0[i] - d[i] for i in range(5)]
    hi = [x0[i] + d[i] for i in range(5)]
    lo[3], hi[3] = (0.0, 1.0) if d[3] >= 0.5 else (max(0.0, lo[3]), min(1.0, hi[3]))
    free = [i for i in range(5) if not (circle and i in (3, 4))]

    def res(y):
        x = list(x0)
        for j, i in enumerate(free):
            x[i] = y[j]
        fi, gi, si = P.T
        f0, g0, a, obla, phi = x
        b = a - a * obla
        th = np.arctan2(gi - g0, fi - f0)
        ang = th + np.radians(phi)
        r = a * b / np.sqrt((a * np.sin(ang)) ** 2 + (b * np.cos(ang)) ** 2)
        return (np.hypot(fi - f0 - r * np.cos(th), gi - g0 - r * np.sin(th)) / np.sqrt(si ** 2 + err ** 2)) ** power

    sol = least_squares(res, [x0[i] for i in free], bounds=([lo[i] for i in free], [hi[i] for i in free]), method='trf')
    x = list(x0)
    for j, i in enumerate(free):
        x[i] = sol.x[j]
    return chi2(P, *x, err=err), x


def best(chi):
    i = int(np.argmin(chi.data['chi2']))
    return float(chi.data['chi2'][i]), [float(chi.data[k][i]) for k in PARS]


def run(name):
    d, ev, cfg = read_case(name)
    circle = cfg['model'] == 'circle'
    err = cfg['ellipse_error'][0]
    oj = [cfg['oj_' + k][0] for k in PARS]
    ojs = [cfg['oj_' + k][1] for k in PARS]
    occ = make_occ(d, ev)
    P = points(occ)
    fix = dict(doblateness=0, dposition_angle=0) if circle else {}
    # (1) rough start, wide windows
    lmax = max(np.hypot(*(P[i, :2] - P[i + 1, :2])) for i in range(0, len(P) - 1, 2)) if len(P) > 1 else 100.0
    rough = dict(center_f=P[:, 0].mean(), center_g=P[:, 1].mean(), equatorial_radius=0.55 * lmax,
                 dcenter_f=lmax, dcenter_g=lmax, dequatorial_radius=0.5 * lmax,
                 oblateness=0.0 if circle else 0.1, doblateness=0.5, position_angle=0.0, dposition_angle=90.0)
    rough.update(fix)
    t0 = time.perf_counter()
    c1 = fit_ellipse(occ, **rough, loop=10, number_chi=1, ellipse_error=err, verbose=False, method='least_squares')
    t_ls = time.perf_counter() - t0
    x1, p1 = best(c1)
    xr2, pr2 = minimise(P, rough, err, circle, 1)     # χ² from the same start: the actual minimum
    xr4, pr4 = minimise(P, rough, err, circle, 2)     # Σ r⁴: reproduces SORA's least squares
    # (2) from the Occultations.jl solution, same windows
    start = dict(rough, center_f=oj[0], center_g=oj[1], equatorial_radius=oj[2],
                 oblateness=oj[3] if not circle else 0.0, position_angle=oj[4] % 180 if not circle else 0.0)
    c2 = fit_ellipse(occ, **start, loop=10, number_chi=1, ellipse_error=err, verbose=False, method='least_squares')
    x2, p2 = best(c2)
    xb, pb = (x1, p1) if x1 <= x2 else (x2, p2)
    # (3) Monte Carlo region around the better solution: ±4σ (Occultations.jl), floors
    w = [max(4 * ojs[0], 2.0), max(4 * ojs[1], 2.0), max(4 * ojs[2], 2.0), min(max(4 * ojs[3], 0.02), 0.5),
         min(max(4 * ojs[4], 2.0), 90.0)]
    box = dict(center_f=pb[0], center_g=pb[1], equatorial_radius=pb[2], oblateness=pb[3], position_angle=pb[4],
               dcenter_f=w[0], dcenter_g=w[1], dequatorial_radius=w[2], doblateness=w[3], dposition_angle=w[4])
    box.update(fix)
    t0 = time.perf_counter()
    c3 = fit_ellipse(occ, **box, loop=20000, number_chi=2000, dchi_min=1, ellipse_error=err, verbose=False,
                     method='least_squares')
    t_mc = time.perf_counter() - t0
    r = c3.get_nsigma(sigma=1)
    x3, _ = best(c3)
    # (4) SORA's χ² at the Occultations.jl solution
    x4 = chi2(P, *oj, err=err)
    with open(os.path.join(d, 'sora.txt'), 'w') as io:
        io.write(f"npoints {len(P)}\n")
        io.write(f"ls_rough {x1:.6f} " + ' '.join(f'{v:.6f}' for v in p1) + f"\nls_rough_time {t_ls:.4f}\n")
        io.write(f"ls_oj {x2:.6f} " + ' '.join(f'{v:.6f}' for v in p2) + "\n")
        io.write(f"scipy_r2 {xr2:.6f} " + ' '.join(f'{v:.6f}' for v in pr2) + "\n")
        io.write(f"scipy_r4 {xr4:.6f} " + ' '.join(f'{v:.6f}' for v in pr4) + "\n")
        io.write(f"mc_chi2min {x3:.6f}\nmc_n1sigma {r['n_points']}\nmc_nsamples {len(c3.data['chi2'])}\nmc_time {t_mc:.4f}\n")
        for k in PARS:
            io.write(f"mc_{k} {r[k][0]:.6f} {r[k][1]:.6f}\n")
        io.write(f"sora_chi2_at_oj {x4:.6f}\n")
    print(f"{name:22s} N={len(P):3d}  OJ χ² {cfg['oj_chi2'][0]:9.3f} (SORA formula {x4:9.3f})  "
          f"SORA LS rough {x1:9.3f}  from OJ {x2:9.3f}  scipy Σr² {xr2:9.3f} Σr⁴ {xr4:9.3f}  MC min {x3:9.3f} [{r['n_points']} in 1σ]  "
          f"a {r['equatorial_radius'][0]:.1f}±{r['equatorial_radius'][1]:.1f}  ε {r['oblateness'][0]:.3f}±{r['oblateness'][1]:.3f}  "
          f"PA {r['position_angle'][0]:.1f}±{r['position_angle'][1]:.1f}  ({t_ls:.2f} s, {t_mc:.1f} s)", flush=True)


cases = sys.argv[1:] or sorted(n for n in os.listdir(DIR) if os.path.isdir(os.path.join(DIR, n)))
for n in cases:
    try:
        run(n)
    except Exception as e:
        print(f"{n:22s} FAILED: {type(e).__name__}: {e}", flush=True)
