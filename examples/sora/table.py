# Table: published | Occultations.jl | SORA (radial residuals, same data and σ_limb)
import os
DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'output', 'sora', 'cases')
# published: a′, σ; ε′, σ; PA(minor, N→E, mod 180), σ   (1σ; MS4 3σ/3; UK126 full axes/2)
PAPER = {
 'varda_2018': ('Souami et al. 2020', 383.0, 3.0, 0.066, 0.047, 67.0, 8.0),
 '1995sm55_2024': ('Ortiz et al. 2025', 104.3, 0.4, 0.200, 0.002, 44.1, 0.4),
 '1995sm55_2024_slimb': ('Ortiz et al. 2025', 104.3, 0.4, 0.200, 0.002, 44.1, 0.4),
 'huya_2023': ('Rommel et al. 2025, free', 222.5, 9.1, 0.107, 0.075, 31.5, 15.5),
 '2002tc302_2018': ('Ortiz et al. 2020', 271.6, 9.0, 0.154, 0.035, 3.0, 1.0),
 '2007uk126_2014': ('Schindler et al. 2017', 322.90, 2.84, 0.074, 0.021, 158.75, 5.65),
 '2002gz32_2017': ('Santos-Sanz et al. 2021, orig. chords', 146.0, 12.0, 0.500, 0.045, 138.0, 3.0),
 'makemake_2011': ('Ortiz et al. 2012, circle', 715.0, 4.5, 0.0, 0.0, 0.0, 0.0),
 '2002ms4_2020': ('Rommel et al. 2023 (3σ/3)', 412.0, 3.3, 0.066, 0.011, 121.0, 5.3),
}
def kv(path):
    d = {}
    for l in open(path):
        k, *v = l.split(); d[k] = [float(x) if x.replace('.','',1).replace('-','',1).replace('e','',1).isdigit() or x.lstrip('-')[0].isdigit() else x for x in v]
    return d
NAMES = {'varda_2018': 'Varda 2018', '1995sm55_2024': '1995 SM55 2024', '1995sm55_2024_slimb': '1995 SM55 2024, σ_limb',
         'huya_2023': 'Huya 2023', '2002tc302_2018': '2002 TC302 2018', '2007uk126_2014': '2007 UK126 2014',
         '2002gz32_2017': '2002 GZ32 2017', 'makemake_2011': 'Makemake 2011 (circle)', '2002ms4_2020': '2002 MS4 2020'}
pm = lambda v, s, d: f"{v:.{d}f} ± {s:.{d}f}"
print("Parameters (published | Occultations.jl | SORA):\n")
print("| Event | σ_limb [km] | N | a′ [km] | ε′ | PA [°] |")
print("|---|---|---|---|---|---|")
diag = []
for n, (ref, a, sa, e, se, pa, spa) in PAPER.items():
    c = kv(os.path.join(DIR, n, 'config.txt')); s = kv(os.path.join(DIR, n, 'sora.txt'))
    circ = c['model'][0] == 'circle'
    col = lambda q, d: (pm(q[0], q[1], d))
    A = " / ".join([pm(a, sa, 1), col(c['oj_equatorial_radius'], 1), col(s['mc_equatorial_radius'], 1)])
    E = "–" if circ else " / ".join([pm(e, se, 3), col(c['oj_oblateness'], 3), col(s['mc_oblateness'], 3)])
    P = "–" if circ else " / ".join([pm(pa, spa, 1), pm(c['oj_position_angle'][0] % 180, c['oj_position_angle'][1], 1),
                                       pm(s['mc_position_angle'][0] % 180, s['mc_position_angle'][1], 1)])
    print(f"| {NAMES[n]} | {c['ellipse_error'][0]:.1f} | {int(s['npoints'][0])} | {A} | {E} | {P} |")
    diag.append((NAMES[n], c['oj_chi2'][0], s['sora_chi2_at_oj'][0], s['ls_rough'][0], s['scipy_r2'][0], s['scipy_r4'][0],
                 s['mc_chi2min'][0], int(s['mc_n1sigma'][0]), int(s['mc_nsamples'][0])))
print("\nχ² (radial residuals):\n")
print("| Event | Occultations.jl minimum | SORA's formula there | SORA least squares | scipy, Σ r² | scipy, Σ r⁴ | SORA Monte Carlo minimum | samples with Δχ² < 1 |")
print("|---|---|---|---|---|---|---|---|")
for d in diag:
    print(f"| {d[0]} | {d[1]:.2f} | {d[2]:.2f} | {d[3]:.2f} | {d[4]:.2f} | {d[5]:.2f} | {d[6]:.2f} | {d[7]} of {d[8]} |")
