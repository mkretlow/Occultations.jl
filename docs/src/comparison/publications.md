# Occultations.jl compared with published reductions

Twenty-three published occultation studies were re-reduced with Occultations.jl from the timings,
stations and stars given in the papers (2002 XV93 and 2007 OC10: from the authors' light curves) (DOIs in
[`references/`](https://github.com/mkretlow/Occultations.jl/blob/main/references/README.md)). The examples, their data and their complete printed
output are:

| Object | Event(s) | Publication | Files |
|---|---|---|---|
| (38628) Huya | 2019-03-18 | [SantosSanz2022](@citet), A&A 664, A130 | [script](https://github.com/mkretlow/Occultations.jl/blob/main/examples/huya_2019.jl) · [data](https://github.com/mkretlow/Occultations.jl/blob/main/examples/data/huya_20190318.jl) · [output](https://github.com/mkretlow/Occultations.jl/blob/main/examples/output/huya_2019/output.txt) |
| (119951) 2002 KX14 | 5 events 2020–2023 | [Rizos2025](@citet), A&A 697, A62 + corrigendum A&A 708, C1 (2026) | [script](https://github.com/mkretlow/Occultations.jl/blob/main/examples/2002kx14_2020_2023.jl) · [data](https://github.com/mkretlow/Occultations.jl/blob/main/examples/data/2002kx14_2020_2023.jl) · [output](https://github.com/mkretlow/Occultations.jl/blob/main/examples/output/2002kx14_2020_2023/output.txt) |
| (143707) 2003 UY117 | 2020-10-23 | [Kretlow2024](@citet), A&A 691, A31 | [script](https://github.com/mkretlow/Occultations.jl/blob/main/examples/2003uy117_2020.jl) · [data](https://github.com/mkretlow/Occultations.jl/blob/main/examples/data/2003uy117_20201023.jl) · [output](https://github.com/mkretlow/Occultations.jl/blob/main/examples/output/2003uy117_2020/output.txt) |
| (612533) 2002 XV93 | 2024-01-10 | [Arimatsu2026](@citet), Nature Astronomy | [script](https://github.com/mkretlow/Occultations.jl/blob/main/examples/2002xv93_2024.jl) · [data](https://github.com/mkretlow/Occultations.jl/blob/main/examples/data/2002xv93_20240110.jl) · [data 2](https://github.com/mkretlow/Occultations.jl/blob/main/examples/data/2002xv93_zenodo/README.md) · [output](https://github.com/mkretlow/Occultations.jl/blob/main/examples/output/2002xv93_2024/output.txt) |
| (24835) 1995 SM55 | 2024-02-25 | [Ortiz2025](@citet), A&A 703, A147 | [script](https://github.com/mkretlow/Occultations.jl/blob/main/examples/1995sm55_2024.jl) · [data](https://github.com/mkretlow/Occultations.jl/blob/main/examples/data/1995sm55_20240225.jl) · [output](https://github.com/mkretlow/Occultations.jl/blob/main/examples/output/1995sm55_2024/output.txt) |
| (38628) Huya | 2023-06-24 | [Rommel2025](@citet), Planet. Sci. J. 6, 48 | [script](https://github.com/mkretlow/Occultations.jl/blob/main/examples/huya_2023.jl) · [data](https://github.com/mkretlow/Occultations.jl/blob/main/examples/data/huya_20230624.jl) · [output](https://github.com/mkretlow/Occultations.jl/blob/main/examples/output/huya_2023/output.txt) |
| (28978) Ixion | 6 events 2020–2023 | [Kilic2026](@citet), A&A 707, A70 | [script](https://github.com/mkretlow/Occultations.jl/blob/main/examples/ixion_2020_2023.jl) · [data](https://github.com/mkretlow/Occultations.jl/blob/main/examples/data/ixion_2020_2023.jl) · [output](https://github.com/mkretlow/Occultations.jl/blob/main/examples/output/ixion_2020_2023/output.txt) |
| (229762) 2007 UK126 | 2014-11-15 | [Schindler2017](@citet), A&A 600, A12 | [script](https://github.com/mkretlow/Occultations.jl/blob/main/examples/2007uk126_2014.jl) · [data](https://github.com/mkretlow/Occultations.jl/blob/main/examples/data/2007uk126_20141115.jl) · [output](https://github.com/mkretlow/Occultations.jl/blob/main/examples/output/2007uk126_2014/output.txt) |
| (174567) Varda | 2018-09-10 | [Souami2020](@citet), A&A 643, A125 | [script](https://github.com/mkretlow/Occultations.jl/blob/main/examples/varda_2018.jl) · [data](https://github.com/mkretlow/Occultations.jl/blob/main/examples/data/varda_20180910.jl) · [output](https://github.com/mkretlow/Occultations.jl/blob/main/examples/output/varda_2018/output.txt) |
| (84522) 2002 TC302 | 2018-01-28 | [Ortiz2020](@citet), A&A 639, A134 | [script](https://github.com/mkretlow/Occultations.jl/blob/main/examples/2002tc302_2018.jl) · [data](https://github.com/mkretlow/Occultations.jl/blob/main/examples/data/2002tc302_20180128.jl) · [output](https://github.com/mkretlow/Occultations.jl/blob/main/examples/output/2002tc302_2018/output.txt) |
| (208996) 2003 AZ84 | 2012-02-03 and 2014-11-15 (+ single chords 2011, 2013) | [DiasOliveira2017](@citet), AJ 154, 22 | [script](https://github.com/mkretlow/Occultations.jl/blob/main/examples/2003az84_2012_2014.jl) · [data](https://github.com/mkretlow/Occultations.jl/blob/main/examples/data/2003az84_2012_2014.jl) · [output](https://github.com/mkretlow/Occultations.jl/blob/main/examples/output/2003az84_2012_2014/output.txt) |
| (84922) 2003 VS2 | 2019-10-22 | [VaraLubiano2022](@citet), A&A 663, A121 | [script](https://github.com/mkretlow/Occultations.jl/blob/main/examples/2003vs2_2019.jl) · [data](https://github.com/mkretlow/Occultations.jl/blob/main/examples/data/2003vs2_20191022.jl) · [output](https://github.com/mkretlow/Occultations.jl/blob/main/examples/output/2003vs2_2019/output.txt) |
| (136108) Haumea | 2017-01-21 | [Ortiz2017](@citet), Nature 550, 219 | [script](https://github.com/mkretlow/Occultations.jl/blob/main/examples/haumea_2017.jl) · [data](https://github.com/mkretlow/Occultations.jl/blob/main/examples/data/haumea_20170121.jl) · [output](https://github.com/mkretlow/Occultations.jl/blob/main/examples/output/haumea_2017/output.txt) |
| (95626) 2002 GZ32 | 2017-05-20 | [SantosSanz2021](@citet), MNRAS 501, 6062 | [script](https://github.com/mkretlow/Occultations.jl/blob/main/examples/2002gz32_2017.jl) · [data](https://github.com/mkretlow/Occultations.jl/blob/main/examples/data/2002gz32_20170520.jl) · [output](https://github.com/mkretlow/Occultations.jl/blob/main/examples/output/2002gz32_2017/output.txt) |
| (136472) Makemake | 2011-04-23 | [Ortiz2012](@citet), Nature 491, 566 | [script](https://github.com/mkretlow/Occultations.jl/blob/main/examples/makemake_2011.jl) · [data](https://github.com/mkretlow/Occultations.jl/blob/main/examples/data/makemake_20110423.jl) · [output](https://github.com/mkretlow/Occultations.jl/blob/main/examples/output/makemake_2011/output.txt) |
| (54598) Bienor | 4 events 2019–2023 | [FernandezValenzuela2023](@citet), A&A 669, A112; [Rizos2024](@citet), A&A 689, A82 | [script](https://github.com/mkretlow/Occultations.jl/blob/main/examples/bienor_2019_2023.jl) · [data](https://github.com/mkretlow/Occultations.jl/blob/main/examples/data/bienor_2019_2023.jl) · [output](https://github.com/mkretlow/Occultations.jl/blob/main/examples/output/bienor_2019_2023/output.txt) |
| (10199) Chariklo | 7 events 2017–2020 | [Morgado2021](@citet), A&A 652, A141 (timings at the CDS) | [script](https://github.com/mkretlow/Occultations.jl/blob/main/examples/chariklo_2017_2020.jl) · [data](https://github.com/mkretlow/Occultations.jl/blob/main/examples/data/chariklo_2017_2020.jl) · [output](https://github.com/mkretlow/Occultations.jl/blob/main/examples/output/chariklo_2017_2020/output.txt) |
| (307261) 2002 MS4 | 2020-08-08 | [Rommel2023](@citet), A&A 678, A167 (stations at the CDS) | [script](https://github.com/mkretlow/Occultations.jl/blob/main/examples/2002ms4_2020.jl) · [data](https://github.com/mkretlow/Occultations.jl/blob/main/examples/data/2002ms4_20200808.jl) · [output](https://github.com/mkretlow/Occultations.jl/blob/main/examples/output/2002ms4_2020/output.txt) |
| (50000) Quaoar | 2011-05-04 and 2022-08-09 | [BragaRibas2013](@citet), ApJ 773, 26; [Pereira2023](@citet), A&A 673, L4 | [script](https://github.com/mkretlow/Occultations.jl/blob/main/examples/quaoar_2011_2022.jl) · [data](https://github.com/mkretlow/Occultations.jl/blob/main/examples/data/quaoar_2011_2022.jl) · [output](https://github.com/mkretlow/Occultations.jl/blob/main/examples/output/quaoar_2011_2022/output.txt) |
| (2060) Chiron | 2018-11-28 and 2019-09-08 | [BragaRibas2023](@citet), A&A 676, A72 | [script](https://github.com/mkretlow/Occultations.jl/blob/main/examples/chiron_2018_2019.jl) · [data](https://github.com/mkretlow/Occultations.jl/blob/main/examples/data/chiron_2018_2019.jl) · [output](https://github.com/mkretlow/Occultations.jl/blob/main/examples/output/chiron_2018_2019/output.txt) |
| (470316) 2007 OC10 | 2022-08-22 (double star) | [GomezLimon2025](@citet), A&A 697, A157 (light curves at the CDS) | [script](https://github.com/mkretlow/Occultations.jl/blob/main/examples/2007oc10_2022.jl) · [data](https://github.com/mkretlow/Occultations.jl/blob/main/examples/data/2007oc10_20220822.jl) · [output](https://github.com/mkretlow/Occultations.jl/blob/main/examples/output/2007oc10_2022/output.txt) |

Each case is also a regression test in `test/runtests.jl`.

## Common setup and caveats

* **Input**: the published ingress/egress times and their 1σ errors, the station coordinates
  and the occulted star (position at the event epoch, angular diameter). The original light
  curves are not public, so the timing step of the reduction is not repeated.
* **Ephemerides**: the current JPL Horizons solution (stored in `examples/data`), while the
  papers used older JPL orbits (JPL#28, #18, #20). The fitted centres (offsets from the
  ephemeris) therefore differ from the published ones; shapes, sizes and the absolute
  astrometric positions do not depend on this.
* **Negative chords**: the observed intervals are not tabulated in the papers; ±5–6 min
  around the event were assumed. Negative chords are hard constraints (the package default);
  the exposure (plus dead) time is stored as the chord's sensitivity for the optional soft
  constraint.
* **Fits**: residuals along the chords (the package default) and radial residuals (as SORA
  and most papers). Uncertainties from the covariance matrix, χ² profiles and posterior
  sampling; where χ²/dof ≫ 1 they are scaled or a limb allowance `σ_limb` is fitted.
* **Position angles**: the package gives the position angle of the minor axis, measured from
  north towards +f (east). The UY117 paper uses the same convention (from north, clockwise
  with f to the right); the angle of the KX14 paper agrees with that of the major axis.

## Geometry

The event geometry computed by the package agrees with the papers (the first three studies are shown here; the other sections give their values):

| Event | Distance [au] (paper) | Shadow velocity [km/s] (paper) | Star diameter [km] (paper) |
|---|---|---|---|
| Huya | 28.350 (28.35) | 8.16 (8.07) | 0.70 (0.70) |
| KX14 A | 38.01 (38.01) | 23.41 (23.41) | 1.27 (1.26) |
| KX14 B | 37.88 (37.88) | 24.33 (24.33) | 0.52 (0.58) |
| KX14 C | 37.97 (37.97) | 20.38 (20.38) | 1.27 (1.25) |
| KX14 D | 37.96 (37.96) | 20.83 (20.83) | 0.88 (0.89) |
| KX14 E | 37.85 (37.85) | 23.74 (23.74) | 1.10 (0.65)¹ |
| UY117 | 33.47 (33.47) | 21.80 (21.80) | 0.43 (0.4) |

The Huya velocity is the geocentric one at closest approach; the paper's value is "as seen
from Earth". ¹ The KX14 paper gives 0.04 mas for star E, which is 1.10 km at 37.85 au; its
0.65 km is inconsistent with that.

The Fresnel scales (1.13 km for Huya at 0.6 µm, 1.32 km for UY117 at 0.7 µm) agree as well.

## (38628) Huya, 2019 March 18

19 positive chords (18 sites), 2 negative chords.

| Fit | a [km] | b [km] | PA minor [°] | D_eq [km] | χ²/dof |
|---|---|---|---|---|---|
| **Paper** (Table A.4) | 217.6 ± 3.5 | 194.1 ± 6.1 | 55.2 ± 9.1 | 411.0 ± 7.3 | – |
| Timing errors only | 217.2 ± 0.1 | 192.4 ± 0.3 | 60.1 ± 0.4 | 408.8 ± 0.3 | 227.8 |
| Radial, σ_limb = 7.3 km | 218.1 ± 2.1 | 192.7 ± 5.9 | 62.9 ± 7.3 | 410.0 ± 6.2 | 1.00 |
| Along chords, σ_limb = 8.4 km | 218.2 ± 2.1 | 194.7 ± 4.0 | 59.8 ± 6.8 | 412.2 ± 4.6 | 1.00 |
| Clock priors N(0, 1 s), σ_limb = 10.1 km | 217.8 ± 3.0 | 195.6 ± 5.3 | 61.0 ± 9.8 | 412.8 ± 5.7 | 1.00 |
| Posterior (along chords, σ_limb) | | 196.1 −3.8 +4.0 | 1σ profile 52–66 | 413.7 −4.4 +4.6 | |

Findings:

* The published ellipse is reproduced within its uncertainties by every variant.
* With the published timing errors alone the chords are inconsistent (χ²/dof ≈ 230). Many
  stations show residuals of the **same sign at immersion and emersion** (1–2 s), i.e.
  clock errors, which matches the paper's remark that the chord centres are not aligned.
  Two remedies give the same shape: a limb allowance chosen for χ²/dof = 1
  (`σ_limb = :auto`, 7–8 km) or a Bayesian clock offset per station with a 1 s prior. Fitting
  free clock offsets for all stations is degenerate (the ellipse can shear).
* Both negative chords are consistent with the limb (18.7 and 31.4 km outside).
* Astrometric position at 00:42:35.2 UTC: RA 16 41 06.41994, Dec −06 43 34.6071
  (σ 0.21, 0.14 mas including the star).

## (119951) 2002 KX14, five events 2020–2023

15 positive chords, 2 negative chords. Corrigendum: the Occ. D timings of GOA and Sonear3
(Table B.2) were typographical errors; the corrected values are used (they are also in
arXiv:2503.19179v3).

| Fit | a [km] | b [km] | PA major [°] | D_eq [km] | p_V |
|---|---|---|---|---|---|
| **Paper**, Occ. A alone (Table 2) | 241.5 ± 11.0 | 159.0 ± 8.1 | 105.8 ± 5.8 | – | |
| Occ. A alone (errors × √(χ²/dof), χ²/dof = 13.4) | 248.8 ± 18.0 | 167.5 ± 11.9 | 118.6 ± 7.9 | 408.3 ± 24.8 | |
| **Paper**, all events (Fig. 1) | 241.0 ± 7.2 | 157.1 ± 5.2 | – | 389.2 ± 8.7 | 0.119 ± 0.007 |
| Joint fit, one ellipse, centre per event (scaled, χ²/dof = 7.8) | 245.0 ± 12.0 | 169.0 ± 9.2 | 117.9 ± 6.3 | 406.9 ± 17.9 | |
| Posterior (timing errors only) | | | | 409.2 −6.4 +7.2 | 0.108 ± 0.005 |

Findings:

* The semi-major axis agrees; the semi-minor axis is larger (169 vs 157 km), so the
  diameter is ≈ 2σ (of the published error) larger and the albedo correspondingly lower.
  The methods differ: the paper fits event A, ties the other events to its ellipse with a
  penalty, re-centres all chords and refits an algebraic ellipse with Monte-Carlo errors;
  here all events are fitted at once with a separate centre per event.
* **Event D has a mirror solution**: with two chords the centre can lie on either side. The
  global minimum (χ² 117.8) puts it on the far side; a prior that D's ephemeris offset is
  within ~30 km of event C's (two days earlier) selects the other solution (χ² 124.6), in
  which D's centre lies 73 km south of C's — the paper has 68 km.
* Event E's centre is ~1500 km off, because of the star (Gaia RUWE 1.5), as noted in the
  paper; the per-event centres absorb it.
* The timings scatter more than their errors (χ²/dof ≈ 8); the posterior based on the
  timing errors alone is therefore narrower than the scaled covariance errors.
* Transcription: Apache Point Observatory is listed at 32°46′ **S** in Table B.1; the
  observatory is in the northern hemisphere.

Astrometric offsets from the Horizons ephemeris (Δα cos δ, Δδ in mas): A (+0.5, +9.2),
B (+0.2, +6.3), C (+1.2, +5.4), D (−0.6, +2.8), E (+49.9, +33.8, star problem).

## (143707) 2003 UY117, 2020 October 23

4 positive chords (3 sites), 13 negative chords, rotational light curve with
P = 12.376 h and Δm = 0.36 ± 0.13 mag.

| Fit | a′ [km] | b′ [km] | φ′ [°] | ε′ | D_eq [km] | χ²/dof |
|---|---|---|---|---|---|---|
| **Paper** (Table 3) | 141 ± 9 | 94 ± 16 | −31 ± 11 | 0.33 ± 0.11 | 228 ± 21 | 0.43 |
| Radial residuals | 139.2 ± 9.6 | 92.9 ± 12.9 | −27.3 ± 11.3 | 0.33 ± 0.12 | 227.4 ± 13.5 | 0.41 |
| Along chords | 138.2 ± 5.2 | 94.2 ± 7.2 | −29.8 ± 6.4 | 0.32 ± 0.05 | 228.2 ± 10.4 | 2.10 |
| Posterior (radial) | | | | | 236.6 −15.0 +23.3 | |

Geometric albedo (posterior, H_V = 5.97 ± 0.07): 12.8 −2.2 +2.1 % (paper 13.9 ± 2.7 %).

Findings:

* The limb fit, including χ²/dof, is reproduced. The closest negative chords (CAHA220,
  Albox, Ibiza) pass 8–9 km outside the limb.
* **Astrometry**: with the Gaia DR3 star position the RA differs from the published one by
  ≈ 7 mas; with the Gaia DR2 position (Table 1) the result is RA 03 23 26.21408,
  Dec +22 47 19.3189 against the published 03 23 26.21404, +22 47 19.3184 — the published
  position was evidently computed with the DR2 star.
* **3-D shape**: a triaxial ellipsoid (rotating about its shortest axis, principal = true)
  fitted to the occultation and the light-curve amplitude with the pole free. One
  apparent ellipse and one amplitude do not determine a, b, c and the pole: least-squares
  fits from different starts all reach the same χ² along a **valley of solutions**:

  | a [km] | b [km] | c [km] | aspect [°] | D_eq,V [km] |
  |---|---|---|---|---|
  | 159.4 | 104.5 | 91.6 | 70 | 230.3 |
  | 164.2 | 101.7 | 90.6 | 65 | 229.6 |
  | 166.2 | 98.6 | 91.1 | 62 | 228.6 |
  | 190.8 | 133.7 | 91.1 | 81 | 264.9 |
  | 190.9 | 109.2 | 81.4 | 65 | 238.5 |
  | **Paper**: 166 ± 12 | 108 ± 12 | 90 +14 −12 | 70 +20 −12 | 235 ± 25 |

  The paper's solution lies in this valley (its errors are conditional, one parameter at a
  time). The posterior over the valley has an isotropic prior on the pole. Adding the
  rotational phase at the occultation time (0.63 after maximum light,
  `RotationalPhase`) narrows it:

  | Posterior | a [km] | b [km] | c [km] | D_eq,V [km] | aspect, folded [°] |
  |---|---|---|---|---|---|
  | occultation + amplitude | 168.0 −21.9 +30.5 | 118.8 −16.5 +18.6 | 93.1 −18.1 +12.7 | 242.8 −23.9 +33.3 | 70.7 −19.0 +13.8 |
  | + rotational phase | 168.1 −20.3 +35.4 | 113.9 −12.2 +17.3 | 91.6 −20.1 +10.3 | 240.9 −20.8 +19.9 | 69.0 −15.6 +12.8 |
  | **Paper** | 166 ± 12 | 108 ± 12 | 90 +14 −12 | 235 ± 25 | 70 +20 −12 |

  Both posteriors agree with the paper's ellipsoid. The negative chords are hard
  constraints, and the three close to the southern limb (CAHA220, Albox, Ibiza, 8–9 km away)
  already exclude most of the long, nearly end-on solutions. The phase narrows the
  distribution of D_eq,V further, mainly its upper tail (+33 → +20 km). The aspect is
  determined only up to a pole flip (ψ ↔ 180° − ψ, which leaves the occultation and the
  amplitude unchanged) and is reported folded. The phase does not narrow it.

## (612533) 2002 XV93, 2024 January 10

Two positive chords (Kyoto, Kiso) and a near miss (Fukushima), with the authors' normalised
light curves and best-fit model profiles (Zenodo doi:10.5281/zenodo.18976431, CC BY 4.0).
Geometry: distance 37.04 au (37.0), shadow velocity 25.17 km/s (25.2), Fresnel scale
1.29 km at 0.6 µm.

| Quantity | Paper | Occultations.jl |
|---|---|---|
| Kiso duration | 16.53 ± 0.08 s | 16.528 ± 0.043 s |
| Kyoto duration | 18.0 ± 0.9 s | 17.904 ± 0.451 s |
| Circle from the chord ends | 235 +22 −15 km, χ² = 1.40 | 230.8 ± 8.7 km, χ² = 5.6 (1 dof) |
| Closest approach Kyoto / Kiso / Fukushima | 41 / 103 / 255 km | 32 / 94 / 245 km |
| Model profile vs authors' (ρ = 1500 kg/m³) | – | rms 0.009–0.010, max 0.035 |
| p(CH₄) | 124 +15 −12 nbar | 129.3 ± 12.3 nbar |
| p(N₂) | 177 +18 −19 nbar | 185.3 ± 19.9 nbar |
| p(CO) | 159 +15 −17 nbar | 166.1 ± 20.1 nbar |
| Surface radius | 244 ± 1 km (Kiso alone) | 235.2 ± 2.0 km (CH₄, three light curves); χ² profile flat from ≈ 230 to 250 km |
| Δχ² airless − atmosphere | 89.2 − 28.7 = 60.5 (Kiso) | 57.7 (free airless fit), 64.9 (profile at p = 0): 8.1σ |
| CH₄ pressure range | – | 95 % one-sided limits 104–157 nbar, 3σ 81–176 nbar |

Findings:

* The timings, the geometry and the circle agree within the errors. The circle through
  four chord ends has a single degree of freedom, so its χ² carries little information.
* **Atmosphere model**: the same model of the paper was used:
  * 47 K at the surface and +2.7 K/km up to 20 km;
  * free molecular flow above 25 km, with smooth 2 km transitions.

  Ray-traced with this package, it reproduces the authors' profiles to 1 % rms, with the
  shadow edge at 229.9 vs 230.1 km. For ρ = 1000 kg/m³ the agreement is worse (rms 0.04, edge 1.5 km apart).
  Presumably the authors' ρ = 1000 profiles use a different radius from the 244 km assumed
  here.
* **Pressures** for all three gases agree with the paper within about half their errors
  (CH₄ 129 vs 124, N₂ 185 vs 177, CO 166 vs 159 nbar). Here all three light curves are
  fitted jointly, with a free centre and radius.
* **Composition**: the three gases fit equally well (χ² 171.3, 171.2, 171.1 for 178 dof).
  As in the paper, the light curves do not distinguish them; they constrain the
  refractivity, and with it a different surface pressure for each gas.
* **Radius**: the best joint fit has the surface at 235 km, against 244 km in the paper.
  The χ² profile in the radius is nearly flat between ≈ 230 and 250 km (Δχ² < 1 at 234,
  238, 242 and 246 km), because the radius trades off against the centre across the chords.
  So the paper's 244 km is compatible, but its ± 1 km holds only for a geometry fixed in
  advance; the paper fitted the Kiso light curve alone. The random multistarts of an
  earlier version found different minima in this valley; the example now searches a fixed
  grid of starting radii and centres.
* **Test for an atmosphere**: the profile likelihood gives a 8.1σ detection. The airless
  limb fits the chord durations but not the gradual ingress/egress at Kiso, nor the dip in
  the Fukushima light curve (flux ≈ 0.5 at closest approach, a few km outside the
  surface shadow of the fit).

## (24835) 1995 SM55, 2024 February 25

Seven positive chords from five sites, 31 negative chords; Gaia DR3 224489389287197184 at the
occultation epoch. Shadow velocity 10.42 km/s (10.42), stellar diameter 0.63 km (0.6).

| Fit | a′ [km] | b′ [km] | φ′ [°] | D_eq,A [km] | χ²/dof |
|---|---|---|---|---|---|
| **Paper** (Table 3, N = 14) | 104.3 ± 0.4 | 83.5 ± 0.5 | 44.1 ± 0.4 | 186.7 ± 1.8 | 1.55 |
| Radial residuals | 102.9 ± 1.4 | 82.2 ± 1.1 | 41.0 ± 1.8 | 183.9 ± 2.4 | 39.3 |
| Along the chords | 105.4 ± 1.3 | 84.9 ± 0.9 | 44.9 ± 1.5 | 189.2 ± 2.1 | 46.3 |
| Radial, σ_limb = 13.5 km | 106 ± 28 | 88 ± 8 | – | 193 ± 26 | 1.00 |
| Posterior (σ_limb) | | | | 182.9 −11.1 +10.7 | |

Findings:

* The ellipse is reproduced, but **not its χ²**. The two Wise telescopes stand at the same
  coordinates (Table A.1). Their immersions (18:17:41.091 ± 0.014 s and 18:17:40.344 ±
  0.035 s) differ by 0.75 s, a 20σ inconsistency (normalised residuals −6.5 and +16.9). With
  the published errors no ellipse reaches χ²/dof ≈ 1.5. The quoted errors (± 0.4 km in a′)
  are therefore far too small.
* The Table 3 errors are also internally inconsistent:
  * ε′ = 0.200 ± 0.002 cannot follow from a′ ± 0.4 and b′ ± 0.5 km, which give ± 0.006;
  * D_eq,A = 186.7 ± 1.8 km does not follow from them either (± 0.7 km);
  * a′ is ± 0.4 km in Table 3 and ± 0.3 km in Sect. 4.1.
* The text (Sect. 2.2) says all positive detections used exposures ≥ 2 s. Table A.1 lists
  0.1 s (RTT150, Chalin) and 0.3 s (T100).
* **Astrometry** (timing errors only): RA 03 40 30.534398, Dec +39 08 06.70796, against the
  published 03 40 30.5344080, +39 08 06.708383. That is a difference of 0.1 and 0.4 mas.
* **Albedo**: 0.82 −0.09 +0.11 (paper 0.80 ± 0.04). The uncertainty here includes the Wise
  inconsistency, through σ_limb.
* **3-D shape** (amplitude 0.05 ± 0.02 mag, rotational phase 0.09, isotropic pole):
  posterior a = 102.6 −6.5 +9.7, b = 92.2 −6.9 +5.4, c = 76.5 −27.9 +11.4 km,
  D_eq,V = 176.9 −18.6 +12.2 km. The paper gives 106 ± 2, 92 ± 12, 76 ± 10 km and 181 ± 12 km
  from a grid search with conditional errors.

## (38628) Huya, 2023 June 24

Eleven positive chords, fourteen negative chords (including La Palma, which recorded only the
satellite); Gaia DR3 4360429542557512064. Shadow velocity 22.43 km/s (22.43).

| Fit | a′ [km] | b′ [km] | PA [°] | R_eq [km] | χ²/dof |
|---|---|---|---|---|---|
| **Paper**, 2023 Free | 222.5 ± 9.1 | 198.7 ± 15.2 | 31.5 ± 15.5 | 210.3 ± 11.0 | 0.68 |
| **Paper**, 2023 Restrict (PA 53.7 ± 2.2°) | 218.7 ± 8.1 | 200.3 ± 14.9 | 53.7 | 209.3 ± 10.3 | 0.75 |
| Free, negative chords ignored | 224.3 ± 12.1 | 198.3 ± 8.8 | 45 ± 28 | 210.9 ± 7.0 | 0.75 |
| Free, negative chords hard | 223.2 ± 10.0 | 185.8 ± 11.5 | 12 ± 13 | 203.6 ± 5.7 | 0.95 |
| PA prior 53.7 ± 2.2°, negatives hard | 247.8 ± 11.0 | 174.9 ± 8.7 | 26 ± 2 | 208.2 ± 5.3 | 10.6 |

Findings:

* Without the negative chords, the paper's free solution is reproduced. That limb, however,
  **crosses three negative chords**: Guirguillano by 41 km, San Marcello Pistoiese by 11 km
  and Belesta by 5 km. The paper filtered its solutions with the Belesta negative only.
* With all negative chords as constraints, the limb is smaller: b′ = 186 km and
  R_eq = 204 km. The PA range from the satellite orbit (53.7 ± 2.2°) is then incompatible
  with the data (χ²/dof = 10.6; the fit pulls the PA to 26°).
* Possible explanations:
  * the Guirguillano data (cycle 10 s, GPS) are faulty, which the paper does not say;
  * topography;
  * the satellite does not orbit in Huya's equator.
* **Astrometry**: within 0.2 and 0.5 mas of Table 4.
* **Satellite**:
  * chord length 182.9 ± 9.5 km (paper 179 ± 14 km);
  * chord midpoint at Δα cos δ = +0.0653″, Δδ = −0.0439″ from Huya's centre, 1599 km. The
    paper gives +0.06546″, −0.04440″ and 1603 km for the satellite centre.
* **Huya and the satellite fitted together** ([`huya_2023_satellite.jl`](../generated/huya_2023_satellite.md);
  `BinaryLimb`, the La Palma contacts assigned to the satellite, its radius fixed at 106.5 km as
  in the paper):
  * two equally good solutions (χ² = 16.16 both), the satellite on either side of the chord;
    the limb of Huya is that of the fit with all negatives (223.2 × 185.8 km);
  * satellite relative to Huya: Δα cos δ = +65.3 ± 0.5 mas, Δδ = −46.5 or −41.2 ± 0.6 mas,
    separation 1631 or 1570 ± 11 km (paper, mean of both: +65.46 ± 3.00, −44.40 ± 5.00 mas,
    1603 km);
  * the satellite positions agree with Table 4 within 0.05 mas, **with the labels interchanged**:
    the solution with the satellite south of the chord (smaller declination,
    −07 00 20.16539 here) is the one labelled "Northern" in Table 4 (−07 00 20.16537);
  * the radius: χ² is flat from half the chord (91 km) up to at least 200 km: the chord gives
    the lower limit D ≥ 183 km and no upper limit (no negative chord near La Palma).
* **2019 Restrict**: the paper's errors are from the timing errors alone (± 0.11 km in a′),
  although χ²/dof ≈ 210 (194 here). Scaled by √(χ²/dof), they are ≈ 15 times larger.
* The caption of Fig. 9 dates the 2019 event to "2019 March 28"; it was March 18.

## (28978) Ixion, six events 2020–2023

Six multi-chord events (23 positive chords with 55 contact times, 13 negative chords) with a
shared ellipse and a centre per event. The distances and velocities of Table 3 are
reproduced to 0.01 au and 0.01 km/s.

| Global fit | a [km] | b [km] | PA of the minor axis [°] | R_eq [km] | χ²ν |
|---|---|---|---|---|---|
| **Paper** (Eq. 4, σ_model = 4 km) | 363.42 +3.53 −3.85 | 333.98 +7.07 −4.96 | 110.97 − 90 = 20.97 +7.00 −6.54 | 348.39 +5.37 −4.43 | 1.17 (N = 56, ν = 41) |
| σ_limb = 4 km | 362.58 ± 2.55 | 337.46 ± 4.52 | 21.78 ± 5.81 | 349.80 ± 2.81 | 1.10 (N = 55, ν = 40) |
| Posterior | 362.7 −2.4 +2.6 | 338.1 −4.7 +4.3 | | 350.1 −2.6 +2.7 | |

Geometric albedo 0.104 ± 0.002 (paper 0.106 ± 0.003).

Findings:

* The global solution is reproduced within its errors, and so are the per-event fits (with
  σ_model = 4 km) where they are constrained.
* **Position-angle conventions**:
  * the per-event ϕ of Table 4 is the position angle of the minor axis (SORA's convention;
    e.g. 18° vs 15°, 38° vs 35°, 27° vs 31°);
  * the global ϕ = 110.97° of Eq. (4) is that of the major axis (= 20.97° for the minor);
  * the paper uses both conventions without saying so.
* Table 4, 28 April 2021: R_equiv = 386.3 km and ε′ = 0.231 do not follow from a = 391.9 and
  b = 322.2 km, which give 355.3 km and 0.178. The other columns are consistent.
* The abstract speaks of five multi-chord events. The global fit (N = 56, ν = 41,
  M = 3 + 2N_events) implies six. Table A.1 lists 55 contact times for the six events.
* **Astrometry** (Table 5):
  * 2020-10-13 agrees to 1.5 mas and 2020-08-17 to 8 mas;
  * the other four events differ by 30–51 mas, i.e. 800–1400 km across the track, at any
    epoch near the published one;
  * the cause is unclear: Table 2 star positions, Table 5 epochs or positions, or the
    reduction here.
* The station table's continuations are labelled "Table D.2" and "Table D.3"; they are
  Table D.1.

## (229762) 2007 UK126, 2014 November 15

Three positive chords (GPS-timed), no negative chords tabulated. The star is UCAC4 448-006503
in the paper; its Gaia DR3 counterpart 3206533286529807488 is used here. Distance
42.573 au (42.572); the chord lengths of Table 2 are reproduced within 0.2 km.

| Fit | 2a [km] | 2b [km] | θ (major axis from E to N) [°] | a/b |
|---|---|---|---|---|
| **Paper** (Table 5, error-in-variables) | 645.80 ± 5.68 | 597.81 ± 12.74 | 21.25 ± 5.65 | 1.080 ± 0.025 |
| Along the chords | 645.9 ± 7.2 | 595.4 ± 36.9 | 20.7 ± 14.6 | 1.085 |
| Radial residuals | 645.8 ± 8.2 | 595.3 ± 47.0 | 20.5 ± 19.2 | 1.085 |
| Posterior | 652.0 −9.3 +20.6 | 623.0 −39.2 +22.4 | 40 −24 +72 | 1.063 −0.044 +0.059 |

Geometric albedo 14.3 −1.1 +1.2 % (paper 15.0 ± 1.6 %).

Findings:

* The ellipse is reproduced, but **not its uncertainties**.
  * Six contacts and five ellipse parameters leave one degree of freedom. With the
    published timing errors, the least-squares covariance gives errors in 2b and θ about
    three times larger than Table 5.
  * The posterior is strongly non-Gaussian: the minor axis and the orientation are only
    weakly constrained by three nearly parallel chords.
  * The paper's error-in-variables fit [Szpak2015](@citep) evidently underestimates them.
* The paper's "axes" a_El, b_El are **full axes**, not the usual semi-axes. The central chord
  alone is 641 km, and the albedo of Sect. 4.3 follows from a_El · b_El as full axes. The
  term is not defined in the paper.
* The paper's Fig. 4 shows a "planar 68.3 % confidence region" of the ellipse, i.e. a region
  containing the true ellipse with 68.3 % probability. This corresponds to `:region` here,
  not to SORA's pointwise envelope (`examples/output/2007uk126_2014/2007uk126_limb_fit_paper.png`).

## (174567) Varda, 2018 September 10

Five positive and six negative chords; Gaia DR2 4367203805493754752 at the epoch. Distance
46.638 au (46.634), shadow velocity 10.86 km/s (10.86), stellar diameter 0.82 km (0.82).

| Fit | a′ [km] | ε′ | P′ [°] | R′_equiv [km] | χ²/dof |
|---|---|---|---|---|---|
| **Paper** (Table 5, N = 10) | 383 ± 3 | 0.066 ± 0.047 | 67 ± 8 | 370 ± 7 | 1.53 |
| Radial residuals | 383.4 ± 3.6 | 0.071 ± 0.054 | 72.4 ± 4.7 | 369.6 ± 7.7 | 1.49 |
| Along the chords | 383.2 ± 3.4 | 0.070 ± 0.050 | 72.5 ± 4.4 | 369.7 ± 6.9 | 1.69 |
| Posterior | 382.6 −2.4 +4.2 | 0.027 −0.019 +0.047 | – | 377.6 −8.2 +4.1 | |

Findings:

* The limb fit, its χ² and the **astrometry** are reproduced: RA 17 18 25.12973, Dec
  −02 05 14.3055 against 17 18 25.12972, −02 05 14.3055.
* The **chord lengths** of Table 4 are 0.8 % longer than here (670.4 vs 664.8 km, …). They
  equal duration × the geocentric shadow velocity (10.86 km/s), not the topocentric
  velocity at each site. The fit uses the timings, so it is unaffected.
* **Albedo**: p_V = 0.099 ± 0.002 in Table 5, but R′_equiv ± 7 km alone gives ± 0.004. Here:
  0.095 −0.002 +0.004.
* The posterior of ε′ piles up towards 0, a nearly circular limb. Its position angle is then
  undetermined. The least-squares ε′ = 0.07 ± 0.05 hides this.

## (84522) 2002 TC302, 2018 January 28

Twelve positive chords in Italy, France, Switzerland and Slovenia and four near misses
[Ortiz2020](@citep). Distance 44.48 au, shadow velocity 4.77 km/s (4.77).

| Fit | a′ [km] | b′ [km] | PA [°] | D_eq [km] | χ²/dof |
|---|---|---|---|---|---|
| **Paper** (Sect. 5) | 271.6 ± 9.0 | 229.8 ± 5.5 | 3.0 ± 1.0 | 499.6 | |
| Timing errors only | 260.0 ± 1.4 | 228.5 ± 1.4 | 12.1 ± 2.0 | 487.6 ± 2.4 | 12.80 |
| Along the chords, σ_limb = 10.9 km | 272.1 ± 11.9 | 232.7 ± 7.3 | 4.2 ± 11.2 | 503.2 ± 17.1 | 1.00 |
| Radial, σ_limb = 10.0 km | 270.8 ± 12.7 | 231.0 ± 9.0 | 7.6 ± 13.8 | 500.2 ± 19.4 | 1.00 |
| Posterior (σ_limb) | 271.0 −8.4 +9.1 | 232.2 −6.2 +5.8 | 4.4 −7.2 +9.5 | 502.9 −13.6 +11.2 | |

Findings:

* The ellipse is reproduced once the scatter of the chords (χ²/dof = 12.8 with the timing
  errors alone) is absorbed by a limb allowance of ≈ 11 km, as topography or timing problems.
* The **position angle** is not determined to ± 1°: with σ_limb its error is ± 11°.
* Clock offsets with N(0, 0.5 s) priors single out **Tavolaia** (−4.85 s); the other offsets
  stay within ± 1.5 s.
* The four near misses lie 31–43 km outside the limb.
* The **star** position of Table 1 is the Gaia DR2 position propagated to J2000.0, 0.64″
  from the star at the event (21 000 km at the body); the Gaia DR2 position at the event is
  used here. The Gaia source ID quoted in the abstract does not exist.
* Astrometry: RA 02 21 49.44338, Dec +28 24 13.4122 at 21:53:36.6 UTC (σ 0.5 mas).

## (208996) 2003 AZ84, 2012 February 3 and 2014 November 15

Three positive chords per event, with negative chords [DiasOliveira2017](@citep). Shadow
velocities 24.82 and 9.41 km/s (25.0, 9.4). Each event has one degree of freedom for an ellipse.

| Fit | a′ [km] | ε′ | PA [°] | R_eq [km] |
|---|---|---|---|---|
| **Paper** 2012 | 426 −16 +26 | 0.34 −0.09 +0.10 | 41 −5 +8 | 346 ± 6 |
| Least squares 2012 | 412.4 ± 17.4 | 0.27 ± 0.12 | 44.0 ± 10.5 | 352.0 ± 17.0 |
| Posterior 2012 | 415 −12 +24 | 0.19 −0.12 +0.09 | 67 −23 +21 | 380 −20 +15 |
| **Paper** 2014 | 393 ± 1 | 0.054 ± 0.003 | 36 ± 4 | 382 ± 3 |
| Least squares 2014 | 394.3 ± 0.6 | 0.049 ± 0.004 | 46.5 ± 1.8 | 384.4 ± 0.2 |
| Posterior 2014 | 394 −0.6 +0.7 | 0.049 ± 0.004 | 46.3 −1.9 +1.8 | 384.4 ± 0.2 |

Findings:

* The 2014 ellipse is reproduced, except that its position angle is 10° larger (2.3σ).
* The 2012 posterior is far from Gaussian: the area-equivalent radius is 380 −20 +15 km, not
  346 ± 6 km.
* The paper's **change of the apparent oblateness** between the events ("more than 3σ") is not
  significant with these posteriors: ε′(2012) ≤ ε′(2014) in 10 % of the draws.
* **Astrometry**: in 2014 the body position agrees with the paper's corrected position to
  2 mas. In 2012 the declination agrees, but the right ascension differs by −425 mas, twice
  the paper's RA offset f_c (6800 km = 212 mas): a **sign error** in the paper's 2012 offset.
* **Single-chord events** of 2011 January 8 and 2013 December 2: the paper does not identify
  the stars. Each is the only Gaia DR3 source within 5″ of the predicted position
  (3151679712277335168, G = 18.4; 3147373651082806272, G = 14.9). With them, the shadow
  velocities (26.00 and 16.51 km/s; paper 26 and about 17) and chord lengths (572.8 and
  669.7 km; paper 573 ± 21 and 657 ± 35) are reproduced. Each chord gives two astrometric
  positions (centre on either side of the chord, circle of the 2014 size), 2–25 mas from the
  current JPL orbit.
* Several station longitudes are given in mixed conventions in the paper and were corrected.

## (84922) 2003 VS2, 2019 October 22

Eleven positive and thirteen negative chords in south-eastern Europe [VaraLubiano2022](@citep).
Shadow velocity 13.01 km/s (13.01), Fresnel scale 1.33 km (1.33).

| Fit | a′ [km] | b′ [km] | PA [°] | D_eq [km] | χ²/dof |
|---|---|---|---|---|---|
| **Paper**, chord 8 shifted | 292 ± 3 | 231 ± 6 | 11 ± 2 | 519 ± 10 | |
| **Paper**, chord centres aligned | 293 ± 3 | 230 ± 6 | 6 ± 2 | 519 ± 10 | |
| Chord 8 clock offset free | 291.2 ± 2.8 | 231.8 ± 7.0 | 7.0 ± 1.1 | 519.6 ± 5.5 | 6.37 |
| All clock offsets, N(0, 1 s) | 293.6 ± 3.0 | 225.9 ± 6.5 | 9.1 ± 3.4 | 515.1 ± 4.9 | 2.98 |
| Chord 8 free, σ_limb = 8.3 km | 292.4 ± 4.5 | 225.9 ± 9.6 | 8.3 ± 1.9 | 514.0 ± 8.1 | 1.00 |

Findings:

* The fitted clock offset of chord 8 is −14.02 ± 0.53 s; the paper shifted it by −14 s. The
  offsets of the other chords follow the paper's alignment shifts within their errors.
* Area-equivalent diameter (posterior) 516.4 −6.9 +6.3 km (519 ± 10), geometric albedo
  0.133 ± 0.009 (0.134 ± 0.010).
* **3-D shape** from the occultation, the light-curve amplitude and the rotational phase:
  least squares a, b, c = 356, 227, 223 km, D_eq,V = 524 km, aspect 54° (paper 339, 235, 226,
  524 km, 59°). The posterior is much wider (a = 365 −17 +33, b = 266 −25 +14 km): one event
  does not constrain the 3-D shape as tightly as the paper's errors suggest.
* Astrometry agrees with the current JPL orbit to 0.2 mas; the paper's offsets refer to JPL#33.
* The star's declination in the paper is misprinted (01.0748″ for 01.748″).

## (95626) 2002 GZ32, 2017 May 20

Five chords used for the limb, a sixth (La Hita 40 cm) with an unknown time offset, and two
close negative chords [SantosSanz2021](@citep). Shadow velocity 20.95 km/s (21.19).

| Fit | a′ [km] | b′ [km] | PA [°] | D_eq [km] | χ² |
|---|---|---|---|---|---|
| **Paper**, original chords | 146 ± 12 | 73.0 ± 4.6 | 318 ± 3 | 206 ± 15 | 29.5 |
| **Paper**, chords aligned | 152.6 ± 8.4 | 73.0 ± 4.2 | 322 ± 3 | 211 ± 12 | 4.0 |
| Radial residuals | 184.3 ± 24.8 | 40.4 ± 1.9 | 163.7 ± 1.0 | 172.5 ± 8.4 | 41.6 |
| Along the chords | 163.7 ± 22.9 | 71.1 ± 4.1 | 139.1 ± 4.7 | 215.7 ± 20.6 | 42.9 |
| Clock offsets, N(0, 0.5 s) | 160.7 ± 29.7 | 67.9 ± 5.9 | 140.1 ± 6.1 | 208.9 ± 25.3 | 4.1 |
| Along the chords, σ_limb = 8.3 km | 151.8 ± 20.8 | 70.7 ± 5.0 | 142.1 ± 6.1 | 207.2 ± 18.2 | 5.0 |

The paper's P′ of 318–322° is the same orientation as 138–142° here.

Findings:

* With the timing errors only, **radial residuals** run into a different, flatter solution
  (b′ = 40 km, PA 164°): for a very elongated limb the radial distance is a poor measure of
  the misfit. Residuals along the chords give the paper's shape.
* Area-equivalent diameter (posterior) 205.2 −9.5 +8.6 km (211 ± 12), albedo 0.044 ± 0.005
  (0.043 ± 0.007), axis ratio 2.1 ± 0.2 (2.09 ± 0.24).
* The negative chord at Sant Esteve passes 2.5 km from the limb and bounds the major axis.
* The fitted clock offsets do not reproduce the paper's +2 s shift of the Athens chord: with
  its large timing errors it needs no shift.
* The shadow velocity from the current ephemeris is 1 % lower than the paper's.

## (136472) Makemake, 2011 April 23

Eight chords from five sites in Chile and Brazil [Ortiz2012](@citep), including the GPS-timed
ULTRACAM light curves of the NTT. The star (NOMAD 1181-0235723) is Gaia DR3
4010602734012333056; the position in the Supplementary Information (SI), from ground-based
astrometry, is 0.10″ from it, within its stated dispersion. The eight negative stations of SI
Table S1 recorded no data (clouds, technical problems).

| Fit | Result | Paper |
|---|---|---|
| Circle | D = 1429.5 ± 9.5 km, χ²/dof = 0.63 | 1430 ± 9 km, χ²/dof = 1.03 |
| Ellipse, posterior | minor axis 1440 −18 +26 km, axial ratio 1.56 ± 0.3 | 1428 ± 17 km, 1.15 ± 0.17 |

Findings:

* The circular fit is reproduced. The lower χ²/dof suggests that the published timing
  errors are conservative.
* All chords cross the central part of the body: the **long axis of an ellipse is not
  constrained**. The minor axis is well determined, the axial ratio is set by the prior and
  its bound; the paper's 1.15 ± 0.17 from one-dimensional χ² cuts is not a constraint of the
  chords. The SI (Sect. 8) reaches the same conclusion and prefers the nearly circular
  solution.
* Astrometry agrees with the current JPL orbit to 1 mas.

## (136108) Haumea and its ring, 2017 January 21

Twelve chords of the body from ten sites in central Europe, and ring events at eight sites
[Ortiz2017](@citep); Mount Agliale saw only the ring. The ring is fitted to the mid-times of the
ring events (`Contact(:mid, …)`), its width (70 km) and opacity (0.5) fixed at the paper's
values.

| Quantity | Here | Paper |
|---|---|---|
| Limb a′, b′ [km] | 853.7 ± 1.7, 579.6 ± 8.7 | 852 ± 2, 569 ± 13 |
| Limb PA (minor axis) [°] | −77.2 ± 0.8 | −76.3 ± 1.2 |
| Ring radius [km] | 2310 ± 5 | 2287 −45 +75 |
| Ring b′ [km], PA [°] | 552, −74.7 | 541 ± 15, −74.3 ± 1.3 |
| Opening angle \|B\| [°] | 13.84 ± 0.09 | 13.8 ± 0.5 |
| Ring pole, solution A / 1 (antipode) [°] | (285.1, −10.2) | (285.1, −10.6) |
| Ring pole, solution B / 2 (antipode) [°] | (312.5, −18.3) | (312.3, −18.6) |
| 3-D, solution A: a, b, c [km] | 1153 ± 24, 852 ± 7, 503 ± 24 | 1161 ± 30, 852 ± 4, 513 ± 16 |
| D_eq,V [km] | 1582 | 1595 ± 11 |

Findings:

* Body, ring and both ring poles are reproduced. The two mirror planes fit the mid-times with
  the same χ² (73.12); a fit started from one plane can stop in a local minimum, so the example
  refits each from the mirror of the other.
* With the ring plane as the equator, solution A reproduces the 3-D ellipsoid. Solution B gives
  an aspect of 76°, at which a and c are not determined; the occultation alone does not decide
  between the poles.
* The negative chord of the body at Mount Agliale passes 9 km from the limb.
* The geocentric shadow velocity from the current ephemeris is 14.36 km/s, against 13.1 km/s
  in the paper. The sizes agree nevertheless, so the paper's chords must rest on a velocity
  close to 14.4 km/s; its quoted value is unexplained.

## (54598) Bienor, 2019–2023

The 2019 January 11 event [FernandezValenzuela2023](@citep): five chords in Spain and Portugal,
four negative chords. The events of 2022 February 6 (A), 2022 December 26 (B) and 2023
February 14 (C) [Rizos2024](@citep): two, seven and two chords. Shadow velocity 2019: 13.61 km/s
(13.47).

| Fit | a′ [km] | b′ [km] | PA [°] | D_eq [km] | χ² (dof) |
|---|---|---|---|---|---|
| **Paper** 2019, S1 (Manzanares shifted) | 94 ± 27 | 66 ± 5 | 52 ± 13 | 160 ± 20 | |
| **Paper** 2019, S2 (centres aligned) | 97 ± 24 | 64 ± 4 | 61 ± 7 | 160 ± 20 | |
| 2019, timing errors only | 122.6 ± 19.6 | 55.1 ± 3.2 | 77.2 ± 4.2 | 164.4 | 56.3 (5) |
| 2019, Manzanares offset fitted (+2.70 ± 0.40 s) | 102.6 ± 21.3 | 65.1 ± 3.1 | 38.0 ± 9.9 | 163.5 | 11.0 (4) |
| 2019, Manzanares +4.06 s (S1) | 110.2 ± 23.4 | 60.4 ± 2.4 | 14.0 ± 6.8 | 163.2 | 22.5 (5) |
| **Paper** Occ. B | 107.0 ± 3.7 | 51.9 ± 2.1 | 109.3 ± 3.4 | 149 ± 4 | |
| Occ. B, timing errors only | 104.2 ± 2.6 | 51.0 ± 0.6 | 110.3 ± 0.7 | 145.8 | 134.7 (9) |
| Occ. B, σ_limb = 9.4 km | 118.3 ± 11.3 | 46.4 ± 3.7 | 106.5 ± 2.7 | 148.1 | 8.1 (9) |

The paper's position angle of Occ. B (19.3°) is that of the major axis.

| 3-D ellipsoid | a [km] | b [km] | c [km] |
|---|---|---|---|
| **Paper** (Occ. B and the light-curve amplitude) | 127 ± 5 | 55 ± 4 | 45 ± 4 |
| Occ. A, B, C jointly, rotational phases | 136.8 ± 15.7 | 52.5 ± 4.0 | 36.2 ± 4.8 |

Findings:

* The ellipse of Occ. B and the 3-D shape are reproduced. The joint fit uses the pole derived
  from the paper's aspect and position angles (α, δ) = (3.58°, +58.23°), which reproduces those
  of all three events within 0.1°, and the rotational phases of Table 5. Its light-curve
  amplitude (0.19 mag) is below the observed 0.23–0.30 mag.
* The **2019 Manzanares shift** of the paper (−4.06 s "towards the north-east") is +4.06 s
  along the motion; a fitted offset gives +2.70 ± 0.40 s and a lower χ². The 2019 shape and
  orientation depend on how the clock problems are treated.
* The **binary scenario** S3: an ellipse and a circle fit the 2019 chords as well as a single
  ellipse only with the separation nearly along the projected pole (PA 127° against 317°),
  which a synchronous binary at maximum light cannot have. With the separation perpendicular to
  the pole the secondary shrinks to a 14 km appendage (χ² = 3.6): the chords give no support for
  two comparable bodies; the paper's S3 needs the thermal diameter as an extra constraint.
* Occ. C is dated 2023 February 14 in Tables B.1 and B.2, 2024 in Table 5 (a misprint). The
  Van Vleck chord of Occ. A needs a clock offset of −4.6 ± 1.4 s (paper: shifted by 4.0 s).
* Two longitudes of Table 2 of the 2019 paper (Sabadell, TJO) are east longitudes in a column
  labelled "W".

## (10199) Chariklo and its rings, 2017–2020

Seven events with timings of the main body and of the ring detections (Table C.1 at the CDS)
[Morgado2021](@citep); the rings are fitted to the mid-times of their detections.

| Quantity | Here | Paper |
|---|---|---|
| 2017-07-23: C1R radius [km] | 385.2 ± 0.5 | 385.9 ± 0.4 |
| 2017-07-23: B, P [°] | 50.00 ± 0.13, −60.94 | 49.91 ± 0.11, −60.85 ± 0.11 |
| 2017-07-23: pole [°] | (150.91 ± 0.17, +41.76 ± 0.09) | (151.03 ± 0.14, +41.81 ± 0.07) |
| 2017-07-23: C2R radius [km] | 401.7 ± 0.6 | 399.8 ± 0.6 |
| 2017-06-22: C1R radius [km], pole [°] | 390.7 ± 2.6, (151.76, +42.27) | 386.6 ± 12.3, (149.46, +40.98) |
| Five ring events jointly: C1R radius [km], pole [°] | 386.9 ± 1.7, (151.13, +41.94) | 385.9 ± 0.4, (151.03, +41.81) |

Findings:

* The ring geometry of each event and the pole are reproduced (errors scaled by √(χ²/dof);
  the rings are not perfectly circular and of variable width, the mid-times have errors of
  a few ms).
* A single event cannot distinguish the ring plane from its mirror image (same χ²); fitting
  the five events with ring detections together rejects the mirror plane by a large margin.
* The astrometry of all seven events agrees with the paper's Table 5 within 1 mas (2019
  August 8 within 3 mas).

## (307261) 2002 MS4, 2020 August 8

59 positive and 40 negative chords [Rommel2023](@citep); stations from the CDS tables.

| Quantity | Here | Paper |
|---|---|---|
| Global ellipse (13 chords, σ_limb = 7 km): a′ [km] | 415.5 ± 2.1 | 412 ± 10 (3σ) |
| ε′ | 0.080 ± 0.006 | 0.066 ± 0.034 (3σ; 0.034 in Sect. 3.4) |
| PA [°] | 126.4 ± 2.7 | 121 ± 16 (3σ) |
| R_equiv [km] | 398.5 ± 1.1 | 398 ± 12 (3σ) |
| Deepest limb point, PA 30–60° [km] | −58 (median of the bin −31) | depression 45.1 ± 1.5 deep, 322 ± 39 km long |
| Astrometry − Table 5 [mas] | 0.0, 0.0 | – |

Findings:

* The global ellipse of the 13 selected chords and the astrometry are reproduced. With σ_limb
  fitted (12.5 km) and the negative chords as soft constraints, the limb crosses the grazing
  Montsec chord by 5 km, within the paper's 7 km tolerance.
* The deviations of the limb points of all 59 chords from the global ellipse show the
  paper's large depression between position angles of about 30° and 60°; the Khmelnytskyi
  ingress (+277 km) is an outlier, as in the paper. Elsewhere the median deviation stays
  within about 10 km.
* The paper gives ε′ = 0.066 in Table 4 but 0.034 in Sect. 3.4 for the same ellipse.

## (50000) Quaoar and its rings, 2011 and 2022

The 2011 May 4 event [BragaRibas2013](@citep) and the 2022 August 9 event with both rings
[Pereira2023](@citep). [Arimatsu2019](@citet) give only the duration of their single chord
of 2019 and are not included.

| Quantity | Here | Paper |
|---|---|---|
| 2022 limb, σ_limb = 6.6 km: a′ [km], ε′, PA [°] | 579.4 ± 4.3, 0.124 ± 0.014, 165.4 ± 1.2 | 579.5 ± 4.0, 0.12 ± 0.01, 345.2 ± 1.2 |
| 2022 R_equiv [km] | 542.4 ± 1.4 | 543 ± 2 |
| 2022 astrometry − Table 2 [mas] | +0.03, −0.16 | – |
| Q1R radius [km] (this event alone) | 4132 ± 30 | 4057 ± 6 (2018–2022) |
| Q1R B, P [°] | −16.0 ± 1.5, 350.5 ± 0.3 | −20.0 ± 0.3, 350.2 ± 0.2 |
| Q2R radius [km] | 2619 ± 54 | 2520 ± 20 |
| 2011 limb, σ_limb = 61.6 km: a′ [km], ε′ | 594 ± 50, 0.10 ± 0.08 | 568–656, 0.14–0.29 |

Findings:

* The 2022 limb and astrometry are reproduced.
* Q1R is far from circular in the 2022 event: the paper's own radial distances are
  3995–4010 km at ingress and 4107–4126 km at egress. With the timing errors alone (3 ms at
  Mauna Kea) a circular ring does not fit (χ²/dof ≈ 4 and a biased opening angle); with an
  allowance `σ_ring` of 18 km for the ring contacts (new) all detections count. The paper's
  pole also uses the events of 2018–2021; from one event the two mirror planes are equally good.
* 2011: the six chords have exposures of 3–14 s and show a large misfit (χ²/dof ≈ 9 for an
  ellipse); the paper explains it by a crater or a mountain. With a limb allowance the
  ellipse is poorly determined.

## (2060) Chiron, 2018 and 2019

[BragaRibas2023](@citep): four chords on 2019 September 8 (three effective; Paris without a
reliable time source), one chord on 2018 November 28. The 2022 event of [Ortiz2023](@citet)
has no tabulated timings.

| Quantity | Here | Paper |
|---|---|---|
| Paris clock offset [s] | +0.99 ± 0.59 | +0.9 |
| 2019 a′ [km] (posterior), ε′ | 119 −12 +21, 0.3 ± 0.2 | 126 ± 22, 0.27 ± 0.18 |
| 2019 R_equiv [km] (posterior) | 102.7 −3.8 +4.0 | 105 −7 +6 (imposed: 98–111) |
| 2018 chord length [km] | 186.9 | 180.5 ± 0.9 |
| Astrometry − Table 7 [mas] | ≤ 0.2 (both events) | – |

Findings:

* The 2019 limb and the Paris offset are reproduced; without the thermal size as a constraint,
  the posterior of the three effective chords still gives an area-equivalent radius close to
  the thermal one.
* The 2018 chord length of the paper is the duration times the geocentric velocity
  (5.77 km/s); at this slow event the station's velocity is 0.2 km/s higher.
* The Boyden negative chord selects the centre north of the SAAO chord, as in the paper.

## (470316) 2007 OC10, 2022 August 22: a double star

The light curves of 13 stations [GomezLimon2025](@citep), published at the CDS: drops of 75 % in
Sicily (GAL Hassin, SLN), of 25 % on Tenerife (Artemis, TAR2) and none at Liverpool on La Palma,
between them. The star is double; the two groups saw the occultation of its two components. As
in the paper, the light curves are fitted directly, with an ellipse and two point-like stars
(`double_star = true`: the flux fraction f₁ of the primary and its offset (x₁, y₁) from the
photocentre; the secondary at −f₁/(1 − f₁)·(x₁, y₁)), here with diffraction, the stellar disc
and the exposure. The starting point is the paper's solution transferred to the JPL ephemeris.
Distance 35.532 au (35.532), shadow velocity 23.62 km/s (23.6).

| Quantity | Here (posterior) | Paper (sampling 1) |
|---|---|---|
| single star: χ²/dof | 2.84 (4101 dof) | – |
| double star: χ²/dof | 1.04 (4098 dof) | – |
| a′ [km] | 215.5 −3.5 +6.2 | 215 −7 +10 |
| b′/a′ | 0.58 −0.13 +0.14 | 0.58 ± 0.16 |
| tilt ϕ (major axis from east towards north) [°] | 17.2 −8.7 +7.6 | 16 −11 +12 |
| D_eq [km] | 329.4 −40.2 +37.2 | 330 −55 +56 |
| f₁ | 0.748 −0.007 +0.006 | 0.748 −0.012 +0.011 |
| x₁, y₁ [km] | −307 −9 +7, −208 ± 7 | −305 −38 +31, −209 −146 +22 |
| separation [mas] | 57.0 −0.5 +0.6 | 57 −11 +4 |
| position angle of the secondary [°] | 56.0 ± 0.8 | 56 −17 +3 |
| Δmag | 1.18 ± 0.04 | 1.18 ± 0.07 |
| 2007 OC10 − Table 5 [mas] (Δα cos δ, Δδ) | −0.4, −0.2 | – |
| primary star − Table 5 [mas] | −0.4, −0.1 | – |
| secondary star − Table 5 [mas] | −0.4, +0.3 | – |

Findings:

* The double star explains the light curves to their noise; a single star cannot (χ²/dof 2.84):
  it cannot produce the two different depths with no drop between them.
* All published values are reproduced: the limb, its area-equivalent diameter, and the flux
  ratio, separation, position angle and magnitude difference of the stars. The positions of the
  body and of the two stars agree with Table 5 to 0.5 mas (at the epoch of the table, the
  predicted closest approach, 10 s after the one computed here).
* The credible intervals of the ellipse are somewhat narrower than the paper's, those of the
  star geometry much narrower (separation ±0.5 against −11 +4 mas; y₁ ±7 against −146 +22 km).
  The paper's intervals have long tails, in y₁ and in the centre (y_c −180 km), that do not
  appear here. The posterior here comes from 32 walkers started at the least-squares solution
  (autocorrelation times ≤ 250 steps, 2500 steps sampled); the paper used 512 walkers started
  from its uniform priors. Whether the tails are real (a second, weaker mode) or an effect of
  the paper's sampling is not resolved here.
* The least-squares covariance is not usable: with exposures of up to 10 s the model is
  piecewise linear in the contact times and its Jacobian by finite differences is unreliable.
  The posterior gives the uncertainties.

## Summary

* Nine of these fits were repeated with SORA on the same data; the two packages agree (see
  [Occultations.jl and SORA](sora.md#The-same-events-reduced-with-both-packages)).

* The package reproduces the published geometry (distances, velocities, stellar diameters,
  Fresnel scales) of the events, and the published limbs within their uncertainties in most
  cases. Exceptions: the 2002 KX14 minor axis and diameter (≈ 2σ larger with a joint
  multi-event fit) and the 2014 position angle of 2003 AZ84 (10° larger).
* Most data sets show more scatter than the timing errors (χ²/dof 8–230 for KX14 and Huya,
  13 for 2002 TC302); limb allowances, clock-offset priors, per-event centres and posterior
  sampling make this explicit instead of hiding it in the fitting method.
* Degeneracies are exposed rather than resolved by a choice of starting point: the mirror
  solution of a two-chord event (KX14 D), the shape–pole valley of a single event with a
  light-curve amplitude (UY117), the mirror ring planes of Haumea, and the unconstrained long
  axes of Makemake and the three-chord ellipses (2003 AZ84, 2007 UK126). Posteriors show
  errors several times those from one-dimensional χ² cuts in these cases.
* Issues found in the publications: the Occ. D timings of KX14 (since corrected by the
  corrigendum), an inconsistent stellar diameter for KX14 event E, a sign error in the
  latitude of Apache Point, the use of the Gaia DR2 star position for the UY117 astrometry, a
  sign error in the 2012 RA offset of 2003 AZ84, the J2000 star position of 2002 TC302 (0.64″
  from the star at the event), and misprinted station or star coordinates (2003 AZ84,
  2003 VS2).
