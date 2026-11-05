# Examples

The tutorial and the synthetic examples simulate occultations with a known truth. Each of the
other examples reduces a published occultation again, from the data given in the paper, and
compares the results with the publication. The transcribed data and the stored JPL Horizons
ephemerides are in [`examples/data/`](https://github.com/mkretlow/Occultations.jl/tree/main/examples/data),
so all examples run offline. Each writes its figures and its printed results to
`examples/output/<example>/`.

The examples are grouped as in the paper: the synthetic events first, then the published
events by the kind of analysis – first the cases the package reproduces, then those with
differences from the publication (discussed in
[Comparison with published reductions](comparison/publications.md)).

The [quick-start notebook](https://github.com/mkretlow/Occultations.jl/blob/main/examples/quickstart.ipynb)
(Jupyter; Varda 2018, [Souami2020](@citet)) shows the basic usage in a few commands: prediction,
chords, limb fit, astrometry, posterior and maps.

## Synthetic events

*Known truth; the stations, their equipment and timing errors are simulated.*

| Example | Event | Publication | Shows |
|---|---|---|---|
| [Tutorial](generated/tutorial.md) | synthetic | – | all features, offline |
| [Diffraction and timing](generated/diffraction_timing.md) | synthetic | – | Fresnel fringes, passband/star/exposure smoothing and a timing fit |
| [A satellite and a binary](generated/sim_binary.md) | synthetic | – | `BinaryLimb`: satellite as circle or ellipse, two drops per chord, barycentre, coverage |
| [A contact binary](generated/sim_contact_binary.md) | synthetic | – | `ContactBinaryLimb`: light curves with two drops, visible segments, joint light-curve fit |
| [A contact binary in 3-D](generated/sim_contact_binary_3d.md) | synthetic, 4 events | – | `ContactBinaryEllipsoids`: two rotating ellipsoids, light-curve amplitude, unknown pole |
| [Two narrow rings](generated/sim_rings.md) | synthetic | – | `RingedLimb`: ring dips in light curves, their timing, ring pole, widths, optical depths |
| [A distant ring, three events](generated/sim_ring_pole.md) | synthetic, 3 events | – | `RingedLimb`: the ring pole from several events, the mirror ambiguity |

## Reproduced: one event, ellipse or circle

*Limb allowance, clock offsets, posterior.*

| Example | Event | Publication | Shows |
|---|---|---|---|
| [Makemake](generated/makemake_2011.md) | (136472) Makemake, 2011-04-23 | [Ortiz2012](@citet) | star from Gaia, circle and ellipse, an unconstrained long axis |
| [Varda](generated/varda_2018.md) | (174567) Varda, 2018-09-10 | [Souami2020](@citet) | limb and astrometry |
| [Huya 2019](generated/huya_2019.md) | (38628) Huya, 2019-03-18 | [SantosSanz2022](@citet) | σ_limb, clock priors, posterior |
| [2003 VS2](generated/2003vs2_2019.md) | (84922) 2003 VS2, 2019-10-22 | [VaraLubiano2022](@citet) | a clock off by a minute, 3-D ellipsoid with amplitude and phase |

## Reproduced: a large campaign

*Many chords, selected chords, topography.*

| Example | Event | Publication | Shows |
|---|---|---|---|
| [2002 MS4](generated/2002ms4_2020.md) | (307261) 2002 MS4, 2020-08-08 | [Rommel2023](@citet) | 59 chords, the global ellipse, a large depression |

## Reproduced: several events

*One shape with a centre per event, clock offsets, a triaxial ellipsoid from several events.*

| Example | Event | Publication | Shows |
|---|---|---|---|
| [Bienor](generated/bienor_2019_2023.md) | (54598) Bienor, 4 events 2019–2023 | [FernandezValenzuela2023](@citet), [Rizos2024](@citet) | clock offsets, a binary alternative, a triaxial ellipsoid from three events |
| [Ixion](generated/ixion_2020_2023.md) | (28978) Ixion, 6 events | [Kilic2026](@citet) | global fit, astrometry |

## Reproduced: a three-dimensional shape

*Triaxial ellipsoid, rotational amplitude and phase, posterior.*

| Example | Event | Publication | Shows |
|---|---|---|---|
| [2003 UY117](generated/2003uy117_2020.md) | (143707) 2003 UY117, 2020-10-23 | [Kretlow2024](@citet) | 13 negative chords, 3-D ellipsoid with light-curve amplitude and phase |

## Reproduced: rings

*Body and rings together, ring mid-times, ring pole and its mirror.*

| Example | Event | Publication | Shows |
|---|---|---|---|
| [Haumea and its ring](generated/haumea_2017.md) | (136108) Haumea, 2017-01-21 | [Ortiz2017](@citet) | `RingedLimb` from ring mid-times, mirror ring poles, 3-D shape with the ring plane as equator |
| [Chariklo](generated/chariklo_2017_2020.md) | (10199) Chariklo, 7 events 2017–2020 | [Morgado2021](@citet) | two rings from mid-times, the mirror pole resolved by several events, astrometry |
| [Quaoar](generated/quaoar_2011_2022.md) | (50000) Quaoar, 2011 and 2022 | [BragaRibas2013](@citet), [Pereira2023](@citet) | limb, two rings far from the body, ring allowance |

## Reproduced: a double star

*Light curves fitted directly, posterior.*

| Example | Event | Publication | Shows |
|---|---|---|---|
| [2007 OC10](generated/2007oc10_2022.md) | (470316) 2007 OC10, 2022 | [GomezLimon2025](@citet) | a double star: the light curves of 13 stations fitted directly, posterior |

## Reproduced: an atmosphere

*Light curves fitted directly, ray-traced refraction.*

| Example | Event | Publication | Shows |
|---|---|---|---|
| [2002 XV93](generated/2002xv93_2024.md) | (612533) 2002 XV93, 2024-01-10 | [Arimatsu2026](@citet) | light-curve fits, atmosphere detection |

## With differences: one event, few chords

*Error estimates, residuals along the chords, star positions.*

| Example | Event | Publication | Shows |
|---|---|---|---|
| [2007 UK126](generated/2007uk126_2014.md) | (229762) 2007 UK126, 2014-11-15 | [Schindler2017](@citet) | three chords, non-Gaussian posterior |
| [2002 GZ32](generated/2002gz32_2017.md) | (95626) 2002 GZ32, 2017-05-20 | [SantosSanz2021](@citet) | an elongated centaur, aligned chords as clock offsets |
| [2002 TC302](generated/2002tc302_2018.md) | 2002 TC302, 2018-01-28 | [Ortiz2020](@citet) | star position of the paper at J2000.0, σ_limb, clock offsets |
| [Chiron](generated/chiron_2018_2019.md) | (2060) Chiron, 2018 and 2019 | [BragaRibas2023](@citet) | three effective chords, a clock offset, a single chord |
| [1995 SM55](generated/1995sm55_2024.md) | (24835) 1995 SM55, 2024-02-25 | [Ortiz2025](@citet) | albedo, 3-D ellipsoid, inconsistent timings |

## With differences: negative chords and a satellite

*Negative chords as constraints, binary fit.*

| Example | Event | Publication | Shows |
|---|---|---|---|
| [Huya 2023](generated/huya_2023.md) | (38628) Huya, 2023-06-24 | [Rommel2025](@citet) | negative chords, satellite |
| [Huya and its satellite](generated/huya_2023_satellite.md) | (38628) Huya, 2023-06-24 | [Rommel2025](@citet) | `BinaryLimb`: Huya and the satellite fitted together, satellite astrometry |

## With differences: several events

*Joint fit, mirror solutions, astrometry.*

| Example | Event | Publication | Shows |
|---|---|---|---|
| [2003 AZ84](generated/2003az84_2012_2014.md) | (208996) 2003 AZ84, 2011–2014 | [DiasOliveira2017](@citet) | three chords per event, posteriors, astrometry against the paper |
| [2002 KX14](generated/2002kx14_2020_2023.md) | (119951) 2002 KX14, 5 events | [Rizos2025](@citet) | multi-event fit, mirror solution, albedo |

The pages are made from the example scripts with [Literate.jl](https://github.com/fredrikekre/Literate.jl):
the comments of a script are the text of its page. The code is not run when the manual is
built; the output after each section is that of the last run of the script.

## Standard figures

Every example writes the same standard figures
([`examples/standard_figures.jl`](https://github.com/mkretlow/Occultations.jl/blob/main/examples/standard_figures.jl)), next to its own:

| File | Content |
|---|---|
| `<name>_limb_fit.png` | the limb fit with the default display: 1σ shape band and centre ellipse |
| `<name>_limb_fit_region.png` | the limb fit with the 68.3 % (1σ) region of whole limbs |
| `<name>_limb_fit_paper.png` | the limb in the style of the paper's figure (band type and level) |
| `<name>_posterior.png` | the limbs of the posterior |
| `<name>_geomap_prediction.png` | the predicted shadow path (ephemeris only, nominal radius) with its cross-track uncertainty |
| `<name>_geomap_postfit.png` | the shadow path reconstructed from the fit |
| `<name>_geomap_comparison.png` | both paths together |
| `<name>_globe.png` | the post-fit path on the globe: Earth limb, night side shaded, terminator |

The results of all examples are compared with the publications in
[Comparison with published reductions](comparison/publications.md); the fitting methods of
SORA, which most of these papers used, are compared with those of this package in
[Occultations.jl and SORA](comparison/sora.md).
