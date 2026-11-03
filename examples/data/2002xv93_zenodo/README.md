# 2002 XV93 occultation of 2024-01-10: light curves and model profiles

Files from the dataset of Arimatsu et al. (2026), Nature Astronomy, doi:10.1038/s41550-026-02846-1:

K. Arimatsu et al., *Raw movie data, reduced light curves, and model profiles for the
2024-01-10 occultation by (612533) 2002 XV93*, Zenodo,
[doi:10.5281/zenodo.18976431](https://doi.org/10.5281/zenodo.18976431),
licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

Only the reduced light curves (`lightcurve_data/`) and the best-fit model profiles
(`model_profiles/`) are included, unchanged; the raw movies (5.5 GB) are not. The original
README files of these two directories are `README_lightcurve_data.dat` and
`README_model_profiles.dat`.

* `lightcurve_<station>.csv` – `MJD` (UTC), `Scaled_flux` (normalised), `Error` (1σ)
* `bestfitmodel_<gas>_rho<ρ>.csv` – distance to the body centre [km], normalised intensity of
  the best-fit refraction model (CH₄, N₂, CO; bulk density 1500 or 1000 kg/m³)
