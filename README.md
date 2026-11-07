# Nested bilateral filtering of geodetic velocity fields (release 1.0)

MATLAB code accompanying

> Chang, F., Hooper, A. J., Dong, S., Yin, H., Fang, J., & Elliott, J. R.
> Robust Strain-Rate Mapping from Geodetic Velocity Fields Using Refined
> Bilateral Filtering. *Journal of Geophysical Research: Solid Earth*

The code filters gridded InSAR (or InSAR + GNSS) velocity fields before
strain-rate estimation. Three filters are compared: a sliding Gaussian
filter, a sliding median filter and the **nested bilateral filter** proposed
in the paper. The nested bilateral filter weights each pixel in the filtering
window by its distance (spatial kernel, sigma_s) and by its velocity
difference (velocity kernel, sigma_v) from a robust reference velocity, the
local median within a smaller center window nested within the filtering
window. All filters are NaN-aware, so data gaps can be handled.

The development version of the code is on GitHub:
https://github.com/fengnian-chang/InSAR_velocity_filter_method

## Contents

| File | Description |
|---|---|
| `bilateral_nanconv.m` | Nested bilateral filter (Equations 1-3 and Text S2). With `center_winsize = 1` it reduces to the conventional bilateral filter. |
| `nanconv.m` | NaN-aware convolution, used for the Gaussian filter (third-party code, see below). |
| `nanmedfilt2.m` | NaN-aware median filter with a rectangular or circular window. |
| `generate_synthetic_insar_velocities.m` | Synthetic fault-parallel velocity field from a 2D screw dislocation model with optional shallow creep (Equations 7-8), noise and data gaps. |
| `generate_correlated_noise_fft_Lc.m` | Spatially correlated noise with a Lorentzian power spectrum (Text S1). |
| `filter_comparsion_synthetic_profile.m` | Synthetic test with a straight fault (Figure 2 and Figure S4). |
| `filter_comparsion_synthetic_fault_step.m` | Synthetic test with a fault step (Figure S5). |
| `filter_comparsion_AHB.m` | Filtering and strain rates for the Alpine-Himalayan Belt (Section 4, Figures 4-10). |
| `haversine.m` | Great-circle distance, used to convert pixel differences to gradients. |
| `loess_nanfit2_plane.m` | NaN-aware local planar (LOESS) smoother. Not used for the results in the paper. |
| `vik.mat` | Diverging colormap used in the figures (*vik*, Crameri, 2018). |
| `AHB/ve_AHB_final.tif`, `AHB/vn_AHB_final.tif` | East and north velocities of the AHB (Elliott et al., 2026), see below. |

## Requirements

- MATLAB R2020a or later (tested with R2023a).
- Image Processing Toolbox (`fspecial`, `padarray`).
- Mapping Toolbox (`readgeoraster`, `geocrop`, `geotiffwrite`), needed only
  for `filter_comparsion_AHB.m`.

No compilation is required. Add the folder to the MATLAB path, or run the
scripts from within it.

## Using the filter on your own data

```matlab
% V: 2-D velocity grid (mm/yr) with NaN for missing pixels, ~1 km spacing
% Nested bilateral filter with the parameters used for most AHB regions
Vf = bilateral_nanconv(V, 15, 3, 91, 40, 'median', 'nonanout');
%                         |   |   |   |      |         |
%                 sigma_s-+   |   |   |      |         +- estimate NaN pixels
%                 sigma_v ----+   |   |      +- reference = local median
%                 window size ----+   +- center window size
```

`sigma_s`, the window size and the center window size are given in pixels;
`sigma_v` is in the units of `V`. The window size is about 6 * sigma_s
(three-sigma rule). Use `'nanout'` to keep NaN at the NaN pixels of the input.
See `help bilateral_nanconv` for all options.

Guidance on parameter choice (Section 3.2 and Figure 3): sigma_s = 15 km
works for a wide range of slip rates and locking depths, whereas sigma_v
should increase with slip rate and locking depth. In the paper, sigma_v =
3 mm/yr was used for most AHB regions and 5 mm/yr for Anatolia, where slip
rates are higher. A center window of 40 km was used throughout.

Strain rates are then computed from the gradients of the filtered velocities,
as in `filter_comparsion_AHB.m`.

## Reproducing the figures

### Synthetic tests (Section 3.1)

Noise for all synthetic tests (Text S1): 1 mm/yr Gaussian white noise plus
0.4 mm/yr spatially correlated noise with a Lorentzian power spectrum
(Lc = 500 km). The random seeds are fixed, so the results are reproducible.
Slip rate 5 mm/yr, locking depth 20 km, 400 x 400 km grid with 1 km spacing.
Filters: 100 km window for all methods; sigma_s = 15 km for the Gaussian and
bilateral filters; sigma_v = 2 mm/yr and a 40 km median center window for the
bilateral filter.

| Figure | Script | Settings in the script |
|---|---|---|
| 2a | `filter_comparsion_synthetic_profile.m` | `creep_rate = 0`, `nan_type = 'blocks'` |
| 2b | `filter_comparsion_synthetic_profile.m` | `creep_rate = 2`, `nan_type = 'blocks'` |
| S4a | `filter_comparsion_synthetic_profile.m` | `creep_rate = 0`, `nan_type = 'distributed'` |
| S4b | `filter_comparsion_synthetic_profile.m` | `creep_rate = 2`, `nan_type = 'distributed'` |
| S5a | `filter_comparsion_synthetic_fault_step.m` | `creep_rate = 0` |
| S5b | `filter_comparsion_synthetic_fault_step.m` | `creep_rate = 2` |

Each script takes about one minute on a desktop PC and plots the results in
MATLAB (velocity maps, dVe/dy maps and swath profiles); the layout differs
from the published figures.

### Alpine-Himalayan Belt (Section 4)

`filter_comparsion_AHB.m` filters the east velocity with the three filters
and writes, for each filter, GeoTIFFs of the east-velocity gradients
(`dvedx`, `dvedy`), the maximum shear strain rate (`Max_shear`, Equation 6),
the dilatation rate (`Dil`, Equation 5, positive = extension) and the
vorticity (`Vort`) in nstrain/yr to the `AHB/` folder. File names contain
the filter parameters, e.g.
`Max_shear_bilateral_win91_sig_s15_sig_v3_median40.tif`.

Parameters used in the paper (set in Section 1 of the script):

| Region | Figures | winsize | sigma_s | sigma_v | center window |
|---|---|---|---|---|---|
| Pamirs, Western and Eastern Tibet, Makran to Hindu Kush | 4-7, 9 | 91 | 15 | 3 | 40 |
| Anatolia | 8 | 91 | 15 | 5 | 40 |

The Gaussian and median filters use the same window size as the bilateral
filter. Window sizes and sigma_s are in pixels (1 pixel = 0.01 deg, ~1 km).

The bilateral and median filters loop over every pixel, so processing the
full AHB grid (3801 x 13101 pixels) takes many hours. To test the script or
process one area, set `region = [lon_min lon_max lat_min lat_max]`; a
2 x 2 deg region takes about 15 s. The maps in the paper were plotted with
GMT from the GeoTIFF outputs.

## AHB velocity data

`AHB/ve_AHB_final.tif` and `AHB/vn_AHB_final.tif` are the Eurasia-fixed east
and north velocities (mm/yr) of Elliott et al. (2026) on a 0.01 deg grid. The
east velocity is mainly constrained by Sentinel-1 InSAR; the north velocity
is derived from smoothed GNSS velocities (VELMAP). Only the east velocity is
filtered; the north velocity is used as provided. Please cite Elliott et al.
(2026) when using these data.

## Notes

- Even window sizes are effectively rounded up to the next odd size,
  because the half-width is `floor(winsize/2)`.
- For the AHB data, pixel-center coordinates are taken from the raster
  reference, and gradients are converted to nstrain/yr with the local
  east-west and north-south pixel spacing (haversine distance).
- Outputs at pixels without an east velocity are set to NaN.

## Third-party code and data

- `nanconv.m` was written by Benjamin Kraus (Copyright (c) 2013) and is
  distributed on MATLAB Central File Exchange under the BSD license. It is
  included unchanged; its license text is reproduced in `LICENSE`.
- `vik.mat` is the *vik* colormap from Scientific colour maps (Crameri, 2018).
- The AHB velocity data are from Elliott et al. (2026).

## References

- Crameri, F. (2018). Scientific colour maps. Zenodo.
  https://doi.org/10.5281/zenodo.1243862
- Elliott, J. R., Fang, J., Lazecký, M., Maghsoudi, Y., Ou, Q., Payne, J. A.,
  et al. (2026). Deformation, strains and velocities for the Alpine Himalayan
  Belt from trans-continental Sentinel-1 InSAR & GNSS. *Remote Sensing of
  Environment*, 338, 115320. https://doi.org/10.1016/j.rse.2026.115320

## License

The code is released under the BSD 3-Clause License, see `LICENSE`. This
license does not cover the third-party code and data listed above, which are
distributed under their own terms.

## Contact

Fengnian Chang (code)(cfn@smail.nju.edu.cn); corresponding author: Shaochun Dong (dsc@nju.edu.cn)
