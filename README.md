This package provides four filtering methods for InSAR velocity fields: bilateral, Gaussian, median, and LOESS (currently a simplified version that fits a low-order polynomial within a moving window). The calling syntax for each method follows the corresponding function or the examples provided here.

### Examples:

```V_bilateral = bilateral_nanconv(V, spatial_sigma, range_sigma, spatial_winsize, center_winsize, center_mode, varargin);```

```V_gaussian = nanconv(V, gauss_kernel, edge_nan_mode);```

```V_median = nanmedfilt2(V, [winsize winsize], nan_mode);```

```V_loess = loess_nanfit2_plane(V, winsize, nan_mode);```

Users can directly run  ```filter_comparison_synthetic_profile.m ``` to compare the performance of different methods on synthetic data. The parameters include the slip rate of an infinitely long strike-slip screw dislocation, locking depth, shallow creep rate, shallow creep depth, and noise (Gaussian white noise and spatially correlated noise). Users can also configure different fault parameters (e.g., fault step) and NaN settings by calling  ```generate_synthetic_insar_velocities.m ```, including NaN blocks and distributed NaNs, to simulate common conditions in real InSAR velocity fields.

For real data applications, different filtering methods can be compared using the script  ```filter_comparison_AHB.m ```. Only the eastward and northward velocity components (Ve.tif and Vn.tif) are required as inputs, and these must be resampled onto a consistent grid with identical spatial resolution and extent. Note that the filtering methods here are pixel-based, and the choice of filtering parameters should be adjusted according to the resolution of the input velocity fields.

### Citing:

Chang, F., Hooper, A. J., Dong, S., Yin, H., Fang, J., & Elliott, J. R. (2026). Robust Strain-Rate Mapping from Geodetic Velocity Fields Using Refined Bilateral Filtering. ESS Open Archive. https://doi.org/10.22541/essoar.15002643/v1

For questions, please contact cfn@smail.nju.edu.cn or xhfk9298@leeds.ac.uk
