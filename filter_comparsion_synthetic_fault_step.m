%% FILTER_COMPARSION_SYNTHETIC_FAULT_STEP
% Compare Gaussian, median and nested bilateral filtering of a synthetic
% velocity field with a fault step (Section 3.1, Figure S5).
%
% The fault lies at y = -50 km for x < 0 and at y = 50 km for x >= 0.
% Settings used in the paper:
%   Figure S5a : creep_rate = 0
%   Figure S5b : creep_rate = 2
%
% Noise (Text S1): 1 mm/yr Gaussian white noise + 0.4 mm/yr spatially
% correlated noise with a Lorentzian power spectrum (Lc = 500 km); no data gaps.
% Filters: 100 km window for all methods; sigma_s = 15 km for the Gaussian
% and bilateral filters; sigma_v = 2 mm/yr and a 40 km median center window
% for the bilateral filter (same as Figure 2).
%
% Each row shows the velocity field and the gradients dVe/dy and dVe/dx
% (nstrain/yr). dVe/dx highlights artifacts at the step.
%
% Required: generate_synthetic_insar_velocities.m,
%           generate_correlated_noise_fft_Lc.m, bilateral_nanconv.m,
%           nanconv.m, nanmedfilt2.m, vik.mat
% Runtime: a few minutes (the bilateral filter is a pixel-by-pixel loop).

clear; close all;

%% ---------- 1. Synthetic model ----------
slip_rate            = 5;       % slip rate S (mm/yr)
locking_depth_km     = 20;      % locking depth D1 (km)
creep_rate           = 0;       % shallow creep rate C (mm/yr): 0 (Fig. S5a) or 2 (Fig. S5b)
creep_depth_km       = 5;       % creeping depth D2 (km)
white_noise_std      = 1;       % Gaussian white noise (mm/yr)
correlated_noise_std = 0.4;     % correlated noise, Lc = 500 km (mm/yr)
grid_extent_km       = [-200 200 -200 200];   % [xmin xmax ymin ymax] (km)
grid_res_km          = 1;       % grid spacing (km)

fault_type = "step";
fault_params.x0_km      = 0;    % x position of the step (km)
fault_params.y_left_km  = -50;  % fault position for x < x0 (km)
fault_params.y_right_km =  50;  % fault position for x >= x0 (km)

nan_type = 'none';
nan_params.do_plot = false;     % true: plot the noise components

[V_clean, V_noisy, X, Y] = generate_synthetic_insar_velocities( ...
    slip_rate, locking_depth_km, creep_rate, creep_depth_km, ...
    white_noise_std, correlated_noise_std, grid_extent_km, grid_res_km, ...
    fault_type, fault_params, nan_type, nan_params);

%% ---------- 2. Filtering ----------
spatial_winsize = 100;   % filtering window (pixels = km)
spatial_sigma   = 15;    % sigma_s (pixels = km)
velocity_sigma  = 2;     % sigma_v (mm/yr)
center_winsize  = 40;    % center window for the reference velocity (pixels = km)

gauss_kernel = fspecial('gaussian', [spatial_winsize spatial_winsize], spatial_sigma);

V.Gaussian  = nanconv(V_noisy, gauss_kernel, 'edge', 'nanout');
V.Median    = nanmedfilt2(V_noisy, [spatial_winsize spatial_winsize], 'nanout');
V.Bilateral = bilateral_nanconv(V_noisy, spatial_sigma, velocity_sigma, ...
    spatial_winsize, center_winsize, 'median', 'nanout');

%% ---------- 3. Velocity gradients (nstrain/yr) ----------
% (mm/yr)/km = 1e-6 /yr = 1e3 nstrain/yr
[dVc_dx, dVc_dy] = gradient(V_clean, grid_res_km);
dVc_dx = dVc_dx * 1e3;
dVc_dy = dVc_dy * 1e3;

methods = ["Median", "Gaussian", "Bilateral"];
for name = methods
    [gx, gy] = gradient(V.(name), grid_res_km);
    dVdx.(name) = gx * 1e3;
    dVdy.(name) = gy * 1e3;
end

%% ---------- 4. Plot ----------
vik = importdata('vik.mat');

figure('Color', 'w', 'Position', [50 50 900 1000]);
set(gcf, 'DefaultAxesFontName', 'Calibri', 'DefaultTextFontName', 'Calibri');

% Row 1: noisy velocity and noise-free gradients
subplot(4,3,1)
map_plot(X, Y, V_noisy, [-4 4], vik);
title("Synthetic Ve with noise"); colorbar;

subplot(4,3,2)
map_plot(X, Y, dVc_dy, [-100 100], vik);
title("\partialVe/\partialy (noise-free)"); colorbar;

subplot(4,3,3)
map_plot(X, Y, dVc_dx, [-100 100], vik);
title("\partialVe/\partialx (noise-free)"); colorbar;

% Rows 2-4: median, Gaussian and bilateral filtering
row_titles = [sprintf("Median (%d\\times%d)", spatial_winsize, spatial_winsize), ...
              sprintf("Gaussian (%d\\times%d, \\sigma_s=%d)", spatial_winsize, spatial_winsize, spatial_sigma), ...
              sprintf("Bilateral (%d\\times%d, \\sigma_s=%d, \\sigma_v=%d)", spatial_winsize, spatial_winsize, spatial_sigma, velocity_sigma)];

for m = 1:numel(methods)
    name = methods(m);
    base = 3*m;

    subplot(4,3,base+1)
    map_plot(X, Y, V.(name), [-4 4], vik);
    title(row_titles(m));

    subplot(4,3,base+2)
    map_plot(X, Y, dVdy.(name), [-100 100], vik);
    title(sprintf("\\partialVe/\\partialy %s", name));

    subplot(4,3,base+3)
    map_plot(X, Y, dVdx.(name), [-100 100], vik);
    title(sprintf("\\partialVe/\\partialx %s", name));
end

%% ---------- Local functions ----------
function map_plot(X, Y, Z, clim, cmap)
    h = imagesc(X(1,:), Y(:,1), Z, clim);
    axis image xy;
    set(gca, 'Color', [0.7 0.7 0.7]);   % NaN areas in gray
    set(h, 'AlphaData', ~isnan(Z));
    colormap(gca, cmap);
    hold on;
end
