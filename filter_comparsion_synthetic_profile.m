%% FILTER_COMPARSION_SYNTHETIC_PROFILE
% Compare Gaussian, median and nested bilateral filtering of a synthetic
% velocity field from a 2D screw dislocation model (Section 3.1).
%
% Settings used in the paper:
%   Figure 2a  : creep_rate = 0, nan_type = 'blocks'
%   Figure 2b  : creep_rate = 2, nan_type = 'blocks'
%   Figure S4a : creep_rate = 0, nan_type = 'distributed'
%   Figure S4b : creep_rate = 2, nan_type = 'distributed'
%
% Noise (Text S1): 1 mm/yr Gaussian white noise + 0.4 mm/yr spatially
% correlated noise with a Lorentzian power spectrum (Lc = 500 km).
% Filters: 100 km window for all methods; sigma_s = 15 km for the Gaussian
% and bilateral filters; sigma_v = 2 mm/yr and a 40 km median center window
% for the bilateral filter.
%
% Strain rates are represented by dVe/dy (nstrain/yr), the gradient of the
% east velocity in the fault-perpendicular direction.
%
% Required: generate_synthetic_insar_velocities.m,
%           generate_correlated_noise_fft_Lc.m, bilateral_nanconv.m,
%           nanconv.m, nanmedfilt2.m, vik.mat
% Runtime: a few minutes (the bilateral filter is a pixel-by-pixel loop).

clear; close all;

%% ---------- 1. Synthetic model ----------
slip_rate            = 5;       % slip rate S (mm/yr)
locking_depth_km     = 20;      % locking depth D1 (km)
creep_rate           = 0;       % shallow creep rate C (mm/yr): 0 (Fig. 2a/S4a) or 2 (Fig. 2b/S4b)
creep_depth_km       = 5;       % creeping depth D2 (km)
white_noise_std      = 1;       % Gaussian white noise (mm/yr)
correlated_noise_std = 0.4;     % correlated noise, Lc = 500 km (mm/yr)
grid_extent_km       = [-200 200 -200 200];   % [xmin xmax ymin ymax] (km)
grid_res_km          = 1;       % grid spacing (km)

fault_type   = "straight";      % east-west fault along y = 0
fault_params = struct();

% Data gaps
nan_type = 'blocks';            % 'blocks' (Figure 2) or 'distributed' (Figure S4)
nan_params.do_plot = false;     % true: plot the noise components and NaN mask
switch nan_type
    case 'blocks'
        % One rectangular gap given by row and column indices
        nan_params.blocks = {{220, 260, 260, 340}};
    case 'distributed'
        % 40% of pixels removed at random within x = 0-200 km, y = 15-200 km
        nan_params.bbox_km = [0 200 15 200];
        nan_params.p_in    = 0.40;
        nan_params.p_out   = 0;
end

[V_clean, V_noisy, X, Y] = generate_synthetic_insar_velocities( ...
    slip_rate, locking_depth_km, creep_rate, creep_depth_km, ...
    white_noise_std, correlated_noise_std, grid_extent_km, grid_res_km, ...
    fault_type, fault_params, nan_type, nan_params);

%% ---------- 2. Filtering ----------
spatial_winsize = 100;   % filtering window (pixels = km)
spatial_sigma   = 15;    % sigma_s (pixels = km), ~1/6 of the window
velocity_sigma  = 2;     % sigma_v (mm/yr), ~2x the noise level
center_winsize  = 40;    % center window for the reference velocity (pixels = km)

gauss_kernel = fspecial('gaussian', [spatial_winsize spatial_winsize], spatial_sigma);

V.Gaussian  = nanconv(V_noisy, gauss_kernel, 'edge', 'nanout');
V.Median    = nanmedfilt2(V_noisy, [spatial_winsize spatial_winsize], 'nanout');
V.Bilateral = bilateral_nanconv(V_noisy, spatial_sigma, velocity_sigma, ...
    spatial_winsize, center_winsize, 'median', 'nanout');

%% ---------- 3. Strain rate dVe/dy (nstrain/yr) ----------
% (mm/yr)/km = 1e-6 /yr = 1e3 nstrain/yr
[~, dVc_dy] = gradient(V_clean, grid_res_km);
dVc_dy = dVc_dy * 1e3;

methods = ["Median", "Gaussian", "Bilateral"];
for name = methods
    [~, d] = gradient(V.(name), grid_res_km);
    dV.(name) = d * 1e3;
end

%% ---------- 4. Swath profiles ----------
% Profiles are averaged over 20 km-wide swaths across the fault
xmask1 = X(1,:) >= -110 & X(1,:) <= -90;    % profile 1 (blue)
xmask2 = X(1,:) >=   90 & X(1,:) <= 110;    % profile 2 (red), crosses the data gap
xmask_clean = X(1,:) >= -10 & X(1,:) <= 10; % continuous noise-free velocity profile

prof_y   = Y(:,1);
prof_Vc  = mean(V_clean(:, xmask_clean), 2, "omitnan");
prof_dVc = {mean(dVc_dy(:, xmask1), 2, "omitnan"), mean(dVc_dy(:, xmask2), 2, "omitnan")};
masks    = {xmask1, xmask2};
colors   = {'b', 'r'};

%% ---------- 5. Plot ----------
vik = importdata('vik.mat');
vlim = [floor(min(prof_Vc)) - 0.2, ceil(max(prof_Vc)) + 0.2];

figure('Color', 'w', 'Position', [50 50 1500 900]);
set(gcf, 'DefaultAxesFontName', 'Calibri', 'DefaultTextFontName', 'Calibri');

% Row 1: noisy velocity and noise-free strain rate
subplot(4,6,1)
map_plot(X, Y, V_noisy, [-4 4], vik);
rectangle('Position', [-110,-200,20,400], 'EdgeColor', 'b', 'LineWidth', 1.5);
rectangle('Position', [90,-200,20,400],  'EdgeColor', 'r', 'LineWidth', 1.5);
rectangle('Position', [-150,50,101,101], 'EdgeColor', 'g', 'LineWidth', 1.5);   % filter window size
title("Synthetic Ve with noise"); colorbar;

for p = 1:2
    subplot(4,6,1+p)
    Vpts = V_noisy(:, masks{p});
    ypts = repmat(prof_y, 1, size(Vpts, 2));
    ok = ~isnan(Vpts);
    plot(ypts(ok), Vpts(ok), '.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 4); hold on;
    plot(prof_y, prof_Vc, 'k:', 'LineWidth', 1);
    plot(prof_y, mean(V_noisy(:, masks{p}), 2, "omitnan"), colors{p}, 'LineWidth', 0.8);
    profile_axes([-6 6]);
    title(sprintf("Profile %d Ve (noisy)", p));
end

subplot(4,6,4)
map_plot(X, Y, dVc_dy, [-100 100], vik);
rectangle('Position', [-110,-200,20,400], 'EdgeColor', 'b', 'LineWidth', 1.5);
rectangle('Position', [90,-200,20,400],  'EdgeColor', 'r', 'LineWidth', 1.5);
title("\partialVe/\partialy (noise-free)"); colorbar;

for p = 1:2
    subplot(4,6,4+p)
    plot(prof_y, prof_dVc{p}, colors{p}, 'LineWidth', 0.8);
    profile_axes([0 100]);
    title(sprintf("Profile %d \\partialVe/\\partialy (noise-free)", p));
end

% Rows 2-4: median, Gaussian and bilateral filtering
row_titles = [sprintf("Median (%d\\times%d)", spatial_winsize, spatial_winsize), ...
              sprintf("Gaussian (%d\\times%d, \\sigma_s=%d)", spatial_winsize, spatial_winsize, spatial_sigma), ...
              sprintf("Bilateral (%d\\times%d, \\sigma_s=%d, \\sigma_v=%d)", spatial_winsize, spatial_winsize, spatial_sigma, velocity_sigma)];

for m = 1:numel(methods)
    name = methods(m);
    base = 6*m;

    subplot(4,6,base+1)
    map_plot(X, Y, V.(name), [-4 4], vik);
    title(row_titles(m));

    for p = 1:2
        subplot(4,6,base+1+p)
        plot(prof_y, prof_Vc, 'k:', 'LineWidth', 1); hold on;
        plot(prof_y, mean(V.(name)(:, masks{p}), 2, "omitnan"), colors{p}, 'LineWidth', 0.8);
        profile_axes(vlim);
        title(sprintf("Profile %d Ve %s", p, name));
    end

    subplot(4,6,base+4)
    map_plot(X, Y, dV.(name), [-100 100], vik);
    title(sprintf("\\partialVe/\\partialy %s", name));

    for p = 1:2
        subplot(4,6,base+4+p)
        plot(prof_y, prof_dVc{p}, 'k:', 'LineWidth', 1); hold on;
        plot(prof_y, mean(dV.(name)(:, masks{p}), 2, "omitnan"), colors{p}, 'LineWidth', 0.8);
        profile_axes([0 100]);
        title(sprintf("Profile %d \\partialVe/\\partialy %s", p, name));
    end
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

function profile_axes(ylims)
    grid on;
    xlim([-200 200]);
    ylim(ylims);
    xlabel('y (km)');
end
