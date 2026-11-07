%% FILTER_COMPARSION_AHB
% Apply Gaussian, median and nested bilateral filters to the east velocity
% field of the Alpine-Himalayan Belt (AHB; Elliott et al., 2026) and compute
% horizontal strain rates (Section 4).
%
% Inputs (folder AHB/, Eurasia-fixed velocities in mm/yr on a 0.01 deg grid):
%   ve_AHB_final.tif - east velocity, mainly constrained by InSAR
%   vn_AHB_final.tif - north velocity from smoothed GNSS (VELMAP)
% Only the east velocity is filtered; the north velocity is used as provided.
% Because the grid spacing is ~1 km, window sizes and sigma_s in pixels are
% approximately in km.
%
% Filter parameters used in the paper (center_winsize = 40 for all results;
% Gaussian and median filters use the same winsize as the bilateral filter):
%   Pamirs, Western and Eastern Tibet, Makran-Hindu Kush:
%       winsize = 91, spatial_sigma = 15, velocity_sigma = 3
%   Anatolia:
%       winsize = 91, spatial_sigma = 15, velocity_sigma = 5
% The window size is ~6*sigma_s (three-sigma rule; Section 3.2).
%
% Outputs (GeoTIFF, written to AHB/), for each filter:
%   dvedx, dvedy - east-velocity gradients (nstrain/yr)
%   Max_shear    - maximum shear strain rate (Equation 6)
%   Dil          - horizontal dilatation rate (Equation 5), positive = extension
%   Vort         - vorticity, dVn/dx - dVe/dy (Figures S7-S11)
% Pixels with no east velocity are set to NaN in all outputs.
%
% Required: bilateral_nanconv.m, nanconv.m, nanmedfilt2.m, haversine.m,
%           Mapping Toolbox (readgeoraster, geotiffwrite, geocrop) and
%           Image Processing Toolbox (fspecial).
% Runtime: the bilateral and median filters loop over every pixel and take
% many hours for the full AHB grid; set 'region' to process a subregion.

clear; close all;

%% ---------- 1. Settings ----------
% Filter parameters (pixels; 1 pixel ~ 1 km)
winsize        = 91;    % filtering window (odd, ~6*sigma_s)
spatial_sigma  = 15;    % sigma_s
velocity_sigma = 3;     % sigma_v (mm/yr)
center_winsize = 40;    % center window for the reference velocity

% Region to process: [lon_min lon_max lat_min lat_max] (deg), [] = full AHB grid
% e.g. Pamirs: [66 76 35 41]
region = [];

outdir = 'AHB';
make_plot = true;       % quick-look map of the maximum shear strain rates

%% ---------- 2. Read the velocity fields ----------
[Ve_insar, R_e] = readgeoraster(fullfile('AHB', 've_AHB_final.tif'), 'CoordinateSystemType', 'geographic');
[Vn_gnss,  R_n] = readgeoraster(fullfile('AHB', 'vn_AHB_final.tif'), 'CoordinateSystemType', 'geographic');

if ~isempty(region)
    [Ve_insar, R_e] = geocrop(Ve_insar, R_e, region(3:4), region(1:2));
    [Vn_gnss,  R_n] = geocrop(Vn_gnss,  R_n, region(3:4), region(1:2));
end
assert(isequal(size(Ve_insar), size(Vn_gnss)), 'East and north velocity grids must have the same size.');

% Pixel-center coordinates (deg)
[nrow, ncol] = size(Ve_insar);
[~, xcoords] = intrinsicToGeographic(R_e, 1:ncol, ones(1, ncol));
[ycoords, ~] = intrinsicToGeographic(R_e, ones(1, nrow), 1:nrow);

%% ---------- 3. Filtering of the east velocity (NaN-aware) ----------
gauss_kernel = fspecial('gaussian', [winsize winsize], spatial_sigma);
V_gaussian  = nanconv(Ve_insar, gauss_kernel, 'edge', 'nonanout');
V_median    = nanmedfilt2(Ve_insar, [winsize winsize], 'nonanout');
V_bilateral = bilateral_nanconv(Ve_insar, spatial_sigma, velocity_sigma, ...
    winsize, center_winsize, 'median', 'nonanout');

%% ---------- 4. Velocity gradients (nstrain/yr) ----------
% gradient() returns differences per pixel (mm/yr per pixel). Dividing by the
% pixel spacing in km gives (mm/yr)/km = 1e-6 /yr; the factor 1000 converts
% this to nstrain/yr. The sign makes y point north and x point east (for a
% north-up grid, rows run from north to south and the y scale is negative).
scaley = sign(ycoords(2) - ycoords(1)) * 1000 / haversine(ycoords(1), xcoords(1), ycoords(2), xcoords(1));
scalex = zeros(length(ycoords), 1);
for i = 1:length(ycoords)
    scalex(i) = sign(xcoords(2) - xcoords(1)) * 1000 / haversine(ycoords(i), xcoords(1), ycoords(i), xcoords(2));
end

[dvndx, dvndy] = gradient(Vn_gnss);
s_dvndx = dvndx .* scalex;
s_dvndy = dvndy * scaley;

methods = {'gaussian', 'median', 'bilateral'};
Vf = {V_gaussian, V_median, V_bilateral};

mask_valid = ~isnan(Ve_insar);   % output mask from the original east velocity
for k = 1:numel(methods)
    [dvedx, dvedy] = gradient(Vf{k});
    s_dvedx = dvedx .* scalex;
    s_dvedy = dvedy * scaley;

    %% ---------- 5. Strain-rate tensor and derived quantities ----------
    Exx = s_dvedx;
    Eyy = s_dvndy;
    Exy = 0.5 * (s_dvedy + s_dvndx);

    S.(methods{k}).dvedx     = s_dvedx;
    S.(methods{k}).dvedy     = s_dvedy;
    S.(methods{k}).Max_shear = sqrt(Exy.^2 + (Exx - Eyy).^2 / 4);   % maximum shear strain rate
    S.(methods{k}).Dil       = Exx + Eyy;                           % dilatation rate
    S.(methods{k}).Vort      = s_dvndx - s_dvedy;                   % vorticity

    % Re-apply the data mask of the original east velocity
    fields = fieldnames(S.(methods{k}));
    for f = 1:numel(fields)
        S.(methods{k}).(fields{f})(~mask_valid) = NaN;
    end
end

%% ---------- 6. Write GeoTIFFs ----------
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

tag.gaussian  = sprintf('win%d_sig_s%d', winsize, spatial_sigma);
tag.median    = sprintf('win%d', winsize);
tag.bilateral = sprintf('win%d_sig_s%d_sig_v%d_median%d', winsize, spatial_sigma, velocity_sigma, center_winsize);
if ~isempty(region)
    region_tag = sprintf('_lon%g-%g_lat%g-%g', region);
else
    region_tag = '';
end

quantities = {'dvedx', 'dvedy', 'Max_shear', 'Dil', 'Vort'};
for k = 1:numel(methods)
    for q = 1:numel(quantities)
        fname = sprintf('%s_%s_%s%s.tif', quantities{q}, methods{k}, tag.(methods{k}), region_tag);
        geotiffwrite(fullfile(outdir, fname), S.(methods{k}).(quantities{q}), R_e);
    end
end
disp(['Results saved to ', outdir]);

%% ---------- 7. Quick-look plot ----------
if make_plot
    figure('Color', 'w', 'Position', [50 50 1500 450]);
    titles = {'Gaussian', 'Median', 'Bilateral'};
    for k = 1:numel(methods)
        subplot(1, 3, k);
        Z = S.(methods{k}).Max_shear;
        imagesc(xcoords, ycoords, Z, 'AlphaData', ~isnan(Z));
        axis image xy;
        caxis([0 150]);
        colormap(gca, flipud(hot));
        colorbar;
        title(sprintf('Max shear, %s (nstrain/yr)', titles{k}), 'FontWeight', 'normal');
        xlabel('Longitude (\circ)'); ylabel('Latitude (\circ)');
    end
end
