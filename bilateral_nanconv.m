function V_filtered = bilateral_nanconv( ...
    V, spatial_sigma, velocity_sigma, ...
    spatial_winsize, center_winsize, ...
    center_mode, varargin)
%BILATERAL_NANCONV  NaN-aware nested bilateral filter for gridded velocity fields.
%
%   V_filtered = bilateral_nanconv(V, spatial_sigma, velocity_sigma, ...
%                                  spatial_winsize, center_winsize, ...
%                                  center_mode, ...)
%
%   For each target pixel x, the filtered velocity is a weighted mean of the
%   original velocities v(i) within a square filtering window N(x)
%   (Equations 1-3 and Supporting Information Text S2, Equations S4-S5):
%
%       v_f(x) = 1/W * sum_i v(i) * Gs(||x - i||) * Gv(|v_ref(x) - v(i)|)
%
%   Gs is a Gaussian spatial kernel with standard deviation spatial_sigma and
%   Gv is a Gaussian velocity kernel with standard deviation velocity_sigma.
%   W is the sum of all weights.
%
%   In the nested bilateral filter, the reference velocity v_ref(x) is not the
%   (possibly noisy) velocity of the central pixel, but a local median (or mean,
%   or Gaussian-weighted average) of the original velocities within a smaller
%   center window Omega_x that is nested within N(x). v_ref is used only to
%   compute the velocity weights; the velocities that enter the weighted mean
%   are the original, unfiltered values. The filter is applied in a single,
%   non-iterative pass. With center_winsize = 1 and center_mode = 'median',
%   the function reduces to the conventional bilateral filter.
%
%   INPUTS
%     V               - 2-D velocity grid (e.g., mm/yr); NaN marks missing data
%     spatial_sigma   - standard deviation of the spatial kernel, sigma_s (pixels)
%     velocity_sigma  - standard deviation of the velocity kernel, sigma_v
%                       (same units as V, e.g., mm/yr)
%     spatial_winsize - size of the square filtering window N(x) (pixels).
%                       The half-width is floor(spatial_winsize/2), so an even
%                       size is effectively rounded up to the next odd size.
%                       In this study the window is ~6*sigma_s.
%     center_winsize  - size of the square center window Omega_x (pixels),
%                       40 (km) in this study
%     center_mode     - how v_ref is computed within Omega_x:
%                         'median'   - local median (used in this study)
%                         'mean'     - local mean
%                         'gaussian' - Gaussian-weighted local average; the
%                                      first extra argument is center_sigma
%     varargin        - for 'gaussian': center_sigma (pixels), followed by
%                       optional flags; for 'median'/'mean': optional flags
%
%   OPTIONAL FLAGS
%     'nanout'   - output is NaN wherever the input V is NaN
%     'nonanout' - (default) output is also estimated at NaN pixels of V where
%                  valid neighbours exist (NaN only where no estimate is possible)
%
%   OUTPUT
%     V_filtered - filtered grid, same size as V
%
%   EXAMPLES
%     % Nested bilateral filter as used in this study (1 km pixels)
%     Vf = bilateral_nanconv(V, 15, 3, 91, 40, 'median', 'nonanout');
%     % Conventional bilateral filter (reference = central pixel)
%     Vf = bilateral_nanconv(V, 15, 3, 91, 1, 'median', 'nanout');
%     % Gaussian-weighted reference velocity (center_sigma = 5 pixels)
%     Vf = bilateral_nanconv(V, 15, 3, 91, 25, 'gaussian', 5, 'nanout');
%
%   See also NANCONV, NANMEDFILT2.

% ---------- Parse center_mode and extra arguments ----------
center_mode = lower(center_mode);

if strcmp(center_mode, 'gaussian')
    % The first extra argument must be center_sigma
    if isempty(varargin) || ~isnumeric(varargin{1})
        error('For center_mode = "gaussian", the first extra argument must be the numeric center_sigma.');
    end
    center_sigma = varargin{1};
    flag_args    = varargin(2:end);
elseif strcmp(center_mode, 'median') || strcmp(center_mode, 'mean')
    flag_args    = varargin;
else
    error('center_mode must be "median", "mean", or "gaussian".');
end

% ---------- Parse flags ----------
nanout = false;
for k = 1:numel(flag_args)
    switch lower(flag_args{k})
        case 'nanout'
            nanout = true;
        case 'nonanout'
            nanout = false;
        otherwise
            error('Unknown flag "%s". Use "nanout" or "nonanout".', flag_args{k});
    end
end

% ---------- Initialization ----------
[m, n] = size(V);
V_filtered = nan(m, n);
valid_mask_global = ~isnan(V);

% ---------- Spatial kernel Gs ----------
half_win_spatial = floor(spatial_winsize/2);
[Xg, Yg] = meshgrid(-half_win_spatial:half_win_spatial, -half_win_spatial:half_win_spatial);
spatial_kernel = exp(-(Xg.^2 + Yg.^2) / (2 * spatial_sigma^2));

% ---------- Reference velocity setup ----------
half_win_center = floor(center_winsize/2);
if strcmp(center_mode, 'gaussian')
    % Pre-compute a NaN-aware Gaussian-smoothed velocity map
    [Xc, Yc] = meshgrid(-half_win_center:half_win_center, -half_win_center:half_win_center);
    center_kernel = exp(-(Xc.^2 + Yc.^2) / (2 * center_sigma^2));
    center_kernel = center_kernel / sum(center_kernel(:));
    V_center = nanconv(V, center_kernel, 'edge', 'nanout');
end

% ---------- Main loop over target pixels ----------
for i = 1:m
    for j = 1:n

        % 1. Filtering window N(x), clipped at the grid edges
        r1 = max(i - half_win_spatial, 1);
        r2 = min(i + half_win_spatial, m);
        c1 = max(j - half_win_spatial, 1);
        c2 = min(j + half_win_spatial, n);

        local_patch = V(r1:r2, c1:c2);
        local_valid = valid_mask_global(r1:r2, c1:c2);

        % Corresponding part of the spatial kernel
        sk = spatial_kernel( ...
            (r1 - i + half_win_spatial + 1):(r2 - i + half_win_spatial + 1), ...
            (c1 - j + half_win_spatial + 1):(c2 - j + half_win_spatial + 1) );

        % 2. Reference velocity v_ref(x) from the center window Omega_x
        if strcmp(center_mode, 'gaussian')
            center_val = V_center(i,j);
            if isnan(center_val)
                continue;   % no reference velocity: output stays NaN
            end
        else
            rr1 = max(i - half_win_center, 1);
            rr2 = min(i + half_win_center, m);
            cc1 = max(j - half_win_center, 1);
            cc2 = min(j + half_win_center, n);

            center_patch = V(rr1:rr2, cc1:cc2);
            center_valid = ~isnan(center_patch);
            if ~any(center_valid(:))
                continue;   % no valid pixel in Omega_x: output stays NaN
            end

            vals = center_patch(center_valid);
            if strcmp(center_mode, 'median')
                center_val = median(vals);
            else
                center_val = mean(vals);
            end
        end

        % 3. Velocity kernel Gv relative to v_ref(x)
        velocity_weight = zeros(size(local_patch));
        velocity_weight(local_valid) = exp( -((local_patch(local_valid) - center_val).^2) / (2 * velocity_sigma^2) );

        % Total weight = spatial kernel x velocity kernel; NaN pixels get zero weight
        total_weight = sk .* velocity_weight;
        total_weight(~local_valid) = 0;

        denom = sum(total_weight(:));
        if denom > 0
            V_filtered(i,j) = sum(local_patch(local_valid) .* total_weight(local_valid)) / denom;
        end

        % 4. 'nanout': keep NaN wherever the input is NaN
        if nanout && isnan(V(i,j))
            V_filtered(i,j) = NaN;
        end
    end
end

end
