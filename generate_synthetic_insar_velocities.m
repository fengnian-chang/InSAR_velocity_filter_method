function [V_clean, V_noisy, X, Y, nan_mask] = generate_synthetic_insar_velocities( ...
    slip_rate, locking_depth_km, creep_rate, creep_depth_km, white_noise_std, ...
    correlated_noise_std, grid_extent_km, grid_res_km, fault_type, fault_params, nan_type, nan_params)
%GENERATE_SYNTHETIC_INSAR_VELOCITIES  Synthetic fault-parallel velocity field with noise and data gaps.
%
%   [V_clean, V_noisy, X, Y, nan_mask] = generate_synthetic_insar_velocities( ...
%       slip_rate, locking_depth_km, creep_rate, creep_depth_km, ...
%       white_noise_std, correlated_noise_std, grid_extent_km, grid_res_km, ...
%       fault_type, fault_params, nan_type, nan_params)
%
%   Generates the east-west (fault-parallel) velocity of an east-west
%   strike-slip fault from a 2D screw dislocation model (Equation 7), with
%   optional shallow creep modeled by the back-slip approach (Equation 8):
%
%       Ve(y) = S/pi * atan(y/D1) - C * [ 1/pi * atan(y/D2) - H(y) ],
%
%   where y is the fault-perpendicular distance, S the slip rate, D1 the
%   locking depth, C the shallow creep rate, D2 the creeping depth and H the
%   Heaviside function. Gaussian white noise and spatially correlated noise
%   with a Lorentzian power spectrum (Lc = 500 km) are added (Supporting
%   Information Text S1), and NaN regions can be imposed to mimic data gaps.
%
%   INPUTS
%     slip_rate            - slip rate below the locking depth, S (mm/yr), e.g. 5
%     locking_depth_km     - locking depth, D1 (km), e.g. 20
%     creep_rate           - shallow creep rate, C (mm/yr), e.g. 0 or 2
%     creep_depth_km       - depth extent of shallow creep, D2 (km), e.g. 5
%     white_noise_std      - standard deviation of Gaussian white noise (mm/yr), e.g. 1
%     correlated_noise_std - standard deviation of correlated noise (mm/yr), e.g. 0.4
%     grid_extent_km       - grid extent [xmin xmax ymin ymax] (km), e.g. [-200 200 -200 200]
%     grid_res_km          - grid spacing (km), e.g. 1
%     fault_type           - 'straight' (default): fault along y = y0
%                            'step'            : fault at y = y_left for x < x0
%                                                and at y = y_right for x >= x0
%     fault_params         - structure with optional fields
%                              y0_km      (straight; default 0)
%                              x0_km      (step; default 0)
%                              y_left_km  (step; default -50)
%                              y_right_km (step; default 50)
%     nan_type             - 'none' (default), 'blocks', 'distributed', or
%                            'both' (distributed + blocks)
%     nan_params           - structure with optional fields
%                              do_plot  : true/false (default false), plot the
%                                         velocities, noise and NaN mask
%                              For 'blocks' / 'both':
%                                blocks : cell array of {r1, r2, c1, c2}
%                                         (row and column index ranges), e.g.
%                                         {{220, 260, 260, 340}}
%                              For 'distributed' / 'both', either
%                                A) bbox mode:
%                                   bbox_km : [xmin xmax ymin ymax] (km)
%                                   p_in    : NaN fraction inside the box (default 0.30)
%                                   p_out   : NaN fraction outside the box (default 0.05)
%                                B) split mode (used when bbox_km is absent):
%                                   y_split_km    : dividing line y = y_split (default 0)
%                                   p_above       : NaN fraction for y > y_split (default 0.40)
%                                   p_below       : NaN fraction for y < y_split (default 0.10)
%                                   include_split : assign y == y_split to the
%                                                   upper side (default false)
%                                seed             : random seed (default 2026)
%                                respect_existing : do not count pixels that are
%                                                   already NaN (default true)
%
%   OUTPUTS
%     V_clean  - noise-free velocity (mm/yr), NaN in the masked pixels
%     V_noisy  - noisy velocity (mm/yr), NaN in the masked pixels
%     X, Y     - coordinate grids (km)
%     nan_mask - logical mask of the NaN pixels
%
%   The noise is reproducible: the white noise uses rng(42) and the
%   correlated noise uses rng(38).
%
%   EXAMPLES
%     % Figure 2a: block data gap, no creep
%     nan_params.blocks = {{220, 260, 260, 340}};
%     [Vc, Vn, X, Y] = generate_synthetic_insar_velocities(5, 20, 0, 5, 1, 0.4, ...
%         [-200 200 -200 200], 1, 'straight', struct(), 'blocks', nan_params);
%
%     % Figure S4a: 40% of pixels removed within x = 0-200 km, y = 15-200 km
%     nan_params = struct('bbox_km', [0 200 15 200], 'p_in', 0.40, 'p_out', 0);
%     [Vc, Vn, X, Y] = generate_synthetic_insar_velocities(5, 20, 0, 5, 1, 0.4, ...
%         [-200 200 -200 200], 1, 'straight', struct(), 'distributed', nan_params);
%
%     % Figure S5a: fault step
%     fp = struct('x0_km', 0, 'y_left_km', -50, 'y_right_km', 50);
%     [Vc, Vn, X, Y] = generate_synthetic_insar_velocities(5, 20, 0, 5, 1, 0.4, ...
%         [-200 200 -200 200], 1, 'step', fp, 'none', struct());
%
%   See also GENERATE_CORRELATED_NOISE_FFT_LC.

    %% ---------- 0. Defaults ----------
    if nargin < 9 || isempty(fault_type)
        fault_type = 'straight';
    end
    if nargin < 10 || isempty(fault_params)
        fault_params = struct();
    end
    if nargin < 11 || isempty(nan_type)
        nan_type = 'none';
    end
    if nargin < 12 || isempty(nan_params)
        nan_params = struct();
    end

    % Correlation length scale of the correlated noise (km), Text S1
    Lc_km = 500;

    %% ---------- 1. Grid ----------
    x_vec = grid_extent_km(1):grid_res_km:grid_extent_km(2);
    y_vec = grid_extent_km(3):grid_res_km:grid_extent_km(4);
    [X, Y] = meshgrid(x_vec, y_vec);

    %% ---------- 2. Screw dislocation with optional shallow creep ----------
    d1 = locking_depth_km * 1000;   % m
    d2 = creep_depth_km   * 1000;   % m
    s  = slip_rate;
    C  = creep_rate;

    % Fault trace y_fault(x) (km)
    fault_type = lower(string(fault_type));
    switch fault_type
        case "straight"
            y0 = getfieldwithdefault(fault_params, 'y0_km', 0);
            y_fault_km = y0 * ones(size(X));
        case "step"
            x0      = getfieldwithdefault(fault_params, 'x0_km', 0);
            y_left  = getfieldwithdefault(fault_params, 'y_left_km', -50);
            y_right = getfieldwithdefault(fault_params, 'y_right_km', 50);
            y_fault_km = y_left * ones(size(X));
            y_fault_km(X >= x0) = y_right;
        otherwise
            error("Unknown fault_type: %s. Use 'straight' or 'step'.", fault_type);
    end

    % Signed fault-perpendicular distance (m) and Heaviside function
    y_perp = (Y - y_fault_km) * 1000;
    H = double(y_perp >= 0);

    V_clean = (s/pi) * atan2(y_perp, d1) ...
              - C * ( (1/pi) * atan2(y_perp, d2) - H );

    %% ---------- 3. Noise (Text S1, Equation S1) ----------
    rng(42);
    white_noise = white_noise_std * randn(size(V_clean));

    rng(38);
    corr_noise = generate_correlated_noise_fft_Lc(X, Y, Lc_km, correlated_noise_std);

    V_noisy = V_clean + white_noise + corr_noise;

    %% ---------- 4. NaN mask ----------
    nan_mask = false(size(V_clean));
    nan_type = lower(string(nan_type));

    switch nan_type
        case "none"
            % no data gaps
        case "blocks"
            apply_blocks();
        case "distributed"
            apply_distributed();
        case "both"
            apply_distributed();
            apply_blocks();
        otherwise
            error("Unknown nan_type: %s. Use 'none', 'blocks', 'distributed', or 'both'.", nan_type);
    end

    V_clean(nan_mask) = NaN;
    V_noisy(nan_mask) = NaN;

    %% ---------- 5. Optional plot ----------
    if getfieldwithdefault(nan_params, 'do_plot', false)
        if exist('vik.mat', 'file')
            cmap = importdata('vik.mat');
        else
            cmap = parula;
        end

        figure('Color', 'w', 'Position', [100 100 1200 700]);
        panels = {V_clean,     [-5 5], 'Noise-free velocity (mm/yr)';
                  V_noisy,     [-5 5], 'Noisy velocity with NaNs (mm/yr)';
                  white_noise, [-3 3], sprintf('White noise (%.2f mm/yr)', white_noise_std);
                  corr_noise,  [-1 1], sprintf('Correlated noise (%.2f mm/yr, L_c = %d km)', correlated_noise_std, Lc_km);
                  white_noise + corr_noise, [-3 3], 'Total noise (mm/yr)';
                  double(nan_mask), [0 1], 'NaN mask'};
        for k = 1:size(panels, 1)
            subplot(2, 3, k);
            h = imagesc(x_vec, y_vec, panels{k, 1}, panels{k, 2});
            axis image xy;
            colormap(gca, cmap);
            set(h, 'AlphaData', ~isnan(panels{k, 1}));
            colorbar;
            title(panels{k, 3}, 'FontWeight', 'normal');
        end
        sgtitle('Synthetic velocity field, noise components and NaN mask');
    end

    %% ---------- Nested functions ----------
    function apply_blocks()
        % Rectangular NaN blocks given by row/column index ranges
        if ~isfield(nan_params, 'blocks') || isempty(nan_params.blocks)
            return;
        end
        blocks = nan_params.blocks;
        if numel(blocks) == 4 && ~iscell(blocks{1})
            blocks = {blocks};
        end
        for kb = 1:numel(blocks)
            b = blocks{kb};
            nan_mask(b{1}:b{2}, b{3}:b{4}) = true;
        end
    end

    function apply_distributed()
        % Randomly distributed NaNs with different densities
        rng(getfieldwithdefault(nan_params, 'seed', 2026));
        respect_existing = getfieldwithdefault(nan_params, 'respect_existing', true);

        if isfield(nan_params, 'bbox_km') && ~isempty(nan_params.bbox_km)
            % Mode A: NaN fraction p_in inside the box and p_out outside it
            bbox = nan_params.bbox_km;
            if numel(bbox) ~= 4
                error('nan_params.bbox_km must be a 1x4 vector [xmin xmax ymin ymax] (km).');
            end
            p_in  = getfieldwithdefault(nan_params, 'p_in', 0.30);
            p_out = getfieldwithdefault(nan_params, 'p_out', 0.05);
            if p_in < 0 || p_in > 1 || p_out < 0 || p_out > 1
                error('NaN fractions p_in and p_out must be within [0, 1].');
            end

            mask_a = (X >= bbox(1)) & (X <= bbox(2)) & (Y >= bbox(3)) & (Y <= bbox(4));
            mask_b = ~mask_a;
            p_a = p_in;
            p_b = p_out;
        else
            % Mode B: NaN fraction p_above for y > y_split and p_below for y < y_split
            y_split = getfieldwithdefault(nan_params, 'y_split_km', 0);
            p_a = getfieldwithdefault(nan_params, 'p_above', 0.40);
            p_b = getfieldwithdefault(nan_params, 'p_below', 0.10);
            if p_a < 0 || p_a > 1 || p_b < 0 || p_b > 1
                error('NaN fractions p_above and p_below must be within [0, 1].');
            end

            if getfieldwithdefault(nan_params, 'include_split', false)
                mask_a = (Y >= y_split);
            else
                mask_a = (Y > y_split);
            end
            mask_b = (Y < y_split);
        end

        if respect_existing
            valid  = ~nan_mask & ~isnan(V_noisy);
            mask_a = mask_a & valid;
            mask_b = mask_b & valid;
        end

        u = rand(size(V_clean));
        nan_mask = nan_mask | (mask_a & (u < p_a)) | (mask_b & (u < p_b));
    end

end

%% ---------- Utility ----------
function v = getfieldwithdefault(s, f, d)
    if isfield(s, f)
        v = s.(f);
    else
        v = d;
    end
end
