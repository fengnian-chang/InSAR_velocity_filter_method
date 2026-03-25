function [V_clean, V_noisy, X, Y, nan_mask] = generate_synthetic_insar_velocities( ...
    slip_rate, locking_depth_km, creep_rate, creep_depth_km, white_noise_std, ...
    correlated_noise_std, grid_extent_km, grid_res_km, fault_type, fault_params, nan_type, nan_params)

% GENERATE_SYNTHETIC_INSAR_VELOCITIES
% Generates synthetic fault-parallel InSAR velocity fields for a 2-D
% strike-slip fault represented by a screw dislocation model, with optional
% shallow creep, additive white noise, spatially correlated noise, and
% flexible NaN masking schemes.
%
% The function is designed for testing filtering, strain-rate estimation,
% and inversion sensitivity under realistic noise and data-gap conditions.
%
% ----------------------------
% INPUTS
% ----------------------------
% slip_rate               - Deep fault slip rate (mm/yr), e.g. 5
% locking_depth_km        - Locking depth of the fault (km), e.g. 20
% creep_rate              - Shallow creep rate (mm/yr), e.g. 1
% creep_depth_km          - Depth extent of shallow creep (km),
%                           e.g. 5 means creep occurs from 0–5 km depth
% white_noise_std         - Standard deviation of Gaussian white noise (mm/yr)
% correlated_noise_std    - Standard deviation of spatially correlated noise
%                           (mm/yr)
% grid_extent_km          - Spatial extent of the grid [xmin xmax ymin ymax]
%                           in km, e.g. [-100 100 -50 50]
% grid_res_km             - Grid resolution in km, e.g. 1
%
% nan_type                - Type of NaN mask to apply:
%                             'none'        : no NaNs
%                             'blocks'      : rectangular NaN blocks
%                             'distributed' : randomly distributed NaNs with
%                                              different densities across
%                                              the fault (y = 0)
%                             'both'        : distributed + blocks
%
% nan_params              - nan_params.do_plot        : true/false(default: false)
%                           
%                           Structure containing NaN parameters, depending
%                           on nan_type:
%
%                           For nan_type = 'blocks' or 'both':
%                             nan_params.blocks = {
%                                 {r1,r2,c1,c2}, ...
%                                 {r1b,r2b,c1b,c2b}, ...
%                             };
%
%                           For nan_type = 'distributed' or 'both':
%
%   A) bbox 模式：若提供 nan_params.bbox_km = [xmin xmax ymin ymax]
%      则 bbox 内 NaN 比例  nan_params.p_in，bbox 外 nan_params.p_out
%   B) split 模式：否则使用 y_split_km，将 Y>y_split_km 与 Y<y_split_km 赋不同 NaN 比例
%      (p_above / p_below)
%                             nan_params.y_split_km = 20;    % 分界线 y = 20 km
%                             nan_params.p_above    = 0.40;  % y > 20 的 NaN 比例
%                             nan_params.p_below    = 0.10;  % y < 20 的 NaN 比例
%                             nan_params.include_split = false; % y==20 是否归到 above

%                             nan_params.seed           : RNG seed                          
%                             nan_params.respect_existing : avoid overwriting
%                                                           existing NaNs
%
% ----------------------------
% OUTPUTS
% ----------------------------
% V_clean                 - Clean synthetic velocity field (mm/yr)
% V_noisy                 - Noisy velocity field with NaNs applied (mm/yr)
% X, Y                    - Coordinate grids in km
% nan_mask                - Logical mask indicating NaN locations
%
% ----------------------------
% EXAMPLES
% ----------------------------
% Example 1: Distributed NaNs with asymmetric density across the fault
%
% nan_type = 'distributed';
% nan_params.y_split_km = 0;    % 分界线 y = 0 km
% nan_params.p_above    = 0.40;  % y > 0 的 NaN 比例
% nan_params.p_below    = 0.10;  % y < 0 的 NaN 比例
% nan_params.seed  = 42;
%
% [V_clean, V_noisy, X, Y] = generate_synthetic_insar_velocities( ...
%     5, 20, 1, 5, 1, 1, [-100 100 -50 50], 1, nan_type, nan_params);
%
%
% Example 2: Distributed NaNs + rectangular NaN gaps
%
% nan_type = 'both';
% nan_params.y_split_km = 20;    % 分界线 y = 20 km
% nan_params.p_above    = 0.40;  % y > 20 的 NaN 比例
% nan_params.p_below    = 0.10;  % y < 20 的 NaN 比例
% nan_params.blocks = { {10,30,20,40}, {60,80,100,140} };
% 
% [V_clean, V_noisy] = generate_synthetic_insar_velocities( ...
%     5, 20, 1, 5, 1, 1, [-100 100 -50 50], 1, nan_type, nan_params);
%
% 
% Example 3: Distributed NaNs + rectangular NaN gaps
%
% fault_type = "step";
% fault_params.x0_km = 0;
% fault_params.y_left_km  = -50;
% fault_params.y_right_km =  50;
% 
% nan_type = "distributed";
% nan_params.bbox_km = [-50 50 -20 20];
% nan_params.p_in = 0.4;
% nan_params.p_out = 0.1;
% 
% [V_clean, V_noisy, X, Y] = generate_synthetic_insar_velocities( ...
%     5, 20, 0, 5, 1, 0.2, [-200 200 -200 200], 1, fault_type, fault_params, nan_type, nan_params);
% ----------------------------
% NOTES
% ----------------------------
% - The fault is assumed to lie at y = 0, with positive/negative y defining
%   the two sides of the fault.
% - All random processes are reproducible via explicit RNG seeds.
% - This function is suitable for synthetic tests of filtering, strain-rate
%   recovery, and slip-deficit inversion.

    %% ---------- 0. Defaults ----------
    
    % ---- fault geometry defaults ----
    if nargin < 9 || isempty(fault_type)
        fault_type = 'straight';   % 默认：直线断层
    end
    
    if nargin < 10 || isempty(fault_params)
        fault_params = struct();
    end
    
    % ---- NaN defaults ----
    if nargin < 11 || isempty(nan_type)
        nan_type = 'none';         % 默认：不加 NaN
    end
    
    if nargin < 12 || isempty(nan_params)
        nan_params = struct();
    end

    %% ---------- 1. Grid ----------
    x_vec = grid_extent_km(1):grid_res_km:grid_extent_km(2);
    y_vec = grid_extent_km(3):grid_res_km:grid_extent_km(4);
    [X, Y] = meshgrid(x_vec, y_vec);
    Y_m = Y * 1000;

    %% ---------- 2. Dislocation + creep ----------
    d1 = locking_depth_km * 1000;
    d2 = creep_depth_km   * 1000;

    s = slip_rate;
    C = creep_rate;

    % ---------- fault geometry: y_fault(X) in km ----------
    fault_type = lower(string(fault_type));

    switch fault_type
        case "straight"
            % default: fault at y = y0 (km), default 0
            y0 = getfieldwithdefault(fault_params,'y0_km',0);
            y_fault_km = y0 * ones(size(X));

        case "step"
            % step fault:
            % X < x0  -> y = y_left
            % X >= x0 -> y = y_right
            x0     = getfieldwithdefault(fault_params,'x0_km',0);
            y_left = getfieldwithdefault(fault_params,'y_left_km',-50);
            y_right= getfieldwithdefault(fault_params,'y_right_km', 50);

            y_fault_km = y_left * ones(size(X));
            y_fault_km(X >= x0) = y_right;

        otherwise
            error("Unknown fault_type: %s. Use 'straight' or 'step'.", fault_type);
    end

    % signed perpendicular distance to the fault (meters)
    x_perp = (Y - y_fault_km) * 1000;

    % Heaviside for creep term
    H = double(x_perp >= 0);

    V_clean =  (s/pi) * atan2(x_perp, d1) ...
              - C * ( (1/pi) * atan2(x_perp, d2) - H );

    %% ---------- 3. Noise ----------
    rng(42);
    white_noise = white_noise_std * randn(size(V_clean));

    rng(38);
    noise_corr_1 = generate_correlated_noise_fft_Lc(X, Y, 30,  correlated_noise_std);
    rng(35);
    noise_corr_2 = generate_correlated_noise_fft_Lc(X, Y, 100, correlated_noise_std);

    V_noisy = V_clean + white_noise + noise_corr_1 + noise_corr_2;

    %% ---------- 4. NaN mask ----------
    nan_mask = false(size(V_clean));
    nan_type = lower(string(nan_type));

    % -------- blocks --------
    function apply_blocks()
        if ~isfield(nan_params,'blocks') || isempty(nan_params.blocks)
            return;
        end
        blocks = nan_params.blocks;
        if numel(blocks) == 4 && ~iscell(blocks{1})
            blocks = {blocks};
        end
        for k = 1:numel(blocks)
            b = blocks{k};
            nan_mask(b{1}:b{2}, b{3}:b{4}) = true;
        end
    end

% -------- distributed (UPDATED) --------
% 支持两种 distributed 模式（自动判别）：
%   A) bbox 模式：若提供 nan_params.bbox_km = [xmin xmax ymin ymax]
%      则 bbox 内 NaN 比例 = p_in，bbox 外 = p_out
%   B) split 模式：否则使用 y_split_km，将 Y>y0 与 Y<y0 赋不同 NaN 比例
%      (p_above / p_below)

function apply_distributed()

    % ---------- common defaults ----------
    seed = getfieldwithdefault(nan_params,'seed',2026);
    respect_existing = getfieldwithdefault(nan_params,'respect_existing',true);
    rng(seed);

    % ---------- Mode A: bbox in/out ----------
    if isfield(nan_params,'bbox_km') && ~isempty(nan_params.bbox_km)

        bbox = nan_params.bbox_km;  % [xmin xmax ymin ymax] in km
        if numel(bbox) ~= 4
            error('nan_params.bbox_km must be a 1x4 vector: [xmin xmax ymin ymax] (km).');
        end
        xmin = bbox(1); xmax = bbox(2);
        ymin = bbox(3); ymax = bbox(4);

        p_in  = getfieldwithdefault(nan_params,'p_in', 0.30);
        p_out = getfieldwithdefault(nan_params,'p_out',0.05);

        if p_in < 0 || p_in > 1 || p_out < 0 || p_out > 1
            error('distributed NaN fractions (p_in/p_out) must be within [0,1].');
        end

        mask_in  = (X >= xmin) & (X <= xmax) & (Y >= ymin) & (Y <= ymax);
        mask_out = ~mask_in;

        if respect_existing
            valid = ~nan_mask & ~isnan(V_noisy);
            mask_in  = mask_in  & valid;
            mask_out = mask_out & valid;
        end

        u = rand(size(V_clean));
        nan_mask = nan_mask | (mask_in  & (u < p_in)) | (mask_out & (u < p_out));

    % ---------- Mode B: y-split above/below ----------
    else
        y0 = getfieldwithdefault(nan_params,'y_split_km',0);      % km
        pA = getfieldwithdefault(nan_params,'p_above',0.40);      % Y > y0
        pB = getfieldwithdefault(nan_params,'p_below',0.10);      % Y < y0
        include_split = getfieldwithdefault(nan_params,'include_split',false);

        if pA < 0 || pA > 1 || pB < 0 || pB > 1
            error('distributed NaN fractions (p_above/p_below) must be within [0,1].');
        end

        if include_split
            mask_above = (Y >= y0);
            mask_below = (Y <  y0);
        else
            mask_above = (Y >  y0);
            mask_below = (Y <  y0);
        end

        if respect_existing
            valid = ~nan_mask & ~isnan(V_noisy);
            mask_above = mask_above & valid;
            mask_below = mask_below & valid;
        end

        u = rand(size(V_clean));
        nan_mask = nan_mask | (mask_above & (u < pA)) | (mask_below & (u < pB));
    end
end

% -------- apply by type --------
    switch nan_type
        case "none"
            % nothing
        case "blocks"
            apply_blocks();
        case "distributed"
            apply_distributed();
        case "both"
            apply_distributed();
            apply_blocks();
        otherwise
            error("Unknown nan_type: %s", nan_type);
    end

    %% ---------- 5. Apply ----------
    V_clean(nan_mask) = NaN;
    V_noisy(nan_mask) = NaN;

    %% ---------- 6. Optional plotting ----------
    % Controlled by nan_params.do_plot (default: false)

    do_plot = false;
    if isfield(nan_params,'do_plot')
        do_plot = nan_params.do_plot;
    end

    if do_plot

        % ---- load colormap (vik) if available ----
        if exist('vik.mat','file')
            vik = importdata('vik.mat');
            cmap = vik;
        else
            warning('vik.mat not found, using parula instead.');
            cmap = parula;
        end

        figure('Color','w','Position',[100 100 1200 700]);

        % ---------- Clean velocity ----------
        subplot(2,3,1);
        h1 = imagesc(x_vec, y_vec, V_clean, [-5 5]);
        axis image xy;
        colormap(gca, cmap);
        set(gca,'Color',[1 1 1]);
        set(h1,'AlphaData',~isnan(V_clean));
        colorbar;
        title('Clean velocity (lock + creep)','FontWeight','normal');

        % ---------- Noisy velocity ----------
        subplot(2,3,2);
        h2 = imagesc(x_vec, y_vec, V_noisy, [-5 5]);
        axis image xy;
        colormap(gca, cmap);
        set(gca,'Color',[1 1 1]);
        set(h2,'AlphaData',~isnan(V_noisy));
        colorbar;
        title('Noisy velocity + NaNs','FontWeight','normal');

        % ---------- White noise ----------
        subplot(2,3,3);
        h3 = imagesc(x_vec, y_vec, white_noise, [-3 3]);
        axis image xy;
        colormap(gca, cmap);
        set(gca,'Color',[1 1 1]);
        set(h3,'AlphaData',~isnan(white_noise));
        colorbar;
        title(sprintf('White noise (%.2f mm/yr)',white_noise_std), ...
              'FontWeight','normal');

        % ---------- Correlated noise (30 km) ----------
        subplot(2,3,4);
        h4 = imagesc(x_vec, y_vec, noise_corr_1, [-1 1]);
        axis image xy;
        colormap(gca, cmap);
        set(gca,'Color',[1 1 1]);
        set(h4,'AlphaData',~isnan(noise_corr_1));
        colorbar;
        title(sprintf('Correlated noise (L_c = 30 km)'), ...
              'FontWeight','normal');

        % ---------- Correlated noise (100 km) ----------
        subplot(2,3,5);
        h5 = imagesc(x_vec, y_vec, noise_corr_2, [-1 1]);
        axis image xy;
        colormap(gca, cmap);
        set(gca,'Color',[1 1 1]);
        set(h5,'AlphaData',~isnan(noise_corr_2));
        colorbar;
        title(sprintf('Correlated noise (L_c = 100 km)'), ...
              'FontWeight','normal');

        % ---------- Total noise ----------
        subplot(2,3,6);
        h6 = imagesc(x_vec, y_vec, white_noise + noise_corr_1 + noise_corr_2, [-3 3]);
        axis image xy;
        colormap(gca, cmap);
        set(gca,'Color',[1 1 1]);
        set(h6,'AlphaData',~isnan(V_noisy));
        colorbar;
        title('Total noise','FontWeight','normal');

        sgtitle('Synthetic InSAR velocity & noise components','FontWeight','bold');
    end

end

%% ---------- utility ----------
function v = getfieldwithdefault(s, f, d)
    if isfield(s,f)
        v = s.(f);
    else
        v = d;
    end
end

