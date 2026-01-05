function [V_clean, V_noisy, X, Y] = generate_synthetic_insar_velocities( ...
    slip_rate, locking_depth_km, creep_rate, creep_depth_km, noise_std, ...
    grid_extent_km, grid_res_km, nan_blocks)
% GENERATE_SYNTHETIC_INSAR_VELOCITIES
% Generates synthetic fault-parallel velocities for a screw dislocation
% representing a strike-slip fault in 2D, with added white/correlated Gaussian
% noise and NaNs.
%
% Inputs:
%   slip_rate         - fault slip rate in mm/yr (e.g., 5)
%   locking_depth_km  - locking depth in km (e.g., 20)
%   creep_rate        - shallow creep rate (e.g., 1)
%   creep_depth_km    - shallow creep depth in km (e.g., 5 indicates that shallow creep occurs within 0 to 5 km near the surface)
%   noise_std         - std deviation of Gaussian noise in mm/yr (e.g., 1)
%   grid_extent_km    - [xmin xmax ymin ymax] in km (e.g., [-100 100 -50 50])
%   grid_res_km       - grid resolution in km (e.g., 1)
%   nan_blocks        - 可以是：
%                         {r1, r2, c1, c2}
%                       或
%                         {{r1, r2, c1, c2}, {r1b, r2b, c1b, c2b}, ...}
%
% Outputs:
%   V_clean           - clean velocity field (mm/yr)
%   V_noisy           - noisy velocity field with NaNs (mm/yr)
%   X, Y              - coordinate grids in km

    % ---------- 1. Grid setup ----------
    x_vec = grid_extent_km(1):grid_res_km:grid_extent_km(2);   % 1D for columns
    y_vec = grid_extent_km(3):grid_res_km:grid_extent_km(4);   % 1D for rows
    [X, Y] = meshgrid(x_vec, y_vec);                           % 2D grids (km)
    Y_m = Y * 1000;                                            % convert to meters

    % ---------- 2. Dislocation + shallow creep model ----------

    locking_depth_m = locking_depth_km * 1000;   % d1
    d1 = locking_depth_m;

    % 浅层蠕滑到 5 km 深度
    d2 = creep_depth_km * 1000;                  % d2, 例子：0–5 km 蠕滑

    % 深部锁固滑动速率 s（mm/yr）
    s = slip_rate;

    % 浅层蠕滑速率 C（例如固定为 2 mm/yr；也可以改成比例）
    C = creep_rate;                                       % 蠕滑速率（mm/yr）

    % 断层位置偏移（公里），如果剖面 0 刚好在断层上，可以设为 0
    xc_km = 0;
    xc_m  = xc_km * 1000;

    % fault-perpendicular 坐标（矩阵），这里用 Y 方向作为剖面
    x_perp = Y_m + xc_m;                         % 与公式中的 x + xc 对应

    % Heaviside 函数 H(x + xc)，x+xc>=0 为 1，其余为 0
    H = double(x_perp >= 0);

    % 组合：深部锁固 + 浅层蠕滑
    % v(x) = -(s/pi)*atan((x+xc)/d1) + C[ (1/pi)*atan((x+xc)/d2) - H(x+xc) ]
    V_clean =  (s/pi) * atan2(x_perp, d1) ...           % 深部锁固对应的螺位错
              - C * ( (1/pi) * atan2(x_perp, d2) ...     % 全深度蠕滑
                     - H );                              % 减去地表以上的滑动
    % 单位依然是 mm/yr

    % ---------- 3. Correlated Gaussian noise (FFT) ----------
    rng(42);  % reproducibility
    white_noise = noise_std * randn(size(V_clean));
    V_noisy = V_clean + white_noise;
    % rng(38);  % reproducibility
    % corr_Lc_km = 50;   % 手动控制大气平滑程度
    % noise_corr = generate_correlated_noise_fft_Lc(X, Y, corr_Lc_km, noise_std);
    % V_noisy = V_clean + white_noise + noise_corr;

    % ---------- 4. Apply NaN masks ----------
    if nargin < 6 || isempty(nan_blocks)
        nan_blocks = {};
    end

    % 统一转换成 cell-of-cells 形式
    if ~isempty(nan_blocks)
        if ~iscell(nan_blocks)
            error('nan_blocks 必须是 cell。');
        end

        % 情况 1：直接传 {r1, r2, c1, c2}
        if numel(nan_blocks) == 4 && ~iscell(nan_blocks{1})
            nan_blocks = {nan_blocks};
        % 情况 2：传 {{...}, {...}, ...}，直接用
        elseif numel(nan_blocks) >= 1 && iscell(nan_blocks{1})
            % do nothing
        else
            error('nan_blocks 格式不符合预期。请用 {r1,r2,c1,c2} 或 {{r1,r2,c1,c2},...}。');
        end

        for k = 1:numel(nan_blocks)
            blk = nan_blocks{k};   % 这里 blk 一定是 1x4 cell
            r1 = blk{1}; r2 = blk{2};
            c1 = blk{3}; c2 = blk{4};
            V_clean(r1:r2, c1:c2) = NaN;
            V_noisy(r1:r2, c1:c2) = NaN;
        end
    end

    % ---------- 5. Load colormap (optional) ----------
    vik = importdata('vik.mat');
    cpt.vik = vik; %#ok<NASGU>

    % % ---------- Optional plotting ----------
    % figure;
    % 
    % subplot(2,2,1);
    % h1 = imagesc(x_vec, y_vec, V_clean, [-5 5]); axis image xy;
    % colormap(gca, cpt.vik);
    % set(gca, 'Color', [1 1 1]);
    % set(h1, 'AlphaData', ~isnan(V_clean));
    % colorbar;
    % title('Synth. Fault-Parallel Velocity (lock+creep)');
    % 
    % subplot(2,2,2);
    % h2 = imagesc(x_vec, y_vec, V_noisy, [-5 5]); axis image xy;
    % colormap(gca, cpt.vik);
    % set(gca, 'Color', [1 1 1]);
    % set(h2, 'AlphaData', ~isnan(V_noisy));
    % colorbar;
    % title('Noisy Velocity with White/Correlated Noise & NaNs');
    % 
    % subplot(2,2,3);
    % h3 = imagesc(x_vec, y_vec, white_noise, [-1 1]); axis image xy;
    % colormap(gca, cpt.vik);
    % set(gca, 'Color', [1 1 1]);
    % set(h3, 'AlphaData', ~isnan(white_noise));
    % colorbar;
    % title(sprintf('White Noise (%.2f mm/yr)', noise_std));
    % 
    % subplot(2,2,4);
    % h4 = imagesc(x_vec, y_vec, noise_corr, [-1 1]); axis image xy;
    % colormap(gca, cpt.vik);
    % set(gca, 'Color', [1 1 1]);
    % set(h4, 'AlphaData', ~isnan(noise_corr));
    % colorbar;
    % title(sprintf('Correlated Noise (%.2f mm/yr, 50 km L_c)', noise_std));

end
