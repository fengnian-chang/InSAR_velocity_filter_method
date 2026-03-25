clear all 
% 走滑断层模型参数
slip_rate = 5;                     % 自由滑动速率 mm/yr
locking_depth_km = 20;             % 锁固深度 km
creep_rate = 2;                    % 浅层蠕滑速率 mm/yr
creep_depth_km = 5;                % 浅层蠕滑深度 km
white_noise_std = 1;               % 白噪声（Gaussian）强度 mm/yr
correlated_noise_std = 0.2;        % 空间相关噪声（Gaussian）强度 mm/yr（correlated length 30km + 100 km）
grid_extent_km = [-200 200 -200 200]; % 区域范围 X: [-200,200] km, Y: [-200,200] km
grid_res_km = 1;                  % 网格分辨率 km

% 断层类型
fault_type = "step";
% fault_type = "straight";
fault_params.x0_km = 0;
fault_params.y_left_km  = -50;
fault_params.y_right_km =  50;

% 设置两个 NaN 区域（用于测试 NaN-aware 滤波器）
% nan_type = 'both';
% nan_type = 'distributed';
% nan_type = 'blocks';
nan_type = 'none';

nan_params.do_plot = false; 
nan_params.bbox_km = [-100 0 25 200];
nan_params.p_in = 0.30;
nan_params.p_out = 0;
% nan_params.y_split_km = 50;
% nan_params.p_above = 0.60;  
% nan_params.p_below = 0;

nan_params.blocks = {
    {220, 260, 260, 340};    % 第一块遮挡区域
};

% 调用函数生成数据
[V_clean, V_noisy, X, Y] = generate_synthetic_insar_velocities( ...
    slip_rate, locking_depth_km, creep_rate, creep_depth_km, white_noise_std, correlated_noise_std, grid_extent_km, grid_res_km, fault_type, fault_params, nan_type, nan_params);

%%
% 滤波参数
spatial_sigma = 10;
range_sigma = 2;
spatial_winsize = 60;

center_winsize = 25;
% center_sigma = 5;   % 大概2σ对应half_winsize(95.45%),3σ 99.8%, 当center_mode = 'gaussian'时使用

gauss_kernel = fspecial('gaussian', [spatial_winsize spatial_winsize], spatial_sigma);

% 高斯滤波（nan-aware）
V_gaussian = nanconv(V_noisy, gauss_kernel, 'edge', 'nanout');

% 中值滤波（nan-aware）
V_median = nanmedfilt2(V_noisy, [spatial_winsize spatial_winsize], 'nanout');

% LOESS滤波（nan-aware）
% V_loess = loess_nanfit2_plane(V_noisy, spatial_winsize, 'nanout');

% 双边滤波（NaN-aware）
V_bilateral = bilateral_nanconv(V_noisy, spatial_sigma, range_sigma, spatial_winsize, center_winsize, 'median', 'nanout');
%V_bilateral = bilateral_nanconv(V_noisy, spatial_sigma, range_sigma, spatial_winsize, center_winsize, 'gaussian', center_sigma, 'nanout');

% 计算 dv/dy（向北方向的速度梯度）
dy = grid_res_km;
[~, dVc_dy] = gradient(V_clean, dy);   % 合成速率的应变率
[~, dVg_dy] = gradient(V_gaussian, dy);   % 高斯滤波后的应变率
[~, dVm_dy] = gradient(V_median, dy);     % 中值滤波后的应变率
% [~, dVl_dy] = gradient(V_loess, 1);    % LOESS滤波后的应变率
[~, dVb_dy] = gradient(V_bilateral, dy);  % 双边滤波后的应变率

% 换算为 nstrain/yr
dVc_dy = dVc_dy * 1e3; 
dVg_dy = dVg_dy * 1e3; 
dVm_dy = dVm_dy * 1e3; 
% dVl_dy = dVl_dy * 1e3;
dVb_dy = dVb_dy * 1e3;

% 计算 dv/dx（向东方向的速度梯度）
dx = grid_res_km;
[dVc_dx, ~] = gradient(V_clean, dx);   % 合成速率的应变率
[dVg_dx, ~] = gradient(V_gaussian, dx);   % 高斯滤波后的应变率
[dVm_dx, ~] = gradient(V_median, dx);     % 中值滤波后的应变率
[dVb_dx, ~] = gradient(V_bilateral, dx);  % 双边滤波后的应变率

% 换算为 nstrain/yr
dVc_dx = dVc_dx * 1e3; 
dVg_dx = dVg_dx * 1e3; 
dVm_dx = dVm_dx * 1e3; 
dVb_dx = dVb_dx * 1e3;

% Load colormap
vik = importdata('vik.mat');
cpt.vik = vik;

% %% 提取剖面 (x = -10~10 km; y 全部)
% xmask1 = X(1,:) >= -110 & X(1,:) <= -90;    % green profile
% xmask2 = X(1,:) >= 90 & X(1,:) <= 110;  % red profile
% xmask_clean = X(1,:) >= -10 & X(1,:) <= 10; %用于绘制完整无间断剖面（合成数据基准）
% 
% prof_y = Y(:,1);  % y 方向
% prof_Vc     = mean(V_clean(:, xmask_clean), 2, "omitnan"); %完整clean剖面
% % profile 1
% prof_Vn1     = mean(V_noisy(:, xmask1), 2, "omitnan");
% prof_Vm1     = mean(V_median(:, xmask1), 2, "omitnan");
% prof_Vg1     = mean(V_gaussian(:, xmask1), 2, "omitnan");
% % prof_Vl1     = mean(V_loess(:, xmask1), 2, "omitnan");
% prof_Vb1     = mean(V_bilateral(:, xmask1), 2, "omitnan");
% 
% prof_dVc1    = mean(dVc_dy(:, xmask1), 2, "omitnan");
% prof_dVm1    = mean(dVm_dy(:, xmask1), 2, "omitnan");
% prof_dVg1    = mean(dVg_dy(:, xmask1), 2, "omitnan");
% % prof_dVl1    = mean(dVl_dy(:, xmask1), 2, "omitnan");
% prof_dVb1    = mean(dVb_dy(:, xmask1), 2, "omitnan");
% 
% % profile 2
% prof_Vn2     = mean(V_noisy(:, xmask2), 2, "omitnan");
% prof_Vm2     = mean(V_median(:, xmask2), 2, "omitnan");
% prof_Vg2     = mean(V_gaussian(:, xmask2), 2, "omitnan");
% % prof_Vl2     = mean(V_loess(:, xmask2), 2, "omitnan");
% prof_Vb2     = mean(V_bilateral(:, xmask2), 2, "omitnan");
% 
% prof_dVc2    = mean(dVc_dy(:, xmask2), 2, "omitnan");
% prof_dVm2    = mean(dVm_dy(:, xmask2), 2, "omitnan");
% prof_dVg2    = mean(dVg_dy(:, xmask2), 2, "omitnan");
% % prof_dVl2    = mean(dVl_dy(:, xmask2), 2, "omitnan");
% prof_dVb2    = mean(dVb_dy(:, xmask2), 2, "omitnan");

%% 可视化开始
figure;

set(groot, 'defaultAxesFontName', 'Calibri');
set(groot, 'defaultTextFontName', 'Calibri');

%%%% ========================
%%%% 1. 第一行（Noisy）
%%%% ========================

%% ---- (1,1) Noisy velocity ----
subplot(4,3,1)
hV1 = imagesc(X(1,:),Y(:,1),V_noisy,[-4 4]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hV1,'AlphaData',~isnan(V_noisy));
title("Synth Ve with noise");
colormap(vik); colorbar;
hold on;
% rectangle('Position',[-110,-200,20,400],'EdgeColor','b','LineWidth',1.5);
% rectangle('Position',[90,-200,20,400],'EdgeColor','r','LineWidth',1.5);
% rectangle('Position',[-150,50,101,101],'EdgeColor','g','LineWidth',1.5);

%% ---- (1,2) Clean strain ----
subplot(4,3,2)
hS1 = imagesc(X(1,:),Y(:,1),dVc_dy,[-100 100]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hS1,'AlphaData',~isnan(dVc_dy));
title("∂Ve/∂y Synth (without noise)"); colormap(vik); colorbar;
hold on;

subplot(4,3,3)
hS1 = imagesc(X(1,:),Y(:,1),dVc_dx,[-100 100]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hS1,'AlphaData',~isnan(dVc_dx));
title("∂Ve/∂x Synth (without noise)"); colormap(vik); colorbar;
hold on;

%%%% ========================
%%%% 2. Median
%%%% ========================

subplot(4,3,4)
hV2 = imagesc(X(1,:),Y(:,1),V_median,[-4 4]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hV2,'AlphaData',~isnan(V_median));
title(sprintf("Median (%d×%d)",spatial_winsize,spatial_winsize));

subplot(4,3,5)
hS2 = imagesc(X(1,:),Y(:,1),dVm_dy,[-100 100]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hS2,'AlphaData',~isnan(dVm_dy));
title("∂Ve/∂y Median");

subplot(4,3,6)
hS2 = imagesc(X(1,:),Y(:,1),dVm_dx,[-100 100]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hS2,'AlphaData',~isnan(dVm_dx));
title("∂Ve/∂x Median");

%%%% ========================
%%%% 3. Gaussian
%%%% ========================

subplot(4,3,7)
hV3 = imagesc(X(1,:),Y(:,1),V_gaussian,[-4 4]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hV3,'AlphaData',~isnan(V_gaussian));
title(sprintf("Gaussian (%d×%d,σ_s=%d)",spatial_winsize,spatial_winsize,spatial_sigma));

subplot(4,3,8)
hS3 = imagesc(X(1,:),Y(:,1),dVg_dy,[-100 100]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hS3,'AlphaData',~isnan(dVg_dy));
title("∂Ve/∂y Gaussian");

subplot(4,3,9)
hS3 = imagesc(X(1,:),Y(:,1),dVg_dx,[-100 100]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hS3,'AlphaData',~isnan(dVg_dx));
title("∂Ve/∂x Gaussian");

%%%% ========================
%%%% 4. Bilateral
%%%% ========================

subplot(4,3,10)
hV5 = imagesc(X(1,:),Y(:,1),V_bilateral,[-4 4]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hV5,'AlphaData',~isnan(V_bilateral));
title(sprintf("Bilateral (%d×%d,σ_s=%d,σ_v=%d)",...
       spatial_winsize,spatial_winsize,spatial_sigma,range_sigma));

subplot(4,3,11)
hS5 = imagesc(X(1,:),Y(:,1),dVb_dy,[-100 100]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hS5,'AlphaData',~isnan(dVb_dy));
title("∂Ve/∂y Bilateral");

subplot(4,3,12)
hS5 = imagesc(X(1,:),Y(:,1),dVb_dx,[-100 100]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hS5,'AlphaData',~isnan(dVb_dx));
title("∂Ve/∂x Bilateral");
