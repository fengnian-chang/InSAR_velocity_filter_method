clear all 
% 走滑断层模型参数
slip_rate = 5;                     % 自由滑动速率 mm/yr
locking_depth_km = 20;             % 锁固深度 km
creep_rate = 0;                    % 浅层蠕滑速率 mm/yr
creep_depth_km = 5;                % 浅层蠕滑深度 km
white_noise_std = 1;               % 白噪声（Gaussian）强度 mm/yr
correlated_noise_std = 0.2;        % 空间相关噪声（Gaussian）强度 mm/yr（correlated length 30km + 100 km）
grid_extent_km = [-200 200 -200 200]; % 区域范围 X: [-200,200] km, Y: [-200,200] km
grid_res_km = 1;                  % 网格分辨率 km

% 断层类型
% fault_type = "step";      % used for test fault step
fault_params.x0_km = 0;     % Invalid when fault_type = "straight"
fault_params.y_left_km  = -50;
fault_params.y_right_km =  50;

fault_type = "straight";

% 设置两个 NaN 区域（用于测试 NaN-aware 滤波器）
% nan_type = 'both';
% nan_type = 'distributed';
nan_type = 'blocks';
% nan_type = 'none';

nan_params.do_plot = true;      % plot velocity and NaN settings

% nan_params.bbox_km = [0 200 15 200];    % pixel density inside and outside the box
% nan_params.p_in = 0.40;
% nan_params.p_out = 0;

% nan_params.y_split_km = 50;       % pixel density above and below the y_split
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
spatial_sigma = 15;
range_sigma = 2;
spatial_winsize = 100;

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

% Load colormap
vik = importdata('vik.mat');
cpt.vik = vik;

%% 提取剖面 (x = -10~10 km; y 全部)
xmask1 = X(1,:) >= -110 & X(1,:) <= -90;    % green profile
xmask2 = X(1,:) >= 90 & X(1,:) <= 110;  % red profile
xmask_clean = X(1,:) >= -10 & X(1,:) <= 10; %用于绘制完整无间断剖面（合成数据基准）

prof_y = Y(:,1);  % y 方向
prof_Vc     = mean(V_clean(:, xmask_clean), 2, "omitnan"); %完整clean剖面
% profile 1
prof_Vn1     = mean(V_noisy(:, xmask1), 2, "omitnan");
prof_Vm1     = mean(V_median(:, xmask1), 2, "omitnan");
prof_Vg1     = mean(V_gaussian(:, xmask1), 2, "omitnan");
% prof_Vl1     = mean(V_loess(:, xmask1), 2, "omitnan");
prof_Vb1     = mean(V_bilateral(:, xmask1), 2, "omitnan");

prof_dVc1    = mean(dVc_dy(:, xmask1), 2, "omitnan");
prof_dVm1    = mean(dVm_dy(:, xmask1), 2, "omitnan");
prof_dVg1    = mean(dVg_dy(:, xmask1), 2, "omitnan");
% prof_dVl1    = mean(dVl_dy(:, xmask1), 2, "omitnan");
prof_dVb1    = mean(dVb_dy(:, xmask1), 2, "omitnan");

% profile 2
prof_Vn2     = mean(V_noisy(:, xmask2), 2, "omitnan");
prof_Vm2     = mean(V_median(:, xmask2), 2, "omitnan");
prof_Vg2     = mean(V_gaussian(:, xmask2), 2, "omitnan");
% prof_Vl2     = mean(V_loess(:, xmask2), 2, "omitnan");
prof_Vb2     = mean(V_bilateral(:, xmask2), 2, "omitnan");

prof_dVc2    = mean(dVc_dy(:, xmask2), 2, "omitnan");
prof_dVm2    = mean(dVm_dy(:, xmask2), 2, "omitnan");
prof_dVg2    = mean(dVg_dy(:, xmask2), 2, "omitnan");
% prof_dVl2    = mean(dVl_dy(:, xmask2), 2, "omitnan");
prof_dVb2    = mean(dVb_dy(:, xmask2), 2, "omitnan");

%% 可视化开始
figure;

set(groot, 'defaultAxesFontName', 'Calibri');
set(groot, 'defaultTextFontName', 'Calibri');

%%%% ========================
%%%% 1. 第一行（Noisy）
%%%% ========================

%% ---- (1,1) Noisy velocity ----
subplot(4,6,1)
hV1 = imagesc(X(1,:),Y(:,1),V_noisy,[-4 4]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hV1,'AlphaData',~isnan(V_noisy));
title("Synth Ve with noise");
colormap(vik); colorbar;
hold on;
rectangle('Position',[-110,-200,20,400],'EdgeColor','b','LineWidth',1.5);
rectangle('Position',[90,-200,20,400],'EdgeColor','r','LineWidth',1.5);
rectangle('Position',[-150,50,101,101],'EdgeColor','g','LineWidth',1.5);

%% ---- (1,2)(1,3) Profile ----
% 准备采样点（灰色散点）
Yvec = prof_y(:);  % Ny×1

% profile 1 window: X in [-110,-90]
Vpts1 = V_noisy(:, xmask1);                 % Ny×Nx1
ypts1 = repmat(Yvec, 1, size(Vpts1,2));     % Ny×Nx1

% profile 2 window: X in [90,110]
Vpts2 = V_noisy(:, xmask2);                 % Ny×Nx2
ypts2 = repmat(Yvec, 1, size(Vpts2,2));     % Ny×Nx2

% ---- (1,2) ----
subplot(4,6,2)

% 灰色采样点（去掉NaN）
m1 = ~isnan(Vpts1);
plot(ypts1(m1), Vpts1(m1), '.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 6); hold on;

plot(prof_y, prof_Vc,'k:','LineWidth',1); hold on;

% 均值曲线
plot(prof_y, prof_Vn1,'b','LineWidth',0.6); grid on;
title("Profile_3 Ve (noisy)");

% ---- (1,3) ----
subplot(4,6,3)

m2 = ~isnan(Vpts2);
plot(ypts2(m2), Vpts2(m2), '.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 6); hold on;

plot(prof_y, prof_Vc,'k:','LineWidth',1); hold on;

plot(prof_y, prof_Vn2,'r','LineWidth',0.6); grid on;
title("Profile_4 Ve (noisy)");


%% ---- (1,4) Clean strain ----
subplot(4,6,4)
hS1 = imagesc(X(1,:),Y(:,1),dVc_dy,[-100 100]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hS1,'AlphaData',~isnan(dVc_dy));
title("∂Ve/∂y Synth (without noise)"); colormap(vik); colorbar;
hold on;
rectangle('Position',[-110,-200,20,400],'EdgeColor','b','LineWidth',1.5);
rectangle('Position',[90,-200,20,400],'EdgeColor','r','LineWidth',1.5);

%% ---- (1,5)(1,6) Profile ----
subplot(4,6,5)
plot(prof_y, prof_dVc1,'b'); grid on;
title("Profile_3 ∂Ve/∂y Synth");
ylim([0 100]);
xlim([-100 100]);

subplot(4,6,6)
plot(prof_y, prof_dVc2,'r'); grid on;
title("Profile_4 ∂Ve/∂y Synth");
ylim([0 100]);
xlim([-100 100]);
%%%% ========================
%%%% 2. Median
%%%% ========================

subplot(4,6,7)
hV2 = imagesc(X(1,:),Y(:,1),V_median,[-4 4]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hV2,'AlphaData',~isnan(V_median));
title(sprintf("Median (%d×%d)",spatial_winsize,spatial_winsize));

subplot(4,6,8)
plot(prof_y, prof_Vc,'k:','LineWidth',1); hold on;
plot(prof_y, prof_Vm1,'b'); grid on;
title("Profile_3 Ve Median");

subplot(4,6,9)
plot(prof_y, prof_Vc,'k:','LineWidth',1); hold on;
plot(prof_y, prof_Vm2,'r'); grid on;
title("Profile_4 Ve Median");

subplot(4,6,10)
hS2 = imagesc(X(1,:),Y(:,1),dVm_dy,[-100 100]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hS2,'AlphaData',~isnan(dVm_dy));
title("∂Ve/∂y Median");

subplot(4,6,11)
plot(prof_y, prof_dVc1,'k:','LineWidth',1); hold on;
plot(prof_y, prof_dVm1,'b'); grid on;
title("Profile_3 ∂Ve/∂y Median");
ylim([0 100]);
xlim([-100 100]);

subplot(4,6,12)
plot(prof_y, prof_dVc2,'k:','LineWidth',1); hold on;
plot(prof_y, prof_dVm2,'r'); grid on;
title("Profile_4 ∂Ve/∂y Median");
ylim([0 100]);
xlim([-100 100]);
%%%% ========================
%%%% 3. Gaussian
%%%% ========================

subplot(4,6,13)
hV3 = imagesc(X(1,:),Y(:,1),V_gaussian,[-4 4]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hV3,'AlphaData',~isnan(V_gaussian));
title(sprintf("Gaussian (%d×%d,σ_s=%d)",spatial_winsize,spatial_winsize,spatial_sigma));

subplot(4,6,14)
plot(prof_y, prof_Vc,'k:','LineWidth',1); hold on;
plot(prof_y, prof_Vg1,'b'); grid on;
title("Profile_3 Ve Gaussian");

subplot(4,6,15)
plot(prof_y, prof_Vc,'k:','LineWidth',1); hold on;
plot(prof_y, prof_Vg2,'r'); grid on;
title("Profile_4 Ve Gaussian");

subplot(4,6,16)
hS3 = imagesc(X(1,:),Y(:,1),dVg_dy,[-100 100]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hS3,'AlphaData',~isnan(dVg_dy));
title("∂Ve/∂y Gaussian");

subplot(4,6,17)
plot(prof_y, prof_dVc1,'k:','LineWidth',1); hold on;
plot(prof_y, prof_dVg1,'b'); grid on;
title("Profile_3 ∂Ve/∂y Gaussian");
ylim([0 100]);
xlim([-100 100]);

subplot(4,6,18)
plot(prof_y, prof_dVc2,'k:','LineWidth',1); hold on;
plot(prof_y, prof_dVg2,'r'); grid on;
title("Profile_4 ∂Ve/∂y Gaussian");
ylim([0 100]);
xlim([-100 100]);
%%%% ========================
%%%% 4. Bilateral
%%%% ========================

subplot(4,6,19)
hV5 = imagesc(X(1,:),Y(:,1),V_bilateral,[-4 4]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hV5,'AlphaData',~isnan(V_bilateral));
title(sprintf("Bilateral (%d×%d,σ_s=%d,σ_v=%d)",...
       spatial_winsize,spatial_winsize,spatial_sigma,range_sigma));

subplot(4,6,20)
plot(prof_y, prof_Vc,'k:','LineWidth',1); hold on;
plot(prof_y, prof_Vb1,'b'); grid on;
title("Profile_3 Ve Bilateral");

subplot(4,6,21)
plot(prof_y, prof_Vc,'k:','LineWidth',1); hold on;
plot(prof_y, prof_Vb2,'r'); grid on;
title("Profile_4 Ve Bilateral");

subplot(4,6,22)
hS5 = imagesc(X(1,:),Y(:,1),dVb_dy,[-100 100]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hS5,'AlphaData',~isnan(dVb_dy));
title("∂Ve/∂y Bilateral");

subplot(4,6,23)
plot(prof_y, prof_dVc1,'k:','LineWidth',1); hold on;
plot(prof_y, prof_dVb1,'b'); grid on;
title("Profile_3 ∂Ve/∂y Bilateral");
ylim([0 100]);
xlim([-100 100]);

subplot(4,6,24)
plot(prof_y, prof_dVc2,'k:','LineWidth',1); hold on;
plot(prof_y, prof_dVb2,'r'); grid on;
title("Profile_4 ∂Ve/∂y Bilateral");
ylim([0 100]);
xlim([-100 100]);
