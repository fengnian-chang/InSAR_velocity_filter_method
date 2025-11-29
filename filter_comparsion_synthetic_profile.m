clear all 
% 走滑断层模型参数
slip_rate = 5;                     % mm/yr
locking_depth_km = 20;            % 锁固深度 km
noise_std = 0.5;                    % 噪声强度 mm/yr
grid_extent_km = [-200 200 -200 200]; % 区域范围 X: [-200,200] km, Y: [-200,200] km
grid_res_km = 1;                  % 网格分辨率 km

% 设置两个 NaN 区域（用于测试 NaN-aware 滤波器）
nan_blocks = {
    {220, 260, 160, 240};    % 第一块遮挡区域
    %{40, 80, 120, 160};    % 第二块遮挡区域
    %{40, 41, 120, 121} 
};

% 调用函数生成数据
[V_clean, V_noisy, X, Y] = generate_synthetic_insar_velocities( ...
    slip_rate, locking_depth_km, noise_std, grid_extent_km, grid_res_km, nan_blocks);

%%
% % 添加短波长水平速度扰动（Mogi源）
% mogi_center = [-80, 100];  % 中心位置 (X, Y)，单位 km
% mogi_depth = 10;            % 源深度 km
% nu = 0.25;                 % 泊松比
% deltaV = 5* pi * 10^3;           % 控制体积变化率（调节强度，默认最大速度 ~2mm/yr）
% 
% [XX, YY] = meshgrid(X(1,:), Y(:,1));
% dx = XX - mogi_center(1);
% dy = YY - mogi_center(2);
% r = sqrt(dx.^2 + dy.^2);
% R2 = r.^2 + mogi_depth^2;
% 
% % 径向水平速度 (mm/yr)
% vr = (3 * (1 - nu) * deltaV .* r) ./ (pi * R2.^2);
% 
% % 投影到东西方向
% v_mogi_E = -vr .* (dx ./ r);
% v_mogi_E(r == 0) = 0;  % 中心点处理
% % 添加局部沉降信号到原始带噪速度场
% V_noisy = V_noisy + v_mogi_E;

%%
% 添加二维高斯噪声区域
% noise_block_extent = [-120, 80, -80, 120];  % xmin, ymin, xmax, ymax
% noise_mean = 1;     % mm/yr
% noise_std = 1;      % mm/yr
% 
% % 生成逻辑掩码
% xmask = X(1,:) >= noise_block_extent(1) & X(1,:) <= noise_block_extent(3);
% ymask = Y(:,1) >= noise_block_extent(2) & Y(:,1) <= noise_block_extent(4);
% mask = ymask * xmask;  % 外积生成二维逻辑掩码
% 
% % 添加二维高斯白噪声（相同大小，仅在掩码处有效）
% rand_noise = noise_mean + noise_std * randn(size(V_noisy));
% V_noisy(mask == 1) = V_noisy(mask == 1) + rand_noise(mask == 1);

%%
% 滤波参数
spatial_sigma = 15;
range_sigma = 2;
spatial_winsize = 101;

center_winsize = 31;
center_sigma = 5;   % 大概2σ对应half_winsize(95.45%),3σ 99.8%, 当center_mode = 'gaussian'时使用

gauss_kernel = fspecial('gaussian', [spatial_winsize spatial_winsize], spatial_sigma);

% 高斯滤波（nan-aware）
V_gaussian = nanconv(V_noisy, gauss_kernel, 'edge', 'nanout');

% 中值滤波（nan-aware）
V_median = nanmedfilt2(V_noisy, [spatial_winsize spatial_winsize], 'nanout');

% LOESS滤波（nan-aware）
V_loess = loess_nanfit2_plane(V_noisy, spatial_winsize, 'nanout');

% 双边滤波（NaN-aware）
V_bilateral = bilateral_nanconv(V_noisy, spatial_sigma, range_sigma, spatial_winsize, center_winsize, 'median', 'nanout');
%V_bilateral = bilateral_nanconv(V_noisy, spatial_sigma, range_sigma, spatial_winsize, center_winsize, 'gaussian', center_sigma, 'nanout');

% 计算 dv/dy（向北方向的速度梯度）
dy = grid_res_km;
[~, dVc_dy] = gradient(V_clean, dy);   % 合成速率的应变率
[~, dVg_dy] = gradient(V_gaussian, dy);   % 高斯滤波后的应变率
[~, dVm_dy] = gradient(V_median, dy);     % 中值滤波后的应变率
[~, dVl_dy] = gradient(V_loess, 1);    % LOESS滤波后的应变率
[~, dVb_dy] = gradient(V_bilateral, dy);  % 双边滤波后的应变率

% 换算为 nstrain/yr
dVc_dy = dVc_dy * 1e3; 
dVg_dy = dVg_dy * 1e3; 
dVm_dy = dVm_dy * 1e3; 
dVl_dy = dVl_dy * 1e3;
dVb_dy = dVb_dy * 1e3;

% Load colormap
vik = importdata('vik.mat');
cpt.vik = vik;

%% 提取剖面 (x = -10~10 km; y 全部)
xmask = X(1,:) >= -10 & X(1,:) <= 10;
xmask2 = X(1,:) >= 90 & X(1,:) <= 110; %用于绘制完整无间断剖面（合成数据基准）

prof_y = Y(:,1);  % y 方向
prof_V_noisy = mean(V_noisy(:, xmask), 2, "omitnan");
prof_Vc     = mean(V_clean(:, xmask2), 2, "omitnan"); %完整剖面
prof_Vm     = mean(V_median(:, xmask), 2, "omitnan");
prof_Vg     = mean(V_gaussian(:, xmask), 2, "omitnan");
prof_Vl     = mean(V_loess(:, xmask), 2, "omitnan");
prof_Vb     = mean(V_bilateral(:, xmask), 2, "omitnan");

prof_dVc    = mean(dVc_dy(:, xmask), 2, "omitnan");
prof_dVc2    = mean(dVc_dy(:, xmask2), 2, "omitnan"); %完整剖面用于对比
prof_dVm    = mean(dVm_dy(:, xmask), 2, "omitnan");
prof_dVg    = mean(dVg_dy(:, xmask), 2, "omitnan");
prof_dVl    = mean(dVl_dy(:, xmask), 2, "omitnan");
prof_dVb    = mean(dVb_dy(:, xmask), 2, "omitnan");

%% 可视化开始
figure;

%%%% ========================
%%%% 1. 第一行（Noisy）
%%%% ========================

%% ---- (1,1) Noisy velocity ----
subplot(5,4,1)
hV1 = imagesc(X(1,:),Y(:,1),V_noisy,[-3 3]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hV1,'AlphaData',~isnan(V_noisy));
title("Synth Ve (SR=5mm/yr;LD=20km) with white&correlated noise (sig=0.5+0.5mm/yr)");
colormap(vik); colorbar;
hold on;
rectangle('Position',[-10,-200,20,400],'EdgeColor','r','LineWidth',1.5);
rectangle('Position',[-150,50,101,101],'EdgeColor','g','LineWidth',1.5);

%% ---- (1,2) Profile ----
subplot(5,4,2)
plot(prof_y, prof_Vc,'r','LineWidth',1.5); hold on;
plot(prof_y, prof_V_noisy,'k'); grid on;
title("Profile Ve (noisy)");

%% ---- (1,3) Noisy strain ----
subplot(5,4,3)
hS1 = imagesc(X(1,:),Y(:,1),dVc_dy,[-100 100]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hS1,'AlphaData',~isnan(dVc_dy));
title("∂Ve/∂y Synth (without noise)"); colormap(vik); colorbar;
hold on;
rectangle('Position',[-10,-200,20,400],'EdgeColor','r','LineWidth',1.5);

%% ---- (1,4) Profile ----
subplot(5,4,4)
plot(prof_y, prof_dVc,'k'); grid on;
title("Profile ∂Ve/∂y Synth");
ylim([0 80]);

%%%% ========================
%%%% 2. Median
%%%% ========================

subplot(5,4,5)
hV2 = imagesc(X(1,:),Y(:,1),V_median,[-3 3]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hV2,'AlphaData',~isnan(V_median));
title(sprintf("Median (%d×%d)",spatial_winsize,spatial_winsize));

subplot(5,4,6)
plot(prof_y, prof_Vc,'r','LineWidth',1.5); hold on;
plot(prof_y, prof_Vm,'k'); grid on;
title("Profile Ve Median");

subplot(5,4,7)
hS2 = imagesc(X(1,:),Y(:,1),dVm_dy,[-100 100]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hS2,'AlphaData',~isnan(dVm_dy));
title("∂Ve/∂y Median");

subplot(5,4,8)
plot(prof_y, prof_dVc2,'r','LineWidth',1.5); hold on;
plot(prof_y, prof_dVm,'k'); grid on;
title("Profile ∂Ve/∂y Median");
ylim([0 80]);


%%%% ========================
%%%% 3. Gaussian
%%%% ========================

subplot(5,4,9)
hV3 = imagesc(X(1,:),Y(:,1),V_gaussian,[-3 3]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hV3,'AlphaData',~isnan(V_gaussian));
title(sprintf("Gaussian (%d×%d,σ_s=%d)",spatial_winsize,spatial_winsize,spatial_sigma));

subplot(5,4,10)
plot(prof_y, prof_Vc,'r','LineWidth',1.5); hold on;
plot(prof_y, prof_Vg,'k'); grid on;
title("Profile Ve Gaussian");

subplot(5,4,11)
hS3 = imagesc(X(1,:),Y(:,1),dVg_dy,[-100 100]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hS3,'AlphaData',~isnan(dVg_dy));
title("∂Ve/∂y Gaussian");

subplot(5,4,12)
plot(prof_y, prof_dVc2,'r','LineWidth',1.5); hold on;
plot(prof_y, prof_dVg,'k'); grid on;
title("Profile ∂Ve/∂y Gaussian");
ylim([0 80]);


%%%% ========================
%%%% 4. LOESS
%%%% ========================

subplot(5,4,13)
hV4 = imagesc(X(1,:),Y(:,1),V_loess,[-3 3]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hV4,'AlphaData',~isnan(V_loess));
title(sprintf("LOESS (%d×%d)",spatial_winsize,spatial_winsize));

subplot(5,4,14)
plot(prof_y, prof_Vc,'r','LineWidth',1.5); hold on;
plot(prof_y, prof_Vl,'k'); grid on;
title("Profile Ve LOESS");

subplot(5,4,15)
hS4 = imagesc(X(1,:),Y(:,1),dVl_dy,[-100 100]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hS4,'AlphaData',~isnan(dVl_dy));
title("∂Ve/∂y LOESS");

subplot(5,4,16)
plot(prof_y, prof_dVc2,'r','LineWidth',1.5); hold on;
plot(prof_y, prof_dVl,'k'); grid on;
title("Profile ∂Ve/∂y LOESS");
ylim([0 80]);


%%%% ========================
%%%% 5. Bilateral
%%%% ========================

subplot(5,4,17)
hV5 = imagesc(X(1,:),Y(:,1),V_bilateral,[-3 3]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hV5,'AlphaData',~isnan(V_bilateral));
title(sprintf("Bilateral (%d×%d,σ_s=%d,σ_v=%d)",...
       spatial_winsize,spatial_winsize,spatial_sigma,range_sigma));

subplot(5,4,18)
plot(prof_y, prof_Vc,'r','LineWidth',1.5); hold on;
plot(prof_y, prof_Vb,'k'); grid on;
title("Profile Ve Bilateral");

subplot(5,4,19)
hS5 = imagesc(X(1,:),Y(:,1),dVb_dy,[-100 100]); axis image xy;
set(gca,'Color',[0.7 0.7 0.7]);
set(hS5,'AlphaData',~isnan(dVb_dy));
title("∂Ve/∂y Bilateral");

subplot(5,4,20)
plot(prof_y, prof_dVc2,'r','LineWidth',1.5); hold on;
plot(prof_y, prof_dVb,'k'); grid on;
title("Profile ∂Ve/∂y Bilateral");
ylim([0 80]);
