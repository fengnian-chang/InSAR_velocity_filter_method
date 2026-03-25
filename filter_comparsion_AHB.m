clear all 
%% READ IN GEOTIFFS
[Vn_gnss, R_n] = readgeoraster('AHB/vn_AHB_final.tif','CoordinateSystemType','geographic');
[Ve_insar, R_e] = readgeoraster('AHB/ve_AHB_final.tif','CoordinateSystemType','geographic');

xcoords = linspace(R_e.LongitudeLimits(1),R_e.LongitudeLimits(2),R_e.RasterSize(2));
ycoords = linspace(R_e.LatitudeLimits(2),R_e.LatitudeLimits(1),R_e.RasterSize(1));

%% Filter parameters (in px)
spatial_sigma = 15;
winsize = 91;

% for bilateral
value_sigma = 3;    % velocity threshold

% Local small window used to calculate the reference velocity
center_winsize = 25;    % Generally good for center_winsize between 20 and 50

% center_sigma = 1;  % if use local gaussuan

%% Different Filter method (Bilateral, Median, Gaussian, and LOESS)
% Gaussian（nan-aware）
gauss_kernel = fspecial('gaussian', [winsize winsize], spatial_sigma);
V_gaussian = nanconv(Ve_insar, gauss_kernel, 'edge', 'nonanout');

% Median（nan-aware）
V_median = nanmedfilt2(Ve_insar, [winsize winsize], 'nonanout');

% LOESS（nan-aware）
% V_loess = loess_nanfit2_plane(V_noisy, winsize, 'nanout');

% Bilateral（NaN-aware）
V_bilateral = bilateral_nanconv(Ve_insar, spatial_sigma, value_sigma, winsize, center_winsize, 'median', 'nonanout');

%% Calculate velocity gradients
% GNSS velocity
[dvndx,dvndy] = gradient(Vn_gnss);
% InSAR velocity
[dvedx_gaussian,dvedy_gaussian] = gradient(V_gaussian);
[dvedx_median,dvedy_median] = gradient(V_median);
[dvedx_bilateral,dvedy_bilateral] = gradient(V_bilateral);

scaley = -1000/haversine(ycoords(1),xcoords(1),ycoords(2),xcoords(1));
scalex = zeros(length(ycoords),1);
for i = 1:length(ycoords)
    scalex(i) = 1000/haversine(ycoords(i),xcoords(1),ycoords(i),xcoords(2));
end

s_dvndy = dvndy*scaley;
s_dvndx = dvndx.*scalex;

s_dvedx_gaussian = dvedx_gaussian.*scalex;
s_dvedy_gaussian = dvedy_gaussian*scaley;
s_dvedx_median = dvedx_median.*scalex;
s_dvedy_median = dvedy_median*scaley;
s_dvedx_bilateral = dvedx_bilateral.*scalex;
s_dvedy_bilateral = dvedy_bilateral*scaley;

%% Read the fault data and prepare to plot
% gmt_file = 'fault/gem_active_faults.gmt';
% 
% fid = fopen(gmt_file,'r');
% faults = {};   % 用 cell 数组存放每一条断层
% current_x = [];
% current_y = [];
% 
% while true
%     line = fgetl(fid);
%     if ~ischar(line), break; end
% 
%     % --- 每段开始标识符 ---
%     if startsWith(line, '>')
%         % 保存上一段
%         if ~isempty(current_x)
%             faults{end+1} = [current_x(:), current_y(:)];
%         end
%         current_x = [];
%         current_y = [];
%         continue;
%     end
% 
%     % --- 解析 lon lat ---
%     vals = sscanf(line, '%f %f');
%     if numel(vals) == 2
%         current_x(end+1) = vals(1);
%         current_y(end+1) = vals(2);
%     end
% end
% fclose(fid);
% 
% % 最后一段加入
% if ~isempty(current_x)
%     faults{end+1} = [current_x(:), current_y(:)];
% end
% 
% % ==== 由速度场确定绘图经纬度范围 ====
% xmin = min(xcoords);
% xmax = max(xcoords);
% ymin = min(ycoords);
% ymax = max(ycoords);

%% Calculate and plot strain rate tensor components plus derived products
Eyy = s_dvndy;
Exx_gaussian = s_dvedx_gaussian;
Exx_median = s_dvedx_median;
Exx_bilateral = s_dvedx_bilateral;
Exy_gaussian = 1/2*(s_dvedy_gaussian+s_dvndx);
Exy_median = 1/2*(s_dvedy_median+s_dvndx);
Exy_bilateral = 1/2*(s_dvedy_bilateral+s_dvndx);

% InSAR GNSS Combine
Evort_gaussian = s_dvndx - s_dvedy_gaussian ;  % vorticity
Edil_gaussian = Exx_gaussian + Eyy; % dilatation
Eshear_gaussian = sqrt(Exy_gaussian.^2 + (Exx_gaussian - Eyy).^2/4); % Max shearstrain rate
EII_gaussian = sqrt(Exx_gaussian.^2 + 2*Exy_gaussian.^2 +Eyy.^2); % 2nd Invariant strain rate 

Evort_median = s_dvndx - s_dvedy_median ;  % vorticity
Edil_median = Exx_median + Eyy; % dilatation
Eshear_median = sqrt(Exy_median.^2 + (Exx_median - Eyy).^2/4); % Max shearstrain rate
EII_median = sqrt(Exx_median.^2 + 2*Exy_median.^2 +Eyy.^2); % 2nd Invariant strain rate 

Evort_bilateral = s_dvndx - s_dvedy_bilateral ;  % vorticity
Edil_bilateral = Exx_bilateral + Eyy; % dilatation
Eshear_bilateral = sqrt(Exy_bilateral.^2 + (Exx_bilateral - Eyy).^2/4); % Max shearstrain rate
EII_bilateral = sqrt(Exx_bilateral.^2 + 2*Exy_bilateral.^2 +Eyy.^2); % 2nd Invariant strain rate 

%% reblank out NaNs if needed
% nnz(isnan(Ve_insar));         % Check the number of NaN values
mask_valid = ~isnan(Ve_insar);   % mask using Original Ve

V_gaussian(~mask_valid)=nan; 
V_median(~mask_valid)=nan; 
V_bilateral(~mask_valid)=nan; 

s_dvedx_gaussian(~mask_valid)=nan; 
s_dvedx_median(~mask_valid)=nan; 
s_dvedx_bilateral(~mask_valid)=nan; 
s_dvedy_gaussian(~mask_valid)=nan; 
s_dvedy_median(~mask_valid)=nan; 
s_dvedy_bilateral(~mask_valid)=nan; 

Evort_gaussian(~mask_valid)=nan; 
Edil_gaussian(~mask_valid)=nan; 
Eshear_gaussian(~mask_valid)=nan; 
EII_gaussian(~mask_valid)=nan; 
Evort_median(~mask_valid)=nan; 
Edil_median(~mask_valid)=nan; 
Eshear_median(~mask_valid)=nan; 
EII_median(~mask_valid)=nan; 
Evort_bilateral(~mask_valid)=nan; 
Edil_bilateral(~mask_valid)=nan; 
Eshear_bilateral(~mask_valid)=nan; 
EII_bilateral(~mask_valid)=nan; 

%% plot velocity gradients
% % Load colormap
% vik = importdata('vik.mat');
% cpt.vik = vik;
% 
% figure(1)
% 
% subplot(3,3,1)
% imagesc(xcoords,ycoords,Eshear_gaussian, ...
%         'AlphaData', ~isnan(Eshear_gaussian));   % NaN 透明
% axis image
% axis xy
% caxis([-100 100])
% title('Shear gaussian (win:100km, sig:15km)')
% colormap(vik)
% set(gca,'Color',[1 1 1])                          % 背景（NaN 区域）纯白
% colorbar
% hold on
% for i = 1:length(faults)
%     lon = faults{i}(:,1);
%     lat = faults{i}(:,2);
% 
%     % 只保留在 bbox 内的点（外面的设为 NaN，这样 plot 会断开）
%     in = lon >= xmin & lon <= xmax & lat >= ymin & lat <= ymax;
%     if any(in)
%         lon_plot = lon;
%         lat_plot = lat;
%         lon_plot(~in) = NaN;
%         lat_plot(~in) = NaN;
%         plot(lon_plot, lat_plot, 'k--', 'LineWidth', 0.1);
%     end
% end
% 
% 
% subplot(3,3,2)
% imagesc(xcoords,ycoords,Eshear_median, ...
%         'AlphaData', ~isnan(Eshear_median));   % NaN 透明
% axis image
% axis xy
% caxis([-100 100])
% title('Shear median square (win:100km)')
% colormap(vik)
% set(gca,'Color',[1 1 1])                          % 背景（NaN 区域）纯白
% %colorbar
% hold on
% for i = 1:length(faults)
%     lon = faults{i}(:,1);
%     lat = faults{i}(:,2);
% 
%     % 只保留在 bbox 内的点（外面的设为 NaN，这样 plot 会断开）
%     in = lon >= xmin & lon <= xmax & lat >= ymin & lat <= ymax;
%     if any(in)
%         lon_plot = lon;
%         lat_plot = lat;
%         lon_plot(~in) = NaN;
%         lat_plot(~in) = NaN;
%         plot(lon_plot, lat_plot, 'k--', 'LineWidth', 0.1);
%     end
% end
% 
% % subplot(3,4,3)
% % imagesc(xcoords,ycoords,Eshear_median, ...
% %         'AlphaData', ~isnan(Eshear_median));   % NaN 透明
% % axis image
% % axis xy
% % caxis([-200 200])
% % title('Eshear median circle (win:150km)')
% % colormap(vik)
% % set(gca,'Color',[1 1 1])                          % 背景（NaN 区域）纯白
% % %colorbar
% % hold on
% % for i = 1:length(faults)
% %     lon = faults{i}(:,1);
% %     lat = faults{i}(:,2);
% % 
% %     % 只保留在 bbox 内的点（外面的设为 NaN，这样 plot 会断开）
% %     in = lon >= xmin & lon <= xmax & lat >= ymin & lat <= ymax;
% %     if any(in)
% %         lon_plot = lon;
% %         lat_plot = lat;
% %         lon_plot(~in) = NaN;
% %         lat_plot(~in) = NaN;
% %         plot(lon_plot, lat_plot, 'k--', 'LineWidth', 0.1);
% %     end
% % end
% 
% subplot(3,3,3)
% imagesc(xcoords,ycoords,Eshear_bilateral, ...
%         'AlphaData', ~isnan(Eshear_bilateral));   % NaN 透明
% axis image
% axis xy
% caxis([-100 100])
% title('Shear bilateral (win:100km, sig_s:15km, sig_v:3mm/yr)')
% colormap(vik)
% set(gca,'Color',[1 1 1])                          % 背景（NaN 区域）纯白
% %colorbar
% hold on
% for i = 1:length(faults)
%     lon = faults{i}(:,1);
%     lat = faults{i}(:,2);
% 
%     % 只保留在 bbox 内的点（外面的设为 NaN，这样 plot 会断开）
%     in = lon >= xmin & lon <= xmax & lat >= ymin & lat <= ymax;
%     if any(in)
%         lon_plot = lon;
%         lat_plot = lat;
%         lon_plot(~in) = NaN;
%         lat_plot(~in) = NaN;
%         plot(lon_plot, lat_plot, 'k--', 'LineWidth', 0.1);
%     end
% end
% 
% 
% subplot(3,3,4)
% imagesc(xcoords,ycoords,Edil_gaussian, ...
%         'AlphaData', ~isnan(Edil_gaussian));   % NaN 透明
% axis image
% axis xy
% caxis([-100 100])
% title('Dil gaussian')
% colormap(vik)
% set(gca,'Color',[1 1 1])                          % 背景（NaN 区域）纯白
% %colorbar
% hold on
% for i = 1:length(faults)
%     lon = faults{i}(:,1);
%     lat = faults{i}(:,2);
% 
%     % 只保留在 bbox 内的点（外面的设为 NaN，这样 plot 会断开）
%     in = lon >= xmin & lon <= xmax & lat >= ymin & lat <= ymax;
%     if any(in)
%         lon_plot = lon;
%         lat_plot = lat;
%         lon_plot(~in) = NaN;
%         lat_plot(~in) = NaN;
%         plot(lon_plot, lat_plot, 'k--', 'LineWidth', 0.1);
%     end
% end
% 
% 
% subplot(3,3,5)
% imagesc(xcoords,ycoords,Edil_median, ...
%         'AlphaData', ~isnan(Edil_median));   % NaN 透明
% axis image
% axis xy
% caxis([-100 100])
% title('Dil median square')
% colormap(vik)
% set(gca,'Color',[1 1 1])                          % 背景（NaN 区域）纯白
% %colorbar
% hold on
% for i = 1:length(faults)
%     lon = faults{i}(:,1);
%     lat = faults{i}(:,2);
% 
%     % 只保留在 bbox 内的点（外面的设为 NaN，这样 plot 会断开）
%     in = lon >= xmin & lon <= xmax & lat >= ymin & lat <= ymax;
%     if any(in)
%         lon_plot = lon;
%         lat_plot = lat;
%         lon_plot(~in) = NaN;
%         lat_plot(~in) = NaN;
%         plot(lon_plot, lat_plot, 'k--', 'LineWidth', 0.1);
%     end
% end
% 
% % subplot(3,4,7)
% % imagesc(xcoords,ycoords,Edil_median, ...
% %         'AlphaData', ~isnan(Edil_median));   % NaN 透明
% % axis image
% % axis xy
% % caxis([-200 200])
% % title('Edil median circle')
% % colormap(vik)
% % set(gca,'Color',[1 1 1])                          % 背景（NaN 区域）纯白
% % %colorbar
% % hold on
% % for i = 1:length(faults)
% %     lon = faults{i}(:,1);
% %     lat = faults{i}(:,2);
% % 
% %     % 只保留在 bbox 内的点（外面的设为 NaN，这样 plot 会断开）
% %     in = lon >= xmin & lon <= xmax & lat >= ymin & lat <= ymax;
% %     if any(in)
% %         lon_plot = lon;
% %         lat_plot = lat;
% %         lon_plot(~in) = NaN;
% %         lat_plot(~in) = NaN;
% %         plot(lon_plot, lat_plot, 'k--', 'LineWidth', 0.1);
% %     end
% % end
% 
% subplot(3,3,6)
% imagesc(xcoords,ycoords,Edil_bilateral, ...
%         'AlphaData', ~isnan(Edil_bilateral));   % NaN 透明
% axis image
% axis xy
% caxis([-100 100])
% title('Dil bilateral')
% colormap(vik)
% set(gca,'Color',[1 1 1])                          % 背景（NaN 区域）纯白
% %colorbar
% hold on
% for i = 1:length(faults)
%     lon = faults{i}(:,1);
%     lat = faults{i}(:,2);
% 
%     % 只保留在 bbox 内的点（外面的设为 NaN，这样 plot 会断开）
%     in = lon >= xmin & lon <= xmax & lat >= ymin & lat <= ymax;
%     if any(in)
%         lon_plot = lon;
%         lat_plot = lat;
%         lon_plot(~in) = NaN;
%         lat_plot(~in) = NaN;
%         plot(lon_plot, lat_plot, 'k--', 'LineWidth', 0.1);
%     end
% end
% 
% subplot(3,3,7)
% imagesc(xcoords,ycoords,Evort_gaussian, ...
%         'AlphaData', ~isnan(Evort_gaussian));   % NaN 透明
% axis image
% axis xy
% caxis([-100 100])
% title('Vort gaussian')
% colormap(vik)
% set(gca,'Color',[1 1 1])                          % 背景（NaN 区域）纯白
% %colorbar
% hold on
% for i = 1:length(faults)
%     lon = faults{i}(:,1);
%     lat = faults{i}(:,2);
% 
%     % 只保留在 bbox 内的点（外面的设为 NaN，这样 plot 会断开）
%     in = lon >= xmin & lon <= xmax & lat >= ymin & lat <= ymax;
%     if any(in)
%         lon_plot = lon;
%         lat_plot = lat;
%         lon_plot(~in) = NaN;
%         lat_plot(~in) = NaN;
%         plot(lon_plot, lat_plot, 'k--', 'LineWidth', 0.1);
%     end
% end
% 
% 
% subplot(3,3,8)
% imagesc(xcoords,ycoords,Evort_median, ...
%         'AlphaData', ~isnan(Evort_median));   % NaN 透明
% axis image
% axis xy
% caxis([-100 100])
% title('Vort median square')
% colormap(vik)
% set(gca,'Color',[1 1 1])                          % 背景（NaN 区域）纯白
% %colorbar
% hold on
% for i = 1:length(faults)
%     lon = faults{i}(:,1);
%     lat = faults{i}(:,2);
% 
%     % 只保留在 bbox 内的点（外面的设为 NaN，这样 plot 会断开）
%     in = lon >= xmin & lon <= xmax & lat >= ymin & lat <= ymax;
%     if any(in)
%         lon_plot = lon;
%         lat_plot = lat;
%         lon_plot(~in) = NaN;
%         lat_plot(~in) = NaN;
%         plot(lon_plot, lat_plot, 'k--', 'LineWidth', 0.1);
%     end
% end
% 
% % subplot(3,4,11)
% % imagesc(xcoords,ycoords,Evort_median, ...
% %         'AlphaData', ~isnan(Evort_median));   % NaN 透明
% % axis image
% % axis xy
% % caxis([-200 200])
% % title('Evort median circle')
% % colormap(vik)
% % set(gca,'Color',[1 1 1])                          % 背景（NaN 区域）纯白
% % %colorbar
% % hold on
% % for i = 1:length(faults)
% %     lon = faults{i}(:,1);
% %     lat = faults{i}(:,2);
% % 
% %     % 只保留在 bbox 内的点（外面的设为 NaN，这样 plot 会断开）
% %     in = lon >= xmin & lon <= xmax & lat >= ymin & lat <= ymax;
% %     if any(in)
% %         lon_plot = lon;
% %         lat_plot = lat;
% %         lon_plot(~in) = NaN;
% %         lat_plot(~in) = NaN;
% %         plot(lon_plot, lat_plot, 'k--', 'LineWidth', 0.1);
% %     end
% % end
% 
% subplot(3,3,9)
% imagesc(xcoords,ycoords,Evort_bilateral, ...
%         'AlphaData', ~isnan(Evort_bilateral));   % NaN 透明
% axis image
% axis xy
% caxis([-100 100])
% title('Vort bilateral')
% colormap(vik)
% set(gca,'Color',[1 1 1])                          % 背景（NaN 区域）纯白
% %colorbar
% hold on
% for i = 1:length(faults)
%     lon = faults{i}(:,1);
%     lat = faults{i}(:,2);
% 
%     % 只保留在 bbox 内的点（外面的设为 NaN，这样 plot 会断开）
%     in = lon >= xmin & lon <= xmax & lat >= ymin & lat <= ymax;
%     if any(in)
%         lon_plot = lon;
%         lat_plot = lat;
%         lon_plot(~in) = NaN;
%         lat_plot(~in) = NaN;
%         plot(lon_plot, lat_plot, 'k--', 'LineWidth', 0.1);
%     end
% end
%% write to tiff
outdir = 'AHB';
% dvedx
geotiffwrite(fullfile(outdir, strcat('dvedx_gaussian_win',num2str(winsize),'_sig_s',num2str(spatial_sigma),'.tif')), s_dvedx_gaussian, R_e);
geotiffwrite(fullfile(outdir, strcat('dvedx_bilateral_win',num2str(winsize),'_sig_s',num2str(spatial_sigma),'_sig_v',num2str(value_sigma),'_median',num2str(center_winsize),'.tif')), s_dvedx_bilateral, R_e);
geotiffwrite(fullfile(outdir, strcat('dvedx_median_win',num2str(winsize),'.tif')), s_dvedx_median, R_e);
% dvedy
geotiffwrite(fullfile(outdir, strcat('dvedy_gaussian_win',num2str(winsize),'_sig_s',num2str(spatial_sigma),'.tif')), s_dvedy_gaussian, R_e);
geotiffwrite(fullfile(outdir, strcat('dvedy_bilateral_win',num2str(winsize),'_sig_s',num2str(spatial_sigma),'_sig_v',num2str(value_sigma),'_median',num2str(center_winsize),'.tif')), s_dvedy_bilateral, R_e);
geotiffwrite(fullfile(outdir, strcat('dvedy_median_win',num2str(winsize),'.tif')), s_dvedy_median, R_e);
% Max_shear
geotiffwrite(fullfile(outdir, strcat('Max_shear_gaussian_win',num2str(winsize),'_sig_s',num2str(spatial_sigma),'.tif')), Eshear_gaussian, R_e);
geotiffwrite(fullfile(outdir, strcat('Max_shear_bilateral_win',num2str(winsize),'_sig_s',num2str(spatial_sigma),'_sig_v',num2str(value_sigma),'_median',num2str(center_winsize),'.tif')), Eshear_bilateral, R_e);
geotiffwrite(fullfile(outdir, strcat('Max_shear_median_win',num2str(winsize),'.tif')), Eshear_median, R_e);
% Dil
geotiffwrite(fullfile(outdir, strcat('Dil_gaussian_win',num2str(winsize),'_sig_s',num2str(spatial_sigma),'.tif')), Edil_gaussian, R_e);
geotiffwrite(fullfile(outdir, strcat('Dil_bilateral_win',num2str(winsize),'_sig_s',num2str(spatial_sigma),'_sig_v',num2str(value_sigma),'_median',num2str(center_winsize),'.tif')), Edil_bilateral, R_e);
geotiffwrite(fullfile(outdir, strcat('Dil_median_win',num2str(winsize),'.tif')), Edil_median, R_e);
% Vort
geotiffwrite(fullfile(outdir, strcat('Vort_gaussian_win',num2str(winsize),'_sig_s',num2str(spatial_sigma),'.tif')), Evort_gaussian, R_e);
geotiffwrite(fullfile(outdir, strcat('Vort_bilateral_win',num2str(winsize),'_sig_s',num2str(spatial_sigma),'_sig_v',num2str(value_sigma),'_median',num2str(center_winsize),'.tif')), Evort_bilateral, R_e);
geotiffwrite(fullfile(outdir, strcat('Vort_median_win',num2str(winsize),'.tif')), Evort_median, R_e);
% II
geotiffwrite(fullfile(outdir, strcat('II_gaussian_win',num2str(winsize),'_sig_s',num2str(spatial_sigma),'.tif')), EII_gaussian, R_e);
geotiffwrite(fullfile(outdir, strcat('II_bilateral_win',num2str(winsize),'_sig_s',num2str(spatial_sigma),'_sig_v',num2str(value_sigma),'_median',num2str(center_winsize),'.tif')), EII_bilateral, R_e);
geotiffwrite(fullfile(outdir, strcat('II_median_win',num2str(winsize),'.tif')), EII_median, R_e);
disp(['Results saved to path']);