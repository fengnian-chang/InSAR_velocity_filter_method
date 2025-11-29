function direct_strain_calc (sigma, winsize)
%% Script to calculate strain rates (for Tianshan) directly from the InSAR resolution files
%
% TJW 17 June 2023
% J. Fang 13 Oct 2023
% J. Fang 6 Mar 2025
% FN Chang 16 Sept 2025

% 先执行filter_velmap_vn_ve.sh滤波

%close all
%clear all

% direct_strain_calc(20,100);
% direct_strain_calc(30,200)
%% READ IN GEOTIFFS

[Vn_gnss, R_n] = readgeoraster('Tianshan_strain/GNSS_N_vel_resample.tif','CoordinateSystemType','geographic');
[Ve_gnss, R_e] = readgeoraster('Tianshan_strain/GNSS_E_vel_resample.tif','CoordinateSystemType','geographic');
[Ve_insar, R_e] = readgeoraster('Tianshan_strain/InSAR_ve_1km.tif','CoordinateSystemType','geographic');

%[Ve, R_e] = readgeoraster('Tianshan_strain/Tianshan_ref2-2_decomp0_vE.geo.tif','CoordinateSystemType','geographic');
%Ve=Ve(1:size(Vn,1),1:size(Vn,2)); % cut to same size as Vn

xcoords = linspace(R_e.LongitudeLimits(1),R_e.LongitudeLimits(2),R_e.RasterSize(2));
ycoords = linspace(R_e.LatitudeLimits(2),R_e.LatitudeLimits(1),R_e.RasterSize(1));

%% Gaussian Filter (ignore Nan)

%sigma = [2.5,5,7.5,10];  %in pixels
%sigma = 7.5;  %in pixels
%winsize = 50; %winsize

disp(['Sigma: ', num2str(sigma)]);
disp(['WinSize: ', num2str(winsize)]);

gauss_kernel = fspecial('gaussian',[winsize,winsize], sigma);
%gauss_kernel2 = fspecial('gaussian',[winsize,winsize], sigma(2));
%gauss_kernel3 = fspecial('gaussian',[winsize,winsize], sigma(3));
%gauss_kernel4 = fspecial('gaussian',[winsize,winsize], sigma(4));

Ve_insar_filt = nanconv(Ve_insar,gauss_kernel);
%Ve_filt2 = nanconv(Ve,gauss_kernel2);
%Ve_filt3 = nanconv(Ve,gauss_kernel3);
%Ve_filt4 = nanconv(Ve,gauss_kernel4);

%% Calculate velocity gradients
[dvedx_gnss,dvedy_gnss] = gradient(Ve_gnss);
[dvedx_insar,dvedy_insar] = gradient(Ve_insar_filt);
%[dvedx2,dvedy2] = gradient(Ve_filt2);
%[dvedx3,dvedy3] = gradient(Ve_filt3);
%[dvedx4,dvedy4] = gradient(Ve_filt4);

%Using v3 (Currently sigma = 7.5) as final version
%dvedx = dvedx3;
%dvedy = dvedy3;

[dvndx,dvndy] = gradient(Vn_gnss);

%% Scaling of velocity gradients
scaley = -1000/haversine(ycoords(1),xcoords(1),ycoords(2),xcoords(1));
scalex = zeros(length(ycoords),1);
for i = 1:length(ycoords)
    scalex(i) = 1000/haversine(ycoords(i),xcoords(1),ycoords(i),xcoords(2));
end

s_dvedy_insar = dvedy_insar*scaley;
s_dvedx_insar = dvedx_insar.*scalex;
s_dvedy_gnss = dvedy_gnss*scaley;
s_dvedx_gnss = dvedx_gnss.*scalex;
s_dvndy = dvndy*scaley;
s_dvndx = dvndx.*scalex;

%% plot velocity gradients
figure(1)
subplot(2,2,1)
 imagesc(xcoords,ycoords,s_dvedy_insar)
 axis image
 axis xy
 caxis([-100 100])
 title('InSAR dVe/dy')
subplot(2,2,2)
 imagesc(xcoords,ycoords,s_dvedx_insar)
 axis image
 axis xy
 caxis([-100 100])
 title('InSAR dVe/dx')
subplot(2,2,3)
 imagesc(xcoords,ycoords,s_dvndx)
 axis image
 axis xy
 caxis([-100 100])
 title('dVn/dx')
subplot(2,2,4)
 imagesc(xcoords,ycoords,s_dvndy)
 axis image
 axis xy
 caxis([-100 100])
 title('dVn/dy')

 %% Calculate and plot strain rate tensor components plus derived products
Exx_insar = s_dvedx_insar;
Exx_gnss = s_dvedx_gnss;
Eyy = s_dvndy;
Exy_insar = 1/2*(s_dvedy_insar+s_dvndx);
Exy_gnss = 1/2*(s_dvedy_gnss+s_dvndx);

% InSAR GNSS Combine
Evort_insar = s_dvndx - s_dvedy_insar ;  % vorticity
Edil_insar = Exx_insar + Eyy; % dilatation
Eshear_insar = sqrt(Exy_insar.^2 + (Exx_insar - Eyy).^2/4); % Max shearstrain rate
EII_insar = sqrt(Exx_insar.^2 + 2*Exy_insar.^2 +Eyy.^2); % 2nd Invariant strain rate 

% GNSS only
Evort_gnss = s_dvndx - s_dvedy_gnss ;  % vorticity
Edil_gnss = Exx_gnss + Eyy; % dilatation
Eshear_gnss = sqrt(Exy_gnss.^2 + (Exx_gnss - Eyy).^2/4); % Max shearstrain rate
EII_gnss = sqrt(Exx_gnss.^2 + 2*Exy_gnss.^2 +Eyy.^2); % 2nd Invariant strain rate 

Eshear_insar(isnan(Ve_insar))=nan;  % reblank out NaNs if needed
Edil_insar(isnan(Ve_insar))=nan;
EII_insar(isnan(Ve_insar))=nan;
Evort_insar(isnan(Ve_insar))=nan;
Eshear_gnss(isnan(Ve_insar))=nan;
Edil_gnss(isnan(Ve_insar))=nan;
EII_gnss(isnan(Ve_insar))=nan;
Evort_gnss(isnan(Ve_insar))=nan;

% Difference
Eshear_diff = Eshear_insar - Eshear_gnss;
Edil_diff = Edil_insar - Edil_gnss;
EII_diff = EII_insar - EII_gnss;

% InSAR strain plot
figure(2)
subplot(3,2,1)
 imagesc(xcoords,ycoords,Exx_insar)
 axis image
 axis xy
 caxis([-100 100])
 title('InSAR Exx')
subplot(3,2,3)
 imagesc(xcoords,ycoords,Eyy)
 axis image
 axis xy
 caxis([-100 100])
 title('Eyy')
subplot(3,2,5)
 imagesc(xcoords,ycoords,Exy_insar)
 axis image
 axis xy
 caxis([-100 100])
 title('InSAR Exy')
subplot(3,2,2)
 imagesc(xcoords,ycoords,Edil_insar)
 axis image
 axis xy
 caxis([-100 100])
 title('InSAR Edil')
subplot(3,2,4)
 imagesc(xcoords,ycoords,Eshear_insar)
 axis image
 axis xy
 caxis([0 100])
 title('InSAR Emaxshear')
subplot(3,2,6)
 imagesc(xcoords,ycoords,EII_insar)
 axis image
 axis xy
 caxis([0 100])
 title('InSAR EII')

% GNSS strain plot
figure(3)
subplot(3,2,1)
 imagesc(xcoords,ycoords,Exx_gnss)
 axis image
 axis xy
 caxis([-100 100])
 title('GNSS Exx')
subplot(3,2,3)
 imagesc(xcoords,ycoords,Eyy)
 axis image
 axis xy
 caxis([-100 100])
 title('Eyy')
subplot(3,2,5)
 imagesc(xcoords,ycoords,Exy_gnss)
 axis image
 axis xy
 caxis([-100 100])
 title('GNSS Exy')
subplot(3,2,2)
 imagesc(xcoords,ycoords,Edil_gnss)
 axis image
 axis xy
 caxis([-100 100])
 title('GNSS Edil')
subplot(3,2,4)
 imagesc(xcoords,ycoords,Eshear_gnss)
 axis image
 axis xy
 caxis([0 100])
 title('GNSS Emaxshear')
subplot(3,2,6)
 imagesc(xcoords,ycoords,EII_gnss)
 axis image
 axis xy
 caxis([0 100])
 title('GNSS EII')

% compare
figure(4)
subplot(2,2,1)
 imagesc(xcoords,ycoords,Eshear_gnss)
 axis image
 axis xy
 caxis([0 100])
 title('GNSS Emaxshear')
subplot(2,2,2)
 imagesc(xcoords,ycoords,Eshear_insar)
 axis image
 axis xy
 caxis([0 100])
 title('InSAR Emaxshear')
subplot(2,2,3)
 imagesc(xcoords,ycoords,Eshear_diff)
 axis image
 axis xy
 caxis([0 100])
 title('Diff Emaxshear')

figure(5)
subplot(2,2,1)
 imagesc(xcoords,ycoords,Edil_gnss)
 axis image
 axis xy
 caxis([0 100])
 title('GNSS Edil')
subplot(2,2,2)
 imagesc(xcoords,ycoords,Edil_insar)
 axis image
 axis xy
 caxis([0 100])
 title('InSAR Edil')
subplot(2,2,3)
 imagesc(xcoords,ycoords,Edil_diff)
 axis image
 axis xy
 caxis([0 100])
 title('Diff Edil')

 %%
outdir = 'Tianshan_strain/';

geotiffwrite(fullfile(outdir, strcat('InSAR_Exx_sig',num2str(sigma),'_win',num2str(winsize),'.tif')), Exx_insar, R_n);
disp(['Results saved to  ', fullfile(outdir, strcat('Exx_sig',num2str(sigma),'_win',num2str(winsize),'.tif'))]);

geotiffwrite(fullfile(outdir, strcat('InSAR_Exy_sig',num2str(sigma),'_win',num2str(winsize),'.tif')), Exy_insar, R_n);
disp(['Results saved to  ', fullfile(outdir, strcat('Exy_sig',num2str(sigma),'_win',num2str(winsize),'.tif'))]);

geotiffwrite(fullfile(outdir, strcat('InSAR_Edil_sig',num2str(sigma),'_win',num2str(winsize),'.tif')), Edil_insar, R_n);
disp(['Results saved to  ', fullfile(outdir, strcat('Edil_sig',num2str(sigma),'_win',num2str(winsize),'.tif'))]);

geotiffwrite(fullfile(outdir, strcat('InSAR_Eshear_sig',num2str(sigma),'_win',num2str(winsize),'.tif')), Eshear_insar, R_n);
disp(['Results saved to  ', fullfile(outdir, strcat('Eshear_sig',num2str(sigma),'_win',num2str(winsize),'.tif'))]);

geotiffwrite(fullfile(outdir, strcat('InSAR_EII_sig',num2str(sigma),'_win',num2str(winsize),'.tif')), EII_insar, R_n);
disp(['Results saved to  ', fullfile(outdir, strcat('EII_sig',num2str(sigma),'_win',num2str(winsize),'.tif'))]);

geotiffwrite(fullfile(outdir, strcat('InSAR_Evort_sig',num2str(sigma),'_win',num2str(winsize),'.tif')), Evort_insar, R_n);
disp(['Results saved to  ', fullfile(outdir, strcat('Evort_sig',num2str(sigma),'_win',num2str(winsize),'.tif'))]);

% gnss output

geotiffwrite(fullfile(outdir, 'GNSS_Exx.tif'), Exx_gnss, R_n);
geotiffwrite(fullfile(outdir, 'GNSS_Exy.tif'), Exy_gnss, R_n);
geotiffwrite(fullfile(outdir, 'GNSS_Eyy.tif'), Eyy, R_n);
geotiffwrite(fullfile(outdir, 'GNSS_Edil.tif'), Edil_gnss, R_n);
geotiffwrite(fullfile(outdir, 'GNSS_Eshear.tif'), Eshear_gnss, R_n);
geotiffwrite(fullfile(outdir, 'GNSS_EII.tif'), EII_gnss, R_n);
geotiffwrite(fullfile(outdir, 'GNSS_Evort.tif'), Evort_gnss, R_n);
disp(['GNSS Results saved']);

geotiffwrite(fullfile(outdir, 'Edil_diff.tif'), Edil_diff, R_n);
geotiffwrite(fullfile(outdir, 'Eshear_diff.tif'), Eshear_diff, R_n);
geotiffwrite(fullfile(outdir, 'EII_diff.tif'), EII_diff, R_n);
disp(['Differences saved']);

end
