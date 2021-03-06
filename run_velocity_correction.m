clear all
%close all
var_save_flag = 1;
plot_save_flag = 1;

%%%%% Things to try???
% -stricter thresholds for debris to get more rain flags
%   -add Zh/bulk rhoHV threshold/dependency for filtering
%   -add SNR threshold for filtering
%       -take max value of PSD and correct points within 10/20/?? dB
% -look at more individual spectra --- ginput to click point on figure and
% retrieve x,y coords [x,y] = ginput(1)
%
% -test real data from David


% Good correction examples: twocell d3,n10000 el=2.0
%   r=15,az=21
%   r=12,az=21
%   r=12,az=19

% Plot idea:
% combine fig 8 plots to more easily compare size, concentration, elevation
% d3 1000000 vs d3 10000 at 0.5
% d3 1000000 vs d3 1000000 at 0.5
% d3 1000000 at 0.5 vs d3 1000000 at 5.0
% do this for PPIs and histograms


% % Load IQ data
LES = 'twocell';
base_dir = ['~/Documents/sims/' LES];
sim_dir = uigetdir(base_dir); % Location of SimRadar output files
filename = blib('choosefile', sim_dir, 'sim-*.mat');
load(filename);

% Separate path name and file name for saving files later
fname = erase(filename, sim_dir);
fname = fname(2:end-4);
dind = regexp(fname,'d');
nind = regexp(fname,'n');
dtype = fname(dind+1: nind-1);
dnum = fname(nind+1: end);
eind = regexp(fname,'PPI');
elev = fname(eind+3:eind+5);

img_dir = ['~/Documents/imgs/' LES '/d' dtype '/DPSD/'];
if ~exist(img_dir,'dir')
    mkdir(img_dir)
    addpath(genpath(img_dir))
end
savepath


% Put all the bulk variables from checkiq.m into one structure
params = struct('zdr', zdr, 'phv', rhohv, 'xx', xx, 'yy', yy, 'zz', zz,...
    'va', dat.params.va);

% if simulation contains multiple sweeps
if size(iqh,4) > 1
    iqh = iqh(:,:,:,1);
    iqv = iqv(:,:,:,1);
end

% Rearrange dimensions to run through DPSD code
iqh = permute(iqh, [3 1 2]);
iqv = permute(iqv, [3 1 2]);
M = size(iqh, 1);
va = dat.params.va;
vvx = flip(linspace(-va+1, va, M));
d = nuttallwin(M); % data window

% Uncorrected PSD and velocity
psd_old = fftshift(abs(fft(iqh.*d, M, 1)).^2, 1);
vr_old = vr;

% Calculate DPSDs
[szdr, sphv, svar, spsd] = dpsd_calc(iqh, iqv, d, M, 20, 1);

%%
% HCA debris classification

[obj_class, agg, memf] = hca_class(szdr, sphv, svar);
%obj_class = movmean(obj_class, 3, 1, 'omitnan', 'Endpoints', 'shrink');


% Filter debris from PSD
psd_filt = psd_old;
obj_class(1) = 0; % can't remember why I did this tbh but I didn't want endpoints marked as debris?
obj_class(end) = 0;
psd_filt(obj_class>0.33) = 0; % set all debris coefficients to 0

% Velocity recalculation -- need to make this more efficient than a loop!
var_norm = svar.phv ./ max(svar.phv,1); % normalize variance
b1 = 6; % aggressiveness of filter
var_inv = (1 - var_norm).^b1; % invert normalized variance
psd_new_var = var_inv .* psd_old;

b2 = 4; % aggressiveness of filter
psd_new_agg = (agg.rain.^b2) .* psd_old;

vr_new = zeros(size(iqh,2), size(iqh,3));
vr_new_var = vr_new;
vr_new_agg = vr_new;
for i = 1:size(iqh,2)
    for j = 1:size(iqh,3)
        vr_temp = calc_velocity_spectral(psd_filt(:,i,j), M, dat.params.prt, dat.params.lambda);
        if length(vr_temp) > 1
            vr_new(i,j) = vr_old(i,j);
        else
            vr_new(i,j) = vr_temp;
        end
        
        vr_new_var(i,j) = calc_velocity_spectral(psd_new_var(:,i,j), M, dat.params.prt, dat.params.lambda);
        vr_new_agg(i,j) = calc_velocity_spectral(psd_new_agg(:,i,j), M, dat.params.prt, dat.params.lambda);
    end
end

clear vr

% Load rain file for comparison
% sim_dir = uigetdir(base_dir); % Location of SimRadar files
% filename = blib('choosefile', sim_dir, '*.mat');
rain_dir = strrep(sim_dir, ['debris' dtype], 'nodebris');
rain_file = [rain_dir '/' strrep(fname, ['d' dtype 'n' dnum], 'nodebris')];
load([rain_file '.mat']);

vr_truth = vr;

%% Plot Doppler velocities


% Compare velocity recalculations
% figure(1)
% clf
% 
% ha = subplot(2,2,1);
% hs = pcolor(xx, yy, vr_truth(:,:,1));
% caxis([-1 1] * va)
% colormap(ha, blib('rgmap2'))
% shading flat
% c = colorbar;
% c.Label.String = 'm s^{-1}';
% c.Label.FontSize = 13;
% c.Label.VerticalAlignment = 'middle';
% title('(a) v_{truth}', 'FontSize', 14)
% xlabel('x (m)', 'FontSize', 14)
% ylabel('y (m)', 'FontSize', 14)
% 
% ha(2) = subplot(2,2,2);
% hs(2) = pcolor(xx, yy, vr_new(:,:,1));
% caxis([-1 1] * va)
% colormap(ha(2), blib('rgmap2'))
% shading flat
% c = colorbar;
% c.Label.String = 'm s^{-1}';
% c.Label.FontSize = 13;
% c.Label.VerticalAlignment = 'middle';
% title('(b) v_{corr} (DCA-based)', 'FontSize', 14)
% xlabel('x (m)', 'FontSize', 14)
% ylabel('y (m)', 'FontSize', 14)
% 
% ha(3) = subplot(2,2,3);
% hs(3) = pcolor(xx, yy, vr_new_var(:,:,1));
% caxis([-1 1] * va)
% colormap(ha(3), blib('rgmap2'))
% shading flat
% c = colorbar;
% c.Label.String = 'm s^{-1}';
% c.Label.FontSize = 13;
% c.Label.VerticalAlignment = 'middle';
% title('(c) v_{corr} (variance-based)', 'FontSize', 14)
% xlabel('x (m)', 'FontSize', 14)
% ylabel('y (m)', 'FontSize', 14)
% 
% ha(4) = subplot(2,2,4);
% hs(4) = pcolor(xx, yy, vr_new_agg(:,:,1));
% caxis([-1 1] * va)
% colormap(ha(4), blib('rgmap2'))
% shading flat
% c = colorbar;
% c.Label.String = 'm s^{-1}';
% c.Label.FontSize = 13;
% c.Label.VerticalAlignment = 'middle';
% title('(d) v_{corr} (aggregation-based)', 'FontSize', 14)
% xlabel('x (m)', 'FontSize', 14)
% ylabel('y (m)', 'FontSize', 14)
% 
% set(gcf, 'Units', 'inches', 'Position', [10 10 12 8])
% axes('Unit', 'Normalized', 'Position', [0.5 0.5 0.01 0.01])
% title(['Debris type ' dtype ', n=' dnum ', tilt=' elev '^o'], 'FontSize', 12);
% axis off
% 
% if plot_save_flag
%     print([img_dir 'vcorr-all-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
% end


figure(1)
clf

ha = subplot(3,2,[1,2]);
hs = pcolor(xx, yy, vr_truth(:,:,1));
caxis([-1 1] * va)
colormap(ha, blib('rgmap2'))
shading flat
c = colorbar;
c.Label.String = 'm s^{-1}';
c.Label.FontSize = 13;
c.Label.VerticalAlignment = 'middle';
title('(a) v_{truth}', 'FontSize', 14)
xlabel('x (m)', 'FontSize', 14)
ylabel('y (m)', 'FontSize', 14)
ha.Position = [0.3575 0.7093 0.2850 0.2157];

ha(2) = subplot(3,2,3);
hs(2) = pcolor(xx, yy, vr_old(:,:,1));
caxis([-1 1] * va)
colormap(ha(2), blib('rgmap2'))
shading flat
c = colorbar;
c.Label.String = 'm s^{-1}';
c.Label.FontSize = 13;
c.Label.VerticalAlignment = 'middle';
title('(c) Uncorrected v_r', 'FontSize', 14)
xlabel('x (m)', 'FontSize', 14)
ylabel('y (m)', 'FontSize', 14)

ha(3) = subplot(3,2,4);
hs(3) = pcolor(xx, yy, vr_new(:,:,1));
caxis([-1 1] * va)
colormap(ha(3), blib('rgmap2'))
shading flat
c = colorbar;
c.Label.String = 'm s^{-1}';
c.Label.FontSize = 13;
c.Label.VerticalAlignment = 'middle';
title('(c) DCA-corrected v_r', 'FontSize', 14)
xlabel('x (m)', 'FontSize', 14)
ylabel('y (m)', 'FontSize', 14)

ha(4) = subplot(3,2,5);
hs(4) = pcolor(xx, yy, vr_new_var(:,:,1));
caxis([-1 1] * va)
colormap(ha(4), blib('rgmap2'))
shading flat
c = colorbar;
c.Label.String = 'm s^{-1}';
c.Label.FontSize = 13;
c.Label.VerticalAlignment = 'middle';
title('(d) Variance-corrected v_r', 'FontSize', 14)
xlabel('x (m)', 'FontSize', 14)
ylabel('y (m)', 'FontSize', 14)

ha(5) = subplot(3,2,6);
hs(5) = pcolor(xx, yy, vr_new_agg(:,:,1));
caxis([-1 1] * va)
colormap(ha(5), blib('rgmap2'))
shading flat
c = colorbar;
c.Label.String = 'm s^{-1}';
c.Label.FontSize = 13;
c.Label.VerticalAlignment = 'middle';
title('(e) Aggregation-corrected v_r', 'FontSize', 14)
xlabel('x (m)', 'FontSize', 14)
ylabel('y (m)', 'FontSize', 14)

set(gcf, 'Units', 'inches', 'Position', [10 10 12 12])
axes('Unit', 'Normalized', 'Position', [0.5 0.03 0.01 0.01])
title(['Debris type ' dtype ', n=' dnum ', tilt=' elev '^o'], 'FontSize', 14);
axis off

if plot_save_flag
    print([img_dir 'vcorr-all-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
end



% 
% figure(2)
% 
% ha = subplot(2,2,1);
% hs = pcolor(xx, yy, vr_new(:,:,1));
% caxis([-1 1] * va)
% colormap(ha, blib('rgmap2'))
% shading flat
% axis square
% c = colorbar;
% c.Label.String = 'm s^{-1}';
% c.Label.FontSize = 13;
% c.Label.VerticalAlignment = 'middle';
% title('(a) v_{corrected}', 'FontSize', 14)
% xlabel('x (m)', 'FontSize', 14)
% ylabel('y (m)', 'FontSize', 14)
% 
% ha(2) = subplot(2,2,2);
% hs(2) = pcolor(xx, yy, vr_old(:,:,1));
% caxis([-1 1] * va)
% colormap(ha(2), blib('rgmap2'))
% shading flat
% axis square
% c = colorbar;
% c.Label.String = 'm s^{-1}';
% c.Label.FontSize = 13;
% c.Label.VerticalAlignment = 'middle';
% title('(b) v_{original}', 'FontSize', 14)
% xlabel('x (m)', 'FontSize', 14)
% ylabel('y (m)', 'FontSize', 14)
% 
% ha(3) = subplot(2,2,[3,4]);
% hs(3) = pcolor(xx, yy, vr_truth(:,:,1));
% caxis([-1 1] * va)
% colormap(ha(3), blib('rgmap2'))
% shading flat
% axis square
% c = colorbar;
% c.Label.String = 'm s^{-1}';
% c.Label.FontSize = 13;
% c.Label.VerticalAlignment = 'middle';
% title('(c) v_{truth}', 'FontSize', 14)
% xlabel('x (m)', 'FontSize', 14)
% ylabel('y (m)', 'FontSize', 14)
% 
% sgtitle('Debris-removed correction', 'FontSize', 14)
% 
% if plot_save_flag
%     print([img_dir 'vcorr-filter'], '-dpng')
% end



% figure(3)
% 
% ha = subplot(2,2,1);
% hs = pcolor(xx, yy, vr_new_var(:,:,1));
% caxis([-1 1] * va)
% colormap(ha, blib('rgmap2'))
% shading flat
% axis square
% c = colorbar;
% c.Label.String = 'm s^{-1}';
% c.Label.FontSize = 13;
% c.Label.VerticalAlignment = 'middle';
% title('(a) v_{corrected}', 'FontSize', 14)
% xlabel('x (m)', 'FontSize', 14)
% ylabel('y (m)', 'FontSize', 14)
% 
% ha(2) = subplot(2,2,2);
% hs(2) = pcolor(xx, yy, vr_old(:,:,1));
% caxis([-1 1] * va)
% colormap(ha(2), blib('rgmap2'))
% shading flat
% axis square
% c = colorbar;
% c.Label.String = 'm s^{-1}';
% c.Label.FontSize = 13;
% c.Label.VerticalAlignment = 'middle';
% title('(b) v_{original}', 'FontSize', 14)
% xlabel('x (m)', 'FontSize', 14)
% ylabel('y (m)', 'FontSize', 14)
% 
% ha(3) = subplot(2,2,[3,4]);
% hs(3) = pcolor(xx, yy, vr_truth(:,:,1));
% caxis([-1 1] * va)
% colormap(ha(3), blib('rgmap2'))
% shading flat
% axis square
% c = colorbar;
% c.Label.String = 'm s^{-1}';
% c.Label.FontSize = 13;
% c.Label.VerticalAlignment = 'middle';
% title('(c) v_{truth}', 'FontSize', 14)
% xlabel('x (m)', 'FontSize', 14)
% ylabel('y (m)', 'FontSize', 14)
% 
% sgtitle('Variance-weighted correction', 'FontSize', 14)
% 
% if plot_save_flag
%     print([img_dir 'vcorr-variance'], '-dpng')
% end



% figure(4)
% 
% ha = subplot(2,2,1);
% hs = pcolor(xx, yy, vr_new_agg(:,:,1));
% caxis([-1 1] * va)
% colormap(ha, blib('rgmap2'))
% shading flat
% axis square
% c = colorbar;
% c.Label.String = 'm s^{-1}';
% c.Label.FontSize = 13;
% c.Label.VerticalAlignment = 'middle';
% title('(a) v_{corrected}', 'FontSize', 14)
% xlabel('x (m)', 'FontSize', 14)
% ylabel('y (m)', 'FontSize', 14)
% 
% ha(2) = subplot(2,2,2);
% hs(2) = pcolor(xx, yy, vr_old(:,:,1));
% caxis([-1 1] * va)
% colormap(ha(2), blib('rgmap2'))
% shading flat
% axis square
% c = colorbar;
% c.Label.String = 'm s^{-1}';
% c.Label.FontSize = 13;
% c.Label.VerticalAlignment = 'middle';
% title('(b) v_{original}', 'FontSize', 14)
% xlabel('x (m)', 'FontSize', 14)
% ylabel('y (m)', 'FontSize', 14)
% 
% ha(3) = subplot(2,2,[3,4]);
% hs(3) = pcolor(xx, yy, vr_truth(:,:,1));
% caxis([-1 1] * va)
% colormap(ha(3), blib('rgmap2'))
% shading flat
% axis square
% c = colorbar;
% c.Label.String = 'm s^{-1}';
% c.Label.FontSize = 13;
% c.Label.VerticalAlignment = 'middle';
% title('(c) v_{truth}', 'FontSize', 14)
% xlabel('x (m)', 'FontSize', 14)
% ylabel('y (m)', 'FontSize', 14)
% 
% sgtitle('Aggregation-weighted correction', 'FontSize', 14)
% 
% if plot_save_flag
%     print([img_dir 'vcorr-aggreg'], '-dpng')
% end


% Plot error of uncorrected and corrected velocities
dv_uncorr = vr_old - vr_truth;
dv_corr = vr_new - vr_truth;
dv_corr_var = vr_new_var - vr_truth;
dv_corr_agg = vr_new_agg - vr_truth;

dv_uncorr_mean = mean(dv_uncorr, 'all');
dv_corr_mean = mean(dv_corr, 'all');
dv_corr_var_mean = mean(dv_corr_var, 'all');
dv_corr_agg_mean = mean(dv_corr_agg, 'all');

dv_uncorr_25 = prctile(dv_uncorr, 25, 'all');
dv_corr_25 = prctile(dv_corr, 25, 'all');
dv_corr_var_25 = prctile(dv_corr_var, 25, 'all');
dv_corr_agg_25 = prctile(dv_corr_agg, 25, 'all');

dv_uncorr_75 = prctile(dv_uncorr, 75, 'all');
dv_corr_75 = prctile(dv_corr, 75, 'all');
dv_corr_var_75 = prctile(dv_corr_var, 75, 'all');
dv_corr_agg_75 = prctile(dv_corr_agg, 75, 'all');


figure(5)
clf

ha = subplot(2,2,1);
hs = pcolor(xx, yy, dv_uncorr(:,:,1));
caxis([-1 1] * va)
colormap(ha, blib('rbmap'))
shading flat
c = colorbar;
c.Label.String = 'm s^{-1}';
c.Label.FontSize = 13;
c.Label.VerticalAlignment = 'middle';
title('(a) Pre-correction', 'FontSize', 14)
xlabel('x (m)', 'FontSize', 14)
ylabel('y (m)', 'FontSize', 14)

ha(2) = subplot(2,2,2);
hs(2) = pcolor(xx, yy, dv_corr(:,:,1));
caxis([-1 1] * va)
colormap(ha(2), blib('rbmap'))
shading flat
c = colorbar;
c.Label.String = 'm s^{-1}';
c.Label.FontSize = 13;
c.Label.VerticalAlignment = 'middle';
title('(b) Post-correction', 'FontSize', 14)
xlabel('x (m)', 'FontSize', 14)
ylabel('y (m)', 'FontSize', 14)

ha(3) = subplot(2,2,3);
histogram(dv_uncorr, -va:5:va)
hold on
plot(ones(1,2)*dv_uncorr_mean, [0 350], 'k', 'LineWidth', 1.2)
plot(ones(1,2)*dv_uncorr_25, [0 350], '--k', 'LineWidth', 1.2)
plot(ones(1,2)*dv_uncorr_75, [0 350], '--k', 'LineWidth', 1.2)
hold off
title('(c)', 'FontSize', 14)
xlabel('v_{bias} (m s^{-1})', 'FontSize', 14)
ylim([0 350])
text(-100, 330, ['\mu = ' num2str(double(dv_uncorr_mean),'%.2f')], 'FontSize', 13)
text(-100, 300, ['Q_1 = ' num2str(double(dv_uncorr_25),'%.2f')], 'FontSize', 13)
text(-100, 270, ['Q_3 = ' num2str(double(dv_uncorr_75),'%.2f')], 'FontSize', 13)

ha(4) = subplot(2,2,4);
histogram(dv_corr, -va:5:va)
hold on
plot(ones(1,2)*dv_corr_mean, [0 350], 'k', 'LineWidth', 1.2)
plot(ones(1,2)*dv_corr_25, [0 350], '--k', 'LineWidth', 1.2)
plot(ones(1,2)*dv_corr_75, [0 350], '--k', 'LineWidth', 1.2)
hold off
title('(d)', 'FontSize', 14)
xlabel('v_{bias} (m s^{-1})', 'FontSize', 14)
ylim([0 350])
text(-100, 330, ['\mu = ' num2str(double(dv_corr_mean),'%.2f')], 'FontSize', 13)
text(-100, 300, ['Q_1 = ' num2str(double(dv_corr_25),'%.2f')], 'FontSize', 13)
text(-100, 270, ['Q_3 = ' num2str(double(dv_corr_75),'%.2f')], 'FontSize', 13)

set(gcf, 'Units', 'inches', 'Position', [10 10 12 8])
axes('Unit', 'Normalized', 'Position', [0.5 0.48 0.01 0.01])
title(['Debris type ' dtype ', n=' dnum ', tilt=' elev '^o'], 'FontSize', 14);
axis off

%sgtitle('v_{bias} - Debris-removed correction','FontSize',16);

if plot_save_flag
    print([img_dir 'vbias-DCA-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
end



figure(6)
clf

ha = subplot(2,2,1);
hs = pcolor(xx, yy, dv_uncorr(:,:,1));
caxis([-1 1] * va)
colormap(ha, blib('rbmap'))
shading flat
c = colorbar;
c.Label.String = 'm s^{-1}';
c.Label.FontSize = 13;
c.Label.VerticalAlignment = 'middle';
title('(a) Pre-correction', 'FontSize', 14)
xlabel('x (m)', 'FontSize', 14)
ylabel('y (m)', 'FontSize', 14)

ha(2) = subplot(2,2,2);
hs(2) = pcolor(xx, yy, dv_corr_var(:,:,1));
caxis([-1 1] * va)
colormap(ha(2), blib('rbmap'))
shading flat
c = colorbar;
c.Label.String = 'm s^{-1}';
c.Label.FontSize = 13;
c.Label.VerticalAlignment = 'middle';
title('(b) Post-correction', 'FontSize', 14)
xlabel('x (m)', 'FontSize', 14)
ylabel('y (m)', 'FontSize', 14)

ha(3) = subplot(2,2,3);
histogram(dv_uncorr, -va:5:va)
hold on
plot(ones(1,2)*dv_uncorr_mean, [0 350], 'k', 'LineWidth', 1.2)
plot(ones(1,2)*dv_uncorr_25, [0 350], '--k', 'LineWidth', 1.2)
plot(ones(1,2)*dv_uncorr_75, [0 350], '--k', 'LineWidth', 1.2)
hold off
title('(c)', 'FontSize', 14)
xlabel('v_{bias} (m s^{-1})', 'FontSize', 14)
ylim([0 350])
text(-100, 330, ['\mu = ' num2str(double(dv_uncorr_mean),'%.2f')], 'FontSize', 13)
text(-100, 300, ['Q_1 = ' num2str(double(dv_uncorr_25),'%.2f')], 'FontSize', 13)
text(-100, 270, ['Q_3 = ' num2str(double(dv_uncorr_75),'%.2f')], 'FontSize', 13)

ha(4) = subplot(2,2,4);
histogram(dv_corr_var, -va:5:va)
hold on
plot(ones(1,2)*dv_corr_var_mean, [0 350], 'k', 'LineWidth', 1.2)
plot(ones(1,2)*dv_corr_var_25, [0 350], '--k', 'LineWidth', 1.2)
plot(ones(1,2)*dv_corr_var_75, [0 350], '--k', 'LineWidth', 1.2)
hold off
title('(d)', 'FontSize', 14)
xlabel('v_{bias} (m s^{-1})', 'FontSize', 14)
ylim([0 350])
text(-100, 330, ['\mu = ' num2str(double(dv_corr_var_mean),'%.2f')], 'FontSize', 13)
text(-100, 300, ['Q_1 = ' num2str(double(dv_corr_var_25),'%.2f')], 'FontSize', 13)
text(-100, 270, ['Q_3 = ' num2str(double(dv_corr_var_75),'%.2f')], 'FontSize', 13)

set(gcf, 'Units', 'inches', 'Position', [10 10 12 8])
axes('Unit', 'Normalized', 'Position', [0.5 0.48 0.01 0.01])
title(['Debris type ' dtype ', n=' dnum ', tilt=' elev '^o'], 'FontSize', 14);
axis off

%sgtitle('v_{bias} - Variance-weighted correction','FontSize',16);

if plot_save_flag
    print([img_dir 'vbias-var-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
end



figure(7)
clf

ha = subplot(2,2,1);
hs = pcolor(xx, yy, dv_uncorr(:,:,1));
caxis([-1 1] * va)
colormap(ha, blib('rbmap'))
shading flat
c = colorbar;
c.Label.String = 'm s^{-1}';
c.Label.FontSize = 13;
c.Label.VerticalAlignment = 'middle';
title('(a) Pre-correction', 'FontSize', 14)
xlabel('x (m)', 'FontSize', 14)
ylabel('y (m)', 'FontSize', 14)

ha(2) = subplot(2,2,2);
hs(2) = pcolor(xx, yy, dv_corr_agg(:,:,1));
caxis([-1 1] * va)
colormap(ha(2), blib('rbmap'))
shading flat
c = colorbar;
c.Label.String = 'm s^{-1}';
c.Label.FontSize = 13;
c.Label.VerticalAlignment = 'middle';
title('(b) Post-correction', 'FontSize', 14)
xlabel('x (m)', 'FontSize', 14)
ylabel('y (m)', 'FontSize', 14)

ha(3) = subplot(2,2,3);
histogram(dv_uncorr, -va:5:va)
hold on
plot(ones(1,2)*dv_uncorr_mean, [0 350], 'k', 'LineWidth', 1.2)
plot(ones(1,2)*dv_uncorr_25, [0 350], '--k', 'LineWidth', 1.2)
plot(ones(1,2)*dv_uncorr_75, [0 350], '--k', 'LineWidth', 1.2)
hold off
title('(c)', 'FontSize', 14)
xlabel('v_{bias} (m s^{-1})', 'FontSize', 14)
ylim([0 350])
text(-100, 330, ['\mu = ' num2str(double(dv_uncorr_mean),'%.2f')], 'FontSize', 13)
text(-100, 300, ['Q_1 = ' num2str(double(dv_uncorr_25),'%.2f')], 'FontSize', 13)
text(-100, 270, ['Q_3 = ' num2str(double(dv_uncorr_75),'%.2f')], 'FontSize', 13)

ha(4) = subplot(2,2,4);
histogram(dv_corr_agg, -va:5:va)
hold on
plot(ones(1,2)*dv_corr_agg_mean, [0 350], 'k', 'LineWidth', 1.2)
plot(ones(1,2)*dv_corr_agg_25, [0 350], '--k', 'LineWidth', 1.2)
plot(ones(1,2)*dv_corr_agg_75, [0 350], '--k', 'LineWidth', 1.2)
hold off
title('(d)', 'FontSize', 14)
xlabel('v_{bias} (m s^{-1})', 'FontSize', 14)
ylim([0 350])
text(-100, 330, ['\mu = ' num2str(double(dv_corr_agg_mean),'%.2f')], 'FontSize', 13)
text(-100, 300, ['Q_1 = ' num2str(double(dv_corr_agg_25),'%.2f')], 'FontSize', 13)
text(-100, 270, ['Q_3 = ' num2str(double(dv_corr_agg_75),'%.2f')], 'FontSize', 13)

set(gcf, 'Units', 'inches', 'Position', [10 10 12 8])
axes('Unit', 'Normalized', 'Position', [0.5 0.48 0.01 0.01])
title(['Debris type ' dtype ', n=' dnum ', tilt=' elev '^o'], 'FontSize', 14);
axis off

%sgtitle('v_{bias} - Aggregation-weighted correction','FontSize',16);

if plot_save_flag
    print([img_dir 'vbias-agg-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
end




% figure(8)
% clf
% 
% ha = subplot(3,2,1);
% hs = pcolor(xx, yy, dv_corr(:,:,1));
% caxis([-1 1] * va)
% colormap(ha, blib('rbmap'))
% shading flat
% c = colorbar;
% c.Label.String = 'm s^{-1}';
% c.Label.FontSize = 13;
% c.Label.VerticalAlignment = 'middle';
% %title('(a) Debris-removed correction', 'FontSize', 14)
% title('(a) v_{bias}', 'FontSize', 14)
% xlabel('x (m)', 'FontSize', 14)
% ylabel('y (m)', 'FontSize', 14)
% 
% subplot(3,2,2)
% histogram(dv_corr, -va:5:va)
% hold on
% plot(ones(1,2)*dv_corr_mean, [0 350], 'k', 'LineWidth', 1.2)
% plot(ones(1,2)*dv_corr_25, [0 350], '--k', 'LineWidth', 1.2)
% plot(ones(1,2)*dv_corr_75, [0 350], '--k', 'LineWidth', 1.2)
% hold off
% title('(b)', 'FontSize', 14)
% xlabel('v_{bias} (m s^{-1})', 'FontSize', 14)
% ylim([0 350])
% text(-100, 330, ['\mu = ' num2str(double(dv_corr_mean),'%.2f')], 'FontSize', 13)
% text(-100, 300, ['Q_1 = ' num2str(double(dv_corr_25),'%.2f')], 'FontSize', 13)
% text(-100, 270, ['Q_3 = ' num2str(double(dv_corr_75),'%.2f')], 'FontSize', 13)
% 
% ha(2) = subplot(3,2,3);
% hs(2) = pcolor(xx, yy, dv_corr_var(:,:,1));
% caxis([-1 1] * va)
% colormap(ha(2), blib('rbmap'))
% shading flat
% c = colorbar;
% c.Label.String = 'm s^{-1}';
% c.Label.FontSize = 13;
% c.Label.VerticalAlignment = 'middle';
% %title('(c) Variance-weighted correction', 'FontSize', 14)
% title('(c) v_{bias}', 'FontSize', 14)
% xlabel('x (m)', 'FontSize', 14)
% ylabel('y (m)', 'FontSize', 14)
% 
% subplot(3,2,4)
% histogram(dv_corr_var, -va:5:va)
% hold on
% plot(ones(1,2)*dv_corr_var_mean, [0 350], 'k', 'LineWidth', 1.2)
% plot(ones(1,2)*dv_corr_var_25, [0 350], '--k', 'LineWidth', 1.2)
% plot(ones(1,2)*dv_corr_var_75, [0 350], '--k', 'LineWidth', 1.2)
% hold off
% title('(d)', 'FontSize', 14)
% xlabel('v_{bias} (m s^{-1})', 'FontSize', 14)
% ylim([0 350])
% text(-100, 330, ['\mu = ' num2str(double(dv_corr_var_mean),'%.2f')], 'FontSize', 13)
% text(-100, 300, ['Q_1 = ' num2str(double(dv_corr_var_25),'%.2f')], 'FontSize', 13)
% text(-100, 270, ['Q_3 = ' num2str(double(dv_corr_var_75),'%.2f')], 'FontSize', 13)
% 
% ha(3) = subplot(3,2,5);
% hs(3) = pcolor(xx, yy, dv_corr_agg(:,:,1));
% caxis([-1 1] * va)
% colormap(ha(3), blib('rbmap'))
% shading flat
% c = colorbar;
% c.Label.String = 'm s^{-1}';
% c.Label.FontSize = 13;
% c.Label.VerticalAlignment = 'middle';
% %title('(e) Aggregation-weighted correction', 'FontSize', 14)
% title('(e) v_{bias}', 'FontSize', 14)
% xlabel('x (m)', 'FontSize', 14)
% ylabel('y (m)', 'FontSize', 14)
% 
% subplot(3,2,6)
% histogram(dv_corr_agg, -va:5:va)
% hold on
% plot(ones(1,2)*dv_corr_agg_mean, [0 350], 'k', 'LineWidth', 1.2)
% plot(ones(1,2)*dv_corr_agg_25, [0 350], '--k', 'LineWidth', 1.2)
% plot(ones(1,2)*dv_corr_agg_75, [0 350], '--k', 'LineWidth', 1.2)
% hold off
% title('(f)', 'FontSize', 14)
% xlabel('v_{bias} (m s^{-1})', 'FontSize', 14)
% ylim([0 350])
% text(-100, 330, ['\mu = ' num2str(double(dv_corr_agg_mean),'%.2f')], 'FontSize', 13)
% text(-100, 300, ['Q_1 = ' num2str(double(dv_corr_agg_25),'%.2f')], 'FontSize', 13)
% text(-100, 270, ['Q_3 = ' num2str(double(dv_corr_agg_75),'%.2f')], 'FontSize', 13)
% 
% set(gcf, 'Units', 'inches', 'Position', [10 10 12 12])
% axes('Unit', 'Normalized', 'Position', [0.5 0.03 0.01 0.01])
% title(['Debris type ' dtype ', n=' dnum ', tilt=' elev '^o'], 'FontSize', 14);
% axis off
% 
% % FUCKING HALLELUJAH
% annotation('textbox', [0.35 0.91 0.35 0.05], 'String', 'DCA-based correction', 'FontSize', 18, 'EdgeColor', 'none')
% annotation('textbox', [0.35 0.61 0.35 0.05], 'String', 'Variance-based correction', 'FontSize', 18, 'EdgeColor', 'none')
% annotation('textbox', [0.35 0.31 0.35 0.05], 'String', 'Aggregation-based correction', 'FontSize', 18, 'EdgeColor', 'none')
% 
% if plot_save_flag
%     print([img_dir 'vcorr-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
% end

figure(8)
clf

ha = subplot(4,2,1);
hs = pcolor(xx, yy, dv_uncorr(:,:,1));
caxis([-1 1] * va)
colormap(ha, blib('rbmap'))
shading flat
c = colorbar;
c.Label.String = 'm s^{-1}';
c.Label.FontSize = 13;
c.Label.VerticalAlignment = 'middle';
title('(a) v_{bias}', 'FontSize', 16)
%xlabel('x (m)', 'FontSize', 14)
ylabel('y (m)', 'FontSize', 14)

subplot(4,2,2);
histogram(dv_uncorr, -va:5:va)
hold on
plot(ones(1,2)*dv_uncorr_mean, [0 350], 'k', 'LineWidth', 1.2)
plot(ones(1,2)*dv_uncorr_25, [0 350], '--k', 'LineWidth', 1.2)
plot(ones(1,2)*dv_uncorr_75, [0 350], '--k', 'LineWidth', 1.2)
hold off
title('(b) v_{bias} distribution', 'FontSize', 16)
%xlabel('v_{bias} (m s^{-1})', 'FontSize', 14)
ylim([0 350])
text(-100, 330, ['\mu = ' num2str(double(dv_uncorr_mean),'%.2f')], 'FontSize', 13)
text(-100, 290, ['Q_1 = ' num2str(double(dv_uncorr_25),'%.2f')], 'FontSize', 13)
text(-100, 250, ['Q_3 = ' num2str(double(dv_uncorr_75),'%.2f')], 'FontSize', 13)

ha(2) = subplot(4,2,3);
hs(2) = pcolor(xx, yy, dv_corr(:,:,1));
caxis([-1 1] * va)
colormap(ha(2), blib('rbmap'))
shading flat
c = colorbar;
c.Label.String = 'm s^{-1}';
c.Label.FontSize = 13;
c.Label.VerticalAlignment = 'middle';
%title('(c) v_{bias}', 'FontSize', 14)
title('(c)', 'FontSize', 16)
%xlabel('x (m)', 'FontSize', 14)
ylabel('y (m)', 'FontSize', 14)

subplot(4,2,4)
histogram(dv_corr, -va:5:va)
hold on
plot(ones(1,2)*dv_corr_mean, [0 350], 'k', 'LineWidth', 1.2)
plot(ones(1,2)*dv_corr_25, [0 350], '--k', 'LineWidth', 1.2)
plot(ones(1,2)*dv_corr_75, [0 350], '--k', 'LineWidth', 1.2)
hold off
%title('(d) v_{bias} distribution', 'FontSize', 14)
title('(d)', 'FontSize', 16)
%xlabel('v_{bias} (m s^{-1})', 'FontSize', 14)
ylim([0 350])
text(-100, 330, ['\mu = ' num2str(double(dv_corr_mean),'%.2f')], 'FontSize', 13)
text(-100, 290, ['Q_1 = ' num2str(double(dv_corr_25),'%.2f')], 'FontSize', 13)
text(-100, 250, ['Q_3 = ' num2str(double(dv_corr_75),'%.2f')], 'FontSize', 13)

ha(3) = subplot(4,2,5);
hs(3) = pcolor(xx, yy, dv_corr_var(:,:,1));
caxis([-1 1] * va)
colormap(ha(3), blib('rbmap'))
shading flat
c = colorbar;
c.Label.String = 'm s^{-1}';
c.Label.FontSize = 13;
c.Label.VerticalAlignment = 'middle';
%title('(e) v_{bias}', 'FontSize', 14)
title('(e)', 'FontSize', 16)
%xlabel('x (m)', 'FontSize', 14)
ylabel('y (m)', 'FontSize', 14)

subplot(4,2,6)
histogram(dv_corr_var, -va:5:va)
hold on
plot(ones(1,2)*dv_corr_var_mean, [0 350], 'k', 'LineWidth', 1.2)
plot(ones(1,2)*dv_corr_var_25, [0 350], '--k', 'LineWidth', 1.2)
plot(ones(1,2)*dv_corr_var_75, [0 350], '--k', 'LineWidth', 1.2)
hold off
%title('(f) v_{bias} distribution', 'FontSize', 14)
title('(f)', 'FontSize', 16)
%xlabel('v_{bias} (m s^{-1})', 'FontSize', 14)
ylim([0 350])
text(-100, 330, ['\mu = ' num2str(double(dv_corr_var_mean),'%.2f')], 'FontSize', 13)
text(-100, 290, ['Q_1 = ' num2str(double(dv_corr_var_25),'%.2f')], 'FontSize', 13)
text(-100, 250, ['Q_3 = ' num2str(double(dv_corr_var_75),'%.2f')], 'FontSize', 13)

ha(4) = subplot(4,2,7);
hs(4) = pcolor(xx, yy, dv_corr_agg(:,:,1));
caxis([-1 1] * va)
colormap(ha(4), blib('rbmap'))
shading flat
c = colorbar;
c.Label.String = 'm s^{-1}';
c.Label.FontSize = 13;
c.Label.VerticalAlignment = 'middle';
%title('(g) v_{bias}', 'FontSize', 14)
title('(g)', 'FontSize', 16)
xlabel('x (m)', 'FontSize', 14)
ylabel('y (m)', 'FontSize', 14)

subplot(4,2,8)
histogram(dv_corr_agg, -va:5:va)
hold on
plot(ones(1,2)*dv_corr_agg_mean, [0 350], 'k', 'LineWidth', 1.2)
plot(ones(1,2)*dv_corr_agg_25, [0 350], '--k', 'LineWidth', 1.2)
plot(ones(1,2)*dv_corr_agg_75, [0 350], '--k', 'LineWidth', 1.2)
hold off
%title('(h) v_{bias} distribution', 'FontSize', 14)
title('(h)', 'FontSize', 16)
xlabel('v_{bias} (m s^{-1})', 'FontSize', 14)
ylim([0 350])
text(-100, 330, ['\mu = ' num2str(double(dv_corr_agg_mean),'%.2f')], 'FontSize', 13)
text(-100, 290, ['Q_1 = ' num2str(double(dv_corr_agg_25),'%.2f')], 'FontSize', 13)
text(-100, 250, ['Q_3 = ' num2str(double(dv_corr_agg_75),'%.2f')], 'FontSize', 13)

set(gcf, 'Units', 'inches', 'Position', [10 10 12 14])
axes('Unit', 'Normalized', 'Position', [0.5 0.03 0.01 0.01])
title(['Debris type ' dtype ', n=' dnum ', tilt=' elev '^o'], 'FontSize', 14);
axis off

% FUCKING HALLELUJAH
annotation('textbox', [0.43 0.91 0.35 0.05], 'String', 'No correction', 'FontSize', 16, 'EdgeColor', 'none')
annotation('textbox', [0.4 0.69 0.35 0.05], 'String', 'DCA-based correction', 'FontSize', 16, 'EdgeColor', 'none')
annotation('textbox', [0.38 0.47 0.35 0.05], 'String', 'Variance-based correction', 'FontSize', 16, 'EdgeColor', 'none')
annotation('textbox', [0.37 0.25 0.35 0.05], 'String', 'Aggregation-based correction', 'FontSize', 16, 'EdgeColor', 'none')

if plot_save_flag
    print([img_dir 'vcorr-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
end



% figure(9)
% 
% ha = subplot(2,1,1);
% hs = pcolor(xx, yy, dv_uncorr(:,:,1));
% caxis([-1 1] * va)
% colormap(ha, blib('rbmap'))
% shading flat
% c = colorbar;
% c.Label.String = 'm s^{-1}';
% c.Label.FontSize = 13;
% c.Label.VerticalAlignment = 'middle';
% title('(a) v_{bias} (uncorrected)', 'FontSize', 14)
% xlabel('x (m)', 'FontSize', 14)
% ylabel('y (m)', 'FontSize', 14)
% 
% 
% subplot(2,1,2)
% histogram(dv_uncorr, -va:5:va)
% hold on
% plot(ones(1,2)*dv_uncorr_mean, [0 350], 'k', 'LineWidth', 1.2)
% plot(ones(1,2)*dv_uncorr_25, [0 350], '--k', 'LineWidth', 1.2)
% plot(ones(1,2)*dv_uncorr_75, [0 350], '--k', 'LineWidth', 1.2)
% hold off
% title('(b)', 'FontSize', 14)
% xlabel('v_{bias} (m s^{-1})', 'FontSize', 14)
% ylim([0 350])
% text(-100, 330, ['\mu = ' num2str(double(dv_uncorr_mean),'%.2f')], 'FontSize', 13)
% text(-100, 300, ['Q_1 = ' num2str(double(dv_uncorr_25),'%.2f')], 'FontSize', 13)
% text(-100, 270, ['Q_3 = ' num2str(double(dv_uncorr_75),'%.2f')], 'FontSize', 13)


% figure(9)
% 
% ha = subplot(3,1,1);
% hs = pcolor(xx, yy, vr_new(:,:,1) - vr_old(:,:,1));
% caxis([-1 1] * va)
% colormap(ha, blib('rbmap'))
% shading flat
% c = colorbar;
% c.Label.String = 'm s^{-1}';
% c.Label.FontSize = 13;
% c.Label.VerticalAlignment = 'middle';
% title('(a) Debris-removed v correction', 'FontSize', 14)
% xlabel('x (m)', 'FontSize', 14)
% ylabel('y (m)', 'FontSize', 14)
% 
% ha(2) = subplot(3,1,2);
% hs(2) = pcolor(xx, yy, vr_new_var(:,:,1) - vr_old(:,:,1));
% caxis([-1 1] * va)
% colormap(ha(2), blib('rbmap'))
% shading flat
% c = colorbar;
% c.Label.String = 'm s^{-1}';
% c.Label.FontSize = 13;
% c.Label.VerticalAlignment = 'middle';
% title('(b) Variance-weighted v correction', 'FontSize', 14)
% xlabel('x (m)', 'FontSize', 14)
% ylabel('y (m)', 'FontSize', 14)
% 
% ha(3) = subplot(3,1,3);
% hs(3) = pcolor(xx, yy, vr_new_agg(:,:,1) - vr_old(:,:,1));
% caxis([-1 1] * va)
% colormap(ha(3), blib('rbmap'))
% shading flat
% c = colorbar;
% c.Label.String = 'm s^{-1}';
% c.Label.FontSize = 13;
% c.Label.VerticalAlignment = 'middle';
% title('(c) Aggregation-weighted v correction', 'FontSize', 14)
% xlabel('x (m)', 'FontSize', 14)
% ylabel('y (m)', 'FontSize', 14)



%%


rind = 15;
azind = 21;

figure(10)
clf

ax1 = subplot(3,1,1);
% DPSD
yyaxis left
plot(vvx(:), squeeze(szdr(:,rind,azind)), 'k', 'LineWidth', 1.1)
ylabel('sZ_{DR} (dB)', 'FontSize', 14)
xlim([-va, va])
ylim([-15, 15])
% PSD
yyaxis right
hs1 = semilogy(vvx(:), squeeze(abs(spsd.h(:,rind,azind))), 'LineWidth', 1);
hs1.Color = [0.6 0.6 0.6];
hs1.LineStyle = '-.';
ylabel('sS_H', 'FontSize', 14)
ax1.YAxis(1).Color = 'k';
ax1.YAxis(2).Color = [0.6 0.6 0.6];
grid on
title('(a) sZ_{DR}, sS_H', 'FontSize', 14)
xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
%legend('sZ_{DR}', 'sS_V', 'Location', 'northwest')
% text(18, 1, 'Debris', 'FontSize', 20)
% text(-30, 1, 'Rain', 'FontSize', 20)
% text(-48, 3.4, 'sZ_{DR}', 'FontSize', 14)
% text(-48, -3.5, 'sS_V', 'FontSize', 14, 'Color', [0.6 0.6 0.6])


ax2 = subplot(3,1,2);
% DPSD
yyaxis left
plot(vvx(:), squeeze(sphv(:,rind,azind)), 'k', 'LineWidth', 1.1)
ylabel('s\rho_{HV}', 'FontSize', 14)
xlim([-va, va])
ylim([0, 1])
% PSD
yyaxis right
hs2 = semilogy(vvx(:), squeeze(abs(spsd.h(:,rind,azind))), 'LineWidth', 1);
hs2.Color = [0.6 0.6 0.6];
hs2.LineStyle = '-.';
ylabel('sS_H', 'FontSize', 14)
ax2.YAxis(1).Color = 'k';
ax2.YAxis(2).Color = [0.6 0.6 0.6];
grid on
title('(b) s\rho_{HV}, sS_H', 'FontSize', 14)
xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
% text(18, 0.6, 'Debris', 'FontSize', 20)
% text(-30, 0.6, 'Rain', 'FontSize', 20)
% text(-48, 0.95, 's\rho_{HV}', 'FontSize', 14)
% text(-48, 0.15, 'sS_V', 'FontSize', 14, 'Color', [0.6 0.6 0.6])


ax3 = subplot(3,1,3);
plot(vvx(:), squeeze(obj_class(:,rind,azind)), 'k', 'LineWidth', 1.1)
grid on
title('(c) Scatterer classification', 'FontSize', 14)
xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
xlim([-va, va])
ylim([-0.5,1.5])
yticks([0 1])
yticklabels({'Rain', 'Debris'})
%text(-48, -1.8, 'sZ_{DR}', 'FontSize', 14)
%text(-47, -13.5, 'sS_H', 'FontSize', 14, 'Color', [0.6 0.6 0.6])

if plot_save_flag
    print([img_dir 'hcaclass-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
end


figure(11)
clf

subplot(3,1,1)
plot(vvx(:), squeeze(obj_class(:,rind,azind)), 'k', 'LineWidth', 1.1)
title('(a) Scatterer classification', 'FontSize', 14)
%xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
xlim([-va, va])
ylim([-0.5,1.5])
yticks([0 1])
yticklabels({'Rain', 'Debris'})

subplot(3,1,2)
semilogy(vvx(:), squeeze(abs(psd_old(:,rind,azind))), 'k', 'LineWidth', 1)
hold on
semilogy(ones(1,2)*vr_old(rind,azind), [1e0 1e12], '--k', 'LineWidth', 1.3)
semilogy(ones(1,2)*vr_truth(rind,azind), [1e0 1e12], ':k', 'LineWidth', 2)
hold off
%xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
%ylabel('sS_H', 'FontSize', 14)
xlim([-va, va])
ylim([1e0 1e12])
title('(b) Original sS_H(v) and v_r', 'FontSize', 14)
if vr_old(rind,azind) >= vr_truth(rind,azind)
    text_pos1 = double(vr_old(rind,azind)+3);
    text_pos2 = double(vr_truth(rind,azind)-51);
else
    text_pos1 = double(vr_old(rind,azind)-51);
    text_pos2 = double(vr_truth(rind,azind)+3);
end
text(text_pos1, 1e2, ['v_{old} = ' num2str(round(vr_old(rind,azind),1)) ' m s^{-1}'], 'FontSize', 12)
text(text_pos2, 1e2, ['v_{truth} = ' num2str(round(vr_truth(rind,azind),1)) ' m s^{-1}'], 'FontSize', 12)

subplot(3,1,3)
semilogy(vvx(:), squeeze(abs(psd_filt(:,rind,azind))), 'k', 'LineWidth', 1.5)
hold on
semilogy(ones(1,2)*vr_new(rind,azind), [1e0 1e12], '--k', 'LineWidth', 1.3)
semilogy(ones(1,2)*vr_truth(rind,azind), [1e0 1e12], ':k', 'LineWidth', 2)
hold off
xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
%ylabel('sS_H', 'FontSize', 14)
xlim([-va, va])
ylim([1e0 1e12])
title('(c) Filtered sS_H(v) and corrected v_r', 'FontSize', 14)
if vr_new(rind,azind) >= vr_truth(rind,azind)
    text_pos1 = double(vr_new(rind,azind)+3);
    text_pos2 = double(vr_truth(rind,azind)-51);
else
    text_pos1 = double(vr_new(rind,azind)-51);
    text_pos2 = double(vr_truth(rind,azind)+3);
end
text(text_pos1, 1e1, ['v_{new} = ' num2str(round(vr_new(rind,azind),1)) ' m s^{-1}'], 'FontSize', 12)
text(text_pos2, 1e1, ['v_{truth} = ' num2str(round(vr_truth(rind,azind),1)) ' m s^{-1}'], 'FontSize', 12)

if plot_save_flag
    print([img_dir 'vnew-DCA-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
end



figure(12)
clf

subplot(3,1,1)
plot(vvx(:), squeeze(svar.phv(:,rind,azind)), 'k', 'LineWidth', 1.1)
title('(a) \sigma^2_{s\rho_{HV}}', 'FontSize', 14)
%xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
xlim([-va, va])
ylim([0, 0.3])


subplot(3,1,2)
semilogy(vvx(:), squeeze(abs(psd_old(:,rind,azind))), 'k', 'LineWidth', 1)
hold on
semilogy(ones(1,2)*vr_old(rind,azind), [1e0 1e12], '--k', 'LineWidth', 1.3)
semilogy(ones(1,2)*vr_truth(rind,azind), [1e0 1e12], ':k', 'LineWidth', 2)
hold off
%xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
%ylabel('sS_H', 'FontSize', 14)
xlim([-va, va])
ylim([1e0 1e12])
title('(b) Original sS_H(v) and v_r', 'FontSize', 14)
if vr_old(rind,azind) >= vr_truth(rind,azind)
    text_pos1 = double(vr_old(rind,azind)+3);
    text_pos2 = double(vr_truth(rind,azind)-51);
else
    text_pos1 = double(vr_old(rind,azind)-51);
    text_pos2 = double(vr_truth(rind,azind)+3);
end
text(text_pos1, 1e1, ['v_{old} = ' num2str(round(vr_old(rind,azind),1)) ' m s^{-1}'], 'FontSize', 12)
text(text_pos2, 1e1, ['v_{truth} = ' num2str(round(vr_truth(rind,azind),1)) ' m s^{-1}'], 'FontSize', 12)


subplot(3,1,3)
semilogy(vvx(:), squeeze(abs(psd_new_var(:,rind,azind))), 'k', 'LineWidth', 1)
hold on
semilogy(ones(1,2)*vr_new_var(rind,azind), [1e0 1e12], '--k', 'LineWidth', 1.3)
semilogy(ones(1,2)*vr_truth(rind,azind), [1e0 1e12], ':k', 'LineWidth', 2)
hold off
xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
%ylabel('sS_H', 'FontSize', 14)
xlim([-va, va])
ylim([1e0 1e12])
title('(c) Filtered sS_H(v) and corrected v_r', 'FontSize', 14)
if vr_new_var(rind,azind) >= vr_truth(rind,azind)
    text_pos1 = double(vr_new_var(rind,azind)+3);
    text_pos2 = double(vr_truth(rind,azind)-51);
else
    text_pos1 = double(vr_new_var(rind,azind)-51);
    text_pos2 = double(vr_truth(rind,azind)+3);
end
text(text_pos1, 1e1, ['v_{new} = ' num2str(round(vr_new_var(rind,azind),1)) ' m s^{-1}'], 'FontSize', 12)
text(text_pos2, 1e1, ['v_{truth} = ' num2str(round(vr_truth(rind,azind),1)) ' m s^{-1}'], 'FontSize', 12)

if plot_save_flag
    print([img_dir 'vnew-var-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
end



figure(13)
clf

subplot(3,1,1)
plot(vvx(:), squeeze(agg.rain(:,rind,azind)), 'k', 'LineWidth', 1.1)
title('(a) A_{rain}', 'FontSize', 14)
%xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
xlim([-va, va])
ylim([0, 1])


subplot(3,1,2)
semilogy(vvx(:), squeeze(abs(psd_old(:,rind,azind))), 'k', 'LineWidth', 1)
hold on
semilogy(ones(1,2)*vr_old(rind,azind), [1e0 1e12], '--k', 'LineWidth', 1.3)
semilogy(ones(1,2)*vr_truth(rind,azind), [1e0 1e12], ':k', 'LineWidth', 2)
hold off
%xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
%ylabel('sS_H', 'FontSize', 14)
xlim([-va, va])
ylim([1e0 1e12])
title('(b) Original sS_H(v) and v_r', 'FontSize', 14)
if vr_old(rind,azind) >= vr_truth(rind,azind)
    text_pos1 = double(vr_old(rind,azind)+3);
    text_pos2 = double(vr_truth(rind,azind)-51);
else
    text_pos1 = double(vr_old(rind,azind)-51);
    text_pos2 = double(vr_truth(rind,azind)+3);
end
text(text_pos1, 1e1, ['v_{old} = ' num2str(round(vr_old(rind,azind),1)) ' m s^{-1}'], 'FontSize', 12)
text(text_pos2, 1e1, ['v_{truth} = ' num2str(round(vr_truth(rind,azind),1)) ' m s^{-1}'], 'FontSize', 12)


subplot(3,1,3)
semilogy(vvx(:), squeeze(abs(psd_new_agg(:,rind,azind))), 'k', 'LineWidth', 1)
hold on
semilogy(ones(1,2)*vr_new_agg(rind,azind), [1e0 1e12], '--k', 'LineWidth', 1.3)
semilogy(ones(1,2)*vr_truth(rind,azind), [1e0 1e12], ':k', 'LineWidth', 2)
hold off
xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
%ylabel('sS_H', 'FontSize', 14)
xlim([-va, va])
ylim([1e0 1e12])
title('(c) Filtered sS_H(v) and corrected v_r', 'FontSize', 14)
if vr_new_agg(rind,azind) >= vr_truth(rind,azind)
    text_pos1 = double(vr_new_agg(rind,azind)+3);
    text_pos2 = double(vr_truth(rind,azind)-51);
else
    text_pos1 = double(vr_new_agg(rind,azind)-51);
    text_pos2 = double(vr_truth(rind,azind)+3);
end
text(text_pos1, 1e1, ['v_{new} = ' num2str(round(vr_new_agg(rind,azind),1)) ' m s^{-1}'], 'FontSize', 12)
text(text_pos2, 1e1, ['v_{truth} = ' num2str(round(vr_truth(rind,azind),1)) ' m s^{-1}'], 'FontSize', 12)

if plot_save_flag
    print([img_dir 'vnew-agg-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
end



figure(14)
clf

ax = subplot(3,1,1);
% DPSD
yyaxis left
plot(vvx(:), squeeze(szdr(:,rind,azind)), 'k', 'LineWidth', 1.1)
ylabel('sZ_{DR} (dB)', 'FontSize', 14)
xlim([-va, va])
ylim([-15, 15])
% PSD
yyaxis right
hs = semilogy(vvx(:), squeeze(abs(spsd.h(:,rind,azind))), 'LineWidth', 1);
hs.Color = [0.6 0.6 0.6];
hs.LineStyle = '-.';
ylabel('sS_H', 'FontSize', 14)
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = [0.6 0.6 0.6];
grid on
title('(a) sZ_{DR}, sS_H', 'FontSize', 14)
%xlabel('{\it v} (m s^{-1})', 'FontSize', 14)

subplot(3,1,2)
plot(vvx(:), squeeze(memf.zdr_rain(:,rind,azind)), 'k', 'LineWidth', 1)
xlim([-va va])
ylim([-0.2 1.2])
yticks(0:0.2:1)
grid on
title('(b) sZ_{DR} rain membership')
%xlabel('{\it v} (m s^{-1})', 'FontSize', 14)

subplot(3,1,3)
plot(vvx(:), squeeze(memf.zdr_debr(:,rind,azind)), 'k', 'LineWidth', 1)
xlim([-va va])
ylim([-0.2 1.2])
yticks(0:0.2:1)
grid on
title('(c) sZ_{DR} debris membership')
xlabel('{\it v} (m s^{-1})', 'FontSize', 14)

if plot_save_flag
    print([img_dir 'memf-szdr-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
end



figure(15)
clf

ax = subplot(3,1,1);
% DPSD
yyaxis left
plot(vvx(:), squeeze(sphv(:,rind,azind)), 'k', 'LineWidth', 1.1)
ylabel('s\rho_{HV}', 'FontSize', 14)
xlim([-va, va])
ylim([0, 1])
% PSD
yyaxis right
hs = semilogy(vvx(:), squeeze(abs(spsd.h(:,rind,azind))), 'LineWidth', 1);
hs.Color = [0.6 0.6 0.6];
hs.LineStyle = '-.';
ylabel('sS_H', 'FontSize', 14)
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = [0.6 0.6 0.6];
grid on
title('(a) s\rho_{HV}, sS_H', 'FontSize', 14)
%xlabel('{\it v} (m s^{-1})', 'FontSize', 14)

subplot(3,1,2)
plot(vvx(:), squeeze(memf.phv_rain(:,rind,azind)), 'k', 'LineWidth', 1)
xlim([-va va])
ylim([-0.2 1.2])
yticks(0:0.2:1)
grid on
title('(b) s\rho_{HV} rain membership')
%xlabel('{\it v} (m s^{-1})', 'FontSize', 14)

subplot(3,1,3)
plot(vvx(:), squeeze(memf.phv_debr(:,rind,azind)), 'k', 'LineWidth', 1)
xlim([-va va])
ylim([-0.2 1.2])
yticks(0:0.2:1)
grid on
title('(c) s\rho_{HV} debris membership')
xlabel('{\it v} (m s^{-1})', 'FontSize', 14)

if plot_save_flag
    print([img_dir 'memf-sphv-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
end



figure(16)
clf

ax = subplot(3,1,1);
% DPSD
yyaxis left
plot(vvx(:), squeeze(svar.phv(:,rind,azind)), 'k', 'LineWidth', 1.1)
ylabel('\sigma^2_{s\rho_{HV}}', 'FontSize', 14)
xlim([-va, va])
ylim([0, 0.3])
% PSD
yyaxis right
hs = semilogy(vvx(:), squeeze(abs(spsd.h(:,rind,azind))), 'LineWidth', 1);
hs.Color = [0.6 0.6 0.6];
hs.LineStyle = '-.';
ylabel('sS_H', 'FontSize', 14)
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = [0.6 0.6 0.6];
grid on
title('(a) \sigma^2_{s\rho_{HV}}, sS_H', 'FontSize', 14)
%xlabel('{\it v} (m s^{-1})', 'FontSize', 14)

subplot(3,1,2)
plot(vvx(:), squeeze(memf.pvar_rain(:,rind,azind)), 'k', 'LineWidth', 1)
xlim([-va va])
ylim([-0.2 1.2])
yticks(0:0.2:1)
grid on
title('(b) \sigma^2_{s\rho_{HV}} rain membership')
%xlabel('{\it v} (m s^{-1})', 'FontSize', 14)

subplot(3,1,3)
plot(vvx(:), squeeze(memf.pvar_debr(:,rind,azind)), 'k', 'LineWidth', 1)
xlim([-va va])
ylim([-0.2 1.2])
yticks(0:0.2:1)
grid on
title('(c) \sigma^2_{s\rho_{HV}} debris membership')
xlabel('{\it v} (m s^{-1})', 'FontSize', 14)

if plot_save_flag
    print([img_dir 'memf-pvar-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
end



figure(17)

subplot(3,1,1)
semilogy(ones(1,2)*vr_truth(rind,azind), [1e0 1e15], ':k', 'LineWidth', 2)
hold on
semilogy(ones(1,2)*vr_new(rind,azind), [1e0 1e15], ':r', 'LineWidth', 2)
semilogy(ones(1,2)*vr_old(rind,azind), [1e0 1e15], ':b', 'LineWidth', 2)
semilogy(vvx(:), squeeze(abs(psd_filt(:,rind,azind))), 'k', 'LineWidth', 2)
hold off
%xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
ylabel('sS_H', 'FontSize', 14)
xlim([-va, va])
ylim([1e0 1e15])
title('(a) DCA-based correction', 'FontSize', 14)
legend(['v_{truth} = ' num2str(round(vr_truth(rind,azind),1))],...
    ['v_{new} = ' num2str(round(vr_new(rind,azind),1))],...
    ['v_{old} = ' num2str(round(vr_old(rind,azind),1))], 'Location', 'northwest')


subplot(3,1,2)
semilogy(ones(1,2)*vr_truth(rind,azind), [1e0 1e15], ':k', 'LineWidth', 2)
hold on
semilogy(ones(1,2)*vr_new_var(rind,azind), [1e0 1e15], ':r', 'LineWidth', 2)
semilogy(ones(1,2)*vr_old(rind,azind), [1e0 1e15], ':b', 'LineWidth', 2)
semilogy(vvx(:), squeeze(abs(psd_new_var(:,rind,azind))), 'k', 'LineWidth', 1)
hold off
%xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
ylabel('sS_H', 'FontSize', 14)
xlim([-va, va])
ylim([1e0 1e15])
title('(b) Variance-based correction', 'FontSize', 14)
legend(['v_{truth} = ' num2str(round(vr_truth(rind,azind),1))],...
    ['v_{new} = ' num2str(round(vr_new_var(rind,azind),1))],...
    ['v_{old} = ' num2str(round(vr_old(rind,azind),1))], 'Location', 'northwest')

subplot(3,1,3)
semilogy(ones(1,2)*vr_truth(rind,azind), [1e0 1e15], ':k', 'LineWidth', 2)
hold on
semilogy(ones(1,2)*vr_new_agg(rind,azind), [1e0 1e15], ':r', 'LineWidth', 2)
semilogy(ones(1,2)*vr_old(rind,azind), [1e0 1e15], ':b', 'LineWidth', 2)
semilogy(vvx(:), squeeze(abs(psd_new_agg(:,rind,azind))), 'k', 'LineWidth', 1)
hold off
xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
ylabel('sS_H', 'FontSize', 14)
xlim([-va, va])
ylim([1e0 1e15])
title('(c) Aggregation-based correction', 'FontSize', 14)
legend(['v_{truth} = ' num2str(round(vr_truth(rind,azind),1))],...
    ['v_{new} = ' num2str(round(vr_new_agg(rind,azind),1))],...
    ['v_{old} = ' num2str(round(vr_old(rind,azind),1))], 'Location', 'northwest')

set(gcf, 'Units', 'inches', 'Position', [10 10 9 10])

if plot_save_flag
    print([img_dir 'vnew-all-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
end

%%

figure(18)

subplot(2,2,1)
%scatter(reshape(params.phv,[],1), reshape(dv_uncorr,[],1), '.')
scatter(reshape(params.phv,[],1), abs(reshape(dv_uncorr,[],1)), '.')
title('(a) No correction', 'FontSize', 14)
xlabel('\rho_{HV}', 'FontSize', 14)
ylabel('v_{bias} (m s^{-1})', 'FontSize', 14)
xlim([0 1])

subplot(2,2,2)
%scatter(reshape(params.phv,[],1), reshape(dv_corr,[],1), '.')
scatter(reshape(params.phv,[],1), abs(reshape(dv_corr,[],1)), '.')
title('(b) DCA correction', 'FontSize', 14)
xlabel('\rho_{HV}', 'FontSize', 14)
ylabel('v_{bias} (m s^{-1})', 'FontSize', 14)

subplot(2,2,3)
%scatter(reshape(params.phv,[],1), reshape(dv_corr_var,[],1), '.')
scatter(reshape(params.phv,[],1), abs(reshape(dv_corr_var,[],1)), '.')
title('(c) Variance correction', 'FontSize', 14)
xlabel('\rho_{HV}', 'FontSize', 14)
ylabel('v_{bias} (m s^{-1})', 'FontSize', 14)

subplot(2,2,4)
%scatter(reshape(params.phv,[],1), reshape(dv_corr_agg,[],1), '.')
scatter(reshape(params.phv,[],1), abs(reshape(dv_corr_agg,[],1)), '.')
title('(d) Aggregation correction', 'FontSize', 14)
xlabel('\rho_{HV}', 'FontSize', 14)
ylabel('v_{bias} (m s^{-1})', 'FontSize', 14)

axes('Unit', 'Normalized', 'Position', [0.5 0.95 0.01 0.01])
title(['Debris type ' dtype ', n=' dnum ', tilt=' elev '^o'], 'FontSize', 14);
axis off



figure(19)

subplot(2,2,1)
%scatter(reshape(vr_truth,[],1), reshape(dv_uncorr,[],1), '.')
scatter(reshape(vr_truth,[],1), abs(reshape(dv_uncorr,[],1)), '.')
title('(a) No correction', 'FontSize', 14)
xlabel('v_{truth} (m s^{-1})', 'FontSize', 14)
ylabel('v_{bias} (m s^{-1})', 'FontSize', 14)

subplot(2,2,2)
%scatter(reshape(vr_truth,[],1), reshape(dv_corr,[],1), '.')
scatter(reshape(vr_truth,[],1), abs(reshape(dv_corr,[],1)), '.')
title('(b) DCA correction', 'FontSize', 14)
xlabel('v_{truth} (m s^{-1})', 'FontSize', 14)
ylabel('v_{bias} (m s^{-1})', 'FontSize', 14)

subplot(2,2,3)
%scatter(reshape(vr_truth,[],1), reshape(dv_corr_var,[],1), '.')
scatter(reshape(vr_truth,[],1), abs(reshape(dv_corr_var,[],1)), '.')
title('(c) Variance correction', 'FontSize', 14)
xlabel('v_{truth} (m s^{-1})', 'FontSize', 14)
ylabel('v_{bias} (m s^{-1})', 'FontSize', 14)

subplot(2,2,4)
%scatter(reshape(vr_truth,[],1), reshape(dv_corr_agg,[],1), '.')
scatter(reshape(vr_truth,[],1), abs(reshape(dv_corr_agg,[],1)), '.')
title('(d) Aggregation correction', 'FontSize', 14)
xlabel('v_{truth} (m s^{-1})', 'FontSize', 14)
ylabel('v_{bias} (m s^{-1})', 'FontSize', 14)

axes('Unit', 'Normalized', 'Position', [0.5 0.95 0.01 0.01])
title(['Debris type ' dtype ', n=' dnum ', tilt=' elev '^o'], 'FontSize', 14);
axis off








%% Save variables

if var_save_flag
    save([sim_dir '/dpsd-' erase(fname,'sim-'), '.mat'], 'szdr', 'sphv', 'svar', 'spsd', 'iqh',...
        'iqv', 'vr_old', 'vr_new', 'vr_truth', 'obj_class', 'vvx', 'params');
end