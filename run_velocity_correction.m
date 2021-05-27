clear all
close all
save_flag = 0;

%%%%% Things to try???
% -stricter thresholds for debris to get more rain flags
%   -add Zh/bulk rhoHV threshold/dependency for filtering
%   -add SNR threshold for filtering
%       -take max value of PSD and correct points within 10/20/?? dB
% -look at more individual spectra --- ginput to click point on figure and
% retrieve x,y coords [x,y] = ginput(1)
%
% -test real data from David

% Load IQ data
base_dir = '~/Documents/sims/';
sim_dir = uigetdir(base_dir); % Location of IQ files
filename = blib('choosefile', sim_dir, '*.mat');
load(filename);

% Separate path name and file name for saving files later
inds = strfind(filename, '/');
fpath = filename(1 : inds(end));
fname = filename(inds(end)+1 : end-4);

% Put all the bulk variables from checkiq.m into one structure
params = struct('zdr', zdr, 'phv', rhohv,'xx', xx, 'yy', yy, 'zz', zz,...
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

% HCA debris classification
[debris_flags, agg] = hca_class(szdr.cf, sphv.cf, svar);

% debris_flags = zeros(M,size(iqh,2),size(iqh,3));
% agg_rain = zeros(M,size(iqh,2),size(iqh,3));
% for i = 1:size(iqh,2)
%     for j = 1:size(iqh,3)
%         [obj_class, agg] = hca_class(szdr.cf(:,i,j), sphv.cf(:,i,j), svar(:,i,j));
%         debris_flags(:,i,j) = obj_class;
%         agg_rain(:,i,j) = agg.rain;
%     end
% end


% Filter debris from PSD
psd_filt = psd_old;
obj_class(1) = 0; % can't remember why I did this tbh but I didn't want endpoints marked as debris?
obj_class(end) = 0;
psd_filt(obj_class==1) = 0; % set all debris coefficients to 0

% Velocity recalculation -- need to make this more efficient than a loop!
var_norm = svar ./ max(svar,1); % normalize variance
b1 = 4; % aggressiveness of filter
var_inv = (1 - var_norm).^b1; % invert normalized variance
psd_new2 = var_inv .* psd_old;

b2 = 4; % aggressiveness of filter
psd_new3 = (agg.rain.^b2) .* psd_old;

vr_new = zeros(size(iqh,2), size(iqh,3));
vr_new2 = vr_new;
vr_new3 = vr_new;
for i = 1:size(iqh,2)
    for j = 1:size(iqh,3)
        vr_temp = calc_velocity_spectral(psd_filt(:,i,j), M, dat.params.prt, dat.params.lambda);
        if length(vr_temp) > 1
            vr_new(i,j) = vr_old(i,j);
        else
            vr_new(i,j) = vr_temp;
        end
        
        vr_new2(i,j) = calc_velocity_spectral(psd_new2(:,i,j), M, dat.params.prt, dat.params.lambda);
        vr_new3(i,j) = calc_velocity_spectral(psd_new3(:,i,j), M, dat.params.prt, dat.params.lambda);
    end
end

clear vr

% Load rain file for comparison
base_dir = '~/Documents/sims/';
sim_dir = uigetdir(base_dir); % Location of IQ files
filename = blib('choosefile', sim_dir, '*.mat');
load(filename);

vr_truth = vr;

%%
% Plot Doppler velocities
figure(1)

ha = subplot(3,1,1);
hs = pcolor(xx, yy, vr_old(:,:,1));
caxis([-1 1] * va)
colormap(ha, blib('rgmap2'))
shading flat
colorbar
title('Uncorrected Velocity (m/s)')

ha(2) = subplot(3,1,2);
hs(2) = pcolor(xx, yy, vr_new(:,:,1));
caxis([-1 1] * va)
colormap(ha(2), blib('rgmap2'))
shading flat
colorbar
title('Corrected Velocity (m/s)')

ha(3) = subplot(3,1,3);
hs(3) = pcolor(xx, yy, vr_truth(:,:,1));
caxis([-1 1] * va)
colormap(ha(3), blib('rgmap2'))
shading flat
colorbar
title('True Velocity (m/s)')


% Plot error of uncorrected and corrected velocities
dv_uncorr = vr_old - vr_truth;
dv_corr = vr_new - vr_truth;

figure(2)

ha = subplot(2,1,1);
hs = pcolor(xx, yy, dv_uncorr(:,:,1));
caxis([-1 1] * va)
colormap(ha, blib('rbmap'))
shading flat
colorbar
title('Error pre-correction')

ha(2) = subplot(2,1,2);
hs(2) = pcolor(xx, yy, dv_corr(:,:,1));
caxis([-1 1] * va)
colormap(ha(2), blib('rbmap'))
shading flat
colorbar
title('Error post-correction')


figure(3)

ha = subplot(2,2,1);
hs = pcolor(xx, yy, vr_new(:,:,1));
caxis([-1 1] * va)
colormap(ha, blib('rgmap2'))
shading flat
colorbar
title('Corrected velocity, debris removed')

ha(2) = subplot(2,2,2);
hs(2) = pcolor(xx, yy, vr_new2(:,:,1));
caxis([-1 1] * va)
colormap(ha(2), blib('rgmap2'))
shading flat
colorbar
title('Corrected velocity, variance-weighted')

ha(3) = subplot(2,2,3);
hs(3) = pcolor(xx, yy, vr_new3(:,:,1));
caxis([-1 1] * va)
colormap(ha(3), blib('rgmap2'))
shading flat
colorbar
title('Corrected velocity, aggregation-weighted')

ha(4) = subplot(2,2,4);
hs(4) = pcolor(xx, yy, vr_truth(:,:,1));
caxis([-1 1] * va)
colormap(ha(4), blib('rgmap2'))
shading flat
colorbar
title('True velocity')

%% Save variables

if save_flag
    dpsd_dir = [fpath 'DPSD'];
    if ~exist(dpsd_dir, 'dir')
        mkdir(dpsd_dir)
        addpath(genpath(dpsd_dir))
    end
    
    save([dpsd_dir '/dpsd-' fname '.mat'], 'szdr', 'sphv', 'svar', 'spsd', 'iqh',...
        'iqv', 'vr_old', 'vr_new', 'vr_truth', 'obj_class', 'vvx', 'params');
end