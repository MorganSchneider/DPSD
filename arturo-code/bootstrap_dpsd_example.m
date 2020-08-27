%% Load example KOUN 5/20/2013 data
load('ex_data.mat');

dr = 250;
va = 27.5;

t_az = data.az_set;
t_el = data.el_set;
t_iqh = data.H;
t_iqv = data.V;
time_string = data.time_set;
el_deg_file = median(t_el);

% KOUN Noise Values
N0.H = 10^(-67*0.1);
N0.V = 10^(-67*0.1);

% Select range gates to compute
gt_range = 3:210;
n_gt = numel(gt_range);
r_unique = gt_range * dr * 1e-3;

t_az(t_az >= 180) = t_az(t_az >= 180) - 360;

% Adjust for desired angular resolution and spacing
az_spacing = 0.5;
az_swath = 2.0;

% Keeps only the last samples 
pulse_mask = 1:(find(t_az > t_az(end), 1, 'first') - 1);
t_az(pulse_mask) = [];
t_el(pulse_mask) = [];
t_iqh(:, pulse_mask) = [];
t_iqv(:, pulse_mask) = [];

% Azimuth binning
az_discrete = round(t_az / az_spacing) * az_spacing;
az_unique = unique(az_discrete);
n_az = numel(az_unique);

% Computes the number of samples per radial
az_info = repmat(struct('pulse_ids',[],'az',[],'n_pulses',[]), [1, n_az]);
for i = 1:n_az
    az_info(i).pulse_ids = find(t_az < az_unique(i) + 0.5 * az_swath & t_az > az_unique(i) - 0.5 * az_swath);
    az_info(i).az = az_unique(i);
    az_info(i).n_pulses = numel(az_info(i).pulse_ids);
end

R0.H = nan(n_gt, n_az);
R0.V = nan(n_gt, n_az);
R1.H = nan(n_gt, n_az);
R0.X = nan(n_gt, n_az);

for i = 1:n_az
    tV.H = t_iqh(gt_range, az_info(i).pulse_ids);
    tV.V = t_iqv(gt_range, az_info(i).pulse_ids);
    R0.H(:,i) = mean(bsxfun(@times, tV.H, conj(tV.H)), 2);
    R0.V(:,i) = mean(bsxfun(@times, tV.V, conj(tV.V)), 2);
    R1.H(:,i) = mean(bsxfun(@times, tV.H(:, 2:end), conj(tV.H(:, 1:end-1))), 2);
    R0.X(:,i) = mean(bsxfun(@times, tV.H, conj(tV.V)), 2);
end

% DPSD Processing
% Choose an azimuth (in deg)
az_deg = 27;
i_az = find([az_info(:).az] == az_deg, 1, 'first');
% az_info(i_az).az % Verify selected azimuth

% Select number of spectra to average for DPSD estimation
K = 1;
dpsd_strat = 1; % 0 for range, 1 for az averaging
NFT = [];

switch dpsd_strat
    case 0
        % % Range averaging
        % This portion of the code reshapes the I/Q matrix to include
        % additional signals for averaging. For example: if K = 3, it
        % groups range gates (1:3) for the first range gate, (2:4) for the
        % second range gate, and so on.
        range_set = hankel(fliplr(gt_range));
        range_set(range_set == 0) = 1;
        range_set = range_set(1:K, :);
        range_set = range_set(end:-1:1);
        
        % Reshape I/Q to include 
        tV.H = t_iqh(range_set, az_info(i_az).pulse_ids);
        tV.V = t_iqv(range_set, az_info(i_az).pulse_ids);
        
        % Set data windowing function
        w = blackman(az_info(i_az).n_pulses).';
    case 1
        % % Azimuth averaging
        % This does the same as previous portion but with using data from
        % different radials instead of range.
        az_set = floor((1:K) - K/2) + i_az;
        az_set(az_set < 1) = 1;
        az_set(az_set > n_az) = n_az;
        
        % Pad to min pulses
        min_pulses = min([az_info(az_set).n_pulses]);
        
        tV.H = zeros(n_gt * K, min_pulses);
        tV.V = zeros(n_gt * K, min_pulses);
        for i = 1:K
            t_cropped_ids = az_info(az_set(i)).pulse_ids.';
            tV.H(i:K:end, 1:min_pulses) = t_iqh(gt_range, t_cropped_ids(1:min_pulses));
            tV.V(i:K:end, 1:min_pulses) = t_iqv(gt_range, t_cropped_ids(1:min_pulses));
        end
        
        % Set data windowing function
        w = blackman(min_pulses).';
end

% Get the Bootstrap DPSD estimate
E = bootstrap_dpsd(tV, w, N0, NFT, 20, K, n_gt);

sp1 = 10*log10(E.sSNR.H);
sp2 = 10*log10(E.sSNR.V);
sp3 = 10*log10(E.sD);
sp4 = E.sR;

snrmask = sp1 > 20 & sp2 > 20;
sp1(~snrmask) = nan;
sp2(~snrmask) = nan;
sp3(~snrmask) = nan;
sp4(~snrmask) = nan;

v_axis = get_velocity_axis(va, az_info(i_az).n_pulses);

% Plot range-Doppler of sSNR H/V, sZDR, and sRHV
figure(1);
imagesc(v_axis, r_unique, sp1);
set(gca,'ydir','normal');
colormap(boonlib('zmap'));
caxis([0 80]);
freezeColors;
set(gca, 'fontsize', 14);
xlabel('Radial velocity (m/s)', 'fontsize', 14);
ylabel('Range (km)', 'fontsize', 14);
title('sSNR_H (dB)')

figure(2);
imagesc(v_axis, r_unique, sp2);
set(gca,'ydir','normal');
colormap(boonlib('zmap'));
caxis([0 80]);
freezeColors;
set(gca, 'fontsize', 14);
xlabel('Radial velocity (m/s)', 'fontsize', 14);
ylabel('Range (km)', 'fontsize', 14);
title('sSNR_V (dB)')

figure(3);
imagesc(v_axis, r_unique, sp3);
set(gca,'ydir','normal');
colormap(dmap_wsr(256));
caxis([-5 8]);
freezeColors;
set(gca, 'fontsize', 14);
xlabel('Radial velocity (m/s)', 'fontsize', 14);
ylabel('Range (km)', 'fontsize', 14);
title('sZ_{DR} (dB)')

figure(4);
imagesc(v_axis, r_unique, sp4);
set(gca,'ydir','normal');
colormap(rmap_wsr(256));
caxis([0.1 1.05]);
set(gca, 'fontsize', 14);
xlabel('Radial velocity (m/s)', 'fontsize', 14);
ylabel('Range (km)', 'fontsize', 14);
title('s\rho_{HV}')