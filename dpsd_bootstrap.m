% why is this happening to me
close all
clear iqh1 iqh2 iqv1 iqv2 iqh iqv

plot_save_flag = 0; % set plot_save_flag = 1 to save plots from this script

K = 20; % number of bootstrapped pseudorealizations

[iq_source, signal_type] = choose_source();

%-----------If using iq_emulate.m------------%
if strcmp(iq_source, 'emulator')
    
    %---Set desired parameter values---%
    %---CHANGE THESE AS NEEDED!!!---%
    rlzs = 1000; % Number of time series to emulate
    M = 64; % Number of samples per realization
    S = [30, 30]; % Signal power in dB
    vr = [-35, 25]; % mean Doppler velocity in m/s
    sw = [5, 3]; % Do5pler spectrum width in m/s
    SNR = 30; % SNR in dB
    zdr = [0, 3]; % Bulk ZDR in dBZ
    phv = [0.98, 0.20]; % Bulk rhoHV
    
    va = 50; % Nyquist velocity
    lambda = 0.1; % wavelength
    Ts = lambda/4/va;
    
    if strcmp(signal_type, 'rain')
        [iqh, iqv, ZDR_bulk, PHV_bulk] = iq_emulate(M, S(1), vr(1), sw(1), SNR, va, rlzs, zdr(1), phv(1));
    elseif strcmp(signal_type, 'debris')
        [iqh, iqv, ZDR_bulk, PHV_bulk] = iq_emulate(M, S(2), vr(2), sw(2), SNR, va, rlzs, zdr(2), phv(2));
    elseif strcmp(signal_type, 'multi')
        [iqh1, iqv1, ZDR_bulk1, PHV_bulk1] = iq_emulate(M, S(1), vr(1), sw(1), SNR, va, rlzs, zdr(1), phv(1));
        [iqh2, iqv2, ZDR_bulk2, PHV_bulk2] = iq_emulate(M, S(2), vr(2), sw(2), SNR, va, rlzs, zdr(2), phv(2));
        iqh = iqh1 + iqh2;
        iqv = iqv1 + iqv2;
    end

%---------If using SimRadar I/Q data---------%
elseif strcmp(iq_source, 'simradar')
    
    if strcmp(signal_type, 'rain')
        filename = '~/Documents/sims/suctvort/200326/nodebris/sim-PPI0.5-DU-nodebris';
        load([filename '.mat']);
    elseif strcmp(signal_type, 'debris')
        filename = '~/Documents/sims/suctvort/200326/debris3/sim-PPI0.5-TCU-d3n100000';
        load([filename '.mat']);
    elseif strcmp(signal_type, 'multi')
        filename = '~/Documents/sims/suctvort/200326/debris3/sim-PPI0.5-DCU-d3n100000';
        load([filename '.mat']);
    end
    
    M = size(iqh, 3);
    lambda = dat.params.lambda;
    Ts = dat.params.prt;
    va = dat.params.va;
    rlzs = 1;
    
    iqh_full = iqh;
    iqv_full = iqv;
    r_ind = randsample(size(iqh,1), 1); % randomly pull from SimRadar output
    az_ind = randsample(size(iqh,2), 1);
    iqh = squeeze(iqh(r_ind, az_ind, :));
    iqv = squeeze(iqv(r_ind, az_ind, :));
    
    
    figure(101)
    ha = subplot(1,2,1);
    hs = pcolor(xx, yy, zh(:,:,1));
    hold on
        plot(xx(r_ind, az_ind), yy(r_ind, az_ind), 'sk', 'MarkerSize', 11, 'MarkerFaceColor', 'w', 'LineWidth', 2);
    hold off
    set(gca, 'DataAspect', [1 1 1])
    caxis([0 80])
    colormap(ha, blib('zmap'))
    colorbar
    shading flat
    set(gca, 'YDir', 'Normal')
    title('Reflectivity (dBZ)')
    
    ha(2) = subplot(1,2,2);
    hs(2) = pcolor(xx, yy, vr(:,:,1));
    hold on
        plot(xx(r_ind, az_ind), yy(r_ind, az_ind), 'sk', 'MarkerSize', 11, 'MarkerFaceColor', 'w', 'LineWidth', 2);
    hold off
    set(gca, 'DataAspect', [1 1 1])
    caxis([-1 1] * round(max(max(abs(vr(:,:,1)))), -1))
    colormap(ha(2), blib('rgmap2'))
    colorbar
    shading flat
    title('Velocity (m/s)')
    
    axes('Unit', 'Normalized', 'Position', [0.5 0.8 0.01 0.01])
    title(['SimRadar-generated spectrum (r=', num2str(round(r(r_ind)/1000,2)),...
        'km, \phi_{az}=' num2str(round(180/pi*az_rad(az_ind),1)) '^o)'], 'FontSize', 12)
    axis off
end

%------------------------

% Choose data window
% d = hamming(M);
% d = blackman(M);
d = hann(M);
% d = rectwin(M);

% Figure out how to ?????????

Rxx_mat = zeros(2*M-1, rlzs);
Sxx_mat = Rxx_mat;
for n = 1:rlzs
    [Rxx_mat(:,n), ~] = xcorr(squeeze(iqh(:,n))); % ACF
    Sxx_mat(:,n) = fftshift(fft(Rxx_mat(:,n))); % PSD
end
inds_2s = -(M-1):(M-1);
lags = inds_2s * Ts;

fs = 2*va/lambda;
f0 = fs/M;
f = f0 * inds_2s;
vv = -f * lambda / 2;

if strcmp(iq_source, 'emulator')
    Rxx = mean(Rxx_mat,2);
    Sxx = mean(Sxx_mat,2);
    R0h = mean(iqh.*conj(iqh),1); % Horiz. lag 0 ACF
    R0v = mean(iqv.*conj(iqv),1); % Vert. lag 0 ACF
elseif strcmp(iq_source, 'simradar')
    Rxx = Rxx_mat;
    Sxx = Sxx_mat;
    R0h = mean(iqh.*conj(iqh),1);
    R0v = mean(iqv.*conj(iqv),1);
end



% ---Signal coherency correction---

CX_HL = 0.5 * iqh(1,:) ./ iqh(M,:);
CX_HR = 0.5 * iqh(M,:) ./ iqh(1,:);
CX_VL = 0.5 * iqv(1,:) ./ iqv(M,:);
CX_VR = 0.5 * iqv(M,:) ./ iqv(1,:);
CX_L = CX_HL + CX_VL; % total left coherency correction factor
CX_R = CX_HR + CX_VR; % total right coherency correction factor
CX_L = repmat(CX_L, [M-1 1]);
CX_R = repmat(CX_R, [M-1 1]);

iqh_L = CX_L .* iqh(1:M-1, :);
iqh_R = CX_R .* iqh(2:M, :);
iqv_L = CX_L .* iqv(1:M-1, :);
iqv_R = CX_R .* iqv(2:M, :);

XH = [iqh_L; iqh; iqh_R];
XV = [iqv_L; iqv; iqv_R];

% XH = [sqrt(R0h ./ mean(iqh_L.*conj(iqh_L),1)) .* iqh_L;
%     iqh;
%     sqrt(R0h ./ mean(iqh_R.*conj(iqh_R),1)) .* iqh_R]; % combine into extended H and V time series
% XV = [sqrt(R0v ./ mean(iqv_L.*conj(iqv_L),1)) .* iqv_L;
%     iqv;
%     sqrt(R0v ./ mean(iqv_R.*conj(iqv_R),1)) .* iqv_R];

ts_plot_ind = randsample(rlzs, 1);

figure(1)
subplot(2,2,1)
plot(1:M, iqh(:,ts_plot_ind), '-b', 1:M, iqv(:,ts_plot_ind), '-r')
axis square
title('Time series')
legend('H', 'V', 'Location', 'northeast')

subplot(2,2,2)
plot(1:3*M-2, XH(:,ts_plot_ind), '-b', 1:3*M-2, XV(:,ts_plot_ind), '-r')
axis square
title('Extended time series')
legend('H', 'V', 'Location', 'northeast')

subplot(2,2,3)
semilogy(lags, abs(Rxx))
axis square
if strcmp(iq_source, 'emulator')
    title(['Mean ACF (' num2str(rlzs) ' rlzs)'])
elseif strcmp(iq_source, 'simradar')
    title('ACF')
end
xlabel('Lag {\it l}')

subplot(2,2,4)
semilogy(vv, abs(Sxx))
axis square
xlim([-va va])
if strcmp(iq_source, 'emulator')
    title(['Mean PSD (' num2str(rlzs) ' rlzs)'])
elseif strcmp(iq_source, 'simradar')
    title('PSD')
end
xlabel('Doppler velocity {\it v_r}')



% Maximum ratio of corrected samples
a = 1/M * sum( (d/max(d)).^2 );
rmax = (1 - sqrt(a)) / 2;

cor = [ones(M-1,1); zeros(M,1); ones(M-1,1)];
% cor = 1 for indices in X_H and X_V where coherency correction was applied,
% 0 for indices without correction (i.e. the original signal)
r = zeros(2*M-1,1);
for m = 1:2*M-1
    % for every possible M-length block of data from extended time series,
    % calculate the number of corrected samples in the block
    r(m) = sum(cor(m: m+M-1))/M;
end

block_start_inds = find(r <= rmax); % find the starting indices of all desirable blocks
% These are the indices from which the bootstrapped pseudorealizations will be generated.

figure(2)
plot(1:3*M-2, XH(:,ts_plot_ind), 'Color', [0.10 0.10 0.75])
hold on
plot(1:3*M-2, XV(:,ts_plot_ind), 'Color', [0.75 0.25 0.25])
plot([block_start_inds(1) block_start_inds(1)], ylim,...
    '-k', 'LineWidth', 1.5)
plot([block_start_inds(end)+M-1, block_start_inds(end)+M-1],...
    ylim, '-k', 'LineWidth', 1.5)
plot([M, M], ylim, 'Color', [0.25 0.50 0.25], 'LineWidth', 1.5)
plot([2*M, 2*M], ylim, 'Color', [0.25 0.50 0.25], 'LineWidth', 1.5)
hold off
title('Extended time series bootstrapping domain')
xlabel('Original time series', 'Color', [0.25 0.50 0.25])
legend('H', 'V', 'Location', 'northeast')


% Bootstrap
j = randsample(block_start_inds, K, true); % bootstrap!

ZH = zeros(M,K,rlzs);
ZV = ZH;
sSH = ZH;
sSV = ZH;
sSX = ZH;
ks = 0:M-1;
for n = 1:K
    % Power correction
    BH = XH(j(n):j(n)+M-1, :);
    BV = XV(j(n):j(n)+M-1, :);
    RH = mean(BH.*conj(BH), 1);
    RV = mean(BV.*conj(BV), 1);
    BH = repmat(sqrt(R0h ./ RH), [M 1]) .* BH;
    BV = repmat(sqrt(R0v ./ RV), [M 1]) .* BV;
    
    % Compute PSDs
    for k = 1:M
        ZH(k,n,:) = sum( repmat(d,[1 rlzs]) .* BH .* exp(-1j*2*pi*(k-1)*repmat(ks(:),[1 rlzs])/M) , 1 );
        ZV(k,n,:) = sum( repmat(d,[1 rlzs]) .* BV .* exp(-1j*2*pi*(k-1)*repmat(ks(:),[1 rlzs])/M) , 1 );
    end
    sSH(:,n,:) = ZH(:,n,:) .* conj(ZH(:,n,:)) / M;
    sSV(:,n,:) = ZV(:,n,:) .* conj(ZV(:,n,:)) / M;
    sSX(:,n,:) = ZH(:,n,:) .* conj(ZV(:,n,:)) / M;
end

% Compute DPSD estimates
sZDR = zeros(M,rlzs);
sPHV = sZDR;
for k = 1:M
    sZDR(k,:) = squeeze(mean(sSH(k,:,:),2)) ./ squeeze(mean(sSV(k,:,:),2));
    sPHV(k,:) = abs(squeeze(mean(sSX(k,:,:),2))) ./ sqrt(squeeze(mean(sSH(k,:,:),2)) .* squeeze(mean(sSV(k,:,:),2)));
end

% Bias correction
if K == 1
    b = ((1-rmax)^(-3.3)) - 2*((1-rmax)^1.1);
elseif K > 1
    b = ((1-rmax)^(-4.5)) - 2*((1-rmax)^(-2.1));
end

sZDR_corr = sZDR .* (1 - ((1/b/K) * (1 - sPHV.^2)));
sPHV_corr = sPHV .* (1 - ((1/b/K) * (1 - sPHV.^2).^2 ./ (4*sPHV.^2)));

sSH = mean(sSH, 3);
sSV = mean(sSV, 3);
sZDR = mean(sZDR, 2);
sPHV = mean(sPHV, 2);
sZDR_corr = mean(sZDR_corr, 2);
sPHV_corr = mean(sPHV_corr, 2);

% Velocity vector for plotting
%vv = -va * linspace(-1, 1, M);
vv = -2 * flip(rot90((0:M-1))) * va / M;
vv(vv < -va) = vv(vv < -va) + (2*va);


sort_mat = sortrows([vv, sSH, sSV, sZDR, sPHV, sZDR_corr, sPHV_corr]);
vv = sort_mat(:,1);
sSH = sort_mat(:,2:K+1);
sSV = sort_mat(:,K+2:2*K+1);
sZDR = sort_mat(:,2*K+2);
sPHV = sort_mat(:,2*K+3);
sZDR_corr = sort_mat(:,2*K+4);
sPHV_corr = sort_mat(:,2*K+5);


% plot ZDR spectrum
figure(3)
subplot(2,2,1)
semilogy(vv, mean(abs(sSH),2))
xlim([-va va])
axis square
grid on
title(['Mean H-channel PSD ({\it K}=', num2str(K), ')'])
xlabel('Doppler velocity {\it v_r}')

subplot(2,2,2)
semilogy(vv, mean(abs(sSV),2))
xlim([-va va])
axis square
grid on
title(['Mean V-channel PSD ({\it K}=', num2str(K), ')'])
xlabel('Doppler velocity {\it v_r}')

subplot(2,2,3)
plot(vv, 10*log10(sZDR))
xlim([-va, va])
ylim([-6 6])
axis square
grid on
title('{\it s}Z_{DR}({\it v_r}) in dBZ')
xlabel('Doppler velocity {\it v_r}')

subplot(2,2,4)
plot(vv, 10*log10(sZDR_corr))
xlim([-va, va])
ylim([-6 6])
axis square
grid on
title('Bias-corrected {\it s}Z_{DR}({\it v_r}) in dBZ')
xlabel('Doppler velocity {\it v_r}')


% plot rhoHV spectrum
figure(4)
subplot(2,2,1)
semilogy(vv, mean(abs(sSH),2))
xlim([-va va])
axis square
grid on
title(['Mean H-channel PSD ({\it K}=', num2str(K), ')'])
xlabel('Doppler velocity {\it v_r}')

subplot(2,2,2)
semilogy(vv, mean(abs(sSV),2))
xlim([-va va])
axis square
grid on
title(['Mean V-channel PSD ({\it K}=', num2str(K), ')'])
xlabel('Doppler velocity {\it v_r}')

subplot(2,2,3)
plot(vv, sPHV)
xlim([-va, va])
ylim([0 1])
axis square
grid on
title('{\it s}\rho_{HV}({\it v_r})')
xlabel('Doppler velocity {\it v_r}')

subplot(2,2,4)
plot(vv, sPHV_corr)
xlim([-va, va])
ylim([0 1])
axis square
grid on
title('Bias-corrected {\it s}\rho_{HV}({\it v_r})')
xlabel('Doppler velocity {\it v_r}')

