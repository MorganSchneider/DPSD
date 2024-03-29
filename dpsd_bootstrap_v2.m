% I fucking BROKE it AND THIS IS WHY WE DO VERSION CONTROL
clear;
%close all
clear iqh1 iqh2 iqv1 iqv2 iqh iqv

plot_save_flag = 0; % set plot_save_flag = 1 to save plots from this script
K = 20; % number of bootstrapped pseudorealizations

[iq_source, signal_type] = choose_source();
% iq_source = 'simradar';
% signal_type = 'multi';
% if strcmp(signal_type, 'rain') || strcmp(signal_type, 'multi')
%     cx = 'DCU';
% else
%     cx = 'TC';
% end
% el = {'2.0', '2.5', '3.0', '3.5', '4.0', '4.5', '5.0'};

zmap = feval('boonlib', 'zmap', 21);
rgmap = feval('boonlib', 'rgmap', 21);


%---------If using iq_emulate_yu.m---------%
if strcmp(iq_source, 'yu')
    
    % Index 1 is rain setting
    % Index 2 is debris setting
    S_dB = [30 30];
    S = 10.^(S_dB/10);
    sw = [5 10];
    vr = [-25 25];
    zdr = [3 -2];
    phv = [0.98 0.5];
    rlzs = 1; % Number of spectra to generate
    M = 100;
    va = 50;
%     % vv2 = linspace(-2*va, 2*va, 2*M-1); % For Folding spectra
    vv = linspace(-va+1, va, M);
    lambda = 0.1; % wavelength
    Ts = lambda/4/va;
    
    Sv(1,:) = S(1) ./ sqrt(2*pi*sw(1)^2) .* exp(-(vv-vr(1)).^2 / (2*sw(1)^2));
    Sv(2,:) = S(2) ./ sqrt(2*pi*sw(2)^2) .* exp(-(vv-vr(2)).^2 / (2*sw(2)^2));
    
    if strcmp(signal_type, 'rain')
        Sv = Sv(1,:);
        sphv = phv(1) * ones(1,M);
        szdr = zdr(1) * ones(1,M);
    elseif strcmp(signal_type, 'debris')
        Sv = Sv(2,:);
        sphv = phv(2) * ones(1,M);
        szdr = zdr(2) * ones(1,M);
    else
        Sv = sum(Sv,1);
        
        % Adjust for Fourier component order
        [c, ind] = min(abs(vv - (vr(1)+vr(2)) / 2));
        if vr(1) < vr(2)
            sphv(1:ind) = phv(1);
            sphv(ind+1: M) = phv(2);
            szdr(1:ind) = zdr(1);
            szdr(ind+1: M) = zdr(2);
        else
            sphv(1:ind) = phv(2);
            sphv(ind+1: M) = phv(1);
            szdr(1:ind) = zdr(2);
            szdr(ind+1: M) = zdr(1);
        end
    end
    
    
    figure(20)
    subplot(1,2,1)
    yyaxis left
    semilogy(vv, Sv, 'k', 'LineWidth', 1)
    ylabel('S_V')
    yyaxis right
    plot(vv, szdr, '--b', 'LineWidth', 1)
    ylabel('sZ_{DR}')
    ylim([min(zdr)-1 max(zdr)+1])
    grid on
    axis square
    
    subplot(1,2,2)
    yyaxis left
    semilogy(vv, Sv, 'k', 'LineWidth', 1)
    ylabel('S_V')
    yyaxis right
    plot(vv, sphv, '--r', 'LineWidth', 1)
    ylabel('s\rho_{HV}')
    ylim([0 1])
    grid on
    axis square
    
    szdr_linear = 10 .^ (szdr./10);
    [iqh, iqv, ZDR_bulk, PHV_bulk] = iq_emulate_yu(Sv, sphv, szdr_linear, rlzs);

    % Calculate error (RMSE) between specified and estimated DPSDs & play
    % with parameters like rlzs, etc
    % Check high and low rhoHV
    % Averaging over time to smooth (averaging "chunks" of time series
    % points)
    % Small azimuthal change (-0.5-0.5 degrees) and lots of pulses (like
    % 10000 pulses) to keep signal from decorrelating, then average over
    % dwells
    
    szdr_save = szdr;
    sphv_save = sphv; 
    clear sphv szdr % Clear these for later
    
    params = struct('signal_type', signal_type, 'S_dB', S_dB, 'S', S, 'Sv', Sv, 'vr', vr, 'sw', sw,...
        'zdr', zdr, 'phv', phv, 'szdr', szdr_save, 'sphv', sphv_save, 'K', K);
    

%---------If using SimRadar I/Q data---------%
elseif strcmp(iq_source, 'simradar')
    
    sim_dir = ['~/Documents/code/DPSD/test_sims/' signal_type];
    filename = blib('choosefile', sim_dir, '*.mat');
    load(filename);
%     filename = ['sim-PPI' el{elx} '-' cx '-d' num2str(dt) 'n' num2str(dn) '.mat'];
%     %filename = ['sim-PPI' el{elx} '-' cx '-nodebris.mat'];
%     load([sim_dir '/' filename]);
    
%     if strcmp(signal_type, 'rain')
%         sim_dir = '~/Documents/code/DPSD/test_sims/rain';
%         filename = blib('choosefile', sim_dir, '*.mat');
%         load(filename);
%     elseif strcmp(signal_type, 'debris')
%         sim_dir = '~/Documents/code/DPSD/test_sims/debris';
%         filename = blib('choosefile', sim_dir, '*.mat');
%         load(filename);
%     elseif strcmp(signal_type, 'realistic')
%         sim_dir = '~/Documents/code/DPSD/test_sims/multi';
%         filename = blib('choosefile', sim_dir, '*.mat');
%         load(filename);
%     end
    
    params = struct('signal_type', signal_type, 'zdr', zdr, 'phv', rhohv, 'K', K,...
        'xx', xx, 'yy', yy, 'zz', zz, 'va', dat.params.va, 'vr', vr);
    
    if size(iqh,4) > 1
        iqh = iqh(:,:,:,1);
        iqv = iqv(:,:,:,1);
    end
    
    iqh = permute(iqh, [3 1 2]);
    iqv = permute(iqv, [3 1 2]);
    
    M = size(iqh, 1);
    rlzs = 1;
    
    lambda = dat.params.lambda;
    Ts = dat.params.prt;
    va = dat.params.va;
    vv = linspace(-va+1, va, M);
end



%% ------------------------

vvx = flip(vv);

% Choose data window
d = hann(M); % hamming(M), hann(M), blackman(M), rectwin(M)
R0_d = mean(abs(d).^2); % window power

dm2 = size(iqh,2);
dm3 = size(iqh,3);
% Random indices for plotting
%ind1 = randsample(dm2,1);
%ind2 = randsample(dm3,1);

if strcmp(iq_source, 'simradar')
    ind1 = 12;
    ind2 = 19;
else
    ind1 = 1;
    ind2 = 1;
end

% Data window
d_rep = repmat(d, [1 dm2 dm3]);
mm = repmat((0:M-1)', [1 dm2 dm3]);

% Calculate DPSD the Original Way prior to bootstrap
DFT_H = fftshift(fft(d_rep .* iqh, M, 1));
DFT_V = fftshift(fft(d_rep .* iqv, M, 1));
% If emulator_yu, 2 for average in dwell
% If simradar,  2 for average in range
%               3 for average in azimuth
SHF = nanmean(abs(DFT_H).^2 / M, 2);
SVF = nanmean(abs(DFT_V).^2 / M, 2);
SXF = nanmean(DFT_H .* conj(DFT_V) / M, 2);

figure(152)

subplot(2,1,1)
yyaxis left
plot(vvx, 10*log10(SHF(:,ind2)), '-k', vvx, 10*log10(SVF(:,ind2)), '-b')
xlabel('V_{r} (m s^{-1})', 'FontSize', 14)
ylabel('S_{H,V}', 'FontSize', 14)
yyaxis right
plot(vvx, 10*log10(SHF(:,ind2) ./ SVF(:,ind2)), '--k')
ylabel('S(Z_{DR}) dB', 'FontSize', 14)

subplot(2,1,2)
yyaxis left
plot(vvx, 10*log10(SHF(:,ind2)), '-k', vvx, 10*log10(SVF(:,ind2)), '-b')
xlabel('V_{r} (m s^{-1})', 'FontSize', 14)
yyaxis right
plot(vvx, abs(SXF(:,ind2)) ./ sqrt(SHF(:,ind2) .* SVF(:,ind2)), '--k')
ylabel('S(\rho_{HV})', 'FontSize', 14)

if plot_save_flag
    print(['~/Documents/imgs/DPSD/' iq_source '_' signal_type '_oldDPSD'], '-dpng')
end

Rxx_mat = zeros(2*M-1, dm2, dm3);
for m = 1:dm2
    for n = 1:dm3
        [Rxx_mat(:,m,n), ~] = xcorr(squeeze(iqh(:,m,n))); % ACF
    end
end
Sxx_mat = fftshift(fft(Rxx_mat, 2*M-1, 1), 1);

inds_2s = -(M-1):(M-1);
lags = inds_2s * Ts;

fs = 2 * va / lambda;
f0 = fs / M;
f = f0 * inds_2s;
vv1 = -f * lambda / 2;

if dm3 == 1
    Rxx = mean(Rxx_mat, 2);
    Sxx = mean(Sxx_mat, 2);
else
    Rxx = Rxx_mat;
    Sxx = Sxx_mat;
end
R0h = mean(iqh.*conj(iqh), 1); % Horiz. lag 0 ACF
R0v = mean(iqv.*conj(iqv), 1); % Vert. lag 0 ACF


%---Signal coherency correction---%
CX_L = repmat(0.5 * (iqh(1,:,:)./iqh(M,:,:) + iqv(1,:,:)./iqv(M,:,:)), [M-1 1 1]);
CX_R = repmat(0.5 * (iqh(M,:,:)./iqh(1,:,:) + iqv(M,:,:)./iqv(1,:,:)), [M-1 1 1]);
iqh_L = CX_L .* iqh(1:M-1,:,:);
iqh_R = CX_R .* iqh(2:M,:,:);
iqv_L = CX_L .* iqv(1:M-1,:,:);
iqv_R = CX_R .* iqv(2:M,:,:);

XH = [iqh_L; iqh; iqh_R];
XV = [iqv_L; iqv; iqv_R];


figure(1)
subplot(2,2,1)
plot(1:M, iqh(:,ind1,ind2), '-b', 1:M, iqv(:,ind1,ind2), '-r')
axis square
title('Time series')
legend('H', 'V', 'Location', 'northeast')

subplot(2,2,2)
plot(1:3*M-2, XH(:,ind1,ind2), '-b', 1:3*M-2, XV(:,ind1,ind2), '-r')
axis square
title('Extended time series')
legend('H', 'V', 'Location', 'northeast')

subplot(2,2,3)
semilogy(lags, abs(Rxx(:,ind1,ind2)))
axis square
if strcmp(iq_source, 'yu')
    title(['Mean ACF (' num2str(rlzs) ' rlzs)'])
elseif strcmp(iq_source, 'simradar')
    title('Sample ACF')
end
xlabel('Lag {\it l}')

subplot(2,2,4)
semilogy(vv1, abs(Sxx(:,ind1,ind2)))
axis square
xlim([-va va])
if strcmp(iq_source, 'yu')
    title(['Mean PSD (' num2str(rlzs) ' rlzs)'])
elseif strcmp(iq_source, 'simradar')
    title('Sample PSD')
end
xlabel('Doppler velocity {\it v_r}')
if plot_save_flag
    print(['~/Documents/imgs/DPSD/' iq_source '_' signal_type '_TS-ACF-PSD'], '-dpng')
end

% Maximum ratio of corrected samples
a = mean((d/max(d)).^2); % this is slightly different in Arturo's code!!!
rmax = (1 - sqrt(a)) / 2;

cor = [ones(M-1,1); zeros(M,1); ones(M-1,1)];
% cor = 1 for indices in X_H and X_V where coherency correction was applied,
% 0 for indices without correction (i.e. the original signal)
r = zeros(2*M-1, 1);
for m = 1: 2*M-1
    % for every possible M-length block of data from extended time series,
    % calculate the number of corrected samples in the block
    r(m) = sum(cor(m : m+M-1)) / M;
end

block_start_inds = find(r <= rmax); % find the starting indices of all desirable blocks
% These are the indices from which the bootstrapped pseudorealizations will be generated.

figure(2)
plot(1:3*M-2, XH(:,ind1,ind2), 'Color', [0.10 0.10 0.75])
hold on
plot(1:3*M-2, XV(:,ind1,ind2), 'Color', [0.75 0.25 0.25])
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

if plot_save_flag
    print(['~/Documents/imgs/DPSD/' iq_source '_' signal_type '_bootstrapdomain'], '-dpng')
end

% Bootstrap
j = randsample(block_start_inds, K, true); % bootstrap!

% ZH has dimensions (length of orig. TS, # pseudorlzs, [# independent rlzs])
ZH = zeros(M,K,dm2,dm3);
ZV = ZH;
ks = 0:M-1;
for n = 1:K
    % Power correction
    BH = XH(j(n):j(n)+M-1, :, :);
    BV = XV(j(n):j(n)+M-1, :, :);
    RH = mean(BH.*conj(BH), 1);
    RV = mean(BV.*conj(BV), 1);
    BH = repmat(sqrt(R0h ./ RH), [M 1 1]) .* BH;
    BV = repmat(sqrt(R0v ./ RV), [M 1 1]) .* BV;
    
    % Compute PSDs
%     ZH(:,n,:,:) = fftshift(fft(BH .* d_rep, M, 1), 1);
%     ZV(:,n,:,:) = fftshift(fft(BV .* d_rep, M, 1), 1);
    for k = 1:M
        ZH(k,n,:,:) = sum( d_rep .* BH .* exp(-1j*2*pi*(k-1)*repmat(ks(:),[1 dm2 dm3])/M) , 1 );
        ZV(k,n,:,:) = sum( d_rep .* BV .* exp(-1j*2*pi*(k-1)*repmat(ks(:),[1 dm2 dm3])/M) , 1 );
    end
end
ZH = fftshift(ZH,1);
ZV = fftshift(ZV,1);
sSH.ms = squeeze(mean(ZH.*conj(ZH) / (M*R0_d), 2));
sSV.ms = squeeze(mean(ZV.*conj(ZV) / (M*R0_d), 2));
sSX.ms = squeeze(mean(ZH.*conj(ZV) / (M*R0_d), 2));

% Compute DPSD estimates
sZDR.ms = sSH.ms ./ sSV.ms;
sPHV.ms = abs(sSX.ms) ./ sqrt(sSH.ms .* sSV.ms);

% Bias correction
% WILL NEED TO CHANGE THIS TO RUN WITH SIMRADAR!!!!
if rlzs == 1
    b = (1-rmax)^-3.3 - 2*(1-rmax)^1.1;
elseif rlzs > 1
    b = (1-rmax)^-4.5 - (1-rmax)^-2.1;
end
sZDR.c = sZDR.ms .* (1 - ((1/b/rlzs) * (1 - (sPHV.ms).^2)));
sPHV.c = sPHV.ms .* (1 - ((1/b/rlzs) * (1 - (sPHV.ms).^2).^2 ./ (4*(sPHV.ms).^2)));
sPHV.c(sPHV.c < 0) = 0;
sZDR.c(sZDR.c < 0) = real(sZDR.c(sZDR.c < 0));

% Average over independent rlzs
if dm3 == 1
    sSH.f = squeeze(mean(sSH.ms, 2));
    sSV.f = squeeze(mean(sSV.ms, 2));
    sZDR.f = squeeze(mean(sZDR.ms, 2));
    sPHV.f = squeeze(mean(sPHV.ms, 2));
    sZDR.cf = squeeze(mean(sZDR.c, 2));
    sPHV.cf = squeeze(mean(sPHV.c, 2));
else
    sSH.f = sSH.ms;
    sSV.f = sSV.ms;
    sZDR.f = sZDR.ms;
    sPHV.f = sPHV.ms;
    sZDR.cf = sZDR.c;
    sPHV.cf = sPHV.c;
end


% % Compare to Arturo's code
% if dm3 > 1
%     iqh_tmp = reshape(iqh, [M dm2*dm3]);
%     iqv_tmp = reshape(iqv, [M dm2*dm3]);
% else
%     iqh_tmp = iqh;
%     iqv_tmp = iqv;
% end
% iqh_tmp = permute(iqh_tmp, [2 1]);
% iqv_tmp = permute(iqv_tmp, [2 1]);
% V = struct('H', iqh_tmp, 'V', iqv_tmp);
% N0 = struct('H', 1, 'V', 1);
% 
% E = bootstrap_dpsd(V, d, N0, [], K, 1, dm2*dm3);
% E.sD = permute(E.sD, [2 1]);
% E.sR = permute(E.sR, [2 1]);
% if dm3 > 1
%     sZDR.au = reshape(E.sD, [M dm2 dm3]);
%     sPHV.au = reshape(E.sR, [M dm2 dm3]);
%     sZDR.auf = sZDR.au;
%     sPHV.auf = sPHV.au;
% else
%     sZDR.au = E.sD;
%     sPHV.au = E.sR;
%     sZDR.auf = squeeze(mean(sZDR.au, 2));
%     sPHV.auf = squeeze(mean(sPHV.au, 2));
% end

%
%
%
%
%

% plot ZDR spectrum
figure(7)
subplot(2,2,1)
yyaxis left
semilogy(vvx, abs(sSH.f(:,ind1,ind2)), '-.k', 'LineWidth', 1.5)
xlabel('Doppler velocity {\it v_r}', 'FontSize', 20)
ylabel('PSD', 'FontSize', 20)
yyaxis right
plot(vvx, 10*log10(sZDR.f(:,ind1,ind2)), '-r', 'LineWidth', 1.2)
% hold on
% plot(vvx, 10*log10(sZDR.auf(:,ind1,ind2)), '-r')
% hold off
% legend('PSD', 'morgan DPSD', 'arturo DPSD', 'Location', 'southwest')
ylim([-10 10])
ylabel('sZ_{DR}', 'FontSize', 20)
xlim([-va va])
%axis square
grid on
title(['Mean H-channel PSD, sZ_{DR} ({\it K}=', num2str(K), ')'], 'FontSize', 20)
%xlabel('Doppler velocity {\it v_r}', 'FontSize', 20)


subplot(2,2,2)
yyaxis left
semilogy(vvx, abs(sSH.f(:,ind1,ind2)), '-.k', 'LineWidth', 1.5)
xlabel('Doppler velocity {\it v_r}', 'FontSize', 20)
ylabel('PSD', 'FontSize', 20)
yyaxis right
plot(vvx, 10*log10(sZDR.cf(:,ind1,ind2)), '-r', 'LineWidth', 1.2)
% hold on
% plot(vvx, 10*log10(sZDR.auf(:,ind1,ind2)), '-r')
% hold off
ylim([-10 10])
ylabel('sZ_{DR}', 'FontSize', 20)
xlim([-va va])
%axis square
grid on
title(['Mean H-channel PSD, sZ_{DR} ({\it K}=', num2str(K), ')'], 'FontSize', 20)
xlabel('Doppler velocity {\it v_r}', 'FontSize', 20)

subplot(2,2,3)
plot(vvx, 10*log10(sZDR.f(:,ind1,ind2)), 'b')
% hold on
% plot(vvx, 10*log10(sZDR.auf(:,ind1,ind2)), 'r')
% hold off
xlim([-va, va])
ylim([-10 10])
%axis square
grid on
title('{\it s}Z_{DR}({\it v_r}) in dBZ')
xlabel('Doppler velocity {\it v_r}')

subplot(2,2,4)
plot(vvx, 10*log10(sZDR.cf(:,ind1,ind2)), 'b')
% hold on
% plot(vvx, 10*log10(sZDR.auf(:,ind1,ind2)), 'r')
% hold off
xlim([-va, va])
ylim([-10 10])
%axis square
grid on
title('Bias-corrected {\it s}Z_{DR}({\it v_r}) in dBZ')
xlabel('Doppler velocity {\it v_r}')

set(gcf, 'Units', 'Inches', 'Position', [2 2 12 10])
if plot_save_flag
    print(['~/Documents/imgs/DPSD/' iq_source '_' signal_type '_sZDR'], '-dpng')
end

% plot rhoHV spectrum
figure(8)
subplot(2,2,1)
yyaxis left
semilogy(vvx, abs(sSH.f(:,ind1,ind2)), '-.k', 'LineWidth', 1.5)
xlabel('Doppler velocity {\it v_r}', 'FontSize', 20)
ylabel('PSD', 'FontSize', 20)
yyaxis right
plot(vvx, sPHV.f(:,ind1,ind2), '-r', 'LineWidth', 1.2)
% hold on
% plot(vvx, sPHV.auf(:,ind1,ind2), '-r')
% hold off
% legend('PSD', 'morgan DPSD', 'arturo DPSD', 'Location', 'southwest')
ylim([0 1])
ylabel('s\rho_{HV}', 'FontSize', 20)
%axis square
grid on
title(['Mean H-channel PSD, s\rho_{HV} ({\it K}=', num2str(K), ')'], 'FontSize', 20)

subplot(2,2,2)
yyaxis left
semilogy(vvx, abs(sSH.f(:,ind1,ind2)), '-.k', 'LineWidth', 1.5)
ylabel('PSD', 'FontSize', 20)
yyaxis right
plot(vvx, sPHV.cf(:,ind1,ind2), '-r', 'LineWidth', 1.2)
% hold on
% plot(vvx, sPHV.auf(:,ind1,ind2), '-r')
% hold off
ylim([0 1])
ylabel('s\rho_{HV}', 'FontSize', 20)
xlim([-va va])
%axis square
grid on
title(['Mean H-channel PSD, s\rho_{HV} ({\it K}=', num2str(K), ')'], 'FontSize', 20)
xlabel('Doppler velocity {\it v_r}', 'FontSize', 20)

subplot(2,2,3)
plot(vvx, sPHV.f(:,ind1,ind2), 'b')
% hold on
% plot(vvx, sPHV.auf(:,ind1,ind2), 'r')
% hold off
xlim([-va, va])
ylim([0 1])
%axis square
grid on
title('{\it s}\rho_{HV}({\it v_r})')
xlabel('Doppler velocity {\it v_r}')

subplot(2,2,4)
plot(vvx, sPHV.cf(:,ind1,ind2), 'b')
% hold on
% plot(vvx, sPHV.auf(:,ind1,ind2), 'r')
% hold off
xlim([-va, va])
ylim([0 1])
%axis square
grid on
title('Bias-corrected {\it s}\rho_{HV}({\it v_r})')
xlabel('Doppler velocity {\it v_r}')

set(gcf, 'Units', 'Inches', 'Position', [2 2 12 10])

if plot_save_flag
    print(['~/Documents/imgs/DPSD/' iq_source '_' signal_type '_sPHV'], '-dpng')
end


%%

if strcmp(iq_source, 'yu')
    figure(5)
    
    ax1 = subplot(2,2,1);
    % DPSD
    yyaxis left
    plot(vv, szdr_save, 'k', 'LineWidth', 1.1)
    ylabel('sZ_{DR} (dB)', 'FontSize', 14)
    xlim([-va, va])
    ylim([-5, 5])
    text(18, 1, 'Debris', 'FontSize', 20)
    text(-30, 1, 'Rain', 'FontSize', 20)
    text(-48, 3.4, 'sZ_{DR}', 'FontSize', 14)
    text(-48, -3.5, 'sS_V', 'FontSize', 14, 'Color', [0.6 0.6 0.6])
    % PSD
    yyaxis right
    hs1 = semilogy(vv, Sv, 'LineWidth', 1);
    hs1.Color = [0.6 0.6 0.6];
    hs1.LineStyle = '-.';
    ylabel('sS_V', 'FontSize', 14)
    ax1.YAxis(1).Color = 'k';
    ax1.YAxis(2).Color = [0.6 0.6 0.6];
    grid on
    title('(a) Input sZ_{DR}, sS_V', 'FontSize', 14)
    xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
    %legend('sZ_{DR}', 'sS_V', 'Location', 'northwest')
    
    
    ax2 = subplot(2,2,2);
    % DPSD
    yyaxis left
    plot(vv, sphv_save, 'k', 'LineWidth', 1.1)
    ylabel('s\rho_{HV}', 'FontSize', 14)
    xlim([-va, va])
    ylim([0, 1])
    text(18, 0.6, 'Debris', 'FontSize', 20)
    text(-30, 0.6, 'Rain', 'FontSize', 20)
    text(-48, 0.95, 's\rho_{HV}', 'FontSize', 14)
    text(-48, 0.15, 'sS_V', 'FontSize', 14, 'Color', [0.6 0.6 0.6])
    % PSD
    yyaxis right
    hs2 = semilogy(vv, Sv, 'LineWidth', 1);
    hs2.Color = [0.6 0.6 0.6];
    hs2.LineStyle = '-.';
    ylabel('sS_V', 'FontSize', 14)
    ax2.YAxis(1).Color = 'k';
    ax2.YAxis(2).Color = [0.6 0.6 0.6];
    grid on
    title('(b) Input s\rho_{HV}, sS_V', 'FontSize', 14)
    xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
    %legend('s\rho_{HV}', 'sS_V', 'Location', 'northwest')
    
    
    ax3 = subplot(2,2,3);
    % DPSD
    yyaxis left
    plot(vvx, 10*log10(sZDR.cf(:,ind1,ind2)), 'k', 'LineWidth', 1.1)
    ylabel('sZ_{DR} (dB)', 'FontSize', 14)
    xlim([-va, va])
    ylim([-15, 15])
    text(-48, -1.8, 'sZ_{DR}', 'FontSize', 14)
    text(-47, -13.5, 'sS_H', 'FontSize', 14, 'Color', [0.6 0.6 0.6])
    % PSD
    yyaxis right
    hs3 = semilogy(vvx, abs(sSH.f(:,ind1,ind2)), 'LineWidth', 1);
    hs3.Color = [0.6 0.6 0.6];
    hs3.LineStyle = '-.';
    ylabel('sS_H', 'FontSize', 14)
    ax3.YAxis(1).Color = 'k';
    ax3.YAxis(2).Color = [0.6 0.6 0.6];
    grid on
    title('(c) Output sZ_{DR}, sS_H', 'FontSize', 14)
    xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
    %legend('sZ_{DR}', 'sS_H', 'Location', 'northwest')
    
    
    ax4 = subplot(2,2,4);
    % DPSD
    yyaxis left
    plot(vvx, sPHV.cf(:,ind1,ind2), 'k', 'LineWidth', 1.1)
    ylabel('s\rho_{HV}', 'FontSize', 14)
    xlim([-va, va])
    ylim([0, 1])
    text(-48, 0.89, 's\rho_{HV}', 'FontSize', 14)
    text(-47, 0.06, 'sS_H', 'FontSize', 14, 'Color', [0.6 0.6 0.6])
    % PSD
    yyaxis right
    hs4 = semilogy(vvx, abs(sSH.f(:,ind1,ind2)), 'LineWidth', 1);
    hs4.Color = [0.6 0.6 0.6];
    hs4.LineStyle = '-.';
    ylabel('sS_H', 'FontSize', 14)
    ax4.YAxis(1).Color = 'k';
    ax4.YAxis(2).Color = [0.6 0.6 0.6];
    grid on
    title('(d) Output s\rho_{HV}, sS_H', 'FontSize', 14)
    xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
    
    print('~/Documents/emulator-ex', '-dpng')
    
end

if strcmp(iq_source, 'simradar')
    
    figure(5)
    clf
    
    ax1 = subplot(2,1,1);
    % DPSD
    yyaxis left
    hs = plot(vvx, real(10*log10(sZDR.cf(:,ind1,ind2))), 'LineWidth', 1.4);
    hs.Color = [0.0 0.0 0.7];
    hs.LineStyle = '-';
    ylabel('sZ_{DR} (dB)', 'FontSize', 16)
    xlim([-va, va])
    ylim([-20, 20])
%     text(-95, 11, 'sZ_{DR}', 'FontSize', 17, 'Color', [0.0 0.0 0.7])
%     text(-87, -17, 'sS_H', 'FontSize', 18, 'Color', [0.3 0.3 0.3])
    % PSD
    yyaxis right
    hs1 = semilogy(vvx, real(abs(sSH.f(:,ind1,ind2))), 'LineWidth', 2);
    %hs1.Color = [0.6 0.6 0.6];
    %hs1.LineStyle = '-.';
    hs1.Color = [0.5 0.5 0.5];
    hs1.LineStyle = ':';
    ylabel('sS_H', 'FontSize', 14)
    ax1.YAxis(1).Color = [0.0 0.0 0.7];
    ax1.YAxis(2).Color = 'r';
    %grid on
    %title('(a) sZ_{DR}, sS_H', 'FontSize', 14)
    title('Spectral differential reflectivity (sZ_{DR})', 'FontSize', 16)
    %xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
    %legend('sZ_{DR}', 'sS_H', 'Location', 'northwest')
    
    ax2 = subplot(2,1,2);
    % DPSD
    yyaxis left
    hs = plot(vvx, real(sPHV.cf(:,ind1,ind2)), 'LineWidth', 1.4);
    hs.Color = [0.7 0.0 0.0];
    ylabel('s\rho_{HV}', 'FontSize', 16)
    xlim([-va, va])
    ylim([-0.2, 1.2])
    yticks(0:0.2:1)
    text(-92, 1.1, 's\rho_{HV}', 'FontSize', 18, 'Color', [0.6 0.0 0.0])
    text(-87, -0.09, 'sS_H', 'FontSize', 18, 'Color', [0.3 0.3 0.3])
    % PSD
    yyaxis right
    hs2 = semilogy(vvx, real(abs(sSH.f(:,ind1,ind2))), 'LineWidth', 2);
    %hs2.Color = [0.6, 0.6, 0.6];
    %hs2.LineStyle = '-.';
    hs2.Color = [0.5 0.5 0.5];
    hs2.LineStyle = ':';
    ylabel('sS_H', 'FontSize', 14)
    ax2.YAxis(1).Color = [0.7 0.0 0.0];
    ax2.YAxis(2).Color = [0.3 0.3 0.3];
    %grid on
    %title('(b) s\rho_{HV}, sS_H', 'FontSize', 14)
    title('Spectral correlation coefficient (s\rho_{HV})', 'FontSize', 16)
    xlabel('{\it v} (m s^{-1})', 'FontSize', 16)
    %legend('s\rho_{HV}', 'sS_V', 'Location', 'northwest')
    
    if strcmp(signal_type, 'multi')
        signal_type = 'realistic';
    end
    sig_type = [upper(signal_type(1)), signal_type(2:end)];
    
    annotation('textbox', [0.32 0.94 0.4 0.05], 'String', [sig_type ' signal DPSDs'],...
       'FontSize',18,'FontWeight','bold','EdgeColor','k','LineWidth',1,...
       'HorizontalAlignment','center','VerticalAlignment','middle')
    
    set(gcf, 'Units', 'Inches', 'Position', [6 8 8 9])
    ax1.Position(2) = 0.54;
    
    print(['~/Documents/DPSD-' signal_type '-ex'], '-dpng')
    
end


% Save variables into .mat file
% if strcmp(iq_source, 'simradar')
%     filename = erase(filename, [sim_dir '/']);
%     save_fname = ['DPSD_' filename(1:end-4)];
% else
%     save_fname = ['DPSD_emulator-K' num2str(K)];
% end
% 
% save(['~/Documents/code/DPSD/dpsd_outputs/' iq_source '/' signal_type '/x100/' save_fname '.mat'],...
%     'sZDR', 'sPHV', 'sSH', 'sSV', 'sSX', 'iqh', 'iqv', 'vvx', 'd', 'params');

