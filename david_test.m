% This code computes DPSDs using the Bootstrap technique discussed in 
% Umeyama et al. (2017). It can generate arbitrary DPSDs and spectra and 
% also read in SimRadar dual-pol IQ. Traditional DPSD methods, such as Yu et 
% al. (2012)'s averaging technique are also computed.

clear
close all

david_file = true;
force_debris = true;
multi_type = 2; % 1: Umeyama case, 2: Debris
M = 128 - 1; % Number of Samples
choose_gate = true;

% Range gate choice for SimRadar if choose_gate = true
% rind = 11;
% azind = 15;
rind = 15;
azind = 31;

% Thresholds for force_debris
rho_thres = 0.8;
z_thres = 40;

num_bootstrap = 50; % Number of bootstrapped psuedorealizations
rlzs = 1; % Number of IQ samples

% Perform range averaging
range_avg = true;
azi_avg = true;
avg_int_range = 3;
avg_int_azi = 3;

plot_save_flag = 0; % set plot_save_flag = 1 to save plots from this script

% K = 20; % number of bootstrapped pseudorealizations

[iq_source, signal_type] = choose_source();

zmap = feval('boonlib', 'zmap', 21);
rgmap = feval('boonlib', 'rgmap', 21);

if strcmp(iq_source, 'emulator_yu')
%     va = 50;

    if multi_type == 1 && strcmp(signal_type, 'multi')
        va = 15;
    else
        va = 50;
    end
    vv = linspace(-va,va,M);
    vv2 = vv;
    S1 = 10 ^ (30/10);
    S2 = 10 ^ (50/10);
    sw1 = 10;
    sw2 = 5;
    vr1 = -25;
    vr2 = 25;

    lambda = 0.1; % wavelength
    Ts = lambda/4/va;
    if strcmp(signal_type, 'rain')
        rhv = [0.98 0.98];
        zdr = [3 3];
        Sv = S1 / sqrt(2*pi*sw1^2) .* exp(-(vv - vr1).^2 / (2*sw1^2));
        sphv = rhv(1) * ones(size(Sv));
        szdr = zdr(1) * ones(size(Sv));
    elseif strcmp(signal_type, 'debris')
        rhv = [0.5 0.5];
        zdr = [0 0];
        Sv = S2 / sqrt(2*pi*sw2^2) .* exp(-(vv - vr2).^2 / (2*sw2^2));
        sphv = rhv(2) * ones(size(Sv));
        szdr = zdr(2) * ones(size(Sv));
    elseif strcmp(signal_type, 'multi')
        if multi_type == 1
            % Overwrite
            rhv = [0.99 0.85];
            zdr = [0.5 2];
            S1 = 10 ^ (20/10);
            S2 = 10 ^ (35/10);
            vr1 = -3;
            vr2 = 3;
            sw1 = 2.5;
            sw2 = 0.75;
        else
            % Sets DPSD rhohv and ZDR values to the values below. The
            % separation between the two occurs at the midpoint between the
            % two radial velocities
            rhv = [0.5 0.98];
            zdr = [-2 5];
        end
 
        % Create two Gaussians
        Sv1 = S1 / sqrt(2*pi*sw1^2) .* exp(-(vv2 - vr1).^2 / (2*sw1^2));
        Sv2 = S2 / sqrt(2*pi*sw2^2) .* exp(-(vv2 - vr2).^2 / (2*sw2^2));
        Sv = Sv1 + Sv2;
%         
%       % Set DPSDs (can modify this to any arbitrary function
        [c,ind] = min(abs(vv-(vr1+vr2)/2));
        sphv(1:ind) = rhv(1) * linspace(1,0,ind);
        sphv(ind+1:max(size(vv))) = rhv(2) * linspace(0,1,ind-1);
        szdr(1:ind) = zdr(1);
        szdr(ind+1:max(size(vv))) = zdr(2);
        
        % Calculate corresponding velocities
        prf = 4 * va / lambda;
        fs = prf / 2; 
        f0 = fs / M;
        vv = -f0 * lambda * (0:M-1);
        vv(vv > va) = vv(vv > va) - 2*va;
        vv(vv < -va) = vv(vv < -va) + 2*va;
        vv = fftshift(vv);
    end
    
    tmp_sphv = sphv;
    tmp_szdr = szdr;
    
    figure(15)
    
    subplot(2,2,1)
    plot(vv, Sv)
    grid on
    
    subplot(2,2,2)
    plot(vv, tmp_sphv)
    grid on
    set(gca, 'ylim', [0 1])
    
    subplot(2,2,3)
    plot(vv, tmp_szdr)
    grid on
    szdr_linear = 10 .^ (szdr./10);
    set(gca, 'ylim', [nanmin(zdr)-1 nanmax(zdr)+1])
    [V_h, V_v, ZDR_bulk, PHV_bulk] = iq_emulate_yu(Sv, sphv, szdr_linear, rlzs);
    % Time series output needs to be of format (N_samples, N_realizations)
    
    iqh=V_h;
    iqv=V_v;
    szdr_save=szdr;
    sphv_save=sphv; 
    clear sphv szdr; % Clear these for later    
    
%---------If using SimRadar I/Q data---------%
elseif strcmp(iq_source, 'simradar')
    
    if strcmp(signal_type, 'rain')
        filename = 'suctvort_rain'; %Same simulation as deb4_10000 except no debris
        load([filename '.mat']);
        force_debris = 0; %Turn force debris off if selected
    elseif strcmp(signal_type, 'debris')
        filename = 'sim-PPI0.5-TCU-d3n100000';
        load([filename '.mat']);
    elseif strcmp(signal_type, 'multi')
        if david_file
            filename = 'deb4_10000';
        else
            filename = 'sim-PPI0.5-DCU-d3n100000';
        end
        load([filename '.mat']);
    end
    
    M = size(iqh, 3);
    lambda = dat.params.lambda;
    Ts = dat.params.prt;
    prf = 1/Ts;
    va = dat.params.va;
    rlzs = 1;
    %Calculate corresponding velocities
    % prf = 4 * va / lambda;
    fs = prf / 2; 
    f0 = fs / M;

%         fs = 2*va/lambda;
%         f0 = fs/M;
%         f = f0 * inds_2s;
    vv = -f0 * lambda * (0:M-1);
    vv(vv > va) = vv(vv > va) - 2*va;
    vv(vv < -va) = vv(vv < -va) + 2*va;
    vv = fftshift(vv);
    
    iqh_full = iqh;
    iqv_full = iqv;
    rhohv_test = 1;
    if choose_gate
         r_ind = rind;
         az_ind = azind;
    else
        % Force the choice of a gate with low rhohv and high Zh
        if force_debris
            while rhohv_test > rho_thres || r_ind == 1 || z_test < z_thres
                r_ind = randsample(size(iqh,1), 1); % randomly pull from SimRadar output
                az_ind = randsample(size(iqh,2), 1);
                rhohv_test = rhohv(r_ind, az_ind);
                z_test = zh(r_ind, az_ind);
            end
        else
             r_ind = randsample(size(iqh,1), 1); % randomly pull from SimRadar output
             az_ind = randsample(size(iqh,2), 1);
        end
    end
    
    % Range averaging option to use multiple gates for additional
    % realizations
    if range_avg && azi_avg
        half_dr = floor(avg_int_range / 2);
        half_az = floor(avg_int_azi / 2);
        iqh = double(squeeze(iqh(r_ind-half_dr:r_ind+half_dr, az_ind-half_az:az_ind+half_az, :)));
        iqv = double(squeeze(iqv(r_ind-half_dr:r_ind+half_dr, az_ind-half_az:az_ind+half_az, :)));
        iqh_tmp = permute(iqh, [3 1 2]);
        iqv_tmp = permute(iqv, [3 1 2]);
        clear iqh iqv
        for idx = 1:size(iqh_tmp,2)
            for jdx = 1:size(iqh_tmp,3)
                iqh(:,(idx-1)*size(iqh_tmp,2)+jdx) = squeeze(iqh_tmp(:,idx,jdx));
                iqv(:,(idx-1)*size(iqv_tmp,2)+jdx) = squeeze(iqv_tmp(:,idx,jdx));
            end
        end
        rlzs = avg_int_range * avg_int_azi;
    elseif range_avg
        half_dr = floor(avg_int_range / 2);
        iqh = double(squeeze(iqh(r_ind-half_dr:r_ind+half_dr, az_ind, :)))';
        iqv = double(squeeze(iqv(r_ind-half_dr:r_ind+half_dr, az_ind, :)))';
        rlzs = avg_int_range;
    elseif azi_avg
        half_az = floor(avg_int_azi / 2);
        iqh = double(squeeze(iqh(r_ind, az_ind-half_az:az_ind+half_az, :)))';
        iqv = double(squeeze(iqv(r_ind, az_ind-half_az:az_ind+half_az, :)))';
        rlzs = avg_int_azi;
    else
        iqh = double(squeeze(iqh(r_ind, az_ind, :)));
        iqv = double(squeeze(iqv(r_ind, az_ind, :)));
    end
    fprintf(['rhohv=', num2str(rhohv(r_ind, az_ind)), '\n']);
    fprintf(['ZDR=', num2str(zdr(r_ind, az_ind)), '\n']);
    
    figure(101)
    
    ha = subplot(2,2,1);
    hs = pcolor(xx, yy, zh(:,:,1));
    hold on
        plot(xx(r_ind, az_ind), yy(r_ind, az_ind), 'sk', 'MarkerSize', 11, 'MarkerFaceColor', 'w', 'LineWidth', 2);
    hold off
    set(gca, 'DataAspect', [1 1 1])
    caxis([0 80])
    colormap(ha, zmap)
    colorbar
    shading flat
    set(gca, 'YDir', 'Normal')
    title('Reflectivity (dBZ)')
    
    ha(2) = subplot(2,2,2);
    hs(2) = pcolor(xx, yy, vr(:,:,1));
    hold on
        plot(xx(r_ind, az_ind), yy(r_ind, az_ind), 'sk', 'MarkerSize', 11, 'MarkerFaceColor', 'w', 'LineWidth', 2);
    hold off
    set(gca, 'DataAspect', [1 1 1])
    caxis([-1 1] * round(max(max(abs(vr(:,:,1)))), -1))
    colormap(ha(2), rgmap)
    colorbar
    shading flat
    title('Velocity (m/s)')
    
    ha(3) = subplot(2,2,3);
    hs(3) = pcolor(xx, yy, rhohv(:,:,1));
    hold on
        plot(xx(r_ind, az_ind), yy(r_ind, az_ind), 'sk', 'MarkerSize', 11, 'MarkerFaceColor', 'w', 'LineWidth', 2);
    hold off
    set(gca, 'DataAspect', [1 1 1])
    caxis([0.2 1])
    colormap(ha(3), zmap)
    colorbar
    shading flat
    title('Rhohv')
    
    axes('Unit', 'Normalized', 'Position', [0.5 0.90 0.01 0.01])
    title(['SimRadar-generated spectrum (r=', num2str(round(r(r_ind)/1000,2)),...
        'km, \phi_{az}=' num2str(round(180/pi*az_rad(az_ind),1)) '^o)'], 'FontSize', 12)
    axis off
end

% Choose data window
% d = hamming(M);
% d = blackman(M);
d = hann(M);
% d = rectwin(M);

if rlzs > 1
    d_rep = repmat(d, [1 size(iqh,2)]);
    mm = repmat((0:M-1)', [1 size(iqh,2)]);

    % Calculate DPSD the Original Way prior to bootstrap
    DFT_H = fftshift(fft(d_rep.*iqh, M, 1));
    DFT_V = fftshift(fft(d_rep.*iqv, M, 1));
    Ph = nanmean(iqh.*conj(iqh), 1);
    Pv = nanmean(iqv.*conj(iqv), 1);

    % Periodogram estimate using Yu et al. (2012)
    SHF_trad = nanmean(abs(DFT_H) .^ 2/M, 2);
    SVF_trad = nanmean(abs(DFT_V) .^ 2/M, 2);
    SXF_trad = nanmean(DFT_H.*conj(DFT_V) / M, 2);
else
%     d_rep = repmat(d, [1 size(iqh,2)]);
%     mm = repmat([0:M-1]', [1 size(iqh,2)]);

    % Calculate DPSD the Original Way prior to bootstrap
    DFT_H = fftshift(fft(d.*iqh, M, 1), 1);
    DFT_V = fftshift(fft(d.*iqv, M, 1), 1);
    Ph = nanmean(iqh.*conj(iqh), 1);
    Pv = nanmean(iqv.*conj(iqv), 1);

    % Periodogram estimate using Yu et al. (2012)
    SHF_trad = abs(DFT_H) .^ 2/M;
    SVF_trad = abs(DFT_V) .^ 2/M;
    SXF_trad = DFT_H .* conj(DFT_V) / M;
end

figure(152)

subplot(2,1,1)
yyaxis left
plot(vv, 10*log10(SHF_trad), '-k', vv, 10*log10(SVF_trad), '-b')
xlabel('V_{r} (m s^{-1})', 'FontSize', 14)
ylabel('S_{H,V}', 'FontSize', 14)
yyaxis right
plot(vv, 10*log10(SHF_trad ./ SVF_trad), '--k')
ylabel('S(Z_{DR}) dB', 'FontSize', 14)

subplot(2,1,2)
yyaxis left
plot(vv, 10*log10(SHF_trad), '-k', vv, 10*log10(SVF_trad), '-b')
xlabel('V_{r} (m s^{-1})', 'FontSize', 14)
yyaxis right
plot(vv, abs(SXF_trad) ./ sqrt(SHF_trad.*SVF_trad), '--k')
ylabel('S(\rho_{HV})', 'FontSize', 14)

% Sv_Save = Sv;
% Sv1_save = Sv1;
% Sv2_save = Sv2;
clear Sv Sv1 Sv2

% Coherency and extended IQ construction
cxp = 1/2 * (iqh(M,:) ./ iqh(1,:) + iqv(M,:) ./ iqv(1,:));
cxm = 1/2 * (iqh(1,:) ./ iqh(M,:) + iqv(1,:) ./ iqv(M,:));
iqlh = cxm .* iqh(1:M-1,:);
iqlv = cxm .* iqv(1:M-1,:);
iqrh = cxp .* iqh(2:M,:);
iqrv = cxp .* iqv(2:M,:);

Xh = [iqlh; iqh; iqrh];
Xv = [iqlv; iqv; iqrv];

% Compute maximum ratio
rmax = (1 - sqrt(1/M * sum(abs(d/max(d)).^2))) / 2;
corr = [ones(size(iqlh)); zeros(size(iqh)); ones(size(iqrh))];

for idx = 1:max(size(corr)) - M + 1
    r(idx) = sum(corr(idx:idx+M-1)) / M;
end

% Find blocks meeting criteria for enough original data points
valid_blocks = find(r < rmax);

rsamp = randsample(valid_blocks, num_bootstrap, 'true');

% Grab random sample and normalized by original signal power
if rlzs > 1
    Bh = zeros(num_bootstrap, M, rlzs); 
    Bv = Bh;
    for jdx = 1:num_bootstrap
        Bh(jdx,:,:) = Xh(rsamp(jdx): rsamp(jdx)+M-1, :);
        Phs = squeeze(nanmean(Bh(jdx,:,:) .* conj(Bh(jdx,:,:))));
        Bh(jdx,:,:) = sqrt(repmat(Ph,[M 1]) ./ repmat(Phs',[M 1])) .* squeeze(Bh(jdx,:,:));
        Bv(jdx,:,:) = Xv(rsamp(jdx): rsamp(jdx)+M-1, :);
        Pvs = squeeze(nanmean(Bv(jdx,:,:) .* conj(Bv(jdx,:,:))));
        Bv(jdx,:,:) = sqrt(repmat(Pv,[M 1]) ./ repmat(Pvs',[M 1])) .* squeeze(Bv(jdx,:,:));
    end
    d_rep = repmat(d', [num_bootstrap 1 rlzs]);
    
    % Calculate DPSD the Original Way prior to bootstrap
    DFT_H_boot = fftshift(fft(d_rep.*Bh, M, 2));
    DFT_V_boot = fftshift(fft(d_rep.*Bv, M, 2));
    SHF_boot = nanmean(nanmean(abs(DFT_H_boot) .^ 2/M, 1), 3);
    SVF_boot = nanmean(nanmean(abs(DFT_V_boot) .^ 2/M, 1), 3);
    SXF_boot = nanmean(nanmean(DFT_H_boot.*conj(DFT_V_boot) / M, 1), 3);
else
    Bh = zeros(num_bootstrap, M); 
    Bv = Bh;
    for jdx = 1:num_bootstrap
        Bh(jdx,:) = Xh(rsamp(jdx): rsamp(jdx)+M-1);
        Phs = squeeze(nanmean(Bh(jdx,:) .* conj(Bh(jdx,:))));
        Bh(jdx,:) = sqrt(repmat(Ph,[M 1]) ./ repmat(Phs',[M 1])) .* squeeze(Bh(jdx,:))';
        Bv(jdx,:) = Xv(rsamp(jdx): rsamp(jdx)+M-1);
        Pvs = squeeze(nanmean(Bv(jdx,:) .* conj(Bv(jdx,:))));
        Bv(jdx,:) = sqrt(repmat(Pv,[M 1]) ./ repmat(Pvs',[M 1])) .* squeeze(Bv(jdx,:))';
    end
    d_rep = repmat(d', [num_bootstrap 1]);
    % Calculate Bootstrap DPSD
    DFT_H_boot = fftshift(fft(d_rep.*Bh, M, 2), 2);
    DFT_V_boot = fftshift(fft(d_rep.*Bv, M, 2), 2);
    SHF_boot = nanmean(abs(DFT_H_boot) .^ 2/M, 1);
    SVF_boot = nanmean(abs(DFT_V_boot) .^ 2/M, 1);
    SXF_boot = nanmean(DFT_H_boot.*conj(DFT_V_boot) / M, 1);
    
    vv = -vv; % Temporary fix since I can't figure out what's wrong (flipped sign somewhere)
end

figure(175)

subplot(2,1,1)
yyaxis left
plot(vv, 10*log10(SHF_boot), '-k', vv, 10*log10(SVF_boot), '-b')
xlabel('V_{r} (m s^{-1})', 'FontSize', 14)
ylabel('S_{H,V}', 'FontSize', 14)
yyaxis right
plot(vv, 10*log10(SHF_boot ./ SVF_boot), '--k')
ylabel('S(Z_{DR}) dB', 'FontSize', 14)

subplot(2,1,2)
yyaxis left
plot(vv, 10*log10(SHF_boot), '-k', vv, 10*log10(SVF_boot), '-b')
xlabel('V_{r} (m s^{-1})', 'FontSize', 14)
yyaxis right
plot(vv, abs(SXF_boot) ./ sqrt(SHF_boot .* SVF_boot), '--k')
ylabel('S(\rho_{HV})', 'FontSize', 14)


