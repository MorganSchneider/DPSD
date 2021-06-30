function [szdr, sphv, svar, spsd] = dpsd_calc(iqh, iqv, d, M, K, rlzs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS
% 
%   iqh:    H channel I/Q data, size(iqh,1) = M
%   iqv:    V channel I/Q data, size(iqv,1) = M
%   d:      Windowing function, size(d) = [M,1]
%   M:      Number of samples per gate
%   K:      Number of pseudorealizations to generate
%   rlzs:   Number of independent realizations per gate (rlzs = 1 if using
%           SimRadar data)
%
% OUTPUTS
%   szdr:   Spectral ZDR (dB), bias-corrected and averaged over rlzs
%   sphv:   Spectral rhoHV, bias-corrected and averaged over rlzs
%   svar
%       svar.zdr:   Moving-block variance of szdr (dB)
%       svar.phv:   Moving block variance of sphv
%   spsd
%       spsd.h:     H channel PSD estimate, sS_H
%       spsd.v:     V channel PSD estimate, sS_V
%       spsd.x:     Cross-channel PSD estimate, sS_X
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nvar = 9; % number of points for moving variance

% Choose data window
%d = hann(M); % hamming(M), hann(M), blackman(M), rectwin(M)
R0_d = mean(abs(d).^2); % window power

dm2 = size(iqh,2);
dm3 = size(iqh,3);

% Data window
d_rep = repmat(d, [1 dm2 dm3]);
mm = repmat((0:M-1)', [1 dm2 dm3]);

Rxx_mat = zeros(2*M-1, dm2, dm3);
for m = 1:dm2
    for n = 1:dm3
        [Rxx_mat(:,m,n), ~] = xcorr(squeeze(iqh(:,m,n))); % ACF
    end
end
Sxx_mat = fftshift(fft(Rxx_mat, 2*M-1, 1), 1);


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
sSH.n = squeeze(mean(ZH.*conj(ZH) / (M*R0_d), 2));
sSV.n = squeeze(mean(ZV.*conj(ZV) / (M*R0_d), 2));
sSX.n = squeeze(mean(ZH.*conj(ZV) / (M*R0_d), 2));

% Compute DPSD estimates
sZDR.n = sSH.n ./ sSV.n;
sPHV.n = abs(sSX.n) ./ sqrt(sSH.n .* sSV.n);

% Bias correction
if rlzs == 1
    b = (1-rmax)^-3.3 - 2*(1-rmax)^1.1;
elseif rlzs > 1
    b = (1-rmax)^-4.5 - (1-rmax)^-2.1;
end
sZDR.c = sZDR.n .* (1 - ((1/b/rlzs) * (1 - (sPHV.n).^2)));
sPHV.c = sPHV.n .* (1 - ((1/b/rlzs) * (1 - (sPHV.n).^2).^2 ./ (4*(sPHV.n).^2)));
sPHV.c(sPHV.c < 0) = 0;
sZDR.c(sZDR.c < 0) = real(sZDR.c(sZDR.c < 0));


% Average over independent rlzs
if dm3 == 1
    sSH.f = squeeze(mean(sSH.n, 2));
    sSV.f = squeeze(mean(sSV.n, 2));
    sSX.f = squeeze(mean(sSX.n, 2));
    sZDR.f = squeeze(mean(sZDR.n, 2));
    sPHV.f = squeeze(mean(sPHV.n, 2));
    sZDR.cf = squeeze(mean(sZDR.c, 2));
    sPHV.cf = squeeze(mean(sPHV.c, 2));
else
    sSH.f = sSH.n;
    sSV.f = sSV.n;
    sSX.f = sSX.n;
    sZDR.f = sZDR.n;
    sPHV.f = sPHV.n;
    sZDR.cf = sZDR.c;
    sPHV.cf = sPHV.c;
end

% Convert sZDR to logarithmic units (dB)
sZDR.n = 10*log10(sZDR.n);
sZDR.f = 10*log10(sZDR.f);
sZDR.c = 10*log10(sZDR.c);
sZDR.cf = 10*log10(sZDR.cf);


% Set output variables
szdr = sZDR.cf;
sphv = sPHV.cf;

spsd.h = sSH.f;
spsd.v = sSV.f;
spsd.x = sSX.f;

svar.phv = movvar(sPHV.cf,nvar,0,1,'omitnan','Endpoints','shrink');
svar.zdr = movvar(sZDR.cf,nvar,0,1,'omitnan','Endpoints','shrink');


end