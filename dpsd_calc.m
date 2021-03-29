function [sZDR, sPHV] = dpsd_calc(VH, VV, d, M, K)

% Coherency correction factor
CX_L = 0.5 * (VH(1)./VH(M) + VV(1)./VV(M)) * ones(M-1,1);
CX_R = 0.5 * (VH(M)./VH(1) + VV(M)./VV(1)) * ones(M-1,1);

VH_L = CX_L .* VH(1:M-1);
VH_R = CX_R .* VH(2:M);
VV_L = CX_L .* VV(1:M-1);
VV_R = CX_R .* VV(2:M);

XH = [VH_L; VH; VH_R];
XV = [VV_L; VV; VV_R];

% Maximum ratio of corrected samples
a = mean((d/max(d)).^2);
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

block_start_inds = find(r <= rmax);
j = randsample(block_start_inds, K, true);


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
    
%     ZH(:,n,:) = fft(d_rep .* BH, M, 1);
%     ZV(:,n,:) = fft(d_rep .* BV, M, 1);
    
    % Compute PSDs
    for k = 1:M
        ZH(k,n,:) = sum( repmat(d, [1 rlzs]) .* BH .* exp(-1j*2*pi*(k-1)*repmat(ks(:),[1 rlzs])/M) , 1 );
        ZV(k,n,:) = sum( repmat(d, [1 rlzs]) .* BV .* exp(-1j*2*pi*(k-1)*repmat(ks(:),[1 rlzs])/M) , 1 );
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
    b = (1-rmax)^-3.3 - 2*(1-rmax)^1.1;
elseif K > 1
    b = (1-rmax)^-4.5 - (1-rmax)^-2.1;
end


end