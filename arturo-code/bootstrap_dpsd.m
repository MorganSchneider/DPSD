function E = bootstrap_dpsd(V, w, N0, NFT, B, K, N)
% BOOTSTRAP DUAL-POLARIMETRIC DENSITY ESTIMATOR FOR POLARIMETRIC WEATHER
% RADAR I/Q SIGNALS
%
% Usage: bootstrap_dpsd(V, w, N0, NFT, B, K, N)
%
% Inputs:
%   V       : I/Q signal (structure) 
%   V.H     : H-channel I/Q signal (vector/matrix)
%             If matrix, each row corresponds to a different data.
%   V.V     : V-channel I/Q signal (vector/matrix)
%             If matrix, each row corresponds to a different data.
%   w       : Data windowing function (vector)
%   N0      : Average noise power (structure)
%   N0.H    : H-channel average noise power (scalar)
%   N0.V    : H-channel average noise power (scalar)
%   NFT     : Number of DFT points (scalar)
%   B (K')  : Number of pseudo-realizations (scalar)
%   K       : Number of spectra to average for DPSD estimation (scalar)
%             If K > 1, spectra will be averaged for DPSD estimation.
%             An example of signal arrangement is:
%               -> <a,1>
%               -> <a,2>  This will produce DPSDs by averaging PSDs <a1,a2>, <b1,b2>,
%               -> <b,1>  and <c1,c2>, and so on.
%               -> <b,2>  And the size of the I/Q matrix is (N*K)xM
%               -> <c,1>  
%               -> <c,2>
%   N       : Number of 'independent' range gates in the matrix
%
% Outputs:
%   E       : Estimates (structure)
%   E.sS    : Spectral power (structure)
%   E.sS.H  : H-spectral power (NxNFT matrix)
%   E.sS.V  : V-spectral power (NxNFT matrix)
%   E.sS.X  : Cross-spectral power (NxNFT matrix)
%   E.sSNR  : Spectral SNR (structure)
%   E.sSNR.H: H-spectral SNR (NxNFT matrix)
%   E.sSNR.V: V-spectral SNR (NxNFT matrix)
%   E.sD    : Spectral ZDR (NxNFT matrix)
%   E.sR    : Spectral RHV (NxNFT matrix)
%
% % % %
% Author: Arturo Umeyama (c) 2016
% % % %

% Compute the maximum ratio
r = 0.5 - sqrt(mean(w.^2))*0.5;

NK = size(V.H, 1); % Number of rows in I/Q - each row can be an independent
                  % signal, a different range gate, etc.
M = size(V.H, 2); % Number of samples in I/Q 

if isempty(NFT), NFT = M; end;

% Set variables to store the PSDs
S.H = nan(NK, NFT);
S.V = nan(NK, NFT);
S.X = nan(NK, NFT);

% Do this for each row of the I/Q matrix 
for i = 1:NK
    % Set temporary variable
    tV.H = V.H(i,:);
    tV.V = V.V(i,:);
    
    % Compute combined correction factors
    CX_left.H = 0.5*(tV.H(1)/tV.H(end) + tV.V(1)/tV.V(end));
    CX_left.V = 0.5*(tV.H(1)/tV.H(end) + tV.V(1)/tV.V(end));
    CX_right.H = 0.5*(tV.H(end)/tV.H(1) + tV.V(end)/tV.V(1));
    CX_right.V = 0.5*(tV.H(end)/tV.H(1) + tV.V(end)/tV.V(1));
    
    % Set maximum ratio, apply correction factors, and drop first/last
    % samples
    VL.H = tV.H((end - round(M * r) + 1):(end - 1)) * CX_left.H;
    VL.V = tV.V((end - round(M * r) + 1):(end - 1)) * CX_left.V;
    VR.H = tV.H(2:round(M * r)) * CX_right.H;
    VR.V = tV.V(2:round(M * r)) * CX_right.V;
    
    % Get extended signal
    X.H = [VL.H, tV.H, VR.H];
    X.V = [VL.V, tV.V, VR.V];
    
    % Get size of extended signal
    Mx = size(X.H, 2);
    % Bootstrapping: generate the starting indexes of the B (K')
    % pseudo-realizations
    boot_indexes = randi(Mx - M + 1, [B, 1]);
    % Then, complete the index sequence with M samples
    boot_indexes = bsxfun(@plus, boot_indexes - 1, repmat(1:M, [B, 1]));
    
    % Save some memory, re-use tV 
    clear tV;
    % Retrieve the I/Q data of the bootstrapped indexes. These will be the
    % bootstrapped pseudo-realizations
    tV.H = X.H(boot_indexes);
    tV.V = X.V(boot_indexes);
    
    % Compute original signal power for power correction
    R0.H = mean(bsxfun(@times, V.H(i,:), conj(V.H(i,:))), 2);
    R0.V = mean(bsxfun(@times, V.V(i,:), conj(V.V(i,:))), 2);
    
    % Compute the signal power for each pseudo-realization
    tR0.H = mean(bsxfun(@times, tV.H, conj(tV.H)), 2);
    tR0.V = mean(bsxfun(@times, tV.V, conj(tV.V)), 2);
    
    % Apply scaling to correct the power of pseudo-realizations
    tV.H = bsxfun(@times, sqrt(bsxfun(@rdivide, R0.H, tR0.H)), tV.H);
    tV.V = bsxfun(@times, sqrt(bsxfun(@rdivide, R0.V, tR0.V)), tV.V);
    
    % Apply windowing function then the Fourier transform
    z.H = fftshift(fft(bsxfun(@times, tV.H, w), NFT, 2), 2);
    z.V = fftshift(fft(bsxfun(@times, tV.V, w), NFT, 2), 2);
    
    % Compute the PSDs the mean of the PSDs of each pseudo-realization.
    % Also account for window power here.
    alpha = mean(abs(w).^2);
    S.H(i,:) = mean(bsxfun(@power, abs(z.H), 2) / (M * alpha), 1);
    S.V(i,:) = mean(bsxfun(@power, abs(z.V), 2) / (M * alpha), 1);
    S.X(i,:) = mean(bsxfun(@times, z.H, conj(z.V)) / (M * alpha), 1);
end

tsh = nan(N, M);
tsv = nan(N, M);
tsx = nan(N, M);
td = nan(N, M);
tr = nan(N, M);

for i = 1:N
    % If averaging spectra from different sources for DPSD estimation,
    % arrange the signals consecutively, i.e.:
    % -> <a,1>
    % -> <a,2>  This will produce DPSDs by averaging PSDs <a1,a2>, <b1,b2>,
    % -> <b,1>  and <c1,c2>, and so on.
    % -> <b,2>
    % -> <c,1>
    % -> <c,2>
    %     ...
    iK = (1:K) + (i-1)*K;
    
    % Averages the PSDs from different sources for DPSD estimation
    tsh(i,:) = mean(S.H(iK,:), 1);
    tsv(i,:) = mean(S.V(iK,:), 1);
    tsx(i,:) = mean(S.X(iK,:), 1);
    
    % Estimates the DPSDs (before correction)
    td(i,:) = bsxfun(@rdivide, tsh(i,:), tsv(i,:));
    tr(i,:) = bsxfun(@rdivide, abs(tsx(i,:)), sqrt( bsxfun(@times, tsh(i,:), tsv(i,:)) ) );
end

if K == 1
    beta = (1-r)^-3.3 - 2*(1-r)^1.1;
else
    beta = (1-r)^-4.5 - (1-r)^-2.1;
end

% Stores spectral variables in structure E
E.sS.H = tsh;
E.sS.V = tsv;
E.sS.X = tsx;
E.sSNR.H = bsxfun(@rdivide, tsh, N0.H);
E.sSNR.V = bsxfun(@rdivide, tsv, N0.V);

% Apply correction to DPSDs
E.sD = bsxfun(@times, td, (1 - 1 / (beta * K) * (1 - bsxfun(@power, tr, 2))) );
E.sR = bsxfun(@times, tr, (1 - 1 / (beta * K) * (bsxfun(@rdivide, bsxfun(@power, (1 - bsxfun(@power, tr, 2)), 2), 4 * bsxfun(@power, tr, 2)))));

% After the correction, some values can be off the range for low SNR
% coefficients 
E.sD(E.sD < 0) = nan;
E.sR(E.sR < 0) = 0;