function [V_h, V_v, ZDR, PHV] = iq_emulate(M, S_dB, vr, sw, SNR, va, rlzs, zdr, phv)
%
% %-----Uncomment for debugging------
% M = 100;
% S_dB = 45;
% vr = -30;
% sw = 5;
% SNR = 30;
% va = 100;
% rlzs = 1;
% zdr = 2; % dB
% phv = 0.70;
% ------------------------

%% -------- METHOD 1 -------- %%

% x = ts_sp_dftm(M, 0, vr, sw, SNR, va, rlzs, [], 'double');
% y = ts_sp_dftm(M, 0, vr, sw, SNR, va, rlzs, [], 'double');
%
% SHH = 10^(S_dB/10);
% SVV = 10^(-zdr/10) * SHH;
%
% V_h = sqrt(SHH) * x;
% V_v = sqrt(SVV) * (phv * x + sqrt(1 - phv^2) * y);
% %V_v = sqrt(SVV / mean(V_v .* conj(V_v))) * V_v;


%% -------- METHOD 2 -------- %%

V_h = zeros(M, rlzs);
V_v = zeros(M, rlzs);

SHH = 10^(S_dB/10);
SVV = 10^(-zdr/10) * SHH; % this is the power we want for the V channel based on ZDR

n = 2*M-1; % two-sided power spectrum
vv = va * linspace(-1, 1, n);

for a = 1:rlzs
    S_k(1,:) = SHH / sqrt(2*pi*sw^2) .* exp(-(vv - vr).^2 / (2*sw^2));
    S_k(2,:) = SVV / sqrt(2*pi*sw^2) .* exp(-(vv - vr).^2 / (2*sw^2));
    X_k = rand(2,n);
    y_k = 2*pi*rand(4,n) - pi; % uniform random phase
    
    
    SNR_linear = 10^(SNR/10);
    N = 1/n/SNR_linear * sum(S_k,2);
    
    P_k = -(S_k + N*ones(1,n)) .* log(X_k); % total power spectrum
    
    V1 = zeros(M,1);
    V2 = V1;
    % The signal and noise should each have random phases
    for m = 1:M % generate uncorrelated H/V time series
        % V1(m) = 1/n * sum( sqrt(P_k(1,:)) .* exp(1j*y_k(1,:)) .* exp(-1j*2*pi*(-(M-1):M-1)*m/n) );
        V1(m) = 1/n * sum( sqrt(N(1)) * exp(1j*y_k(1,:)) .* exp(-1j*2*pi*(-(M-1):M-1)*m/n) + ...
            sqrt(S_k(1,:)) .* exp(1j*y_k(3,:)) .* exp(-1j*2*pi*(-(M-1):M-1)*m/n) );
        % V2(m) = 1/n * sum( sqrt(P_k(2,:)) .* exp(1j*y_k(2,:)) .* exp(-1j*2*pi*(-(M-1):M-1)*m/n) );
        V2(m) = 1/n * sum(sqrt(N(2)) * exp(1j*y_k(2,:)) .* exp(-1j*2*pi*(-(M-1):M-1)*m/n) + ...
            sqrt(S_k(2,:)) .* exp(1j*y_k(4,:)) .* exp(-1j*2*pi*(-(M-1):M-1)*m/n) );
    end
    
%     rh = abs(mean(V1 .* conj(V2))) ./ sqrt(mean(V1 .* conj(V1)) .* mean(V2 .* conj(V2)));
%     fprintf(['rh: ', num2str(rh), '\n']);
    if a == 1
        figure(175)
        plot(1:M, real(V1), '-k', 1:M, imag(V1), '--k',...
            1:M, real(V2), '-b', 1:M, imag(V2), '--b')
        grid on
    end
    
    
    %------Implement time series partial correlation dependence???-----------
    
    s1 = sqrt(1 / mean(V1.*conj(V1))) * V1; % unit-power
    s2 = sqrt(1 / mean(V2.*conj(V2))) * V2;
    
    x1 = sqrt(SHH) * s1;
    % Calculate x2 to be partially correlated with x1 with desired rhoHV
    x2 = sqrt(SVV) * (phv * s1 + sqrt(1 - phv^2) * s2);
    
%     zdr_tmp = 10 * log10(mean(x1.*conj(x1)) ./ mean(x2.*conj(x2)));
%     fprintf(['zdr: ', num2str(zdr_tmp), '\n']);
    
    V_h(:,a) = x1(:);
    V_v(:,a) = x2(:);
    
    if a == 1
        figure(175)
        plot(1:M, real(V_h(:,a)), '-k', 1:M, imag(V_h(:,a)), '--k',...
            1:M, real(V_v(:,a)), '-b', 1:M, imag(V_v(:,a)), '--b')
        grid on
    end
    
end

%% -------- Bulk parameters -------- %%
% Don't comment this out!

SHH = mean(V_h .* conj(V_h));
SVV = mean(V_v .* conj(V_v));
SXX = mean(V_h .* conj(V_v));
ZDR = 10 * log10(SHH ./ SVV);
PHV = abs(SXX) ./ sqrt(SHH .* SVV);



