clear
close all

M = 64;
va = 50;
vv = linspace(-va, va, M);
S1 = 10 ^ (30/10);
S2 = 10 ^ (50/10);
K = 100; % Number of spectra
sw1 = 5;
sw2 = 5;
vr1 = -10;
vr2 = 15;
rhv = [0.5 0.98];
zdr = [0 5];
snr = 30;

sim_type = input('Enter simulation type (1: Single peak, 2: Bimodal): ');

Sv1 = S1 / sqrt(2*pi*sw1^2) .* exp(-(vv - vr1).^2 / (2*sw1^2));

if sim_type == 2
    Sv2 = S2 / sqrt(2*pi*sw2^2) .* exp(-(vv - vr2).^2 / (2*sw2^2));
    Sv = Sv1 + Sv2;
    [c, ind] = min(abs(vv - (vr1+vr2) / 2));
    sphv(1:ind) = rhv(1);
    sphv(ind+1: max(size(vv))) = rhv(2);
    szdr(1:ind) = zdr(1);
    szdr(ind+1: max(size(vv)))= zdr(2);
else
    Sv = Sv1;
    sphv = rhv(1) * ones(size(Sv));
    szdr = zdr(1) * ones(size(Sv));
end

figure(1)

subplot(2,2,1)
plot(vv, Sv)
grid on

subplot(2,2,2)
plot(vv, sphv)
grid on
set(gca, 'ylim', [0 1])

subplot(2,2,3)
plot(vv, szdr)
grid on
szdr_linear = 10 .^ (szdr./10);
set(gca, 'ylim', [nanmin(zdr)-1 nanmax(zdr)+1])
[V_h, V_v, ZDR, PHV] = iq_emulate_yu(Sv, sphv, szdr_linear, K);

figure(4)

subplot(1,2,1)
hist(ZDR(:), -10:0.5:10)
xlabel('ZDR (dB)', 'FontSize', 14)

subplot(1,2,2)
hist(PHV(:), 0:0.01:1)
xlabel('Rhohv', 'FontSize',14)

% Compute DPSD
d = hann(M);
% d = rectwin(M);
d_rep = repmat(d, [1 size(V_h,2)]);
mm = repmat((0:M-1)', [1 size(V_h,2)]);
% Calculate DPSD the Original Way prior to bootstrap
DFT_H = fftshift(fft(d_rep .* V_h, M, 1));
DFT_V = fftshift(fft(d_rep .* V_v, M, 1));
% DFT_H = fftshift(sum(d_rep .* iqh .* exp(-j*2*pi*mm), 2));
% DFT_V = fftshift(sum(d_rep .* iqv .* exp(-j*2*pi*mm), 2));
SHF = nanmean(abs(DFT_H) .^ 2 / M, 2);
SVF = nanmean(abs(DFT_V) .^ 2 / M, 2);
SXF = nanmean(DFT_H .* conj(DFT_V) / M, 2);

figure(152)

subplot(2,1,1)
yyaxis left
plot(vv, 10*log10(SHF), '-k', vv, 10*log10(SVF), '-b')
xlabel('V_{r} (m s^{-1})', 'FontSize', 14)
ylabel('S_{H,V}', 'FontSize', 14)
yyaxis right
plot(-vv, 10*log10(SHF ./ SVF), '--k', vv, szdr, '--r')
ylabel('S(Z_{DR}) dB', 'FontSize', 14)

subplot(2,1,2)
yyaxis left
plot(vv, 10*log10(SHF), '-k', vv, 10*log10(SVF), '-b')
xlabel('V_{r} (m s^{-1})', 'FontSize', 14)
yyaxis right
plot(-vv, abs(SXF) ./ sqrt(SHF .* SVF), '--k', vv, sphv, '--r')
ylabel('S(\rho_{HV})', 'FontSize', 14)



