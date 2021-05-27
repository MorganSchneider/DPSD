% HCA code :(

% start by removing debris entirely to filter, then make more complicated
% as needed (this is much later!!!)

% rain 0.9-0.95, 0.95-1 phv
% debris anything phv

% debris zdr -10-5 at 1, taper outside to something (10 dB each side?)
% rain zdr 0-3, taper over ~2dB?

% variance rain 0-0.02, play with tapering (0.05, 0.1?)
% variance debris start at 0.05 taper up to 0.1, taper down to 0.2(5),
% start at 0.5 for function value, end point something like 0.3

% keep dummy values at top to test

%%%%%% IT KIND OF WORKS????????? 03/29/21

%%%%% ADD SNR TO HCA AT SOME POINT --- want to KEEP low SNR components in
%%%%% rain signal, don't filter out with debris

%clf

rdx = 16; % good: 16 20 1 2; 15 21 1 2
adx = 20; % bad: 15 19 7 2
edx = 1; % overcorrected: 16 20 7 2
ddx = 2; % no change: 15 20 1 2; 15 20 7 2

psd = mult.psd(:,rdx,adx,edx,ddx);
psd_rain = rain.psd(:,rdx,adx,edx);

szdr = mult.szdr(:,rdx,adx,edx,ddx);
sphv = mult.sphv(:,rdx,adx,edx,ddx);
svar = sphv_var_mult(:,rdx,adx,edx,ddx);
% szdr = linspace(-25,25,100);
% sphv = linspace(0,1,100);
% svar = linspace(0,0.3,100);
vel = vvx;

M = length(vel);

obj_class = zeros(M,1);

%% Membership functions (turn this into a separate matlab function)

tpr = [0.9, 0.95]; % threshold phv rain
tpd = tpr(1); % threshold phv debris
tzr = [-2, 0, 3, 5]; % threshold zdr rain
tzd = [-20, -10, 5, 15]; % threshold zdr debris
tvr = [0.02, 0.05]; % threshold var rain
tvd = [0, 0.05, 0.15, 0.25]; % threshold var debris

 
% test for debris
mem_zdr_debr = zeros(M,1);
mem_phv_debr = ones(M,1);
mem_var_debr = ones(M,1);

mem_phv_debr(sphv < tpd) = 1;
mem_phv_debr(sphv >= tpd) = 0.5;

%inds = find(szdr < thres_zdr_debr(1));
mem_zdr_debr(szdr < tzd(1)) = 0;
%inds = find(szdr < thres_zdr_debr(2) & szdr >= thres_zdr_debr(1));
mem_zdr_debr(szdr < tzd(2) & szdr >= tzd(1)) = ...
    (szdr(szdr < tzd(2) & szdr >= tzd(1)) - tzd(1)) / (tzd(2) - tzd(1));
%inds = find(szdr < thres_zdr_debr(3) & szdr >= thres_zdr_debr(2));
mem_zdr_debr(szdr < tzd(3) & szdr >= tzd(2)) = 1;
%inds = find(szdr < thres_zdr_debr(4) & szdr >= thres_zdr_debr(3));
mem_zdr_debr(szdr < tzd(4) & szdr >= tzd(3)) = ...
    1 - ((szdr(szdr < tzd(4) & szdr >= tzd(3)) - tzd(3)) / (tzd(4) - tzd(3)));
%inds = find(szdr >= thres_zdr_debr(4));
mem_zdr_debr(szdr >= tzd(4)) = 0;


mem_var_debr(svar < tvd(1)) = 0;
mem_var_debr(svar < tvd(2) & svar >= tvd(1)) = ...
    (svar(svar < tvd(2) & svar >= tvd(1)) - tvd(1)) / (tvd(2) - tvd(1));
mem_var_debr(svar < tvd(3) & svar >= tvd(2)) = 1;
mem_var_debr(svar < tvd(4) & svar >= tvd(3)) = ...
    1 - ((svar(svar < tvd(4) & svar >= tvd(3)) - tvd(3)) / (tvd(4) - tvd(3)));
mem_var_debr(svar >= tvd(4)) = 0;


% test for rain
mem_zdr_rain = zeros(M,1);
mem_phv_rain = ones(M,1);
mem_var_rain = ones(M,1);


mem_phv_rain(sphv < tpr(1)) = 0;
mem_phv_rain(sphv < tpr(2) & sphv >= tpr(1)) = ...
    (sphv(sphv < tpr(2) & sphv >= tpr(1)) - tpr(1)) / (tpr(2) - tpr(1));
mem_phv_rain(sphv >= tpr(2)) = 1;

mem_zdr_rain(szdr < tzr(1)) = 0;
mem_zdr_rain(szdr < tzr(2) & szdr >= tzr(1)) = ...
    (szdr(szdr < tzr(2) & szdr >= tzr(1)) - tzr(1)) / (tzr(2) - tzr(1));
mem_zdr_rain(szdr < tzr(3) & szdr >= tzr(2)) = 1;
mem_zdr_rain(szdr < tzr(4) & szdr >= tzr(3)) = ...
    1 - ((szdr(szdr < tzr(4) & szdr >= tzr(3)) - tzr(3)) / (tzr(4) - tzr(3)));
mem_zdr_rain(szdr >= tzr(4)) = 0;

mem_var_rain(svar < tvr(1)) = 1;
mem_var_rain(svar < tvr(2) & svar >= tvr(1)) = ...
    1 - ((svar(svar < tvr(2) & svar >= tvr(1)) - tvr(1)) / (tvr(2) - tvr(1)));
mem_var_rain(svar >= tvr(2)) = 0;



%% Signal classification

% Aggregation parameters
agg_debr = mem_zdr_debr + mem_phv_debr + mem_var_debr;
agg_rain = mem_zdr_rain + mem_phv_rain + mem_var_rain;

obj_class(agg_debr > agg_rain) = 1;


figure(20)
subplot(3,1,1)
plot(vel, szdr)
title('sZDR')

subplot(3,1,2)
plot(vel, mem_zdr_rain)
title('Rain membership fx')

subplot(3,1,3)
plot(vel, mem_zdr_debr)
title('Debris membership fx')


figure(21)
subplot(3,1,1)
plot(vel, sphv)
title('sPHV')

subplot(3,1,2)
plot(vel, mem_phv_rain)
title('Rain membership fx')

subplot(3,1,3)
plot(vel, mem_phv_debr)
title('Debris membership fx')


figure(22)
subplot(3,1,1)
plot(vel, svar)
title('sPHV variance')

subplot(3,1,2)
plot(vel, mem_var_rain)
title('Rain membership fx')

subplot(3,1,3)
plot(vel, mem_var_debr)
title('Debris membership fx')


figure(23)
subplot(3,1,1)
plot(vel, agg_rain)
title('Rain aggregation param')

subplot(3,1,2)
plot(vel, agg_debr)
title('Debris aggregation param')

subplot(3,1,3)
plot(vel, obj_class)
title('HCA class (debr=1 rain=0)')


figure(24)
subplot(2,1,1)
yyaxis left
semilogy(vel, abs(psd), '-k', 'LineWidth', 1.5)
hold on
semilogy(vel, abs(rain.psd(:,rdx,adx,edx)), '-b', 'LineWidth', 0.8)
hold off
ylabel('PSD', 'FontSize', 20)
ylim([1e0 1e8])
xlim([-100 100])
yyaxis right
plot(vel, sphv, '-r', 'LineWidth', 1.5)
ylabel('sPHV', 'FontSize', 20)
title('Doppler spectrum + s\rho_{HV}', 'FontSize', 22)
%xlabel('Velocity')

subplot(2,1,2)
plot(vel, obj_class, '-k', 'LineWidth', 1.5)
title('HCA class (debris=1, rain=0)', 'FontSize', 22)
xlabel('Velocity', 'FontSize', 20)
ylim([-0.5 1.5])


figure(25)
subplot(2,1,1)
yyaxis left
semilogy(vel, abs(psd), '-k', 'LineWidth', 1.5)
hold on
semilogy(vel, abs(psd_rain), '-b', 'LineWidth', 0.8)
hold off
ylabel('PSD', 'FontSize', 20)
xlim([-100 100])
yyaxis right
plot(vel, szdr, '-r', 'LineWidth', 1.5)
ylabel('sZDR', 'FontSize', 20)
%xlabel('Velocity')
title('Doppler spectrum + sZ_{DR}', 'FontSize', 22)

subplot(2,1,2)
plot(vel, obj_class, '-k', 'LineWidth', 1.5)
title('HCA class (debris=1, rain=0)', 'FontSize', 22)
xlabel('Velocity', 'FontSize', 20)
ylim([-0.5 1.5])

%% Weighting by variance/rain aggregation param for debris mitigation
% plan B after recalculating using only rain coefficients (next
% section) it's just easier to code so here we are

var_norm = svar / max(svar); % normalize variance
b1 = 2; % aggressiveness of filter
var_inv = (1 - var_norm).^b1; % invert normalized variance
psd_new1 = var_inv .* psd;

figure()
subplot(2,1,1)
semilogy(vel, psd, '-k', 'LineWidth', 1.5)
title('Original Doppler Spectrum', 'FontSize', 22)
xlim([-100 100])
ylim([1e0 1e9])

subplot(2,1,2)
semilogy(vel, psd_new1, '-k', 'LineWidth', 1.5)
title('Filtered Doppler Spectrum', 'FontSize', 22)
xlim([-100 100])
ylim([1e0 1e9])

% this makes almost no difference unless you exponentiate(?) the inverted
% variance -- 2 or 3 for less aggressive


b2 = 2; % aggressiveness of filter
psd_new2 = (agg_rain.^2) .* psd;
% maybe plot things someday



%% Sebastian this is the part i need help with

psd_filt = psd; % psd = original psd, psd_filt = filtered psd, this is all linear not dB
obj_class(1) = 0; % can't remember why I did this tbh but I didn't want endpoints marked as debris?
obj_class(end) = 0;
psd_filt(obj_class==1) = 0; % set all debris coefficients to 0


figure(26) % random plot not really necessary
subplot(2,1,1)
semilogy(vel, psd_filt, '-k', 'LineWidth', 1.5)
title('Filtered Doppler Spectrum', 'FontSize', 22)
xlim([-100 100])
ylim([1e0 1e9])

subplot(2,1,2)
plot(vel, obj_class, '-k', 'LineWidth', 1.5)
title('HCA class (debris=1, rain=0)', 'FontSize', 22)
xlabel('Velocity', 'FontSize', 20)
ylim([-0.5 1.5])


M = length(vel); % vel is just my velocity vector for plotting spectra
km = find(psd_filt == max(psd_filt)) - M/2;

k = km - M/2 + 1: km + M/2;
Ts = dat.params.prt;

%P = 1/M * sum(psd); % not sure which I'm supposed to use for periodogram power?
%P = 1/M * sum(psd_filt);
P = sum(psd_filt);
lambda = dat.params.lambda;

sum_term = 0;
st2 = 0;
for i = 1:length(k) % is this loop right?? I'm pretty sure it's not
    sum_term = sum_term + ((k(i) - km) * psd_filt(mod(k(i)+M/2,M)+1)); % from D&Z
%    st2 = st2 + (k(i) * psd_filt(mod(k(i),M)+1)); % from Bob's wx radar notes
end

v1 = -lambda/(2*M) * (km/Ts + (1/(P*Ts) * sum_term)); % from D&Z

%v2 = -lambda/(2*M*Ts) * (st2 - M/2); % from Bob's notes

% the calculation from Doviak and Zrnic and the calculation from Bob's
% notes are giving me two very different numbers lol



%% linear interpolation and recalculation

psd_new4 = psd_filt;

nanx = isnan(psd_new4);
t = 1:numel(psd_new4);
psd_new4(nanx) = interp1(t(~nanx), psd_new4(~nanx), t(nanx));



% figure(27)
% subplot(2,1,1)
% semilogy(vel, psd_new, '-k', 'LineWidth', 1.5)
% title('Reconstructed Doppler Spectrum (linear interp)', 'FontSize', 22)
% xlim([-100 100])
% ylim([1e0 1e9])
% 
% subplot(2,1,2)
% semilogy(vel, psd, '-k', 'LineWidth', 1.5)
% title('Original Doppler Spectrum', 'FontSize', 22)
% xlim([-100 100])
% ylim([1e0 1e9])



% I realized this is totally wrong but I actually have no idea how to convert a
% periodogram back into a time series please help me
ts_old = ifft(ifftshift(psd));
ts_new = ifft(ifftshift(psd_new4));
ts_rain = ifft(ifftshift(psd_rain));

%vr = -dat.params.va / pi * angle(mean(iqh(:, :, 2:end,:) .* conj(iqh(:, :, 1:end-1,:)), 3));
vr_old = -100 / pi * angle(mean(ts_old(2:end) .* conj(ts_old(1:end-1))));
vr_new = -100 / pi * angle(mean(ts_new(2:end) .* conj(ts_new(1:end-1))));
vr_rain = -100 / pi * angle(mean(ts_rain(2:end) .* conj(ts_rain(1:end-1))));

figure(27)
subplot(3,1,1)
semilogy(vel, psd_new4, '-k', 'LineWidth', 1.5)
hold on
semilogy([vr_new, vr_new], [1e0, 1e9], '-.k', 'LineWidth', 1.5)
hold off
title('Reconstructed Doppler Spectrum (linear interp)', 'FontSize', 22)
xlim([-100 100])
ylim([1e0 1e9])

subplot(3,1,2)
semilogy(vel, psd, '-k', 'LineWidth', 1.5)
hold on
semilogy([vr_old, vr_old], [1e0, 1e9], '-.k', 'LineWidth', 1.5)
hold off
title('Original Doppler Spectrum', 'FontSize', 22)
xlim([-100 100])
ylim([1e0 1e9])

subplot(3,1,3)
semilogy(vel, psd_rain, '-k', 'LineWidth', 1.5)
hold on
semilogy([vr_rain, vr_rain], [1e0, 1e9], '-.k', 'LineWidth', 1.5)
hold off
title('Truth Doppler Spectrum', 'FontSize', 22)
xlim([-100 100])
ylim([1e0 1e9])
xlabel('Velocity', 'FontSize', 20)


