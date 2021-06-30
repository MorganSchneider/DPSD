
close all

% The data window for the
% analysis is a Blackman?Nuttall window, though it should be
% noted that other windows with sufficient sidelobe levels yield
% similar results.

% iqh.*nuttall(128)
% new_iq_windowed=iqh.*nuttall(128)

%% Load files

% Dimensions: (M, radius, azimuth, elev, dtype)

% set windowing function here!
win = nuttallwin(100);

% number of debris
dn = 100000;
els = {'0.5', '1.0', '1.5', '2.0', '2.5', '3.0', '3.5', '4.0', '4.5', '5.0'};


%%% Rain %%%
rain_dir = '~/Documents/code/DPSD/dpsd_outputs/simradar/rain';
rain = struct('vh', [], 'vv', [], 'szdr', [], 'sphv', [], 'psd', []);

for elx = 1:7
    ff = ['DPSD_sim-PPI' els{elx} '-DCU-nodebris.mat'];
    load([rain_dir '/' ff])
    
    rain.vh(:,:,:,elx) = iqh .* win;
    rain.vv(:,:,:,elx) = iqv .* win;
    rain.szdr(:,:,:,elx) = real(10*log10(sZDR.cf));
    rain.sphv(:,:,:,elx) = sPHV.cf;
end

M = size(iqh,1);
    
%%% Multi %%%
    
mult_dir = '~/Documents/code/DPSD/dpsd_outputs/simradar/multi';
mult = struct('vh', [], 'vv', [], 'szdr', [], 'sphv', [], 'psd', []);
    
for elx = 1:7
    for dt = [1 3] % add windowing function
        % nuttallwin(N)
        % zero padding for FFT
        ff = ['DPSD_sim-PPI' els{elx} '-DCU-d' num2str(dt) 'n' num2str(dn) '.mat'];
        load([mult_dir '/' ff])
        
        mult.vh(:,:,:,elx,(dt+1)/2) = iqh .* win;
        mult.vv(:,:,:,elx,(dt+1)/2) = iqv .* win;
        mult.szdr(:,:,:,elx,(dt+1)/2) = real(10*log10(sZDR.cf));
        mult.sphv(:,:,:,elx,(dt+1)/2) = sPHV.cf;
    end
    
end

%% Find debris velocities

deb = struct('vh', [], 'psd', []);


deb.vh = mult.vh - repmat(rain.vh, [1 1 1 1 2]);
zh = fftshift(fft(deb.vh, M, 1));
deb.psd = flip(abs(zh).^2 / M);

zhd = fftshift(fft(mult.vh, M, 1));
mult.psd = flip(abs(zhd).^2 / M);

zhr = fftshift(fft(rain.vh, M, 1));
rain.psd = flip(abs(zhr).^2 / M);


% CHECK PARAMS for 14,20 in d3 0.5deg
rdx = 14; % rad index
azdx = 20; % az index
eldx = 7; % elev index
ddx = 2; % dtype

figure(1)
subplot(1,2,1)
% rdx=14,azdx=20,eldx=1,ddx=3
yyaxis left
semilogy(vvx, abs(mult.psd(:,rdx,azdx,eldx,ddx)), '-k', 'LineWidth', 0.8)
hold on
semilogy(vvx, abs(deb.psd(:,rdx,azdx,eldx,ddx)), '-r', 'LineWidth', 0.8)
semilogy(vvx, abs(rain.psd(:,rdx,azdx,eldx)), '--b', 'LineWidth', 0.8)
hold off
yyaxis right
plot(vvx, mult.sphv(:,rdx,azdx,eldx,ddx), '-k')
legend('multi', 'debris', 'rain', 'Location', 'northwest')

% do this in dB instead of linear
subplot(1,2,2)
yyaxis left
plot(vvx, 10*log10(mult.psd(:,rdx,azdx,eldx,ddx)) - 10*log10(rain.psd(:,rdx,azdx,eldx)), '-r', 'LineWidth', 0.8)
hold on
plot(vvx, 2*ones(1,length(vvx)), '-k', 'LineWidth', 0.8)
plot(vvx, 5*ones(1,length(vvx)), '-k', 'LineWidth', 0.8)
hold off
yyaxis right
plot(vvx, mult.sphv(:,rdx,azdx,eldx,ddx), '-k');



psd_diff = 10*log10(mult.psd) - repmat(10*log10(rain.psd), [1 1 1 1 size(mult.vh,5)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% debris ID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thres = 5;
pthres = 1e4;


d = struct('el', [], 'szdr', [], 'sphv', []);
d.el = struct('szdr', [], 'sphv', []);
r = struct('el', [], 'szdr', [], 'sphv', []);
r.el = struct('szdr', [], 'sphv', []);

d_szdr = []; d_sphv = []; r_szdr = []; r_sphv = [];

for dd = 1:2
    
    dszdr_tot = []; dsphv_tot = []; rszdr_tot = []; rsphv_tot = [];
    for ee = 1:7
        
        dszdr = []; dsphv = []; rszdr = []; rsphv = [];
        for rr = 1:21
            for aa = 1:48
                d_inds = find(psd_diff(:,rr,aa,ee,dd) >= thres & ...
                    abs(mult.psd(:,rr,aa,ee,dd)) > pthres);
                r_inds = find(psd_diff(:,rr,aa,ee,dd) < thres & ...
                    abs(mult.psd(:,rr,aa,ee,dd)) > pthres);
                
                tmp = squeeze(mult.szdr(d_inds,rr,aa,ee,dd));
                dszdr = [dszdr; tmp];
                dszdr_tot = [dszdr_tot; tmp];
                
                tmp = squeeze(mult.sphv(d_inds,rr,aa,ee,dd));
                dsphv = [dsphv; tmp];
                dsphv_tot = [dsphv_tot; tmp];
                
                tmp = squeeze(mult.szdr(r_inds,rr,aa,ee,dd));
                rszdr = [rszdr; tmp];
                rszdr_tot = [rszdr_tot; tmp];
                
                tmp = squeeze(mult.sphv(r_inds,rr,aa,ee,dd));
                rsphv = [rsphv; tmp];
                rsphv_tot = [rsphv_tot; tmp];
            end
        end
        
        d(dd).el(ee).szdr = dszdr;
        d(dd).el(ee).sphv = dsphv;
        r(dd).el(ee).szdr = rszdr;
        r(dd).el(ee).sphv = rsphv;
        
    end
    
    d(dd).szdr = dszdr_tot;
    d(dd).sphv = dsphv_tot;
    r(dd).szdr = rszdr_tot;
    r(dd).sphv = rsphv_tot;
    
    d_szdr = [d_szdr; dszdr_tot];
    d_sphv = [d_sphv; dsphv_tot];
    r_szdr = [r_szdr; rszdr_tot];
    r_sphv = [r_sphv; rsphv_tot];
end



sphv_bins = 0.025:0.025:1;
szdr_bins = -20:1:20;


%% PLOTS


% Plot SPHV grouped by debris type

% for dd = 1:4
%     figure()
%     for ee = 1:3
%         subplot(2,3,ee)
%         histogram(d(dd).el(ee).sphv, sphv_bins)
%         title(['d' num2str(dd) ' el=' els{ee} '.0 (thres=' num2str(thres) 'dB)'])
%         xlabel('sphv')
%         
%         subplot(2,3,ee+3)
%         histogram(r(dd).el(ee).sphv, sphv_bins)
%         title(['rain el=' els{ee} '.0 (thres=' num2str(thres) 'dB)'])
%         xlabel('sphv')
%     end
%     
%     figure()
%     for ee = 1:3
%         subplot(2,3,ee)
%         histogram(d(dd).el(ee).szdr, szdr_bins)
%         title(['d' num2str(dd) ' el=' els{ee} '.0 (thres=' num2str(thres) 'dB)'])
%         xlabel('szdr')
%         
%         subplot(2,3,ee+3)
%         histogram(r(dd).el(ee).szdr, szdr_bins)
%         title(['rain el=' els{ee} '.0 (thres=' num2str(thres) 'dB)'])
%         xlabel('szdr')
%     end
% end
% 
% figure()
% for dd = 1:4
%     subplot(2,4,dd)
%     histogram(d(dd).sphv, sphv_bins)
%     title(['d' num2str(dd) ' (thres=' num2str(thres) 'dB)'])
%     xlabel('sphv')
%     
%     subplot(2,4,dd+4)
%     histogram(r(dd).sphv, sphv_bins)
%     title(['rain (thres=' num2str(thres) 'dB)'])
%     xlabel('sphv')
% end
% 
% figure()
% for dd = 1:4
%     subplot(2,4,dd)
%     histogram(d(dd).szdr, szdr_bins)
%     title(['d' num2str(dd) ' (thres=' num2str(thres) 'dB)'])
%     xlabel('szdr')
%     
%     subplot(2,4,dd+4)
%     histogram(r(dd).szdr, szdr_bins)
%     title(['rain (thres=' num2str(thres) 'dB)'])
%     xlabel('szdr')
% end


% Plot everything together

figure()

subplot(2,2,1)
histogram(d_szdr, szdr_bins)
title(['Debris SZDR values (thres=' num2str(thres) 'dB)'])
xlabel('szdr')

subplot(2,2,2)
histogram(d_sphv, sphv_bins)
title(['Debris SPHV values (thres=' num2str(thres) 'dB)'])
xlabel('sphv')

subplot(2,2,3)
histogram(r_szdr, szdr_bins)
title(['Rain SZDR values (thres=' num2str(thres) 'dB)'])
xlabel('szdr')

subplot(2,2,4)
histogram(r_sphv, sphv_bins)
title(['Rain SPHV values (thres=' num2str(thres) 'dB)'])
xlabel('sphv')


%% Check against actual debris DPSD

%%% Debris %%%

debr_dir = '~/Documents/code/DPSD/dpsd_outputs/simradar/debris/x100';
debris = struct('vh', [], 'vv', [], 'szdr', [], 'sphv', [], 'psd', []);

for elx = 1:7
    for dt = [1 3]
        ff = ['DPSD_sim-PPI' els{elx} '-TCU-d' num2str(dt) 'n' num2str(dn) '.mat'];
        load([debr_dir '/' ff])
        
        debris.vh(:,:,:,elx,(dt+1)/2) = iqh .* win;
        debris.vv(:,:,:,elx,(dt+1)/2) = iqv .* win;
        debris.szdr(:,:,:,elx,(dt+1)/2) = real(10*log10(sZDR.cf));
        debris.sphv(:,:,:,elx,(dt+1)/2) = sPHV.cf;
    end
end

% actual debris PSD
zh = fftshift(fft(debris.vh, M, 1));
debris.psd = flip(abs(zh).^2 / M);


%% Moving variance

n = 9; % number of pts for moving variance
%var_bins = stuff; % variance bins for plotting

sphv_var_rain = movvar(rain.sphv, n,0,1,'omitnan','Endpoints','shrink');
sphv_var_debr = movvar(debris.sphv, n,0,1,'omitnan','Endpoints','shrink');
sphv_var_mult = movvar(mult.sphv, n,0,1,'omitnan','Endpoints','shrink');



% sphv_var_rain = sphv_var_rain(rain.psd >= pthres);
% sphv_var_debr = sphv_var_debr(mult.psd >= pthres);
% sphv_var_mult = sphv_var_mult(mult.psd >= pthres);

pvar_bins = 0.005:0.005:0.3;

figure(14)

subplot(1,3,1)
histogram(sphv_var_rain, pvar_bins)
title(['Rain ' num2str(n) 'pt running variance'])
%ylim([0 1e4])

subplot(1,3,2)
histogram(sphv_var_debr, pvar_bins)
title(['Debris ' num2str(n) 'pt running variance'])
%ylim([0 10e4])

subplot(1,3,3)
histogram(sphv_var_mult, pvar_bins)
title(['Multi ' num2str(n) 'pt running variance'])
%ylim([0 10e4])




figure(15)

subplot(2,3,1)
%histogram(debris.szdr(abs(mult.psd) > pthres), szdr_bins)
yyaxis left
histogram(debris.szdr, szdr_bins)
ylim([0 2e5])
yyaxis right
plot(szdr, mem_zdr_debr, '-k', 'LineWidth', 2)
ylim([0 2])
%histfit(reshape(debris.szdr, [1612800 1]), 61)
title('Debris sZDR', 'FontSize', 20)

subplot(2,3,4)
yyaxis left
histogram(rain.szdr, szdr_bins)
ylim([0 2e5])
yyaxis right
plot(szdr, mem_zdr_rain, '-k', 'LineWidth', 2)
ylim([0 2])
%histfit(reshape(rain.szdr, [403200 1]), 61)
title('Rain sZDR', 'FontSize', 20)

subplot(2,3,2)
%histogram(debris.sphv(abs(mult.psd) > pthres), sphv_bins)
yyaxis left
histogram(debris.sphv, sphv_bins)
ylim([0 6e5])
yyaxis right
plot(sphv, mem_phv_debr, '-k', 'LineWidth', 2)
ylim([0 2])
%histfit(reshape(debris.sphv(debris.sphv>0), [1612800 1]), 21, 'gamma')
title('Debris sPHV', 'FontSize', 20)

subplot(2,3,5)
yyaxis left
histogram(rain.sphv, sphv_bins)
yyaxis right
plot(sphv, mem_phv_rain, '-k', 'LineWidth', 2)
ylim([0 2])
%histfit(reshape(rain.sphv, [403200 1]), 21, 'gamma')
title('Rain sPHV', 'FontSize', 20)

subplot(2,3,3)
yyaxis left
histogram(sphv_var_debr, pvar_bins)
yyaxis right
plot(svar, mem_var_debr, '-k', 'LineWidth', 2)
ylim([0 2])
ylabel('Membership Function', 'FontSize', 20)
title(['Debris ' num2str(n) 'pt moving var'], 'FontSize', 20)

subplot(2,3,6)
yyaxis left
histogram(sphv_var_rain, pvar_bins)
yyaxis right
plot(svar, mem_var_rain, '-k', 'LineWidth', 2)
ylim([0 2])
ylabel('Membership Function', 'FontSize', 20)
title(['Rain ' num2str(n) 'pt moving var'], 'FontSize', 20)

