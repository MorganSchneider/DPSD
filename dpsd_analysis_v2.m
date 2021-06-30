%% Load DPSD variables

LES = 'twocell';
sim_date = '210614';
dtypes = [1 3];
dnums = [10000 100000 1000000];


base_dir = ['~/Documents/sims/' LES '/' sim_date];

elevs = 0.5:0.5:5;
rain_dir = [base_dir '/nodebris'];
for el = elevs
    dpsd_file = [rain_dir '/dpsd-PPI' num2str(el,'%.1f') '-DCU-nodebris.mat'];
    if ~isfile(dpsd_file)
        sim_file = strrep(dpsd_file, 'dpsd', 'sim');
        load(sim_file)
        
        iqh = permute(iqh, [3 1 2]);
        iqv = permute(iqv, [3 1 2]);
        M = size(iqh, 1);
        d = nuttallwin(M); % data window
        
        [szdr, sphv, svar, ~] = dpsd_calc(iqh, iqv, d, M, 20, 1);
    else
        load(dpsd_file)
    end
    
    if ~exist('tmp1','var')
        tmp1 = szdr;
        tmp2 = sphv;
        tmp3 = svar.phv;
        tmp4 = svar.zdr;
    else
        tmp1 = cat(1, tmp1, szdr);
        tmp2 = cat(1, tmp2, sphv);
        tmp3 = cat(1, tmp3, svar.phv);
        tmp4 = cat(1, tmp4, svar.zdr);
    end
end
rain = struct('szdr', tmp1, 'sphv', tmp2, 'pvar', tmp3, 'zvar', tmp4);
clear tmp1 tmp2 tmp3 tmp4


elevs = [1 2];
for dt = dtypes
    sim_dir = [base_dir '/debris' num2str(dt)];
    for dn = dnums
        for el = elevs
            dpsd_file = [sim_dir '/dpsd-PPI' num2str(el,'%.1f') '-TC-d' num2str(dt) 'n' num2str(dn) '.mat'];
            
            if ~isfile(dpsd_file)
                sim_file = strrep(dpsd_file, 'dpsd', 'sim');
                load(sim_file)
                
                iqh = permute(iqh, [3 1 2]);
                iqv = permute(iqv, [3 1 2]);
                M = size(iqh, 1);
                d = nuttallwin(M); % data window
                
                [szdr, sphv, svar, ~] = dpsd_calc(iqh, iqv, d, M, 20, 1);
            else
                load(dpsd_file)
            end
            
            if ~exist('tmp1','var')
                tmp1 = szdr;
                tmp2 = sphv;
                tmp3 = svar.phv;
                tmp4 = svar.zdr;
            else
                tmp1 = cat(1, tmp1, szdr);
                tmp2 = cat(1, tmp2, sphv);
                tmp3 = cat(1, tmp3, svar.phv);
                tmp4 = cat(1, tmp4, svar.zdr);
            end
        end
    end
end
debris = struct('szdr', tmp1, 'sphv', tmp2, 'pvar', tmp3, 'zvar', tmp4);
clear tmp1 tmp2 tmp3 tmp4


%% Make histograms

debris.sphv_bins = 0.025:0.025:1;
rain.sphv_bins = 0.025:0.025:1;
debris.szdr_bins = -20:1:20;
rain.szdr_bins = -20:0.5:20;
pvar_bins = 0.005:0.005:0.25;
debris.zvar_bins = 0:1:80;
rain.zvar_bins = 0:0.05:3;

szdr_tst = linspace(-20,20,200);
sphv_tst = linspace(0,1,200);
pvar_tst = linspace(0,0.25,200);
zvar_tst = linspace(0,80,200);
svar_tst = struct('phv', pvar_tst, 'zdr', zvar_tst);

[obj_class, agg, memf] = hca_class(szdr_tst, sphv_tst, svar_tst);



figure(1);
clf

ax = subplot(3,2,1);
yyaxis right
histogram(debris.szdr, debris.szdr_bins)
xlim([-20 20])
ylim([0 2e5])
yticks(5e4:5e4:2e5)
%set(gca, 'YScale', 'log')
yyaxis left
plot(szdr_tst, memf.zdr_debr, '-k', 'LineWidth', 2)
ylim([0 1.5])
yticks([0, 0.5, 1])
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = [0.6 0.6 0.6];
grid on
title('(a) sZ_{DR,debris}', 'FontSize', 14)

ax(2) = subplot(3,2,2);
yyaxis right
histogram(rain.szdr, rain.szdr_bins)
xlim([-20 20])
ylim([0 2e5])
yticks(5e4:5e4:2e5)
yyaxis left
plot(szdr_tst, memf.zdr_rain, '-k', 'LineWidth', 2)
ylim([0 1.5])
yticks([0, 0.5, 1])
ax(2).YAxis(1).Color = 'k';
ax(2).YAxis(2).Color = [0.6 0.6 0.6];
grid on
title('(b) sZ_{DR,rain}', 'FontSize', 14)

ax(3) = subplot(3,2,3);
yyaxis right
histogram(debris.sphv, debris.sphv_bins)
ylim([0 4e5])
yticks(1e5:1e5:4e5)
%set(gca, 'YScale', 'log')
yyaxis left
plot(sphv_tst, memf.phv_debr, '-k', 'LineWidth', 2)
xlim([0 1])
ylim([0 1.5])
yticks([0, 0.5, 1])
ax(3).YAxis(1).Color = 'k';
ax(3).YAxis(2).Color = [0.6 0.6 0.6];
grid on
title('(c) s\rho_{HV,debris}', 'FontSize', 14)

ax(4) = subplot(3,2,4);
yyaxis right
histogram(rain.sphv, rain.sphv_bins)
ylim([0 4e5])
yticks(1e5:1e5:4e5)
%set(gca, 'YScale', 'log')
yyaxis left
plot(sphv_tst, memf.phv_rain, '-k', 'LineWidth', 2)
xlim([0 1])
ylim([0 1.5])
yticks([0, 0.5, 1])
ax(4).YAxis(1).Color = 'k';
ax(4).YAxis(2).Color = [0.6 0.6 0.6];
grid on
title('(d) s\rho_{HV,rain}', 'FontSize', 14)

ax(5) = subplot(3,2,5);
yyaxis right
histogram(debris.pvar, pvar_bins)
yticks(2e4:2e4:8e4)
%set(gca, 'YScale', 'log')
yyaxis left
plot(pvar_tst, memf.pvar_debr, '-k', 'LineWidth', 2)
xlim([0 0.25])
ylim([0 1.5])
yticks([0, 0.5, 1])
ax(5).YAxis(1).Color = 'k';
ax(5).YAxis(2).Color = [0.6 0.6 0.6];
grid on
title('(e) \sigma^2_{s\rho_{HV},debris}', 'FontSize', 14)


ax(6) = subplot(3,2,6);
yyaxis right
histogram(rain.pvar, pvar_bins)
yticks(2e3:2e3:1.2e4)
%set(gca, 'YScale', 'log')
yyaxis left
plot(pvar_tst, memf.pvar_rain, '-k', 'LineWidth', 2)
xlim([0 0.25])
ylim([0 1.5])
yticks([0, 0.5, 1])
ax(6).YAxis(1).Color = 'k';
ax(6).YAxis(2).Color = [0.6 0.6 0.6];
grid on
title('(f) \sigma^2_{s\rho_{HV},rain}', 'FontSize', 14)


% ax(7) = subplot(4,2,7);
% yyaxis left
% histogram(debris.zvar, debris.zvar_bins)
% yyaxis right
% plot(zvar_tst, memf.zvar_debr, '-k', 'LineWidth', 2)
% ylim([-0.1 1.5])
% yticks([0 1])
% ax(7).YAxis(1).Color = 'k';
% ax(7).YAxis(2).Color = 'k';
% title('(g) \sigma^2_{sZ_{DR}}', 'FontSize', 14)
% xlabel('Debris', 'FontSize', 20)
% 
% 
% ax(8) = subplot(4,2,8);
% yyaxis left
% histogram(rain.zvar, rain.zvar_bins)
% yyaxis right
% plot(zvar_tst, memf.zvar_rain, '-k', 'LineWidth', 2)
% xlim([0 5])
% ylim([-0.1 1.5])
% yticks([0 1])
% ax(8).YAxis(1).Color = 'k';
% ax(8).YAxis(2).Color = 'k';
% title('(h) \sigma^2_{sZ_{DR}}', 'FontSize', 14)
% xlabel('Rain', 'FontSize', 20)


