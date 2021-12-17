%% Load DPSD variables

LES = 'twocell';
dd_date = '210614';
nd_date = '211005';
dtypes = [1 3];
dnums = [10000 100000];


base_dir = ['~/Documents/sims/' LES];

elevs = 0.5:0.5:5;
rain_dir = [base_dir '/' nd_date '/nodebris'];
for el = elevs
    dpsd_file = [rain_dir '/dpsd-PPI' num2str(el,'%.1f') '-U-nodebris.mat'];
    if ~isfile(dpsd_file)
        sim_file = strrep(dpsd_file, 'dpsd', 'sim');
        if ~isfile(sim_file)
            
        end
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


elevs = 0.5:0.5:2;
for dt = dtypes
    sim_dir = [base_dir '/' dd_date '/debris' num2str(dt)];
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
rain.szdr_bins = -20:0.75:20;
pvar_bins = 0.005:0.005:0.25;
debris.zvar_bins = 0:1:80;
rain.zvar_bins = 0:0.05:3;

szdr_tst = linspace(-20,20,200);
sphv_tst = linspace(0,1,200);
pvar_tst = linspace(0,0.25,200);
zvar_tst = linspace(0,80,200);
svar_tst = struct('phv', pvar_tst, 'zdr', zvar_tst);

[obj_class, agg, memf] = hca_class(szdr_tst, sphv_tst, svar_tst);


histc = [0.5 0.5 0.6];
yt = 0:0.5:1;
yl = [0 1.25];

figure(1);
clf

ax = subplot(2,3,1);
yyaxis left
plot(szdr_tst, memf.zdr_debr, '-k', 'LineWidth', 2)
ylim(yl)
yticks(yt)
xlabel('sZ_{DR}', 'FontSize', 14)
ylabel('Debris membership', 'FontSize', 14)
yyaxis right
histogram(debris.szdr, debris.szdr_bins, 'FaceColor', histc)
xlim([-20 20])
ylim([0 2e5])
yticks(5e4:5e4:3e5)
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = histc;
title('(a) sZ_{DR}', 'FontSize', 14)
ax.SortMethod = 'depth';

ax(2) = subplot(2,3,2);
yyaxis left
plot(sphv_tst, memf.phv_debr, '-k', 'LineWidth', 2)
xlim([0 1])
ylim(yl)
yticks(yt)
xlabel('s\rho_{HV}', 'FontSize', 14)
yyaxis right
histogram(debris.sphv, debris.sphv_bins, 'FaceColor', histc)
ylim([0 4e5])
yticks(1e5:1e5:5e5)
ax(2).YAxis(1).Color = 'k';
ax(2).YAxis(2).Color = histc;
title('(b) s\rho_{HV}', 'FontSize', 14)
ax(2).SortMethod = 'depth';

ax(3) = subplot(2,3,3);
yyaxis left
plot(pvar_tst, memf.pvar_debr, '-k', 'LineWidth', 2)
xlim([0 0.25])
ylim(yl)
yticks(yt)
xlabel('\sigma^2', 'FontSize', 14)
yyaxis right
histogram(debris.pvar, pvar_bins, 'FaceColor', histc)
ylim([0 10e4])
yticks(2e4:2e4:12e4)
ax(3).YAxis(1).Color = 'k';
ax(3).YAxis(2).Color = histc;
title('(c) s\rho_{HV} variance', 'FontSize', 14)
ax(3).SortMethod = 'depth';

ax(4) = subplot(2,3,4);
yyaxis left
plot(szdr_tst, memf.zdr_rain, '-k', 'LineWidth', 2)
ylim(yl)
yticks(yt)
xlabel('sZ_{DR}', 'FontSize', 14)
ylabel('Rain membership', 'FontSize', 14)
yyaxis right
histogram(rain.szdr, rain.szdr_bins, 'FaceColor', histc)
xlim([-20 20])
ylim([0 3e5])
yticks(5e4:5e4:3e5)
ax(4).YAxis(1).Color = 'k';
ax(4).YAxis(2).Color = histc;
title('(d) sZ_{DR}', 'FontSize', 14)
ax(4).SortMethod = 'depth';

ax(5) = subplot(2,3,5);
yyaxis left
plot(sphv_tst, memf.phv_rain, '-k', 'LineWidth', 2)
xlim([0 1])
ylim(yl)
yticks(yt)
xlabel('s\rho_{HV}', 'FontSize', 14)
yyaxis right
histogram(rain.sphv, rain.sphv_bins, 'FaceColor', histc)
ylim([0 4e5])
yticks(1e5:1e5:5e5)
ax(5).YAxis(1).Color = 'k';
ax(5).YAxis(2).Color = histc;
title('(e) s\rho_{HV}', 'FontSize', 14)
ax(5).SortMethod = 'depth';

ax(6) = subplot(2,3,6);
yyaxis left
plot(pvar_tst, memf.pvar_rain, '-k', 'LineWidth', 2)
xlim([0 0.25])
ylim(yl)
yticks(yt)
xlabel('\sigma^2', 'FontSize', 14)
yyaxis right
histogram(rain.pvar, pvar_bins, 'FaceColor', histc)
ylim([0 10e3])
yticks(2e3:2e3:12e3)
ax(6).YAxis(2).Exponent = 3;
ax(6).YAxis(1).Color = 'k';
ax(6).YAxis(2).Color = histc;
title('(f) s\rho_{HV} variance', 'FontSize', 14)
ax(6).SortMethod = 'depth';


set(gcf, 'Units', 'Inches', 'Position', [1 10 14 7])

print('~/Documents/articles/2021/membership', '-dpng')
% annotation('textbox', [0.386 0.93 0.25 0.08], 'String', 'Debris Membership Functions',...
%     'FontSize',16,'FontWeight','bold','EdgeColor','none','HorizontalAlignment','center',...
%     'VerticalAlignment','middle')
% annotation('textbox', [0.386 0.46 0.25 0.08], 'String', 'Rain Membership Functions',...
%     'FontSize',16,'FontWeight','bold','EdgeColor','none','HorizontalAlignment','center',...
%     'VerticalAlignment','middle')



