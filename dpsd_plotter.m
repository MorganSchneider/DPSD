%clear;
%close all

[iq_source, signal_type] = choose_source();

zmap = feval('boonlib', 'zmap', 21);
rgmap = feval('boonlib', 'rgmap', 21);

sim_dir = ['~/Documents/code/DPSD/dpsd_outputs/' iq_source '/' signal_type];
filename = blib('choosefile', sim_dir, '*.mat');
load(filename);

if strcmp(iq_source, 'simradar')
    inds = find(filename=='/');
    filename2 = ['~/Documents/code/DPSD/test_sims/' signal_type '/'...
        filename(inds(end)+6:end)];
    load(filename2, 'xx', 'yy', 'zz', 'zh', 'vr');
end

% DPSD file:
% 
% sZDR, sPHV    : Structure array of spectral ZDR and rhoHV estimates
%  sZDR.ms      - Uncorrected, unaveraged
%  sZDR.c       - Corrected, unaveraged
%  sZDR.f       - Uncorrected, averaged (if multiple realizations)
%  sZDR.cf      - Corrected, averaged (if multiple realizations)
% sSH, sSV, sSX : Horizontal, vertical, and cross-pol power spectral density
%  sSH.ms       - Unaveraged
%  sSH.f        - Averaged (if multiple realizations)
% iqh, iqv      : Horizontal and vertical I/Q time series
% vvx           : Velocity vector for plotting spectra
% d             : Windowing function
% params        : Structure array of other useful quantities
%  signal_type  - Prescribed signal type (i.e., rain, debris, or multi)
%  zdr          - Bulk ZDR value of gate
%  phv          - Bulk rhoHV value of gate
%  K            - Number of bootstrapped pseudorealizations used

% szdr = sZDR.cf;
% sphv = sPHV.cf;
ZDR = params.zdr;
PHV = params.phv;


%% 

r_ind = 14;
az_ind = 20;


%for r_ind = 1:size(iqh,2)
%    for az_ind = 1:size(iqh,3)
        
        figure(1)
        clf
        
        ha = subplot(2,2,1);
        hs = pcolor(xx, yy, zh(:,:,1));
        hold on
        plot(xx(r_ind, az_ind), yy(r_ind, az_ind), 'sk', 'MarkerSize', 11, 'MarkerFaceColor', 'w', 'LineWidth', 2);
        hold off
        set(gca, 'DataAspect', [1 1 1])
        caxis([0 80])
        colormap(ha, blib('zmap'))
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
        colormap(ha(2), blib('rgmap2'))
        colorbar
        shading flat
        title('Velocity (m/s)')
        
        
        subplot(2,2,3)
        yyaxis left
        semilogy(vvx, abs(sSH.f(:,r_ind,az_ind)), '-k', 'LineWidth', 1)
        xlabel('Doppler velocity {\it v_r}')
        ylabel('PSD')
        yyaxis right
        plot(vvx, 10*log10(sZDR.cf(:,r_ind,az_ind)), '-r', 'LineWidth', 0.6)
        ylabel('sZ_{DR}')
        title('Spectral Z_{DR}')
        %xlim([-va va])
        ylim([-10 10])
        %axis square
        grid on
        
        subplot(2,2,4)
        yyaxis left
        semilogy(vvx, abs(sSH.f(:,r_ind,az_ind)), '-k', 'LineWidth', 1)
        xlabel('Doppler velocity {\it v_r}')
        ylabel('PSD')
        yyaxis right
        plot(vvx, sPHV.cf(:,r_ind,az_ind), '-r', 'LineWidth', 0.6)
        ylabel('s\rho_{HV}')
        title('Spectral \rho_{HV}')
        %xlim([-va va])
        ylim([0 1])
        %axis square
        grid on
        
        
        axes('Unit', 'Normalized', 'Position', [0.5 0.5 0.01 0.01])
        title(['Debris type 3, ' signal_type ' signal (Z_{DR} = ',...
            num2str(ZDR(r_ind,az_ind),2) ' dBZ, \rho_{HV} = ' num2str(PHV(r_ind,az_ind),2) ')'],...
            'FontSize', 12)
        axis off
        
        print(['~/Documents/imgs/dpsd_' iq_source '_' signal_type '_d11'], '-dpng')
        
%        pause
%         save_fig = input(['Save figure? (y/n)' newline]);
%         if strcmp(save_fig, 'Y') || strcmp(save_fig, 'y')
%             
%         end

%    end
%end

if strcmp(signal_type, 'rain')
    szdr.rain = sZDR.cf;
    sphv.rain = sPHV.cf;
    zdr.rain = ZDR;
    phv.rain = PHV;
    ssh.rain = sSH.f;
    ssv.rain = sSV.f;
    vh.rain = iqh;
    vv.rain = iqv;
elseif strcmp(signal_type, 'debris')
    szdr.debris = sZDR.cf;
    sphv.debris = sPHV.cf;
    zdr.debris = ZDR;
    phv.debris = PHV;
    ssh.debris = sSH.f;
    ssv.debris = sSV.f;
    vh.debris = iqh;
    vv.debris = iqv;
elseif strcmp(signal_type, 'multi')
    szdr.multi = sZDR.cf;
    sphv.multi = sPHV.cf;
    zdr.multi = ZDR;
    phv.multi = PHV;
    ssh.multi = sSH.f;
    ssv.multi = sSV.f;
    vh.multi = iqh;
    vv.multi = iqv;
end

figure(2)
clf

semilogy(vvx, abs(ssh.rain(:,r_ind,az_ind)), '-k', 'LineWidth', 1)
hold on
semilogy(vvx, abs(ssh.debris(:,r_ind,az_ind)), '--b', 'LineWidth', 1)
semilogy(vvx, abs(ssh.multi(:,r_ind,az_ind)), '--k', 'LineWidth', 1)
semilogy(vvx, abs(ssh.multi(:,r_ind,az_ind) - ssh.debris(:,r_ind,az_ind)), 'r', 'LineWidth', 1)
hold off
grid on
title('sSH by signal type')
legend('Rain only', 'Debris only', 'Multi', 'Multi - Debris', 'Location', 'northeast')
xlabel('Velocity (m/s)')

if strcmp(signal_type, 'multi')
    print(['~/Documents/imgs/dpsd_' iq_source '_PSDs_d11'], '-dpng')
end





% figure(3)
% clf
% 
% ha = subplot(2,4,1);
% hs = pcolor(xx, yy, zh(:,:,1));
% hold on
% plot(xx(r_ind, az_ind), yy(r_ind, az_ind), 'sk', 'MarkerSize', 11, 'MarkerFaceColor', 'w', 'LineWidth', 2);
% hold off
% set(gca, 'DataAspect', [1 1 1])
% caxis([0 80])
% colormap(ha, blib('zmap'))
% colorbar
% shading flat
% set(gca, 'YDir', 'Normal')
% title('Reflectivity (dBZ)')
% 
% ha(2) = subplot(2,4,5);
% hs(2) = pcolor(xx, yy, vr(:,:,1));
% hold on
% plot(xx(r_ind, az_ind), yy(r_ind, az_ind), 'sk', 'MarkerSize', 11, 'MarkerFaceColor', 'w', 'LineWidth', 2);
% hold off
% set(gca, 'DataAspect', [1 1 1])
% caxis([-1 1] * round(max(max(abs(vr(:,:,1)))), -1))
% colormap(ha(2), blib('rgmap2'))
% colorbar
% shading flat
% title('Velocity (m/s)')
% 
% 
% subplot(2,4,2)
% yyaxis left
% semilogy(vvx, abs(sSH.f(:,r_ind,az_ind)), '-k', 'LineWidth', 1)
% xlabel('Doppler velocity {\it v_r}')
% ylabel('PSD')
% yyaxis right
% plot(vvx, 10*log10(sZDR.cf(:,r_ind,az_ind)), '-r', 'LineWidth', 0.6)
% ylabel('sZ_{DR}')
% title('Spectral Z_{DR}')
% %xlim([-va va])
% ylim([-10 10])
% %axis square
% grid on
% 
% subplot(2,2,6)
% yyaxis left
% semilogy(vvx, abs(sSH.f(:,r_ind,az_ind)), '-k', 'LineWidth', 1)
% xlabel('Doppler velocity {\it v_r}')
% ylabel('PSD')
% yyaxis right
% plot(vvx, sPHV.cf(:,r_ind,az_ind), '-r', 'LineWidth', 0.6)
% ylabel('s\rho_{HV}')
% title('Spectral \rho_{HV}')
% %xlim([-va va])
% ylim([0 1])
% %axis square
% grid on
% 
% subplot(2,4,3)
% 
% 
% subplot(2,4,7)
% 
% 
% subplot(2,4,4)
% 
% 
% subplot(2,4,8)
% 
% 
% axes('Unit', 'Normalized', 'Position', [0.5 0.5 0.01 0.01])
% title(['Debris type 14, ' signal_type ' signal (Z_{DR} = ',...
%     num2str(ZDR(r_ind,az_ind),2) ' dBZ, \rho_{HV} = ' num2str(PHV(r_ind,az_ind),2) ')'],...
%     'FontSize', 12)
% axis off


