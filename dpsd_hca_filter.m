%sim_dir = ['~/Documents/code/DPSD/test_sims/' signal_type];
%filename = blib('choosefile', sim_dir, '*.mat');

les = 'twocell';
simdate = '200630';
rain_concept = 'DCU';
debr_concept = 'TCU';
mult_concept = 'DCU';
dtypes = 0:14;
dnums = [10000 100000 1000000];

K = 20;

%%%%% IT'S NOT WORKINGGGGGGGGG %%%%%


for dt = dtypes
    if dt == 0 % no debris
        dn = [];
        concept = rain_concept;
        sim_dir = ['~/Documents/sims/' les '/' simdate '/nodebris/'];
        fname = ['sim-*-', concept, '-nodebris.mat'];
        mat_files = dir([sim_dir fname]);
        nels = length(mat_files);
        
        for i = 1:nels
            filename = mat_files(i).name;
            save_fname = ['dpsd-' filename(5:end)];
            load([sim_dir filename]);
            
            if size(iqh,4) > 1
                iqh = iqh(:,:,:,1);
                iqv = iqv(:,:,:,1);
            end
            
            iqh = permute(iqh, [3 1 2]);
            iqv = permute(iqv, [3 1 2]);
            M = size(iqh, 1);
            rlzs = 1;
            d = hann(M);
            lambda = dat.params.lambda;
            Ts = dat.params.prt;
            va = dat.params.va;
            vv = linspace(-va+1, va, M);
            vvx = flip(vv);
            
            params = struct('zdr', zdr, 'phv', rhohv, 'K', K, 'xx', xx, 'yy', yy, 'zz', zz);
            
            [sZDR, sPHV, sPSD] = dpsd_calc(iqh, iqv, d, M, K, rlzs);
            
            save([sim_dir save_fname],...
                'sZDR', 'sPHV', 'sPSD', 'iqh', 'iqv', 'vvx', 'd', 'params');
            
        end
        
        
        
    else % with debris
        concept = mult_concept;
        for dn = dnums
            sim_dir = ['~/Documents/sims/' les '/' simdate '/debris' num2str(dt) '/'];
            fname = ['sim-*-', concept, '-d', num2str(dt), 'n', num2str(dn), '.mat'];
            mat_files = dir([sim_dir fname]);
            nels = length(mat_files);
            for i = 1:nels
                filename = mat_files(i).name;
                save_fname = ['dpsd-' filename(5:end)];
                load([sim_dir filename]);
                
                if size(iqh,4) > 1
                    iqh = iqh(:,:,:,1);
                    iqv = iqv(:,:,:,1);
                end
                
                iqh = permute(iqh, [3 1 2]);
                iqv = permute(iqv, [3 1 2]);
                M = size(iqh, 1);
                rlzs = 1;
                d = hann(M);
                lambda = dat.params.lambda;
                Ts = dat.params.prt;
                va = dat.params.va;
                vv = linspace(-va+1, va, M);
                vvx = flip(vv);
                
                params = struct('zdr', zdr, 'phv', rhohv, 'K', K, 'xx', xx, 'yy', yy, 'zz', zz);
                
                [sZDR, sPHV, sPSD] = dpsd_calc(iqh, iqv, d, M, K, rlzs);
                
                save([sim_dir save_fname],...
                    'sZDR', 'sPHV', 'sPSD', 'iqh', 'iqv', 'vvx', 'd', 'params');
                
            end
            
            
        end
    end
end

%%% Have to add this in somewhere %%%


% plot ZDR spectrum
figure(3)
subplot(2,2,1)
yyaxis left
semilogy(vvx, abs(sSH.f(:,ind1,ind2)), 'k', 'LineWidth', 1)
xlabel('Doppler velocity {\it v_r}')
ylabel('PSD')
yyaxis right
plot(vvx, 10*log10(sZDR.f(:,ind1,ind2)), '-b')
% hold on
% plot(vvx, 10*log10(sZDR.auf(:,ind1,ind2)), '-r')
% hold off
% legend('PSD', 'morgan DPSD', 'arturo DPSD', 'Location', 'southwest')
ylim([-10 10])
ylabel('sZ_{DR}')
xlim([-va va])
axis square
grid on
title(['Mean H-channel PSD, sZ_{DR} ({\it K}=', num2str(K), ')'])
xlabel('Doppler velocity {\it v_r}')


subplot(2,2,2)
yyaxis left
semilogy(vvx, abs(sSV.f(:,ind1,ind2)), 'k', 'LineWidth', 1)
xlabel('Doppler velocity {\it v_r}')
ylabel('PSD')
yyaxis right
plot(vvx, 10*log10(sZDR.f(:,ind1,ind2)), '-b')
% hold on
% plot(vvx, 10*log10(sZDR.auf(:,ind1,ind2)), '-r')
% hold off
ylim([-10 10])
ylabel('sZ_{DR}')
xlim([-va va])
axis square
grid on
title(['Mean V-channel PSD, sZ_{DR} ({\it K}=', num2str(K), ')'])
xlabel('Doppler velocity {\it v_r}')

subplot(2,2,3)
plot(vvx, 10*log10(sZDR.f(:,ind1,ind2)), 'b')
% hold on
% plot(vvx, 10*log10(sZDR.auf(:,ind1,ind2)), 'r')
% hold off
xlim([-va, va])
ylim([-10 10])
axis square
grid on
title('{\it s}Z_{DR}({\it v_r}) in dBZ')
xlabel('Doppler velocity {\it v_r}')

subplot(2,2,4)
plot(vvx, 10*log10(sZDR.cf(:,ind1,ind2)), 'b')
% hold on
% plot(vvx, 10*log10(sZDR.auf(:,ind1,ind2)), 'r')
% hold off
xlim([-va, va])
ylim([-10 10])
axis square
grid on
title('Bias-corrected {\it s}Z_{DR}({\it v_r}) in dBZ')
xlabel('Doppler velocity {\it v_r}')

set(gcf, 'Units', 'Inches', 'Position', [2 2 12 10])
if plot_save_flag
    print(['~/Documents/imgs/DPSD/' iq_source '_' signal_type '_sZDR'], '-dpng')
end

% plot rhoHV spectrum
figure(4)
subplot(2,2,1)
yyaxis left
semilogy(vvx, abs(sSH.f(:,ind1,ind2)), 'k', 'LineWidth', 1)
xlabel('Doppler velocity {\it v_r}')
ylabel('PSD')
yyaxis right
plot(vvx, sPHV.f(:,ind1,ind2), '-b')
% hold on
% plot(vvx, sPHV.auf(:,ind1,ind2), '-r')
% hold off
% legend('PSD', 'morgan DPSD', 'arturo DPSD', 'Location', 'southwest')
ylim([0 1])
ylabel('s\rho_{HV}')
axis square
grid on
title(['Mean H-channel PSD, s\rho_{HV} ({\it K}=', num2str(K), ')'])

subplot(2,2,2)
yyaxis left
semilogy(vvx, abs(sSV.f(:,ind1,ind2)), 'k', 'LineWidth', 1)
ylabel('PSD')
yyaxis right
plot(vvx, sPHV.f(:,ind1,ind2), '-b')
% hold on
% plot(vvx, sPHV.auf(:,ind1,ind2), '-r')
% hold off
ylim([0 1])
ylabel('s\rho_{HV}')
xlim([-va va])
axis square
grid on
title(['Mean V-channel PSD, s\rho_{HV} ({\it K}=', num2str(K), ')'])
xlabel('Doppler velocity {\it v_r}')

subplot(2,2,3)
plot(vvx, sPHV.f(:,ind1,ind2), 'b')
% hold on
% plot(vvx, sPHV.auf(:,ind1,ind2), 'r')
% hold off
xlim([-va, va])
ylim([0 1])
axis square
grid on
title('{\it s}\rho_{HV}({\it v_r})')
xlabel('Doppler velocity {\it v_r}')

subplot(2,2,4)
plot(vvx, sPHV.cf(:,ind1,ind2), 'b')
% hold on
% plot(vvx, sPHV.auf(:,ind1,ind2), 'r')
% hold off
xlim([-va, va])
ylim([0 1])
axis square
grid on
title('Bias-corrected {\it s}\rho_{HV}({\it v_r})')
xlabel('Doppler velocity {\it v_r}')

set(gcf, 'Units', 'Inches', 'Position', [2 2 12 10])

if plot_save_flag
    print(['~/Documents/imgs/DPSD/' iq_source '_' signal_type '_sPHV'], '-dpng')
end



