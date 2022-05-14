
%close all

%%%%% Things to try???
% -stricter thresholds for debris to get more rain flags
%   -add Zh/bulk rhoHV threshold/dependency for filtering
%   -add SNR threshold for filtering
%       -take max value of PSD and correct points within 10/20/?? dB
% -look at more individual spectra --- ginput to click point on figure and
% retrieve x,y coords [x,y] = ginput(1)
%
% -test real data from David


% Good correction examples: twocell d3,n10000 el=2.0
%   r=15,az=21
%   r=12,az=21
%   r=12,az=19

% Plot idea:
% combine fig 8 plots to more easily compare size, concentration, elevation
% d3 1000000 vs d3 10000 at 0.5
% d3 1000000 vs d3 1000000 at 0.5
% d3 1000000 at 0.5 vs d3 1000000 at 5.0
% do this for PPIs and histograms

[st,~] = dbstack('-completenames');
if length(st) > 1
    external_call = 1;
else
    external_call = 0;
end



% Load IQ data
if ~external_call
    plot_flag = 1;
    var_save_flag = 1;
    plot_save_flag = 0;
    
    x = input(['Choose I/Q data source: ', newline,...
        '  (1) SimRadar', newline,...
        '  (2) KOUN', newline]);
end

if x == 1 %if SimRadar
    if ~external_call
        LES = 'twocell';
        base_dir = ['Documents/sims/' LES];
        data_dir = uigetdir(base_dir); % Location of SimRadar output files
        filename = blib('choosefile', data_dir, 'sim-*.mat');
    end
    load(filename);
    
    % Separate path name and file name for saving files later
    fname = erase(filename, [data_dir '/']);
    fname = erase(fname, '.mat');
    dind = regexp(fname,'d');
    nind = regexp(fname,'n');
    dtype = fname(dind+1: nind-1);
    dnum = fname(nind+1: end);
    eind = regexp(fname,'PPI');
    elev = fname(eind+3:eind+5);
    
    va = dat.params.va;
    Ts = dat.params.prt;
    lambda = dat.params.lambda;
    
    img_dir = ['Documents/imgs/' LES '/d' dtype '/DPSD/'];
    
    % Put all the bulk variables from checkiq.m into one structure
    params = struct('zh', zh, 'zdr', zdr, 'phv', rhohv, 'xx', xx, 'yy', yy, 'zz', zz, 'va', va);
    
elseif x == 2 %if KOUN
    if ~external_call
        data_dir = 'Documents/code/obsdata/';
        filename = blib('choosefile', data_dir, '*.mat');
        
        if ~isfile(filename)
            process_data
        end
    end
    
    load(filename);
    
    fname = erase(filename, data_dir);
    fname = erase(fname, '.mat');
    
    iqh = iqh_tds;
    iqv = iqv_tds;
    xx = xx(tds.xinds, tds.yinds);
    yy = yy(tds.xinds, tds.yinds);
    zz = zz(tds.xinds, tds.yinds);
    zh = zh(tds.xinds, tds.yinds);
    zdr = zdr(tds.xinds, tds.yinds);
    rhohv = rhohv(tds.xinds, tds.yinds);
    vr = vr(tds.xinds, tds.yinds);
    % vr_unfolded = vr_unfolded(tds.xinds, tds.yinds);
    
    img_dir = 'Documents/imgs/obsdata/';
    
    % Put all the bulk variables into one structure
    params = struct('zh', zh, 'zdr', zdr, 'phv', rhohv, 'xx', xx, 'yy', yy, 'zz', zz, 'va', va);
    
end

if ~exist(img_dir,'dir')
    mkdir(img_dir)
    addpath(genpath(img_dir))
    savepath
end


% Put all the bulk variables from checkiq.m into one structure
%params = struct('zh', zh, 'vr', vr, 'zdr', zdr, 'phv', rhohv, 'xx', xx, 'yy', yy, 'zz', zz, 'va', va);

% if simulation contains multiple sweeps
if size(iqh,4) > 1
    iqh = iqh(:,:,:,1);
    iqv = iqv(:,:,:,1);
end

% Rearrange dimensions to run through DPSD code
iqh = permute(iqh, [3 1 2]);
iqv = permute(iqv, [3 1 2]);
M = size(iqh, 1);
%va = dat.params.va;
vvx = flip(linspace(-va+1, va, M));
d = nuttallwin(M); % data window

% Uncorrected PSD and velocity
PSD.old = fftshift(abs(fft(iqh.*d, M, 1)).^2, 1);
vr_old = vr;

% Calculate DPSDs
[szdr, sphv, svar, spsd] = dpsd_calc(iqh, iqv, d, M, 20, 1);

%%
% HCA debris classification

[obj_class, agg, memf] = hca_class(szdr, sphv, svar);
%obj_class = movmean(obj_class, 3, 1, 'omitnan', 'Endpoints', 'shrink');


% Filter debris from PSD
PSD.dca = PSD.old;
obj_class(1) = 0; % can't remember why I did this tbh but I didn't want endpoints marked as debris?
obj_class(end) = 0;
PSD.dca(obj_class>0.33) = 0; % set all debris coefficients to 0

% Velocity recalculation -- need to make this more efficient than a loop!
var_norm = svar.phv ./ max(svar.phv,1); % normalize variance
b1 = 6; % aggressiveness of filter
var_inv = (1 - var_norm).^b1; % invert normalized variance
PSD.var = var_inv .* PSD.old;

b2 = 4; % aggressiveness of filter
PSD.agg = (agg.rain.^b2) .* PSD.old;

vr_new.dca = zeros(size(iqh,2), size(iqh,3));
vr_new.var = vr_new.dca;
vr_new.agg = vr_new.dca;
for i = 1:size(iqh,2)
    for j = 1:size(iqh,3)
        vr_temp = calc_velocity_spectral(PSD.dca(:,i,j), M, Ts, lambda);
        if length(vr_temp) > 1
            vr_new.dca(i,j) = vr_old(i,j);
        else
            vr_new.dca(i,j) = vr_temp;
        end
        
        vr_new.var(i,j) = calc_velocity_spectral(PSD.var(:,i,j), M, Ts, lambda);
        vr_new.agg(i,j) = calc_velocity_spectral(PSD.agg(:,i,j), M, Ts, lambda);
    end
end

clear vr


if x == 1
    rain_dir = strrep(data_dir, ['debris' dtype], 'nodebris');
    rain_file = [rain_dir '/' strrep(fname, ['DCU-d' dtype 'n' dnum], 'U-nodebris')];
    rain_sim = load([rain_file '.mat']);
    
    iqh_rain = permute(rain_sim.iqh, [3 1 2]);
    vr_truth = rain_sim.vr;
    PSD.truth = fftshift(abs(fft(iqh_rain.*d, M, 1)).^2, 1);
end


if x == 1
    dv = struct('uncorr', [], 'dca', [], 'var', [], 'agg', []);
    
    dv.uncorr = struct('data', vr_old - vr_truth,...
        'mean', mean(vr_old - vr_truth, 'all'),...
        'prc25', prctile(vr_old - vr_truth, 25, 'all'),...
        'prc75', prctile(vr_old - vr_truth, 75, 'all'));
    
    dv.dca = struct('data', vr_new.dca - vr_truth,...
        'mean', mean(vr_new.dca - vr_truth, 'all'),...
        'prc25', prctile(vr_new.dca - vr_truth, 25, 'all'),...
        'prc75', prctile(vr_new.dca - vr_truth, 75, 'all'));
    
    dv.var = struct('data', vr_new.var - vr_truth,...
        'mean', mean(vr_new.var - vr_truth, 'all'),...
        'prc25', prctile(vr_new.var - vr_truth, 25, 'all'),...
        'prc75', prctile(vr_new.var - vr_truth, 75, 'all'));
    
    dv.agg = struct('data', vr_new.agg - vr_truth,...
        'mean', mean(vr_new.agg - vr_truth, 'all'),...
        'prc25', prctile(vr_new.agg - vr_truth, 25, 'all'),...
        'prc75', prctile(vr_new.agg - vr_truth, 75, 'all'));
    
elseif x == 2
    vr_new.unfolded.dca = vr_new.dca;
    vr_new.unfolded.var = vr_new.var;
    vr_new.unfolded.agg = vr_new.agg;

    i1 = find(vr_new.dca.*vr_unfolded < 0 & vr_unfolded < 0);
    i2 = find(vr_new.dca.*vr_unfolded < 0 & vr_unfolded > 0);
    vr_new.unfolded.dca(i1) = vr_new.dca(i1) - 2*va;
    vr_new.unfolded.dca(i2) = vr_new.dca(i2) + 2*va;
    
    i1 = find(vr_new.var.*vr_unfolded < 0 & vr_unfolded < 0);
    i2 = find(vr_new.var.*vr_unfolded < 0 & vr_unfolded > 0);
    vr_new.unfolded.var(i1) = vr_new.var(i1) - 2*va;
    vr_new.unfolded.var(i2) = vr_new.var(i2) + 2*va;
    
    i1 = find(vr_new.agg.*vr_unfolded < 0 & vr_unfolded < 0);
    i2 = find(vr_new.agg.*vr_unfolded < 0 & vr_unfolded > 0);
    vr_new.unfolded.agg(i1) = vr_new.agg(i1) - 2*va;
    vr_new.unfolded.agg(i2) = vr_new.agg(i2) + 2*va;
    
    dv = struct('unfolded', [], 'dca', vr_new.dca - vr_old,...
        'var', vr_new.var - vr_old, 'agg', vr_new.agg - vr_old);
    dv.unfolded = struct('dca', vr_new.unfolded.dca - vr_unfolded,...
                'var', vr_new.unfolded.var - vr_unfolded,...
                'agg', vr_new.unfolded.agg - vr_unfolded);
end


%% Save variables


if var_save_flag
    if x == 1 %if SimRadar
        save([data_dir '/dpsd-' erase(fname,'sim-'), '.mat'], 'szdr', 'sphv', 'svar',...
            'spsd', 'iqh', 'iqv', 'vr_old', 'vr_new', 'vr_truth', 'dv', 'obj_class',...
            'agg', 'memf', 'vvx', 'PSD', 'params');
        
        % szdr, sphv: [M,r,az]
        %
        % svar: [1x1 struct]
        %     Fields (phv, zdr): [M,r,az]
        %
        % spsd: [1x1 struct] -> Periodogram PSD estimates
        %     Fields (h, v, x): [M,r,az]
        %
        % iqh, iqv: [M,r,az]
        %
        % vr_old, vr_truth: [r,az]
        %
        % vr_new: [1x1 struct]
        %     Fields (dca, var, agg): [r,az]
        %
        % dv: [1x1] struct
        %     Fields (uncorr, dca, var, agg): [1x1 struct]
        %       Subfield (data): [r,az]
        %       Subfields (mean, prc25, prc75): Single float
        %
        % obj_class: [M,r,az] -> 0=rain, 1=debris
        %
        % agg: [1x1 struct]
        %     Fields (debr, rain): [M,r,az]
        %
        % memf: [1x1 struct]
        %     Fields (debr, rain): [1x1 struct]
        %       Subfields (zdr, phv, pvar, zvar): [M,r,az]
        %
        % vvx: [1,M] -> Velocity plotting axis
        %
        % PSD: [1x1 struct] -> Time-series PSD estimates
        %     Fields (old, dca, var, agg, truth): [M,r,az]
        %
        % params: [1x1 struct]
        %     Fields (xx, yy, zz, zh, zdr, phv): [r,az]
        %     Field (va): Single float
        
        
    elseif x == 2 %if KOUN
        save([data_dir '/dpsd-' fname, '.mat'], 'szdr', 'sphv', 'svar', 'spsd',...
            'iqh', 'iqv', 'vr_old', 'vr_new', 'vr_unfolded', 'dv', 'obj_class',...
            'agg', 'memf', 'vvx', 'PSD', 'params');
    end
end


if ~plot_flag
    return
end

%% Plots for SimRadar

if x == 1
    
    if false
        
        % PPIs of truth, original, and corrected velocity for all methods
        figure(1)
        clf
        
        ha = subplot(3,2,[1,2]);
        hs = pcolor(xx, yy, vr_truth(:,:,1));
        caxis([-1 1] * va)
        colormap(ha, blib('rgmap2'))
        shading flat
        c = colorbar;
        c.Label.String = 'm s^{-1}';
        c.Label.FontSize = 13;
        c.Label.VerticalAlignment = 'middle';
        title('(a) v_{truth}', 'FontSize', 14)
        xlabel('x (m)', 'FontSize', 14)
        ylabel('y (m)', 'FontSize', 14)
        ha.Position = [0.3575 0.7093 0.2850 0.2157];
        
        ha(2) = subplot(3,2,3);
        hs(2) = pcolor(xx, yy, vr_old(:,:,1));
        caxis([-1 1] * va)
        colormap(ha(2), blib('rgmap2'))
        shading flat
        c = colorbar;
        c.Label.String = 'm s^{-1}';
        c.Label.FontSize = 13;
        c.Label.VerticalAlignment = 'middle';
        title('(c) Uncorrected v_r', 'FontSize', 14)
        xlabel('x (m)', 'FontSize', 14)
        ylabel('y (m)', 'FontSize', 14)
        
        ha(3) = subplot(3,2,4);
        hs(3) = pcolor(xx, yy, vr_new.dca(:,:,1));
        caxis([-1 1] * va)
        colormap(ha(3), blib('rgmap2'))
        shading flat
        c = colorbar;
        c.Label.String = 'm s^{-1}';
        c.Label.FontSize = 13;
        c.Label.VerticalAlignment = 'middle';
        title('(c) DCA-corrected v_r', 'FontSize', 14)
        xlabel('x (m)', 'FontSize', 14)
        ylabel('y (m)', 'FontSize', 14)
        
        ha(4) = subplot(3,2,5);
        hs(4) = pcolor(xx, yy, vr_new.var(:,:,1));
        caxis([-1 1] * va)
        colormap(ha(4), blib('rgmap2'))
        shading flat
        c = colorbar;
        c.Label.String = 'm s^{-1}';
        c.Label.FontSize = 13;
        c.Label.VerticalAlignment = 'middle';
        title('(d) Variance-corrected v_r', 'FontSize', 14)
        xlabel('x (m)', 'FontSize', 14)
        ylabel('y (m)', 'FontSize', 14)
        
        ha(5) = subplot(3,2,6);
        hs(5) = pcolor(xx, yy, vr_new.agg(:,:,1));
        caxis([-1 1] * va)
        colormap(ha(5), blib('rgmap2'))
        shading flat
        c = colorbar;
        c.Label.String = 'm s^{-1}';
        c.Label.FontSize = 13;
        c.Label.VerticalAlignment = 'middle';
        title('(e) Aggregation-corrected v_r', 'FontSize', 14)
        xlabel('x (m)', 'FontSize', 14)
        ylabel('y (m)', 'FontSize', 14)
        
        set(gcf, 'Units', 'inches', 'Position', [1 1 12 12])
        axes('Unit', 'Normalized', 'Position', [0.5 0.03 0.01 0.01])
        title(['Debris type ' dtype ', n=' dnum ', tilt=' elev '^o'], 'FontSize', 14);
        axis off
        
        if plot_save_flag
            print([img_dir 'vcorr-all-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
        end
        
    end
    
%%
    % PPIs of truth, original, and aggregation-corrected velocity
    figure(3)
    clf
    
    ha = subplot(2,2,[1,2]);
    hs = pcolor(xx, yy, vr_truth(:,:,1));
    caxis([-1 1] * va)
    colormap(ha, blib('rgmap2'))
    shading flat
    c = colorbar;
    c.Label.String = 'm s^{-1}';
    c.Label.FontSize = 12;
    c.Label.VerticalAlignment = 'middle';
    title('Truth v_r', 'FontSize', 12)
    xlabel('x (m)', 'FontSize', 12)
    ylabel('y (m)', 'FontSize', 12)
    ha.Position = [0.3584 0.5838 0.2832 0.3412];
    
    ha(2) = subplot(2,2,4);
    hs(2) = pcolor(xx, yy, vr_new.agg(:,:,1));
    caxis([-1 1] * va)
    colormap(ha(2), blib('rgmap2'))
    shading flat
    c = colorbar;
    c.Label.String = 'm s^{-1}';
    c.Label.FontSize = 12;
    c.Label.VerticalAlignment = 'middle';
    title(['Corrected v_r (10,000 boards)'], 'FontSize', 12)
    xlabel('x (m)', 'FontSize', 12)
    ylabel('y (m)', 'FontSize', 12)
    
    ha(5) = subplot(2,2,3);
    hs(5) = pcolor(xx, yy, vr_old(:,:,1));
    caxis([-1 1] * va)
    colormap(ha(5), blib('rgmap2'))
    shading flat
    c = colorbar;
    c.Label.String = 'm s^{-1}';
    c.Label.FontSize = 12;
    c.Label.VerticalAlignment = 'middle';
    title(['Uncorrected v_r (10,000 boards)'], 'FontSize', 12)
    xlabel('x (m)', 'FontSize', 12)
    ylabel('y (m)', 'FontSize', 12)
    
    set(gcf, 'Units', 'inches', 'Position', [1 1 8 5])
%%
    if plot_save_flag
        print([img_dir 'vppi-agg-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
    end
    

    % Plot error of uncorrected and corrected velocities
    
%     dv.uncorr.data = vr_old - vr_truth;
%     dv.dca.data = vr_new.dca - vr_truth;
%     dv.var.data = vr_new.var - vr_truth;
%     dv.agg.data = vr_new.agg - vr_truth;
%     
%     dv.uncorr.mean = mean(dv.uncorr.data, 'all');
%     dv.dca.mean = mean(dv.dca.data, 'all');
%     dv.var.mean = mean(dv.var.data, 'all');
%     dv.agg.mean = mean(dv.agg.data, 'all');
%     
%     dv.uncorr.prc25 = prctile(dv.uncorr.data, 25, 'all');
%     dv.dca.prc25 = prctile(dv.dca.data, 25, 'all');
%     dv.var.prc25 = prctile(dv.var.data, 25, 'all');
%     dv.agg.prc25 = prctile(dv.agg.data, 25, 'all');
%     
%     dv.uncorr.prc75 = prctile(dv.uncorr.data, 75, 'all');
%     dv.dca.prc75 = prctile(dv.dca.data, 75, 'all');
%     dv.var.prc75 = prctile(dv.var.data, 75, 'all');
%     dv.agg.prc75 = prctile(dv.agg.data, 75, 'all');
    
    if false
        
        % PPI/histograms of DCA-corrected velocity bias
        figure(5)
        clf
        
        ha = subplot(2,2,1);
        hs = pcolor(xx, yy, dv.uncorr.data(:,:,1));
        caxis([-1 1] * va)
        colormap(ha, blib('rbmap'))
        shading flat
        c = colorbar;
        c.Label.String = 'm s^{-1}';
        c.Label.FontSize = 13;
        c.Label.VerticalAlignment = 'middle';
        title('(a) Pre-correction', 'FontSize', 14)
        xlabel('x (m)', 'FontSize', 14)
        ylabel('y (m)', 'FontSize', 14)
        
        ha(2) = subplot(2,2,2);
        hs(2) = pcolor(xx, yy, dv.dca.data(:,:,1));
        caxis([-1 1] * va)
        colormap(ha(2), blib('rbmap'))
        shading flat
        c = colorbar;
        c.Label.String = 'm s^{-1}';
        c.Label.FontSize = 13;
        c.Label.VerticalAlignment = 'middle';
        title('(b) Post-correction', 'FontSize', 14)
        xlabel('x (m)', 'FontSize', 14)
        ylabel('y (m)', 'FontSize', 14)
        
        ha(3) = subplot(2,2,3);
        histogram(dv.uncorr.data, -va:5:va)
        hold on
        plot(ones(1,2)*dv.uncorr.mean, [0 350], 'k', 'LineWidth', 1.2)
        plot(ones(1,2)*dv.uncorr.prc25, [0 350], '--k', 'LineWidth', 1.2)
        plot(ones(1,2)*dv.uncorr.prc75, [0 350], '--k', 'LineWidth', 1.2)
        hold off
        title('(c)', 'FontSize', 14)
        xlabel('v_{bias} (m s^{-1})', 'FontSize', 14)
        ylim([0 350])
        text(-100, 330, ['\mu = ' num2str(double(dv.uncorr.mean),'%.2f')], 'FontSize', 13)
        text(-100, 300, ['Q_1 = ' num2str(double(dv.uncorr.prc25),'%.2f')], 'FontSize', 13)
        text(-100, 270, ['Q_3 = ' num2str(double(dv.uncorr.prc75),'%.2f')], 'FontSize', 13)
        
        ha(4) = subplot(2,2,4);
        histogram(dv.dca.data, -va:5:va)
        hold on
        plot(ones(1,2)*dv.dca.mean, [0 350], 'k', 'LineWidth', 1.2)
        plot(ones(1,2)*dv.dca.prc25, [0 350], '--k', 'LineWidth', 1.2)
        plot(ones(1,2)*dv.dca.prc75, [0 350], '--k', 'LineWidth', 1.2)
        hold off
        title('(d)', 'FontSize', 14)
        xlabel('v_{bias} (m s^{-1})', 'FontSize', 14)
        ylim([0 350])
        text(-100, 330, ['\mu = ' num2str(double(dv.dca.mean),'%.2f')], 'FontSize', 13)
        text(-100, 300, ['Q_1 = ' num2str(double(dv.dca.prc25),'%.2f')], 'FontSize', 13)
        text(-100, 270, ['Q_3 = ' num2str(double(dv.dca.prc75),'%.2f')], 'FontSize', 13)
        
        set(gcf, 'Units', 'inches', 'Position', [10 10 12 8])
        axes('Unit', 'Normalized', 'Position', [0.5 0.48 0.01 0.01])
        title(['Debris type ' dtype ', n=' dnum ', tilt=' elev '^o'], 'FontSize', 14);
        axis off
        
        %sgtitle('v_{bias} - Debris-removed correction','FontSize',16);
        
        if plot_save_flag
            print([img_dir 'vbias-DCA-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
        end
        
        
        % PPI/histograms of variance-corrected velocity bias
        figure(6)
        clf
        
        ha = subplot(2,2,1);
        hs = pcolor(xx, yy, dv.uncorr.data(:,:,1));
        caxis([-1 1] * va)
        colormap(ha, blib('rbmap'))
        shading flat
        c = colorbar;
        c.Label.String = 'm s^{-1}';
        c.Label.FontSize = 13;
        c.Label.VerticalAlignment = 'middle';
        title('(a) Pre-correction', 'FontSize', 14)
        xlabel('x (m)', 'FontSize', 14)
        ylabel('y (m)', 'FontSize', 14)
        
        ha(2) = subplot(2,2,2);
        hs(2) = pcolor(xx, yy, dv.var.data(:,:,1));
        caxis([-1 1] * va)
        colormap(ha(2), blib('rbmap'))
        shading flat
        c = colorbar;
        c.Label.String = 'm s^{-1}';
        c.Label.FontSize = 13;
        c.Label.VerticalAlignment = 'middle';
        title('(b) Post-correction', 'FontSize', 14)
        xlabel('x (m)', 'FontSize', 14)
        ylabel('y (m)', 'FontSize', 14)
        
        ha(3) = subplot(2,2,3);
        histogram(dv.uncorr.data, -va:5:va)
        hold on
        plot(ones(1,2)*dv.uncorr.mean, [0 350], 'k', 'LineWidth', 1.2)
        plot(ones(1,2)*dv.uncorr.prc25, [0 350], '--k', 'LineWidth', 1.2)
        plot(ones(1,2)*dv.uncorr.prc75, [0 350], '--k', 'LineWidth', 1.2)
        hold off
        title('(c)', 'FontSize', 14)
        xlabel('v_{bias} (m s^{-1})', 'FontSize', 14)
        ylim([0 350])
        text(-100, 330, ['\mu = ' num2str(double(dv.uncorr.mean),'%.2f')], 'FontSize', 13)
        text(-100, 300, ['Q_1 = ' num2str(double(dv.uncorr.prc25),'%.2f')], 'FontSize', 13)
        text(-100, 270, ['Q_3 = ' num2str(double(dv.uncorr.prc75),'%.2f')], 'FontSize', 13)
        
        ha(4) = subplot(2,2,4);
        histogram(dv.var.data, -va:5:va)
        hold on
        plot(ones(1,2)*dv.var.mean, [0 350], 'k', 'LineWidth', 1.2)
        plot(ones(1,2)*dv.var.prc25, [0 350], '--k', 'LineWidth', 1.2)
        plot(ones(1,2)*dv.var.prc75, [0 350], '--k', 'LineWidth', 1.2)
        hold off
        title('(d)', 'FontSize', 14)
        xlabel('v_{bias} (m s^{-1})', 'FontSize', 14)
        ylim([0 350])
        text(-100, 330, ['\mu = ' num2str(double(dv.var.mean),'%.2f')], 'FontSize', 13)
        text(-100, 300, ['Q_1 = ' num2str(double(dv.var.prc25),'%.2f')], 'FontSize', 13)
        text(-100, 270, ['Q_3 = ' num2str(double(dv.var.prc75),'%.2f')], 'FontSize', 13)
        
        set(gcf, 'Units', 'inches', 'Position', [10 10 12 8])
        axes('Unit', 'Normalized', 'Position', [0.5 0.48 0.01 0.01])
        title(['Debris type ' dtype ', n=' dnum ', tilt=' elev '^o'], 'FontSize', 14);
        axis off
        
        %sgtitle('v_{bias} - Variance-weighted correction','FontSize',16);
        
        if plot_save_flag
            print([img_dir 'vbias-var-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
        end
        
        
        % PPI/histograms of aggregation-corrected velocity bias
        figure(7)
        clf
        
        ha = subplot(2,2,1);
        hs = pcolor(xx, yy, dv.uncorr.data(:,:,1));
        caxis([-1 1] * va)
        colormap(ha, blib('rbmap'))
        shading flat
        c = colorbar;
        c.Label.String = 'm s^{-1}';
        c.Label.FontSize = 13;
        c.Label.VerticalAlignment = 'middle';
        %title('(a) Pre-correction', 'FontSize', 14)
        title('Pre-Correction Bias', 'FontSize', 16)
        xlabel('x (m)', 'FontSize', 14)
        ylabel('y (m)', 'FontSize', 14)
        
        ha(2) = subplot(2,2,2);
        hs(2) = pcolor(xx, yy, dv.agg.data(:,:,1));
        caxis([-1 1] * va)
        colormap(ha(2), blib('rbmap'))
        shading flat
        c = colorbar;
        c.Label.String = 'm s^{-1}';
        c.Label.FontSize = 13;
        c.Label.VerticalAlignment = 'middle';
        %title('(b) Post-correction', 'FontSize', 14)
        title('Post-Correction Bias', 'FontSize', 16)
        xlabel('x (m)', 'FontSize', 14)
        ylabel('y (m)', 'FontSize', 14)
        
        ha(3) = subplot(2,2,3);
        histogram(dv.uncorr.data, -va:5:va)
        hold on
        plot(ones(1,2)*dv.uncorr.mean, [0 350], 'k', 'LineWidth', 1.2)
        plot(ones(1,2)*dv.uncorr.prc25, [0 350], '--k', 'LineWidth', 1.2)
        plot(ones(1,2)*dv.uncorr.prc75, [0 350], '--k', 'LineWidth', 1.2)
        hold off
        %title('(c)', 'FontSize', 14)
        xlabel('v_{bias} (m s^{-1})', 'FontSize', 14)
        ylim([0 350])
        text(-100, 330, ['\mu = ' num2str(double(dv.uncorr.mean),'%.2f')], 'FontSize', 13)
        text(-100, 300, ['Q_1 = ' num2str(double(dv.uncorr.prc25),'%.2f')], 'FontSize', 13)
        text(-100, 270, ['Q_3 = ' num2str(double(dv.uncorr.prc75),'%.2f')], 'FontSize', 13)
        
        ha(4) = subplot(2,2,4);
        histogram(dv.agg.data, -va:5:va)
        hold on
        plot(ones(1,2)*dv.agg.mean, [0 350], 'k', 'LineWidth', 1.2)
        plot(ones(1,2)*dv.agg.prc25, [0 350], '--k', 'LineWidth', 1.2)
        plot(ones(1,2)*dv.agg.prc75, [0 350], '--k', 'LineWidth', 1.2)
        hold off
        %title('(d)', 'FontSize', 14)
        xlabel('v_{bias} (m s^{-1})', 'FontSize', 14)
        ylim([0 350])
        text(-100, 330, ['\mu = ' num2str(double(dv.agg.mean),'%.2f')], 'FontSize', 13)
        text(-100, 300, ['Q_1 = ' num2str(double(dv.agg.prc25),'%.2f')], 'FontSize', 13)
        text(-100, 270, ['Q_3 = ' num2str(double(dv.agg.prc75),'%.2f')], 'FontSize', 13)
        
        set(gcf, 'Units', 'inches', 'Position', [10 10 13 8])
        axes('Unit', 'Normalized', 'Position', [0.5 0.48 0.01 0.01])
        title(['Debris type ' dtype ', n=' dnum ', tilt=' elev '^o'], 'FontSize', 14);
        axis off
        
        %sgtitle('v_{bias} - Aggregation-weighted correction','FontSize',16);
        
        if plot_save_flag
            print([img_dir 'vbias-agg-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
        end
        
    end
    
    
    % PPI/histograms of aggregation-corrected velocity bias, but prettier?
    figure(8)
    clf
    
    %ha = subplot(4,2,1);
    ha = subplot(2,2,1);
    hs = pcolor(xx, yy, dv.uncorr.data(:,:,1));
    caxis([-1 1] * va)
    colormap(ha, blib('rbmap'))
    shading flat
    c = colorbar;
    c.Label.String = 'm s^{-1}';
    c.Label.FontSize = 13;
    c.Label.VerticalAlignment = 'middle';
    %title('(a) v_{bias}', 'FontSize', 16)
    title('Velocity Bias', 'FontSize', 16)
    xlabel('x (m)', 'FontSize', 14)
    ylabel('y (m)', 'FontSize', 14)
    
    %subplot(4,2,2)
    subplot(2,2,2)
    histogram(dv.uncorr.data, -va:5:va)
    hold on
    plot(ones(1,2)*dv.uncorr.mean, [0 350], 'k', 'LineWidth', 1.2)
    plot(ones(1,2)*dv.uncorr.prc25, [0 350], '--k', 'LineWidth', 1.2)
    plot(ones(1,2)*dv.uncorr.prc75, [0 350], '--k', 'LineWidth', 1.2)
    hold off
    %title('(b) v_{bias} distribution', 'FontSize', 16)
    title('Velocity Bias Distribution', 'FontSize', 16)
    %xlabel('v_{bias} (m s^{-1})', 'FontSize', 14)
    ylim([0 350])
    text(-100, 330, ['\mu = ' num2str(double(dv.uncorr.mean),'%.2f')], 'FontSize', 13)
    text(-100, 290, ['Q_1 = ' num2str(double(dv.uncorr.prc25),'%.2f')], 'FontSize', 13)
    text(-100, 250, ['Q_3 = ' num2str(double(dv.uncorr.prc75),'%.2f')], 'FontSize', 13)
    
    
    %ha(4) = subplot(4,2,7);
    ha(4) = subplot(2,2,3);
    hs(4) = pcolor(xx, yy, dv.agg.data(:,:,1));
    caxis([-1 1] * va)
    colormap(ha(4), blib('rbmap'))
    shading flat
    c = colorbar;
    c.Label.String = 'm s^{-1}';
    c.Label.FontSize = 13;
    c.Label.VerticalAlignment = 'middle';
    %title('(g) v_{bias}', 'FontSize', 14)
    %title('(g)', 'FontSize', 16)
    xlabel('x (m)', 'FontSize', 14)
    ylabel('y (m)', 'FontSize', 14)
    
    %subplot(4,2,8)
    subplot(2,2,4)
    histogram(dv.agg.data, -va:5:va)
    hold on
    plot(ones(1,2)*dv.agg.mean, [0 350], 'k', 'LineWidth', 1.2)
    plot(ones(1,2)*dv.agg.prc25, [0 350], '--k', 'LineWidth', 1.2)
    plot(ones(1,2)*dv.agg.prc75, [0 350], '--k', 'LineWidth', 1.2)
    hold off
    %title('(h) v_{bias} distribution', 'FontSize', 14)
    %title('(h)', 'FontSize', 16)
    xlabel('v_{bias} (m s^{-1})', 'FontSize', 14)
    ylim([0 350])
    text(-100, 330, ['\mu = ' num2str(double(dv.agg.mean),'%.2f')], 'FontSize', 13)
    text(-100, 290, ['Q_1 = ' num2str(double(dv.agg.prc25),'%.2f')], 'FontSize', 13)
    text(-100, 250, ['Q_3 = ' num2str(double(dv.agg.prc75),'%.2f')], 'FontSize', 13)
    
    %set(gcf, 'Units', 'inches', 'Position', [10 10 12 14])
    set(gcf, 'Units', 'inches', 'Position', [10 10 13 8.3])
    % axes('Unit', 'Normalized', 'Position', [0.5 0.03 0.01 0.01])
    % title(['Debris type ' dtype ', n=' dnum ', tilt=' elev '^o'], 'FontSize', 14);
    % axis off
    
    % FUCKING HALLELUJAH
    %annotation('textbox', [0.43 0.91 0.35 0.05], 'String', 'No correction', 'FontSize', 16, 'EdgeColor', 'none')
    %annotation('textbox', [0.4 0.69 0.35 0.05], 'String', 'DCA-based correction', 'FontSize', 16, 'EdgeColor', 'none')
    %annotation('textbox', [0.38 0.47 0.35 0.05], 'String', 'Variance-based correction', 'FontSize', 16, 'EdgeColor', 'none')
    %annotation('textbox', [0.37 0.25 0.35 0.05], 'String', 'Aggregation-based correction', 'FontSize', 16, 'EdgeColor', 'none')
    annotation('textbox', [0.41 0.94 0.16 0.05], 'String', 'Uncorrected',...
        'FontSize',18,'FontWeight','bold','EdgeColor','k','LineWidth',1,...
        'HorizontalAlignment','center','VerticalAlignment','middle')
    annotation('textbox', [0.34 0.47 0.32 0.05], 'String', 'Aggregation-Based Correction',...
        'FontSize',18,'FontWeight','bold','EdgeColor','k','LineWidth',1,...
        'HorizontalAlignment','center','VerticalAlignment','middle')
    
    if plot_save_flag
        print([img_dir 'vcorr-agg-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
        %print([img_dir 'vcorr-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
    end
    
    
    
    
    
    rind = 15;
    azind = 21;
    
    if false
        
        % Plot of sZDR, sPHV, and DCA classification w/ PSD underlaid
        figure(10)
        clf
        
        ax1 = subplot(3,1,1);
        % DPSD
        yyaxis left
        plot(vvx(:), squeeze(szdr(:,rind,azind)), 'k', 'LineWidth', 1.1)
        ylabel('sZ_{DR} (dB)', 'FontSize', 14)
        xlim([-va, va])
        ylim([-15, 15])
        % PSD
        yyaxis right
        hs1 = semilogy(vvx(:), squeeze(abs(spsd.h(:,rind,azind))), 'LineWidth', 1);
        hs1.Color = [0.6 0.6 0.6];
        hs1.LineStyle = '-.';
        ylabel('sS_H', 'FontSize', 14)
        ax1.YAxis(1).Color = 'k';
        ax1.YAxis(2).Color = [0.6 0.6 0.6];
        grid on
        title('(a) sZ_{DR}, sS_H', 'FontSize', 14)
        xlabel('v (m s^{-1})', 'FontSize', 14)
        %legend('sZ_{DR}', 'sS_V', 'Location', 'northwest')
        % text(18, 1, 'Debris', 'FontSize', 20)
        % text(-30, 1, 'Rain', 'FontSize', 20)
        % text(-48, 3.4, 'sZ_{DR}', 'FontSize', 14)
        % text(-48, -3.5, 'sS_V', 'FontSize', 14, 'Color', [0.6 0.6 0.6])
        
        
        ax2 = subplot(3,1,2);
        % DPSD
        yyaxis left
        plot(vvx(:), squeeze(sphv(:,rind,azind)), 'k', 'LineWidth', 1.1)
        ylabel('s\rho_{HV}', 'FontSize', 14)
        xlim([-va, va])
        ylim([0, 1])
        % PSD
        yyaxis right
        hs2 = semilogy(vvx(:), squeeze(abs(spsd.h(:,rind,azind))), 'LineWidth', 1);
        hs2.Color = [0.6 0.6 0.6];
        hs2.LineStyle = '-.';
        ylabel('sS_H', 'FontSize', 14)
        ax2.YAxis(1).Color = 'k';
        ax2.YAxis(2).Color = [0.6 0.6 0.6];
        grid on
        title('(b) s\rho_{HV}, sS_H', 'FontSize', 14)
        xlabel('v (m s^{-1})', 'FontSize', 14)
        % text(18, 0.6, 'Debris', 'FontSize', 20)
        % text(-30, 0.6, 'Rain', 'FontSize', 20)
        % text(-48, 0.95, 's\rho_{HV}', 'FontSize', 14)
        % text(-48, 0.15, 'sS_V', 'FontSize', 14, 'Color', [0.6 0.6 0.6])
        
        
        ax3 = subplot(3,1,3);
        plot(vvx(:), squeeze(obj_class(:,rind,azind)), 'k', 'LineWidth', 1.1)
        grid on
        title('(c) Scatterer classification', 'FontSize', 14)
        xlabel('v (m s^{-1})', 'FontSize', 14)
        xlim([-va, va])
        ylim([-0.5,1.5])
        yticks([0 1])
        yticklabels({'Rain', 'Debris'})
        %text(-48, -1.8, 'sZ_{DR}', 'FontSize', 14)
        %text(-47, -13.5, 'sS_H', 'FontSize', 14, 'Color', [0.6 0.6 0.6])
        
        if plot_save_flag
            print([img_dir 'hcaclass-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
        end
        
    end
    
    
    % Plot of DCA classification w/ original and filtered PSDs + velocities
    figure(11)
    clf
    
    subplot(3,1,1)
    plot(vvx(:), squeeze(obj_class(:,rind,azind)), 'k', 'LineWidth', 1.1)
    %title('(a) DCA classification', 'FontSize', 14)
    title('DCA Classification', 'FontSize', 18)
    %xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
    xlim([-va, va])
    ylim([-0.5,1.5])
    yticks([0 1])
    yticklabels([])
    annotation('textbox', [0.05 0.83 0.1 0.05], 'String', 'Debris', 'FontSize', 15, 'EdgeColor', 'none')
    annotation('textbox', [0.07 0.73 0.1 0.05], 'String', 'Rain', 'FontSize', 15, 'EdgeColor', 'none')
    
    subplot(3,1,2)
    semilogy(vvx(:), squeeze(abs(PSD.old(:,rind,azind))), 'k', 'LineWidth', 1.1)
    hold on
    semilogy(ones(1,2)*vr_old(rind,azind), [1e0 1e12], ':b', 'LineWidth', 2.1)
    semilogy(ones(1,2)*vr_truth(rind,azind), [1e0 1e12], ':k', 'LineWidth', 2.1)
    hold off
    %xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
    %ylabel('sS_H', 'FontSize', 14)
    xlim([-va, va])
    ylim([1e0 1e12])
    %title('(b) Original sS_H(v) and v_r', 'FontSize', 14)
    title('Original Spectrum', 'FontSize', 18)
    if vr_old(rind,azind) >= vr_truth(rind,azind)
        text_pos1 = double(vr_old(rind,azind)+3);
        text_pos2 = double(vr_truth(rind,azind)-51);
    else
        text_pos1 = double(vr_old(rind,azind)-51);
        text_pos2 = double(vr_truth(rind,azind)+3);
    end
    text(text_pos1, 1e2, ['v_{old} = ' num2str(round(vr_old(rind,azind),1)) ' m s^{-1}'], 'FontSize', 15, 'Color', 'b')
    text(text_pos2, 1e2, ['v_{truth} = ' num2str(round(vr_truth(rind,azind),1)) ' m s^{-1}'], 'FontSize', 15)
    
    subplot(3,1,3)
    semilogy(vvx(:), squeeze(abs(PSD.dca(:,rind,azind))), 'k', 'LineWidth', 1.5)
    hold on
    semilogy(ones(1,2)*vr_new.dca(rind,azind), [1e0 1e12], ':r', 'LineWidth', 2.1)
    semilogy(ones(1,2)*vr_truth(rind,azind), [1e0 1e12], ':k', 'LineWidth', 2.1)
    hold off
    xlabel('{\it v} (m s^{-1})', 'FontSize', 16)
    %ylabel('sS_H', 'FontSize', 14)
    xlim([-va, va])
    ylim([1e0 1e12])
    %title('(c) Filtered sS_H(v) and corrected v_r', 'FontSize', 14)
    title('Filtered Spectrum', 'FontSize', 18)
    if vr_new.dca(rind,azind) >= vr_truth(rind,azind)
        text_pos1 = double(vr_new.dca(rind,azind)+3);
        text_pos2 = double(vr_truth(rind,azind)-51);
    else
        text_pos1 = double(vr_new.dca(rind,azind)-51);
        text_pos2 = double(vr_truth(rind,azind)+3);
    end
    text(text_pos1, 1e2, ['v_{new} = ' num2str(round(vr_new.dca(rind,azind),1)) ' m s^{-1}'], 'FontSize', 15, 'Color', [0.8 0.2 0.2])
    text(text_pos2, 1e2, ['v_{truth} = ' num2str(round(vr_truth(rind,azind),1)) ' m s^{-1}'], 'FontSize', 15)
    
    set(gcf, 'Units', 'inches', 'Position', [6 8 9 7.5])
    
    if plot_save_flag
        print([img_dir 'vnew-DCA-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
    end
    
    
    
    % Plot of sPHV variance w/ original and filtered PSDs + velocities
    figure(12)
    clf
    
    subplot(3,1,1)
    plot(vvx(:), squeeze(svar.phv(:,rind,azind)), 'k', 'LineWidth', 1.1)
    %title('(a) \sigma^2_{s\rho_{HV}}', 'FontSize', 14)
    title('s\rho_{HV} Variance', 'FontSize', 18)
    %xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
    xlim([-va, va])
    ylim([0, 0.3])
    
    
    subplot(3,1,2)
    semilogy(vvx(:), squeeze(abs(PSD.old(:,rind,azind))), 'k', 'LineWidth', 1.1)
    hold on
    semilogy(ones(1,2)*vr_old(rind,azind), [1e0 1e12], ':b', 'LineWidth', 2.1)
    semilogy(ones(1,2)*vr_truth(rind,azind), [1e0 1e12], ':k', 'LineWidth', 2.1)
    hold off
    %xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
    %ylabel('sS_H', 'FontSize', 14)
    xlim([-va, va])
    ylim([1e0 1e12])
    %title('(b) Original sS_H(v) and v_r', 'FontSize', 14)
    title('Original Spectrum', 'FontSize', 18)
    if vr_old(rind,azind) >= vr_truth(rind,azind)
        text_pos1 = double(vr_old(rind,azind)+3);
        text_pos2 = double(vr_truth(rind,azind)-51);
    else
        text_pos1 = double(vr_old(rind,azind)-51);
        text_pos2 = double(vr_truth(rind,azind)+3);
    end
    text(text_pos1, 1e2, ['v_{old} = ' num2str(round(vr_old(rind,azind),1)) ' m s^{-1}'], 'FontSize', 15, 'Color', 'b')
    text(text_pos2, 1e2, ['v_{truth} = ' num2str(round(vr_truth(rind,azind),1)) ' m s^{-1}'], 'FontSize', 15)
    
    
    subplot(3,1,3)
    semilogy(vvx(:), squeeze(abs(PSD.var(:,rind,azind))), 'k', 'LineWidth', 1.1)
    hold on
    semilogy(ones(1,2)*vr_new.var(rind,azind), [1e0 1e12], ':r', 'LineWidth', 2.1)
    semilogy(ones(1,2)*vr_truth(rind,azind), [1e0 1e12], ':k', 'LineWidth', 2.1)
    hold off
    xlabel('{\it v} (m s^{-1})', 'FontSize', 16)
    %ylabel('sS_H', 'FontSize', 14)
    xlim([-va, va])
    ylim([1e0 1e12])
    %title('(c) Filtered sS_H(v) and corrected v_r', 'FontSize', 14)
    title('Filtered Spectrum', 'FontSize', 18)
    if vr_new.var(rind,azind) >= vr_truth(rind,azind)
        text_pos1 = double(vr_new.var(rind,azind)+3);
        text_pos2 = double(vr_truth(rind,azind)-51);
    else
        text_pos1 = double(vr_new.var(rind,azind)-51);
        text_pos2 = double(vr_truth(rind,azind)+3);
    end
    text(text_pos1, 1e2, ['v_{new} = ' num2str(round(vr_new.var(rind,azind),1)) ' m s^{-1}'], 'FontSize', 15, 'Color', [0.8 0.2 0.2])
    text(text_pos2, 1e2, ['v_{truth} = ' num2str(round(vr_truth(rind,azind),1)) ' m s^{-1}'], 'FontSize', 15)
    
    set(gcf, 'Units', 'inches', 'Position', [6 8 9 7.5])
    
    if plot_save_flag
        print([img_dir 'vnew-var-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
    end
    
    
    
    % Plot of rain aggregation w/ original and filtered PSDs + velocities
    figure(13)
    clf
    
    subplot(3,1,1)
    plot(vvx(:), squeeze(agg.rain(:,rind,azind)), 'k', 'LineWidth', 1.1)
    %title('(a) A_{rain}', 'FontSize', 14)
    title('Rain Aggregation Value', 'FontSize', 18)
    %xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
    xlim([-va, va])
    ylim([0, 1])
    
    
    subplot(3,1,2)
    semilogy(vvx(:), squeeze(abs(PSD.old(:,rind,azind))), 'k', 'LineWidth', 1.1)
    hold on
    semilogy(ones(1,2)*vr_old(rind,azind), [1e0 1e12], ':b', 'LineWidth', 2.1)
    semilogy(ones(1,2)*vr_truth(rind,azind), [1e0 1e12], ':k', 'LineWidth', 2.1)
    hold off
    %xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
    %ylabel('sS_H', 'FontSize', 14)
    xlim([-va, va])
    ylim([1e0 1e12])
    %title('(b) Original sS_H(v) and v_r', 'FontSize', 14)
    title('Original Spectrum', 'FontSize', 18)
    if vr_old(rind,azind) >= vr_truth(rind,azind)
        text_pos1 = double(vr_old(rind,azind)+3);
        text_pos2 = double(vr_truth(rind,azind)-51);
    else
        text_pos1 = double(vr_old(rind,azind)-51);
        text_pos2 = double(vr_truth(rind,azind)+3);
    end
    text(text_pos1, 1e2, ['v_{old} = ' num2str(round(vr_old(rind,azind),1)) ' m s^{-1}'], 'FontSize', 15, 'Color', 'b')
    text(text_pos2, 1e2, ['v_{truth} = ' num2str(round(vr_truth(rind,azind),1)) ' m s^{-1}'], 'FontSize', 15)
    
    
    subplot(3,1,3)
    semilogy(vvx(:), squeeze(abs(PSD.agg(:,rind,azind))), 'k', 'LineWidth', 1.1)
    hold on
    semilogy(ones(1,2)*vr_new.agg(rind,azind), [1e0 1e12], ':r', 'LineWidth', 2.1)
    semilogy(ones(1,2)*vr_truth(rind,azind), [1e0 1e12], ':k', 'LineWidth', 2.1)
    hold off
    xlabel('{\it v} (m s^{-1})', 'FontSize', 16)
    %ylabel('sS_H', 'FontSize', 14)
    xlim([-va, va])
    ylim([1e0 1e12])
    %title('(c) Filtered sS_H(v) and corrected v_r', 'FontSize', 14)
    title('Filtered Spectrum', 'FontSize', 18)
    if vr_new.agg(rind,azind) >= vr_truth(rind,azind)
        text_pos1 = double(vr_new.agg(rind,azind)+3);
        text_pos2 = double(vr_truth(rind,azind)-51);
    else
        text_pos1 = double(vr_new.agg(rind,azind)-51);
        text_pos2 = double(vr_truth(rind,azind)+3);
    end
    text(text_pos1, 10^10.5, ['v_{new} = ' num2str(round(vr_new.agg(rind,azind),1)) ' m s^{-1}'], 'FontSize', 15, 'Color', [0.8 0.2 0.2])
    text(text_pos2, 10^10.5, ['v_{truth} = ' num2str(round(vr_truth(rind,azind),1)) ' m s^{-1}'], 'FontSize', 15)
    
    set(gcf, 'Units', 'inches', 'Position', [6 8 9 7.5])
    
    if plot_save_flag
        print([img_dir 'vnew-agg-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
    end
    
    
    if false
        
        % Plots of sZDR + rain and debris membership
        figure(14)
        clf
        
        ax = subplot(3,1,1);
        % DPSD
        yyaxis left
        plot(vvx(:), squeeze(szdr(:,rind,azind)), 'k', 'LineWidth', 1.1)
        ylabel('sZ_{DR} (dB)', 'FontSize', 14)
        xlim([-va, va])
        ylim([-15, 15])
        % PSD
        yyaxis right
        hs = semilogy(vvx(:), squeeze(abs(spsd.h(:,rind,azind))), 'LineWidth', 1);
        hs.Color = [0.6 0.6 0.6];
        hs.LineStyle = '-.';
        ylabel('sS_H', 'FontSize', 14)
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = [0.6 0.6 0.6];
        grid on
        title('(a) sZ_{DR}, sS_H', 'FontSize', 14)
        %xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
        
        subplot(3,1,2)
        plot(vvx(:), squeeze(memf.rain.zdr(:,rind,azind)), 'k', 'LineWidth', 1)
        xlim([-va va])
        ylim([-0.2 1.2])
        yticks(0:0.2:1)
        grid on
        title('(b) sZ_{DR} rain membership')
        %xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
        
        subplot(3,1,3)
        plot(vvx(:), squeeze(memf.debr.zdr(:,rind,azind)), 'k', 'LineWidth', 1)
        xlim([-va va])
        ylim([-0.2 1.2])
        yticks(0:0.2:1)
        grid on
        title('(c) sZ_{DR} debris membership')
        xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
        
        if plot_save_flag
            print([img_dir 'memf-szdr-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
        end
        
        
        % Plots of sPHV + rain and debris membership
        figure(15)
        clf
        
        ax = subplot(3,1,1);
        % DPSD
        yyaxis left
        plot(vvx(:), squeeze(sphv(:,rind,azind)), 'k', 'LineWidth', 1.1)
        ylabel('s\rho_{HV}', 'FontSize', 14)
        xlim([-va, va])
        ylim([0, 1])
        % PSD
        yyaxis right
        hs = semilogy(vvx(:), squeeze(abs(spsd.h(:,rind,azind))), 'LineWidth', 1);
        hs.Color = [0.6 0.6 0.6];
        hs.LineStyle = '-.';
        ylabel('sS_H', 'FontSize', 14)
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = [0.6 0.6 0.6];
        grid on
        title('(a) s\rho_{HV}, sS_H', 'FontSize', 14)
        %xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
        
        subplot(3,1,2)
        plot(vvx(:), squeeze(memf.rain.phv(:,rind,azind)), 'k', 'LineWidth', 1)
        xlim([-va va])
        ylim([-0.2 1.2])
        yticks(0:0.2:1)
        grid on
        title('(b) s\rho_{HV} rain membership')
        %xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
        
        subplot(3,1,3)
        plot(vvx(:), squeeze(memf.debr.phv(:,rind,azind)), 'k', 'LineWidth', 1)
        xlim([-va va])
        ylim([-0.2 1.2])
        yticks(0:0.2:1)
        grid on
        title('(c) s\rho_{HV} debris membership')
        xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
        
        if plot_save_flag
            print([img_dir 'memf-sphv-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
        end
        
        
        % Plots of sPHV variance + rain and debris membership
        figure(16)
        clf
        
        ax = subplot(3,1,1);
        % DPSD
        yyaxis left
        plot(vvx(:), squeeze(svar.phv(:,rind,azind)), 'k', 'LineWidth', 1.1)
        ylabel('\sigma^2_{s\rho_{HV}}', 'FontSize', 14)
        xlim([-va, va])
        ylim([0, 0.3])
        % PSD
        yyaxis right
        hs = semilogy(vvx(:), squeeze(abs(spsd.h(:,rind,azind))), 'LineWidth', 1);
        hs.Color = [0.6 0.6 0.6];
        hs.LineStyle = '-.';
        ylabel('sS_H', 'FontSize', 14)
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = [0.6 0.6 0.6];
        grid on
        title('(a) \sigma^2_{s\rho_{HV}}, sS_H', 'FontSize', 14)
        %xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
        
        subplot(3,1,2)
        plot(vvx(:), squeeze(memf.rain.pvar(:,rind,azind)), 'k', 'LineWidth', 1)
        xlim([-va va])
        ylim([-0.2 1.2])
        yticks(0:0.2:1)
        grid on
        title('(b) \sigma^2_{s\rho_{HV}} rain membership')
        %xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
        
        subplot(3,1,3)
        plot(vvx(:), squeeze(memf.debr.pvar(:,rind,azind)), 'k', 'LineWidth', 1)
        xlim([-va va])
        ylim([-0.2 1.2])
        yticks(0:0.2:1)
        grid on
        title('(c) \sigma^2_{s\rho_{HV}} debris membership')
        xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
        
        if plot_save_flag
            print([img_dir 'memf-pvar-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
        end
        
        
        % IMPORTANT!
        % Plots of filtered PSDs + truth/original/corrected velocities
        figure(17)
        
        subplot(3,1,1)
        semilogy(ones(1,2)*vr_truth(rind,azind), [1e0 1e15], ':k', 'LineWidth', 2)
        hold on
        semilogy(ones(1,2)*vr_new.dca(rind,azind), [1e0 1e15], ':r', 'LineWidth', 2)
        semilogy(ones(1,2)*vr_old(rind,azind), [1e0 1e15], ':b', 'LineWidth', 2)
        semilogy(vvx(:), squeeze(abs(PSD.dca(:,rind,azind))), 'k', 'LineWidth', 2)
        hold off
        %xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
        ylabel('sS_H', 'FontSize', 14)
        xlim([-va, va])
        ylim([1e0 1e15])
        title('(a) DCA-based correction', 'FontSize', 14)
        legend(['v_{truth} = ' num2str(round(vr_truth(rind,azind),1))],...
            ['v_{new} = ' num2str(round(vr_new.dca(rind,azind),1))],...
            ['v_{old} = ' num2str(round(vr_old(rind,azind),1))], 'Location', 'northwest')
        
        
        subplot(3,1,2)
        semilogy(ones(1,2)*vr_truth(rind,azind), [1e0 1e15], ':k', 'LineWidth', 2)
        hold on
        semilogy(ones(1,2)*vr_new.var(rind,azind), [1e0 1e15], ':r', 'LineWidth', 2)
        semilogy(ones(1,2)*vr_old(rind,azind), [1e0 1e15], ':b', 'LineWidth', 2)
        semilogy(vvx(:), squeeze(abs(PSD.var(:,rind,azind))), 'k', 'LineWidth', 1)
        hold off
        %xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
        ylabel('sS_H', 'FontSize', 14)
        xlim([-va, va])
        ylim([1e0 1e15])
        title('(b) Variance-based correction', 'FontSize', 14)
        legend(['v_{truth} = ' num2str(round(vr_truth(rind,azind),1))],...
            ['v_{new} = ' num2str(round(vr_new.var(rind,azind),1))],...
            ['v_{old} = ' num2str(round(vr_old(rind,azind),1))], 'Location', 'northwest')
        
        subplot(3,1,3)
        semilogy(ones(1,2)*vr_truth(rind,azind), [1e0 1e15], ':k', 'LineWidth', 2)
        hold on
        semilogy(ones(1,2)*vr_new.agg(rind,azind), [1e0 1e15], ':r', 'LineWidth', 2)
        semilogy(ones(1,2)*vr_old(rind,azind), [1e0 1e15], ':b', 'LineWidth', 2)
        semilogy(vvx(:), squeeze(abs(PSD.agg(:,rind,azind))), 'k', 'LineWidth', 1)
        hold off
        xlabel('{\it v} (m s^{-1})', 'FontSize', 14)
        ylabel('sS_H', 'FontSize', 14)
        xlim([-va, va])
        ylim([1e0 1e15])
        title('(c) Aggregation-based correction', 'FontSize', 14)
        legend(['v_{truth} = ' num2str(round(vr_truth(rind,azind),1))],...
            ['v_{new} = ' num2str(round(vr_new.agg(rind,azind),1))],...
            ['v_{old} = ' num2str(round(vr_old(rind,azind),1))], 'Location', 'northwest')
        
        set(gcf, 'Units', 'inches', 'Position', [10 10 9 10])
        
        if plot_save_flag
            print([img_dir 'vnew-all-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
        end
        
    end
    
    if strcmp(dtype, '3')
        dt_str = 'Wood Boards';
    elseif strcmp(dtype, '1')
        dt_str = 'Leaves';
    end
    
    if length(dnum) > 6
        dn_str = [dnum(1:end-5) ',' dnum(end-5:end-3) ',' dnum(end-2:end)];
    else
        dn_str = [dnum(1:end-3) ',' dnum(end-2:end)];
    end
    
    
    % Plot of aggregation-filtered PSD + truth/original/corrected velocity
    figure(20)
    clf
    
    semilogy(ones(1,2)*vr_truth(rind,azind), [1e0 1e15], ':k', 'LineWidth', 2.1)
    hold on
    semilogy(ones(1,2)*vr_new.agg(rind,azind), [1e0 1e15], ':r', 'LineWidth', 2.1)
    semilogy(ones(1,2)*vr_old(rind,azind), [1e0 1e15], ':b', 'LineWidth', 2.1)
    semilogy(vvx(:), squeeze(abs(PSD.agg(:,rind,azind))), 'k', 'LineWidth', 1)
    hold off
    xlabel('{\it v} (m s^{-1})', 'FontSize', 16)
    %ylabel('sS_H', 'FontSize', 14)
    xlim([-va, va])
    ylim([1e0 1e15])
    title(['Velocity Correction (' dn_str ' ' dt_str ' at ' elev '\circ)'], 'FontSize', 16)
    legend(['v_{truth} = ' num2str(round(vr_truth(rind,azind),1))],...
        ['v_{new} = ' num2str(round(vr_new.agg(rind,azind),1))],...
        ['v_{old} = ' num2str(round(vr_old(rind,azind),1))], 'Location', 'northwest', 'FontSize', 14)
    
    set(gcf, 'Units', 'inches', 'Position', [4 14 13 4])
    
    if plot_save_flag
        print([img_dir 'vnew-agg-d' dtype 'n' dnum '-el' elev([1 3]) '-1p'], '-dpng')
    end
    
    
    if false
    
    % IMPORTANT!
    % Scatter plots of velocity bias by bulk PHV
    figure(18)
    
    subplot(2,2,1)
    %scatter(reshape(params.phv,[],1), reshape(dv.uncorr.data,[],1), '.')
    scatter(reshape(params.phv,[],1), abs(reshape(dv.uncorr.data,[],1)), '.')
    title('(a) No correction', 'FontSize', 14)
    xlabel('\rho_{HV}', 'FontSize', 14)
    ylabel('v_{bias} (m s^{-1})', 'FontSize', 14)
    xlim([0 1])
    
    subplot(2,2,2)
    %scatter(reshape(params.phv,[],1), reshape(dv.dca.data,[],1), '.')
    scatter(reshape(params.phv,[],1), abs(reshape(dv.dca.data,[],1)), '.')
    title('(b) DCA correction', 'FontSize', 14)
    xlabel('\rho_{HV}', 'FontSize', 14)
    ylabel('v_{bias} (m s^{-1})', 'FontSize', 14)
    
    subplot(2,2,3)
    %scatter(reshape(params.phv,[],1), reshape(dv.var.data,[],1), '.')
    scatter(reshape(params.phv,[],1), abs(reshape(dv.var.data,[],1)), '.')
    title('(c) Variance correction', 'FontSize', 14)
    xlabel('\rho_{HV}', 'FontSize', 14)
    ylabel('v_{bias} (m s^{-1})', 'FontSize', 14)
    
    subplot(2,2,4)
    %scatter(reshape(params.phv,[],1), reshape(dv.agg.data,[],1), '.')
    scatter(reshape(params.phv,[],1), abs(reshape(dv.agg.data,[],1)), '.')
    title('(d) Aggregation correction', 'FontSize', 14)
    xlabel('\rho_{HV}', 'FontSize', 14)
    ylabel('v_{bias} (m s^{-1})', 'FontSize', 14)
    
    axes('Unit', 'Normalized', 'Position', [0.5 0.95 0.01 0.01])
    title(['Debris type ' dtype ', n=' dnum ', tilt=' elev '^o'], 'FontSize', 14);
    axis off
    
    
    % IMPORTANT!
    % Scatter plots of velocity bias by truth velocity
    figure(19)
    
    subplot(2,2,1)
    %scatter(reshape(vr_truth,[],1), reshape(dv_uncorr,[],1), '.')
    scatter(reshape(vr_truth,[],1), abs(reshape(dv.uncorr.data,[],1)), '.')
    title('(a) No correction', 'FontSize', 14)
    xlabel('v_{truth} (m s^{-1})', 'FontSize', 14)
    ylabel('v_{bias} (m s^{-1})', 'FontSize', 14)
    
    subplot(2,2,2)
    %scatter(reshape(vr_truth,[],1), reshape(dv.dca.data,[],1), '.')
    scatter(reshape(vr_truth,[],1), abs(reshape(dv.dca.data,[],1)), '.')
    title('(b) DCA correction', 'FontSize', 14)
    xlabel('v_{truth} (m s^{-1})', 'FontSize', 14)
    ylabel('v_{bias} (m s^{-1})', 'FontSize', 14)
    
    subplot(2,2,3)
    %scatter(reshape(vr_truth,[],1), reshape(dv.var.data,[],1), '.')
    scatter(reshape(vr_truth,[],1), abs(reshape(dv.var.data,[],1)), '.')
    title('(c) Variance correction', 'FontSize', 14)
    xlabel('v_{truth} (m s^{-1})', 'FontSize', 14)
    ylabel('v_{bias} (m s^{-1})', 'FontSize', 14)
    
    subplot(2,2,4)
    %scatter(reshape(vr_truth,[],1), reshape(dv.agg.data,[],1), '.')
    scatter(reshape(vr_truth,[],1), abs(reshape(dv.agg.data,[],1)), '.')
    title('(d) Aggregation correction', 'FontSize', 14)
    xlabel('v_{truth} (m s^{-1})', 'FontSize', 14)
    ylabel('v_{bias} (m s^{-1})', 'FontSize', 14)
    
    axes('Unit', 'Normalized', 'Position', [0.5 0.95 0.01 0.01])
    title(['Debris type ' dtype ', n=' dnum ', tilt=' elev '^o'], 'FontSize', 14);
    axis off
    
    end
    
end

%% Plots for KOUN

if x == 2
    
    
    
    figure(1)
    clf
    
    ha = subplot(2,3,1);
    pcolor(xx, yy, vr_old(:,:,1))
    caxis([-1 1] * va)
    colormap(ha, blib('rgmap2'))
    shading flat
    c = colorbar;
    c.Label.String = 'm s^{-1}';
    c.Label.FontSize = 13;
    c.Label.VerticalAlignment = 'middle';
    title('(a) Uncorrected v_r', 'FontSize', 14)
    xlabel('x (km)', 'FontSize', 14)
    ylabel('y (km)', 'FontSize', 14)
    set(gca, 'DataAspect', [1 1 1])
    ha.Position = [0.2704 0.55 0.1738 0.3412];
    
    ha(2) = subplot(2,3,3);
    pcolor(xx, yy, vr_unfolded(:,:,1))
    caxis([-1 1] * 50)
    colormap(ha(2), blib('rgmap2'))
    shading flat
    c = colorbar;
    c.Label.String = 'm s^{-1}';
    c.Label.FontSize = 13;
    c.Label.VerticalAlignment = 'middle';
    title('(b) Uncorrected v_r (dealiased)', 'FontSize', 14)
    xlabel('x (km)', 'FontSize', 14)
    ylabel('y (km)', 'FontSize', 14)
    set(gca, 'DataAspect', [1 1 1])
    ha(2).Position = [0.5512 0.55 0.1738 0.3412];
    
    ha(3) = subplot(2,3,4);
    pcolor(xx, yy, vr_new.dca(:,:,1))
    caxis([-1 1] * 50)
    colormap(ha(3), blib('rgmap2'))
    shading flat
    c = colorbar;
    c.Label.String = 'm s^{-1}';
    c.Label.FontSize = 13;
    c.Label.VerticalAlignment = 'middle';
    title('(c) DCA-corrected v_r', 'FontSize', 14)
    xlabel('x (km)', 'FontSize', 14)
    ylabel('y (km)', 'FontSize', 14)
    set(gca, 'DataAspect', [1 1 1])
    ha(3).Position = [0.1300 0.1100 0.1738 0.3412];
    
    ha(4) = subplot(2,3,5);
    pcolor(xx, yy, vr_new.var(:,:,1))
    caxis([-1 1] * 50)
    colormap(ha(4), blib('rgmap2'))
    shading flat
    c = colorbar;
    c.Label.String = 'm s^{-1}';
    c.Label.FontSize = 13;
    c.Label.VerticalAlignment = 'middle';
    title('(d) Variance-corrected v_r', 'FontSize', 14)
    xlabel('x (km)', 'FontSize', 14)
    ylabel('y (km)', 'FontSize', 14)
    set(gca, 'DataAspect', [1 1 1])
    ha(4).Position = [0.4108 0.1100 0.1738 0.3412];
    
    ha(5) = subplot(2,3,6);
    pcolor(xx, yy, vr_new.agg(:,:,1))
    caxis([-1 1] * 50)
    colormap(ha(5), blib('rgmap2'))
    shading flat
    c = colorbar;
    c.Label.String = 'm s^{-1}';
    c.Label.FontSize = 13;
    c.Label.VerticalAlignment = 'middle';
    title('(e) Aggregation-corrected v_r', 'FontSize', 14)
    xlabel('x (km)', 'FontSize', 14)
    ylabel('y (km)', 'FontSize', 14)
    set(gca, 'DataAspect', [1 1 1])
    ha(5).Position = [0.6916 0.1100 0.1738 0.3412];
    
    %set(gcf, 'Units', 'inches', 'Position', [5 5 12 8])
    axes('Unit', 'Normalized', 'Position', [0.5 0.95 0.01 0.01])
    title('KOUN 0.5\circ (05-20-2013 20:14 UTC)', 'FontSize', 14);
    axis off
    
    if plot_save_flag
        %print([img_dir 'vcorr-all-d' dtype 'n' dnum '-el' elev([1 3])], '-dpng')
    end
    
    
    figure(2)
    clf
    
    ha = subplot(1,3,1);
    pcolor(xx, yy, vr_new.dca(:,:,1) - vr_unfolded(:,:,1))
    caxis([-1 1] * va)
    colormap(ha, blib('rbmap'))
    shading flat
    c = colorbar;
    c.Label.String = 'm s^{-1}';
    c.Label.FontSize = 13;
    c.Label.VerticalAlignment = 'middle';
    title('(a) DCA correction magnitude', 'FontSize', 14)
    xlabel('x (m)', 'FontSize', 14)
    ylabel('y (m)', 'FontSize', 14)
    set(gca, 'DataAspect', [1 1 1])
    
    ha(2) = subplot(1,3,2);
    pcolor(xx, yy, vr_new.var(:,:,1) - vr_unfolded(:,:,1))
    caxis([-1 1] * va)
    colormap(ha(2), blib('rbmap'))
    shading flat
    c = colorbar;
    c.Label.String = 'm s^{-1}';
    c.Label.FontSize = 13;
    c.Label.VerticalAlignment = 'middle';
    title('(b) Variance correction magnitude', 'FontSize', 14)
    xlabel('x (m)', 'FontSize', 14)
    ylabel('y (m)', 'FontSize', 14)
    set(gca, 'DataAspect', [1 1 1])
    
    ha(3) = subplot(1,3,3);
    pcolor(xx, yy, vr_new.agg(:,:,1) - vr_unfolded(:,:,1))
    caxis([-1 1] * va)
    colormap(ha(3), blib('rbmap'))
    shading flat
    c = colorbar;
    c.Label.String = 'm s^{-1}';
    c.Label.FontSize = 13;
    c.Label.VerticalAlignment = 'middle';
    title('(c) Aggregation correction magnitude', 'FontSize', 14)
    xlabel('x (m)', 'FontSize', 14)
    ylabel('y (m)', 'FontSize', 14)
    set(gca, 'DataAspect', [1 1 1])
    


    figure(3)

    ha = subplot(1,3,1);
    pcolor(xx, yy, vr_old(:,:,1))
    caxis([-1 1] * va)
    colormap(ha, blib('rgmap2'))
    shading flat
    c = colorbar;
    c.Label.String = 'm s^{-1}';
    c.Label.FontSize = 13;
    c.Label.VerticalAlignment = 'middle';
    title('(a) Uncorrected v_r', 'FontSize', 14)
    xlabel('x (km)', 'FontSize', 14)
    ylabel('y (km)', 'FontSize', 14)
    set(gca, 'DataAspect', [1 1 1])

    ha(2) = subplot(1,3,2);
    pcolor(xx, yy, vr_new.agg(:,:,1))
    caxis([-1 1] * va)
    colormap(ha(2), blib('rgmap2'))
    shading flat
    c = colorbar;
    c.Label.String = 'm s^{-1}';
    c.Label.FontSize = 13;
    c.Label.VerticalAlignment = 'middle';
    title('(b) Corrected v_r', 'FontSize', 14)
    xlabel('x (km)', 'FontSize', 14)
    ylabel('y (km)', 'FontSize', 14)
    set(gca, 'DataAspect', [1 1 1])
    
    ha(3) = subplot(1,3,3);
    pcolor(xx, yy, vr_new.agg(:,:,1) - vr_old(:,:,1))
    caxis([-1 1] * va)
    colormap(ha(3), blib('rbmap'))
    shading flat
    c = colorbar;
    c.Label.String = 'm s^{-1}';
    c.Label.FontSize = 13;
    c.Label.VerticalAlignment = 'middle';
    title('(c) Correction magnitude', 'FontSize', 14)
    xlabel('x (m)', 'FontSize', 14)
    ylabel('y (m)', 'FontSize', 14)
    set(gca, 'DataAspect', [1 1 1])
    
    if false
    ha(4) = subplot(2,3,4);
    pcolor(xx, yy, vr_unfolded(:,:,1))
    caxis([-1 1] * va)
    colormap(ha(4), blib('rgmap2'))
    shading flat
    c = colorbar;
    c.Label.String = 'm s^{-1}';
    c.Label.FontSize = 13;
    c.Label.VerticalAlignment = 'middle';
    title('(d) Uncorrected, unfolded v_r', 'FontSize', 14)
    xlabel('x (km)', 'FontSize', 14)
    ylabel('y (km)', 'FontSize', 14)
    set(gca, 'DataAspect', [1 1 1])
    
    ha(5) = subplot(2,3,5);
    pcolor(xx, yy, vr_new.unfolded.agg(:,:,1))
    caxis([-1 1] * va)
    colormap(ha(5), blib('rgmap2'))
    shading flat
    c = colorbar;
    c.Label.String = 'm s^{-1}';
    c.Label.FontSize = 13;
    c.Label.VerticalAlignment = 'middle';
    title('(e) Corrected, unfolded v_r', 'FontSize', 14)
    xlabel('x (km)', 'FontSize', 14)
    ylabel('y (km)', 'FontSize', 14)
    set(gca, 'DataAspect', [1 1 1])
    
    ha(6) = subplot(2,3,6);
    pcolor(xx, yy, vr_new.unfolded.agg(:,:,1) - vr_unfolded(:,:,1))
    caxis([-1 1] * va)
    colormap(ha(6), blib('rbmap'))
    shading flat
    c = colorbar;
    c.Label.String = 'm s^{-1}';
    c.Label.FontSize = 13;
    c.Label.VerticalAlignment = 'middle';
    title('(f) Correction magnitude', 'FontSize', 14)
    xlabel('x (m)', 'FontSize', 14)
    ylabel('y (m)', 'FontSize', 14)
    set(gca, 'DataAspect', [1 1 1])
    end

    
    % aliased inds: 49,56 (8,14) IN; 48,61 (7,19) OUT
    % nonaliased inds: 48,54 (7,12) IN; 48,62 (7,20) or 48,60 (7,18) OUT
    ri = 7;
    azi = 20;
    
    figure(4)
    clf
    
    subplot(2,1,1)
    yyaxis left
    plot(vvx, squeeze(szdr(:,ri,azi)), '-k', 'LineWidth', 1)
    hold on
    plot(ones(1,2)*vr_old(ri,azi), [-15 15], ':b', 'LineWidth', 2)
    plot(ones(1,2)*vr_new.agg(ri,azi), [-15 15], ':r', 'LineWidth', 2)
    hold off
    xlabel('v (m s^{-1})')
    ylabel('sZ_{DR}')
    xlim([-va va])
    ylim([-15 15])
    yyaxis right
    semilogy(vvx, squeeze(PSD.old(:,ri,azi)), ':k', 'LineWidth', 0.7)
    ylabel('sS_H')
    title('(a) sZ_{DR}, sS_H', 'FontSize', 14)
    
    subplot(2,1,2)
    yyaxis left
    plot(vvx, squeeze(sphv(:,ri,azi)), '-k', 'LineWidth', 1)
    hold on
    plot(ones(1,2)*vr_old(ri,azi), [0 1], ':b', 'LineWidth', 2)
    plot(ones(1,2)*vr_new.agg(ri,azi), [0 1], ':r', 'LineWidth', 2)
    hold off
    xlabel('v (m s^{-1})')
    ylabel('s\rho_{HV}')
    xlim([-va va])
    yyaxis right
    semilogy(vvx, squeeze(PSD.old(:,ri,azi)), ':k', 'LineWidth', 0.7)
    ylabel('sS_H')
    title('(b) s\rho_{HV}, sS_H', 'FontSize', 14)

    
    
    
end




