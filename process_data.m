if ~exist('header', 'var')
    read_data
end

gates = size(iqh,1);
delta_R = header.fRangeRes;
va = header.fWavelengthCM*1e-2/(4*header.fPRTUSec*1e-6);
M = size(iqh,3);
Ts = header.fPRTUSec * 1e-6; % convert PRT from us to s
lambda = header.fWavelengthCM * 1e-2; % convert wavelength from cm to m

% From https://www.roc.noaa.gov/WSR88D/Engineering/NEXRADTechInfo.aspx:
antenna_gain_db = 45.5;
tx_power_watt = 7e5; %700 kW

P = squeeze(mean(abs(iqh).^2,3));
R1 = squeeze(mean(conj(iqh(:,:,1:M-1)).*iqh(:,:,2:M),3));

sh = real(mean(iqh .* conj(iqh), 3));
sv = real(mean(iqv .* conj(iqv), 3));
mh = repmat(mean(iqh,3), [1 1 M]);
mv = repmat(mean(iqv,3), [1 1 M]);
sh2 = mean((iqh - mh) .* conj(iqh - mh), 3);
sv2 = mean((iqv - mv) .* conj(iqv - mv), 3);
sx2 = mean((iqh - mh) .* conj(iqv - mv), 3);

r = (1:gates)*0.25;
el = fEL_desired;
[az_mat, r_mat] = meshgrid(az,r);
xx = r_mat .* sind(az_mat) * cosd(el);
yy = r_mat .* cosd(az_mat) * cosd(el);
zz = r_mat * sind(el);
vr = -va/pi*angle(R1);
vr_unfolded = vr;

zdr = 10*log10(squeeze(sh./sv));
rhohv = abs(sx2) ./ sqrt(sh2.*sv2);

sh = 10*log10(squeeze(sh));
sv = 10*log10(squeeze(sv));
zcor = header.fSYCALDB;
rcor = 10*log10((r(:)*1000).^2);
rcor = repmat(rcor, [1 numel(az)]);
zh = sh - zcor + rcor;

% Dealiasing
vr_unfolded(47:49,55:56) = vr(47:49,55:56) - 2*va;
vr_unfolded(46:49,57) = vr(46:49,57) - 2*va;
vr_unfolded(47,58) = vr(47,58) - 2*va;
vr_unfolded(47,59:63) = vr(47,59:63) + 2*va;
vr_unfolded(48,61) = vr(48,61) + 2*va;

%% Plots

figure(6)
pcolor(xx,yy,10*log10(P))
ylabel('y (km)')
xlabel('x (km)')
title('Power (dB)')
colorbar
shading flat
set(gca, 'DataAspect', [1 1 1])

figure(7)
%pcolor(xx,yy,10*log10(P))
pcolor(xx,yy,zh)
ylabel('y (km)')
xlabel('x (km)')
%title('Power (dB)')
title('Reflectivity (dB)')
caxis([0 50])
colorbar
shading flat
set(gca, 'DataAspect', [1 1 1])

figure(8)
pcolor(xx,yy,vr)
ylabel('y (km)')
xlabel('x (km)')
title('Velocity (m/s)')
colorbar
shading flat
set(gca, 'DataAspect', [1 1 1])

tds.xlim = [-9.5, -5.5];
tds.ylim = [8, 11];

tds.xinds = 42:55;
tds.yinds = 43:70;

tvs.xinds = 44:52;
tvs.yinds = 48:68;

rho_clims = [0.3, 1];
zdr_clims = [-5, 5];

figure(1)
pcolor(xx(tds.xinds,tds.yinds), yy(tds.xinds,tds.yinds), 10*log10(P(tds.xinds,tds.yinds)))
colorbar
shading flat
xlabel('x (km)')
ylabel('y (km)')
title('Power (dB)')
% xlim([-9.5, -5.5])
% ylim([8, 11])
set(gca, 'DataAspect', [1 1 1])

figure(2)
pcolor(xx(tds.xinds,tds.yinds), yy(tds.xinds,tds.yinds), vr(tds.xinds,tds.yinds))
colormap(blib('rgmap2'))
colorbar
shading flat
xlabel('x (km)')
ylabel('y (km)')
title('Velocity (m/s)')
% xlim([-9.5, -5.5])
% ylim([8, 11])
set(gca, 'DataAspect', [1 1 1])

figure(3)
pcolor(xx(tds.xinds,tds.yinds), yy(tds.xinds,tds.yinds), zdr(tds.xinds,tds.yinds))
colormap(blib('nwsdmap'))
caxis(zdr_clims)
colorbar
shading flat
xlabel('x (km)')
ylabel('y (km)')
title('ZDR (dB)')
% xlim([-9.5, -5.5])
% ylim([8, 11])
set(gca, 'DataAspect', [1 1 1])

figure(4)
pcolor(xx(tds.xinds,tds.yinds), yy(tds.xinds,tds.yinds), rhohv(tds.xinds,tds.yinds))
colormap(blib('nwsrmap'))
caxis(rho_clims)
colorbar
shading flat
xlabel('x (km)')
ylabel('y (km)')
title('rhoHV')
set(gca, 'DataAspect', [1 1 1])


figure(5)
pcolor(xx(tds.xinds, tds.yinds), yy(tds.xinds, tds.yinds), vr_unfolded(tds.xinds, tds.yinds))
colormap(blib('rgmap2'))
colorbar
shading flat
xlabel('x (km)')
ylabel('y (km)')
title('Dealiased Velocity (m/s)')
set(gca, 'DataAspect', [1 1 1])

figure(6)
pcolor(xx(tvs.xinds, tvs.yinds), yy(tvs.xinds, tvs.yinds), vr_unfolded(tvs.xinds, tvs.yinds))
colormap(blib('rgmap2'))
colorbar
shading flat
xlabel('x (km)')
ylabel('y (km)')
title('Dealiased Velocity (m/s)')
set(gca, 'DataAspect', [1 1 1])

iqh_full = iqh;
iqv_full = iqv;
iqh_tds = iqh(tds.xinds, tds.yinds, :);
iqv_tds = iqv(tds.xinds, tds.yinds, :);

save(['~/Documents/code/obsdata/KOUN_data_', num2str(el,'%.1f'), '.mat'],...
    'iqh_tds', 'iqv_tds', 'vr', 'vr_unfolded', 'xx', 'yy', 'zz', 'zh', 'zdr',...
    'rhohv', 'P', 'va', 'tds', 'tvs', 'iqh_full', 'iqv_full', 'Ts', 'lambda')

clear vr vr_unfolded