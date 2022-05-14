el = 0.5;
obsdate = '20130520';
fname = ['KOUN_' obsdate '_' num2str(el,'%.1f') '.mat'];

if ~isfile(fname)
    fEL_desired = el;
    read_data
end

gates = size(iqh,1);
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
[az_mat, r_mat] = meshgrid(az,r);
xx = r_mat .* sind(az_mat) * cosd(el);
yy = r_mat .* cosd(az_mat) * cosd(el);
zz = r_mat * sind(el);
vr = -va/pi*angle(R1);

zdr = 10*log10(squeeze(sh./sv));
rhohv = abs(sx2) ./ sqrt(sh2.*sv2);

sh = 10*log10(squeeze(sh));
sv = 10*log10(squeeze(sv));
zcor = header.fSYCALDB;
rcor = 10*log10((r(:)*1000).^2);
rcor = repmat(rcor, [1 numel(az)]);
zh = sh - zcor + rcor;

% Dealiasing
% vr_unfolded = vr;
% vr_unfolded(47:49,55:56) = vr(47:49,55:56) - 2*va;
% vr_unfolded(46:49,57) = vr(46:49,57) - 2*va;
% vr_unfolded(47,58) = vr(47,58) - 2*va;
% vr_unfolded(47,59:63) = vr(47,59:63) + 2*va;
% vr_unfolded(48,61) = vr(48,61) + 2*va;

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

if false
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
end

%% 
tmp = struct('P', P, 'zh', zh, 'vr', vr, 'zdr', zdr, 'rhohv', rhohv);
tds = tmp;

% storm thresholds
P_thres = -20;
v_thres = 1;
r_thres = 40;

inds = zeros(size(P));
inds((10*log10(P) > P_thres) & (abs(vr) > v_thres) & (r_mat < r_thres)) = 1;
tmp.P(~inds) = nan;
tmp.zh(~inds) = nan;
tmp.vr(~inds) = nan;
tmp.zdr(~inds) = nan;
tmp.rhohv(~inds) = nan;


% TDS thresholds
zh_thres = 30;
zdr_thres = 1;
rho_thres = 0.82;

itds = inds;
itds((tmp.zh < zh_thres) | (tmp.zdr > zdr_thres) | (tmp.rhohv > rho_thres)) = 0;
tds.P(~itds) = nan;
tds.zh(~itds) = nan;
tds.vr(~itds) = nan;
tds.zdr(~itds) = nan;
tds.rhohv(~itds) = nan;


rc = nanmedian(r_mat(itds==1));
azc = nanmedian(az_mat(itds==1));
xc = rc * sind(azc);
yc = rc * cosd(azc);

rmin = rc - 3;
rmax = rc + 3;
amin = azc - 10;
amax = azc + 10;

% need to add iqh, iqv to this structure
[ri,ai] = find(r_mat>=rmin & r_mat<=rmax & az_mat>=amin & az_mat<=amax);
tor = struct('xx', xx(ri,ai), 'yy', yy(ri,ai), 'zh', zh(ri,ai),...
    'vr', vr(ri,ai), 'zdr', zdr(ri,ai), 'rhohv', rhohv(ri,ai));



%%


figure()

subplot(2,4,1)
pcolor(xx, yy, tmp.zh)
caxis([0 70])
colormap(gca, blib('zmap'))
colorbar
shading flat
xlim([-20 30])
ylim([0 40])
%xlim([-15 0])
%ylim([5 15])
%set(gca, 'DataAspect', [1 1 1])

subplot(2,4,2)
pcolor(xx, yy, tmp.vr)
caxis([-va va])
colormap(gca, blib('rgmap2'))
colorbar
shading flat
xlim([-20 30])
ylim([0 40])
%set(gca, 'DataAspect', [1 1 1])

subplot(2,4,3)
pcolor(xx, yy, tmp.zdr)
caxis([-5 5])
colormap(gca, blib('nwsdmap'))
colorbar
shading flat
xlim([-20 30])
ylim([0 40])
%set(gca, 'DataAspect', [1 1 1])

subplot(2,4,4)
pcolor(xx, yy, tmp.rhohv)
caxis([0 1])
colormap(gca, blib('nwsrmap'))
colorbar
shading flat
xlim([-20 30])
ylim([0 40])
%set(gca, 'DataAspect', [1 1 1])

subplot(2,4,5)
pcolor(xx, yy, tds.zh)
caxis([0 70])
colormap(gca, blib('zmap'))
colorbar
shading flat
hold on
scatter(xc, yc, 80, '.k')
hold off
xlim([-15 0])
ylim([5 15])

subplot(2,4,6)
pcolor(xx, yy, tds.vr)
caxis([-va va])
colormap(gca, blib('rgmap2'))
colorbar
shading flat
hold on
scatter(xc, yc, 80, '.k')
hold off
xlim([-15 0])
ylim([5 15])

subplot(2,4,7)
pcolor(xx, yy, tds.zdr)
caxis([-5 5])
colormap(gca, blib('nwsdmap'))
colorbar
shading flat
hold on
scatter(xc, yc, 80, '.k')
hold off
xlim([-15 0])
ylim([5 15])

subplot(2,4,8)
pcolor(xx, yy, tds.rhohv)
caxis([0 1])
colormap(gca, blib('nwsrmap'))
colorbar
shading flat
hold on
scatter(xc, yc, 80, '.k')
hold off
xlim([-15 0])
ylim([5 15])



figure()

subplot(2,2,1)
pcolor(tor.xx, tor.yy, tor.zh)
caxis([0 70])
colormap(gca, blib('zmap'))
colorbar
shading flat
hold on
scatter(xc, yc, 400, '.k')
hold off

subplot(2,2,2)
pcolor(tor.xx, tor.yy, tor.vr)
caxis([-va va])
colormap(gca, blib('rgmap2'))
colorbar
shading flat
hold on
scatter(xc, yc, 400, '.k')
hold off

subplot(2,2,3)
pcolor(tor.xx, tor.yy, tor.zdr)
caxis([-5 5])
colormap(gca, blib('nwsdmap'))
colorbar
shading flat
hold on
scatter(xc, yc, 400, '.k')
hold off

subplot(2,2,4)
pcolor(tor.xx, tor.yy, tor.rhohv)
caxis([0 1])
colormap(gca, blib('nwsrmap'))
colorbar
shading flat
hold on
scatter(xc, yc, 400, '.k')
hold off



% save(['Documents/code/obsdata/KOUN_data_', num2str(el,'%.1f'), '.mat'],...
%    'iqh_tds', 'iqv_tds', 'vr', 'vr_unfolded', 'xx', 'yy', 'zz', 'zh', 'zdr',...
%    'rhohv', 'P', 'va', 'tds', 'tvs', 'iqh_full', 'iqv_full', 'Ts', 'lambda')

% save(['Documents/code/obsdata/KOUN_' obsdate '_' num2str(el,'%.1f') '.mat'],...
%     'iqh', 'iqv', 'xx', 'yy', 'zz', 'vr', 'zh', 'zdr', 'rhohv', 'va')
