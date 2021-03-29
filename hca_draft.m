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

vel

% szdr = mult.szdr(:,14,20,1,3);
% sphv = mult.sphv(:,14,20,1,3);
% svar = sphv_var_mult(:,14,20,1,3);
% vel = vvx;

M = length(vel);

hca_class = zeros(M,1);

tpr = [0.9, 0.95]; % threshold phv rain
tpd = [0, 1]; % threshold phv debris
tzr = [-2, 0, 3, 5]; % threshold zdr rain
tzd = [-20, -10, 5, 15]; % threshold zdr debris
tvr = [0.02, 0.05]; % threshold var rain
tvd = [0.05, 0.1, 0.15, 0.25]; % threshold var debris
 
% test for debris
mem_zdr_debr = zeros(M,1);
mem_phv_debr = ones(M,1);
mem_var_debr = ones(M,1);

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


mem_var_debr(svar < tvd(1)) = 0.5;
mem_var_debr(svar < tvd(2) & svar >= tvd(1)) = ...
    (svar(svar < tvd(2) & svar >= tvd(1)) - tvd(1)) / (tvd(2) - tvd(1));
mem_var_debr(svar < tvd(3) & svar >= tvd(2)) = 1;
mem_var_debr(svar < tvd(4) & svar >= tvd(3)) = ...
    1 - ((svar(svar < tvd(4) & svar >= tvd(3)) - tvd(3)) / (tvd(4) - tvd(3)));
mem_var_debr(svar >= tvd(4)) = 0.3;


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


agg_debr = mem_zdr_debr + mem_phv_debr + mem_var_debr;
agg_rain = mem_zdr_rain + mem_phv_rain + mem_var_rain;

hca_class(agg_debr > agg_rain) = 1;


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
plot(vel, hca_class)
title('HCA class (debr=1 rain=0)')

