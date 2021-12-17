function [obj_class, agg, memf] = hca_class(szdr, sphv, svar)
%
% INPUTS
%   szdr:   DPSD for differential reflectivity (3D array, size M*r*az)
%   sphv:   DPSD for correlation coefficient (3D array, size M*r*az)
%   svar:   n-point moving variances (struct)
%     svar.zdr = sZDR moving variance (3D array, size M*r*az)
%     svar.phv = sPHV moving variance (3D array, size M*r*az)
%
% OUTPUTS
%   obj_class:  HCA classification vector (3D array, size M*r*az)
%       1 for debris, 0 for rain
%
%   agg:        Aggregation parameters (struct of 3D arrays, size M*r*az)
%       agg.debris:      Debris aggregation parameter
%       agg.rain:        Rain aggregation parameter
%
%   memf:       Membership values (struct of 3D arrays, size M*r*az)
%       memf.zdr_debr:   Debris-szdr membership
%       memf.zdr_rain:   Rain-szdr membership
%       memf.phv_debr:   Debris-sphv membership
%       memf.phv_rain:   Rain-sphv membership
%       memf.pvar_debr:  Debris-pvar membership
%       memf.pvar_rain:  Rain-pvar membership
%       memf.zvar_debr:  Debris-zvar membership
%       memf.zvar_rain:  Rain-zvar membership

%%%%% Weights for calculating aggregation parameters

a = 0.2; % weight for sZDR
b = 0.3; % weight for sPHV
c = 0.5; % weight for sVAR

% a = 1; % weight for sZDR
% b = 1; % weight for sPHV
% c = 1; % weight for sPHV variance
% d = 1; % weight for sZDR variance

%%%%%

M = size(szdr,1);
pvar = svar.phv;
zvar = svar.zdr;

obj_class = zeros(M, size(szdr,2), size(szdr,3));

% thres_phv_rain = [0.9, 0.95]; % sPHV thresholds for rain
% thres_phv_debr = [0.85, 1]; % sPHV thresholds for debris
% thres_zdr_rain = [-2, 0, 3, 5]; % sZDR thresholds for rain
% thres_zdr_debr = [-20, -6, 2, 10]; % sZDR thresholds for debris
% thres_pvar_rain = [0.03, 0.07]; % sPHV variance thresholds for rain
% thres_pvar_debr = [0, 0.08, 0.15, 0.25]; % sPHV variance thresholds for debris
% thres_zvar_rain = [0.5, 1.5]; % sZDR variance thresholds for rain
% thres_zvar_debr = [0, 10]; % sZDR variance thresholds for debris


%%% 0.5, 20, 50, 80, 99.5 percentiles -- [x1, x2, (mu), x3, x4] %%%
% thres_zdr_rain = [-1.5, 1, (1.5), 2, 3.7];
% thres_zdr_debr = [-18, -5, (-0.8), 2, 15];
% thres_phv_rain = [0.8, 0.98, (0.9963), 0.999, 1];
% thres_phv_debr = [0, 0.6, (0.9117), 0.98, 1];
% thres_pvar_rain = [0, 5e-6, (3.17e-5), 2e-4, 0.025];
% thres_pvar_debr = [3e-6, 0.0014, (0.03), 0.1032, 0.1924];
% thres_zvar_rain = [0.015, 0.12, (0.3), 0.65, 4.8];
% thres_zvar_debr = [0.1, 2.2, (9.3), 23, 92];

thres_zdr_rain = [-1.5, 1, 2, 4];
thres_zdr_debr = [-18, -5, 2, 15];
thres_phv_rain = [0.8, 0.95]; % [x1, x2] - 5th percentile 0.9551
thres_phv_debr = [0.6, 0.91]; % [x3, x4] - 50th percentile 0.9117
thres_pvar_rain = [0.002, 0.025]; % [x3, x4] - 95th percentile 0.0018
thres_pvar_debr = [0, 0.03, 0.11, 0.25]; % 50th percentile 0.0308
thres_zvar_rain = [0, 0.12, 0.65, 5];
thres_zvar_debr = [0.1, 2.2, 23, 100];



%% Debris membership


mem_zdr_debr = zeros(M, size(szdr,2), size(szdr,3));
mem_phv_debr = ones(M, size(szdr,2), size(szdr,3));
mem_pvar_debr = ones(M, size(szdr,2), size(szdr,3));
mem_zvar_debr = ones(M, size(szdr,2), size(szdr,3));

% sPHV
mem_phv_debr(sphv<thres_phv_debr(1)) = 1;
mem_phv_debr(sphv<thres_phv_debr(2) & sphv>=thres_phv_debr(1)) = 1 - ((sphv(sphv<thres_phv_debr(2) & sphv>=thres_phv_debr(1)) - thres_phv_debr(1)) / (thres_phv_debr(2) - thres_phv_debr(1)));
mem_phv_debr(sphv>=thres_phv_debr(2)) = 0;

% sZDR
mem_zdr_debr(szdr<thres_zdr_debr(1)) = 0;
mem_zdr_debr(szdr<thres_zdr_debr(2) & szdr>=thres_zdr_debr(1)) = (szdr(szdr<thres_zdr_debr(2) & szdr>=thres_zdr_debr(1)) - thres_zdr_debr(1)) / (thres_zdr_debr(2) - thres_zdr_debr(1));
mem_zdr_debr(szdr<thres_zdr_debr(3) & szdr>=thres_zdr_debr(2)) = 1;
mem_zdr_debr(szdr<thres_zdr_debr(4) & szdr>=thres_zdr_debr(3)) = 1 - ((szdr(szdr<thres_zdr_debr(4) & szdr>=thres_zdr_debr(3)) - thres_zdr_debr(3)) / (thres_zdr_debr(4) - thres_zdr_debr(3)));
mem_zdr_debr(szdr>=thres_zdr_debr(4)) = 0;

% sPHV variance
mem_pvar_debr(pvar<thres_pvar_debr(1)) = 0;
mem_pvar_debr(pvar<thres_pvar_debr(2) & pvar>=thres_pvar_debr(1)) = (pvar(pvar<thres_pvar_debr(2) & pvar>=thres_pvar_debr(1)) - thres_pvar_debr(1)) / (thres_pvar_debr(2) - thres_pvar_debr(1));
mem_pvar_debr(pvar<thres_pvar_debr(3) & pvar>=thres_pvar_debr(2)) = 1;
mem_pvar_debr(pvar<thres_pvar_debr(4) & pvar>=thres_pvar_debr(3)) = 1 - ((pvar(pvar<thres_pvar_debr(4) & pvar>=thres_pvar_debr(3)) - thres_pvar_debr(3)) / (thres_pvar_debr(4) - thres_pvar_debr(3)));
mem_pvar_debr(pvar>=thres_pvar_debr(4)) = 0;

% sZDR variance
mem_zvar_debr(zvar<thres_zvar_debr(1)) = 0;
mem_zvar_debr(zvar<thres_zvar_debr(2) & zvar>=thres_zvar_debr(1)) = (zvar(zvar<thres_zvar_debr(2) & zvar>=thres_zvar_debr(1)) - thres_zvar_debr(1)) / (thres_zvar_debr(2) - thres_zvar_debr(1));
mem_zvar_debr(zvar<thres_zvar_debr(3) & zvar>=thres_zvar_debr(2)) = 1;
mem_zvar_debr(zvar<thres_zvar_debr(4) & zvar>=thres_zvar_debr(3)) = 1 - ((zvar(zvar<thres_zvar_debr(4) & zvar>=thres_zvar_debr(3)) - thres_zvar_debr(3)) / (thres_zvar_debr(4) - thres_zvar_debr(3)));
mem_zvar_debr(zvar>=thres_zvar_debr(4)) = 0;

%% Rain membership

mem_zdr_rain = zeros(M, size(szdr,2), size(szdr,3));
mem_phv_rain = ones(M, size(szdr,2), size(szdr,3));
mem_pvar_rain = ones(M, size(szdr,2), size(szdr,3));
mem_zvar_rain = ones(M, size(szdr,2), size(szdr,3));

% sPHV
mem_phv_rain(sphv<thres_phv_rain(1)) = 0;
mem_phv_rain(sphv<thres_phv_rain(2) & sphv>=thres_phv_rain(1)) = (sphv(sphv<thres_phv_rain(2) & sphv>=thres_phv_rain(1)) - thres_phv_rain(1)) / (thres_phv_rain(2) - thres_phv_rain(1));
mem_phv_rain(sphv>=thres_phv_rain(2)) = 1;

% sZDR
mem_zdr_rain(szdr<thres_zdr_rain(1)) = 0;
mem_zdr_rain(szdr<thres_zdr_rain(2) & szdr>=thres_zdr_rain(1)) = (szdr(szdr<thres_zdr_rain(2) & szdr>=thres_zdr_rain(1)) - thres_zdr_rain(1)) / (thres_zdr_rain(2) - thres_zdr_rain(1));
mem_zdr_rain(szdr<thres_zdr_rain(3) & szdr>=thres_zdr_rain(2)) = 1;
mem_zdr_rain(szdr<thres_zdr_rain(4) & szdr>=thres_zdr_rain(3)) = 1 - ((szdr(szdr<thres_zdr_rain(4) & szdr>=thres_zdr_rain(3)) - thres_zdr_rain(3)) / (thres_zdr_rain(4) - thres_zdr_rain(3)));
mem_zdr_rain(szdr>=thres_zdr_rain(4)) = 0;

% sPHV variance
mem_pvar_rain(pvar<thres_pvar_rain(1)) = 1;
mem_pvar_rain(pvar<thres_pvar_rain(2) & pvar>=thres_pvar_rain(1)) = 1 - ((pvar(pvar<thres_pvar_rain(2) & pvar>=thres_pvar_rain(1)) - thres_pvar_rain(1)) / (thres_pvar_rain(2) - thres_pvar_rain(1)));
mem_pvar_rain(pvar>=thres_pvar_rain(2)) = 0;

% sZDR variance
mem_zvar_rain(zvar<thres_zvar_rain(1)) = 0;
mem_zvar_rain(zvar<thres_zvar_rain(2) & zvar>=thres_zvar_rain(1)) = (zvar(zvar<thres_zvar_rain(2) & zvar>=thres_zvar_rain(1)) - thres_zvar_rain(1)) / (thres_zvar_rain(2) - thres_zvar_rain(1));
mem_zvar_rain(zvar<thres_zvar_rain(3) & zvar>=thres_zvar_rain(2)) = 1;
mem_zvar_rain(zvar<thres_zvar_rain(4) & zvar>=thres_zvar_rain(3)) = 1 - ((zvar(zvar<thres_zvar_rain(4) & zvar>=thres_zvar_rain(3)) - thres_zvar_rain(3)) / (thres_zvar_rain(4) - thres_zvar_rain(3)));
mem_zvar_rain(zvar>=thres_zvar_rain(4)) = 0;


%% Calculate membership

agg_debr = a*mem_zdr_debr + b*mem_phv_debr + c*mem_pvar_debr;
agg_rain = a*mem_zdr_rain + b*mem_phv_rain + c*mem_pvar_rain;
% agg_debr = a*mem_zdr_debr + b*mem_phv_debr + c*mem_pvar_debr + d*mem_zvar_debr;
% agg_rain = a*mem_zdr_rain + b*mem_phv_rain + c*mem_pvar_rain + d*mem_zvar_rain;

agg.debr = agg_debr;
agg.rain = agg_rain;

obj_class(agg_debr > agg_rain) = 1;

memf = struct('debr', [], 'rain', []);
memf.debr = struct('zdr', mem_zdr_debr, 'phv', mem_phv_debr,...
    'pvar', mem_pvar_debr, 'zvar', mem_zvar_debr);
memf.rain = struct('zdr', mem_zdr_rain, 'phv', mem_phv_rain,...
    'pvar', mem_pvar_rain, 'zvar', mem_zvar_rain);

% memf = struct('zdr_debr', mem_zdr_debr, 'zdr_rain', mem_zdr_rain, 'phv_debr', ...
%     mem_phv_debr, 'phv_rain', mem_phv_rain, 'pvar_debr', mem_pvar_debr, 'pvar_rain', ...
%     mem_pvar_rain, 'zvar_debr', mem_zvar_debr, 'zvar_rain', mem_zvar_rain);

end