function [obj_class, agg, varargout] = hca_class(szdr, sphv, svar)
%
% INPUTS
%   szdr:   DPSD for differential reflectivity (column vector)
%   sphv:   DPSD for correlation coefficient (column vector)
%   svar:   n-point moving variance of sphv (column vector)
%
% OUTPUTS
%   obj_class:  HCA classification vector (column vector)
%       obj_class =1 for debris, =0 for rain
%   varargout:  Membership functions (cell array)
%       varargout{1}:   Debris-szdr membership function (column vector)
%       varargout{2}:   Debris-sphv membership function (column vector)
%       varargout{3}:   Debris-svar membership function (column vector)
%       varargout{4}:   Rain-szdr membership function (column vector)
%       varargout{5}:   Rain-sphv membership function (column vector)
%       varargout{6}:   Rain svar membership function (column vector)

%%%%% Weights for calculating aggregation parameter

a = 0.3; % weight for sZDR
b = 0.3; % weight for sPHV
c = 0.4; % weight for sVAR

%%%%%

M = length(szdr);

obj_class = zeros(M,size(szdr,2),size(szdr,3));

tpr = [0.85, 0.9]; % threshold phv rain
tpd = tpr(1); % threshold phv debris
tzr = [-2, 0, 3, 5]; % threshold zdr rain
tzd = [-20, -6, 5, 15]; % threshold zdr debris
tvr = [0.03, 0.07]; % threshold var rain
tvd = [0, 0.05, 0.15, 0.25]; % threshold var debris


% test for debris
mem_zdr_debr = zeros(M,size(szdr,2),size(szdr,3));
mem_phv_debr = ones(M,size(szdr,2),size(szdr,3));
mem_var_debr = ones(M,size(szdr,2),size(szdr,3));

mem_phv_debr(sphv < tpd) = 1;
mem_phv_debr(sphv >= tpd) = 0.5;

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


mem_var_debr(svar < tvd(1)) = 0;
mem_var_debr(svar < tvd(2) & svar >= tvd(1)) = ...
    (svar(svar < tvd(2) & svar >= tvd(1)) - tvd(1)) / (tvd(2) - tvd(1));
mem_var_debr(svar < tvd(3) & svar >= tvd(2)) = 1;
mem_var_debr(svar < tvd(4) & svar >= tvd(3)) = ...
    1 - ((svar(svar < tvd(4) & svar >= tvd(3)) - tvd(3)) / (tvd(4) - tvd(3)));
mem_var_debr(svar >= tvd(4)) = 0;


% test for rain
mem_zdr_rain = zeros(M,size(szdr,2),size(szdr,3));
mem_phv_rain = ones(M,size(szdr,2),size(szdr,3));
mem_var_rain = ones(M,size(szdr,2),size(szdr,3));


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



agg_debr = a*mem_zdr_debr + b*mem_phv_debr + c*mem_var_debr;
agg_rain = a*mem_zdr_rain + b*mem_phv_rain + c*mem_var_rain;

agg.debr = agg_debr;
agg.rain = agg_rain;

obj_class(agg_debr > agg_rain) = 1;

varargout{1} = mem_zdr_debr;
varargout{2} = mem_phv_debr;
varargout{3} = mem_var_debr;
varargout{4} = mem_zdr_rain;
varargout{5} = mem_phv_rain;
varargout{6} = mem_var_rain;

end