function out = rmap(n)
% Differential Reflectivity Map
% WSR-88D Operational
%
% For caxis -> [0.1 1.05]
%
%

colors = nan(16,3);
colors(1,:) = [hex2dec('00') hex2dec('00') hex2dec('00')]/255;
colors(2,:) = [hex2dec('95') hex2dec('94') hex2dec('9C')]/255;
colors(3,:) = [hex2dec('16') hex2dec('14') hex2dec('8C')]/255;
colors(4,:) = [hex2dec('09') hex2dec('02') hex2dec('D9')]/255;
colors(5,:) = [hex2dec('89') hex2dec('87') hex2dec('D6')]/255;
colors(6,:) = [hex2dec('5C') hex2dec('FF') hex2dec('59')]/255;
colors(7,:) = [hex2dec('8B') hex2dec('CF') hex2dec('02')]/255;
colors(8,:) = [hex2dec('FF') hex2dec('FB') hex2dec('00')]/255;
colors(9,:) = [hex2dec('FF') hex2dec('C4') hex2dec('00')]/255;
colors(10,:) = [hex2dec('FF') hex2dec('89') hex2dec('03')]/255;
colors(11,:) = [hex2dec('FF') hex2dec('2B') hex2dec('00')]/255;
colors(12,:) = [hex2dec('E3') hex2dec('00') hex2dec('00')]/255;
colors(13,:) = [hex2dec('A1') hex2dec('00') hex2dec('00')]/255;
colors(14,:) = [hex2dec('97') hex2dec('05') hex2dec('56')]/255;
colors(15,:) = [hex2dec('FA') hex2dec('AC') hex2dec('D1')]/255;
colors(16,:) = [hex2dec('77') hex2dec('00') hex2dec('7D')]/255;

% Corresponding to 0.2 < rho <= 1.05

cmap = [...
 0.10, colors(1,:); 
 0.20, colors(2,:);
 0.45, colors(3,:);
 0.65, colors(4,:);
 0.75, colors(5,:);
 0.80, colors(6,:);
 0.85, colors(7,:);
 0.90, colors(8,:);
 0.93, colors(9,:);
 0.95, colors(10,:);
 0.96, colors(11,:);
 0.97, colors(12,:);
 0.98, colors(13,:);
 0.99, colors(14,:);
1.005, colors(14,:);
 1.01, colors(15,:);
 1.05, colors(15,:);
%   nan, colors(16,:)
    ];

out = interp1(cmap(:,1), cmap(:,2:4), linspace(cmap(1,1), cmap(end,1), n));