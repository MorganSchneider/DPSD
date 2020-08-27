function out = dmap(n)
% Differential Reflectivity Map
% WSR-88D Operational
%
% For caxis -> [-5 8]
%
%

colors = nan(16,3);
colors(1,:) = [hex2dec('00') hex2dec('00') hex2dec('00')]/255;
colors(2,:) = [hex2dec('40') hex2dec('40') hex2dec('40')]/255;
colors(3,:) = [hex2dec('9C') hex2dec('9C') hex2dec('9C')]/255;
colors(4,:) = [hex2dec('C9') hex2dec('C9') hex2dec('C9')]/255;
colors(5,:) = [hex2dec('8C') hex2dec('78') hex2dec('B4')]/255;
colors(6,:) = [hex2dec('00') hex2dec('00') hex2dec('98')]/255;
colors(7,:) = [hex2dec('23') hex2dec('98') hex2dec('D3')]/255;
colors(8,:) = [hex2dec('44') hex2dec('FF') hex2dec('D2')]/255;
colors(9,:) = [hex2dec('57') hex2dec('DB') hex2dec('56')]/255;
colors(10,:) = [hex2dec('FF') hex2dec('FF') hex2dec('60')]/255;
colors(11,:) = [hex2dec('FF') hex2dec('90') hex2dec('45')]/255;
colors(12,:) = [hex2dec('DA') hex2dec('00') hex2dec('00')]/255;
colors(13,:) = [hex2dec('AE') hex2dec('00') hex2dec('00')]/255;
colors(14,:) = [hex2dec('F7') hex2dec('82') hex2dec('BE')]/255;
colors(15,:) = [hex2dec('FF') hex2dec('FF') hex2dec('FF')]/255;
colors(16,:) = [hex2dec('77') hex2dec('00') hex2dec('7D')]/255;

cmap = [...
   -5, colors(1,:); 
   -4, colors(2,:);
   -2, colors(3,:);
 -0.5, colors(4,:);
    0, colors(5,:);
 0.25, colors(6,:);
  0.5, colors(7,:);
    1, colors(8,:);
  1.5, colors(9,:);
    2, colors(10,:);
  2.5, colors(11,:);
    3, colors(12,:);
    4, colors(13,:);
    5, colors(14,:);
    6, colors(15,:);
    8, colors(15,:);
%   nan, colors(16,:)
    ];

out = interp1(cmap(:,1), cmap(:,2:4), linspace(cmap(1,1), cmap(end,1), n));