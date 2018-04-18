% getSourceSink4.m
%
% author: Julia Dohner
% April 17, 2018
%
% provides record of fossil fuel emissions, land use emissions, and
% extratropical land use emissions for the time frame 1765 to present.
%
% note: need to add framework to append on datasets that go to the future.
% contact Fabian to see if his future emission scenarios will include data
% all the way back to 1765 or if just to present.

function [ff,LU,LUhigher] = getSourceSink4(year,ts);

addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/joosModel/co2_forward_data_2016/JDfiles'));

FF_2016 = csvread('GCP_FFann_1751-2016.csv'); % in gigatons/year
LU_2016 = csvread('Houghton_LUann_1850-2016_2.csv'); % in gigatons/year


d = 1/2.31; % 1 ppm CO2 = 2.31 gton CO2

% interpolate to monthly
FF_2016mo(:,1) = (1765:(1/ts):FF_2016(end,1))';
FF_2016mo(:,2) = (interp1(FF_2016(:,1),FF_2016(:,2),FF_2016mo(:,1)));
FF_2016mo(:,2) = FF_2016mo(:,2)*d; %convert to ppm

% set LU = 0 in 1765, interpolated back
LU_2016mo(:,1) = (1765:(1/ts):LU_2016(end,1))';
LU_2016mo(:,2) = interp1(LU_2016(:,1),LU_2016(:,2),LU_2016mo(:,1));
LU_2016mo(:,2) = LU_2016mo(:,2)*d; % convert to ppm

% shorten datasets to match time frame of year vector (currently thru 2016)
FF_start = find(FF_2016mo(:,1) == year(1));
FF_end = find(FF_2016mo(:,1) == year(end));
ff = FF_2016mo(FF_start:FF_end,:);
LU_start = find(LU_2016mo(:,1) == year(1));
LU_end = find(LU_2016mo(:,1) == year(end));
LU = LU_2016mo(LU_start:LU_end,:);

%% LU 2006 (other record)

LU_higher = csvread('LR_landUse_JDedits_1750-2016.csv');

% interpolate to monthly
LU_2006mo(:,1) = (1765:(1/12):2016)';
LU_2006mo(:,2) = interp1(LU_higher(:,1),LU_higher(:,12),LU_2006mo(:,1));
LU_start = find(LU_2006mo(:,1) == year(1)); % start at 1850
LUhigher = LU_2006mo(LU_start:end,:); % end at end year


% leaving extratropical land use out of this
% LUex_2000 = csvread('landUseExtra_1800-2000.csv'); % from Houghton
% % TODO: extend LUex forward and backwards to 1765-2016 by padding w zeros
% LUex_2000mo(:,1) = (1800:(1/12):2000)';
% LUex_2000mo(:,2) = interp1(LUex_2000(:,1),LUex_2000(:,2),LUex_2000mo(:,1));
% LUex_start = find(LUex_2000mo(:,1) == year2(1)); % start at 1850
% LUex_2000mo = LUex_2000mo(LUex_start:end,:); % end at end year
% LUex_end = find(LUex(:,1) == year2(end));
% LUex = LUex(1:LUex_end,:);