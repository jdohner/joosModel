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
 
function [ff,LU] = getSourceData(year,ts)
 
addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/joosModel/co2_forward_data_2016/JDfiles'));
 
FF_2016 = csvread('FF_Boden_2016.csv'); % in gigatons/year
LU_2016 = csvread('LU_Houghton_2016.csv'); % in gigatons/year
 
d = 1/2.31; % 1 ppm CO2 = 2.31 gton CO2
 
% interpolate to monthly
FF_2016mo(:,1) = (1765:(1/ts):FF_2016(end,1))';
FF_2016mo(:,2) = (interp1(FF_2016(:,1),FF_2016(:,2),FF_2016mo(:,1)));
FF_2016mo(:,2) = FF_2016mo(:,2)*d; %convert to ppm
 
LU_2016mo(:,1) = (1765:(1/ts):LU_2016(end,1))';
LU_2016mo(:,2) = interp1(LU_2016(:,1),LU_2016(:,2),LU_2016mo(:,1));
LU_2016mo(:,2) = LU_2016mo(:,2)*d; % convert to ppm
 
% shorten datasets to match time frame of year vector
FF_start = find(FF_2016mo(:,1) == year(1));
FF_end = find(FF_2016mo(:,1) == year(end));
ff = FF_2016mo(FF_start:FF_end,:);
 
LU_start = find(LU_2016mo(:,1) == year(1));
LU_end = find(LU_2016mo(:,1) == year(end));
LU = LU_2016mo(LU_start:LU_end,:);


