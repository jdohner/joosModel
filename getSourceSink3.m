% getSourceSink2.m
%
% author: Lauren Rafelski, modified by Julia Dohner
% January 23, 2018
%
% reorganized version of getsourcesink_scale3 from LR code. Allows for
% forward model running to 2016
%
% this version joins records at latest timepoint and extends through 2016
% but will cut off at end year

% LU in ppm/year
% ff in ppm/year


function [ff,LU,LUex] = getSourceSink3(year2, ts);

    
year3 = (year2(1,1):(1/ts):2016)'; % making full data vector thru to 2016
    
% load data thru 2009, 2006, 2000 - ff is monthly, lu is annual

% both vectors begin in 1800
LU_2006 = csvread('landUse_1800-2006.csv');
LUex_2000 = csvread('landUseExtra_1800-2000.csv');
load fossilFuel_1751-2009.mat;
FF_2009 = ff1; % already monthly resolution

% shortening ff vector to begin at start_year (have data back thru 1700)
FF_start = find(FF_2009(:,1) == year2(1));
FF_2009 = FF_2009(FF_start:end,:);

% interpolate to monthly
LU_2006mo(:,1) = (1800:(1/12):2006)';
LU_2006mo(:,2) = interp1(LU_2006(:,1),LU_2006(:,12),LU_2006mo(:,1));
LU_start = find(LU_2006mo(:,1) == year2(1)); % start at 1850
LU_2006mo = LU_2006mo(LU_start:end,:); % end at end year

LUex_2000mo(:,1) = (1800:(1/12):2000)';
LUex_2000mo(:,2) = interp1(LUex_2000(:,1),LUex_2000(:,2),LUex_2000mo(:,1));
LUex_start = find(LUex_2000mo(:,1) == year2(1)); % start at 1850
LUex_2000mo = LUex_2000mo(LUex_start:end,:); % end at end year


% extend land use with last value thru to latest year
% if want to emulate LR exactly, use this rather than filling in obs to 2016
LU(:,1) = year3;
LU(1:length(LU_2006mo),2) = LU_2006mo(:,2);
LU(length(LU_2006mo)+1:end,2) = LU_2006mo(end,2);

    
% extend extratropical land use to end_year with zeros thru to latest year
LUex(:,1) = year3; 
LUex(1:length(LUex_2000mo),2) = LUex_2000mo(:,2);
LUex(length(LUex_2000mo)+1:end,2) = 0;

%% extend vectors to 2016

% load extended data thru 2016 - all annual
FF_2016 = csvread('fossilFuel_1959-2016.csv'); % in gigatons/year
% 1 ppm CO2 = 2.31 gton CO2
d = 1/2.31; % gigaton to ppm conversion factor

% this one isn't working for csvread, so reading in as text file
%landUse_2016 = csvread('landUse_1959-2016.csv'); 
fid = fopen('landUse_1959-2016.txt');
C = textscan(fid,'%f %f', 'delimiter','\t');
fclose(fid);
LU_2016(:,1) = C{1};
LU_2016(:,2) = C{2};

% fossil fuel 2016 to monthly (arrives as 1959-2016)
month_2016 = 1959:(1/12):2016;
FF_2016mo_0 = (interp1(FF_2016(:,1),FF_2016(:,2),month_2016)).';
FF_2016mo(:,1) = month_2016;
FF_2016mo(:,2) = FF_2016mo_0*d; %convert to ppm

% land use 2016 to monthly
LU_2016mo_0 = (interp1(LU_2016(:,1),LU_2016(:,2),month_2016)).';
LU_2016mo(:,1) = month_2016;
LU_2016mo(:,2) = LU_2016mo_0*d;

% extend (making full length vectors (1800-2016)) and patch

clear ff; % get rid of dimension mismatch since looks like below pulls
% from original data source (FF_2009) anyways
% Houghton record begins in 1959, combine two records at 1959
%i_2009 = find(FF_2009(:,1) == 2009); % same for ff, landuse, extraLU
i_2009 = length(FF_2009);
enddate_2009 = FF_2009(end,1);
startdate_2016 = find(FF_2016mo(:,1) == enddate_2009);
ff(:,1) = year3;
%j_2009 = find(FF_2009(:,1) == 2009);
ff(1:i_2009,2) = FF_2009(:,2);
ff(i_2009+1:end,2) = FF_2016mo(startdate_2016+1:end,2);


%% shorten datasets to match time frame of year vector

FF_end = find(ff(:,1) == year2(end));
ff = ff(1:FF_end,:);
LU_end = find(LU(:,1) == year2(end));
LU = LU(1:LU_end,:);
LUex_end = find(LUex(:,1) == year2(end));
LUex = LUex(1:LUex_end,:);

