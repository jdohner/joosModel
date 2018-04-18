% 3/13/08: changed CO2 dataset to spline fit with sigma =0.6 to capture
% 1940s plateau
% 4/2/08: changed CO2 dataset to updated ice core data, sigma=0.6
% 4/24/08: add in an option to decrease ice core data by 2 ppm; change to
% linear interpolation for this option
% 1/10/11: add updated CO2 dataset through 2010
        
function [annincMLOSPO,dpCO2a,co2_combine_trunc,co2_preind] = getObservedCO2(ts,start_year,end_year)

addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/joosModel/co2_forward_data_2016'));

addpath(genpath(...
    '/Users/juliadohner/Documents/MATLAB/joosModel/co2_forward_data_2016/JDfiles'));

load co2_2011_2.mat

dt = 1/ts; 

% processing data for co2 up to 2011
% loads in mlospo_meure (4427x2), 1640-Feb 2010, monthly
meure_years = mlospo_meure(:,1);
mlostart = mlospo_meure(1,1);
mloend = mlospo_meure(end,1);
meure_CO2 = mlospo_meure(:,2);

% Create new time array
meureInterp_years = mlostart:1/ts:mloend;

% Do interpolation
meureInterp_CO2 = interp1(meure_years,meure_CO2,meureInterp_years,'spline');

MLOSPOiceinterp(:,1) = meureInterp_years; % ends Feb 2010
MLOSPOiceinterp(:,2) = meureInterp_CO2;

% seems super similar to merged co2 2016 record below

%% everything past where co2_2011 ends

% starts at year 1, incremented by year
CO2_2016 = csvread('mergedCO2_2016.csv');
year_2016 = (CO2_2016(1,1):1/ts:CO2_2016(end,1))';
CO2_2016mo(:,1) = year_2016;
CO2_2016mo(:,2) = (interp1(CO2_2016(:,1),CO2_2016(:,2),year_2016)).';

% joinYear is 2010 + 2/12 (i.e. March 2010)
joinYear = MLOSPOiceinterp(end,1)+(1/24); % this brings us from 2010+(3/24) (i.e. Feb 2010) to March 2010 (starting at 2010+2/12)
i = find(CO2_2016mo(:,1) >= joinYear,1);
year_full(:,1) = [MLOSPOiceinterp(:,1) ; CO2_2016mo(i:end,1)];

% starts at beginning of MLOSPO, ends at end of 2016 record
co2_combine(:,1) = year_full; 
co2_combine(:,2) = [MLOSPOiceinterp(1:end,2); CO2_2016mo(i:end,2)];

co2_preind = mean(co2_combine(1:1000,2));

%% Calculate CO2 increment with monthly resolution, in ppm/year
% n = 7 is 7/1958, last value is 7/2005

for n = ((ts/2)+1):(length(co2_combine)-(ts/2))
    annincMLOSPO(n,1) = co2_combine(n,1);
    annincMLOSPO(n,2) = co2_combine(n+(ts/2),2) - co2_combine(n-(ts/2),2);
end

year = start_year:dt:end_year;
i1 = find(co2_combine(:,1) >= start_year,1);
j1 = find(co2_combine(:,1) >= end_year);
co2_trunc = co2_combine(i1:j1,:);

%% Calculate change in atmospheric concentration
dpCO2a(:,1) = co2_trunc(:,1); 
dpCO2a(:,2) = co2_trunc(:,2)-co2_trunc(1,2);

n = find(co2_combine(:,1) >= end_year,1);
m = find(co2_combine(:,1) >= start_year,1);
co2_combine_trunc = co2_combine(m:n,:);