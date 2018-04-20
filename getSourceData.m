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
 
function [ff,LU] = getSourceData(year,ts,FF_data,LU_data)
    
% comment out the two lines below if using files being fed in
FF_data = csvread('FFdata_Boden_2016.csv'); % in gigatons/year
LU_data = csvread('LUdata_Houghton_2016.csv'); % in gigatons/year
 
d = 1/2.31; % 1 ppm CO2 = 2.31 gton CO2
 
% interpolate to monthly
FFmo(:,1) = (1765:(1/ts):FF_data(end,1))';
FFmo(:,2) = (interp1(FF_data(:,1),FF_data(:,2),FFmo(:,1)));
FFmo(:,2) = FFmo(:,2)*d; %convert to ppm
 
LUmo(:,1) = (1765:(1/ts):LU_data(end,1))';
LUmo(:,2) = interp1(LU_data(:,1),LU_data(:,2),LUmo(:,1));
LUmo(:,2) = LUmo(:,2)*d; % convert to ppm
 
% shorten datasets to match time frame of year vector
FF_start = find(FFmo(:,1) == year(1));
FF_end = find(FFmo(:,1) == year(end));
ff = FFmo(FF_start:FF_end,:);
 
LU_start = find(LUmo(:,1) == year(1));
LU_end = find(LUmo(:,1) == year(end));
LU = LUmo(LU_start:LU_end,:);


