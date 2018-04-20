% joosModelDriver.m
% 
% Author: Julia dohner
% April 13, 2018
%
% Forward model of atmospheric CO2 based on Joos et al. (1996) ocean and
% land CO2 uptake models.

clear all

%% adjustables are between the percent sign borders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data files to load in - must be in ppm
FF_data = csvread('dataFF_Boden_2016.csv'); % in gigatons/year
LU_data = csvread('dataLU_Houghton_2016.csv'); % in gigatons/year
% convert to ppm
d = 1/2.31; % 1 ppm CO2 = 2.31 gton CO2
FF_data(:,2) = FF_data(:,2)*d;
LU_data(:,2) = LU_data(:,2)*d;

start_year = 1765;
end_year = 2016;
ts = 12; % timesteps per year
dt = 1/ts;
year = start_year:dt:end_year;

beta = 0.287; % fertilization factor
CO2_preind = 283; % 278 in Joos, but tweaking to match observed record
c1 = 0.85; % sinks scaling factor


Aoc = 3.62E14; % surface area of ocean, m^2, from Joos 1996
c = 1.722E17; % unit converter, umol m^3 ppm^-1 kg^-1, from Joos 1996
h = 75; % mixed layer depth, m, from Joos 1996
T_const = 18.2; % surface temperature, deg C, from Joos 1996
kg = 1/9.06; % gas exchange rate, yr^-1, from Joos 1996

[t,r,rdecay] = HILDAresponse(year);
[ff, LU] = getSourceData(year,ts,FF_data,LU_data);
load dataObservedCO2.mat; % loads in dtdelpCO2a_obs,dpCO2a_obs,CO2a_obs
%[dtdelpCO2a_obs,dpCO2a_obs,CO2a_obs] = getObservedCO2(ts,start_year,end_year);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize vectors
blankVec(:,1) = year;
blankVec(:,2) = 0;
dtdelpCO2a = blankVec; % annual change in CO2atm (ppm/yr)
dpCO2a = blankVec; % change in CO2atm from preindustrial
dpCO2s = blankVec; % change in dissolved CO2 from preindustrial
delDIC = blankVec; % change in dissolved inorganic carbon
fas = blankVec; % air-sea flux
delfnpp = blankVec; % change in net primary production from added CO2
delfdecay = blankVec; % change in organic matter decay from added CO2
ffer = blankVec; % total perturbation (delfnpp + delfdecay) from added CO2
CO2a = blankVec; % atmospheric CO2
CO2a(1,1) = year(1); % intialize first CO2a value
CO2a(1,2) = CO2_preind;
sumCheck = blankVec;


%% loop to calculate change in atmospheric co2 growth rate 
% calculates values at each month by summing sources and sinks

for i = 1:length(year)-1
    
    % calculate air-sea flux
    fas(i,2) = (kg/Aoc)*(dpCO2a(i,2) - dpCO2s(i,2)); 

    w = conv(fas(1:i,2),r(1:i,2)); % Eq. 3 (Joos '96)
    v = conv(delfnpp(1:i,2),rdecay(1:i,2)); % Eq. 16 (Joos '96)

    % calculate change in DIC
    delDIC(i+1,2) = (c/h)*w(i)*dt; % Eq. 3 (Joos '96)

    %Calculate dpCO2s from DIC - Eq. 6b (Joos '96)
    dpCO2s(i+1,2) = (1.5568 - (1.3993E-2)*T_const)*delDIC(i+1,2) + ...
        (7.4706-0.20207*T_const)*10^(-3)*(delDIC(i+1,2))^2 - ...
        (1.2748-0.12015*T_const)*10^(-5)*(delDIC(i+1,2))^3 + ...
        (2.4491-0.12639*T_const)*10^(-7)*(delDIC(i+1,2))^4 - ...
        (1.5468-0.15326*T_const)*10^(-10)*(delDIC(i+1,2))^5;

    % calculate NPP perturbation
    delfnpp(i+1,2) = 60*beta*log(CO2a(i,2)/CO2_preind); % Eq. 17 (Joos '96)
    delfdecay(i+1,2) = v(i)*dt; % Eq. 16 (Joos '96)
    ffer(i+1,2) = delfnpp(i+1,2) - delfdecay(i+1,2); % Eq. 16 (Joos '96)

    % calculate change in atmospheric CO2
    dtdelpCO2a(i,2) =  ff(i,2) + LU(i,2) - c1*(Aoc*fas(i,2) + ffer(i,2)) ; % Eq. 4 (Joos '96)
    dpCO2a(i+1,2) = dpCO2a(i,2) + dtdelpCO2a(i,2)/12; 
    CO2a(i+1,2) = dpCO2a(i,2) + CO2a(1,2); 
    sumCheck(i,2) = dtdelpCO2a(i,2)+Aoc*fas(i,2)+ffer(i,2)-LU(i,2)-ff(i,2);

end

% calculate final time points
fas(end,2) = (kg/Aoc)*(dpCO2a(end,2) - dpCO2s(end,2)); 
dtdelpCO2a(end,2) =  ff(end,2) + LU(end,2) - c1*(Aoc*fas(i,2) + ffer(i,2)); 
     
%% plotting

figure
plot(ff(:,1),ff(:,2),fas(:,1),Aoc*fas(:,2),ffer(:,1),ffer(:,2),LU(:,1),LU(:,2),dtdelpCO2a(:,1),dtdelpCO2a(:,2))
legend('Fossil fuel','Air-sea flux','Land sink','Land use','Change in atmospheric CO2','location','northwest')
ylabel('ppm/yr')
xlabel('year')
grid

figure
plot(CO2a_obs(:,1),CO2a_obs(:,2),CO2a(:,1),CO2a(:,2));
legend('Observed atmospheric CO2','Calculated atmospheric CO2','location','northwest')
ylabel('ppm')
xlabel('year')
grid