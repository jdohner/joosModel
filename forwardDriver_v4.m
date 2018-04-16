% forward v4
% april 13, 2018
% julia dohner
% 
% forward model of atmospheric co2 based on joos et al. (1996) ocean and
% land uptake models

clear all

LU = 1; % 1 for high land use, 2 for low land use
start_year = 1850;
end_year = 2020;
ts = 12;
dt = 1/ts;

year = start_year:dt:end_year;

% ff data - load_fossil2.m
% LU data - BP_extrap_CDIAC_data_2009

Aoc = 3.62E14; % surface area of ocean, m^2, from Joos 1996
c = 1.722E17; % unit converter, umol m^3 ppm^-1 kg^-1, from Joos 1996
h = 75; % mixed layer depth, m, from Joos 1996
T_const = 18.2; % surface temperature, deg C, from Joos 1996
kg = 1/9.06; % gas exchange rate, yr^-1, from Joos 1996
beta = 0.287; % fertilization factor
co2_preind = 280;

[t,r] = HILDAResponse(year);
[ff, LU, LUex] = getSourceSink3(year,ts);


%% the motherloop

for i = 1:length(year)
    
    % ocean calcs
    
    % calculate air-sea flux
    fas(i,1) = year(i);
    fas(i,2) = (kg/Aoc)*(dpCO2a(i,2) - dpCO2s(i,2)); % air-sea flux of CO2
    
    % calculate NPP perturbation
    delfnpp(i,1) = year(i);
    delfnpp(i,2) = 60*beta*log(CO2a(i,2)/278);
    
    % convolve the air-sea flux and the pulse response function (Joos 1996)
    % Note: LR convolves fas and r to do the integral calculation for
    % equation 3 in Joos. I'm trying the same thing for calculating
    % the integral in equation 16, but the worrisome difference is that the
    % integration bounds for equation 3 are t_0 to t, whereas the bounds in
    % equation 16 are negative infinity to t.
    w = conv(fas(1:i,2),r(1:i,2));
    v = conv(delfnpp(1:i,2),rdecay(1:i,2));
    
    if i < length(year2)
        % Calculate delDIC
        % Joos equation (3)
        delDIC(i+1,1) = year(i+1); 
        delDIC(i+1,2) = (c/h)*w(i)*dt; % change in DIC

        %Calculate dpCO2s from DIC - from Joos 1996
        dpCO2s(i+1,2) = (1.5568 - (1.3993E-2)*T_const)*delDIC(i+1,2) + ...
            (7.4706-0.20207*T_const)*10^(-3)*(delDIC(i+1,2))^2 - ...
            (1.2748-0.12015*T_const)*10^(-5)*(delDIC(i+1,2))^3 + ...
            (2.4491-0.12639*T_const)*10^(-7)*(delDIC(i+1,2))^4 - ...
            (1.5468-0.15326*T_const)*10^(-10)*(delDIC(i+1,2))^5;
    end
    
    % land calcs
    % equation 17 (page 412)
    % assumption that NPP depends logarithmically on atmospheric CO2
    % in units of ppm/yr
    
    delfdecay(i+1,1) = year(i+1);
    delfdecay(i+1,2) = delfnpp(i+1)-v(i)*dt; % only issue is that this doesn't go from neg infinity maybe?
    % other options for integral: syms, erf

    
    ffer(i,2) = delfnpp(i,2) + delfdecay(i,2);

    
    dtdelpCO2a(i,2) =  ff(i,2) - Aoc*fas(i,2) - delCdt(i,2) +landuse(i,2);
    
    if i < length(year2)
            dpCO2a(i+1,2) = dpCO2a(i,2) + dtdelpCO2a(i,2)/12; 
    end
    
        % calculate CO2a
    CO2a(i,1) = year(i);
    CO2a(i,2) = dpCO2a(i,2) + co2_preind;

end

% calculate final time points
% make sure sign convention is consistent

    
%% plotting