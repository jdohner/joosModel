% forward v4
% april 13, 2018
% julia dohner
% 
% forward model of atmospheric co2 based on joos et al. (1996) ocean and
% land uptake models

clear all

start_year = 1850;
end_year = 2009+(7/12);
ts = 12;
dt = 1/ts;
year = start_year:dt:end_year;

Aoc = 3.62E14; % surface area of ocean, m^2, from Joos 1996
c = 1.722E17; % unit converter, umol m^3 ppm^-1 kg^-1, from Joos 1996
h = 75; % mixed layer depth, m, from Joos 1996
T_const = 18.2; % surface temperature, deg C, from Joos 1996
kg = 1/9.06; % gas exchange rate, yr^-1, from Joos 1996
beta = 0.287; % fertilization factor
co2_preind = 278;

[t,r,rdecay] = HILDAResponse(year);
%[ff, LU,~] = getSourceSink4(year,ts);
[ff, LU,~] = getSourceSink3(year,ts);
[~,~,CO2a_obs,~] = getObservedCO2(ts,start_year,end_year);

% NOTE: load files for fossil fuel, land use and extratropical land use 
% emissions here if using forward projections

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
CO2a(1,2) = co2_preind;
sumCheck = blankVec;


%% the motherloop

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
    delfnpp(i+1,2) = 60*beta*log(CO2a(i,2)/278); % Eq. 17 (Joos '96)
    delfdecay(i+1,2) = v(i)*dt; % Eq. 16 (Joos '96)
    ffer(i+1,2) = delfnpp(i+1,2) - delfdecay(i+1,2); % Eq. 16 (Joos '96)

    %ff(i,2)
    % calculate change in atmospheric CO2
    dtdelpCO2a(i,2) =  ff(i,2) ;%+ LU(i,2);%- Aoc*fas(i,2) - ffer(i,2) ; % Eq. 4 (Joos '96)
    dpCO2a(i+1,2) = dpCO2a(i,2) + dtdelpCO2a(i,2)/12; 
    CO2a(i+1,2) = dpCO2a(i,2) + CO2a(1,2); 
    sumCheck(i,2) = dtdelpCO2a(i,2)+Aoc*fas(i,2)+ffer(i,2)-LU(i,2)-ff(i,2);

end

% calculate final time points
fas(end,2) = (kg/Aoc)*(dpCO2a(end,2) - dpCO2s(end,2)); 
%dtdelpCO2a(end,2) =  ff(end,2) - Aoc*fas(end,2) - ffer(end,2) + LU(end,2); 
     
%% plotting

figure
plot(ff(:,1),ff(:,2),fas(:,1),Aoc*fas(:,2),ffer(:,1),ffer(:,2),LU(:,1),LU(:,2),dtdelpCO2a(:,1),dtdelpCO2a(:,2))
legend('fossil fuels','air-sea flux','additional land sink','land use','change in atmospheric co2','location','northwest')
ylabel('ppm/yr')
xlabel('year')

figure
plot(CO2a_obs(:,1),CO2a_obs(:,2),CO2a(:,1),CO2a(:,2));
legend('observed atmospheric co2','calculated atmospheric co2','location','northwest')
ylabel('ppm')
xlabel('year')