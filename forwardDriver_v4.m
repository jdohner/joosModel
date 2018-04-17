% forward v4
% april 13, 2018
% julia dohner
% 
% forward model of atmospheric co2 based on joos et al. (1996) ocean and
% land uptake models

clear all

LU = 1; % 1 for high land use, 2 for low land use
start_year = 1800;%1765;%
end_year = 2009+(7/12);%2020; %2300;
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
co2_preind = 278;

[t,r,rdecay] = HILDAResponse(year);
[ff, LU, LUex] = getSourceSink3(year,ts);
[~,~,CO2a_obs,~] = MLOinterpolate_increment2(ts,start_year,end_year);

%% initialize vectors

dtdelpCO2a = zeros(length(year),2); 
dtdelpCO2a(:,1) = year;
dpCO2a = zeros(length(year),2); %change in CO2a from preind
dpCO2a(:,1) = year; 
dpCO2s = zeros(length(dpCO2a),2); % dissolved CO2
dpCO2s(:,1) = dpCO2a(:,1);
delDIC = zeros(length(year),2); 
fas = zeros(length(year),2);
fas(:,1) = year;
delfnpp = zeros(length(year),2);
delfnpp(:,1) = year;
delfdecay = zeros(length(year),2);
delfdecay(:,1) = year;
ffer = zeros(length(year),2);
ffer(:,1) = year;
CO2a = zeros(length(year),2); 
CO2a(:,1) = year; 
CO2a(1,1) = year(1); % intialize first CO2a value
CO2a(1,2) = co2_preind;

%% the motherloop

for i = 1:length(year)-1
    
    % calculate air-sea flux
    fas(i,1) = year(i);
    fas(i,2) = (kg/Aoc)*(dpCO2a(i,2) - dpCO2s(i,2)); % air-sea flux of CO2

    % convolve the air-sea flux and the pulse response function (Joos 1996)
    w = conv(fas(1:i,2),r(1:i,2));
    v = conv(delfnpp(1:i,2),rdecay(1:i,2));
    
    if i < length(year)
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
        
        % calculate NPP perturbation
        delfnpp(i+1,2) = 60*beta*log(CO2a(i,2)/278);
        delfdecay(i+1,2) = v(i)*dt; 
        ffer(i+1,2) = delfnpp(i+1,2) - delfdecay(i+1,2);
    end
    
    % land calcs
    % equation 17 (page 412)
    % assumption that NPP depends logarithmically on atmospheric CO2
    % in units of ppm/yr

    dtdelpCO2a(i,2) =  ff(i,2) - Aoc*fas(i,2) - ffer(i,2) + LU(i,2);
    
    if i < length(year)
            dpCO2a(i+1,2) = dpCO2a(i,2) + dtdelpCO2a(i,2)/12;
            CO2a(i+1,2) = dpCO2a(i,2) + CO2a(1,2);
    end
    
%     if i > 1
%         % calculate CO2a
%         CO2a(i,1) = year(i);
%         %CO2a(i,2) = dpCO2a(i,2) + CO2a(i,2);
%     end

end

% calculate final time points
% make sure sign convention is consistent

%     
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