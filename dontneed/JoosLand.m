% joos land model
%
% author: julia dohner
%
% april 12, 2018
%
% code for land model from Joos et al. (1996)

% equation A.3. Biosphere decay response function (page 416)

start_year = 1850;
end_year = 2020;
dt = 1/12;

year = start_year:dt:end_year;

% sanity check: rdecay(t = 5yrs) = 0.062610
rdecay(:,2) = 0.70211   *exp(-0.35t)+0.013414 ...
                        *exp(-t/20)-0.71846 ...
                        *exp(-55*t/120)+0.0029323 ...
                        *exp(-t/100);

% equation 17 (page 412)
% assumption that NPP depends logarithmically on atmospheric CO2
% in units of ppm/yr
del_fnpp(:,2) = 60*beta*ln(CO2a(:,2)/278);

% equation 16 (page 412)
% additional NPP (del_npp) and decay (del_fdecay) due to fert. alone
del_fnpp; % defined in equation 17
del_fdecay = 1; %integral from neg infinity to t etc etc
ffer(:,2) = del_fnpp - del_fdecay;


% equation 4 (page 402)
% relationship between atmospheric co2 perturbations and carbon emissions
% del_pco2a = atmos. CO2 perturbations
% e = carbon emissions
% units in ppm/yr
diff(del_pco2a) = e(:,2) - Aoc*fas(:,2);

% equation 14 (page 410)
% the emission term e may be subdivided into anthrop emissions (ff, LU) and
% an additional biospheric sink term ffer
e = e_anth - ffer;


