% file HILDAResponse.m
% 
% author Julia Dohner, borrowing from Lauren Rafelski
% 
% note Be sure to run defaults.m in outer folder before running this code
% 
% brief Calculates the response "r" from the HILDA model (A.2.2 in Joos
% 1996)
% 
% code from jooshildascale_annotate2.m (aka ocean model driver)

function [t,r,rdecay]= HILDAresponse(year)

% allocating space for t and r matrices before for loop

% tracer concentration (Joos 1996 A.1. pg 415)
% one column matrix
t = NaN(length(year), 1); 

% pulse response function (Joos 1996 A.1. pg 415)
% two column matrix
r = NaN(length(year), 2); 

% Response function to calculate ocean uptake
% year comes from load_fossil2 here
% A.2.2
for i = 1:length(year)
    % each value in first column of t is the year at index i minus first
    % year
     t(i,1) = year(i) - year(1);
     r(i,1) = t(i,1);
    %Calculate response function based on HILDA equation in Joos 1996
     if t(i,1) == 0
         r(i,2) = 1;
     elseif t(i,1) <= 2
         r(i,2)= (1/0.95873)*(0.12935+0.21898*exp(-t(i,1)/0.034569)+0.17003*exp(-t(i,1)/0.26936)...
             +0.24071*exp(-t(i,1)/0.96083)+0.24093*exp(-t(i,1)/4.9792));
     else
         r(i,2) = (1/0.95873)*(0.022936+0.24278*exp(-t(i,1)/1.2679)+0.13963*exp(-t(i,1)/5.2528)...
                +0.089318*exp(-t(i,1)/18.601)+0.037820*exp(-t(i,1)/68.736)...
                +0.035549*exp(-t(i,1)/232.3));
     end
end

% A.3 biosphere decay response function
for i = 1:length(year)
    rdecay(i,1) = t(i,1);
    
    if t(i,1) == 0
         rdecay(i,2) = 0; % set first value = 0
    else
        rdecay(i,2) = 0.70211   *exp(-0.35*t(i,1))+0.013414 ...
                        *exp(-t(i,1)/20)-0.71846 ...
                        *exp(-55*t(i,1)/120)+0.0029323 ...
                        *exp(-t(i,1)/100);
    end
end

  
end