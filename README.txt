README.txt for joosModel
Model for future CO2 levels to be used to answer "When will the Keeling Curve bend?"
April 19, 2018

Project description:
This code runs the model for oceanic and biospheric CO2 uptake from the 1996 paper by Joos et al. The code uses the Joos model to calculate the change in atmospheric CO2 given fossil fuel and land use emission data. This code will be used to predict future atmospheric CO2 levels given future CO2 emission scenarios. For the current project with Dr. Ralph Keeling, Dr. Armin Schwartzman, Dr. David Victor, Dr. Ahmed Abdulla, Dr. Fabian Telschlow, and Julia Dohner, the future CO2 levels calculated from this code will ultimately be used to determine when we can detect bifurcations in the trajectories of atmospheric CO2 levels of different emission scenarios from the business-as-usual trajectory.


Prerequisites:
To run this code, you will need MATLAB installed on your computer. If choosing to run this code with fossil fuel and land use emission records other than those running from 1751-2016 provided by Boden et al. and Houghton et al., respectively, you will need data files of your emissions data of choice saved in the joosModel folder.


Running tests:
In its original form, the code should yield two plots identical to those saved as MATLAB figures (extension .fig) in the joosModel folder. The original code (and resultant figures) were run with the following conditions:
start_year = 1756;
end_year = 2016;
ts = 12;
CO2_preind = 283;
c1 = 0.85;
FF data = Boden thru 2016
LU data = Houghton thru 2016

Tweaking parameters:
To achieve an optimal fit with the observed CO2 record, the parameters CO2_preind and c1 should be tweaked. C1 will adjust the strength of the land and ocean sinks (a fair tactic because both contain a degree of uncertainty). The Joos et al. (1996) paper sets CO2_preind as 278, but it is clear that adjusting this value results in a better fit to the observed CO2 record.

Miscellaneous notes on the code:

1. Extratropical land use is not included in this model (though it is included as a "low land use" scenario in the Rafelski et al. (2009) paper.

2. To change the timestep or time boundaries on the observed CO2 record, uncomment the line calling the function getObservedCO2. See notes on CO2 record in "Information on data" section below. CO2a_obs record has monthly averages recorded in the middle of each month at increments of 1/12 yr starting at 1/24th of a year, but then switches in 2011 to being recorded on the first of each month. This should not be an issue.

3. The Boden and Houghton FF and LU emissions files are in gigaton of carbon per year, and are converted to ppm/yr using the conversion factor d in the driver. If the input files you feed in are already in ppm/yr, be sure to comment this out. If your input files are not in ppm/yr, be sure to include a conversion factor.


Authors:
Julia Dohner, with code adapted from Lauren Rafelski.


Information on data:

FFdata_Boden_2016.csv
Annual data from 1751-2016 in gigaton carbon/yearFossil fuel combustion and cement production emissions: Boden, T. A., Marland, G., and Andres, R. J.: Global, Regional, and National Fossil-Fuel CO2 Emissions, Oak Ridge National Laboratory, U.S. Department of Energy, Oak Ridge, Tenn., U.S.A., doi 10.3334/CDIAC/00001_V2017, 2017; available at: http://cdiac.ess-dive.lbl.gov/trends/emis/overview_2014.html. 
Downloaded from http://www.globalcarbonproject.org/carbonbudget/17/data.htmLUdata_Houghton_2016.csv
Annual data from 1751-2016 in gigaton carbon/year. Value of 0 added at 1751 to interpolate backwards. Actual data begins in 1850.Land use emissions data: Houghton, R. A. and Nassikas, A. A.: Global and regional fluxes of carbon from land use and land cover change 1850-2015, Global Biogeochemical Cycles, 31, 456-472, 2017.Downloaded from http://www.globalcarbonproject.org/carbonbudget/17/data.htmdataObservedCO2.matCO2 record: Received from Lauren Rafelski in package of code for paper Rafelski et al. (2009). From the paper: "We use a global CO2 record based on a combination of the ice core record from Law Dome before 1958 (Etheridge et al.,1996; MacFarling Meure et al., 2006) and a seasonally detrended arithmetic average of monthly air measurements from Maura Loa and the South Pole from the Scripps CO2 program after 1958. The records were combined without adjustment. The ice core data are approximated to monthly resolution using a spline with a standard error sigma of 0.6 ppm CO2 (Rafelski et al., 2009)."

Note: if you want to have the observed CO2 record in a timestep other than monthly or use a time frame other than 1765-2016, the CO2 record can be re-interpolated using the function getObservedCO2.m included in the joosModel file folder. Simply feed the desired tilmestep (ts), start year and end year. Otherwise, a datafile of the saved output variables from this function for ts = 12, start year = 1765 and end year = 2016 will be used in the model driver code.