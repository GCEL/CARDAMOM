
###
## Minimum script for re-running DALEC outside of the CARDAMOM framework
## Author: T. L. Smallman (t.l.smallman@ed.ac.uk)
## Created: 11/11/2021, last edited: 11/11/2021
###

# Set working directory to the location which the CARDAMOM framework can be found
setwd("<path to CARDAMOM directory>")
# This is so we can read in the R functions for CARDAMOM, after this you can change the directory
source("./R_functions/load_all_cardamom_functions.r")

# Load the info file for the project you will be calling from
load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_LAB_MHMCMC/Miombo_kilwa_nhambita_1km/infofile.RData")
# Load the already processed DALEC outputs - this has the parameters and drivers but also the existing DALEC output against which you can compare any modifications
load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_LAB_MHMCMC/Miombo_kilwa_nhambita_1km/RESULTS_PROCESSED/kilwa.RData")
# If the project has more than one site within set "n" to the correct value, otherwise leave as 1
n = 1

# Set missing observation values to NA as this will be easiler for plotting, -9999 needed for for CARDAMOM itself
drivers$obs[which(drivers$obs == -9999)] = NA

###
## Modify the drivers if you wish

## You can extend the analysis time series while keeping the existing dataset intact

# Copy original drivers into a new object that we will manipulate
new_drivers = drivers

nos_loops = 0
if (nos_loops > 0) {
    # Loop drivers to allow simulation to steady state
    for (r in seq(1,nos_loops)) {
         new_drivers$met = rbind(new_drivers$met,drivers$met)
    }
    # Just looping the existing drivers doesn't work. 
    # At the very least the simulation day need to be updated so continuously increase.
    # We can achieve this by cumulative sum of the day of year variable
    new_drivers$met[,1] =  PROJECT$model$timestep_days ; new_drivers$met[,1] = cumsum(new_drivers$met[,1])
} # nos_loops

## You can manipulate the meteorology

new_drivers$met[,2] = new_drivers$met[,2] + 0 # minimum temperature (C)
new_drivers$met[,3] = new_drivers$met[,3] + 0 # maximum temperature (C)
new_drivers$met[,4] = new_drivers$met[,4] + 0 # mean daily shortwave radiation (MJ/m2/day)
new_drivers$met[,5] = new_drivers$met[,5] + 0 # atmospheric CO2 concentration (ppm)
#new_drivers$met[,6] = new_drivers$met[,6] + 0 # Julian day of year - DON'T MESS WITH THIS
new_drivers$met[,7] = new_drivers$met[,7] + 0 # mean precipitation rate (kgH2O/m2/s)

# ACM2 enabled
new_drivers$met[,14] = new_drivers$met[,14] + 0 # mean temperature (C)
new_drivers$met[,15] = new_drivers$met[,15] + 0 # mean wind speed (m/s)
new_drivers$met[,16] = new_drivers$met[,16] + 0 # mean vapour pressure deficit (Pa)

# GSI models only
#new_drivers$met[,10] = new_drivers$met[,10] + 0 # 21 day rolling mean avg max temperature (C)
#new_drivers$met[,11] = new_drivers$met[,11] + 0 # 21 day rolling mean avg day length (seconds)
#new_drivers$met[,12] = new_drivers$met[,12] + 0 # 21 day rolling mean avg vapour pressure deficit (Pa)

## You can manipulate disturbance drivers

new_drivers$met[,8] = new_drivers$met[,8] + 0 # harvested fraction (0-1)
new_drivers$met[,9] = new_drivers$met[,9] + 0 # burnt fraction (0-1)

# GSI models only
#new_drivers$met[,13] = new_drivers$met[,13] + 0 # forest management type (catagorical)

## Update timing information if needed  

# Determine number of steps per year
steps_per_year = dim(drivers$met)[1] / length(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year))
# Determine the number of years in the simulation in total
nos_years = dim(new_drivers$met)[1] / steps_per_year
# Number years needed for output from dalec.so, make sure that the original PROJECT time period is updated
analysis_years = length(c(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)))
new_PROJECT = PROJECT
if (analysis_years != nos_years) {
    new_PROJECT$end_year = as.numeric(PROJECT$end_year) + (nos_years - analysis_years)
}

##
# Now run the actual model

# run subsample of parameters for full results / propogation
soil_info = c(drivers$top_sand,drivers$bot_sand,drivers$top_clay,drivers$bot_clay)
C_cycle = simulate_all(n,new_PROJECT,PROJECT$model$name,new_drivers$met,parameters[1:PROJECT$model$nopars[n],,],
                          drivers$lat,PROJECT$ctessel_pft[n],PROJECT$parameter_type,
                          PROJECT$exepath,soil_info)

# Then compare original states_all with C_cycle

