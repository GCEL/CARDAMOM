
###
## Minimum script for re-running DALEC outside of the CARDAMOM framework
## Author: T. L. Smallman (t.l.smallman@ed.ac.uk)
## Created: 11/11/2021, last edited: 11/11/2021
###

# Set working directory to the location which the CARDAMOM framework can be found
setwd("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/")
# This is so we can read in the R functions for CARDAMOM, after this you can change the directory
source("./R_functions/load_all_cardamom_functions.r")

# Load the info file for the project you will be calling from
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/CZO_lai_EDC_TendyMet_longterm/infofile.RData")
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A3.C1.D2.F2.H2.P1.#_MHMCMC/CZO_lai_EDC_TendyMet_longterm/infofile.RData")
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A3.C1.D2.F2.H2.P1.#_MHMCMC/CZO_fapar_EDC_TendyMet_longterm/infofile.RData")
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A3.C1.D2.F2.H2.P1.#_MHMCMC/CZO_lai_fapar_EDC_TendyMet_longterm/infofile.RData")
load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A3.C3.H2.M1.#_MHMCMC/ATEC/infofile.RData")
# If the project has more than one site within set "n" to the correct value, otherwise leave as 1
n = 1
# Load the already processed DALEC outputs - this has the parameters and drivers but also the existing DALEC output against which you can compare any modifications
load(paste(PROJECT$results_processedpath,"/",PROJECT$sites[n],".RData", sep=""))

# Set missing observation values to NA as this will be easier for plotting, -9999 needed for for CARDAMOM itself
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

#new_drivers$met[56,8] = new_drivers$met[56,8] + 0.8 # harvested fraction (0-1)
#new_drivers$met[,9] = new_drivers$met[,9] + 0 # burnt fraction (0-1)
#new_drivers$met[56,13] = 4 # forest management type (catagorical)

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
#parameters[11,,] = 0.3 #parameters[11,,]*2
#parameters[4,,] = parameters[4,,]*1.5
#parameters[5,,] = parameters[5,,]*0.5
#parameters[7,,] = 0.03
#parameters[13,,] = 0.21
#parameters[35,,] = 0.99
# run subsample of parameters for full results / propogation
soil_info = c(drivers$top_sand,drivers$bot_sand,drivers$top_clay,drivers$bot_clay)
C_cycle = simulate_all(n,new_PROJECT,PROJECT$model$name,new_drivers$met,parameters[1:PROJECT$model$nopars[n],,],
                       drivers$lat,PROJECT$ctessel_pft[n],PROJECT$parameter_type,
                       PROJECT$exepath,soil_info)

# Read in the overall observations file
#czo = read.csv("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/SENSE_REP_2023/CZO_lai_obs/CZO_timeseries_obs.csv")
#czo = read.csv("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/SENSE_REP_2023/CZO_lai_obs/soaprootRawHalfHourly.csv", header=FALSE)
#czo_par_in = rollapply(czo$V11, by = 48, width=48, FUN=sum) ; czo_par_out = rollapply(czo$V12, by = 48, width=48, FUN=sum)
# Then compare original states_all with C_cycle

par(mfrow=c(2,4))
plot(apply(states_all$gpp_gCm2day,2,median), type="l", lwd=3, ylim=c(0,15), main="GPP")
lines(apply(C_cycle$gpp_gCm2day,2,median), col="green", lwd=3)
#plotCI(czo$GPP_gCm2day, uiw = czo$GPP_unc_gCm2day, add=TRUE, ylim=c(0,9), col="red")
#plot(apply(states_all$ET_kgH2Om2day,2,median), type="l", lwd=3, ylim=c(0,8), main="ET")
#lines(apply(C_cycle$ET_kgH2Om2day,2,median), col="green", lwd=3)
plot(apply(states_all$foliage_gCm2,2,median), type="l", lwd=3, ylim=c(0,250), main="Fol")
lines(apply(C_cycle$foliage_gCm2,2,median), col="green", lwd=3)
#plotCI(czo$Evap_kgH2Om2day, uiw = czo$Evap_unc_kgH2Om2day, add=TRUE, col="red")
plot(apply(states_all$SurfWater_kgH2Om2,2,median), type="l", lwd=3, ylim=c(0,150), main="SurfWater")
lines(apply(C_cycle$SurfWater_kgH2Om2,2,median), col="green", lwd=3)
plot(apply(states_all$lai_m2m2,2,median), type="l", lwd=3, ylim=c(0,10), main="LAI")
lines(apply(C_cycle$lai_m2m2,2,median), col="green", lwd=3)
plotCI(drivers$obs[,3], uiw = drivers$obs[,4], add=TRUE)
#plotCI(czo$LAI_m2m2, uiw = czo$LAI_unc_m2m2, add=TRUE, col="red")
#plot(apply(states_all$rauto_gCm2day+states_all$rhet_litter_gCm2day+states_all$rhet_som_gCm2day,2,median), type="l", lwd=3, ylim=c(0,16))
#lines(apply(C_cycle$rauto_gCm2day+C_cycle$rhet_litter_gCm2day+C_cycle$rhet_som_gCm2day,2,median), col="green", lwd=3)
#plotCI(czo$Reco_gCm2day, uiw = czo$Reco_unc_gCm2day, add=TRUE, ylim=c(0,10), col="red")
plot(apply(states_all$RootDepth_m,2,median), type="l", lwd=3, ylim=c(0,max(states_all$RootDepth_m, na.rm=TRUE)), main="RootDepth")
lines(apply(C_cycle$RootDepth_m,2,median), col="green", lwd=3)
plot(apply(states_all$harvest_gCm2day,2,median), type="l", lwd=3, ylim=c(0,350), main="Harvest")
lines(apply(C_cycle$harvest_gCm2day,2,median), col="green", lwd=3)
plot(apply(states_all$StorageOrgan_gCm2,2,median), type="l", lwd=3, ylim=c(0,350), main="StorageOrgan")
lines(apply(C_cycle$StorageOrgan_gCm2,2,median), col="green", lwd=3)
plot(apply(states_all$som_gCm2,2,median), type="l", lwd=3, main="Soil C")
lines(apply(C_cycle$som_gCm2,2,median), col="green", lwd=3)

par(mfrow=c(1,1))
plot(apply(states_all$alloc_foliage_gCm2day,2,median), type="l", lwd=3, main="Allocation pattern", col="green", ylim=c(0,4))
lines(apply(C_cycle$alloc_wood_gCm2day,2,median), col="brown", lwd=3)
lines(apply(C_cycle$alloc_labile_gCm2day,2,median), col="blue", lwd=3)
lines(apply(C_cycle$alloc_roots_gCm2day,2,median), col="red", lwd=3)
lines(apply(C_cycle$alloc_StorageOrgan_gCm2day,2,median), col="black", lwd=3)
lines(apply(C_cycle$alloc_autotrophic_gCm2day,2,median), col="yellow", lwd=3)
par(new=TRUE) ; plot(apply(C_cycle$DevelopmentStage,2,median), col="black", lwd=3, lty = 2, type="l")

plot(cumsum(apply(C_cycle$alloc_StorageOrgan_gCm2day,2,median)), col="black", lwd=3, lty = 2, type="l")

for (i in seq(1, 300)) { if (i == 1) {plot(states_all$harvest_gCm2day[i,], ylim=c(0,250)) } else {lines(states_all$harvest_gCm2day[i,])} ; print(states_all$harvest_gCm2day[i,which(states_all$harvest_gCm2day[i,] > 0)]) }

#plot(apply(states_all$APAR_MJm2day,2,median)/(drivers$met[,4]*0.5), type="l", lwd=3, ylim=c(0,1))
#lines(apply(C_cycle$APAR_MJm2day,2,median)/(new_drivers$met[,4]*0.5), col="green", lwd=3)
##plotCI(drivers$obs[,23], uiw = drivers$obs[,24], add=TRUE)


