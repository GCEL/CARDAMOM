
###
## This script is intended to take CARDAMOM gridded analyses and 
## post-process them into temporally downscaled estimates (i.e. >= daily -> hourly).
## These estimates continue to have quantile based uncertainty information for each location and time step
## Version 1: T. L. Smallman (22/04/2021)
###

###
## Load R libraries and user defined functions into memory

# set CARDAMOM working directory
setwd("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/")
# load all functions
source("./R_functions/load_all_cardamom_functions.r")
# Load specific library for this script
library(ncdf4)

###
## Define input project information

# load CARDAMOM gridded project file
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_LU_FIRES_MHMCMC/global_1deg_C1/infofile.RData")
load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_MHMCMC/DAREUK_0.5deg_monthly/infofile.RData")
# Then load parameter and C stocks / flux files
load(paste(PROJECT$results_processedpath,"/",PROJECT$name,"_parameter_maps.RData",sep=""))
load(paste(PROJECT$results_processedpath,"/",PROJECT$name,"_stock_flux.RData",sep=""))
# Output file name modifiers
output_prefix = "" # must be followed by "_"
output_suffix = "" # must be preceeded by "_"

###
## Load weather_generator to provide downscaling information on radiation and temperature curves

## NOTE: that the weather generator source code is currently hardcoded to provide hourly information
## Therefore we will hardcode the number of steps_per_day expected for the final output
steps_per_day = 24 # assuming hourly
hrs_in_analysis = sum(PROJECT$model$timestep_days)*steps_per_day
simulation_years = c(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year))
nos_years = length(simulation_years)
steps_per_year = length(PROJECT$model$timestep_days) / nos_years
if (steps_per_year != floor(steps_per_year)) {
    stop("Need to update the code to account for years with different number of time steps")
}

# create the shared object from the weather generator code
system('R CMD SHLIB -o ./LIBRARY/weather_generator_F/weather_generator.so ./LIBRARY/weather_generator_F/weather_generator.f90 ./LIBRARY/weather_generator_F/wg_interface.f90')
# load the shared object into memory for use in downscaling meteorology
dyn.load("./LIBRARY/weather_generator_F/weather_generator.so")

### 
## Declare gridded output variables for those we want to use

# Specify the number of quantiles we will be extracting from the CARDAMOM output.
# We will assume we take the lowest, median and highest quantiles
quantiles_locs = median(c(1:dim(grid_output$nee_gCm2day)[2]))
#quantiles_locs = c(1,tmp,dim(grid_output$nee_gCm2day)[2])
#quantiles_wanted = grid_output$num_quantiles[quantiles_locs]
quantiles_wanted = c(0.5)
#quantiles_wanted = c(0.025,0.5,0.975)

###
## Loop through each year in turn

start = 0 ; finish = 0
for (yr in seq(1,nos_years)) {

     # How many days in the current year
     days_in_year = nos_days_in_year(simulation_years[yr])
     # Determine how many hours in the current year
     hrs_in_year = days_in_year * steps_per_day
     start = finish + 1 ; finish = yr * steps_per_year

     # Create the hourly output variables
     dims = dim(grid_output$mean_nee_gCm2day)
     nee_gCm2hr = array(NA, dim=c(dims[1],dims[2],length(quantiles_locs),hrs_in_year))
     nbe_gCm2hr = array(NA, dim=c(dims[1],dims[2],length(quantiles_locs),hrs_in_year))
     gpp_gCm2hr = array(NA, dim=c(dims[1],dims[2],length(quantiles_locs),hrs_in_year))
     reco_gCm2hr = array(NA, dim=c(dims[1],dims[2],length(quantiles_locs),hrs_in_year))
     fire_gCm2hr = array(NA, dim=c(dims[1],dims[2],length(quantiles_locs),hrs_in_year))
     harvest_gCm2hr = array(NA, dim=c(dims[1],dims[2],length(quantiles_locs),hrs_in_year))

     # Loop through each pixel in turn and downscale
     for (n in seq(1,PROJECT$nosites)) {

          # Update user
          if (n%%100) {print(paste("Current at site ",n," of ",PROJECT$nosites," in year: ",simulation_years[yr], sep=""))}

          # Load the drivers for the current pixel
          drivers = read_binary_file_format(paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[n],".bin",sep=""))

          # Extract the grid location for the current pixel
          # This just helps keep the coding tidy
          slot_j=grid_output$j_location[n]
          slot_i=grid_output$i_location[n]
          
          # Ensure that the model has a run for this location
          if (is.na(slot_i) == FALSE & is.na(slot_j) == FALSE) {

              # Extract location information
              latitude = grid_output$lat[slot_i,slot_j] # degree (-90/90)
              # Extract time information showing the day of year
              time = drivers$met[start:finish,6] # Julian day of year
              # generate but from different file structure
              sat_avg_in=drivers$met[start:finish,14]# C
              sat_min_in=drivers$met[start:finish,2] # C
              sat_max_in=drivers$met[start:finish,3] # C
              swrad_in=drivers$met[start:finish,4]*11.57407 # MJ/m2/day -> W/m2
              ## not used in this analysis
              coa_in = drivers$met[start:finish,5] # CO2 ppm
              ppt_in=drivers$met[start:finish,7] # kgH2O/m2/s
              wind_in=drivers$met[start:finish,15] # m/s
              rh_avg_in=rep(0.5,length(time)) # 0-1
              rh_max_in=rep(0.5,length(time)) # 0-1
              rh_min_in=rep(0.5,length(time)) # 0-1
              sf_pressure_in=rep(101300,length(time)) # Pa

              # The weather generator will downscale each "day" given to it but does not expect every day to be provided
              tmp=.Fortran("weathergeneratorinterface",latitude=as.single(rep(latitude,length.out=length(time))),nos_days=as.integer(length(time)),days_in_year=as.single(days_in_year)
	                                              ,time=as.single(time)
                                                      ,sat_avg_in=as.single(sat_avg_in),sat_max_in=as.single(sat_max_in),sat_min_in=as.single(sat_min_in)
                                                      ,ppt_in=as.single(ppt_in),swrad_in=as.single(swrad_in),coa_in=as.single(coa_in)
                                                      ,rh_avg_in=as.single(rh_avg_in),rh_max_in=as.single(rh_max_in),rh_min_in=as.single(rh_min_in),wind_in=as.single(wind_in)
                                                      ,sat_out=as.single(array(0,dim=c(length(time)*steps_per_day))),ppt_out=as.single(array(0,dim=c(length(time)*steps_per_day)))
                                                      ,swrad_out=as.single(array(0,dim=c(length(time)*steps_per_day))),coa_out=as.single(array(0,dim=c(length(time)*steps_per_day)))
                                                      ,rh_out=as.single(array(0,dim=c(length(time)*steps_per_day))),wind_out=as.single(array(0,dim=c(length(time)*steps_per_day))) )

              # Extract wanted output
              sat_out=tmp$sat_out ; swrad_out=tmp$swrad_out
              #coa_out=tmp$coa_out ; ppt_out=tmp$ppt_out 
              #rh_out=tmp$rh_out   ; wind_out=tmp$wind_out
              # Tidy away output
              rm(tmp) ; gc()

              ## Determine downscaling coefficient for GPP from daily to hourly based on radiation curve
              hold_gpp = swrad_out/rep(rollapply(swrad_out, width = steps_per_day, by = steps_per_day, FUN=sum, na.rm=TRUE), each = steps_per_day)
              ## Determine downscale respiration fluxes based on temperature curve and exponential response function
              hold_resp=exp(grid_parameters$parameters[slot_i,slot_j,10,quantiles_locs[2]]*sat_out)
              hold_resp=hold_resp/rep(rollapply(hold_resp, width = steps_per_day, by = steps_per_day, FUN=sum, na.rm=TRUE), each=steps_per_day)
              # Determine hourly flux rate for variables which we assume constant emission over the day
              hold_const = 1/steps_per_day

              # Now apply the scaler adjustments and expand to fill each day of the month with the same value
              timestep_days_local = PROJECT$model$timestep_days[start:finish]
              for (q in seq(1, length(quantiles_locs))) {
                   yy = 0 ; ww = 0 
                   for (zz in seq(1,length(time))) {
                        # Determine the period to be covered by the same day
                        yy = yy + steps_per_day ; ww = ww + (timestep_days_local[zz]*steps_per_day) ; adjustment=(ww-((timestep_days_local[zz]*steps_per_day)-1))
                        # Apply disaggregation based on radiation curve for gpp
                        gpp_gCm2hr[slot_i,slot_j,q,adjustment:ww] = grid_output$gpp_gCm2day[n,quantiles_locs[q],(start+zz-1)]*hold_gpp[(yy-(steps_per_day-1)):yy]
                        # Apply disaggredation based on temperature curve for ra
                        reco_gCm2hr[slot_i,slot_j,q,adjustment:ww] = grid_output$reco_gCm2day[n,quantiles_locs[q],(start+zz-1)]*hold_resp[(yy-(steps_per_day-1)):yy]
                        # Apply disaggredation based on constant emissions for fire
                        fire_gCm2hr[slot_i,slot_j,q,adjustment:ww] = grid_output$fire_gCm2day[n,quantiles_locs[q],(start+zz-1)]*hold_const
                        # Apply disaggredation based on constant emissions for fire
                        harvest_gCm2hr[slot_i,slot_j,q,adjustment:ww] = grid_output$harvest_gCm2day[n,quantiles_locs[q],(start+zz-1)]*hold_const
                   } # Downscaling and filling
              } # loop each quantile

          } # Is there an analysis for this location

     } # end location loop

     # At the moment we assume to derive the NEE and NBE based on their components.
     # NOTE: that due to the non-gaussian nature of many of our fluxes, and varying over time that these will not be mass balanced
     nee_gCm2hr = reco_gCm2hr - gpp_gCm2hr 
     nbe_gCm2hr = fire_gCm2hr + reco_gCm2hr - gpp_gCm2hr

     ###
     ## Output now

     # create lat / long axes, assumes WGS-84 grid
     latitude = seq(PROJECT$latitude[1]+(PROJECT$resolution*0.5),PROJECT$latitude[2]-(PROJECT$resolution*0.5), length.out = PROJECT$lat_dim)
     longitude = seq(PROJECT$longitude[1]+(PROJECT$resolution*0.5),PROJECT$longitude[2]-(PROJECT$resolution*0.5), length.out = PROJECT$long_dim)

     ## define dimension
     lat_dimen <- ncdim_def( "lat", units="degree north (-90->90)", latitude )
     long_dimen <- ncdim_def( "lon", units="degree east (-180->180)", longitude )
     time_dimen <- ncdim_def( "time", units="", 1:hrs_in_year)
     quantile_dimen <- ncdim_def( "quantile", units="-", quantiles_wanted)
 
     ## define output variable
     var0 = ncvar_def("Time", units = "hr", longname = paste("Hours since 01/01/",simulation_years[yr],sep=""), 
                      dim=list(time_dimen), missval = -99999, prec="double", compression = 9)

     ## FLUXES
     # GPP
     var1  = ncvar_def("gpp_ensemble",   unit="gC.m-2.h-1", longname = "Gross Primary Productivity - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
     # Ecosystem respiration
     var2  = ncvar_def("reco_ensemble",  unit="gC.m-2.h-1", longname = "Ecosystem (Ra + Rh) Respiration - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
     # Net Ecosystem Exchange
     var3  = ncvar_def("nee_ensemble",   unit="gC.m-2.h-1", longname = "Net Ecosystem Exchange - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
     # Net Biome Exchange
     var4  = ncvar_def("nbe_ensemble",   unit="gC.m-2.h-1", longname = "Net Biome Exchange (NEE + Fire) - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
     # Fire emissions
     var5  = ncvar_def("fFire_ensemble", unit="gC.m-2.h-1", longname = "Fire - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
     # Flux from forest loss
     var6 = ncvar_def("fLuc_ensemble",  unit="gC.m-2.h-1", longname = "Forest harvest - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)

     # create the empty file
     output_name = paste(PROJECT$results_processedpath,output_prefix,"CFLUX_hourly_downscaled_",simulation_years[yr],output_suffix,".nc",sep="")
     new_file=nc_create(filename=output_name, vars=list(var0,                                      
                                                        var1,var2,var3,var4,var5,var6),
                                                        force_v4 = TRUE)

## Load data into output variable

## TIMING
ncvar_put(new_file, var0, 1:hrs_in_year)

## FLUXES
# GPP
ncvar_put(new_file, var1, gpp_gCm2hr)
# RECO
ncvar_put(new_file, var2, reco_gCm2hr)
# NEE
ncvar_put(new_file, var3, nee_gCm2hr)
# NBE
ncvar_put(new_file, var4, nbe_gCm2hr)
# FIR
ncvar_put(new_file, var5, fire_gCm2hr)
# Forest Harvest
ncvar_put(new_file, var6, harvest_gCm2hr)

## close the file to write to disk
nc_close(new_file)

} # Year loop

# unload the object
dyn.unload("./LIBRARY/weather_generator_F/weather_generator.so")

