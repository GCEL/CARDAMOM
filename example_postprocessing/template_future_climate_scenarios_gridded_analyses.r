
###
## Script to run a single site under climate change scenarios
## The script has the following assumptions:
##   1) we are using the same executable as that used
##      in the calibration process. 
##   2) simulating to steady state, i.e. no disturbance scenarios
##   3) simulating to 2100 at a monthly time step
##   4) output full ensembles for later plotting
##   5) this script does not plot the outputs
###

# ensure clean R
rm(list=ls())

# Update the user
print("Starting steady_state")

###
## User options
### 

# Set your cardamom working directory
setwd("</your/CARDAMOM/directory/here/>")

# Set the location of cardamom outputs
load("./CARDAMOM_LOCAL/cardamompathdef.RData")
# Specify the infofile for the project to be used
infofile = paste(load_paths$cardamom_outputs,"/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/global_2_2.5deg_C7_GCP_AGB/infofile.RData",sep="")
# Load project for climate
load(infofile)

# Suffix for the output files
output_suffix = "_steadystate" # must include "_" at the beginning
# Location to place the outputs of this script
outdir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/global_2_2.5deg_C7_GCP_AGB/RESULTS_PROCESSED/"

# Which Shared Socioeconomic Pathways (SSPs) do you want to use?
# Currently available in GCEL are: "ssp119","ssp126","ssp434","ssp245","ssp370","ssp585"
ssp_scenarios = c("ssp119","ssp126","ssp434","ssp245","ssp370","ssp585")

# Which ESM to extract climate change from?
# Currently available are: "MOHC" 
ESM = "MOHC"
# How far into the future do we go? 
# NOTE: max 2100
new_end_year = 2100

# Load libraries and CARDAMOM R code
source("./R_functions/load_all_cardamom_functions.r")
library(RcppRoll)
library(abind)


# Local version (reduced complexity) of the run_each_site function to allow for parallel operation
run_each_site_local<-function(site,ssp,PROJECT,outdir,original_end,cardamom_lat,cardamom_long,
                              lat_esm,long_esm,begin_anomaly_esm,end_anomaly_esm,end_esm,
                              grid_rainfall_kgm2s,grid_swrad_MJm2day,grid_maxt_C,grid_mint_C,
                              grid_airt_C,grid_vpd_Pa) { #,grid_wind_spd_ms,grid_harvest)

  # Update the user
  print(paste("At site = ",site," of ",PROJECT$nosites," ",Sys.time(),sep=""))

  # Determine which parameter chains we will be using
  output = determine_parameter_chains_to_run(PROJECT,site)
  parameters = output$parameters ; converged = output$converged ; rm(output)

  # Define the output variables name                  
  outfile_stock_fluxes = paste(outdir,PROJECT$sites[site],"_",ssp_scenarios[ssp],"_",ESM,output_suffix,".RData",sep="")
              
  # Check whether this site likely has a valid set of parameters
  if (length(as.vector(parameters)) != 1 & file.exists(outfile_stock_fluxes) == FALSE) {

      # Read in the current sites driver information
      drivers = read_binary_file_format(paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[site],".bin",sep=""))
              
      # Find the latitude / longitude of the site to be ran
      slot_j = as.numeric(PROJECT$sites[site])/PROJECT$long_dim
      slot_i = as.numeric(PROJECT$sites[site])-(floor(slot_j)*PROJECT$long_dim)
      if(slot_i == 0) {slot_i = PROJECT$long_dim} ; slot_j = ceiling(slot_j)

      # Find site within the climate scenario analysis grid.
      # We will use this to filter the larger grid into something more efficient
      # NOTE: this assumes that each "site" (experiment) in the selected project is the same
      output = closest2d_3(1,lat_esm,long_esm,cardamom_lat[slot_i,slot_j],cardamom_long[slot_i,slot_j])
      i1 = unlist(output)[1] ; j1 = unlist(output)[2]

      # Extract only the location we want
      rainfall_kgm2s = grid_rainfall_kgm2s[i1,j1,]
      swrad_MJm2day  = grid_swrad_MJm2day[i1,j1,]
      #wind_spd_ms    = grid_wind_spd_ms[i1,j1,]
      maxt_C         = grid_maxt_C[i1,j1,]
      mint_C         = grid_mint_C[i1,j1,]
      airt_C         = grid_airt_C[i1,j1,]
      vpd_Pa         = grid_vpd_Pa[i1,j1,]
      # NOTE: co2_ppm is a vector not a grid so a common value is assumed across the globe

#      # Find site within the harvest scenario analysis grid
#      # NOTE: assume that all within the current PROJECT are the same location but different experiments
#      output = closest2d_3(1,lat_luh2,long_luh2,cardamom_lat[slot_i,slot_j],cardamom_long[slot_i,slot_j])
#      i1 = unlist(output)[1] ; j1 = unlist(output)[2]
#     
#      # Extract just the location we want
#      harvest = grid_harvest[i1,j1,]

      ##
      # Create site specific meteorological drivers

      if (exists("steps_per_year") == FALSE) {
          steps_per_year = dim(drivers$met)[1] / length(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year))
      }
    
      # row bind the entire met array to expand by adding the times series on each time
      new_drivers = drivers
      for (r in seq(1,nos_loops)) {
           new_drivers$met = rbind(new_drivers$met,drivers$met)
      }
      
      # Just looping the existing drivers doesn't work. 
      # At the very least the simulation day need to be updated so continuously increase.
      # We can achieve this by cumulative sum of the day of year variable
      new_drivers$met[,1] = PROJECT$model$timestep_days ; new_drivers$met[,1] = cumsum(new_drivers$met[,1])
      # Remove biomass removals (i.e. management) - but leave fire
      new_drivers$met[(dim(drivers$met)[1]+1):dim(new_drivers$met)[1],8] = 0

#      # Determine the start point in the management scenario, i.e. the year after the observed period
#      begin_anomaly = which(years_since_2015 == (max(0,(original_end-2015))+1))
#      harvest = harvest[begin_anomaly:length(years_since_2015)]
#      rm(begin_anomaly)
#      # Impose LUH2 biomass removal scenarios
#      start = (dim(drivers$met)[1]+1) ; finish = dim(new_drivers$met)[1]
#      new_drivers$met[start:finish,8] = rep(harvest, each = steps_per_year)[1:length(c(start:finish))]
#      new_drivers$met[is.na(new_drivers$met[,8]) == TRUE,8] = 0

      # Adjust to make input array match the maximum length of the anomaly
      new_drivers$met = new_drivers$met[1:(dim(drivers$met)[1]+length(co2_ppm)),]
      # Determine the number of years in the simulation in total
      nos_years = dim(new_drivers$met)[1] / steps_per_year
      # Determine the number of years in the anomaly
      nos_years_anomaly = nos_years - length(c(as.numeric(PROJECT$start_year):as.numeric(original_end)))
      # Number years needed for output from dalec.so, make sure that the original PROJECT time period is updated
      analysis_years = length(c(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)))
      if (analysis_years != nos_years) {
          PROJECT$end_year = as.numeric(PROJECT$end_year) + (nos_years - analysis_years)
      }

      # Determine the start point for introducing the ESM time series
      start = dim(drivers$met)[1]+1 ; finish = dim(new_drivers$met)[1]
      # Determine the overlap period from the original time series. 
      # We assume here that both time series are in monthly step.
      f_final_yr = dim(drivers$met)[1] ; s_final_yr = f_final_yr-(length(begin_anomaly_esm:end_anomaly_esm)-1)
      # CO2 is a special case as we used a flat time series previously and the below approach attempt to bring the ESM time series towards
      # our existing time series.
      co2_bias = co2_ppm + new_drivers$met[f_final_yr,5]
      # Determine the bias between each month in the overlap period
      mint_bias     = new_drivers$met[s_final_yr:f_final_yr,2]  - mint_C[begin_anomaly_esm:end_anomaly_esm]
      maxt_bias     = new_drivers$met[s_final_yr:f_final_yr,3]  - maxt_C[begin_anomaly_esm:end_anomaly_esm]
      swrad_bias    = new_drivers$met[s_final_yr:f_final_yr,4]  - swrad_MJm2day[begin_anomaly_esm:end_anomaly_esm]
      rainfall_bias = new_drivers$met[s_final_yr:f_final_yr,7]  - rainfall_kgm2s[begin_anomaly_esm:end_anomaly_esm]
      rollmaxt_bias = new_drivers$met[s_final_yr:f_final_yr,10] - maxt_C[begin_anomaly_esm:end_anomaly_esm]
      rollvpd_bias  = new_drivers$met[s_final_yr:f_final_yr,12] - vpd_Pa[begin_anomaly_esm:end_anomaly_esm]
      airt_bias     = new_drivers$met[s_final_yr:f_final_yr,14] - airt_C[begin_anomaly_esm:end_anomaly_esm]
#      wind_spd_bias = new_drivers$met[s_final_yr:f_final_yr,15] - wind_spd_ms[begin_anomaly_esm:end_anomaly_esm]
      vpd_bias      = new_drivers$met[s_final_yr:f_final_yr,16] - vpd_Pa[begin_anomaly_esm:end_anomaly_esm]
      # Loop through each month and determine the mean bias for each month
      tmp = rep(NA, steps_per_year)
      # Minimum temperature
      for (step in seq(1, steps_per_year)) {tmp[step] = mean(mint_bias[seq(step,length(mint_bias),steps_per_year)])} # month loop
      mint_bias = tmp
      # Maximum temperature
      for (step in seq(1, steps_per_year)) {tmp[step] = mean(maxt_bias[seq(step,length(maxt_bias),steps_per_year)])} # month loop
      maxt_bias = tmp
      # SW radiation
      for (step in seq(1, steps_per_year)) {tmp[step] = mean(swrad_bias[seq(step,length(swrad_bias),steps_per_year)])} # month loop
      swrad_bias = tmp
      # Rainfall
      for (step in seq(1, steps_per_year)) {tmp[step] = mean(rainfall_bias[seq(step,length(rainfall_bias),steps_per_year)])} # month loop
      rainfall_bias = tmp
      # Rolling max temperature
      for (step in seq(1, steps_per_year)) {tmp[step] = mean(rollmaxt_bias[seq(step,length(rollmaxt_bias),steps_per_year)])} # month loop
      rollmaxt_bias = tmp
      # Rolling VPD
      for (step in seq(1, steps_per_year)) {tmp[step] = mean(rollvpd_bias[seq(step,length(rollvpd_bias),steps_per_year)])} # month loop
      rollvpd_bias = tmp
      # Mean temperature
      for (step in seq(1, steps_per_year)) {tmp[step] = mean(airt_bias[seq(step,length(airt_bias),steps_per_year)])} # month loop
      airt_bias = tmp
      ## Wind speed
      #for (step in seq(1, steps_per_year)) {tmp[step] = mean(wind_spd_bias[seq(step,length(wind_spd_bias),steps_per_year)])} # month loop
      #wind_spd_bias = tmp
      # VPD
      for (step in seq(1, steps_per_year)) {tmp[step] = mean(vpd_bias[seq(step,length(vpd_bias),steps_per_year)])} # month loop
      vpd_bias = tmp
      # Tidy
      rm(step,tmp)
      # Apply mean bias to the ESM drivers
      mint_bias = mint_bias + mint_C
      maxt_bias = maxt_bias + maxt_C
      swrad_bias = swrad_bias + swrad_MJm2day
      rainfall_bias = rainfall_bias + rainfall_kgm2s
      rollmaxt_bias = rollmaxt_bias + maxt_C
      rollvpd_bias = rollvpd_bias + vpd_Pa               
      airt_bias = airt_bias + airt_C
#      wind_spd_bias = wind_spd_bias + wind_spd_ms
      vpd_bias = vpd_bias + vpd_Pa
      # Remove the first year, i.e. the overlap year
      mint_bias = mint_bias[-c(begin_anomaly_esm:end_anomaly_esm)]
      maxt_bias = maxt_bias[-c(begin_anomaly_esm:end_anomaly_esm)]
      swrad_bias = swrad_bias[-c(begin_anomaly_esm:end_anomaly_esm)]
      rainfall_bias = rainfall_bias[-c(begin_anomaly_esm:end_anomaly_esm)]
      rollmaxt_bias = rollmaxt_bias[-c(begin_anomaly_esm:end_anomaly_esm)]
      rollvpd_bias = rollvpd_bias[-c(begin_anomaly_esm:end_anomaly_esm)]
      airt_bias = airt_bias[-c(begin_anomaly_esm:end_anomaly_esm)]
      #wind_spd_bias = wind_spd_bias[-c(begin_anomaly_esm:end_anomaly_esm)]
      vpd_bias = vpd_bias[-c(begin_anomaly_esm:end_anomaly_esm)]   
      # Replace the final part of the time series with the bias corrected ESM data
      new_drivers$met[start:finish,2]  = mint_bias
      new_drivers$met[start:finish,3]  = maxt_bias
      new_drivers$met[start:finish,4]  = swrad_bias
      new_drivers$met[start:finish,5]  = co2_bias
      new_drivers$met[start:finish,7]  = rainfall_bias
      new_drivers$met[start:finish,10] = rollmaxt_bias
      new_drivers$met[start:finish,12] = rollvpd_bias
      new_drivers$met[start:finish,14] = airt_bias
#      new_drivers$met[start:finish,15] = wind_spd_bias
      new_drivers$met[start:finish,16] = vpd_bias
      # Ensure realism is obeyed
      new_drivers$met[new_drivers$met[,4] < 0,4] = 0   # Can't have less than zero radiation
      new_drivers$met[new_drivers$met[,7] < 0,7] = 0   # Can't have less than zero rainfall
      new_drivers$met[new_drivers$met[,12] < 0,12] = 0 # Can't have less than zero 21-day mean VPD
      new_drivers$met[new_drivers$met[,16] < 0,16] = 0 # Can't have less than zero VPD

      # Run the model
      # NOTE: we assume that we are using site one as an example of the correct parameter number etc.
      soil_info = c(new_drivers$top_sand,new_drivers$bot_sand,new_drivers$top_clay,new_drivers$bot_clay)
      states_all = simulate_all(1,PROJECT,PROJECT$model$name,new_drivers$met,
                                parameters[1:PROJECT$model$nopars[1],,],
                                new_drivers$lat,PROJECT$ctessel_pft[1],PROJECT$parameter_type,
                                PROJECT$exepath,soil_info)

      # Post-process the DALEC model output for both site and gridded analyses
      states_all = post_process_dalec(states_all,parameters,drivers,PROJECT,site)
                                      
      # we want straight forward reduced dataset of common stocks and fluxes
      #num_quantiles = c(0.025,0.05,0.25,0.5,0.75,0.95,0.975)
      num_quantiles = c(0.025,0.5,0.975)
      na_flag = TRUE

      # Run post-processing for gridded analysis
      dummy = post_process_for_grid(outfile_stock_fluxes,PROJECT,drivers,parameters,num_quantiles,na_flag,converged,states_all)
      if (dummy == 0) { return(outfile_stock_fluxes) }

      ## Return the site_output list object for subsequent slotting into the wider grid
      #if (dummy != 0) { 
      #    print(paste("Error post-processing site = ",PROJECT$site[site])) 
      #} else {
      #    site_output_all = append(site_output_all,outfile_stock_fluxes)
      #}

     } # valid set of parameters?
     
} # end function run_each_site_local

###
## SSP information which should not need changing
###

# Path to CMIP6 files
cmip6dir = "/exports/csce/datastore/geos/groups/gcel/cmip6/TLS/ScenarioMIP/"
## Path to land use scenarios (LUH2)
#luh2dir = "/exports/csce/datastore/geos/groups/gcel/LUH2/"

###
## Load needed libraries and functions
###

## Create any user defined functions
time_sum = function(var, interval) {
   return(rollapply(var, sum, by = interval, width = interval, na.rm=TRUE))
}
time_mean = function(var, interval) {
   return(rollapply(var, mean, by = interval, width = interval, na.rm=TRUE))
}

###
## Determine required timing information for the analysis
###

# We will go beyond the currently calibrated time period. 
# Store the original for later resets
original_end = as.numeric(PROJECT$end_year)

# How many times do we need to loop our meteorology to reach the target year?
# -1 adjusts for the current meteorological dataset
nos_loops = ceiling(((new_end_year-as.numeric(PROJECT$start_year)) 
                     / length(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year))) - 1)

# Calculate the number of years
PROJECT$nos_years = (as.numeric(PROJECT$end_year) - as.numeric(PROJECT$start_year))+1
  
###
## Determine the required spatial information
###

# generate the lat / long grid again
output = generate_wgs84_grid(PROJECT$latitude,PROJECT$longitude,PROJECT$resolution)
cardamom_lat = array(output$lat, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
cardamom_long = array(output$long,dim=c(PROJECT$long_dim,PROJECT$lat_dim))

###
## Any flags needed to allow for the borrowing of CARDAMOM R functions into this analysis script
use_parallel = TRUE
numWorkers = 48
repair = 1
          
###
## Begin projection analysis
###

# Loop through the desired SSPs
for (ssp in seq(1, length(ssp_scenarios))) {

     ###
     ## Load climate + CO2 information for the current SSP
     
     # List all available scenario filesFind all the anomalies files for the current scenario / ESM climate anomaly information
     list_of_files = list.files(paste(cmip6dir,ESM,"/",ssp_scenarios[ssp],"/month",sep=""), full.names=TRUE, recursive = TRUE)
     if (new_end_year > 2049) { 
         list_of_files = list_of_files[grepl("201501",list_of_files) == TRUE | grepl("205001",list_of_files) == TRUE]
         expected_nos_files = 2
     } else {
         list_of_files = list_of_files[grepl("201501",list_of_files) == TRUE]
         expected_nos_files = 1
     } 

     # Check that files exist for each of the variables we want
     if (length(list_of_files[grepl("/pr/",list_of_files)]) == expected_nos_files & 
         length(list_of_files[grepl("/rsds/",list_of_files)]) == expected_nos_files & 
         #length(list_of_files[grepl("/sfcWind/",list_of_files)]) == expected_nos_files & 
         length(list_of_files[grepl("/tasmax/",list_of_files)]) == expected_nos_files &
         length(list_of_files[grepl("/tasmin/",list_of_files)]) == expected_nos_files & 
         length(list_of_files[grepl("/tas/",list_of_files)]) == expected_nos_files & 
         length(list_of_files[grepl("/vpd/",list_of_files)]) == expected_nos_files &
         length(list_of_files[grepl("/co2mass/",list_of_files)]) == expected_nos_files ) {

         # Open the climate variables we want
         sub_list_of_files = list_of_files[grepl("201501",list_of_files) == TRUE]
         rainfall_anomaly_file = nc_open(sub_list_of_files[grepl("/pr/",sub_list_of_files)])
         swrad_anomaly_file    = nc_open(sub_list_of_files[grepl("/rsds/",sub_list_of_files)])
         #wind_spd_anomaly_file = nc_open(sub_list_of_files[grepl("/sfcWind/",sub_list_of_files)])
         maxt_anomaly_file     = nc_open(sub_list_of_files[grepl("/tasmax/",sub_list_of_files)]) 
         mint_anomaly_file     = nc_open(sub_list_of_files[grepl("/tasmin/",sub_list_of_files)]) 
         airt_anomaly_file     = nc_open(sub_list_of_files[grepl("/tas/",sub_list_of_files)])    
         vpd_anomaly_file      = nc_open(sub_list_of_files[grepl("/vpd/",sub_list_of_files)])
         co2_anomaly_file      = nc_open(sub_list_of_files[grepl("/co2mass/",sub_list_of_files)])
         # Now read in each variable, applying relevant unit conversions
         # e.g. rainfall_anomaly_file$var$pr$units
         grid_rainfall_kgm2s = ncvar_get(rainfall_anomaly_file,"pr")
         grid_swrad_MJm2day  = ncvar_get(swrad_anomaly_file,"rsds") * 86400 * 1e-6 # W/m2 -> MJ/m2/day
         #grid_wind_spd_ms    = ncvar_get(wind_spd_anomaly_file,"sfcWind") 
         grid_maxt_C         = ncvar_get(maxt_anomaly_file,"tasmax") - 273.15 # K->C
         grid_mint_C         = ncvar_get(mint_anomaly_file,"tasmin") - 273.15 # K->C
         grid_airt_C         = ncvar_get(airt_anomaly_file,"tas") - 273.15
         grid_vpd_Pa         = ncvar_get(vpd_anomaly_file,"vpd")
         co2_ppm        = ncvar_get(co2_anomaly_file,"co2mass")/1e+12/7.804816 # kg -> ppm 
         # Where:
         # 1 ppmv of CO2= 2.13 Gt of C
         # 44.01 CO2 mass
         # 1 ppm CO2 = 2.13/12.0107*44.01 Gt = 7.804816
         # convert to Gt
   
         # Read in lat / long and time information from one of the files - it should be the same in all cases
         days_since_1850 = ncvar_get(vpd_anomaly_file, "time")
         lat_esm = ncvar_get(vpd_anomaly_file, "lat") # degrees north (-89.5 / 89.5)
         long_esm = ncvar_get(vpd_anomaly_file, "lon") # degrees east (-179.5 / 179.5)
         # Come out 0-360, convert to correct
         long_esm[which(long_esm > 180)] = long_esm[which(long_esm > 180)]-360 
    
         # Tidy up currently open files
         nc_close(rainfall_anomaly_file) ; nc_close(swrad_anomaly_file) #; nc_close(wind_spd_anomaly_file)
         nc_close(maxt_anomaly_file) ; nc_close(mint_anomaly_file) ; nc_close(airt_anomaly_file)
         nc_close(vpd_anomaly_file) ; nc_close(co2_anomaly_file) 

         # If we need to second half of the time series then read them in here
         if (expected_nos_files == 2) {
             # Open the climate variables we want
             sub_list_of_files = list_of_files[grepl("20500",list_of_files)==TRUE]
             rainfall_anomaly_file = nc_open(sub_list_of_files[grepl("/pr/",sub_list_of_files)])
             swrad_anomaly_file    = nc_open(sub_list_of_files[grepl("/rsds/",sub_list_of_files)])
             #wind_spd_anomaly_file = nc_open(sub_list_of_files[grepl("/sfcWind/",sub_list_of_files)])
             maxt_anomaly_file     = nc_open(sub_list_of_files[grepl("/tasmax/",sub_list_of_files)]) 
             mint_anomaly_file     = nc_open(sub_list_of_files[grepl("/tasmin/",sub_list_of_files)]) 
             airt_anomaly_file     = nc_open(sub_list_of_files[grepl("/tas/",sub_list_of_files)])    
             vpd_anomaly_file      = nc_open(sub_list_of_files[grepl("/vpd/",sub_list_of_files)])
             co2_anomaly_file      = nc_open(sub_list_of_files[grepl("/co2mass/",sub_list_of_files)])
             # Now read in each variable, applying relevant unit conversions
             # e.g. rainfall_anomaly_file$var$pr$units
             tmp_rainfall_kgm2s = ncvar_get(rainfall_anomaly_file,"pr")
             tmp_swrad_MJm2day  = ncvar_get(swrad_anomaly_file,"rsds") * 86400 * 1e-6 # W/m2 -> MJ/m2/day
             #tmp_wind_spd_ms    = ncvar_get(wind_spd_anomaly_file,"sfcWind")
             tmp_maxt_C         = ncvar_get(maxt_anomaly_file,"tasmax") - 273.15 # K->C
             tmp_mint_C         = ncvar_get(mint_anomaly_file,"tasmin") - 273.15 # K->C
             tmp_airt_C         = ncvar_get(airt_anomaly_file,"tas") - 273.15
             tmp_vpd_Pa         = ncvar_get(vpd_anomaly_file,"vpd")

             # Read time information and directly append
             days_since_1850 = append(days_since_1850,ncvar_get(vpd_anomaly_file, "time"))
             # Read CO2 information and directly append
             # NOTE: co2_ppm is a vector not a grid so a common value is assumed across the globe
             co2_ppm            = append(co2_ppm,ncvar_get(co2_anomaly_file,"co2mass")/1e+12/7.804816) # kg -> ppm                
             # Where:
             # 1 ppmv of CO2= 2.13 Gt of C
             # 44.01 CO2 mass
             # 1 ppm CO2 = 2.13/12.0107*44.01 Gt = 7.804816
             # convert to Gt

             # Append each dataset from the second file to the first - creating our complete time series
             grid_rainfall_kgm2s = abind(grid_rainfall_kgm2s,tmp_rainfall_kgm2s, along=3)
             grid_swrad_MJm2day  = abind(grid_swrad_MJm2day,tmp_swrad_MJm2day, along=3)
             #grid_wind_spd_ms    = abind(grid_wind_spd_ms,tmp_wind_spd_ms, along=3)
             grid_maxt_C         = abind(grid_maxt_C,tmp_maxt_C, along=3)
             grid_mint_C         = abind(grid_mint_C,tmp_mint_C, along=3)
             grid_airt_C         = abind(grid_airt_C,tmp_airt_C, along=3)
             grid_vpd_Pa         = abind(grid_vpd_Pa,tmp_vpd_Pa, along=3)

             # Tidy up the tmp_* variables
             rm(tmp_rainfall_kgm2s,tmp_swrad_MJm2day,tmp_maxt_C,tmp_mint_C,tmp_airt_C,tmp_vpd_Pa)#,tmp_wind_spd_ms)

             # Tidy up currently open files
             nc_close(rainfall_anomaly_file) ; nc_close(swrad_anomaly_file) #; nc_close(wind_spd_anomaly_file)
             nc_close(maxt_anomaly_file) ; nc_close(mint_anomaly_file) ; nc_close(airt_anomaly_file)
             nc_close(vpd_anomaly_file) ; nc_close(co2_anomaly_file)  
 
         } # open second file if needed

#         # Open harvest scenario files
#         harvest_file = nc_open(paste(luh2dir,"/",ssp_scenarios[ssp],"/LUH2_biomass_removal_scenario.nc",sep=""))
#         grid_harvest = ncvar_get(harvest_file, "harvest_fraction")
#         lat_luh2 = ncvar_get(harvest_file, "lat") ; long_luh2 = ncvar_get(harvest_file, "lon")
#         years_since_2015 = ncvar_get(harvest_file, "time")
#
#         # Tidy up files
#         nc_close(harvest_file)
#
#         # Estimate division needed to the annual fraction of biomass harvest to that needed to be applied for each timestep
#         grid_harvest = grid_harvest / (365.25 / mean(PROJECT$model$timestep_days))

         # Input is monthly but we want to know this information for estimating the anomaly points
         # NOTE: these are only applicable for UKESM models which assume a 360 day year
         months_per_year = 12 ; days_per_month = 30 ; days_per_year = 360 
      
         # What time period is the of the climate scenario is the anomaly from?
         if (original_end < 2015) {stop('There is no overlap between the ESM and current analysis time periods - anomaly calculations are not possible')}
         # We will calculate the mean anomaly between the current and the ESM meteorology from the beginning of the ESM period
         begin_anomaly_esm = 1
         # and end it at the end of the current climate period
         end_anomaly_esm = max(0,(original_end - 2015))+1
         end_anomaly_esm = (end_anomaly_esm * months_per_year) 
         end_esm = dim(grid_rainfall_kgm2s)[3]
 
         # We will not keep the whole time series to that we can deterine a pixel / site specific
         # mean residual to apply a bias correction to the ESM time series
                                   
         # CO2 is a global average value only - i.e. no grid, based on global obserations so is consistent with the 
         # observation mean used in the calibration period here. So we will use a direct residual approach.
         # However, to achieve a smooth increase we bias correct from the final value in overlap period.
         # This approach contrasts the per month bias used in the rest of the meteorology.
         co2_ppm = co2_ppm - co2_ppm[months_per_year]
         # Select the future time period only
         co2_ppm = co2_ppm[(end_anomaly_esm+1):end_esm]

     } # have all of 2015-2049 files

     ###
     ## Loop through sites for the current climate scenario

     # NOTE: that the use of mclapply() is due to reported improved efficiency over creating a virtual cluster.
     # However, mclapply does not (at the time of typing) work on Windows, i.e. Linux and Mac only

     if (use_parallel) {
     
         ## Parallel case
         
         # Now we can deploy in anger
         cl <- min(PROJECT$nosites,numWorkers)
         site_output_all = mclapply(c(1:PROJECT$nosites),FUN=run_each_site_local,ssp=ssp,
                                    PROJECT=PROJECT,outdir=outdir,original_end=original_end,cardamom_lat=cardamom_lat,cardamom_long=cardamom_long,
                                    lat_esm=lat_esm,long_esm=long_esm,begin_anomaly_esm=begin_anomaly_esm,end_anomaly_esm=end_anomaly_esm,
                                    end_esm=end_esm,
                                    grid_rainfall_kgm2s,grid_swrad_MJm2day=grid_swrad_MJm2day,grid_maxt_C=grid_maxt_C,grid_mint_C=grid_mint_C,
                                    grid_airt_C=grid_airt_C,grid_vpd_Pa=grid_vpd_Pa,
                                    #grid_wind_spd_ms=grid_wind_spd_ms,grid_harvest=grid_harvest
                                    mc.cores = cl,mc.preschedule = FALSE)                                 
     } else {
     
         ## Serial case
     
         site_output_all = lapply(c(1:PROJECT$nosites),FUN=run_each_site_local,ssp=ssp,
                                  PROJECT=PROJECT,outdir=outdir,original_end=original_end,cardamom_lat=cardamom_lat,cardamom_long=cardamom_long,
                                  lat_esm=lat_esm,long_esm=long_esm,begin_anomaly_esm=begin_anomaly_esm,end_anomaly_esm=end_anomaly_esm,
                                  end_esm=end_esm,
                                  grid_rainfall_kgm2s,grid_swrad_MJm2day=grid_swrad_MJm2day,grid_maxt_C=grid_maxt_C,grid_mint_C=grid_mint_C,
                                  grid_airt_C=grid_airt_C,grid_vpd_Pa=grid_vpd_Pa)
                                  #grid_wind_spd_ms=grid_wind_spd_ms,grid_harvest=grid_harvest)
     }                        
     # Now assuming this has worked correctly, we can create the grid_output
     # First by loading an example of the site specific file
     n = 0 ; site_output = -1
     while (class(site_output) != "character" & n < length(site_output_all)) {
            n = n + 1
            site_output = site_output_all[[n]]
     }
     load(site_output)

     # Temporally update some timing information to ensure the grid_output object is correct
     PROJECT$end_year = as.character(new_end_year)
     PROJECT$nos_years = length(c(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)))

     # Define the intended output grid file name
     outfile_grid = paste(outdir,PROJECT$name,"_",ssp_scenarios[ssp],"_",ESM,output_suffix,".RData",sep="")
     # Create the list object
     grid_output = define_grid_output(PROJECT,repair,outfile_grid,site_output)
     # tidy up
     rm(site_output)

     # Combine into single common grid for further analysis
     grid_output = post_process_into_grid(grid_output,site_output_all,PROJECT)
      
     # update the user
     print("...writing combined grid_output to file")

     # now save the combined grid file
     save(grid_output, file=outfile_grid, compress = "gzip", compression_level = 9)

     # Tidy up...
     rm(grid_output) ; gc(reset=TRUE)
     # ...and reset
     PROJECT$end_year = as.character(original_end)
     PROJECT$nos_years = length(c(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)))
     
} # loop ssps

