#########################################################################################
# CARbon DAta MOdel fraMework (CARDAMOM) and DALEC terrestrial ecosystem model suite
# CARDAMOM is a Bayesian model-data fusion software framework. CARDAMOM is used to 
# assimilate observations and ecological theory to retrieve parameters for the 
# DALEC suite of intermediate complexity terrestrial ecosystem models. DALEC can be
# used as a fully integrated component of CARDAMOM or independently. 
# Copyright (C) 2024  University of Edinburgh,
#                     Mathew Williams (mat.williams@ed.ac.uk), 
#                     T. Luke Smallman (t.l.smallman@ed.ac.uk)
# UoE = University of Edinburgh

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# ########## File specific description ##########
# Function to run CARDAMOM parameters via the chosen model
# 
# Author: T. Luke Smallman (12/11/2024)
# Exceptions states below in specific functions
#
#########################################################################################

run_each_site<-function(n,PROJECT,repair,grid_override) {

  # Update the user 
  if (use_parallel == FALSE) {print(paste("Site = ",PROJECT$sites[n]," ",n," of ",PROJECT$nosites," ",Sys.time(),sep=""))}

  # Define the output file names
  outfile_site         = paste(PROJECT$results_processedpath,PROJECT$sites[n],".RData",sep="")
  outfile_parameters   = paste(PROJECT$results_processedpath,PROJECT$sites[n],"_parameters.RData",sep="")
  outfile_stock_fluxes = paste(PROJECT$results_processedpath,PROJECT$sites[n],"_stock_fluxes.RData",sep="")

  # Set dummy output variable, the value may be changes by the code below
  dummy = 0

  if (repair == 1 | (file.exists(outfile_stock_fluxes) == FALSE & file.exists(outfile_parameters) == FALSE)) {

      # Determine which parameter chains we will be using
      output = determine_parameter_chains_to_run(PROJECT,n)
      # Check if we likely have an error flag
      if (length(as.vector(output)) == 1) {
          # Return the flag and allow the job to move on
          return(output)
      }
      # Otherwise we should assume these variables exist
      parameters = output$parameters ; converged = output$converged ; kept_chains = output$kept_chains 
      # Then tidy
      rm(output)      

      # load the met data for each site
      drivers = read_binary_file_format(paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[n],".bin",sep=""))
## HACK to remove CO2 effect
#drivers$met[,5] = drivers$met[1,5]
# HACK to create S2 simulations for GCP / Trendy v15
#drivers$met[,8] = 0
# Irrigation of crops experiment, a minimum of at least 10 mm/day
#drivers$met[,7] = pmax(drivers$met[,7],100/86400)
#drivers$met[,7] = drivers$met[,7]*0.1
      # run parameters for full results / propogation
      soil_info = c(drivers$top_sand,drivers$bot_sand,drivers$top_clay,drivers$bot_clay)
      if (use_parallel == FALSE) {print("running model ensemble")}
      states_all = simulate_all(n,PROJECT,PROJECT$model$name,drivers$met,parameters[1:PROJECT$model$nopars[n],,],
                                drivers$lat,PROJECT$ctessel_pft[n],PROJECT$parameter_type,
                                PROJECT$exepath,soil_info)
## EDC hack to reduce the possible range of soil C changes
#filter = which(abs((states_all$som_gCm2[,dim(states_all$som_gCm2)[2]] - states_all$som_gCm2[,1]) / PROJECT$nos_years) < 250)
#if (length(filter) > 0) {
#parameters = array(parameters, dim=c((PROJECT$model$nopars[n]+1),prod(dim(parameters)[2:3])))[,filter]
#parameters = array(parameters, dim=c(dim(parameters)[1],dim(parameters)[2],2))
#states_all = simulate_all(n,PROJECT,PROJECT$model$name,drivers$met,parameters[1:PROJECT$model$nopars[n],,],
#                          drivers$lat,PROJECT$ctessel_pft[n],PROJECT$parameter_type,
#                          PROJECT$exepath,soil_info)
#} else {
#return(-9999)
#}
      if (use_parallel == FALSE) {print("model ensemble ran, on to post-processing")}

      # Avoid running with ACM basically where not all fluxes exist
      if (grepl("DALEC",PROJECT$model$name)) {

          ###
          # Derive stocks and fluxes used in the calculation of gridded aggregates
          # These are variables which for a site analysis would be easy to calculate
          # from the ensembles but difficult if not determined here and now before aggregation
          ###

          # Post-process the DALEC model output for both site and gridded analyses
          states_all = post_process_dalec(states_all,parameters,drivers,PROJECT,n)   
          # Determine how many ensemble members are within the observational uncertainties
          # of the calibration datasets
          states_all = assess_ensemble_fit_to_calibration_data(states_all,parameters,drivers,PROJECT)
      
      } # DALEC model or not?

      check_list = names(states_all)
      # pass to local variable for saving
      site_ctessel_pft = PROJECT$ctessel_pft[n]
      NPP_fraction = list(NPP_foliage_fraction = states_all$NPP_foliage_fraction)
      if (any(check_list == "NPP_roots_fraction") == TRUE) {NPP_fraction$NPP_roots_fraction = states_all$NPP_roots_fraction}
      if (any(check_list == "NPP_wood_fraction") == TRUE) {NPP_fraction$NPP_wood_fraction = states_all$NPP_wood_fraction}
      MTT_years = list(MTT_foliage_years = states_all$MTT_foliage_years)
      if (any(check_list == "MTT_labile_years") == TRUE) {MTT_years$MTT_labile_years = states_all$MTT_labile_years}
      if (any(check_list == "MTT_roots_years") == TRUE) {MTT_years$MTT_roots_years = states_all$MTT_roots_years}
      if (any(check_list == "MTT_wood_years") == TRUE) {MTT_years$MTT_wood_years = states_all$MTT_wood_years}
      if (any(check_list == "MTT_litter_years") == TRUE) {MTT_years$MTT_litter_years = states_all$MTT_litter_years}
      if (any(check_list == "MTT_woodlitter_years") == TRUE) {MTT_years$MTT_woodlitter_years = states_all$MTT_woodlitter_years}
      if (any(check_list == "MTT_som_years") == TRUE) {MTT_years$MTT_som_years = states_all$MTT_som_years}
      SS_gCm2 = list(SS_foliage_gCm2 = states_all$SS_foliage_gCm2)
      if (any(check_list == "SS_labile_gCm2") == TRUE) {SS_gCm2$SS_labile_gCm2 = states_all$SS_labile_gCm2}
      if (any(check_list == "SS_roots_gCm2") == TRUE) {SS_gCm2$SS_roots_gCm2 = states_all$SS_roots_gCm2}
      if (any(check_list == "SS_wood_gCm2") == TRUE) {SS_gCm2$SS_wood_gCm2 = states_all$SS_wood_gCm2}
      if (any(check_list == "SS_litter_gCm2") == TRUE) {SS_gCm2$SS_litter_gCm2 = states_all$SS_litter_gCm2}
      if (any(check_list == "SS_woodlitter_gCm2") == TRUE) {SS_gCm2$SS_woodlitter_gCm2 = states_all$SS_woodlitter_gCm2}
      if (any(check_list == "SS_som_gCm2") == TRUE) {SS_gCm2$SS_som_gCm2 = states_all$SS_som_gCm2}

      # Sanity check
      if (any(is.na(as.vector(NPP_fraction)))) {
      #if (length(which(is.na(as.vector(NPP_fraction))) == TRUE) > 0) {
          print(paste("NA value found in NPP for site ",PROJECT$site[n],sep="")) ; dummy = -4 ; return(dummy)
      }
      #if (use_parallel == FALSE) {print("processing and storing ensemble output")}

      # determine whether this is a gridded run (or one with the override in place)
      if (PROJECT$spatial_type == "site" | grid_override == TRUE) {

          # ...if this is a site run save the full ensemble and everything else...
          save(kept_chains,parameters,drivers,states_all,site_ctessel_pft,file=outfile_site, compress="gzip", compression_level = 6)
          # store the parameters and driver information
          save(kept_chains,parameters,drivers,site_ctessel_pft,NPP_fraction,MTT_years,SS_gCm2,#converged,
               file=outfile_parameters, compress="gzip", compression_level = 6)
          # Return
          dummy = 0 ; return(dummy)

      } else {

          # ...otherwise this is a grid and we want straight forward reduced dataset of common stocks and fluxes
          # The current quantiles allow for calculation of the 95 % CI (0.975-0.025) and the standard deviation equivalent (0.839-0.1607)
          # The other quantiles provide a equal description of the distribution.
          #num_quantiles = c(0.025,0.05,0.16,0.5,0.84,0.95,0.975) #; num_quantiles_agg = seq(0.0,1, length = 100)
          #num_quantiles = c(0.025,0.1607143,0.2964286,0.4321429,0.5,0.5678571,0.7035714,0.8392857,0.975)
          num_quantiles <- c(0,       # ensemble minimum
                             0.010,   # lower tail interpolation
                             0.025,   # 95% CI lower bound
                             0.100,   # lower shoulder resolution
                             0.160,   # 68% CI lower bound
                             0.250,   # lower interior (quartile)
                             0.500,   # median
                             0.750,   # upper interior (quartile)
                             0.840,   # 68% CI upper bound
                             0.900,   # upper shoulder resolution
                             0.975,   # 95% CI upper bound
                             0.990,   # upper tail interpolation
                             1)       # ensemble maximum         
          na_flag = TRUE
     
          # Run post-processing for gridded analysis
          dummy = post_process_for_grid(outfile_stock_fluxes,PROJECT,drivers,parameters,num_quantiles,na_flag,converged,states_all)
         
          # Optionally update the user
          if (use_parallel == FALSE) {print("Postprocessing for grid done, return completed filename to run_each_site")}
          # Tidy local environment
          rm(states_all,drivers) ; gc()
          # Return the site_output list object for subsequent slotting into the wider grid
          if (dummy == 0) { return(outfile_stock_fluxes) }

      } # gridded run?

      # Return, check just in case
      dummy = -5 ; return(dummy)

  } else { # *parameters.RData already exists

      # Report to user
      print('Already extracted result vectors (set repair = 1 if re-run is needed)')
      # Return
      dummy = 0 ; return(dummy)


  } # *parameters.RData already exists

} # end of run_each_site
## Use byte compile
run_each_site<-cmpfun(run_each_site)