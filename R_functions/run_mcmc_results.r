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
# Functions run and manage input / output data for DALEC runs (part of stage 3).
# 
# Author: T. Luke Smallman (12/11/2024)
# Exceptions states below in specific functions
#
#########################################################################################

rollapply_mean_annual_max<-function(var_in, step) {

   # Determine the mean of a maximums extracted from a rolling descrete windows
   # function specific command turned into a function to allow
   # for use within apply()
   return(mean(rollapply(var_in, by = step, width = step, FUN=max), na.rm = TRUE))

} # rollapply_mean_annual_max
## Use byte compile
rollapply_mean_annual_max<-cmpfun(rollapply_mean_annual_max)

rollapply_mean_annual<-function(var_in, step) {

   # Determine the mean of a maximums extracted from a rolling descrete windows
   # function specific command turned into a function to allow
   # for use within apply()
   return(rollapply(var_in, by = step, width = step, FUN=mean))

} # rollapply_mean_annual
## Use byte compile
rollapply_mean_annual<-cmpfun(rollapply_mean_annual)

###
## Function to create the output objected for a gridded CARDAMOM analysis
###

define_grid_output<-function(PROJECT,repair,outfile_grid,site_output){

      # Begin creation of all variables for the output gridded dataset
      # otherwise load the existing but incomplete version from file
      if (file.exists(outfile_grid) == FALSE | repair == 1) {

          print("Creating new grid_output object")

          #
          # Generate the summary information
          # Time invariant but contain uncertainty, shaped into the full spatial grid
          # i.e. including areas not part of the analysis but within the spatial domain
          #

          # Define grid_output object and add some basic information
          grid_output = list(readme = c("For gridded variables (e.g. mean_*, or final_*) the default dimesions are long/lat/quantile",
                                        "For time varying variables (e.g. wood_gCm2) default dimensions are site/quantile/time"),
                             num_quantiles = site_output$num_quantiles,
                             long_dim = PROJECT$long_dim, lat_dim = PROJECT$lat_dim, time_dim = dim(site_output$labile_gCm2)[2],
                             start_year = PROJECT$start_year, end_year = PROJECT$end_year,
                             steps_per_year = site_output$steps_per_year, nos_years = PROJECT$nos_years)

          # Load to local variables
          steps_per_year = site_output$steps_per_year ; nos_years = PROJECT$nos_years
          check_list = names(site_output)
 
          ###
          # Define likelihood / parameter / driver information
          ###

          # loop through parameters + likelihood
          grid_output$parameters = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)+1,dim(site_output$labile_gCm2)[1]))
          # track which parameters have converged + likelihood
          grid_output$parameters_converged = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)+1))
          # Parameter priors
          grid_output$parameter_priors_array = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)))
          grid_output$parameter_priors_uncertainty_array = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)))
          # Mean meteorological conditions
          grid_output$mean_temperature_C = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$mean_radiation_MJm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$mean_vpd_Pa = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$mean_precipitation_kgH2Om2yr = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          # Assimilated leaf area index information
          grid_output$assimilated_lai_max_m2m2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$assimilated_lai_mean_m2m2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$assimilated_lai_sd_m2m2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$assimilated_lai_unc_m2m2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          # Assimilated wood stock / prior information
          grid_output$assimilated_wood_mean_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$assimilated_wood_mean_unc_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          # Assimilated som stock / prior information
          grid_output$assimilated_som_mean_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$assimilated_som_mean_unc_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          # Generic dump of the whole driver$met and drivers$obs arrays
          grid_output$met_array_averages = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(met_array_names)))
          grid_output$obs_array_averages = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(obs_array_names)))
          grid_output$met_array_annual_averages = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_years,length(met_array_names)))
          grid_output$obs_array_annual_averages = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_years,length(obs_array_names)))

          # Net primary production allocation fractions
          if (any(check_list == "NPP_foliage_fraction") == TRUE) {grid_output$NPP_foliage_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (any(check_list == "NPP_roots_fraction") == TRUE) {grid_output$NPP_roots_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (any(check_list == "NPP_wood_fraction") == TRUE) {grid_output$NPP_wood_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          # Analysis mean transit (residence) times (years)
          if (any(check_list == "MTT_labile_years") == TRUE) {grid_output$MTT_labile_years =array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (any(check_list == "MTT_foliage_years") == TRUE) {grid_output$MTT_foliage_years = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (any(check_list == "MTT_roots_years") == TRUE) {grid_output$MTT_roots_years = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (any(check_list == "MTT_wood_years") == TRUE) {grid_output$MTT_wood_years = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (any(check_list == "MTT_litter_years") == TRUE) {grid_output$MTT_litter_years = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (any(check_list == "MTT_woodlitter_years") == TRUE) {grid_output$MTT_woodlitter_years = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (any(check_list == "MTT_som_years") == TRUE) {grid_output$MTT_som_years = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (any(check_list == "MTT_Ctotal_years") == TRUE) {grid_output$MTT_Ctotal_years = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (any(check_list == "MTT_annual_Ctotal_years") == TRUE) {grid_output$MTT_annual_Ctotal_years = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))}
          # Steady state C stock estimates (gC/m2)
          if (any(check_list == "SS_labile_gCm2") == TRUE) {grid_output$SS_labile_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (any(check_list == "SS_foliage_gCm2") == TRUE) {grid_output$SS_foliage_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (any(check_list == "SS_roots_gCm2") == TRUE) {grid_output$SS_roots_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (any(check_list == "SS_wood_gCm2") == TRUE) {grid_output$SS_wood_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (any(check_list == "SS_litter_gCm2") == TRUE) {grid_output$SS_litter_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (any(check_list == "SS_woodlitter_gCm2") == TRUE) {grid_output$SS_woodlitter_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}
          if (any(check_list == "SS_som_gCm2") == TRUE) {grid_output$SS_som_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))}

          # Total ecosystem C is always present
          grid_output$mean_Ctotal_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$final_Ctotal_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$final_dCtotal_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          # Ecosystem leaf area index is always present
          grid_output$mean_lai_m2m2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$final_dlai_m2m2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          # Mean ecosystem bulk fluxes are always present
          grid_output$mean_nee_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_gpp_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_rauto_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_rhet_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_reco_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_npp_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_harvest_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_fire_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_nbe_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_nbp_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          # Always present but at annual time step
          grid_output$mean_cue = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          # Time varying pixel based (i.e. not within the grid) values with quantile based uncertainty for...
          # States
          grid_output$Ctotal_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$dCtotal_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$lai_m2m2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$dlai_m2m2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          # Fluxes
          grid_output$nee_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$gpp_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$rauto_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$rhet_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$reco_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$npp_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$harvest_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$fire_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$nbe_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$nbp_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          # Always present but at annual time step
          grid_output$mean_annual_cue = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          grid_output$mean_annual_Ctotal_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          grid_output$mean_annual_lai_m2m2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          grid_output$mean_annual_nee_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          grid_output$mean_annual_gpp_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          grid_output$mean_annual_rauto_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          grid_output$mean_annual_rhet_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          grid_output$mean_annual_reco_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          grid_output$mean_annual_npp_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          grid_output$mean_annual_harvest_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          grid_output$mean_annual_fire_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          grid_output$mean_annual_nbe_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          grid_output$mean_annual_nbp_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          # Always present but at annual time step differences
          grid_output$mean_annual_dlai_m2m2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          grid_output$mean_annual_dCtotal_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          grid_output$mean_annual_dnee_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          grid_output$mean_annual_dgpp_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          grid_output$mean_annual_drauto_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          grid_output$mean_annual_drhet_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          grid_output$mean_annual_dreco_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          grid_output$mean_annual_dnpp_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          grid_output$mean_annual_dharvest_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          grid_output$mean_annual_dfire_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          grid_output$mean_annual_dnbe_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          grid_output$mean_annual_dnbp_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))

          if (any(check_list == "grazing_gCm2day")) {
              grid_output$grazing_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$mean_grazing_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_annual_grazing_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_dgrazing_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          }

          # Annual changes
          #grid_output$annual_change_Ctotal_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          #grid_output$annual_change_biomass_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          #grid_output$annual_change_dom_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))

          # Periodically request garbage collection to maintain memory management - to avoid R's attempts to claim more memory than exists
          gc() ; gc()

          # Based on the presence of each pool define the grids for the mean and final values.
          # Also, create the time varying but quantile based values and time

          # Gridded biomass information
          if (any(check_list == "biomass_gCm2") == TRUE) {
              # Grid mean / finals for globally available variables
              grid_output$mean_biomass_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_biomass_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_dCbiomass_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_outflux_biomass_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_combined_biomass_to_litter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_biomass_to_litter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$MTT_biomass_years =array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Time varying pixel specific with quantiles
              grid_output$biomass_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$dCbiomass_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$outflux_biomass_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$biomass_to_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$combined_biomass_to_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # Annual information
              grid_output$mean_annual_biomass_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_dCbiomass_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_outflux_biomass_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_combined_biomass_to_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_biomass_to_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$MTT_annual_biomass_years = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              # Fractional partitioning of turnover to different drivers - should they exist
              grid_output$NaturalFractionOfTurnover_biomass = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$FireFractionOfTurnover_biomass = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$HarvestFractionOfTurnover_biomass = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$GRAZINGFractionOfTurnover_biomass = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Conditional variables
              if (any(check_list == "FIREemiss_biomass_gCm2day") == TRUE) {
                  grid_output$FIREemiss_biomass_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIREemiss_biomass_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_FIREemiss_biomass_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "FIRElitter_biomass_gCm2day") == TRUE) {
                  grid_output$FIRElitter_biomass_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIRElitter_biomass_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_FIRElitter_biomass_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "HARVESTextracted_biomass_gCm2day") == TRUE) {
                  grid_output$HARVESTextracted_biomass_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTextracted_biomass_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_HARVESTextracted_biomass_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
                  grid_output$HARVESTlitter_biomass_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTlitter_biomass_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_HARVESTlitter_biomass_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "GRAZINGextracted_biomass_gCm2day") == TRUE) {
                  grid_output$GRAZINGextracted_biomass_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_GRAZINGextracted_biomass_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_GRAZINGextracted_biomass_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
                  grid_output$GRAZINGlitter_biomass_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_GRAZINGlitter_biomass_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_GRAZINGlitter_biomass_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }              
          }
          # Gridded labile information
          if (any(check_list == "labile_gCm2") == TRUE) {
              # Grid averages
              grid_output$mean_labile_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_labile_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_dClabile_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_outflux_labile_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_labile_to_foliage_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_alloc_labile_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_combined_labile_to_litter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Time varying pixel specific with quantiles
              grid_output$labile_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$dClabile_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$outflux_labile_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$labile_to_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$alloc_labile_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$combined_labile_to_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # Annual information
              grid_output$mean_annual_labile_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_dClabile_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))              
              grid_output$mean_annual_outflux_labile_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_labile_to_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_alloc_labile_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_combined_labile_to_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$MTT_annual_labile_years = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              # Fractional partitioning of turnover to different drivers - should they exist
              grid_output$NaturalFractionOfTurnover_labile = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$FireFractionOfTurnover_labile = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$HarvestFractionOfTurnover_labile = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$GRAZINGFractionOfTurnover_labile = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Conditional variables
              if (any(check_list == "FIREemiss_labile_gCm2day") == TRUE) {
                  grid_output$FIREemiss_labile_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIREemiss_labile_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_FIREemiss_labile_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "FIRElitter_labile_gCm2day") == TRUE) {
                  grid_output$FIRElitter_labile_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIRElitter_labile_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_FIRElitter_labile_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "HARVESTextracted_labile_gCm2day") == TRUE) {
                  grid_output$HARVESTextracted_labile_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTextracted_labile_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_HARVESTextracted_labile_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
                  grid_output$HARVESTlitter_labile_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTlitter_labile_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_HARVESTlitter_labile_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "GRAZINGextracted_labile_gCm2day") == TRUE) {
                  grid_output$GRAZINGextracted_labile_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_GRAZINGextracted_labile_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_GRAZINGextracted_labile_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
                  grid_output$GRAZINGlitter_labile_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_GRAZINGlitter_labile_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_GRAZINGlitter_labile_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }              
          }
          # Gridded foliage information
          if (any(check_list == "foliage_gCm2") == TRUE) {
              # Grided mean / final
              grid_output$mean_foliage_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_foliage_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_dCfoliage_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_outflux_foliage_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_foliage_to_litter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_combined_alloc_foliage_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_combined_foliage_to_litter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Time varying pixel specific with quantiles
              grid_output$foliage_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$dCfoliage_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$outflux_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$foliage_to_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$combined_alloc_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$combined_foliage_to_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # Annual information
              grid_output$mean_annual_foliage_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_dCfoliage_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))              
              grid_output$mean_annual_outflux_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_foliage_to_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_combined_alloc_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_combined_foliage_to_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$MTT_annual_foliage_years = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              # Fractional partitioning of turnover to different drivers - should they exist
              grid_output$NaturalFractionOfTurnover_foliage = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$FireFractionOfTurnover_foliage = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$HarvestFractionOfTurnover_foliage = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$GRAZINGFractionOfTurnover_foliage = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Conditional variables
              if (any(check_list == "alloc_foliage_gCm2day")) {
                  grid_output$alloc_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_alloc_foliage_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_alloc_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "FIREemiss_foliage_gCm2day")) {
                  grid_output$FIREemiss_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIREemiss_foliage_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_FIREemiss_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "FIRElitter_foliage_gCm2day") == TRUE) {
                  grid_output$FIRElitter_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIRElitter_foliage_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_FIRElitter_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "HARVESTextracted_foliage_gCm2day") == TRUE) {
                  grid_output$HARVESTextracted_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTextracted_foliage_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_HARVESTextracted_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
                  grid_output$HARVESTlitter_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTlitter_foliage_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_HARVESTlitter_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "GRAZINGextracted_foliage_gCm2day") == TRUE) {
                  grid_output$GRAZINGextracted_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_GRAZINGextracted_foliage_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_GRAZINGextracted_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
                  grid_output$GRAZINGlitter_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_GRAZINGlitter_foliage_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_GRAZINGlitter_foliage_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }              
          }
          # Gridded roots information
          if (any(check_list == "roots_gCm2") == TRUE) {
              # Gridded mean / final
              grid_output$mean_roots_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_roots_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_dCroots_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_outflux_roots_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_roots_to_litter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_alloc_roots_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$annual_max_roots_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_combined_alloc_roots_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_combined_roots_to_litter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Time varying pixel specific with quantiles
              grid_output$roots_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$dCroots_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$outflux_roots_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$roots_to_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$alloc_roots_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$combined_alloc_roots_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))              
              grid_output$combined_roots_to_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # Annual information
              grid_output$mean_annual_roots_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_dCroots_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))               
              grid_output$mean_annual_outflux_roots_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_roots_to_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_alloc_roots_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_combined_alloc_roots_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_combined_roots_to_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$MTT_annual_roots_years = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              # Fractional partitioning of turnover to different drivers - should they exist
              grid_output$NaturalFractionOfTurnover_roots = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$FireFractionOfTurnover_roots = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$HarvestFractionOfTurnover_roots = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$GRAZINGFractionOfTurnover_roots = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Is rooting depth calculate (m) by this model?
              if (any(check_list == "RootDepth_m") == TRUE) {
                  grid_output$RootDepth_m = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$dRootDepth_m = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_RootDepth_m = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$final_dRootDepth_m = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_RootDepth_m = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "FIREemiss_roots_gCm2day") == TRUE) {
                  grid_output$FIREemiss_roots_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIREemiss_roots_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_FIREemiss_roots_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "FIRElitter_roots_gCm2day") == TRUE) {
                  grid_output$FIRElitter_roots_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIRElitter_roots_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_FIRElitter_roots_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "HARVESTextracted_roots_gCm2day") == TRUE) {
                  grid_output$HARVESTextracted_roots_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTextracted_roots_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_HARVESTextracted_roots_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
                  grid_output$HARVESTlitter_roots_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTlitter_roots_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_HARVESTlitter_roots_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "GRAZINGextracted_roots_gCm2day") == TRUE) {
                  grid_output$GRAZINGextracted_roots_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_GRAZINGextracted_roots_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_GRAZINGextracted_roots_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
                  grid_output$GRAZINGlitter_roots_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_GRAZINGlitter_roots_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_GRAZINGlitter_roots_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }              
          }
          # Gridded wood information
          if (any(check_list == "wood_gCm2") == TRUE) {
              # Gridded mean / final variables
              grid_output$mean_wood_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_wood_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_dCwood_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_outflux_wood_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_wood_to_litter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_alloc_wood_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$annual_max_wood_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_combined_alloc_wood_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))             
              grid_output$mean_combined_wood_to_litter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Time varying pixel specific with quantiles
              grid_output$wood_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$dCwood_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$outflux_wood_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$wood_to_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$alloc_wood_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$combined_alloc_wood_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))                            
              grid_output$combined_wood_to_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # Annual information
              grid_output$mean_annual_wood_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_dCwood_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_outflux_wood_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_wood_to_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_alloc_wood_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_combined_alloc_wood_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_combined_wood_to_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$MTT_annual_wood_years = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              # Fractional partitioning of turnover to different drivers - should they exist
              grid_output$NaturalFractionOfTurnover_wood = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$FireFractionOfTurnover_wood = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$HarvestFractionOfTurnover_wood = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$GRAZINGFractionOfTurnover_wood = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Conditional variables
              if (any(check_list == "FIREemiss_wood_gCm2day") == TRUE) {
                  grid_output$FIREemiss_wood_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIREemiss_wood_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_FIREemiss_wood_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "FIRElitter_wood_gCm2day") == TRUE) {
                  grid_output$FIRElitter_wood_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIRElitter_wood_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_FIRElitter_wood_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "HARVESTextracted_wood_gCm2day") == TRUE) {
                  grid_output$HARVESTextracted_wood_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTextracted_wood_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_HARVESTextracted_wood_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
                  grid_output$HARVESTlitter_wood_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTlitter_wood_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_HARVESTlitter_wood_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "GRAZINGextracted_wood_gCm2day") == TRUE) {
                  grid_output$GRAZINGextracted_wood_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_GRAZINGextracted_wood_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_GRAZINGextracted_wood_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
                  grid_output$GRAZINGlitter_wood_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_GRAZINGlitter_wood_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_GRAZINGlitter_wood_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }              
          }
          # Gridded litter information
          if (any(check_list == "litter_gCm2") == TRUE) {
              # Gridded mean / final
              grid_output$mean_litter_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_litter_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_dClitter_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_outflux_litter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_litter_to_som_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_rhet_litter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_combined_litter_to_som_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Time varying pixel specific with quantiles
              grid_output$litter_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$dClitter_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$outflux_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$litter_to_som_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$rhet_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$combined_litter_to_som_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # Annual information
              grid_output$mean_annual_litter_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_dClitter_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_outflux_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_litter_to_som_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_rhet_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_drhet_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_combined_litter_to_som_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$MTT_annual_litter_years = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              # Fractional partitioning of turnover to different drivers - should they exist
              grid_output$NaturalFractionOfTurnover_litter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$FireFractionOfTurnover_litter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$HarvestFractionOfTurnover_litter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$GRAZINGFractionOfTurnover_litter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Conditional variables
              if (any(check_list == "FIREemiss_litter_gCm2day") == TRUE) {
                  grid_output$FIREemiss_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIREemiss_litter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_FIREemiss_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "FIRElitter_litter_gCm2day") == TRUE) {
                  grid_output$FIRElitter_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIRElitter_litter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_FIRElitter_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "HARVESTextracted_litter_gCm2day") == TRUE) {
                  grid_output$HARVESTextracted_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTextracted_litter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_HARVESTextracted_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "GRAZINGextracted_litter_gCm2day") == TRUE) {
                  grid_output$GRAZINGextracted_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_GRAZINGextracted_litter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_GRAZINGextracted_litter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }              
          }
          # Gridded wood litter information
          if (any(check_list == "woodlitter_gCm2") == TRUE) {
              # Gridded mean / final
              grid_output$mean_woodlitter_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_woodlitter_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_dCwoodlitter_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_outflux_woodlitter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_woodlitter_to_som_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_rhet_woodlitter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_combined_woodlitter_to_som_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Time varying pixel specific with quantiles
              grid_output$woodlitter_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$dCwoodlitter_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$outflux_woodlitter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$woodlitter_to_som_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$rhet_woodlitter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$combined_woodlitter_to_som_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # Annual information
              grid_output$mean_annual_woodlitter_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_dCwoodlitter_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_outflux_woodlitter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_woodlitter_to_som_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_rhet_woodlitter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_combined_woodlitter_to_som_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$MTT_annual_woodlitter_years = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              # Fractional partitioning of turnover to different drivers - should they exist
              grid_output$NaturalFractionOfTurnover_woodlitter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$FireFractionOfTurnover_woodlitter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$HarvestFractionOfTurnover_woodlitter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$GRAZINGFractionOfTurnover_woodlitter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Conditional variables
              if (any(check_list == "FIREemiss_woodlitter_gCm2day") == TRUE) {
                  grid_output$FIREemiss_woodlitter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIREemiss_woodlitter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_FIREemiss_woodlitter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "FIRElitter_woodlitter_gCm2day") == TRUE) {
                  grid_output$FIRElitter_woodlitter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIRElitter_woodlitter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_FIRElitter_woodlitter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "HARVESTextracted_woodlitter_gCm2day") == TRUE) {
                  grid_output$HARVESTextracted_woodlitter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTextracted_woodlitter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_HARVESTextracted_woodlitter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "GRAZINGextracted_woodlitter_gCm2day") == TRUE) {
                  grid_output$GRAZINGextracted_woodlitter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_GRAZINGextracted_woodlitter_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_GRAZINGextracted_woodlitter_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }              
          }
          # Gridded som information
          if (any(check_list == "som_gCm2") == TRUE) {
              # Gridded mean / finals
              grid_output$mean_som_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_som_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_dCsom_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_outflux_som_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_rhet_som_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Time varying pixel specific with quantiles
              grid_output$som_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$dCsom_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$outflux_som_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$rhet_som_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # Annual information
              grid_output$mean_annual_som_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_dCsom_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_outflux_som_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_rhet_som_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_drhet_som_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$MTT_annual_som_years = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              # Fractional partitioning of turnover to different drivers - should they exist
              grid_output$NaturalFractionOfTurnover_som = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$FireFractionOfTurnover_som = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$HarvestFractionOfTurnover_som = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$GRAZINGFractionOfTurnover_som = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Time varying
              if (any(check_list == "FIREemiss_som_gCm2day") == TRUE) {
                  grid_output$FIREemiss_som_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIREemiss_som_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_FIREemiss_som_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "HARVESTextracted_som_gCm2day") == TRUE) {
                  grid_output$HARVESTextracted_som_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTextracted_som_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_HARVESTextracted_som_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "GRAZINGextracted_som_gCm2day") == TRUE) {
                  grid_output$GRAZINGextracted_som_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_GRAZINGextracted_som_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_GRAZINGextracted_som_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }              
          }
          # Gridded dom information
          if (any(check_list == "dom_gCm2") == TRUE) {
              # Gridded mean / finals
              grid_output$mean_dom_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_dom_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_dCdom_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_outflux_dom_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_rhet_dom_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$MTT_dom_years = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Time varying pixel specific with quantiles
              grid_output$dom_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$dCdom_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$outflux_dom_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$rhet_dom_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # Annual information
              grid_output$mean_annual_dom_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_dCdom_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_outflux_dom_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_rhet_dom_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$MTT_annual_dom_years = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              # Fractional partitioning of turnover to different drivers - should they exist
              grid_output$NaturalFractionOfTurnover_dom = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$FireFractionOfTurnover_dom = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$HarvestFractionOfTurnover_dom = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$GRAZINGFractionOfTurnover_dom = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # Conditional
              if (any(check_list == "FIREemiss_dom_gCm2day") == TRUE) {
                  grid_output$FIREemiss_dom_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_FIREemiss_dom_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_FIREemiss_dom_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "HARVESTextracted_dom_gCm2day") == TRUE) {
                  grid_output$HARVESTextracted_dom_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_HARVESTextracted_dom_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_HARVESTextracted_dom_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }
              if (any(check_list == "GRAZINGextracted_dom_gCm2day") == TRUE) {
                  grid_output$GRAZINGextracted_dom_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_GRAZINGextracted_dom_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$mean_annual_GRAZINGextracted_dom_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              }              
          }

          # Periodically request garbage collection to maintain memory management - to avoid R's attempts to claim more memory than exists
          gc() ; gc()
          
          # Water cycle specific variables
          if (any(check_list == "ET_kgH2Om2day") == TRUE) {
              # currently water in the soil surface layer (0-30 cm)
              grid_output$mean_SurfWater_kgH2Om2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_annual_SurfWater_kgH2Om2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_dSurfWater_kgH2Om2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$final_SurfWater_kgH2Om2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_dSurfWater_kgH2Om2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$SurfWater_kgH2Om2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$dSurfWater_kgH2Om2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # plant apparent soil water potential (MPa)
              grid_output$mean_wSWP_MPa = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_annual_wSWP_MPa = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_dwSWP_MPa = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$final_wSWP_MPa = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$wSWP_MPa = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$dwSWP_MPa = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # evapotranspiration (Etrans + Esoil + Ewetcanopy)
              grid_output$mean_annual_dET_kgH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_ET_kgH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_ET_kgH2Om2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$ET_kgH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # Ecosystem water use efficiency (GPP/ET)
              grid_output$mean_annual_wue_eco_gCkgH2O = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_wue_eco_gCkgH2O = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$wue_eco_gCkgH2O = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # Check whether the evaporation components exist
              if (any(check_list == "Etrans_kgH2Om2day") == TRUE) {
                  # Transpiration
                  grid_output$mean_annual_dEtrans_kgH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
                  grid_output$mean_annual_Etrans_kgH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
                  grid_output$mean_Etrans_kgH2Om2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$Etrans_kgH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  # Plant water use efficiency (GPP/Etrans)
                  grid_output$mean_annual_wue_plant_gCkgH2O = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
                  grid_output$mean_wue_plant_gCkgH2O = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$wue_plant_gCkgH2O = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              }
              if (any(check_list == "Esoil_kgH2Om2day") == TRUE) {
                  # Soil evaporation
                  grid_output$mean_annual_dEsoil_kgH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
                  grid_output$mean_annual_Esoil_kgH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
                  grid_output$mean_Esoil_kgH2Om2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$Esoil_kgH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              }
              if (any(check_list == "Ewetcanopy_kgH2Om2day") == TRUE) {
                  # Wet canopy evaporation
                  grid_output$mean_annual_dEwetcanopy_kgH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
                  grid_output$mean_annual_Ewetcanopy_kgH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
                  grid_output$mean_Ewetcanopy_kgH2Om2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$Ewetcanopy_kgH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              }
              if (any(check_list == "runoff_kgH2Om2day") == TRUE) {
                  # Surface water runoff
                  grid_output$mean_annual_drunoff_kgH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
                  grid_output$mean_annual_runoff_kgH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
                  grid_output$mean_runoff_kgH2Om2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$runoff_kgH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              }
              if (any(check_list == "underflow_kgH2Om2day") == TRUE) {
                  # Underflow from bottom of soil water column
                  grid_output$mean_annual_dunderflow_kgH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
                  grid_output$mean_annual_underflow_kgH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
                  grid_output$mean_underflow_kgH2Om2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$underflow_kgH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              }
              if (any(check_list == "total_drainage_kgH2Om2day") == TRUE) {
                  # Total drainage from soil surface andn bottom of soil water column
                  grid_output$mean_annual_dtotal_drainage_kgH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
                  grid_output$mean_annual_total_drainage_kgH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
                  grid_output$mean_total_drainage_kgH2Om2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
                  grid_output$total_drainage_kgH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              }
              if (any(check_list == 'LWP_MPa') == TRUE){
                  # Leaf water potential
                  grid_output$LWP_MPa = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
                  grid_output$mean_LWP_MPa = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              }
          }
          # Snow specific
          if (any(check_list == "snow_kgH2Om2") == TRUE) {
              ## snow on soil surface
              grid_output$snow_kgH2Om2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$mean_snow_kgH2Om2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_annual_snow_kgH2Om2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_dsnow_kgH2Om2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
          }
          # Canopy process variables
          if (any(check_list == "APAR_MJm2day") == TRUE) {
              # Absorbed photosynthetically active radation
              grid_output$mean_annual_APAR_MJm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_dAPAR_MJm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_APAR_MJm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$APAR_MJm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          }
          if (any(check_list == "CiCa") == TRUE) {
              # Canopy Ci:Ca
              grid_output$mean_annual_CiCa = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_dCiCa = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_CiCa = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$CiCa = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          }
          if (any(check_list == "gs_demand_supply_ratio") == TRUE) {
              # Ratio of stomatal conductance relative to its maximum value,
              # this metric provides information on the demand vs supply constrains on stomatal conductance
              grid_output$mean_annual_gs_demand_supply_ratio = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_annual_dgs_demand_supply_ratio = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_gs_demand_supply_ratio = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$gs_demand_supply_ratio = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          }
          if (any(check_list == "gs_mmolH2Om2s") == TRUE) {
              # Canopy stomatal conductance
              grid_output$mean_annual_gs_mmolH2Om2s = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_gs_mmolH2Om2s = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$gs_mmolH2Om2s = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          }
          if (any(check_list == "gb_mmolH2Om2s") == TRUE) {
              # Canopy boundary layer conductance
              grid_output$mean_annual_gb_mmolH2Om2s = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_gb_mmolH2Om2s = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$gb_mmolH2Om2s = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          }
          if (any(check_list == "leaf_temperature_celcius") == TRUE) {
              # Canopy temperature 
              grid_output$mean_annual_leaf_temperature_celcius = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_leaf_temperature_celcius = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$leaf_temperature_celcius = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # Soil temperature
              grid_output$mean_annual_soil_temperature_celcius = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],nos_years))
              grid_output$mean_soil_temperature_celcius = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$soil_temperature_celcius = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          }

          # Create overlap statistics variables - may not always get filled in the end
          grid_output$priors_assim_data_overlap_fraction = array(NA, dim=c(max(PROJECT$model$nopars),PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$gpp_assim_data_overlap_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$lai_assim_data_overlap_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$nee_assim_data_overlap_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$wood_assim_data_overlap_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$soil_assim_data_overlap_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$fapar_assim_data_overlap_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$et_assim_data_overlap_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$nbe_assim_data_overlap_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$fire_assim_data_overlap_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))

          # Time and uncertainty invarient information,
          # this is the correlation between ensemble members for parameter and C-cycle flux variables
          grid_output$lai_parameter_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)))
          grid_output$nbp_parameter_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)))
          grid_output$nee_parameter_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)))
          grid_output$gpp_parameter_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)))
          grid_output$rauto_parameter_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)))
          grid_output$rhet_parameter_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)))
          grid_output$fire_parameter_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)))
          # Assess within pixel correlations with LAI
          grid_output$lai_m2m2_to_GPP_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$lai_m2m2_to_NEE_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$lai_m2m2_to_NBP_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$lai_m2m2_to_Rauto_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$lai_m2m2_to_Rhet_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$lai_m2m2_to_wood_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$lai_m2m2_to_som_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$lai_m2m2_to_dCwood_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$lai_m2m2_to_dCsom_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))     
          # Assess within pixel correlations with NBP
          grid_output$NBP_gCm2day_to_GPP_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))     
          grid_output$NBP_gCm2day_to_NEE_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))     
          grid_output$NBP_gCm2day_to_lai_m2m2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))     
          grid_output$NBP_gCm2day_to_Rauto_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))     
          grid_output$NBP_gCm2day_to_Rhet_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))     
          grid_output$NBP_gCm2day_to_wood_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))     
          grid_output$NBP_gCm2day_to_som_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))     
          grid_output$NBP_gCm2day_to_dCwood_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))     
          grid_output$NBP_gCm2day_to_dCsom_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))     
          # If harvest 
          if (any(check_list == "harvest_gCm2day")) {
              grid_output$lai_m2m2_to_harvest_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$NBP_gCm2day_to_harvest_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          }                                        
          # If CiCa
          if (any(check_list == "CiCa_parameter_correlation") == TRUE) {
              grid_output$CiCa_parameter_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)))
          }
          # If Leaf Water Potential
          if (any(check_list == "LWP_parameter_correlation") == TRUE) {
              grid_output$LWP_parameter_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)))
          }
          # If Mean transit time for wood correlation exists, ensure we store it for the gridded run too
          if (any(check_list == "MTT_wood_years_parameter_correlation") == TRUE) {
              grid_output$MTT_wood_years_parameter_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)))
          }
          # If Mean mean allocation to wood correlation exists, ensure we store it for the gridded run too
          if (any(check_list == "NPP_wood_gCm2day_parameter_correlation") == TRUE) {
              grid_output$NPP_wood_gCm2day_parameter_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)))
          }
          # If the combined correlation between wood MTT and wood allocation has been calculated
          if (any(check_list == "MTT_wood_years_to_NPP_wood_gCm2day_correlation") == TRUE) {
              grid_output$MTT_wood_years_to_NPP_wood_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$MTT_wood_years_to_NPP_wood_fraction_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$MTT_wood_years_to_MTT_som_years_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              # Assess within pixel correlations with NPP flux to wood
              grid_output$NPP_wood_gCm2day_to_GPP_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$NPP_wood_gCm2day_to_NEE_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$NPP_wood_gCm2day_to_Rauto_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$NPP_wood_gCm2day_to_Rhet_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$NPP_wood_gCm2day_to_wood_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$NPP_wood_gCm2day_to_som_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$NPP_wood_gCm2day_to_lai_m2m2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$NPP_wood_gCm2day_to_dCwood_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$NPP_wood_gCm2day_to_dCsom_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))              
              # Assess within pixel correlations with NPP fraction to wood
              grid_output$NPP_wood_fraction_to_GPP_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$NPP_wood_fraction_to_NEE_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$NPP_wood_fraction_to_Rauto_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$NPP_wood_fraction_to_Rhet_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$NPP_wood_fraction_to_wood_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$NPP_wood_fraction_to_som_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$NPP_wood_fraction_to_lai_m2m2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$NPP_wood_fraction_to_dCwood_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$NPP_wood_fraction_to_dCsom_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              # Assess within pixel correlations with wood turnover
              grid_output$MTT_wood_years_to_GPP_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$MTT_wood_years_to_NEE_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$MTT_wood_years_to_Rauto_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$MTT_wood_years_to_Rhet_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$MTT_wood_years_to_wood_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$MTT_wood_years_to_som_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$MTT_wood_years_to_lai_m2m2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$MTT_wood_years_to_dCwood_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$MTT_wood_years_to_dCsom_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))    
              grid_output$dCwood_gCm2_to_gpp_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))    
              grid_output$dCwood_gCm2_to_rauto_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))    
              grid_output$dCwood_gCm2_to_nee_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))    
              grid_output$dCwood_gCm2_to_rhet_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))    
              grid_output$dCwood_gCm2_to_wood_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))    
              grid_output$dCwood_gCm2_to_som_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))    
              grid_output$dCwood_gCm2_to_dCsom_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))                  
          }
          # If Mean mean allocation to wood correlation exists, ensure we store it for the gridded run too
          if (any(check_list == "MTT_som_years_parameter_correlation") == TRUE) {
              grid_output$MTT_som_years_parameter_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)))
              # Assess within pixel correlations with soil turnover
              grid_output$MTT_som_years_to_GPP_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$MTT_som_years_to_NEE_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$MTT_som_years_to_Rauto_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$MTT_som_years_to_Rhet_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$MTT_som_years_to_wood_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$MTT_som_years_to_som_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$MTT_som_years_to_lai_m2m2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$MTT_som_years_to_dCwood_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
              grid_output$MTT_som_years_to_dCsom_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))             
              grid_output$dCsom_gCm2_to_gpp_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))    
              grid_output$dCsom_gCm2_to_rauto_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))    
              grid_output$dCsom_gCm2_to_nee_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))    
              grid_output$dCsom_gCm2_to_rhet_gCm2day_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))    
              grid_output$dCsom_gCm2_to_wood_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))    
              grid_output$dCsom_gCm2_to_som_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))    
              grid_output$dCsom_gCm2_to_som_input_gCm2_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))    
          }
          # Quantify the mean absolute magnitude of correlations between parameters
          grid_output$absolute_mean_parameter_correlation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))

          #
          # Generate the spatial information needed to relate i,j within the grid to long / lat and the summary vs detailed output
          #

          # vector into which the i,j position information within the grid will be stored
          grid_output$i_location = rep(NA, length.out = PROJECT$nosites)
          grid_output$j_location = rep(NA, length.out = PROJECT$nosites)

          # generate the lat / long grid again
          output = generate_grid(PROJECT$grid_type,PROJECT$latitude,PROJECT$longitude,PROJECT$resolution)
          grid_output$lat = output$lat
          grid_output$long = output$long
          # Structure area (in m) as a grid
          grid_output$area_m2 = output$area

          # Load the land mask...
          grid_output$landmask = array(0, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          # ...and land fraction
          grid_output$land_fraction = array(PROJECT$landsea, dim=c(PROJECT$long_dim,PROJECT$lat_dim))         

      } else {

          # if file already exists load it up
          load(outfile_grid)

      } # have we already an output file

      # Periodically request garbage collection to maintain memory management - to avoid R's attempts to claim more memory than exists
      gc() ; gc()

      # Return
      return(grid_output)

} # define_grid_output
## Use byte compile
define_grid_output<-cmpfun(define_grid_output)

###
## Function to control the site level running of the re-processing by DALEC
###

run_mcmc_results <- function (PROJECT,repair,grid_override) {

  # Update the user
  print('Welcome to RUN_MCMC_RESULTS!!')

  # start marker
  stime = proc.time()["elapsed"]

  # how many plots in total do we have
  nos_plots = 1:PROJECT$nosites

  # Calculate the number of years
  PROJECT$nos_years = (as.numeric(PROJECT$end_year) - as.numeric(PROJECT$start_year))+1
  Sys.sleep(1) # wait to allow the memory to be updates?
  save(PROJECT, file = paste(PROJECT$localpath,"infofile.RData",sep=""))

  # determine what the output file name is here, so that we can check if one already exists
  outfile_grid = paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep="")

  # now check which ones we need to calculate, but only if override not in play
  existing_files = rep(NA, length(nos_plots)) ; existing_list = 0
  if (repair != 1) {
      # Inform the user
      print("...beginning filterings for sites we have already processed")
      # Generate lists of expected files, which changes depending on whether this is a gridded or a site based run
      if (PROJECT$spatial_type == "grid") {
          outfile_stocks = paste(PROJECT$results_processedpath,PROJECT$sites,"_stock_fluxes.RData",sep="")
      } else if (PROJECT$spatial_type == "site") {
          outfile_stocks = paste(PROJECT$results_processedpath,PROJECT$sites,".RData",sep="")
      } else {
          stop("PROJECT$spatial_type not of valid value, i.e. grid or site")
      }
      # Determine whether each exists
      outfile_stocks_exists = file.exists(outfile_stocks)
      # Determine whether each exists
      outfile_stocks_exists = file.exists(outfile_stocks)
      keep_list = which(outfile_stocks_exists == FALSE) # i.e. keep in the to do list
      existing_list = which(outfile_stocks_exists)
      existing_files[which(outfile_stocks_exists)] = outfile_stocks[which(outfile_stocks_exists)] # already existing output files
      # Update user
      print(paste("......removing ",length(nos_plots)-length(keep_list)," sites out of ",length(nos_plots)," from the analysis",sep=""))
      nos_plots = nos_plots[keep_list]
  } # repair !=1
  existing_list = existing_list[-1] # do not remove - see below useage for why

  # now request the creation of the plots
  if (length(nos_plots) > 1 & request_use_local_slurm) {

      # use parallel in interactive mode

      # Inform the user
      print("...beginning parallel operations using slurm mode")

      # NOTE: by creating submitting slurm jobs to process each site we are firmly restricted
      # a linux operating system, and one using a slurm job scheduler

      # Create a job ID which will be common to all associated tasks
      job_ID = floor(as.vector(Sys.time())) # this will be unique as the number of seconds since some reference start datetime
      # Now we can deploy in anger
      dummy = submit_R_run_each_site_to_local_slurm_machine(PROJECT,repair,job_ID)

      # Create empty output object into which to insert the required file names from the existing files
      site_output_all = vector("list", PROJECT$nosites)
      for (n in 1:PROJECT$nosites) {
           # Search the processed results folder for *_stock_fluxes.RData files
           tmp = paste(PROJECT$results_processedpath,PROJECT$sites[n],"_stock_fluxes.RData",sep="")
           if (file.exists(tmp)) { site_output_all[[n]] = tmp }
      } 
      # site_output_all = a list object of the filepaths to all files successfully ran

      print("...finished parallel operations using slurm mode")

  } else if (use_parallel & length(nos_plots) > 1) {

      # use parallel in interactive mode

      # Inform the user
      print("...beginning parallel operations using interactive mode")

      # NOTE: that the use of mclapply() is due to reported improved efficiency over creating a virtual cluster.
      # However, mclapply does not (at the time of typing) work on Windows, i.e. Linux and Mac only

      # Now we can deploy in anger
      cl <- min(length(nos_plots),numWorkers)
      site_output_all = mclapply(nos_plots,FUN=run_each_site,PROJECT=PROJECT,
                                 repair=repair,grid_override=grid_override, mc.cores = cl)

      print("...finished parallel operations using interactive mode")

  } else if (length(nos_plots) > 0) {

      # or use serial

      # Inform the user
      print("...beginning serial operations")

      # Run analysis
      site_output_all = lapply(nos_plots,FUN=run_each_site,PROJECT=PROJECT,
                               repair=repair,grid_override=grid_override)

      print("...finished serial operations")

  } else {

      # Create empty output object into which to insert the required file names from the existing files
      site_output_all = vector("list", PROJECT$nosites)

  } # parallel option

  # Check whether we have some existing files...
  if (length(existing_list) > 0) {
      if (existing_list[1] > 0) {
          # ...we do, so insert them into the overall output file list
          for (n in 1:length(existing_list)) {
               site_output_all[[existing_list[n]]] = existing_files[existing_list[n]]
          }
      } # More subtle control
  } # does the variable have a value?

  # now if this is a gridded run we want to take out individual site specific summary files and combine them into a single file
  if (PROJECT$spatial_type == "grid" & grid_override == FALSE) {

      # Load individual site information into the combined grid

      if (file.exists(outfile_grid) == FALSE || repair == 1) {
          # Load the first completed site_output file to great the grid_output
          n = 0 ; site_output = -1
          while (class(site_output) != "character" & n < length(site_output_all)) {
                 n = n + 1
                 site_output = site_output_all[[n]]
          }
          # Check if the loop has finished with error
          if (n == length(site_output_all) & class(site_output_all[[n]]) != "character") {
              # Difficult to say what will come out of here, so dump it all!
              print(unlist(site_output_all))
              print("Above error from run_mcmc_results.r L1018, problem with site_output_all")
              return(c(-1))
          }
          # Now assuming this has worked correctly, we can create the grid_output
          # First by loading the selected file
          load(site_output)
          # then creating the list object
          grid_output = define_grid_output(PROJECT,repair,outfile_grid,site_output)
          # tidy up
          rm(site_output)
      } else {
          # Load the existing file
          load(outfile_grid)
      }

      ###
      # Now place within the overall grid
      ###

      # Load into the grid_output object
      grid_output = post_process_into_grid(grid_output,site_output_all,PROJECT)

      # update the user
      print("...writing combined grid_output to file")

      # now save the combined grid file
      save(grid_output, file = outfile_grid, compress = "gzip", compression_level = 9)

      # Tidy up
      rm(grid_output) ; gc(reset=TRUE)

  } # gridded run?

  # tell me whats happening
  print(paste("...time to process ",round((proc.time()["elapsed"]-stime)/60,1)," minutes",sep=""))
  # Tidy up
  rm(site_output_all) ; gc(reset=TRUE)

  # return
  dummy = 0 ; return(dummy)

} # end function run_mcmc_results

## Use byte compile
run_mcmc_results<-cmpfun(run_mcmc_results)
