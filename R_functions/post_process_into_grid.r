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
# Function to insert post-processed CARDAMOM output into a gridded output variable
# 
# Author: T. Luke Smallman (02/05/2024)
#
#########################################################################################

post_process_into_grid<-function(grid_output,site_output_all,PROJECT) {

  # Loop through all sites
  for (n in 1:PROJECT$nosites) {

       # Update the user
       if (n%%10000 == 0) {print(paste("...",floor((n/PROJECT$nosites)*100),"% completed",sep=""))}

       # Extract current site file name from output object
       site_output = site_output_all[[n]]

       # Check that the file name is a character string, and we will assume it
       # exists
       if (class(site_output) == "character") {

           # Load the file
           load(site_output)

           # Update the list variables in states_all which we will be searching
           check_list = names(site_output)

           # determine the lat / long location within the grid
           slot_j = as.numeric(PROJECT$sites[n])/PROJECT$long_dim
           slot_i = as.numeric(PROJECT$sites[n])-(floor(slot_j)*PROJECT$long_dim)
           if(slot_i == 0) {slot_i = PROJECT$long_dim} ; slot_j = ceiling(slot_j)
           # save for later
           grid_output$i_location[n] = slot_i ; grid_output$j_location[n] = slot_j
           # Update the landmask to specify locations with actual anlysis
           grid_output$landmask[slot_i,slot_j] = 1

           # loop through parameters + likelihood
           grid_output$parameters[slot_i,slot_j,,] = site_output$parameters
           # track which parameters have converged + likelihood
           grid_output$parameters_converged[slot_i,slot_j,] = site_output$parameters_converged
           # Parameter priors
           grid_output$parameter_priors_array[slot_i,slot_j,] = site_output$parameter_priors_array
           grid_output$parameter_priors_uncertainty_array[slot_i,slot_j,] = site_output$parameter_priors_uncertainty_array
           # Generic dump of the whole driver$met and drivers$obs arrays
           # Long term average first...
           grid_output$met_array_averages[slot_i,slot_j,] = site_output$met_array_averages
           grid_output$obs_array_averages[slot_i,slot_j,] = site_output$obs_array_averages
           # Annual averages second
           grid_output$met_array_annual_averages[slot_i,slot_j,,] = site_output$met_array_annual_averages
           grid_output$obs_array_annual_averages[slot_i,slot_j,,] = site_output$obs_array_annual_averages           
           # Mean meteorological conditions
           grid_output$mean_temperature_C[slot_i,slot_j] = site_output$mean_temperature_C
           grid_output$mean_radiation_MJm2day[slot_i,slot_j] = site_output$mean_radiation_MJm2day
           grid_output$mean_vpd_Pa[slot_i,slot_j] = site_output$mean_vpd_Pa
           grid_output$mean_precipitation_kgH2Om2yr[slot_i,slot_j] = site_output$mean_precipitation_kgH2Om2yr
           # Assimilated leaf area index information
           grid_output$assimilated_lai_max_m2m2[slot_i,slot_j] = site_output$assimilated_lai_max_m2m2
           grid_output$assimilated_lai_mean_m2m2[slot_i,slot_j] = site_output$assimilated_lai_mean_m2m2
           grid_output$assimilated_lai_sd_m2m2[slot_i,slot_j] = site_output$assimilated_lai_sd_m2m2
           grid_output$assimilated_lai_unc_m2m2[slot_i,slot_j] = site_output$assimilated_lai_unc_m2m2
           # Assimilated wood stock / prior information
           grid_output$assimilated_wood_mean_gCm2[slot_i,slot_j] = site_output$assimilated_wood_mean_gCm2
           grid_output$assimilated_wood_mean_unc_gCm2[slot_i,slot_j] = site_output$assimilated_wood_mean_unc_gCm2
           # Assimilated som stock / prior information
           grid_output$assimilated_som_mean_gCm2[slot_i,slot_j] = site_output$assimilated_som_mean_gCm2
           grid_output$assimilated_som_mean_unc_gCm2[slot_i,slot_j] = site_output$assimilated_som_mean_unc_gCm2

           # Net primary production allocation fractions
           if (any(check_list == "NPP_foliage_fraction")) {grid_output$NPP_foliage_fraction[slot_i,slot_j,] = site_output$NPP_foliage_fraction}
           if (any(check_list == "NPP_roots_fraction")) {grid_output$NPP_roots_fraction[slot_i,slot_j,] = site_output$NPP_roots_fraction}
           if (any(check_list == "NPP_wood_fraction")) {grid_output$NPP_wood_fraction[slot_i,slot_j,] = site_output$NPP_wood_fraction}
           # Analysis mean transit (residence) times (years)
           if (any(check_list == "MTT_labile_years")) {grid_output$MTT_labile_years[slot_i,slot_j,] = site_output$MTT_labile_years}
           if (any(check_list == "MTT_foliage_years")) {grid_output$MTT_foliage_years[slot_i,slot_j,] = site_output$MTT_foliage_years}
           if (any(check_list == "MTT_roots_years")) {grid_output$MTT_roots_years[slot_i,slot_j,] = site_output$MTT_roots_years}
           if (any(check_list == "MTT_wood_years")) {grid_output$MTT_wood_years[slot_i,slot_j,] = site_output$MTT_wood_years}
           if (any(check_list == "MTT_litter_years")) {grid_output$MTT_litter_years[slot_i,slot_j,] = site_output$MTT_litter_years}
           if (any(check_list == "MTT_woodlitter_years")) {grid_output$MTT_woodlitter_years[slot_i,slot_j,] = site_output$MTT_woodlitter_years}
           if (any(check_list == "MTT_som_years")) {grid_output$MTT_som_years[slot_i,slot_j,] = site_output$MTT_som_years}
           if (any(check_list == "MTT_Ctotal_years")) {grid_output$MTT_Ctotal_years[slot_i,slot_j,] = site_output$MTT_Ctotal_years}
           if (any(check_list == "MTT_annual_Ctotal_years")) {grid_output$MTT_annual_Ctotal_years[n,,] = site_output$MTT_annual_Ctotal_years}
           # Steady state C stock estimates (gC/m2)
           if (any(check_list == "SS_labile_gCm2")) {grid_output$SS_labile_gCm2[slot_i,slot_j,] = site_output$SS_labile_gCm2}
           if (any(check_list == "SS_foliage_gCm2")) {grid_output$SS_foliage_gCm2[slot_i,slot_j,] = site_output$SS_foliage_gCm2}
           if (any(check_list == "SS_roots_gCm2")) {grid_output$SS_roots_gCm2[slot_i,slot_j,] = site_output$SS_roots_gCm2}
           if (any(check_list == "SS_wood_gCm2")) {grid_output$SS_wood_gCm2[slot_i,slot_j,] = site_output$SS_wood_gCm2}
           if (any(check_list == "SS_litter_gCm2")) {grid_output$SS_litter_gCm2[slot_i,slot_j,] = site_output$SS_litter_gCm2}
           if (any(check_list == "SS_woodlitter_gCm2")) {grid_output$SS_woodlitter_gCm2[slot_i,slot_j,] = site_output$SS_woodlitter_gCm2}
           if (any(check_list == "SS_som_gCm2")) {grid_output$SS_som_gCm2[slot_i,slot_j,] = site_output$SS_som_gCm2}

           # Total ecosystem C is always present
           grid_output$mean_Ctotal_gCm2[slot_i,slot_j,] = site_output$mean_Ctotal_gCm2
           grid_output$final_Ctotal_gCm2[slot_i,slot_j,] = site_output$Ctotal_gCm2[,grid_output$time_dim]
           grid_output$final_dCtotal_gCm2[slot_i,slot_j,] = site_output$dCtotal_gCm2[,grid_output$time_dim]
           # Ecosystem leaf area index is always present
           grid_output$mean_lai_m2m2[slot_i,slot_j,] = site_output$mean_lai_m2m2
           grid_output$final_dlai_m2m2[slot_i,slot_j,] = site_output$dlai_m2m2[,grid_output$time_dim]
           # Mean ecosystem bulk fluxes are always present
           grid_output$mean_nee_gCm2day[slot_i,slot_j,] = site_output$mean_nee_gCm2day
           grid_output$mean_gpp_gCm2day[slot_i,slot_j,] = site_output$mean_gpp_gCm2day
           grid_output$mean_rauto_gCm2day[slot_i,slot_j,] = site_output$mean_rauto_gCm2day
           grid_output$mean_rhet_gCm2day[slot_i,slot_j,] = site_output$mean_rhet_gCm2day
           grid_output$mean_reco_gCm2day[slot_i,slot_j,] = site_output$mean_reco_gCm2day
           grid_output$mean_npp_gCm2day[slot_i,slot_j,] = site_output$mean_npp_gCm2day
           grid_output$mean_harvest_gCm2day[slot_i,slot_j,] = site_output$mean_harvest_gCm2day
           grid_output$mean_fire_gCm2day[slot_i,slot_j,] = site_output$mean_fire_gCm2day
           grid_output$mean_nbe_gCm2day[slot_i,slot_j,] = site_output$mean_nbe_gCm2day
           grid_output$mean_nbp_gCm2day[slot_i,slot_j,] = site_output$mean_nbp_gCm2day
           # Always present but at annual time step
           grid_output$mean_cue[slot_i,slot_j,] = site_output$mean_cue
           # Time varying pixel based (i.e. not within the grid) values with quantile based uncertainty for...
           # States
           grid_output$Ctotal_gCm2[n,,]  = site_output$Ctotal_gCm2
           grid_output$dCtotal_gCm2[n,,] = site_output$dCtotal_gCm2
           grid_output$lai_m2m2[n,,]     = site_output$lai_m2m2
           grid_output$dlai_m2m2[n,,]    = site_output$dlai_m2m2
           # Fluxes
           grid_output$nee_gCm2day[n,,]     = site_output$nee_gCm2day
           grid_output$gpp_gCm2day[n,,]     = site_output$gpp_gCm2day
           grid_output$rauto_gCm2day[n,,]   = site_output$rauto_gCm2day
           grid_output$rhet_gCm2day[n,,]    = site_output$rhet_gCm2day
           grid_output$reco_gCm2day[n,,]    = site_output$reco_gCm2day
           grid_output$npp_gCm2day[n,,]     = site_output$npp_gCm2day
           grid_output$harvest_gCm2day[n,,] = site_output$harvest_gCm2day
           grid_output$fire_gCm2day[n,,]    = site_output$fire_gCm2day
           grid_output$nbe_gCm2day[n,,]     = site_output$nbe_gCm2day
           grid_output$nbp_gCm2day[n,,]     = site_output$nbp_gCm2day
           # Always present but at annual time step
           grid_output$mean_annual_cue[n,,]             = site_output$mean_annual_cue
           grid_output$mean_annual_Ctotal_gCm2[n,,]     = site_output$mean_annual_Ctotal_gCm2
           grid_output$mean_annual_lai_m2m2[n,,]        = site_output$mean_annual_lai_m2m2
           grid_output$mean_annual_nee_gCm2day[n,,]     = site_output$mean_annual_nee_gCm2day
           grid_output$mean_annual_gpp_gCm2day[n,,]     = site_output$mean_annual_gpp_gCm2day
           grid_output$mean_annual_rauto_gCm2day[n,,]   = site_output$mean_annual_rauto_gCm2day
           grid_output$mean_annual_rhet_gCm2day[n,,]    = site_output$mean_annual_rhet_gCm2day
           grid_output$mean_annual_reco_gCm2day[n,,]    = site_output$mean_annual_reco_gCm2day
           grid_output$mean_annual_npp_gCm2day[n,,]     = site_output$mean_annual_npp_gCm2day
           grid_output$mean_annual_harvest_gCm2day[n,,] = site_output$mean_annual_harvest_gCm2day
           grid_output$mean_annual_fire_gCm2day[n,,]    = site_output$mean_annual_fire_gCm2day
           grid_output$mean_annual_nbe_gCm2day[n,,]     = site_output$mean_annual_nbe_gCm2day
           grid_output$mean_annual_nbp_gCm2day[n,,]     = site_output$mean_annual_nbp_gCm2day
           # Always present but at annual time step differences
           grid_output$mean_annual_dlai_m2m2[n,,]        = site_output$mean_annual_dlai_m2m2
           grid_output$mean_annual_dCtotal_gCm2[n,,]     = site_output$mean_annual_dCtotal_gCm2
           grid_output$mean_annual_dnee_gCm2day[n,,]     = site_output$mean_annual_dnee_gCm2day
           grid_output$mean_annual_dgpp_gCm2day[n,,]     = site_output$mean_annual_dgpp_gCm2day
           grid_output$mean_annual_drauto_gCm2day[n,,]   = site_output$mean_annual_drauto_gCm2day
           grid_output$mean_annual_drhet_gCm2day[n,,]    = site_output$mean_annual_drhet_gCm2day
           grid_output$mean_annual_dreco_gCm2day[n,,]    = site_output$mean_annual_dreco_gCm2day
           grid_output$mean_annual_dnpp_gCm2day[n,,]     = site_output$mean_annual_dnpp_gCm2day
           grid_output$mean_annual_dharvest_gCm2day[n,,] = site_output$mean_annual_dharvest_gCm2day
           grid_output$mean_annual_dfire_gCm2day[n,,]    = site_output$mean_annual_dfire_gCm2day
           grid_output$mean_annual_dnbe_gCm2day[n,,]     = site_output$mean_annual_dnbe_gCm2day
           grid_output$mean_annual_dnbp_gCm2day[n,,]     = site_output$mean_annual_dnbp_gCm2day

           if (any(check_list == "grazing_gCm2day")) {
               grid_output$grazing_gCm2day[n,,]                 = site_output$grazing_gCm2day
               grid_output$mean_grazing_gCm2day[slot_i,slot_j,] = site_output$mean_grazing_gCm2day
               grid_output$mean_annual_grazing_gCm2day[n,,]     = site_output$mean_annual_grazing_gCm2day
               grid_output$mean_annual_dgrazing_gCm2day[n,,]    = site_output$mean_annual_dgrazing_gCm2day
           }
           # Extract quantiles
           #grid_output$annual_change_Ctotal_gCm2[n,,] = site_output$annual_change_Ctotal_gCm2
           #grid_output$annual_change_biomass_gCm2[n,,] = site_output$annual_change_biomass_gCm2
           #grid_output$annual_change_dom_gCm2[n,,] = site_output$annual_change_dom_gCm2

           # Based on the presence of each pool define the grids for the mean and final values.
           # Also, create the time varying but quantile based values and time

           # Gridded biomass information
           if (any(check_list == "biomass_gCm2")) {
               # Grid mean / finals for globally available variables
               grid_output$mean_biomass_gCm2[slot_i,slot_j,] = site_output$mean_biomass_gCm2
               grid_output$final_biomass_gCm2[slot_i,slot_j,] = site_output$biomass_gCm2[,grid_output$time_dim]
               grid_output$final_dCbiomass_gCm2[slot_i,slot_j,] = site_output$dCbiomass_gCm2[,grid_output$time_dim]
               grid_output$mean_outflux_biomass_gCm2day[slot_i,slot_j,] = site_output$mean_outflux_biomass_gCm2day
               grid_output$mean_biomass_to_litter_gCm2day[slot_i,slot_j,] = site_output$mean_biomass_to_litter_gCm2day
               grid_output$mean_combined_biomass_to_litter_gCm2day[slot_i,slot_j,] = site_output$mean_combined_biomass_to_litter_gCm2day
               grid_output$MTT_biomass_years[slot_i,slot_j,] = site_output$MTT_biomass_years
               # Time varying pixel specific with quantiles
               grid_output$biomass_gCm2[n,,] = site_output$biomass_gCm2
               grid_output$dCbiomass_gCm2[n,,] = site_output$dCbiomass_gCm2
               grid_output$outflux_biomass_gCm2day[n,,] = site_output$outflux_biomass_gCm2day
               grid_output$biomass_to_litter_gCm2day[n,,] = site_output$biomass_to_litter_gCm2day
               grid_output$combined_biomass_to_litter_gCm2day[n,,] = site_output$combined_biomass_to_litter_gCm2day
               # Annual information
               grid_output$mean_annual_biomass_gCm2[n,,] = site_output$mean_annual_biomass_gCm2
               grid_output$mean_annual_dCbiomass_gCm2[n,,] = site_output$mean_annual_dCbiomass_gCm2
               grid_output$mean_annual_outflux_biomass_gCm2day[n,,] = site_output$mean_annual_outflux_biomass_gCm2day
               grid_output$mean_annual_combined_biomass_to_litter_gCm2day[n,,] = site_output$mean_annual_combined_biomass_to_litter_gCm2day
               grid_output$mean_annual_biomass_to_litter_gCm2day[n,,] = site_output$mean_annual_biomass_to_litter_gCm2day
               grid_output$MTT_annual_biomass_years[n,,] = site_output$MTT_annual_biomass_years
               # Fractional partitioning of tunover to different drivers - should they exist
               grid_output$NaturalFractionOfTurnover_biomass[slot_i,slot_j,] = site_output$NaturalFractionOfTurnover_biomass
               grid_output$FireFractionOfTurnover_biomass[slot_i,slot_j,] = site_output$FireFractionOfTurnover_biomass
               grid_output$HarvestFractionOfTurnover_biomass[slot_i,slot_j,] = site_output$HarvestFractionOfTurnover_biomass
               grid_output$GRAZINGFractionOfTurnover_biomass[slot_i,slot_j,] = site_output$GRAZINGFractionOfTurnover_biomass
               # Conditional variables
               if (any(check_list == "FIREemiss_biomass_gCm2day")) {
                   grid_output$FIREemiss_biomass_gCm2day[n,,] = site_output$FIREemiss_biomass_gCm2day
                   grid_output$mean_FIREemiss_biomass_gCm2day[slot_i,slot_j,] = site_output$mean_FIREemiss_biomass_gCm2day
                   grid_output$mean_annual_FIREemiss_biomass_gCm2day[n,,] = site_output$mean_annual_FIREemiss_biomass_gCm2day
               }
               if (any(check_list == "FIRElitter_biomass_gCm2day")) {
                   grid_output$FIRElitter_biomass_gCm2day[n,,] = site_output$FIRElitter_biomass_gCm2day
                   grid_output$mean_FIRElitter_biomass_gCm2day[slot_i,slot_j,] = site_output$mean_FIRElitter_biomass_gCm2day
                   grid_output$mean_annual_FIRElitter_biomass_gCm2day[n,,] = site_output$mean_annual_FIRElitter_biomass_gCm2day
               }
               if (any(check_list == "HARVESTextracted_biomass_gCm2day")) {
                   grid_output$HARVESTextracted_biomass_gCm2day[n,,] = site_output$HARVESTextracted_biomass_gCm2day
                   grid_output$mean_HARVESTextracted_biomass_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTextracted_biomass_gCm2day
                   grid_output$mean_annual_HARVESTextracted_biomass_gCm2day[n,,] = site_output$mean_annual_HARVESTextracted_biomass_gCm2day
                   grid_output$HARVESTlitter_biomass_gCm2day[n,,] = site_output$HARVESTlitter_biomass_gCm2day
                   grid_output$mean_HARVESTlitter_biomass_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTlitter_biomass_gCm2day
                   grid_output$mean_annual_HARVESTlitter_biomass_gCm2day[n,,] = site_output$mean_annual_HARVESTlitter_biomass_gCm2day
               }
               if (any(check_list == "GRAZINGextracted_biomass_gCm2day")) {
                   grid_output$GRAZINGextracted_biomass_gCm2day[n,,] = site_output$GRAZINGextracted_biomass_gCm2day
                   grid_output$mean_GRAZINGextracted_biomass_gCm2day[slot_i,slot_j,] = site_output$mean_GRAZINGextracted_biomass_gCm2day
                   grid_output$mean_annual_GRAZINGextracted_biomass_gCm2day[n,,] = site_output$mean_annual_GRAZINGextracted_biomass_gCm2day
                   grid_output$GRAZINGlitter_biomass_gCm2day[n,,] = site_output$GRAZINGlitter_biomass_gCm2day
                   grid_output$mean_GRAZINGlitter_biomass_gCm2day[slot_i,slot_j,] = site_output$mean_GRAZINGlitter_biomass_gCm2day
                   grid_output$mean_annual_GRAZINGlitter_biomass_gCm2day[n,,] = site_output$mean_annual_GRAZINGlitter_biomass_gCm2day
               }               
           }
           # Gridded labile information
           if (any(check_list == "labile_gCm2")) {
               # Grid averages
               grid_output$mean_labile_gCm2[slot_i,slot_j,] = site_output$mean_labile_gCm2
               grid_output$final_labile_gCm2[slot_i,slot_j,] = site_output$labile_gCm2[,grid_output$time_dim]
               grid_output$final_dClabile_gCm2[slot_i,slot_j,] = site_output$dClabile_gCm2[,grid_output$time_dim]
               grid_output$mean_outflux_labile_gCm2day[slot_i,slot_j,] = site_output$mean_outflux_labile_gCm2day
               grid_output$mean_labile_to_foliage_gCm2day[slot_i,slot_j,] = site_output$mean_labile_to_foliage_gCm2day
               grid_output$mean_alloc_labile_gCm2day[slot_i,slot_j,] = site_output$mean_alloc_labile_gCm2day
               grid_output$mean_combined_labile_to_litter_gCm2day[slot_i,slot_j,] = site_output$mean_combined_labile_to_litter_gCm2day
               # Time varying
               grid_output$labile_gCm2[n,,] = site_output$labile_gCm2
               grid_output$dClabile_gCm2[n,,] = site_output$dClabile_gCm2
               grid_output$outflux_labile_gCm2day[n,,] = site_output$outflux_labile_gCm2day
               grid_output$labile_to_foliage_gCm2day[n,,] = site_output$labile_to_foliage_gCm2day
               grid_output$alloc_labile_gCm2day[n,,] = site_output$alloc_labile_gCm2day
               grid_output$combined_labile_to_litter_gCm2day[n,,] = site_output$combined_labile_to_litter_gCm2day
               # Annual information
               grid_output$mean_annual_labile_gCm2[n,,] = site_output$mean_annual_labile_gCm2
               grid_output$mean_annual_dClabile_gCm2[n,,] = site_output$mean_annual_dClabile_gCm2
               grid_output$mean_annual_outflux_labile_gCm2day[n,,] = site_output$mean_annual_outflux_labile_gCm2day
               grid_output$mean_annual_labile_to_foliage_gCm2day[n,,] = site_output$mean_annual_labile_to_foliage_gCm2day
               grid_output$mean_annual_alloc_labile_gCm2day[n,,] = site_output$mean_annual_alloc_labile_gCm2day
               grid_output$mean_annual_combined_labile_to_litter_gCm2day[n,,] = site_output$mean_annual_combined_labile_to_litter_gCm2day
               grid_output$MTT_annual_labile_years[n,,] = site_output$MTT_annual_labile_years
               # Fractional partitioning of tunover to different drivers - should they exist
               grid_output$NaturalFractionOfTurnover_labile[slot_i,slot_j,] = site_output$NaturalFractionOfTurnover_labile
               grid_output$FireFractionOfTurnover_labile[slot_i,slot_j,] = site_output$FireFractionOfTurnover_labile
               grid_output$HarvestFractionOfTurnover_labile[slot_i,slot_j,] = site_output$HarvestFractionOfTurnover_labile
               grid_output$GRAZINGFractionOfTurnover_labile[slot_i,slot_j,] = site_output$GRAZINGFractionOfTurnover_labile
               # Conditional variables
               if (any(check_list == "FIREemiss_labile_gCm2day")) {
                   grid_output$FIREemiss_labile_gCm2day[n,,] = site_output$FIREemiss_labile_gCm2day
                   grid_output$mean_FIREemiss_labile_gCm2day[slot_i,slot_j,] = site_output$mean_FIREemiss_labile_gCm2day
                   grid_output$mean_annual_FIREemiss_labile_gCm2day[n,,] = site_output$mean_annual_FIREemiss_labile_gCm2day
               }
               if (any(check_list == "FIRElitter_labile_gCm2day")) {
                   grid_output$FIRElitter_labile_gCm2day[n,,] = site_output$FIRElitter_labile_gCm2day
                   grid_output$mean_FIRElitter_labile_gCm2day[slot_i,slot_j,] = site_output$mean_FIRElitter_labile_gCm2day
                   grid_output$mean_annual_FIRElitter_labile_gCm2day[n,,] = site_output$mean_annual_FIRElitter_labile_gCm2day
               }
               if (any(check_list == "HARVESTextracted_labile_gCm2day")) {
                   grid_output$HARVESTextracted_labile_gCm2day[n,,] = site_output$HARVESTextracted_labile_gCm2day
                   grid_output$mean_HARVESTextracted_labile_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTextracted_labile_gCm2day
                   grid_output$mean_annual_HARVESTextracted_labile_gCm2day[n,,] = site_output$mean_annual_HARVESTextracted_labile_gCm2day
                   grid_output$HARVESTlitter_labile_gCm2day[n,,] = site_output$HARVESTlitter_labile_gCm2day
                   grid_output$mean_HARVESTlitter_labile_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTlitter_labile_gCm2day
                   grid_output$mean_annual_HARVESTlitter_labile_gCm2day[n,,] = site_output$mean_annual_HARVESTlitter_labile_gCm2day
               }
               if (any(check_list == "GRAZINGextracted_labile_gCm2day")) {
                   grid_output$GRAZINGextracted_labile_gCm2day[n,,] = site_output$GRAZINGextracted_labile_gCm2day
                   grid_output$mean_GRAZINGextracted_labile_gCm2day[slot_i,slot_j,] = site_output$mean_GRAZINGextracted_labile_gCm2day
                   grid_output$mean_annual_GRAZINGextracted_labile_gCm2day[n,,] = site_output$mean_annual_GRAZINGextracted_labile_gCm2day
                   grid_output$GRAZINGlitter_labile_gCm2day[n,,] = site_output$GRAZINGlitter_labile_gCm2day
                   grid_output$mean_GRAZINGlitter_labile_gCm2day[slot_i,slot_j,] = site_output$mean_GRAZINGlitter_labile_gCm2day
                   grid_output$mean_annual_GRAZINGlitter_labile_gCm2day[n,,] = site_output$mean_annual_GRAZINGlitter_labile_gCm2day
               }               
           }
           # Gridded foliage information
           if (any(check_list == "foliage_gCm2")) {
               # Grided mean / final
               grid_output$mean_foliage_gCm2[slot_i,slot_j,] = site_output$mean_foliage_gCm2
               grid_output$final_foliage_gCm2[slot_i,slot_j,] = site_output$foliage_gCm2[,grid_output$time_dim]
               grid_output$final_dCfoliage_gCm2[slot_i,slot_j,] = site_output$dCfoliage_gCm2[,grid_output$time_dim]
               grid_output$mean_outflux_foliage_gCm2day[slot_i,slot_j,] = site_output$mean_outflux_foliage_gCm2day
               grid_output$mean_foliage_to_litter_gCm2day[slot_i,slot_j,] = site_output$mean_foliage_to_litter_gCm2day
               grid_output$mean_combined_alloc_foliage_gCm2day[slot_i,slot_j,] = site_output$mean_combined_alloc_foliage_gCm2day
               grid_output$mean_combined_foliage_to_litter_gCm2day[slot_i,slot_j,] = site_output$mean_combined_foliage_to_litter_gCm2day
               # Time varying pixel specific
               grid_output$foliage_gCm2[n,,] = site_output$foliage_gCm2
               grid_output$dCfoliage_gCm2[n,,] = site_output$dCfoliage_gCm2
               grid_output$outflux_foliage_gCm2day[n,,] = site_output$outflux_foliage_gCm2day
               grid_output$foliage_to_litter_gCm2day[n,,] = site_output$foliage_to_litter_gCm2day
               grid_output$combined_alloc_foliage_gCm2day[n,,] = site_output$combined_alloc_foliage_gCm2day
               grid_output$combined_foliage_to_litter_gCm2day[n,,] = site_output$combined_foliage_to_litter_gCm2day
               # Annual information
               grid_output$mean_annual_foliage_gCm2[n,,] = site_output$mean_annual_foliage_gCm2
               grid_output$mean_annual_dCfoliage_gCm2[n,,] = site_output$mean_annual_dCfoliage_gCm2               
               grid_output$mean_annual_outflux_foliage_gCm2day[n,,] = site_output$mean_annual_outflux_foliage_gCm2day
               grid_output$mean_annual_foliage_to_litter_gCm2day[n,,] = site_output$mean_annual_foliage_to_litter_gCm2day
               grid_output$mean_annual_combined_alloc_foliage_gCm2day[n,,] = site_output$mean_annual_combined_alloc_foliage_gCm2day
               grid_output$mean_annual_combined_foliage_to_litter_gCm2day[n,,] = site_output$mean_annual_combined_foliage_to_litter_gCm2day
               grid_output$MTT_annual_foliage_years[n,,] = site_output$MTT_annual_foliage_years
               # Fractional partitioning of tunover to different drivers - should they exist
               grid_output$NaturalFractionOfTurnover_foliage[slot_i,slot_j,] = site_output$NaturalFractionOfTurnover_foliage
               grid_output$FireFractionOfTurnover_foliage[slot_i,slot_j,] = site_output$FireFractionOfTurnover_foliage
               grid_output$HarvestFractionOfTurnover_foliage[slot_i,slot_j,] = site_output$HarvestFractionOfTurnover_foliage
               grid_output$GRAZINGFractionOfTurnover_foliage[slot_i,slot_j,] = site_output$GRAZINGFractionOfTurnover_foliage
               # Conditional variables
               if (any(check_list == "alloc_foliage_gCm2day")) {
                   grid_output$alloc_foliage_gCm2day[n,,] = site_output$alloc_foliage_gCm2day
                   grid_output$mean_alloc_foliage_gCm2day[slot_i,slot_j,] = site_output$mean_alloc_foliage_gCm2day
                   grid_output$mean_annual_alloc_foliage_gCm2day[n,,] = site_output$mean_annual_alloc_foliage_gCm2day
               }
               if (any(check_list == "FIREemiss_foliage_gCm2day")) {
                   grid_output$FIREemiss_foliage_gCm2day[n,,] = site_output$FIREemiss_foliage_gCm2day
                   grid_output$mean_FIREemiss_foliage_gCm2day[slot_i,slot_j,] = site_output$mean_FIREemiss_foliage_gCm2day
                   grid_output$mean_annual_FIREemiss_foliage_gCm2day[n,,] = site_output$mean_annual_FIREemiss_foliage_gCm2day
               }
               if (any(check_list == "FIRElitter_foliage_gCm2day")) {
                   grid_output$FIRElitter_foliage_gCm2day[n,,] = site_output$FIRElitter_foliage_gCm2day
                   grid_output$mean_FIRElitter_foliage_gCm2day[slot_i,slot_j,] = site_output$mean_FIRElitter_foliage_gCm2day
                   grid_output$mean_annual_FIRElitter_foliage_gCm2day[n,,] = site_output$mean_annual_FIRElitter_foliage_gCm2day
               }
               if (any(check_list == "HARVESTextracted_foliage_gCm2day")) {
                   grid_output$HARVESTextracted_foliage_gCm2day[n,,] = site_output$HARVESTextracted_foliage_gCm2day
                   grid_output$mean_HARVESTextracted_foliage_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTextracted_foliage_gCm2day
                   grid_output$mean_annual_HARVESTextracted_foliage_gCm2day[n,,] = site_output$mean_annual_HARVESTextracted_foliage_gCm2day
                   grid_output$HARVESTlitter_foliage_gCm2day[n,,] = site_output$HARVESTlitter_foliage_gCm2day
                   grid_output$mean_HARVESTlitter_foliage_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTlitter_foliage_gCm2day
                   grid_output$mean_annual_HARVESTlitter_foliage_gCm2day[n,,] = site_output$mean_annual_HARVESTlitter_foliage_gCm2day
               }
               if (any(check_list == "GRAZINGextracted_foliage_gCm2day")) {
                   grid_output$GRAZINGextracted_foliage_gCm2day[n,,] = site_output$GRAZINGextracted_foliage_gCm2day
                   grid_output$mean_GRAZINGextracted_foliage_gCm2day[slot_i,slot_j,] = site_output$mean_GRAZINGextracted_foliage_gCm2day
                   grid_output$mean_annual_GRAZINGextracted_foliage_gCm2day[n,,] = site_output$mean_annual_GRAZINGextracted_foliage_gCm2day
                   grid_output$GRAZINGlitter_foliage_gCm2day[n,,] = site_output$GRAZINGlitter_foliage_gCm2day
                   grid_output$mean_GRAZINGlitter_foliage_gCm2day[slot_i,slot_j,] = site_output$mean_GRAZINGlitter_foliage_gCm2day
                   grid_output$mean_annual_GRAZINGlitter_foliage_gCm2day[n,,] = site_output$mean_annual_GRAZINGlitter_foliage_gCm2day
               }               
           }
           #  Gridded roots information
           if (any(check_list == "roots_gCm2")) {
               # Gridded mean / final
               grid_output$mean_roots_gCm2[slot_i,slot_j,] = site_output$mean_roots_gCm2
               grid_output$final_roots_gCm2[slot_i,slot_j,] = site_output$roots_gCm2[,grid_output$time_dim]
               grid_output$final_dCroots_gCm2[slot_i,slot_j,] = site_output$dCroots_gCm2[,grid_output$time_dim]
               grid_output$mean_outflux_roots_gCm2day[slot_i,slot_j,] = site_output$mean_outflux_roots_gCm2day
               grid_output$mean_roots_to_litter_gCm2day[slot_i,slot_j,] = site_output$mean_roots_to_litter_gCm2day
               grid_output$mean_alloc_roots_gCm2day[slot_i,slot_j,] = site_output$mean_alloc_roots_gCm2day
               grid_output$annual_max_roots_gCm2[slot_i,slot_j,] = site_output$annual_max_roots_gCm2
               grid_output$mean_combined_alloc_roots_gCm2day[slot_i,slot_j,] = site_output$mean_combined_alloc_roots_gCm2day                              
               grid_output$mean_combined_roots_to_litter_gCm2day[slot_i,slot_j,] = site_output$mean_combined_roots_to_litter_gCm2day
               # Pixel specific time varying
               grid_output$roots_gCm2[n,,] = site_output$roots_gCm2
               grid_output$dCroots_gCm2[n,,] = site_output$dCroots_gCm2
               grid_output$outflux_roots_gCm2day[n,,] = site_output$outflux_roots_gCm2day
               grid_output$roots_to_litter_gCm2day[n,,] = site_output$roots_to_litter_gCm2day
               grid_output$alloc_roots_gCm2day[n,,] = site_output$alloc_roots_gCm2day
               grid_output$combined_alloc_roots_gCm2day[n,,] = site_output$combined_alloc_roots_gCm2day                              
               grid_output$combined_roots_to_litter_gCm2day[n,,] = site_output$combined_roots_to_litter_gCm2day
               # Annual information
               grid_output$mean_annual_roots_gCm2[n,,] = site_output$mean_annual_roots_gCm2
               grid_output$mean_annual_dCroots_gCm2[n,,] = site_output$mean_annual_dCroots_gCm2
               grid_output$mean_annual_outflux_roots_gCm2day[n,,] = site_output$mean_annual_outflux_roots_gCm2day
               grid_output$mean_annual_roots_to_litter_gCm2day[n,,] = site_output$mean_annual_roots_to_litter_gCm2day
               grid_output$mean_annual_alloc_roots_gCm2day[n,,] = site_output$mean_annual_alloc_roots_gCm2day
               grid_output$mean_annual_combined_alloc_roots_gCm2day[n,,] = site_output$mean_annual_combined_alloc_roots_gCm2day               
               grid_output$mean_annual_combined_roots_to_litter_gCm2day[n,,] = site_output$mean_annual_combined_roots_to_litter_gCm2day
               grid_output$MTT_annual_roots_years[n,,] = site_output$MTT_annual_roots_years
               # Fractional partitioning of tunover to different drivers - should they exist
               grid_output$NaturalFractionOfTurnover_roots[slot_i,slot_j,] = site_output$NaturalFractionOfTurnover_roots
               grid_output$FireFractionOfTurnover_roots[slot_i,slot_j,] = site_output$FireFractionOfTurnover_roots
               grid_output$HarvestFractionOfTurnover_roots[slot_i,slot_j,] = site_output$HarvestFractionOfTurnover_roots
               grid_output$GRAZINGFractionOfTurnover_roots[slot_i,slot_j,] = site_output$GRAZINGFractionOfTurnover_roots
               # Is rooting depth calculate (m) by this model?
               if (any(check_list == "RootDepth_m")) {
                   grid_output$RootDepth_m[n,,] = site_output$RootDepth_m
                   grid_output$dRootDepth_m[n,,] = site_output$dRootDepth_m
                   grid_output$mean_RootDepth_m[slot_i,slot_j,] = site_output$mean_RootDepth_m
                   grid_output$final_dRootDepth_m[slot_i,slot_j,] = site_output$dRootDepth_m[,grid_output$time_dim]
                   grid_output$mean_annual_RootDepth_m[n,,] = site_output$mean_annual_RootDepth_m
               }
               # Conditional variables
               if (any(check_list == "FIREemiss_roots_gCm2day")) {
                   grid_output$FIREemiss_roots_gCm2day[n,,] = site_output$FIREemiss_roots_gCm2day
                   grid_output$mean_FIREemiss_roots_gCm2day[slot_i,slot_j,] = site_output$mean_FIREemiss_roots_gCm2day
                   grid_output$mean_annual_FIREemiss_roots_gCm2day[n,,] = site_output$mean_annual_FIREemiss_roots_gCm2day
               }
               if (any(check_list == "FIRElitter_roots_gCm2day")) {
                   grid_output$FIRElitter_roots_gCm2day[n,,] = site_output$FIRElitter_roots_gCm2day
                   grid_output$mean_FIRElitter_roots_gCm2day[slot_i,slot_j,] = site_output$mean_FIRElitter_roots_gCm2day
                   grid_output$mean_annual_FIRElitter_roots_gCm2day[n,,] = site_output$mean_annual_FIRElitter_roots_gCm2day
               }
               if (any(check_list == "HARVESTextracted_roots_gCm2day")) {
                   grid_output$HARVESTextracted_roots_gCm2day[n,,] = site_output$HARVESTextracted_roots_gCm2day
                   grid_output$mean_HARVESTextracted_roots_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTextracted_roots_gCm2day
                   grid_output$mean_annual_HARVESTextracted_roots_gCm2day[n,,] = site_output$mean_annual_HARVESTextracted_roots_gCm2day
                   grid_output$HARVESTlitter_roots_gCm2day[n,,] = site_output$HARVESTlitter_roots_gCm2day
                   grid_output$mean_HARVESTlitter_roots_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTlitter_roots_gCm2day
                   grid_output$mean_annual_HARVESTlitter_roots_gCm2day[n,,] = site_output$mean_annual_HARVESTlitter_roots_gCm2day
               }
               if (any(check_list == "GRAZINGextracted_roots_gCm2day")) {
                   grid_output$GRAZINGextracted_roots_gCm2day[n,,] = site_output$GRAZINGextracted_roots_gCm2day
                   grid_output$mean_GRAZINGextracted_roots_gCm2day[slot_i,slot_j,] = site_output$mean_GRAZINGextracted_roots_gCm2day
                   grid_output$mean_annual_GRAZINGextracted_roots_gCm2day[n,,] = site_output$mean_annual_GRAZINGextracted_roots_gCm2day
                   grid_output$GRAZINGlitter_roots_gCm2day[n,,] = site_output$GRAZINGlitter_roots_gCm2day
                   grid_output$mean_GRAZINGlitter_roots_gCm2day[slot_i,slot_j,] = site_output$mean_GRAZINGlitter_roots_gCm2day
                   grid_output$mean_annual_GRAZINGlitter_roots_gCm2day[n,,] = site_output$mean_annual_GRAZINGlitter_roots_gCm2day
               }               
           }
           # Gridded wood information
           if (any(check_list == "wood_gCm2")) {
               # Gridded mean / final variables
               grid_output$mean_wood_gCm2[slot_i,slot_j,] = site_output$mean_wood_gCm2
               grid_output$final_wood_gCm2[slot_i,slot_j,] = site_output$wood_gCm2[,grid_output$time_dim]
               grid_output$final_dCwood_gCm2[slot_i,slot_j,] = site_output$dCwood_gCm2[,grid_output$time_dim]
               grid_output$mean_outflux_wood_gCm2day[slot_i,slot_j,] = site_output$mean_outflux_wood_gCm2day
               grid_output$mean_wood_to_litter_gCm2day[slot_i,slot_j,] = site_output$mean_wood_to_litter_gCm2day
               grid_output$mean_alloc_wood_gCm2day[slot_i,slot_j,] = site_output$mean_alloc_wood_gCm2day
               grid_output$annual_max_wood_gCm2[slot_i,slot_j,] = site_output$annual_max_wood_gCm2
               grid_output$mean_combined_alloc_wood_gCm2day[slot_i,slot_j,] = site_output$mean_combined_alloc_wood_gCm2day
               grid_output$mean_combined_wood_to_litter_gCm2day[slot_i,slot_j,] = site_output$mean_combined_wood_to_litter_gCm2day
               # Pixel specific time varying
               grid_output$wood_gCm2[n,,] = site_output$wood_gCm2
               grid_output$dCwood_gCm2[n,,] = site_output$dCwood_gCm2
               grid_output$outflux_wood_gCm2day[n,,] = site_output$outflux_wood_gCm2day
               grid_output$wood_to_litter_gCm2day[n,,] = site_output$wood_to_litter_gCm2day
               grid_output$alloc_wood_gCm2day[n,,] = site_output$alloc_wood_gCm2day
               grid_output$combined_alloc_wood_gCm2day[n,,] = site_output$combined_alloc_wood_gCm2day               
               grid_output$combined_wood_to_litter_gCm2day[n,,] = site_output$combined_wood_to_litter_gCm2day
               # Annual information
               grid_output$mean_annual_wood_gCm2[n,,] = site_output$mean_annual_wood_gCm2
               grid_output$mean_annual_dCwood_gCm2[n,,] = site_output$mean_annual_dCwood_gCm2
               grid_output$mean_annual_outflux_wood_gCm2day[n,,] = site_output$mean_annual_outflux_wood_gCm2day
               grid_output$mean_annual_wood_to_litter_gCm2day[n,,] = site_output$mean_annual_wood_to_litter_gCm2day
               grid_output$mean_annual_alloc_wood_gCm2day[n,,] = site_output$mean_annual_alloc_wood_gCm2day
               grid_output$mean_annual_combined_alloc_wood_gCm2day[n,,] = site_output$mean_annual_combined_alloc_wood_gCm2day               
               grid_output$mean_annual_combined_wood_to_litter_gCm2day[n,,] = site_output$mean_annual_combined_wood_to_litter_gCm2day
               grid_output$MTT_annual_wood_years[n,,] = site_output$MTT_annual_wood_years
               # Fractional partitioning of tunover to different drivers - should they exist
               grid_output$NaturalFractionOfTurnover_wood[slot_i,slot_j,] = site_output$NaturalFractionOfTurnover_wood
               grid_output$FireFractionOfTurnover_wood[slot_i,slot_j,] = site_output$FireFractionOfTurnover_wood
               grid_output$HarvestFractionOfTurnover_wood[slot_i,slot_j,] = site_output$HarvestFractionOfTurnover_wood
               grid_output$GRAZINGFractionOfTurnover_wood[slot_i,slot_j,] = site_output$GRAZINGFractionOfTurnover_wood
               # Conditional variables
               if (any(check_list == "FIREemiss_wood_gCm2day")) {
                   grid_output$FIREemiss_wood_gCm2day[n,,] = site_output$FIREemiss_wood_gCm2day
                   grid_output$mean_FIREemiss_wood_gCm2day[slot_i,slot_j,] = site_output$mean_FIREemiss_wood_gCm2day
                   grid_output$mean_annual_FIREemiss_wood_gCm2day[n,,] = site_output$mean_annual_FIREemiss_wood_gCm2day
               }
               if (any(check_list == "FIRElitter_wood_gCm2day")) {
                   grid_output$FIRElitter_wood_gCm2day[n,,] = site_output$FIRElitter_wood_gCm2day
                   grid_output$mean_FIRElitter_wood_gCm2day[slot_i,slot_j,] = site_output$mean_FIRElitter_wood_gCm2day
                   grid_output$mean_annual_FIRElitter_wood_gCm2day[n,,] = site_output$mean_annual_FIRElitter_wood_gCm2day
               }
               if (any(check_list == "HARVESTextracted_wood_gCm2day")) {
                   grid_output$HARVESTextracted_wood_gCm2day[n,,] = site_output$HARVESTextracted_wood_gCm2day
                   grid_output$mean_HARVESTextracted_wood_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTextracted_wood_gCm2day
                   grid_output$mean_annual_HARVESTextracted_wood_gCm2day[n,,] = site_output$mean_annual_HARVESTextracted_wood_gCm2day
                   grid_output$HARVESTlitter_wood_gCm2day[n,,] = site_output$HARVESTlitter_wood_gCm2day
                   grid_output$mean_HARVESTlitter_wood_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTlitter_wood_gCm2day
                   grid_output$mean_annual_HARVESTlitter_wood_gCm2day[n,,] = site_output$mean_annual_HARVESTlitter_wood_gCm2day
               }
               if (any(check_list == "GRAZINGextracted_wood_gCm2day")) {
                   grid_output$GRAZINGextracted_wood_gCm2day[n,,] = site_output$GRAZINGextracted_wood_gCm2day
                   grid_output$mean_GRAZINGextracted_wood_gCm2day[slot_i,slot_j,] = site_output$mean_GRAZINGextracted_wood_gCm2day
                   grid_output$mean_annual_GRAZINGextracted_wood_gCm2day[n,,] = site_output$mean_annual_GRAZINGextracted_wood_gCm2day
                   grid_output$GRAZINGlitter_wood_gCm2day[n,,] = site_output$GRAZINGlitter_wood_gCm2day
                   grid_output$mean_GRAZINGlitter_wood_gCm2day[slot_i,slot_j,] = site_output$mean_GRAZINGlitter_wood_gCm2day
                   grid_output$mean_annual_GRAZINGlitter_wood_gCm2day[n,,] = site_output$mean_annual_GRAZINGlitter_wood_gCm2day
               }               
           }
           # Gridded litter information
           if (any(check_list == "litter_gCm2")) {
               # Gridded mean / final
               grid_output$mean_litter_gCm2[slot_i,slot_j,] = site_output$mean_litter_gCm2
               grid_output$final_litter_gCm2[slot_i,slot_j,] = site_output$litter_gCm2[,grid_output$time_dim]
               grid_output$final_dClitter_gCm2[slot_i,slot_j,] = site_output$dClitter_gCm2[,grid_output$time_dim]
               grid_output$mean_outflux_litter_gCm2day[slot_i,slot_j,] = site_output$mean_outflux_litter_gCm2day
               grid_output$mean_litter_to_som_gCm2day[slot_i,slot_j,] = site_output$mean_litter_to_som_gCm2day
               grid_output$mean_rhet_litter_gCm2day[slot_i,slot_j,] = site_output$mean_rhet_litter_gCm2day
               grid_output$mean_combined_litter_to_som_gCm2day[slot_i,slot_j,] = site_output$mean_combined_litter_to_som_gCm2day
               # Pixel specific time varying
               grid_output$litter_gCm2[n,,] = site_output$litter_gCm2
               grid_output$dClitter_gCm2[n,,] = site_output$dClitter_gCm2
               grid_output$outflux_litter_gCm2day[n,,] = site_output$outflux_litter_gCm2day
               grid_output$litter_to_som_gCm2day[n,,] = site_output$litter_to_som_gCm2day
               grid_output$rhet_litter_gCm2day[n,,] = site_output$rhet_litter_gCm2day
               grid_output$combined_litter_to_som_gCm2day[n,,] = site_output$combined_litter_to_som_gCm2day
               # Annual information
               grid_output$mean_annual_litter_gCm2[n,,] = site_output$mean_annual_litter_gCm2
               grid_output$mean_annual_dClitter_gCm2[n,,] = site_output$mean_annual_dClitter_gCm2
               grid_output$mean_annual_outflux_litter_gCm2day[n,,] = site_output$mean_annual_outflux_litter_gCm2day
               grid_output$mean_annual_litter_to_som_gCm2day[n,,] = site_output$mean_annual_litter_to_som_gCm2day
               grid_output$mean_annual_rhet_litter_gCm2day[n,,] = site_output$mean_annual_rhet_litter_gCm2day
               grid_output$mean_annual_drhet_litter_gCm2day[n,,] = site_output$mean_annual_drhet_litter_gCm2day
               grid_output$mean_annual_combined_litter_to_som_gCm2day[n,,] = site_output$mean_annual_combined_litter_to_som_gCm2day
               grid_output$MTT_annual_litter_years[n,,] = site_output$MTT_annual_litter_years
               # Fractional partitioning of tunover to different drivers - should they exist
               grid_output$NaturalFractionOfTurnover_litter[slot_i,slot_j,] = site_output$NaturalFractionOfTurnover_litter
               grid_output$FireFractionOfTurnover_litter[slot_i,slot_j,] = site_output$FireFractionOfTurnover_litter
               grid_output$HarvestFractionOfTurnover_litter[slot_i,slot_j,] = site_output$HarvestFractionOfTurnover_litter
               grid_output$GRAZINGFractionOfTurnover_litter[slot_i,slot_j,] = site_output$GRAZINGFractionOfTurnover_litter
               # Conditional variables
               if (any(check_list == "FIREemiss_litter_gCm2day")) {
                   grid_output$FIREemiss_litter_gCm2day[n,,] = site_output$FIREemiss_litter_gCm2day
                   grid_output$mean_FIREemiss_litter_gCm2day[slot_i,slot_j,] = site_output$mean_FIREemiss_litter_gCm2day
                   grid_output$mean_annual_FIREemiss_litter_gCm2day[n,,] = site_output$mean_annual_FIREemiss_litter_gCm2day
               }
               if (any(check_list == "FIRElitter_litter_gCm2day")) {
                   grid_output$FIRElitter_litter_gCm2day[n,,] = site_output$FIRElitter_litter_gCm2day
                   grid_output$mean_FIRElitter_litter_gCm2day[slot_i,slot_j,] = site_output$mean_FIRElitter_litter_gCm2day
                   grid_output$mean_annual_FIRElitter_litter_gCm2day[n,,] = site_output$mean_annual_FIRElitter_litter_gCm2day
               }
               if (any(check_list == "HARVESTextracted_litter_gCm2day")) {
                   grid_output$HARVESTextracted_litter_gCm2day[n,,] = site_output$HARVESTextracted_litter_gCm2day
                   grid_output$mean_HARVESTextracted_litter_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTextracted_litter_gCm2day
                   grid_output$mean_annual_HARVESTextracted_litter_gCm2day[n,,] = site_output$mean_annual_HARVESTextracted_litter_gCm2day
               }
               if (any(check_list == "GRAZINGextracted_litter_gCm2day")) {
                   grid_output$GRAZINGextracted_litter_gCm2day[n,,] = site_output$GRAZINGextracted_litter_gCm2day
                   grid_output$mean_GRAZINGextracted_litter_gCm2day[slot_i,slot_j,] = site_output$mean_GRAZINGextracted_litter_gCm2day
                   grid_output$mean_annual_GRAZINGextracted_litter_gCm2day[n,,] = site_output$mean_annual_GRAZINGextracted_litter_gCm2day
               }               
           }
           # Gridded wood litter information
           if (any(check_list == "woodlitter_gCm2")) {
               # Gridded mean / final
               grid_output$mean_woodlitter_gCm2[slot_i,slot_j,] = site_output$mean_woodlitter_gCm2
               grid_output$final_woodlitter_gCm2[slot_i,slot_j,] = site_output$woodlitter_gCm2[,grid_output$time_dim]
               grid_output$final_dCwoodlitter_gCm2[slot_i,slot_j,] = site_output$dCwoodlitter_gCm2[,grid_output$time_dim]
               grid_output$mean_outflux_woodlitter_gCm2day[slot_i,slot_j,] = site_output$mean_outflux_woodlitter_gCm2day
               grid_output$mean_woodlitter_to_som_gCm2day[slot_i,slot_j,] = site_output$mean_woodlitter_to_som_gCm2day
               grid_output$mean_rhet_woodlitter_gCm2day[slot_i,slot_j,] = site_output$mean_rhet_woodlitter_gCm2day
               grid_output$mean_combined_woodlitter_to_som_gCm2day[slot_i,slot_j,] = site_output$mean_combined_woodlitter_to_som_gCm2day
               # Pixel specific time varying
               grid_output$woodlitter_gCm2[n,,] = site_output$woodlitter_gCm2
               grid_output$dCwoodlitter_gCm2[n,,] = site_output$dCwoodlitter_gCm2
               grid_output$outflux_woodlitter_gCm2day[n,,] = site_output$outflux_woodlitter_gCm2day
               grid_output$woodlitter_to_som_gCm2day[n,,] = site_output$woodlitter_to_som_gCm2day
               grid_output$rhet_woodlitter_gCm2day[n,,] = site_output$rhet_woodlitter_gCm2day
               grid_output$combined_woodlitter_to_som_gCm2day[n,,] = site_output$combined_woodlitter_to_som_gCm2day
               # Annual information
               grid_output$mean_annual_woodlitter_gCm2[n,,] = site_output$mean_annual_woodlitter_gCm2
               grid_output$mean_annual_dCwoodlitter_gCm2[n,,] = site_output$mean_annual_dCwoodlitter_gCm2
               grid_output$mean_annual_outflux_woodlitter_gCm2day[n,,] = site_output$mean_annual_outflux_woodlitter_gCm2day
               grid_output$mean_annual_woodlitter_to_som_gCm2day[n,,] = site_output$mean_annual_woodlitter_to_som_gCm2day
               grid_output$mean_annual_rhet_woodlitter_gCm2day[n,,] = site_output$mean_annual_rhet_woodlitter_gCm2day
               grid_output$mean_annual_combined_woodlitter_to_som_gCm2day[n,,] = site_output$mean_annual_combined_woodlitter_to_som_gCm2day
               grid_output$MTT_annual_woodlitter_years[n,,] = site_output$MTT_annual_woodlitter_years
               # Fractional partitioning of tunover to different drivers - should they exist
               grid_output$NaturalFractionOfTurnover_woodlitter[slot_i,slot_j,] = site_output$NaturalFractionOfTurnover_woodlitter
               grid_output$FireFractionOfTurnover_woodlitter[slot_i,slot_j,] = site_output$FireFractionOfTurnover_woodlitter
               grid_output$HarvestFractionOfTurnover_woodlitter[slot_i,slot_j,] = site_output$HarvestFractionOfTurnover_woodlitter
               grid_output$GRAZINGFractionOfTurnover_woodlitter[slot_i,slot_j,] = site_output$GRAZINGFractionOfTurnover_woodlitter
               # Conditional variables
               if (any(check_list == "FIREemiss_woodlitter_gCm2day")) {
                   grid_output$FIREemiss_woodlitter_gCm2day[n,,] = site_output$FIREemiss_woodlitter_gCm2day
                   grid_output$mean_FIREemiss_woodlitter_gCm2day[slot_i,slot_j,] = site_output$mean_FIREemiss_woodlitter_gCm2day
                   grid_output$mean_annual_FIREemiss_woodlitter_gCm2day[n,,] = site_output$mean_annual_FIREemiss_woodlitter_gCm2day
               }
               if (any(check_list == "FIRElitter_woodlitter_gCm2day")) {
                   grid_output$FIRElitter_woodlitter_gCm2day[n,,] = site_output$FIRElitter_woodlitter_gCm2day
                   grid_output$mean_FIRElitter_woodlitter_gCm2day[slot_i,slot_j,] = site_output$mean_FIRElitter_woodlitter_gCm2day
                   grid_output$mean_annual_FIRElitter_woodlitter_gCm2day[n,,] = site_output$mean_annual_FIRElitter_woodlitter_gCm2day
               }
               if (any(check_list == "HARVESTextracted_woodlitter_gCm2day")) {
                   grid_output$HARVESTextracted_woodlitter_gCm2day[n,,] = site_output$HARVESTextracted_woodlitter_gCm2day
                   grid_output$mean_HARVESTextracted_woodlitter_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTextracted_woodlitter_gCm2day
                   grid_output$mean_annual_HARVESTextracted_woodlitter_gCm2day[n,,] = site_output$mean_annual_HARVESTextracted_woodlitter_gCm2day
               }
               if (any(check_list == "GRAZINGextracted_woodlitter_gCm2day")) {
                   grid_output$GRAZINGextracted_woodlitter_gCm2day[n,,] = site_output$GRAZINGextracted_woodlitter_gCm2day
                   grid_output$mean_GRAZINGextracted_woodlitter_gCm2day[slot_i,slot_j,] = site_output$mean_GRAZINGextracted_woodlitter_gCm2day
                   grid_output$mean_annual_GRAZINGextracted_woodlitter_gCm2day[n,,] = site_output$mean_annual_GRAZINGextracted_woodlitter_gCm2day
               }               
           }
           # Gridded som information
           if (any(check_list == "som_gCm2")) {
               # Gridded mean / finals
               grid_output$mean_som_gCm2[slot_i,slot_j,] = site_output$mean_som_gCm2
               grid_output$final_som_gCm2[slot_i,slot_j,] = site_output$som_gCm2[,grid_output$time_dim]
               grid_output$final_dCsom_gCm2[slot_i,slot_j,] = site_output$dCsom_gCm2[,grid_output$time_dim]
               grid_output$mean_outflux_som_gCm2day[slot_i,slot_j,] = site_output$mean_outflux_som_gCm2day
               grid_output$mean_rhet_som_gCm2day[slot_i,slot_j,] = site_output$mean_rhet_som_gCm2day
               # Pixel specific time varying
               grid_output$som_gCm2[n,,] = site_output$som_gCm2
               grid_output$dCsom_gCm2[n,,] = site_output$dCsom_gCm2
               grid_output$outflux_som_gCm2day[n,,] = site_output$outflux_som_gCm2day
               grid_output$rhet_som_gCm2day[n,,] = site_output$rhet_som_gCm2day
               # Annual information
               grid_output$mean_annual_som_gCm2[n,,] = site_output$mean_annual_som_gCm2
               grid_output$mean_annual_dCsom_gCm2[n,,] = site_output$mean_annual_dCsom_gCm2
               grid_output$mean_annual_outflux_som_gCm2day[n,,] = site_output$mean_annual_outflux_som_gCm2day
               grid_output$mean_annual_rhet_som_gCm2day[n,,] = site_output$mean_annual_rhet_som_gCm2day
               grid_output$mean_annual_drhet_som_gCm2day[n,,] = site_output$mean_annual_drhet_som_gCm2day               
               grid_output$MTT_annual_som_years[n,,] = site_output$MTT_annual_som_years
               # Fractional partitioning of tunover to different drivers - should they exist
               grid_output$NaturalFractionOfTurnover_som[slot_i,slot_j,] = site_output$NaturalFractionOfTurnover_som
               grid_output$FireFractionOfTurnover_som[slot_i,slot_j,] = site_output$FireFractionOfTurnover_som
               grid_output$HarvestFractionOfTurnover_som[slot_i,slot_j,] = site_output$HarvestFractionOfTurnover_som
               grid_output$GRAZINGFractionOfTurnover_som[slot_i,slot_j,] = site_output$GRAZINGFractionOfTurnover_som
               # Time varying
               if (any(check_list == "FIREemiss_som_gCm2day")) {
                   grid_output$FIREemiss_som_gCm2day[n,,] = site_output$FIREemiss_som_gCm2day
                   grid_output$mean_FIREemiss_som_gCm2day[slot_i,slot_j,] = site_output$mean_FIREemiss_som_gCm2day
                   grid_output$mean_annual_FIREemiss_som_gCm2day[n,,] = site_output$mean_annual_FIREemiss_som_gCm2day
               }
               if (any(check_list == "HARVESTextracted_som_gCm2day")) {
                   grid_output$HARVESTextracted_som_gCm2day[n,,] = site_output$HARVESTextracted_som_gCm2day
                   grid_output$mean_HARVESTextracted_som_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTextracted_som_gCm2day
                   grid_output$mean_annual_HARVESTextracted_som_gCm2day[n,,] = site_output$mean_annual_HARVESTextracted_som_gCm2day
               }
               if (any(check_list == "GRAZINGextracted_som_gCm2day")) {
                   grid_output$GRAZINGextracted_som_gCm2day[n,,] = site_output$GRAZINGextracted_som_gCm2day
                   grid_output$mean_GRAZINGextracted_som_gCm2day[slot_i,slot_j,] = site_output$mean_GRAZINGextracted_som_gCm2day
                   grid_output$mean_annual_GRAZINGextracted_som_gCm2day[n,,] = site_output$mean_annual_GRAZINGextracted_som_gCm2day
               }               
           }
           # Gridded dom information
           if (any(check_list == "dom_gCm2")) {
               # Gridded mean / finals
               grid_output$mean_dom_gCm2[slot_i,slot_j,] = site_output$mean_dom_gCm2
               grid_output$final_dom_gCm2[slot_i,slot_j,] = site_output$dom_gCm2[,grid_output$time_dim]
               grid_output$final_dCdom_gCm2[slot_i,slot_j,] = site_output$dCdom_gCm2[,grid_output$time_dim]
               grid_output$mean_outflux_dom_gCm2day[slot_i,slot_j,] = site_output$mean_outflux_dom_gCm2day
               grid_output$mean_rhet_dom_gCm2day[slot_i,slot_j,] = site_output$mean_rhet_dom_gCm2day
               grid_output$MTT_dom_years[slot_i,slot_j,] = site_output$MTT_dom_years               
               # Pixel specific time varying
               grid_output$dom_gCm2[n,,] = site_output$dom_gCm2
               grid_output$dCdom_gCm2[n,,] = site_output$dCdom_gCm2
               grid_output$outflux_dom_gCm2day[n,,] = site_output$outflux_dom_gCm2day
               grid_output$rhet_dom_gCm2day[n,,] = site_output$rhet_dom_gCm2day
               # Annual information
               grid_output$mean_annual_dom_gCm2[n,,] = site_output$mean_annual_dom_gCm2
               grid_output$mean_annual_dCdom_gCm2[n,,] = site_output$mean_annual_dCdom_gCm2
               grid_output$mean_annual_outflux_dom_gCm2day[n,,] = site_output$mean_annual_outflux_dom_gCm2day
               grid_output$mean_annual_rhet_dom_gCm2day[n,,] = site_output$mean_annual_rhet_dom_gCm2day
               grid_output$MTT_annual_dom_years[n,,] = site_output$MTT_annual_dom_years
               # Fractional partitioning of tunover to different drivers - should they exist
               grid_output$NaturalFractionOfTurnover_dom[slot_i,slot_j,] = site_output$NaturalFractionOfTurnover_dom
               grid_output$FireFractionOfTurnover_dom[slot_i,slot_j,] = site_output$FireFractionOfTurnover_dom
               grid_output$HarvestFractionOfTurnover_dom[slot_i,slot_j,] = site_output$HarvestFractionOfTurnover_dom
               grid_output$GRAZINGFractionOfTurnover_dom[slot_i,slot_j,] = site_output$GRAZINGFractionOfTurnover_dom
               # Conditional
               if (any(check_list == "FIREemiss_dom_gCm2day")) {
                   grid_output$FIREemiss_dom_gCm2day[n,,] = site_output$FIREemiss_dom_gCm2day
                   grid_output$mean_FIREemiss_dom_gCm2day[slot_i,slot_j,] = site_output$mean_FIREemiss_dom_gCm2day
                   grid_output$mean_annual_FIREemiss_dom_gCm2day[n,,] = site_output$mean_annual_FIREemiss_dom_gCm2day
               }
               if (any(check_list == "HARVESTextracted_dom_gCm2day")) {
                   grid_output$HARVESTextracted_dom_gCm2day[n,,] = site_output$HARVESTextracted_dom_gCm2day
                   grid_output$mean_HARVESTextracted_dom_gCm2day[slot_i,slot_j,] = site_output$mean_HARVESTextracted_dom_gCm2day
                   grid_output$mean_annual_HARVESTextracted_dom_gCm2day[n,,] = site_output$mean_annual_HARVESTextracted_dom_gCm2day
               }
               if (any(check_list == "GRAZINGextracted_dom_gCm2day")) {
                   grid_output$GRAZINGextracted_dom_gCm2day[n,,] = site_output$GRAZINGextracted_dom_gCm2day
                   grid_output$mean_GRAZINGextracted_dom_gCm2day[slot_i,slot_j,] = site_output$mean_GRAZINGextracted_dom_gCm2day
                   grid_output$mean_annual_GRAZINGextracted_dom_gCm2day[n,,] = site_output$mean_annual_GRAZINGextracted_dom_gCm2day
               }               
           }

           # Water cycle specific variables
           if (any(check_list == "ET_kgH2Om2day")) {
               # currently water in the soil surface layer (0-30 cm)
               grid_output$mean_SurfWater_kgH2Om2[slot_i,slot_j,] = site_output$mean_SurfWater_kgH2Om2
               grid_output$mean_annual_SurfWater_kgH2Om2[n,,] = site_output$mean_annual_SurfWater_kgH2Om2
               grid_output$mean_annual_dSurfWater_kgH2Om2[n,,] = site_output$mean_annual_dSurfWater_kgH2Om2
               grid_output$final_SurfWater_kgH2Om2[slot_i,slot_j,] = site_output$SurfWater_kgH2Om2[,grid_output$time_dim]
               grid_output$final_dSurfWater_kgH2Om2[slot_i,slot_j,] = site_output$dSurfWater_kgH2Om2[,grid_output$time_dim]
               grid_output$SurfWater_kgH2Om2[n,,] = site_output$SurfWater_kgH2Om2
               grid_output$dSurfWater_kgH2Om2[n,,] = site_output$dSurfWater_kgH2Om2
               # plant apparent soil water potential (MPa)
               grid_output$mean_wSWP_MPa[slot_i,slot_j,] = site_output$mean_wSWP_MPa
               grid_output$mean_annual_wSWP_MPa[n,,] = site_output$mean_annual_wSWP_MPa
               grid_output$mean_annual_dwSWP_MPa[n,,] = site_output$mean_annual_dwSWP_MPa
               grid_output$final_wSWP_MPa[slot_i,slot_j,] = site_output$wSWP_MPa[,grid_output$time_dim]
               grid_output$wSWP_MPa[n,,] = site_output$wSWP_MPa
               grid_output$dwSWP_MPa[n,,] = site_output$dwSWP_MPa
               # evapotranspiration (Etrans + Esoil + Ewetcanopy)
               grid_output$mean_annual_dET_kgH2Om2day[n,,] = site_output$mean_annual_dET_kgH2Om2day               
               grid_output$mean_annual_ET_kgH2Om2day[n,,] = site_output$mean_annual_ET_kgH2Om2day
               grid_output$mean_ET_kgH2Om2day[slot_i,slot_j,] = site_output$mean_ET_kgH2Om2day
               grid_output$ET_kgH2Om2day[n,,] = site_output$ET_kgH2Om2day
               # Ecosystem water use efficiency (GPP/ET)
               grid_output$mean_annual_wue_eco_gCkgH2O[n,,] = site_output$mean_annual_wue_eco_gCkgH2O
               grid_output$mean_wue_eco_gCkgH2O[slot_i,slot_j,] = site_output$mean_wue_eco_gCkgH2O
               grid_output$wue_eco_gCkgH2O[n,,] = site_output$wue_eco_gCkgH2O
               # Check whether the evaporation components exist
               if (any(check_list == "Etrans_kgH2Om2day")) {
                   # Transpiration
                   grid_output$mean_annual_dEtrans_kgH2Om2day[n,,] = site_output$mean_annual_dEtrans_kgH2Om2day               
                   grid_output$mean_annual_Etrans_kgH2Om2day[n,,] = site_output$mean_annual_Etrans_kgH2Om2day
                   grid_output$mean_Etrans_kgH2Om2day[slot_i,slot_j,] = site_output$mean_Etrans_kgH2Om2day
                   grid_output$Etrans_kgH2Om2day[n,,] = site_output$Etrans_kgH2Om2day
                   # Plant water use efficiency (GPP/Etrans)
                   grid_output$mean_annual_wue_plant_gCkgH2O[n,,] = site_output$mean_annual_wue_plant_gCkgH2O
                   grid_output$mean_wue_plant_gCkgH2O[slot_i,slot_j,] = site_output$mean_wue_plant_gCkgH2O
                   grid_output$wue_plant_gCkgH2O[n,,] = site_output$wue_plant_gCkgH2O
               }
               if (any(check_list == "Esoil_kgH2Om2day")) {
                   # Soil evaporation
                   grid_output$mean_annual_dEsoil_kgH2Om2day[n,,] = site_output$mean_annual_dEsoil_kgH2Om2day               
                   grid_output$mean_annual_Esoil_kgH2Om2day[n,,] = site_output$mean_annual_Esoil_kgH2Om2day
                   grid_output$mean_Esoil_kgH2Om2day[slot_i,slot_j,] = site_output$mean_Esoil_kgH2Om2day
                   grid_output$Esoil_kgH2Om2day[n,,] = site_output$Esoil_kgH2Om2day
               }
               if (any(check_list == "Ewetcanopy_kgH2Om2day")) {
                   # Wet canopy evaporation
                   grid_output$mean_annual_dEwetcanopy_kgH2Om2day[n,,] = site_output$mean_annual_dEwetcanopy_kgH2Om2day               
                   grid_output$mean_annual_Ewetcanopy_kgH2Om2day[n,,] = site_output$mean_annual_Ewetcanopy_kgH2Om2day
                   grid_output$mean_Ewetcanopy_kgH2Om2day[slot_i,slot_j,] = site_output$mean_Ewetcanopy_kgH2Om2day
                   grid_output$Ewetcanopy_kgH2Om2day[n,,] = site_output$Ewetcanopy_kgH2Om2day
               }
               if (any(check_list == "runoff_kgH2Om2day")) {
                   # Surface water 
                   grid_output$mean_annual_drunoff_kgH2Om2day[n,,] = site_output$mean_annual_drunoff_kgH2Om2day               
                   grid_output$mean_annual_runoff_kgH2Om2day[n,,] = site_output$mean_annual_runoff_kgH2Om2day
                   grid_output$mean_runoff_kgH2Om2day[slot_i,slot_j,] = site_output$mean_runoff_kgH2Om2day
                   grid_output$runoff_kgH2Om2day[n,,] = site_output$runoff_kgH2Om2day
               }
               if (any(check_list == "underflow_kgH2Om2day")) {
                   # Underflow from bottom of soil column
                   grid_output$mean_annual_dunderflow_kgH2Om2day[n,,] = site_output$mean_annual_dunderflow_kgH2Om2day               
                   grid_output$mean_annual_underflow_kgH2Om2day[n,,] = site_output$mean_annual_underflow_kgH2Om2day
                   grid_output$mean_underflow_kgH2Om2day[slot_i,slot_j,] = site_output$mean_underflow_kgH2Om2day
                   grid_output$underflow_kgH2Om2day[n,,] = site_output$underflow_kgH2Om2day
               }
               if (any(check_list == "total_drainage_kgH2Om2day")) {
                   # Total drainage from soil surface and bottom of soil 
                   grid_output$mean_annual_dtotal_drainage_kgH2Om2day[n,,] = site_output$mean_annual_dtotal_drainage_kgH2Om2day               
                   grid_output$mean_annual_total_drainage_kgH2Om2day[n,,] = site_output$mean_annual_total_drainage_kgH2Om2day
                   grid_output$mean_total_drainage_kgH2Om2day[slot_i,slot_j,] = site_output$mean_total_drainage_kgH2Om2day
                   grid_output$total_drainage_kgH2Om2day[n,,] = site_output$total_drainage_kgH2Om2day
               }
               if (any(check_list == "LWP_MPa")){ 
                   # Leaf water potential
                   grid_output$mean_LWP_MPa[slot_i,slot_j,] = site_output$mean_LWP_MPa
                   grid_output$LWP_MPa[n,,] = site_output$LWP_MPa
               }
           } # ET_kgH2Om2day exists
           # Snow specific
           if (any(check_list == "snow_kgH2Om2")) {
               ## snow on soil surface
               grid_output$mean_annual_snow_kgH2Om2[n,,] = site_output$mean_annual_snow_kgH2Om2
               grid_output$mean_annual_dsnow_kgH2Om2[n,,] = site_output$mean_annual_dsnow_kgH2Om2
               grid_output$mean_snow_kgH2Om2[slot_i,slot_j,] = site_output$mean_snow_kgH2Om2
               grid_output$snow_kgH2Om2[n,,] = site_output$snow_kgH2Om2
           }
           # Canopy process variables
           if (any(check_list == "APAR_MJm2day")) {
               # Absorbed photosynthetically active radation
               grid_output$mean_annual_APAR_MJm2day[n,,] = site_output$mean_annual_APAR_MJm2day
               grid_output$mean_annual_dAPAR_MJm2day[n,,] = site_output$mean_annual_dAPAR_MJm2day
               grid_output$mean_APAR_MJm2day[slot_i,slot_j,] = site_output$mean_APAR_MJm2day
               grid_output$APAR_MJm2day[n,,] = site_output$APAR_MJm2day
           }
           if (any(check_list == "CiCa")) {
               # Canopy Ci:Ca
               grid_output$mean_annual_CiCa[n,,] = site_output$mean_annual_CiCa
               grid_output$mean_annual_dCiCa[n,,] = site_output$mean_annual_dCiCa
               grid_output$mean_CiCa[slot_i,slot_j,] = site_output$mean_CiCa
               grid_output$CiCa[n,,] = site_output$CiCa
           }
           if (any(check_list == "gs_demand_supply_ratio")) {
               # Ratio of stomatal conductance relative to its maximum value,
               # this metric provides information on the demand vs supply constrains on stomatal conductance
               grid_output$mean_annual_gs_demand_supply_ratio[n,,] = site_output$mean_annual_gs_demand_supply_ratio
               grid_output$mean_annual_dgs_demand_supply_ratio[n,,] = site_output$mean_annual_dgs_demand_supply_ratio
               grid_output$mean_gs_demand_supply_ratio[slot_i,slot_j,] = site_output$mean_gs_demand_supply_ratio
               grid_output$gs_demand_supply_ratio[n,,] = site_output$gs_demand_supply_ratio
           }
           if (any(check_list == "gs_mmolH2Om2s")) {
               # Canopy stomatal conductance
               grid_output$mean_annual_gs_mmolH2Om2s[n,,] = site_output$mean_annual_gs_mmolH2Om2s
               grid_output$mean_gs_mmolH2Om2s[slot_i,slot_j,] = site_output$mean_gs_mmolH2Om2s
               grid_output$gs_mmolH2Om2s[n,,] = site_output$gs_mmolH2Om2s
           }
           if (any(check_list == "gb_mmolH2Om2s")) {
               # Canopy boundary layer conductance
               grid_output$mean_annual_gb_mmolH2Om2s[n,,] = site_output$mean_annual_gb_mmolH2Om2s
               grid_output$mean_gb_mmolH2Om2s[slot_i,slot_j,] = site_output$mean_gb_mmolH2Om2s
               grid_output$gb_mmolH2Om2s[n,,] = site_output$gb_mmolH2Om2s
           }
           if (any(check_list == "leaf_temperature_celcius") == TRUE) {
              # Canopy temperature 
              grid_output$mean_annual_leaf_temperature_celcius[n,,] = site_output$mean_annual_leaf_temperature_celcius
              grid_output$mean_leaf_temperature_celcius[slot_i,slot_j,] = site_output$mean_leaf_temperature_celcius 
              grid_output$leaf_temperature_celcius[n,,] = site_output$leaf_temperature_celcius
              # Soil temperature
              grid_output$mean_annual_soil_temperature_celcius[n,,] = site_output$mean_annual_soil_temperature_celcius
              grid_output$mean_soil_temperature_celcius[slot_i,slot_j,] = site_output$mean_soil_temperature_celcius
              grid_output$soil_temperature_celcius[n,,] = site_output$soil_temperature_celcius
           }

           # Parameter prior assimilated overlaps
           grid_output$priors_assim_data_overlap_fraction[,slot_i,slot_j] = site_output$priors_assim_data_overlap_fraction
           # Any time series assimilated data overlaps?
           if (any(check_list == "gpp_assim_data_overlap_fraction")) {
               grid_output$gpp_assim_data_overlap_fraction[slot_i,slot_j] = site_output$gpp_assim_data_overlap_fraction
           }
           if (any(check_list == "lai_assim_data_overlap_fraction")) {
               grid_output$lai_assim_data_overlap_fraction[slot_i,slot_j] = site_output$lai_assim_data_overlap_fraction
           }
           if (any(check_list == "nee_assim_data_overlap_fraction")) {
               grid_output$nee_assim_data_overlap_fraction[slot_i,slot_j] = site_output_all$nee_assim_data_overlap_fraction
           }
           if (any(check_list == "wood_assim_data_overlap_fraction")) {
               grid_output$wood_assim_data_overlap_fraction[slot_i,slot_j] = site_output$wood_assim_data_overlap_fraction
           }
           if (any(check_list == "soil_assim_data_overlap_fraction")) {
               grid_output$soil_assim_data_overlap_fraction[slot_i,slot_j] = site_output$soil_assim_data_overlap_fraction
           }
           if (any(check_list == "fapar_assim_data_overlap_fraction")) {
               grid_output$fapar_assim_data_overlap_fraction[slot_i,slot_j] = site_output$fapar_assim_data_overlap_fraction
           }             
           if (any(check_list == "et_assim_data_overlap_fraction")) {
               grid_output$et_assim_data_overlap_fraction[slot_i,slot_j] = site_output$et_assim_data_overlap_fraction
           }
           if (any(check_list == "nbe_assim_data_overlap_fraction")) {
               grid_output$nbe_assim_data_overlap_fraction[slot_i,slot_j] = site_output$nbe_assim_data_overlap_fraction
           }
           if (any(check_list == "fire_assim_data_overlap_fraction")) {
               grid_output$fire_assim_data_overlap_fraction[slot_i,slot_j] = site_output$fire_assim_data_overlap_fraction
           }

           # Store mean absolute parameter correlation information
           grid_output$absolute_mean_parameter_correlation[slot_i,slot_j] = site_output$absolute_mean_parameter_correlation
           # Parameter vs C-cycle flux correlation across ensemble member
           grid_output$lai_parameter_correlation[slot_i,slot_j,] = site_output$lai_parameter_correlation
           grid_output$nbp_parameter_correlation[slot_i,slot_j,] = site_output$nbp_parameter_correlation           
           grid_output$nee_parameter_correlation[slot_i,slot_j,] = site_output$nee_parameter_correlation
           grid_output$gpp_parameter_correlation[slot_i,slot_j,] = site_output$gpp_parameter_correlation
           grid_output$rauto_parameter_correlation[slot_i,slot_j,] = site_output$rauto_parameter_correlation
           grid_output$rhet_parameter_correlation[slot_i,slot_j,] = site_output$rhet_parameter_correlation
           grid_output$fire_parameter_correlation[slot_i,slot_j,] = site_output$fire_parameter_correlation
           # Correlations between LAI and key gross and net fluxes
           grid_output$lai_m2m2_to_GPP_gCm2day_correlation[slot_i,slot_j] = site_output$lai_m2m2_to_GPP_gCm2day_correlation
           grid_output$lai_m2m2_to_NEE_gCm2day_correlation[slot_i,slot_j] = site_output$lai_m2m2_to_NEE_gCm2day_correlation
           grid_output$lai_m2m2_to_NBP_gCm2day_correlation[slot_i,slot_j] = site_output$lai_m2m2_to_NBP_gCm2day_correlation           
           grid_output$lai_m2m2_to_Rauto_gCm2day_correlation[slot_i,slot_j] = site_output$lai_m2m2_to_Rauto_gCm2day_correlation
           grid_output$lai_m2m2_to_Rhet_gCm2day_correlation[slot_i,slot_j] = site_output$lai_m2m2_to_Rhet_gCm2day_correlation
           grid_output$lai_m2m2_to_som_gCm2_correlation[slot_i,slot_j] = site_output$lai_m2m2_to_som_gCm2_correlation
           grid_output$lai_m2m2_to_dCsom_gCm2_correlation[slot_i,slot_j] = site_output$lai_m2m2_to_dCsom_gCm2_correlation
           # Correlations between NBP and key gross and net fluxes
           grid_output$NBP_gCm2day_to_GPP_gCm2day_correlation[slot_i,slot_j] = site_output$NBP_gCm2day_to_GPP_gCm2day_correlation
           grid_output$NBP_gCm2day_to_NEE_gCm2day_correlation[slot_i,slot_j] = site_output$NBP_gCm2day_to_NEE_gCm2day_correlation
           grid_output$NBP_gCm2day_to_lai_m2m2_correlation[slot_i,slot_j] = site_output$NBP_gCm2day_to_lai_m2m2_correlation           
           grid_output$NBP_gCm2day_to_Rauto_gCm2day_correlation[slot_i,slot_j] = site_output$NBP_gCm2day_to_Rauto_gCm2day_correlation
           grid_output$NBP_gCm2day_to_Rhet_gCm2day_correlation[slot_i,slot_j] = site_output$NBP_gCm2day_to_Rhet_gCm2day_correlation
           grid_output$NBP_gCm2day_to_som_gCm2_correlation[slot_i,slot_j] = site_output$NBP_gCm2day_to_som_gCm2_correlation
           grid_output$NBP_gCm2day_to_dCsom_gCm2_correlation[slot_i,slot_j] = site_output$NBP_gCm2day_to_dCsom_gCm2_correlation  

           # If wood is estimated
           if (any(check_list == "wood_gCm2")) {           
               grid_output$lai_m2m2_to_wood_gCm2_correlation[slot_i,slot_j] = site_output$lai_m2m2_to_wood_gCm2_correlation
               grid_output$lai_m2m2_to_dCwood_gCm2_correlation[slot_i,slot_j] = site_output$lai_m2m2_to_dCwood_gCm2_correlation
               grid_output$NBP_gCm2day_to_wood_gCm2_correlation[slot_i,slot_j] = site_output$NBP_gCm2day_to_wood_gCm2_correlation
               grid_output$NBP_gCm2day_to_dCwood_gCm2_correlation[slot_i,slot_j] = site_output$NBP_gCm2day_to_dCwood_gCm2_correlation
           }   
           # If harvest is estimated
           if (any(check_list == "harvest_gCm2day")) {
               grid_output$lai_m2m2_to_harvest_gCm2day_correlation[slot_i,slot_j] = site_output$lai_m2m2_to_harvest_gCm2day_correlation
               grid_output$NBP_gCm2day_to_harvest_gCm2day_correlation[slot_i,slot_j] = site_output$NBP_gCm2day_to_harvest_gCm2day_correlation
           }           
           # If Mean transit time for wood correlation exists, ensure we store it for the gridded run too
           if (any(check_list == "CiCa_parameter_correlation")) {
               grid_output$CiCa_parameter_correlation[slot_i,slot_j,] = site_output$CiCa_parameter_correlation
           }
           # If Mean transit time for wood correlation exists, ensure we store it for the gridded run too
           if (any(check_list == "LWP_parameter_correlation")) {
               grid_output$LWP_parameter_correlation[slot_i,slot_j,] = site_output$LWP_parameter_correlation
           }           
           # If Mean transit time for wood correlation exists, ensure we store it for the gridded run too
           if (any(check_list == "MTT_wood_years_parameter_correlation")) {
               grid_output$MTT_wood_years_parameter_correlation[slot_i,slot_j,] = site_output$MTT_wood_years_parameter_correlation
           }
           # If Mean mean allocation to wood correlation exists, ensure we store it for the gridded run too
           if (any(check_list == "NPP_wood_gCm2day_parameter_correlation")) {
               grid_output$NPP_wood_gCm2day_parameter_correlation[slot_i,slot_j,] = site_output$NPP_wood_gCm2day_parameter_correlation
           }
           # If the correlation between wood MTT and wood allocation have been determined
           if (any(check_list == "MTT_wood_years_to_NPP_wood_gCm2day_correlation")) {
               grid_output$MTT_wood_years_to_NPP_wood_gCm2day_correlation[slot_i,slot_j] = site_output$MTT_wood_years_to_NPP_wood_gCm2day_correlation
               grid_output$MTT_wood_years_to_NPP_wood_fraction_correlation[slot_i,slot_j] = site_output$MTT_wood_years_to_NPP_wood_fraction_correlation
               grid_output$MTT_wood_years_to_MTT_som_years_correlation[slot_i,slot_j] = site_output$MTT_wood_years_to_MTT_som_years_correlation
               # ...against key other observables
               grid_output$NPP_wood_gCm2day_to_GPP_gCm2day_correlation[slot_i,slot_j] = site_output$NPP_wood_gCm2day_to_GPP_gCm2day_correlation
               grid_output$NPP_wood_gCm2day_to_NEE_gCm2day_correlation[slot_i,slot_j] = site_output$NPP_wood_gCm2day_to_NEE_gCm2day_correlation
               grid_output$NPP_wood_gCm2day_to_Rauto_gCm2day_correlation[slot_i,slot_j] = site_output$NPP_wood_gCm2day_to_Rauto_gCm2day_correlation
               grid_output$NPP_wood_gCm2day_to_Rhet_gCm2day_correlation[slot_i,slot_j] = site_output$NPP_wood_gCm2day_to_Rhet_gCm2day_correlation
               grid_output$NPP_wood_gCm2day_to_wood_gCm2_correlation[slot_i,slot_j] = site_output$NPP_wood_gCm2day_to_wood_gCm2_correlation
               grid_output$NPP_wood_gCm2day_to_som_gCm2_correlation[slot_i,slot_j] = site_output$NPP_wood_gCm2day_to_som_gCm2_correlation
               grid_output$NPP_wood_gCm2day_to_lai_m2m2_correlation[slot_i,slot_j] = site_output$NPP_wood_gCm2day_to_lai_m2m2_correlation            
               grid_output$NPP_wood_gCm2day_to_dCwood_gCm2_correlation[slot_i,slot_j] = site_output$NPP_wood_gCm2day_to_dCwood_gCm2_correlation  
               grid_output$NPP_wood_gCm2day_to_dCsom_gCm2_correlation[slot_i,slot_j] = site_output$NPP_wood_gCm2day_to_dCsom_gCm2_correlation
               grid_output$NPP_wood_fraction_to_GPP_gCm2day_correlation[slot_i,slot_j] = site_output$NPP_wood_fraction_to_GPP_gCm2day_correlation
               grid_output$NPP_wood_fraction_to_NEE_gCm2day_correlation[slot_i,slot_j] = site_output$NPP_wood_fraction_to_NEE_gCm2day_correlation
               grid_output$NPP_wood_fraction_to_Rauto_gCm2day_correlation[slot_i,slot_j] = site_output$NPP_wood_fraction_to_Rauto_gCm2day_correlation
               grid_output$NPP_wood_fraction_to_Rhet_gCm2day_correlation[slot_i,slot_j] = site_output$NPP_wood_fraction_to_Rhet_gCm2day_correlation
               grid_output$NPP_wood_fraction_to_wood_gCm2_correlation[slot_i,slot_j] = site_output$NPP_wood_fraction_to_wood_gCm2_correlation
               grid_output$NPP_wood_fraction_to_som_gCm2_correlation[slot_i,slot_j] = site_output$NPP_wood_fraction_to_som_gCm2_correlation
               grid_output$NPP_wood_fraction_to_lai_m2m2_correlation[slot_i,slot_j] = site_output$NPP_wood_fraction_to_lai_m2m2_correlation            
               grid_output$NPP_wood_fraction_to_dCwood_gCm2_correlation[slot_i,slot_j] = site_output$NPP_wood_fraction_to_dCwood_gCm2_correlation
               grid_output$NPP_wood_fraction_to_dCsom_gCm2_correlation[slot_i,slot_j] = site_output$NPP_wood_fraction_to_dCsom_gCm2_correlation
               grid_output$MTT_wood_years_to_GPP_gCm2day_correlation[slot_i,slot_j] = site_output$MTT_wood_years_to_GPP_gCm2day_correlation
               grid_output$MTT_wood_years_to_NEE_gCm2day_correlation[slot_i,slot_j] = site_output$MTT_wood_years_to_NEE_gCm2day_correlation
               grid_output$MTT_wood_years_to_Rauto_gCm2day_correlation[slot_i,slot_j] = site_output$MTT_wood_years_to_Rauto_gCm2day_correlation 
               grid_output$MTT_wood_years_to_Rhet_gCm2day_correlation[slot_i,slot_j] = site_output$MTT_wood_years_to_Rhet_gCm2day_correlation
               grid_output$MTT_wood_years_to_wood_gCm2_correlation[slot_i,slot_j] = site_output$MTT_wood_years_to_wood_gCm2_correlation
               grid_output$MTT_wood_years_to_som_gCm2_correlation[slot_i,slot_j] = site_output$MTT_wood_years_to_som_gCm2_correlation
               grid_output$MTT_wood_years_to_lai_m2m2_correlation[slot_i,slot_j] = site_output$MTT_wood_years_to_lai_m2m2_correlation     
               grid_output$MTT_wood_years_to_dCwood_gCm2_correlation[slot_i,slot_j] = site_output$MTT_wood_years_to_dCwood_gCm2_correlation
               grid_output$MTT_wood_years_to_dCsom_gCm2_correlation[slot_i,slot_j] = site_output$MTT_wood_years_to_dCsom_gCm2_correlation
               grid_output$dCwood_gCm2_to_gpp_gCm2day_correlation[slot_i,slot_j] = site_output$dCwood_gCm2_to_gpp_gCm2day_correlation 
               grid_output$dCwood_gCm2_to_rauto_gCm2day_correlation[slot_i,slot_j] = site_output$dCwood_gCm2_to_rauto_gCm2day_correlation 
               grid_output$dCwood_gCm2_to_nee_gCm2day_correlation[slot_i,slot_j]  = site_output$dCwood_gCm2_to_nee_gCm2day_correlation 
               grid_output$dCwood_gCm2_to_rhet_gCm2day_correlation[slot_i,slot_j] = site_output$dCwood_gCm2_to_rhet_gCm2day_correlation
               grid_output$dCwood_gCm2_to_wood_gCm2_correlation[slot_i,slot_j] = site_output$dCwood_gCm2_to_wood_gCm2_correlation 
               grid_output$dCwood_gCm2_to_som_gCm2_correlation[slot_i,slot_j] = site_output$dCwood_gCm2_to_som_gCm2_correlation
               grid_output$dCwood_gCm2_to_dCsom_gCm2_correlation[slot_i,slot_j] = site_output$dCwood_gCm2_to_dCsom_gCm2_correlation
           }
           if (any(check_list == "MTT_som_years_to_dCsom_gCm2_correlation") == TRUE) {
               grid_output$MTT_som_years_parameter_correlation[slot_i,slot_j,] = site_output$MTT_som_years_parameter_correlation
               # Assess within pixel correlations with soil turnover
               grid_output$MTT_som_years_to_GPP_gCm2day_correlation[slot_i,slot_j] = site_output$MTT_som_years_to_GPP_gCm2day_correlation
               grid_output$MTT_som_years_to_NEE_gCm2day_correlation[slot_i,slot_j] = site_output$MTT_som_years_to_NEE_gCm2day_correlation
               grid_output$MTT_som_years_to_Rauto_gCm2day_correlation[slot_i,slot_j] = site_output$MTT_som_years_to_Rauto_gCm2day_correlation 
               grid_output$MTT_som_years_to_Rhet_gCm2day_correlation[slot_i,slot_j] = site_output$MTT_som_years_to_Rhet_gCm2day_correlation
               grid_output$MTT_som_years_to_som_gCm2_correlation[slot_i,slot_j] = site_output$MTT_som_years_to_som_gCm2_correlation
               grid_output$MTT_som_years_to_lai_m2m2_correlation[slot_i,slot_j] = site_output$MTT_som_years_to_lai_m2m2_correlation
               grid_output$MTT_som_years_to_dCsom_gCm2_correlation[slot_i,slot_j] = site_output$MTT_som_years_to_dCsom_gCm2_correlation
               grid_output$dCsom_gCm2_to_gpp_gCm2day_correlation[slot_i,slot_j] = site_output$dCsom_gCm2_to_gpp_gCm2day_correlation 
               grid_output$dCsom_gCm2_to_rauto_gCm2day_correlation[slot_i,slot_j] = site_output$dCsom_gCm2_to_rauto_gCm2day_correlation 
               grid_output$dCsom_gCm2_to_nee_gCm2day_correlation[slot_i,slot_j]  = site_output$dCsom_gCm2_to_nee_gCm2day_correlation 
               grid_output$dCsom_gCm2_to_rhet_gCm2day_correlation[slot_i,slot_j] = site_output$dCsom_gCm2_to_rhet_gCm2day_correlation
               grid_output$dCsom_gCm2_to_som_gCm2_correlation[slot_i,slot_j] = site_output$dCsom_gCm2_to_som_gCm2_correlation   
               grid_output$dCsom_gCm2_to_som_input_gCm2_correlation[slot_i,slot_j] = site_output$dCsom_gCm2_to_som_input_gCm2_correlation   

               if (any(check_list == "MTT_wood_years_to_NPP_wood_gCm2day_correlation")) {
                   grid_output$MTT_som_years_to_wood_gCm2_correlation[slot_i,slot_j] = site_output$MTT_som_years_to_wood_gCm2_correlation
                   grid_output$MTT_som_years_to_dCwood_gCm2_correlation[slot_i,slot_j] = site_output$MTT_som_years_to_dCwood_gCm2_correlation                    
                   grid_output$dCsom_gCm2_to_wood_gCm2_correlation[slot_i,slot_j]  = site_output$dCsom_gCm2_to_wood_gCm2_correlation                
               }

           }            

           # Tidy up
           rm(site_output) ; file.remove(site_output_all[[n]])

       } # Does the current site have processed output?

  } # loop sites to load into grid_output

  # Return back to user
  return(grid_output)
  
} # end function post_process_into_grid

## Use byte compile
post_process_into_grid<-cmpfun(post_process_into_grid)
