
###
## Function to load potential biomass maps to be used as attractor
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_potential_biomass_maps_for_extraction<-function(latlon_in,Cwood_potential_source,start,finish,timestep_days) {

   if (Cwood_potential_source == "UoE_potAGB") {

       # Read in desired potential AGB map
       infile = nc_open(paste(path_to_Cwood_potential,"BRA_020_AGB_potential_RFR_avitabile_worldclim_soilgrids_dist2edges_final.nc", sep=""))
       lat = ncvar_get(infile, "lat") ; long = ncvar_get(infile, "lon")
       # Estimates are in Mg/ha, convert to gC/m2 later
       biomass_gCm2 = ncvar_get(infile, "AGBpot")
       biomass_uncertainty_gCm2 = ncvar_get(infile, "AGBpot_max")
       tmp = ncvar_get(infile, "AGBpot_min")
       biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 - tmp ; rm(tmp)

       # Create lat / long grid from vectors
       idim = length(long) ; jdim = length(lat)
       lat = array(lat, dim=c(jdim,idim)) ; lat = t(lat)
       long = array(long, dim=c(idim,jdim))

       # Use allometry to estimate below ground biomass stock and
       # combined with the above ground (Mg/ha) to give a total woody biomass estimate
       # Saatchi et al., (2011), PNAS, 108, 9899-9904, https://www.pnas.org/content/108/24/9899
       biomass_gCm2 = biomass_gCm2 + (0.489 * biomass_gCm2 ** 0.89)
       biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 + (0.489 * biomass_uncertainty_gCm2 ** 0.89)
       # Convert units of biomass and its uncertainty from MgCha -> gC/m2
       biomass_gCm2 = biomass_gCm2 * 1e2 * 0.48
       biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 1e2 * 0.48

       # Output variables
       return(list(lat = lat, long = long,
                   biomass_gCm2 = biomass_gCm2, biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2))

   } else {

       # Output variables
       return(list(lat = -9999, long = -9999,
                   biomass_gCm2 = -9999, biomass_uncertainty_gCm2 = -9999))

   } # Cwood_potential_source


} # function end
