
###
## Function to post-process CARDAMOM output for a gridded analysis
###

# This function was created by T. L Smallman (t.l.smallman@ed.ac.uk, UoE)

post_process_for_grid<-function(outfile_stock_fluxes,PROJECT,drivers,parameters,num_quantiles,na_flag,converged,states_all) {

  # Determine some useful information for the analysis below
  nos_years = PROJECT$nos_years
  steps_per_year = floor(dim(drivers$met)[1] / nos_years)

  # Declare the site level output list object
  site_output = list(num_quantiles = num_quantiles, steps_per_year = steps_per_year, nos_years = nos_years)

  # Check the list variables in states_all which we will be searching
  check_list = names(states_all)

  ###
  # Derive stocks and fluxes used in the calculation of gridded aggregates
  # These are variables which for a site analysis would be easy to calcuate
  # from the ensembles but difficult if not determined here and now before aggregation
  ###

  # The total allocation of C to the foliage pool can, depending on model,
  # be the combined total of direct allocation and that via a labile pool.
  # For many comparison we will need their combined total.
  if (any(check_list == "alloc_foliage_gCm2day") &
      any(check_list == "labile_to_foliage_gCm2day")) {
      states_all$combined_alloc_foliage_gCm2day = states_all$alloc_foliage_gCm2day + states_all$labile_to_foliage_gCm2day
  } else {
      if (any(check_list == "labile_to_foliage_gCm2day")) {
          states_all$combined_alloc_foliage_gCm2day = states_all$labile_to_foliage_gCm2day
      } else if (any(check_list == "alloc_foliage_gCm2day")) {
          states_all$combined_alloc_foliage_gCm2day = states_all$alloc_foliage_gCm2day
      } else {
          stop("Error, CARDAMOM cannnot determine where C allocation foliage has come from")
      }
  }
  # Update the list variables in states_all which we will be searching
  check_list = names(states_all)


  # Determine pool specific

  # Determine the total biomass within the system
  # Concurrently determine the combined natural, fire and harvest litter fluxes.
  # All models have a foliage pool
  states_all$biomass_gCm2 = states_all$foliage_gCm2
  states_all$biomass_to_litter_gCm2day = states_all$foliage_to_litter_gCm2day
  # Now look for accumulating optional disturbance drivers
  if (any(check_list == "FIREemiss_foliage_gCm2day")) {
      states_all$FIREemiss_biomass_gCm2day = states_all$FIREemiss_foliage_gCm2day
  }
  if (any(check_list == "FIRElitter_foliage_gCm2day")) {
      states_all$FIRElitter_biomass_gCm2day = states_all$FIRElitter_foliage_gCm2day
  }
  if (any(check_list == "HARVESTextracted_foliage_gCm2day")) {
      states_all$HARVESTextracted_biomass_gCm2day = states_all$HARVESTextracted_foliage_gCm2day
      states_all$HARVESTlitter_biomass_gCm2day = states_all$HARVESTlitter_foliage_gCm2day
  }
  if (any(check_list == "labile_gCm2")) {
      states_all$biomass_gCm2 = states_all$biomass_gCm2 + states_all$labile_gCm2
      # Now look for accumulating optional disturbance drivers
      if (any(check_list == "FIREemiss_labile_gCm2day")) {
          states_all$FIREemiss_biomass_gCm2day = states_all$FIREemiss_biomass_gCm2day + states_all$FIREemiss_labile_gCm2day
      }
      if (any(check_list == "FIRElitter_labile_gCm2day")) {
          states_all$FIRElitter_biomass_gCm2day = states_all$FIRElitter_biomass_gCm2day + states_all$FIRElitter_labile_gCm2day
      }
      if (any(check_list == "HARVESTextracted_labile_gCm2day")) {
          states_all$HARVESTextracted_biomass_gCm2day = states_all$HARVESTextracted_biomass_gCm2day + states_all$HARVESTextracted_labile_gCm2day
          states_all$HARVESTlitter_biomass_gCm2day = states_all$HARVESTlitter_biomass_gCm2day + states_all$HARVESTlitter_labile_gCm2day
      }
  }
  if (any(check_list == "roots_gCm2")) {
      states_all$biomass_gCm2 = states_all$biomass_gCm2 + states_all$roots_gCm2
      states_all$biomass_to_litter_gCm2day = states_all$biomass_to_litter_gCm2day + states_all$roots_to_litter_gCm2day
      # Now look for accumulating optional disturbance drivers
      if (any(check_list == "FIREemiss_roots_gCm2day")) {
          states_all$FIREemiss_biomass_gCm2day = states_all$FIREemiss_biomass_gCm2day + states_all$FIREemiss_roots_gCm2day
      }
      if (any(check_list == "FIRElitter_roots_gCm2day")) {
          states_all$FIRElitter_biomass_gCm2day = states_all$FIRElitter_biomass_gCm2day + states_all$FIRElitter_roots_gCm2day
      }
      if (any(check_list == "HARVESTextracted_roots_gCm2day")) {
          states_all$HARVESTextracted_biomass_gCm2day = states_all$HARVESTextracted_biomass_gCm2day + states_all$HARVESTextracted_roots_gCm2day
          states_all$HARVESTlitter_biomass_gCm2day = states_all$HARVESTlitter_biomass_gCm2day + states_all$HARVESTlitter_roots_gCm2day
      }
  }
  if (any(check_list == "wood_gCm2")) {
      states_all$biomass_gCm2 = states_all$biomass_gCm2 + states_all$wood_gCm2
      states_all$biomass_to_litter_gCm2day = states_all$biomass_to_litter_gCm2day + states_all$wood_to_litter_gCm2day
      # Account for wood lost as emission due to fire
      if (any(check_list == "FIREemiss_wood_gCm2day")) {
          states_all$FIREemiss_biomass_gCm2day = states_all$FIREemiss_biomass_gCm2day + states_all$FIREemiss_wood_gCm2day
      }
      # Account for wood lost as mortality induced litter / non-combusted residues due to fire
      if (any(check_list == "FIRElitter_wood_gCm2day")) {
          states_all$FIRElitter_biomass_gCm2day = states_all$FIRElitter_biomass_gCm2day + states_all$FIRElitter_wood_gCm2day
      }
      # Account for wood extracted by harvest activities
      if (any(check_list == "HARVESTextracted_wood_gCm2day")) {
          states_all$HARVESTextracted_biomass_gCm2day = states_all$HARVESTextracted_biomass_gCm2day + states_all$HARVESTextracted_wood_gCm2day
          states_all$HARVESTlitter_biomass_gCm2day = states_all$HARVESTlitter_biomass_gCm2day + states_all$HARVESTlitter_wood_gCm2day
      }
  } # wood
  # Update the list variables in states_all which we will be searching
  check_list = names(states_all)

  # Determine the total dead organic matter within the system
  # plus their associated loss terms.
  # All models have a som pool so start there
  if (any(check_list == "dom_gCm2") == FALSE) {
      states_all$dom_gCm2 = states_all$som_gCm2
      states_all$rhet_dom_gCm2day = states_all$rhet_som_gCm2day
      # Create corresponding variables for the disturbance drivers.
      # NOTE: that we do not need to account for disturbance induced
      # litter fluxes as these always remain within the DOM pool
      if (any(check_list == "FIREemiss_som_gCm2day")) {
          states_all$FIREemiss_dom_gCm2day = states_all$FIREemiss_som_gCm2day
      }
      if (any(check_list == "HARVESTextracted_som_gCm2day")) {
          states_all$HARVESTextracted_dom_gCm2day = states_all$HARVESTextracted_som_gCm2day
      }
      # Check for presence of litter or wood litter pools
      if (any(check_list == "litter_gCm2")) {
          states_all$dom_gCm2 = states_all$dom_gCm2 + states_all$litter_gCm2
          states_all$rhet_dom_gCm2day = states_all$rhet_dom_gCm2day + states_all$rhet_litter_gCm2day
          if (any(check_list == "FIREemiss_litter_gCm2day")) {
              states_all$FIREemiss_dom_gCm2day = states_all$FIREemiss_dom_gCm2day + states_all$FIREemiss_litter_gCm2day
          }
          if (any(check_list == "HARVESTextracted_litter_gCm2day")) {
              states_all$HARVESTextracted_dom_gCm2day = states_all$HARVESTextracted_dom_gCm2day + states_all$HARVESTextracted_litter_gCm2day
          }
      } # litter
      if (any(check_list == "woodlitter_gCm2")) {
          states_all$dom_gCm2 = states_all$dom_gCm2 + states_all$woodlitter_gCm2
          states_all$rhet_dom_gCm2day = states_all$rhet_dom_gCm2day + states_all$rhet_woodlitter_gCm2day
          if (any(check_list == "FIREemiss_woodlitter_gCm2day")) {
              states_all$FIREemiss_dom_gCm2day = states_all$FIREemiss_dom_gCm2day + states_all$FIREemiss_woodlitter_gCm2day
          }
          if (any(check_list == "HARVESTextracted_woodlitter_gCm2day")) {
              states_all$HARVESTextracted_dom_gCm2day = states_all$HARVESTextracted_dom_gCm2day + states_all$HARVESTextracted_woodlitter_gCm2day
          }
      } # wood litter
  } # DOM pool exists

  # Determine the total C within the system
  # Biomass and dead organic matter have already been determined,
  # so combine.
  states_all$Ctotal_gCm2 = states_all$biomass_gCm2 + states_all$dom_gCm2
  # Update the list variables in states_all which we will be searching
  check_list = names(states_all)

  ###
  # Extract likelihood / parameter / driver information
  ###

  # loop through parameters + likelihood
  site_output$parameters = array(NA, dim=c(dim(parameters)[1],length(num_quantiles)))
  for (p in seq(1, dim(parameters)[1])) {
       site_output$parameters[p,] = quantile(as.vector(parameters[p,,]), prob=num_quantiles)
  }
  # track which parameters have converged + likelihood
  site_output$parameters_converged = rep(0, dim(parameters)[1])
  site_output$parameters_converged[which(converged == "PASS")] = 1
  # Mean meteorological conditions
  site_output$mean_temperature_C = mean((drivers$met[,3]+drivers$met[,2])*0.5)
  site_output$mean_radiation_MJm2day = mean(drivers$met[,4])
  site_output$mean_vpd_Pa = mean(drivers$met[,16])
  site_output$mean_precipitation_kgH2Om2yr = mean(drivers$met[,7])*86400*365.25
  # Assimilated LAI information
  if (max(drivers$obs[,3]) > 0) {
      filter = which(drivers$obs[,3] != -9999)
      site_output$assimilated_lai_max_m2m2 = max(drivers$obs[filter,3])
      site_output$assimilated_lai_mean_m2m2 = mean(drivers$obs[filter,3])
      site_output$assimilated_lai_sd_m2m2 = sd(drivers$obs[filter,3])
      site_output$assimilated_lai_unc_m2m2 = mean(drivers$obs[filter,4])
  } else {
      site_output$assimilated_lai_max_m2m2 = NA
      site_output$assimilated_lai_mean_m2m2 = NA
      site_output$assimilated_lai_sd_m2m2 = NA
      site_output$assimilated_lai_unc_m2m2 = NA
  }
  # Assimilated wood stock / prior information
  # Do we have one or both wood stock prior and time series inforamtion
  if (drivers$parpriors[21] > 0 |
      length(which(drivers$obs[,13] > -9999)) > 0) {
      # Initialise the variable with zero to allow averaging
      site_output$assimilated_wood_mean_gCm2 = 0
      site_output$assimilated_wood_mean_unc_gCm2 = 0
      # Do we have time series information?
      if (length(which(drivers$obs[,13] > -9999)) > 0) { # Does time series info exist
          site_output$assimilated_wood_mean_gCm2 = append(site_output$assimilated_wood_mean_gCm2,
                                                          drivers$obs[which(drivers$obs[,13] != -9999),13])
          site_output$assimilated_wood_mean_unc_gCm2 = append(site_output$assimilated_wood_mean_unc_gCm2,
                                                              drivers$obs[which(drivers$obs[,14] != -9999),14])
      }
      # Do we have a prior estimate on initial wood stocks?
      if (drivers$parpriors[21] > 0) { # Does a prior exist
          site_output$assimilated_wood_mean_gCm2 = append(site_output$assimilated_wood_mean_gCm2,
                                                          drivers$parpriors[21])
          site_output$assimilated_wood_mean_unc_gCm2 = append(site_output$assimilated_wood_mean_unc_gCm2,
                                                              drivers$parpriorunc[21])
      }
      # Average available observation and uncertainty information
      # remembering to remove the initial value
      site_output$assimilated_wood_mean_gCm2 = mean(site_output$assimilated_wood_mean_gCm2[-1], na.rm = na_flag)
      site_output$assimilated_wood_mean_unc_gCm2 = mean(site_output$assimilated_wood_mean_unc_gCm2[-1], na.rm = na_flag)
  } else {
      # assign missing value flag for consistency
      site_output$assimilated_wood_mean_gCm2 = NA
      site_output$assimilated_wood_mean_unc_gCm2 = NA
  } # wood stock or prior information was assimilated

  # Assimilated som stock / prior information
  # Do we have one or both som stock prior and time series inforamtion
  if (drivers$parpriors[23] > 0 |
      length(which(drivers$obs[,19] > -9999)) > 0) {
      # Initialise the variable with zero to allow averaging
      site_output$assimilated_som_mean_gCm2 = 0
      site_output$assimilated_som_mean_unc_gCm2 = 0
      # Do we have time series information?
      if (length(which(drivers$obs[,19] > -9999)) > 0) { # Does time series info exist
          site_output$assimilated_som_mean_gCm2 = append(site_output$assimilated_som_mean_gCm2,
                                                         drivers$obs[which(drivers$obs[,19] != -9999),19])
          site_output$assimilated_som_mean_unc_gCm2 = append(site_output$assimilated_som_mean_unc_gCm2,
                                                             drivers$obs[which(drivers$obs[,20] != -9999),20])
      }
      # Do we have a prior estimate on initial wood stocks?
      if (drivers$parpriors[23] > 0) { # Does a prior exist
          site_output$assimilated_som_mean_gCm2 = append(site_output$assimilated_som_mean_gCm2,
                                                         drivers$parpriors[23])
          site_output$assimilated_som_mean_unc_gCm2 = append(site_output$assimilated_som_mean_unc_gCm2,
                                                             drivers$parpriorunc[23])
      }
      # Average available observation and uncertainty information
      # remembering to remove the initial value
      site_output$assimilated_som_mean_gCm2 = mean(site_output$assimilated_som_mean_gCm2[-1], na.rm = na_flag)
      site_output$assimilated_som_mean_unc_gCm2 = mean(site_output$assimilated_som_mean_unc_gCm2[-1], na.rm = na_flag)
  } else {
      # assign missing value flag for consistency
      site_output$assimilated_som_mean_gCm2 = NA
      site_output$assimilated_som_mean_unc_gCm2 = NA
  } # som stock or prior information was assimilated

  ###
  # Begin aggregating the states_all information into
  # quantile based site_output
  ###

  # Net primary production allocation fractions
  if (any(check_list == "NPP_foliage_fraction")) {site_output$NPP_foliage_fraction = quantile(states_all$NPP_foliage_fraction, prob=num_quantiles, na.rm = na_flag)}
  if (any(check_list == "NPP_roots_fraction")) {site_output$NPP_roots_fraction = quantile(states_all$NPP_roots_fraction, prob=num_quantiles, na.rm = na_flag)}
  if (any(check_list == "NPP_wood_fraction")) {site_output$NPP_wood_fraction = quantile(states_all$NPP_wood_fraction, prob=num_quantiles, na.rm = na_flag)}
  # Analysis mean transit (residence) times (years)
  if (any(check_list == "MTT_labile_years")) {site_output$MTT_labile_years = quantile(states_all$MTT_labile_years, prob=num_quantiles, na.rm = na_flag)}
  if (any(check_list == "MTT_foliage_years")) {site_output$MTT_foliage_years = quantile(states_all$MTT_foliage_years, prob=num_quantiles, na.rm = na_flag)}
  if (any(check_list == "MTT_roots_years")) {site_output$MTT_roots_years = quantile(states_all$MTT_roots_years, prob=num_quantiles, na.rm = na_flag)}
  if (any(check_list == "MTT_wood_years")) {site_output$MTT_wood_years = quantile(states_all$MTT_wood_years, prob=num_quantiles, na.rm = na_flag)}
  if (any(check_list == "MTT_litter_years")) {site_output$MTT_litter_years = quantile(states_all$MTT_litter_years, prob=num_quantiles, na.rm = na_flag)}
  if (any(check_list == "MTT_woodlitter_years")) {site_output$MTT_woodlitter_years = quantile(states_all$MTT_woodlitter_years, prob=num_quantiles, na.rm = na_flag)}
  if (any(check_list == "MTT_som_years")) {site_output$MTT_som_years = quantile(states_all$MTT_som_years, prob=num_quantiles, na.rm = na_flag)}
  # Steady state C stock estimates (gC/m2)
  if (any(check_list == "SS_labile_gCm2")) {site_output$SS_labile_gCm2 = quantile(states_all$SS_labile_gCm2, prob=num_quantiles, na.rm = na_flag)}
  if (any(check_list == "SS_foliage_gCm2")) {site_output$SS_foliage_gCm2 = quantile(states_all$SS_foliage_gCm2, prob=num_quantiles, na.rm = na_flag)}
  if (any(check_list == "SS_roots_gCm2")) {site_output$SS_roots_gCm2 = quantile(states_all$SS_roots_gCm2, prob=num_quantiles, na.rm = na_flag)}
  if (any(check_list == "SS_wood_gCm2")) {site_output$SS_wood_gCm2 = quantile(states_all$SS_wood_gCm2, prob=num_quantiles, na.rm = na_flag)}
  if (any(check_list == "SS_litter_gCm2")) {site_output$SS_litter_gCm2 = quantile(states_all$SS_litter_gCm2, prob=num_quantiles, na.rm = na_flag)}
  if (any(check_list == "SS_woodlitter_gCm2")) {site_output$SS_woodlitter_gCm2 = quantile(states_all$SS_woodlitter_gCm2, prob=num_quantiles, na.rm = na_flag)}
  if (any(check_list == "SS_som_gCm2")) {site_output$SS_som_gCm2 = quantile(states_all$SS_som_gCm2, prob=num_quantiles, na.rm = na_flag)}

  # State variables - NOTE: extraction of pixel specific means done here to account for different ensemble trajectories, i.e. correlation in time.
  site_output$lai_m2m2                = apply(states_all$lai_m2m2,2,quantile,prob=num_quantiles,na.rm = na_flag)
  site_output$mean_lai_m2m2           = quantile(rowMeans(states_all$lai_m2m2,na.rm = na_flag) ,prob=num_quantiles, na.rm = TRUE)
  site_output$mean_annual_lai_m2m2    = apply(t(apply(states_all$lai_m2m2,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
  site_output$Ctotal_gCm2             = apply(states_all$Ctotal_gCm2,2,quantile,prob=num_quantiles,na.rm = na_flag)
  site_output$mean_Ctotal_gCm2        = quantile(rowMeans(states_all$Ctotal_gCm2,na.rm = na_flag) ,prob=num_quantiles, na.rm = TRUE)
  site_output$mean_annual_Ctotal_gCm2 = apply(t(apply(states_all$Ctotal_gCm2,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
  # Fluxes second
  site_output$gpp_gCm2day                 = apply(states_all$gpp_gCm2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
  site_output$mean_gpp_gCm2day            = quantile(rowMeans(states_all$gpp_gCm2day,na.rm = na_flag) ,prob=num_quantiles, na.rm = TRUE)
  site_output$mean_annual_gpp_gCm2day     = apply(t(apply(states_all$gpp_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
  site_output$rauto_gCm2day               = apply(states_all$rauto_gCm2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
  site_output$mean_rauto_gCm2day          = quantile(rowMeans(states_all$rauto_gCm2day,na.rm = na_flag) ,prob=num_quantiles, na.rm = TRUE)
  site_output$mean_annual_rauto_gCm2day   = apply(t(apply(states_all$rauto_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
  site_output$rhet_gCm2day                = apply(states_all$rhet_gCm2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
  site_output$mean_rhet_gCm2day           = quantile(rowMeans(states_all$rhet_gCm2day,na.rm = na_flag) ,prob=num_quantiles, na.rm = TRUE)
  site_output$mean_annual_rhet_gCm2day    = apply(t(apply(states_all$rhet_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
  site_output$nee_gCm2day                 = apply(states_all$nee_gCm2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
  site_output$mean_nee_gCm2day            = quantile(rowMeans(states_all$nee_gCm2day,na.rm = na_flag) ,prob=num_quantiles, na.rm = TRUE)
  site_output$mean_annual_nee_gCm2day     = apply(t(apply(states_all$nee_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
  site_output$npp_gCm2day                 = apply(states_all$npp_gCm2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
  site_output$mean_npp_gCm2day            = quantile(rowMeans(states_all$npp_gCm2day,na.rm = na_flag) ,prob=num_quantiles, na.rm = TRUE)
  site_output$mean_annual_npp_gCm2day     = apply(t(apply(states_all$npp_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
  site_output$nbe_gCm2day                 = apply(states_all$nbe_gCm2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
  site_output$mean_nbe_gCm2day            = quantile(rowMeans(states_all$nbe_gCm2day,na.rm = na_flag) ,prob=num_quantiles, na.rm = TRUE)
  site_output$mean_annual_nbe_gCm2day     = apply(t(apply(states_all$nbe_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
  site_output$nbp_gCm2day                 = apply(states_all$nbp_gCm2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
  site_output$mean_nbp_gCm2day            = quantile(rowMeans(states_all$nbp_gCm2day,na.rm = na_flag) ,prob=num_quantiles, na.rm = TRUE)
  site_output$mean_annual_nbp_gCm2day     = apply(t(apply(states_all$nbp_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
  site_output$reco_gCm2day                = apply(states_all$reco_gCm2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
  site_output$mean_reco_gCm2day           = quantile(rowMeans(states_all$reco_gCm2day,na.rm = na_flag) ,prob=num_quantiles, na.rm = TRUE)
  site_output$mean_annual_reco_gCm2day    = apply(t(apply(states_all$reco_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
  site_output$harvest_gCm2day             = apply(states_all$harvest_gCm2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
  site_output$mean_harvest_gCm2day        = quantile(rowMeans(states_all$harvest_gCm2day,na.rm = na_flag) ,prob=num_quantiles, na.rm = TRUE)
  site_output$mean_annual_harvest_gCm2day = apply(t(apply(states_all$harvest_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
  site_output$fire_gCm2day                = apply(states_all$fire_gCm2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
  site_output$mean_fire_gCm2day           = quantile(rowMeans(states_all$fire_gCm2day,na.rm = na_flag) ,prob=num_quantiles, na.rm = TRUE)
  site_output$mean_annual_fire_gCm2day    = apply(t(apply(states_all$fire_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
  # C-cycle diagnostics which are mean annuals apply vs rowMeans
  site_output$mean_annual_cue      = apply(states_all$mean_annual_cue, 2, quantile, prob = num_quantiles, na.rm=TRUE)
  site_output$mean_cue             = quantile(rowMeans(states_all$mean_annual_cue, na.rm=na_flag), prob = num_quantiles, na.rm=TRUE)

  ###
  # Track net pool change over time
  ###

  # Start with the common variables
  dCbio = states_all$Ctotal_gCm2 - states_all$Ctotal_gCm2[,1] # difference in total C from initial
  site_output$dCtotal_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)

  dCbio = states_all$lai_m2m2 - states_all$lai_m2m2[,1] #  difference in lai from initial
  site_output$dlai_m2m2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)

  # Then do model specific pool combinations
  # Include the pool, its net change, the allocation of C from GPP / NPP direct
  # and its output flow(s).

  # Labile related pool, change, input and output variables
  if (any(check_list == "labile_gCm2")) {
      # Assign pool to site_output
      site_output$labile_gCm2 = apply(states_all$labile_gCm2,2,quantile,prob=num_quantiles, na.rm = na_flag)
      site_output$mean_labile_gCm2 = quantile(rowMeans(states_all$labile_gCm2, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_labile_gCm2 = apply(t(apply(states_all$labile_gCm2,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Determine net pool change over time
      dCbio = states_all$labile_gCm2 - states_all$labile_gCm2[,1] # difference in labile from initial
      site_output$dClabile_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
      # Determine the allocation to labile - in all cases this must be a direct variable
      site_output$alloc_labile_gCm2day = apply(states_all$alloc_labile_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
      site_output$mean_alloc_labile_gCm2day = quantile(rowMeans(states_all$alloc_labile_gCm2day, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_alloc_labile_gCm2day = apply(t(apply(states_all$alloc_labile_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Declare combined natural, fire and harvest driven creation of litter
      site_output$combined_labile_to_litter_gCm2day = array(0, dim=dim(states_all$labile_gCm2))
      # Check for the possible loss pathways
      if (any(check_list == "labile_to_foliage_gCm2day")) {
          site_output$labile_to_foliage_gCm2day = apply(states_all$labile_to_foliage_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_labile_to_foliage_gCm2day = quantile(rowMeans(states_all$labile_to_foliage_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_labile_to_foliage_gCm2day = apply(t(apply(states_all$labile_to_foliage_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$NaturalFractionOfTurnover_labile = states_all$labile_to_foliage_gCm2day
          # Begin accumulating the total output
          site_output$outflux_labile_gCm2day = states_all$labile_to_foliage_gCm2day
      }
      # Other natural flux pathways should really go here before disturbance related
      if (any(check_list == "FIREemiss_labile_gCm2day")) {
          site_output$FIREemiss_labile_gCm2day = apply(states_all$FIREemiss_labile_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_FIREemiss_labile_gCm2day = quantile(rowMeans(states_all$FIREemiss_labile_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_FIREemiss_labile_gCm2day = apply(t(apply(states_all$FIREemiss_labile_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$FireFractionOfTurnover_labile = states_all$FIREemiss_labile_gCm2day
          site_output$outflux_labile_gCm2day = site_output$outflux_labile_gCm2day + states_all$FIREemiss_labile_gCm2day
      }
      if (any(check_list == "FIRElitter_labile_gCm2day")) {
          site_output$FIRElitter_labile_gCm2day = apply(states_all$FIRElitter_labile_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_FIRElitter_labile_gCm2day = quantile(rowMeans(states_all$FIRElitter_labile_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_FIRElitter_labile_gCm2day = apply(t(apply(states_all$FIRElitter_labile_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$outflux_labile_gCm2day = site_output$outflux_labile_gCm2day + states_all$FIRElitter_labile_gCm2day
          site_output$FireFractionOfTurnover_labile = site_output$FireFractionOfTurnover_labile + states_all$FIRElitter_labile_gCm2day
          # Accumulate the combined total of pool lost as litter,
          # NOTE depending on DALEC version may go to either a litter pool or som
          site_output$combined_labile_to_litter_gCm2day = site_output$combined_labile_to_litter_gCm2day + states_all$FIRElitter_labile_gCm2day
      }
      if (any(check_list == "HARVESTextracted_labile_gCm2day")) {
          site_output$HARVESTextracted_labile_gCm2day = apply(states_all$HARVESTextracted_labile_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_HARVESTextracted_labile_gCm2day = quantile(rowMeans(states_all$HARVESTextracted_labile_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_HARVESTextracted_labile_gCm2day = apply(t(apply(states_all$HARVESTextracted_labile_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$HARVESTlitter_labile_gCm2day = apply(states_all$HARVESTlitter_labile_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_HARVESTlitter_labile_gCm2day = quantile(rowMeans(states_all$HARVESTlitter_labile_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_HARVESTlitter_labile_gCm2day = apply(t(apply(states_all$HARVESTlitter_labile_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$outflux_labile_gCm2day = site_output$outflux_labile_gCm2day + states_all$HARVESTextracted_labile_gCm2day + states_all$HARVESTlitter_labile_gCm2day
          site_output$HarvestFractionOfTurnover_labile = states_all$HARVESTextracted_labile_gCm2day + states_all$HARVESTlitter_labile_gCm2day
          # Accumulate the combined total of pool lost as litter,
          # NOTE depending on DALEC version may go to either a litter pool or som
          site_output$combined_labile_to_litter_gCm2day = site_output$combined_labile_to_litter_gCm2day + states_all$HARVESTlitter_labile_gCm2day
      }
      # Use this information to determine the mean residence times as it evolves over time
      # NOTE: rollapply inverts the dimensions from that wanted
      site_output$MTT_annual_labile_years = t( apply(states_all$labile_gCm2,1, rollapply_mean_annual, step = steps_per_year)
                                             / (apply(site_output$outflux_labile_gCm2day,1, rollapply_mean_annual, step = steps_per_year) * 365.25) )
      site_output$MTT_annual_labile_years = apply(site_output$MTT_annual_labile_years,2,quantile,prob=num_quantiles, na.rm = na_flag)
      # Convert to fractions
      tmp = rowMeans(site_output$outflux_labile_gCm2day, na.rm = na_flag)
      site_output$NaturalFractionOfTurnover_labile = quantile(rowMeans(site_output$NaturalFractionOfTurnover_labile, na.rm = na_flag) / tmp,
                                                              prob = num_quantiles, na.rm = na_flag)
      site_output$FireFractionOfTurnover_labile = quantile(rowMeans(site_output$FireFractionOfTurnover_labile, na.rm = na_flag) / tmp,
                                                           prob = num_quantiles, na.rm = na_flag)
      site_output$HarvestFractionOfTurnover_labile = quantile(rowMeans(site_output$HarvestFractionOfTurnover_labile, na.rm = na_flag) / tmp,
                                                              prob = num_quantiles, na.rm = na_flag)
      # Mean outflux from labile
      site_output$mean_outflux_labile_gCm2day = quantile(rowMeans(site_output$outflux_labile_gCm2day, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_combined_labile_to_litter_gCm2day = quantile(rowMeans(site_output$combined_labile_to_litter_gCm2day, na.rm = na_flag), prob=num_quantiles)
      # Mean annual outflux from labile
      site_output$mean_annual_outflux_labile_gCm2day = apply(t(apply(site_output$outflux_labile_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      site_output$mean_annual_combined_labile_to_litter_gCm2day = apply(t(apply(site_output$combined_labile_to_litter_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Now aggregate across quantiles
      site_output$outflux_labile_gCm2day = apply(site_output$outflux_labile_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
      site_output$combined_labile_to_litter_gCm2day = apply(site_output$combined_labile_to_litter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
  } # exists("labile")

  # Foliage related pool, change, input and output variables
  if (any(check_list == "foliage_gCm2")) {
      # A combined total of C to foliage must always exist
      site_output$combined_alloc_foliage_gCm2day = apply(states_all$combined_alloc_foliage_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
      site_output$mean_combined_alloc_foliage_gCm2day = quantile(rowMeans(states_all$combined_alloc_foliage_gCm2day, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_combined_alloc_foliage_gCm2day = apply(t(apply(states_all$combined_alloc_foliage_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Assign pool to site_output
      site_output$foliage_gCm2 = apply(states_all$foliage_gCm2,2,quantile,prob=num_quantiles, na.rm = na_flag)
      site_output$mean_foliage_gCm2 = quantile(rowMeans(states_all$foliage_gCm2, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_foliage_gCm2 = apply(t(apply(states_all$foliage_gCm2,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Determine net pool change over time
      dCbio = states_all$foliage_gCm2 - states_all$foliage_gCm2[,1] # difference in root from initial
      site_output$dCfoliage_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
      # Check for the possible pathways
      if (any(check_list == "alloc_foliage_gCm2day")) {
          site_output$alloc_foliage_gCm2day = apply(states_all$alloc_foliage_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_alloc_foliage_gCm2day = quantile(rowMeans(states_all$alloc_foliage_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_alloc_foliage_gCm2day = apply(t(apply(states_all$alloc_foliage_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      }
      if (any(check_list == "foliage_to_litter_gCm2day")) {
          site_output$foliage_to_litter_gCm2day = apply(states_all$foliage_to_litter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_foliage_to_litter_gCm2day = quantile(rowMeans(states_all$foliage_to_litter_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_foliage_to_litter_gCm2day = apply(t(apply(states_all$foliage_to_litter_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$NaturalFractionOfTurnover_foliage = states_all$foliage_to_litter_gCm2day
          # Begin accumulating the total output fluxes here
          site_output$outflux_foliage_gCm2day = states_all$foliage_to_litter_gCm2day
          # Accumulate the combined total of pool lost as litter,
          # NOTE depending on DALEC version may go to either a litter pool or som
          site_output$combined_foliage_to_litter_gCm2day = states_all$foliage_to_litter_gCm2day
      }
      if (any(check_list == "FIREemiss_foliage_gCm2day")) {
          site_output$FIREemiss_foliage_gCm2day = apply(states_all$FIREemiss_foliage_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_FIREemiss_foliage_gCm2day = quantile(rowMeans(states_all$FIREemiss_foliage_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_FIREemiss_foliage_gCm2day = apply(t(apply(states_all$FIREemiss_foliage_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$FireFractionOfTurnover_foliage = states_all$FIREemiss_foliage_gCm2day
          site_output$outflux_foliage_gCm2day = site_output$outflux_foliage_gCm2day + states_all$FIREemiss_foliage_gCm2day
      }
      if (any(check_list == "FIRElitter_foliage_gCm2day")) {
          site_output$FIRElitter_foliage_gCm2day = apply(states_all$FIRElitter_foliage_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_FIRElitter_foliage_gCm2day = quantile(rowMeans(states_all$FIRElitter_foliage_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_FIRElitter_foliage_gCm2day = apply(t(apply(states_all$FIRElitter_foliage_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$outflux_foliage_gCm2day = site_output$outflux_foliage_gCm2day + states_all$FIRElitter_foliage_gCm2day
          site_output$FireFractionOfTurnover_foliage = site_output$FireFractionOfTurnover_foliage + states_all$FIRElitter_foliage_gCm2day
          # Accumulate the combined total of pool lost as litter,
          # NOTE depending on DALEC version may go to either a litter pool or som
          site_output$combined_foliage_to_litter_gCm2day = site_output$combined_foliage_to_litter_gCm2day + states_all$FIRElitter_foliage_gCm2day
      }
      if (any(check_list == "HARVESTextracted_foliage_gCm2day")) {
          site_output$HARVESTextracted_foliage_gCm2day = apply(states_all$HARVESTextracted_foliage_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_HARVESTextracted_foliage_gCm2day = quantile(rowMeans(states_all$HARVESTextracted_foliage_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_HARVESTextracted_foliage_gCm2day = apply(t(apply(states_all$HARVESTextracted_foliage_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$HARVESTlitter_foliage_gCm2day = apply(states_all$HARVESTlitter_foliage_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_HARVESTlitter_foliage_gCm2day = quantile(rowMeans(states_all$HARVESTlitter_foliage_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_HARVESTlitter_foliage_gCm2day = apply(t(apply(states_all$HARVESTlitter_foliage_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$outflux_foliage_gCm2day = site_output$outflux_foliage_gCm2day + states_all$HARVESTextracted_foliage_gCm2day + states_all$HARVESTlitter_foliage_gCm2day
          site_output$HarvestFractionOfTurnover_foliage = states_all$HARVESTextracted_foliage_gCm2day + states_all$HARVESTlitter_foliage_gCm2day
          # Accumulate the combined total of pool lost as litter,
          # NOTE depending on DALEC version may go to either a litter pool or som
          site_output$combined_foliage_to_litter_gCm2day = site_output$combined_foliage_to_litter_gCm2day + states_all$HARVESTlitter_foliage_gCm2day
      }
      # Use this information to determine the mean residence times as it evolves over time
      # NOTE: rollapply inverts the dimensions from that wanted
      site_output$MTT_annual_foliage_years = t( apply(states_all$foliage_gCm2,1, rollapply_mean_annual, step = steps_per_year)
                                              / (apply(site_output$outflux_foliage_gCm2day,1, rollapply_mean_annual, step = steps_per_year) * 365.25) )
      site_output$MTT_annual_foliage_years = apply(site_output$MTT_annual_foliage_years,2,quantile,prob=num_quantiles, na.rm = na_flag)
      # Convert to fractions
      tmp = rowMeans(site_output$outflux_foliage_gCm2day, na.rm = na_flag)
      site_output$NaturalFractionOfTurnover_foliage = quantile(rowMeans(site_output$NaturalFractionOfTurnover_foliage, na.rm = na_flag) / tmp,
                                                               prob = num_quantiles, na.rm = na_flag)
      site_output$FireFractionOfTurnover_foliage = quantile(rowMeans(site_output$FireFractionOfTurnover_foliage, na.rm = na_flag) / tmp,
                                                            prob = num_quantiles, na.rm = na_flag)
      site_output$HarvestFractionOfTurnover_foliage = quantile(rowMeans(site_output$HarvestFractionOfTurnover_foliage, na.rm = na_flag) / tmp,
                                                               prob = num_quantiles, na.rm = na_flag)
      # Mean outflux from foliage
      site_output$mean_outflux_foliage_gCm2day = quantile(rowMeans(site_output$outflux_foliage_gCm2day, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_combined_foliage_to_litter_gCm2day = quantile(rowMeans(site_output$combined_foliage_to_litter_gCm2day, na.rm = na_flag), prob=num_quantiles)
      # Mean annual outflux from foliage
      site_output$mean_annual_outflux_foliage_gCm2day = apply(t(apply(site_output$outflux_foliage_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      site_output$mean_annual_combined_foliage_to_litter_gCm2day = apply(t(apply(site_output$combined_foliage_to_litter_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Now aggregate across quantiles
      site_output$outflux_foliage_gCm2day = apply(site_output$outflux_foliage_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
      site_output$combined_foliage_to_litter_gCm2day = apply(site_output$combined_foliage_to_litter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
  } # exists("foliage")

  # Fine roots related pool, change, input and output variables
  if (any(check_list == "roots_gCm2")) {
      # Assign pool to site_output
      site_output$roots_gCm2 = apply(states_all$roots_gCm2,2,quantile,prob=num_quantiles, na.rm = na_flag)
      site_output$mean_roots_gCm2 = quantile(rowMeans(states_all$roots_gCm2, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_roots_gCm2 = apply(t(apply(states_all$roots_gCm2,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Determine net pool change over time
      dCbio = states_all$roots_gCm2 - states_all$roots_gCm2[,1] # difference in root from initial
      site_output$dCroots_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
      # Calculate the mean annual maximums
      dCbio = apply(states_all$roots_gCm2, 1, rollapply_mean_annual_max, step = steps_per_year)
      site_output$annual_max_roots_gCm2 = quantile(dCbio, prob=num_quantiles, na.rm = na_flag)
      # Declare combined natural, fire and harvest driven creation of litter
      site_output$combined_roots_to_litter_gCm2day = array(NA, dim=dim(states_all$roots_gCm2))
      # Is rooting depth calculate (m) by this model?
      if (any(check_list == "RootDepth_m")) {
          site_output$RootDepth_m = apply(states_all$RootDepth_m,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_RootDepth_m = quantile(rowMeans(states_all$RootDepth_m, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_RootDepth_m = apply(t(apply(states_all$RootDepth_m,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          dCbio = states_all$RootDepth_m - states_all$RootDepth_m[,1] # difference in root from initial
          site_output$dRootDepth_m = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
      }
      # Check for the possible pathways
      if (any(check_list == "alloc_roots_gCm2day")) {
          site_output$alloc_roots_gCm2day = apply(states_all$alloc_roots_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_alloc_roots_gCm2day = quantile(rowMeans(states_all$alloc_roots_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_alloc_roots_gCm2day = apply(t(apply(states_all$alloc_roots_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      }
      if (any(check_list == "roots_to_litter_gCm2day")) {
          site_output$roots_to_litter_gCm2day = apply(states_all$roots_to_litter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_roots_to_litter_gCm2day = quantile(rowMeans(states_all$roots_to_litter_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_roots_to_litter_gCm2day = apply(t(apply(states_all$roots_to_litter_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$NaturalFractionOfTurnover_roots = states_all$roots_to_litter_gCm2day
          # Begin accumulation of output fluxes
          site_output$outflux_roots_gCm2day = states_all$roots_to_litter_gCm2day
          # Accumulate the combined total of pool lost as litter,
          # NOTE depending on DALEC version may go to either a litter pool or som
          site_output$combined_roots_to_litter_gCm2day = states_all$roots_to_litter_gCm2day
      }
      # If this one exists then maybe some other fluxes do
      if (any(check_list == "FIREemiss_roots_gCm2day")) {
          site_output$FIREemiss_roots_gCm2day = apply(states_all$FIREemiss_roots_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_FIREemiss_roots_gCm2day = quantile(rowMeans(states_all$FIREemiss_roots_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_FIREemiss_roots_gCm2day = apply(t(apply(states_all$FIREemiss_roots_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$FireFractionOfTurnover_roots = states_all$FIREemiss_roots_gCm2day
          site_output$outflux_roots_gCm2day = site_output$outflux_roots_gCm2day + states_all$FIREemiss_roots_gCm2day
      }
      if (any(check_list == "FIRElitter_roots_gCm2day")) {
          site_output$FIRElitter_roots_gCm2day = apply(states_all$FIRElitter_roots_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_FIRElitter_roots_gCm2day = quantile(rowMeans(states_all$FIRElitter_roots_gCm2day,1, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_FIRElitter_roots_gCm2day = apply(t(apply(states_all$FIRElitter_roots_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$outflux_roots_gCm2day = site_output$outflux_roots_gCm2day + states_all$FIRElitter_roots_gCm2day
          site_output$FireFractionOfTurnover_roots = site_output$FireFractionOfTurnover_roots + states_all$FIRElitter_roots_gCm2day
          # Accumulate the combined total of pool lost as litter,
          # NOTE depending on DALEC version may go to either a litter pool or som
          site_output$combined_roots_to_litter_gCm2day = site_output$combined_roots_to_litter_gCm2day + states_all$FIRElitter_roots_gCm2day
      }
      if (any(check_list == "HARVESTextracted_roots_gCm2day")) {
          site_output$HARVESTextracted_roots_gCm2day = apply(states_all$HARVESTextracted_roots_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_HARVESTextracted_roots_gCm2day = quantile(rowMeans(states_all$HARVESTextracted_roots_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_HARVESTextracted_roots_gCm2day = apply(t(apply(states_all$HARVESTextracted_roots_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$HARVESTlitter_roots_gCm2day = apply(states_all$HARVESTlitter_roots_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_HARVESTlitter_roots_gCm2day = quantile(rowMeans(states_all$HARVESTlitter_roots_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_HARVESTlitter_roots_gCm2day = apply(t(apply(states_all$HARVESTlitter_roots_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$outflux_roots_gCm2day = site_output$outflux_roots_gCm2day + states_all$HARVESTextracted_roots_gCm2day + states_all$HARVESTlitter_roots_gCm2day
          site_output$HarvestFractionOfTurnover_roots = states_all$HARVESTextracted_roots_gCm2day + states_all$HARVESTlitter_roots_gCm2day
          # Accumulate the combined total of pool lost as litter,
          # NOTE depending on DALEC version may go to either a litter pool or som
          site_output$combined_roots_to_litter_gCm2day = site_output$combined_roots_to_litter_gCm2day + states_all$HARVESTlitter_roots_gCm2day
      }
      # Use this information to determine the mean residence times as it evolves over time
      # NOTE: rollapply inverts the dimensions from that wanted
      site_output$MTT_annual_roots_years = t( apply(states_all$roots_gCm2,1, rollapply_mean_annual, step = steps_per_year)
                                            / (apply(site_output$outflux_roots_gCm2day,1, rollapply_mean_annual, step = steps_per_year) * 365.25) )
      site_output$MTT_annual_roots_years = apply(site_output$MTT_annual_roots_years,2,quantile,prob=num_quantiles, na.rm = na_flag)
      # Convert to fractions
      tmp = rowMeans(site_output$outflux_roots_gCm2day, na.rm = na_flag)
      site_output$NaturalFractionOfTurnover_roots = quantile(rowMeans(site_output$NaturalFractionOfTurnover_roots, na.rm = na_flag) / tmp,
                                                             prob = num_quantiles, na.rm = na_flag)
      site_output$FireFractionOfTurnover_roots = quantile(rowMeans(site_output$FireFractionOfTurnover_roots, na.rm = na_flag) / tmp,
                                                          prob = num_quantiles, na.rm = na_flag)
      site_output$HarvestFractionOfTurnover_roots = quantile(rowMeans(site_output$HarvestFractionOfTurnover_roots, na.rm = na_flag) / tmp,
                                                             prob = num_quantiles, na.rm = na_flag)
      # Mean outflux from roots
      site_output$mean_outflux_roots_gCm2day = quantile(rowMeans(site_output$outflux_roots_gCm2day, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_combined_roots_to_litter_gCm2day = quantile(rowMeans(site_output$combined_roots_to_litter_gCm2day, na.rm = na_flag), prob=num_quantiles)
      # Mean annual outflux from roots
      site_output$mean_annual_outflux_roots_gCm2day = apply(t(apply(site_output$outflux_roots_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      site_output$mean_annual_combined_roots_to_litter_gCm2day = apply(t(apply(site_output$combined_roots_to_litter_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Now aggregate across quantiles
      site_output$outflux_roots_gCm2day = apply(site_output$outflux_roots_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
      site_output$combined_roots_to_litter_gCm2day = apply(site_output$combined_roots_to_litter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
  } # exists("roots")

  # Wood related pool, change, input and output variables
  if (any(check_list == "wood_gCm2")) {
      # Assign pool to site_output
      site_output$wood_gCm2 = apply(states_all$wood_gCm2,2,quantile,prob=num_quantiles, na.rm = na_flag)
      site_output$mean_wood_gCm2 = quantile(rowMeans(states_all$wood_gCm2, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_wood_gCm2 = apply(t(apply(states_all$wood_gCm2,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Determine net pool change over time
      dCbio = states_all$wood_gCm2 - states_all$wood_gCm2[,1] # difference in wood from initial
      site_output$dCwood_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
      # Calculate the mean annual maximums
      dCbio = apply(states_all$wood_gCm2, 1, rollapply_mean_annual_max, step = steps_per_year)
      site_output$annual_max_wood_gCm2 = quantile(dCbio, prob=num_quantiles, na.rm = na_flag)
      # Declare combined natural, fire and harvest driven creation of litter
      site_output$combined_wood_to_litter_gCm2day = array(NA, dim=dim(states_all$wood_gCm2))
      # Check for the possible pathways
      if (any(check_list == "alloc_wood_gCm2day")) {
          site_output$alloc_wood_gCm2day = apply(states_all$alloc_wood_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_alloc_wood_gCm2day = quantile(rowMeans(states_all$alloc_wood_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_alloc_wood_gCm2day = apply(t(apply(states_all$alloc_wood_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      }
      if (any(check_list == "wood_to_litter_gCm2day")) {
          site_output$wood_to_litter_gCm2day = apply(states_all$wood_to_litter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_wood_to_litter_gCm2day = quantile(rowMeans(states_all$wood_to_litter_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_wood_to_litter_gCm2day = apply(t(apply(states_all$wood_to_litter_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$NaturalFractionOfTurnover_wood = states_all$wood_to_litter_gCm2day
          # Begin accumulating output fluxes
          site_output$outflux_wood_gCm2day = states_all$wood_to_litter_gCm2day
          # Accumulate the combined total of pool lost as litter,
          # NOTE depending on DALEC version may go to either a litter pool or som
          site_output$combined_wood_to_litter_gCm2day = states_all$wood_to_litter_gCm2day
      }
      # If this one exists then maybe some other fluxes do
      if (any(check_list == "FIREemiss_wood_gCm2day")) {
          site_output$FIREemiss_wood_gCm2day = apply(states_all$FIREemiss_wood_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_FIREemiss_wood_gCm2day = quantile(rowMeans(states_all$FIREemiss_wood_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_FIREemiss_wood_gCm2day = apply(t(apply(states_all$FIREemiss_wood_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$FireFractionOfTurnover_wood = states_all$FIREemiss_wood_gCm2day
          site_output$outflux_wood_gCm2day = site_output$outflux_wood_gCm2day + states_all$FIREemiss_wood_gCm2day
      }
      if (any(check_list == "FIRElitter_wood_gCm2day")) {
          site_output$FIRElitter_wood_gCm2day = apply(states_all$FIRElitter_wood_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_FIRElitter_wood_gCm2day = quantile(rowMeans(states_all$FIRElitter_wood_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_FIRElitter_wood_gCm2day = apply(t(apply(states_all$FIRElitter_wood_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$outflux_wood_gCm2day = site_output$outflux_wood_gCm2day + states_all$FIRElitter_wood_gCm2day
          site_output$FireFractionOfTurnover_wood = site_output$FireFractionOfTurnover_wood + states_all$FIRElitter_wood_gCm2day
          # Accumulate the combined total of pool lost as litter,
          # NOTE depending on DALEC version may go to either a litter pool or som
          site_output$combined_wood_to_litter_gCm2day = site_output$combined_wood_to_litter_gCm2day + states_all$FIRElitter_wood_gCm2day
      }
      if (any(check_list == "HARVESTextracted_wood_gCm2day")) {
          site_output$HARVESTextracted_wood_gCm2day = apply(states_all$HARVESTextracted_wood_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_HARVESTextracted_wood_gCm2day = quantile(rowMeans(states_all$HARVESTextracted_wood_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_HARVESTextracted_wood_gCm2day = apply(t(apply(states_all$HARVESTextracted_wood_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$HARVESTlitter_wood_gCm2day = apply(states_all$HARVESTlitter_wood_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_HARVESTlitter_wood_gCm2day = quantile(rowMeans(states_all$HARVESTlitter_wood_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_HARVESTlitter_wood_gCm2day = apply(t(apply(states_all$HARVESTlitter_wood_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$outflux_wood_gCm2day = site_output$outflux_wood_gCm2day + states_all$HARVESTextracted_wood_gCm2day + states_all$HARVESTlitter_wood_gCm2day
          site_output$HarvestFractionOfTurnover_wood = states_all$HARVESTextracted_wood_gCm2day + states_all$HARVESTlitter_wood_gCm2day
          # Accumulate the combined total of pool lost as litter,
          # NOTE depending on DALEC version may go to either a litter pool or som
          site_output$combined_wood_to_litter_gCm2day = site_output$combined_wood_to_litter_gCm2day + states_all$HARVESTlitter_wood_gCm2day
      }
      # Use this information to determine the mean residence times as it evolves over time
      # NOTE: rollapply inverts the dimensions from that wanted
      site_output$MTT_annual_wood_years = t( apply(states_all$wood_gCm2,1, rollapply_mean_annual, step = steps_per_year)
                                           / (apply(site_output$outflux_wood_gCm2day,1, rollapply_mean_annual, step = steps_per_year) * 365.25) )
      site_output$MTT_annual_wood_years = apply(site_output$MTT_annual_wood_years,2,quantile,prob=num_quantiles, na.rm = na_flag)
      # Convert to fractions
      tmp = rowMeans(site_output$outflux_wood_gCm2day, na.rm = na_flag)
      site_output$NaturalFractionOfTurnover_wood = quantile(rowMeans(site_output$NaturalFractionOfTurnover_wood, na.rm = na_flag) / tmp,
                                                            prob = num_quantiles, na.rm = na_flag)
      site_output$FireFractionOfTurnover_wood = quantile(rowMeans(site_output$FireFractionOfTurnover_wood, na.rm = na_flag) / tmp,
                                                         prob = num_quantiles, na.rm = na_flag)
      site_output$HarvestFractionOfTurnover_wood = quantile(rowMeans(site_output$HarvestFractionOfTurnover_wood, na.rm = na_flag) / tmp,
                                                            prob = num_quantiles, na.rm = na_flag)
      # Mean outflux from wood
      site_output$mean_outflux_wood_gCm2day = quantile(rowMeans(site_output$outflux_wood_gCm2day, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_combined_wood_to_litter_gCm2day = quantile(rowMeans(site_output$combined_wood_to_litter_gCm2day, na.rm = na_flag), prob=num_quantiles)
      # Mean annual outflux from wood
      site_output$mean_annual_outflux_wood_gCm2day = apply(t(apply(site_output$outflux_wood_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      site_output$mean_annual_combined_wood_to_litter_gCm2day = apply(t(apply(site_output$combined_wood_to_litter_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Now aggregate across quantiles
      site_output$outflux_wood_gCm2day = apply(site_output$outflux_wood_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
      site_output$combined_wood_to_litter_gCm2day = apply(site_output$combined_wood_to_litter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
  } # exists("wood")

  # Fine litter pool, change and output variables
  if (any(check_list == "litter_gCm2")) {
      # Assign pool to site_output
      site_output$litter_gCm2 = apply(states_all$litter_gCm2,2,quantile,prob=num_quantiles, na.rm = na_flag)
      site_output$mean_litter_gCm2 = quantile(rowMeans(states_all$litter_gCm2, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_litter_gCm2 = apply(t(apply(states_all$litter_gCm2,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Determine net pool change over time
      dCbio = states_all$litter_gCm2 - states_all$litter_gCm2[,1] # difference in litter from initial
      site_output$dClitter_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
      # Accumulate total output
      site_output$outflux_litter_gCm2day = states_all$rhet_litter_gCm2day + states_all$litter_to_som_gCm2day
      site_output$NaturalFractionOfTurnover_litter = site_output$outflux_litter_gCm2day
      # Heterotrophic respiration
      site_output$rhet_litter_gCm2day = apply(states_all$rhet_litter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
      site_output$mean_rhet_litter_gCm2day = quantile(rowMeans(states_all$rhet_litter_gCm2day, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_rhet_litter_gCm2day = apply(t(apply(states_all$rhet_litter_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Decomposition
      site_output$litter_to_som_gCm2day = apply(states_all$litter_to_som_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
      site_output$mean_litter_to_som_gCm2day = quantile(rowMeans(states_all$litter_to_som_gCm2day, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_litter_to_som_gCm2day = apply(t(apply(states_all$litter_to_som_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Accumulate combined natural, fire and harvest related fluxes to som
      site_output$combined_litter_to_som_gCm2day = states_all$litter_to_som_gCm2day
      # If this one exists then maybe some other fluxes do
      if (any(check_list == "FIREemiss_litter_gCm2day")) {
          site_output$FIREemiss_litter_gCm2day = apply(states_all$FIREemiss_litter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_FIREemiss_litter_gCm2day = quantile(rowMeans(states_all$FIREemiss_litter_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_FIREemiss_litter_gCm2day = apply(t(apply(states_all$FIREemiss_litter_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$FireFractionOfTurnover_litter = states_all$FIREemiss_litter_gCm2day
          site_output$outflux_litter_gCm2day = site_output$outflux_litter_gCm2day + states_all$FIREemiss_litter_gCm2day
      }
      if (any(check_list == "FIRElitter_litter_gCm2day")) {
          site_output$FIRElitter_litter_gCm2day = apply(states_all$FIRElitter_litter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_FIRElitter_litter_gCm2day = quantile(rowMeans(states_all$FIRElitter_litter_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_FIRElitter_litter_gCm2day = apply(t(apply(states_all$FIRElitter_litter_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$outflux_litter_gCm2day = site_output$outflux_litter_gCm2day + states_all$FIRElitter_litter_gCm2day
          site_output$FireFractionOfTurnover_litter = site_output$FireFractionOfTurnover_litter + states_all$FIRElitter_litter_gCm2day
          # Accumulate combined natural, fire and harvest related fluxes to som
          site_output$combined_litter_to_som_gCm2day = site_output$combined_litter_to_som_gCm2day + states_all$FIRElitter_litter_gCm2day
      }
      if (any(check_list == "HARVESTextracted_litter_gCm2day")) {
          site_output$HARVESTextracted_litter_gCm2day = apply(states_all$HARVESTextracted_litter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_HARVESTextracted_litter_gCm2day = quantile(rowMeans(states_all$HARVESTextracted_litter_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_HARVESTextracted_litter_gCm2day = apply(t(apply(states_all$HARVESTextracted_litter_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$outflux_litter_gCm2day = site_output$outflux_litter_gCm2day + states_all$HARVESTextracted_litter_gCm2day
          site_output$HarvestFractionOfTurnover_litter = states_all$HARVESTextracted_litter_gCm2day
      }
      # Use this information to determine the mean residence times as it evolves over time
      # NOTE: rollapply inverts the dimensions from that wanted
      site_output$MTT_annual_litter_years = t( apply(states_all$litter_gCm2,1, rollapply_mean_annual, step = steps_per_year)
                                             / (apply(site_output$outflux_litter_gCm2day,1, rollapply_mean_annual, step = steps_per_year) * 365.25) )
      site_output$MTT_annual_litter_years = apply(site_output$MTT_annual_litter_years,2,quantile,prob=num_quantiles, na.rm = na_flag)
      # Convert to fractions
      tmp = rowMeans(site_output$outflux_litter_gCm2day, na.rm = na_flag)
      site_output$NaturalFractionOfTurnover_litter = quantile(rowMeans(site_output$NaturalFractionOfTurnover_litter, na.rm = na_flag) / tmp,
                                                              prob = num_quantiles, na.rm = na_flag)
      site_output$FireFractionOfTurnover_litter = quantile(rowMeans(site_output$FireFractionOfTurnover_litter, na.rm = na_flag) / tmp,
                                                           prob = num_quantiles, na.rm = na_flag)
      site_output$HarvestFractionOfTurnover_litter = quantile(rowMeans(site_output$HarvestFractionOfTurnover_litter, na.rm = na_flag) / tmp,
                                                              prob = num_quantiles, na.rm = na_flag)
      # Mean outflux from litter
      site_output$mean_outflux_litter_gCm2day = quantile(rowMeans(site_output$outflux_litter_gCm2day, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_combined_litter_to_som_gCm2day = quantile(rowMeans(site_output$combined_litter_to_som_gCm2day, na.rm = na_flag), prob=num_quantiles)
      # Mean annual outflux from litter
      site_output$mean_annual_outflux_litter_gCm2day = apply(t(apply(site_output$outflux_litter_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      site_output$mean_annual_combined_litter_to_som_gCm2day = apply(t(apply(site_output$combined_litter_to_som_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Now aggregate across quantiles
      site_output$outflux_litter_gCm2day = apply(site_output$outflux_litter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
      site_output$combined_litter_to_som_gCm2day = apply(site_output$combined_litter_to_som_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
  } # exists("litter")

  # Wood litter pool, change and output variables
  if (any(check_list == "woodlitter_gCm2")) {
      # Assign pool to site_output
      site_output$woodlitter_gCm2 = apply(states_all$woodlitter_gCm2,2,quantile,prob=num_quantiles, na.rm = na_flag)
      site_output$mean_woodlitter_gCm2 = quantile(rowMeans(states_all$woodlitter_gCm2, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_woodlitter_gCm2 = apply(t(apply(states_all$woodlitter_gCm2,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Determine net pool change over time
      dCbio = states_all$woodlitter_gCm2 - states_all$woodlitter_gCm2[,1] # difference in wood litter from initial
      site_output$dCwoodlitter_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
      # Begin accumulating losses from wood litter
      site_output$outflux_woodlitter_gCm2day = states_all$rhet_woodlitter_gCm2day + states_all$woodlitter_to_som_gCm2day
      site_output$NaturalFractionOfTurnover_woodlitter = site_output$outflux_woodlitter_gCm2day
      # Heterotrophic respiration
      site_output$rhet_woodlitter_gCm2day = apply(states_all$rhet_woodlitter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
      site_output$mean_rhet_woodlitter_gCm2day = quantile(rowMeans(states_all$rhet_woodlitter_gCm2day, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_rhet_woodlitter_gCm2day = apply(t(apply(states_all$rhet_woodlitter_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Decomposition
      site_output$woodlitter_to_som_gCm2day = apply(states_all$woodlitter_to_som_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
      site_output$mean_woodlitter_to_som_gCm2day = quantile(rowMeans(states_all$woodlitter_to_som_gCm2day, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_woodlitter_to_som_gCm2day = apply(t(apply(states_all$woodlitter_to_som_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Accumulate combined natural, fire and harvest related fluxes to som
      site_output$combined_woodlitter_to_som_gCm2day = states_all$woodlitter_to_som_gCm2day
      # If this one exists then maybe some other fluxes do
      if (any(check_list == "FIREemiss_woodlitter_gCm2day")) {
          site_output$FIREemiss_woodlitter_gCm2day = apply(states_all$FIREemiss_woodlitter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_FIREemiss_woodlitter_gCm2day = quantile(rowMeans(states_all$FIREemiss_woodlitter_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_FIREemiss_woodlitter_gCm2day = apply(t(apply(states_all$FIREemiss_woodlitter_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$FireFractionOfTurnover_woodlitter = states_all$FIREemiss_woodlitter_gCm2day
          site_output$outflux_woodlitter_gCm2day = site_output$outflux_woodlitter_gCm2day + states_all$FIREemiss_woodlitter_gCm2day
      }
      if (any(check_list == "FIRElitter_woodlitter_gCm2day")) {
          site_output$FIRElitter_woodlitter_gCm2day = apply(states_all$FIRElitter_woodlitter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_FIRElitter_woodlitter_gCm2day = quantile(rowMeans(states_all$FIRElitter_woodlitter_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_FIRElitter_woodlitter_gCm2day = apply(t(apply(states_all$FIRElitter_woodlitter_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$outflux_woodlitter_gCm2day = site_output$outflux_woodlitter_gCm2day + states_all$FIRElitter_woodlitter_gCm2day
          site_output$FireFractionOfTurnover_woodlitter = site_output$FireFractionOfTurnover_woodlitter + states_all$FIRElitter_woodlitter_gCm2day
          # Accumulate combined natural, fire and harvest related fluxes to som
          site_output$combined_woodlitter_to_som_gCm2day = site_output$combined_woodlitter_to_som_gCm2day + states_all$FIRElitter_woodlitter_gCm2day
      }
      if (any(check_list == "HARVESTextracted_woodlitter_gCm2day")) {
          site_output$HARVESTextracted_woodlitter_gCm2day = apply(states_all$HARVESTextracted_woodlitter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_HARVESTextracted_woodlitter_gCm2day = quantile(rowMeans(states_all$HARVESTextracted_woodlitter_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_HARVESTextracted_woodlitter_gCm2day = apply(t(apply(states_all$HARVESTextracted_woodlitter_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$outflux_woodlitter_gCm2day = site_output$outflux_woodlitter_gCm2day + states_all$HARVESTextracted_woodlitter_gCm2day
          site_output$HarvestFractionOfTurnover_woodlitter = states_all$HARVESTextracted_woodlitter_gCm2day
      }
      # Use this information to determine the mean residence times as it evolves over time
      # NOTE: rollapply inverts the dimensions from that wanted
      site_output$MTT_annual_woodlitter_years = t( apply(states_all$woodlitter_gCm2,1, rollapply_mean_annual, step = steps_per_year)
                                                 / (apply(site_output$outflux_woodlitter_gCm2day,1, rollapply_mean_annual, step = steps_per_year) * 365.25) )
      site_output$MTT_annual_woodlitter_years = apply(site_output$MTT_annual_woodlitter_years,2,quantile,prob=num_quantiles, na.rm = na_flag)
      # Convert to fractions
      tmp = rowMeans(site_output$outflux_woodlitter_gCm2day, na.rm = na_flag)
      site_output$NaturalFractionOfTurnover_woodlitter = quantile(rowMeans(site_output$NaturalFractionOfTurnover_woodlitter, na.rm = na_flag) / tmp,
                                                                  prob = num_quantiles, na.rm = na_flag)
      site_output$FireFractionOfTurnover_woodlitter = quantile(rowMeans(site_output$FireFractionOfTurnover_woodlitter, na.rm = na_flag) / tmp,
                                                               prob = num_quantiles, na.rm = na_flag)
      site_output$HarvestFractionOfTurnover_woodlitter = quantile(rowMeans(site_output$HarvestFractionOfTurnover_woodlitter, na.rm = na_flag) / tmp,
                                                                  prob = num_quantiles, na.rm = na_flag)
      # Mean outflux from wood
      site_output$mean_outflux_woodlitter_gCm2day = quantile(rowMeans(site_output$outflux_woodlitter_gCm2day, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_combined_woodlitter_to_som_gCm2day = quantile(rowMeans(site_output$combined_woodlitter_to_som_gCm2day, na.rm = na_flag), prob=num_quantiles)
      # Mean annual outflux from wood
      site_output$mean_annual_outflux_woodlitter_gCm2day = apply(t(apply(site_output$outflux_woodlitter_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      site_output$mean_annual_combined_woodlitter_to_som_gCm2day = apply(t(apply(site_output$combined_woodlitter_to_som_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Now aggregate across quantiles
      site_output$outflux_woodlitter_gCm2day = apply(site_output$outflux_woodlitter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
      site_output$combined_woodlitter_to_som_gCm2day = apply(site_output$combined_woodlitter_to_som_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
  } # exists("woodlitter")

  # Soil organic matter pool, change and output variables
  if (any(check_list == "som_gCm2")) {
      # Assign pool to site_output
      site_output$som_gCm2 = apply(states_all$som_gCm2,2,quantile,prob=num_quantiles, na.rm = na_flag)
      site_output$mean_som_gCm2 = quantile(rowMeans(states_all$som_gCm2, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_som_gCm2 = apply(t(apply(states_all$som_gCm2,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Determine net pool change over time
      dCbio = states_all$som_gCm2 - states_all$som_gCm2[,1] # difference in som from initial
      site_output$dCsom_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
      # Begin accumulating losses from som
      site_output$outflux_som_gCm2day = states_all$rhet_som_gCm2day
      site_output$NaturalFractionOfTurnover_som = site_output$outflux_som_gCm2day
      # Heterotrphic respiration
      site_output$rhet_som_gCm2day = apply(states_all$rhet_som_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
      site_output$mean_rhet_som_gCm2day = quantile(rowMeans(states_all$rhet_som_gCm2day, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_rhet_som_gCm2day = apply(t(apply(states_all$rhet_som_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # If this one exists then maybe some other fluxes do
      if (any(check_list == "FIREemiss_som_gCm2day")) {
          site_output$FIREemiss_som_gCm2day = apply(states_all$FIREemiss_som_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_FIREemiss_som_gCm2day = quantile(rowMeans(states_all$FIREemiss_som_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_FIREemiss_som_gCm2day = apply(t(apply(states_all$FIREemiss_som_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$FireFractionOfTurnover_som = states_all$FIREemiss_som_gCm2day
          site_output$outflux_som_gCm2day = site_output$outflux_som_gCm2day + states_all$FIREemiss_som_gCm2day
      }
      # NOTE: FIRElitter does not exist as there is not litter which leaves the som pool
      if (any(check_list == "HARVESTextracted_som_gCm2day")) {
          site_output$HARVESTextracted_som_gCm2day = apply(states_all$HARVESTextracted_som_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_HARVESTextracted_som_gCm2day = quantile(rowMeans(states_all$HARVESTextracted_som_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_HARVESTextracted_som_gCm2day = apply(t(apply(states_all$HARVESTextracted_som_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$outflux_som_gCm2day = site_output$outflux_som_gCm2day + states_all$HARVESTextracted_som_gCm2day
          site_output$HarvestFractionOfTurnover_som = states_all$HARVESTextracted_som_gCm2day
      }
      # Use this information to determine the mean residence times as it evolves over time
      # NOTE: rollapply inverts the dimensions from that wanted
      site_output$MTT_annual_som_years = t( apply(states_all$som_gCm2,1, rollapply_mean_annual, step = steps_per_year)
                                          / (apply(site_output$outflux_som_gCm2day,1, rollapply_mean_annual, step = steps_per_year) * 365.25) )
      site_output$MTT_annual_som_years = apply(site_output$MTT_annual_som_years,2,quantile,prob=num_quantiles, na.rm = na_flag)
      # Convert to fractions
      tmp = rowMeans(site_output$outflux_som_gCm2day, na.rm = na_flag)
      site_output$NaturalFractionOfTurnover_som = quantile(rowMeans(site_output$NaturalFractionOfTurnover_som, na.rm = na_flag) / tmp,
                                                           prob = num_quantiles, na.rm = na_flag)
      site_output$FireFractionOfTurnover_som = quantile(rowMeans(site_output$FireFractionOfTurnover_som, na.rm = na_flag) / tmp,
                                                        prob = num_quantiles, na.rm = na_flag)
      site_output$HarvestFractionOfTurnover_som = quantile(rowMeans(site_output$HarvestFractionOfTurnover_som, na.rm = na_flag) / tmp,
                                                           prob = num_quantiles, na.rm = na_flag)
      # Mean outflux from som
      site_output$mean_outflux_som_gCm2day = quantile(rowMeans(site_output$outflux_som_gCm2day, na.rm = na_flag), prob=num_quantiles)
      # Mean annual outflux from som
      site_output$mean_annual_outflux_som_gCm2day = apply(t(apply(site_output$outflux_som_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Now aggregate across quantiles
      site_output$outflux_som_gCm2day = apply(site_output$outflux_som_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
  } # exists("som")

  # Biomass related pool, change and output variables
  if (any(check_list == "biomass_gCm2")) {
      # Assign pool to site_output
      site_output$biomass_gCm2 = apply(states_all$biomass_gCm2,2,quantile,prob=num_quantiles, na.rm = na_flag)
      site_output$mean_biomass_gCm2 = quantile(rowMeans(states_all$biomass_gCm2, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_biomass_gCm2 = apply(t(apply(states_all$biomass_gCm2,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Determine net pool change over time
      dCbio = states_all$biomass_gCm2 - states_all$biomass_gCm2[,1] # difference in labile from initial
      site_output$dCbiomass_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
      # Track natural biomass losses
      site_output$biomass_to_litter_gCm2day = apply(states_all$biomass_to_litter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
      site_output$mean_biomass_to_litter_gCm2day = quantile(rowMeans(states_all$biomass_to_litter_gCm2day, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_biomass_to_litter_gCm2day = apply(t(apply(states_all$biomass_to_litter_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Begin accumulating the total output
      site_output$outflux_biomass_gCm2day = states_all$biomass_to_litter_gCm2day
      site_output$NaturalFractionOfTurnover_biomass = states_all$biomass_to_litter_gCm2day
      # Begin accumulating the combined natural, fire and harvest litter creation
      site_output$combined_biomass_to_litter_gCm2day = states_all$biomass_to_litter_gCm2day
      # Other natural flux pathways should really go here before disturbance related
      if (any(check_list == "FIREemiss_biomass_gCm2day")) {
          site_output$FIREemiss_biomass_gCm2day = apply(states_all$FIREemiss_biomass_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_FIREemiss_biomass_gCm2day = quantile(rowMeans(states_all$FIREemiss_biomass_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_FIREemiss_biomass_gCm2day = apply(t(apply(states_all$FIREemiss_biomass_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$outflux_biomass_gCm2day = site_output$outflux_biomass_gCm2day + states_all$FIREemiss_biomass_gCm2day
          site_output$FireFractionOfTurnover_biomass = states_all$FIREemiss_biomass_gCm2day
      }
      if (any(check_list == "FIRElitter_biomass_gCm2day")) {
          site_output$FIRElitter_biomass_gCm2day = apply(states_all$FIRElitter_biomass_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_FIRElitter_biomass_gCm2day = quantile(rowMeans(states_all$FIRElitter_biomass_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_FIRElitter_biomass_gCm2day = apply(t(apply(states_all$FIRElitter_biomass_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$outflux_biomass_gCm2day = site_output$outflux_biomass_gCm2day + states_all$FIRElitter_biomass_gCm2day
          site_output$FireFractionOfTurnover_biomass = site_output$FireFractionOfTurnover_biomass + states_all$FIRElitter_biomass_gCm2day
          # Accumulte the combined fire, harvest and natural biomass flux to litter
          site_output$combined_biomass_to_litter_gCm2day = site_output$combined_biomass_to_litter_gCm2day + states_all$FIRElitter_biomass_gCm2day
      }
      if (any(check_list == "HARVESTextracted_biomass_gCm2day")) {
          site_output$HARVESTextracted_biomass_gCm2day = apply(states_all$HARVESTextracted_biomass_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_HARVESTextracted_biomass_gCm2day = quantile(rowMeans(states_all$HARVESTextracted_biomass_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_HARVESTextracted_biomass_gCm2day = apply(t(apply(states_all$HARVESTextracted_biomass_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$HARVESTlitter_biomass_gCm2day = apply(states_all$HARVESTlitter_biomass_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_HARVESTlitter_biomass_gCm2day = quantile(rowMeans(states_all$HARVESTlitter_biomass_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_HARVESTlitter_biomass_gCm2day = apply(t(apply(states_all$HARVESTlitter_biomass_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$outflux_biomass_gCm2day = site_output$outflux_biomass_gCm2day + states_all$HARVESTextracted_biomass_gCm2day + states_all$HARVESTlitter_biomass_gCm2day
          site_output$HarvestFractionOfTurnover_biomass = states_all$HARVESTextracted_biomass_gCm2day + states_all$HARVESTlitter_biomass_gCm2day
      }
      # Estimate the ecosystem mean transit (residence) times as a function of natural and disturbance processes
      # NOTE: rollapply inverts the dimensions from that wanted
      site_output$MTT_annual_biomass_years = t( apply(states_all$biomass_gCm2,1, rollapply_mean_annual, step = steps_per_year)
                                              / (apply(site_output$outflux_biomass_gCm2day,1, rollapply_mean_annual, step = steps_per_year) * 365.25) )
      site_output$MTT_annual_biomass_years = apply(site_output$MTT_annual_biomass_years,2,quantile,prob=num_quantiles, na.rm = na_flag)
      # Convert to fractions
      tmp = rowMeans(site_output$outflux_biomass_gCm2day, na.rm = na_flag)
      site_output$NaturalFractionOfTurnover_biomass = quantile(rowMeans(site_output$NaturalFractionOfTurnover_biomass, na.rm = na_flag) / tmp,
                                                               prob = num_quantiles, na.rm = na_flag)
      site_output$FireFractionOfTurnover_biomass = quantile(rowMeans(site_output$FireFractionOfTurnover_biomass, na.rm = na_flag) / tmp,
                                                            prob = num_quantiles, na.rm = na_flag)
      site_output$HarvestFractionOfTurnover_biomass = quantile(rowMeans(site_output$HarvestFractionOfTurnover_biomass, na.rm = na_flag) /tmp,
                                                               prob = num_quantiles, na.rm = na_flag)
      # Mean outflux from biomass
      site_output$mean_outflux_biomass_gCm2day = quantile(rowMeans(site_output$outflux_biomass_gCm2day, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_combined_biomass_to_litter_gCm2day = quantile(rowMeans(site_output$combined_biomass_to_litter_gCm2day, na.rm = na_flag), prob=num_quantiles)
      # Mean annual outflux from biomass
      site_output$mean_annual_outflux_biomass_gCm2day = apply(t(apply(site_output$outflux_biomass_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      site_output$mean_annual_combined_biomass_to_litter_gCm2day = apply(t(apply(site_output$combined_biomass_to_litter_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Now aggregate across quantiles
      site_output$outflux_biomass_gCm2day = apply(site_output$outflux_biomass_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
      site_output$combined_biomass_to_litter_gCm2day = apply(site_output$combined_biomass_to_litter_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
  } # exists("biomass_gCm2")

  # Dead organic matter pool, change and output variables
  # NOTE: this is potentially the combination of som, litter and wood litter
  if (any(check_list == "dom_gCm2")) {
      # Assign pool to site_output
      site_output$dom_gCm2 = apply(states_all$dom_gCm2,2,quantile,prob=num_quantiles, na.rm = na_flag)
      site_output$mean_dom_gCm2 = quantile(rowMeans(states_all$dom_gCm2, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_dom_gCm2 = apply(t(apply(states_all$dom_gCm2,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Determine net pool change over time
      dCbio = states_all$dom_gCm2 - states_all$dom_gCm2[,1] # difference in dom from initial
      site_output$dCdom_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
      # Accumulating output fluxes from dead organic matter
      site_output$outflux_dom_gCm2day = states_all$rhet_dom_gCm2day
      site_output$NaturalFractionOfTurnover_dom = site_output$outflux_dom_gCm2day
      # Heterotrophic respiration
      site_output$rhet_dom_gCm2day = apply(states_all$rhet_dom_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
      site_output$mean_rhet_dom_gCm2day = quantile(rowMeans(states_all$rhet_dom_gCm2day, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_rhet_dom_gCm2day = apply(t(apply(states_all$rhet_dom_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # If this one exists then maybe some other fluxes do
      if (any(check_list == "FIREemiss_dom_gCm2day")) {
          site_output$FIREemiss_dom_gCm2day = apply(states_all$FIREemiss_dom_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_FIREemiss_dom_gCm2day = quantile(rowMeans(states_all$FIREemiss_dom_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_FIREemiss_dom_gCm2day = apply(t(apply(states_all$FIREemiss_dom_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$FireFractionOfTurnover_dom = states_all$FIREemiss_dom_gCm2day
          site_output$outflux_dom_gCm2day = site_output$outflux_dom_gCm2day + states_all$FIREemiss_dom_gCm2day
      }
      # NOTE: FIRElitter does not exist as there is not litter which leaves the dom pool
      if (any(check_list == "HARVESTextracted_dom_gCm2day")) {
          site_output$HARVESTextracted_dom_gCm2day = apply(states_all$HARVESTextracted_dom_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
          site_output$mean_HARVESTextracted_dom_gCm2day = quantile(rowMeans(states_all$HARVESTextracted_dom_gCm2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_HARVESTextracted_dom_gCm2day = apply(t(apply(states_all$HARVESTextracted_dom_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2, quantile, prob=num_quantiles, na.rm = TRUE)
          site_output$outflux_dom_gCm2day = site_output$outflux_dom_gCm2day + states_all$HARVESTextracted_som_gCm2day
          site_output$HarvestFractionOfTurnover_dom = states_all$HARVESTextracted_dom_gCm2day
      }
      # Use this information to determine the mean residence times as it evolves over time
      # NOTE: rollapply inverts the dimensions from that wanted
      site_output$MTT_annual_dom_years = t( apply(states_all$dom_gCm2,1, rollapply_mean_annual, step = steps_per_year)
                                          / (apply(site_output$outflux_dom_gCm2day,1, rollapply_mean_annual, step = steps_per_year) * 365.25) )
      site_output$MTT_annual_dom_years = apply(site_output$MTT_annual_dom_years,2,quantile,prob=num_quantiles, na.rm = na_flag)
      # Convert to fractions
      tmp = rowMeans(site_output$outflux_dom_gCm2day, na.rm = na_flag)
      site_output$NaturalFractionOfTurnover_dom = quantile(rowMeans(site_output$NaturalFractionOfTurnover_dom, na.rm = na_flag) / tmp,
                                                           prob = num_quantiles, na.rm = na_flag)
      site_output$FireFractionOfTurnover_dom = quantile(rowMeans(site_output$FireFractionOfTurnover_dom, na.rm = na_flag) / tmp,
                                                        prob = num_quantiles, na.rm = na_flag)
      site_output$HarvestFractionOfTurnover_dom = quantile(rowMeans(site_output$HarvestFractionOfTurnover_dom, na.rm = na_flag) / tmp,
                                                           prob = num_quantiles, na.rm = na_flag)
      # Mean outflux from dom
      site_output$mean_outflux_dom_gCm2day = quantile(rowMeans(site_output$outflux_dom_gCm2day, na.rm = na_flag), prob=num_quantiles)
      # Mean annual outflux from dom
      site_output$mean_annual_outflux_dom_gCm2day = apply(t(apply(site_output$outflux_dom_gCm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Now aggregate across quantiles
      site_output$outflux_dom_gCm2day = apply(site_output$outflux_dom_gCm2day,2,quantile,prob=num_quantiles, na.rm = na_flag)
  }

  ###
  # Aggregate Water cycle related fluxes
  ###

  # Water cycle specific if available
  if (any(check_list == "ET_kgH2Om2day")) {
      ## current water in the soil surface layer (0-30 cm)
      site_output$SurfWater_kgH2Om2 = apply(states_all$SurfWater_kgH2Om2,2,quantile,prob=num_quantiles,na.rm = na_flag)
      site_output$mean_SurfWater_kgH2Om2 = quantile(rowMeans(states_all$SurfWater_kgH2Om2, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_SurfWater_kgH2Om2 = apply(t(apply(states_all$SurfWater_kgH2Om2,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Calculate change over time
      dCbio = states_all$SurfWater_kgH2Om2 - states_all$SurfWater_kgH2Om2[,1] # difference in surface water from initial
      site_output$dSurfWater_kgH2Om2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)

      # plant apparent soil water potential (MPa)
      site_output$wSWP_MPa = apply(states_all$wSWP_MPa,2,quantile,prob=num_quantiles,na.rm = na_flag)
      site_output$mean_wSWP_MPa = quantile(rowMeans(states_all$wSWP_MPa, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_wSWP_MPa = apply(t(apply(states_all$wSWP_MPa,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Calculate change over time
      dCbio = states_all$wSWP_MPa - states_all$wSWP_MPa[,1] # difference from initial
      site_output$dwSWP_MPa = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)

      # evapotranspiration (Etrans + Esoil + Ewetcanopy)
      site_output$ET_kgH2Om2day = apply(states_all$ET_kgH2Om2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
      site_output$mean_ET_kgH2Om2day = quantile(rowMeans(states_all$ET_kgH2Om2day, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_ET_kgH2Om2day = apply(t(apply(states_all$ET_kgH2Om2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)

      # Calculate the ecosystem water use efficiency
      site_output$wue_eco_gCkgH2O = apply(states_all$gpp_gCm2day/states_all$ET_kgH2Om2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
      site_output$mean_wue_eco_gCkgH2O = quantile(rowSums(states_all$gpp_gCm2day)/rowSums(states_all$ET_kgH2Om2day), prob=num_quantiles, na.rm = na_flag)
      site_output$mean_annual_wue_eco_gCkgH2O = apply(t(apply(states_all$gpp_gCm2day,1, rollapply_mean_annual, step = steps_per_year) /
                                                        apply(states_all$ET_kgH2Om2day,1, rollapply_mean_annual, step = steps_per_year)),
                                                      2, quantile, prob=num_quantiles, na.rm = TRUE)

      # Check whether the evaporation components exist
      if (any(check_list == "Etrans_kgH2Om2day")) {
          # Transpiration
          site_output$Etrans_kgH2Om2day = apply(states_all$Etrans_kgH2Om2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
          site_output$mean_Etrans_kgH2Om2day = quantile(rowMeans(states_all$Etrans_kgH2Om2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_Etrans_kgH2Om2day = apply(t(apply(states_all$Etrans_kgH2Om2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          # Calculate the plant water use efficiency
          site_output$wue_plant_gCkgH2O = apply(states_all$gpp_gCm2day/states_all$Etrans_kgH2Om2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
          site_output$mean_wue_plant_gCkgH2O = quantile(rowSums(states_all$gpp_gCm2day)/rowSums(states_all$Etrans_kgH2Om2day), prob=num_quantiles, na.rm = na_flag)
          site_output$mean_annual_wue_plant_gCkgH2O = apply(t(apply(states_all$gpp_gCm2day,1, rollapply_mean_annual, step = steps_per_year) /
                                                              apply(states_all$Etrans_kgH2Om2day,1, rollapply_mean_annual, step = steps_per_year)),
                                                            2, quantile, prob=num_quantiles, na.rm = TRUE)
      }
      if (any(check_list == "Esoil_kgH2Om2day")) {
          # Soil evaporation
          site_output$Esoil_kgH2Om2day = apply(states_all$Esoil_kgH2Om2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
          site_output$mean_Esoil_kgH2Om2day = quantile(rowMeans(states_all$Esoil_kgH2Om2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_Esoil_kgH2Om2day = apply(t(apply(states_all$Esoil_kgH2Om2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      }
      if (any(check_list == "Ewetcanopy_kgH2Om2day")) {
          # Wet canopy evaporation
          site_output$Ewetcanopy_kgH2Om2day = apply(states_all$Ewetcanopy_kgH2Om2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
          site_output$mean_Ewetcanopy_kgH2Om2day = quantile(rowMeans(states_all$Ewetcanopy_kgH2Om2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_Ewetcanopy_kgH2Om2day = apply(t(apply(states_all$Ewetcanopy_kgH2Om2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      }
      # Check whether surface runoff exists?
      if (any(check_list == "runoff_kgH2Om2day")) {
          # Surface water drainage
          site_output$runoff_kgH2Om2day = apply(states_all$runoff_kgH2Om2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
          site_output$mean_runoff_kgH2Om2day = quantile(rowMeans(states_all$runoff_kgH2Om2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_runoff_kgH2Om2day = apply(t(apply(states_all$runoff_kgH2Om2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          # Accumulating total drainage (runoff and underflow)
          site_output$total_drainage_kgH2Om2day = states_all$runoff_kgH2Om2day
      }
      # Check whether underflow exists, i.e. drainage out of the bottom of the soil water column?
      if (any(check_list == "underflow_kgH2Om2day")) {
          # Drainage from bottom of soil column
          site_output$underflow_kgH2Om2day = apply(states_all$underflow_kgH2Om2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
          site_output$mean_underflow_kgH2Om2day = quantile(rowMeans(states_all$underflow_kgH2Om2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_underflow_kgH2Om2day = apply(t(apply(states_all$underflow_kgH2Om2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
          # Accumulating total drainage (runoff and underflow)
          site_output$total_drainage_kgH2Om2day = site_output$total_drainage_kgH2Om2day + states_all$underflow_kgH2Om2day
      }
      # Update the total drainage variable
      if (any(names(site_output) == "total_drainage_kgH2Om2day")) {
          # Total drainage from surface and soil bottom
          site_output$total_drainage_kgH2Om2day = apply(site_output$total_drainage_kgH2Om2day,2,quantile,prob=num_quantiles,na.rm = na_flag)
          site_output$mean_total_drainage_kgH2Om2day = quantile(rowMeans(site_output$total_drainage_kgH2Om2day, na.rm = na_flag), prob=num_quantiles)
          site_output$mean_annual_total_drainage_kgH2Om2day = apply(t(apply(site_output$total_drainage_kgH2Om2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      }
      # Update the mean Leaf Water Potential
      if (any(check_list == "LWP_MPa")){ # if (exists(x = "LWP_MPa", where = site_output)){
          site_output$LWP_MPa = apply(states_all$LWP_MPa,2,quantile,prob=num_quantiles,na.rm = na_flag)
          site_output$LWP_MPa = quantile(rowMeans(states_all$LWP_MPa, na.rm = na_flag), prob=num_quantiles)
      }
  } # water cycle variables

  # Snow related
  if (any(check_list == "snow_kgH2Om2")) {
      ## Snow on soil surface
      site_output$snow_kgH2Om2 = apply(states_all$snow_kgH2Om2,2,quantile,prob=num_quantiles,na.rm = na_flag)
      site_output$mean_snow_kgH2Om2 = quantile(rowMeans(states_all$snow_kgH2Om2, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_snow_kgH2Om2 = apply(t(apply(states_all$snow_kgH2Om2,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
  }
  
  ###
  # Aggregate ACM diagnositic information
  ###

  # Canopy process information if available
  if (any(check_list == "APAR_MJm2day")) {
      # Extract the absorbed photosynthetically active radiation by the canopy
      site_output$APAR_MJm2day = apply(states_all$APAR_MJm2day,2,quantile, prob=num_quantiles, na.rm = na_flag)
      site_output$mean_APAR_MJm2day = quantile(rowMeans(states_all$APAR_MJm2day, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_APAR_MJm2day = apply(t(apply(states_all$APAR_MJm2day,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Calculate change over time
      dCbio = states_all$APAR_MJm2day - states_all$APAR_MJm2day[,1] # difference in dom from initial
      site_output$dAPAR_MJm2day = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
  }
  if (any(check_list == "CiCa")) {
      # Extract the internal vs ambient CO2 ratio
      site_output$CiCa = apply(states_all$CiCa,2,quantile, prob=num_quantiles, na.rm = na_flag)
      site_output$mean_CiCa = quantile(rowMeans(states_all$CiCa, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_CiCa = apply(t(apply(states_all$CiCa,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Calculate change over time
      dCbio = states_all$CiCa - states_all$CiCa[,1] # difference in dom from initial
      site_output$dCiCa = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
  }
  if (any(check_list == "gs_demand_supply_ratio")) {
      # Extract the ratio of stomatal conductance relative to its maximum value,
      # this metric provides information on the demand vs supply constrains on stomatal conductance
      site_output$gs_demand_supply_ratio = apply(states_all$gs_demand_supply_ratio,2,quantile, prob=num_quantiles, na.rm = na_flag)
      site_output$mean_gs_demand_supply_ratio = quantile(rowMeans(states_all$gs_demand_supply_ratio, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_gs_demand_supply_ratio = apply(t(apply(states_all$gs_demand_supply_ratio,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Calculate change over time
      dCbio = states_all$gs_demand_supply_ratio - states_all$gs_demand_supply_ratio[,1] # difference in dom from initial
      site_output$dgs_demand_supply_ratio = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
  }
  if (any(check_list == "gs_mmolH2Om2s")) {
      # Extract the canopy stomatal conductance
      site_output$gs_mmolH2Om2s = apply(states_all$gs_mmolH2Om2s,2,quantile, prob=num_quantiles, na.rm = na_flag)
      site_output$mean_gs_mmolH2Om2s = quantile(rowMeans(states_all$gs_mmolH2Om2s, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_gs_mmolH2Om2s = apply(t(apply(states_all$gs_mmolH2Om2s,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Calculate change over time
      dCbio = states_all$gs_mmolH2Om2s - states_all$gs_mmolH2Om2s[,1] # difference in dom from initial
      site_output$dgs_mmolH2Om2s = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
  }
  if (any(check_list == "gb_mmolH2Om2s")) {
      # Extract the canopy boundary layer conductance
      site_output$gb_mmolH2Om2s = apply(states_all$gb_mmolH2Om2s,2,quantile, prob=num_quantiles, na.rm = na_flag)
      site_output$mean_gb_mmolH2Om2s = quantile(rowMeans(states_all$gb_mmolH2Om2s, na.rm = na_flag), prob=num_quantiles)
      site_output$mean_annual_gb_mmolH2Om2s = apply(t(apply(states_all$gb_mmolH2Om2s,1, rollapply_mean_annual, step = steps_per_year)), 2,quantile, prob=num_quantiles, na.rm = TRUE)
      # Calculate change over time
      dCbio = states_all$gb_mmolH2Om2s - states_all$gb_mmolH2Om2s[,1] # difference in dom from initial
      site_output$dgb_mmolH2Om2s = apply(dCbio,2,quantile,prob=num_quantiles,na.rm = na_flag)
  }

  ###
  # Aggregate model ensemble - observation uncertainty consistency
  # information
  ###

  # Any time series assimilated data overlaps?
  if (any(check_list == "gpp_assim_data_overlap_fraction")) {
      site_output$gpp_assim_data_overlap_fraction = states_all$gpp_assim_data_overlap_fraction
  }
  if (any(check_list == "lai_assim_data_overlap_fraction")) {
      site_output$lai_assim_data_overlap_fraction = states_all$lai_assim_data_overlap_fraction
  }
  if (any(check_list == "nee_assim_data_overlap_fraction")) {
      site_output$nee_assim_data_overlap_fraction = states_all$nee_assim_data_overlap_fraction
  }
  if (any(check_list == "wood_assim_data_overlap_fraction")) {
      site_output$wood_assim_data_overlap_fraction = states_all$wood_assim_data_overlap_fraction
  }
  if (any(check_list == "soil_assim_data_overlap_fraction")) {
      site_output$soil_assim_data_overlap_fraction = states_all$soil_assim_data_overlap_fraction
  }
  if (any(check_list == "et_assim_data_overlap_fraction")) {
      site_output$et_assim_data_overlap_fraction = states_all$et_assim_data_overlap_fraction
  }
  if (any(check_list == "nbe_assim_data_overlap_fraction")) {
      site_output$nbe_assim_data_overlap_fraction = states_all$nbe_assim_data_overlap_fraction
  }
  if (any(check_list == "fire_assim_data_overlap_fraction")) {
      site_output$fire_assim_data_overlap_fraction = states_all$fire_assim_data_overlap_fraction
  }

  # Store mean absolute parameter correlation information
  site_output$absolute_mean_parameter_correlation = states_all$absolute_mean_parameter_correlation
  # C-cycle flux correlation with parameters
  site_output$nee_parameter_correlation = states_all$nee_parameter_correlation
  site_output$gpp_parameter_correlation = states_all$gpp_parameter_correlation
  site_output$rauto_parameter_correlation = states_all$rauto_parameter_correlation
  site_output$rhet_parameter_correlation = states_all$rhet_parameter_correlation
  site_output$fire_parameter_correlation = states_all$fire_parameter_correlation
  # If Mean transit time for wood correlation exists, ensure we store it for the gridded run too
  if (any(check_list == "MTT_wood_years_parameter_correlation")) {
      site_output$MTT_wood_years_parameter_correlation = states_all$MTT_wood_years_parameter_correlation
  }
  # If Mean mean allocation to wood correlation exists, ensure we store it for the gridded run too
  if (any(check_list == "NPP_wood_gCm2day_parameter_correlation")) {
      site_output$NPP_wood_gCm2day_parameter_correlation = states_all$NPP_wood_gCm2day_parameter_correlation
  }
  # If the correlation between wood MTT and wood allocation have been determined
  if (any(check_list == "MTT_wood_years_to_NPP_wood_gCm2day_correlation")) {
      site_output$MTT_wood_years_to_NPP_wood_gCm2day_correlation = states_all$MTT_wood_years_to_NPP_wood_gCm2day_correlation
  }

  # save to pixel specific file for the moment... in "run_mcmc_results" these will be combined into a single grid
  save(site_output,file=outfile_stock_fluxes, compress = "gzip", compression_level = 6)

  # Return all clear message to ensure that the file name is added to the list reported back up
  dummy = 0 ; return(dummy)

} # end function post_process_for_grid
## Use byte compile
post_process_for_grid<-cmpfun(post_process_for_grid)


