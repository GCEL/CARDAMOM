
###
# Derive stocks and fluxes used in the calculation of gridded aggregates
# These are variables which for a site analysis would be easy to calculate
# from the ensembles but difficult if not determined here and now before aggregation
###

# This function was created by T. L Smallman (t.l.smallman@ed.ac.uk, UoE)

post_process_dalec<-function(states_all,parameters,drivers,PROJECT,n) {

  # Determine some useful information for the analysis below
  nos_years = PROJECT$nos_years
  steps_per_year = floor(dim(drivers$met)[1] / nos_years)

  # Check the list variables in states_all which we will be searching
  check_list = names(states_all)

  # If a combined ecosystem heterotrophic respiration flux does not
  # exist we shall calculate it
  if (any(check_list == "rhet_gCm2day") == FALSE) {
      if (any(check_list == "rhet_dom_gCm2day")) {
          states_all$rhet_gCm2day = states_all$rhet_dom_gCm2day
      } else {
          # Calculate the combined ecosystem heterotrophic respiration.
          # All models have a som pool, so start with that
          states_all$rhet_gCm2day = states_all$rhet_som_gCm2day
          # If the model has a litter pool (foliar + fine root) add this
          if (any(check_list == "rhet_litter_gCm2day")) {
              states_all$rhet_gCm2day = states_all$rhet_gCm2day + states_all$rhet_litter_gCm2day
          }
          # If the model has a wood litter pool add this
          if (any(check_list == "rhet_woodlitter_gCm2day")) {
              states_all$rhet_gCm2day = states_all$rhet_gCm2day + states_all$rhet_woodlitter_gCm2day
          }
      } # does rhet_dom_gCm2day exist?
  } # does rhet_gCm2day exist?
  # Combine autotrophic and heterotrophic respiration into ecosystem respiration
  states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  # Calculate the net ecosystem exchange of CO2
  states_all$nee_gCm2day = states_all$reco_gCm2day - states_all$gpp_gCm2day
  # Calculate net primary productivity
  states_all$npp_gCm2day = states_all$gpp_gCm2day - states_all$rauto_gCm2day
  # Begin calculation of net biome exchange and net biome productivity
  states_all$nbe_gCm2day = states_all$nee_gCm2day  # negative = sink
  states_all$nbp_gCm2day = -states_all$nee_gCm2day # positive = sink
  # If fire exists then update the NBE and NBP accordingly
  if (any(check_list == "fire_gCm2day")) {
      states_all$nbe_gCm2day = states_all$nbe_gCm2day + states_all$fire_gCm2day
      states_all$nbp_gCm2day = states_all$nbp_gCm2day - states_all$fire_gCm2day
  }
  # If a harvest flux exists update the NBP. NOTE: that this harvest flux
  # specifically accouts for C removed, there may be mortality due to harvest
  # but remains in system as residues.
  if (any(check_list == "harvest_gCm2day")) {
      states_all$nbp_gCm2day = states_all$nbp_gCm2day - states_all$harvest_gCm2day
  }
  # Now calculate the mean annual carbon use efficiency (NPP:GPP) as some models do now have a parameter for this
  # NOTE: rollapply inverts the dimensions from that wanted, hence t()
  states_all$mean_annual_cue = t( apply(states_all$npp_gCm2day,1, rollapply_mean_annual, step = steps_per_year)
                                / apply(states_all$gpp_gCm2day,1, rollapply_mean_annual, step = steps_per_year) )

  ###
  ## Post-hoc calculation of parameter correlations with key C-cycle variables

  # Construct and rearrange array of parameter suitable for correlation determination
  tmp = t(array(as.vector(parameters[1:PROJECT$model$nopars[n],,]),dim=c(PROJECT$model$nopars[n],prod(dim(parameters)[2:3]))))
  # Determine the correlation matrix between all parameters
  states_all$absolute_mean_parameter_correlation = cor(tmp)
  # Determine the mean of the absolute correlations from the matrix
  states_all$absolute_mean_parameter_correlation = mean(abs(states_all$absolute_mean_parameter_correlation[lower.tri(states_all$absolute_mean_parameter_correlation,diag=FALSE)]))

  # Determine correlations between parameter values and various state variables
  states_all$nee_parameter_correlation = cor(tmp,rowMeans(states_all$nee_gCm2day))
  states_all$gpp_parameter_correlation = cor(tmp,rowMeans(states_all$gpp_gCm2day))
  states_all$rauto_parameter_correlation = cor(tmp,rowMeans(states_all$rauto_gCm2day))
  states_all$rhet_parameter_correlation = cor(tmp,rowMeans(states_all$rhet_gCm2day))
  # Avoid error flag when no fire
  if (any(check_list == "fire_gCm2day")) {
      if (max(as.vector(states_all$fire_gCm2day)) > 0) {
          states_all$fire_parameter_correlation = cor(tmp,rowMeans(states_all$fire_gCm2day))
      } else {
          states_all$fire_parameter_correlation = array(0, dim = c(PROJECT$model$nopars[n],1))
      }
  }
  # Determine whether have have both mean transit time and allocation to wood
  if (any(check_list == "MTT_wood_years") & any(check_list == "alloc_wood_gCm2day")) {
      # As both exist determine their correlations with parameters...
      states_all$MTT_wood_years_parameter_correlation = cor(tmp,states_all$MTT_wood_years)
      states_all$NPP_wood_gCm2day_parameter_correlation = cor(tmp,rowMeans(states_all$alloc_wood_gCm2day))
      # ...and with each other
      states_all$MTT_wood_years_to_NPP_wood_gCm2day_correlation = cor(states_all$MTT_wood_years,rowMeans(states_all$alloc_wood_gCm2day))
  } else {
      # Both are not present, so we will determine whether we can generate one of the correlation estimates

      # If Mean transit time for wood is provided generate a correlation estimate
      if (any(check_list == "MTT_wood_years")) {
          states_all$MTT_wood_years_parameter_correlation = cor(tmp,states_all$MTT_wood_years)
      }
      # If Mean mean allocation to wood is provided generate a correlation estimate
      if (any(check_list == "alloc_wood_gCm2day")) {
          states_all$NPP_wood_gCm2day_parameter_correlation = cor(tmp,rowMeans(states_all$alloc_wood_gCm2day))
          }
  } # Both MTT wood and alloc_wood present?

  # Return back to user
  return(states_all)
          
} # end function post_process_dalec
## Use byte compile
post_process_dalec<-cmpfun(post_process_dalec)

###
# Quantify the proportion of ensemble members within the uncertainty bounds of the calibration datasets
###

# This function was created by T. L Smallman (t.l.smallman@ed.ac.uk, UoE)

assess_ensemble_fit_to_calibration_data<-function(states_all,drivers,PROJECT) {

  ###
  ## Comparison with assimilated observation - to what extent does the ensemble overlap?

  ## GPP (gC/m2/day)
  obs_id = 1 ; unc_id = obs_id+1
  if (length(which(drivers$obs[,obs_id] != -9999)) > 0) {
      # Loop through time to assess model overlap with observations
      nobs = 0 ; states_all$gpp_assim_data_overlap_fraction = 0
      to_do = which(drivers$obs[,obs_id] != -9999)
      for (a in seq(1, length(to_do))) {
           # Assign correct time step
           t = to_do[a]
           # Estimate the min / max values for the observations
           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
           # Create list object containing each observations distributions
           hist_list = list(o = c(obs_min,obs_max), m = states_all$gpp_gCm2day[,t])
           # Estimate average model ensemble within observated range
           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
           states_all$gpp_assim_data_overlap_fraction = states_all$gpp_assim_data_overlap_fraction + tmp2
           nobs = nobs + 1
      } # time loop
      # Average the overlap
      if (nobs > 0) {
          states_all$gpp_assim_data_overlap_fraction = states_all$gpp_assim_data_overlap_fraction / nobs
      } else {
          states_all$gpp_assim_data_overlap_fraction = 0
      }
  } # was the obs assimilated?

  ## LAI (m2/m2)
  obs_id = 3 ; unc_id = obs_id+1
  if (length(which(drivers$obs[,obs_id] != -9999)) > 0) {
      # Loop through time to assess model overlap with observations
      nobs = 0 ; states_all$lai_assim_data_overlap_fraction = 0
      to_do = which(drivers$obs[,obs_id] != -9999)
      for (a in seq(1, length(to_do))) {
           # Assign correct time step
           t = to_do[a]
           # Estimate the min / max values for the observations
           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
           # Create list object containing each observations distributions
           hist_list = list(o = c(obs_min,obs_max), m = states_all$lai_m2m2[,t])
           # Estimate average model ensemble within observated range
           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
           states_all$lai_assim_data_overlap_fraction = states_all$lai_assim_data_overlap_fraction + tmp2
           nobs = nobs + 1
           } # time loop
      # Average the overlap
      if (nobs > 0) {
          states_all$lai_assim_data_overlap_fraction = states_all$lai_assim_data_overlap_fraction / nobs
      } else {
          states_all$lai_assim_data_overlap_fraction = 0
      }
  } # was the obs assimilated?

  ## NEE (gC/m2/day)
  obs_id = 5 ; unc_id = obs_id+1
  if (length(which(drivers$obs[,obs_id] != -9999)) > 0) {
      # Loop through time to assess model overlap with observations
      nobs = 0 ; states_all$nee_assim_data_overlap_fraction = 0
      to_do = which(drivers$obs[,obs_id] != -9999)
      for (a in seq(1, length(to_do))) {
           # Assign correct time step
           t = to_do[a]
           # Estimate the min / max values for the observations
           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
           # Create list object containing each observations distributions
           hist_list = list(o = c(obs_min,obs_max), m = states_all$nee_gCm2day[,t])
           # Estimate average model ensemble within observated range
           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
           states_all$nee_assim_data_overlap_fraction = states_all$nee_assim_data_overlap_fraction + tmp2
           nobs = nobs + 1
      } # time loop
      # Average the overlap
      if (nobs > 0) {
          states_all$nee_assim_data_overlap_fraction = states_all$nee_assim_data_overlap_fraction / nobs
      } else {
          states_all$nee_assim_data_overlap_fraction = 0
      }
  } # was the obs assimilated?

  ## Reco (gC/m2/day)
  obs_id = 9 ; unc_id = obs_id+1
  if (length(which(drivers$obs[,obs_id] != -9999)) > 0) {
      # Loop through time to assess model overlap with observations
      nobs = 0 ; states_all$reco_assim_data_overlap_fraction = 0
      to_do = which(drivers$obs[,obs_id] != -9999)
      for (a in seq(1, length(to_do))) {
           # Assign correct time step
           t = to_do[a]
           # Estimate the min / max values for the observations
           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
           # Create list object containing each observations distributions
           hist_list = list(o = c(obs_min,obs_max), m = states_all$nee_gCm2day[,t])
           # Estimate average model ensemble within observated range
           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
           states_all$reco_assim_data_overlap_fraction = states_all$reco_assim_data_overlap_fraction + tmp2
           nobs = nobs + 1
      } # time loop
      # Average the overlap
      if (nobs > 0) {
          states_all$reco_assim_data_overlap_fraction = states_all$reco_assim_data_overlap_fraction / nobs
      } else {
          states_all$reco_assim_data_overlap_fraction = 0
      }
  } # was the obs assimilated?

  ## Wood (gC/m2)
  obs_id = 13 ; unc_id = obs_id+1
  # If there is a prior assign it to the first timestep of the observation timeseries
  if (drivers$parpriors[21] > 0) { drivers$obs[1,obs_id] = drivers$parpriors[21] ; drivers$obs[1,unc_id] = drivers$parpriorunc[21] }
  if (length(which(drivers$obs[,obs_id] != -9999)) > 0) {
      # Loop through time to assess model overlap with observations
      nobs = 0 ; states_all$wood_assim_data_overlap_fraction = 0
      to_do = which(drivers$obs[,obs_id] != -9999)
      for (a in seq(1, length(to_do))) {
           # Assign correct time step
           t = to_do[a]
           # Estimate the min / max values for the observations
           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
           # Create list object containing each observations distributions
           hist_list = list(o = c(obs_min,obs_max), m = states_all$wood_gCm2[,t])
           # Estimate average model ensemble within observated range
           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
           states_all$wood_assim_data_overlap_fraction = states_all$wood_assim_data_overlap_fraction + tmp2
           nobs = nobs + 1
      } # time loop
      # Average the overlap
      if (nobs > 0) {
          states_all$wood_assim_data_overlap_fraction = states_all$wood_assim_data_overlap_fraction / nobs
      } else {
          states_all$wood_assim_data_overlap_fraction = 0
      }
  } # was the obs assimilated?

  ## Soil (gC/m2)
  obs_id = 19 ; unc_id = obs_id+1
  # If there is a prior assign it to the first timestep of the observation timeseries
  if (drivers$parpriors[23] > 0) { drivers$obs[1,obs_id] = drivers$parpriors[23] ; drivers$obs[1,unc_id] = drivers$parpriorunc[23] }
  if (length(which(drivers$obs[,obs_id] != -9999)) > 0) {
      # Loop through time to assess model overlap with observations
      nobs = 0 ; states_all$soil_assim_data_overlap_fraction = 0
      to_do = which(drivers$obs[,obs_id] != -9999)
      for (a in seq(1, length(to_do))) {
           # Assign correct time step
           t = to_do[a]
           # Estimate the min / max values for the observations
           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
           # Create list object containing each observations distributions
           hist_list = list(o = c(obs_min,obs_max), m = states_all$som_gCm2[,t])
           # Estimate average model ensemble within observated range
           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
           states_all$soil_assim_data_overlap_fraction = states_all$soil_assim_data_overlap_fraction + tmp2
           nobs = nobs + 1
      } # time loop
      # Average the overlap
      if (nobs > 0) {
          states_all$soil_assim_data_overlap_fraction = states_all$soil_assim_data_overlap_fraction / nobs
      } else {
          states_all$soil_assim_data_overlap_fraction = 0
      }
  } # was the obs assimilated?

  ## ET (kgH2O/m2/day)
  obs_id = 31 ; unc_id = obs_id+1
  if (length(which(drivers$obs[,obs_id] != -9999)) > 0) {
      # Loop through time to assess model overlap with observations
      nobs = 0 ; states_all$et_assim_data_overlap_fraction = 0
      to_do = which(drivers$obs[,obs_id] != -9999)
      for (a in seq(1, length(to_do))) {
           # Assign correct time step
           t = to_do[a]
           # Estimate the min / max values for the observations
           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
           # Create list object containing each observations distributions
           hist_list = list(o = c(obs_min,obs_max), m = states_all$ET_kgH2Om2day[,t])
           # Estimate average model ensemble within observated range
           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
           states_all$et_assim_data_overlap_fraction = states_all$et_assim_data_overlap_fraction + tmp2
           nobs = nobs + 1
      } # time loop
      # Average the overlap
      if (nobs > 0) {
          states_all$et_assim_data_overlap_fraction = states_all$et_assim_data_overlap_fraction / nobs
      } else {
          states_all$et_assim_data_overlap_fraction = 0
      }
  } # was the obs assimilated?

  ## NBE (gC/m2/day)
  obs_id = 35 ; unc_id = obs_id+1 
  if (length(which(drivers$obs[,obs_id] != -9999)) > 0) {
      # Loop through time to assess model overlap with observations
      nobs = 0 ; states_all$nbe_assim_data_overlap_fraction = 0
      to_do = which(drivers$obs[,obs_id] != -9999)
      for (a in seq(1, length(to_do))) {
           # Assign correct time step
           t = to_do[a]
           # Estimate the min / max values for the observations
           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
           # Create list object containing each observations distributions
           hist_list = list(o = c(obs_min,obs_max), m = states_all$nbe_gCm2day[,t])
           # Estimate average model ensemble within observated range
           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
           states_all$nbe_assim_data_overlap_fraction = states_all$nbe_assim_data_overlap_fraction + tmp2
           nobs = nobs + 1
      } # time loop
      # Average the overlap
      if (nobs > 0) {
          states_all$nbe_assim_data_overlap_fraction = states_all$nbe_assim_data_overlap_fraction / nobs
      } else {
          states_all$nbe_assim_data_overlap_fraction = 0
      }
  } # was the obs assimilated?

  ## Fire (gC/m2/day)
  obs_id = 7 ; unc_id = obs_id+1
  if (length(which(drivers$obs[,obs_id] != -9999)) > 0) {
      # Loop through time to assess model overlap with observations
      nobs = 0 ; states_all$fire_assim_data_overlap_fraction = 0
      to_do = which(drivers$obs[,obs_id] != -9999)
      for (a in seq(1, length(to_do))) {
           # Assign correct time step
           t = to_do[a]
           # Estimate the min / max values for the observations
           obs_max = drivers$obs[t,obs_id] + drivers$obs[t,unc_id]
           obs_min = drivers$obs[t,obs_id] - drivers$obs[t,unc_id]
           # Create list object containing each observations distributions
           hist_list = list(o = c(obs_min,obs_max), m = states_all$fire_gCm2day[,t])
           # Estimate average model ensemble within observated range
           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
           states_all$fire_assim_data_overlap_fraction = states_all$fire_assim_data_overlap_fraction + tmp2
           nobs = nobs + 1
      } # time loop
      # Average the overlap
      if (nobs > 0) {
          states_all$fire_assim_data_overlap_fraction = states_all$fire_assim_data_overlap_fraction / nobs
      } else {
          states_all$fire_assim_data_overlap_fraction = 0
      }
  } # was the obs assimilated?

  # Return back to user
  return(states_all)

} # end function assess_ensemble_fit_to_calibration_data
## Use byte compile
assess_ensemble_fit_to_calibration_data<-cmpfun(assess_ensemble_fit_to_calibration_data)

