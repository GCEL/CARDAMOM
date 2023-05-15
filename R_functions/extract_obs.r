
###
## Function to extract obs needed for CARDAMOM
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

extract_obs<-function(grid_long_loc,grid_lat_loc,latlon_wanted,lai_all,Csom_all,forest_all
                     ,Cwood_initial_all,Cwood_stock_all,Cwood_potential_all
                     ,sand_clay_all,crop_man_all,burnt_all,soilwater_all,nbe_all
                     ,lca_all,lifespan_all,leaf_fall_period_all,gpp_all,Cwood_inc_all,Cwood_mortality_all
					 ,fire_all
                     ,ctessel_pft,site_name,start_year,end_year
                     ,timestep_days,spatial_type,resolution,grid_type,modelname) {

    # Create useful timing information for multiple functions
    # Years to be simulated
    years_to_load = as.numeric(start_year):as.numeric(end_year)
    # Determine how many days are in each year
    doy_obs = 0
    for (i in seq(1, length(years_to_load))) {
         nos_days = nos_days_in_year(years_to_load[i])
         # count up days needed
         doy_obs = append(doy_obs,1:nos_days)
         # Days per time step if not already provided
         if (length(timestep_days) == 1) {timestep_days = rep(timestep_days, nos_days, by = timestep_days)}
    }
    doy_obs = doy_obs[-1]

    # Create useful timing information for multiple functions
    if (length(timestep_days) == 1) {
        analysis_years = seq(as.numeric(start_year),as.numeric(end_year))
        nos_days = 0
        for (t in seq(1,length(analysis_years))) {nos_days = nos_days + nos_days_in_year(analysis_years[t])}
        timestep_days = rep(timestep_days, nos_days, by = timestep_days)
    }

    ###
    ## Get some NBE information (gC/m2/day); negative is sink
    ###

    if (nbe_source == "GEOSCHEM" | nbe_source == "Global_Combined") {

      # Extract NBE and uncertainty information
      # NOTE: assume default uncertainty (+/- scale)
      output = extract_nbe(grid_long_loc,grid_lat_loc,timestep_days,
                           spatial_type,resolution,grid_type,latlon_wanted,
                           nbe_all,years_to_load,doy_obs)
      nbe = output$nbe ; nbe_unc = output$nbe_unc

    } else if (nbe_source == "site_specific") {

      # read from .csv or netcdf
      infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
      nbe = read_site_specific_obs("NBE_gCm2day",infile) ; nbe_unc = read_site_specific_obs("NBE_unc_gCm2day",infile)
      if (max(nbe_unc) == -9999) {
          nbe_unc = rep(-9999,times = length(nbe))
          # apply default uncertainty consistent with Eddy covariance estimates
          nbe_unc[which(nbe != -9999)] = 1.0
      }

    } else {

      nbe = -9999
      nbe_unc = -9999

    }
    # Add model structural uncertainty to the uncertainty estimate if present
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    nbe_unc[nbe_unc >= 0] = sqrt(nbe_unc[nbe_unc >= 0]**2 + 1.0**2)

    ###
    ## Get some LAI information (m2/m2)
    ###

    if (lai_source == "MODIS" | lai_source == "COPERNICUS") {

        # Extract lai and uncertainty information
        # NOTE: assume default uncertainty (+/- scale)
        output = extract_lai_timeseries(grid_long_loc,grid_lat_loc,timestep_days,
                                        spatial_type,resolution,grid_type,
                                        latlon_wanted,lai_all,years_to_load,doy_obs)
        lai = output$lai ; lai_unc = output$lai_unc

    } else if (lai_source == "site_specific") {

        # read from .csv or netcdf
        infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        lai = read_site_specific_obs("LAI_m2m2",infile) ; lai_unc = read_site_specific_obs("LAI_unc_m2m2",infile)
        if (max(lai_unc) == -9999) {
            lai_unc = rep(-9999,times = length(lai))
            # apply default uncertainty
            lai_unc[which(lai != -9999)] = 0.5
        }

    } else {

        lai = -9999
        lai_unc = -9999

    }
    # Assume minimum uncertainty to reflect model structural uncertainty
    # Estimates from comparison of LAI uncertainties trials at 0.5 and 0.25,
    # resultant CI in both instances is range of ~0.50. Therefore CI of +/- 0.25
    #lai_unc[lai_unc >= 0] = sqrt(lai_unc[lai_unc >= 0]**2 + 0.25**2)
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    if (length(which(lai_unc >= 0)) > 0) {
        lai_unc[lai_unc >= 0] = pmax(0.25,sqrt(lai_unc[lai_unc >= 0]**2 + (0.1*mean(lai[lai >= 0]))**2))
    }

    ###
    ## Get some Cfoliage information (stock; gC/m2)
    ###

    if (Cfol_stock_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Cfol_stock=read_site_specific_obs("Cfol_stock_gCm2",infile)
        Cfol_stock_unc=read_site_specific_obs("Cfol_stock_unc_gCm2",infile)
        if (length(Cfol_stock_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cfol_stock_unc = rep(-9999,times = length(Cfol_stock))
            # See Smallman et al., (2017) for uncertainty estimate
            Cfol_stock_unc[which(Cfol_stock > 0)] = 0.38 * Cfol_stock[which(Cfol_stock > 0)]
        }
    } else {
        # assume no data available
        Cfol_stock = -9999 ; Cfol_stock_unc = -9999
    }

    ###
    ## Get some initial Csom (gC/m2) information
    ###

    if (Csom_source == "HWSD" | Csom_source == "SoilGrids"  | Csom_source == "SoilGrids_v2" | Csom_source == "NCSCD" | Csom_source == "NCSCD3m") {
        Csom_info = extract_Csom_prior(grid_long_loc,grid_lat_loc,spatial_type,
                                       resolution,grid_type,latlon_wanted,Csom_all)
        Csom_initial = Csom_info$Csom_initial ; Csom_initial_unc = Csom_info$Csom_initial_unc
    } else if (Csom_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        Csom_initial = read_site_specific_obs("Csom_initial_gCm2",infile)
        Csom_initial_unc = read_site_specific_obs("Csom_initial_unc_gCm2",infile)
        if (Csom_initial_unc == -9999 & Csom_initial > 0) {
            # on the other hand if not then we have no uncertainty info, so use default
            Csom_initial_unc = 0.24 * Csom_initial
        }
    } else {
        # assume no data available
        Csom_initial = -9999 ; Csom_initial_unc = -9999
    }
    # Now assuming we have actual information we need to add the model structural uncertainty.
    # A structural of uncertainty has been estimates at ~ 1000 gC/m2 based on Smallman et al., (2017)
    #Csom_initial_unc[Csom_initial_unc >= 0] = sqrt(Csom_initial_unc[Csom_initial_unc >= 0]**2 + 1000**2)
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    #Csom_initial_unc[Csom_initial_unc > 0] = pmax(1000,sqrt(Csom_initial_unc[Csom_initial_unc > 0]**2 + (0.1*mean(Csom_initial[Csom_initial > 0]))**2))
    # Assume structural uncertainty does not apply to initial conditions
    Csom_initial_unc[Csom_initial_unc > 0] = pmax(1000,Csom_initial_unc[Csom_initial_unc > 0])

    ###
    ## Get some sand / clay information (%)
    ###

    if (sand_clay_source == "HWSD" | sand_clay_source == "SoilGrids" | sand_clay_source == "SoilGrids_v2") {
        sand_clay=extract_sand_clay(grid_long_loc,grid_lat_loc,spatial_type,
                                    resolution,grid_type,latlon_wanted,sand_clay_all)
        top_sand = sand_clay$top_sand ; bot_sand = sand_clay$bot_sand
        top_clay = sand_clay$top_clay ; bot_clay = sand_clay$bot_clay
    } else if (sand_clay_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        top_sand = read_site_specific_obs("top_sand_initial_percent",infile)
        bot_sand = read_site_specific_obs("bot_sand_initial_percent",infile)
        top_clay = read_site_specific_obs("top_clay_initial_percent",infile)
        bot_clay = read_site_specific_obs("bot_clay_initial_percent",infile)
    } else {
        # assume no data available
        top_sand = 40 ; bot_sand = 40
        top_clay = 15 ; bot_clay = 15
    }

    ###
    ## Get some crop management information (day)
    ###

    if (ctessel_pft == 1 & crop_management_source == "sacks_crop_calendar") {
        # could add other variables such as SOM (gC.m-2)
        crop_dates=extract_sacks_crop_info(spatial_type,resolution,grid_type,latlon_wanted,crop_man_all)
        plant = crop_dates$plant ; plant_range = crop_dates$plant_range
        harvest = crop_dates$harvest ; harvest_range = crop_dates$harvest_range
    } else if (crop_management_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        plant = read_site_specific_obs("plant_initial",infile)
        plant_range = read_site_specific_obs("plant_range_initial",infile)
        harvest = read_site_specific_obs("harvest_initial",infile)
        harvest_range = read_site_specific_obs("harvest_range_initial",infile)
    } else {
        # assume no data available
        plant = 304 ; plant_range = 15 # days
        harvest = 208 ; harvest_range = 15 # days
    }

    ###
    ## Get some Wood increment information (gC/m2/day; time series)
    ###

    if (Cwood_inc_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Cwood_inc=read_site_specific_obs("Cwood_increment_gCm2day",infile)
        Cwood_inc_unc=read_site_specific_obs("Cwood_increment_uncertainty_gCm2day",infile)
        Cwood_inc_lag=read_site_specific_obs("Cwood_increment_lag_step",infile) # in model time steps
        # Has uncertainty information been provided?
        if (length(Cwood_inc_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cwood_inc_unc = rep(-9999,times = length(Cwood_inc))
            Cwood_inc_unc[which(Cwood_inc > 0)] = 0.25 * Cwood_inc[which(Cwood_inc > 0)]
        }
        # Has lag information been provided
        if (length(Cwood_inc_lag) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cwood_inc_lag = rep(-9999,times = length(Cwood_inc))
            Cwood_inc_lag[which(Cwood_inc > 0)] = 1 # assume applies to current time step only
        }
    } else if (Cwood_inc_source == "Rainfor") {
        # If there are any values in the analysis window
        if (max(Cwood_inc_all$place_obs_in_step) > 0) {
            # Extract the current location
            output = extract_wood_productivity(grid_long_loc,grid_lat_loc,timestep_days,
                                               spatial_type,resolution,grid_type,
                                               latlon_wanted,Cwood_inc_all)
            Cwood_inc = output$Cwood_inc ; Cwood_inc_unc = output$Cwood_inc_unc ; Cwood_inc_lag = output$Cwood_inc_lag
            # Tidy up
            rm(output)
        } else {
            # assume no data available
            Cwood_inc = -9999 ; Cwood_inc_unc = -9999 ; Cwood_inc_lag = -9999
        }
    } else {
        # assume no data available
        Cwood_inc = -9999 ; Cwood_inc_unc = -9999 ; Cwood_inc_lag = -9999
    }
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    Cwood_inc_unc[Cwood_inc_unc >= 0] = pmax(0.1,sqrt(Cwood_inc_unc[Cwood_inc_unc >= 0]**2 + (0.1*mean(Cwood_inc_unc[Cwood_inc_unc >= 0]))**2))

    ###
    ## Get some Wood natural mortality information (gC/m2/day; time series)
    ###

    if (Cwood_mortality_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Cwood_mortality=read_site_specific_obs("Cwood_mortality_gCm2day",infile)
        Cwood_mortality_unc=read_site_specific_obs("Cwood_mortality_uncertainty_gCm2day",infile)
        Cwood_mortality_lag=read_site_specific_obs("Cwood_mortality_lag_step",infile) # in model time steps
        # Has uncertainty information been provided?
        if (length(Cwood_mortality_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cwood_mortality_unc = rep(-9999,times = length(Cwood_mortality))
            Cwood_mortality_unc[which(Cwood_mortality > 0)] = 0.25 * Cwood_mortality[which(Cwood_mortality > 0)]
        }
        # Has lag information been provided
        if (length(Cwood_mortality_lag) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cwood_mortality_lag = rep(-9999,times = length(Cwood_mortality))
            Cwood_mortality_lag[which(Cwood_mortality > 0)] = 1 # assume applies to current time step only
        }
    } else if (Cwood_mortality_source == "Rainfor") {
        # If there are any values in the analysis window
        if (max(Cwood_mortality_all$place_obs_in_step) > 0) {
            # Extract the current location
            output = extract_wood_mortality(grid_long_loc,grid_lat_loc,timestep_days,
                                            spatial_type,resolution,grid_type,
                                            latlon_wanted,Cwood_mortality_all)
            Cwood_mortality = output$Cwood_mortality
            Cwood_mortality_unc = output$Cwood_mortality_unc
            Cwood_mortality_lag = output$Cwood_mortality_lag
            # Tidy up
            rm(output)
        } else {
            # assume no data available
            Cwood_mortality = -9999 ; Cwood_mortality_unc = -9999 ; Cwood_mortality_lag = -9999
        }
    } else {
        # assume no data available
        Cwood_mortality = -9999 ; Cwood_mortality_unc = -9999 ; Cwood_mortality_lag = -9999
    }
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    Cwood_mortality_unc[Cwood_mortality_unc >= 0] = pmax(0.1,sqrt(Cwood_mortality_unc[Cwood_mortality_unc >= 0]**2 + (0.1*mean(Cwood_mortality_unc[Cwood_mortality_unc >= 0]))**2))

    ###
    ## Get some GPP information (time series; gC/m2/day)
    ###

    if (GPP_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        if (modelname == "ACM") {infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")}
        #if (modelname == "ACM") {infile=paste(path_to_site_obs,site_name,"_timeseries_obs_iWUE_trunk_nowater.csv",sep="")}
        #if (modelname == "ACM") {infile=paste(path_to_site_obs,site_name,"_timeseries_obs_iWUE_trunk_nowater_copy.csv",sep="")}
#        if (modelname == "ACM") {infile=paste(path_to_site_obs,site_name,"_timeseries_obs_iWUE_trunk_nowater_avgNlessthan4.csv",sep="")}
        if (modelname == "ACM") {infile=paste(path_to_site_obs,site_name,"_timeseries_obs_iWUE_trunk_nowater_avgNlessthan4_subsample.csv.csv",sep="")}
        GPP = read_site_specific_obs("GPP_gCm2day",infile)
        GPP_unc = read_site_specific_obs("GPP_unc_gCm2day",infile)
        if (length(GPP_unc) == 1) {
            GPP_unc = rep(-9999,times = length(GPP))
            # Composed of NEE 0.58 gC/m2/day (Hill et al., 2012) plus mass balance mismatch of
            # 0.16 gC/m2/day, therefore 0.74 gC/m2/day
            GPP_unc[which(GPP > 0)] = 0.74
            if (modelname == "ACM") {GPP_unc = rep(mean(GPP)*0.40,times=length(GPP))}
        }
    } else if (GPP_source == "Global_Combined") {

        # Extract GPP and uncertainty information
        # NOTE: assume default uncertainty (+/- scale)
        output = extract_gpp(grid_long_loc,grid_lat_loc,timestep_days,spatial_type,
                             resolution,grid_type,latlon_wanted,gpp_all,years_to_load,doy_obs)
        GPP = output$GPP ; GPP_unc = output$GPP_unc

    } else {
        # assume no data available
        GPP = -9999 ; GPP_unc = -9999
    }
    # Combine with an estimate of model structural error.
    # A mean rmse of ~2 gC/m2/day was estimated when evaluating ACM-GPP-ET (ACM2)
    # GPP against fluxnet observations (Smallman & Williams 2019)
    # NOTE: this is a gross overestimate as ACM-GPP-ET was not calibrated against the data,
    # i.e. we don't know that it couldn't fit it
    #GPP_unc[GPP_unc >= 0] = sqrt(GPP_unc[GPP_unc >= 0]**2 + 2**2)
    # Assumed uncertainty structure as agreed with Anthony Bloom,
    # NOTE minimum bound also applied
    GPP_unc[GPP_unc >= 0] = pmax(1.0,sqrt(GPP_unc[GPP_unc >= 0]**2 + (0.1*mean(GPP[GPP >= 0]))**2))

    ###
    ## Get some fire C emission information (time series; gC/m2/day)
    ###

    if (fire_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Fire = read_site_specific_obs("fire_gCm2day",infile)
        Fire_unc = read_site_specific_obs("fire_unc_gCm2day",infile)
        if (length(Fire_unc) == 1) {
            Fire_unc = rep(-9999,times = length(Fire))
            # Ill defined assumption
            Fire_unc[which(Fire > 0)] = 0.1
        }
    } else if (fire_source == "Global_Combined") {

        # Extract Fire and uncertainty information
        # NOTE: assume default uncertainty (+/- scale)
        output = extract_fire(grid_long_loc,grid_lat_loc,timestep_days,spatial_type,
                              resolution,grid_type,latlon_wanted,fire_all,years_to_load,doy_obs)
        Fire = output$Fire ; Fire_unc = output$Fire_unc

    } else {
        # assume no data available
        Fire = -9999 ; Fire_unc = -9999
    }
    # Combine with an estimate of model structural error.
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates, based on median difference between GFED and GFAS
    #Fire_unc[Fire_unc >= 0] = pmax(0.1,sqrt(Fire_unc[Fire_unc >= 0]**2 + (0.1*mean(Fire[Fire >= 0]))**2))
    #Fire_unc[Fire_unc >= 0] = pmax(0.1*mean(Fire[Fire >= 0]),pmax(0.01,Fire_unc[Fire_unc >= 0]))
    Fire_unc[Fire_unc >= 0] = pmax(0.1*mean(Fire[Fire >= 0]),pmax(0.01,pmin(Fire[Fire >= 0],Fire_unc[Fire_unc >= 0])))

    # Fire emissions are dominated by zero (or effectively zero) values. This bias' the APMCMC (or MHMCMC)
    # away from fitting the more important emission values in favour of the more common zero (or near zero values.)
    # To address this we remove fire emission observations less than 0.01 gC/m2/day
    if (use_parallel == FALSE) {print("Fire observations < 0.01 gC/m2/day have been removed to prevent bias to MDF calibration")}
    Fire_unc[Fire < 0.01] = -9999 ; Fire[Fire < 0.01] = -9999

    ###
    ## Get some Evapotranspiration information (time series; kgH2O/m2/day)
    ###

    if (Evap_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        if (modelname == "ACM") {infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")}
        #if (modelname == "ACM") {infile=paste(path_to_site_obs,site_name,"_timeseries_obs_iWUE_trunk_nowater.csv",sep="")}
        #if (modelname == "ACM") {infile=paste(path_to_site_obs,site_name,"_timeseries_obs_iWUE_trunk_nowater_copy.csv",sep="")}
#        if (modelname == "ACM") {infile=paste(path_to_site_obs,site_name,"_timeseries_obs_iWUE_trunk_nowater_avgNlessthan4.csv",sep="")}
        if (modelname == "ACM") {infile=paste(path_to_site_obs,site_name,"_timeseries_obs_iWUE_trunk_nowater_avgNlessthan4_subsample.csv.csv",sep="")}
        Evap = read_site_specific_obs("Evap_kgH2Om2day",infile)
        Evap_unc = read_site_specific_obs("Evap_unc_kgH2Om2day",infile)
        if (length(Evap_unc) == 1) {
            Evap_unc = rep(-9999,times = length(Evap))
            Evap_unc[which(Evap > -9999)] = 0.77 # Assuming Hollinger & Richardson (2005) Tree Physiology, 25, 873-885
        }
        if (modelname == "ACM") {
            # borrow woody increment for soil evaporation in ACM_ET recalibration
            Cwood_inc = read_site_specific_obs("soilevap_kgH2Om2day",infile)
            # borrow Cfol_stock for wet canopy evaporation in ACM_ET recalibration
            Cfol_stock = read_site_specific_obs("wetevap_kgH2Om2day",infile)
            # actually lets make uncertainty half mean of total ET
            Evap_unc = rep(abs(mean(Evap))*0.40, length.out = length(Evap))
            Cwood_inc_unc = rep(abs(mean(Cwood_inc))*0.40, length.out = length(Evap))
            Cfol_stock_unc = rep(abs(mean(Cfol_stock))*0.40, length.out = length(Evap))
        }
    } else {
        # assume no data available
        Evap = -9999 ; Evap_unc = -9999
    }
    # Combine with an estimate of model structural error.
    # A mean rmse of ~1 kgH2O/m2/day was estimated when evaluating ACM-GPP-ET (ACM2)
    # ET against fluxnet observations (Smallman & Williams 2019)
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    Evap_unc[Evap_unc >= 0] = pmax(0.1,sqrt(Evap_unc[Evap_unc >= 0]**2 + (0.1*mean(Evap[Evap >= 0]))**2))

    ###
    ## Get some Reco information (time series; gC/m2/day)
    ###

    if (Reco_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Reco=read_site_specific_obs("Reco_gCm2day",infile)
        Reco_unc=read_site_specific_obs("Reco_unc_gCm2day",infile)
        if (length(Reco_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Reco_unc = rep(-9999,times = length(Reco))
            # Composed of NEE 0.58 gC/m2/day (Hill et al., 2012) plus mass balance mismatch of
            # 0.16 gC/m2/day, therefore 0.74 gC/m2/day
            Reco_unc[which(Reco > 0)] = 0.74
        }
    } else {
        # assume no data available
        Reco = -9999 ; Reco_unc = -9999
    }
    # apply lower bound in all cases to the uncertainty
    Reco_unc[Reco_unc >= 0] = sqrt(Reco_unc[Reco_unc >= 0]**2 + 0.50**2)

    ###
    ## Get some NEE information (time series; gC/m2/day)
    ###

    if (NEE_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        NEE=read_site_specific_obs("NEE_gCm2day",infile)
        NEE_unc=read_site_specific_obs("NEE_unc_gCm2day",infile)
        if (length(NEE_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            NEE_unc = rep(-9999,times = length(NEE))
            NEE_unc[which(NEE != -9999)] = 0.58 # Hill et al., (2012)
        }
    } else {
        # assume no data available
        NEE = -9999 ; NEE_unc = -9999
    }
    # The largest site level mean rmse achieved across the COMPLEX (Famiglietti et al., 2021) sites was 0.99 gC/m2/day.
    # This was achieved in the reduced uncertainty analysis providing and indication of the model structural error
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    NEE_unc[NEE_unc >= 0] = sqrt(NEE_unc[NEE_unc >= 0]**2 + 1**2)

    ###
    ## Get some Cfoliage information (initial conditions)
    ###

    if (Cfol_initial_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        Cfol_initial=read_site_specific_obs("Cfol_initial_gCm2",infile)
        Cfol_initial_unc=read_site_specific_obs("Cfol_initial_unc_gCm2",infile)
        if (Cfol_initial_unc == -9999 & Cfol_initial > 0) {
          # on the other hand if not then we have no uncertainty info, so use default
          Cfol_initial_unc = 0.25 * Cfol_initial
        }
    } else {
        # assume no data available
        Cfol_initial = -9999 ; Cfol_initial_unc = -9999
    }

    ###
    ## Get some Cwood information (initial conditions)
    ###

    if (Cwood_initial_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        Cwood_initial=read_site_specific_obs("Cwood_initial_gCm2",infile)
        Cwood_initial_unc=read_site_specific_obs("Cwood_initial_unc_gCm2",infile)
        if (Cwood_initial_unc == -9999 & Cwood_initial > 0) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cwood_initial_unc = 0.25 * Cwood_initial
        }

    } else if (Cwood_initial_source == "Avitabile" | Cwood_initial_source == "mpi_biomass" |
               Cwood_initial_source == "UoL_stable_forest" |Cwood_initial_source == "Rainfor" |
               Cwood_initial_source == "UoL_stable_savannah") {

        # All maps converted into common format, therefore a common extraction subroutine can be used
        output = extract_Cwood_initial(grid_long_loc,grid_lat_loc,spatial_type,
                                       resolution,grid_type,latlon_wanted,Cwood_initial_all)
        Cwood_initial = output$Cwood_stock ; Cwood_initial_unc = output$Cwood_stock_unc
    } else {
        # assume no data available
        Cwood_initial=-9999 ; Cwood_initial_unc=-9999
    }

    ###
    ## Get some Croots information (initial conditions)
    ###

    if (Croots_initial_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        Croots_initial=read_site_specific_obs("Croots_initial_gCm2",infile)
        Croots_initial_unc=read_site_specific_obs("Croots_initial_unc_gCm2",infile)
        if (Croots_initial_unc == -9999 & Croots_initial > 0) {
          # on the other hand if not then we have no uncertainty info, so use default
          Croots_initial_unc = 0.44 * Croots_initial
        }
    } else {
        # assume no data available
        Croots_initial = -9999 ; Croots_initial_unc = -9999
    }

    ###
    ## Get some Clitter information (initial conditions)
    ###

    if (Clit_initial_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        Clit_initial=read_site_specific_obs("Clit_initial_gCm2",infile)
        Clit_initial_unc=read_site_specific_obs("Clit_initial_unc_gCm2",infile)
        if (Clit_initial_unc == -9999 & Clit_initial > 0) {
            # on the other hand if not then we have no uncertainty info, so use default
            Clit_initial_unc = 0.25 * Clit_initial
        }
    } else {
        # assume no data available
        Clit_initial = -9999 ; Clit_initial_unc = -9999
    }

    ###
    ## Get some Cwood information (stock)
    ###

    if (Cwood_stock_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Cwood_stock = read_site_specific_obs("Cwood_stock_gCm2",infile)
        Cwood_stock_unc = read_site_specific_obs("Cwood_stock_unc_gCm2",infile)
        if (length(Cwood_stock_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cwood_stock_unc = rep(-9999,times = length(Cwood_stock))
            Cwood_stock_unc[which(Cwood_stock != -9999)] = abs(0.25 * Cwood_stock[which(Cwood_stock != -9999)])
        }
    } else if (Cwood_stock_source == "mpi_biomass" |
               Cwood_stock_source == "INPE_Avitabile" | Cwood_stock_source == "Avitabile" |
               Cwood_stock_source == "Rainfor" | Cwood_stock_source == "Rainfor_annual" |
               Cwood_stock_source == "McNicol" | Cwood_stock_source == "Biomass_maps_Africa_UoL" |
               Cwood_stock_source == "ESA_CCI_Biomass") {

        # All maps converted into common format, therefore a common extraction subroutine can be used
        if (max(Cwood_stock_all$place_obs_in_step) > 0) {
            output = extract_Cwood_stocks(grid_long_loc,grid_lat_loc,timestep_days,
                                          spatial_type,resolution,grid_type,latlon_wanted,
                                          Cwood_stock_all)
            Cwood_stock = output$Cwood_stock ; Cwood_stock_unc = output$Cwood_stock_unc
#            tmp = which(Cwood_stock > 0) # first AGB only
#            if (length(tmp) > 1) {
#                Cwood_stock[tmp[-1]] = -9999 ; Cwood_stock_unc[tmp[-1]] = -9999
#            }
        } else {
            Cwood_stock = rep(-9999, length(timestep_days))
            Cwood_stock_unc = rep(-9999, length(timestep_days))
        }

        # INPE Map is a merge of two separate ones, therefore it is a bad idea to have two time steps,
        # one from each map as these data are inconsistent.
        # Therefore, check whether we have two data points and take an average + average location
        if (Cwood_stock_source == "INPE_Avitabile") {
            tmp = which(Cwood_stock > 0)
            if (length(tmp) == 2) {
                tmp2 = Cwood_stock[tmp] ; tmp3 = Cwood_stock_unc[tmp]
                Cwood_stock[tmp] = -9999 ; Cwood_stock_unc[tmp] = -9999
                Cwood_stock[floor(mean(tmp))] = mean(tmp2) ; Cwood_stock_unc[floor(mean(tmp))] = mean(tmp3)
                rm(tmp,tmp2,tmp3)
            }
        } # INPE_Avitabile map adjustment
    } else {
        # assume no data available
        Cwood_stock = -9999 ; Cwood_stock_unc = -9999
    }
    # apply lower bound in all cases to the uncertainty
    #Cwood_stock_unc[Cwood_stock_unc >= 0] = sqrt(Cwood_stock_unc[Cwood_stock_unc >= 0]**2 + 100**2)
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    Cwood_stock_unc[Cwood_stock_unc >= 0] = pmax(100,sqrt(Cwood_stock_unc[Cwood_stock_unc >= 0]**2 + (0.1*mean(Cwood_stock[Cwood_stock >= 0]))**2))

    ###
    ## Get some Cagb information (stock)
    ###

    if (Cagb_stock_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Cagb_stock=read_site_specific_obs("Cagb_stock_gCm2",infile)
        Cagb_stock_unc=read_site_specific_obs("Cagb_stock_unc_gCm2",infile)
        if (length(Cagb_stock_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cagb_stock_unc = rep(-9999,times = length(Cagb_stock))
            Cagb_stock_unc[which(Cagb_stock != -9999)] = abs(0.25 * Cagb_stock[which(Cagb_stock != -9999)])
        }
    } else {
        # assume no data available
        Cagb_stock = -9999 ; Cagb_stock_unc = -9999
    }
    # apply lower bound in all cases to the uncertainty
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    Cagb_stock_unc[Cagb_stock_unc >= 0] = pmax(100,sqrt(Cagb_stock_unc[Cagb_stock_unc >= 0]**2 + (0.1*mean(Cagb_stock[Cagb_stock >= 0]))**2))

    ###
    ## Get some Croots information (stock)
    ###

    if (Croots_stock_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Croots_stock=read_site_specific_obs("Croots_stock_gCm2",infile)
        Croots_stock_unc=read_site_specific_obs("Croots_stock_unc_gCm2",infile)
        if (length(Croots_stock_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Croots_stock_unc = rep(-9999,times = length(Croots_stock))
            Croots_stock_unc[which(Croots_stock != -9999)] = abs(0.44 * Croots_stock[which(Croots_stock != -9999)])
        }
    } else {
        # assume no data available
        Croots_stock = -9999 ; Croots_stock_unc = -9999
    }

    ###
    ## Get some Clitter information (stock)
    ###
    #    print("checking Clitter")
    if (Clit_stock_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Clit_stock=read_site_specific_obs("Clit_stock_gCm2",infile)
        Clit_stock_unc=read_site_specific_obs("Clit_stock_unc_gCm2",infile)
        if (length(Clit_stock_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
           Clit_stock_unc = rep(-9999,times = length(Clit_stock))
           Clit_stock_unc[which(Clit_stock != -9999)] = abs(0.38 * Clit_stock[which(Clit_stock != -9999)])
        }
    } else {
        # assume no data available
        Clit_stock = -9999 ; Clit_stock_unc = -9999
    }

    ###
    ## Get some Csom information (stock)
    ###

    if (Csom_stock_source == "site_specific") {
      infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
      Csom_stock=read_site_specific_obs("Csom_stock_gCm2",infile)
      Csom_stock_unc=read_site_specific_obs("Csom_stock_unc_gCm2",infile)
      if (length(Csom_stock_unc) == 1) {
        # on the other hand if not then we have no uncertainty info, so use default
        Csom_stock_unc = rep(-9999,times = length(Csom_stock))
        Csom_stock_unc[which(Csom_stock != -9999)] = abs(0.24 * Csom_stock[which(Csom_stock != -9999)])
      }
    } else {
      # assume no data available
      Csom_stock = -9999 ; Csom_stock_unc = -9999
    }
    # Now assuming we have actual information we need to add the model structural uncertainty.
    # A structural of uncertainty has been estimates at ~ 1000 gC/m2 based on Smallman et al., (2017)
    #Csom_stock_unc[Csom_stock_unc >= 0] = sqrt(Csom_stock_unc[Csom_stock_unc >= 0]**2 + 1000**2)
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    Csom_stock_unc[Csom_stock_unc >= 0] = pmax(1000,sqrt(Csom_stock_unc[Csom_stock_unc >= 0]**2 + (0.1*mean(Csom_stock[Csom_stock > 0]))**2))

    ###
    ## Get some Ccoarseroot information (stock)
    ###

    if (Ccoarseroot_stock_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Ccoarseroot_stock=read_site_specific_obs("Ccoarseroot_stock_gCm2",infile)
        Ccoarseroot_stock_unc=read_site_specific_obs("Ccoarseroot_stock_unc_gCm2",infile)
      if (length(Ccoarseroot_stock_unc) == 1) {
          # on the other hand if not then we have no uncertainty info, so use default
          Ccoarseroot_stock_unc = rep(-9999,times = length(Ccoarseroot_stock))
          Ccoarseroot_stock_unc[which(Ccoarseroot_stock != -9999)] = abs(0.24 * Ccoarseroot_stock[which(Ccoarseroot_stock != -9999)])
      }
    } else {
        # assume no data available
        Ccoarseroot_stock = -9999 ; Ccoarseroot_stock_unc = -9999
    }

    ###
    ## Get some Cfolmax information (stock)
    ###

    if (Cfolmax_stock_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Cfolmax_stock=read_site_specific_obs("Cfolmax_stock_gCm2",infile)
        Cfolmax_stock_unc=read_site_specific_obs("Cfolmax_stock_unc_gCm2",infile)
        if (length(Cfolmax_stock_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cfolmax_stock_unc = rep(-9999,times = length(Cfolmax_stock))
            Cfolmax_stock_unc[which(Cfolmax_stock != -9999)] = abs(0.24 * Cfolmax_stock[which(Cfolmax_stock != -9999)])
        }
    } else {
        # assume no data available
        Cfolmax_stock = -9999 ; Cfolmax_stock_unc = -9999
    }

    ###
    ## Get some deforestation information (fraction time series)
    ###

    if (deforestation_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        deforestation = read_site_specific_obs("deforestation_fraction",infile)
        forest_management = read_site_specific_obs("management_type",infile)
        if (length(forest_management) == 1) {forest_management = rep(2, times = length(deforestation))}
        yield_class = -9999 #read_site_specific_obs("yield_class",infile)
        age = read_site_specific_obs("age",infile)
        if (length(age) > 1) {age = age[1]} # we only want the age at the beginning of the simulation
    } else if (deforestation_source == "GFW") {
        output = extract_forestry_information(grid_long_loc,grid_lat_loc,timestep_days,
                                              spatial_type,resolution,grid_type,latlon_wanted,
                                              forest_all,start_year,end_year,ctessel_pft,
                                              years_to_load,doy_obs)
        ctessel_pft = output$ctessel_pft
        deforestation = output$deforestation
        yield_class = output$yield_class
        age = output$age
        forest_management = 2
    } else {
        # assume no data available
        deforestation = 0
        forest_management = 2
        yield_class = 0
        age = -9999
    }

    ###
    ## Get some burnt area information (fraction time series)
    ###

    if (burnt_area_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        burnt_area = read_site_specific_obs("burnt_area_fraction",infile)
    } else if (burnt_area_source == " "){
        # assume no data available
        burnt_area = 0
    } else {
        # Assume all burnt area product follow the same structure
        burnt_area = extract_burnt_area_information(grid_long_loc,grid_lat_loc,latlon_wanted,
                                                    timestep_days,spatial_type,grid_type,resolution,
                                                    start_year,end_year,burnt_all,years_to_load,doy_obs)
    }

    ###
    ## Get some snow water equivalent (mm/day)
    ###

    if (snow_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        SWE=read_site_specific_obs("snow_water_mmday",infile)
        SWE_unc=read_site_specific_obs("snow_water_unc_mmday",infile)
        if (length(which(SWE_unc != -9999)) == 0) {
            # on the other hand if not then we have no uncertainty info, so use default
           SWE_unc=rep(sd(SWE[which(SWE != -9999)],na.rm=TRUE),length.out=length(SWE))
        }
    } else {
        # assume no data available
        SWE = -9999 ; SWE_unc = -9999
    }

    ###
    ## Get initial soil water fraction prior (initial conditions)
    ###

    if (soilwater_initial_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        soilwater=read_site_specific_obs("soil_water",infile)
        soilwater_unc=read_site_specific_obs("soil_water_unc",infile)
        if (soilwater_unc == -9999 & soilwater > 0) {
          # on the other hand if not then we have no uncertainty info, so use default
          soilwater_unc = 0.10 * soilwater
        }
    } else if (soilwater_initial_source == "GLEAM") {
        output = extract_soilwater_initial(spatial_type,resolution,grid_type,latlon_wanted,soilwater_all)
        soilwater = output$soil_water ; soilwater_unc = output$soil_water_unc
    } else {
        # assume no data available
        soilwater = -9999 ; soilwater_unc = -9999
    }

    ###
    ## Get some Cwood information (potential stock)
    ###

    if (Cwood_potential_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        Cwood_potential=read_site_specific_obs("Cwood_potential_gCm2",infile)
        Cwood_potential_unc=read_site_specific_obs("Cwood_potential_unc_gCm2",infile)
    } else if (Cwood_potential_source == "UoE_potAGB") {
        # get Cwood
        output = extract_Cwood_potential(grid_long_loc,grid_lat_loc,timestep_days,
                                         spatial_type,resolution,grid_type,latlon_wanted,
                                         Cwood_potential_all)
        Cwood_potential = output$Cwood_stock
        Cwood_potential_unc = output$Cwood_stock_unc
    } else {
        # assume no data available
        Cwood_potential=-9999 ; Cwood_potential_unc=-9999
    }

    ###
    ## Get some leaf carbon per unit leaf area (gC/m2) information (potential stock)
    ###

    if (lca_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        lca=read_site_specific_obs("LCA_gCm2",infile)
        lca_unc=read_site_specific_obs("LCA_unc_gCm2",infile)
    } else if (lca_source == "Butler") {
        # get Cwood
        output = extract_lca_prior(grid_long_loc,grid_lat_loc,spatial_type,resolution,
                                   grid_type,latlon_wanted,lca_all)
        lca = output$lca_gCm2
        lca_unc = output$lca_unc_gCm2
    } else {
        # assume no data available
        lca = -9999 ; lca_unc = -9999
    }
	
	 ###
    ## Get leaf_fall_period information
    ###

    if (leaf_fall_period_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        leaf_fall_period=read_site_specific_obs("leaf_fall_period",infile)
        leaf_fall_period_unc=read_site_specific_obs("leaf_fall_period_unc",infile)
    } else if (leaf_fall_period_source == "TG") {
        # 
        output = extract_leaf_fall_period_prior(grid_long_loc,grid_lat_loc,spatial_type,resolution,
                                   grid_type,latlon_wanted,leaf_fall_period_all)
        leaf_fall_period = output$leaf_fall_period
        #leaf_fall_period_unc = output$leaf_fall_period_unc_days
    } else {
        # assume no data available
        leaf_fall_period = -9999 #; leaf_fall_period_unc = -9999
    }
	
	
	
	    ###
    ## Get lifespan/MRT ("natural") fol information
    ###

    if (lifespan_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        lifespan=read_site_specific_obs("lifespan",infile)
        #lifespan_unc=read_site_specific_obs("lifespan_unc",infile)
    } else if (lifespan_source == "Tupek") {
        # 
        output = extract_lifespan_prior(grid_long_loc,grid_lat_loc,spatial_type,resolution,
                                   grid_type,latlon_wanted,lifespan_all)
        lifespan = output$lifespan
        #lifespan_unc = output$lifespan_unc_yrs
    } else {
        # assume no data available
        lifespan = -9999 #; lifespan_unc = -9999
    }
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    #lca_unc[lca_unc >= 0] = pmax(10,sqrt(lca_unc[lca_unc >= 0]**2 + (0.1*mean(lca[lca > 0]))**2))

    # return output now
    return(list(LAT = latlon_wanted[1], LAI = lai, LAI_unc = lai_unc, GPP = GPP, GPP_unc = GPP_unc, Fire = Fire, Fire_unc = Fire_unc
               ,Evap = Evap, Evap_unc = Evap_unc, NEE = NEE, NEE_unc = NEE_unc, Reco = Reco, Reco_unc = Reco_unc
               ,Cfol_stock = Cfol_stock, Cfol_stock_unc = Cfol_stock_unc
               ,Cwood_stock = Cwood_stock, Cwood_stock_unc = Cwood_stock_unc, Cagb_stock=Cagb_stock, Cagb_stock_unc = Cagb_stock_unc
               ,Croots_stock = Croots_stock, Croots_stock_unc = Croots_stock_unc, Clit_stock = Clit_stock, Clit_stock_unc = Clit_stock_unc
               ,Csom_stock = Csom_stock, Csom_stock_unc = Csom_stock_unc, Ccoarseroot_stock = Ccoarseroot_stock
               ,Ccoarseroot_stock_unc = Ccoarseroot_stock_unc, Cfolmax_stock = Cfolmax_stock, Cfolmax_stock_unc = Cfolmax_stock_unc
               ,Csom_initial = Csom_initial, Csom_initial_unc = Csom_initial_unc, Cfol_initial = Cfol_initial, Cfol_initial_unc = Cfol_initial_unc
               ,Cwood_initial = Cwood_initial, Cwood_initial_unc = Cwood_initial_unc, Croots_initial = Croots_initial
               ,Croots_initial_unc = Croots_initial_unc, Clit_initial = Clit_initial, Clit_initial_unc = Clit_initial_unc
               ,deforestation = deforestation, burnt_area = burnt_area, ctessel_pft = ctessel_pft, yield_class = yield_class
               ,age = age, forest_management = forest_management, top_sand = top_sand, bot_sand = bot_sand, top_clay = top_clay
               ,bot_clay = bot_clay, plant = plant, plant_range = plant_range, harvest = harvest, harvest_range = harvest_range
               ,SWE = SWE, SWE_unc = SWE_unc, soilwater = soilwater, soilwater_unc = soilwater_unc, nbe = nbe, nbe_unc = nbe_unc
               ,Cwood_potential = Cwood_potential, Cwood_potential_unc = Cwood_potential_unc, lca = lca, lca_unc = lca_unc
               ,Cwood_inc = Cwood_inc, Cwood_inc_unc = Cwood_inc_unc, Cwood_inc_lag = Cwood_inc_lag
               ,Cwood_mortality = Cwood_mortality, Cwood_mortality_unc = Cwood_mortality_unc, Cwood_mortality_lag = Cwood_mortality_lag
			   , lifespan = lifespan, leaf_fall_period = leaf_fall_period))#, lifespan_unc = lifespan_unc))


} # end function extract_obs

## Use byte compile
extract_obs<-cmpfun(extract_obs)
