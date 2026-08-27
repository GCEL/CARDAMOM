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
# Controls the extraction of location specific observations from different datasets
# 
# This function is based on an original Matlab function development by A. A. Bloom 
# (UoE, now at the Jet Propulsion Laboratory). Translation to R and subsequent 
# modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).
#
#########################################################################################

extract_obs<-function(grid_long_loc,grid_lat_loc,latlon_wanted,lai_all,Csom_all,forest_all
                     ,Cwood_initial_all,Cwood_stock_all,Cwood_potential_all
                     ,sand_clay_all,crop_man_all,burnt_all,soilwater_all,nbe_all
                     ,lca_all,gpp_all,Cwood_change_all,Cwood_growth_all,Cwood_loss_all
                     ,fire_all,dlai_all,fapar_all,et_all,RhetQ10_all,MTTsom_all,MaxRootDepth_all
                     ,Cagb_stock_all,Cagb_change_all,Cagb_growth_all,leaflifespan_all
                     ,labile_release_timing_all,labile_release_period_all,frac_Cwood_coarse_root_all
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
    ## Extract the local information for timeseries information with uncertainty,
    ## i.e. those values which are assimilated
    ###

    ###
    ## Get some NBE information (gC/m2/day); negative is sink

    if (nbe_source == "Gridded_nc" | nbe_source == "Gridded_tif") {

        if (nbe_all$data_available) {
            # Extract NBE and uncertainty information
            # NOTE: assume default uncertainty (+/- scale)
            # Extract the current location from the gridded dataset
            output = extract_timeseries_observations_with_uncertainty(grid_long_loc,grid_lat_loc,timestep_days,years_to_load,doy_obs,
                                                                      nbe_all,agg_func = "mean",na_flag = -9999,
                                                                      est_var_name_in = "nbe_gCm2day",
                                                                      unc_var_name_in = "nbe_unc_gCm2day",
                                                                      lag_var_name_in = "nbe_lag_day",
                                                                      est_var_name_out = "nbe",
                                                                      unc_var_name_out = "nbe_unc",
                                                                      lag_var_name_out = "nbe_lag")
            # Assign to local variables                           
            nbe = output$nbe ; nbe_unc = output$nbe_unc ; nbe_lag = output$nbe_lag
        } else {
            # Set missing data value
            nbe = -9999 ; nbe_unc = -9999 ; nbe_lag = -9999
        }
        
    } else if (nbe_source == "site_specific") {

        # read from .csv or netcdf
        infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        nbe = read_site_specific_obs("NBE_gCm2day",infile) 
        nbe_unc = read_site_specific_obs("NBE_unc_gCm2day",infile)
        nbe_lag = read_site_specific_obs("NBE_lag_step",infile)
        if (max(nbe_unc) == -9999) {
            nbe_unc = rep(-9999,times = length(nbe))
            # apply default uncertainty consistent with Eddy covariance estimates of NEE, 
            # Composed of NEE 0.58 gC/m2/day (Hill et al., 2012) plus mass balance mismatch of
            nbe_unc[which(nbe != -9999)] = 0.58
        }
        if (max(nbe_lag) == -9999) {
            nbe_lag = rep(0,times = length(nbe))
        }

    } else {

        nbe = -9999 ; nbe_unc = -9999 ; nbe_lag = -9999

    }
    # Add model structural uncertainty to the uncertainty estimate if present
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    nbe_unc[nbe_unc >= 0] = sqrt(nbe_unc[nbe_unc >= 0]**2 + 1.0**2)

    ###
    ## Get some LAI information (m2/m2)

    if (lai_source == "Gridded_tif" | lai_source == "Gridded_nc") {
        if (lai_all$data_available) {
            # Extract lai and uncertainty information
            # NOTE: assume default uncertainty (+/- scale)
            # Extract the current location from the gridded dataset
            output = extract_timeseries_observations_with_uncertainty(grid_long_loc,grid_lat_loc,timestep_days,years_to_load,doy_obs,
                                                                      lai_all,agg_func = "mean", na_flag = -9999,
                                                                      est_var_name_in = "lai_m2m2",
                                                                      unc_var_name_in = "lai_unc_m2m2",
                                                                      lag_var_name_in = "lai_lag_day",
                                                                      est_var_name_out = "lai",
                                                                      unc_var_name_out = "lai_unc",
                                                                      lag_var_name_out = "lai_lag")
            # Assign to local variables                           
            lai = output$lai ; lai_unc = output$lai_unc ; lai_lag = output$lai_lag
        } else {
            # Set missing data value
            lai = -9999 ; lai_unc = -9999 ; lai_lag = -9999
        }
    } else if (lai_source == "site_specific") {
        # read from .csv or netcdf
        infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        lai = read_site_specific_obs("LAI_m2m2",infile) 
        lai_unc = read_site_specific_obs("LAI_unc_m2m2",infile)
        lai_lag = read_site_specific_obs("LAI_lag_step",infile)
        if (max(lai_unc) == -9999) {
            lai_unc = rep(-9999,times = length(lai))
            # apply default uncertainty
            lai_unc[which(lai != -9999)] = 0.25
        }
        if (max(lai_lag) == -9999) {
            lai_lag = rep(1,times = length(lai))
        }        
    } else {
        # Set missing data value
        lai = -9999 ; lai_unc = -9999 ; lai_lag = -9999
    }
    # Assume minimum uncertainty to reflect model structural uncertainty
    # Estimates from comparison of LAI uncertainties trials at 0.5 and 0.25,
    # resultant CI in both instances is range of ~0.50. Therefore CI of +/- 0.25
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    if (length(which(lai_unc >= 0)) > 0) {
        lai_unc[lai_unc >= 0] = pmax(0.25,sqrt(lai_unc[lai_unc >= 0]**2 + (0.1*mean(lai[lai >= 0]))**2))
    }

    ###
    ## Get some fAPAR information (0-1)

    if (fapar_source == "Gridded_nc" | fapar_source == "Gridded_tif") {
        if (fapar_all$data_available) {
            # Extract fapar and uncertainty information
            # NOTE: assume default uncertainty (+/- scale)
            # Extract the current location from the gridded dataset
            output = extract_timeseries_observations_with_uncertainty(grid_long_loc,grid_lat_loc,timestep_days,years_to_load,doy_obs,
                                                                      fapar_all,agg_func = "mean", na_flag = -9999,
                                                                      est_var_name_in = "fapar",
                                                                      unc_var_name_in = "fapar_unc",
                                                                      lag_var_name_in = "fapar_lag_day",
                                                                      est_var_name_out = "fapar",
                                                                      unc_var_name_out = "fapar_unc",
                                                                      lag_var_name_out = "fapar_lag")
            # Assign to local variables                           
            fapar = output$fapar ; fapar_unc = output$fapar_unc ; fapar_lag = output$fapar_lag
        } else {
            # Set missing data value
            fapar = -9999 ; fapar_unc = -9999 ; fapar_lag = -9999
        }
    } else if (fapar_source == "site_specific") {
        # read from .csv or netcdf
        infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        fapar = read_site_specific_obs("fAPAR_fraction",infile) 
        fapar_unc = read_site_specific_obs("fAPAR_unc_fraction",infile)
        fapar_lag = read_site_specific_obs("fAPAR_lag_step",infile)
        if (max(fapar_unc) == -9999) {
            fapar_unc = rep(-9999,times = length(fapar))
            # apply default uncertainty
            fapar_unc[which(fapar != -9999)] = 0.05
        }
        if (max(fapar_lag) == -9999) {
            fapar_lag = rep(0,times = length(fapar))
        }        
    } else {
        # Set missing data value
        fapar = -9999 ; fapar_unc = -9999 ; fapar_lag = -9999
    }
    # Assume minimum uncertainty to reflect model structural uncertainty.
    # NOTE minimum uncertainty bound irrespective of the dataset estimates,
    # further more we currently have no good estimate of a realistic minimum 
    # uncertainty for FAPAR. The value used here is arbitary.
    if (length(which(fapar_unc >= 0)) > 0) {
        fapar_unc[fapar_unc >= 0] = pmax(0.05,sqrt(fapar_unc[fapar_unc >= 0]**2 + (0.1*mean(fapar[fapar >= 0]))**2))
    }

    ###
    ## Get some Cfoliage information (stock; gC/m2)

    if (Cfol_stock_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Cfol_stock = read_site_specific_obs("Cfol_stock_gCm2",infile)
        Cfol_stock_unc = read_site_specific_obs("Cfol_stock_unc_gCm2",infile)
        Cfol_stock_lag = read_site_specific_obs("Cfol_stock_lag_step",infile)
        if (length(Cfol_stock_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cfol_stock_unc = rep(-9999,times = length(Cfol_stock))
            # See Smallman et al., (2017) for uncertainty estimate
            Cfol_stock_unc[which(Cfol_stock > 0)] = 0.38 * Cfol_stock[which(Cfol_stock > 0)]
        }
        if (length(Cfol_stock_lag) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cfol_stock_lag = rep(1,times = length(Cfol_stock))
        }        
    } else {
        # assume no data available
        Cfol_stock = -9999 ; Cfol_stock_unc = -9999 ; Cfol_stock_lag = -9999
    }

    ###
    ## Get some information on C extracted due to harvest
    ## This can be either crop yield, grassland cutting or forest loss
    ## Specificially related to C removed from the site (horizontal transfer), 
    ## not that which remains as litter.
    ## (gC/m2/day; time series)

    if (harvest_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        harvest = read_site_specific_obs("harvest_gCm2day",infile)
        harvest_unc = read_site_specific_obs("harvest_uncertainty_gCm2day",infile)
        harvest_lag = read_site_specific_obs("harvest_lag_step",infile) # in model time steps
        if (length(harvest) == 1) {
            stop("Timeseries of harvest information was expected (harvest_source == 'site_specific') but not provided")
        }
        # Has uncertainty information been provided?
        if (length(harvest_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            harvest_unc = rep(-9999,times = length(harvest))
            harvest_unc[which(harvest > 0)] = 0.25 * harvest[which(harvest > 0)]
        }
        # Has lag information been provided
        if (length(harvest_lag) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            harvest_lag = rep(-9999,times = length(harvest))
            harvest_lag[which(harvest > 0)] = 0 # assume applies to current time step only
        }
    } else {
        harvest = -9999              # Extracted C due to harvest over lag period (gC/m2/day)
        harvest_unc = -9999          # Extracted C due to harvest varince
        harvest_lag = -9999          # Lag period over which to average (steps)
    }

    ###
    ## Get some Wood net change information (gC/m2/day; time series)

    if (Cwood_change_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Cwood_change = read_site_specific_obs("Cwood_change_gCm2day",infile)
        Cwood_change_unc = read_site_specific_obs("Cwood_change_uncertainty_gCm2day",infile)
        Cwood_change_lag = read_site_specific_obs("Cwood_change_lag_step",infile) # in model time steps
        # Has uncertainty information been provided?
        if (length(Cwood_change_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cwood_change_unc = rep(-9999,times = length(Cwood_change))
            Cwood_change_unc[which(Cwood_change > 0)] = 0.25 * Cwood_change[which(Cwood_change > 0)]
        }
        # Has lag information been provided
        if (length(Cwood_change_lag) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cwood_change_lag = rep(-9999,times = length(Cwood_change))
            Cwood_change_lag[which(Cwood_change > 0)] = 0 # assume applies to current time step only
        }
    } else if (Cwood_change_source == "Gridded_nc" | Cwood_change_source == "Gridded_tif") {
        # If there are any values in the analysis window
        if (Cwood_change_all$data_available) {        
            # Extract gross wood production / growth and uncertainty information
            # NOTE: assume default uncertainty (+/- scale)
            # Extract the current location from the gridded dataset
            output = extract_timeseries_observations_with_uncertainty(grid_long_loc,grid_lat_loc,timestep_days,years_to_load,doy_obs,
                                                                      Cwood_change_all,agg_func = "mean", na_flag = -9999,
                                                                      est_var_name_in = "Cwood_change_gCm2day",
                                                                      unc_var_name_in = "Cwood_change_uncertainty_gCm2day",
                                                                      lag_var_name_in = "Cwood_change_lag",
                                                                      est_var_name_out = "Cwood_change",
                                                                      unc_var_name_out = "Cwood_change_unc",
                                                                      lag_var_name_out = "Cwood_change_lag")                                                  
            Cwood_change = output$Cwood_change ; Cwood_change_unc = output$Cwood_change_unc ; Cwood_change_lag = output$Cwood_change_lag
            # Tidy up
            rm(output)
        } else {
            # assume no data available
            Cwood_change = -9999 ; Cwood_change_unc = -9999 ; Cwood_change_lag = -9999
        }
    } else {
        # assume no data available
        Cwood_change = -9999 ; Cwood_change_unc = -9999 ; Cwood_change_lag = -9999
    }
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    Cwood_change_unc[Cwood_change_unc >= 0] = pmax(0.1,sqrt(Cwood_change_unc[Cwood_change_unc >= 0]**2 + (0.1*mean(Cwood_change[Cwood_change_unc >= 0]))**2))

    ###
    ## Get some gross woody production information (gC/m2/day; time series)

    if (Cwood_growth_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Cwood_growth = read_site_specific_obs("Cwood_production_gCm2day",infile)
        Cwood_growth_unc = read_site_specific_obs("Cwood_production_uncertainty_gCm2day",infile)
        Cwood_growth_lag = read_site_specific_obs("Cwood_production_lag_step",infile) # in model time steps
        # Has uncertainty information been provided?
        if (length(Cwood_growth_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cwood_growth_unc = rep(-9999,times = length(Cwood_growth))
            Cwood_growth_unc[which(Cwood_growth > 0)] = 0.25 * Cwood_growth[which(Cwood_growth > 0)]
        }
        # Has lag information been provided
        if (length(Cwood_growth_lag) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cwood_growth_lag = rep(-9999,times = length(Cwood_growth))
            Cwood_growth_lag[which(Cwood_growth > 0)] = 0 # assume applies to current time step only
        }
    } else if (Cwood_growth_source == "Gridded_nc" | Cwood_growth_source == "Gridded_tif") {
        # If there are any values in the analysis window
        if (Cwood_growth_all$data_available) {      
            # Extract gross wood production / growth and uncertainty information
            # NOTE: assume default uncertainty (+/- scale)
            # Extract the current location from the gridded dataset
            output = extract_timeseries_observations_with_uncertainty(grid_long_loc,grid_lat_loc,timestep_days,years_to_load,doy_obs,
                                                                      Cwood_growth_all,agg_func = "mean", na_flag = -9999,
                                                                      est_var_name_in = "Cwood_growth_gCm2day",
                                                                      unc_var_name_in = "Cwood_growth_uncertainty_gCm2day",
                                                                      lag_var_name_in = "Cwood_growth_lag",
                                                                      est_var_name_out = "Cwood_growth",
                                                                      unc_var_name_out = "Cwood_growth_unc",
                                                                      lag_var_name_out = "Cwood_growth_lag")   
            # Assign to local variables                           
            Cwood_growth = output$Cwood_growth ; Cwood_growth_unc = output$Cwood_growth_unc ; Cwood_growth_lag = output$Cwood_growth_lag
            # Tidy up
            rm(output)
        } else {
            # assume no data available
            Cwood_growth = -9999 ; Cwood_growth_unc = -9999 ; Cwood_growth_lag = -9999
        }
    } else {
        # assume no data available
        Cwood_growth = -9999 ; Cwood_growth_unc = -9999 ; Cwood_growth_lag = -9999
    }
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    Cwood_growth_unc[Cwood_growth_unc >= 0] = pmax(0.1,sqrt(Cwood_growth_unc[Cwood_growth_unc >= 0]**2 + (0.1*mean(Cwood_growth[Cwood_growth_unc >= 0]))**2))

    ###
    ## Get some Wood loss information (gC/m2/day; time series)

    if (Cwood_mortality_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Cwood_loss = read_site_specific_obs("Cwood_loss_gCm2day",infile)
        Cwood_loss_unc = read_site_specific_obs("Cwood_loss_uncertainty_gCm2day",infile)
        Cwood_loss_lag = read_site_specific_obs("Cwood_loss_lag_step",infile) # in model time steps
        # Has uncertainty information been provided?
        if (length(Cwood_loss_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cwood_loss_unc = rep(-9999,times = length(Cwood_loss))
            Cwood_loss_unc[which(Cwood_loss > 0)] = 0.25 * Cwood_loss[which(Cwood_loss > 0)]
        }
        # Has lag information been provided
        if (length(Cwood_loss_lag) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cwood_loss_lag = rep(-9999,times = length(Cwood_loss))
            Cwood_loss_lag[which(Cwood_loss > 0)] = 0 # assume applies to current time step only
        }
    } else if (Cwood_mortality_source == "Gridded_nc" | Cwood_mortality_source == "Gridded_tif") {
        # If there are any values in the analysis window
        if (Cwood_loss_all$data_available) {      
            # Extract wood mortality and uncertainty information
            # NOTE: assume default uncertainty (+/- scale)
            # Extract the current location from the gridded dataset
            output = extract_timeseries_observations_with_uncertainty(grid_long_loc,grid_lat_loc,timestep_days,years_to_load,doy_obs,
                                                                      Cwood_loss_all,agg_func = "mean", na_flag = -9999,
                                                                      est_var_name_in = "Cwood_loss_gCm2day",
                                                                      unc_var_name_in = "Cwood_loss_uncertainty_gCm2day",
                                                                      lag_var_name_in = "Cwood_loss_lag",
                                                                      est_var_name_out = "Cwood_loss",
                                                                      unc_var_name_out = "Cwood_loss_unc",
                                                                      lag_var_name_out = "Cwood_loss_lag")   
            # Assign to local variables                           
            Cwood_loss = output$Cwood_loss ; Cwood_loss_unc = output$Cwood_loss_unc ; Cwood_loss_lag = output$Cwood_loss_lag

            # Tidy up
            rm(output)
        } else {
            # assume no data available
            Cwood_loss = -9999 ; Cwood_loss_unc = -9999 ; Cwood_loss_lag = -9999
        }
    } else {
        # assume no data available
        Cwood_loss = -9999 ; Cwood_loss_unc = -9999 ; Cwood_loss_lag = -9999
    }
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    Cwood_loss_unc[Cwood_loss_unc >= 0] = pmax(0.1,sqrt(Cwood_loss_unc[Cwood_loss_unc >= 0]**2 + 
                                                                  (0.1*mean(Cwood_loss[Cwood_loss_unc >= 0]))**2))

    ###
    ## Get some foliage to litter flux information (gC/m2/day; time series)

    if (foliage_to_litter_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        foliage_to_litter = read_site_specific_obs("foliage_to_litter_gCm2day",infile)
        foliage_to_litter_unc = read_site_specific_obs("foliage_to_litter_unc_gCm2day",infile)
        foliage_to_litter_lag = read_site_specific_obs("foliage_to_litter_lag_step",infile) # in model time steps
        # Has uncertainty information been provided?
        if (length(foliage_to_litter_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            foliage_to_litter_unc = rep(-9999,times = length(foliage_to_litter))
            foliage_to_litter_unc[which(foliage_to_litter > 0)] = 0.25 * foliage_to_litter[which(foliage_to_litter > 0)]
        }
        # Has lag information been provided
        if (length(foliage_to_litter_lag) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            foliage_to_litter_lag = rep(-9999,times = length(foliage_to_litter))
            foliage_to_litter_lag[which(foliage_to_litter > 0)] = 0 # assume applies to current time step only
        }
    } else {
        # assume no data available
        foliage_to_litter = -9999 ; foliage_to_litter_unc = -9999 ; foliage_to_litter_lag = -9999
    }
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    foliage_to_litter_unc[foliage_to_litter_unc >= 0] = pmax(0.1,sqrt(foliage_to_litter_unc[foliage_to_litter_unc >= 0]**2 + 
                                                                      (0.1*mean(foliage_to_litter[foliage_to_litter_unc >= 0]))**2))

    ###
    ## Get some GPP information (time series; gC/m2/day)

    if (gpp_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        GPP = read_site_specific_obs("GPP_gCm2day",infile)
        GPP_unc = read_site_specific_obs("GPP_unc_gCm2day",infile)
        GPP_lag = read_site_specific_obs("GPP_lag_step",infile)
        if (length(GPP_unc) == 1) {
            GPP_unc = rep(-9999,times = length(GPP))
            # Composed of NEE 0.58 gC/m2/day (Hill et al., 2012) plus mass balance mismatch of
            # 0.16 gC/m2/day, therefore 0.74 gC/m2/day
            GPP_unc[which(GPP > 0)] = 0.74
        }
        if (length(GPP_lag) == 1) {
            GPP_lag = rep(0,times = length(GPP))
        }
    } else if (gpp_source == "Gridded_nc" | gpp_source == "Gridded_tif") {

        if (gpp_all$data_available) {
            # Extract GPP and uncertainty information
            # NOTE: assume default uncertainty (+/- scale)
            # Extract the current location from the gridded dataset
            output = extract_timeseries_observations_with_uncertainty(grid_long_loc,grid_lat_loc,timestep_days,years_to_load,doy_obs,
                                                                      gpp_all,agg_func = "mean", na_flag = -9999,
                                                                      est_var_name_in = "gpp_gCm2day",
                                                                      unc_var_name_in = "gpp_unc_gCm2day",
                                                                      lag_var_name_in = "gpp_lag_day",
                                                                      est_var_name_out = "GPP",
                                                                      unc_var_name_out = "GPP_unc",
                                                                      lag_var_name_out = "GPP_lag")
            # Assign to local variables                           
            GPP = output$GPP ; GPP_unc = output$GPP_unc ; GPP_lag = output$GPP_lag
        } else {
            # Set missing data value
            GPP = -9999 ; GPP_unc = -9999 ; GPP_lag = -9999
        }

    } else {
        # assume no data available
        GPP = -9999 ; GPP_unc = -9999 ; GPP_lag = -9999
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

    if (fire_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Fire = read_site_specific_obs("fire_gCm2day",infile)
        Fire_unc = read_site_specific_obs("fire_unc_gCm2day",infile)
        Fire_lag = read_site_specific_obs("fire_lag_step",infile)
        if (length(Fire_unc) == 1) {
            Fire_unc = rep(-9999,times = length(Fire))
            # Ill defined assumption
            Fire_unc[which(Fire > 0)] = 0.1
        }
        if (length(Fire_lag) == 1) {
            Fire_lag = rep(0,times = length(Fire))
        }        
    } else if (fire_source == "Gridded_nc" | fire_source == "Gridded_tif") {
        if (fire_all$data_available) {
            # Extract fire and uncertainty information
            # NOTE: assume default uncertainty (+/- scale)
            # Extract the current location from the gridded dataset
            output = extract_timeseries_observations_with_uncertainty(grid_long_loc,grid_lat_loc,timestep_days,years_to_load,doy_obs,
                                                                      fire_all,agg_func = "mean", na_flag = -9999,
                                                                      est_var_name_in = "fire_gCm2day",
                                                                      unc_var_name_in = "fire_unc_gCm2day",
                                                                      lag_var_name_in = "fire_lag_day",
                                                                      est_var_name_out = "Fire",
                                                                      unc_var_name_out = "Fire_unc",
                                                                      lag_var_name_out = "Fire_lag")
            # Assign to local variables                           
            Fire = output$Fire ; Fire_unc = output$Fire_unc ; Fire_lag = output$Fire_lag
        } else {
            # Set missing data value
            Fire = -9999 ; Fire_unc = -9999 ; Fire_lag = -9999
        }
    } else {
        # assume no data available
        Fire = -9999 ; Fire_unc = -9999 ; Fire_lag = -9999
    }
    # Combine with an estimate of model structural error.
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates, based on median difference between GFED and GFAS
    #Fire_unc[Fire_unc >= 0] = pmax(0.1,sqrt(Fire_unc[Fire_unc >= 0]**2 + (0.1*mean(Fire[Fire >= 0]))**2))
    #Fire_unc[Fire_unc >= 0] = pmax(0.1*mean(Fire[Fire >= 0]),pmax(0.01,pmin(Fire[Fire >= 0],Fire_unc[Fire_unc >= 0])))
    Fire_unc[Fire_unc >= 0] = pmax(0.1,sqrt(Fire_unc[Fire_unc >= 0]**2 + (0.1*mean(Fire[Fire >= 0]))**2))
    # Fire emissions are dominated by zero (or effectively zero) values. This bias' the AP-MCMC (or MHMCMC)
    # away from fitting the more important emission values in favour of the more common zero (or near zero values.)
    # To address this we remove fire emission observations less than 0.01 gC/m2/day
    if (use_parallel == FALSE) {print("Fire observations (if present) < 0.01 gC/m2/day have been removed to prevent bias to MDF calibration")}
    Fire_unc[Fire < 0.01] = -9999 ; Fire[Fire < 0.01] = -9999

    ###
    ## Get some Evapotranspiration information (time series; kgH2O/m2/day)

    if (et_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        if (modelname == "ACM") {infile=paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")}
        ET = read_site_specific_obs("Evap_kgH2Om2day",infile)
        ET_unc = read_site_specific_obs("Evap_unc_kgH2Om2day",infile)
        ET_lag = read_site_specific_obs("Evap_lag_step",infile)
        if (length(ET_unc) == 1) {
            ET_unc = rep(-9999,times = length(ET))
            ET_unc[which(ET > -9999)] = 0.77 # Assuming Hollinger & Richardson (2005) Tree Physiology, 25, 873-885
        }
        if (length(ET_lag) == 1) {
            ET_lag = rep(0,times = length(ET))
        }        
    } else if (et_source == "Gridded_nc" | et_source == "Gridded_tif") {

        if (et_all$data_available) {
            # Extract ET and uncertainty information
            # NOTE: assume default uncertainty (+/- scale)
            # Extract the current location from the gridded dataset
            output = extract_timeseries_observations_with_uncertainty(grid_long_loc,grid_lat_loc,timestep_days,years_to_load,doy_obs,
                                                                      et_all,agg_func = "mean", na_flag = -9999,
                                                                      est_var_name_in = "et_kgH2Om2day",
                                                                      unc_var_name_in = "et_unc_kgH2Om2day",
                                                                      lag_var_name_in = "et_lag_day",
                                                                      est_var_name_out = "ET",
                                                                      unc_var_name_out = "ET_unc",
                                                                      lag_var_name_out = "ET_lag")
            # Assign to local variables                           
            ET = output$ET ; ET_unc = output$ET_unc ; ET_lag = output$ET_lag
        } else {
            # assume no data available
            ET = -9999 ; ET_unc = -9999 ; ET_lag = -9999
        }      
    } else {
        # assume no data available
        ET = -9999 ; ET_unc = -9999 ; ET_lag = -9999
    }
    # Combine with an estimate of model structural error.
    # A mean rmse of ~1 kgH2O/m2/day was estimated when evaluating ACM-GPP-ET (ACM2)
    # ET against fluxnet observations (Smallman & Williams 2019)
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    ET_unc[ET_unc >= 0] = pmax(1,sqrt(ET_unc[ET_unc >= 0]**2 + (0.1*mean(ET[ET >= 0]))**2))

    ###
    ## Get some Reco information (time series; gC/m2/day)

    if (Reco_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Reco = read_site_specific_obs("Reco_gCm2day",infile)
        Reco_unc = read_site_specific_obs("Reco_unc_gCm2day",infile)
        Reco_lag = read_site_specific_obs("Reco_lag_step",infile)
        if (length(Reco_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Reco_unc = rep(-9999,times = length(Reco))
            # Composed of NEE 0.58 gC/m2/day (Hill et al., 2012) plus mass balance mismatch of
            # 0.16 gC/m2/day, therefore 0.74 gC/m2/day
            Reco_unc[which(Reco > 0)] = 0.74
        }
        if (length(Reco_lag) == 1) {
            Reco_lag = rep(0,length(Reco))
        }
    } else {
        # assume no data available
        Reco = -9999 ; Reco_unc = -9999 ; Reco_lag = -9999
    }
    # Apply model structural and minimum uncertainty
    Reco_unc[Reco_unc >= 0] = pmax(1.0,sqrt(Reco_unc[Reco_unc >= 0]**2 + (0.1*mean(Reco[Reco >= 0]))**2))

    ###
    ## Get some NEE information (time series; gC/m2/day)

    if (nee_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        NEE = read_site_specific_obs("NEE_gCm2day",infile)
        NEE_unc = read_site_specific_obs("NEE_unc_gCm2day",infile)
        NEE_lag = read_site_specific_obs("NEE_lag_step",infile)
        if (length(NEE_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            NEE_unc = rep(-9999,times = length(NEE))
            NEE_unc[which(NEE != -9999)] = 0.58 # Hill et al., (2012)
        }
        if (length(NEE_lag) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            NEE_lag = rep(0,times = length(NEE))
        }        
    } else {
        # assume no data available
        NEE = -9999 ; NEE_unc = -9999 ; NEE_lag = -9999
    }
    # The largest site level mean rmse achieved across the COMPLEX (Famiglietti et al., 2021) sites was 0.99 gC/m2/day.
    # This was achieved in the reduced uncertainty analysis providing and indication of the model structural error
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    NEE_unc[NEE_unc >= 0] = sqrt(NEE_unc[NEE_unc >= 0]**2 + 1**2)   

    ###
    ## Get some Cwood information (stock)

    if (Cwood_stock_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Cwood_stock = read_site_specific_obs("Cwood_stock_gCm2",infile)
        Cwood_stock_unc = read_site_specific_obs("Cwood_stock_unc_gCm2",infile)
        Cwood_stock_lag = read_site_specific_obs("Cwood_stock_lag_step",infile)
        if (length(Cwood_stock_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cwood_stock_unc = rep(-9999,times = length(Cwood_stock))
            Cwood_stock_unc[which(Cwood_stock != -9999)] = abs(0.25 * Cwood_stock[which(Cwood_stock != -9999)])
        }
        if (length(Cwood_stock_lag) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cwood_stock_lag = rep(1,times = length(Cwood_stock))
        }        
    } else if (Cwood_stock_source == "Gridded_nc" | Cwood_stock_source == "Gridded_tif") {
        if (Cwood_stock_all$data_available) {
            # Extract the current location from the gridded dataset
            output = extract_timeseries_observations_with_uncertainty(grid_long_loc,grid_lat_loc,timestep_days,years_to_load,doy_obs,
                                                                      Cwood_stock_all,agg_func = "mean", na_flag = -9999,
                                                                      est_var_name_in = "biomass_gCm2",
                                                                      unc_var_name_in = "biomass_uncertainty_gCm2",
                                                                      lag_var_name_in = "biomass_lag_day",
                                                                      est_var_name_out = "Cwood_stock",
                                                                      unc_var_name_out = "Cwood_stock_unc",
                                                                      lag_var_name_out = "Cwood_stock_lag")
            # Assign to local variables                                                             
            Cwood_stock = output$Cwood_stock ; Cwood_stock_unc = output$Cwood_stock_unc ; Cwood_stock_lag = output$Cwood_stock_lag
            # All maps converted into common format, therefore a common extraction subroutine can be used
#            if (max(Cwood_stock_all$doy_obs) > 0) {
#                output = extract_Cwood_stocks(grid_long_loc,grid_lat_loc,timestep_days,
#                                              spatial_type,resolution,grid_type,latlon_wanted,
#                                              Cwood_stock_all)
#                Cwood_stock = output$Cwood_stock ; Cwood_stock_unc = output$Cwood_stock_unc
#                tmp = which(Cwood_stock > 0) # first AGB only
#                if (length(tmp) > 1) {
#                    Cwood_stock[tmp[-1]] = -9999 ; Cwood_stock_unc[tmp[-1]] = -9999
#                }
#            } else {
#                Cwood_stock = rep(-9999, length(timestep_days))
#                Cwood_stock_unc = rep(-9999, length(timestep_days))
#            }
        } else {
                Cwood_stock = rep(-9999, length(timestep_days))
                Cwood_stock_unc = rep(-9999, length(timestep_days))           
                Cwood_stock_lag = rep(-9999, length(timestep_days))    
        }
    } else {
        # assume no data available
        Cwood_stock = -9999 ; Cwood_stock_unc = -9999 ; Cwood_stock_lag = -9999
    }
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    Cwood_stock_unc[Cwood_stock_unc >= 0] = pmax(500,sqrt(Cwood_stock_unc[Cwood_stock_unc >= 0]**2 + (0.1*mean(Cwood_stock[Cwood_stock >= 0]))**2),na.rm=TRUE)

    ###
    ## Get some Cagb information (stock)

    if (Cagb_stock_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Cagb_stock = read_site_specific_obs("Cagb_stock_gCm2",infile)
        Cagb_stock_unc = read_site_specific_obs("Cagb_stock_unc_gCm2",infile)
        Cagb_stock_lag = read_site_specific_obs("Cagb_stock_lag_step",infile)        
        if (length(Cagb_stock_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cagb_stock_unc = rep(-9999,times = length(Cagb_stock))
            Cagb_stock_unc[which(Cagb_stock != -9999)] = abs(0.25 * Cagb_stock[which(Cagb_stock != -9999)])
        }
        if (length(Cagb_stock_lag) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cagb_stock_lag = rep(1,times = length(Cagb_stock))
        }
    } else if (Cagb_change_source == "Gridded_nc" | Cagb_change_source == "Gridded_tif") {
        # If there are any values in the analysis window
        if (Cagb_change_all$data_available) {        
            # Extract gross wood production / growth and uncertainty information
            # NOTE: assume default uncertainty (+/- scale)
            # Extract the current location from the gridded dataset
            output = extract_timeseries_observations_with_uncertainty(grid_long_loc,grid_lat_loc,timestep_days,years_to_load,doy_obs,
                                                                      Cagb_change_all,agg_func = "mean", na_flag = -9999,
                                                                      est_var_name_in = "Cagb_stock_gCm2day",
                                                                      unc_var_name_in = "Cagb_stock_uncertainty_gCm2day",
                                                                      lag_var_name_in = "Cagb_stock_lag",
                                                                      est_var_name_out = "Cagb_stock",
                                                                      unc_var_name_out = "Cagb_stock_unc",
                                                                      lag_var_name_out = "Cagb_stock_lag")                                                  
            Cagb_change = output$Cagb_change ; Cagb_change_unc = output$Cagb_change_unc ; Cagb_change_lag = output$Cagb_change_lag
            # Tidy up
            rm(output)
        } else {
            # assume no data available
            Cagb_change = -9999 ; Cagb_change_unc = -9999 ; Cagb_change_lag = -9999
        }
    } else {
        # assume no data available
        Cagb_stock = -9999 ; Cagb_stock_unc = -9999 ; Cagb_stock_lag = -9999
    }
    # apply lower bound in all cases to the uncertainty
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    Cagb_stock_unc[Cagb_stock_unc >= 0] = pmax(100,sqrt(Cagb_stock_unc[Cagb_stock_unc >= 0]**2 + (0.1*mean(Cagb_stock[Cagb_stock >= 0]))**2))

    ###
    ## Get some AGB net change information (gC/m2/day; time series)

    if (Cagb_change_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Cagb_change = read_site_specific_obs("Cagb_change_gCm2day",infile)
        Cagb_change_unc = read_site_specific_obs("Cagb_change_uncertainty_gCm2day",infile)
        Cagb_change_lag = read_site_specific_obs("Cagb_change_lag_step",infile) # in model time steps
        # Has uncertainty information been provided?
        if (length(Cagb_change_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cagb_change_unc = rep(-9999,times = length(Cagb_change))
            Cagb_change_unc[which(Cagb_change > 0)] = 0.25 * Cagb_change[which(Cagb_change > 0)]
        }
        # Has lag information been provided
        if (length(Cagb_change_lag) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cagb_change_lag = rep(-9999,times = length(Cagb_change))
            Cagb_change_lag[which(Cagb_change > 0)] = 0 # assume applies to current time step only
        }
    } else if (Cagb_change_source == "Gridded_nc" | Cagb_change_source == "Gridded_tif") {
        # If there are any values in the analysis window
        if (Cagb_change_all$data_available) {        
            # Extract gross wood production / growth and uncertainty information
            # NOTE: assume default uncertainty (+/- scale)
            # Extract the current location from the gridded dataset
            output = extract_timeseries_observations_with_uncertainty(grid_long_loc,grid_lat_loc,timestep_days,years_to_load,doy_obs,
                                                                      Cagb_change_all,agg_func = "mean", na_flag = -9999,
                                                                      est_var_name_in = "Cagb_change_gCm2day",
                                                                      unc_var_name_in = "Cagb_change_uncertainty_gCm2day",
                                                                      lag_var_name_in = "Cagb_change_lag",
                                                                      est_var_name_out = "Cagb_change",
                                                                      unc_var_name_out = "Cagb_change_unc",
                                                                      lag_var_name_out = "Cagb_change_lag")                                                  
            Cagb_change = output$Cagb_change ; Cagb_change_unc = output$Cagb_change_unc ; Cagb_change_lag = output$Cagb_change_lag
            # Tidy up
            rm(output)
        } else {
            # assume no data available
            Cagb_change = -9999 ; Cagb_change_unc = -9999 ; Cagb_change_lag = -9999
        }
    } else {
        # assume no data available
        Cagb_change = -9999 ; Cagb_change_unc = -9999 ; Cagb_change_lag = -9999
    }
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    Cagb_change_unc[Cagb_change_unc >= 0] = pmax(0.1,sqrt(Cagb_change_unc[Cagb_change_unc >= 0]**2 + (0.1*mean(Cagb_change[Cagb_change_unc >= 0]))**2))

    ###
    ## Get some gross woody production information (gC/m2/day; time series)

    if (Cagb_growth_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Cagb_growth = read_site_specific_obs("Cagb_production_gCm2day",infile)
        Cagb_growth_unc = read_site_specific_obs("Cagb_production_uncertainty_gCm2day",infile)
        Cagb_growth_lag = read_site_specific_obs("Cagb_production_lag_step",infile) # in model time steps
        # Has uncertainty information been provided?
        if (length(Cagb_growth_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cagb_growth_unc = rep(-9999,times = length(Cagb_growth))
            Cagb_growth_unc[which(Cagb_growth > 0)] = 0.25 * Cagb_growth[which(Cagb_growth > 0)]
        }
        # Has lag information been provided
        if (length(Cagb_growth_lag) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cagb_growth_lag = rep(-9999,times = length(Cagb_growth))
            Cagb_growth_lag[which(Cagb_growth > 0)] = 0 # assume applies to current time step only
        }
    } else if (Cagb_growth_source == "Gridded_nc" | Cagb_growth_source == "Gridded_tif") {
        # If there are any values in the analysis window
        if (Cagb_growth_all$data_available) {      
            # Extract gross wood production / growth and uncertainty information
            # NOTE: assume default uncertainty (+/- scale)
            # Extract the current location from the gridded dataset
            output = extract_timeseries_observations_with_uncertainty(grid_long_loc,grid_lat_loc,timestep_days,years_to_load,doy_obs,
                                                                      Cagb_growth_all,agg_func = "mean", na_flag = -9999,
                                                                      est_var_name_in = "Cagb_growth_gCm2day",
                                                                      unc_var_name_in = "Cagb_growth_uncertainty_gCm2day",
                                                                      lag_var_name_in = "Cagb_growth_lag",
                                                                      est_var_name_out = "Cagb_growth",
                                                                      unc_var_name_out = "Cagb_growth_unc",
                                                                      lag_var_name_out = "Cagb_growth_lag")   
            # Assign to local variables                           
            Cagb_growth = output$Cagb_growth ; Cagb_growth_unc = output$Cagb_growth_unc ; Cagb_growth_lag = output$Cagb_growth_lag
            # Tidy up
            rm(output)
        } else {
            # assume no data available
            Cagb_growth = -9999 ; Cagb_growth_unc = -9999 ; Cagb_growth_lag = -9999
        }
    } else {
        # assume no data available
        Cagb_growth = -9999 ; Cagb_growth_unc = -9999 ; Cagb_growth_lag = -9999
    }
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    Cagb_growth_unc[Cagb_growth_unc >= 0] = pmax(0.1,sqrt(Cagb_growth_unc[Cagb_growth_unc >= 0]**2 + (0.1*mean(Cagb_growth[Cagb_growth_unc >= 0]))**2))

    ###
    ## Get some Wood loss information (gC/m2/day; time series)

    if (Cagb_mortality_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Cagb_loss = read_site_specific_obs("Cagb_loss_gCm2day",infile)
        Cagb_loss_unc = read_site_specific_obs("Cagb_loss_uncertainty_gCm2day",infile)
        Cagb_loss_lag = read_site_specific_obs("Cagb_loss_lag_step",infile) # in model time steps
        # Has uncertainty information been provided?
        if (length(Cagb_loss_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cagb_loss_unc = rep(-9999,times = length(Cagb_loss))
            Cagb_loss_unc[which(Cagb_loss > 0)] = 0.25 * Cagb_loss[which(Cagb_loss > 0)]
        }
        # Has lag information been provided
        if (length(Cagb_loss_lag) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cagb_loss_lag = rep(-9999,times = length(Cagb_loss))
            Cagb_loss_lag[which(Cagb_loss > 0)] = 0 # assume applies to current time step only
        }
    } else if (Cagb_mortality_source == "Gridded_nc" | Cagb_mortality_source == "Gridded_tif") {
        # If there are any values in the analysis window
        if (Cagb_loss_all$data_available) {      
            # Extract wood mortality and uncertainty information
            # NOTE: assume default uncertainty (+/- scale)
            # Extract the current location from the gridded dataset
            output = extract_timeseries_observations_with_uncertainty(grid_long_loc,grid_lat_loc,timestep_days,years_to_load,doy_obs,
                                                                      Cagb_loss_all,agg_func = "mean", na_flag = -9999,
                                                                      est_var_name_in = "Cagb_loss_gCm2day",
                                                                      unc_var_name_in = "Cagb_loss_uncertainty_gCm2day",
                                                                      lag_var_name_in = "Cagb_loss_lag",
                                                                      est_var_name_out = "Cagb_loss",
                                                                      unc_var_name_out = "Cagb_loss_unc",
                                                                      lag_var_name_out = "Cagb_loss_lag")   
            # Assign to local variables                           
            Cagb_loss = output$Cagb_loss ; Cagb_loss_unc = output$Cagb_loss_unc ; Cagb_loss_lag = output$Cagb_loss_lag

            # Tidy up
            rm(output)
        } else {
            # assume no data available
            Cagb_loss = -9999 ; Cagb_loss_unc = -9999 ; Cagb_loss_lag = -9999
        }
    } else {
        # assume no data available
        Cagb_loss = -9999 ; Cagb_loss_unc = -9999 ; Cagb_loss_lag = -9999
    }
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    Cagb_loss_unc[Cagb_loss_unc >= 0] = pmax(0.1,sqrt(Cagb_loss_unc[Cagb_loss_unc >= 0]**2 + 
                                                                  (0.1*mean(Cagb_loss[Cagb_loss_unc >= 0]))**2))



    ###
    ## Get some Croots information (stock)

    if (Croots_stock_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Croots_stock = read_site_specific_obs("Croots_stock_gCm2",infile)
        Croots_stock_unc = read_site_specific_obs("Croots_stock_unc_gCm2",infile)
        Croots_stock_lag = read_site_specific_obs("Croots_stock_lag_step",infile)
        if (length(Croots_stock_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Croots_stock_unc = rep(-9999,times = length(Croots_stock))
            Croots_stock_unc[which(Croots_stock != -9999)] = abs(0.44 * Croots_stock[which(Croots_stock != -9999)])
        }
        if (length(Croots_stock_lag) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
            Croots_stock_lag = rep(1,times = length(Croots_stock))
        }
    } else {
        # assume no data available
        Croots_stock = -9999 ; Croots_stock_unc = -9999 ; Croots_stock_lag = -9999
    }

    ###
    ## Get some Clitter information (stock)

    if (Clit_stock_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Clit_stock = read_site_specific_obs("Clit_stock_gCm2",infile)
        Clit_stock_unc = read_site_specific_obs("Clit_stock_unc_gCm2",infile)
        Clit_stock_lag = read_site_specific_obs("Clit_stock_lag_step",infile)
        if (length(Clit_stock_unc) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
           Clit_stock_unc = rep(-9999,times = length(Clit_stock))
           Clit_stock_unc[which(Clit_stock != -9999)] = abs(0.38 * Clit_stock[which(Clit_stock != -9999)])
        }
        if (length(Clit_stock_lag) == 1) {
            # on the other hand if not then we have no uncertainty info, so use default
           Clit_stock_lag = rep(1,times = length(Clit_stock))
        }
    } else {
        # assume no data available
        Clit_stock = -9999 ; Clit_stock_unc = -9999 ; Clit_stock_lag = -9999
    }

    ###
    ## Get some Csom information (stock)

    if (Csom_stock_source == "site_specific") {
      infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
      Csom_stock = read_site_specific_obs("Csom_stock_gCm2",infile)
      Csom_stock_unc = read_site_specific_obs("Csom_stock_unc_gCm2",infile)
      Csom_stock_lag = read_site_specific_obs("Csom_stock_lag_step",infile)
      if (length(Csom_stock_unc) == 1) {
        # on the other hand if not then we have no uncertainty info, so use default
        Csom_stock_unc = rep(-9999,times = length(Csom_stock))
        Csom_stock_unc[which(Csom_stock != -9999)] = abs(0.24 * Csom_stock[which(Csom_stock != -9999)])
      }
      if (length(Csom_stock_lag) == 1) {
        # on the other hand if not then we have no uncertainty info, so use default
        Csom_stock_lag = rep(1,times = length(Csom_stock))
      }
    } else {
      # assume no data available
      Csom_stock = -9999 ; Csom_stock_unc = -9999 ; Csom_stock_lag = -9999
    }
    # Now assuming we have actual information we need to add the model structural uncertainty.
    # A structural of uncertainty has been estimates at ~ 1000 gC/m2 based on Smallman et al., (2017)
    #Csom_stock_unc[Csom_stock_unc >= 0] = sqrt(Csom_stock_unc[Csom_stock_unc >= 0]**2 + 1000**2)
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    Csom_stock_unc[Csom_stock_unc >= 0] = pmax(100,sqrt(Csom_stock_unc[Csom_stock_unc >= 0]**2 + (0.1*mean(Csom_stock[Csom_stock > 0]))**2))

    ###
    ## Get some Ccoarseroot information (stock)

    if (Ccoarseroot_stock_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        Ccoarseroot_stock = read_site_specific_obs("Ccoarseroot_stock_gCm2",infile)
        Ccoarseroot_stock_unc = read_site_specific_obs("Ccoarseroot_stock_unc_gCm2",infile)
        Ccoarseroot_stock_lag = read_site_specific_obs("Ccoarseroot_stock_lag_step",infile)
      if (length(Ccoarseroot_stock_unc) == 1) {
          # on the other hand if not then we have no uncertainty info, so use default
          Ccoarseroot_stock_unc = rep(-9999,times = length(Ccoarseroot_stock))
          Ccoarseroot_stock_unc[which(Ccoarseroot_stock != -9999)] = abs(0.24 * Ccoarseroot_stock[which(Ccoarseroot_stock != -9999)])
      }
      if (length(Ccoarseroot_stock_lag) == 1) {
          # on the other hand if not then we have no uncertainty info, so use default
          Ccoarseroot_stock_lag = rep(1,times = length(Ccoarseroot_stock))
      }
    } else {
        # assume no data available
        Ccoarseroot_stock = -9999 ; Ccoarseroot_stock_unc = -9999 ; Ccoarseroot_stock_lag = -9999
    }

    ###
    ## Get some snow water equivalent (kgH2O/m2 or mm)

    if (snow_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        SWE = read_site_specific_obs("snow_water_kgH2Om2",infile)
        SWE_unc = read_site_specific_obs("snow_water_unc_kgH2Om2",infile)
        SWE_lag = read_site_specific_obs("snow_water_lag_step",infile)
        if (length(which(SWE_unc != -9999)) == 0) {
            # on the other hand if not then we have no uncertainty info, so use default
           SWE_unc=rep(sd(SWE[which(SWE != -9999)],na.rm=TRUE),length.out=length(SWE))
        }
        if (length(which(SWE_lag != -9999)) == 0) {
            # on the other hand if not then we have no uncertainty info, so use default
           SWE_lag=rep(1,length.out=length(SWE))
        }
    } else {
        # assume no data available
        SWE = -9999 ; SWE_unc = -9999 ; SWE_lag = -9999
    }

    ###
    ## Extract the local information for static information with uncertainty,
    ## i.e. those values which are assimilated
    ###

    ###
    ## Get some initial Csom (gC/m2) information

    if (Csom_source == "Gridded_nc" | Csom_source == "Gridded_tif") {
        # Get an initial estimate for soil C
        output = extract_static_observations_with_uncertainty(grid_long_loc,grid_lat_loc,Csom_all,
                                                              na_flag = -9999,
                                                              est_var_name_in="Csom", 
                                                              unc_var_name_in="Csom_unc",
                                                              est_var_name_out="Csom_initial",
                                                              unc_var_name_out="Csom_initial_unc") 
        Csom_initial = output$Csom_initial ; Csom_initial_unc = output$Csom_initial_unc
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
    if (Csom_initial_unc > 0) { Csom_initial_unc = max(100,sqrt(Csom_initial_unc**2 + (0.1*Csom_initial)**2)) } 

    ###
    ## Get some crop management information (day)

    if (crop_management_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        planting_doy = read_site_specific_obs("planting_doy_initial",infile)
        planting_doy_unc = read_site_specific_obs("planting_doy_unc_initial",infile)
        growing_season_doy = read_site_specific_obs("growing_season_doy_initial",infile)
        growing_season_doy_unc = read_site_specific_obs("growing_season_doy_unc_initial",infile)
        # Prior parameter range for sowing date span 365.25-> but rescale to 1-365.25 by taking the modulus.
        # This means that we must put this prior into the parameter prior range space
        planting_doy = planting_doy + 365.25
    } else {
        # assume no data available
        planting_doy = -9999 ; planting_doy_unc = -9999 # days
        growing_season_doy = -9999  ; growing_season_doy_unc = -9999 # days # note +365.25 to account for the parameter range
        #planting_doy = 273 + 365.25 ; planting_doy_unc = 15 # days
        #growing_season_doy = 330 ; growing_season_doy_unc = 15 # days

    }

    ###
    ## Get some Cfoliage information (initial conditions)

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

    if (Cwood_initial_source == "site_specific") {
        infile=paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        Cwood_initial=read_site_specific_obs("Cwood_initial_gCm2",infile)
        Cwood_initial_unc=read_site_specific_obs("Cwood_initial_unc_gCm2",infile)
        if (Cwood_initial_unc == -9999 & Cwood_initial > 0) {
            # on the other hand if not then we have no uncertainty info, so use default
            Cwood_initial_unc = 0.25 * Cwood_initial
        }
    } else if (Cwood_initial_source == "Gridded_nc" | Cwood_initial_source == "Gridded_tif") {
        # Extract the initial wood stock estimate
        output = extract_static_observations_with_uncertainty(grid_long_loc,grid_lat_loc,Cwood_initial_all,
                                                              na_flag = -9999,
                                                              est_var_name_in="biomass_gCm2",
                                                              unc_var_name_in="biomass_uncertainty_gCm2",
                                                              est_var_name_out="Cwood_stock",
                                                              unc_var_name_out="Cwood_stock_unc") 
        Cwood_initial = output$Cwood_stock ; Cwood_initial_unc = output$Cwood_stock_unc
    } else {
        # assume no data available
        Cwood_initial=-9999 ; Cwood_initial_unc=-9999
    }
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    Cwood_initial_unc[Cwood_initial_unc >= 0] = pmax(100,sqrt(Cwood_initial_unc[Cwood_initial_unc >= 0]**2 + 
                                                              (0.1*mean(Cwood_initial[Cwood_initial_unc >= 0]))**2),na.rm=TRUE)

    ###
    ## Get some Croots information (initial conditions)
    
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
    ## Get some leaf carbon per unit leaf area (gC/m2) information 

    if (lca_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        lca = read_site_specific_obs("LCA_gCm2",infile)
        lca_unc = read_site_specific_obs("LCA_unc_gCm2",infile)
    } else if (lca_source == "Gridded_nc" | lca_source == "Gridded_tif") {
        # get leaf carbon per unit leaf area from gridded dataset  
        output = extract_static_observations_with_uncertainty(grid_long_loc,grid_lat_loc,lca_all,
                                                              na_flag = -9999,
                                                              est_var_name_in="lca_gCm2",
                                                              unc_var_name_in="lca_uncertainty_gCm2",
                                                              est_var_name_out="lca_gCm2",
                                                              unc_var_name_out="lca_unc_gCm2") 
        # Load into local variables
        lca = output$lca_gCm2
        lca_unc = output$lca_unc_gCm2
    } else {
        # assume no data available
        lca = -9999 ; lca_unc = -9999
    }
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    if (lca > -9999) { lca_unc = max(10,sqrt(lca_unc**2 + (0.1*lca)**2)) }

    ###
    ## Get some leaf lifespan information (y)

    if (leaflifespan_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        leaflifespan = read_site_specific_obs("leaflifespan_y",infile)
        leaflifespan_unc = read_site_specific_obs("leaflifespan_unc_y",infile)
    } else if (lca_source == "Gridded_nc" | lca_source == "Gridded_tif") {
        # get leaf carbon per unit leaf area from gridded dataset  
        output = extract_static_observations_with_uncertainty(grid_long_loc,grid_lat_loc,lca_all,
                                                              na_flag = -9999,
                                                              est_var_name_in="leaflifespan_y",
                                                              unc_var_name_in="leaflifespan_uncertainty_y",
                                                              est_var_name_out="leaflifespan_y",
                                                              unc_var_name_out="leaflifespan_unc_y") 
        # Load into local variables
        leaflifespan = output$leaflifespan_y
        leaflifespan_unc = output$leaflifespan_unc_y
    } else {
        # assume no data available
        leaflifespan = -9999 ; leaflifespan_unc = -9999
    }
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    if (leaflifespan > -9999) { leaflifespan_unc = max(0.1,sqrt(leaflifespan_unc**2 + (0.1*leaflifespan)**2)) }

    ###
    ## Get some labile release timing information (d)

    if (labile_release_timing_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        labile_release_timing = read_site_specific_obs("labile_release_timing_d",infile)
        labile_release_timing_unc = read_site_specific_obs("labile_release_timing_unc_d",infile)
    } else if (lca_source == "Gridded_nc" | lca_source == "Gridded_tif") {
        # get leaf carbon per unit leaf area from gridded dataset  
        output = extract_static_observations_with_uncertainty(grid_long_loc,grid_lat_loc,lca_all,
                                                              na_flag = -9999,
                                                              est_var_name_in="labile_release_timing_d",
                                                              unc_var_name_in="labile_release_timing_uncertainty_d",
                                                              est_var_name_out="labile_release_timing_d",
                                                              unc_var_name_out="labile_release_timing_unc_d") 
        # Load into local variables
        labile_release_timing = output$labile_release_timing_d
        labile_release_timing_unc = output$labile_release_timing_unc_d
    } else {
        # assume no data available
        labile_release_timing = -9999 ; labile_release_timing_unc = -9999
    }
    # Assumed uncertainty structure
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    if (labile_release_timing > -9999) { labile_release_timing_unc = max(7,labile_release_timing_unc**2) }


    ###
    ## Get some labile release period information (d)

    if (labile_release_period_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        labile_release_period = read_site_specific_obs("labile_release_period_d",infile)
        labile_release_period_unc = read_site_specific_obs("labile_release_period_unc_d",infile)
    } else if (lca_source == "Gridded_nc" | lca_source == "Gridded_tif") {
        # get leaf carbon per unit leaf area from gridded dataset  
        output = extract_static_observations_with_uncertainty(grid_long_loc,grid_lat_loc,lca_all,
                                                              na_flag = -9999,
                                                              est_var_name_in="labile_release_period_d",
                                                              unc_var_name_in="labile_release_period_uncertainty_d",
                                                              est_var_name_out="labile_release_period_d",
                                                              unc_var_name_out="labile_release_period_unc_d") 
        # Load into local variables
        labile_release_period = output$labile_release_period_d
        labile_release_period_unc = output$labile_release_period_unc_d
    } else {
        # assume no data available
        labile_release_period = -9999 ; labile_release_period_unc = -9999
    }
    # Assumed uncertainty structure
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    if (labile_release_period > -9999) { labile_release_period_unc = max(7,labile_release_period_unc**2) }



    ###
    ## Get some maximum rooting depth information (m)

    if (MaxRootDepth_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        MaxRootDepth = read_site_specific_obs("MaxRootDepth_m",infile)
        MaxRootDepth_unc = read_site_specific_obs("MaxRootDepth_unc_m",infile)
    } else if (MaxRootDepth_source == "Gridded_nc" | MaxRootDepth_source == "Gridded_tif") {
        # get leaf carbon per unit leaf area from gridded dataset  
        output = extract_static_observations_with_uncertainty(grid_long_loc,grid_lat_loc,MaxRootDepth_all,
                                                              na_flag = -9999,
                                                              est_var_name_in="MaxRootDepth_m",
                                                              unc_var_name_in="MaxRootDepth_uncertainty_m",
                                                              est_var_name_out="MaxRootDepth_m",
                                                              unc_var_name_out="MaxRootDepth_unc_m") 
        # Load into local variables
        MaxRootDepth = output$MaxRootDepth_m
        MaxRootDepth_unc = output$MaxRootDepth_unc_m
    } else {
        # assume no data available
        MaxRootDepth = -9999 ; MaxRootDepth_unc = -9999
    }
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    if (MaxRootDepth > -9999) { MaxRootDepth_unc = max(0.3,sqrt(MaxRootDepth_unc**2 + (0.1*MaxRootDepth)**2)) }

    ###
    ## Get some prior info on fraction of Cwood belowground as course roots (fraction) 

    if (frac_Cwood_coarse_root_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        frac_Cwood_coarse_root_prior = read_site_specific_obs("frac_Cwood_coarse_root_prior",infile)
        frac_Cwood_coarse_root_prior_unc = read_site_specific_obs("frac_Cwood_coarse_root_prior_unc",infile)
    } else if (MaxRootDepth_source == "Gridded_nc" | frac_Cwood_coarse_root_source == "Gridded_tif") {
        # get leaf carbon per unit leaf area from gridded dataset  
        output = extract_static_observations_with_uncertainty(grid_long_loc,grid_lat_loc,frac_Cwood_coarse_root_all,
                                                              na_flag = -9999,
                                                              est_var_name_in="frac_Cwood_coarse_root_m",
                                                              unc_var_name_in="frac_Cwood_coarse_root_uncertainty_m",
                                                              est_var_name_out="frac_Cwood_coarse_root_m",
                                                              unc_var_name_out="frac_Cwood_coarse_root_unc_m") 
        # Load into local variables
        frac_Cwood_coarse_root = output$frac_Cwood_coarse_root_m
        frac_Cwood_coarse_root_unc = output$frac_Cwood_coarse_root_unc_m
    }  else {
        # assume no data available
        frac_Cwood_coarse_root_prior = -9999 ; frac_Cwood_coarse_root_prior_unc = -9999
        # If we have some wood stock information we can do better with a prior value
        if (max(Cwood_stock) > 0) {
            # Based on the creation of a log~log fit to the below ground stock estimates 
            # using allometry from # Saatchi et al., (2011), PNAS, 108, 9899-9904, https://www.pnas.org/content/108/24/9899
            # BGB/AGB = exp(log(TWB) * -9.356e-2 + -6.536e-1 ), R2 = 0.9997, +0.02 is the max error between fits
            frac_Cwood_coarse_root_prior = exp(log(max(Cwood_stock)) * -9.356e-2 - 6.536e-1)
            frac_Cwood_coarse_root_prior_unc = exp(log(max(Cwood_stock_unc)) * -9.356e-2 - 6.536e-1) + 0.02 
        }       
    }
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    if (frac_Cwood_coarse_root_prior > -9999) { 
        frac_Cwood_coarse_root_prior_unc = max(0.1,sqrt(frac_Cwood_coarse_root_prior_unc**2 + (0.1*frac_Cwood_coarse_root_prior)**2)) 
    }

    ###
    ## Get soil water fraction (m3/m3 vol.)

    if (soilwater_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_timeseries_obs.csv",sep="")
        soilwater = read_site_specific_obs("soil_water_m3m3",infile)
        soilwater_unc = read_site_specific_obs("soil_water_unc_m3m3",infile)
        soilwater_lag = read_site_specific_obs("soil_water_lag_step",infile)
        if (soilwater_unc == -9999 & soilwater > 0) {
          # on the other hand if not then we have no uncertainty info, so use default
          soilwater_unc = 0.10 * soilwater
        }
        if (soilwater_lag == -9999 & soilwater > 0) {
          # on the other hand if not then we have no uncertainty info, so use default
          soilwater_lag = rep(1, length(soilwater))
        }        
    } else {
        # assume no data available
        soilwater = -9999 ; soilwater_unc = -9999 ; soilwater_lag = -9999
    }

    ###
    ## Get some Cwood information (potential stock)

    if (Cwood_potential_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        Cwood_potential = read_site_specific_obs("Cwood_potential_gCm2",infile)
        Cwood_potential_unc = read_site_specific_obs("Cwood_potential_unc_gCm2",infile)
    } else if (Cwood_potential_source == "Gridded_nc" | Cwood_potential_source == "Gridded_tif") {
        # Extract potential wood stock information
        output = extract_static_observations_with_uncertainty(grid_long_loc,grid_lat_loc,Cwood_potential_all,
                                                              na_flag = -9999,
                                                              est_var_name_in="biomass_gCm2",
                                                              unc_var_name_in="biomass_uncertainty_gCm2",
                                                              est_var_name_out="Cwood_stock",
                                                              unc_var_name_out="Cwood_stock_unc")         
        Cwood_potential = output$Cwood_stock
        Cwood_potential_unc = output$Cwood_stock_unc
    } else {
        # assume no data available
        Cwood_potential = -9999 ; Cwood_potential_unc = -9999
    }

    ###
    ## Get minimum LWP (MPa) information 

    if (minLWP_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        minLWP = read_site_specific_obs("minLWP_MPa",infile)
        minLWP_unc = read_site_specific_obs("minLWP_unc_MPa",infile)
    } else {
        # assume no data available
        minLWP = -9999 ; minLWP_unc = -9999
    }

    ###
    ## Get some information on the temperature sensitivity of heterotrophic respiration.
    ## The contents of the file has the same name, but may be expressed as a Q10 or 
    ## exp(c*T) depending on the DALEC model version. Please be very careful that these match

    if (RhetQ10_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        RhetQ10 = read_site_specific_obs("RhetQ10",infile)
        RhetQ10_unc = read_site_specific_obs("RhetQ10_unc",infile)
    } else if (RhetQ10_source == "Gridded_nc" | RhetQ10_source == "Gridded_tif") {
        # get Rhet temperature sensitivity information
        output = extract_static_observations_with_uncertainty(grid_long_loc,grid_lat_loc,RhetQ10_all,
                                                              na_flag = -9999,
                                                              est_var_name_in="RhetQ10",
                                                              unc_var_name_in="RhetQ10_uncertainty",
                                                              est_var_name_out="RhetQ10",
                                                              unc_var_name_out="RhetQ10_unc") 
        # Load into local variables
        RhetQ10 = output$RhetQ10
        RhetQ10_unc = output$RhetQ10_unc
    } else {
        # assume no data available
        RhetQ10 = -9999 ; RhetQ10_unc = -9999
    }
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    if (RhetQ10 > -9999) { RhetQ10_unc = max(0.005,sqrt(RhetQ10_unc**2 + (0.1*RhetQ10)**2)) }

    ###
    ## Mean transit time (MTT, years) for the soil organic matter (aka soil, som).
    ## This will be applied as an other prior, not directly on the turnover parameter
    ## as these are not equivalent

    if (MTTsom_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
        MTTsom = read_site_specific_obs("MTTsom",infile)
        MTTsom_unc = read_site_specific_obs("MTTsom_unc",infile)
    } else if (MTTsom_source == "Gridded_nc" | MTTsom_source == "Gridded_tif") {
        # get mean transis time of soil information
        output = extract_static_observations_with_uncertainty(grid_long_loc,grid_lat_loc,MTTsom_all,
                                                              na_flag = -9999,
                                                              est_var_name_in="MTTsom",
                                                              unc_var_name_in="MTTsom_uncertainty",
                                                              est_var_name_out="MTTsom",
                                                              unc_var_name_out="MTTsom_unc") 
        # Load into local variables
        MTTsom = output$MTTsom
        MTTsom_unc = output$MTTsom_unc
    } else {
        # assume no data available
        MTTsom = -9999 ; MTTsom_unc = -9999
    }
    # Assumed uncertainty structure as agreed with Anthony Bloom
    # NOTE minimum uncertainty bound irrespective of the dataset estimates
    if (MTTsom > -9999) { MTTsom_unc = max(10,sqrt(MTTsom_unc**2 + (0.1*MTTsom)**2)) }

    ###
    ## Extract the local information for timeseries information without observations,
    ## i.e. those values which are forcings
    ###

    ###
    ## Get some deforestation information (fraction time series)

    if (deforestation_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_timeseries_met.csv",sep="")
        deforestation = read_site_specific_obs("deforestation_fraction",infile)
        # If the default variable can't be found, check for the lai_loss variable
        # used in the managed grassland system
        if (length(deforestation) == 1 && deforestation == -9999) {
            deforestation = read_site_specific_obs("lai_loss",infile)
        }
        forest_management = read_site_specific_obs("management_type",infile)
        if (length(forest_management) == 1) {forest_management = rep(2, times = length(deforestation))}
        yield_class = -9999 #read_site_specific_obs("yield_class",infile)
        age = read_site_specific_obs("age",infile)
        if (length(age) > 1) {age = age[1]} # we only want the age at the beginning of the simulation
    } else if (deforestation_source == "Gridded_nc" | deforestation_source == "Gridded_tif") {
        # Extract from the gridded array
#        output = extract_timeseries_observations_without_uncertainty(grid_long_loc,grid_lat_loc,timestep_days,years_to_load,doy_obs,
#                                                                     forest_all,agg_func = "sum", na_flag = 0,
#                                                                     est_var_name_in="loss_fraction",lag_var_name_in="loss_fraction_lag",
#                                                                     est_var_name_out="deforestation",lag_var_name_out="deforestation_lag")      
#        deforestation = output$deforestation ; deforestation_lag = output$deforestation_lag
        output = extract_timeseries_forcing(grid_long_loc,grid_lat_loc,timestep_days,years_to_load,doy_obs,
                                            forest_all,agg_func = "sum", fraction = TRUE,na_flag = 0,
                                            est_var_name_in="loss_fraction",lag_var_name_in="loss_fraction_lag",
                                            est_var_name_out="deforestation")
        deforestation = output$deforestation 
        yield_class = -9999
        age = -9999
        forest_management = 2 # Default option, check model specific code for their actual effects
    } else {
        # assume no data available
        deforestation = 0
        forest_management = 2
        yield_class = 0
        age = -9999
    }

    ###
    ## Get some burnt area information (fraction time series)

    if (burnt_area_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_timeseries_met.csv",sep="")
        burnt_area = read_site_specific_obs("burnt_area_fraction",infile)
    } else if (burnt_area_source == " "){
        # assume no data available
        burnt_area = 0 
    } else {
        # Extract from the gridded array
#        output = extract_timeseries_observations_without_uncertainty(grid_long_loc,grid_lat_loc,timestep_days,years_to_load,doy_obs,
#                                                                     burnt_all,agg_func = "sum",na_flag = 0,
#                                                                     est_var_name_in="burnt_area",lag_var_name_in="burnt_area_lag",
#                                                                     est_var_name_out="burnt_area",lag_var_name_out="burnt_area_lag")
        output = extract_timeseries_forcing(grid_long_loc,grid_lat_loc,timestep_days,years_to_load,doy_obs,
                                            burnt_all,agg_func = "sum",fraction = TRUE, na_flag = 0,
                                            est_var_name_in="burnt_area",lag_var_name_in="burnt_area_lag",
                                            est_var_name_out="burnt_area")
        # Extract out of the output object
        burnt_area = output$burnt_area 
    }

    ###
    ## Get some LAI change forcing information (m2/m2 time series)
    ## Currently used for managed grasslands to enforce periods of loss due to cutting or grazing.
    ## NOTE: the time series here will have both LAI losses and gains provided, the bounding will
    ##       occur specifically in binary_data.r for the grassland models only (DALEC.16/17)

    if (lai_change_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_timeseries_met.csv",sep="")
        lai_change = read_site_specific_obs("lai_change_m2m2",infile) # NOTE "lai_loss" should be changed to "lai_change"
    } else if (lai_change_source == "Gridded_nc" | lai_change_source == "Gridded_tif") {
        # Extract from the gridded array
        output = extract_timeseries_forcing(grid_long_loc,grid_lat_loc,timestep_days,years_to_load,doy_obs,
                                            dlai_all,agg_func = "sum", fraction = FALSE, na_flag = 0,
                                            est_var_name_in="dlai_m2m2",lag_var_name_in="dlai_lag_day",
                                            est_var_name_out="lai_change")
        lai_change = output$lai_change 
    } else {
        # assume no data available
        lai_change = 0
    }

    ###
    ## Extract the local information for static information without observations,
    ## i.e. those values which are forcings
    ###

    ###
    ## Get some sand / clay information (%)

    if (sand_clay_source == "Gridded_nc" | sand_clay_source == "Gridded_tif") {
        ## Extract each layer and type in turn
        # Extract local sand content (top soil, 0-30cm)
        output = extract_static_observations_without_uncertainty(grid_long_loc,grid_lat_loc,sand_clay_all,
                                                                 na_flag = NA,
                                                                 est_var_name_in="top_sand",
                                                                 est_var_name_out="top_sand") 
        top_sand = output$top_sand 
        # Extract local sand content (bottom soil, 31-100cm)
        output = extract_static_observations_without_uncertainty(grid_long_loc,grid_lat_loc,sand_clay_all,
                                                                 na_flag = NA,
                                                                 est_var_name_in="bot_sand",
                                                                 est_var_name_out="bot_sand")           
        bot_sand = output$bot_sand
        # Extract local clay content (top soil, 0-30cm)
        output = extract_static_observations_without_uncertainty(grid_long_loc,grid_lat_loc,sand_clay_all,
                                                                 na_flag = NA,
                                                                 est_var_name_in="top_clay",
                                                                 est_var_name_out="top_clay")     
        top_clay = output$top_clay 
        # Extract local clay content (bottom soil, 31-100cm)
        output = extract_static_observations_without_uncertainty(grid_long_loc,grid_lat_loc,sand_clay_all,
                                                                 na_flag = NA,
                                                                 est_var_name_in="bot_clay",
                                                                 est_var_name_out="bot_clay")     
        bot_clay = output$bot_clay

        ## Sanity and mass balance checks
        # Guard against NaN values
        if (is.na(top_sand) | is.infinite(top_sand)) {top_sand = 40}
        if (is.na(bot_sand) | is.infinite(bot_sand)) {bot_sand = 40}
        if (is.na(top_clay) | is.infinite(top_clay)) {top_clay = 15}
        if (is.na(bot_clay) | is.infinite(bot_clay)) {bot_clay = 15}
        # ML based approaches, such as those typically available, generate the different layers 
        # independently. As a result correlations and mass balance will not be preserved. 
        # i.e. the sand / clay combinations can be > 100 %
        # 94 % chosesn as this is the highest total % found in the HWSD dataset
        if ((top_sand+top_clay) > 94) {
             tmp1 = top_sand / (top_sand + top_clay + 6) # 6 % is implicit in the 94 % max value for silt / gravel
             tmp2 = top_clay / (top_sand + top_clay + 6) # 6 % is implicit in the 94 % max value for silt / gravel
             top_sand = tmp1*100 ; top_clay = tmp2*100
        }
        if ((bot_sand+bot_clay) > 94) {
             tmp1 = bot_sand / (bot_sand + bot_clay + 6) # 6 % is implicit in the 94 % max value for silt / gravel
             tmp2 = bot_clay / (bot_sand + bot_clay + 6) # 6 % is implicit in the 94 % max value for silt / gravel
             bot_sand = tmp1*100 ; bot_clay = tmp2*100
        }

    } else if (sand_clay_source == "site_specific") {
        infile = paste(path_to_site_obs,site_name,"_initial_obs.csv",sep="")
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
    ## Prepare the final output object and return
    ###

    # return output now
    return(list(LAT = latlon_wanted[1], ctessel_pft = ctessel_pft, 
                top_sand = top_sand, bot_sand = bot_sand, 
                top_clay = top_clay, bot_clay = bot_clay, 
                LAI = lai, LAI_unc = lai_unc, LAI_lag = lai_lag, 
                GPP = GPP, GPP_unc = GPP_unc, GPP_lag = GPP_lag, 
                Fire = Fire, Fire_unc = Fire_unc, Fire_lag = Fire_lag,
                ET = ET, ET_unc = ET_unc, ET_lag = ET_lag, 
                NEE = NEE, NEE_unc = NEE_unc, NEE_lag = NEE_lag, 
                Reco = Reco, Reco_unc = Reco_unc, Reco_lag = Reco_lag,
                fAPAR = fapar, fAPAR_unc = fapar_unc, fAPAR_lag = fapar_lag,
                Cfol_stock = Cfol_stock, Cfol_stock_unc = Cfol_stock_unc, Cfol_stock_lag = Cfol_stock_lag,
                Cwood_stock = Cwood_stock, Cwood_stock_unc = Cwood_stock_unc, Cwood_stock_lag = Cwood_stock_lag, 
                Cagb_stock = Cagb_stock, Cagb_stock_unc = Cagb_stock_unc, Cagb_stock_lag = Cagb_stock_lag,
                Croots_stock = Croots_stock, Croots_stock_unc = Croots_stock_unc, Croots_stock_lag = Croots_stock_lag, 
                Clit_stock = Clit_stock, Clit_stock_unc = Clit_stock_unc, Clit_stock_lag = Clit_stock_lag,
                Csom_stock = Csom_stock, Csom_stock_unc = Csom_stock_unc, Csom_stock_lag = Csom_stock_lag, 
                Ccoarseroot_stock = Ccoarseroot_stock, Ccoarseroot_stock_unc = Ccoarseroot_stock_unc, Ccoarseroot_stock_lag = Ccoarseroot_stock_lag,
                nbe = nbe, nbe_unc = nbe_unc, nbe_lag = nbe_lag,
                SWE = SWE, SWE_unc = SWE_unc, SWE_lag = SWE_lag,
                soilwater = soilwater, soilwater_unc = soilwater_unc, soilwater_lag = soilwater_lag,
                Cwood_change = Cwood_change, Cwood_change_unc = Cwood_change_unc, Cwood_change_lag = Cwood_change_lag,
                Cwood_growth = Cwood_growth, Cwood_growth_unc = Cwood_growth_unc, Cwood_growth_lag = Cwood_growth_lag,
                Cwood_loss = Cwood_loss, Cwood_loss_unc = Cwood_loss_unc, Cwood_loss_lag = Cwood_loss_lag,
                harvest = harvest, harvest_unc = harvest_unc, harvest_lag = harvest_lag,
                Csom_initial = Csom_initial, Csom_initial_unc = Csom_initial_unc, 
                Cfol_initial = Cfol_initial, Cfol_initial_unc = Cfol_initial_unc,
                Cwood_initial = Cwood_initial, Cwood_initial_unc = Cwood_initial_unc, 
                Croots_initial = Croots_initial, Croots_initial_unc = Croots_initial_unc, 
                Clit_initial = Clit_initial, Clit_initial_unc = Clit_initial_unc, 
                deforestation = deforestation, 
                age = age, forest_management = forest_management, yield_class = yield_class,
                burnt_area = burnt_area, lai_change = lai_change,
                planting_doy = planting_doy, planting_doy_unc = planting_doy, 
                growing_season_doy = growing_season_doy, growing_season_doy_unc = growing_season_doy_unc,              
                Cwood_potential = Cwood_potential, Cwood_potential_unc = Cwood_potential_unc, 
                lca = lca, lca_unc = lca_unc,
                foliage_to_litter = foliage_to_litter, foliage_to_litter_unc = foliage_to_litter_unc, foliage_to_litter_lag = foliage_to_litter_lag,
                frac_Cwood_coarse_root_prior = frac_Cwood_coarse_root_prior, frac_Cwood_coarse_root_prior_unc = frac_Cwood_coarse_root_prior_unc,
                minLWP = minLWP, minLWP_unc = minLWP_unc, RhetQ10 = RhetQ10, RhetQ10_unc = RhetQ10_unc, MTTsom = MTTsom, MTTsom_unc = MTTsom_unc,
                MaxRootDepth = MaxRootDepth, MaxRootDepth_unc = MaxRootDepth_unc,
                Cagb_change = Cagb_change, Cagb_change_unc = Cagb_change_unc, Cagb_change_lag = Cagb_change_lag,
                Cagb_growth = Cagb_growth, Cagb_growth_unc = Cagb_growth_unc, Cagb_growth_lag = Cagb_growth_lag,
                Cagb_loss = Cagb_loss, Cagb_loss_unc = Cagb_loss_unc, Cagb_loss_lag = Cagb_loss_lag,
                leaflifespan = leaflifespan,
                labile_release_timing = labile_release_timing, labile_release_period = labile_release_period,
                frac_Cwood_coarse_root = frac_Cwood_coarse_root))

} # end function extract_obs

## Use byte compile
extract_obs<-cmpfun(extract_obs)
