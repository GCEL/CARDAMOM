
###
## Standard functions to extract obs data from:
## 1) time varying 
## 2) Static static datasets
## 3) With uncertainty or 
## 4) Without uncertainty
## 5) By averaging or
## 6) By sum
## 7) over varied lag periods
### 

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
# Four function used to extract location specific information on from gridded dataset, 
# not including meteorological forcings.
# Functions are:
# extract_timeseries_observations_with_uncertainty
# extract_timeseries_observations_without_uncertainty
# extract_static_observations_with_uncertainty
# extract_static_observations_without_uncertainty
# extract_timeseries_forcing
# 
# Author: T. Luke Smallman (24/01/2025)
#
#########################################################################################

extract_timeseries_observations_with_uncertainty<- function(i1,j1,timestep_days,years_to_load,doy_obs,
                                                            data_all,agg_func, na_flag,
                                                            est_var_name_in,unc_var_name_in,
                                                            lag_var_name_in,est_var_name_out,
                                                            unc_var_name_out,lag_var_name_out) {
                                                  
   # Update the user
   if (use_parallel == FALSE) {print(paste("extracting ",est_var_name_in," for current location ",Sys.time(),sep=""))}

   # Extract current location to local variable
   obs = data_all[[est_var_name_in]][i1,j1,]
   unc = data_all[[unc_var_name_in]][i1,j1,]
   lag = data_all[[lag_var_name_in]][i1,j1,]

   # declare output variable
   obs_out = array(NA, dim=length(doy_obs))
   obs_unc_out = array(NA, dim=length(doy_obs))
   obs_lag_out = array(NA, dim=length(doy_obs))

   ## Line up the days of year which have observations into a complete timeseries of days...
   b = 1 ; i = 1 ; a = 1 ; start_year = as.numeric(years_to_load[1])
   #print("...begin inserting LAI observations into model time steps")
   while (b <= length(data_all$doy_obs)) {

      # if we are in a year which is missing then we do not allow consideration of DOY
      if (start_year != data_all$missing_years[a]) {
          if (doy_obs[i] == data_all$doy_obs[b]) {
              obs_out[i] = obs[b] ; obs_unc_out[i] = unc[b] ; obs_lag_out[i] = lag[b] ; b = b + 1
          } # end if doy matches
      } # end if missing year

      # but we do keep counting through the total vector length which we expect
      i = i + 1

      # each time we come back to doy_obs[i]==1 we need to count on the year
      #if (length(data_all$doy_obs) < 10) {print(data_all$doy_obs)}

      if (doy_obs[i] == 1 & b <= length(data_all$doy_obs)) {
          # and if we have just been in a missing year we need to count on the missing years vector to
          if (start_year == data_all$missing_years[a]) { a = min(length(data_all$missing_years),a+1) }
          start_year = start_year + 1
      } # end if doy_obs[i] == 1

   } # end while condition

   ## Aggregate to the model time steps by sum or mean
   if (length(timestep_days) == 1 & timestep_days[1] == 1) {

       # well actually we do nothing

   } else {

       # Aggregate when the time step is not daily

       # Sanity check in case we have not been 
       # given a complete timeseries of step sizes
       if (length(timestep_days) == 1) {
           run_day_selector = seq(1,length(obs_out),timestep_days)
           timestep_days = rep(timestep_days, length.out=length(obs_out))
       }

       # Determine the actual cumulative number of days to have passed
       run_day_selector = cumsum(timestep_days)
       # create needed variables
       obs_agg = array(NA,dim=length(run_day_selector))
       obs_unc_agg = array(NA,dim=length(run_day_selector))
       obs_lag_agg = array(NA,dim=length(run_day_selector))

       # Determine whether we are aggregating by mean or by sum
       if (agg_func == "mean") {
           # Loop through timeseries and aggregate
           for (y in seq(1,length(run_day_selector))) {
                # Only run for time steps when there is a value
                pick = (run_day_selector[y]-timestep_days[y]+1):run_day_selector[y]
                if (length(which(is.na(obs_out[pick]) == FALSE)) > 0) {       
                    obs_agg[y] = weighted.mean(x = obs_out[pick], w = pmax(1,obs_lag_out[pick], na.rm=TRUE), na.rm=TRUE)
                    obs_unc_agg[y] = weighted.mean(x = obs_unc_out[pick], w = pmax(1,obs_lag_out[pick], na.rm=TRUE), na.rm=TRUE)
                    obs_lag_agg[y] = sum(obs_lag_out[pick], na.rm=TRUE)
                }
                   
           }
       } else if (agg_func == "sum") {
           # Loop through timeseries and aggregate
           for (y in seq(1,length(run_day_selector))) {
                # Only run for time steps when there is a value           
                pick = (run_day_selector[y]-timestep_days[y]+1):run_day_selector[y]
                if (length(which(is.na(obs_out[pick]) == FALSE)) > 0) {                           
                    # Take the mean first to allow for weighting based on the lags
                    obs_agg[y] = weighted.mean(x = obs_out[pick], w = obs_lag_out[pick]+1, na.rm=TRUE)
                    obs_unc_agg[y] = weighted.mean(x = obs_unc_out[pick], w = obs_lag_out[pick]+1, na.rm=TRUE)
                    obs_lag_agg[y] = sum(obs_lag_out[pick]+1, na.rm=TRUE)
                    # Then reaccumulate based on the total number of lags
                    obs_agg[y] = obs_agg[y] * obs_lag_agg[y]
                    obs_unc_agg[y] = obs_unc_agg[y] * obs_lag_agg[y]           
                }
           }       
       } else {
            stop("A non-valid function has been specified for the extract_timeseries_observations_with_uncertainty()")
       }

       # Convert the lag periods into model time steps
       obs_lag_agg = ceiling(obs_lag_agg / mean(timestep_days))

       # update with new output information
       obs_out = obs_agg ; obs_unc_out = obs_unc_agg ; obs_lag_out = obs_lag_agg
       # clean up
       rm(obs_agg,obs_unc_agg,obs_lag_agg,y) ; gc()

   } # temporal aggregation etc

   # convert missing data to -9999
   na_loc = which(is.na(obs_out) | is.na(obs_unc_out) | is.na(obs_lag_out))
   obs_out[na_loc] = na_flag ; obs_unc_out[na_loc] = na_flag ; obs_lag_out[na_loc] = na_flag

   # clean up
   rm(i1,j1,obs,unc,lag,i,a) ; gc(reset=TRUE,verbose=FALSE)

   # Create output object
   output = list(obs_out, obs_unc_out, obs_lag_out,na_loc = na_loc)
   # Update with the correct variable names
   names(output)[1:3]<-c(est_var_name_out,unc_var_name_out,lag_var_name_out)
   # Return function
   return(output)

} # end function extract_timeseries_observations_with_uncertainty
## Use byte compile
extract_timeseries_observations_with_uncertainty<-cmpfun(extract_timeseries_observations_with_uncertainty)

extract_timeseries_observations_without_uncertainty<- function(i1,j1,timestep_days,years_to_load,doy_obs,
                                                               data_all,agg_func, na_flag,
                                                               est_var_name_in,lag_var_name_in,
                                                               est_var_name_out,lag_var_name_out) {

   # Update the user
   if (use_parallel == FALSE) {print(paste("extracting ",est_var_name_in," for current location ",Sys.time(),sep=""))}

   # Extract current location to local variable
   obs = data_all[[est_var_name_in]][i1,j1,]
   lag = data_all[[lag_var_name_in]][i1,j1,]

   # declare output variable
   obs_out = array(NA, dim=length(doy_obs))
   obs_lag_out = array(NA, dim=length(doy_obs))

   ## Line up the days of year which have observations into a complete timeseries of days...
   b = 1 ; i = 1 ; a = 1 ; start_year = as.numeric(years_to_load[1])
   #print("...begin inserting LAI observations into model time steps")
   while (b <= length(data_all$doy_obs)) {

      # if we are in a year which is missing then we do not allow consideration of DOY
      if (start_year != data_all$missing_years[a]) {
          if (doy_obs[i] == data_all$doy_obs[b]) {
              obs_out[i] = obs[b] ; obs_lag_out[i] = lag[b] ; b = b + 1
          } # end if doy matches
      } # end if missing year

      # but we do keep counting through the total vector length which we expect
      i = i + 1

      # each time we come back to doy_obs[i]==1 we need to count on the year
      if (doy_obs[i] == 1 & b <= length(data_all$doy_obs)) {
          # and if we have just been in a missing year we need to count on the missing years vector to
          if (start_year == data_all$missing_years[a]) { a = min(length(data_all$missing_years),a+1) }
          start_year = start_year + 1
      } # end if doy_obs[i] == 1

   } # end while condition

   ## Aggregate to the model time steps by sum or mean
   if (length(timestep_days) == 1 & timestep_days[1] == 1) {

       # well actually we do nothing

   } else {

       # Aggregate when the time step is not daily

       # Sanity check in case we have not been 
       # given a complete timeseries of step sizes
       if (length(timestep_days) == 1) {
           run_day_selector = seq(1,length(obs_out),timestep_days)
           timestep_days = rep(timestep_days, length.out=length(obs_out))
       }

       # Determine the actual cumulative number of days to have passed
       run_day_selector = cumsum(timestep_days)
       # create needed variables
       obs_agg = array(NA,dim=length(run_day_selector))
       obs_lag_agg = array(NA,dim=length(run_day_selector))

       # Determine whether we are aggregating by mean or by sum
       if (agg_func == "mean") {
           # Loop through timeseries and aggregate
           for (y in seq(1,length(run_day_selector))) {
                # Only run for time steps when there is a value
                pick = (run_day_selector[y]-timestep_days[y]+1):run_day_selector[y]
                if (length(which(is.na(obs_out[pick]) == FALSE)) > 0) {       
                    obs_agg[y] = weighted.mean(x = obs_out[pick], w = pmax(1,obs_lag_out[pick], na.rm=TRUE), na.rm=TRUE)
                    obs_lag_agg[y] = sum(pmax(1,obs_lag_out[pick],na.rm=TRUE), na.rm=TRUE)
                }
           }
       } else if (agg_func == "sum") {
           # Loop through timeseries and aggregate
           for (y in seq(1,length(run_day_selector))) {
                # Only run for time steps when there is a value
                pick = (run_day_selector[y]-timestep_days[y]+1):run_day_selector[y]
                if (length(which(is.na(obs_out[pick]) == FALSE)) > 0) {       
                    obs_agg[y] = weighted.mean(x = obs_out[pick], w = obs_lag_out[pick]+1, na.rm=TRUE)                          
                    obs_lag_agg[y] = sum(obs_lag_out[pick]+1, na.rm=TRUE)
                    # Then reaccumulate based on the total number of lags
                    obs_agg[y] = obs_agg[y] * obs_lag_agg[y]
                } 
           }    
       } else {
            stop("A non-valid function has been specified for the extract_timeseries_observations_without_uncertainty()")
       }
       # Convert the lag periods into model time steps
       obs_lag_agg = ceiling(obs_lag_agg / mean(timestep_days))

       # update with new output information
       obs_out = obs_agg ; obs_lag_out = obs_lag_agg
       # clean up
       rm(obs_agg,y) ; gc()

   } # temporal aggregation etc

   # convert missing data to -9999
   na_loc = which(is.na(obs_out) | is.na(obs_lag_out)) 
   obs_out[na_loc] = na_flag ; obs_lag_out[na_loc] = na_flag

   # clean up
   rm(i1,j1,obs,lag,i,a) ; gc(reset=TRUE,verbose=FALSE)

   # Create output object
   output = list(obs_out,obs_lag_out,na_loc = na_loc)
   # Update with the correct variable names
   names(output)[1:2]<-c(est_var_name_out,lag_var_name_out)
   # Return function
   return(output)

} # end function extract_timeseries_observations_without_uncertainty
## Use byte compile
extract_timeseries_observations_without_uncertainty<-cmpfun(extract_timeseries_observations_without_uncertainty)

extract_static_observations_with_uncertainty<- function(i1,j1,data_all,na_flag,
                                                        est_var_name_in,unc_var_name_in,
                                                        est_var_name_out,unc_var_name_out) {

   # Update the user
   if (use_parallel == FALSE) {print(paste("extracting ",est_var_name_in," for current location ",Sys.time(),sep=""))}

   # Extract current location to local variable
   obs = data_all[[est_var_name_in]][i1,j1]
   unc = data_all[[unc_var_name_in]][i1,j1]

   # convert missing data to na_flag (typically -9999)
   obs[which(is.na(obs))] = na_flag ; unc[which(is.na(unc))] = na_flag

   # Create output object
   output = list(obs, unc)
   # Update with the correct variable names
   names(output)[1:2]<-c(est_var_name_out,unc_var_name_out)
   # Return function
   return(output)

} # end function extract_observations_with_uncertainty
## Use byte compile
extract_static_observations_with_uncertainty<-cmpfun(extract_static_observations_with_uncertainty)

extract_static_observations_without_uncertainty<- function(i1,j1,data_all,na_flag,
                                                           est_var_name_in,est_var_name_out) {

   # Update the user
   if (use_parallel == FALSE) {print(paste("extracting ",est_var_name_in," for current location ",Sys.time(),sep=""))}

   # Extract current location to local variable
   obs = data_all[[est_var_name_in]][i1,j1]

   # convert missing data to -9999
   obs[which(is.na(obs))] = na_flag 

   # Create output object
   output = list(obs)
   # Update with the correct variable names
   names(output)[1]<-c(est_var_name_out)
   # Return function
   return(output)

} # end function extract_observations_with_uncertainty
## Use byte compile
extract_static_observations_without_uncertainty<-cmpfun(extract_static_observations_without_uncertainty)

extract_timeseries_forcing<- function(i1,j1,timestep_days,years_to_load,doy_obs,
                                      data_all,agg_func,fraction,na_flag,
                                      est_var_name_in,lag_var_name_in,
                                      est_var_name_out) {

   # This subroutine is specifically intended for extracting non-meteorological forcings.
   # We assume no uncertainties.
   # Lags are now assumed to be used to distribute the forcings across a number of days preceeding the specified date.
   # NOTE / WARNING: We assume that the lags are non-overlapping. 

   # Update the user
   if (use_parallel == FALSE) {print(paste("extracting ",est_var_name_in," for current location ",Sys.time(),sep=""))}

   # Extract current location to local variable
   obs = data_all[[est_var_name_in]][i1,j1,]
   lag = data_all[[lag_var_name_in]][i1,j1,]

   # declare output variable
   obs_out = array(NA, dim=length(doy_obs))
   obs_lag_out = array(NA, dim=length(doy_obs))

   ## Line up the days of year which have observations into a complete timeseries of days...
   b = 1 ; i = 1 ; a = 1 ; start_year = as.numeric(years_to_load[1])
   #print("...begin inserting LAI observations into model time steps")
   while (b <= length(data_all$doy_obs)) {

      # if we are in a year which is missing then we do not allow consideration of DOY
      if (start_year != data_all$missing_years[a]) {
          if (doy_obs[i] == data_all$doy_obs[b]) {
              obs_out[i] = obs[b] ; obs_lag_out[i] = lag[b] ; b = b + 1
          } # end if doy matches
      } # end if missing year

      # but we do keep counting through the total vector length which we expect
      i = i + 1

      # each time we come back to doy_obs[i]==1 we need to count on the year
      if (doy_obs[i] == 1 & b <= length(data_all$doy_obs)) {
          # and if we have just been in a missing year we need to count on the missing years vector to
          if (start_year == data_all$missing_years[a]) { a = min(length(data_all$missing_years),a+1) }
          start_year = start_year + 1
      } # end if doy_obs[i] == 1

   } # end while condition

   ## Aggregate to the model time steps by sum or mean
   if (length(timestep_days) == 1 & timestep_days[1] == 1) {

       ## Lags for distributing a variable backwards in time must now be applied before any temporal aggregation
       for (y in seq(1, length(obs_out))) {
            # If the lag is == zero then we don't need to do anything as its is instantaneously applied
            if (obs_lag_out[y] != 0) {
                if (fraction) {
                    # If this is a fraction, i.e. not an absolute rate, 
                    # then we do a time varied average based on compound 
                    # interest calculation to preserve the mass balance
                    # If lag is not zero then we should go through and distribute
                    obs_out[max(1,y-obs_lag_out[y]):y] = 1-(1-obs_out[y])**(1/obs_lag_out[y])
                    # We also must update the lag now to be instantanuous
                    obs_lag_out[max(1,y-obs_lag_out[y]):y] = 0
                } else {
                    # If this is an absolute, i.e. not a fractional rate, then we do a direct averaging
                    # to preserve the mass balance
                    # If lag is not zero then we should go through and distribute
                    obs_out[max(1,y-obs_lag_out[y]):y] = obs_out[y] / obs_lag_out[y]
                    # We also must update the lag now to be instantanuous
                    obs_lag_out[max(1,y-obs_lag_out[y]):y] = 0
                } 
           }
       } #  lag distributing loop

   } else {

       # Aggregate when the time step is not daily

       # Sanity check in case we have not been 
       # given a complete timeseries of step sizes
       if (length(timestep_days) == 1) {
           run_day_selector = seq(1,length(obs_out),timestep_days)
           timestep_days = rep(timestep_days, length.out=length(obs_out))
       }

       # Determine the actual cumulative number of days to have passed
       run_day_selector = cumsum(timestep_days)
       # create needed variables
       obs_agg = array(NA,dim=length(run_day_selector))
       obs_lag_agg = array(NA,dim=length(run_day_selector))

       # Determine whether we are aggregating by mean or by sum
       if (agg_func == "mean") {
           # Loop through timeseries and aggregate
           for (y in seq(1,length(run_day_selector))) {
                # Only run for time steps when there is a value
                pick = (run_day_selector[y]-timestep_days[y]+1):run_day_selector[y]
                if (length(which(is.na(obs_out[pick]) == FALSE)) > 0) {                           
                    obs_agg[y] = weighted.mean(x = obs_out[pick], w = pmax(1,obs_lag_out[pick], na.rm=TRUE), na.rm=TRUE)
                    obs_lag_agg[y] = sum(obs_lag_out[pick], na.rm=TRUE)
                }
           }
       } else if (agg_func == "sum") {
           # Loop through timeseries and aggregate
           for (y in seq(1,length(run_day_selector))) {
                # Only run for time steps when there is a value
                pick = (run_day_selector[y]-timestep_days[y]+1):run_day_selector[y]
                if (length(which(is.na(obs_out[pick]) == FALSE)) > 0) {                  
                    # Ensure a minimum weighting value of 1 to ensure that instantanous values are counted.         
                    obs_agg[y] = weighted.mean(x = obs_out[pick], w = pmax(1,obs_lag_out[pick], na.rm=TRUE), na.rm=TRUE)                          
                    # Accumuate the weightings to allow the mass balance to be estimated below.
                    obs_lag_agg[y] = sum(pmax(1,obs_lag_out[pick], na.rm=TRUE), na.rm=TRUE)
                    # Then reaccumulate based on the total number of lags
                    obs_agg[y] = obs_agg[y] * obs_lag_agg[y]
                } 
           }    
           # Final sanity check to ensure that time steps without values do not erroneous get turned into a -1 which is an impossible lag
           #obs_lag_agg = pmax(0,obs_lag_agg)
       } else {
            stop("A non-valid function has been specified for the extract_timeseries_forcing()")
       }

       # Convert the lag periods into model time steps
       obs_lag_agg = ceiling(obs_lag_agg / mean(timestep_days))        

       # As these are intended to be forcings we need to now distribute the values over the lagged periods.
       # We will use the inverse of the compound interest approach to achieve this.

       # For the following loop all NA lags must now be set to zero
       obs_lag_agg[is.na(obs_lag_agg)] = 0

       ## Lags for distributing a variable backwards in time
       if (sum(obs_lag_agg, na.rm=TRUE) != 0) {
           for (y in seq(1, length(obs_agg))) {
                # If the lag is == zero then we don't need to do anything as its is instantaneously applied
                if (obs_lag_agg[y] != 0) {
                    if (fraction) {
                        # Sanity check 
                        if (obs_agg[y] > 1) {
                            print("WARNING: the extraction of a fractional forcing estimate which lagged sum is greater than 1.\\
                                            This will likely result in the value being lost to NaN.")
                        }
                        # If this is a fraction, i.e. not an absolute rate, 
                        # then we do a time varied average based on compound 
                        # interest calculation to preserve the mass balance
                        # If lag is not zero then we should go through and distribute
                        obs_agg[max(1,y-obs_lag_agg[y]):y] = 1-(1-obs_agg[y])**(1/obs_lag_agg[y])
                        # We also must update the lag now to be instantanuous
                        obs_lag_agg[max(1,y-obs_lag_agg[y]):y] = 0
                    } else {
                        # If this is an absolute, i.e. not a fractional rate, then we do a direct averaging
                        # to preserve the mass balance
                        # If lag is not zero then we should go through and distribute
                        obs_agg[max(1,y-obs_lag_agg[y]):y] = obs_agg[y] / obs_lag_agg[y]
                        # We also must update the lag now to be instantanuous
                        obs_lag_agg[max(1,y-obs_lag_agg[y]):y] = 0
                    } 
                }
           } #  lag distributing loop  
       } # any lags at all?     

       # update with new output information
       obs_out = obs_agg ; obs_lag_out = obs_lag_agg
       # clean up
       rm(obs_agg,y) ; gc()

   } # temporal aggregation etc

   # convert missing data to -9999
   na_loc = which(is.na(obs_out)) 
   obs_out[na_loc] = na_flag 

   # clean up
   rm(i1,j1,obs,lag,i,a) ; gc(reset=TRUE,verbose=FALSE)

   # Create output object
   output = list(obs_out,na_loc = na_loc)
   # Update with the correct variable names
   names(output)[1]<-c(est_var_name_out)
   # Return function
   return(output)

} # end function extract_timeseries_forcing
## Use byte compile
extract_timeseries_forcing<-cmpfun(extract_timeseries_forcing)