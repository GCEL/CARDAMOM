
###
## Function extracts locatio specific information on Burned area
## from already loaded gridded datasets
###

# This function is based by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

extract_burnt_area_information<- function(i1,j1,latlon_in,timestep_days,spatial_type,
                                          grid_type,resolution,start_year,end_year,burnt_all,
                                          years_to_load,doy_obs) {

  # Update the user
  if (use_parallel == FALSE) {print(paste("Beginning burned fraction data extraction for current location ",Sys.time(),sep=""))}

  # Extract location specific information to a local variable
  burnt_area = burnt_all$burnt_area[i1,j1,]

  # just incase there is no missing data we best make sure there is a value which can be assessed
  if (length(burnt_all$missing_years) == 0) { burnt_all$missing_years=1066 }

  # declare output variable
  burnt_area_out = array(0, dim=length(doy_obs))
  # now line up the obs days with all days
  b = 1 ; i = 1 ; a = 1 ; start_year = as.numeric(years_to_load[1])
  while (b <= length(burnt_all$doy_obs)) {
         # if we are in a year which is missing then we do not allow consideration of DOY
         if (start_year != burnt_all$missing_years[a]) {
             if (doy_obs[i] == burnt_all$doy_obs[b]) {
                 burnt_area_out[i] = burnt_area[b] ; b = b+1
             } # end if doy matches
         } # end if missing year
         # but we do keep counting through the total vector length which we expect
         i = i + 1
         # have we just looped round the year?
         if (i != 1 & doy_obs[i-1] > doy_obs[i]) {
             # and if we have just been in a missing year we need to count on the missing years vector to
             if (start_year == burnt_all$missing_years[a]) {a = min(length(burnt_all$missing_years),a+1)}
             start_year = start_year+1
         } # end if doy_obs[i] == 1
  } # end while condition

  if (length(timestep_days) == 1 & timestep_days[1] == 1) {
      # well actually we do nothing
  } else {
      # generally this now deals with time steps which are not daily.
      # However if not monthly special case
      if (length(timestep_days) == 1) {
          run_day_selector = seq(1,length(burnt_area_out),timestep_days)
          timestep_days = rep(timestep_days, length.out=length(burnt_area_out))
      }
      #print("...calculating monthly averages for burnt area")
      # determine the actual daily positions
      run_day_selector = cumsum(timestep_days)
      # create needed variables
      burnt_area_agg = array(NA,dim=length(run_day_selector))
      # As far is an absolute event we sum all events within the time period not average!
      for (y in seq(1,length(run_day_selector))) {
           burnt_area_agg[y] = sum(burnt_area_out[(run_day_selector[y]-timestep_days[y]+1):run_day_selector[y]],na.rm=TRUE)
      }
      # update with new output information
      burnt_area_out = burnt_area_agg
      # clean up
      rm(burnt_area_agg) ; gc()
  } # monthly aggregation etc

  # convert missing data to -9999
  burnt_area_out[is.na(burnt_area_out) == TRUE] = -9999

  # pass the information back
  return(burnt_area=burnt_area_out)

} # end function extract_burnt_area_information

## Use byte compile
extract_burnt_area_information<-cmpfun(extract_burnt_area_information)
