
###
## Function used to extract location specific information on GPP from already loaded gridded dataset
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

extract_gpp<- function(i1,j1,timestep_days,spatial_type,resolution,grid_type,
                       latlon_in,gpp_all,years_to_load,doy_obs) {

  # Update the user
  if (use_parallel == FALSE) {print(paste("GPP data extracted for current location ",Sys.time(),sep=""))}

#  # find the nearest location
#  output = closest2d_2(1,gpp_all$lat,gpp_all$long,latlon_in[1],latlon_in[2])
#  i1 = unlist(output, use.names=FALSE)[1] ; j1 = unlist(output, use.names=FALSE)[2]

  # Extract to local variable
  gpp = gpp_all$gpp_gCm2day[i1,j1,]
  gpp_unc = gpp_all$gpp_unc_gCm2day[i1,j1,]

  # just incase there is no missing data we best make sure there is a value which can be assessed
  if (length(gpp_all$missing_years) == 0) { gpp_all$missing_years=1066 }

  # declare output variable
  gpp_out = array(NA, dim=length(doy_obs))
  gpp_unc_out = array(NA, dim=length(doy_obs))
  # now line up the obs days with all days
  b = 1 ; i = 1 ; a = 1 ; start_year=as.numeric(years_to_load[1])
  while (b <= length(gpp_all$doy_obs)) {

      # if we are in a year which is missing then we do not allow consideration of DOY
      if (start_year != gpp_all$missing_years[a]) {
          if (doy_obs[i] == gpp_all$doy_obs[b]) {
              gpp_out[i] = gpp[b] ; gpp_unc_out[i] = gpp_unc[b] ; b = b+1
          } # end if doy matches
      } # end if missing year
      # but we do keep counting through the total vector length which we expect
      i = i+1

      # each time we come back to doy_obs[i]==1 we need to count on the year
      if (doy_obs[i-1] > doy_obs[i] & b <= length(gpp_all$doy_obs)) {
          # and if we have just been in a missing year we need to count on the missing years vector to
          if (start_year == gpp_all$missing_years[a]) {a = min(length(gpp_all$missing_years),a+1)}
          start_year=start_year+1
      } # end if doy_obs[i] == 1

  } # end while condition

  if (length(timestep_days) == 1 & timestep_days[1] == 1) {

      # well actually we do nothing

  } else {

      # Generally this now deals with time steps which are not daily.
      # However if not monthly special case
      if (length(timestep_days) == 1) {
          run_day_selector = seq(1,length(gpp_out),timestep_days)
          timestep_days = rep(timestep_days, length.out=length(gpp_out))
      }
      #print("...calculating monthly averages for GPP")
      # determine the actual daily positions
      run_day_selector = cumsum(timestep_days)
      # create needed variables
      gpp_agg = array(NA,dim=length(run_day_selector))
      gpp_unc_agg = array(NA, dim=length(run_day_selector))
      for (y in seq(1,length(run_day_selector))) {
           pick = (run_day_selector[y]-timestep_days[y]+1):run_day_selector[y]
           gpp_agg[y] = mean(gpp_out[pick],na.rm=TRUE)
           gpp_unc_agg[y] = mean(gpp_unc_out[pick],na.rm=TRUE)
      }
      # update with new output information
      gpp_out = gpp_agg ; gpp_unc_out = gpp_unc_agg
      # clean up
      rm(gpp_agg,gpp_unc_agg,y) ; gc()

  } # monthly aggregation etc

  # convert missing data back to -9999
  gpp_out[which(is.na(gpp_out))] = -9999
  gpp_unc_out[which(is.na(gpp_unc_out))] = -9999

  # pass the information back
  output = list(GPP = gpp_out, GPP_unc = gpp_unc_out)
  return(output)

} # end function extract_gpp

## Use byte compile
extract_gpp<-cmpfun(extract_gpp)
