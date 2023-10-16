
###
## Function extracts fraction absorbed photosynthetically active radiation (fAPAR)
## information from a pre-loaded gridded dataset
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

extract_fapar_timeseries<- function(i1,j1,timestep_days,spatial_type,resolution,
                                    grid_type,latlon_in,fapar_all,years_to_load,doy_obs) {

   # Update the user
   if (use_parallel == FALSE) {print(paste("fAPAR data extracted for current location ",Sys.time(),sep=""))}

   # Extract current location to local variable
   fapar = fapar_all$fapar_all[i1,j1,]
   fapar_unc = fapar_all$fapar_unc_all[i1,j1,]

   # Just incase there is no missing data we best make sure there is a value which can be assessed
   if (length(fapar_all$missing_years) == 0) { fapar_all$missing_years=1066 }

   # declare output variable
   fapar_out = array(NA, dim=length(doy_obs))
   fapar_unc_out = array(NA, dim=length(doy_obs))
   # now line up the obs days with all days
   b = 1 ; i = 1 ; a = 1 ; start_year = as.numeric(years_to_load[1])
   #print("...begin inserting fAPAR observations into model time steps")
   while (b <= length(fapar_all$doy_obs)) {

      # if we are in a year which is missing then we do not allow consideration of DOY
      if (start_year != fapar_all$missing_years[a]) {
          if (doy_obs[i] == fapar_all$doy_obs[b]) {
              fapar_out[i] = fapar[b] ; fapar_unc_out[i] = fapar_unc[b] ; b = b + 1
          } # end if doy matches
      } # end if missing year

      # but we do keep counting through the total vector length which we expect
      i = i + 1

      # each time we come back to doy_obs[i]==1 we need to count on the year
      if (doy_obs[i] == 1 & b <= length(fapar_all$doy_obs)) {
          # and if we have just been in a missing year we need to count on the missing years vector to
          if (start_year == fapar_all$missing_years[a]) { a = min(length(fapar_all$missing_years),a+1) }
          start_year = start_year + 1
      } # end if doy_obs[i] == 1

   } # end while condition

   if (length(timestep_days) == 1 & timestep_days[1] == 1) {

       # well actually we do nothing

   } else {
       # generally this now deals with time steps which are not daily.
       # However if not monthly special case
       if (length(timestep_days) == 1) {
           run_day_selector=seq(1,length(fapar_out),timestep_days)
           timestep_days=rep(timestep_days, length.out=length(fapar_out))
       }
       #print("...calculating monthly averages for lai")
       # determine the actual daily positions
       run_day_selector=cumsum(timestep_days)
       # create needed variables
       fapar_agg = array(NA,dim=length(run_day_selector))
       fapar_unc_agg = array(NA,dim=length(run_day_selector))
       # Loop through
       for (y in seq(1,length(run_day_selector))) {
            pick = (run_day_selector[y]-timestep_days[y]+1):run_day_selector[y]
            fapar_agg[y] = mean(fapar_out[pick], na.rm=TRUE)
            fapar_unc_agg[y] = mean(fapar_unc_out[pick], na.rm=TRUE)
       }
       # update with new output information
       fapar_out = fapar_agg ; fapar_unc_out = fapar_unc_agg
       # clean up
       rm(fapar_agg,fapar_unc_agg,y) ; gc()

   } # monthly aggregation etc

   # convert missing data to -9999
   fapar_out[which(is.na(fapar_out))] = -9999 ; fapar_unc_out[which(is.na(fapar_unc_out))] = -9999

   # clean up
   rm(i1,j1,fapar,i,a) ; gc(reset=TRUE,verbose=FALSE)

   # CARDAMOM works best if the uncertainties are the same across each LAI observation as the framework tends towards lower LAI values
   # Therefore, to make use of the uncertainty information we take the mean for this site and apply it across each value.
   # NOTE: we put a book end the upper uncertainty linked to the max LAI estimate to ensure that there is some constraint
   #fapar_unc_out[lai_out >= 0] = max(0.25,min(mean(lai_unc_out[lai_out >= 0]), max(lai_out[lai_out >= 0])))
   #fapar_unc_out[fapar_out >= 0] = pmax(0.25,fapar_unc_out[fapar_out >= 0])

   # pass the information back
   output = list(fapar = fapar_out, fapar_unc = fapar_unc_out)
   return(output)

} # end function extract_lai_timeseries

## Use byte compile
extract_fapar_timeseries<-cmpfun(extract_fapar_timeseries)
