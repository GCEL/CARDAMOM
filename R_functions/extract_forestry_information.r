###
## Function deals with the extraction of relevant forest clearance and growth information for UK forestry only!
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

extract_forestry_information<-function(i1,j1,timestep_days,spatial_type,resolution,grid_type,
                                       latlon_in,forest_all,start_year,end_year,ctessel_pft_in,
                                       years_to_load,doy_obs) {

#   # find the nearest location
#   output=closest2d_2(1,forest_all$lat,forest_all$long,latlon_in[1],latlon_in[2])
#   i1=unlist(output, use.names=FALSE)[1] ; j1=unlist(output, use.names=FALSE)[2]

   # Assume this location does not have forest commission information
   # The age, yield class and pft override will be removed at a later date
   # TLS: 11/01/2022
   age=-9999
   yield_class=-9999
   ctessel_pft=ctessel_pft_in

   # declare output variable
   deforestation = rep(0, times=length(doy_obs))
   ## normal assumption
   start_of_years = which(doy_obs == 340)
   # which year is the one in which deforestation occurs?
   # then find the appropriate beginning of a year and make deforestation
   for (aa in seq(1,length(forest_all$year_of_loss))) {
        start_point = start_of_years[which(as.numeric(years_to_load) == forest_all$year_of_loss[aa])]  
       # end_point = start_point 
        deforestation[start_point] = (forest_all$loss_fraction[i1,j1,aa])  
   }

   # generally this now deals with time steps which are not daily.
   # However if not monthly special case
   if (length(timestep_days) == 1) {
       run_day_selector=seq(1,length(deforestation),timestep_days)
       timestep_days=rep(timestep_days, length.out=length(deforestation))
   }

   if (use_parallel == FALSE) {print("...calculating weekly / monthly averages for deforestation info")}
   # determine the actual daily positions
   run_day_selector = cumsum(timestep_days)
   # create needed variables
   deforestation_agg = array(NA,dim=length(run_day_selector))
   for (y in seq(1,length(run_day_selector))) {
        deforestation_agg[y] = sum(deforestation[(run_day_selector[y]-timestep_days[y]):run_day_selector[y]],na.rm=TRUE)
        # having picked from this period, ensure no overlap by clearing it!
        deforestation[(run_day_selector[y]-timestep_days[y]):run_day_selector[y]] = NA
   }
   deforestation_agg[is.na(deforestation_agg)] = -9999
   # update with new output information
   deforestation = deforestation_agg
   # clean up
   rm(deforestation_agg) ; gc()

   # return time series and updates pft information
   return(list(ctessel_pft=ctessel_pft,deforestation=deforestation,yield_class=yield_class,age=age))

} # end function extract_forestry_information

## Use byte compile
extract_forestry_information<-cmpfun(extract_forestry_information)
