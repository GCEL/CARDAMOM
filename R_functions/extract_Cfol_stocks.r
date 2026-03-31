
###
## Function extracts location specific information on the timeseries of Cfol stock information
## drawn from an already loaded gridded dataset
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

extract_Cfol_stocks<- function(i1,j1,timestep_days,spatial_type,resolution,grid_type,
                                latlon_in,Cfol_stock_all) {

   # Update the user
   if (use_parallel == FALSE) {print(paste("Cfol stocks extracted for current location ",Sys.time(),sep=""))}

#   # find the nearest location
#   output = closest2d_2(1,Cfol_stock_all$lat,Cfol_stock_all$long,latlon_in[1],latlon_in[2])
#   i1 = unlist(output, use.names=FALSE)[1] ; j1 = unlist(output, use.names=FALSE)[2]

   # Create time series output variables
   Cfol_stock = rep(-9999, length(timestep_days))
   Cfol_stock_unc = rep(-9999, length(timestep_days))

   # Loop through each time step of the Cfol time series obs and
   # estimate average value
   for (t in seq(1, length(Cfol_stock_all$place_obs_in_step))) {
        Cfol_stock[Cfol_stock_all$place_obs_in_step[t]+6] = Cfol_stock_all$foliage_gCm2[i1,j1,t]
        tmp = min(Cfol_stock[Cfol_stock_all$place_obs_in_step[t]+6], Cfol_stock_all$foliage_uncertainty_gCm2[i1,j1,t])
        Cfol_stock_unc[Cfol_stock_all$place_obs_in_step[t]+6] = tmp
   }

   # Set any time series values with NaN to missing data flag (-9999)
   Cfol_stock[which(is.na(Cfol_stock))] = -9999
   Cfol_stock_unc[which(is.na(Cfol_stock_unc))] = -9999

   # pass the information back
   return(list(Cfol_stock = Cfol_stock, Cfol_stock_unc = Cfol_stock_unc))
   #return(list(Cfol_stock = Cfol_stock, Cfol_stock_unc = 250))

} # end function extract_Cfol_stocks

## Use byte compile
extract_Cfol_stocks<-cmpfun(extract_Cfol_stocks)
