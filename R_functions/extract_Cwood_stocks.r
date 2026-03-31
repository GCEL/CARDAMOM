
###
## Function extracts location specific information on the timeseries of wood stock information
## drawn from an already loaded gridded dataset
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

extract_Cwood_stocks<- function(i1,j1,timestep_days,spatial_type,resolution,grid_type,
                                latlon_in,Cwood_stock_all) {

   # Update the user
   if (use_parallel == FALSE) {print(paste("Cwood stocks extracted for current location ",Sys.time(),sep=""))}

#   # find the nearest location
#   output = closest2d_2(1,Cwood_stock_all$lat,Cwood_stock_all$long,latlon_in[1],latlon_in[2])
#   i1 = unlist(output, use.names=FALSE)[1] ; j1 = unlist(output, use.names=FALSE)[2]

   # Create time series output variables
   Cwood_stock = rep(-9999, length(timestep_days))
   Cwood_stock_unc = rep(-9999, length(timestep_days))

   # Loop through each time step of the Cwood time series obs and
   # estimate average value
   for (t in seq(1, length(Cwood_stock_all$place_obs_in_step))) {
        Cwood_stock[Cwood_stock_all$place_obs_in_step[t]+6] = Cwood_stock_all$biomass_gCm2[i1,j1,t]
        tmp = min(Cwood_stock[Cwood_stock_all$place_obs_in_step[t]+6], Cwood_stock_all$biomass_uncertainty_gCm2[i1,j1,t])
        Cwood_stock_unc[Cwood_stock_all$place_obs_in_step[t]+6] = tmp
   }

   # Set any time series values with NaN to missing data flag (-9999)
   Cwood_stock[which(is.na(Cwood_stock))] = -9999
   Cwood_stock_unc[which(is.na(Cwood_stock_unc))] = -9999

   # pass the information back
   return(list(Cwood_stock = Cwood_stock, Cwood_stock_unc = Cwood_stock_unc))
   #return(list(Cwood_stock = Cwood_stock, Cwood_stock_unc = 250))

} # end function extract_Cwood_stocks

## Use byte compile
extract_Cwood_stocks<-cmpfun(extract_Cwood_stocks)
