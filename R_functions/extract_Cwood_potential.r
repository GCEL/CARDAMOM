

###
## Function extracts information on the potential steady state value of the wood stock pool
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

extract_Cwood_potential<- function(timestep_days,spatial_type,resolution,grid_type,latlon_in,Cwood_potential_all) {

   # Update the user
   print(paste("Cwood potential extracted for current location ",Sys.time(),sep=""))

   # find the nearest location
   output = closest2d(1,Cwood_potential_all$lat,Cwood_potential_all$long,latlon_in[1],latlon_in[2],2)
   i1 = unlist(output)[1] ; j1 = unlist(output)[2]

   # If resolution has been provides as single value then adjust this here
   if (length(resolution) == 1) {tmp_res = resolution * c(1,1)} else {tmp_res = resolution}

   # work out number of pixels to average over
   if (spatial_type == "grid") {
       # resolution of the product
       product_res = c(abs(Cwood_potential_all$long[2,1]-Cwood_potential_all$long[1,1]),abs(Cwood_potential_all$lat[1,2]-Cwood_potential_all$lat[1,1]))
       if (grid_type == "wgs84") {
           # radius is ceiling of the ratio of the product vs analysis ratio
           radius = floor(0.5*(resolution / product_res))
        } else if (grid_type == "UK") {
           # Estimate radius for UK grid assuming radius is determine by the longitude size
           # 6371e3 = mean earth radius (m)
           radius = round(rad2deg(sqrt((resolution / 6371e3**2))) / product_res, digits=0)
           #radius = max(0,floor(1*resolution*1e-3*0.5))
        } else {
           stop("have not specified the grid used in this analysis")
        }
   } else {
        radius = c(0,0)
   }

   # Work out average areas
   average_i = (i1-radius[1]):(i1+radius[1]) ; average_j = (j1-radius[2]):(j1+radius[2])
   average_i = max(1,(i1-radius[1])):min(dim(Cwood_potential_all$biomass_gCm2)[1],(i1+radius[1]))
   average_j = max(1,(j1-radius[2])):min(dim(Cwood_potential_all$biomass_gCm2)[2],(j1+radius[2]))
   # Carry out averaging
   Cwood = mean(Cwood_potential_all$biomass_gCm2[average_i,average_j], na.rm=TRUE)
   Cwood_unc = mean(Cwood_potential_all$biomass_uncertainty_gCm2[average_i,average_j], na.rm=TRUE)

   # Convert any NaN to missing data flag -9999
   Cwood[which(is.na(Cwood))] = -9999 ; Cwood_unc[which(is.na(Cwood_unc))] = -9999

   # pass the information back
   return(list(Cwood_stock = Cwood, Cwood_stock_unc = Cwood_unc))

} # end function extract_Cwood_potential

## Use byte compile
extract_Cwood_potential<-cmpfun(extract_Cwood_potential)