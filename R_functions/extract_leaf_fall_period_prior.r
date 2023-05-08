
###
## Extracts location specific information on the leaf fall period from a loaded gridded dataset
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

extract_leaf_fall_period_prior<- function(i1,j1,spatial_type,resolution,grid_type,latlon_in,leaf_fall_period_all) {

   # Update the user
   print(paste("Leaf fall period prior extracted for current location ",Sys.time(),sep=""))

#   # find the nearest location
#   output = closest2d_2(1,lifespan_all$lat,lifespan_all$long,latlon_in[1],latlon_in[2])
#   i1 = unlist(output, use.names=FALSE)[1] ; j1 = unlist(output, use.names=FALSE)[2]

   # Extract target location
   leaf_fall_period = leaf_fall_period_all$leaf_fall_period[i1,j1]
   #leaf_fall_period_unc_days = leaf_fall_period_all$leaf_fall_period_uncertainty_days[i1,j1]

   # Convert any NaN to missing data flag -9999
   leaf_fall_period[is.na(leaf_fall_period)] = -9999 #; leaf_fall_period_unc[is.na(leaf_fall_period_unc_days)] = -9999

   # pass the information back
   return(list(leaf_fall_period = leaf_fall_period))#, leaf_fall_period_unc_days = leaf_fall_period_unc_days))

} # end function extract_leaf_fall_period_prior

## Use byte compile
extract_leaf_fall_period_prior<-cmpfun(extract_leaf_fall_period_prior)