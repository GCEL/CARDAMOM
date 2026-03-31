
###
## Extracts location specific information on the leaf max_root (yrs) from a loaded gridded dataset
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

extract_max_root_prior<- function(i1,j1,spatial_type,resolution,grid_type,latlon_in,max_root_all) {

   # Update the user
   #print(paste("max_root prior extracted for current location ",Sys.time(),sep=""))

#   # find the nearest location
#   output = closest2d_2(1,max_root_all$lat,max_root_all$long,latlon_in[1],latlon_in[2])
#   i1 = unlist(output, use.names=FALSE)[1] ; j1 = unlist(output, use.names=FALSE)[2]

   # Extract target location
   max_root = max_root_all$max_root[i1,j1]
   #max_root_unc_yrs = max_root_all$max_root_uncertainty_yrs[i1,j1]

   # Convert any NaN to missing data flag -9999
   max_root[is.na(max_root)] = -9999 #; max_root_unc[is.na(max_root_unc_yrs)] = -9999

   # pass the information back
   return(list(max_root = max_root))#, max_root_unc_yrs = max_root_unc_yrs))

} # end function extract_max_root_prior

## Use byte compile
extract_max_root_prior<-cmpfun(extract_max_root_prior)