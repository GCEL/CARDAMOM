
###
## Extracts location specific information on the leaf lifespan (yrs) from a loaded gridded dataset
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

extract_lifespan_prior<- function(i1,j1,spatial_type,resolution,grid_type,latlon_in,lifespan_all) {

   # Update the user
   print(paste("Lifespan prior extracted for current location ",Sys.time(),sep=""))

#   # find the nearest location
#   output = closest2d_2(1,lifespan_all$lat,lifespan_all$long,latlon_in[1],latlon_in[2])
#   i1 = unlist(output, use.names=FALSE)[1] ; j1 = unlist(output, use.names=FALSE)[2]

   # Extract target location
   lifespan = lifespan_all$lifespan[i1,j1]
   #lifespan_unc_yrs = lifespan_all$lifespan_uncertainty_yrs[i1,j1]

   # Convert any NaN to missing data flag -9999
   lifespan[is.na(lifespan)] = -9999 #; lifespan_unc[is.na(lifespan_unc_yrs)] = -9999

   # pass the information back
   return(list(lifespan = lifespan))#, lifespan_unc_yrs = lifespan_unc_yrs))

} # end function extract_lifespan_prior

## Use byte compile
extract_lifespan_prior<-cmpfun(extract_lifespan_prior)