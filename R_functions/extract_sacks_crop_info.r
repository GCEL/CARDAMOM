
####
## Function to extract the a prior estimate of crop harvest and sowing data from the SACKS database
####

# This function is based on an original Python function development by D. Slevin (UoE, now at Forestry Commission, Scotland).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

extract_sacks_crop_info<- function(spatial_type,resolution,grid_type,latlon_in,crop_man_all) {

	# convert input data long to conform to what we need
	check1=which(crop_man_all$long > 180) ; if (length(check1) > 0) { crop_man_all$long[check1]=crop_man_all$long[check1]-360 }

	# find the nearest location
	output=closest2d(1,crop_man_all$lat,crop_man_all$long,latlon_in[1],latlon_in[2],3)
	#i1=unlist(output)[2] ; j1=unlist(output)[1]
	i1=unlist(output)[1] ; j1=unlist(output)[2]
	print(paste("Crop management data extracted for current location ",Sys.time(),sep=""))

	# return long to 0-360
	if (length(check1) > 0) { crop_man_all$long[check1]=crop_man_all$long[check1]+360 }

	# work out number of pixels to average over
  if (spatial_type == "grid") {
      # resolution of the product
      product_res = abs(lai_all$lat[1,2]-lai_all$lat[1,1])+abs(lai_all$long[2,1]-lai_all$long[1,1])
      product_res = product_res * 0.5 # NOTE: averaging needed for line above
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
      radius = 0
  }


	answer=NA ; plant = -9999 ; harvest = -9999 ; plant_range = -9999 ; harvest_range = -9999
	while (is.na(answer) == TRUE) {
	    # work out average areas
	    average_i=max(1,i1-radius):min(length(crop_man_all$long),i1+radius) ; average_j=max(1,j1-radius):min(length(crop_man_all$lat),j1+radius)
	    # carry out averaging
	    tmp=crop_man_all$plant[average_i,average_j] ; tmp[which(tmp == -9999)]=NA
	    plant=mean(tmp, na.rm=TRUE)
	    tmp=crop_man_all$harvest[average_i,average_j] ; tmp[which(tmp == -9999)]=NA
	    harvest=mean(tmp, na.rm=TRUE)
#	    tmp=crop_man_all$plant_range[average_i,average_j] ; tmp[which(tmp == -9999)]=NA
#	    plant_range=mean(tmp, na.rm=TRUE)
#	    tmp=crop_man_all$harvest_range[average_i,average_j] ; tmp[which(tmp == -9999)]=NA
#	    harvest_range=mean(tmp, na.rm=TRUE)
	    # error checking
	    if (is.na(plant) | plant == 0 | harvest == 0 | is.na(harvest)) {radius=radius+1 ; answer=NA} else {answer=0}
	}
	print(paste("NOTE Sowing and Harvest dates averaged over a pixel radius (i.e. centre + radius) of ",radius," points",sep=""))
	# pass the information back
	return(list(plant=plant,plant_range=plant_range,harvest=harvest,harvest_range=harvest_range))

} # end of function
