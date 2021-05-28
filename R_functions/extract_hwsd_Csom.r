
####
## Function to extract the HWSD soil carbon estimate for initial conditions
####

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

extract_hwsd_Csom<- function(spatial_type,resolution,grid_type,latlon_in,Csom_all) {

	# convert input data long to conform to what we need
	check1=which(Csom_all$long > 180) ; if (length(check1) > 0) { Csom_all$long[check1]=Csom_all$long[check1]-360 }

	# find the nearest location
	output=closest2d(1,Csom_all$lat,Csom_all$long,latlon_in[1],latlon_in[2],2)
	j1=unlist(output)[2];i1=unlist(output)[1]
	print(paste("Csom data extracted for current location ",Sys.time(),sep=""))

	# return long to 0-360
	if (length(check1) > 0) { Csom_all$long[check1]=Csom_all$long[check1]+360 }

	# work out number of pixels to average over
	if (spatial_type == "grid") {
      # resolution of the product
      product_res = abs(Csom_all$lat[1,2]-Csom_all$lat[1,1])+abs(Csom_all$long[2,1]-Csom_all$long[1,1])
      product_res = product_res * 0.5 # NOTE: averaging needed for line above
	    if (grid_type == "wgs84") {
          # radius is ceiling of the ratio of the product vs analysis ratio
#          radius = round(resolution / product_res, digits=0)
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
	    radius=0
	}

	answer=NA
	while (is.na(answer) == TRUE) {
	    # work out average areas
	    average_i=(i1-radius):(i1+radius) ; average_j=(j1-radius):(j1+radius)
	    average_i=max(1,(i1-radius)):min(dim(Csom_all$Csom)[1],(i1+radius)) ; average_j=max(1,(j1-radius)):min(dim(Csom_all$Csom)[2],(j1+radius))
	    # carry out averaging
	    tmp=Csom_all$Csom[average_i,average_j] ; tmp[which(tmp == -9999)]=NA
	    Csom=mean(tmp, na.rm=TRUE)
	    # error checking
	    if (is.na(Csom) | Csom == 0) {radius=radius+1 ; answer=NA} else {answer=0}

	}
	print(paste("NOTE Csom averaged over a pixel radius (i.e. centre + radius) of ",radius," points",sep=""))
	# pass the information back
	return(Csom)

} # end of function
