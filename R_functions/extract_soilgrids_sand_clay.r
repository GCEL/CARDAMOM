
###
## Function to extract location specific information on soil texture from the gridded SoilGrids database
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

extract_soilgrid_sand_clay<- function(spatial_type,resolution,grid_type,latlon_in,sand_clay_all) {

  # Update the user
	print(paste("Sand/clay data extracted for current location ",Sys.time(),sep=""))

	# convert input data long to conform to what we need
	check1=which(sand_clay_all$long > 180) ; if (length(check1) > 0) { sand_clay_all$long[check1]=sand_clay_all$long[check1]-360 }

	# find the nearest location
	output=closest2d(1,sand_clay_all$lat,sand_clay_all$long,latlon_in[1],latlon_in[2],2)
	j1=unlist(output)[2];i1=unlist(output)[1]

	# return long to 0-360
	if (length(check1) > 0) { sand_clay_all$long[check1]=sand_clay_all$long[check1]+360 }
  # If resolution has been provides as single value then adjust this here
  if (length(resolution) == 1 & spatial_type == "grid") {tmp_res = resolution * c(1,1)} else {tmp_res = resolution}

	# work out number of pixels to average over
#	print("NOTE all Csom values are a minimum average of 9 pixels (i.e centre+1 )")
	if (spatial_type == "grid") {
      # Extract the product resolution assumed to be degrees (x,y)
      product_res = c(abs(sand_clay_all$long[2,1]-sand_clay_all$long[1,1]),abs(sand_clay_all$lat[1,2]-sand_clay_all$lat[1,1]))
	    if (grid_type == "wgs84") {
          # radius is floor of the ratio of the product vs analysis ratio
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

	answer = NA
	while (is.na(answer) == TRUE) {
	    # work out average areas
	    average_i = (i1-radius[1]):(i1+radius[1]) ; average_j = (j1-radius[2]):(j1+radius[2])
	    average_i = max(1,(i1-radius[1])):min(dim(sand_clay_all$top_sand)[1],(i1+radius[1]))
      average_j = max(1,(j1-radius[2])):min(dim(sand_clay_all$top_sand)[2],(j1+radius[2]))
	    # carry out averaging
	    tmp1 = sand_clay_all$top_sand[average_i,average_j] ; tmp1[which(tmp1 == -9999)] = NA
	    tmp2 = sand_clay_all$bot_sand[average_i,average_j] ; tmp2[which(tmp2 == -9999)] = NA
	    tmp3 = sand_clay_all$top_clay[average_i,average_j] ; tmp3[which(tmp3 == -9999)] = NA
	    tmp4 = sand_clay_all$bot_clay[average_i,average_j] ; tmp4[which(tmp4 == -9999)] = NA
	    top_sand = mean(tmp1, na.rm=TRUE) ; bot_sand=mean(tmp2, na.rm=TRUE)
	    top_clay = mean(tmp3, na.rm=TRUE) ; bot_clay=mean(tmp4, na.rm=TRUE)
	    # error checking
	    if (is.na(top_sand) | top_sand == 0) {radius = radius+1 ; answer = NA} else {answer = 0}
	}
  # Inform the user
	print(paste("NOTE sand/clay averaged over a pixel radius (i.e. centre + radius) of ",radius," points",sep=""))

	# just to check because when averaging sometimes the sand / clay combinations can be > 100 %
	# 94 % chosesn as this is the highest total % found in the HWSD dataset
	if ((top_sand+top_clay) > 94) {
	    tmp1 = top_sand / (top_sand + top_clay + 6) # 6 % is implicit in the 94 % max value for silt / gravel
	    tmp2 = top_clay / (top_sand + top_clay + 6) # 6 % is implicit in the 94 % max value for silt / gravel
	    top_sand = tmp1*100 ; top_clay = tmp2*100
	}
	if ((bot_sand+bot_clay) > 94) {
	    tmp1 = bot_sand / (bot_sand + bot_clay + 6) # 6 % is implicit in the 94 % max value for silt / gravel
	    tmp2 = bot_clay / (bot_sand + bot_clay + 6) # 6 % is implicit in the 94 % max value for silt / gravel
	    top_sand = tmp1*100 ; top_clay = tmp2*100
	}

	# pass the information back
	return(list(top_sand=top_sand,bot_sand=bot_sand,top_clay=top_clay,bot_clay=bot_clay))

} # end function extract_soilgrid_sand_clay

## Use byte compile
extract_soilgrid_sand_clay<-cmpfun(extract_soilgrid_sand_clay)
