
###
## Function to load met data from global field ECMWF data
## subsequently extracted in extract_met_drivers.txt
###

load_hwsd_sand_clay_fields_for_extraction<-function(latlon_in,sand_clay_source) {

    if (sand_clay_source == "SoilGrids") {

        # Read in the data for both the sand and clay
        # Sand
        top_sand = raster(paste(path_to_sand_clay,"clay_percent_mean_0to30cm.tif", sep=""))
        bot_sand = raster(paste(path_to_sand_clay,"clay_percent_mean_30to100cm.tif", sep=""))
        # Clay
        top_clay = raster(paste(path_to_sand_clay,"sand_percent_mean_0to30cm.tif", sep=""))
        bot_clay = raster(paste(path_to_sand_clay,"sand_percent_mean_30to100cm.tif", sep=""))

        # Extract dimension information for the grid.
        # Note 1) the axis switching between raster and actual array
        #      2) we only do this once as the lat / long grid for both maps is identical
        xdim = dim(top_sand)[2] ; ydim = dim(top_sand)[1]
        # extract the lat / long information needed
        long = coordinates(top_sand)[,1] ; lat = coordinates(top_sand)[,2]
        # restructure into correct orientation
        long = array(long, dim=c(xdim,ydim))
        lat = array(lat, dim=c(xdim,ydim))

        # Break out from the rasters into arrays which we can manipulate
        # Sand
        top_sand = array(as.vector(unlist(top_sand)), dim=c(xdim,ydim))
        bot_sand = array(as.vector(unlist(bot_sand)), dim=c(xdim,ydim))
         # Clay
        top_clay = array(as.vector(unlist(top_clay)), dim=c(xdim,ydim))
        bot_clay = array(as.vector(unlist(bot_clay)), dim=c(xdim,ydim))

        # output variables
        return(list(top_sand = top_sand, top_clay = top_clay, bot_sand = bot_sand, bot_clay = bot_clay,lat = lat,long = long))

    } else if(Csom_source == "SoilGrids2") {

          # file names - note that 0-30cm and 30-100cm in different files
          sand_file1=paste(path_to_sand_clay,"/sand_000-030cm_mean(.)*.nc",sep="")
        	sand1=nc_open(sand_file_1)
          sand_file2=paste(path_to_sand_clay,"/sand_030-100cm_mean(.)*.nc",sep="")
        	sand2=nc_open(sand_file_2)
          clay_file1=paste(path_to_sand_clay,"/clay_000-030cm_mean(.)*.nc",sep="")
        	clay1=nc_open(clay_file_1)
          clay_file2=paste(path_to_sand_clay,"/clay_030-100cm_mean(.)*.nc",sep="")
        	clay2=nc_open(clay_file_2)

        	# extract location variables
        	lat=ncvar_get(sand1, "latitude") ; long=ncvar_get(sand1, "longitude")
        	# read the texxture data from the relevant files
          top_sand=ncvar_get(sand1, "sand")
          bot_sand=ncvar_get(sand2, "sand")
          top_clay=ncvar_get(clay1, "clay")
          bot_clay=ncvar_get(clay2, "clay")

        	if (length(dim(latlon_in)) > 1) {
        	    max_lat=max(latlon_in[,1])+0.5 ; max_long=max(latlon_in[,2])+0.5
        	    min_lat=min(latlon_in[,1])-0.5 ; min_long=min(latlon_in[,2])-0.5
        	} else {
        	    max_lat=max(latlon_in[1])+0.5 ; max_long=max(latlon_in[2])+0.5
        	    min_lat=min(latlon_in[1])-0.5 ; min_long=min(latlon_in[2])-0.5
        	}
        	keep_lat=which(lat[1,] > min_lat & lat[1,] < max_lat)
        	keep_long=which(long[,1] > min_long & long[,1] < max_long)
          top_sand=top_sand[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
          bot_sand=bot_sand[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
          top_clay=top_clay[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
          bot_clay=bot_clay[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]

          lat=lat[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
        	long=long[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]

        	# close files after use
        	nc_close(sand1)
          nc_close(sand2)
        	nc_close(clay1)
          nc_close(clay2)
          return(list(top_sand = top_sand, top_clay = top_clay, bot_sand = bot_sand, bot_clay = bot_clay, lat = lat,long = long))

    } else if (sand_clay_source == "HWSD") {

	# let the user know this might take some time
	print("Loading processed HWSD sand clay fields for subsequent sub-setting ...")

	# open processed modis files
	input_file_1=paste(path_to_sand_clay,"/HWSD_sand_clay_with_lat_long.nc",sep="")
	data1=nc_open(input_file_1)

	# extract location variables
	lat=ncvar_get(data1, "lat") ; long=ncvar_get(data1, "long")
	# read the modis lai drivers
	hwsd_top_sand=ncvar_get(data1, "HWSD_top_sand") ; hwsd_top_clay=ncvar_get(data1, "HWSD_top_clay")
	hwsd_bot_sand=ncvar_get(data1, "HWSD_bot_sand") ; hwsd_bot_clay=ncvar_get(data1, "HWSD_bot_clay")

	if (length(dim(latlon_in)) > 1) {
	    max_lat=max(latlon_in[,1])+0.5 ; max_long=max(latlon_in[,2])+0.5
	    min_lat=min(latlon_in[,1])-0.5 ; min_long=min(latlon_in[,2])-0.5
	} else {
	    max_lat=max(latlon_in[1])+0.5 ; max_long=max(latlon_in[2])+0.5
	    min_lat=min(latlon_in[1])-0.5 ; min_long=min(latlon_in[2])-0.5
	}
	keep_lat=which(lat[1,] > min_lat & lat[1,] < max_lat)
	keep_long=which(long[,1] > min_long & long[,1] < max_long)
	hwsd_top_sand=hwsd_top_sand[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
	hwsd_bot_sand=hwsd_bot_sand[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
	hwsd_top_clay=hwsd_top_clay[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
	hwsd_bot_clay=hwsd_bot_clay[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
	lat=lat[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
	long=long[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
	# close files after use
	nc_close(data1)

	# clean
	rm(keep_lat,keep_long,max_lat,max_long,min_lat,min_long) ; gc(reset=TRUE,verbose=FALSE)

	# output variables
	return(list(top_sand=hwsd_top_sand,top_clay=hwsd_top_clay,bot_sand=hwsd_bot_sand,bot_clay=hwsd_bot_clay,lat=lat,long=long))

    } else {
	# output variables
	return(list(top_sand=40,top_clay=15,bot_sand=40,bot_clay=15,lat=-9999,long=-9999))
    }

} # function end
