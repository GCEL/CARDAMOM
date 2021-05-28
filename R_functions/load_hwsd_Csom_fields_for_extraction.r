
###
## Function to load gridded dataset of soil prior information from HWSD
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_hwsd_Csom_fields_for_extraction<-function(latlon_in,Csom_source) {

    if (Csom_source == "SoilGrids") {

        # this is a very bespoke modification so leave it here to avoid getting lost
        Csom = raster(paste(path_to_Csom,"Csom_gCm2_mean_0to1m.tif", sep=""))
        Csom_unc = raster(paste(path_to_Csom,"Csom_gCm2_sd_0to1m.tif", sep=""))
        # extract dimension information for the grid, note the axis switching between raster and actual array
        xdim = dim(Csom)[2] ; ydim = dim(Csom)[1]
        # extract the lat / long information needed
        long = coordinates(Csom)[,1] ; lat = coordinates(Csom)[,2]
        # restructure into correct orientation
        long = array(long, dim=c(xdim,ydim))
        lat = array(lat, dim=c(xdim,ydim))
        # break out from the rasters into arrays which we can manipulate
        Csom = array(as.vector(unlist(Csom)), dim=c(xdim,ydim))
        Csom_unc = array(as.vector(unlist(Csom_unc)), dim=c(xdim,ydim))


        return(list(Csom = Csom, Csom_unc = Csom_unc, lat = lat,long = long))

    } else if(Csom_source == "SoilGrids2") {

        # file names - note that 0-30cm and 30-100cm in different files
        soc_file1=paste(path_to_Csom,"/SOC_000-030cm_mean(.)*.nc",sep="")
      	soc1=nc_open(soc_file_1)
        soc_file2=paste(path_to_Csom,"/SOC_030-100cm_mean(.)*.nc",sep="")
      	soc2=nc_open(soc_file_2)

      	# extract location variables
      	lat=ncvar_get(soc1, "latitude") ; long=ncvar_get(soc1, "longitude")
      	# read the SOC from the two files
        Csom1=ncvar_get(soc1, "SOC")
        Csom2=ncvar_get(soc2, "SOC")

      	if (length(dim(latlon_in)) > 1) {
      	    max_lat=max(latlon_in[,1])+0.5 ; max_long=max(latlon_in[,2])+0.5
      	    min_lat=min(latlon_in[,1])-0.5 ; min_long=min(latlon_in[,2])-0.5
      	} else {
      	    max_lat=max(latlon_in[1])+0.5 ; max_long=max(latlon_in[2])+0.5
      	    min_lat=min(latlon_in[1])-0.5 ; min_long=min(latlon_in[2])-0.5
      	}
      	keep_lat=which(lat[1,] > min_lat & lat[1,] < max_lat)
      	keep_long=which(long[,1] > min_long & long[,1] < max_long)
        Csom1=Csom1[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
        Csom2=Csom2[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
        Csom = Csom1 + Csom2
        lat=lat[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
      	long=long[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]

        # repeat for uncertainty (uncertainty added)
        # file names - note that 0-30cm and 30-100cm in different files
        unc_file1=paste(path_to_Csom,"/SOC_000-030cm_uncertainty(.)*.nc",sep="")
      	unc1=nc_open(unc_file_1)
        unc_file2=paste(path_to_Csom,"/SOC_030-100cm_uncertainty(.)*.nc",sep="")
      	unc2=nc_open(unc_file_2)

      	# read the SOC from the two files
        Cunc1=ncvar_get(unc1, "SOC_uncertainty")
        Cunc2=ncvar_get(unc2, "SOC_uncertainty")
        Cunc1=Cunc1[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
        Cunc2=Cunc2[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
        Csom_unc = Cunc1 + Cunc2
      	# close files after use
      	nc_close(soc1)
        nc_close(soc2)
      	nc_close(unc1)
        nc_close(unc2)

        return(list(Csom = Csom, Csom_unc = Csom_unc, lat = lat,long = long))

    } else if (Csom_source == "HWSD") {

	# let the user know this might take some time
	print("Loading processed HWSD Csom fields for subsequent sub-setting ...")

	# open processed modis files
	input_file_1=paste(path_to_Csom,"/HWSD_Csom_with_lat_long.nc",sep="")
	data1=nc_open(input_file_1)

	# extract location variables
	lat=ncvar_get(data1, "lat") ; long=ncvar_get(data1, "long")
	# read the modis lai drivers
	hwsd_Csom=ncvar_get(data1, "HWSD_Csom")

	if (length(dim(latlon_in)) > 1) {
	    max_lat=max(latlon_in[,1])+0.5 ; max_long=max(latlon_in[,2])+0.5
	    min_lat=min(latlon_in[,1])-0.5 ; min_long=min(latlon_in[,2])-0.5
	} else {
	    max_lat=max(latlon_in[1])+0.5 ; max_long=max(latlon_in[2])+0.5
	    min_lat=min(latlon_in[1])-0.5 ; min_long=min(latlon_in[2])-0.5
	}
	keep_lat=which(lat[1,] > min_lat & lat[1,] < max_lat)
	keep_long=which(long[,1] > min_long & long[,1] < max_long)
	hwsd_Csom=hwsd_Csom[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
	lat=lat[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
	long=long[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
	# close files after use
	nc_close(data1)

	# clean
	rm(keep_lat,keep_long,max_lat,max_long,min_lat,min_long) ; gc(reset=TRUE,verbose=FALSE)

	# output variables
	return(list(Csom=hwsd_Csom,Csom_unc = -9999,lat=lat,long=long))

} else {
	# output variables
	return(list(Csom=-9999, Csom_unc = -9999, lat=-9999,long=-9999))

    }

} # function end
