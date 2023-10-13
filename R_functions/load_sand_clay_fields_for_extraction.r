
###
## Function to load soil texture information from global gridded HWSD
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_sand_clay_fields_for_extraction<-function(latlon_in,sand_clay_source,cardamom_ext,spatial_type) {

    if (sand_clay_source == "SoilGrids" | sand_clay_source == "SoilGrids_v2") {

        print("Loading sand / clay fractions from SoilGrids")

        # Read in the data for both the sand and clay
        # Sand
        top_sand = rast(paste(path_to_sand_clay,"sand_percent_mean_0to30cm.tif", sep=""))
        bot_sand = rast(paste(path_to_sand_clay,"sand_percent_mean_30to100cm.tif", sep=""))
        # Clay
        top_clay = rast(paste(path_to_sand_clay,"clay_percent_mean_0to30cm.tif", sep=""))
        bot_clay = rast(paste(path_to_sand_clay,"clay_percent_mean_30to100cm.tif", sep=""))

        # Create raster with the target crs
        target = rast(crs = ("+init=epsg:4326"), ext = ext(top_sand), resolution = res(top_sand))
        # Check whether the target and actual analyses have the same CRS
        if (compareGeom(top_sand,target) == FALSE) {
            # Resample to correct grid
            top_sand = resample(top_sand, target, method="ngb") ; gc() 
            bot_sand = resample(bot_sand, target, method="ngb") ; gc()
            top_clay = resample(top_clay, target, method="ngb") ; gc()
            bot_clay = resample(bot_clay, target, method="ngb") ; gc()
        }
        # Extend the extent of the overall grid to the analysis domain
        top_sand = extend(top_sand,cardamom_ext) ; bot_sand = extend(bot_sand,cardamom_ext)
        top_clay = extend(top_clay,cardamom_ext) ; bot_clay = extend(bot_clay,cardamom_ext)
        # Trim the extent of the overall grid to the analysis domain
        top_sand = crop(top_sand,cardamom_ext) ; bot_sand = crop(bot_sand,cardamom_ext)
        top_clay = crop(top_clay,cardamom_ext) ; bot_clay = crop(bot_clay,cardamom_ext)
        # Adjust spatial resolution of the datasets, this occurs in all cases
        if (res(top_sand)[1] != res(cardamom_ext)[1] | res(top_sand)[2] != res(cardamom_ext)[2]) {

            # Create raster with the target resolution
            target = rast(crs = crs(cardamom_ext), ext = ext(cardamom_ext), resolution = res(cardamom_ext))

            # Resample to correct grid
            top_sand = resample(top_sand, target, method="bilinear") ; gc() 
            bot_sand = resample(bot_sand, target, method="bilinear") ; gc() 
            top_clay = resample(top_clay, target, method="bilinear") ; gc() 
            bot_clay = resample(bot_clay, target, method="bilinear") ; gc() 

        } # Aggrgeate to resolution

        # Extract dimension information for the grid.
        # Note 1) the axis switching between raster and actual array
        #      2) we only do this once as the lat / long grid for both maps is identical
        xdim = dim(top_sand)[2] ; ydim = dim(top_sand)[1]
        # extract the lat / long information needed
        long = crds(top_sand,df=TRUE, na.rm=FALSE)
        lat  = long$y ; long = long$x        # restructure into correct orientation
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

     }  else if(Csom_source == "SoilGrids_v2_stratified") {
        # This has been added by DTM for the DARE-UK project
        # file names - note that 0-30cm and 30-100cm in different files
        top_sand_file=paste(path_to_Csom,"/sand_000-030cm_mean*.nc",sep="")
      	top_sand_nc=nc_open(Sys.glob(top_sand_file))
        bottom_sand_file=paste(path_to_Csom,"/sand_030-100cm_mean*.nc",sep="")
      	bottom_sand_nc=nc_open(Sys.glob(bottom_sand_file))
        
        top_clay_file=paste(path_to_Csom,"/clay_000-030cm_mean*.nc",sep="")
      	top_clay_nc=nc_open(Sys.glob(top_sand_file))
        bottom_clay_file=paste(path_to_Csom,"/clay_030-100cm_mean*.nc",sep="")
      	bottom_clay_nc=nc_open(Sys.glob(bottom_clay_file))
        
      	# extract location variables
      	lat=ncvar_get(top_sand_nc, "latitude") ; long=ncvar_get(top_sand_nc, "longitude")
        res = abs(diff(lat[1:2]))
        max_lat = max(latlon_in[,1])+res/2. ; max_long = max(latlon_in[,2])+res/2.
        min_lat = min(latlon_in[,1])-res/2. ; min_long = min(latlon_in[,2])-res/2.
        keep_lat = which((lat > min_lat) & (lat < max_lat))
        keep_long = which((long > min_long) & (long < max_long))

        lat = lat[keep_lat]
        long = long[keep_long]
        xdim = length(long) ; ydim = length(lat)
        long = array(long, dim=c(xdim,ydim))
        lat = t(array(lat, dim=c(ydim,xdim)))

        # read the sand from the two files
        top_sand=ncvar_get(top_sand_nc, "sand")
        bottom_sand=ncvar_get(bottom_sand_nc, "sand")
	nc_close(top_sand_nc)
	nc_close(bottom_sand_nc)
	    
      	top_sand=top_sand[keep_long,keep_lat]
        bottom_sand=bottom_sand[keep_long,keep_lat]

	# read the clay from the two files
        top_clay=ncvar_get(top_clay_nc, "clay")
        bottom_clay=ncvar_get(bottom_clay_nc, "clay")
	nc_close(top_clay_nc)
	nc_close(bottom_clay_nc)
	    
      	top_clay=top_clay[keep_long,keep_lat]
        bottom_clay=bottom_clay[keep_long,keep_lat]
        
	
	return(list(top_sand = top_sand, top_clay = top_clay, bot_sand = bot_sand, bot_clay = bot_clay, lat = lat,long = long))

   } else if (sand_clay_source == "HWSD") {

        # let the user know this might take some time
        print("Loading processed HWSD sand clay fields for subsequent sub-setting ...")

        # open processed modis files
        input_file_1 = paste(path_to_sand_clay,"/HWSD_sand_clay_with_lat_long.nc",sep="")
        data1 = nc_open(input_file_1)

        # extract location variables
        lat = ncvar_get(data1, "lat") ; long = ncvar_get(data1, "long")
        # read the HWSD datasets
	      top_sand = ncvar_get(data1, "HWSD_top_sand") ; top_clay = ncvar_get(data1, "HWSD_top_clay")
        bot_sand = ncvar_get(data1, "HWSD_bot_sand") ; bot_clay = ncvar_get(data1, "HWSD_bot_clay")
        nc_close(data1)

        # Convert to a raster, assuming standad WGS84 grid
        top_sand = data.frame(x = as.vector(long), y = as.vector(lat), z = as.vector(top_sand))
        top_sand = rast(top_sand, crs = ("+init=epsg:4326"), type="xyz")
        bot_sand = data.frame(x = as.vector(long), y = as.vector(lat), z = as.vector(bot_sand))
        bot_sand = rast(bot_sand, crs = ("+init=epsg:4326"), type="xyz")
        top_clay = data.frame(x = as.vector(long), y = as.vector(lat), z = as.vector(top_clay))
        top_clay = rast(top_clay, crs = ("+init=epsg:4326"), type="xyz")
        bot_clay = data.frame(x = as.vector(long), y = as.vector(lat), z = as.vector(bot_clay))
        bot_clay = rast(bot_clay, crs = ("+init=epsg:4326"), type="xyz")

        # Create raster with the target crs
        target = rast(crs = ("+init=epsg:4326"), ext = ext(top_sand), resolution = res(top_sand))
        # Check whether the target and actual analyses have the same CRS
        if (compareGeom(top_sand,target) == FALSE) {
            # Resample to correct grid
            top_sand = resample(top_sand, target, method="ngb") ; gc() 
            bot_sand = resample(bot_sand, target, method="ngb") ; gc() 
            top_clay = resample(top_clay, target, method="ngb") ; gc() 
            bot_clay = resample(bot_clay, target, method="ngb") ; gc() 
        }
        # Extend the extent of the overall grid to the analysis domain
        top_sand = extend(top_sand,cardamom_ext) ; bot_sand = extend(bot_sand,cardamom_ext)
        top_clay = extend(top_clay,cardamom_ext) ; bot_clay = extend(bot_clay,cardamom_ext)
        # Trim the extent of the overall grid to the analysis domain
        top_sand = crop(top_sand,cardamom_ext) ; bot_sand = crop(bot_sand,cardamom_ext)
        top_clay = crop(top_clay,cardamom_ext) ; bot_clay = crop(bot_clay,cardamom_ext)

        # Adjust spatial resolution of the datasets, this occurs in all cases
        if (res(top_sand)[1] != res(cardamom_ext)[1] | res(top_sand)[2] != res(cardamom_ext)[2]) {

            # Create raster with the target resolution
            target = rast(crs = crs(cardamom_ext), ext = ext(cardamom_ext), resolution = res(cardamom_ext))

            # Resample to correct grid
            top_sand = resample(top_sand, target, method="bilinear") ; gc() 
            bot_sand = resample(bot_sand, target, method="bilinear") ; gc()
            top_clay = resample(top_clay, target, method="bilinear") ; gc()
            bot_clay = resample(bot_clay, target, method="bilinear") ; gc()

        } # Aggrgeate to resolution

        # Extract dimension information for the grid.
        # Note 1) the axis switching between raster and actual array
        #      2) we only do this once as the lat / long grid for both maps is identical
        xdim = dim(top_sand)[2] ; ydim = dim(top_sand)[1]
        # extract the lat / long information needed
        long = crds(top_sand,df=TRUE, na.rm=FALSE)
        lat  = long$y ; long = long$x
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
        return(list(top_sand=top_sand,top_clay=top_clay,bot_sand=bot_sand,bot_clay=bot_clay,lat=lat,long=long))

    } else {
        # output variables
        return(list(top_sand=40,top_clay=15,bot_sand=40,bot_clay=15,lat=-9999,long=-9999))
    }

} # function end load_sand_clay_fields_for_extraction

## Use byte compile
load_sand_clay_fields_for_extraction<-cmpfun(load_sand_clay_fields_for_extraction)
