
###
## Function to load gridded dataset of soil prior information from HWSD
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_Csom_fields_for_extraction<-function(latlon_in,Csom_source,cardamom_ext,spatial_type) {

    if (Csom_source == "SoilGrids") {

        # let the user know this might take some time
        print("Loading processed SoilGrids Csom fields for subsequent sub-setting ...")

        # This is a very bespoke modification so leave it here to avoid getting lost
        Csom = rast(paste(path_to_Csom,"Csom_gCm2_mean_0to1m.tif", sep=""))
        Csom_unc = rast(paste(path_to_Csom,"Csom_gCm2_sd_0to1m.tif", sep=""))

        # Create raster with the target crs
        target = rast(crs = ("+init=epsg:4326"), ext = ext(Csom), resolution = res(Csom))
        # Check whether the target and actual analyses have the same CRS
        if (compareGeom(Csom,target) == FALSE) {
            # Resample to correct grid
            Csom = resample(Csom, target, method="ngb") ; gc() 
            Csom_unc = resample(Csom_unc, target, method="ngb") ; gc() 
        }
        # Extend the extent of the overall grid to the analysis domain
        Csom = extend(Csom,cardamom_ext) ; Csom_unc = extend(Csom_unc,cardamom_ext)
        # Trim the extent of the overall grid to the analysis domain
        Csom = crop(Csom,cardamom_ext) ; Csom_unc = crop(Csom_unc,cardamom_ext)
        # Adjust spatial resolution of the datasets, this occurs in all cases
        if (res(Csom)[1] != res(cardamom_ext)[1] | res(Csom)[2] != res(cardamom_ext)[2]) {

            # Create raster with the target resolution
            target = rast(crs = crs(cardamom_ext), ext = ext(cardamom_ext), resolution = res(cardamom_ext))

            # Resample to correct grid
            Csom = resample(Csom, target, method="bilinear") ; gc() 
            Csom_unc = resample(Csom_unc, target, method="bilinear") ; gc() 

        } # Aggrgeate to resolution

        # extract dimension information for the grid, note the axis switching between raster and actual array
        xdim = dim(Csom)[2] ; ydim = dim(Csom)[1]
        # extract the lat / long information needed
        long = crds(Csom,df=TRUE, na.rm=FALSE)
        lat  = long$y ; long = long$x
        # restructure into correct orientation
        long = array(long, dim=c(xdim,ydim))
        lat = array(lat, dim=c(xdim,ydim))
        # break out from the rasters into arrays which we can manipulate
        Csom = array(as.vector(unlist(Csom)), dim=c(xdim,ydim))
        Csom_unc = array(as.vector(unlist(Csom_unc)), dim=c(xdim,ydim))

        # Assume in all cases than a zero prior value should be classed as missing data
        Csom[Csom < 1] = NA ; Csom_unc[is.na(Csom)] = NA

        return(list(Csom = Csom, Csom_unc = Csom_unc, lat = lat, long = long))

    } else if (Csom_source == "SoilGrids_v2") {

        # let the user know this might take some time
        print("Loading processed SoilGrids_v2 Csom fields for subsequent sub-setting ...")

        # This is a very bespoke modification so leave it here to avoid getting lost
        Csom = rast(paste(path_to_Csom,"Csom_gCm2_mean_0to100cm.tif", sep=""))
        Csom_unc = rast(paste(path_to_Csom,"Csom_gCm2_uncertainty_0to100cm.tif", sep=""))

        # Create raster with the target crs
        target = rast(crs = ("+init=epsg:4326"), ext = ext(Csom), resolution = res(Csom))
        # Check whether the target and actual analyses have the same CRS
        if (compareGeom(Csom,target) == FALSE) {
            # Resample to correct grid
            Csom = resample(Csom, target, method="ngb") ; gc() 
            Csom_unc = resample(Csom_unc, target, method="ngb") ; gc() 
        }
        # Extend the extent of the overall grid to the analysis domain
        Csom = extend(Csom,cardamom_ext) ; Csom_unc = extend(Csom_unc,cardamom_ext)
        # Trim the extent of the overall grid to the analysis domain
        Csom = crop(Csom,cardamom_ext) ; Csom_unc = crop(Csom_unc,cardamom_ext)
        # Adjust spatial resolution of the datasets, this occurs in all cases
        if (res(Csom)[1] != res(cardamom_ext)[1] | res(Csom)[2] != res(cardamom_ext)[2]) {

            # Create raster with the target resolution
            target = rast(crs = crs(cardamom_ext), ext = ext(cardamom_ext), resolution = res(cardamom_ext))

            # Resample to correct grid
            Csom = resample(Csom, target, method="bilinear") ; gc() 
            Csom_unc = resample(Csom_unc, target, method="bilinear") ; gc() 

        } # Aggrgeate to resolution

        # extract dimension information for the grid, note the axis switching between raster and actual array
        xdim = dim(Csom)[2] ; ydim = dim(Csom)[1]
        # extract the lat / long information needed
        long = crds(Csom,df=TRUE, na.rm=FALSE)
        lat  = long$y ; long = long$x
        # restructure into correct orientation
        long = array(long, dim=c(xdim,ydim))
        lat = array(lat, dim=c(xdim,ydim))
        # break out from the rasters into arrays which we can manipulate
        Csom = array(as.vector(unlist(Csom)), dim=c(xdim,ydim))
        Csom_unc = array(as.vector(unlist(Csom_unc)), dim=c(xdim,ydim))

        # Assume in all cases than a zero prior value should be classed as missing data
        Csom[Csom < 1] = NA ; Csom_unc[is.na(Csom)] = NA

        return(list(Csom = Csom, Csom_unc = Csom_unc, lat = lat, long = long))

      }  else if(Csom_source == "SoilGrids_v2_stratified") {
        # This has been added by DTM for the DARE-UK project
        # file names - note that 0-30cm and 30-100cm in different files
        soc_file1=paste(path_to_Csom,"/SOC_000-030cm_mean*.nc",sep="")
      	soc1=nc_open(Sys.glob(soc_file1))
        soc_file2=paste(path_to_Csom,"/SOC_030-100cm_mean*.nc",sep="")
      	soc2=nc_open(Sys.glob(soc_file2))

      	# extract location variables
      	lat=ncvar_get(soc1, "latitude") ; long=ncvar_get(soc1, "longitude")
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

        # read the SOC from the two files
        Csom1=ncvar_get(soc1, "SOC")
        Csom2=ncvar_get(soc2, "SOC")

      	Csom1=Csom1[keep_long,keep_lat]
        Csom2=Csom2[keep_long,keep_lat]
        Csom = array(Csom1 + Csom2, dim=c(xdim,ydim))

        # repeat for uncertainty (uncertainty added)
        # file names - note that 0-30cm and 30-100cm in different files
        unc_file1=paste(path_to_Csom,"/SOC_000-030cm_uncertainty*.nc",sep="")
      	unc1=nc_open(Sys.glob(unc_file1))
        unc_file2=paste(path_to_Csom,"/SOC_030-100cm_uncertainty*.nc",sep="")
      	unc2=nc_open(Sys.glob(unc_file2))

      	# read the SOC from the two files
        Cunc1=ncvar_get(unc1, "SOC_uncertainty")
        Cunc2=ncvar_get(unc2, "SOC_uncertainty")

        Cunc1=Cunc1[keep_long,keep_lat]
        Cunc2=Cunc2[keep_long,keep_lat]

        Csom_unc = array(Cunc1 + Cunc2, dim=c(xdim,ydim))
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
        # read the HWSD soil C prior
        Csom=ncvar_get(data1, "HWSD_Csom")

        # Convert to a raster, assuming standad WGS84 grid
        Csom = data.frame(x = as.vector(long), y = as.vector(lat), z = as.vector(Csom))
        Csom = rast(Csom, crs = ("+init=epsg:4326"), type="xyz")

        # Create raster with the target crs (technically this bit is not required)
        target = rast(crs = ("+init=epsg:4326"), ext = ext(Csom), resolution = res(Csom))
        # Check whether the target and actual analyses have the same CRS
        if (compareGeom(Csom,target) == FALSE) {
            # Resample to correct grid
            Csom = resample(Csom, target, method="ngb") ; gc() 
        }
        # Extend the extent of the overall grid to the analysis domain
        Csom = extend(Csom,cardamom_ext)
        # Trim the extent of the overall grid to the analysis domain
        Csom = crop(Csom,cardamom_ext)
        # Adjust spatial resolution of the datasets, this occurs in all cases
        if (res(Csom)[1] != res(cardamom_ext)[1] | res(Csom)[2] != res(cardamom_ext)[2]) {

            # Create raster with the target resolution
            target = rast(crs = crs(cardamom_ext), ext = ext(cardamom_ext), resolution = res(cardamom_ext))
            # Resample to correct grid
            Csom = resample(Csom, target, method="bilinear") ; gc() 

        } # Aggrgeate to resolution

        # extract dimension information for the grid, note the axis switching between raster and actual array
        xdim = dim(Csom)[2] ; ydim = dim(Csom)[1]
        # extract the lat / long information needed
        long = crds(Csom,df=TRUE, na.rm=FALSE)
        lat  = long$y ; long = long$x
        # restructure into correct orientation
        long = array(long, dim=c(xdim,ydim))
        lat = array(lat, dim=c(xdim,ydim))
        # break out from the rasters into arrays which we can manipulate
        Csom = array(as.vector(unlist(Csom)), dim=c(xdim,ydim))

        # Assume in all cases than a zero prior value should be classed as missing data
        Csom[Csom < 1] = NA
        # see papers assessing uncertainty of HWSD, ~47 %
        Csom_unc = array(Csom * 0.47, dim=c(xdim,ydim))
        # With a minium bound assumption
        Csom_unc[Csom_unc < 100] = 100
        # Ensure consistency for missing values
        Csom_unc[Csom == -9999 | is.na(Csom) == TRUE] = -9999

        # Return the loaded dataset
        return(list(Csom = Csom, Csom_unc = Csom_unc, lat = lat, long = long))

    } else if (Csom_source == "NCSCD") {

        # let the user know this might take some time
        print("Loading processed NCSCD Csom fields for subsequent sub-setting ...")

        # open processed modis files
        input_file_1 = paste(path_to_Csom,"/NCSCD_Csom_0-1m_with_lat_long.nc",sep="")
        data1 = nc_open(input_file_1)

        # extract location variables
        lat = ncvar_get(data1, "Latitude") ; long = ncvar_get(data1, "Longitude")
        # read the NCSCD soil C prior
        Csom = ncvar_get(data1, "Csom")

        # expand the one directional values here into 2 directional
        lat_dim = length(lat) ; long_dim = length(long)
        long = array(long,dim=c(long_dim,lat_dim))
        lat = array(lat,dim=c(lat_dim,long_dim)) ; lat=t(lat)

        # Convert to a raster, assuming standad WGS84 grid
        Csom = data.frame(x = as.vector(long), y = as.vector(lat), z = as.vector(Csom))
        Csom = rast(Csom, crs = ("+init=epsg:4326"), type="xyz")

        # Create raster with the target crs (technically this bit is not required)
        target = rast(crs = ("+init=epsg:4326"), ext = ext(Csom), resolution = res(Csom))
        # Check whether the target and actual analyses have the same CRS
        if (compareGeom(Csom,target) == FALSE) {
            # Resample to correct grid
            Csom = resample(Csom, target, method="ngb") ; gc() 
        }
        # Extend the extent of the overall grid to the analysis domain
        Csom = extend(Csom,cardamom_ext)
        # Trim the extent of the overall grid to the analysis domain
        Csom = crop(Csom,cardamom_ext)
        # Adjust spatial resolution of the datasets, this occurs in all cases
        if (res(Csom)[1] != res(cardamom_ext)[1] | res(Csom)[2] != res(cardamom_ext)[2]) {

            # Create raster with the target resolution
            target = rast(crs = crs(cardamom_ext), ext = ext(cardamom_ext), resolution = res(cardamom_ext))
            # Resample to correct grid
            Csom = resample(Csom, target, method="bilinear") ; gc() 

        } # Aggrgeate to resolution

        # extract dimension information for the grid, note the axis switching between raster and actual array
        xdim = dim(Csom)[2] ; ydim = dim(Csom)[1]
        # extract the lat / long information needed
        long = crds(Csom,df=TRUE, na.rm=FALSE)
        lat  = long$y ; long = long$x
        # restructure into correct orientation
        long = array(long, dim=c(xdim,ydim))
        lat = array(lat, dim=c(xdim,ydim))
        # break out from the rasters into arrays which we can manipulate
        Csom = array(as.vector(unlist(Csom)), dim=c(xdim,ydim))

        # Assume in all cases than a zero prior value should be classed as missing data
        Csom[Csom < 1] = NA
        # assume uncertainty half that of HWSD as more targetted analysis, 0.5 * ~47 %
        Csom_unc = array(Csom * 0.47 * 0.5, dim=c(xdim,ydim))
        # With a minimum bound assumption
        Csom_unc[Csom_unc < 100] = 100
        # Ensure consistency for missing values
        Csom_unc[Csom == -9999 | is.na(Csom) == TRUE] = -9999

        # Return the loaded dataset
        return(list(Csom = Csom, Csom_unc = Csom_unc, lat = lat, long = long))

    } else if (Csom_source == "NCSCD3m") {

        # let the user know this might take some time
        print("Loading processed NCSCD Csom fields for subsequent sub-setting ...")

        # open processed modis files
        input_file_1 = paste(path_to_Csom,"/NCSCD_Csom_0-3m_with_lat_long.nc",sep="")
        data1 = nc_open(input_file_1)

        # extract location variables
        lat = ncvar_get(data1, "Latitude") ; long = ncvar_get(data1, "Longitude")
        # read the NCSCD soil C prior
        Csom = ncvar_get(data1, "Csom")

        # expand the one directional values here into 2 directional
        lat_dim = length(lat) ; long_dim = length(long)
        long = array(long,dim=c(long_dim,lat_dim))
        lat = array(lat,dim=c(lat_dim,long_dim)) ; lat=t(lat)

        # Convert to a raster, assuming standad WGS84 grid
        Csom = data.frame(x = as.vector(long), y = as.vector(lat), z = as.vector(Csom))
        Csom = rast(Csom, crs = ("+init=epsg:4326"), type="xyz")

        # Create raster with the target crs (technically this bit is not required)
        target = rast(crs = ("+init=epsg:4326"), ext = ext(Csom), resolution = res(Csom))
        # Check whether the target and actual analyses have the same CRS
        if (compareGeom(Csom,target) == FALSE) {
            # Resample to correct grid
            Csom = resample(Csom, target, method="ngb") ; gc() 
        }
        # Extend the extent of the overall grid to the analysis domain
        Csom = extend(Csom,cardamom_ext)
        # Trim the extent of the overall grid to the analysis domain
        Csom = crop(Csom,cardamom_ext)
        # Adjust spatial resolution of the datasets, this occurs in all cases
        if (res(Csom)[1] != res(cardamom_ext)[1] | res(Csom)[2] != res(cardamom_ext)[2]) {

            # Create raster with the target resolution
            target = rast(crs = crs(cardamom_ext), ext = ext(cardamom_ext), resolution = res(cardamom_ext))
            # Resample to correct grid
            Csom = resample(Csom, target, method="bilinear") ; gc() 

        } # Aggrgeate to resolution

        # extract dimension information for the grid, note the axis switching between raster and actual array
        xdim = dim(Csom)[2] ; ydim = dim(Csom)[1]
        # extract the lat / long information needed
        long = crds(Csom,df=TRUE, na.rm=FALSE)
        lat  = long$y ; long = long$x
        # restructure into correct orientation
        long = array(long, dim=c(xdim,ydim))
        lat = array(lat, dim=c(xdim,ydim))
        # break out from the rasters into arrays which we can manipulate
        Csom = array(as.vector(unlist(Csom)), dim=c(xdim,ydim))

        # Assume in all cases than a zero prior value should be classed as missing data
        Csom[Csom < 1] = NA
        # assume uncertainty, ~47 %
        Csom_unc = array(Csom * 0.47, dim=c(xdim,ydim))
        # With a minimum bound assumption
        Csom_unc[Csom_unc < 100] = 100
        # Ensure consistency for missing values
        Csom_unc[Csom == -9999 | is.na(Csom) == TRUE] = -9999

        # Return the loaded dataset
        return(list(Csom = Csom, Csom_unc = Csom_unc, lat = lat, long = long))

    } else {
        # output variables
	      return(list(Csom=-9999, Csom_unc = -9999, lat=-9999,long=-9999))
    }

} # function end load_Csom_fields_for_extraction

## Use byte compile
load_Csom_fields_for_extraction<-cmpfun(load_Csom_fields_for_extraction)
