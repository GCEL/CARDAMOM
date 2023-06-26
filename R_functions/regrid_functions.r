
###
## Functions to aggregate spatial information to coarser resolutions
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

# Define local function
regrid_func<-function(var1_in, lat_in, long_in, cardamom_ext, landmask) {

   # Set flags
   lat_done = FALSE

   # Check dimensions of the input variable
   if (length(dim(var1_in)) == 2) {
       var1_in = array(var1_in, dim=c(dim(var1_in),1))
   }

   # Loop through each timestep in the year
   for (t in seq(1, dim(var1_in)[3])) {
        # Convert to a raster, assuming standard WGS84 grid
        var1 = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(var1_in[,,t]))
        if (length(unique(diff(input_long[,1]))) > 1 | length(unique(diff(input_lat[1,]))) > 1) {
            var1 = griddify(var1, dim(var1_in)[1], dim(var1_in)[2])
        } else {
            var1 = rasterFromXYZ(var1, crs = ("+init=epsg:4326"))
        }

        # Create raster with the target crs (technically this bit is not required)
        target = raster(crs = ("+init=epsg:4326"), ext = extent(var1), resolution = res(var1))
        # Check whether the target and actual analyses have the same CRS
        if (compareCRS(var1,target) == FALSE) {
            # Resample to correct grid
            var1 = resample(var1, target, method="ngb") ; gc() ; removeTmpFiles()
        }
        # Extend the extent of the overall grid to the analysis domain
        var1 = extend(var1,cardamom_ext, snap="near", value = NA)
        # Trim the extent of the overall grid to the analysis domain
        var1 = crop(var1,cardamom_ext)

        # If this is a gridded analysis and the desired CARDAMOM resolution is coarser than the currently provided then aggregate here.
        # Despite creation of a cardamom_ext for a site run do not allow aggregation here as this will damage the fine resolution datasets
        if (res(var1)[1] != res(cardamom_ext)[1] | res(var1)[2] != res(cardamom_ext)[2]) {

            # Create raster with the target resolution
            target = raster(crs = crs(cardamom_ext), ext = extent(cardamom_ext), resolution = res(cardamom_ext))
            # Resample to correct grid
            var1 = resample(var1, target, method="bilinear") ; gc() ; removeTmpFiles()

        } # Aggregate to resolution

        # If a land mask is present then also restrict to this target domain
        if (missing("landmask") == FALSE) {var1 = crop(var1,landmask)}

        if (lat_done == FALSE) {
            # extract dimension information for the grid, note the axis switching between raster and actual array
            xdim = dim(var1)[2] ; ydim = dim(var1)[1]
            # extract the lat / long information needed
            long = coordinates(var1)[,1] ; lat = coordinates(var1)[,2]
            # restructure into correct orientation
            long = array(long, dim=c(xdim,ydim))
            lat = array(lat, dim=c(xdim,ydim))
        }
        # break out from the raster into arrays which we can manipulate
        var1 = array(as.vector(unlist(var1)), dim=c(xdim,ydim))

        # vectorise at this time
        if (lat_done == FALSE) {
            var_out = as.vector(var1)
        } else {
            var_out = append(var_out,as.vector(var1))
        }

        # update flag for lat / long load
        if (lat_done == FALSE) {lat_done = TRUE}

   } # Within variable time loop

   # restructure
   var_out = array(var_out, dim=c(xdim,ydim,dim(var1_in)[3]))

   # Return to user
   return(list(var = var_out, lat = lat, long = long))

} # end function regrid_func

# Function to use gdal libraries to aggregate to target spatial resolution
regrid_gdal_func<-function(tmp_dir, var1_in, lat_in, long_in, cardamom_ext, landmask) {

   # Load libraries needed
   require(gdalUtils)
   # Help avoid build up of raster tmp files
   rasterOptions(overwrite = TRUE)

   # Set flags
   lat_done = FALSE

   # Check dimensions of the input variable
   if (length(dim(var1_in)) == 2) {
       var1_in = array(var1_in, dim=c(dim(var1_in),1))
   }

   # Store current working directory
   current_wd = getwd()
   # Change to tmp working directory
   setwd(tmp_dir)
   # Create file names for intermediate files
   outfile_tmp = paste(tmp_dir,"/tmp.tif",sep="")
   outfile_agg = paste(tmp_dir,"/tmp_agg.tif",sep="")

   # Loop through each timestep in the year
   for (t in seq(1, dim(var1_in)[3])) {

        # Convert to a raster, assuming standard WGS84 grid
        var1 = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(var1_in[,,t]))
        if (length(unique(diff(input_long[,1]))) > 1 | length(unique(diff(input_lat[1,]))) > 1) {
            var1 = griddify(var1, dim(var1_in)[1], dim(var1_in)[2])
        } else {
            var1 = rasterFromXYZ(var1, crs = ("+init=epsg:4326"))
        }

        # Create raster with the target crs (technically this bit is not required)
        target = raster(crs = ("+init=epsg:4326"), ext = extent(var1), resolution = res(var1))
        # Check whether the target and actual analyses have the same CRS
        if (compareCRS(var1,target) == FALSE) {
            # Resample to correct grid
            var1 = resample(var1, target, method="ngb") ; gc() ; removeTmpFiles()
        }
        # Extend the extent of the overall grid to the analysis domain
        var1 = extend(var1,cardamom_ext, snap="near", value = NA)
        # Trim the extent of the overall grid to the analysis domain
        var1 = crop(var1,cardamom_ext)

        # Aggregate to coarser spatial resolution using gdal libraries as these are faster
        if (res(var1)[1] != res(cardamom_ext)[1] | res(var1)[2] != res(cardamom_ext)[2]) {
            # Write out file to temporary raster file
            writeRaster(var1, outfile_tmp, format = "GTiff", overwrite=TRUE)
            # Carry out aggregation using gdal libraries
            gdal_translate(src_dataset = outfile_tmp, dst_dataset = outfile_agg,
                           a_srs = "EPSG:4326", of = "GTiff", tr = res(cardamom_ext),
                           r = "average")
#            gdalwarp(srcfile = outfile_tmp, dstfile = outfile_agg,
#                     a_srs = "EPSG:4326", of = "GTiff", tr = res(cardamom_ext),
#                     r = "average")
            # Load back aggregated variable into memory
            var1 = raster(outfile_agg)
            # Tidy away the intermediate files
            file.remove(outfile_tmp)
        }

        # If a land mask is present then also restrict to this target domain
        if (missing("landmask") == FALSE) {var1 = crop(var1,landmask)}

        if (lat_done == FALSE) {
            # extract dimension information for the grid, note the axis switching between raster and actual array
            xdim = dim(var1)[2] ; ydim = dim(var1)[1]
            # extract the lat / long information needed
            long = coordinates(var1)[,1] ; lat = coordinates(var1)[,2]
            # restructure into correct orientation
            long = array(long, dim=c(xdim,ydim))
            lat = array(lat, dim=c(xdim,ydim))
        }
        # break out from the raster into arrays which we can manipulate
        var1 = array(as.vector(unlist(var1)), dim=c(xdim,ydim))

        # vectorise at this time
        if (lat_done == FALSE) {
            var_out = as.vector(var1)
        } else {
            var_out = append(var_out,as.vector(var1))
        }

        # Ensure var1 and assocated tmp files have been removed
        rm(var1) ; if (file.exists(outfile_agg)) {file.remove(outfile_agg)}

        # update flag for lat / long load
        if (lat_done == FALSE) {lat_done = TRUE}

   } # Within variable time loop

   # restructure
   var_out = array(var_out, dim=c(xdim,ydim,dim(var1_in)[3]))

   # Return to original working directory
   setwd(current_wd)

   # Return to user
   return(list(var = var_out, lat = lat, long = long))

} # end function regrid_gdal_func
