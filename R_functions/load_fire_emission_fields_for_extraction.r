
###
## Function to load Fire C emisions from various gridded datasets
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_fire_emission_fields_for_extraction<-function(latlon_in,fire_source,start_year,end_year,cardamom_ext,spatial_type) {

    # Determine the years for which the analysis will occur
    years_to_do = as.character(seq(start_year,end_year))

    if (fire_source == "Global_Combined") {

        # let the user know this might take some time
        print("Loading Global_Combined Fire estimates for subsequent sub-setting ...")

        lat_done = FALSE ; missing_years = 0 ; keepers = 0 ; yrs = 1
        # loop for year here
        for (yr in seq(1, length(years_to_do))) {
             print(paste("... ",round((yr/length(years_to_do))*100,0),"% completed ",Sys.time(),sep=""))

             if (yr == 1) {
                 # list the available files
                 available_files = list.files(path_to_fire,full.names=TRUE)
                 # first check how many files we have
                 for (yrr in seq(1, length(years_to_do))) {
                      if (length(which(grepl(paste("Combined_FIRE_OBS_",years_to_do[yrr],sep=""),available_files))) > 0) {
                          keepers = keepers+1
                      } else {
                          missing_years = append(missing_years,years_to_do[yrr])
                      } # missing years
                 } # loop years to check which we have
	            } # first year

              # Determine the name of the first file first file files
  	          input_file_1 = paste(path_to_fire,"/Combined_FIRE_OBS_",years_to_do[yr],".nc",sep="")

              # check to see if file exists if it does then we read it in,
              # if not then we assume its a year we don't have data for and move on
	            if (file.exists(input_file_1) == TRUE) {

                  # open the file
                  data1 = nc_open(input_file_1)

                  # extract location variables
                  lat_in = ncvar_get(data1, "lat_axis") ; long_in = ncvar_get(data1, "long_axis")
                  # read the Fire emission estimate (units are gC/m2/day)
                  var1_in = ncvar_get(data1, "Fire")
                  # get some uncertainty information - in this case the max / min which will be the basis of our uncertainty
                  var2_in = ncvar_get(data1, "Fire_min") ; var3_in = ncvar_get(data1, "Fire_max")
                  var2_in = (var3_in - var2_in) * 0.5 ; rm(var3_in)
                  # extract timing information
                  if (lat_done == FALSE) {
                      # Months in a year
                      doy_obs = ncvar_get(data1, "time_axis")
                  } else {
                      # Months in a year
                      tmp1 = ncvar_get(data1, "time_axis")
                      doy_obs = append(doy_obs,tmp1)
                  }
                  # close files after use
                  nc_close(data1)

                  # Turn lat_in / long_in from vectors to arrays
                  lat_in = t(array(lat_in, dim=c(dim(var1_in)[2],dim(var1_in)[1])))
                  long_in = array(long_in, dim=c(dim(var1_in)[1],dim(var1_in)[2]))

                  # Loop through each timestep in the year
                  for (t in seq(1, dim(var1_in)[3])) {
                       # Convert to a raster, assuming standad WGS84 grid
                       var1 = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(var1_in[,,t]))
                       var1 = rasterFromXYZ(var1, crs = ("+init=epsg:4326"))
                       var2 = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(var2_in[,,t]))
                       var2 = rasterFromXYZ(var2, crs = ("+init=epsg:4326"))

                       # Create raster with the target crs (technically this bit is not required)
                       target = raster(crs = ("+init=epsg:4326"), ext = extent(var1), resolution = res(var1))
                       # Check whether the target and actual analyses have the same CRS
                       if (compareCRS(var1,target) == FALSE) {
                           # Resample to correct grid
                           var1 = resample(var1, target, method="ngb") ; gc() ; removeTmpFiles()
                           var2 = resample(var2, target, method="ngb") ; gc() ; removeTmpFiles()
                       }
                       # Extend the extent of the overall grid to the analysis domain
                       var1 = extend(var1,cardamom_ext) ; var2 = extend(var2,cardamom_ext)
                       # Trim the extent of the overall grid to the analysis domain
                       var1 = crop(var1,cardamom_ext) ; var2 = crop(var2,cardamom_ext)

                       # Adjust spatial resolution of the datasets, this occurs in all cases
                       if (res(var1)[1] != res(cardamom_ext)[1] | res(var1)[2] != res(cardamom_ext)[2]) {

                           # Create raster with the target resolution
                           target = raster(crs = crs(cardamom_ext), ext = extent(cardamom_ext), resolution = res(cardamom_ext))
                           # Resample to correct grid
                           var1 = resample(var1, target, method="bilinear") ; gc() ; removeTmpFiles()
                           var2 = resample(var2, target, method="bilinear") ; gc() ; removeTmpFiles()
                       } # Aggrgeate to resolution

                       if (lat_done == FALSE) {
                           # extract dimension information for the grid, note the axis switching between raster and actual array
                           xdim = dim(var1)[2] ; ydim = dim(var1)[1]
                           # extract the lat / long information needed
                           long = coordinates(var1)[,1] ; lat = coordinates(var1)[,2]
                           # restructure into correct orientation
                           long = array(long, dim=c(xdim,ydim))
                           lat = array(lat, dim=c(xdim,ydim))
                       }
                       # break out from the rasters into arrays which we can manipulate
                       var1 = array(as.vector(unlist(var1)), dim=c(xdim,ydim))
                       var2 = array(as.vector(unlist(var2)), dim=c(xdim,ydim))

                       # vectorise at this time
                       if (lat_done == FALSE) {
                           fire_gCm2day = as.vector(var1)
                           fire_unc_gCm2day = as.vector(var2)
                       } else {
                           fire_gCm2day = append(fire_gCm2day,as.vector(var1))
                           fire_unc_gCm2day = append(fire_unc_gCm2day,as.vector(var2))
                       }

                       # update flag for lat / long load
                       if (lat_done == FALSE) {lat_done = TRUE}

                  } # loop steps within year

                  # keep track of years actually ran
                  yrs = yrs+1

	            } # end of does file exist

  	     } # year loop

         # remove initial value
         missing_years=missing_years[-1]

         # clean up variables
         gc(reset=TRUE,verbose=FALSE)

         # restructure
         fire_gCm2day = array(fire_gCm2day, dim=c(xdim,ydim,length(doy_obs)))
         fire_unc_gCm2day = array(fire_unc_gCm2day, dim=c(xdim,ydim,length(doy_obs)))

         # output variables
         return(list(fire_gCm2day = fire_gCm2day, fire_unc_gCm2day = fire_unc_gCm2day,
                     doy_obs = doy_obs, lat = lat, long = long, missing_years = missing_years))

    } else if (fire_source == " " | fire_source == "site_specific"){

	        # Do nothing as this should be read directly from files or not needed
          return(list(fire_gCm2day = -9999, fire_unc_gCm2day = -9999,
                      doy_obs = -9999, lat = -9999, long = -9999, missing_years = -9999))

    } else {

	        stop(paste("FIRE option (",fire_source,") not valid"))

    } # if Fire type

} # function end

## Use byte compile
load_fire_emission_fields_for_extraction<-cmpfun(load_fire_emission_fields_for_extraction)
