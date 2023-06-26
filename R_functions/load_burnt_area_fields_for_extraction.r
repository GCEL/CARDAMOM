
###
## Function to load burned area data from various gridded datasets
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_burnt_area_fields_for_extraction<-function(latlon_in,burnt_area_source,path_to_burnt_area,start_year,end_year,cardamom_ext,spatial_type) {

    # Determine the years for which the analysis will occur
    years_to_do = as.character(seq(start_year,end_year))

    if (burnt_area_source == "GFED4") {

        # let the user know this might take some time
        print("Loading GFED4 processed burnt area estimates for subsequent sub-setting ...")

        lat_done = FALSE ; missing_years = 0 ; keepers = 0 ; yrs = 1
        # loop for year here
        for (yr in seq(1, length(years_to_do))) {
             print(paste("... ",round((yr/length(years_to_do))*100,0),"% completed ",Sys.time(),sep=""))

             if (yr == 1) {
                 # list the available files
                 available_files = list.files(path_to_burnt_area,full.names=TRUE)
                 # first check how many files we have
                 for (yrr in seq(1, length(years_to_do))) {
                      if (length(which(grepl(years_to_do[yrr],available_files))) > 0) {
                          keepers = keepers+1
                      } else {
                          missing_years = append(missing_years,years_to_do[yrr])
                      } # missing years
                 } # loop years to check which we have
	            } # first year
              # open processed modis files
  	          input_file_1 = paste(path_to_burnt_area,"/GFED4_",years_to_do[yr],".nc",sep="")

              # check to see if file exists if it does then we read it in,
              # if not then we assume its a year we don't have data for and move on
	            if (file.exists(input_file_1) == TRUE) {

                  # open the file
                  data1 = nc_open(input_file_1)

                  # extract location variables
                  lat_in = ncvar_get(data1, "latitude") ; long_in = ncvar_get(data1, "longitude")
                  # read the burnt fraction estimate (units are 0-1)
                  var1 = ncvar_get(data1, "BurnedFraction")
                  # get time information (month in this case)
                  var2 = ncvar_get(data1, "time") ; time_steps_per_year = 12
                  # approximate doy of the mid-month and allocate fire to that point
                  if (lat_done == FALSE) {
                      doy_obs = floor((var2*(365.25/12))-(365.25/24))
                  } else {
                      doy_obs = append(doy_obs,floor((var2*(365.25/12))-(365.25/24)))
                  }

                  # close files after use
                  nc_close(data1) ; rm(var2)

                  # Turn lat_in / long_in from vectors to arrays
                  lat_in = t(array(lat_in, dim=c(dim(var1)[2],dim(var1)[1])))
                  long_in = array(long_in, dim=c(dim(var1)[1],dim(var1)[2]))

                  # Convert to a raster, assuming standad WGS84 grid
                  var1 = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(var1))
                  var1 = rasterFromXYZ(var1, crs = ("+init=epsg:4326"))
                  # Remove the input lat / long information
                  rm(lat_in,long_in)

                  # extend the extent of the overall grid to the analysis domain
                  var1 = extend(var1,cardamom_ext)
                  # Trim the extent of the overall grid to the analysis domain
                  var1 = crop(var1,cardamom_ext)
                  # set actual missing data to 0 as missing data is actually no fire
                  var1[which(is.na(as.vector(var1)))] = 0
                   # Adjust spatial resolution of the datasets, this occurs in all cases
                   if (res(var1)[1] != res(cardamom_ext)[1] | res(var1)[2] != res(cardamom_ext)[2]) {
                       # Create raster with the target resolution
                       target = raster(crs = crs(cardamom_ext), ext = extent(cardamom_ext), resolution = res(cardamom_ext))
                       # Resample to correct grid.
                       # Probably should be done via aggregate function to allow for correct error propogation
                       var1 = resample(var1, target, method="bilinear") ; gc() ; removeTmpFiles()
                   } # Aggrgeate to resolution

                  # Extract spatial information just the once
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

                  # vectorise at this time
                  if (lat_done == FALSE) {
                      burnt_area = as.vector(var1)
                  } else {
                      burnt_area = append(burnt_area,as.vector(var1))
                  }

                  # update flag for lat / long load
                  if (lat_done == FALSE) {lat_done = TRUE}
                  # keep track of years actually ran
                  yrs = yrs+1

	            } # end of does file exist

  	     } # year loop

         # remove initial value
         missing_years=missing_years[-1]

         # clean up variables
         rm(var1,var2) ; gc(reset=TRUE,verbose=FALSE)

         # restructure
         burnt_area=array(burnt_area, dim=c(xdim,ydim,length(doy_obs)))

         # output variables
         return(list(burnt_area=burnt_area,doy_obs=doy_obs,lat=lat,long=long,missing_years=missing_years))

    } else if (burnt_area_source == "MCD64A1") {

        # let the user know this might take some time
        print("Loading processed MCD61A1 burnt fraction fields for subsequent sub-setting ...")

        # check which file prefix we are using today
        # list all available files which we will then search
        avail_files = list.files(path_to_burnt_area,full.names=TRUE)
        #prefix = "MCD15A2H_LAI_(.)*" # (.)* wildcard characters for unix standard MCD15A2H_LAI_*
        prefix = "MCD64A1_BurnedFraction_"

        # timing information on the number of day in a month
        month_days = rep(31,length.out=12)
        month_days[2] = 28 ; month_days[c(4,6,9,11)] = 30

        lat_done = FALSE ; missing_years = 0 ; keepers = 0 ; yrs = 1 ; doy_out = 0
        # loop for year here
        for (yr in seq(1, length(years_to_do))) {

             # Update the user as to our progress
             print(paste("... ",round((yr/length(years_to_do))*100,0),"% completed ",Sys.time(),sep=""))

             # If we are just starting check how many files we have
             if (yr == 1) {
                 nsteps = 0
                 # Loop through all the analyses years and check whether files exist for it
                 for (yrr in seq(1, length(years_to_do))) {
                      # create the prefix for the observation file files we will want for a given year
                      input_file_1 = paste(prefix,years_to_do[yrr],sep="")
                      # then check whether this pattern is found in the available files
                      this_year = grepl(input_file_1, avail_files) ; this_year = which(this_year == TRUE)
                      # If we have at least one timestep for this year then we have some information otherwise it is missing!
                      if (length(this_year) > 0) {
                          # Track the number of obserations
                          keepers = keepers+1 ; nsteps = max(nsteps,length(this_year))
                      } else {
                          missing_years = append(missing_years,years_to_do[yrr])
                      }
                 } # loop through possible years
                 rm(yrr)
             } # first year?

             # Begin reading the files in now for real

             # Determine the unique file name pattern
             input_file_1 = paste(prefix,years_to_do[yr],sep="")

             # Then check whether this pattern is found in the available files
             this_year = avail_files[grepl(input_file_1, avail_files)]
             if (length(this_year) > 0) {

                 # The files should be in order due to the YYYYDOY format used
                 this_year = this_year[order(this_year)]

                 # Loop through the available files for the current year
                 for (t in seq(1, length(this_year))) {

                      # Inform user
                      #print(paste("...reading the following data file = ",this_year[t],sep=""))
                      # open the file
                      data1 = nc_open(this_year[t])

                      # Get timing variable
                      doy_in = ncvar_get(data1, "doy")
                      # Extract spatial information
                      lat_in = ncvar_get(data1, "lat") ; long_in = ncvar_get(data1, "lon")
                      # read the LAI observations
                      var1 = ncvar_get(data1, "BurnedFraction") # Fraction of pixel burned
                      # Close the current file
                      nc_close(data1)

                      # Convert to a raster, assuming standad WGS84 grid
                      var1 = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(var1))
                      var1 = rasterFromXYZ(var1, crs = ("+init=epsg:4326")) ; f1 = filename(var1)
                      # Remove the input lat / long information
                      rm(lat_in,long_in,lat_in_sd)

                      # Extend the extent of the overall grid to the analysis domain
                      var1 = extend(var1,cardamom_ext)
                      ff1 = filename(var1)

                      # Check and remove unwanted tmp files
                      if (f1 != ff1 & f1 != "") {
                          if (file.exists(f1)) { file.remove(f1) }
                          if (file.exists(gsub(".grd",".gri",f1))) { file.remove(gsub(".grd",".gri",f1)) }
                      } ; f1 = filename(var1)
                      # Trim the extent of the overall grid to the analysis domain
                      var1 = crop(var1,cardamom_ext)
                      ff1 = filename(var1)

                      # Check and remove unwanted tmp files
                      if (f1 != ff1 & f1 != "") {
                          if (file.exists(f1)) { file.remove(f1) }
                          if (file.exists(gsub(".grd",".gri",f1))) { file.remove(gsub(".grd",".gri",f1)) }
                      } ; f1 = filename(var1)
                      # Adjust spatial resolution of the datasets, this occurs in all cases
                      if (res(var1)[1] != res(cardamom_ext)[1] | res(var1)[2] != res(cardamom_ext)[2]) {
                          # Create raster with the target resolution
                          target = raster(crs = crs(cardamom_ext), ext = extent(cardamom_ext), resolution = res(cardamom_ext))
                          # Resample to correct grid.
                          # Probably should be done via aggregate function to allow for correct error propogation
                          var1 = resample(var1, target, method="bilinear") ; gc() ; removeTmpFiles()
                          ff1 = filename(var1)
                          # Check and remove unwanted tmp files
                          if (f1 != ff1 & f1 != "") {
                              if (file.exists(f1)) { file.remove(f1) }
                              if (file.exists(gsub(".grd",".gri",f1))) { file.remove(gsub(".grd",".gri",f1)) }
                          } ; f1 = filename(var1)

                      } # Aggrgeate to resolution

                      # Extract spatial information just the once
                      if (lat_done == FALSE) {
                          # Set flag to true
                          lat_done = TRUE
                          # extract dimension information for the grid, note the axis switching between raster and actual array
                          xdim = dim(var1)[2] ; ydim = dim(var1)[1]
                          # extract the lat / long information needed
                          long = coordinates(var1)[,1] ; lat = coordinates(var1)[,2]
                          # restructure into correct orientation
                          long = array(long, dim=c(xdim,ydim))
                          lat = array(lat, dim=c(xdim,ydim))
                          # create holding arrays for the burned information...
                          burnt_area_hold = array(NA, dim=c(xdim*ydim,keepers*nsteps))
                      }
                      # break out from the rasters into arrays which we can manipulate
                      var1 = array(as.vector(unlist(var1)), dim=c(xdim,ydim))

                      # Check and remove unwanted tmp files
                      if (f1 != "") {
                          if (file.exists(f1)) { file.remove(f1) }
                          if (file.exists(gsub(".grd",".gri",f1))) { file.remove(gsub(".grd",".gri",f1)) }
                      } # no longer using rasters at this point, so no file name to track
                      # set actual missing data to -9999
                      var1[which(is.na(as.vector(var1)))] = -9999

                      # begin populating the various outputs
                      burnt_area_hold[1:length(as.vector(var1)),(t+((yrs-1)*nsteps))] = as.vector(var1)
                      doy_out = append(doy_out,doy_in)

                 } # loop through available time steps in the current year

                 # keep track of years actually ran
                 yrs = yrs + 1
                 # clean up allocated memeory
                 rm(var1) ; gc()

             } # is there information for the current year?

        } # year loop

        # Correct for initialisation
        doy_out = doy_out[-1]

        # Sanity check for LAI
        if (lat_done == FALSE) {stop('No burned fraction information could be found...')}

        # remove initial value
        missing_years = missing_years[-1]

        # check which ones are NA because I made them up
        not_na = is.na(as.vector(burnt_area_hold))
        not_na = which(not_na == FALSE)

        filter = as.vector(burnt_area_hold) == -9999
        # now remove the ones that are actual missing data
        burnt_area_hold[filter] = NA
        # return spatial structure to data
        burnt_area_out = array(as.vector(burnt_area_hold)[not_na], dim=c(xdim,ydim,length(doy_out)))

        # output variables
        burn_area_all = list(burnt_area = burnt_area_out,
                             doy_obs = doy_out,
                             lat = lat, long = long, missing_years=missing_years)
        # clean up variables
        rm(doy_in,burnt_area_hold,not_na,lai_out,doy_out,lat,long,missing_years)
        gc(reset=TRUE,verbose=FALSE)
        return(burn_area_all)

    } else if (burnt_area_source == " " | burnt_area_source == "site_specific"){

        # Do nothing as this should be read directly from files or not needed

    } else {

        stop(paste("Burnt area option (",burnt_area_source,") not valid"))

    } # if data_source

} # function end

## Use byte compile
load_burnt_area_fields_for_extraction<-cmpfun(load_burnt_area_fields_for_extraction)
