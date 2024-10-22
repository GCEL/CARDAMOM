
###
## Function to load Net Biome Exchange for subsquent subsetting
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_nbe_fields_for_extraction<-function(latlon_in,nbe_source,years_to_load,cardamom_ext,spatial_type) {

    if (nbe_source == "Global_Combined") {

        # let the user know this might take some time
        print("Loading Global_Combined NBE estimates for subsequent sub-setting ...")

        lat_done = FALSE ; missing_years = 0 ; keepers = 0 ; yrs = 1
        # loop for year here
        for (yr in seq(1, length(years_to_do))) {
             print(paste("... ",round((yr/length(years_to_do))*100,0),"% completed ",Sys.time(),sep=""))

             if (yr == 1) {
                 # list the available files
                 available_files = list.files(path_to_nbe,full.names=TRUE)
                 # first check how many files we have
                 for (yrr in seq(1, length(years_to_do))) {
                      if (length(which(grepl(paste("Combined_NBE_OBS_",years_to_do[yrr],sep=""),available_files))) > 0) {
                          keepers = keepers+1
                      } else {
                          missing_years = append(missing_years,years_to_do[yrr])
                      } # missing years
                 } # loop years to check which we have
              } # first year

              # Determine the name of the first file first file files
              input_file_1 = paste(path_to_nbe,"/Combined_NBE_OBS_",years_to_do[yr],".nc",sep="")

              # check to see if file exists if it does then we read it in,
              # if not then we assume its a year we don't have data for and move on
              if (file.exists(input_file_1) == TRUE) {

                  # open the file
                  data1 = nc_open(input_file_1)

                  # extract location variables
                  lat_in = ncvar_get(data1, "lat_axis") ; long_in = ncvar_get(data1, "long_axis")
                  # read the NBE estimate (units are gC/m2/day)
                  var1_in = ncvar_get(data1, "NBE")
                  # get some uncertainty information - in this case the max / min which will be the basis of our uncertainty
                  var2_in = ncvar_get(data1, "NBE_min") ; var3_in = ncvar_get(data1, "NBE_max")
                  var2_in = (var3_in - var2_in) * 0.5 ; rm(var3_in)
                  # Months in a year
                  time_steps_per_year = 12
                  # approximate doy of the mid-month and allocate nbe to that point
                  if (lat_done == FALSE) {
                      doy_obs = floor((c(1:time_steps_per_year)*(365.25/12))-(365.25/24))
                  } else {
                      doy_obs = append(doy_obs,floor((c(1:time_steps_per_year)*(365.25/12))-(365.25/24)))
                  }

                  # close files after use
                  nc_close(data1)

                  # Turn lat_in / long_in from vectors to arrays
                  lat_in = t(array(lat_in[length(lat_in):1], dim=c(dim(var1_in)[2],dim(var1_in)[1])))
                  long_in = array(long_in, dim=c(dim(var1_in)[1],dim(var1_in)[2]))

                  # Loop through each timestep in the year
                  for (t in seq(1, dim(var1_in)[3])) {
                       # Convert to a raster, assuming standad WGS84 grid
                       var1 = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(var1_in[,,t]))
                       var1 = rast(var1, crs = ("+init=epsg:4326"), type="xyz")
                       var2 = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(var2_in[,,t]))
                       var2 = rast(var2, crs = ("+init=epsg:4326"), type="xyz")

                       # Create raster with the target crs (technically this bit is not required)
                       target = rast(crs = ("+init=epsg:4326"), ext = ext(var1), resolution = res(var1))
                       # Check whether the target and actual analyses have the same CRS
                       if (compareGeom(var1,target) == FALSE) {
                           # Resample to correct grid
                           var1 = resample(var1, target, method="ngb") ; gc() 
                           var2 = resample(var2, target, method="ngb") ; gc() 
                       }
                       # Extend the extent of the overall grid to the analysis domain
                       var1 = extend(var1,cardamom_ext) ; var2 = extend(var2,cardamom_ext)
                       # Trim the extent of the overall grid to the analysis domain
                       var1 = crop(var1,cardamom_ext) ; var2 = crop(var2,cardamom_ext)

                       # Adjust spatial resolution of the datasets, this occurs in all cases
                       if (res(var1)[1] != res(cardamom_ext)[1] | res(var1)[2] != res(cardamom_ext)[2]) {

                           # Create raster with the target resolution
                           target = rast(crs = crs(cardamom_ext), ext = ext(cardamom_ext), resolution = res(cardamom_ext))
                           # Resample to correct grid
                           var1 = resample(var1, target, method="bilinear") ; gc()
                           var2 = resample(var2, target, method="bilinear") ; gc() 

                       } # Aggrgeate to resolution

                       if (lat_done == FALSE) {
                           # extract dimension information for the grid, note the axis switching between raster and actual array
                           xdim = dim(var1)[2] ; ydim = dim(var1)[1]
                           # extract the lat / long information needed
                           long = crds(var1,df=TRUE, na.rm=FALSE)
                           lat  = long$y ; long = long$x
                           # restructure into correct orientation
                           long = array(long, dim=c(xdim,ydim))
                           lat = array(lat, dim=c(xdim,ydim))
                       }
                       # break out from the rasters into arrays which we can manipulate
                       var1 = array(as.vector(unlist(var1)), dim=c(xdim,ydim))
                       var2 = array(as.vector(unlist(var2)), dim=c(xdim,ydim))

                       # vectorise at this time
                       if (lat_done == FALSE) {
                           nbe_gCm2day = as.vector(var1)
                           nbe_unc_gCm2day = as.vector(var2)
                       } else {
                           nbe_gCm2day = append(nbe_gCm2day,as.vector(var1))
                           nbe_unc_gCm2day = append(nbe_unc_gCm2day,as.vector(var2))
                       }

                       # update flag for lat / long load
                       if (lat_done == FALSE) {lat_done = TRUE}

                  } # within year step

                  # keep track of years actually ran
                  yrs = yrs+1

              } # end of does file exist

        } # year loop

        # remove initial value
        missing_years=missing_years[-1]

        # clean up variables
        gc(reset=TRUE,verbose=FALSE)

        # restructure
        nbe_gCm2day = array(nbe_gCm2day, dim=c(xdim,ydim,length(doy_obs)))
        nbe_unc_gCm2day = array(nbe_unc_gCm2day, dim=c(xdim,ydim,length(doy_obs)))

        # output variables
        return(list(retrieval_valid = TRUE, nbe_gCm2day = nbe_gCm2day, nbe_unc_gCm2day = nbe_unc_gCm2day,
                    doy_obs = doy_obs, lat = lat, long = long, missing_years = missing_years))

  } else if (nbe_source == "GEOSCHEM_GCP") {

      # let the user know this might take some time
      print("Loading processed GEOCHEM_GCP NBE fields for subsequent sub-setting ...")

      # check which file prefix we are using today
      # list all available files which we will then search
      avail_files = list.files(path_to_nbe,full.names=TRUE, pattern = "gcp2023_v3_uoe_1x1_2001b")
      if (length(avail_files) != 1) {avail_files = list.files(path_to_nbe,full.names=TRUE, pattern = "merra2_oco2_insitu_v10r_odiac_uoe")}
      if (length(avail_files) != 1) {avail_files = list.files(path_to_nbe,full.names=TRUE, pattern = "merra2_gosat_insitu_uoe_2009_2019_full")}
      if (length(avail_files) != 1) {avail_files = list.files(path_to_nbe,full.names=TRUE, pattern = "oco2_v11r_lnlgogis_uoe_2014_full")}
      #prefix = "gcp2023_v3_uoe_1x1_2001b" # (.)* wildcard characters for unix standard c_gls*_
      # How many files do we have to work with? Hopefully just the one
      if (length(avail_files) != 1) {stop("We do not have a single GEOSCHEM_GCP NBE estimate. Correct and re-run.")}
      
      # Open link to the file
      data1 = nc_open(avail_files)
      # Load the time variable from the file
      # This should come in format [YYYY,MM,DD]
      # NOTE the script assumes that the file is at a monthly timestep and that only completed years are provided
      time = ncvar_get(data1, "start_date")
      # Check that the file has whole years matching our monthly assumption
      if (dim(time)[2]/12 != floor(dim(time)[2]/12)) {
          stop("The number of time steps in GEOSCHEM_GCP file is not consistent with our assumptions of monthly timesteps and whole years only")
      }

      # Determimine which years of the analysis are included in the file
      obs_years = unique(time[1,])
      # Determine the common years
      common_years = intersect(obs_years,years_to_load)
      
      # If there is no overlap then we will just move on from this dataset
      if (length(common_years) == 0) {
          return(retrieval_valid = FALSE)
      }
      
      # Do all required years have observations?
      if (length(common_years) != length(years_to_load)) {
          # No, so we will work out which are missing
          missing_years = 0
          for (yr in seq(1, length(years_to_load))) {
               if (length(which(common_years == years_to_load[yr])) > 0) {
                   # For the moment do nothing!               
               } else { 
                   # Track which year is missing
                   missing_years = append(missing_years,years_to_load[yr])
               }
          }
      } else {
          # Yes, so update appropriate variables
          missing_years = 0
      }
      # remove initial value
      missing_years = missing_years[-1]

      # Based on the common years create a filter for the timesteps we will be keeping
      keeper_steps = 0 
      for (yr in seq(1, length(common_years))) {
           keeper_steps = append(keeper_steps,which(time[1,] == common_years[yr]))
      }
      keeper_steps = keeper_steps[-1]

      # timing information on the number of day in a month
      month_days = rep(31,length.out=12)
      month_days[2] = 28 ; month_days[c(4,6,9,11)] = 30
      # get timing variable
      doy_obs = time[3,]
      for (d in seq(1,length(time[3,]))) {
           # January is correct already, so only adjust if month is >= February
           if (time[2,d] > 1) { # Month counter
               doy_obs[d] = doy_obs[d] + sum(month_days[1:(time[2,d]-1)])
           }
      }
      # Apply keeper filer
      doy_obs = doy_obs[keeper_steps]
      
      # Extract latitude and longitude information
      lat_in = ncvar_get(data1, "lat") ; long_in = ncvar_get(data1, "lon")
      # Store their dimensions
      lat_dim = length(lat_in) ; long_dim = length(long_in)
      # Turn lat_in / long_in from vectors to arrays
      lat_in = t(array(lat_in, dim=c(lat_dim,long_dim)))
      long_in = array(long_in, dim=c(long_dim,lat_dim))

      # Read in biospheric flux NBE
      var1 = ncvar_get(data1, "bio") # net biome exchange of CO2 (kgC/m2/s)
      # Read in biospheric flux NBE uncertainty
      if (length(which(grepl("post_flux_uncert",names(data1$var)) == TRUE)) > 0) {
          var2 = ncvar_get(data1, "post_flux_uncert") # NBE error estimate (kgC/m2/s)
      } else if (length(which(grepl("uncertainty",names(data1$var)) == TRUE)) > 0) {
          var2 = ncvar_get(data1, "uncertainty") # NBE error estimate (kgC/m2/s)
      } else {
          
      }
      # Close file connection
      nc_close(data1)
      
      # Filter based on the keeper_steps
      var1 = var1[,,keeper_steps] ; var2 = var2[,,keeper_steps]
      # Reset missing data to NA
      # NOTE these get converted to -9999 in extract_obe.r
      var1[which(var1 == -999)] = NA ; var2[which(var2 == -999)] = NA
      # Apply unit correction kgC/m2/s -> gC/m2/day
      var1 = var1 * 86400 * 1e3
      var2 = var2 * 86400 * 1e3

      # Aggregate to target resolution and extent
      output = regrid_func(var1, lat_in, long_in, cardamom_ext)
      nbe_out = output$var ; lat = output$lat ; long = output$long ; rm(output,var1)
      output = regrid_func(var2, lat_in, long_in, cardamom_ext)
      nbe_unc_out = output$var ; rm(output,var2)

      # enforce minimum uncertainty value
      nbe_unc_out[nbe_unc_out < 0.01] = 0.01 # (gC/m2/day)

      # output variables
      nbe_all = list(retrieval_valid = TRUE, nbe_gCm2day = nbe_out, nbe_unc_gCm2day = nbe_unc_out,
                     doy_obs = doy_obs, lat = lat, long = long, missing_years = missing_years)

      # Tidy up
      rm(lat_in,long_in,time,nbe_out,nbe_unc_out,lat,long,doy_obs,missing_years,lat_in,long_in)
      gc(reset=TRUE,verbose=FALSE)

      # Return to function
      return(nbe_all)

  } else if (nbe_source == "GEOSCHEM") {

      # let the user know this might take some time
      print("Loading processed NBE fields for subsequent sub-setting ...")

      # check which file prefix we are using today
      # list all available files which we will then search
      avail_files = list.files(path_to_nbe,full.names=TRUE)
      prefix = "geoschem_nbe_" # (.)* wildcard characters for unix standard c_gls*_

      # timing information on the number of day in a month
      month_days = rep(31,length.out=12)
      month_days[2] = 28 ; month_days[c(4,6,9,11)] = 30

      lat_done = FALSE ; missing_years = 0 ; keepers = 0 ; yrs = 1
      # loop for year here
      for (yr in seq(1, length(years_to_load))) {
           print(paste("... ",round((yr/length(years_to_load))*100,0),"% completed ",Sys.time(),sep=""))

           # first check how many files we have
           if (yr == 1) {
               nsteps = 0
               for (yrr in seq(1, length(years_to_load))) {
                    # create the prefix to the files we will want for a given year
                    input_file_1=paste(prefix,years_to_load[yrr],sep="")
                    # then check whether this pattern is found in the available files
                    this_year = grepl(input_file_1, avail_files) ; this_year = which(this_year == TRUE)
                    # if we have at least one timestep for this year then we have some information otherwise it is missing!
                    if (length(this_year) > 0) {
                        keepers = keepers+1 ; nsteps = max(nsteps,length(this_year))
                    } else {
                        missing_years = append(missing_years,years_to_load[yrr])
                    }
               } # loop through possible years
               rm(yrr)
           } # first year?

           # open processed files
           input_file_1=paste(prefix,years_to_load[yr],sep="")
           # then check whether this pattern is found in the available files
           this_year = avail_files[grepl(input_file_1, avail_files)]

           if (length(this_year) > 0) {

               # Strip out the time component entirely
               tmp = unlist(strsplit(this_year,input_file_1))
               tmp = tmp[seq(2,length(tmp),2)]
               # The re-order the files correctly
               this_year = this_year[order(tmp)]

               # now loop through the available files for the current year
               for (t in seq(1, length(this_year))) {

                    # open the file
                    data1 = nc_open(this_year[t])

                    # get timing variable
                    tmp = strsplit(this_year[t],input_file_1)[[1]][2]
                    month = as.numeric(substring(tmp,1,2)) #; doy_in = as.numeric(substring(tmp,3,4))
                    doy_in = ncvar_get(data1,"day")
                    # January is correct already, so only adjust if month is >= February
                    if (month > 1) {
                        doy_in = doy_in + sum(month_days[1:(month-1)])
                    }

                    # extract location variables
                    lat_in = ncvar_get(data1, "latitude") ; long_in = ncvar_get(data1, "longitude")
                    # read the NBE observations
                    var1 = ncvar_get(data1, "NBE") # net biome exchange of CO2 (gC/m2/day)
                    # check for error variable
                    var2 = ncvar_get(data1, "NBE_unc") # NBE error estimate(gC/m2/day)

                    # close files after use
                    nc_close(data1)

                    # Turn lat_in / long_in from vectors to arrays
                    lat_in = t(array(lat_in, dim=c(dim(var1)[2],dim(var1)[1])))
                    long_in = array(long_in, dim=c(dim(var1)[1],dim(var1)[2]))

                    # Convert to a raster, assuming standad WGS84 grid
                    var1 = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(var1))
                    var1 = rast(var1, crs = ("+init=epsg:4326"), type="xyz")
                    var2 = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(var2))
                    var2 = rast(var2, crs = ("+init=epsg:4326"), type="xyz")

                    # Create raster with the target crs (technically this bit is not required)
                    target = rast(crs = ("+init=epsg:4326"), ext = ext(var1), resolution = res(var1))
                    # Check whether the target and actual analyses have the same CRS
                    if (compareGeom(var1,target) == FALSE) {
                        # Resample to correct grid
                        var1 = resample(var1, target, method="ngb") ; gc() 
                        var2 = resample(var2, target, method="ngb") ; gc() 
                    }
                    # Extend the extent of the overall grid to the analysis domain
                    var1 = extend(var1,cardamom_ext) ; var2 = extend(var2,cardamom_ext)
                    # Trim the extent of the overall grid to the analysis domain
                    var1 = crop(var1,cardamom_ext) ; var2 = crop(var2,cardamom_ext)
                    var1[which(as.vector(var1) == -9999)] = NA ; var2[which(as.vector(var2) == -9999)] = NA
                    # Adjust spatial resolution of the datasets, this occurs in all cases
                    if (res(var1)[1] != res(cardamom_ext)[1] | res(var1)[2] != res(cardamom_ext)[2]) {

                        # Create raster with the target resolution
                        target = rast(crs = crs(cardamom_ext), ext = ext(cardamom_ext), resolution = res(cardamom_ext))
                        # Resample to correct grid
                        var1 = resample(var1, target, method="bilinear") ; gc() 
                        var2 = resample(var2, target, method="bilinear") ; gc() 

                    } # Aggrgeate to resolution

                    if (lat_done == FALSE) {
                        # extract dimension information for the grid, note the axis switching between raster and actual array
                        xdim = dim(var1)[2] ; ydim = dim(var1)[1]
                        # extract the lat / long information needed
                        long = crds(var1,df=TRUE, na.rm=FALSE)
                        lat  = long$y ; long = long$x
                        # restructure into correct orientation
                        long = array(long, dim=c(xdim,ydim))
                        lat = array(lat, dim=c(xdim,ydim))
                    }
                    # break out from the rasters into arrays which we can manipulate
                    var1 = array(as.vector(unlist(var1)), dim=c(xdim,ydim))
                    var2 = array(as.vector(unlist(var2)), dim=c(xdim,ydim))

                    # remove additional spatial information
                    if (lat_done == FALSE) {
                        # create holding arrays for the nbe information...
                        nbe_hold = as.vector(var1)
                        # ...and its uncertainty information...
                        nbe_unc_hold = as.vector(var2)
                        # ...and timing
                        doy_obs = doy_in
                    } else {
                        # begin populating the various outputs
                        nbe_hold = append(nbe_hold,as.vector(var1))
                        nbe_unc_hold = append(nbe_unc_hold,as.vector(var2))
                        doy_obs = append(doy_obs,doy_in)
                    }

                    # update flag for lat / long load
                    if (lat_done == FALSE) {lat_done = TRUE}

               } # loop through available time steps in the current year

               # keep track of years actually ran
               yrs = yrs + 1
               # clean up allocated memeory
               gc()

           } # is there information for the current year?

      } # year loop

      # Sanity check for NBE
      if (lat_done == FALSE) {stop('No NBE information could be found...')}

      # remove initial value
      missing_years = missing_years[-1]

      # enforce minimum uncertainty value
      nbe_unc_hold[abs(nbe_unc_hold) < 0.01] = 0.01

      # return spatial structure to data
      nbe_out = array(as.vector(nbe_hold), dim=c(xdim,ydim,length(doy_obs)))
      nbe_unc_out = array(as.vector(nbe_unc_hold), dim=c(xdim,ydim,length(doy_obs)))

      # output variables
      nbe_all = list(retrieval_valid = TRUE, nbe_gCm2day = nbe_out, nbe_unc_gCm2day = nbe_unc_out,
                     doy_obs = doy_obs, lat = lat, long = long, missing_years=missing_years)
      # clean up variables
      rm(doy_in,nbe_hold,nbe_unc_hold,nbe_out,doy_obs,lat,long,missing_years) ; gc(reset=TRUE,verbose=FALSE)
      return(nbe_all)

  } else if (nbe_source == "OCO2MIP") {# if GEOSCHEM or other

      # let the user know this might take some time
      print("Loading processed NBE fields for subsequent sub-setting ...")

      # check which file prefix we are using today
      # list all available files which we will then search
      avail_files = list.files(path_to_nbe,full.names=TRUE)
      prefix_est = "EnsMean_gridded" # (.)* wildcard characters for unix standard c_gls*_
      prefix_unc = "EnsStd_gridded" # (.)* wildcard characters for unix standard c_gls*_

      # Only a single file is provided here, from this the correct time period
      # must be extracted

      # Exsure that both files exist
      est_file = avail_files[grepl(prefix_est, avail_files)]
      unc_file = avail_files[grepl(prefix_unc, avail_files)]
      if (length(est_file) != 1 | length(unc_file) != 1) {
          print(paste("Incorrect number of NBE estimte and uncertainty files found ", sep=""))
          print(paste("est_file = ",est_file, sep=""))
          print(paste("unc_file = ",unc_file, sep=""))
          stop()
      }
      # Open both files
      data_est = nc_open(est_file)
      data_unc = nc_open(unc_file)
      # Extract datetime information
      time_in = ncvar_get(data_est, "start_date")
      # Extract estimate and uncertainty information
      nbe_in = ncvar_get(data_est, "land") # gC/m2/year
      nbe_unc_in = ncvar_get(data_unc, "land") # gC/m2/year

      # Extract latitude and longitude
      lat_in = ncvar_get(data_est, "latitude")
      long_in = ncvar_get(data_est, "longitude")
      # Turn lat_in / long_in from vectors to arrays
      lat_in = t(array(lat_in, dim=c(dim(nbe_in)[2],dim(nbe_in)[1])))
      long_in = array(long_in, dim=c(dim(nbe_in)[1],dim(nbe_in)[2]))

      # Close both files
      nc_close(data_est) ; nc_close(data_unc)

      # timing information on the number of day in a month
      month_days = rep(31,length.out=12)
      month_days[2] = 28 ; month_days[c(4,6,9,11)] = 30

      lat_done = FALSE ; missing_years = 0 ; keepers = 0 ; yrs = 1
      # loop for year here
      for (yr in seq(1, length(years_to_load))) {
           print(paste("... ",round((yr/length(years_to_load))*100,0),"% completed ",Sys.time(),sep=""))

           # first check how many files we have
           if (yr == 1) {
               for (yrr in seq(1, length(years_to_load))) {
                    # Check whether the desired year is in the file
                    this_year = which(time_in[1,] == years_to_load[yrr] & time_in[2,] == 1)
                    # There should be a single value
                    if (length(this_year) > 0) {
                        keepers = keepers+1
                    } else {
                        missing_years = append(missing_years,years_to_load[yrr])
                    }
               } # loop through possible years
               rm(yrr)
           } # first year?

           # Check where the start point and end points are for the desired year
           year_start = which(as.numeric(time_in[1,]) == years_to_load[yr] & as.numeric(time_in[2,]) == 1)
           year_end = which(as.numeric(time_in[1,]) == years_to_load[yr] & as.numeric(time_in[2,]) == 12)

           # Assuming we have right year, begin running
           if (length(year_start) > 0) {

               # Ensure these are the first and final values to cover all time steps
               year_start = year_start[1] ; year_end = year_end[length(year_end)]

               # Now loop through the available time steps
               for (t in seq(year_start, year_end)) {

                    # Determine the day of year variable for each time step
                    month = time_in[2,t] ; doy_in = time_in[3,t]
                    # January is correct already, so only adjust if month is >= February
                    if (month > 1) {
                        doy_in = doy_in + sum(month_days[1:(month-1)])
                    }

                    # read the NBE observations
                    var1 = nbe_in[,,t] # Land based net biome exchange of CO2 (gC/m2/yr)
                    # check for error variable
                    var2 = nbe_unc_in[,,t] # NBE error estimate(gC/m2/yr)
                    # Convert units into gC/m2/day
                    var1 = var1 / 365.25 ; var2 = var2 / 365.25

                    # Convert to a raster, assuming standad WGS84 grid
                    var1 = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(var1))
                    var1 = rast(var1, crs = ("+init=epsg:4326"), type="xyz")
                    var2 = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(var2))
                    var2 = rast(var2, crs = ("+init=epsg:4326"), type="xyz")

                    # Create raster with the target crs (technically this bit is not required)
                    target = rast(crs = ("+init=epsg:4326"), ext = ext(var1), resolution = res(var1))
                    # Check whether the target and actual analyses have the same CRS
                    if (compareGeom(var1,target) == FALSE) {
                        # Resample to correct grid
                        var1 = resample(var1, target, method="ngb") ; gc() 
                        var2 = resample(var2, target, method="ngb") ; gc() 
                    }
                    # Extend the extent of the overall grid to the analysis domain
                    var1 = extend(var1,cardamom_ext) ; var2 = extend(var2,cardamom_ext)
                    # Trim the extent of the overall grid to the analysis domain
                    var1 = crop(var1,cardamom_ext) ; var2 = crop(var2,cardamom_ext)
                    var1[which(as.vector(var1) == -9999)] = NA ; var2[which(as.vector(var2) == -9999)] = NA
                    # If this is a gridded analysis and the desired CARDAMOM resolution is coarser than the currently provided then aggregate here.
                    # Despite creation of a cardamom_ext for a site run do not allow aggragation here as tis will damage the fine resolution datasets
                    if (res(var1)[1] != res(cardamom_ext)[1] | res(var1)[2] != res(cardamom_ext)[2]) {

                        # Create raster with the target resolution
                        target = rast(crs = crs(cardamom_ext), ext = ext(cardamom_ext), resolution = res(cardamom_ext))
                        # Resample to correct grid
                        var1 = resample(var1, target, method="bilinear") ; gc()
                        var2 = resample(var2, target, method="bilinear") ; gc() 

                    } # Aggrgeate to resolution

                    if (lat_done == FALSE) {
                        # extract dimension information for the grid, note the axis switching between raster and actual array
                        xdim = dim(var1)[2] ; ydim = dim(var1)[1]
                        # extract the lat / long information needed
                        long = crds(var1,df=TRUE, na.rm=FALSE)
                        lat  = long$y ; long = long$x
                        # restructure into correct orientation
                        long = array(long, dim=c(xdim,ydim))
                        lat = array(lat, dim=c(xdim,ydim))
                    }
                    # break out from the rasters into arrays which we can manipulate
                    var1 = array(as.vector(unlist(var1)), dim=c(xdim,ydim))
                    var2 = array(as.vector(unlist(var2)), dim=c(xdim,ydim))

                    # remove additional spatial information
                    if (lat_done == FALSE) {
                        # create holding arrays for the nbe information...
                        nbe_hold = as.vector(var1)
                        # ...and its uncertainty information...
                        nbe_unc_hold = as.vector(var2)
                        # ...and timing
                        doy_obs = doy_in
                    } else {
                        # begin populating the various outputs
                        nbe_hold = append(nbe_hold,as.vector(var1))
                        nbe_unc_hold = append(nbe_unc_hold,as.vector(var2))
                        doy_obs = append(doy_obs,doy_in)
                    }

                    # update flag for lat / long load
                    if (lat_done == FALSE) {lat_done = TRUE}

               } # loop through available time steps in the current year

               # keep track of years actually ran
               yrs = yrs + 1
               # clean up allocated memeory
               gc()

           } # is there information for the current year?

      } # year loop

      # Sanity check for NBE
      if (lat_done == FALSE) {stop('No NBE information could be found...')}

      # remove initial value
      missing_years = missing_years[-1]

      # enforce minimum uncertainty value
      nbe_unc_hold[abs(nbe_unc_hold) < 0.01] = 0.01

      # return spatial structure to data
      nbe_out = array(as.vector(nbe_hold), dim=c(xdim,ydim,length(doy_obs)))
      nbe_unc_out = array(as.vector(nbe_unc_hold), dim=c(xdim,ydim,length(doy_obs)))

      # output variables
      nbe_all = list(retrieval_valid = TRUE, nbe_gCm2day = nbe_out, nbe_unc_gCm2day = nbe_unc_out,
                     doy_obs = doy_obs, lat = lat, long = long, missing_years=missing_years)
      # clean up variables
      rm(doy_in,nbe_hold,nbe_unc_hold,nbe_out,doy_obs,lat,long,missing_years) ; gc(reset=TRUE,verbose=FALSE)
      return(nbe_all)

  } # select correct NBE estimate

} # function end load_nbe_fields_for_extraction

## Use byte compile
load_nbe_fields_for_extraction<-cmpfun(load_nbe_fields_for_extraction)
