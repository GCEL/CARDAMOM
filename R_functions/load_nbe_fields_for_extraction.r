
###
## Function to load Net Biome Exchange for subsquent subsetting
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_nbe_fields_for_extraction<-function(latlon_in,nbe_source,years_to_load) {

  if (nbe_source == "GEOSCHEM") {

      # let the user know this might take some time
      print("Loading processed NBE fields for subsequent sub-setting ...")

      # check which file prefix we are using today
      # list all available files which we will then search
      avail_files = list.files(path_to_nbe,full.names=TRUE)
      if (length(which(grepl("uoe_ln_v1.0.1x1.",avail_files) == TRUE)) > 0) {prefix = "uoe_ln_v1.0.1x1."}
      if (length(which(grepl("uoe_ln_v1.0.4x5.",avail_files) == TRUE)) > 0) {prefix = "uoe_ln_v1.0.4x5."}

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
                    # create the prefix to the files we will want for a given year
                    input_file_1 = paste(prefix,years_to_load[yrr],sep="")
                    # then check whether this pattern is found in the available files
                    this_year = grepl(input_file_1, avail_files) ; this_year = which(this_year == TRUE)
                    # if we have at least one timestep for this year then we have some information otherwise it is missing!
                    if (length(this_year) > 0) {
                         keepers = keepers+1
                    } else {
                         missing_years = append(missing_years,years_to_load[yrr])
                    }
               } # loop through possible years
               rm(yrr)
           } # first year?

           # open processed modis files
           input_file_1=paste(prefix,years_to_load[yr],sep="")
           # then check whether this pattern is found in the available files
           this_year = avail_files[grepl(input_file_1, avail_files)]

           if (length(this_year) > 0) {
               # now loop through the available files for the current year
               nos_days = nos_days_in_year(years_to_load[yr])
               for (t in seq(1, length(this_year))) {

                    # open the file
                    data1 = nc_open(this_year[t])

                    # get timing variable
                    tmp = strsplit(this_year[t],input_file_1)[[1]][2]
                    month = as.numeric(substring(tmp,1,2))
                    day_of_month = ncvar_get(data1, "Day")
                    if (nos_days == 366) {
                        month_days[2] = 29
                    } else {
                        month_days[2] = 28
                    }
                    # January is correct already, so only adjust if month is >= February
                    if (month > 1) {
                        doy_in = day_of_month + sum(month_days[1:(month-1)])
                    } else {
                        doy_in = day_of_month
                    }

                    # extract location variables
                    if (lat_done == FALSE) {
                        lat = ncvar_get(data1, "lat") ; long = ncvar_get(data1, "lon")
                        lat = array(lat, dim=c(length(lat),length(long))) ; lat = t(lat)
                        long = array(long, dim=dim(lat))
                        # restrict the spatial extent based on latlong ranges provided
                        remove_lat = intersect(which(lat < (max(latlon_in[,1])+2)),which(lat > (min(latlon_in[,1])-2)))
                        remove_long = intersect(which(long < (max(latlon_in[,2])+2)),which(long > (min(latlon_in[,2])-2)))
                        # now find common where out in both contexts
                        remove_lat = intersect(remove_lat,remove_long)
                        # update both variables because of common matrix
                        remove_long = remove_lat
                        # adjust for matrix rather than vector arrangement
                        remove_lat = remove_lat/dim(lat)[1]
                        remove_long = (remove_long-(floor(remove_lat)*dim(lat)[1]))+1
                        remove_lat = ceiling(remove_lat)
                        # update new dimensions
                        lat_dim = length(min(remove_lat):max(remove_lat)) ; long_dim = length(min(remove_long):max(remove_long))
                        lat = lat[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat)]
                        long = long[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat)]
                    }

                    # read the observations
                    var1 = ncvar_get(data1, "NBE") # Net Biome Exchange (gC/m2/day) - negative is sink
                    var2 = ncvar_get(data1, "NBE_unc") # Uncertainty (gC/m2/day)
                    # reduce spatial cover to the desired area only
                    var1 = var1[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat),]
                    var2 = var2[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat),]
                    # set actual missing data to -9999
                    var1[is.na(as.vector(var1))] = -9999 ; var2[is.na(as.vector(var2))] = -9999

                    # close files after use
                    nc_close(data1)

                    # remove additional spatial information
                    if (lat_done == FALSE) {
                        # create holding arrays for the lai information...
                        nbe_gCm2day = array(var1, dim=c(long_dim*lat_dim,length(doy_in)))
                        # ...and its uncertainty information...
                        nbe_unc_gCm2day = array(var2, dim=c(long_dim*lat_dim,length(doy_in)))
                        # ...and timing
                        doy_out = doy_in
                    } else {
                        # begin populating the various outputs
                        nbe_gCm2day = cbind(nbe_gCm2day,array(var1,dim=c(long_dim*lat_dim,length(doy_in))))
                        nbe_unc_gCm2day = cbind(nbe_unc_gCm2day,array(var2,dim=c(long_dim*lat_dim,length(doy_in))))
                        doy_out = append(doy_out,doy_in)
                    }

                    # update flag for lat / long load
                    if (lat_done == FALSE) {lat_done = TRUE}

               } # loop through available time steps in the current year

               # keep track of years actually ran
               yrs = yrs + 1
               # clean up allocated memeory
               rm(var1,var2) ; gc()

           } # is there information for the current year?

      } # year loop

      # Sanity check for NBE
      if (lat_done == FALSE) {stop('No NBE information could be found...')}

      # remove initial value
      missing_years = missing_years[-1] ; missing_years = as.numeric(missing_years)
      # now remove the ones that are actual missing data
      nbe_gCm2day[nbe_gCm2day == -9999] = NA ; nbe_unc_gCm2day[nbe_unc_gCm2day == -9999] = NA
      # enforce minimum uncertainty value
      nbe_unc_gCm2day[abs(nbe_unc_gCm2day) < 0.01] = 0.01
      # return spatial structure to data
      nbe_gCm2day = array(as.vector(nbe_gCm2day), dim=c(long_dim,lat_dim,length(doy_out)))
      nbe_unc_gCm2day = array(as.vector(nbe_unc_gCm2day), dim=c(long_dim,lat_dim,length(doy_out)))

      # output variables
      nbe_all = list(nbe_gCm2day = nbe_gCm2day, nbe_unc_gCm2day = nbe_unc_gCm2day,
                     doy_obs = doy_out, lat = lat, long = long, missing_years=missing_years)
      return(nbe_all)

  } # if GEOSCHEM

} # function end
