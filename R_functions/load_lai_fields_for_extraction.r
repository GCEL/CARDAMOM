
###
## Function to load leaf area index from global databases
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_lai_fields_for_extraction<-function(latlon_in,lai_source,years_to_load) {

  if (lai_source == "MODIS") {

    # let the user know this might take some time
    print("Loading processed lai fields for subsequent sub-setting ...")

    lat_done = FALSE ; missing_years=0 ; keepers=0 ; yrs=1
    # loop for year here
    for (yr in seq(1, length(years_to_load))) {
      print(paste("... ",round((yr/length(years_to_load))*100,0),"% completed ",Sys.time(),sep=""))

      if (yr == 1) {
        # first check how many files we have
        for (yrr in seq(1, length(years_to_load))) {
          # open processed modis files
          input_file_1=paste(path_to_lai,"modis_lai_with_lat_long_",years_to_load[yrr],".nc",sep="")
          if (file.exists(input_file_1) == TRUE) {keepers=keepers+1} else {missing_years=append(missing_years,years_to_load[yrr])}
        }
      }
      # open processed modis files
      input_file_1=paste(path_to_lai,"modis_lai_with_lat_long_",years_to_load[yr],".nc",sep="")

      # check to see if file exists if it does then we read it in, if not then we assume its a year we don't have data for and move on
      if (file.exists(input_file_1) == TRUE) {
        # open the file
        data1=nc_open(input_file_1)

        # get timing variable
        doy_in=ncvar_get(data1,"doy")

        # extract location variables
        if (lat_done == FALSE) {
          lat=ncvar_get(data1, "lat") ; long=ncvar_get(data1, "long")
          # restrict the spatial extent based on latlong ranges provided
          remove_lat=intersect(which(lat < (max(latlon_in[,1])+2)),which(lat > (min(latlon_in[,1])-2)))
          remove_long=intersect(which(long < (max(latlon_in[,2])+2)),which(long > (min(latlon_in[,2])-2)))
          # now find common where out in both contexts
          remove_lat=intersect(remove_lat,remove_long)
          # update both variables because of common matrix
          remove_long=remove_lat
          # adjust for matrix rather than vector arrangement
          remove_lat=remove_lat/dim(lat)[1]
          remove_long=(remove_long-(floor(remove_lat)*dim(lat)[1]))+1
          remove_lat=ceiling(remove_lat)
          # update new dimensions
          lat_dim=length(min(remove_lat):max(remove_lat)) ; long_dim=length(min(remove_long):max(remove_long))
          lat=lat[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat)] ; long=long[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat)]
        }

        # read the modis lai drivers
        var1=ncvar_get(data1, "modis_lai")
        # reduce spatial cover to the desired area only
        var1=var1[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat),]
        # set actual missing data to -9999
        var1[which(is.na(as.vector(var1)))] = -9999

        # close files after use
        nc_close(data1)

        # remove additional spatial information
        if (lat_done == FALSE) {
          lai_hold=array(NA, dim=c(long_dim*lat_dim*48,keepers))
          lai_hold[1:length(as.vector(var1)),yrs]=as.vector(var1)
          doy_out=doy_in
        } else {
          lai_hold[1:length(as.vector(var1)),yrs]=as.vector(var1)
          doy_out=append(doy_out,doy_in)
        }

        # update flag for lat / long load
        if (lat_done == FALSE) {lat_done = TRUE}
        # keep track of years actually ran
        yrs=yrs+1
        # clean up allocated memeory
        rm(var1) ; gc()
      } # end of does file exist

    } # year loop

    # Sanity check for LAI
    if (lat_done == FALSE) {stop('No LAI information could be found...')}

    # remove initial value
    missing_years=missing_years[-1]

    # check which ones are NA because I made them up
    not_na=is.na(as.vector(lai_hold))
    not_na=which(not_na == FALSE)

    # now remove the ones that are actual missing data
    lai_hold[which(as.vector(lai_hold) == -9999)] = NA
    # return spatial structure to data
    lai_out=array(as.vector(lai_hold)[not_na], dim=c(long_dim,lat_dim,length(doy_out)))
    lai_unc_out = array(-9999, dim=dim(lai_out))
    # output variables
    lai_all = list(lai_all = lai_out, lai_unc_all = lai_unc_out,
                  doy_obs = doy_out, lat = lat, long = long, missing_years=missing_years)
    # clean up variables
    rm(doy_in,lai_hold,not_na,lai_out,doy_out,lat,long,missing_years) ; gc(reset=TRUE,verbose=FALSE)
    return(lai_all)

  } else if (lai_source == "COPERNICUS") {

    # let the user know this might take some time
    print("Loading processed lai fields for subsequent sub-setting ...")

    # check which file prefix we are using today
    # list all available files which we will then search
    avail_files = list.files(path_to_lai,full.names=TRUE)
    prefix = "c_gls(.)*_" # (.)* wildcard characters for unix standard c_gls*_
#    if (length(which(grepl("c_gls_LAI_", avail_files))) > 0) {
#        prefix = "c_gls_LAI_"
#    } else {
#        prefix = "c_gls_LAI300_"
#    }

    # timing information on the number of day in a month
    month_days = rep(31,length.out=12)
    month_days[2] = 28 ; month_days[c(4,6,9,11)] = 30

    lat_done = FALSE ; missing_years = 0 ; keepers = 0 ; yrs = 1
    # loop for year here
    for (yr in seq(1, length(years_to_load))) {
      print(paste("... ",round((yr/length(years_to_load))*100,0),"% completed ",Sys.time(),sep=""))

      # first check how many files we have
      if (yr == 1) {
        # list all available files which we will then search
#        avail_files = list.files(path_to_lai,full.names=TRUE)
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
      # Strip out the time component entirely
      tmp = unlist(strsplit(this_year,input_file_1))
      tmp = tmp[seq(2,length(tmp),2)]
      # The re-order the files correctly
      this_year = this_year[order(tmp)]

      if (length(this_year) > 0) {
          # now loop through the available files for the current year
          for (t in seq(1, length(this_year))) {

            # open the file
            data1 = nc_open(this_year[t])

            # get timing variable
            tmp = strsplit(this_year[t],input_file_1)[[1]][2]
            month = as.numeric(substring(tmp,1,2)) ; doy_in = as.numeric(substring(tmp,3,4))
            # January is correct already, so only adjust if month is >= February
            if (month > 1) {
              doy_in = doy_in + sum(month_days[1:(month-1)])
            }

            # extract location variables
            if (lat_done == FALSE) {
                lat = ncvar_get(data1, "lat") ; long = ncvar_get(data1, "lon")
                lat = array(rev(lat), dim=c(length(lat),length(long))) ; lat = t(lat)
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

            # read the LAI observations
            var1 = ncvar_get(data1, "LAI") # leaf area index (m2/m2)
            # check for erro variable
            if (length(which(grepl("LAI_ERR",names(data1$var)) == TRUE)) > 0) {
                var2 = ncvar_get(data1, "LAI_ERR") # standard error (m2/m2)
            } else if (length(which(grepl("RMSE",names(data1$var)) == TRUE)) > 0) {
                var2 = ncvar_get(data1, "RMSE") # standard error (m2/m2)
            } else {
                stop("LAI error variable cannot be found for copernicus...")
            }
            # re-structure to matching orientation with the lat / long information
            var1 = var1[,dim(var1)[2]:1] ; var2 = var2[,dim(var2)[2]:1]
            # reduce spatial cover to the desired area only
            var1 = var1[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat)]
            var2 = var2[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat)]
            # set actual missing data to -9999
            var1[which(is.na(as.vector(var1)))] = -9999
            var2[which(is.na(as.vector(var2)))] = -9999

            # close files after use
            nc_close(data1)

            # remove additional spatial information
            if (lat_done == FALSE) {
              # create holding arrays for the lai information...
              lai_hold = array(NA, dim=c(long_dim*lat_dim,keepers*nsteps))
              lai_hold[1:length(as.vector(var1)),(t+((yrs-1)*nsteps))] = as.vector(var1)
              # ...and its uncertainty information...
              lai_unc_hold = array(NA, dim=c(long_dim*lat_dim,keepers*nsteps))
              lai_unc_hold[1:length(as.vector(var2)),(t+((yrs-1)*nsteps))] = as.vector(var2)
              # ...and timing
              doy_out = doy_in
            } else {
              # begin populating the various outputs
              lai_hold[1:length(as.vector(var1)),(t+((yrs-1)*nsteps))] = as.vector(var1)
              lai_unc_hold[1:length(as.vector(var2)),(t+((yrs-1)*nsteps))] = as.vector(var2)
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

    # Sanity check for LAI
    if (lat_done == FALSE) {stop('No LAI information could be found...')}

    # remove initial value
    missing_years = missing_years[-1]

    # check which ones are NA because I made them up
    not_na = is.na(as.vector(lai_hold))
    not_na = which(not_na == FALSE)

    filter = as.vector(lai_hold) == -9999 | as.vector(lai_unc_hold) == -9999
    # now remove the ones that are actual missing data
    lai_hold[filter] = NA ; lai_unc_hold[filter] = NA
    # return spatial structure to data
    lai_out = array(as.vector(lai_hold)[not_na], dim=c(long_dim,lat_dim,length(doy_out)))
    lai_unc_out = array(as.vector(lai_unc_hold)[not_na], dim=c(long_dim,lat_dim,length(doy_out)))

    # output variables
    lai_all = list(lai_all = lai_out, lai_unc_all = lai_unc_out,
                  doy_obs = doy_out, lat = lat, long = long, missing_years=missing_years)
    # clean up variables
    rm(doy_in,lai_hold,lai_unc_hold,not_na,lai_out,doy_out,lat,long,missing_years) ; gc(reset=TRUE,verbose=FALSE)
    return(lai_all)

  } # if MODIS

} # function end
## Use byte compile
load_lai_fields_for_extraction<-cmpfun(load_lai_fields_for_extraction)
