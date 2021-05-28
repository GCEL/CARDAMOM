
###
## Function to load met data from global field ECMWF data
## subsequently extracted in extract_met_drivers.txt
###

load_burnt_area_fields_for_extraction<-function(latlon_in,burnt_area_source,path_to_burnt_area,start_year,end_year) {

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
                  if (lat_done == FALSE) {lat=ncvar_get(data1, "latitude") ; long=ncvar_get(data1, "longitude")}

                  # read the burnt fraction estimate (units are 0-1)
                  var1 = ncvar_get(data1, "BurnedFraction")
                  # set actual missing data to 0 as missing data is actually no fire
                  var1[which(is.na(as.vector(var1)))] = 0
                  # get time information (month in this case)
                  var2 = ncvar_get(data1, "time") ; time_steps_per_year = 12
                  # approximate doy of the mid-month and allocate fire to that point
                  if (lat_done == FALSE) {
                      doy_obs = floor((var2*(365.25/12))-(365.25/24))
                  } else {
                      doy_obs = append(doy_obs,floor((var2*(365.25/12))-(365.25/24)))
                  }

                  # close files after use
                  nc_close(data1)

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
         burnt_area=array(burnt_area, dim=c(length(long),length(lat),length(doy_obs)))

         # output variables
         return(list(burnt_area=burnt_area,doy_obs=doy_obs,lat=lat,long=long,missing_years=missing_years))

    } else if (burnt_area_source == "MCD64A1") {

         # let the user know this might take some time
         print("Loading MCD64A1 processed burnt area estimates for subsequent sub-setting ...")

         lat_done = FALSE ; missing_years = 0 ; keepers = 0 ; yrs = 1
         # loop for year here
         for (yr in seq(1, length(years_to_do))) {
              print(paste("... ",round((yr/length(years_to_do))*100,0),"% completed ",Sys.time(),sep=""))

              if (yr == 1) {
                  # list the available files
                  available_files = list.files(path_to_burnt_area,full.names=TRUE)
                  # first check how many files we have
                  for (yrr in seq(1, length(years_to_do))) {
                      if (length(which(grepl(years_to_do[yrr],available_files))) > 0) {keepers = keepers+1} else {missing_years = append(missing_years,years_to_do[yrr])}
                  }
              }
              # open processed modis files
              input_file_1 = paste(path_to_burnt_area,"/MCD64A1_",years_to_do[yr],"(.)*.nc",sep="")

              # check to see if file exists if it does then we read it in, if not then we assume its a year we don't have data for and move on
              if (file.exists(input_file_1) == TRUE) {

                  # open the file
                  data1 = nc_open(input_file_1)

                  # extract location variables
                  if (lat_done == FALSE) {
                      lat = ncvar_get(data1, "latitude") ; long = ncvar_get(data1, "longitude")
                      # Restrict the spatial extent based on latlong ranges provided
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

                  # read the burnt fraction estimate (units are 0-1)
                  var1 = ncvar_get(data1, "BurnedFraction")
                  # restrict spatial area
                  var1 = var1[min(remove_long):max(remove_long),min(remove_lat):max(remove_lat),]
                  # set actual missing data to 0 as missing data is actually no fire
                  var1[which(is.na(as.vector(var1)))] = 0

                  # get time information (day of year)
                  var2 = ncvar_get(data1, "time") ; time_steps_per_year = 12
                  # approximate doy of the mid-month and allocate fire to that point
                  if (lat_done == FALSE) {
                      doy_obs = var2
                  } else {
                      doy_obs = append(doy_obs,var2)
                  }

                  # close files after use
                  nc_close(data1)

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
          missing_years = missing_years[-1]

          # clean up variables
          rm(var1,var2) ; gc(reset=TRUE,verbose=FALSE)

          # restructure
          burnt_area=array(burnt_area, dim=c(dim(long)[1],dim(lat)[2],length(doy_obs)))

          # output variables
          return(list(burnt_area=burnt_area,doy_obs=doy_obs,lat=lat,long=long,missing_years=missing_years))

    } else if (burnt_area_source == " " | burnt_area_source == "site_specific"){

	        # Do nothing as this should be read directly from files or not needed

    } else {

	        stop(paste("Burnt area option (",burnt_area_source,") not valid"))

    } # if MPI biomass

} # function end
