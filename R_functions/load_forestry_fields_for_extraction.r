load_forestry_fields_for_extraction<-function(latlon_in,forestry_source,years_to_load) {

  UK_forest_hack = FALSE

  if (forestry_source == "combined_dataset") {
      # let the user know this might take some time
      print("Loading processed 1 km forestry fields for subsequent sub-setting ...")

      # open processed modis files
      #input_file_1=paste(path_to_forestry,"UK_forestry_planting.nc",sep="")
      input_file_1=paste(path_to_forestry,"UK_forestry_planting_public_and_private.nc",sep="")
      #input_file_1=paste(path_to_forestry,"UK_forestry_planting_FC_public_and_private.nc",sep="")
      input_file_2=paste(path_to_forestry,"year_of_forest_loss.nc",sep="")
      data1=nc_open(input_file_1) ; data2=nc_open(input_file_2)

      # extract location variables
      lat=ncvar_get(data1, "lat") ; long=ncvar_get(data1, "long")

      # read pft information
      primary_pft=ncvar_get(data1, "primary_pft")
      #	secondary_pft=ncvar_get(data1, "secondary_pft")
      #	tertiary_pft=ncvar_get(data1, "tertiary_pft")
      # read cover to determine the max cover
      primary_cover=ncvar_get(data1, "primary_cover")
      #	secondary_cover=ncvar_get(data1, "secondary_cover")
      #	tertiary_cover=ncvar_get(data1, "tertiary_cover")
      # generate empty space for this
      dims=dim(primary_cover)

      # first make areas with no pft information lost to the analysis
      primary_cover[which(is.na(as.vector(primary_pft)))]=-9999
      # secondary_cover[which(is.na(as.vector(secondary_pft)))]=-9999
      # tertiary_cover[which(is.na(as.vector(tertiary_pft)))]=-9999

      # might not be the most efficient way of doing things but here we go
      #	keep_primary=which(as.vector(primary_cover) > 0 &
      #                    as.vector(primary_cover) > as.vector(secondary_cover) &
      #                    as.vector(primary_cover) > as.vector(tertiary_cover))
      # keep_secondary=which(as.vector(secondary_cover) > 0 &
      #                      as.vector(secondary_cover) > as.vector(primary_cover) &
      #                      as.vector(secondary_cover) > as.vector(tertiary_cover))
      # keep_tertiary=which(as.vector(tertiary_cover) > 0 &
      #                     as.vector(tertiary_cover) > as.vector(primary_cover) &
      #                     as.vector(tertiary_cover) > as.vector(secondary_cover))

      # create the empty spaces with default value -9999 present
      planting_year=array(-9999, dim=prod(dims))
      planting_yield=array(-9999, dim=prod(dims))
      planting_pft=array(-9999, dim=prod(dims))

      # re-read the data
      primary_cover=ncvar_get(data1, "primary_pft")
      #	secondary_cover=ncvar_get(data1, "secondary_pft")
      #	tertiary_cover=ncvar_get(data1, "tertiary_pft")
      # ensure that all missing data are set as -9999 first
      primary_cover[which(is.na(as.vector(primary_cover)))]=-9999
      #	secondary_cover[which(is.na(as.vector(secondary_cover)))]=-9999
      #	tertiary_cover[which(is.na(as.vector(tertiary_cover)))]=-9999
      # now extract the bits we want
      planting_pft=as.vector(primary_cover)
      #	planting_pft[keep_primary]=as.vector(primary_cover)[keep_primary]
      #	planting_pft[keep_secondary]=as.vector(secondary_cover)[keep_secondary]
      #	planting_pft[keep_tertiary]=as.vector(tertiary_cover)[keep_tertiary]

      # read planting year information (could inform time since variables)
      primary_cover=ncvar_get(data1, "primary_plant")
      #	secondary_cover=ncvar_get(data1, "secondary_plant")
      #	tertiary_cover=ncvar_get(data1, "tertiary_plant")
      # ensure that all missing data are set as -9999 first
      primary_cover[which(is.na(as.vector(primary_cover)))]=-9999
      #	secondary_cover[which(is.na(as.vector(secondary_cover)))]=-9999
      #	tertiary_cover[which(is.na(as.vector(tertiary_cover)))]=-9999
      # now allocate the appropriate information
      planting_year=as.vector(primary_cover)
      #	planting_year[keep_primary]=as.vector(primary_cover)[keep_primary]
      #	planting_year[keep_secondary]=as.vector(secondary_cover)[keep_secondary]
      #	planting_year[keep_tertiary]=as.vector(tertiary_cover)[keep_tertiary]

      # read planting year information (could inform time since variables)
      primary_cover=ncvar_get(data1, "primary_yield")
      #	secondary_cover=ncvar_get(data1, "secondary_yield")
      #	tertiary_cover=ncvar_get(data1, "tertiary_yield")
      # ensure that all missing data are set as -9999 first
      primary_cover[which(is.na(as.vector(primary_cover)))]=-9999
      #	secondary_cover[which(is.na(as.vector(secondary_cover)))]=-9999
      #	tertiary_cover[which(is.na(as.vector(tertiary_cover)))]=-9999
      # now allocate the appropriate information
      planting_yield=as.vector(primary_cover)
      #	planting_yield[keep_primary]=as.vector(primary_cover)[keep_primary]
      #	planting_yield[keep_secondary]=as.vector(secondary_cover)[keep_secondary]
      #	planting_yield[keep_tertiary]=as.vector(tertiary_cover)[keep_tertiary]

      # reconstruct the data
      planting_year=array(planting_year, dim=dims)
      planting_yield=array(planting_yield, dim=dims)
      planting_pft=array(planting_pft, dim=dims)
      # keep some spatial information
      planting_lat=lat ; planting_long=long

      # close files after use
      nc_close(data1)

      # begin reading the forest loss information instead
      loss_lat=ncvar_get(data2, "lat") ; loss_long=ncvar_get(data2, "long")
      # read year of forest loss informatin
      year_of_loss=ncvar_get(data2, "yearloss")
      dims=dim(year_of_loss)
      year_of_loss[which(as.vector(is.na(year_of_loss)))] = -9999
      year_of_loss=array(year_of_loss, dim=dims)
      loss_fraction=array(0, dim=dims)
      loss_fraction[which(year_of_loss != -9999)]=1

      # tidy up
      nc_close(data2)

      #rm(primary_cover,secondary_cover,tertiary_cover,dims,lat,long,input_file_1,input_file_2) ; gc(reset=TRUE,verbose=FALSE)

      # output variables
      return(list(planting_lat=planting_lat,planting_long=planting_long,planting_year=planting_year,planting_yield=planting_yield,planting_pft=planting_pft
                 ,loss_lat=loss_lat,loss_long=loss_long,year_of_loss=year_of_loss,loss_fraction=loss_fraction))

  } else if (forestry_source == "GFW" & UK_forest_hack) {

      # let the user know this might take some time
      print("Loading processed 1 km forestry fields for subsequent sub-setting ...")

      # open processed UK FC planting information files
      input_file_1=paste(path_to_forestry,"UK_forestry_planting_public_and_private.nc",sep="")
      data1=nc_open(input_file_1)

      # extract location variables
      lat=ncvar_get(data1, "lat") ; long=ncvar_get(data1, "long")

      # read pft information
      primary_pft=ncvar_get(data1, "primary_pft")

      # read cover to determine the max cover
      primary_cover=ncvar_get(data1, "primary_cover")

      # generate empty space for this
      dims=dim(primary_cover)

      # first make areas with no pft information lost to the analysis
      primary_cover[which(is.na(as.vector(primary_pft)))]=-9999

      # create the empty spaces with default value -9999 present
      planting_year=array(-9999, dim=prod(dims))
      planting_yield=array(-9999, dim=prod(dims))
      planting_pft=array(-9999, dim=prod(dims))

      # re-read the data
      primary_cover=ncvar_get(data1, "primary_pft")
      # ensure that all missing data are set as -9999 first
      primary_cover[which(is.na(as.vector(primary_cover)))]=-9999
      # now extract the bits we want
      planting_pft=as.vector(primary_cover)

      # read planting year information (could inform time since variables)
      primary_cover=ncvar_get(data1, "primary_plant")

      # ensure that all missing data are set as -9999 first
      primary_cover[which(is.na(as.vector(primary_cover)))]=-9999
      # remove information for sites planted after 2000
      primary_cover[which(primary_cover < 0)] = -9999

      # now allocate the appropriate information
      planting_year=as.vector(primary_cover)

      # read planting year information (could inform time since variables)
      primary_cover=ncvar_get(data1, "primary_yield")
      # ensure that all missing data are set as -9999 first
      primary_cover[which(is.na(as.vector(primary_cover)))]=-9999
      # now allocate the appropriate information
      planting_yield=as.vector(primary_cover)

      # reconstruct the data
      planting_year=array(planting_year, dim=dims)
      planting_yield=array(planting_yield, dim=dims)
      planting_pft=array(planting_pft, dim=dims)
      # keep some spatial information
      planting_lat=lat ; planting_long=long

      # close files after use
      nc_close(data1)

      for (yrr in seq(1,length(years_to_load))){

           input_file_2=paste(path_to_forestry,"GFW_forest_loss_",years_to_load[yrr],".nc",sep="")
           data2=nc_open(input_file_2)

           if (yrr == 1) {
               # begin reading the forest loss information instead
               loss_lat=ncvar_get(data2, "latitude") ; loss_long=ncvar_get(data2, "longitude")
               dims=c(dim(loss_lat),length(years_to_load))
               loss_fraction=array(NA, dim=dims)
           }
           # read year of forest loss informatin
           loss_fraction_tmp = ncvar_get(data2, "forest_loss")
           loss_fraction_tmp[which(as.vector(is.na(loss_fraction_tmp)))] = 0
           loss_fraction=array(loss_fraction, dim=dims)
           # place new clearance information into the output array
           loss_fraction[,,yrr]=loss_fraction_tmp

           # tidy up
           nc_close(data2)

      } # looping years

      # output variables
      return(list(planting_lat=planting_lat,planting_long=planting_long,planting_year=planting_year,planting_yield=planting_yield,planting_pft=planting_pft
                 ,loss_lat=loss_lat,loss_long=loss_long,year_of_loss=years_to_load,loss_fraction=loss_fraction))

  } else if (forestry_source == "GFW" & UK_forest_hack == FALSE) {

      # let the user know this might take some time
      print("Loading processed Global Forest Watch clearance information for subsequent sub-setting ...")

      # Set flag
      lat_done = FALSE
      # Loop through years
      for (yrr in seq(1,length(years_to_load))){

          # Create file name # two possible versions of  naming syntax
          input_file_2a = paste(path_to_forestry,"GFW_forest_loss_",years_to_load[yrr],".nc",sep="")
          input_file_2b = paste(path_to_forestry,"tree_cover_loss_fraction",years_to_load[yrr],"(.)*.nc",sep="")
          # Check it exists
          if (any(file.exists(input_file_2a),file.exists(input_file_2b))) {

              # Open the current file
              if (file.exists(input_file_2a)){
                  data2 = nc_open(input_file_2a)
                  loss_variable = "forest_loss"
              } else{
                  data_2 = nc_open(input_file_2b)
                  loss_variable = "tree_cover_loss"
              }

              # Determine  lat / long bounds
              if (length(dim(latlon_in)) > 1) {
                  max_lat=max(latlon_in[,1])+1.0 ; max_long=max(latlon_in[,2])+1.0
                  min_lat=min(latlon_in[,1])-1.0 ; min_long=min(latlon_in[,2])-1.0
              } else {
                  max_lat=max(latlon_in[1])+1.0 ; max_long=max(latlon_in[2])+1.0
                  min_lat=min(latlon_in[1])-1.0 ; min_long=min(latlon_in[2])-1.0
              }

              # Create the grid once only
              if (lat_done == FALSE) {

                  # begin reading the forest loss information instead
                  loss_lat=ncvar_get(data2, "latitude") ; loss_long=ncvar_get(data2, "longitude")
                  keep_lat=which(loss_lat[1,] > min_lat & loss_lat[1,] < max_lat)
                  keep_long=which(loss_long[,1] > min_long & loss_long[,1] < max_long)
                  lat_dim = length(min(keep_lat):max(keep_lat))
                  long_dim = length(min(keep_long):max(keep_long))
                  loss_lat=loss_lat[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
                  loss_long=loss_long[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
                  # rebuild dimensions
                  loss_lat = array(loss_lat, dim=c(long_dim,lat_dim)) ; loss_long = array(loss_long, dim=c(long_dim,lat_dim))
                  dims=c(dim(loss_lat),length(years_to_load))
                  loss_fraction=array(NA, dim=dims)
                  # set flag
                  lat_done = TRUE
              } # lat_done

              # read year of forest loss informatin
              loss_fraction_tmp = ncvar_get(data2, loss_variable)
              loss_fraction_tmp = loss_fraction_tmp[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
              loss_fraction_tmp[which(as.vector(is.na(loss_fraction_tmp)))] = 0
              loss_fraction = array(loss_fraction, dim=dims)
              # place new clearance information into the output array
              loss_fraction[,,yrr] = loss_fraction_tmp

              # tidy up
              nc_close(data2)

          } # File exists
      } # looping years

      # output variables
      return(list(planting_lat=-9999,planting_long=-9999,planting_year=-9999,planting_yield=-9999,planting_pft=-9999
                 ,loss_lat=loss_lat,loss_long=loss_long,year_of_loss=years_to_load,loss_fraction=loss_fraction))

  } else {

      # output variables
      return(list(planting_lat=-9999,planting_long=-9999,planting_year=-9999,planting_yield=-9999,planting_pft=-9999
                 ,loss_lat=-9999,loss_long=-9999,year_of_loss=-9999,loss_fraction=-9999))

  }

} # end of function
