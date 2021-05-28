###
## Function deals with the extraction of relevant forest clearance and growth information for UK forestry only!
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

extract_forestry_information<-function(timestep_days,spatial_type,resolution,grid_type,latlon_in,forest_all,start_year,end_year,ctessel_pft_in) {

  # extract information from all maps if there is any
  if (length(dim(forest_all$planting_year)) > 1) {

    # find locations information for UK forest commission
    # convert input data long to conform to what we need
    check1=which(forest_all$planting_long > 180) ; if (length(check1) > 0) { forest_all$planting_long[check1]=forest_all$planting_long[check1]-360 }
    # find the nearest location
    output=closest2d(1,forest_all$planting_lat,forest_all$planting_long,latlon_in[1],latlon_in[2],2)
    i1=unlist(output)[1] ; j1=unlist(output)[2]

    if (forest_all$planting_year[i1,j1] > 0) {
      age=as.numeric(start_year)-as.numeric(forest_all$planting_year[i1,j1])
      yield_class=forest_all$planting_yield[i1,j1]
      ctessel_pft=forest_all$planting_pft[i1,j1]
    } else {
      # this location does not have forest commission informatin
      age=-9999
      yield_class=-9999
      ctessel_pft=ctessel_pft_in
    }

  } else {

    # this location does not have forest commission informatin
    age=-9999
    yield_class=-9999
    ctessel_pft=ctessel_pft_in

  }

  # next move on to working out when there is any deforestation
  # convert input data long to conform to what we need
  check1=which(forest_all$loss_long > 180) ; if (length(check1) > 0) { forest_all$loss_long[check1]=forest_all$loss_long[check1]-360 }

  # find the nearest location
  output=closest2d(1,forest_all$loss_lat,forest_all$loss_long,latlon_in[1],latlon_in[2],2)
  i1=unlist(output)[1] ; j1=unlist(output)[2]

  # decide between data types available
  if (length(dim(forest_all$year_of_loss)) == 0 & length(dim(forest_all$loss_fraction)) == 3) {

    ## using default GFW assumptions

    # now determine if the site has any actual information
    if (max(forest_all$loss_fraction[i1,j1,],na.rm=TRUE) > 0) {
      # we will not insert the deforestation event into the correct location
      doy_out = 0
      for (i in seq(1, length(years_to_do))) {
        # is current year a leap or not
        nos_days = 365
        mod=as.numeric(years_to_do[i])-round((as.numeric(years_to_do[i])/4))*4
        if (mod == 0) {
          nos_days = 366
          mod=as.numeric(years_to_do[i])-round((as.numeric(years_to_do[i])/100))*100
          if (mod == 0) {
            nos_days  = 365
            mod=as.numeric(years_to_do[i])-round((as.numeric(years_to_do[i])/400))*400
            if (mod == 0) {
              nos_days  = 366
            }
          }
        }
        # count up days needed
        doy_out = append(doy_out,1:nos_days)
      }
      doy_out = doy_out[-1]

      # declare output variable
      #deforestation = array(0, dim=length(doy_out))
      deforestation = rep(0, times=length(doy_out))
      ## Assumption used in UK forestry
      #	    # find locations of beginnings of year
      #	    start_of_years=which(doy_out == 120)
      #	    # which year is the one in which deforestation occurs?
      #	    # then find the appropriate beginning of a year and make deforestation
      #            for (aa in seq(1,length(forest_all$year_of_loss))) {
      #	         deforestation[start_of_years[which(as.numeric(years_to_do) == forest_all$year_of_loss[aa])]] = forest_all$loss_fraction[i1,j1,aa]
      #            }
      ## normal assumption
      start_of_years = which(doy_out == 1)
      # which year is the one in which deforestation occurs?
      # then find the appropriate beginning of a year and make deforestation
      for (aa in seq(1,length(forest_all$year_of_loss))) {
        start_point = start_of_years[which(as.numeric(years_to_do) == forest_all$year_of_loss[aa])]
        end_point = start_point + 364
        deforestation[start_point:end_point] = (forest_all$loss_fraction[i1,j1,aa]) / 364
      }

      #	    if (length(timestep_days) == 1 & timestep_days[1] == 1) {
      # well actually we do nothing

      #	    } else {
      # generally this now deals with time steps which are not daily.
      # However if not monthly special case
      if (length(timestep_days) == 1) {
        run_day_selector=seq(1,length(deforestation),timestep_days)
        timestep_days=rep(timestep_days, length.out=length(deforestation))
      }

      print("...calculating weekly / monthly averages for deforestation info")
      # determine the actual daily positions
      run_day_selector = cumsum(timestep_days)
      # create needed variables
      deforestation_agg = array(NA,dim=length(run_day_selector))
      for (y in seq(1,length(run_day_selector))) {
        pick = deforestation[(run_day_selector[y]-timestep_days[y]):run_day_selector[y]]
        deforestation_agg[y] = sum(pick[which(pick != -9999)],na.rm=TRUE)
        # having picked from this period, ensure no overlap by clearing it!
        deforestation[(run_day_selector[y]-timestep_days[y]):run_day_selector[y]] = -9999
      }
      deforestation_agg[which(is.na(deforestation_agg))] = -9999
      # update with new output information
      deforestation = deforestation_agg
      # clean up
      rm(deforestation_agg) ; gc()
      #	    } # monthly aggregation etc

    } else {

      # else we will go with no deforestation at all
      deforestation = 0

    } # missing data condition

  } else {

    ## using UK only planting information and dominant clearance behaviour

    # now determine if the site has any actual information
    deforestation_date=forest_all$year_of_loss[i1,j1]
    if (deforestation_date != -9999) {deforestation_date=deforestation_date+2000}

    # actually in preference to this data set if the forest management location suggest that this is a negative value then we will use this instead
    # i.e. a negative value from the forestry commission map means that is was planted after the start of the simulation
    # also make the age unknown now
    # BASED on information from the NFI the average number of years between clearance and restocking is 2 years.
    if (age < 1 & age != -9999) {deforestation_date=as.numeric(start_year)-(age+2) ; age=-9999}

    if (deforestation_date != -9999) {
      # we will not insert the deforestation event into the correct location
      doy_out=0
      for (i in seq(1, length(years_to_do))) {
        # is current year a leap or not
        nos_days = 365
        mod=as.numeric(years_to_do[i])-round((as.numeric(years_to_do[i])/4))*4
        if (mod == 0) {
          nos_days = 366
          mod=as.numeric(years_to_do[i])-round((as.numeric(years_to_do[i])/100))*100
          if (mod == 0) {
            nos_days  = 365
            mod=as.numeric(years_to_do[i])-round((as.numeric(years_to_do[i])/400))*400
            if (mod == 0) {
              nos_days  = 366
            }
          }
        }
        # count up days needed
        doy_out=append(doy_out,1:nos_days)
      }
      doy_out=doy_out[-1]

      # declare output variable
      deforestation=array(0, dim=length(doy_out))
      # find locations of beginnings of year
      start_of_years=which(doy_out == 120)
      # which year is the one in which deforestation occurs?
      # then find the appropriate beginning of a year and make deforestation
      deforestation[start_of_years[which(as.numeric(years_to_do) == deforestation_date)]] = 1.0

      if (length(timestep_days) == 1 & timestep_days[1] == 1) {
        # well actually we do nothing
      } else {
        # generally this now deals with time steps which are not daily.
        # However if not monthly special case
        if (length(timestep_days) == 1) {
          run_day_selector=seq(1,length(deforestation),timestep_days)
          timestep_days=rep(timestep_days, length.out=length(deforestation))
        }
        print("...calculating monthly averages for deforestation info")

        # determine the actual daily positions
        run_day_selector=cumsum(timestep_days)
        # create needed variables
        deforestation_agg=array(NA,dim=length(run_day_selector))
        for (y in seq(1,length(run_day_selector))) {
          pick=deforestation[(run_day_selector[y]-timestep_days[y]):run_day_selector[y]]
          deforestation_agg[y]=max(pick[which(pick != -9999)],na.rm=TRUE)
          # having picked from this period, ensure no overlap by clearing it!
          deforestation[(run_day_selector[y]-timestep_days[y]):run_day_selector[y]]=-9999
        }
        deforestation_agg[which(is.na(deforestation_agg))] = -9999
        # update with new output information
        deforestation=deforestation_agg
        # clean up
        rm(deforestation_agg) ; gc()
      } # monthly aggregation etc

    } else {

      # else we will go with no deforestation at all
      deforestation=0

    } # missing data condition

  } # type of information

  # return time series and updates pft information
  return(list(ctessel_pft=ctessel_pft,deforestation=deforestation,yield_class=yield_class,age=age))

} # end function
