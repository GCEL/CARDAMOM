#########################################################################################
# CARbon DAta MOdel fraMework (CARDAMOM) and DALEC terrestrial ecosystem model suite
# CARDAMOM is a Bayesian model-data fusion software framework. CARDAMOM is used to 
# assimilate observations and ecological theory to retrieve parameters for the 
# DALEC suite of intermediate complexity terrestrial ecosystem models. DALEC can be
# used as a fully integrated component of CARDAMOM or independently. 
# Copyright (C) 2024  University of Edinburgh,
#                     Mathew Williams (mat.williams@ed.ac.uk), 
#                     T. Luke Smallman (t.l.smallman@ed.ac.uk)
# UoE = University of Edinburgh

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# ########## File specific description ##########
# Function to load meteorological forcing from various gridded datasets.
# 
# Author: T. Luke Smallman (02/05/2024)
# Modification History:
# 1) T. Luke Smallman (11/03/2025) 
#    ERA5 workflow updated to read whole year files and to only read required spatial domain rather than whole domain
# 2) 
#########################################################################################

load_met_fields_for_extraction<-function(latlon_in,met_source,modelname,startyear,endyear,spatial_type,cardamom_ext) {
  
  # let the user know this might take some time
  print("Loading global met fields for subsequent sub-setting...")
  # declare timing variables
  t_grid = 0
  
  if (met_source == "site_specific") {
    
    # contruct output
    met_all = list(site_specific=TRUE, wheat = c(1:dim(latlon_in)[1]))
    
  } else {
    
    if (met_source == "trendy") {
      
      # declare variable ids needed to select files / infile variables
      varid = c("dswrf","tmx","pre","vpd","tmn","wsp")
      infile_varid = c("dswrf","tmx","pre","vpd","tmn","wsp")
      
      # open first ecmwf file to extract needed information
      input_file_1 = paste(path_to_met_source,"/",varid[1],"_",startyear,"_monthly.nc",sep="")
      data1 = nc_open(input_file_1)
      
      # Get timing information
      steps_in_day = 1/30.4375
      
      # extract location variables
      lat = ncvar_get(data1, "lat") ; long = ncvar_get(data1, "lon")
      #            # convert input data long to conform to what we need
      #            check1 = which(long > 180) ; if (length(check1) > 0) { long[check1] = long[check1]-360 }
      
      # expand the one directional values here into 2 directional
      lat_dim = length(lat) ; long_dim = length(long)
      long = array(long,dim=c(long_dim,lat_dim))
      lat = array(lat,dim=c(lat_dim,long_dim)) ; lat = t(lat)
      
      # Read in an example of the first variable - this should be Shortwave radiation
      tmp1 = ncvar_get(data1,infile_varid[1])
      # We only need the first time step for this
      tmp1 = tmp1[,,1]
      
      # tidy
      nc_close(data1)
      
      # Unlike ERA5, or rather the standard CARDAMOM meterological input format,
      # for trendy we read in the whole domain at all times. So the xy_bounds variable 
      # has a default filling, where 1 is the start point and -1 is the read until end of file.
      xy_bounds = c(1,-1,1,-1)
      
      # If we are using the GSI model we need the 30 days (or month) before the start date of the simulation, so we need to check if we have this information
      # NOTE that this section of code is duplicated for each of the available datasets because of differences in storage and file name
      extra_year = FALSE
      present = 0
      for (lag in seq(1,length(varid))) {
        if (varid[lag] != "") {
          input_file_1 = paste(path_to_met_source,"/",varid[lag],"_",as.character(as.numeric(startyear)-1),"_monthly.nc",sep="")
        }
        # check whether the files exist or not
        if (file.exists(input_file_1)) {present = present+1}
      }
      # if all the files exist then we will use them
      if (present == length(which(varid != ""))) {extra_year=TRUE}
      
    } else if (met_source == "ERA" | met_source == "isimip3a") {
      
      # declare variable ids needed to select files / infile variables
      varid = c("sw_radiation_daily_mean","airt_daily_max","precipitation_daily_mean","vpd_daily_mean","airt_daily_min","wind_spd_daily_mean")
      infile_varid = c("daily_swrad","airt_max","daily_precip","vpd_mean","airt_min","wind_spd")
      
      # Open the ERA5 file for the shortwave radiation.
      # Could have been a different variable or file, but this is easy to work with.
      input_file_1 = paste(path_to_met_source,varid[1],"_",startyear,".nc",sep="")
      data1 = nc_open(input_file_1)
      
      # get timing variable
      doy = ncvar_get(data1, "doy")
      steps_in_day = 1 / (doy[2]-doy[1])
      # extract location variables
      lat = ncvar_get(data1, "lat") ; long = ncvar_get(data1, "lon")
      
      # Determine the spatial domain we want to read for the cardamom domain
      # NOTE: this code assumes a regular grid only
      xy_bounds = unlist(crs(cardamom_ext, describe=TRUE)$extent) # first convert the project domain back into lat / long values
      xy_bounds[1] = max(1,min(which(long[,1] >= xy_bounds[1]))) ; xy_bounds[2] = min(dim(long)[1],max(which(long[,1] <= xy_bounds[2])))
      xy_bounds[3] = max(1,min(which(lat[1,] >= xy_bounds[3])))  ; xy_bounds[4] = min(dim(lat)[2],max(which(lat[1,] <= xy_bounds[4])))
      # The netcdf function requires a starting point and a number of reads.
      # The starting points are xy_bounds[1] and xy_bounds[3], 
      # the counts need to be filled into xy_bounds[2] and xy_bounds[4]
      xy_bounds[2] = length(c(xy_bounds[1]:xy_bounds[2]))
      xy_bounds[4] = length(c(xy_bounds[3]:xy_bounds[4]))
      
      # Re-read the location information to ensure match with the variable
      long = ncvar_get(data1, "lon", start=c(xy_bounds[1],xy_bounds[3]), count=c(xy_bounds[2],xy_bounds[4]))
      lat = ncvar_get(data1, "lat", start=c(xy_bounds[1],xy_bounds[3]), count=c(xy_bounds[2],xy_bounds[4])) 
      
      # Read in an example of the first variable - this should be Shortwave radiation
      tmp1 = ncvar_get(data1,infile_varid[1], 
                       start=c(xy_bounds[1],xy_bounds[3],1), 
                       count=c(xy_bounds[2],xy_bounds[4],1))
      # tidy
      nc_close(data1)
      
      # Several model structures require at least the preceeding months meteorology
      # for rolling averages. Check here if we have this information.
      extra_year = FALSE ; present = 0
      for (lag in seq(1,length(varid))) {
        input_file_1 = paste(path_to_met_source,varid[lag],"_",as.character(as.numeric(startyear)-1),".nc",sep="")
        # check whether the files exist or not
        if (file.exists(input_file_1)) {present = present+1}
      }
      # if all the files exist then we will use them
      if (present == length(varid)) {extra_year=TRUE}
      
    } # met_source
    
    # Convert to a raster, assuming standad WGS84 grid
    # This depends on the lat / long / tmp1 spatially matching each other AND
    # latitude ranging -90/90 and longitude ranging -180/180 degrees
    tmp1 = data.frame(x = as.vector(long), y = as.vector(lat), z = as.vector(tmp1))
    tmp1 = rast(tmp1, crs = ("epsg:4326"), type="xyz")
    
    # Extract the epsg from the file
    epsg = crs(tmp1, describe = TRUE)$code
    # If we have an epsg then we want to know if it differs from the one desired by the analysis
    if (epsg != gsub("epsg:","",cardamom_grid_type)) {
      # Ensure that the extent of the input object is consistent 
      # with the possible extent of the selected epsg
      tmp1 = crop(tmp1, ext(unlist(crs(cardamom_grid_type, describe=TRUE)$extent)))
      # If it does not match we need to reproject it
      tmp1 = project(tmp1, cardamom_grid_type, method="near", align = FALSE) ; gc()
    }
    # Extend the extent of the overall grid to the analysis domain
    tmp1 = extend(tmp1,cardamom_ext)
    # Trim the extent of the overall grid to the analysis domain
    tmp1 = crop(tmp1,cardamom_ext)
    
    # Set error flags to NA
    tmp1[which(as.vector(tmp1) == -9999)] = NA
    # We assume where that the first variable is shortwave radiation or
    # other positive definite variable. Thus, locations which are < 0
    # Are assumed to be non-valid locations
    tmp1[tmp1 < 0] = NA
    
    # Match resolutions
    if (res(tmp1)[1] != res(cardamom_ext)[1] | res(tmp1)[2] != res(cardamom_ext)[2]) {
      # Resample to correct grid
      tmp1 = resample(tmp1, cardamom_ext, method="average") ; gc() 
    } # Aggrgeate to resolution
    
    # Filter through the reduced dataset for the specific locations
    # NOTE: cell number (or pixel number in the grid), this variable is not flipped to human viewing
    # whereas the i_loc and j_loc have been
    output = closest2d_tif(tmp1,latlon_in[,1],latlon_in[,2])
    output_n = output$n_loc
    
    # Extract dimension information for the aggregated grid
    # NOTE the axis switching between raster and actual array
    long_dim = dim(tmp1)[2] ; lat_dim = dim(tmp1)[1]
    
    # Extract from raster for easier manipulation
    tmp1 = values(tmp1)
    # Filter for those points which match out desired land locations
    tmp1 = tmp1[output_n]
    # Of these location we now want to track the locations which contain valid
    # meteorology
    wheat_from_chaff = output_n[which(is.finite(tmp1))]
    rm(output_n,tmp1)
    
    # define timing variables
    years_to_load = as.numeric(startyear):as.numeric(endyear)
    if (extra_year) {
      load_years = c(as.character(as.numeric(startyear)-1),as.numeric(startyear):as.numeric(endyear))
    } else {
      load_years = as.numeric(startyear):as.numeric(endyear)
    }
    print("have finished determining met extraction parameters, now for the main load")
    # extract the data into lists which we will contruct.
    if (use_parallel) {
      
      # NOTE: that the use of mclapply() is due to reported improved efficiency over creating a virtual cluster.
      # However, mclapply does not (at the time of typing) work on Windows, i.e. Linux and Mac only
      cl <- min(length(load_years),numWorkers)
      if (varid[1] != "") {
        swrad_out_list=mclapply(load_years,FUN=load_met_function,varid=varid[1],infile_varid=infile_varid[1],
                                spatial_type=spatial_type,cardamom_ext=cardamom_ext,path_to_met_source=path_to_met_source,
                                met_source=met_source,wheat=wheat_from_chaff,xy_bounds = xy_bounds,mc.cores = cl)
      }
      print("...met load 17 %")
      if (varid[2] != "") {
        maxt_out_list=mclapply(load_years,FUN=load_met_function,varid=varid[2],infile_varid=infile_varid[2],
                               spatial_type=spatial_type,cardamom_ext=cardamom_ext,path_to_met_source=path_to_met_source,
                               met_source=met_source,wheat=wheat_from_chaff,xy_bounds = xy_bounds,mc.cores = cl)
      }
      print("...met load 33 %")
      if (varid[3] != "") {
        precip_out_list=mclapply(load_years,FUN=load_met_function,varid=varid[3],infile_varid=infile_varid[3],
                                 spatial_type=spatial_type,cardamom_ext=cardamom_ext,path_to_met_source=path_to_met_source,
                                 met_source=met_source,wheat=wheat_from_chaff,xy_bounds = xy_bounds,mc.cores = cl)
      }
      print("...met load 50 %")
      if (varid[4] != "") {
        vpd_out_list=mclapply(load_years,FUN=load_met_function,varid=varid[4],infile_varid=infile_varid[4],
                              spatial_type=spatial_type,cardamom_ext=cardamom_ext,path_to_met_source=path_to_met_source,
                              met_source=met_source,wheat=wheat_from_chaff,xy_bounds = xy_bounds,mc.cores = cl)
      }
      print("...met load 67 %")
      if (varid[5] != "") {
        mint_out_list=mclapply(load_years,FUN=load_met_function,varid=varid[5],infile_varid=infile_varid[5],
                               spatial_type=spatial_type,cardamom_ext=cardamom_ext,path_to_met_source=path_to_met_source,
                               met_source=met_source,wheat=wheat_from_chaff,xy_bounds = xy_bounds,mc.cores = cl)
      }
      print("...met load 83 %")
      if (varid[6] != "") {
        wind_out_list=mclapply(load_years,FUN=load_met_function,varid=varid[6],infile_varid=infile_varid[6],
                               spatial_type=spatial_type,cardamom_ext=cardamom_ext,path_to_met_source=path_to_met_source,
                               met_source=met_source,wheat=wheat_from_chaff,xy_bounds = xy_bounds,mc.cores = cl)
      }
      print("...met load 100 %")
      
    } else {
      
      # or use serial
      if (varid[1] != "") {
        swrad_out_list=lapply(load_years,FUN=load_met_function,varid=varid[1],infile_varid=infile_varid[1],
                              spatial_type=spatial_type,cardamom_ext=cardamom_ext,path_to_met_source=path_to_met_source,
                              met_source=met_source,wheat=wheat_from_chaff,xy_bounds = xy_bounds)
      }
      print("...met load 17 %")
      if (varid[2] != "") {
        maxt_out_list=lapply(load_years,FUN=load_met_function,varid=varid[2],infile_varid=infile_varid[2],
                             spatial_type=spatial_type,cardamom_ext=cardamom_ext,path_to_met_source=path_to_met_source,
                             met_source=met_source,wheat=wheat_from_chaff,xy_bounds = xy_bounds)
      }
      print("...met load 33 %")
      if (varid[3] != "") {
        precip_out_list=lapply(load_years,FUN=load_met_function,varid=varid[3],infile_varid=infile_varid[3],
                               spatial_type=spatial_type,cardamom_ext=cardamom_ext,path_to_met_source=path_to_met_source,
                               met_source=met_source,wheat=wheat_from_chaff,xy_bounds = xy_bounds)
      }
      print("...met load 50 %")
      if (varid[4] != "") {
        vpd_out_list=lapply(load_years,FUN=load_met_function,varid=varid[4],infile_varid=infile_varid[4],
                            spatial_type=spatial_type,cardamom_ext=cardamom_ext,path_to_met_source=path_to_met_source,
                            met_source=met_source,wheat=wheat_from_chaff,xy_bounds = xy_bounds)
      }
      print("...met load 67 %")
      if (varid[5] != "") {
        mint_out_list=lapply(load_years,FUN=load_met_function,varid=varid[5],infile_varid=infile_varid[5],
                             spatial_type=spatial_type,cardamom_ext=cardamom_ext,path_to_met_source=path_to_met_source,
                             met_source=met_source,wheat=wheat_from_chaff,xy_bounds = xy_bounds)
      }
      print("...met load 83 %")
      if (varid[6] != "") {
        wind_out_list=lapply(load_years,FUN=load_met_function,varid=varid[6],infile_varid=infile_varid[6],
                             spatial_type=spatial_type,cardamom_ext=cardamom_ext,path_to_met_source=path_to_met_source,
                             met_source=met_source,wheat=wheat_from_chaff,xy_bounds = xy_bounds)
      }
      print("...met load 100 %")
      
    } # parallel option
    
    # user update
    print("...beginning restructuring of meteorological datasets")
    # Declare variables with an initial value which will be removed later
    swrad_out = 0 ; maxt_out = 0 ; precip_out = 0 ; vpd_out = 0 ; mint_out = 0 ; wind_out = 0 ; tmp_out = 0; t_grid = 0
    # Start extracting from the output list objects the data we actually want
    if (varid[1] != "") {
      for (i in seq(1, length(swrad_out_list))) {
        swrad_out=append(swrad_out,swrad_out_list[[i]]$var_out)
        t_grid=append(t_grid,swrad_out_list[[i]]$t_grid)
      }
      rm(swrad_out_list)
    }
    if (varid[2] != "") {
      for (i in seq(1, length(maxt_out_list))) {maxt_out=append(maxt_out,maxt_out_list[[i]]$var_out)}
      rm(maxt_out_list)
    }
    if (varid[3] != "") {
      for (i in seq(1, length(precip_out_list))) {precip_out = append(precip_out,precip_out_list[[i]]$var_out)}
      rm(precip_out_list)
    }
    if (varid[4] != "") {
      for (i in seq(1, length(vpd_out_list))) {vpd_out = append(vpd_out,vpd_out_list[[i]]$var_out)}
      rm(vpd_out_list)
    }
    if (varid[5] != "") {
      for (i in seq(1, length(mint_out_list))) {mint_out = append(mint_out,mint_out_list[[i]]$var_out)}
      rm(mint_out_list)
    }
    if (varid[6] != "") {
      for (i in seq(1, length(wind_out_list))) {wind_out = append(wind_out,wind_out_list[[i]]$var_out)}
      rm(wind_out_list)
    }
    
    # remove initial value
    swrad_out = swrad_out[-1] ; maxt_out = maxt_out[-1]
    precip_out = precip_out[-1] ; vpd_out = vpd_out[-1]
    mint_out = mint_out[-1] ; wind_out = wind_out[-1]
    # we want total t_grid only
    t_grid = sum(t_grid)
    # Force tidy of memory
    gc()
    
    #
    # Unit conversions based on specific datasets
    #
    
    # convert Trendy air temperature of oC to K
    if (met_source == "trendy") {
      maxt_out = maxt_out + 273.15
      mint_out = mint_out + 273.15
    }
    
    # convert ERA from W/m2 to MJ/m2/day
    if (met_source == "ERA" | met_source == "isimip3a") {
      swrad_out = swrad_out * 86400 * 1e-6
    }
    
    # convert precipitation from kgH2O/m2/day -> per second
    if (met_source == "ERA" | met_source == "isimip3a") {
      precip_out = precip_out / 86400
    }
    
    # generate some additional timing and data
    # NOTE that if we are using an extra of data for the purposes of a GSI model run then the number of days should be the number intended for simulation
    # not to match with the met files just at the moment.
    # This will be corrected later on when the location specific data are extracted
    print("...generating day of year variables")
    
    # Read in Global CO2 source (ppm)
    # Default is the Mauna Loa background. The default dataset is provided with CARDAMOM source code and covers 1959-2021 at monthly time step.
    # The data are sourced from https://www.esrl.noaa.gov/gmd/ccgg/trends/data.html .
    # The interpolated (gap-filled) observations were extracted and processed into a simple file for our use.
    # Assuming the same name and file convension are used, alternate sources can be easily used.
    # Assumed format, csv col 1 = year, 2 = month, 3 co2_ppm
    co2_background = read.csv(paste(path_to_co2,"/co2_monthly.csv",sep=""), header=TRUE)
    if (years_to_load[1] < co2_background$year[1] | years_to_load[length(years_to_load)] > co2_background$year[dim(co2_background)[1]]) {
      stop("Available CO2 information in ./<path_to_co2>/co2_monthly.csv does not cover project time frame...")
    }
    for (yr in seq(1,length(years_to_load))) {
      # if this includes the extra year add it on to the beginning
      if (extra_year) {extra_nos_days = nos_days_in_year(load_years[1])}
      # is current year a leap or not
      nos_days = nos_days_in_year(years_to_load[yr])
      # Extract this years CO2
      co2_annual = co2_background$co2_ppm[which(co2_background$year == years_to_load[yr])]
      if (nos_days == 366) {
        days_per_month=c(31,29,31,30,31,30,31,31,30,31,30,31)
      } else {
        days_per_month=c(31,28,31,30,31,30,31,31,30,31,30,31)
      }
      # reconstruct arrays before doing site level
      if (yr == 1) {
        doy = seq(1,nos_days)
        # mass of dry air = 28.97(g/mol) ; mass of co2 = 44.01 (g/mol); *1e-6 scale from umol -> mol
        co2 = rep(co2_annual,times = days_per_month) # specifically generate the right number of days per month
        co2 = rep(co2, each = max(1,steps_in_day))
        #co2 = rep(380,length.out=(nos_days*steps_in_day)) # default ppm 570 / 370 face
        if (extra_year) {
          if (extra_nos_days == 366) {
            days_per_month=c(31,29,31,30,31,30,31,31,30,31,30,31)
          } else {
            days_per_month=c(31,28,31,30,31,30,31,31,30,31,30,31)
          }
          extra_co2 = rep(co2_annual,times = days_per_month) # specifically generate the right number of days per month
          extra_co2 = rep(extra_co2, each = max(1,steps_in_day))
          co2 = append(extra_co2,co2)
        }
      } else {
        doy = append(doy,seq(1,nos_days))
        extra_co2 = rep(co2_annual,times = days_per_month) # specifically generate the right number of days per month
        extra_co2 = rep(extra_co2, each = max(1,steps_in_day))
        co2 = append(co2, extra_co2)
        #                co2 = append(co2,rep(380,length.out=(nos_days*steps_in_day)))
      }
    } # end of years loop
    
    # create day of run variable
    run_day = seq(1:length(doy))
    
    # user update
    print("...preparing final met output variables")
    # restructure to output variables
    swrad_out  = array(swrad_out,  dim=c(length(wheat_from_chaff),t_grid)) # SW Rad (MJ/m2/day)
    print(dim(swrad_out))
    maxt_out   = array(maxt_out,   dim=c(length(wheat_from_chaff),t_grid)) # Max T (K)
    precip_out = array(precip_out, dim=c(length(wheat_from_chaff),t_grid)) # Precipitation (kgH2O/m2/s)
    vpd_out    = array(vpd_out,    dim=c(length(wheat_from_chaff),t_grid)) # VPD (Pa)
    mint_out   = array(mint_out,   dim=c(length(wheat_from_chaff),t_grid)) # Min T (K)
    wind_out   = array(wind_out,   dim=c(length(wheat_from_chaff),t_grid)) # Wind (m/s)
    
    # output variables
    met_all = list(lat=lat,long=long,run_day=run_day,wheat=wheat_from_chaff,t_grid=t_grid
                   ,maxt=maxt_out,swrad=swrad_out,co2=co2,doy=doy,precip=precip_out,vpd=vpd_out
                   ,mint=mint_out,wind_spd=wind_out,extra_year=extra_year)
    
    # quick sanity check
    if (any(is.na(as.vector(met_all$swrad)))) {stop(paste("SW_RAD contains NaN summary: ",summary(as.vector(met_all$swrad)),sep="")) }
    if (min(as.vector(met_all$swrad)) < -1) {stop(paste("SW_RAD summary: ",summary(as.vector(met_all$swrad)),sep="")) }
    met_all$swrad[which(met_all$swrad < 0)] = 0
    
    # clean up loose memory
    rm(lat_dim,long_dim,years_to_load,load_years,nos_days,
       lat,long,run_day,wheat_from_chaff,t_grid,swrad_out,mint_out,
       co2,doy,precip_out,vpd_out,extra_year,maxt_out,wind_out)#,check1)
    
  } # site_specific or not
  
  # final tidy up and return
  gc(reset=TRUE,verbose=FALSE)
  return(met_all)
  
} # function end
## Use byte compile
load_met_fields_for_extraction<-cmpfun(load_met_fields_for_extraction)