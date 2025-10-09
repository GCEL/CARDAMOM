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
# Function to carry out stage -5 processes,
# i.e. create a CARDAMOM project
# Author: T. Luke Smallman (02/05/2024)
#
#########################################################################################

cardamom_stage_minus_1<-function(PROJECTfile,PROJECTtype,paths,model,method,projname) { 
  
  # load data for passing to list
  PROJECT=list(name=projname,
               type=PROJECTtype,
               source=language,
               paths=paths,
               request_cost_function_scaling = request_cost_function_scaling)
  
  # create sites names file name
  if (cardamom_type == "site") {
    PROJECT$nosites=length(sites_cardamom)
    PROJECT$waterpixels = 0
    PROJECT$landsea = 1 # 1=land, 0 = sea # NOTE not used for site runs
    PROJECT$sites = sites_cardamom
    if (pft_wanted) {
      site_info = find_pft(sites_cardamom_lat,sites_cardamom_long)
      PROJECT$ctessel_pft = site_info$ctessel_pft
    } else {
      PROJECT$ctessel_pft = rep(0, length.out=PROJECT$nosites)
    } # pft_wanted
  } else { # if (cardamom_type == "site")
    site_info = how_many_points(path_to_landsea,sites_cardamom_lat,sites_cardamom_long,
                                cardamom_resolution,cardamom_grid_type,sites_cardamom)
    PROJECT$nosites = site_info$nosites
    PROJECT$sites = site_info$sites
    PROJECT$landsea = site_info$landsea
    PROJECT$waterpixels = site_info$waterpixels
    if (pft_wanted) {
      PROJECT$ctessel_pft = site_info$ctessel_pft
    } else {
      PROJECT$ctessel_pft = rep(0,length.out=length(site_info$ctessel_pft))
    } # pft_wanted
    PROJECT$lat_dim = site_info$lat_dim
    PROJECT$long_dim = site_info$long_dim
    if (exists(x = "landsea_frac", where = site_info)) {
      PROJECT$landsea_frac = site_info$landsea_frac
    }
    if (exists(x = "area_m2", where = site_info)) {
      PROJECT$area_m2 = site_info$area_m2
    }
  } # if (cardamom_type == "site")
  
  # whether we use PFT specific or global parameter ranges
  if (pft_specific_parameters) {
    PROJECT$parameter_type = "pft_specific"
  } else if (pft_specific_parameters == FALSE) {
    PROJECT$parameter_type = "global"
  } else {
    stop(paste("user has not defined logistic 'pft_specific_parameters' option",sep=""))
  }
  
  # Load additional model information
  PROJECT$model = cardamom_model_details(model,PROJECT$parameter_type,PROJECT$ctessel_pft)
  # define PROJECT set up
  PROJECT = cardamom_project_setup(paths,PROJECT)
  
  # use EDCs
  failed=TRUE
  while (failed){
    if (exists("request_use_EDCs")) {
      if (request_use_EDCs == TRUE | request_use_EDCs == FALSE) {
        tmp = request_use_EDCs
      } else {
        tmp = readline("Use EDCs (TRUE/FALSE)?")
      }
    } else {
      tmp = readline("Use EDCs (TRUE/FALSE)?")
    }
    if (tmp == TRUE) {PROJECT$edc = 1} else {PROJECT$edc = 0}
    if (tmp != TRUE & tmp != FALSE) {failed = TRUE} else {failed = FALSE}
  } # whilte(failed)
  # load start and end year
  PROJECT$start_year = years_to_do[1]
  PROJECT$end_year = years_to_do[length(years_to_do)]
  
  # remember LCM used
  PROJECT$land_cover_option=use_lcm
  
  # Create parameter file name
  PROJECT$Rparfile=paste(paths$cardamom_PROJECTs,"/",PROJECT$type,"/","parfile.csv",sep="")
  
  # Save spatial type
  if (cardamom_type == "grid") {
    PROJECT$spatial_type = "grid"
    PROJECT$resolution = cardamom_resolution
    PROJECT$grid_type = cardamom_grid_type
  } else if (cardamom_type == "site") {
    PROJECT$spatial_type = "site"
    PROJECT$resolution = " "
    PROJECT$grid_type = "epsg:4326" # default assumption that the lat / long locations are given in the WGS-84 grid
  } else {
    stop("missing cardamom_type variable")
  } #
  # in both cases we will actually want to have the lat / long information
  PROJECT$latitude = sites_cardamom_lat ; PROJECT$longitude = sites_cardamom_long
  PROJECT$model$timestep = timestep_type
  
  # Determine he number of days per time step, options are monthly, weekly and daily
  if (PROJECT$model$timestep == "monthly") {
    print("...model will use calender monthly timestep ")
    # How many years are we running?
    all_years = as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)
    # How many days in the first year
    nos_days = nos_days_in_year(all_years[1])
    if (nos_days == 366) {timestep_days = c(31,29,31,30,31,30,31,31,30,31,30,31)} else {timestep_days = c(31,28,31,30,31,30,31,31,30,31,30,31)}
    # Determine how many days in each year of the analysis
    for (y in seq(2, length(all_years))) {
      # calculate increment
      nos_days = nos_days_in_year(all_years[y])
      if (nos_days == 366) {
        timestep_days = append(timestep_days,c(31,29,31,30,31,30,31,31,30,31,30,31))
      } else {
        timestep_days = append(timestep_days,c(31,28,31,30,31,30,31,31,30,31,30,31))
      }
    } # loop through days
    # clean up
    rm(all_years)
  } else if (PROJECT$model$timestep == "weekly") {
    print("...model will use weekly timestep ")
    all_years = as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)
    nos_days = nos_days_in_year(all_years[1])
    if (nos_days == 366) {timestep_days = c(rep(7,times=51),9)} else {timestep_days = c(rep(7,times=51),8)}
    for (y in seq(2, length(all_years))) {
      # calculate increment
      nos_days=nos_days_in_year(all_years[y])
      if (nos_days == 366) {
        timestep_days = append(timestep_days,c(rep(7,times=51),9))
      } else {
        timestep_days = append(timestep_days,c(rep(7,times=51),8))
      }
    } # loop through days
    # clean up
    rm(all_years)
  } else if (PROJECT$model$timestep == "daily") {
    print("...model will use daily timestep")
    timestep_days=1
  } else {
    timestep_days=1
    print("WARNING: user has not specified a time step option ('daily', 'weekly' or 'monthly'), default assumption is daily")
  } # time step diagnosis loop
  
  # load to PROJECT for perminent use
  PROJECT$model$timestep_days = timestep_days
  
  # output PROJECT info to R binary file
  save(PROJECT, file=PROJECTfile)
  #if (cardamom_type == "site") {print(PROJECT)}
  print("Project file saves as R object")
  print(paste("file path = ",PROJECTfile, sep=""))
  # return state of the operation if successful
  return(paste("CARDAMOM Report: -1 completed", sep=""))
  
} # end function cardamom_stage_minus_1

## Use byte compile
cardamom_stage_minus_1<-cmpfun(cardamom_stage_minus_1)