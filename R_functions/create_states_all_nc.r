
###
## Two functions to create and load site CARDAMOM runs with netcdf file formats
###

create_states_all_nc<-function(PROJECT) {

  # Loop through each site in turn
  for (n in seq(1, PROJECT$nosites)) {

       ###
       ## Load the states_all object
       ###

       load(paste(PROJECT$results_processedpath,PROJECT$sites[n],".RData",sep=""))

       ###
       ## Define any bespoke additions
       ###

       ###
       ## Create file and its basic common variables
       ###

       # Create timing variables needed
       years = as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)
       PROJECT$nos_years = length(years)
       steps_per_year = floor(dim(drivers$met)[1] / PROJECT$nos_years)
       time_dim = length(drivers$met[,1])
       # If timestep_days is of length == 1 then we are using a daily time step, but timstep_days is
       # required to be a vector of length equal to the analysis...make it so...
       if (length(PROJECT$model$timestep_days) == 1) {
           PROJECT$model$timestep_days = rep(PROJECT$model$timestep_days,
                                             length.out = steps_per_year*nos_years)
       }
       # Determine the number of days past since the beginning of the analysis
       time = cumsum(PROJECT$model$timestep_days)

       # We have three list objects we want to preserve the information of
       # if these files are read back in again. The states_all and drivers objects
       # are pretty complex so these will be kept seperate, while parameters
       # will be folded into the states_all file. Below is a list of all potential
       # variables.

       # Define the possible dimensions
       time_dimen <- ncdim_def( "time", units="", 1:time_dim)
       ensemble_dimen <- ncdim_def( "ensemble", units="-", 1:length(states_all$lai_m2m2[,1]))
       year_dimen <- ncdim_def( "year", units="", years)
       nchains_dimen <- ncdim_def("nochains", units="", 1:dim(parameters)[2])
       ensemble_per_chain_dimen <- ncdim_def("ensemble_per_chain_dimen", units="", 1:dim(parameters)[3])
       nprior_dimen <- ncdim_def( "prior_dim", units="", 1:length(drivers$parpriors))
       notherprior_dimen <- ncdim_def( "otherprior_dim", units="", 1:length(drivers$otherpriors))
       npar_dimen <- ncdim_def( "nos_parameters", units="", 1:(max(PROJECT$model$nopars))) # NOTE: +1 is to account for the log-likelihood
       nparLL_dimen <- ncdim_def( "nos_parameters_plus_loglikelihood", units="", 1:(max(PROJECT$model$nopars)+1)) # NOTE: +1 is to account for the log-likelihood
       nmet_dimen <- ncdim_def( "nos_met_drivers", units="", 1:dim(drivers$met)[2]) #
       nobs_dimen <- ncdim_def( "nos_obs_drivers", units="", 1:dim(drivers$obs)[2]) #
       scalar_dimen <- ncdim_def("scalar", units="",1) # used for any single value dimension
       dimnchar <- ncdim_def("nchar", "", 1:100, create_dimvar=FALSE ) # Maximum number of characters used in string vector

       # List all valid dimension for later searching
       # NOTE: that this order has been selected to guard against errors such as same lat / long dimension.
       available_dimen = c("time_dimen","ensemble_dimen","year_dimen","nchains_dimen",
                           "ensemble_per_chain_dimen","nprior_dimen","notherprior_dimen",
                           "npar_dimen","nparLL_dimen","scalar_dimen","nmet_dimen",
                           "nobs_dimen")

       ###
       ## Write combined parameter and states_all file
       ###

       ## define output variable
       var0 = ncvar_def("latitude", units = "", longname = "Latitude of the site",
                        dim=list(scalar_dimen),missval = -9999, prec="single", compression = 9)
       var1 = ncvar_def("longitude", units = "", longname = "Longitude of the site",
                        dim=list(scalar_dimen),missval = -9999, prec="single", compression = 9)
       var2 = ncvar_def("Time", units = "d", longname = paste("Monthly time step given in days since 01/01/",PROJECT$start_year,sep=""),
                        dim=list(time_dimen), missval = -9999, prec="single", compression = 9)
       var3 = ncvar_def("steps_per_year", units = "", longname = "Mean number of model time steps per year",
                        dim=list(scalar_dimen), missval = -9999, prec="single", compression = 9)
       var4 = ncvar_def("parameters", units = "1", longname = "Parameters ensembles",
                        dim=list(nparLL_dimen,nchains_dimen,ensemble_per_chain_dimen), missval = -9999, prec="double", compression = 9)

       ## Write the basic file
       # Define the output file name
       # determine what the output file name is here, so that we can check if one already exists
       output_name = paste(PROJECT$results_processedpath,PROJECT$sites[n],"_states_all.nc",sep="")

       # Delete if the file currently exists
       if (file.exists(output_name)) {file.remove(output_name)}
       # Create the empty file space
       new_file=nc_create(filename=output_name, force_v4 = TRUE,
                          vars=list(var0,var1,var2,var3,var4))

       # Load first variable into the file
       # Site latitude
       ncvar_put(new_file, var0, PROJECT$latitude)
       # Site longitude
       ncvar_put(new_file, var1, PROJECT$longitude)
       # Day of analysis
       ncvar_put(new_file, var2, time)
       # Time steps per simulated year
       ncvar_put(new_file, var3, steps_per_year)
       # Posterior parameter ensembles
       ncvar_put(new_file, var4, parameters)

       # Close the existing file to ensure its written to file
       nc_close(new_file)

       ###
       ## Re-open the file so that we can add to it a variable at a time

       new_file <- nc_open( output_name, write=TRUE )

       ###
       ## Determine which variables we will dump into the netcdf file

       # Make list of all variables
       available_variables = names(states_all)
       # Sort remaining into alphabetical order
       available_variables = sort(available_variables)

       ###
       ## ADD VARIABLES

       # Loop through all variable to create a netcdf dump file
       for (v in seq(1, length(available_variables))) {
            # Check that the desired variable exists
            if (exists(x = available_variables[v], where = states_all)) {
                # Check the dimensions of the variable
                tmp_dim = dim(get(available_variables[v], pos = states_all))
                # Check we have dimensions
                if (length(tmp_dim) == 0) {
                    # Then we have only one dimenions which needs to be extracted by length
                    tmp_dim = length(get(available_variables[v], pos = states_all))
                }
                # Create empty list object
                dimension_list = list()
                # Search for matching netcdf dimension for this variable
                for (d in seq(1, length(tmp_dim))) {
                     # Looping all possibilities here
                     found = FALSE ; dd = 0
                     while (found == FALSE) {
                        dd = dd + 1
                        # Check whether the length of the currently proposed dimension
                        # matches that needed for the current dimension
                        if (get(available_dimen[dd])$len == tmp_dim[d]) {
                            found = TRUE
                        }
                        # Ensure we don't get stuck in a loop
                        if (dd > length(available_dimen)) {
                            stop(paste("Could not find valid dimension for ",available_variables[v]," and dimension ",d,sep=""))
                        }
                     } # Looping through available dimensions
                     # We found the correct dimension we can assign it to the list
                     dimension_list[[available_dimen[dd]]] = get(available_dimen[dd])
                } # Loop through the needed dimensions

                # We should now have all the pieces we need to make a dump file of the cardamom output
                var_new  = ncvar_def(available_variables[v], unit="", longname = "", dim=dimension_list, missval = -9999, prec="single",compression = 9)
                new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
                ncvar_put(new_file, var_new,  get(available_variables[v], pos=states_all))

            } # Desired variable exists
       } # Loop through all variables

       # close the file to write to disk
       nc_close(new_file)
       # Update the user
       print("Have created a netcdf file dump of the states_all + parameters object")

       ###
       ## Write drivers file
       ###

       ## define output variable
       var0 = ncvar_def("id", units = "", longname = "Internal CARDAMOM model ID code",
                        dim=list(scalar_dimen), missval = -9999, prec="single", compression = 9)

       ## Write the basic file
       # Define the output file name
       # determine what the output file name is here, so that we can check if one already exists
       output_name = paste(PROJECT$results_processedpath,PROJECT$sites[n],"_drivers.nc",sep="")

       # Delete if the file currently exists
       if (file.exists(output_name)) {file.remove(output_name)}
       # Create the empty file space
       new_file=nc_create(filename=output_name, force_v4 = TRUE,
                          vars=list(var0))

       # Load first variable into the file
       # Internal CARDAMOM model ID code
       ncvar_put(new_file, var0, drivers$id)

       # Close the existing file to ensure its written to file
       nc_close(new_file)

       ###
       ## Re-open the file so that we can add to it a variable at a time

       new_file <- nc_open( output_name, write=TRUE )

       ###
       ## Determine which variables we will dump into the netcdf file

       # Make list of all variables
       available_variables = names(drivers)
       # List the variables we need to remove
       remove_list = c("id")
       # Loop through and remove each in turn
       for (v in seq(1, length(remove_list))) {
            tmp = which(available_variables == remove_list[v])
            if (length(tmp) > 0) {
                available_variables = available_variables[-tmp]
            } else {
                print(paste("The ",remove_list[v]," variable was requested to be removed from available_variables, but could not be found.",sep=""))
            } # Check the desired variable actually exists
       } # Loop each variable to be removed
       # Sort remaining into alphabetical order
       available_variables = sort(available_variables)

       ###
       ## ADD VARIABLES

       # Loop through all variable to create a netcdf dump file
       for (v in seq(1, length(available_variables))) {
            # Check that the desired variable exists
            if (exists(x = available_variables[v], where = drivers)) {
                # Check the dimensions of the variable
                tmp_dim = dim(get(available_variables[v], pos = drivers))
                # Check we have dimensions
                if (length(tmp_dim) == 0) {
                    # Then we have only one dimenions which needs to be extracted by length
                    tmp_dim = length(get(available_variables[v], pos = drivers))
                }
                # Create empty list object
                dimension_list = list()
                # Search for matching netcdf dimension for this variable
                for (d in seq(1, length(tmp_dim))) {
                     # Looping all possibilities here
                     found = FALSE ; dd = 0
                     while (found == FALSE) {
                        dd = dd + 1
                        # Check whether the length of the currently proposed dimension
                        # matches that needed for the current dimension
                        if (get(available_dimen[dd])$len == tmp_dim[d]) {
                            found = TRUE
                        }
                        # Ensure we don't get stuck in a loop
                        if (dd > length(available_dimen)) {
                            stop(paste("Could not find valid dimension for ",available_variables[v]," and dimension ",d,sep=""))
                        }
                     } # Looping through available dimensions
                     # We found the correct dimension we can assign it to the list
                     dimension_list[[available_dimen[dd]]] = get(available_dimen[dd])
                } # Loop through the needed dimensions

                # We should now have all the pieces we need to make a dump file of the cardamom output
                var_new  = ncvar_def(available_variables[v], unit="", longname = "", dim=dimension_list, missval = -9999, prec="double",compression = 9)
                new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
                ncvar_put(new_file, var_new,  get(available_variables[v], pos=drivers))

            } # Desired variable exists
       } # Loop through all variables

       # close the file to write to disk
       nc_close(new_file)
       # Update the user
       print("Have created a netcdf file dump of the drivers object")

  } # loop through sites

  # Return to user
  return(0)

} # end function create_states_all_nc

load_states_all_from_nc<-function(PROJECT,n) {

  ###
  ## Load the states_all and parameters arrays
  ###

  # Create the expected file name
  output_name = paste(PROJECT$results_processedpath,PROJECT$sites[n],"_states_all.nc",sep="")

  # Check the file exists
  if (file.exists(output_name) == FALSE) {
      # Inform the user
      print(paste("The expected file ",output_name," cannot be found",sep=""))
      return(-1)
  }

  # Open the file
  data_in = nc_open(output_name)

  # Specifically extract posterior parameter ensembles
  parameters = ncvar_get(data_in, "parameters")

  # Make a list of all variables
  available_variables = names(data_in$var)
  # Remove parameters from this list, as it has already been extracted
  available_variables = available_variables[-which(available_variables == "parameters")]

  # Create an empty list object
  states_all = list()
  # Begin looping through all variables and adding them to the states_all object
  for (v in seq(1, length(available_variables))) {
       states_all[[available_variables[v]]] = ncvar_get(data_in, available_variables[v])
  } # Looping all variables

  # Close the netcdf file
  nc_close(data_in)

  ###
  ## Load the drivers
  ###

  # Create the expected file name
  output_name = paste(PROJECT$results_processedpath,PROJECT$sites[n],"_drivers.nc",sep="")

  # Check the file exists
  if (file.exists(output_name) == FALSE) {
      # Inform the user
      print(paste("The expected file ",output_name," cannot be found",sep=""))
      return(-1)
  }

  # Open the file
  data_in = nc_open(output_name)

  # Make a list of all variables
  available_variables = names(data_in$var)

  # Create an empty list object
  drivers = list()
  # Begin looping through all variables and adding them to the states_all object
  for (v in seq(1, length(available_variables))) {
       drivers[[available_variables[v]]] = ncvar_get(data_in, available_variables[v])
  } # Looping all variables

  # Close the netcdf file
  nc_close(data_in)

  # Return to user
  return(list(states_all = states_all, parameters = parameters, drivers = drivers))

} # end function load_states_all_from_nc
