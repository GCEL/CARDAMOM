
###
## Two functions to create and load gridded CARDAMOM runs with netcdf file formats
###

create_grid_output_nc<-function(PROJECT) {

  ###
  ## Load the grid_output object
  ###

  load(paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep=""))

  ###
  ## Define any bespoke additions
  ###

  # Specify any extra information for the filename
  output_prefix = "" # follow with "_"
  output_suffix = "" # begin with "_"

  ###
  ## Create file and its basic common variables
  ###

  # Create vectors for latitude / longitude
  longitude = grid_output$long[,1] ; latitude = grid_output$lat[1,]
  # Create timing variables needed
  years = as.numeric(grid_output$start_year):as.numeric(grid_output$end_year)
  nos_years = length(years)
  # If timestep_days is of length == 1 then we are using a daily time step, but timstep_days is
  # required to be a vector of length equal to the analysis...make it so...
  if (length(PROJECT$model$timestep_days) == 1) {PROJECT$model$timestep_days = rep(PROJECT$model$timestep_days, length.out = grid_output$steps_per_year*nos_years)}
  # Determine the number of days past since the beginning of the analysis
  time = cumsum(PROJECT$model$timestep_days)

  # Define the possible dimensions
  lat_dimen <- ncdim_def( "lat_vector", units="degree north (-90->90)", latitude )
  long_dimen <- ncdim_def( "lon_vector", units="degree east (-180->180)", longitude )
  time_dimen <- ncdim_def( "time", units="", 1:grid_output$time_dim)
  nsites_dimen <- ncdim_def( "nsites", units="", 1:PROJECT$nosites)
  quantile_dimen <- ncdim_def( "quantile", units="-", grid_output$num_quantiles)
  year_dimen <- ncdim_def( "year", units="", years)
  npar_dimen <- ncdim_def( "nos_parameters", units="", 1:(max(PROJECT$model$nopars))) # NOTE: +1 is to account for the log-likelihood
  nparLL_dimen <- ncdim_def( "nos_parameters_plus_loglikelihood", units="", 1:(max(PROJECT$model$nopars)+1)) # NOTE: +1 is to account for the log-likelihood
  scalar_dimen <- ncdim_def("scalar", units="",1) # used for any single value dimension
  readme_dimen <- ncdim_def("readme_length", units="",length(grid_output$readme)) #
  dimnchar <- ncdim_def("nchar", "", 1:100, create_dimvar=FALSE ) # Maximum number of characters used in string vector
  # List all valid dimension for later searching
  # NOTE: that this order has been selected to guard against errors such as same lat / long dimension.
  available_dimen = c("long_dimen","lat_dimen","nsites_dimen","time_dimen",
                      "quantile_dimen","year_dimen","npar_dimen","nparLL_dimen",
                      "scalar_dimen")

  ## define output variable
  var0  = ncvar_def("readme", units = "-", longname = "",
                    dim=list(dimnchar,readme_dimen), prec="char", compression = 9)
  var1  = ncvar_def("Time", units = "d", longname = paste("Monthly time step given in days since 01/01/",PROJECT$start_year,sep=""),
                    dim=list(time_dimen), missval = -9999, prec="single", compression = 9)
  var2  = ncvar_def("start_year", units = "", longname = "First year of the analysis",
                    dim=list(dimnchar,scalar_dimen), prec="char", compression = 9)
  var3  = ncvar_def("end_year", units = "", longname = "Final year of the analysis",
                    dim=list(dimnchar,scalar_dimen), prec="char", compression = 9)
  var4  = ncvar_def("grid_area", units = "m2", longname = "Pixel area",
                    dim=list(long_dimen,lat_dimen), missval = -9999, prec="single", compression = 9)
  var5  = ncvar_def("land_fraction", units = "1", longname = "Fraction of pixel which is land",
                    dim=list(long_dimen,lat_dimen), missval = -9999, prec="single", compression = 9)
  var6  = ncvar_def("landmask", units = "1", longname = "Binary mask of land water classification for each pixel",
                    dim=list(long_dimen,lat_dimen), missval = -9999, prec="single", compression = 9)
  var7  = ncvar_def("i_location", units = "1", longname = "Converts site number into the corresponding i or x coordinate for gridded variables",
                    dim=list(nsites_dimen), missval = -9999, prec="single", compression = 9)
  var8  = ncvar_def("j_location", units = "1", longname = "Converts site number into the corresponding j or y coordinate for gridded variables",
                    dim=list(nsites_dimen), missval = -9999, prec="single", compression = 9)
  var9  = ncvar_def("latitude", units = "degree", longname = "Latitude grid (-90/90)",
                    dim=list(long_dimen,lat_dimen), missval = -9999, prec="single", compression = 9)
  var10 = ncvar_def("longitude", units = "degree", longname = "Longigude grid (-180/180)",
                    dim=list(long_dimen,lat_dimen), missval = -9999, prec="single", compression = 9)
  var11 = ncvar_def("steps_per_year", units = "", longname = "Mean number of model time steps per year",
                    dim=list(scalar_dimen), missval = -9999, prec="single", compression = 9)
  var12 = ncvar_def("parameters", units = "1", longname = "Parameters estimates across space and quantiles",
                    dim=list(long_dimen,lat_dimen,nparLL_dimen,quantile_dimen), missval = -9999, prec="double", compression = 9)
  var13 = ncvar_def("parameters_converged", units = "1", longname = "Parameters pass / fail for Gelman-Rubin convergence criterion",
                    dim=list(long_dimen,lat_dimen,nparLL_dimen), missval = -9999, prec="single", compression = 9)

  ## Write the basic file
  # Define the output file name
  # determine what the output file name is here, so that we can check if one already exists
  output_name = paste(PROJECT$results_processedpath,output_prefix,PROJECT$name,"_stock_flux",output_suffix,".nc",sep="")

  # Delete if the file currently exists
  if (file.exists(output_name)) {file.remove(output_name)}
  # Create the empty file space
  new_file=nc_create(filename=output_name, force_v4 = TRUE,
                     vars=list(var0,var1,var2,var3,var4,var5,var6,var7,var8,
                               var9,var10,var11,var12,var13))

  # Load first variable into the file
  # Brief readme
  ncvar_put(new_file, var0, grid_output$readme)
  # Day of analysis
  ncvar_put(new_file, var1, time)
  # First year of the analysis
  ncvar_put(new_file, var2, grid_output$start_year)
  # Final year of the analysis
  ncvar_put(new_file, var3, grid_output$end_year)
  # Grid area
  ncvar_put(new_file, var4, grid_output$area_m2)
  # Land fraction
  ncvar_put(new_file, var5, grid_output$land_fraction)
  # Binary land mask
  ncvar_put(new_file, var6, grid_output$landmask)
  # i_location
  ncvar_put(new_file, var7, grid_output$i_location)
  # j_location
  ncvar_put(new_file, var8, grid_output$j_location)
  # Latitude as grid (-90/90)
  ncvar_put(new_file, var9, grid_output$lat)
  # Longitude as grid (-180/180)
  ncvar_put(new_file, var10, grid_output$long)
  # Number of time steps per year
  ncvar_put(new_file, var11, grid_output$steps_per_year)
  # Array defining the parameter quantiles estimates across the grid
  ncvar_put(new_file, var12, grid_output$parameters)
  # Array defining which parameters and where passed the
  # Gelman-Rubin convergence criterion.
  ncvar_put(new_file, var13, grid_output$parameters_converged)

  # Close the existing file to ensure its written to file
  nc_close(new_file)

  ###
  ## Re-open the file so that we can add to it a variable at a time
  ###

  new_file <- nc_open( output_name, write=TRUE )

  ###
  ## Determine which variables we will dump into the netcdf file
  ###

  # Make list of all variables
  available_variables = names(grid_output)
  # List the variables we need to remove
  remove_list = c("area_m2","land_fraction","landmask","i_location","j_location",
                  "steps_per_year","parameters","parameters_converged","readme",
                  "start_year","end_year","lat","long","lat_dim","long_dim",
                  "time_dim","nos_years")
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
  ###

  # Loop through all variable to create a netcdf dump file
  for (v in seq(1, length(available_variables))) {
       # Check that the desired variable exists
       if (exists(x = available_variables[v], where = grid_output)) {
           # Check the dimensions of the variable
           tmp_dim = dim(get(available_variables[v], pos = grid_output))
           # Check we have dimensions
           if (length(tmp_dim) == 0) {
               # Then we have only one dimenions which needs to be extracted by length
               tmp_dim = length(get(available_variables[v], pos = grid_output))
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
           ncvar_put(new_file, var_new,  get(available_variables[v], pos=grid_output))

       } # Desired variable exists
  } # Loop through all variables

  # close the file to write to disk
  nc_close(new_file)

  # Update the user
  print("Have created a netcdf file dump of the grid_output object")

  # Return to user
  return(0)

} # end function create_grid_output_nc

load_grid_output_from_nc<-function(PROJECT) {

  # Specify any extra information for the filename
  output_prefix = "" # follow with "_"
  output_suffix = "" # begin with "_"

  # Create the expected file name
  output_name = paste(PROJECT$results_processedpath,output_prefix,PROJECT$name,"_stock_flux",output_suffix,".nc",sep="")
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
  grid_output = list()

  # Specifically extract / determine the following variables
  grid_output$time_dim = length(ncvar_get(data_in, "time"))
  grid_output$years = ncvar_get(data_in,"year")
  grid_output$nos_years = length(grid_output$years)
  grid_output$lat_dim = length(ncvar_get(data_in,"lat_vector"))
  grid_output$long_dim = length(ncvar_get(data_in,"lon_vector"))
  grid_output$num_quantiles = ncvar_get(data_in,"quantile")

  # Begin looping through all variables and adding them to the grid_output object
  for (v in seq(1, length(available_variables))) {
       grid_output[[available_variables[v]]] = ncvar_get(data_in, available_variables[v])
  } # Looping all variables

  # Close the netcdf file
  nc_close(data_in)

  # Return to user
  return(grid_output)

} # end function load_grid_output_from_nc
