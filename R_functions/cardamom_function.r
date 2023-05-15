
###
## CARDAMOM function
## from here all other components are called
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

cardamom <-function (projname,model,method,stage) {
#stage = 2 ; repair = 1 ; use_parallel = TRUE
  ## load needed functions into R environment
  paths = load_paths()

  # Set defaults incase missing, NOTE: <<- to assign global
  if (exists("select_country") == FALSE) {select_country <<- FALSE}
  if (exists("path_to_landsea") == FALSE) {path_to_landsea <<- "default"}
  if (exists("request_compile_local") == FALSE) {request_compile_local <<- FALSE}
  if (exists("path_to_co2") == FALSE) {path_to_co2 <<- "./R_functions/"}
  if (exists("request_extended_mcmc") == FALSE) {request_extended_mcmc <<- FALSE}

  # Use this function to ensure that if the short model name has been provided that we translate
  # this into the full internal code version
  tmp = cardamom_model_details(model,"global",1)
  model = tmp$name

  # define file name for PROJECT file
  # this file will contain all information relating the the PROJECT
  PROJECTfile = paste(paths$cardamom_outputs,model,"_",method,"/",projname,"/infofile.RData",sep="")
  PROJECTtype = paste(model,"_",method,sep="")
  # information to the user
  print(paste("When this all began ",Sys.time(),sep=""))

  if (stage == -1 & model == "ACM") {
      # check the user understand what is about to happen
      failed=TRUE
      while (failed){
        understands=readline("Does the user understand that using the ACM model all datastreams except GPP MUST be set to '  ' and that 'path_to_site_obs' file has all needed information (yes/no)?")
        if (understands != "yes") {failed = TRUE} else {failed = FALSE}
        if (understands == "no") {stop('then you need to read the code to figure out what it assumes....')}
      } # while loop
  } # model ACM

  ###
  ## Begin Stage -1
  ## create or repair PROJECT info file, PROJECT initialiation

  if (file.exists(PROJECTfile) == FALSE | stage == -1){

      # load data for passing to list
      PROJECT=list(name=projname,type=PROJECTtype,source=language,paths=paths)

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

      # write out PROJECT info
#      cardamom_write_project_info(PROJECT)
#      print("Project info written to text document")

#      bob=readline("Freeze / backup project (y/n)")
#      if (bob != "y" | bob != "n") {bob=readline("Freeze / backup project (y/n)")}
      # create a tar.gz of the project for later use
#      if (bob == "y") {cardamom_freeze_code(PROJECT)}

      # save spatial type
      if (cardamom_type == "grid") {
          PROJECT$spatial_type = "grid"
          PROJECT$resolution = cardamom_resolution
          PROJECT$grid_type = cardamom_grid_type
      } else if (cardamom_type == "site") {
          PROJECT$spatial_type = "site"
          PROJECT$resolution = " "
          PROJECT$grid_type = " "
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
      if (cardamom_type == "site") {print(PROJECT)}
      print("Project file saves as R object")
      print(paste("file path = ",PROJECTfile, sep=""))
      # return state of the operation if successful
      return(paste("CARDAMOM Report: ",stage," completed", sep=""))

  } else {

      # Inform the user
      print(paste("Loading the infofile = ",PROJECTfile,sep=""))

      # load PROJECT file from binary file
      load(PROJECTfile)
      if (cardamom_type == "site") {print(PROJECT)}

  } # if (stage == -1)

  ###
  ## Stage 0 is a compile only option. This assumes that the PROJECT has already been created

  if (stage == 0) {

      # call compile of the model
      # define PROJECT set up
      PROJECT = cardamom_project_setup(paths,PROJECT)

  } # if (stage == 0)

  ###
  ## Begin Stage 1

  if (stage == 1) {

      # Check for existing binary input files
      print("Check for existance of binary input files")

      # Check whether there are any files which still need creating
      if (repair == 0) {
          n = 0 ; check = TRUE
          while (n < PROJECT$nosites & check) {
               # Increment
               n = n + 1
               # Check whether the current file exists
               if (exists(paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[n],".bin",sep="")) == FALSE) {
                   # Escape loop now we know some sites need making
                   check = FALSE
               }
          }
          # If all files are present we can stop the current proess
          if (check) {return(paste("CARDAMOM Report: ",stage," completed", sep=""))}
      } # repair == 0

      # Update the user
      print("Beginning creation of binary input files")

      # flag for met drivers load
      met_all = 0 ; lai_all = 0 ; Csom_all = 0 ; forest_all = 0 ; Cwood_all = 0
      # load from PROJECT time step information
      timestep_days = PROJECT$model$timestep_days
      noyears = length(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year))
      # Determine location information
      if (cardamom_type == "grid") {
          print("Determining number / locations of grid points for this run ...")
          output = determine_lat_long_needed(PROJECT$latitude,PROJECT$longitude,PROJECT$resolution,PROJECT$grid_type,PROJECT$waterpixels)
          print("Have now determined grid point locations")
          # Bind together the latitude / longitudes for the grid, extract grid information needed to aid further processing (cardamom_ext)
          latlon = cbind(output$lat,output$long) ; cardamom_ext = output$cardamom_ext
          # Extract the lat / long grids designed for determining the extraction locations for the gridded observation datasets
          obs_long_grid = output$obs_long_grid ; obs_lat_grid = output$obs_lat_grid
          # Tidy up
          rm(output) ; gc(reset=TRUE,verbose=FALSE)
      } else if (cardamom_type != "grid") {
          print("Determining number / locations of grid points for this run ...")
          # Combine the latitude / longitude from the site list
          latlon = cbind(PROJECT$latitude,PROJECT$longitude)
          # However we still need a reduced area domain for extracting the site level analyses. +c(-0.5,0.5) allows buffer
          output = determine_lat_long_needed(lat = range(PROJECT$latitude)+c(-0.5,0.5), long = range(PROJECT$longitude)+c(-0.5,0.5)
                                            ,resolution = 0.125, grid_type = "wgs84", remove = 0)
          cardamom_ext = output$cardamom_ext ; gc(reset=TRUE,verbose=FALSE)
          # Extract the lat / long grids designed for determining the extraction locations for the gridded observation datasets
          obs_long_grid = output$obs_long_grid ; obs_lat_grid = output$obs_lat_grid
          # Tidy up
          rm(output) ; gc(reset=TRUE,verbose=FALSE)
          print("Have now determined grid point locations")
      }

      # Load available data from gridded datasets?
      if (PROJECT$model$name != "ACM") {
          # if this is the first time of all creating new met files this time round load the whole dataset for rapid access
          met_all = load_met_fields_for_extraction(latlon,met_source,PROJECT$model$name,PROJECT$start_year,PROJECT$end_year,PROJECT$spatial_type,cardamom_ext)
          lai_all = load_lai_fields_for_extraction(latlon,lai_source,as.character(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)),cardamom_ext,PROJECT$spatial_type)
          nbe_all = load_nbe_fields_for_extraction(latlon,nbe_source,as.character(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)),cardamom_ext,PROJECT$spatial_type)
          gpp_all = load_gpp_fields_for_extraction(latlon,GPP_source,as.numeric(PROJECT$start_year),as.numeric(PROJECT$end_year),cardamom_ext,PROJECT$spatial_type)
          fire_all = load_fire_emission_fields_for_extraction(latlon,fire_source,as.numeric(PROJECT$start_year),as.numeric(PROJECT$end_year),cardamom_ext,PROJECT$spatial_type)
          Csom_all = load_Csom_fields_for_extraction(latlon,Csom_source,cardamom_ext,PROJECT$spatial_type)
          crop_man_all = load_sacks_calendar_fields_for_extraction(latlon,crop_management_source)
          sand_clay_all = load_sand_clay_fields_for_extraction(latlon,sand_clay_source,cardamom_ext,PROJECT$spatial_type)
          forest_all = load_forestry_fields_for_extraction(latlon,deforestation_source,as.character(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)),cardamom_ext,PROJECT$spatial_type)
          Cwood_initial_all = load_initial_biomass_maps_for_extraction(latlon,Cwood_initial_source,as.numeric(PROJECT$start_year),as.numeric(PROJECT$end_year),timestep_days,cardamom_ext,PROJECT$spatial_type)
          Cwood_stock_all = load_biomass_stocks_maps_for_extraction(latlon,Cwood_stock_source,as.numeric(PROJECT$start_year),as.numeric(PROJECT$end_year),timestep_days,cardamom_ext,PROJECT$spatial_type)
          Cwood_potential_all = load_potential_biomass_maps_for_extraction(latlon,Cwood_potential_source,as.numeric(PROJECT$start_year),as.numeric(PROJECT$end_year),timestep_days,cardamom_ext,PROJECT$spatial_type)
          burnt_all = load_burnt_area_fields_for_extraction(latlon,burnt_area_source,path_to_burnt_area,as.numeric(PROJECT$start_year),as.numeric(PROJECT$end_year),cardamom_ext,PROJECT$spatial_type)
          soilwater_all = load_soilwater_fields_for_extraction(latlon,soilwater_initial_source)
          lca_all = load_lca_maps_for_extraction(latlon,lca_source,cardamom_ext,PROJECT$spatial_type)
          Cwood_inc_all = load_wood_productivity_maps_for_extraction(Cwood_inc_source,cardamom_ext,PROJECT$spatial_type,latlon,as.numeric(PROJECT$start_year),as.numeric(PROJECT$end_year),timestep_days)
          Cwood_mortality_all = load_wood_mortality_maps_for_extraction(Cwood_mortality_source,cardamom_ext,PROJECT$spatial_type,latlon,as.numeric(PROJECT$start_year),as.numeric(PROJECT$end_year),timestep_days)
		  lifespan_all = load_lifespan_maps_for_extraction(latlon,lifespan_source,cardamom_ext,PROJECT$spatial_type)
		  leaf_fall_period_all = load_leaf_fall_period_maps_for_extraction(latlon,leaf_fall_period_source,cardamom_ext,PROJECT$spatial_type)
      } # # if (PROJECT$model$name != "ACM")

      # Update user
      print("Loading completed, beginning writing out file write out")

      if (use_parallel) {
          # Create function needed to process the site specific creation
          write_bin_files<-function(n) {

              # create the file name for the met/obs binary
              filename = paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[n],".bin",sep="")

              # All CARDAMOM read gridded datasets now map onto the same projection, extent and resolution.
              # This means that we can extract the location of the current site within any grid just the
              # once and pass around the solution to all the extraction functions.
              # NOTE: that met_all is an exception to this as the for the analysis to work we must always have
              # meteorology driving the model, so we assume a nearest neighbour approach of valid
              # locations, rather than accepting data gaps as done in observations / disturbance drivers.
              output = closest2d_2(1,obs_lat_grid,obs_long_grid,latlon[n,1],latlon[n,2])
              grid_long_loc = unlist(output, use.names=FALSE)[1] ; grid_lat_loc = unlist(output, use.names=FALSE)[2]
              rm(output)

              # load met drivers and obs from CTESSEL and convert to daily if needed
              # these data are site specific so
              if (file.exists(filename) == FALSE | repair == 1){
                  # Load met drivers for ACM or other models
                  if (PROJECT$model$name != "ACM") {
                      met = extract_met_drivers(n,timestep_days,PROJECT$start_year,PROJECT$end_year,latlon[n,],met_all,met_source,PROJECT$sites[n])
                  } else { # if (PROJECT$model$name != "ACM")
                      # assume ACM special case
                      met = extract_acm_met_drivers(PROJECT,latlon[n,],PROJECT$sites[n])
                  }
                  # Load observations
                  obs = extract_obs(grid_long_loc,grid_lat_loc,latlon[n,],lai_all,Csom_all,forest_all
                                   ,Cwood_initial_all,Cwood_stock_all,Cwood_potential_all
                                   ,sand_clay_all,crop_man_all,burnt_all,soilwater_all
                                   ,nbe_all, lca_all, lifespan_all, leaf_fall_period_all,gpp_all,Cwood_inc_all,Cwood_mortality_all, fire_all
                                   ,PROJECT$ctessel_pft[n],PROJECT$sites[n],PROJECT$start_year,PROJECT$end_year
                                   ,timestep_days,PROJECT$spatial_type,PROJECT$resolution,PROJECT$grid_type,PROJECT$model$name)

                  # update ctessel pft in the project and potentially the model information
                  PROJECT$ctessel_pft[n] = obs$ctessel_pft
                  # Load additional model information
                  PROJECT$model = cardamom_model_details(PROJECT$model$name,pft_specific_parameters,PROJECT$ctessel_pft)
                  # write out the relevant binary files
                  binary_data(met,obs,filename,PROJECT$edc,latlon[n,],PROJECT$ctessel_pft[n],
                              PROJECT$model$name,PROJECT$parameter_type,PROJECT$model$nopars[n],noyears)
              } # if (file.exists(filename) == FALSE | repair == 1)
          } # end function

          # NOTE: that the use of mclapply() is due to reported improved efficiency over creating a virtual cluster.
          # However, mclapply does not (at the time of typing) work on Windows, i.e. Linux and Mac only
          cl <- min(PROJECT$nosites,numWorkers)
          dummy = mclapply(c(1:PROJECT$nosites), FUN = write_bin_files, mc.cores = cl)
      } else { # use parallel
          # start looping through sites to create site specific files of obs and met
          for (n in seq(1, PROJECT$nosites)) {

               print(paste("Site ",n," of ",PROJECT$nosites," ",Sys.time(),sep=""))
               # create the file name for the met/obs binary
               filename=paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[n],".bin",sep="")

               # All CARDAMOM read gridded datasets now map onto the same projection, extent and resolution.
               # This means that we can extract the location of the current site within any grid just the
               # once and pass around the solution to all the extraction functions.
               # NOTE: that met_all is an exception to this as the for the analysis to work we must always have
               # meteorology driving the model, so we assume a nearest neighbour approach of valid
               # locations, rather than accepting data gaps as done in observations / disturbance drivers.
               output = closest2d_2(1,obs_lat_grid,obs_long_grid,latlon[n,1],latlon[n,2])
               grid_long_loc = unlist(output, use.names=FALSE)[1] ; grid_lat_loc = unlist(output, use.names=FALSE)[2]
               rm(output)

               # load met drivers and obs from CTESSEL and convert to daily if needed
               # these data are site specific so
               if (file.exists(filename) == FALSE | repair == 1){
                   # Load met drivers for ACM or other models
                   if (PROJECT$model$name != "ACM") {
                       met = extract_met_drivers(n,timestep_days,PROJECT$start_year,PROJECT$end_year,latlon[n,],met_all,met_source,PROJECT$sites[n])
                   } else { # if (PROJECT$model$name != "ACM")
                       # assume ACM special case
                       met = extract_acm_met_drivers(PROJECT,latlon[n,],PROJECT$sites[n])
                   }
                   # Load observations
                   obs = extract_obs(grid_long_loc,grid_lat_loc,latlon[n,],lai_all,Csom_all,forest_all
                                    ,Cwood_initial_all,Cwood_stock_all,Cwood_potential_all
                                    ,sand_clay_all,crop_man_all,burnt_all,soilwater_all
                                    ,nbe_all, lca_all, lifespan_all,leaf_fall_period_all,gpp_all,Cwood_inc_all,Cwood_mortality_all, fire_all
                                    ,PROJECT$ctessel_pft[n],PROJECT$sites[n],PROJECT$start_year,PROJECT$end_year
                                    ,timestep_days,PROJECT$spatial_type,PROJECT$resolution,PROJECT$grid_type,PROJECT$model$name)

                   # update ctessel pft in the project and potentially the model information
                   PROJECT$ctessel_pft[n] = obs$ctessel_pft
                   # Load additional model information
                   PROJECT$model = cardamom_model_details(PROJECT$model$name,pft_specific_parameters,PROJECT$ctessel_pft)
                   # write out the relevant binary files
                   binary_data(met,obs,filename,PROJECT$edc,latlon[n,],PROJECT$ctessel_pft[n],
                               PROJECT$model$name,PROJECT$parameter_type,PROJECT$model$nopars[n],noyears)
               } # if (file.exists(filename) == FALSE | repair == 1)
          } # site loop
      } # use_parallel

      # clean up to remove large drains on computer memeory
      rm(met_all,lai_all) ; gc(reset=TRUE,verbose=FALSE)

      # copy files to eddie?
      if (PROJECT$ecdf) {
          failed = TRUE ; copy_to = "y" ; failed = FALSE
          while(failed) {
             copy_to = readline("Copy binary driver files to cluster? (y/n)")
              if (copy_to != "y" & copy_to != "n") {failed = TRUE} else {failed = FALSE}
          }
          if (copy_to == "y") {
             #home_computer=Sys.info()["nodename"]
              # Check whether an existing cardamom_inputs.zip exists
              if (file.exists(paste(PROJECT$datapath,"cardamom_inputs.zip",sep=""))) {
                  # Delete this file before creating a new one of the latest input files
                  system(paste("rm ",PROJECT$datapath,"cardamom_inputs.zip", sep=""))
              }
              # Compress all input files into zip directory
              system(paste("zip -j -r -q ",PROJECT$datapath,"cardamom_inputs.zip ",PROJECT$datapath," -i '*.bin'",sep=""))
              # Copy the zip directory to the remote server
              command = paste("scp -r -q ",username,"@",home_computer,":",PROJECT$datapath,"cardamom_inputs.zip ",PROJECT$edatapath,sep="")
              # Unzip on remote server
              command = c(command,paste("unzip -o -qq ",PROJECT$edatapath,"cardamom_inputs.zip -d ",PROJECT$edatapath, sep=""))
              # Remove the zip directory on remote server
              command = c(command,paste("rm ",PROJECT$edatapath,"cardamom_inputs.zip" ,sep=""))
              #command = paste("scp -r ",username,"@",home_computer,":",PROJECT$datapath,"* ",PROJECT$edatapath,sep="")
              print(command)
              # Execute command on remote server
              ecdf_execute(command,PROJECT$paths$cardamom_cluster)
              # Delete local copy of the zip directory
              system(paste("rm ",PROJECT$datapath,"cardamom_inputs.zip", sep=""))
          }
      } # copy to Eddie

      # report to the user
      return(paste("CARDAMOM Report: ",stage," completed", sep=""))
  } # if stage == 1

  ###
  ## Begin Stage 2

  if (stage == 2) {

      print('Welcome to Stage 2 - running CARDAMOM')
      print('The code will be run on cluster or your local machine');

      if (PROJECT$ecdf) {
          # submit files to eddie
          submit_processes_to_cluster(PROJECT)
      } else {
          # submit to local machine
          submit_processes_to_local_machine(PROJECT)
      }
      # report to the user
      return(paste("CARDAMOM Report: ",stage," completed", sep=""))

  } # stage == 2

  ###
  ## Begin Stage 3

  if (stage == 3) {

      print("Stage 3 will copy files back from cluster and begin postprocessing")
      print("NOTE: this will only be effective if cluster has completed its tasks")

      if (PROJECT$ecdf) {
          failed = TRUE
          while(failed) {
             # do we copy back the files?
             copy_back = readline("Copy results back from cluster? (y/n)")
             if (copy_back != "y" & copy_back != "n") { failed=TRUE } else { failed=FALSE }
          }
          # Are we copying back the files
          if (copy_back == "y") {
              #home_computer=Sys.info()["nodename"]
              # If yes, then we mist delete the existing files to ensure we do not mix analysis versions
              if (length(list.files(paste(PROJECT$resultspath,"/",sep=""))) > 0) {
                  system(paste("find ",PROJECT$resultspath," -type f -name '*' -delete",sep=""))
              }
              # Ensure any existing zip directory has been deleted before creating a new one.
              command = paste("rm ",PROJECT$eresultspath,"cardamom_outputs*.zip",sep="")
              # Prepare copy back command
              # Compress all existing files into zip directory
              # There is a limit on how many files (based on the command length) that can be added at once using zip alone.
              # However, we can get around this by listing all files using find and then piping these into zip
              for (i in seq(1,PROJECT$nochains)) {
                   command = c(command,paste("zip -j -r -q ",PROJECT$eresultspath,"cardamom_outputs_",i,".zip ",PROJECT$eresultspath," -i '*_",i,"_PARS'",sep=""))
              }
              command = c(command,paste("scp -r -q ",PROJECT$eresultspath,"cardamom_outputs*.zip ",username,"@",home_computer,":",PROJECT$resultspath,sep=""))
              #command = c(command,paste("rm ",PROJECT$eresultspath,"cardamom_outputs.zip",sep=""))
              #command = paste("scp -r ",PROJECT$eresultspath,"* ",username,"@",home_computer,":",PROJECT$resultspath,sep="")
              # Execute on remote server
              ecdf_execute(command,PROJECT$paths$cardamom_cluster)
          } # copy back
          # Locally check that we have the cardamom_outputs.zip copied back from remote server.
          # Assuming they are present unzip and delete the zip directory
          for (i in seq(1,PROJECT$nochains)) {
               if (file.exists(paste(PROJECT$resultspath,"cardamom_outputs_",i,".zip",sep=""))) {
                   # Unzip
                   system(paste("unzip -qq -o ",PROJECT$resultspath,"cardamom_outputs_",i,".zip -d ",PROJECT$resultspath, sep=""))
                   # Delete file now
                   system(paste("rm ",PROJECT$resultspath,"cardamom_outputs_",i,".zip", sep=""))
               }
          }
      } # ecdf condition
      # do we run the parameters yet for analysis
      # Changed to a hardcoded run of the analysis
      run_all = "y"#readline("Run all parameter vectors to generate confidence intervals? (y/n)")
      failed = TRUE
      while(failed) {
         if (run_all != "y" & run_all != "n") {run_all = readline("Run all parameter vectors to generate confidence intervals? (y/n)") ; failed=TRUE} else {failed = FALSE}
      }
      # If we are running
      if (run_all == "y") {
          # Specify how much of the storged parameter set to run
          # NOTE: this is deprecated in favour of hardcoded last 100 of each parameter file.
          # This is more resilient to running incomplete chains
          PROJECT$latter_sample_frac = 0.75 #0.5 # 0.75 #readline("What (latter) fraction of accepted parameters to use (e.g. 0.5)?")
          # Run the parameter back through DALEC
          run_mcmc_results(PROJECT,stage,repair,grid_override)
      }

      # now save the project
      save(PROJECT,file=PROJECTfile)

      # report to the user
      return(paste("CARDAMOM Report: ",stage," completed", sep=""))

  } # stage == 3

  ###
  ## Begin Stage 4

  if (stage == 4 | stage == 4.5) {

      print("Beginning stage 4: generating stardard outputs")

      # Generating site level plots or gridded
      if (PROJECT$spatial_type == "site" | grid_override) {

          # Generate figures of parameters and model values including
          # uncertainty information
          generate_uncertainty_figures(PROJECT)

      } else if (PROJECT$spatial_type == "grid") {

          # will generate spatial maps instead
          generate_parameter_maps(PROJECT)
          if (stage == 4.5) {
              # Overly complicated piece of code superceeded by simplifed approach.
              # This function and associated code could be removed safely
              generate_stocks_and_fluxes_maps(PROJECT)
          } else {
              generate_simplified_stock_and_flux_maps(PROJECT)
          }

      } else {
          stop('missing spatial_type definition (i.e. grid or site)')
      } # grid or site run

      # report to the user
      return(paste("CARDAMOM Report: ",stage," completed", sep=""))

  } # stage == 4

  ###
  ## Begin Stage 5

  # Currently empty
  if (stage == 5) {

      print("Stage 5 current uncoded and open for new uses")

      # report to the user
      return(paste("CARDAMOM Report: ",stage," completed", sep=""))

  } # stage == 5

} # end function cardamom

## Use byte compile
cardamom<-cmpfun(cardamom)
