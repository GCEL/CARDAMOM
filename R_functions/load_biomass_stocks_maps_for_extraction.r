
###
## Function to load biomass maps which apply to gridded domain
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_biomass_stocks_maps_for_extraction<-function(latlon_in,Cwood_stock_source,start,finish,timestep_days,cardamom_ext,spatial_type) {

    # Generate timing information need in most cases
    analysis_years = seq(as.numeric(start),as.numeric(finish))
    if (length(timestep_days) == 1 & timestep_days[1] == 1) {
        nos_days = 0
        for (t in seq(1,length(analysis_years))) {nos_days = nos_days + nos_days_in_year(analysis_years[t])}
        timestep_days = rep(timestep_days,nos_days)
    }
    # Cumulative of days, used when finding the correct analysis / observtion time step
    run_day_selector = cumsum(timestep_days)
    # How many steps per year
    steps_per_year = length(timestep_days) / length(analysis_years)

    ###
    ## Select the correct Cwood source for specific time points

    if (Cwood_stock_source == "mpi_biomass") {

        # let the user know this might take some time
        print("Loading MPI - >30N Forest Biomass map...")

        # Create the full file paths to both 2010 and 2017 AGB estimates
        input_file = paste(path_to_Cwood,"2014121116258biomass_v3_total.nc",sep="")
        years_with_obs = c(2010)

        # Extract if appropriate
        done_lat = FALSE

        # determine whether the first year is within the analysis period
        if (years_with_obs >= as.numeric(start) & years_with_obs <= as.numeric(finish)) {

            # Open the first file
            data1 = nc_open(input_file)

            # Set flag to TRUE, impacts what will be returned from this function
            done_lat = TRUE
            # Read in lat / long information
            lat = ncvar_get(data1, "latitude") ; long = ncvar_get(data1, "longitude")
            # Create lat / long grid from vectors
            idim = length(long) ; jdim = length(lat)
            lat = array(lat, dim=c(jdim,idim)) ; lat = t(lat)
            long = array(long, dim=c(idim,jdim))

            # Read the biomass estimates and uncertainty
            # NOTE: Units are kgC/m2, conversion to gC/m2
            biomass_gCm2 = ncvar_get(data1, "biomass_total") ; biomass_uncertainty_gCm2 = ncvar_get(data1, "uncertainty_biomass_total")
            # close files after use
            nc_close(data1)

            # Convert to a raster, assuming standad WGS84 grid
            biomass_gCm2 = data.frame(x = as.vector(long), y = as.vector(lat), z = as.vector(biomass_gCm2))
            biomass_gCm2 = rasterFromXYZ(biomass_gCm2, crs = ("+init=epsg:4326"))
            biomass_uncertainty_gCm2 = data.frame(x = as.vector(long), y = as.vector(lat), z = as.vector(biomass_uncertainty_gCm2))
            biomass_uncertainty_gCm2 = rasterFromXYZ(biomass_uncertainty_gCm2, crs = ("+init=epsg:4326"))

            # Create raster with the target crs (technically this bit is not required)
            target = raster(crs = ("+init=epsg:4326"), ext = extent(biomass_gCm2), resolution = res(biomass_gCm2))
            # Check whether the target and actual analyses have the same CRS
            if (compareCRS(biomass_gCm2,target) == FALSE) {
                # Resample to correct grid
                biomass_gCm2 = resample(biomass_gCm2, target, method="ngb") ; gc() ; removeTmpFiles()
                biomass_uncertainty_gCm2 = resample(biomass_uncertainty_gCm2, target, method="ngb") ; gc() ; removeTmpFiles()
            }
            # Extend if required to the target area
            biomass_gCm2 = extend(biomass_gCm2,cardamom_ext) ; biomass_uncertainty_gCm2 = extend(biomass_uncertainty_gCm2,cardamom_ext)
            # Trim the extent of the overall grid to the analysis domain
            biomass_gCm2 = crop(biomass_gCm2,cardamom_ext) ; biomass_uncertainty_gCm2 = crop(biomass_uncertainty_gCm2,cardamom_ext)
            # Remove any missing or un-realistic data points
            biomass_gCm2[which(as.vector(biomass_gCm2) < 0)] = NA
            biomass_uncertainty_gCm2[which(as.vector(biomass_uncertainty_gCm2) < 0)] = NA

            # Adjust spatial resolution of the datasets, this occurs in all cases
            if (res(biomass_gCm2)[1] != res(cardamom_ext)[1] | res(biomass_gCm2)[2] != res(cardamom_ext)[2]) {

                # Create raster with the target resolution
                target = raster(crs = crs(cardamom_ext), ext = extent(cardamom_ext), resolution = res(cardamom_ext))
                # Resample to correct grid
                biomass_gCm2 = resample(biomass_gCm2, target, method="bilinear") ; gc() ; removeTmpFiles()
                biomass_uncertainty_gCm2 = resample(biomass_uncertainty_gCm2, target, method="bilinear") ; gc() ; removeTmpFiles()

            } # Aggrgeate to resolution

            # extract dimension information for the grid, note the axis switching between raster and actual array
            xdim = dim(biomass_gCm2)[2] ; ydim = dim(biomass_gCm2)[1]
            # extract the lat / long information needed
            long = coordinates(biomass_gCm2)[,1] ; lat = coordinates(biomass_gCm2)[,2]
            # restructure into correct orientation
            long = array(long, dim=c(xdim,ydim))
            lat = array(lat, dim=c(xdim,ydim))
            # break out from the rasters into arrays which we can manipulate
            biomass_gCm2 = array(as.vector(unlist(biomass_gCm2)), dim=c(xdim,ydim))
            biomass_uncertainty_gCm2 = array(as.vector(unlist(biomass_uncertainty_gCm2)), dim=c(xdim,ydim))

            # Determine when in the analysis time series the observations should go
            # NOTE: We assume the biomass estimate is placed at the beginning of the year

            # What year of the analysis does the data fall?
            place_obs_in_step = which(run_day_selector >= floor(which(analysis_years == years_with_obs) * 365.25))[1]
            place_obs_in_step = place_obs_in_step - (steps_per_year-1)

        } # Is dataset within the analysis time period?

        # Convert kgCm-2-> gCm-2
        biomass_gCm2 = biomass_gCm2 * 1e3
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 1e3

        # Re-construct arrays for output
        idim = dim(lat)[1] ; jdim = dim(long)[2] ; tdim = length(biomass_gCm2) / (idim * jdim)
        biomass_gCm2 = array(biomass_gCm2, dim=c(idim,jdim,tdim))
        biomass_uncertainty_gCm2 = array(biomass_uncertainty_gCm2, dim=c(idim,jdim,tdim))

        # clean up variables
        gc(reset=TRUE,verbose=FALSE)

        if (done_lat) {
            # Output variables
            return(list(place_obs_in_step = place_obs_in_step, lat = lat, long = long,
                        biomass_gCm2 = biomass_gCm2, biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2))
        } else {
            # Output dummy variables
            return(list(place_obs_in_step = -9999, lat = -9999, long = -9999,
                        biomass_gCm2 = -9999, biomass_uncertainty_gCm2 = -9999))
        } # done_lat

    } else if (Cwood_stock_source == "Avitabile") {

        # let the user know this might take some time
        print("Loading processed Avitabile AGB...")

        # Create the full file paths to both 2010 and 2017 AGB estimates
        input_file = paste(path_to_Cwood,"Biomass_stocks_with_lat_long.nc",sep="")
        years_with_obs = c(2007)

        # Extract if appropriate
        done_lat = FALSE

        # determine whether the first year is within the analysis period
        if (years_with_obs >= as.numeric(start) & years_with_obs <= as.numeric(finish)) {

            # Open the first file
            data1 = nc_open(input_file)

            # Set flag to TRUE, impacts what will be returned from this function
            done_lat = TRUE
            # Read in lat / long information
            lat = ncvar_get(data1, "lat") ; long = ncvar_get(data1, "long")

            # Read the biomass estimates and uncertainty
            # NOTE: Units are MgC/ha, above ground biomass, conversion to gC/m2 and total biomass done later
            biomass_gCm2 = ncvar_get(data1, "Biomass") ; biomass_uncertainty_gCm2 = ncvar_get(data1, "Biomass_Uncertainty")
            # close files after use
            nc_close(data1)

            # Convert to a raster, assuming standad WGS84 grid
            biomass_gCm2 = data.frame(x = as.vector(long), y = as.vector(lat), z = as.vector(biomass_gCm2))
            biomass_gCm2 = rasterFromXYZ(biomass_gCm2, crs = ("+init=epsg:4326"))
            biomass_uncertainty_gCm2 = data.frame(x = as.vector(long), y = as.vector(lat), z = as.vector(biomass_uncertainty_gCm2))
            biomass_uncertainty_gCm2 = rasterFromXYZ(biomass_uncertainty_gCm2, crs = ("+init=epsg:4326"))

            # Create raster with the target crs (technically this bit is not required)
            target = raster(crs = ("+init=epsg:4326"), ext = extent(biomass_gCm2), resolution = res(biomass_gCm2))
            # Check whether the target and actual analyses have the same CRS
            if (compareCRS(biomass_gCm2,target) == FALSE) {
                # Resample to correct grid
                biomass_gCm2 = resample(biomass_gCm2, target, method="ngb") ; gc() ; removeTmpFiles()
                biomass_uncertainty_gCm2 = resample(biomass_uncertainty_gCm2, target, method="ngb") ; gc() ; removeTmpFiles()
            }
            # extend the extent of the overall grid to the analysis domain
            biomass_gCm2 = extend(biomass_gCm2,cardamom_ext) ; biomass_uncertainty_gCm2 = extend(biomass_uncertainty_gCm2,cardamom_ext)
            # Trim the extent of the overall grid to the analysis domain
            biomass_gCm2 = crop(biomass_gCm2,cardamom_ext) ; biomass_uncertainty_gCm2 = crop(biomass_uncertainty_gCm2,cardamom_ext)
            # Remove any missing or un-realistic data points
            biomass_gCm2[which(as.vector(biomass_gCm2) < 0)] = NA
            biomass_uncertainty_gCm2[which(as.vector(biomass_uncertainty_gCm2) < 0)] = NA

            # Adjust spatial resolution of the datasets, this occurs in all cases
            if (res(biomass_gCm2)[1] != res(cardamom_ext)[1] | res(biomass_gCm2)[2] != res(cardamom_ext)[2]) {

                # Create raster with the target resolution
                target = raster(crs = crs(cardamom_ext), ext = extent(cardamom_ext), resolution = res(cardamom_ext))
                # Resample to correct grid
                biomass_gCm2 = resample(biomass_gCm2, target, method="bilinear") ; gc() ; removeTmpFiles()
                biomass_uncertainty_gCm2 = resample(biomass_uncertainty_gCm2, target, method="bilinear") ; gc() ; removeTmpFiles()

            } # Aggrgeate to resolution

            # extract dimension information for the grid, note the axis switching between raster and actual array
            xdim = dim(biomass_gCm2)[2] ; ydim = dim(biomass_gCm2)[1]
            # extract the lat / long information needed
            long = coordinates(biomass_gCm2)[,1] ; lat = coordinates(biomass_gCm2)[,2]
            # restructure into correct orientation
            long = array(long, dim=c(xdim,ydim))
            lat = array(lat, dim=c(xdim,ydim))
            # break out from the rasters into arrays which we can manipulate
            biomass_gCm2 = array(as.vector(unlist(biomass_gCm2)), dim=c(xdim,ydim))
            biomass_uncertainty_gCm2 = array(as.vector(unlist(biomass_uncertainty_gCm2)), dim=c(xdim,ydim))

            # Determine when in the analysis time series the observations should go
            # NOTE: We assume the biomass estimate is placed at the beginning of the year

            # What year of the analysis does the data fall?
            place_obs_in_step = which(run_day_selector >= floor(which(analysis_years == years_with_obs) * 365.25))[1]
            place_obs_in_step = place_obs_in_step - (steps_per_year-1)

        } # Is dataset within the analysis time period?

        # Convert to MgC/ha -> Mg/ha needed for Saatchi et al (2011)
        biomass_gCm2 = biomass_gCm2 * 2.083333
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 2.083333 # ! Based on 250 gC/m2 CI per-com Mat
        # Use allometry to estimate below ground biomass stock and
        # combined with the above ground (Mg/ha) to give a total woody biomass estimate
        # Saatchi et al., (2011), PNAS, 108, 9899-9904, https://www.pnas.org/content/108/24/9899
        biomass_gCm2 = biomass_gCm2 + (0.489 * (biomass_gCm2) ** 0.89)
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 + (0.489 * biomass_uncertainty_gCm2 ** 0.89)
        # Convert units of biomass and its uncertainty from MgCha -> gC/m2
        biomass_gCm2 = biomass_gCm2 * 1e2 * 0.48
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 1e2 * 0.48

        # Re-construct arrays for output
        idim = dim(lat)[1] ; jdim = dim(long)[2] ; tdim = length(biomass_gCm2) / (idim * jdim)
        biomass_gCm2 = array(biomass_gCm2, dim=c(idim,jdim,tdim))
        biomass_uncertainty_gCm2 = array(biomass_uncertainty_gCm2, dim=c(idim,jdim,tdim))

        # clean up variables
        gc(reset=TRUE,verbose=FALSE)

        if (done_lat) {
            # Output variables
            return(list(place_obs_in_step = place_obs_in_step, lat = lat, long = long,
                        biomass_gCm2 = biomass_gCm2, biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2))
        } else {
            # Output dummy variables
            return(list(place_obs_in_step = -9999, lat = -9999, long = -9999,
                        biomass_gCm2 = -9999, biomass_uncertainty_gCm2 = -9999))
        } # done_lat

    } else if (Cwood_stock_source == "INPE_Avitabile") {

        # let the user know this might take some time
        print("Loading processed combined INPE (Amazon) Avitabile (pan-tropics) biomass...")

        # Create the full file paths to both 2007 and 2013 Total Biomass estimates
        input_file = paste(path_to_Cwood,"Cagb_2007_gCm2.tif",sep="")
        input_file = append(input_file,paste(path_to_Cwood,"Cagb_2013_gCm2.tif",sep=""))
        input_file_uncertainty = paste(path_to_Cwood,"Cagb_2007_uncertainty_gCm2.tif",sep="")
        input_file_uncertainty = append(input_file_uncertainty,paste(path_to_Cwood,"Cagb_2013_uncertainty_gCm2.tif",sep=""))
        years_with_obs = c(2007,2013)

        # Loop through each year and extract if appropriate
        done_lat = FALSE
        for (t in seq(1, length(years_with_obs))) {

             # determine whether the first year is within the analysis period
             if (years_with_obs[t] >= as.numeric(start) & years_with_obs[t] <= as.numeric(finish)) {

                 # Read in the estimate and uncertainty rasters
                 biomass = raster(input_file[t])
                 biomass_uncertainty = raster(input_file_uncertainty[t])

                 # Create raster with the target crs
                 target = raster(crs = ("+init=epsg:4326"), ext = extent(biomass), resolution = res(biomass))
                 # Check whether the target and actual analyses have the same CRS
                 if (compareCRS(biomass,target) == FALSE) {
                     # Resample to correct grid
                     biomass = resample(biomass, target, method="ngb") ; gc() ; removeTmpFiles()
                     biomass_uncertainty = resample(biomass_uncertainty, target, method="ngb") ; gc() ; removeTmpFiles()
                 }
                 # Extend the extent of the overall grid to the analysis domain
                 biomass = extend(biomass,cardamom_ext) ; biomass_uncertainty = extend(biomass_uncertainty,cardamom_ext)
                 # Trim the extent of the overall grid to the analysis domain
                 biomass = crop(biomass,cardamom_ext) ; biomass_uncertainty = crop(biomass_uncertainty,cardamom_ext)
                 # now remove the ones that are actual missing data
                 biomass[which(as.vector(biomass) < 0)] = NA
                 biomass_uncertainty[which(as.vector(biomass_uncertainty) < 0)] = NA
                 # Adjust spatial resolution of the datasets, this occurs in all cases
                 if (res(biomass)[1] != res(cardamom_ext)[1] | res(biomass)[2] != res(cardamom_ext)[2]) {

                     # Create raster with the target resolution
                     target = raster(crs = crs(cardamom_ext), ext = extent(cardamom_ext), resolution = res(cardamom_ext))

                     # Resample to correct grid
                     biomass = resample(biomass, target, method="bilinear") ; gc() ; removeTmpFiles()
                     biomass_uncertainty = resample(biomass_uncertainty, target, method="bilinear") ; gc() ; removeTmpFiles()

                 } # Aggrgeate to resolution

                 # If the first file to be read extract the lat / long information
                 if (done_lat == FALSE) {
                     # Set flag to TRUE, impacts what will be returned from this function
                     done_lat = TRUE

                     # extract dimension information for the grid, note the axis switching between raster and actual array
                     xdim = dim(biomass)[2] ; ydim = dim(biomass)[1]
                     # extract the lat / long information needed
                     long = coordinates(biomass)[,1] ; lat = coordinates(biomass)[,2]
                     # restructure into correct orientation
                     long = array(long, dim=c(xdim,ydim))
                     lat = array(lat, dim=c(xdim,ydim))

                 } # extract lat / long...just the once

                 # break out from the rasters into arrays which we can manipulate
                 biomass = array(as.vector(unlist(biomass)), dim=c(xdim,ydim))
                 biomass_uncertainty = array(as.vector(unlist(biomass_uncertainty)), dim=c(xdim,ydim))

                 # Determine when in the analysis time series the observations should go
                 # NOTE: We assume the biomass estimate is placed at the beginning of the year

                 # What year of the analysis does the data fall?
                 obs_step = which(run_day_selector >= floor(which(analysis_years == years_with_obs[t]) * 365.25))[1]
                 obs_step = obs_step - (steps_per_year-1)

                 # Combine with the other time step
                 if (exists("place_obs_in_step")) {
                     # Output variables already exits to append them
                     place_obs_in_step = append(place_obs_in_step, obs_step)
                     biomass_gCm2 = append(biomass_gCm2, as.vector(biomass)) ; rm(biomass)
                     biomass_uncertainty_gCm2 = append(biomass_uncertainty_gCm2, as.vector(biomass_uncertainty)) ; rm(biomass_uncertainty)
                 } else {
                     # Output variables do not already exist, assign them
                     place_obs_in_step = obs_step ; rm(obs_step)
                     biomass_gCm2 = as.vector(biomass) ; rm(biomass)
                     biomass_uncertainty_gCm2 = as.vector(biomass_uncertainty) ; rm(biomass_uncertainty)
                 } # obs_step exists

             } # Is dataset within the analysis time period?

        } # looping available years

        # Convert MgC/ha -> Mg/ha needed for Saatchi et al (2011)
        biomass_gCm2 = biomass_gCm2 * 1e-2 * 2.083333
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 1e-2 * 2.083333
        # Use allometry to estimate below ground biomass stock and
        # combined with the above ground (Mg/ha) to give a total woody biomass estimate
        # Saatchi et al., (2011), PNAS, 108, 9899-9904, https://www.pnas.org/content/108/24/9899
        biomass_gCm2 = biomass_gCm2 + (0.489 * biomass_gCm2 ** 0.89)
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 + (0.489 * biomass_uncertainty_gCm2 ** 0.89)
        # Convert units of biomass and its uncertainty from MgCha -> gC/m2
        biomass_gCm2 = biomass_gCm2 * 1e2 * 0.48
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 1e2 * 0.48

        # Re-construct arrays for output
        idim = dim(lat)[1] ; jdim = dim(long)[2] ; tdim = length(biomass_gCm2) / (idim * jdim)
        biomass_gCm2 = array(biomass_gCm2, dim=c(idim,jdim,tdim))
        biomass_uncertainty_gCm2 = array(biomass_uncertainty_gCm2, dim=c(idim,jdim,tdim))

        # clean up variables
        gc(reset=TRUE,verbose=FALSE)

        if (done_lat) {
            # Output variables
            return(list(place_obs_in_step = place_obs_in_step, lat = lat, long = long,
                        biomass_gCm2 = biomass_gCm2, biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2))
        } else {
            # Output dummy variables
            return(list(place_obs_in_step = -9999, lat = -9999, long = -9999,
                        biomass_gCm2 = -9999, biomass_uncertainty_gCm2 = -9999))
        } # done_lat

    } else if (Cwood_stock_source == "ESA_CCI_Biomass") {

        # this is a very bespoke modification so leave it here to avoid getting lost
        print("Loading ESA CCI AGB maps")

        # Create the full file paths estimates and their uncertainty (MgC/ha)
        input_file = list.files(path_to_Cwood)
        # extract only .tif files, $ symbol asks for strings that end in the given pattern
        # The \\ also specifies that the . is not to be considered a wildcard
        input_file = input_file[grepl("\\.tif$",input_file) == TRUE]
        # Extract the uncertainty files from the original list
        unc_input_file = input_file[grepl("SD_AGB",input_file) == TRUE]
        input_file = input_file[grepl("SD_AGB",input_file) == FALSE]
        # Check that we have the same number of files for both biomass and uncertainty
        if (length(input_file) != length(unc_input_file)) {stop("Different number of observation and uncertainty files found...")}
        # Determine the number of years found
        years_with_obs = gsub("AGB_map_MgCha_","",input_file)
        years_with_obs = as.numeric(gsub("\\.tif$","",years_with_obs))

        # Loop through each year and extract if appropriate
        done_lat = FALSE
        for (t in seq(1, length(years_with_obs))) {

             # determine whether the first year is within the analysis period
             if (years_with_obs[t] >= as.numeric(start) & years_with_obs[t] <= as.numeric(finish)) {

                 # Read in the estimate and uncertainty rasters
                 biomass = raster(paste(path_to_Cwood,input_file[t],sep=""))
                 biomass_uncertainty = raster(paste(path_to_Cwood,unc_input_file[t],sep=""))

                 # Create raster with the target crs
                 target = raster(crs = ("+init=epsg:4326"), ext = extent(biomass), resolution = res(biomass))
                 # Check whether the target and actual analyses have the same CRS
                 if (compareCRS(biomass,target) == FALSE) {
                     # Resample to correct grid
                     biomass = resample(biomass, target, method="ngb") ; gc() ; removeTmpFiles()
                     biomass_uncertainty = resample(biomass_uncertainty, target, method="ngb") ; gc() ; removeTmpFiles()
                 }
                 # Extend the extent of the overall grid to the analysis domain
                 biomass = extend(biomass,cardamom_ext) ; biomass_uncertainty = extend(biomass_uncertainty,cardamom_ext)
                 # Trim the extent of the overall grid to the analysis domain
                 biomass = crop(biomass,cardamom_ext) ; biomass_uncertainty = crop(biomass_uncertainty,cardamom_ext)
                 # now remove the ones that are actual missing data
                 biomass[which(as.vector(biomass) < 0)] = NA
                 biomass_uncertainty[which(as.vector(biomass_uncertainty) < 0)] = NA
                 # Adjust spatial resolution of the datasets, this occurs in all cases
                 if (res(biomass)[1] != res(cardamom_ext)[1] | res(biomass)[2] != res(cardamom_ext)[2]) {

                     # Create raster with the target resolution
                     target = raster(crs = crs(cardamom_ext), ext = extent(cardamom_ext), resolution = res(cardamom_ext))

                     # Resample to correct grid
                     biomass = resample(biomass, target, method="bilinear") ; gc() ; removeTmpFiles()
                     biomass_uncertainty = resample(biomass_uncertainty, target, method="bilinear") ; gc() ; removeTmpFiles()

                 } # Aggrgeate to resolution

                 # If the first file to be read extract the lat / long information
                 if (done_lat == FALSE) {
                     # Set flag to TRUE, impacts what will be returned from this function
                     done_lat = TRUE

                     # extract dimension information for the grid, note the axis switching between raster and actual array
                     xdim = dim(biomass)[2] ; ydim = dim(biomass)[1]
                     # extract the lat / long information needed
                     long = coordinates(biomass)[,1] ; lat = coordinates(biomass)[,2]
                     # restructure into correct orientation
                     long = array(long, dim=c(xdim,ydim))
                     lat = array(lat, dim=c(xdim,ydim))

                 } # extract lat / long...just the once

                 # break out from the rasters into arrays which we can manipulate
                 biomass = array(as.vector(unlist(biomass)), dim=c(xdim,ydim))
                 biomass_uncertainty = array(as.vector(unlist(biomass_uncertainty)), dim=c(xdim,ydim))

                 # Determine when in the analysis time series the observations should go
                 # NOTE: We assume the biomass estimate is placed at the beginning of the year

                 # What year of the analysis does the data fall?
                 obs_step = which(run_day_selector >= floor(which(analysis_years == years_with_obs[t]) * 365.25))[1]
                 obs_step = obs_step - (steps_per_year-1)
                 # Combine with the other time step
                 if (exists("place_obs_in_step")) {
                     # Output variables already exits to append them
                     place_obs_in_step = append(place_obs_in_step, obs_step)
                     biomass_gCm2 = append(biomass_gCm2, as.vector(biomass)) ; rm(biomass)
                     biomass_uncertainty_gCm2 = append(biomass_uncertainty_gCm2, as.vector(biomass_uncertainty)) ; rm(biomass_uncertainty)
                 } else {
                     # Output variables do not already exist, assign them
                     place_obs_in_step = obs_step ; rm(obs_step)
                     biomass_gCm2 = as.vector(biomass) ; rm(biomass)
                     biomass_uncertainty_gCm2 = as.vector(biomass_uncertainty) ; rm(biomass_uncertainty)
                 } # obs_step exists

             } # Is dataset within the analysis time period?

        } # looping available years

        # Convert MgC/ha -> Mg/ha needed for Saatchi et al (2011)
        biomass_gCm2 = biomass_gCm2 * 2.083333
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 2.083333
        # Use allometry to estimate below ground biomass stock and
        # combined with the above ground (Mg/ha) to give a total woody biomass estimate
        # Saatchi et al., (2011), PNAS, 108, 9899-9904, https://www.pnas.org/content/108/24/9899
        biomass_gCm2 = biomass_gCm2 + (0.489 * biomass_gCm2 ** 0.89)
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 + (0.489 * biomass_uncertainty_gCm2 ** 0.89)
        # Now back to desired units gC/m2
        biomass_gCm2 = biomass_gCm2 * 0.48 * 1e2
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 0.48 * 1e2

        # Re-construct arrays for output
        idim = dim(lat)[1] ; jdim = dim(long)[2] ; tdim = length(biomass_gCm2) / (idim * jdim)
        biomass_gCm2 = array(biomass_gCm2, dim=c(idim,jdim,tdim))
        biomass_uncertainty_gCm2 = array(biomass_uncertainty_gCm2, dim=c(idim,jdim,tdim))
        if (done_lat) {
            # Output variables
            return(list(place_obs_in_step = place_obs_in_step, lat = lat, long = long,
                        biomass_gCm2 = biomass_gCm2, biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2))
        } else {
            # Output dummy variables
            return(list(place_obs_in_step = -9999, lat = -9999, long = -9999,
                        biomass_gCm2 = -9999, biomass_uncertainty_gCm2 = -9999))
        } # done_lat

    } else if (Cwood_stock_source == "Saatchi_2021") {

        # this is a very bespoke modification so leave it here to avoid getting lost
        print("Loading Saatchi 2021 wood stock maps")

        # Create the full file paths estimates and their uncertainty (MgC/ha)
        input_file = list.files(path_to_Cwood)
        # extract only .tif files, $ symbol asks for strings that end in the given pattern
        # The \\ also specifies that the . is not to be considered a wildcard
        input_file = input_file[grepl("\\.tif$",input_file) == TRUE]
        # Extract the specific files from the original list
        input_unc_file = input_file[grepl("saatchi_wood_uncertainty_MgCha",input_file) == TRUE]
        # Extract the specific files from the original list
        input_file = input_file[grepl("saatchi_wood_MgCha",input_file) == TRUE]

        # Sense check
        if (length(input_unc_file) != length(input_file)) {stop("number of uncertainty and data files differ for Cwood_stock_source = Saatchi_2021")}

        # Determine the number of years found
        years_with_obs = gsub("saatchi_wood_MgCha_","",input_file)
        years_with_obs = as.numeric(gsub("\\.tif$","",years_with_obs))

        # Loop through each year and extract if appropriate
        done_lat = FALSE
        for (t in seq(1, length(years_with_obs))) {

             # determine whether the first year is within the analysis period
             if (years_with_obs[t] >= as.numeric(start) & years_with_obs[t] <= as.numeric(finish)) {

                 # Read in the estimate and uncertainty rasters
                 biomass = raster(paste(path_to_Cwood,input_file[t],sep=""))
                 biomass_unc = raster(paste(path_to_Cwood,input_unc_file[t],sep=""))

                 # Create raster with the target crs
                 target = raster(crs = ("+init=epsg:4326"), ext = extent(biomass), resolution = res(biomass))
                 # Check whether the target and actual analyses have the same CRS
                 if (compareCRS(biomass,target) == FALSE) {
                     # Resample to correct grid
                     biomass = resample(biomass, target, method="ngb") ; gc() ; removeTmpFiles()
                     biomass_unc = resample(biomass_unc, target, method="ngb") ; gc() ; removeTmpFiles()
                 }
                 # Extend the extent of the overall grid to the analysis domain
                 biomass = extend(biomass,cardamom_ext)
                 biomass_unc = extend(biomass_unc,cardamom_ext)
                 # Trim the extent of the overall grid to the analysis domain
                 biomass = crop(biomass,cardamom_ext)
                 biomass_unc = crop(biomass_unc,cardamom_ext)
                 # now remove the ones that are actual missing data
                 biomass[which(as.vector(biomass) < 0)] = NA
                 biomass_unc[which(as.vector(biomass_unc) < 0)] = NA

                 # Adjust spatial resolution of the datasets, this occurs in all cases
                 if (res(biomass)[1] != res(cardamom_ext)[1] | res(biomass)[2] != res(cardamom_ext)[2]) {

                     # Create raster with the target resolution
                     target = raster(crs = crs(cardamom_ext), ext = extent(cardamom_ext), resolution = res(cardamom_ext))
                     # Resample to correct grid
                     biomass = resample(biomass, target, method="bilinear") ; gc() ; removeTmpFiles()
                     biomass_unc = resample(biomass_unc, target, method="bilinear") ; gc() ; removeTmpFiles()

                 } # Aggrgeate to resolution


                 # If the first file to be read extract the lat / long information
                 if (done_lat == FALSE) {
                     # Set flag to TRUE, impacts what will be returned from this function
                     done_lat = TRUE

                     # extract dimension information for the grid, note the axis switching between raster and actual array
                     xdim = dim(biomass)[2] ; ydim = dim(biomass)[1]
                     # extract the lat / long information needed
                     long = coordinates(biomass)[,1] ; lat = coordinates(biomass)[,2]
                     # restructure into correct orientation
                     long = array(long, dim=c(xdim,ydim))
                     lat = array(lat, dim=c(xdim,ydim))

                 } # extract lat / long...just the once

                 # break out from the rasters into arrays which we can manipulate
                 biomass = array(as.vector(unlist(biomass)), dim=c(xdim,ydim))
                 biomass_unc = array(as.vector(unlist(biomass_unc)), dim=c(xdim,ydim))

                 # Determine when in the analysis time series the observations should go
                 # NOTE: We assume the biomass estimate is placed at the beginning of the year

                 # What year of the analysis does the data fall?
                 obs_step = which(run_day_selector >= floor(which(analysis_years == years_with_obs[t]) * 365.25))[1]
                 obs_step = obs_step - (steps_per_year-1)
                 # Combine with the other time step
                 if (exists("place_obs_in_step")) {
                     # Output variables already exits to append them
                     place_obs_in_step = append(place_obs_in_step, obs_step)
                     biomass_gCm2 = append(biomass_gCm2, as.vector(biomass)) ; rm(biomass)
                     biomass_uncertainty_gCm2 = append(biomass_uncertainty_gCm2, as.vector(biomass_unc)) ; rm(biomass_unc)
                 } else {
                     # Output variables do not already exist, assign them
                     place_obs_in_step = obs_step ; rm(obs_step)
                     biomass_gCm2 = as.vector(biomass) ; rm(biomass)
                     biomass_uncertainty_gCm2 = as.vector(biomass_unc) ; rm(biomass_unc)
                 } # obs_step exists

             } # Is dataset within the analysis time period?

        } # looping available years

        # Convert MgC/ha -> gCm2 needed for Saatchi et al (2011)
        biomass_gCm2 = biomass_gCm2 * 1e2
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 1e2

        # Re-construct arrays for output
        idim = dim(lat)[1] ; jdim = dim(long)[2] ; tdim = length(biomass_gCm2) / (idim * jdim)
        biomass_gCm2 = array(biomass_gCm2, dim=c(idim,jdim,tdim))
        biomass_uncertainty_gCm2 = array(biomass_uncertainty_gCm2, dim=c(idim,jdim,tdim))
        #biomass_uncertainty_gCm2 = biomass_gCm2 * 0.18 # assume uncertainty is 18 %
        if (done_lat) {
            # Output variables
            return(list(place_obs_in_step = place_obs_in_step, lat = lat, long = long,
                        biomass_gCm2 = biomass_gCm2, biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2))
        } else {
            # Output dummy variables
            return(list(place_obs_in_step = -9999, lat = -9999, long = -9999,
                        biomass_gCm2 = -9999, biomass_uncertainty_gCm2 = -9999))
        } # done_lat

    } else if (Cwood_stock_source == "UoL_stable_forest") {

        # this is a very bespoke modification so leave it here to avoid getting lost
        print("Loading UoL Stable Forests AGB...")

        # Data represent the following year
        years_with_obs = c(2014) ; done_lat = FALSE
        # determine whether the first year is within the analysis period
        if (years_with_obs >= as.numeric(start) & years_with_obs <= as.numeric(finish)) {

            done_lat = TRUE

            # Read in estimate and uncertainty rasters
            # NOTE: that map is above ground, total biomass will be estimated later
            biomass_gCm2 = raster(paste(path_to_Cwood,"Kenya_0.25deg_AGB_stable_forest_2015_2017.tif", sep=""))
            biomass_uncertainty_gCm2 = raster(paste(path_to_Cwood,"Kenya_0.25deg_AGB_std_stable_forest_2015_2017.tif", sep=""))

            # Create raster with the target crs
            target = raster(crs = ("+init=epsg:4326"), ext = extent(biomass_gCm2), resolution = res(biomass_gCm2))
            # Check whether the target and actual analyses have the same CRS
            if (compareCRS(biomass_gCm2,target) == FALSE) {
                # Resample to correct grid
                biomass_gCm2 = resample(biomass_gCm2, target, method="ngb") ; gc() ; removeTmpFiles()
                biomass_uncertainty_gCm2 = resample(biomass_uncertainty_gCm2, target, method="ngb") ; gc() ; removeTmpFiles()
            }
            # Extend the extent of the overall grid to the analysis domain
            biomass_gCm2 = extend(biomass_gCm2,cardamom_ext) ; biomass_uncertainty_gCm2 = extend(biomass_uncertainty_gCm2,cardamom_ext)
            # Trim the extent of the overall grid to the analysis domain
            biomass_gCm2 = crop(biomass_gCm2,cardamom_ext) ; biomass_uncertainty_gCm2 = crop(biomass_uncertainty_gCm2,cardamom_ext)
            # now remove the ones that are actual missing data
            biomass_gCm2[which(as.vector(biomass_gCm2) < 0)] = NA
            biomass_uncertainty_gCm2[which(as.vector(biomass_uncertainty_gCm2) < 0)] = NA
            # Adjust spatial resolution of the datasets, this occurs in all cases
            if (res(biomass_gCm2)[1] != res(cardamom_ext)[1] | res(biomass_gCm2)[2] != res(cardamom_ext)[2]) {

                # Create raster with the target resolution
                target = raster(crs = crs(cardamom_ext), ext = extent(cardamom_ext), resolution = res(cardamom_ext))
                # Resample to correct grid
                biomass_gCm2 = resample(biomass_gCm2, target, method="bilinear") ; gc() ; removeTmpFiles()
                biomass_uncertainty_gCm2 = resample(biomass_uncertainty_gCm2, target, method="bilinear") ; gc() ; removeTmpFiles()

            } # Aggrgeate to resolution

            # extract dimension information for the grid, note the axis switching between raster and actual array
            xdim = dim(biomass_gCm2)[2] ; ydim = dim(biomass_gCm2)[1]
            # extract the lat / long information needed
            long = coordinates(biomass_gCm2)[,1] ; lat = coordinates(biomass_gCm2)[,2]
            # restructure into correct orientation
            long = array(long, dim=c(xdim,ydim))
            lat = array(lat, dim=c(xdim,ydim))
            # break out from the rasters into arrays which we can manipulate
            biomass_gCm2 = array(as.vector(unlist(biomass_gCm2)), dim=c(xdim,ydim))
            biomass_uncertainty_gCm2 = array(as.vector(unlist(biomass_uncertainty_gCm2)), dim=c(xdim,ydim))

            # Determine when in the analysis time series the observations should go
            # NOTE: We assume the biomass estimate is placed at the beginning of the year

            # What year of the analysis does the data fall?
            place_obs_in_step = which(run_day_selector >= floor(which(analysis_years == years_with_obs) * 365.25))[1]
            place_obs_in_step = place_obs_in_step - (steps_per_year-1)

            # Re-construct arrays for output
            idim = dim(lat)[1] ; jdim = dim(long)[2] ; tdim = length(biomass_gCm2) / (idim * jdim)
            biomass_gCm2 = array(biomass_gCm2, dim=c(idim,jdim,tdim))
            biomass_uncertainty_gCm2 = array(biomass_uncertainty_gCm2, dim=c(idim,jdim,tdim))

       } # data within analysis period

        # Convert gC/m2 -> Mg/ha needed for Saatchi et al (2011)
        biomass_gCm2 = biomass_gCm2 * 2.083333 * 1e-2
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 2.083333 * 1e-2
        # Use allometry to estimate below ground biomass stock and
        # combined with the above ground (Mg/ha) to give a total woody biomass estimate
        # Saatchi et al., (2011), PNAS, 108, 9899-9904, https://www.pnas.org/content/108/24/9899
        biomass_gCm2 = biomass_gCm2 + (0.489 * biomass_gCm2 ** 0.89)
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 + (0.489 * biomass_uncertainty_gCm2 ** 0.89)
        # Now back to desired units gC/m2
        biomass_gCm2 = biomass_gCm2 * 0.48 * 1e2
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 0.48 * 1e2

        if (done_lat) {
            # Output variables
            return(list(place_obs_in_step = place_obs_in_step, lat = lat, long = long,
                        biomass_gCm2 = biomass_gCm2, biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2))
        } else {
            # Output dummy variables
            return(list(place_obs_in_step = -9999, lat = -9999, long = -9999,
                        biomass_gCm2 = -9999, biomass_uncertainty_gCm2 = -9999))
        } # done_lat

    } else if (Cwood_stock_source == "UoL_stable_savannah") {

        # this is a very bespoke modification so leave it here to avoid getting lost
        print("Loading UoL Stable Savannah AGB...")

        years_with_obs = c(2014) ; done_lat = FALSE
        # determine whether the first year is within the analysis period
        if (years_with_obs >= as.numeric(start) & years_with_obs <= as.numeric(finish)) {

            done_lat = TRUE

            # Read in estimate and uncertainty rasters
            # NOTE: that map is above ground, total biomass will be estimated later
            biomass_gCm2 = raster(paste(path_to_Cwood,"Kenya_0.25deg_AGB_stable_savannah_2015_2017.tif", sep=""))
            biomass_uncertainty_gCm2 = raster(paste(path_to_Cwood,"Kenya_0.25deg_AGB_std_stable_savannah_2015_2017.tif", sep=""))

            # Create raster with the target crs
            target = raster(crs = ("+init=epsg:4326"), ext = extent(biomass_gCm2), resolution = res(biomass_gCm2))
            # Check whether the target and actual analyses have the same CRS
            if (compareCRS(biomass_gCm2,target) == FALSE) {
                # Resample to correct grid
                biomass_gCm2 = resample(biomass_gCm2, target, method="ngb") ; gc() ; removeTmpFiles()
                biomass_uncertainty_gCm2 = resample(biomass_uncertainty_gCm2, target, method="ngb") ; gc() ; removeTmpFiles()
            }
            # Extend the extent of the overall grid to the analysis domain
            biomass_gCm2 = extend(biomass_gCm2,cardamom_ext) ; biomass_uncertainty_gCm2 = extend(biomass_uncertainty_gCm2,cardamom_ext)
            # Trim the extent of the overall grid to the analysis domain
            biomass_gCm2 = crop(biomass_gCm2,cardamom_ext) ; biomass_uncertainty_gCm2 = crop(biomass_uncertainty_gCm2,cardamom_ext)
            # now remove the ones that are actual missing data
            biomass_gCm2[which(as.vector(biomass_gCm2) < 0)] = NA
            biomass_uncertainty_gCm2[which(as.vector(biomass_uncertainty_gCm2) < 0)] = NA
            # Adjust spatial resolution of the datasets, this occurs in all cases
            if (res(biomass_gCm2)[1] != res(cardamom_ext)[1] | res(biomass_gCm2)[2] != res(cardamom_ext)[2]) {

                 # Create raster with the target resolution
                 target = raster(crs = crs(cardamom_ext), ext = extent(cardamom_ext), resolution = res(cardamom_ext))

                 # Resample to correct grid
                 biomass_gCm2 = resample(biomass_gCm2, target, method="bilinear") ; gc() ; removeTmpFiles()
                 biomass_uncertainty_gCm2 = resample(biomass_uncertainty_gCm2, target, method="bilinear") ; gc() ; removeTmpFiles()

            } # Aggrgeate to resolution

            # extract dimension information for the grid, note the axis switching between raster and actual array
            xdim = dim(biomass_gCm2)[2] ; ydim = dim(biomass_gCm2)[1]
            # extract the lat / long information needed
            long = coordinates(biomass_gCm2)[,1] ; lat = coordinates(biomass_gCm2)[,2]
            # restructure into correct orientation
            long = array(long, dim=c(xdim,ydim))
            lat = array(lat, dim=c(xdim,ydim))
            # break out from the rasters into arrays which we can manipulate
            biomass_gCm2 = array(as.vector(unlist(biomass_gCm2)), dim=c(xdim,ydim))
            biomass_uncertainty_gCm2 = array(as.vector(unlist(biomass_uncertainty_gCm2)), dim=c(xdim,ydim))

            # Determine when in the analysis time series the observations should go
            # NOTE: We assume the biomass estimate is placed at the beginning of the year

            # What year of the analysis does the data fall?
            place_obs_in_step = which(run_day_selector >= floor(which(analysis_years == years_with_obs) * 365.25))[1]
            place_obs_in_step = place_obs_in_step - (steps_per_year-1)

            # Re-construct arrays for output
            idim = dim(lat)[1] ; jdim = dim(long)[2] ; tdim = length(biomass_gCm2) / (idim * jdim)
            biomass_gCm2 = array(biomass_gCm2, dim=c(idim,jdim,tdim))
            biomass_uncertainty_gCm2 = array(biomass_uncertainty_gCm2, dim=c(idim,jdim,tdim))

       } # data within analysis period

        # Convert gC/m2 -> Mg/ha needed for Saatchi et al (2011)
        biomass_gCm2 = biomass_gCm2 * 2.083333 * 1e-2
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 2.083333 * 1e-2
        # Use allometry to estimate below ground biomass stock and
        # combined with the above ground (Mg/ha) to give a total woody biomass estimate
        # Saatchi et al., (2011), PNAS, 108, 9899-9904, https://www.pnas.org/content/108/24/9899
        biomass_gCm2 = biomass_gCm2 + (0.489 * biomass_gCm2 ** 0.89)
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 + (0.489 * biomass_uncertainty_gCm2 ** 0.89)
        # Now back to desired units gC/m2
        biomass_gCm2 = biomass_gCm2 * 0.48 * 1e2
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 0.48 * 1e2

       if (done_lat) {
           # Output variables
           return(list(place_obs_in_step = place_obs_in_step, lat = lat, long = long,
                       biomass_gCm2 = biomass_gCm2, biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2))
       } else {
           # Output dummy variables
           return(list(place_obs_in_step = -9999, lat = -9999, long = -9999,
                       biomass_gCm2 = -9999, biomass_uncertainty_gCm2 = -9999))
       } # done_lat

    } else if (Cwood_stock_source == "Rainfor") {

      # this is a very bespoke modification so leave it here to avoid getting lost
      print("Loading UoL Rainfor AGB...")

      # Data represent the following year
      years_with_obs = c(2010) ; done_lat = FALSE
      # determine whether the first year is within the analysis period
      if (years_with_obs >= as.numeric(start) & years_with_obs <= as.numeric(finish)) {

        done_lat = TRUE

        # Read in estimate and uncertainty rasters
        # NOTE: this data assimilates total biomass
        biomass_gCm2 = raster(paste(path_to_Cwood,"wood_biomass_gCm2_2010.tif", sep=""))
        biomass_uncertainty_gCm2 = raster(paste(path_to_Cwood,"unc_wood_biomass_gCm2_2010.tif", sep=""))

        # Create raster with the target crs
        target = raster(crs = ("+init=epsg:4326"), ext = extent(biomass_gCm2), resolution = res(biomass_gCm2))
        # Check whether the target and actual analyses have the same CRS
        if (compareCRS(biomass_gCm2,target) == FALSE) {
          # Resample to correct grid
          biomass_gCm2 = resample(biomass_gCm2, target, method="ngb") ; gc() ; removeTmpFiles()
          biomass_uncertainty_gCm2 = resample(biomass_uncertainty_gCm2, target, method="ngb") ; gc() ; removeTmpFiles()
        }
        # Extend the extent of the overall grid to the analysis domain
        biomass_gCm2 = extend(biomass_gCm2,cardamom_ext) ; biomass_uncertainty_gCm2 = extend(biomass_uncertainty_gCm2,cardamom_ext)
        # Trim the extent of the overall grid to the analysis domain
        biomass_gCm2 = crop(biomass_gCm2,cardamom_ext) ; biomass_uncertainty_gCm2 = crop(biomass_uncertainty_gCm2,cardamom_ext)
        # now remove the ones that are actual missing data
        biomass_gCm2[which(as.vector(biomass_gCm2) < 0)] = NA
        biomass_uncertainty_gCm2[which(as.vector(biomass_uncertainty_gCm2) < 0)] = NA
        # Adjust spatial resolution of the datasets, this occurs in all cases
        if (res(biomass_gCm2)[1] != res(cardamom_ext)[1] | res(biomass_gCm2)[2] != res(cardamom_ext)[2]) {

            # Create raster with the target resolution
            target = raster(crs = crs(cardamom_ext), ext = extent(cardamom_ext), resolution = res(cardamom_ext))

            # Resample to correct grid
            biomass_gCm2 = resample(biomass_gCm2, target, method="bilinear") ; gc() ; removeTmpFiles()
            biomass_uncertainty_gCm2 = resample(biomass_uncertainty_gCm2, target, method="bilinear") ; gc() ; removeTmpFiles()

        } # Aggrgeate to resolution

        # extract dimension information for the grid, note the axis switching between raster and actual array
        xdim = dim(biomass_gCm2)[2] ; ydim = dim(biomass_gCm2)[1]
        # extract the lat / long information needed
        long = coordinates(biomass_gCm2)[,1] ; lat = coordinates(biomass_gCm2)[,2]
        # restructure into correct orientation
        long = array(long, dim=c(xdim,ydim))
        lat = array(lat, dim=c(xdim,ydim))
        # break out from the rasters into arrays which we can manipulate
        biomass_gCm2 = array(as.vector(unlist(biomass_gCm2)), dim=c(xdim,ydim))
        biomass_uncertainty_gCm2 = array(as.vector(unlist(biomass_uncertainty_gCm2)), dim=c(xdim,ydim))

        # Determine when in the analysis time series the observations should go
        # NOTE: We assume the biomass estimate is placed at the beginning of the year

        # What year of the analysis does the data fall?
        place_obs_in_step = which(run_day_selector >= floor(which(analysis_years == years_with_obs) * 365.25))[1]
        place_obs_in_step = place_obs_in_step - (steps_per_year-1)

        # Re-construct arrays for output
        idim = dim(lat)[1] ; jdim = dim(long)[2] ; tdim = length(biomass_gCm2) / (idim * jdim)
        biomass_gCm2 = array(biomass_gCm2, dim=c(idim,jdim,tdim))
        biomass_uncertainty_gCm2 = array(biomass_uncertainty_gCm2, dim=c(idim,jdim,tdim))

      } # data within analysis period

      # # Convert gC/m2 -> Mg/ha needed for Saatchi et al (2011)
      # biomass_gCm2 = biomass_gCm2 * 2.083333 * 1e-2
      # biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 2.083333 * 1e-2
      # # Use allometry to estimate below ground biomass stock and
      # # combined with the above ground (Mg/ha) to give a total woody biomass estimate
      # # Saatchi et al., (2011), PNAS, 108, 9899-9904, https://www.pnas.org/content/108/24/9899
      # biomass_gCm2 = biomass_gCm2 + (0.489 * biomass_gCm2 ** 0.89)
      # biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 + (0.489 * biomass_uncertainty_gCm2 ** 0.89)
      # # Now back to desired units gC/m2
      # biomass_gCm2 = biomass_gCm2 * 0.48 * 1e2
      # biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 0.48 * 1e2

      if (done_lat) {
        # Output variables
        return(list(place_obs_in_step = place_obs_in_step, lat = lat, long = long,
                    biomass_gCm2 = biomass_gCm2, biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2))
      } else {
        # Output dummy variables
        return(list(place_obs_in_step = -9999, lat = -9999, long = -9999,
                    biomass_gCm2 = -9999, biomass_uncertainty_gCm2 = -9999))
      } # done_lat

    } else if (Cwood_stock_source == "Rainfor_annual") {

        # this is a very bespoke modification so leave it here to avoid getting lost
        print("Loading UoL Rainfor biomass annually...")

        # Create the full file paths estimates and their uncertainty (MgC/ha)
        input_file = list.files(path_to_Cwood)
        # extract only .tif files, $ symbol asks for strings that end in the given pattern
        # The \\ also specifies that the . is not to be considered a wildcard
        input_file = input_file[grepl("\\.tif$",input_file) == TRUE]
        # Extract the uncertainty files from the original list
        unc_input_file = input_file[grepl("unc_wood",input_file) == TRUE]
        input_file = input_file[grepl("unc_wood",input_file) == FALSE]
        # Check that we have the same number of files for both biomass and uncertainty
        if (length(input_file) != length(unc_input_file)) {stop("Different number of observation and uncertainty files found...")}
        # Determine the number of years found
        years_with_obs = gsub("wood_biomass_gCm2_","",input_file)
        years_with_obs = as.numeric(gsub("\\.tif$","",years_with_obs))
        print(length(years_with_obs))
        if (length(years_with_obs) == length(unc_input_file)) {print("Same number of years and uncertainty files...")}
        # Loop through each year and extract if appropriate
        done_lat = FALSE
        for (t in seq(1, length(years_with_obs))) {

        # determine whether the first year is within the analysis period
        if (years_with_obs[t] >= as.numeric(start) & years_with_obs[t] <= as.numeric(finish)) {

            # Read in the estimate and uncertainty rasters
            biomass = raster(paste(path_to_Cwood,input_file[t],sep=""))
            biomass_uncertainty = raster(paste(path_to_Cwood,unc_input_file[t],sep=""))
            print(summary(biomass)[3])
            # Create raster with the target crs
            target = raster(crs = ("+init=epsg:4326"), ext = extent(biomass), resolution = res(biomass))
            # Check whether the target and actual analyses have the same CRS
            if (compareCRS(biomass,target) == FALSE) {
                # Resample to correct grid
                biomass = resample(biomass, target, method="ngb") ; gc() ; removeTmpFiles()
                biomass_uncertainty = resample(biomass_uncertainty, target, method="ngb") ; gc() ; removeTmpFiles()
            }
            # Extend the extent of the overall grid to the analysis domain
            biomass = extend(biomass,cardamom_ext) ; biomass_uncertainty = extend(biomass_uncertainty,cardamom_ext)
            # Trim the extent of the overall grid to the analysis domain
            biomass = crop(biomass,cardamom_ext) ; biomass_uncertainty = crop(biomass_uncertainty,cardamom_ext)
            # now remove the ones that are actual missing data
            biomass[which(as.vector(biomass) < 0)] = NA
            biomass_uncertainty[which(as.vector(biomass_uncertainty) < 0)] = NA

            # Adjust spatial resolution of the datasets, this occurs in all cases
            if (res(biomass)[1] != res(cardamom_ext)[1] | res(biomass)[2] != res(cardamom_ext)[2]) {
                # Create raster with the target resolution
                target = raster(crs = crs(cardamom_ext), ext = extent(cardamom_ext), resolution = res(cardamom_ext))
                # Resample to correct grid
                biomass = resample(biomass, target, method="bilinear") ; gc() ; removeTmpFiles()
                biomass_uncertainty = resample(biomass_uncertainty, target, method="bilinear") ; gc() ; removeTmpFiles()
            } # Aggrgeate to resolution

            # If the first file to be read extract the lat / long information
            if (done_lat == FALSE) {
                # Set flag to TRUE, impacts what will be returned from this function
                done_lat = TRUE

                # extract dimension information for the grid, note the axis switching between raster and actual array
                xdim = dim(biomass)[2] ; ydim = dim(biomass)[1]
                # extract the lat / long information needed
                long = coordinates(biomass)[,1] ; lat = coordinates(biomass)[,2]
                # restructure into correct orientation
                long = array(long, dim=c(xdim,ydim))
                lat = array(lat, dim=c(xdim,ydim))
            } # extract lat / long...just the once

            # break out from the rasters into arrays which we can manipulate
            biomass = array(as.vector(unlist(biomass)), dim=c(xdim,ydim))
            biomass_uncertainty = array(as.vector(unlist(biomass_uncertainty)), dim=c(xdim,ydim))
            # Determine when in the analysis time series the observations should go
            # NOTE: We assume the biomass estimate is placed at the beginning of the year

            # What year of the analysis does the data fall?
            obs_step = which(run_day_selector >= floor(which(analysis_years == years_with_obs[t]) * 365.25))[1]
            obs_step = obs_step - (steps_per_year-1)
            # Combine with the other time step
            if (exists("place_obs_in_step")) {
                # Output variables already exits to append them
                place_obs_in_step = append(place_obs_in_step, obs_step)
                biomass_gCm2 = append(biomass_gCm2, as.vector(biomass)) ; rm(biomass)
                biomass_uncertainty_gCm2 = append(biomass_uncertainty_gCm2, as.vector(biomass_uncertainty)) ; rm(biomass_uncertainty)
            } else {
                # Output variables do not already exist, assign them
                place_obs_in_step = obs_step ; rm(obs_step)
                biomass_gCm2 = as.vector(biomass) ; rm(biomass)
                biomass_uncertainty_gCm2 = as.vector(biomass_uncertainty) ; rm(biomass_uncertainty)
            } # obs_step exists

         } # Is dataset within the analysis time period?

        } # looping available years

        # Re-construct arrays for output
        idim = dim(lat)[1] ; jdim = dim(long)[2] ; tdim = length(biomass_gCm2) / (idim * jdim)
        biomass_gCm2 = array(biomass_gCm2, dim=c(idim,jdim,tdim))
        biomass_uncertainty_gCm2 = array(biomass_uncertainty_gCm2, dim=c(idim,jdim,tdim))

        if (done_lat) {
            # Output variables
            return(list(place_obs_in_step = place_obs_in_step, lat = lat, long = long,
                        biomass_gCm2 = biomass_gCm2, biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2))
        } else {
            # Output dummy variables
            return(list(place_obs_in_step = -9999, lat = -9999, long = -9999,
                        biomass_gCm2 = -9999, biomass_uncertainty_gCm2 = -9999))
        } # done_lat

    } else if (Cwood_stock_source == "Biomass_maps_Africa_UoL") {

        # this is a very bespoke modification so leave it here to avoid getting lost
        print("Loading UoL Africa maps")

        # Create the full file paths estimates and their uncertainty (MgC/ha)
        input_file = list.files(path_to_Cwood)
        # extract only .tif files
        input_file = input_file[grepl(".tif",input_file) == TRUE]
        # Check that we have found some files
        if (length(input_file) == 0) {stop("No files have been found in target directory")}
        # Extract the uncertainty files from the original list
        unc_input_file = input_file[grepl("SD_AGB",input_file) == TRUE]
        input_file = input_file[grepl("SD_AGB",input_file) == FALSE]
        # Check that we have the same number of files for both biomass and uncertainty
        if (length(input_file) != length(unc_input_file)) {stop("Different number of observation and uncertainty files found...")}
        # Determine the number of years found
        years_with_obs = gsub("AGB_map_combined_","",input_file)
        years_with_obs = as.numeric(gsub(".tif","",years_with_obs))

        # Loop through each year and extract if appropriate
        done_lat = FALSE
        for (t in seq(1, length(years_with_obs))) {

             # determine whether the first year is within the analysis period
             if (years_with_obs[t] >= as.numeric(start) & years_with_obs[t] <= as.numeric(finish)) {

                 # Read in the estimate and uncertainty rasters
                 biomass = raster(paste(path_to_Cwood,input_file[t],sep=""))
                 biomass_uncertainty = raster(paste(path_to_Cwood,unc_input_file[t],sep=""))

                 # Create raster with the target crs
                 target = raster(crs = ("+init=epsg:4326"), ext = extent(biomass), resolution = res(biomass))
                 # Check whether the target and actual analyses have the same CRS
                 if (compareCRS(biomass,target) == FALSE) {
                     # Resample to correct grid
                     biomass = resample(biomass, target, method="ngb") ; gc() ; removeTmpFiles()
                     biomass_uncertainty = resample(biomass_uncertainty, target, method="ngb") ; gc() ; removeTmpFiles()
                 }
                 # Extend the extent of the overall grid to the analysis domain
                 biomass = extend(biomass,cardamom_ext) ; biomass_uncertainty = extend(biomass_uncertainty,cardamom_ext)
                 # Trim the extent of the overall grid to the analysis domain
                 biomass = crop(biomass,cardamom_ext) ; biomass_uncertainty = crop(biomass_uncertainty,cardamom_ext)
                 # now remove the ones that are actual missing data
                 biomass[which(as.vector(biomass) < 0)] = NA
                 biomass_uncertainty[which(as.vector(biomass_uncertainty) < 0)] = NA
                 # Adjust spatial resolution of the datasets, this occurs in all cases
                 if (res(biomass)[1] != res(cardamom_ext)[1] | res(biomass)[2] != res(cardamom_ext)[2]) {

                     # Create raster with the target resolution
                     target = raster(crs = crs(cardamom_ext), ext = extent(cardamom_ext), resolution = res(cardamom_ext))

                     # Resample to correct grid
                     biomass = resample(biomass, target, method="bilinear") ; gc() ; removeTmpFiles()
                     biomass_uncertainty = resample(biomass_uncertainty, target, method="bilinear") ; gc() ; removeTmpFiles()

                 } # Aggrgeate to resolution

                 # If the first file to be read extract the lat / long information
                 if (done_lat == FALSE) {
                     # Set flag to TRUE, impacts what will be returned from this function
                     done_lat = TRUE

                     # extract dimension information for the grid, note the axis switching between raster and actual array
                     xdim = dim(biomass)[2] ; ydim = dim(biomass)[1]
                     # extract the lat / long information needed
                     long = coordinates(biomass)[,1] ; lat = coordinates(biomass)[,2]
                     # restructure into correct orientation
                     long = array(long, dim=c(xdim,ydim))
                     lat = array(lat, dim=c(xdim,ydim))

                 } # extract lat / long...just the once

                 # break out from the rasters into arrays which we can manipulate
                 biomass = array(as.vector(unlist(biomass)), dim=c(xdim,ydim))
                 biomass_uncertainty = array(as.vector(unlist(biomass_uncertainty)), dim=c(xdim,ydim))

                 # Determine when in the analysis time series the observations should go
                 # NOTE: We assume the biomass estimate is placed at the beginning of the year

                 # What year of the analysis does the data fall?
                 obs_step = which(run_day_selector >= floor(which(analysis_years == years_with_obs[t]) * 365.25))[1]
                 obs_step = obs_step - (steps_per_year-1)

                 # Combine with the other time step
                 if (exists("place_obs_in_step")) {
                     # Output variables already exits to append them
                     place_obs_in_step = append(place_obs_in_step, obs_step)
                     biomass_gCm2 = append(biomass_gCm2, as.vector(biomass)) ; rm(biomass)
                     biomass_uncertainty_gCm2 = append(biomass_uncertainty_gCm2, as.vector(biomass_uncertainty)) ; rm(biomass_uncertainty)
                 } else {
                     # Output variables do not already exist, assign them
                     place_obs_in_step = obs_step ; rm(obs_step)
                     biomass_gCm2 = as.vector(biomass) ; rm(biomass)
                     biomass_uncertainty_gCm2 = as.vector(biomass_uncertainty) ; rm(biomass_uncertainty)
                 } # obs_step exists

             } # Is dataset within the analysis time period?

        } # looping available years

        # Convert MgC/ha -> Mg/ha needed for Saatchi et al (2011)
        biomass_gCm2 = biomass_gCm2 * 2.083333
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 2.083333
        # Use allometry to estimate below ground biomass stock and
        # combined with the above ground (Mg/ha) to give a total woody biomass estimate
        # Saatchi et al., (2011), PNAS, 108, 9899-9904, https://www.pnas.org/content/108/24/9899
        biomass_gCm2 = biomass_gCm2 + (0.489 * biomass_gCm2 ** 0.89)
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 + (0.489 * biomass_uncertainty_gCm2 ** 0.89)
        # Now back to desired units gC/m2
        biomass_gCm2 = biomass_gCm2 * 0.48 * 1e2
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 0.48 * 1e2

        # Re-construct arrays for output
        idim = dim(lat)[1] ; jdim = dim(long)[2] ; tdim = length(biomass_gCm2) / (idim * jdim)
        biomass_gCm2 = array(biomass_gCm2, dim=c(idim,jdim,tdim))
        biomass_uncertainty_gCm2 = array(biomass_uncertainty_gCm2, dim=c(idim,jdim,tdim))

        if (done_lat) {
            # Output variables
            return(list(place_obs_in_step = place_obs_in_step, lat = lat, long = long,
                        biomass_gCm2 = biomass_gCm2, biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2))
        } else {
            # Output dummy variables
            return(list(place_obs_in_step = -9999, lat = -9999, long = -9999,
                        biomass_gCm2 = -9999, biomass_uncertainty_gCm2 = -9999))
        } # done_lat

    } else if (Cwood_stock_source == "McNicol") {

        print("Loading McNichol AGB maps...")

        # Create the full file paths estimates...uncertainty is assumed to be 250gC/m2
#        input_file = paste(path_to_Cwood,"mcnicol_AGC2007_0.25d.tif",sep="")
##        years_with_obs = c(2007)
#        input_file = append(input_file,paste(path_to_Cwood,"mcnicol_AGC2008_0.25d.tif",sep=""))
#        input_file = append(input_file,paste(path_to_Cwood,"mcnicol_AGC2009_0.25d.tif",sep=""))
#        input_file = append(input_file,paste(path_to_Cwood,"mcnicol_AGC2010_0.25d.tif",sep=""))
#        years_with_obs = c(2007:2010)
        input_file = paste(path_to_Cwood,"mcnicol_AGC2007_1km.tif",sep="")
#        years_with_obs = c(2007)
        input_file = append(input_file,paste(path_to_Cwood,"mcnicol_AGC2008_1km.tif",sep=""))
        input_file = append(input_file,paste(path_to_Cwood,"mcnicol_AGC2009_1km.tif",sep=""))
        input_file = append(input_file,paste(path_to_Cwood,"mcnicol_AGC2010_1km.tif",sep=""))
        years_with_obs = c(2007:2010)

        # Loop through each year and extract if appropriate
        done_lat = FALSE
        for (t in seq(1, length(years_with_obs))) {

             # determine whether the first year is within the analysis period
             if (years_with_obs[t] >= as.numeric(start) & years_with_obs[t] <= as.numeric(finish)) {

                 # Read in the estimate and uncertainty rasters
                 biomass = raster(input_file[t])

                 # Create raster with the target crs
                 target = raster(crs = ("+init=epsg:4326"), ext = extent(biomass), resolution = res(biomass))
                 # Check whether the target and actual analyses have the same CRS
                 if (compareCRS(biomass,target) == FALSE) {
                     # Resample to correct grid
                     biomass = resample(biomass, target, method="ngb") ; gc() ; removeTmpFiles()
                 }
                 # Extend the extent of the overall grid to the analysis domain
                 biomass = extend(biomass,cardamom_ext)
                 # Trim the extent of the overall grid to the analysis domain
                 biomass = crop(biomass,cardamom_ext)
                 # now remove the ones that are actual missing data
                 biomass[which(as.vector(biomass) < 0)] = NA
                 # Adjust spatial resolution of the datasets, this occurs in all cases
                 if (res(biomass)[1] != res(cardamom_ext)[1] | res(biomass)[2] != res(cardamom_ext)[2]) {
                     # Create raster with the target resolution
                     target = raster(crs = crs(cardamom_ext), ext = extent(cardamom_ext), resolution = res(cardamom_ext))
                     # Resample to correct grid
                     biomass = resample(biomass, target, method="bilinear") ; gc() ; removeTmpFiles()
                 } # Aggrgeate to resolution

                 # If the first file to be read extract the lat / long information
                 if (done_lat == FALSE) {
                     # Set flag to TRUE, impacts what will be returned from this function
                     done_lat = TRUE

                     # extract dimension information for the grid, note the axis switching between raster and actual array
                     xdim = dim(biomass)[2] ; ydim = dim(biomass)[1]
                     # extract the lat / long information needed
                     long = coordinates(biomass)[,1] ; lat = coordinates(biomass)[,2]
                     # restructure into correct orientation
                     long = array(long, dim=c(xdim,ydim))
                     lat = array(lat, dim=c(xdim,ydim))

                 } # extract lat / long...just the once

                 # break out from the rasters into arrays which we can manipulate
                 biomass = array(as.vector(unlist(biomass)), dim=c(xdim,ydim))

                 # Determine when in the analysis time series the observations should go
                 # NOTE: We assume the biomass estimate is placed at the beginning of the year

                 # What year of the analysis does the data fall?
                 obs_step = which(run_day_selector >= floor(which(analysis_years == years_with_obs[t]) * 365.25))[1]
                 obs_step = obs_step - (steps_per_year-1)

                 # Combine with the other time step
                 if (exists("place_obs_in_step")) {
                     # Output variables already exits to append them
                     place_obs_in_step = append(place_obs_in_step, obs_step)
                     biomass_gCm2 = append(biomass_gCm2, as.vector(biomass)) ; rm(biomass)
                 } else {
                     # Output variables do not already exist, assign them
                     place_obs_in_step = obs_step ; rm(obs_step)
                     biomass_gCm2 = as.vector(biomass) ; rm(biomass)
                 } # obs_step exists

             } # Is dataset within the analysis time period?

        } # looping available years

# Below ground estimatation following Saatchi et al (2011)
#        # Convert MgC/ha -> Mg/ha needed for Saatchi et al (2011)
#        biomass_gCm2 = biomass_gCm2 * 2.083333
#        # Use allometry to estimate below ground biomass stock and
#        # combined with the above ground (Mg/ha) to give a total woody biomass estimate
#        # Saatchi et al., (2011), PNAS, 108, 9899-9904, https://www.pnas.org/content/108/24/9899
#        biomass_gCm2 = biomass_gCm2 + (0.489 * biomass_gCm2 ** 0.89)
#        # Convert units of biomass and its uncertainty from MgCha -> gC/m2
#        biomass_gCm2 = biomass_gCm2 * 1e2 * 0.48
# Below round estimation following Ryan et al., (2011)
        # Ryan, C. M., M. Williams and J. Grace (2011).
        # "Above and Below Ground Carbon Stocks in a Miombo Woodland Landscape of Mozambique." Biotropica 43: 423-432
        biomass_gCm2 = biomass_gCm2 + (0.42 * biomass_gCm2)
        # Convert units of biomass and its uncertainty from MgCha -> gC/m2
        biomass_gCm2 = biomass_gCm2 * 1e2

        # Re-construct arrays for output
        idim = dim(lat)[1] ; jdim = dim(long)[2] ; tdim = length(biomass_gCm2) / (idim * jdim)
        biomass_gCm2 = array(biomass_gCm2, dim=c(idim,jdim,tdim))
        # Now create the uncertainty arrays, currently assuming constant 250 gC/m2 uncertainty
        # Ask Mat Williams
        biomass_uncertainty_gCm2 = array(250, dim=c(idim,jdim,tdim))

        # clean up variables
        gc(reset=TRUE,verbose=FALSE)

        if (done_lat) {
            # Output variables
            return(list(place_obs_in_step = place_obs_in_step, lat = lat, long = long,
                        biomass_gCm2 = biomass_gCm2, biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2))
        } else {
            # Output dummy variables
            return(list(place_obs_in_step = -9999, lat = -9999, long = -9999,
                        biomass_gCm2 = -9999, biomass_uncertainty_gCm2 = -9999))
        } # done_lat

    } else {
        # Output dummy variables
        return(list(place_obs_in_step = -9999, lat = -9999, long = -9999,
                    biomass_gCm2 = -9999, biomass_uncertainty_gCm2 = -9999))
    } # which biomass source?

} # function end load_biomass_stocks_maps_for_extraction

## Use byte compile
load_biomass_stocks_maps_for_extraction<-cmpfun(load_biomass_stocks_maps_for_extraction)
