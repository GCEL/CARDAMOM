
###
## Function to load wood mortality maps which apply to gridded domain
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).
# Created: 28/12/2021
# Last edited: 28/12/2021 (T. L. Smallman)

load_wood_mortality_maps_for_extraction<-function(Cwood_mortality_source,cardamom_ext,spatial_type,latlon_in,start,finish,timestep_days) {

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
    ## Select the correct Cwood increment source for specific time points

    if (exists("Cwood_mortality_source") == FALSE) {
        # Output dummy variables
        return(list(place_obs_in_step = -9999, lat = -9999, long = -9999,
                    Cwood_mortality_gCm2 = -9999, Cwood_mortality_uncertainty_gCm2 = -9999, Cwood_mortality_lag = -9999))
    }

    # If we got here then this must exist
    if (Cwood_mortality_source == "Rainfor") {

        # this is a very bespoke modification so leave it here to avoid getting lost
        print("Loading Rainfor wood mortality maps maps")

        # Create the full file paths estimates and their uncertainty (MgC/ha)
        input_file = list.files(path_to_Cwood_mortality)
        # extract only .tif files, $ symbol asks for strings that end in the given pattern
        # The \\ also specifies that the . is not to be considered a wildcard
        input_file = input_file[grepl("\\.tif$",input_file) == TRUE]
        # Extract the uncertainty files from the original list
        unc_input_file = input_file[grepl("unc_wood_mortality_gCm2",input_file) == TRUE]
        input_file = input_file[grepl("unc_wood_mortality_gCm2",input_file) == FALSE]
        # Check that we have the same number of files for both biomass and uncertainty
        if (length(input_file) != length(unc_input_file)) {stop("Different number of observation and uncertainty files found...")}

        # Extract timing information
        years_with_obs = gsub("wood_mortality_gCm2_","",input_file)
        years_with_obs = gsub("\\.tif$","",years_with_obs)
        # These files have two years provided, defining the beginning and ending of the averaging periods
        endings = rep(NA,length(years_with_obs)) ; avg_period_of_obs = rep(NA,length(years_with_obs))
        for (t in seq(1, length(years_with_obs))) {
             tmp = unlist(strsplit(years_with_obs[t],"_"))
             endings[t] = as.numeric(tmp[2])
             # Use the beginnings and endings to etermine the averaging period
             avg_period_of_obs[t] = length(c(as.numeric(tmp[1]):endings[t]))
        }
        # Overwrite years with obs with the endings, as these are the points which the information will be loaded
        years_with_obs = endings
        # Correct averaging period to that of the model timestep
        avg_period_of_obs = avg_period_of_obs * steps_per_year

        # Loop through each year and extract if appropriate
        done_lat = FALSE
        for (t in seq(1, length(years_with_obs))) {

             # determine whether the first year is within the analysis period
             if (years_with_obs[t] >= as.numeric(start) & years_with_obs[t] <= as.numeric(finish)) {

                 # Read in the estimate and uncertainty rasters
                 Cwood_mortality = raster(paste(path_to_Cwood_mortality,input_file[t],sep=""))
                 Cwood_mortality_uncertainty = raster(paste(path_to_Cwood_mortality,unc_input_file[t],sep=""))

                 # Create raster with the target crs
                 target = raster(crs = ("+init=epsg:4326"), ext = extent(Cwood_mortality), resolution = res(Cwood_mortality))
                 # Check whether the target and actual analyses have the same CRS
                 if (compareCRS(Cwood_mortality,target) == FALSE) {
                     # Resample to correct grid
                     Cwood_mortality = resample(Cwood_mortality, target, method="ngb") ; gc() ; removeTmpFiles()
                     Cwood_mortality_uncertainty = resample(Cwood_mortality_uncertainty, target, method="ngb") ; gc() ; removeTmpFiles()
                 }
                 # Extend the extent of the overall grid to the analysis domain
                 Cwood_mortality = extend(Cwood_mortality,cardamom_ext) ; Cwood_mortality_uncertainty = extend(Cwood_mortality_uncertainty,cardamom_ext)
                 # Trim the extent of the overall grid to the analysis domain
                 Cwood_mortality = crop(Cwood_mortality,cardamom_ext) ; Cwood_mortality_uncertainty = crop(Cwood_mortality_uncertainty,cardamom_ext)
                 # now remove the ones that are actual missing data
                 Cwood_mortality[which(as.vector(Cwood_mortality) < 0)] = NA
                 Cwood_mortality_uncertainty[which(as.vector(Cwood_mortality_uncertainty) < 0)] = NA

                 # Adjust spatial resolution of the datasets, this occurs in all cases
                 if (res(Cwood_mortality)[1] != res(cardamom_ext)[1] | res(Cwood_mortality)[2] != res(cardamom_ext)[2]) {

                     # Create raster with the target resolution
                     target = raster(crs = crs(cardamom_ext), ext = extent(cardamom_ext), resolution = res(cardamom_ext))

                     # Resample to correct grid
                     Cwood_mortality = resample(Cwood_mortality, target, method="bilinear") ; gc() ; removeTmpFiles()
                     Cwood_mortality_uncertainty = resample(Cwood_mortality_uncertainty, target, method="bilinear") ; gc() ; removeTmpFiles()

                 } # Aggrgeate to resolution

                 # If the first file to be read extract the lat / long information
                 if (done_lat == FALSE) {
                     # Set flag to TRUE, impacts what will be returned from this function
                     done_lat = TRUE

                     # extract dimension information for the grid, note the axis switching between raster and actual array
                     xdim = dim(Cwood_mortality)[2] ; ydim = dim(Cwood_mortality)[1]
                     # extract the lat / long information needed
                     long = coordinates(Cwood_mortality)[,1] ; lat = coordinates(Cwood_mortality)[,2]
                     # restructure into correct orientation
                     long = array(long, dim=c(xdim,ydim))
                     lat = array(lat, dim=c(xdim,ydim))

                 } # extract lat / long...just the once

                 # break out from the rasters into arrays which we can manipulate
                 Cwood_mortality = array(as.vector(unlist(Cwood_mortality)), dim=c(xdim,ydim))
                 Cwood_mortality_uncertainty = array(as.vector(unlist(Cwood_mortality_uncertainty)), dim=c(xdim,ydim))

                 # Determine when in the analysis time series the observations should go
                 # NOTE: We assume the biomass estimate is placed at the beginning of the year

                 # Put the observation at the end of the desired year
                 obs_step = which(run_day_selector >= floor(which(analysis_years == years_with_obs[t]) * 365.25))[1]
                 # Combine with the other time step
                 if (exists("place_obs_in_step")) {
                     # Output variables already exits to append them
                     Cwood_mortality_lag_step = append(Cwood_mortality_lag_step,rep(avg_period_of_obs[t], length(as.vector(Cwood_mortality))))
                     place_obs_in_step = append(place_obs_in_step, obs_step)
                     Cwood_mortality_gCm2 = append(Cwood_mortality_gCm2, as.vector(Cwood_mortality)) ; rm(Cwood_mortality)
                     Cwood_mortality_uncertainty_gCm2 = append(Cwood_mortality_uncertainty_gCm2, as.vector(Cwood_mortality_uncertainty)) ; rm(Cwood_mortality_uncertainty)
                 } else {
                     # Output variables do not already exist, assign them
                     Cwood_mortality_lag_step = rep(avg_period_of_obs[t], length(as.vector(Cwood_mortality)))
                     place_obs_in_step = obs_step ; rm(obs_step)
                     Cwood_mortality_gCm2 = as.vector(Cwood_mortality) ; rm(Cwood_mortality)
                     Cwood_mortality_uncertainty_gCm2 = as.vector(Cwood_mortality_uncertainty) ; rm(Cwood_mortality_uncertainty)
                 } # obs_step exists

             } # Is dataset within the analysis time period?

        } # looping available years

#        # Convert gC/m2 -> Mg/ha needed for Saatchi et al (2011)
#        Cwood_mortality_gCm2 = Cwood_mortality_gCm2 * 2.083333 * 1e2
#        Cwood_mortality_uncertainty_gCm2 = Cwood_mortality_uncertainty_gCm2 * 2.083333 * 1e2
#        # Use allometry to estimate below ground biomass stock and
#        # combined with the above ground (Mg/ha) to give a total woody biomass estimate
#        # Saatchi et al., (2011), PNAS, 108, 9899-9904, https://www.pnas.org/content/108/24/9899
#        Cwood_mortality_gCm2 = Cwood_mortality_gCm2 + (0.489 * Cwood_mortality_gCm2 ** 0.89)
#        Cwood_mortality_uncertainty_gCm2 = Cwood_mortality_uncertainty_gCm2 + (0.489 * Cwood_mortality_uncertainty_gCm2 ** 0.89)
#        # Now back to desired units gC/m2
#        Cwood_mortality_gCm2 = Cwood_mortality_gCm2 * 0.48 * 1e2
#        Cwood_mortality_uncertainty_gCm2 = Cwood_mortality_uncertainty_gCm2 * 0.48 * 1e2

        # Re-construct arrays for output
        idim = dim(lat)[1] ; jdim = dim(long)[2] ; tdim = length(Cwood_mortality_gCm2) / (idim * jdim)
        Cwood_mortality_gCm2 = array(Cwood_mortality_gCm2, dim=c(idim,jdim,tdim))
        Cwood_mortality_uncertainty_gCm2 = array(Cwood_mortality_uncertainty_gCm2, dim=c(idim,jdim,tdim))
        Cwood_mortality_lag_step = array(Cwood_mortality_lag_step, dim=c(idim,jdim,tdim))

        if (done_lat) {
            # Output variables
            return(list(place_obs_in_step = place_obs_in_step, lat = lat, long = long,
                        Cwood_mortality_gCm2 = Cwood_mortality_gCm2, Cwood_mortality_uncertainty_gCm2 = Cwood_mortality_uncertainty_gCm2,
                        Cwood_mortality_lag_step = Cwood_mortality_lag_step))
        } else {
            # Output dummy variables
            return(list(place_obs_in_step = -9999, lat = -9999, long = -9999,
                        Cwood_mortality_gCm2 = -9999, Cwood_mortality_uncertainty_gCm2 = -9999, Cwood_mortality_lag = -9999))
        } # done_lat

    } else {

       # Output dummy variables
       return(list(place_obs_in_step = -9999, lat = -9999, long = -9999,
                   Cwood_mortality_gCm2 = -9999, Cwood_mortality_uncertainty_gCm2 = -9999, Cwood_mortality_lag = -9999))

    } # which Cwood increment source?

} # function end load_biomass_stocks_maps_for_extraction

## Use byte compile
load_wood_mortality_maps_for_extraction<-cmpfun(load_wood_mortality_maps_for_extraction)
