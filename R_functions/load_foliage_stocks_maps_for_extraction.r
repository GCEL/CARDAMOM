###
## Function to load biomass maps which apply to gridded domain
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_foliage_stocks_maps_for_extraction<-function(latlon_in,Cfol_stock_source,start,finish,timestep_days,cardamom_ext,spatial_type) {

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
	
	
	if (Cfol_stock_source == "MSNFI") {

        # this is a very bespoke modification so leave it here to avoid getting lost
        print("Loading Finland MSNFI Cfol")

        # Create the full file paths estimates and their uncertainty (MgC/ha)
        input_file = list.files(path_to_Cfol)
        # extract only .tif files, $ symbol asks for strings that end in the given pattern
        # The \\ also specifies that the . is not to be considered a wildcard
        input_file = input_file[grepl("\\.tif$",input_file) == TRUE]
        # Extract the specific files from the original list
        input_unc_file = input_file[grepl("Cfol_uncertainty_gCm2",input_file) == TRUE]
        # Extract the specific files from the original list
        input_file = input_file[grepl("Cfol_gCm2",input_file) == TRUE]

        # Sense check
        if (length(input_unc_file) != length(input_file)) {stop("number of uncertainty and data files differ for Cfol_stock_source = MSNFI")}

        # Determine the number of years found
        years_with_obs = gsub("Cfol_gCm2_","",input_file)
        years_with_obs = as.numeric(gsub("\\.tif$","",years_with_obs))

        # Loop through each year and extract if appropriate
        done_lat = FALSE
        for (t in seq(1, length(years_with_obs))) {

             # determine whether the first year is within the analysis period
             if (years_with_obs[t] >= as.numeric(start) & years_with_obs[t] <= as.numeric(finish)) {

                 # Read in the estimate and uncertainty rasters
                 foliage = rast(paste(path_to_Cfol,input_file[t],sep=""))
                 foliage_unc = rast(paste(path_to_Cfol,input_unc_file[t],sep=""))

                 # Create raster with the target crs
                 target = rast(crs = ("+init=epsg:4326"), extent = ext(foliage), resolution = res(foliage))
                 # Check whether the target and actual analyses have the same CRS
                 if (compareGeom(foliage,target) == FALSE) {
                     # Resample to correct grid
                     foliage = resample(foliage, target, method="near") ; gc()
                     foliage_unc = resample(foliage_unc, target, method="near") ; gc() 
                 }
                 # Extend the extent of the overall grid to the analysis domain
                 foliage = extend(foliage,cardamom_ext)
                 foliage_unc = extend(foliage_unc,cardamom_ext)
                 # Trim the extent of the overall grid to the analysis domain
                 foliage = crop(foliage,cardamom_ext)
                 foliage_unc = crop(foliage_unc,cardamom_ext)
                 # now remove the ones that are actual missing data
                 foliage[which(as.vector(foliage) < 0)] = NA
                 foliage_unc[which(as.vector(foliage_unc) < 0)] = NA

                 # Adjust spatial resolution of the datasets, this occurs in all cases
                 if (res(foliage)[1] != res(cardamom_ext)[1] | res(foliage)[2] != res(cardamom_ext)[2]) {

                     # Create raster with the target resolution
                     target = rast(crs = crs(cardamom_ext), extent = ext(cardamom_ext), resolution = res(cardamom_ext))
                     # Resample to correct grid
                     foliage = resample(foliage, target, method="bilinear") ; gc()
                     foliage_unc = resample(foliage_unc, target, method="bilinear") ; gc() 

                 } # Aggrgeate to resolution


                 # If the first file to be read extract the lat / long information
                 if (done_lat == FALSE) {
                     # Set flag to TRUE, impacts what will be returned from this function
                     done_lat = TRUE

                     # extract dimension information for the grid, note the axis switching between raster and actual array
                     xdim = dim(foliage)[2] ; ydim = dim(foliage)[1]
                     # extract the lat / long information needed
                     long = crds(foliage,df=TRUE, na.rm=FALSE)
                     lat  = long$y ; long = long$x
                     # restructure into correct orientation
                     long = array(long, dim=c(xdim,ydim))
                     lat = array(lat, dim=c(xdim,ydim))

                 } # extract lat / long...just the once

                 # break out from the rasters into arrays which we can manipulate
                 foliage = array(as.vector(unlist(foliage)), dim=c(xdim,ydim))
                 foliage_unc = array(as.vector(unlist(foliage_unc)), dim=c(xdim,ydim))

                 # Determine when in the analysis time series the observations should go
                 # NOTE: We assume the foliage estimate is placed at the beginning of the year

                 # What year of the analysis does the data fall?
                 obs_step = which(run_day_selector >= floor(which(analysis_years == years_with_obs[t]) * 365.25))[1]
                 obs_step = obs_step - (steps_per_year-1)
                 # Combine with the other time step
                 if (exists("place_obs_in_step")) {
                     # Output variables already exits to append them
                     place_obs_in_step = append(place_obs_in_step, obs_step)
                     foliage_gCm2 = append(foliage_gCm2, as.vector(foliage)) ; rm(foliage)
                     foliage_uncertainty_gCm2 = append(foliage_uncertainty_gCm2, as.vector(foliage_unc)) ; rm(foliage_unc)
                 } else {
                     # Output variables do not already exist, assign them
                     place_obs_in_step = obs_step ; rm(obs_step)
                     foliage_gCm2 = as.vector(foliage) ; rm(foliage)
                     foliage_uncertainty_gCm2 = as.vector(foliage_unc) ; rm(foliage_unc)
                 } # obs_step exists

             } # Is dataset within the analysis time period?

        } # looping available years


        # Re-construct arrays for output
        idim = dim(lat)[1] ; jdim = dim(long)[2] ; tdim = length(foliage_gCm2) / (idim * jdim)
        foliage_gCm2 = array(foliage_gCm2, dim=c(idim,jdim,tdim))
        foliage_uncertainty_gCm2 = array(foliage_uncertainty_gCm2, dim=c(idim,jdim,tdim))
        #foliage_uncertainty_gCm2 = foliage_gCm2 * 0.18 # assume uncertainty is 18 %
        if (done_lat) {
            # Output variables
            return(list(place_obs_in_step = place_obs_in_step, lat = lat, long = long,
                        foliage_gCm2 = foliage_gCm2, foliage_uncertainty_gCm2 = foliage_uncertainty_gCm2))
        } else {
            # Output dummy variables
            return(list(place_obs_in_step = -9999, lat = -9999, long = -9999,
                        foliage_gCm2 = -9999, foliage_uncertainty_gCm2 = -9999))
        } # done_lat


    } else {
        # Output dummy variables
        return(list(place_obs_in_step = -9999, lat = -9999, long = -9999,
                    foliage_gCm2 = -9999, foliage_uncertainty_gCm2 = -9999))
    } # which foliage source?

} # function end load_foliage_stocks_maps_for_extraction

## Use byte compile
load_foliage_stocks_maps_for_extraction<-cmpfun(load_foliage_stocks_maps_for_extraction)
