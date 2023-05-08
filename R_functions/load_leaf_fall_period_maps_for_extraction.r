
###
## Function to load leaf fall period map
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_leaf_fall_period_maps_for_extraction<-function(latlon_in,leaf_fall_period_source,cardamom_ext,spatial_type) {

    ###
    ## Select the correct LCA source for specific time points

    if (leaf_fall_period_source == "TG") {
	
		# let the user know this might take some time
        print("Loading leaf fall period map...")

        # Create the full file paths estimates and their uncertainty (gC/m2)
        input_file = list.files(path_to_leaf_fall_period)
        # extract only .tif files, $ symbol asks for strings that end in the given pattern
        # The \\ also specifies that the . is not to be considered a wildcard
        input_file = input_file[grepl("\\.tif$",input_file) == TRUE]
        # Extract the uncertainty files from the original list
        #unc_input_file = input_file[grepl("leaf_fall_period_SD",input_file) == TRUE]
        input_file = input_file[grepl("leaf_fall_period_SD",input_file) == FALSE]
        # Check that we have the same number of files for both leaf fall period and uncertainty
        #if (length(input_file) != length(unc_input_file)) {stop("Different number of observation and uncertainty files found...")}
        #if (length(input_file) > 1 | length(unc_input_file) > 1) {stop("More than one file has been found for the estimate and its uncertainty, there should only be one")}
		
		 # Read in the estimate and uncertainty rasters
        leaf_fall_period = raster(paste(path_to_leaf_fall_period,input_file,sep=""))
        #leaf_fall_period_uncertainty = raster(paste(path_to_leaf_fall_period,unc_input_file,sep=""))

        # Create raster with the target crs
        target = raster(crs = ("+init=epsg:4326"), ext = extent(leaf_fall_period), resolution = res(leaf_fall_period))
        # Check whether the target and actual analyses have the same CRS
        if (compareCRS(leaf_fall_period,target) == FALSE) {
            # Resample to correct grid
            leaf_fall_period = resample(leaf_fall_period, target, method="ngb") ; gc() ; removeTmpFiles()
            #leaf_fall_period_uncertainty_days = resample(leaf_fall_period_uncertainty_days, target, method="ngb") ; gc() ; removeTmpFiles()
        }
		# Extend the extent of the overall grid to the analysis domain
        leaf_fall_period = extend(leaf_fall_period,cardamom_ext) #; leaf_fall_period_uncertainty_days = extend(leaf_fall_period_uncertainty_days,cardamom_ext)
        # Trim the extent of the overall grid to the analysis domain
        leaf_fall_period = crop(leaf_fall_period,cardamom_ext) #; leaf_fall_period_uncertainty_days = crop(leaf_fall_period_uncertainty_days,cardamom_ext)
        # now remove the ones that are actual missing data
        leaf_fall_period[which(as.vector(leaf_fall_period) < 0)] = NA
        #leaf_fall_period_uncertainty_days[which(as.vector(leaf_fall_period_uncertainty_days) < 0)] = NA
        # If this is a gridded analysis and the desired CARDAMOM resolution is coarser than the currently provided then aggregate here
        # Despite creation of a cardamom_ext for a site run do not allow aggragation here as tis will damage the fine resolution datasets
        if (spatial_type == "grid") {
            if (res(leaf_fall_period)[1] < res(cardamom_ext)[1] | res(leaf_fall_period)[2] < res(cardamom_ext)[2]) {

                # Create raster with the target resolution
                target = raster(crs = crs(cardamom_ext), ext = extent(cardamom_ext), resolution = res(cardamom_ext))

                # Resample to correct grid
                leaf_fall_period = resample(leaf_fall_period, target, method="bilinear") ; gc() ; removeTmpFiles()
                #leaf_fall_period_uncertainty_days = resample(leaf_fall_period_uncertainty_days, target, method="bilinear") ; gc() ; removeTmpFiles()

            } # Aggrgeate to resolution
        } # spatial_type == "grid"
		
		# extract dimension information for the grid, note the axis switching between raster and actual array
        xdim = dim(leaf_fall_period)[2] ; ydim = dim(leaf_fall_period)[1]
        # extract the lat / long information needed
        long = coordinates(leaf_fall_period)[,1] ; lat = coordinates(leaf_fall_period)[,2]
        # restructure into correct orientation
        long = array(long, dim=c(xdim,ydim))
        lat = array(lat, dim=c(xdim,ydim))
        # break out from the rasters into arrays which we can manipulate
        leaf_fall_period = array(as.vector(unlist(leaf_fall_period)), dim=c(xdim,ydim))
        #leaf_fall_period_uncertainty_days = array(as.vector(unlist(leaf_fall_period_uncertainty_days)), dim=c(xdim,ydim))

        # Output variables
        return(list(lat = lat, long = long, leaf_fall_period = leaf_fall_period))#, leaf_fall_period_uncertainty_days = leaf_fall_period_uncertainty_days))

    } else {
	
		  # Output dummy variables
        return(list(lat = -9999, long = -9999, leaf_fall_period = -9999)) #, leaf_fall_period_uncertainty_days = -9999))

    } # which leaf fall period source?

} # function end load_leaf_fall_period_maps_for_extraction

## Use byte compile
load_leaf_fall_period_maps_for_extraction <-cmpfun(load_leaf_fall_period_maps_for_extraction)
