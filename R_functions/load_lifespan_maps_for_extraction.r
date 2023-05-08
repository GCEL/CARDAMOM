
###
## Function to load leaf lifespan maps
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_lifespan_maps_for_extraction<-function(latlon_in,lifespan_source,cardamom_ext,spatial_type) {

    ###
    ## Select the correct LCA source for specific time points

    if (lifespan_source == "Tupek") {
	
		# let the user know this might take some time
        print("Loading Tupek lifespan map...")

        # Create the full file paths estimates and their uncertainty (gC/m2)
        input_file = list.files(path_to_lifespan)
        # extract only .tif files, $ symbol asks for strings that end in the given pattern
        # The \\ also specifies that the . is not to be considered a wildcard
        input_file = input_file[grepl("\\.tif$",input_file) == TRUE]
        # Extract the uncertainty files from the original list
        #unc_input_file = input_file[grepl("lifespan_SD",input_file) == TRUE]
        input_file = input_file[grepl("lifespan_SD",input_file) == FALSE]
        # Check that we have the same number of files for both lifespan and uncertainty
        #if (length(input_file) != length(unc_input_file)) {stop("Different number of observation and uncertainty files found...")}
        #if (length(input_file) > 1 | length(unc_input_file) > 1) {stop("More than one file has been found for the estimate and its uncertainty, there should only be one")}
		
		 # Read in the estimate and uncertainty rasters
        lifespan = raster(paste(path_to_lifespan,input_file,sep=""))
        #lifespan_uncertainty = raster(paste(path_to_lifespan,unc_input_file,sep=""))

        # Create raster with the target crs
        target = raster(crs = ("+init=epsg:4326"), ext = extent(lifespan), resolution = res(lifespan))
        # Check whether the target and actual analyses have the same CRS
        if (compareCRS(lifespan,target) == FALSE) {
            # Resample to correct grid
            lifespan = resample(lifespan, target, method="ngb") ; gc() ; removeTmpFiles()
            #lifespan_uncertainty_yrs = resample(lifespan_uncertainty_yrs, target, method="ngb") ; gc() ; removeTmpFiles()
        }
		# Extend the extent of the overall grid to the analysis domain
        lifespan = extend(lifespan,cardamom_ext) #; lifespan_uncertainty_yrs = extend(lifespan_uncertainty_yrs,cardamom_ext)
        # Trim the extent of the overall grid to the analysis domain
        lifespan = crop(lifespan,cardamom_ext) #; lifespan_uncertainty_yrs = crop(lifespan_uncertainty_yrs,cardamom_ext)
        # now remove the ones that are actual missing data
        lifespan[which(as.vector(lifespan) < 0)] = NA
        #lifespan_uncertainty_yrs[which(as.vector(lifespan_uncertainty_yrs) < 0)] = NA
        # If this is a gridded analysis and the desired CARDAMOM resolution is coarser than the currently provided then aggregate here
        # Despite creation of a cardamom_ext for a site run do not allow aggragation here as tis will damage the fine resolution datasets
        if (spatial_type == "grid") {
            if (res(lifespan)[1] < res(cardamom_ext)[1] | res(lifespan)[2] < res(cardamom_ext)[2]) {

                # Create raster with the target resolution
                target = raster(crs = crs(cardamom_ext), ext = extent(cardamom_ext), resolution = res(cardamom_ext))

                # Resample to correct grid
                lifespan = resample(lifespan, target, method="bilinear") ; gc() ; removeTmpFiles()
                #lifespan_uncertainty_yrs = resample(lifespan_uncertainty_yrs, target, method="bilinear") ; gc() ; removeTmpFiles()

            } # Aggrgeate to resolution
        } # spatial_type == "grid"
		
		# extract dimension information for the grid, note the axis switching between raster and actual array
        xdim = dim(lifespan)[2] ; ydim = dim(lifespan)[1]
        # extract the lat / long information needed
        long = coordinates(lifespan)[,1] ; lat = coordinates(lifespan)[,2]
        # restructure into correct orientation
        long = array(long, dim=c(xdim,ydim))
        lat = array(lat, dim=c(xdim,ydim))
        # break out from the rasters into arrays which we can manipulate
        lifespan = array(as.vector(unlist(lifespan)), dim=c(xdim,ydim))
        #lifespan_uncertainty_yrs = array(as.vector(unlist(lifespan_uncertainty_yrs)), dim=c(xdim,ydim))

        # Output variables
        return(list(lat = lat, long = long, lifespan = lifespan))#, lifespan_uncertainty_yrs = lifespan_uncertainty_yrs))

    } else {
	
		  # Output dummy variables
        return(list(lat = -9999, long = -9999, lifespan = -9999)) #, lifespan_uncertainty_yrs = -9999))

    } # which lifespan source?

} # function end load_lifespan_maps_for_extraction

## Use byte compile
load_lifespan_maps_for_extraction<-cmpfun(load_lifespan_maps_for_extraction)
