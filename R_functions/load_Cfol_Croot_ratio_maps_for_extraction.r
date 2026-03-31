
###
## Function to load Cfol:Croot maps
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_Cfol_Croot_ratio_maps_for_extraction<-function(latlon_in,Cfol_Croot_ratio_source,cardamom_ext,spatial_type) {

    ###
    ## Select the correct LL source for specific time points

    if (Cfol_Croot_ratio_source == "Helmisaari") {
	
		# let the user know this might take some time
        print("Loading Cfol_Croot_ratio map...")

        # Create the full file paths estimates and their uncertainty (gC/m2)
        #input_file = list.files(path_to_max_root)
        # extract only .tif files, $ symbol asks for strings that end in the given pattern
        # The \\ also specifies that the . is not to be considered a wildcard
        #input_file = input_file[grepl("\\.tif$",input_file) == TRUE]
        # Extract the uncertainty files from the original list
        #unc_input_file = input_file[grepl("max_root_SD",input_file) == TRUE]
        #input_file = input_file[grepl("max_root_SD",input_file) == FALSE]
        # Check that we have the same number of files for both max_root and uncertainty
        #if (length(input_file) != length(unc_input_file)) {stop("Different number of observation and uncertainty files found...")}
        #if (length(input_file) > 1 | length(unc_input_file) > 1) {stop("More than one file has been found for the estimate and its uncertainty, there should only be one")}
		
		 # Read in the estimate and uncertainty rasters
        Cfol_Croot_ratio = rast(paste(path_to_Cfol_Croot_ratio,"Cfol_Croots_ratio.tif",sep=""))
        #lifespan_uncertainty = rast(paste(path_to_lifespan,unc_input_file,sep=""))

        # Create raster with the target crs
        target = rast(crs = ("+init=epsg:4326"), extent = ext(Cfol_Croot_ratio), resolution = res(Cfol_Croot_ratio))
        # Check whether the target and actual analyses have the same CRS
        if (compareGeom(Cfol_Croot_ratio,target) == FALSE) {
            # Resample to correct grid
            Cfol_Croot_ratio = resample(Cfol_Croot_ratio, target, method="bilinear") ; gc() 
            #lifespan_uncertainty_yrs = resample(lifespan_uncertainty_yrs, target, method="near") ; gc() 
        }
		# Extend the extent of the overall grid to the analysis domain
        Cfol_Croot_ratio = extend(Cfol_Croot_ratio,cardamom_ext) #; max_root_uncertainty_yrs = extend(max_root_uncertainty_yrs,cardamom_ext)
        # Trim the extent of the overall grid to the analysis domain
        Cfol_Croot_ratio = crop(Cfol_Croot_ratio,cardamom_ext) #; max_root_uncertainty_yrs = crop(max_root_uncertainty_yrs,cardamom_ext)
        # now remove the ones that are actual missing data
        Cfol_Croot_ratio[which(as.vector(Cfol_Croot_ratio) < 0)] = NA
        #max_root_uncertainty_yrs[which(as.vector(max_root_uncertainty_yrs) < 0)] = NA
        # If this is a gridded analysis and the desired CARDAMOM resolution is coarser than the currently provided then aggregate here
        # Despite creation of a cardamom_ext for a site run do not allow aggragation here as tis will damage the fine resolution datasets
        if (spatial_type == "grid") {
            if (res(Cfol_Croot_ratio)[1] < res(cardamom_ext)[1] | res(Cfol_Croot_ratio)[2] < res(cardamom_ext)[2]) {

                # Create raster with the target resolution
                target = rast(crs = crs(cardamom_ext), extent = ext(cardamom_ext), resolution = res(cardamom_ext))

                # Resample to correct grid
                Cfol_Croot_ratio = resample(Cfol_Croot_ratio, target, method="bilinear") ; gc() 
                #max_root_uncertainty_yrs = resample(max_root_uncertainty_yrs, target, method="bilinear") ; gc() 

            } # Aggrgeate to resolution
        } # spatial_type == "grid"
		
		# extract dimension information for the grid, note the axis switching between raster and actual array
        xdim = dim(Cfol_Croot_ratio)[2] ; ydim = dim(Cfol_Croot_ratio)[1]
        # extract the lat / long information needed
        long = crds(Cfol_Croot_ratio,df=T,na.rm=F) 
		lat = long$y ; long = long$x
        # restructure into correct orientation
        long = array(long, dim=c(xdim,ydim))
        lat = array(lat, dim=c(xdim,ydim))
        # break out from the rasters into arrays which we can manipulate
        Cfol_Croot_ratio = array(as.vector(unlist(Cfol_Croot_ratio)), dim=c(xdim,ydim))
        #max_root_uncertainty_yrs = array(as.vector(unlist(max_root_uncertainty_yrs)), dim=c(xdim,ydim))

        # Output variables
        return(list(lat = lat, long = long, Cfol_Croot_ratio = Cfol_Croot_ratio))#, max_root_uncertainty_yrs = max_root_uncertainty_yrs))

    } else {
	
		  # Output dummy variables
        return(list(lat = -9999, long = -9999, Cfol_Croot_ratio = -9999)) #, max_root_uncertainty_yrs = -9999))

    } # which max_root source?

} # function end load_max_root_maps_for_extraction

## Use byte compile
load_Cfol_Croot_ratio_maps_for_extraction<-cmpfun(load_Cfol_Croot_ratio_maps_for_extraction)
