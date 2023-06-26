
###
## Function to load leaf carbon per unit leaf area (gC/m2) maps
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_lca_maps_for_extraction<-function(latlon_in,lca_source,cardamom_ext,spatial_type) {

    ###
    ## Select the correct LCA source for specific time points

    if (lca_source == "Butler") {

        # let the user know this might take some time
        print("Loading Butler LCA map...")

        # Create the full file paths estimates and their uncertainty (gC/m2)
        input_file = list.files(path_to_lca)
        # extract only .tif files, $ symbol asks for strings that end in the given pattern
        # The \\ also specifies that the . is not to be considered a wildcard
        input_file = input_file[grepl("\\.tif$",input_file) == TRUE]
        # Extract the uncertainty files from the original list
        unc_input_file = input_file[grepl("LCA_SD",input_file) == TRUE]
        input_file = input_file[grepl("LCA_SD",input_file) == FALSE]
        # Check that we have the same number of files for both biomass and uncertainty
        if (length(input_file) != length(unc_input_file)) {stop("Different number of observation and uncertainty files found...")}
        if (length(input_file) > 1 | length(unc_input_file) > 1) {stop("More than one file has been found for the estimate and its uncertainty, there should only be one")}

        # Read in the estimate and uncertainty rasters
        lca_gCm2 = raster(paste(path_to_lca,input_file,sep=""))
        lca_uncertainty_gCm2 = raster(paste(path_to_lca,unc_input_file,sep=""))

        # Create raster with the target crs
        target = raster(crs = ("+init=epsg:4326"), ext = extent(lca_gCm2), resolution = res(lca_gCm2))
        # Check whether the target and actual analyses have the same CRS
        if (compareCRS(lca_gCm2,target) == FALSE) {
            # Resample to correct grid
            lca_gCm2 = resample(lca_gCm2, target, method="ngb") ; gc() ; removeTmpFiles()
            lca_uncertainty_gCm2 = resample(lca_uncertainty_gCm2, target, method="ngb") ; gc() ; removeTmpFiles()
        }
        # Extend the extent of the overall grid to the analysis domain
        lca_gCm2 = extend(lca_gCm2,cardamom_ext) ; lca_uncertainty_gCm2 = extend(lca_uncertainty_gCm2,cardamom_ext)
        # Trim the extent of the overall grid to the analysis domain
        lca_gCm2 = crop(lca_gCm2,cardamom_ext) ; lca_uncertainty_gCm2 = crop(lca_uncertainty_gCm2,cardamom_ext)
        # now remove the ones that are actual missing data
        lca_gCm2[which(as.vector(lca_gCm2) < 0)] = NA
        lca_uncertainty_gCm2[which(as.vector(lca_uncertainty_gCm2) < 0)] = NA
        # Adjust spatial resolution of the datasets, this occurs in all cases
        if (res(lca_gCm2)[1] != res(cardamom_ext)[1] | res(lca_gCm2)[2] != res(cardamom_ext)[2]) {

            # Create raster with the target resolution
            target = raster(crs = crs(cardamom_ext), ext = extent(cardamom_ext), resolution = res(cardamom_ext))

            # Resample to correct grid
            lca_gCm2 = resample(lca_gCm2, target, method="bilinear") ; gc() ; removeTmpFiles()
            lca_uncertainty_gCm2 = resample(lca_uncertainty_gCm2, target, method="bilinear") ; gc() ; removeTmpFiles()

        } # Aggrgeate to resolution

        # extract dimension information for the grid, note the axis switching between raster and actual array
        xdim = dim(lca_gCm2)[2] ; ydim = dim(lca_gCm2)[1]
        # extract the lat / long information needed
        long = coordinates(lca_gCm2)[,1] ; lat = coordinates(lca_gCm2)[,2]
        # restructure into correct orientation
        long = array(long, dim=c(xdim,ydim))
        lat = array(lat, dim=c(xdim,ydim))
        # break out from the rasters into arrays which we can manipulate
        lca_gCm2 = array(as.vector(unlist(lca_gCm2)), dim=c(xdim,ydim))
        lca_uncertainty_gCm2 = array(as.vector(unlist(lca_uncertainty_gCm2)), dim=c(xdim,ydim))

        # Output variables
        return(list(lat = lat, long = long, lca_gCm2 = lca_gCm2, lca_uncertainty_gCm2 = lca_uncertainty_gCm2))

    } else {

        # Output dummy variables
        return(list(lat = -9999, long = -9999, lca_gCm2 = -9999, lca_uncertainty_gCm2 = -9999))

    } # which LCA source?

} # function end load_lca_maps_for_extraction

## Use byte compile
load_lca_maps_for_extraction<-cmpfun(load_lca_maps_for_extraction)
