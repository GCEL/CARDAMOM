
###
## Function to load potential biomass maps to be used as attractor
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_potential_biomass_maps_for_extraction<-function(latlon_in,Cwood_potential_source,start,finish,timestep_days,cardamom_ext,spatial_type) {

   if (Cwood_potential_source == "UoE_potAGB") {

       # Read in desired potential AGB map
       infile = nc_open(paste(path_to_Cwood_potential,"BRA_020_AGB_potential_RFR_avitabile_worldclim_soilgrids_dist2edges_final.nc", sep=""))
       lat = ncvar_get(infile, "lat") ; long = ncvar_get(infile, "lon")
       # Estimates are in Mg/ha, convert to gC/m2 later
       biomass_gCm2 = ncvar_get(infile, "AGBpot")
       biomass_uncertainty_gCm2 = ncvar_get(infile, "AGBpot_max")
       tmp = ncvar_get(infile, "AGBpot_min")
       biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 - tmp ; rm(tmp)

       # Create lat / long grid from vectors
       idim = length(long) ; jdim = length(lat)
       lat = array(lat, dim=c(jdim,idim)) ; lat = t(lat)
       long = array(long, dim=c(idim,jdim))

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
       # Extend the extent of the overall grid to the analysis domain
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

       # Use allometry to estimate below ground biomass stock and
       # combined with the above ground (Mg/ha) to give a total woody biomass estimate
       # Saatchi et al., (2011), PNAS, 108, 9899-9904, https://www.pnas.org/content/108/24/9899
       biomass_gCm2 = biomass_gCm2 + (0.489 * biomass_gCm2 ** 0.89)
       biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 + (0.489 * biomass_uncertainty_gCm2 ** 0.89)
       # Convert units of biomass and its uncertainty from MgCha -> gC/m2
       biomass_gCm2 = biomass_gCm2 * 1e2 * 0.48
       biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 1e2 * 0.48

       # Output variables
       return(list(lat = lat, long = long,
                   biomass_gCm2 = biomass_gCm2, biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2))

   } else {

       # Output variables
       return(list(lat = -9999, long = -9999,
                   biomass_gCm2 = -9999, biomass_uncertainty_gCm2 = -9999))

   } # Cwood_potential_source


} # function end load_potential_biomass_maps_for_extraction

## Use byte compile
load_potential_biomass_maps_for_extraction<-cmpfun(load_potential_biomass_maps_for_extraction)
