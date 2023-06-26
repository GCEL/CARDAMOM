
###
## Function to load biomass data to be applied as an initial condition
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_initial_biomass_maps_for_extraction<-function(latlon_in,Cwood_initial_source,start,finish,timestep_days,cardamom_ext,spatial_type) {

      # Generate timing information need in most cases
      analysis_years = seq(as.numeric(start),as.numeric(finish))

      ###
      ## Select the correct Cwood source for initial conditions prior

      if (Cwood_initial_source == "mpi_biomass") {

          # let the user know this might take some time
          print("Loading MPI - >30N Forest Biomass map...")

          # Create the full file paths
          input_file = paste(path_to_Cwood_initial,"2014121116258biomass_v3_total.nc",sep="")

          # Open the first file
          data1 = nc_open(input_file)

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

          # Convert kgCm-2-> gCm-2
          biomass_gCm2 = biomass_gCm2 * 1e3
          biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 1e3

          # clean up variables
          gc(reset=TRUE,verbose=FALSE)

          # Output variables
          return(list(lat = lat, long = long,
                      biomass_gCm2 = biomass_gCm2, biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2))

      } else if (Cwood_initial_source == "Avitabile") {

          # let the user know this might take some time
          print("Loading processed Avitabile AGB...")

          # Create the full file paths
          input_file = paste(path_to_Cwood_initial,"Biomass_stocks_with_lat_long.nc",sep="")

          # Open the first file
          data1 = nc_open(input_file)

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

          # Convert to MgC/ha -> Mg/ha needed for Saatchi et al (2011)
          biomass_gCm2 = biomass_gCm2 * 2.083333
          biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 2.083333

          # Use allometry to estimate below ground biomass stock and
          # combined with the above ground (Mg/ha) to give a total woody biomass estimate
          # Saatchi et al., (2011), PNAS, 108, 9899-9904, https://www.pnas.org/content/108/24/9899
          biomass_gCm2 = biomass_gCm2 + (0.489 * biomass_gCm2 ** 0.89)
          biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 + (0.489 * biomass_uncertainty_gCm2 ** 0.89)
          # Convert units of biomass and its uncertainty from MgCha -> gC/m2
          biomass_gCm2 = biomass_gCm2 * 1e2 * 0.48
          biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 1e2 * 0.48

          # clean up variables
          gc(reset=TRUE,verbose=FALSE)

          # Output variables
          return(list(lat = lat, long = long,
                      biomass_gCm2 = biomass_gCm2, biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2))

    } else if (Cwood_initial_source == "UoL_stable_forest") {

        # this is a very bespoke modification so leave it here to avoid getting lost
        print("Loading UoL Stable Forests AGB...")

        # Read in estimate and uncertainty rasters
        # NOTE: that map is above ground, total biomass will be estimated later
        biomass_gCm2 = raster(paste(path_to_Cwood_initial,"Kenya_0.25deg_AGB_stable_forest_2015_2017.tif", sep=""))
        biomass_uncertainty_gCm2 = raster(paste(path_to_Cwood_initial,"Kenya_0.25deg_AGB_std_stable_forest_2015_2017.tif", sep=""))

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

        # Convert gC/m2 -> Mg/ha needed for Saatchi et al (2011)
        biomass_gCm2 = biomass_gCm2 * 2.083333 * 1e-2
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 2.083333 * 1e-2
        # Use allometry to estimate below ground biomass stock and
        # combined with the above ground (Mg/ha) to give a total woody biomass estimate
        # Saatchi et al., (2011), PNAS, 108, 9899-9904, https://www.pnas.org/content/108/24/9899
        biomass_gCm2 = biomass_gCm2 + (0.489 * biomass_gCm2 ** 0.89)
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 + (0.489 * biomass_uncertainty_gCm2 ** 0.89)
        # Now back to desired units gC/m2
        biomass_gC = biomass_gCm2 * 0.48 * 1e2
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 0.48 * 1e2

        # Output variables
        return(list(lat = lat, long = long,
                    biomass_gCm2 = biomass_gCm2, biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2))

    } else if (Cwood_initial_source == "UoL_stable_savannah") {

        # this is a very bespoke modification so leave it here to avoid getting lost
        print("Loading UoL Stable Savannah AGB...")

        # Read in estimate and uncertainty rasters
        # NOTE: that map is above ground, total biomass will be estimated later
        biomass_gCm2 = raster(paste(path_to_Cwood_initial,"Kenya_0.25deg_AGB_stable_savannah_2015_2017.tif", sep=""))
        biomass_uncertainty_gCm2 = raster(paste(path_to_Cwood_initial,"Kenya_0.25deg_AGB_std_stable_savannah_2015_2017.tif", sep=""))

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

        # Convert gC/m2 -> Mg/ha needed for Saatchi et al (2011)
        biomass_gCm2 = biomass_gCm2 * 2.083333 * 1e-2
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 2.083333 * 1e-2
        # Use allometry to estimate below ground biomass stock and
        # combined with the above ground (Mg/ha) to give a total woody biomass estimate
        # Saatchi et al., (2011), PNAS, 108, 9899-9904, https://www.pnas.org/content/108/24/9899
        biomass_gCm2 = biomass_gCm2 + (0.489 * biomass_gCm2 ** 0.89)
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 + (0.489 * biomass_uncertainty_gCm2 ** 0.89)
        # Now back to desired units gC/m2
        biomass_gC = biomass_gCm2 * 0.48 * 1e2
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 0.48 * 1e2

        # Output variables
        return(list(lat = lat, long = long,
                    biomass_gCm2 = biomass_gCm2, biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2))

    } else if (Cwood_initial_source == "Rainfor") {

        # this is a very bespoke modification so leave it here to avoid getting lost
        print("Loading UoL Rainfor AGB...")

        # Read in estimate and uncertainty rasters
        # NOTE: this data assimilates total biomass
        biomass_gCm2 = raster(paste(path_to_Cwood_initial,"wood_biomass_gCm2_2000.tif", sep=""))
        biomass_uncertainty_gCm2 = raster(paste(path_to_Cwood_initial,"unc_wood_biomass_gCm2_2000.tif", sep=""))

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

        # Output variables
        return(list(lat = lat, long = long,
                    biomass_gCm2 = biomass_gCm2, biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2))

    } else {
         # Output variables
         return(list(lat = -9999, long = -9999,
                     biomass_gCm2 = -9999, biomass_uncertainty_gCm2 = -9999))
    } # which biomass source?

} # function end load_initial_biomass_maps_for_extraction

## Use byte compile
load_initial_biomass_maps_for_extraction<-cmpfun(load_initial_biomass_maps_for_extraction)
