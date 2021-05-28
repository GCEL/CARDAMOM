
###
## Function to load biomass data to be applied as an initial condition
###

load_initial_biomass_maps_for_extraction<-function(latlon_in,Cwood_initial_source,start,finish,timestep_days) {

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

          # Remove any missing or un-realistic data points
          biomass_gCm2[which(as.vector(biomass_gCm2) < 0)] = NA
          biomass_uncertainty_gCm2[which(as.vector(biomass_uncertainty_gCm2) < 0)] = NA

          # Filter around target area
          max_lat = max(latlon_in[,1])+1.0 ; max_long=max(latlon_in[,2])+1.0
          min_lat = min(latlon_in[,1])-1.0 ; min_long=min(latlon_in[,2])-1.0
          keep_lat_min = min(which(lat[1,] > min_lat))
          keep_lat_max = max(which(lat[1,] < max_lat))
          keep_long_min = min(which(long[,1] > min_long))
          keep_long_max = max(which(long[,1] < max_long))
          # Remove data outside of target area
          biomass_gCm2 = biomass_gCm2[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]
          biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]
          lat = lat[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]
          long = long[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]

          # Convert kgCm-2-> gCm-2
          biomass_gCm2 = biomass_gCm2 * 1e3
          biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2 * 1e3

          # Re-construct arrays for output
          idim = dim(lat)[1] ; jdim = dim(long)[2] ; tdim = length(biomass_gCm2) / (idim * jdim)
          biomass_gCm2 = array(biomass_gCm2, dim=c(idim,jdim,tdim))
          biomass_uncertainty_gCm2 = array(biomass_uncertainty_gCm2, dim=c(idim,jdim,tdim))

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

          # Remove any missing or un-realistic data points
          biomass_gCm2[which(as.vector(biomass_gCm2) < 0)] = NA
          biomass_uncertainty_gCm2[which(as.vector(biomass_uncertainty_gCm2) < 0)] = NA

          # Filter around target area
          max_lat = max(latlon_in[,1])+1.0 ; max_long=max(latlon_in[,2])+1.0
          min_lat = min(latlon_in[,1])-1.0 ; min_long=min(latlon_in[,2])-1.0
          keep_lat_min = min(which(lat[1,] > min_lat))
          keep_lat_max = max(which(lat[1,] < max_lat))
          keep_long_min = min(which(long[,1] > min_long))
          keep_long_max = max(which(long[,1] < max_long))
          # Remove data outside of target area
          biomass_gCm2 = biomass_gCm2[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]
          biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]
          lat = lat[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]
          long = long[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]

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

          # Re-construct arrays for output
          idim = dim(lat)[1] ; jdim = dim(long)[2] ; tdim = length(biomass_gCm2) / (idim * jdim)
          biomass_gCm2 = array(biomass_gCm2, dim=c(idim,jdim,tdim))
          biomass_uncertainty_gCm2 = array(biomass_uncertainty_gCm2, dim=c(idim,jdim,tdim))

          # clean up variables
          gc(reset=TRUE,verbose=FALSE)

          # Output variables
          return(list(lat = lat, long = long,
                      biomass_gCm2 = biomass_gCm2, biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2))

      } else if (Cwood_initial_source == "GlobBIOMASS") {

         # let the user know this might take some time
         print("Loading processed GlobBIOMASS for subsequent sub-setting ...")

         # Create the full file paths to both 2010 and 2017 AGB estimates
         # Units MgC/ha conversion to gC/m2 will occur later
         input_file = paste(path_to_Cwood_initial,"AGBiomass_stocks_2010_with_lat_long.nc",sep="")
         input_file = append(input_file,paste(path_to_Cwood_initial,"AGBiomass_stocks_2017_with_lat_long.nc",sep=""))
         years_with_obs = c(2010,2017)

         # Determine which year is closest to the start point
         t = years_with_obs - as.numeric(start)
         t = which(t == min(t))[1]

         # Open the first file
         data1 = nc_open(input_file[t])

         # Read in lat / long information
         lat = ncvar_get(data1, "lat") ; long = ncvar_get(data1, "long")

         # Determine our target area
         max_lat = max(latlon_in[,1])+1.0 ; max_long = max(latlon_in[,2])+1.0
         min_lat = min(latlon_in[,1])-1.0 ; min_long = min(latlon_in[,2])-1.0
         keep_lat_min = min(which(lat[1,] > min_lat)) ; keep_lat_max = max(which(lat[1,] < max_lat))
         keep_long_min = min(which(long[,1] > min_long)) ; keep_long_max = max(which(long[,1] < max_long))
         lat = lat[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]
         long = long[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]

         # Read above ground biomass estimate and uncertainty
         biomass_gCm2 = ncvar_get(data1, "AGBiomass")
         biomass_uncertainty_gCm2 = ncvar_get(data1, "AGBiomass_Uncertainty")
         # Close files and tidy up
         nc_close(data1) ; gc(reset=TRUE,verbose=FALSE)

         # Remove data outside of target area
         biomass_gCm2 = biomass_gCm2[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]
         biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]
         # Remove missing data flags or un-realistic values
         biomass_gCm2[which(as.vector(biomass_gCm2) < 0)] = NA
         biomass_uncertainty_gCm2[which(as.vector(biomass_uncertainty_gCm2) < 0)] = NA

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

         # Re-construct arrays for output
         idim = dim(lat)[1] ; jdim = dim(long)[2] ; tdim = length(biomass_gCm2) / (idim * jdim)
         biomass_gCm2 = array(biomass_gCm2, dim=c(idim,jdim,tdim))
         biomass_uncertainty_gCm2 = array(biomass_uncertainty_gCm2, dim=c(idim,jdim,tdim))

         # clean up variables
         gc(reset=TRUE,verbose=FALSE)

         # Output variables
         return(list(lat = lat, long = long,
                     biomass_gCm2 = biomass_gCm2, biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2))

    } else if (Cwood_initial_source == "ESA-CCI") {

       # let the user know this might take some time
       print("Loading processed GlobBIOMASS for subsequent sub-setting ...")

       # Create the full file paths to both 2010 and 2017 AGB estimates
       # Units MgC/ha conversion to gC/m2 will occur later
       input_file = paste(path_to_Cwood_initial,"AGBiomass_stocks_2010_with_lat_long(.)*.nc",sep="")
       input_file = append(input_file,paste(path_to_Cwood_initial,"AGBiomass_stocks_2017_with_lat_long(.)*.nc",sep=""))
       years_with_obs = c(2017,2018)

       # Determine which year is closest to the start point
       t = years_with_obs - as.numeric(start)
       t = which(t == min(t))[1]

       # Open the first file
       data1 = nc_open(input_file[t])

       # Read in lat / long information
       lat = ncvar_get(data1, "latitude") ; long = ncvar_get(data1, "longitude")

       # Determine our target area
       max_lat = max(latlon_in[,1])+1.0 ; max_long = max(latlon_in[,2])+1.0
       min_lat = min(latlon_in[,1])-1.0 ; min_long = min(latlon_in[,2])-1.0
       keep_lat_min = min(which(lat[1,] > min_lat)) ; keep_lat_max = max(which(lat[1,] < max_lat))
       keep_long_min = min(which(long[,1] > min_long)) ; keep_long_max = max(which(long[,1] < max_long))
       lat = lat[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]
       long = long[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]

       # Read above ground biomass estimate and uncertainty
       biomass_gCm2 = ncvar_get(data1, "AGBiomass")
       biomass_uncertainty_gCm2 = ncvar_get(data1, "AGBiomass_Uncertainty")
       # Close files and tidy up
       nc_close(data1) ; gc(reset=TRUE,verbose=FALSE)

       # Remove data outside of target area
       biomass_gCm2 = biomass_gCm2[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]
       biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]
       # Remove missing data flags or un-realistic values
       biomass_gCm2[which(as.vector(biomass_gCm2) < 0)] = NA
       biomass_uncertainty_gCm2[which(as.vector(biomass_uncertainty_gCm2) < 0)] = NA

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

       # Re-construct arrays for output
       idim = dim(lat)[1] ; jdim = dim(long)[2] ; tdim = length(biomass_gCm2) / (idim * jdim)
       biomass_gCm2 = array(biomass_gCm2, dim=c(idim,jdim,tdim))
       biomass_uncertainty_gCm2 = array(biomass_uncertainty_gCm2, dim=c(idim,jdim,tdim))

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

        # Store dimension information
        dims = dim(biomass_gCm2)[1:2]
        # Extract latitude / longitude information
        lat = coordinates(biomass_gCm2)
        # Split between long and lat
        long = lat[,1] ; lat = lat[,2]
        # Reconstruct the full lat / long grid and flip dimensions as needed
        long = array(long, dim=c(dims[2],dims[1]))
        lat = array(lat, dim=c(dims[2],dims[1]))
        long = long[,dim(long)[2]:1]
        lat = lat[,dim(lat)[2]:1]

        # filter around target area
        max_lat = max(latlon_in[,1])+1.0 ; max_long=max(latlon_in[,2])+1.0
        min_lat = min(latlon_in[,1])-1.0 ; min_long=min(latlon_in[,2])-1.0
        keep_lat_min = min(which(lat[1,] > min_lat))
        keep_lat_max = max(which(lat[1,] < max_lat))
        keep_long_min = min(which(long[,1] > min_long))
        keep_long_max = max(which(long[,1] < max_long))
        lat = lat[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]
        long = long[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]

        # Similarly break apart the raster and re-construct into the correct orientation
        biomass_gCm2 = array(as.vector(biomass_gCm2), dim=c(dims[2],dims[1]))
        biomass_uncertainty_gCm2 = array(as.vector(biomass_uncertainty_gCm2), dim=c(dims[2],dims[1]))
        biomass_gCm2 = biomass_gCm2[,dim(biomass_gCm2)[2]:1]
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2[,dim(biomass_uncertainty_gCm2)[2]:1]

        # remove data outside of target area
        biomass_gCm2 = biomass[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]
        biomass_uncertainty_gCm2 = biomass_uncertainty[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]

        # now remove the ones that are actual missing data
        biomass_gCm2[which(as.vector(biomass_gCm2) < 0)] = NA
        biomass_uncertainty_gCm2[which(as.vector(biomass_uncertainty_gCm2) < 0)] = NA

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

        # Store dimension information
        dims = dim(biomass_gCm2)[1:2]
        # Extract latitude / longitude information
        lat = coordinates(biomass_gCm2)
        # Split between long and lat
        long = lat[,1] ; lat = lat[,2]
        # Reconstruct the full lat / long grid and flip dimensions as needed
        long = array(long, dim=c(dims[2],dims[1]))
        lat = array(lat, dim=c(dims[2],dims[1]))
        long = long[,dim(long)[2]:1]
        lat = lat[,dim(lat)[2]:1]

        # filter around target area
        max_lat = max(latlon_in[,1])+1.0 ; max_long=max(latlon_in[,2])+1.0
        min_lat = min(latlon_in[,1])-1.0 ; min_long=min(latlon_in[,2])-1.0
        keep_lat_min = min(which(lat[1,] > min_lat))
        keep_lat_max = max(which(lat[1,] < max_lat))
        keep_long_min = min(which(long[,1] > min_long))
        keep_long_max = max(which(long[,1] < max_long))
        lat = lat[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]
        long = long[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]

        # Similarly break apart the raster and re-construct into the correct orientation
        biomass_gCm2 = array(as.vector(biomass_gCm2), dim=c(dims[2],dims[1]))
        biomass_uncertainty_gCm2 = array(as.vector(biomass_uncertainty_gCm2), dim=c(dims[2],dims[1]))
        biomass_gCm2 = biomass_gCm2[,dim(biomass_gCm2)[2]:1]
        biomass_uncertainty_gCm2 = biomass_uncertainty_gCm2[,dim(biomass_uncertainty_gCm2)[2]:1]

        # remove data outside of target area
        biomass_gCm2 = biomass[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]
        biomass_uncertainty_gCm2 = biomass_uncertainty[keep_long_min:keep_long_max,keep_lat_min:keep_lat_max]

        # now remove the ones that are actual missing data
        biomass_gCm2[which(as.vector(biomass_gCm2) < 0)] = NA
        biomass_uncertainty_gCm2[which(as.vector(biomass_uncertainty_gCm2) < 0)] = NA

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

    } else {
         # Output variables
         return(list(lat = -9999, long = -9999,
                     biomass_gCm2 = -9999, biomass_uncertainty_gCm2 = -9999))
    } # which biomass source?

} # function end
