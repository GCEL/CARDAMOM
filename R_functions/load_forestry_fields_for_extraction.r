
###
## Function to load forestry planting and disturbance information from a gridded datset
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_forestry_fields_for_extraction<-function(latlon_in,forestry_source,years_to_load,cardamom_ext,spatial_type) {

  if (forestry_source == "GFW") {

      # let the user know this might take some time
      print("Loading processed Global Forest Watch clearance information for subsequent sub-setting ...")

      # Set flag
      lat_done = FALSE
      # Loop through years
      for (yrr in seq(1,length(years_to_load))){

          # Create file name # two possible versions of  naming syntax
          input_file_2a = paste(path_to_forestry,"GFW_forest_loss_",years_to_load[yrr],".nc",sep="")
          input_file_2b = Sys.glob(paste(path_to_forestry,"/tree_cover_loss_fraction_",years_to_load[yrr],"*.nc",sep=""))

          # Check it exists
          if (any(file.exists(input_file_2a),file.exists(input_file_2b))) {

              # Open the current file
              if (file.exists(input_file_2a)){
                  input_file_2 = input_file_2a
                  loss_variable = "forest_loss"
              } else{
                  input_file_2 = input_file_2b
                  loss_variable = "tree_cover_loss"
              }

              # Open the current file
              data2 = nc_open(input_file_2)
              # begin reading the forest loss information instead
              lat_in = ncvar_get(data2, "latitude") ; long_in = ncvar_get(data2, "longitude")
              # read year of forest loss informatin
              var1 = ncvar_get(data2, loss_variable)
              # tidy up
              nc_close(data2)

              # Check that lat / long are 2D arrays, if not we can try to force them
              if (length(dim(lat_in)) == 1 | length(dim(long_in)) == 1) {
                  # Ok, one (or hopefully both) of the lat / long variables is not in the expected grid format
                  if (length(dim(lat_in)) == 1 & length(dim(long_in)) == 2) {
                      # Not good
                      stop("GFW latitude is a vector but longitude is an array")
                  }
                  if (length(dim(lat_in)) == 2 & length(dim(long_in)) == 1) {
                      # Not good
                      stop("GFW longitude is a vector but latitude is an array")
                  }
                  # From here assume they are both vectors
                  # Turn lat_in / long_in from vectors to arrays
                  lat_in = t(array(lat_in, dim=c(dim(var1)[2],dim(var1)[1])))
                  long_in = array(long_in, dim=c(dim(var1)[1],dim(var1)[2]))
              }

              # Convert to a raster, assuming standad WGS84 grid
              var1 = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(var1))
              var1 = rast(var1, crs = ("+init=epsg:4326"), type="xyz")

              # Create raster with the target crs (technically this bit is not required)
              target = rast(crs = ("+init=epsg:4326"), ext = ext(var1), resolution = res(var1))
              # Check whether the target and actual analyses have the same CRS
              if (compareGeom(var1,target) == FALSE) {
                  # Resample to correct grid
                  var1 = resample(var1, target, method="ngb") ; gc() 
              }
              # Extend the extent of the overall grid to the analysis domain
              var1 = extend(var1,cardamom_ext)
              # Trim the extent of the overall grid to the analysis domain
              var1 = crop(var1,cardamom_ext)
              var1[which(as.vector(var1) < 0)] = NA
              # Adjust spatial resolution of the datasets, this occurs in all cases
              if (res(var1)[1] != res(cardamom_ext)[1] | res(var1)[2] != res(cardamom_ext)[2]) {

                  # Create raster with the target resolution
                  target = rast(crs = crs(cardamom_ext), ext = ext(cardamom_ext), resolution = res(cardamom_ext))
                  # Resample to correct grid
                  var1 = resample(var1, target, method="bilinear") ; gc() 

              } # Aggrgeate to resolution

              if (lat_done == FALSE) {
                  # extract dimension information for the grid, note the axis switching between raster and actual array
                  xdim = dim(var1)[2] ; ydim = dim(var1)[1]
                  # extract the lat / long information needed
                  long = crds(var1,df=TRUE, na.rm=FALSE)
                  lat  = long$y ; long = long$x
                  # restructure into correct orientation
                  long = array(long, dim=c(xdim,ydim))
                  lat = array(lat, dim=c(xdim,ydim))
                  loss_fraction = array(NA, dim=c(xdim,ydim,length(years_to_load)))
                  lat_done = TRUE
              }
              # break out from the rasters into arrays which we can manipulate
              var1 = array(as.vector(unlist(var1)), dim=c(xdim,ydim))

              # place new clearance information into the output array
              loss_fraction[,,yrr] = var1

          } # File exists
      } # looping years

      # output variables
      return(list(lat=lat,long=long,year_of_loss=years_to_load,loss_fraction=loss_fraction))

  } else {

      # output variables
      return(list(lat=-9999,long=-9999,year_of_loss=-9999,loss_fraction=-9999))

  }

} # end function load_forestry_fields_for_extraction

## Use byte compile
load_forestry_fields_for_extraction<-cmpfun(load_forestry_fields_for_extraction)
