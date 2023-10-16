
###
## Load met function, called by met_fields_for_extraction.r
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_met_function<- function (year_to_do,varid,infile_varid,spatial_type,cardamom_ext,
                              path_to_met_source,met_source,wheat) {

    # Create target grid for aggregation if needed
    target = rast(crs = crs(cardamom_ext), ext = ext(cardamom_ext), resolution = res(cardamom_ext))

    if (met_source == "ERA" | met_source == "isimip3a") {

        # Loop through each monthly file
        t_grid = 0 ; var1_out = 0
        for (m in seq(1, 12)) {

             # Read in first file
             if (m > 9) {
                 input_file_1=paste(path_to_met_source,varid[1],"_",year_to_do,m,".nc",sep="")
             } else {
                 input_file_1=paste(path_to_met_source,varid[1],"_",year_to_do,"0",m,".nc",sep="")
             }
             # open netcdf files
             data1 = nc_open(input_file_1)
             # read the met drivers
             var_in = ncvar_get(data1, infile_varid[1])
             # keep count of time steps
             t_grid = t_grid + dim(var_in)[3]
             # read in the location information
             lat = ncvar_get(data1, "Latitude") ; long = ncvar_get(data1, "Longitude")
#             # convert input data long to conform to what we need
#             check1 = which(long > 180) ; if (length(check1) > 0) { long[check1] = long[check1]-360 }

             # expand the one directional values here into 2 directional
             lat_dim = length(lat) ; long_dim = length(long)
             long = array(long,dim=c(long_dim,lat_dim))
             lat = array(lat,dim=c(lat_dim,long_dim)) ; lat=t(lat)
             # close files after use
             nc_close(data1)

             # Loop each time step and aggregte before placing into to an out variables
             for (t in seq(1, dim(var_in)[3])) {
                  # Convert to a raster, assuming standad WGS84 grid
                  # This dependes on the lat / long / tmp1 spatially matching each other AND
                  # latitude ranging -90/90 and longitude ranging -180/180 degrees
                  var1 = data.frame(x = as.vector(long), y = as.vector(lat), z = as.vector(var_in[,,t]))
                  var1 = rast(var1, crs = ("+init=epsg:4326"), type="xyz")

                  # Extend the extent of the overall grid to the analysis domain
                  var1 = extend(var1,cardamom_ext)
                  # Trim the extent of the overall grid to the analysis domain
                  var1 = crop(var1,cardamom_ext)
                  var1[which(as.vector(var1) == -9999)] = NA
                  # Match resolutions
                  if (res(var1)[1] != res(cardamom_ext)[1] | res(var1)[2] != res(cardamom_ext)[2]) {
                      # Resample to correct grid
                      var1 = resample(var1, target, method="bilinear") ; gc() 
                  } # Aggrgeate to resolution

                  # break out from the rasters so can manipulate
                  var1_out = append(var1_out, as.vector(unlist(var1))[wheat])
             } # step within month loop

        }  # month loop
        
        # Remove initial value
        var1_out = var1_out[-1]

        # clean up
        rm(var1,m) ; gc(reset=TRUE,verbose=FALSE)

    } else if (met_source == "trendy_v9" | met_source == "trendy_v11" | met_source == "trendy_v12") {

        # Check whether this is a leap year or not
        nos_days = nos_days_in_year(year_to_do)
        # Determine what the number of days per month are
        # Define days in month
        days_per_month = rep(31,12) ; days_per_month[c(9,4,6,11)] = 30
        if (nos_days == 365) {days_per_month[2] = 28} else {days_per_month[2] = 29}

        # open first file in the sequence
        input_file_1 = paste(path_to_met_source,"/",varid[1],"_",year_to_do,"_monthly.nc",sep="")
        # open netcdf files
        data1 = nc_open(input_file_1)
        # read the met drivers
        var_in = ncvar_get(data1, infile_varid[1])
        # read in the location information
        lat = ncvar_get(data1, "lat") ; long = ncvar_get(data1, "lon")
#        # convert input data long to conform to what we need
#        check1 = which(long > 180) ; if (length(check1) > 0) { long[check1] = long[check1]-360 }

        # expand the one directional values here into 2 directional
        lat_dim = length(lat) ; long_dim = length(long)
        long = array(long,dim=c(long_dim,lat_dim))
        lat = array(lat,dim=c(lat_dim,long_dim)) ; lat=t(lat)
        # close files after use
        nc_close(data1)

        # Loop each time step and aggregte before placing into to an out variables
        var1_out = 0
        for (t in seq(1, dim(var_in)[3])) {
             # Convert to a raster, assuming standad WGS84 grid
             # This dependes on the lat / long / tmp1 spatially matching each other AND
             # latitude ranging -90/90 and longitude ranging -180/180 degrees
             var1 = data.frame(x = as.vector(long), y = as.vector(lat), z = as.vector(var_in[,,t]))
             var1 = rast(var1, crs = ("+init=epsg:4326"), type="xyz")

             # Extend the extent of the overall grid to the analysis domain
             var1 = extend(var1,cardamom_ext)
             # Trim the extent of the overall grid to the analysis domain
             var1 = crop(var1,cardamom_ext)
             var1[which(as.vector(var1) == -9999)] = NA
             # Match resolutions of the datasets
             if (res(var1)[1] != res(cardamom_ext)[1] | res(var1)[2] != res(cardamom_ext)[2]) {
                 # Resample to correct grid
                 var1 = resample(var1, target, method="bilinear") ; gc() 
             } # Aggrgeate to resolution

             # break out from the rasters so can manipulate
             var1_out = append(var1_out, rep(as.vector(unlist(var1))[wheat], times = days_per_month[t]))

        } # end of time loop

        # Remove initial value
        var1_out = var1_out[-1]

        # keep count of time steps
        t_grid = sum(days_per_month)

        # clean up
        rm(var_in,var1,t) ; gc(reset=TRUE,verbose=FALSE)

    } # end data source selection

    # return back to the user
    return(list(var_out=var1_out,t_grid=t_grid))

} # end function
## Use byte compile
load_met_function<-cmpfun(load_met_function)
