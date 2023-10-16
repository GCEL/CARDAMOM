
###
## Creates lat / long grid on WGS-84 framework
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

generate_wgs84_grid<-function(lat,long,resolution) {

    # Check whether we have two values for the resolution (longitude, latitude)
    if (length(resolution) == 1) {tmp = rep(resolution, length.out = 2)} else {tmp = resolution}

    # pt create raster of our desired resolution (in degrees) with the spatial extent set in the control file.
    # This function makes use of the WGS-84 coordinate system
    pt = rast(vals = 1, resolution = tmp, xmin = long[1], xmax = long[2], ymin = lat[1], ymax = lat[2], crs = "+init=epsg:4326")

    # extract the spatial information needed else where
    dims = dim(pt) ; lat_dim = dims[1] ; long_dim = dims[2]
    # extract the lat / long information needed
    long = crds(pt,df=TRUE, na.rm=FALSE)
    lat  = long$y ; long = long$x
    # reverse latitude vector to correct with R plotting orientation (i.e. make the world N/S rather than S/N)
    # It is very important to remember this has happened as the output from this fuction will not be consistent
    # with much of what is done with the calibration datasets.
    lat = lat[length(lat):1]

    # combine into output object
    output = list(cardamom_exts = pt, lat = lat, long = long, lat_dim = lat_dim, long_dim = long_dim)
    # tidy up
    rm(lat,long,lat_dim,long_dim,pt) ; gc(reset=TRUE,verbose=FALSE)

    # return back the solution
    return(output)

} # end function

## Use byte compile
generate_wgs84_grid<-cmpfun(generate_wgs84_grid)
