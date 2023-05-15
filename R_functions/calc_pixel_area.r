
###
## Function to determine the area of a pixel at any location in meter squared
###

# This function is based on an Python function development by J. F. Exbrayat (UoE).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

## function to convert degrees to radians
deg2rad<-function(degree_in) {return(degree_in*(pi/180))}
rad2deg<-function(radian_in) {return(radian_in*(180/pi))}

#calc_pixel_area<-function(latitude,longitude,resolution) {
#
#    # Resolution in degrees
#    # latitude (-90,90) degrees
#    # longitude (-180,180) degrees
#
#    # Ensure that resolution has two dimensions.
#    # If only 1 provided assume equal size lat / long
#    if (length(resolution) == 1) {resolution = rep(resolution,2)}
#
#    # mean earth radius (m)
#    R = 6371e3
#
#    # Estimate pixel area (m2)
#    pixel_area = ( R**2 * ( deg2rad(longitude+resolution[1]*0.5)-deg2rad(longitude-resolution[1]*0.5) )
#                        * (sin(deg2rad(latitude+resolution[2]*0.5))-sin(deg2rad(latitude-resolution[2]*0.5))) )
#
#    # return to user in meters
#    return(pixel_area)
#
#} # end function
### Use byte compile
#calc_pixel_area<-cmpfun(calc_pixel_area)

#' Calculate the grid cell area for a centered lat/lon grid.
#'
#' @param lat latitude coordinates of the grid centers
#' @param lon longitude coordinates of the grid centers
#' @return The grid cell area in m^2 (meter*meter)
#' @details Currently the lon must be uniform but the lat does not need to be.
#' @keywords internal
#' @note This is an internal RCMIP5 function and not exported.
#' Coded into the CARDAMOM framrwork by T. L. Smallman (t.l.smallman@ed.ac.uk)
#' 30/03/2023
calc_pixel_area<- function(lon,lat) {

    # Deal with backwards compatibility for old 1D arrays of lon and lat
    # Assuming each provides single axis coordintes
    if(length(dim(lon)) == 0) {
        lon <- matrix(lon, nrow=length(lon), ncol=length(lat))
    }
    if(length(dim(lat)) == 0) {
        lat <- matrix(lat, nrow=dim(lon)[1], ncol=dim(lon)[2], byrow=TRUE)
    }

    # Determine the dimensions of the overall grid
    numLon <- dim(lon)[1]
    numLat <- dim(lat)[2]

    # If for some reason we have a -180:180 lon base, reset to span 0:360
    lon[lon < 0] <- 360 + lon[lon < 0]
    # Calculate the longitude degrees spanned by a grid cell
    # ... modulo 360 to deal with wrapping boundries
    deltaLon <- (lon[c(2:numLon,1),] - lon[1:numLon,]) %% 360

    # Calculate the min/max latitude for each grid cell
    edgeLat <- (lat[,2:numLat]+lat[,2:numLat-1])/2
    minLat <- cbind(-90, edgeLat)
    maxLat <- cbind(edgeLat, 90)
    # Check that the latitudes are centered in the grids
    if(any(abs(lat) == 90)) {
        warning('Grid cells centered at poles will have zero area.')
        minLat[minLat < -90] <- -90
        maxLat[maxLat > 90] <- 90
    }

    # Convert from degree to radius
    deltaLon <- deg2rad(deltaLon)
    minLat <- deg2rad(minLat)
    maxLat <- deg2rad(maxLat)
    lat <- deg2rad(lat)

    # Assume the radius of the earth: 6371e3 meter
    R <- 6371e3 # meters

    # Calculate the east/west edges by assuming the earth is spherical and
    # ...east/west edges are defined by latitude arc lengths
    # ... => R*(maxLat-minLat)
    # Calculate the north/south edges by assuming the arc length of longitude
    # ...is the lattitude corrected radius (R*cos(lat)) times the change in lon
    # ... => (R*cos(lat))*deltaLon
    return( abs(R*(maxLat-minLat) * (R*cos(lat))*deltaLon))

    # Old formulation for reference (updated 29 September 2014)
    # ...no significant difference but harder to explain
    #R^2*(sin(maxLat)-sin(minLat))*deltaLon

} # calcGridArea
calc_pixel_area<-cmpfun(calc_pixel_area)
