#########################################################################################
# CARbon DAta MOdel fraMework (CARDAMOM) and DALEC terrestrial ecosystem model suite
# CARDAMOM is a Bayesian model-data fusion software framework. CARDAMOM is used to 
# assimilate observations and ecological theory to retrieve parameters for the 
# DALEC suite of intermediate complexity terrestrial ecosystem models. DALEC can be
# used as a fully integrated component of CARDAMOM or independently. 
# Copyright (C) 2024  University of Edinburgh,
#                     Mathew Williams (mat.williams@ed.ac.uk), 
#                     T. Luke Smallman (t.l.smallman@ed.ac.uk)
# UoE = University of Edinburgh

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# ########## File specific description ##########
# This function is to determine the area of a pixel at any location in meter squared.
# This function is based on an Python function development by J. F. Exbrayat (UoE).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).
# Note assume regular grid, though allows for differing lat and long resolutions
#
#########################################################################################

###
## Function to determine the area of a pixel at any location in meter squared
###

# This function is based on an Python function development by J. F. Exbrayat (UoE).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).
# Note assume regular grid, though allows for differing lat and long resolutions

## function to convert degrees to radians
deg2rad<-function(degree_in) {return(degree_in*(pi/180))}
rad2deg<-function(radian_in) {return(radian_in*(180/pi))}

calc_pixel_area<-function(longitude,latitude) {

    # Resolution in degrees
    # latitude (-90,90) degrees
    # longitude (-180,180) degrees

    # Deal with backwards compatibility for old 1D arrays of lon and lat
    # Assuming each provides single axis coordintes
    if(length(dim(longitude)) == 0) {
        longitude <- matrix(longitude, nrow=length(longitude), ncol=length(latitude))
    }
    if(length(dim(latitude)) == 0) {
        latitude <- matrix(latitude, nrow=dim(longitude)[1], ncol=dim(longitude)[2], byrow=TRUE)
    }
   
    # Assumptions:
    # 1) The grid is regular
    # 2) The provided lat / long are pixel centers
    # 3) That area at 90/-90 latitude is zero

    # Determine the longitude / latitude difference
    dlon = diff(longitude[,1]) ; dlat = diff(latitude[1,])
    # Ensure that the grid is regular
    #if (length(unique(dlon)) != 1) {stop("calc_pixel_area: longitude dimension not on regular grid")} TG
    #if (length(unique(dlat)) != 1) {stop("calc_pixel_area: latitude dimension not on regular grid")} TG
    # Extract resolution (lon/lat)
    resolution = c(dlon[1],dlat[1])
    
    # Ensure that resolution has two dimensions.
    # If only 1 provided assume equal size lat / long
    #if (length(resolution) == 1) {resolution = rep(resolution,2)}

    # mean earth radius (m)
    R = 6371e3 # 6378137 m
             
    # Estimate pixel area (m2)
    pixel_area = ( R**2 * deg2rad(resolution[1]) 
                        * (sin(deg2rad(latitude+resolution[2]*0.5))-sin(deg2rad(latitude-resolution[2]*0.5))) )

    # Correct for poles
    pixel_area[abs(latitude) > 90] = 0

    # return to user in meters
    return(pixel_area)

} # end function

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
#' Coded into the CARDAMOM framework by T. L. Smallman (t.l.smallman@ed.ac.uk)
#' 30/03/2023
#' Note function assumes a global grid is being fed in
#calc_pixel_area<- function(lon,lat) {
#
#    # Deal with backwards compatibility for old 1D arrays of lon and lat
#    # Assuming each provides single axis coordintes
#    if(length(dim(lon)) == 0) {
#        lon <- matrix(lon, nrow=length(lon), ncol=length(lat))
#    }
#    if(length(dim(lat)) == 0) {
#        lat <- matrix(lat, nrow=dim(lon)[1], ncol=dim(lon)[2], byrow=TRUE)
#    }
#
#    # Determine the dimensions of the overall grid
#    numLon <- dim(lon)[1]
#    numLat <- dim(lat)[2]
#
#    # If for some reason we have a -180:180 lon base, reset to span 0:360
#    lon[lon < 0] <- 360 + lon[lon < 0]
#    # Calculate the longitude degrees spanned by a grid cell
#    # ... modulo 360 to deal with wrapping boundries
#    deltaLon <- (lon[c(2:numLon,1),] - lon[1:(numLon-1),]) %% 360
#
#    # Calculate the min/max latitude for each grid cell
#    edgeLat <- (lat[,2:numLat]+lat[,2:numLat-1])/2
#    minLat <- cbind(-90, edgeLat)
#    maxLat <- cbind(edgeLat, 90)
#    # Check that the latitudes are centered in the grids
#    if(any(abs(lat) == 90)) {
#        warning('Grid cells centered at poles will have zero area.')
#        minLat[minLat < -90] <- -90
#        maxLat[maxLat > 90] <- 90
#    }
#
#    # Convert from degree to radius
#    deltaLon <- deg2rad(deltaLon)
#    minLat <- deg2rad(minLat)
#    maxLat <- deg2rad(maxLat)
#    lat <- deg2rad(lat)
#
#    # Assume the radius of the earth: 6371e3 meter
#    R <- 6371e3 # meters
#
#    # Calculate the east/west edges by assuming the earth is spherical and
#    # ...east/west edges are defined by latitude arc lengths
#    # ... => R*(maxLat-minLat)
#    # Calculate the north/south edges by assuming the arc length of longitude
#    # ...is the lattitude corrected radius (R*cos(lat)) times the change in lon
#    # ... => (R*cos(lat))*deltaLon
#    return( abs(R*(maxLat-minLat) * (R*cos(lat))*deltaLon))
#
#    # Old formulation for reference (updated 29 September 2014)
#    # ...no significant difference but harder to explain
#    #R^2*(sin(maxLat)-sin(minLat))*deltaLon
#
#} # calcGridArea
calc_pixel_area<-cmpfun(calc_pixel_area)