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
# Function to create regular latitude / longitude grids for global WGS-84.
# 
# Author: T. Luke Smallman (12/11/2024)
#
#########################################################################################

generate_grid<-function(cardamom_grid_type,lat,long,resolution) {

    # Check whether we have two values for the resolution (longitude, latitude)
    if (length(resolution) == 1) {tmp = rep(resolution, length.out = 2)} else {tmp = resolution}

    # pt create raster of our desired resolution (e.g. degrees) with the spatial extent set in the control file.
    # For example, espg:4326 is the WGS-84 coordinate system
    pt = rast(vals = 1, resolution = tmp, xmin = long[1], xmax = long[2], ymin = lat[1], ymax = lat[2], crs = cardamom_grid_type)
    # Estimate the pixel areas
    area = cellSize(pt, unit="m", transform=FALSE, rcx=100)

    # extract the spatial information needed else where
    dims = dim(pt) ; lat_dim = dims[1] ; long_dim = dims[2]
    # extract the lat / long information needed
    long = crds(pt,df=TRUE, na.rm=FALSE)
    lat  = long$y ; long = long$x
    # restructure into correct orientation
    long = array(long, dim=c(long_dim,lat_dim))
    lat = array(lat, dim=c(long_dim,lat_dim))
    # reverse latitude vector to correct with R plotting orientation (i.e. make the world N/S rather than S/N)
    # It is very important to remember this has happened as the output from this fuction will not be consistent
    # with much of what is done with the calibration datasets.
    long = long[,lat_dim:1] ; lat = lat[,lat_dim:1]

    # Extract the fractional cover, 
    # assumed to be unchanged for any given year.
    area = array(values(area), dim=c(long_dim,lat_dim))
    # Flip the north-south axis to invert for human eyes
    area = area[,lat_dim:1]

    # combine into output object
    output = list(cardamom_ext = pt, lat = lat, long = long, lat_dim = lat_dim, long_dim = long_dim, area = area)
    # tidy up
    rm(lat,long,lat_dim,long_dim,pt) ; gc(reset=TRUE,verbose=FALSE)

    # return back the solution
    return(output)

} # end function

## Use byte compile
generate_grid<-cmpfun(generate_grid)