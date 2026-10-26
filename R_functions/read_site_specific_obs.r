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
# Default function to read in site specific information for assimilation
# 
# Author: T. Luke Smallman (12/11/2024)
#
#########################################################################################

read_site_specific_obs <- function(variable,infile) {

    # read in the data assuming hearder are present and row do not have numbers
    input = read.csv(infile, header=TRUE)
    # check the variable exists in file
    if (length(which(names(input) == variable)) > 0) {
        # read from the table the desired informatin
        output = input[,variable]
#if (variable == "Evap_kgH2Om2day") {
#    tmp1=input[,"soilevap_kgH2Om2day"]
#    tmp2=input[,"wetevap_kgH2Om2day"]
#output = output - tmp1 - tmp2
#}
    } else {
        # if not then return missing value variable
        output = -9999
    }
    # return the desired variable
    return(output)

} # end of function read_site_specific_obs

## Use byte compile
read_site_specific_obs<-cmpfun(read_site_specific_obs)
