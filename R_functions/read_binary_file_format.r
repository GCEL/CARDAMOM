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
# Function to extract met drivers from binary file
# 
# This function is based on an original Matlab function development by A. A. Bloom 
# (UoE, now at the Jet Propulsion Laboratory). Translation to R and subsequent 
# modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).
# Exceptions are given within specific functions.
#
#########################################################################################

read_binary_file_format<- function(infile) {

      #print("Beginning read of binary input files...")
      # Open and read the DALEC binary driver file
      # open this chains binary file into R, instructing 'r' to read and 'b' for binary
      bob=file(infile,'rb') ; nos_var=1e6
      bd=readBin(bob, double(),nos_var)
      # keep reading until we have read all that can be read
      set1=NA
      while (length(set1) > 0) {
	      set1=readBin(bob, double(),nos_var)
	      bd=append(bd,set1)
      }
      # now close this chain
      close(bob)

      # begin preparing to disagregate the different sections of the file
      k=0
      # extract static data (50 places)
      static=bd[(k+1):(k+50)]
      k=k+50
      # extract priors
      pr=bd[(k+1):(k+100)]
      k=k+100
      # extract prior uncertainties (100 places)
      pru=bd[(k+1):(k+100)]
      k=k+100
      # extract prior weighting (100 places)
      prw=bd[(k+1):(k+100)]
      k=k+100
      # other priors (50 places)
      opr=bd[(k+1):(k+50)]
      k=k+50
      # other prior uncertainties (50 places)
      opru=bd[(k+1):(k+50)]
      k=k+50
      # other prior weighting (50 places)
      oprw=bd[(k+1):(k+50)]
      k=k+50

      # store prior information
      md=list(parpriors=pr,parpriorunc=pru,parpriorweight=prw,otherpriors=opr,otherpriorunc=opru,otherpriorweight=oprw)

      # id code (not currently used)
      md$id = static[1]
      # latitude
      md$lat = static[2]
      # number of days in simulation
      md$nodays = static[3]
      # number of met fields
      md$nomet = static[4]
      # number of observation streams
      md$noobs = static[5]
      # implement reality checks (1=yes,0=no)
      md$RC = static[6]
      # ctessel pft
      md$ctessel_pft = static[7]
      # yield if forest
      md$yield = static[8]
      # age if forest
      md$age = static[9]
      # 10 will be nos_pars
      ###
      # start searching EDCs from anywhere (1) or from prescribed starting point (0)
      md$rc_random_search = static[11] == 1
      # Top sand %
      md$top_sand = static[12]
      # bot sand %
      md$bot_sand = static[13]
      # Top clay %
      md$top_clay = static[14]
      # Bot clay %
      md$bot_clay = static[15]
      # Top residual water fraction (m3m3) from Van-Genuchten Mualem
      md$top_residual_waterfrac = static[16]
      # Bottom residual water fraction (m3m3) from Van-Genuchten Mualem
      md$bot_residual_waterfrac = static[17]
      # Top porosity (m3m3) from Van-Genuchten Mualem
      md$top_porosity = static[18]
      # Bottom porosity (m3m3) from Van-Genuchten Mualem
      md$bot_porosity = static[19]
      # Top pore size distribution (-) from Van-Genuchten Mualem
      md$top_pore_size_dist = static[20]
      # Bottom pore size distribution (-) from Van-Genuchten Mualem
      md$bot_pore_size_dist = static[21]
      # Top air entry pressure (m-1) from Van-Genuchten Mualem
      md$top_air_entry = static[22]
      # Bottom air entry pressure (m-1) from Van-Genuchten Mualem
      md$bot_air_entry = static[23]
      # Top saturated hydraulic conductivity (ms-1) from Van-Genuchten Mualem
      md$top_sat_conductivity = static[24]
      # Bottom saturated hydraulic conductivity (ms-1) from Van-Genuchten Mualem
      md$bot_sat_conductivity = static[25]


      # extract temporal data (met and obs)
      tempdata = bd[(k+1):(k+((md$nomet+md$noobs)*md$nodays))]
      # restructure correctly
      tempdata = array(tempdata,dim=c(md$nomet+md$noobs,md$nodays))

      # pass all met data
      md$met = t(tempdata[1:md$nomet,1:md$nodays])
      # pass all observations (if any)
      # current defaults are:
      # 1st column = GPP
      # 2nd column = LAI
      # 3rd column = NEE
      md$obs = t(tempdata[(md$nomet+1):(md$nomet+md$noobs),1:md$nodays])

      # pass back out information
      return(md)

} # end of function
## Use byte compile
read_binary_file_format<-cmpfun(read_binary_file_format)
