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
# Function to generate mean state variable information by running the parameters and model choice
# 
# Author: T. Luke Smallman (12/11/2024)
# Exceptions states below in specific functions
#
#########################################################################################

# Function to actually restrict values of a double precision numeric into single precision allowed range
# Note: modified from "readBrukerFlexData" library
double2single <-function(x) {

    # Ensure the incoming variable is a double precision variable
    stopifnot(is.double(x))
    # Create a virtual connection which we will use to write / read the single precision variable version
    virtualCon = raw()
    # Write out to virtual connection in 4L (i.e. 4 byte single precision)
    virtualCon = writeBin(object = x, con = virtualCon, size = 4L) # 4L defines single precision
    # Read back into the new variable. what = double() as R does not have a format for single precision
    y = readBin(con = virtualCon, what = double(), size = 4L, n = length(x))
    # Return back to user
    return(y)

} # end function double2single

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

simulate_all<- function (site,PROJECT,model_name,met,pars,lat,pft,parameter_type,exepath,soil_info) {

  output_dim=17 ; aNPP_dim = 3 ; MTT_dim = 5 ; SS_dim = 5 ; fire_dim = 6
  noyears = length(c(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)))

  # restructure pars
  if (length(dim(pars)) > 2) {
      pars_in = array(0,dim=c(dim(pars)[1],dim(pars)[2]*dim(pars)[3]))
      nos_iter = dim(pars)[2]*dim(pars)[3]
      for (n in seq(1,dim(pars)[1])){
           pars_in[n,] = pars[n,,]
      }
  } else if (length(dim(pars)) == 2) {
      nos_iter = dim(pars)[2]
      pars_in = pars
  } else {
      nos_iter = 1
      pars_in = array(pars, dim=c(length(pars),1))
  }

  # loop through combinations
  if (model_name == "ACM") {
      # load the function code from the dalec shared object.
      # NOTE: that the name of the shared object is hardcoded to dalec.so
      # while the function call within will be specific to the actual model being called
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      output_dim = 11
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "racm",output_dim=as.integer(output_dim),met=as.double(t(met)),pars=as.double(pars_in)
                          ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                          ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                          ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                          ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                          ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                          ,soil_frac_clay=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                          ,soil_frac_sand=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3))) )
      output=tmp$out_var
      output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # construct output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      ET_kgH2Om2day = output[,,3],soilET_kgH2Om2day = output[,,4],
                      Rtot_MPasm2mmol = output[,,5],wetcanET_kgH2Om2day = output[,,6],
                      gs_demand_supply = output[,,7], gs_total_canopy = output[,,8],
                      APAR_MJm2day = output[,,9], gb_total_canopy = output[,,10],
                      CiCa = output[,,11])
  } else if (model_name == "DALEC.C3.M1.014") {
      output_dim = 58 ; MTT_dim = 8 ; SS_dim = 8
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      #crop_file_location=paste(PROJECT$exepath,"winter_wheat_development.csv", sep="")
      crop_type = 1 # Winter Wheat
      wd_old = getwd() ; setwd(PROJECT$exepath)
      tmp=.Fortran( "rdalec14",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                              ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                             
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,nos_years=as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                              ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2]))
                              ,pathlength=as.integer(crop_type))
                              #,exepath=as.character(crop_file_location),pathlength=as.integer(nchar(crop_file_location)))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc() ; setwd(wd_old)
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      harvest_gCm2day = output[,,5],
                      mean_harvest_gCm2day = output_mean[,5],
                      mean_annual_harvest_gCm2day = output_annual[,,5],
                      extracted_residue_gCm2day = output[,,6],
                      mean_extracted_residue_gCm2day = output_mean[,6],
                      mean_annual_extracted_residue_gCm2day = output_annual[,,6],
                      fire_gCm2day = array(0, dim(output)[1:2]), # kept because fire is a common management routine outside of europe
                      mean_fire_gCm2day = array(0, dim(output_mean)[1]), # kept because fire is a common management routine outside of europe
                      mean_annual_fire_gCm2day = array(0, dim(output_annual)[1:2]), # kept because fire is a common management routine outside of europe
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      alloc_autotrophic_gCm2day = output[,,12],
                      mean_alloc_autotrophic_gCm2day = output_mean[,12],
                      mean_annual_alloc_autotrophic_gCm2day = output_annual[,,12],
                      alloc_StorageOrgan_gCm2day = output[,,13],
                      mean_alloc_StorageOrgan_gCm2day = output_mean[,13],
                      mean_annual_alloc_StorageOrgan_gCm2day = output_annual[,,13],
                      foliage_to_litter_gCm2day = output[,,14],
                      mean_foliage_to_litter_gCm2day = output_mean[,14],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,14],
                      roots_to_litter_gCm2day = output[,,15],
                      mean_roots_to_litter_gCm2day = output_mean[,15],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,15],
                      wood_to_litter_gCm2day = output[,,16],
                      mean_wood_to_litter_gCm2day = output_mean[,16],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,16],
                      litter_to_som_gCm2day = output[,,17],
                      mean_litter_to_som_gCm2day = output_mean[,17],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,17],
                      rauto_maintenance_gCm2day = output[,,18],
                      mean_rauto_maintenance_gCm2day = output_mean[,18],
                      mean_annual_rauto_maintenance_gCm2day = output_annual[,,18],
                      rauto_labile_to_foliage_gCm2day = output[,,19],
                      mean_rauto_labile_to_foliage_gCm2day = output_mean[,19],
                      mean_annual_rauto_labile_to_foliage_gCm2day = output_annual[,,19],
                      rauto_npp_to_labile_gCm2day = output[,,20],
                      mean_rauto_npp_to_labile_gCm2day = output_mean[,20],
                      mean_annual_rauto_npp_to_labile_gCm2day = output_annual[,,20],
                      rauto_foliage_to_litter_gCm2day = output[,,21],
                      mean_rauto_foliage_to_litter_gCm2day = output_mean[,21],
                      mean_annual_rauto_foliage_to_litter_gCm2day = output_annual[,,21],
                      rauto_wood_to_litter_gCm2day = output[,,22],
                      mean_rauto_wood_to_litter_gCm2day = output_mean[,22],
                      mean_annual_rauto_wood_to_litter_gCm2day = output_annual[,,22],
                      HARVESTextracted_foliage_gCm2day = output[,,23],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,23],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,23],
                      HARVESTextracted_wood_gCm2day = output[,,24],
                      mean_HARVESTextracted_wood_gCm2day = output_annual[,,24],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,24],
                      HARVESTextracted_DeadFoliage_gCm2day = output[,,25],
                      mean_HARVESTextracted_DeadFoliage_gCm2day = output_mean[,25],
                      mean_annual_HARVESTextracted_DeadFoliage_gCm2day = output_annual[,,25],
                      HARVESTextracted_labile_gCm2day = output[,,26],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,26],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,26],
                      HARVESTlitter_foliage_gCm2day = output[,,27],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,27],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,27],
                      HARVESTlitter_wood_gCm2day = output[,,28],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,28],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,28],
                      HARVESTlitter_DeadFoliage_gCm2day = output[,,29],
                      mean_HARVESTlitter_DeadFoliage_gCm2day = output_mean[,29],
                      mean_annual_HARVESTlitter_DeadFoliage_gCm2day = output_annual[,,29],
                      HARVESTlitter_autotrophic_gCm2day = output[,,30],
                      mean_HARVESTlitter_autotrophic_gCm2day = output_mean[,30],
                      mean_annual_HARVESTlitter_autotrophic_gCm2day = output_annual[,,30],
                      HARVESTlitter_labile_gCm2day = output[,,31],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,31],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,31],
                      PLOUGHlitter_roots_gCm2day = output[,,32],
                      mean_PLOUGHlitter_roots_gCm2day = output_mean[,32],
                      mean_annual_PLOUGHlitter_roots_gCm2day = output_annual[,,32],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,33],
                      mean_labile_gCm2 = output_mean[,33],
                      mean_annual_labile_gCm2 = output_annual[,,33],
                      foliage_gCm2 = output[,,34],
                      mean_foliage_gCm2 = output_mean[,34],
                      mean_annual_foliage_gCm2 = output_annual[,,34],
                      roots_gCm2 = output[,,35],
                      mean_roots_gCm2 = output_mean[,35],
                      mean_annual_roots_gCm2 = output_annual[,,35],
                      wood_gCm2 = output[,,36],
                      mean_wood_gCm2 = output_mean[,36],
                      mean_annual_wood_gCm2 = output_annual[,,36],
                      litter_gCm2 = output[,,37],
                      mean_litter_gCm2 = output_mean[,37],
                      mean_annual_litter_gCm2 = output_annual[,,37],
                      som_gCm2 = output[,,38],
                      mean_som_gCm2 = output_mean[,38],
                      mean_annual_som_gCm2 = output_annual[,,38],
                      autotrophic_gCm2 = output[,,39],
                      mean_autotrophic_gCm2 = output_mean[,39],
                      mean_annual_autotrophic_gCm2 = output_annual[,,39],
                      StorageOrgan_gCm2 = output[,,40],
                      mean_StorageOrgan_gCm2 = output_mean[,40],
                      mean_annual_StorageOrgan_gCm2 = output_annual[,,40],
                      DeadFoliage_gCm2 = output[,,41],
                      mean_DeadFoliage_gCm2 = output_mean[,41],
                      mean_annual_DeadFoliage_gCm2 = output_annual[,,41],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,42],
                      mean_lai_m2m2 = output_mean[,42],
                      mean_annual_lai_m2m2 = output_annual[,,42],
                      # Photosynthesis / C~water coupling related
                      CiCa = output[,,43],
                      mean_CiCa = output_mean[,43],
                      mean_annual_CiCa = output_annual[,,43],
                      # Misc
                      DevelopmentStage = output[,,44],
                      mean_DevelopmentStage = output_mean[,44],
                      mean_annual_DevelopmentStage = output_annual[,,44],
                      FoliarN_gNm2 = output[,,45],
                      mean_FoliarN_gNm2  = output_mean[,45],
                      mean_annual_FoliarN_gNm2  = output_annual[,,45],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      MTT_autotrophic_years = MTT_years[,7],
                      MTT_DeadFoliage_years = MTT_years[,8],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6],
                      SS_autotrophic_gCm2 = SS_gCm2[,7],
                      SS_DeadFoliage_gCm2 = SS_gCm2[,8])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day +
                           states_all$alloc_autotrophic_gCm2day +
                           states_all$alloc_StorageOrgan_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean),
                           apply(states_all$alloc_autotrophic_gCm2day,1,mean),
                           apply(states_all$alloc_StorageOrgan_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      states_all$NPP_autotrophic_fraction = NPP_fraction[,4]
      states_all$NPP_StorageOrgan_fraction = NPP_fraction[,5]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A3.C3.H2.M1.015") {
      output_dim = 62 ; MTT_dim = 8 ; SS_dim = 8
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      #crop_file_location=paste(PROJECT$exepath,"winter_wheat_development.csv", sep="")
      crop_type = 1 # Winter Wheat
      wd_old = getwd() ; setwd(PROJECT$exepath)
      tmp=.Fortran( "rdalec15",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                             ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                             
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodiags=as.integer(PROJECT$model$nodiags[site]),nodays=as.integer(dim(met)[1])
                             ,nos_years=as.integer(noyears)
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                             ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                             ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2]))
                             ,residual_waterfrac_in=as.double(c(soil_info[5],soil_info[6],soil_info[6])) 
                             ,porosity_in = as.double(c(soil_info[7],soil_info[8],soil_info[8]))
                             ,pore_size_dist_in = as.double(c(soil_info[9],soil_info[10],soil_info[10]))
                             ,air_entry_in = as.double(c(soil_info[11],soil_info[12],soil_info[12]))
                             ,sat_conductivity_in = as.double(c(soil_info[13],soil_info[14],soil_info[14]))                             
                             ,pathlength=as.integer(crop_type))
                             #,exepath=as.character(crop_file_location),pathlength=as.integer(nchar(crop_file_location)))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc() ; setwd(wd_old)
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      harvest_gCm2day = output[,,5],
                      mean_harvest_gCm2day = output_mean[,5],
                      mean_annual_harvest_gCm2day = output_annual[,,5],
                      extracted_residue_gCm2day = output[,,6],
                      mean_extracted_residue_gCm2day = output_mean[,6],
                      mean_annual_extracted_residue_gCm2day = output_annual[,,6],
                      fire_gCm2day = array(0, dim(output)[1:2]), # kept because fire is a common management routine outside of europe
                      mean_fire_gCm2day = array(0, dim(output_mean)[1]), # kept because fire is a common management routine outside of europe
                      mean_annual_fire_gCm2day = array(0, dim(output_annual)[1:2]), # kept because fire is a common management routine outside of europe
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      alloc_autotrophic_gCm2day = output[,,12],
                      mean_alloc_autotrophic_gCm2day = output_mean[,12],
                      mean_annual_alloc_autotrophic_gCm2day = output_annual[,,12],
                      alloc_StorageOrgan_gCm2day = output[,,13],
                      mean_alloc_StorageOrgan_gCm2day = output_mean[,13],
                      mean_annual_alloc_StorageOrgan_gCm2day = output_annual[,,13],
                      foliage_to_litter_gCm2day = output[,,14],
                      mean_foliage_to_litter_gCm2day = output_mean[,14],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,14],
                      roots_to_litter_gCm2day = output[,,15],
                      mean_roots_to_litter_gCm2day = output_mean[,15],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,15],
                      wood_to_litter_gCm2day = output[,,16],
                      mean_wood_to_litter_gCm2day = output_mean[,16],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,16],
                      litter_to_som_gCm2day = output[,,17],
                      mean_litter_to_som_gCm2day = output_mean[,17],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,17],
                      rauto_maintenance_gCm2day = output[,,18],
                      mean_rauto_maintenance_gCm2day = output_mean[,18],
                      mean_annual_rauto_maintenance_gCm2day = output_annual[,,18],
                      rauto_labile_to_foliage_gCm2day = output[,,19],
                      mean_rauto_labile_to_foliage_gCm2day = output_mean[,19],
                      mean_annual_rauto_labile_to_foliage_gCm2day = output_annual[,,19],
                      rauto_npp_to_labile_gCm2day = output[,,20],
                      mean_rauto_npp_to_labile_gCm2day = output_mean[,20],
                      mean_annual_rauto_npp_to_labile_gCm2day = output_annual[,,20],
                      rauto_foliage_to_litter_gCm2day = output[,,21],
                      mean_rauto_foliage_to_litter_gCm2day = output_mean[,21],
                      mean_annual_rauto_foliage_to_litter_gCm2day = output_annual[,,21],
                      rauto_wood_to_litter_gCm2day = output[,,22],
                      mean_rauto_wood_to_litter_gCm2day = output_mean[,22],
                      mean_annual_rauto_wood_to_litter_gCm2day = output_annual[,,22],
                      HARVESTextracted_foliage_gCm2day = output[,,23],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,23],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,23],
                      HARVESTextracted_wood_gCm2day = output[,,24],
                      mean_HARVESTextracted_wood_gCm2day = output_annual[,,24],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,24],
                      HARVESTextracted_DeadFoliage_gCm2day = output[,,25],
                      mean_HARVESTextracted_DeadFoliage_gCm2day = output_mean[,25],
                      mean_annual_HARVESTextracted_DeadFoliage_gCm2day = output_annual[,,25],
                      HARVESTextracted_labile_gCm2day = output[,,26],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,26],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,26],
                      HARVESTlitter_foliage_gCm2day = output[,,27],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,27],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,27],
                      HARVESTlitter_wood_gCm2day = output[,,28],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,28],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,28],
                      HARVESTlitter_DeadFoliage_gCm2day = output[,,29],
                      mean_HARVESTlitter_DeadFoliage_gCm2day = output_mean[,29],
                      mean_annual_HARVESTlitter_DeadFoliage_gCm2day = output_annual[,,29],
                      HARVESTlitter_autotrophic_gCm2day = output[,,30],
                      mean_HARVESTlitter_autotrophic_gCm2day = output_mean[,30],
                      mean_annual_HARVESTlitter_autotrophic_gCm2day = output_annual[,,30],
                      HARVESTlitter_labile_gCm2day = output[,,31],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,31],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,31],
                      PLOUGHlitter_roots_gCm2day = output[,,32],
                      mean_PLOUGHlitter_roots_gCm2day = output_mean[,32],
                      mean_annual_PLOUGHlitter_roots_gCm2day = output_annual[,,32],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,33],
                      mean_labile_gCm2 = output_mean[,33],
                      mean_annual_labile_gCm2 = output_annual[,,33],
                      foliage_gCm2 = output[,,34],
                      mean_foliage_gCm2 = output_mean[,34],
                      mean_annual_foliage_gCm2 = output_annual[,,34],
                      roots_gCm2 = output[,,35],
                      mean_roots_gCm2 = output_mean[,35],
                      mean_annual_roots_gCm2 = output_annual[,,35],
                      wood_gCm2 = output[,,36],
                      mean_wood_gCm2 = output_mean[,36],
                      mean_annual_wood_gCm2 = output_annual[,,36],
                      litter_gCm2 = output[,,37],
                      mean_litter_gCm2 = output_mean[,37],
                      mean_annual_litter_gCm2 = output_annual[,,37],
                      som_gCm2 = output[,,38],
                      mean_som_gCm2 = output_mean[,38],
                      mean_annual_som_gCm2 = output_annual[,,38],
                      autotrophic_gCm2 = output[,,39],
                      mean_autotrophic_gCm2 = output_mean[,39],
                      mean_annual_autotrophic_gCm2 = output_annual[,,39],
                      StorageOrgan_gCm2 = output[,,40],
                      mean_StorageOrgan_gCm2 = output_mean[,40],
                      mean_annual_StorageOrgan_gCm2 = output_annual[,,40],
                      DeadFoliage_gCm2 = output[,,41],
                      mean_DeadFoliage_gCm2 = output_mean[,41],
                      mean_annual_DeadFoliage_gCm2 = output_annual[,,41],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,42],
                      mean_ET_kgH2Om2day = output_mean[,42],
                      mean_annual_ET_kgH2Om2day = output_annual[,,42],
                      Etrans_kgH2Om2day = output[,,43],
                      mean_Etrans_kgH2Om2day = output_mean[,43],
                      mean_annual_Etrans_kgH2Om2day = output_annual[,,43],
                      Esoil_kgH2Om2day = output[,,44],
                      mean_Esoil_kgH2Om2day = output_mean[,44],
                      mean_annual_Esoil_kgH2Om2day = output_annual[,,44],
                      Ewetcanopy_kgH2Om2day = output[,,45],
                      mean_Ewetcanopy_kgH2Om2day = output_mean[,45],
                      mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,45],
                      runoff_kgH2Om2day = output[,,46],
                      mean_runoff_kgH2Om2day = output_mean[,46],
                      mean_annual_runoff_kgH2Om2day = output_annual[,,46],
                      underflow_kgH2Om2day = output[,,47],
                      mean_underflow_kgH2Om2day = output_mean[,47],
                      mean_annual_underflow_kgH2Om2day = output_annual[,,47],
                      SurfDrainage_kgH2Om2day = output[,,48],
                      mean_SurfDrainage_kgH2Om2day = output_mean[,48],
                      mean_annual_SurfDrainage_kgH2Om2day = output_annual[,,48],
                      SurfInfiltrated_kgH2Om2day = output[,,49],
                      mean_SurfInfiltrated_kgH2Om2day = output_mean[,49],
                      mean_annual_SurfInfiltrated_kgH2Om2day = output_annual[,,49],
                      Etrans_1st_root_layer_uptake_fraction = output[,,50],
                      mean_Etrans_1st_root_layer_uptake_fraction = output_mean[,50],
                      mean_annual_Etrans_1st_root_layer_uptake_fraction = output_annual[,,50],
                      Etrans_2nd_root_layer_uptake_fraction = output[,,51],
                      mean_Etrans_2nd_root_layer_uptake_fraction = output_mean[,51],
                      mean_annual_Etrans_2nd_root_layer_uptake_fraction = output_annual[,,51],
                      SurfWater_kgH2Om2 = output[,,52],
                      mean_SurfWater_kgH2Om2 = output_mean[,52],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,52],
                      wSWP_MPa = output[,,53],
                      mean_wSWP_MPa = output_mean[,53],
                      mean_annual_wSWP_MPa = output_annual[,,53],
                      snow_kgH2Om2 = output[,,54],
                      mean_snow_kgH2Om2 = output_mean[,54],
                      mean_annual_snow_kgH2Om2 = output_annual[,,54],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,55],
                      mean_lai_m2m2 = output_mean[,55],
                      mean_annual_lai_m2m2 = output_annual[,,55],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,56],
                      mean_gs_demand_supply_ratio = output_mean[,56],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,56],
                      gs_mmolH2Om2s = output[,,57],
                      mean_gs_mmolH2Om2s = output_mean[,57],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,57],
                      APAR_MJm2day = output[,,58],
                      mean_APAR_MJm2day = output_mean[,58],
                      mean_annual_APAR_MJm2day = output_annual[,,58],
                      gb_mmolH2Om2s = output[,,59],
                      mean_gb_mmolH2Om2s = output_mean[,59],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,59],
                      CiCa = output[,,60],
                      mean_CiCa = output_mean[,60],
                      mean_annual_CiCa = output_annual[,,60],
                      # Misc
                      RootDepth_m = output[,,61],
                      mean_RootDepth_m = output_mean[,61],
                      mean_annual_RootDepth_m = output_annual[,,61],
                      DevelopmentStage = output[,,62],
                      mean_DevelopmentStage = output_mean[,62],
                      mean_annual_DevelopmentStage = output_annual[,,62],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      MTT_autotrophic_years = MTT_years[,7],
                      MTT_DeadFoliage_years = MTT_years[,8],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6],
                      SS_autotrophic_gCm2 = SS_gCm2[,7],
                      SS_DeadFoliage_gCm2 = SS_gCm2[,8])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day +
                           states_all$alloc_autotrophic_gCm2day +
                           states_all$alloc_StorageOrgan_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean),
                           apply(states_all$alloc_autotrophic_gCm2day,1,mean),
                           apply(states_all$alloc_StorageOrgan_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      states_all$NPP_autotrophic_fraction = NPP_fraction[,4]
      states_all$NPP_StorageOrgan_fraction = NPP_fraction[,5]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P4.R2.011") {
    output_dim = 72 ; MTT_dim = 7 ; SS_dim = 7
    # Load the required dalec shared object
    dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
    tmp=.Fortran( "rdalec11",output_dim=as.integer(output_dim)
                            ,aNPP_dim=as.integer(aNPP_dim)
                            ,MTT_dim=as.integer(MTT_dim)
                            ,SS_dim = as.integer(SS_dim)
                            ,fire_dim=as.integer(fire_dim)
                            ,met=as.double(t(met))
                            ,pars=as.double(pars_in)
                            ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                            ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                            ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                            ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                            ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                            
                            ,lat=as.double(lat)
                            ,nopars=as.integer(PROJECT$model$nopars[site])
                            ,nomet=as.integer(dim(met)[2])
                            ,nofluxes=as.integer(PROJECT$model$nofluxes[site])
                            ,nopools=as.integer(PROJECT$model$nopools[site])
                            ,pft=as.integer(pft)
                            ,nodays=as.integer(dim(met)[1])
                            ,nos_years=as.integer(noyears)
                            ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1]))))
                            ,nos_iter=as.integer(nos_iter)
                            ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                            ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3)))
                            ,residual_waterfrac_in=as.double(c(soil_info[5],soil_info[6],soil_info[6])) 
                            ,porosity_in = as.double(c(soil_info[7],soil_info[8],soil_info[8]))
                            ,pore_size_dist_in = as.double(c(soil_info[9],soil_info[10],soil_info[10]))
                            ,air_entry_in = as.double(c(soil_info[11],soil_info[12],soil_info[12]))
                            ,sat_conductivity_in = as.double(c(soil_info[13],soil_info[14],soil_info[14])))
    # Extract the different output variables
    output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
    output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
    output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))
    MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
    SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
    # Unload the current dalec shared object
    dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
    rm(tmp) ; gc()
    # create output object
    states_all=list(# Ecosystem fluxes
                    gpp_gCm2day = output[,,1],
                    mean_gpp_gCm2day = output_mean[,1],
                    mean_annual_gpp_gCm2day = output_annual[,,1],
                    rauto_gCm2day = output[,,2],
                    mean_rauto_gCm2day = output_mean[,2],
                    mean_annual_rauto_gCm2day = output_annual[,,2],
                    rg_foliage_gCm2day = output[,,3],
                    mean_rg_foliage_gCm2day = output_mean[,3],
                    mean_annual_rg_foliage_gCm2day = output_annual[,,3],
                    rhet_litter_gCm2day = output[,,4],
                    mean_rhet_litter_gCm2day = output_mean[,4],
                    mean_annual_rhet_litter_gCm2day = output_annual[,,4],
                    rhet_som_gCm2day = output[,,5],
                    mean_rhet_som_gCm2day = output_mean[,5],
                    mean_annual_rhet_som_gCm2day = output_annual[,,5],
                    rhet_woodlitter_gCm2day = output[,,6],
                    mean_rhet_woodlitter_gCm2day = output_mean[,6],
                    mean_annual_rhet_woodlitter_gCm2day = output_annual[,,6],
                    fire_gCm2day = output[,,7],
                    mean_fire_gCm2day = output_mean[,7],
                    mean_annual_fire_gCm2day = output_annual[,,7],
                    harvest_gCm2day = output[,,8],
                    mean_harvest_gCm2day = output_mean[,8],
                    mean_annual_harvest_gCm2day = output_annual[,,8],
                    # Internal fluxes
                    alloc_labile_gCm2day = output[,,9],
                    mean_alloc_labile_gCm2day = output_mean[,9],
                    mean_annual_alloc_labile_gCm2day = output_annual[,,9],
                    alloc_roots_gCm2day = output[,,10],
                    mean_alloc_roots_gCm2day = output_mean[,10],
                    mean_annual_alloc_roots_gCm2day = output_annual[,,10],
                    alloc_wood_gCm2day = output[,,11],
                    mean_alloc_wood_gCm2day = output_mean[,11],
                    mean_annual_alloc_wood_gCm2day = output_annual[,,11],
                    labile_to_foliage_gCm2day = output[,,12],
                    mean_labile_to_foliage_gCm2day = output_mean[,12],
                    mean_annual_labile_to_foliage_gCm2day = output_annual[,,12],
                    foliage_to_litter_gCm2day = output[,,13],
                    mean_foliage_to_litter_gCm2day = output_mean[,13],
                    mean_annual_foliage_to_litter_gCm2day = output_annual[,,13],
                    roots_to_litter_gCm2day = output[,,14],
                    mean_roots_to_litter_gCm2day = output_mean[,14],
                    mean_annual_roots_to_litter_gCm2day = output_annual[,,14],
                    wood_to_litter_gCm2day = output[,,15],
                    mean_wood_to_litter_gCm2day = output_mean[,15],
                    mean_annual_wood_to_litter_gCm2day = output_annual[,,15],
                    litter_to_som_gCm2day = output[,,16],
                    mean_litter_to_som_gCm2day = output_mean[,16],
                    mean_annual_litter_to_som_gCm2day = output_annual[,,16],
                    woodlitter_to_som_gCm2day = output[,,17],
                    mean_woodlitter_to_som_gCm2day = output_mean[,17],
                    mean_annual_woodlitter_to_som_gCm2day = output_annual[,,17],
                    # Disturbance fluxes
                    FIREemiss_labile_gCm2day = output[,,18],
                    mean_FIREemiss_labile_gCm2day = output_mean[,18],
                    mean_annual_FIREemiss_labile_gCm2day = output_annual[,,18],
                    FIRElitter_labile_gCm2day = output[,,19],
                    mean_FIRElitter_labile_gCm2day = output_mean[,19],
                    mean_annual_FIRElitter_labile_gCm2day = output_annual[,,19],
                    FIREemiss_foliage_gCm2day = output[,,20],
                    mean_FIREemiss_foliage_gCm2day = output_mean[,20],
                    mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,20],
                    FIRElitter_foliage_gCm2day = output[,,21],
                    mean_FIRElitter_foliage_gCm2day = output_mean[,21],
                    mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,21],
                    FIREemiss_roots_gCm2day = output[,,22],
                    mean_FIREemiss_roots_gCm2day = output_mean[,22],
                    mean_annual_FIREemiss_roots_gCm2day = output_annual[,,22],
                    FIRElitter_roots_gCm2day = output[,,23],
                    mean_FIRElitter_roots_gCm2day = output_mean[,23],
                    mean_annual_FIRElitter_roots_gCm2day = output_annual[,,23],
                    FIREemiss_wood_gCm2day = output[,,24],
                    mean_FIREemiss_wood_gCm2day = output_mean[,24],
                    mean_annual_FIREemiss_wood_gCm2day = output_annual[,,24],
                    FIRElitter_wood_gCm2day = output[,,25],
                    mean_FIRElitter_wood_gCm2day = output_mean[,25],
                    mean_annual_FIRElitter_wood_gCm2day = output_annual[,,25],
                    FIREemiss_litter_gCm2day = output[,,26],
                    mean_FIREemiss_litter_gCm2day = output_mean[,26],
                    mean_annual_FIREemiss_litter_gCm2day = output_annual[,,26],
                    FIRElitter_litter_gCm2day = output[,,27],
                    mean_FIRElitter_litter_gCm2day = output_mean[,27],
                    mean_annual_FIRElitter_litter_gCm2day = output_annual[,,27],
                    FIREemiss_woodlitter_gCm2day = output[,,28],
                    mean_FIREemiss_woodlitter_gCm2day = output_mean[,28],
                    mean_annual_FIREemiss_woodlitter_gCm2day = output_annual[,,28],
                    FIRElitter_woodlitter_gCm2day = output[,,29],
                    mean_FIRElitter_woodlitter_gCm2day = output_mean[,29],
                    mean_annual_FIRElitter_woodlitter_gCm2day = output_annual[,,29],
                    FIREemiss_som_gCm2day = output[,,30],
                    mean_FIREemiss_som_gCm2day = output_mean[,30],
                    mean_annual_FIREemiss_som_gCm2day = output_annual[,,30],
                    HARVESTextracted_labile_gCm2day = output[,,31],
                    mean_HARVESTextracted_labile_gCm2day = output_mean[,31],
                    mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,31],
                    HARVESTextracted_foliage_gCm2day = output[,,32],
                    mean_HARVESTextracted_foliage_gCm2day = output_mean[,32],
                    mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,32],
                    HARVESTextracted_roots_gCm2day = output[,,33],
                    mean_HARVESTextracted_roots_gCm2day = output_mean[,33],
                    mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,33],
                    HARVESTextracted_wood_gCm2day = output[,,34],
                    mean_HARVESTextracted_wood_gCm2day = output_mean[,34],
                    mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,34],
                    HARVESTextracted_litter_gCm2day = output[,,35],
                    mean_HARVESTextracted_litter_gCm2day = output_mean[,35],
                    mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,35],
                    HARVESTextracted_woodlitter_gCm2day = output[,,36],
                    mean_HARVESTextracted_woodlitter_gCm2day = output_mean[,36],
                    mean_annual_HARVESTextracted_woodlitter_gCm2day = output_annual[,,36],
                    HARVESTextracted_som_gCm2day = output[,,37],
                    mean_HARVESTextracted_som_gCm2day = output_mean[,37],
                    mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,37],
                    HARVESTlitter_labile_gCm2day = output[,,38],
                    mean_HARVESTlitter_labile_gCm2day = output_mean[,38],
                    mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,38],
                    HARVESTlitter_foliage_gCm2day = output[,,39],
                    mean_HARVESTlitter_foliage_gCm2day = output_mean[,39],
                    mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,39],
                    HARVESTlitter_roots_gCm2day = output[,,40],
                    mean_HARVESTlitter_roots_gCm2day = output_mean[,40],
                    mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,40],
                    HARVESTlitter_wood_gCm2day = output[,,41],
                    mean_HARVESTlitter_wood_gCm2day = output_mean[,41],
                    mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,41],
                    # C pools (gC/m2)
                    labile_gCm2 = output[,,42],
                    mean_labile_gCm2 = output_mean[,42],
                    mean_annual_labile_gCm2 = output_annual[,,42],
                    foliage_gCm2 = output[,,43],
                    mean_foliage_gCm2 = output_mean[,43],
                    mean_annual_foliage_gCm2 = output_annual[,,43],
                    roots_gCm2 = output[,,44],
                    mean_roots_gCm2 = output_mean[,44],
                    mean_annual_roots_gCm2 = output_annual[,,44],
                    wood_gCm2 = output[,,45],
                    mean_wood_gCm2 = output_mean[,45],
                    mean_annual_wood_gCm2 = output_annual[,,45],
                    litter_gCm2 = output[,,46],
                    mean_litter_gCm2 = output_mean[,46],
                    mean_annual_litter_gCm2 = output_annual[,,46],
                    woodlitter_gCm2 = output[,,47],
                    mean_woodlitter_gCm2 = output_mean[,47],
                    mean_annual_woodlitter_gCm2 = output_annual[,,47],
                    som_gCm2 = output[,,48],
                    mean_som_gCm2 = output_mean[,48],
                    mean_annual_som_gCm2 = output_annual[,,48],
                    # Water cycle related
                    ET_kgH2Om2day = output[,,49],
                    mean_ET_kgH2Om2day = output_mean[,49],
                    mean_annual_ET_kgH2Om2day = output_annual[,,49],
                    Etrans_kgH2Om2day = output[,,50],
                    mean_Etrans_kgH2Om2day = output_mean[,50],
                    mean_annual_Etrans_kgH2Om2day = output_annual[,,50],
                    Esoil_kgH2Om2day = output[,,51],
                    mean_Esoil_kgH2Om2day = output_mean[,51],
                    mean_annual_Esoil_kgH2Om2day = output_annual[,,51],
                    Ewetcanopy_kgH2Om2day = output[,,52],
                    mean_Ewetcanopy_kgH2Om2day = output_mean[,52],
                    mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,52],
                    runoff_kgH2Om2day = output[,,53],
                    mean_runoff_kgH2Om2day = output_mean[,53],
                    mean_annual_runoff_kgH2Om2day = output_annual[,,53],
                    underflow_kgH2Om2day = output[,,54],
                    mean_underflow_kgH2Om2day = output_mean[,54],
                    mean_annual_underflow_kgH2Om2day = output_annual[,,54],
                    SurfDrainage_kgH2Om2day = output[,,55],
                    mean_SurfDrainage_kgH2Om2day = output_mean[,55],
                    mean_annual_SurfDrainage_kgH2Om2day = output_annual[,,55],
                    SurfInfiltrated_kgH2Om2day = output[,,56],
                    mean_SurfInfiltrated_kgH2Om2day = output_mean[,56],
                    mean_annual_SurfInfiltrated_kgH2Om2day = output_annual[,,56],
                    Etrans_1st_root_layer_uptake_fraction = output[,,57],
                    mean_Etrans_1st_root_layer_uptake_fraction = output_mean[,57],
                    mean_annual_Etrans_1st_root_layer_uptake_fraction = output_annual[,,57],
                    Etrans_2nd_root_layer_uptake_fraction = output[,,58],
                    mean_Etrans_2nd_root_layer_uptake_fraction = output_mean[,58],
                    mean_annual_Etrans_2nd_root_layer_uptake_fraction = output_annual[,,58],
                    SurfWater_kgH2Om2 = output[,,59],
                    mean_SurfWater_kgH2Om2 = output_mean[,59],
                    mean_annual_SurfWater_kgH2Om2 = output_annual[,,59],
                    wSWP_MPa = output[,,60],
                    mean_wSWP_MPa = output_mean[,60],
                    mean_annual_wSWP_MPa = output_annual[,,60],
                    snow_kgH2Om2 = output[,,61],
                    mean_snow_kgH2Om2 = output_mean[,61],
                    mean_annual_snow_kgH2Om2 = output_annual[,,61],
                    # Canopy (phenology) properties
                    lai_m2m2 = output[,,62],
                    mean_lai_m2m2 = output_mean[,62],
                    mean_annual_lai_m2m2 = output_annual[,,62],
                    gsi = output[,,63],
                    mean_gsi = output_mean[,63],
                    mean_annual_gsi = output_annual[,,63],
                    gsi_itemp = output[,,64],
                    mean_gsi_itemp = output_mean[,64],
                    mean_annual_gsi_itemp = output_annual[,,64],
                    gsi_iphoto = output[,,65],
                    mean_gsi_iphoto = output_mean[,65],
                    mean_annual_gsi_iphoto = output_annual[,,65],
                    gsi_ivpd = output[,,66],
                    mean_gsi_ivpd = output_mean[,66],
                    mean_annual_gsi_ivpd = output_annual[,,66],
                    # Photosynthesis / C~water coupling related
                    gs_demand_supply_ratio = output[,,67],
                    mean_gs_demand_supply_ratio = output_mean[,67],
                    mean_annual_gs_demand_supply_ratio = output_annual[,,67],
                    gs_mmolH2Om2s = output[,,68],
                    mean_gs_mmolH2Om2s = output_mean[,68],
                    mean_annual_gs_mmolH2Om2s = output_annual[,,68],
                    APAR_MJm2day = output[,,69],
                    mean_APAR_MJm2day = output_mean[,69],
                    mean_annual_APAR_MJm2day = output_annual[,,69],
                    gb_mmolH2Om2s = output[,,70],
                    mean_gb_mmolH2Om2s = output_mean[,70],
                    mean_annual_gb_mmolH2Om2s = output_annual[,,70],
                    CiCa = output[,,71],
                    mean_CiCa = output_mean[,71],
                    mean_annual_CiCa = output_annual[,,71],
                    # Misc
                    RootDepth_m = output[,,72],
                    mean_RootDepth_m = output_mean[,72],
                    mean_annual_RootDepth_m = output_annual[,,72],
                    ## Aggregated variables
                    # Mean Transit times
                    MTT_labile_years = MTT_years[,1],
                    MTT_foliage_years = MTT_years[,2],
                    MTT_roots_years = MTT_years[,3],
                    MTT_wood_years = MTT_years[,4],
                    MTT_litter_years = MTT_years[,5],
                    MTT_woodlitter_years = MTT_years[,6],
                    MTT_som_years = MTT_years[,7],
                    # Steady state estimates
                    SS_labile_gCm2 = SS_gCm2[,1],
                    SS_foliage_gCm2 = SS_gCm2[,2],
                    SS_roots_gCm2 = SS_gCm2[,3],
                    SS_wood_gCm2 = SS_gCm2[,4],
                    SS_litter_gCm2 = SS_gCm2[,5],
                    SS_woodlitter_gCm2 = SS_gCm2[,6],
                    SS_som_gCm2 = SS_gCm2[,7])
    # Determine the NPP fraction of expressed NPP
    # i.e. actual growth not GPP-Ra
    NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                         states_all$alloc_roots_gCm2day +
                         states_all$alloc_wood_gCm2day,1,mean)
    NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day,1,mean),
                         apply(states_all$alloc_roots_gCm2day,1,mean),
                         apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
    states_all$NPP_foliage_fraction = NPP_fraction[,1]
    states_all$NPP_roots_fraction = NPP_fraction[,2]
    states_all$NPP_wood_fraction = NPP_fraction[,3]
    # Tidy up variables
    rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P7.R2.023") {
      output_dim = 69 ; MTT_dim = 7 ; SS_dim = 7
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec23",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                              ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                              
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,nos_years=as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                              ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2]))
                              ,residual_waterfrac_in=as.double(c(soil_info[5],soil_info[6],soil_info[6])) 
                              ,porosity_in = as.double(c(soil_info[7],soil_info[8],soil_info[8]))
                              ,pore_size_dist_in = as.double(c(soil_info[9],soil_info[10],soil_info[10]))
                              ,air_entry_in = as.double(c(soil_info[11],soil_info[12],soil_info[12]))
                              ,sat_conductivity_in = as.double(c(soil_info[13],soil_info[14],soil_info[14])))
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))      
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      rhet_woodlitter_gCm2day = output[,,5],
                      mean_rhet_woodlitter_gCm2day = output_mean[,5],
                      mean_annual_rhet_woodlitter_gCm2day = output_annual[,,5],
                      fire_gCm2day = output[,,6],
                      mean_fire_gCm2day = output_mean[,6],
                      mean_annual_fire_gCm2day = output_annual[,,6],
                      harvest_gCm2day = output[,,7],
                      mean_harvest_gCm2day = output_mean[,7],
                      mean_annual_harvest_gCm2day = output_annual[,,7],
                      # Internal fluxes
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      woodlitter_to_som_gCm2day = output[,,16],
                      mean_woodlitter_to_som_gCm2day = output_mean[,16],
                      mean_annual_woodlitter_to_som_gCm2day = output_annual[,,16],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,17],
                      mean_FIREemiss_labile_gCm2day = output_mean[,17],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,17],
                      FIRElitter_labile_gCm2day = output[,,18],
                      mean_FIRElitter_labile_gCm2day = output_mean[,18],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,18],
                      FIREemiss_foliage_gCm2day = output[,,19],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,19],
                      FIRElitter_foliage_gCm2day = output[,,20],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,20],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,20],
                      FIREemiss_roots_gCm2day = output[,,21],
                      mean_FIREemiss_roots_gCm2day = output_mean[,21],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,21],
                      FIRElitter_roots_gCm2day = output[,,22],
                      mean_FIRElitter_roots_gCm2day = output_mean[,22],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,22],
                      FIREemiss_wood_gCm2day = output[,,23],
                      mean_FIREemiss_wood_gCm2day = output_mean[,23],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,23],
                      FIRElitter_wood_gCm2day = output[,,24],
                      mean_FIRElitter_wood_gCm2day = output_mean[,24],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,24],
                      FIREemiss_litter_gCm2day = output[,,25],
                      mean_FIREemiss_litter_gCm2day = output_mean[,25],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,25],
                      FIRElitter_litter_gCm2day = output[,,26],
                      mean_FIRElitter_litter_gCm2day = output_mean[,26],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,26],
                      FIREemiss_woodlitter_gCm2day = output[,,27],
                      mean_FIREemiss_woodlitter_gCm2day = output_mean[,27],
                      mean_annual_FIREemiss_woodlitter_gCm2day = output_annual[,,27],
                      FIRElitter_woodlitter_gCm2day = output[,,28],
                      mean_FIRElitter_woodlitter_gCm2day = output_mean[,28],
                      mean_annual_FIRElitter_woodlitter_gCm2day = output_annual[,,28],
                      FIREemiss_som_gCm2day = output[,,29],
                      mean_FIREemiss_som_gCm2day = output_mean[,29],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,29],
                      HARVESTextracted_labile_gCm2day = output[,,30],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,30],
                      HARVESTextracted_foliage_gCm2day = output[,,31],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,31],
                      HARVESTextracted_roots_gCm2day = output[,,32],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,32],
                      HARVESTextracted_wood_gCm2day = output[,,33],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,33],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,33],
                      HARVESTextracted_litter_gCm2day = output[,,34],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,34],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,34],
                      HARVESTextracted_woodlitter_gCm2day = output[,,35],
                      mean_HARVESTextracted_woodlitter_gCm2day = output_mean[,35],
                      mean_annual_HARVESTextracted_woodlitter_gCm2day = output_annual[,,35],
                      HARVESTextracted_som_gCm2day = output[,,36],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,36],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,36],
                      HARVESTlitter_labile_gCm2day = output[,,37],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,37],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,37],
                      HARVESTlitter_foliage_gCm2day = output[,,38],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,38],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,38],
                      HARVESTlitter_roots_gCm2day = output[,,39],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,39],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,39],
                      HARVESTlitter_wood_gCm2day = output[,,40],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,40],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,40],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,41],
                      mean_labile_gCm2 = output_mean[,41],
                      mean_annual_labile_gCm2 = output_annual[,,41],
                      foliage_gCm2 = output[,,42],
                      mean_foliage_gCm2 = output_mean[,42],
                      mean_annual_foliage_gCm2 = output_annual[,,42],
                      roots_gCm2 = output[,,43],
                      mean_roots_gCm2 = output_mean[,43],
                      mean_annual_roots_gCm2 = output_annual[,,43],
                      wood_gCm2 = output[,,44],
                      mean_wood_gCm2 = output_mean[,44],
                      mean_annual_wood_gCm2 = output_annual[,,44],
                      litter_gCm2 = output[,,45],
                      mean_litter_gCm2 = output_mean[,45],
                      mean_annual_litter_gCm2 = output_annual[,,45],
                      woodlitter_gCm2 = output[,,46],
                      mean_woodlitter_gCm2 = output_mean[,46],
                      mean_annual_woodlitter_gCm2 = output_annual[,,46],
                      som_gCm2 = output[,,47],
                      mean_som_gCm2 = output_mean[,47],
                      mean_annual_som_gCm2 = output_annual[,,47],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,48],
                      mean_ET_kgH2Om2day = output_mean[,48],
                      mean_annual_ET_kgH2Om2day = output_annual[,,48],
                      Etrans_kgH2Om2day = output[,,49],
                      mean_Etrans_kgH2Om2day = output_mean[,49],
                      mean_annual_Etrans_kgH2Om2day = output_annual[,,49],
                      Esoil_kgH2Om2day = output[,,50],
                      mean_Esoil_kgH2Om2day = output_mean[,50],
                      mean_annual_Esoil_kgH2Om2day = output_annual[,,50],
                      Ewetcanopy_kgH2Om2day = output[,,51],
                      mean_Ewetcanopy_kgH2Om2day = output_mean[,51],
                      mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,51],
                      runoff_kgH2Om2day = output[,,52],
                      mean_runoff_kgH2Om2day = output_mean[,52],
                      mean_annual_runoff_kgH2Om2day = output_annual[,,52],
                      underflow_kgH2Om2day = output[,,53],
                      mean_underflow_kgH2Om2day = output_mean[,53],
                      mean_annual_underflow_kgH2Om2day = output_annual[,,53],
                      SurfDrainage_kgH2Om2day = output[,,54],
                      mean_SurfDrainage_kgH2Om2day = output_mean[,54],
                      mean_annual_SurfDrainage_kgH2Om2day = output_annual[,,54],
                      SurfInfiltrated_kgH2Om2day = output[,,55],
                      mean_SurfInfiltrated_kgH2Om2day = output_mean[,55],
                      mean_annual_SurfInfiltrated_kgH2Om2day = output_annual[,,55],
                      Etrans_1st_root_layer_uptake_fraction = output[,,56],
                      mean_Etrans_1st_root_layer_uptake_fraction = output_mean[,56],
                      mean_annual_Etrans_1st_root_layer_uptake_fraction = output_annual[,,56],
                      Etrans_2nd_root_layer_uptake_fraction = output[,,57],
                      mean_Etrans_2nd_root_layer_uptake_fraction = output_mean[,57],
                      mean_annual_Etrans_2nd_root_layer_uptake_fraction = output_annual[,,57],
                      SurfWater_kgH2Om2 = output[,,58],
                      mean_SurfWater_kgH2Om2 = output_mean[,58],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,58],
                      wSWP_MPa = output[,,59],
                      mean_wSWP_MPa = output_mean[,59],
                      mean_annual_wSWP_MPa = output_annual[,,59],
                      snow_kgH2Om2 = output[,,60],
                      mean_snow_kgH2Om2 = output_mean[,60],
                      mean_annual_snow_kgH2Om2 = output_annual[,,60],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,61],
                      mean_lai_m2m2 = output_mean[,56],
                      mean_annual_lai_m2m2 = output_annual[,,61],
                      cgi = output[,,62],
                      mean_cgi = output_mean[,62],
                      mean_annual_cgi = output_annual[,,62],
                      ncce_gCm2day = output[,,63],
                      mean_ncce_gCm2day = output_mean[,63],
                      mean_annual_ncce_gCm2day = output_annual[,,63],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,64],
                      mean_gs_demand_supply_ratio = output_mean[,64],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,64],
                      gs_mmolH2Om2s = output[,,65],
                      mean_gs_mmolH2Om2s = output_mean[,65],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,65],
                      APAR_MJm2day = output[,,66],
                      mean_APAR_MJm2day = output_mean[,66],
                      mean_annual_APAR_MJm2day = output_annual[,,66],
                      gb_mmolH2Om2s = output[,,67],
                      mean_gb_mmolH2Om2s = output_mean[,67],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,67],
                      CiCa = output[,,68],
                      mean_CiCa = output_mean[,68],
                      mean_annual_CiCa = output_annual[,,68],
                      # Misc
                      RootDepth_m = output[,,69],
                      mean_RootDepth_m = output_mean[,69],
                      mean_annual_RootDepth_m = output_annual[,,69],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_woodlitter_years = MTT_years[,6],
                      MTT_som_years = MTT_years[,7],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_woodlitter_gCm2 = SS_gCm2[,6],
                      SS_som_gCm2 = SS_gCm2[,7])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P8.R2.024") {
      output_dim = 70 ; MTT_dim = 7 ; SS_dim = 7
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec24",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                              ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                              
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,nos_years=as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                              ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2]))
                              ,residual_waterfrac_in=as.double(c(soil_info[5],soil_info[6],soil_info[6])) 
                              ,porosity_in = as.double(c(soil_info[7],soil_info[8],soil_info[8]))
                              ,pore_size_dist_in = as.double(c(soil_info[9],soil_info[10],soil_info[10]))
                              ,air_entry_in = as.double(c(soil_info[11],soil_info[12],soil_info[12]))
                              ,sat_conductivity_in = as.double(c(soil_info[13],soil_info[14],soil_info[14])))
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      rhet_woodlitter_gCm2day = output[,,5],
                      mean_rhet_woodlitter_gCm2day = output_mean[,5],
                      mean_annual_rhet_woodlitter_gCm2day = output_annual[,,5],
                      fire_gCm2day = output[,,6],
                      mean_fire_gCm2day = output_mean[,6],
                      mean_annual_fire_gCm2day = output_annual[,,6],
                      harvest_gCm2day = output[,,7],
                      mean_harvest_gCm2day = output_mean[,7],
                      mean_annual_harvest_gCm2day = output_annual[,,7],
                      # Internal fluxes
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      woodlitter_to_som_gCm2day = output[,,16],
                      mean_woodlitter_to_som_gCm2day = output_mean[,16],
                      mean_annual_woodlitter_to_som_gCm2day = output_annual[,,16],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,17],
                      mean_FIREemiss_labile_gCm2day = output_mean[,17],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,17],
                      FIRElitter_labile_gCm2day = output[,,18],
                      mean_FIRElitter_labile_gCm2day = output_mean[,18],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,18],
                      FIREemiss_foliage_gCm2day = output[,,19],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,19],
                      FIRElitter_foliage_gCm2day = output[,,20],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,20],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,20],
                      FIREemiss_roots_gCm2day = output[,,21],
                      mean_FIREemiss_roots_gCm2day = output_mean[,21],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,21],
                      FIRElitter_roots_gCm2day = output[,,22],
                      mean_FIRElitter_roots_gCm2day = output_mean[,22],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,22],
                      FIREemiss_wood_gCm2day = output[,,23],
                      mean_FIREemiss_wood_gCm2day = output_mean[,23],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,23],
                      FIRElitter_wood_gCm2day = output[,,24],
                      mean_FIRElitter_wood_gCm2day = output_mean[,24],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,24],
                      FIREemiss_litter_gCm2day = output[,,25],
                      mean_FIREemiss_litter_gCm2day = output_mean[,25],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,25],
                      FIRElitter_litter_gCm2day = output[,,26],
                      mean_FIRElitter_litter_gCm2day = output_mean[,26],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,26],
                      FIREemiss_woodlitter_gCm2day = output[,,27],
                      mean_FIREemiss_woodlitter_gCm2day = output_mean[,27],
                      mean_annual_FIREemiss_woodlitter_gCm2day = output_annual[,,27],
                      FIRElitter_woodlitter_gCm2day = output[,,28],
                      mean_FIRElitter_woodlitter_gCm2day = output_mean[,28],
                      mean_annual_FIRElitter_woodlitter_gCm2day = output_annual[,,28],
                      FIREemiss_som_gCm2day = output[,,29],
                      mean_FIREemiss_som_gCm2day = output_mean[,29],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,29],
                      HARVESTextracted_labile_gCm2day = output[,,30],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,30],
                      HARVESTextracted_foliage_gCm2day = output[,,31],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,31],
                      HARVESTextracted_roots_gCm2day = output[,,32],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,32],
                      HARVESTextracted_wood_gCm2day = output[,,33],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,33],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,33],
                      HARVESTextracted_litter_gCm2day = output[,,34],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,34],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,34],
                      HARVESTextracted_woodlitter_gCm2day = output[,,35],
                      mean_HARVESTextracted_woodlitter_gCm2day = output_mean[,35],
                      mean_annual_HARVESTextracted_woodlitter_gCm2day = output_annual[,,35],
                      HARVESTextracted_som_gCm2day = output[,,36],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,36],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,36],
                      HARVESTlitter_labile_gCm2day = output[,,37],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,37],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,37],
                      HARVESTlitter_foliage_gCm2day = output[,,38],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,38],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,38],
                      HARVESTlitter_roots_gCm2day = output[,,39],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,39],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,39],
                      HARVESTlitter_wood_gCm2day = output[,,40],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,40],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,40],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,41],
                      mean_labile_gCm2 = output_mean[,41],
                      mean_annual_labile_gCm2 = output_annual[,,41],
                      foliage_gCm2 = output[,,42],
                      mean_foliage_gCm2 = output_mean[,42],
                      mean_annual_foliage_gCm2 = output_annual[,,42],
                      roots_gCm2 = output[,,43],
                      mean_roots_gCm2 = output_mean[,43],
                      mean_annual_roots_gCm2 = output_annual[,,43],
                      wood_gCm2 = output[,,44],
                      mean_wood_gCm2 = output_mean[,44],
                      mean_annual_wood_gCm2 = output_annual[,,44],
                      litter_gCm2 = output[,,45],
                      mean_litter_gCm2 = output_mean[,45],
                      mean_annual_litter_gCm2 = output_annual[,,45],
                      woodlitter_gCm2 = output[,,46],
                      mean_woodlitter_gCm2 = output_mean[,46],
                      mean_annual_woodlitter_gCm2 = output_annual[,,46],
                      som_gCm2 = output[,,47],
                      mean_som_gCm2 = output_mean[,47],
                      mean_annual_som_gCm2 = output_annual[,,47],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,48],
                      mean_ET_kgH2Om2day = output_mean[,48],
                      mean_annual_ET_kgH2Om2day = output_annual[,,48],
                      Etrans_kgH2Om2day = output[,,49],
                      mean_Etrans_kgH2Om2day = output_mean[,49],
                      mean_annual_Etrans_kgH2Om2day = output_annual[,,49],
                      Esoil_kgH2Om2day = output[,,50],
                      mean_Esoil_kgH2Om2day = output_mean[,50],
                      mean_annual_Esoil_kgH2Om2day = output_annual[,,50],
                      Ewetcanopy_kgH2Om2day = output[,,51],
                      mean_Ewetcanopy_kgH2Om2day = output_mean[,51],
                      mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,51],
                      runoff_kgH2Om2day = output[,,52],
                      mean_runoff_kgH2Om2day = output_mean[,52],
                      mean_annual_runoff_kgH2Om2day = output_annual[,,52],
                      underflow_kgH2Om2day = output[,,53],
                      mean_underflow_kgH2Om2day = output_mean[,53],
                      mean_annual_underflow_kgH2Om2day = output_annual[,,53],
                      SurfDrainage_kgH2Om2day = output[,,54],
                      mean_SurfDrainage_kgH2Om2day = output_mean[,54],
                      mean_annual_SurfDrainage_kgH2Om2day = output_annual[,,54],
                      SurfInfiltrated_kgH2Om2day = output[,,55],
                      mean_SurfInfiltrated_kgH2Om2day = output_mean[,55],
                      mean_annual_SurfInfiltrated_kgH2Om2day = output_annual[,,55],
                      Etrans_1st_root_layer_uptake_fraction = output[,,56],
                      mean_Etrans_1st_root_layer_uptake_fraction = output_mean[,56],
                      mean_annual_Etrans_1st_root_layer_uptake_fraction = output_annual[,,56],
                      Etrans_2nd_root_layer_uptake_fraction = output[,,57],
                      mean_Etrans_2nd_root_layer_uptake_fraction = output_mean[,57],
                      mean_annual_Etrans_2nd_root_layer_uptake_fraction = output_annual[,,57],
                      SurfWater_kgH2Om2 = output[,,58],
                      mean_SurfWater_kgH2Om2 = output_mean[,58],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,58],
                      wSWP_MPa = output[,,59],
                      mean_wSWP_MPa = output_mean[,59],
                      mean_annual_wSWP_MPa = output_annual[,,59],
                      snow_kgH2Om2 = output[,,60],
                      mean_snow_kgH2Om2 = output_mean[,60],
                      mean_annual_snow_kgH2Om2 = output_annual[,,60],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,61],
                      mean_lai_m2m2 = output_mean[,61],
                      mean_annual_lai_m2m2 = output_annual[,,61],
                      cgi = output[,,62],
                      mean_cgi = output_mean[,62],
                      mean_annual_cgi = output_annual[,,62],
                      cmi = output[,,63],
                      mean_cmi = output_mean[,63],
                      mean_annual_cmi = output_annual[,,63],
                      ncce_gCm2day = output[,,64],
                      mean_ncce_gCm2day = output_mean[,64],
                      mean_annual_ncce_gCm2day = output_annual[,,64],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,65],
                      mean_gs_demand_supply_ratio = output_mean[,65],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,65],
                      gs_mmolH2Om2s = output[,,66],
                      mean_gs_mmolH2Om2s = output_mean[,66],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,66],
                      APAR_MJm2day = output[,,67],
                      mean_APAR_MJm2day = output_mean[,67],
                      mean_annual_APAR_MJm2day = output_annual[,,67],
                      gb_mmolH2Om2s = output[,,68],
                      mean_gb_mmolH2Om2s = output[,,68],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,68],
                      CiCa = output[,,69],
                      mean_CiCa = output_mean[,69],
                      mean_annual_CiCa = output_annual[,,69],
                      # Misc
                      RootDepth_m = output[,,70],
                      mean_RootDepth_m = output_mean[,70],
                      mean_annual_RootDepth_m = output_annual[,,70],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_woodlitter_years = MTT_years[,6],
                      MTT_som_years = MTT_years[,7],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_woodlitter_gCm2 = SS_gCm2[,6],
                      SS_som_gCm2 = SS_gCm2[,7])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P9.R2.025") {
      output_dim = 70 ; MTT_dim = 7 ; SS_dim = 7
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec25",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                              ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                              
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,nos_years=as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                              ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2]))
                              ,residual_waterfrac_in=as.double(c(soil_info[5],soil_info[6],soil_info[6])) 
                              ,porosity_in = as.double(c(soil_info[7],soil_info[8],soil_info[8]))
                              ,pore_size_dist_in = as.double(c(soil_info[9],soil_info[10],soil_info[10]))
                              ,air_entry_in = as.double(c(soil_info[11],soil_info[12],soil_info[12]))
                              ,sat_conductivity_in = as.double(c(soil_info[13],soil_info[14],soil_info[14])))
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      rhet_woodlitter_gCm2day = output[,,5],
                      mean_rhet_woodlitter_gCm2day = output_mean[,5],
                      mean_annual_rhet_woodlitter_gCm2day = output_annual[,,5],
                      fire_gCm2day = output[,,6],
                      mean_fire_gCm2day = output_mean[,6],
                      mean_annual_fire_gCm2day = output_annual[,,6],
                      harvest_gCm2day = output[,,7],
                      mean_harvest_gCm2day = output_mean[,7],
                      mean_annual_harvest_gCm2day = output_annual[,,7],
                      # Internal fluxes
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      woodlitter_to_som_gCm2day = output[,,16],
                      mean_woodlitter_to_som_gCm2day = output_mean[,16],
                      mean_annual_woodlitter_to_som_gCm2day = output_annual[,,16],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,17],
                      mean_FIREemiss_labile_gCm2day = output_mean[,17],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,17],
                      FIRElitter_labile_gCm2day = output[,,18],
                      mean_FIRElitter_labile_gCm2day = output_mean[,18],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,18],
                      FIREemiss_foliage_gCm2day = output[,,19],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,19],
                      FIRElitter_foliage_gCm2day = output[,,20],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,20],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,20],
                      FIREemiss_roots_gCm2day = output[,,21],
                      mean_FIREemiss_roots_gCm2day = output_mean[,21],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,21],
                      FIRElitter_roots_gCm2day = output[,,22],
                      mean_FIRElitter_roots_gCm2day = output_mean[,22],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,22],
                      FIREemiss_wood_gCm2day = output[,,23],
                      mean_FIREemiss_wood_gCm2day = output_mean[,23],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,23],
                      FIRElitter_wood_gCm2day = output[,,24],
                      mean_FIRElitter_wood_gCm2day = output_mean[,24],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,24],
                      FIREemiss_litter_gCm2day = output[,,25],
                      mean_FIREemiss_litter_gCm2day = output_mean[,25],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,25],
                      FIRElitter_litter_gCm2day = output[,,26],
                      mean_FIRElitter_litter_gCm2day = output_mean[,26],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,26],
                      FIREemiss_woodlitter_gCm2day = output[,,27],
                      mean_FIREemiss_woodlitter_gCm2day = output_mean[,27],
                      mean_annual_FIREemiss_woodlitter_gCm2day = output_annual[,,27],
                      FIRElitter_woodlitter_gCm2day = output[,,28],
                      mean_FIRElitter_woodlitter_gCm2day = output_mean[,28],
                      mean_annual_FIRElitter_woodlitter_gCm2day = output_annual[,,28],
                      FIREemiss_som_gCm2day = output[,,29],
                      mean_FIREemiss_som_gCm2day = output_mean[,29],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,29],
                      HARVESTextracted_labile_gCm2day = output[,,30],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,30],
                      HARVESTextracted_foliage_gCm2day = output[,,31],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,31],
                      HARVESTextracted_roots_gCm2day = output[,,32],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,32],
                      HARVESTextracted_wood_gCm2day = output[,,33],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,33],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,33],
                      HARVESTextracted_litter_gCm2day = output[,,34],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,34],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,34],
                      HARVESTextracted_woodlitter_gCm2day = output[,,35],
                      mean_HARVESTextracted_woodlitter_gCm2day = output_mean[,35],
                      mean_annual_HARVESTextracted_woodlitter_gCm2day = output_annual[,,35],
                      HARVESTextracted_som_gCm2day = output[,,36],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,36],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,36],
                      HARVESTlitter_labile_gCm2day = output[,,37],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,37],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,37],
                      HARVESTlitter_foliage_gCm2day = output[,,38],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,38],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,38],
                      HARVESTlitter_roots_gCm2day = output[,,39],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,39],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,39],
                      HARVESTlitter_wood_gCm2day = output[,,40],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,40],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,40],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,41],
                      mean_labile_gCm2 = output_mean[,41],
                      mean_annual_labile_gCm2 = output_annual[,,41],
                      foliage_gCm2 = output[,,42],
                      mean_foliage_gCm2 = output_mean[,42],
                      mean_annual_foliage_gCm2 = output_annual[,,42],
                      roots_gCm2 = output[,,43],
                      mean_roots_gCm2 = output_mean[,43],
                      mean_annual_roots_gCm2 = output_annual[,,43],
                      wood_gCm2 = output[,,44],
                      mean_wood_gCm2 = output_mean[,44],
                      mean_annual_wood_gCm2 = output_annual[,,44],
                      litter_gCm2 = output[,,45],
                      mean_litter_gCm2 = output_mean[,45],
                      mean_annual_litter_gCm2 = output_annual[,,45],
                      woodlitter_gCm2 = output[,,46],
                      mean_woodlitter_gCm2 = output_mean[,46],
                      mean_annual_woodlitter_gCm2 = output_annual[,,46],
                      som_gCm2 = output[,,47],
                      mean_som_gCm2 = output_mean[,47],
                      mean_annual_som_gCm2 = output_annual[,,47],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,48],
                      mean_ET_kgH2Om2day = output_mean[,48],
                      mean_annual_ET_kgH2Om2day = output_annual[,,48],
                      Etrans_kgH2Om2day = output[,,49],
                      mean_Etrans_kgH2Om2day = output_mean[,49],
                      mean_annual_Etrans_kgH2Om2day = output_annual[,,49],
                      Esoil_kgH2Om2day = output[,,50],
                      mean_Esoil_kgH2Om2day = output_mean[,50],
                      mean_annual_Esoil_kgH2Om2day = output_annual[,,50],
                      Ewetcanopy_kgH2Om2day = output[,,51],
                      mean_Ewetcanopy_kgH2Om2day = output_mean[,51],
                      mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,51],
                      runoff_kgH2Om2day = output[,,52],
                      mean_runoff_kgH2Om2day = output_mean[,52],
                      mean_annual_runoff_kgH2Om2day = output_annual[,,52],
                      underflow_kgH2Om2day = output[,,53],
                      mean_underflow_kgH2Om2day = output_mean[,53],
                      mean_annual_underflow_kgH2Om2day = output_annual[,,53],
                      SurfDrainage_kgH2Om2day = output[,,54],
                      mean_SurfDrainage_kgH2Om2day = output_mean[,54],
                      mean_annual_SurfDrainage_kgH2Om2day = output_annual[,,54],
                      SurfInfiltrated_kgH2Om2day = output[,,55],
                      mean_SurfInfiltrated_kgH2Om2day = output_mean[,55],
                      mean_annual_SurfInfiltrated_kgH2Om2day = output_annual[,,55],
                      Etrans_1st_root_layer_uptake_fraction = output[,,56],
                      mean_Etrans_1st_root_layer_uptake_fraction = output_mean[,56],
                      mean_annual_Etrans_1st_root_layer_uptake_fraction = output_annual[,,56],
                      Etrans_2nd_root_layer_uptake_fraction = output[,,57],
                      mean_Etrans_2nd_root_layer_uptake_fraction = output_mean[,57],
                      mean_annual_Etrans_2nd_root_layer_uptake_fraction = output_annual[,,57],
                      SurfWater_kgH2Om2 = output[,,58],
                      mean_SurfWater_kgH2Om2 = output_mean[,58],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,58],
                      wSWP_MPa = output[,,59],
                      mean_wSWP_MPa = output_mean[,59],
                      mean_annual_wSWP_MPa = output_annual[,,59],
                      snow_kgH2Om2 = output[,,60],
                      mean_snow_kgH2Om2 = output_mean[,60],
                      mean_annual_snow_kgH2Om2 = output_annual[,,60],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,61],
                      mean_lai_m2m2 = output_mean[,61],
                      mean_annual_lai_m2m2 = output_annual[,,61],
                      cgi = output[,,62],
                      mean_cgi = output_mean[,62],
                      mean_annual_cgi = output_annual[,,62],
                      cmi = output[,,63],
                      mean_cmi = output_mean[,63],
                      mean_annual_cmi = output_annual[,,63],
                      ncce_gCm2day = output[,,64],
                      mean_ncce_gCm2day = output_mean[,64],
                      mean_annual_ncce_gCm2day = output_annual[,,64],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,65],
                      mean_gs_demand_supply_ratio = output_mean[,65],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,65],
                      gs_mmolH2Om2s = output[,,66],
                      mean_gs_mmolH2Om2s = output_mean[,66],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,66],
                      APAR_MJm2day = output[,,67],
                      mean_APAR_MJm2day = output_mean[,67],
                      mean_annual_APAR_MJm2day = output_annual[,,67],
                      gb_mmolH2Om2s = output[,,68],
                      mean_gb_mmolH2Om2s = output[,,68],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,68],
                      CiCa = output[,,69],
                      mean_CiCa = output_mean[,69],
                      mean_annual_CiCa = output_annual[,,69],
                      # Misc
                      RootDepth_m = output[,,70],
                      mean_RootDepth_m = output_mean[,70],
                      mean_annual_RootDepth_m = output_annual[,,70],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_woodlitter_years = MTT_years[,6],
                      MTT_som_years = MTT_years[,7],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_woodlitter_gCm2 = SS_gCm2[,6],
                      SS_som_gCm2 = SS_gCm2[,7])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P10.R2.026") {
# Interface needs updating
      output_dim=28
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      crop_file_location=paste(PROJECT$exepath,"winter_wheat_development.csv", sep="")
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdalec26",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met)),pars=as.double(pars_in)
                              ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var4=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,out_var5=as.double(array(0,dim=c(nos_iter,MTT_dim,noyears)))
                              ,out_var6=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site])
                              ,nomet=as.integer(dim(met)[2]),nofluxes=as.integer(PROJECT$model$nofluxes[site])
                              ,nopools=as.integer(PROJECT$model$nopools[site]),pft=as.integer(pft)
                              ,pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                              ,nos_years=as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                              ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3)))
                              ,residual_waterfrac_in=as.double(c(soil_info[5],soil_info[6],soil_info[6])) 
                              ,porosity_in = as.double(c(soil_info[7],soil_info[8],soil_info[8]))
                              ,pore_size_dist_in = as.double(c(soil_info[9],soil_info[10],soil_info[10]))
                              ,air_entry_in = as.double(c(soil_info[11],soil_info[12],soil_info[12]))
                              ,sat_conductivity_in = as.double(c(soil_info[13],soil_info[14],soil_info[14]))
                              ,exepath=as.character(crop_file_location),pathlength=as.integer(nchar(crop_file_location)))
      output = tmp$out_var   ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))
      aNPP = tmp$out_var2    ; aNPP = array(aNPP, dim=c(nos_iter,aNPP_dim))
      MTT = tmp$out_var3     ; MTT = array(MTT, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var4 ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      aMTT = tmp$out_var5    ; aMTT = array(aMTT, dim=c(nos_iter,MTT_dim,noyears))
      MTTnat = tmp$out_var6  ; MTTnat = array(MTTnat, dim=c(nos_iter,MTT_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], litter_gCm2 = output[,,10],
                      labile_gCm2 = output[,,11], foliage_gCm2 = output[,,12],
                      harvest_gCm2day = output[,,13], gsi = output[,,14],
                      gsi_itemp = output[,,15], gsi_iphoto = output[,,16],
                      gsi_ivpd = output[,,17], ET_kgH2Om2day = output[,,18],
                      SurfWater_kgH2Om2 = output[,,19], wSWP_MPa = output[,,20],
                      woodlitter_gCm2 = output[,,21], fire_gCm2day = output[,,23],
                      gs_demand_supply = output[,,24], gs_total_canopy = output[,,25],
                      APAR_MJm2day = output[,,26], gb_total_canopy = output[,,27],
                      CiCa = output[,,28],
                      aNPP = aNPP, MTT = MTT, SS_gCm2 = SS_gCm2, aMTT = aMTT, natMTT = MTTnat)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
      # Final tidy
      rm(output,MTT_gCm2,SS_gCm2)
  } else if (model_name == "DALEC.A1.C2.D2.F2.H1.P4.R2.010") {
    output_dim = 59 ; MTT_dim = 7 ; SS_dim = 7
    # Load the required dalec shared object
    dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
    tmp=.Fortran( "rdalec10",output_dim=as.integer(output_dim)
                            ,aNPP_dim=as.integer(aNPP_dim)
                            ,MTT_dim=as.integer(MTT_dim)
                            ,SS_dim = as.integer(SS_dim)
                            ,fire_dim=as.integer(fire_dim)
                            ,met=as.double(t(met))
                            ,pars=as.double(pars_in)
                            ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                            ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                            ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                            ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                            ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                           
                            ,lat=as.double(lat)
                            ,nopars=as.integer(PROJECT$model$nopars[site])
                            ,nomet=as.integer(dim(met)[2])
                            ,nofluxes=as.integer(PROJECT$model$nofluxes[site])
                            ,nopools=as.integer(PROJECT$model$nopools[site])
                            ,pft=as.integer(pft)
                            ,nodays=as.integer(dim(met)[1])
                            ,nos_years=as.integer(noyears)
                            ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1]))))
                            ,nos_iter=as.integer(nos_iter))
    # Extract the different output variables
    output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
    output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
    output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))
    MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
    SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
    # Unload the current dalec shared object
    dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
    rm(tmp) ; gc()
    # create output object
    states_all=list(# Ecosystem fluxes
                    gpp_gCm2day = output[,,1],
                    mean_gpp_gCm2day = output_mean[,1],
                    mean_annual_gpp_gCm2day = output_annual[,,1],
                    rauto_gCm2day = output[,,2],
                    mean_rauto_gCm2day = output_mean[,2],
                    mean_annual_rauto_gCm2day = output_annual[,,2],
                    rg_foliage_gCm2day = output[,,3],
                    mean_rg_foliage_gCm2day = output_mean[,3],
                    mean_annual_rg_foliage_gCm2day = output_annual[,,3],
                    rhet_litter_gCm2day = output[,,4],
                    mean_rhet_litter_gCm2day = output_mean[,4],
                    mean_annual_rhet_litter_gCm2day = output_annual[,,4],
                    rhet_som_gCm2day = output[,,5],
                    mean_rhet_som_gCm2day = output_mean[,5],
                    mean_annual_rhet_som_gCm2day = output_annual[,,5],
                    rhet_woodlitter_gCm2day = output[,,6],
                    mean_rhet_woodlitter_gCm2day = output_mean[,6],
                    mean_annual_rhet_woodlitter_gCm2day = output_annual[,,6],
                    fire_gCm2day = output[,,7],
                    mean_fire_gCm2day = output_mean[,7],
                    mean_annual_fire_gCm2day = output_annual[,,7],
                    harvest_gCm2day = output[,,8],
                    mean_harvest_gCm2day = output_mean[,8],
                    mean_annual_harvest_gCm2day = output_annual[,,8],
                    # Internal fluxes
                    alloc_labile_gCm2day = output[,,9],
                    mean_alloc_labile_gCm2day = output_mean[,9],
                    mean_annual_alloc_labile_gCm2day = output_annual[,,9],
                    alloc_roots_gCm2day = output[,,10],
                    mean_alloc_roots_gCm2day = output_mean[,10],
                    mean_annual_alloc_roots_gCm2day = output_annual[,,10],
                    alloc_wood_gCm2day = output[,,11],
                    mean_alloc_wood_gCm2day = output_mean[,11],
                    mean_annual_alloc_wood_gCm2day = output_annual[,,11],
                    labile_to_foliage_gCm2day = output[,,12],
                    mean_labile_to_foliage_gCm2day = output_mean[,12],
                    mean_annual_labile_to_foliage_gCm2day = output_annual[,,12],
                    foliage_to_litter_gCm2day = output[,,13],
                    mean_foliage_to_litter_gCm2day = output_mean[,13],
                    mean_annual_foliage_to_litter_gCm2day = output_annual[,,13],
                    roots_to_litter_gCm2day = output[,,14],
                    mean_roots_to_litter_gCm2day = output_mean[,14],
                    mean_annual_roots_to_litter_gCm2day = output_annual[,,14],
                    wood_to_litter_gCm2day = output[,,15],
                    mean_wood_to_litter_gCm2day = output_mean[,15],
                    mean_annual_wood_to_litter_gCm2day = output_annual[,,15],
                    litter_to_som_gCm2day = output[,,16],
                    mean_litter_to_som_gCm2day = output_mean[,16],
                    mean_annual_litter_to_som_gCm2day = output_annual[,,16],
                    woodlitter_to_som_gCm2day = output[,,17],
                    mean_woodlitter_to_som_gCm2day = output_mean[,17],
                    mean_annual_woodlitter_to_som_gCm2day = output_annual[,,17],
                    # Disturbance fluxes
                    FIREemiss_labile_gCm2day = output[,,18],
                    mean_FIREemiss_labile_gCm2day = output_mean[,18],
                    mean_annual_FIREemiss_labile_gCm2day = output_annual[,,18],
                    FIRElitter_labile_gCm2day = output[,,19],
                    mean_FIRElitter_labile_gCm2day = output_mean[,19],
                    mean_annual_FIRElitter_labile_gCm2day = output_annual[,,19],
                    FIREemiss_foliage_gCm2day = output[,,20],
                    mean_FIREemiss_foliage_gCm2day = output_mean[,20],
                    mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,20],
                    FIRElitter_foliage_gCm2day = output[,,21],
                    mean_FIRElitter_foliage_gCm2day = output_mean[,21],
                    mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,21],
                    FIREemiss_roots_gCm2day = output[,,22],
                    mean_FIREemiss_roots_gCm2day = output_mean[,22],
                    mean_annual_FIREemiss_roots_gCm2day = output_annual[,,22],
                    FIRElitter_roots_gCm2day = output[,,23],
                    mean_FIRElitter_roots_gCm2day = output_mean[,23],
                    mean_annual_FIRElitter_roots_gCm2day = output_annual[,,23],
                    FIREemiss_wood_gCm2day = output[,,24],
                    mean_FIREemiss_wood_gCm2day = output_mean[,24],
                    mean_annual_FIREemiss_wood_gCm2day = output_annual[,,24],
                    FIRElitter_wood_gCm2day = output[,,25],
                    mean_FIRElitter_wood_gCm2day = output_mean[,25],
                    mean_annual_FIRElitter_wood_gCm2day = output_annual[,,25],
                    FIREemiss_litter_gCm2day = output[,,26],
                    mean_FIREemiss_litter_gCm2day = output_mean[,26],
                    mean_annual_FIREemiss_litter_gCm2day = output_annual[,,26],
                    FIRElitter_litter_gCm2day = output[,,27],
                    mean_FIRElitter_litter_gCm2day = output_mean[,27],
                    mean_annual_FIRElitter_litter_gCm2day = output_annual[,,27],
                    FIREemiss_woodlitter_gCm2day = output[,,28],
                    mean_FIREemiss_woodlitter_gCm2day = output_mean[,28],
                    mean_annual_FIREemiss_woodlitter_gCm2day = output_annual[,,28],
                    FIRElitter_woodlitter_gCm2day = output[,,29],
                    mean_FIRElitter_woodlitter_gCm2day = output_mean[,29],
                    mean_annual_FIRElitter_woodlitter_gCm2day = output_annual[,,29],
                    FIREemiss_som_gCm2day = output[,,30],
                    mean_FIREemiss_som_gCm2day = output_mean[,30],
                    mean_annual_FIREemiss_som_gCm2day = output_annual[,,30],
                    HARVESTextracted_labile_gCm2day = output[,,31],
                    mean_HARVESTextracted_labile_gCm2day = output_mean[,31],
                    mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,31],
                    HARVESTextracted_foliage_gCm2day = output[,,32],
                    mean_HARVESTextracted_foliage_gCm2day = output_mean[,32],
                    mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,32],
                    HARVESTextracted_roots_gCm2day = output[,,33],
                    mean_HARVESTextracted_roots_gCm2day = output_mean[,33],
                    mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,33],
                    HARVESTextracted_wood_gCm2day = output[,,34],
                    mean_HARVESTextracted_wood_gCm2day = output_mean[,34],
                    mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,34],
                    HARVESTextracted_litter_gCm2day = output[,,35],
                    mean_HARVESTextracted_litter_gCm2day = output_mean[,35],
                    mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,35],
                    HARVESTextracted_woodlitter_gCm2day = output[,,36],
                    mean_HARVESTextracted_woodlitter_gCm2day = output_mean[,36],
                    mean_annual_HARVESTextracted_woodlitter_gCm2day = output_annual[,,36],
                    HARVESTextracted_som_gCm2day = output[,,37],
                    mean_HARVESTextracted_som_gCm2day = output_mean[,37],
                    mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,37],
                    HARVESTlitter_labile_gCm2day = output[,,38],
                    mean_HARVESTlitter_labile_gCm2day = output_mean[,38],
                    mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,38],
                    HARVESTlitter_foliage_gCm2day = output[,,39],
                    mean_HARVESTlitter_foliage_gCm2day = output_mean[,39],
                    mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,39],
                    HARVESTlitter_roots_gCm2day = output[,,40],
                    mean_HARVESTlitter_roots_gCm2day = output_mean[,40],
                    mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,40],
                    HARVESTlitter_wood_gCm2day = output[,,41],
                    mean_HARVESTlitter_wood_gCm2day = output_mean[,41],
                    mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,41],
                    # C pools (gC/m2)
                    labile_gCm2 = output[,,42],
                    mean_labile_gCm2 = output_mean[,42],
                    mean_annual_labile_gCm2 = output_annual[,,42],
                    foliage_gCm2 = output[,,43],
                    mean_foliage_gCm2 = output_mean[,43],
                    mean_annual_foliage_gCm2 = output_annual[,,43],
                    roots_gCm2 = output[,,44],
                    mean_roots_gCm2 = output_mean[,44],
                    mean_annual_roots_gCm2 = output_annual[,,44],
                    wood_gCm2 = output[,,45],
                    mean_wood_gCm2 = output_mean[,45],
                    mean_annual_wood_gCm2 = output_annual[,,45],
                    litter_gCm2 = output[,,46],
                    mean_litter_gCm2 = output_mean[,46],
                    mean_annual_litter_gCm2 = output_annual[,,46],
                    woodlitter_gCm2 = output[,,47],
                    mean_woodlitter_gCm2 = output_mean[,47],
                    mean_annual_woodlitter_gCm2 = output_annual[,,47],
                    som_gCm2 = output[,,48],
                    mean_som_gCm2 = output_mean[,48],
                    mean_annual_som_gCm2 = output_annual[,,48],
                    # Canopy (phenology) properties
                    lai_m2m2 = output[,,49],
                    mean_lai_m2m2 = output_mean[,49],
                    mean_annual_lai_m2m2 = output_annual[,,49],
                    gsi = output[,,50],
                    mean_gsi = output_mean[,50],
                    mean_annual_gsi = output_annual[,,50],
                    gsi_itemp = output[,,51],
                    mean_gsi_itemp = output_mean[,51],
                    mean_annual_gsi_itemp = output_annual[,,51],
                    gsi_iphoto = output[,,52],
                    mean_gsi_iphoto = output_mean[,52],
                    mean_annual_gsi_iphoto = output_annual[,,52],
                    gsi_ivpd = output[,,53],
                    mean_gsi_ivpd = output_mean[,53],
                    mean_annual_gsi_ivpd = output_annual[,,53],
                    # Photosynthesis / C~water coupling related
                    gs_demand_supply_ratio = output[,,54],
                    mean_gs_demand_supply_ratio = output_mean[,54],
                    mean_annual_gs_demand_supply_ratio = output_annual[,,54],
                    gs_mmolH2Om2s = output[,,55],
                    mean_gs_mmolH2Om2s = output_mean[,55],
                    mean_annual_gs_mmolH2Om2s = output_annual[,,55],
                    APAR_MJm2day = output[,,56],
                    mean_APAR_MJm2day = output_mean[,56],
                    mean_annual_APAR_MJm2day = output_annual[,,56],
                    gb_mmolH2Om2s = output[,,57],
                    mean_gb_mmolH2Om2s = output_mean[,57],
                    mean_annual_gb_mmolH2Om2s = output_annual[,,57],
                    CiCa = output[,,58],
                    mean_CiCa = output_mean[,58],
                    mean_annual_CiCa = output_annual[,,58],
                    # Misc
                    RootDepth_m = output[,,59],
                    mean_RootDepth_m = output_mean[,59],
                    mean_annual_RootDepth_m = output_annual[,,59],
                    ## Aggregated variables
                    # Mean Transit times
                    MTT_labile_years = MTT_years[,1],
                    MTT_foliage_years = MTT_years[,2],
                    MTT_roots_years = MTT_years[,3],
                    MTT_wood_years = MTT_years[,4],
                    MTT_litter_years = MTT_years[,5],
                    MTT_woodlitter_years = MTT_years[,6],
                    MTT_som_years = MTT_years[,7],
                    # Steady state estimates
                    SS_labile_gCm2 = SS_gCm2[,1],
                    SS_foliage_gCm2 = SS_gCm2[,2],
                    SS_roots_gCm2 = SS_gCm2[,3],
                    SS_wood_gCm2 = SS_gCm2[,4],
                    SS_litter_gCm2 = SS_gCm2[,5],
                    SS_woodlitter_gCm2 = SS_gCm2[,6],
                    SS_som_gCm2 = SS_gCm2[,7])
    # Determine the NPP fraction of expressed NPP
    # i.e. actual growth not GPP-Ra
    NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                         states_all$alloc_roots_gCm2day +
                         states_all$alloc_wood_gCm2day,1,mean)
    NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day,1,mean),
                         apply(states_all$alloc_roots_gCm2day,1,mean),
                         apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
    states_all$NPP_foliage_fraction = NPP_fraction[,1]
    states_all$NPP_roots_fraction = NPP_fraction[,2]
    states_all$NPP_wood_fraction = NPP_fraction[,3]
    # Tidy up variables
    rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P3.R1.009") {
      output_dim = 72 ; MTT_dim = 7 ; SS_dim = 7
      # Load the required dalec shared object
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec9",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,fire_dim=as.integer(fire_dim)
                             ,met=as.double(t(met)),pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                             ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                            
                             ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site])
                             ,nomet=as.integer(dim(met)[2]),nofluxes=as.integer(PROJECT$model$nofluxes[site])
                             ,nopools=as.integer(PROJECT$model$nopools[site]),pft=as.integer(pft)
                             ,nodays=as.integer(dim(met)[1])
                             ,nos_years=as.integer(noyears)
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                             ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                             ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3)))
                             ,residual_waterfrac_in=as.double(c(soil_info[5],soil_info[6],soil_info[6])) 
                             ,porosity_in = as.double(c(soil_info[7],soil_info[8],soil_info[8]))
                             ,pore_size_dist_in = as.double(c(soil_info[9],soil_info[10],soil_info[10]))
                             ,air_entry_in = as.double(c(soil_info[11],soil_info[12],soil_info[12]))
                             ,sat_conductivity_in = as.double(c(soil_info[13],soil_info[14],soil_info[14])))
      # Extract the different output variables
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rg_foliage_gCm2day = output[,,3],
                      mean_rg_foliage_gCm2day = output_mean[,3],
                      mean_annual_rg_foliage_gCm2day = output_annual[,,3],
                      rhet_litter_gCm2day = output[,,4],
                      mean_rhet_litter_gCm2day = output_mean[,4],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,4],
                      rhet_som_gCm2day = output[,,5],
                      mean_rhet_som_gCm2day = output_mean[,5],
                      mean_annual_rhet_som_gCm2day = output_annual[,,5],
                      rhet_woodlitter_gCm2day = output[,,6],
                      mean_rhet_woodlitter_gCm2day = output_mean[,6],
                      mean_annual_rhet_woodlitter_gCm2day = output_annual[,,6],
                      fire_gCm2day = output[,,7],
                      mean_fire_gCm2day = output_mean[,7],
                      mean_annual_fire_gCm2day = output_annual[,,7],
                      harvest_gCm2day = output[,,8],
                      mean_harvest_gCm2day = output_mean[,8],
                      mean_annual_harvest_gCm2day = output_annual[,,8],
                      # Internal fluxes
                      alloc_labile_gCm2day = output[,,9],
                      mean_alloc_labile_gCm2day = output_mean[,9],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,9],
                      alloc_roots_gCm2day = output[,,10],
                      mean_alloc_roots_gCm2day = output_mean[,10],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,10],
                      alloc_wood_gCm2day = output[,,11],
                      mean_alloc_wood_gCm2day = output_mean[,11],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,11],
                      labile_to_foliage_gCm2day = output[,,12],
                      mean_labile_to_foliage_gCm2day = output_mean[,12],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,12],
                      foliage_to_litter_gCm2day = output[,,13],
                      mean_foliage_to_litter_gCm2day = output_mean[,13],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,13],
                      roots_to_litter_gCm2day = output[,,14],
                      mean_roots_to_litter_gCm2day = output_mean[,14],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,14],
                      wood_to_litter_gCm2day = output[,,15],
                      mean_wood_to_litter_gCm2day = output_mean[,15],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,15],
                      litter_to_som_gCm2day = output[,,16],
                      mean_litter_to_som_gCm2day = output_mean[,16],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,16],
                      woodlitter_to_som_gCm2day = output[,,17],
                      mean_woodlitter_to_som_gCm2day = output_mean[,17],
                      mean_annual_woodlitter_to_som_gCm2day = output_annual[,,17],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,18],
                      mean_FIREemiss_labile_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,18],
                      FIRElitter_labile_gCm2day = output[,,19],
                      mean_FIRElitter_labile_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,19],
                      FIREemiss_foliage_gCm2day = output[,,20],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,20],
                      FIRElitter_foliage_gCm2day = output[,,21],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,21],
                      FIREemiss_roots_gCm2day = output[,,22],
                      mean_FIREemiss_roots_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,22],
                      FIRElitter_roots_gCm2day = output[,,23],
                      mean_FIRElitter_roots_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,23],
                      FIREemiss_wood_gCm2day = output[,,24],
                      mean_FIREemiss_wood_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,24],
                      FIRElitter_wood_gCm2day = output[,,25],
                      mean_FIRElitter_wood_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,25],
                      FIREemiss_litter_gCm2day = output[,,26],
                      mean_FIREemiss_litter_gCm2day = output_mean[,26],
                      Fmean_annual_IREemiss_litter_gCm2day = output_annual[,,26],
                      FIRElitter_litter_gCm2day = output[,,27],
                      mean_FIRElitter_litter_gCm2day = output_mean[,27],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,27],
                      FIREemiss_woodlitter_gCm2day = output[,,28],
                      mean_FIREemiss_woodlitter_gCm2day = output_mean[,28],
                      mean_annual_FIREemiss_woodlitter_gCm2day = output_annual[,,28],
                      FIRElitter_woodlitter_gCm2day = output[,,29],
                      mean_FIRElitter_woodlitter_gCm2day = output_mean[,29],
                      mean_annual_FIRElitter_woodlitter_gCm2day = output_annual[,,29],
                      FIREemiss_som_gCm2day = output[,,30],
                      mean_FIREemiss_som_gCm2day = output_mean[,30],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,30],
                      HARVESTextracted_labile_gCm2day = output[,,31],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,31],
                      HARVESTextracted_foliage_gCm2day = output[,,32],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,32],
                      HARVESTextracted_roots_gCm2day = output[,,33],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,33],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,33],
                      HARVESTextracted_wood_gCm2day = output[,,34],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,34],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,34],
                      HARVESTextracted_litter_gCm2day = output[,,35],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,35],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,35],
                      HARVESTextracted_woodlitter_gCm2day = output[,,36],
                      mean_HARVESTextracted_woodlitter_gCm2day = output_mean[,36],
                      mean_annual_HARVESTextracted_woodlitter_gCm2day = output_annual[,,36],
                      HARVESTextracted_som_gCm2day = output[,,37],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,37],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,37],
                      HARVESTlitter_labile_gCm2day = output[,,38],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,38],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,38],
                      HARVESTlitter_foliage_gCm2day = output[,,39],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,39],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,39],
                      HARVESTlitter_roots_gCm2day = output[,,40],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,40],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,40],
                      HARVESTlitter_wood_gCm2day = output[,,41],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,41],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,41],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,42],
                      mean_labile_gCm2 = output_mean[,42],
                      mean_annual_labile_gCm2 = output_annual[,,42],
                      foliage_gCm2 = output[,,43],
                      mean_foliage_gCm2 = output_mean[,43],
                      mean_annual_foliage_gCm2 = output_annual[,,43],
                      roots_gCm2 = output[,,44],
                      mean_roots_gCm2 = output_mean[,44],
                      mean_annual_roots_gCm2 = output_annual[,,44],
                      wood_gCm2 = output[,,45],
                      mean_wood_gCm2 = output_mean[,45],
                      mean_annual_wood_gCm2 = output_annual[,,45],
                      litter_gCm2 = output[,,46],
                      mean_litter_gCm2 = output_mean[,46],
                      mean_annual_litter_gCm2 = output_annual[,,46],
                      woodlitter_gCm2 = output[,,47],
                      mean_woodlitter_gCm2 = output_mean[,47],
                      mean_annual_woodlitter_gCm2 = output_annual[,,47],
                      som_gCm2 = output[,,48],
                      mean_som_gCm2 = output_mean[,48],
                      mean_annual_som_gCm2 = output_annual[,,48],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,49],
                      mean_ET_kgH2Om2day = output_mean[,49],
                      mean_annual_ET_kgH2Om2day = output_annual[,,49],
                      Etrans_kgH2Om2day = output[,,50],
                      mean_Etrans_kgH2Om2day = output_mean[,50],
                      mean_annual_Etrans_kgH2Om2day = output_annual[,,50],
                      Esoil_kgH2Om2day = output[,,51],
                      mean_Esoil_kgH2Om2day = output_mean[,51],
                      mean_annual_Esoil_kgH2Om2day = output_annual[,,51],
                      Ewetcanopy_kgH2Om2day = output[,,52],
                      mean_Ewetcanopy_kgH2Om2day = output_mean[,52],
                      mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,52],
                      runoff_kgH2Om2day = output[,,53],
                      mean_runoff_kgH2Om2day = output_mean[,53],
                      mean_annual_runoff_kgH2Om2day = output_annual[,,53],
                      underflow_kgH2Om2day = output[,,54],
                      mean_underflow_kgH2Om2day = output_mean[,54],
                      mean_annual_underflow_kgH2Om2day = output_annual[,,54],
                      SurfDrainage_kgH2Om2day = output[,,55],
                      mean_SurfDrainage_kgH2Om2day = output_mean[,55],
                      mean_annual_SurfDrainage_kgH2Om2day = output_annual[,,55],
                      SurfInfiltrated_kgH2Om2day = output[,,56],
                      mean_SurfInfiltrated_kgH2Om2day = output_mean[,56],
                      mean_annual_SurfInfiltrated_kgH2Om2day = output_annual[,,56],
                      Etrans_1st_root_layer_uptake_fraction = output[,,57],
                      mean_Etrans_1st_root_layer_uptake_fraction = output_mean[,57],
                      mean_annual_Etrans_1st_root_layer_uptake_fraction = output_annual[,,57],
                      Etrans_2nd_root_layer_uptake_fraction = output[,,58],
                      mean_Etrans_2nd_root_layer_uptake_fraction = output_mean[,58],
                      mean_annual_Etrans_2nd_root_layer_uptake_fraction = output_annual[,,58],
                      SurfWater_kgH2Om2 = output[,,59],
                      mean_SurfWater_kgH2Om2 = output_mean[,59],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,59],
                      wSWP_MPa = output[,,60],
                      mean_wSWP_MPa = output_mean[,60],
                      mean_annual_wSWP_MPa = output_annual[,,60],
                      snow_kgH2Om2 = output[,,61],
                      mean_snow_kgH2Om2 = output_mean[,61],
                      mean_annual_snow_kgH2Om2 = output_annual[,,61],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,62],
                      mean_lai_m2m2 = output_mean[,62],
                      mean_annual_lai_m2m2 = output_annual[,,62],
                      gsi = output[,,63],
                      mean_gsi = output_mean[,63],
                      mean_annual_gsi = output_annual[,,63],
                      gsi_itemp = output[,,64],
                      mean_gsi_itemp = output_mean[,64],
                      mean_annual_gsi_itemp = output_annual[,,64],
                      gsi_iphoto = output[,,65],
                      mean_gsi_iphoto = output_mean[,65],
                      mean_annual_gsi_iphoto = output_annual[,,65],
                      gsi_ivpd = output[,,66],
                      mean_gsi_ivpd = output_mean[,66],
                      mean_annual_gsi_ivpd = output_annual[,,66],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,67],
                      mean_gs_demand_supply_ratio = output_mean[,67],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,67],
                      gs_mmolH2Om2s = output[,,68],
                      mean_gs_mmolH2Om2s = output_mean[,68],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,68],
                      APAR_MJm2day = output[,,69],
                      mean_APAR_MJm2day = output_mean[,69],
                      mean_annual_APAR_MJm2day = output_annual[,,69],
                      gb_mmolH2Om2s = output[,,70],
                      mean_gb_mmolH2Om2s = output_mean[,70],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,70],
                      CiCa = output[,,71],
                      mean_CiCa = output_mean[,71],
                      mean_annual_CiCa = output_annual[,,71],
                      # Misc
                      RootDepth_m = output[,,72],
                      mean_RootDepth_m = output_mean[,72],
                      mean_annual_RootDepth_m = output_annual[,,72],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_woodlitter_years = MTT_years[,6],
                      MTT_som_years = MTT_years[,7],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_woodlitter_gCm2 = SS_gCm2[,6],
                      SS_som_gCm2 = SS_gCm2[,7])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.C1.D1.F2.P1.002") {
      output_dim = 44 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec2",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                             ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                                 
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodays=as.integer(dim(met)[1])
                             ,nos_years=as.integer(noyears)
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter) )
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,43],
                      mean_lai_m2m2 = output_mean[,43],
                      mean_annual_lai_m2m2 = output_annual[,,43],
                      # Photosynthesis / C~water coupling related
                      CiCa = output[,,44],
                      mean_CiCa = output_mean[,44],
                      mean_annual_CiCa = output_annual[,,44],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC_1005") {
      output_dim = 46 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec1005",output_dim=as.integer(output_dim)
                                ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                                ,met=as.double(t(met))
                                ,pars=as.double(pars_in)
                                ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                                ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                                ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                                ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                       
                                ,lat=as.double(lat)
                                ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                ,nodays=as.integer(dim(met)[1])
                                ,nos_years=as.integer(noyears)
                                ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
      # extract output variables
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output_mean[,11],
                      mean_labile_to_foliage_gCm2day = output[,,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      mean_ET_kgH2Om2day = output_mean[,43],
                      mean_annual_ET_kgH2Om2day = output_annual[,,43],
                      SurfWater_kgH2Om2 = output[,,44],
                      mean_SurfWater_kgH2Om2 = output_mean[,44],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,44],
                      # wSWP_MPa = output[,,45],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,45],
                      mean_lai_m2m2 = output_mean[,45],
                      mean_annual_lai_m2m2 = output_annual[,,45],
                      # Photosynthesis / C~water coupling related
                      CiCa = output[,,46],
                      mean_CiCa = output_mean[,46],
                      mean_annual_CiCa = output_annual[,,46],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC_1005a") {
      output_dim = 46 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      # if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdalec1005a",output_dim=as.integer(output_dim)
                                 ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                                 ,met=as.double(t(met))
                                 ,pars=as.double(pars_in)
                                 ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                                 ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                 ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                                 ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                                 ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                        
                                 ,lat=as.double(lat)
                                 ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                 ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                 ,nodays=as.integer(dim(met)[1])
                                 ,nos_years=as.integer(noyears)
                                 ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
      # extract output variables
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output_mean[,11],
                      mean_labile_to_foliage_gCm2day = output[,,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      mean_ET_kgH2Om2day = output_mean[,43],
                      mean_annual_ET_kgH2Om2day = output_annual[,,43],
                      SurfWater_kgH2Om2 = output[,,44],
                      mean_SurfWater_kgH2Om2 = output_mean[,44],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,44],
                      # wSWP_MPa = output[,,45],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,45],
                      mean_lai_m2m2 = output_mean[,45],
                      mean_annual_lai_m2m2 = output_annual[,,45],
                      # Photosynthesis / C~water coupling related
                      CiCa = output[,,46],
                      mean_CiCa = output_mean[,46],
                      mean_annual_CiCa = output_annual[,,46],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # add newly calculated variables
      # states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Final tidy
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C1.D2.F2.H1.P1.003") {
      output_dim = 48 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec3",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                             ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                    
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodays=as.integer(dim(met)[1])
                             ,nos_years=as.integer(noyears)
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                             ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                             ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3)))
                             ,residual_waterfrac_in=as.double(c(soil_info[5],soil_info[6],soil_info[6])) 
                             ,porosity_in = as.double(c(soil_info[7],soil_info[8],soil_info[8]))
                             ,pore_size_dist_in = as.double(c(soil_info[9],soil_info[10],soil_info[10]))
                             ,air_entry_in = as.double(c(soil_info[11],soil_info[12],soil_info[12]))
                             ,sat_conductivity_in = as.double(c(soil_info[13],soil_info[14],soil_info[14])))
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))      
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output_mean[,11],
                      mean_labile_to_foliage_gCm2day = output[,,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter__annualgCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,43],
                      mean_lai_m2m2 = output_mean[,43],
                      mean_annual_lai_m2m2 = output_annual[,,43],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,44],
                      mean_gs_demand_supply_ratio = output_mean[,44],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,44],
                      gs_mmolH2Om2s = output[,,45],
                      mean_gs_mmolH2Om2s = output_mean[,45],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,45],
                      APAR_MJm2day = output[,,46],
                      mean_APAR_MJm2day = output_mean[,46],
                      mean_annual_APAR_MJm2day = output_annual[,,46],
                      gb_mmolH2Om2s = output[,,47],
                      mean_gb_mmolH2Om2s = output_mean[,47],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,47],
                      CiCa = output[,,48],
                      mean_CiCa = output_mean[,48],
                      mean_annual_CiCa = output_annual[,,48],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # Determine the NPP fraction of expressed NPPHARVESTextracted_roots_gCm2day = output[,,29],
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C1.D2.F2.H2.P1.004") {
      output_dim = 68 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec4",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                             ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                    
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodiags=as.integer(PROJECT$model$nodiags[site]),nodays=as.integer(dim(met)[1])
                             ,nos_years=as.integer(noyears)
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                             ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                             ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2]))
                             ,residual_waterfrac_in=as.double(c(soil_info[5],soil_info[6],soil_info[6])) 
                             ,porosity_in = as.double(c(soil_info[7],soil_info[8],soil_info[8]))
                             ,pore_size_dist_in = as.double(c(soil_info[9],soil_info[10],soil_info[10]))
                             ,air_entry_in = as.double(c(soil_info[11],soil_info[12],soil_info[12]))
                             ,sat_conductivity_in = as.double(c(soil_info[13],soil_info[14],soil_info[14])))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))      
      MTT_years = tmp$out_var2     ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3       ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      mean_ET_kgH2Om2day = output_mean[,43],
                      mean_annual_ET_kgH2Om2day = output_annual[,,43],
                      Etrans_kgH2Om2day = output[,,44],
                      mean_Etrans_kgH2Om2day = output_mean[,44],
                      mean_annual_Etrans_kgH2Om2day = output_annual[,,44],
                      Esoil_kgH2Om2day = output[,,45],
                      mean_Esoil_kgH2Om2day = output_mean[,45],
                      mean_annual_Esoil_kgH2Om2day = output_annual[,,45],
                      Ewetcanopy_kgH2Om2day = output[,,46],
                      mean_Ewetcanopy_kgH2Om2day = output_mean[,46],
                      mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,46],
                      runoff_kgH2Om2day = output[,,47],
                      mean_runoff_kgH2Om2day = output_mean[,47],
                      mean_annual_runoff_kgH2Om2day = output_annual[,,47],
                      underflow_kgH2Om2day = output[,,48],
                      mean_underflow_kgH2Om2day = output_mean[,48],
                      mean_annual_underflow_kgH2Om2day = output_annual[,,48],
                      SurfDrainage_kgH2Om2day = output[,,49],
                      mean_SurfDrainage_kgH2Om2day = output_mean[,49],
                      mean_annual_SurfDrainage_kgH2Om2day = output_annual[,,49],
                      SurfInfiltrated_kgH2Om2day = output[,,50],
                      mean_SurfInfiltrated_kgH2Om2day = output_mean[,50],
                      mean_annual_SurfInfiltrated_kgH2Om2day = output_annual[,,50],
                      Etrans_1st_root_layer_uptake_fraction = output[,,51],
                      mean_Etrans_1st_root_layer_uptake_fraction = output_mean[,51],
                      mean_annual_Etrans_1st_root_layer_uptake_fraction = output_annual[,,51],
                      Etrans_2nd_root_layer_uptake_fraction = output[,,52],
                      mean_Etrans_2nd_root_layer_uptake_fraction = output_mean[,52],
                      mean_annual_Etrans_2nd_root_layer_uptake_fraction = output_annual[,,52],
                      SurfWater_kgH2Om2 = output[,,53],
                      mean_SurfWater_kgH2Om2 = output_mean[,53],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,53],
                      wSWP_MPa = output[,,54],
                      mean_wSWP_MPa = output_mean[,54],
                      mean_annual_wSWP_MPa = output_annual[,,54],
                      snow_kgH2Om2 = output[,,55],
                      mean_snow_kgH2Om2 = output_mean[,55],
                      mean_annual_snow_kgH2Om2 = output_annual[,,55],
                      hydraulic_conductivity_ms = output[,,66],
                      soil_waterpot_MPa = output[,,67],
                      field_capacity = output[,,68],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,56],
                      mean_lai_m2m2 = output_mean[,56],
                      mean_annual_lai_m2m2 = output_annual[,,56],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,57],
                      mean_gs_demand_supply_ratio = output_mean[,57],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,57],
                      gs_mmolH2Om2s = output[,,58],
                      mean_gs_mmolH2Om2s = output_mean[,58],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,58],
                      APAR_MJm2day = output[,,59],
                      mean_APAR_MJm2day = output_mean[,59],
                      mean_annual_APAR_MJm2day = output_annual[,,59],
                      gb_mmolH2Om2s = output[,,60],
                      mean_gb_mmolH2Om2s = output_mean[,60],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,60],
                      CiCa = output[,,61],
                      mean_CiCa = output_mean[,61],
                      mean_annual_CiCa = output_annual[,,61],
                      # Misc
                      RootDepth_m = output[,,62],
                      mean_RootDepth_m = output_mean[,62],
                      mean_annual_RootDepth_m = output_annual[,,62],
                      # Leaf water potential
                      LWP_MPa = output[,,63],
                      mean_LWP_MPa = output_mean[,63],
                      mean_annual_LWP_MPa = output_annual[,,63], 
                      # Canopy aerodynamic diagnostics
                      canopy_area_scaling_wind = output[,,64],
                      mean_canopy_area_scaling_wind = output_mean[,64],
                      mean_annual_canopy_area_scaling_wind = output_annual[,,64],
                      canopy_area_scaling_light = output[,,65],
                      mean_canopy_area_scaling_light = output_mean[,65],
                      mean_annual_canopy_area_scaling_light = output_annual[,,65],    
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C1.D2.F2.H3.P1.029") {
      output_dim = 63 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec29",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                             ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                    
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodays=as.integer(dim(met)[1])
                             ,nos_years=as.integer(noyears)
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                             ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                             ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2]))
                             ,residual_waterfrac_in=as.double(c(soil_info[5],soil_info[6],soil_info[6])) 
                             ,porosity_in = as.double(c(soil_info[7],soil_info[8],soil_info[8]))
                             ,pore_size_dist_in = as.double(c(soil_info[9],soil_info[10],soil_info[10]))
                             ,air_entry_in = as.double(c(soil_info[11],soil_info[12],soil_info[12]))
                             ,sat_conductivity_in = as.double(c(soil_info[13],soil_info[14],soil_info[14])))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))      
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],  
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output[,,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mea_meann[,,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      mean_ET_kgH2Om2day = output_mean[,43],
                      mean_annual_ET_kgH2Om2day = output_annual[,,43],
                      Etrans_kgH2Om2day = output[,,44],
                      mean_Etrans_kgH2Om2day = output_mean[,44],
                      mean_annual_Etrans_kgH2Om2day = output_annual[,,44],
                      Esoil_kgH2Om2day = output[,,45],
                      mean_Esoil_kgH2Om2day = output_mean[,45],
                      mean_annual_Esoil_kgH2Om2day = output_annual[,,45],
                      Ewetcanopy_kgH2Om2day = output[,,46],
                      mean_Ewetcanopy_kgH2Om2day = output_mean[,46],
                      mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,46],
                      runoff_kgH2Om2day = output[,,47],
                      mean_runoff_kgH2Om2day = output_mean[,47],
                      mean_annual_runoff_kgH2Om2day = output_annual[,,47],
                      underflow_kgH2Om2day = output[,,48],
                      mean_underflow_kgH2Om2day = output_mean[,48],
                      mean_annual_underflow_kgH2Om2day = output_annual[,,48],
                      SurfDrainage_kgH2Om2day = output[,,49],
                      mean_SurfDrainage_kgH2Om2day = output_mean[,49],
                      mean_annual_SurfDrainage_kgH2Om2day = output_annual[,,49],
                      SurfInfiltrated_kgH2Om2day = output[,,50],
                      mean_SurfInfiltrated_kgH2Om2day = output_mean[,50],
                      mean_annual_SurfInfiltrated_kgH2Om2day = output_annual[,,50],
                      Etrans_1st_root_layer_uptake_fraction = output[,,51],
                      mean_Etrans_1st_root_layer_uptake_fraction = output_mean[,51],
                      mean_annual_Etrans_1st_root_layer_uptake_fraction = output_annual[,,51],
                      Etrans_2nd_root_layer_uptake_fraction = output[,,52],
                      mean_Etrans_2nd_root_layer_uptake_fraction = output_mean[,52],
                      mean_annual_Etrans_2nd_root_layer_uptake_fraction = output_annual[,,52],
                      SurfWater_kgH2Om2 = output[,,53],
                      mean_SurfWater_kgH2Om2 = output_mean[,53],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,53],
                      wSWP_MPa = output[,,54],
                      mean_wSWP_MPa = output_mean[,54],
                      mean_annual_wSWP_MPa = output_annual[,,54],
                      snow_kgH2Om2 = output[,,55],
                      mean_snow_kgH2Om2 = output_mean[,55],
                      mean_annual_snow_kgH2Om2 = output_annual[,,55],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,56],
                      mean_lai_m2m2 = output_mean[,56],
                      mean_annual_lai_m2m2 = output_annual[,,56],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,57],
                      mean_gs_demand_supply_ratio = output_mean[,57],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,57],
                      gs_mmolH2Om2s = output[,,58],
                      mean_gs_mmolH2Om2s = output_mean[,58],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,58],
                      APAR_MJm2day = output[,,59],
                      mean_APAR_MJm2day = output_mean[,59],
                      mean_annual_APAR_MJm2day = output_annual[,,59],
                      gb_mmolH2Om2s = output[,,60],
                      mean_gb_mmolH2Om2s = output_mean[,60],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,60],
                      CiCa = output[,,61],
                      mean_CiCa = output_mean[,61],
                      mean_annual_CiCa = output_annual[,,61],
                      # Misc
                      RootDepth_m = output[,,62],
                      mean_RootDepth_m = output_mean[,62],
                      mean_annual_RootDepth_m = output_annual[,,62],                    
                      LWP_MPa = output[,,63],
                      mean_LWP_MPa = output_mean[,63],
                      mean_annual_LWP_MPa = output_annual[,,63],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A2.C1.D2.F2.H2.P1.020") {
      output_dim = 62 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec20",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                              ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                     
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,nos_years=as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                              ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2]))
                              ,residual_waterfrac_in=as.double(c(soil_info[5],soil_info[6],soil_info[6])) 
                              ,porosity_in = as.double(c(soil_info[7],soil_info[8],soil_info[8]))
                              ,pore_size_dist_in = as.double(c(soil_info[9],soil_info[10],soil_info[10]))
                              ,air_entry_in = as.double(c(soil_info[11],soil_info[12],soil_info[12]))
                              ,sat_conductivity_in = as.double(c(soil_info[13],soil_info[14],soil_info[14])))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))      
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],  
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output[,,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      mean_ET_kgH2Om2day = output_mean[,43],
                      mean_annual_ET_kgH2Om2day = output_annual[,,43],
                      Etrans_kgH2Om2day = output[,,44],
                      mean_Etrans_kgH2Om2day = output_mean[,44],
                      mean_annual_Etrans_kgH2Om2day = output_annual[,,44],
                      Esoil_kgH2Om2day = output[,,45],
                      mean_Esoil_kgH2Om2day = output_mean[,45],
                      mean_annual_Esoil_kgH2Om2day = output_annual[,,45],
                      Ewetcanopy_kgH2Om2day = output[,,46],
                      mean_Ewetcanopy_kgH2Om2day = output_mean[,46],
                      mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,46],
                      runoff_kgH2Om2day = output[,,47],
                      mean_runoff_kgH2Om2day = output_mean[,47],
                      mean_annual_runoff_kgH2Om2day = output_annual[,,47],
                      underflow_kgH2Om2day = output[,,48],
                      mean_underflow_kgH2Om2day = output_mean[,48],
                      mean_annual_underflow_kgH2Om2day = output_annual[,,48],
                      SurfDrainage_kgH2Om2day = output[,,49],
                      mean_SurfDrainage_kgH2Om2day = output_mean[,49],
                      mean_annual_SurfDrainage_kgH2Om2day = output_annual[,,49],
                      SurfInfiltrated_kgH2Om2day = output[,,50],
                      mean_SurfInfiltrated_kgH2Om2day = output_mean[,50],
                      mean_annual_SurfInfiltrated_kgH2Om2day = output_annual[,,50],
                      Etrans_1st_root_layer_uptake_fraction = output[,,51],
                      mean_Etrans_1st_root_layer_uptake_fraction = output_mean[,51],
                      mean_annual_Etrans_1st_root_layer_uptake_fraction = output_annual[,,51],
                      Etrans_2nd_root_layer_uptake_fraction = output[,,52],
                      mean_Etrans_2nd_root_layer_uptake_fraction = output_mean[,52],
                      mean_annual_Etrans_2nd_root_layer_uptake_fraction = output_annual[,,52],
                      SurfWater_kgH2Om2 = output[,,53],
                      mean_SurfWater_kgH2Om2 = output_mean[,53],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,53],
                      wSWP_MPa = output[,,54],
                      mean_wSWP_MPa = output_mean[,54],
                      mean_annual_wSWP_MPa = output_annual[,,54],
                      snow_kgH2Om2 = output[,,55],
                      mean_snow_kgH2Om2 = output_mean[,55],
                      mean_annual_snow_kgH2Om2 = output_annual[,,55],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,56],
                      mean_lai_m2m2 = output_mean[,56],
                      mean_annual_lai_m2m2 = output_annual[,,56],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,57],
                      mean_gs_demand_supply_ratio = output_mean[,57],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,57],
                      gs_mmolH2Om2s = output[,,58],
                      mean_gs_mmolH2Om2s = output_mean[,58],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,58],
                      APAR_MJm2day = output[,,59],
                      mean_APAR_MJm2day = output_mean[,59],
                      mean_annual_APAR_MJm2day = output_annual[,,59],
                      gb_mmolH2Om2s = output[,,60],
                      mean_gb_mmolH2Om2s = output_mean[,60],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,60],
                      CiCa = output[,,61],
                      mean_CiCa = output_mean[,61],
                      mean_annual_CiCa = output_annual[,,61],
                      # Misc
                      RootDepth_m = output[,,62],
                      mean_RootDepth_m = output_mean[,62],
                      mean_annual_RootDepth_m = output_annual[,,62],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A4.C6.D2.F2.H2.P11.031") {
      output_dim = 75 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec31",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                              ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodiags=as.integer(PROJECT$model$nodiags[site]),nodays=as.integer(dim(met)[1])
                              ,nos_years = as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                              ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2]))
                              ,residual_waterfrac_in=as.double(c(soil_info[5],soil_info[6],soil_info[6])) 
                              ,porosity_in = as.double(c(soil_info[7],soil_info[8],soil_info[8]))
                              ,pore_size_dist_in = as.double(c(soil_info[9],soil_info[10],soil_info[10]))
                              ,air_entry_in = as.double(c(soil_info[11],soil_info[12],soil_info[12]))
                              ,sat_conductivity_in = as.double(c(soil_info[13],soil_info[14],soil_info[14])))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))
      MTT_years = tmp$out_var2     ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3       ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      mean_ET_kgH2Om2day = output_mean[,43],
                      mean_annual_ET_kgH2Om2day = output_annual[,,43],
                      Etrans_kgH2Om2day = output[,,44],
                      mean_Etrans_kgH2Om2day = output_mean[,44],
                      mean_annual_Etrans_kgH2Om2day = output_annual[,,44],
                      Esoil_kgH2Om2day = output[,,45],
                      mean_Esoil_kgH2Om2day = output_mean[,45],
                      mean_annual_Esoil_kgH2Om2day = output_annual[,,45],
                      Ewetcanopy_kgH2Om2day = output[,,46],
                      mean_Ewetcanopy_kgH2Om2day = output_mean[,46],
                      mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,46],
                      runoff_kgH2Om2day = output[,,47],
                      mean_runoff_kgH2Om2day = output_mean[,47],
                      mean_annual_runoff_kgH2Om2day = output_annual[,,47],
                      underflow_kgH2Om2day = output[,,48],
                      mean_underflow_kgH2Om2day = output_mean[,48],
                      mean_annual_underflow_kgH2Om2day = output_annual[,,48],
                      SurfDrainage_kgH2Om2day = output[,,49],
                      mean_SurfDrainage_kgH2Om2day = output_mean[,49],
                      mean_annual_SurfDrainage_kgH2Om2day = output_annual[,,49],
                      SurfInfiltrated_kgH2Om2day = output[,,50],
                      mean_SurfInfiltrated_kgH2Om2day = output_mean[,50],
                      mean_annual_SurfInfiltrated_kgH2Om2day = output_annual[,,50],
                      Etrans_1st_root_layer_uptake_fraction = output[,,51],
                      mean_Etrans_1st_root_layer_uptake_fraction = output_mean[,51],
                      mean_annual_Etrans_1st_root_layer_uptake_fraction = output_annual[,,51],
                      Etrans_2nd_root_layer_uptake_fraction = output[,,52],
                      mean_Etrans_2nd_root_layer_uptake_fraction = output_mean[,52],
                      mean_annual_Etrans_2nd_root_layer_uptake_fraction = output_annual[,,52],
                      SurfWater_kgH2Om2 = output[,,53],
                      mean_SurfWater_kgH2Om2 = output_mean[,53],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,53],
                      wSWP_MPa = output[,,54],
                      mean_wSWP_MPa = output_mean[,54],
                      mean_annual_wSWP_MPa = output_annual[,,54],
                      snow_kgH2Om2 = output[,,55],
                      mean_snow_kgH2Om2 = output_mean[,55],
                      mean_annual_snow_kgH2Om2 = output_annual[,,55],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,56],
                      mean_lai_m2m2 = output_mean[,56],
                      mean_annual_lai_m2m2 = output_annual[,,56],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,57],
                      mean_gs_demand_supply_ratio = output_mean[,57],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,57],
                      gs_mmolH2Om2s = output[,,58],
                      mean_gs_mmolH2Om2s = output_mean[,58],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,58],
                      APAR_MJm2day = output[,,59],
                      mean_APAR_MJm2day = output_mean[,59],
                      mean_annual_APAR_MJm2day = output_annual[,,59],
                      gb_mmolH2Om2s = output[,,60],
                      mean_gb_mmolH2Om2s = output_mean[,60],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,60],
                      CiCa = output[,,61],
                      mean_CiCa = output_mean[,61],
                      mean_annual_CiCa = output_annual[,,61],
                      # Misc
                      RootDepth_m = output[,,62],
                      mean_RootDepth_m = output_mean[,62],
                      mean_annual_RootDepth_m = output_annual[,,62],
                      leaf_temperature_celcius = output[,,63],
                      mean_leaf_temperature_celcius = output_mean[,63],
                      mean_annual_leaf_temperature_celcius = output_annual[,,63],                      
                      soil_temperature_celcius = output[,,64],
                      mean_soil_temperature_celcius = output_mean[,64],
                      mean_annual_soil_temperature_celcius = output_annual[,,64],   
                      # Plant water stress
                      LWP_MPa = output[,,65],
                      mean_LWP_MPa = output_mean[,65],
                      mean_annual_LWP_MPa = output_annual[,,65],
                      # C allocation diagnositics
                      LabBio_limitation = output[,,66],
                      mean_LabBio_limitation = output_mean[,66],
                      mean_annual_LabBio_limitation = output_annual[,,66],
                      foliage_leafT_limitation = output[,,67],
                      mean_foliage_leafT_limitation = output_mean[,67],
                      mean_annual_foliage_leafT_limitation = output_annual[,,67],
                      roots_leafT_limitation = output[,,68],
                      mean_roots_leafT_limitation = output_mean[,68],
                      mean_annual_roots_leafT_limitation = output_annual[,,68],
                      wood_leafT_limitation = output[,,69],
                      mean_wood_leafT_limitation = output_mean[,69],
                      mean_annual_wood_leafT_limitation = output_annual[,,69],
                      foliage_wSWP_limitation = output[,,70],
                      mean_foliage_wSWP_limitation = output_mean[,70],
                      mean_annual_foliage_wSWP_limitation = output_annual[,,70],
                      roots_wSWP_limitation = output[,,71],
                      mean_roots_wSWP_limitation = output_mean[,71],
                      mean_annual_roots_wSWP_limitation = output_annual[,,71],
                      wood_wSWP_limitation = output[,,72],
                      mean_wood_wSWP_limitation = output_mean[,72],
                      mean_annual_wood_wSWP_limitation = output_annual[,,72],
                      gpp_return_gCm2day = output[,,73],
                      mean_gpp_return_gCm2day = output_mean[,73],
                      mean_annual_gpp_return_gCm2day = output_annual[,,73],
                      canopy_area_scaling_wind = output[,,74],
                      mean_canopy_area_scaling_wind = output_mean[,74],
                      mean_annual_canopy_area_scaling_wind = output_annual[,,74],
                      canopy_area_scaling_light = output[,,75],
                      mean_canopy_area_scaling_light = output_mean[,75],
                      mean_annual_canopy_area_scaling_light = output_annual[,,75],                      
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,output_mean,output_annual,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A4.C6.D2.F2.H3.P12.033") {
      output_dim = 83 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec33",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                              ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodiags=as.integer(PROJECT$model$nodiags[site]),nodays=as.integer(dim(met)[1])
                              ,nos_years = as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                              ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2]))
                              ,residual_waterfrac_in=as.double(c(soil_info[5],soil_info[6],soil_info[6])) 
                              ,porosity_in = as.double(c(soil_info[7],soil_info[8],soil_info[8]))
                              ,pore_size_dist_in = as.double(c(soil_info[9],soil_info[10],soil_info[10]))
                              ,air_entry_in = as.double(c(soil_info[11],soil_info[12],soil_info[12]))
                              ,sat_conductivity_in = as.double(c(soil_info[13],soil_info[14],soil_info[14])))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))
      MTT_years = tmp$out_var2     ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3       ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      rgrow_gCm2day = output[,,11],
                      mean_rgrow_gCm2day = output_mean[,11],
                      mean_annual_rgrow_gCm2day = output_annual[,,11],
                      rmain_from_labile_gCm2day = output[,,12],
                      mean_rmain_from_labile_gCm2day = output_mean[,12],
                      mean_annual_rmain_from_labile_gCm2day = output_annual[,,12],                      
                      foliage_to_litter_gCm2day = output[,,13],
                      mean_foliage_to_litter_gCm2day = output_mean[,13],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,13],
                      roots_to_litter_gCm2day = output[,,14],
                      mean_roots_to_litter_gCm2day = output_mean[,14],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,15],
                      mean_wood_to_litter_gCm2day = output_mean[,15],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,15],
                      litter_to_som_gCm2day = output[,,16],
                      mean_litter_to_som_gCm2day = output_mean[,16],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,16],
                      rmain_wood_root_gCm2day = output[,,17],
                      mean_rmain_wood_root_gCm2day = output_mean[,17],
                      mean_annual_rmain_wood_root_gCm2day = output_annual[,,17],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,18],
                      mean_FIREemiss_labile_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,18],
                      FIRElitter_labile_gCm2day = output[,,19],
                      mean_FIRElitter_labile_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,19],
                      FIREemiss_foliage_gCm2day = output[,,20],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,20],
                      FIRElitter_foliage_gCm2day = output[,,21],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,21],
                      FIREemiss_roots_gCm2day = output[,,22],
                      mean_FIREemiss_roots_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,22],
                      FIRElitter_roots_gCm2day = output[,,23],
                      mean_FIRElitter_roots_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,23],
                      FIREemiss_wood_gCm2day = output[,,24],
                      mean_FIREemiss_wood_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,24],
                      FIRElitter_wood_gCm2day = output[,,25],
                      mean_FIRElitter_wood_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,25],
                      FIREemiss_litter_gCm2day = output[,,26],
                      mean_FIREemiss_litter_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,26],
                      FIRElitter_litter_gCm2day = output[,,27],
                      mean_FIRElitter_litter_gCm2day = output_mean[,27],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,27],
                      FIREemiss_som_gCm2day = output[,,28],
                      mean_FIREemiss_som_gCm2day = output_mean[,28],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,28],
                      HARVESTextracted_labile_gCm2day = output[,,29],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,29],
                      HARVESTextracted_foliage_gCm2day = output[,,30],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,30],
                      HARVESTextracted_roots_gCm2day = output[,,31],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,31],
                      HARVESTextracted_wood_gCm2day = output[,,32],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,32],
                      HARVESTextracted_litter_gCm2day = output[,,33],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,33],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,33],
                      HARVESTextracted_som_gCm2day = output[,,34],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,34],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,34],
                      HARVESTlitter_labile_gCm2day = output[,,35],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,35],
                      HARVESTlitter_foliage_gCm2day = output[,,36],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,36],
                      HARVESTlitter_roots_gCm2day = output[,,37],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,37],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,37],
                      HARVESTlitter_wood_gCm2day = output[,,38],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,38],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,38],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,39],
                      mean_labile_gCm2 = output_mean[,39],
                      mean_annual_labile_gCm2 = output_annual[,,39],
                      foliage_gCm2 = output[,,40],
                      mean_foliage_gCm2 = output_mean[,40],
                      mean_annual_foliage_gCm2 = output_annual[,,40],
                      roots_gCm2 = output[,,41],
                      mean_roots_gCm2 = output_mean[,41],
                      mean_annual_roots_gCm2 = output_annual[,,41],
                      wood_gCm2 = output[,,42],
                      mean_wood_gCm2 = output_mean[,42],
                      mean_annual_wood_gCm2 = output_annual[,,42],
                      litter_gCm2 = output[,,43],
                      mean_litter_gCm2 = output_mean[,43],
                      mean_annual_litter_gCm2 = output_annual[,,43],
                      som_gCm2 = output[,,44],
                      mean_som_gCm2 = output_mean[,44],
                      mean_annual_som_gCm2 = output_annual[,,44],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,45],
                      mean_ET_kgH2Om2day = output_mean[,45],
                      mean_annual_ET_kgH2Om2day = output_annual[,,45],
                      Etrans_kgH2Om2day = output[,,46],
                      mean_Etrans_kgH2Om2day = output_mean[,46],
                      mean_annual_Etrans_kgH2Om2day = output_annual[,,46],
                      Esoil_kgH2Om2day = output[,,47],
                      mean_Esoil_kgH2Om2day = output_mean[,47],
                      mean_annual_Esoil_kgH2Om2day = output_annual[,,47],
                      Ewetcanopy_kgH2Om2day = output[,,48],
                      mean_Ewetcanopy_kgH2Om2day = output_mean[,48],
                      mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,48],
                      runoff_kgH2Om2day = output[,,49],
                      mean_runoff_kgH2Om2day = output_mean[,49],
                      mean_annual_runoff_kgH2Om2day = output_annual[,,49],
                      underflow_kgH2Om2day = output[,,50],
                      mean_underflow_kgH2Om2day = output_mean[,50],
                      mean_annual_underflow_kgH2Om2day = output_annual[,,50],
                      SurfDrainage_kgH2Om2day = output[,,51],
                      mean_SurfDrainage_kgH2Om2day = output_mean[,51],
                      mean_annual_SurfDrainage_kgH2Om2day = output_annual[,,51],
                      SurfInfiltrated_kgH2Om2day = output[,,52],
                      mean_SurfInfiltrated_kgH2Om2day = output_mean[,52],
                      mean_annual_SurfInfiltrated_kgH2Om2day = output_annual[,,52],
                      Etrans_1st_root_layer_uptake_fraction = output[,,53],
                      mean_Etrans_1st_root_layer_uptake_fraction = output_mean[,53],
                      mean_annual_Etrans_1st_root_layer_uptake_fraction = output_annual[,,53],
                      Etrans_2nd_root_layer_uptake_fraction = output[,,54],
                      mean_Etrans_2nd_root_layer_uptake_fraction = output_mean[,54],
                      mean_annual_Etrans_2nd_root_layer_uptake_fraction = output_annual[,,54],
                      SurfWater_kgH2Om2 = output[,,55],
                      mean_SurfWater_kgH2Om2 = output_mean[,55],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,55],
                      wSWP_MPa = output[,,56],
                      mean_wSWP_MPa = output_mean[,56],
                      mean_annual_wSWP_MPa = output_annual[,,56],
                      snow_kgH2Om2 = output[,,57],
                      mean_snow_kgH2Om2 = output_mean[,57],
                      mean_annual_snow_kgH2Om2 = output_annual[,,57],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,58],
                      mean_lai_m2m2 = output_mean[,58],
                      mean_annual_lai_m2m2 = output_annual[,,58],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,59],
                      mean_gs_demand_supply_ratio = output_mean[,59],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,59],
                      gs_mmolH2Om2s = output[,,60],
                      mean_gs_mmolH2Om2s = output_mean[,60],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,60],
                      APAR_MJm2day = output[,,61],
                      mean_APAR_MJm2day = output_mean[,61],
                      mean_annual_APAR_MJm2day = output_annual[,,61],
                      gb_mmolH2Om2s = output[,,62],
                      mean_gb_mmolH2Om2s = output_mean[,62],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,62],
                      CiCa = output[,,63],
                      mean_CiCa = output_mean[,63],
                      mean_annual_CiCa = output_annual[,,63],
                      # Misc
                      RootDepth_m = output[,,64],
                      mean_RootDepth_m = output_mean[,64],
                      mean_annual_RootDepth_m = output_annual[,,64],
                      leaf_temperature_celcius = output[,,65],
                      mean_leaf_temperature_celcius = output_mean[,65],
                      mean_annual_leaf_temperature_celcius = output_annual[,,65],                      
                      soil_temperature_celcius = output[,,66],
                      mean_soil_temperature_celcius = output_mean[,66],
                      mean_annual_soil_temperature_celcius = output_annual[,,66],   
                      # Plant water stress
                      LWP_MPa = output[,,67],
                      mean_LWP_MPa = output_mean[,67],
                      mean_annual_LWP_MPa = output_annual[,,67],
                      # C allocation diagnositics
                      LabBio_limitation = output[,,68],
                      mean_LabBio_limitation = output_mean[,68],
                      mean_annual_LabBio_limitation = output_annual[,,68],
                      foliage_leafT_limitation = output[,,69],
                      mean_foliage_leafT_limitation = output_mean[,69],
                      mean_annual_foliage_leafT_limitation = output_annual[,,69],
                      roots_leafT_limitation = output[,,70],
                      mean_roots_leafT_limitation = output_mean[,70],
                      mean_annual_roots_leafT_limitation = output_annual[,,70],
                      wood_leafT_limitation = output[,,71],
                      mean_wood_leafT_limitation = output_mean[,71],
                      mean_annual_wood_leafT_limitation = output_annual[,,71],
                      foliage_wSWP_limitation = output[,,72],
                      mean_foliage_wSWP_limitation = output_mean[,72],
                      mean_annual_foliage_wSWP_limitation = output_annual[,,72],
                      roots_wSWP_limitation = output[,,73],
                      mean_roots_wSWP_limitation = output_mean[,73],
                      mean_annual_roots_wSWP_limitation = output_annual[,,73],
                      wood_wSWP_limitation = output[,,74],
                      mean_wood_wSWP_limitation = output_mean[,74],
                      mean_annual_wood_wSWP_limitation = output_annual[,,74],
                      ncce_gCm2day = output[,,75],
                      mean_ncce_gCm2day = output_mean[,75],
                      mean_annual_ncce_gCm2day = output_annual[,,75],
                      canopy_area_scaling_wind = output[,,76],
                      mean_canopy_area_scaling_wind = output_mean[,76],
                      mean_annual_canopy_area_scaling_wind = output_annual[,,76],
                      canopy_area_scaling_light = output[,,77],
                      mean_canopy_area_scaling_light = output_mean[,77],
                      mean_annual_canopy_area_scaling_light = output_annual[,,77],     
                      foliar_growth_limitation = output[,,78],
                      mean_foliar_growth_limitation = output_mean[,78],
                      mean_annual_foliar_growth_limitation = output_annual[,,78], 
                      MTT_foliage_days = output[,,79],
                      mean_MTT_foliage_days = output_mean[,79],
                      mean_annual_MTT_foliage_days = output_annual[,,79],
                      avg_whole_plant_balance_gCm2day = output[,,80],
                      mean_avg_whole_plant_balance_gCm2day = output_mean[,80],
                      mean_annual_foliage_leafT_turnover_limitation = output_annual[,,80],    
                      dummy = output[,,81],
                      mean_dummy = output_mean[,81],
                      mean_annual_dummy = output_annual[,,81],
                      ncce_grow_gCgC = output[,,82],
                      mean_ncce_grow_gCgC = output_mean[,82],
                      mean_annual_ncce_grow_gCgC = output_annual[,,82],
                      ncce_loss_gCgC = output[,,83],
                      mean_ncce_loss_gCgC = output_mean[,83],
                      mean_annual_ncce_loss_gCgC = output_annual[,,83],       
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,output_mean,output_annual,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A3.C1.D2.F2.H2.P1.030") {
      output_dim = 62 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec30",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                             ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                    
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodays=as.integer(dim(met)[1])
                             ,nos_years=as.integer(noyears)
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                             ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                             ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2]))
                             ,residual_waterfrac_in=as.double(c(soil_info[5],soil_info[6],soil_info[6])) 
                             ,porosity_in = as.double(c(soil_info[7],soil_info[8],soil_info[8]))
                             ,pore_size_dist_in = as.double(c(soil_info[9],soil_info[10],soil_info[10]))
                             ,air_entry_in = as.double(c(soil_info[11],soil_info[12],soil_info[12]))
                             ,sat_conductivity_in = as.double(c(soil_info[13],soil_info[14],soil_info[14])))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))      
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      mean_ET_kgH2Om2day = output_mean[,43],
                      mean_annual_ET_kgH2Om2day = output_annual[,,43],
                      Etrans_kgH2Om2day = output[,,44],
                      mean_Etrans_kgH2Om2day = output_mean[,44],
                      mean_annual_Etrans_kgH2Om2day = output_annual[,,44],
                      Esoil_kgH2Om2day = output[,,45],
                      mean_Esoil_kgH2Om2day = output_mean[,45],
                      mean_annual_Esoil_kgH2Om2day = output_annual[,,45],
                      Ewetcanopy_kgH2Om2day = output[,,46],
                      mean_Ewetcanopy_kgH2Om2day = output_mean[,46],
                      mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,46],
                      runoff_kgH2Om2day = output[,,47],
                      mean_runoff_kgH2Om2day = output_mean[,47],
                      mean_annual_runoff_kgH2Om2day = output_annual[,,47],
                      underflow_kgH2Om2day = output[,,48],
                      mean_underflow_kgH2Om2day = output_mean[,48],
                      mean_annual_underflow_kgH2Om2day = output_annual[,,48],
                      SurfDrainage_kgH2Om2day = output[,,49],
                      mean_SurfDrainage_kgH2Om2day = output_mean[,49],
                      mean_annual_SurfDrainage_kgH2Om2day = output_annual[,,49],
                      SurfInfiltrated_kgH2Om2day = output[,,50],
                      mean_SurfInfiltrated_kgH2Om2day = output_mean[,50],
                      mean_annual_SurfInfiltrated_kgH2Om2day = output_annual[,,50],
                      Etrans_1st_root_layer_uptake_fraction = output[,,51],
                      mean_Etrans_1st_root_layer_uptake_fraction = output_mean[,51],
                      mean_annual_Etrans_1st_root_layer_uptake_fraction = output_annual[,,51],
                      Etrans_2nd_root_layer_uptake_fraction = output[,,52],
                      mean_Etrans_2nd_root_layer_uptake_fraction = output_mean[,52],
                      mean_annual_Etrans_2nd_root_layer_uptake_fraction = output_annual[,,52],
                      SurfWater_kgH2Om2 = output[,,53],
                      mean_SurfWater_kgH2Om2 = output_mean[,53],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,53],
                      wSWP_MPa = output[,,54],
                      mean_wSWP_MPa = output_mean[,54],
                      mean_annual_wSWP_MPa = output_annual[,,54],
                      snow_kgH2Om2 = output[,,55],
                      mean_snow_kgH2Om2 = output_mean[,55],
                      mean_annual_snow_kgH2Om2 = output_annual[,,55],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,56],
                      mean_lai_m2m2 = output_mean[,56],
                      mean_annual_lai_m2m2 = output_annual[,,56],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,57],
                      mean_gs_demand_supply_ratio = output_mean[,57],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,57],
                      gs_mmolH2Om2s = output[,,58],
                      mean_gs_mmolH2Om2s = output_mean[,58],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,58],
                      APAR_MJm2day = output[,,59],
                      mean_APAR_MJm2day = output_mean[,59],
                      mean_annual_APAR_MJm2day = output_annual[,,59],
                      gb_mmolH2Om2s = output[,,60],
                      mean_gb_mmolH2Om2s = output_mean[,60],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,60],
                      CiCa = output[,,61],
                      mean_CiCa = output_mean[,61],
                      mean_annual_CiCa = output_annual[,,61],
                      # Misc
                      RootDepth_m = output[,,62],
                      mean_RootDepth_m = output_mean[,62],
                      mean_annual_RootDepth_m = output_annual[,,62],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C1.D2.F2.H2.P2.018") {
      output_dim = 62 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec18",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                              ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                     
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,nos_years=as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                              ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3)))
                              ,residual_waterfrac_in=as.double(c(soil_info[5],soil_info[6],soil_info[6])) 
                              ,porosity_in = as.double(c(soil_info[7],soil_info[8],soil_info[8]))
                              ,pore_size_dist_in = as.double(c(soil_info[9],soil_info[10],soil_info[10]))
                              ,air_entry_in = as.double(c(soil_info[11],soil_info[12],soil_info[12]))
                              ,sat_conductivity_in = as.double(c(soil_info[13],soil_info[14],soil_info[14])))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))      
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      mean_ET_kgH2Om2day = output_mean[,43],
                      mean_annual_ET_kgH2Om2day = output_annual[,,43],
                      Etrans_kgH2Om2day = output[,,44],
                      mean_Etrans_kgH2Om2day = output_mean[,44],
                      mean_annual_Etrans_kgH2Om2day = output_annual[,,44],
                      Esoil_kgH2Om2day = output[,,45],
                      mean_Esoil_kgH2Om2day = output_mean[,45],
                      mean_annual_Esoil_kgH2Om2day = output_annual[,,45],
                      Ewetcanopy_kgH2Om2day = output[,,46],
                      mean_Ewetcanopy_kgH2Om2day = output_mean[,46],
                      mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,46],
                      runoff_kgH2Om2day = output[,,47],
                      mean_runoff_kgH2Om2day = output_mean[,47],
                      mean_annual_runoff_kgH2Om2day = output_annual[,,47],
                      underflow_kgH2Om2day = output[,,48],
                      mean_underflow_kgH2Om2day = output_mean[,48],
                      mean_annual_underflow_kgH2Om2day = output_annual[,,48],
                      SurfDrainage_kgH2Om2day = output[,,49],
                      mean_SurfDrainage_kgH2Om2day = output_mean[,49],
                      mean_annual_SurfDrainage_kgH2Om2day = output_annual[,,49],
                      SurfInfiltrated_kgH2Om2day = output[,,50],
                      mean_SurfInfiltrated_kgH2Om2day = output_mean[,50],
                      mean_annual_SurfInfiltrated_kgH2Om2day = output_annual[,,50],
                      Etrans_1st_root_layer_uptake_fraction = output[,,51],
                      mean_Etrans_1st_root_layer_uptake_fraction = output_mean[,51],
                      mean_annual_Etrans_1st_root_layer_uptake_fraction = output_annual[,,51],
                      Etrans_2nd_root_layer_uptake_fraction = output[,,52],
                      mean_Etrans_2nd_root_layer_uptake_fraction = output_mean[,52],
                      mean_annual_Etrans_2nd_root_layer_uptake_fraction = output_annual[,,52],
                      SurfWater_kgH2Om2 = output[,,53],
                      mean_SurfWater_kgH2Om2 = output_mean[,53],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,53],
                      wSWP_MPa = output[,,54],
                      mean_wSWP_MPa = output_mean[,54],
                      mean_annual_wSWP_MPa = output_annual[,,54],
                      snow_kgH2Om2 = output[,,55],
                      mean_snow_kgH2Om2 = output_mean[,55],
                      mean_annual_snow_kgH2Om2 = output_annual[,,55],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,56],
                      mean_lai_m2m2 = output_mean[,56],
                      mean_annual_lai_m2m2 = output_annual[,,56],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,57],
                      mean_gs_demand_supply_ratio = output_mean[,57],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,57],
                      gs_mmolH2Om2s = output[,,58],
                      mean_gs_mmolH2Om2s = output_mean[,58],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,58],
                      APAR_MJm2day = output[,,59],
                      mean_APAR_MJm2day = output_mean[,59],
                      mean_annual_APAR_MJm2day = output_annual[,,59],
                      gb_mmolH2Om2s = output[,,60],
                      mean_gb_mmolH2Om2s = output_mean[,60],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,60],
                      CiCa = output[,,61],
                      mean_CiCa = output_mean[,61],
                      mean_annual_CiCa = output_annual[,,61],
                      # Misc
                      RootDepth_m = output[,,62],
                      mean_RootDepth_m = output_mean[,62],
                      mean_annual_RootDepth_m = output_annual[,,62],

                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C1.D2.F2.H2.P5.021") {
      output_dim = 52 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec21",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                              ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                     
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,nos_years=as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                              ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3)))
                              ,residual_waterfrac_in=as.double(c(soil_info[5],soil_info[6],soil_info[6])) 
                              ,porosity_in = as.double(c(soil_info[7],soil_info[8],soil_info[8]))
                              ,pore_size_dist_in = as.double(c(soil_info[9],soil_info[10],soil_info[10]))
                              ,air_entry_in = as.double(c(soil_info[11],soil_info[12],soil_info[12]))
                              ,sat_conductivity_in = as.double(c(soil_info[13],soil_info[14],soil_info[14])))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))      
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      mean_ET_kgH2Om2day = output_mean[,43],
                      mean_annual_ET_kgH2Om2day = output_annual[,,43],
                      SurfWater_kgH2Om2 = output[,,44],
                      mean_SurfWater_kgH2Om2 = output_mean[,44],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,44],
                      wSWP_MPa = output[,,45],
                      mean_wSWP_MPa = output_mean[,45],
                      mean_annual_wSWP_MPa = output_annual[,,45],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,46],
                      mean_lai_m2m2 = output_mean[,46],
                      mean_annual_lai_m2m2 = output_annual[,,46],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,47],
                      mean_gs_demand_supply_ratio = output_mean[,47],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,47],
                      gs_mmolH2Om2s = output[,,48],
                      mean_gs_mmolH2Om2s = output_mean[,48],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,48],
                      APAR_MJm2day = output[,,49],
                      mean_APAR_MJm2day = output_mean[,49],
                      mean_annual_APAR_MJm2day = output_annual[,,49],
                      gb_mmolH2Om2s = output[,,50],
                      mean_gb_mmolH2Om2s = output_mean[,50],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,50],
                      CiCa = output[,,51],
                      mean_CiCa = output_mean[,51],
                      mean_annual_CiCa = output_annual[,,51],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C1.D2.F2.H2.P6.022") {
      output_dim = 62 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec22",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                              ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                     
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,nos_years=as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                              ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3)))
                              ,residual_waterfrac_in=as.double(c(soil_info[5],soil_info[6],soil_info[6])) 
                              ,porosity_in = as.double(c(soil_info[7],soil_info[8],soil_info[8]))
                              ,pore_size_dist_in = as.double(c(soil_info[9],soil_info[10],soil_info[10]))
                              ,air_entry_in = as.double(c(soil_info[11],soil_info[12],soil_info[12]))
                              ,sat_conductivity_in = as.double(c(soil_info[13],soil_info[14],soil_info[14])))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))      
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      mean_ET_kgH2Om2day = output_mean[,43],
                      mean_annual_ET_kgH2Om2day = output_annual[,,43],
                      Etrans_kgH2Om2day = output[,,44],
                      mean_Etrans_kgH2Om2day = output_mean[,44],
                      mean_annual_Etrans_kgH2Om2day = output_annual[,,44],
                      Esoil_kgH2Om2day = output[,,45],
                      mean_Esoil_kgH2Om2day = output_mean[,45],
                      mean_annual_Esoil_kgH2Om2day = output_annual[,,45],
                      Ewetcanopy_kgH2Om2day = output[,,46],
                      mean_Ewetcanopy_kgH2Om2day = output_mean[,46],
                      mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,46],
                      runoff_kgH2Om2day = output[,,47],
                      mean_runoff_kgH2Om2day = output_mean[,47],
                      mean_annual_runoff_kgH2Om2day = output_annual[,,47],
                      underflow_kgH2Om2day = output[,,48],
                      mean_underflow_kgH2Om2day = output_mean[,48],
                      mean_annual_underflow_kgH2Om2day = output_annual[,,48],
                      SurfDrainage_kgH2Om2day = output[,,49],
                      mean_SurfDrainage_kgH2Om2day = output_mean[,49],
                      mean_annual_SurfDrainage_kgH2Om2day = output_annual[,,49],
                      SurfInfiltrated_kgH2Om2day = output[,,50],
                      mean_SurfInfiltrated_kgH2Om2day = output_mean[,50],
                      mean_annual_SurfInfiltrated_kgH2Om2day = output_annual[,,50],
                      Etrans_1st_root_layer_uptake_fraction = output[,,51],
                      mean_Etrans_1st_root_layer_uptake_fraction = output_mean[,51],
                      mean_annual_Etrans_1st_root_layer_uptake_fraction = output_annual[,,51],
                      Etrans_2nd_root_layer_uptake_fraction = output[,,52],
                      mean_Etrans_2nd_root_layer_uptake_fraction = output_mean[,52],
                      mean_annual_Etrans_2nd_root_layer_uptake_fraction = output_annual[,,52],
                      SurfWater_kgH2Om2 = output[,,53],
                      mean_SurfWater_kgH2Om2 = output_mean[,53],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,53],
                      wSWP_MPa = output[,,54],
                      mean_wSWP_MPa = output_mean[,54],
                      mean_annual_wSWP_MPa = output_annual[,,54],
                      snow_kgH2Om2 = output[,,55],
                      mean_snow_kgH2Om2 = output_mean[,55],
                      mean_annual_snow_kgH2Om2 = output_annual[,,55],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,56],
                      mean_lai_m2m2 = output_mean[,56],
                      mean_annual_lai_m2m2 = output_annual[,,56],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,57],
                      mean_gs_demand_supply_ratio = output_mean[,57],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,57],
                      gs_mmolH2Om2s = output[,,58],
                      mean_gs_mmolH2Om2s = output_mean[,58],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,58],
                      APAR_MJm2day = output[,,59],
                      mean_APAR_MJm2day = output_mean[,59],
                      mean_annual_APAR_MJm2day = output_annual[,,59],
                      gb_mmolH2Om2s = output[,,60],
                      mean_gb_mmolH2Om2s = output_mean[,60],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,60],
                      CiCa = output[,,61],
                      mean_CiCa = output_mean[,61],
                      mean_annual_CiCa = output_annual[,,61],
                      # Misc
                      RootDepth_m = output[,,62],
                      mean_RootDepth_m = output_mean[,62],
                      mean_annual_RootDepth_m = output_annual[,,62],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C1.D2.F2.H2.P1.R1.05") {
      output_dim = 62 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec5",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                             ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                    
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodays=as.integer(dim(met)[1])
                             ,nos_years=as.integer(noyears)
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                             ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                             ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3)))
                             ,residual_waterfrac_in=as.double(c(soil_info[5],soil_info[6],soil_info[6])) 
                             ,porosity_in = as.double(c(soil_info[7],soil_info[8],soil_info[8]))
                             ,pore_size_dist_in = as.double(c(soil_info[9],soil_info[10],soil_info[10]))
                             ,air_entry_in = as.double(c(soil_info[11],soil_info[12],soil_info[12]))
                             ,sat_conductivity_in = as.double(c(soil_info[13],soil_info[14],soil_info[14])))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))      
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      mean_alloc_labile_gCm2day = output_mean[,8],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      mean_alloc_roots_gCm2day = output_mean[,9],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      mean_alloc_wood_gCm2day = output_mean[,10],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      mean_labile_to_foliage_gCm2day = output_mean[,11],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      mean_foliage_to_litter_gCm2day = output_mean[,12],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      mean_roots_to_litter_gCm2day = output_mean[,13],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      mean_wood_to_litter_gCm2day = output_mean[,14],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      mean_litter_to_som_gCm2day = output_mean[,15],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      mean_FIREemiss_labile_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      mean_FIRElitter_labile_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      mean_FIREemiss_roots_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      mean_FIRElitter_roots_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      mean_FIREemiss_wood_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      mean_FIRElitter_wood_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      mean_FIREemiss_litter_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      mean_FIRElitter_litter_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      mean_FIREemiss_som_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      mean_labile_gCm2 = output_mean[,37],
                      mean_annual_labile_gCm2 = output_annual[,,37],
                      foliage_gCm2 = output[,,38],
                      mean_foliage_gCm2 = output_mean[,38],
                      mean_annual_foliage_gCm2 = output_annual[,,38],
                      roots_gCm2 = output[,,39],
                      mean_roots_gCm2 = output_mean[,39],
                      mean_annual_roots_gCm2 = output_annual[,,39],
                      wood_gCm2 = output[,,40],
                      mean_wood_gCm2 = output_mean[,40],
                      mean_annual_wood_gCm2 = output_annual[,,40],
                      litter_gCm2 = output[,,41],
                      mean_litter_gCm2 = output_mean[,41],
                      mean_annual_litter_gCm2 = output_annual[,,41],
                      som_gCm2 = output[,,42],
                      mean_som_gCm2 = output_mean[,42],
                      mean_annual_som_gCm2 = output_annual[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      mean_ET_kgH2Om2day = output_mean[,43],
                      mean_annual_ET_kgH2Om2day = output_annual[,,43],
                      Etrans_kgH2Om2day = output[,,44],
                      mean_Etrans_kgH2Om2day = output_mean[,44],
                      mean_annual_Etrans_kgH2Om2day = output_annual[,,44],
                      Esoil_kgH2Om2day = output[,,45],
                      mean_Esoil_kgH2Om2day = output_mean[,45],
                      mean_annual_Esoil_kgH2Om2day = output_annual[,,45],
                      Ewetcanopy_kgH2Om2day = output[,,46],
                      mean_Ewetcanopy_kgH2Om2day = output_mean[,46],
                      mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,46],
                      runoff_kgH2Om2day = output[,,47],
                      mean_runoff_kgH2Om2day = output_mean[,47],
                      mean_annual_runoff_kgH2Om2day = output_annual[,,47],
                      underflow_kgH2Om2day = output[,,48],
                      mean_underflow_kgH2Om2day = output_mean[,48],
                      mean_annual_underflow_kgH2Om2day = output_annual[,,48],
                      SurfDrainage_kgH2Om2day = output[,,49],
                      mean_SurfDrainage_kgH2Om2day = output_mean[,49],
                      mean_annual_SurfDrainage_kgH2Om2day = output_annual[,,49],
                      SurfInfiltrated_kgH2Om2day = output[,,50],
                      mean_SurfInfiltrated_kgH2Om2day = output_mean[,50],
                      mean_annual_SurfInfiltrated_kgH2Om2day = output_annual[,,50],
                      Etrans_1st_root_layer_uptake_fraction = output[,,51],
                      mean_Etrans_1st_root_layer_uptake_fraction = output_mean[,51],
                      mean_annual_Etrans_1st_root_layer_uptake_fraction = output_annual[,,51],
                      Etrans_2nd_root_layer_uptake_fraction = output[,,52],
                      mean_Etrans_2nd_root_layer_uptake_fraction = output_mean[,52],
                      mean_annual_Etrans_2nd_root_layer_uptake_fraction = output_annual[,,52],
                      SurfWater_kgH2Om2 = output[,,53],
                      mean_SurfWater_kgH2Om2 = output_mean[,53],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,53],
                      wSWP_MPa = output[,,54],
                      mean_wSWP_MPa = output_mean[,54],
                      mean_annual_wSWP_MPa = output_annual[,,54],
                      snow_kgH2Om2 = output[,,55],
                      mean_snow_kgH2Om2 = output_mean[,55],
                      mean_annual_snow_kgH2Om2 = output_annual[,,55],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,56],
                      mean_lai_m2m2 = output_mean[,56],
                      mean_annual_lai_m2m2 = output_annual[,,56],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,57],
                      mean_gs_demand_supply_ratio = output_mean[,57],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,57],
                      gs_mmolH2Om2s = output[,,58],
                      mean_gs_mmolH2Om2s = output_mean[,58],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,58],
                      APAR_MJm2day = output[,,59],
                      mean_APAR_MJm2day = output_mean[,59],
                      mean_annual_APAR_MJm2day = output_annual[,,59],
                      gb_mmolH2Om2s = output[,,60],
                      mean_gb_mmolH2Om2s = output_mean[,60],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,60],
                      CiCa = output[,,61],
                      mean_CiCa = output_mean[,61],
                      mean_annual_CiCa = output_annual[,,61],
                      # Misc
                      RootDepth_m = output[,,62],
                      mean_RootDepth_m = output_mean[,62],
                      mean_annual_RootDepth_m = output_annual[,,62],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_som_years = MTT_years[,6],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_som_gCm2 = SS_gCm2[,6])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                     states_all$alloc_foliage_gCm2day +
                     states_all$alloc_roots_gCm2day +
                     states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                     apply(states_all$alloc_roots_gCm2day,1,mean),
                     apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P1.R1.006") {
      output_dim = 68 ; MTT_dim = 7 ; SS_dim = 7
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec6",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                             ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                    
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodays=as.integer(dim(met)[1])
                             ,nos_years=as.integer(noyears)
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                             ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                             ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2]))
                             ,residual_waterfrac_in=as.double(c(soil_info[5],soil_info[6],soil_info[6])) 
                             ,porosity_in = as.double(c(soil_info[7],soil_info[8],soil_info[8]))
                             ,pore_size_dist_in = as.double(c(soil_info[9],soil_info[10],soil_info[10]))
                             ,air_entry_in = as.double(c(soil_info[11],soil_info[12],soil_info[12]))
                             ,sat_conductivity_in = as.double(c(soil_info[13],soil_info[14],soil_info[14])))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))      
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      rhet_woodlitter_gCm2day = output[,,5],
                      mean_rhet_woodlitter_gCm2day = output_mean[,5],
                      mean_annual_rhet_woodlitter_gCm2day = output_annual[,,5],
                      fire_gCm2day = output[,,6],
                      mean_fire_gCm2day = output_mean[,6],
                      mean_annual_fire_gCm2day = output_annual[,,6],
                      harvest_gCm2day = output[,,7],
                      mean_harvest_gCm2day = output_mean[,7],
                      mean_annual_harvest_gCm2day = output_annual[,,7],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,8],
                      mean_alloc_foliage_gCm2day = output_mean[,8],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,8],
                      alloc_labile_gCm2day = output[,,9],
                      mean_alloc_labile_gCm2day = output_mean[,9],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,9],
                      alloc_roots_gCm2day = output[,,10],
                      mean_alloc_roots_gCm2day = output_mean[,10],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,10],
                      alloc_wood_gCm2day = output[,,11],
                      mean_alloc_wood_gCm2day = output_mean[,11],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,11],
                      labile_to_foliage_gCm2day = output[,,12],
                      mean_labile_to_foliage_gCm2day = output_mean[,12],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,12],
                      foliage_to_litter_gCm2day = output[,,13],
                      mean_foliage_to_litter_gCm2day = output_mean[,13],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,13],
                      roots_to_litter_gCm2day = output[,,14],
                      mean_roots_to_litter_gCm2day = output_mean[,14],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,14],
                      wood_to_litter_gCm2day = output[,,15],
                      mean_wood_to_litter_gCm2day = output_mean[,15],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,15],
                      litter_to_som_gCm2day = output[,,16],
                      mean_litter_to_som_gCm2day = output_mean[,16],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,16],
                      woodlitter_to_som_gCm2day = output[,,17],
                      mean_woodlitter_to_som_gCm2day = output_mean[,17],
                      mean_annual_woodlitter_to_som_gCm2day = output_annual[,,17],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,18],
                      mean_FIREemiss_labile_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,18],
                      FIRElitter_labile_gCm2day = output[,,19],
                      mean_FIRElitter_labile_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,19],
                      FIREemiss_foliage_gCm2day = output[,,20],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,20],
                      FIRElitter_foliage_gCm2day = output[,,21],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,21],
                      FIREemiss_roots_gCm2day = output[,,22],
                      mean_FIREemiss_roots_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,22],
                      FIRElitter_roots_gCm2day = output[,,23],
                      mean_FIRElitter_roots_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,23],
                      FIREemiss_wood_gCm2day = output[,,24],
                      mean_FIREemiss_wood_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,24],
                      FIRElitter_wood_gCm2day = output[,,25],
                      mean_FIRElitter_wood_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,25],
                      FIREemiss_litter_gCm2day = output[,,26],
                      mean_FIREemiss_litter_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,26],
                      FIRElitter_litter_gCm2day = output[,,27],
                      mean_FIRElitter_litter_gCm2day = output_mean[,27],
                      mean_annual_FIRElitter_litter_gCm2day = output[,,27],
                      FIREemiss_woodlitter_gCm2day = output[,,28],
                      mean_FIREemiss_woodlitter_gCm2day = output_mean[,28],
                      mean_annual_FIREemiss_woodlitter_gCm2day = output_annual[,,28],
                      FIRElitter_woodlitter_gCm2day = output[,,29],
                      mean_FIRElitter_woodlitter_gCm2day = output_mean[,29],
                      mean_annual_FIRElitter_woodlitter_gCm2day = output_annual[,,29],
                      FIREemiss_som_gCm2day = output[,,30],
                      mean_FIREemiss_som_gCm2day = output_mean[,30],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,30],
                      HARVESTextracted_labile_gCm2day = output[,,31],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,31],
                      HARVESTextracted_foliage_gCm2day = output[,,32],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,32],
                      HARVESTextracted_roots_gCm2day = output[,,33],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,33],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,33],
                      HARVESTextracted_wood_gCm2day = output[,,34],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,34],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,34],
                      HARVESTextracted_litter_gCm2day = output[,,35],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,35],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,35],
                      HARVESTextracted_woodlitter_gCm2day = output[,,36],
                      mean_HARVESTextracted_woodlitter_gCm2day = output_mean[,36],
                      mean_annual_HARVESTextracted_woodlitter_gCm2day = output_annual[,,36],
                      HARVESTextracted_som_gCm2day = output[,,37],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,37],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,37],
                      HARVESTlitter_labile_gCm2day = output[,,38],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,38],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,38],
                      HARVESTlitter_foliage_gCm2day = output[,,39],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,39],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,39],
                      HARVESTlitter_roots_gCm2day = output[,,40],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,40],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,40],
                      HARVESTlitter_wood_gCm2day = output[,,41],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,41],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,41],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,42],
                      mean_labile_gCm2 = output_mean[,42],
                      mean_annual_labile_gCm2 = output_annual[,,42],
                      foliage_gCm2 = output[,,43],
                      mean_foliage_gCm2 = output_mean[,43],
                      mean_annual_foliage_gCm2 = output_annual[,,43],
                      roots_gCm2 = output[,,44],
                      mean_roots_gCm2 = output_mean[,44],
                      mean_annual_roots_gCm2 = output_annual[,,44],
                      wood_gCm2 = output[,,45],
                      mean_wood_gCm2 = output_mean[,45],
                      mean_annual_wood_gCm2 = output_annual[,,45],
                      litter_gCm2 = output[,,46],
                      mean_litter_gCm2 = output_mean[,46],
                      mean_annual_litter_gCm2 = output_annual[,,46],
                      woodlitter_gCm2 = output[,,47],
                      mean_woodlitter_gCm2 = output_mean[,47],
                      mean_annual_woodlitter_gCm2 = output_annual[,,47],
                      som_gCm2 = output[,,48],
                      mean_som_gCm2 = output_mean[,48],
                      mean_annual_som_gCm2 = output_annual[,,48],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,49],
                      mean_ET_kgH2Om2day = output_mean[,49],
                      mean_annual_ET_kgH2Om2day = output_annual[,,49],
                      Etrans_kgH2Om2day = output[,,50],
                      mean_Etrans_kgH2Om2day = output_mean[,50],
                      mean_annual_Etrans_kgH2Om2day = output_annual[,,50],
                      Esoil_kgH2Om2day = output[,,51],
                      mean_Esoil_kgH2Om2day = output_mean[,51],
                      mean_annual_Esoil_kgH2Om2day = output_annual[,,51],
                      Ewetcanopy_kgH2Om2day = output[,,52],
                      mean_Ewetcanopy_kgH2Om2day = output_mean[,52],
                      mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,52],
                      runoff_kgH2Om2day = output[,,53],
                      mean_runoff_kgH2Om2day = output_mean[,53],
                      mean_annual_runoff_kgH2Om2day = output_annual[,,53],
                      underflow_kgH2Om2day = output[,,54],
                      mean_underflow_kgH2Om2day = output_mean[,54],
                      mean_annual_underflow_kgH2Om2day = output_annual[,,54],
                      SurfDrainage_kgH2Om2day = output[,,55],
                      mean_SurfDrainage_kgH2Om2day = output_mean[,55],
                      mean_annual_SurfDrainage_kgH2Om2day = output_annual[,,55],
                      SurfInfiltrated_kgH2Om2day = output[,,56],
                      mean_SurfInfiltrated_kgH2Om2day = output_mean[,56],
                      mean_annual_SurfInfiltrated_kgH2Om2day = output_annual[,,56],
                      Etrans_1st_root_layer_uptake_fraction = output[,,57],
                      mean_Etrans_1st_root_layer_uptake_fraction = output_mean[,57],
                      mean_annual_Etrans_1st_root_layer_uptake_fraction = output_annual[,,57],
                      Etrans_2nd_root_layer_uptake_fraction = output[,,58],
                      mean_Etrans_2nd_root_layer_uptake_fraction = output_mean[,58],
                      mean_annual_Etrans_2nd_root_layer_uptake_fraction = output_annual[,,58],
                      SurfWater_kgH2Om2 = output[,,59],
                      mean_SurfWater_kgH2Om2 = output_mean[,59],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,59],
                      wSWP_MPa = output[,,60],
                      mean_wSWP_MPa = output_mean[,60],
                      mean_annual_wSWP_MPa = output_annual[,,60],
                      snow_kgH2Om2 = output[,,61],
                      mean_snow_kgH2Om2 = output_mean[,61],
                      mean_annual_snow_kgH2Om2 = output_annual[,,61],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,62],
                      mean_lai_m2m2 = output_mean[,62],
                      mean_annual_lai_m2m2 = output_annual[,,62],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,63],
                      mean_gs_demand_supply_ratio = output_mean[,63],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,63],
                      gs_mmolH2Om2s = output[,,64],
                      mean_gs_mmolH2Om2s = output_mean[,64],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,64],
                      APAR_MJm2day = output[,,65],
                      mean_APAR_MJm2day = output_mean[,65],
                      mean_annual_APAR_MJm2day = output_annual[,,65],
                      gb_mmolH2Om2s = output[,,66],
                      mean_gb_mmolH2Om2s = output_mean[,66],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,66],
                      CiCa = output[,,67],
                      mean_CiCa = output_mean[,67],
                      mean_annual_CiCa = output_annual[,,67],
                      # Misc
                      RootDepth_m = output[,,68],
                      mean_RootDepth_m = output_mean[,68],
                      mean_annual_RootDepth_m = output_annual[,,68],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_woodlitter_years = MTT_years[,6],
                      MTT_som_years = MTT_years[,7],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_woodlitter_gCm2 = SS_gCm2[,6],
                      SS_som_gCm2 = SS_gCm2[,7])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P2.R1.007") {
      output_dim = 68 ; MTT_dim = 7 ; SS_dim = 7
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec7",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                             ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                    
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodays=as.integer(dim(met)[1])
                             ,nos_years=as.integer(noyears)
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                             ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                             ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2]))
                             ,residual_waterfrac_in=as.double(c(soil_info[5],soil_info[6],soil_info[6])) 
                             ,porosity_in = as.double(c(soil_info[7],soil_info[8],soil_info[8]))
                             ,pore_size_dist_in = as.double(c(soil_info[9],soil_info[10],soil_info[10]))
                             ,air_entry_in = as.double(c(soil_info[11],soil_info[12],soil_info[12]))
                             ,sat_conductivity_in = as.double(c(soil_info[13],soil_info[14],soil_info[14])))
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))            
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      rhet_woodlitter_gCm2day = output[,,5],
                      mean_rhet_woodlitter_gCm2day = output_mean[,5],
                      mean_annual_rhet_woodlitter_gCm2day = output_annual[,,5],
                      fire_gCm2day = output[,,6],
                      mean_fire_gCm2day = output_mean[,6],
                      mean_annual_fire_gCm2day = output_annual[,,6],
                      harvest_gCm2day = output[,,7],
                      mean_harvest_gCm2day = output_mean[,7],
                      mean_annual_harvest_gCm2day = output_annual[,,7],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,8],
                      mean_alloc_foliage_gCm2day = output_mean[,8],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,8],
                      alloc_labile_gCm2day = output[,,9],
                      mean_alloc_labile_gCm2day = output_mean[,9],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,9],
                      alloc_roots_gCm2day = output[,,10],
                      mean_alloc_roots_gCm2day = output_mean[,10],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,10],
                      alloc_wood_gCm2day = output[,,11],
                      mean_alloc_wood_gCm2day = output_mean[,11],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,11],
                      labile_to_foliage_gCm2day = output[,,12],
                      mean_labile_to_foliage_gCm2day = output_mean[,12],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,12],
                      foliage_to_litter_gCm2day = output[,,13],
                      mean_foliage_to_litter_gCm2day = output_mean[,13],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,13],
                      roots_to_litter_gCm2day = output[,,14],
                      mean_roots_to_litter_gCm2day = output_mean[,14],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,14],
                      wood_to_litter_gCm2day = output[,,15],
                      mean_wood_to_litter_gCm2day = output_mean[,15],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,15],
                      litter_to_som_gCm2day = output[,,16],
                      mean_litter_to_som_gCm2day = output_mean[,16],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,16],
                      woodlitter_to_som_gCm2day = output[,,17],
                      mean_woodlitter_to_som_gCm2day = output_mean[,17],
                      mean_annual_woodlitter_to_som_gCm2day = output_annual[,,17],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,18],
                      mean_FIREemiss_labile_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,18],
                      FIRElitter_labile_gCm2day = output[,,19],
                      mean_FIRElitter_labile_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,19],
                      FIREemiss_foliage_gCm2day = output[,,20],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,20],
                      FIRElitter_foliage_gCm2day = output[,,21],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,21],
                      FIREemiss_roots_gCm2day = output[,,22],
                      mean_FIREemiss_roots_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,22],
                      FIRElitter_roots_gCm2day = output[,,23],
                      mean_FIRElitter_roots_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,23],
                      FIREemiss_wood_gCm2day = output[,,24],
                      mean_FIREemiss_wood_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,24],
                      FIRElitter_wood_gCm2day = output[,,25],
                      mean_FIRElitter_wood_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,25],
                      FIREemiss_litter_gCm2day = output[,,26],
                      mean_FIREemiss_litter_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,26],
                      FIRElitter_litter_gCm2day = output[,,27],
                      mean_FIRElitter_litter_gCm2day = output_mean[,27],
                      mean_annual_FIRElitter_litter_gCm2day = output[,,27],
                      FIREemiss_woodlitter_gCm2day = output[,,28],
                      mean_FIREemiss_woodlitter_gCm2day = output_mean[,28],
                      mean_annual_FIREemiss_woodlitter_gCm2day = output_annual[,,28],
                      FIRElitter_woodlitter_gCm2day = output[,,29],
                      mean_FIRElitter_woodlitter_gCm2day = output_mean[,29],
                      mean_annual_FIRElitter_woodlitter_gCm2day = output_annual[,,29],
                      FIREemiss_som_gCm2day = output[,,30],
                      mean_FIREemiss_som_gCm2day = output_mean[,30],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,30],
                      HARVESTextracted_labile_gCm2day = output[,,31],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,31],
                      HARVESTextracted_foliage_gCm2day = output[,,32],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,32],
                      HARVESTextracted_roots_gCm2day = output[,,33],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,33],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,33],
                      HARVESTextracted_wood_gCm2day = output[,,34],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,34],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,34],
                      HARVESTextracted_litter_gCm2day = output[,,35],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,35],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,35],
                      HARVESTextracted_woodlitter_gCm2day = output[,,36],
                      mean_HARVESTextracted_woodlitter_gCm2day = output_mean[,36],
                      mean_annual_HARVESTextracted_woodlitter_gCm2day = output_annual[,,36],
                      HARVESTextracted_som_gCm2day = output[,,37],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,37],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,37],
                      HARVESTlitter_labile_gCm2day = output[,,38],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,38],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,38],
                      HARVESTlitter_foliage_gCm2day = output[,,39],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,39],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,39],
                      HARVESTlitter_roots_gCm2day = output[,,40],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,40],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,40],
                      HARVESTlitter_wood_gCm2day = output[,,41],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,41],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,41],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,42],
                      mean_labile_gCm2 = output_mean[,42],
                      mean_annual_labile_gCm2 = output_annual[,,42],
                      foliage_gCm2 = output[,,43],
                      mean_foliage_gCm2 = output_mean[,43],
                      mean_annual_foliage_gCm2 = output_annual[,,43],
                      roots_gCm2 = output[,,44],
                      mean_roots_gCm2 = output_mean[,44],
                      mean_annual_roots_gCm2 = output_annual[,,44],
                      wood_gCm2 = output[,,45],
                      mean_wood_gCm2 = output_mean[,45],
                      mean_annual_wood_gCm2 = output_annual[,,45],
                      litter_gCm2 = output[,,46],
                      mean_litter_gCm2 = output_mean[,46],
                      mean_annual_litter_gCm2 = output_annual[,,46],
                      woodlitter_gCm2 = output[,,47],
                      mean_woodlitter_gCm2 = output_mean[,47],
                      mean_annual_woodlitter_gCm2 = output_annual[,,47],
                      som_gCm2 = output[,,48],
                      mean_som_gCm2 = output_mean[,48],
                      mean_annual_som_gCm2 = output_annual[,,48],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,49],
                      mean_ET_kgH2Om2day = output_mean[,49],
                      mean_annual_ET_kgH2Om2day = output_annual[,,49],
                      Etrans_kgH2Om2day = output[,,50],
                      mean_Etrans_kgH2Om2day = output_mean[,50],
                      mean_annual_Etrans_kgH2Om2day = output_annual[,,50],
                      Esoil_kgH2Om2day = output[,,51],
                      mean_Esoil_kgH2Om2day = output_mean[,51],
                      mean_annual_Esoil_kgH2Om2day = output_annual[,,51],
                      Ewetcanopy_kgH2Om2day = output[,,52],
                      mean_Ewetcanopy_kgH2Om2day = output_mean[,52],
                      mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,52],
                      runoff_kgH2Om2day = output[,,53],
                      mean_runoff_kgH2Om2day = output_mean[,53],
                      mean_annual_runoff_kgH2Om2day = output_annual[,,53],
                      underflow_kgH2Om2day = output[,,54],
                      mean_underflow_kgH2Om2day = output_mean[,54],
                      mean_annual_underflow_kgH2Om2day = output_annual[,,54],
                      SurfDrainage_kgH2Om2day = output[,,55],
                      mean_SurfDrainage_kgH2Om2day = output_mean[,55],
                      mean_annual_SurfDrainage_kgH2Om2day = output_annual[,,55],
                      SurfInfiltrated_kgH2Om2day = output[,,56],
                      mean_SurfInfiltrated_kgH2Om2day = output_mean[,56],
                      mean_annual_SurfInfiltrated_kgH2Om2day = output_annual[,,56],
                      Etrans_1st_root_layer_uptake_fraction = output[,,57],
                      mean_Etrans_1st_root_layer_uptake_fraction = output_mean[,57],
                      mean_annual_Etrans_1st_root_layer_uptake_fraction = output_annual[,,57],
                      Etrans_2nd_root_layer_uptake_fraction = output[,,58],
                      mean_Etrans_2nd_root_layer_uptake_fraction = output_mean[,58],
                      mean_annual_Etrans_2nd_root_layer_uptake_fraction = output_annual[,,58],
                      SurfWater_kgH2Om2 = output[,,59],
                      mean_SurfWater_kgH2Om2 = output_mean[,59],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,59],
                      wSWP_MPa = output[,,60],
                      mean_wSWP_MPa = output_mean[,60],
                      mean_annual_wSWP_MPa = output_annual[,,60],
                      snow_kgH2Om2 = output[,,61],
                      mean_snow_kgH2Om2 = output_mean[,61],
                      mean_annual_snow_kgH2Om2 = output_annual[,,61],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,62],
                      mean_lai_m2m2 = output_mean[,62],
                      mean_annual_lai_m2m2 = output_annual[,,62],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,63],
                      mean_gs_demand_supply_ratio = output_mean[,63],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,63],
                      gs_mmolH2Om2s = output[,,64],
                      mean_gs_mmolH2Om2s = output_mean[,64],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,64],
                      APAR_MJm2day = output[,,65],
                      mean_APAR_MJm2day = output_mean[,65],
                      mean_annual_APAR_MJm2day = output_annual[,,65],
                      gb_mmolH2Om2s = output[,,66],
                      mean_gb_mmolH2Om2s = output_mean[,66],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,66],
                      CiCa = output[,,67],
                      mean_CiCa = output_mean[,67],
                      mean_annual_CiCa = output_annual[,,67],
                      # Misc
                      RootDepth_m = output[,,68],
                      mean_RootDepth_m = output_mean[,68],
                      mean_annual_RootDepth_m = output_annual[,,68],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_woodlitter_years = MTT_years[,6],
                      MTT_som_years = MTT_years[,7],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_woodlitter_gCm2 = SS_gCm2[,6],
                      SS_som_gCm2 = SS_gCm2[,7])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P2.R3.019") {
      output_dim = 68 ; MTT_dim = 7 ; SS_dim = 7
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec19",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                              ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                     
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,nos_years=as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                              ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2]))
                              ,residual_waterfrac_in=as.double(c(soil_info[5],soil_info[6],soil_info[6])) 
                              ,porosity_in = as.double(c(soil_info[7],soil_info[8],soil_info[8]))
                              ,pore_size_dist_in = as.double(c(soil_info[9],soil_info[10],soil_info[10]))
                              ,air_entry_in = as.double(c(soil_info[11],soil_info[12],soil_info[12]))
                              ,sat_conductivity_in = as.double(c(soil_info[13],soil_info[14],soil_info[14])))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))            
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      rhet_woodlitter_gCm2day = output[,,5],
                      mean_rhet_woodlitter_gCm2day = output_mean[,5],
                      mean_annual_rhet_woodlitter_gCm2day = output_annual[,,5],
                      fire_gCm2day = output[,,6],
                      mean_fire_gCm2day = output_mean[,6],
                      mean_annual_fire_gCm2day = output_annual[,,6],
                      harvest_gCm2day = output[,,7],
                      mean_harvest_gCm2day = output_mean[,7],
                      mean_annual_harvest_gCm2day = output_annual[,,7],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,8],
                      mean_alloc_foliage_gCm2day = output_mean[,8],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,8],
                      alloc_labile_gCm2day = output[,,9],
                      mean_alloc_labile_gCm2day = output_mean[,9],
                      mean_annual_alloc_labile_gCm2day = output_annual[,,9],
                      alloc_roots_gCm2day = output[,,10],
                      mean_alloc_roots_gCm2day = output_mean[,10],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,10],
                      alloc_wood_gCm2day = output[,,11],
                      mean_alloc_wood_gCm2day = output_mean[,11],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,11],
                      labile_to_foliage_gCm2day = output[,,12],
                      mean_labile_to_foliage_gCm2day = output_mean[,12],
                      mean_annual_labile_to_foliage_gCm2day = output_annual[,,12],
                      foliage_to_litter_gCm2day = output[,,13],
                      mean_foliage_to_litter_gCm2day = output_mean[,13],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,13],
                      roots_to_litter_gCm2day = output[,,14],
                      mean_roots_to_litter_gCm2day = output_mean[,14],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,14],
                      wood_to_litter_gCm2day = output[,,15],
                      mean_wood_to_litter_gCm2day = output_mean[,15],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,15],
                      litter_to_som_gCm2day = output[,,16],
                      mean_litter_to_som_gCm2day = output_mean[,16],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,16],
                      woodlitter_to_som_gCm2day = output[,,17],
                      mean_woodlitter_to_som_gCm2day = output_mean[,17],
                      mean_annual_woodlitter_to_som_gCm2day = output_annual[,,17],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,18],
                      mean_FIREemiss_labile_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,18],
                      FIRElitter_labile_gCm2day = output[,,19],
                      mean_FIRElitter_labile_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,19],
                      FIREemiss_foliage_gCm2day = output[,,20],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,20],
                      FIRElitter_foliage_gCm2day = output[,,21],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,21],
                      FIREemiss_roots_gCm2day = output[,,22],
                      mean_FIREemiss_roots_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,22],
                      FIRElitter_roots_gCm2day = output[,,23],
                      mean_FIRElitter_roots_gCm2day = output_mean[,23],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,23],
                      FIREemiss_wood_gCm2day = output[,,24],
                      mean_FIREemiss_wood_gCm2day = output_mean[,24],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,24],
                      FIRElitter_wood_gCm2day = output[,,25],
                      mean_FIRElitter_wood_gCm2day = output_mean[,25],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,25],
                      FIREemiss_litter_gCm2day = output[,,26],
                      mean_FIREemiss_litter_gCm2day = output_mean[,26],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,26],
                      FIRElitter_litter_gCm2day = output[,,27],
                      mean_FIRElitter_litter_gCm2day = output_mean[,27],
                      mean_annual_FIRElitter_litter_gCm2day = output[,,27],
                      FIREemiss_woodlitter_gCm2day = output[,,28],
                      mean_FIREemiss_woodlitter_gCm2day = output_mean[,28],
                      mean_annual_FIREemiss_woodlitter_gCm2day = output_annual[,,28],
                      FIRElitter_woodlitter_gCm2day = output[,,29],
                      mean_FIRElitter_woodlitter_gCm2day = output_mean[,29],
                      mean_annual_FIRElitter_woodlitter_gCm2day = output_annual[,,29],
                      FIREemiss_som_gCm2day = output[,,30],
                      mean_FIREemiss_som_gCm2day = output_mean[,30],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,30],
                      HARVESTextracted_labile_gCm2day = output[,,31],
                      mean_HARVESTextracted_labile_gCm2day = output_mean[,31],
                      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,31],
                      HARVESTextracted_foliage_gCm2day = output[,,32],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,32],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,32],
                      HARVESTextracted_roots_gCm2day = output[,,33],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,33],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,33],
                      HARVESTextracted_wood_gCm2day = output[,,34],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,34],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,34],
                      HARVESTextracted_litter_gCm2day = output[,,35],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,35],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,35],
                      HARVESTextracted_woodlitter_gCm2day = output[,,36],
                      mean_HARVESTextracted_woodlitter_gCm2day = output_mean[,36],
                      mean_annual_HARVESTextracted_woodlitter_gCm2day = output_annual[,,36],
                      HARVESTextracted_som_gCm2day = output[,,37],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,37],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,37],
                      HARVESTlitter_labile_gCm2day = output[,,38],
                      mean_HARVESTlitter_labile_gCm2day = output_mean[,38],
                      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,38],
                      HARVESTlitter_foliage_gCm2day = output[,,39],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,39],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,39],
                      HARVESTlitter_roots_gCm2day = output[,,40],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,40],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,40],
                      HARVESTlitter_wood_gCm2day = output[,,41],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,41],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,41],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,42],
                      mean_labile_gCm2 = output_mean[,42],
                      mean_annual_labile_gCm2 = output_annual[,,42],
                      foliage_gCm2 = output[,,43],
                      mean_foliage_gCm2 = output_mean[,43],
                      mean_annual_foliage_gCm2 = output_annual[,,43],
                      roots_gCm2 = output[,,44],
                      mean_roots_gCm2 = output_mean[,44],
                      mean_annual_roots_gCm2 = output_annual[,,44],
                      wood_gCm2 = output[,,45],
                      mean_wood_gCm2 = output_mean[,45],
                      mean_annual_wood_gCm2 = output_annual[,,45],
                      litter_gCm2 = output[,,46],
                      mean_litter_gCm2 = output_mean[,46],
                      mean_annual_litter_gCm2 = output_annual[,,46],
                      woodlitter_gCm2 = output[,,47],
                      mean_woodlitter_gCm2 = output_mean[,47],
                      mean_annual_woodlitter_gCm2 = output_annual[,,47],
                      som_gCm2 = output[,,48],
                      mean_som_gCm2 = output_mean[,48],
                      mean_annual_som_gCm2 = output_annual[,,48],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,49],
                      mean_ET_kgH2Om2day = output_mean[,49],
                      mean_annual_ET_kgH2Om2day = output_annual[,,49],
                      Etrans_kgH2Om2day = output[,,50],
                      mean_Etrans_kgH2Om2day = output_mean[,50],
                      mean_annual_Etrans_kgH2Om2day = output_annual[,,50],
                      Esoil_kgH2Om2day = output[,,51],
                      mean_Esoil_kgH2Om2day = output_mean[,51],
                      mean_annual_Esoil_kgH2Om2day = output_annual[,,51],
                      Ewetcanopy_kgH2Om2day = output[,,52],
                      mean_Ewetcanopy_kgH2Om2day = output_mean[,52],
                      mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,52],
                      runoff_kgH2Om2day = output[,,53],
                      mean_runoff_kgH2Om2day = output_mean[,53],
                      mean_annual_runoff_kgH2Om2day = output_annual[,,53],
                      underflow_kgH2Om2day = output[,,54],
                      mean_underflow_kgH2Om2day = output_mean[,54],
                      mean_annual_underflow_kgH2Om2day = output_annual[,,54],
                      SurfDrainage_kgH2Om2day = output[,,55],
                      mean_SurfDrainage_kgH2Om2day = output_mean[,55],
                      mean_annual_SurfDrainage_kgH2Om2day = output_annual[,,55],
                      SurfInfiltrated_kgH2Om2day = output[,,56],
                      mean_SurfInfiltrated_kgH2Om2day = output_mean[,56],
                      mean_annual_SurfInfiltrated_kgH2Om2day = output_annual[,,56],
                      Etrans_1st_root_layer_uptake_fraction = output[,,57],
                      mean_Etrans_1st_root_layer_uptake_fraction = output_mean[,57],
                      mean_annual_Etrans_1st_root_layer_uptake_fraction = output_annual[,,57],
                      Etrans_2nd_root_layer_uptake_fraction = output[,,58],
                      mean_Etrans_2nd_root_layer_uptake_fraction = output_mean[,58],
                      mean_annual_Etrans_2nd_root_layer_uptake_fraction = output_annual[,,58],
                      SurfWater_kgH2Om2 = output[,,59],
                      mean_SurfWater_kgH2Om2 = output_mean[,59],
                      mean_annual_SurfWater_kgH2Om2 = output_annual[,,59],
                      wSWP_MPa = output[,,60],
                      mean_wSWP_MPa = output_mean[,60],
                      mean_annual_wSWP_MPa = output_annual[,,60],
                      snow_kgH2Om2 = output[,,61],
                      mean_snow_kgH2Om2 = output_mean[,61],
                      mean_annual_snow_kgH2Om2 = output_annual[,,61],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,62],
                      mean_lai_m2m2 = output_mean[,62],
                      mean_annual_lai_m2m2 = output_annual[,,62],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,63],
                      mean_gs_demand_supply_ratio = output_mean[,63],
                      mean_annual_gs_demand_supply_ratio = output_annual[,,63],
                      gs_mmolH2Om2s = output[,,64],
                      mean_gs_mmolH2Om2s = output_mean[,64],
                      mean_annual_gs_mmolH2Om2s = output_annual[,,64],
                      APAR_MJm2day = output[,,65],
                      mean_APAR_MJm2day = output_mean[,65],
                      mean_annual_APAR_MJm2day = output_annual[,,65],
                      gb_mmolH2Om2s = output[,,66],
                      mean_gb_mmolH2Om2s = output_mean[,66],
                      mean_annual_gb_mmolH2Om2s = output_annual[,,66],
                      CiCa = output[,,67],
                      mean_CiCa = output_mean[,67],
                      mean_annual_CiCa = output_annual[,,67],
                      # Misc
                      RootDepth_m = output[,,68],
                      mean_RootDepth_m = output_mean[,68],
                      mean_annual_RootDepth_m = output_annual[,,68],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_labile_years = MTT_years[,1],
                      MTT_foliage_years = MTT_years[,2],
                      MTT_roots_years = MTT_years[,3],
                      MTT_wood_years = MTT_years[,4],
                      MTT_litter_years = MTT_years[,5],
                      MTT_woodlitter_years = MTT_years[,6],
                      MTT_som_years = MTT_years[,7],
                      # Steady state estimates
                      SS_labile_gCm2 = SS_gCm2[,1],
                      SS_foliage_gCm2 = SS_gCm2[,2],
                      SS_roots_gCm2 = SS_gCm2[,3],
                      SS_wood_gCm2 = SS_gCm2[,4],
                      SS_litter_gCm2 = SS_gCm2[,5],
                      SS_woodlitter_gCm2 = SS_gCm2[,6],
                      SS_som_gCm2 = SS_gCm2[,7])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Tidy up variables
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.C5.D1.F2.P1.013") {
    output_dim = 22 ; MTT_dim = 4 ; SS_dim = 4
    dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
    tmp=.Fortran( "rdalec13",output_dim=as.integer(output_dim)
                            ,MTT_dim=as.integer(MTT_dim)
                            ,SS_dim = as.integer(SS_dim)
                            ,met=as.double(t(met))
                            ,pars=as.double(pars_in)
                            ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                            ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                            ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                            ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                            ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                   
                            ,lat=as.double(lat)
                            ,nopars=as.integer(PROJECT$model$nopars[site])
                            ,nomet=as.integer(dim(met)[2])
                            ,nofluxes=as.integer(PROJECT$model$nofluxes[site])
                            ,nopools=as.integer(PROJECT$model$nopools[site])
                            ,nodays=as.integer(dim(met)[1])
                            ,nos_years=as.integer(noyears)
                            ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1]))))
                            ,nos_iter=as.integer(nos_iter) )
    output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
    output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
    output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))          
    MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
    SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
    # Unload the current dalec shared object
    dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
    rm(tmp) ; gc()
    # create output object
    states_all=list(# Ecosystem fluxes
                    gpp_gCm2day = output[,,1],
                    mean_gpp_gCm2day = output_mean[,1],
                    mean_annual_gpp_gCm2day = output_annual[,,1],
                    rauto_gCm2day = output[,,2],
                    mean_rauto_gCm2day = output_mean[,2],
                    mean_annual_rauto_gCm2day = output_annual[,,2],
                    rhet_dom_gCm2day = output[,,3],
                    mean_rhet_dom_gCm2day = output_mean[,3],
                    mean_annual_rhet_dom_gCm2day = output_annual[,,3],
                    fire_gCm2day = output[,,4],
                    mean_fire_gCm2day = output_mean[,4],
                    mean_annual_fire_gCm2day = output_annual[,,4],
                    # Internal fluxes
                    alloc_foliage_gCm2day = output[,,5],
                    mean_alloc_foliage_gCm2day = output_mean[,5],
                    mean_annual_alloc_foliage_gCm2day = output_annual[,,5],
                    alloc_labile_gCm2day = output[,,6],
                    mean_alloc_labile_gCm2day = output_mean[,6],
                    mean_annual_alloc_labile_gCm2day = output_annual[,,6],
                    alloc_roots_wood_gCm2day = output[,,7],
                    mean_alloc_roots_wood_gCm2day = output_mean[,7],
                    mean_annual_alloc_roots_wood_gCm2day = output_annual[,,7],
                    labile_to_foliage_gCm2day = output[,,8],
                    mean_labile_to_foliage_gCm2day = output_mean[,8],
                    mean_annual_labile_to_foliage_gCm2day = output_annual[,,8],
                    foliage_to_litter_gCm2day = output[,,9],
                    mean_foliage_to_litter_gCm2day = output_mean[,9],
                    mean_annual_foliage_to_litter_gCm2day = output_annual[,,9],
                    roots_wood_to_litter_gCm2day = output[,,10],
                    mean_roots_wood_to_litter_gCm2day = output_mean[,10],
                    mean_annual_roots_wood_to_litter_gCm2day = output_annual[,,10],
                    # Disturbance fluxes
                    FIREemiss_labile_gCm2day = output[,,11],
                    mean_FIREemiss_labile_gCm2day = output_mean[,11],
                    mean_annual_FIREemiss_labile_gCm2day = output_annual[,,11],
                    FIRElitter_labile_gCm2day = output[,,12],
                    mean_FIRElitter_labile_gCm2day = output_mean[,12],
                    mean_annual_FIRElitter_labile_gCm2day = output_annual[,,12],
                    FIREemiss_foliage_gCm2day = output[,,13],
                    mean_FIREemiss_foliage_gCm2day = output_mean[,13],
                    mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,13],
                    FIRElitter_foliage_gCm2day = output[,,14],
                    mean_FIRElitter_foliage_gCm2day = output_mean[,14],
                    mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,14],
                    FIREemiss_roots_wood_gCm2day = output[,,15],
                    mean_FIREemiss_roots_wood_gCm2day = output_mean[,15],
                    mean_annual_FIREemiss_roots_wood_gCm2day = output_annual[,,15],
                    FIRElitter_roots_wood_gCm2day = output[,,16],
                    mean_FIRElitter_roots_wood_gCm2day = output_mean[,16],
                    mean_annual_FIRElitter_roots_wood_gCm2day = output_annual[,,16],
                    FIREemiss_dom_gCm2day = output[,,17],
                    mean_FIREemiss_dom_gCm2day = output_mean[,17],
                    mean_annual_FIREemiss_dom_gCm2day = output_annual[,,17],
                    # C pools (gC/m2)
                    labile_gCm2 = output[,,18],
                    mean_labile_gCm2 = output_mean[,18],
                    mean_annual_labile_gCm2 = output_annual[,,18],
                    foliage_gCm2 = output[,,19],
                    mean_foliage_gCm2 = output_mean[,19],
                    mean_annual_foliage_gCm2 = output_annual[,,19],
                    roots_wood_gCm2 = output[,,20],
                    mean_roots_wood_gCm2 = output_mean[,20],
                    mean_annual_roots_wood_gCm2 = output_annual[,,20],
                    dom_gCm2 = output[,,21],
                    mean_dom_gCm2 = output_mean[,21],
                    mean_annual_dom_gCm2 = output_annual[,,21],
                    # Canopy (phenology) properties
                    lai_m2m2 = output[,,22],
                    mean_lai_m2m2 = output_mean[,22],
                    mean_annual_lai_m2m2 = output_annual[,,22],
                    ## Aggregated variables
                    # Mean Transit times
                    MTT_labile_years = MTT_years[,1],
                    MTT_foliage_years = MTT_years[,2],
                    MTT_roots_wood_years = MTT_years[,3],
                    MTT_dom_years = MTT_years[,4],
                    # Steady state estimates
                    SS_labile_gCm2 = SS_gCm2[,1],
                    SS_foliage_gCm2 = SS_gCm2[,2],
                    SS_roots_wood_gCm2 = SS_gCm2[,3],
                    SS_dom_gCm2 = SS_gCm2[,4])
    # Determine the NPP fraction of expressed NPP
    # i.e. actual growth not GPP-Ra
    NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                         states_all$alloc_foliage_gCm2day +
                         states_all$alloc_roots_wood_gCm2day,1,mean)
    NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                         apply(states_all$alloc_roots_wood_gCm2day,1,mean)) / NPP_fraction
    states_all$NPP_foliage_fraction = NPP_fraction[,1]
    states_all$NPP_roots_wood_fraction = NPP_fraction[,2]
    # Tidy up variables
    rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.D1.F2.001") {
      output_dim = 36 ; MTT_dim = 5 ; SS_dim = 5
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec1",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                             ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                    
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodays=as.integer(dim(met)[1])
                             ,nos_years=as.integer(noyears)
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
      output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
      output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))            
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # Create the output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      mean_gpp_gCm2day = output_mean[,1],
                      mean_annual_gpp_gCm2day = output_annual[,,1],
                      rauto_gCm2day = output[,,2],
                      mean_rauto_gCm2day = output_mean[,2],
                      mean_annual_rauto_gCm2day = output_annual[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      mean_rhet_litter_gCm2day = output_mean[,3],
                      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                      rhet_som_gCm2day = output[,,4],
                      mean_rhet_som_gCm2day = output_mean[,4],
                      mean_annual_rhet_som_gCm2day = output_annual[,,4],
                      fire_gCm2day = output[,,5],
                      mean_fire_gCm2day = output_mean[,5],
                      mean_annual_fire_gCm2day = output_annual[,,5],
                      harvest_gCm2day = output[,,6],
                      mean_harvest_gCm2day = output_mean[,6],
                      mean_annual_harvest_gCm2day = output_annual[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      mean_alloc_foliage_gCm2day = output_mean[,7],
                      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
                      alloc_roots_gCm2day = output[,,8],
                      mean_alloc_roots_gCm2day = output_mean[,8],
                      mean_annual_alloc_roots_gCm2day = output_annual[,,8],
                      alloc_wood_gCm2day = output[,,9],
                      mean_alloc_wood_gCm2day = output_mean[,9],
                      mean_annual_alloc_wood_gCm2day = output_annual[,,9],
                      foliage_to_litter_gCm2day = output[,,10],
                      mean_foliage_to_litter_gCm2day = output_mean[,10],
                      mean_annual_foliage_to_litter_gCm2day = output_annual[,,10],
                      roots_to_litter_gCm2day = output[,,11],
                      mean_roots_to_litter_gCm2day = output_mean[,11],
                      mean_annual_roots_to_litter_gCm2day = output_annual[,,11],
                      wood_to_litter_gCm2day = output[,,12],
                      mean_wood_to_litter_gCm2day = output_mean[,12],
                      mean_annual_wood_to_litter_gCm2day = output_annual[,,12],
                      litter_to_som_gCm2day = output[,,13],
                      mean_litter_to_som_gCm2day = output_mean[,13],
                      mean_annual_litter_to_som_gCm2day = output_annual[,,13],
                      # Disturbance fluxes
                      FIREemiss_foliage_gCm2day = output[,,14],
                      mean_FIREemiss_foliage_gCm2day = output_mean[,14],
                      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,14],
                      FIRElitter_foliage_gCm2day = output[,,15],
                      mean_FIRElitter_foliage_gCm2day = output_mean[,15],
                      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,15],
                      FIREemiss_roots_gCm2day = output[,,16],
                      mean_FIREemiss_roots_gCm2day = output_mean[,16],
                      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,16],
                      FIRElitter_roots_gCm2day = output[,,17],
                      mean_FIRElitter_roots_gCm2day = output_mean[,17],
                      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,17],
                      FIREemiss_wood_gCm2day = output[,,18],
                      mean_FIREemiss_wood_gCm2day = output_mean[,18],
                      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,18],
                      FIRElitter_wood_gCm2day = output[,,19],
                      mean_FIRElitter_wood_gCm2day = output_mean[,19],
                      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,19],
                      FIREemiss_litter_gCm2day = output[,,20],
                      mean_FIREemiss_litter_gCm2day = output_mean[,20],
                      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,20],
                      FIRElitter_litter_gCm2day = output[,,21],
                      mean_FIRElitter_litter_gCm2day = output_mean[,21],
                      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,21],
                      FIREemiss_som_gCm2day = output[,,22],
                      mean_FIREemiss_som_gCm2day = output_mean[,22],
                      mean_annual_FIREemiss_som_gCm2day = output_annual[,,22],
                      HARVESTextracted_foliage_gCm2day = output[,,23],
                      mean_HARVESTextracted_foliage_gCm2day = output_mean[,23],
                      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,23],
                      HARVESTextracted_roots_gCm2day = output[,,24],
                      mean_HARVESTextracted_roots_gCm2day = output_mean[,24],
                      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,24],
                      HARVESTextracted_wood_gCm2day = output[,,25],
                      mean_HARVESTextracted_wood_gCm2day = output_mean[,25],
                      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,25],
                      HARVESTextracted_litter_gCm2day = output[,,26],
                      mean_HARVESTextracted_litter_gCm2day = output_mean[,26],
                      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,26],
                      HARVESTextracted_som_gCm2day = output[,,27],
                      mean_HARVESTextracted_som_gCm2day = output_mean[,27],
                      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,27],
                      HARVESTlitter_foliage_gCm2day = output[,,28],
                      mean_HARVESTlitter_foliage_gCm2day = output_mean[,28],
                      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,28],
                      HARVESTlitter_roots_gCm2day = output[,,29],
                      mean_HARVESTlitter_roots_gCm2day = output_mean[,29],
                      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,29],
                      HARVESTlitter_wood_gCm2day = output[,,30],
                      mean_HARVESTlitter_wood_gCm2day = output_mean[,30],
                      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,30],
                      # C pools (gC/m2)
                      foliage_gCm2 = output[,,31],
                      mean_foliage_gCm2 = output_mean[,31],
                      mean_annual_foliage_gCm2 = output_annual[,,31],
                      roots_gCm2 = output[,,32],
                      mean_roots_gCm2 = output_mean[,32],
                      mean_annual_roots_gCm2 = output_annual[,,32],
                      wood_gCm2 = output[,,33],
                      mean_wood_gCm2 = output_mean[,33],
                      mean_annual_wood_gCm2 = output_annual[,,33],
                      litter_gCm2 = output[,,34],
                      mean_litter_gCm2 = output_mean[,34],
                      mean_annual_litter_gCm2 = output_annual[,,34],
                      som_gCm2 = output[,,35],
                      mean_som_gCm2 = output_mean[,35],
                      mean_annual_som_gCm2 = output_annual[,,35],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,36],
                      mean_lai_m2m2 = output_mean[,36],
                      mean_annual_lai_m2m2 = output_annual[,,36],
                      ## Aggregated variables
                      # Mean Transit times
                      MTT_foliage_years = MTT_years[,1],
                      MTT_roots_years = MTT_years[,2],
                      MTT_wood_years = MTT_years[,3],
                      MTT_litter_years = MTT_years[,4],
                      MTT_som_years = MTT_years[,5],
                      # Steady state estimates
                      SS_foliage_gCm2 = SS_gCm2[,1],
                      SS_roots_gCm2 = SS_gCm2[,2],
                      SS_wood_gCm2 = SS_gCm2[,3],
                      SS_litter_gCm2 = SS_gCm2[,4],
                      SS_som_gCm2 = SS_gCm2[,5])
      # Determine the NPP fraction of expressed NPP
      # i.e. actual growth not GPP-Ra
      NPP_fraction = apply(states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
      NPP_fraction = cbind(apply(states_all$alloc_foliage_gCm2day,1,mean),
                           apply(states_all$alloc_roots_gCm2day,1,mean),
                           apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
      states_all$NPP_foliage_fraction = NPP_fraction[,1]
      states_all$NPP_roots_fraction = NPP_fraction[,2]
      states_all$NPP_wood_fraction = NPP_fraction[,3]
      # Final tidy
      rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.C4.D1.F2.012") {
    output_dim = 23 ; MTT_dim = 3 ; SS_dim = 3
    dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
    tmp=.Fortran( "rdalec12",output_dim=as.integer(output_dim)
                            ,MTT_dim=as.integer(MTT_dim)
                            ,SS_dim = as.integer(SS_dim)
                            ,met=as.double(t(met))
                            ,pars=as.double(pars_in)
                            ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                            ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                            ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                            ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                            ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                   
                            ,lat=as.double(lat)
                            ,nopars=as.integer(PROJECT$model$nopars[site])
                            ,nomet=as.integer(dim(met)[2])
                            ,nofluxes=as.integer(PROJECT$model$nofluxes[site])
                            ,nopools=as.integer(PROJECT$model$nopools[site])
                            ,nodays=as.integer(dim(met)[1])
                            ,nos_years=as.integer(noyears)
                            ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1]))))
                            ,nos_iter=as.integer(nos_iter) )
    output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
    output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
    output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))          
    MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
    SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
    # Unload the current dalec shared object
    dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
    rm(tmp) ; gc()
    # create output object
    states_all=list(# Ecosystem fluxes
                    gpp_gCm2day = output[,,1],
                    mean_gpp_gCm2day = output_mean[,1],
                    mean_annual_gpp_gCm2day = output_annual[,,1],
                    rauto_gCm2day = output[,,2],
                    mean_rauto_gCm2day = output_mean[,2],
                    mean_annual_rauto_gCm2day = output_annual[,,2],
                    rhet_dom_gCm2day = output[,,3],
                    mean_rhet_dom_gCm2day = output_mean[,3],
                    mean_annual_rhet_dom_gCm2day = output_annual[,,3],
                    fire_gCm2day = output[,,4],
                    mean_fire_gCm2day = output_mean[,4],
                    mean_annual_fire_gCm2day = output_annual[,,4],
                    harvest_gCm2day = output[,,5],
                    mean_harvest_gCm2day = output_mean[,5],
                    mean_annual_harvest_gCm2day = output_annual[,,5],
                    # Internal fluxes
                    alloc_foliage_gCm2day = output[,,6],
                    mean_alloc_foliage_gCm2day = output_mean[,6],
                    mean_annual_alloc_foliage_gCm2day = output_annual[,,6],
                    alloc_roots_wood_gCm2day = output[,,7],
                    mean_alloc_roots_wood_gCm2day = output_mean[,7],
                    mean_annual_alloc_roots_wood_gCm2day = output_annual[,,7],
                    foliage_to_litter_gCm2day = output[,,8],
                    mean_foliage_to_litter_gCm2day = output_mean[,8],
                    mean_annual_foliage_to_litter_gCm2day = output_annual[,,8],
                    roots_wood_to_litter_gCm2day = output[,,9],
                    mean_roots_wood_to_litter_gCm2day = output_mean[,9],
                    mean_annual_roots_wood_to_litter_gCm2day = output_annual[,,9],
                    # Disturbance fluxes
                    FIREemiss_foliage_gCm2day = output[,,10],
                    mean_FIREemiss_foliage_gCm2day = output_mean[,10],
                    mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,10],
                    FIRElitter_foliage_gCm2day = output[,,11],
                    mean_FIRElitter_foliage_gCm2day = output_mean[,11],
                    mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,11],
                    FIREemiss_roots_wood_gCm2day = output[,,12],
                    mean_FIREemiss_roots_wood_gCm2day = output_mean[,12],
                    mean_annual_FIREemiss_roots_wood_gCm2day = output_annual[,,12],
                    FIRElitter_roots_wood_gCm2day = output[,,13],
                    mean_FIRElitter_roots_wood_gCm2day = output_mean[,13],
                    mean_annual_FIRElitter_roots_wood_gCm2day = output_annual[,,13],
                    FIREemiss_dom_gCm2day = output[,,14],
                    mean_FIREemiss_dom_gCm2day = output_mean[,14],
                    mean_annual_FIREemiss_dom_gCm2day = output_annual[,,14],
                    HARVESTextracted_foliage_gCm2day = output[,,15],
                    mean_HARVESTextracted_foliage_gCm2day = output_mean[,15],
                    mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,15],
                    HARVESTextracted_roots_wood_gCm2day = output[,,16],
                    mean_HARVESTextracted_roots_wood_gCm2day = output_mean[,16],
                    mean_annual_HARVESTextracted_roots_wood_gCm2day = output_annual[,,16],
                    HARVESTextracted_som_gCm2day = output[,,17],
                    mean_HARVESTextracted_som_gCm2day = output_mean[,17],
                    mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,17],
                    HARVESTlitter_foliage_gCm2day = output[,,18],
                    mean_HARVESTlitter_foliage_gCm2day = output_mean[,18],
                    mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,18],
                    HARVESTlitter_roots_wood_gCm2day = output[,,19],
                    mean_HARVESTlitter_roots_wood_gCm2day = output_mean[,19],
                    mean_annual_HARVESTlitter_roots_wood_gCm2day = output_annual[,,19],
                    # C pools (gC/m2)
                    foliage_gCm2 = output[,,20],
                    mean_foliage_gCm2 = output[,,20],
                    mean_annual_foliage_gCm2 = output_annual[,,20],
                    roots_wood_gCm2 = output[,,21],
                    mean_roots_wood_gCm2 = output[,,21],
                    mean_annual_roots_wood_gCm2 = output_annual[,,21],
                    dom_gCm2 = output[,,22],
                    mean_dom_gCm2 = output[,,22],
                    mean_annual_dom_gCm2 = output_annual[,,22],
                    # Canopy (phenology) properties
                    lai_m2m2 = output[,,23],
                    mean_lai_m2m2 = output[,,23],
                    mean_annual_lai_m2m2 = output_annual[,,23],
                    ## Aggregated variables
                    # Mean Transit times
                    MTT_foliage_years = MTT_years[,1],
                    MTT_roots_wood_years = MTT_years[,2],
                    MTT_dom_years = MTT_years[,3],
                    # Steady state estimates
                    SS_foliage_gCm2 = SS_gCm2[,1],
                    SS_roots_wood_gCm2 = SS_gCm2[,2],
                    SS_dom_gCm2 = SS_gCm2[,3])
    # Determine the NPP fraction of expressed NPP
    # i.e. actual growth not GPP-Ra
    NPP_fraction = apply(states_all$alloc_foliage_gCm2day +
                         states_all$alloc_roots_wood_gCm2day,1,mean)
    NPP_fraction = cbind(apply(states_all$alloc_foliage_gCm2day,1,mean),
                         apply(states_all$alloc_roots_wood_gCm2day,1,mean)) / NPP_fraction
    states_all$NPP_foliage_fraction = NPP_fraction[,1]
    states_all$NPP_roots_wood_fraction = NPP_fraction[,2]
    # Tidy up variables
    rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C2.D2.F2.H1.P3.R1.008") {
    output_dim = 59 ; MTT_dim = 7 ; SS_dim = 7
    # Load the required dalec shared object
    dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
    tmp=.Fortran( "rdalec8",output_dim=as.integer(output_dim)
                           ,aNPP_dim=as.integer(aNPP_dim)
                           ,MTT_dim=as.integer(MTT_dim)
                           ,SS_dim = as.integer(SS_dim)
                           ,fire_dim=as.integer(fire_dim)
                           ,met=as.double(t(met))
                           ,pars=as.double(pars_in)
                           ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                           ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                           ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                           ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                           ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                  
                           ,lat=as.double(lat)
                           ,nopars=as.integer(PROJECT$model$nopars[site])
                           ,nomet=as.integer(dim(met)[2])
                           ,nofluxes=as.integer(PROJECT$model$nofluxes[site])
                           ,nopools=as.integer(PROJECT$model$nopools[site])
                           ,pft=as.integer(pft)
                           ,nodays=as.integer(dim(met)[1])
                           ,nos_years=as.integer(noyears)
                           ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1]))))
                           ,nos_iter=as.integer(nos_iter))
    # Extract the different output variables
    output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
    output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
    output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))          
    MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
    SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
    # Unload the current dalec shared object
    dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
    rm(tmp) ; gc()
    # create output object
    states_all=list(# Ecosystem fluxes
                    gpp_gCm2day = output[,,1],
                    mean_gpp_gCm2day = output_mean[,1],
                    mean_annual_gpp_gCm2day = output_annual[,,1],
                    rauto_gCm2day = output[,,2],
                    mean_rauto_gCm2day = output_mean[,2],
                    mean_annual_rauto_gCm2day = output_annual[,,2],
                    rgrow_foliage_gCm2day = output[,,3],
                    mean_rgrow_foliage_gCm2day = output_mean[,3],
                    mean_annual_rgrow_foliage_gCm2day = output_annual[,,3],
                    rhet_litter_gCm2day = output[,,4],
                    mean_rhet_litter_gCm2day = output_mean[,4],
                    mean_annual_rhet_litter_gCm2day = output_annual[,,4],
                    rhet_som_gCm2day = output[,,5],
                    mean_rhet_som_gCm2day = output_mean[,5],
                    mean_annual_rhet_som_gCm2day = output_annual[,,5],
                    rhet_woodlitter_gCm2day = output[,,6],
                    mean_rhet_woodlitter_gCm2day = output_mean[,6],
                    mean_annual_rhet_woodlitter_gCm2day = output_annual[,,6],
                    fire_gCm2day = output[,,7],
                    mean_fire_gCm2day = output_mean[,7],
                    mean_annual_fire_gCm2day = output_annual[,,7],
                    harvest_gCm2day = output[,,8],
                    mean_harvest_gCm2day = output_mean[,8],
                    mean_annual_harvest_gCm2day = output_annual[,,8],
                    # Internal fluxes
                    alloc_labile_gCm2day = output[,,9],
                    mean_alloc_labile_gCm2day = output_mean[,9],
                    mean_annual_alloc_labile_gCm2day = output_annual[,,9],
                    alloc_roots_gCm2day = output[,,10],
                    mean_alloc_roots_gCm2day = output_mean[,10],
                    mean_annual_alloc_roots_gCm2day = output_annual[,,10],
                    alloc_wood_gCm2day = output[,,11],
                    mean_alloc_wood_gCm2day = output_mean[,11],
                    mean_annual_alloc_wood_gCm2day = output_annual[,,11],
                    labile_to_foliage_gCm2day = output[,,12],
                    mean_labile_to_foliage_gCm2day = output_mean[,12],
                    mean_annual_labile_to_foliage_gCm2day = output_annual[,,12],
                    foliage_to_litter_gCm2day = output[,,13],
                    mean_foliage_to_litter_gCm2day = output_mean[,13],
                    mean_annual_foliage_to_litter_gCm2day = output_annual[,,13],
                    roots_to_litter_gCm2day = output[,,14],
                    mean_roots_to_litter_gCm2day = output_mean[,14],
                    mean_annual_roots_to_litter_gCm2day = output_annual[,,14],
                    wood_to_litter_gCm2day = output[,,15],
                    mean_wood_to_litter_gCm2day = output_mean[,15],
                    mean_annual_wood_to_litter_gCm2day = output_annual[,,15],
                    litter_to_som_gCm2day = output[,,16],
                    mean_litter_to_som_gCm2day = output_mean[,16],
                    mean_annual_litter_to_som_gCm2day = output_annual[,,16],
                    woodlitter_to_som_gCm2day = output[,,17],
                    mean_woodlitter_to_som_gCm2day = output_mean[,17],
                    mean_annual_woodlitter_to_som_gCm2day = output_annual[,,17],
                    # Disturbance fluxes
                    FIREemiss_labile_gCm2day = output[,,18],
                    mean_FIREemiss_labile_gCm2day = output_mean[,18],
                    mean_annual_FIREemiss_labile_gCm2day = output_annual[,,18],
                    FIRElitter_labile_gCm2day = output[,,19],
                    mean_FIRElitter_labile_gCm2day = output_mean[,19],
                    mean_annual_FIRElitter_labile_gCm2day = output_annual[,,19],
                    FIREemiss_foliage_gCm2day = output[,,20],
                    mean_FIREemiss_foliage_gCm2day = output_mean[,20],
                    mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,20],
                    FIRElitter_foliage_gCm2day = output[,,21],
                    mean_FIRElitter_foliage_gCm2day = output_mean[,21],
                    mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,21],
                    FIREemiss_roots_gCm2day = output[,,22],
                    mean_FIREemiss_roots_gCm2day = output[,,22],
                    mean_annual_FIREemiss_roots_gCm2day = output_annual[,,22],
                    FIRElitter_roots_gCm2day = output[,,23],
                    mean_FIRElitter_roots_gCm2day = output_mean[,23],
                    mean_annual_FIRElitter_roots_gCm2day = output_annual[,,23],
                    FIREemiss_wood_gCm2day = output[,,24],
                    mean_FIREemiss_wood_gCm2day = output_mean[,24],
                    mean_annual_FIREemiss_wood_gCm2day = output_annual[,,24],
                    FIRElitter_wood_gCm2day = output[,,25],
                    mean_FIRElitter_wood_gCm2day = output_mean[,25],
                    mean_annual_FIRElitter_wood_gCm2day = output_annual[,,25],
                    FIREemiss_litter_gCm2day = output[,,26],
                    mean_FIREemiss_litter_gCm2day = output_mean[,26],
                    mean_annual_FIREemiss_litter_gCm2day = output_annual[,,26],
                    FIRElitter_litter_gCm2day = output[,,27],
                    mean_FIRElitter_litter_gCm2day = output[,,27],
                    mean_annual_FIRElitter_litter_gCm2day = output_annual[,,27],
                    FIREemiss_woodlitter_gCm2day = output[,,28],
                    mean_FIREemiss_woodlitter_gCm2day = output_mean[,28],
                    mean_annual_FIREemiss_woodlitter_gCm2day = output_annual[,,28],
                    FIRElitter_woodlitter_gCm2day = output[,,29],
                    mean_FIRElitter_woodlitter_gCm2day = output_mean[,29],
                    mean_annual_FIRElitter_woodlitter_gCm2day = output_annual[,,29],
                    FIREemiss_som_gCm2day = output[,,30],
                    mean_FIREemiss_som_gCm2day = output_mean[,30],
                    mean_annual_FIREemiss_som_gCm2day = output_annual[,,30],
                    HARVESTextracted_labile_gCm2day = output[,,31],
                    mean_HARVESTextracted_labile_gCm2day = output_mean[,31],
                    mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,31],
                    HARVESTextracted_foliage_gCm2day = output[,,32],
                    mean_HARVESTextracted_foliage_gCm2day = output_mean[,32],
                    mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,32],
                    HARVESTextracted_roots_gCm2day = output[,,33],
                    mean_HARVESTextracted_roots_gCm2day = output_mean[,33],
                    mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,33],
                    HARVESTextracted_wood_gCm2day = output[,,34],
                    mean_HARVESTextracted_wood_gCm2day = output_mean[,34],
                    mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,34],
                    HARVESTextracted_litter_gCm2day = output[,,35],
                    mean_HARVESTextracted_litter_gCm2day = output_mean[,35],
                    mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,35],
                    HARVESTextracted_woodlitter_gCm2day = output[,,36],
                    mean_HARVESTextracted_woodlitter_gCm2day = output_mean[,36],
                    mean_annual_HARVESTextracted_woodlitter_gCm2day = output_annual[,,36],
                    HARVESTextracted_som_gCm2day = output[,,37],
                    mean_HARVESTextracted_som_gCm2day = output_mean[,37],
                    mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,37],
                    HARVESTlitter_labile_gCm2day = output[,,38],
                    mean_HARVESTlitter_labile_gCm2day = output_mean[,38],
                    mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,38],
                    HARVESTlitter_foliage_gCm2day = output[,,39],
                    mean_HARVESTlitter_foliage_gCm2day = output_mean[,39],
                    mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,39],
                    HARVESTlitter_roots_gCm2day = output[,,40],
                    mean_HARVESTlitter_roots_gCm2day = output_mean[,40],
                    mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,40],
                    HARVESTlitter_wood_gCm2day = output[,,41],
                    mean_HARVESTlitter_wood_gCm2day = output_mean[,41],
                    mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,41],
                    # C pools (gC/m2)
                    labile_gCm2 = output[,,42],
                    mean_labile_gCm2 = output[,,42],
                    mean_annual_labile_gCm2 = output_annual[,,42],
                    foliage_gCm2 = output[,,43],
                    mean_foliage_gCm2 = output[,,43],
                    mean_annual_foliage_gCm2 = output_annual[,,43],
                    roots_gCm2 = output[,,44],
                    mean_roots_gCm2 = output[,,44],
                    mean_annual_roots_gCm2 = output_annual[,,44],
                    wood_gCm2 = output[,,45],
                    mean_wood_gCm2 = output[,,45],
                    mean_annual_wood_gCm2 = output_annual[,,45],
                    litter_gCm2 = output[,,46],
                    mean_litter_gCm2 = output[,,46],
                    mean_annual_litter_gCm2 = output_annual[,,46],
                    woodlitter_gCm2 = output[,,47],
                    mean_woodlitter_gCm2 = output[,,47],
                    mean_annual_woodlitter_gCm2 = output_annual[,,47],
                    som_gCm2 = output[,,48],
                    mean_som_gCm2 = output[,,48],
                    mean_annual_som_gCm2 = output_annual[,,48],
                    # Canopy (phenology) properties
                    lai_m2m2 = output[,,49],
                    mean_lai_m2m2 = output[,,49],
                    mean_annual_lai_m2m2 = output_annual[,,49],
                    gsi = output[,,50],
                    mean_gsi = output[,,50],
                    mean_annual_gsi = output_annual[,,50],
                    gsi_itemp = output[,,51],
                    mean_gsi_itemp = output[,,51],
                    mean_annual_gsi_itemp = output_annual[,,51],
                    gsi_iphoto = output[,,52],
                    mean_gsi_iphoto = output[,,52],
                    mean_annual_gsi_iphoto = output_annual[,,52],
                    gsi_ivpd = output[,,53],
                    mean_gsi_ivpd = output[,,53],
                    mean_annual_gsi_ivpd = output_annual[,,53],
                    # Photosynthesis / C~water coupling related
                    gs_demand_supply_ratio = output[,,54],
                    mean_gs_demand_supply_ratio = output[,,54],
                    mean_annual_gs_demand_supply_ratio = output_annual[,,54],
                    gs_mmolH2Om2s = output[,,55],
                    mean_gs_mmolH2Om2s = output[,,55],
                    mean_annual_gs_mmolH2Om2s = output_annual[,,55],
                    APAR_MJm2day = output[,,56],
                    mean_APAR_MJm2day = output[,,56],
                    mean_annual_APAR_MJm2day = output_annual[,,56],
                    gb_mmolH2Om2s = output[,,57],
                    mean_gb_mmolH2Om2s = output[,,57],
                    mean_annual_gb_mmolH2Om2s = output_annual[,,57],
                    CiCa = output[,,58],
                    mean_CiCa = output[,,58],
                    mean_annual_CiCa = output_annual[,,58],
                    # Misc
                    RootDepth_m = output[,,59],
                    mean_RootDepth_m = output[,,59],
                    mean_annual_RootDepth_m = output_annual[,,59],
                    ## Aggregated variables
                    # Mean Transit times
                    MTT_labile_years = MTT_years[,1],
                    MTT_foliage_years = MTT_years[,2],
                    MTT_roots_years = MTT_years[,3],
                    MTT_wood_years = MTT_years[,4],
                    MTT_litter_years = MTT_years[,5],
                    MTT_woodlitter_years = MTT_years[,6],
                    MTT_som_years = MTT_years[,7],
                    # Steady state estimates
                    SS_labile_gCm2 = SS_gCm2[,1],
                    SS_foliage_gCm2 = SS_gCm2[,2],
                    SS_roots_gCm2 = SS_gCm2[,3],
                    SS_wood_gCm2 = SS_gCm2[,4],
                    SS_litter_gCm2 = SS_gCm2[,5],
                    SS_woodlitter_gCm2 = SS_gCm2[,6],
                    SS_som_gCm2 = SS_gCm2[,7])
    # Determine the NPP fraction of expressed NPP
    # i.e. actual growth not GPP-Ra
    NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                         states_all$alloc_roots_gCm2day +
                         states_all$alloc_wood_gCm2day,1,mean)
    NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day,1,mean),
                         apply(states_all$alloc_roots_gCm2day,1,mean),
                         apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
    states_all$NPP_foliage_fraction = NPP_fraction[,1]
    states_all$NPP_roots_fraction = NPP_fraction[,2]
    states_all$NPP_wood_fraction = NPP_fraction[,3]
    # Tidy up variables
    rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.M2.016") { 
    output_dim = 49 ; MTT_dim = 5 ; SS_dim = 5
    # Load the required dalec shared object
    dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
    tmp=.Fortran( "rdalec16",output_dim=as.integer(output_dim)
                            ,MTT_dim=as.integer(MTT_dim)
                            ,SS_dim = as.integer(SS_dim)
                            ,met=as.double(t(met))
                            ,pars=as.double(pars_in)
                            ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                            ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                            ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                            ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                            ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                  
                            ,lat=as.double(lat)
                            ,nopars=as.integer(PROJECT$model$nopars[site])
                            ,nomet=as.integer(dim(met)[2])
                            ,nofluxes=as.integer(PROJECT$model$nofluxes[site])
                            ,nopools=as.integer(PROJECT$model$nopools[site])
                            ,nodays=as.integer(dim(met)[1])
                            ,nos_years=as.integer(noyears)
                            ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1]))))
                            ,nos_iter=as.integer(nos_iter))
    # Extract the different output variables
    output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
    output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
    output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))              
    MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
    SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
    # Unload the current dalec shared object
    dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
    rm(tmp) ; gc()
    # create output object
    states_all=list(# Ecosystem fluxes
                    gpp_gCm2day = output[,,1],
                    mean_gpp_gCm2day = output_mean[,1],
                    mean_annual_gpp_gCm2day = output_annual[,,1],
                    rauto_gCm2day = output[,,2],
                    mean_rauto_gCm2day = output_mean[,2],
                    mean_annual_rauto_gCm2day = output_annual[,,2],
                    rhet_litter_gCm2day = output[,,3],
                    mean_rhet_litter_gCm2day = output_mean[,3],
                    mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                    rhet_som_gCm2day = output[,,4],
                    mean_rhet_som_gCm2day = output_mean[,4],
                    mean_annual_rhet_som_gCm2day = output_annual[,,4],
                    fire_gCm2day = output[,,5],
                    mean_fire_gCm2day = output_mean[,5],
                    mean_annual_fire_gCm2day = output_annual[,,5],
                    harvest_gCm2day = output[,,6],
                    mean_harvest_gCm2day = output_mean[,6],
                    mean_annual_harvest_gCm2day = output_annual[,,6],
                    grazing_gCm2day = output[,,7],
                    mean_grazing_gCm2day = output_mean[,7],
                    mean_annual_grazing_gCm2day = output_annual[,,7],
                    # Internal fluxes
                    alloc_foliage_gCm2day = output[,,8],
                    mean_alloc_foliage_gCm2day = output_mean[,8],
                    mean_annual_alloc_foliage_gCm2day = output_annual[,,8],
                    alloc_labile_gCm2day = output[,,9],
                    mean_alloc_labile_gCm2day = output_mean[,9],
                    mean_annual_alloc_labile_gCm2day = output_annual[,,9],
                    alloc_roots_gCm2day = output[,,10],
                    mean_alloc_roots_gCm2day = output_mean[,10],
                    mean_annual_alloc_roots_gCm2day = output_annual[,,10],
                    labile_to_foliage_gCm2day = output[,,11],
                    mean_labile_to_foliage_gCm2day = output_mean[,11],
                    mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                    foliage_to_litter_gCm2day = output[,,12],
                    mean_foliage_to_litter_gCm2day = output_mean[,12],
                    mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                    roots_to_litter_gCm2day = output[,,13],
                    mean_roots_to_litter_gCm2day = output_mean[,13],
                    mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                    litter_to_som_gCm2day = output[,,14],
                    mean_litter_to_som_gCm2day = output_mean[,14],
                    mean_annual_litter_to_som_gCm2day = output_annual[,,14],
                    # Disturbance fluxes
                    FIREemiss_labile_gCm2day = output[,,15],
                    mean_FIREemiss_labile_gCm2day = output_mean[,15],
                    mean_annual_FIREemiss_labile_gCm2day = output_annual[,,15],
                    FIRElitter_labile_gCm2day = output[,,16],
                    mean_FIRElitter_labile_gCm2day = output_mean[,16],
                    mean_annual_FIRElitter_labile_gCm2day = output_annual[,,16],
                    FIREemiss_foliage_gCm2day = output[,,17],
                    mean_FIREemiss_foliage_gCm2day = output_mean[,17],
                    mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,17],
                    FIRElitter_foliage_gCm2day = output[,,18],
                    mean_FIRElitter_foliage_gCm2day = output_mean[,18],
                    mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,18],
                    FIREemiss_roots_gCm2day = output[,,19],
                    mean_FIREemiss_roots_gCm2day = output_mean[,19],
                    mean_annual_FIREemiss_roots_gCm2day = output_annual[,,19],
                    FIRElitter_roots_gCm2day = output[,,20],
                    mean_FIRElitter_roots_gCm2day = output_mean[,20],
                    mean_annual_FIRElitter_roots_gCm2day = output_annual[,,20],
                    FIREemiss_litter_gCm2day = output[,,21],
                    mean_FIREemiss_litter_gCm2day = output_mean[,21],
                    mean_annual_FIREemiss_litter_gCm2day = output_annual[,,21],
                    FIRElitter_litter_gCm2day = output[,,22],
                    mean_FIRElitter_litter_gCm2day = output_mean[,22],
                    mean_annual_FIRElitter_litter_gCm2day = output_annual[,,22],
                    FIREemiss_som_gCm2day = output[,,23],
                    mean_FIREemiss_som_gCm2day = output_mean[,23],
                    mean_annual_FIREemiss_som_gCm2day = output_annual[,,23],
                    HARVESTextracted_labile_gCm2day = output[,,24],
                    mean_HARVESTextracted_labile_gCm2day = output_mean[,24],
                    mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,24],
                    HARVESTextracted_foliage_gCm2day = output[,,25],
                    mean_HARVESTextracted_foliage_gCm2day = output_mean[,25],
                    mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,25],
                    HARVESTextracted_roots_gCm2day = output[,,26],
                    mean_HARVESTextracted_roots_gCm2day = output_mean[,26],
                    mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,26],
                    HARVESTlitter_labile_gCm2day = output[,,27],
                    mean_HARVESTlitter_labile_gCm2day = output_mean[,27],
                    mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,27],
                    HARVESTlitter_foliage_gCm2day = output[,,28],
                    mean_HARVESTlitter_foliage_gCm2day = output_mean[,28],
                    mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,28],
                    HARVESTlitter_roots_gCm2day = output[,,29],
                    mean_HARVESTlitter_roots_gCm2day = output_mean[,29],
                    mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,29],
                    GRAZINGextracted_labile_gCm2day = output[,,30],
                    mean_GRAZINGextracted_labile_gCm2day = output_mean[,30],
                    mean_annual_GRAZINGextracted_labile_gCm2day = output_annual[,,30],
                    GRAZINGextracted_foliage_gCm2day = output[,,31],
                    mean_GRAZINGextracted_foliage_gCm2day = output_mean[,31],
                    mean_annual_GRAZINGextracted_foliage_gCm2day = output_annual[,,31],
                    GRAZINGextracted_roots_gCm2day = output[,,32],
                    mean_GRAZINGextracted_roots_gCm2day = output_mean[,32],
                    mean_annual_GRAZINGextracted_roots_gCm2day = output_annual[,,32],
                    GRAZINGlitter_labile_gCm2day = output[,,33],
                    mean_GRAZINGlitter_labile_gCm2day = output_mean[,33],
                    mean_annual_GRAZINGlitter_labile_gCm2day = output_annual[,,33],
                    GRAZINGlitter_foliage_gCm2day = output[,,34],
                    mean_GRAZINGlitter_foliage_gCm2day = output_mean[,34],
                    mean_annual_GRAZINGlitter_foliage_gCm2day = output_annual[,,34],
                    GRAZINGlitter_roots_gCm2day = output[,,35],
                    mean_GRAZINGlitter_roots_gCm2day = output_mean[,35],
                    mean_annual_GRAZINGlitter_roots_gCm2day = output_annual[,,35],
                    # Animal output (gC/m2/day)
                    animal_manure_to_soil_gCm2day = output[,,36],
                    mean_animal_manure_to_soil_gCm2day = output_mean[,36],
                    mean_annual_animal_manure_to_soil_gCm2day = output_annual[,,36],
                    animal_respiration_gCm2day = output[,,37],
                    mean_animal_respiration_gCm2day = output_mean[,37],
                    mean_annual_animal_respiration_gCm2day = output_annual[,,37],
                    animal_methane_gCm2day = output[,,38],
                    mean_animal_methane_gCm2day = output_mean[,38],
                    mean_annual_animal_methane_gCm2day = output_annual[,,38],
                    # C pools (gC/m2)
                    labile_gCm2 = output[,,39],
                    mean_labile_gCm2 = output_mean[,39],
                    mean_annual_labile_gCm2 = output_annual[,,39],
                    foliage_gCm2 = output[,,40],
                    mean_foliage_gCm2 = output_mean[,40],
                    mean_annual_foliage_gCm2 = output_annual[,,40],
                    roots_gCm2 = output[,,41],
                    mean_roots_gCm2 = output[,,41],
                    mean_annual_roots_gCm2 = output_annual[,,41],
                    litter_gCm2 = output[,,42],
                    mean_litter_gCm2 = output_mean[,42],
                    mean_annual_litter_gCm2 = output_annual[,,42],
                    som_gCm2 = output[,,43],
                    mean_som_gCm2 = output_mean[,43],
                    mean_annual_som_gCm2 = output_annual[,,43],
                    # Canopy (phenology) properties
                    lai_m2m2 = output[,,44],
                    mean_lai_m2m2 = output_mean[,44],
                    mean_annual_lai_m2m2 = output_annual[,,44],
                    # Photosynthesis diagnostic
                    CiCa = output[,,45], 
                    mean_CiCa = output_mean[,45], 
                    mean_annual_CiCa = output_annual[,,45], 
                    # Canopy Phenology 
                    gsi = output[,,46],
                    mean_gsi = output_mean[,46],
                    mean_annual_gsi = output_annual[,,46],
                    gsi_itemp = output[,,47],
                    mean_gsi_itemp = output_mean[,47],
                    mean_annual_gsi_itemp = output_annual[,,47],
                    gsi_iphoto = output[,,48],
                    mean_gsi_iphoto = output_mean[,48],
                    mean_annual_gsi_iphoto = output_annual[,,48],
                    gsi_ivpd = output[,,49],
                    mean_gsi_ivpd = output_mean[,49],
                    mean_annual_gsi_ivpd = output_annual[,,49],
                    ## Aggregated variables
                    # Mean Transit times
                    MTT_labile_years = MTT_years[,1],
                    MTT_foliage_years = MTT_years[,2],
                    MTT_roots_years = MTT_years[,3],
                    MTT_litter_years = MTT_years[,4],
                    MTT_som_years = MTT_years[,5],
                    # Steady state estimates
                    SS_labile_gCm2 = SS_gCm2[,1],
                    SS_foliage_gCm2 = SS_gCm2[,2],
                    SS_roots_gCm2 = SS_gCm2[,3],
                    SS_litter_gCm2 = SS_gCm2[,4],
                    SS_som_gCm2 = SS_gCm2[,5])
    # Determine the NPP fraction of expressed NPP
    # i.e. actual growth not GPP-Ra                    
    NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                         states_all$alloc_roots_gCm2day +
                         states_all$alloc_foliage_gCm2day,1,mean)
    NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                         apply(states_all$alloc_roots_gCm2day,1,mean)) / NPP_fraction
    states_all$NPP_foliage_fraction = NPP_fraction[,1]
    states_all$NPP_roots_fraction = NPP_fraction[,2]
    # Tidy up variables
    rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A3.H2.M2.017") { 
    output_dim = 67 ; MTT_dim = 5 ; SS_dim = 5
    # Load the required dalec shared object
    dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
    tmp=.Fortran( "rdalec17",output_dim=as.integer(output_dim)
                            ,MTT_dim=as.integer(MTT_dim)
                            ,SS_dim = as.integer(SS_dim)
                            ,met=as.double(t(met))
                            ,pars=as.double(pars_in)
                            ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                            ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                            ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                            ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                            ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                  
                            ,lat=as.double(lat)
                            ,nopars=as.integer(PROJECT$model$nopars[site])
                            ,nomet=as.integer(dim(met)[2])
                            ,nofluxes=as.integer(PROJECT$model$nofluxes[site])
                            ,nopools=as.integer(PROJECT$model$nopools[site])
                            ,nodays=as.integer(dim(met)[1])
                            ,nos_years=as.integer(noyears)
                            ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1]))))
                            ,nos_iter=as.integer(nos_iter))
    # Extract the different output variables
    output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
    output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
    output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))              
    MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
    SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
    # Unload the current dalec shared object
    dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
    rm(tmp) ; gc()
    # create output object
    states_all=list(# Ecosystem fluxes
                    gpp_gCm2day = output[,,1],
                    mean_gpp_gCm2day = output_mean[,1],
                    mean_annual_gpp_gCm2day = output_annual[,,1],
                    rauto_gCm2day = output[,,2],
                    mean_rauto_gCm2day = output_mean[,2],
                    mean_annual_rauto_gCm2day = output_annual[,,2],
                    rhet_litter_gCm2day = output[,,3],
                    mean_rhet_litter_gCm2day = output_mean[,3],
                    mean_annual_rhet_litter_gCm2day = output_annual[,,3],
                    rhet_som_gCm2day = output[,,4],
                    mean_rhet_som_gCm2day = output_mean[,4],
                    mean_annual_rhet_som_gCm2day = output_annual[,,4],
                    fire_gCm2day = output[,,5],
                    mean_fire_gCm2day = output_mean[,5],
                    mean_annual_fire_gCm2day = output_annual[,,5],
                    harvest_gCm2day = output[,,6],
                    mean_harvest_gCm2day = output_mean[,6],
                    mean_annual_harvest_gCm2day = output_annual[,,6],
                    grazing_gCm2day = output[,,7],
                    mean_grazing_gCm2day = output_mean[,7],
                    mean_annual_grazing_gCm2day = output_annual[,,7],
                    # Internal fluxes
                    alloc_foliage_gCm2day = output[,,8],
                    mean_alloc_foliage_gCm2day = output_mean[,8],
                    mean_annual_alloc_foliage_gCm2day = output_annual[,,8],
                    alloc_labile_gCm2day = output[,,9],
                    mean_alloc_labile_gCm2day = output_mean[,9],
                    mean_annual_alloc_labile_gCm2day = output_annual[,,9],
                    alloc_roots_gCm2day = output[,,10],
                    mean_alloc_roots_gCm2day = output_mean[,10],
                    mean_annual_alloc_roots_gCm2day = output_annual[,,10],
                    labile_to_foliage_gCm2day = output[,,11],
                    mean_labile_to_foliage_gCm2day = output_mean[,11],
                    mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
                    foliage_to_litter_gCm2day = output[,,12],
                    mean_foliage_to_litter_gCm2day = output_mean[,12],
                    mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
                    roots_to_litter_gCm2day = output[,,13],
                    mean_roots_to_litter_gCm2day = output_mean[,13],
                    mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
                    litter_to_som_gCm2day = output[,,14],
                    mean_litter_to_som_gCm2day = output_mean[,14],
                    mean_annual_litter_to_som_gCm2day = output_annual[,,14],
                    # Disturbance fluxes
                    FIREemiss_labile_gCm2day = output[,,15],
                    mean_FIREemiss_labile_gCm2day = output_mean[,15],
                    mean_annual_FIREemiss_labile_gCm2day = output_annual[,,15],
                    FIRElitter_labile_gCm2day = output[,,16],
                    mean_FIRElitter_labile_gCm2day = output_mean[,16],
                    mean_annual_FIRElitter_labile_gCm2day = output_annual[,,16],
                    FIREemiss_foliage_gCm2day = output[,,17],
                    mean_FIREemiss_foliage_gCm2day = output_mean[,17],
                    mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,17],
                    FIRElitter_foliage_gCm2day = output[,,18],
                    mean_FIRElitter_foliage_gCm2day = output_mean[,18],
                    mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,18],
                    FIREemiss_roots_gCm2day = output[,,19],
                    mean_FIREemiss_roots_gCm2day = output_mean[,19],
                    mean_annual_FIREemiss_roots_gCm2day = output_annual[,,19],
                    FIRElitter_roots_gCm2day = output[,,20],
                    mean_FIRElitter_roots_gCm2day = output_mean[,20],
                    mean_annual_FIRElitter_roots_gCm2day = output_annual[,,20],
                    FIREemiss_litter_gCm2day = output[,,21],
                    mean_FIREemiss_litter_gCm2day = output_mean[,21],
                    mean_annual_FIREemiss_litter_gCm2day = output_annual[,,21],
                    FIRElitter_litter_gCm2day = output[,,22],
                    mean_FIRElitter_litter_gCm2day = output_mean[,22],
                    mean_annual_FIRElitter_litter_gCm2day = output_annual[,,22],
                    FIREemiss_som_gCm2day = output[,,23],
                    mean_FIREemiss_som_gCm2day = output_mean[,23],
                    mean_annual_FIREemiss_som_gCm2day = output_annual[,,23],
                    HARVESTextracted_labile_gCm2day = output[,,24],
                    mean_HARVESTextracted_labile_gCm2day = output_mean[,24],
                    mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,24],
                    HARVESTextracted_foliage_gCm2day = output[,,25],
                    mean_HARVESTextracted_foliage_gCm2day = output_mean[,25],
                    mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,25],
                    HARVESTextracted_roots_gCm2day = output[,,26],
                    mean_HARVESTextracted_roots_gCm2day = output_mean[,26],
                    mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,26],
                    HARVESTlitter_labile_gCm2day = output[,,27],
                    mean_HARVESTlitter_labile_gCm2day = output_mean[,27],
                    mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,27],
                    HARVESTlitter_foliage_gCm2day = output[,,28],
                    mean_HARVESTlitter_foliage_gCm2day = output_mean[,28],
                    mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,28],
                    HARVESTlitter_roots_gCm2day = output[,,29],
                    mean_HARVESTlitter_roots_gCm2day = output_mean[,29],
                    mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,29],
                    GRAZINGextracted_labile_gCm2day = output[,,30],
                    mean_GRAZINGextracted_labile_gCm2day = output_mean[,30],
                    mean_annual_GRAZINGextracted_labile_gCm2day = output_annual[,,30],
                    GRAZINGextracted_foliage_gCm2day = output[,,31],
                    mean_GRAZINGextracted_foliage_gCm2day = output_mean[,31],
                    mean_annual_GRAZINGextracted_foliage_gCm2day = output_annual[,,31],
                    GRAZINGextracted_roots_gCm2day = output[,,32],
                    mean_GRAZINGextracted_roots_gCm2day = output_mean[,32],
                    mean_annual_GRAZINGextracted_roots_gCm2day = output_annual[,,32],
                    GRAZINGlitter_labile_gCm2day = output[,,33],
                    mean_GRAZINGlitter_labile_gCm2day = output_mean[,33],
                    mean_annual_GRAZINGlitter_labile_gCm2day = output_annual[,,33],
                    GRAZINGlitter_foliage_gCm2day = output[,,34],
                    mean_GRAZINGlitter_foliage_gCm2day = output_mean[,34],
                    mean_annual_GRAZINGlitter_foliage_gCm2day = output_annual[,,34],
                    GRAZINGlitter_roots_gCm2day = output[,,35],
                    mean_GRAZINGlitter_roots_gCm2day = output_mean[,35],
                    mean_annual_GRAZINGlitter_roots_gCm2day = output_annual[,,35],
                    # Animal output (gC/m2/day)
                    animal_manure_to_soil_gCm2day = output[,,36],
                    mean_animal_manure_to_soil_gCm2day = output_mean[,36],
                    mean_annual_animal_manure_to_soil_gCm2day = output_annual[,,36],
                    animal_respiration_gCm2day = output[,,37],
                    mean_animal_respiration_gCm2day = output_mean[,37],
                    mean_annual_animal_respiration_gCm2day = output_annual[,,37],
                    animal_methane_gCm2day = output[,,38],
                    mean_animal_methane_gCm2day = output_mean[,38],
                    mean_annual_animal_methane_gCm2day = output_annual[,,38],
                    # C pools (gC/m2)
                    labile_gCm2 = output[,,39],
                    mean_labile_gCm2 = output_mean[,39],
                    mean_annual_labile_gCm2 = output_annual[,,39],
                    foliage_gCm2 = output[,,40],
                    mean_foliage_gCm2 = output_mean[,40],
                    mean_annual_foliage_gCm2 = output_annual[,,40],
                    roots_gCm2 = output[,,41],
                    mean_roots_gCm2 = output[,,41],
                    mean_annual_roots_gCm2 = output_annual[,,41],
                    litter_gCm2 = output[,,42],
                    mean_litter_gCm2 = output_mean[,42],
                    mean_annual_litter_gCm2 = output_annual[,,42],
                    som_gCm2 = output[,,43],
                    mean_som_gCm2 = output_mean[,43],
                    mean_annual_som_gCm2 = output_annual[,,43],
                    # Canopy (phenology) properties
                    lai_m2m2 = output[,,44],
                    mean_lai_m2m2 = output_mean[,44],
                    mean_annual_lai_m2m2 = output_annual[,,44],
                    # Photosynthesis diagnostic
                    CiCa = output[,,45], 
                    mean_CiCa = output_mean[,45], 
                    mean_annual_CiCa = output_annual[,,45], 
                    gs_demand_supply_ratio = output[,,46],
                    mean_gs_demand_supply_ratio = output_mean[,46],
                    mean_annual_gs_demand_supply_ratio = output_annual[,,46],
                    gs_mmolH2Om2s = output[,,47],
                    mean_gs_mmolH2Om2s = output_mean[,47],
                    mean_annual_gs_mmolH2Om2s = output_annual[,,47],
                    APAR_MJm2day = output[,,48],
                    mean_APAR_MJm2day = output_mean[,48],
                    mean_annual_APAR_MJm2day = output_annual[,,48],
                    gb_mmolH2Om2s = output[,,49],
                    mean_gb_mmolH2Om2s = output_mean[,49],
                    mean_annual_gb_mmolH2Om2s = output_annual[,,49],
                    # Water cycle related
                    ET_kgH2Om2day = output[,,50],
                    mean_ET_kgH2Om2day = output_mean[,50],
                    mean_annual_ET_kgH2Om2day = output_annual[,,50],
                    Etrans_kgH2Om2day = output[,,51],
                    mean_Etrans_kgH2Om2day = output_mean[,51],
                    mean_annual_Etrans_kgH2Om2day = output_annual[,,51],
                    Esoil_kgH2Om2day = output[,,52],
                    mean_Esoil_kgH2Om2day = output_mean[,52],
                    mean_annual_Esoil_kgH2Om2day = output_annual[,,52],
                    Ewetcanopy_kgH2Om2day = output[,,53],
                    mean_Ewetcanopy_kgH2Om2day = output_mean[,53],
                    mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,53],
                    runoff_kgH2Om2day = output[,,54],
                    mean_runoff_kgH2Om2day = output_mean[,54],
                    mean_annual_runoff_kgH2Om2day = output_annual[,,54],
                    underflow_kgH2Om2day = output[,,55],
                    mean_underflow_kgH2Om2day = output_mean[,55],
                    mean_annual_underflow_kgH2Om2day = output_annual[,,55],
                    SurfDrainage_kgH2Om2day = output[,,56],
                    mean_SurfDrainage_kgH2Om2day = output_mean[,56],
                    mean_annual_SurfDrainage_kgH2Om2day = output_annual[,,56],
                    SurfInfiltrated_kgH2Om2day = output[,,57],
                    mean_SurfInfiltrated_kgH2Om2day = output_mean[,57],
                    mean_annual_SurfInfiltrated_kgH2Om2day = output_annual[,,57],
                    Etrans_1st_root_layer_uptake_fraction = output[,,58],
                    mean_Etrans_1st_root_layer_uptake_fraction = output_mean[,58],
                    mean_annual_Etrans_1st_root_layer_uptake_fraction = output_annual[,,58],
                    Etrans_2nd_root_layer_uptake_fraction = output[,,59],
                    mean_Etrans_2nd_root_layer_uptake_fraction = output_mean[,59],
                    mean_annual_Etrans_2nd_root_layer_uptake_fraction = output_annual[,,59],
                    SurfWater_kgH2Om2 = output[,,60],
                    mean_SurfWater_kgH2Om2 = output_mean[,60],
                    mean_annual_SurfWater_kgH2Om2 = output_annual[,,60],
                    wSWP_MPa = output[,,61],
                    mean_wSWP_MPa = output_mean[,61],
                    mean_annual_wSWP_MPa = output_annual[,,61],
                    snow_kgH2Om2 = output[,,62],
                    mean_snow_kgH2Om2 = output_mean[,62],
                    mean_annual_snow_kgH2Om2 = output_annual[,,62],
                    # Misc
                    RootDepth_m = output[,,63],
                    mean_RootDepth_m = output_mean[,63],
                    mean_annual_RootDepth_m = output_annual[,,63],
                    # Canopy Phenology 
                    gsi = output[,,64],
                    mean_gsi = output_mean[,64],
                    mean_annual_gsi = output_annual[,,64],
                    gsi_itemp = output[,,65],
                    mean_gsi_itemp = output_mean[,65],
                    mean_annual_gsi_itemp = output_annual[,,65],
                    gsi_iphoto = output[,,66],
                    mean_gsi_iphoto = output_mean[,66],
                    mean_annual_gsi_iphoto = output_annual[,,66],
                    gsi_ivpd = output[,,67],
                    mean_gsi_ivpd = output_mean[,67],
                    mean_annual_gsi_ivpd = output_annual[,,67],
                    ## Aggregated variables
                    # Mean Transit times
                    MTT_labile_years = MTT_years[,1],
                    MTT_foliage_years = MTT_years[,2],
                    MTT_roots_years = MTT_years[,3],
                    MTT_litter_years = MTT_years[,4],
                    MTT_som_years = MTT_years[,5],
                    # Steady state estimates
                    SS_labile_gCm2 = SS_gCm2[,1],
                    SS_foliage_gCm2 = SS_gCm2[,2],
                    SS_roots_gCm2 = SS_gCm2[,3],
                    SS_litter_gCm2 = SS_gCm2[,4],
                    SS_som_gCm2 = SS_gCm2[,5])
    # Determine the NPP fraction of expressed NPP
    # i.e. actual growth not GPP-Ra                    
    NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                         states_all$alloc_roots_gCm2day +
                         states_all$alloc_foliage_gCm2day,1,mean)
    NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                         apply(states_all$alloc_roots_gCm2day,1,mean)) / NPP_fraction
    states_all$NPP_foliage_fraction = NPP_fraction[,1]
    states_all$NPP_roots_fraction = NPP_fraction[,2]
    # Tidy up variables
    rm(output,MTT_years,SS_gCm2)
  } else if (model_name == "DALEC.A1.C1.D2.F2.H4.P1.032") {
    output_dim = 69 ; MTT_dim = 6 ; SS_dim = 6
    dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
    tmp=.Fortran( "rdalec32",output_dim=as.integer(output_dim)
                  ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                  ,met=as.double(t(met))
                  ,pars=as.double(pars_in)
                  ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                  ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                  ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                  ,out_var4=as.double(array(0,dim=c(nos_iter,output_dim)))
                  ,out_var5=as.double(array(0,dim=c(nos_iter,noyears,output_dim)))                                    
                  ,lat=as.double(lat)
                  ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                  ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                  ,nodiags=as.integer(PROJECT$model$nodiags[site]),nodays=as.integer(dim(met)[1])
                  ,nos_years=as.integer(noyears)
                  ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                  ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                  ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2]))
                  ,residual_waterfrac_in=as.double(c(soil_info[5],soil_info[6],soil_info[6])) 
                  ,porosity_in = as.double(c(soil_info[7],soil_info[8],soil_info[8]))
                  ,pore_size_dist_in = as.double(c(soil_info[9],soil_info[10],soil_info[10]))
                  ,air_entry_in = as.double(c(soil_info[11],soil_info[12],soil_info[12]))
                  ,sat_conductivity_in = as.double(c(soil_info[13],soil_info[14],soil_info[14])))
    output = tmp$out_var1        ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
    output_mean = tmp$out_var4   ; output_mean = array(output_mean, dim=c(nos_iter,output_dim))
    output_annual = tmp$out_var5 ; output_annual = array(output_annual, dim=c(nos_iter,noyears,output_dim))      
    MTT_years = tmp$out_var2     ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
    SS_gCm2 = tmp$out_var3       ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
    # Unload the current dalec shared object
    dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
    rm(tmp) ; gc()
    # create output object
    states_all=list(# Ecosystem fluxes
      gpp_gCm2day = output[,,1],
      mean_gpp_gCm2day = output_mean[,1],
      mean_annual_gpp_gCm2day = output_annual[,,1],
      rauto_gCm2day = output[,,2],
      mean_rauto_gCm2day = output_mean[,2],
      mean_annual_rauto_gCm2day = output_annual[,,2],
      rhet_litter_gCm2day = output[,,3],
      mean_rhet_litter_gCm2day = output_mean[,3],
      mean_annual_rhet_litter_gCm2day = output_annual[,,3],
      rhet_som_gCm2day = output[,,4],
      mean_rhet_som_gCm2day = output_mean[,4],
      mean_annual_rhet_som_gCm2day = output_annual[,,4],
      fire_gCm2day = output[,,5],
      mean_fire_gCm2day = output_mean[,5],
      mean_annual_fire_gCm2day = output_annual[,,5],
      harvest_gCm2day = output[,,6],
      mean_harvest_gCm2day = output_mean[,6],
      mean_annual_harvest_gCm2day = output_annual[,,6],
      # Internal fluxes
      alloc_foliage_gCm2day = output[,,7],
      mean_alloc_foliage_gCm2day = output_mean[,7],
      mean_annual_alloc_foliage_gCm2day = output_annual[,,7],
      alloc_labile_gCm2day = output[,,8],
      mean_alloc_labile_gCm2day = output_mean[,8],
      mean_annual_alloc_labile_gCm2day = output_annual[,,8],
      alloc_roots_gCm2day = output[,,9],
      mean_alloc_roots_gCm2day = output_mean[,9],
      mean_annual_alloc_roots_gCm2day = output_annual[,,9],
      alloc_wood_gCm2day = output[,,10],
      mean_alloc_wood_gCm2day = output_mean[,10],
      mean_annual_alloc_wood_gCm2day = output_annual[,,10],
      labile_to_foliage_gCm2day = output[,,11],
      mean_labile_to_foliage_gCm2day = output_mean[,11],
      mean_annual_labile_to_foliage_gCm2day = output_annual[,,11],
      foliage_to_litter_gCm2day = output[,,12],
      mean_foliage_to_litter_gCm2day = output_mean[,12],
      mean_annual_foliage_to_litter_gCm2day = output_annual[,,12],
      roots_to_litter_gCm2day = output[,,13],
      mean_roots_to_litter_gCm2day = output_mean[,13],
      mean_annual_roots_to_litter_gCm2day = output_annual[,,13],
      wood_to_litter_gCm2day = output[,,14],
      mean_wood_to_litter_gCm2day = output_mean[,14],
      mean_annual_wood_to_litter_gCm2day = output_annual[,,14],
      litter_to_som_gCm2day = output[,,15],
      mean_litter_to_som_gCm2day = output_mean[,15],
      mean_annual_litter_to_som_gCm2day = output_annual[,,15],
      # Disturbance fluxes
      FIREemiss_labile_gCm2day = output[,,16],
      mean_FIREemiss_labile_gCm2day = output_mean[,16],
      mean_annual_FIREemiss_labile_gCm2day = output_annual[,,16],
      FIRElitter_labile_gCm2day = output[,,17],
      mean_FIRElitter_labile_gCm2day = output_mean[,17],
      mean_annual_FIRElitter_labile_gCm2day = output_annual[,,17],
      FIREemiss_foliage_gCm2day = output[,,18],
      mean_FIREemiss_foliage_gCm2day = output_mean[,18],
      mean_annual_FIREemiss_foliage_gCm2day = output_annual[,,18],
      FIRElitter_foliage_gCm2day = output[,,19],
      mean_FIRElitter_foliage_gCm2day = output_mean[,19],
      mean_annual_FIRElitter_foliage_gCm2day = output_annual[,,19],
      FIREemiss_roots_gCm2day = output[,,20],
      mean_FIREemiss_roots_gCm2day = output_mean[,20],
      mean_annual_FIREemiss_roots_gCm2day = output_annual[,,20],
      FIRElitter_roots_gCm2day = output[,,21],
      mean_FIRElitter_roots_gCm2day = output_mean[,21],
      mean_annual_FIRElitter_roots_gCm2day = output_annual[,,21],
      FIREemiss_wood_gCm2day = output[,,22],
      mean_FIREemiss_wood_gCm2day = output_mean[,22],
      mean_annual_FIREemiss_wood_gCm2day = output_annual[,,22],
      FIRElitter_wood_gCm2day = output[,,23],
      mean_FIRElitter_wood_gCm2day = output_mean[,23],
      mean_annual_FIRElitter_wood_gCm2day = output_annual[,,23],
      FIREemiss_litter_gCm2day = output[,,24],
      mean_FIREemiss_litter_gCm2day = output_mean[,24],
      mean_annual_FIREemiss_litter_gCm2day = output_annual[,,24],
      FIRElitter_litter_gCm2day = output[,,25],
      mean_FIRElitter_litter_gCm2day = output_mean[,25],
      mean_annual_FIRElitter_litter_gCm2day = output_annual[,,25],
      FIREemiss_som_gCm2day = output[,,26],
      mean_FIREemiss_som_gCm2day = output_mean[,26],
      mean_annual_FIREemiss_som_gCm2day = output_annual[,,26],
      HARVESTextracted_labile_gCm2day = output[,,27],
      mean_HARVESTextracted_labile_gCm2day = output_mean[,27],
      mean_annual_HARVESTextracted_labile_gCm2day = output_annual[,,27],
      HARVESTextracted_foliage_gCm2day = output[,,28],
      mean_HARVESTextracted_foliage_gCm2day = output_mean[,28],
      mean_annual_HARVESTextracted_foliage_gCm2day = output_annual[,,28],
      HARVESTextracted_roots_gCm2day = output[,,29],
      mean_HARVESTextracted_roots_gCm2day = output_mean[,29],
      mean_annual_HARVESTextracted_roots_gCm2day = output_annual[,,29],
      HARVESTextracted_wood_gCm2day = output[,,30],
      mean_HARVESTextracted_wood_gCm2day = output_mean[,30],
      mean_annual_HARVESTextracted_wood_gCm2day = output_annual[,,30],
      HARVESTextracted_litter_gCm2day = output[,,31],
      mean_HARVESTextracted_litter_gCm2day = output_mean[,31],
      mean_annual_HARVESTextracted_litter_gCm2day = output_annual[,,31],
      HARVESTextracted_som_gCm2day = output[,,32],
      mean_HARVESTextracted_som_gCm2day = output_mean[,32],
      mean_annual_HARVESTextracted_som_gCm2day = output_annual[,,32],
      HARVESTlitter_labile_gCm2day = output[,,33],
      mean_HARVESTlitter_labile_gCm2day = output_mean[,33],
      mean_annual_HARVESTlitter_labile_gCm2day = output_annual[,,33],
      HARVESTlitter_foliage_gCm2day = output[,,34],
      mean_HARVESTlitter_foliage_gCm2day = output_mean[,34],
      mean_annual_HARVESTlitter_foliage_gCm2day = output_annual[,,34],
      HARVESTlitter_roots_gCm2day = output[,,35],
      mean_HARVESTlitter_roots_gCm2day = output_mean[,35],
      mean_annual_HARVESTlitter_roots_gCm2day = output_annual[,,35],
      HARVESTlitter_wood_gCm2day = output[,,36],
      mean_HARVESTlitter_wood_gCm2day = output_mean[,36],
      mean_annual_HARVESTlitter_wood_gCm2day = output_annual[,,36],
      # C pools (gC/m2)
      labile_gCm2 = output[,,37],
      mean_labile_gCm2 = output_mean[,37],
      mean_annual_labile_gCm2 = output_annual[,,37],
      foliage_gCm2 = output[,,38],
      mean_foliage_gCm2 = output_mean[,38],
      mean_annual_foliage_gCm2 = output_annual[,,38],
      roots_gCm2 = output[,,39],
      mean_roots_gCm2 = output_mean[,39],
      mean_annual_roots_gCm2 = output_annual[,,39],
      wood_gCm2 = output[,,40],
      mean_wood_gCm2 = output_mean[,40],
      mean_annual_wood_gCm2 = output_annual[,,40],
      litter_gCm2 = output[,,41],
      mean_litter_gCm2 = output_mean[,41],
      mean_annual_litter_gCm2 = output_annual[,,41],
      som_gCm2 = output[,,42],
      mean_som_gCm2 = output_mean[,42],
      mean_annual_som_gCm2 = output_annual[,,42],
      # Water cycle related
      ET_kgH2Om2day = output[,,43],
      mean_ET_kgH2Om2day = output_mean[,43],
      mean_annual_ET_kgH2Om2day = output_annual[,,43],
      Etrans_kgH2Om2day = output[,,44],
      mean_Etrans_kgH2Om2day = output_mean[,44],
      mean_annual_Etrans_kgH2Om2day = output_annual[,,44],
      Esoil_kgH2Om2day = output[,,45],
      mean_Esoil_kgH2Om2day = output_mean[,45],
      mean_annual_Esoil_kgH2Om2day = output_annual[,,45],
      Ewetcanopy_kgH2Om2day = output[,,46],
      mean_Ewetcanopy_kgH2Om2day = output_mean[,46],
      mean_annual_Ewetcanopy_kgH2Om2day = output_annual[,,46],
      runoff_kgH2Om2day = output[,,47],
      mean_runoff_kgH2Om2day = output_mean[,47],
      mean_annual_runoff_kgH2Om2day = output_annual[,,47],
      underflow_kgH2Om2day = output[,,48],
      mean_underflow_kgH2Om2day = output_mean[,48],
      mean_annual_underflow_kgH2Om2day = output_annual[,,48],
      SurfDrainage_kgH2Om2day = output[,,49],
      mean_SurfDrainage_kgH2Om2day = output_mean[,49],
      mean_annual_SurfDrainage_kgH2Om2day = output_annual[,,49],
      SurfInfiltrated_kgH2Om2day = output[,,50],
      mean_SurfInfiltrated_kgH2Om2day = output_mean[,50],
      mean_annual_SurfInfiltrated_kgH2Om2day = output_annual[,,50],
      Etrans_1st_root_layer_uptake_fraction = output[,,51],
      mean_Etrans_1st_root_layer_uptake_fraction = output_mean[,51],
      mean_annual_Etrans_1st_root_layer_uptake_fraction = output_annual[,,51],
      Etrans_2nd_root_layer_uptake_fraction = output[,,52],
      mean_Etrans_2nd_root_layer_uptake_fraction = output_mean[,52],
      mean_annual_Etrans_2nd_root_layer_uptake_fraction = output_annual[,,52],
      SurfWater_kgH2Om2 = output[,,53],
      mean_SurfWater_kgH2Om2 = output_mean[,53],
      mean_annual_SurfWater_kgH2Om2 = output_annual[,,53],
      wSWP_MPa = output[,,54],
      mean_wSWP_MPa = output_mean[,54],
      mean_annual_wSWP_MPa = output_annual[,,54],
      snow_kgH2Om2 = output[,,55],
      mean_snow_kgH2Om2 = output_mean[,55],
      mean_annual_snow_kgH2Om2 = output_annual[,,55],
      hydraulic_conductivity_ms = output[,,66],
      relative_waterfrac_m3m3 = output[,,67],
      soil_waterpot_MPa = output[,,68],
      field_capacity = output[,,69],
      # Canopy (phenology) properties
      lai_m2m2 = output[,,56],
      mean_lai_m2m2 = output_mean[,56],
      mean_annual_lai_m2m2 = output_annual[,,56],
      # Photosynthesis / C~water coupling related
      gs_demand_supply_ratio = output[,,57],
      mean_gs_demand_supply_ratio = output_mean[,57],
      mean_annual_gs_demand_supply_ratio = output_annual[,,57],
      gs_mmolH2Om2s = output[,,58],
      mean_gs_mmolH2Om2s = output_mean[,58],
      mean_annual_gs_mmolH2Om2s = output_annual[,,58],
      APAR_MJm2day = output[,,59],
      mean_APAR_MJm2day = output_mean[,59],
      mean_annual_APAR_MJm2day = output_annual[,,59],
      gb_mmolH2Om2s = output[,,60],
      mean_gb_mmolH2Om2s = output_mean[,60],
      mean_annual_gb_mmolH2Om2s = output_annual[,,60],
      CiCa = output[,,61],
      mean_CiCa = output_mean[,61],
      mean_annual_CiCa = output_annual[,,61],
      # Misc
      RootDepth_m = output[,,62],
      mean_RootDepth_m = output_mean[,62],
      mean_annual_RootDepth_m = output_annual[,,62],
      # Leaf water potential
      LWP_MPa = output[,,63],
      mean_LWP_MPa = output_mean[,63],
      mean_annual_LWP_MPa = output_annual[,,63], 
      # Canopy aerodynamic diagnostics
      canopy_area_scaling_wind = output[,,64],
      mean_canopy_area_scaling_wind = output_mean[,64],
      mean_annual_canopy_area_scaling_wind = output_annual[,,64],
      canopy_area_scaling_light = output[,,65],
      mean_canopy_area_scaling_light = output_mean[,65],
      mean_annual_canopy_area_scaling_light = output_annual[,,65],    
      ## Aggregated variables
      # Mean Transit times
      MTT_labile_years = MTT_years[,1],
      MTT_foliage_years = MTT_years[,2],
      MTT_roots_years = MTT_years[,3],
      MTT_wood_years = MTT_years[,4],
      MTT_litter_years = MTT_years[,5],
      MTT_som_years = MTT_years[,6],
      # Steady state estimates
      SS_labile_gCm2 = SS_gCm2[,1],
      SS_foliage_gCm2 = SS_gCm2[,2],
      SS_roots_gCm2 = SS_gCm2[,3],
      SS_wood_gCm2 = SS_gCm2[,4],
      SS_litter_gCm2 = SS_gCm2[,5],
      SS_som_gCm2 = SS_gCm2[,6])
    # Determine the NPP fraction of expressed NPP
    # i.e. actual growth not GPP-Ra
    NPP_fraction = apply(states_all$labile_to_foliage_gCm2day +
                           states_all$alloc_foliage_gCm2day +
                           states_all$alloc_roots_gCm2day +
                           states_all$alloc_wood_gCm2day,1,mean)
    NPP_fraction = cbind(apply(states_all$labile_to_foliage_gCm2day+states_all$alloc_foliage_gCm2day,1,mean),
                         apply(states_all$alloc_roots_gCm2day,1,mean),
                         apply(states_all$alloc_wood_gCm2day,1,mean)) / NPP_fraction
    states_all$NPP_foliage_fraction = NPP_fraction[,1]
    states_all$NPP_roots_fraction = NPP_fraction[,2]
    states_all$NPP_wood_fraction = NPP_fraction[,3]
    # Tidy up variables
    rm(output,MTT_years,SS_gCm2)
  } else {
    stop(paste("Model choice (",model_name,") does not have corresponding R interface",sep=""))
  }
  

  # return state variable means
  return(states_all) ; gc(verbose=FALSE)

} # end of function
## Use byte compile
simulate_all<-cmpfun(simulate_all)
