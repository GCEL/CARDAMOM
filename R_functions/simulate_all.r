
###
## Function to generate mean state variable information by running the parameters and model choice
###

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
  } else if (model_name == "DALEC.C3.M1.#") {
# THIS CODE AND THE R INTERFACE NEED UPDATING TO MAKE OPTIMAL USE OF THE VARIABLES HERE
      output_dim=28
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      crop_file_location=paste(PROJECT$exepath,"winter_wheat_development.csv", sep="")
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdalec14",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim)
                                      ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                                      ,fire_dim=as.integer(fire_dim)
                                      ,met=as.double(t(met)),pars=as.double(pars_in)
                                      ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                                      ,out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                      ,out_var3=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                      ,out_var4=as.double(array(0,dim=c(nos_iter,SS_dim)))
                                      ,out_var5=as.double(array(0,dim=c(nos_iter,MTT_dim,noyears)))
                                      ,out_var6=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                      ,out_var7=as.double(array(0,dim=c(nos_iter,fire_dim,noyears)))
                                      ,out_var8=as.double(array(0,dim=c(nos_iter,fire_dim,noyears)))
                                      ,out_var9=as.double(array(0,dim=c(nos_iter,fire_dim,noyears)))
                                      ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site])
                                      ,nomet=as.integer(dim(met)[2]),nofluxes=as.integer(PROJECT$model$nofluxes[site])
                                      ,nopools=as.integer(PROJECT$model$nopools[site]),pft=as.integer(pft)
                                      ,pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                      ,noyears=as.integer(noyears)
                                      ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                                      ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                                      ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3)))
                                      ,exepath=as.character(crop_file_location),pathlength=as.integer(nchar(crop_file_location)))
      output = tmp$out_var   ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP = tmp$out_var2    ; aNPP = array(aNPP, dim=c(nos_iter,aNPP_dim))
      MTT = tmp$out_var3     ; MTT = array(MTT, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var4 ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      aMTT = tmp$out_var5    ; aMTT = array(aMTT, dim=c(nos_iter,MTT_dim,noyears))
      MTTnat = tmp$out_var6  ; MTTnat = array(MTTnat, dim=c(nos_iter,MTT_dim))
      FIREemiss = tmp$out_var7; FIREemiss = array(FIREemiss, dim=c(nos_iter,fire_dim,noyears))
      FIRElit = tmp$out_var8  ; FIRElit = array(FIRElit, dim=c(nos_iter,fire_dim,noyears))
      outflux_nat = tmp$out_var9  ; outflux_nat = array(outflux_nat, dim=c(nos_iter,fire_dim,noyears))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object (14,15,16,17 unused)
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], litter_gCm2 = output[,,10],
                      labile_gCm2 = output[,,11], foliage_gCm2 = output[,,12],
                      harvest_gCm2day = output[,,13], ET_kgH2Om2day = output[,,18],
                      SurfWater_kgH2Om2 = output[,,19], wSWP_MPa = output[,,20],
                      woodlitter_gCm2 = output[,,21], auto_gCm2 = output[,,22],
                      fire_gCm2day = output[,,23],
                      gs_demand_supply = output[,,24], gs_total_canopy = output[,,25],
                      APAR_MJm2day = output[,,26], gb_total_canopy = output[,,27],
                      CiCa = output[,,28],
                      aNPP = aNPP, MTT = MTT, SS_gCm2 = SS_gCm2, aMTT = aMTT, natMTT = MTTnat,
                      FIREemiss_gCm2yr = FIREemiss, FIRElit_gCm2yr = FIRElit,
                      NAToutflux_gCm2yr = outflux_nat)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
      # Final tidy
      rm(output,MTT_gCm2,SS_gCm2)
  } else if (model_name == "DALEC.A1.C3.H2.M1.#") {
# THIS CODE AND THE R INTERFACE NEED UPDATING TO MAKE OPTIMAL USE OF THE VARIABLES HERE
      output_dim=28
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      crop_file_location=paste(PROJECT$exepath,"winter_wheat_development.csv", sep="")
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdalec15",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim)
                                      ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                                      ,fire_dim=as.integer(fire_dim)
                                      ,met=as.double(t(met)),pars=as.double(pars_in)
                                      ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                                      ,out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                      ,out_var3=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                      ,out_var4=as.double(array(0,dim=c(nos_iter,SS_dim)))
                                      ,out_var5=as.double(array(0,dim=c(nos_iter,MTT_dim,noyears)))
                                      ,out_var6=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                      ,out_var7=as.double(array(0,dim=c(nos_iter,fire_dim,noyears)))
                                      ,out_var8=as.double(array(0,dim=c(nos_iter,fire_dim,noyears)))
                                      ,out_var9=as.double(array(0,dim=c(nos_iter,fire_dim,noyears)))
                                      ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site])
                                      ,nomet=as.integer(dim(met)[2]),nofluxes=as.integer(PROJECT$model$nofluxes[site])
                                      ,nopools=as.integer(PROJECT$model$nopools[site]),pft=as.integer(pft)
                                      ,pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                      ,noyears=as.integer(noyears)
                                      ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                                      ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                                      ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3)))
                                      ,exepath=as.character(crop_file_location),pathlength=as.integer(nchar(crop_file_location)))
      output = tmp$out_var   ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP = tmp$out_var2    ; aNPP = array(aNPP, dim=c(nos_iter,aNPP_dim))
      MTT = tmp$out_var3     ; MTT = array(MTT, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var4 ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      aMTT = tmp$out_var5    ; aMTT = array(aMTT, dim=c(nos_iter,MTT_dim,noyears))
      MTTnat = tmp$out_var6  ; MTTnat = array(MTTnat, dim=c(nos_iter,MTT_dim))
      FIREemiss = tmp$out_var7; FIREemiss = array(FIREemiss, dim=c(nos_iter,fire_dim,noyears))
      FIRElit = tmp$out_var8  ; FIRElit = array(FIRElit, dim=c(nos_iter,fire_dim,noyears))
      outflux_nat = tmp$out_var9  ; outflux_nat = array(outflux_nat, dim=c(nos_iter,fire_dim,noyears))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object (14,15,16,17 unused)
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], litter_gCm2 = output[,,10],
                      labile_gCm2 = output[,,11], foliage_gCm2 = output[,,12],
                      harvest_gCm2day = output[,,13], ET_kgH2Om2day = output[,,18],
                      SurfWater_kgH2Om2 = output[,,19], wSWP_MPa = output[,,20],
                      woodlitter_gCm2 = output[,,21], auto_gCm2 = output[,,22],
                      fire_gCm2day = output[,,23],
                      gs_demand_supply = output[,,24], gs_total_canopy = output[,,25],
                      APAR_MJm2day = output[,,26], gb_total_canopy = output[,,27],
                      CiCa = output[,,28],
                      aNPP = aNPP, MTT = MTT, SS_gCm2 = SS_gCm2, aMTT = aMTT, natMTT = MTTnat,
                      FIREemiss_gCm2yr = FIREemiss, FIRElit_gCm2yr = FIRElit,
                      NAToutflux_gCm2yr = outflux_nat)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
      # Final tidy
      rm(output,MTT_gCm2,SS_gCm2)
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P4.R2.#") {
    output_dim = 62 ; MTT_dim = 7 ; SS_dim = 7
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
                            ,lat=as.double(lat)
                            ,nopars=as.integer(PROJECT$model$nopars[site])
                            ,nomet=as.integer(dim(met)[2])
                            ,nofluxes=as.integer(PROJECT$model$nofluxes[site])
                            ,nopools=as.integer(PROJECT$model$nopools[site])
                            ,pft=as.integer(pft)
                            ,nodays=as.integer(dim(met)[1])
                            ,noyears=as.integer(noyears)
                            ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1]))))
                            ,nos_iter=as.integer(nos_iter)
                            ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                            ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3))))
    # Extract the different output variables
    output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
    MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
    SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
    # Unload the current dalec shared object
    dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
    rm(tmp) ; gc()
    # create output object
    states_all=list(# Ecosystem fluxes
                    gpp_gCm2day = output[,,1],
                    rauto_gCm2day = output[,,2],
                    rg_foliage_gCm2day = output[,,3],
                    rhet_litter_gCm2day = output[,,4],
                    rhet_som_gCm2day = output[,,5],
                    rhet_woodlitter_gCm2day = output[,,6],
                    fire_gCm2day = output[,,7],
                    harvest_gCm2day = output[,,8],
                    # Internal fluxes
                    alloc_labile_gCm2day = output[,,9],
                    alloc_roots_gCm2day = output[,,10],
                    alloc_wood_gCm2day = output[,,11],
                    labile_to_foliage_gCm2day = output[,,12],
                    foliage_to_litter_gCm2day = output[,,13],
                    roots_to_litter_gCm2day = output[,,14],
                    wood_to_litter_gCm2day = output[,,15],
                    litter_to_som_gCm2day = output[,,16],
                    woodlitter_to_som_gCm2day = output[,,17],
                    # Disturbance fluxes
                    FIREemiss_labile_gCm2day = output[,,18],
                    FIRElitter_labile_gCm2day = output[,,19],
                    FIREemiss_foliage_gCm2day = output[,,20],
                    FIRElitter_foliage_gCm2day = output[,,21],
                    FIREemiss_roots_gCm2day = output[,,22],
                    FIRElitter_roots_gCm2day = output[,,23],
                    FIREemiss_wood_gCm2day = output[,,24],
                    FIRElitter_wood_gCm2day = output[,,25],
                    FIREemiss_litter_gCm2day = output[,,26],
                    FIRElitter_litter_gCm2day = output[,,27],
                    FIREemiss_woodlitter_gCm2day = output[,,28],
                    FIRElitter_woodlitter_gCm2day = output[,,29],
                    FIREemiss_som_gCm2day = output[,,30],
                    HARVESTextracted_labile_gCm2day = output[,,31],
                    HARVESTextracted_foliage_gCm2day = output[,,32],
                    HARVESTextracted_roots_gCm2day = output[,,33],
                    HARVESTextracted_wood_gCm2day = output[,,34],
                    HARVESTextracted_litter_gCm2day = output[,,35],
                    HARVESTextracted_woodlitter_gCm2day = output[,,36],
                    HARVESTextracted_som_gCm2day = output[,,37],
                    HARVESTlitter_labile_gCm2day = output[,,38],
                    HARVESTlitter_foliage_gCm2day = output[,,39],
                    HARVESTlitter_roots_gCm2day = output[,,40],
                    HARVESTlitter_wood_gCm2day = output[,,41],
                    # C pools (gC/m2)
                    labile_gCm2 = output[,,42],
                    foliage_gCm2 = output[,,43],
                    roots_gCm2 = output[,,44],
                    wood_gCm2 = output[,,45],
                    litter_gCm2 = output[,,46],
                    woodlitter_gCm2 = output[,,47],
                    som_gCm2 = output[,,48],
                    # Water cycle related
                    ET_kgH2Om2day = output[,,49],
                    SurfWater_kgH2Om2 = output[,,50],
                    wSWP_MPa = output[,,51],
                    # Canopy (phenology) properties
                    lai_m2m2 = output[,,52],
                    gsi = output[,,53],
                    gsi_itemp = output[,,54],
                    gsi_iphoto = output[,,55],
                    gsi_ivpd = output[,,56],
                    # Photosynthesis / C~water coupling related
                    gs_demand_supply_ratio = output[,,57],
                    gs_mmolH2Om2s = output[,,58],
                    APAR_MJm2day = output[,,59],
                    gb_mmolH2Om2s = output[,,60],
                    CiCa = output[,,61],
                    # Misc
                    RootDepth_m = output[,,62],
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
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P7.R2.#") {
      output_dim = 59 ; MTT_dim = 7 ; SS_dim = 7
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec23",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                              ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2])))
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      rauto_gCm2day = output[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      rhet_som_gCm2day = output[,,4],
                      rhet_woodlitter_gCm2day = output[,,5],
                      fire_gCm2day = output[,,6],
                      harvest_gCm2day = output[,,7],
                      # Internal fluxes
                      alloc_labile_gCm2day = output[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      woodlitter_to_som_gCm2day = output[,,16],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,17],
                      FIRElitter_labile_gCm2day = output[,,18],
                      FIREemiss_foliage_gCm2day = output[,,19],
                      FIRElitter_foliage_gCm2day = output[,,20],
                      FIREemiss_roots_gCm2day = output[,,21],
                      FIRElitter_roots_gCm2day = output[,,22],
                      FIREemiss_wood_gCm2day = output[,,23],
                      FIRElitter_wood_gCm2day = output[,,24],
                      FIREemiss_litter_gCm2day = output[,,25],
                      FIRElitter_litter_gCm2day = output[,,26],
                      FIREemiss_woodlitter_gCm2day = output[,,27],
                      FIRElitter_woodlitter_gCm2day = output[,,28],
                      FIREemiss_som_gCm2day = output[,,29],
                      HARVESTextracted_labile_gCm2day = output[,,30],
                      HARVESTextracted_foliage_gCm2day = output[,,31],
                      HARVESTextracted_roots_gCm2day = output[,,32],
                      HARVESTextracted_wood_gCm2day = output[,,33],
                      HARVESTextracted_litter_gCm2day = output[,,34],
                      HARVESTextracted_woodlitter_gCm2day = output[,,35],
                      HARVESTextracted_som_gCm2day = output[,,36],
                      HARVESTlitter_labile_gCm2day = output[,,37],
                      HARVESTlitter_foliage_gCm2day = output[,,38],
                      HARVESTlitter_roots_gCm2day = output[,,39],
                      HARVESTlitter_wood_gCm2day = output[,,40],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,41],
                      foliage_gCm2 = output[,,42],
                      roots_gCm2 = output[,,43],
                      wood_gCm2 = output[,,44],
                      litter_gCm2 = output[,,45],
                      woodlitter_gCm2 = output[,,46],
                      som_gCm2 = output[,,47],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,48],
                      SurfWater_kgH2Om2 = output[,,49],
                      wSWP_MPa = output[,,50],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,51],
                      cgi = output[,,52],
                      ncce_gCm2day = output[,,53],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,54],
                      gs_mmolH2Om2s = output[,,55],
                      APAR_MJm2day = output[,,56],
                      gb_mmolH2Om2s = output[,,57],
                      CiCa = output[,,58],
                      # Misc
                      RootDepth_m = output[,,59],
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
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P8.R2.#") {
      output_dim = 60 ; MTT_dim = 7 ; SS_dim = 7
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec24",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                              ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2])))
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      rauto_gCm2day = output[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      rhet_som_gCm2day = output[,,4],
                      rhet_woodlitter_gCm2day = output[,,5],
                      fire_gCm2day = output[,,6],
                      harvest_gCm2day = output[,,7],
                      # Internal fluxes
                      alloc_labile_gCm2day = output[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      woodlitter_to_som_gCm2day = output[,,16],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,17],
                      FIRElitter_labile_gCm2day = output[,,18],
                      FIREemiss_foliage_gCm2day = output[,,19],
                      FIRElitter_foliage_gCm2day = output[,,20],
                      FIREemiss_roots_gCm2day = output[,,21],
                      FIRElitter_roots_gCm2day = output[,,22],
                      FIREemiss_wood_gCm2day = output[,,23],
                      FIRElitter_wood_gCm2day = output[,,24],
                      FIREemiss_litter_gCm2day = output[,,25],
                      FIRElitter_litter_gCm2day = output[,,26],
                      FIREemiss_woodlitter_gCm2day = output[,,27],
                      FIRElitter_woodlitter_gCm2day = output[,,28],
                      FIREemiss_som_gCm2day = output[,,29],
                      HARVESTextracted_labile_gCm2day = output[,,30],
                      HARVESTextracted_foliage_gCm2day = output[,,31],
                      HARVESTextracted_roots_gCm2day = output[,,32],
                      HARVESTextracted_wood_gCm2day = output[,,33],
                      HARVESTextracted_litter_gCm2day = output[,,34],
                      HARVESTextracted_woodlitter_gCm2day = output[,,35],
                      HARVESTextracted_som_gCm2day = output[,,36],
                      HARVESTlitter_labile_gCm2day = output[,,37],
                      HARVESTlitter_foliage_gCm2day = output[,,38],
                      HARVESTlitter_roots_gCm2day = output[,,39],
                      HARVESTlitter_wood_gCm2day = output[,,40],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,41],
                      foliage_gCm2 = output[,,42],
                      roots_gCm2 = output[,,43],
                      wood_gCm2 = output[,,44],
                      litter_gCm2 = output[,,45],
                      woodlitter_gCm2 = output[,,46],
                      som_gCm2 = output[,,47],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,48],
                      SurfWater_kgH2Om2 = output[,,49],
                      wSWP_MPa = output[,,50],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,51],
                      cgi = output[,,52],
                      cmi = output[,,53],
                      ncce_gCm2day = output[,,54],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,55],
                      gs_mmolH2Om2s = output[,,56],
                      APAR_MJm2day = output[,,57],
                      gb_mmolH2Om2s = output[,,58],
                      CiCa = output[,,59],
                      # Misc
                      RootDepth_m = output[,,60],
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
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P9.R2.#") {
      output_dim = 60 ; MTT_dim = 7 ; SS_dim = 7
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec25",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                              ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2])))
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      rauto_gCm2day = output[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      rhet_som_gCm2day = output[,,4],
                      rhet_woodlitter_gCm2day = output[,,5],
                      fire_gCm2day = output[,,6],
                      harvest_gCm2day = output[,,7],
                      # Internal fluxes
                      alloc_labile_gCm2day = output[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      woodlitter_to_som_gCm2day = output[,,16],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,17],
                      FIRElitter_labile_gCm2day = output[,,18],
                      FIREemiss_foliage_gCm2day = output[,,19],
                      FIRElitter_foliage_gCm2day = output[,,20],
                      FIREemiss_roots_gCm2day = output[,,21],
                      FIRElitter_roots_gCm2day = output[,,22],
                      FIREemiss_wood_gCm2day = output[,,23],
                      FIRElitter_wood_gCm2day = output[,,24],
                      FIREemiss_litter_gCm2day = output[,,25],
                      FIRElitter_litter_gCm2day = output[,,26],
                      FIREemiss_woodlitter_gCm2day = output[,,27],
                      FIRElitter_woodlitter_gCm2day = output[,,28],
                      FIREemiss_som_gCm2day = output[,,29],
                      HARVESTextracted_labile_gCm2day = output[,,30],
                      HARVESTextracted_foliage_gCm2day = output[,,31],
                      HARVESTextracted_roots_gCm2day = output[,,32],
                      HARVESTextracted_wood_gCm2day = output[,,33],
                      HARVESTextracted_litter_gCm2day = output[,,34],
                      HARVESTextracted_woodlitter_gCm2day = output[,,35],
                      HARVESTextracted_som_gCm2day = output[,,36],
                      HARVESTlitter_labile_gCm2day = output[,,37],
                      HARVESTlitter_foliage_gCm2day = output[,,38],
                      HARVESTlitter_roots_gCm2day = output[,,39],
                      HARVESTlitter_wood_gCm2day = output[,,40],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,41],
                      foliage_gCm2 = output[,,42],
                      roots_gCm2 = output[,,43],
                      wood_gCm2 = output[,,44],
                      litter_gCm2 = output[,,45],
                      woodlitter_gCm2 = output[,,46],
                      som_gCm2 = output[,,47],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,48],
                      SurfWater_kgH2Om2 = output[,,49],
                      wSWP_MPa = output[,,50],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,51],
                      cgi = output[,,52],
                      cmi = output[,,53],
                      ncce_gCm2day = output[,,54],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,55],
                      gs_mmolH2Om2s = output[,,56],
                      APAR_MJm2day = output[,,57],
                      gb_mmolH2Om2s = output[,,58],
                      CiCa = output[,,59],
                      # Misc
                      RootDepth_m = output[,,60],
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

  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P10.R2.#") {
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
                              ,noyears=as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                              ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3)))
                              ,exepath=as.character(crop_file_location),pathlength=as.integer(nchar(crop_file_location)))
      output = tmp$out_var   ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
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
  } else if (model_name == "DALEC.A1.C2.D2.F2.H1.P4.R2.#") {
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
                            ,lat=as.double(lat)
                            ,nopars=as.integer(PROJECT$model$nopars[site])
                            ,nomet=as.integer(dim(met)[2])
                            ,nofluxes=as.integer(PROJECT$model$nofluxes[site])
                            ,nopools=as.integer(PROJECT$model$nopools[site])
                            ,pft=as.integer(pft)
                            ,nodays=as.integer(dim(met)[1])
                            ,noyears=as.integer(noyears)
                            ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1]))))
                            ,nos_iter=as.integer(nos_iter))
    # Extract the different output variables
    output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
    MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
    SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
    # Unload the current dalec shared object
    dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
    rm(tmp) ; gc()
    # create output object
    states_all=list(# Ecosystem fluxes
                    gpp_gCm2day = output[,,1],
                    rauto_gCm2day = output[,,2],
                    rg_foliage_gCm2day = output[,,3],
                    rhet_litter_gCm2day = output[,,4],
                    rhet_som_gCm2day = output[,,5],
                    rhet_woodlitter_gCm2day = output[,,6],
                    fire_gCm2day = output[,,7],
                    harvest_gCm2day = output[,,8],
                    # Internal fluxes
                    alloc_labile_gCm2day = output[,,9],
                    alloc_roots_gCm2day = output[,,10],
                    alloc_wood_gCm2day = output[,,11],
                    labile_to_foliage_gCm2day = output[,,12],
                    foliage_to_litter_gCm2day = output[,,13],
                    roots_to_litter_gCm2day = output[,,14],
                    wood_to_litter_gCm2day = output[,,15],
                    litter_to_som_gCm2day = output[,,16],
                    woodlitter_to_som_gCm2day = output[,,17],
                    # Disturbance fluxes
                    FIREemiss_labile_gCm2day = output[,,18],
                    FIRElitter_labile_gCm2day = output[,,19],
                    FIREemiss_foliage_gCm2day = output[,,20],
                    FIRElitter_foliage_gCm2day = output[,,21],
                    FIREemiss_roots_gCm2day = output[,,22],
                    FIRElitter_roots_gCm2day = output[,,23],
                    FIREemiss_wood_gCm2day = output[,,24],
                    FIRElitter_wood_gCm2day = output[,,25],
                    FIREemiss_litter_gCm2day = output[,,26],
                    FIRElitter_litter_gCm2day = output[,,27],
                    FIREemiss_woodlitter_gCm2day = output[,,28],
                    FIRElitter_woodlitter_gCm2day = output[,,29],
                    FIREemiss_som_gCm2day = output[,,30],
                    HARVESTextracted_labile_gCm2day = output[,,31],
                    HARVESTextracted_foliage_gCm2day = output[,,32],
                    HARVESTextracted_roots_gCm2day = output[,,33],
                    HARVESTextracted_wood_gCm2day = output[,,34],
                    HARVESTextracted_litter_gCm2day = output[,,35],
                    HARVESTextracted_woodlitter_gCm2day = output[,,36],
                    HARVESTextracted_som_gCm2day = output[,,37],
                    HARVESTlitter_labile_gCm2day = output[,,38],
                    HARVESTlitter_foliage_gCm2day = output[,,39],
                    HARVESTlitter_roots_gCm2day = output[,,40],
                    HARVESTlitter_wood_gCm2day = output[,,41],
                    # C pools (gC/m2)
                    labile_gCm2 = output[,,42],
                    foliage_gCm2 = output[,,43],
                    roots_gCm2 = output[,,44],
                    wood_gCm2 = output[,,45],
                    litter_gCm2 = output[,,46],
                    woodlitter_gCm2 = output[,,47],
                    som_gCm2 = output[,,48],
                    # Canopy (phenology) properties
                    lai_m2m2 = output[,,49],
                    gsi = output[,,50],
                    gsi_itemp = output[,,51],
                    gsi_iphoto = output[,,52],
                    gsi_ivpd = output[,,53],
                    # Photosynthesis / C~water coupling related
                    gs_demand_supply_ratio = output[,,54],
                    gs_mmolH2Om2s = output[,,55],
                    APAR_MJm2day = output[,,56],
                    gb_mmolH2Om2s = output[,,57],
                    CiCa = output[,,58],
                    # Misc
                    RootDepth_m = output[,,59],
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
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P3.R1.#") {
      output_dim = 62 ; MTT_dim = 7 ; SS_dim = 7
      # Load the required dalec shared object
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec9",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim)
                                     ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                                     ,fire_dim=as.integer(fire_dim)
                                     ,met=as.double(t(met)),pars=as.double(pars_in)
                                     ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                                     ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                     ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                                     ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site])
                                     ,nomet=as.integer(dim(met)[2]),nofluxes=as.integer(PROJECT$model$nofluxes[site])
                                     ,nopools=as.integer(PROJECT$model$nopools[site]),pft=as.integer(pft)
                                     ,nodays=as.integer(dim(met)[1])
                                     ,noyears=as.integer(noyears)
                                     ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                                     ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                                     ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3))))
      # Extract the different output variables
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      rauto_gCm2day = output[,,2],
                      rg_foliage_gCm2day = output[,,3],
                      rhet_litter_gCm2day = output[,,4],
                      rhet_som_gCm2day = output[,,5],
                      rhet_woodlitter_gCm2day = output[,,6],
                      fire_gCm2day = output[,,7],
                      harvest_gCm2day = output[,,8],
                      # Internal fluxes
                      alloc_labile_gCm2day = output[,,9],
                      alloc_roots_gCm2day = output[,,10],
                      alloc_wood_gCm2day = output[,,11],
                      labile_to_foliage_gCm2day = output[,,12],
                      foliage_to_litter_gCm2day = output[,,13],
                      roots_to_litter_gCm2day = output[,,14],
                      wood_to_litter_gCm2day = output[,,15],
                      litter_to_som_gCm2day = output[,,16],
                      woodlitter_to_som_gCm2day = output[,,17],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,18],
                      FIRElitter_labile_gCm2day = output[,,19],
                      FIREemiss_foliage_gCm2day = output[,,20],
                      FIRElitter_foliage_gCm2day = output[,,21],
                      FIREemiss_roots_gCm2day = output[,,22],
                      FIRElitter_roots_gCm2day = output[,,23],
                      FIREemiss_wood_gCm2day = output[,,24],
                      FIRElitter_wood_gCm2day = output[,,25],
                      FIREemiss_litter_gCm2day = output[,,26],
                      FIRElitter_litter_gCm2day = output[,,27],
                      FIREemiss_woodlitter_gCm2day = output[,,28],
                      FIRElitter_woodlitter_gCm2day = output[,,29],
                      FIREemiss_som_gCm2day = output[,,30],
                      HARVESTextracted_labile_gCm2day = output[,,31],
                      HARVESTextracted_foliage_gCm2day = output[,,32],
                      HARVESTextracted_roots_gCm2day = output[,,33],
                      HARVESTextracted_wood_gCm2day = output[,,34],
                      HARVESTextracted_litter_gCm2day = output[,,35],
                      HARVESTextracted_woodlitter_gCm2day = output[,,36],
                      HARVESTextracted_som_gCm2day = output[,,37],
                      HARVESTlitter_labile_gCm2day = output[,,38],
                      HARVESTlitter_foliage_gCm2day = output[,,39],
                      HARVESTlitter_roots_gCm2day = output[,,40],
                      HARVESTlitter_wood_gCm2day = output[,,41],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,42],
                      foliage_gCm2 = output[,,43],
                      roots_gCm2 = output[,,44],
                      wood_gCm2 = output[,,45],
                      litter_gCm2 = output[,,46],
                      woodlitter_gCm2 = output[,,47],
                      som_gCm2 = output[,,48],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,49],
                      SurfWater_kgH2Om2 = output[,,50],
                      wSWP_MPa = output[,,51],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,52],
                      gsi = output[,,53],
                      gsi_itemp = output[,,54],
                      gsi_iphoto = output[,,55],
                      gsi_ivpd = output[,,56],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,57],
                      gs_mmolH2Om2s = output[,,58],
                      APAR_MJm2day = output[,,59],
                      gb_mmolH2Om2s = output[,,60],
                      CiCa = output[,,61],
                      # Misc
                      RootDepth_m = output[,,62],
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
  } else if (model_name == "DALEC.C1.D1.F2.P1.#") {
      output_dim = 44 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec2",output_dim=as.integer(output_dim)
                                          ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                                          ,met=as.double(t(met))
                                          ,pars=as.double(pars_in)
                                          ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                                          ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                          ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                                          ,lat=as.double(lat)
                                          ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                          ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                          ,nodays=as.integer(dim(met)[1])
                                          ,noyears=as.integer(noyears)
                                          ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter) )
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      rauto_gCm2day = output[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      rhet_som_gCm2day = output[,,4],
                      fire_gCm2day = output[,,5],
                      harvest_gCm2day = output[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      foliage_gCm2 = output[,,38],
                      roots_gCm2 = output[,,39],
                      wood_gCm2 = output[,,40],
                      litter_gCm2 = output[,,41],
                      som_gCm2 = output[,,42],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,43],
                      # Photosynthesis / C~water coupling related
                      CiCa = output[,,44],
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
                                ,lat=as.double(lat)
                                ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                ,nodays=as.integer(dim(met)[1])
                                ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
      # extract output variables
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      rauto_gCm2day = output[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      rhet_som_gCm2day = output[,,4],
                      fire_gCm2day = output[,,5],
                      harvest_gCm2day = output[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      foliage_gCm2 = output[,,38],
                      roots_gCm2 = output[,,39],
                      wood_gCm2 = output[,,40],
                      litter_gCm2 = output[,,41],
                      som_gCm2 = output[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      SurfWater_kgH2Om2 = output[,,44],
                      # wSWP_MPa = output[,,45],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,45],
                      # Photosynthesis / C~water coupling related
                      CiCa = output[,,46],
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
                                 ,lat=as.double(lat)
                                 ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                 ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                 ,nodays=as.integer(dim(met)[1])
                                 ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
      # extract output variables
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      rauto_gCm2day = output[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      rhet_som_gCm2day = output[,,4],
                      fire_gCm2day = output[,,5],
                      harvest_gCm2day = output[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      foliage_gCm2 = output[,,38],
                      roots_gCm2 = output[,,39],
                      wood_gCm2 = output[,,40],
                      litter_gCm2 = output[,,41],
                      som_gCm2 = output[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      SurfWater_kgH2Om2 = output[,,44],
                      # wSWP_MPa = output[,,45],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,45],
                      # Photosynthesis / C~water coupling related
                      CiCa = output[,,46],
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
  } else if (model_name == "DALEC.A1.C1.D2.F2.H1.P1.#") {
      output_dim = 48 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec3",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodays=as.integer(dim(met)[1])
                             ,noyears=as.integer(noyears)
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                             ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                             ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3))))
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      rauto_gCm2day = output[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      rhet_som_gCm2day = output[,,4],
                      fire_gCm2day = output[,,5],
                      harvest_gCm2day = output[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      foliage_gCm2 = output[,,38],
                      roots_gCm2 = output[,,39],
                      wood_gCm2 = output[,,40],
                      litter_gCm2 = output[,,41],
                      som_gCm2 = output[,,42],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,43],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,44],
                      gs_mmolH2Om2s = output[,,45],
                      APAR_MJm2day = output[,,46],
                      gb_mmolH2Om2s = output[,,47],
                      CiCa = output[,,48],
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
  } else if (model_name == "DALEC.A1.C1.D2.F2.H2.P1.#") {
      output_dim = 58 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec4",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodays=as.integer(dim(met)[1])
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                             ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                             ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2])))
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      rauto_gCm2day = output[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      rhet_som_gCm2day = output[,,4],
                      fire_gCm2day = output[,,5],
                      harvest_gCm2day = output[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      foliage_gCm2 = output[,,38],
                      roots_gCm2 = output[,,39],
                      wood_gCm2 = output[,,40],
                      litter_gCm2 = output[,,41],
                      som_gCm2 = output[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      Etrans_kgH2Om2day = output[,,44],
                      Esoil_kgH2Om2day = output[,,45],
                      Ewetcanopy_kgH2Om2day = output[,,46],
                      runoff_kgH2Om2day = output[,,47],
                      underflow_kgH2Om2day = output[,,48],
                      SurfWater_kgH2Om2 = output[,,49],
                      wSWP_MPa = output[,,50],
                      snow_kgH2Om2 = output[,,51],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,52],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,53],
                      gs_mmolH2Om2s = output[,,54],
                      APAR_MJm2day = output[,,55],
                      gb_mmolH2Om2s = output[,,56],
                      CiCa = output[,,57],
                      # Misc
                      RootDepth_m = output[,,58],
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
  } else if (model_name == "DALEC.A2.C1.D2.F2.H2.P1.#") {
      output_dim = 58 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec20",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                              ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2])))
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      rauto_gCm2day = output[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      rhet_som_gCm2day = output[,,4],
                      fire_gCm2day = output[,,5],
                      harvest_gCm2day = output[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      foliage_gCm2 = output[,,38],
                      roots_gCm2 = output[,,39],
                      wood_gCm2 = output[,,40],
                      litter_gCm2 = output[,,41],
                      som_gCm2 = output[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      Etrans_kgH2Om2day = output[,,44],
                      Esoil_kgH2Om2day = output[,,45],
                      Ewetcanopy_kgH2Om2day = output[,,46],
                      runoff_kgH2Om2day = output[,,47],
                      underflow_kgH2Om2day = output[,,48],
                      SurfWater_kgH2Om2 = output[,,49],
                      wSWP_MPa = output[,,50],
                      snow_kgH2Om2 = output[,,51],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,52],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,53],
                      gs_mmolH2Om2s = output[,,54],
                      APAR_MJm2day = output[,,55],
                      gb_mmolH2Om2s = output[,,56],
                      CiCa = output[,,57],
                      # Misc
                      RootDepth_m = output[,,58],
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
  } else if (model_name == "DALEC.A1.C1.D2.F2.H2.P2.#") {
      output_dim = 51 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec18",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,noyears=as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                              ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3))))
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      rauto_gCm2day = output[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      rhet_som_gCm2day = output[,,4],
                      fire_gCm2day = output[,,5],
                      harvest_gCm2day = output[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      foliage_gCm2 = output[,,38],
                      roots_gCm2 = output[,,39],
                      wood_gCm2 = output[,,40],
                      litter_gCm2 = output[,,41],
                      som_gCm2 = output[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      SurfWater_kgH2Om2 = output[,,44],
                      wSWP_MPa = output[,,45],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,46],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,47],
                      gs_mmolH2Om2s = output[,,48],
                      APAR_MJm2day = output[,,49],
                      gb_mmolH2Om2s = output[,,50],
                      CiCa = output[,,51],
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
  } else if (model_name == "DALEC.A1.C1.D2.F2.H2.P5.#") {
      output_dim = 52 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec21",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,noyears=as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                              ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3))))
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      rauto_gCm2day = output[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      rhet_som_gCm2day = output[,,4],
                      fire_gCm2day = output[,,5],
                      harvest_gCm2day = output[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      foliage_gCm2 = output[,,38],
                      roots_gCm2 = output[,,39],
                      wood_gCm2 = output[,,40],
                      litter_gCm2 = output[,,41],
                      som_gCm2 = output[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      SurfWater_kgH2Om2 = output[,,44],
                      wSWP_MPa = output[,,45],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,46],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,47],
                      gs_mmolH2Om2s = output[,,48],
                      APAR_MJm2day = output[,,49],
                      gb_mmolH2Om2s = output[,,50],
                      CiCa = output[,,51],
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
  } else if (model_name == "DALEC.A1.C1.D2.F2.H2.P6.#") {
      output_dim = 51 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec22",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,noyears=as.integer(noyears)
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                              ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3))))
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      rauto_gCm2day = output[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      rhet_som_gCm2day = output[,,4],
                      fire_gCm2day = output[,,5],
                      harvest_gCm2day = output[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      foliage_gCm2 = output[,,38],
                      roots_gCm2 = output[,,39],
                      wood_gCm2 = output[,,40],
                      litter_gCm2 = output[,,41],
                      som_gCm2 = output[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      SurfWater_kgH2Om2 = output[,,44],
                      wSWP_MPa = output[,,45],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,46],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,47],
                      gs_mmolH2Om2s = output[,,48],
                      APAR_MJm2day = output[,,49],
                      gb_mmolH2Om2s = output[,,50],
                      CiCa = output[,,51],
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
  } else if (model_name == "DALEC.A1.C1.D2.F2.H2.P1.R1.#") {
      output_dim = 51 ; MTT_dim = 6 ; SS_dim = 6
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec5",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodays=as.integer(dim(met)[1])
                             ,noyears=as.integer(noyears)
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                             ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                             ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3))))
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      rauto_gCm2day = output[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      rhet_som_gCm2day = output[,,4],
                      fire_gCm2day = output[,,5],
                      harvest_gCm2day = output[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      alloc_labile_gCm2day = output[,,8],
                      alloc_roots_gCm2day = output[,,9],
                      alloc_wood_gCm2day = output[,,10],
                      labile_to_foliage_gCm2day = output[,,11],
                      foliage_to_litter_gCm2day = output[,,12],
                      roots_to_litter_gCm2day = output[,,13],
                      wood_to_litter_gCm2day = output[,,14],
                      litter_to_som_gCm2day = output[,,15],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,16],
                      FIRElitter_labile_gCm2day = output[,,17],
                      FIREemiss_foliage_gCm2day = output[,,18],
                      FIRElitter_foliage_gCm2day = output[,,19],
                      FIREemiss_roots_gCm2day = output[,,20],
                      FIRElitter_roots_gCm2day = output[,,21],
                      FIREemiss_wood_gCm2day = output[,,22],
                      FIRElitter_wood_gCm2day = output[,,23],
                      FIREemiss_litter_gCm2day = output[,,24],
                      FIRElitter_litter_gCm2day = output[,,25],
                      FIREemiss_som_gCm2day = output[,,26],
                      HARVESTextracted_labile_gCm2day = output[,,27],
                      HARVESTextracted_foliage_gCm2day = output[,,28],
                      HARVESTextracted_roots_gCm2day = output[,,29],
                      HARVESTextracted_wood_gCm2day = output[,,30],
                      HARVESTextracted_litter_gCm2day = output[,,31],
                      HARVESTextracted_som_gCm2day = output[,,32],
                      HARVESTlitter_labile_gCm2day = output[,,33],
                      HARVESTlitter_foliage_gCm2day = output[,,34],
                      HARVESTlitter_roots_gCm2day = output[,,35],
                      HARVESTlitter_wood_gCm2day = output[,,36],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,37],
                      foliage_gCm2 = output[,,38],
                      roots_gCm2 = output[,,39],
                      wood_gCm2 = output[,,40],
                      litter_gCm2 = output[,,41],
                      som_gCm2 = output[,,42],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,43],
                      SurfWater_kgH2Om2 = output[,,44],
                      wSWP_MPa = output[,,45],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,46],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,47],
                      gs_mmolH2Om2s = output[,,48],
                      APAR_MJm2day = output[,,49],
                      gb_mmolH2Om2s = output[,,50],
                      CiCa = output[,,51],
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
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P1.R1.#") {
      output_dim = 57 ; MTT_dim = 7 ; SS_dim = 7
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec6",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodays=as.integer(dim(met)[1])
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                             ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                             ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2])))
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      rauto_gCm2day = output[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      rhet_som_gCm2day = output[,,4],
                      rhet_woodlitter_gCm2day = output[,,5],
                      fire_gCm2day = output[,,6],
                      harvest_gCm2day = output[,,7],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,8],
                      alloc_labile_gCm2day = output[,,9],
                      alloc_roots_gCm2day = output[,,10],
                      alloc_wood_gCm2day = output[,,11],
                      labile_to_foliage_gCm2day = output[,,12],
                      foliage_to_litter_gCm2day = output[,,13],
                      roots_to_litter_gCm2day = output[,,14],
                      wood_to_litter_gCm2day = output[,,15],
                      litter_to_som_gCm2day = output[,,16],
                      woodlitter_to_som_gCm2day = output[,,17],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,18],
                      FIRElitter_labile_gCm2day = output[,,19],
                      FIREemiss_foliage_gCm2day = output[,,20],
                      FIRElitter_foliage_gCm2day = output[,,21],
                      FIREemiss_roots_gCm2day = output[,,22],
                      FIRElitter_roots_gCm2day = output[,,23],
                      FIREemiss_wood_gCm2day = output[,,24],
                      FIRElitter_wood_gCm2day = output[,,25],
                      FIREemiss_litter_gCm2day = output[,,26],
                      FIRElitter_litter_gCm2day = output[,,27],
                      FIREemiss_woodlitter_gCm2day = output[,,28],
                      FIRElitter_woodlitter_gCm2day = output[,,29],
                      FIREemiss_som_gCm2day = output[,,30],
                      HARVESTextracted_labile_gCm2day = output[,,31],
                      HARVESTextracted_foliage_gCm2day = output[,,32],
                      HARVESTextracted_roots_gCm2day = output[,,33],
                      HARVESTextracted_wood_gCm2day = output[,,34],
                      HARVESTextracted_litter_gCm2day = output[,,35],
                      HARVESTextracted_woodlitter_gCm2day = output[,,36],
                      HARVESTextracted_som_gCm2day = output[,,37],
                      HARVESTlitter_labile_gCm2day = output[,,38],
                      HARVESTlitter_foliage_gCm2day = output[,,39],
                      HARVESTlitter_roots_gCm2day = output[,,40],
                      HARVESTlitter_wood_gCm2day = output[,,41],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,42],
                      foliage_gCm2 = output[,,43],
                      roots_gCm2 = output[,,44],
                      wood_gCm2 = output[,,45],
                      litter_gCm2 = output[,,46],
                      woodlitter_gCm2 = output[,,47],
                      som_gCm2 = output[,,48],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,49],
                      SurfWater_kgH2Om2 = output[,,50],
                      wSWP_MPa = output[,,51],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,52],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,53],
                      gs_mmolH2Om2s = output[,,54],
                      APAR_MJm2day = output[,,55],
                      gb_mmolH2Om2s = output[,,56],
                      CiCa = output[,,57],
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
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P2.R1.#") {
      output_dim = 57 ; MTT_dim = 7 ; SS_dim = 7
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec7",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodays=as.integer(dim(met)[1])
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                             ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                             ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2])))
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      rauto_gCm2day = output[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      rhet_som_gCm2day = output[,,4],
                      rhet_woodlitter_gCm2day = output[,,5],
                      fire_gCm2day = output[,,6],
                      harvest_gCm2day = output[,,7],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,8],
                      alloc_labile_gCm2day = output[,,9],
                      alloc_roots_gCm2day = output[,,10],
                      alloc_wood_gCm2day = output[,,11],
                      labile_to_foliage_gCm2day = output[,,12],
                      foliage_to_litter_gCm2day = output[,,13],
                      roots_to_litter_gCm2day = output[,,14],
                      wood_to_litter_gCm2day = output[,,15],
                      litter_to_som_gCm2day = output[,,16],
                      woodlitter_to_som_gCm2day = output[,,17],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,18],
                      FIRElitter_labile_gCm2day = output[,,19],
                      FIREemiss_foliage_gCm2day = output[,,20],
                      FIRElitter_foliage_gCm2day = output[,,21],
                      FIREemiss_roots_gCm2day = output[,,22],
                      FIRElitter_roots_gCm2day = output[,,23],
                      FIREemiss_wood_gCm2day = output[,,24],
                      FIRElitter_wood_gCm2day = output[,,25],
                      FIREemiss_litter_gCm2day = output[,,26],
                      FIRElitter_litter_gCm2day = output[,,27],
                      FIREemiss_woodlitter_gCm2day = output[,,28],
                      FIRElitter_woodlitter_gCm2day = output[,,29],
                      FIREemiss_som_gCm2day = output[,,30],
                      HARVESTextracted_labile_gCm2day = output[,,31],
                      HARVESTextracted_foliage_gCm2day = output[,,32],
                      HARVESTextracted_roots_gCm2day = output[,,33],
                      HARVESTextracted_wood_gCm2day = output[,,34],
                      HARVESTextracted_litter_gCm2day = output[,,35],
                      HARVESTextracted_woodlitter_gCm2day = output[,,36],
                      HARVESTextracted_som_gCm2day = output[,,37],
                      HARVESTlitter_labile_gCm2day = output[,,38],
                      HARVESTlitter_foliage_gCm2day = output[,,39],
                      HARVESTlitter_roots_gCm2day = output[,,40],
                      HARVESTlitter_wood_gCm2day = output[,,41],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,42],
                      foliage_gCm2 = output[,,43],
                      roots_gCm2 = output[,,44],
                      wood_gCm2 = output[,,45],
                      litter_gCm2 = output[,,46],
                      woodlitter_gCm2 = output[,,47],
                      som_gCm2 = output[,,48],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,49],
                      SurfWater_kgH2Om2 = output[,,50],
                      wSWP_MPa = output[,,51],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,52],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,53],
                      gs_mmolH2Om2s = output[,,54],
                      APAR_MJm2day = output[,,55],
                      gb_mmolH2Om2s = output[,,56],
                      CiCa = output[,,57],
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
  } else if (model_name == "DALEC.A1.C2.D2.F2.H2.P2.R3.#") {
      output_dim = 57 ; MTT_dim = 7 ; SS_dim = 7
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec19",output_dim=as.integer(output_dim)
                              ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                              ,met=as.double(t(met))
                              ,pars=as.double(pars_in)
                              ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                              ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                              ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                              ,lat=as.double(lat)
                              ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                              ,nodays=as.integer(dim(met)[1])
                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                              ,soil_frac_clay_in=as.double(c(soil_info[3],soil_info[4],soil_info[4]))
                              ,soil_frac_sand_in=as.double(c(soil_info[1],soil_info[2],soil_info[2])))
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      rauto_gCm2day = output[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      rhet_som_gCm2day = output[,,4],
                      rhet_woodlitter_gCm2day = output[,,5],
                      fire_gCm2day = output[,,6],
                      harvest_gCm2day = output[,,7],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,8],
                      alloc_labile_gCm2day = output[,,9],
                      alloc_roots_gCm2day = output[,,10],
                      alloc_wood_gCm2day = output[,,11],
                      labile_to_foliage_gCm2day = output[,,12],
                      foliage_to_litter_gCm2day = output[,,13],
                      roots_to_litter_gCm2day = output[,,14],
                      wood_to_litter_gCm2day = output[,,15],
                      litter_to_som_gCm2day = output[,,16],
                      woodlitter_to_som_gCm2day = output[,,17],
                      # Disturbance fluxes
                      FIREemiss_labile_gCm2day = output[,,18],
                      FIRElitter_labile_gCm2day = output[,,19],
                      FIREemiss_foliage_gCm2day = output[,,20],
                      FIRElitter_foliage_gCm2day = output[,,21],
                      FIREemiss_roots_gCm2day = output[,,22],
                      FIRElitter_roots_gCm2day = output[,,23],
                      FIREemiss_wood_gCm2day = output[,,24],
                      FIRElitter_wood_gCm2day = output[,,25],
                      FIREemiss_litter_gCm2day = output[,,26],
                      FIRElitter_litter_gCm2day = output[,,27],
                      FIREemiss_woodlitter_gCm2day = output[,,28],
                      FIRElitter_woodlitter_gCm2day = output[,,29],
                      FIREemiss_som_gCm2day = output[,,30],
                      HARVESTextracted_labile_gCm2day = output[,,31],
                      HARVESTextracted_foliage_gCm2day = output[,,32],
                      HARVESTextracted_roots_gCm2day = output[,,33],
                      HARVESTextracted_wood_gCm2day = output[,,34],
                      HARVESTextracted_litter_gCm2day = output[,,35],
                      HARVESTextracted_woodlitter_gCm2day = output[,,36],
                      HARVESTextracted_som_gCm2day = output[,,37],
                      HARVESTlitter_labile_gCm2day = output[,,38],
                      HARVESTlitter_foliage_gCm2day = output[,,39],
                      HARVESTlitter_roots_gCm2day = output[,,40],
                      HARVESTlitter_wood_gCm2day = output[,,41],
                      # C pools (gC/m2)
                      labile_gCm2 = output[,,42],
                      foliage_gCm2 = output[,,43],
                      roots_gCm2 = output[,,44],
                      wood_gCm2 = output[,,45],
                      litter_gCm2 = output[,,46],
                      woodlitter_gCm2 = output[,,47],
                      som_gCm2 = output[,,48],
                      # Water cycle related
                      ET_kgH2Om2day = output[,,49],
                      SurfWater_kgH2Om2 = output[,,50],
                      wSWP_MPa = output[,,51],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,52],
                      # Photosynthesis / C~water coupling related
                      gs_demand_supply_ratio = output[,,53],
                      gs_mmolH2Om2s = output[,,54],
                      APAR_MJm2day = output[,,55],
                      gb_mmolH2Om2s = output[,,56],
                      CiCa = output[,,57],
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
  } else if (model_name == "DALEC.C5.D1.F2.P1.#") {
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
                            ,lat=as.double(lat)
                            ,nopars=as.integer(PROJECT$model$nopars[site])
                            ,nomet=as.integer(dim(met)[2])
                            ,nofluxes=as.integer(PROJECT$model$nofluxes[site])
                            ,nopools=as.integer(PROJECT$model$nopools[site])
                            ,nodays=as.integer(dim(met)[1])
                            ,noyears=as.integer(noyears)
                            ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1]))))
                            ,nos_iter=as.integer(nos_iter) )
    output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
    MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
    SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
    # Unload the current dalec shared object
    dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
    rm(tmp) ; gc()
    # create output object
    states_all=list(# Ecosystem fluxes
                    gpp_gCm2day = output[,,1],
                    rauto_gCm2day = output[,,2],
                    rhet_dom_gCm2day = output[,,3],
                    fire_gCm2day = output[,,4],
                    # Internal fluxes
                    alloc_foliage_gCm2day = output[,,5],
                    alloc_labile_gCm2day = output[,,6],
                    alloc_roots_wood_gCm2day = output[,,7],
                    labile_to_foliage_gCm2day = output[,,8],
                    foliage_to_litter_gCm2day = output[,,9],
                    roots_wood_to_litter_gCm2day = output[,,10],
                    # Disturbance fluxes
                    FIREemiss_labile_gCm2day = output[,,11],
                    FIRElitter_labile_gCm2day = output[,,12],
                    FIREemiss_foliage_gCm2day = output[,,13],
                    FIRElitter_foliage_gCm2day = output[,,14],
                    FIREemiss_roots_wood_gCm2day = output[,,15],
                    FIRElitter_roots_wood_gCm2day = output[,,16],
                    FIREemiss_dom_gCm2day = output[,,17],
                    # C pools (gC/m2)
                    labile_gCm2 = output[,,18],
                    foliage_gCm2 = output[,,19],
                    roots_wood_gCm2 = output[,,20],
                    dom_gCm2 = output[,,21],
                    # Canopy (phenology) properties
                    lai_m2m2 = output[,,22],
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
  } else if (model_name == "DALEC.D1.F2.#") {
      output_dim = 36 ; MTT_dim = 5 ; SS_dim = 5
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      tmp=.Fortran( "rdalec1",output_dim=as.integer(output_dim)
                             ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                             ,met=as.double(t(met))
                             ,pars=as.double(pars_in)
                             ,out_var1=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                             ,out_var2=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                             ,out_var3=as.double(array(0,dim=c(nos_iter,SS_dim)))
                             ,lat=as.double(lat)
                             ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                             ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                             ,nodays=as.integer(dim(met)[1])
                             ,noyears=as.integer(noyears)
                             ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
      output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
      SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
      # Unload the current dalec shared object
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # Create the output object
      states_all=list(# Ecosystem fluxes
                      gpp_gCm2day = output[,,1],
                      rauto_gCm2day = output[,,2],
                      rhet_litter_gCm2day = output[,,3],
                      rhet_som_gCm2day = output[,,4],
                      fire_gCm2day = output[,,5],
                      harvest_gCm2day = output[,,6],
                      # Internal fluxes
                      alloc_foliage_gCm2day = output[,,7],
                      alloc_roots_gCm2day = output[,,8],
                      alloc_wood_gCm2day = output[,,9],
                      foliage_to_litter_gCm2day = output[,,10],
                      roots_to_litter_gCm2day = output[,,11],
                      wood_to_litter_gCm2day = output[,,12],
                      litter_to_som_gCm2day = output[,,13],
                      # Disturbance fluxes
                      FIREemiss_foliage_gCm2day = output[,,14],
                      FIRElitter_foliage_gCm2day = output[,,15],
                      FIREemiss_roots_gCm2day = output[,,16],
                      FIRElitter_roots_gCm2day = output[,,17],
                      FIREemiss_wood_gCm2day = output[,,18],
                      FIRElitter_wood_gCm2day = output[,,19],
                      FIREemiss_litter_gCm2day = output[,,20],
                      FIRElitter_litter_gCm2day = output[,,21],
                      FIREemiss_som_gCm2day = output[,,22],
                      HARVESTextracted_foliage_gCm2day = output[,,23],
                      HARVESTextracted_roots_gCm2day = output[,,24],
                      HARVESTextracted_wood_gCm2day = output[,,25],
                      HARVESTextracted_litter_gCm2day = output[,,26],
                      HARVESTextracted_som_gCm2day = output[,,27],
                      HARVESTlitter_foliage_gCm2day = output[,,28],
                      HARVESTlitter_roots_gCm2day = output[,,29],
                      HARVESTlitter_wood_gCm2day = output[,,30],
                      # C pools (gC/m2)
                      foliage_gCm2 = output[,,31],
                      roots_gCm2 = output[,,32],
                      wood_gCm2 = output[,,33],
                      litter_gCm2 = output[,,34],
                      som_gCm2 = output[,,35],
                      # Canopy (phenology) properties
                      lai_m2m2 = output[,,36],
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
  } else if (model_name == "DALEC.C4.D1.F2.#") {
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
                            ,lat=as.double(lat)
                            ,nopars=as.integer(PROJECT$model$nopars[site])
                            ,nomet=as.integer(dim(met)[2])
                            ,nofluxes=as.integer(PROJECT$model$nofluxes[site])
                            ,nopools=as.integer(PROJECT$model$nopools[site])
                            ,nodays=as.integer(dim(met)[1])
                            ,noyears=as.integer(noyears)
                            ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1]))))
                            ,nos_iter=as.integer(nos_iter) )
    output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
    MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
    SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
    # Unload the current dalec shared object
    dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
    rm(tmp) ; gc()
    # create output object
    states_all=list(# Ecosystem fluxes
                    gpp_gCm2day = output[,,1],
                    rauto_gCm2day = output[,,2],
                    rhet_dom_gCm2day = output[,,3],
                    fire_gCm2day = output[,,4],
                    harvest_gCm2day = output[,,5],
                    # Internal fluxes
                    alloc_foliage_gCm2day = output[,,6],
                    alloc_roots_wood_gCm2day = output[,,7],
                    foliage_to_litter_gCm2day = output[,,8],
                    roots_wood_to_litter_gCm2day = output[,,9],
                    # Disturbance fluxes
                    FIREemiss_foliage_gCm2day = output[,,10],
                    FIRElitter_foliage_gCm2day = output[,,11],
                    FIREemiss_roots_wood_gCm2day = output[,,12],
                    FIRElitter_roots_wood_gCm2day = output[,,13],
                    FIREemiss_dom_gCm2day = output[,,14],
                    HARVESTextracted_foliage_gCm2day = output[,,15],
                    HARVESTextracted_roots_wood_gCm2day = output[,,16],
                    HARVESTextracted_som_gCm2day = output[,,17],
                    HARVESTlitter_foliage_gCm2day = output[,,18],
                    HARVESTlitter_roots_wood_gCm2day = output[,,19],
                    # C pools (gC/m2)
                    foliage_gCm2 = output[,,20],
                    roots_wood_gCm2 = output[,,21],
                    dom_gCm2 = output[,,22],
                    # Canopy (phenology) properties
                    lai_m2m2 = output[,,23],
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
  } else if (model_name == "DALEC.A1.C2.D2.F2.H1.P3.R1.#") {
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
                           ,lat=as.double(lat)
                           ,nopars=as.integer(PROJECT$model$nopars[site])
                           ,nomet=as.integer(dim(met)[2])
                           ,nofluxes=as.integer(PROJECT$model$nofluxes[site])
                           ,nopools=as.integer(PROJECT$model$nopools[site])
                           ,pft=as.integer(pft)
                           ,nodays=as.integer(dim(met)[1])
                           ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1]))))
                           ,nos_iter=as.integer(nos_iter))
    # Extract the different output variables
    output = tmp$out_var1    ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
    MTT_years = tmp$out_var2 ; MTT_years = array(MTT_years, dim=c(nos_iter,MTT_dim))
    SS_gCm2 = tmp$out_var3   ; SS_gCm2 = array(SS_gCm2, dim=c(nos_iter,SS_dim))
    # Unload the current dalec shared object
    dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
    rm(tmp) ; gc()
    # create output object
    states_all=list(# Ecosystem fluxes
                    gpp_gCm2day = output[,,1],
                    rauto_gCm2day = output[,,2],
                    rg_foliage_gCm2day = output[,,3],
                    rhet_litter_gCm2day = output[,,4],
                    rhet_som_gCm2day = output[,,5],
                    rhet_woodlitter_gCm2day = output[,,6],
                    fire_gCm2day = output[,,7],
                    harvest_gCm2day = output[,,8],
                    # Internal fluxes
                    alloc_labile_gCm2day = output[,,9],
                    alloc_roots_gCm2day = output[,,10],
                    alloc_wood_gCm2day = output[,,11],
                    labile_to_foliage_gCm2day = output[,,12],
                    foliage_to_litter_gCm2day = output[,,13],
                    roots_to_litter_gCm2day = output[,,14],
                    wood_to_litter_gCm2day = output[,,15],
                    litter_to_som_gCm2day = output[,,16],
                    woodlitter_to_som_gCm2day = output[,,17],
                    # Disturbance fluxes
                    FIREemiss_labile_gCm2day = output[,,18],
                    FIRElitter_labile_gCm2day = output[,,19],
                    FIREemiss_foliage_gCm2day = output[,,20],
                    FIRElitter_foliage_gCm2day = output[,,21],
                    FIREemiss_roots_gCm2day = output[,,22],
                    FIRElitter_roots_gCm2day = output[,,23],
                    FIREemiss_wood_gCm2day = output[,,24],
                    FIRElitter_wood_gCm2day = output[,,25],
                    FIREemiss_litter_gCm2day = output[,,26],
                    FIRElitter_litter_gCm2day = output[,,27],
                    FIREemiss_woodlitter_gCm2day = output[,,28],
                    FIRElitter_woodlitter_gCm2day = output[,,29],
                    FIREemiss_som_gCm2day = output[,,30],
                    HARVESTextracted_labile_gCm2day = output[,,31],
                    HARVESTextracted_foliage_gCm2day = output[,,32],
                    HARVESTextracted_roots_gCm2day = output[,,33],
                    HARVESTextracted_wood_gCm2day = output[,,34],
                    HARVESTextracted_litter_gCm2day = output[,,35],
                    HARVESTextracted_woodlitter_gCm2day = output[,,36],
                    HARVESTextracted_som_gCm2day = output[,,37],
                    HARVESTlitter_labile_gCm2day = output[,,38],
                    HARVESTlitter_foliage_gCm2day = output[,,39],
                    HARVESTlitter_roots_gCm2day = output[,,40],
                    HARVESTlitter_wood_gCm2day = output[,,41],
                    # C pools (gC/m2)
                    labile_gCm2 = output[,,42],
                    foliage_gCm2 = output[,,43],
                    roots_gCm2 = output[,,44],
                    wood_gCm2 = output[,,45],
                    litter_gCm2 = output[,,46],
                    woodlitter_gCm2 = output[,,47],
                    som_gCm2 = output[,,48],
                    # Canopy (phenology) properties
                    lai_m2m2 = output[,,49],
                    gsi = output[,,50],
                    gsi_itemp = output[,,51],
                    gsi_iphoto = output[,,52],
                    gsi_ivpd = output[,,53],
                    # Photosynthesis / C~water coupling related
                    gs_demand_supply_ratio = output[,,54],
                    gs_mmolH2Om2s = output[,,55],
                    APAR_MJm2day = output[,,56],
                    gb_mmolH2Om2s = output[,,57],
                    CiCa = output[,,58],
                    # Misc
                    RootDepth_m = output[,,59],
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
  } else {
      stop(paste("Model choice (",model_name,") does not have corresponding R interface",sep=""))
  }

  # return state variable means
  return(states_all) ; gc(verbose=FALSE)

} # end of function
## Use byte compile
simulate_all<-cmpfun(simulate_all)
