
###
## Function to generate mean state variable information by running the parameters and model choice
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

simulate_all<- function (site,PROJECT,model_name,met,pars,lat,pft,parameter_type,exepath,soil_info) {

  output_dim=17 ; aNPP_dim = 3 ; MTT_dim = 5 ; SS_dim = 5
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
                          ,soil_frac_clay=as.double(array(c(soil_info[3],soil_info[3],soil_info[4],soil_info[4]),dim=c(4)))
                          ,soil_frac_sand=as.double(array(c(soil_info[1],soil_info[1],soil_info[2],soil_info[2]),dim=c(4))) )
      output=tmp$out_var
      output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # construct output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      evap_kgH2Om2day = output[,,3],soilevap_kgH2Om2day = output[,,4],
                      Rtot_MPasm2mmol = output[,,5],wetcanevap_kgH2Om2day = output[,,6],
                      gs_demand_supply = output[,,7], gs_total_canopy = output[,,8],
                      APAR_MJm2day = output[,,9], gb_total_canopy = output[,,10],
                      CiCa = output[,,11])
  } else if (model_name == "DALEC_BUCKET") {
      output_dim=28
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      crop_file_location=paste(PROJECT$exepath,"winter_wheat_development.csv", sep="")
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdalecbucket",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim)
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
                                     ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[3],soil_info[4],soil_info[4]),dim=c(4)))
                                     ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[1],soil_info[2],soil_info[2]),dim=c(4)))
                                     ,exepath=as.character(crop_file_location),pathlength=as.integer(nchar(crop_file_location)))
      output = tmp$out_var ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP = tmp$out_var2  ; aNPP = array(aNPP, dim=c(nos_iter,aNPP_dim))
      MTT = tmp$out_var3   ; MTT = array(MTT, dim=c(nos_iter,MTT_dim))
      SS = tmp$out_var4    ; SS = array(SS, dim=c(nos_iter,SS_dim))
      aMTT = tmp$out_var5  ; aMTT = array(aMTT, dim=c(nos_iter,MTT_dim,noyears))
      MTTnat = tmp$out_var6; MTTnat = array(MTTnat, dim=c(nos_iter,MTT_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], lit_gCm2 = output[,,10],
                      lab_gCm2 = output[,,11], fol_gCm2 = output[,,12],
                      harvest_C_gCm2day = output[,,13], gsi = output[,,14],
                      gsi_itemp = output[,,15], gsi_iphoto = output[,,16],
                      gsi_ivpd = output[,,17], evap_kgH2Om2day = output[,,18],
                      sfc_water_mm = output[,,19], wSWP_MPa = output[,,20],
                      litwood_gCm2 = output[,,21], fire_gCm2day = output[,,23],
                      gs_demand_supply = output[,,24], gs_total_canopy = output[,,25],
                      APAR_MJm2day = output[,,26], gb_total_canopy = output[,,27],
                      CiCa = output[,,28],
                      aNPP = aNPP, MTT = MTT, SS = SS, aMTT = aMTT, natMTT = MTTnat)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  } else if (model_name == "DALEC_G5") {
      output_dim=28
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      crop_file_location=paste(PROJECT$exepath,"winter_wheat_development.csv", sep="")
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdalecg5",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim)
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
                                     ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[3],soil_info[4],soil_info[4]),dim=c(4)))
                                     ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[1],soil_info[2],soil_info[2]),dim=c(4)))
                                     ,exepath=as.character(crop_file_location),pathlength=as.integer(nchar(crop_file_location)))
      output = tmp$out_var ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP = tmp$out_var2  ; aNPP = array(aNPP, dim=c(nos_iter,aNPP_dim))
      MTT = tmp$out_var3   ; MTT = array(MTT, dim=c(nos_iter,MTT_dim))
      SS = tmp$out_var4    ; SS = array(SS, dim=c(nos_iter,SS_dim))
      aMTT = tmp$out_var5  ; aMTT = array(aMTT, dim=c(nos_iter,MTT_dim,noyears))
      MTTnat = tmp$out_var6; MTTnat = array(MTTnat, dim=c(nos_iter,MTT_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], lit_gCm2 = output[,,10],
                      lab_gCm2 = output[,,11], fol_gCm2 = output[,,12],
                      harvest_C_gCm2day = output[,,13], gsi = output[,,14],
                      gsi_itemp = output[,,15], gsi_iphoto = output[,,16],
                      gsi_ivpd = output[,,17], evap_kgH2Om2day = output[,,18],
                      sfc_water_mm = output[,,19], wSWP_MPa = output[,,20],
                      litwood_gCm2 = output[,,21], fire_gCm2day = output[,,23],
                      gs_demand_supply = output[,,24], gs_total_canopy = output[,,25],
                      APAR_MJm2day = output[,,26], gb_total_canopy = output[,,27],
                      CiCa = output[,,28],
                      aNPP = aNPP, MTT = MTT, SS = SS, aMTT = aMTT, natMTT = MTTnat)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  } else if (model_name == "DALEC_G6") {
      output_dim=28
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      crop_file_location=paste(PROJECT$exepath,"winter_wheat_development.csv", sep="")
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdalecg6",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim)
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
      output = tmp$out_var ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP = tmp$out_var2  ; aNPP = array(aNPP, dim=c(nos_iter,aNPP_dim))
      MTT = tmp$out_var3   ; MTT = array(MTT, dim=c(nos_iter,MTT_dim))
      SS = tmp$out_var4    ; SS = array(SS, dim=c(nos_iter,SS_dim))
      aMTT = tmp$out_var5  ; aMTT = array(aMTT, dim=c(nos_iter,MTT_dim,noyears))
      MTTnat = tmp$out_var6; MTTnat = array(MTTnat, dim=c(nos_iter,MTT_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], lit_gCm2 = output[,,10],
                      lab_gCm2 = output[,,11], fol_gCm2 = output[,,12],
                      harvest_C_gCm2day = output[,,13], cgi = output[,,14],
                      ncce_gCm2day = output[,,15], rSWP_MPa = output[,,16],
                      potH2O_supply_mmolH2Om2s = output[,,17], evap_kgH2Om2day = output[,,18],
                      sfc_water_mm = output[,,19], wSWP_MPa = output[,,20],
                      litwood_gCm2 = output[,,21], fire_gCm2day = output[,,23],
                      gs_demand_supply = output[,,24], gs_total_canopy = output[,,25],
                      APAR_MJm2day = output[,,26], gb_total_canopy = output[,,27],
                      CiCa = output[,,28],
                      aNPP = aNPP, MTT = MTT, SS = SS, aMTT = aMTT, natMTT = MTTnat)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  } else if (model_name == "DALEC_BUCKET_CanAGE") {
      output_dim=28
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      crop_file_location=paste(PROJECT$exepath,"winter_wheat_development.csv", sep="")
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdalecbucketcanage",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim)
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
                                     ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[3],soil_info[4],soil_info[4]),dim=c(4)))
                                     ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[1],soil_info[2],soil_info[2]),dim=c(4)))
                                     ,exepath=as.character(crop_file_location),pathlength=as.integer(nchar(crop_file_location)))
      output = tmp$out_var ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP = tmp$out_var2  ; aNPP = array(aNPP, dim=c(nos_iter,aNPP_dim))
      MTT = tmp$out_var3   ; MTT = array(MTT, dim=c(nos_iter,MTT_dim))
      SS = tmp$out_var4    ; SS = array(SS, dim=c(nos_iter,SS_dim))
      aMTT = tmp$out_var5  ; aMTT = array(aMTT, dim=c(nos_iter,MTT_dim,noyears))
      MTTnat = tmp$out_var6; MTTnat = array(MTTnat, dim=c(nos_iter,MTT_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], lit_gCm2 = output[,,10],
                      lab_gCm2 = output[,,11], fol_gCm2 = output[,,12],
                      harvest_C_gCm2day = output[,,13], gsi = output[,,14],
                      gsi_itemp = output[,,15], gsi_iphoto = output[,,16],
                      gsi_ivpd = output[,,17], evap_kgH2Om2day = output[,,18],
                      sfc_water_mm = output[,,19], wSWP_MPa = output[,,20],
                      litwood_gCm2 = output[,,21], fire_gCm2day = output[,,23],
                      gs_demand_supply = output[,,24], gs_total_canopy = output[,,25],
                      APAR_MJm2day = output[,,26], gb_total_canopy = output[,,27],
                      CiCa = output[,,28],
                      aNPP = aNPP, MTT = MTT, SS = SS, aMTT = aMTT, natMTT = MTTnat)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  } else if (model_name == "DALEC") {
      output_dim=25
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      crop_file_location=paste(PROJECT$exepath,"winter_wheat_development.csv", sep="")
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdalec",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim)
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
                            ,exepath=as.character(crop_file_location),pathlength=as.integer(nchar(crop_file_location)))
      output = tmp$out_var ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP = tmp$out_var2  ; aNPP = array(aNPP, dim=c(nos_iter,aNPP_dim))
      MTT = tmp$out_var3   ; MTT = array(MTT, dim=c(nos_iter,MTT_dim))
      SS = tmp$out_var4    ; SS = array(SS, dim=c(nos_iter,SS_dim))
      aMTT = tmp$out_var5  ; aMTT = array(aMTT, dim=c(nos_iter,MTT_dim,noyears))
      MTTnat = tmp$out_var6; MTTnat = array(MTTnat, dim=c(nos_iter,MTT_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], lit_gCm2 = output[,,10],
                      lab_gCm2 = output[,,11], fol_gCm2 = output[,,12],
                      harvest_C_gCm2day = output[,,13], gsi = output[,,14],
                      gsi_itemp = output[,,15], gsi_iphoto = output[,,16],
                      gsi_ivpd = output[,,17], gs_demand_supply = output[,,18],
                      gs_total_canopy = output[,,19], gb_total_canopy = output[,,20],
                      litwood_gCm2 = output[,,21], fire_gCm2day = output[,,23],
                      APAR_MJm2day = output[,,24], CiCa = output[,,25],
                      aNPP = aNPP, MTT = MTT, SS = SS, aMTT = aMTT, natMTT = MTTnat)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  } else if (model_name == "DALEC_GSI_BUCKET") {
      output_dim=28
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      crop_file_location=paste(PROJECT$exepath,"winter_wheat_development.csv", sep="")
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdalecgsibucket",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim)
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
                                     ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[3],soil_info[4],soil_info[4]),dim=c(4)))
                                     ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[1],soil_info[2],soil_info[2]),dim=c(4)))
                                     ,exepath=as.character(crop_file_location),pathlength=as.integer(nchar(crop_file_location)))
      output = tmp$out_var ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP = tmp$out_var2  ; aNPP = array(aNPP, dim=c(nos_iter,aNPP_dim))
      MTT = tmp$out_var3   ; MTT = array(MTT, dim=c(nos_iter,MTT_dim))
      SS = tmp$out_var4    ; SS = array(SS, dim=c(nos_iter,SS_dim))
      aMTT = tmp$out_var5  ; aMTT = array(aMTT, dim=c(nos_iter,MTT_dim,noyears))
      MTTnat = tmp$out_var6; MTTnat = array(MTTnat, dim=c(nos_iter,MTT_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], lit_gCm2 = output[,,10],
                      lab_gCm2 = output[,,11], fol_gCm2 = output[,,12],
                      harvest_C_gCm2day = output[,,13], gsi = output[,,14],
                      gsi_itemp = output[,,15], gsi_iphoto = output[,,16],
                      gsi_ivpd = output[,,17], evap_kgH2Om2day = output[,,18],
                      sfc_water_mm = output[,,19], wSWP_MPa = output[,,20],
                      litwood_gCm2 = output[,,21], fire_gCm2day = output[,,23],
                      gs_demand_supply = output[,,24], gs_total_canopy = output[,,25],
                      APAR_MJm2day = output[,,26], gb_total_canopy = output[,,27],
                      CiCa = output[,,28],
                      aNPP = aNPP, MTT = MTT, SS = SS, aMTT = aMTT, natMTT = MTTnat)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  } else if (model_name == "DALECN_GSI_BUCKET") {
      output_dim=23
      if (is.loaded("rdalecngsibucket") == FALSE) { dyn.load(paste(PROJECT$exepath,"/dalec.so", sep="")) }
      crop_file_location=paste(PROJECT$exepath,"winter_wheat_development.csv", sep="")
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran("rdalecngsibucket",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
                                     ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                                     ,out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                     ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                     ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                     ,pft=as.integer(pft),nodays=as.integer(dim(met)[1]),deltat=as.double(array(0,dim=c(as.integer(dim(met)[1]))))
                                     ,nos_iter=as.integer(nos_iter)
                                     ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[3],soil_info[4],soil_info[4]),dim=c(4)))
                                     ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[1],soil_info[2],soil_info[2]),dim=c(4)))
                                     ,exepath=as.character(crop_file_location),pathlength=as.integer(nchar(crop_file_location)))
      output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
      if (site == PROJECT$sites[length(PROJECT$sites)]) {dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))}
      rm(tmp) ; gc()
      stop('Code to assign variables to output has not been re-written to current code standard - ooops')
  } else if (model_name == "DALECN_BUCKET") {
      output_dim=23
      if (is.loaded("rdalecnbucket") == FALSE) { dyn.load(paste(PROJECT$exepath,"/dalec.so", sep="")) }
      crop_file_location=as.character(paste(PROJECT$exepath,"winter_wheat_development.csv", sep=""))
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran("rdalecnbucket",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
                                  ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                  ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                  ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                  ,pft=as.integer(pft),nodays=as.integer(dim(met)[1]),deltat=as.double(array(0,dim=c(as.integer(dim(met)[1]))))
                                  ,nos_iter=as.integer(nos_iter)
                                  ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[3],soil_info[4],soil_info[4]),dim=c(4)))
                                  ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[1],soil_info[2],soil_info[2]),dim=c(4)))
                                  ,exepath=as.character(crop_file_location),pathlength=as.integer(nchar(crop_file_location)))
      output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
      if (site == PROJECT$sites[length(PROJECT$sites)]) {dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))}
      rm(tmp) ; gc()
      # create output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], lit_gCm2 = output[,,10],
                      lab_gCm2 = output[,,11], fol_gCm2 = output[,,12],
                      harvest_C_gCm2day = output[,,13], canopyage_days = output[,,14],
                      evap_kgH2Om2day = output[,,18],
                      sfc_water_mm = output[,,19], wSWP_MPa = output[,,20],
                      litwood_gCm2 = output[,,21], fire_gCm2day = output[,,23], aNPP = aNPP)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  } else if (model_name == "DALEC_CDEA") {
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      pft_specific = 0
      tmp=.Fortran( "rdaleccdea",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
                                ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
      output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
  } else if (model_name == "DALEC_CDEA_FR") {
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdaleccdeafr",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
                                  ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                  ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                  ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                  ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                  ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
      output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      stop('Code to assign variables to output has not been re-written to current code standard - ooops')
  } else if (model_name == "DALEC_CDEA_LU_FIRES") {
      output_dim=15
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdaleccdealufires",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim)
                                       ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                                       ,met=as.double(t(met))
                                       ,pars=as.double(pars_in),out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                                       ,out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                       ,out_var3=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                       ,out_var4=as.double(array(0,dim=c(nos_iter,SS_dim)))
                                       ,out_var5=as.double(array(0,dim=c(nos_iter,MTT_dim,noyears)))
                                       ,out_var6=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                       ,lat=as.double(lat)
                                       ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                       ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                       ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                       ,noyears=as.integer(noyears)
                                       ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
      output = tmp$out_var ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP = tmp$out_var2  ; aNPP = array(aNPP, dim=c(nos_iter,aNPP_dim))
      MTT = tmp$out_var3   ; MTT = array(MTT, dim=c(nos_iter,MTT_dim))
      SS = tmp$out_var4    ; SS = array(SS, dim=c(nos_iter,SS_dim))
      aMTT = tmp$out_var5  ; aMTT = array(aMTT, dim=c(nos_iter,MTT_dim,noyears))
      MTTnat = tmp$out_var6; MTTnat = array(MTTnat, dim=c(nos_iter,MTT_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # Create the output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], lit_gCm2 = output[,,10],
                      lab_gCm2 = output[,,11], fol_gCm2 = output[,,12],
                      harvest_C_gCm2day = output[,,13], fire_gCm2day = output[,,14],
                      CiCa = output[,,15],
                      aNPP = aNPP, MTT = MTT, SS = SS, aMTT = aMTT, natMTT = MTTnat)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  } else if (model_name == "DALEC_CDEA_ACM2") {
      output_dim=19
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdaleccdeaacm2",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim)
                                    ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                                    ,met=as.double(t(met))
                                    ,pars=as.double(pars_in),out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                                    ,out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                    ,out_var3=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                    ,out_var4=as.double(array(0,dim=c(nos_iter,SS_dim)))
                                    ,out_var5=as.double(array(0,dim=c(nos_iter,MTT_dim,noyears)))
                                    ,out_var6=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                    ,lat=as.double(lat)
                                    ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                    ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                    ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                    ,noyears=as.integer(noyears)
                                    ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
      output = tmp$out_var ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP = tmp$out_var2  ; aNPP = array(aNPP, dim=c(nos_iter,aNPP_dim))
      MTT = tmp$out_var3   ; MTT = array(MTT, dim=c(nos_iter,MTT_dim))
      SS = tmp$out_var4    ; SS = array(SS, dim=c(nos_iter,SS_dim))
      aMTT = tmp$out_var5  ; aMTT = array(aMTT, dim=c(nos_iter,MTT_dim,noyears))
      MTTnat = tmp$out_var6; MTTnat = array(MTTnat, dim=c(nos_iter,MTT_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # Create the output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], lit_gCm2 = output[,,10],
                      lab_gCm2 = output[,,11], fol_gCm2 = output[,,12],
                      harvest_C_gCm2day = output[,,13], fire_gCm2day = output[,,14],
                      gs_demand_supply = output[,,15], gs_total_canopy = output[,,16],
                      APAR_MJm2day = output[,,17], gb_total_canopy = output[,,18],
                      CiCa = output[,,19],
                      aNPP = aNPP, MTT = MTT, SS = SS, aMTT = aMTT, natMTT = MTTnat)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  } else if (model_name == "DALEC_CDEA_ACM2_BUCKET") {
      output_dim=25
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdaleccdeaacm2bucket",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim)
                                          ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                                          ,met=as.double(t(met))
                                          ,pars=as.double(pars_in),out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                                          ,out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                          ,out_var3=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                          ,out_var4=as.double(array(0,dim=c(nos_iter,SS_dim)))
                                          ,out_var5=as.double(array(0,dim=c(nos_iter,MTT_dim,noyears)))
                                          ,out_var6=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                          ,lat=as.double(lat)
                                          ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                          ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                          ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                          ,noyears=as.integer(noyears)
                                          ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                                          ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[3],soil_info[4],soil_info[4]),dim=c(4)))
                                          ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[1],soil_info[2],soil_info[2]),dim=c(4))))
      output = tmp$out_var ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP = tmp$out_var2  ; aNPP = array(aNPP, dim=c(nos_iter,aNPP_dim))
      MTT = tmp$out_var3   ; MTT = array(MTT, dim=c(nos_iter,MTT_dim))
      SS = tmp$out_var4    ; SS = array(SS, dim=c(nos_iter,SS_dim))
      aMTT = tmp$out_var5  ; aMTT = array(aMTT, dim=c(nos_iter,MTT_dim,noyears))
      MTTnat = tmp$out_var6; MTTnat = array(MTTnat, dim=c(nos_iter,MTT_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # Create the output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], lit_gCm2 = output[,,10],
                      lab_gCm2 = output[,,11], fol_gCm2 = output[,,12],
                      harvest_C_gCm2day = output[,,13], fire_gCm2day = output[,,14],
                      evap_kgH2Om2day = output[,,18],sfc_water_mm = output[,,19],
                      wSWP_MPa = output[,,20],gs_demand_supply = output[,,21],
                      gs_total_canopy = output[,,22],APAR_MJm2day = output[,,23],
                      gb_total_canopy = output[,,24],CiCa = output[,,25],
                      aNPP = aNPP, MTT = MTT, SS = SS, aMTT = aMTT, natMTT = MTTnat)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  } else if (model_name == "DALEC_CDEA_ACM2_BUCKET_RmRg") {
      output_dim=25
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdaleccdeaacm2bucketRmRg",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim)
                                          ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                                          ,met=as.double(t(met))
                                          ,pars=as.double(pars_in),out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                                          ,out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                          ,out_var3=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                          ,out_var4=as.double(array(0,dim=c(nos_iter,SS_dim)))
                                          ,out_var5=as.double(array(0,dim=c(nos_iter,MTT_dim,noyears)))
                                          ,out_var6=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                          ,lat=as.double(lat)
                                          ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                          ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                          ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                          ,noyears=as.integer(noyears)
                                          ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                                          ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[3],soil_info[4],soil_info[4]),dim=c(4)))
                                          ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[1],soil_info[2],soil_info[2]),dim=c(4))))
      output = tmp$out_var ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP = tmp$out_var2  ; aNPP = array(aNPP, dim=c(nos_iter,aNPP_dim))
      MTT = tmp$out_var3   ; MTT = array(MTT, dim=c(nos_iter,MTT_dim))
      SS = tmp$out_var4    ; SS = array(SS, dim=c(nos_iter,SS_dim))
      aMTT = tmp$out_var5  ; aMTT = array(aMTT, dim=c(nos_iter,MTT_dim,noyears))
      MTTnat = tmp$out_var6; MTTnat = array(MTTnat, dim=c(nos_iter,MTT_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # Create the output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], lit_gCm2 = output[,,10],
                      lab_gCm2 = output[,,11], fol_gCm2 = output[,,12],
                      harvest_C_gCm2day = output[,,13], fire_gCm2day = output[,,14],
                      evap_kgH2Om2day = output[,,18],sfc_water_mm = output[,,19],
                      wSWP_MPa = output[,,20],gs_demand_supply = output[,,21],
                      gs_total_canopy = output[,,22],APAR_MJm2day = output[,,23],
                      gb_total_canopy = output[,,24],CiCa = output[,,25],
                      aNPP = aNPP, MTT = MTT, SS = SS, aMTT = aMTT, natMTT = MTTnat)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  } else if (model_name == "DALEC_CDEA_ACM2_BUCKET_RmRg_CWD") {
      output_dim=26
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdaleccdeaacm2bucketRmRgcwd",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim)
                                          ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                                          ,met=as.double(t(met))
                                          ,pars=as.double(pars_in),out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                                          ,out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                          ,out_var3=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                          ,out_var4=as.double(array(0,dim=c(nos_iter,SS_dim)))
                                          ,out_var5=as.double(array(0,dim=c(nos_iter,MTT_dim,noyears)))
                                          ,out_var6=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                          ,lat=as.double(lat)
                                          ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                          ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                          ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                          ,noyears=as.integer(noyears)
                                          ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                                          ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[3],soil_info[4],soil_info[4]),dim=c(4)))
                                          ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[1],soil_info[2],soil_info[2]),dim=c(4))))
      output = tmp$out_var ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP = tmp$out_var2  ; aNPP = array(aNPP, dim=c(nos_iter,aNPP_dim))
      MTT = tmp$out_var3   ; MTT = array(MTT, dim=c(nos_iter,MTT_dim))
      SS = tmp$out_var4    ; SS = array(SS, dim=c(nos_iter,SS_dim))
      aMTT = tmp$out_var5  ; aMTT = array(aMTT, dim=c(nos_iter,MTT_dim,noyears))
      MTTnat = tmp$out_var6; MTTnat = array(MTTnat, dim=c(nos_iter,MTT_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # Create the output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], lit_gCm2 = output[,,10],
                      lab_gCm2 = output[,,11], fol_gCm2 = output[,,12],
                      harvest_C_gCm2day = output[,,13], fire_gCm2day = output[,,14],
                      evap_kgH2Om2day = output[,,18],sfc_water_mm = output[,,19],
                      wSWP_MPa = output[,,20],gs_demand_supply = output[,,21],
                      gs_total_canopy = output[,,22],APAR_MJm2day = output[,,23],
                      gb_total_canopy = output[,,24],litwood_gCm2 = output[,,25],
                      CiCa = output[,,26],
                      aNPP = aNPP, MTT = MTT, SS = SS, aMTT = aMTT, natMTT = MTTnat)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  } else if (model_name == "DALEC_CDEA_ACM2_BUCKET_RmRg_CWD_wMRT") {
      output_dim=26
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdaleccdeaacm2bucketRmRgcwdwmrt",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim)
                                          ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                                          ,met=as.double(t(met))
                                          ,pars=as.double(pars_in),out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                                          ,out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                          ,out_var3=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                          ,out_var4=as.double(array(0,dim=c(nos_iter,SS_dim)))
                                          ,out_var5=as.double(array(0,dim=c(nos_iter,MTT_dim,noyears)))
                                          ,out_var6=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                          ,lat=as.double(lat)
                                          ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                          ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                          ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                          ,noyears=as.integer(noyears)
                                          ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                                          ,soil_frac_clay_in=as.double(array(c(soil_info[3],soil_info[3],soil_info[4],soil_info[4]),dim=c(4)))
                                          ,soil_frac_sand_in=as.double(array(c(soil_info[1],soil_info[1],soil_info[2],soil_info[2]),dim=c(4))))
      output = tmp$out_var ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP = tmp$out_var2  ; aNPP = array(aNPP, dim=c(nos_iter,aNPP_dim))
      MTT = tmp$out_var3   ; MTT = array(MTT, dim=c(nos_iter,MTT_dim))
      SS = tmp$out_var4    ; SS = array(SS, dim=c(nos_iter,SS_dim))
      aMTT = tmp$out_var5  ; aMTT = array(aMTT, dim=c(nos_iter,MTT_dim,noyears))
      MTTnat = tmp$out_var6; MTTnat = array(MTTnat, dim=c(nos_iter,MTT_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # Create the output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], lit_gCm2 = output[,,10],
                      lab_gCm2 = output[,,11], fol_gCm2 = output[,,12],
                      harvest_C_gCm2day = output[,,13], fire_gCm2day = output[,,14],
                      evap_kgH2Om2day = output[,,18],sfc_water_mm = output[,,19],
                      wSWP_MPa = output[,,20],gs_demand_supply = output[,,21],
                      gs_total_canopy = output[,,22],APAR_MJm2day = output[,,23],
                      gb_total_canopy = output[,,24],litwood_gCm2 = output[,,25],
                      CiCa = output[,,26],
                      aNPP = aNPP, MTT = MTT, SS = SS, aMTT = aMTT, natMTT = MTTnat)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  } else if (model_name == "DALEC_CDEA_no_lit_root") {
      output_dim=19
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdaleccdeanolitroot",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim)
                                         ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                                         ,met=as.double(t(met))
                                         ,pars=as.double(pars_in),out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                                         ,out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                         ,out_var3=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                         ,out_var4=as.double(array(0,dim=c(nos_iter,SS_dim)))
                                         ,out_var5=as.double(array(0,dim=c(nos_iter,MTT_dim,noyears)))
                                         ,out_var6=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                         ,lat=as.double(lat)
                                         ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                         ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                         ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                         ,noyears=as.integer(noyears)
                                         ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
      output = tmp$out_var ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP = tmp$out_var2  ; aNPP = array(aNPP, dim=c(nos_iter,aNPP_dim))
      MTT = tmp$out_var3   ; MTT = array(MTT, dim=c(nos_iter,MTT_dim))
      SS = tmp$out_var4    ; SS = array(SS, dim=c(nos_iter,SS_dim))
      aMTT = tmp$out_var5  ; aMTT = array(aMTT, dim=c(nos_iter,MTT_dim,noyears))
      MTTnat = tmp$out_var6; MTTnat = array(MTTnat, dim=c(nos_iter,MTT_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # Create the output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], lit_gCm2 = output[,,10],
                      lab_gCm2 = output[,,11], fol_gCm2 = output[,,12],
                      harvest_C_gCm2day = output[,,13], fire_gCm2day = output[,,14],
                      aNPP = aNPP, MTT = MTT, SS = SS, aMTT = aMTT, natMTT = MTTnat)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  } else if (model_name == "DALEC_EVERGREEN") {
      output_dim=19
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdalecevergreen",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim)
                                     ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                                     ,met=as.double(t(met))
                                     ,pars=as.double(pars_in),out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                                     ,out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                     ,out_var3=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                     ,out_var4=as.double(array(0,dim=c(nos_iter,SS_dim)))
                                     ,out_var5=as.double(array(0,dim=c(nos_iter,MTT_dim,noyears)))
                                     ,out_var6=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                     ,lat=as.double(lat)
                                     ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                     ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                     ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                     ,noyears=as.integer(noyears)
                                     ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
      output = tmp$out_var ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP = tmp$out_var2  ; aNPP = array(aNPP, dim=c(nos_iter,aNPP_dim))
      MTT = tmp$out_var3   ; MTT = array(MTT, dim=c(nos_iter,MTT_dim))
      SS = tmp$out_var4    ; SS = array(SS, dim=c(nos_iter,SS_dim))
      aMTT = tmp$out_var5  ; aMTT = array(aMTT, dim=c(nos_iter,MTT_dim,noyears))
      MTTnat = tmp$out_var6; MTTnat = array(MTTnat, dim=c(nos_iter,MTT_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # Create the output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], lit_gCm2 = output[,,10],
                      lab_gCm2 = output[,,11], fol_gCm2 = output[,,12],
                      harvest_C_gCm2day = output[,,13], fire_gCm2day = output[,,14],
                      aNPP = aNPP, MTT = MTT, SS = SS, aMTT = aMTT, natMTT = MTTnat)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  } else if (model_name == "DALEC_EVERGREEN_no_lit_root") {
      output_dim=19
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdalecevergreennolitroot",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim)
                                     ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                                     ,met=as.double(t(met))
                                     ,pars=as.double(pars_in),out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                                     ,out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                     ,out_var3=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                     ,out_var4=as.double(array(0,dim=c(nos_iter,SS_dim)))
                                     ,out_var5=as.double(array(0,dim=c(nos_iter,MTT_dim,noyears)))
                                     ,out_var6=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                     ,lat=as.double(lat)
                                     ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                     ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                     ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                     ,noyears=as.integer(noyears)
                                     ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
      output = tmp$out_var ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP = tmp$out_var2  ; aNPP = array(aNPP, dim=c(nos_iter,aNPP_dim))
      MTT = tmp$out_var3   ; MTT = array(MTT, dim=c(nos_iter,MTT_dim))
      SS = tmp$out_var4    ; SS = array(SS, dim=c(nos_iter,SS_dim))
      aMTT = tmp$out_var5  ; aMTT = array(aMTT, dim=c(nos_iter,MTT_dim,noyears))
      MTTnat = tmp$out_var6; MTTnat = array(MTTnat, dim=c(nos_iter,MTT_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # Create the output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], lit_gCm2 = output[,,10],
                      lab_gCm2 = output[,,11], fol_gCm2 = output[,,12],
                      harvest_C_gCm2day = output[,,13], fire_gCm2day = output[,,14],
                      aNPP = aNPP, MTT = MTT, SS = SS, aMTT = aMTT, natMTT = MTTnat)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  } else if (model_name == "DALEC_GSI_FR_LABILE") {
      output_dim=18
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdalecgsifr",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
                                 ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                 ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                 ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                 ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                 ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
      output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], lit_gCm2 = output[,,10],
                      lab_gCm2 = output[,,11], fol_gCm2 = output[,,12],
                      harvest_C_gCm2day = output[,,13], gsi = output[,,14],
                      gsi_itemp = output[,,15], gsi_iphoto = output[,,16],
                      gsi_ivpd = output[,,17], labile_slow_gCm2 = output[,,18], aNPP = aNPP)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  } else if (model_name == "DALECN_GSI_FR") {
      output_dim=22
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdalecngsifr",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
                                  ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                  ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                  ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                  ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                  ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
      output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], lit_gCm2 = output[,,10],
                      lab_gCm2 = output[,,11], fol_gCm2 = output[,,12],
                      harvest_C_gCm2day = output[,,13], gsi = output[,,14],
                      gsi_itemp = output[,,15], gsi_iphoto = output[,,16],
                      gsi_ivpd = output[,,17], litwood_gCm2 = output[,,18],
                      litN_gNm2 = output[,,19], labN_gNm2 = output[,,20],
                      DIN_gNm2 = output[,,21], N_mineralisation_gNm2 = output[,,22], aNPP = aNPP)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  } else if (model_name == "DALEC_GSI_FR") {
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdalecgsifr",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
                                 ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                 ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                 ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                 ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                 ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
    output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
    aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
    dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
    rm(tmp) ; gc()
    stop('Code to assign variables to output has not been re-written to current code standard - ooops')
  } else if (model_name == "DALEC_GSI_DFOL_FR") {
    dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
    if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
    tmp=.Fortran( "rdalecgsidfolfr",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
                                   ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                   ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                   ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                   ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                   ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
      output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      stop('Code to assign variables to output has not been re-written to current code standard - ooops')
  } else if (model_name == "DALEC_GSI_DFOL_CWD_FR") {
      output_dim=24
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      crop_file_location=paste(PROJECT$exepath,"winter_wheat_development.csv", sep="")
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdalecgsidfolcwdfr",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim)
                                        ,MTT_dim=as.integer(MTT_dim),SS_dim = as.integer(SS_dim)
                                        ,met=as.double(t(met))
                                        ,pars=as.double(pars_in),out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                                        ,out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                        ,out_var3=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                        ,out_var4=as.double(array(0,dim=c(nos_iter,SS_dim)))
                                        ,out_var5=as.double(array(0,dim=c(nos_iter,MTT_dim,noyears)))
                                        ,out_var6=as.double(array(0,dim=c(nos_iter,MTT_dim)))
                                        ,lat=as.double(lat)
                                        ,nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                        ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                        ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                        ,noyears=as.integer(noyears)
                                        ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                                        ,exepath=as.character(crop_file_location),pathlength=as.integer(nchar(crop_file_location)))
      output = tmp$out_var ; output = array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP = tmp$out_var2  ; aNPP = array(aNPP, dim=c(nos_iter,aNPP_dim))
      MTT = tmp$out_var3   ; MTT = array(MTT, dim=c(nos_iter,MTT_dim))
      SS = tmp$out_var4    ; SS = array(SS, dim=c(nos_iter,SS_dim))
      aMTT = tmp$out_var5  ; aMTT = array(aMTT, dim=c(nos_iter,MTT_dim,noyears))
      MTTnat = tmp$out_var6; MTTnat = array(MTTnat, dim=c(nos_iter,MTT_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], lit_gCm2 = output[,,10],
                      lab_gCm2 = output[,,11], fol_gCm2 = output[,,12],
                      harvest_C_gCm2day = output[,,13], gsi = output[,,14],
                      gsi_itemp = output[,,15], gsi_iphoto = output[,,16],
                      gsi_ivpd = output[,,17], litwood_gCm2 = output[,,18],
                      fire_gCm2day = output[,,19], gs_demand_supply = output[,,20],
                      gs_total_canopy = output[,,21], APAR_MJm2day = output[,,22],
                      gb_total_canopy = output[,,23], CiCa = output[,,24],
                      aNPP = aNPP, MTT = MTT, SS = SS, aMTT = aMTT, natMTT = MTTnat)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  } else if (model_name == "DALEC_GSI_DFOL_FROOT_FR") {
      output_dim=19
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdalecgsidfolfrootfr",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
                                          ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                                          ,out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                          ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                          ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                          ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                          ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
      output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], lit_gCm2 = output[,,10],
                      lab_gCm2 = output[,,11], fol_gCm2 = output[,,12],
                      harvest_C_gCm2day = output[,,13], gsi = output[,,14],
                      gsi_itemp = output[,,15], gsi_iphoto = output[,,16],
                      gsi_ivpd = output[,,17], labroot_gCm2 = output[,,18],
                      labwood_gCm2 = output[,,19], aNPP = aNPP)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  } else if (model_name == "DALECN_GSI_DFOL_LABILE_FR") {
      output_dim=24
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdalecngsidfollabfr",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
                                         ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                         ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                         ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                         ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                         ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
      output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], lit_gCm2 = output[,,10],
                      lab_gCm2 = output[,,11], fol_gCm2 = output[,,12],
                      harvest_C_gCm2day = output[,,13], gsi = output[,,14],
                      gsi_itemp = output[,,15], gsi_iphoto = output[,,16],
                      gsi_ivpd = output[,,17], labroot_gCm2 = output[,,18],
                      labwood_gCm2 = output[,,19], litwood_gCm2 = output[,,20],
                      litN_gNm2 = output[,,21], labN_gNm2 = output[,,22],
                      DIN_gNm2 = output[,,23], N_mineralisation_gNm2day = output[,,24], aNPP = aNPP)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  } else if (model_name == "DALECN_GSI_DFOL_LABILE_FROOT_FR") {
      output_dim=24
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdalecngsidfollabfrootfr",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
                                              ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                              ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                              ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                              ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                              ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
      output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], lit_gCm2 = output[,,10],
                      lab_gCm2 = output[,,11], fol_gCm2 = output[,,12],
                      harvest_C_gCm2day = output[,,13], gsi = output[,,14],
                      gsi_itemp = output[,,15], gsi_iphoto = output[,,16],
                      gsi_ivpd = output[,,17], labroot_gCm2 = output[,,18],
                      labwood_gCm2 = output[,,19], litwood_gCm2 = output[,,20],
                      litN_gNm2 = output[,,21], labN_gNm2 = output[,,22],
                      DIN_gNm2 = output[,,23], N_mineralisation_gNm2day = output[,,24], aNPP = aNPP)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  } else if (model_name == "DALEC_GSI_DFOL_LABILE_FR") {
      output_dim=20
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdalecgsidfollabfr",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
                                        ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                        ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                        ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                        ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                        ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
      output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], lit_gCm2 = output[,,10],
                      lab_gCm2 = output[,,11], fol_gCm2 = output[,,12],
                      harvest_C_gCm2day = output[,,13], gsi = output[,,14],
                      gsi_itemp = output[,,15], gsi_iphoto = output[,,16],
                      gsi_ivpd = output[,,17], labroot_gCm2 = output[,,18],
                      labwood_gCm2 = output[,,19], litwood_gCm2 = output[,,20], aNPP = aNPP)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  } else if (model_name == "DALEC_GSI_MFOL_FR") {
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdalecgsimfolfr",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
                                     ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                     ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                     ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                     ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                     ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
      output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      stop('Code to assign variables to output has not been re-written to current code standard - ooops')
  } else if (model_name == "DALEC_GSI_DBio_FR") {
      output_dim=22
      dyn.load(paste(PROJECT$exepath,"/dalec.so", sep=""))
      if (parameter_type == "pft_specific") {pft_specific = 1} else {pft_specific = 0}
      tmp=.Fortran( "rdalecgsibiofr",output_dim=as.integer(output_dim),aNPP_dim=as.integer(aNPP_dim),met=as.double(t(met)),pars=as.double(pars_in)
                                    ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim))),out_var2=as.double(array(0,dim=c(nos_iter,aNPP_dim)))
                                    ,lat=as.double(lat),nopars=as.integer(PROJECT$model$nopars[site]),nomet=as.integer(dim(met)[2])
                                    ,nofluxes=as.integer(PROJECT$model$nofluxes[site]),nopools=as.integer(PROJECT$model$nopools[site])
                                    ,pft=as.integer(pft),pft_specific=as.integer(pft_specific),nodays=as.integer(dim(met)[1])
                                    ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter))
      output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
      aNPP=tmp$out_var2  ; aNPP=array(aNPP, dim=c(nos_iter,aNPP_dim))
      dyn.unload(paste(PROJECT$exepath,"/dalec.so", sep=""))
      rm(tmp) ; gc()
      # create output object
      states_all=list(lai_m2m2 = output[,,1], gpp_gCm2day = output[,,2],
                      rauto_gCm2day = output[,,3], rhet_gCm2day = output[,,4],
                      nee_gCm2day = output[,,5], wood_gCm2 = output[,,6],
                      som_gCm2 = output[,,7], bio_gCm2 = output[,,8],
                      root_gCm2 = output[,,9], lit_gCm2 = output[,,10],
                      lab_gCm2 = output[,,11], fol_gCm2 = output[,,12],
                      harvest_C_gCm2day = output[,,13], gsi = output[,,14],
                      gsi_itemp = output[,,15], gsi_iphoto = output[,,16],
                      gsi_ivpd = output[,,17], somfast_gCm2 = output[,,18],
                      litroot_gCm2 = output[,,19], litwood_gCm2 = output[,,20],
                      microact = output[,,21], microbial_gCm2 = output[,,22], aNPP = aNPP)
      # add newly calculated variables
      states_all$reco_gCm2day = states_all$rauto_gCm2day + states_all$rhet_gCm2day
  } else {
      stop(paste("Model choice (",model_name,") does not have corresponding R interface",sep=""))
  }

  # return state variable means
  return(states_all) ; gc(verbose=FALSE)

  } # end of function
  ## Use byte compile
  simulate_all<-cmpfun(simulate_all)
