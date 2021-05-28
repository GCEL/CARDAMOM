
###
## Function to run CARDAMOM parameters via the chosen model
###

run_each_site<-function(n,PROJECT,stage,repair,grid_override,stage5modifiers) {

  if (stage == 5) {
      stage5name = as.vector(unlist(stage5modifiers))
      stage5name = paste(stage5name[1],stage5name[2],stage5name[3],stage5name[4],stage5name[5],stage5name[6],stage5name[7],stage5name[8],sep="_")
      outfile = paste(PROJECT$results_processedpath,PROJECT$sites[n],"_",stage5name,".RData",sep="")
      outfile1 = paste(PROJECT$results_processedpath,PROJECT$sites[n],"_",stage5name,"_parameters.RData",sep="")
  } else {
      outfile = paste(PROJECT$results_processedpath,PROJECT$sites[n],".RData",sep="")
      outfile1 = paste(PROJECT$results_processedpath,PROJECT$sites[n],"_parameters.RData",sep="")
      outfile2 = paste(PROJECT$results_processedpath,PROJECT$sites[n],"_stock_fluxes.RData",sep="")
  }

  if (file.exists(outfile1) == FALSE | repair == 1) {
      # load only the desired latter fraction of the parameter vectors
      # output is order dimensions(npar+1,iter,chain)
      parameters = read_parameter_chains(PROJECT,n,3)
      parameter_covariance = read_parameter_covariance(PROJECT,n)

      # determine whether we have any actual completed chains and whether they include EDC consistent value only
      error_check = FALSE
      if (parameters[1] == -9999) {
          error_check = TRUE
      } else {
          if (length(which(as.vector(is.na(parameters)))) > 0 ) {
              error_check = TRUE
              print("NA found in likelihood score")
           } else if (min(as.vector(parameters)) == -Inf) {
              error_check = TRUE
              print("Inf found in likelihood score")
          } # NaN / Inf check
      } # error check

      # ok so if we actually have some parameters we will
      if (error_check == FALSE) {

          # test for convergence and whether or not there is any single chain which can be removed in they do not converge
          notconv = TRUE
          while (dim(parameters)[3] > 2 & notconv) {
              print("begin convergence checking")
              converged = have_chains_converged(parameters)
              # if log-likelihood has passed then we are not interested
              if (converged[length(converged)] == "FAIL") {
                  print("...not converged begin removing potential parameter vectors")
                  i = 1 ; max_likelihood = rep(NA, length.out=dim(parameters)[3]) ; CI90 = rep(NA,length.out=c(2))
                  while (notconv){
                     print(paste("......trying removal of parameter vector ",i,sep=""))
                     # Track the maximum likelihood across each chain.
                     max_likelihood[i] = max(parameters[dim(parameters)[1],,i])
                     converged = have_chains_converged(parameters[,,-i]) ; i = i + 1
                     # if removing one of the chains get convergence then great
                     if (converged[length(converged)] == "PASS") {
                         print(".........convergence found on chain removal")
                         # likelihoods converge now but we need to check for the possibility that the chain we have removed is actually better than the others
                         CI90[1] = quantile(parameters[dim(parameters)[1],,(i-1)], prob=c(0.10)) ; CI90[2] = quantile(parameters[dim(parameters)[1],,-(i-1)], prob=c(0.90))
                         # if the rejected chain is significantly better (at 90 % CI) than the converged chains then we have a problem
                         if (CI90[1] > CI90[2]) {
                             # rejected chain (while others converge) is actually better and the others have gotten stuck in a local minima.
                             # we will now assume that we use the single good chain instead...
                             parameters = array(parameters[,,(i-1)],dim=c(dim(parameters)[1:2],2))
                             notconv = FALSE ; i = (i-1) * -1
                             print(paste("............chain ",i*-1," only has been accepted",sep=""))
                         } else {
                             # if the non-converged chain is worse or just the same in likelihood terms as the others then we will ditch it
                             notconv = FALSE ; i = i-1 # converged now?
                             parameters = parameters[,,-i]
                             print(paste("............chain rejected = ",i,sep=""))
                         }
                     }
                     # If removing one chain does not lead to convergence then lowest average likelihood chain
                     if (i > dim(parameters)[3] & notconv) {
                         # Which is lowest likelihood
                         i = which(max_likelihood == min(max_likelihood))
                         # Remove from array
                         parameters = parameters[,,-i]
                         # Update the maximum likelihood vector also
                         max_likelihood = max_likelihood[-i]
                         # Update the user
                         print(paste(".........single chain removal couldn't find convergence; removing lowest likelihood chain = ",i,sep=""))
                         # reset counter
                         i = 1
                         # If we have removed chains down to 2 (or kept just one) we need to abort
                         if (dim(parameters)[3] < 3 & notconv) {
                             notconv = FALSE ; print(".........have removed all low likelihood chains without successful convergence")
                         }
                     }
                  } # while to removing chains
              } else {
                  print("All chains converge")
                  # we have conveged
                  notconv = FALSE
              } # if likelihood not converged
          } # if more than 2 chains

          # load the met data for each site
          drivers = read_binary_file_format(paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[n],".bin",sep=""))

          # if Stage 5 sensitivity analysis is enacted - apply factor scaled adjustment to the drivers
          if (stage == 5) {

              # do weather drivers first
              drivers$met[,2]  = drivers$met[,2]  * stage5modifiers$airt_factor # min temperature
              drivers$met[,3]  = drivers$met[,3]  * stage5modifiers$airt_factor # max temperature
              drivers$met[,4]  = drivers$met[,4]  * stage5modifiers$swrad_factor
              drivers$met[,5]  = drivers$met[,5]  * stage5modifiers$co2_factor
              drivers$met[,7]  = drivers$met[,7]  * stage5modifiers$rainfall_factor
              drivers$met[,15] = drivers$met[,15] * stage5modifiers$wind_spd_factor
              drivers$met[,16] = drivers$met[,16] * stage5modifiers$vpd_factor
              # do disturbance drivers next
              disturbed_locations = which(drivers$met[,8] > 0)
              drivers$met[disturbed_locations,8] = drivers$met[disturbed_locations,8] * stage5modifiers$deforestation_factor
              disturbed_locations = which(drivers$met[,9] > 0)
              drivers$met[disturbed_locations,9] = drivers$met[disturbed_locations,9] * stage5modifiers$burnt_area_factor
              # do GSI related parameters next
              drivers$met[,10] = drivers$met[,10] * stage5modifiers$airt_factor # 21 day rolling average daily minimum temperature

          }

          # run parameters for full results / propogation
          soil_info = c(drivers$top_sand,drivers$bot_sand,drivers$top_clay,drivers$bot_clay)
          print("running model ensemble")
          states_all = simulate_all(n,PROJECT,PROJECT$model$name,drivers$met,parameters[1:PROJECT$model$nopars[n],,],
                                    drivers$lat,PROJECT$ctessel_pft[n],PROJECT$parameter_type,
                                    PROJECT$exepath,soil_info)

          # Avoid running with ACM basically where not all fluxes exist
          if (grepl("DALEC",PROJECT$model$name)) {
              # Post-hoc calculation of parameter correlations with key C-cycle variables
              tmp = t(array(as.vector(parameters[1:PROJECT$model$nopars[n],,]),dim=c(PROJECT$model$nopars[n],prod(dim(parameters)[2:3]))))
              states_all$nee_par_cor = cor(tmp,apply(states_all$nee_gCm2day,1,mean))
              states_all$gpp_par_cor = cor(tmp,apply(states_all$gpp_gCm2day,1,mean))
              states_all$rauto_par_cor = cor(tmp,apply(states_all$rauto_gCm2day,1,mean))
              states_all$rhet_par_cor = cor(tmp,apply(states_all$rhet_gCm2day,1,mean))
              # Avoid error flag when no fire
              if (max(as.vector(states_all$fire_gCm2day), na.rm=TRUE) > 0) {
                  states_all$fire_par_cor = cor(tmp,apply(states_all$fire_gCm2day,1,mean))
              } else {
                  states_all$fire_par_cor = array(0, dim = c(PROJECT$model$nopars[n],1))
              }
          }

          # pass to local variable for saving
          site_ctessel_pft = PROJECT$ctessel_pft[n]
          aNPP = states_all$aNPP ; natMTT = states_all$natMTT; MTT = states_all$MTT ; aMTT = states_all$aMTT ; SS = states_all$SS
          # Sanity check
          if (length(which(is.na(as.vector(aNPP))) == TRUE) > 0) {
              print(paste("NA value found in NPP for site ",PROJECT$site[n],sep="")) ; dummy = -1 ; return(dummy)
          }
          print("processing and storing ensemble output")
          # store the results now in binary file
          save(parameter_covariance,parameters,drivers,site_ctessel_pft,aNPP,natMTT,MTT,aMTT,SS,file=outfile1, compress="gzip")
          # determine whether this is a gridded run (or one with the override in place)
          if (PROJECT$spatial_type == "site" | grid_override == TRUE) {
              # ...if this is a site run save the full ensemble and everything else...
              save(parameters,drivers,states_all,site_ctessel_pft,file=outfile, compress="gzip")
          } else {
              # ...otherwise this is a grid and we want straight forward reduced dataset of common stocks and fluxes
              num_quantiles = c(0.025,0.05,0.25,0.5,0.75,0.95,0.975) ; num_quantiles_agg = seq(0.0,1, length = 1000)
              na_flag = TRUE
              # Estimate multiple use fluxes
              npp = states_all$gpp_gCm2day - states_all$rauto_gCm2day
              # Stocks first
              site_output = list(num_quantiles = num_quantiles,
                                 labile_gCm2 = apply(states_all$lab_gCm2,2,quantile,prob=num_quantiles,na.rm=na_flag),
#follab_gCm2 = apply(states_all$fol_gCm2+states_all$lab_gCm2,2,quantile,prob=num_quantiles,na.rm=na_flag),
                                 biomass_gCm2 = apply(states_all$bio_gCm2,2,quantile,prob=num_quantiles,na.rm=na_flag),
                                 lai_m2m2 = apply(states_all$lai_m2m2,2,quantile,prob=num_quantiles,na.rm=na_flag),
                                 foliage_gCm2 = apply(states_all$fol_gCm2,2,quantile,prob=num_quantiles,na.rm=na_flag),
                                 roots_gCm2 = apply(states_all$root_gCm2,2,quantile,prob=num_quantiles,na.rm=na_flag),
                                 wood_gCm2 = apply(states_all$wood_gCm2,2,quantile,prob=num_quantiles,na.rm=na_flag),
                                 lit_gCm2 = apply(states_all$lit_gCm2,2,quantile,prob=num_quantiles,na.rm=na_flag),
                                 som_gCm2 = apply(states_all$som_gCm2,2,quantile,prob=num_quantiles,na.rm=na_flag),
              # Fluxes second
                                 nee_gCm2day = apply(states_all$nee_gCm2day,2,quantile,prob=num_quantiles,na.rm=na_flag),
                                 gpp_gCm2day = apply(states_all$gpp_gCm2day,2,quantile,prob=num_quantiles,na.rm=na_flag),
                                 rauto_gCm2day = apply(states_all$rauto_gCm2day,2,quantile,prob=num_quantiles,na.rm=na_flag),
                                 rhet_gCm2day = apply(states_all$rhet_gCm2day,2,quantile,prob=num_quantiles,na.rm=na_flag),
                                 reco_gCm2day = apply(states_all$rhet_gCm2day+states_all$rauto_gCm2day,2,quantile,prob=num_quantiles,na.rm=na_flag),
                                 npp_gCm2day = apply(npp,2,quantile,prob=num_quantiles,na.rm=na_flag),
                                 fnpp_gCm2day = apply(npp * states_all$aNPP[,1],2,quantile,prob=num_quantiles,na.rm=na_flag),
                                 rnpp_gCm2day = apply(npp * states_all$aNPP[,2],2,quantile,prob=num_quantiles,na.rm=na_flag),
                                 wnpp_gCm2day = apply(npp * states_all$aNPP[,3],2,quantile,prob=num_quantiles,na.rm=na_flag),
                                 harvest_gCm2day = apply(states_all$harvest_C_gCm2day,2,quantile,prob=num_quantiles,na.rm=na_flag),
                                 fire_gCm2day = apply(states_all$fire_gCm2day,2,quantile,prob=num_quantiles,na.rm=na_flag))
              # Some derived fluxes / states
              site_output$nbe_gCm2day = apply(states_all$nee_gCm2day + states_all$fire_gCm2day,2,quantile,prob=num_quantiles,na.rm=na_flag)
              site_output$nbp_gCm2day = apply(-states_all$nee_gCm2day - states_all$fire_gCm2day - states_all$harvest_C_gCm2day,2,quantile,prob=num_quantiles,na.rm=na_flag)
              dCbio = states_all$bio_gCm2 - states_all$bio_gCm2[,1] # difference in biomass from initial
              site_output$dCbio_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm=na_flag)
              dCbio = states_all$fol_gCm2 - states_all$fol_gCm2[,1] # difference in biomass from initial
              site_output$dCfoliage_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm=na_flag)
              dCbio = states_all$root_gCm2 - states_all$root_gCm2[,1] # difference in root from initial
              site_output$dCroots_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm=na_flag)
              dCbio = states_all$wood_gCm2 - states_all$wood_gCm2[,1] # difference in wood from initial
              site_output$dCwood_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm=na_flag)
              dCbio = states_all$lit_gCm2 - states_all$lit_gCm2[,1] # difference in lit from initial
              site_output$dClit_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm=na_flag)
              dCbio = states_all$som_gCm2 - states_all$som_gCm2[,1] # difference in som from initial
              site_output$dCsom_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm=na_flag)
              if (length(which(names(states_all) == "litwood_gCm2")) > 0) {
                  # Total C including wood litter
                  site_output$totalC_gCm2 = apply(states_all$bio_gCm2+states_all$litwood_gCm2+states_all$lit_gCm2+states_all$som_gCm2,2,quantile,prob=num_quantiles,na.rm=na_flag)
                  dCbio = (states_all$bio_gCm2+states_all$litwood_gCm2+states_all$lit_gCm2+states_all$som_gCm2)
                  dCbio = dCbio - (states_all$bio_gCm2[,1]+states_all$litwood_gCm2[,1]+states_all$lit_gCm2[,1]+states_all$som_gCm2[,1])
                  site_output$dCtotalC_gCm2 = apply(dCbio,2,quantile, prob=num_quantiles, na.rm=TRUE)
                  # DOM (soil + litter + wood litter)
                  site_output$dom_gCm2 = apply(states_all$litwood_gCm2+states_all$lit_gCm2+states_all$som_gCm2,2,quantile,prob=num_quantiles,na.rm=na_flag)
                  dCbio = (states_all$litwood_gCm2+states_all$lit_gCm2+states_all$som_gCm2)
                  dCbio = dCbio - (states_all$litwood_gCm2[,1]+states_all$lit_gCm2[,1]+states_all$som_gCm2[,1])
                  site_output$dCdom_gCm2 = apply(dCbio,2,quantile, prob=num_quantiles, na.rm=TRUE)
                  # wood litter
                  site_output$litwood_gCm2 = apply(states_all$litwood_gCm2,2,quantile,prob=num_quantiles,na.rm=na_flag)
                  dCbio = states_all$litwood_gCm2 - states_all$litwood_gCm2[,1] # difference in cwd from initial
                  site_output$dClitwood_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm=na_flag)
                  # Combined litter and wood litter
                  site_output$deadorg_gCm2 = apply(states_all$litwood_gCm2+states_all$lit_gCm2,2,quantile,prob=num_quantiles,na.rm=na_flag)
                  dCbio = (states_all$litwood_gCm2+states_all$lit_gCm2) - (states_all$litwood_gCm2[,1]+states_all$lit_gCm2[,1]) # difference in deadorg from initial
                  site_output$dCdeadorg_gCm2 = apply(dCbio,2,quantile,prob=num_quantiles,na.rm=na_flag)

              } else {
                  # Total C without wood litter
                  site_output$totalC_gCm2 = apply(states_all$bio_gCm2+states_all$lit_gCm2+states_all$som_gCm2,2,quantile,prob=num_quantiles,na.rm=na_flag)
                  dCbio = (states_all$bio_gCm2+states_all$lit_gCm2+states_all$som_gCm2)
                  dCbio = dCbio - (states_all$bio_gCm2[,1]+states_all$lit_gCm2[,1]+states_all$som_gCm2[,1])
                  site_output$dCtotalC_gCm2 = apply(dCbio,2,quantile, prob=num_quantiles, na.rm=TRUE)
                  # DOM (soil + litter)
                  site_output$dom_gCm2 = apply(states_all$lit_gCm2+states_all$som_gCm2,2,quantile,prob=num_quantiles,na.rm=na_flag)
                  dCbio = (states_all$lit_gCm2+states_all$som_gCm2)
                  dCbio = dCbio - (states_all$lit_gCm2[,1]+states_all$som_gCm2[,1])
                  site_output$dCdom_gCm2 = apply(dCbio,2,quantile, prob=num_quantiles, na.rm=TRUE)
              }
              # Water cycle specific if available
              if (length(which(names(states_all) == "evap_kgH2Om2day")) > 0) {
                  # current water in the soil surface layer (0-30 cm)
                  site_output$SurfWater_kgH2Om2 = apply(states_all$sfc_water_mm,2,quantile,prob=num_quantiles,na.rm=na_flag)
                  # plant apparent soil water potential (MPa)
                  site_output$wSWP_MPa = apply(states_all$wSWP_MPa,2,quantile,prob=num_quantiles,na.rm=na_flag)
                  # evapotranspiration (Etrans + Esoil + Ewetcanopy)
                  site_output$evap_kgH2Om2day = apply(states_all$evap_kgH2Om2day,2,quantile,prob=num_quantiles,na.rm=na_flag)
              }
              # Canopy process information if available
              if (length(which(names(states_all) == "APAR_MJm2day")) > 0) {
                  # Extract the absorbed photosynthetically active radiation by the canopy
                  site_output$APAR_MJm2day = apply(states_all$APAR_MJm2day,2,quantile, prob=num_quantiles, na.rm=na_flag)
              }
              if (length(which(names(states_all) == "CiCa")) > 0) {
                  # Extract the internal vs ambient CO2 ratio
                  site_output$CiCa = apply(states_all$CiCa,2,quantile, prob=num_quantiles, na.rm=na_flag)
              }

              ## Now keeping the whole ensemble extract the pixel level values needed for grid scale aggregates
              ## All units remain at this point as the are output by DALEC
              steps_per_year = floor(dim(drivers$met)[1] / ((as.numeric(PROJECT$end_year) - as.numeric(PROJECT$start_year))+1))
              ss = dim(drivers$met)[1]-steps_per_year ;  ff = dim(drivers$met)[1]
              # Mean of final year
              if (length(which(names(states_all) == "litwood_gCm2")) > 0) {
                  ## With wood litter
                  # Total C and change
                  site_output$agg_totalC = quantile(apply(states_all$bio_gCm2[,ss:ff] +
                                                          states_all$lit_gCm2[,ss:ff] +
                                                          states_all$litwood_gCm2[,ss:ff] +
                                                          states_all$som_gCm2[,ss:ff],1,mean, na.rm=na_flag), prob = num_quantiles_agg, na.rm=na_flag)
                  dCbio = states_all$bio_gCm2[,ff] + states_all$lit_gCm2[,ff] + states_all$litwood_gCm2[,ff] + states_all$som_gCm2[,ff]
                  dCbio = dCbio - (states_all$bio_gCm2[,1] + states_all$lit_gCm2[,1] + states_all$litwood_gCm2[,1] + states_all$som_gCm2[,1])
                  site_output$agg_dCtotalC = quantile(dCbio, prob = num_quantiles_agg, na.rm=na_flag)
                  # DOM and change
                  site_output$agg_dom = quantile(apply(states_all$som_gCm2[,ss:ff] +
                                                       states_all$lit_gCm2[,ss:ff] +
                                                       states_all$litwood_gCm2[,ss:ff],1,mean, na.rm=na_flag), prob = num_quantiles_agg, na.rm=na_flag)
                  dCbio = states_all$lit_gCm2[,ff] + states_all$litwood_gCm2[,ff] + states_all$som_gCm2[,ff]
                  dCbio = dCbio - (states_all$lit_gCm2[,1] + states_all$litwood_gCm2[,1] + states_all$som_gCm2[,1])
                  site_output$agg_dCdom = quantile(dCbio, prob = num_quantiles_agg, na.rm=na_flag)
                  # Wood litter and change
                  site_output$agg_litwood = quantile(apply(states_all$litwood_gCm2[,ss:ff],1,mean, na.rm=na_flag), prob = num_quantiles_agg, na.rm=na_flag)
                  dCbio = states_all$litwood_gCm2[,ff] - states_all$litwood_gCm2[,1]
                  site_output$agg_dClitwood = quantile(dCbio, prob = num_quantiles_agg, na.rm=na_flag)
              } else {
                  ## Without wood litter
                  # Total C and change
                  site_output$agg_totalC = quantile(apply(states_all$bio_gCm2[,ss:ff] +
                                                          states_all$lit_gCm2[,ss:ff] +
                                                          states_all$som_gCm2[,ss:ff],1,mean, na.rm=na_flag), prob = num_quantiles_agg, na.rm=na_flag)
                  dCbio = states_all$bio_gCm2[,ff] + states_all$lit_gCm2[,ff] + states_all$som_gCm2[,ff]
                  dCbio = dCbio - (states_all$bio_gCm2[,1] + states_all$lit_gCm2[,1] + states_all$som_gCm2[,1])
                  site_output$agg_dCtotalC = quantile(dCbio, prob = num_quantiles_agg, na.rm=na_flag)
                  # DOM and change
                  site_output$agg_dom = quantile(apply(states_all$lit_gCm2[,ss:ff]+states_all$som_gCm2[,ss:ff],1,mean), prob = num_quantiles_agg, na.rm=na_flag)
                  dCbio = states_all$lit_gCm2[,ff] + states_all$som_gCm2[,ff]
                  dCbio = dCbio - (states_all$lit_gCm2[,1] + states_all$som_gCm2[,1])
                  site_output$agg_dCdom = quantile(dCbio, prob = num_quantiles_agg, na.rm=na_flag)

              }

              # Biomass and change
              site_output$agg_biomass = quantile(apply(states_all$bio_gCm2[,ss:ff],1,mean, na.rm=na_flag), prob = num_quantiles_agg, na.rm=na_flag)
              dCbio = states_all$bio_gCm2[,ff] - states_all$bio_gCm2[,1] # difference in biomass from initial
              site_output$agg_dCbio = quantile(dCbio, prob = num_quantiles_agg, na.rm=na_flag)
              # Labile and change
              site_output$agg_labile = quantile(apply(states_all$lab_gCm2[,ss:ff],1,mean, na.rm=na_flag), prob = num_quantiles_agg, na.rm=na_flag)
              dCbio = states_all$lab_gCm2[,ff] - states_all$lab_gCm2[,1] # difference in biomass from initial
              site_output$agg_dClabile = quantile(dCbio, prob = num_quantiles_agg, na.rm=na_flag)
              # Foliage and change
              site_output$agg_foliage = quantile(apply(states_all$fol_gCm2[,ss:ff],1,mean, na.rm=na_flag), prob = num_quantiles_agg, na.rm=na_flag)
              dCbio = states_all$fol_gCm2[,ff] - states_all$fol_gCm2[,1] # difference in biomass from initial
              site_output$agg_dCfoliage = quantile(dCbio, prob = num_quantiles_agg, na.rm=na_flag)
              # Fine root and change
              site_output$agg_root = quantile(apply(states_all$root_gCm2[,ss:ff],1,mean, na.rm=na_flag), prob = num_quantiles_agg, na.rm=na_flag)
              dCbio = states_all$root_gCm2[,ff] - states_all$root_gCm2[,1] # difference in biomass from initial
              site_output$agg_dCroot = quantile(dCbio, prob = num_quantiles_agg, na.rm=na_flag)
              # Wood and change
              site_output$agg_wood = quantile(apply(states_all$wood_gCm2[,ss:ff],1,mean, na.rm=na_flag), prob = num_quantiles_agg, na.rm=na_flag)
              dCbio = states_all$wood_gCm2[,ff] - states_all$wood_gCm2[,1] # difference in biomass from initial
              site_output$agg_dCwood = quantile(dCbio, prob = num_quantiles_agg, na.rm=na_flag)
              # Foliar and fine root litter and chaneg
              site_output$agg_lit = quantile(apply(states_all$lit_gCm2[,ss:ff],1,mean, na.rm=na_flag), prob = num_quantiles_agg, na.rm=na_flag)
              dCbio = states_all$lit_gCm2[,ff] - states_all$lit_gCm2[,1] # difference in biomass from initial
              site_output$agg_dClitter = quantile(dCbio, prob = num_quantiles_agg, na.rm=na_flag)
              # Soil C and change
              site_output$agg_som = quantile(apply(states_all$som_gCm2[,ss:ff],1,mean, na.rm=na_flag), prob = num_quantiles_agg, na.rm=na_flag)
              dCbio = states_all$som_gCm2[,ff] - states_all$som_gCm2[,1] # difference in biomass from initial
              site_output$agg_dCsom = quantile(dCbio, prob = num_quantiles_agg, na.rm=na_flag)

              ## Mean of whole time period
              site_output$agg_nbp = quantile(apply(-states_all$nee_gCm2day - states_all$fire_gCm2day - states_all$harvest_C_gCm2day,1,mean, na.rm=na_flag), prob=num_quantiles_agg, na.rm=na_flag)
              site_output$agg_nbe = quantile(apply(states_all$nee_gCm2day + states_all$fire_gCm2day,1,mean, na.rm=na_flag), prob=num_quantiles_agg, na.rm=na_flag)
              site_output$agg_nee = quantile(apply(states_all$nee_gCm2day,1,mean, na.rm=na_flag), prob = num_quantiles_agg, na.rm=na_flag)
              site_output$agg_gpp = quantile(apply(states_all$gpp_gCm2day,1,mean, na.rm=na_flag), prob = num_quantiles_agg, na.rm=na_flag)
              site_output$agg_rauto = quantile(apply(states_all$rauto_gCm2day,1,mean, na.rm=na_flag), prob = num_quantiles_agg, na.rm=na_flag)
              site_output$agg_rhet = quantile(apply(states_all$rhet_gCm2day,1,mean, na.rm=na_flag), prob = num_quantiles_agg, na.rm=na_flag)
              site_output$agg_reco = quantile(apply(states_all$reco_gCm2day,1,mean, na.rm=na_flag), prob = num_quantiles_agg, na.rm=na_flag)
              site_output$agg_harvest = quantile(apply(states_all$harvest_C_gCm2day,1,mean, na.rm=na_flag), prob = num_quantiles_agg, na.rm=na_flag)
              site_output$agg_fire = quantile(apply(states_all$fire_gCm2day,1,mean, na.rm=na_flag), prob = num_quantiles_agg, na.rm=na_flag)
              site_output$agg_npp = quantile(apply(npp,1,mean, na.rm=na_flag), prob = num_quantiles_agg, na.rm=na_flag)
              site_output$agg_fnpp = quantile(apply(npp * states_all$aNPP[,1],1,mean, na.rm=na_flag), prob = num_quantiles_agg, na.rm=na_flag)
              site_output$agg_rnpp = quantile(apply(npp * states_all$aNPP[,2],1,mean, na.rm=na_flag), prob = num_quantiles_agg, na.rm=na_flag)
              site_output$agg_wnpp = quantile(apply(npp * states_all$aNPP[,3],1,mean, na.rm=na_flag), prob = num_quantiles_agg, na.rm=na_flag)
              # Water cycle specific if available
              if (length(which(names(states_all) == "evap_kgH2Om2day")) > 0) {
                  # evapotranspiration (Etrans + Esoil + Ewetcanopy)
                  site_output$agg_evap = quantile(apply(states_all$evap_kgH2Om2day,1,mean, na.rm=na_flag), prob = num_quantiles_agg, na.rm=na_flag)
              }

              # C-cycle flux correlation with parameters
              site_output$nee_par_cor = states_all$nee_par_cor
              site_output$gpp_par_cor = states_all$gpp_par_cor
              site_output$rauto_par_cor = states_all$rauto_par_cor
              site_output$rhet_par_cor = states_all$rhet_par_cor
              site_output$fire_par_cor = states_all$fire_par_cor

              # save to pixel specific file for the moment... in "run_mcmc_results" these will be combined into a single grid
              save(site_output,file=outfile2, compress = "gzip")
        }

        dummy = 0

      } else {

        dummy = -1

      } # error_check == FALSE

  } else {

      dummy = -1
      print('Already extracted result vectors (set repair = 1 if re-run is needed)')

  } # *parameters.RData already exists

  return(dummy)

} # end of run_each_site

## Use byte compile
run_each_site<-cmpfun(run_each_site)

###
## Function to control the site level running of the re-processing by DALEC
###

run_mcmc_results <- function (PROJECT,stage,repair,grid_override) {

  print('Welcome to RUN_MCMC_RESULTS!!')

  # bundle needed functions down the chain
  functions_list=c("read_parameter_chains","read_binary_file_format","simulate_all",
                   "read_binary_response_surface","crop_development_parameters",
                   "have_chains_converged","psrf","read_parameter_covariance")
  # start marker
  stime = proc.time()["elapsed"]

  # how many plots in total do we have
  nos_plots = 1:PROJECT$nosites

  # now check which ones we need to calculate, but only if override not in play
  if (repair != 1) {
      print("...beginning filterings for sites we have already processed")
      keep_list = 0
      for (i in seq(1, length(nos_plots))) {
           outfile1 = paste(PROJECT$results_processedpath,PROJECT$sites[i],"_parameters.RData",sep="")
           if (file.exists(outfile1) == FALSE) {keep_list=append(keep_list,i)}
      }
      # filter out the sites we already have then
      keep_list = keep_list[-1] ; print(paste("......removing ",length(nos_plots)-length(keep_list)," sites out of ",length(nos_plots)," from the analysis",sep=""))
      nos_plots = nos_plots[keep_list]
  }

  # Combine driver modifiers into a list to be passed to the run_each_site()
  stage5modifiers=list(airt_factor = airt_factor, swrad_factor = swrad_factor, co2_factor = co2_factor,
                       rainfall_factor = rainfall_factor, wind_spd_factor = wind_spd_factor, vpd_factor = vpd_factor,
                       deforestation_factor = deforestation_factor, burnt_area_factor = burnt_area_factor)

  # now request the creation of the plots
  if (use_parallel & length(nos_plots) > 1) {
      print("...beginning parallel operations")
      cl <- makeCluster(min(length(nos_plots),numWorkers), type = "PSOCK")
      clusterExport(cl,functions_list)
      # load R libraries in cluster
      clusterExport(cl,"load_r_libraries") ; clusterEvalQ(cl, load_r_libraries())
      dummy=parLapply(cl,nos_plots,fun=run_each_site,PROJECT=PROJECT,stage=stage,
                      repair=repair,grid_override=grid_override,stage5modifiers=stage5modifiers)
      stopCluster(cl)
  } else {
      print("...beginning serial operations")
      # or use serial
      dummy=lapply(nos_plots,FUN=run_each_site,PROJECT=PROJECT,stage=stage,
                   repair=repair,grid_override=grid_override,stage5modifiers=stage5modifiers)
  } # parallel option

  # now if this is a gridded run we want to take out individual site specific summary files and combine them into a single file
  if (PROJECT$spatial_type == "grid" & grid_override == FALSE) {

      # update the user
      print("...combining pixel level stock / flux output into single array")

      # determine what the output file name is here, so that we can check if one already exists
      outfile_grid = paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep="")
      # read in example input files for using in some calculation
      drivers = read_binary_file_format(paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[1],".bin",sep=""))
      # Determine some useful information for the analysis below
      nos_years = (as.numeric(PROJECT$end_year) - as.numeric(PROJECT$start_year))+1
      steps_per_year = floor(dim(drivers$met)[1] / nos_years)
      # Set number of final aggregated quantiles wanted
      agg_quantiles_final = seq(0,1,length.out=100)

      if (file.exists(outfile_grid) == FALSE | repair == 1) {

          # make a list of all the files we will be reading in
          to_do = list.files(PROJECT$results_processedpath, full.names=TRUE)
          to_do = to_do[grepl("_stock_fluxes",to_do)]
          # read in the first file so that we can then set up all the grids for combined summary output
          load(to_do[1])

          #
          # Generate the summary information
          # Time invariant but contain uncertainty, shaped into the full spatial grid
          # i.e. including areas not part of the analysis but within the spatail domain
          #
          # Mean stocks first
          grid_output = list(num_quantiles = site_output$num_quantiles)
          grid_output$mean_labile_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_totalC_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_dom_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_biomass_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_foliage_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_roots_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_wood_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_lit_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_som_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_dCtotalC_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_dCdom_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_dCbio_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_dCfoliage_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_dCroots_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_dCwood_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_dClit_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_dCsom_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_lai_m2m2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          # Final stocks second
          grid_output$final_labile_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$final_totalC_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$final_dom_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$final_biomass_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$final_foliage_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$final_roots_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$final_wood_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$final_lit_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$final_som_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$final_dCtotalC_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$final_dCdom_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$final_dCbio_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$final_dCfoliage_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$final_dCroots_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$final_dCwood_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$final_dClit_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$final_dCsom_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$final_lai_m2m2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          # Fluxes third
          grid_output$mean_nee_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_gpp_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_rauto_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_rhet_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_reco_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_npp_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_fnpp_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_rnpp_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_wnpp_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_harvest_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_fire_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_nbe_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$mean_nbp_gCm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          # For those which we currently have need, estimate the mean annual maximum
          grid_output$annual_max_roots_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          grid_output$annual_max_wood_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          # Models where we have a wood litter pool and therefore a total dead organic matter combination also
          if (length(which(names(site_output) == "litwood_gCm2")) > 0) {
              grid_output$mean_litwood_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_dClitwood_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_litwood_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_dClitwood_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_deadorg_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$mean_dCdeadorg_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_deadorg_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_dCdeadorg_gCm2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          }
          # Finally water cycle specific if available
          if (length(which(names(site_output) == "evap_kgH2Om2day")) > 0) {
              # currently water in the soil surface layer (0-10 cm)
              grid_output$mean_SurfWater_kgH2Om2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_SurfWater_kgH2Om2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # plant apparent soil water potential (MPa)
              grid_output$mean_wSWP_MPa = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              grid_output$final_wSWP_MPa = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
              # evapotranspiration (Etrans + Esoil + Ewetcanopy)
              grid_output$mean_evap_kgH2Om2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          }
          # Canopy process variables
          if (length(which(names(site_output) == "APAR_MJm2day")) > 0) {
              # Absorbed photosynthetically active radation
              grid_output$mean_APAR_MJm2day = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          }
          if (length(which(names(site_output) == "CiCa")) > 0) {
              # Canopy Ci:Ca
              grid_output$mean_CiCa = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,dim(site_output$labile_gCm2)[1]))
          }

          # Time and uncertainty invarient information,
          # this is the correlation between ensemble members for parameter and C-cycle flux variables
          grid_output$nee_par_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)))
          grid_output$gpp_par_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)))
          grid_output$rauto_par_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)))
          grid_output$rhet_par_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)))
          grid_output$fire_par_cor = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)))

          #
          # Generate the variables needed to contain just the analysis locations
          # Full time series provided along with uncertainty information but for the analysis pixels only
          # Each pixel will have a corresponding i and j location which related to its lat/long within the overall domain grid
          #

          # Mean stocks first
          grid_output$labile_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$totalC_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$dom_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$biomass_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$foliage_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
#grid_output$follab_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$roots_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$wood_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$lit_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$som_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$dCtotalC_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$dCdom_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$dCbio_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$dCfoliage_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$dCroots_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$dCwood_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$dClit_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$dCsom_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$lai_m2m2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          # Fluxes seconds
          grid_output$nee_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$gpp_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$rauto_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$rhet_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$reco_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$npp_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$fnpp_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$rnpp_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$wnpp_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$harvest_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$fire_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$nbe_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          grid_output$nbp_gCm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          # Models where we have a CWD pool and therefore a total dead organic matter combination also
          if (length(which(names(site_output) == "litwood_gCm2")) > 0) {
              grid_output$litwood_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$dClitwood_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$dCdeadorg_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              grid_output$deadorg_gCm2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          }
          # Finally water cycle specific if available
          if (length(which(names(site_output) == "evap_kgH2Om2day")) > 0) {
              # currently water in the soil surface layer (0-10 cm)
              grid_output$SurfWater_kgH2Om2 = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # plant apparent soil water potential (MPa)
              grid_output$wSWP_MPa = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
              # evapotranspiration (Etrans + Esoil + Ewetcanopy)
              grid_output$evap_kgH2Om2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          }
          # Canopy process information
          if (length(which(names(site_output) == "APAR_MJm2day")) > 0) {
              # Canopy absorbed photosynthetically active radiation
              grid_output$APAR_MJm2day = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          }
          if (length(which(names(site_output) == "CiCa")) > 0) {
              # Canopy CiCa
              grid_output$CiCa = array(NA, dim=c(PROJECT$nosites,dim(site_output$labile_gCm2)[1],dim(site_output$labile_gCm2)[2]))
          }

          #
          # Generate the variables needed to contain grid aggregate values
          # These will be 'global' totals in units of TgC/yr, TgC, TgH20/yr or unit specific
          # Stocks variables will be estimated as the mean stock in the final year
          # Fluxes variables will be estimated as the mean flux over time
          #

          # Number of iterations available for resample
          agg_iter = length(site_output$agg_labile)

          # Mean stocks first
          grid_output$agg_labile_TgC = array(0, dim=c(agg_iter))
          grid_output$agg_totalC_TgC = array(0, dim=c(agg_iter))
          grid_output$agg_dom_TgC = array(0, dim=c(agg_iter))
          grid_output$agg_biomass_TgC = array(0, dim=c(agg_iter))
          grid_output$agg_foliage_TgC = array(0, dim=c(agg_iter))
          grid_output$agg_roots_TgC = array(0, dim=c(agg_iter))
          grid_output$agg_wood_TgC = array(0, dim=c(agg_iter))
          grid_output$agg_lit_TgC = array(0, dim=c(agg_iter))
          grid_output$agg_som_TgC = array(0, dim=c(agg_iter))
          grid_output$agg_dCtotalC_TgC = array(0, dim=c(agg_iter))
          grid_output$agg_dCdom_TgC = array(0, dim=c(agg_iter))
          grid_output$agg_dCbio_TgC = array(0, dim=c(agg_iter))
          grid_output$agg_dCfoliage_TgC = array(0, dim=c(agg_iter))
          grid_output$agg_dCroots_TgC = array(0, dim=c(agg_iter))
          grid_output$agg_dCwood_TgC = array(0, dim=c(agg_iter))
          grid_output$agg_dClit_TgC = array(0, dim=c(agg_iter))
          grid_output$agg_dCsom_TgC = array(0, dim=c(agg_iter))
          # Fluxes seconds
          grid_output$agg_nee_TgCyr = array(0, dim=c(agg_iter))
          grid_output$agg_gpp_TgCyr = array(0, dim=c(agg_iter))
          grid_output$agg_rauto_TgCyr = array(0, dim=c(agg_iter))
          grid_output$agg_rhet_TgCyr = array(0, dim=c(agg_iter))
          grid_output$agg_reco_TgCyr = array(0, dim=c(agg_iter))
          grid_output$agg_npp_TgCyr = array(0, dim=c(agg_iter))
          grid_output$agg_fnpp_TgCyr = array(0, dim=c(agg_iter))
          grid_output$agg_rnpp_TgCyr = array(0, dim=c(agg_iter))
          grid_output$agg_wnpp_TgCyr = array(0, dim=c(agg_iter))
          grid_output$agg_harvest_TgCyr = array(0, dim=c(agg_iter))
          grid_output$agg_fire_TgCyr = array(0, dim=c(agg_iter))
          grid_output$agg_nbe_TgCyr = array(0, dim=c(agg_iter))
          grid_output$agg_nbp_TgCyr = array(0, dim=c(agg_iter))
          # Models where we have a CWD pool and therefore a total dead organic matter combination also
          if (length(which(names(site_output) == "litwood_gCm2")) > 0) {
              grid_output$agg_litwood_TgC = array(0, dim=c(agg_iter))
              grid_output$agg_dClitwood_TgC = array(0, dim=c(agg_iter))
          }
          # Finally water cycle specific if available
          if (length(which(names(site_output) == "evap_kgH2Om2day")) > 0) {
              # evapotranspiration (Etrans + Esoil + Ewetcanopy)
              grid_output$agg_evap_PgH2Oyr = array(0, dim=c(agg_iter))
          }

          #
          # Generate the spatial information needed to relate i,j within the grid to long / lat and the summary vs detailed output
          #

          # vector into which the i,j position information within the grid will be stored
          grid_output$i_location = rep(NA, length.out = PROJECT$nosites)
          grid_output$j_location = rep(NA, length.out = PROJECT$nosites)

          # generate the lat / long grid again
          output = generate_wgs84_grid(PROJECT$latitude,PROJECT$longitude,PROJECT$resolution)
          grid_output$lat = array(output$lat, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
          grid_output$long = array(output$long,dim=c(PROJECT$long_dim,PROJECT$lat_dim))

          # Determine grid area (m2)
          grid_output$area = calc_pixel_area(grid_output$lat,grid_output$long,PROJECT$resolution)
          # this output is in vector form and we need matching array shapes so...
          grid_output$area = array(grid_output$area, dim=c(PROJECT$long_dim,PROJECT$lat_dim))

      } else {

          # if file already exists load it up
          load(outfile_grid)

      } # have we already an output file

      for (n in seq(1,PROJECT$nosites)) {
           # load the current file
           outfile2 = paste(PROJECT$results_processedpath,PROJECT$sites[n],"_stock_fluxes.RData",sep="")
           if (file.exists(outfile2)) {
               # load the site
               load(outfile2)
               # determine the lat / long location within the grid
               slot_j=as.numeric(PROJECT$sites[n])/PROJECT$long_dim
               slot_i=as.numeric(PROJECT$sites[n])-(floor(slot_j)*PROJECT$long_dim)
               if(slot_i == 0) {slot_i = PROJECT$long_dim} ; slot_j=ceiling(slot_j)
               # save for later
               grid_output$i_location[n] = slot_i ; grid_output$j_location[n] = slot_j
               # Store a record of the quantiles extracted for the site data
               grid_output$quantiles_wanted

               ## Begin aggregation to grid totals

               # Number of iterations available for resample
               agg_iter = length(site_output$agg_labile)
               # Unit adjustment
               unit_adj = 1e-12 # gC --> TgC ; or kgH2O -> PgH2O

               # Mean stocks first
               grid_output$agg_labile_TgC  = grid_output$agg_labile_TgC  + sample(site_output$agg_labile*grid_output$area[slot_i,slot_j]*unit_adj, size = agg_iter, replace = FALSE)
               grid_output$agg_totalC_TgC  = grid_output$agg_totalC_TgC  + sample(site_output$agg_totalC*grid_output$area[slot_i,slot_j]*unit_adj, size = agg_iter, replace = FALSE)
               grid_output$agg_dom_TgC     = grid_output$agg_dom_TgC     + sample(site_output$agg_dom*grid_output$area[slot_i,slot_j]*unit_adj, size = agg_iter, replace = FALSE)
               grid_output$agg_biomass_TgC = grid_output$agg_biomass_TgC + sample(site_output$agg_biomass*grid_output$area[slot_i,slot_j]*unit_adj, size = agg_iter, replace = FALSE)
               grid_output$agg_foliage_TgC = grid_output$agg_foliage_TgC + sample(site_output$agg_foliage*grid_output$area[slot_i,slot_j]*unit_adj, size = agg_iter, replace = FALSE)
               grid_output$agg_roots_TgC   = grid_output$agg_roots_TgC   + sample(site_output$agg_root*grid_output$area[slot_i,slot_j]*unit_adj, size = agg_iter, replace = FALSE)
               grid_output$agg_wood_TgC    = grid_output$agg_wood_TgC    + sample(site_output$agg_wood*grid_output$area[slot_i,slot_j]*unit_adj, size = agg_iter, replace = FALSE)
               grid_output$agg_lit_TgC     = grid_output$agg_lit_TgC     + sample(site_output$agg_lit*grid_output$area[slot_i,slot_j]*unit_adj, size = agg_iter, replace = FALSE)
               grid_output$agg_som_TgC     = grid_output$agg_som_TgC     + sample(site_output$agg_som*grid_output$area[slot_i,slot_j]*unit_adj, size = agg_iter, replace = FALSE)
               # Stock changes second
               grid_output$agg_dCtotalC_TgC  = grid_output$agg_dCtotalC_TgC + sample(site_output$agg_dCtotalC*grid_output$area[slot_i,slot_j]*unit_adj, size = agg_iter, replace = FALSE)
               grid_output$agg_dCdom_TgC     = grid_output$agg_dCdom_TgC + sample(site_output$agg_dCdom*grid_output$area[slot_i,slot_j]*unit_adj, size = agg_iter, replace = FALSE)
               grid_output$agg_dCbio_TgC     = grid_output$agg_dCbio_TgC + sample(site_output$agg_dCbio*grid_output$area[slot_i,slot_j]*unit_adj, size = agg_iter, replace = FALSE)
               grid_output$agg_dCfoliage_TgC = grid_output$agg_dCfoliage_TgC + sample(site_output$agg_dCfoliage*grid_output$area[slot_i,slot_j]*unit_adj, size = agg_iter, replace = FALSE)
               grid_output$agg_dCroots_TgC   = grid_output$agg_dCroots_TgC + sample(site_output$agg_dCroot*grid_output$area[slot_i,slot_j]*unit_adj, size = agg_iter, replace = FALSE)
               grid_output$agg_dCwood_TgC    = grid_output$agg_dCwood_TgC + sample(site_output$agg_dCwood*grid_output$area[slot_i,slot_j]*unit_adj, size = agg_iter, replace = FALSE)
               grid_output$agg_dClit_TgC     = grid_output$agg_dClit_TgC + sample(site_output$agg_dClitter*grid_output$area[slot_i,slot_j]*unit_adj, size = agg_iter, replace = FALSE)
               grid_output$agg_dCsom_TgC     = grid_output$agg_dCsom_TgC + sample(site_output$agg_dCsom*grid_output$area[slot_i,slot_j]*unit_adj, size = agg_iter, replace = FALSE)
               # Fluxes third
               grid_output$agg_nee_TgCyr = grid_output$agg_nee_TgC + sample(site_output$agg_nee*grid_output$area[slot_i,slot_j]*unit_adj*365.25, size = agg_iter, replace = FALSE)
               grid_output$agg_gpp_TgCyr = grid_output$agg_gpp_TgC + sample(site_output$agg_gpp*grid_output$area[slot_i,slot_j]*unit_adj*365.25, size = agg_iter, replace = FALSE)
               grid_output$agg_rauto_TgCyr = grid_output$agg_rauto_TgC + sample(site_output$agg_rauto*grid_output$area[slot_i,slot_j]*unit_adj*365.25, size = agg_iter, replace = FALSE)
               grid_output$agg_rhet_TgCyr = grid_output$agg_rhet_TgC + sample(site_output$agg_rhet*grid_output$area[slot_i,slot_j]*unit_adj*365.25, size = agg_iter, replace = FALSE)
               grid_output$agg_reco_TgCyr = grid_output$agg_reco_TgC + sample(site_output$agg_reco*grid_output$area[slot_i,slot_j]*unit_adj*365.25, size = agg_iter, replace = FALSE)
               grid_output$agg_npp_TgCyr = grid_output$agg_npp_TgC + sample(site_output$agg_npp*grid_output$area[slot_i,slot_j]*unit_adj*365.25, size = agg_iter, replace = FALSE)
               grid_output$agg_fnpp_TgCyr = grid_output$agg_fnpp_TgC + sample(site_output$agg_fnpp*grid_output$area[slot_i,slot_j]*unit_adj*365.25, size = agg_iter, replace = FALSE)
               grid_output$agg_rnpp_TgCyr = grid_output$agg_rnpp_TgC + sample(site_output$agg_rnpp*grid_output$area[slot_i,slot_j]*unit_adj*365.25, size = agg_iter, replace = FALSE)
               grid_output$agg_wnpp_TgCyr = grid_output$agg_wnpp_TgC + sample(site_output$agg_wnpp*grid_output$area[slot_i,slot_j]*unit_adj*365.25, size = agg_iter, replace = FALSE)
               grid_output$agg_harvest_TgCyr = grid_output$agg_harvest_TgC + sample(site_output$agg_harvest*grid_output$area[slot_i,slot_j]*unit_adj*365.25, size = agg_iter, replace = FALSE)
               grid_output$agg_fire_TgCyr = grid_output$agg_fire_TgC + sample(site_output$agg_fire*grid_output$area[slot_i,slot_j]*unit_adj*365.25, size = agg_iter, replace = FALSE)
               grid_output$agg_nbe_TgCyr = grid_output$agg_nbe_TgC + sample(site_output$agg_nbe*grid_output$area[slot_i,slot_j]*unit_adj*365.25, size = agg_iter, replace = FALSE)
               grid_output$agg_nbp_TgCyr = grid_output$agg_nbp_TgC + sample(site_output$agg_nbp*grid_output$area[slot_i,slot_j]*unit_adj*365.25, size = agg_iter, replace = FALSE)
               # Models where we have a CWD pool and therefore a total dead organic matter combination also
               if (length(which(names(site_output) == "litwood_gCm2")) > 0) {
                   grid_output$agg_litwood_TgC = grid_output$agg_litwood_TgC + sample(site_output$agg_litwood*grid_output$area[slot_i,slot_j]*unit_adj, size = agg_iter, replace = FALSE)
                   grid_output$agg_dClitwood_TgC = grid_output$agg_dClitwood_TgC + sample(site_output$agg_dClitwood*grid_output$area[slot_i,slot_j]*unit_adj, size = agg_iter, replace = FALSE)
               }
               # Finally water cycle specific if available
               if (length(which(names(site_output) == "evap_kgH2Om2day")) > 0) {
                   # evapotranspiration (Etrans + Esoil + Ewetcanopy)
                   grid_output$agg_evap_PgH2Oyr = grid_output$agg_evap_PgH2Oyr + sample(site_output$agg_evap*grid_output$area[slot_i,slot_j]*unit_adj*365.25, size = agg_iter, replace = FALSE)
               }
               # now assign to correct location in array
               # Stocks first
               grid_output$labile_gCm2[n,,] = site_output$labile_gCm2
               grid_output$totalC_gCm2[n,,] = site_output$totalC_gCm2
               grid_output$dom_gCm2[n,,] = site_output$dom_gCm2
               grid_output$biomass_gCm2[n,,] = site_output$biomass_gCm2
#grid_output$follab_gCm2[n,,] = site_output$follab_gCm2
               grid_output$foliage_gCm2[n,,] = site_output$foliage_gCm2
               grid_output$roots_gCm2[n,,] = site_output$roots_gCm2
               grid_output$wood_gCm2[n,,] = site_output$wood_gCm2
               grid_output$lit_gCm2[n,,] = site_output$lit_gCm2
               grid_output$som_gCm2[n,,] = site_output$som_gCm2
               grid_output$dCtotalC_gCm2[n,,] = site_output$dCtotalC_gCm2
               grid_output$dCdom_gCm2[n,,] = site_output$dCdom_gCm2
               grid_output$dCbio_gCm2[n,,] = site_output$dCbio_gCm2
               grid_output$dCfoliage_gCm2[n,,] = site_output$dCfoliage_gCm2
               grid_output$dCroots_gCm2[n,,] = site_output$dCroots_gCm2
               grid_output$dCwood_gCm2[n,,] = site_output$dCwood_gCm2
               grid_output$dClit_gCm2[n,,] = site_output$dClit_gCm2
               grid_output$dCsom_gCm2[n,,] = site_output$dCsom_gCm2
               grid_output$lai_m2m2[n,,] = site_output$lai_m2m2
               # Fluxes second
               grid_output$nee_gCm2day[n,,] = site_output$nee_gCm2day
               grid_output$gpp_gCm2day[n,,] = site_output$gpp_gCm2day
               grid_output$rauto_gCm2day[n,,] = site_output$rauto_gCm2day
               grid_output$rhet_gCm2day[n,,] = site_output$rhet_gCm2day
               grid_output$reco_gCm2day[n,,] = site_output$reco_gCm2day
               grid_output$npp_gCm2day[n,,] = site_output$npp_gCm2day
               grid_output$fnpp_gCm2day[n,,] = site_output$fnpp_gCm2day
               grid_output$rnpp_gCm2day[n,,] = site_output$rnpp_gCm2day
               grid_output$wnpp_gCm2day[n,,] = site_output$wnpp_gCm2day
               grid_output$harvest_gCm2day[n,,] = site_output$harvest_gCm2day
               grid_output$fire_gCm2day[n,,] = site_output$fire_gCm2day
               grid_output$nbe_gCm2day[n,,] = site_output$nbe_gCm2day
               grid_output$nbp_gCm2day[n,,] = site_output$nbp_gCm2day
               # Models where we have a CWD pool and therefore a total dead organic matter combination also
               if (length(which(names(site_output) == "litwood_gCm2")) > 0) {
                   grid_output$litwood_gCm2[n,,] = site_output$litwood_gCm2
                   grid_output$dClitwood_gCm2[n,,] = site_output$dClitwood_gCm2
                   grid_output$dCdeadorg_gCm2[n,,] = site_output$dCdeadorg_gCm2
                   grid_output$deadorg_gCm2[n,,] = site_output$deadorg_gCm2
               }
               # Finally water cycle specific if available
               if (length(which(names(site_output) == "evap_kgH2Om2day")) > 0) {
                   # currently water in the soil surface layer (0-10 cm)
                   grid_output$SurfWater_kgH2Om2[n,,] = site_output$SurfWater_kgH2Om2
                   # plant apparent soil water potential (MPa)
                   grid_output$wSWP_MPa[n,,] = site_output$wSWP_MPa
                   # evapotranspiration (Etrans + Esoil + Ewetcanopy)
                   grid_output$evap_kgH2Om2day[n,,] = site_output$evap_kgH2Om2day
               }
               if (length(which(names(site_output) == "APAR_MJm2day")) > 0) {
                  grid_output$APAR_MJm2day[n,,] = site_output$APAR_MJm2day
               }
               if (length(which(names(site_output) == "CiCa")) > 0) {
                  grid_output$CiCa[n,,] = site_output$CiCa
               }

               # now assign to correct location in array
               # Mean stocks first
               grid_output$mean_labile_gCm2[slot_i,slot_j,] = apply(site_output$labile_gCm2,1,mean)
               grid_output$mean_totalC_gCm2[slot_i,slot_j,] = apply(site_output$totalC_gCm2,1,mean)
               grid_output$mean_dom_gCm2[slot_i,slot_j,] = apply(site_output$dom_gCm2,1,mean)
               grid_output$mean_biomass_gCm2[slot_i,slot_j,] = apply(site_output$biomass_gCm2,1,mean)
               grid_output$mean_foliage_gCm2[slot_i,slot_j,] = apply(site_output$foliage_gCm2,1,mean)
               grid_output$mean_roots_gCm2[slot_i,slot_j,] = apply(site_output$roots_gCm2,1,mean)
               grid_output$mean_wood_gCm2[slot_i,slot_j,] = apply(site_output$wood_gCm2,1,mean)
               grid_output$mean_lit_gCm2[slot_i,slot_j,] = apply(site_output$lit_gCm2,1,mean)
               grid_output$mean_som_gCm2[slot_i,slot_j,] = apply(site_output$som_gCm2,1,mean)
               grid_output$mean_dCtotalC_gCm2[slot_i,slot_j,] = apply(site_output$dCtotalC_gCm2,1,mean)
               grid_output$mean_dCdom_gCm2[slot_i,slot_j,] = apply(site_output$dCdom_gCm2,1,mean)
               grid_output$mean_dCbio_gCm2[slot_i,slot_j,] = apply(site_output$dCbio_gCm2,1,mean)
               grid_output$mean_dCfoliage_gCm2[slot_i,slot_j,] = apply(site_output$dCfoliage_gCm2,1,mean)
               grid_output$mean_dCroots_gCm2[slot_i,slot_j,] = apply(site_output$dCroots_gCm2,1,mean)
               grid_output$mean_dCwood_gCm2[slot_i,slot_j,] = apply(site_output$dCwood_gCm2,1,mean)
               grid_output$mean_dClit_gCm2[slot_i,slot_j,] = apply(site_output$dClit_gCm2,1,mean)
               grid_output$mean_dCsom_gCm2[slot_i,slot_j,] = apply(site_output$dCsom_gCm2,1,mean)
               grid_output$mean_lai_m2m2[slot_i,slot_j,] = apply(site_output$lai_m2m2,1,mean)
               # Final stocks second
               final_step = dim(site_output$labile_gCm2)[2]
               grid_output$final_labile_gCm2[slot_i,slot_j,] = site_output$labile_gCm2[,final_step]
               grid_output$final_totalC_gCm2[slot_i,slot_j,] = site_output$totalC_gCm2[,final_step]
               grid_output$final_dom_gCm2[slot_i,slot_j,] = site_output$dom_gCm2[,final_step]
               grid_output$final_biomass_gCm2[slot_i,slot_j,] = site_output$biomass_gCm2[,final_step]
               grid_output$final_foliage_gCm2[slot_i,slot_j,] = site_output$foliage_gCm2[,final_step]
               grid_output$final_roots_gCm2[slot_i,slot_j,] = site_output$roots_gCm2[,final_step]
               grid_output$final_wood_gCm2[slot_i,slot_j,] = site_output$wood_gCm2[,final_step]
               grid_output$final_lit_gCm2[slot_i,slot_j,] = site_output$lit_gCm2[,final_step]
               grid_output$final_som_gCm2[slot_i,slot_j,] = site_output$som_gCm2[,final_step]
               grid_output$final_dCtotalC_gCm2[slot_i,slot_j,] = site_output$dCtotalC_gCm2[,final_step]
               grid_output$final_dCdom_gCm2[slot_i,slot_j,] = site_output$dCdom_gCm2[,final_step]
               grid_output$final_dCbio_gCm2[slot_i,slot_j,] = site_output$dCbio_gCm2[,final_step]
               grid_output$final_dCfoliage_gCm2[slot_i,slot_j,] = site_output$dCfoliage_gCm2[,final_step]
               grid_output$final_dCroots_gCm2[slot_i,slot_j,] = site_output$dCroots_gCm2[,final_step]
               grid_output$final_dCwood_gCm2[slot_i,slot_j,] = site_output$dCwood_gCm2[,final_step]
               grid_output$final_dClit_gCm2[slot_i,slot_j,] = site_output$dClit_gCm2[,final_step]
               grid_output$final_dCsom_gCm2[slot_i,slot_j,] = site_output$dCsom_gCm2[,final_step]
               grid_output$final_lai_m2m2[slot_i,slot_j,] = site_output$lai_m2m2[,final_step]
               # Fluxes third
               grid_output$mean_gpp_gCm2day[slot_i,slot_j,] = apply(site_output$gpp_gCm2day,1,mean)
               grid_output$mean_nee_gCm2day[slot_i,slot_j,] = apply(site_output$nee_gCm2day,1,mean)
               grid_output$mean_rauto_gCm2day[slot_i,slot_j,] = apply(site_output$rauto_gCm2day,1,mean)
               grid_output$mean_rhet_gCm2day[slot_i,slot_j,] = apply(site_output$rhet_gCm2day,1,mean)
               grid_output$mean_reco_gCm2day[slot_i,slot_j,] = apply(site_output$reco_gCm2day,1,mean)
               grid_output$mean_npp_gCm2day[slot_i,slot_j,] = apply(site_output$npp_gCm2day,1,mean)
               grid_output$mean_fnpp_gCm2day[slot_i,slot_j,] = apply(site_output$fnpp_gCm2day,1,mean)
               grid_output$mean_rnpp_gCm2day[slot_i,slot_j,] = apply(site_output$rnpp_gCm2day,1,mean)
               grid_output$mean_wnpp_gCm2day[slot_i,slot_j,] = apply(site_output$wnpp_gCm2day,1,mean)
               grid_output$mean_harvest_gCm2day[slot_i,slot_j,] = apply(site_output$harvest_gCm2day,1,mean)
               grid_output$mean_fire_gCm2day[slot_i,slot_j,] = apply(site_output$fire_gCm2day,1,mean)
               grid_output$mean_nbe_gCm2day[slot_i,slot_j,] = apply(site_output$nbe_gCm2day,1,mean)
               grid_output$mean_nbp_gCm2day[slot_i,slot_j,] = apply(site_output$nbp_gCm2day,1,mean)
               # For those which we currently have need, estimate the mean annual maximum
               grid_output$annual_max_roots_gCm2[slot_i,slot_j,] = 0
               grid_output$annual_max_wood_gCm2[slot_i,slot_j,] = 0
               for (y in seq(1,nos_years)) {
                    s = 1 + (steps_per_year * (y - 1)) ; f = s + (steps_per_year-1)
                    grid_output$annual_max_roots_gCm2[slot_i,slot_j,] = grid_output$annual_max_roots_gCm2[slot_i,slot_j,] + apply(site_output$roots_gCm2[,s:f],1,max)
                    grid_output$annual_max_wood_gCm2[slot_i,slot_j,] = grid_output$annual_max_wood_gCm2[slot_i,slot_j,] + apply(site_output$wood_gCm2[,s:f],1,max)
               }
               grid_output$annual_max_roots_gCm2[slot_i,slot_j,] = grid_output$annual_max_roots_gCm2[slot_i,slot_j,] / nos_years
               grid_output$annual_max_wood_gCm2[slot_i,slot_j,] = grid_output$annual_max_wood_gCm2[slot_i,slot_j,] / nos_years
               # Models where we have a CWD pool and therefore a total dead organic matter combination also
               if (length(which(names(site_output) == "litwood_gCm2")) > 0) {
                   grid_output$mean_litwood_gCm2[slot_i,slot_j,] = apply(site_output$litwood_gCm2,1,mean)
                   grid_output$mean_dClitwood_gCm2[slot_i,slot_j,] = apply(site_output$dClitwood_gCm2,1,mean)
                   grid_output$mean_dCdeadorg_gCm2[slot_i,slot_j,] = apply(site_output$dCdeadorg_gCm2,1,mean)
                   grid_output$mean_deadorg_gCm2[slot_i,slot_j,] = apply(site_output$deadorg_gCm2,1,mean)
                   grid_output$final_litwood_gCm2[slot_i,slot_j,] = site_output$litwood_gCm2[,final_step]
                   grid_output$final_dClitwood_gCm2[slot_i,slot_j,] = site_output$dClitwood_gCm2[,final_step]
                   grid_output$final_dCdeadorg_gCm2[slot_i,slot_j,] = site_output$dCdeadorg_gCm2[,final_step]
                   grid_output$final_deadorg_gCm2[slot_i,slot_j,] = site_output$deadorg_gCm2[,final_step]
               }
               # Finally water cycle specific if available
               if (length(which(names(site_output) == "evap_kgH2Om2day")) > 0) {
                   # currently water in the soil surface layer (0-10 cm)
                   grid_output$mean_SurfWater_kgH2Om2[slot_i,slot_j,] = apply(site_output$SurfWater_kgH2Om2,1,mean)
                   grid_output$final_SurfWater_kgH2Om2[slot_i,slot_j,] = site_output$SurfWater_kgH2Om2[,final_step]
                   # plant apparent soil water potential (MPa)
                   grid_output$mean_wSWP_MPa[slot_i,slot_j,] = apply(site_output$wSWP_MPa,1,mean)
                   grid_output$final_wSWP_MPa[slot_i,slot_j,] = site_output$wSWP_MPa[,final_step]
                   # evapotranspiration (Etrans + Esoil + Ewetcanopy)
                   grid_output$mean_evap_kgH2Om2day[slot_i,slot_j,] = apply(site_output$evap_kgH2Om2day,1,mean)
               }
               # Canopy process information
               if (length(which(names(site_output) == "APAR_MJm2day")) > 0) {
                   grid_output$mean_APAR_MJm2day[slot_i,slot_j,] = apply(site_output$APAR_MJm2day,1,mean)
               }
               if (length(which(names(site_output) == "CiCa")) > 0) {
                   grid_output$mean_CiCa[slot_i,slot_j,] = apply(site_output$CiCa,1,mean)
               }

               # Parameter vs C-cycle flux correlation across ensemble member
               grid_output$nee_par_cor[slot_i,slot_j,] = site_output$nee_par_cor
               grid_output$gpp_par_cor[slot_i,slot_j,] = site_output$gpp_par_cor
               grid_output$rauto_par_cor[slot_i,slot_j,] = site_output$rauto_par_cor
               grid_output$rhet_par_cor[slot_i,slot_j,] = site_output$rhet_par_cor
               grid_output$fire_par_cor[slot_i,slot_j,] = site_output$fire_par_cor

               # now tidy away the file
               file.remove(outfile2)

          } # if site has been done and can now be made part of the overall...

      } # looping through sites

      # Mean stocks first
      grid_output$agg_labile_TgC = quantile(grid_output$agg_labile_TgC, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_totalC_TgC = quantile(grid_output$agg_totalC_TgC, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_dom_TgC = quantile(grid_output$agg_dom_TgC, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_biomass_TgC = quantile(grid_output$agg_biomass_TgC, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_foliage_TgC = quantile(grid_output$agg_foliage_TgC, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_roots_TgC = quantile(grid_output$agg_roots_TgC, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_wood_TgC = quantile(grid_output$agg_wood_TgC, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_lit_TgC = quantile(grid_output$agg_lit_TgC, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_som_TgC = quantile(grid_output$agg_som_TgC, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_dCtotalC_TgC = quantile(grid_output$agg_dCtotalC_TgC, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_dCdom_TgC = quantile(grid_output$agg_dCdom_TgC, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_dCbio_TgC = quantile(grid_output$agg_dCbio_TgC, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_dCfoliage_TgC = quantile(grid_output$agg_dCfoliage_TgC, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_dCroots_TgC = quantile(grid_output$agg_dCroots_TgC, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_dCwood_TgC = quantile(grid_output$agg_dCwood_TgC, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_dClit_TgC = quantile(grid_output$agg_dClit_TgC, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_dCsom_TgC = quantile(grid_output$agg_dCsom_TgC, prob = agg_quantiles_final, na.rm=TRUE)
      # Fluxes seconds
      grid_output$agg_nee_TgCyr = quantile(grid_output$agg_nee_TgCyr, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_gpp_TgCyr = quantile(grid_output$agg_gpp_TgCyr, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_rauto_TgCyr = quantile(grid_output$agg_rauto_TgCyr, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_rhet_TgCyr = quantile(grid_output$agg_rhet_TgCyr, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_reco_TgCyr = quantile(grid_output$agg_reco_TgCyr, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_npp_TgCyr = quantile(grid_output$agg_npp_TgCyr, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_fnpp_TgCyr = quantile(grid_output$agg_fnpp_TgCyr, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_rnpp_TgCyr = quantile(grid_output$agg_rnpp_TgCyr, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_wnpp_TgCyr = quantile(grid_output$agg_wnpp_TgCyr, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_harvest_TgCyr = quantile(grid_output$agg_harvest_TgCyr, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_fire_TgCyr = quantile(grid_output$agg_fire_TgCyr, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_nbe_TgCyr = quantile(grid_output$agg_nbe_TgCyr, prob = agg_quantiles_final, na.rm=TRUE)
      grid_output$agg_nbp_TgCyr = quantile(grid_output$agg_nbp_TgCyr, prob = agg_quantiles_final, na.rm=TRUE)
      # Models where we have a CWD pool and therefore a total dead organic matter combination also
      if (length(which(names(grid_output) == "agg_litwood_TgC")) > 0) {
          grid_output$agg_litwood_TgC = quantile(grid_output$agg_litwood_TgC, prob = agg_quantiles_final, na.rm=TRUE)
          grid_output$agg_dClitwood_TgC = quantile(grid_output$agg_dClitwood_TgC, prob = agg_quantiles_final, na.rm=TRUE)
      }
      # Finally water cycle specific if available
      if (length(which(names(grid_output) == "agg_evap_PgH2Oyr")) > 0) {
          # evapotranspiration (Etrans + Esoil + Ewetcanopy)
          grid_output$agg_evap_PgH2Oyr = quantile(grid_output$agg_evap_PgH2Oyr, prob = agg_quantiles_final, na.rm=TRUE)
      }

      # now save the combined grid file
      save(grid_output, file=outfile_grid, compress = "gzip")

  } # gridded run?

  # tell me whats happening
  print(paste("...time to process ",round((proc.time()["elapsed"]-stime)/60,1)," minutes",sep=""))

} # end function run_mcmc_results

## Use byte compile
run_mcmc_results<-cmpfun(run_mcmc_results)
