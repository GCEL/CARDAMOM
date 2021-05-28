
###
## Function to create generic plots of gridded CARDAMOM analysis parameters and emergent traits
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

generate_parameter_maps<-function(PROJECT) {

  # how many years of the analysis
  nos_years = length(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year))
  # load timesteps to local variable
  timestep_days = PROJECT$model$timestep_days ; seconds_per_day = 86400
  if (length(timestep_days) == 1) {
      eh = TRUE ; n = 0
      while(eh) {
        n = n + 1
        # generate file name of the output file created in stage 3
        loadfile = paste(PROJECT$results_processedpath,PROJECT$sites[n],"_parameters.RData",sep="")
        if (file.exists(loadfile)) {load(loadfile) ; eh = FALSE}
    }
    timestep_days = rep(timestep_days, length.out=drivers$nodays)
  } # length(timestep_days) == 1

  # determine new output file for aggregated values
  outfile = paste(PROJECT$results_processedpath,PROJECT$name,"_parameter_maps.RData",sep="")

  # Which quantiles will we extract, these should be kept the same as those for the stock / flux outputs
  num_quantiles = c(0.025,0.05,0.25,0.5,0.75,0.95,0.975) ; na_flag = TRUE
  median_loc = 4 ; upper_loc = 7 ; lower_loc = 1

  ## Parameters needed to estimate the fire related residence time.
  # Combustion completeness parameters
  cf = rep(0.1,7)  # 0.1 applies to labile, roots and wood
  cf[2] = 0.9      # Update foliar
  cf[c(5,7)] = 0.7 # Update litter and wood litter
  cf[6] = 0.01     # Update soil
  # Resilience factor for non-combusted tissue
  rfac = rep(0.5,7) ; rfac[5] = 0.1 ; rfac[6] = 0 ; rfac[7] = 0.1

  if (file.exists(outfile) == FALSE | repair == 1) {
      # Create output object
      grid_parameters = list(pft_array=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim)))
      # Parameter information
      grid_parameters$parameters = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,(max(PROJECT$model$nopars)+1),length(num_quantiles)))
      grid_parameters$parameters_converged=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,max(PROJECT$model$nopars)+1))
      # Derived Ecosystem traits (time invarient)
      grid_parameters$NPP_foliar_fraction=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(num_quantiles)))
      grid_parameters$NPP_wood_fraction=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(num_quantiles)))
      grid_parameters$NPP_root_fraction=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(num_quantiles)))
      grid_parameters$MTT_foliar_years=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(num_quantiles)))
      grid_parameters$MTT_wood_years=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(num_quantiles)))
      grid_parameters$MTT_root_years=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(num_quantiles)))
      grid_parameters$MTT_som_years=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(num_quantiles)))
      grid_parameters$MTT_DeadOrg_years=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(num_quantiles)))
      grid_parameters$SS_foliar_gCm2=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(num_quantiles)))
      grid_parameters$SS_wood_gCm2=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(num_quantiles)))
      grid_parameters$SS_root_gCm2=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(num_quantiles)))
      grid_parameters$SS_DeadOrg_gCm2=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(num_quantiles)))
      grid_parameters$SS_som_gCm2=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(num_quantiles)))
      # Derived Ecoystem traits (time varient)
      grid_parameters$aMTT_foliar_years=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_years,length(num_quantiles)))
      grid_parameters$aMTT_wood_years=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_years,length(num_quantiles)))
      grid_parameters$aMTT_root_years=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_years,length(num_quantiles)))
      grid_parameters$aMTT_som_years=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_years,length(num_quantiles)))
      grid_parameters$aMTT_DeadOrg_years=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_years,length(num_quantiles)))
      # Derived Ecosystem traits (turnover partitioning)
      # Natural mean annual transit time
      grid_parameters$MTTnatural_foliar_years=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(num_quantiles)))
      grid_parameters$MTTnatural_wood_years=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(num_quantiles)))
      grid_parameters$MTTnatural_root_years=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(num_quantiles)))
      grid_parameters$MTTnatural_som_years=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(num_quantiles)))
      grid_parameters$MTTnatural_DeadOrg_years=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(num_quantiles)))
      # Fire
      grid_parameters$MTTfire_foliar_years=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(num_quantiles)))
      grid_parameters$MTTfire_wood_years=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(num_quantiles)))
      grid_parameters$MTTfire_root_years=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(num_quantiles)))
      # Biomass removal (harvest)
      grid_parameters$MTTharvest_foliar_years=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(num_quantiles)))
      grid_parameters$MTTharvest_wood_years=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(num_quantiles)))
      grid_parameters$MTTharvest_root_years=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(num_quantiles)))
      # Model Driver information
      grid_parameters$mean_temperature_C=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
      grid_parameters$mean_vpd_Pa=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
      grid_parameters$mean_radiation_MJm2day=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
      grid_parameters$mean_precipitation_kgm2yr=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
      # Observational information
      grid_parameters$obs_lai_max_m2m2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
      grid_parameters$obs_lai_mean_m2m2=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
      grid_parameters$obs_lai_sigma_m2m2=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
      grid_parameters$obs_wood_gCm2=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
      grid_parameters$obs_wood_unc_gCm2=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
      grid_parameters$obs_som_gCm2=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
      grid_parameters$obs_som_unc_gCm2=array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))

      # loop through all the files to eventually build up a complete vector
      for (n in seq(1, PROJECT$nosites)) {
           # generate file name of the output file created in stage 3
           loadfile = paste(PROJECT$results_processedpath,PROJECT$sites[n],"_parameters.RData",sep="")
           if (file.exists(loadfile) == TRUE) {

               # Load the specific file
               load(loadfile)

               # Sanity check
               if (length(which(is.na(as.vector(aNPP)) == TRUE)) > 0) {print(paste("NA found in aNPP site = ",loadfile,sep=""))}
               # Update user
               if (n < 100 | n%%100 == 0) {print(paste("...have loaded ",round((n/PROJECT$nosites)*100, digits=0),"% of pixels",sep=""))}

               # now load parameter information into the array
               slot_j=as.numeric(PROJECT$sites[n])/PROJECT$long_dim
               slot_i=as.numeric(PROJECT$sites[n])-(floor(slot_j)*PROJECT$long_dim)
               if (slot_i == 0){slot_i = PROJECT$long_dim} ; slot_j=ceiling(slot_j)
               grid_parameters$mean_temperature_C[slot_i,slot_j]=mean((drivers$met[,2]+drivers$met[,3])*0.5)
               grid_parameters$mean_radiation_MJm2day[slot_i,slot_j]=mean(drivers$met[,4])
               grid_parameters$mean_vpd_Pa[slot_i,slot_j]=mean(drivers$met[,16])
               grid_parameters$mean_precipitation_kgm2yr[slot_i,slot_j]=sum(drivers$met[,7]*timestep_days*seconds_per_day)/nos_years
               grid_parameters$obs_lai_max_m2m2[slot_i,slot_j]=max(drivers$obs[,3])
               grid_parameters$obs_lai_mean_m2m2[slot_i,slot_j]=mean(drivers$obs[which(drivers$obs[,3] != -9999),3])
               grid_parameters$obs_lai_sigma_m2m2[slot_i,slot_j]=sd(drivers$obs[which(drivers$obs[,3] != -9999),3])
               if (length(which(drivers$obs[,13] > -9999)) > 0) {
                   grid_parameters$obs_wood_gCm2[slot_i,slot_j]=max(drivers$obs[which(drivers$obs[,13] != -9999),13],na.rm=TRUE)
                   grid_parameters$obs_wood_unc_gCm2[slot_i,slot_j]=max(drivers$obs[which(drivers$obs[,14] != -9999),14],na.rm=TRUE)
               } else if (drivers$parpriors[21] > 0) {
                   grid_parameters$obs_wood_gCm2[slot_i,slot_j]=drivers$parpriors[21]
                   grid_parameters$obs_wood_unc_gCm2[slot_i,slot_j]=drivers$parpriorunc[21]
               }
               if (length(which(drivers$obs[,19] > -9999)) > 0) {
                   grid_parameters$obs_som_gCm2[slot_i,slot_j] = max (drivers$obs[which(drivers$obs[,19] != -9999),19],na.rm=TRUE)
                   grid_parameters$obs_som_unc_gCm2[slot_i,slot_j] = max(drivers$obs[which(drivers$obs[,20] != -9999),20],na.rm=TRUE)
               } else if (drivers$parpriors[23] > 0) {
                   grid_parameters$obs_som_gCm2[slot_i,slot_j] = drivers$parpriors[23]
                   grid_parameters$obs_som_unc_gCm2[slot_i,slot_j] = drivers$parpriorunc[23]
               }
               grid_parameters$parameters_converged[slot_i,slot_j,] = 0
               converged=have_chains_converged(parameters)
               grid_parameters$parameters_converged[slot_i,slot_j,which(converged=="PASS")]=1
               grid_parameters$pft_array[slot_i,slot_j]=drivers$ctessel_pft
               # calculate NPP allocation fraction (fol,root,wood)
               grid_parameters$NPP_foliar_fraction[slot_i,slot_j,]=quantile(aNPP[,1], prob=num_quantiles,na.rm=TRUE)
               grid_parameters$NPP_wood_fraction[slot_i,slot_j,]=quantile(aNPP[,2], prob=num_quantiles,na.rm=TRUE)
               grid_parameters$NPP_root_fraction[slot_i,slot_j,]=quantile(aNPP[,3], prob=num_quantiles,na.rm=TRUE)
               # calculate mean residence time variables (fol,root,wood,lit+litwood,som)
               # NOTE: that depending on the model DeadOrg may be litter or litter + cwd
               grid_parameters$MTT_foliar_years[slot_i,slot_j,]=quantile(MTT[,1], prob=num_quantiles,na.rm=TRUE)
               grid_parameters$MTT_root_years[slot_i,slot_j,]=quantile(MTT[,2], prob=num_quantiles,na.rm=TRUE)
               grid_parameters$MTT_wood_years[slot_i,slot_j,]=quantile(MTT[,3], prob=num_quantiles,na.rm=TRUE)
               grid_parameters$MTT_DeadOrg_years[slot_i,slot_j,]=quantile(MTT[,4], prob=num_quantiles,na.rm=TRUE)
               grid_parameters$MTT_som_years[slot_i,slot_j,]=quantile(MTT[,5], prob=num_quantiles,na.rm=TRUE)
               # calculate natural (i.e. without disturbance) mean residence time variables (fol,root,wood,lit+litwood,som)
               # NOTE: that depending on the model DeadOrg may be litter or litter + cwd
               grid_parameters$MTTnatural_foliar_years[slot_i,slot_j,]=quantile(natMTT[,1], prob=num_quantiles,na.rm=TRUE)
               grid_parameters$MTTnatural_root_years[slot_i,slot_j,]=quantile(natMTT[,2], prob=num_quantiles,na.rm=TRUE)
               grid_parameters$MTTnatural_wood_years[slot_i,slot_j,]=quantile(natMTT[,3], prob=num_quantiles,na.rm=TRUE)
               grid_parameters$MTTnatural_DeadOrg_years[slot_i,slot_j,]=quantile(natMTT[,4], prob=num_quantiles,na.rm=TRUE)
               grid_parameters$MTTnatural_som_years[slot_i,slot_j,]=quantile(natMTT[,5], prob=num_quantiles,na.rm=TRUE)
               # calculate annual residence time variables (fol,root,wood,lit+litwood,som)
               # NOTE: that depending on the model DeadOrg may be litter or litter + cwd
               grid_parameters$aMTT_foliar_years[slot_i,slot_j,,]=t(apply(aMTT[,1,],2,quantile, prob=num_quantiles,na.rm=TRUE))
               grid_parameters$aMTT_root_years[slot_i,slot_j,,]=t(apply(aMTT[,2,],2,quantile, prob=num_quantiles,na.rm=TRUE))
               grid_parameters$aMTT_wood_years[slot_i,slot_j,,]=t(apply(aMTT[,3,],2,quantile, prob=num_quantiles,na.rm=TRUE))
               grid_parameters$aMTT_DeadOrg_years[slot_i,slot_j,,]=t(apply(aMTT[,4,],2,quantile, prob=num_quantiles,na.rm=TRUE))
               grid_parameters$aMTT_som_years[slot_i,slot_j,,]=t(apply(aMTT[,5,],2,quantile, prob=num_quantiles,na.rm=TRUE))
               # calculate disturbance specific residence time variables (fol,root,wood,lit+litwood,som)
               # NOTE: that depending on the model DeadOrg may be litter or litter + cwd
               # Mean annual fire fraction
               tmp = sum(drivers$met[which(drivers$met[,9] > 0),9]) / nos_years
               # If there is fire we will estimate its impact on residence time
               if (tmp > 0) {
                   # (fire*cc) + (fire*(1-cc)*(1-rfac))
                   if (PROJECT$model$name == "DALEC_CDEA_ACM2_BUCKET_RmRg_CWD_wMRT") {
                       # Use the model calibrated resiliance factors
                       # Foliage
                       tmp2 = ((tmp * parameters[32,,]) + (tmp * (1-parameters[32,,]) * (1-parameters[31,,]))) ** -1
                       tmp2 = quantile(tmp2, prob = num_quantiles, na.rm=TRUE)
                       grid_parameters$MTTfire_foliar_years[slot_i,slot_j,] = tmp2
                       # Currently fine roots and wood have a commmon rfac and combustion completeness
                       tmp2 = ((tmp * parameters[33,,]) + (tmp * (1-parameters[33,,]) * (1-parameters[31,,]))) ** -1
                       tmp2 = quantile(tmp2, prob = num_quantiles, na.rm=TRUE)
                       grid_parameters$MTTfire_root_years[slot_i,slot_j,] = tmp2
                       grid_parameters$MTTfire_wood_years[slot_i,slot_j,] = tmp2
                   } else {
                       # Use default assumptions
                       tmp2 = ((tmp * cf[2]) + (tmp * (1-cf[2]) * (1-rfac[2]))) ** -1
                       tmp2 = quantile(tmp2, prob = num_quantiles, na.rm=TRUE)
                       grid_parameters$MTTfire_foliar_years[slot_i,slot_j,] = tmp2
                       tmp2 = ((tmp * cf[3]) + (tmp * (1-cf[3]) * (1-rfac[3]))) ** -1
                       tmp2 = quantile(tmp2, prob = num_quantiles, na.rm=TRUE)
                       grid_parameters$MTTfire_root_years[slot_i,slot_j,] = tmp2
                       tmp2 = ((tmp * cf[4]) + (tmp * (1-cf[4]) * (1-rfac[4]))) ** -1
                       tmp2 = quantile(tmp2, prob = num_quantiles, na.rm=TRUE)
                       grid_parameters$MTTfire_wood_years[slot_i,slot_j,] = tmp2
                   } # generic or model specific parameters
               } # is there any fire
               # Forest biomass removal / harvest
               tmp = sum(drivers$met[which(drivers$met[,8] > 0),8]) / nos_years
               if (tmp > 0) {
                   tmp = tmp ** -1
                   grid_parameters$MTTharvest_foliar_years[slot_i,slot_j,] = tmp
                   grid_parameters$MTTharvest_root_years[slot_i,slot_j,] = tmp
                   grid_parameters$MTTharvest_wood_years[slot_i,slot_j,] = tmp
               } # is there any disturbance
               # Calculate C stock steady states, as function of natural, fire and biomass extraction
               grid_parameters$SS_foliar_gCm2[slot_i,slot_j,]=quantile(SS[,1], prob=num_quantiles,na.rm=TRUE)
               grid_parameters$SS_root_gCm2[slot_i,slot_j,]=quantile(SS[,2], prob=num_quantiles,na.rm=TRUE)
               grid_parameters$SS_wood_gCm2[slot_i,slot_j,]=quantile(SS[,3], prob=num_quantiles,na.rm=TRUE)
               grid_parameters$SS_DeadOrg_gCm2[slot_i,slot_j,]=quantile(SS[,4], prob=num_quantiles,na.rm=TRUE)
               grid_parameters$SS_som_gCm2[slot_i,slot_j,]=quantile(SS[,5], prob=num_quantiles,na.rm=TRUE)
               # loop through parameters + likelihood
               for (p in seq(1, dim(parameters)[1])) {
                    grid_parameters$parameters[slot_i,slot_j,p,] = quantile(as.vector(parameters[p,,]), prob=num_quantiles)
               }

           } #
      } #  site loop

      # Convert avgN log10-normal to gN/m2
      if (PROJECT$model$name == "DALECN_GSI_FR" | PROJECT$model$name == "DALECN_GSI_DFOL_FR"
        | PROJECT$model$name == "DALEC_GSI_BUCKET"| PROJECT$model$name == "DALEC_GSI_DFOL_CWD_FR"
        | PROJECT$model$name == "DALECN_GSI_BUCKET" | PROJECT$model$name == "DALECN_BUCKET"
        | PROJECT$model$name == "DALEC_BUCKET" | PROJECT$model$name == "DALEC"
        | PROJECT$model$name == "DALEC_G5" | PROJECT$model$name == "DALEC_G6"
        | PROJECT$model$name == "DALEC_BUCKET_CanAGE") {
          grid_parameters$parameters[,,11,]=10**grid_parameters$parameters[,,11,]
      }

  } else {

      # if we have loaded these before just load the final product
      load(outfile)

  } # have these previously been loaded

  # inform the user
  print("......have finished loading - now beginning cluster analysis")

  if (file.exists(outfile) == FALSE | repair == 1) {
      nos_uk_clusters = 1 ; uk_cluster = 1 ; uk_cluster_pft = 1
      if (PROJECT$model$name == "DALEC_EVERGREEN" | PROJECT$model$name == "DALEC_CDEA_LU_FIRES" | PROJECT$model$name == "DALEC_CDEA_ACM2" |
          PROJECT$model$name == "DALEC_CDEA_ACM2_BUCKET" | PROJECT$model$name == "DALEC_GSI_DFOL_CWD_FR" |
          PROJECT$model$name == "DALEC_CDEA_ACM2_BUCKET_RmRg" | PROJECT$model$name == "DALEC_CDEA_ACM2_BUCKET_RmRg_CWD" |
          PROJECT$model$name == "DALEC_CDEA_ACM2_BUCKET_RmRg_CWD_wMRT" |
          PROJECT$model$name == "DALEC_GSI_BUCKET" | PROJECT$model$name == "DALEC" |
          PROJECT$model$name == "DALEC_BUCKET" | PROJECT$model$name == "DALEC_BUCKET_CanAGE" |
          PROJECT$model$name == "DALEC_G5" | PROJECT$model$name == "DALEC_G6") {

          # remove non-constrained parameters (i.e. those not actually used in this analysis)
          initial_conditions=c(18:23)
          par_array_median_normalised = grid_parameters$parameters[,,-c(initial_conditions),median_loc]
          if (PROJECT$model$name == "DALEC_GSI_DFOL_CWD_FR") {
              par_array_median_normalised = grid_parameters$parameters[,,-c(initial_conditions,30,31,32,33,37),median_loc]
          } else if (PROJECT$model$name == "DALEC_GSI_BUCKET" | PROJECT$model$name == "DALEC_BUCKET" |
                     PROJECT$model$name == "DALEC" | PROJECT$model$name == "DALEC_BUCKET_CanAGE") {
              par_array_median_normalised = grid_parameters$parameters[,,-c(initial_conditions,30,31,32,33,37),median_loc]
          } else if (PROJECT$model$name == "DALEC_G5" | PROJECT$model$name == "DALEC_G6") {
              par_array_median_normalised = grid_parameters$parameters[,,-c(initial_conditions,37),median_loc]
          } # PROJECT$model$name == "DALEC_GSI_DFOL_CWD_FR"

          # now normalise the dataset
          for (i in seq(1,dim(par_array_median_normalised)[3])) {
               min_par_val=min(grid_parameters$parameters[,,i,median_loc],na.rm=TRUE)
               max_par_val=max(grid_parameters$parameters[,,i,median_loc],na.rm=TRUE)
               par_array_median_normalised[,,i]=((grid_parameters$parameters[,,i,median_loc]-min_par_val)/(max_par_val-min_par_val))
          }
          par_array_tmp=array(NA,dim=c(prod(dim(grid_parameters$parameters)[1:2]),dim(par_array_median_normalised)[3]))
          par_array_tmp[1:prod(dim(par_array_median_normalised)[1:2]),1:dim(par_array_median_normalised)[3]]=par_array_median_normalised
          actual_forests=which(is.na(par_array_tmp[,1]) == FALSE)
          par_array_tmp=par_array_tmp[actual_forests,]
          par_array_tmp=array(par_array_tmp,dim=c((length(par_array_tmp)/dim(par_array_median_normalised)[3]),dim(par_array_median_normalised)[3]))

          tmp = 0
          # Looping to find preference_input, the preferenceRange() returns 2 values, the first of which minimises the number of clusters,
          # while the seconds would return as many clusters as there are observations.
          # It is the responsibility of the user to ensure the most appropriate use of these information to result in an appropriate number of clusters for error propagation
          for (i in seq(1,10)) {
               tmp=append(tmp,preferenceRange(negDistMat(par_array_tmp[sample(1:dim(par_array_tmp)[1],0.05*dim(par_array_tmp)[1], replace=FALSE),],r=2))[1])
          } ; preference_input=max(mean(tmp[-1]),median(tmp[-1]))
          grid_parameters$uk_clusters=apclusterL(negDistMat(r=2),par_array_tmp,frac=0.1,sweeps=10, p=preference_input,maxits=1000, convits=100)
          #grid_parameters$uk_clusters=apclusterL(negDistMat(r=2),par_array_tmp,frac=0.1, sweeps=10, q=0.05, maxits=1000, convits=100)
          grid_parameters$nos_uk_clusters=length(grid_parameters$uk_clusters@clusters) ; uk_clusters_exemplars=grid_parameters$uk_clusters@exemplars
          grid_parameters$uk_cluster_pft=array(NA,dim=c(dim(par_array_median_normalised)[1:2]))
          for (i in seq(1,length(grid_parameters$uk_clusters@clusters))) {
               grid_parameters$uk_cluster_pft[actual_forests[grid_parameters$uk_clusters@clusters[[i]]]] = i
          }
          grid_parameters$uk_cluster_pft=array(grid_parameters$uk_cluster_pft,dim=c(dim(par_array_median_normalised)[1:2]))
      }
  }

  # inform the user
  print("......now generating parameter maps")

  # describe parameter numbers with log liklihood added on the end
  character_bit=rep("p",times=(PROJECT$model$nopars[1]))
  number_bit=1:(PROJECT$model$nopars[1])
  # merge the letter and numbers together
  par_names=c(paste(character_bit,number_bit,sep=""),"log-likelihood")

  # determine correct height and widths
  hist_height=4000 ; hist_width=7200
  fig_height=7000 ; fig_width = ((PROJECT$long_dim/PROJECT$lat_dim)+0.25) * fig_height #7200
  if (PROJECT$grid_type == "UK") { fig_height=8000 ; fig_width=7200 }
  # load colour palette
  colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral")))

  # calculate land mask
  grid_parameters$landmask=array(PROJECT$landsea, dim=c(PROJECT$long_dim,PROJECT$lat_dim))

  # assuming we have generated one lets create the Cluster analysis map
  if (PROJECT$model$name == "DALEC_EVERGREEN" | PROJECT$model$name == "DALEC_CDEA_LU_FIRES" | PROJECT$model$name == "DALEC_CDEA_ACM2" |
      PROJECT$model$name == "DALEC_CDEA_ACM2_BUCKET" | PROJECT$model$name == "DALEC_CDEA_ACM2_BUCKET_RmRg" |
      PROJECT$model$name == "DALEC_CDEA_ACM2_BUCKET_RmRg_CWD" | PROJECT$model$name == "DALEC_CDEA_ACM2_BUCKET_RmRg_CWD_wMRT" |
      PROJECT$model$name == "DALEC_GSI_DFOL_CWD_FR" | PROJECT$model$name == "DALEC_GSI_BUCKET" |
      PROJECT$model$name == "DALEC" | PROJECT$model$name == "DALEC_BUCKET" | PROJECT$model$name == "DALEC_BUCKET_CanAGE" |
      PROJECT$model$name == "DALEC_G5" | PROJECT$model$name == "DALEC_G6") {
      jpeg(file=paste(PROJECT$figpath,"Cluster_map_of_median_parameters_",PROJECT$name,".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
      image.plot(grid_parameters$uk_cluster_pft, main=paste("Cluster analysis potential PFT map",sep=""),axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
      contour(grid_parameters$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
      dev.off()
  }

  # create tempory array of parameters for the covariance analysis
  jpeg(file=paste(PROJECT$figpath,"parameter_correlations_median_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
  tmp=array(grid_parameters$parameters[,,1:dim(grid_parameters$parameters)[3],median_loc],dim=c(prod(dim(grid_parameters$parameters)[1:2]),(dim(grid_parameters$parameters)[3])))
  image.plot(cor(tmp,use="na.or.complete", method="spearman"), main="Median Parameter Correlation (Spearmans)")
  dev.off()

  ###
  ## Some things that are really parameters but are dependent on the evolution of state variables so need more information to be calculated
  # assign this value once as aall arrays have the same number of values present
  colour_choices=colour_choices_upper(length(grid_parameters$NPP_foliar_fraction[,,median_loc]))

  # create map of NPP allocations
  jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_","Cfoliar_NPP_alloc","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
  image.plot(grid_parameters$NPP_foliar_fraction[,,median_loc], col=colour_choices, main=paste("Cfoliar NPP allocation median estimate",sep=""), axes=FALSE
            ,cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
  contour(grid_parameters$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
  dev.off()
  jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_","Cwood_NPP_alloc","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
  image.plot(grid_parameters$NPP_wood_fraction[,,median_loc],col=colour_choices, main=paste("Cwood NPP allocation median estimate",sep=""), axes=FALSE
            ,cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
  contour(grid_parameters$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
  dev.off()
  jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_","Croot_NPP_alloc","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
  image.plot(grid_parameters$NPP_root_fraction[,,median_loc], col=colour_choices, main=paste("Croot NPP allocation median estimate",sep=""), axes=FALSE
            ,cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
  contour(grid_parameters$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
  dev.off()
  jpeg(file=paste(PROJECT$figpath,"parameter_maps_95CI_","Cfoliar_NPP_alloc","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
  image.plot(grid_parameters$NPP_foliar_fraction[,,upper_loc]-grid_parameters$NPP_foliar_fraction[,,lower_loc]
            ,col=colour_choices, main=paste("Cfoliar NPP allocation uncertainty",sep=""),axes=FALSE
            ,cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
  contour(grid_parameters$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
  dev.off()
  jpeg(file=paste(PROJECT$figpath,"parameter_maps_95CI_","Cwood_NPP_alloc","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
  image.plot(grid_parameters$NPP_wood_fraction[,,upper_loc]-grid_parameters$NPP_wood_fraction[,,lower_loc]
            ,col=colour_choices, main=paste("Cwood NPP allocation uncertainty",sep=""),axes=FALSE
            ,cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
  contour(grid_parameters$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
  dev.off()
  jpeg(file=paste(PROJECT$figpath,"parameter_maps_95CI_","Croot_NPP_alloc","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
  image.plot(grid_parameters$NPP_root_fraction[,,upper_loc]-grid_parameters$NPP_root_fraction[,,lower_loc]
            ,col=colour_choices, main=paste("Croot NPP allocation uncertainty",sep=""),axes=FALSE
            ,cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
  contour(grid_parameters$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
  dev.off()

  # create map of residence times
  jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_","Cfoliar_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
  image.plot(grid_parameters$MTT_foliar_years[,,median_loc], col=colour_choices, main=paste("Cfoliar residence time median estimate",sep="")
            ,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
  contour(grid_parameters$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
  dev.off()
  jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_","Cwood_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
  image.plot(grid_parameters$MTT_wood_years[,,median_loc], col=colour_choices, main=paste("Cwood residence time median estimate",sep="")
            ,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
  contour(grid_parameters$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
  dev.off()
  jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_","Croot_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
  image.plot(grid_parameters$MTT_root_years[,,median_loc], col=colour_choices, main=paste("Croot residence time median estimate",sep="")
            ,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
  contour(grid_parameters$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
  dev.off()
  jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_","Csom_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
  image.plot(grid_parameters$MTT_som_years[,,median_loc], col=colour_choices, main=paste("Csom residence time median estimate",sep="")
            ,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
  contour(grid_parameters$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
  dev.off()
  jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_","CDeadOrg_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
  image.plot(grid_parameters$MTT_DeadOrg_years[,,median_loc], col=colour_choices, main=paste("CDeadOrg residence time median estimate",sep="")
            ,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
  contour(grid_parameters$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
  dev.off()

  jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_uncertainty_","Cfoliar_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
  image.plot(grid_parameters$MTT_foliar_years[,,upper_loc]-grid_parameters$MTT_foliar_years[,,lower_loc]
            ,col=colour_choices, main=paste("Cfoliar residence time uncertainty",sep="")
            ,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
  contour(grid_parameters$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
  dev.off()
  jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_uncertainty_","Cwood_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
  image.plot(grid_parameters$MTT_wood_years[,,upper_loc]-grid_parameters$MTT_wood_years[,,lower_loc]
            ,col=colour_choices, main=paste("Cwood residence time uncertainty",sep="")
            ,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
  contour(grid_parameters$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
  dev.off()
  jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_uncertainty_","Croot_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
  image.plot(grid_parameters$MTT_root_years[,,upper_loc]-grid_parameters$MTT_root_years[,,lower_loc]
            ,col=colour_choices, main=paste("Croot residence time uncertainty",sep="")
            ,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
  contour(grid_parameters$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
  dev.off()
  jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_uncertainty_","CDeadOrg_residence_time","_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
  image.plot(grid_parameters$MTT_DeadOrg_years[,,upper_loc]-grid_parameters$MTT_DeadOrg_years[,,lower_loc]
            ,col=colour_choices, main=paste("CDeadOrg residence time uncertainty",sep="")
            ,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
  contour(grid_parameters$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
  dev.off()

  # now everything has been loaded into a nice array we begin plotting
  for (p in seq(1,dim(grid_parameters$parameters)[3])) {
       ###
       ## spatial maps
       jpeg(file=paste(PROJECT$figpath,"parameter_maps_median_",par_names[p],"_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
       par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
       image.plot(grid_parameters$parameters[,,p,median_loc], col=colour_choices, main=paste(par_names[p]," median estimate",sep=""),axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
       contour(grid_parameters$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
       dev.off()
       jpeg(file=paste(PROJECT$figpath,"parameter_maps_95CI_",par_names[p],"_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
       par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
       image.plot(grid_parameters$parameters[,,p,upper_loc]-grid_parameters$parameters[,,p,lower_loc], col=colour_choices, main=paste(par_names[p]," uncertainty range",sep=""),axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
       contour(grid_parameters$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
       dev.off()
       jpeg(file=paste(PROJECT$figpath,"parameter_maps_converged_",par_names[p],"_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
       par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
       image.plot(grid_parameters$parameters_converged[,,p], col=colour_choices, main=paste(par_names[p]," Gelmen-Rubens convergence (1 = PASS / 0 = FALSE)",sep="")
                 ,axes=FALSE, zlim=c(0,1), cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
       contour(grid_parameters$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
       dev.off()
  }

  # Plant function type map
  if (length(which(is.na((as.vector(grid_parameters$pft_array))) == FALSE & as.vector(grid_parameters$pft_array) > 0)) > 0) {
      jpeg(file=paste(PROJECT$figpath,"pft_maps_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
      image.plot(grid_parameters$pft_array, main="ECMWF PFT"
                ,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5
                ,axis.args=list(cex.axis=1.8,hadj=0.1),zlim=c(0,20))
      contour(grid_parameters$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
      dev.off()
  }
  # mean temperature
  jpeg(file=paste(PROJECT$figpath,"mean_temperature_maps_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
  image.plot(grid_parameters$mean_temperature_C, col=rev(colour_choices), main="Mean air temperature (oC)",axes=FALSE, cex.main=2.4,legend.width=3.0
            ,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1)
            ,zlim=c(min(as.vector(grid_parameters$mean_temperature_C),na.rm=TRUE),max(as.vector(grid_parameters$mean_temperature_C),na.rm=TRUE)))
  contour(grid_parameters$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
  dev.off()
  # mean radiation
  jpeg(file=paste(PROJECT$figpath,"mean_radiation_maps_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
  image.plot(grid_parameters$mean_radiation_MJm2day, col=rev(colour_choices), main="Mean radiation (MJ/m2/day)"
            ,axes=FALSE, cex.main=2.4,legend.width=3.0
            ,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1),zlim=c(0,36))
  contour(grid_parameters$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
  dev.off()
  # mean vpd
  jpeg(file=paste(PROJECT$figpath,"mean_vpd_maps_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
  image.plot(grid_parameters$mean_vpd_Pa, col=rev(colour_choices), main="Mean VPD (Pa)"
            ,axes=FALSE, cex.main=2.4, legend.width=3.0, cex=1.5
            ,axis.args=list(cex.axis=1.8,hadj=0.1),zlim=c(0,max(as.vector(grid_parameters$mean_vpd_Pa),na.rm=TRUE)))
  contour(grid_parameters$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
  dev.off()
  # mean precipitation
  jpeg(file=paste(PROJECT$figpath,"mean_precipitation_maps_",PROJECT$name,".jpg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
  image.plot(grid_parameters$mean_precipitation_kgm2yr, col=colour_choices, main="Mean precipitation (mm/yr)",axes=FALSE, cex.main=2.4,legend.width=3.0
            ,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1),zlim=c(0,max(as.vector(grid_parameters$mean_precipitation_kgm2yr),na.rm=TRUE)))
  contour(grid_parameters$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
  dev.off()

  # Add readme information
  grid_parameters$readme = data.frame(Note_1 = "Fire MTT is calculated assuming hardcoded combustion completeness values which may from those used in a specific version of DALEC",
                                      Note_2 = "Biomass removal / harvest MTT calculation is calculated assuming a proportional amount of fine root dies, this is not the assumption used in CDEA model versions.")
  # output some aggragated values
  # probably best to add some aggregated met drivers to this concoction here
  save(grid_parameters,file=outfile, compress = "gzip")

  # tidy before leaving
  gc(reset=TRUE, verbose=FALSE)

  print("...done with parameters...")

  } # end function generate_parameter_maps

  ## Use byte compile
  generate_parameter_maps<-cmpfun(generate_parameter_maps)
