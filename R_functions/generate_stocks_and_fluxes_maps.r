
###
## Function to create gridded, uncertainty propogated plots for the gridded CARDAMOM analyses.
## Includes functions to determine the linear components of GSI's associations with DALEC simulated LAI.
## NOTE: these functions are largely superceded by those found in "generate_simplified_stock_and_flux_maps.r"
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).
# The GSI concept is derived from the Growing Season Index described by Jolley et al., (2005)

gsi_controlling<-function(gsi_components,gsi,gsi_wanted,timesteps,m) {

  ###
  ## function to determine the fraction of time a given GSI component is providing the dominent control on the gradient for GSI control

  # gsi_wanted = scaler i.e. which temp, vpd, photoperiod to assess
  # gsi_components[iter,time,component] all the components to determine ht minimum
  # gsi[iter,time] = the actual gsi estimate

  # check dimensions
  if (gsi_wanted > dim(gsi_components)[3]) {stop('Desired GSI component is greater than available drivers')}

  # calculate gradients over averaging period
  x_all=1:22
  gsi_gradient=array(0,dim=dim(gsi)) ; components_gradient=array(0,dim=c(dim(gsi),dim(gsi_components)[3]))
  # loop through time steps
  for (j in seq(2,length(m))) {
    # calculate time step averaging period
    adjust_to_daily = round(sum(timesteps[(j-m[j]+1):j]/m[j]),digits=0)
    # calculate gradient of the overall GSI variable
    x = x_all[1:(m[j]+1)] ; interval=length(x_all[1:(m[j]+1)])
    # calculate the sum of x ; # calculate the sum of y
    sum_x = sum(x) ; sum_y = rowSums(gsi[,(j-m[j]):j])
    # calculate the sum of squares of x ; # calculate the sum of the product of xy
    sumsq_x = sum(x*x)  ; sum_product_xy = rowSums(t(((x)*t(gsi[,(j-m[j]):j]))))
    # calculate the gradient
    gsi_gradient[,j] = ( (interval*sum_product_xy) - (sum_x*sum_y) ) / ( (interval*sumsq_x) - (sum_x*sum_x) )
    gsi_gradient[,j] = gsi_gradient[,j] / adjust_to_daily
    # now loop through and calculate the gradients for each GSI component
    for (comp in seq(1, dim(gsi_components)[3])) {
      # calculate the sum of x ; # calculate the sum of y
      sum_y = rowSums(gsi_components[1:dim(gsi)[1],(j-m[j]):j,comp])
      # calculate the sum of squares of x ; # calculate the sum of the product of xy
      sum_product_xy = rowSums(t((x)*t(gsi_components[1:dim(gsi)[1],(j-m[j]):j,comp])))
      # calculate the gradient
      components_gradient[,j,comp] = ( (interval*sum_product_xy) - (sum_x*sum_y) ) / ( (interval*sumsq_x) - (sum_x*sum_x) )
      components_gradient[,j,comp] = components_gradient[,j,comp] / adjust_to_daily
    }
  }

  # declare output variable
  tmp=array(0,dim=c(dim(gsi)[1],dim(components_gradient)[3]))
  output=array(0,dim=c(dim(gsi)[1],dim(components_gradient)[3]))
  # loop through components to determine solution
  for (t in seq(1, dim(gsi_gradient)[1])) {
    vars = 1:3 ; vars=vars[which(apply(components_gradient[t,,],2,max) != 0)]
    aa = lm(gsi_gradient[t,]~components_gradient[t,,1]+components_gradient[t,,2]+components_gradient[t,,3])
    aa = (anova(aa)$"Sum Sq")[1:length(vars)]
    tmp[t,vars]=aa/sum(aa)
  }

  # calculate % time control is by the chosen variable
  output=(tmp/rowSums(tmp))*100

  # back to user
  return(output)

} # end function gsi_controlling
# Use byte compile
gsi_controlling<-cmpfun(gsi_controlling)

gsi_limiting<-function(gsi_components,gsi_wanted) {

  ###
  ## function to determine the fraction of time a given GSI component is limiting

  # gsi_wanted = scaler i.e. which temp, vpd, photoperiod to assess
  # gsi_components[iter,time,component] all the components to determine the minimum

  # check dimensions
  if (gsi_wanted > dim(gsi_components)[3]) {stop('Desired GSI component is greater than available drivers')}

  # determine which of the GSI components are the ones I am comparing against
  to_do=1:dim(gsi_components)[3] ; to_do = to_do[which(to_do != gsi_wanted)]

  # declare output variable
  output=array(FALSE,dim=c(dim(gsi_components)[1:2]))

  # loop through components to determine solution
  for (t in seq(1,dim(gsi_components)[2])) {
    # true / false and I most limiting
    output[,t]=gsi_components[,t,gsi_wanted] < pmin(gsi_components[,t,to_do[1]],gsi_components[,t,to_do[2]])
  }

  # tidy memory
  gc()

  # back to user
  return(output)

} #  end function gsi_limiting
## Use byte compile
gsi_limiting<-cmpfun(gsi_limiting)

run_mcmc_results_for_grid <-function (n,PROJECT,parameters,sub_parameter,drivers) {

  # To avoid additional io requirement when running the grid analysis we will not save the states and fluxes to a file.
  # Insteady we generate them as an when into local memory

  # run subsample of parameters for full results / propogation
  soil_info=c(drivers$top_sand,drivers$bot_sand,drivers$top_clay,drivers$bot_clay)
  states_all=simulate_all(n,PROJECT,PROJECT$model$name,drivers$met,sub_parameter[1:PROJECT$model$nopars[n],,],drivers$lat,PROJECT$ctessel_pft[n],PROJECT$parameter_type,PROJECT$exepath,soil_info)

  # back to user
  return(states_all)

} # end function run_mcmc_results_for_grid
## Use byte compile
run_mcmc_results_for_grid<-cmpfun(run_mcmc_results_for_grid)

parallel_cluster_processing<-function(c,analysis_info,uk_cluster_pft,
                                      outfile_tmp,area,area_with_g_Tg,
                                      total_area,total_area_1,PROJECT) {

    # Function is called to loop through an entire cluster

    # determine whcih sites are part of the cluster called for
    sites_to_do = analysis_info$site_locations[which(analysis_info$site_cluster == c)]
    print(paste("......there are ",length(sites_to_do)," sites in cluster ",c,sep=""))

    # now begin looping through the sites in this cluster
    for (nn in seq(1,length(sites_to_do))) {

         print(paste(".........Cluster: ",c," beginning site ",nn," of ",length(sites_to_do),sep=""))

         # assign correct site ID
         n = sites_to_do[nn]

         # calculate pixel location information
         slot_j=as.numeric(PROJECT$sites[n])/PROJECT$long_dim
         slot_i=as.numeric(PROJECT$sites[n])-(floor(slot_j)*PROJECT$long_dim)
         if(slot_i == 0) {slot_i = PROJECT$long_dim} ; slot_j=ceiling(slot_j)

         # generate file name of the output file created in stage 3
         loadfile = paste(PROJECT$results_processedpath,PROJECT$sites[n],"_parameters.RData",sep="")
         #print(paste(".........processing site:",PROJECT$sites[n],sep=""))
         # load driver / parameter information
         load(loadfile)
         # run the parameter chains to generate states but without actually saving them this time
         states_all = run_mcmc_results_for_grid(n,PROJECT,parameters,sub_parameter,drivers)
         # remove variable that we don't need here
         states_all = states_all[which(names(states_all) != "aNPP")] ; gc()
         pixel_dims = dim(states_all[[1]]) ; current_pars = pixel_dims[1]

         if (nn == 1) {
             # define the local variables needed for output
             states_array_median = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,analysis_info$pixel_vars+analysis_info$nos_extra))
             states_array_unc = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,analysis_info$pixel_vars+analysis_info$nos_extra))
             gsi_control_median = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,3))
             gsi_control_unc = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,3))
             final_stocks_median = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,analysis_info$no_pools+1))
             final_stocks_uncertainty = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,analysis_info$no_pools+1))
             stock_change_rate_median = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,analysis_info$no_pools))
             stock_change_rate_uncertainty = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,analysis_info$no_pools))
             pixel_management_array = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
             cumulative_management_array = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
             # cluster / timeseries related
             states_mean_timeseries = array(0, dim=c(analysis_info$max_pars,pixel_dims[2],analysis_info$pixel_vars))
             states_sum_timeseries = array(0, dim=c(analysis_info$max_pars,pixel_dims[2],analysis_info$pixel_vars))
             biomass_change_sum_timeseries = array(0, dim=c(analysis_info$max_pars,pixel_dims[2],length(analysis_info$biomass_change_names)))
             biomass_change_mean_timeseries = array(0, dim=c(analysis_info$max_pars,pixel_dims[2],length(analysis_info$biomass_change_names)))
         }
         if (max(drivers$met[,8]) > 0) {
            pixel_management_array[slot_i,slot_j] = ceiling(max(which(drivers$met[,8] > 0))/analysis_info$steps_per_year)
            cumulative_management_array[slot_i,slot_j] = sum(drivers$met[which(drivers$met[,8] > 0),8])
         }
         # calculate the final stock medians and uncertainty
         bob=quantile(states_all$bio_gCm2[,pixel_dims[2]], prob=c(0.975,0.5,0.025))
         final_stocks_median[slot_i,slot_j,1]=bob[2] ; final_stocks_uncertainty[slot_i,slot_j,1]=abs(bob[1]-bob[3])
         bob=quantile(states_all$lab_gCm2[,pixel_dims[2]], prob=c(0.975,0.5,0.025))
         final_stocks_median[slot_i,slot_j,2]=bob[2] ; final_stocks_uncertainty[slot_i,slot_j,2]=abs(bob[1]-bob[3])
         bob=quantile(states_all$fol_gCm2[,pixel_dims[2]], prob=c(0.975,0.5,0.025))
         final_stocks_median[slot_i,slot_j,3]=bob[2] ; final_stocks_uncertainty[slot_i,slot_j,3]=abs(bob[1]-bob[3])
         bob=quantile(states_all$root_gCm2[,pixel_dims[2]], prob=c(0.975,0.5,0.025))
         final_stocks_median[slot_i,slot_j,4]=bob[2] ; final_stocks_uncertainty[slot_i,slot_j,4]=abs(bob[1]-bob[3])
         bob=quantile(states_all$wood_gCm2[,pixel_dims[2]], prob=c(0.975,0.5,0.025))
         final_stocks_median[slot_i,slot_j,5]=bob[2] ; final_stocks_uncertainty[slot_i,slot_j,5]=abs(bob[1]-bob[3])
         bob=quantile(states_all$lit_gCm2[,pixel_dims[2]], prob=c(0.975,0.5,0.025))
         final_stocks_median[slot_i,slot_j,6]=bob[2] ; final_stocks_uncertainty[slot_i,slot_j,6]=abs(bob[1]-bob[3])
         bob=quantile(states_all$som_gCm2[,pixel_dims[2]], prob=c(0.975,0.5,0.025))
         final_stocks_median[slot_i,slot_j,7]=bob[2] ; final_stocks_uncertainty[slot_i,slot_j,7]=abs(bob[1]-bob[3])
         if (length(which(names(states_all) == "litwood_gCm2")) > 0) {
             bob=quantile(states_all$litwood_gCm2[,pixel_dims[2]], prob=c(0.975,0.5,0.025))
             final_stocks_median[slot_i,slot_j,8]=bob[2] ; final_stocks_uncertainty[slot_i,slot_j,8]=abs(bob[1]-bob[3])
         }

         # calculate stock change over time without accounting for clearances
         bob = quantile((states_all$lab_gCm2)[,pixel_dims[2]]-(states_all$lab_gCm2)[,1], prob=c(0.975,0.50,0.025))
         stock_change_rate_median[slot_i,slot_j,1] = bob[2]/analysis_info$nos_years
         stock_change_rate_uncertainty[slot_i,slot_j,1] = abs(bob[3]-bob[1])/analysis_info$nos_years
         bob = quantile((states_all$fol_gCm2)[,pixel_dims[2]]-(states_all$fol_gCm2)[,1], prob=c(0.975,0.50,0.025))
         stock_change_rate_median[slot_i,slot_j,2] = bob[2]/analysis_info$nos_years
         stock_change_rate_uncertainty[slot_i,slot_j,2] = abs(bob[3]-bob[1])/analysis_info$nos_years
         bob = quantile((states_all$root_gCm2)[,pixel_dims[2]]-(states_all$root_gCm2)[,1], prob=c(0.975,0.50,0.025))
         stock_change_rate_median[slot_i,slot_j,3] = bob[2]/analysis_info$nos_years
         stock_change_rate_uncertainty[slot_i,slot_j,3] = abs(bob[3]-bob[1])/analysis_info$nos_years
         bob = quantile((states_all$wood_gCm2)[,pixel_dims[2]]-(states_all$wood_gCm2)[,1], prob=c(0.975,0.50,0.025))
         stock_change_rate_median[slot_i,slot_j,4] = bob[2]/analysis_info$nos_years
         stock_change_rate_uncertainty[slot_i,slot_j,4] = abs(bob[3]-bob[1])/analysis_info$nos_years
         bob = quantile((states_all$lit_gCm2)[,pixel_dims[2]]-(states_all$lit_gCm2)[,1], prob=c(0.975,0.50,0.025))
         stock_change_rate_median[slot_i,slot_j,5] = bob[2]/analysis_info$nos_years
         stock_change_rate_uncertainty[slot_i,slot_j,5] = abs(bob[3]-bob[1])/analysis_info$nos_years
         bob = quantile((states_all$som_gCm2)[,pixel_dims[2]]-(states_all$som_gCm2)[,1], prob=c(0.975,0.50,0.025))
         stock_change_rate_median[slot_i,slot_j,6] = bob[2]/analysis_info$nos_years
         stock_change_rate_uncertainty[slot_i,slot_j,6] = abs(bob[3]-bob[1])/analysis_info$nos_years
         if (length(which(names(states_all) == "litwood_gCm2")) > 0) {
             bob = quantile((states_all$litwood_gCm2)[,pixel_dims[2]]-(states_all$litwood_gCm2)[,1], prob=c(0.975,0.50,0.025))
             stock_change_rate_median[slot_i,slot_j,7] = bob[2]/analysis_info$nos_years
             stock_change_rate_uncertainty[slot_i,slot_j,7] = abs(bob[3]-bob[1])/analysis_info$nos_years
         }

         # move on to biomass related / spatial maps of states and fluxes, these will have a cluster factor associated with them which is reconstructed at a higher level.

         # stock change in time series
         tmp = (t(t(apply((states_all$wood_gCm2-states_all$wood_gCm2[,1]),2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))))*area_with_g_Tg[slot_i,slot_j])
         biomass_change_sum_timeseries[1:analysis_info$max_pars,,1]=biomass_change_sum_timeseries[1:analysis_info$max_pars,,1]+tmp
         tmp=apply((states_all$wood_gCm2-states_all$wood_gCm2[,1]),2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))*(area[slot_i,slot_j]*total_area_1)
         biomass_change_mean_timeseries[1:analysis_info$max_pars,,1]=biomass_change_mean_timeseries[1:analysis_info$max_pars,,1]+tmp
         tmp=(t(t(apply((states_all$som_gCm2-states_all$som_gCm2[,1]),2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))))*area_with_g_Tg[slot_i,slot_j])
         biomass_change_sum_timeseries[1:analysis_info$max_pars,,2]=biomass_change_sum_timeseries[1:analysis_info$max_pars,,2]+tmp
         tmp=apply((states_all$som_gCm2-states_all$som_gCm2[,1]),2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))*(area[slot_i,slot_j]*total_area_1)
         biomass_change_mean_timeseries[1:analysis_info$max_pars,,2]=biomass_change_mean_timeseries[1:analysis_info$max_pars,,2]+tmp
         tmp=(t(t(apply((states_all$bio_gCm2-states_all$bio_gCm2[,1]),2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))))*area_with_g_Tg[slot_i,slot_j])
         biomass_change_sum_timeseries[1:analysis_info$max_pars,,3]=biomass_change_sum_timeseries[1:analysis_info$max_pars,,3]+tmp
         tmp=apply((states_all$bio_gCm2-states_all$bio_gCm2[,1]),2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))*(area[slot_i,slot_j]*total_area_1)
         biomass_change_mean_timeseries[1:analysis_info$max_pars,,3]=biomass_change_mean_timeseries[1:analysis_info$max_pars,,3]+tmp
         tmp=(t(t(apply((states_all$root_gCm2-states_all$root_gCm2[,1]),2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))))*area_with_g_Tg[slot_i,slot_j])
         biomass_change_sum_timeseries[1:analysis_info$max_pars,,4]=biomass_change_sum_timeseries[1:analysis_info$max_pars,,4]+tmp
         tmp=apply((states_all$root_gCm2-states_all$root_gCm2[,1]),2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))*(area[slot_i,slot_j]*total_area_1)
         biomass_change_mean_timeseries[1:analysis_info$max_pars,,4]=biomass_change_mean_timeseries[1:analysis_info$max_pars,,4]+tmp
         tmp=(t(t(apply((states_all$lit_gCm2-states_all$lit_gCm2[,1]),2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))))*area_with_g_Tg[slot_i,slot_j])
         biomass_change_sum_timeseries[1:analysis_info$max_pars,,5]=biomass_change_sum_timeseries[1:analysis_info$max_pars,,5]+tmp
         tmp=apply((states_all$lit_gCm2-states_all$lit_gCm2[,1]),2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))*(area[slot_i,slot_j]*total_area_1)
         biomass_change_mean_timeseries[1:analysis_info$max_pars,,5]=biomass_change_mean_timeseries[1:analysis_info$max_pars,,5]+tmp
         tmp=(t(t(apply((states_all$lab_gCm2-states_all$lab_gCm2[,1]),2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))))*area_with_g_Tg[slot_i,slot_j])
         biomass_change_sum_timeseries[1:analysis_info$max_pars,,6]=biomass_change_sum_timeseries[1:analysis_info$max_pars,,6]+tmp
         tmp=apply((states_all$lab_gCm2-states_all$lab_gCm2[,1]),2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))*(area[slot_i,slot_j]*total_area_1)
         biomass_change_mean_timeseries[1:analysis_info$max_pars,,6]=biomass_change_mean_timeseries[1:analysis_info$max_pars,,6]+tmp
         tmp=(t(t(apply((states_all$fol_gCm2-states_all$fol_gCm2[,1]),2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))))*area_with_g_Tg[slot_i,slot_j])
         biomass_change_sum_timeseries[1:analysis_info$max_pars,,7]=biomass_change_sum_timeseries[1:analysis_info$max_pars,,7]+tmp
         tmp=apply((states_all$fol_gCm2states_all$fol_gCm2[,1]),2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))*(area[slot_i,slot_j]*total_area_1)
         biomass_change_mean_timeseries[1:analysis_info$max_pars,,7]=biomass_change_mean_timeseries[1:analysis_info$max_pars,,7]+tmp
         # some optional values
         if (length(which(names(states_all) == "litwood_gCm2")) > 0) {
             tmp=(t(t(apply((states_all$litwood_gCm2-states_all$litwood_gCm2[,1]),2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))))*area_with_g_Tg[slot_i,slot_j])
             biomass_change_sum_timeseries[1:analysis_info$max_pars,,8]=biomass_change_sum_timeseries[1:analysis_info$max_pars,,8]+tmp
             tmp=apply((states_all$litwood_gCm2-states_all$litwood_gCm2[,1]),2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))*(area[slot_i,slot_j]*total_area_1)
             biomass_change_mean_timeseries[1:analysis_info$max_pars,,8]=biomass_change_mean_timeseries[1:analysis_info$max_pars,,8]+tmp
         } else if (length(which(names(states_all) == "sfc_water_mm")) > 0) {
             tmp=(t(t(apply((states_all$sfc_water_mm-states_all$sfc_water_mm[,1]),2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))))*area[slot_i,slot_j])
             biomass_change_sum_timeseries[1:analysis_info$max_pars,,9]=biomass_change_sum_timeseries[1:analysis_info$max_pars,,9]+tmp
             tmp=apply((states_all$sfc_water_mm-states_all$sfc_water_mm[,1]),2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))*(area[slot_i,slot_j]*total_area_1)
             biomass_change_mean_timeseries[1:analysis_info$max_pars,,9]=biomass_change_mean_timeseries[1:analysis_info$max_pars,,9]+tmp
         }

         # loop through all processes states for generic output
         for (p in seq(1, length(analysis_info$in_file_names))) {
              infile_location=which(names(states_all) == analysis_info$in_file_names[p])
              if (length(infile_location) > 0) {
                  if (analysis_info$in_file_names[p] == "gsi") {
                      # calculate percentage of time that a given GSI parameter smallest
                      gsi_counter = array(0,dim=c(current_pars,3)) ; gsi_input = array(0,dim=c(current_pars,pixel_dims[2],3))
                      gsi_input[1:pixel_dims[1],,1] = states_all$gsi_itemp ; gsi_input[1:pixel_dims[1],,2]=states_all$gsi_iphoto
                      gsi_input[1:pixel_dims[1],,3] = states_all$gsi_ivpd
                      # extract mean temperature components
                      bob = gsi_limiting(gsi_input,gsi_wanted=1)
                      gsi_counter[,1] = (apply(bob,1,sum,na.rm=TRUE)/pixel_dims[2])*1e2
                      bob = quantile(gsi_counter[,1],prob=c(0.975,0.50,0.025),na.rm=TRUE)
                      states_array_median[slot_i,slot_j,which(analysis_info$par_names == "GSI_itemp")] = bob[2]
                      states_array_unc[slot_i,slot_j,which(analysis_info$par_names == "GSI_itemp")] = abs(bob[3]-bob[1])
                      # extract mean photoperiod component
                      bob = gsi_limiting(gsi_input,gsi_wanted=2)
                      gsi_counter[,2] = (apply(bob,1,sum,na.rm=TRUE)/pixel_dims[2])*1e2
                      bob = quantile(gsi_counter[,2],prob=c(0.975,0.50,0.025),na.rm=TRUE)
                      states_array_median[slot_i,slot_j,which(analysis_info$par_names == "GSI_iphoto")] = bob[2]
                      states_array_unc[slot_i,slot_j,which(analysis_info$par_names == "GSI_iphoto")] = abs(bob[3]-bob[1])
                      # extract mean VPD  / soil moisture component
                      bob = gsi_limiting(gsi_input,gsi_wanted=3)
                      gsi_counter[,3] = (apply(bob,1,sum,na.rm=TRUE)/pixel_dims[2])*1e2
                      bob = quantile(gsi_counter[,3],prob=c(0.975,0.50,0.025),na.rm=TRUE)
                      states_array_median[slot_i,slot_j,which(analysis_info$par_names == "GSI_ivpd")] = bob[2]
                      states_array_unc[slot_i,slot_j,which(analysis_info$par_names == "GSI_ivpd")] = abs(bob[3]-bob[1])
                      # now calculate the time when the GSI value is having the greatest impact on the gradient that determines leaf growth or senescence
                      bob = apply(gsi_controlling(gsi_input,states_all$gsi,1,analysis_info$timestep_days,analysis_info$tmp_m),2,quantile,prob=c(0.975,0.5,0.025),na.rm=TRUE)
                      gsi_control_median[slot_i,slot_j,1:3]=bob[2,]
                      gsi_control_unc[slot_i,slot_j,1:3]=abs(bob[1,]-bob[3,])
                      # now calculate time series information for gsi and each of its components
                      # sample distribution to max_pars samples (mostly this is to save compute time)
                      # GSI
                      tmp = t(t(apply(states_all$gsi,2,quantile,prob = seq(0,1,length.out=analysis_info$max_pars)))*analysis_info$timestep_days)*area[slot_i,slot_j]
                      states_sum_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],p] = states_sum_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],p]+tmp
                      tmp = apply(states_all$gsi,2,quantile,prob = seq(0,1,length.out=analysis_info$max_pars))*(area[slot_i,slot_j]*total_area_1)
                      states_mean_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],p] = states_mean_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],p] + tmp
                      # gsi_itemp
                      tmp = t(t(apply(states_all$gsi_itemp,2,quantile,prob = seq(0,1,length.out=analysis_info$max_pars))))*area[slot_i,slot_j]
                      states_sum_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],which(analysis_info$in_file_names == "gsi_itemp")] = states_sum_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],which(analysis_info$in_file_names == "gsi_itemp")]+tmp
                      tmp = apply(states_all$gsi_itemp,2,quantile,prob = seq(0,1,length.out=analysis_info$max_pars))*(area[slot_i,slot_j]*total_area_1)
                      states_mean_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],which(analysis_info$in_file_names == "gsi_itemp")] = states_mean_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],which(analysis_info$in_file_names == "gsi_itemp")]+tmp
                      # gsi_iphoto
                      tmp=t(t(apply(states_all$gsi_iphoto,2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))))*area[slot_i,slot_j]
                      states_sum_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],which(analysis_info$in_file_names == "gsi_iphoto")] = states_sum_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],which(analysis_info$in_file_names == "gsi_iphoto")]+tmp
                      tmp=apply(states_all$gsi_iphoto,2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))*(area[slot_i,slot_j]*total_area_1)
                      states_mean_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],which(analysis_info$in_file_names == "gsi_iphoto")] = states_mean_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],which(analysis_info$in_file_names == "gsi_iphoto")]+tmp
                      # gsi_ivpd / soil moisture
                      tmp=t(t(apply(states_all$gsi_ivpd,2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))))*area[slot_i,slot_j]
                      states_sum_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],which(analysis_info$in_file_names == "gsi_ivpd")] = states_sum_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],which(analysis_info$in_file_names == "gsi_ivpd")]+tmp
                      tmp=apply(states_all$gsi_ivpd,2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))*(area[slot_i,slot_j]*total_area_1)
                      states_mean_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],which(analysis_info$in_file_names == "gsi_ivpd")] = states_mean_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],which(analysis_info$in_file_names == "gsi_ivpd")]+tmp
                      rm(gsi_counter,gsi_input)
                  } else if (analysis_info$in_file_names[p] == "lai_m2m2" | analysis_info$in_file_names[p] == "wSWP_MPa") {
                      # non-C pools these do not have a unit or time correction applied to them
                      bob = apply(states_all[[infile_location]],2,quantile, prob=c(0.975,0.50,0.025))
                      states_array_median[slot_i,slot_j,p] = mean(bob[2,])
                      states_array_unc[slot_i,slot_j,p] = mean(abs(bob[1,]-bob[3,]))
                      # sample distribution to max_pars samples (mostly this is to save compute time)
                      tmp=t(t(apply(states_all[[infile_location]],2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))))*area[slot_i,slot_j]
                      states_sum_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],p] = states_sum_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],p] + tmp
                      tmp=apply(states_all[[infile_location]],2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))*(area[slot_i,slot_j]*total_area_1)
                      states_mean_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],p] = states_mean_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],p] + tmp
                  } else if (analysis_info$in_file_names[p] == "fire_gCm2day"  | analysis_info$in_file_names[p] == "harvest_C_gCm2day" | analysis_info$in_file_names[p] == "evap_kgH2Om2day" |
                             analysis_info$in_file_names[p] == "gpp_gCm2day"   | analysis_info$in_file_names[p] == "nee_gCm2day"       | analysis_info$in_file_names[p] == "reco_gCm2day" |
                             analysis_info$in_file_names[p] == "rauto_gCm2day" | analysis_info$in_file_names[p] == "rhet_gCm2day"      | analysis_info$in_file_names[p] == "wood_gCm2" |
                             analysis_info$in_file_names[p] == "som_gCm2"      | analysis_info$in_file_names[p] == "bio_gCm2"          | analysis_info$in_file_names[p] == "root_gCm2" |
                             analysis_info$in_file_names[p] == "lit_gCm2"      | analysis_info$in_file_names[p] == "lab_gCm2"          | analysis_info$in_file_names[p] == "fol_gCm2"  |
                             analysis_info$in_file_names[p] == "litwood_gCm2") {
                      # remove time step multiplication from time step but do include unit correction with area
                      bob = apply(states_all[[infile_location]],2,quantile, prob=c(0.975,0.50,0.025))
                      states_array_median[slot_i,slot_j,p] = mean(bob[2,])
                      states_array_unc[slot_i,slot_j,p] = mean(abs(bob[1,]-bob[3,]))
                      # sample distribution to max_pars samples (mostly this is to save compute time)
                      tmp=t(t(apply(states_all[[infile_location]],2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))))*area_with_g_Tg[slot_i,slot_j]
                      states_sum_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],p] = states_sum_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],p]+tmp
                      tmp=apply(states_all[[infile_location]],2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))*(area[slot_i,slot_j]*total_area_1)
                      states_mean_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],p] = states_mean_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],p]+tmp
                  } else if (analysis_info$in_file_names[p] == "sfc_water_mm") {
                      # remove time step multiplication from time step but do include unit correction with area
                      bob = apply(states_all[[infile_location]],2,quantile, prob=c(0.975,0.50,0.025))
                      states_array_median[slot_i,slot_j,p] = mean(bob[2,])
                      states_array_unc[slot_i,slot_j,p] = mean(abs(bob[1,]-bob[3,]))
                      # sample distribution to max_pars samples (mostly this is to save compute time)
                      tmp=t(t(apply(states_all[[infile_location]],2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))))*area[slot_i,slot_j]
                      states_sum_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],p] = states_sum_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],p] + tmp
                      tmp=apply(states_all[[infile_location]],2,quantile,prob=seq(0,1,length.out=analysis_info$max_pars))*(area[slot_i,slot_j]*total_area_1)
                      states_mean_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],p] = states_mean_timeseries[1:analysis_info$max_pars,1:pixel_dims[2],p]+tmp
                  } # variable choice
              } # is the variable present
         } # end of for loop for variables

         # tidy up attempted
         rm(states_all,bob,tmp) ; gc()

    } # end of site loop

    # prepare output to be sent back to the master
    output=list(states_array_median = states_array_median,states_array_unc = states_array_unc,states_sum_timeseries = states_sum_timeseries,
                states_mean_timeseries = states_mean_timeseries,biomass_change_mean_timeseries = biomass_change_mean_timeseries,
                biomass_change_sum_timeseries = biomass_change_sum_timeseries,gsi_control_median = gsi_control_median,gsi_control_unc = gsi_control_unc,
                final_stocks_median = final_stocks_median,final_stocks_uncertainty = final_stocks_uncertainty,
                pixel_management_array = pixel_management_array,cumulative_management_array=cumulative_management_array,
                stock_change_rate_median = stock_change_rate_median,stock_change_rate_uncertainty = stock_change_rate_uncertainty)
    # back top the user
    return(output)

} # parallel_cluster_processing

###
## Function to extract state variables and direct the production of uncertainty plots for the key states and fluxes
###

generate_stocks_and_fluxes_maps<-function(PROJECT) {

    print("...beginning stocks and fluxes...")

    # how many years of the analysis
    nos_years=length(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year))
    # load timesteps to local variable
    timestep_days=PROJECT$model$timestep_days
    if (length(timestep_days) == 1) {
        eh = TRUE ; n = 0
        while(eh) {
          n = n + 1
          # generate file name of the output file created in stage 3
          loadfile=paste(PROJECT$results_processedpath,PROJECT$sites[n],"_parameters.RData",sep="")
          if (file.exists(loadfile)) {load(loadfile) ; eh = FALSE}
        }
        timestep_days=rep(timestep_days, length.out=drivers$nodays)
    }
    # caculate number of model steps in the year
    steps_per_year=floor(365.25/mean(timestep_days))
    # extra variables for flux list
    nos_extra=0#6

    # determine averaging period for timeseries
    tmp_m=rep(0, times=length(timestep_days))
    for (n in seq(1,length(timestep_days))) {
         # calculate the gradient / trend of GSI
         if (sum(timestep_days[1:n]) < 21) {
             tmp_m[n] = n-1
         } else {
             # else we will try and work out the gradient to see what is
             # happening to the system over all. The default assumption will be to
             # consider the averaging period of GSI model (i.e. 21 days). If this is not
             # possible either the time step of the system is used (if step
             # greater
             # than 21 days) or all available steps (if n < 21).
             m = 0 ; test = 0
            while (test < 21) {
              m = m+1 ; test = sum(timestep_days[(n-m):n])
              if (m > (n-1)) {test = 21}
            } # while
            tmp_m[n] = m
         } # for calculating gradient
    } # timestep_days

    # work out area matrix for the pixels in meters
    # include adjustment for g-> Tg (*1e-12)
    if (PROJECT$grid_type == "UK") {
        area_with_g_Tg = array(PROJECT$resolution**2, dim=c(PROJECT$long_dim,PROJECT$lat_dim))*1e-12
        area = array(PROJECT$resolution**2, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
    } else if (PROJECT$grid_type == "wgs84") {
        # generate the lat / long grid again
        output = generate_wgs84_grid(PROJECT$latitude,PROJECT$longitude,PROJECT$resolution)
        # then generate the area estimates for each pixel
        area_with_g_Tg = calc_pixel_area(output$lat,output$long,PROJECT$resolution)*1e-12
        area = calc_pixel_area(output$lat,output$long,PROJECT$resolution)
        # this output is in vector form and we need matching array shapes so...
        area = array(area, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
        area_with_g_Tg = array(area_with_g_Tg, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
    } else {
        stop("valid spatial grid option not selected (UK, or wgs84)")
    }

    # load the parameter file as we may need some of the stored information
    load(paste(PROJECT$results_processedpath,PROJECT$name,"_parameter_maps.RData",sep=""))

    # determine new output file for aggregated values
    outfile=paste(PROJECT$results_processedpath,PROJECT$name,"_flux_maps.RData",sep="")

    # variable names
    par_names=c("LAI","GPP","NEE","Reco","Rauto","Rhet","Cwood","Csom","CTotal","Croot","Clitter","Clabile","Cfoliage","Charvested","GSI","GSI_itemp","GSI_iphoto","GSI_ivpd","Ccwd","Evap","RootWater","wSWP","Fire","GSI_itemp_limiting","GSI_iphoto_limiting","GSI_ivpd_limiting")
    in_file_names=c("lai_m2m2","gpp_gCm2day","nee_gCm2day","reco_gCm2day","rauto_gCm2day","rhet_gCm2day","wood_gCm2","som_gCm2","bio_gCm2","root_gCm2","lit_gCm2","lab_gCm2","fol_gCm2","harvest_C_gCm2day","gsi","gsi_itemp","gsi_iphoto","gsi_ivpd","litwood_gCm2","evap_kgH2Om2day","sfc_water_mm","wSWP_MPa","fire_gCm2day")

    # variable units
    par_names_units=c("m2/m2","gC.m-2.day-1","gC.m-2.day-1","gC.m-2.day-1","gC.m-2.day-1","gC.m-2.day-1","gC.m-2","gC.m-2","gC.m-2","gC.m-2","gC.m-2","gC.m-2","gC.m-2","gC.m-2","0-1","0-1","0-1","0-1","gC.m-2","kgH2O.m-2.day-1","mm","MPa","gC.m-2.day-1","%","%","%")
    # for area sum values
    par_names_units_sum=c("m2/m2","TgC.day-1","TgC.day-1","TgC.day-1","TgC.day-1","TgC.day-1","TgC","TgC","TgC","TgC","TgC","TgC","TgC","TgC","0-1","0-1","0-1","0-1","TgC","TgH2O","mm","MPa","TgC.day-1")
    # final stock / stock change names
    stock_change_names=c("Clabile","Cfoliage","Croots","Cwood","Clitter","Csom","Ccwd","RootWater")
    stock_change_names_mean=c("gC.m-2","gC.m-2","gC.m-2","gC.m-2","gC.m-2","gC.m-2","gC.m-2","mm")
    stock_change_names_sum=c("TgC","TgC","TgC","TgC","TgC","TgC","TgC")
    # final stock / stock change names
    final_stock_names=c("Clabile","Cfoliage","Croots","Cwood","Clitter","Csom","Ccwd","RootWater")
    final_stock_names_mean=c("gC.m-2","gC.m-2","gC.m-2","gC.m-2","gC.m-2","gC.m-2","gC.m-2","mm")
    final_stock_names_sum=c("TgC","TgC","TgC","TgC","TgC","TgC","TgC")
    # biomass change names
    biomass_change_names=c("Cwood","Csom","CTotal","Croot","Clitter","Clabile","Cfoliage","Ccwd","RootWater")
    biomass_change_names_mean=c("gC.m-2","gC.m-2","gC.m-2","gC.m-2","gC.m-2","gC.m-2","gC.m-2","gC.m-2","mm")
    biomass_change_names_sum=c("TgC","TgC","TgC","TgC","TgC","TgC","TgC","TgC")

    if (file.exists(outfile) == FALSE | repair == 1) {

        # counter for sites actually in analysis
        sites_processed=0 ; total_area=0 ; site_locations = 0 ; site_cluster = 0
        for (n in seq(1, PROJECT$nosites)) {
             # generate file name of the output file created in stage 3
             loadfile=paste(PROJECT$results_processedpath,PROJECT$sites[n],"_parameters.RData",sep="")
             if (file.exists(loadfile)) {
                 site_locations=append(site_locations,n)
                 sites_processed=sites_processed+1
                 slot_j=as.numeric(PROJECT$sites[n])/PROJECT$long_dim
                 slot_i=as.numeric(PROJECT$sites[n])-(floor(slot_j)*PROJECT$long_dim)
                 if(slot_i == 0){slot_i = PROJECT$long_dim} ; slot_j=ceiling(slot_j)
                 total_area=total_area+area[slot_i,slot_j]
                 site_cluster = append(site_cluster,uk_cluster_pft[slot_i,slot_j])
             }
        }
        # calculate inverse of total area too
        total_area_1 = 1/total_area
        # remove initial value
        site_locations=site_locations[-1] ; site_cluster = site_cluster[-1]

        # create array we will be filling here if we have not already
        if (exists("states_array_median") == FALSE) {
            # dimension info
            max_pars = 200 ; no_pools = 7
            pixel_dims = c(max_pars,length(timestep_days)) ; pixel_vars = length(in_file_names)
            # spatial / not cluster related values
            states_array_median = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,pixel_vars+nos_extra))
            states_array_unc = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,pixel_vars+nos_extra))
            gsi_control_median = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,3))
            gsi_control_unc = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,3))
            final_stocks_median = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,no_pools+1))
            final_stocks_uncertainty = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,no_pools+1))
            stock_change_rate_median = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,no_pools))
            stock_change_rate_uncertainty = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,no_pools))
            pixel_management_array = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
            cumulative_management_array = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
            # cluster / timeseries related
            states_mean_timeseries = array(0, dim=c(max_pars,pixel_dims[2],pixel_vars,nos_uk_clusters))
            states_sum_timeseries = array(0, dim=c(max_pars,pixel_dims[2],pixel_vars,nos_uk_clusters))
            biomass_change_sum_timeseries = array(0, dim=c(max_pars,pixel_dims[2],length(biomass_change_names),nos_uk_clusters))
            biomass_change_mean_timeseries = array(0, dim=c(max_pars,pixel_dims[2],length(biomass_change_names),nos_uk_clusters))
        }

        for (c in seq(1,nos_uk_clusters)) {print(paste("Cluster ",c," has ",length(which(site_cluster == c))," sites",sep=""))}

        # loop through the clusters generating all flux values
        if (use_parallel & numWorkers > nos_uk_clusters) {
            # bundle needed functions down the chain
            functions_list=c("simulate_all","run_mcmc_results_for_grid","gsi_limiting","gsi_controlling")
            analysis_info=list(site_locations = site_locations,site_cluster = site_cluster,max_pars = max_pars, no_pools = no_pools, timestep_days = timestep_days, tmp_m = tmp_m,
            nos_years = nos_years, steps_per_year = steps_per_year,pixel_vars = pixel_vars,nos_extra = nos_extra,
            in_file_names = in_file_names, par_names = par_names, biomass_change_names = biomass_change_names)
            cl <- makeCluster(nos_uk_clusters, type = "PSOCK",outfile="./parallel_reporting.txt")
            # load R libraries in cluster
            clusterExport(cl,"load_r_libraries") ; clusterEvalQ(cl, load_r_libraries()) ; clusterExport(cl,functions_list)
            out_list=parLapply(cl,c(1:nos_uk_clusters),fun=parallel_cluster_processing,analysis_info=analysis_info,
            uk_cluster_pft = uk_cluster_pft, area = area, area_with_g_Tg = area_with_g_Tg,
            total_area = total_area, total_area_1 = total_area_1,PROJECT = PROJECT)  ; gc() ; gc()
            stopCluster(cl) ; closeAllConnections() ; gc(); gc()
            #; file.remove("./parallel_reporting.txt") # close cluster and delete the reporting file assuming that things were successful...
            print(".........beginning re-arranging of output")

            # extract output from each core and merge into single grid for each cluster
            for (c in seq(1,nos_uk_clusters)) {

                 # do cluster averaged material first
                 for (p in seq(1,pixel_vars)) {
                      states_mean_timeseries[1:max_pars,,p,c]=states_mean_timeseries[1:max_pars,,p,c]+out_list[[c]]$states_mean_timeseries[1:max_pars,,p]
                      states_sum_timeseries[1:max_pars,,p,c]=states_sum_timeseries[1:max_pars,,p,c]+out_list[[c]]$states_sum_timeseries[1:max_pars,,p]
                 }
                 for (p in seq(1,length(biomass_change_names))) {
                      biomass_change_mean_timeseries[1:max_pars,,p,c]=biomass_change_mean_timeseries[1:max_pars,,p,c]+out_list[[c]]$biomass_change_mean_timeseries[1:max_pars,,p]
                      biomass_change_sum_timeseries[1:max_pars,,p,c]=biomass_change_sum_timeseries[1:max_pars,,p,c]+out_list[[c]]$biomass_change_sum_timeseries[1:max_pars,,p]
                 }

                 # determine sites within this cluster
                 sites_to_do = site_locations[which(site_cluster == c)]
                 # loop through every location now to position all the data.
                 for (nn in seq(1,length(sites_to_do))) {
                      # calculate pixel location information
                      n = sites_to_do[nn]
                      slot_j=as.numeric(PROJECT$sites[n])/PROJECT$long_dim
                      slot_i=as.numeric(PROJECT$sites[n])-(floor(slot_j)*PROJECT$long_dim)
                      if(slot_i == 0) {slot_i = PROJECT$long_dim}
                      slot_j=ceiling(slot_j)

                      # model outputs
                      states_array_median[slot_i,slot_j,]=out_list[[c]]$states_array_median[slot_i,slot_j,]
                      states_array_unc[slot_i,slot_j,]=out_list[[c]]$states_array_unc[slot_i,slot_j,]
                      gsi_control_median[slot_i,slot_j,]=out_list[[c]]$gsi_control_median[slot_i,slot_j,]
                      gsi_control_unc[slot_i,slot_j,]=out_list[[c]]$gsi_control_unc[slot_i,slot_j,]
                      final_stocks_median[slot_i,slot_j,]=out_list[[c]]$final_stocks_median[slot_i,slot_j,]
                      final_stocks_uncertainty[slot_i,slot_j,]=out_list[[c]]$final_stocks_uncertainty[slot_i,slot_j,]
                      stock_change_rate_median[slot_i,slot_j,]=out_list[[c]]$stock_change_rate_median[slot_i,slot_j,]
                      stock_change_rate_uncertainty[slot_i,slot_j,]=out_list[[c]]$stock_change_rate_uncertainty[slot_i,slot_j,]
                      # driver / observations outputs
                      pixel_management_array[slot_i,slot_j]=out_list[[c]]$pixel_management_array[slot_i,slot_j]
                      cumulative_management_array[slot_i,slot_j]=out_list[[c]]$cumulative_management_array[slot_i,slot_j]
                 } # loop through all sites not collecting their information
            } # looping clusters
            rm(out_list) ; gc() ; gc()
            print("......done re-arranging")
            gc() ; file.remove("./parallel_reporting.txt")

        } else {

            analysis_info=list(site_locations = site_locations,site_cluster = site_cluster,max_pars = max_pars,
            no_pools = no_pools, timestep_days = timestep_days, tmp_m = tmp_m,
            nos_years = nos_years, steps_per_year = steps_per_year,pixel_vars = pixel_vars,nos_extra = nos_extra,
            in_file_names = in_file_names, par_names = par_names, biomass_change_names = biomass_change_names)
            # or use serial
            out_list=lapply(c(1:nos_uk_clusters),FUN=parallel_cluster_processing,analysis_info=analysis_info,
            uk_cluster_pft = uk_cluster_pft, area = area, area_with_g_Tg = area_with_g_Tg, total_area = total_area,
            total_area_1 = total_area_1,PROJECT = PROJECT) ; gc() ; gc()
            print(".........beginning re-arranging of output")
            # extract output from each core and merge into single grid for each cluster
            for (c in seq(1,nos_uk_clusters)) {

                 # do cluster averaged material first
                 for (p in seq(1,pixel_vars)) {
                      states_mean_timeseries[1:max_pars,,p,c]=states_mean_timeseries[1:max_pars,,p,c]+out_list[[c]]$states_mean_timeseries[1:max_pars,,p]
                      states_sum_timeseries[1:max_pars,,p,c]=states_sum_timeseries[1:max_pars,,p,c]+out_list[[c]]$states_sum_timeseries[1:max_pars,,p]
                 }
                 for (p in seq(1,length(biomass_change_names))) {
                      biomass_change_mean_timeseries[1:max_pars,,p,c]=biomass_change_mean_timeseries[1:max_pars,,p,c]+out_list[[c]]$biomass_change_mean_timeseries[1:max_pars,,p]
                      biomass_change_sum_timeseries[1:max_pars,,p,c]=biomass_change_sum_timeseries[1:max_pars,,p,c]+out_list[[c]]$biomass_change_sum_timeseries[1:max_pars,,p]
                 }

                 # determine sites within this cluster
                 sites_to_do = site_locations[which(site_cluster == c)]
                 # loop through every location now to position all the data.
                 for (nn in seq(1,length(sites_to_do))) {
                      # calculate pixel location information
                      n = sites_to_do[nn]
                      slot_j=as.numeric(PROJECT$sites[n])/PROJECT$long_dim
                      slot_i=as.numeric(PROJECT$sites[n])-(floor(slot_j)*PROJECT$long_dim)
                      if(slot_i == 0) {slot_i = PROJECT$long_dim}
                      slot_j=ceiling(slot_j)

                       # model outputs
                      states_array_median[slot_i,slot_j,]=out_list[[c]]$states_array_median[slot_i,slot_j,]
                      states_array_unc[slot_i,slot_j,]=out_list[[c]]$states_array_unc[slot_i,slot_j,]
                      gsi_control_median[slot_i,slot_j,]=out_list[[c]]$gsi_control_median[slot_i,slot_j,]
                      gsi_control_unc[slot_i,slot_j,]=out_list[[c]]$gsi_control_unc[slot_i,slot_j,]
                      final_stocks_median[slot_i,slot_j,]=out_list[[c]]$final_stocks_median[slot_i,slot_j,]
                      final_stocks_uncertainty[slot_i,slot_j,]=out_list[[c]]$final_stocks_uncertainty[slot_i,slot_j,]
                      stock_change_rate_median[slot_i,slot_j,]=out_list[[c]]$stock_change_rate_median[slot_i,slot_j,]
                      stock_change_rate_uncertainty[slot_i,slot_j,]=out_list[[c]]$stock_change_rate_uncertainty[slot_i,slot_j,]
                      # driver / observations outputs
                      pixel_management_array[slot_i,slot_j]=out_list[[c]]$pixel_management_array[slot_i,slot_j]
                      cumulative_management_array[slot_i,slot_j]=out_list[[c]]$cumulative_management_array[slot_i,slot_j]
                 } # loop through all sites not collecting their information
            } # looping clusters
            rm(out_list) ; gc() ; gc()
            print("......done re-arranging")

        } # parallel or serial option

    } else {

        # if we have loaded these before just load the final product
        load(outfile)

    } # have these previously been loaded

    # inform the user
    print(paste("......in total ",round(((sites_processed/PROJECT$nosites)*100),digits=1),"% of site processed",sep=""))
    print("......now generating flux / state maps")

    # calculate some timing information
    timestep=1
    if (PROJECT$model$timestep != "daily") {timestep=mean(timestep_days)}
    time_vector=1:dim(states_sum_timeseries)[2]
    year_vector=time_vector/(365.25/timestep)
    year_vector=year_vector+as.numeric(PROJECT$start_year)
    interval=floor(length(year_vector)/10)

    # calculate land mask
    landmask=array(PROJECT$landsea, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
    # load colour palette
    colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral")))
    # array sizes are always the same so
    colour_choices = colour_choices_upper(length(par_array_median[,,1]))

    # determine correct height and widths
    fig_height=4000 ; fig_width=7200
    if (PROJECT$grid_type == "UK") { fig_height=8000 ; fig_width=7200 }

    mean_rooting_depth = NA
    if (PROJECT$model$name == "DALEC_GSI_BUCKET" | PROJECT$model$name == "DALECN_GSI_BUCKET") {
        jpeg(file=paste(PROJECT$figpath,"median_root_depth_maps_",PROJECT$name,".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
        par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
        mean_rooting_depth = par_array_median[,,40] * (states_array_median[,,10]*2) / (par_array_median[,,39] + (states_array_median[,,10]*2))
        z_axis=c(min(as.vector(mean_rooting_depth),na.rm=TRUE),max(as.vector(mean_rooting_depth),na.rm=TRUE))
        image.plot(mean_rooting_depth,col=colour_choices, main=paste("Median root depth (m)",sep=""),zlim=z_axis,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
        contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
        dev.off()
    } else if (PROJECT$model$name == "DALEC" | PROJECT$model$name == "DALEC_BUCKET" | PROJECT$model$name == "DALEC_BUCKET_CanAGE") {
        jpeg(file=paste(PROJECT$figpath,"median_root_depth_maps_",PROJECT$name,".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
        par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
        mean_rooting_depth = par_array_median[,,40] * (states_array_median[,,7]*par_array_median[,,29]*2) / (par_array_median[,,39] + (states_array_median[,,7]*par_array_median[,,29]*2))
        z_axis=c(min(as.vector(mean_rooting_depth),na.rm=TRUE),max(as.vector(mean_rooting_depth),na.rm=TRUE))
        image.plot(mean_rooting_depth,col=colour_choices, main=paste("Median root depth (m)",sep=""),zlim=z_axis,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
        contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
        dev.off()
    } else if (PROJECT$model$name == "DALECN_BUCKET") {
        jpeg(file=paste(PROJECT$figpath,"median_root_depth_maps_",PROJECT$name,".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
        par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
        mean_rooting_depth = par_array_median[,,35] * (states_array_median[,,10]*2) / (par_array_median[,,34] + (states_array_median[,,10]*2))
        z_axis=c(min(as.vector(mean_rooting_depth),na.rm=TRUE),max(as.vector(mean_rooting_depth),na.rm=TRUE))
        image.plot(mean_rooting_depth,col=colour_choices, main=paste("Median root depth (m)",sep=""),zlim=z_axis,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
        contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
        dev.off()
    }

    if (grepl("GSI",PROJECT$model$name)) {
        # generate GSI dominent control plots in bespoke fashion
        ## GSI_itemp component with dominent control
        z_axis=c(0,100)
        jpeg(file=paste(PROJECT$figpath,"median_gsi_itemp_control_maps_",PROJECT$name,".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
        par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
        image.plot(gsi_control_median[,,1],col=colour_choices, main=paste("Median GSI min temp dominant control",sep=""),zlim=z_axis,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
        contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
        dev.off()
        jpeg(file=paste(PROJECT$figpath,"median_uncertainty_gsi_itemp_control_maps_",PROJECT$name,".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
        par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
        image.plot(gsi_control_unc[,,1],col=colour_choices, main=paste("Median uncertainty GSI min temp dominant control",sep=""),axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
        contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
        dev.off()
        jpeg(file=paste(PROJECT$figpath,"median_gsi_itemp_control_hist_",PROJECT$name,".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
        par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
        hist(gsi_control_median[,,1], main=paste("Median GSI min temp dominant control",sep=""), cex.main=2.4,cex=1.5)
        dev.off()
        jpeg(file=paste(PROJECT$figpath,"median_uncertainty_gsi_itemp_control_hist_",PROJECT$name,".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
        par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
        hist(gsi_control_unc[,,1], main=paste("Median uncertainty GSI min temp dominant control",sep=""), cex.main=2.4,cex=1.5)
        dev.off()
        ## GSI iphoto
        jpeg(file=paste(PROJECT$figpath,"median_gsi_iphoto_control_maps_",PROJECT$name,".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
        par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
        image.plot(gsi_control_median[,,2],col=colour_choices, main=paste("Median GSI photoperiod dominant control",sep=""),zlim=z_axis,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
        contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
        dev.off()
        jpeg(file=paste(PROJECT$figpath,"median_uncertainty_gsi_iphoto_control_maps_",PROJECT$name,".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
        par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
        image.plot(gsi_control_unc[,,2],col=colour_choices, main=paste("Median uncertainty GSI photoperiod dominant control",sep=""),axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
        contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
        dev.off()
        jpeg(file=paste(PROJECT$figpath,"median_gsi_iphoto_control_hist_",PROJECT$name,".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
        par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
        hist(gsi_control_median[,,2], main=paste("Median GSI photoperiod dominant control",sep=""), cex.main=2.4,cex=1.5)
        dev.off()
        jpeg(file=paste(PROJECT$figpath,"median_uncertainty_gsi_iphoto_control_hist_",PROJECT$name,".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
        par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
        hist(gsi_control_unc[,,2], main=paste("Median uncertainty GSI photoperiod dominant control",sep=""), cex.main=2.4,cex=1.5)
        dev.off()
        ## GSI ivpd
        jpeg(file=paste(PROJECT$figpath,"median_gsi_ivpd_control_maps_",PROJECT$name,".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
        par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
        image.plot(gsi_control_median[,,3],col=colour_choices, main=paste("Median GSI VPD dominant control",sep=""),zlim=z_axis,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
        contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
        dev.off()
        jpeg(file=paste(PROJECT$figpath,"median_uncertainty_gsi_ivpd_control_maps_",PROJECT$name,".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
        par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
        image.plot(gsi_control_unc[,,3],col=colour_choices, main=paste("Median uncertainty GSI VPD dominant control",sep=""),axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
        contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
        dev.off()
        jpeg(file=paste(PROJECT$figpath,"median_gsi_ivpd_control_hist_",PROJECT$name,".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
        par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
        hist(gsi_control_median[,,3], main=paste("Median GSI VPD dominant control",sep=""), cex.main=2.4,cex=1.5)
        dev.off()
        jpeg(file=paste(PROJECT$figpath,"median_uncertainty_gsi_ivpd_control_hist_",PROJECT$name,".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
        par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
        hist(gsi_control_unc[,,3], main=paste("Median uncertainty GSI VPD dominant control",sep=""), cex.main=2.4,cex=1.5)
        dev.off()
    } # GSI

    # now everything has been loaded into a nice array we begin plotting
    # maps
    for (p in seq(1,dim(states_array_median)[3])) {
         if (length(which(as.vector(states_array_median[,,p]) != 0)) > 0) {
              # define z axis for colour consistency
              z_axis=c(min(states_array_median[,,p],na.rm=TRUE),max((states_array_median[,,p])[which( abs(states_array_median[,,p]) < abs(log(0)))],na.rm=TRUE))
              # special case for GSI components
              if ( grepl("GSI",par_names[p]) ) {z_axis=c(0,1)}
              if ( grepl("_limiting",par_names[p]) ) {z_axis=c(0,100)}

              jpeg(file=paste(PROJECT$figpath,"state_flux_maps_mean_",par_names[p],"_",PROJECT$name,".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
              par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
              image.plot(states_array_median[,,p],col=colour_choices, main=paste("Mean ",par_names[p]," (",par_names_units[p],")",sep=""),zlim=z_axis,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
              contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
              dev.off()
              jpeg(file=paste(PROJECT$figpath,"state_flux_uncertainty_maps_",par_names[p],"_",PROJECT$name,".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
              # define z axis for colour consistency
              z_axis=c(min(states_array_unc[,,p],na.rm=TRUE),max((states_array_unc[,,p])[which( abs(states_array_unc[,,p]) < abs(log(0)))],na.rm=TRUE))
              par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
              image.plot(states_array_unc[,,p],col=colour_choices, main=paste("Mean ",par_names[p]," CI range  (",par_names_units[p],")",sep=""),zlim=z_axis,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
              contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
              dev.off()

              if (grepl("GSI",par_names[p])) {
                  # limiting GSI scaler
                  jpeg(file=paste(PROJECT$figpath,"state_flux_hist_mean_",par_names[p],"_",PROJECT$name,".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
                  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
                  hist(states_array_median[,,p], main=paste("Mean ",par_names[p]," (",par_names_units[p],")",sep=""), cex.main=2.4,cex=1.5)
                  dev.off()
                  jpeg(file=paste(PROJECT$figpath,"state_flux_hist_uncertainty_maps_",par_names[p],"_",PROJECT$name,".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
                  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
                  hist(states_array_unc[,,p], main=paste("Mean ",par_names[p]," CI range  (",par_names_units[p],")",sep=""), cex.main=2.4,cex=1.5)
                  dev.off()
              } # if GSI variable
         } # is variable present
    } # loop through state variables

    # time series data
    for (p in seq(1,dim(states_array_median)[3]-nos_extra)) {
         if (length(which(as.vector(states_sum_timeseries[,,p,]) != 0)) > 0) {
              jpeg(file=paste(PROJECT$figpath,"grid_sum_timeseries_",par_names[p],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=300, quality=100)
              info = " " # assume default is no header, but sometimes we add something extra...
              if (par_names[p] == "NEE"){
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_sum_timeseries[sample(1:max_pars,max_pars),,p,c]*timestep
                       } else {
                           tmp_var=tmp_var+states_sum_timeseries[sample(1:max_pars,max_pars),,p,c]*timestep
                       }
                  }
                  var1=sum(apply(t(tmp_var),1,median,na.rm=TRUE))/nos_years
                  var2=sum(apply(t(tmp_var),1,quantile, prob=0.975,na.rm=TRUE))/nos_years
                  var3=sum(apply(t(tmp_var),1,quantile, prob=0.025,na.rm=TRUE))/nos_years
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Cumulative ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (TgC.yr-1)",sep="")
              } # NEE
              if (par_names[p] == "GPP"){
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_sum_timeseries[sample(1:max_pars,max_pars),,p,c]*timestep
                       } else {
                           tmp_var=tmp_var+states_sum_timeseries[sample(1:max_pars,max_pars),,p,c]*timestep
                       }
                  }
                  var1=sum(apply(t(tmp_var),1,median,na.rm=TRUE))/nos_years
                  var2=sum(apply(t(tmp_var),1,quantile, prob=0.975,na.rm=TRUE))/nos_years
                  var3=sum(apply(t(tmp_var),1,quantile, prob=0.025,na.rm=TRUE))/nos_years
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Cumulative ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (TgC.yr-1)",sep="")
              } # GPP
              if (par_names[p] == "Evap"){
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_sum_timeseries[sample(1:max_pars,max_pars),,p,c]*timestep
                       } else {
                           tmp_var=tmp_var+states_sum_timeseries[sample(1:max_pars,max_pars),,p,c]*timestep
                       }
                  }
                  var1=sum(apply(t(tmp_var),1,median,na.rm=TRUE))/nos_years
                  var2=sum(apply(t(tmp_var),1,quantile, prob=0.975,na.rm=TRUE))/nos_years
                  var3=sum(apply(t(tmp_var),1,quantile, prob=0.025,na.rm=TRUE))/nos_years
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Cumulative ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (PgH2O.yr-1)",sep="")
              } # Evap
              if (par_names[p] == "Reco"){
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_sum_timeseries[sample(1:max_pars,max_pars),,p,c]*timestep
                       } else {
                           tmp_var=tmp_var+states_sum_timeseries[sample(1:max_pars,max_pars),,p,c]*timestep
                       }
                  }
                  var1=sum(apply(t(tmp_var),1,median,na.rm=TRUE))/nos_years
                  var2=sum(apply(t(tmp_var),1,quantile, prob=0.975,na.rm=TRUE))/nos_years
                  var3=sum(apply(t(tmp_var),1,quantile, prob=0.025,na.rm=TRUE))/nos_years
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Cumulative ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (TgC.yr-1)",sep="")
              } # Reco
              if (par_names[p] == "Rauto"){
                  for (c in seq(1,nos_uk_clusters)) {
                      if (c == 1) {
                          tmp_var=states_sum_timeseries[sample(1:max_pars,max_pars),,p,c]*timestep
                      } else {
                          tmp_var=tmp_var+states_sum_timeseries[sample(1:max_pars,max_pars),,p,c]*timestep
                      }
                  }
                  var1=sum(apply(t(tmp_var),1,median,na.rm=TRUE))/nos_years
                  var2=sum(apply(t(tmp_var),1,quantile, prob=0.975,na.rm=TRUE))/nos_years
                  var3=sum(apply(t(tmp_var),1,quantile, prob=0.025,na.rm=TRUE))/nos_years
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Cumulative ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (TgC.yr-1)",sep="")
              } # Rauto
              if (par_names[p] == "Rhet"){
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_sum_timeseries[sample(1:max_pars,max_pars),,p,c]*timestep
                       } else {
                           tmp_var=tmp_var+states_sum_timeseries[sample(1:max_pars,max_pars),,p,c]*timestep
                       }
                  }
                  var1=sum(apply(t(tmp_var),1,median,na.rm=TRUE))/nos_years
                  var2=sum(apply(t(tmp_var),1,quantile, prob=0.975,na.rm=TRUE))/nos_years
                  var3=sum(apply(t(tmp_var),1,quantile, prob=0.025,na.rm=TRUE))/nos_years
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Cumulative ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (TgC.yr-1)",sep="")
              } # Rhet
              if (par_names[p] == "Fire"){
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_sum_timeseries[sample(1:max_pars,max_pars),,p,c]*timestep
                       } else {
                           tmp_var=tmp_var+states_sum_timeseries[sample(1:max_pars,max_pars),,p,c]*timestep
                       }
                  }
                  var1=sum(apply(t(tmp_var),1,median,na.rm=TRUE))/nos_years
                  var2=sum(apply(t(tmp_var),1,quantile, prob=0.975,na.rm=TRUE))/nos_years
                  var3=sum(apply(t(tmp_var),1,quantile, prob=0.025,na.rm=TRUE))/nos_years
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Cumulative ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (TgC.yr-1)",sep="")
              } # Fire
              if (par_names[p] == "Charvested"){
                  for (c in seq(1,nos_uk_clusters)) {
                      if (c == 1) {
                          tmp_var=states_sum_timeseries[sample(1:max_pars,max_pars),,p,c]*timestep
                      } else {
                          tmp_var=tmp_var+states_sum_timeseries[sample(1:max_pars,max_pars),,p,c]*timestep
                      }
                  }
                  var1=sum(apply(t(tmp_var),1,median,na.rm=TRUE))/nos_years
                  var2=sum(apply(t(tmp_var),1,quantile, prob=0.975,na.rm=TRUE))/nos_years
                  var3=sum(apply(t(tmp_var),1,quantile, prob=0.025,na.rm=TRUE))/nos_years
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Cumulative ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (TgC.yr-1)",sep="")
              } # Charvest
              if (par_names[p] == "Cfoliage"){
                  # at the end of the analysis
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_sum_timeseries[sample(1:max_pars,max_pars),dim(states_sum_timeseries)[2],p,c]
                       } else {
                           tmp_var=tmp_var+states_sum_timeseries[sample(1:max_pars,max_pars),dim(states_sum_timeseries)[2],p,c]
                       }
                  }
                  var1=median(tmp_var) ; var2=quantile(tmp_var, prob=0.975,na.rm=TRUE) ; var3=quantile(tmp_var, prob=0.025,na.rm=TRUE)
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Final ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (TgC)",sep="")
              } # Cfoliage
              if (par_names[p] == "Croot"){
                  # at the end of the analysis
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_sum_timeseries[sample(1:max_pars,max_pars),dim(states_sum_timeseries)[2],p,c]
                       } else {
                           tmp_var=tmp_var+states_sum_timeseries[sample(1:max_pars,max_pars),dim(states_sum_timeseries)[2],p,c]
                       }

                  }
                  var1=median(tmp_var) ; var2=quantile(tmp_var, prob=0.975,na.rm=TRUE) ; var3=quantile(tmp_var, prob=0.025,na.rm=TRUE)
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Final ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (TgC)",sep="")
              } # Croot
              if (par_names[p] == "Clabile"){
                  # at the end of the analysis
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_sum_timeseries[sample(1:max_pars,max_pars),dim(states_sum_timeseries)[2],p,c]
                       } else {
                           tmp_var=tmp_var+states_sum_timeseries[sample(1:max_pars,max_pars),dim(states_sum_timeseries)[2],p,c]
                       }
                  }
                  var1=median(tmp_var) ; var2=quantile(tmp_var, prob=0.975,na.rm=TRUE) ; var3=quantile(tmp_var, prob=0.025,na.rm=TRUE)
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Final ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (TgC)",sep="")
              } # Clabile
              if (par_names[p] == "Cwood"){
                  # at the end of the analysis
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_sum_timeseries[sample(1:max_pars,max_pars),dim(states_sum_timeseries)[2],p,c]
                       } else {
                           tmp_var=tmp_var+states_sum_timeseries[sample(1:max_pars,max_pars),dim(states_sum_timeseries)[2],p,c]
                       }
                  }
                  var1=median(tmp_var) ; var2=quantile(tmp_var, prob=0.975,na.rm=TRUE) ; var3=quantile(tmp_var, prob=0.025,na.rm=TRUE)
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Final ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (TgC)",sep="")
              } # Cwood
              if (par_names[p] == "Ccwd"){
                  # at the end of the analysis
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_sum_timeseries[sample(1:max_pars,max_pars),dim(states_sum_timeseries)[2],p,c]
                       } else {
                           tmp_var=tmp_var+states_sum_timeseries[sample(1:max_pars,max_pars),dim(states_sum_timeseries)[2],p,c]
                       }
                  }
                  var1=median(tmp_var) ; var2=quantile(tmp_var, prob=0.975,na.rm=TRUE) ; var3=quantile(tmp_var, prob=0.025,na.rm=TRUE)
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Final ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (TgC)",sep="")
              } # Ccwd
              if (par_names[p] == "Clitter"){
                  # at the end of the analysis
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_sum_timeseries[sample(1:max_pars,max_pars),dim(states_sum_timeseries)[2],p,c]
                       } else {
                           tmp_var=tmp_var+states_sum_timeseries[sample(1:max_pars,max_pars),dim(states_sum_timeseries)[2],p,c]
                       }
                  }
                  var1=median(tmp_var) ; var2=quantile(tmp_var, prob=0.975,na.rm=TRUE) ; var3=quantile(tmp_var, prob=0.025,na.rm=TRUE)
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Final ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (TgC)",sep="")
              } # Clitter
              if (par_names[p] == "Csom"){
                  # at the end of the analysis
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_sum_timeseries[sample(1:max_pars,max_pars),dim(states_sum_timeseries)[2],p,c]
                       } else {
                           tmp_var=tmp_var+states_sum_timeseries[sample(1:max_pars,max_pars),dim(states_sum_timeseries)[2],p,c]
                       }
                  }
                  var1=median(tmp_var) ; var2=quantile(tmp_var, prob=0.975,na.rm=TRUE) ; var3=quantile(tmp_var, prob=0.025,na.rm=TRUE)
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Final ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (TgC)",sep="")
              } # Csom
              if (par_names[p] == "CTotal"){
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_sum_timeseries[sample(1:max_pars,max_pars),dim(states_sum_timeseries)[2],p,c]
                       } else {
                           tmp_var=tmp_var+states_sum_timeseries[sample(1:max_pars,max_pars),dim(states_sum_timeseries)[2],p,c]
                       }
                  }
                  var1=median(tmp_var) ; var2=quantile(tmp_var, prob=0.975,na.rm=TRUE) ; var3=quantile(tmp_var, prob=0.025,na.rm=TRUE)
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Final ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (TgC)",sep="")
              } # CTotal
              par(mfrow=c(1,1), mar=c(5,5,3,1))
              for (c in seq(1,nos_uk_clusters)) {if (c == 1) {ymin=states_sum_timeseries[sample(1:max_pars,max_pars),,p,c]} else {ymin=ymin+states_sum_timeseries[sample(1:max_pars,max_pars),,p,c]}}
              ymax=apply(ymin,2,quantile,prob=c(0.975)); ymax=quantile(ymax[which(abs(ymax) != abs(log(0)))],na.rm=TRUE,prob=c(0.975)) ; ymin=min(apply(ymin,2,quantile,prob=c(0.001)))
              plot(rep(-9999,dim(states_sum_timeseries)[2]), pch=16,xaxt="n",ylim=c(ymin,ymax),
              cex=0.8,ylab=paste(par_names[p]," (",par_names_units_sum[p],")",sep=""),xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=info)
              axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
              # add the confidence intervals
              plotconfidence(states_sum_timeseries[,,p,])
              # calculate and draw the median values, could be mean instead or other
              for (c in seq(1,nos_uk_clusters)) {
                   if (c == 1) {
                       tmp_var=states_sum_timeseries[sample(1:max_pars,max_pars),,p,c]
                   } else {
                       tmp_var=tmp_var+states_sum_timeseries[sample(1:max_pars,max_pars),,p,c]
                   }
              }
              lines(apply(tmp_var,2,median), lwd=1,col="blue")
              dev.off()
              jpeg(file=paste(PROJECT$figpath,"grid_averaged_timeseries_",par_names[p],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=300, quality=100)
              info = " " # assume default is no header, but sometimes we add something extra...
              if (par_names[p] == "NEE"){
                  for (c in seq(1,nos_uk_clusters)) {
                      if (c == 1) {
                          tmp_var=states_mean_timeseries[sample(1:max_pars,max_pars),,p,c]
                      } else {
                          tmp_var=tmp_var+states_mean_timeseries[sample(1:max_pars,max_pars),,p,c]
                      }
                  }
                  var1=sum(apply(t(tmp_var*mean(timestep_days)),1,median,na.rm=TRUE))/nos_years
                  var2=sum(apply(t(tmp_var*mean(timestep_days)),1,quantile, prob=0.975,na.rm=TRUE))/nos_years
                  var3=sum(apply(t(tmp_var*mean(timestep_days)),1,quantile, prob=0.025,na.rm=TRUE))/nos_years
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Mean annual ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (gC.m-2.yr-1)",sep="")
              } # NEE
              if (par_names[p] == "CTotal"){
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_mean_timeseries[sample(1:max_pars,max_pars),,p,c]
                       } else {
                           tmp_var=tmp_var+states_mean_timeseries[sample(1:max_pars,max_pars),,p,c]
                       }
                  }
                  var1=sum(apply(t(tmp_var*mean(timestep_days)),1,median,na.rm=TRUE))/nos_years
                  var2=sum(apply(t(tmp_var*mean(timestep_days)),1,quantile, prob=0.975,na.rm=TRUE))/nos_years
                  var3=sum(apply(t(tmp_var*mean(timestep_days)),1,quantile, prob=0.025,na.rm=TRUE))/nos_years
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Mean annual ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (gC.m-2.yr-1)",sep="")
              } # CTotal
              if (par_names[p] == "GPP"){
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_mean_timeseries[sample(1:max_pars,max_pars),,p,c]
                       } else {
                           tmp_var=tmp_var+states_mean_timeseries[sample(1:max_pars,max_pars),,p,c]
                       }
                  }
                  var1=sum(apply(t(tmp_var*mean(timestep_days)),1,median,na.rm=TRUE))/nos_years
                  var2=sum(apply(t(tmp_var*mean(timestep_days)),1,quantile, prob=0.975,na.rm=TRUE))/nos_years
                  var3=sum(apply(t(tmp_var*mean(timestep_days)),1,quantile, prob=0.025,na.rm=TRUE))/nos_years
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Mean annual ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (gC.m-2.yr-1)",sep="")
              } # GPP
              if (par_names[p] == "Rhet"){
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_mean_timeseries[sample(1:max_pars,max_pars),,p,c]
                       } else {
                           tmp_var=tmp_var+states_mean_timeseries[sample(1:max_pars,max_pars),,p,c]
                       }
                  }
                  var1=sum(apply(t(tmp_var*mean(timestep_days)),1,median,na.rm=TRUE))/nos_years
                  var2=sum(apply(t(tmp_var*mean(timestep_days)),1,quantile, prob=0.975,na.rm=TRUE))/nos_years
                  var3=sum(apply(t(tmp_var*mean(timestep_days)),1,quantile, prob=0.025,na.rm=TRUE))/nos_years
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Mean annual ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (gC.m-2.yr-1)",sep="")
              } # Rhet
              if (par_names[p] == "Reco"){
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_mean_timeseries[sample(1:max_pars,max_pars),,p,c]
                       } else {
                           tmp_var=tmp_var+states_mean_timeseries[sample(1:max_pars,max_pars),,p,c]
                       }
                  }
                  var1=sum(apply(t(tmp_var*mean(timestep_days)),1,median,na.rm=TRUE))/nos_years
                  var2=sum(apply(t(tmp_var*mean(timestep_days)),1,quantile, prob=0.975,na.rm=TRUE))/nos_years
                  var3=sum(apply(t(tmp_var*mean(timestep_days)),1,quantile, prob=0.025,na.rm=TRUE))/nos_years
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Mean annual ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (gC.m-2.yr-1)",sep="")
              } # Reco
              if (par_names[p] == "Rauto"){
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_mean_timeseries[sample(1:max_pars,max_pars),,p,c]
                       } else {
                           tmp_var=tmp_var+states_mean_timeseries[sample(1:max_pars,max_pars),,p,c]
                       }
                  }
                  var1=sum(apply(t(tmp_var*mean(timestep_days)),1,median,na.rm=TRUE))/nos_years
                  var2=sum(apply(t(tmp_var*mean(timestep_days)),1,quantile, prob=0.975,na.rm=TRUE))/nos_years
                  var3=sum(apply(t(tmp_var*mean(timestep_days)),1,quantile, prob=0.025,na.rm=TRUE))/nos_years
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Mean annual ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (gC.m-2.yr-1)",sep="")
              } # Rauto
              if (par_names[p] == "Evap"){
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_mean_timeseries[sample(1:max_pars,max_pars),,p,c]
                       } else {
                           tmp_var=tmp_var+states_mean_timeseries[sample(1:max_pars,max_pars),,p,c]
                       }
                  }
                  var1=sum(apply(t(tmp_var*mean(timestep_days)),1,median,na.rm=TRUE))/nos_years
                  var2=sum(apply(t(tmp_var*mean(timestep_days)),1,quantile, prob=0.975,na.rm=TRUE))/nos_years
                  var3=sum(apply(t(tmp_var*mean(timestep_days)),1,quantile, prob=0.025,na.rm=TRUE))/nos_years
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Mean annual ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (kgH2O.m-2.yr-1)",sep="")
              }  # Evap
              if (par_names[p] == "Fire"){
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_mean_timeseries[sample(1:max_pars,max_pars),,p,c]
                       } else {
                           tmp_var=tmp_var+states_mean_timeseries[sample(1:max_pars,max_pars),,p,c]
                       }
                  }
                  var1=sum(apply(t(tmp_var*mean(timestep_days)),1,median,na.rm=TRUE))/nos_years
                  var2=sum(apply(t(tmp_var*mean(timestep_days)),1,quantile, prob=0.975,na.rm=TRUE))/nos_years
                  var3=sum(apply(t(tmp_var*mean(timestep_days)),1,quantile, prob=0.025,na.rm=TRUE))/nos_years
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Mean annual ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (gC.m-2.yr-1)",sep="")
              } # Fire
              if (par_names[p] == "Charvested"){
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_mean_timeseries[sample(1:max_pars,max_pars),,p,c]
                       } else {
                           tmp_var=tmp_var+states_mean_timeseries[sample(1:max_pars,max_pars),,p,c]
                       }
                  }
                  var1=sum(apply(t(tmp_var*mean(timestep_days)),1,median,na.rm=TRUE))/nos_years
                  var2=sum(apply(t(tmp_var*mean(timestep_days)),1,quantile, prob=0.975,na.rm=TRUE))/nos_years
                  var3=sum(apply(t(tmp_var*mean(timestep_days)),1,quantile, prob=0.025,na.rm=TRUE))/nos_years
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Mean annual ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (gC.m-2.yr-1)",sep="")
              } # Charvested
              if (par_names[p] == "Cfoliage"){
                  # at the end of the analysis
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_mean_timeseries[sample(1:max_pars,max_pars),dim(states_mean_timeseries)[2],p,c]
                       } else {
                           tmp_var=tmp_var+states_mean_timeseries[sample(1:max_pars,max_pars),dim(states_mean_timeseries)[2],p,c]
                       }
                  }
                  var1=median(tmp_var) ; var2=quantile(tmp_var, prob=0.975,na.rm=TRUE) ; var3=quantile(tmp_var, prob=0.025,na.rm=TRUE)
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Mean Final ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (gC.m-2)",sep="")
              } # Cfoliage
              if (par_names[p] == "Croot"){
                  # at the end of the analysis
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_mean_timeseries[sample(1:max_pars,max_pars),dim(states_mean_timeseries)[2],p,c]
                       } else {
                           tmp_var=tmp_var+states_mean_timeseries[sample(1:max_pars,max_pars),dim(states_mean_timeseries)[2],p,c]
                       }
                  }
                  var1=median(tmp_var) ; var2=quantile(tmp_var, prob=0.975,na.rm=TRUE) ; var3=quantile(tmp_var, prob=0.025,na.rm=TRUE)
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Mean Final ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (gC.m-2)",sep="")
              } # Croot
              if (par_names[p] == "Clabile"){
                # at the end of the analysis
                for (c in seq(1,nos_uk_clusters)) {
                  if (c == 1) {
                    tmp_var=states_mean_timeseries[sample(1:max_pars,max_pars),dim(states_mean_timeseries)[2],p,c]
                  } else {
                    tmp_var=tmp_var+states_mean_timeseries[sample(1:max_pars,max_pars),dim(states_mean_timeseries)[2],p,c]
                  }

                }
                var1=median(tmp_var) ; var2=quantile(tmp_var, prob=0.975,na.rm=TRUE) ; var3=quantile(tmp_var, prob=0.025,na.rm=TRUE)
                var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                info = paste("Mean Final ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (gC.m-2)",sep="")
              } # Clabile
              if (par_names[p] == "Clitter"){
                  # at the end of the analysis
                  for (c in seq(1,nos_uk_clusters)) {
                      if (c == 1) {
                          tmp_var=states_mean_timeseries[sample(1:max_pars,max_pars),dim(states_mean_timeseries)[2],p,c]
                      } else {
                          tmp_var=tmp_var+states_mean_timeseries[sample(1:max_pars,max_pars),dim(states_mean_timeseries)[2],p,c]
                      }
                  }
                  var1=median(tmp_var) ; var2=quantile(tmp_var, prob=0.975,na.rm=TRUE) ; var3=quantile(tmp_var, prob=0.025,na.rm=TRUE)
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Mean Final ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (gC.m-2)",sep="")
              } # Clitter
              if (par_names[p] == "Ccwd"){
                  # at the end of the analysis
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_mean_timeseries[sample(1:max_pars,max_pars),dim(states_mean_timeseries)[2],p,c]
                       } else {
                           tmp_var=tmp_var+states_mean_timeseries[sample(1:max_pars,max_pars),dim(states_mean_timeseries)[2],p,c]
                       }
                  }
                  var1=median(tmp_var) ; var2=quantile(tmp_var, prob=0.975,na.rm=TRUE) ; var3=quantile(tmp_var, prob=0.025,na.rm=TRUE)
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Mean Final ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (gC.m-2)",sep="")
              } # Ccwd
              if (par_names[p] == "Csom"){
                  # at the end of the analysis
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_mean_timeseries[sample(1:max_pars,max_pars),dim(states_mean_timeseries)[2],p,c]
                       } else {
                           tmp_var=tmp_var+states_mean_timeseries[sample(1:max_pars,max_pars),dim(states_mean_timeseries)[2],p,c]
                       }
                  }
                  var1=median(tmp_var) ; var2=quantile(tmp_var, prob=0.975,na.rm=TRUE) ; var3=quantile(tmp_var, prob=0.025,na.rm=TRUE)
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Mean Final ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (gC.m-2)",sep="")
              } # Csom
              if (par_names[p] == "Cwood"){
                  # at the end of the analysis
                  for (c in seq(1,nos_uk_clusters)) {
                       if (c == 1) {
                           tmp_var=states_mean_timeseries[sample(1:max_pars,max_pars),dim(states_mean_timeseries)[2],p,c]
                       } else {
                           tmp_var=tmp_var+states_mean_timeseries[sample(1:max_pars,max_pars),dim(states_mean_timeseries)[2],p,c]
                       }
                  }
                  var1=median(tmp_var) ; var2=quantile(tmp_var, prob=0.975,na.rm=TRUE) ; var3=quantile(tmp_var, prob=0.025,na.rm=TRUE)
                  var1=round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
                  info = paste("Mean Final ",par_names[p]," median estimate = ",var1,"; 97.5 % = ",var2,"; 2.5 % = ",var3," (gC.m-2)",sep="")
              } # Cwood
              for (c in seq(1,nos_uk_clusters)) {if (c == 1) {ymin=states_mean_timeseries[sample(1:max_pars,max_pars),,p,c]} else {ymin=ymin+states_mean_timeseries[sample(1:max_pars,max_pars),,p,c]}}
              ymax=apply(ymin,2,quantile,prob=c(0.975)); ymax=quantile(ymax[which(abs(ymax) != abs(log(0)))],na.rm=TRUE,prob=c(0.975)) ; ymin=min(apply(ymin,2,quantile,prob=c(0.001)))
              par(mfrow=c(1,1), mar=c(5,5,3,1))
              plot(rep(-9999,dim(states_mean_timeseries)[2]), pch=16,xaxt="n",ylim=c(ymin,ymax),
              cex=0.8,ylab=paste(par_names[p]," (",par_names_units[p],")",sep=""),xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=info)
              axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
              # add the confidence intervals
              plotconfidence(states_mean_timeseries[,,p,])
              # calculate and draw the median values, could be mean instead or other
              for (c in seq(1,nos_uk_clusters)) {
                   if (c == 1) {
                       tmp_var=states_mean_timeseries[sample(1:max_pars,max_pars),,p,c]
                   } else {
                       tmp_var=tmp_var+states_mean_timeseries[sample(1:max_pars,max_pars),,p,c]
                   }
              }
              lines(apply(tmp_var,2,median), lwd=1,col="blue")
              dev.off()
         }
    }

    # output some aggragated values
    # probably best to add some aggregated met drivers to this concoction here
    save(landmask,area,states_array_median,states_array_unc,states_sum_timeseries,states_mean_timeseries,
         biomass_change_mean_timeseries,biomass_change_sum_timeseries,gsi_control_median,gsi_control_unc,
         biomass_change_names,biomass_change_names_mean,biomass_change_names_sum,final_stocks_median,
         final_stocks_uncertainty,pixel_management_array,cumulative_management_array,
         stock_change_rate_median,stock_change_rate_uncertainty,
         mean_rooting_depth,sites_processed,par_names,max_pars,file=outfile)#, compress="gzip")

         # tidy before leaving
         gc(reset=TRUE, verbose=FALSE)

} # end function generate_stocks_and_fluxes_maps

## Use byte compile
#generate_stocks_and_fluxes_maps<-cmpfun(generate_stocks_and_fluxes_maps)
