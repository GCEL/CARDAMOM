
###
## Function to create generic plots of gridded CARDAMOM analysis parameters and emergent traits
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

generate_parameter_maps<-function(PROJECT) {

   # Determine the lat / long for the grid
   if (PROJECT$grid_type == "UK") {
       output = generate_uk_grid(PROJECT$latitude,PROJECT$longitude,PROJECT$resolution)
       grid_lat = array(output$lat, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
       grid_long = array(output$long,dim=c(PROJECT$long_dim,PROJECT$lat_dim))
       rm(output)
   } else if (PROJECT$grid_type == "wgs84") {
       # generate the lat / long grid again
       output = generate_wgs84_grid(PROJECT$latitude,PROJECT$longitude,PROJECT$resolution)
       grid_lat = array(output$lat, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
       grid_long = array(output$long,dim=c(PROJECT$long_dim,PROJECT$lat_dim))
       rm(output)
   } else {
       stop("valid spatial grid option not selected (UK, or wgs84)")
   }

   # Move working directory
   old_wd = getwd() ; setwd(PROJECT$figpath)

   # Which quantiles will we extract
   na_flag = TRUE
   median_loc = 4 ; upper_loc = 7 ; lower_loc = 1 # 0.50, 0.025, 0.975 assumed

   # Loaded the grid aggregated dataset into memory
   infile = paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep="")
   load(infile)
   # Ensure that any Inf values are removed to NA
   # This is a hack implemented during a development phase and should be removed
   # or commented out under general operation
   for (i in seq(1, length(grid_output))) {
        grid_output[[i]][which(is.infinite(grid_output[[i]]))] = NA
   }

   # load timesteps to local variable
   timestep_days = PROJECT$model$timestep_days ; seconds_per_day = 86400
   timestep_days = rep(timestep_days, length.out=grid_output$time_dim)

   # Convert avgN log10-normal to gN/m2, in models which use foliar N in gN/m2
   if (PROJECT$model$name == "DALEC_GSI_BUCKET"| PROJECT$model$name == "DALEC_GSI_DFOL_CWD_FR" |
       PROJECT$model$name == "DALECN_GSI_BUCKET" | PROJECT$model$name == "DALEC_BUCKET" |
       PROJECT$model$name == "DALEC" | PROJECT$model$name == "DALEC_G5" |
       PROJECT$model$name == "DALEC_G6" | PROJECT$model$name == "DALEC_BUCKET_CanAGE") {
       grid_output$parameters[,,11,] = 10**grid_output$parameters[,,11,]
   }

  # Determine correct height and widths
  hist_height = 4000*0.65 ; hist_width = 7200*0.65
  fig_height = 3000*0.65 ; fig_width = ((PROJECT$long_dim/PROJECT$lat_dim)+0.25) * fig_height
  if (PROJECT$grid_type == "UK") { fig_height = 8000*0.65 ; fig_width = 7200*0.65 }
  # load colour palette
  colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral")))

  # inform the user
  print("......have finished loading - now beginning cluster analysis")

  if (exists(x = "clusters", where = grid_output) == FALSE | repair == 1) {
#  if (exists(x = "clusters", where = grid_output) == FALSE) {
      if (grepl("DALEC",PROJECT$model$name)) {

          # remove non-constrained parameters (i.e. those not actually used in this analysis)
          initial_conditions=c(18:23)
          par_array_median_normalised = grid_output$parameters[,,-c(initial_conditions),median_loc]
          if (PROJECT$model$name == "DALEC_GSI_DFOL_CWD_FR") {
              par_array_median_normalised = grid_output$parameters[,,-c(initial_conditions,30,31,32,33,37),median_loc]
          } else if (PROJECT$model$name == "DALEC_BUCKET") {
              par_array_median_normalised = grid_output$parameters[,,-c(initial_conditions,30,31,32,33,37,41),median_loc]
          } else if (PROJECT$model$name == "DALEC") {
              par_array_median_normalised = grid_output$parameters[,,-c(initial_conditions,30,31,32,33,37),median_loc]
          } else if (PROJECT$model$name == "DALEC_GSI_BUCKET") {
              par_array_median_normalised = grid_output$parameters[,,-c(initial_conditions,30,31,32,33,35,37),median_loc]
          } else if (PROJECT$model$name == "DALEC_GSI_BUCKET_CanAGE") {
              par_array_median_normalised = grid_output$parameters[,,-c(initial_conditions,30,31,32,33,35,37,41),median_loc]
          } else if (PROJECT$model$name == "DALEC_CDEA_ACM2_BUCKET") {
              par_array_median_normalised = grid_output$parameters[,,-c(initial_conditions,24),median_loc]
          } else if (PROJECT$model$name == "DALEC_CDEA_ACM_FARQUHAR_BUCKET") {
              par_array_median_normalised = grid_output$parameters[,,-c(initial_conditions,24),median_loc]
          } else if (PROJECT$model$name == "DALEC_CDEA_ACM2_BUCKET_wMRT") {
              par_array_median_normalised = grid_output$parameters[,,-c(initial_conditions,24),median_loc]
          } else if (PROJECT$model$name == "DALEC_CDEA_ACM2_BUCKET_LAB") {
              par_array_median_normalised = grid_output$parameters[,,-c(initial_conditions,24),median_loc]
          } else if (PROJECT$model$name == "DALEC_CDEA_ACM2_BUCKET_LAB_wMRT") {
              par_array_median_normalised = grid_output$parameters[,,-c(initial_conditions,24),median_loc]
          } else if (PROJECT$model$name == "DALEC_CDEA_ACM2_BUCKET_RmRg") {
              par_array_median_normalised = grid_output$parameters[,,-c(initial_conditions,24),median_loc]
          } else if (PROJECT$model$name == "DALEC_CDEA_ACM2_BUCKET_RmRg_CWD") {
              par_array_median_normalised = grid_output$parameters[,,-c(initial_conditions,24,28),median_loc]
          } else if (PROJECT$model$name == "DALEC_CDEA_ACM2_BUCKET_RmRg_CWD_wMRT") {
              par_array_median_normalised = grid_output$parameters[,,-c(initial_conditions,24,28),median_loc]
          } else if (PROJECT$model$name == "DALEC_G5" | PROJECT$model$name == "DALEC_G6") {
              par_array_median_normalised = grid_output$parameters[,,-c(initial_conditions,37,41),median_loc]
          } else if (PROJECT$model$name == "DALEC_1005" | PROJECT$model$name == "DALEC_1005a") {
              par_array_median_normalised = grid_output$parameters[,,-c(initial_conditions,27,36),median_loc]
          } # PROJECT$model$name == "DALEC_GSI_DFOL_CWD_FR"

          ###
          ## Begin Affinity propagation clustering - for process parameters only

          # Normalise the parameter values as this helps with the cluster
          # analyses to account for very different parameter ranges
          for (i in seq(1,dim(par_array_median_normalised)[3])) {
               min_par_val=min(par_array_median_normalised[,,i],na.rm=TRUE)
               max_par_val=max(par_array_median_normalised[,,i],na.rm=TRUE)
               par_array_median_normalised[,,i] = ((par_array_median_normalised[,,i]-min_par_val)/(max_par_val-min_par_val))
          }

          # Create temporary arrays needed to allow removing of NAs and convert array into (space,par)
          par_array_tmp=array(NA,dim=c(prod(dim(grid_output$parameters)[1:2]),dim(par_array_median_normalised)[3]))
          par_array_tmp[1:prod(dim(par_array_median_normalised)[1:2]),1:dim(par_array_median_normalised)[3]]=par_array_median_normalised
          actual_forests=which(is.na(par_array_tmp[,1]) == FALSE)
          par_array_tmp=par_array_tmp[actual_forests,]
          par_array_tmp=array(par_array_tmp,dim=c((length(par_array_tmp)/dim(par_array_median_normalised)[3]),dim(par_array_median_normalised)[3]))

          # Looping to find preference_input, the preferenceRange() returns 2 values,
          # the first of which minimises the number of clusters, while the seconds
          # would return as many clusters as there are observations. It is the
          # responsibility of the user to ensure the most appropriate use of these
          # information to result in an appropriate number of clusters for error propagation
          tmp = 0
          for (i in seq(1,3)) {
               tmp=append(tmp,preferenceRange(negDistMat(par_array_tmp[sample(1:dim(par_array_tmp)[1],0.05*dim(par_array_tmp)[1], replace=FALSE),],r=2))[1])
          } ; preference_input=max(mean(tmp[-1]),median(tmp[-1]))

          grid_output$cluster_analysis=apclusterL(negDistMat(r=2),par_array_tmp,frac=0.1,sweeps=4, p=preference_input,maxits=1000, convits=100)
          #grid_output$cluster_analysis=apclusterL(negDistMat(r=2),par_array_tmp,frac=0.1, sweeps=10, q=0.05, maxits=1000, convits=100)
          grid_output$nos_pars_clusters=length(grid_output$cluster_analysis@clusters) ; clusters_exemplars=grid_output$cluster_analysis@exemplars
          grid_output$pars_clusters=array(NA,dim=c(dim(par_array_median_normalised)[1:2]))
          for (i in seq(1,length(grid_output$cluster_analysis@clusters))) {
               grid_output$pars_clusters[actual_forests[grid_output$cluster_analysis@clusters[[i]]]] = i
          }
          grid_output$pars_clusters=array(grid_output$pars_clusters,dim=c(dim(par_array_median_normalised)[1:2]))
          # Tidy away the overall analysis in faviour of what we have extracted
          grid_output = within(grid_output, rm(cluster_analysis))

          ###
          ## Begin Affinity propagation clustering - for process parameters + initial conditions

          # Now generate cluster map using all parameters including the initial conditions
          par_array_median_normalised = grid_output$parameters[,,,median_loc]
          # now normalise the parameter values
          for (i in seq(1,dim(par_array_median_normalised)[3])) {
               min_par_val=min(par_array_median_normalised[,,i],na.rm=TRUE)
               max_par_val=max(par_array_median_normalised[,,i],na.rm=TRUE)
               par_array_median_normalised[,,i]=((par_array_median_normalised[,,i]-min_par_val)/(max_par_val-min_par_val))
          }
          # Create temporary arrays needed to allow removing of NAs and convert array into (space,par)
          par_array_tmp=array(NA,dim=c(prod(dim(grid_output$parameters)[1:2]),dim(par_array_median_normalised)[3]))
          par_array_tmp[1:prod(dim(par_array_median_normalised)[1:2]),1:dim(par_array_median_normalised)[3]]=par_array_median_normalised
          actual_forests=which(is.na(par_array_tmp[,1]) == FALSE)
          par_array_tmp=par_array_tmp[actual_forests,]
          par_array_tmp=array(par_array_tmp,dim=c((length(par_array_tmp)/dim(par_array_median_normalised)[3]),dim(par_array_median_normalised)[3]))

          # Looping to find preference_input, the preferenceRange() returns 2 values,
          # the first of which minimises the number of clusters, while the seconds
          # would return as many clusters as there are observations. It is the
          # responsibility of the user to ensure the most appropriate use of these
          # information to result in an appropriate number of clusters for error propagation
          tmp = 0
          for (i in seq(1,3)) {
               tmp=append(tmp,preferenceRange(negDistMat(par_array_tmp[sample(1:dim(par_array_tmp)[1],0.05*dim(par_array_tmp)[1], replace=FALSE),],r=2))[1])
          } ; preference_input=max(mean(tmp[-1]),median(tmp[-1]))

          grid_output$cluster_analysis=apclusterL(negDistMat(r=2),par_array_tmp,frac=0.1,sweeps=4, p=preference_input,maxits=1000, convits=100)
          #grid_output$cluster_analysis=apclusterL(negDistMat(r=2),par_array_tmp,frac=0.1, sweeps=10, q=0.05, maxits=1000, convits=100)
          grid_output$nos_clusters=length(grid_output$cluster_analysis@clusters) ; clusters_exemplars=grid_output$cluster_analysis@exemplars
          grid_output$clusters=array(NA,dim=c(dim(par_array_median_normalised)[1:2]))
          for (i in seq(1,length(grid_output$cluster_analysis@clusters))) {
               grid_output$clusters[actual_forests[grid_output$cluster_analysis@clusters[[i]]]] = i
          }
          grid_output$clusters=array(grid_output$clusters,dim=c(dim(par_array_median_normalised)[1:2]))
          # Tidy away the overall analysis in faviour of what we have extracted
          grid_output = within(grid_output, rm(cluster_analysis))

          # Save the now updated cluster analysis into the grid_output
          save(grid_output, file = infile)

          # Now plot both possible cluster maps
          figname = paste("Cluster_map_of_median_parameters_",gsub("%","_",PROJECT$name),".jpeg",sep="")
          jpeg(file=figname, width=fig_width, height=fig_height, res=300, quality=100)
          par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
          image.plot(x = grid_long, y = grid_lat, z = grid_output$pars_clusters, main=paste("Parameter based cluster maps",sep=""),axes=FALSE, cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
          map(add=TRUE, lwd = 2)
          dev.off()

          figname = paste("Cluster_map_of_median_parameters_with_initial_",gsub("%","_",PROJECT$name),".jpeg",sep="")
          jpeg(file=figname, width=fig_width, height=fig_height, res=300, quality=100)
          par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
          image.plot(x = grid_long, y = grid_lat, z = grid_output$clusters, main=paste("Parameter + initial based cluster map",sep=""),axes=FALSE, cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
          map(add=TRUE, lwd = 2)
          dev.off()

      } # make cluster map

  } # reprocessing or not?

  # inform the user
  print("......now generating parameter maps")

  # describe parameter numbers with log liklihood added on the end
  character_bit=rep("p",times=(PROJECT$model$nopars[1]))
  number_bit=1:(PROJECT$model$nopars[1])
  # merge the letter and numbers together
  par_names=c(paste(character_bit,number_bit,sep=""),"log-likelihood")

  # create tempory array of parameters for the covariance analysis
  figname = paste("parameter_correlations_median_",gsub("%","_",PROJECT$name),".jpeg",sep="")
  jpeg(file=figname, width=fig_width, height=fig_height, res=300, quality=100)
  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
  tmp=array(grid_output$parameters[,,1:dim(grid_output$parameters)[3],median_loc],dim=c(prod(dim(grid_output$parameters)[1:2]),(dim(grid_output$parameters)[3])))
  image.plot(cor(tmp,use="na.or.complete", method="spearman"), main="Median Parameter Correlation (Spearmans)")
  dev.off()

  ###
  ## Some things that are really parameters but are dependent on the evolution of state variables so need more information to be calculated
  # assign this value once as aall arrays have the same number of values present
  colour_choices=colour_choices_upper(prod(c(PROJECT$long_dim,PROJECT$lat_dim)))

  # create map of NPP allocation fractions
  if (exists(x = "NPP_foliage_fraction", where = grid_output)) {
      # Plot median
      figname = paste("parameter_maps_median_","Cfoliar_NPP_alloc","_",gsub("%","_",PROJECT$name),".jpeg",sep="")
      jpeg(file=figname, width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
      zrange=range(pretty(c(min(grid_output$NPP_foliage_fraction[,,median_loc], na.rm=TRUE),max(grid_output$NPP_foliage_fraction[,,median_loc],na.rm=TRUE))))
      image.plot(x = grid_long, y = grid_lat, z = grid_output$NPP_foliage_fraction[,,median_loc], col=colour_choices
                ,main=paste("NPP fractional allocation to foliage",sep=""), axes=FALSE, zlim=zrange
                ,cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
      map(add=TRUE, lwd = 2)
      #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
      dev.off()
      # Plot uncertainty
      figname = paste("parameter_maps_95CI_","Cfoliar_NPP_alloc","_",gsub("%","_",PROJECT$name),".jpeg",sep="")
      jpeg(file=figname, width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
      var = grid_output$NPP_foliage_fraction[,,upper_loc]-grid_output$NPP_foliage_fraction[,,lower_loc]
      zrange=range(pretty(c(min(var, na.rm=TRUE),max(var,na.rm=TRUE))))
      image.plot(x = grid_long, y = grid_lat, z = grid_output$NPP_foliage_fraction[,,upper_loc]-grid_output$NPP_foliage_fraction[,,lower_loc]
                ,col=colour_choices, main=paste("Uncertainty on NPP fractional allocation to foliage",sep=""),axes=FALSE
                ,cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1), zlim=zrange)
      map(add=TRUE, lwd = 2)
      #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
      dev.off()
  }
  if (exists(x = "NPP_wood_fraction", where = grid_output)) {
      # Plot median
      figname = paste("parameter_maps_median_","Cwood_NPP_alloc","_",gsub("%","_",PROJECT$name),".jpeg",sep="")
      jpeg(file=figname, width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
      var = grid_output$NPP_wood_fraction[,,median_loc]
      zrange=range(pretty(c(min(var, na.rm=TRUE),max(var,na.rm=TRUE))))
      image.plot(x = grid_long, y = grid_lat, z = var,col=colour_choices
                ,main=paste("NPP fractional allocation to wood",sep=""), axes=FALSE
                ,cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1), zlim=zrange)
      map(add=TRUE, lwd = 2)
      #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
      dev.off()
      # plot uncertainty
      figname = paste("parameter_maps_95CI_","Cwood_NPP_alloc","_",gsub("%","_",PROJECT$name),".jpeg",sep="")
      jpeg(file=figname, width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
      var = grid_output$NPP_wood_fraction[,,upper_loc]-grid_output$NPP_wood_fraction[,,lower_loc]
      zrange=range(pretty(c(min(var, na.rm=TRUE),max(var,na.rm=TRUE))))
      image.plot(x = grid_long, y = grid_lat, z = var
                ,col=colour_choices, main=paste("Uncertainty on NPP fractional allocation to wood",sep=""),axes=FALSE
                ,cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1), zlim=zrange)
      map(add=TRUE, lwd = 2)
      #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
      dev.off()
  }
  if (exists(x = "NPP_roots_fraction", where = grid_output)) {
      # Plot median
      figname = paste("parameter_maps_median_","Croot_NPP_alloc","_",gsub("%","_",PROJECT$name),".jpeg",sep="")
      jpeg(file=figname, width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
      var = grid_output$NPP_roots_fraction[,,median_loc]
      zrange=range(pretty(c(min(var, na.rm=TRUE),max(var,na.rm=TRUE))))
      image.plot(x = grid_long, y = grid_lat, z = var, col=colour_choices
                ,main=paste("NPP fractional allocation to fine roots",sep=""), axes=FALSE
                ,cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1), zlim=zrange)
      map(add=TRUE, lwd = 2)
      #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
      dev.off()
      # Plot uncertainty
      figname = paste("parameter_maps_95CI_","Croot_NPP_alloc","_",gsub("%","_",PROJECT$name),".jpeg",sep="")
      jpeg(file=figname, width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
      var = grid_output$NPP_roots_fraction[,,upper_loc]-grid_output$NPP_roots_fraction[,,lower_loc]
      zrange=range(pretty(c(min(var, na.rm=TRUE),max(var,na.rm=TRUE))))
      image.plot(x = grid_long, y = grid_lat, z = var
                ,col=colour_choices, main=paste("Uncertainty on NPP fractional allocation to fine roots",sep=""),axes=FALSE
                ,cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1), zlim=zrange)
      map(add=TRUE, lwd = 2)
      #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
      dev.off()
  }

  # create map of residence times
  if (exists(x = "MTT_foliage_years", where = grid_output)) {
      # Plot medians
      jpeg(file=paste("parameter_maps_median_","Cfoliar_residence_time","_",gsub("%","_",PROJECT$name),".jpeg",sep=""),
           width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
      var = grid_output$MTT_foliage_years[,,median_loc]
      zrange=range(pretty(c(min(var, na.rm=TRUE),max(var,na.rm=TRUE))))
      image.plot(x = grid_long, y = grid_lat, z = var, col=colour_choices
                ,main=paste("Mean foliage transit time (years)",sep="")
                ,axes=FALSE, cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1), zlim=zrange)
      map(add=TRUE, lwd = 2)
      #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
      dev.off()
      # Plot uncertainty
      jpeg(file=paste("parameter_maps_median_uncertainty_","Cfoliar_residence_time","_",gsub("%","_",PROJECT$name),".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
      var = grid_output$MTT_foliage_years[,,upper_loc]-grid_output$MTT_foliage_years[,,lower_loc]
      zrange=range(pretty(c(min(var, na.rm=TRUE),max(var,na.rm=TRUE))))
      image.plot(x = grid_long, y = grid_lat, z = var
                ,col=colour_choices, main=paste("Uncertainty on foliage transit time (years)",sep="")
                ,axes=FALSE, cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1), zlim=zrange)
      map(add=TRUE, lwd = 2)
      #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
      dev.off()
  }

  if (exists(x = "MTT_wood_years", where = grid_output)) {
      # Plot median
      jpeg(file=paste("parameter_maps_median_","Cwood_residence_time","_",gsub("%","_",PROJECT$name),".jpeg",sep=""),
           width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
      var = grid_output$MTT_wood_years[,,median_loc]
      zrange=range(pretty(c(min(var, na.rm=TRUE),max(var,na.rm=TRUE))))
      image.plot(x = grid_long, y = grid_lat, z = var, col=colour_choices
                ,main=paste("Mean wood transit time (years)",sep="")
                ,axes=FALSE, cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1), zlim=zrange)
      map(add=TRUE, lwd = 2)
      #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
      dev.off()
      # Plot uncertainty
      jpeg(file=paste("parameter_maps_median_uncertainty_","Cwood_residence_time","_",gsub("%","_",PROJECT$name),".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
      var = grid_output$MTT_wood_years[,,upper_loc]-grid_output$MTT_wood_years[,,lower_loc]
      zrange=range(pretty(c(min(var, na.rm=TRUE),max(var,na.rm=TRUE))))
      image.plot(x = grid_long, y = grid_lat, z = var
               ,col=colour_choices, main=paste("Uncertainty on wood transit time (years)",sep="")
               ,axes=FALSE, cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1), zlim=zrange)
      map(add=TRUE, lwd = 2)
      #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
      dev.off()
  }

  if (exists(x = "MTT_roots_years", where = grid_output)) {
      # Plot median
      jpeg(file=paste("parameter_maps_median_","Croot_residence_time","_",gsub("%","_",PROJECT$name),".jpeg",sep=""),
           width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
      var = grid_output$MTT_roots_years[,,median_loc]
      zrange=range(pretty(c(min(var, na.rm=TRUE),max(var,na.rm=TRUE))))
      image.plot(x = grid_long, y = grid_lat, z = var, col=colour_choices
                ,main=paste("Mean fine roots transit time (years)",sep="")
                ,axes=FALSE, cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1), zlim=zrange)
      map(add=TRUE, lwd = 2)
      #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
      dev.off()
      # Plot uncertainty
      jpeg(file=paste("parameter_maps_median_uncertainty_","Croot_residence_time","_",gsub("%","_",PROJECT$name),".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
      var = grid_output$MTT_roots_years[,,upper_loc]-grid_output$MTT_roots_years[,,lower_loc]
      zrange=range(pretty(c(min(var, na.rm=TRUE),max(var,na.rm=TRUE))))
      image.plot(x = grid_long, y = grid_lat, z = var
                ,col=colour_choices, main=paste("Uncertainty on fine roots transit time (years)",sep="")
                ,axes=FALSE, cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1), zlim=zrange)
      map(add=TRUE, lwd = 2)
      #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
      dev.off()
  }

  if (exists(x = "MTT_som_years", where = grid_output)) {
      # Plot medians
      jpeg(file=paste("parameter_maps_median_","Csom_residence_time","_",gsub("%","_",PROJECT$name),".jpeg",sep=""),
           width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
      var = grid_output$MTT_som_years[,,median_loc]
      zrange=range(pretty(c(min(var, na.rm=TRUE),max(var,na.rm=TRUE))))
      image.plot(x = grid_long, y = grid_lat, z = var, col=colour_choices
                ,main=paste("Csom residence time median estimate",sep="")
                ,axes=FALSE, cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1), zlim=zrange)
      map(add=TRUE, lwd = 2)
      #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
      dev.off()
      # Plot uncertainty
      jpeg(file=paste("parameter_maps_median_uncertainty_","Csom_residence_time","_",gsub("%","_",PROJECT$name),".jpeg",sep=""),
           width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
      var = grid_output$MTT_som_years[,,upper_loc]-grid_output$MTT_som_years[,,lower_loc]
      zrange=range(pretty(c(min(var, na.rm=TRUE),max(var,na.rm=TRUE))))
      image.plot(x = grid_long, y = grid_lat, z = var
                ,col=colour_choices, main=paste("Uncertainty on som transit time (years)",sep="")
                ,axes=FALSE, cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1), zlim=zrange)
      map(add=TRUE, lwd = 2)
      #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
      dev.off()
  }

  if (exists(x = "MTT_litter_years", where = grid_output)) {
      # Plot median
      jpeg(file=paste("parameter_maps_median_","Clitter_residence_time","_",gsub("%","_",PROJECT$name),".jpeg",sep=""),
           width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
      var = grid_output$MTT_litter_years[,,median_loc]
      zrange=range(pretty(c(min(var, na.rm=TRUE),max(var,na.rm=TRUE))))
      image.plot(x = grid_long, y = grid_lat, z = var, col=colour_choices
                ,main=paste("Clitter residence time median estimate",sep="")
                ,axes=FALSE, cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1), zlim=zrange)
      map(add=TRUE, lwd = 2)
      #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
      dev.off()
      # Plot uncertainty
      jpeg(file=paste("parameter_maps_median_uncertainty_","Clitter_residence_time","_",gsub("%","_",PROJECT$name),".jpeg",sep=""),
           width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
      var = grid_output$MTT_litter_years[,,upper_loc]-grid_output$MTT_litter_years[,,lower_loc]
      zrange=range(pretty(c(min(var, na.rm=TRUE),max(var,na.rm=TRUE))))
      image.plot(x = grid_long, y = grid_lat, z = var
                ,col=colour_choices, main=paste("Uncertainty on litter transit time (years)",sep="")
                ,axes=FALSE, cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1), zlim=zrange)
      map(add=TRUE, lwd = 2)
      #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
      dev.off()
  }

  if (exists(x = "MTT_woodlitter_years", where = grid_output)) {
      # Plot medians
      jpeg(file=paste("parameter_maps_median_","Cwoodlitter_residence_time","_",gsub("%","_",PROJECT$name),".jpeg",sep=""),
           width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
      var = grid_output$MTT_woodlitter_years[,,median_loc]
      zrange=range(pretty(c(min(var, na.rm=TRUE),max(var,na.rm=TRUE))))
      image.plot(x = grid_long, y = grid_lat, z = var
                ,main=paste("Cwoodlitter residence time median estimate",sep="")
                ,axes=FALSE, cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1), zlim=zrange)
      map(add=TRUE, lwd = 2)
      #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
      dev.off()
      # Plot uncertainty
      jpeg(file=paste("parameter_maps_median_uncertainty_","Cwoodlitter_residence_time","_",gsub("%","_",PROJECT$name),".jpeg",sep=""),
           width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
      var = grid_output$MTT_woodlitter_years[,,upper_loc]-grid_output$MTT_woodlitter_years[,,lower_loc]
      zrange=range(pretty(c(min(var, na.rm=TRUE),max(var,na.rm=TRUE))))
      image.plot(x = grid_long, y = grid_lat, z = var
                ,col=colour_choices, main=paste("Uncertainty on wood litter transit time (years)",sep="")
                ,axes=FALSE, cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1), zlim=zrange)
      map(add=TRUE, lwd = 2)
      #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
      dev.off()
  }

  # now everything has been loaded into a nice array we begin plotting
  for (p in seq(1,dim(grid_output$parameters)[3])) {
       ###
       ## Histrograms of parameters
       jpeg(file=paste("parameter_hist_median_",par_names[p],"_",gsub("%","_",PROJECT$name),".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
       par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
       hist(grid_output$parameters[,,p,median_loc], main=paste(par_names[p]," median estimate",sep=""), cex.main=1.1, cex=1.5, cex.axis=1.8)
       dev.off()
       ###
       ## spatial maps
       jpeg(file=paste("parameter_maps_median_",par_names[p],"_",gsub("%","_",PROJECT$name),".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
       par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
       var = grid_output$parameters[,,p,median_loc]
       zrange=range(pretty(c(min(var, na.rm=TRUE),max(var,na.rm=TRUE))))
       image.plot(x = grid_long, y = grid_lat, z = var, col=colour_choices, zlim=zrange,
                  main=paste(par_names[p]," median estimate",sep=""),axes=FALSE, cex.main=1.1,legend.width=3.0,
                  cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
       map(add=TRUE, lwd = 2)
       #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
       dev.off()
       jpeg(file=paste("parameter_maps_95CI_",par_names[p],"_",gsub("%","_",PROJECT$name),".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
       par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
       var = grid_output$parameters[,,p,upper_loc]-grid_output$parameters[,,p,lower_loc]
       zrange=range(pretty(c(min(var, na.rm=TRUE),max(var,na.rm=TRUE))))
       image.plot(x = grid_long, y = grid_lat, z = var,
                  col=colour_choices, main=paste(par_names[p]," uncertainty range",sep=""), axes=FALSE, zlim=zrange,
                  cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
       map(add=TRUE, lwd = 2)
       #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
       dev.off()
       jpeg(file=paste("parameter_maps_converged_",par_names[p],"_",gsub("%","_",PROJECT$name),".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
       par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
       image.plot(x = grid_long, y = grid_lat, z = grid_output$parameters_converged[,,p], col=colour_choices
                 ,main=paste(par_names[p]," Gelmen-Rubens convergence (1 = PASS / 0 = FALSE)",sep="")
                 ,axes=FALSE, zlim=c(0,1), cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
       map(add=TRUE, lwd = 2)
       #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
       dev.off()
  } # parameter loop

  # mean temperature
  jpeg(file=paste("mean_temperature_maps_",gsub("%","_",PROJECT$name),".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
  image.plot(x = grid_long, y = grid_lat, z = grid_output$mean_temperature_C, col=rev(colour_choices)
            ,main="Mean air temperature (oC)",axes=FALSE, cex.main=1.1,legend.width=3.0
            ,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1)
            ,zlim=c(min(as.vector(grid_output$mean_temperature_C),na.rm=TRUE),max(as.vector(grid_output$mean_temperature_C),na.rm=TRUE)))
  map(add=TRUE, lwd = 2)
  #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
  dev.off()
  # mean radiation
  jpeg(file=paste("mean_radiation_maps_",gsub("%","_",PROJECT$name),".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
  image.plot(x = grid_long, y = grid_lat, z = grid_output$mean_radiation_MJm2day, col=rev(colour_choices)
            ,main="Mean radiation (MJ/m2/day)"
            ,axes=FALSE, cex.main=1.1,legend.width=3.0
            ,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1),zlim=c(0,36))
  map(add=TRUE, lwd = 2)
  #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
  dev.off()
  # mean vpd
  jpeg(file=paste("mean_vpd_maps_",gsub("%","_",PROJECT$name),".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
  image.plot(x = grid_long, y = grid_lat, z = grid_output$mean_vpd_Pa, col=rev(colour_choices)
            ,main="Mean VPD (Pa)"
            ,axes=FALSE, cex.main=1.1, legend.width=3.0, cex=1.5
            ,axis.args=list(cex.axis=1.8,hadj=0.1),zlim=c(0,max(as.vector(grid_output$mean_vpd_Pa),na.rm=TRUE)))
  map(add=TRUE, lwd = 2)
  #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
  dev.off()
  # mean precipitation
  jpeg(file=paste("mean_precipitation_maps_",gsub("%","_",PROJECT$name),".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
  par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
  image.plot(x = grid_long, y = grid_lat, z = grid_output$mean_precipitation_kgH2Om2yr, col=colour_choices,
             main="Mean precipitation (kgH2O/m2/yr)",axes=FALSE, cex.main=1.1,legend.width=3.0
            ,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1),zlim=c(0,max(as.vector(grid_output$mean_precipitation_kgH2Om2yr),na.rm=TRUE)))
  map(add=TRUE, lwd = 2)
  #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
  dev.off()

  # tidy before leaving
  gc(reset=TRUE, verbose=FALSE)

  # Move working directory
  setwd(old_wd)

  print("...done plotting parameter / observation / driver maps...")

  } # end function generate_parameter_maps

  ## Use byte compile
  generate_parameter_maps<-cmpfun(generate_parameter_maps)
