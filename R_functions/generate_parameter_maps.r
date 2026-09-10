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
# Function to create generic plots of gridded CARDAMOM analysis parameters and emergent traits
# 
# Author: T. Luke Smallman (12/11/2024)
#
#########################################################################################

generate_parameter_maps<-function(PROJECT) {

   # Move working directory
   old_wd = getwd() ; setwd(PROJECT$figpath)

   # Which quantiles will we extract
   na_flag = TRUE

   # Loaded the grid aggregated dataset into memory
   infile = paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep="")
   load(infile)

   # Determine the correct quantiles
   ee = length(grid_output$num_quantiles) ; ss = 1 ; ci95 = rep(NA, 2) ; ci68 = rep(NA, 2)    
   for (i in seq(1, floor(length(grid_output$num_quantiles)*0.5))) {
        tmp1 = grid_output$num_quantiles[ss] ; tmp2 = grid_output$num_quantiles[ee]         
        if (round(tmp2-tmp1, digits=2) == 0.95) { ci95[1] = ss ; ci95[2] = ee }
        if (round(tmp2-tmp1, digits=2) == 0.68) { ci68[1] = ss ; ci68[2] = ee }
        ss = ss + 1 ; ee = ee - 1
   }
   # Median estimate array index
   median_loc = which(round(grid_output$num_quantiles, digits=2) == 0.5)
   # Lower and lower quantiles array index
   #low_quant = ci95[1]  ; high_quant = ci95[2] # 95 %
   lower_loc = ci68[1] ; upper_loc = ci68[2] # 68, i.e. 1-SD %

   # Extract the lat / long information
   grid_lat = grid_output$lat ; grid_long = grid_output$long

   # Ensure that any Inf values are removed to NA
   # This is a hack implemented during a development phase and should be removed
   # or commented out under general operation
   for (i in seq(1, length(grid_output))) {
        grid_output[[i]][which(is.infinite(grid_output[[i]]))] = NA
   }

   # load timesteps to local variable
   timestep_days = PROJECT$model$timestep_days ; seconds_per_day = 86400
   timestep_days = rep(timestep_days, length.out = grid_output$time_dim)

   # Convert avgN log10-normal to gN/m2, in models which use foliar N in gN/m2
   if (PROJECT$model$name == "DALEC.A1.C2.D2.F2.H2.P3.R1.008"| PROJECT$model$name == "DALEC.A1.C2.D2.F2.H1.P3.R1.009" |
       PROJECT$model$name == "DALEC.A1.C2.D2.F2.H2.P4.R2.010" | 
       PROJECT$model$name == "DALEC.A1.C2.D2.F2.H1.P4.R2.011" | PROJECT$model$name == "DALEC.A1.C2.D2.F2.H2.P7.R2.023" |
       PROJECT$model$name == "DALEC.A1.C2.D2.F2.H2.P8.R2.024" | PROJECT$model$name == "DALEC.A1.C2.D2.F2.H2.P10.R2.026") {
       grid_output$parameters[,,11,] = 10**grid_output$parameters[,,11,]
   }

  # Determine correct height and widths
  hist_height = 4000*0.65 ; hist_width = 7200*0.65
  fig_height = 3000*0.65 ; fig_width = ((PROJECT$long_dim/PROJECT$lat_dim)+0.25) * fig_height
  if (grepl("27700",PROJECT$grid_type)) { fig_height = 8000*0.65 ; fig_width = 7200*0.65 }
  # Set up colour schemes
  colour_choices_default = colorRampPalette(brewer.pal(11,"Spectral")) 
  colour_choices_sign    = colorRampPalette(brewer.pal(11,"PRGn"))
  colour_choices_gain    = colorRampPalette(brewer.pal(9,"YlGnBu"))
  colour_choices_loss    = colorRampPalette(brewer.pal(9,"YlOrRd"))
  colour_choices_CI      = colorRampPalette(brewer.pal(9,"Purples"))
  # Now extract out the final colours we will use
  colour_choices_default = colour_choices_default(100)
  colour_choices_sign    = colour_choices_sign(100)
  colour_choices_gain    = colour_choices_gain(100)
  colour_choices_loss    = colour_choices_loss(100)
  colour_choices_CI      = colour_choices_CI(100)


  # inform the user
  print("......have finished loading - now beginning cluster analysis")

  if (exists(x = "clusters", where = grid_output) == FALSE | repair == 1) {

      if (grepl("DALEC",PROJECT$model$name)) {

          ###
          ## Begin Affinity propagation clustering - for process parameters + initial conditions

          # Now generate cluster map using all parameters including the initial conditions
          par_array_median_normalised = grid_output$parameters[,,,median_loc]
          # now normalise the parameter values
          for (i in seq(1,dim(par_array_median_normalised)[3])) {
               min_par_val = min(as.vector(par_array_median_normalised[,,i]),na.rm=TRUE)
               max_par_val = max(as.vector(par_array_median_normalised[,,i]),na.rm=TRUE)
               par_array_median_normalised[,,i] = ((par_array_median_normalised[,,i]-min_par_val)/(max_par_val-min_par_val))
          }
          # Create temporary arrays needed to allow removing of NAs and convert array into (space,par)
          par_array_tmp = array(NA,dim=c(prod(dim(grid_output$parameters)[1:2]),dim(par_array_median_normalised)[3]))
          par_array_tmp[1:prod(dim(par_array_median_normalised)[1:2]),1:dim(par_array_median_normalised)[3]] = par_array_median_normalised
          actual_forests = which(is.na(par_array_tmp[,1]) == FALSE)
          par_array_tmp = par_array_tmp[actual_forests,]
          par_array_tmp = array(par_array_tmp,dim=c((length(par_array_tmp)/dim(par_array_median_normalised)[3]),dim(par_array_median_normalised)[3]))

          # Determine the number of pixels to be in the subsample.
          # 2000 sites is an arbitary selection.
          subsample_frac = min(0.1,2000/PROJECT$nosites)
          # Looping to find preference_input, the preferenceRange() returns 2 values,
          # the first of which minimises the number of clusters, while the seconds
          # would return as many clusters as there are observations. It is the
          # responsibility of the user to ensure the most appropriate use of these
          # information to result in an appropriate number of clusters for error propagation
          tmp = 0 
          for (i in seq(1,3)) {
               tmp = append(tmp,preferenceRange(negDistMat(par_array_tmp[sample(1:dim(par_array_tmp)[1],subsample_frac*dim(par_array_tmp)[1], replace=FALSE),],r=2))[1])
          } ; preference_input = max(mean(tmp[-1]), median(tmp[-1]))
          #print(paste("q = ",preference_input))
          # Conduct the affinity propagation clustering analysis.
          # convits = how many iterations to wait for a change in clustering before cancelling process.
          # maxits  = maximum number of iteration to perform regardless or clustering still changing
          # frac    = if set, clustering is done on a randomly selected sub-sample of pixels. The sub-sample is the fraction specified.
          # sweeps  = if frac set, how many times to repeat the subsampling processes
          grid_output$cluster_analysis = apclusterL(negDistMat(r=2), par_array_tmp, frac = subsample_frac, sweeps = 5, p = preference_input, maxits=500, convits=50)
          #grid_output$cluster_analysis=apclusterL(negDistMat(r=2),par_array_tmp,frac=0.1, sweeps=10, q=0.05, maxits=1000, convits=100)
          grid_output$nos_clusters = length(grid_output$cluster_analysis@clusters) ; grid_output$clusters_exemplars=grid_output$cluster_analysis@exemplars
          grid_output$clusters = array(NA,dim=c(dim(par_array_median_normalised)[1:2]))
          for (i in seq(1,length(grid_output$cluster_analysis@clusters))) {
               grid_output$clusters[actual_forests[grid_output$cluster_analysis@clusters[[i]]]] = i
          }
          grid_output$clusters=array(grid_output$clusters,dim=c(dim(par_array_median_normalised)[1:2]))
          # Tidy away the overall analysis in faviour of what we have extracted
          grid_output = within(grid_output, rm(cluster_analysis))

          # Save the now updated cluster analysis into the grid_output
          save(grid_output, file = infile)

#          # Now plot both possible cluster maps
#          figname = paste("Cluster_map_of_median_parameters_",gsub("%","_",PROJECT$name),".jpeg",sep="")
#          jpeg(file=figname, width=fig_width, height=fig_height, res=300, quality=100)
#          par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
#          image.plot(x = grid_long, y = grid_lat, z = grid_output$pars_clusters, main=paste("Parameter based cluster maps",sep=""),axes=FALSE, cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
#          map(add=TRUE, lwd = 2)
#          dev.off()

          figname = paste("Cluster_map_of_median_parameters_with_initial_",gsub("%","_",PROJECT$name),".jpeg",sep="")
          jpeg(file=figname, width=fig_width, height=fig_height, res=300, quality=100)
          par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
          image.plot(x = grid_long, y = grid_lat, z = grid_output$clusters, main=paste("Parameter + initial based cluster map",sep=""),
                     axes=FALSE, cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
          map(add=TRUE, lwd = 2)
          dev.off()

      } # make cluster map

  } # reprocessing or not?

  # inform the user
  print("......now generating parameter correlations matrix")

  # describe parameter numbers with log liklihood added on the end
  character_bit=rep("p",times=(PROJECT$model$nopars[1]))
  number_bit=1:(PROJECT$model$nopars[1])
  # merge the letter and numbers together
  par_names=c(paste(character_bit,number_bit,sep=""),"log-likelihood")

  # create tempory array of parameters for the covariance analysis
  tmp = array(grid_output$parameters[,,1:dim(grid_output$parameters)[3],median_loc],dim=c(prod(dim(grid_output$parameters)[1:2]),(dim(grid_output$parameters)[3])))
  tmp = cor(tmp,use="na.or.complete", method="spearman")
  if (any(is.finite(tmp))) {
      figname = paste("parameter_correlations_median_",gsub("%","_",PROJECT$name),".jpeg",sep="")
      jpeg(file=figname, width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
      image.plot(cor(tmp,use="na.or.complete", method="spearman"), main="Median Parameter Correlation (Spearmans)")
      dev.off()
  } else {
      print("WARNING: problem with the grid_output$parameters array as correlation analysis returned NaN or Inf")
  }

  ###
  ## Some things that are really parameters but are dependent on the evolution of state variables so need more information to be calculated
  ## assign this value once as all arrays have the same number of values present

  # inform the user
  print("......now generating NPP allocation figures")

  # create map of NPP allocation fractions
  if (exists(x = "NPP_foliage_fraction", where = grid_output)) {
      # Plot median
      figname = paste("parameter_maps_median_","Cfoliar_NPP_alloc","_",gsub("%","_",PROJECT$name),".jpeg",sep="")
      jpeg(file=figname, width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
      zrange=range(pretty(c(min(grid_output$NPP_foliage_fraction[,,median_loc], na.rm=TRUE),max(grid_output$NPP_foliage_fraction[,,median_loc],na.rm=TRUE))))
      image.plot(x = grid_long, y = grid_lat, z = grid_output$NPP_foliage_fraction[,,median_loc], col=colour_choices_gain
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
                ,col=colour_choices_CI, main=paste("Uncertainty on NPP fractional allocation to foliage",sep=""),axes=FALSE
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
      image.plot(x = grid_long, y = grid_lat, z = var,col=colour_choices_gain
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
                ,col=colour_choices_CI, main=paste("Uncertainty on NPP fractional allocation to wood",sep=""),axes=FALSE
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
      image.plot(x = grid_long, y = grid_lat, z = var, col=colour_choices_gain
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
                ,col=colour_choices_CI, main=paste("Uncertainty on NPP fractional allocation to fine roots",sep=""),axes=FALSE
                ,cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1), zlim=zrange)
      map(add=TRUE, lwd = 2)
      #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
      dev.off()
  }

  # inform the user
  print("......now generating mean transit time figures")

  # create map of residence times
  if (exists(x = "MTT_foliage_years", where = grid_output)) {
      # Plot medians
      jpeg(file=paste("parameter_maps_median_","Cfoliar_residence_time","_",gsub("%","_",PROJECT$name),".jpeg",sep=""),
           width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
      var = grid_output$MTT_foliage_years[,,median_loc]
      zrange=range(pretty(c(min(var, na.rm=TRUE),max(var,na.rm=TRUE))))
      image.plot(x = grid_long, y = grid_lat, z = var, col=colour_choices_gain
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
                ,col=colour_choices_CI, main=paste("Uncertainty on foliage transit time (years)",sep="")
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
      image.plot(x = grid_long, y = grid_lat, z = var, col=colour_choices_gain
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
               ,col=colour_choices_CI, main=paste("Uncertainty on wood transit time (years)",sep="")
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
      image.plot(x = grid_long, y = grid_lat, z = var, col=colour_choices_gain
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
                ,col=colour_choices_CI, main=paste("Uncertainty on fine roots transit time (years)",sep="")
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
      image.plot(x = grid_long, y = grid_lat, z = var, col=colour_choices_gain
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
                ,col=colour_choices_CI, main=paste("Uncertainty on som transit time (years)",sep="")
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
      image.plot(x = grid_long, y = grid_lat, z = var, col=colour_choices_gain
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
                ,col=colour_choices_CI, main=paste("Uncertainty on litter transit time (years)",sep="")
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
                ,col=colour_choices_gain, main=paste("Cwoodlitter residence time median estimate",sep="")
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
                ,col=colour_choices_CI, main=paste("Uncertainty on wood litter transit time (years)",sep="")
                ,axes=FALSE, cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1), zlim=zrange)
      map(add=TRUE, lwd = 2)
      #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
      dev.off()
  }

  # inform the user
  print("......now generating parameter figures")

  # now everything has been loaded into a nice array we begin plotting
  for (p in seq(1,dim(grid_output$parameters)[3])) {
       if (any(is.finite(as.vector(grid_output$parameters[,,p,median_loc])))) {
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
       zrange = range(pretty(c(min(var, na.rm=TRUE),max(var,na.rm=TRUE))))
       if (zrange[1] >= 0 & zrange[2] > 0) {
           colour_choices = colour_choices_gain
       } else if (zrange[1] < 0 & zrange[2] > 0) {
           colour_choices = colour_choices_sign
       } else if (zrange[1] < 0 & zrange[2] <= 0) {
           colour_choices = rev(colour_choices_loss)
       } else {
           colour_choices = colour_choices_default
       }
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
                  col=colour_choices_CI, main=paste(par_names[p]," uncertainty range",sep=""), axes=FALSE, zlim=zrange,
                  cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
       map(add=TRUE, lwd = 2)
       #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
       dev.off()
       jpeg(file=paste("parameter_maps_converged_",par_names[p],"_",gsub("%","_",PROJECT$name),".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
       par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
       image.plot(x = grid_long, y = grid_lat, z = grid_output$parameters_converged[,,p], col=colour_choices_sign
                 ,main=paste(par_names[p]," Gelmen-Rubens convergence (1 = PASS / 0 = FALSE)",sep="")
                 ,axes=FALSE, zlim=c(0,1), cex.main=1.1,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
       map(add=TRUE, lwd = 2)
       #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
       dev.off()
       }

       # Generate generic maps of spatial aggregates of drivers$parpriors
       if (p <= max(PROJECT$model$nopars)) {
           if (length(which(is.na(grid_output$parameter_priors_array[,,p]) != TRUE)) > 0) {
               zrange = c(min(as.vector(grid_output$parameter_priors_array[,,p]),na.rm=TRUE),max(as.vector(grid_output$parameter_priors_array[,,p]),na.rm=TRUE))
               zrange = sort(zrange + (c(-0.01,0.01) * zrange))
               if (zrange[1] >= 0 & zrange[2] > 0) {
                   colour_choices = colour_choices_gain
               } else if (zrange[1] < 0 & zrange[2] > 0) {
                   colour_choices = colour_choices_sign
               } else if (zrange[1] < 0 & zrange[2] <= 0) {
                   colour_choices = rev(colour_choices_loss)
               } else {
                   colour_choices = colour_choices_default
               }               
               fig_name = paste("mean_parameter_priors_array_maps_",gsub(" ","_",par_names[p]),"_",gsub("%","_",PROJECT$name),".jpeg",sep="")
               fig_name = gsub("\\(","", fig_name) ; fig_name = gsub("\\)","", fig_name)
               fig_name = gsub("/","", fig_name) ; fig_name = gsub("/","", fig_name)
               jpeg(file=fig_name, width=fig_width, height=fig_height, res=300, quality=100)
               par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
               image.plot(x = grid_long, y = grid_lat, z = grid_output$parameter_priors_array[,,p], col = colour_choices
                         ,main=paste(par_names[p],sep=""),axes=FALSE, cex.main=1.1,legend.width=3.0
                         ,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1)
                         ,zlim=zrange)
               map(add=TRUE, lwd = 2)
               #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
               dev.off()
               zrange = c(min(as.vector(grid_output$parameter_priors_uncertainty_array[,,p]),na.rm=TRUE),max(as.vector(grid_output$parameter_priors_uncertainty_array[,,p]),na.rm=TRUE))
               zrange = sort(zrange + (c(-0.01,0.01) * zrange))
               fig_name = paste("mean_parameter_priors_uncertainty_array_maps_",gsub(" ","_",par_names[p]),"_",gsub("%","_",PROJECT$name),".jpeg",sep="")
               fig_name = gsub("\\(","", fig_name) ; fig_name = gsub("\\)","", fig_name)
               fig_name = gsub("/","", fig_name) ; fig_name = gsub("/","", fig_name)
               jpeg(file=fig_name, width=fig_width, height=fig_height, res=300, quality=100)
               par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
               image.plot(x = grid_long, y = grid_lat, z = grid_output$parameter_priors_uncertainty_array[,,p], col=colour_choices_CI
                         ,main=paste(par_names[p],sep=""),axes=FALSE, cex.main=1.1,legend.width=3.0
                         ,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1)
                         ,zlim=zrange)
               map(add=TRUE, lwd = 2)
               #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
               dev.off()
           } # has a value worth reporting
       } # Is a parameter not log-likelihood

  } # parameter loop
  
  # inform the user
  print("......now generating forcings figures")

  # Generate generic maps of spatial aggregates of drivers$met
  for (m in seq(1, length(met_array_names))) {
       zrange = c(min(as.vector(grid_output$met_array_averages[,,m]),na.rm=TRUE),max(as.vector(grid_output$met_array_averages[,,m]),na.rm=TRUE))
       zrange = zrange + (c(-0.05,0.05) * zrange)
       if (diff(zrange) == 0) {zrange = c(-0.05,0.05)}
       zrange = sort(zrange)
       fig_name = paste("mean_met_array_maps_",gsub(" ","_",met_array_names[m]),"_",gsub("%","_",PROJECT$name),".jpeg",sep="")
       fig_name = gsub("\\(","", fig_name) ; fig_name = gsub("\\)","", fig_name)
       fig_name = gsub("/","", fig_name) ; fig_name = gsub("/","", fig_name)
       jpeg(file=fig_name, width=fig_width, height=fig_height, res=300, quality=100)
       par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
       image.plot(x = grid_long, y = grid_lat, z = grid_output$met_array_averages[,,m], col=rev(colour_choices_default)
                 ,main=met_array_names[m],axes=FALSE, cex.main=1.1,legend.width=3.0
                 ,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1)
                 ,zlim=zrange)
       map(add=TRUE, lwd = 2)
       #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
       dev.off()
  }

  # inform the user
  print("......now generating assimilated observations figures")

  # Generate generic maps of spatial aggregates of drivers$obs
  for (m in seq(1, length(obs_array_names))) {
       if (any(is.na(grid_output$obs_array_averages[,,m]) == FALSE)) {
           zrange = c(min(as.vector(grid_output$obs_array_averages[,,m]),na.rm=TRUE),max(as.vector(grid_output$obs_array_averages[,,m]),na.rm=TRUE))
           zrange = zrange + (c(-0.05,0.05) * zrange)     
           if (diff(zrange) == 0) {zrange = c(-0.05,0.05)}           
           fig_name = paste("mean_obs_array_maps_",gsub(" ","_",obs_array_names[m]),"_",gsub("%","_",PROJECT$name),".jpeg",sep="")
           fig_name = gsub("\\(","", fig_name) ; fig_name = gsub("\\)","", fig_name)
           fig_name = gsub("/","", fig_name) ; fig_name = gsub("/","", fig_name)
           jpeg(file=fig_name, width=fig_width, height=fig_height, res=300, quality=100)
           par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
           image.plot(x = grid_long, y = grid_lat, z = grid_output$obs_array_averages[,,m], col=rev(colour_choices_default)
                     ,main=obs_array_names[m],axes=FALSE, cex.main=1.1,legend.width=3.0
                     ,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1)
                     ,zlim=zrange)
           map(add=TRUE, lwd = 2)
           #contour(grid_output$landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
           dev.off()
       }
  }

  # tidy before leaving
  gc(reset=TRUE, verbose=FALSE)

  # Move working directory
  setwd(old_wd)

  print("...done plotting parameter / observation / driver maps...")

  } # end function generate_parameter_maps

  ## Use byte compile
  generate_parameter_maps<-cmpfun(generate_parameter_maps)
