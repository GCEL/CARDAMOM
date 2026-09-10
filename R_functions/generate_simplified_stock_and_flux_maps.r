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
# Function to create generic stock and flux plots for gridded CARDAMOM analyses
# 
# Author: T. Luke Smallman (12/11/2024)
#
#########################################################################################

generate_simplified_stock_and_flux_maps<-function(PROJECT) {

  print("...beginning stocks and fluxes...")

  # Move working directory
  old_wd = getwd()
  setwd(PROJECT$figpath)

  # read in the grid stock information
  infile = paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep="")
  if (file.exists(infile) == FALSE) {stop("grid_outputs for 'generate_simplified_stock_and_flux_maps' missing")}
  load(paste(infile))
   
  # generate the lat / long grid again
  output = generate_grid(PROJECT$grid_type,PROJECT$latitude,PROJECT$longitude,PROJECT$resolution)
  area = output$area ; grid_lat = output$lat ; grid_long = output$long
  # include adjustment for g-> Tg (*1e-12)  
  area_with_g_Tg = area*1e-12
  rm(output)

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

  # calculate land mask
  landmask = array(PROJECT$landsea, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
  # Set up colour scheme
  smoothScatter_colours  = colorRampPalette(c("white",rep(rev(brewer.pal(11,"Spectral")),each=3)))
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

  # determine correct height and widths
  fig_height = 3000*0.65 ; fig_width = ((PROJECT$long_dim/PROJECT$lat_dim)+0.25) * fig_height
  if (grepl("27700",PROJECT$grid_type)) { fig_height = 8000*0.65 ; fig_width = 7200*0.65 }

  # If root depth information has been provided plot some of it up here.
  if (exists(x = "mean_RootDepth_m", where = grid_output)) {
      jpeg(file=paste("median_root_depth_maps_",gsub("%","_",PROJECT$name),".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.6), omi=c(0.2, 0.2, 0.2, 0.40))
      z_axis = c(min(as.vector(grid_output$mean_RootDepth_m),na.rm=TRUE),max(as.vector(grid_output$mean_RootDepth_m),na.rm=TRUE))
      image.plot(x = grid_long, y = grid_lat, z = grid_output$mean_RootDepth_m[,,median_loc],col=colour_choices_gain
                ,main=paste("Median root depth (m)",sep=""),zlim=z_axis,axes=FALSE
                ,cex.main=0.9,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
      map(add=TRUE, lwd = 2)
      #contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
      dev.off()
  } # plot rooting depth

  # extract a list of all the variables stored in the output object
  par_names = names(grid_output)
  # filter for those related to the 'mean' status
  par_names = par_names[grepl("mean", par_names)]
  # Special case removal for mean annual variables which will not be gridded
  par_names = par_names[grepl("mean_annual", par_names) == FALSE]

  # loop through these mean variables, output median and CI range for each of these variables
  for (p in seq(1,length(par_names))) {

       # determine position in grid_output list which contains the variable of interest
       pp = which(names(grid_output) == par_names[p])
       # Check that the number of dimensions matches that expected.
       # Strictly speaking this is a bit of a hack to account for not being able to identify the model variables only.
       # Could store non-model output variables into a different array to solve this.
       if (length(dim(grid_output[[pp]])) == 3) {
           # create maps
           jpeg(file=paste("grid_mean_map_",par_names[p],"_",gsub("%","_",PROJECT$name),".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
           par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.8), omi=c(0.2, 0.2, 0.2, 0.40))
           info = " " # assume default is no header, but sometimes we add something extra...
           var1 = mean(grid_output[[pp]][,,median_loc], na.rm=TRUE)
           var2 = mean(grid_output[[pp]][,,upper_loc], na.rm=TRUE)
           var3 = mean(grid_output[[pp]][,,lower_loc], na.rm=TRUE)
           var1 = round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2) 
           info = paste("Mean estimate: ",par_names[p]," (97.5 % = ",var2,"; 50 % = ",var1,"; 2.5 % = ",var3,")", sep="")
           zrange = range(pretty(c(min(grid_output[[pp]][,,median_loc], na.rm=TRUE),max(grid_output[[pp]][,,median_loc],na.rm=TRUE))))
           if (zrange[1] >= 0 & zrange[2] > 0) {
               colour_choices = colour_choices_gain
           } else if (zrange[1] < 0 & zrange[2] > 0) {
               colour_choices = colour_choices_sign
           } else if (zrange[1] < 0 & zrange[2] <= 0) {
               colour_choices = rev(colour_choices_loss)
           } else {
               colour_choices = colour_choices_default
           }
           image.plot(x = grid_long, y = grid_lat, z = grid_output[[pp]][,,median_loc], zlim=zrange, main=info, col = colour_choices,
                      axes=FALSE, cex.main=0.9, legend.width=3.0, cex=1.5, axis.args=list(cex.axis=1.8, hadj=0.1))
           map(add=TRUE, lwd = 2)
           #contour(landmask, add = TRUE, lwd = 1.0, nlevels = 1,axes = FALSE,drawlabels = FALSE,col = "black")
           dev.off()

           # Histrograms of fluxes
           jpeg(file=paste("grid_mean_hist_median_",par_names[p],"_",gsub("%","_",PROJECT$name),".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
           par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
           hist(grid_output[[pp]][,,median_loc], main=info, cex.main=0.9, cex=1.5, cex.axis=1.8)
           dev.off()

           # create maps
           jpeg(file=paste("grid_95CI_map_",par_names[p],"_",gsub("%","_",PROJECT$name),".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
           par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.8), omi=c(0.2, 0.2, 0.2, 0.40))
           info = " " # assume default is no header, but sometimes we add something extra...
           var2 = grid_output[[pp]][,,upper_loc] - grid_output[[pp]][,,lower_loc]
           var1 = mean(var2, na.rm=TRUE)
           var1 = round(var1,digit=2) 
           info = paste("Mean 95CI estimate: ",par_names[p]," = ",var1, sep="")
           zrange = range(pretty(c(min(var2, na.rm=TRUE),max(var2,na.rm=TRUE))))
           image.plot(x = grid_long, y = grid_lat, z = var2, zlim=zrange, main=info, col = colour_choices_CI,
                      axes=FALSE, cex.main=0.9, legend.width=3.0, cex=1.5, axis.args=list(cex.axis=1.8, hadj=0.1))
           map(add=TRUE, lwd = 2)
           #contour(landmask, add = TRUE, lwd = 1.0, nlevels = 1,axes = FALSE,drawlabels = FALSE,col = "black")
           dev.off()

           # Histrograms of fluxes
           jpeg(file=paste("grid_95CI_hist_median_",par_names[p],"_",gsub("%","_",PROJECT$name),".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
           par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.3), omi=c(0.2, 0.2, 0.2, 0.40))
           hist(var2, main=info, cex.main=0.9, cex=1.5, cex.axis=1.8)
           dev.off()

       } # dimension check
  } # loop through all "grid_output" objects

  # extract a list of all the variables stored in the output object
  par_names = names(grid_output)
  # filter for those related to the 'mean' status
  par_names = par_names[grepl("final", par_names)]

  # loop through these mean variables, output median and CI range for each of these variables
  for (p in seq(1,length(par_names))) {

       # determine position in grid_output list which contains the variable of interest
       pp = which(names(grid_output) == par_names[p])

       jpeg(file=paste("grid_final_map_",par_names[p],"_",gsub("%","_",PROJECT$name),".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
       par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.8), omi=c(0.2, 0.2, 0.2, 0.40))
       info = " " # assume default is no header, but sometimes we add something extra...
       var1 = mean(grid_output[[pp]][,,median_loc], na.rm=TRUE)
       var2 = mean(grid_output[[pp]][,,upper_loc], na.rm=TRUE)
       var3 = mean(grid_output[[pp]][,,lower_loc], na.rm=TRUE)
       var1 = round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2) #; var4 = round(var4,digit=2) ; var5 = round(var5,digit=2)
       info = paste("Final estimate: ",par_names[p]," (97.5 % = ",var2,"; 50 % = ",var1,"; 2.5 % = ",var3,")", sep="")
       zrange = range(pretty(c(min(grid_output[[pp]][,,median_loc], na.rm=TRUE),max(grid_output[[pp]][,,median_loc],na.rm=TRUE))))
       if (zrange[1] >= 0 & zrange[2] > 0) {
           colour_choices = colour_choices_gain
       } else if (zrange[1] < 0 & zrange[2] > 0) {
           colour_choices = colour_choices_sign
       } else if (zrange[1] < 0 & zrange[2] <= 0) {
           colour_choices = rev(colour_choices_loss)
       } else {
           colour_choices = colour_choices_default
       }
       image.plot(x = grid_long, y = grid_lat, z = grid_output[[pp]][,,median_loc], main=info, col = colour_choices, zlim=zrange,
                  axes=FALSE, cex.main=0.9, legend.width=3.0, cex=1.5, axis.args=list(cex.axis=1.8, hadj=0.1))
       map(add=TRUE, lwd = 2)
       #contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
       dev.off()

  } # loop through all "grid_output" objects

  # extract a list of all the variables stored in the output object
  par_names = names(grid_output)
  # filter for those related to the 'mean' status
  par_names = par_names[grepl("SS_", par_names)]

  # loop through these mean variables, output median and CI range for each of these variables
  for (p in seq(1,length(par_names))) {

       # determine position in grid_output list which contains the variable of interest
       pp = which(names(grid_output) == par_names[p])

       jpeg(file=paste("grid_SteadyState_map_",par_names[p],"_",gsub("%","_",PROJECT$name),".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
       par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.8), omi=c(0.2, 0.2, 0.2, 0.40))
       info = " " # assume default is no header, but sometimes we add something extra...
       var1 = mean(grid_output[[pp]][,,median_loc], na.rm=TRUE)
       var2 = mean(grid_output[[pp]][,,upper_loc], na.rm=TRUE)
       var3 = mean(grid_output[[pp]][,,lower_loc], na.rm=TRUE)
       var1 = round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2)
       info = paste("Steady State: ",par_names[p]," (97.5 % = ",var2,"; 50 % = ",var1,"; 2.5 % = ",var3,")", sep="")
       zrange = range(pretty(c(min(grid_output[[pp]][,,median_loc], na.rm=TRUE),max(grid_output[[pp]][,,median_loc],na.rm=TRUE))))
       if (zrange[1] >= 0 & zrange[2] > 0) {
           colour_choices = colour_choices_gain
       } else if (zrange[1] < 0 & zrange[2] > 0) {
           colour_choices = colour_choices_sign
       } else if (zrange[1] < 0 & zrange[2] <= 0) {
           colour_choices = rev(colour_choices_loss)
       } else {
           colour_choices = colour_choices_default
       }
       image.plot(x = grid_long, y = grid_lat, z = grid_output[[pp]][,,median_loc], main=info, col = colour_choices, zlim=zrange,
                  axes=FALSE, cex.main=0.9, legend.width=3.0, cex=1.5, axis.args=list(cex.axis=1.8, hadj=0.1))
       map(add=TRUE, lwd = 2)
       #contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
       dev.off()

  } # loop through all "grid_output" objects

  print("...done plotting stocks and fluxes maps...")

  # tidy before leaving
  gc(reset=TRUE, verbose=FALSE)

  # Move working directory
  setwd(old_wd)

} # end function generate_simplified_stock_and_fluxe_maps

## Use byte compile
generate_simplified_stock_and_flux_maps<-cmpfun(generate_simplified_stock_and_flux_maps)
