
###
## Function to create generic stock and flux plots for gridded CARDAMOM analyses
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

generate_simplified_stock_and_flux_maps<-function(PROJECT) {

  print("...beginning stocks and fluxes...")

  # read in the grid stock information
  infile = paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep="")
  if (file.exists(infile) == FALSE) {stop("grid_outputs for 'generate_simplified_stock_and_flux_maps' missing")}
  load(paste(infile))
  # read in the grid parameter information
  infile=paste(PROJECT$results_processedpath,PROJECT$name,"_parameter_maps.RData",sep="")
  if (file.exists(infile) == FALSE & grepl("BUCKET",PROJECT$model$name)) {stop("parameter_maps for 'generate_simplified_stock_and_flux_maps' missing")}
  load(paste(infile))

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

  # determine the array value for the median,
  num_quantiles = dim(grid_output$mean_labile_gCm2)[3]
  if (num_quantiles == 7) {
      # then we assume we are dealing with 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975 quantiles
      median_loc = 4 ; loc_25 = 3 ; loc_75 = 5; lower_loc = 1 ; upper_loc = 7
  } else {
      # otherwise we need to approximate it...
      median_loc = round(num_quantiles / 2,digits=0)
      lower_loc = ceiling(num_quantiles * 0.025)
      upper_loc = floor(num_quantiles * 0.975)
      loc_25 = ceiling(num_quantiles * 0.25)
      loc_75 = floor(num_quantiles * 0.75)
  }

  # calculate land mask
  landmask = array(PROJECT$landsea, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
  # load colour palette
  colour_choices_upper = colorRampPalette(brewer.pal(11,"Spectral"))
  # array sizes are always the same so
  colour_choices = colour_choices_upper(length(area))

  # determine correct height and widths
  fig_height = 7000 ; fig_width = ((PROJECT$long_dim/PROJECT$lat_dim)+0.25) * fig_height

  # create a map summarising the rooting depth information
  mean_rooting_depth = NA
  if (PROJECT$model$name ==  "DALEC" | PROJECT$model$name == "DALEC_BUCKET" |
      PROJECT$model$name == "DALEC_G5" | PROJECT$model$name == "DALEC_G6" |
      PROJECT$model$name == "DALEC_GSI_BUCKET" | PROJECT$model$name == "DALEC_BUCKET_CanAGE") {

      jpeg(file=paste(PROJECT$figpath,"median_root_depth_maps_",PROJECT$name,".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.6), omi=c(0.2, 0.2, 0.2, 0.40))
      root_biomass = (grid_output$mean_roots_gCm2[,,median_loc]+(grid_output$mean_wood_gCm2[,,median_loc]*grid_parameters$parameters[,,29,median_loc]))*2
      mean_rooting_depth = grid_parameters$parameters[,,40,median_loc] * root_biomass / (grid_parameters$parameters[,,39,median_loc] + root_biomass)
      z_axis=c(min(as.vector(mean_rooting_depth),na.rm=TRUE),max(as.vector(mean_rooting_depth),na.rm=TRUE))
      image.plot(mean_rooting_depth,col=colour_choices, main=paste("Median root depth (m)",sep=""),zlim=z_axis,axes=FALSE
                ,cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
      contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
      dev.off()

      jpeg(file=paste(PROJECT$figpath,"median_max_root_depth_maps_",PROJECT$name,".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.6), omi=c(0.2, 0.2, 0.2, 0.40))
      root_biomass = (grid_output$annual_max_roots_gCm2[,,median_loc]+(grid_output$annual_max_wood_gCm2[,,median_loc]*grid_parameters$parameters[,,29,median_loc]))*2
      mean_rooting_depth = grid_parameters$parameters[,,40,median_loc] * root_biomass / (grid_parameters$parameters[,,39,median_loc] + root_biomass)
      z_axis=c(min(as.vector(mean_rooting_depth),na.rm=TRUE),max(as.vector(mean_rooting_depth),na.rm=TRUE))
      image.plot(mean_rooting_depth,col=colour_choices, main=paste("Annual max root depth (m)",sep=""),zlim=z_axis,axes=FALSE
                ,cex.main=2.4,legend.width=3.0,cex=1.5,axis.args=list(cex.axis=1.8,hadj=0.1))
      contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
      dev.off()

  } # rooting depth

  # extract a list of all the variables stored in the output object
  par_names = names(grid_output)
  # filter for those related to the 'mean' status
  par_names = par_names[grepl("mean", par_names)]

  # loop through these mean variables, output median and CI range for each of these variables
  for (p in seq(1,length(par_names))) {

       # determine position in grid_output list which contains the variable of interest
       pp = which(names(grid_output) == par_names[p])

       # create maps
       jpeg(file=paste(PROJECT$figpath,"grid_mean_map_",par_names[p],"_",PROJECT$name,".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
       par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.8), omi=c(0.2, 0.2, 0.2, 0.40))
       info = " " # assume default is no header, but sometimes we add something extra...
       var1 = mean(grid_output[[pp]][,,median_loc], na.rm=TRUE)
       var2 = mean(grid_output[[pp]][,,upper_loc], na.rm=TRUE)
       var3 = mean(grid_output[[pp]][,,lower_loc], na.rm=TRUE)
       var4 = mean(grid_output[[pp]][,,loc_25], na.rm=TRUE)
       var5 = mean(grid_output[[pp]][,,loc_75], na.rm=TRUE)
       var1 = round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2) ; var4 = round(var4,digit=2) ; var5 = round(var5,digit=2)
       info = paste("Mean estimate: ",par_names[p]," (97.5 % = ",var2,"; 75 % = ",var5,"; 50 % = ",var1,"; 25 % = ",var4,"; 2.5 % = ",var3,")", sep="")
       image.plot(grid_output[[pp]][,,median_loc], main=info, col = colour_choices, axes=FALSE, cex.main=2.4, legend.width=3.0, cex=1.5, axis.args=list(cex.axis=1.8, hadj=0.1))
       contour(landmask, add = TRUE, lwd = 1.0, nlevels = 1,axes = FALSE,drawlabels = FALSE,col = "black")
       dev.off()

       # create distribution information
       jpeg(file=paste(PROJECT$figpath,"grid_mean_hist_",par_names[p],"_",PROJECT$name,".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
       par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.8), omi=c(0.2, 0.2, 0.2, 0.40))
       info = paste("Distribution of median estimates: ",par_names[p], sep="")
       hist(as.vector(grid_output[[pp]][,,median_loc]), main=info, col = "grey", cex.main=2.4, cex=1.5)
       dev.off()

  } # loop through all "grid_output" objects

  # extract a list of all the variables stored in the output object
  par_names = names(grid_output)
  # filter for those related to the 'mean' status
  par_names = par_names[grepl("final", par_names)]

  # loop through these mean variables, output median and CI range for each of these variables
  for (p in seq(1,length(par_names))) {

       # determine position in grid_output list which contains the variable of interest
       pp = which(names(grid_output) == par_names[p])

       jpeg(file=paste(PROJECT$figpath,"grid_final_map_",par_names[p],"_",PROJECT$name,".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
       par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.8), omi=c(0.2, 0.2, 0.2, 0.40))
       info = " " # assume default is no header, but sometimes we add something extra...
       var1 = mean(grid_output[[pp]][,,median_loc], na.rm=TRUE)
       var2 = mean(grid_output[[pp]][,,upper_loc], na.rm=TRUE)
       var3 = mean(grid_output[[pp]][,,lower_loc], na.rm=TRUE)
       var4 = mean(grid_output[[pp]][,,loc_25], na.rm=TRUE)
       var5 = mean(grid_output[[pp]][,,loc_75], na.rm=TRUE)
       var1 = round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2) ; var4 = round(var4,digit=2) ; var5 = round(var5,digit=2)
       info = paste("Final estimate: ",par_names[p]," (97.5 % = ",var2,"; 75 % = ",var5,"; 50 % = ",var1,"; 25 % = ",var4,"; 2.5 % = ",var3,")", sep="")
       image.plot(grid_output[[pp]][,,median_loc], main=info, col = colour_choices, axes=FALSE, cex.main=2.4, legend.width=3.0, cex=1.5, axis.args=list(cex.axis=1.8, hadj=0.1))
       contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
       dev.off()

  } # loop through all "grid_output" objects

  # tidy before leaving
  gc(reset=TRUE, verbose=FALSE)

} # end function generate_simplified_stock_and_fluxe_maps
