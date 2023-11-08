
###
## Function to create generic stock and flux plots for gridded CARDAMOM analyses
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

generate_simplified_stock_and_flux_maps<-function(PROJECT) {

  print("...beginning stocks and fluxes...")

  # Move working directory
  old_wd = getwd()
  setwd(PROJECT$figpath)

  # read in the grid stock information
  infile = paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep="")
  if (file.exists(infile) == FALSE) {stop("grid_outputs for 'generate_simplified_stock_and_flux_maps' missing")}
  load(paste(infile))

  # work out area matrix for the pixels in meters
  # include adjustment for g-> Tg (*1e-12)
  if (PROJECT$grid_type == "UK") {
      output = generate_uk_grid(PROJECT$latitude,PROJECT$longitude,PROJECT$resolution)
      grid_lat = array(output$lat, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
      grid_long = array(output$long,dim=c(PROJECT$long_dim,PROJECT$lat_dim))
      area_with_g_Tg = array(PROJECT$resolution**2, dim=c(PROJECT$long_dim,PROJECT$lat_dim))*1e-12
      area = array(PROJECT$resolution**2, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
      rm(output)
  } else if (PROJECT$grid_type == "wgs84") {
      # generate the lat / long grid again
      output = generate_wgs84_grid(PROJECT$latitude,PROJECT$longitude,PROJECT$resolution)
      grid_lat = array(output$lat, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
      grid_long = array(output$long,dim=c(PROJECT$long_dim,PROJECT$lat_dim))
      # then generate the area estimates for each pixel
      area_with_g_Tg = calc_pixel_area(grid_long,grid_lat)*1e-12
      area = calc_pixel_area(grid_long,grid_lat)
      # this output is in vector form and we need matching array shapes so...
      area = array(area, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
      area_with_g_Tg = array(area_with_g_Tg, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
      rm(output)
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
  fig_height = 3000*0.65 ; fig_width = ((PROJECT$long_dim/PROJECT$lat_dim)+0.25) * fig_height

  # If root depth information has been provided plot some of it up here.
  if (exists(x = "mean_RootDepth_m", where = grid_output)) {
      jpeg(file=paste("median_root_depth_maps_",gsub("%","_",PROJECT$name),".jpeg",sep=""), width=fig_width, height=fig_height, res=300, quality=100)
      par(mfrow=c(1,1), mar=c(1.2, 1.0, 2.2, 6.6), omi=c(0.2, 0.2, 0.2, 0.40))
      z_axis=c(min(as.vector(grid_output$mean_RootDepth_m),na.rm=TRUE),max(as.vector(grid_output$mean_RootDepth_m),na.rm=TRUE))
      image.plot(x = grid_long, y = grid_lat, z = grid_output$mean_RootDepth_m[,,median_loc],col=colour_choices
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
           #var4 = mean(grid_output[[pp]][,,loc_25], na.rm=TRUE)
           #var5 = mean(grid_output[[pp]][,,loc_75], na.rm=TRUE)
           var1 = round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2) #; var4 = round(var4,digit=2) ; var5 = round(var5,digit=2)
           #var1 = round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2) ; var4 = round(var4,digit=2) ; var5 = round(var5,digit=2)
           #info = paste("Mean estimate: ",par_names[p]," (97.5 % = ",var2,"; 75 % = ",var5,"; 50 % = ",var1,"; 25 % = ",var4,"; 2.5 % = ",var3,")", sep="")
           info = paste("Mean estimate: ",par_names[p]," (97.5 % = ",var2,"; 50 % = ",var1,"; 2.5 % = ",var3,")", sep="")
           zrange=range(pretty(c(min(grid_output[[pp]][,,median_loc], na.rm=TRUE),max(grid_output[[pp]][,,median_loc],na.rm=TRUE))))
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
           #var2 = mean(grid_output[[pp]][,,upper_loc], na.rm=TRUE)
           #var3 = mean(grid_output[[pp]][,,lower_loc], na.rm=TRUE)
           #var4 = mean(grid_output[[pp]][,,loc_25], na.rm=TRUE)
           #var5 = mean(grid_output[[pp]][,,loc_75], na.rm=TRUE)
           var1 = round(var1,digit=2) #; var2=round(var2,digit=2) ; var3=round(var3,digit=2) #; var4 = round(var4,digit=2) ; var5 = round(var5,digit=2)
           #var1 = round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2) ; var4 = round(var4,digit=2) ; var5 = round(var5,digit=2)
           #info = paste("Mean estimate: ",par_names[p]," (97.5 % = ",var2,"; 75 % = ",var5,"; 50 % = ",var1,"; 25 % = ",var4,"; 2.5 % = ",var3,")", sep="")
           info = paste("Mean 95CI estimate: ",par_names[p]," = ",var1, sep="")
           zrange=range(pretty(c(min(var2, na.rm=TRUE),max(var2,na.rm=TRUE))))
           image.plot(x = grid_long, y = grid_lat, z = var2, zlim=zrange, main=info, col = colour_choices,
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
       #var4 = mean(grid_output[[pp]][,,loc_25], na.rm=TRUE)
       #var5 = mean(grid_output[[pp]][,,loc_75], na.rm=TRUE)
       var1 = round(var1,digit=2) ; var2=round(var2,digit=2) ; var3=round(var3,digit=2) #; var4 = round(var4,digit=2) ; var5 = round(var5,digit=2)
#       info = paste("Final estimate: ",par_names[p]," (97.5 % = ",var2,"; 75 % = ",var5,"; 50 % = ",var1,"; 25 % = ",var4,"; 2.5 % = ",var3,")", sep="")
       info = paste("Final estimate: ",par_names[p]," (97.5 % = ",var2,"; 50 % = ",var1,"; 2.5 % = ",var3,")", sep="")
       zrange=range(pretty(c(min(grid_output[[pp]][,,median_loc], na.rm=TRUE),max(grid_output[[pp]][,,median_loc],na.rm=TRUE))))
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
       zrange=range(pretty(c(min(grid_output[[pp]][,,median_loc], na.rm=TRUE),max(grid_output[[pp]][,,median_loc],na.rm=TRUE))))
       image.plot(x = grid_long, y = grid_lat, z = grid_output[[pp]][,,median_loc], main=info, col = colour_choices, zlim=zrange,
                  axes=FALSE, cex.main=0.9, legend.width=3.0, cex=1.5, axis.args=list(cex.axis=1.8, hadj=0.1))
       map(add=TRUE, lwd = 2)
       #contour(landmask, add = TRUE, lwd=1.0, nlevels=1,axes=FALSE,drawlabels=FALSE,col="black")
       dev.off()

  } # loop through all "grid_output" objects

  # tidy before leaving
  gc(reset=TRUE, verbose=FALSE)

  # Move working directory
  setwd(old_wd)

} # end function generate_simplified_stock_and_fluxe_maps

## Use byte compile
generate_simplified_stock_and_flux_maps<-cmpfun(generate_simplified_stock_and_flux_maps)
