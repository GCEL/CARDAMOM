
###
## Function which allows figure generation to be spread across the cores
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).
# Exceptions are given within specific functions.

#Function to determine rmse
rmse <- function(obs, pred) sqrt(mean((obs-pred)^2, na.rm=TRUE))
## Use byte compile
rmse<-cmpfun(rmse)

uncertainty_figures<-function(n,PROJECT,load_file) {

   # This function contains all the site level state and flux plots.
   # Their selection is based on the conditional presence of particular variables
   # or specific model versions.

   # Load the current file
   load(load_file)

	 # Calculate some timing information
	 timestep = 1
	 if (PROJECT$model$timestep == "monthly") {timestep = mean(PROJECT$model$timestep_days)}
	 if (PROJECT$model$timestep == "weekly") {timestep = mean(PROJECT$model$timestep_days)}
	 time_vector = 1:dim(states_all$gpp_gCm2day)[2]
	 year_vector = time_vector/(365.25/timestep)
	 year_vector = year_vector+as.numeric(PROJECT$start_year)
	 interval = floor(length(year_vector)/10)

   # Plot rooting depth information
   plot_root_depth = FALSE
   if (PROJECT$model$name == "DALEC.A1.C3.H2.M1.#") {
       # These models assume rooting depth is controlled by coarse root, which is a fraction of the woody pool!
       var = t(states_all$roots_gCm2)
       # parameter numbers adjusted for crop model
       var = as.vector(parameters[37,,]) * (var*2) / (as.vector(parameters[36,,]) + (var*2))
       plot_root_depth = TRUE
   }  else if (PROJECT$model$name == "DALEC.A1.C2.D2.F2.H2.P4.R2.#" | PROJECT$model$name == "DALEC.A1.C2.D2.F2.H1.P4.R2.#" |
       PROJECT$model$name == "DALEC.A1.C2.D2.F2.H2.P7.R2.#" | PROJECT$model$name == "DALEC.A1.C2.D2.F2.H2.P8.R2.#" |
       PROJECT$model$name == "DALEC.A1.C2.D2.F2.H2.P3.R1.#" | PROJECT$model$name == "DALEC.A1.C2.D2.F2.H2.P10.R2.#" |
       PROJECR$model$name == "DALEC.A1.C1.D2.F2.H2.P2.#"){
       # These models assume rooting depth is controlled by coarse root, which is a fraction of the woody pool!
       tmp = t(states_all$wood_gCm2)*as.vector(parameters[29,,])
       var = t(states_all$roots_gCm2) + tmp
       # Now estimate the rooting depth based on the equation imbedded in DALEC.A1.C2.D2.F2.H2.P3.R1.
       var = as.vector(parameters[40,,]) * (var*2) / (as.vector(parameters[39,,]) + (var*2))
       plot_root_depth = TRUE
   }  else if (PROJECT$model$name == "DALEC.A1.C1.D2.F2.H2.P1.#"){
       # These models assume rooting depth is controlled by coarse root, which is a fraction of the woody pool!
       tmp = t(states_all$wood_gCm2)*as.vector(parameters[25,,])
       var = t(states_all$roots_gCm2) + tmp
       # Now estimate the rooting depth based on the equation imbedded in DALEC.A1.C2.D2.F2.H2.P3.R1.
       var = as.vector(parameters[27,,]) * (var*2) / (as.vector(parameters[26,,]) + (var*2))
       plot_root_depth = TRUE
   }  else if (PROJECT$model$name == "DALEC.A2.C1.D2.F2.H2.P1.#"){
       # These models assume rooting depth is controlled by coarse root, which is a fraction of the woody pool!
       tmp = t(states_all$wood_gCm2)*as.vector(parameters[25,,])
       var = t(states_all$roots_gCm2) + tmp
       # Now estimate the rooting depth based on the equation imbedded in DALEC.A1.C2.D2.F2.H2.P3.R1.
       var = as.vector(parameters[27,,]) * (var*2) / (as.vector(parameters[26,,]) + (var*2))
       plot_root_depth = TRUE
   }  else if (PROJECT$model$name == "DALEC.A1.C1.D2.F2.H2.P1.R1.#"){
       # These models assume rooting depth is controlled by coarse root, which is a fraction of the woody pool!
       tmp = t(states_all$wood_gCm2)*as.vector(parameters[25,,])
       var = t(states_all$roots_gCm2) + tmp
       # Now estimate the rooting depth based on the equation imbedded in DALEC.A1.C2.D2.F2.H2.P3.R1.
       var = as.vector(parameters[27,,]) * (var*2) / (as.vector(parameters[26,,]) + (var*2))
       plot_root_depth = TRUE
   }  else if (PROJECT$model$name == "DALEC.A1.C2.D2.F2.H2.P1.R1.#"){
       # These models assume rooting depth is controlled by coarse root, which is a fraction of the woody pool!
       tmp = t(states_all$wood_gCm2)*as.vector(parameters[25,,])
       var = t(states_all$roots_gCm2) + tmp
       # Now estimate the rooting depth based on the equation imbedded in DALEC.A1.C2.D2.F2.H2.P3.R1.
       var = as.vector(parameters[27,,]) * (var*2) / (as.vector(parameters[26,,]) + (var*2))
       plot_root_depth = TRUE
   } else if (PROJECT$model$name == "DALEC.A1.C2.D2.F2.H2.P2.R1.#") {
       # These models assume rooting depth is controlled by coarse root, which is a fraction of the woody pool!
       tmp = t(states_all$wood_gCm2)*as.vector(parameters[25,,])
       var = t(states_all$roots_gCm2) + tmp
       # Now estimate the rooting depth based on the equation imbedded in DALEC.A1.C2.D2.F2.H2.P3.R1.
       var = as.vector(parameters[27,,]) * (var*2) / (as.vector(parameters[26,,]) + (var*2))
       plot_root_depth = TRUE
   } else if (PROJECT$model$name == "DALEC.A1.C2.D2.F2.H2.P2.R3.#") {
       # These models assume rooting depth is controlled by coarse root, which is a fraction of the woody pool!
       tmp = t(states_all$wood_gCm2)*as.vector(parameters[25,,])
       var = t(states_all$roots_gCm2) + tmp
       # Now estimate the rooting depth based on the equation imbedded in DALEC.A1.C2.D2.F2.H2.P3.R1.
       var = as.vector(parameters[27,,]) * (var*2) / (as.vector(parameters[26,,]) + (var*2))
       plot_root_depth = TRUE
   } else if (PROJECT$model$name == "DALEC.A1.C1.D2.F2.H2.P5.#"){
       # These models assume rooting depth is controlled by coarse root, which is a fraction of the woody pool!
       tmp = t(states_all$wood_gCm2)*as.vector(parameters[25,,])
       var = t(states_all$roots_gCm2) + tmp
       # Now estimate the rooting depth based on the equation imbedded in DALEC.A1.C2.D2.F2.H2.P3.R1.
       var = as.vector(parameters[27,,]) * (var*2) / (as.vector(parameters[26,,]) + (var*2))
       plot_root_depth = TRUE
   }  else if (PROJECT$model$name == "DALEC.A1.C1.D2.F2.H2.P6.#"){
       # These models assume rooting depth is controlled by coarse root, which is a fraction of the woody pool!
       tmp = t(states_all$wood_gCm2)*as.vector(parameters[25,,])
       var = t(states_all$roots_gCm2) + tmp
       # Now estimate the rooting depth based on the equation imbedded in DALEC.A1.C2.D2.F2.H2.P3.R1.
       var = as.vector(parameters[27,,]) * (var*2) / (as.vector(parameters[26,,]) + (var*2))
       plot_root_depth = TRUE
   } # model specific plotting for root depth
   # Plot the actual time series now, assuming one has been calculated
   if (plot_root_depth) {
       jpeg(file=paste(PROJECT$figpath,"timeseries_RootDepth_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(var)[which(var != Inf)], prob=c(0.75), na.rm=TRUE)), cex=0.8,ylab="Rooting Depth (m)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")
       dev.off()
   } # plot_root_depth

   # Soil surface snow cover
   if (exists(x = "CiCa", where = states_all)) {

       # incoming data from states_all is dim=c(iter, chain, time)
       # structure needed by function is dim=c(time,iter)

       # flip it to get the right shape
       var = t(states_all$CiCa)

       jpeg(file=paste(PROJECT$figpath,"timeseries_CiCa_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
            width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(quantile(as.vector(var)[which(var != Inf)], prob=c(0.001,0.999), na.rm=TRUE)),
            cex=0.8,ylab="Ci:Ca (0-1)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),
            tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")
       dev.off()

   } # CiCa

   # Soil surface snow cover
   if (exists(x = "snow_kgH2Om2", where = states_all)) {

       # incoming data from states_all is dim=c(iter, chain, time)
       # structure needed by function is dim=c(time,iter)

       # flip it to get the right shape
       var = t(states_all$snow_kgH2Om2)

       jpeg(file=paste(PROJECT$figpath,"timeseries_snow_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
            width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(quantile(as.vector(var)[which(var != Inf)], prob=c(0.001,0.999), na.rm=TRUE)),
            cex=0.8,ylab="Soil surface snow cover (kgH2O/m2)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),
            tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")
       dev.off()

   } # snow_kgH2Om2

   # Weighted soil water potential
   if (exists(x = "wSWP_MPa", where = states_all)) {

       # incoming data from states_all is dim=c(iter, chain, time)
       # structure needed by function is dim=c(time,iter)

       # flip it to get the right shape
       var = t(states_all$wSWP_MPa)

       jpeg(file=paste(PROJECT$figpath,"timeseries_wSWP_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
            width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(quantile(as.vector(var)[which(var != Inf)], prob=c(0.001,0.999), na.rm=TRUE)),
            cex=0.8,ylab="Weighted Soil Water Potential (MPa)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),
            tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")
       dev.off()

   } # wSWP_MPa

   # Soil Surface Water (mm or kgH2O/m2) in 0-30cm depth
   if (exists(x = "SurfWater_kgH2Om2", where = states_all)) {

       # flip it to get the right shape
       var = t(states_all$SurfWater_kgH2Om2)

       jpeg(file=paste(PROJECT$figpath,"timeseries_SurfWater_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
            width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(var), prob=c(0.999), na.rm=TRUE)),
            cex = 0.8,ylab="Soil surface water (0-30 cm; kgH2O.m-2)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8,
            cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],
            labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),
            tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")
       dev.off()

   }  # SurfWater_kgH2Om2

   # Evapotranspiration (kgH2O/m2/day)
   if (exists(x = "ET_kgH2Om2day", where = states_all)) {

       # incoming data from states_all is dim=c(iter, time)
       # structure needed by function is dim=c(time,iter)
       # flip it to get the right shape
       var = t(states_all$ET_kgH2Om2day)

       obs = drivers$obs[,31] ; obs_unc = drivers$obs[,32]
       # filter -9999 to NA
       filter = which(obs == -9999) ; obs[filter] = NA ; obs_unc[filter] = NA

       jpeg(file=paste(PROJECT$figpath,"timeseries_evap_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
       if (PROJECT$model$name == "ACM" & max(as.vector(var),na.rm=TRUE) > 0) {
           par(mfrow=c(1,1), omi=c(0.1,0.1,0.1,0.1), mai=c(1,1,1,1))
           maxl=which(as.vector(parameters[dim(parameters)[1],,]) == max(as.vector(parameters[dim(parameters)[1],,])))
           maxl=maxl[1] # if we happen to have more than one values with the same likelihood we will just pick the first one....
           plot(states_all$ET_kgH2Om2day[maxl,],evap_obs, ylab="SPA", xlab="ACM", main="Evap (kgH2O.m-2.day-1)",
                pch=16,cex=0.8,cex.main=1.8,cex.lab=1.8,cex.axis=1.8)
           abline(0,1,col="red", lwd=4)
           hey=lm(obs~states_all$ET_kgH2Om2day[maxl,]) ; beta1=coef(hey)[2] ; intercept=coef(hey)[1]
           explained=summary(hey)$adj.r.squared
           trend = mean(states_all$ET_kgH2Om2day[maxl,]-obs)
           error=rmse(obs,states_all$ET_kgH2Om2day[maxl,])
           prop_error=mean(abs((obs-states_all$ET_kgH2Om2day[maxl,])/obs)*100)
           text(max(states_all$ET_kgH2Om2day[maxl,])*0.12,max(obs*0.97),
                label=bquote( SPA == .(round(beta1,3)) * ACM + .(round(intercept,3))), cex=2.)
           text(max(states_all$ET_kgH2Om2day[maxl,])*0.182,max(obs*0.92),
                label=paste("R2 = ",round(explained,3)," rmse = ",round(error,3)," bias = ",round(trend,3), sep=""), cex=2.)
           text(max(states_all$ET_kgH2Om2day[maxl,])*0.074,max(obs*0.87),
                label=paste("% error = ",round(prop_error,3), sep=""), cex=2.)
           text(max(states_all$ET_kgH2Om2day[maxl,])*0.146,max(obs*0.81),
                label=paste("ACM WUE (gC/kgH2O) = ",round(mean(states_all$gpp_gCm2day[maxl,]/states_all$ET_kgH2Om2day[maxl,]),3)," (",round(quantile(states_all$gpp_gCm2day[maxl,]/states_all$ET_kgH2Om2day[maxl,],prob=c(0.025)),3),"/",round(quantile(states_all$gpp_gCm2day[maxl,]/states_all$ET_kgH2Om2day[maxl,],prob=c(0.975)),3),")", sep=""), cex=2.)
           text(max(states_all$ET_kgH2Om2day[maxl,])*0.143,max(obs*0.76),
                label=paste("SPA WUE (gC/kgH2O) = ",round(mean(drivers$obs[,1]/obs),3)," (",round(quantile(drivers$obs[,1]/obs,prob=c(0.025)),3),"/",round(quantile(drivers$obs[,1]/obs,prob=c(0.975)),3),")", sep=""), cex=2.)
       } else {
            # now create the plotting area
            par(mfrow=c(1,1), mar=c(5,5,3,1))
            plot(obs, pch=16,xaxt="n", ylim=c(0,quantile(as.vector(var), prob=c(0.999),na.rm=TRUE)), cex=0.8,
                 ylab="Evap (kgH2O.m-2.day-1)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
                 main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
            axis(1, at=time_vector[seq(1,length(time_vector),interval)],
                 labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
            # add the confidence intervals
            plotconfidence(var)
            # calculate and draw the median values, could be mean instead or other
            lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")
            # add the data on top if there is any
            if (length(which(is.na(obs))) != length(obs) ) {
                points(obs, pch=16, cex=0.8)
                plotCI(obs,gap=0,uiw=obs_unc, col="black", add=TRUE, cex=1,lwd=2,sfrac=0.01,lty=1,pch=16)
            }
       } # acm or not
       dev.off()

   }

   # Canopy scale aerodynamic conductance (mmH2O/m2/s)
   if (exists(x = "gb_mmolH2Om2s", where = states_all)) {

       # incoming data from states_all is dim=c(iter, time)
       # structure needed by function is dim=c(time,iter)
       # flip it to get the right shape
       var = t(states_all$gs_mmolH2Om2s)
       obs = rep(NA, dim(var)[1])

       jpeg(file=paste(PROJECT$figpath,"timeseries_aerodynamic_conductance_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(obs, pch=16,xaxt="n", ylim=c(0,quantile(as.vector(var), prob=c(0.999),na.rm=TRUE)), cex=0.8,
            ylab="Canopy-scale aerodynamic conductance (mmH2O.m-2.s-1)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],
            labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")

       dev.off()

   }

   # Canopy scale stomatal conductance (mmH2O/m2/s)
   if (exists(x = "gs_mmolH2Om2s", where = states_all)) {

       # incoming data from states_all is dim=c(iter, time)
       # structure needed by function is dim=c(time,iter)
       # flip it to get the right shape
       var = t(states_all$gs_mmolH2Om2s)
       obs = rep(NA, dim(var)[1])

       jpeg(file=paste(PROJECT$figpath,"timeseries_stomatal_conductance_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(obs, pch=16,xaxt="n", ylim=c(0,quantile(as.vector(var), prob=c(0.999),na.rm=TRUE)), cex=0.8,
            ylab="Canopy-scale stomatal conductance (mmH2O.m-2.s-1)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],
            labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")

       dev.off()

   }

   # Transpiration (kgH2O/m2/day)
   if (exists(x = "Etrans_kgH2Om2day", where = states_all)) {

       # incoming data from states_all is dim=c(iter, time)
       # structure needed by function is dim=c(time,iter)
       # flip it to get the right shape
       var = t(states_all$Etrans_kgH2Om2day)
       obs = rep(NA, dim(var)[1])

       jpeg(file=paste(PROJECT$figpath,"timeseries_Etrans_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(obs, pch=16,xaxt="n", ylim=c(0,quantile(as.vector(var), prob=c(0.999),na.rm=TRUE)), cex=0.8,
            ylab="Transpiration (kgH2O.m-2.day-1)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],
            labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")

       dev.off()

   }

   # Soil evaporation (kgH2O/m2/day)
   if (exists(x = "Esoil_kgH2Om2day", where = states_all)) {

       # incoming data from states_all is dim=c(iter, time)
       # structure needed by function is dim=c(time,iter)
       # flip it to get the right shape
       var = t(states_all$Esoil_kgH2Om2day)
       obs = rep(NA, dim(var)[1])

       jpeg(file=paste(PROJECT$figpath,"timeseries_Esoil_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(obs, pch=16,xaxt="n", ylim=c(0,quantile(as.vector(var), prob=c(0.999),na.rm=TRUE)), cex=0.8,
            ylab="Soil evaporation (kgH2O.m-2.day-1)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],
            labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")

       dev.off()

   }

   # SoiWet canopy evaporation (kgH2O/m2/day)
   if (exists(x = "Ewetcanopy_kgH2Om2day", where = states_all)) {

       # incoming data from states_all is dim=c(iter, time)
       # structure needed by function is dim=c(time,iter)
       # flip it to get the right shape
       var = t(states_all$Ewetcanopy_kgH2Om2day)
       obs = rep(NA, dim(var)[1])

       jpeg(file=paste(PROJECT$figpath,"timeseries_Ewetcanopy_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(obs, pch=16,xaxt="n", ylim=c(0,quantile(as.vector(var), prob=c(0.999),na.rm=TRUE)), cex=0.8,
            ylab="Wet canopy evaporation (kgH2O.m-2.day-1)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],
            labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")

       dev.off()

   }

   # Foliage C stock (gC/m2)
   if (exists(x = "foliage_gCm2", where = states_all)) {

	     # flip it to get the right shape
		   var = t(states_all$foliage_gCm2)

		   # pass observations driver
		   obs = drivers$obs[,11] ; obs_unc = drivers$obs[,12]
		   # filter -9999 to NA
		   filter = which(obs == -9999) ; obs[filter] = NA ; obs_unc[filter] = NA

		   ymax = quantile(as.vector(var), prob=c(0.999), na.rm=TRUE)
		   ymin = quantile(as.vector(var), prob=c(0.001), na.rm=TRUE)
		   xloc = 0.15*dim(var)[1] ; yloc=(1-0.05)*ymax

		   jpeg(file=paste(PROJECT$figpath,"timeseries_foliage_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
            width=7200, height=4000, res=280, quality=100)
		   # now create the plotting area
		   par(mfrow=c(1,1), mar=c(5,5,3,1))
		   plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(ymin,ymax), cex=0.8,ylab="Foliage (gC/m2)",xlab="Time (Year)",
            cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
		   axis(1, at=time_vector[seq(1,length(time_vector),interval)],
            labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
		   # add the confidence intervals
		   plotconfidence(var)
   	   # calculate and draw the median values, could be mean instead or other
     	 lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")
		   # add the data on top if there is any
   		 if (length(which(is.na(obs))) != length(obs) ) {
		     	 points(obs, pch=16, cex=0.8)
			     plotCI(obs,gap=0,uiw=Cfol_obs_unc, col="black", add=TRUE, cex=1,lwd=2,sfrac=0.01,lty=1,pch=16)
		   }
		   dev.off()

	  } # foliage_gCm2

    # leaf area index (m2m2)
    if (exists(x = "lai_m2m2", where = states_all)) {

        # flip it to get the right shape
		    var = t(states_all$lai_m2m2)
		    obs = drivers$obs[,3]

  		  jpeg(file=paste(PROJECT$figpath,"timeseries_lai_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
             width=7200, height=4000, res=280, quality=100)
		    # now create the plotting area
		    par(mfrow=c(1,1), mar=c(5,5,3,1))
	   	  plot(obs, pch=16,xaxt="n", ylim=c(0,max(max(obs),quantile(as.vector(var), prob=c(0.999), na.rm=TRUE))),
             cex=0.8,ylab="LAI (m2/m2)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
             main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
		    axis(1, at=time_vector[seq(1,length(time_vector),interval)],
             labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
		    # add the confidence intervals
		    plotconfidence(var)
		    # calculate and draw the median values, could be mean instead or other
		    lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")
		    # add the data on top
		    points(obs, pch=16, cex=0.8)
    		dev.off()

   } # lai_m2m2

   # GPP (gC/m2/day)
   if (exists(x = "gpp_gCm2day", where = states_all)) {

       # flip it to get the right shape
       var = t(states_all$gpp_gCm2day)
       obs = drivers$obs[,1] ; obs_unc = drivers$obs[,2]
       # filter -9999 to NA
       filter = which(obs == -9999) ; obs[filter] = NA ; obs_unc[filter] = NA

       jpeg(file=paste(PROJECT$figpath,"timeseries_gpp_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
       if (PROJECT$model$name == "ACM") {
           par(mfrow=c(1,1), omi=c(0.1,0.1,0.1,0.1), mai=c(1,1,1,1))
           maxl = which(as.vector(parameters[dim(parameters)[1],,]) == max(as.vector(parameters[dim(parameters)[1],,])))
           maxl = maxl[1] # if we happen to have more than one values with the same likelihood we will just pick the first one....
           plot(states_all$gpp_gCm2day[maxl,],drivers$obs[,1], ylab="SPA", xlab="ACM", main="GPP (gC/m2/day)",
                pch=16,cex=0.8,cex.main=1.8,cex.lab=1.8,cex.axis=1.8) ; abline(0,1,col="red", lwd=4)
           hey = lm(drivers$obs[,1]~states_all$gpp_gCm2day[maxl,]) ; beta1=coef(hey)[2] ; intercept=coef(hey)[1]
           explained = summary(hey)$adj.r.squared ; trend = mean(states_all$gpp_gCm2day[maxl,]-drivers$obs[,1]) ; error=rmse(drivers$obs[,1],states_all$gpp_gCm2day[maxl,])
           prop_error = mean(abs((drivers$obs[,1]-states_all$gpp_gCm2day[maxl,])/drivers$obs[,1])*100)
           text(max(states_all$gpp_gCm2day[maxl,])*0.12,max(drivers$obs[,1]*0.97),
                label=bquote( SPA == .(round(beta1,3)) * ACM + .(round(intercept,3))), cex=2.)
           text(max(states_all$gpp_gCm2day[maxl,])*0.182,max(drivers$obs[,1]*0.92),
                label=paste("R2 = ",round(explained,3)," rmse = ",round(error,3)," bias = ",round(trend,3), sep=""), cex=2.)
           text(max(states_all$gpp_gCm2day[maxl,])*0.074,max(drivers$obs[,1]*0.87),
                label=paste("% error = ",round(prop_error,3), sep=""), cex=2.)
           text(max(states_all$gpp_gCm2day[maxl,])*0.146,max(drivers$obs[,1]*0.81),
                label=paste("ACM NUE = ",round(mean(states_all$gpp_gCm2day[maxl,]/drivers$met[,14]),3)," (",round(quantile(states_all$gpp_gCm2day[maxl,]/drivers$met[,14],prob=c(0.025)),3),"/",round(quantile(states_all$gpp_gCm2day[maxl,]/drivers$met[,14],prob=c(0.975)),3),")", sep=""), cex=2.)
           text(max(states_all$gpp_gCm2day[maxl,])*0.143,max(drivers$obs[,1]*0.76),
                label=paste("SPA NUE = ",round(mean(obs/drivers$met[,14]),3)," (",round(quantile(obs/drivers$met[,14],prob=c(0.025)),3),"/",round(quantile(gpp_obs/drivers$met[,14],prob=c(0.975)),3),")", sep=""), cex=2.)
       } else {
           # create the plotting area
           par(mfrow=c(1,1), mar=c(5,5,3,1))
           plot(obs, pch=16,xaxt="n", ylim=c(0,quantile(as.vector(var), prob=c(0.999),na.rm=TRUE)),
                cex=0.8,ylab="GPP (gC/m2/day)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
                main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
           axis(1, at=time_vector[seq(1,length(time_vector),interval)],
                labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
           # add the confidence intervals
           plotconfidence(var)
           # calculate and draw the median values, could be mean instead or other
           lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")
           # add the data on top if there is any
           if (length(which(is.na(obs))) != length(obs) ) {
               points(obs, pch=16, cex=0.8)
               plotCI(obs, gap=0, uiw=obs_unc, col="black", add=TRUE, cex=1,lwd=2,sfrac=0.01,lty=1,pch=16)
           }
       } # acm or not
       dev.off()

   } # gpp_gCm2day

   # Net Ecosystem Exchange of CO2 (gC/m2day)
   if (exists(x = "nee_gCm2day", where = states_all)) {

   	   # flip it to get the right shape
       var = t(states_all$nee_gCm2day)
       # pass observations driver
       obs = drivers$obs[,5] ; obs_unc = drivers$obs[,6]
       # filter -9999 to NA
       filter = which(obs == -9999) ; obs[filter] = NA ; obs_unc[filter] = NA

       jpeg(file=paste(PROJECT$figpath,"timeseries_nee_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
            width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(obs, pch=16,xaxt="n", ylim=c(quantile(as.vector(var),prob=c(0.001),na.rm=TRUE),quantile(as.vector(var),prob=c(0.999),na.rm=TRUE)),
            cex=0.8,ylab="NEE (gC/m2/day)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],
            labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")
       # add the data on top if there is any
       if (length(which(is.na(obs))) != length(obs) ) {
           points(obs, pch=16, cex=0.8)
           if (length(which(is.na(obs)))/length(obs) > 0.01) {
               plotCI(obs,gap=0,uiw=obs_unc, col="black", add=TRUE, cex=1,lwd=2,sfrac=0.01,lty=1,pch=16)
           }
       }

       dev.off()

   } # nee_gCm2day

   # Ecosystem respiration (gC/m2/day)
   if (exists(x = "reco_gCm2day", where = states_all)) {

       # flip it to get the right shape
       var = t(states_all$reco_gCm2day)
       # pass observations driver
       obs=drivers$obs[,9] ; obs_unc=drivers$obs[,10]
       # filter -9999 to NA
       filter = which(obs == -9999) ; obs[filter] = NA ; obs_unc[filter] = NA

       jpeg(file=paste(PROJECT$figpath,"timeseries_eco_resp_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
            width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(var), prob=c(0.999),na.rm=TRUE)),
            cex=0.8,ylab="Ecosystem respiration (gC/m2/day)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],
            labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")
       # add the data on top if there is any
       if (length(which(is.na(obs))) != length(obs) ) {
           points(obs, pch=16, cex=0.8)
           plotCI(obs,gap=0,uiw=obs_unc, col="black", add=TRUE, cex=1,lwd=2,sfrac=0.01,lty=1,pch=16)
       }
       dev.off()

   } # reco_gCm2day

   # Total Heterotrophic respiration
   if (exists(x = "rhet_gCm2day", where = states_all)) {

       # flip it to get the right shape
       var = t(states_all$rhet_gCm2day)

       jpeg(file=paste(PROJECT$figpath,"timeseries_heterotrophic_respiration_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
            width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(var), prob=c(0.999),na.rm=TRUE)),
            cex=0.8,ylab="Heterotrophic respiration (gC/m2/day)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8,
            cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)],
            digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")

       dev.off()

   } # rhet_gCm2day

   # labile (gC/m2)
   if (exists(x = "labile_gCm2", where = states_all)) {

       # flip it to get the right shape
       var=t(states_all$labile_gCm2)

       jpeg(file=paste(PROJECT$figpath,"timeseries_labile_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
            width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(var), prob=c(0.999), na.rm=TRUE)),
            cex=0.8,ylab="labile (gC/m2)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],
            labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")
       dev.off()

   } # labile_gCm2

   # Foliage + fine root litter (gCm2)
   if (exists(x = "litter_gCm2", where = states_all)) {

       # flip it to get the right shape
       var = t(states_all$litter_gCm2)

       # pass observations driver
       obs = drivers$obs[,17] ; obs_unc = drivers$obs[,18]
       # filter -9999 to NA
       filter = which(obs == -9999) ; obs[filter] = NA ; obs_unc[filter] = NA

       jpeg(file=paste(PROJECT$figpath,"timeseries_litter_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
            width=7200, height=4000, res=280, quality=100)

       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(var[1:(dim(var)[1]-1),]), prob=c(0.999), na.rm=TRUE)),
            cex=0.8,ylab="litter (gC/m2)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],
            labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")
       # add the data on top if there is any
       if (length(which(is.na(obs))) != length(obs) ) {
           points(obs, pch=16, cex=0.8)
           plotCI(obs,gap=0,uiw=obs_unc, col="black", add=TRUE, cex=1,lwd=2,sfrac=0.01,lty=1,pch=16)
       }
       dev.off()

   } # litter_gCm2

   # Foliage + fine root litter (gCm2)
   if (exists(x = "woodlitter_gCm2", where = states_all)) {

       # flip it to get the right shape
       var = t(states_all$woodlitter_gCm2)

#       # pass observations driver
#       obs = drivers$obs[,17] ; obs_unc = drivers$obs[,18]
#       # filter -9999 to NA
#       filter = which(obs == -9999) ; obs[filter] = NA ; obs_unc[filter] = NA

       jpeg(file=paste(PROJECT$figpath,"timeseries_woodlitter_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
            width=7200, height=4000, res=280, quality=100)

       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(var[1:(dim(var)[1]-1),]), prob=c(0.999), na.rm=TRUE)),
            cex=0.8,ylab="wood litter (gC/m2)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],
            labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")
       # add the data on top if there is any
#       if (length(which(is.na(obs))) != length(obs) ) {
#           points(obs, pch=16, cex=0.8)
#           plotCI(obs,gap=0,uiw=obs_unc, col="black", add=TRUE, cex=1,lwd=2,sfrac=0.01,lty=1,pch=16)
#       }
       dev.off()

   } # litter_gCm2

   # Fine roots (gC/m2)
   if (exists(x = "roots_gCm2", where = states_all)) {

       # flip it to get the right shape
       var=t(states_all$roots_gCm2)

       # pass observations driver
       obs=drivers$obs[,15] ; obs_unc=drivers$obs[,16]
       # filter -9999 to NA
       filter = which(obs == -9999) ; obs[filter] = NA ; obs_unc[filter] = NA

       jpeg(file=paste(PROJECT$figpath,"timeseries_roots_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
            width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(var), prob=c(0.999), na.rm=TRUE)),
            cex=0.8,ylab="Roots (gC/m2)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],
            labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")
       # add the data on top if there is any
       if (length(which(is.na(obs))) != length(obs) ) {
           points(obs, pch=16, cex=0.8)
           plotCI(obs,gap=0,uiw=obs_unc, col="black", add=TRUE, cex=1,lwd=2,sfrac=0.01,lty=1,pch=16)
       }
       dev.off()

   } # roots_gCm2


   # Wood stocks (gC/m2)
   if (exists(x = "wood_gCm2", where = states_all)) {

       # flip it to get the right shape
       var=t(states_all$wood_gCm2)
       # pass observations driver
       obs=drivers$obs[,13] ; obs_unc=drivers$obs[,14]
       # filter -9999 to NA
       filter = which(obs == -9999) ; obs[filter] = NA ; obs_unc[filter] = NA

       jpeg(file=paste(PROJECT$figpath,"timeseries_wood_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
            width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(var), prob=c(0.999), na.rm=TRUE)),
            cex=0.8,ylab="Wood (gC/m2)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")
       # add the data on top if there is any
       if (length(which(is.na(obs))) != length(obs) ) {
           points(obs, pch=16, cex=0.8)
           plotCI(obs,gap=0,uiw=obs_unc, col="black", add=TRUE, cex=1,lwd=2,sfrac=0.01,lty=1,pch=16)
       }
       dev.off()

   } # wood_gCm2

   # soil organic carbon (gC/m2)
   if (exists(x = "som_gCm2", where = states_all)) {

       # flip it to get the right shape
       var = t(states_all$som_gCm2)

       # pass observations driver
       obs = drivers$obs[,19] ; obs_unc = drivers$obs[,20]
       # filter -9999 to NA
       filter = which(obs == -9999) ; obs[filter] = NA ; obs_unc[filter] = NA

       jpeg(file=paste(PROJECT$figpath,"timeseries_som_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
            width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(var), prob=c(0.999), na.rm=TRUE)),
            cex=0.8,ylab="som (gC/m2)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],
            labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")
       # add the data if there is any which is not missing
       if (length(which(is.na(obs))) != length(obs) ) {
           # add the data on top
           points(obs, pch=16, cex=0.8)
           plotCI(obs,gap=0,uiw=obs_unc, col="black", add=TRUE, cex=1,lwd=2,sfrac=0.01,lty=1,pch=16)
       }
       dev.off()

    } # som_gCm2

    # Biomass (gC/m2)
    if (exists(x = "biomass_gCm2", where = states_all)) {

        # structure needed by function is dim=c(time,iter)
        # flip it to get the right shape
        var=t(states_all$biomass_gCm2)

        jpeg(file=paste(PROJECT$figpath,"timeseries_Cbiomass_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
             width=7200, height=4000, res=280, quality=100)
        # now create the plotting area
        par(mfrow=c(1,1), mar=c(5,5,3,1))
        plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(var), prob=c(0.999), na.rm=TRUE)),
             cex=0.8,ylab="Biomass (gC/m2)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
             main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
        axis(1, at=time_vector[seq(1,length(time_vector),interval)],
             labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
        # add the confidence intervals
        plotconfidence(var)
        # calculate and draw the median values, could be mean instead or other
        lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")

        dev.off()

   }  # biomass_gCm2

   # Canopy growth index (CGI; 0-1)
   if (exists(x = "ncce_gCm2day", where = states_all)) {

       # structure needed by function is dim=c(time,iter)
       # flip it to get the right shape
       var=t(states_all$ncce_gCm2day)

       jpeg(file=paste(PROJECT$figpath,"timeseries_NCCE_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
            width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(var), prob=c(0.999), na.rm=TRUE)),
            cex=0.8,ylab="NCCE (gC/m2/day)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],
            labels=round(year_vector[seq(1,length(time_vector),interval)], digits=1),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")
       dev.off()

   } # ncce_gCm2day

   # Canopy growth index (CGI; 0-1)
   if (exists(x = "cgi", where = states_all)) {

       # structure needed by function is dim=c(time,iter)
       # flip it to get the right shape
       var=t(states_all$cgi)

       jpeg(file=paste(PROJECT$figpath,"timeseries_CGI_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
            width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(var), prob=c(0.999), na.rm=TRUE)),
            cex=0.8,ylab="CGI",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],
            labels=round(year_vector[seq(1,length(time_vector),interval)], digits=1),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")
       dev.off()

   } # cgi

   # Canopy mortality index (CMI; 0-1)
   if (exists(x = "cmi", where = states_all)) {

       # structure needed by function is dim=c(time,iter)
       # flip it to get the right shape
       var=t(states_all$cmi)

       jpeg(file=paste(PROJECT$figpath,"timeseries_CMI_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
            width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(var), prob=c(0.999), na.rm=TRUE)),
            cex=0.8,ylab="CMI",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],
            labels=round(year_vector[seq(1,length(time_vector),interval)], digits=1),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")
       dev.off()

   } # cmi

   # Growing season index (GSI; 0-1)
   if (exists(x = "gsi", where = states_all)) {

       # structure needed by function is dim=c(time,iter)
       # flip it to get the right shape
       var=t(states_all$gsi)

       jpeg(file=paste(PROJECT$figpath,"timeseries_GSI_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
            width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(var), prob=c(0.999), na.rm=TRUE)),
            cex=0.8,ylab="GSI",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],
            labels=round(year_vector[seq(1,length(time_vector),interval)], digits=1),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")
       dev.off()

   } # gsi

   # Growing season index - photoperiod component (GSI; 0-1)
   if (exists(x = "gsi_iphoto", where = states_all)) {

       # structure needed by function is dim=c(time,iter)
       # flip it to get the right shape
       var=t(states_all$gsi_iphoto)

       jpeg(file=paste(PROJECT$figpath,"timeseries_GSI_photoperiod_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
            width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(var), prob=c(0.999), na.rm=TRUE)),
            cex=0.8,ylab="GSI-iphoto",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],
            labels=round(year_vector[seq(1,length(time_vector),interval)], digits=1),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")
       dev.off()

   } # gsi_iphoto

   # Growing season index - maximum temperature component (GSI; 0-1)
   if (exists(x = "gsi_itemp", where = states_all)) {

       # structure needed by function is dim=c(time,iter)
       # flip it to get the right shape
       var=t(states_all$gsi_itemp)

       jpeg(file=paste(PROJECT$figpath,"timeseries_GSI_temperature_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
            width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(var), prob=c(0.999), na.rm=TRUE)),
            cex=0.8,ylab="GSI-itemp",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],
            labels=round(year_vector[seq(1,length(time_vector),interval)], digits=1),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")
       dev.off()

   } # gsi_itemp

   # Growing season index - vapour pressure component (GSI; 0-1)
   if (exists(x = "gsi_ivpd", where = states_all)) {

       # structure needed by function is dim=c(time,iter)
       # flip it to get the right shape
       var=t(states_all$gsi_ivpd)

       jpeg(file=paste(PROJECT$figpath,"timeseries_GSI_vpd_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
            width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(var), prob=c(0.999), na.rm=TRUE)),
            cex=0.8,ylab="GSI-ivpd",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],
            labels=round(year_vector[seq(1,length(time_vector),interval)], digits=1),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")
       dev.off()

   } # gsi_ivpd

   # Growing season index - wSWP component (GSI; 0-1)
   if (exists(x = "gsi_iwswp", where = states_all)) {

       # structure needed by function is dim=c(time,iter)
       # flip it to get the right shape
       var=t(states_all$gsi_iwswp)

       jpeg(file=paste(PROJECT$figpath,"timeseries_GSI_wSWP_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
            width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(var), prob=c(0.999), na.rm=TRUE)),
            cex=0.8,ylab="GSI-iwSWP",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],
            labels=round(year_vector[seq(1,length(time_vector),interval)], digits=1),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")
       dev.off()

   } # gsi_iwSWP

   # Total ecosystem harvested C (gC/m2/day)
   if (exists(x = "harvest_gCm2day", where = states_all)) {

       # structure needed by function is dim=c(time,iter)
       # flip it to get the right shape
       var=t(states_all$harvest_gCm2day)

       ymax=quantile(as.vector(var), prob=c(0.999), na.rm=TRUE)
       jpeg(file=paste(PROJECT$figpath,"timeseries_harvestedC_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
            width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(0,ymax), cex=0.8,ylab="Harvested C",xlab="Time (Year)",
       cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],
            labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")

       dev.off()

   } # harvest

   # Total ecosystem harvested C (gC/m2/day)
   if (exists(x = "grid_output", where = states_all)) {

       # structure needed by function is dim=c(time,iter)
       # flip it to get the right shape
       var=t(states_all$canopyage_days)

       ymax=quantile(as.vector(var), prob=c(0.999), na.rm=TRUE)
       jpeg(file=paste(PROJECT$figpath,"timeseries_CanopyAge_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
            width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(0,ymax), cex=0.8,
            ylab="Mean Canopy Age (days)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")

       dev.off()

   } # canopy age

   # Ecosystem fire C emissions (gC/m2/day)
   if (exists(x = "fire_gCm2day", where = states_all)) {

       # structure needed by function is dim=c(time,iter)
       # flip it to get the right shape
       var = t(states_all$fire_gCm2day)

       ymax = quantile(as.vector(var), prob=c(0.999), na.rm=TRUE)
       jpeg(file=paste(PROJECT$figpath,"timeseries_fire_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""),
            width=7200, height=4000, res=280, quality=100)
       # now create the plotting area
       par(mfrow=c(1,1), mar=c(5,5,3,1))
       plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(0,ymax), cex=0.8,
            ylab="Fire (gC/m2/day)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8,
            main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
       axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
       # add the confidence intervals
       plotconfidence(var)
       # calculate and draw the median values, could be mean instead or other
       lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="red")

       dev.off()

   } # fire_gCm2day

   # tidy before leaving
   gc(reset=TRUE, verbose=FALSE)

} # end of function
## Use byte compile
uncertainty_figures<-cmpfun(uncertainty_figures)
