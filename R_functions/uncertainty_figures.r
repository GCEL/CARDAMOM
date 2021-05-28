
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

uncertainty_figures<-function(which_plot,PROJECT,states_all,drivers,parameters,n,plotconfidence) {

	# calculate some timing information
	timestep=1
	if (PROJECT$model$timestep == "monthly") {timestep=mean(PROJECT$model$timestep_days)}
	if (PROJECT$model$timestep == "weekly") {timestep=mean(PROJECT$model$timestep_days)}
	time_vector=1:dim(states_all$gpp_gCm2day)[2]
	year_vector=time_vector/(365.25/timestep)
	year_vector=year_vector+as.numeric(PROJECT$start_year)
	interval=floor(length(year_vector)/10)

        if (which_plot == -5) {

            if (PROJECT$ctessel_pft[n] == 1) {
                # These models assume rooting depth is controlled by coarse root, which is a fraction of the woody pool!
                var=t(states_all$root_gCm2)
                # parameter numbers adjusted for crop model
                var=as.vector(parameters[37,,]) * (var*2) / (as.vector(parameters[36,,]) + (var*2))
            } else if (PROJECT$model$name == "DALECN_BUCKET") {
                # Now estimate the rooting depth based on the equation imbedded in DALECN_BUCKET
                var=as.vector(parameters[35,,]) * (var*2) / (as.vector(parameters[34,,]) + (var*2))
            }  else if (PROJECT$model$name == "DALEC_BUCKET" | PROJECT$model$name == "DALEC" |
                PROJECT$model$name == "DALEC_G5" | PROJECT$model$name == "DALEC_G6" |
                PROJECT$model$name == "DALEC_GSI_BUCKET" | PROJECT$model$name == "DALEC_BUCKET_CanAGE"){
                # These models assume rooting depth is controlled by coarse root, which is a fraction of the woody pool!
                tmp = t(states_all$wood_gCm2)*as.vector(parameters[29,,])
                var=t(states_all$root_gCm2) + tmp
                # Now estimate the rooting depth based on the equation imbedded in DALEC_GSI_BUCKET
                var=as.vector(parameters[40,,]) * (var*2) / (as.vector(parameters[39,,]) + (var*2))
            }  else if (PROJECT$model$name == "DALEC_CDEA_ACM2_BUCKET"){
                # These models assume rooting depth is controlled by coarse root, which is a fraction of the woody pool!
                tmp = t(states_all$wood_gCm2)*as.vector(parameters[25,,])
                var=t(states_all$root_gCm2) + tmp
                # Now estimate the rooting depth based on the equation imbedded in DALEC_GSI_BUCKET
                var=as.vector(parameters[27,,]) * (var*2) / (as.vector(parameters[26,,]) + (var*2))
            }  else if (PROJECT$model$name == "DALEC_CDEA_ACM2_BUCKET_RmRg"){
                # These models assume rooting depth is controlled by coarse root, which is a fraction of the woody pool!
                tmp = t(states_all$wood_gCm2)*as.vector(parameters[25,,])
                var=t(states_all$root_gCm2) + tmp
                # Now estimate the rooting depth based on the equation imbedded in DALEC_GSI_BUCKET
                var=as.vector(parameters[27,,]) * (var*2) / (as.vector(parameters[26,,]) + (var*2))
            }  else if (PROJECT$model$name == "DALEC_CDEA_ACM2_BUCKET_RmRg_CWD"){
                # These models assume rooting depth is controlled by coarse root, which is a fraction of the woody pool!
                tmp = t(states_all$wood_gCm2)*as.vector(parameters[25,,])
                var=t(states_all$root_gCm2) + tmp
                # Now estimate the rooting depth based on the equation imbedded in DALEC_GSI_BUCKET
                var=as.vector(parameters[27,,]) * (var*2) / (as.vector(parameters[26,,]) + (var*2))
            } else if (PROJECT$model$name == "DALEC_CDEA_ACM2_BUCKET_RmRg_CWD_wMRT") {
                # These models assume rooting depth is controlled by coarse root, which is a fraction of the woody pool!
                tmp = t(states_all$wood_gCm2)*as.vector(parameters[25,,])
                var=t(states_all$root_gCm2) + tmp
                # Now estimate the rooting depth based on the equation imbedded in DALEC_GSI_BUCKET
                var=as.vector(parameters[27,,]) * (var*2) / (as.vector(parameters[26,,]) + (var*2))
            }

            jpeg(file=paste(PROJECT$figpath,"timeseries_RootDepth_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
            # now create the plotting area
            par(mfrow=c(1,1), mar=c(5,5,3,1))
            plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(var)[which(var != Inf)], prob=c(0.75), na.rm=TRUE)), cex=0.8,ylab="Rooting Depth (m)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
            axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
            # add the confidence intervals
            plotconfidence(var)
            # calculate and draw the median values, could be mean instead or other
            #lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")

            dev.off()

        } else if (which_plot == -4) {

            # incoming data from states_all is dim=c(iter, chain, time)
            # structure needed by function is dim=c(time,iter)

            # flip it to get the right shape
            var=t(states_all$wSWP_MPa)

            jpeg(file=paste(PROJECT$figpath,"timeseries_wSWP_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
            # now create the plotting area
            par(mfrow=c(1,1), mar=c(5,5,3,1))
            plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(quantile(as.vector(var)[which(var != Inf)], prob=c(0.001,0.999), na.rm=TRUE)), cex=0.8,ylab="Weighted Soil Water Potential (MPa)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
            axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
            # add the confidence intervals
            plotconfidence(var)
            # calculate and draw the median values, could be mean instead or other
            #lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")
            dev.off()

        } else if (which_plot == -3) {

            # flip it to get the right shape
            var=t(states_all$sfc_water_mm)

            jpeg(file=paste(PROJECT$figpath,"timeseries_sfc_water_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
            # now create the plotting area
            par(mfrow=c(1,1), mar=c(5,5,3,1))
            plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(var), prob=c(0.999), na.rm=TRUE)), cex=0.8,ylab="Soil surface water (0-30 cm; kgH2O.m-2)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
            axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
            # add the confidence intervals
            plotconfidence(var)
            # calculate and draw the median values, could be mean instead or other
            #lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")
            dev.off()

        } else if (which_plot == -2) {

            if (length(which(is.na(as.vector(states_all$evap_kgH2Om2day)) == FALSE)) > 0) {

                # incoming data from states_all is dim=c(iter, time)
                # structure needed by function is dim=c(time,iter)
                # flip it to get the right shape
                evap_var=t(states_all$evap_kgH2Om2day)

                evap_obs=drivers$obs[,31] ; evap_obs_unc = drivers$obs[,32]
                # filter -9999 to NA
                filter = which(evap_obs == -9999) ; evap_obs[filter] = NA ; evap_obs_unc[filter] = NA

                jpeg(file=paste(PROJECT$figpath,"timeseries_evap_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
                if (PROJECT$model$name == "ACM" & max(as.vector(evap_var),na.rm=TRUE) > 0) {
                    par(mfrow=c(1,1), omi=c(0.1,0.1,0.1,0.1), mai=c(1,1,1,1))
                    maxl=which(as.vector(parameters[dim(parameters)[1],,]) == max(as.vector(parameters[dim(parameters)[1],,])))
                    maxl=maxl[1] # if we happen to have more than one values with the same likelihood we will just pick the first one....
                    plot(states_all$evap_kgH2Om2day[maxl,],evap_obs, ylab="SPA", xlab="ACM", main="Evap (kgH2O.m-2.day-1)",pch=16,cex=0.8,cex.main=1.8,cex.lab=1.8,cex.axis=1.8) ; abline(0,1,col="red", lwd=4)
                    hey=lm(evap_obs~states_all$evap_kgH2Om2day[maxl,]) ; beta1=coef(hey)[2] ; intercept=coef(hey)[1]
                    explained=summary(hey)$adj.r.squared ; trend = mean(states_all$evap_kgH2Om2day[maxl,]-evap_obs) ; error=rmse(evap_obs,states_all$evap_kgH2Om2day[maxl,])
                    prop_error=mean(abs((evap_obs-states_all$evap_kgH2Om2day[maxl,])/evap_obs)*100)
                    text(max(states_all$evap_kgH2Om2day[maxl,])*0.12,max(evap_obs*0.97),label=bquote( SPA == .(round(beta1,3)) * ACM + .(round(intercept,3))), cex=2.)
                    text(max(states_all$evap_kgH2Om2day[maxl,])*0.182,max(evap_obs*0.92),label=paste("R2 = ",round(explained,3)," rmse = ",round(error,3)," bias = ",round(trend,3), sep=""), cex=2.)
                    text(max(states_all$evap_kgH2Om2day[maxl,])*0.074,max(evap_obs*0.87),label=paste("% error = ",round(prop_error,3), sep=""), cex=2.)
                    text(max(states_all$evap_kgH2Om2day[maxl,])*0.146,max(evap_obs*0.81),label=paste("ACM WUE (gC/kgH2O) = ",round(mean(states_all$gpp_gCm2day[maxl,]/states_all$evap_kgH2Om2day[maxl,]),3)," (",round(quantile(states_all$gpp_gCm2day[maxl,]/states_all$evap_kgH2Om2day[maxl,],prob=c(0.025)),3),"/",round(quantile(states_all$gpp_gCm2day[maxl,]/states_all$evap_kgH2Om2day[maxl,],prob=c(0.975)),3),")", sep=""), cex=2.)
                    text(max(states_all$evap_kgH2Om2day[maxl,])*0.143,max(evap_obs*0.76),label=paste("SPA WUE (gC/kgH2O) = ",round(mean(drivers$obs[,1]/evap_obs),3)," (",round(quantile(drivers$obs[,1]/evap_obs,prob=c(0.025)),3),"/",round(quantile(drivers$obs[,1]/evap_obs,prob=c(0.975)),3),")", sep=""), cex=2.)
                } else {
                    # now create the plotting area
                    par(mfrow=c(1,1), mar=c(5,5,3,1))
                    plot(evap_obs, pch=16,xaxt="n", ylim=c(0,quantile(as.vector(evap_var), prob=c(0.999),na.rm=TRUE)), cex=0.8,ylab="Evap (kgH2O.m-2.day-1)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
                    axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
                    # add the confidence intervals
                    plotconfidence(evap_var)
                    # calculate and draw the median values, could be mean instead or other
                    #lines(apply(evap_var[1:(dim(evap_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")
                    # add the data on top if there is any
                    if (length(which(is.na(evap_obs))) != length(evap_obs) ) {
                        points(evap_obs, pch=16, cex=0.8)
                        plotCI(evap_obs,gap=0,uiw=evap_obs_unc, col="black", add=TRUE, cex=1,lwd=2,sfrac=0.01,lty=1,pch=16)
                    }
                } # acm or not
                dev.off()
            }

        } else if (which_plot == 0) {

		# flip it to get the right shape
		Cfol_var=t(states_all$fol_gCm2)

		# pass observations driver
		Cfol_obs=drivers$obs[,11] ; Cfol_obs_unc=drivers$obs[,12]
		# filter -9999 to NA
		filter = which(Cfol_obs == -9999) ; Cfol_obs[filter] = NA ; Cfol_obs_unc[filter] = NA

		# calculate foliarCN
		if (PROJECT$model$name == "DALECN_GSI_FR" | PROJECT$model$name == "DALECN_GSI_DFOL_LABILE_FR" | PROJECT$model$name == "DALECN_GSI_DFOL_LABILE_FROOT_FR") {
			info = round(quantile(as.vector(parameters[17,,]/(10**parameters[11,,])),prob=c(0.025,0.5,0.975)),digits=2)
		}
		ymax=quantile(as.vector(Cfol_var), prob=c(0.999), na.rm=TRUE)
		ymin=quantile(as.vector(Cfol_var), prob=c(0.001), na.rm=TRUE)
		xloc=0.15*dim(Cfol_var)[1] ; yloc=(1-0.05)*ymax

		jpeg(file=paste(PROJECT$figpath,"timeseries_Cfol_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
		# now create the plotting area
		par(mfrow=c(1,1), mar=c(5,5,3,1))
		plot(rep(-9999,dim(Cfol_var)[1]),xaxt="n", pch=16, ylim=c(ymin,ymax), cex=0.8,ylab="Cfol (gC.m-2)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
		axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
		# add the confidence intervals
		plotconfidence(Cfol_var)
		if (PROJECT$model$name == "DALECN_GSI_FR" | PROJECT$model$name == "DALECN_GSI_DFOL_LABILE_FR" | PROJECT$model$name == "DALECN_GSI_DFOL_LABILE_FROOT_FR") {
			text(xloc,yloc, paste("foliarCN = ",info[2],"(",info[1],"/",info[3],")",sep=""),cex=2)
		}
		# calculate and draw the median values, could be mean instead or other
		#lines(apply(Cfol_var[1:(dim(Cfol_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")
		# add the data on top if there is any
		if (length(which(is.na(Cfol_obs))) != length(Cfol_obs) ) {
			points(Cfol_obs, pch=16, cex=0.8)
			plotCI(Cfol_obs,gap=0,uiw=Cfol_obs_unc, col="black", add=TRUE, cex=1,lwd=2,sfrac=0.01,lty=1,pch=16)
		}

		dev.off()

	} else if (which_plot == 1) {

                # flip it to get the right shape
		lai_var=t(states_all$lai_m2m2)
		lai_obs=drivers$obs[,3]

		jpeg(file=paste(PROJECT$figpath,"timeseries_lai_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
		# now create the plotting area
		par(mfrow=c(1,1), mar=c(5,5,3,1))
		plot(lai_obs, pch=16,xaxt="n", ylim=c(0,max(max(lai_obs),quantile(as.vector(lai_var), prob=c(0.999), na.rm=TRUE))), cex=0.8,ylab="LAI (m2/m2)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
		axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
		# add the confidence intervals
		plotconfidence(lai_var)
		# calculate and draw the median values, could be mean instead or other
		#lines(apply(lai_var[1:(dim(lai_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")
		# add the data on top
		points(lai_obs, pch=16, cex=0.8)

		dev.off()

	} else if (which_plot == 2) {

		# flip it to get the right shape
 		gpp_var=t(states_all$gpp_gCm2day)
                gpp_obs=drivers$obs[,1] ; gpp_obs_unc = drivers$obs[,2]
		# filter -9999 to NA
		filter = which(gpp_obs == -9999) ; gpp_obs[filter] = NA ; gpp_obs_unc[filter] = NA

		jpeg(file=paste(PROJECT$figpath,"timeseries_gpp_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
		if (PROJECT$model$name == "ACM") {
			par(mfrow=c(1,1), omi=c(0.1,0.1,0.1,0.1), mai=c(1,1,1,1))
			maxl=which(as.vector(parameters[dim(parameters)[1],,]) == max(as.vector(parameters[dim(parameters)[1],,])))
			maxl=maxl[1] # if we happen to have more than one values with the same likelihood we will just pick the first one....
			plot(states_all$gpp_gCm2day[maxl,],drivers$obs[,1], ylab="SPA", xlab="ACM", main="GPP (gC.m-2.day-1)",pch=16,cex=0.8,cex.main=1.8,cex.lab=1.8,cex.axis=1.8) ; abline(0,1,col="red", lwd=4)
			hey=lm(drivers$obs[,1]~states_all$gpp_gCm2day[maxl,]) ; beta1=coef(hey)[2] ; intercept=coef(hey)[1]
			explained=summary(hey)$adj.r.squared ; trend = mean(states_all$gpp_gCm2day[maxl,]-drivers$obs[,1]) ; error=rmse(drivers$obs[,1],states_all$gpp_gCm2day[maxl,])
			prop_error=mean(abs((drivers$obs[,1]-states_all$gpp_gCm2day[maxl,])/drivers$obs[,1])*100)
			text(max(states_all$gpp_gCm2day[maxl,])*0.12,max(drivers$obs[,1]*0.97),label=bquote( SPA == .(round(beta1,3)) * ACM + .(round(intercept,3))), cex=2.)
			text(max(states_all$gpp_gCm2day[maxl,])*0.182,max(drivers$obs[,1]*0.92),label=paste("R2 = ",round(explained,3)," rmse = ",round(error,3)," bias = ",round(trend,3), sep=""), cex=2.)
			text(max(states_all$gpp_gCm2day[maxl,])*0.074,max(drivers$obs[,1]*0.87),label=paste("% error = ",round(prop_error,3), sep=""), cex=2.)
			text(max(states_all$gpp_gCm2day[maxl,])*0.146,max(drivers$obs[,1]*0.81),label=paste("ACM NUE = ",round(mean(states_all$gpp_gCm2day[maxl,]/drivers$met[,14]),3)," (",round(quantile(states_all$gpp_gCm2day[maxl,]/drivers$met[,14],prob=c(0.025)),3),"/",round(quantile(states_all$gpp_gCm2day[maxl,]/drivers$met[,14],prob=c(0.975)),3),")", sep=""), cex=2.)
			text(max(states_all$gpp_gCm2day[maxl,])*0.143,max(drivers$obs[,1]*0.76),label=paste("SPA NUE = ",round(mean(gpp_obs/drivers$met[,14]),3)," (",round(quantile(gpp_obs/drivers$met[,14],prob=c(0.025)),3),"/",round(quantile(gpp_obs/drivers$met[,14],prob=c(0.975)),3),")", sep=""), cex=2.)
		} else {
			# now create the plotting area
			par(mfrow=c(1,1), mar=c(5,5,3,1))
			plot(gpp_obs, pch=16,xaxt="n", ylim=c(0,quantile(as.vector(gpp_var), prob=c(0.999),na.rm=TRUE)), cex=0.8,ylab="GPP (gC.m-2.day-1)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
			axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
			# add the confidence intervals
			plotconfidence(gpp_var)
			# calculate and draw the median values, could be mean instead or other
			#lines(apply(gpp_var[1:(dim(gpp_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")
			# add the data on top if there is any
			if (length(which(is.na(gpp_obs))) != length(gpp_obs) ) {
				points(gpp_obs, pch=16, cex=0.8)
				plotCI(gpp_obs,gap=0,uiw=gpp_obs_unc, col="black", add=TRUE, cex=1,lwd=2,sfrac=0.01,lty=1,pch=16)
			}
		} # acm or not
		dev.off()

	} else if (which_plot == 3) {

		# flip it to get the right shape
		nee_var=t(states_all$nee_gCm2day)
		# pass observations driver
		nee_obs=drivers$obs[,5] ; nee_obs_unc=drivers$obs[,6]
		# filter -9999 to NA
		filter = which(nee_obs == -9999) ; nee_obs[filter] = NA ; nee_obs_unc[filter] = NA

		jpeg(file=paste(PROJECT$figpath,"timeseries_nee_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
		# now create the plotting area
		par(mfrow=c(1,1), mar=c(5,5,3,1))
		plot(nee_obs, pch=16,xaxt="n", ylim=c(quantile(as.vector(nee_var),prob=c(0.001),na.rm=TRUE),quantile(as.vector(nee_var),prob=c(0.999),na.rm=TRUE)), cex=0.8,ylab="NEE (gC.m-2.day-1)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
		axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
		# add the confidence intervals
		plotconfidence(nee_var)
		# calculate and draw the median values, could be mean instead or other
		#lines(apply(nee_var[1:(dim(nee_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")
		# add the data on top if there is any
		if (length(which(is.na(nee_obs))) != length(nee_obs) ) {
			points(nee_obs, pch=16, cex=0.8)
			if (length(which(is.na(nee_obs)))/length(nee_obs) > 0.01) {
				plotCI(nee_obs,gap=0,uiw=nee_obs_unc, col="black", add=TRUE, cex=1,lwd=2,sfrac=0.01,lty=1,pch=16)
			}
		}

		dev.off()

	} else if (which_plot == 4) {

                # flip it to get the right shape
		Reco_var=t(states_all$reco_gCm2day)
		# pass observations driver
		eco_resp_obs=drivers$obs[,9] ; eco_resp_obs_unc=drivers$obs[,10]
		# filter -9999 to NA
		filter = which(eco_resp_obs == -9999) ; eco_resp_obs[filter] = NA ; eco_resp_obs_unc[filter] = NA

		jpeg(file=paste(PROJECT$figpath,"timeseries_eco_resp_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
		# now create the plotting area
		par(mfrow=c(1,1), mar=c(5,5,3,1))
		plot(rep(-9999,dim(Reco_var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(Reco_var), prob=c(0.999),na.rm=TRUE)), cex=0.8,ylab="Eco Resp (gC.m-2.day-1)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
		axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
		# add the confidence intervals
		plotconfidence(Reco_var)
		# calculate and draw the median values, could be mean instead or other
		#lines(apply(Reco_var[1:(dim(Reco_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")
		# add the data on top if there is any
		if (length(which(is.na(eco_resp_obs))) != length(eco_resp_obs) ) {
			points(eco_resp_obs, pch=16, cex=0.8)
			plotCI(eco_resp_obs,gap=0,uiw=eco_resp_obs_unc, col="black", add=TRUE, cex=1,lwd=2,sfrac=0.01,lty=1,pch=16)
		}
		dev.off()

	} else if (which_plot == 5) {

		# flip it to get the right shape
		Rhet_var=t(states_all$rhet_gCm2day)

		jpeg(file=paste(PROJECT$figpath,"timeseries_het_resp_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
		# now create the plotting area
		par(mfrow=c(1,1), mar=c(5,5,3,1))
		plot(rep(-9999,dim(Rhet_var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(Rhet_var), prob=c(0.999),na.rm=TRUE)), cex=0.8,ylab="Het Resp (gC.m-2.day-1)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
		axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
		# add the confidence intervals
		plotconfidence(Rhet_var)
		# calculate and draw the median values, could be mean instead or other
		#lines(apply(Rhet_var[1:(dim(Rhet_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")

		dev.off()

	} else if (which_plot == 6) {

                # flip it to get the right shape
   		Clab_var=t(states_all$lab_gCm2)

		jpeg(file=paste(PROJECT$figpath,"timeseries_Clab_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
		# now create the plotting area
		par(mfrow=c(1,1), mar=c(5,5,3,1))
		plot(rep(-9999,dim(Clab_var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(Clab_var), prob=c(0.999), na.rm=TRUE)), cex=0.8,ylab="Clab (gC.m-2)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
		axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
		# add the confidence intervals
		plotconfidence(Clab_var)
		# calculate and draw the median values, could be mean instead or other
		#lines(apply(Clab_var[1:(dim(Clab_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")
		dev.off()

		if (grepl("FROOT",PROJECT$model$name) | grepl("LABILE",PROJECT$model$name)) {

			# flip it to get the right shape
			Clab_var=t(states_all$labroot_gCm2)


			jpeg(file=paste(PROJECT$figpath,"timeseries_Clabroot_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
			# now create the plotting area
			par(mfrow=c(1,1), mar=c(5,5,3,1))
			plot(rep(-9999,dim(Clab_var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(Clab_var), prob=c(0.999), na.rm=TRUE)), cex=0.8,ylab="Clab_root (gC.m-2)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
			axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
			# add the confidence intervals
			plotconfidence(Clab_var)
			# calculate and draw the median values, could be mean instead or other
			#lines(apply(Clab_var[1:(dim(Clab_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")
			dev.off()

			# flip it to get the right shape
			Clab_var=t(states_all$labwood_gCm2)

			jpeg(file=paste(PROJECT$figpath,"timeseries_Clabwood_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
			# now create the plotting area
			par(mfrow=c(1,1), mar=c(5,5,3,1))
			plot(rep(-9999,dim(Clab_var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(Clab_var), prob=c(0.999), na.rm=TRUE)), cex=0.8,ylab="Clab_wood (gC.m-2)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
			axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
			# add the confidence intervals
			plotconfidence(Clab_var)
			# calculate and draw the median values, could be mean instead or other
			#lines(apply(Clab_var[1:(dim(Clab_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")
			dev.off()

		}

		if (PROJECT$model$name=="DALEC_GSI_FR_LABILE") {

			# flip it to get the right shape
			Clab_var=t(states_all$Clabslow_gCm2)

			jpeg(file=paste(PROJECT$figpath,"timeseries_Clabslow_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
			# now create the plotting area
			par(mfrow=c(1,1), mar=c(5,5,3,1))
			plot(rep(-9999,dim(Clab_var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(Clab_var), prob=c(0.999), na.rm=TRUE)), cex=0.8,ylab="Clab_slow (gC.m-2)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
			axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
			# add the confidence intervals
			plotconfidence(Clab_var)
			# calculate and draw the median values, could be mean instead or other
			#lines(apply(Clab_var[1:(dim(Clab_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")
			dev.off()

		}

	} else if (which_plot == 7) {

		# flip it to get the right shape
		Clit_var=t(states_all$lit_gCm2)

		# pass observations driver
		Clit_obs=drivers$obs[,17] ; Clit_obs_unc=drivers$obs[,18]
		# filter -9999 to NA
		filter = which(Clit_obs == -9999) ; Clit_obs[filter] = NA ; Clit_obs_unc[filter] = NA

		filename_alt="timeseries_Clit_" ; var_name_alt="Clit (gC.m-2)"
		if (PROJECT$model$name == "DALEC_GSI_DBio_FR") {filename_alt="timeseries_Clitfoliage_" ; var_name_alt="Clit_foliage (gC.m-2)"}

		jpeg(file=paste(PROJECT$figpath,filename_alt,PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)

		# now create the plotting area
		par(mfrow=c(1,1), mar=c(5,5,3,1))
		plot(rep(-9999,dim(Clit_var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(Clit_var[1:(dim(Clit_var)[1]-1),]), prob=c(0.999), na.rm=TRUE)), cex=0.8,ylab=var_name_alt,xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
		axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
		# add the confidence intervals
		plotconfidence(Clit_var)
		# calculate and draw the median values, could be mean instead or other
		#lines(apply(Clit_var[1:(dim(Clit_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")
		# add the data on top if there is any
		if (length(which(is.na(Clit_obs))) != length(Clit_obs) ) {
			points(Clit_obs, pch=16, cex=0.8)
			plotCI(Clit_obs,gap=0,uiw=Clit_obs_unc, col="black", add=TRUE, cex=1,lwd=2,sfrac=0.01,lty=1,pch=16)
		}
		dev.off()

		if (PROJECT$model$name == "DALEC_GSI_DBio_FR") {

			# structure needed by function is dim=c(time,iter)
			# flip it to get the right shape
			Clit_var = Clit_var + t(states_all$litroot_gCm2)

			jpeg(file=paste(PROJECT$figpath,"timeseries_Clit_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
			# now create the plotting area
			par(mfrow=c(1,1), mar=c(5,5,3,1))
			plot(rep(-9999,dim(Clit_var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(Clit_var), prob=c(0.999), na.rm=TRUE)), cex=0.8,ylab="Clit (gC m-2)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
			axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
			# add the confidence intervals
			plotconfidence(Clit_var)
			# calculate and draw the median values, could be mean instead or other
			#lines(apply(Clit_var[1:(dim(Clit_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")

			dev.off()
		}

	} else if (which_plot == 8) {

                # flip it to get the right shape
		Cr_var=t(states_all$root_gCm2)

		# pass all observations (if any)
		# current defaults are:
		# 1st column = GPP
		# 2nd column = LAI
		# 3rd column = NEE
		# pass observations driver
		Croots_obs=drivers$obs[,15] ; Croots_obs_unc=drivers$obs[,16]
		# filter -9999 to NA
		filter = which(Croots_obs == -9999) ; Croots_obs[filter] = NA ; Croots_obs_unc[filter] = NA

		jpeg(file=paste(PROJECT$figpath,"timeseries_Croots_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
		# now create the plotting area
		par(mfrow=c(1,1), mar=c(5,5,3,1))
		plot(rep(-9999,dim(Cr_var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(Cr_var), prob=c(0.999), na.rm=TRUE)), cex=0.8,ylab="Croots (gC.m-2)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
		axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
		# add the confidence intervals
		plotconfidence(Cr_var)
		# calculate and draw the median values, could be mean instead or other
		#lines(apply(Cr_var[1:(dim(Cr_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")
		# add the data on top if there is any
		if (length(which(is.na(Croots_obs))) != length(Croots_obs) ) {
			points(Croots_obs, pch=16, cex=0.8)
			plotCI(Croots_obs,gap=0,uiw=Croots_obs_unc, col="black", add=TRUE, cex=1,lwd=2,sfrac=0.01,lty=1,pch=16)
		}
		dev.off()

	} else if (which_plot == 9) {

                # flip it to get the right shape
		Cw_var=t(states_all$wood_gCm2)
		# pass observations driver
		Cwood_obs=drivers$obs[,13] ; Cwood_obs_unc=drivers$obs[,14]
		# filter -9999 to NA
		filter = which(Cwood_obs == -9999) ; Cwood_obs[filter] = NA ; Cwood_obs_unc[filter] = NA

		jpeg(file=paste(PROJECT$figpath,"timeseries_Cwood_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
		# now create the plotting area
		par(mfrow=c(1,1), mar=c(5,5,3,1))
		plot(rep(-9999,dim(Cw_var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(Cw_var), prob=c(0.999), na.rm=TRUE)), cex=0.8,ylab="Cwood (gC.m-2)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
		axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
		# add the confidence intervals
		plotconfidence(Cw_var)
		# calculate and draw the median values, could be mean instead or other
		#lines(apply(Cw_var[1:(dim(Cw_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")
		# add the data on top if there is any
		if (length(which(is.na(Cwood_obs))) != length(Cwood_obs) ) {
			points(Cwood_obs, pch=16, cex=0.8)
			plotCI(Cwood_obs,gap=0,uiw=Cwood_obs_unc, col="black", add=TRUE, cex=1,lwd=2,sfrac=0.01,lty=1,pch=16)
		}
		dev.off()

	} else if (which_plot == 10) {

		# flip it to get the right shape
		Csom_var=t(states_all$som_gCm2)

		# pass observations driver
		Csom_obs=drivers$obs[,19] ; Csom_obs_unc=drivers$obs[,20]
		# filter -9999 to NA
		filter = which(Csom_obs == -9999) ; Csom_obs[filter] = NA ; Csom_obs_unc[filter] = NA

		jpeg(file=paste(PROJECT$figpath,"timeseries_Csom_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
		# now create the plotting area
		par(mfrow=c(1,1), mar=c(5,5,3,1))
		plot(rep(-9999,dim(Csom_var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(Csom_var), prob=c(0.999), na.rm=TRUE)), cex=0.8,ylab="Csom (gC.m-2)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
		axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
		# add the confidence intervals
		plotconfidence(Csom_var)
		# calculate and draw the median values, could be mean instead or other
		#lines(apply(Csom_var[1:(dim(Csom_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")
		# add the data if there is any which is not missing
		if (length(which(is.na(Csom_obs))) != length(Csom_obs) ) {
			# add the data on top
			points(Csom_obs, pch=16, cex=0.8)
			plotCI(Csom_obs,gap=0,uiw=Csom_obs_unc, col="black", add=TRUE, cex=1,lwd=2,sfrac=0.01,lty=1,pch=16)
		}
		dev.off()

	} else if (which_plot == 11) {

		# structure needed by function is dim=c(time,iter)
		# flip it to get the right shape
		Call_var=t(states_all$bio_gCm2)

		jpeg(file=paste(PROJECT$figpath,"timeseries_Cbio_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
		# now create the plotting area
		par(mfrow=c(1,1), mar=c(5,5,3,1))
		plot(rep(-9999,dim(Call_var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(Call_var), prob=c(0.999), na.rm=TRUE)), cex=0.8,ylab="Cbio (gC.m-2)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
		axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
		# add the confidence intervals
		plotconfidence(Call_var)
		# calculate and draw the median values, could be mean instead or other
		#lines(apply(Call_var[1:(dim(Call_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")

		dev.off()

	} else if (which_plot == 12) {

		if (PROJECT$ctessel_pft[n] != 1) {
			if (grepl("DALECN",PROJECT$model$name) || grepl("DFOL",PROJECT$model$name) || grepl("DBio",PROJECT$model$name) || grepl("BUCKET",PROJECT$model$name)) {
				# load timesteps to local variable
				timestep_days=PROJECT$model$timestep_days
				if (length(timestep_days) == 1) {
					timestep_days=rep(timestep_days, times=dim(states_all[[1]])[2])
				}
				# determine averaging period for timeseries
				tmp_m=rep(0, times=length(timestep_days))
				for (step in seq(1,length(timestep_days))) {
					# calculate the gradient / trend of GSI
					if (sum(timestep_days[1:step]) < 21) {
						tmp_m[step] = step-1
					} else {
						# else we will try and work out the gradient to see what is
						# happening
						# to the system over all. The default assumption will be to
						# consider
						# the averaging period of GSI model (i.e. 21 days). If this is not
						# possible either the time step of the system is used (if step
						# greater
						# than 21 days) or all available steps (if n < 21).
						m = 0 ; test = 0
						while (test < 21) {
							m=m+1 ; test = sum(timestep_days[(step-m):step])
							if (m == (step-1)) {test = 21}
						} # while
						tmp_m[step] = m
					} # for calculating gradient
				} # timestep_days
			} # *DFOL*

			# structure needed by function is dim=c(time,iter)
			# flip it to get the right shape
			GSI_var=t(states_all$gsi)
			# bulk value
			jpeg(file=paste(PROJECT$figpath,"timeseries_GSI_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
			# now create the plotting area
			par(mfrow=c(1,1), mar=c(5,5,3,1))
			plot(rep(-9999,dim(GSI_var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(GSI_var), prob=c(0.999), na.rm=TRUE)), cex=0.8,ylab="GSI",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
			axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=1),tck=-0.02, padj=+0.15, cex.axis=1.9)
			# add the confidence intervals
			plotconfidence(GSI_var)
			# calculate and draw the median values, could be mean instead or other
			#lines(apply(GSI_var[1:(dim(GSI_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")
			dev.off()

			# structure needed by function is dim=c(time,iter)
			# flip it to get the right shape
			GSI_var=t(states_all$gsi_iphoto)

			bob=rep(0, times=3)
#			if (grepl("DALECN",PROJECT$model$name) || grepl("DFOL",PROJECT$model$name) || grepl("DBio",PROJECT$model$name) || grepl("BUCKET",PROJECT$model$name)) {
#				gsi_input = array(0,dim=c(dim(states_all[[1]])[1],dim(states_all[[1]])[2],3))
#                                gsi_input[1:dim(states_all[[16]])[1],,1] = states_all$gsi_itemp ; gsi_input[1:dim(states_all[[16]])[1],,2]=states_all$gsi_iphoto
#                                gsi_input[1:dim(states_all[[16]])[1],,3] = states_all$gsi_ivpd
#           			# now calculate the time when the GSI value is having the greatest impact on the gradient that determines leaf growth or senescence
#				bob = gsi_controlling(gsi_input,states_all$gsi,2,timestep_days,tmp_m)
#				bob = quantile(bob[,2],prob=c(0.975,0.50,0.025),na.rm=TRUE)
#			}

			# photoperiod
			jpeg(file=paste(PROJECT$figpath,"timeseries_GSI_photoperiod_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
			# now create the plotting area
			par(mfrow=c(1,1), mar=c(5,5,3,1))
			plot(rep(-9999,dim(GSI_var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(GSI_var), prob=c(0.999), na.rm=TRUE)), cex=0.8,ylab="GSI",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name," ( median = ",round(bob[2],digits=2),", 97.5% = ",round(bob[1],digits=2),", 2.5% = ",round(bob[3],digits=2),")", sep=""))
			axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=1),tck=-0.02, padj=+0.15, cex.axis=1.9)
			# add the confidence intervals
			plotconfidence(GSI_var)
			# calculate and draw the median values, could be mean instead or other
			#lines(apply(GSI_var[1:(dim(GSI_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")
			dev.off()

			# structure needed by function is dim=c(time,iter)
			# flip it to get the right shape
			GSI_var=t(states_all$gsi_itemp)
#			if (grepl("DALECN",PROJECT$model$name) || grepl("DFOL",PROJECT$model$name) || grepl("DBio",PROJECT$model$name) || grepl("BUCKET",PROJECT$model$name)) {
#				# now calculate the time when the GSI value is having the greatest impact on the gradient that determines leaf growth or senescence
#				bob = gsi_controlling(gsi_input,states_all$gsi,1,timestep_days,tmp_m)
#				bob = quantile(bob[,1],prob=c(0.975,0.50,0.025),na.rm=TRUE)
#			}
			# temperature
			jpeg(file=paste(PROJECT$figpath,"timeseries_GSI_temperature_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
			# now create the plotting area
			par(mfrow=c(1,1), mar=c(5,5,3,1))
			plot(rep(-9999,dim(GSI_var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(GSI_var), prob=c(0.999), na.rm=TRUE)), cex=0.8,ylab="GSI",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name," ( median = ",round(bob[2],digits=2),", 97.5% = ",round(bob[1],digits=2),", 2.5% = ",round(bob[3],digits=2),")", sep=""))
			axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=1),tck=-0.02, padj=+0.15, cex.axis=1.9)
			# add the confidence intervals
			plotconfidence(GSI_var)
			# calculate and draw the median values, could be mean instead or other
			#lines(apply(GSI_var[1:(dim(GSI_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")
			dev.off()

			# structure needed by function is dim=c(time,iter)
			# flip it to get the right shape
			GSI_var=t(states_all$gsi_ivpd)
#			if (grepl("DALECN",PROJECT$model$name) || grepl("DFOL",PROJECT$model$name) || grepl("DBio",PROJECT$model$name) || grepl("BUCKET",PROJECT$model$name)) {
#				# now calculate the time when the GSI value is having the greatest impact on the gradient that determines leaf growth or senescence
#				bob = gsi_controlling(gsi_input,states_all$gsi,3,timestep_days,tmp_m)
#				bob = quantile(bob[,3],prob=c(0.975,0.50,0.025),na.rm=TRUE)
#			}
			# VPD
			jpeg(file=paste(PROJECT$figpath,"timeseries_GSI_vpd_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
			# now create the plotting area
			par(mfrow=c(1,1), mar=c(5,5,3,1))
			plot(rep(-9999,dim(GSI_var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(GSI_var), prob=c(0.999), na.rm=TRUE)), cex=0.8,ylab="GSI",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name," ( median = ",round(bob[2],digits=2),", 97.5% = ",round(bob[1],digits=2),", 2.5% = ",round(bob[3],digits=2),")", sep=""))
			axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=1),tck=-0.02, padj=+0.15, cex.axis=1.9)
			# add the confidence intervals
			plotconfidence(GSI_var)
			# calculate and draw the median values, could be mean instead or other
			#lines(apply(GSI_var[1:(dim(GSI_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")
			dev.off()
		} #pft not == 1 (crop)
	} else if (which_plot == 13) {

		# structure needed by function is dim=c(time,iter)
		# flip it to get the right shape
		somfast_var=t(states_all$somfast_gCm2)

		jpeg(file=paste(PROJECT$figpath,"timeseries_Csomfast_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
		# now create the plotting area
		par(mfrow=c(1,1), mar=c(5,5,3,1))
		plot(rep(-9999,dim(somfast_var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(somfast_var), prob=c(0.999), na.rm=TRUE)), cex=0.8,ylab="Csom_fast",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
		axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
		# add the confidence intervals
		plotconfidence(somfast_var)
		# calculate and draw the median values, could be mean instead or other
		#lines(apply(somfast_var[1:(dim(somfast_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")
		dev.off()

		# structure needed by function is dim=c(time,iter)
		# flip it to get the right shape
		somfast_var=t(states_all$microbial_gCm2)

		jpeg(file=paste(PROJECT$figpath,"timeseries_Cmicrobial_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
		# now create the plotting area
		par(mfrow=c(1,1), mar=c(5,5,3,1))
		plot(rep(-9999,dim(somfast_var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(somfast_var), prob=c(0.999), na.rm=TRUE)), cex=0.8,ylab="Cmicrobial",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
		axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
		# add the confidence intervals
		plotconfidence(somfast_var)
		# calculate and draw the median values, could be mean instead or other
		#lines(apply(somfast_var[1:(dim(somfast_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")
		dev.off()

	} else if (which_plot == 14) {

		# structure needed by function is dim=c(time,iter)
		# flip it to get the right shape
		litroot_var=t(states_all$litroot_gCm2)

		jpeg(file=paste(PROJECT$figpath,"timeseries_Crootlitter_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
		# now create the plotting area
		par(mfrow=c(1,1), mar=c(5,5,3,1))
		plot(rep(-9999,dim(litroot_var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(litroot_var), prob=c(0.999), na.rm=TRUE)), cex=0.8,ylab="Clitter_root",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
		axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
		# add the confidence intervals
		plotconfidence(litroot_var)
		# calculate and draw the median values, could be mean instead or other
		#lines(apply(litroot_var[1:(dim(litroot_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")

		dev.off()

	} else if (which_plot == 15) {

		# structure needed by function is dim=c(time,iter)
		# flip it to get the right shape
		litwood_var=t(states_all$litwood_gCm2)

		jpeg(file=paste(PROJECT$figpath,"timeseries_Cwoodlitter_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
		# now create the plotting area
		par(mfrow=c(1,1), mar=c(5,5,3,1))
		plot(rep(-9999,dim(litwood_var)[1]),xaxt="n", pch=16, ylim=c(0,quantile(as.vector(litwood_var), prob=c(0.999), na.rm=TRUE)), cex=0.8,ylab="Clitter_wood",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
		axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
		# add the confidence intervals
		plotconfidence(litwood_var)
		# calculate and draw the median values, could be mean instead or other
		#lines(apply(litwood_var[1:(dim(litwood_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")

		dev.off()

	} else if (which_plot == 17) {

		# structure needed by function is dim=c(time,iter)
		# flip it to get the right shape
		litN_var=t(states_all$litN_gNm2)

		# calculate mean litCN
		info = round(quantile(apply((t(states_all$lit_gCm2)/litN_var),2,mean),prob=c(0.025,0.5,0.975)),digits=2)
		ymax=quantile(as.vector(litN_var), prob=c(0.999), na.rm=TRUE)
		xloc=0.15*dim(litN_var)[1] ; yloc=(1-0.05)*ymax
		jpeg(file=paste(PROJECT$figpath,"timeseries_Nlitter_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
		# now create the plotting area
		par(mfrow=c(1,1), mar=c(5,5,3,1))
		plot(rep(-9999,dim(litN_var)[1]),xaxt="n", pch=16, ylim=c(0,ymax), cex=0.8,ylab="Nlitter",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
		axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
		text(xloc,yloc, paste("mean litC:N = ",info[2],"(",info[1],"/",info[3],")",sep=""),cex=2)
		# add the confidence intervals
		plotconfidence(litN_var)
		# calculate and draw the median values, could be mean instead or other
		#lines(apply(litN_var[1:(dim(litN_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")

		dev.off()

	} else if (which_plot == 18) {

		# structure needed by function is dim=c(time,iter)
		# flip it to get the right shape
		labN_var=t(states_all$labN_gNm2)

		ymax=quantile(as.vector(labN_var), prob=c(0.999), na.rm=TRUE)
		jpeg(file=paste(PROJECT$figpath,"timeseries_Nlabile_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
		# now create the plotting area
		par(mfrow=c(1,1), mar=c(5,5,3,1))
		plot(rep(-9999,dim(labN_var)[1]),xaxt="n", pch=16, ylim=c(0,ymax), cex=0.8,ylab="Nlabile",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
		axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
		# add the confidence intervals
		plotconfidence(labN_var)
		# calculate and draw the median values, could be mean instead or other
		#lines(apply(labN_var[1:(dim(labN_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")

		dev.off()

	} else if (which_plot == 19) {

		# structure needed by function is dim=c(time,iter)
		# flip it to get the right shape
		DIN_var=t(states_all$DIN_gNm2)
		# calculate mean DIN
		info = round(quantile(apply(DIN_var,2,mean),prob=c(0.025,0.5,0.975)),digits=2)
		ymax=quantile(as.vector(DIN_var), prob=c(0.999), na.rm=TRUE)
		xloc=0.15*dim(DIN_var)[1] ; yloc=(1-0.05)*ymax

		jpeg(file=paste(PROJECT$figpath,"timeseries_DIN_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
		# now create the plotting area
		par(mfrow=c(1,1), mar=c(5,5,3,1))
		plot(rep(-9999,dim(DIN_var)[1]),xaxt="n", pch=16, ylim=c(0,ymax), cex=0.8,ylab="Dissolved Inorganic Nitrogen",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
		axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
		# add the confidence intervals
		plotconfidence(DIN_var)
		text(xloc,yloc, paste("mean DIN = ",info[2],"(",info[1],"/",info[3],")",sep=""),cex=2)
		# calculate and draw the median values, could be mean instead or other
		#lines(apply(DIN_var[1:(dim(DIN_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")

		dev.off()

	} else if (which_plot == 20) {

		# structure needed by function is dim=c(time,iter)
		# flip it to get the right shape
		DIN_var=t(states_all$N_mineralisation_gNm2day)
		# calculate mean DIN
		info = round(quantile(apply(DIN_var,2,sum)/(dim(DIN_var)[1]/(365.25/timestep)),prob=c(0.025,0.5,0.975)),digits=2)
		ymax=quantile(as.vector(DIN_var), prob=c(0.999), na.rm=TRUE)
		ymin=quantile(as.vector(DIN_var), prob=c(0.001), na.rm=TRUE)
		xloc=0.15*dim(DIN_var)[1] ; yloc=(1-0.05)*ymax

		jpeg(file=paste(PROJECT$figpath,"timeseries_Nmineralisation_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
		# now create the plotting area
		par(mfrow=c(1,1), mar=c(5,5,3,1))
		plot(rep(-9999,dim(DIN_var)[1]),xaxt="n", pch=16, ylim=c(ymin,ymax), cex=0.8,ylab="N mineralisation rate",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
		axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
		# add the confidence intervals
		plotconfidence(DIN_var)
		text(xloc,yloc, paste("mean annual N mineralisation = ",info[2],"(",info[1],"/",info[3],")",sep=""),cex=2)
		# calculate and draw the median values, could be mean instead or other
		#lines(apply(DIN_var[1:(dim(DIN_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")

		dev.off()

	} else if (which_plot == 21) {

		# structure needed by function is dim=c(time,iter)
		# flip it to get the right shape
		DIN_var=states_all$rauto_gCm2day/states_all$gpp_gCm2day
		DIN_var[which(states_all$gpp_gCm2day == 0)] = NA
		DIN_var=t(DIN_var)

		# calculate meanRa:GPP
		info = round(quantile(apply(DIN_var,2,mean,na.rm=TRUE),prob=c(0.025,0.5,0.975)),digits=2)
		ymax=quantile(as.vector(DIN_var), prob=c(0.99), na.rm=TRUE)
		ymin=quantile(as.vector(DIN_var), prob=c(0.01), na.rm=TRUE)
		xloc=0.15*dim(DIN_var)[1] ; yloc=(1-0.05)*ymax

		jpeg(file=paste(PROJECT$figpath,"timeseries_RaGPP_actual_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
		# now create the plotting area
		par(mfrow=c(1,1), mar=c(5,5,3,1))
		plot(rep(-9999,dim(DIN_var)[1]),xaxt="n", pch=16, ylim=c(ymin,ymax), cex=0.8,ylab="Actual Ra:GPP",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
		axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
		# add the confidence intervals
		plotconfidence(DIN_var)
		text(xloc,yloc, paste("mean Ra:GPP = ",info[2],"(",info[1],"/",info[3],")",sep=""),cex=2)
		# calculate and draw the median values, could be mean instead or other
		#lines(apply(DIN_var[1:(dim(DIN_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")

		dev.off()

	} else if (which_plot == 22) {

		# structure needed by function is dim=c(time,iter)
		# flip it to get the right shape
		harvestC_var=t(states_all$harvest_C_gCm2day)

		ymax=quantile(as.vector(harvestC_var), prob=c(0.999), na.rm=TRUE)
		jpeg(file=paste(PROJECT$figpath,"timeseries_harvestedC_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
		# now create the plotting area
		par(mfrow=c(1,1), mar=c(5,5,3,1))
		plot(rep(-9999,dim(harvestC_var)[1]),xaxt="n", pch=16, ylim=c(0,ymax), cex=0.8,ylab="Harvested C",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
		axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
		# add the confidence intervals
		plotconfidence(harvestC_var)
		# calculate and draw the median values, could be mean instead or other
		#lines(apply(harvestC_var[1:(dim(harvestC_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")

		dev.off()

	} else if (which_plot == 23) {

		# structure needed by function is dim=c(time,iter)
		# flip it to get the right shape
		canopyage_var=t(states_all$canopyage_days)

		ymax=quantile(as.vector(canopyage_var), prob=c(0.999), na.rm=TRUE)
		jpeg(file=paste(PROJECT$figpath,"timeseries_CanopyAge_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
		# now create the plotting area
		par(mfrow=c(1,1), mar=c(5,5,3,1))
		plot(rep(-9999,dim(canopyage_var)[1]),xaxt="n", pch=16, ylim=c(0,ymax), cex=0.8,ylab="Mean Canopy Age (days)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
		axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
		# add the confidence intervals
		plotconfidence(canopyage_var)
		# calculate and draw the median values, could be mean instead or other
		#lines(apply(canopyage_var[1:(dim(canopyage_var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")

		dev.off()

	} else if (which_plot == 24) {

		# structure needed by function is dim=c(time,iter)
		# flip it to get the right shape
		var=t(states_all$fire_gCm2day)

		ymax=quantile(as.vector(var), prob=c(0.999), na.rm=TRUE)
		jpeg(file=paste(PROJECT$figpath,"timeseries_fire_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=280, quality=100)
		# now create the plotting area
		par(mfrow=c(1,1), mar=c(5,5,3,1))
		plot(rep(-9999,dim(var)[1]),xaxt="n", pch=16, ylim=c(0,ymax), cex=0.8,ylab="Fire (gC/m2/day)",xlab="Time (Year)", cex.lab=1.8, cex.axis=1.8, cex.main=1.8, main=paste(PROJECT$sites[n]," - ",PROJECT$name, sep=""))
		axis(1, at=time_vector[seq(1,length(time_vector),interval)],labels=round(year_vector[seq(1,length(time_vector),interval)], digits=0),tck=-0.02, padj=+0.15, cex.axis=1.9)
		# add the confidence intervals
		plotconfidence(var)
		# calculate and draw the median values, could be mean instead or other
		#lines(apply(var[1:(dim(var)[1]-1),],1,median,na.rm=TRUE), pch=1, col="blue")

		dev.off()

	} else {
		print("have requested a figure which we have not actually scripted yet")
	}# end of conditional

	# tidy before leaving
	gc(reset=TRUE, verbose=FALSE)

} # end of function
## Use byte compile
uncertainty_figures<-cmpfun(uncertainty_figures)
