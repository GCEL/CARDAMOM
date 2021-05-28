
###
## Function to create and save plots of CARDAMOM parameter vectors
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

plot_parameters<- function(PROJECT,parameters,converged,n) {

      # input is order for parameters dimensions(npar+1,iter,chain)

      # now check out the pdfs
      # construct the parameter name info
      character_bit=rep("p",times=(dim(parameters)[1]-1))
      number_bit=1:(dim(parameters)[1]-1)
      # merge the letter and numbers together
      par_names=c(paste(character_bit,number_bit,sep=""),"log-likelihood")
      # now add whether the parameter has converged
      par_names=paste(par_names," (",converged,") ",sep="")

      jpeg(file=paste(PROJECT$figpath,"random_walk_of_parameters_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=400, quality=100)
      if (PROJECT$parameter_type == "pft_specific" & PROJECT$ctessel_pft[n] == 1) {
          par(mfrow=c(7,7),mar=c(3, 3, 3, 1), oma=c(0,0,1,0))
      } else if (PROJECT$model$name == "ACM") {
          par(mfrow=c(5,6),mar=c(5, 4, 4, 2), oma=c(0,0,1,0))
      } else if (PROJECT$model$name == "DALEC_GSI_DBio_FR") {
          par(mfrow=c(7,8),mar=c(5, 4, 4, 2), oma=c(0,0,1,0))
      } else if (PROJECT$model$name == "DALECN_GSI_FR" ) {
          par(mfrow=c(7,8),mar=c(5, 4, 4, 2), oma=c(0,0,1,0))
      } else if (PROJECT$model$name == "DALECN_GSI_DFOL_LABILE_FR") {
          par(mfrow=c(8,8),mar=c(3, 1.5, 1, 1), oma=c(0,0,1,0))
      } else if (PROJECT$model$name == "DALECN_GSI_DFOL_LABILE_FROOT_FR") {
          par(mfrow=c(8,8),mar=c(3, 1.5, 1, 1), oma=c(0,0,1,0))
      } else {
          par(mfrow=c(7,7),mar=c(2.4, 3, 3.75, 1), oma=c(0,0,1,0))
      }
      for (i in seq(1,dim(parameters)[1])) {
           plot(as.vector(parameters[i,,]), main=par_names[i],cex.lab=1.4, cex.axis=1.4, cex.main=1.4, ylab="", xlab="")
           if (i == 3) {
                mtext(paste(PROJECT$sites[n]," ",PROJECT$name), padj=-2.3, cex=1.4)
           }
      }
      dev.off()

      jpeg(file=paste(PROJECT$figpath,"/histogram_of_parameters_",PROJECT$sites[n],"_",PROJECT$name,".jpeg",sep=""), width=7200, height=4000, res=400, quality=100)
      if (PROJECT$parameter_type == "pft_specific" & PROJECT$ctessel_pft[n] == 1) {
          par(mfrow=c(7,7),mar=c(3, 3, 3, 1), oma=c(0,0,1,0))
      } else if (PROJECT$model$name == "ACM") {
          par(mfrow=c(5,6),mar=c(5, 4, 4, 2), oma=c(0,0,1,0))
      } else if (PROJECT$model$name == "DALEC_GSI_DBio_FR") {
          par(mfrow=c(7,8),mar=c(5, 4, 4, 2), oma=c(0,0,1,0))
      } else if (PROJECT$model$name == "DALECN_GSI_FR") {
          par(mfrow=c(7,8),mar=c(5, 4, 4, 2), oma=c(0,0,1,0))
      } else if (PROJECT$model$name == "DALECN_GSI_DFOL_LABILE_FR") {
          par(mfrow=c(8,8),mar=c(3, 1.5, 1, 1), oma=c(0,0,1,0))
      } else if (PROJECT$model$name == "DALECN_GSI_DFOL_LABILE_FROOT_FR") {
          par(mfrow=c(8,8),mar=c(3, 1.5, 1, 1), oma=c(0,0,1,0))
      } else {
          par(mfrow=c(7,7),mar=c(2.4, 3, 3.75, 1), oma=c(0,0,1,0))
      }
      for (i in seq(1,dim(parameters)[1])) {
           hist(as.vector(parameters[i,,]), main=par_names[i], cex.lab=1.4, cex.axis=1.4, cex.main=1.4, ylab="", xlab="")
           if (i == 3) {
               mtext(paste(PROJECT$sites[n]," ",PROJECT$name), padj=-2.3, cex=1.4)
           }
      }
      dev.off()

} # function end
