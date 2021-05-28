
###
## Function to generate mean state variable information by running the vs and model choice
## NOTE: this file is no longer in use
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

simulate_statevars<- function (model_name,met,pars,lat) {

      noedc=100

      # update user
      print("Beginning simulation to determine mean state variables")

      # declare output variables
      # order is npar,iter,chain
      laiall=array(0,dim=c(dim(pars)[2],dim(pars)[3]))
      gppall=array(0,dim=c(dim(pars)[2],dim(pars)[3]))
      neeall=array(0,dim=c(dim(pars)[2],dim(pars)[3]))
      somallin=array(0,dim=c(dim(pars)[2],dim(pars)[3]))
      somallout=array(0,dim=c(dim(pars)[2],dim(pars)[3]))
      agball=array(0,dim=c(dim(pars)[2],dim(pars)[3]))
      somall=array(0,dim=c(dim(pars)[2],dim(pars)[3]))
      bioall=array(0,dim=c(dim(pars)[2],dim(pars)[3]))
      rooall=array(0,dim=c(dim(pars)[2],dim(pars)[3]))
      litall=array(0,dim=c(dim(pars)[2],dim(pars)[3]))
      laball=array(0,dim=c(dim(pars)[2],dim(pars)[3]))
      folall=array(0,dim=c(dim(pars)[2],dim(pars)[3]))
      # 100 gives plenty of space for lots of ecological and dynamical constraints
      edcdall=array(0,dim=c(dim(pars)[2],dim(pars)[3],100))

      # loop through combinations
      for (n in seq(1,dim(pars)[2])) {
           # loop through chains
           for (c in seq(1, dim(pars)[3])) {
                if (model_name == "DALEC_CDEA") {
                    output=dalec_cdea_r(met,pars[,n,c],lat)
                } else {
                    stop(paste("Model choice (",model_name,") does not have corresponding R version",sep=""))
                }
                # calculate means
                laiall[n,c]=mean(output$lai)
                gppall[n,c]=mean(output$gpp)
                neeall[n,c]=mean(output$nee)
                somallin[n,c]=mean((output$litter2som+output$woodlitter_production))
                somallout[n,c]=mean(output$respiration_het_som)
                agball[n,c]=mean(output$C_wood)
                somall[n,c]=mean(output$C_som)
                bioall[n,c]=mean(output$Biomass)
                rooall[n,c]=mean(output$C_root)
                litall[n,c]=mean(output$C_litter)
                laball[n,c]=mean(output$C_labile)
                folall[n,c]=mean(output$C_foliar)
#                edcdall[n,c]=mean(output$edcdiags)
	              if (n%%100 == 0 & c == 1) {print(paste(round((n/dim(pars)[2])*100,0)," % completed of mean state calculations",sep=""))}
           } # loop for chains
      } # loop for combinations

      # combine output
      mean_states=list(lai=laiall
                      ,gpp=gppall
                      ,nee=neeall
                      ,somin=somallin
                      ,somout=somallout
                      ,reco=neeall+gppall
                      ,agb=agball
                      ,som=somall
                      ,bio=bioall
                      ,root=rooall
                      ,lit=litall
                      ,lab=laball
                      ,fol=folall)

      # return state variable means
      return(mean_states)

} # end of function
