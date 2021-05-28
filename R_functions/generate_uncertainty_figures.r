
###
## Function to direct the creation of site level timeseries plots of CARDAMOM analyes
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

generate_uncertainty_figures<-function(PROJECT,n) {

	# load locally needed library
	require(gplots)

	# generate file name of the output file created in stage 3
	loadfile=paste(PROJECT$results_processedpath,PROJECT$sites[n],".RData",sep="")

	if (file.exists(loadfile) == TRUE) {
	    #stime=proc.time()["elapsed"]
	    load(loadfile) ; print(paste("DALEC simulations will be loaded from ",loadfile,sep=""))
	    #print(paste("load dalec in ",proc.time()["elapsed"]-stime," seconds",sep=""))
	} else {
	    # do we run the parameters yet for analysis
	    run_all=readline("Raw results have not been processed therefore we will do it now. Do you want to run all parameter vectors to generate confidence intervals? (y/n)")
	    if (run_all == "y") {
          PROJECT$latter_sample_frac = 0.5 # readline("What (latter) fraction of accepted parameters to use (e.g. 0.5)?")
          # Run the parameters
          run_mcmc_results(PROJECT,stage,repair,grid_override)
          # and read in the results
          load(loadfile)
	    } # if condition
	} # file exists statement

	# how many plots in total do we have
	nos_plots=0:11
  if (PROJECT$model$name == "DALEC_CDEA_ACM2_BUCKET") {nos_plots=c(-5,-4,-3,-2,nos_plots,24)}
  if (PROJECT$model$name == "DALEC_CDEA_ACM2_BUCKET_RmRg") {nos_plots=c(-5,-4,-3,-2,nos_plots,24)}
  if (PROJECT$model$name == "DALEC_CDEA_ACM2_BUCKET_RmRg_CWD") {nos_plots=c(-5,-4,-3,-2,nos_plots,15,24)}
  if (PROJECT$model$name == "DALEC_CDEA_ACM2_BUCKET_RmRg_CWD_wMRT") {nos_plots = c(-5,-4,-3,-2,nos_plots,15,24)}
  if (PROJECT$model$name == "DALEC_GSI_BUCKET") {nos_plots=c(-5,-4,-3,-2,nos_plots,12,15,22,24)}
  if (PROJECT$model$name == "DALEC_BUCKET") {nos_plots=c(-5,-4,-3,-2,nos_plots,12,15,22,24)}
  if (PROJECT$model$name == "DALEC_G5") {nos_plots=c(-5,-4,-3,-2,nos_plots,12,15,22,24)}
  if (PROJECT$model$name == "DALEC_G6") {nos_plots=c(-5,-4,-3,-2,nos_plots,15,22,24)}
  if (PROJECT$model$name == "DALEC_BUCKET_CanAGE") {nos_plots=c(-5,-4,-3,-2,nos_plots,12,15,22,24)}
  if (PROJECT$model$name == "DALEC") {nos_plots=c(-5,nos_plots,12,15,22,24)}
  if (PROJECT$model$name == "DALECN_GSI_BUCKET") {nos_plots=c(-5,-4,-3,-2,nos_plots,12,15,21,22)}
  if (PROJECT$model$name == "DALECN_BUCKET") {nos_plots=c(-5,-4,-3,-2,nos_plots,15,21,22,23)}
  if (PROJECT$model$name == "DALEC_GSI_DFOL_CWD_FR") {nos_plots=c(nos_plots,12,15,22,24)}
  if (PROJECT$model$name == "DALEC_GSI_DBio_FR") {nos_plots=0:16}
  if (PROJECT$model$name == "DALECN_GSI_FR" | PROJECT$model$name == "DALECN_GSI_DFOL_LABILE_FR" | PROJECT$model$name == "DALECN_GSI_DFOL_LABILE_FROOT_FR") {nos_plots=c(0:12,15,17:21)}
	if (PROJECT$model$name == "ACM") {nos_plots=c(2,-2)}
	# now request the creation of the plots
	if (use_parallel & length(nos_plots) > 1) {
	    cl <- makeCluster(min(length(nos_plots),numWorkers), type = "PSOCK")
	    # load R libraries in cluster
	    clusterExport(cl,c("load_r_libraries","rmse","gsi_controlling"))
	    clusterEvalQ(cl, load_r_libraries())
	    dummy=parLapply(cl,nos_plots,fun=uncertainty_figures,PROJECT=PROJECT,states_all=states_all,drivers=drivers,parameters=parameters,n=n,plotconfidence=plotconfidence)
	    stopCluster(cl)
	} else {
	    # or use serial
	    dummy=lapply(nos_plots,FUN=uncertainty_figures,PROJECT=PROJECT,states_all=states_all,drivers=drivers,parameters=parameters,n=n,plotconfidence=plotconfidence)
	} # parallel option

	# tidy before leaving
	gc(reset=TRUE, verbose=FALSE)

}
## Use byte compile
generate_uncertainty_figures<-cmpfun(generate_uncertainty_figures)
