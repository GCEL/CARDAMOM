
###
## Function to submit processes to eddie
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE) & J. F. Exbrayat (UoE).

submit_processes_to_local_machine<-function (PROJECT_in) {

    print('PREPARING TO SUBMIT MCMC TO LOCAL MACHINE')
    print('This function should be valid for all CARDAMOM compatible DALEC MCMC functions')

    # do we want to remove any previous output files?
    delete_old = readline("Delete any previous output files for this project name?(y/n)")
    if (delete_old == "y") {
        system(paste("rm ",PROJECT_in$resultspath,"/*",sep=""))
    }
    # do we want to run tasks in the background?
    background = readline("Run jobs in background?(y/n)")

    # Determine number of sites to submit concurrently
    if (use_parallel) {
        # In parallel use number of cores into indicate the override
        concurrent_sites = ceiling(numWorkers / PROJECT_in$nochains)
    } else {
        # In serial do not over right
        concurrent_sites = 1
    }

    # Assume MCMC withh use pre-mcmc where lielihoods are normalised by sample size
    pre_mcmc = 1   
    # Combined the default set of parameter proposals with the extended run number
    # if this is an extened run
    nsamples = as.integer(PROJECT_in$nsamples)
    if (request_extended_mcmc) {
        nsamples = nsamples + request_nos_extended_samples
        # Assume pre-mcmc is not to be used if this is an extended run
        pre_mcmc = 0
    }

    # Check presence of PROJECT_in$cost_function_scaling
    if (exists(x = "cost_function_scaling", where = PROJECT_in) == FALSE) {
        # If not, assume default cost function
        PROJECT_in$cost_function_scaling = 0
    }
    
    # begin submitting the different tasks
    cwd = getwd()
    setwd(PROJECT_in$exepath)
    for (n in seq(1, PROJECT_in$nosites)) {
         # Override background request?
         bg_override = TRUE ; if (n%%concurrent_sites == 0) {bg_override = FALSE}
         for (c in seq(1, PROJECT_in$nochains)) {
              # Create the input / output file names for the current job
              infile=paste(PROJECT_in$datapath,PROJECT_in$name,"_",PROJECT_in$sites[n],".bin",sep="")
              output=paste(PROJECT_in$resultspath,PROJECT_in$name,"_",PROJECT_in$sites[n],"_",c,"_",sep="")
              # Depending on whether or not a background run submit the job
              if (background == "y" | c != PROJECT_in$nochains | bg_override) {
                  system(paste(PROJECT_in$exepath,PROJECT_in$exe," ",
                               infile," ",
                               output," ",
                               as.integer(nsamples),
                               " 0 ",
                               as.integer(PROJECT_in$samplerate)," ",
                               as.integer(pre_mcmc)," ",
                               as.integer(PROJECT_in$cost_function_scaling)," & ",sep=""))
              } else {
                  system(paste(PROJECT_in$exepath,PROJECT_in$exe," ",
                               infile," ",
                               output," ",
                               as.integer(nsamples),
                               " 0 ",
                               as.integer(PROJECT_in$samplerate)," ",
                               as.integer(pre_mcmc)," ",
                               as.integer(PROJECT_in$cost_function_scaling),sep=""))
              }
              # To ensure that each chain is submitted at a unique time we want to delay the code - this impacts the seed value used in the random number generator
              Sys.sleep(1) # wait for 1 seconds
         } # chain loop
    } # site loop
    setwd(cwd) ; rm(cwd)
    print("Command issued to local machine")

} # end of function submit_processes_to_local_machine

## Use byte compile
submit_processes_to_local_machine<-cmpfun(submit_processes_to_local_machine)
