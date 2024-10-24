
###
## Function to submit processes to eddie
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE) & J. F. Exbrayat (UoE).

submit_processes_to_local_slurm_machine<-function (PROJECT_in) {

    print('PREPARING TO SUBMIT MCMC TO LOCAL CLUSTER MACHINE (SLURM scheduler assumed)')
    print('This function should be valid for all CARDAMOM compatible DALEC MCMC functions')

    ## Some housekeeping

    # Remove any previous output files?
    delete_old = readline("Delete any previous output files for this project name?(y/n)")
    if (delete_old == "y") {
        system(paste("rm ",PROJECT_in$resultspath,"/*",sep=""))
    }

    # CARDAMOM typically uses a multi-phase MCMC process. First, an EDC searching phase to
    # to find parameters with an EDC compliant starting point. Second, a pre-mcmc during which
    # the likelihood scores for each data stream are normalised by their sample size. Third, the 
    # main analysis during which the likelihood scores are weighted based on the 'cost_function_scaling'
    # However, if this is an extended run, i.e. going beyond the parameter proposals originally requested
    # the pre-mcmc must be turned off to maintain consistency in the likelihood scores being assessed.
    if (request_extended_mcmc) {
        # In an extended run the number of proposals (samples) is added to by the requested extension
        nsamples = as.integer(PROJECT_in$nsamples) + request_nos_extended_samples
        # Assume pre-mcmc is not to be used if this is an extended run
        pre_mcmc = 0
    } else {
        # In a normal run the number of proposals is passed into the local variable unchanged
        nsamples = as.integer(PROJECT_in$nsamples)
        # Asse the pre-mcmc is used as default
        pre_mcmc = 1
    }

    # Check presence of PROJECT_in$cost_function_scaling
    if (exists(x = "cost_function_scaling", where = PROJECT_in) == FALSE) {
        # If not, assume default cost function
        PROJECT_in$cost_function_scaling = 0
    }

    ## Create the two files needed, one which contains the list of jobs to be ran
    ## and a second which write the correct slurm job submission script

    # create the new file name in the correct location
    outfile = paste(PROJECT_in$exepath,"CARDAMOM_ECDF_EXECUTABLES_LIST.txt",sep="")
    # begin writing out the file contents
    # construct the file now
    first_pass=TRUE
    for (c in seq(1, PROJECT_in$nochains)) {
         for (n in seq(1, PROJECT_in$nosites)) {
              infile = paste(PROJECT_in$datapath,PROJECT_in$name,"_",PROJECT_in$sites[n],".bin",sep="")
              output = paste(PROJECT_in$resultspath,PROJECT_in$name,"_",PROJECT_in$sites[n],"_",c,"_",sep="")
              if (first_pass) {
                  write(paste(PROJECT_in$exepath,PROJECT_in$exe," ",
                              infile," ",
                              output," ",
                              as.integer(nsamples),
                              " 0 ",
                              as.integer(PROJECT_in$samplerate)," ",
                              as.integer(pre_mcmc)," ",
                              as.integer(PROJECT_in$cost_function_scaling),sep=""),sep=" ", ncolumn=1,file=outfile,append="F")
                  first_pass=FALSE
              } else {
                  write(paste(PROJECT_in$exepath,PROJECT_in$exe," ",
                              infile," ",
                              output," ",
                              as.integer(nsamples),
                              " 0 ",
                              as.integer(PROJECT_in$samplerate)," ",
                              as.integer(pre_mcmc)," ",
                              as.integer(PROJECT_in$cost_function_scaling),sep=""),sep=" ", ncolumn=1,file=outfile,append="T")
              } # first pass or not
         } # chain no
    } # nosite

    ## default information for cluster submission
    # number of tasks per array
    max_nbundle_size = 1000
    # If the slurm account does not already exist in memory assume this one
    if (exists("slurm_account") == FALSE) { slurm_account = "geos_research" }

    # number of tasks required
    ntasks=PROJECT_in$nochains*PROJECT_in$nosites
    # number of bundles needed for tasks 
    nbundle=ceiling(ntasks/max_nbundle_size)
    # number of tasks per bundle
    ntaskbundles=ceiling(ntasks/nbundle)
    # make the size bundle specific to adjust for hangers on
    ntaskbundles=rep(ntaskbundles, times=nbundle)
    # place any hangers on into the last bundle
    ntaskbundles[nbundle]=ntaskbundles[nbundle]+(ntasks%%nbundle)

    print(paste('Number of tasks to be submitted = ',ntasks,sep=""))
    print(paste('Maximum number of tasks allowed = ',max_nbundle_size,sep=""))
    print(paste('Tasks will be bundled in groups of  ~',mean(ntaskbundles),sep=""))
    print(paste('Number of bundles to be submitted = ',nbundle,sep=""))

    # Loop through each bundle and submit to the local slurm cluster
    for (b in seq(1, nbundle)) {
         # Determine the start and end points of the bundles to be submitted
         if (b == 1) { bundle_start = 1 } else { bundle_start = sum(ntaskbundles[1:(b-1)] + 1) }
         bundle_end = sum(ntaskbundles[1:b])

         # Create the shell script for submitting the job to slurm on the local cluster
         slurm_file = paste(PROJECT_in$exepath,"/slurm_submission.sh",sep="")
         col_sep = "" ; nos_cols = 20
         write(    c("#!/bin/bash"), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = FALSE)
         write(    c(" "), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c("# Slurm directives"), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c(paste('#SBATCH --account=',slurm_account,sep="")), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c(paste("#SBATCH --job-name=Bundle_",b,sep="")), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c("#SBATCH --ntasks=1"), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c("#SBATCH --cpus-per-task=1"), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c("#SBATCH --mem=1G "), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c(paste('#SBATCH --output="',PROJECT_in$oestreampath,'/slurm-%A_%a.out"',sep="")), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c(paste("#SBATCH --time=",as.numeric(PROJECT_in$chain_runtime),":00:00",sep="")), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c(paste("#SBATCH --array=",bundle_start,"-",bundle_end,sep="")), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c(" "), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c("# THIS SCRIPT MUST BE ACCOMPANIED BY CARDAMOM_ECDF_EXECUTABLES_LIST.txt IN THE SAME DIRECTORY"), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c("# arguments are start and end lines!"), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c(" "), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c(paste("task=$( cat $1CARDAMOM_ECDF_EXECUTABLES_LIST.txt | sed $SLURM_ARRAY_TASK_ID\\!d )",sep="")), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c("command ${task}"), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)

         # Record directory to change back in a moment
         cwd = getwd()
         # Set working directory to the location of the executable we want to run
         setwd(PROJECT_in$exepath)
         # Submit jobs to the local slurm cluster
         system(paste("sbatch ",slurm_file,sep=""))
         # Return back to normal working directory
         setwd(cwd) ; rm(cwd)

    } # loop for submission of batches
    
    # Inform the user
    print("Command issued to local slurm machine")

} # end of function submit_processes_to_local_slurm_machine

## Use byte compile
submit_processes_to_local_slurm_machine<-cmpfun(submit_processes_to_local_slurm_machine)


###
## Function to submit processes to eddie
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE) & J. F. Exbrayat (UoE).

submit_R_run_each_site_to_local_slurm_machine<-function(PROJECT_in,repair,job_ID) {

    # Update the user
    print('PREPARING TO SUBMIT STAGE 3 REPROCESSING TO LOCAL CLUSTER MACHINE (SLURM scheduler assumed)')
    print('This function should be valid for all CARDAMOM compatible DALEC models functions')

    ## Create the correct slurm job submission script for run_each_site()

    ## default information for cluster submission
    # number of tasks per array
    max_nbundle_size = 5000
    # If the slurm account does not already exist in memory assume this one
    if (exists("slurm_account") == FALSE) { slurm_account = "geos_research" }

    # number of tasks required
    ntasks = PROJECT_in$nosites
    # number of bundles needed for tasks 
    nbundle = ceiling(ntasks/(max_nbundle_size-1))
    # number of tasks per bundle
    ntaskbundles = ceiling(ntasks/nbundle)
    # make the size bundle specific to adjust for hangers on
    ntaskbundles = rep(ntaskbundles, times=nbundle)
    # place any hangers on into the last bundle
    ntaskbundles[nbundle] = ntaskbundles[nbundle]+(ntasks%%nbundle)

    # Update the user
    print(paste('Number of tasks to be submitted = ',ntasks,sep=""))
    print(paste('Maximum number of tasks allowed = ',max_nbundle_size,sep=""))
    print(paste('Tasks will be bundled in groups of  ~',mean(ntaskbundles),sep=""))
    print(paste('Number of bundles to be submitted = ',nbundle,sep=""))
# MAKE SCRIPT NOW RUN MULTIPLE SITES PER TASKS
    # Determine the file path for the PROJECT infofile.RData
    project_path = paste(PROJECT_in$localpath,"/infofile.RData",sep="")

    # Loop through each bundle and submit to the local slurm cluster
    for (b in seq(1, nbundle)) {
         # Locally extract the size of the current bundle
         bundle_end = ntaskbundles[b]
         # and determine the offset required to get to the right site number.
         # This is done to avoid the SLURM max array counter size limit.
         if (b == 1) { bundle_offset = 0 } else { bundle_offset = sum(ntaskbundles[1:(b-1)]) }

         # Create the shell script for submitting the job to slurm on the local cluster
         slurm_file = paste(PROJECT_in$exepath,"/slurm_stage_3_submission.sh",sep="")
         col_sep = "" ; nos_cols = 20
         write(    c("#!/bin/bash"), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = FALSE)
         write(    c(" "), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c("# Slurm directives"), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c(paste('#SBATCH --account=',slurm_account,sep="")), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c(paste("#SBATCH --job-name=",job_ID,"_Bundle_",b,sep="")), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c("#SBATCH --ntasks=1"), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c("#SBATCH --cpus-per-task=1"), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c("#SBATCH --mem=1G "), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c(paste('#SBATCH --output="',PROJECT_in$oestreampath,'/slurm-%A_%a.out"',sep="")), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c(paste("#SBATCH --time=00:30:00",sep="")), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c(paste("#SBATCH --array=[1-",bundle_end,"]",sep="")), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c(" "), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c("sitenum=$((SLURM_ARRAY_TASK_ID+",bundle_offset,"))"), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)
         write(    c(paste("R --no-save < ",PROJECT_in$paths$cardamom,"/R_functions/run_each_site_slurm.r --args ",PROJECT_in$localpath,"/infofile.RData ",repair," $sitenum",sep="")), file = slurm_file, ncolumns = nos_cols, sep=col_sep, append = TRUE)

         # Record directory to change back in a moment
         cwd = getwd()
         # Set working directory to the location of the executable we want to run
         setwd(PROJECT_in$exepath)
         # Submit jobs to the local slurm cluster
         print(paste("...bundle ",b," of ",nbundle,sep=""))
         system(paste("sbatch ",slurm_file,sep=""))
         # Return back to normal working directory
         setwd(cwd) ; rm(cwd)

         ## Hack to work around maximum number of array jobs being submitted
         # Check whether the slurm scheduler has finished all jobs
         Sys.sleep(1) ; ongoing = TRUE
         while(ongoing) {
            # Query ongoing jobs, assumes only your user name is returned
            system(paste('squeue -u ',username,' --format="%20j %2t" > q',sep="")) ; q = read.table("q", header=TRUE)
            # If this job ID is still found in the queue we probably want to continue waiting...
            q_filter = grepl(job_ID,q[,1])
            if (length(which(q_filter)) > 0) {
                # ...but we want to check whether the tasks listed are anything other than 
                # completing (CG)
                if (length(which(grepl("CG",q[q_filter,2]))) == dim(q[q_filter,])[1]) {
                    # Just completing, therefore we will submit the next batch in the meantime
                    file.remove("q") ; ongoing = FALSE 
                } else {
                    # Otherwise, we will wait 60 seconds and check again
                    file.remove("q") ; print("...waiting") ; Sys.sleep(60)                
                }
            } else {
                # All tasks with this ID have now completely left the cluster
                # we should move on
                file.remove("q") ; ongoing = FALSE 
            } # check for active jobs
         } # while loop

    } # loop for submission of batches

    # Even after all tasks have been submitted we want to 
    # make sure that absolutely everything has passed through 
    # the completing phase
    ongoing = TRUE
    while(ongoing) {
       # Query ongoing jobs, assumes only your user name is returned
       system(paste('squeue -u ',username,' --format="%20j %2t" > q',sep="")) ; q = read.table("q", header=TRUE)
       # If anything is still on the cluster we should continue to wait
       if (length(which(grepl(job_ID,q[,1]))) > 0) {
           # Otherwise, we will wait 60 seconds and check again
           file.remove("q") ; print("...waiting") ; Sys.sleep(60)                
       } else {
           # All tasks with this ID have now completely left the cluster
           # we should move on
           file.remove("q") ; ongoing = FALSE
       } # queue check for active jobs
    } # while loop

    # Inform the user
    print("Command issued to local slurm machine")

    # Return back run_mcmc_results()
    return(0)

} # end of function submit_R_run_each_site_to_local_slurm_machine

## Use byte compile
submit_R_run_each_site_to_local_slurm_machine<-cmpfun(submit_R_run_each_site_to_local_slurm_machine)
