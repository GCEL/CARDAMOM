
###
## Function to submit processes to eddie
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE) & J. F. Exbrayat (UoE).

submit_processes_to_cluster<-function (PROJECT_in) {

    print('PREPARING TO SUBMIT MCMC TO ECDF')
    print('This function should be valid for all CARDAMOM compatible DALEC MCMC functions')

    # Assume MCMC withh use pre-mcmc where likelihoods are normalised by sample size
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

    # create the new file name in the correct location
    outfile=paste(PROJECT_in$exepath,"CARDAMOM_ECDF_EXECUTABLES_LIST.txt",sep="")
    # begin writing out the file contents
    # construct the file now
    first_pass=TRUE
    for (c in seq(1, PROJECT_in$nochains)) {
         for (n in seq(1, PROJECT_in$nosites)) {
              infile=paste(PROJECT_in$edatapath,PROJECT_in$name,"_",PROJECT_in$sites[n],".bin",sep="")
              output=paste(PROJECT_in$eresultspath,PROJECT_in$name,"_",PROJECT_in$sites[n],"_",c,"_",sep="")
              if (first_pass) {
                  write(paste(PROJECT_in$eexepath,PROJECT_in$exe," ",
                              infile," ",
                              output," ",
                              as.integer(nsamples),
                              " 0 ",
                              as.integer(PROJECT_in$samplerate)," ",
                              as.integer(pre_mcmc)," ",
                              as.integer(PROJECT_in$cost_function_scaling),sep=""),sep=" ", ncolumn=1,file=outfile,append="F")
                  first_pass=FALSE
              } else {
                  write(paste(PROJECT_in$eexepath,PROJECT_in$exe," ",
                              infile," ",
                              output," ",
                              as.integer(nsamples),
                              " 0 ",
                              as.integer(PROJECT_in$samplerate)," ",
                              as.integer(pre_mcmc)," ",
                              as.integer(PROJECT_in$cost_function_scaling),sep=""),sep=" ", ncolumn=1,file=outfile,append="T")
              }
         } # chain no
    } # nosite

    # do we want to remove any previous output files?
    delete_old=readline("Delete any previous output files for this project name?(y/n)")
    while(delete_old != "y" & delete_old != "n") {
          # do we want to remove the existing files
          delete_old=readline("Delete any previous output files for this project name?(y/n)")
    }
    # write commands to clear previous runs and copy new commands
    commands=c(paste("rm ",PROJECT_in$eresultspath,"*",sep="")
              ,paste("rm ",PROJECT_in$eoestreampath,"*",sep="")
              ,paste("scp ",username,"@",home_computer,":",outfile," ",PROJECT_in$eexepath,sep="")
              ,paste("cd ",PROJECT_in$eexepath, sep=""))
    if (delete_old == "n") {commands=commands[-1]}

    ## default information for eddie submission
    # number of bundles allowed
    nbundles=9000

    # number of tasks required
    ntasks=PROJECT_in$nochains*PROJECT_in$nosites
    # number of bundles needed for tasks (assuming max 5000 bundle limit)
    bundlesize=ceiling(ntasks/nbundles)
    # number of tasks per bundle
    ntaskbundles=ceiling((PROJECT_in$nochains*PROJECT_in$nosites)/bundlesize)
    # make the size bundle specific to adjust for hangers on
    ntaskbundles=rep(ntaskbundles, times=bundlesize)
    # place any hangers on into the last bundle
    ntaskbundles[bundlesize]=ntaskbundles[bundlesize]+(ntasks%%bundlesize)
    # task run time
    runtimestr=paste(" -l h_rt=",as.numeric(PROJECT_in$chain_runtime),":00:00 ",sep="")
    # task memory needs
    memorystr=paste(" -l h_vmem=1G",sep="")

    # eddie email link
    if (grepl("@",PROJECT_in$email)) {
        emailstr=paste(" -m beas -M ",PROJECT_in$email,sep="")
    } else {
        emailstr=""
    }
    # eddie output and error txt files dump
    oestream=paste(" -o ",PROJECT_in$eoestreampath," -e ",PROJECT_in$eoestreampath,sep="")

    print(paste('Number of tasks to be submitted = ',ntasks,sep=""))
    print(paste('Maximum number of tasks allowed = ',nbundles,sep=""))
    print(paste('Tasks will be bundled in groups of  ~ ',mean(ntaskbundles),sep=""))
    print(paste('Number of bundles to be submitted = ',bundlesize,sep=""))

    # bundle submission options
    bundle_choice=readline("Submit all (a), or submit specific bundles (b)?")
    if (bundle_choice != "a" & bundle_choice != "b") {bundle_choice=readline("Did not correctly select all (a), or specific bundles (b)?")}

    # now further information
    if (bundle_choice == "a") {
        start_point=1 ; end_point=bundlesize
    } else if (bundle_choice == "b") {
        start_point=as.numeric(readline("Please provide the number of the first bundle to submit"))
        end_point=as.numeric(readline("Please provide the number of the last bundle to submit"))
        if (start_point < 0 | start_point > bundlesize | end_point < 0 | end_point > bundlesize | end_point < start_point) {
            start_point=as.numeric(readline("Invalid start and / end point bundle provided try again, FIRST bundle to submit now please"))
            end_point=as.numeric(readline("now the last bundle to submit"))
        }
    }

    # now carry out the submission of the different arrays
    for (current_bundle in seq(start_point,end_point)) {
         # adjust to ntask number
         if (current_bundle == 1) {
             bundle_start=1
             bundle_end=as.numeric(ntaskbundles[current_bundle])
         } else {
            bundle_start=sum(as.numeric(ntaskbundles[1:(as.numeric(current_bundle)-1)]))+1
            bundle_end=sum(as.numeric(ntaskbundles[1:as.numeric(current_bundle)]))
         } # current = 1 condition
         # now submit them all
         jobnamestr=paste(" -N ",PROJECT_in$name,"_bundle_",bundle_start,"_",bundle_end,sep="")
         commands=append(commands,paste("qsub -t ",bundle_start,"-",bundle_end,jobnamestr,emailstr,memorystr,runtimestr,oestream," ",PROJECT_in$eexepath,"CARDAMOM_ECDF_SUBMIT_BUNDLES.sh ",PROJECT_in$eexepath,sep=""))
         #print(commands)
    } # bundle looping

    # issue commands to eddie
    ecdf_execute(commands,PROJECT_in$paths$cardamom_cluster)
    print("Command issued to cluster")

} # end of function
## Use byte compile
submit_processes_to_cluster<-cmpfun(submit_processes_to_cluster)
