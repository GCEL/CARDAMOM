
###
## Function to describe the PROJECT requirements
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

cardamom_project_setup <- function (paths,PROJECT) {

  # local paths to be set
  typepath=paste(paths$cardamom_output,"/",PROJECT$type,"/",sep="")
  localpath=paste(paths$cardamom_output,"/",PROJECT$type,"/",PROJECT$name,"/",sep="")
  datapath=paste(localpath,"DATA/",sep="")
  resultspath=paste(localpath,"RESULTS/",sep="")
  results_processedpath=paste(localpath,"RESULTS_PROCESSED/",sep="")
  figpath=paste(localpath,"FIGURES/",sep="")
  exepath=paste(localpath,"EXECUTABLE/",sep="")

  # some useful variables
  modelname=PROJECT$model$name
  parameter_type=PROJECT$parameter_type
  project_src=PROJECT$source
  project_type=PROJECT$type

  # create local paths if they do no exist already
  if (dir.exists(typepath) == FALSE ){system(paste("mkdir ",typepath,sep=""))}
  if (dir.exists(localpath) == FALSE ){system(paste("mkdir ",localpath,sep=""))}
  if (dir.exists(datapath) == FALSE ){system(paste("mkdir ",datapath,sep=""))}
  if (dir.exists(resultspath) == FALSE ){system(paste("mkdir ",resultspath,sep=""))}
  if (dir.exists(results_processedpath) == FALSE ){system(paste("mkdir ",results_processedpath,sep=""))}
  if (dir.exists(figpath) == FALSE ){system(paste("mkdir ",figpath,sep=""))}
  if (dir.exists(exepath) == FALSE ){system(paste("mkdir ",exepath,sep=""))}

  # number of chains desired?
  failed=TRUE
  # Check whether we have one defined already
  if (exists("request_nos_chains")) {
      nochains = request_nos_chains
      if (nochains > 2 & nochains < 13) {failed=FALSE} else {failed=TRUE}
  }
  # If not or not defined correctly ask the user
  while (failed) {
     nochains=as.integer(readline("How many chains?"))
     if (nochains > 1 & nochains < 11) {failed=FALSE} else {failed=TRUE}
  }

  # number of iterations
  failed=TRUE
  # Check whether we have one defined already
  if (exists("request_nos_samples")) {
      nsamples = request_nos_samples
      if (nsamples >= 1e4 & nsamples <= 1000e6) {failed=FALSE} else {failed=TRUE}
  }
  # If not or not defined correctly ask the user
  while(failed) {
    nsamples=as.integer(readline("Confirm how many parameters to sample (max = 500e6)?"))
    if (nsamples >= 1e4 & nsamples <= 1000e6) {failed=FALSE} else {failed=TRUE}
  }

  # number of parameters to sub-sample per chain
  failed=TRUE
  # Check whether we have one defined already
  if (exists("request_nos_subsamples")) {
      nsubsamples = request_nos_subsamples
      if (nsubsamples > 1e2 & nsubsamples < 1e6 & nsubsamples < nsamples) {failed=FALSE} else {failed=TRUE}
  }
  # If not or not defined correctly ask the user
  while(failed) {
    nsubsamples=as.integer(readline("How many parameters to keep per chain (recommend = 1000)?"))
    if (nsubsamples > 1e2 & nsubsamples < 1e6 & nsubsamples < nsamples) {failed=FALSE} else {failed=TRUE}
  }

  # how much of sample chain to be kept
  latter_sample_frac=0.5
  # approximate chain runtime estimate in hours (approx 1 our for 1 million iterations)
  cre = nsamples/1.0e6

  # PROJECT discription
  description = ""#readline("Any other comments?")

  # calculate the sample rate
  samplerate = nsamples/nsubsamples

  # creation date
  date = Sys.time()
  # remove the spaces and other characters
  date = gsub("-", "",date)
  date = gsub(" ", "_",date)
  date = gsub(":", "_",date)

  # define executable name
  exe = paste(PROJECT$name,".exe",sep="")

  # Are we using the remote server?
  # Has this information already been provided...
  if (exists("request_use_server")) {
      use_eddie = request_use_server
      # If the value is not valid then ask the user
      if (use_eddie != TRUE & use_eddie != FALSE) {
          use_eddie = readline("Will you run this PROJECT on remote server (TRUE/FALSE)")
      }
  } else {
      # ... or ask the user
      use_eddie=readline("Will you run this PROJECT on remote server (TRUE/FALSE)")
  }

  # If so then so something
  if (use_eddie == TRUE) {
      use_eddie = TRUE

      # Do we already know how long to run for on server?
      if (exists("request_runtime")) {
          chain_runtime = request_runtime
      } else {
          # ask the user how long they want to set the simulation to run for
          chain_runtime = readline(paste("How much run time per chain do you want to request (whole hours)? (NOTE: Given ",nsamples," required parameter vectors per chain, approximately ",cre," hours needed."))
      }
      # If response is not acceptable try again
      if (as.numeric(chain_runtime) > 48 | as.numeric(chain_runtime) < 1) {
          chain_runtime=readline(paste("Maximum number of hours to be submitted to Eddie is 48 hours, please re-select the number of hours"))
      }
      # do we want an email to notify you of eddie works
      email=""#readline("Enter your email address for remote server notification (if you want)")
  } else {
      use_eddie = FALSE
      chain_runtime = 48
  }

  # do I compile on eddie?
  if (use_eddie) {
    # are we using eddie or not
    eddiepath=paste(paths$cardamom_ecdf,project_type,PROJECT$name,sep="/")
    ecdf_source=paste(paths$cardamom_ecdf,"LIBRARY/",sep="/")
    # declare the eddie specific paths
    edatapath=paste(eddiepath,"DATA/",sep="/")
    eresultspath=paste(eddiepath,"RESULTS/",sep="/")
    eoestreampath=paste(eddiepath,"OUTPUT_ERROR_STREAM/",sep="/")
    eexepath=paste(eddiepath,"EXECUTABLES/",sep="/")

    # check current host address
    #home_computer=Sys.info()["nodename"]

    # generate cardamom submit scripts
    #	  generate_eddie_submit_script(paths)

    # create directories on eddie and copy some important shell scripts
    commands=c(paste("mkdir ",paths$cardamom_ecdf,sep="")
              ,paste("mkdir ",ecdf_source,sep="")
              ,paste("mkdir ",paths$cardamom_ecdf,"/",project_type,sep="")
              ,paste("mkdir ",paths$cardamom_ecdf,"/",project_type,"/",PROJECT$name,sep="")
              ,paste("mkdir ",edatapath,sep="")
              ,paste("mkdir ",eresultspath,sep="")
              ,paste("mkdir ",eoestreampath,sep="")
              ,paste("mkdir ",eexepath,sep="")
              ,paste("scp ",username,"@",home_computer,":",paths$cardamom,"/R_functions/CARDAMOM_ECDF_SUBMIT_BUNDLES.sh ",eexepath,sep="")
              ,paste("chmod +x ",eexepath,"/CARDAMOM_ECDF_SUBMIT_BUNDLES.sh",sep=""))

    # Have we been given this information already?
    if (exists("request_compile_server")) {
        comline = request_compile_server
    } else {
        comline = readline("Copy and compile any source code updates to Eddie (TRUE/FALSE)?")
    }
    if (comline) {
        print("Backup source code currently on eddie first")
        print("Then copying source code to eddie")
        print("Finally compile source code on eddie")
        if (project_src == "C") {
            commands = append(commands,c(paste("mv ",ecdf_source,"CARDAMOM_C ",ecdf_source,"CARDAMOM_C_BKP",sep="")
                             ,paste("scp -r ",username,"@",home_computer,":",paths$cardamom,"LIBRARY/CARDAMOM_C ",ecdf_source,sep="")
                             ,paste("gcc ",ecdf_source,"CARDAMOM_C/projects/DALEC_CDEA_TEMPLATE/DALEC_CDEA_TEMPLATE.c -o ",ecdf_source,
                                    "CARDAMOM_C/projects/DALEC_CDEA_TEMPLATE/a.out -lm",sep="")
                             ,paste("cp ",ecdf_source,"CARDAMOM_C/projects/DALEC_CDEA_TEMPLATE/a.out ",eexepath,"/",exe,sep="")))
        } else if (project_src == "Fortran") {
            # compiler options
            compiler_options=""#"-xhost -ipo -no-ftz"
            if (timing) {compiler_options=paste(compiler_options," -pg",sep="")}
            if (debug) {compiler_options=paste(compiler_options," -debug -backtrace",sep="")}
            commands=append(commands,c(paste("rm -r ",ecdf_source,"CARDAMOM_F_BKP",sep="")
                       ,paste("mv ",ecdf_source,"CARDAMOM_F ",ecdf_source,"CARDAMOM_F_BKP",sep="")
                       ,paste("scp -r ",username,"@",home_computer,":",paths$cardamom,"LIBRARY/CARDAMOM_F ",ecdf_source,sep="")
                       ,paste("cd ",ecdf_source,"CARDAMOM_F/executable",sep="")
                       ,paste("rm cardamom.exe") # depends on working directory "executable"
                       ,paste("rm *.mod")        # depends on working directory "executable"
                       ,paste(compiler," -O2 ",compiler_options," ../misc/math_functions.f90 ../misc/oksofar.f90 ../model/",modelname,"/src/",modelname,".f90",
                              " ../general/cardamom_structures.f90 ../method/MHMCMC/MCMC_FUN/MHMCMC_STRUCTURES.f90 ../method/MHMCMC/MCMC_FUN/MHMCMC_StressTests.f90",
                              " ../model/",modelname,"/src/",modelname,"_PARS.f90 ../general/cardamom_io.f90 ../method/MHMCMC/MCMC_FUN/MHMCMC.f90",
                              " ../model/",modelname,"/likelihood/MODEL_LIKELIHOOD.f90 ../general/cardamom_main.f90 -o cardamom.exe",sep="")
                       ,paste("cp ",ecdf_source,"CARDAMOM_F/executable/cardamom.exe ",eexepath,"/",exe,sep="")))
            # If a crop model the copy the crop development files into place too
            if (modelname == "DALEC.A1.C3.H2.M1" | modelname == "DALEC.C3.M1") {
                commands=append(commands,paste("cp ",ecdf_source,"CARDAMOM_F/model/",modelname,"/src/winter_wheat_development.csv ",eexepath,"/",sep=""))
                system(paste("cp ",paths$cardamom,"LIBRARY/CARDAMOM_F/model/",modelname,"/src/winter_wheat_development.csv ",exepath,"/",sep=""))
            } #

        } else {
          stop('Source code language has not been specified')

        } # Language choice

    } # Compile on remote server

    # issue commands to eddie
    #print(commands)
    ecdf_execute(commands,PROJECT$paths$cardamom_cluster)

  }
  # then compile locally
  comline="y"
  if (comline == "y") {
    print("Finally compile source code locally")
    if (project_src == "C") {
        system(paste("gcc ",paths$cardamom,"LIBRARY/CARDAMOM_C/projects/DALEC_CDEA_TEMPLATE/DALEC_CDEA_TEMPLATE.c -o ",
                     paths$cardamom,"LIBRARY/CARDAMOM_C/projects/DALEC_CDEA_TEMPLATE/a.out -lm",sep=""))
        system(paste("cp ",paths$cardamom,"LIBRARY/CARDAMOM_C/projects/DALEC_CDEA_TEMPLATE/a.out ",exepath,"/",exe,sep=""))
    } else if (project_src == "Fortran") {

        # store current working directory so that we can leave it briefly but return later
        cwd=getwd()
        # Move to the directory containing the source ode
        setwd(paste(paths$cardamom,"LIBRARY/CARDAMOM_F/executable/",sep=""))
        # Remove the evidence of the previuos compilation
        system("rm *.mod") # depends on being in executable directory

        if (request_compile_server == FALSE | request_compile_local == TRUE) {
            # compiler options
            compiler_options=""#"-xhost -ipo -no-ftz"
            if (timing) {compiler_options=paste(compiler_options," -pg",sep="")}
            if (debug) {compiler_options=paste(compiler_options," -debug -traceback",sep="")}
            # if executables are present either in source library or in project folder remove them
            #if (file.exists(paste(paths$cardamom,"LIBRARY/CARDAMOM_F/executable/cardamom.exe",sep=""))) { system(paste("rm cardamom.exe")) }
            if (file.exists(paste(exepath,"/",exe,sep=""))) {system(paste("rm ",exepath,"/",exe,sep=""))}
            # issue compile commands
            system(paste(compiler," -O2 ",compiler_options," ../misc/math_functions.f90 ../misc/oksofar.f90 ../model/",modelname,"/src/",modelname,".f90",
                         " ../general/cardamom_structures.f90 ../method/MHMCMC/MCMC_FUN/MHMCMC_STRUCTURES.f90 ../method/MHMCMC/MCMC_FUN/MHMCMC_StressTests.f90",
                         " ../model/",modelname,"/src/",modelname,"_PARS.f90 ../general/cardamom_io.f90 ../method/MHMCMC/MCMC_FUN/MHMCMC.f90",
                         " ../model/",modelname,"/likelihood/MODEL_LIKELIHOOD.f90 ../general/cardamom_main.f90 -o cardamom.exe",sep=""))
#            print(paste(compiler," -O2 ",compiler_options," ../misc/math_functions.f90 ../misc/oksofar.f90 ../model/",modelname,"/src/",modelname,".f90",
#                         " ../general/cardamom_structures.f90 ../method/MHMCMC/MCMC_FUN/MHMCMC_STRUCTURES.f90 ../method/MHMCMC/MCMC_FUN/MHMCMC_StressTests.f90",
#                         " ../model/",modelname,"/src/",modelname,"_PARS.f90 ../general/cardamom_io.f90 ../method/MHMCMC/MCMC_FUN/MHMCMC.f90",
#                         " ../model/",modelname,"/likelihood/MODEL_LIKELIHOOD.f90 ../general/cardamom_main.f90 -o cardamom.exe",sep=""))
            system(paste("cp ",paths$cardamom,"LIBRARY/CARDAMOM_F/executable/cardamom.exe ",exepath,"/",exe,sep=""))
        } # compile executable on local machine too?

        # Generate the shared library needed later by R
        #print(paste("gfortran -fcheck=all -O2 -shared ../model/",modelname,"/src/",modelname,".f90 ",
        #             "../model/",modelname,"/src/",modelname,"_R_interface.f90 ","-o dalec.so -fPIC",sep=""))
        system(paste("gfortran -fcheck=all -O2 -shared ../model/",modelname,"/src/",modelname,".f90 ",
                     "../model/",modelname,"/src/",modelname,"_R_interface.f90 ","-o dalec.so -fPIC",sep=""))
        system(paste("cp ",paths$cardamom,"LIBRARY/CARDAMOM_F/executable/dalec.so ",exepath,"/dalec.so",sep=""))
        # Copy crop development file into position
        if (modelname == "DALEC.A1.C3.H2.M1." | modelname == "DALEC.C3.M1.") {
             system(paste("cp ",paths$cardamom,"LIBRARY/CARDAMOM_F/model/",modelname,"/src/winter_wheat_development.csv ",exepath,"/",sep=""))
        } #
        # return to original working directory
        setwd(cwd)

    } else {

        stop('Source code language has not been specified')

    } # Compile locally

  } # copy and compile to eddie

  # prepare output
  PROJECT$ecdf=use_eddie
  PROJECT$localpath=localpath
  PROJECT$datapath=datapath
  PROJECT$resultspath=resultspath
  PROJECT$results_processedpath=results_processedpath
  PROJECT$figpath=figpath
  PROJECT$exepath=exepath
  PROJECT$nochains=nochains
  PROJECT$nsamples=nsamples
  PROJECT$latter_sample_frac=latter_sample_frac
  PROJECT$nsubsamples=nsubsamples
  PROJECT$chain_runtime=chain_runtime
  PROJECT$description=description
  PROJECT$samplerate=samplerate
  PROJECT$date=date
  PROJECT$exe=exe
  # eddie specific information
  if (use_eddie) {
    PROJECT$eddiepath=eddiepath
    PROJECT$edatapath=edatapath
    PROJECT$eresultspath=eresultspath
    PROJECT$eoestreampath=eoestreampath
    PROJECT$eexepath=eexepath
    PROJECT$email=email
  }
  return(PROJECT)
} # function end cardamom_project_setup

## Use byte compile
cardamom_project_setup<-cmpfun(cardamom_project_setup)
