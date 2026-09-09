#########################################################################################
# CARbon DAta MOdel fraMework (CARDAMOM) and DALEC terrestrial ecosystem model suite
# CARDAMOM is a Bayesian model-data fusion software framework. CARDAMOM is used to 
# assimilate observations and ecological theory to retrieve parameters for the 
# DALEC suite of intermediate complexity terrestrial ecosystem models. DALEC can be
# used as a fully integrated component of CARDAMOM or independently. 
# Copyright (C) 2024  University of Edinburgh,
#                     Mathew Williams (mat.williams@ed.ac.uk), 
#                     T. Luke Smallman (t.l.smallman@ed.ac.uk)
# UoE = University of Edinburgh

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# ########## File specific description ##########
# Functions within are used to determine the commands needed to create a CARDAMOM project
# and compile CARDAMOM and DALEC codes.
# This function is based on an original Matlab function development by A. A. Bloom 
# (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).
#
#########################################################################################

cardamom_project_setup <- function (paths,PROJECT) {

  # local paths to be set
  typepath=paste(paths$cardamom_output,"/",PROJECT$type,"/",sep="")
  localpath=paste(paths$cardamom_output,"/",PROJECT$type,"/",PROJECT$name,"/",sep="")
  datapath=paste(localpath,"DATA/",sep="")
  oestreampath=paste(localpath,"OUTPUT_ERROR_STREAM/",sep="/")  
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
  if (dir.exists(oestreampath) == FALSE ){system(paste("mkdir ",oestreampath,sep=""))}
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

      ## Use of remote server (Eddie) has been chosen

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

      ## Use of local machine has been selected

      use_eddie = FALSE
      chain_runtime = 48
      if (request_use_local_slurm) {
          chain_runtime = request_runtime
      }
  }

  # do I compile on eddie?
  if (use_eddie) {
      # are we using eddie or not
      eddiepath = paste(paths$cardamom_ecdf,project_type,PROJECT$name,sep="/")
      ecdf_source = paste(paths$cardamom_ecdf,"LIBRARY/",sep="/")
      # declare the eddie specific paths
      edatapath = paste(eddiepath,"DATA/",sep="/")
      eresultspath = paste(eddiepath,"RESULTS/",sep="/")
      eoestreampath = paste(eddiepath,"OUTPUT_ERROR_STREAM/",sep="/")
      eexepath = paste(eddiepath,"EXECUTABLES/",sep="/")

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
                #,paste("scp -i ",sshpass_key_home," ",username,"@",home_computer,":",paths$cardamom,"/R_functions/CARDAMOM_ECDF_SUBMIT_BUNDLES.sh ",eexepath,sep="")
                ,paste("scp ",username,"@",home_computer,":",paths$cardamom,"/R_functions/CARDAMOM_ECDF_SUBMIT_BUNDLES.sh ",eexepath,sep="")
                ,paste("chmod +x ",eexepath,"/CARDAMOM_ECDF_SUBMIT_BUNDLES.sh",sep=""))

      # Have we been given this information already?
      if (exists("request_compile_server")) {
          comline = request_compile_server
      } else {
          comline = readline("Copy and compile any source code updates to Eddie (TRUE/FALSE)?")
      }
      if (comline) {
          print("Compiling instructions to backup source code currently on eddie first...")
          print("...then copying source code to eddie and compile")
          if (project_src == "C") {
              commands = append(commands,c(paste("mv ",ecdf_source,"CARDAMOM_C ",ecdf_source,"CARDAMOM_C_BKP",sep="")
                               #,paste("scp -r -i ",sshpass_key_home," ",username,"@",home_computer,":",paths$cardamom,"LIBRARY/CARDAMOM_C ",ecdf_source,sep="")
                               ,paste("scp -r ",username,"@",home_computer,":",paths$cardamom,"LIBRARY/CARDAMOM_C ",ecdf_source,sep="")                               
                               ,paste("gcc ",ecdf_source,"CARDAMOM_C/projects/DALEC_CDEA_TEMPLATE/DALEC_CDEA_TEMPLATE.c -o ",ecdf_source,
                                      "CARDAMOM_C/projects/DALEC_CDEA_TEMPLATE/a.out -lm",sep="")
                               ,paste("cp ",ecdf_source,"CARDAMOM_C/projects/DALEC_CDEA_TEMPLATE/a.out ",eexepath,"/",exe,sep="")))
          } else if (project_src == "Fortran") {
              # Map the debug / timing switches onto the CMake build type
              build_type = "RELEASE" ; if (timing | debug) {build_type = "DEBUG"}
              # Select the target and executable name
              target = "cardamom"        ; built_exe = "cardamom.exe"
              # Remote build via CMake, mirroring the local branch. The build tree
              # is LIBRARY/CARDAMOM_F/build; executables land in
              # LIBRARY/CARDAMOM_F/executable. 'module load cmake' may be required
              # depending on the cluster environment.
              commands=append(commands,c(paste("rm -r ",ecdf_source,"CARDAMOM_F_BKP",sep="")
                                        ,paste("mv ",ecdf_source,"CARDAMOM_F ",ecdf_source,"CARDAMOM_F_BKP",sep="")
                                        #,paste("scp -r -i ",sshpass_key_home," ",username,"@",home_computer,":",paths$cardamom,"LIBRARY/CARDAMOM_F ",ecdf_source,sep="")
                                        ,paste("scp -r ",username,"@",home_computer,":",paths$cardamom,"LIBRARY/CARDAMOM_F ",ecdf_source,sep="")
                                        # The CMake project root is one level above LIBRARY; copy the
                                        # top-level CMakeLists.txt so 'cmake -S <cardamom_ecdf>' can find it
                                        # (it does project() + add_subdirectory(LIBRARY/CARDAMOM_F)).
                                        ,paste("scp ",username,"@",home_computer,":",paths$cardamom,"CMakeLists.txt ",paths$cardamom_ecdf,"/",sep="")
                                        ,paste("rm -rf ",ecdf_source,"CARDAMOM_F/build",sep="")
                                        ,paste("cmake --build . --clean-first",sep="")
                                        ,paste("cmake -S ",paths$cardamom_ecdf," -B ",ecdf_source,"CARDAMOM_F/build",
                                               " -DCMAKE_BUILD_TYPE=",build_type,
                                               " -DMODEL=",modelname,
                                               " -DCMAKE_Fortran_COMPILER=",compiler,sep="")
                                        ,paste("cmake --build ",ecdf_source,"CARDAMOM_F/build --clean-first",sep="")                                        
                                        ,paste("cmake --build ",ecdf_source,"CARDAMOM_F/build --target ",target," -j",sep="")
                                        ,paste("cp ",ecdf_source,"CARDAMOM_F/executable/",built_exe," ",eexepath,"/",exe,sep="")
                                        ,paste("cp ",ecdf_source,"CARDAMOM_F/build","/LIBRARY/CARDAMOM_F/libCARDAMOM.so* ",eexepath,"/",sep="")
                                        ,paste("cp ",ecdf_source,"CARDAMOM_F/build","/LIBRARY/CARDAMOM_F/libcardamom-samplers-lib.so* ",eexepath,"/",sep="")))
              # If a crop model the copy the crop development files into place too
              if (modelname == "DALEC.A3.C3.H2.M1.015" | modelname == "DALEC.C3.M1.014") {
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

  } # on eddie

  # then compile locally
  comline="y"
  if (comline == "y") {
      print("Finally compile source code locally")
      if (project_src == "C") {
          system(paste("gcc ",paths$cardamom,"LIBRARY/CARDAMOM_C/projects/DALEC_CDEA_TEMPLATE/DALEC_CDEA_TEMPLATE.c -o ",
                       paths$cardamom,"LIBRARY/CARDAMOM_C/projects/DALEC_CDEA_TEMPLATE/a.out -lm",sep=""))
          system(paste("cp ",paths$cardamom,"LIBRARY/CARDAMOM_C/projects/DALEC_CDEA_TEMPLATE/a.out ",exepath,"/",exe,sep=""))
      } else if (project_src == "Fortran") {

          # The parallel-samplers code base is built with CMake (see
          # LIBRARY/CARDAMOM_F/CMakeLists.txt). We configure an out-of-source
          # build under LIBRARY/CARDAMOM_F/build and build the executable for the
          # requested sampler. The CMake project writes the executables and
          # dalec.so into LIBRARY/CARDAMOM_F/executable.
          if (Sys.which("cmake") == "") {stop("cmake not found on PATH (CMake >= 3.22 required to build CARDAMOM_F)")}

          # Fixed source, build and output locations
          srcroot   = paths$cardamom
          buildpath = paste(paths$cardamom,"LIBRARY/CARDAMOM_F/build",sep="")
          exedir    = paste(paths$cardamom,"LIBRARY/CARDAMOM_F/executable",sep="")

          # Map the debug / timing switches onto the CMake build type
          build_type = "RELEASE" ; if (timing | debug) {build_type = "DEBUG"}

          # Select the sampler executable to build:
          target = "cardamom"        ; built_exe = "cardamom.exe"

          if (request_compile_server == FALSE | request_compile_local == TRUE) {
              # The compiler, model and build type are baked into the CMake cache
              # at configure time and cannot be changed in place, so force a clean
              # configure by removing any existing build tree.
              if (dir.exists(buildpath)) {system(paste("rm -rf ",buildpath,sep=""))}
              # Remove any stale executable in the project folder
              if (file.exists(paste(exepath,"/",exe,sep=""))) {system(paste("rm ",exepath,"/",exe,sep=""))}
              if (file.exists(paste(exepath,"/libCARDAMOM.so",sep=""))) {system(paste("rm ",exepath,"/libCARDAMOM.so*",sep=""))}              
              # Configure (compiler / model / build type set here)
              system(paste("cmake -S ",srcroot," -B ",buildpath,
                           " -DCMAKE_BUILD_TYPE=",build_type,
                           " -DMODEL=",modelname,
                           " -DCMAKE_Fortran_COMPILER=",compiler,sep=""))
              system(paste("cmake --build ",buildpath," --clean-first",sep=""))                           
              # Build the chosen sampler executable
              system(paste("cmake --build ",buildpath," --target ",target," -j",sep=""))
              # Copy the built executable (source name = CMake OUTPUT_NAME) to the
              # project executable directory, renamed to the project executable name
              system(paste("cp ",exedir,"/",built_exe," ",exepath,"/",exe,sep=""))
              # Copy the shared library the executable is linked against
              system(paste("cp ",buildpath,"/LIBRARY/CARDAMOM_F/libCARDAMOM.so* ",exepath,"/",sep=""))    
              system(paste("cp ",buildpath,"/LIBRARY/CARDAMOM_F/libcardamom-samplers-lib.so** ",exepath,"/",sep=""))    
              # Build the shared library needed later by R for forward DALEC runs
              system(paste("cmake --build ",buildpath," --target dalec -j",sep=""))
              system(paste("cp ",exedir,"/dalec.so ",exepath,"/dalec.so",sep=""))
          } # compile executable on local machine too?

          # Copy crop development file into position
          if (modelname == "DALEC.A3.C3.H2.M1.015" | modelname == "DALEC.C3.M1.014") {
               system(paste("cp ",paths$cardamom,"LIBRARY/CARDAMOM_F/model/",modelname,"/src/winter_wheat_development.csv ",exepath,"/",sep=""))
          } #

    } else {

        stop('Source code language has not been specified')

    } # Compile locally

  } # copy and compile to eddie

  # prepare output
  PROJECT$ecdf=use_eddie
  PROJECT$request_use_local_slurm=request_use_local_slurm
  PROJECT$localpath=localpath
  PROJECT$datapath=datapath
  PROJECT$oestreampath=oestreampath
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
