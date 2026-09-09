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
# Function to write commands need to run and submit processes to local interactive machine
# 
# This function is based on an original Matlab function development by A. A. Bloom 
# (UoE, now at the Jet Propulsion Laboratory). Translation to R and subsequent 
# modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).
# Exceptions are given within specific functions.
#
#########################################################################################

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

    # Check presence of PROJECT_in$request_cost_function_scaling
    if (exists(x = "request_cost_function_scaling", where = PROJECT_in) == FALSE) {
        # If not, assume default cost function
        PROJECT_in$request_cost_function_scaling = 0
    }

    # Each site is now a single process that runs all PROJECT_in$nochains chains
    # internally via OpenMP. Request that thread count for the executable.
    Sys.setenv(OMP_NUM_THREADS = as.integer(PROJECT_in$nochains))

    # begin submitting the different tasks
    cwd = getwd()
    setwd(PROJECT_in$exepath)
    for (n in seq(1, PROJECT_in$nosites)) {
         # Throttle: every concurrent_sites-th job runs in the foreground so we do
         # not launch more concurrent sites than the core budget allows.
         run_background = TRUE ; if (n%%concurrent_sites == 0) {run_background = FALSE}
         # Create the input / output file names for the current job.
         infile=paste(PROJECT_in$datapath,PROJECT_in$name,"_",PROJECT_in$sites[n],".bin",sep="")
         # Output base only; the chain id and PARS/STEP/COV/COVINFO suffix are
         # appended by the executable (one set of files per internal chain).
         output=paste(PROJECT_in$resultspath,PROJECT_in$name,"_",PROJECT_in$sites[n],sep="")
         # Build the command (8 positional arguments, the last being the chain count)
         command=paste(PROJECT_in$exepath,PROJECT_in$exe," ",
                       infile," ",
                       output," ",
                       as.integer(nsamples),
                       " 0 ",
                       as.integer(PROJECT_in$samplerate)," ",
                       as.integer(pre_mcmc)," ",
                       as.integer(PROJECT_in$request_cost_function_scaling)," ",
                       as.integer(PROJECT_in$nochains),sep="")
         # Depending on whether or not a background run submit the job
         if (background == "y" | run_background) {
             system(paste(command," & ",sep=""))
         } else {
             system(command)
         }
         # To ensure that each site is submitted at a unique time we want to delay the code - this impacts the seed value used in the random number generator
         Sys.sleep(1) # wait for 1 seconds
    } # site loop
    setwd(cwd) ; rm(cwd)
    print("Command issued to local machine")

} # end of function submit_processes_to_local_machine

## Use byte compile
submit_processes_to_local_machine<-cmpfun(submit_processes_to_local_machine)
