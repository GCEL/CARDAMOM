
###
## Function to read parameter chains info
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).
# Exceptions are given within specific functions.

dump_binary_files <-function(infile) {

  # Function will read in the contents of a raw fortran binary file format and
  # dump into a big vector for the user to sort out themselves

  # Check file exists
  if (file.exists(infile)) {

      # open connection to the file 
      # 'r' to read and 'b' for binary
      bob = file(infile,'rb') ; nos_var = 1e6
      # Read from the file, double() is the data type expected
      set1 = readBin(bob, double(),nos_var) ; temp = 0
      # keep reading until we have read all that can be read
      while (length(temp) > 0) {
         temp = readBin(bob, double(),nos_var)
         set1 = append(set1,temp)
      }
      # now close this chain
      close(bob)

      # Tidy up
      rm(temp)

      # Back to user
      return(set1)

  } # file exists

} # end function

read_parameter_chains<- function(PROJECT_in,n,ndim) {

  # search for all output files
  pfile = list.files(paste(PROJECT_in$resultspath,sep=""), full.names=TRUE)
  # select the correct project
  is_it = grepl(PROJECT_in$name,pfile) ; pfile = pfile[is_it]
  # select the PARS files only
  is_it = grepl("PARS",pfile) ; pfile = pfile[is_it]
  # need to duplicate the list at this point to ensure that we can be certain we do not confuse the chain number and site numbers
  pfile_tmp = gsub(c("_PARS"),"",pfile)
  # select the correct site
  is_it = grepl(paste(PROJECT_in$name,"_",PROJECT_in$sites[n],"_",sep=""),pfile_tmp) ; pfile = pfile[is_it] ; rm(pfile_tmp)
  # Find and remove any files which have no data in them
  is_it = file.size(pfile) ; is_it = which(is_it > 0) ; pfile = pfile[is_it]

  # just in case
  if (length(pfile) < 1) {return(-9999)}

  # search for all output files
  sfile = list.files(paste(PROJECT_in$resultspath,sep=""), full.names=TRUE)
  # select the correct project
  is_it = grepl(PROJECT_in$name,pfile) ; pfile = pfile[is_it]
  # select the correct site
  is_it = grepl(PROJECT_in$sites[n],pfile) ; pfile = pfile[is_it]
  # select the STEP files only
  sfiles = paste(PROJECT_in$resultspath,PROJECT_in$name,"_",PROJECT_in$sites[n],"_*_STEP",sep="")

  # calculate the number of chains
  chains = seq(1,length(pfile))
  # load the fraction of samples to lose
  frac = as.numeric(PROJECT_in$latter_sample_frac)
  # calculate the number of parameter vectors this is
  par_vector_length = 100
  # which site are we on now
  print("Beginning parameter extraction and chain merge")
  print(paste("Site = ",PROJECT_in$sites[n]," ",n," of ",PROJECT_in$nosites," ",Sys.time(),sep=""))
  # create error flag, initial value zero
  status = array(0,dim=c(length(chains)))

  # define output variable
  param_sets_out = array(NA,dim=c((PROJECT_in$model$nopars[n]+1),par_vector_length,length(chains)))
  # loop through each chain
  for (c in seq(1,length(chains))) {
       print(paste("...chain ",c," of ",length(chains),sep=""))
       # open this chains binary file into R, instructing 'r' to read and 'b' for binary
       bob = file(paste(pfile[c],sep=""),'rb') ; nos_var = 1e6
       set1 = readBin(bob, double(),nos_var) ; temp = 0
       # keep reading until we have read all that can be read
       while (length(temp) > 0) {
          temp = readBin(bob, double(),nos_var)
          set1 = append(set1,temp)
       }
       # now close this chain
       close(bob)

       # re-arrange into correct structure
       param_sets = array(set1,dim=c((PROJECT_in$model$nopars[n]+1),(length(set1)/(PROJECT_in$model$nopars[n]+1))))
       set1 = 0 ; rm(set1)

       # check for inconsistencies
       if (abs(dim(param_sets)[2] - PROJECT_in$nsubsamples) > 1) {
           if (dim(param_sets)[2] > PROJECT_in$nsubsamples ) {
               print('*************************************************************')
               print(paste('Warning! Too many parameter vectors in ',pfile[c],sep=""))
               print('*************************************************************')
               print('Likely cause is appended solutions in previously existing file')
               print(paste('To solve this, only the last ',PROJECT_in$nsubsamples,' will be used',sep=""))
               # keep only the end of the parameter sets
               param_sets = param_sets[,((dim(param_sets)[2]-PROJECT_in$nsubsamples):dim(param_sets)[2])]
               status[c] = 2
           } else if (dim(param_sets)[2] < PROJECT_in$nsubsamples) {
               print('*************************************************************')
               print(paste('Warning! Missing parameter vectors in ',pfile[c],sep=""))
               print('*************************************************************')
               print('Likely cause is incomplete chain: it may have been stopped by ECDF, or')
               print('MCMC chain got stuck (happens sometimes). This may also result as an')
               print('inconsistency between requested number of samples vs PROJECT$nsubsamples')
               print('CONSIDER DELETING (or re-running) THIS CHAIN!!')
               print('Likely error will occur next!')
               status[c] = 1
               # just in case
               if (dim(param_sets)[2] < par_vector_length) {return(-9999)}
           } else {
               status[c] = 0
           }
       } # mismatch between expected and actual parameter outputs > 1

       # keep a sample form the end of the chain
       param_sets = param_sets[,(((dim(param_sets)[2]-par_vector_length)+1):dim(param_sets)[2])]
       # add these output to the final output variable
       param_sets_out[1:(PROJECT_in$model$nopars[n]+1),,c]=param_sets
       # clean up
       param_sets = 0 ; rm(param_sets)

  } # end of chains loop

  if (PROJECT_in$model$name == "DALEC_CDEA" || PROJECT_in$model$name == "DALEC_CDEA_LU_FIRES" ||
      PROJECT_in$model$name == "DALEC_CDEA_ACM2" || PROJECT_in$model$name == "DALEC_CDEA_ACM2_BUCKET" ||
      PROJECT_in$model$name == "DALEC_CDEA_ACM2_BUCKET_RmRg" || PROJECT_in$model$name == "DALEC_CDEA_ACM2_BUCKET_RmRg" ||
      PROJECT_in$model$name == "DALEC_CDEA_ACM2_BUCKET_RmRg_CWD_wMRT") {
      param_sets_out[c(12,15),,] = ((param_sets_out[c(12,15),,]-1)%%365.25)+1
  }
  if (PROJECT_in$model$name == "DALEC_CDEA_no_lit_root") {
      param_sets_out[c(8,11),,] = ((param_sets_out[c(8,11),,]-1)%%365.25)+1
  }

  # return the parameter solutions
  return(param_sets_out)

} # end of function
## Use byte compile
read_parameter_chains<-cmpfun(read_parameter_chains)
