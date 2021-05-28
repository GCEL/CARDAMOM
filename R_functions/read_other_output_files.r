
###
## Function to read final parameter covariance matrix from file
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).
# Exceptions are given within specific functions.

read_other_output_files<- function(PROJECT_in,n,pattern_wanted) {

  # search for all output files
  cfile = list.files(paste(PROJECT_in$resultspath,sep=""), full.names=TRUE)
  # select the correct project
  is_it = grepl(PROJECT_in$name,cfile) ; cfile = cfile[is_it]
  # select the COV files only
  is_it = grepl(paste("_",pattern_wanted,sep=""),cfile) ; cfile = cfile[is_it]
  # finally select the specific site we are interested in
  is_it = grepl(PROJECT_in$sites[n],cfile) ; cfile = cfile[is_it]

  # just in case
  if (length(cfile) <= 1) {return(-9999)}

  # calculate the number of chains
  chains = seq(1,length(cfile))
  # which site are we on now
  print("Beginning parameter covariance extraction")
  print(paste("Site = ",PROJECT_in$sites[n]," ",n," of ",PROJECT_in$nosites," ",Sys.time(),sep=""))
  # create error flag, initial value zero
  status = array(0,dim=c(length(chains)))

  # define the output varable
  output = array(NA, dim=c(PROJECT_in$model$nopars[n]+1,PROJECT_in$nsubsamples,length(chains)))

  # loop through each chain
  for (c in seq(1,length(chains))) {
    print(paste("...chain ",c," of ",length(chains),sep=""))
    # open this chains binary file into R, instructing 'r' to read and 'b' for binary
    bob = file(paste(cfile[c],sep=""),'rb') ; nos_var = 1e6
    set1 = readBin(bob, double(),nos_var) ; temp = 0
    if (length(set1) > 0) {
       # keep reading until we have read all that can be read
       while (length(temp) > 0) {
         temp = readBin(bob, double(),nos_var)
         set1 = append(set1,temp)
       }
       # We have both, re-arrange into correct structure
       tmp = length(set1)/(PROJECT_in$model$nopars[n]+1)
       temp = array(set1, dim=c(PROJECT_in$model$nopars[n]+1,tmp))
       tmp = min(tmp, PROJECT_in$nsubsamples)
       output[,1:tmp,c] = temp[,1:tmp]
    }
    # now close this chain
    close(bob) ; set1 = 0 ; rm(set1)

  } # end of chains loop

  # return the parameter solutions
  return(output)

} # end of function
