
###
## Function to determine convergence of chains
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

have_chains_converged<-function (param_sets) {

      # input is order dimensions(npar+1,iter,chain)

      # lets not beat about things lets see if the models have converged
      # need to re-arrange the array to the same as the function needs
      var = array(NA,dim=c(dim(param_sets)[2],dim(param_sets)[1],dim(param_sets)[3]))
      for (i in seq(1, dim(param_sets)[1])) {
	         var[1:dim(param_sets)[2],i,1:dim(param_sets)[3]] = param_sets[i,1:dim(param_sets)[2],1:dim(param_sets)[3]]
      }

      # pass parameters to convergence function
      converged = psrf(var) #; print(converged$R[length(converged$R)])
      # assume critical value of 1 (default is 1.1 but we want to be a little conservative)
      bob = which(converged$R > 1.2)
      # start with all pass assumption
      converged = array("PASS",length(converged$R))
      # changed passes to fails if needed
      converged[bob] = "FAIL"
      # return the result
      return(converged)

} # end of function
