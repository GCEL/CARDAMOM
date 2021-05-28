
###
## Function determines whether or not parametr chains have statistically converged
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).
# Exceptions are given within specific functions.

psrf<- function (X) {
      #PSRF Potential Scale Reduction Factor
      #
      #   [R,NEFF,V,W,B] = PSRF(X) or
      #   [R,NEFF,V,W,B] = PSRF(x1,x2,...,xs)
      #   returns "Potential Scale Reduction Factor" (PSRF) for
      #   collection of MCMC-simulations. X is a NxDxM matrix
      #   which contains M MCMC simulations of length N, each with
      #   dimension D. MCMC-simulations can be given as separate
      #   arguments x1,x2,... which should have the same length.
      #
      #   Returns
      #     R     PSRF in a row vector of length D
      #     neff  estimated effective number of samples M*N*V/B
      #     V     estimated mixture-of-sequences variances
      #     W     estimated within sequence variances
      #     B     estimated between sequence variances
      #
      #   The idea of the PSRF is that if R is not near 1 (below 1.1 for
      #   example) one may conclude that the tested samples were not from
      #   the same distribution (chain might not have been converged
      #   yet).
      #
      #   If only one simulation is given, the factor is calculated
      #   between first and last third of the chain. Note that use of
      #   only one chain will produce over-optimistic result.
      #
      #   Method is from:
      #      Brooks, S.P. and Gelman, A. (1998) General methods for
      #      monitoring convergence of iterative simulations. Journal of
      #      Computational and Graphical Statistics. 7, 434-455. Note that
      #      this function returns square-root definiton of R (see Gelman
      #      et al (2003), Bayesian Data Analsyis p. 297).
      #
      #   See also
      #     CPSRF, MPSRF, IPSRF

      # Copyright (C) 1999 Simo Särkkä
      # Copyright (C) 2003 Aki Vehtari
      #
      # This software is distributed under the GNU General Public
      # Licence (version 2 or later); please refer to the file
      # Licence.txt, included with the software, for details.

      # 2004-01-22 Aki.Vehtari@hut.fi Added neff, R^2->R, and cleaning

      # M = nos_chains
      # N = nos_timesteps
      # D = nos_parameters

      # Therefore X = N,D,M
      N=dim(X)[1] ; D=dim(X)[2] ; M=dim(X)[3]

      # must have more than 1 time step for variance to be assessed
      if (N < 1) {stop('Too few samples')}

      # Calculate means W of the variances
      W = array(0,dim=c(1,D))
      for (n in seq(1,M)) { ####arrays don't match here try in matlab first
	x = X[,,n] - matrix(colMeans(X[,,n]),nrow=N,ncol=length(colMeans(X[,,n])),byrow=TRUE)
	W = W + colSums(x*x)
      }
      W = W / ((N-1) * M)

      # Calculate variances B (in fact B/n) of the means.
      Bpn = array(0,dim=c(1,D))
      m = rowMeans(colMeans(X))
      for (n in seq(1,M)) {
	x = colMeans(X[,,n]) - m
	Bpn = Bpn + x*x
      }
      Bpn = Bpn / (M-1)

      # Calculate reduction factors
      S = (N-1)/N * W + Bpn
      R = (M+1)/M * S / W - (N-1)/M/N
      V = R * W
      R = sqrt(R)
      B = Bpn*N
      neff = min(M*N*V/B,M*N)
      output=list(R=R,neff=neff,V=V,W=W,B=B)
      return(output)
}
## Use byte compile
psrf<-cmpfun(psrf)
