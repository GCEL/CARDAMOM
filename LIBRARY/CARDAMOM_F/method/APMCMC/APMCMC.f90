!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CARbon DAta MOdel fraMework (CARDAMOM) and DALEC terrestrial ecosystem model suite
! CARDAMOM is a Bayesian model-data fusion software framework. CARDAMOM is used to
! assimilate observations and ecological theory to retrieve parameters for the
! DALEC suite of intermediate complexity terrestrial ecosystem models. DALEC can be
! used as a fully integrated component of CARDAMOM or independently.
! Copyright (C) 2024  University of Edinburgh,
!                     Mathew Williams (mat.williams@ed.ac.uk),
!                     T. Luke Smallman (t.l.smallman@ed.ac.uk)
! UoE = University of Edinburgh
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!!!!!!!!!!!! File specific description !!!!!!!!!!
! Module contains all subroutine and functions relevant specifically to the
! AP-MCMC method. The choice of EDC, likelihood and model are made else where and
! are thus contains within a seperate module
!
! On normalized parameters:
! "raw" values, as convenient for loglikelihood calculation and file writing,
! are the default  For function arguments, internal saved data.  Parameters
! are converted to normalized values as needed by adaptive statistics functions.  History matrix
!  PARSALL and covariance matrix work refer to the normalized parameter space.
!
! How to use this sampler:
!  Create an object of type PARINFO : number and bounds of parameters
!  and DEMCzOPT-options, containing as many or few of the fields as needed, the rest default to the
!                 default values in type definiton here
!  Create an object of type MCMC_OUTPUT to write reults to.
!  Create a real function loglikelihood taking a vector of npars (same as in PARINFO) parameters and
!                 returning rel:: loglikelihood.
!  (not implemented yet) Optionally set OMP_NUM_THREADS
!  Call subroutine DEMCz(fct, parinfo, mcopt, mcmcout)
!
! Relevant source references:
! Haario et al., (2001) An adaptive Metropolis algorithm. Bernoulli 7.2: 223-242.
! Haario et al., (2006) Stat. Comput., 16:339–354, DOI 10.1007/s11222-006-9438-0,
! Roberts and Rosenthal (2009), Examples of Adaptive MCMC, J. Comp. Graph. Stat. 18:349-367
!
! This code is based on the original C verion of the University of Edinburgh
! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
! All code translation into Fortran, integration into the University of
! Edinburgh CARDAMOM code and subsequent modifications by:
! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
! J. F. Exbrayat (University of Edinburgh)
! refactored jklebes 2024-25
! See function / subroutine specific comments for exceptions and contributors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module cardamom_MHMCMC

   use samplers_shared, only: PARINFO, MCMC_output, MCMC_options, filenames_insert_threadid, neg_inf
   use samplers_io, only: io_buffer_space, initialize_buffers, open_output_files, close_output_files
   use OMP_LIB

   implicit none(type, external)

   public
contains
   !
   !--------------------------------------------------------------------
   !
   subroutine run_parallel_mcmc(model_likelihood, PI, MCO, MCOUT_list, model_likelihood_write_in, nchains, seed, standard)
      !- Run multiple parallel MCMC simulations (adaptive MCMC algorithm with CARDAMOM-specific quirks)
      !
      !/* ***********INPUTS************
      ! *
      ! * model_likelihood : A function wholly responsible for
      ! * (a) running the model given the DATA and parameters,
      ! * (b) comparing it to observations, and
      ! * (c) returning  the (log) likelihood.
      ! * The subroutine will be run as MODEL_LIKELIHOOD(PARS, npars, loglikelihood, chainid)
      ! * and outputs a single loglikelihood value
      ! *
      ! * PI: This structure contains information on parameter number and bounds
      ! *
      ! * MCO: This structure contains option values for the MCMC run.
      ! * These will be set to default values if empty, restart, nchains . Options include:
      ! * (a) number of runs
      ! * (b) filename for writing file with results
      ! * (c) step adaptation frequency
      ! * (d) initial step size
      ! * */
      !
      !/* **************OUTPUTS*************
      ! *
      ! * RESULTS FILE: File includes (a) results (b) likelihood and (c) final step
      ! size
      ! *
      ! * */
      ! Write to MCMCOUT
      !
      !-

      implicit none(type, external)

      !! input and output structs
      ! read-only, shared beteen chains:
      type(PARINFO), intent(in):: PI
      !!PARINFO struct from model giving number, bounds of parameters
      type(MCMC_OPTIONS), intent(inout):: MCO
      !! struct of options for the run, shared between all threads
      type(MCMC_OUTPUT), dimension(:), allocatable, intent(inout):: MCOUT_list  ! array of MCOUT objects
      !! Array of MCMC_OUTPUT structs for each thread's results
      integer, optional, intent(in):: nchains
      integer, intent(in) :: seed
      logical, intent(in), optional :: standard
	!! Run standard adaptive MCMC algorithm, or include cardamom-specific modifications
        !! Default standard=false => use cardamom quirks

      logical :: standard_

      !! number chains optional, default 1
      integer:: i
      !! internal loop index

      !> the function to maximize.
      !> Completely agnostic, samples any functions vector -> double
      !> Typically a loglikelihood evaluation of a model against observation data
      !> given the inputted parameter values.
      interface
         subroutine model_likelihood(param_vector, n, ML, id) bind(c)
            use iso_c_binding, only:  c_int, c_double
            implicit none(type, external)
            real(c_double), dimension(n), intent(inout):: param_vector
			!! used as read-only, but has to be inout for compatibility with R via C
            integer(c_int), intent(in):: n, id
            real(c_double), intent(out):: ML
         end subroutine model_likelihood
      end interface

      !> optionally  give a second function with same shape as model_likelihood,
      !> for writing to file.  model_likelihood_write_in is the input arg, which may not be present.
      procedure(model_likelihood), optional:: model_likelihood_write_in
      !> A second function with same shape as model_likelihood,
      !> for writing to file.  Internal variable equal to model_likelihood_write_in if present
      !> or (default) same as main model_likelihood function
      procedure(model_likelihood), pointer:: model_likelihood_write

      ! Argument processing  !!!!!!!!!!!!!!!

      if (.not. present(standard)) then
         standard_  = .false.
      else
         standard_ = standard
      end if


      if (.not. present(nchains)) then
         MCO%nchains = 1
         ! or try to infer from OMP_THREADS or columns in MCOUT ...?
      else
         MCO%nchains = nchains
      end if

      ! process function arguments
      if (present(model_likelihood_write_in)) then
         ! if second function given use it for writing to file
         model_likelihood_write => model_likelihood_write_in
      else
         ! default, if no argument given: use same function for calculation and printing
         model_likelihood_write => model_likelihood
      end if

      ! Outputs
      ! each one an MCOUT struct with unitialized output data fields
      if (.not. allocated(MCOUT_list)) allocate (MCOUT_list(MCO%nchains))

      ! This adaptive MCMC is "parallelized" to do N chains at the same time, but they
      ! each do a complete independent run.  The parallelization structure of this module is
      ! trivial.  It exists mainly as template for samplers with more crossover and more complex
      ! structure.
      ! `shared` is already omp default behvaior, just making it explicit.  This is ok because
      ! all arguments are read-only / intent(in) except MCOUT_list(i), which has its own array entry.
      !$OMP parallel do default(shared)
      do i = 1, MCO%nchains
         ! saves latest, best loglikelihood and associated parameters to MCOUT_list(i)
         call run_mcmc(model_likelihood, PI, MCO, MCOUT_list(i), model_likelihood_write, i, seed+i, standard_)
      end do
      !$OMP end parallel do

   end subroutine run_parallel_mcmc

   subroutine run_mcmc(model_likelihood, PI, MCO, MCOUT, model_likelihood_write_in, chainid, seed, standard)
   !! Main function for a single-thread adaptive MCMC simulation
      use samplers_math, only: log_par2nor, log_nor2par, par2nor, nor2par
      use samplers_shared, only: init_pars_random, bounds_check,  metropolis_choice
      use samplers_io, only: write_parameters, write_variances, write_covariance_matrix &
                             , write_covariance_info, write_mcmc_output, open_output_files
      use random_uniform, ONLY: UNIF_VECTOR, initialize_random

      ! declare any local variables
      ! all variables in here are local to the single chain and the duration of its run
      ! input and output structs
      type(PARINFO), intent(in):: PI
      type(MCMC_OPTIONS), intent(in):: MCO
      type(MCMC_OUTPUT), intent(inout):: MCOUT

      integer, intent(in), optional:: chainid
      integer, intent(in) :: seed
      logical, intent(in), optional :: standard
      integer:: chainid_
      logical:: standard_

      type(io_buffer_space):: io_space
      !! buffer for writing to out files, private to this chain
      character(350):: outfile, stepfile, covfile, covifile
      !! filenames
      double precision, dimension(PI%npars) :: PARS_previous & ! parameter values for current state
                                              ,PARS_proposed & ! parameter values for current proposal
                                              ,BESTPARS        ! best set of parameters so far

      double precision, allocatable, dimension(:,:) :: PARSALL  ! The history, accepted normalised parameter values in intervals nadapt .  npars x nadapt
      double precision :: loglikelihood_previous, &
                          loglikelihood_proposed, & ! loglikelihood of a set of parameters
                            output_loglikelihood, & ! loglikelihood according to  alternative loglikelihood calculation for writing to file
                                           llmax, & ! best loglikelihood seen so far
                                  burn_in_period, & ! global TODO module data ?! for how many proposals will we adapt the covariance matrix as a minimum
                                     opt_scaling, & ! = opt_scaling_const/npars
                                            beta, &
                                 par_minstepsize, &
                                        P_target
      type(UNIF_VECTOR):: uniform_random_vector
      !! object holding array of pre-generated random values - (supposedly faster to pregenerate) - local to this
      !! chain
      logical:: multivariate
      integer:: MAXITER, nchains, npars
      !> counters-local to this chain's run
      integer:: ITER, ACC, ACC_FIRST, ACCLOC, N_before_mv_target
      integer :: pfileunit, sfileunit, cfileunit, cifileunit ! file unit numbers
      double precision:: ACCRATE, ACCRATE_GLOBAL

      ! declare interface for the model likelihood function.
      ! A function pars vector -> loglikelihood
      ! Do not modify here, sampler is generic.
      ! Wrap the model loglikelihood function in another function elsewhere to
      ! make it conform to this form.
      interface
         subroutine model_likelihood(param_vector, n, ML, id) bind(c)
            use iso_c_binding, only:  c_int, c_double
            implicit none(type, external)
            real(c_double), dimension(n), intent(inout):: param_vector
            integer(c_int), intent(in):: n, id
            real(c_double), intent(out):: ML
         end subroutine model_likelihood
      end interface

      ! optionally give a second function with same shape as model_likelihood,
      ! for writing to file;
      procedure(model_likelihood), optional:: model_likelihood_write_in
      ! Going forwards this pointer is the alternateive likelihood function for writing to file
      procedure(model_likelihood), pointer:: model_likelihood_write

      logical:: accept

      ! process arguments...

      if (.not. present(standard)) then
         standard_  = .false.
      else
         standard_ = standard
      end if


      if (present(model_likelihood_write_in)) then
         ! if second function given use it for writing to file
         model_likelihood_write => model_likelihood_write_in
      else
         ! default, if no argument given: use same function for calculation and printing
         model_likelihood_write => model_likelihood
      end if

      ! read settings from input struct ...
      npars = PI%npars
      nchains = MCO%Nchains
      MAXITER = MCO%nout
      P_target = MCO%P_target
      !MCOUT%use_multivariate  ! default .false.
      N_before_mv_target = MCO%N_before_mv*PI%npars
      beta = MCO%beta
      par_minstepsize = MCO%par_minstepsize

      if (present(chainid)) then
         chainid_ = chainid
      else
         chainid_ = 1
      end if

      ! Allocate PARSALL, the history of the chain
      allocate(PARSALL(PI%npars, MCO%nadapt))

      ! Initialize pregenerated random numbers, if using-local to this chain
      call uniform_random_vector%initialize_random(seed)

      ! process file names
      outfile = MCO%outfile
      stepfile = MCO%stepfile
      covfile = MCO%covfile
      covifile = MCO%covifile
      if (MCO%nchains > 1 .and. present(chainid)) then
         call filenames_insert_threadid(outfile, stepfile, covfile, covifile, chainid_)
      end if

    !!! prepare file writing
      if (MCO%nwrite > 0) then
         ! allocate buffers (different one for each chain)
         call initialize_buffers(npars, MAXITER/MCO%nwrite, io_space)
         call open_output_files(outfile, stepfile, covfile, covifile, chainid_, pfileunit, sfileunit, cfileunit, cifileunit)
      end if

    !! Calculate derived  settings of the run...
      ! Determine how long we will continue to adapt our proposal covariance
      ! matrix and use of Delayed Rejection
      burn_in_period = MCO%fADAPT*dble(MAXITER)
      ! See step() for relevant references.
      ! scd = 2.381204 the optimal scaling parameter for MCMC search, when applied
      ! to multivariate proposal.
      ! NOTE 1: 2.38/sqrt(npars) sometimes used when applied to the Cholesky factor
      ! NOTE 2: 2.381204**2 = 5.670132
      opt_scaling = sqrt(MCO%opt_scaling_const/dble(PI%npars))

      if (MCO%restart) then
         ! For aborted simulations
         ! keep MCOUT, starting point, counters statistics collection
         ! implies start from last state
         ITER = 0
         MCOUT%complete = .false.
         ACC = 0
         ACC_first = 0
         ACCLOC = 0
         ACCRATE = 0d0
         ACCRATE_GLOBAL = 0d0
      else
         ! Not a restart from aborted simulation
         ! start new counters

         ! init MCOUT
         MCOUT%complete = .false.
         MCOUT%nos_iterations = 0
         ! Initialize basic MCMC counters
         ITER = 0
         ! Initialize further counters for adaptive
         ACC = 0
         ACC_first = 0
         ACCLOC = 0
         ACCRATE = 0d0
         ACCRATE_GLOBAL = 0d0

         if (.not. MCO%fixedpars) then
            ! completely new simulation from new starting point

            if (.not. allocated(MCOUT%pars)) allocate (MCOUT%pars(npars))
            if (.not. allocated(MCOUT%bestpars)) allocate (MCOUT%bestpars(npars))

            ! initialize statistics collection output fields
            if (.not. allocated(MCOUT%parvar)) then
               ! we recieved blank new MCOUT, start new stats collection
               allocate (MCOUT%parvar(npars))
               allocate (MCOUT%meanpar(npars))
               allocate (MCOUT%covariance(npars, npars))
            end if

         else
            ! ((.not. restart) .and. fixedpars) : keep MCOUT%pars and statistics from previous phase
            ! write "initial" covariance matrix that we inherited from previous phase to file
            if (MCO%nwrite > 0) then
               call write_covariance_matrix(MCOUT%covariance, npars, .true., chainid_)
            end if
         end if

      end if

      !!!! init params
      ! initalize bestpars to current pars
      if (.not. MCO%restart .and. .not. MCO%fixedpars) then
         call init_pars_random(PI, PARS_previous, PI%fix_pars, uniform_random_vector)
         ! Inform the user
         write (*, *) "Have generated random initial pars"
      else  ! restart or fixedpars case-from aborted simulation or from previous phase-use pars in MCOUT as stating point
         PARS_previous = MCOUT%PARS
         write (*, *) "Starting from pars from previous phase"
         ! start statistics from previous statistics in MCOUT
      end if

      if (.not. MCO%restart) then  ! if new simulation, whether from fixedpars or random pars
         ! calculate initial ll
         call model_likelihood(PARS_previous, npars, loglikelihood_previous, chainid_)
         ! best value so far = starting value
         BESTPARS = PARS_previous
         llmax = loglikelihood_previous
      else
         ! TODO recalculatig because not read from file in case of restart , although
         ! present in column npars+1
         ! loglikelihood_previous = MCOUT%ll
         call model_likelihood(PARS_previous, npars, loglikelihood_previous, chainid_)
         ! BESTPARS = MCOUT%bestll  ! not available from file read
         BESTPARS = PARS_previous
         ! recalculating this one too
         call model_likelihood(BESTPARS, npars, llmax, chainid_)
      end if

      if (loglikelihood_previous == neg_inf) then
         write (*, *) "WARNING  ! loglikelihood = ", loglikelihood_previous, " - &
         & AP-MCMC will get stuck, if so please check initial conditions"
         error stop 1
      end if

      ! Begin the main AP-MCMC loop
      do while (ITER < MAXITER .and. loglikelihood_previous < (P_target - 10*epsilon(1.0d0)))

         ! set flag for mutlivariate phase
         multivariate = MCOUT%use_multivariate .and. (MCOUT%Nparvar > N_before_mv_target)

         ! take a step in parameter space: generate proposed
         ! new parameters PARS
         call step_pars_real(PARS_previous, PARS_proposed, PI, multivariate, MCOUT%covariance, beta, opt_scaling, par_minstepsize, &
                             uniform_random_vector)

         ! if parameter proposal in bounds check the model
         if (bounds_check(PI, PARS_proposed)) then
            ! calculate the model likelihood
            call model_likelihood(PARS_proposed, npars, loglikelihood_proposed, chainid_)
            accept = metropolis_choice(loglikelihood_proposed, loglikelihood_previous)
         else
            accept = .false.
         end if  ! in bound

         if (accept) then
            ! Keep count of the number of accepted proposals in this local period
            ACCLOC = ACCLOC + 1
            ! Store accepted parameter proposals (unnormalized values)
            ! keep record of all parameters accepted since step adaption
            ! (this chain)
            ! Because this history matrix is used for (normalized) statistics for adaptiveness,
            ! store normalized version of pars
            PARSALL(1:npars, ACCLOC) = log_par2nor(PARS_proposed, PI%parmin, PI%parmax, PI%paradj)  ! add row in history matrix
            ! Accepted first proposal from multivariate
            if (multivariate) ACC_first = ACC_first + 1

            ! Save accepted pars as previous pars, save previous log-likelihood
            PARS_previous(1:npars) = PARS_proposed(1:npars)
            loglikelihood_previous = loglikelihood_proposed
            ! store the best parameter set
            if (loglikelihood_previous >= llmax) then
               BESTPARS = PARS_previous
                llmax = loglikelihood_previous
            end if
         else
            ! here we write the same state to history again, in a standard MCMC
            if (standard_) then
                ACCLOC = ACCLOC + 1
                PARSALL(1:npars, ACCLOC) = log_par2nor(PARS_previous, PI%parmin, PI%parmax, PI%paradj)
            endif
            ! Cardamom quirk 1:  not done in CARDAMOM-MHMCMC version to match original
            ! Likely the samplers lose ergodicity.  Helps samplers not get stuck, from experience.
         end if  ! accept or reject proposed pars

         ! count iteration
         ITER = ITER + 1

         if (MCO%nwrite > 0) then
            if (mod(ITER, MCO%nwrite) == 0) then
               ! calculate the likelihood for the actual uncertainties-this avoid
               ! issues with different phases of the MCMC which may use sub-samples
               ! of observations or inflated uncertainties to aid parameter
               ! searching
               call model_likelihood_write(PARS_previous, npars, output_loglikelihood, chainid_)
               ! Now write out to files
               call write_mcmc_output(MCOUT%parvar, ACCRATE, &
                                      MCOUT%covariance, &
                                      MCOUT%meanpar, MCOUT%Nparvar, &
                                      PARS_previous, output_loglikelihood, npars, ITER == MAXITER, &
                                      io_space, chainid_)
            end if
         end if  ! write or not to write

         ! time to adapt?
         if (mod(ITER, MCO%nadapt) == 0) then

           !! update the acceptance counters and acceptance ratios
            ! Total accepted values
            ACC = ACC + ACCLOC
            ! Calculate global acceptance rate
            ACCRATE_GLOBAL = ACC/ITER
            ! Calculate local acceptance rate (i.e. since last adapt)
            ACCRATE = ACCLOC/dble(MCO%nadapt)

            ! Second, are we still in the adaption phase?
            if (burn_in_period > ITER .or. .not. MCOUT%use_multivariate) then
               if (standard_) then
               ! PARSALL and ACCLOC stay as they are - full history and count of all states sampled
               ! in this phase i.e. since last (mod(ITER, MCO%nadapt) == 0)
               ! to be added to running statistacs calculations
               ! there should be MCO%nadapt new rows in this matrix
               else

               ! Cardamom quirk 2: initial estimate comes from first, middle, and last 3 states out of the last
               ! period of nadapt steps only.
               ! then (new 2025) continue with all accepted samples from last period.
               if (.not. MCOUT%cov .and. ACCLOC > 3) then
                  PARSALL(1:npars, 2) = PARSALL(1:npars, ceiling(ACCLOC*0.5d0))
                  PARSALL(1:npars, 3) = PARSALL(1:npars, ACCLOC)
                  ACCLOC = 3
                  ! leads to attempting to make new covariance matrix from the first, middle, and last rows of local history
               end if
               endif

               ! update the covariance matrix for multivariate proposal
               call update_statistics(PARSALL, npars, MCOUT, MCOUT%use_multivariate, ACCLOC, N_before_mv_target, standard_)

               if (MCO%nwrite > 0 .and. ACCLOC == 3 .and. MCOUT%cov) then
                  ! we just switched from not having a covraiance matrix to having an initial covariance matrix
                  ! write the first covariance matrix to file
                  call write_covariance_matrix(MCOUT%covariance, npars, .true., chainid_)
               end if

            end if !  have enough parameter been accepted

            ! resets the local acceptance counter
            ACCLOC = 0

         end if  ! time to adapt?

         ! Should I be write(*,*)ing to screen or not?
         if (MCO%nPRINT > 0) then
            if (mod(ITER, MCO%nPRINT) == 0) then
               write (*, *) "Chain ", chainid, "of", mco%nchains
               write (*, *) "Using multivariate sampling = ", MCOUT%use_multivariate
               write (*, *) "Total proposal = ", ITER, " out of ", MAXITER
               write (*, *) "Total accepted = ", ACC
               write (*, *) "Overall acceptance rate    = ", dble(ACC)/dble(ITER)
               write (*, *) "Local   acceptance rate    = ", ACCRATE
               write (*, *) "Current obs   = ", loglikelihood_previous, "proposed = ", loglikelihood_proposed, " log-likelihood"
               write (*, *) "Maximum likelihood = ", llmax
               ! NOTE: that-infinity in current obs only indicates failure of EDCs
               ! but-infinity in both obs and parameter likelihood scores indicates
               ! that proposed parameters are out of bounds
            end if
         end if  ! write(*,*) to screen or not

      end do  ! end MHMCMC sampling loop conditions

    !!!  Finalize

      ! write out final covariance matrix for the analysis
      if (MCO%nwrite > 0) call write_covariance_matrix(MCOUT%covariance, npars, .false., chainid_)

      ! record the best single set of parameters
      MCOUT%bestpars = BESTPARS
      MCOUT%bestll = llmax
      ! Output the current state, so that next simulations can potentially resume here
      MCOUT%pars = PARS_previous  !! (!) TODO should be done outside
      MCOUT%ll = loglikelihood_previous
      ! record how many iterations were taken to complete
      MCOUT%nos_iterations = MCOUT%nos_iterations + ITER
      ! set flag MCMC completed
      MCOUT%complete = .true.
      ! tidy up
      deallocate(PARSALL)

      ! completed AP-MCMC loop
      write (*, *) "ITER ", ITER
      write (*, *) "AP-MCMC loop completed"
      write (*, *) "Overall acceptance rate     = ", dble(ACC)/dble(ITER)
      write (*, *) "Final local acceptance rate = ", ACCRATE
      write (*, *) "Best log-likelihood = ", llmax
      !write (*, *) "Best parameters = ", MCOUT%bestpars

      if (MCO%nwrite > 0) then
         call close_output_files(pfileunit, sfileunit, cfileunit, cifileunit)
      endif

   end subroutine run_mcmc
   !
   !------------------------------------------------------------------
   !
   subroutine update_statistics(PARSALL, npars, MCOUT, use_multivariate, nadapt, N_before_mv_target, standard)
      !! Calculate/update running mean, variance, and covariance.
      !! Formerly adapt_step_size, this function is central to the adaptive method from Roberts and Rosenthal 2009
      !! which adjusts proposal step size proportional to observed variances/covariances.
      !! three CARDAMOM quirks here : 1) only accepted steps were recorded, not repeat entries for non-accepted steps
      !!                              2) initial update_statistics is passed only first or first, middle, and last rows
      !!                              3) Modified increment_covariance_matrix call with artificially low 4th argument means more
      !!                                 distant history is downweighted
      use samplers_math, only: cholesky_factor, covariance_matrix, &
                               increment_covariance_matrix

      implicit none(type, external)

      ! declare input types
      type(MCMC_OUTPUT), intent(inout):: MCOUT  ! incl statistics collection
      logical, intent(inout):: use_multivariate
      !! covariance matrix exists and is useable, i.e. positive definite
      ! declare inputs variables
      integer, intent(in):: npars
      !! total length of history of simulation (including the nadapt recent steps)
      integer, intent(in):: nadapt
      !! nadapt = ACCLOC "accepted local", number of new states to add into running statistics
      double precision, intent(in):: PARSALL(npars, nadapt)  ! collection of recently accepted normalised parameter combinations

      ! declare local variables
      integer :: p, info  ! counters
      double precision, dimension(npars, npars):: cov_backup
      double precision, dimension(npars, npars):: cholesky
      double precision, dimension(npars):: meanpar_backup
      integer:: Nparvar_backup
	!! number of accepted steps as of before the current round of nadapt , = ITER - nadapt
        !!     Read it from MCOUT%Nparvar where it was saved at the end of last roung
      integer:: Nparvar_local
        !! min(Nparvar_backup, N_bfore_mv_target) ,
        !! a version of Nparvar_backup that is artificially capped at a constant
      integer, intent(in):: N_before_mv_target
      logical, intent(in) :: standard ! use standard literature adaptive MCMC (true)
                          ! or apply downweighting of history (false) for increased adaptiveness

      ! if we have a covariance matrix then we want to update it, if not then we need to create one
      if (MCOUT%cov) then

         ! Increment the variance-covariance matrix with new accepted parameter sets
         ! NOTE: that this also increments the total accepted counter (PI%Nparvar)

         cov_backup = MCOUT%covariance; meanpar_backup = MCOUT%meanpar; Nparvar_backup = MCOUT%Nparvar
         if (standard) then
		 ! update statistics :
                 !  Requesting running update of covariance and mean with 4th argument correctly
                 ! Giving number of states that entered into the calculation so far,
		 ! Nparvar_backup
		 ! = ITER-nadapt instead of being artificially capped at 100
		 Nparvar_local = Nparvar_backup
		 call increment_covariance_matrix(PARSALL, MCOUT%meanpar, npars &
		                                  , Nparvar_local, nadapt, MCOUT%covariance)
         else
		 ! 3rd difference between cardamom APMCMC and standard MHMCMC :
		 ! Have started hardcoding a maximum number of observations to be N_before_mv_target.
		 ! While not strictly following Haario et al., (2001) or Roberts and Rosenthal, (2009)
		 ! this allows for the covariance matrix to be more responsive to its local environment.
		 ! in fact this is having the effect of extremely downweighting history
		 ! in the running calculations of mean and covariance matrix in new covariance matrix .
		 ! They will not converge, and represent mean and covariance for the local neighborhood.

		 Nparvar_local = min(N_before_mv_target, Nparvar_backup)

         endif
		 ! caution : subroutine increment_covariance_matrix changes not just its last argument 'covariance',
		 ! but also its second argument 'mean' and its
		 ! 4th argument 'Nparvar_local' .
         call increment_covariance_matrix(PARSALL(1:npars, 1:nadapt), MCOUT%meanpar, npars &
		                                  , Nparvar_local, nadapt, MCOUT%covariance)
         ! Calculate the cholesky factor as this includes a determination of
         ! whether the covariance matrix is positive definite.
         cholesky = MCOUT%covariance
         ! caution: writes to its second argument, in addition to checking positive definiteness
         call cholesky_factor(npars, cholesky, info)
         ! If the updated covariance matrix is not positive definite we should
         ! reject the update in favour of the existing matrix
         if (info == 0) then
            ! Set multivariate sampling to true
            use_multivariate = .true.
            MCOUT%Nparvar = Nparvar_local
         else
            ! The current addition of a parameter leads to a matrix which is not
            ! positive definite. If we previously had a matrix which is positive
            ! definite then we should reject totally the new matrix, if not we
            ! should keep it and accumulate the information
            if (use_multivariate) then
               ! return original matrix to place
               MCOUT%covariance = cov_backup
               MCOUT%meanpar = meanpar_backup
               MCOUT%Nparvar = Nparvar_backup
            else
               ! Keep accumulating use_multivariatethe information
               use_multivariate = .false.
            end if
         end if

      else  ! PI%cov == .false.

         ! we have not yet created a covariance matrix based on accepted
         ! parameters. Assuming we have some then create one...
         if (nadapt > 2) then

            ! estimate covariance matrix
            call covariance_matrix(PARSALL(1:npars, 1:nadapt), MCOUT%meanpar, &
            & npars, nadapt, MCOUT%covariance)
            MCOUT%cov = .true.; MCOUT%Nparvar = nadapt

            ! Calculate the cholesky factor as this includes a determination of
            ! whether the covariance matrix is positive definite.
            ! Caution: cholesky_factor alters its second argument
            ! that's why we only input a copy of the matrix
            cholesky = MCOUT%covariance
            call cholesky_factor(npars, cholesky, info)

            if (info /= 0) then
               ! Keep accumulating information until positive definite matrix
               ! calculated
               use_multivariate = .false.
            else
               ! Positive definite found straight away-might as well use it!
               use_multivariate = .true.
            end if

         end if  ! nadapt > 2

      end if  ! PI%cov == .true.

      ! variance is the diagonal of the covariance matrix
      do p = 1, npars
         MCOUT%parvar(p) = MCOUT%covariance(p, p)
      end do

      return

   end subroutine update_statistics
   !
   !--------------------------------------------------------------------
   !
   subroutine step_pars_real(PARS0, PARS, PI, multivariate, covariance, beta, opt_scaling, &
                             par_minstepsize, random_uniform_vector)
      use samplers_math, only: log_par2nor, log_nor2par
      use random_uniform, only: UNIF_VECTOR
      implicit none(type, external)

      !> Generates new proposed state from currect state in real parameter space.
      !> Wraps step_pars

      ! Arguements
      double precision, dimension(:), intent(in):: pars0    ! current parameters
      double precision, dimension(:), intent(out):: pars    ! proposal
      type(UNIF_VECTOR), intent(inout):: random_uniform_vector
      !integer, intent(in):: npars
      type(PARINFO), intent(in):: PI
      logical, intent(in):: multivariate
      double precision, dimension(:,:), intent(in):: covariance
      double precision, intent(in):: beta, opt_scaling, par_minstepsize

      ! Local variables
      double precision, dimension(PI%npars) :: pars0_norm, pars_norm

      ! Normalise the current parameters
      pars0_norm = log_par2nor(pars0, PI%parmin, PI%parmax, PI%paradj)
      ! Step in normalised space
      call step_pars(pars0_norm, pars_norm, PI%npars, multivariate, covariance, beta, opt_scaling, par_minstepsize, &
                     random_uniform_vector)
      ! Return to real parameter values for the assessment
      pars = log_nor2par(pars_norm, PI%parmin, PI%parmax, PI%paradj)

   end subroutine step_pars_real
   !
   !------------------------------------------------------------------
   !
   subroutine step_pars(PARS0, PARS, npars, multivariate, covariance, beta, opt_scaling, &
                        par_minstepsize, random_uniform_vector)  ! TODO check against original !!
      use samplers_math, only: random_normal, random_multivariate
      use random_uniform, only: UNIF_VECTOR

      ! Applies Roberts and Rosenthal 2009-Eq 3 to generate new proposed state in (normalized)
      ! parameter space
      ! IN : Pars0 current state (normalized)
      ! IN : Stats MCOSTATS object holding info like covariance matrix of the run so far (calculated on
      ! normalized space )
      ! OUT: PARS new proposed state (normalized)
      ! plus take beta from module data

      implicit none(type, external)

      ! declare input variables
      !double precision, dimension(PI%npars), intent(inout):: !norpars0 & ! normalised current parameters
      !,norpars  & ! normalised proposal
      double precision, dimension(:), intent(in):: pars0 ! current parameters
      double precision, dimension(:), intent(out):: pars ! proposed new parameters to generate
      integer, intent(in):: npars ! length of pars vectors
      type(UNIF_VECTOR), intent(inout):: random_uniform_vector ! this chain's uniform random number object
      logical, intent(in):: multivariate ! Are we in multivariate sampling phase, i.e. guided by complete covariance matrix
                                         ! or just variances ?
      double precision, dimension(:, :), intent(in):: covariance ! covariance matrix as measured by simulation so far
      double precision, intent(in):: beta, opt_scaling, par_minstepsize ! parameters of the adaptive algorithm

      ! declare local variables
      integer:: p ! loop counter
      double precision:: rn(npars), mu(npars), rn2(npars) ! internal random vectors

      ! mean of distributions
      mu = 0d0

      ! Splitting step calculation based on number of parameter vectors accepted
      ! is linked to the need build a covariance matrix prior to multivariate
      ! sampling.
      ! See Roberts and Rosenthal, Examples of Adaptive MCMC, J. Comp. Graph. Stat. 18:349-367, 2009.

      ! Sample random normal distribution (mean = 0, sd = 1)
      ! TODO vectorize
      do p = 1, npars
         call random_normal(random_uniform_vector, rn2(p))
      end do

      if (multivariate) then

         ! Draw a vector from multivariate distribution
         ! NOTE: if covariance matrix provided is not positive definite
         !       a sample from normal distribution is returned
         call random_multivariate(npars, 1, covariance, mu, rn, random_uniform_vector)

         ! Estimate the step to be applied to the current parameter vector to
         ! create the new proposal. scd = a scaling parameter linking searching
         ! stepping to the number of parameters being retrieved by the analysis.
         ! See Haario et al., (2001) An adaptive Metropolis algorithm. Bernoulli 7.2: 223-242.
         ! and references therein. See also, Roberts & Rosenthal (2009) for beta scaling.
         pars = pars0 + (rn*opt_scaling*(1d0 - beta)) + (par_minstepsize*rn2*beta)

      else

         !MCOUT%multivariate_proposal = .false.
         pars = pars0 + (par_minstepsize*rn2)

      end if

   end subroutine step_pars
   !
   !------------------------------------------------------------------
   !
end module cardamom_MHMCMC
