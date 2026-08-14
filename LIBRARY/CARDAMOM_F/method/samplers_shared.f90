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

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

!!!!!!!!!!!! File specific description !!!!!!!!!!
! This file contains the source code of various different 
! declarable types used in the MCMC solvers.
!
! This code was implemented by Jason Klebes (Jason.Klebes@ed.ac.uk)
! Subsequent modifications by:
! T. L. Smallman (University of Edinburgh, t.l.smallman@ed.ac.uk)
!
! See function / subroutine specific comments for exceptions and contributors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module samplers_shared
   
   implicit none(type, external)

   public

   private filename_insert_threadid_single

   !> A collection of info about the model's parameters
   !> npars and min, max bounds as two arrays
   !> Closely related to model fct; model fct must take this number
   !> and type of paramters
   !> Unlike in previous versions, intended to be intent(in) only
   !> TODO we could bundle this in a type with the function:
   !> the fortran quasi-object
   type PARINFO
      integer:: npars
      double precision, allocatable, dimension(:):: parmin, parmax, paradj
      logical, allocatable, dimension(:):: fix_pars
      logical, allocatable, dimension(:):: parfix
   end type PARINFO

   !> A collection of input options to the MCMC sampler run
   !> contains default values
   type SAMPLER_options
      integer:: MAXITER = 10000, & ! overall steps, if convergence not reached
                 nadapt = 1000,  & ! steps per "local" sampling period, between adaptation steps
                Nchains = 1,     & ! consider setting OMP env to something compatible
                 nwrite = 1000,  & ! Frequency (steps) of writing current parameters to file
                 nprint = 1000     ! Frequency (steps) of prining current solver information to screen
      integer:: nout !
      double precision:: P_target = 0d0 ! termination criterion-a loglikelihood to stop at (optional)

      !> file names for outputs, note output format is a raw binary format.
      character(350)::  outfile = "parout.bin", &
                      stepfile = "stepout.bin", &
                        covfile = "covout.bin", &
                   covifile = "covinfoout.bin"
      real:: fadapt  ! TODO fraction adapt-move to outside
      logical:: append, & !
            randparini, & !
            returnpars    ! a variable that is never used and has no effect, needs deleting in all model likelihood files
      logical:: restart = .false., & ! is it a restart ?
              fixedpars = .false.    ! Continue from last state in MCOUT (don't initialize to random points) ?
   end type SAMPLER_OPTIONS

   type , extends (SAMPLER_OPTIONS) :: MCMC_OPTIONS ! Additional settings for adaptive-MCMC family samplers
      !> setting for adaptive AP-MCMC step size
      double precision:: par_minstepsize = 0.001d0 & ! 0.0005 -> 0.001 -> 0.01 -> 0.1 -> 0.005
                        ,par_maxstepsize = 0.01d0  & !
                       ,par_initstepsize = 0.005d0 & !
                                   ,beta = 0.05d0    ! weighting for gaussian step in multivariate proposals
      !> Optimal scaling variable for parameter searching
      double precision:: opt_scaling_const = 2.381204**2 ! scd = 2.381204 the optimal scaling parameter
                                                         ! for MCMC search, when applied to  multivariate proposal.
                                                         ! NOTE 1: 2.38/sqrt(npars) sometimes used when applied to the Cholesky factor.
                                                         ! NOTE 2: 2.381204**2 = 5.670132
      integer :: N_before_mv = 10 ! Number of accepted proposals before attempting to build multi-variate sampler
   end type MCMC_OPTIONS


   !> Collection of info for output of the sampling run
   !> , can also be passed to next run to continue from the last state
   !> Note output is mainly via file writing
   type MCMC_OUTPUT
      double precision:: bestll ! best (maximum) loglikelihood value found so far
      double precision:: ll ! latest loglikelihood value
      double precision, allocatable, dimension(:):: bestpars ! best (loglikelihood-maximizing) parameter values found so far
      double precision, allocatable, dimension(:):: pars     ! latest parameter values
      double precision:: acceptance_rate ! acceptance rate
      logical:: complete ! Did the main loop finish?
      integer:: nos_iterations = 0 ! number main loop interations run so far
      integer:: Nparvar, Nparvar_local ! Number of states in history that have gone into running
                                        ! mean, variance, and covariance calculation.
      double precision, allocatable, dimension(:):: parvar ! variance
      double precision, allocatable, dimension(:):: meanpar ! mean of parameters
      double precision, allocatable, dimension(:, :):: covariance ! covariance matrix measured during sampling
      logical:: cov = .false. ! Does the covariance matrix exist yet?
      logical:: use_multivariate = .false. ! Covariance matrix exists and is positive definite
   end type MCMC_OUTPUT

   double precision:: neg_inf = -huge(1d0)

contains
   !
   !--------------------------------------------------------------------
   !
   subroutine init_infinity()
      ! TODO not used.  `infini=0` is still used in model/.
      use, intrinsic :: ieee_arithmetic, only: ieee_negative_inf, ieee_support_inf, ieee_value

      if (ieee_support_inf(1d0)) then
           neg_inf = ieee_value(1d0, ieee_negative_inf)
      else
           neg_inf = -huge(1d0)
      end if
   end subroutine init_infinity
   !
   !--------------------------------------------------------------------
   !
   logical function metropolis_choice(new_loglikelihood, old_loglikelihood)

      ! core sampler math
      !> Accept or reject
      !> First argument is new/proposed log(!)likelihood,
      !> second is old log likelihood
      !> return logical

      ! Arguments
      double precision, intent(in):: new_loglikelihood, old_loglikelihood

      ! Local variables
      double precision:: r  ! draw random number 0 to 1
      call random_number(r)

      ! TODO add optional pregen random
      ! l1/l2 > r  <=> logl1-logl2 > log(r)
      ! TODO very small chance of r = exactly 0 .  (is this true with our generator)  
      ! should catch and supply log(r) = -inf .  Performance  impact of check?
      metropolis_choice = ((new_loglikelihood - old_loglikelihood) > log(r))
 
   end function metropolis_choice
   !
   !--------------------------------------------------------------------
   !
   ! routines for random initialization
   !
   !--------------------------------------------------------------------
   !
   subroutine init_pars_random(PI, pars0, fix_pars_flag, uniform_random_vector)
      use random_uniform, only: UNIF_VECTOR, next_random_uniform
      use samplers_math, only: log_nor2par
      implicit none(type, external)
      
      ! Arguments
      type(PARINFO), intent(in):: PI  ! give number, bounds of params
      double precision, dimension(PI%npars), intent(inout):: pars0  ! return random initial values-nonnormalized
      logical, dimension(PI%npars), optional, intent(in):: fix_pars_flag  ! flags .true. to keep inidividual pars
      type(UNIF_VECTOR), optional, intent(inout):: uniform_random_vector  ! object supplying pre-generated randoms 0 to 1  ! TODO make optional

      ! Local variables
      logical, dimension(PI%npars):: fix_pars_flag_  ! internal version
      integer:: i

      if (.not. (present(fix_pars_flag))) then
         fix_pars_flag_ = .false.
      else
         fix_pars_flag_ = fix_pars_flag
      end if

      do i = 1, PI%npars
         ! parfix = 1 stay at prior value, parfix = 0 randomly search
         ! TODO condense this horrible logic to restart and fix_pars_flag somewhere outside of this fct
         !if (MCO%fixedpars .and. PI%parini(i) /= -9999d0) PI%parfix(i) = 1d0
         ! only assign random parameters if (a) randparini == .true. or (b) PI$parini(n) == -9999)
         !if (MCO%randparini .and. PI%parfix(i) == 0d0 .and. .not.restart_flag) then
!
         !TODO vectorize
         if (.not. (fix_pars_flag_(i))) then
            ! make sure to give it a random number vector unique to the chain
            ! scale each to the parameter's range
            pars0(i) = log_nor2par(uniform_random_vector%next_random_uniform(), PI%parmin(i), PI%parmax(i), PI%paradj(i))
         end if

      end do  ! for PI%npar loop
      !TODO test
      !TODO merge, optional uniform_random_vector arg
   end subroutine init_pars_random
   !
   !--------------------------------------------------------------------
   !
   subroutine init_latin_square(PI, pars0, n_chains)  ! TODO
      use samplers_math, only: nor2par
      implicit none(type, external)

      ! Arguements
      type(PARINFO), intent(in):: PI  ! give number, bounds, and potentially current value of params
      integer, intent(in):: n_chains
      double precision, dimension(PI%npars, n_chains), intent(out):: pars0  ! return initial values-nonnormalized

      ! Local variables
      integer:: i, j
      double precision, dimension(PI%npars, n_chains):: points

      ! generate N initial points on the (0..1)^N space in latin hypercube distribution...
      ! Must be run for all chains at once outside of parallel regions

      ! convert to real parameter values space
      do j = 1, n_chains
         do i = 1, PI%npars
            pars0(i, j) = nor2par(points(i, j), PI%parmin(i), PI%parmax(i))
         end do
      end do

   end subroutine init_latin_square
   !
   !--------------------------------------------------------------------
   !
   pure logical function bounds_check(PI, PARS)

      ! Arguments
      type(PARINFO), intent(in):: PI
      double precision, dimension(:), intent(in):: PARS

      ! given real-space params ... convert to lognorm space, check if all in (0, 1)
      ! or check directly against real boundary values
      bounds_check = all((PARS > PI%parmin) .and. (PARS < PI%parmax))

   end function bounds_check
   ! 
   !--------------------------------------------------------------------
   !
   subroutine filenames_insert_threadid(outfile, stepfile, covfile, covifile, chainid)
      !! Amends the given 4 filenames by inserting chainid in the appropriate place
      !! e.g. stem_COV -> stem_1_COV
      
      ! Arguments
      character(len=*), intent(inout):: outfile, stepfile, covfile, covifile
      integer, intent(in):: chainid
      
      call filename_insert_threadid_single(outfile, chainid)
      call filename_insert_threadid_single(stepfile, chainid)
      call filename_insert_threadid_single(covfile, chainid)
      call filename_insert_threadid_single(covifile, chainid)

   end subroutine filenames_insert_threadid
   !
   !--------------------------------------------------------------------
   !
   subroutine filename_insert_threadid_single(filename, chainid) 

      ! inserts thread id into a single filename,
      ! STEM_FILETYPE -> STEM_NUMBER_FILETYPE

      ! Arguements
      character(len=*), intent(inout):: filename  ! expected format STEM_FILETYPE eg "UK_baseline_sites_AliceHolt_COV" 
      integer, intent(in):: chainid

      ! Local variables
      integer :: index_split
      character(4):: chainid_str ! internal char version of chainid number, for filenames

      ! internal write to convert int -> str
      write (chainid_str, '(i0)') chainid

      ! separate filename into stem and suffix again at last _
      index_split = scan(filename, '_', back=.true.) !location of last _     
      ! assemble new filename
      filename =  trim(filename(1:index_split))//trim(chainid_str)//"_"//trim(filename((index_split+1):))
 
   end subroutine filename_insert_threadid_single
   !
   !--------------------------------------------------------------------
   !
end module samplers_shared
