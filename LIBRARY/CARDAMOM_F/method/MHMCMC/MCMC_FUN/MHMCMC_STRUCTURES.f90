
module MCMCOPT

  ! module contains the main derived types required by the MCMC itself. Input met
  ! drivers and observational data are however stored else where in a cardamom
  ! related module

  !!!!!!!!!!!
  ! Authorship contributions
  !
  ! This code is based on the original C verion of the University of Edinburgh
  ! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
  ! All code translation into Fortran, integration into the University of
  ! Edinburgh CARDAMOM code and subsequent modifications by:
  ! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
  ! J. F. Exbrayat (University of Edinburgh)
  ! See function / subroutine specific comments for exceptions and contributors
  !!!!!!!!!!!

  implicit none

  ! assume all private
  private

  ! explicit public statements
  public :: MCMC_OPTIONS, MCMC_OUTPUT ,COUNTERS, PARAMETER_INFO &
           ,PI,MCOUT,MCO, initialise_mcmc_output

  !
  ! module contains the declarable types used in the MCMC
  !

  ! contains MHMCMC options
  type MCMC_OPTIONS

    double precision :: sub_fraction = 0.2d0, inflation_factor = 6d0

    logical :: sub_sample_complete = .false. & ! Has a sub-sample / inflated uncertainty simulation taken place?
              ,returnpars = .true. & ! return best fit parameters or not
              ,randparini = .true. & ! use random initial values parameters
              ,fixedpars  = .true.   ! use fixed initial values where inputs are not = -9999
    character(350) :: outfile   & ! parameter output file name
                     ,stepfile  & ! step file name
                     ,covfile   & ! covariance file name
                     ,covifile    ! covariance information file
    double precision :: fADAPT    ! adapt step size for a given fraction of the full run
    integer :: nADAPT     & ! adapt step size after every N iterations
              ,nOUT       & ! number of requested output parameter sets
              ,nPRINT     & ! print info to screen every N solutions (0 to silent)
              ,nWRITE     & ! write to file every N solutions
              ,append       ! append to existing output files (0 = delete existing file or 1 = append)

  end type ! MCMC_OPTIONS
  ! create options type
  type(MCMC_OPTIONS), save :: MCO

  ! contains output information
  type MCMC_OUTPUT
    integer :: complete ! is MHMCMC completed (1 = yes, 0 = no)
    ! further metrics could be added here
    integer :: nos_iterations
    double precision :: acceptance_rate
    double precision, allocatable, dimension(:) :: best_pars ! store current best parameter set
  end type ! MCMC_OUTPUT
  ! create output type
  type(MCMC_OUTPUT), save :: MCOUT

  ! information which is needed determine progress
  type COUNTERS

    double precision :: ACC            & ! total number of accepted solutions
                       ,ACC_first      & ! total number of pre-DR accepted solution
                       ,ACC_beta       & ! total number of beta accepted solutions
                       ,ACCLOC         & ! number of recently accepted solutions
                       ,ACCLOC_beta    & ! number of recently accepted solutions from beta proposals
                       ,ITER_beta      & ! number of iterations using beta proposal
                       ,ITER           & ! number of iterations attempted
                       ,ACCEDC         & ! number of EDC complient iterations
                       ,ACCRATE_beta   & ! local beta step acceptance rate
                       ,ACCRATE        & ! local acceptance rate
                       ,ACCRATE_GLOBAL   ! global acceptance rate

  end type ! COUNTERS

  ! parameter structure defined here
  type PARAMETER_INFO
    logical :: cov = .false., use_multivariate = .false.
    double precision, allocatable, dimension(:,:) :: covariance, & ! parameter covariance matrix
                                                             iC    ! inverse covariance matrix
    double precision, allocatable, dimension(:) :: mean_par & ! mean parameter value
                                                  ,parmax   & ! maximum parameter values
                                                  ,parmin   & ! minimum parameter values
                                                  ,parini   & ! initial parameter values
                                                  ,paradj   & ! adjustment to allow log-normalised stepping
                                                  ,parfix   & ! do they need fixing (i.e. randomly generated)
                                                  ,parvar     ! variance of accepted parameter

    double precision :: Nparvar ! Number of samples forming variance of accepted parameters

    integer :: npars  ! number of parameters to be solved

    ! crop specific variables
    double precision :: stock_seed_labile
    double precision, allocatable, dimension(:)  ::    DS_shoot, & !
                                                        DS_root, & !
                                                       fol_frac, & !
                                                      stem_frac, & !
                                                      root_frac, & !
                                                        DS_LRLV, & !
                                                           LRLV, & !
                                                        DS_LRRT, & !
                                                           LRRT

  end type ! PARAMETER_INFO
  ! create parameter info type
  type (PARAMETER_INFO), save :: PI

  contains
  !
  !------------------------------------------------------------------
  !
  subroutine initialise_mcmc_output

    ! subroutine allocated memory to the output arrays

    implicit none

    ! define dimensions for output structure
    allocate(MCOUT%best_pars(PI%npars))

    ! set to zero, will become 1 when MCMCM complete
    MCOUT%complete = 0
    ! and clear initial memory
    MCOUT%best_pars = 0d0

  end subroutine initialise_mcmc_output
  !
  !------------------------------------------------------------------
  !
end module
