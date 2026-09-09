module test_APMCMC
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use test_functions
  use test_math, only: approx
  use samplers_math, only: random_int
  use random_uniform
  use APMCMC
  use OMP_LIB
  implicit none
  private

  public:: collect_APMCMCtests

contains

!> Collect all exported unit tests
subroutine collect_APMCMCtests(testsuite)
  implicit none
  !> Collection of tests
  type(unittest_type), allocatable, intent(out):: testsuite(:)

  testsuite = [ &
    new_unittest("step_pars", test_step_pars), &
    new_unittest("step_pars_real", test_step_pars_real), &
    new_unittest("mcmc_output_type", test_mcmc_output_type), &
    new_unittest("mcmc_len0", test_run_mcmc_len0), &
    new_unittest("mcmc_not_static", test_not_static), &
    new_unittest("mcmc_len1000", test_run_mcmc_len1000), &
    new_unittest("mcmc_nchains1_len0", test_run_parallel_mcmc_nchains1_len0), &
    new_unittest("mcmc_nchains4_len0", test_run_parallel_mcmc_nchains4_len0), &
    new_unittest("mcmc_nchains4omp_len0", test_run_parallel_mcmc_nchains4_enforceomp_len0), &
    new_unittest("mcmc_nchains4omp_len100000", test_run_parallel_mcmc_nchains4_enforceomp_len100000), &
    new_unittest("mcmc_stop_condition", test_mcmc_stop_condition), &
    new_unittest("mcmc_stop_condition_exact", test_mcmc_stop_condition_exact), &
    new_unittest("mcmc_stop_condition_initial", test_mcmc_stop_condition_initial), &
    new_unittest("mcmc_two_phase", test_mcmc_two_phase), &
    new_unittest("mcmc_two_phase_parallel", test_mcmc_two_phase_parallel) &
    ]

end subroutine collect_APMCMCtests


subroutine test_step_pars(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  ! generate new proposal PARS given current point PARS0 and other specs
  ! example point x, y = 0, 0 with independent variances
  integer, parameter:: npars = 2
  double precision, dimension(npars):: pars0, pars
  logical:: multivariate = .true.
  double precision, dimension(npars, npars):: covariance
  double precision:: beta = 0.05d0  ! recommended constant, only relevant in multivariate phase
  double precision:: opt_scaling = 5.67/dble(npars)  ! recommended constant
  double precision:: par_minstepsize = 0.001d0 ! ?? what does this do ?
  type(UNIF_VECTOR):: random_uniform
  integer:: seed
  seed = random_int()
  call random_uniform%initialize_random(seed)
  call init_PI()
  pars0 = [0d0, 0d0]
  covariance = reshape(source = [1d0, 0d0, 0d0, 1d0], shape = [2, 2])  ! uncorellated with large variance
  call step_pars(pars0, pars, npars, multivariate, covariance, beta, opt_scaling, par_minstepsize, random_uniform)
  ! expect pars is close to, but not equal pars0
  call check(error, (pars(1)/=0d0) .and. (pars(2)/=0d0))
  ! Note step_pars does not respect bouns 0..1 of normed space !!  as per originnal cardamom, checked and rejected later
  ! in bounds
  ! call check(error, (pars(1)>=0d0) .and. (pars(1)<=1d0))
  ! call check(error, (pars(2)>=0d0) .and. (pars(2)<=1d0))
end subroutine test_step_pars

subroutine test_step_pars_real(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  ! generate new proposal PARS given current point PARS0 and other specs
  ! example point x, y = 0, 0 with independent variances
  double precision, dimension(PI_xy%npars):: pars0, pars
  logical:: multivariate = .true.
  double precision, dimension(PI_xy%npars, PI_xy%npars):: covariance
  double precision:: beta = 0.05d0  ! recommended constant, only relevant in multivariate phase
  double precision:: opt_scaling  ! recommended constant
  double precision:: par_minstepsize = 0.001d0 ! ?? what does this do ?
  type(UNIF_VECTOR):: random_uniform
  integer:: seed
  seed = random_int()
  opt_scaling = 5.67/dble(PI_xy%npars)
  call random_uniform%initialize_random(seed)
  call init_PI()
  pars0 = [0d0, 3.0d0]  ! within bounds of PI_xy
  covariance = reshape(source = [1d0, 0d0, 0d0, 1d0], shape = [2, 2])  ! uncorellated with large variance
  call step_pars_real(pars0, pars, PI_xy, multivariate, covariance, beta, opt_scaling, par_minstepsize, random_uniform)
  ! expect pars is close to, but not equal pars0
  call check(error, (pars(1)/=pars0(1)) .and. (pars(2)/=pars0(2)))
  ! Note step_pars does not respect bouns 0..1 of normed space !!  as per originnal cardamom, checked and rejected later
  ! in bounds
  !call check(error, (pars(1)>=PI_xy%parmin(1)) .and. (pars(1)<=PI_xy%parmax(1)))
  !call check(error, (pars(2)>=PI_xy%parmin(2)) .and. (pars(2)<=PI_xy%parmax(2)))
end subroutine test_step_pars_real

subroutine test_adapt_step_size(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  ! ...
  ! expect either a covariance matrix exists or pi%use_multivariate is set to false
  ! TODO should be called update_covariance
end subroutine test_adapt_step_size

subroutine test_increment_covariance_matrix(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  ! call increment_covariance_matrix(history, mean, N_local, N_history, cov)
  ! expect changes to mean, cov
  ! and n_local +=1
end subroutine test_increment_covariance_matrix

subroutine test_cholesky_factor(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  !call cholesky_factor(N, cov, posdef)
end subroutine test_cholesky_factor

subroutine test_covariance_matrix(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  !call covariance_matrix(history, mean, npars, N_history, cov)
  ! expect output at mean, cov
end subroutine test_covariance_matrix

subroutine test_mcmc_output_type(error)
  use APMCMC, only: MCMC_OUTPUT
  implicit none
  type(error_type), allocatable, intent(out):: error
  ! test declare an object of type MCMC_OUTPUT
  type(MCMC_OUTPUT):: MCOUT
end subroutine test_mcmc_output_type

! TODO code to run before tests, such as init_pi

subroutine test_run_mcmc_len0(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  type(mcmc_output):: mcout
  type(mcmc_options):: mcopt  ! filled with defaults only
  ! zero length run : takes expected input arguments, setup works,
  ! outputs/writes unchanged state
  ! all on defaults, without optional arguments
  integer:: seed
  seed = random_int()
  call init_pi()
  mcopt%nout = 0
  call run_mcmc(ll_normal, pi_xy, mcopt, mcout, seed=seed)
  ! expect values in mcout : a random initial state (within given parameter bounds)
  ! and its loglikelihood
  ! call check(error, mcout%ll > 0d0, .true. )
  write(*,*) mcout%ll, mcout%pars
  call check(error, mcout%pars(1) >= pi_xy%parmin(1) .and. mcout%pars(1) <= pi_xy%parmax(1))
  call check(error, mcout%pars(2) >= pi_xy%parmin(2) .and. mcout%pars(2) <= pi_xy%parmax(2) )
  ! given starting values which happen to be the correct (min loglik) solution,
  ! expect after zero run steps mcout has these initial values and
  ! lolglik is low.
  !call check(error, mcout%pars(1), x_ideal )
  !call check(error, mcout%pars(2), y_ideal)
  !call check(error, MCOUT% ll, 0.0 )
end subroutine test_run_mcmc_len0

subroutine test_not_static(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  type(mcmc_output):: mcout
  type(mcmc_options):: mcopt  ! filled with defaults only
  double precision, dimension(2):: state0, state1, state2
    !! Evolving values of x, y estimate of normal distribution ll_normal
  integer:: seed
  seed = random_int()
  call init_pi()
  mcopt%nout = 0
  ! a Zero-length run to intialize to random state
  call run_mcmc(ll_normal, pi_xy, mcopt, mcout, seed=seed)
  state0 = mcout%pars
  mcopt%fixedpars = .true. ! start from same state
  mcopt%nout = 100
  ! A short run
  call run_mcmc(ll_normal, pi_xy, mcopt, mcout, seed=seed)
  state1 = mcout%pars
  ! another 100 sampling steps
  call run_mcmc(ll_normal, pi_xy, mcopt, mcout, seed=seed)
  state2 = mcout%pars
  ! check that there were fluctuations
  call check(error, state1(1) /= state0(1) .and. state1(2) /= state0(2))
  call check(error, state2(1) /= state1(1) .and. state2(2) /= state1(2))
end subroutine test_not_static

subroutine test_run_mcmc_len1000(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  type(mcmc_output):: mcout
  type(mcmc_options):: mcopt  ! filled with defaults only
  ! all on defaults, without optional arguments
  integer:: seed
  seed = random_int()
  call init_pi()
  mcopt%nout = 1000
  call run_mcmc(ll_normal, pi_xy, mcopt, mcout, seed=seed)
  ! expect values in mcout : loglikelihood and parameters in bounds
  ! call check(error, mcout%ll > 0d0 )
  call check(error, mcout%pars(1) >= pi_xy%parmin(1) .and. mcout%pars(1) <= pi_xy%parmax(1))
  call check(error, mcout%pars(2) >= pi_xy%parmin(2) .and. mcout%pars(2) <= pi_xy%parmax(2) )
  ! this is an easy problem (quadratic potential), expect convergence on solution in a short run
  !write(*,*) mcout%ll, mcout%pars
  !call check(error, mcout%pars(1), x_ideal )
  !call check(error, mcout%pars(2), y_ideal)
  !call check(error, MCOUT%ll, 50.0 )
end subroutine test_run_mcmc_len1000

subroutine test_run_parallel_mcmc_nchains1_len0(error)
  implicit none
  integer, parameter:: nchains = 1
  type(error_type), allocatable, intent(out):: error
  type(mcmc_output), dimension(:), allocatable:: mcout
  type(mcmc_output):: mcout1
  type(mcmc_options):: mcopt  ! filled with defaults only
  ! zero length run : takes expected input arguments, setup works,
  ! outputs/writes unchanged state
  ! all on defaults, without optional arguments
  integer:: seed
  seed = random_int()
  call init_pi()
  mcopt%nout = 0
  call run_parallel_mcmc(ll_normal, pi_xy, mcopt, mcout, seed=seed)
  ! expect values in mcout : a random initial state (within given parameter bounds)
  ! and its loglikelihood
  mcout1 = mcout(1)
  call check(error, mcout1%pars(1) >= pi_xy%parmin(1) .and. mcout1%pars(1) <= pi_xy%parmax(1))
  call check(error, mcout1%pars(2) >= pi_xy%parmin(2) .and. mcout1%pars(2) <= pi_xy%parmax(2) )
end subroutine test_run_parallel_mcmc_nchains1_len0

subroutine test_run_parallel_mcmc_nchains4_len0(error)
  implicit none
  integer, parameter:: nchains = 4
  type(error_type), allocatable, intent(out):: error
  type(mcmc_output), dimension(:), allocatable:: mcout
  type(mcmc_output):: mcout1
  type(mcmc_options):: mcopt  ! filled with defaults only
  integer:: i
  ! zero length run : takes expected input arguments, setup works,
  ! outputs/writes unchanged state
  ! all on defaults, without optional arguments
  integer:: seed
  seed = random_int()
  call init_pi()
  mcopt%nout = 0
  write(*,*) "calling"
  call run_parallel_mcmc(ll_normal, pi_xy, mcopt, mcout, nchains = nchains, seed=seed)
  ! expect values in mcout : a random initial state (within given parameter bounds)
  ! and its loglikelihood
  do i = 1, nchains
    mcout1 = mcout(i)
    write(*,*) i, mcout1%pars(1), mcout1%pars(2)
    call check(error, mcout1%pars(1) >= pi_xy%parmin(1) .and. mcout1%pars(1) <= pi_xy%parmax(1))
    call check(error, mcout1%pars(2) >= pi_xy%parmin(2) .and. mcout1%pars(2) <= pi_xy%parmax(2) )
  end do
end subroutine test_run_parallel_mcmc_nchains4_len0


subroutine test_run_parallel_mcmc_nchains4_enforceomp_len0(error)
  implicit none
  integer, parameter:: nchains = 4
  type(error_type), allocatable, intent(out):: error
  type(mcmc_output), dimension(:), allocatable:: mcout
  type(mcmc_options):: mcopt  ! filled with defaults only
  type(mcmc_output):: mcout1
  integer:: i
  ! zero length run : takes expected input arguments, setup works,
  ! outputs/writes unchanged state
  ! all on defaults, without optional arguments
  integer:: seed
  seed = random_int()
  call omp_set_num_threads(4)
  call init_pi()
  mcopt%nout = 0
  call run_parallel_mcmc(ll_normal, pi_xy, mcopt, mcout, nchains = nchains, seed=seed)
  ! expect values in mcout : a random initial state (within given parameter bounds)
  ! and its loglikelihood
  do i = 1, nchains
    mcout1 = mcout(i)
    write(*,*) i, mcout1%pars(1), mcout1%pars(2)
    call check(error, mcout1%pars(1) >= pi_xy%parmin(1) .and. mcout1%pars(1) <= pi_xy%parmax(1))
    call check(error, mcout1%pars(2) >= pi_xy%parmin(2) .and. mcout1%pars(2) <= pi_xy%parmax(2) )
  end do
end subroutine test_run_parallel_mcmc_nchains4_enforceomp_len0

subroutine test_run_parallel_mcmc_len1000(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  type(MCMC_OUTPUT):: MCOUT
  type(MCMC_OPTIONS):: MCOPT
  ! all on defaults, without optional arguments
  integer:: seed
  seed = random_int()
  MCOPT%nout = 1000
  call run_mcmc(ll_normal, PI_xy, MCOPT, MCOUT, seed=seed)
  ! Expect x and y close to true values were found
  ! and best loglik is close to 0
  ! TODO "is close"  helper
  call check(error, MCOUT%pars(1)== x_ideal )
  call check(error, MCOUT%pars(2)== y_ideal)
  call check(error, MCOUT%ll == 0.0 )
  ! This is an easy problem, expect convergence was reached long before MAXITER.
end subroutine test_run_parallel_mcmc_len1000

subroutine test_run_parallel_mcmc_nchains4_enforceomp_len100000(error)
  implicit none
  integer, parameter:: nchains = 4
  type(error_type), allocatable, intent(out):: error
  type(mcmc_output), dimension(:), allocatable:: mcout
  type(mcmc_options):: mcopt  ! filled with defaults only
  type(mcmc_output):: mcout1
  integer:: i
  ! zero length run : takes expected input arguments, setup works,
  ! outputs/writes unchanged state
  ! all on defaults, without optional arguments
  integer:: seed
  seed = random_int()
  call omp_set_num_threads(4)  ! error if not compiled with omp library
  call init_pi()
  mcopt%nout = 10000
  call run_parallel_mcmc(ll_normal, pi_xy, mcopt, mcout, nchains = nchains, seed=seed)
  ! expect values in mcout : a random initial state (within given parameter bounds)
  ! and its loglikelihood
  do i = 1, nchains
    mcout1 = mcout(i)
    write(*,*) i, mcout1%pars
    call check(error, mcout1%pars(1) >= pi_xy%parmin(1) .and. mcout1%pars(1) <= pi_xy%parmax(1))
    call check(error, mcout1%pars(2) >= pi_xy%parmin(2) .and. mcout1%pars(2) <= pi_xy%parmax(2) )
    ! expect bestll and bestpars better than final one (unless they happen to be the same one, unlikely)
    call check(error, mcout1%bestll > mcout1%ll )
    call check(error, (abs(mcout1%pars(1) - x_ideal) > abs(mcout1%bestpars(1)-x_ideal)) &
    & .or.  (abs(mcout1%pars(2) - y_ideal) > abs(mcout1%bestpars(2)-y_ideal) ))
    ! outputs complete and nos_iterations
    call check(error, mcout1%complete)  ! expect .true.
    call check(error, mcout1%nos_iterations, mcopt%nout)  ! because no convergence checks implemented at the moment
  end do
end subroutine test_run_parallel_mcmc_nchains4_enforceomp_len100000

subroutine test_mcmc_stop_condition(error)
  !! stops on MCMC stop criterion of reaching loglikelihood threshold
  !! and %bestll is expected to be same as latest ll
  implicit none
  type(error_type), allocatable, intent(out):: error
  type(mcmc_output):: mcout
  type(mcmc_options):: mcopt  ! filled with defaults only
  integer:: maxsteps
  integer:: seed
  seed = random_int()
  call init_pi()
  maxsteps = 1000000  ! don't expect to actually run for this long before convergence ll = 0.0
  mcopt%nout = maxsteps
  mcopt%P_target = -0.0d0  ! convergence criteria : loglikelood reached 0
  call run_mcmc(ll_step, pi_xy, mcopt, mcout, seed=seed)
  ! Expect we have optimized the values to the correct ranges
  write(*,*) mcout%pars, x_lower, x_upper, y_lower, y_upper
  call check(error, mcout%pars(1) >= x_lower .and. mcout%pars(1) <= x_upper )
  call check(error, mcout%pars(2) >= y_lower  .and. mcout%pars(2) <= y_upper )
  ! Expect the optimization stopped early
  call check(error, mcout%nos_iterations < maxsteps)
  ! Expect bestll is exactly 0 and current ll is exactly zero
  call check(error, approx(mcout%bestll, 0d0))
  call check(error, approx(mcout%ll, 0d0))
end subroutine test_mcmc_stop_condition

subroutine test_mcmc_stop_condition_exact(error)
  !! stops on MCMC stop criterion of reaching loglikelihood threshold
  !! with only exact equality to likelihood threshold, check 0d0 >= 0d0
  !! and %bestll is expected to be same as latest ll
  implicit none
  type(error_type), allocatable, intent(out):: error
  type(mcmc_output):: mcout
  type(mcmc_options):: mcopt  ! filled with defaults only
  integer:: maxsteps
  integer:: seed
  seed = random_int()
  call init_pi()
  maxsteps = 1000000  ! don't expect to actually run for this long before convergence ll = 0.0
  mcopt%nout = maxsteps
  mcopt%P_target = -0.0d0  ! convergence criteria : loglikelood reached 0
  call run_mcmc(ll_step, pi_xy, mcopt, mcout, seed=seed)
  ! Expect we have optimized the values to the correct ranges
  call check(error, mcout%pars(1) >= x_lower .and. mcout%pars(1) <= x_upper )
  call check(error, mcout%pars(2) >= y_lower  .and. mcout%pars(2) <= y_upper )
  ! Expect the optimization stopped early
  call check(error, mcout%nos_iterations < maxsteps)
  ! Expect bestll is exactly 0 and current ll is exactly zero
  call check(error, approx(mcout%bestll, 0d0))
  call check(error, approx(mcout%ll, 0d0))
end subroutine test_mcmc_stop_condition_exact

subroutine test_mcmc_stop_condition_initial(error)
  !! stops on MCMC stop criterion of reaching loglikelihood threshold
  !! with a threshold such that the initial state likely already fulfils
  !! the criterion and no sampling iterationsare run.  Still have complete output
  !! and %bestll is expected to be same as latest ll
  implicit none
  type(error_type), allocatable, intent(out):: error
  type(mcmc_output):: mcout
  type(mcmc_options):: mcopt  ! filled with defaults only
  integer:: maxsteps
  integer:: seed
  seed = random_int()
  call init_pi()
  maxsteps = 1000000  ! don't expect to actually run for this long before convergence ll = 0.0
  mcopt%nout = maxsteps
  mcopt%P_target = -10000.0d0  ! Extremely broad convergence criterion, already fulfilled
  call run_mcmc(ll_step, pi_xy, mcopt, mcout, seed=seed)
  ! Expect the optimization stopped early, immediately
  call check(error, mcout%nos_iterations < maxsteps)
  call check(error, mcout%nos_iterations <= 0)
  ! Expect bestll is current ll
  call check(error, approx(mcout%bestll, mcout%ll, 0d0))
end subroutine test_mcmc_stop_condition_initial

subroutine test_mcmc_two_phase(error)
  !! Reproducing find_edc+main MCMC cardamom run.
  !! First find a point in ll = 0 region of ll_step potential, then
  !! start at this point with ll_bounded potential, which returns
  !! -Infinity if started at points outside of this region
  implicit none
  type(error_type), allocatable, intent(out):: error
  type(mcmc_output):: mcout
  type(mcmc_options):: mcopt  ! filled with defaults only
  integer:: maxsteps
  double precision, dimension(2):: startingpars
  integer:: seed
  seed = random_int()
  call init_pi()
  maxsteps = 1000000  ! don't expect to actually run for this long before convergence ll = 0.0
  mcopt%nout = maxsteps
  mcopt%P_target = 0.0d0
  call run_mcmc(ll_step, pi_xy, mcopt, mcout, seed=seed)
  ! Expect we have optimized the values to the correct ranges
  call check(error, mcout%pars(1) >= x_lower .and. mcout%pars(1) <= x_upper )
  call check(error, mcout%pars(2) >= y_lower  .and. mcout%pars(2) <= y_upper )
  ! remember point found by the first round
  startingpars = mcout%pars
  ! Now put the previously found point mcout$pars as starting point for another run with
  ! ll_bounded potential
  ! Keep McOUT, but reset step counter and statistics
  mcopt%nout = 0  ! run for zero steps
  ! pass optional restart=.true. argument to not generate new random starting point
  mcopt%restart = .true.
  call run_mcmc(ll_bounded, pi_xy, mcopt, mcout, seed=seed)
  ! Expect previously found point has been passed to new mcmc run
  call check(error, mcout%pars(1) == startingpars(1) )
  call check(error, mcout%pars(2) == startingpars(2) )
  ! expect the loglikelihood is never-Infinity
  call check(error, mcout%ll > -999999 )
  ! Run for a few more rounds
  mcopt%nout = 100000
  call run_mcmc(ll_bounded, pi_xy, mcopt, mcout, seed=seed)
  ! Expect we continue to be in the "allowed" region
  call check(error, mcout%ll > -999999 )
end subroutine test_mcmc_two_phase

subroutine test_mcmc_two_phase_parallel(error)
  !! Reproducing find_edc+main MCMC cardamom run.
  !! Repoducing parallel structure in main:  find_edc is called in top level omp parallel do,
  !! then run_parallel_mcmc is called.
  !! First find starting points point in ll = 0 region of ll_step potential, then
  !! start at these points with ll_bounded potential, which returns
  !! -Infinity if started at points outside of this region
  implicit none
  integer, parameter:: nchains = 4
  type(error_type), allocatable, intent(out):: error
  type(mcmc_output), dimension(:), allocatable:: mcout_list
  type(mcmc_options):: mcopt
  integer:: maxsteps
  double precision, dimension(nchains, 2):: startingpars
  integer:: i
  integer:: seed
  seed = random_int()
  call omp_set_num_threads(nchains)
  allocate(mcout_list(nchains))
  call init_pi()
  maxsteps = 1000000
  ! shared MCopt struct
  mcopt%nout = maxsteps
  mcopt%P_target = 0.0d0

  !$omp parallel do default(shared)
  do i = 1, nchains
  call run_mcmc(ll_step, pi_xy, mcopt, mcout_list(i), chainid = i, seed=seed+i)
  end do
  !$omp end parallel do

  do i = 1, nchains
  ! Expect we have optimized the values to the correct ranges
  call check(error, mcout_list(i)%pars(1) >= x_lower .and. mcout_list(i)%pars(1) <= x_upper )
  call check(error, mcout_list(i)%pars(2) >= y_lower .and. mcout_list(i)%pars(2) <= y_upper )
  ! remember points found by the first round
  startingpars(i, :) = mcout_list(i)%pars
  end do
  call check(error, mcout_list(1)%pars(2) /= mcout_list(2)%pars(2))

  ! Now put the previously found point mcout$pars as starting point for another run with
  ! ll_bounded potential
  ! Keep McOUT, but reset step counter and statistics
  mcopt%nout = 0  ! run for zero steps
  ! pass optional restart=.true. argument to not generate new random starting points
  mcopt%restart = .true.
  call run_parallel_mcmc(ll_bounded, pi_xy, mcopt, mcout_list, seed=seed)

  ! loop of checks
  do i = 1, nchains
  ! Expect previously found point has been passed to new mcmc run
  call check(error, mcout_list(i)%pars(1) == startingpars(i, 1) )
  call check(error, mcout_list(i)%pars(2) == startingpars(i, 2) )
  ! expect the loglikelihood is never-Infinity
  call check(error, mcout_list(i)%ll > -999999 )
  end do
  ! Run for a few more rounds
  mcopt%nout = 100000
  call run_parallel_mcmc(ll_bounded, pi_xy, mcopt, mcout_list, seed=seed)
  ! Expect we continue to be in the "allowed" region
  do i = 1, nchains
  write(*,*) mcout_list(i)%pars
  call check(error, mcout_list(i)%pars(1) >= x_lower .and. mcout_list(i)%pars(1) <= x_upper )
  call check(error, mcout_list(i)%pars(2) >= y_lower .and. mcout_list(i)%pars(2) <= y_upper )
  call check(error, mcout_list(i)%ll > -999999 )
  end do
end subroutine test_mcmc_two_phase_parallel

end module test_APMCMC
