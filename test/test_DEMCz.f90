Module test_DEMCz
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use DEMCz
   use test_functions
   implicit none

   public:: collect_DEMCztests

   integer, parameter:: dp = kind(0.0d0)

contains

!> Collect all exported unit tests
   subroutine collect_DEMCztests(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out):: testsuite(:)

      testsuite = [ &
                  new_unittest("Options_type", test_MCO), &
                  new_unittest("random_int", test_random_int), &
                  new_unittest("DEMcz runs", test_DEMcz_runs), &
                  new_unittest("DEMcz runs with omp threads", test_DEMcz_runs_enforce_omp) &
                  ]

   end subroutine collect_DEMCztests

   subroutine test_MCO(error)
      type(error_type), allocatable, intent(out):: error
      !Having imported DEMCz_module, we should have a type(DEMCzOpt)
      ! and access a module-level object holding default options
      type(DEMCzOPT):: options
   end subroutine test_MCO

! TODO move to common
   subroutine test_random_int(error)
      use DEMCz, only: random_int
      implicit none
      type(error_type), allocatable, intent(out):: error
      integer:: r
      r = random_int(2)
  !! should choose from (1, 2) (inclusive)
      call check(error, r == 1 .or. r == 2)
   end subroutine test_random_int

   subroutine test_DEMCz_runs(error)
      use DEMCz, only: run_DEMCz, PARINFO, random_int
      implicit none
      type(error_type), allocatable, intent(out):: error
      ! test the main DEMCz function just runs when given a function

      type(DEMCzOPT):: options
     !! new DEMCZ options struct with default values
      type(MCMC_OUTPUT), dimension(:), allocatable:: DEMCzOUT
     !! new (blank) struct to write results to

      integer :: seed

      ! PI: use the PI_xy struct from test_functions quadratic potential
      call init_pi()

      seed=random_int(300000)

      ! with MAXITER < nadapt
      options%nout = 10
      call run_DEMCz(ll_normal, PI_xy, options, DEMCzOUT, seed=seed)

      ! with MAXITER = N*nadapt
      options%nout = 50
      options%nadapt = 10
      call run_DEMCz(ll_normal, PI_xy, options, DEMCzOUT, seed=seed)

      ! with MAXITER /= N*nadapt
      options%nout = 53
      options%nadapt = 10
      call run_DEMCz(ll_normal, PI_xy, options, DEMCzOUT, seed=seed)

   end subroutine test_DEMCz_runs

   subroutine test_DEMCz_runs_enforce_omp(error)
      use DEMCz, only: run_DEMCz, PARINFO, random_int
      implicit none
      type(error_type), allocatable, intent(out):: error
      ! test the main DEMCz function just runs when given a function

      type(DEMCzOPT):: options
     !! new DEMCZ options struct with default values
      type(MCMC_OUTPUT), dimension(:), allocatable:: DEMCzOUT
     !! new (blank) struct to write results to
      integer:: nchains
      integer :: seed
      nchains = 4

      ! PI: use the PI_xy struct from test_functions quadratic potential
      call init_pi()

      seed=random_int(300000)

      call omp_set_num_threads(nchains)
      options%nout = 10

      call run_DEMCz(ll_normal, PI_xy, options, DEMCzOUT, seed=seed)

   end subroutine test_DEMCz_runs_enforce_omp

end module test_DEMCz
