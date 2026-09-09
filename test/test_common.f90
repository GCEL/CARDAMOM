module test_common
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use random_uniform
   use test_functions
   use samplers_shared
   use samplers_math, only: random_int
   implicit none
   private

   public:: collect_commontests

   integer, parameter:: dp = kind(0.0d0)

contains

!> Collect all exported unit tests
   subroutine collect_commontests(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out):: testsuite(:)

      testsuite = [ &
                  new_unittest("infini", test_is_infinity), &
                  new_unittest("test init pars random", test_init_pars_random), &
                  new_unittest("metropolis_choice", test_metropolis_choice), &
                  new_unittest("metropolis_increase", test_metropolis_increase), &
                  new_unittest("metropolis_stochastic", test_metropolis_stochastic) &
                  ]

   end subroutine collect_commontests

   subroutine test_is_infinity(error)
      type(error_type), allocatable, intent(out):: error
      double precision:: P, log_P
      P = 0.0_dp
      log_P = log(P)
      call init_infinity()
      call check(error, log_P == neg_inf, .true.)
      P = 0.01_dp
      log_P = log(P)
      call check(error, log_P == neg_inf, .false.)
   end subroutine test_is_infinity

   subroutine test_init_pars_random(error)
      type(error_type), allocatable, intent(out):: error
      logical, dimension(:), allocatable:: fix_pars_flag
      double precision, dimension(:), allocatable:: pars0
      type(UNIF_VECTOR):: random_uniform
      integer:: seed
      seed = random_int()  ! this test with a different seed each time
      call random_uniform%initialize_random(seed)

      if (.not. allocated(pars0)) allocate (pars0(PI_xy%npars))
      ! with no fix_pars_flag
      call init_pars_random(PI_xy, pars0, uniform_random_vector=random_uniform)
      if (.not. allocated(fix_pars_flag)) allocate (fix_pars_flag(PI_xy%npars))
      ! with fix_pars_flag all false
      fix_pars_flag = .false.
      call init_pars_random(PI_xy, pars0, fix_pars_flag, random_uniform)
      ! With fix_pars_flag all true
      fix_pars_flag = .true.
      call init_pars_random(PI_xy, pars0, fix_pars_flag, random_uniform)
   end subroutine test_init_pars_random

   subroutine test_metropolis_choice(error)
      !> Metropolis chioce function takes two log(!) likelihoods
      !> and returns logical
      !> The first argument is the new/proposed loglikelihood
      use samplers_shared, only: metropolis_choice
      implicit none
      type(error_type), allocatable, intent(out):: error
      ! certain acceptance of state with probability 1 vs 0
      call check(error, metropolis_choice(log(1.0_dp), log(1e-15_dp)), .true.)
      ! certain rejection of state with probability 0 vs 1
      call check(error, metropolis_choice(log(1e-15_dp), log(1.0_dp)), .false.)

   end subroutine test_metropolis_choice

   subroutine test_metropolis_increase(error)
      !> We are maximizing ll
      !> unconditional acceptance if new ll is bigger
      use samplers_shared, only: metropolis_choice
      implicit none
      type(error_type), allocatable, intent(out):: error
      ! certain accept if new (first argument) ll is bigger
      call check(error, metropolis_choice(1.1_dp, 1.0_dp), .true.)
   end subroutine test_metropolis_increase

   subroutine test_metropolis_stochastic(error)
      !> Metropolis chioce function takes two log(!) likelihoods
      !> and returns logical
      !> The first argument is the new/proposed loglikelihood
      use samplers_shared, only: metropolis_choice
      implicit none
      type(error_type), allocatable, intent(out):: error
      double precision:: new_loglikelihood, old_loglikelihood, accept_ratio
      integer:: i, N, accept_count
      ! Something with a likelihood l1 = 1/2 l2 should be acceped
      ! 50% of the time.
      new_loglikelihood = log(.3)
      old_loglikelihood = log(.6)
      ! get acceptance N times  - expect about 50% true
      N = 500
      accept_count = 0
      do i = 1, N
         if (metropolis_choice(new_loglikelihood, old_loglikelihood)) then
            accept_count = accept_count + 1
         end if
      end do
      accept_ratio = accept_count/real(N)
      write (*, *) accept_count
      write (*, *) accept_ratio
      call check(error, accept_ratio > .4 .and. accept_ratio < .6)
   end subroutine test_metropolis_stochastic

end module test_common
