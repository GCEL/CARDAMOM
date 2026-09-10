module test_random
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use random_uniform
  use samplers_math, only: random_int
  implicit none
  private

  public:: collect_randomtests

  integer, parameter:: dp = kind(0.0d0)

contains

!> Collect all exported unit tests
subroutine collect_randomtests(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out):: testsuite(:)

  testsuite = [ &
    new_unittest("test random_int(3)", test_random_int ),  &
    new_unittest("test random_int()", test_random_int_noarg ),  &
    new_unittest("fill random_uniform", test_fill_random_uniform),  &
    new_unittest("initialize random uniform", test_initialize),  &
    new_unittest("get_random_uniform", test_get_random_uniform),  &
    new_unittest("next_random_uniform", test_next_random_uniform), &
    new_unittest("set seed", test_set_seed), &
    new_unittest("seed consistency", test_set_seed_consistency), &
    new_unittest("seed independence", test_set_threadsafe_seed), &
    new_unittest("refill independence", test_threadsafe_refill) &
    ]

end subroutine collect_randomtests

subroutine test_random_int(error)
  !! Check behavior of random_int(N) on a small range 1 to N
  implicit none
  type(error_type), allocatable, intent(out):: error
  integer, parameter :: n_samples = 100
  integer :: upper_bound = 3
  integer, dimension(n_samples) :: results
  integer :: i
  do i=1, n_samples
    results(i) = random_int(upper_bound)
  end do
  ! expect integers in range 1 to 3 inclusive.
  ! check 1 occurs in results
  call check(error, any(results == 1 ) )
  ! check 3 occurs in results
  call check(error, any(results == upper_bound ) )
  ! check only 1,2,3 occur in results
  call check(error, all(results > 0 .and. results <= upper_bound  ) )
end subroutine test_random_int

subroutine test_random_int_noarg(error)
  !! Check behavior of random_int() with no arg, ints 1 to HUGE
  implicit none
  type(error_type), allocatable, intent(out):: error
  integer, parameter :: n_samples = 1000
  integer :: expected_upper_bound = HUGE(1)
  integer, dimension(n_samples) :: results
  double precision :: avg
  integer :: i
  do i=1, n_samples
    results(i) = random_int()
  end do
  ! expect integers in range 1 to HUGE inclusive.
  ! check all positive
  call check(error, all(results > 0) )
  ! Check values are approx evenly distributed in quarters of the range
  call check(error, count(results <= (expected_upper_bound * 0.25d0) ) >= n_samples * 0.2  )
  call check(error, count(results > (expected_upper_bound * 0.25d0) .and. (results <= expected_upper_bound * 0.5d0) ) >= n_samples * 0.2  )
  call check(error, count(results > (expected_upper_bound * 0.5d0) .and. (results <= expected_upper_bound * 0.75d0) ) >= n_samples * 0.2  )
  call check(error, count(results > expected_upper_bound * 0.75d0 ) >= n_samples * 0.2  )
  ! check average is in the ballpark of HUGE/2
  avg =sum( results / dble(n_samples))  ! scaling happens before summing to avoid integer overflow
  call check(error, avg >= expected_upper_bound * 0.4d0 .and. avg <= expected_upper_bound * 0.6d0 )
end subroutine test_random_int_noarg

subroutine test_fill_random_uniform(error)
  !! Run fill_random_uniform and check that the array
  !! on the object contains values in range 0 to 1
  implicit none
  type(error_type), allocatable, intent(out):: error
  type(UNIF_VECTOR):: random_uniform
  integer:: n
  double precision, dimension(:), allocatable:: arr
  integer:: seed
  seed = random_int()
  call random_uniform%initialize_random(seed)
  n = random_uniform%length
  allocate(arr(n))
  arr = 0d0
  call fill_random_uniform(arr, n, random_uniform%ranx)
  !write(*,*) arr
  call check(error, arr(1) > 0d0 .and. arr(1) <= 1.0 )
  call check(error, arr(5) > 0d0 .and. arr(5) <= 1.0 )
  call check(error, arr(n) > 0d0 .and. arr(n) <= 1.0 )
  call check(error, arr(2) /= arr(1) )
end subroutine test_fill_random_uniform

subroutine test_initialize(error)
  !! Call type bound initialize() of UNIF_VECTOR
  implicit none
  type(error_type), allocatable, intent(out):: error
  type(UNIF_VECTOR):: random_uniform
  integer:: seed
  seed = random_int()
  call random_uniform%initialize_random(seed)
  call check(error, random_uniform%index, 1 )
  call check(error, allocated(random_uniform%u))
  call check(error, random_uniform%u(2) > 0d0 .and. random_uniform%u(2) <= 1.0 )
  call check(error, random_uniform%u(2) /= random_uniform%u(1) )
  call check(error, random_uniform%u(random_uniform%length) > 0d0 .and. random_uniform%u(random_uniform%length) <= 1.0 )
end subroutine test_initialize

subroutine test_get_random_uniform(error)
  !! Get array of random values from UNIF_VECTOR object
  implicit none
  type(error_type), allocatable, intent(out):: error
  type(UNIF_VECTOR):: random_uniform
  integer:: n
  double precision, dimension(:), allocatable:: x
  integer:: seed
  seed = random_int()
  call random_uniform%initialize_random(seed)
  n = 1
  x = random_uniform%get_random_uniform(n)
  call check(error, x(1) > 0d0 .and. x(1) <= 1.0 )
  n = 12
  x = random_uniform%get_random_uniform(n)
  call check(error, x(1) > 0d0 .and. x(1) <= 1.0 )
  call check(error, x(n) > 0d0 .and. x(n) <= 1.0 )
  !n = random_uniform%length+1
  !x = random_uniform%get_random_uniform(n)  ! TODO not handled,  error/crash expected
end subroutine test_get_random_uniform

subroutine test_next_random_uniform(error)
  !! Get single value from UNIF_VECTOR object
  implicit none
  type(error_type), allocatable, intent(out):: error
  double precision  :: x
  integer ::  index
  type(UNIF_VECTOR):: random_uniform
  integer:: seed
  seed = random_int()
  call random_uniform%initialize_random(seed)
  index = random_uniform%index
  x = random_uniform%next_random_uniform()
  call check(error, x > 0d0 .and. x <= 1.0 )
  call check(error, random_uniform%index, index+1)
end subroutine test_next_random_uniform

subroutine test_set_seed(error)
  !! initialize UNIF_VECTOR object with a seed
  implicit none
  type(error_type), allocatable, intent(out):: error
  type(UNIF_VECTOR):: random_uniform
  integer:: n
  double precision, dimension(:), allocatable:: x
  integer:: seed = 155
  call random_uniform%initialize_random(seed)
  call check(error, random_uniform%seed, seed)
  n = 110
  x = random_uniform%get_random_uniform(n)
  call check(error, x(1) > 0d0 .and. x(1) <= 1.0 )
  call check(error, x(103) > 0d0 .and. x(103) <= 1.0 )
end subroutine test_set_seed

subroutine test_set_seed_consistency(error)
  !! Check that same explicit seed leads to same random values
  implicit none
  type(error_type), allocatable, intent(out):: error
  type(UNIF_VECTOR):: random_uniform
  integer:: seed1 = 155
  integer:: seed2 = 80
  double precision:: value1
  double precision, dimension(11):: values
  double precision:: value12
  call random_uniform%initialize_random(seed1)
  value1 = random_uniform%next_random_uniform()
  values = random_uniform%get_random_uniform(11)
  value12 = values(11)
  call random_uniform%initialize_random(seed2)
  ! check values are different with different seed
  call check(error, value1 /= random_uniform%next_random_uniform())
  values = random_uniform%get_random_uniform(11)
  call check(error, value12 /= values(11))
  call random_uniform%initialize_random(seed1)
  ! check values are the same as previously from the same seed 155
  call check(error, value1, random_uniform%next_random_uniform())
  values = random_uniform%get_random_uniform(11)
  call check(error, value12, values(11))
end subroutine test_set_seed_consistency

subroutine test_set_threadsafe_seed(error)
  !! Initialize multiple UNIF_VECTOR objects with different explicit seeds
  !! and check that they received different seeds, values.
  type(error_type), allocatable, intent(out):: error
  type(UNIF_VECTOR):: random_uniform1
  type(UNIF_VECTOR):: random_uniform2
  integer:: seed1 = 155
  integer:: seed2 = 298
  call random_uniform1%initialize_random(seed1)
  call random_uniform2%initialize_random(seed2)
  call check(error, random_uniform1%next_random_uniform() /= random_uniform2%next_random_uniform())
  call check(error, random_uniform1%next_random_uniform() /= random_uniform2%next_random_uniform())
  call check(error, random_uniform1%seed /= random_uniform2%seed)
end subroutine test_set_threadsafe_seed


subroutine test_threadsafe_refill(error)
  !! Initialize multiple UNIF_VECTOR objects with different explicit seeds
  !! and check that they received different seeds, values.
  type(error_type), allocatable, intent(out):: error
  type(UNIF_VECTOR):: random_uniform1
  type(UNIF_VECTOR):: random_uniform2
  type(UNIF_VECTOR):: random_uniform3
  double precision:: value1, value2, value3
  integer:: seed1
  integer:: seed2
  seed1 = random_int()
  seed2 = random_int()

  call random_uniform1%initialize_random(seed1)

  ! possible, this sets some shared values like ranx to 'wrong' values
  ! matching seed 2
  call random_uniform2%initialize_random(seed2)

  ! simulate internally-triggered refill-possibly from ranx influenced by seed2
  call fill_random_uniform(random_uniform1%u, random_uniform1%length, random_uniform1%ranx)

  random_uniform1%index = 1
  value1 = random_uniform1%next_random_uniform()
  value2 = random_uniform1%next_random_uniform()
  value3 = random_uniform1%next_random_uniform()

  !vs same sequence from seed 1 without possible seed2 contamination
  call random_uniform3%initialize_random(seed1)
  call fill_random_uniform(random_uniform3%u, random_uniform3%length, random_uniform3%ranx)

  random_uniform3%index = 1
  call check(error, value1, random_uniform3%next_random_uniform())
  call check(error, value2, random_uniform3%next_random_uniform())
  call check(error, value3, random_uniform3%next_random_uniform())
  call check(error, random_uniform1%seed, random_uniform3%seed)
  call check(error, random_uniform1%seed /= random_uniform2%seed)
end subroutine test_threadsafe_refill

subroutine test_random_multivariate(error)
  !! TODO
  implicit none
  type(error_type), allocatable, intent(out):: error
  ! make covariance matrix
  ! make npars
  ! make vector mu (mean) to center the step on
  ! get vector rn
  !call random_multivariate(npars, 1, covariance_matrix, mu, rn)
end subroutine test_random_multivariate


subroutine test_random_normal(error)
  !! TODO
  implicit none
  type(error_type), allocatable, intent(out):: error
  ! make covariance matrix
  ! make npars
  ! make vector mu (mean) to center the step on
  ! get output value result
  !call random_normal(result)
end subroutine test_random_normal


end module test_random
