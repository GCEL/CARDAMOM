module test_math
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use samplers_math
  implicit none
  private

  public:: collect_mathtests, approx

  integer, parameter:: dp = kind(0.0d0)

  ! some example arguments for conversion
  integer, parameter:: npars = 5
  double precision, dimension(npars):: parmax = [1d0, 0d0, 1.20d0, -13.5d0, 0.00000043d0]
  double precision, dimension(npars):: parmin = [0.00d0, -0.3d0, -1.20d0,  -190.5d0, 0.00000042d0]
  double precision, dimension(npars):: paradj = [1.0d0, 1.3d0, 2.20d0,  191.5d0, 0.0d0] ! abs(parmin1 if nonpositive
  double precision, dimension(npars):: pars1 = [1.0d0, -.01d0, 0.0d0, -100d0, 0.000000421d0]
  double precision, dimension(npars):: pars_norm1 = [0.0d0, 0.99d0, 0.32d0, 1.0d0, 0.000001d0]
contains

!> Collect all exported unit tests
subroutine collect_mathtests(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out):: testsuite(:)

  testsuite = [ &
    new_unittest("approx", test_approx),  &
    new_unittest("par2nor", test_par2nor),  &
    new_unittest("nor2par", test_nor2par),  &
    new_unittest("log_par2nor", test_log_par2nor),  &
    new_unittest("log_nor2par", test_log_nor2par),  &
    new_unittest("running_avg", test_running_avg_scalar),  &
    new_unittest("increment covariance matrix", test_running_covariance),  &
    new_unittest("running_covariance_of_stationary", test_zero_variance)  &
    ]

end subroutine collect_mathtests

! testing util approx equal
function approx(a, b, rel_tol) result(eq)
  implicit none
  double precision, intent(in):: a, b
  double precision, intent(in), optional:: rel_tol
  double precision:: rel_tol_default = 1.0e-8
  double precision:: rel_tol_
  logical:: eq
  if (.not. present(rel_tol)) then
    rel_tol_ = rel_tol_default
  else
    rel_tol_ = rel_tol
  endif
  ! check close to zero
  if (abs(a) <= epsilon(0d0)) then
    eq = (abs(b) <= epsilon(0d0))
  else if (abs(b) <= epsilon(0d0)) then
    eq = (abs(a) <= epsilon(0d0))
  else
  ! else check relative difference
    eq = (abs(a-b) <= rel_tol_*abs(a))  .or.   (abs(a-b) <= rel_tol*abs(b))
  endif
end function approx


subroutine test_approx(error)
  !! test my 'approx equals' testing utility
  implicit none
  type(error_type), allocatable, intent(out):: error
  call check(error, approx(1d0, 0d0), .false.)
  call check(error, approx(0d0, -1d0), .false.)
  call check(error, approx(0d0, 0d0), .true.)
  call check(error, approx(27d0/9d0, 3d0), .true.)
end subroutine test_approx

subroutine test_par2nor(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  double precision, dimension(npars):: pars_norm, pars2
  integer:: i
  pars_norm = par2nor( pars1, parmin, parmax)
  ! expect values between 0 and 1
  call check(error, all(pars_norm >= 0d0) .and. all(pars_norm <= 1d0))
  ! try converting back
  pars2 = nor2par(pars_norm, parmin, parmax)
  ! expect approximately equal to the original
  call check(error, all ([(approx(pars1(i), pars2(i)), i = 1, npars)] ))
end subroutine test_par2nor


subroutine test_nor2par(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  double precision, dimension(npars):: pars_real
  integer:: i
  pars_real = nor2par(pars_norm1, parmin, parmax)
  ! expect values in bounds parmin, parmax
  call check(error, all(pars_real >= parmin) .and. all(pars_real <= parmax))
end subroutine test_nor2par


subroutine test_log_par2nor(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  double precision, dimension(npars):: pars_norm, pars2
  integer:: i
  pars_norm = log_par2nor(pars1, parmin, parmax, paradj)
  ! expect values between 0 and 1
  call check(error, all(pars_norm >= 0d0) .and. all(pars_norm <= 1d0))
  ! try converting back
  pars2 = log_nor2par(pars_norm, parmin, parmax, paradj)
  ! expect approximately equal to the original
  call check(error, all ([(approx(pars1(i), pars2(i)), i = 1, npars)] ))
end subroutine test_log_par2nor


subroutine test_log_nor2par(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  double precision, dimension(npars):: pars_real
  pars_real = log_nor2par(pars_norm1, parmin, parmax, paradj)
  ! expect values in bounds parmin, parmax
  call check(error, all(pars_real >= parmin) .and. all(pars_real <= parmax))
end subroutine test_log_nor2par


subroutine test_running_avg_scalar(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  double precision, dimension(100):: series
  double precision, dimension(1, 1):: covariance   ! irrelevant covariance matrix
  double precision, dimension(1):: mean, mean75, mean85  ! correct answers
  double precision, dimension(1):: running_mean85  ! answers from increment_covariance_matrix update to mean
  integer:: cur
  covariance = 1d0
  ! ransom values 0..1
  call random_number(series)
  ! insert some higher values towards the end
  series(80:90) = series(80:90)* 5.0d0
  series(90:100) = series(80:90)* 10.0d0
  ! true average of the numbers
  mean(1) = sum(series)/size(series)
  ! mean up to 75
  mean75(1) = sum(series(1:75))/size(series(1:75))
  ! mean up to 85
  mean85(1) = sum(series(1:85))/size(series(1:85))
  ! try getting mean up to 85 from mean up to 75+running average fct
  running_mean85(1) =  mean75(1)
  cur = 75  ! is also incremented by the call
  call increment_covariance_matrix(series(76:85), running_mean85, 1, cur, 10, covariance)
  call check(error, cur, 85)
  call check(error, approx(mean85(1), running_mean85(1)))
end subroutine test_running_avg_scalar


subroutine test_running_covariance(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  double precision, dimension(2, 100):: series  ! random series
  double precision, dimension(2):: mean, mean10
  double precision, dimension(2, 2):: covariance, covariance_10   ! 2 x 2 covariance matrix
  double precision, dimension(2, 2):: covariance_running
  integer:: cur
  covariance = 0
  call random_number(series)

  ! true covariance of whole series
  call covariance_matrix(series, mean, 2, 100, covariance)
  ! covariance up to 10
  call covariance_matrix(series(1:2, 1:10), mean, 2, 10, covariance_10)
  ! mean up to 10
  mean10(1) = sum(series(1, 1:10))/size(series(1, 1:10))
  mean10(2) = sum(series(2, 1:10))/size(series(2, 1:10))
  covariance_running = covariance_10
  ! try getting covariance by incrementing
  cur = 10
  call increment_covariance_matrix(series(1:2, 11:100), mean10, 2, cur, 90, covariance_running)
  ! expect values in
  write(*,*) covariance
  write(*,*) covariance_10
  write(*,*) covariance_running
  call check(error, approx(covariance(1, 1), covariance_running(1, 1)))
  call check(error, approx(covariance(2, 1), covariance_running(2, 1)))
end subroutine test_running_covariance

subroutine test_zero_variance(error)
  implicit none
  type(error_type), allocatable, intent(out):: error
  double precision, dimension(2, 500000):: series  ! a very static series
  double precision, dimension(2):: mean  ! irrelevant mean
  double precision, dimension(2, 2):: covariance, covariance_10   ! 2 x 2 covariance matrix
  double precision, dimension(2, 2):: covariance_running
  integer:: cur
  covariance = 0
  series = 0.5d0  ! then constant
  series(1:2, 2) = 0.1d0  ! some initial movement
  series(1:2, 4) = 0.9d0
  series(1:2, 1) = 0.1d0

  ! true covariance of whole series
  call covariance_matrix(series, mean, 2, 500000, covariance)
  ! covariance up to 10
  call covariance_matrix(series(1:2, 1:10), mean, 2, 10, covariance_10)

  ! try getting covariance by incrementing
  covariance_running = covariance_10
  cur = 10
  call increment_covariance_matrix(series(1:2, 11:500000), mean, 2, cur, 499990, covariance_running)
  ! expect values in
  write(*,*) covariance
  write(*,*) covariance_10
  write(*,*) covariance_running
  call check(error,  covariance(2, 2) < .000001d0)
  call check(error,  covariance_running(1, 1) < .000001d0)
  call check(error, approx(covariance(1, 1), covariance_running(1, 1)))
  call check(error, approx(covariance(2, 1), covariance_running(2, 1)))
end subroutine test_zero_variance


end module test_math
