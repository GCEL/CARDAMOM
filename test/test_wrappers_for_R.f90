module test_wrappers_for_R
  !! Helper functions for R scripts
  !! if they don't run from fortran they definitely won't run from R!
  !! Test here to see info on segfaults
  !! Behavior from R is hard to debug: fails intermittently and no info
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use model_shared, only :PI
  use cardamom_structures, only : DATAin
  use iso_c_binding
  use cardamom_Rinterfaces
  implicit none
  private

  public:: collect_test_wrappers_for_R

contains

subroutine collect_test_wrappers_for_R(testsuite)
  type(unittest_type), allocatable, intent(out):: testsuite(:)

  testsuite = [ &
    new_unittest("stresstest setup", test_initialize_stresstest_circle), &
    new_unittest("stresstest initialize twice", test_repeat_initialize_stresstest), &
    new_unittest("model setup from data file", test_initialize_model), &
    new_unittest("model initialize twice", test_repeat_initialize_model), &
    new_unittest("Report npars for R", test_get_model_npars), &
    new_unittest("Report PI%parmax for R", test_get_model_parmax), &
    new_unittest("Report PI%parmin for R", test_get_model_parmin), &
    new_unittest("Report DATAin%parini for R", test_get_initial_pars) &
    ]

end subroutine collect_test_wrappers_for_R

subroutine test_initialize_stresstest_circle(error)
  !! this should be able to run first in a program, set up and allocate
  !! all objects needed in stresstest (PI, DATAin, MCOUT, ...)
  use model_shared, only: PI
  use cardamom_structures, only: DATAin
  type(error_type), allocatable, intent(out):: error
  call initialize_stresstest_circle()
  !! expect PI is initialized to expected values
  call check(error, PI%parmin(1), -10d0)
  call check(error, PI%parmin(10), 1d0)
  call check(error, PI%parmax(3), 300d0)
  call check(error, PI%parmax(9), 300d0)
  ! datain was set, and "Circle" routine was selected
  call check(error, DATAin%noobs, 9)
end subroutine test_initialize_stresstest_circle


subroutine test_repeat_initialize_stresstest(error)
  !! check we don't get crashes when called initialize_stresstest_circle twice in a R session
  type(error_type), allocatable, intent(out):: error
  call initialize_stresstest_circle()
  call initialize_stresstest_circle()
end subroutine test_repeat_initialize_stresstest

subroutine test_initialize_model(error)
  type(error_type), allocatable, intent(out):: error
  !call initialize_cardamom("../../test/data/UK_baseline_sites_AliceHolt.bin")  ! assuming ctest is run from CARDAMOM/build/
  call initialize_cardamom()
    ! path relative to test exectuable at CARDAMOM/build/test/CARADAMOM-tester

  ! check we now have allocated and filled PI, DATAin
  ! from model file, input data
  call check(error, PI%parmin(1), -10d0)
  call check(error, PI%parmin(10), 1d0)
  call check(error, PI%parmax(3), 300d0)
  call check(error, PI%parmax(9), 300d0)
  call check(error, DATAin%noobs, 9)
end subroutine test_initialize_model

subroutine test_repeat_initialize_model(error)
  !! check this won't cause unexplained R crashes when called initialize_model() twice in a R session
  type(error_type), allocatable, intent(out):: error
  !call initialize_cardamom("../../test/data/UK_baseline_sites_AliceHolt.bin")
  call initialize_cardamom()
  !call initialize_cardamom("../../test/data/UK_baseline_sites_AliceHolt.bin")
  call initialize_cardamom()
end subroutine test_repeat_initialize_model

subroutine test_get_model_npars(error)
  !! test fortran side of function that should report compiled model's PI%npars to R,
  !! when its C name is called from R and with an R integer to write to
  type(error_type), allocatable, intent(out):: error
  integer:: npars
  !call initialize_cardamom("../../test/data/UK_baseline_sites_AliceHolt.bin")
  call initialize_cardamom()
  call get_npars(npars)
  call check(error, npars, 32)
end subroutine test_get_model_npars

subroutine test_get_model_parmin(error)
  !! function that should write model PI%parmin to real c_double intent(inout) array,
  !! getter, for C/R compatility.
  type(error_type), allocatable, intent(out):: error
  integer:: npars
  real(c_double), dimension(:), allocatable:: parmin
  !call initialize_cardamom("../../test/data/UK_baseline_sites_AliceHolt.bin")
  call initialize_cardamom()
  call get_npars(npars)
  allocate(parmin(npars))
  call get_parmin(npars, parmin)
  call check(error, size(parmin), npars)
  call check(error, parmin(1) == 0_c_double)
  call check(error, parmin(32) == 0_c_double)
end subroutine test_get_model_parmin

subroutine test_get_model_parmax(error)
  !! function that should write model PI%parmin to real c_double intent(inout) array,
  !! getter, for C/R compatility.
  type(error_type), allocatable, intent(out):: error
  integer:: npars
  real(c_double), dimension(:), allocatable:: parmax
  !call initialize_cardamom("../../test/data/UK_baseline_sites_AliceHolt.bin")
  call initialize_cardamom()
  call get_npars(npars)
  allocate(parmax(npars))
  call get_parmax(npars, parmax)
  call check(error, size(parmax), npars)
  call check(error, parmax(1)== 0_c_double)
  call check(error, parmax(32)== 0_c_double)
end subroutine test_get_model_parmax

subroutine test_get_initial_pars(error)
  !! function to report DATAin%parini to R, after initializtion from data file
  type(error_type), allocatable, intent(out):: error
  integer:: npars
  real(c_double), dimension(:), allocatable:: parmin, parmax, parini
  !call initialize_cardamom("../../test/data/UK_baseline_sites_AliceHolt.bin")
  call initialize_cardamom()
  call get_npars(npars)
  allocate(parmin(npars), parmax(npars), parini(npars))
  call get_parmin(npars, parmin)
  call get_parmax(npars, parmax)
  !call get_initial_pars(npars, parini)
  call check(error, size(parini), npars)
  call check(error, parini(2) >= parmin(2))
  call check(error, parini(2) <= parmax(2))
  call check(error, parini(19) >= parmin(19))
  call check(error, parini(19) <= parmax(19))
end subroutine test_get_initial_pars

subroutine test_get_model_loglikelihood(error)
  type(error_type), allocatable, intent(out):: error
end subroutine test_get_model_loglikelihood

end module test_wrappers_for_R
