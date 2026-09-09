module test_random
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use CARBON_MODEL_MOD
  implicit none
  private

  public:: collect_modeltests

  integer, parameter:: dp = kind(0.0d0)

contains

!> Collect all exported unit tests
subroutine collect_modeltests(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out):: testsuite(:)

  testsuite = [ &
    new_unittest("Same output for same input", test_consistency) &
    ]

end subroutine collect_modeltests

subroutine test_consistency(error)
  !! As in _likelihood.f90 model_sanity_check  :
  !! Calling carbon_model twice with same input,
  !! we expect same output

  ! variables for carbon_model intent(in)

  ! arrays for carbon_model intent(out)

  ! second set for second run
end subroutine test_consistency


end module test_random
