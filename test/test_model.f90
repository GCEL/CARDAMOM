module test_model
  !!Testing the example model DALEC.A1.C1.D2.F2.H2.P1
  !! These tests are flaky !
  !! - they always print warnings 'Integer overflow when calculating the amount of memory to allocate'
  !!   Attempt to DEALLOCATE unallocated '%s' , but these don't seem to cause tests to fail
  !! - Line numbers of failing checks are sometimes not pointing to the true problem
  !! - They pass or fail for random reasons such as presence of a write(*,*) statement in
  !!    (not even run) model_sanity_check, or allocation of 'test' array .
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use test_functions
  use test_math, only: approx
  use random_uniform
  use samplers_math, only: random_int
  use model_shared, only: initialize_carbon_model, destroy_carbon_model
  use APMCMC
  use model_shared, only: PI
  use cardamom_io, only: initialize
  use cardamom_main_utils
  use OMP_LIB
  use ieee_arithmetic
  implicit none
  private

  public:: collect_modeltests

  character(len=*), parameter:: infile = "../../test/data/UK_baseline_sites_AliceHolt.bin"

  ! relative to ctest working directory, by default CARDAMOM/build/test/ , where add_test was called

contains

!> Collect all exported unit tests
subroutine collect_modeltests(testsuite)
  implicit none
  !> Collection of tests
  type(unittest_type), allocatable, intent(out):: testsuite(:)

  testsuite = [ &
    new_unittest("model_initialize", test_model_initialize), &
    new_unittest("carbon_model_sanity_check", test_model_sanity_check), &
    new_unittest("carbon_model_not_nan", test_carbon_model_not_nan), &
    new_unittest("model_repeat_evaluation", test_model_repeat_evaluation) &
    !new_unittest("model_group_repeat_evaluation", test_model_group_repeat_evaluation), &
    !new_unittest("model_group_repeat_evaluation_3", test_model_group_repeat_evaluation_3) &
    ]

end subroutine collect_modeltests


subroutine test_model_initialize(error)
  use CARBON_MODEL_MEMORY, only: mvs
  implicit none
  type(error_type), allocatable, intent(out):: error
  integer:: nchains
  type(MCMC_OUTPUT):: MCOUT
  nchains = 1
    call initialize(infile)
    call initialize_carbon_model(nchains)
    call initialize_stats(MCOUT, PI%npars)
    ! check initialized?
    call check(error, allocated(mVs(1)%rainfall_time))
    call destroy_carbon_model()
    call check(error, .not. allocated(mVs))
end subroutine test_model_initialize

subroutine test_model_loglikelihood(error)
  !! A number should come back as loglikelihood,
  !! expect not NaN and negative
  type(error_type), allocatable, intent(out):: error
end subroutine test_model_loglikelihood

subroutine test_model_sanity_check(error)
    !! check model repeatabilty via inbuilt "sanity check" function,
    !! should have same result as tests here
    use model_likelihood_module, only: model_sanity_check, sanity_check
    use cardamom_structures, only: DATAin
    use samplers_shared, only : init_pars_random
    type(error_type), allocatable, intent(out):: error
    double precision, dimension(:), allocatable:: PARS
    integer:: nchains, seed
    type(UNIF_VECTOR):: random_uniform
    nchains = 1
    call initialize(infile)
    call initialize_carbon_model(nchains)
    PARS = DATAin%parpriors(1:PI%npars)
    seed = random_int()
    call random_uniform%initialize_random(seed)
    call init_pars_random(PI, PARS, PI%fix_pars, random_uniform)
    call model_sanity_check(PARS, 1)
    call check(error, sanity_check )
    call destroy_carbon_model()
end subroutine test_model_sanity_check


subroutine test_carbon_model_not_nan(error)
  !! Test Model-internal subroutine carbon_model :
  !! POOLS, FLUXES, DIAGS arrays to compare to data should
  !! come back filled and not NaN
  use CARBON_MODEL_MOD, only: carbon_model
  use CARBON_MODEL_MEMORY, only: mvs
    use cardamom_structures, only: DATAin
    use samplers_shared, only : init_pars_random
    implicit none
  type(error_type), allocatable, intent(out):: error
  double precision, dimension(:), allocatable:: PARS
  !dimensions are not known before they're read from file at initialize()
    double precision, dimension(:,:), allocatable:: pools
    double precision, dimension(:,:), allocatable:: fluxes
    double precision, dimension(:,:), allocatable:: diags
    double precision, dimension(:,:), allocatable:: test
    integer:: nchains, seed
    type(UNIF_VECTOR):: random_uniform
    allocate(test(97, 49))
    nchains = 1
    call initialize(infile)
    allocate(pools(DATAin%nodays+1, DATAin%nopools))
    ! TODO here we get ' Integer overflow when calculating the amount of memory to allocate '
    allocate(fluxes(DATAin%nodays, DATAin%nofluxes))
    allocate(diags(DATAin%nodays, DATAin%nodiags))
    call initialize_carbon_model(nchains)
    PARS = DATAin%parpriors(1:PI%npars)
    seed = random_int()
    call random_uniform%initialize_random(seed)
    call init_pars_random(PI, PARS, PI%fix_pars, random_uniform)
    call carbon_model(1, DATAin%nodays, DATAin%MET, PARS, DATAin%deltat &
                     ,DATAin%nodays, DATAin%LAT &
                     ,fluxes, pools, diags,  DATAin%nopars &
                     ,DATAin%nomet, DATAin%nopools, DATAin%nofluxes  &
                     ,DATAin%nodiags, mVs(1))
    call check(error, .not. ieee_is_nan(sum(fluxes)) )
    call check(error, .not. ieee_is_nan(sum(pools)) )
    call check(error, .not. ieee_is_nan(sum(diags)) )
    call destroy_carbon_model()
end subroutine test_carbon_model_not_nan


subroutine test_model_repeat_evaluation(error)
  !! Test that repeat calls to carbon_model with same PARS
  !! gives same results.
  !! Similar to original cardamom "Sanity check" function.
  use CARBON_MODEL_MOD, only: carbon_model
  use CARBON_MODEL_MEMORY, only: mvs
    use cardamom_structures, only: DATAin
    use samplers_shared, only : init_pars_random
  implicit none
  type(error_type), allocatable, intent(out):: error
  integer:: nchains
    double precision, dimension(:), allocatable:: PARS
    double precision, dimension(:,:), allocatable:: pools1, pools2
    double precision, dimension(:,:), allocatable:: fluxes1, fluxes2
    double precision, dimension(:,:), allocatable:: diags1, diags2
    double precision:: pool_error, flux_error, diag_error
    double precision, dimension(:,:), allocatable:: test
    ! These are here because they're kept at this level in model_likelihood.f90 files,
    ! and they are there in DALEC models because this text insertion was the simplest way to edit 37 models
    integer:: seed
    type(UNIF_VECTOR):: random_uniform
    allocate(test(97, 49))
    nchains = 1
    call initialize(infile)
    allocate(pools1(DATAin%nodays+1, DATAin%nopools), pools2(DATAin%nodays+1, DATAin%nopools))
    allocate(fluxes1(DATAin%nodays, DATAin%nofluxes), fluxes2(DATAin%nodays, DATAin%nofluxes))
    allocate(diags1(DATAin%nodays, DATAin%nodiags), diags2(DATAin%nodays, DATAin%nodiags))
    call initialize_carbon_model(nchains)
    PARS = DATAin%parpriors(1:PI%npars)
    seed = random_int()  ! TODO record later
    call random_uniform%initialize_random(seed)
    call init_pars_random(PI, PARS, PI%fix_pars, random_uniform)
    call carbon_model(1, DATAin%nodays, DATAin%MET, PARS, DATAin%deltat &
                     ,DATAin%nodays, DATAin%LAT &
                     ,fluxes1, pools1, diags1,  DATAin%nopars &
                     ,DATAin%nomet, DATAin%nopools, DATAin%nofluxes  &
                     ,DATAin%nodiags, mVs(1))
    call check(error, .not. ieee_is_nan(sum(pools1)) )
    call check(error, .not. ieee_is_nan(sum(fluxes1)) )
    call check(error, .not. ieee_is_nan(sum(diags1)) )
    call carbon_model(1, DATAin%nodays, DATAin%MET, PARS, DATAin%deltat &
                     ,DATAin%nodays, DATAin%LAT &
                     ,fluxes2, pools2, diags2,  DATAin%nopars &
                     ,DATAin%nomet, DATAin%nopools, DATAin%nofluxes  &
                     ,DATAin%nodiags, mVs(1))
    call check(error, .not. ieee_is_nan(sum(fluxes2)) )
    call check(error, .not. ieee_is_nan(sum(pools2)) )
    call check(error, .not. ieee_is_nan(sum(diags2)) )
    flux_error = sum(abs(fluxes1-fluxes2))
    pool_error = sum(abs(pools1-pools2))
    diag_error = sum(abs(diags1-diags2))
    call check(error, flux_error < .0000000001)
    call check(error, pool_error < .0000000001)
    call check(error, diag_error < .0000000001)
    call destroy_carbon_model()
end subroutine test_model_repeat_evaluation

subroutine test_model_repeat_and_independent_evaluation(error)
  !! Similar to original cardamom "Sanity check" function
  implicit none
  type(error_type), allocatable, intent(out):: error
end subroutine test_model_repeat_and_independent_evaluation


subroutine test_model_group_repeat_evaluation(error)
  use CARBON_MODEL_MOD, only: carbon_model
  use CARBON_MODEL_MEMORY, only: mvs
    use cardamom_structures, only: DATAin
    use samplers_shared, only : init_pars_random
  !! Similar to original cardamom "Sanity check" function
  implicit none
  type(error_type), allocatable, intent(out):: error
  integer:: nchains
    double precision, dimension(PI%npars):: PARS
    double precision, dimension(:,:), allocatable:: pools1, pools2
    double precision, dimension(:,:), allocatable:: fluxes1, fluxes2
    double precision, dimension(:,:), allocatable:: diags1, diags2
    double precision:: pool_error, flux_error, diag_error
integer:: seed, i, clock
    type(UNIF_VECTOR):: random_uniform
    nchains = 1
    call initialize(infile)
    allocate(pools1((DATAin%nodays+1), DATAin%nopools), pools2((DATAin%nodays+1), DATAin%nopools))
    allocate(fluxes1(DATAin%nodays, DATAin%nofluxes), fluxes2(DATAin%nodays, DATAin%nofluxes))
    allocate(diags1(DATAin%nodays, DATAin%nodiags), diags2(DATAin%nodays, DATAin%nodiags))
    PARS = DATAin%parpriors(1:PI%npars)
    seed = random_int()
    call random_uniform%initialize_random(seed)
    call init_pars_random(PI, PARS, PI%fix_pars, random_uniform)
    call initialize_carbon_model(nchains)
    do i = 1, nchains
    call carbon_model(1, DATAin%nodays, DATAin%MET, PARS, DATAin%deltat &
                     ,DATAin%nodays, DATAin%LAT &
                     ,fluxes1, pools1, diags1,  DATAin%nopars &
                     ,DATAin%nomet, DATAin%nopools, DATAin%nofluxes  &
                     ,DATAin%nodiags, mVs(i))
    call check(error, .not. ieee_is_nan(sum(fluxes1)) )
    call check(error, .not. ieee_is_nan(sum(pools1)) )
    call check(error, .not. ieee_is_nan(sum(diags1)) )
    call carbon_model(1, DATAin%nodays, DATAin%MET, PARS, DATAin%deltat &
                     ,DATAin%nodays, DATAin%LAT &
                     ,fluxes2, pools2, diags2,  DATAin%nopars &
                     ,DATAin%nomet, DATAin%nopools, DATAin%nofluxes  &
                     ,DATAin%nodiags, mVs(i))
    call check(error, .not. ieee_is_nan(sum(fluxes1)) )
    call check(error, .not. ieee_is_nan(sum(pools1)) )
    call check(error, .not. ieee_is_nan(sum(diags1)) )
    call check(error, .not. ieee_is_nan(sum(fluxes2)) )
    call check(error, .not. ieee_is_nan(sum(pools2)) )
    call check(error, .not. ieee_is_nan(sum(diags2)) )
    flux_error = sum(abs(fluxes1-fluxes2))
    pool_error = sum(abs(pools1-pools2))
    diag_error = sum(abs(diags1-diags2))
    call check(error, flux_error < .0000000001)
    call check(error, pool_error < .0000000001)
    call check(error, diag_error < .0000000001)
    end do
    call destroy_carbon_model()
end subroutine test_model_group_repeat_evaluation

subroutine test_model_group_repeat_evaluation_3(error)
  !! Similar to original cardamom "Sanity check" function
  !! This variant with 4 x 3 evaluations instead of 4 x 2 passes
  !! when group_repeat_evaluation doesn't, mysterious
  use CARBON_MODEL_MOD, only: carbon_model
  use CARBON_MODEL_MEMORY, only: mvs
    use cardamom_structures, only: DATAin
    use samplers_shared, only : init_pars_random
  implicit none
  type(error_type), allocatable, intent(out):: error
  integer:: nchains
    double precision, dimension(PI%npars):: PARS
    double precision, dimension(:,:), allocatable:: pools1, pools2, pools3
    double precision, dimension(:,:), allocatable:: fluxes1, fluxes2, fluxes3
    double precision, dimension(:,:), allocatable:: diags1, diags2, diags3
    double precision:: pool_error, flux_error, diag_error
integer:: seed, i, nodays
    type(UNIF_VECTOR):: random_uniform
    nchains = 1
    call initialize(infile)
    nodays = DATAin%nodays
    allocate(pools1((nodays+1), DATAin%nopools), pools2((nodays+1), DATAin%nopools), pools3((nodays+1), DATAin%nopools))
    allocate(fluxes1(nodays, DATAin%nofluxes), fluxes2(nodays, DATAin%nofluxes), fluxes3(nodays, DATAin%nofluxes))
    allocate(diags1(nodays, DATAin%nodiags), diags2(nodays, DATAin%nodiags), diags3(nodays, DATAin%nodiags))
    call initialize_carbon_model(nchains)
    PARS = DATAin%parpriors(1:PI%npars)
    seed = random_int()  ! TODO record later
    call random_uniform%initialize_random(seed)
    call init_pars_random(PI, PARS, PI%fix_pars, random_uniform)
    do i = 1, nchains
    call carbon_model(1, DATAin%nodays, DATAin%MET, PARS, DATAin%deltat &
                     ,DATAin%nodays, DATAin%LAT &
                     ,fluxes1, pools1, diags1,  DATAin%nopars &
                     ,DATAin%nomet, DATAin%nopools, DATAin%nofluxes  &
                     ,DATAin%nodiags, mVs(i))
    call carbon_model(1, DATAin%nodays, DATAin%MET, PARS, DATAin%deltat &
                     ,DATAin%nodays, DATAin%LAT &
                     ,fluxes2, pools2, diags2,  DATAin%nopars &
                     ,DATAin%nomet, DATAin%nopools, DATAin%nofluxes  &
                     ,DATAin%nodiags, mVs(i))
    flux_error = sum(abs(fluxes1-fluxes2))
    pool_error = sum(abs(pools1-pools2))
    diag_error = sum(abs(diags1-diags2))
    call check(error, flux_error < .0000000001)
    call check(error, pool_error < .0000000001)
    call check(error, diag_error < .0000000001)
    call carbon_model(1, DATAin%nodays, DATAin%MET, PARS, DATAin%deltat &
                     ,DATAin%nodays, DATAin%LAT &
                     ,fluxes3, pools3, diags3,  DATAin%nopars &
                     ,DATAin%nomet, DATAin%nopools, DATAin%nofluxes  &
                     ,DATAin%nodiags, mVs(i))
    flux_error = sum(abs(fluxes3-fluxes2))
    pool_error = sum(abs(pools3-pools2))
    diag_error = sum(abs(diags1-diags2))
    call check(error, flux_error < .0000000001)
    call check(error, pool_error < .0000000001)
    end do
    call destroy_carbon_model()
end subroutine test_model_group_repeat_evaluation_3

subroutine test_model_cross_repeat_evaluation(error)
  !! Multiple threads that change threadid and evaluate the state again with mV arrays
  !! previously used by a different thread
  implicit none
  type(error_type), allocatable, intent(out):: error
end subroutine test_model_cross_repeat_evaluation

end module test_model
