! TODO move to tests
module MHMCMC_StressTests
   use samplers_shared, only: PARINFO
   use model_shared, only: PI

   ! Module contains a number of diagnostic tests used to ensure that the MCMC
   ! is able to retrieve a known distribution of parameters for simple models.

  !!!!!!!!!!!
   ! Authorship contributions
   !
   ! Created: 12/05/2021, T. L. Smallman (UoE, t.l.smallman@ed.ac.uk)
   ! Subsequent contributions by:
   ! T. L. Smallman (UoE)
   ! A. A. Bloom (JPL, USA)
   ! Version history:
   ! Version 1.0: A single parameter retrieval with known PDF.
   !              Estimating pi and radius from a circle with multiple 'observation'
   !              Estimating pi and radii from 9 circles, a single observtion per circle

   implicit none

   ! Assume all contents private unless explicitly states
   private

   ! Explicit statement of public variables or functions
   public:: prepare_for_stress_test, StressTest_likelihood_fct, StressTest_sublikelihood_fct

   ! Declare any module level variables

   ! Stress Test 1-estimate parameters for multiple circle
   ! Parameter 1 = pi, parameter 2:10 = radi
   ! DO NOT USE SOME REFINEMENT NEEDED
   double precision, parameter:: circle_par_1 = 3.141d0, &
                                 circle_par_2 = 1.2d0, &
                                 circle_par_3 = 3d0, &
                                 circle_par_4 = 8d0, &
                                 circle_par_5 = 10d0, &
                                 circle_par_6 = 15d0, &
                                 circle_par_7 = 200d0, &
                                 circle_par_8 = 193d0, &
                                 circle_par_9 = 88d0, &
                                circle_par_10 = 291d0, &
                               circle_obs_unc = 1d0
   double precision, dimension(9):: circle_obs

   ! Stress Test 2-estimate known PDF for single parameter
   double precision, parameter:: single_obs_mean = 0d0, &
                                  single_obs_unc = 1d0

   ! Stress Test 3-estimate parameters for a single circle
   ! Parameter 1 = pi, parameter 2 = radius
   double precision, parameter:: single_circle_par_1 = 3.141d0, &
                                 single_circle_par_2 = 88d0, &
                               single_circle_obs_unc = single_circle_par_2*0.2d0, &
                               single_circle_par_unc = single_circle_par_1*0.5d0
   double precision, dimension(29):: single_circle_obs

contains
   !
   !--------------------------------------------------------------------
   !
   subroutine circle(pars, nopars, output)

      implicit none

      ! Subroutine estimates the area of a circle using two parameters
      ! which define pi and the radius. A prior estimate of pi has also been applied.

      ! arguments
      integer, intent(in):: nopars
      double precision, dimension(nopars), intent(in):: pars
      double precision, intent(out):: output
      ! local variables
      integer:: i
      double precision, dimension(nopars - 1):: area

      ! Determine the area of the circle for the current parameters
      do i = 2, nopars
         area(i - 1) = pars(1)*pars(i)**2d0
      end do
      if (isnan(pars(1))) then
         stop 1
      end if
      ! Convert into log-likelihood
      output = sum(-0.5d0*(((area - circle_obs)/circle_obs_unc)**2))

   end subroutine circle
   !
   !--------------------------------------------------------------------
   !
   subroutine circle_parameter_prior_ranges
      use model_shared, only: PI

      ! define the parameter prior ranges

      implicit none

      ! Pi
      PI%parmin(1) = -10d0
      PI%parmax(1) = 10d0

      ! Radius-1
      PI%parmin(2) = 1d0
      PI%parmax(2) = 300.0d0

      ! Radius-2
      PI%parmin(3) = 1d0
      PI%parmax(3) = 300.0d0

      ! Radius-3
      PI%parmin(4) = 1d0
      PI%parmax(4) = 300.0d0

      ! Radius-4
      PI%parmin(5) = 1d0
      PI%parmax(5) = 300.0d0

      ! Radius-5
      PI%parmin(6) = 1d0
      PI%parmax(6) = 300.0d0

      ! Radius-6
      PI%parmin(7) = 1d0
      PI%parmax(7) = 300.0d0

      ! Radius-7
      PI%parmin(8) = 1d0
      PI%parmax(8) = 300.0d0

      ! Radius-8
      PI%parmin(9) = 1d0
      PI%parmax(9) = 300.0d0

      ! Radius-9
      PI%parmin(10) = 1d0
      PI%parmax(10) = 300.0d0

      ! Assign observations values
      circle_obs(1) = circle_par_1*circle_par_2**2d0
      circle_obs(2) = circle_par_1*circle_par_3**2d0
      circle_obs(3) = circle_par_1*circle_par_4**2d0
      circle_obs(4) = circle_par_1*circle_par_5**2d0
      circle_obs(5) = circle_par_1*circle_par_6**2d0
      circle_obs(6) = circle_par_1*circle_par_7**2d0
      circle_obs(7) = circle_par_1*circle_par_8**2d0
      circle_obs(8) = circle_par_1*circle_par_9**2d0
      circle_obs(9) = circle_par_1*circle_par_10**2d0

   end subroutine circle_parameter_prior_ranges
   !
   !--------------------------------------------------------------------
   !
   subroutine single_circle(pars, nopars, output)

      implicit none

      ! Subroutine estimates the area of a circle using two parameters
      ! which define pi and the radius. A prior estimate of pi has also been applied.

      ! arguments
      integer, intent(in):: nopars
      double precision, dimension(nopars), intent(in):: pars
      double precision, intent(out):: output

      ! Determine the area of the circle for the current parameters
      output = pars(1)*pars(2)**2d0

      ! Convert into log-likelihood
      output = sum(-0.5d0*(((output - single_circle_obs)/single_circle_obs_unc)**2))

   end subroutine single_circle
   !
   !--------------------------------------------------------------------
   !
   subroutine single_circle_parameter_prior_ranges

      ! define the parameter prior ranges

      implicit none

      ! Pi
      PI%parmin(1) = -10d0
      PI%parmax(1) = 10d0

      ! Radius-1
      PI%parmin(2) = 1d0
      PI%parmax(2) = 300.0d0

      ! Assign observations values, taken from a Gaussian distribution with mean of 88 and variance of 1
      single_circle_obs(1) = 24329.30d0
      single_circle_obs(2) = 24328.31d0
      single_circle_obs(3) = 24328.23d0
      single_circle_obs(4) = 24328.89d0
      single_circle_obs(5) = 24326.58d0
      single_circle_obs(6) = 24329.30d0
      single_circle_obs(7) = 24327.14d0
      single_circle_obs(8) = 24330.30d0
      single_circle_obs(9) = 24329.10d0
      single_circle_obs(10) = 24330.63d0
      single_circle_obs(11) = 24328.21d0
      single_circle_obs(12) = 24327.76d0
      single_circle_obs(13) = 24327.22d0
      single_circle_obs(14) = 24329.36d0
      single_circle_obs(15) = 24329.99d0
      single_circle_obs(16) = 24327.93d0
      single_circle_obs(17) = 24328.18d0
      single_circle_obs(18) = 24327.78d0
      single_circle_obs(19) = 24328.47d0
      single_circle_obs(20) = 24328.25d0
      single_circle_obs(21) = 24327.39d0
      single_circle_obs(22) = 24328.63d0
      single_circle_obs(23) = 24328.01d0
      single_circle_obs(24) = 24329.58d0
      single_circle_obs(25) = 24328.56d0
      single_circle_obs(26) = 24329.65d0
      single_circle_obs(27) = 24326.77d0
      single_circle_obs(28) = 24328.45d0
      single_circle_obs(29) = 24327.42d0

   end subroutine single_circle_parameter_prior_ranges
   !
   !--------------------------------------------------------------------
   !
   subroutine single_parameter_prior_ranges

      ! define the parameter prior ranges

      ! Estimate the PDF of a single value with known variance
      ! Test concept from A. A. Bloom (JPL, USA)
      ! Implemented by T. L. Smallman (UoE, t.l.smallman@ed.ac.uk)

      implicit none

      ! Very large uniform range
      PI%parmin(1) = -100000d0
      PI%parmax(1) = 100000d0

   end subroutine single_parameter_prior_ranges
   !
   !--------------------------------------------------------------------
   !
   subroutine prepare_for_stress_test(infile, outfile)
      use apmcmc, only: MCMC_OUTPUT, MCMC_OPTIONS
      use cardamom_structures, only: DATA_type, set_datain
      use cardamom_main_utils, only: initialize_stats

      ! Function by-passes the main CARDAMOM i/o code to allow
      ! for a non-standard operation of the model stress test

      implicit none

      ! Arguments
      character(350), intent(inout):: infile, outfile
      type(Data_type):: DATAin  ! local tmp copy

      ! local variables
      integer:: i
      type(MCMC_OUTPUT):: MCOUT  ! TODO no effect
      type(MCMC_OPTIONS):: MCO  ! TODO no effect

      ! Set internal parameters in the absence of an input file
      ! allocate the default run information

      if (outfile == "Circle") then
          ! ID = -1 StressTest-Circle
          DATAin%ID = -1
          DATAin%nodays = 1
          DATAin%nomet = 1
          DATAin%noobs = 9
          DATAin%nopools = 1
          DATAin%nopars = 10
          DATAin%nofluxes = 1
      else if (outfile == "Single") then
          ! ID = -2 StressTest-Single parameter
          DATAin%ID = -2
          DATAin%nodays = 1
          DATAin%nomet = 1
          DATAin%noobs = 1
          DATAin%nopools = 1
          DATAin%nopars = 1  ! 2
          DATAin%nofluxes = 1
      else if (outfile == "SingleCircle") then
          ! ID = -3 StressTest-Single Circle
          DATAin%ID = -3
          DATAin%nodays = 1
          DATAin%nomet = 1
          DATAin%noobs = 29
          DATAin%nopools = 1
          DATAin%nopars = 2
          DATAin%nofluxes = 1
      else
          print *, "Valid Stress Test has not been specified"
          stop 1
      end if

      ! Now we have used the infile to determine that this is going to be stress test,
      ! and the specific one has been determined from the outfile,
      ! we will now overwrite the outfile to give a default output location
      outfile = "stress_test_output_"

      ! need to allocate memory to the model output variables
      ! ->They are now local variables of model likelihood fcts
      if (.not. allocated(DATAin%M_FLUXES)) allocate (DATAin%M_FLUXES(DATAin%nodays, DATAin%nofluxes))
      if (.not. allocated(DATAin%M_POOLS)) allocate (DATAin%M_POOLS(DATAin%nodays + 1, DATAin%nopools))

      ! alert the user
      write (*, *) "Created fields for model output"

      ! Begin allocating parameter info
      PI%npars = DATAin%nopars
      if (.not. allocated(MCOUT%parvar)) call initialize_stats(MCOUT, PI%npars)

      if (.not. allocated(PI%parmin)) then
         allocate (PI%parmin(PI%npars), PI%parmax(PI%npars), &
                   PI%parfix(PI%npars), PI%paradj(PI%npars))
      end if
      if (.not. allocated(MCOUT%pars)) allocate (MCOUT%pars(PI%npars))

      ! force zero
      PI%parmin = 0d0; PI%parmax = 0d0; MCOUT%pars = 0d0
      PI%parfix = .false.

      ! load parameter max/min information
      if (DATAin%ID == -1) then
          call circle_parameter_prior_ranges
      else if (DATAin%ID == -2) then
          call single_parameter_prior_ranges
      else if (DATAin%ID == -3) then
          call single_circle_parameter_prior_ranges
      end if

      ! For log-normalisation procedure, no parameter can be <= 0.
      ! To facilitate easy of setting parameter ranges to real values
      ! we here instead calculate the adjustment need to ensure positive only values
      where (PI%parmin <= 0d0) PI%paradj = abs(PI%parmin) + 1d0

      call set_datain(DATAin)

   end subroutine prepare_for_stress_test
   !
   !------------------------------------------------------------------
   !
   subroutine stresstest_likelihood_fct(params, npars, loglikelihood, id) bind(c, name="C_stresstest_likelihood")
    !! wrapper to make stresstest_likelihood function compatible with cardamom-samplers lib, C, R
      use iso_c_binding
      implicit none
      integer(c_int), intent(in)  :: npars
      real(c_double), intent(inout), dimension(npars):: params
      real(c_double), intent(out):: loglikelihood
      integer(c_int), intent(in), optional:: id

      real(c_double):: ML_obs_out, ML_prior_out
      call stresstest_likelihood(params, ML_obs_out, ML_prior_out)

      loglikelihood = ML_obs_out + ML_prior_out
   end subroutine stresstest_likelihood_fct

   !
   !------------------------------------------------------------------
   !
   subroutine StressTest_likelihood(PARS, ML_obs_out, ML_prior_out, thread_id)
      use cardamom_structures, only: DATAin

      ! this subroutine is responsible, under normal circumstances for the running
      ! of the DALEC model, calculation of the log-likelihood for comparison
      ! assessment of parameter performance and use of the EDCs if they are
      ! present/selected

      implicit none

      ! declare inputs
      double precision, dimension(PI%npars), intent(inout):: PARS  ! current parameter vector
      ! output
      double precision, intent(inout):: ML_obs_out, &  ! observation+EDC log-likelihood
         ML_prior_out   ! prior log-likelihood
      integer, intent(in), optional:: thread_id

      ! local variables
      double precision:: output
      ! initial values
      ML_obs_out = 0d0; ML_prior_out = 0d0

      if (DATAin%ID == -1) then
         ! run the circle model
         call circle(PARS, DATAin%nopars, output)
         ! Estimate the likelihood score
         ML_obs_out = output + (-0.5d0*((((pars(1) - circle_par_1))/circle_obs_unc)**2))
      else if (DATAin%ID == -2) then
         ! Estimate likelihood for a single parameter retrieval
         ML_obs_out = -0.5d0*(((PARS(1) - single_obs_mean)/single_obs_unc)**2)
      else if (DATAin%ID == -3) then
         ! run the circle model
         call single_circle(PARS, DATAin%nopars, output)
         ! Estimate the likelihood score
         ML_obs_out = output + (-0.5d0*((((pars(1) - single_circle_par_1))/single_circle_par_unc)**2))
      end if

   end subroutine StressTest_likelihood
   !
   !------------------------------------------------------------------
   !
   subroutine stresstest_sublikelihood_fct(params, npars, loglikelihood, id)
    !! wrapper to make stresstest_likelihood function compatible with cardamom-samplers lib, C, R
      use iso_c_binding
      implicit none
      integer(c_int), intent(in)  :: npars
      real(c_double), intent(inout), dimension(npars):: params
      real(c_double), intent(out):: loglikelihood

      real(c_double):: ML_obs_out, ML_prior_out
      integer(c_int), intent(in), optional:: id

      call stresstest_sublikelihood(params, ML_obs_out, ML_prior_out)

      loglikelihood = ML_obs_out + ML_prior_out

   end subroutine stresstest_sublikelihood_fct
   !
   !------------------------------------------------------------------
   !
   subroutine StressTest_sublikelihood(PARS, ML_obs_out, ML_prior_out, thread_id)
      use cardamom_structures, only: DATAin
      ! this subroutine is responsible, under normal circumstances for the running
      ! of the DALEC model, calculation of the log-likelihood for comparison
      ! assessment of parameter performance and use of the EDCs if they are
      ! present/selected

      implicit none

      ! declare inputs
      double precision, dimension(PI%npars), intent(inout):: PARS  ! current parameter vector
      double precision, dimension(DATAin%nodays):: M_LAI, M_NEE, M_GPP
      double precision, dimension(DATAin%nodays, DATAin%nofluxes):: M_FLUXES
      double precision, dimension((DATAin%nodays + 1), DATAin%nopools):: M_POOLS
      ! output
      double precision, intent(inout):: ML_obs_out, &  ! observation+EDC log-likelihood
         ML_prior_out   ! prior log-likelihood
      integer, intent(in), optional:: thread_id

      ! local variables
      double precision:: output
      ! initial values
      ML_obs_out = 0d0; ML_prior_out = 0d0

      if (DATAin%ID == -1) then
         ! run the circle model
         call circle(PARS, DATAin%nopars, output)
         ! Estimate the likelihood score, normalise by sample size in the sub-case
         ML_obs_out = (output/DATAin%noobs) + (-0.5d0*((((pars(1) - circle_par_1))/circle_obs_unc)**2))
      else if (DATAin%ID == -2) then
         ! Estimate likelihood for a single parameter retrieval
         ML_obs_out = -0.5d0*(((PARS(1) - single_obs_mean)/single_obs_unc)**2)
      else if (DATAin%ID == -3) then
         ! run the circle model
         call single_circle(PARS, DATAin%nopars, output)
         ! Estimate the likelihood score
         ML_obs_out = (output/DATAin%noobs) + (-0.5d0*((((pars(1) - single_circle_par_1))/single_circle_par_unc)**2))
      end if

   end subroutine StressTest_sublikelihood
   !
   !--------------------------------------------------------------------
   !
end module MHMCMC_StressTests  ! MHMCMC_StressTests
