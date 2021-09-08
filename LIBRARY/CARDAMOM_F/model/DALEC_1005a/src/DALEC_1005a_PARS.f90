module MODEL_PARAMETERS

  implicit none

  !!!!!!!!!!!
  ! Authorship contributions
  !
  ! This code is based on the original C verion of the University of Edinburgh
  ! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
  ! All code translation into Fortran, integration into the University of
  ! Edinburgh CARDAMOM code and subsequent modifications by:
  ! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
  ! See function / subroutine specific comments for exceptions and contributors
  !!!!!!!!!!!

  ! make all private
  private

  ! specify explicitly the public
  public :: pars_info

  contains

  !
  !------------------------------------------------------------------
  !
  subroutine pars_info
    use MCMCOPT, only: PI

    ! Subroutine contains a list of parameter ranges for the model.
    ! These could or possibly should go into an alternate file which can be read in.
    ! This may improve the usability when it comes to reading these information
    ! in for different PFTs

    implicit none

    !
    ! declare parameters
    !

    ! Decomposition litter -> som (day-1 at mean temperature)
    PI%parmin(1) = 0.00001d0
    PI%parmax(1) = 0.01d0

    ! Fraction of GPP respired as autotrophic
    PI%parmin(2) = 0.2d0
    PI%parmax(2) = 0.8d0

    ! Fraction of (1-fgpp) to foliage
    PI%parmin(3) = 0.01d0
    PI%parmax(3) = 0.5d0

    ! Fraction of (1-fgpp) to roots*/
    PI%parmin(4) = 0.01d0
    PI%parmax(4) = 1.0d0

    ! Leaf Lifespan (yr)
    ! Wright et al. 2004
    PI%parmin(5) = 1.001d0
    PI%parmax(5) = 8d0

    ! TOR wood* - 1% loss per year value
    !PI%parmin(6) = 0.000009d0 ! 304  years
    PI%parmin(6) = 0.000025d0 ! 109  years
    PI%parmax(6) = 0.001d0    ! 2.74 years

    ! TOR roots
    !PI%parmin(7) = 0.001368925d0 ! 2    years
    !PI%parmax(7) = 0.02d0        ! 0.13 years
    PI%parmin(7) = 0.0001d0 ! 27    years
    PI%parmax(7) = 0.01d0   ! 0.27 years

    ! Turnover of litter (fraction; temperature adjusted)
    !PI%parmin(8) = 0.0001141d0 ! 24   years at 0oC
    !PI%parmax(8) = 0.02d0      ! 0.13 years at 0oC
    PI%parmin(8) = 0.0001d0 ! 27   years at 0oC
    PI%parmax(8) = 0.01d0   ! 0.27 years at 0oC

    ! Turnover of som to Rhet (fraction; temperature adjusted)
!    PI%parmin(9) = 1.368925d-06   ! 2000 years at 0oC
!    PI%parmax(9) = 9.126169d-05   !   30 years at 0oC
    PI%parmin(9) = 1.0d-07  ! 2737851   years at 0oC
    PI%parmax(9) = 1.0d-3   !       2.7 years at 0oC

    ! Temp factor* = Q10 = 1.2-1.6
    PI%parmin(10) = 0.018d0
    PI%parmax(10) = 0.08d0

    ! Canopy Efficiency
    ! NUE and avN combination give a Vcmax equivalent, the canopy efficiency.
    ! Kattge et al (2011) offers a prior of 3.4 - 30.7 gC/m2leaf/day.
    ! Here, to be cautious we will expand accepted range
    ! Thus CUE = NUE * avN -> 1.64 / 42.0
    PI%parmin(11) = 1.64d0 !5d0
    PI%parmax(11) = 42d0 !50d0

    ! max bud burst day
    PI%parmin(12) = 365.25d0
    PI%parmax(12) = 365.25d0*4d0

    ! Fraction to Clab*/
    PI%parmin(13) = 0.01d0
    PI%parmax(13) = 0.5d0

    ! Clab Release period
    !PI%parmin(14) = 10d0
    PI%parmin(14) = 30d0
    PI%parmax(14) = 100d0

    ! max leaf fall day
    PI%parmin(15) = 365.25d0
    PI%parmax(15) = 365.25d0*4d0

    ! Leaf fall period
    !PI%parmin(16) = 20d0
    PI%parmin(16) = 30d0
    PI%parmax(16) = 150d0

    ! LMA (gC.m-2)
    ! Kattge et al. 2011
    PI%parmin(17) = 20d0
    PI%parmax(17) = 180d0

    ! uWUE: GPP*sqrt(VPD)/ ET
    ! gC/kgH2O per hPa
    ! Boese et al., 2017
    PI%parmin(24) = 0.5d0
    PI%parmax(24) = 30d0

    ! Inverse of second order runoff constant for plant available water pool
    ! See Bloom et al., (2020) https://doi.org/10.5194/bg-17-6393-2020
    ! (mm day)
    PI%parmin(25) = 1d0
    PI%parmax(25) = 100000d0

    ! Plant wilting point (mm)
    ! Here this is considered as the point at which water limitation begins
    PI%parmin(26) = 1d0
    PI%parmax(26) = 10000d0

    ! Combustion completeness factor for foliage
    PI%parmin(28) = 0.01d0
    PI%parmax(28) = 1d0
    ! Combustion completeness factor for fine root and wood
    PI%parmin(29) = 0.01d0
    PI%parmax(29) = 1d0
    ! Combustion completeness factor for soil
    PI%parmin(30) = 0.01d0
    PI%parmax(30) = 1d0
    ! Resilience factor for burned but not combusted C stocks
    PI%parmin(31) = 0.01d0
    PI%parmax(31) = 1d0

    ! Labile pool lifespan (years)
    PI%parmin(32) = 1.001d0
    PI%parmax(32) = 8d0

    ! Moisture response coefficient for heterotrophic activity
    PI%parmin(33) = 0.01d0
    PI%parmax(33) = 1d0

    ! Fraction of draingage from plant available water that drains to
    ! plant unavailable
    PI%parmin(34) = 0.01d0
    PI%parmax(34) = 1d0

    ! Inverse of second order runoff constant for plant unavailable water pool
    ! See Bloom et al., (2020) https://doi.org/10.5194/bg-17-6393-2020
    ! (mm day)
    PI%parmin(35) = 1d0
    PI%parmax(35) = 100000d0

    ! Radiation coefficient modifying the uWUE calculation
    ! Boese et al., (2017)
    PI%parmin(37) = 0.01d0
    PI%parmax(37) = 0.3d0

    !
    ! INITIAL VALUES DECLARED HERE
    !

    ! C labile
    PI%parmin(18) = 1d0
    PI%parmax(18) = 2000d0

    ! C foliar
    PI%parmin(19) = 1d0
    PI%parmax(19) = 2000d0

    ! C roots
    PI%parmin(20) =  1.0d0
    PI%parmax(20) = 2000d0

    ! C_wood
    PI%parmin(21) = 1d0
    !PI%parmax(21) = 30000d0
    PI%parmax(21) = 100000d0

    ! C litter
    PI%parmin(22) = 1d0
    PI%parmax(22) = 2000d0

    ! C_som
    PI%parmin(23) = 200d0
    PI%parmax(23) = 250000d0

    ! Plant available water (mm)
    PI%parmin(27) = 1d0
    PI%parmax(27) = 10000d0

    ! Plant available water (mm)
    PI%parmin(36) = 1d0
    PI%parmax(36) = 10000d0

  end subroutine pars_info
  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
