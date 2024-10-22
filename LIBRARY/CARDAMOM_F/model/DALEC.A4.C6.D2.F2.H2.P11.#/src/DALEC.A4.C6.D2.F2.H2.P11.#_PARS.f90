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

    ! NOTE: that these parameter ranges have been matched with Bloom's C code
    ! 22/11/2019 - try not to lose this information as it is needed for comparability

    !
    ! declare parameters
    !

    ! Decomposition litter -> som (day-1)
    PI%parmin(1) = 0.00001d0
    PI%parmax(1) = 0.01d0

    ! Fraction of GPP respired as autotrophic
    PI%parmin(2) = 0.2d0
    PI%parmax(2) = 0.8d0

    ! Potential rate of direct, i.e. without CDEA control labile to foliage (gC/m2/day)
    PI%parmin(3) = 0.001d0
    PI%parmax(3) = 8d0

    ! Potential rate of direct labile to fine root (gC/m2/day)
    PI%parmin(4) = 0.1d0
    PI%parmax(4) = 8d0

    ! Leaf Lifespan (yr)
    ! Wright et al. 2004
    PI%parmin(5) = 1.001d0
    PI%parmax(5) = 6d0 !8d0

    ! TOR wood* - 1% loss per year value
    PI%parmin(6) = 0.000009d0 ! 304  years
    PI%parmax(6) = 0.001d0    ! 2.74 years

    ! TOR roots
    PI%parmin(7) = 0.001368925d0 ! 2    years !0.0006844627d0 ! 4 years
    PI%parmax(7) = 0.02d0        ! 0.13 years

    ! Turnover of litter (fraction; temperature adjusted)
    PI%parmin(8) = 0.0001141d0 ! 24   years at 0oC
    PI%parmax(8) = 0.02d0      ! 0.13 years at 0oC

    ! Turnover of som to Rhet (fraction; temperature adjusted)
    PI%parmin(9) = 1.368925d-06   ! 2000 years at 0oC
    PI%parmax(9) = 9.126169d-05   !   30 years at 0oC !0.0001368926d0 !   20 years at 0oC
!    PI%parmin(9) = 0.0000001d0 ! 27378.0 years at 0oC
!    PI%parmax(9) = 0.001d0     !     2.7 years at 0oC

    ! Temp factor* = Q10 = 1.2-1.6
    PI%parmin(10) = 0.019d0
    PI%parmax(10) = 0.08d0

    ! Vcmax, the maximum rate of carboxylation at the canopy top
    ! umolC/m2/s
    PI%parmin(11) = 10d0
    PI%parmax(11) = 100d0

    ! max bud burst day
    PI%parmin(12) = 365.25d0
    PI%parmax(12) = 365.25d0*4d0

    ! Potential rate of seasonal labile to foliage (gC/m2/day)
    PI%parmin(13) = 0.001d0
    PI%parmax(13) = 10d0

    ! Clab Release period
    PI%parmin(14) = 10d0
    PI%parmax(14) = 100d0

    ! max leaf fall day
    PI%parmin(15) = 365.25d0
    PI%parmax(15) = 365.25d0*4d0

    ! Leaf fall period
    PI%parmin(16) = 20d0
    PI%parmax(16) = 150d0

    ! LMA (gC.m-2)
    ! Kattge et al. 2011
    PI%parmin(17) = 20d0
    PI%parmax(17) = 180d0

    ! fraction of Cwood which is coarse root
    PI%parmin(25) = 0.15d0
    PI%parmax(25) = 0.50d0

    ! BUCKET - coarse root biomass (i.e. gbio/m2 not gC/m2) needed to reach 50 %
    ! of max depth
    PI%parmin(26) = 100d0
    PI%parmax(26) = 2500d0 !500d0

    ! BUCKET - maximum rooting depth
    PI%parmin(27) = 0.35d0
    PI%parmax(27) = 20d0

    ! Resilience factor for burned but not combusted C stocks
    PI%parmin(28) = 0.01d0
    PI%parmax(28) = 0.99d0
    ! Combustion completeness factor for foliage
    PI%parmin(29) = 0.01d0
    PI%parmax(29) = 0.99d0
    ! Combustion completeness factor for fine root and wood
    PI%parmin(30) = 0.01d0
    PI%parmax(30) = 0.99d0
    ! Combustion completeness factor for soil
    PI%parmin(31) = 0.01d0
    PI%parmax(31) = 0.1d0
    ! Combustion completeness factor for foliage + fine root litter
    PI%parmin(32) = 0.01d0
    PI%parmax(32) = 0.99d0

    ! labile:biomass at which growth is limited by 50 %
    PI%parmin(33) = 0.005d0 ! 0.5 %
    PI%parmax(33) = 0.1d0   ! 10 %

    ! Temperature (oC) above p36 at which fine root growth is limited by 50 %
    PI%parmin(34) = 1d0
    PI%parmax(34) = 10d0
    ! Temperature (oC) above p37 at which wood growth is limited by 50 %
    PI%parmin(35) = 1d0
    PI%parmax(35) = 10d0
    ! Temperature (oC) at which fine root growth is prevented
    PI%parmin(36) =  0.01d0 !-8d0
    PI%parmax(36) =  8d0
    ! Temperature (oC) at which wood growth is prevented
    PI%parmin(37) = 1d0
    PI%parmax(37) = 8d0

    ! Potential growth rate of wood (gC/m2/day)
    PI%parmin(38) = 0.05d0
    PI%parmax(38) = 8d0

    ! Potential supply of water from roots (mmolH2O/m2/s) at which wood growth is limited by 50 %
    PI%parmin(39) = 0.001d0
    PI%parmax(39) = 5d0
    ! Potential supply of water from roots (mmolH2O/m2/s) at which wood growth is prevented
    PI%parmin(40) = 0.001d0
    PI%parmax(40) = 5d0

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
    PI%parmin(20) = 1.0d0
    PI%parmax(20) = 2000d0

    ! C_wood
    PI%parmin(21) = 1d0
    PI%parmax(21) = 30000d0

    ! C litter
    PI%parmin(22) = 1d0
    PI%parmax(22) = 2000d0

    ! C_som
    PI%parmin(23) = 200d0
    PI%parmax(23) = 250000d0 !90000d0

    ! Initial soil water
    ! a fraction of field capacity
    PI%parmin(24) = 0.50d0
    PI%parmax(24) = 1.00d0

  end subroutine pars_info

  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
