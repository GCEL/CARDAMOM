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

!    ! Decomposition litter -> som (day-1)
!    PI%parmin(1) = 0.00001d0
!    PI%parmax(1) = 0.01d0
    ! Decomposition efficiency of litter/CWD to som (fraction)
    PI%parmin(1) = 0.25d0
    PI%parmax(1) = 0.75d0

    ! Fraction of GPP respired as Rm(fol,root,wood)
    PI%parmin(2) = 0.1d0
    PI%parmax(2) = 0.8d0

    ! Fraction of (1-fgpp) to foliage
    PI%parmin(3) = 0.01d0
    PI%parmax(3) = 0.5d0

    ! Fraction of (1-fgpp) to roots*/
    PI%parmin(4) = 0.1d0
    PI%parmax(4) = 0.80d0

    ! Leaf Lifespan (yr)
    ! Wright et al. 2004
    PI%parmin(5) = 1.001d0
    PI%parmax(5) = 8d0

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

    ! Canopy Efficiency
    ! NUE and avN combination give a Vcmax equivalent, the canopy efficiency.
    ! Kattge et al (2011) offers a prior of 3.4 - 30.7 gC/m2leaf/day.
    ! Here, to be cautious we will expand accepted range
    ! Thus CUE = NUE * avN -> 1.64 / 42.0
    PI%parmin(11) = 20!1.64d0 !5d0
    PI%parmax(11) = 42d0 !50d0
    ! log10 avg foliar N (gN.m-2)
    ! Kattge et al., (2011) (Quantiles 25% / 75%)
    ! and Thomas et al., (2019) (Aconite canopy paper)
!    PI%parmin(11) = 0.07918125d0!0d0 !-0.2218487d0 !TLS: restricted to 1.2 gN/m2leaf
!    PI%parmax(11) = 0.4771213d0 ! 0.5563025d0 ! TLS: restricted to 3 gC/m2leaf

    ! max bud burst day
    PI%parmin(12) = 365.25d0
    PI%parmax(12) = 365.25d0*4d0

    ! Fraction to Clab*/
    PI%parmin(13) = 0.01d0
    PI%parmax(13) = 0.5d0

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

!    ! Optimum Nitrogen use efficiency (gC/gN/m2/day)
!    PI%parmin(28) =  1.6d0
!    PI%parmax(28) = 40.0d0

    ! Turnover rate for CWD
    PI%parmin(29) = 1.368925d-05 ! 200.00 years at 0oC
    PI%parmax(29) = 0.001d0      !   2.74 years at 0oC

    ! Half saturation coefficient for self-thinning supression on wood turnover
    PI%parmin(30) = 1000d0
    PI%parmax(30) = 20000d0

    ! Resilience factor for burned but not combusted C stocks
    PI%parmin(31) = 0.1d0
    PI%parmax(31) = 1d0
    ! Combustion completeness factor for foliage
    PI%parmin(32) = 0.01d0
    PI%parmax(32) = 1d0
    ! Combustion completeness factor for fine root and wood
    PI%parmin(33) = 0.01d0
    PI%parmax(33) = 1d0
    ! Combustion completeness factor for soil
    PI%parmin(34) = 0.01d0
    PI%parmax(34) = 1d0

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

    ! C CWD
    PI%parmin(28) = 1d0
    PI%parmax(28) = 10000d0

  end subroutine pars_info

  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
