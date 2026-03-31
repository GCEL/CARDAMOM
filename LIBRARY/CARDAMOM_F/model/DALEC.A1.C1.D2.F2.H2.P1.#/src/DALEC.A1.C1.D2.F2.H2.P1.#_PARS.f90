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

    ! Fraction of (1-fgpp) to foliage
    PI%parmin(3) = 0.1d0
    PI%parmax(3) = 0.5d0

    ! Fraction of (1-fgpp) to roots*/
    PI%parmin(4) = 0.1d0
    PI%parmax(4) = 0.80d0

    ! Leaf Lifespan (yr)
    ! Wright et al. 2004
    PI%parmin(5) = 1.001d0
    PI%parmax(5) = 13d0 !8d0

    ! TOR wood* - 1% loss per year value
    PI%parmin(6) = 0.000009d0 ! 304  years
    PI%parmax(6) = 0.001d0    ! 2.74 years

    ! TOR roots
    PI%parmin(7) = 0.00089181d0             !0.002191781   =1.25 years   |   0.001368925d0  = 2    years (default)  | 0.0006844627d0 = 4 years |  0.00130374  =2.1 years #Solly and Yuan papers | 0.00089181 = 3.07 years # Hansson 2013
    PI%parmax(7) = 0.02d0             !0.003652968   =0.75 years   |   0.02d0        = 0.13 years (default)    |                         |  0.00391122  =0.7 years # Solly, 2018 & Yuan, 2010

    ! Turnover of litter (fraction; temperature adjusted)
    PI%parmin(8) = 0.0001141d0 ! 24   years at 0oC
    PI%parmax(8) = 0.02d0     ! 0.033 = 30 day turnover             ! default = 0.02d0      ! 0.13 years at 0oC

    ! Turnover of som to Rhet (fraction; temperature adjusted)
    PI%parmin(9) = 1.368925d-06   ! 2000 years at 0oC
    PI%parmax(9) = 9.126169d-05!9.126169d-05 -> DEFAULT  !   30 years at 0oC !0.0001368926d0 !   20 years at 0oC
!    PI%parmin(9) = 0.0000001d0 ! 27378.0 years at 0oC
!    PI%parmax(9) = 0.001d0     !     2.7 years at 0oC

    ! Temp factor* = Q10 = 1.2-1.6
    PI%parmin(10) =  0.019d0 ! 0.19868d0 =  Q10 = 2.7 (Karhu, 2010) | default =  0.019d0
    PI%parmax(10) = 0.08d0 ! 0.23265d0 = Q10 = 3.2 (Karhu, 2010) | default = 0.08d0

    ! Canopy Efficiency
    ! NUE and avN combination give a Vcmax equivalent, the canopy efficiency.
    ! Kattge et al (2011) offers a potential prior range of 3.4 - 30.7 gC/m2leaf/day.
    ! Here, to be cautious we will expand accepted range
    ! Thus CUE = NUE * avN -> 1.64 / 42.0
    ! TLS: 27/10/2021 restricted again based now on 95 %CI (12.61 / 29.68) from TRY
    PI%parmin(11) = 1.64d0 !5d0 ! current default (27/1/2024) is 10 TG, has been lowered to 1.64 for N site runs....
    PI%parmax(11) = 100d0 !42d0 !50d0

    ! max bud burst day
    PI%parmin(12) = 60.0d0 !60.0d0!365.25d0 !90.0d0 TG - no longer necessary with new max LAI summer EDC
    PI%parmax(12) = 183.0d0!365.25d0*4d0 ! 180.0d0 TG

    ! Fraction to Clab*/
    PI%parmin(13) = 0.01d0
    PI%parmax(13) = 0.5d0

    ! Clab Release period
    PI%parmin(14) = 10d0
    PI%parmax(14) = 100d0

    ! max leaf fall day
    PI%parmin(15) = 183d0 !365.25d0  !183.0d0 TG - no longer necessary with new max LAI summer EDC
    PI%parmax(15) = 365.25d0   !365.25d0 TG

    ! Leaf fall period
    PI%parmin(16) = 20d0 !20d0 DEF
    PI%parmax(16) = 150d0 !150d0 DEF

    ! LMA (gC.m-2)
    ! Kattge et al. 2011
    PI%parmin(17) = 20d0
    PI%parmax(17) = 180d0

    ! fraction of Cwood which is coarse root
    PI%parmin(25) = 0.15d0
    PI%parmax(25) = 0.50d0 !default = 0.5

    ! BUCKET - coarse root biomass (i.e. gbio/m2 not gC/m2) needed to reach 50 %
    ! of max depth
    PI%parmin(26) = 100d0
    PI%parmax(26) = 2500d0 ! default 2500d0

    ! BUCKET - maximum rooting depth
    PI%parmin(27) = 0.35d0
    PI%parmax(27) =  4d0 !default = 20d0

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
