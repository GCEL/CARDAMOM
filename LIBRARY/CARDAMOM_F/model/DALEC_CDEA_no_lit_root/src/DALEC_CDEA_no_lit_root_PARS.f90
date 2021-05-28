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

!    PI%npars=23;

    ! NOTE: that these parameter ranges have been matched with Bloom's C code
    ! 22/11/2019 - try not to lose this information as it is needed for comparability

    !
    ! declare parameters
    !

    ! Fraction of GPP respired as autotrophic
    PI%parmin(1) = 0.2d0
    PI%parmax(1) = 0.8d0

    ! Fraction of (1-fgpp) to foliage
    PI%parmin(2) = 0.01d0
    PI%parmax(2) = 0.5d0

    ! Leaf Lifespan (yr)
    ! Wright et al. 2004
    PI%parmin(3) = 1.001d0
    PI%parmax(3) = 8d0

    ! TOR wood* - 1% loss per year value
    PI%parmin(4) = 0.000025d0 ! 109  years
    PI%parmax(4) = 0.01d0     ! 0.27 years

    ! Turnover of som to Rhet (fraction; temperature adjusted)
    PI%parmin(5) = 1.368925d-06   ! 2000 years at 0oC
    PI%parmax(5) = 0.01d0      !     0.27 years at 0oC

    ! Temp factor* = Q10 = 1.2-1.6
    PI%parmin(6) = 0.019d0
    PI%parmax(6) = 0.08d0

    ! Canopy Efficiency
    ! NUE and avN combination give a Vcmax equivalent, the canopy efficiency.
    ! Kattge et al (2011) offers a prior of 3.4 - 30.7 gC/m2leaf/day.
    ! Here, to be cautious we will expand accepted range
    ! Thus CUE = NUE * avN -> 1.64 / 42.0
    PI%parmin(7) = 1.64d0 !5d0
    PI%parmax(7) = 42d0 !50d0

    ! max bud burst day
    PI%parmin(8) = 365.25d0
    PI%parmax(8) = 365.25d0*4d0

    ! Fraction to Clab*/
    PI%parmin(9) = 0.01d0
    PI%parmax(9) = 0.5d0

    ! Clab Release period
    PI%parmin(10) = 10d0
    PI%parmax(10) = 100d0

    ! max leaf fall day
    PI%parmin(11) = 365.25d0
    PI%parmax(11) = 365.25d0*4d0

    ! Leaf fall period
    PI%parmin(12) = 20d0
    PI%parmax(12) = 150d0

    ! LMA (gC.m-2)
    ! Kattge et al. 2011
    PI%parmin(13) = 5d0
    PI%parmax(13) = 200d0

    !
    ! INITIAL VALUES DECLARED HERE
    !

    ! C labile
    PI%parmin(14) = 1d0
    PI%parmax(14) = 2000d0

    ! C foliar
    PI%parmin(15) = 1d0
    PI%parmax(15) = 2000d0

    ! C_wood + roots
    PI%parmin(16) = 1d0
    PI%parmax(16) = 100000d0

    ! C_som + litter
    PI%parmin(17) = 1d0
    PI%parmax(17) = 200000d0

  end subroutine pars_info

  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
