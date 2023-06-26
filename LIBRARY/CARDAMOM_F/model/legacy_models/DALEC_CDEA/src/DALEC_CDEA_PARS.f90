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
    ! These could or
    ! possibly should go into an alternate file which can be read in.
    ! This may
    ! improve the usability when it comes to reading these information
    ! in for
    ! different PFTs

    implicit none

    PI%npars=23;

    ! contains 6 fields with min max log for par and par

    !
    ! declare parameters
    !

    ! Decomposition litter -> som (day-1)
    PI%parmin(1)=0.00001
    PI%parmax(1)=0.01

    ! Fraction of GPP respired as autotrophic
    PI%parmin(2)=0.3
    PI%parmax(2)=0.7

    ! Fraction of (1-fgpp) to foliage
    PI%parmin(3)=0.01
    PI%parmax(3)=0.5

    ! Fraction of (1-fgpp) to roots*/
    PI%parmin(4)=0.01
    PI%parmax(4)=1.

    ! Leaf Lifespan (yr)
    ! Wright et al. 2004
    PI%parmin(5)=1.001
    PI%parmax(5)=8.

    ! TOR wood* - 1% loss per year value
    PI%parmin(6)=0.000025
    PI%parmax(6)=0.001

    ! TOR roots
    PI%parmin(7)=0.0001
    PI%parmax(7)=0.01

    ! TOR litter
    PI%parmin(8)=0.0001
    PI%parmax(8)=0.01

    ! TOR SOM
    PI%parmin(9)=0.0000001
    PI%parmax(9)=0.001

    ! Temp factor* = Q10 = 1.2-1.6
    PI%parmin(10)=0.019
    PI%parmax(10)=0.08

    ! Canopy Efficiency
    ! set to parmin=1 for FLUXCOM only
    ! e.g. for wetlands etc.
    PI%parmin(11)=10.
    PI%parmax(11)=100.

    ! max bud burst day
    PI%parmin(12)=365.25
    PI%parmax(12)=365.25*4

    ! Fraction to Clab*/
    PI%parmin(13)=0.01
    PI%parmax(13)=0.5

    ! Clab Release period
    PI%parmin(14)=10.
    PI%parmax(14)=100.

    ! max leaf fall day
    PI%parmin(15)=365.25
    PI%parmax(15)=365.25*4

    ! Leaf fall period
    PI%parmin(16)=20.
    PI%parmax(16)=150.

    ! LMA (gC.m-2)
    ! Kattge et al. 2011
    PI%parmin(17)=10d0
    PI%parmax(17)=200d0

    !
    ! INITIAL VALUES DECLARED HERE
    !

    ! C labile
    PI%parmin(18)=20.0
    PI%parmax(18)=2000.0

    ! C foliar
    PI%parmin(19)=20.0
    PI%parmax(19)=2000.0

    ! C roots
    PI%parmin(20)=20.0
    PI%parmax(20)=2000.0

    ! C_wood
    PI%parmin(21)=100.0
    PI%parmax(21)=100000.0

    ! C litter
    PI%parmin(22)=20.0
    PI%parmax(22)=2000.0

    ! C_som
    PI%parmin(23)=100.0
    PI%parmax(23)=200000.0

  end subroutine pars_info

  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
