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
    use cardamom_structures, only: DATAin

    ! Subroutine contains a list of parameter ranges for the model.
    ! These could or
    ! possibly should go into an alternate file which can be read in.
    ! This may
    ! improve the usability when it comes to reading these information
    ! in for
    ! different PFTs

    implicit none

!    PI%npars=26 ! dont forget to change in cardamom_io.f90

    ! contains 6 fields with min max log for par and par

    !
    ! declare parameters
    !

    ! Decomposition rate
    PI%parmin(1)=0.00001
    PI%parmax(1)=0.01

    ! Fraction of GPP respired
    PI%parmin(2)=0.3
    PI%parmax(2)=0.7

    ! Fraction of (1-fgpp) to foliage
    PI%parmin(3)=0.01
    PI%parmax(3)=0.5

    ! Fraction of (1-fgpp) to roots*/
    PI%parmin(4)=0.01
    PI%parmax(4)=1.

    ! Leaf Lifespan
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
    PI%parmin(10)=0.018
    PI%parmax(10)=0.08

    ! Canopy Efficiency
    ! set to parmin=1 for FLUXCOM only
    ! e.g. for wetlands etc.
    PI%parmin(11)=10.
    PI%parmax(11)=100.

    ! Bday
    PI%parmin(12)=365.25
    PI%parmax(12)=365.25*4

    ! Fraction to Clab*/
    PI%parmin(13)=0.01
    PI%parmax(13)=0.5

    ! Clab Release period
    PI%parmin(14)=10.
    PI%parmax(14)=100.

    ! Fday
    PI%parmin(15)=365.25
    PI%parmax(15)=365.25*4

    ! Leaf fall period
    PI%parmin(16)=20.
    PI%parmax(16)=150.

    ! LMA
    ! Kattge et al. 2011,
    ! Restricted based on Duke forest site information
    PI%parmin(17)=10.
    PI%parmax(17)=400.

    !
    ! INITIAL VALUES DECLARED HERE
    !

    ! initial conditions updated based on
    ! Adegbidi et al (2005)

    ! C labile
    PI%parmin(18)=1.0
    PI%parmax(18)=3000.0 ! sum of all other living max

    ! C foliar
    PI%parmin(19)=1.0
    PI%parmax(19)=1000.0

    ! C roots
    PI%parmin(20)=1.0
    PI%parmax(20)=1000.0

    ! C_wood
    PI%parmin(21)=1.0
    PI%parmax(21)=1000.0

    ! C litter
    PI%parmin(22)=1.0
    PI%parmax(22)=1000.0

    ! C_som
    PI%parmin(23)=100.0
    PI%parmax(23)=200000.0

  end subroutine pars_info

  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
