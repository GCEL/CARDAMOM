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

    ! contains 6 fields with min max log for par and par

    !
    ! declare parameters
    !

    ! Fraction of GPP respired as autotrophic
    PI%parmin(1) = 0.20d0
    PI%parmax(1) = 0.80d0

    ! Fraction of (1-fgpp) to foliage
    PI%parmin(2) = 0.01d0
    PI%parmax(2) = 0.5d0

    ! Leaf Lifespan (yr)
    ! Wright et al. (2004) 55 - 2922 days
    PI%parmin(3) = 0.15d0
    PI%parmax(3) = 8d0

    ! TOR wood + root
    PI%parmin(4) = 0.000009d0 ! 304  years
    PI%parmax(4) = 0.01d0     ! 0.27 years

    ! TOR litter + som
    PI%parmin(5) = 1.368925d-06   ! 2000 years at 0oC
    PI%parmax(5) = 9.126169d-05   !   30 years at 0oC !0.0001368926d0 !   20 years at 0oC

    ! Temp factor* = Q10 = 1.2-1.6
    PI%parmin(6) = 0.019d0
    PI%parmax(6) = 0.08d0

    ! Canopy Efficiency (gC/m2leaf/day)
    PI%parmin(7) = 5d0
    PI%parmax(7) = 50d0

    ! LMA (gC.m-2)
    ! Kattge et al. 2011
    PI%parmin(8) = 5d0
    PI%parmax(8) = 200d0

    !
    ! INITIAL VALUES DECLARED HERE
    !

    ! C foliar
    PI%parmin(9) = 1d0
    PI%parmax(9) = 2000d0

    ! C root + wood
    PI%parmin(10) = 1d0
    PI%parmax(10) = 100000d0

    ! C_som
    PI%parmin(11) = 1d0
    PI%parmax(11) = 200000d0

  end subroutine pars_info
  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
