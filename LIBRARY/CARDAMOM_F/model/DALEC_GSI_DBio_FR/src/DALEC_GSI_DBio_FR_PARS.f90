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
  ! G. Xenakis (University of Edinburgh, now at Forest Research, UK)
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

!    PI%npars=43 ! dont forget to change in cardamom_io.f90

    ! contains 6 fields with min max log for par and par

    !
    ! declare parameters
    !

    ! Efficiency of decomposition of the first stage decomposition
    PI%parmin(1)=0.45
    PI%parmax(1)=0.55

    ! Fraction of GPP respired
    PI%parmin(2)=0.2
    PI%parmax(2)=0.7

    ! GSI sensitivity for leaf growth
    PI%parmin(3)=1.00
    PI%parmax(3)=1.03  !1.05

    ! Fraction of (1-fgpp) to roots*/
    PI%parmin(4)=0.01
    PI%parmax(4)=1.0

    ! GSI max leaf turnover
    PI%parmin(5)=0.0001
    PI%parmax(5)=0.2

    ! TOR wood* - 1% loss per year value
    PI%parmin(6)=0.00001
    PI%parmax(6)=0.001

    ! TOR roots
    PI%parmin(7)=0.0001
    PI%parmax(7)=0.01

    ! TOR foliar litter
    PI%parmin(8)=0.0001
    PI%parmax(8)=0.01

    ! Microbial decomposition efficiency
    PI%parmin(9)=0.01
    PI%parmax(9)=0.04

    ! Temp factor* = Q10 = 1.2-1.6
    PI%parmin(10)=1.2 ! 1.6
    PI%parmax(10)=1.6 ! 2.4

    ! Canopy Efficiency
    ! set to parmin=1 for FLUXCOM only
    ! e.g. for wetlands etc.
    PI%parmin(11)=-0.50
    PI%parmax(11)= 0.70

    ! GSI max labile turnover
    PI%parmin(12)=0.0001
    PI%parmax(12)=0.2

    ! Fraction to Clab*/
    PI%parmin(13)=0.01
    PI%parmax(13)=0.5

    ! GSI min temperature threshold (oC)
    PI%parmin(14)=225d0
    PI%parmax(14)=330d0

    ! GSI max temperature threshold (oC)
    PI%parmin(15)=225d0
    PI%parmax(15)=330d0

    ! GSI min photoperiod threshold (sec)
    PI%parmin(16)=3600d0 ! 21600d0 ! 6 hours
    PI%parmax(16)=3600d0*10d0 ! 18 hours

    ! LMA
    ! Kattge et al. 2011,
    PI%parmin(17)=10d0
    PI%parmax(17)=200d0

    ! GSI max photoperiod threshold (sec)
    PI%parmin(24)=3600d0 !21600d0 ! 6 hours
    PI%parmax(24)=64800d0 ! 18 hours

    ! GSI min VPD threshold (Pa)
    PI%parmin(25)=1d0
    PI%parmax(25)=5500d0

    ! GSI max VPD threshold (Pa)
    PI%parmin(26)=1d0
    PI%parmax(26)=5500d0

    ! Root litter turnover rate (DBio)
    PI%parmin(27)=0.0001
    PI%parmax(27)=0.01

    ! Wood litter turnover rate (DBio)
    PI%parmin(28)=0.0001
    PI%parmax(28)=0.01

    ! Efficiency of substrate uptake by microbes (DBio)
    ! original value from Xenakis & Williams (2014)
    PI%parmin(29)=0.62-(0.62*0.5)
    PI%parmax(29)=0.62+(0.62*0.5)

    ! 2nd order rate constant for microbial uptake (DBio)
    PI%parmin(30)=0.1129056-(0.1129056*0.5)
    PI%parmax(30)=0.1129056+(0.1129056*0.5)

    ! Maximum microbial death rate (DBio)
    ! original value from Xenakis & Williams (2014)
    PI%parmin(31)=0.24-(0.24*0.5)
    PI%parmax(31)=0.24+(0.24*0.5)

    ! Inhibition constant for microbial death (DBio)
    ! original value from Xenakis & Williams (2014)
    PI%parmin(32)=0.213-(0.213*0.5)
    PI%parmax(32)=0.213+(0.213*0.5)

    ! Maintenance respiration coefficient (DBio)
    PI%parmin(33)=0.45
    PI%parmax(33)=0.55

    ! GSI senstivity for leaf senescence
    PI%parmin(34)=0.96 !0.95
    PI%parmax(34)=1.00

    ! Decomposition of som_slow (DBio)
    ! original value from Xenakis & Williams (2014)
    PI%parmin(35)=0.455868-(0.455868*0.5)
    PI%parmax(35)=0.455868+(0.455868*0.5)

    ! Initial microbial activity (DBio)
    PI%parmin(36)=0.01
    PI%parmax(36)=0.1

    ! Foliar lignin fraction (DBio)
    PI%parmin(41)=0.01
    PI%parmax(41)=0.40

    ! Fine root lignin fraction (DBio)
    PI%parmin(42)=0.01
    PI%parmax(42)=0.40

    ! Wood lignin fraction (DBio)
    ! lignin fractions based on Cornwell et al (2009) GCB
    PI%parmin(43)=0.15
    PI%parmax(43)=0.40

    ! critical GPP for LAI increase (fraction)
    PI%parmin(44)=1e-10
    PI%parmax(44)=0.30

    ! fraction of Cwood which is branch
    PI%parmin(45)=0.05
    PI%parmax(45)=0.40 !0.65

    ! fraction of Cwood which is coarse root
    PI%parmin(46)=0.15
    PI%parmax(46)=0.30 !0.45

    ! Inhibition constant for C dependant microbial activity (DBio)
    PI%parmin(51)=75.0
    PI%parmax(51)=250.0

    ! GSI - have I just left a growing state (>1)
    PI%parmin(52)=0.50
    PI%parmax(52)=1.50

    ! GSI - initial GSI value
    PI%parmin(53)=1.0
    PI%parmax(53)=2.0

    !
    ! INITIAL VALUES DECLARED HERE
    !

    ! initial conditions updated based on
    ! Adegbidi et al (2005)

    ! C labile
    PI%parmin(18)=1d0
    PI%parmax(18)=1000d0

    ! C foliar
    PI%parmin(19)=1d0
    PI%parmax(19)=1000d0

    ! C roots
    PI%parmin(20)=1d0
    PI%parmax(20)=1000d0

    ! C_wood
    PI%parmin(21)=1d0
    PI%parmax(21)=50000d0

    ! C som (fast)
    PI%parmin(22)=1d0
    PI%parmax(22)=100d0

    ! C_som (slow)
    PI%parmin(23)=100d0
    PI%parmax(23)=200000d0

    ! C litter (foliar)
    PI%parmin(37)=1d0
    PI%parmax(37)=5000d0

    ! C litter (roots)
    PI%parmin(38)=1d0
    PI%parmax(38)=5000d0

    ! C litter (wood)
    PI%parmin(39)=1d0
    PI%parmax(39)=10000d0

    ! C microbial
    PI%parmin(40)=1d0
    PI%parmax(40)=100d0

    !
    ! Replanting pools values
    !

    ! C labile
    PI%parmin(47)=1.0
    PI%parmax(47)=500.0

    ! C foliar
    PI%parmin(48)=1.0
    PI%parmax(48)=500.0

    ! C roots
    PI%parmin(49)=1.0
    PI%parmax(49)=500.0

    ! C_wood
    PI%parmin(50)=1.0
    PI%parmax(50)=1000.0

  end subroutine pars_info

  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
