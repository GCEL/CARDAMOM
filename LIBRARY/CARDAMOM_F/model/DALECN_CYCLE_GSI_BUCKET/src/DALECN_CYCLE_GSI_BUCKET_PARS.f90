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

!    PI%npars=29 ! dont forget to change in cardamom_io.f90

    if (DATAin%PFT == 1) then
       ! crop model will be ran and therefore needs specific parameters to be
       ! called
       call crop_parameters
       call crop_development_parameters(PI%stock_seed_labile,PI%DS_shoot &
                                       ,PI%DS_root,PI%fol_frac,PI%stem_frac &
                                       ,PI%root_frac,PI%DS_LRLV,PI%LRLV &
                                       ,PI%DS_LRRT,PI%LRRT)

    else

       ! generic model

       !
       ! declare parameters
       !

       ! Decomposition rate
       PI%parmin(1)=0.0001 ! 27 years
       PI%parmax(1)=0.01   ! 0.27 years

       ! CN_Root (Kattge et al., 2011)
       PI%parmin(2)=13.85
       PI%parmax(2)=192.31

       ! GSI sensitivity for leaf growth
       PI%parmin(3)=1.00
       PI%parmax(3)=1.03 !1.05

       ! fraction of labile to root
       PI%parmin(4)=3e-5 ! 1e-4
       PI%parmax(4)=0.1

       ! GSI max leaf turnover
       PI%parmin(5)=0.00027 ! 10 years
       PI%parmax(5)=0.10    ! 10 days

       ! TOR wood* - 1% loss per year value
       PI%parmin(6)=0.000009 ! 300 years
       PI%parmax(6)=0.001    ! 2.73 years

       ! TOR roots
       PI%parmin(7)=0.00025  ! 11 years
       PI%parmax(7)=0.01     ! 0.27 years

       ! TOR litter
       PI%parmin(8)=0.0001 ! 27 years
       PI%parmax(8)=0.01   ! 0.27 years

       ! TOR SOM
       PI%parmin(9)=1.368925e-6 ! 2000 years !0.000001 ! ~ 2737 years
       PI%parmax(9)=0.0005      ! ~ 5 years

       ! Temp factor* = Q10 = 1.2-1.6
       PI%parmin(10)=0.018
       PI%parmax(10)=0.08

       ! log10 avg foliar N (gN.m-2)
       ! set to parmin=1 for FLUXCOM only
       ! e.g. for wetlands etc.
       PI%parmin(11)=-0.50
       PI%parmax(11)= 0.70

       ! max labile turnover to foliage
       PI%parmin(12)=3e-5 ! 1e-3
       PI%parmax(12)=0.20 !

       ! fraction of labile to wood*/
       PI%parmin(13)=1e-4 !
       PI%parmax(13)=0.20 ! 0.15

       ! GSI min temperature threshold (oC)
       PI%parmin(14)=225d0
       PI%parmax(14)=330d0

       ! GSI max temperature threshold (oC)
       PI%parmin(15)=225d0
       PI%parmax(15)=330d0

       ! GSI min photoperiod threshold (sec)
       PI%parmin(16)=3600d0 !21600d0 ! 6 hours
       PI%parmax(16)=3600d0*10d0 ! 64800d0 ! 18 hours

       ! LMA
       ! Kattge et al. 2011,
       PI%parmin(17)=10d0
       PI%parmax(17)=200d0

       ! GSI max photoperiod threshold (sec)
       PI%parmin(24)=3600d0 !21600d0 ! 6 hours
       PI%parmax(24)=64800d0 ! 18 hours

!       ! GSI min Rtot threshold (MPa)
!       PI%parmin(25)=0.01
!       PI%parmax(25)=25.0 !50.0
!
!       ! GSI max Rtot threshold (MPa)
!       PI%parmin(26)=0.1
!       PI%parmax(26)=75.0 !100.0

       ! GSI min deltaWP threshold (MPa)
       PI%parmin(25)=-2.0
       PI%parmax(25)= -1e-4

       ! GSI max deltaWP threshold (MPa)
       PI%parmin(26)=-2.0
       PI%parmax(26)= -1e-4

       ! CN_wood (Kattge et al., 2011)
       PI%parmin(27)=169.5
       PI%parmax(27)=909.1

!       ! CN_wood baseline value (gC/gN) for logarithmic relation with woody biomass
!       PI%parmin(27)=100
!       PI%parmax(27)=300

       ! fraction of Cwood which is branch
       PI%parmin(28)=0.05
       PI%parmax(28)=0.40 !0.65

       ! fraction of Cwood which is coarse root
       PI%parmin(29)=0.15
       PI%parmax(29)=0.30 !0.45

       ! GSI senstivity for leaf senescence
       PI%parmin(34)=0.96
       PI%parmax(34)=1.00

       ! GSI - have I just left a growing state (>1)
       PI%parmin(35)=0.50
       PI%parmax(35)=1.5

       ! GSI - initial GSI value
       PI%parmin(36)=1.0
       PI%parmax(36)=2.0

       ! Turnover rate for CWD
       PI%parmin(38)=0.0001 ! 0.00001
       PI%parmax(38)=1d0/365.25 !0.005  ! 0.01

       ! BUCKET - root biomass needed to reach 50 % of max depth
       PI%parmin(39)=50.0
       PI%parmax(39)=500.0

       ! BUCKET - maximum rooting depth
       PI%parmin(40)=0.35
       PI%parmax(40)=15.0 !20.0

       ! Reich - Leaf N linked respiration exponential coefficient
       PI%parmin(41)=1.639-0.01 !0.935 ! 1.639-0.01
       PI%parmax(41)=1.639+0.01 !1.774 ! 1.639+0.01

       ! Reich - Leaf N linked respiration intercept
       PI%parmin(42)=0.01 !0.645
       PI%parmax(42)=1.25 !0.911

       ! Reich - root N linked respiration exponential coefficient
       PI%parmin(43)=1.352-0.01 !1.012 ! 1.352-0.01
       PI%parmax(43)=1.352+0.01 !1.478 ! 1.352+0.01

       ! Reich - root N linked respiration intercept
       PI%parmin(44)=0.01 !0.915
       PI%parmax(44)=1.25 !1.079

       ! Reich - wood N linked respiration exponential coefficient
       PI%parmin(45)=1.344-0.01 !1.170 ! 1.344-0.01
       PI%parmax(45)=1.344+0.01 !1.478 ! 1.344+0.01

       ! Reich - wood N linked respiration intercept
       PI%parmin(46)=0.01 !0.839
       PI%parmax(46)=1.25 !1.053

       ! Initial leaf life span (days)
       PI%parmin(47)=30.0
       PI%parmax(47)=365.25*8

       ! baseline NUE (gC/gN/m2/day-1)
       PI%parmin(48)=15.0
       PI%parmax(48)=30.0

!       ! CN_wood coefficient for increase due to C_wood (deltalCN per gC.m-2)
!       ! NOTE: values in log scale
!       PI%parmin(49)=0.01
!       PI%parmax(49)=0.25

       !
       ! INITIAL VALUES DECLARED HERE
       !

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

       ! C litter
       PI%parmin(22)=1d0
       PI%parmax(22)=5000d0

       ! C_som
       PI%parmin(23)=100d0
       PI%parmax(23)=100000d0

       ! C CWD
       PI%parmin(37)=1d0
       PI%parmax(37)=10000d0

       !
       ! Replanting pools values
       !

       ! C labile
       PI%parmin(30)=1.0
       PI%parmax(30)=500.0

       ! C foliar
       PI%parmin(31)=1.0
       PI%parmax(31)=500.0

       ! C roots
       PI%parmin(32)=1.0
       PI%parmax(32)=500.0

       ! C_wood
       PI%parmin(33)=1.0
       PI%parmax(33)=1000.0

    endif ! crop / default split

  end subroutine pars_info
  !
  !------------------------------------------------------------------
  !
  subroutine crop_parameters

    ! Subroutine reads specific parameter ranges for the
    ! generic AT_DALEC model

    implicit none

    ! declare inputs
    type ( parameter_info ), intent(inout) :: PI

    !
    ! declare parameters
    !

!    PI%npars=34;

    ! Decomposition rate (frac/day)
    PI%parmin(1)=0.000001*24.0 ; PI%parmax(1)=0.0001*24.0

    ! Fraction of GPP to autotrophic pool
    PI%parmin(2)=0.2 ; PI%parmax(2)=0.7

    ! max development rate (day-1) DS (0->1)
    PI%parmin(3)=0.030 ; PI%parmax(3)=0.050
    ! max development rate (day-1) DS (1->2)
    PI%parmin(4)=0.030 ; PI%parmax(4)=0.050

    ! turnover rate foliage (frac/day)
    PI%parmin(5)=1.0e-4*24.0 ; PI%parmax(5)=0.01*24.0

    ! TOR wood* - 1% loss per year value (day-1)
    PI%parmin(6)=1e-4*24.0 ; PI%parmax(6)=0.01*24.0
    ! maximum rate of foliar turnover (hr-1) due to self-shading
    PI%parmin(7)=1e-5*24.0 ; PI%parmax(7)=0.01*24.0

    ! effective vernalisation days when plants are 50 % vernalised
    PI%parmin(8)=12.0 ; PI%parmax(8)=32.0

    ! mineralisation rate of litter (hr-1)
    PI%parmin(9)=1e-5*24.0 ; PI%parmax(9)=1e-2*24.0
    ! mineralisation rate of SOM (hr-1)
    PI%parmin(10)=1e-6*24.0 ; PI%parmax(10)=1e-3*24.0

    ! log10 avg foliar N (gN.m-2)
    ! set to parmin=1 for FLUXCOM only
    ! e.g. for wetlands etc.
    PI%parmin(11)=-0.50 ; PI%parmax(11)=0.7

    ! sow day
    PI%parmin(12)=115 ; PI%parmax(12)=350

    ! respiratory cost of labile transfer (per gC.m-2 labile)
    PI%parmin(13)=0.05 ; PI%parmax(13)=0.4

    ! phenological heat units required for emergence
    PI%parmin(14)=100.0 ; PI%parmax(14)=150.0

    ! harvest day
    PI%parmin(15)=15 ; PI%parmax(15)=350
    ! Plough day
    PI%parmin(16)=365.25 ; PI%parmax(16)=365.25*4.0

    ! LMA
    PI%parmin(17)=10.0 ; PI%parmax(17)=100.0

    !
    ! NOTE number order not consistent
    !

    ! minimum temperature for development (oC)
    PI%parmin(26)=(-1.0+273.15) ; PI%parmax(26)=(10.0+273.15)  ! -10,10
    ! maximum temperature for development (oC)
    PI%parmin(27)=(10.0+273.15) ; PI%parmax(27)=(36.0+273.15)   ! 20,42
    ! optimum temperature for development (oC)
    PI%parmin(28)=(10.0+273.15) ; PI%parmax(28)=(30.0+273.15)   ! 10,35

    ! minimum temperature for vernalisation (oC)
    PI%parmin(29)=(-5.3+273.15) ; PI%parmax(29)=(-0.3+273.15)   ! -15,10
    ! maximum temperature for vernalisation (oC)
    PI%parmin(30)=(12.7+273.15) ; PI%parmax(30)=(18.7+273.15)    ! 5,30
    ! optimum temperature for vernalisation (oC)
    PI%parmin(31)=(2.9+273.15) ; PI%parmax(31)=(6.9+273.15)   ! -5,15

    ! critical photoperiod for development (hrs)
    PI%parmin(32)=6.0 ; PI%parmax(32)=12.0
    ! photoperiod sensitivity
    PI%parmin(33)=0.10 ; PI%parmax(33)=0.35

    ! turnover rate of labile
    PI%parmin(34)=1e-5*24.0  ; PI%parmax(34)=0.00625*24.0
    ! turnover rate of autotrophic pool
    PI%parmin(35)=0.001*24.0 ; PI%parmax(35)=0.07*24.0

    ! BUCKET - root biomass needed to reach 50 % of max depth
    PI%parmin(36)=10.0
    PI%parmax(36)=500.0

    ! BUCKET - maximum rooting depth
    PI%parmin(37)=0.35
    PI%parmax(37)=5.0 !20.0

    !
    ! INITIAL VALUES (gC.m-2) DECLARED HERE
    !

    ! C labile
    PI%parmin(18)=1.0 ; PI%parmax(18)=10.0
    ! C foliar
    PI%parmin(19)=1.0 ; PI%parmax(19)=5.0
    ! C roots
    PI%parmin(20)=1.0 ; PI%parmax(20)=5.0
    ! C_wood
    PI%parmin(21)=1.0 ; PI%parmax(21)=5.0
    ! C litter
    PI%parmin(22)=1.0 ; PI%parmax(22)=10
    ! C_som
    PI%parmin(23)=100.0 ; PI%parmax(23)=200000.0
    ! C autotrophic pool
    PI%parmin(24)=0.1 ; PI%parmax(24)=5.0
    ! C storage organ
    PI%parmin(25)=0.1 ; PI%parmax(25)=1.0

  end subroutine crop_parameters
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine crop_development_parameters(stock_seed_labile,DS_shoot,DS_root,fol_frac &
                                        ,stem_frac,root_frac,DS_LRLV,LRLV,DS_LRRT,LRRT)

    ! subroutine reads in the fixed crop development files which are linked the
    ! the development state of the crops. The development model varies between
    ! which species. e.g. winter wheat and barley, spring wheat and barley

    implicit none

    ! declare inputs
    ! crop specific variables
    double precision,intent(inout) :: stock_seed_labile
    double precision, allocatable, dimension(:),intent(inout)  :: DS_shoot, & !
                                                                   DS_root, & !
                                                                  fol_frac, & !
                                                                 stem_frac, & !
                                                                 root_frac, & !
                                                                   DS_LRLV, & !
                                                                      LRLV, & !
                                                                   DS_LRRT, & !
                                                                      LRRT

    ! local variables..
    integer                 :: columns, i, rows, input_crops_unit, ios
    character(100) :: variables,filename

    ! for the moment hard code the file name
    filename="winter_wheat_development.csv"
    input_crops_unit = 20 ; ios = 0

    ! crop development file
    open(unit = input_crops_unit, file=trim(filename),iostat=ios, status='old', action='read')

    ! ensure we are definitely at the beginning
    rewind(input_crops_unit)

    ! read in the amount of carbon available (as labile) in each seed..
    read(unit=input_crops_unit,fmt=*)variables,stock_seed_labile,variables,variables

    ! read in C partitioning/fraction data and corresponding developmental
    ! stages (DS)
    ! shoot
    read(unit=input_crops_unit,fmt=*) variables
    read(unit=input_crops_unit,fmt=*) rows , columns
    allocate( DS_shoot(rows) , fol_frac(rows) , stem_frac(rows)  )
    do i = 1 , rows
      read(unit=input_crops_unit,fmt=*) DS_shoot(i), fol_frac(i), stem_frac(i)
    enddo

    ! root
    read(unit=input_crops_unit,fmt=*) variables
    read(unit=input_crops_unit,fmt=*) rows , columns
    allocate( DS_root(rows) , root_frac(rows) )
    do i = 1 , rows
      read(unit=input_crops_unit,fmt=*) DS_root(i), root_frac(i)
    enddo

    ! loss rates of leaves and roots
    ! leaves
    read(unit=input_crops_unit,fmt=*) variables
    read(unit=input_crops_unit,fmt=*) rows , columns
    allocate( DS_LRLV(rows) , LRLV(rows) )
    do i = 1 , rows
      read(unit=input_crops_unit,fmt=*) DS_LRLV(i), LRLV(i)
    enddo

    ! roots
    read(unit=input_crops_unit,fmt=*) variables
    read(unit=input_crops_unit,fmt=*) rows , columns
    allocate( DS_LRRT(rows) , LRRT(rows) )
    do i = 1 , rows
      read(unit=input_crops_unit,fmt=*) DS_LRRT(i), LRRT(i)
    enddo

    ! rewind and close
    rewind(input_crops_unit) ; close(input_crops_unit)

  end subroutine crop_development_parameters
  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
