module MODEL_PARAMETERS

  implicit none

  ! make all private
  private

  ! specify explicitly the public
  public :: pars_info

  contains
  
  !
  !------------------------------------------------------------------
  !
  subroutine pars_info(PI)
    use MCMCOPT, only: parameter_info
    use cardamom_structures, only: DATAin

    ! Subroutine contains a list of parameter ranges for the model.
    ! These could or possibly should go into an alternate file which can be read in.
    ! This may improve the usability when it comes to reading these information
    ! in for different PFTs

    implicit none

    ! declare inputs
    type ( parameter_info ), intent(inout) :: PI

!    PI%npars=29 ! dont forget to change in cardamom_io.f90

    if (DATAin%PFT == 1) then
       ! crop model will be ran and therefore needs specific parameters to be
       ! called
       call crop_parameters(PI)
       call crop_development_parameters(PI)
 
    else 
       ! generic model

       ! 
       ! declare parameters
       ! 

       ! Decomposition rate
       PI%parmin(1)=0.00001
       PI%parmax(1)=0.01

       ! Fraction of GPP respired
       PI%parmin(2)=0.3
       PI%parmax(2)=0.7

       ! GSI sensitivity for leaf growth
       PI%parmin(3)=1.00
       PI%parmax(3)=1.025 !1.05

       ! Fraction of (1-fgpp) to roots*/
       PI%parmin(4)=0.01
       PI%parmax(4)=1.0

       ! GSI max leaf turnover
       PI%parmin(5)=0.000001
       PI%parmax(5)=0.2

       ! TOR wood* - 1% loss per year value
       PI%parmin(6)=0.00001
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

       ! log10 avg foliar N (gN.m-2) ! Canopy Efficiency
       ! set to parmin=1 for FLUXCOM only
       ! e.g. for wetlands etc.
       PI%parmin(11)=10d0 ! -0.50
       PI%parmax(11)=100d0! 1d0  

       ! GSI max labile turnover
       PI%parmin(12)=0.000001
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
       PI%parmin(16)=3600d0 !21600d0 ! 6 hours
       PI%parmax(16)=3600d0*10d0 ! 64800d0 ! 18 hours

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

       ! critical GPP for LAI increase (fraction)
       PI%parmin(27)=1e-10
       PI%parmax(27)=0.30

       ! fraction of Cwood which is branch
       PI%parmin(28)=0.05
       PI%parmax(28)=0.65

       ! fraction of Cwood which is coarse root
       PI%parmin(29)=0.15
       PI%parmax(29)=0.45

       ! GSI senstivity for leaf senescence 
       PI%parmin(34)=0.96
       PI%parmax(34)=1.00

       ! GSI - have I just left a growing state (>1)
       PI%parmin(35)=0.50
       PI%parmax(35)=1.5

       ! GSI - initial GSI value 
       PI%parmin(36)=1.0
       PI%parmax(36)=2.0

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
       PI%parmax(21)=20000d0

       ! C litter
       PI%parmin(22)=1d0
       PI%parmax(22)=10000d0

       ! C_som
       PI%parmin(23)=100d0
       PI%parmax(23)=200000d0

       !
       ! Replanting pools values 
       !

       ! C labile
       PI%parmin(30)=1.0
       PI%parmax(30)=1000.0

       ! C foliar
       PI%parmin(31)=1.0
       PI%parmax(31)=1000.0

       ! C roots
       PI%parmin(32)=1.0
       PI%parmax(32)=1000.0

       ! C_wood
       PI%parmin(33)=1.0
       PI%parmax(33)=1000.0

    endif ! crop / default split 
  end subroutine pars_info
  !
  !------------------------------------------------------------------
  !
  subroutine crop_parameters(PI)

    ! Subroutine reads specific parameter ranges for the 
    ! generic AT_DALEC model

    use MCMCOPT, only: parameter_info

    implicit none

    ! declare inputs
    type ( parameter_info ), intent(inout) :: PI

    ! 
    ! declare parameters
    ! 

!    PI%npars=34;

    ! Decomposition rate (frac/hr)
    PI%parmin(1)=0.000001 ; PI%parmax(1)=0.0001

    ! Fraction of GPP to autotrophic pool
    PI%parmin(2)=0.3 ; PI%parmax(2)=0.7

    ! max development rate (day-1) DS (0->1)
    PI%parmin(3)=0.030 ; PI%parmax(3)=0.050
    ! max development rate (day-1) DS (1->2)
    PI%parmin(4)=0.030 ; PI%parmax(4)=0.050

    ! turnover rate foliage (frac/hour)
    PI%parmin(5)=1.0e-4 ; PI%parmax(5)=0.01

    ! TOR wood* - 1% loss per year value (day-1)
    PI%parmin(6)=1e-4 ; PI%parmax(6)=0.01
    ! maximum rate of foliar turnover (hr-1) due to self-shading
    PI%parmin(7)=1e-5 ; PI%parmax(7)=0.01

    ! effective vernalisation days when plants are 50 % vernalised
    PI%parmin(8)=12.0 ; PI%parmax(8)=32.0

    ! mineralisation rate of SOM (hr-1)
    PI%parmin(9)=1e-6 ; PI%parmax(9)=1e-3
    ! mineralisation rate of litter (hr-1)
    PI%parmin(10)=1e-5 ; PI%parmax(10)=1e-2

    ! log10 avg foliar N (gN.m-2) ! Canopy Efficiency
    ! set to parmin=1 for FLUXCOM only
    ! e.g. for wetlands etc.
    PI%parmin(11)=10.0  ; PI%parmax(11)=100.0

    ! sow day
    PI%parmin(12)=240.0 ; PI%parmax(12)=365.25+150.0

    ! respiratory cost of labile transfer (per gC.m-2 labile)
    PI%parmin(13)=0.05 ; PI%parmax(13)=0.4

    ! phenological heat units required for emergence
    PI%parmin(14)=100.0 ; PI%parmax(14)=150.0

    ! harvest day
    PI%parmin(15)=150.0 ; PI%parmax(15)=280.0
    ! Plough day
    PI%parmin(16)=365.25 ; PI%parmax(16)=365.25*4.0

    ! LMA
    PI%parmin(17)=10.0 ; PI%parmax(17)=50.0

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
    PI%parmin(34)=1e-5  ; PI%parmax(34)=0.1
    ! turnover rate of autotrophic pool
    PI%parmin(35)=0.001 ; PI%parmax(35)=0.1

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
  subroutine crop_development_parameters(PI)

    use MCMCOPT, only: parameter_info

    ! subroutine reads in the fixed crop development files which are linked the
    ! the development state of the crops. The development model varies between
    ! which species. e.g. winter wheat and barley, spring wheat and barley

    implicit none

    ! declare inputs
    type ( parameter_info ), intent(inout) :: PI  

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
    read(unit=input_crops_unit,fmt=*)variables,PI%stock_seed_labile,variables,variables

    ! read in C partitioning/fraction data and corresponding developmental
    ! stages (DS)
    ! shoot
    read(unit=input_crops_unit,fmt=*) variables
    read(unit=input_crops_unit,fmt=*) rows , columns
    allocate( PI%DS_shoot(rows) , PI%fol_frac(rows) , PI%stem_frac(rows)  )
    do i = 1 , rows
      read(unit=input_crops_unit,fmt=*) PI%DS_shoot(i), PI%fol_frac(i), PI%stem_frac(i)
    enddo

    ! root
    read(unit=input_crops_unit,fmt=*) variables
    read(unit=input_crops_unit,fmt=*) rows , columns
    allocate( PI%DS_root(rows) , PI%root_frac(rows) )
    do i = 1 , rows
      read(unit=input_crops_unit,fmt=*) PI%DS_root(i), PI%root_frac(i)
    enddo

    ! loss rates of leaves and roots
    ! leaves
    read(unit=input_crops_unit,fmt=*) variables
    read(unit=input_crops_unit,fmt=*) rows , columns
    allocate( PI%DS_LRLV(rows) , PI%LRLV(rows) )
    do i = 1 , rows
      read(unit=input_crops_unit,fmt=*) PI%DS_LRLV(i), PI%LRLV(i)
    enddo

    ! roots
    read(unit=input_crops_unit,fmt=*) variables
    read(unit=input_crops_unit,fmt=*) rows , columns
    allocate( PI%DS_LRRT(rows) , PI%LRRT(rows) )
    do i = 1 , rows
      read(unit=input_crops_unit,fmt=*) PI%DS_LRRT(i), PI%LRRT(i)
    enddo

    ! rewind and close
    rewind(input_crops_unit) ; close(input_crops_unit)

  end subroutine crop_development_parameters
  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
