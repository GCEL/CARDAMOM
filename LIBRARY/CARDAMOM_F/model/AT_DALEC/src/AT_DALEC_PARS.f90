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
    ! These could or possibly should go into an alternate file which can be read in.
    ! This may improve the usability when it comes to reading these information
    ! in for different PFTs

    implicit none

    ! default case is to use global parameters for DALEC_CDEA, however there are
    ! model / PFT specifics available to reduce uncertainties
    ! contains 6 fields with min max log for par and par

    if (DATAin%ID == 1 .or. DATAin%ID == 2) then
       ! DALEC_CDEA and DALEC_BUCKET require alternate source
       write(*,*)"Incorrect source code has been compiled"
       write(*,*)"DALEC_CDEA or DALEC_BUCKET parameters requested in source code for AT_DALEC and variants"
       stop
    else if (DATAin%ID == 3) then
       ! AT_DALEC generic
       call at_dalec_parameters
    else if (DATAin%ID == 4) then
       ! AT_DALEC with PFT specific parameters
       if (DATAin%PFT == 1) then
          ! crops
          call at_dalec_pft1
          call crop_development_parameters
       else if (DATAin%PFT == 2) then
          ! short grass

       else if (DATAin%PFT == 3) then
          ! evergreen forest

       else if (DATAin%PFT == 5) then
          ! deciduous forest

       else if (DATAin%PFT == 11) then
          ! sparse vegetation

       else if (DATAin%PFT == 13) then
          ! Bog

       else if (DATAin%PFT == 18) then
          ! mixed forest

       else if (DATAin%PFT == 19) then
          ! interupted forest

       else
          !  missing values
          write(*,*)"PFT for AT_DALEC with PFT specific parameters not found"
          write(*,*)"PFT asked for =",DATAin%PFT
          stop
       end if! pft condition
    else if (DATAin%ID == 5) then
       ! Forest rotation
       write(*,*)"Holding space for forest rotation model" ; stop
    else
       write(*,*) "Cock model cannot be found"
       stop
    endif ! MODEL type condition

  end subroutine pars_info

  !
  !------------------------------------------------------------------
  !
  subroutine at_dalec_parameters

    ! Subroutine reads specific parameter ranges for the
    ! generic AT_DALEC model

    use MCMCOPT, only: PI

    implicit none

    !
    ! declare parameters
    !

!    PI%npars=22;

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
    PI%parmax(4)=1.0

    ! Leaf Lifespan
    ! Wright et al. 2004
    PI%parmin(5)=1.001
    PI%parmax(5)=8.0

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

    ! Bday
    PI%parmin(11)=365.25
    PI%parmax(11)=365.25*4

    ! Fraction to Clab*/
    PI%parmin(12)=0.01
    PI%parmax(12)=0.5

    ! Clab Release period
    PI%parmin(13)=10.
    PI%parmax(13)=100.

    ! Fday
    PI%parmin(14)=365.25
    PI%parmax(14)=365.25*4

    ! Leaf fall period
    PI%parmin(15)=20.
    PI%parmax(15)=150.

    ! LMA
    ! Kattge et al. 2011
    PI%parmin(16)=10.
    PI%parmax(16)=200.

    !
    ! INITIAL VALUES DECLARED HERE
    !

    ! C labile
    PI%parmin(17)=20.0
    PI%parmax(17)=2000.0

    ! C foliar
    PI%parmin(18)=20.0
    PI%parmax(18)=2000.0

    ! C roots
    PI%parmin(19)=20.0
    PI%parmax(19)=2000.0

    ! C_wood
    PI%parmin(20)=100.0
    PI%parmax(20)=100000.0

    ! C litter
    PI%parmin(21)=20.0
    PI%parmax(21)=2000.0

    ! C_som
    PI%parmin(22)=100.0
    PI%parmax(22)=200000.0

  end subroutine at_dalec_parameters

  !
  !------------------------------------------------------------------
  !
  subroutine at_dalec_pft1

    ! Subroutine reads specific parameter ranges for the
    ! generic AT_DALEC model

    use MCMCOPT, only: PI

    implicit none

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
    PI%parmin(8)=12d0 ; PI%parmax(8)=32d0

    ! mineralisation rate of SOM (hr-1)
    PI%parmin(9)=1e-6 ; PI%parmax(9)=1e-3
    ! mineralisation rate of litter (hr-1)
    PI%parmin(10)=1e-5 ; PI%parmax(10)=1e-2

    ! turnover rate of autotrophic pool
    PI%parmin(11)=0.001 ; PI%parmax(11)=0.1

    ! sow day
    PI%parmin(12)=240d0 ; PI%parmax(12)=365.25+150d0

    ! respiratory cost of labile transfer (per gC.m-2 labile)
    PI%parmin(13)=0.05 ; PI%parmax(13)=0.4

    ! phenological heat units required for emergence
    PI%parmin(14)=100d0 ; PI%parmax(14)=150d0

    ! harvest day
    PI%parmin(15)=150d0 ; PI%parmax(15)=280d0
    ! Plough day
    PI%parmin(16)=365.25 ; PI%parmax(16)=365.25*4d0

    ! LMA
    PI%parmin(17)=10d0 ; PI%parmax(17)=50d0

    !
    ! NOTE number order not consistent
    !

    ! minimum temperature for development (oC)
!    PI%parmin(26)=(-1d0+273.15) ; PI%parmax(26)=(8d0+273.15)  ! -10,10
!    ! maximum temperature for development (oC)
!    PI%parmin(27)=(28d0+273.15) ; PI%parmax(27)=(36d0+273.15)   ! 20,42
!    ! optimum temperature for development (oC)
!    PI%parmin(28)=(23d0+273.15) ; PI%parmax(28)=(30d0+273.15)   ! 10,35
    ! minimum temperature for development (oC)
    PI%parmin(26)=(-1d0+273.15) ; PI%parmax(26)=(10d0+273.15)  ! -10,10
    ! maximum temperature for development (oC)
    PI%parmin(27)=(10d0+273.15) ; PI%parmax(27)=(36d0+273.15)   ! 20,42
    ! optimum temperature for development (oC)
    PI%parmin(28)=(10d0+273.15) ; PI%parmax(28)=(30d0+273.15)   ! 10,35

    ! minimum temperature for vernalisation (oC)
    PI%parmin(29)=(-5.3+273.15) ; PI%parmax(29)=(-0.3+273.15)   ! -15,10
    ! maximum temperature for vernalisation (oC)
    PI%parmin(30)=(12.7+273.15) ; PI%parmax(30)=(18.7+273.15)    ! 5,30
    ! optimum temperature for vernalisation (oC)
    PI%parmin(31)=(2.9+273.15) ; PI%parmax(31)=(6.9+273.15)   ! -5,15

    ! critical photoperiod for development (hrs)
    PI%parmin(32)=6d0 ; PI%parmax(32)=12d0
    ! photoperiod sensitivity
    PI%parmin(33)=0.10 ; PI%parmax(33)=0.35

    ! turnover rate of labile
    PI%parmin(34)=1e-5  ; PI%parmax(34)=0.1

    !
    ! INITIAL VALUES (gC.m-2) DECLARED HERE
    !

    ! C labile
    PI%parmin(18)=0.1 ; PI%parmax(18)=10d0
    ! C foliar
    PI%parmin(19)=0.1 ; PI%parmax(19)=5d0
    ! C roots
    PI%parmin(20)=0.1 ; PI%parmax(20)=5d0
    ! C_wood
    PI%parmin(21)=0.1 ; PI%parmax(21)=5d0
    ! C litter
    PI%parmin(22)=0.1 ; PI%parmax(22)=10
    ! C_som
    PI%parmin(23)=100.0 ; PI%parmax(23)=100000.0
    ! C autotrophic pool
    PI%parmin(24)=0.1 ; PI%parmax(24)=5d0
    ! C storage organ
    PI%parmin(25)=0.1 ; PI%parmax(25)=1d0

  end subroutine at_dalec_pft1
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine crop_development_parameters

    use MCMCOPT, only: PI

    ! subroutine reads in the fixed crop development files which are linked the
    ! the development state of the crops. The development model varies between
    ! which species. e.g. winter wheat and barley, spring wheat and barley

    implicit none

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
    read(unit=input_crops_unit,fmt=*) variables,PI%stock_seed_labile,variables,variables

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
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
end module MODEL_PARAMETERS
