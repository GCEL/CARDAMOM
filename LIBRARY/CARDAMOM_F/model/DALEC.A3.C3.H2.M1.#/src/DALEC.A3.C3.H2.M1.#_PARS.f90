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

    ! crop model will be ran and therefore needs specific parameters to be
    ! called

    !
    ! declare parameters
    !

    ! Decomposition rate (frac/day)
    PI%parmin(1) = 0.00001d0
    PI%parmax(1) = 0.01d0

    ! Fraction of GPP to autotrophic pool
    PI%parmin(2) = 0.2d0
    PI%parmax(2) = 0.8d0

    ! max development rate (day-1) DS (0->1)
    PI%parmin(3) = 0.01d0 
    PI%parmax(3) = 0.05d0
    ! max development rate (day-1) DS (1->2)
    PI%parmin(4) = 0.01d0 
    PI%parmax(4) = 0.05d0

    ! turnover rate foliage (frac/day)
    PI%parmin(5) = 0.002737851d0 ! 365.25 days
    PI%parmax(5) = 0.1428571d0   ! 7 days 

    ! TOR wood* - 1% loss per year value (day-1)
    PI%parmin(6) = 0.002737851d0 ! 365.25 days
    PI%parmax(6) = 0.1428571d0   ! 7 days 0.24d0
    ! maximum rate of foliar turnover (day-1) due to self-shading
    PI%parmin(7) = 0.02380952d0  ! 42 days
    PI%parmax(7) = 0.07142857d0  ! 17 days

    ! effective vernalisation days when plants are 50 % vernalised
    PI%parmin(8) = 12d0 
    PI%parmax(8) = 32d0

    ! Turnover of litter (fraction; temperature adjusted)
    PI%parmin(9) = 0.0001141d0 ! 24   years at 0oC
    PI%parmax(9) = 0.02d0      ! 0.13 years at 0oC

    ! Turnover of som to Rhet (fraction; temperature adjusted)
    PI%parmin(10) = 1.368925d-06   ! 2000 years at 0oC
    PI%parmax(10) = 9.126169d-05   !   30 years at 0oC !0.0001368926d0 !   20 years at 0oC

    ! log10 avg foliar N (gN.m-2)
    PI%parmin(11) = -0.2218487d0 
    PI%parmax(11) = 0.6382028d0

    ! sow day
    PI%parmin(12) = 100d0 
    PI%parmax(12) = 3650d0

    ! Respiratory cost of labile transfer (per gC.m-2 labile)
    ! Nominal assumption here is 21 %, i.e. this is otherwise 
    ! referred to as the growth respiration
    PI%parmin(13) = 0.15d0 
    PI%parmax(13) = 0.25d0

    ! phenological heat units required for emergence
    PI%parmin(14) = 100d0 
    PI%parmax(14) = 150d0

    ! harvest day
    PI%parmin(15) = 15d0 
    PI%parmax(15) = 365d0

    ! plough day
    ! NOT CURRENTLY IN USE
    PI%parmin(16) = 365.25d0 
    PI%parmax(16) = 365.25d0*4d0

    ! LMA
    PI%parmin(17) = 10d0 
    PI%parmax(17) = 100d0

    !
    ! NOTE number order not consistent
    !

    ! minimum temperature for development (oC)
    PI%parmin(26) = (-1d0+273.15d0) 
    PI%parmax(26) = (10d0+273.15d0)  ! -10,10
    ! maximum temperature for development (oC)
    PI%parmin(27) = (10d0+273.15d0) 
    PI%parmax(27) = (36d0+273.15d0)   ! 20,42
    ! optimum temperature for development (oC)
    PI%parmin(28) = (10d0+273.15d0) 
    PI%parmax(28) = (30d0+273.15d0)   ! 10,35

    ! minimum temperature for vernalisation (oC)
    PI%parmin(29) = (-5.3d0+273.15d0)
    PI%parmax(29) = (-0.3d0+273.15d0)   ! -15,10
    ! maximum temperature for vernalisation (oC)
    PI%parmin(30) = (12.7d0+273.15d0) 
    PI%parmax(30) = (18.7d0+273.15d0)    ! 5,30
    ! optimum temperature for vernalisation (oC)
    PI%parmin(31) = (2.9d0+273.15d0) 
    PI%parmax(31) = (6.9d0+273.15d0)   ! -5,15

    ! critical photoperiod for development (hrs)
    PI%parmin(32) = 6d0 
    PI%parmax(32) = 13d0
    ! photoperiod sensitivity
    PI%parmin(33) = 0.10d0 
    PI%parmax(33) = 0.35d0

    ! turnover rate of labile
    PI%parmin(34) = 1d-5
    PI%parmax(34) = 0.15d0
    ! turnover rate of autotrophic pool
    PI%parmin(35) = 0.032d0 ! 31.25 days 
    PI%parmax(35) = 0.99d0  ! ~ 1 day

    ! BUCKET - root biomass needed to reach 50 % of max depth
    PI%parmin(36) = 10d0
    PI%parmax(36) = 150d0

    ! BUCKET - maximum rooting depth
    PI%parmin(37) = 0.35d0
    PI%parmax(37) = 5d0 !20.0

!    ! Canopy NIR reflectance
!    PI%parmin(39) = 0.05d0
!    PI%parmax(39) = 0.60d0
!    ! Canopy PAR reflectance
!    PI%parmin(40) = 0.05d0
!    PI%parmax(40) = 0.60d0
!    ! Canopy NIR transmittance
!    PI%parmin(41) = 0.05d0
!    PI%parmax(41) = 0.60d0
!    ! Canopy PAR transmittance
!    PI%parmin(42) = 0.05d0
!    PI%parmax(42) = 0.60d0
!    ! Soil reflectance to near infrared radiation
!    PI%parmin(43) = 0.01d0
!    PI%parmax(43) = 0.20d0
!    ! Soil reflectance to photosynthetically active radiation
!    PI%parmin(44) = 0.01d0
!    PI%parmax(44) = 0.20d0

    !
    ! INITIAL VALUES DECLARED HERE
    !

    ! C labile (gC.m-2)
    PI%parmin(18) = 1d0 
    PI%parmax(18) = 10d0
    ! C foliar (gC.m-2)
    PI%parmin(19) = 1d0 
    PI%parmax(19) = 5d0
    ! C roots (gC.m-2)
    PI%parmin(20) = 1d0 
    PI%parmax(20) = 5d0
    ! C_wood (gC.m-2)
    PI%parmin(21) = 1d0 
    PI%parmax(21) = 5d0
    ! C litter (gC.m-2)
    PI%parmin(22) = 1d0 
    PI%parmax(22) = 10d0
    ! C_som (gC.m-2)
    PI%parmin(23) = 100d0 
    PI%parmax(23) = 200000d0
    ! C autotrophic pool (gC.m-2)
    PI%parmin(24) = 0.1d0 
    PI%parmax(24) = 5d0
    ! C storage organ (gC.m-2)
    PI%parmin(25) = 0.1d0 
    PI%parmax(25) = 1d0

    ! Soil water fraction (m3/m3)
    PI%parmin(38) = 0.1d0 
    PI%parmax(38) = 0.9d0

    ! Read in the crop type specific development file
    call crop_development_parameters(PI%stock_seed_labile,PI%DS_shoot &
                                    ,PI%DS_root,PI%fol_frac,PI%stem_frac &
                                    ,PI%root_frac,PI%DS_LRLV,PI%LRLV &
                                    ,PI%DS_LRRT,PI%LRRT)

  end subroutine pars_info
  !
  !------------------------------------------------------------------
  !
  subroutine crop_development_parameters(stock_seed_labile,DS_shoot,DS_root,fol_frac &
                                        ,stem_frac,root_frac,DS_LRLV,LRLV,DS_LRRT,LRRT)

    ! subroutine reads in the fixed crop development files which are linked the
    ! the development state of the crops. The development model varies between
    ! which species. e.g. winter wheat and barley, spring wheat and barley
    ! NOTE: duplicate function in the R_interface.f90

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
    integer        :: columns, i, rows, input_crops_unit, ios
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
