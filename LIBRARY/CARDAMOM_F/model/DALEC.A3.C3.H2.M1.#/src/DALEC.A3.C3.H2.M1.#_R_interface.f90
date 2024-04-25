

subroutine rdalec15(output_dim,MTT_dim,SS_dim &
                   ,met,pars &
                   ,out_var1,out_var2,out_var3,out_var4,out_var5 &
                   ,lat,nopars,nomet &
                   ,nofluxes,nopools,nodays,nos_years,deltat &
                   ,nos_iter,soil_frac_clay_in,soil_frac_sand_in &
                   ,pathlength)

  use CARBON_MODEL_MOD, only: CARBON_MODEL, wSWP_time, soil_frac_clay, &
                              soil_frac_sand, nos_soil_layers, &
                              gs_demand_supply_ratio, cica_time, &
                              gs_total_canopy, gb_total_canopy, &
                              canopy_par_MJday_time, root_depth_time, &
                              snow_storage_time, DS_time

  ! subroutine specificially deals with the calling of the fortran code model by
  ! R

  !!!!!!!!!!!
  ! Authorship contributions
  !
  ! This code is by:
  ! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
  ! See function / subroutine specific comments for exceptions and contributors
  !!!!!!!!!!!

  implicit none
  interface
    subroutine crop_development_parameters(stock_seed_labile,DS_shoot,DS_root,fol_frac &
                                          ,stem_frac,root_frac,DS_LRLV,LRLV,DS_LRRT,LRRT &
                                          ,exepath)
      implicit none
      ! declare inputs
      ! crop specific variables
      character(350),intent(in) :: exepath
      double precision :: stock_seed_labile
      double precision, allocatable, dimension(:) :: DS_shoot, & !
                                                      DS_root, & !
                                                     fol_frac, & !
                                                    stem_frac, & !
                                                    root_frac, & !
                                                      DS_LRLV, & !
                                                         LRLV, & !
                                                      DS_LRRT, & !
                                                         LRRT
      ! local variables..
      integer :: columns, i, rows, input_crops_unit, ios
      character(225) :: variables,filename
    end subroutine crop_development_parameters
  end interface

  ! declare input variables
  integer, intent(in) :: pathlength
  !character(pathlength), intent(in) :: exepath
  integer, intent(in) :: nopars         & ! number of paremeters in vector
                        ,output_dim     & !
                        ,MTT_dim        & ! number of pools mean transit time estimates
                        ,SS_dim         & ! number of pools the steady state will be output for
                        ,nos_iter       & !
                        ,nomet          & ! number of meteorological fields
                        ,nofluxes       & ! number of model fluxes
                        ,nopools        & ! number of model pools
                        ,nodays         & ! number of time steps in simulation
                        ,nos_years        ! number of years in simulation

  double precision, intent(inout) :: deltat(nodays)     ! time step in decimal days
  double precision, intent(in) :: met(nomet,nodays)   & ! met drivers, note reverse of needed
                  ,soil_frac_clay_in(nos_soil_layers) & ! clay in soil (%)
                  ,soil_frac_sand_in(nos_soil_layers) & ! sand in soil (%)
                       ,pars(nopars,nos_iter)         & ! number of parameters
                       ,lat                 ! site latitude (degrees)

  ! output declaration
  double precision, intent(out), dimension(nos_iter,nodays,output_dim) :: out_var1
  double precision, intent(out), dimension(nos_iter,MTT_dim) :: out_var2  ! Mean annual MRT (years)
  double precision, intent(out), dimension(nos_iter,SS_dim) :: out_var3  ! Steady State (gC/m2)
  double precision, intent(out), dimension(nos_iter,output_dim) :: out_var4 ! Long term mean of out_var1
  double precision, intent(out), dimension(nos_iter,nos_years,output_dim) :: out_var5 ! Mean annual of out_var1

  ! local variables
  ! vector of ecosystem pools
  character(350) :: exepath
  integer :: a, e, i, s, v, steps_per_year!, nos_years
  integer, dimension(nodays) :: pool_hak
  double precision, dimension((nodays+1),nopools) :: POOLS
  ! vector of ecosystem fluxes
  double precision, dimension(nodays,nofluxes) :: FLUXES
  double precision, dimension(nodays) :: tmp &
                                        ,lai & ! leaf area index
                                        ,GPP & ! Gross primary productivity
                                        ,NEE   ! net ecosystem exchange of CO2


  ! crop development parameters declared here. These are also found in
  ! MHMCMC_STRUCTURES PI%
  ! crop specific variables
  double precision :: stock_seed_labile
  double precision, allocatable, dimension(:)  ::  DS_shoot, & !
                                                    DS_root, & !
                                                   fol_frac, & !
                                                  stem_frac, & !
                                                  root_frac, & !
                                                    DS_LRLV, & !
                                                       LRLV, & !
                                                    DS_LRRT, & !
                                                       LRRT

  ! zero initial conditions
  lai = 0d0 ; GPP = 0d0 ; NEE = 0d0 ; POOLS = 0d0 ; FLUXES = 0d0
  out_var1 = 0d0 ; out_var2 = 0d0 ; out_var3 = 0d0

  ! update soil parameters
  soil_frac_clay(1:nos_soil_layers) = soil_frac_clay_in(1:nos_soil_layers)
  soil_frac_sand(1:nos_soil_layers) = soil_frac_sand_in(1:nos_soil_layers)

  ! generate deltat step from input data
  deltat(1) = met(1,1)
  do i = 2, nodays
     deltat(i) = met(1,i)-met(1,(i-1))
  end do
  ! number of time steps per year
  steps_per_year = nodays/nos_years

  ! Determine which crop development file are we looking for
  if (pathlength == 1) then
      exepath = "winter_wheat_development.csv"
  else 
      ! Not currently defined
      print*,"No valid crop type has been specified"
      stop
  end if 

  ! Load crop development parameters here
  ! TLS: should the exepath and pathlength be made hardcoded in the assumption that 
  ! the needed file always has the same name (i.e. copy and rename to EXE to match specific crop types)
  ! and that the file is in the working directory of the PROJECT
  call crop_development_parameters(stock_seed_labile,DS_shoot,DS_root,fol_frac &
                                  ,stem_frac,root_frac,DS_LRLV,LRLV,DS_LRRT,LRRT &
                                  ,exepath)

  ! begin iterations
  do i = 1, nos_iter
     ! call the model
     call CARBON_MODEL(1,nodays,met,pars(1:nopars,i),deltat,nodays,lat &
                      ,lai,NEE,FLUXES,POOLS,nopars,nomet,nopools,nofluxes &
                      ,GPP,stock_seed_labile,DS_shoot,DS_root,fol_frac &
                      ,stem_frac,root_frac,DS_LRLV,LRLV,DS_LRRT,LRRT)
!if (i == 1) then
!    open(unit=666,file="/home/lsmallma/out.csv", &
!         status='replace',action='readwrite' )
!write(666,*)"deltat",deltat
!    write(666,*),"GSI",FLUXES(:,14)(1:365)
!    close(666)
!endif

     !
     ! Allocate the output the our 'output' variable
     !

     ! C ecosystem fluxes (gC/m2/day)
     out_var1(i,1:nodays,1)  = FLUXES(1:nodays,1)          ! GPP (gC/m2/day)
     out_var1(i,1:nodays,2)  = FLUXES(1:nodays,3)          ! Rauto (gC/m2/day)
     out_var1(i,1:nodays,3)  = FLUXES(1:nodays,13)         ! Rhet_litter (gC/m2/day)
     out_var1(i,1:nodays,4)  = FLUXES(1:nodays,14)         ! Rhet_som (gC/m2/day)
     out_var1(i,1:nodays,5)  = FLUXES(1:nodays,21)         ! harvested material (gC/m2/day)
     out_var1(i,1:nodays,6)  = FLUXES(1:nodays,22)         ! non-yield harvested material (gC/m2/day)
     ! C internal fluxes (gC/m2/day)
     out_var1(i,1:nodays,7)  = FLUXES(1:nodays,4)          ! allocation to foliage (gC/m2/day)
     out_var1(i,1:nodays,8)  = FLUXES(1:nodays,5)          ! allocation to labile (gC/m2/day)
     out_var1(i,1:nodays,9)  = FLUXES(1:nodays,6)          ! allocation to fine roots (gC/m2/day)
     out_var1(i,1:nodays,10) = FLUXES(1:nodays,7)          ! allocation to stem (gC/m2/day)
     out_var1(i,1:nodays,11) = FLUXES(1:nodays,8)          ! labile from foliage (gC/m2/day)
     out_var1(i,1:nodays,12) = FLUXES(1:nodays,16)         ! alloc to autotrophic (gC/m2/day)
     out_var1(i,1:nodays,13) = FLUXES(1:nodays,9)          ! alloc to storage organ (gC/m2/day)
     out_var1(i,1:nodays,14) = FLUXES(1:nodays,10)         ! foliage scenesence (gC/m2/day)
     out_var1(i,1:nodays,15) = FLUXES(1:nodays,12)         ! fine root turnover (gC/m2/day)
     out_var1(i,1:nodays,16) = FLUXES(1:nodays,11)         ! stem turnover (gC /m2/day)
     out_var1(i,1:nodays,17) = FLUXES(1:nodays,15)         ! Decomp_litter (gC/m2/day)
     out_var1(i,1:nodays,18) = FLUXES(1:nodays,23)         ! Respiration from autotrophic allocation (gC.m-2.day-1)
     out_var1(i,1:nodays,19) = FLUXES(1:nodays,24)         ! Respiration from labile to foliage translocation (gC.m-2.day-1)
     out_var1(i,1:nodays,20) = FLUXES(1:nodays,25)         ! Respiration from npp to labile translocation (gC.m-2.day-1)
     out_var1(i,1:nodays,21) = FLUXES(1:nodays,26)         ! Respiration from foliage remobilisation (gC.m-2.day-1)
     out_var1(i,1:nodays,22) = FLUXES(1:nodays,27)         ! Respiration from stem remobilisation (gC.m-2.day-1)
     out_var1(i,1:nodays,23) = FLUXES(1:nodays,28)         ! Foliage extracted from harvest (gC.m-2.d-1)
     out_var1(i,1:nodays,24) = FLUXES(1:nodays,29)         ! Stem extracted from harvest (gC.m-2.d-1)
     out_var1(i,1:nodays,25) = FLUXES(1:nodays,30)         ! Dead still standing foliage extracted from harvest (gC.m-2.d-1)
     out_var1(i,1:nodays,26) = FLUXES(1:nodays,31)         ! Labile extracted from harvest (gC.m-2.d-1)
     out_var1(i,1:nodays,27) = FLUXES(1:nodays,32)         ! Foliage added to litter due to harvest (gC.m-2.d-1)
     out_var1(i,1:nodays,28) = FLUXES(1:nodays,33)         ! Stem added to litter due to harvest (gC.m-2.d-1)
     out_var1(i,1:nodays,29) = FLUXES(1:nodays,34)         ! Dead still standing foliage added to litter due to harvest (gC.m-2.d-1)
     out_var1(i,1:nodays,30) = FLUXES(1:nodays,35)         ! Autotrophic pool added to litter due to harvest (gC.m-2.d-1)
     out_var1(i,1:nodays,31) = FLUXES(1:nodays,36)         ! Labile added to litter due to harvest (gC.m-2.d-1)
     out_var1(i,1:nodays,32) = FLUXES(1:nodays,37)         ! Root added to litter due to plough (gC.m-2.d-1)
     ! C pools (gC/m2)
     out_var1(i,1:nodays,33) = POOLS(1:nodays,1)           ! labile (gC/m2)
     out_var1(i,1:nodays,34) = POOLS(1:nodays,2)           ! foliage (gC/m2)
     out_var1(i,1:nodays,35) = POOLS(1:nodays,3)           ! fine root (gC/m2)
     out_var1(i,1:nodays,36) = POOLS(1:nodays,4)           ! wood (gC/m2)
     out_var1(i,1:nodays,37) = POOLS(1:nodays,5)           ! litter (gC/m2)
     out_var1(i,1:nodays,38) = POOLS(1:nodays,6)           ! som (gC/m2)
     out_var1(i,1:nodays,39) = POOLS(1:nodays,7)           ! autotrophic (gC/m2)
     out_var1(i,1:nodays,40) = POOLS(1:nodays,9)           ! storage organ (gC/m2)
     out_var1(i,1:nodays,41) = POOLS(1:nodays,10)          ! dead but still standing foliage (gC/m2)
     ! Water cycle related
     out_var1(i,1:nodays,42) = FLUXES(1:nodays,19)         ! Evapotranspiration (kgH2O.m-2.day-1)
     out_var1(i,1:nodays,43) = FLUXES(1:nodays,38)         ! transpiration (kgH2O.m-2.day-1)
     out_var1(i,1:nodays,44) = FLUXES(1:nodays,39)         ! soil evaporation (kgH2O.m-2.day-1)
     out_var1(i,1:nodays,45) = FLUXES(1:nodays,40)         ! wet canopy evaporation (kgH2O.m-2.day-1)
     out_var1(i,1:nodays,46) = FLUXES(1:nodays,41)         ! runoff (kgH2O.m-2.day-1)
     out_var1(i,1:nodays,47) = FLUXES(1:nodays,42)         ! underflow (kgH2O.m-2.day-1)
     out_var1(i,1:nodays,48) = POOLS(1:nodays,8)           ! surface water (kgH2O.m-2.30cmdepth)
     out_var1(i,1:nodays,49) = wSWP_time(1:nodays)         ! Weighted Soil Water Potential (MPa)
     out_var1(i,1:nodays,50) = snow_storage_time(1:nodays) ! Snow storage (kgH2O/m2)
     ! Canopy (phenology) properties
     out_var1(i,1:nodays,51) = lai                         ! LAI (m2/m2)
     ! Photosynthesis / C~water coupling related
     out_var1(i,1:nodays,52) = gs_demand_supply_ratio      ! ratio of evaporative demand over supply
     out_var1(i,1:nodays,53) = gs_total_canopy             ! Canopy scale stomatal conductance during day light (mmolH2O/m2ground/s)
     out_var1(i,1:nodays,54) = canopy_par_MJday_time       ! Canopy absorbed PAR (MJ/m2ground/day)
     out_var1(i,1:nodays,55) = gb_total_canopy             ! Canopy scale aerodynamic conductance (mmolH2O/m2ground/s)
     out_var1(i,1:nodays,56) = cica_time                   ! ratio of leaf internal to external CO2
     ! misc
     out_var1(i,1:nodays,57) = root_depth_time             ! rooting depth (m)
     out_var1(i,1:nodays,58) = DS_time                     ! Development stage (0-2)

     !
     ! Calculate long-term mean of out_var1
     !
     
     ! Loop across each variable
     do v = 1, output_dim
        ! Calculate mean value
        out_var4(i,v) = sum(out_var1(i,1:nodays,v)) / dble(nodays)
     end do

     !
     ! Calculate the mean annual of out_var1
     !

     ! Calculate mean annual
     s = 1 ; e = steps_per_year
     do a = 1, nos_years
        do v = 1, output_dim
           out_var5(i,a,v) = sum(out_var1(i,s:e,v)) / dble(steps_per_year)
        end do
        ! Iterate counters
        s = s + steps_per_year ; e = s + steps_per_year - 1
     end do
     
     !!!
     ! Estimate residence time information
     !!!

     ! Labile
     ! Estimate MRT (years)
     pool_hak = 1 ; tmp = 0d0
     where (POOLS(1:nodays,1) > 0d0) ! protection against NaN from division by zero
            pool_hak = 0 
            tmp = ((FLUXES(1:nodays,8) + FLUXES(1:nodays,31) + FLUXES(1:nodays,36)) &
                  / POOLS(1:nodays,1))
     end where
     out_var2(i,1) = sum(tmp) / dble(nodays-sum(pool_hak))
     ! Foliage
     ! Estimate MRT (years)
     pool_hak = 1 ; tmp = 0d0
     where (POOLS(1:nodays,2) > 0d0) ! protection against NaN from division by zero
            pool_hak = 0 
            tmp = ((FLUXES(1:nodays,10) + FLUXES(1:nodays,28) + FLUXES(1:nodays,32)) &
                  / POOLS(1:nodays,2))
     end where
     out_var2(i,2) = sum(tmp) / dble(nodays-sum(pool_hak))
     ! Fine roots
     ! Estimate MRT (years)
     pool_hak = 1 ; tmp = 0d0
     where (POOLS(1:nodays,3) > 0d0) ! protection against NaN from division by zero
            pool_hak = 0 
            tmp = ((FLUXES(1:nodays,12) + FLUXES(1:nodays,36)) &
                  / POOLS(1:nodays,3))
     end where
     out_var2(i,3) = sum(tmp) / dble(nodays-sum(pool_hak))
     ! Wood
     ! Estimate MRT (years)
     pool_hak = 1 ; tmp = 0d0
     where (POOLS(1:nodays,4) > 0d0) ! protection against NaN from division by zero
            pool_hak = 0 
            tmp = ((FLUXES(1:nodays,11) + FLUXES(1:nodays,29) + FLUXES(1:nodays,33)) &
                  / POOLS(1:nodays,4))
     end where
     out_var2(i,4) = sum(tmp) / dble(nodays-sum(pool_hak))
     ! Foliage + fine root litter
     ! Estimate MRT (years)
     pool_hak = 1 ; tmp = 0d0
     where (POOLS(1:nodays,5) > 0d0) ! protection against NaN from division by zero
            pool_hak = 0 
            tmp = ((FLUXES(1:nodays,13) + FLUXES(1:nodays,15)) &
                  / POOLS(1:nodays,5))
     end where
     out_var2(i,5) = sum(tmp) / dble(nodays-sum(pool_hak))
     ! Soil
     ! Estimate MRT (years)
     pool_hak = 1 ; tmp = 0d0
     where (POOLS(1:nodays,6) > 0d0) ! protection against NaN from division by zero
            pool_hak = 0 
            tmp = ((FLUXES(1:nodays,14)) &
                  / POOLS(1:nodays,6))
     end where
     out_var2(i,6) = sum(tmp) / dble(nodays-sum(pool_hak))
     ! Autotrophic
     ! Estimate MRT (years)
     pool_hak = 1 ; tmp = 0d0
     where (POOLS(1:nodays,7) > 0d0) ! protection against NaN from division by zero
            pool_hak = 0 
            tmp = ((FLUXES(1:nodays,23) + FLUXES(1:nodays,35)) &
                  / POOLS(1:nodays,7))
     end where
     out_var2(i,7) = sum(tmp) / dble(nodays-sum(pool_hak))
     ! Dead still standing foliage
     ! Estimate MRT (years)
     pool_hak = 1 ; tmp = 0d0
     where (POOLS(1:nodays,10) > 0d0) ! protection against NaN from division by zero
            pool_hak = 0
            tmp = ((FLUXES(1:nodays,30) + FLUXES(1:nodays,34)) &
                  / POOLS(1:nodays,10))
     end where
     out_var2(i,8) = sum(tmp) / dble(nodays-sum(pool_hak))

     !
     ! Estimate pool inputs needed for steady state calculation
     !

     ! Once the canopy has closes the inputs to the live biomass are stable
     ! and can thus be estimated from the simulated inputs
     out_var3(i,1) = sum(FLUXES(:,5)) ! Labile
     out_var3(i,2) = sum(FLUXES(:,4)) ! Foliage
     out_var3(i,3) = sum(FLUXES(:,6)) ! Fine root
     out_var3(i,4) = sum(FLUXES(:,7)) ! Wood
     out_var3(i,5) = sum(FLUXES(:,12) + FLUXES(:,32) + FLUXES(:,33) + &
                         FLUXES(:,34) + FLUXES(:,35) + FLUXES(:,36) + FLUXES(:,37)) ! litter (foliage + roots)
     ! In the normal DALEC models where wood can be a major input to the soil pool,
     ! at this point we would account for disturbance inputs (including wood)
     ! but NOT natural wood. However, as crops are inherently annual in their time
     ! scales we combine all inputs here.
     out_var3(i,6) = sum(FLUXES(:,15)) ! som
     out_var3(i,7) = sum(FLUXES(:,16)) ! autotrophic
     out_var3(i,8) = sum(FLUXES(:,10)*0.5d0) ! dead still standing foliage

  end do ! nos_iter loop

  ! MTT - Convert daily fractional loss to years
  out_var2 = (out_var2*365.25d0)**(-1d0) ! iter,(lab,fol,root,wood,lit,som)

  ! Steady state gC/m2 estimation
  ! Determine the mean annual input (gC/m2/yr) based on current inputs for all pool,
  ! litter and soil pools updated below...
  out_var3 = (out_var3 / dble(nodays)) * 365.25d0 ! convert to annual mean input
  ! Then estimate the labile, foliar, fine root, wood, litter and soil,
  ! autotrophic and dead standing foliage steady states.
  out_var3(1:nos_iter,1:8) = out_var3(1:nos_iter,1:8) * out_var2(1:nos_iter,1:8) ! multiply by residence time in years

  ! return back to the subroutine then
  return

end subroutine rdalec15
!
!------------------------------------------------------------------
!
  subroutine crop_development_parameters(stock_seed_labile,DS_shoot,DS_root,fol_frac &
                                        ,stem_frac,root_frac,DS_LRLV,LRLV,DS_LRRT,LRRT &
                                        ,exepath)

    ! subroutine reads in the fixed crop development files which are linked the
    ! the development state of the crops. The development model varies between
    ! which species. e.g. winter wheat and barley, spring wheat and barley
    ! NOTE: duplicate function in the *_PARS.f90

    implicit none

    ! declare inputs
    ! crop specific variables
    !integer,intent(in) :: pathlength
    character(350),intent(in) :: exepath
    double precision :: stock_seed_labile
    double precision, allocatable, dimension(:) :: DS_shoot, & !
                                                    DS_root, & !
                                                   fol_frac, & !
                                                  stem_frac, & !
                                                  root_frac, & !
                                                    DS_LRLV, & !
                                                       LRLV, & !
                                                    DS_LRRT, & !
                                                       LRRT

    ! local variables..
    integer :: columns, i, rows, input_crops_unit, ios
    character(225) :: variables,filename

    ! file info needed
    input_crops_unit = 20 ; ios = 0

    ! crop development file passed in from the R code (this is different from
    ! *_PARS.f90 where this subroutine is hardcoded)
    open(unit = input_crops_unit, file=trim(exepath),iostat=ios, status='old', action='read')

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
