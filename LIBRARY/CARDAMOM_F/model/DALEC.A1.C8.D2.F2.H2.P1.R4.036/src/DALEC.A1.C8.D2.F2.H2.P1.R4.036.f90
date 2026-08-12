!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CARbon DAta MOdel fraMework (CARDAMOM) and DALEC terrestrial ecosystem model suite
! CARDAMOM is a Bayesian model-data fusion software framework. CARDAMOM is used to 
! assimilate observations and ecological theory to retrieve parameters for the 
! DALEC suite of intermediate complexity terrestrial ecosystem models. DALEC can be
! used as a fully integrated component of CARDAMOM or independently. 
! Copyright (C) 2024  University of Edinburgh,
!                     Mathew Williams (mat.williams@ed.ac.uk), 
!                     T. Luke Smallman (t.l.smallman@ed.ac.uk)
! UoE = University of Edinburgh

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

!!!!!!!!!!!! File specific description !!!!!!!!!!
! This file contains the source code of DALEC.A1.C8.D2.F2.H2.P1.R4.036
!
! This code contains a variant of the Data Assimilation Linked ECosystem (DALEC) model.
! This version of DALEC is derived from the following primary references:
! Xenakis & Williams (2014) https://doi.org/10.5194/gmd-7-1519-2014
! Bloom & Williams (2015), https://doi.org/10.5194/bg-12-1299-2015.
! Smallman & Williams (2019) https://doi.org/10.5194/gmd-12-2227-2019.
! Thomas et al., (2019), https://doi.org/10.1029/2019MS001679
! This code is based on that created by A. A. Bloom (UoE, now at JPL, USA).
! Subsequent modifications by:
! T. L. Smallman (University of Edinburgh, t.l.smallman@ed.ac.uk)
! J. F. Exbrayat (University of Edinburgh)
! See function / subroutine specific comments for exceptions and contributors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module CARBON_MODEL_MOD

  implicit none

  ! make all private
  private

  ! explicit publics
  public :: CARBON_MODEL  &
           ,initialize_mv


  contains
  !
  !--------------------------------------------------------------------
  !
  subroutine initialize_mv(mV, nodays, nomet, nopars, deltat, soil_frac_sand, soil_frac_clay, met, lat)
    use carbon_model_memory

    !! For a single chain's model_working_varibles type object mV, allocate arrays
    !! and calculate initial values.
    implicit none

    type(model_working_variables), intent(out):: mV
    integer, intent(in):: nodays, nomet, nopars
    double precision, intent(in) :: deltat(nodays)     ! time step in decimal days
    double precision, intent(in), dimension(:) :: soil_frac_sand, soil_frac_clay
    double precision, intent(in) :: met(nomet, nodays)  ! met drivers
    double precision, intent(in) :: lat

    integer:: n

    mV%soil_frac_sand = soil_frac_sand
    mV%soil_frac_clay = soil_frac_clay

    ! allocate variables dimension which are fixed per site only the once
    allocate(mV%deltat_1(nodays),  &
             mV%daylength_hours(nodays), mV%daylength_seconds(nodays), mV%daylength_seconds_1(nodays), &
             mV%rainfall_time(nodays), mV%airt_zero_fraction_time(nodays))

    !
    ! Timing variables which are needed first
    !

    mV%deltat_1 = deltat**(-1d0) 

    !
    ! Iteration independent variables using functions and thus need to be in a loop
    !

    ! first those linked to the time period of the analysis
    do n = 1, nodays
       ! check positive values only for rainfall input
       mV%rainfall_time(n) = max(0d0,met(7,n))
       ! calculate daylength in hours and seconds
       call calculate_daylength((met(6,n)-(deltat(n)*0.5d0)),lat, mV)
       mV%daylength_hours(n) = mV%dayl_hours ; mV%daylength_seconds(n) = mV%dayl_seconds
    end do

    ! calculate inverse for each time step in seconds
    mV%daylength_seconds_1 = mV%daylength_seconds ** (-1d0)
    ! fraction of temperture period above freezing
    mV%airt_zero_fraction_time = 0d0
    where (met(2,:) > 0d0) mV%airt_zero_fraction_time = 1d0
    where (met(3,:) > 0d0 .and. met(2,:) < 0d0) mV%airt_zero_fraction_time = (met(3,:)-0d0) / (met(3,:)-met(2,:))

    ! number of time steps per year
    mV%steps_per_year = nint(dble(nodays)/(sum(deltat)*0.002737851d0))
    ! mean days per step
    mV%mean_days_per_step = sum(deltat) / dble(nodays)

    !
    ! Initialise the water model
    !

    ! zero variables not done elsewhere
    mV%total_water_flux = 0d0 ; mV%water_flux_mmolH2Om2s = 0d0
    ! initialise some time invarient parameters
    call saxton_parameters(mV%soil_frac_clay,mV%soil_frac_sand, mV)
    call initialise_soils(mV%soil_frac_clay,mV%soil_frac_sand, mV)

    !call update_soil_initial_conditions(pars(24), mV) 
    ! can't be done here because of dependence on PARS
    ! that's ok because it's done in every carbon_model() call at else block
    ! does not affect field_capacity or porosity, or anything used before
    ! the call in CARBON_MODEL else block

    ! save the initial conditions for later
    mV%field_capacity_initial = mV%field_capacity
    mV%porosity_initial = mV%porosity

  end subroutine initialize_mv
  !
  !--------------------------------------------------------------------
  !
  subroutine destroy_mv(mV)
    use carbon_model_memory
    !! deallocate arrays in mV
    type(model_working_variables):: mV
    integer:: n
        ! allocate variables dimension which are fixed per site only the once
        deallocate(mV%deltat_1, &
                   mV%daylength_hours, mV%daylength_seconds, mV%daylength_seconds_1, &
                   mV%rainfall_time, mV%airt_zero_fraction_time)
  end subroutine
  !
  !--------------------------------------------------------------------
  !
  subroutine CARBON_MODEL(start,finish,met,pars,deltat,nodays,lat,FLUXES,POOLS,DIAGS &
                         ,nopars,nomet,nopools,nofluxes,nodiags, mV)
    use carbon_model_memory

    ! The Data Assimilation Linked Ecosystem Carbon (DALEC.A1.C8.D2.F2.H2.P1.R4.036) model.
    ! GPP is estimated using the Aggregated Canopy Model (ACM; Williams et al. 1997) driven by
    ! canopy efficiency, radiation, temperature and CO2. Maintenance autotrophic respiration is
    ! a fixed fraction of GPP; remaining NPP is allocated to foliage, labile, fine roots and wood
    ! using fixed fractions with growth respiration deducted at Rg_fraction = 0.21875.
    ! Labile C is released to foliage via CDEA sinusoidal phenology (bud-burst and leaf-fall
    ! controlled by DOY and duration parameters). Dead organic matter cycles through six pools:
    ! foliar litter, fine root litter, wood litter, fast SOM, slow SOM and a microbial pool.
    ! Decomposition of litter to fast or slow SOM is partitioned by lignin fractions (p41-p43).
    ! Microbial dynamics follow Xenakis & Williams (2014). Soil hydrology uses a 3-layer BUCKET
    ! model with snow accumulation and melt. Fire combustion and deforestation disturbance are
    ! implemented with pool-specific combustion completeness and five management scenarios (R4).

    implicit none

    ! declare input variables
    type(model_working_variables) :: mV
    integer, intent(in) :: start    &
                          ,finish   &
                          ,nopars   & ! number of parameters in vector
                          ,nomet    & ! number of meteorological fields
                          ,nofluxes & ! number of model fluxes
                          ,nopools  & ! number of model pools
                          ,nodays   & ! number of days in simulation
                          ,nodiags    ! number of model diagnositic variables

    double precision, intent(in) :: met(nomet,nodays) & ! met drivers
                         ,deltat(nodays)    & ! time step in decimal days
                         ,pars(nopars)      & ! number of parameters
                         ,lat                 ! site latitude (degrees)

    double precision, dimension((nodays+1),nopools), intent(inout) :: POOLS ! vector of ecosystem pools
    double precision, dimension(nodays,nofluxes), intent(inout) :: FLUXES ! vector of ecosystem fluxes
    double precision, dimension(nodays,nodiags), intent(inout) :: DIAGS ! vector of ecosystem diagnostics

    ! declare local variables
    double precision ::  tmp,infi &
              ,microbial_activity & 
                 ,microbial_death &
                   ,transpiration & ! kgH2O/m2/day
                 ,soilevaporation & ! kgH2O/m2/day
                  ,wetcanopy_evap & ! kgH2O/m2/day
                 ,snowsublimation & ! kgH2O/m2/day
       ,wf,wl,ff,fl,osf,osl,sf,ml & ! phenological controls
       ,sla                       & ! specific leaf area: inverse of LMA (m2/gC); precomputed from 1/pars(17)
       ,ff_wf,fl_wl               & ! precomputed phenology amplitude coefficients
       ,sf_over_wf,sf_over_wl       ! precomputed phenology scaling ratios

    ! JFE added 4 May 2018 - combustion efficiencies and fire resilience
    double precision :: burnt_area
    double precision, dimension(9) :: cf,rfac
    ! local deforestation related variables
    double precision, dimension(5) :: post_harvest_burn      & ! how much burning to occur after
                                     ,foliage_frac_res       &
                                     ,roots_frac_res         &
                                     ,rootcr_frac_res        &
                                     ,stem_frac_res          &
                                     ,roots_frac_removal     &
                                     ,rootcr_frac_removal    &
                                     ,Crootcr_part           &
                                     ,soil_loss_frac
    double precision :: labile_loss,foliar_loss      &
                       ,roots_loss,wood_loss         &
                       ,rootcr_loss,stem_loss        &
                       ,labile_residue,foliar_residue&
                       ,roots_residue,wood_residue   &
                       ,C_total,labile_frac_res      &
                       ,labile_frac_removal          &
                       ,Cstem,Crootcr,stem_residue   &
                       ,coarse_root_residue          &
                       ,soil_loss_with_roots

    integer :: harvest_management,n

    ! met drivers are:
    ! 1st run day
    ! 2nd min daily temp (oC)
    ! 3rd max daily temp (oC)
    ! 4th Radiation (MJ.m-2.day-1)
    ! 5th CO2 (ppm)
    ! 6th DOY
    ! 7th precipitation (kgH2O.m-2.s-1)
    ! 8th deforestation fraction
    ! 9th burnt area fraction
    ! 10th 21 day average min temperature (oC)
    ! 11th 21 day average photoperiod (seconds)
    ! 12th 21 day average VPD (Pa)
    ! 13th Forest management practice to accompany any clearing
    ! 14th avg daily temperature (oC)
    ! 15th avg daily wind speed (m.s-1)
    ! 16th vapour pressure deficit (Pa)

    ! POOLS are:
    ! 1  = labile (gC/m2) (initial value: p18)
    ! 2  = foliar (gC/m2) (initial value: p19)
    ! 3  = fine root (gC/m2) (initial value: p20)
    ! 4  = wood (gC/m2) (initial value: p21)
    ! 5  = foliar litter (gC/m2) (initial value: p22)
    ! 6  = fine root litter (gC/m2) (initial value: p36)
    ! 7  = wood litter (gC/m2) (initial value: p34)
    ! 8  = fast SOM (gC/m2) (initial value: p37)
    ! 9  = slow SOM (gC/m2) (initial value: p23)
    ! 10 = microbial (gC/m2) (initial value: p38)
    ! 11 = surface soil water 0-10 cm (mm) (initial value: p24)

    ! FLUXES are:
    ! 1  = GPP (gC/m2/day)
    ! 2  = temperature rate modifier (unitless)
    ! 3  = total autotrophic respiration — maintenance Ra (fraction of GPP) +
    !       growth respiration from all allocations (gC/m2/day)
    ! 4  = direct allocation to foliage from GPP, pre-growth-respiration (gC/m2/day)
    ! 5  = direct allocation to labile from GPP (gC/m2/day)
    ! 6  = direct allocation to fine roots from GPP, pre-growth-respiration (gC/m2/day)
    ! 7  = direct allocation to wood from GPP, pre-growth-respiration (gC/m2/day)
    ! 8  = labile release to foliage, after growth respiration deduction (gC/m2/day)
    ! 9  = leaf fall factor (unitless, CDEA phenology)
    ! 10 = total leaf litter production (gC/m2/day)
    ! 11 = wood turnover to wood litter (gC/m2/day)
    ! 12 = fine root litter production (gC/m2/day)
    ! 13 = heterotrophic respiration from foliar litter (gC/m2/day)
    ! 14 = heterotrophic respiration from slow SOM (gC/m2/day)
    ! 15 = foliar litter decomposition to SOM (gC/m2/day)
    ! 16 = labile release factor (unitless, CDEA phenology)
    ! 17 = total ecosystem fire emission — sum(F18:F23, F32, F63, F64) (gC/m2/day)
    ! 18 = fire emission from labile (gC/m2/day)
    ! 19 = fire emission from foliar (gC/m2/day)
    ! 20 = fire emission from fine roots (gC/m2/day)
    ! 21 = fire emission from wood (gC/m2/day)
    ! 22 = fire emission from foliar litter (gC/m2/day)
    ! 23 = fire emission from slow SOM (gC/m2/day)
    ! 24 = fire mortality: labile → fast SOM (gC/m2/day)
    ! 25 = fire mortality: foliar → foliar litter (gC/m2/day)
    ! 26 = fire mortality: fine roots → fine root litter (gC/m2/day)
    ! 27 = fire mortality: wood → wood litter (gC/m2/day)
    ! 28 = fire mortality: foliar litter → SOM (gC/m2/day)
    ! 29 = evapotranspiration (kgH2O/m2/day)
    ! 30 = heterotrophic respiration from wood litter (gC/m2/day)
    ! 31 = wood litter decomposition to SOM (gC/m2/day)
    ! 32 = fire emission from wood litter (gC/m2/day)
    ! 33 = fire mortality: wood litter → SOM (gC/m2/day)
    ! 34 = harvest extraction from labile (gC/m2/day)
    ! 35 = harvest extraction from foliage (gC/m2/day)
    ! 36 = harvest extraction from fine roots (gC/m2/day)
    ! 37 = harvest extraction from wood (gC/m2/day)
    ! 38 = harvest extraction from foliar litter (gC/m2/day)
    ! 39 = harvest extraction from wood litter (gC/m2/day)
    ! 40 = harvest: mechanical soil C loss with coarse root extraction (gC/m2/day)
    ! 41 = harvest residue: labile → fast SOM (gC/m2/day)
    ! 42 = harvest residue: foliage → foliar litter (gC/m2/day)
    ! 43 = harvest residue: fine roots → fine root litter (gC/m2/day)
    ! 44 = harvest residue: wood → wood litter (gC/m2/day)
    ! 45 = total harvest C extraction — sum(F34:F44) (gC/m2/day)
    ! 46 = transpiration (kgH2O/m2/day)
    ! 47 = soil evaporation (kgH2O/m2/day)
    ! 48 = wet canopy evaporation (kgH2O/m2/day)
    ! 49 = surface runoff (kgH2O/m2/day)
    ! 50 = drainage from bottom of soil column (kgH2O/m2/day)
    ! 51 = gravitational drainage from surface soil layer to 2nd layer (kgH2O/m2/day)
    ! 52 = infiltration into top soil layer by rain (kgH2O/m2/day)
    ! 53 = transpiration fraction extracted from 1st rooting layer (fraction)
    ! 54 = transpiration fraction extracted from 2nd rooting layer (fraction)
    ! 55 = growth respiration from labile release to foliage (gC/m2/day)
    ! 56 = heterotrophic respiration from fine root litter (gC/m2/day)
    ! 57 = fine root litter decomposition to SOM (gC/m2/day)
    ! 58 = heterotrophic respiration from fast SOM (gC/m2/day)
    ! 59 = heterotrophic respiration from microbial pool (gC/m2/day)
    ! 60 = microbial death flux → slow SOM (gC/m2/day)
    ! 61 = microbial-mediated transfer: slow SOM → fast SOM (gC/m2/day)
    ! 62 = fast SOM uptake into microbial pool (gC/m2/day)
    ! 63 = fire emission from fine root litter (gC/m2/day)
    ! 64 = fire emission from fast SOM (gC/m2/day)
    ! 65 = fire mortality: fine root litter → SOM (gC/m2/day)
    ! 66 = infiltration into middle soil layer by rain (kgH2O/m2/day)
    ! 67 = infiltration into bottom soil layer by rain (kgH2O/m2/day)

    ! PARAMETERS
    ! p(1)  = decomposition efficiency of litter to SOM (fraction)
    ! p(2)  = fraction of GPP as maintenance autotrophic respiration (fraction)
    ! p(3)  = fraction of (GPP-Ra) allocated to foliage (fraction)
    ! p(4)  = fraction of (GPP-Ra-fol-lab) allocated to fine roots (fraction)
    ! p(5)  = leaf lifespan (yr)
    ! p(6)  = wood turnover rate (fraction/day)
    ! p(7)  = fine root turnover rate (fraction/day)
    ! p(8)  = foliar litter turnover rate, temperature adjusted (fraction/day)
    ! p(9)  = fine root litter turnover rate, temperature adjusted (fraction/day)
    ! p(10) = temperature sensitivity of heterotrophic respiration (oC-1)
    ! p(11) = canopy efficiency — ACM parameter (umolC/m2/s)
    ! p(12) = day of year for peak labile release / bud-burst — CDEA phenology (day)
    ! p(13) = fraction of (GPP-Ra-fol) allocated to labile (fraction)
    ! p(14) = labile release duration period (days)
    ! p(15) = day of year for peak leaf fall — CDEA phenology (day)
    ! p(16) = leaf fall duration period (days)
    ! p(17) = LMA (gC/m2)
    ! p(18) = initial labile C (gC/m2)
    ! p(19) = initial foliar C (gC/m2)
    ! p(20) = initial fine root C (gC/m2)
    ! p(21) = initial wood C (gC/m2)
    ! p(22) = initial foliar litter C (gC/m2)
    ! p(23) = initial slow SOM C (gC/m2)
    ! p(24) = initial soil water fraction (fraction)
    ! p(25) = fraction of Cwood as coarse root (fraction)
    ! p(26) = coarse root biomass for 50% of maximum rooting depth (gBiomass/m2)
    ! p(27) = maximum rooting depth (m)
    ! p(28) = fire resilience factor for live pools (fraction)
    ! p(29) = combustion completeness for foliage (fraction)
    ! p(30) = combustion completeness for labile, fine roots and wood (fraction)
    ! p(31) = combustion completeness for SOM (fraction)
    ! p(32) = combustion completeness for foliar and fine root litter (fraction)
    ! p(33) = combustion completeness for wood litter (fraction)
    ! p(34) = initial wood litter C (gC/m2)
    ! p(35) = wood litter turnover rate, temperature adjusted (fraction/day)
    ! p(36) = initial fine root litter C (gC/m2)
    ! p(37) = initial fast SOM C (gC/m2)
    ! p(38) = initial microbial C (gC/m2)
    ! p(39) = initial microbial activity (unitless)
    ! p(40) = inhibition constant for C-dependent microbial activity (gC/m2)
    ! p(41) = foliar lignin fraction controlling litter decomposition partitioning to slow SOM (fraction)
    ! p(42) = fine root lignin fraction controlling litter decomposition partitioning to slow SOM (fraction)
    ! p(43) = wood lignin fraction controlling litter decomposition partitioning to slow SOM (fraction)
    ! p(44) = microbial substrate uptake efficiency from fast SOM (fraction)
    ! p(45) = maximum microbial death rate (fraction/day)
    ! p(46) = inhibition constant for microbial death (unitless)
    ! p(47) = microbial maintenance respiration coefficient (fraction/day)
    ! p(48) = microbial decomposition efficiency — fraction of slow SOM turnover to fast SOM (fraction)
    ! p(49) = turnover constant for slow SOM via microbial activity (day-1)
    ! p(50) = 2nd order rate constant for microbial uptake from fast SOM (day-1)

    ! Set some initial states
    infi = 0d0 ; FLUXES = 0d0 ; POOLS = 0d0 ; DIAGS = 0d0
    ! Reset hydrology variables
    mV%intercepted_rainfall = 0d0 ; mV%canopy_storage = 0d0 ; mV%snow_storage = 0d0
    transpiration = 0d0 ; soilevaporation = 0d0 ; wetcanopy_evap = 0d0 ; snowsublimation = 0d0
    ! Reset radiation variabes
    mV%canopy_swrad_MJday = 0d0 ; mV%canopy_par_MJday = 0d0 ; mV%soil_swrad_MJday = 0d0
    mV%canopy_lwrad_Wm2 = 0d0 ; mV%soil_lwrad_Wm2 = 0d0 ; mV%sky_lwrad_Wm2 = 0d0
    ! Reset conductance variables
    mV%soil_conductance = 0d0

    ! load ACM-GPP-ET parameters
    mV%ceff = pars(11) ! Canopy efficiency (umolC/m2/s)
                       ! This is in the full model the product of Nitrogen use efficiency (umolC/gN/m2leaf)
                       ! and average foliar nitrogen gN/m2leaf
    ! Rooting parameters
    mV%root_k = pars(26) ; mV%max_depth = pars(27)
    ! Specific leaf area: precomputed reciprocal of LMA (pars(17)) to replace per-step division
    sla = 1d0 / pars(17)

    ! assigning initial conditions
    POOLS(1,1) = pars(18) ! labile
    POOLS(1,2) = pars(19) ! foliar
    POOLS(1,3) = pars(20) ! roots
    POOLS(1,4) = pars(21) ! wood
    POOLS(1,5) = pars(22) ! foliage litter
    POOLS(1,6) = pars(36) ! root litter
    POOLS(1,7) = pars(34) ! wood litter
    POOLS(1,8) = pars(37) ! fast som
    POOLS(1,9) = pars(23) ! slow som
    POOLS(1,10) = pars(38)! microbial  
    !POOLS(1,11) = assigned later ! soil water (0-10cm)

    ! Initialise the microbial activity 
    microbial_activity = pars(39)

    ! Some time consuming variables we only want to set once
    if (.not.allocated(mV%deltat_1)) then !never used - initialize_model() is always called before carbon_model()
        write(*,*) "Error - arrays not allocated - probably carbon_model() was called without initialize_model()"
        STOP 1
    else ! deltat_1 allocated?

        !
        ! Load initial soil water conditions from memory
        !

        mV%total_water_flux = 0d0; mV%water_flux_mmolH2Om2s = 0d0
        mV%field_capacity = mV%field_capacity_initial
        mV%porosity = mV%porosity_initial

        ! input initial soil water fraction then
        ! update SWP and soil conductivity accordingly
        call update_soil_initial_conditions(pars(24), mV)

    endif  ! deltat_1 allocated

    ! Defining phenological variables
    ! release period coefficient, based on duration of labile turnover or leaf
    ! fall durations
    wf = pars(16)*sqrt(2d0) * 0.5d0
    wl = pars(14)*sqrt(2d0) * 0.5d0
    ! magnitude coefficient
    ff = (log(pars(5))-log(pars(5)-1d0)) * 0.5d0
    fl = 3.45437738965761021d0!(log(1.001d0)-log(0.001d0)) * 0.5d0
    ! set minium labile life span to one year
    ml = 1.001d0
    ! offset for labile and leaf turnovers
    osf = ospolynomial(pars(5),wf)
    osl = ospolynomial(ml,wl)
    ! scaling to biyearly sine curve
    sf = 116.262685928629551d0 !365.25d0/pi
    ! precompute constant phenology amplitude coefficients and scaling ratios
    ff_wf     = (2d0/sqrt(pi)) * (ff/wf) ! leaf fall amplitude coefficient
    fl_wl     = (2d0/sqrt(pi)) * (fl/wl) ! labile release amplitude coefficient
    sf_over_wf = sf/wf                   ! leaf fall scaling ratio
    sf_over_wl = sf/wl                   ! labile release scaling ratio

    ! now load the hardcoded forest management parameters into their scenario locations

    ! Deforestation process functions in a sequential way.
    ! Thus, the pool_loss is first determined as a function of met(8,n) and
    ! for fine and coarse roots whether this felling is associated with a mechanical
    ! removal from the ground. As the canopy and stem is removed (along with a proportion of labile)
    ! fine and coarse roots may subsequently undergo mortality from which they do not recover
    ! but allows for management activities such as grazing, mowing and coppice.
    ! The pool_loss is then partitioned between the material which is left within the system
    ! as a residue and thus directly placed within one of the dead organic matter pools.

    !! Parameter values for deforestation variables
    !! Scenario 1
    ! Define 'removal' for coarse and fine roots, i.e. fraction of imposed
    ! removal which is imposed directly on these pools. These fractions vary
    ! the assumption that the fine and coarse roots are mechanically removed.
    ! 1 = all removed, 0 = all remains.
    roots_frac_removal(1)  = 0d0
    rootcr_frac_removal(1) = 0d0
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(1) = 1d0
    roots_frac_res(1)   = 1d0
    rootcr_frac_res(1)  = 1d0
    stem_frac_res(1)    = 0.20d0 !
    ! wood partitioning (fraction)
    Crootcr_part(1) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
    ! Csom loss due to phyical removal with roots
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(1) = 0.02d0 ! actually between 1-3 %
    ! was the forest burned after deforestation (0-1)
    ! NOTE: that we refer here to the fraction of the cleared land to be burned
    post_harvest_burn(1) = 1d0

    !! Scenario 2
    ! Define 'removal' for coarse and fine roots, i.e. fraction of imposed
    ! removal which is imposed directly on these pools. These fractions vary
    ! the assumption that the fine and coarse roots are mechanically removed.
    ! 1 = all removed, 0 = all remains.
    roots_frac_removal(2)  = 0d0
    rootcr_frac_removal(2) = 0d0
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(2) = 1d0
    roots_frac_res(2)   = 1d0
    rootcr_frac_res(2) = 1d0
    stem_frac_res(2)   = 0.20d0 !
    ! wood partitioning (fraction)
    Crootcr_part(2) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
    ! Csom loss due to phyical removal with roots
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(2) = 0.02d0 ! actually between 1-3 %
    ! was the forest burned after deforestation (0-1)
    ! NOTE: that we refer here to the fraction of the cleared land to be burned
    post_harvest_burn(2) = 0d0

    !! Scenario 3
    ! Define 'removal' for coarse and fine roots, i.e. fraction of imposed
    ! removal which is imposed directly on these pools. These fractions vary
    ! the assumption that the fine and coarse roots are mechanically removed.
    ! 1 = all removed, 0 = all remains.
    roots_frac_removal(3)  = 0d0
    rootcr_frac_removal(3) = 0d0
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(3) = 0.5d0
    roots_frac_res(3)   = 1d0
    rootcr_frac_res(3) = 1d0
    stem_frac_res(3)   = 0d0 !
    ! wood partitioning (fraction)
    Crootcr_part(3) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
    ! Csom loss due to phyical removal with roots
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(3) = 0.02d0 ! actually between 1-3 %
    ! was the forest burned after deforestation (0-1)
    ! NOTE: that we refer here to the fraction of the cleared land to be burned
    post_harvest_burn(3) = 0d0

    !! Scenario 4
    ! Define 'removal' for coarse and fine roots, i.e. fraction of imposed
    ! removal which is imposed directly on these pools. These fractions vary
    ! the assumption that the fine and coarse roots are mechanically removed.
    ! 1 = all removed, 0 = all remains.
    roots_frac_removal(4)  = 1d0
    rootcr_frac_removal(4) = 1d0
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(4) = 0.5d0
    roots_frac_res(4)   = 1d0
    rootcr_frac_res(4) = 0d0
    stem_frac_res(4)   = 0d0
    ! wood partitioning (fraction)
    Crootcr_part(4) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
    ! Csom loss due to phyical removal with roots
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(4) = 0.02d0 ! actually between 1-3 %
    ! was the forest burned after deforestation (0-1)
    ! NOTE: that we refer here to the fraction of the cleared land to be burned
    post_harvest_burn(4) = 0d0

    !## Scenario 5 (grassland grazing / cutting)
    ! Define 'removal' for coarse and fine roots, i.e. fraction of imposed
    ! removal which is imposed directly on these pools. These fractions vary
    ! the assumption that the fine and coarse roots are mechanically removed.
    ! 1 = all removed, 0 = all remains.
    roots_frac_removal(5)  = 0d0
    rootcr_frac_removal(5) = 0d0
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(5) = 0.1d0
    roots_frac_res(5)   = 0d0
    rootcr_frac_res(5)  = 0d0
    stem_frac_res(5)    = 0.12d0
    ! wood partitioning (fraction)
    Crootcr_part(5) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
    ! Csom loss due to phyical removal with roots
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(5) = 0d0 ! actually between 1-3 %
    ! was the forest burned after deforestation (0-1)
    ! NOTE: that we refer here to the fraction of the cleared land to be burned
    post_harvest_burn(5) = 0d0

    ! Override all paritioning parameters with those coming from
    ! CARDAMOM
    Crootcr_part = pars(25)

!    ! Declare combustion efficiency (labile, foliar, roots, wood, litter, soil, woodlitter)
!    cf(1) = 0.1d0 ; cf(2) = 0.9d0
!    cf(3) = 0.1d0 ; cf(4) = 0.1d0
!    cf(5) = 0.7d0 ; cf(6) = 0.01d0
!    ! Resilience factor for non-combusted tissue
!    rfac = 0.5d0 ; rfac(5) = 0.1d0 ; rfac(6) = 0d0

    ! JFE added 4 May 2018 - define fire constants
    ! Update fire parameters derived from
    ! Yin et al., (2020), doi: 10.1038/s414647-020-15852-2
    ! Subsequently expanded by T. L. Smallman & Mat Williams (UoE, 03/09/2021)
    ! to provide specific CC for litter and wood litter.
    ! NOTE: changes also result in the addition of further EDCs

    ! Assign proposed resilience factor
    rfac(1:4) = pars(28) ! live pools
    rfac(5) = 0.1d0 ! foliar litter 
    rfac(6) = 0.1d0 ! fine root litter
    rfac(7) = 0.1d0 ! wood litter 
    rfac(8) = 0d0   ! fast som
    rfac(9) = 0d0   ! slow som
    ! Assign combustion completeness to foliage
    cf(2) = pars(29) ! foliage
    ! Assign combustion completeness to non-photosynthetic
    cf(1) = pars(30) ; cf(3) = pars(30) ; cf(4) = pars(30)
    cf(8) = pars(31) ; cf(9) = pars(31) ! soil
    ! derived values for foliar, fine root and wood litter
    cf(5) = pars(32) ; cf(6) = pars(32) ; cf(7) = pars(33)

    !
    ! Begin looping through each time step
    !

    ! load some needed module level values
    mV%lai = POOLS(1,2)*sla
    mV%mint = met(2,1)  ! minimum temperature (oC)
    mV%maxt = met(3,1)  ! maximum temperature (oC)
    mV%swrad = met(4,1) ! incoming short wave radiation (MJ/m2/day)
    mV%co2 = met(5,1)   ! CO2 (ppm)
    mV%doy = met(6,1)   ! Day of year
    mV%rainfall = mV%rainfall_time(1) ! rainfall (kgH2O/m2/s)
    mV%wind_spd = met(15,1) ! wind speed (m/s)
    mV%vpd_kPa = met(16,1)*1d-3 ! vapour pressure deficit (Pa->kPa)
    mV%meant = (mV%mint + mV%maxt) * 0.5d0 ! mean air temperature (oC)
    mV%leafT = (mV%maxt*0.75d0) + (mV%mint*0.25d0)   ! initial day time canopy temperature (oC)
    mV%soilT = mV%meant
    mV%seconds_per_step = deltat(1) * seconds_per_day
    mV%days_per_step =  deltat(1)
    mV%days_per_step_1 =  mV%deltat_1(1)
    mV%dayl_seconds_1 = mV%daylength_seconds_1(1) 

    ! calculate some temperature dependent meteorologial properties
    call meteorological_constants(mV%leafT,mV%leafT+freeze,mV%vpd_kPa, mV)

    ! Initialise root reach based on initial coarse root biomass
    mV%fine_root_biomass = max(min_root,POOLS(1,3)*2d0)
    mV%root_biomass = mV%fine_root_biomass + max(min_root,POOLS(1,4)*pars(25)*2d0)
    ! calculate soil depth to which roots reach - needed here to set up
    ! layer_thickness correctly!
    mV%root_reach = mV%max_depth * mV%root_biomass / (mV%root_k + mV%root_biomass)
    ! Determine initial soil layer thickness
    mV%layer_thickness(1) = top_soil_depth ; mV%layer_thickness(2) = max(min_layer,mV%root_reach-top_soil_depth)
    mV%layer_thickness(3) = mV%max_depth - sum(mV%layer_thickness(1:2))
    mV%layer_thickness(4) = top_soil_depth
    mV%previous_depth = sum(mV%layer_thickness(1:2))
    ! Needed to initialise soils
    call calculate_Rtot(mV)
    call calculate_update_soil_water(transpiration,soilevaporation,snowsublimation,&
                                     0d0,FLUXES(1,29), mV) ! assume no evap or rainfall

    ! Store soil water content of the surface zone (mm)
    POOLS(1,11) = 1d3 * mV%soil_waterfrac(1) * mV%layer_thickness(1)

    do n = start, finish

       !!!!!!!!!!
       ! assign drivers and update some prognostic variables
       !!!!!!!!!!

       ! Incoming drivers
       mV%mint = met(2,n)  ! minimum temperature (oC)
       mV%maxt = met(3,n)  ! maximum temperature (oC)
       mV%swrad = met(4,n) ! incoming short wave radiation (MJ/m2/day)
       mV%co2 = met(5,n)   ! CO2 (ppm)
       mV%doy = met(6,n)   ! Day of year
       mV%rainfall = mV%rainfall_time(n)
       mV%meant = (mV%mint + mV%maxt) * 0.5d0 ! mean air temperature (oC)
       mV%leafT = (mV%meant + mV%maxt) * 0.5d0 ! estimate mean daytime air temperature (oC)
       mV%soilT = mV%meant ! Estimate mean day time soil temperature (oC)
       mV%wind_spd = met(15,n) ! wind speed (m/s)
       mV%vpd_kPa = met(16,n)*1d-3  ! Vapour pressure deficit (Pa -> kPa)
       mV%airt_zero_fraction = mV%airt_zero_fraction_time(n) ! fraction of above / below freezing temperature

       ! calculate LAI value
       mV%lai = POOLS(n,2)*sla
       DIAGS(n,1) = mV%lai

       ! extract timing related values
       mV%dayl_hours = mV%daylength_hours(n)
       mV%dayl_hours_fraction = mV%dayl_hours * 0.04166667d0 ! 1/24 = 0.04166667
       mV%dayl_seconds = mV%daylength_seconds(n) ; mV%dayl_seconds_1 = mV%daylength_seconds_1(n)
       mV%days_per_step = deltat(n) ; mV%days_per_step_1 = mV%deltat_1(n)
       mV%seconds_per_step = seconds_per_day * mV%days_per_step

       !!!!!!!!!!
       ! Adjust snow balance balance based on temperature
       !!!!!!!!!!

       ! snowing or not...?
       if (((mV%mint + mV%maxt) * 0.5d0) > 0d0) then
           ! on average above freezing so no snow
           mV%snowfall = 0d0
       else
           ! on average below freezing, so some snow based on proportion of temperature
           ! below freezing
           mV%snowfall = mV%rainfall * (1d0 - mV%airt_zero_fraction) ; mV%rainfall = mV%rainfall - mV%snowfall
           ! Add rainfall to the snowpack and clear rainfall variable
           mV%snow_storage = mV%snow_storage + (mV%snowfall*mV%seconds_per_step)
       end if

       ! melting or not...?
       if (mV%mint < 0d0 .and. mV%maxt > 0d0) then
           ! Also melt some of the snow based on airt_zero_fraction
           ! default assumption is that snow is melting at 10 % per day hour above freezing
           mV%snow_melt = min(mV%snow_storage, mV%airt_zero_fraction * mV%snow_storage * 0.1d0 * mV%days_per_step)
           mV%snow_storage = mV%snow_storage - mV%snow_melt
           ! adjust to rate for later addition to rainfall
           mV%snow_melt = mV%snow_melt / mV%seconds_per_step
       elseif (mV%maxt < 0d0) then
           mV%snow_melt = 0d0
       else if (mV%mint > 0d0 .and. mV%snow_storage > 0d0) then
           ! otherwise we assume snow is melting at 10 % per day above hour
           mV%snow_melt = min(mV%snow_storage, mV%snow_storage * 0.1d0 * mV%days_per_step)
           mV%snow_storage = mV%snow_storage - mV%snow_melt
           ! adjust to rate for later addition to rainfall
           mV%snow_melt = mV%snow_melt / mV%seconds_per_step
       else
           mV%snow_melt = 0d0
       end if
       DIAGS(n,2) = mV%snow_storage

       !!!!!!!!!!
       ! Calculate surface exchange coefficients
       !!!!!!!!!!

       ! calculate some temperature dependent meteorologial properties
       call meteorological_constants(mV%leafT,mV%leafT+freeze,mV%vpd_kPa, mV)
       ! pass variables from memory objects
       mV%convert_ms1_mmol_1 = mV%convert_ms1_mol_1 * 1d3
       ! calculate aerodynamic using consistent approach with SPA
       call calculate_aerodynamic_conductance(mV)
       ! Units converted from canopy top m/s to canopy scale (mmolH2O/m2ground/s)
       DIAGS(n,6) = mV%aerodynamic_conductance * mV%convert_ms1_mmol_1 * &
                    mV%leaf_canopy_wind_scaling
       DIAGS(n,15) = mV%leaf_canopy_wind_scaling ! canopy area scaling as a function of wind profiles

       !!!!!!!!!!
       ! Determine net shortwave and isothermal longwave energy balance
       !!!!!!!!!!

       call calculate_radiation_balance(mV)
       DIAGS(n,3) = mV%canopy_par_MJday ! Absorbed PAR by canopy (MJ/m2ground/day)
       DIAGS(n,14) = mV%leaf_canopy_light_scaling ! canopy area scaling as a function of light profiles
       DIAGS(n,13) = mV%soil_par_MJday

       !!!!!!!!!!
       ! Calculate physically constrained evaporation and
       ! soil water potential and total hydraulic resistance
       !!!!!!!!!!

       ! Canopy intercepted rainfall evaporation (kgH2O/m2/day)
       if (mV%lai > 0d0) then ! is this conditional needed?
           call calculate_wetcanopy_evaporation(wetcanopy_evap,mV%canopy_storage, mV)
       else
           ! reset pools
           mV%intercepted_rainfall = 0d0 ; mV%canopy_storage = 0d0 ; wetcanopy_evap = 0d0
       endif

       ! calculate the minimum soil & root hydraulic resistance based on total
       ! fine root mass ! *2*2 => *RS*C->Bio
       mV%fine_root_biomass = max(min_root,POOLS(n,3)*2d0)
       mV%root_biomass = mV%fine_root_biomass + max(min_root,POOLS(n,4)*pars(25)*2d0)
       call calculate_Rtot(mV)
       ! Pass root~water~soil variables to output variable
       DIAGS(n,8) = mV%root_reach ! Rooting depth (m)
       DIAGS(n,10) = mV%wSWP      ! Soil water potential weighted by supply of water
       DIAGS(n,11) = mV%rSWP      ! Soil water potential weighted by access water
       DIAGS(n,12) = mV%Reff      ! Effective hydraulic resistance MPa.s.m2.mmol-1 H20

       ! calculate radiation absorption and estimate stomatal conductance
       call calculate_stomatal_conductance(mV)
       ! Estimate stomatal conductance relative to its minimum / maximum, i.e. how
       ! close are we to maxing out supply
       DIAGS(n,7) = (mV%stomatal_conductance  - mV%minimum_conductance) &
                  / (mV%potential_conductance - mV%minimum_conductance)
       ! Store the canopy level stomatal conductance (mmolH2O/m2ground/s)
       DIAGS(n,5) = mV%stomatal_conductance

       ! Note that soil mass balance will be calculated after phenology
       ! adjustments

       ! Reset output variable
       if (mV%stomatal_conductance > vsmall) then
           ! Gross primary productivity (gC/m2/day)
           ! Assumes acm_gpp_stage_1 ran as part of stomatal conductance calculation
           FLUXES(n,1) = acm_gpp_stage_2(mV%stomatal_conductance, mV) * umol_to_gC * mV%dayl_seconds
           ! Estimate the ratio of leaf internal to ambient CO2 concentrations
           DIAGS(n,4) = mV%ci / mV%co2
           ! Canopy transpiration (kgH2O/m2/day)
           call calculate_transpiration(transpiration, mV)
           ! restrict transpiration to positive only
           transpiration = max(0d0,transpiration)
       else
           ! assume zero fluxes
           FLUXES(n,1) = 0d0 ; transpiration = 0d0 ; DIAGS(n,4) = 0d0
       endif

       ! Estimate average leaf water potential (MPa) based on effective hydraulic resistance, wSWP and transpiration.
       ! Positive LWPs can be estimated given very small gs and cold temperatures.
       ! Debugging print statements
       !print*,"Estimate LWP"
       !LWP = SWP(1:nos_root_layers) - head*canopy_height &
       !     - transpiration*uptake_fraction(1:nos_root_layers) &
       !     * (dayl_seconds_1/mmol_to_kg_water)/Rcond_layer(1:nos_root_layers)
       DIAGS(n,9) =  min(0d0, mV%wSWP - (head*canopy_height) - (((transpiration*mV%dayl_seconds_1)/mmol_to_kg_water) * mV%Reff))

       ! temprate (i.e. temperature modified rate of metabolic activity))
       FLUXES(n,2) = exp(pars(10)*0.5d0*(met(3,n)+met(2,n)))

       ! Update the microbial activity and mortality rate
       call update_microbial_activity_and_death(pars(40),pars(46),pars(45),mV%days_per_step, &
                                                FLUXES(n,2),POOLS(n,8),microbial_activity,microbial_death)
       DIAGS(n,16) = microbial_death ; DIAGS(n,17) = microbial_activity
!if (DIAGS(n,16) /= DIAGS(n,16)) then
!print*,"D16 ",DIAGS(n,16), DIAGS(n,17), FLUXES(n,2),POOLS(n,8)
!stop
!end if
!if (DIAGS(n,17) /= DIAGS(n,17)) then
!print*,"D17 ",DIAGS(n,16), DIAGS(n,17), FLUXES(n,2),POOLS(n,8)
! !D17          0.00E+000    NaN          1.036       Infinity
!stop
!end if
       ! 
       ! Allocate photosynthate
       ! 

       ! Estimate the maintenance respiration component of autotrophic respiration (gC.m-2.day-1)
       FLUXES(n,3) = pars(2)*FLUXES(n,1)
       ! leaf production rate (gC.m-2.day-1)
       FLUXES(n,4) = (FLUXES(n,1)-FLUXES(n,3))*pars(3)
       ! labile production (gC.m-2.day-1)
       FLUXES(n,5) = (FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,4))*pars(13)
       ! root production (gC.m-2.day-1)
       FLUXES(n,6) = (FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,4)-FLUXES(n,5))*pars(4)
       ! wood production
       FLUXES(n,7) = FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,4)-FLUXES(n,5)-FLUXES(n,6)

       ! Labile release and leaffall factors
       FLUXES(n,9)  = ff_wf*exp(-(sin((mV%doy-pars(15)+osf)/sf)*sf_over_wf)**2d0)
       FLUXES(n,16) = fl_wl*exp(-(sin((mV%doy-pars(12)+osl)/sf)*sf_over_wl)**2d0)

       !
       ! those with time dependancies
       !

       ! total labile release
       FLUXES(n,8) = POOLS(n,1)*(1d0-(1d0-FLUXES(n,16))**mV%days_per_step)*mV%days_per_step_1
       ! total leaf litter production
       FLUXES(n,10) = POOLS(n,2)*(1d0-(1d0-FLUXES(n,9))**mV%days_per_step)*mV%days_per_step_1
       ! total wood production
       FLUXES(n,11) = POOLS(n,4)*(1d0-(1d0-pars(6))**mV%days_per_step)*mV%days_per_step_1
       ! total root litter production
       FLUXES(n,12) = POOLS(n,3)*(1d0-(1d0-pars(7))**mV%days_per_step)*mV%days_per_step_1

       !
       ! those with temperature AND time dependancies
       !

       ! Turnover of foliage litter 
       FLUXES(n,13) = POOLS(n,5)*(1d0-(1d0-FLUXES(n,2)*pars(8))**mV%days_per_step)*mV%days_per_step_1
       ! Turnover of wood litter
       FLUXES(n,30) = POOLS(n,7)*(1d0-(1d0-FLUXES(n,2)*pars(35))**mV%days_per_step)*mV%days_per_step_1
       ! Turnover of fine root litter
       FLUXES(n,56) = POOLS(n,6)*(1d0-(1d0-FLUXES(n,2)*pars(9))**mV%days_per_step)*mV%days_per_step_1

       ! Decomposition of foliage litter to som and update heterotrophic flux
       FLUXES(n,15) = FLUXES(n,13) * pars(1) ; FLUXES(n,13) = FLUXES(n,13) - FLUXES(n,15)
       ! Decomposition of wood litter to som and update heterotrophic flux
       FLUXES(n,31) = FLUXES(n,30) * pars(1) ; FLUXES(n,30) = FLUXES(n,30) - FLUXES(n,31)
       ! Decomposition of fine root litter to som and update heterotrophic flux
       FLUXES(n,57) = FLUXES(n,56) * pars(1) ; FLUXES(n,56) = FLUXES(n,56) - FLUXES(n,57)

       ! Soil and microbial dynamics are estimated mostly as absolute fluxes via functions of the pools themselves. 
       ! To reduce the non-linear effects and associated mass balance challenges, the absolute fluxes are converted 
       ! into a fraction of the target pool and then the actual flux is recalculated using our standard compound 
       ! interest approach.

       ! Respiration heterotrophic slow som
       FLUXES(n,14) = max(0d0,min(1d0,(POOLS(n,10)*(1d0-pars(48))*microbial_activity*pars(49)) / max(POOLS(n,9),vsmall)))
       FLUXES(n,14) = POOLS(n,9) * (1d0-(1d0-FLUXES(n,14))**mV%days_per_step)*mV%days_per_step_1
       ! Respiration heterotrophic fast som
       FLUXES(n,58) = max(0d0,min(1d0,(POOLS(n,10)*(1d0-pars(44))*pars(50)*microbial_activity)))
       FLUXES(n,58) = POOLS(n,8) * (1d0-(1d0-FLUXES(n,58))**mV%days_per_step)*mV%days_per_step_1
       ! Respiration heterotrophic microbial
       FLUXES(n,59) = POOLS(n,10) * (1d0-(1d0-(pars(47)*microbial_activity))**mV%days_per_step)*mV%days_per_step_1
       ! Microbial death allocation to slow som
       FLUXES(n,60) = max(0d0,min(1d0,(FLUXES(n,2)*microbial_death*microbial_activity)))
       FLUXES(n,60) = POOLS(n,10) * (1d0-(1d0-FLUXES(n,60))**mV%days_per_step)*mV%days_per_step_1
       ! Microbial mediated transfer of carbon from slow to fast
       FLUXES(n,61) = POOLS(n,10) * (1d0-(1d0-(pars(48)*microbial_activity*pars(49)))**mV%days_per_step)*mV%days_per_step_1
       ! Accumulation of fast som into microbial carbon
       FLUXES(n,62) = max(0d0,min(1d0,(POOLS(n,10)*(pars(44)*pars(50)*microbial_activity))))
       FLUXES(n,62) = POOLS(n,8) *  (1d0-(1d0-FLUXES(n,62))**mV%days_per_step)*mV%days_per_step_1
!if (FLUXES(n,62) < 0d0) print*,"62 ", FLUXES(n,62), POOLS(n,8),POOLS(n,10),microbial_activity,pars(44),pars(50),mV%days_per_step,mV%days_per_step_1

       !!!!!!!!!!
       ! calculate growth respiration and adjust allocation to pools assuming
       ! 0.21875 of total C allocation towards each pool (i.e. 0.28 .eq. xNPP)
       !!!!!!!!!!

       ! Determine growth respiration from labile to foliage.
       ! Update Ra flux and then correct the flux making it to the foliage pool
       FLUXES(n,55) = FLUXES(n,8)*Rg_fraction ; FLUXES(n,3)  = FLUXES(n,3) + FLUXES(n,55)
       FLUXES(n,8) = FLUXES(n,8) * one_Rg_fraction
       ! Direct GPP allocation to foliage
       FLUXES(n,3) = FLUXES(n,3) + (FLUXES(n,4)*Rg_fraction) ; FLUXES(n,4) = FLUXES(n,4) * one_Rg_fraction
       ! Direct GPP allocation to fine roots
       FLUXES(n,3) = FLUXES(n,3) + (FLUXES(n,6)*Rg_fraction) ; FLUXES(n,6) = FLUXES(n,6) * one_Rg_fraction
       ! Direct GPP allocation to wood 
       FLUXES(n,3) = FLUXES(n,3) + (FLUXES(n,7)*Rg_fraction) ; FLUXES(n,7) = FLUXES(n,7) * one_Rg_fraction

       !
       ! update pools for next timestep
       !

       ! labile pool
       POOLS(n+1,1)  = POOLS(n,1)  + (FLUXES(n,5)-FLUXES(n,8)-FLUXES(n,55))*mV%days_per_step
       ! foliar pool
       POOLS(n+1,2)  = POOLS(n,2)  + (FLUXES(n,4)+FLUXES(n,8)-FLUXES(n,10))*mV%days_per_step
       ! root pool
       POOLS(n+1,3)  = POOLS(n,3)  + (FLUXES(n,6)-FLUXES(n,12))*mV%days_per_step
       ! wood pool
       POOLS(n+1,4)  = POOLS(n,4)  + (FLUXES(n,7)-FLUXES(n,11))*mV%days_per_step
       ! foliar litter pool
       POOLS(n+1,5)  = POOLS(n,5)  + (FLUXES(n,10)-FLUXES(n,13)-FLUXES(n,15))*mV%days_per_step
       ! fine root litter pool
       POOLS(n+1,6)  = POOLS(n,6)  + (FLUXES(n,12)-FLUXES(n,56)-FLUXES(n,57))*mV%days_per_step
       ! wood litter
       POOLS(n+1,7)  = POOLS(n,7)  + (FLUXES(n,11)-FLUXES(n,30)-FLUXES(n,31))*mV%days_per_step
       ! fast som pool
       POOLS(n+1,8)  = POOLS(n,8)  + ( (FLUXES(n,15)*(1d0-pars(41)))+ & ! Foliage litter
                                       (FLUXES(n,31)*(1d0-pars(43)))+ & ! Wood litter 
                                       (FLUXES(n,57)*(1d0-pars(42)))+ & ! Fine root litter
                                       FLUXES(n,61)-FLUXES(n,58)-FLUXES(n,62))*mV%days_per_step
       ! slow som pool
       POOLS(n+1,9)  = POOLS(n,9)  + ( (FLUXES(n,15)*pars(41))+ &
                                       (FLUXES(n,31)*pars(43))+ &
                                       (FLUXES(n,57)*pars(42))+ &
                                       FLUXES(n,60)-FLUXES(n,14)-FLUXES(n,61))*mV%days_per_step
       ! microbial pool
       POOLS(n+1,10) = POOLS(n,10) + (FLUXES(n,62)-FLUXES(n,59)-FLUXES(n,60))*mV%days_per_step

!if (POOLS(n+1,8) /= POOLS(n+1,8) .or. POOLS(n+1,8) > 1d6) then
!print*,"P8  ",POOLS(n+1,8),POOLS(n,8),POOLS(n+1,9),POOLS(n,9),POOLS(n+1,10),POOLS(n,10) 
!print*,"    ",FLUXES(n,15),FLUXES(n,31),FLUXES(n,57),FLUXES(n,61)
!print*,"    ",FLUXES(n,58),FLUXES(n,62),FLUXES(n,60),FLUXES(n,14),FLUXES(n,61),FLUXES(n,59)
!stop
!end if
!if (POOLS(n+1,9) /= POOLS(n+1,9)) then
!print*,"P9  ",POOLS(n+1,8),POOLS(n,8),POOLS(n+1,9),POOLS(n,9),POOLS(n+1,10),POOLS(n,10)
!print*,"    ",FLUXES(n,15),FLUXES(n,31),FLUXES(n,57),FLUXES(n,61)
!print*,"    ",FLUXES(n,58),FLUXES(n,62),FLUXES(n,60),FLUXES(n,14),FLUXES(n,61),FLUXES(n,59)
!stop
!end if
!if (POOLS(n+1,10) /= POOLS(n+1,10)) then
!print*,"P10 ",POOLS(n+1,8),POOLS(n,8),POOLS(n+1,9),POOLS(n,9),POOLS(n+1,10),POOLS(n,10) 
!print*,"    ",FLUXES(n,15),FLUXES(n,31),FLUXES(n,57),FLUXES(n,61)
!print*,"    ",FLUXES(n,58),FLUXES(n,62),FLUXES(n,60),FLUXES(n,14),FLUXES(n,61),FLUXES(n,59)
!stop
!end if

       ! Enforce mass reality - note this breaks mass balance as fluxes are not proportionally updated
       POOLS(n+1,8)  = max(vsmall,POOLS(n+1,8))
       POOLS(n+1,9)  = max(vsmall,POOLS(n+1,9))
       POOLS(n+1,10) = max(vsmall,POOLS(n+1,10))

       !!!!!!!!!!
       ! Update soil water balance
       !!!!!!!!!!

       ! add any snow melt to the rainfall now that we have already dealt with the canopy interception
       mV%rainfall = mV%rainfall + mV%snow_melt
       ! do mass balance (i.e. is there enough water to support ET)
       call calculate_update_soil_water(transpiration,soilevaporation,snowsublimation, &
                                        ((mV%rainfall-mV%intercepted_rainfall)*seconds_per_day) &
                                       ,FLUXES(n,29), mV)
       ! now that soil mass balance has been updated we can add the wet canopy
       ! evaporation (kgH2O.m-2.day-1)
       FLUXES(n,29) = FLUXES(n,29) + wetcanopy_evap
       ! store soil water content of the surface zone (mm)
       POOLS(n+1,11) = 1d3 * mV%soil_waterfrac(1) * mV%layer_thickness(1)
       ! Assign all water variables to output variables (kgH2O/m2/day)
       FLUXES(n,46) = transpiration   ! transpiration (kgH2O/m2/day)
       FLUXES(n,47) = soilevaporation ! soil evaporation (kgH2O/m2/day)
       FLUXES(n,48) = wetcanopy_evap  ! wet canopy evaporation (kgH2O/m2/day)
       FLUXES(n,49) = mV%runoff          ! soil surface runoff (kgH2O/m2/day)
       FLUXES(n,50) = mV%underflow       ! drainage from bottom of soil column (kgH2O/m2/day)
       FLUXES(n,51) = mV%water_grav_flow(1) ! drainage from the surface soil layer to 2nd (kgH2O/m2/day)
       FLUXES(n,52) = mV%infiltrated(1)  ! top soil surface infiltration by rain (kgH2O/m2/day)
       FLUXES(n,66) = mV%infiltrated(2)  ! middle soil surface infiltration by rain (kgH2O/m2/day)
       FLUXES(n,67) = mV%infiltrated(3)  ! bottom soil surface infiltration by rain (kgH2O/m2/day)       
       FLUXES(n,53) = mV%uptake_fraction(1) ! transpiration fraction extracted from 1st rooting layer (the soil surface)
       FLUXES(n,54) = mV%uptake_fraction(2) ! transpiration fraction extracted from 2nd rooting layer (dynamic 2nd layer)       

       !!!!!!!!!!
       ! Extract biomass - e.g. deforestation / degradation
       !!!!!!!!!!

       ! reset values
       harvest_management = 0 ; burnt_area = 0d0

       ! Does harvest activities occur?
       if (met(8,n) > 0d0) then

           ! Load the management type / scenario into local variable
           harvest_management = int(met(13,n))

           ! Determine the fraction of cut labile C which remains in system as residue.
           ! We assume that labile is proportionally distributed through the plants
           ! root and wood (structural C).
           C_total = POOLS(n+1,3) + POOLS(n+1,4)
           ! Ensure there is available C for extraction
           if (C_total > 0d0) then
               ! Harvest activities on the wood / structural pool varies depending on
               ! whether it is above or below ground. As such, partition the wood pool
               ! between above ground stem(+branches) and below ground coarse root.
               Crootcr = POOLS(n+1,4)*Crootcr_part(harvest_management)
               Cstem   = POOLS(n+1,4)-Crootcr
               ! Calculate the fraction of harvested labile which remains in system as residue
               labile_frac_res = ((POOLS(n+1,3)/C_total) * roots_frac_res(harvest_management)  ) &
                               + ((Cstem/C_total)        * stem_frac_res(harvest_management)   ) &
                               + ((Crootcr/C_total)      * rootcr_frac_res(harvest_management) )
               ! Calculate the management scenario specific resistance fraction
               labile_frac_removal = ((POOLS(n+1,3)/C_total) * roots_frac_removal(harvest_management)  ) &
                                     + ((Cstem/C_total)        * 1d0   ) &
                                     + ((Crootcr/C_total)      * rootcr_frac_removal(harvest_management) )

               ! Calculate the total loss from biomass pools
               ! We assume that fractional clearing always equals the fraction
               ! of foliage and above ground (stem) wood removal. However, we assume
               ! that coarse root and fine root extractions are dependent on the
               ! management activity type, e.g. in coppice below ground remains.
               ! Thus, labile extractions are also dependent.
               labile_loss = POOLS(n+1,1) * labile_frac_removal * met(8,n)
               foliar_loss = POOLS(n+1,2) * met(8,n)
               roots_loss  = POOLS(n+1,3) * roots_frac_removal(harvest_management) * met(8,n)
               stem_loss   = (Cstem * met(8,n))
               rootcr_loss = (Crootcr * rootcr_frac_removal(harvest_management) * met(8,n))
               wood_loss   =  stem_loss + rootcr_loss

               ! Transfer fraction of harvest waste to litter, wood litter or som pools.
               ! This includes explicit calculation of the stem and coarse root residues due
               ! to their potentially different treatments under management scenarios
               labile_residue = labile_loss*labile_frac_res
               foliar_residue = foliar_loss*foliage_frac_res(harvest_management)
               roots_residue  = roots_loss*roots_frac_res(harvest_management)
               coarse_root_residue = rootcr_loss*rootcr_frac_res(harvest_management)
               stem_residue = stem_loss*stem_frac_res(harvest_management)
               wood_residue = stem_residue + coarse_root_residue
               ! Mechanical loss of Csom due to coarse root extraction,
               ! less the loss remaining as residue
               soil_loss_with_roots = (rootcr_loss-coarse_root_residue) &
                                    * soil_loss_frac(harvest_management)

               ! Update pools
               POOLS(n+1,1) = POOLS(n+1,1) - labile_loss          ! labile
               POOLS(n+1,2) = POOLS(n+1,2) - foliar_loss          ! foliage
               POOLS(n+1,3) = POOLS(n+1,3) - roots_loss           ! fine roots
               POOLS(n+1,4) = POOLS(n+1,4) - wood_loss            ! wood
               POOLS(n+1,5) = POOLS(n+1,5) + foliar_residue       ! foliar litter
               POOLS(n+1,6) = POOLS(n+1,6) + roots_residue        ! root litter
               POOLS(n+1,7) = POOLS(n+1,7) + wood_residue         ! wood litter
               POOLS(n+1,8) = POOLS(n+1,8) + labile_residue       ! fast som
               POOLS(n+1,9) = POOLS(n+1,9) - soil_loss_with_roots ! slow som
               ! mass balance check
               where (POOLS(n+1,1:9) < 0d0) POOLS(n+1,1:9) = 0d0

               ! Convert harvest related extractions to daily rate for output
               ! For dead organic matter pools, in most cases these will be zeros.
               ! But these variables allow for subseqent management where surface litter
               ! pools are removed (i.e. we make no provision for root litter removal) 
               ! or mechanical extraction from soil occurs.

               FLUXES(n,34) = (labile_loss-labile_residue) * mV%days_per_step_1
               FLUXES(n,35) = (foliar_loss-foliar_residue) * mV%days_per_step_1
               FLUXES(n,36) = (roots_loss-roots_residue) * mV%days_per_step_1
               FLUXES(n,37) = (wood_loss-wood_residue) * mV%days_per_step_1
               FLUXES(n,38) = 0d0 ! foliage litter extraction by harvest
               FLUXES(n,39) = 0d0 ! wood litter extraction by harvest
               FLUXES(n,40) = soil_loss_with_roots * mV%days_per_step_1
               ! Convert harvest related residue generations to daily rate for output
               FLUXES(n,41) = labile_residue * mV%days_per_step_1
               FLUXES(n,42) = foliar_residue * mV%days_per_step_1
               FLUXES(n,43) = roots_residue * mV%days_per_step_1
               FLUXES(n,44) = wood_residue * mV%days_per_step_1

               ! Total C extraction, this includes any som or litter clearing
               FLUXES(n,45) = sum(FLUXES(n,34:44))

           end if ! C_total > 0d0

       endif ! end deforestation info

       !!!!!!!!!!
       ! Impose fire
       !!!!!!!!!!

       if (met(9,n) > 0d0 .or.(met(8,n) > 0d0 .and. harvest_management > 0)) then

           ! Adjust burnt area to account for the managment decisions which may not be
           ! reflected in the burnt area drivers
           burnt_area = met(9,n)
           if (met(8,n) > 0d0 .and. burnt_area > 0d0) then
               ! pass harvest management to local integer
               burnt_area = min(1d0,burnt_area + post_harvest_burn(harvest_management))
           else if (met(8,n) > 0d0 .and. burnt_area <= 0d0) then
               burnt_area = post_harvest_burn(harvest_management)
           endif

           ! Determine the corrected burnt area
           if (burnt_area > 0d0) then

               ! first calculate combustion / emissions fluxes in g C m-2 d-1
               FLUXES(n,18) = POOLS(n+1,1)*burnt_area*cf(1) * mV%days_per_step_1 ! labile
               FLUXES(n,19) = POOLS(n+1,2)*burnt_area*cf(2) * mV%days_per_step_1 ! foliar
               FLUXES(n,20) = POOLS(n+1,3)*burnt_area*cf(3) * mV%days_per_step_1 ! roots
               FLUXES(n,21) = POOLS(n+1,4)*burnt_area*cf(4) * mV%days_per_step_1 ! wood
               FLUXES(n,22) = POOLS(n+1,5)*burnt_area*cf(5) * mV%days_per_step_1 ! foliar litter
               FLUXES(n,63) = POOLS(n+1,6)*burnt_area*cf(6) * mV%days_per_step_1 ! root litter               
               FLUXES(n,32) = POOLS(n+1,7)*burnt_area*cf(7) * mV%days_per_step_1 ! wood litter
               FLUXES(n,64) = POOLS(n+1,8)*burnt_area*cf(8) * mV%days_per_step_1 ! fast som
               FLUXES(n,23) = POOLS(n+1,9)*burnt_area*cf(9) * mV%days_per_step_1 ! slow som

               ! second calculate litter transfer fluxes in g C m-2 d-1, all pools except som
               FLUXES(n,24) = POOLS(n+1,1)*burnt_area*(1d0-cf(1))*(1d0-rfac(1)) * mV%days_per_step_1 ! labile into fast som
               FLUXES(n,25) = POOLS(n+1,2)*burnt_area*(1d0-cf(2))*(1d0-rfac(2)) * mV%days_per_step_1 ! foliar into litter
               FLUXES(n,26) = POOLS(n+1,3)*burnt_area*(1d0-cf(3))*(1d0-rfac(3)) * mV%days_per_step_1 ! roots into litter
               FLUXES(n,27) = POOLS(n+1,4)*burnt_area*(1d0-cf(4))*(1d0-rfac(4)) * mV%days_per_step_1 ! wood into litter
               FLUXES(n,28) = POOLS(n+1,5)*burnt_area*(1d0-cf(5))*(1d0-rfac(5)) * mV%days_per_step_1 ! foliar litter into som
               FLUXES(n,65) = POOLS(n+1,6)*burnt_area*(1d0-cf(6))*(1d0-rfac(6)) * mV%days_per_step_1 ! root litter into som
               FLUXES(n,33) = POOLS(n+1,7)*burnt_area*(1d0-cf(7))*(1d0-rfac(7)) * mV%days_per_step_1 ! wood litter into som

               ! update pools - first remove burned vegetation
               POOLS(n+1,1) = POOLS(n+1,1) - (FLUXES(n,18) + FLUXES(n,24)) * mV%days_per_step ! labile
               POOLS(n+1,2) = POOLS(n+1,2) - (FLUXES(n,19) + FLUXES(n,25)) * mV%days_per_step ! foliar
               POOLS(n+1,3) = POOLS(n+1,3) - (FLUXES(n,20) + FLUXES(n,26)) * mV%days_per_step ! roots
               POOLS(n+1,4) = POOLS(n+1,4) - (FLUXES(n,21) + FLUXES(n,27)) * mV%days_per_step ! wood
               ! update pools - add litter transfer
               POOLS(n+1,5) = POOLS(n+1,5) + (FLUXES(n,25) - FLUXES(n,22) - FLUXES(n,28)) * mV%days_per_step
               POOLS(n+1,6) = POOLS(n+1,6) + (FLUXES(n,26) - FLUXES(n,63)) * mV%days_per_step
               POOLS(n+1,7) = POOLS(n+1,7) + (FLUXES(n,27) - FLUXES(n,32)) * mV%days_per_step
               POOLS(n+1,8) = POOLS(n+1,8) + (FLUXES(n,24) + &
                                              (FLUXES(n,28)*(1d0-pars(41))) + &
                                              (FLUXES(n,65)*(1d0-pars(42))) + &
                                              (FLUXES(n,33)*(1d0-pars(43))) - FLUXES(n,64)) * mV%days_per_step
               POOLS(n+1,9) = POOLS(n+1,9) + ((FLUXES(n,28)*pars(41)) + &
                                              (FLUXES(n,65)*pars(42)) + &
                                              (FLUXES(n,33)*pars(43)) - FLUXES(n,23)) * mV%days_per_step

               ! calculate ecosystem emissions
               FLUXES(n,17) = sum(FLUXES(n,18:23)) + FLUXES(n,32)+FLUXES(n,63)+FLUXES(n,64)

           end if ! burnt_area > 0
       else
           ! set fluxes to zero
           FLUXES(n,17:28) = 0d0 ; FLUXES(n,32:33) = 0d0 ; FLUXES(n,63:65) = 0d0
       end if

    end do ! nodays loop

  end subroutine CARBON_MODEL
  !
  !------------------------------------------------------------------
  !
  subroutine acm_gpp_stage_1(mV)
    use carbon_model_memory

    ! Estimate the light and temperature limited photosynthesis components.
    ! See acm_gpp_stage_2() for estimation of CO2 supply limitation and
    ! combination of light, temperature and CO2 co-limitation

    implicit none

    ! Argument
    type(model_working_variables) :: mV

    ! Declare local variables
    double precision :: a, b, c, Pl_max, PAR_m2, airt_ad

    !
    ! Metabolic limited photosynthesis
    !

    ! maximum rate of temperature and nitrogen (canopy efficiency) limited
    ! photosynthesis (gC.m-2.day-1 -> umolC/m2/s). Scaling from leaf to canopy
    ! scaled assumed to follow integral of light environment.
    mV%metabolic_limited_photosynthesis = gC_to_umol*mV%leaf_canopy_light_scaling*mV%ceff*seconds_per_day_1 &
                                     * ((mV%leafT - Vc_minT) / ((mV%leafT - Vc_minT) + Vc_coef))         &
                                     * opt_max_scaling(pn_max_temp,pn_min_temp,pn_opt_temp,pn_kurtosis,mV%leafT)

    !
    ! Light limited photosynthesis
    !

    ! Calculate light limted rate of photosynthesis (umolC.m-2.s-1, daylight) as a function
    ! light capture and leaf to canopy scaling on quantum yield (e0).
    !light_limited_photosynthesis = e0 * canopy_par_MJday
    mV%light_limited_photosynthesis = e0 * mV%canopy_par_MJday * mV%dayl_seconds_1 * gC_to_umol

    !
    ! Stomatal conductance independent variables for diffusion limited
    ! photosynthesis
    !

    ! Canopy level boundary layer conductance unit change
    ! (m.s-1 -> mol.m-2.s-1) assuming sea surface pressure only.
    ! Note the ratio of H20:CO2 diffusion through leaf level boundary layer is
    ! 1.37 (Jones appendix 2). Note conversion to resistance for easiler merging
    ! with stomatal conductance in acm_gpp_stage_2).
    mV%rb_mol_1 = (mV%aerodynamic_conductance * mV%convert_ms1_mol_1 * gb_H2O_CO2 * &
                   mV%leaf_canopy_wind_scaling) ** (-1d0)

    ! Arrhenious Temperature adjustments for Michaelis-Menten coefficients
    ! for CO2 (kc) and O2 (ko) and CO2 compensation point
    ! See McMurtrie et al., (1992) Australian Journal of Botany, vol 40, 657-677
    mV%co2_half_sat   = arrhenious(kc_half_sat_25C,kc_half_sat_gradient,mV%leafT)
    mV%co2_comp_point = arrhenious(co2comp_sat_25C,co2comp_gradient,mV%leafT)

    ! don't forget to return
    return

  end subroutine acm_gpp_stage_1
  !
  !------------------------------------------------------------------
  !
  double precision function acm_gpp_stage_2(gs, mV)
    use carbon_model_memory

    ! Combine the temperature (pn) and light (pl) limited gross primary productivity
    ! estimates with CO2 supply limited via stomatal conductance (gs).
    ! See acm_gpp_stage_1() for additional details on pn and pl calculation.

    implicit none

      type(model_working_variables) :: mV

    ! declare input variables
    double precision, intent(in) :: gs

    ! declare local variables
    double precision :: pp, qq, mult, rc, pd

    !
    ! Combined diffusion limitation and carboxylation limited photosynthesis
    !

    ! Estimation of ci is based on the assumption that metabilic limited
    ! photosynthesis is equal to diffusion limited. For details
    ! see Williams et al, (1997), Ecological Applications,7(3), 1997, pp. 882–894

    ! Daily canopy conductance dertmined through combination of aerodynamic and
    ! stomatal conductances. Both conductances are scaled to canopy aggregate.
    ! aerodynamic conductance already in units of molCO2.m-2.s-1 (see acm_gpp_stage_1).
    ! Stomatal conductance scaled from mmolH2O to molCO2.
    ! The ratio of H20:CO2 diffusion is 1.646259 (Jones appendix 2).
    !
    ! Combining in series the stomatal and boundary layer conductances
    ! to make canopy resistence (s/m2/molCO2)
    rc = (gs*gs_H2Ommol_CO2mol) ** (-1d0) + mV%rb_mol_1

    ! pp and qq represent limitation by metabolic (temperature & N) and
    ! diffusion (co2 supply) respectively
    pp = mV%metabolic_limited_photosynthesis*rc ; qq = mV%co2_comp_point-mV%co2_half_sat
    mult = mV%co2+qq-pp
    ! calculate internal CO2 concentration (ppm or umol/mol)
    mV%ci = 0.5d0*(mult+sqrt((mult*mult)-4d0*(mV%co2*qq-pp*mV%co2_comp_point)))

    ! calculate CO2 limited rate of photosynthesis (umolC.m-2.s-1)
    ! Then scale to day light period as this is then consistent with the light
    ! capture period (1/24 = 0.04166667)
    pd = ((mV%co2-mV%ci)/rc)

    !
    ! Estimate CO2 and light co-limitation
    !

    ! calculate combined light and CO2 limited photosynthesis (umolC/m2/s)
    acm_gpp_stage_2 = mV%light_limited_photosynthesis*pd/(mV%light_limited_photosynthesis+pd)

    ! Estimate ci as a function of the final combined GPP estimate
    pp = acm_gpp_stage_2*rc ; mult = mV%co2+qq-pp
    ! calculate internal CO2 concentration (ppm or umol/mol)
    mV%ci = 0.5d0*(mult+sqrt((mult*mult)-4d0*(mV%co2*qq-pp*mV%co2_comp_point)))

    ! sanity check
    if (acm_gpp_stage_2 /= acm_gpp_stage_2) then
        acm_gpp_stage_2 = 0d0 ; mV%ci = 0d0
    end if

    ! don't forget to return
    return

  end function acm_gpp_stage_2
  !
  !----------------------------------------------------------------------
  !
  double precision function find_gs_iWUE(gs_in, mV)
    use carbon_model_memory

    ! Calculate CO2 limited photosynthesis as a function of metabolic limited
    ! photosynthesis (pn), atmospheric CO2 concentration and stomatal
    ! conductance (gs_in). Photosynthesis is calculated twice to allow for
    ! testing of senstivity to iWUE.

    ! arguments
    double precision, intent(in) :: gs_in

      type(model_working_variables), optional :: mV

    !!!!!!!!!!
    ! Optimise intrinsic water use efficiency
    !!!!!!!!!!

    ! Determine impact of gs increment on pd and how far we are from iWUE
    find_gs_iWUE = mV%iWUE_step - (acm_gpp_stage_2(gs_in + mV%delta_gs, mV) - acm_gpp_stage_2(gs_in, mV))

    ! Remember to return back to the user
    return

  end function find_gs_iWUE
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_stomatal_conductance (mV)
    use carbon_model_memory

    use brent_zero, only: zbrent

    ! Determines an approximation of canopy scale stomatal conductance (gc)
    ! mmolH2O.m-2.s-1 based on potential hydraulic flow, air temperature and absorbed radiation.

    implicit none

      type(model_working_variables) :: mV

    ! local variables
    double precision :: denom, iWUE_upper!, iWUE_lower
    double precision, parameter :: max_gs = 1000d0, &  ! mmolH2O.m-2.s-1 (leaf area)
                                   min_gs = 0.01d0, &  ! mmolH2O.m-2.s-1 (leaf area)
                                   tol_gs = 0.01d0     ! mmolH2O.m-2.s-1 (leaf area)

    !!!!!!!!!!
    ! Calculate stomatal conductance under H2O and CO2 limitations
    !!!!!!!!!!

    if (mV%aerodynamic_conductance > vsmall .and. mV%total_water_flux > vsmall .and. &
        mV%leafT > Vc_minT .and. mV%leaf_canopy_light_scaling > vsmall) then

        ! Pass minimum conductance from local parameter to global value
        mV%minimum_conductance = min_gs * mV%leaf_canopy_light_scaling

        ! Invert Penman-Monteith equation to give gs (m.s-1) needed to meet
        ! maximum possible evaporation for the day.
        ! This will then be reduced based on CO2 limits for diffusion based
        ! photosynthesis
        denom = mV%slope * (((mV%canopy_swrad_MJday * 1d6 * mV%dayl_seconds_1) + mV%canopy_lwrad_Wm2)) &
              + (mV%ET_demand_coef * mV%aerodynamic_conductance * mV%leaf_canopy_wind_scaling)
        denom = (denom / (mV%lambda * mV%total_water_flux * mmol_to_kg_water)) - mV%slope
        mV%potential_conductance = (mV%aerodynamic_conductance * mV%leaf_canopy_wind_scaling) / (denom / mV%psych)

        ! convert m.s-1 to mmolH2O.m-2.d-1, per unit ground area, note that this
        ! is implicitly the canopy scaled value
        mV%potential_conductance = mV%potential_conductance * mV%convert_ms1_mmol_1
        ! if conditions are dew forming then set conductance to maximum as we
        ! are not going to be limited by water demand
        if (mV%potential_conductance <= 0d0 .or. mV%potential_conductance > max_gs*mV%leaf_canopy_light_scaling) then
            mV%potential_conductance = max_gs*mV%leaf_canopy_light_scaling
        end if

        ! If there is a positive demand for water then we will solve for
        ! photosynthesis limits on gs through iterative solution

        ! Determine the appropriate canopy scaled gs increment and return threshold
        mV%delta_gs = 1d0 * mV%leaf_canopy_light_scaling ! mmolH2O/m2leaf/s
        mV%iWUE_step = iWUE * mV%leaf_canopy_light_scaling ! umolC/mmolH2Ogs/s

        ! Calculate stage one acm, temperature and light limitation which
        ! are independent of stomatal conductance effects
        call acm_gpp_stage_1(mV)

        ! In all other cases iterate
        !potential_stomatal_conductance = zbrent('calculate_gs:find_gs_iWUE', &
        !                                        find_gs_iWUE,minimum_conductance,max_gs*leaf_canopy_light_scaling, &
        !                                        tol_gs*lai,iWUE_step*0.10d0)

!        if (do_iWUE) then
            ! Intrinsic WUE optimisation
            ! Check that the water restricted water range brackets the root solution for the bisection
            iWUE_upper = find_gs_iWUE(mV%potential_conductance, mV) !; iWUE_lower = find_gs_iWUE(min_gs)
            if ( iWUE_upper * find_gs_iWUE(min_gs, mV) > 0d0 ) then
                ! Then both proposals indicate that photosynthesis
                ! would be increased by greater opening of the stomata
                ! and is therefore water is limiting!
                mV%stomatal_conductance = mV%potential_conductance
                ! Exception being if both are positive - therefore assume
                ! lowest
                if (iWUE_upper > 0d0) mV%stomatal_conductance = mV%minimum_conductance
            else if (mV%potential_conductance < mV%minimum_conductance) then
                ! If the potential conductance is less than the hardcoded minimum
                ! assume stomatal conductance is the minimum and move on.
                mV%stomatal_conductance = mV%minimum_conductance
            else

                ! In all other cases iterate
                mV%stomatal_conductance = zbrent('calculate_gs:find_gs_iWUE', &
                                                 find_gs_iWUE,mV,mV%minimum_conductance,mV%potential_conductance, &
                                                 tol_gs*mV%lai,mV%iWUE_step*0.10d0)

            end if

    else

        ! if no LAI then there can be no stomatal conductance
        mV%potential_conductance = max_gs ; mV%minimum_conductance = vsmall
        mV%stomatal_conductance = vsmall

    endif ! if aerodynamic conductance > vsmall


  end subroutine calculate_stomatal_conductance
  !
  !------------------------------------------------------------------
  !
  subroutine meteorological_constants(input_temperature,input_temperature_K,input_vpd_kPa, mV)
    use carbon_model_memory

    ! Determine some multiple use constants used by a wide range of functions
    ! All variables here are linked to air temperature and thus invarient between
    ! iterations and can be stored in memory...

    implicit none

      type(model_working_variables) :: mV

    ! arguments
    double precision, intent(in) :: input_temperature, input_temperature_K, &
                                    input_vpd_kPa

    ! local variables
    double precision :: mult, &
              dynamic_viscosity    ! dynamic viscosity (kg.m-2.s-1)

    !
    ! Used for soil, canopy evaporation and transpiration
    !

    ! Density of air (kg/m3)
    mV%air_density_kg = 353d0/input_temperature_K
    ! Conversion ratio for m.s-1 -> mol.m-2.s-1
    mV%convert_ms1_mol_1 = const_sfc_pressure / (input_temperature_K*Rcon)
    ! latent heat of vapourisation,
    ! function of air temperature (J.kg-1)
    mV%lambda = 2501000d0-2364d0*input_temperature

    ! psychrometric constant (kPa K-1)
    mV%psych = (0.0646d0*exp(0.00097d0*input_temperature))
    ! Straight line approximation of the true slope; used in determining
    ! relationship slope
    mult = input_temperature+237.3d0
    ! 2502.935945 = 0.61078*17.269*237.3
    ! Rate of change of saturation vapour pressure with temperature (kPa.K-1)
    mV%slope = (2502.935945d0*exp(17.269d0*input_temperature/mult)) / (mult*mult)

    ! estimate frequently used atmmospheric demand component
    mV%ET_demand_coef = mV%air_density_kg*cpair*input_vpd_kPa

    !
    ! Used for soil evaporation and leaf level conductance
    !

    ! Determine diffusion coefficient (m2.s-1), temperature dependant (pressure dependence neglected). Jones p51; appendix 2
    ! Temperature adjusted from standard 20oC (293.15 K), NOTE that 1/293.15 = 0.003411223
    ! 0.0000242 = conversion to make diffusion specific for water vapor (um2.s-1)
    mV%water_vapour_diffusion = 0.0000242d0*((input_temperature_K/293.15d0)**1.75d0)

    !
    ! Used for calculation of leaf level conductance
    !

    ! Calculate the dynamic viscosity of air (kg.m-2.s-1)
    dynamic_viscosity = ((input_temperature_K**1.5d0)/(input_temperature_K+120d0))*1.4963d-6
    ! and kinematic viscosity (m2.s-1)
    mV%kinematic_viscosity = dynamic_viscosity/mV%air_density_kg

  end subroutine meteorological_constants
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_transpiration(transpiration, mV)
    use carbon_model_memory

    ! Models leaf cnaopy transpiration based on the Penman-Monteith model of
    ! evapotranspiration used to estimate SPA's daily evapotranspiration flux
    ! (kgH20.m-2.day-1).

    implicit none

      type(model_working_variables) :: mV

    ! arguments
    double precision, intent(out) :: transpiration ! kgH2O.m-2.day-1

    ! local variables
    double precision :: canopy_radiation & ! isothermal net radiation (W/m2)
                                  ,gs,gb   ! stomatal and boundary layer conductance (m.s-1)

    !!!!!!!!!!
    ! Estimate energy radiation balance (W.m-2)
    !!!!!!!!!!

    ! Absorbed shortwave radiation MJ.m-2.day-1 -> J.m-2.s-1
    canopy_radiation = mV%canopy_lwrad_Wm2 + (mV%canopy_swrad_MJday * 1d6 * mV%dayl_seconds_1)

    !!!!!!!!!!
    ! Calculate canopy conductance (to water vapour)
    !!!!!!!!!!

    ! Change units of potential stomatal conductance
    ! (mmolH2O.m-2.d-1 -> m.s-1).
    ! Note assumption of sea surface pressure only
    gs = mV%stomatal_conductance / mV%convert_ms1_mmol_1
    ! Scale aerodynamic conductance to canopy scale
    gb = mV%aerodynamic_conductance * mV%leaf_canopy_wind_scaling

    !!!!!!!!!!
    ! Calculate canopy evaporative fluxes (kgH2O/m2/day)
    !!!!!!!!!!

    ! Calculate numerator of Penman Montheith (kgH2O.m-2.day-1)
    ! NOTE: that restriction within water supply restriction is determined
    ! during stomatal conductance level.
    transpiration = ( ( (mV%slope*canopy_radiation) + (mV%ET_demand_coef*gb) ) &
                      / (mV%lambda*(mV%slope+(mV%psych*(1d0+gb/gs)))) )*mV%dayl_seconds

  end subroutine calculate_transpiration
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_wetcanopy_evaporation(wetcanopy_evap,storage, mV)
    use carbon_model_memory

    ! Estimates evaporation of canopy intercepted rainfall based on the Penman-Monteith model of
    ! evapotranspiration used to estimate SPA's daily evapotranspiration flux
    ! (kgH20.m-2.day-1).

    implicit none

      type(model_working_variables) :: mV

    ! arguments
    double precision, intent(inout) :: storage      ! canopy water storage kgH2O/m2
    double precision, intent(out) :: wetcanopy_evap ! kgH2O.m-2.day-1

    ! local variables
    double precision :: canopy_radiation, & ! isothermal net radiation (W/m2)
                                      gb    ! stomatal and boundary layer conductance (m.s-1)

    ! Assuming there is any rainfall, currently water on the canopy or dew formation
    if (mV%rainfall > 0d0 .or. storage > 0d0) then

        !!!!!!!!!!
        ! Calculate canopy conductance (to water vapour)
        !!!!!!!!!!

        ! Combine in series stomatal conductance with boundary layer
        gb = mV%aerodynamic_conductance * mV%leaf_canopy_wind_scaling

        !!!!!!!!!!
        ! Estimate energy radiation balance (W.m-2)
        !!!!!!!!!!

        ! Absorbed shortwave radiation MJ.m-2.day-1 -> J.m-2.s-1
        canopy_radiation = mV%canopy_lwrad_Wm2 + (mV%canopy_swrad_MJday * 1d6 * mV%dayl_seconds_1)

        !!!!!!!!!!
        ! Calculate canopy evaporative fluxes (kgH2O/m2/day)
        !!!!!!!!!!

        ! Calculate potential Penman Montheith (kgH2O.m-2.day-1)
        wetcanopy_evap = max(0d0,(((mV%slope*canopy_radiation) + (mV%ET_demand_coef*gb)) &
                                 / (mV%lambda*(mV%slope+mV%psych))) * mV%dayl_seconds)

        ! Update based on canopy water storage
        call canopy_interception_and_storage(wetcanopy_evap,storage, mV)

    else

        ! there is no water movement possible
        mV%intercepted_rainfall = 0d0 ; wetcanopy_evap = 0d0

    endif

  end subroutine calculate_wetcanopy_evaporation
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_potential_evaporation(potential_evap, mV)
    use carbon_model_memory

    ! Estimates potential surface evapotransporation based on the Penman-Monteith model
    ! (kgH20.m-2.day-1). FAO Chapter 3 Determination of ETo, see chapter 2 for derivation.

    implicit none

      type(model_working_variables) :: mV

    ! arguments
    double precision, intent(out) :: potential_evap ! kgH2O.m-2.day-1

    ! local variables
    double precision :: canopy_radiation  ! isothermal net radiation (W/m2)

    !!!!!!!!!!
    ! Estimate energy radiation balance (W.m-2)
    !!!!!!!!!!

    ! Absorbed shortwave radiation MJ.m-2.day-1
    canopy_radiation = mV%canopy_swrad_MJday + mV%soil_swrad_MJday &
                     + (mV%canopy_lwrad_Wm2 * 1d-6 * seconds_per_day) &
                     + (mV%soil_lwrad_Wm2 * 1d-6 * seconds_per_day)

    !!!!!!!!!!
    ! Calculate canopy evaporative fluxes (kgH2O/m2/day)
    !!!!!!!!!!

    ! Calculate numerator of Penman Montheith (kgH2O.m-2.day-1)
    ! NOTE: Rn - G, neglected as G (ground heat) near zero on daily scales
    ! 0.34 estimates the ratio of canopy and stomatal conductance
    ! 0.408 is the inverse of lambda as described in this code.
    potential_evap = ((0.408d0*mV%slope*canopy_radiation) + &
                      (mV%psych*(900d0 / (mV%meant + 273d0)) * mV%wind_spd * mV%vpd_kPa)) &
                   / (mV%slope + mV%psych * (1d0+0.34d0*mV%wind_spd))

  end subroutine calculate_potential_evaporation
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_soil_evaporation(soilevap, mV)
    use carbon_model_memory

    ! Estimate soil surface evaporation based on the Penman-Monteith model of
    ! evapotranspiration used to estimate SPA's daily evapotranspiration flux
    ! (kgH20.m-2.day-1).

    implicit none

      type(model_working_variables) :: mV

    ! arguments
    double precision, intent(out) :: soilevap ! kgH2O.m-2.day-1

    ! local variables
    double precision :: local_temp &
                   ,soil_radiation & ! isothermal net radiation (W/m2)
                            ,esurf & ! see code below
                             ,esat & ! soil air space saturation vapour pressure
                              ,gws   ! water vapour conductance through soil air space (m.s-1)

    ! oC -> K for local temperature value
    local_temp = mV%soilT + freeze

    !!!!!!!!!!
    ! Estimate energy radiation balance (W.m-2)
    !!!!!!!!!!

    ! Absorbed shortwave radiation MJ.m-2.day-1 -> J.m-2.s-1
    soil_radiation = mV%soil_lwrad_Wm2 + (mV%soil_swrad_MJday * 1d6 * mV%dayl_seconds_1)

    !!!!!!!!!!
    ! Calculate soil evaporative fluxes (kgH2O/m2/day)
    !!!!!!!!!!

    ! calculate saturated vapour pressure (kPa), function of temperature.
    esat = 0.1d0 * exp( 1.80956664d0 + ( 17.2693882d0 * local_temp - 4717.306081d0 ) / ( local_temp - 35.86d0 ) )
    mV%air_vapour_pressure = esat - mV%vpd_kPa

    ! Soil conductance to water vapour diffusion (m s-1)...
    gws = mV%porosity(1) * mV%water_vapour_diffusion / (tortuosity*mV%drythick)

    ! vapour pressure in soil airspace (kPa), dependent on soil water potential
    ! - Jones p.110. partial_molar_vol_water. Less vapour pressure of the air to
    ! estimate the deficit between soil and canopy air spaces
    esurf = (esat * exp( 1d6 * mV%SWP(1) * partial_molar_vol_water / (Rcon * local_temp) )) - mV%air_vapour_pressure

    ! Estimate potential soil evaporation flux (kgH2O.m-2.day-1)
    soilevap = ( ((mV%slope*soil_radiation) + (mV%air_density_kg*cpair*esurf*mV%soil_conductance)) &
               / (mV%lambda*(mV%slope+(mV%psych*(1d0+mV%soil_conductance/gws)))) ) * mV%dayl_seconds

    return

  end subroutine calculate_soil_evaporation
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_aerodynamic_conductance (mV)
    use carbon_model_memory

    !
    ! Calculates the aerodynamic or bulk canopy conductance (m.s-1). Here we
    ! assume neutral conditions due to the lack of an energy balance calculation
    ! in either ACM or DALEC. The equations used here are extracted from SPA.
    !

    implicit none

      type(model_working_variables) :: mV

    ! local variables
    double precision :: local_lai, &
           mixing_length_momentum, & ! mixing length parameter for momentum (m)
            length_scale_momentum, & ! length scale parameter for momentum (m)
                     canopy_decay    ! within canopy decay coefficient for wind

    ! parameters
    double precision, parameter :: foliage_drag = 0.2d0 ! foliage drag coefficient

    ! Restrict LAI used here to greater than a minium value which prevents un-realistic outputs
    local_lai = max(min_lai,mV%lai)

    ! Calculate the zero plane displacement and roughness length
    call z0_displacement(mV%ustar_Uh,local_lai, mV)
    ! Calculate friction velocity at tower height (reference height) (m.s-1)
    ! WARNING neutral conditions only; see WRF module_sf_sfclay.F for 'with
    ! stability versions'
    !    ustar = (wind_spd / log((tower_height-displacement)/roughl)) * vonkarman
    mV%ustar = mV%wind_spd * mV%ustar_Uh

    ! Both the length scale and mixing length are assumed to be constant within
    ! the canopy (under dense canopy conditions).
    ! Calculate length scale (lc) for momentum absorption within the canopy; Harman & Finnigan (2007)
    ! Calculate mixing length (lm) for vertical momentum within the canopy Harman & Finnigan (2008)
    length_scale_momentum = (4d0*canopy_height) / local_lai
    mixing_length_momentum = 2d0*(mV%ustar_Uh**3)*length_scale_momentum

    ! Estimate above canopy wind speed decay followin Harman & Finnigan (2008)
    ! and assuming neutral conditions only.
    call log_law_decay(mV)

    ! Calculate canopy decay coefficient with stability correction
    ! NOTE this is not consistent with canopy momentum decay done by Harman &
    ! Finnigan (2008) instead follows Nui & Yang 2004; Qin et al 2002.
    canopy_decay = sqrt((foliage_drag*canopy_height*local_lai)/mixing_length_momentum)

    ! Estimating the within canopy wind speed, we assume that the wind speed
    ! just inside of the canopy is most important.
    !canopy_wind = canopy_wind*exp((ustar_Uh*((canopy_height*1d0)-canopy_height))/mixing_length_momentum)

    ! Calculate_soil_conductance
    call calculate_soil_conductance(mixing_length_momentum,local_lai,canopy_decay, mV)
    ! Calculate top of canopy leaf conductance (m/s) for water vapour under forced convective conditions
    ! The leaf-to-canopy scaling is achieved with the leaf_canopy_wind_scaling scalar.
    call average_leaf_conductance(mV%aerodynamic_conductance, mV)

    ! Estimate leaf to canopy scaling factor for use with aerodynamic conductance.
    ! Based on the canopy scaling of photosynthetic capacity due to light from
    ! Sellers et al., (1992), Remote Sensing Environment, 42(3), 187-216.
    ! But now applied on the within canopy decay gradient.
    mV%leaf_canopy_wind_scaling = exp(canopy_decay*(1d0-(soil_roughl/canopy_height))) &
                             - exp(canopy_decay*(1d0-((mV%roughl+mV%displacement)/canopy_height)))
    mV%leaf_canopy_wind_scaling = mV%leaf_canopy_wind_scaling / exp(canopy_decay)

  end subroutine calculate_aerodynamic_conductance
  !
  !------------------------------------------------------------------
  !
  subroutine average_leaf_conductance(gv_forced, mV)
    use carbon_model_memory

    !
    ! Subroutine calculates the forced conductance of water vapour for non-cylinder within canopy leaves (i.e. broadleaf)
    ! Free convection (i.e. that driven by energy balance) is negelected here due to the lack of an energy balance
    ! calculation in DALEC. Should a energy balance be added then this code could be expanded include free conductance
    ! Follows a simplified approach to that used in SPA (Smallman et al 2013).
    !

    implicit none

      type(model_working_variables) :: mV

    ! arguments
    double precision, intent(out) :: gv_forced ! canopy conductance (m/s) for water vapour under forced convection

    ! local parameters
    double precision, parameter :: leaf_width_coef = 25d0, & ! (1/leaf_width) * 0.5,
                                                               ! where 0.5 accounts for one half
                                                               ! of the leaf used in water exchange
                                        leaf_width = 0.02d0    ! leaf width (m) (alternates 0.04, 0.08)
!                                                Pr = 0.72d0, & ! Prandtl number
!                                           Pr_coef = 1.05877d0 !1.18d0*(Pr**(0.33d0))
    ! local variables
    double precision :: &
              Sh_forced & ! Sherwood number under forced convection
                    ,Re   ! Reynolds number

    ! Sherwood number under forced convection. NOTE: 0.962 * Pr_coef = 1.018537
!    Sh_forced = 0.962d0*Pr_coef*(sqrt((leaf_width*canopy_wind)/kinematic_viscosity))
    Sh_forced = 1.018537d0*(sqrt((leaf_width*mV%canopy_wind)/mV%kinematic_viscosity))
    ! Estimate the the forced conductance of water vapour
    gv_forced = mV%water_vapour_diffusion*Sh_forced*leaf_width_coef

  end subroutine average_leaf_conductance
  !
  !------------------------------------------------------------------
  !
  subroutine log_law_decay (mV)
    use carbon_model_memory

    ! Standard log-law above canopy wind speed (m.s-1) decay under neutral
    ! conditions.
    ! See Harman & Finnigan 2008; Jones 1992 etc for details.

    implicit none

      type(model_working_variables) :: mV

    ! log law decay, NOTE: given canopy height (9 m) the log function reduces
    ! to a constant value down to ~ 7 decimal place (0.3161471806). Therefore
    ! 1/vonkarman * 0.31 = 0.7710906
    mV%canopy_wind = mV%ustar * vonkarman_1 * log((canopy_height-mV%displacement) / mV%roughl)

  end subroutine log_law_decay
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_field_capacity (mV)
    use carbon_model_memory

    use brent_zero, only: zbrent

    ! field capacity calculations for saxton eqns !

    implicit none

      type(model_working_variables) :: mV

    ! local variables..
    integer :: i
    double precision :: x1, x2

    x1 = 0.1d0 ; x2 = 0.7d0 ! low/high guess
    do i = 1 , nos_soil_layers+1
       mV%water_retention_pass = i
       ! field capacity is water content at which SWP = -10 kPa
       mV%field_capacity(i) = zbrent('water_retention:water_retention_saxton_eqns', &
                                   water_retention_saxton_eqns,mV, x1 , x2 , 0.001d0, 0d0 )
    enddo

  end subroutine calculate_field_capacity
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_daylength(doy,lat, mV)
    use carbon_model_memory

    ! Subroutine uses day of year and latitude (-90 / 90 degrees) as inputs,
    ! combined with trigonomic functions to calculate day length in hours and seconds

    implicit none

    ! arguments
    double precision, intent(in) :: doy, lat
    type(model_working_variables) :: mV
    
    ! local variables
    double precision :: dec, mult, sinld, cosld, aob

    !
    ! Estimate solar geometry variables needed
    !

    ! Declination
    ! NOTE: 0.002739726d0 = 1/365
    !    dec = - asin( sin( 23.45d0 * deg_to_rad ) * cos( 2d0 * pi * ( doy + 10d0 ) / 365d0 ) )
    !    dec = - asin( sin_dayl_deg_to_rad * cos( two_pi * ( doy + 10d0 ) / 365d0 ) )
    dec = - asin( sin_dayl_deg_to_rad * cos( two_pi * ( doy + 10d0 ) * 0.002739726d0 ) )

    ! latitude in radians
    mult = lat * deg_to_rad
    ! day length is estimated as the ratio of sin and cos of the product of declination an latitude in radiation
    sinld = sin( mult ) * sin( dec )
    cosld = cos( mult ) * cos( dec )
    aob = max(-1d0,min(1d0,sinld / cosld))

    ! estimate day length in hours and seconds and upload to module variables
    mV%dayl_hours = 12d0 * ( 1d0 + 2d0 * asin( aob ) * pi_1 )
    mV%dayl_seconds = mV%dayl_hours * seconds_per_hour

    ! return to user
    return

  end subroutine calculate_daylength
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_longwave_isothermal(canopy_temperature,soil_temperature, mV)
    use carbon_model_memory

    ! Subroutine estimates the isothermal net longwave radiation (W.m-2) for
    ! the canopy and soil surface. SPA uses a complex multi-layer radiative
    ! transfer scheme including reflectance, transmittance any absorption.
    ! However, for a given canopy vertical profiles, the LAI absorption
    ! relationship is readily predicted via Michaelis-Menten or
    ! non-rectangular hyperbola as done here.

    implicit none

      type(model_working_variables) :: mV

    ! arguments
    double precision, intent(in) :: canopy_temperature, soil_temperature ! oC

    ! local variables
    double precision :: lwrad, & ! downward long wave radiation from sky (W.m-2)
         transmitted_fraction, & ! fraction of LW which is not incident on the canopy
  canopy_transmitted_fraction, & !
                    delta_iso, &
        longwave_release_soil, & ! emission of long wave radiation from surfaces per m2
      longwave_release_canopy, & ! assuming isothermal condition (W.m-2)
            trans_lw_fraction, &
        reflected_lw_fraction, &
         absorbed_lw_fraction, &
      canopy_release_fraction, & ! fraction of longwave emitted from within the canopy to ultimately be released
   canopy_absorption_from_sky, & ! canopy absorbed radiation from downward LW (W.m-2)
  canopy_absorption_from_soil, & ! canopy absorbed radiation from soil surface (W.m-2)
                  canopy_loss, & ! longwave radiation released from canopy surface (W.m-2).
                                 ! i.e. this value is released from the top and
                                 ! the bottom
       soil_incident_from_sky, &
     soil_absorption_from_sky, & ! soil absorbed radiation from sky (W.m-2)
  soil_absorption_from_canopy    ! soil absorbed radiation emitted from canopy (W.m-2)

    ! local parameters
    double precision, parameter :: nos_layers = 4d0 & ! Number of canopy layers in source model
                                  ,nos_layers_1 = nos_layers ** (-1d0) &
                                  ,clump = 0.75d0   & ! Clumping factor (1 = uniform, 0 totally clumped, mean = 0.75)
                                                      ! He et al., (2012) http://dx.doi.org/10.1016/j.rse.2011.12.008
                                  ,decay = -0.5d0     ! decay coefficient for incident radiation

    ! estimate long wave radiation from atmosphere (W.m-2)
    lwrad = emiss_boltz * (mV%maxt+freeze-20d0) ** 4
    ! estimate isothermal long wave emission per unit area
    longwave_release_soil = emiss_boltz * (soil_temperature+freeze) ** 4
    ! estimate isothermal long wave emission per unit area
    longwave_release_canopy = emiss_boltz * (canopy_temperature+freeze) ** 4

    !!!!!!!!!!
    ! Determine fraction of longwave absorbed by canopy and returned to the sky
    !!!!!!!!!!

    ! First, we consider how much radiation is likely to be incident on the
    ! canopy, or put another way what fraction passes straight through the
    ! canopy?
    transmitted_fraction = exp(decay * mV%lai * clump)

    ! second, we partition the radiation which is incident on the canopy into
    ! that which is transmitted, reflected or absorbed.

    ! Likewise we assume that the reflectance and transmittance are equal
    ! However, the non-linear interception under Beer's Law means that actual
    ! canopy transmittenace to the soil surface and reflectance back to the sky
    ! skews towards reduced transmittance at higher LAI. Both transmittance and
    ! reflectance follow linear a relationship with a common intercept
    ! NOTE: 0.02 = (1-emissivity) * 0.5.
    ! NOTE: lai*0.5 reflects that interacting LAI will be somewhere within the canopy
    !       and this its transmittance or reflectance will not be subject to the entire canopy
    canopy_transmitted_fraction = exp(decay * mV%lai * 0.5d0 * clump)
    trans_lw_fraction     = 0.02d0 * canopy_transmitted_fraction
    reflected_lw_fraction = 0.02d0 * canopy_transmitted_fraction
    ! Absorption is the residual
    absorbed_lw_fraction = 1d0 - trans_lw_fraction - reflected_lw_fraction

    ! Calculate the potential absorption of longwave radiation lost from the
    ! canopy to soil / sky. The first line below estimates the actual release fraction,
    ! second line accounts for changing LAI scaling within a vertical canopy.
    ! NOTE: That assuming the Beer's law emission from a single canopy layer leads to a rough 50 % underestimate of LW emission.
    !       This is why there is the nos_layer correction here
    canopy_release_fraction = (1d0 - (max_lai_lwrad_release*mV%lai) / (mV%lai+lai_half_lwrad_release)) &
                            * (1d0 - exp(decay * mV%lai * nos_layers_1 * clump)) * nos_layers

    !!!!!!!!!!
    ! Distribute longwave from sky
    !!!!!!!!!!

    ! Estimate the radiation which directly bypasses the canopy...
    soil_incident_from_sky = lwrad * transmitted_fraction
    ! ...and update the canopy intercepted radiation
    lwrad = lwrad - soil_incident_from_sky

    ! long wave absorbed by the canopy from the sky
    canopy_absorption_from_sky = lwrad * absorbed_lw_fraction
    ! Long wave absorbed by soil from the sky, soil absorption assumed to be
    ! equal to emissivity
    soil_incident_from_sky = soil_incident_from_sky + (trans_lw_fraction * lwrad)
    soil_absorption_from_sky = soil_incident_from_sky * emissivity
    ! Long wave reflected directly back into sky
    mV%sky_lwrad_Wm2 = lwrad * reflected_lw_fraction

    !!!!!!!!!!
    ! Distribute longwave from soil
    !!!!!!!!!!

    ! Calculate longwave radiation coming up from the soil plus the radiation
    ! which is reflected
    canopy_absorption_from_soil = longwave_release_soil + (soil_incident_from_sky * (1d0-emissivity))
    ! First how much directly bypasses the canopy...
    mV%sky_lwrad_Wm2 = mV%sky_lwrad_Wm2 + (canopy_absorption_from_soil * transmitted_fraction)
    canopy_absorption_from_soil = canopy_absorption_from_soil * (1d0 - transmitted_fraction)
    ! Second, use this to estimate the longwave returning to the sky
    mV%sky_lwrad_Wm2 = mV%sky_lwrad_Wm2 + (canopy_absorption_from_soil * trans_lw_fraction)
    ! Third, now calculate the longwave from the soil surface absorbed by the
    ! canopy
    canopy_absorption_from_soil = canopy_absorption_from_soil * absorbed_lw_fraction

    !!!!!!!!!!
    ! Distribute longwave originating from the canopy itself
    !!!!!!!!!!

    ! calculate two-sided long wave radiation emitted from canopy which is
    ! ultimately lost from to soil or sky (i.e. this value is used twice, once
    ! to soil once to sky)
    canopy_loss = longwave_release_canopy * canopy_release_fraction
    ! Calculate longwave absorbed by soil which is released by the canopy itself
    soil_absorption_from_canopy = canopy_loss * emissivity
    ! Canopy released longwave returned to the sky
    mV%sky_lwrad_Wm2 = mV%sky_lwrad_Wm2 + canopy_loss

    !!!!!!!!!!
    ! Isothermal net long wave canopy and soil balance (W.m-2)
    !!!!!!!!!!

    ! determine isothermal net canopy. Note two canopy_loss used to account for
    ! upwards and downwards emissions
    mV%canopy_lwrad_Wm2 = (canopy_absorption_from_sky + canopy_absorption_from_soil) - (canopy_loss + canopy_loss)
    ! determine isothermal net soil
    mV%soil_lwrad_Wm2 = (soil_absorption_from_sky + soil_absorption_from_canopy) - longwave_release_soil

    !!!!!!!!!!
    ! Convert isothermal to net radiation
    !!!!!!!!!!

    ! Apply linear correction to soil surface isothermal->net longwave radiation
    ! balance based on absorbed shortwave radiation
    delta_iso = (soil_iso_to_net_coef_LAI * mV%lai) + &
                (soil_iso_to_net_coef_SW * (mV%soil_swrad_MJday * 1d6 * seconds_per_day_1)) + &
                 soil_iso_to_net_const
    ! In addition to the iso to net adjustment, SPA analysis shows that soil net never gets much below zero
    mV%soil_lwrad_Wm2 = mV%soil_lwrad_Wm2 + delta_iso
    ! Estimate the mean soil surface temperature as a result of net radiation update
    mV%soilT = (((longwave_release_soil - delta_iso) / emiss_boltz) ** (0.25d0)) - freeze

    ! Apply linear correction to canopy isothermal->net longwave radiation
    ! balance based on absorbed shortwave radiation
    delta_iso = (canopy_iso_to_net_coef_LAI * mV%lai) + &
                (canopy_iso_to_net_coef_SW * (mV%canopy_swrad_MJday * 1d6 * seconds_per_day_1)) + &
                canopy_iso_to_net_const
    mV%canopy_lwrad_Wm2 = mV%canopy_lwrad_Wm2 + delta_iso
    ! Estimate the mean leaf temperature as a result of net radiation update.
    ! This can only be attempted if the canopy release fraction (an empirical fit) is greater than
    ! a minimum value. Otherwise the leafT defaults to infinity which is unrealistic.
    if (canopy_release_fraction > 1d-4) then
        mV%leafT = (((((canopy_loss + canopy_loss) - delta_iso) / (canopy_release_fraction * 2d0)) &
                 / emiss_boltz) ** (0.25d0)) - freeze
    end if

  end subroutine calculate_longwave_isothermal
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_radiation_balance (mV)
    use carbon_model_memory

    implicit none

      type(model_working_variables) :: mV

    ! subroutine call ensures that both shortwave and longwave radiation balance
    ! are calculated at the same time but with the more readable code split
    ! between a shortwave and longwave specific subroutines.

    ! NOTE: that this code provides a daily timescale linear correction on
    ! isothermal longwave balance to net based on soil surface incident shortwave
    ! radiation

    ! Estimate shortwave radiation balance
    call calculate_shortwave_balance(mV)
    ! Estimate isothermal long wave radiation balance
    !call calculate_longwave_isothermal(meant,meant)
    call calculate_longwave_isothermal(mV%leafT,mV%soilT, mV)

  end subroutine calculate_radiation_balance
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_shortwave_balance (mV)
    use carbon_model_memory

    ! Subroutine estimates the canopy and soil absorbed shortwave radiation
    ! (MJ/m2/day).
    ! Radiation absorption is paritioned into NIR and PAR for canopy, and NIR +
    ! PAR for soil.

    ! SPA uses a complex multi-layer radiative transfer scheme including
    ! reflectance, transmittance any absorption. However, for a given
    ! canopy vertical profiles, the LAI absorption relationship is readily
    ! predicted via Michaelis-Menten or non-rectangular hyperbola as done here.

    implicit none

      type(model_working_variables) :: mV

    ! local variables
    double precision :: balance                     &
                       ,transmitted_fraction        &
                       ,canopy_transmitted_fraction &
                       ,absorbed_nir_fraction_soil  &
                       ,absorbed_par_fraction_soil  &
                       ,fsnow,par,nir               &
                       ,soil_nir_MJday              &
                       ,trans_nir_MJday             &
                       ,trans_par_MJday             &
                       ,canopy_nir_MJday            &
                       ,refl_par_MJday              &
                       ,refl_nir_MJday              &
                       ,reflected_nir_fraction      & !
                       ,reflected_par_fraction      & !
                       ,absorbed_nir_fraction       & !
                       ,absorbed_par_fraction       & !
                       ,trans_nir_fraction          & !
                       ,trans_par_fraction

    ! local parameters
    double precision, parameter :: clump = 0.75d0 & ! Clumping factor (1 = uniform, 0 totally clumped, mean = 0.75)
                                                    ! He et al., (2012) http://dx.doi.org/10.1016/j.rse.2011.12.008
                                  ,decay = -0.5d0 & ! decay coefficient for incident radiation
                                  ,newsnow_nir_abs = 0.27d0 & ! NIR absorption fraction
                                  ,newsnow_par_abs = 0.05d0   ! PAR absorption fraction

    !!!!!!!!!!
    ! Determine canopy absorption, reflectance and transmittance as function of
    ! LAI
    !!!!!!!!!!

    ! First, we consider how much radiation is likely to be incident on the
    ! canopy, or put another way what fraction passes straight through the
    ! canopy?
    transmitted_fraction = exp(decay * mV%lai * clump)
    ! Estimate the integral of light interception for use as a leaf to canopy
    ! scaler for photosynthesis, transpiration, and gs
    ! Based on the canopy scaling of photosynthetic capacity due to light from
    ! Sellers et al., (1992), Remote Sensing Environment, 42(3), 187-216.
    mV%leaf_canopy_light_scaling = (1d0-transmitted_fraction) / (-decay*clump)

    ! Second, of the radiation which is incident on the canopy what fractions
    ! are transmitted through, reflected from or absorbed by the canopy

    canopy_transmitted_fraction = exp(decay * mV%lai * 0.5d0 * clump)

    ! Canopy transmitted of PAR & NIR radiation towards the soil
    trans_par_fraction = canopy_transmitted_fraction * max_par_transmitted
    trans_nir_fraction = canopy_transmitted_fraction * max_nir_transmitted
    ! Canopy reflected of near infrared and photosynthetically active radiation
    reflected_nir_fraction = canopy_transmitted_fraction * max_nir_reflected
    reflected_par_fraction = canopy_transmitted_fraction * max_par_reflected
    ! Canopy absorption of near infrared and photosynthetically active radiation
    absorbed_nir_fraction = 1d0 - reflected_nir_fraction - trans_nir_fraction
    absorbed_par_fraction = 1d0 - reflected_par_fraction - trans_par_fraction

    !!!!!!!!!!
    ! Estimate canopy absorption of incoming shortwave radiation
    !!!!!!!!!!

    ! Estimate multiple use par and nir components
    par = sw_par_fraction * mV%swrad
    nir = (1d0 - sw_par_fraction) * mV%swrad

    ! Estimate the radiation which directly bypasses the canopy...
    trans_par_MJday = par * transmitted_fraction
    trans_nir_MJday = nir * transmitted_fraction
    ! ...and update the canopy intercepted radiation
    par = par - trans_par_MJday
    nir = nir - trans_nir_MJday

! NOTE: Can something about snow reflectance be added here to reduce radiation available for subsequent absorption / reflectance, and just reflect back to sky?

    ! Estimate incoming shortwave radiation absorbed, transmitted and reflected
    ! by the canopy (MJ.m-2.day-1)
    mV%canopy_par_MJday = par * absorbed_par_fraction
    canopy_nir_MJday = nir * absorbed_nir_fraction
    trans_par_MJday = trans_par_MJday + (par * trans_par_fraction)
    trans_nir_MJday = trans_nir_MJday + (nir * trans_nir_fraction)
    refl_par_MJday = par * reflected_par_fraction
    refl_nir_MJday = nir * reflected_nir_fraction

    !!!!!!!!!
    ! Estimate soil absorption of shortwave passing through the canopy
    !!!!!!!!!

    ! Update soil reflectance based on snow cover
    if (mV%snow_storage > 0d0) then
        fsnow = 1d0 - exp( - mV%snow_storage * 0.1d0 )  ! fraction of snow cover on the ground
        absorbed_par_fraction_soil = ((1d0 - fsnow) * soil_swrad_absorption) + (fsnow * newsnow_par_abs)
        absorbed_nir_fraction_soil = ((1d0 - fsnow) * soil_swrad_absorption) + (fsnow * newsnow_nir_abs)
    else
        absorbed_par_fraction_soil = soil_swrad_absorption
        absorbed_nir_fraction_soil = soil_swrad_absorption
    endif

    ! Then the radiation incident and ultimately absorbed by the soil surface
    ! itself (MJ.m-2.day-1)
    mV%soil_par_MJday = trans_par_MJday * absorbed_par_fraction_soil
    soil_nir_MJday = trans_nir_MJday * absorbed_nir_fraction_soil
    ! combine totals for use is soil evaporation
    mV%soil_swrad_MJday = soil_nir_MJday + mV%soil_par_MJday

    !!!!!!!!!
    ! Estimate canopy absorption of soil reflected shortwave radiation
    ! This additional reflection / absorption cycle is needed to ensure > 0.99
    ! of incoming radiation is explicitly accounted for in the energy balance.
    !!!!!!!!!

    ! calculate multiple use variables
    par = trans_par_MJday-mV%soil_par_MJday
    nir = trans_nir_MJday-soil_nir_MJday
    ! how much of the reflected radiation directly bypasses the canopy...
    refl_par_MJday = refl_par_MJday + (par * transmitted_fraction)
    refl_nir_MJday = refl_nir_MJday + (nir * transmitted_fraction)
    ! ...and update the canopy on this basis
    par = par * (1d0-transmitted_fraction)
    nir = nir * (1d0-transmitted_fraction)

    ! Update the canopy radiation absorption based on the reflected radiation
    ! (MJ.m-2.day-1)
    mV%canopy_par_MJday = mV%canopy_par_MJday + (par * absorbed_par_fraction)
    canopy_nir_MJday = canopy_nir_MJday + (nir * absorbed_nir_fraction)
    ! Update the total radiation reflected back into the sky, i.e. that which is
    ! now transmitted through the canopy
    refl_par_MJday = refl_par_MJday + (par * trans_par_fraction)
    refl_nir_MJday = refl_nir_MJday + (nir * trans_nir_fraction)

    ! Combine to estimate total shortwave canopy absorbed radiation
    mV%canopy_swrad_MJday = mV%canopy_par_MJday + canopy_nir_MJday

!    ! check energy balance
!    balance = swrad - canopy_par_MJday - canopy_nir_MJday - refl_par_MJday - refl_nir_MJday - soil_swrad_MJday
!    if (((balance - swrad) / swrad) > 0.01) then
!        print*,"SW residual frac = ",(balance - swrad) / swrad,"SW residual = ",balance,"SW in = ",swrad
!    endif

  end subroutine calculate_shortwave_balance
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_Rtot (mV)
    use carbon_model_memory

    ! Purpose of this subroutine is to calculate the minimum soil-root hydraulic
    ! resistance input into ACM. The approach used here is identical to that
    ! found in SPA.

    ! local variables
    integer :: i, rooted_layer

      type(model_working_variables) :: mV
    double precision :: transpiration_resistance,root_reach_local, &
                        slpa, mult, prev, exp_func!, root_depth_50, bonus
    double precision, dimension(nos_root_layers) :: Rcond_layer, &
                                                    root_mass,  &
                                                    root_length
    double precision, parameter :: rootdist_tol = 13.81551d0 ! log(1d0/rootdist_tol - 1d0) were rootdist_tol = 1d-6
                                  !rootdist_tol = 1d-6!, & ! Root density assessed for the max rooting depth
                                   !root_depth_frac_50 = 0.25d0 ! fractional soil depth above which 50 %
                                                               ! of the root mass is assumed to be located

    ! reset water flux
    mV%total_water_flux = 0d0 ; mV%water_flux_mmolH2Om2s = 0d0 ; mV%wSWP = 0d0 ; mV%rSWP = 0d0 ; mV%Reff = 0d0
    slpa = 0d0 ; root_length = 0d0 ; root_mass = 0d0 ; Rcond_layer = 0d0 ; mV%conductance_mmolH2OMPam2s = 0d0
    ! calculate soil depth to which roots reach
    mV%root_reach = mV%max_depth * mV%root_biomass / (mV%root_k + mV%root_biomass)
    ! calculate the plant hydraulic resistance component. Currently unclear
    ! whether this actually varies with height or whether tall trees have a
    ! xylem architecture which keeps the whole plant conductance (gplant) 1-10 (ish).
    !    transpiration_resistance = (gplant * lai)**(-1d0)
    ! Following Weiburg function, the potential conductance is reduced under increasing
    ! potential differences between the canopy and soil
 !SOME THOUGHT NEEDED HERE ON HOW TO QUANTIFY THIS IMPACT ON GPLANT, AS STRESS VARIES DEPENDING ON THE ROOT ZONE LAYERS...
    !transpiration_resistance = gplant * exp(-(abs(SWP(1:nos_root_layers) - minlwp) / s1)**(s2))
    !transpiration_resistance = canopy_height / (transpiration_resistance * max(min_lai,lai))
    transpiration_resistance = canopy_height / (gplant * max(min_lai,mV%lai))

    !!!!!!!!!!!
    ! calculate current steps soil hydraulic conductivity
    !!!!!!!!!!!

    ! seperately calculate the soil conductivity as this applies to each layer
    do i = 1, nos_soil_layers
       call calculate_soil_conductivity(i,mV%soil_waterfrac(i),mV%soil_conductivity(i), mV)
    end do ! soil layers

    !!!!!!!!!!!
    ! Calculate root profile
    !!!!!!!!!!!

    ! The original SPA src generates an exponential distribution which aims
    ! to maintain 50 % of root biomass in the top 25 % of the rooting depth.
    ! In a simple 3 root layer system this can be estimates more simply

!    ! top 25 % of root profile
!    root_depth_50 = root_reach * root_depth_frac_50
!    if (root_depth_50 <= layer_thickness(1)) then
!
!        ! Greater than 50 % of the fine root biomass can be found in the top
!        ! soil layer
!
!        ! Start by assigning all 50 % of root biomass to the top soil layer
!        root_mass(1) = fine_root_biomass * 0.5d0
!        ! Then quantify how much additional root is found in the top soil layer
!        ! assuming that the top 25 % depth is found somewhere within the top
!        ! layer
!        bonus = (fine_root_biomass-root_mass(1)) &
!              * (layer_thickness(1)-root_depth_50) / (root_reach - root_depth_50)
!        root_mass(1) = root_mass(1) + bonus
!        ! partition the remaining root biomass between the seconds and third
!        ! soil layers
!        if (root_reach > sum(layer_thickness(1:2))) then
!            root_mass(2) = (fine_root_biomass - root_mass(1)) &
!                         * (layer_thickness(2)/(root_reach-layer_thickness(1)))
!            root_mass(3) = fine_root_biomass - sum(root_mass(1:2))
!        else
!            root_mass(2) = fine_root_biomass - root_mass(1)
!        endif
!
!    else if (root_depth_50 > layer_thickness(1) .and. root_depth_50 <= sum(layer_thickness(1:2))) then
!
!        ! Greater than 50 % of fine root biomass found in the top two soil
!        ! layers. We will divide the root biomass uniformly based on volume,
!        ! plus bonus for the second layer (as done above)
!        root_mass(1) = fine_root_biomass * (layer_thickness(1)/root_depth_50)
!        root_mass(2) = fine_root_biomass * ((root_depth_50-layer_thickness(1))/root_depth_50)
!        root_mass(1:2) = root_mass(1:2) * 0.5d0
!
!        ! determine bonus for the seconds layer
!        bonus = (fine_root_biomass-sum(root_mass(1:2))) &
!              * ((sum(layer_thickness(1:2))-root_depth_50)/(root_reach-root_depth_50))
!        root_mass(2) = root_mass(2) + bonus
!        root_mass(3) = fine_root_biomass - sum(root_mass(1:2))
!
!    else
!
!        ! Greater than 50 % of fine root biomass stock spans across all three
!        ! layers
!        root_mass(1:2) = fine_root_biomass * 0.5d0 * (layer_thickness(1:2)/root_depth_50)
!        root_mass(3) = fine_root_biomass - sum(root_mass(1:2))
!
!    endif
!    ! now convert root mass into lengths
!    root_length = root_mass * root_mass_length_coef_1
!!    root_length = root_mass / (root_density * root_cross_sec_area)

    !!!!!!!!!!!
    ! Calculate hydraulic properties and each rooted layer
    !!!!!!!!!!!

    ! calculate and accumulate steady state water flux in mmol.m-2.s-1
    ! NOTE: Depth correction already accounted for in soil resistance
    ! calculations and this is the maximum potential rate of transpiration
    ! assuming saturated soil and leaves at their minimum water potential.
    mV%demand = max(0d0, (mV%SWP(1:nos_root_layers) - (head*canopy_height)) - mV%minlwp )
    ! now loop through soil layers, where root is present
    rooted_layer = 1
    ! Determine the exponential coefficient needed for an exponential decay to the current root reach
    ! Exponential decay profile following:
    ! Y = 1 / (1 + exp(-B * Z)), where Y = density at Z, B = gradient, Z = depth
    ! To determine gradient for current maximum root depth assuming density reaches rootdist_tol value, rearranges to:
    ! B = ln(1/Y - 1) / Z
!d = seq(0,2, 0.01) ; c = -2.6 ; rmax = 1 ; d50 = 0.25 ; rd = rmax / (1+ (d/d50)**c)
!rmax = rd * (1 + (d/d50)**c)
!(((rmax / rd) - 1)**(1/c)) * d50 = d ! Depth at which rd = 99 %
    !slpa = log(1d0/rootdist_tol - 1d0) / root_reach
    slpa = rootdist_tol / mV%root_reach
    prev = 1d0
    do i = 1, nos_root_layers
       ! Determine the exponential function for the current cumulative depth
       exp_func = exp(-slpa * sum(mV%layer_thickness(1:i)))
       ! Calculate the difference in the integral between depths, i.e. the proportion of root in the current volume
       mult = prev - (1d0 - (1d0/(1d0+exp_func)) + (0.5d0 * exp_func))
       ! Assign fine roo the the current layer...
       root_mass(i) = mV%fine_root_biomass * mult
       ! and determine the associated amount of root
       root_length(i) = root_mass(i) * root_mass_length_coef_1
       prev = prev - mult
       ! If there is root in the current layer then we should calculate the resistances
       if (root_mass(i) > 0d0) then
           ! Track the deepest root layer assessed
           rooted_layer = i
           ! if there is root then there is a water flux potential...
           root_reach_local = min(mV%root_reach,mV%layer_thickness(i))
           ! calculate and accumulate steady state water flux in mmol.m-2.s-1
           call plant_soil_flow(i,root_length(i),root_mass(i) &
                               ,mV%demand(i),root_reach_local &
                               ,transpiration_resistance,Rcond_layer(i), mV)
       else
           ! ...if there is not then we wont have any below...
           exit
       end if ! root present in current layer?
    end do ! nos_root_layers
    ! Turn the output resistance into conductance
    Rcond_layer = Rcond_layer**(-1d0)

    ! If freezing then assume soil surface is frozen, therefore no water flux
    if (mV%soilT < 1d0) then
        mV%water_flux_mmolH2Om2s(1) = 0d0
        Rcond_layer(1) = 0d0
    end if

    ! Calculate sum value (mmolH2O.m-2.s-1)
    mV%total_water_flux = sum(mV%water_flux_mmolH2Om2s)
    ! Calculate effective resistance
    ! NOTE: minimum condition used to guard against zero conductance and propagation of Inf / NaN
    ! through the model structure/
    mV%Reff = min(1d6,sum(mV%conductance_mmolH2OMPam2s)**(-1d0))
    if (mV%total_water_flux <= vsmall) then
        ! Set values for no water flow situation
        mV%uptake_fraction = (mV%layer_thickness(1:nos_root_layers) / sum(mV%layer_thickness(1:nos_root_layers)))
        ! Estimate weighted soil water potential based on fractional extraction from soil layers
        mV%wSWP = sum(mV%SWP(1:nos_root_layers) * mV%uptake_fraction(1:nos_root_layers))
        ! rSWP based on the conductance due to the roots themselves.
        ! However, similar to the wSWP we need a special case calculation
        ! when there is no extraction from the soil. Here we use the ratio of root mass itself.
        mV%rSWP = sum(mV%SWP(1:rooted_layer) * (root_mass(1:rooted_layer) / sum(root_mass(1:rooted_layer))))
        mV%total_water_flux = 0d0
      else
        ! calculate weighted SWP and uptake fraction
        mV%uptake_fraction(1:nos_root_layers) = mV%water_flux_mmolH2Om2s(1:nos_root_layers) / mV%total_water_flux
        ! Estimate weighted soil water potential based on fractional extraction from soil layers
        mV%wSWP = sum(mV%SWP(1:nos_root_layers) * mV%uptake_fraction(1:nos_root_layers))
        ! rSWP based on the conductance due to the roots themselves.
        ! The idea being that the plant may hedge against growth based on the majority of the
        ! profile being dry while not losing leaves within some toleration.
        mV%rSWP = sum(mV%SWP(1:rooted_layer) * (Rcond_layer(1:rooted_layer) / sum(Rcond_layer(1:rooted_layer))))
    endif

    ! and return
    return

  end subroutine calculate_Rtot
  !
  !-----------------------------------------------------------------
  !
  subroutine canopy_interception_and_storage(potential_evaporation,storage, mV)
    use carbon_model_memory

    ! Simple daily time step integration of canopy rainfall interception, runoff
    ! and rainfall (kgH2O.m-2.s-1). NOTE: it is possible for intercepted rainfall to be
    ! negative if stored water running off into the soil is greater than
    ! rainfall (i.e. when leaves have died between steps)

    implicit none

      type(model_working_variables) :: mV

    ! arguments
    double precision, intent(inout) :: storage, & ! canopy water storage (kgH2O/m2)
                         potential_evaporation    ! wet canopy evaporation (kgH2O.m-2.day-1),
                                                  ! enters as potential but leaves as water balance adjusted.
                                                  ! Note that this assumes a completely wet leaf surface
    ! local variables
    double precision :: a, through_fall, max_storage, max_storage_1, daily_addition, wetcanopy_evaporation &
                       ,potential_drainage_rate ,drain_rate, evap_rate, initial_canopy, co_mass_balance, dx, dz, tmp(3)
    ! local parameters
    double precision, parameter :: CanIntFrac = -0.5d0,     & ! Coefficient scaling rainfall interception fraction with LAI
                                        clump = 0.75d0,     & ! Clumping factor (1 = uniform, 0 totally clumped)
                                                              ! He et al., (2012) http://dx.doi.org/10.1016/j.rse.2011.12.008
                                  CanStorFrac = 0.2d0,      & ! Coefficient scaling canopy water storage with LAI
                                 RefDrainRate = 0.002d0,    & ! Reference drainage rate (mm/min; Rutter et al 1975)
                                  RefDrainLAI = 0.952381d0, & ! Reference drainage 1/LAI (m2/m2; Rutter et al 1975, 1/1.05)
                                 RefDrainCoef = 3.7d0,      & ! Reference drainage Coefficient (Rutter et al 1975)
                               RefDrainCoef_1 = RefDrainCoef ** (-1d0)

    ! hold initial canopy storage in memory
    initial_canopy = storage
    ! determine maximum canopy storage & through fall fraction
    through_fall = exp(CanIntFrac*mV%lai*clump)
    ! maximum canopy storage (mm); minimum is applied to prevent errors in
    ! drainage calculation. Assume minimum capacity due to wood.
    max_storage = max(min_storage,CanStorFrac*mV%lai)
    ! caclulate inverse for efficient calculations below
    max_storage_1 = max_storage**(-1d0)
    ! potential intercepted rainfall (kgH2O.m-2.s-1)
    mV%intercepted_rainfall = mV%rainfall * (1d0 - through_fall)

    ! calculate drainage coefficients (Rutter et al 1975); Corsican Pine
    ! 0.002 is canopy specific coefficient modified by 0.002*(max_storage/1.05)
    ! where max_storage is the canopy maximum capacity (mm) (LAI based) and
    ! 1.05 is the original canopy capacitance
    a = log( RefDrainRate * ( max_storage * RefDrainLAI ) ) - RefDrainCoef * max_storage

    ! average rainfall intercepted by canopy (kgH2O.m-2.day-1)
    daily_addition = mV%intercepted_rainfall * seconds_per_day

    ! reset cumulative variables
    through_fall = 0d0 ; wetcanopy_evaporation = 0d0
    drain_rate = 0d0 ; evap_rate = 0d0

    ! add rain to the canopy and overflow as needed
    storage = storage + daily_addition

    if (storage > max_storage) then

        if (potential_evaporation > 0d0) then

            ! Assume water drainage will always occur an order magnitude above evaporation
            ! so water above canopy capacity is drained. Water below max_storage is accessable by evaporation only.

            ! Assume drainage is all water above the maximum canopy storage (kg/m2/day)
            drain_rate = storage - max_storage

            ! Estimate evaporation from remaining water (i.e. that left after
            ! initial co-access of evaporation and drainage).
            ! Assume evaporation is now restricted by:
            ! 1) energy already spent on evaporation (the -evap_rate) and
            ! 2) linear increase in surface resistance as the leaf surface
            ! dries (i.e. the 0.5).
            evap_rate = min(potential_evaporation * 0.5d0 * storage * max_storage_1, storage - drain_rate)

        else

            ! Load dew formation to the current local evap_rate variable
            evap_rate = potential_evaporation
            ! Restrict drainage the quantity above max_storage, adding dew formation too
            drain_rate = (storage - evap_rate) - max_storage

        endif

    else

        ! no drainage just apply evaporation / dew formation fluxes directly
        drain_rate = 0d0 ; evap_rate = potential_evaporation
        if (evap_rate > 0d0) then
            ! evaporation restricted by fraction of surface actually covered
            ! in water and integrated over period to bare leaf (i.e. the *0.5)
            evap_rate = evap_rate * storage * max_storage_1 * 0.5d0
            ! and the total amount of water
            evap_rate = min(evap_rate,storage)
        else
            ! then dew formation has occurred, if this pushes storage > max_storage add it to drainage
            drain_rate = max(0d0,(storage - evap_rate) - max_storage)
        endif ! evap_rate > 0

    endif ! storage > max_storage

    ! update canopy storage with water flux
    storage = storage - evap_rate - drain_rate
    wetcanopy_evaporation = wetcanopy_evaporation + evap_rate
    through_fall = through_fall + drain_rate

    ! correct intercepted rainfall rate to kgH2O.m-2.s-1
    mV%intercepted_rainfall = mV%intercepted_rainfall - (through_fall * seconds_per_day_1)

!    ! sanity checks; note 1e-8 prevents precision errors causing flags
!    if (intercepted_rainfall > rainfall .or. storage < -1d-8 .or. &
!       (wetcanopy_evaporation * days_per_step_1) > (1d-8 + initial_canopy + (rainfall*seconds_per_day)) ) then
!        print*,"Condition 1",intercepted_rainfall > rainfall
!        print*,"Condition 2",storage < -1d-8
!        print*,"Condition 3",(wetcanopy_evaporation * days_per_step_1) > (1d-8 + initial_canopy + (rainfall*seconds_per_day))
!        print*,"storage (kgH2O/m2)",storage,"max_storage (kgH2O/m2)",max_storage,"initial storage (kgH2O/m2)", initial_canopy
!        print*,"rainfall (kgH2O/m2/day)", rainfall*seconds_per_day, "through_fall (kgH2O/m2/day)", (through_fall * days_per_step_1)
!        print*,"through_fall_total (kgH2O/m2/step)",through_fall
!        print*,"potential_evaporation (kgH2O/m2/day)",potential_evaporation
!        print*,"actual evaporation    (kgH2O/m2/day)",wetcanopy_evaporation * days_per_step_1
!        stop
!    endif

    ! average evaporative flux to daily rate (kgH2O/m2/day)
    potential_evaporation = wetcanopy_evaporation

    ! final clearance of canopy storage of version small values at the level of system precision
    if (storage < 10d0*vsmall) storage = 0d0

  end subroutine canopy_interception_and_storage
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_update_soil_water(Eleaf,Esoil,Esnow,rainfall_in,corrected_ET, mV)
    use carbon_model_memory

    !
    ! Function updates the soil water status and layer thickness
    ! Soil water profile is updated in turn with evaporative losses,
    ! rainfall infiltration and gravitational drainage
    ! Root layer thickness is updated based on changes in the rooting depth from
    ! the previous step
    !

    implicit none

      type(model_working_variables) :: mV

    ! arguments
    double precision, intent(in) :: rainfall_in   ! rainfall (kgH2O.m-2.day-1)
    double precision, intent(inout) :: Eleaf, Esoil, Esnow ! evapotranspiration estimate (kgH2O.m-2.day-1)
    double precision, intent(out) :: corrected_ET      ! water balance corrected evapotranspiration (kgH2O/m2/day)

    ! local variables
    integer :: day, a
    double precision :: depth_change, water_change, initial_soilwater, balance, mass_check, &
                        Esoil_local, Esnow_local
    double precision, dimension(nos_root_layers) :: avail_flux, evaporation_losses, pot_evap_losses
    !logical :: iter_soil = .true.

    ! set soil water exchanges
    Esoil = 0d0 ; Esnow = 0d0 ; corrected_ET = 0d0 ; evaporation_losses = 0d0
    mV%underflow = 0d0 ; mV%runoff = 0d0 ; mV%infiltrated = 0d0 ; mV%water_grav_flow = 0d0 ; pot_evap_losses = 0d0
    initial_soilwater = 1d3 * sum(mV%soil_waterfrac(1:nos_soil_layers) * mV%layer_thickness(1:nos_soil_layers))

    !! Assume leaf transpiration is drawn from the soil based on the
    !! update_fraction estimated in calculate_Rtot
    pot_evap_losses = Eleaf * mV%uptake_fraction
    !! Assume all soil evaporation comes from the soil surface only
    !pot_evap_losses(1) = pot_evap_losses(1) + Esoil

! Conditions under which iterative solution should not be needed...
! Scenario 1
!          (i) Soil layers at or below field capacity, therefore no drainage
!         (ii) The existing water supply and rainfall can support evaporative demanded by the canopy and soil
!     Outcome: Extract all needed water, potentially leaving soil in negetative status, followed by infilatration.
!              Allow for drainage if soil is above field capacity as a result of this proecss
! Scenario 2
!          (i) Soil layers ABOVE field capacity, therefore THERE is drainage
!         (ii) The existing water supply and rainfall can support evaporative demanded by the canopy and soil
!     Outcome: Extract allow water and add all infiltration into the soil.
!              Allow for drainage in the final instance as strongly exponential drainage flow should negate time difference.
!              NOTE: that this may bias between runoff and underflow estimation
! UPDATE (19/04/2024): Soil layer above or near field capacity do require daily iteration or other solution at a later time
!                      due to the non-linear effects on soil surface evaporation as a result of drythick development

    ! to allow for smooth water balance integration carry this out at daily time step
    ! NOTE: Should inestigate conditions where we don't need to do all the days.
    !       For example, can we check how the various terms are changing and whether
    !       we can assume an average for the rest of the time step?
    do day = 1, nint(mV%days_per_step)
       ! Possible conditions for avoiding looping all days
       ! 1) When drythick == min_drythick and rainfall_in > Esoil in first day
       ! 2) If initially, drythick > min_drythick, but rainfall_in > Esoil, iteration can stop once drythick == min_drythick
       ! 3) If drythick == min_drythick and rainfall_in < Esoil in first day, iteration can still be avoided if the time step multiple does not result in drythick >  min_drythick
       ! 4)

       !!!!!!!!!!
       ! Evaporative losses
       !!!!!!!!!!

       ! Estimate drythick for the current step
       !drythick = max(min_drythick, top_soil_depth * max(0d0,(1d0 - (soil_waterfrac(1) / field_capacity(1)))))
       mV%drythick = max(min_drythick, top_soil_depth * max(0d0,(1d0 - (mV%soil_waterfrac(1) / mV%porosity(1)))))
       ! Soil surface (kgH2O.m-2.day-1)
       call calculate_soil_evaporation(Esoil_local, mV)

       ! If snow present assume that soil evaporation is sublimation of soil first
       if (Esoil_local > 0d0 .and. mV%snow_storage > 0d0) then
           if (mV%snow_storage > Esoil_local) then
               Esnow_local = Esoil_local
               mV%snow_storage = mV%snow_storage - Esnow_local
               Esoil_local = 0d0
           else
               Esnow_local = mV%snow_storage
               Esoil_local = Esoil_local - Esnow_local
               mV%snow_storage = 0d0
           end if
       else
           Esnow_local = 0d0
       end if
       ! Accumulate overall time step soil and snow evaporation
       Esoil = Esoil + Esoil_local
       Esnow = Esnow + Esnow_local

       ! load transpiration losses from soil profile
       evaporation_losses = pot_evap_losses
       ! Update with current daily estimate of soil evaporation losses, surface only
       evaporation_losses(1) = evaporation_losses(1) + Esoil_local

       ! can not evaporate from soil more than is available (m -> mm)
       ! NOTE: This is due to the fact that both soil evaporation and transpiration
       !       are drawing from the same water supply.
       avail_flux = mV%soil_waterfrac(1:nos_root_layers) * mV%layer_thickness(1:nos_root_layers) * 1d3
       do a = 1, nos_root_layers ! note: timed comparison between "where" and do loop supports do loop for smaller vectors
          if (evaporation_losses(a) > avail_flux(a)) evaporation_losses(a) = avail_flux(a) * 0.9999d0
       end do
       ! this will update the ET estimate outside of the function
       ! days_per_step corrections happens outside of the loop below
       corrected_ET = corrected_ET + sum(evaporation_losses)

       ! adjust water already committed to evaporation
       ! convert kg.m-2 (or mm) -> Mg.m-2 (or m)
       mV%soil_waterfrac(1:nos_root_layers) = mV%soil_waterfrac(1:nos_root_layers) &
                                         + ((-evaporation_losses(1:nos_root_layers)*1d-3) / mV%layer_thickness(1:nos_root_layers))

       ! Correct for dew formation; any water above porosity in the top layer is assumed runoff
       ! NOTE: layer_thickness * 1d3 scales between m3/m3 to kg/m2
       ! Worth investigating whether this term is actually important, dew formation in unlikely
       ! to be large. Therefore, maybe a conditional statment and calculation unrequired
       if (mV%soil_waterfrac(1) > mV%porosity(1)) then
           mV%runoff = mV%runoff + ((mV%soil_waterfrac(1)-mV%porosity(1)) * mV%layer_thickness(1) * 1d3)
           mV%soil_waterfrac(1) = mV%porosity(1)
       endif

       !!!!!!!!!!
       ! Rainfall infiltration drainage
       !!!!!!!!!!

       ! Determine infiltration from rainfall (kgH2O/m2/day),
       ! if rainfall is probably liquid / soil surface is probably not frozen
       ! reset soil water change variable
       mV%waterchange = 0d0
       call infiltrate(rainfall_in, mV)
       ! update soil profiles. Convert fraction into depth specific values
       ! (rather than m3/m3) then update fluxes
       mV%soil_waterfrac(1:nos_soil_layers) = mV%soil_waterfrac(1:nos_soil_layers) &
                                         + (mV%waterchange(1:nos_soil_layers) / mV%layer_thickness(1:nos_soil_layers))
       ! soil waterchange variable reset in gravitational_drainage()

       !!!!!!!!!!
       ! Gravitational drainage
       !!!!!!!!!!

       ! Determine drainage flux between surface -> sub surface
       call gravitational_drainage(1, mV)

    end do ! days_per_step

    ! apply time step correction kgH2O/m2/step -> kgH2O/m2/day
    mV%water_grav_flow = mV%water_grav_flow * mV%days_per_step_1
    mV%infiltrated = mV%infiltrated * mV%days_per_step_1
    corrected_ET = corrected_ET * mV%days_per_step_1
    mV%underflow = mV%underflow * mV%days_per_step_1
    mV%runoff = mV%runoff * mV%days_per_step_1
    Esoil = Esoil * mV%days_per_step_1
    Esnow = Esnow * mV%days_per_step_1

    ! Based on the soil mass balance corrected_ET, make assumptions to correct Eleaf and Esnow
    balance = (corrected_ET-Eleaf) / Esoil
    if (balance < 1d0 .and. balance > 0d0) Esoil = Esoil * balance

    ! Update corrected_ET with snow sublimation
    corrected_ET = corrected_ET + Esnow

    !!!!!!!!!!
    ! Update soil layer thickness
    !!!!!!!!!!

    ! Determine the root depth below which we will consider layer thickness changes,
    ! and reset the water_change variable, used for maintaining mass balance.
    depth_change = top_soil_depth+min_layer ; water_change = 0
    ! if roots extent down into the bucket
    if (mV%root_reach > depth_change) then

        !!!!!!!!!!
        ! Soil profile is within the bucket layer (layer 3)
        ! therefore we must expand the soil layer (layer 2).
        !!!!!!!!!!

        if (mV%previous_depth > depth_change) then
            ! how much has root depth extended since last step?
            depth_change = mV%root_reach - mV%previous_depth
        else
            ! how much has root depth extended since last step?
            depth_change = mV%root_reach - depth_change
        endif

        ! if there has been an increase
        if (depth_change > 0.05d0) then

            ! determine how much water (mm) is within the new volume of soil
            water_change = mV%soil_waterfrac(nos_soil_layers) * depth_change
            ! now assign that new volume of water to the deep rooting layer
            mV%soil_waterfrac(nos_root_layers) = ((mV%soil_waterfrac(nos_root_layers)*mV%layer_thickness(nos_root_layers))+water_change) &
                                            / (mV%layer_thickness(nos_root_layers)+depth_change)

            ! explicitly update the soil profile if there has been rooting depth
            ! changes
            mV%layer_thickness(1) = top_soil_depth
            mV%layer_thickness(2) = mV%root_reach - top_soil_depth
            mV%layer_thickness(3) = mV%max_depth - sum(mV%layer_thickness(1:2))

            ! keep track of the previous rooting depth
            mV%previous_depth = mV%root_reach

        else if (depth_change < -0.05d0) then

            ! make positive to ensure easier calculations
            depth_change = -depth_change

            ! determine how much water is lost from the old volume of soil
            water_change = mV%soil_waterfrac(nos_root_layers) * depth_change
            ! now assign that new volume of water to the deep rooting layer
            mV%soil_waterfrac(nos_soil_layers) = ((mV%soil_waterfrac(nos_soil_layers)*mV%layer_thickness(nos_soil_layers))+water_change) &
                                            / (mV%layer_thickness(nos_soil_layers)+depth_change)

            ! explicitly update the soil profile if there has been rooting depth
            ! changes
            mV%layer_thickness(1) = top_soil_depth
            mV%layer_thickness(2) = mV%root_reach - top_soil_depth
            mV%layer_thickness(3) = mV%max_depth - sum(mV%layer_thickness(1:2))

            ! keep track of the previous rooting depth
            mV%previous_depth = mV%root_reach

        else

            ! keep track of the previous rooting depth
            mV%previous_depth = mV%previous_depth

        end if ! depth change

    else if (mV%root_reach < depth_change .and. mV%previous_depth > depth_change) then

        !!!!!!!!!!
        ! Model has explicitly contracted from the bucket layer
        !!!!!!!!!!

        ! In this circumstance we want to return the soil profile to it's
        ! default structure with a minimum sized third layer
        depth_change = mV%previous_depth - depth_change

        ! determine how much water is lost from the old volume of soil
        water_change = mV%soil_waterfrac(nos_root_layers) * depth_change
        ! now assign that new volume of water to the deep rooting layer
        mV%soil_waterfrac(nos_soil_layers) = ((mV%soil_waterfrac(nos_soil_layers)*mV%layer_thickness(nos_soil_layers))+water_change) &
                                        / (mV%layer_thickness(nos_soil_layers)+depth_change)

        ! explicitly update the soil profile if there has been rooting depth
        ! changes
        mV%layer_thickness(1) = top_soil_depth
        mV%layer_thickness(2) = min_layer
        mV%layer_thickness(3) = mV%max_depth - sum(mV%layer_thickness(1:2))

        ! keep track of the previous rooting depth
        mV%previous_depth = min_layer

    else ! root_reach > (top_soil_depth + min_layer)

        ! if we are outside of the range when we need to consider rooting depth changes keep track in case we move into a zone when we do
        mV%previous_depth = mV%previous_depth

    endif ! root reach beyond top layer

    ! Update soil water potential
    call soil_water_potential(mV)

!    ! check water balance
!    balance = (rainfall_in - corrected_ET - underflow - runoff) * days_per_step
!    balance = balance &
!            - (sum(soil_waterfrac(1:nos_soil_layers) * layer_thickness(1:nos_soil_layers) * 1d3) &
!            - initial_soilwater)
!
!    if (abs(balance) > 1d-6 .or. soil_waterfrac(1) < -1d-6) then
!        print*,"Soil water miss-balance (mm)",balance
!        print*,"Initial_soilwater (mm) = ",initial_soilwater
!        print*,"Final_soilwater (mm) = ",sum(soil_waterfrac(1:nos_soil_layers) * layer_thickness(1:nos_soil_layers) * 1d3)
!        print*,"State balance = ",sum(soil_waterfrac(1:nos_soil_layers)*layer_thickness(1:nos_soil_layers)*1d3)-initial_soilwater
!        print*,"Flux balance = ",(rainfall_in - corrected_ET - underflow - runoff) * days_per_step
!        print*,"Top soilwater (fraction)",soil_waterfrac(1)
!        print*,"Rainfall (mm/step)",rainfall_in,"ET",corrected_ET,"underflow",underflow,"runoff",runoff
!        print*,"Rainfall (kgH2O/m2/s)",rainfall
!        print*,"Soil Water Fraction = ",soil_waterfrac
!    end if ! abs(balance) > 1d-10

    ! explicit return needed to ensure that function runs all needed code

    return

  end subroutine calculate_update_soil_water
  !
  !-----------------------------------------------------------------
  !
  subroutine infiltrate(rainfall_in, mV)
    use carbon_model_memory

    ! Takes surface_watermm and distributes it among top !
    ! layers. Assumes total infilatration in timestep.   !
    ! NOTE: Assumes that any previous water movement due to infiltration and evaporation
    !       has already been updated in soil mass balance

    implicit none

      type(model_working_variables) :: mV

    ! arguments
    double precision, intent(in) :: rainfall_in ! rainfall (kgH2O.m-2.day-1)

    ! local argumemts
    integer :: i
    double precision :: &
            pot_runoff, & ! Potential rate of runoff (MgH2O/m2/day)
                   add, & ! surface water available for infiltration (m)
                 wdiff    ! available space in a given soil layer for water to fill (m)

    ! Parameters
    double precision, parameter :: beta_s = 0.5d0 ! Infiltration enhancement factor (See Best et al., 2011, Table 6)
                                                  ! i.e. if this was == 1 then there would be more infiltration
    ! convert rainfall water from mm -> m (or kgH2O.m-2.day-1 -> MgH2O.m-2.day-1)
    add = rainfall_in * 1d-3

    ! Estimate the potential infiltration rate, assumed as half the conductivity rate at porosity.
    ! This is consistent with JULES, ORCHIDEE, CLM models as examples.
    ! The *0.5 is the Beta_s factor taken from JULES, but other values are used in the other models.
    call calculate_soil_conductivity(1,mV%porosity(1),pot_runoff, mV)
    ! Scaled to 2.5 hours, on the mean number of hours over which rainfall typically occurs.
    ! This is based on reported precipitation time series from global fluxnet2025 database.
    pot_runoff = pot_runoff * seconds_per_day * (2.5d0/24d0) * beta_s
    pot_runoff = add * exp(-(pot_runoff / add)) ! potential runoff
    ! Update cumulative runoff and substract from the rainfall addition
    mV%runoff = mV%runoff + (pot_runoff*1d3) ; add = add - pot_runoff ; pot_runoff = 0d0

    ! Loop through soil layers to drain
    do i = 1 , nos_soil_layers

       ! is the input of water greater than available space
       ! if so fill and subtract from input and move on to the next
       ! layer determine the available pore space in current soil layer
       wdiff = max(0d0,(mV%porosity(i)-mV%soil_waterfrac(i)) * mV%layer_thickness(i))

       if (add > wdiff) then
           ! if so fill and subtract from input and move on to the next layer
           mV%waterchange(i) = mV%waterchange(i) + wdiff
           add = add - wdiff
       else
           ! otherwise infiltate all in the current layer
           mV%waterchange(i) = mV%waterchange(i) + add
           add = 0d0 ; exit
       end if

    end do ! nos_soil_layers

    ! if after all of this we have some water left assume it is runoff (kgH2O.m-2.day-1)
    ! NOTE that runoff is reset outside of the daily soil loop
    mV%runoff = mV%runoff + (add * 1d3)
    mV%infiltrated = mV%infiltrated + (mV%waterchange(1:nos_soil_layers) * 1d3)

  end subroutine infiltrate
  !
  !-----------------------------------------------------------------
  !
!  subroutine gravitational_drainage(time_period_days)
!
!    ! Integrator for soil gravitational drainage.
!    ! Due to the longer time steps undertake by ACM / DALEC and the fact that
!    ! drainage is a concurrent processes we assume that drainage occurs at
!    ! the bottom of the column first creating space into which water can drain
!    ! from the top down. Therefore we draing from the bottom first and then the top.
!    ! NOTE: Assumes that any previous water movement due to infiltration and evaporation
!    !       has already been updated in soil mass balance
!
!    implicit none
!
!    ! arguments
!    integer, intent(in) :: time_period_days
!
!    ! local variables..
!    integer :: t
!    double precision, dimension(nos_soil_layers) :: dx, & ! range between the start and end points of the integration
!                                               halfway, & ! half way point between start and end point of integration
!                                                liquid, & ! liquid water in local soil layer (m3/m3)
!                                         avail_to_flow, & ! liquid content above field capacity (m3/m3)
!                                               iceprop, & ! fraction of soil layer which is ice
!                                          pot_drainage    ! estimats of time step potential drainage rate (m/s)
!    double precision  :: tmp1,tmp2,tmp3 &
!                                 ,unsat & ! unsaturated pore space in soil_layer below the current (m3/m3)
!                                ,change   ! absolute volume of water drainage in current layer (m3/day)
!
!    ! calculate soil ice proportion; at the moment
!    ! assume everything liquid
!    iceprop = 0d0
!
!    ! except the surface layer in the mean daily temperature is < 0oC
!    if (meant < 1d0) iceprop(1) = 1d0
!
!    ! zero water fluxes
!    waterchange = 0d0
!
!    ! underflow and water_grav_flow are tracked in kgH2O/m2/day but estimated here in MgH2O/m2/day
!    ! therefore we must convert
!    underflow = underflow * 1d-3
!    water_grav_flow = water_grav_flow * 1d-3
!
!    ! estimate potential drainage rate for the current time period
!    liquid = soil_waterfrac(1:nos_soil_layers) * ( 1d0 - iceprop(1:nos_soil_layers) )
!    ! estimate how much liquid is available to flow
!    avail_to_flow = liquid - field_capacity(1:nos_soil_layers)
!    ! trapezium rule scaler and the half-way point between current and field capacity
!    dx = avail_to_flow*0.5d0 ; halfway = liquid - dx
!    do t = 1, nos_soil_layers
!       if (avail_to_flow(t) > 0d0) then
!           ! Trapezium rule for approximating integral of drainage rate
!           call calculate_soil_conductivity(t,liquid(t),tmp1)
!           call calculate_soil_conductivity(t,field_capacity(t),tmp2)
!           call calculate_soil_conductivity(t,halfway(t),tmp3)
!           pot_drainage(t) = 0.5d0 * dx(t) * ((tmp1 + tmp2) + 2d0 * tmp3)
!       else
!           ! We are at field capacity currently even after rainfall has been infiltrated.
!           ! Assume that the potential drainage rate is that at field capacity
!           call calculate_soil_conductivity(t,field_capacity(t),pot_drainage(t))
!       endif ! water above field capacity to flow?
!    end do ! soil layers
!    ! Scale potential drainage from per second to per day
!    pot_drainage = pot_drainage * seconds_per_day
!
!    ! Integrate drainage over each day until time period has been reached or
!    ! each soil layer has reached field capacity
!    t = 1
!    do while (t < (time_period_days+1) .and. maxval(soil_waterfrac - field_capacity) > vsmall)
!
!       ! Estimate liquid content and how much is available to flow / drain
!       avail_to_flow = ( soil_waterfrac(1:nos_soil_layers) * (1d0 - iceprop(1:nos_soil_layers)) ) &
!                     - field_capacity(1:nos_soil_layers)
!
!       ! ...then from the top down
!       do soil_layer = 1, nos_soil_layers
!
!          ! initial conditions; i.e. is there liquid water and more water than
!          ! layer can hold
!          if (avail_to_flow(soil_layer) > 0d0 .and. soil_waterfrac(soil_layer+1) < porosity(soil_layer+1)) then
!
!              ! Unsaturated volume of layer below (m3 m-2)
!              unsat = ( porosity(soil_layer+1) - soil_waterfrac(soil_layer+1) ) &
!                    * layer_thickness(soil_layer+1) / layer_thickness(soil_layer)
!              ! Restrict potential rate calculate above for the available water
!              ! and available space in the layer below.
!              ! NOTE: * layer_thickness(soil_layer) converts units from m3/m2 -> (m3)
!              change = min(unsat,min(pot_drainage(soil_layer),avail_to_flow(soil_layer))) * layer_thickness(soil_layer)
!              ! update soil layer below with drained liquid
!              waterchange( soil_layer + 1 ) = waterchange( soil_layer + 1 ) + change
!              waterchange( soil_layer     ) = waterchange( soil_layer     ) - change
!              ! Also track only the positive flows from one layer to another (MgH2O/m2/day)
!              water_grav_flow(soil_layer) = water_grav_flow(soil_layer) + change
!
!          end if ! some liquid water and drainage possible
!
!       end do ! soil layers
!
!       ! update soil water profile
!       soil_waterfrac(1:nos_soil_layers) = soil_waterfrac(1:nos_soil_layers) &
!                                         + (waterchange(1:nos_soil_layers)/layer_thickness(1:nos_soil_layers))
!       ! estimate drainage from bottom of soil column (MgH2O/m2/day)
!       ! NOTES: that underflow is reset outside of the daily soil loop
!       underflow = underflow + waterchange(nos_soil_layers+1)
!
!       ! Reset now we have moves that liquid
!       waterchange = 0d0
!       ! integerate through time period
!       t = t + 1
!
!    end do ! while condition
!
!    ! convert underflow and water_grav_flow from MgH2O/m2/day -> kgH2O/m2/day
!    underflow = underflow * 1d3
!    water_grav_flow = water_grav_flow * 1d3
!
!  end subroutine gravitational_drainage
  !
  !-----------------------------------------------------------------
  !
  subroutine gravitational_drainage(time_period_days, mV)
    use carbon_model_memory

    ! Integrator for soil gravitational drainage.
    ! Due to the longer time steps undertake by ACM / DALEC and the fact that
    ! drainage is a concurrent processes we assume that drainage occurs at
    ! the bottom of the column first creating space into which water can drain
    ! from the top down. Therefore we draing from the bottom first and then the top.
    ! NOTE: Assumes that any previous water movement due to infiltration and evaporation
    !       has already been updated in soil mass balance

    implicit none

      type(model_working_variables) :: mV

    ! arguments
    integer, intent(in) :: time_period_days

    ! local variables..
    integer :: t, d, s
    double precision, dimension(nos_soil_layers+1) :: soil_waterfrac_local ! local copy of soil water content (m3.m-3)
    double precision, dimension(nos_soil_layers) :: dx, & ! range between the start and end points of the integration
                                               halfway, & ! half way point between start and end point of integration
                                         avail_to_flow, & ! liquid content above field capacity (m3/m3)
                                       liquid_fraction    ! fraction of soil layer which is liquid
    double precision, parameter :: integration_step_seconds = 1800d0*16d0, &
                                   sixth = (1d0/6d0)
    double precision  :: tmp1,tmp2,tmp3, tmp4 &
                             ,pot_drainage_k0 & ! estimates of time step potential drainage rate (m/s)
                             ,pot_drainage_k1 &
                             ,pot_drainage_k2 & 
                             ,pot_drainage_k3 &
                             ,pot_drainage_k4 &
                                      ,liquid & ! liquid water in local soil layer (m3/m3)                             
                                       ,unsat & ! unsaturated pore space in soil layer below the current (m3/m3)
                                      ,change   ! absolute volume of water drainage in current layer (m3/day)

    ! Reset
    liquid_fraction = 1d0 ; mV%waterchange = 0d0 ; d = 1

    ! calculate soil ice proportion;
    ! except the surface layer in the mean daily temperature is < 0oC
    if (mV%meant < 1d0) liquid_fraction(1) = 0d0

    ! underflow and water_grav_flow are tracked in kgH2O/m2/day but estimated here in MgH2O/m2/day
    ! therefore we must convert
    mV%underflow = mV%underflow * 1d-3
    mV%water_grav_flow = mV%water_grav_flow * 1d-3

    ! Estimate liquid content and how much is available to flow / drain
    avail_to_flow = (mV%soil_waterfrac(1:nos_soil_layers) * liquid_fraction(1:nos_soil_layers) ) &
                  - mV%field_capacity(1:nos_soil_layers)

    ! Integrate drainage over each 30 min within day until time period has been reached or
    ! each soil layer has reached field capacity
    do while (d < 4 .and. maxval(avail_to_flow) > vsmall) 

        ! ...then from the top down
        do s = 1, nos_soil_layers

           ! Determine whethere we have any liquid water in the current layer available to flow and the layer below is 
           ! able to accept any water (i.e. is less than porosity).
           if (avail_to_flow(s) > 0d0 .and. mV%soil_waterfrac(s+1) < mV%porosity(s+1)) then

               ! Reset the soil water change variable, this must be done at the beginning of the loop
               ! to ensure that its content are left for extracting the underflow flux after the loop.
               mV%waterchange = 0d0

               !! Implement an explicit 4th order Runge-Kutta approach, to estimate the aggregated flux
               !
               ! Mathematical basis
               ! ------------------
               ! RK4 approximates the solution of  dS/dt = f(t, S)  over a step [t, t+dt]
               ! by evaluating f at four points within the interval and forming a weighted
               ! average of the resulting slopes:
               !
               !   k1 = f(t,          S)                   slope at the start
               !   k2 = f(t + dt/2,   S + dt/2 * k1)       slope at the midpoint (Euler)
               !   k3 = f(t + dt/2,   S + dt/2 * k2)       slope at the midpoint (corrected)
               !   k4 = f(t + dt,     S + dt   * k3)       slope at the end
               !
               !   S_new = S + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
               !
               ! The weighting (1, 2, 2, 1)/6 gives Simpson's-rule-like accuracy.

               ! Load the current soil water content into a local variable to be updated
               soil_waterfrac_local = mV%soil_waterfrac
               ! Estimate the local liquid content in the current layer
               liquid = soil_waterfrac_local(s) * liquid_fraction(s) 

               !! Estimate K1 - rate at the start
               ! Estimate the soil water conductance (k1) at the start of the step
               call calculate_soil_conductivity(s,liquid,pot_drainage_k1, mV)

               !! Estimate K2 - rate at the mid-point
               ! Increment the soil layer using the initial rate (k1) for half the time step
               pot_drainage_k0 = pot_drainage_k1 * (integration_step_seconds * 0.5d0)
               call gravitational_drainage_local_update(s,soil_waterfrac_local,pot_drainage_k0,liquid_fraction(s), &
                                                        mV%layer_thickness,mV%field_capacity(s),mV%porosity)
               ! Estimate the soil water conductance at this new state
               liquid = (soil_waterfrac_local(s) * liquid_fraction(s) )      
               call calculate_soil_conductivity(s,liquid,pot_drainage_k2, mV)

               !! Estimate K3 - rate at the corrected mid-point
               ! Load the current soil water content into a local variable to be updated
               soil_waterfrac_local = mV%soil_waterfrac
               ! Increment the soil layer using the mid-point rate (k2) for half the time step
               pot_drainage_k0 = pot_drainage_k2 * (integration_step_seconds * 0.5d0)
               call gravitational_drainage_local_update(s,soil_waterfrac_local,pot_drainage_k0,liquid_fraction(s), &
                                                        mV%layer_thickness,mV%field_capacity(s),mV%porosity)
               ! Estimate the soil water conductance at this new state
               liquid = (soil_waterfrac_local(s) * liquid_fraction(s) )         
               call calculate_soil_conductivity(s,liquid,pot_drainage_k3, mV)

               !! Estimate K4 - rate at the end
               ! Load the current soil water content into a local variable to be updated
               soil_waterfrac_local = mV%soil_waterfrac
               ! Increment the soil layer using the mid-point rate (k3) for the full time step
               pot_drainage_k0 = pot_drainage_k3 * integration_step_seconds
               call gravitational_drainage_local_update(s,soil_waterfrac_local,pot_drainage_k0,liquid_fraction(s), &
                                                        mV%layer_thickness,mV%field_capacity(s),mV%porosity)
               ! Estimate the soil water conductance at this new state
               liquid = (soil_waterfrac_local(s) * liquid_fraction(s) )         
               call calculate_soil_conductivity(s,liquid,pot_drainage_k4, mV)

               ! Calculate the Simpson's rule weighted average of the rates to estimate the effective average
               ! Note the multiplication to the step time period.
               pot_drainage_k0 = integration_step_seconds * sixth &
                               * (pot_drainage_k1 + 2d0*pot_drainage_k2 + 2d0*pot_drainage_k3 + pot_drainage_k4)

               !! Now update soils for real

               ! Unsaturated volume of layer below (m3 m-2)
               unsat = ( mV%porosity(s+1) - mV%soil_waterfrac(s+1) ) * mV%layer_thickness(s+1) / mV%layer_thickness(s)
               ! Restrict potential rate calculate above for the available water
               ! and available space in the layer below.
               ! NOTE: * layer_thickness(s) converts units from m3/m2 -> (m3)
               change = min(unsat,min(pot_drainage_k0,avail_to_flow(s))) * mV%layer_thickness(s)
               ! update soil layer below with drained liquid
               mV%waterchange( s + 1 ) = mV%waterchange( s + 1 ) + change
               mV%waterchange( s     ) = mV%waterchange( s     ) - change
               ! Also track only the positive flows from one layer to another (MgH2O/m2/day)
               mV%water_grav_flow(s) = mV%water_grav_flow(s) + change

               ! Update the current and below layer, note to avoid a min() bound being used we are allowing the core layer to be updated to. 
               ! This MUST be corrected outside of this loop back to the field capacity
               mV%soil_waterfrac(s:(s+1)) = mV%soil_waterfrac(s:(s+1)) + (mV%waterchange(s:(s+1))/mV%layer_thickness(s:(s+1)))

               ! Update the available liquid content for drainage
               avail_to_flow = (mV%soil_waterfrac(1:nos_soil_layers) * liquid_fraction(1:nos_soil_layers) ) &
                             - mV%field_capacity(1:nos_soil_layers)

           end if ! some liquid water and drainage possible

        end do ! soil layers

        ! Return the core layer back to field capacity
        mV%soil_waterfrac(nos_soil_layers+1) = mV%field_capacity(nos_soil_layers)
        ! Extract drainage from bottom of soil column (MgH2O/m2/day)
        ! NOTES: that underflow is reset outside of the daily soil loop
        mV%underflow = mV%underflow + mV%waterchange(nos_soil_layers+1)

        ! Reset now we have moves that liquid
        mV%waterchange = 0d0
        ! integerate through time period
        d = d + 1

    end do ! while condition

    ! convert underflow and water_grav_flow from MgH2O/m2/day -> kgH2O/m2/day
    mV%underflow = mV%underflow * 1d3
    mV%water_grav_flow = mV%water_grav_flow * 1d3

  end subroutine gravitational_drainage
  !
  !-----------------------------------------------------------------
  !
  subroutine gravitational_drainage_local_update(s,soil_waterfrac_local,pot_drainage,liquid_fraction, &
                                                 layer_thickness_local,field_capacity_local,porosity_local)
    use carbon_model_memory

     ! Subroutine will update a local copy of the soil water fraction for each soil layer

     ! Arguments
     integer, intent(in) :: s
     double precision, intent(in) :: liquid_fraction, & ! Fraction of current layer assumed to be liquid (0-1)
                                        pot_drainage, & ! Proprosed drainage for current layer (m/step)
                                field_capacity_local    ! Field capacity of the current layer (m3/m3)
     double precision, dimension(nos_soil_layers+1), intent(inout) :: soil_waterfrac_local ! Local copy of the soil water fraction (m3/m3)
     double precision, dimension(nos_soil_layers+1), intent(in) :: layer_thickness_local, & ! Local copy of the soil layer thickness (m)
                                                                          porosity_local    ! Local copy of the soil layer porosities (m3/m3)

     ! Local variables
     double precision :: unsat, change, avail_to_flow
     double precision, dimension(nos_soil_layers+1) :: waterchange_local

     ! Initialise
     waterchange_local = 0d0 ; change = 0d0 ; unsat = 0d0 ; avail_to_flow = 0d0

     ! Determine how much liquid water is available to flow in the current profile
     avail_to_flow = (soil_waterfrac_local(s) * liquid_fraction ) - field_capacity_local

     ! Determine whethere we have any liquid water in the current layer available to flow and the layer below is 
     ! able to accept any water (i.e. is less than porosity).
     if (avail_to_flow > 0d0 .and. soil_waterfrac_local(s+1) < porosity_local(s+1)) then

          ! Determine the unsaturated volume of layer below (m3 m-2)
          unsat = ( porosity_local(s+1) - soil_waterfrac_local(s+1) ) &
                * layer_thickness_local(s+1) / layer_thickness_local(s)
          ! Restrict potential rate calculate above for the available water
          ! and available space in the layer below.
          ! NOTE: * layer_thickness_local(s) converts units from m3/m2 -> (m3)
          change = min(unsat,min(pot_drainage,avail_to_flow)) * layer_thickness_local(s)
          ! update soil layer below with drained liquid
          waterchange_local( s + 1 ) = waterchange_local( s + 1 ) + change
          waterchange_local( s     ) = waterchange_local( s     ) - change

     end if ! some liquid water and drainage possible

    ! Update soil water profile, just for the current layer and below
    soil_waterfrac_local(s:nos_soil_layers) = soil_waterfrac_local(s:nos_soil_layers) &
                                            + (waterchange_local(s:nos_soil_layers)/layer_thickness_local(s:nos_soil_layers))

    ! Return back to user
    return

  end subroutine gravitational_drainage_local_update  
  !
  !-----------------------------------------------------------------
  !
  subroutine soil_porosity(soil_frac_clay,soil_frac_sand, mV)
    use carbon_model_memory

    ! Porosity is estimated from Saxton equations. !

    implicit none

      type(model_working_variables) :: mV

    ! arguments
    double precision, dimension(nos_soil_layers) :: soil_frac_clay &
                                                   ,soil_frac_sand
    ! local variables..
    double precision, parameter :: H = 0.332d0, &
                                 J = -7.251d-4, &
                                  K = 0.1276d0

    ! loop over soil layers..
    mV%porosity(1:nos_soil_layers) = H + J * soil_frac_sand(1:nos_soil_layers) + &
                                  K * log10(soil_frac_clay(1:nos_soil_layers))
    ! then assign same to core layer
    mV%porosity(nos_soil_layers+1) = mV%porosity(nos_soil_layers)

  end subroutine soil_porosity
  !
  !---------------------------------------------------------------------
  !
  subroutine initialise_soils(soil_frac_clay,soil_frac_sand, mV)
    use carbon_model_memory

    !
    ! Subroutine calculate the soil layers field capacities and sets the initial
    ! soil water potential set to field capacity
    !

    implicit none

      type(model_working_variables) :: mV

    ! arguments
    double precision, dimension(nos_soil_layers) :: soil_frac_clay &
                                                   ,soil_frac_sand

    ! local variables
    integer :: i

    ! Include some hardcoded boundaries for the Saxton equations
    ! NOTE: do loop was found to be faster than 'where' for small vectors
    do i = 1, nos_soil_layers
       if (soil_frac_sand(i) < 5d0) soil_frac_sand(i) = 5d0
       if (soil_frac_clay(i) < 5d0) soil_frac_clay(i) = 5d0
       !if (soil_frac_sand(i) > 60d0) soil_frac_sand(i) = 60d0
       if (soil_frac_sand(i) > 70d0) soil_frac_sand(i) = 70d0
    end do
    ! calculate soil porosity (m3/m3)
    call soil_porosity(soil_frac_clay,soil_frac_sand, mV)
    ! calculate field capacity (m3/m3)
    call calculate_field_capacity(mV)

    ! final sanity check for porosity
    do i = 1, nos_soil_layers+1
       if (mV%porosity(i) < (mV%field_capacity(i)+0.05d0)) mV%porosity(i) = mV%field_capacity(i) + 0.05d0
    end do

  end subroutine initialise_soils
  !
  !---------------------------------------------------------------------
  !
  subroutine update_soil_initial_conditions(input_soilwater_frac, mV)
    use carbon_model_memory

    !
    ! Subroutine calculate the soil layers field capacities and sets the initial
    ! soil water potential set to field capacity
    !

    implicit none

      type(model_working_variables) :: mV

    ! arguments
    double precision, intent(in) :: input_soilwater_frac ! initial soil water status as fraction of field capacity

    ! local variables
    integer :: i

    ! Load initial soil water fraction to the dynamic layers
    !soil_waterfrac(1:nos_soil_layers) = input_soilwater_frac
    ! Load initial soil water fraction into the top soil layer...
    mV%soil_waterfrac(1) = input_soilwater_frac
    ! ...and assume that the below dynamic layers have the same water content as a fraction of porosity.
    ! This is to avoid the adverse changes to soil water potential estimates in deeper layers due to
    ! differences in hydraulic characteristics at depth
    mV%soil_waterfrac(2:nos_soil_layers) = (input_soilwater_frac / mV%field_capacity(1)) * mV%field_capacity(2:nos_soil_layers)
    ! Assume that the 'core' soil layer is field capacity
    mV%soil_waterfrac(nos_soil_layers+1) = mV%field_capacity(nos_soil_layers)
    ! calculate initial soil water potential
    call soil_water_potential(mV)

    ! Seperately calculate the soil conductivity as this applies to each layer
    do i = 1, nos_soil_layers
       call calculate_soil_conductivity(i,mV%soil_waterfrac(i),mV%soil_conductivity(i), mV)
    end do ! soil layers
    ! but apply the lowest soil layer to the core as well in initial conditions
    mV%soil_conductivity(nos_soil_layers+1) = mV%soil_conductivity(nos_soil_layers)

  end subroutine update_soil_initial_conditions
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_soil_conductivity(soil_layer,waterfrac,conductivity, mV)
    use carbon_model_memory

    ! Calculate the soil conductivity (m s-1) of water based on soil
    ! characteristics and current water content

    implicit none

      type(model_working_variables) :: mV

    ! arguments
    integer, intent(in) :: soil_layer
    double precision, intent(in) :: waterfrac
    double precision, intent(out) :: conductivity

    ! soil conductivity for the dynamic soil layers (i.e. not including core)
    conductivity = mV%cond1(soil_layer) * exp(mV%cond2(soil_layer)+mV%cond3(soil_layer)/waterfrac)

    ! protection against floating point error
    if (waterfrac < 0.05d0) conductivity = 1d-30

  end subroutine calculate_soil_conductivity
  !
  !------------------------------------------------------------------
  !
  subroutine saxton_parameters(soil_frac_clay,soil_frac_sand, mV)
    use carbon_model_memory

    ! Calculate the key parameters of the Saxton, that is cond1,2,3 !
    ! and potA,B                                                    !

    implicit none

      type(model_working_variables) :: mV

    ! arguments
    double precision, dimension(nos_soil_layers) :: soil_frac_clay &
                                                   ,soil_frac_sand

    ! local variables
    double precision, parameter :: A = -4.396d0,  B = -0.0715d0,   CC = -4.880d-4, D = -4.285d-5, &
                                   E = -3.140d0,  F = -2.22d-3,     G = -3.484d-5, H = 0.332d0,   &
                                   J = -7.251d-4, K = 0.1276d0,     P = 12.012d0,  Q = -7.551d-2, &
                                   R = -3.895d0,  T = 3.671d-2,     U = -0.1103d0, V = 8.7546d-4, &
                                   mult1 = 100d0, mult2 = 2.778d-6

    ! layed out in this manor to avoid memory management issues in module
    ! variables
    mV%potA(1:nos_soil_layers) = A + (B * soil_frac_clay) + &
                             (CC * soil_frac_sand * soil_frac_sand) + &
                              (D * soil_frac_sand * soil_frac_sand * soil_frac_clay)
    mV%potA(1:nos_soil_layers) = exp(mV%potA(1:nos_soil_layers))
    mV%potA(1:nos_soil_layers) = mV%potA(1:nos_soil_layers) * mult1

    mV%potB(1:nos_soil_layers) = E + (F * soil_frac_clay * soil_frac_clay) + &
                                  (G * soil_frac_sand * soil_frac_sand * soil_frac_clay)

    mV%cond1(1:nos_soil_layers) = mult2
    mV%cond2(1:nos_soil_layers) = P + (Q * soil_frac_sand)
    mV%cond3(1:nos_soil_layers) = R + (T * soil_frac_sand) + (U * soil_frac_clay) + &
                                   (V * soil_frac_clay * soil_frac_clay)

    ! assign bottom of soil column value to core
    mV%potA(nos_soil_layers+1)  = mV%potA(nos_soil_layers)
    mV%potB(nos_soil_layers+1)  = mV%potB(nos_soil_layers)
    mV%cond1(nos_soil_layers+1) = mult2
    mV%cond2(nos_soil_layers+1) = mV%cond2(nos_soil_layers)
    mV%cond3(nos_soil_layers+1) = mV%cond3(nos_soil_layers)

  end subroutine saxton_parameters
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_soil_conductance(lm,local_lai,canopy_decay, mV)
    use carbon_model_memory

    ! proceedsure to solve for soil surface resistance based on Monin-Obukov
    ! similarity theory stability correction momentum & heat are integrated
    ! through the under canopy space and canopy air space to the surface layer
    ! references are Nui & Yang 2004; Qin et al 2002
    ! NOTE: conversion to conductance at end

    implicit none

      type(model_working_variables) :: mV

    ! declare arguments
    double precision, intent(in) :: lm, &
                             local_lai, &
                          canopy_decay    ! & ! canopy decay coefficient for soil exchange

    ! local variables
    double precision :: Kh_canht       ! eddy diffusivity at canopy height (m2.s-1)

    ! parameters
    double precision, parameter :: foliage_drag = 0.2d0 ! foliage drag coefficient

    ! calculate eddy diffusivity at the top of the canopy (m2.s-1)
    ! Kaimal & Finnigan 1994; for near canopy approximation
    Kh_canht = vonkarman*mV%ustar*(canopy_height-mV%displacement)

    ! approximation of integral for soil resistance (s/m) and conversion to
    ! conductance (m/s)
!    soil_conductance = ( canopy_height/(canopy_decay*Kh_canht) &
!                       * (exp(canopy_decay*(1d0-(soil_roughl/canopy_height)))- &
!                          exp(canopy_decay*(1d0-((roughl+displacement)/canopy_height)))) ) ** (-1d0)
    mV%soil_conductance = 1d0 / ( canopy_height/(canopy_decay*Kh_canht) &
                             * (exp(canopy_decay*(1d0-(soil_roughl/canopy_height)))- &
                                exp(canopy_decay*(1d0-((mV%roughl+mV%displacement)/canopy_height)))) )

    return

  end subroutine calculate_soil_conductance
  !
  !----------------------------------------------------------------------
  !
  subroutine soil_water_potential (mV)
    use carbon_model_memory

    ! Find SWP without updating waterfrac yet (we do that in !
    ! waterthermal). Waterfrac is m3 m-3, soilwp is MPa.     !

    implicit none

      type(model_working_variables) :: mV

    integer :: i

    ! reformulation aims to remove if statement within loop to hopefully improve
    ! optimisation
    mV%SWP(1:nos_soil_layers) = -0.001d0 * mV%potA(1:nos_soil_layers) &
                           * mV%soil_waterfrac(1:nos_soil_layers)**mV%potB(1:nos_soil_layers)
    ! NOTE: profiling indiates that 'where' is slower for very short vectors
    do i = 1, nos_soil_layers
       if (mV%SWP(i) < -20d0 .or. mV%SWP(i) /= mV%SWP(i)) mV%SWP(i) = -20d0
    end do

  end subroutine soil_water_potential
  !
  !------------------------------------------------------------------
  !
  subroutine z0_displacement(ustar_Uh,local_lai, mV)
    use carbon_model_memory

    ! dynamic calculation of roughness length and zero place displacement (m)
    ! based on canopy height and lai. Raupach (1994)

    implicit none

      type(model_working_variables) :: mV

    ! arguments
    double precision, intent(out) :: ustar_Uh ! ratio of friction velocity over wind speed at canopy top
    double precision, intent(in) :: local_lai
    ! local variables
    double precision  sqrt_cd1_lai
    double precision, parameter :: cd1 = 7.5d0,   & ! Canopy drag parameter; fitted to data
                                    Cs = 0.003d0, & ! Substrate drag coefficient
                                    Cr = 0.3d0,   & ! Roughness element drag coefficient
                          ustar_Uh_max = 0.35d0,  &
                          ustar_Uh_min = 0.05d0,  &
                                    Cw = 2d0,     &  ! Characterises roughness sublayer depth (m)
                                 phi_h = 0.19314718056d0 ! Roughness sublayer influence function;

    ! describes the departure of the velocity profile from just above the
    ! roughness from the intertial sublayer log law


    ! Estimate canopy drag coefficient
    sqrt_cd1_lai = sqrt(cd1 * local_lai)

    ! calculate estimate of ratio of friction velocity / canopy wind speed.
    ! NOTE: under current min LAI and fixed canopy height (9 m) this ratio is
    ! fixed at 0.3
    ustar_Uh = 0.3d0
!    ustar_Uh = max(ustar_Uh_min,min(sqrt(Cs+Cr*local_lai*0.5d0),ustar_Uh_max))
!    ustar_Uh = sqrt(Cs+Cr*local_lai*0.5d0)

    ! calculate displacement (m); assume minimum lai 1.0 or 1.5 as height is not
    ! varied
    mV%displacement = (1d0-((1d0-exp(-sqrt_cd1_lai))/sqrt_cd1_lai))*canopy_height

    ! calculate roughness sublayer influence function;
    ! this describes the departure of the velocity profile from just above the
    ! roughness from the intertial sublayer log law
    ! phi_h = log(Cw)-1d0+Cw**(-1d0) ! DO NOT FORGET TO UPDATE IF Cw CHANGES

    ! finally calculate roughness length, dependant on displacement, friction
    ! velocity and lai.
    mV%roughl = ((1d0-mV%displacement/canopy_height)*exp(-vonkarman*ustar_Uh-phi_h))*canopy_height

    ! sanity check
!    if (roughl /= roughl) then
!        write(*,*)"TLS:  ERROR roughness length calculations"
!        write(*,*)"Roughness lenght", roughl, "Displacement", displacement
!        write(*,*)"canopy height", canopy_height, "lai", lai
!    endif

  end subroutine z0_displacement
  !
  !------------------------------------------------------------------
  !
  subroutine plant_soil_flow(root_layer,root_length,root_mass &
                            ,demand,root_reach_in,transpiration_resistance &
                            ,Rtot_layer, mV)
    use carbon_model_memory

    !
    ! Calculate soil layer specific water flow form the soil to canopy (mmolH2O.m-2.s-1)
    ! Accounting for soil, root and plant resistance, and canopy demand
    !

    ! calculate and accumulate steady state water flux in mmol.m-2.s-1
    ! From the current soil layer given an amount of root within the soil layer.

    implicit none

      type(model_working_variables) :: mV

    ! arguments
    integer, intent(in) :: root_layer
    double precision, intent(in) :: root_length, &
                                      root_mass, &
                                         demand, &
                                  root_reach_in, &
                       transpiration_resistance
    double precision, intent(out) :: Rtot_layer

    ! local arguments
    double precision :: soilR1, soilR2

    ! Estimate soil hydraulic resistance to water flow (MPa m2 s mmol-1)
    ! Note: 1) soil conductivity converted from m.s-1 -> m2.s-1.MPa-1 by head.
    !       2) soil resistance calculation in single line to reduce assignment costs
    soilR1 = ( log(root_radius_1*(root_length*pi)**(-0.5d0)) &
               /(two_pi*root_length*root_reach_in*(mV%soil_conductivity(root_layer)*head_1))) &
             * 1d-9 * mol_to_g_water
    ! Calculates root hydraulic resistance (MPa m2 s mmol-1) in a soil-root zone
    soilR2 = root_resist / (root_mass*root_reach_in)
    ! Estimate the total hydraulic resistance for the layer
    Rtot_layer = transpiration_resistance + soilR1 + soilR2
    ! Track for later diagnosis of the effective canopy leaf water potential
    if (demand > 0d0) mV%conductance_mmolH2OMPam2s(root_layer) = 1d0/Rtot_layer
    ! Estimate the soil to plant flow of water mmolH2O/m2/s
    mV%water_flux_mmolH2Om2s(root_layer) = demand/Rtot_layer

    ! return
    return

  end subroutine plant_soil_flow
  !
  !------------------------------------------------------------------
  !
  pure function arrhenious( a , b , t )
    use carbon_model_memory

    ! The equation is simply...                        !
    !    a * exp( b * ( t - 25.0 ) / ( t + 273.15 ) )  !
    ! However, precision in this routine matters as it !
    ! affects many others. To maximise precision, the  !
    ! calculations have been split.                    !

    implicit none

    ! arguments..
    double precision,intent(in) :: a , b , t
    double precision            :: arrhenious

    arrhenious = a * exp( b * (t - 25d0) / (t + freeze) )

  end function arrhenious
  !
  !----------------------------------------------------------------------
  !
  double precision function opt_max_scaling( max_val, min_val , optimum , kurtosis , current )
    use carbon_model_memory

    ! Estimates a 0-1 scaling based on a skewed guassian distribution with a
    ! given optimum, maximum and kurtosis. Minimum is assumed to be at infinity
    ! (or near enough)

    implicit none

    ! arguments..
    double precision, intent(in) :: max_val, min_val, optimum, kurtosis, current

!    ! local variables..
!    double precision, parameter :: min_val = -1d6

    ! Code with implicit assumption of min bound at infinity
!    if ( current >= max_val ) then
!         opt_max_scaling = 0d0
!    else
!         dummy     = exp( log((max_val - current) / (max_val - optimum)) * kurtosis * (max_val - optimum) )
!         opt_max_scaling = dummy * exp( kurtosis * ( current - optimum ) )
!    end if

    ! Code with explicit min bound
    opt_max_scaling = exp( kurtosis * log((max_val-current)/(max_val-optimum)) * (max_val-optimum) ) &
                    * exp( kurtosis * log((current-min_val)/(optimum-min_val)) * (optimum-min_val) )
    ! Sanity check, allows for overlapping parameter ranges
    if (opt_max_scaling /= opt_max_scaling) opt_max_scaling = 0d0

  end function opt_max_scaling
  !
  !------------------------------------------------------------------
  !
  double precision function water_retention_saxton_eqns( xin , mV)
    use carbon_model_memory

    ! field capacity calculations for saxton eqns !

    implicit none

      type(model_working_variables), optional :: mV

    ! arguments..
    double precision, intent(in) :: xin

    ! local variables..
    double precision :: soil_wp

    ! calculate the soil water potential (kPa)..
    soil_wp = -mV%potA(mV%water_retention_pass) * xin**mV%potB(mV%water_retention_pass)
    water_retention_saxton_eqns = soil_wp + 10d0    ! 10 kPa represents air-entry swp

    return

  end function water_retention_saxton_eqns
  !
  !------------------------------------------------------------------
  !
  double precision function ospolynomial(L,w)
    use carbon_model_memory

    ! Function calculates the day offset for Labile release and leaf turnover
    ! functions

    implicit none

    ! declare input variables
    double precision, intent(in) ::  L, w ! polynomial coefficients and scaling factor

    ! declare local variables
    double precision :: LLog, mxc(7) ! polynomial coefficients and scaling factor

    ! assign polynomial terms
    mxc(1) = (0.000023599784710d0)
    mxc(2) = (0.000332730053021d0)
    mxc(3) = (0.000901865258885d0)
    mxc(4) = (-0.005437736864888d0)
    mxc(5) = (-0.020836027517787d0)
    mxc(6) = (0.126972018064287d0)
    mxc(7) = (-0.188459767342504d0)

    ! load log of leaf / labile turnovers
    LLog = log(L-1d0)

    ! calculate the polynomial function
    ospolynomial = (mxc(1)*LLog**6d0 + mxc(2)*LLog**5d0 + &
                    mxc(3)*LLog**4d0 + mxc(4)*LLog**3d0 + &
                    mxc(5)*LLog**2d0 + mxc(6)*LLog      + mxc(7))*w

    ! back to the user...
    return

  end function ospolynomial
  !
  !---------------------------------------------------------------------
  !
  subroutine update_microbial_activity_and_death(p40,p46,p45,deltat, &
                                                 T_scalar,fast_som,microbial_activity,microbial_death)
    use carbon_model_memory

     ! DecoBio v1.0 (Xenakis & Williams 2014, Comparing microbial and chemical
     ! kinetics for modelling soil organic carbon decomposition using the
     ! DecoChem 1.0 and DecoBio v1.0 models, Geoscientific Model Development,
     ! 7, 1519-1533)

     ! This subroutine estimates the new microbial activity (0-1) 
     ! and mortality rates (day-1) as a function of labile C in the soil.
 
     implicit none

     ! Arguments
     double precision, intent(in) :: p40,p46,p45,deltat, &
                                     T_scalar,fast_som
     double precision, intent(inout) :: microbial_activity ! Microbial activity (0-1)
     double precision, intent(out) :: microbial_death      ! Microbial mortality rate (day-1)

     ! Calculate microbial activity (0-1) in relation labile soil C.
     ! Note this is a simplification of the published model due to DALEC operating at >=daily
     ! timesteps. We assume that activity is in steady state with its potential activity.
     microbial_activity = min(1d0,max(0d0, T_scalar * ( fast_som  / ( fast_som + p40 ) ) ))
     ! Calculate the microbial death rate
     microbial_death = p45 / (1d0 + (p46 * fast_som))

     return

  end subroutine update_microbial_activity_and_death
  !
  !------------------------------------------------------------------
  !
!
!--------------------------------------------------------------------
!
end module CARBON_MODEL_MOD

