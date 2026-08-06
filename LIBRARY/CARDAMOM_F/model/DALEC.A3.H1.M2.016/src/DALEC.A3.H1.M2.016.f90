!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CARbon DAta MOdel fraMework (CARDAMOM) and DALEC terrestrial ecosystem model suite
! CARDAMOM is a Bayesian model-data fusion software framework. CARDAMOM is used to 
! assimilate observations and ecological theory to retrieve parameters for the 
! DALEC suite of intermediate complexity terrestrial ecosystem models. DALEC can be
! used as a fully integrated component of CARDAMOM or independently. 
! Copyright (C) 2024  University of Edinburgh,
!                     Mathew Williams (mat.williams@ed.ac.uk), 
!                     T. Luke Smallman (t.l.smallman@ed.ac.uk), 
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
! This file contains the source code of DALEC.A3.H1.M2.016
!
! This code contains a variant of the Data Assimilation Linked ECosystem (DALEC) model.
! This version of DALEC is derived from the following primary references:
! Bloom & Williams (2015), https://doi.org/10.5194/bg-12-1299-2015.
! Smallman & Williams (2019) https://doi.org/10.5194/gmd-12-2227-2019.
! Thomas et al., (2019), https://doi.org/10.1029/2019MS001679
! Myrgiotis et al., (2022). https://doi.org/10.5194/bg-19-4147-2022
! This code is based on that created by A. A. Bloom (UoE, now at JPL, USA).
! Subsequent modifications by:
! T. L. Smallman (University of Edinburgh, t.l.smallman@ed.ac.uk)
! J. F. Exbrayat (University of Edinburgh)
! V. Myrgiotis (UK Centre for Ecology & Hydrology)
! S. Zhu (University of Edinburgh, now University of Southampton)
! See function / subroutine specific comments for exceptions and contributors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module CARBON_MODEL_MOD

  implicit none

  ! grassland model (DALEC.A3.H1.M2.016)
  ! ----------------------------------------------------------------------------------------------------------------
  ! POOLS:   1.labile 2.foliar 3.root                        ! PARAMETERS:
  !          4.litter 5.som                                  !
  ! ------------------------------------------               ! 1.  Decomposition rate
  ! FLUXES:  1.GPP                                           ! 2.  Fraction of GPP respired 
  ! (daily)  2.temprate                                      ! 3.  GSI sens for leaf growth
  !          3.respiration_auto                              ! 4.  NPP belowground allocation parameter
  !          4.leaf production                               ! 5.  GSI max leaf turnover 
  !          5.labile production                             ! 6.  TOR roots
  !          6.root production                               ! 7.  TOR litter
  !          7.aboveground production                        ! 8.  TOR SOM
  !          8.labile consumption -> leaves                  ! 9.  Temp factor Q10 (1.2-1.6)
  !          9.leaffall factor                               ! 10. Photosynthetic N use efficiency
  !          10.leaf litter production                       ! 11. GSI max labile turnover
  !         *11.woodlitter production                        ! 12. GSI min temperature threshold (K)
  !          12.rootlitter production                        ! 13. GSI max temperature threshold (K)
  !          13.respiration het litter                       ! 14. GSI min photoperiod threshold (sec)
  !          14.respiration het som                          ! 15. LCA - g.C.leaf_m-2
  !          15.litter2som                                   ! 16. C labile (initialization)
  !          16.labrelease factor(leaf growth)               ! 17. C foliar (initialization)
  !         *17.carbon flux due to fire                      ! 18. C roots  (initialization)
  !          18.growing season index                         ! 19. C litter (initialization)
  !          19.animal manure C soil input (per time step)   ! 20. GSI max photoperiod threshold (sec)
  !          20.animal resp co2 (per time step)              ! 21. GSI min VPD threshold (Pa) 
  !          21.animal ch4 (per time step)                   ! 22. GSI max VPD threshold (Pa)
  !          22. labile loss (per time step)                 ! 23. critical GPP for LAI increase (gC.m-2.day-1)
  !          23. foliage loss (per time step)                ! 24. GSI senstivity for leaf senescence 
  ! ------------------------------------------               ! 25. GSI - have I just left a growing state (>1)
  ! MET:     1.run day                                       ! 26. GSI - initial GSI value
  !          2.min T (C)                                     ! 27. DM min lim for grazing (kg.DM.ha-1)
  !          3.max T (C)                                     ! 28. DM min lim for cutting (kg.DM.ha-1) 
  !          4.Radiation (MJ.m-2)                            ! 29. leaf-vs-stem allocation factor
  !          5.CO2 (ppm)                                     ! 30. C SOM (initialization)
  !          6.DOY                                           ! 31. DM demand of animal weight (fraction) 
  !         *7.lagged precip                                 ! 32. Post-grazing labile loss (fraction)
  !          8.cutting/grazing :                             ! 33. Post-cut labile loss (fraction)
  !            - spatial mode = lai removed (m2.m-2)         ! 34. Minimum grazed biomass to allow grazing
  !            - field mode = LSU.ha-1                       
  !         *9.burnt area fraction                           
  !          10.21-day avg min T (K)                        
  !          11.21-day avg photoperiod (sec) 
  !          12.21-day avg VPD (Pa)         
  !         *13.Forest mgmt after clearing
  !         *14.Mean T
  ! ----------------------------------------------------------------------------------------------------------------
  ! NOTES : '*' above means not used/applicable for grasslands 
  !         1 LSU per ha = 1 cow that weighs 650kg and grazes on 1 ha of grassland 
  !         carbon = 0.475 * dry matter 
  !         1 g.C.m-2 = 1 * 0.021 t.DM.ha-1
  !         to compile this .f90 into a python shared object (.so) run: f2py -c DALEC_GRASS.f90 -m DALEC_GRASS
  ! ----------------------------------------------------------------------------------------------------------------
  !                  autotrophic      heterotrophic     loss due to     --->    manure from          
  !                  respiration      respiration       grazing/cutting         grazing livestock       
  !                       ^            ^                   ^                       |       
  !                       |            |                   |                       V      
  !                                                                                   
  !PHOTOSYNTHESIS -----> [0] -------> [0] --------------> [0] <-----------------> [0]     
  !                 GPP         NPP            NEE                   NBE                 
  ! ----------------------------------------------------------------------------------------------------------------

  !!!!!!!!!!!
  ! Authorship contributions
  !
  ! This code contains a variant of the Data Assimilation Linked ECosystem (DALEC) model.
  ! This version of DALEC is derived from the following primary references:
  ! Bloom & Williams (2015), https://doi.org/10.5194/bg-12-1299-2015.
  ! Smallman & Williams (2019) https://doi.org/10.5194/gmd-12-2227-2019.
  ! Thomas et al., (2019), https://doi.org/10.1029/2019MS001679
  ! This code is based on that created by A. A. Bloom (UoE, now at JPL, USA).
  ! Subsequent modifications by:
  ! T. L. Smallman (University of Edinburgh, t.l.smallman@ed.ac.uk)
  ! J. F. Exbrayat (University of Edinburgh)
  ! V. Myrgiotis (UK Centre for Ecology & Hydrology)
  ! S. Zhu (University of Edinburgh)
  ! See function / subroutine specific comments for exceptions and contributors
  !!!!!!!!!!!

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
    allocate(mV%deltat_1(nodays),mV%daylength_hours(nodays),mV%daylength_seconds(nodays), &
             mV%daylength_seconds_1(nodays),mV%rainfall_time(nodays),mV%airt_zero_fraction_time(nodays))

    !
    ! Timing variables which are needed first
    !

    mV%deltat_1 = deltat**(-1d0)

    !
    ! Iteration independent variables using functions and thus need to be in a loop
    !

    ! Generate some generic location specific variables for radiation balance
    !call calculate_radiation_commons(lat,pars(33:38))
    call calculate_radiation_commons(lat, mV)

    ! first those linked to the time period of the analysis
    do n = 1, nodays
       ! check positive values only for rainfall input
       mV%rainfall_time(n) = max(0d0,met(7,n))
       ! Calculate solar declination for the current time step
       mV%declination = calculate_declination(met(6,n))           
       ! calculate daylength in hours and seconds
       call calculate_daylength(mV)       
       mV%daylength_hours(n) = mV%dayl_hours ; mV%daylength_seconds(n) = mV%dayl_seconds 
    end do
    ! Create inverse for day length seconds
    mV%daylength_seconds_1 = mV%daylength_seconds**(-1d0)

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
    ! call update_soil_initial_conditions(pars(24), mV)
    ! save the initial conditions for later
    mV%field_capacity_initial = mV%field_capacity
    mV%porosity_initial = mV%porosity

    ! Initialise the Growing Season Index (GSI)
    ! canopy phenology model

    ! Determine the number of model time steps over which GSI will be lagged
    ! NOTE: 21 days is the default assumption from the original source paper (Jolly et al., 2005)
    mV%gsi_lag_steps = max(2,nint(21d0/mV%mean_days_per_step))
    allocate(mV%gsi_lag_days(mV%gsi_lag_steps),mV%gsi_lag_history(mV%gsi_lag_steps))
    call initialise_gsi(mV%mean_days_per_step,mV%gsi_lag_steps,mV%gsi_lag_days)
     
  end subroutine initialize_mv
  !
  !--------------------------------------------------------------------
  !
  subroutine destroy_mv(mV)
    use carbon_model_memory
    !! deallocate arrays in mV
    type(model_working_variables):: mV
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

    ! The Data Assimilation Linked ECosystem model - DALEC.A3.H1.M2.016.
    ! A water enables version of DALEC, specifically designed for understanding managed grasslands.

    ! This version include grassland management specific functions and lacks a wood / structural C pool

    implicit none

    type(model_working_variables) :: mV

    ! declare input variables
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
                          ,f_root & 
                    ,labile_ratio & 
               ,gsi_lai_reduction & 
                ,available_labile & ! available labile for allocation (gC/m2)
                   ,transpiration & ! kgH2O/m2/day
                 ,soilevaporation & ! kgH2O/m2/day
                  ,wetcanopy_evap & ! kgH2O/m2/day
                 ,snowsublimation   ! kgH2O/m2/day

    ! JFE added 4 May 2018 - combustion efficiencies and fire resilience
    double precision :: burnt_area
    double precision, dimension(7) :: cf,rfac
    ! local deforestation related variables        
    integer :: harvest_management,n

    ! met drivers are:
    ! 1st run day
    ! 2nd min daily temp (oC)
    ! 3rd max daily temp (oC)
    ! 4th Radiation (MJ.m-2.day-1)
    ! 5th CO2 (ppm)
    ! 6th DOY
    ! 7th precipitation (kgH2O.m-2.s-1)
    ! 8th lai loss m2m2/step
    ! 9th burnt area fraction
    ! 10th 21 day average min temperature (oC)
    ! 11th 21 day average photoperiod (seconds)
    ! 12th 21 day average VPD (Pa)
    ! 13th Management practice to accompany any clearing
    ! 14th avg daily temperature (oC)
    ! 15th avg daily wind speed (m.s-1)
    ! 16th vapour pressure deficit (Pa)

    ! POOLS are:
    ! 1 = labile (gC/m2) (initial value: p16)
    ! 2 = foliar (gC/m2) (initial value: p17)
    ! 3 = root   (gC/m2) (initial value: p18)
    ! 4 = litter (gC/m2) (initial value: p19)
    ! 5 = som    (gC/m2) (initial value: p23)

    ! FLUXES are:
    ! 1  = GPP (gC/m2/day)
    ! 2  = temperature rate modifier (unitless)
    ! 3  = autotrophic respiration (gC/m2/day)
    ! 4  = GPP allocation to foliage (gC/m2/day)
    ! 5  = GPP allocation to labile (gC/m2/day)
    ! 6  = GPP allocation to roots (gC/m2/day)
    ! 7  = labile->leaf transfer (gC/m2/day)
    ! 8  = leaf fall factor (fraction/day)
    ! 9  = leaf litter production (gC/m2/day)
    ! 10 = root litter production (gC/m2/day)
    ! 11 = heterotrophic respiration from litter (gC/m2/day)
    ! 12 = heterotrophic respiration from som (gC/m2/day)
    ! 13 = litter decomposition to som (gC/m2/day)
    ! 14 = GSI labile release factor (fraction/day)
    ! 15 = NOT IN USE
    ! 16 = NOT IN USE
    ! 17 = NOT IN USE
    ! 18 = NOT IN USE
    ! 19 = animal manure deposited to soil from grazing (gC/m2/day)
    ! 20 = animal respiration from grazing (gC/m2/day)
    ! 21 = animal methane emissions from grazing (gC/m2/day)
    ! 22 = total biomass extracted per cutting event (gC/m2/day)
    ! 23 = total biomass extracted per grazing event (gC/m2/day)
    ! 24 = NOT IN USE
    ! 25 = labile extracted from cutting (gC/m2/day)
    ! 26 = foliage extracted from cutting (gC/m2/day)
    ! 27 = roots extracted from cutting (gC/m2/day)
    ! 28 = labile added to litter from cutting (gC/m2/day)
    ! 29 = foliage added to litter from cutting (gC/m2/day)
    ! 30 = roots added to litter from cutting (gC/m2/day)
    ! 31 = labile extracted from grazing (gC/m2/day)
    ! 32 = foliage extracted from grazing (gC/m2/day)
    ! 33 = roots extracted from grazing (gC/m2/day)
    ! 34 = labile added to litter from grazing (gC/m2/day)
    ! 35 = foliage added to litter from grazing (gC/m2/day)
    ! 36 = roots added to litter from grazing (gC/m2/day)
    ! 37 = NOT IN USE
    ! 38 = NOT IN USE
    ! 39 = NOT IN USE
    ! 40 = NOT IN USE
    ! 41 = NOT IN USE
    ! 42 = NOT IN USE
    ! 43 = NOT IN USE
    ! 44 = NOT IN USE
    ! 45 = NOT IN USE
    ! 46 = evapotranspiration (kgH2O/m2/day)
    ! 47 = transpiration (kgH2O/m2/day)
    ! 48 = soil evaporation (kgH2O/m2/day)
    ! 49 = wet canopy evaporation (kgH2O/m2/day)
    ! 50 = surface runoff (kgH2O/m2/day)
    ! 51 = drainage from bottom of soil column (kgH2O/m2/day)
    ! 52 = drainage from surface to 2nd soil layer (kgH2O/m2/day)
    ! 53 = infiltration into top soil layer (kgH2O/m2/day)
    ! 54 = fraction of transpiration from 1st rooting layer (0-1)
    ! 55 = fraction of transpiration from 2nd rooting layer (0-1)
    ! 56 = infiltration into middle soil layer (kgH2O/m2/day)
    ! 57 = infiltration into bottom soil layer (kgH2O/m2/day)

    ! PARAMETERS are:
    ! p(1)  = litter decomposition rate (fraction/day)
    ! p(2)  = fraction of GPP as autotrophic respiration (fraction)
    ! p(3)  = canopy GSI phenology gradient threshold
    ! p(4)  = NPP belowground allocation exponential parameter
    ! p(5)  = potential leaf turnover rate (fraction/day)
    ! p(6)  = fine root turnover rate (fraction/day)
    ! p(7)  = litter turnover rate (fraction/day)
    ! p(8)  = som turnover rate, temperature adjusted (fraction/day)
    ! p(9)  = temperature sensitivity of heterotrophic respiration (oC-1)
    ! p(10) = max labile turnover rate to foliage (fraction/day)
    ! p(11) = canopy efficiency (umolC/m2leaf/s)
    ! p(12) = GSI minimum temperature threshold (K)
    ! p(13) = GSI maximum temperature threshold (K)
    ! p(14) = GSI minimum photoperiod threshold (seconds)
    ! p(15) = leaf mass per area LMA (gC/m2)
    ! p(16) = initial labile C pool (gC/m2)
    ! p(17) = initial foliar C pool (gC/m2)
    ! p(18) = initial root C pool (gC/m2)
    ! p(19) = initial litter C pool (gC/m2)
    ! p(20) = GSI maximum photoperiod threshold (seconds)
    ! p(21) = GSI minimum VPD threshold (Pa)
    ! p(22) = GSI maximum VPD threshold (Pa)
    ! p(23) = initial som C pool (gC/m2)
    ! p(24) = root biomass for 50% of maximum rooting depth (gBiomass/m2)
    ! p(25) = maximum rooting depth (m)
    ! p(26) = initial canopy GSI value (0-1)
    ! p(27) = minimum above-ground biomass for grazing to occur (gC/m2)
    ! p(28) = minimum above-ground biomass for cutting to occur (gC/m2)
    ! p(29) = NOT IN USE
    ! p(30) = GPP return on new foliage investment (gC/gC)
    ! p(31) = NOT IN USE
    ! p(32) = post-grazing labile loss fraction (fraction)
    ! p(33) = post-cutting labile loss fraction (fraction)
    ! p(34) = minimum biomass removal for a grazing instance to occur (gC/m2/day)

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

    ! post-removal residues and root death | 0:none 1:all
    mV%foliage_frac_res  = 0.05d0  !Â fraction of removed foliage that goes to litter
    mV%labile_frac_res   = 0.05d0  !Â fraction of removed labile that goes to litter
    mV%roots_frac_res    = 1d0     ! fraction of roots that die which go to litter 
    mV%roots_frac_death  = 0.01d0  !Â fraction of roots that die in response to management
    
    ! How many steps in 2 weeks
    mV%two_week_lag = ceiling(14d0/deltat(1))
    ! How many steps in 4 weeks
    mV%four_week_lag = ceiling(28d0/deltat(1))

    ! Generate some generic location specific variables for radiation balance
    !call calculate_radiation_commons(lat,pars(33:38))
    !call calculate_radiation_commons(lat, mV)

    ! load ACM-GPP-ET parameters
    mV%ceff = pars(11) ! Canopy efficiency (umolC/m2/s)
                    ! This is in the full model the product of Nitrogen use efficiency (umolC/gN/m2leaf)
                    ! and average foliar nitrogen gN/m2leaf
    ! Rooting parameters
    mV%root_k = pars(24) ; mV%max_depth = pars(25)

    ! assigning initial conditions
    POOLS(1,1) = pars(16)
    POOLS(1,2) = pars(17)
    POOLS(1,3) = pars(18)
    POOLS(1,4) = pars(19)
    POOLS(1,5) = pars(23)

    if (.not.allocated(mV%deltat_1)) then 
        write(*,*) "Error - arrays not allocated - probably carbon_model() was called without initialize_mv()"
        STOP 1
    else ! deltat_1 allocated?

        !
        ! Load initial soil water conditions from memory
        !

        mV%total_water_flux = 0d0 ; mV%water_flux_mmolH2Om2s = 0d0
        mV%field_capacity = mV%field_capacity_initial
        mV%porosity = mV%porosity_initial

        ! input initial soil water fraction then
        ! update SWP and soil conductivity accordingly
        call update_soil_initial_conditions( mV)

        ! Initialise the Growing Season Index (GSI)
        ! canopy phenology model

        ! Determine the number of model time steps over which GSI will be lagged
        ! NOTE: 21 days is the default assumption from the original source paper (Jolly et al., 2005)
        mV%gsi_lag_steps = max(2,nint(21d0/mV%mean_days_per_step))
        ! We assume the other initialised variables are stored in memory.

    endif ! deltat_1 allocated

    ! assign our starting value
    mV%gsi_lag_history = pars(26)

    ! 
    ! Begin looping through each time step
    ! 

    ! load some needed module level values
    mV%lai = POOLS(1,2)/pars(15)
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
    mV%root_biomass = mV%fine_root_biomass 
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
                                     0d0,FLUXES(1,46), mV) ! assume no evap or rainfall

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
       mV%lai = POOLS(n,2)/pars(15)
       DIAGS(n,1) = mV%lai

       ! extract timing related values
       mV%dayl_hours = mV%daylength_hours(n)
       mV%dayl_hours_fraction = mV%dayl_hours * 0.04166667d0 ! 1/24 = 0.04166667
       mV%dayl_seconds = mV%daylength_seconds(n) ; mV%dayl_seconds_1 = mV%daylength_seconds_1(n)
       mV%days_per_step = deltat(n) ; mV%days_per_step_1 = mV%deltat_1(n)
       mV%seconds_per_step = seconds_per_day * mV%days_per_step

       !!!!!!!!!!
       ! Adjust snow balance balance based on temperture
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
       mV%root_biomass = mV%fine_root_biomass 
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

       ! Estimate photosynthesis and transpiration
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
!print*,lai,stomatal_conductance,FLUXES(n,1)

       ! Estimate average leaf water potential (MPa) based on effective hydraulic resistance, wSWP and transpiration.
       ! Positive LWPs can be estimated given very small gs and cold temperatures.
       ! Debugging print statements
       !print*,"Estimate LWP"
       !LWP = SWP(1:nos_root_layers) - head*canopy_height &
       !     - transpiration*uptake_fraction(1:nos_root_layers) &
       !     * (dayl_seconds_1/mmol_to_kg_water)/Rcond_layer(1:nos_root_layers)
       DIAGS(n,9) =  min(0d0, mV%wSWP - (head*canopy_height) - (((transpiration*mV%dayl_seconds_1)/mmol_to_kg_water) * mV%Reff))

       ! temprate (i.e. T modified rate of metabolic activity))
       FLUXES(n,2) = exp(pars(9)*0.5d0*(met(3,n)+met(2,n)))
       ! Allocate to autotrophic respiration (gC.m-2.day-1)
       ! TLS: Question, should Ra be split between Rg + Rm?
       FLUXES(n,3) = FLUXES(n,1) * pars(2)

       !! Determine direct allocation of GPP to plant tissues   
       ! Dynamic allocation to roots vs aboveground biomass after Reyes.et.al.2017 (10.1002/2017MS001022)
       ! min/max allocation to roots as fraction of NPP
       f_root = max(min_f_root,min(max_f_root,1d0-exp(-1d0*pars(4)*mV%lai)))
       ! allocation to fine roots 
       FLUXES(n,6) = ( FLUXES(n,1)-FLUXES(n,3) ) * f_root
       ! Direct allocation of allocation of ABG C to leaves
       ! 0.1666667 = 1/6
       ! TLS: Is this line actually valid given the GSI phenology? Should we actually have a stem pool for this?
       ! Then the split is between labile (that supports foliage) and a stem.
       FLUXES(n,4) = 0d0 !( FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,6) ) * (1d0 - (pars(29)*min(1d0,lai*0.1666667d0)))
       ! allocation of ABG C to labile/stem using pars(29)
       ! Ostrem.et.al.2013 (10.1080/09064710.2013.819440)            
       FLUXES(n,5) = FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,6)-FLUXES(n,4)

       ! Accumulate this time steps labile C (gC.m-2.day-1)
       available_labile = POOLS(n,1) + (FLUXES(n,5) * mV%days_per_step)                          
!print*,available_labile,POOLS(n,1),FLUXES(n,5),FLUXES(n,1),lai
       ! Do plant allocation to leaves
       call plant_canopy_phenology(nodays, mV%gsi_lag_steps, n, mV%days_per_step,              & ! Timing
                                   met(10,n),met(11,n),met(12,n),                    & ! GSI forcings
                                   pars(12), pars(13), pars(14), pars(20),           & ! GSI parameters 
                                   pars(21), pars(22), pars(3),                      & !
                                   pars(10), pars(5),                                & ! 
                                   pars(15), pars(30), available_labile, POOLS(n,2), & ! To calculate GPP return
                                   FLUXES(n,1), DIAGS(n,21),                         & ! 
                                   FLUXES(n,14), FLUXES(n,7),                        & ! Output variables
                                   FLUXES(n,8), FLUXES(n,9),                         &
                                   DIAGS(n,18), DIAGS(n,19), DIAGS(n,20),            &
                                   DIAGS(n,17), DIAGS(:,16), mV) 

       ! FLUXES WITH TIME DEPENDENCIES

       ! root litter production = P_root * (1-(1-rootTOR)**deltat)/deltat  
       FLUXES(n,10) = POOLS(n,3)*(1d0-(1d0-pars(6))**mV%days_per_step)*mV%deltat_1(n)

       ! FLUXES WITH TEMP AND TIME DEPENDENCIES

       ! resp het litter = P_litter * (1-(1-GPP_respired*litterTOR)**deltat)/deltat  
       FLUXES(n,11) = POOLS(n,4)*(1d0-(1d0-FLUXES(n,2)*pars(7))**mV%days_per_step)*mV%deltat_1(n)
       ! resp het som = P_som * (1-(1-GPP_respired*somTOR)**deltat)/deltat
       FLUXES(n,12) = POOLS(n,5)*(1d0-(1d0-FLUXES(n,2)*pars(8))**mV%days_per_step)*mV%deltat_1(n)
       ! litter to som = P_litter * (1-(1-dec_rate*temprate)**deltat)/deltat
       FLUXES(n,13) = POOLS(n,4)*(1d0-(1d0-FLUXES(n,2)*pars(1))**mV%days_per_step)*mV%deltat_1(n)

       ! update pools for next timestep

       ! labile pool = labile_pool[â€ -1] + (lab_prod - lab_cons)*deltat
       POOLS(n+1,1) = POOLS(n,1) + (FLUXES(n,5)-FLUXES(n,7))*mV%days_per_step
       ! foliar pool = foliar_pool[â€ -1] + (leaf_prod + lab_prod2 - leaf_litter_prod)*deltat
       POOLS(n+1,2) = POOLS(n,2) + (FLUXES(n,4)+FLUXES(n,7)-FLUXES(n,9))*mV%days_per_step
       ! root pool = root_pool[â€ -1] + (root_prod - root_litter_prod)*deltat
       POOLS(n+1,3) = POOLS(n,3) + (FLUXES(n,6)-FLUXES(n,10))*mV%days_per_step
       ! litter pool = litter_pool[â€ -1] + (leaf_litter_prod + root_litter_prod - resp_het_litter - litter2som)*deltat
       POOLS(n+1,4) = POOLS(n,4) + (FLUXES(n,9)+FLUXES(n,10)-FLUXES(n,11)-FLUXES(n,13))*mV%days_per_step
       ! som pool = som_pool[â€ -1] + (litter2som - resp_het_som)
       POOLS(n+1,5) = POOLS(n,5) + (FLUXES(n,13)-FLUXES(n,12))*mV%days_per_step

       !!!!!!!!!!
       ! Update soil water balance
       !!!!!!!!!!

       ! add any snow melt to the rainfall now that we have already dealt with the canopy interception
       mV%rainfall = mV%rainfall + mV%snow_melt
       ! do mass balance (i.e. is there enough water to support ET)
       call calculate_update_soil_water(transpiration,soilevaporation,snowsublimation, &
                                        ((mV%rainfall-mV%intercepted_rainfall)*seconds_per_day) &
                                       ,FLUXES(n,46), mV)
       ! now that soil mass balance has been updated we can add the wet canopy
       ! evaporation (kgH2O.m-2.day-1)
       FLUXES(n,46) = FLUXES(n,46) + wetcanopy_evap
       ! Assign all water variables to output variables (kgH2O/m2/day)
       FLUXES(n,47) = transpiration   ! transpiration  (kgH2O/m2/day)
       FLUXES(n,48) = soilevaporation ! soil evaporation  (kgH2O/m2/day)
       FLUXES(n,49) = wetcanopy_evap  ! wet canopy evaporation  (kgH2O/m2/day)
       FLUXES(n,50) = mV%runoff          ! surface run off  (kgH2O/m2/day)
       FLUXES(n,51) = mV%underflow       ! drainage from the bottom of the soil column (kgH2O/m2/day)
       FLUXES(n,52) = mV%water_grav_flow(1) ! drainage from the surface soil layer to 2nd  (kgH2O/m2/day)
       FLUXES(n,53) = mV%infiltrated(1)  ! top soil surface infiltration by rain (kgH2O/m2/day)
       FLUXES(n,56) = mV%infiltrated(2)  ! middle soil surface infiltration by rain (kgH2O/m2/day)
       FLUXES(n,57) = mV%infiltrated(3)  ! bottom soil surface infiltration by rain (kgH2O/m2/day)       
       FLUXES(n,54) = mV%uptake_fraction(1) ! transpiration fraction extracted from 1st rooting layer (the soil surface)
       FLUXES(n,55) = mV%uptake_fraction(2) ! transpiration fraction extracted from 2nd rooting layer (dynamic 2nd layer)       

      ! CUTTING 
      ! ------------------------------------------------------------------------------------------------------------- ! 

      ! LAI losses which have precise value of 1 are assumed to be a cutting
      if (met(8,n) == 1d0) then

          ! Estimate labile stored above ground assuming uniform distribution across biomass
          labile_ratio = POOLS(n+1,2) / (POOLS(n+1,2)+POOLS(n+1,3))
          call grass_cutting(labile_ratio,POOLS(n+1,1),POOLS(n+1,2),POOLS(n+1,3), & 
                             POOLS(n+1,4),POOLS(n+1,5),met(6,n), &
                             FLUXES(:,22),FLUXES(n,25),FLUXES(n,26),FLUXES(n,27), &
                             FLUXES(n,28),FLUXES(n,29),FLUXES(n,30), &
                             n,nodays,mV%days_per_step, & 
                             pars(28),pars(33), mV)

      else 

          ! GRAZING 
          ! ------------------------------------------------------------------------------------------------------------- ! 

          ! Estimate LAI change for the next time step. 
          ! This structure allows the grazing to include consumption of new growth.
          ! NOTE: met(8,n) sign change to make the lai losses now positive
          gsi_lai_reduction = -met(8,n) - ((POOLS(n,2)-POOLS(n+1,2)) / pars(15))
          ! ...if LAI has increased but the driver suggests a loss 
          !    or which is greater than already simulated consider grazing
          if (met(8,n) < 0d0 .and. gsi_lai_reduction > 0d0 .and. FLUXES(n,4)+FLUXES(n,7) > pars(34)) then
              ! Estimate labile stored above ground assuming uniform distribution across biomass
              labile_ratio = POOLS(n+1,2) / (POOLS(n+1,2)+POOLS(n+1,3))
              call grass_grazing(labile_ratio,POOLS(n+1,1),POOLS(n+1,2),POOLS(n+1,3),  & 
                                 POOLS(n+1,4),POOLS(n+1,5),met(6,n),gsi_lai_reduction, &
                                 FLUXES(n,19),FLUXES(n,20),FLUXES(n,21), &
                                 FLUXES(:,22),FLUXES(:,23),FLUXES(n,31),FLUXES(n,32), &
                                 FLUXES(n,33),FLUXES(n,34),FLUXES(n,35),FLUXES(n,36), &
                                 n,nodays,mV%days_per_step, & 
                                 pars(15),pars(27),pars(32),pars(34), mV)
          end if                      
      end if 

      ! NO FIRE MODEL HAS BEEN IMPLEMENTED AT THIS TIME - THIS IS SOMETHING WE SHOULD CONSIDER...
      
    end do ! nodays loop

  end subroutine CARBON_MODEL
  !
  !------------------------------------------------------------------
  !
  ! Subroutines below this line
  !
  !-----------------------------------------------------------------
  !
  !
  !------------------------------------------------------------------
  !
  subroutine acm_gpp_stage_1 (mV)
    use carbon_model_memory

    ! Estimate the light and temperature limited photosynthesis components.
    ! See acm_gpp_stage_2() for estimation of CO2 supply limitation and
    ! combination of light, temperature and CO2 co-limitation

    implicit none

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
    mV%rb_mol_1 = 1d0 / (mV%aerodynamic_conductance * mV%convert_ms1_mol_1 * gb_H2O_CO2 * &
              mV%leaf_canopy_wind_scaling)

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
    ! see Williams et al, (1997), Ecological Applications,7(3), 1997, pp. 882â€“894

    ! Daily canopy conductance dertmined through combination of aerodynamic and
    ! stomatal conductances. Both conductances are scaled to canopy aggregate.
    ! aerodynamic conductance already in units of molCO2.m-2.s-1 (see acm_gpp_stage_1).
    ! Stomatal conductance scaled from mmolH2O to molCO2.
    ! The ratio of H20:CO2 diffusion is 1.646259 (Jones appendix 2).
    !
    ! Combining in series the stomatal and boundary layer conductances
    ! to make canopy resistence (s/m2/molCO2)
    rc = 1d0 / (gs*gs_H2Ommol_CO2mol) + mV%rb_mol_1

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
            iWUE_upper = find_gs_iWUE(mV%potential_conductance, mV) !; iWUE_lower = find_gs_iWUE(min_gs, mV)
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
                                              find_gs_iWUE,mV,mV%minimum_conductance,mV%potential_conductance,tol_gs*mV%lai,mV%iWUE_step*0.10d0)

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
                                   water_retention_saxton_eqns, mV, x1 , x2 , 0.001d0, 0d0 )
    enddo

  end subroutine calculate_field_capacity  
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
    double precision ::    dT, & ! Canopy transmittance for long wave radiation
                        lwrad, & ! downward long wave radiation from sky (W.m-2)
        longwave_release_soil, & ! emission of long wave radiation from surfaces per m2
      longwave_release_canopy    ! assuming isothermal condition (W.m-2)

    ! estimate long wave radiation from atmosphere (W.m-2)
    lwrad = emiss_boltz * (mV%maxt+freeze-20d0) ** 4
    ! estimate isothermal long wave emission per unit area
    longwave_release_soil = emiss_boltz * (soil_temperature+freeze) ** 4
    ! estimate isothermal long wave emission per unit area
    longwave_release_canopy = emiss_boltz * (canopy_temperature+freeze) ** 4
    ! Canopy transmittance for thermal radiation
    dT = 1d0-exp(-mV%lai/mV%Vc*mV%mu_obar)

    !!!!!!!!!!
    ! Isothermal net long wave canopy and soil balance (W.m-2)
    !!!!!!!!!!

    ! Diffuse longwave absorbed by the canopy
    mV%canopy_lwrad_Wm2 = (lwrad*mV%Vc*dT) - (mV%Vc*dT*2d0*longwave_release_canopy) + (mV%Vc*dT*longwave_release_soil)
    ! Diffuse longwave absorbed by the soil
    mV%soil_lwrad_Wm2 = (lwrad*(1d0-(mV%Vc*dT))) + (mV%Vc*dT*longwave_release_canopy) - longwave_release_soil

  end subroutine calculate_longwave_isothermal
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_radiation_balance (mV)
    use carbon_model_memory

    implicit none

      type(model_working_variables) :: mV

    !
    ! Parametric approximation for 2-stream, multi-layer radiative transfer
    !

    !
    ! References
    !
    ! Sellers (1985) Canopy reflectance, photosynthesis and transpiration.
    !                International Journal of Remote Sensing, 6(8), 1335-1772, doi: 10.1080/01431168508948283
    ! Sellers et al., (1986) A simple biosphere model (SiB) for use in general circulation models.
    !                        Journal of Atmospheric Sciences, 43(6), 505-531
    ! Sellers et al., (1996) A revised land surface parameterisation (SiB2) for atmospheric GCMs. Part 1: Model formualtion.
    !                        Journal of Climate, 9(4), 676-705, doi: 10.1175/1520_0442(1996)009<0676:ARLSPF>2.0.CO;2

    !
    ! Description
    !
    ! The purpose of these equations is the approximation of the behaviuour of a 2-stream multiple canopy layer model
    ! accounting for solar angle, within canopy interception, reflectance and transmittance.
    !
    ! As a result the actual parametric inputs are few, but requires a large number of empiracal relationships to describe
    ! the non-linear dynamics of radiative transfer

    ! NOTE: that this code provides a daily timescale linear correction on
    ! isothermal longwave balance to net based on soil surface incident shortwave
    ! radiation. This correction is drawn from the standard ACM-GPP-ET-v1 approach

    ! Declare local variables
    double precision :: delta_iso

    ! Calculate the declination
    mV%declination = calculate_declination(mV%doy)
    ! Calculate cosine zenith angle
    call calculate_cosine_solar_zenith_angle(mV)
    ! Estimate shortwave radiation balance
    call calculate_shortwave_balance(mV)
    ! Estimate isothermal long wave radiation balance
    call calculate_longwave_isothermal(mV%meant,mV%meant, mV)
    ! Apply linear correction to soil surface isothermal->net longwave radiation
    ! balance based on absorbed shortwave radiation
    delta_iso = (soil_iso_to_net_coef_LAI * mV%lai) + &
                (soil_iso_to_net_coef_SW * (mV%soil_swrad_MJday * 1d6 * seconds_per_day_1)) + &
                 soil_iso_to_net_const
    ! In addition to the iso to net adjustment, SPA analysis shows that soil net never gets much below zero
    !soil_lwrad_Wm2 = max(-0.1d0,soil_lwrad_Wm2 + delta_iso)
    mV%soil_lwrad_Wm2 = mV%soil_lwrad_Wm2 + delta_iso
    ! Apply linear correction to canopy isothermal->net longwave radiation
    ! balance based on absorbed shortwave radiation
    delta_iso = (canopy_iso_to_net_coef_LAI * mV%lai) + &
                (canopy_iso_to_net_coef_SW * (mV%canopy_swrad_MJday * 1d6 * seconds_per_day_1)) + &
                canopy_iso_to_net_const
    mV%canopy_lwrad_Wm2 = mV%canopy_lwrad_Wm2 + delta_iso

  end subroutine calculate_radiation_balance
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_shortwave_balance (mV)
    use carbon_model_memory

    ! Subroutine estimates the canopy and soil absorbed shortwave radiation
    ! (MJ/m2/day). Radiation absorption is paritioned into NIR and PAR for
    ! canopy, and NIR + PAR for soil.
    ! Follows an implementation of the Sellers (1985) approximation

    implicit none

      type(model_working_variables) :: mV

    ! Declare local parameters
    double precision, dimension(no_wavelength), parameter :: &
                                   newsnow_reflectance = (/0.73d0,0.95d0/) ! NIR/PAR new snow reflectance fraction
    ! local variables
    double precision :: Gu, K, mu, S2, S3, fsnow, diffuse_fraction
    ! Local variables with different values per wavelength
    double precision, dimension(no_wavelength) :: &
                      as_mu, dd, ff, sigma, u1, u2, u3, &
                      S1, p1, p2, p3, p4, D1, D2, &
                      h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, &
                      Iup, Idown, soil_albedo, &
                      canopy_absorption_fraction_diffuse, &
                      canopy_absorption_fraction_direct, &
                      soil_absorption_fraction_diffuse, &
                      soil_absorption_fraction_direct, &
                      soil_nir_par_MJday, canopy_nir_par_MJday, &
                      swrad_direct, swrad_diffuse

    ! Estimate the diffuse fraction of shortwave radiation
    call calculate_diffuse_fraction(diffuse_fraction, mV)
    ! Estimate multiple use par and nir components (Note: units are MJ/m2/d)
    swrad_diffuse(1) = (1d0 - sw_par_fraction) * mV%swrad * diffuse_fraction      ! NIR
    swrad_diffuse(2) = sw_par_fraction * mV%swrad * diffuse_fraction              ! PAR
    swrad_direct(1) = (1d0 - sw_par_fraction) * mV%swrad * (1d0-diffuse_fraction) ! NIR
    swrad_direct(2) = sw_par_fraction * mV%swrad * (1d0-diffuse_fraction)         ! PAR

    ! Assign cosine_solar_zenith_angle to a local variable for easier readability
    mu = mV%cosine_solar_zenith_angle

    ! Relative projected area of leaf elements in direction of the cosine_solar_zenith_angle (mu).
    ! This variable is determined as the result of two empirical functions related to the leaf_distribution_deviance
    ! Note the notation used here Gu is varied, for clarity, from the actual used in Sellers (1985) which is G(mu).
    ! Similarly, the coefficients used in calculating Gu indicate "empty set".
    ! For clarity and visual similarity, I've used the capital letter O
    Gu = mV%O1 + (mV%O2 * mu)

    ! Extinction coefficient for direct radiation, related to Gu and the cosine zenith angle
    K = Gu / mu
    ! Single leaf scattering albedo within the canopy, varied by mu, leaf distribution and wavelength
    as_mu = ((mV%canopy_scattering * 0.5d0) * (Gu / (Gu+(mu*mV%O2)))) &
          * (1d0 - (mu*(mV%O1/(Gu+(mu*mV%O2)))*log((Gu+(mu*mV%O2)+(mu*mV%O1))/(mu*mV%O1)) ) )

    ! Upscatting coefficient for direct radiation, varied by mu and wavelength
    mV%beta0 = ((1d0+(mV%mu_obar*K)) / (mV%canopy_scattering*mV%mu_obar*K)) * as_mu

    ! Various terms, yet to have their specific functions determined
    ! Note that notations from Sellers (1985) have been given double letters if only single character was used,
    ! or written word for greek notation
    dd = mV%canopy_scattering * mV%mu_obar * K * mV%beta0
    ff = mV%canopy_scattering * mV%mu_obar * K * (1-mV%beta0)
    sigma = mV%cc**2 + mV%bb**2 + (mV%mu_obar*K)**2
    u1 = mv%bb - (mV%cc/mV%soil_reflectance) ; u2 = mV%bb - (mV%cc*mV%soil_reflectance) ; u3 = ff + (mV%cc*mV%soil_reflectance)
    S1 = exp(-mV%hh*mV%lai) ; S2 = exp(-K*mV%lai)

    ! Related to diffuse radiation
    p1 = mV%bb + (mV%mu_obar*mV%hh) ; p2 = mV%bb - (mV%mu_obar*mV%hh)
    ! Related to direct radiation
    p3 = mV%bb + (mV%mu_obar*K) ; p4 = mV%bb - (mV%mu_obar*K)
    !
    D1 = (p1 * (u1 - (mV%mu_obar*mV%hh)) * (1d0/S1)) - (p2*(u2+(mV%mu_obar*mV%hh))*S1)
    D2 = ((u2+(mV%mu_obar*mV%hh))*(1d0/S1)) - ((u2-(mV%mu_obar*mV%hh))*S1)

    !
    ! Direct radiation specific components
    !

    h1 = (-dd*p4) - (mV%cc*ff)
    h2 =  (1d0/D1) * ( ((dd-((h1/sigma)*p3))*((u1-(mV%mu_obar*mV%hh))*(1d0/S1))) &
                     - (p2*(dd-mV%cc-((h1/sigma)*(u1+(mV%mu_obar*K))))*S2) )
    h3 = (-1d0/D1) * ( ((dd-((h1/sigma)*p3))*(u1+(mV%mu_obar*mV%hh))*S1) &
                     - (p1*(dd-mV%cc-((h1/sigma)*(u1+(mV%mu_obar*K))))*S2) )
    h4 = (-dd*p3) - (mV%cc*ff) ! NOTE: "-" at the beginning is a correction identified in Sellers et al., (1996)
    h5 = (-1d0/D2) * ( ((h4/sigma)*(u2+(mV%mu_obar*mV%hh))*(1d0/S1)) &
                     + (u3-((h4/sigma)*(u2-(mV%mu_obar*K))*S2)) )
    h6 = (1d0/D2) * ( ((h4/sigma)*(u2-(mV%mu_obar*mV%hh))*S1) &
                     + (u3-((h4/sigma)*(u2-(mV%mu_obar*K))*S2)) )

    ! Fraction of direct radiation which leaves the canopy top as diffuse
    Iup = ((h1*exp(-K*mV%lai))/sigma) + (h2*exp(-mV%hh*mV%lai)) + (h3*exp(mV%hh*mV%lai))
    ! Fraction of direct radiation which leaves the canopy base as diffuse
    Idown = ((h4*exp(-K*mV%lai))/sigma) + (h5*exp(-mV%hh*mV%lai)) + (h6*exp(mV%hh*mV%lai))

    ! Update soil reflectance based on snow cover
    if (mV%snow_storage > 0d0) then
        fsnow = 1d0 - exp( - mV%snow_storage * 0.1d0 )  ! fraction of snow cover on the ground
        soil_albedo = ((1d0 - fsnow) * mV%soil_reflectance) + (fsnow * newsnow_reflectance) ! NIR & PAR
    else
        ! Estimate the combined soil and surface layer reflectances, assume direct and diffuse reflectances are isotropic
        ! NOTE; this variable is also used in the diffuse calculation too
        !soil_albedo = (soil_surface_reflectance * Vg) + ((1d0-Vg) * soil_reflectance)
        soil_albedo = mV%soil_reflectance
    endif

    S3 = exp(-K*mV%lai/mV%Vc)
    ! Fraction of direct radiation absorbed by the canopy
    canopy_absorption_fraction_direct = mV%Vc * (1d0 - Iup - (Idown*(1d0-soil_albedo)) &
                                             - (S3*(1d0-soil_albedo)))
    ! Fraction of direct radiation absorbed by the soil
    soil_absorption_fraction_direct = ((1d0-mV%Vc)*(1d0-soil_albedo)) &
                                    + (mV%Vc*((Idown*(1d0-soil_albedo)) + (S3*(1d0-soil_albedo))))

    !
    ! Diffuse radiation specific components
    !

    h7  = (mV%cc/D1) * (u1-(mV%mu_obar*mV%hh)) * (1d0/S1)
    h8  = (-mV%cc/D1) * (u1+(mV%mu_obar*mV%hh)) * S1
    h9  = (1d0/D2) * (u2+(mV%mu_obar*mV%hh)) * (1d0/S1)
    h10 = (-1d0/D2) * (u2-(mV%mu_obar*mV%hh)) * S1

    ! Fraction of diffuse radiation which leaves the canopy top as diffuse
    Iup = (h7*exp(-mV%hh*mV%lai)) + (h8*exp(mV%hh*mV%lai))
    ! Fraction of diffuse radiation which leaves the canopy base as diffuse
    Idown = (h9*exp(-mV%hh*mV%lai)) + (h10*exp(mV%hh*mV%lai))

    ! Fraction of diffuse radiation absorbed by the canopy
    canopy_absorption_fraction_diffuse = mV%Vc * (1d0 - Iup - (Idown*(1d0-soil_albedo)))
    ! Fraction of diffuse radiation absorbed by the soil
    soil_absorption_fraction_diffuse = ((1d0-mV%Vc)*(1d0-soil_albedo)) + (mV%Vc*((Idown*(1d0-soil_albedo))))

    !
    ! Combine direct and diffuse, convert into actual units of energy (MJ/m2/d)
    !

    ! Determine the combined direct and diffuse absorptions across wavelength
    soil_nir_par_MJday = (soil_absorption_fraction_diffuse * swrad_diffuse) &
                       + (soil_absorption_fraction_direct * swrad_direct)
    canopy_nir_par_MJday = (canopy_absorption_fraction_diffuse * swrad_diffuse) &
                       + (canopy_absorption_fraction_direct * swrad_direct)
    ! Assign to specific output variables for use elsewhere in the model
    mV%canopy_par_MJday = canopy_nir_par_MJday(2)
    mV%soil_par_MJday = soil_nir_par_MJday(2)
    ! combine totals for use is soil evaporation
    mV%soil_swrad_MJday = sum(soil_nir_par_MJday)
    ! Combine to estimate total shortwave canopy absorbed radiation
    mV%canopy_swrad_MJday = sum(canopy_nir_par_MJday)

    ! Estimate the integral of light interception for use as a leaf to canopy
    ! scaler for photosynthesis, transpiration, and gs
    ! Based on the canopy scaling of photosynthetic capacity due to light from 
    ! Sellers et al., (1992), Remote Sensing Environment, 42(3), 187-216.
    mV%leaf_canopy_light_scaling = (1d0-S2) / K

    ! check energy balance
!    balance = swrad - canopy_par_MJday - canopy_nir_MJday - refl_par_MJday - refl_nir_MJday - soil_swrad_MJday
!    if ((balance - swrad) / swrad > 0.01) then
!        print*,"SW residual frac = ",(balance - swrad) / swrad,"SW residual = ",balance,"SW in = ",swrad
!   endif

  end subroutine calculate_shortwave_balance
  !
  !-----------------------------------------------------------------
  !
  !subroutine calculate_radiation_commons(lat,rad_pars)
  subroutine calculate_radiation_commons(lat, mV)
    use carbon_model_memory

    implicit none

      type(model_working_variables) :: mV

    ! Description

    ! Declare arguments
    double precision, intent(in) :: lat!, & ! site latitude in degrees
                               !rad_pars(6) !

    ! Calculate some common variables and place into memory
    mV%latitude = lat
    mV%latitude_radians = lat * deg_to_rad
    mV%sin_latitude_radians = sin(mV%latitude_radians)
    mV%cos_latitude_radians = cos(mV%latitude_radians)

    ! Load canopy optical properties to their module variables
!    canopy_reflectance(1)   = rad_pars(1) ! canopy_nir_reflectance
!    canopy_reflectance(2)   = rad_pars(2) ! canopy_par_reflectance
!    canopy_transmittance(1) = rad_pars(3) ! canopy_nir_transmittance
!    canopy_transmittance(2) = rad_pars(4) ! canopy_par_transmittance
!    soil_reflectance(1)     = rad_pars(5) ! soil_nir_reflectance
!    soil_reflectance(2)     = rad_pars(6) ! soil_par_reflectance

    ! Canopy scattering of incident light, varied by wavelength
    ! NOTE: if we want to put snow fall on the canopy in the model, then this
    ! line will need to move into shortwave_balance() or calculate_radiation_balance()
    ! to be recalculated with each update of the canopy reflectance and transmittances.
    mV%canopy_scattering = mV%canopy_reflectance + mv%canopy_transmittance

    ! Two empirical functions related to the leaf_distribution_deviance
    ! used in the calcuation of several variable varying by day of year
    ! The coefficient notation used in Sellers (1985) indicate "empty set".
    ! For clarity and visual similarity, I've used the capital letter O
    mV%O1 = 0.5d0-(0.633d0*leaf_distribution_deviance)-(0.33d0*leaf_distribution_deviance**2d0)
    mV%O2 = 0.877d0*(1d0-(2d0*mV%O1))

    ! The average inverse diffuse optical depth per unit leaf area.
    ! mu_obar indicates mu with over bar
    ! See notes or references for the original integral equation
    mV%mu_obar = (1d0/mV%O2) * ( 1d0-(mV%O1/mV%O2)*log((mV%O1+mV%O2)/mV%O1) )
    ! Mean leaf angle deviation from the horizontal (radians)
    ! SPA default assumption is 30 degrees, where (pi/180) is the conversion to radians
    mV%leaf_angle = 0.5235988d0 !30d0 * (pi/180d0)
    ! Analytical correction for leaf angle (radians) on light scattering within the canopy
    ! Notation is the verbal description of that used in Sellers (1985)
    mV%cos2theta = (1d0 + cos(2d0*mV%leaf_angle)) * 0.5d0

    ! Upward scattering as diffuse radiation, a function of canopy_transmittance, canopy_reflectance and leaf angle.
    ! Notation in Sellers (1985) is "c"
    mV%cc = 0.5 * (mV%canopy_reflectance + mV%canopy_transmittance + (mV%canopy_reflectance - mV%canopy_transmittance) * mV%cos2theta)
    ! The upward scattering fraction / coefficient for diffuse radiation.
    ! NOTE there is a direct radiation equivalent found below, noted as beta0
    mV%beta = mV%cc / mV%canopy_scattering
    ! Downward scatting of diffuse radiation
    ! Note that notations from Sellers (1985) have been given double letters if only single character was used,
    ! or written word for greek notation
    mV%bb = (1d0-(1d0-mV%beta)*mV%canopy_scattering)
    ! Extinction coefficient for diffuse radiation, related to absorption normalised by mu_obar
    ! Note that notations from Sellers (1985) have been given double letters if only single character was used,
    ! or written word for greek notation
    mV%hh = (((mV%bb**2d0 - mV%cc**2d0))**(0.5d0)) / mV%mu_obar

    ! Define the vegetated and covered soil (i.e. by litter) fractions
    ! Vc, could also be considered a canopy clumping factor
    mV%Vc = 0.75d0 ; mV%Vg = 0d0

    ! Return back to user
    return

  end subroutine calculate_radiation_commons  
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_cosine_solar_zenith_angle (mV)
    use carbon_model_memory

    implicit none

      type(model_working_variables) :: mV

    ! Calculate some common variables needed for the Sellers (1985)
    ! parameteric appoximation

    ! Declare local parameters

    ! Solar angle for the hour of the day, assumed in daily model to be 10 am (or 2 pm)
    ! therefore, 2 = 14 - 12, where 14 = timestep of day and 12 = half the number
    ! hour-angle = ( hr - 12 ) * 15. * ( pi / 180.d0 )
    ! of time steps per day.
    double precision, parameter :: cos_hour_angle = cos(0.5235988d0)

    ! Declare local variables

    ! Now calculate the solar-zenith-angle, as per..
    !  sin(latitude)sin(declination) + cos(latitude)cos(declination)cos(hour_angle)
    ! where hour-angle = ( hr - 12 ) * 15. * ( pi / 180.d0 )
    ! and from that the solar flux.
    ! eqn 2.15 in Hartmann's Global Physical Climatology...
    ! Minimum allowed value constraint from SIB3 implementation of Sellers (1985)
    mV%cosine_solar_zenith_angle = max(0.01747d0, mV%sin_latitude_radians * sin(mV%declination) + &
                                               mV%cos_latitude_radians * cos(mV%declination) * &
                                               cos_hour_angle)

  end subroutine calculate_cosine_solar_zenith_angle
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_daylength (mV)
    use carbon_model_memory

    ! Subroutine uses day of year and latitude (-90 / 90 degrees) as inputs,
    ! combined with trigonomic functions to calculate day length in hours and seconds

    implicit none

      type(model_working_variables) :: mV

    ! local variables
    double precision :: dec, sinld, cosld, aob

    !
    ! Estimate solar geometry variables needed
    !

    ! day length is estimated as the ratio of sin and cos of the product of declination an latitude in radiation
    sinld = mV%sin_latitude_radians * sin( mV%declination )
    cosld = mV%cos_latitude_radians * cos( mV%declination )
    aob = max(-1d0,min(1d0,sinld / cosld))

    ! estimate day length in hours and seconds and upload to module variables
    mV%dayl_hours = 12d0 * ( 1d0 + 2d0 * asin( aob ) * pi_1 )
    mV%dayl_seconds = mV%dayl_hours * seconds_per_hour

    ! return to user
    return

  end subroutine calculate_daylength
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_diffuse_fraction(diffuse_fraction, mV)
    use carbon_model_memory

    implicit none

      type(model_working_variables) :: mV

    ! Description estimate ratio of actual to extra-solar radiation
    ! for day length period following Erbs et al., (1982)

    ! Declare arguments
    double precision, intent(out) :: diffuse_fraction

    ! Declare local parameters
    double precision, parameter :: So = 117.504d0 ! solar constant (1360 Wm-2)
                                                  ! unit scaled to MJ/m2/day
    ! Declare local variables
    double precision :: clear_day_swrad, Kt

    ! Estimate sunset solar angle
    mV%sunset_solar_angle = acos(-tan(mV%latitude_radians)*tan(mV%declination))
    ! Estimate clear day light solar radiation (MJ/m2/d)
    clear_day_swrad = So * pi_1 * (1d0+0.033d0*cos(two_pi*mV%doy/365d0)) &
                    * (mV%cos_latitude_radians*cos(mV%declination)*sin(mV%sunset_solar_angle) &
                      + mV%sin_latitude_radians*sin(mV%declination))
    ! Estimate the ratio of actual to clear day radiation (MJ/m2/d over MJ/m2/d)
    Kt = min(1d0,mV%swrad / clear_day_swrad)

    ! Calculate diffuse ratio
    if (mV%sunset_solar_angle < 1.4208d0) then
        if (Kt < 0.715d0) then
            diffuse_fraction = 1d0 - 0.2727d0*Kt + 2.4495d0*Kt**2 - 11.9514d0*Kt**3 + 9.3879d0*Kt**4
        else
            diffuse_fraction = 0.143d0
        end if
    else
        if (Kt < 0.722d0) then
            diffuse_fraction = 1d0 + 0.28832d0*Kt - 2.5557d0*Kt**2 + 0.8448d0*Kt**3
        else
            diffuse_fraction = 0.175d0
        end if
    end if

    ! Sanity check
    diffuse_fraction = min(1d0, diffuse_fraction)

    ! return
    return

  end subroutine calculate_diffuse_fraction
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
    ! calculate the plant hydraulic resistance component.
    transpiration_resistance = canopy_height / (gplant * max(min_lai,mV%lai))

    !!!!!!!!!!!
    ! calculate current steps soil hydraulic conductivity
    !!!!!!!!!!!

    ! seperately calculate the soil conductivity as this applies to each layer
    do i = 1, nos_soil_layers
       call calculate_soil_conductivity(i,mV%soil_waterfrac(i),mV%soil_conductivity(i), mV)
    end do ! soil layers

    !!!!!!!!!!!
    ! Calculate hydraulic properties and each rooted layer
    !!!!!!!!!!!

    ! calculate and accumulate steady state water flux in mmol.m-2.s-1
    ! NOTE: Depth correction already accounted for in soil resistance
    ! calculations and this is the maximum potential rate of transpiration
    ! assuming saturated soil and leaves at their minimum water potential.
    mV%demand = max(0d0, (mV%SWP(1:nos_root_layers) - (head*canopy_height)) - minlwp )
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

    ! if freezing then assume soil surface is frozen, therefore no water flux
    if (mV%meant < 1d0) then
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

    depth_change = (top_soil_depth+min_layer) ; water_change = 0
    ! if roots extent down into the bucket
    if (mV%root_reach > depth_change) then

        !!!!!!!!!!
        ! Soil profile is within the bucket layer (layer 3)
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

    ! Return soil water back to its initial field capacity value
    mV%soil_waterfrac(1:nos_soil_layers) = mV%field_capacity(nos_soil_layers)

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
  subroutine gravitational_drainage(time_period_days, mV)
    use carbon_model_memory

    ! Integrator for soil gravitational drainage.
    ! Due to the longer time steps undertake by ACM / DALEC and the fact that
    ! drainage is a concurrent processes we assume that drainage occurs at
    ! the bottom of the column first creating space into which water can drain
    ! from the top down. Therefore we drain from the bottom first and then the top.
    ! NOTE: Assumes that any previous water movement due to infiltration and evaporation
    !       has already been updated in soil mass balance

    implicit none

      type(model_working_variables) :: mV

    ! arguments
    integer, intent(in) :: time_period_days

    ! local variables..
    integer ::  soil_layer
    integer :: t
    double precision, dimension(nos_soil_layers) :: dx, & ! range between the start and end points of the integration
                                               halfway, & ! half way point between start and end point of integration
                                                liquid, & ! liquid water in local soil layer (m3/m3)
                                         avail_to_flow, & ! liquid content above field capacity (m3/m3)
                                               iceprop, & ! fraction of soil layer which is ice
                                          pot_drainage    ! estimats of time step potential drainage rate (m/s)
    double precision  :: tmp1,tmp2,tmp3 &
                                 ,unsat & ! unsaturated pore space in soil_layer below the current (m3/m3)
                                ,change   ! absolute volume of water drainage in current layer (m3/day)

    ! calculate soil ice proportion; at the moment
    ! assume everything liquid
    iceprop = 0d0

    ! except the surface layer in the mean daily temperature is < 0oC
    if (mV%meant < 1d0) iceprop(1) = 1d0

    ! zero water fluxes
    mV%waterchange = 0d0

    ! underflow is tracked in kgH2O/m2/day but estimated here in MgH2O/m2/day
    ! therefore we must convert
    mV%underflow = mV%underflow * 1d-3
    mV%water_grav_flow = mV%water_grav_flow * 1d-3

    ! estimate potential drainage rate for the current time period
    liquid = mV%soil_waterfrac(1:nos_soil_layers) * ( 1d0 - iceprop(1:nos_soil_layers) )
    ! estimate how much liquid is available to flow
    avail_to_flow = liquid - mV%field_capacity(1:nos_soil_layers)

    ! trapezium rule scaler and the half-way point between current and field capacity
    dx = avail_to_flow*0.5d0 ; halfway = liquid - dx
    do t = 1, nos_soil_layers
       if (avail_to_flow(t) > 0d0) then
           ! Trapezium rule for approximating integral of drainage rate
           call calculate_soil_conductivity(t,liquid(t),tmp1, mV)         ! This is the maximum rate
           call calculate_soil_conductivity(t,mV%field_capacity(t),tmp2, mV) ! This is the final rate
           call calculate_soil_conductivity(t,halfway(t),tmp3, mV)        ! This is half way between
           pot_drainage(t) = 0.5d0 * dx(t) * ((tmp1 + tmp2) + 2d0 * tmp3)
       else
           ! We are at field capacity currently even after rainfall has been infiltrated.
           ! Assume that the potential drainage rate is that at field capacity
           call calculate_soil_conductivity(t,mV%field_capacity(t),pot_drainage(t), mV)
       endif ! water above field capacity to flow?
    end do ! soil layers
    ! Scale potential drainage from per second to per day
    pot_drainage = pot_drainage * seconds_per_day

    ! Integrate drainage over each day until time period has been reached or
    ! each soil layer has reached field capacity
    t = 1
    do while (t < (time_period_days+1) .and. maxval(mV%soil_waterfrac - mV%field_capacity) > vsmall)

       ! Estimate liquid content and how much is available to flow / drain
       avail_to_flow = ( mV%soil_waterfrac(1:nos_soil_layers) * (1d0 - iceprop(1:nos_soil_layers)) ) &
                     - mV%field_capacity(1:nos_soil_layers)

       ! ...then from the top down
       do soil_layer = 1, nos_soil_layers

          ! initial conditions; i.e. is there liquid water and more water than
          ! layer can hold
          if (avail_to_flow(soil_layer) > 0d0 .and. mV%soil_waterfrac(soil_layer+1) < mV%porosity(soil_layer+1)) then

              ! Unsaturated volume of layer below (m3 m-2)
              unsat = ( mV%porosity(soil_layer+1) - mV%soil_waterfrac(soil_layer+1) ) &
                    * mV%layer_thickness(soil_layer+1) / mV%layer_thickness(soil_layer)
              ! Restrict potential rate calculate above for the available water
              ! and available space in the layer below.
              ! NOTE: * layer_thickness(soil_layer) converts units from m3/m2 -> (m3)
              change = min(unsat,min(pot_drainage(soil_layer),avail_to_flow(soil_layer))) * mV%layer_thickness(soil_layer)
              ! update soil layer below with drained liquid
              mV%waterchange( soil_layer + 1 ) = mV%waterchange( soil_layer + 1 ) + change
              mV%waterchange( soil_layer     ) = mV%waterchange( soil_layer     ) - change
              ! Also track only the positive flows from one layer to another (MgH2O/m2/day)
              mV%water_grav_flow(soil_layer) = mV%water_grav_flow(soil_layer) + change

          end if ! some liquid water and drainage possible

       end do ! soil layers

       ! update soil water profile
       mV%soil_waterfrac(1:nos_soil_layers) = mV%soil_waterfrac(1:nos_soil_layers) &
                                         + (mV%waterchange(1:nos_soil_layers)/mV%layer_thickness(1:nos_soil_layers))
       ! estimate drainage from bottom of soil column (MgH2O/m2/day)
       ! NOTES: that underflow is reset outside of the daily soil loop
       mV%underflow = mV%underflow + mV%waterchange(nos_soil_layers+1)

       ! Reset now we have moves that liquid
       mV%waterchange = 0d0
       ! integerate through time period
       t = t + 1

    end do ! while condition

    ! convert underflow and water_grav_flow from MgH2O/m2/day -> kgH2O/m2/day
    mV%underflow = mV%underflow * 1d3
    mV%water_grav_flow = mV%water_grav_flow * 1d3

  end subroutine gravitational_drainage
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
    ! calculate field capacity (m3/m-3)
    call calculate_field_capacity(mV)

    ! final sanity check for porosity
    do i = 1, nos_soil_layers+1
       if (mV%porosity(i) < (mV%field_capacity(i)+0.05d0)) mV%porosity(i) = mV%field_capacity(i) + 0.05d0
    end do

  end subroutine initialise_soils
  !
  !---------------------------------------------------------------------
  !
  subroutine update_soil_initial_conditions (mV)
    use carbon_model_memory

    !
    ! Subroutine calculate the soil layers field capacities and sets the initial
    ! soil water potential set to field capacity
    !

    implicit none

      type(model_working_variables) :: mV

    ! local variables
    integer :: i

    ! Load initial soil water fraction to the dynamic layers
    mV%soil_waterfrac(1:nos_soil_layers) = mV%field_capacity(1:nos_soil_layers)
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
  !----------------------------------------------------------------------
  !
  subroutine initialise_gsi(mean_days_per_step,gsi_lag_steps,gsi_lag_days)
    use carbon_model_memory

    ! Subroutine tidys away the calculation of the 
    ! number of days for the growing season index lag period

    implicit none

    ! Arguments
    integer, intent(in) :: gsi_lag_steps
    double precision, intent(in) :: mean_days_per_step
    double precision, dimension(gsi_lag_steps), intent(out) :: gsi_lag_days
    ! Local variables
    integer :: f          

    ! Determine the number of days equivalent for the lag period
    do f = 1, gsi_lag_steps
       gsi_lag_days(f) = dble(f) * mean_days_per_step
    end do

    ! Return
    return

  end subroutine initialise_gsi
  !
  !------------------------------------------------------------------
  !
  subroutine grass_cutting(labile_ratio,labile,foliage,roots,litter,som,doy,  &
                           harvest,HARVESTextracted_labile,     & 
                           HARVESTextracted_foliage,HARVESTextracted_roots,   &
                           HARVESTlitter_labile,HARVESTlitter_foliage,        &
                           HARVESTlitter_roots,                               &
                           timestep,nodays,step_length,                       &
                           cutting_threshold,post_cutting_labile_loss, mV)
    use carbon_model_memory
  
    !! Determine whether cutting has occured and the resulting impacts on the C-cycle
    
    implicit none

      type(model_working_variables) :: mV
    
    ! Arguments
    integer, intent(in) :: timestep, nodays
    double precision, intent(in) :: doy, &
                            step_length, &
                           labile_ratio, &
                      cutting_threshold, &
               post_cutting_labile_loss
    
    double precision, intent(inout) :: labile, &
                                      foliage, &
                                        roots, &
                                       litter, & 
                                          som, &
                      HARVESTextracted_labile, &
                     HARVESTextracted_foliage, &
                       HARVESTextracted_roots, &
                         HARVESTlitter_labile, &
                        HARVESTlitter_foliage, &
                          HARVESTlitter_roots   
                          
    double precision, dimension(nodays), intent(inout) :: harvest
                           
    ! Local variables
    double precision :: labile_loss    & 
                       ,foliar_loss    &
                       ,roots_loss     &
                       ,labile_residue &
                       ,foliar_residue &
                       ,roots_residue  
           
    ! Determine whether cutting is plausible
    ! 1) Labile+leaf C > cutting threshold 
    ! 2) LAI > 3 (note replaced here with day of year constraints?)
    ! 3) & LAI reduction = -1 & no cut in past month 
    ! TLS: Question, conditions here are very temperate centric, can these be modifed?
    ! TLS: The timestepping of the cutting assumptions need to be dynamics in code to timestep
    if ( (labile+foliage) >= cutting_threshold   & 
          .and. doy >= 91d0 .and. doy <= 304d0 & 
          ! .and. LAI(n) >= 3 & 
          .and. sum(harvest(max(1,timestep-mV%four_week_lag):timestep)) == 0d0 ) then
                  
        ! direct C losses
        labile_loss  = labile * 0.95d0 * post_cutting_labile_loss
        foliar_loss  = foliage * 0.95d0 ! 95% of leaves lost after cutting probably 99% lost in reality 
        roots_loss   = 0d0 ! POOLS(n+1,3) * roots_frac_death ! allocation to roots will be reduced due to reduced LAI 

        ! fraction of harvest wasted 
        labile_residue = labile_loss * mV%labile_frac_res
        foliar_residue = foliar_loss * mV%foliage_frac_res
        roots_residue  = roots_loss  * mV%roots_frac_res

        ! if havest yields > 1500 kg.DM.ha-1 proceed with cut
        ! Note converted to gC/m2 equivalent assuming 47.5 % C content
        ! yields 71.25 gC/m2
        if ( ( (foliar_loss-foliar_residue)+ &
               (labile_loss-labile_residue)+ &
               (roots_loss -roots_residue ) ) >= 71.25d0 ) then
                      
            ! Assign to output the biomass extracted due to cutting
            HARVESTextracted_labile  = (labile_loss-labile_residue) / step_length
            HARVESTextracted_foliage = (foliar_loss-foliar_residue) / step_length
            HARVESTextracted_roots   = (roots_loss -roots_residue) / step_length
            ! Assign the output the biomass entering litter due to cutting
            HARVESTlitter_labile  = labile_residue / step_length
            HARVESTlitter_foliage = foliar_residue / step_length
            HARVESTlitter_roots   = roots_residue / step_length

            ! Combine to total extracted based on cutting
            harvest(timestep) = HARVESTextracted_labile  + &
                                HARVESTextracted_foliage + & 
                                HARVESTextracted_roots 

            ! update pools 
            labile  = max(0d0,labile-labile_loss)
            foliage = max(0d0,foliage-foliar_loss)
            roots   = max(0d0,roots-roots_loss)
            litter  = max(0d0,litter + (labile_residue+foliar_residue+roots_residue))
            som     = max(0d0,som)

        endif ! Determine whether the cut is actually plausible

    endif ! end cutting process   

  end subroutine grass_cutting
  !
  !------------------------------------------------------------------
  !
  subroutine grass_grazing(labile_ratio,labile,foliage,roots,litter,som,doy,        &
                           lai_reduction,animal_manure_to_soil,animal_respiration,  &
                           animal_methane,harvest,grazing,GRAZINGextracted_labile,  & 
                           GRAZINGextracted_foliage,GRAZINGextracted_roots,         &
                           GRAZINGlitter_labile,GRAZINGlitter_foliage,              &
                           GRAZINGlitter_roots,                                     &
                           timestep,nodays,step_length,                             &
                           lca,grazing_threshold,post_grazing_labile_loss,          &
                           min_grazing_removal_threshold, mV)
    use carbon_model_memory
  
    !! Determine whether grazing has occured and the resulting impacts on the C-cycle
    !! Livestock units (LU) per hectare following the assumptions that:
    !! (1) one cattle is 1â€‰LU and one sheep is 0.11â€‰LU, (2) 1â€‰LU weighs 650â€‰kg, 
    !! (3) an animal demands â‰ˆ 2.5â€‰% (p31) of its weight in the form grass dry matter (DM) when grazing, 
    !! and (4) 47.5â€‰% of DM consists of C (VertÃ¨s et al., 2018). 
    !! This default assumption equates to 1 LSU / ha / day requiring ~0.77 gC/m2/day.
    
    implicit none

      type(model_working_variables) :: mV
    
    ! Arguments
    integer, intent(in) :: timestep, nodays
    double precision, intent(in) :: doy, &
                            step_length, &
                           labile_ratio, &
                          lai_reduction, &
                                    lca, &
                      grazing_threshold, &
          min_grazing_removal_threshold, &
               post_grazing_labile_loss
    
    double precision, intent(inout) :: labile, &
                                      foliage, &
                                        roots, &
                                       litter, & 
                                          som, & 
                      GRAZINGextracted_labile, &
                     GRAZINGextracted_foliage, &
                       GRAZINGextracted_roots, &
                         GRAZINGlitter_labile, &
                        GRAZINGlitter_foliage, &
                          GRAZINGlitter_roots, &
                        animal_manure_to_soil, &
                           animal_respiration, &
                               animal_methane
                       
    double precision, dimension(nodays), intent(inout) :: harvest, &
                                                          grazing

    ! Local variables
    double precision :: fraction_loss  &
                       ,labile_loss    & 
                       ,foliar_loss    &
                       ,roots_loss     &
                       ,labile_residue &
                       ,foliar_residue &
                       ,roots_residue 
    double precision, parameter :: &
                 grazing_to_manure = 0.32d0, & ! Fraction of grazed C that ends up as manure
     grazing_to_animal_respiration = 0.54d0, & ! Fraction of grazed C that is respired by animals
                grazing_to_methane = 0.04d0    ! Fraction of grazed C that ends as methane via animal digestion

    ! Determine whether grazing is plausible
    ! 1) An LAI reduction is specified
    ! 2) Labile+leaf C > grazing threshold 
    ! 3) No cutting in the last 2 weeks
!    if ((labile*labile_ratio)+foliage >= grazing_threshold .and. &
!        sum(harvest(max(1,timestep-two_week_lag):timestep)) == 0d0) then
    if ((labile*labile_ratio)+foliage >= grazing_threshold .and. foliage > 0d0 .and. &
        sum(harvest(max(1,timestep-mV%two_week_lag):timestep)) == 0d0) then

        !   
        ! direct C losses
        !

        ! Determine the fraction of foliage being removed by the lai reduction
        fraction_loss = (lai_reduction * lca) / foliage
        ! Assume that a fraction of labile is also storged proportionally within the leaf / above ground
        labile_loss  = labile * labile_ratio * fraction_loss * post_grazing_labile_loss
        ! Apply fractional loss to foliage
        foliar_loss  = foliage * fraction_loss
        roots_loss   = 0d0 ! POOLS(n+1,3) * roots_frac_death
!        labile_loss  = labile * post_grazing_labile_loss
!        foliar_loss  = max(0d0,(lai_reduction * lca) - labile_loss)  
!        roots_loss   = 0d0 ! POOLS(n+1,3) * roots_frac_death

        ! fraction of grazing lost as residue (messy eaters)
        labile_residue = labile_loss * mV%labile_frac_res
        foliar_residue = foliar_loss * mV%foliage_frac_res
        roots_residue  = roots_loss  * mV%roots_frac_res

        ! extracted C via grazing: if remaining AGB > pre-grazing limit DM & grazed biomass > pars(34) gC/m2/step
        if ( (((labile*labile_ratio)+foliage)-(foliar_loss+labile_loss)) >= grazing_threshold .and. & 
               (foliar_loss+labile_loss) >= min_grazing_removal_threshold*step_length ) then

            ! Assign to output the biomass extracted due to cutting
            GRAZINGextracted_labile  = (labile_loss-labile_residue) / step_length
            GRAZINGextracted_foliage = (foliar_loss-foliar_residue) / step_length
            GRAZINGextracted_roots   = (roots_loss -roots_residue)  / step_length
            ! Assign the output the biomass entering litter due to cutting
            GRAZINGlitter_labile  = labile_residue / step_length
            GRAZINGlitter_foliage = foliar_residue / step_length
            GRAZINGlitter_roots   = roots_residue  / step_length

            ! Combine to total extracted based on grazing
            grazing(timestep) = GRAZINGextracted_labile  + &
                                GRAZINGextracted_foliage + &
                                GRAZINGextracted_roots
               
            ! Constants used for animal C fluxes from Vertes.et.al.2019 (10.1016/B978-0-12-811050-8.00002-9)
            ! TLS: Meaning of the constants needs to be defined.

            ! animal manure-C production (gC/m2)
            animal_manure_to_soil = (grazing(timestep) * grazing_to_manure) / step_length
            ! animal respiration CO2-C (gC/m2)
            animal_respiration = (grazing(timestep) * grazing_to_animal_respiration) / step_length
            ! animal CH4-C (gC/m2)
            animal_methane = (grazing(timestep) * grazing_to_methane) / step_length                

            ! update pools 
            labile  = max(0d0,labile-labile_loss)
            foliage = max(0d0,foliage-foliar_loss)
            roots   = max(0d0,roots-roots_loss)
            litter  = max(0d0,litter+((animal_manure_to_soil+GRAZINGlitter_labile + &
                                       GRAZINGlitter_foliage+GRAZINGlitter_roots) * step_length))
            som     = max(0d0,som)

        else if (foliar_loss+labile_loss >= min_grazing_removal_threshold*step_length) then

            ! If losses are > than min grazing removal threshold, then it is implicit that 
            ! the remaining labile+foliage was below the critical above ground threshold.
            ! Therefore, we will attempt to determine whether removals required to reach the minimum 
            ! grazing threshold (which is also the minimum which must remain afterwards).

            ! 
            ! Direct C losses
            !

            ! Determine the adjustment to the existing loss terms, to ensure mass balance with the grazing_threshold value
            fraction_loss = (((labile*labile_ratio) + foliage) - grazing_threshold) / (labile_loss + foliar_loss)
            ! Adjust the existing labile loss based on adjustment fraction
            labile_loss  = labile_loss * fraction_loss
            ! Apply fractional loss to foliage
            foliar_loss  = foliar_loss * fraction_loss
            roots_loss   = 0d0 ! POOLS(n+1,3) * roots_frac_death

            !! Determine losses which result in remaining biomass to be the grazing_threshold value
            !labile_loss  = labile * post_grazing_labile_loss
            !foliar_loss  = foliage - (grazing_threshold - (labile - labile_loss))
            !!foliar_loss  = foliage - (grazing_threshold - labile_loss)
            !roots_loss   = 0d0 ! roots * roots_frac_death

            ! fraction of harvest wasted 
            labile_residue = labile_loss * mV%labile_frac_res
            foliar_residue = foliar_loss * mV%foliage_frac_res
            roots_residue = roots_loss * mV%roots_frac_res

            ! proceed if simulating this grazing will remove > ~0.5 gCm-2 from AGB
! TLS: Commented out due to greater agreement with the grazing density information at North Wyke
! That is not to say there isn't some utility in a 2nd chance component in the future when we have more 
! diverse datasets available.
!            if ((foliar_loss+labile_loss) >= min_grazing_removal_threshold*step_length) then
!
!                ! Assign to output the biomass extracted due to cutting
!                GRAZINGextracted_labile  = (labile_loss-labile_residue) / step_length
!                GRAZINGextracted_foliage = (foliar_loss-foliar_residue) / step_length
!                GRAZINGextracted_roots   = (roots_loss -roots_residue) / step_length
!                ! Assign the output the biomass entering litter due to cutting
!                GRAZINGlitter_labile  = labile_residue / step_length
!                GRAZINGlitter_foliage = foliar_residue / step_length
!                GRAZINGlitter_roots   = roots_residue / step_length
!
!                ! Combine to total extracted based on grazing
!                grazing(timestep) = GRAZINGextracted_labile + &
!                                    GRAZINGextracted_foliage + &
!                                    GRAZINGextracted_roots
!
!                ! animal manure-C production (gC/m2)
!                animal_manure_to_soil = (grazing(timestep) * grazing_to_manure) / step_length
!                ! animal respiration CO2-C (gC/m2)
!                animal_respiration = (grazing(timestep) * grazing_to_animal_respiration) / step_length
!                ! animal CH4-C (gC/m2)
!                animal_methane = (grazing(timestep) * grazing_to_methane) / step_length                
!
!                ! update pools 
!                labile  = max(0d0,labile-labile_loss)
!                foliage = max(0d0,foliage-foliar_loss)
!                roots   = max(0d0,roots-roots_loss)
!                litter  = max(0d0,litter+((animal_manure_to_soil+GRAZINGlitter_labile + &
!                                           GRAZINGlitter_foliage+GRAZINGlitter_roots) * step_length))
!                som     = max(0d0,som)
!
!            endif ! carry out grazing?

        endif ! Is grazing plausible

    endif ! end grazing process

  end subroutine grass_grazing
  !
  !-----------------------------------------------------------------
  !
  subroutine plant_canopy_phenology(nodays, gsi_lag_steps, step, time,                     & ! Timing
                                    avgTmax_C, photoperiod_s, avgVPD_Pa,                   & ! GSI forcings
                                    leafT_min, leafT_max, leafP_min, leafP_max,            & ! GSI parameters 
                                    leafV_min, leafV_max, leaf_phenology_threshold,        & !
                                    potential_labile_turnover, potential_foliage_turnover, & ! 
                                    lca, gpp_return_threshold, available_labile, foliage,  & ! To calculate GPP return
                                    current_gpp, delta_gpp_gCgC,                           & ! 
                                    alloc_leaf_fraction, alloc_leaf_gCm2day,               & ! Output variables
                                    leaf_litter_fraction, leaf_litter_gCm2day,             &
                                    leafT_limit, leafP_limit, leafV_limit, gsi_gradient, gsi, mV) 
    use carbon_model_memory

       ! Subroutine deals with the determining of allocated carbon to foliage from 
       ! the labile / non-structural carbohydrates pool and senescence of foliage 
       ! to fine litter. Growth and loss are determined by a canopy growing season index (GSI)
       ! which is the product of 3-limiting terms represented by linear functions of
       ! temperature, vapour pressure deficit and photoperiod. The linear functions are scaled
       ! 0-1 ensuring the product is also 0-1 and applied to a potential growth and loss parameters.
       ! Growth is a function of GSI, while senescene is a function of 1-GSI with the decision on
       ! which process is being applied is based on the dynamics of the local GSI history (i.e. going up or down).
       ! Where growth is considered a simple economic test is also applied to ensure that the potential 
       ! photosynthetic return is greater than a parameterised value.

       ! Here we calculate the Growing Season Index based on
       ! Jolly et al. A generalized, bioclimatic index to predict foliar
       ! phenology in response to climate Global Change Biology, Volume 11, page 619-632,
       ! 2005 (doi: 10.1111/j.1365-2486.2005.00930.x)
       ! Stoeckli, R., T. Rutishauser, I. Baker, M. A. Liniger, and A. S.
       ! Denning (2011), A global reanalysis of vegetation phenology, J. Geophys. Res.,
       ! 116, G03020, doi:10.1029/2010JG001545.

       implicit none

      type(model_working_variables) :: mV

       ! Arguments
       integer, intent(in) :: step, & ! current model time step
                            nodays, & ! Number of time steps in the analysis
                     gsi_lag_steps    ! Number of model time steps over which GSI calculations will be made.
       double precision, intent(in) :: time, & ! number of days in current time step
                                  avgTmax_C, & ! GSI period average for maximum air temperature (Celcius)
                              photoperiod_s, & ! GSI period average for photo period (seconds)
                                  avgVPD_Pa, & ! GSI period average for vapour pressure deficit (Pa)
                                  leafT_min, & ! Minimum temperaure for canopy growth (Kelvin)
                                  leafT_max, & ! Maximum temperaure for canopy growth (Kelvin)
                                  leafP_min, & ! Minimum photo period for canopy growth (seconds)
                                  leafP_max, & ! Maximum photo period for canopy growth (seconds)
                                  leafV_min, & ! Minimum vapour pressure deficit for canopy growth (Pa)
                                  leafV_max, & ! Maximum vapour pressure deficit for canopy growth (Pa)
                   leaf_phenology_threshold, & ! GSI gradient threshold above which leaf growth is possible 
                  potential_labile_turnover, & ! Potential fractional rate of labile turnover to foliage (0-1)
                 potential_foliage_turnover, & ! Potential fractional rate of foliage turnover to fine litter (0-1)                    
                                        lca, & ! leaf carbon per leaf area (gC/m2)
                       gpp_return_threshold, & ! GPP return for growth to go ahead (gC/gC/m2/day)
                           available_labile, & ! labile C available to spend this time step (gC/m2)
                                    foliage, & ! foliage C in this time step (gC/m2) 
                                current_gpp    ! current time step GPP estimate (gC/m2/day)
       double precision, intent(out) :: &
                             delta_gpp_gCgC, & ! change in gross primary production (gC/gCinvested)
                        alloc_leaf_fraction, & ! allocation to lablile to leaf (0-1)
                         alloc_leaf_gCm2day, & ! allocation to labile to leaf (gC/m2/day)
                       leaf_litter_fraction, & ! leaf litter fall fraction (0-1)
                        leaf_litter_gCm2day, & ! leaf litter fall (gC/m2/day)
                                leafT_limit, & ! temperature limitation on foliage (0-1)
                                leafP_limit, & ! photperiod limitation on foliage (0-1)
                                leafV_limit, & ! vapour pressure deficit limitation on foliage (0-1)
                               gsi_gradient    ! GSI gradient over the lag period (-1->1 / day)
       double precision, dimension(nodays), intent(inout) ::    &
                                        gsi    ! canopy Growing Season Index (gsi, 0-1)   

       ! Local variables
       integer :: interval
       double precision :: rescale, leafT_adj, lai_orig, gs_orig, scaling_orig, leaf_investment

       ! Set intial value of all output variables
       alloc_leaf_fraction = 0d0 ; alloc_leaf_gCm2day = 0d0
       leaf_litter_fraction = 0d0 ; leaf_litter_gCm2day = 0d0
       leafT_limit = 0d0 ; leafP_limit = 0d0 ; leafV_limit = 0d0

       !
       ! Calculate the current time steps GSI value
       !

       ! temperature limitation, then restrict to 0-1; unit change K-> oC
       leafT_limit = (avgTmax_C-(leafT_min-freeze)) / (leafT_max-leafT_min)
       leafT_limit = min(1d0,max(0d0,leafT_limit))
       ! photoperiod limitation, then restrict to 0-1
       leafP_limit = (photoperiod_s-leafP_min) / (leafP_max-leafP_min)
       leafP_limit = min(1d0,max(0d0,leafP_limit))
       ! VPD limitation, then restrict to 0-1
       leafV_limit = 1d0 - ( (avgVPD_Pa-leafV_min) / (leafV_max-leafV_min) )
       leafV_limit = min(1d0,max(0d0,leafV_limit))

       ! Calculate current time steps GSI
       gsi(step) = leafT_limit * leafP_limit * leafV_limit
!print*,gsi(step),leafT_limit,leafP_limit,leafV_limit
       !
       ! Calculate the GSI gradient across the GSI lag period
       !

       ! Determine the current GSI values and period for the linear fit.
       ! This code must account for the beginning of a simulation where 
       ! there is no existing history.
       if (step < gsi_lag_steps) then
           if (step == 1) then
               mV%gsi_lag_history(2) = gsi(step)
               interval = 2
           else
               mV%gsi_lag_history(1:step) = gsi(1:step)
               interval = step
           endif
       else
           mV%gsi_lag_history(1:gsi_lag_steps) = gsi((step-gsi_lag_steps+1):step)
           interval = gsi_lag_steps
       end if
       ! Now calculate the linear gradient
       gsi_gradient = linear_model_gradient(mV%gsi_lag_days(1:interval),mV%gsi_lag_history(1:interval),interval)

       ! We can only allocate if we have labile to spend and a GSI gradient above the leaf phenology threshold.
       if (gsi_gradient > leaf_phenology_threshold .and. gsi(step) > vsmall .and. available_labile > 0d0) then

           !
           ! Labile allocation to plant tissues (gC/m2/day)
           !

           ! Labile to foliage rate (gC.m-2.day-1)
           alloc_leaf_fraction = potential_labile_turnover*gsi(step)

           ! 
           ! Apply (loose-)optimality theory for the proposed allocation
           !

           ! Finally quantify the impact of increasing LAI on GPP, less Rd(24)
           if (alloc_leaf_fraction > 0d0) then
               ! Store the existing LAI and canopy_scaling
               lai_orig = mV%lai ; scaling_orig = mV%leaf_canopy_light_scaling ; gs_orig = mV%stomatal_conductance
               ! Calculate the total time step investment
               alloc_leaf_gCm2day = available_labile * (1d0-(1d0-alloc_leaf_fraction)**time)
               ! Calculate the new LAI based C investment, less allocation to growth respiration
               mV%lai = mV%lai + ((alloc_leaf_gCm2day * one_Rg_fraction) / lca)
               ! Update the shortwave radiation
               call calculate_shortwave_balance(mV)
               ! Update acm_gpp_stage_1
               call acm_gpp_stage_1(mV)      
               ! Update stomatal conductance 
               mV%stomatal_conductance = (mV%stomatal_conductance / scaling_orig) * mV%leaf_canopy_light_scaling
               ! Estimate the change in gross primary production                              
               delta_gpp_gCgC = (acm_gpp_stage_2(mV%stomatal_conductance, mV) * umol_to_gC * mV%dayl_seconds) &
                               - current_gpp
               ! rescale GPP to per gC investment but including the C gone to growth respiration
               delta_gpp_gCgC = delta_gpp_gCgC / alloc_leaf_gCm2day
               ! If non-economical do not grow                               
               if (delta_gpp_gCgC < gpp_return_threshold) then
                   alloc_leaf_fraction = 0d0 ; alloc_leaf_gCm2day = 0d0
               else
                   ! Rescale to the per day flux to match the model timestep
                   alloc_leaf_gCm2day = alloc_leaf_gCm2day / time
               end if
               ! Return initial values
               mV%lai = lai_orig ; mV%stomatal_conductance = gs_orig
               ! Update the shortwave radiation
               call calculate_shortwave_balance(mV)
               ! Update acm_gpp_stage_1
               call acm_gpp_stage_1(mV)      
           end if ! alloc_leaf_fraction > 0

       else if (gsi_gradient <= leaf_phenology_threshold .and. gsi(step) <= vsmall) then

           !
           ! Leaf fall to litter (gC/m2/day)
           !

           ! Estimate the fractional loss rate of foliage to litter
           leaf_litter_fraction = potential_foliage_turnover*(1d0-gsi(step))
           ! Estimate the absolute flux value loss of foliage to litter
           leaf_litter_gCm2day = foliage * (1d0-(1d0-leaf_litter_fraction)**time)/time

       end if ! gsi_gradient > leaf_phenology_threshold .and. available_labile > 0

       ! Return subroutine
       return

  end subroutine plant_canopy_phenology
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
  !---------------------------------------------------------------------
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
  ! Functions below this line
  !
  !------------------------------------------------------------------
  !
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
  double precision function calculate_declination(doy)
    use carbon_model_memory

    implicit none

     ! Declare arguments
     double precision, intent(in) :: doy

     ! Declination calculation
     ! NOTE: 0.002739726d0 = 1/365
     !    dec = - asin( sin( 23.45d0 * deg_to_rad ) * cos( 2d0 * pi * ( doy + 10d0 ) / 365d0 ) )
     !    dec = - asin( sin_dayl_deg_to_rad * cos( two_pi * ( doy + 10d0 ) / 365d0 ) )
     calculate_declination = - asin( sin_dayl_deg_to_rad * cos( two_pi * ( doy + 10d0 ) * 0.002739726d0 ) )

     ! return to user
     return

  end function calculate_declination  
  !
  !------------------------------------------------------------------
  !
  double precision function linear_model_gradient(x,y,interval)
    use carbon_model_memory

    ! Function to calculate the gradient of a linear model for a given dependent
    ! variable (y) based on predictive variable (x). The typical use of this
    ! function will in fact be to assume that x is time.

    implicit none

    ! declare input variables
    integer :: interval
    double precision, dimension(interval) :: x,y

    ! declare local variables
!    double precision :: sum_x, sum_y!, sumsq_x,sum_product_xy
    double precision :: sum_x, sum_y, sumsq_x, sum_product_xy
    integer :: j

    ! single-pass accumulation loop replacing four separate sum() reductions
!    ! calculate the sum of x
!    sum_x = sum(x)
!    ! calculate the sum of y
!    sum_y = sum(y)
    sum_x = 0d0 ; sum_y = 0d0 ; sumsq_x = 0d0 ; sum_product_xy = 0d0
    do j = 1, interval
       sum_x          = sum_x          + x(j)
       sum_y          = sum_y          + y(j)
       sumsq_x        = sumsq_x        + x(j)*x(j)
       sum_product_xy = sum_product_xy + x(j)*y(j)
    end do
    ! calculate the sum of squares of x
    !sumsq_x = sum(x*x)
    ! calculate the sum of the product of xy
    !sum_product_xy = sum(x*y)
    ! calculate the gradient
    !linear_model_gradient = ( (dble(interval)*sum_product_xy) - (sum_x*sum_y) )
    !&
    !                      / ( (dble(interval)*sumsq_x) - (sum_x*sum_x) )
!    ! Linear regression done as single line to reduce assignment requirements
!    linear_model_gradient = ( (dble(interval)*sum(x*y)) - (sum_x*sum_y) ) &
!                          / ( (dble(interval)*sum(x*x)) - (sum_x*sum_x) )
    linear_model_gradient = ( (dble(interval)*sum_product_xy) - (sum_x*sum_y) ) &
                          / ( (dble(interval)*sumsq_x)        - (sum_x*sum_x) )

    ! for future reference here is how to calculate the intercept
!    intercept = ( (sum_y*sumsq_x) - (sum_x*sum_product_xy) ) &
!              / ( (dble(interval)*sumsq_x) - (sum_x*sum_x) )

    ! don't forget to return to the user
    return

  end function linear_model_gradient
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
!
!--------------------------------------------------------------------
!
end module CARBON_MODEL_MOD
