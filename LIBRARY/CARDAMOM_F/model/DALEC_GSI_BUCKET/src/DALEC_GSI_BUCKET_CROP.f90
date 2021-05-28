
module CARBON_MODEL_CROP_MOD

implicit none

! make all private
private

! explicit publics
public :: CARBON_MODEL_CROP    &
         ,resp_rate_temp_coeff &
         ,ts_length            &
         ,sec_in_hour

! declare some module level variables
double precision, parameter :: sec_in_day = 86400d0
double precision, parameter :: sec_in_hour = 3600d0
double precision, parameter :: pi    = 3.14159265d0
double precision, parameter :: deg_to_rad = pi/180d0

! variables local to this module..
integer ::   plough_day, & ! day-of-year when field is ploughed  (default)
                sow_day, & ! day-of-year when field is sown      (default)
            harvest_day, & ! day-of-year when field is harvested (default)
                  stmob, & ! remoblise stem C to labile (1 = on)
 turnover_labile_switch    ! begin turnover of labile C

logical :: vernal_calcs, &  ! do vernalisation calculations?
               ploughed, &  !
        use_seed_labile, & != .False. ! whether to use seed labile for growth
                   sown, & != .False. ! has farmer sown crop yet?
                emerged    != .False. ! has crop emerged yet?

double precision ::           ts_length, & ! time step length in hours
                            step_of_day, & ! current step of the day (default = 1)
                           steps_in_day, & ! number of steps in a day (default = 1)
                                gpp_acm, & ! gross primary productivity (gC.m-2.day-1)
                    stock_storage_organ, & ! storage organ C pool, i.e. the desired crop (gC.m--2)
                     stock_dead_foliage, & ! dead but still standing foliage (gC.m--2)
                        stock_resp_auto, & ! autotrophic respiration pool (gC.m--2)
                           stock_labile, & ! labile C pool (gC.m--2)
                          stock_foliage, & ! foliage C pool (gC.m--2)
                             stock_stem, & ! stem C pool (gC.m--2)
                            stock_roots, & ! roots C pool (gC.m--2)
                           stock_litter, & ! litter C pool (gC.m--2)
                    stock_soilOrgMatter, & ! SOM C pool (gC.m--2)
                              resp_auto, & ! autotrophic respiration (gC.m-2.t-1)
                          resp_h_litter, & ! litter heterotrophic respiration (gC.m-2.t-1)
                   resp_h_soilOrgMatter, & ! SOM heterotrophic respiration (gC.m-2)
                                    npp, & ! net primary productivity (gC.m-2.t-1)
                              nee_dalec, & ! net ecosystem exchange (gC.m-2.t-1)
                                     DS, & ! Developmental state and initial condition
                                    LCA, & ! leaf mass area (gC.m-2)
            mean_alloc_to_storage_organ, & ! rolling average allocation of GPP to storage organ (gC.m-2)
        mean_alloc_to_storage_organ_old, & ! ...same but previous value...
                     decomposition_rate, & ! decomposition rate (frac / hr)
                     frac_GPP_resp_auto, & ! fraction of GPP allocated to autotrophic carbon pool
                  turnover_rate_foliage, & ! turnover rate of foliage (frac/hr)
                     turnover_rate_stem, & ! same for stem
                   turnover_rate_labile, & ! same for labile
                turnover_rate_resp_auto, & ! same for autotrophic C pool
                 resp_cost_labile_trans, & ! labile lost to respiration per gC labile to GPP
             mineralisation_rate_litter, & ! mineralisation rate of litter
      mineralisation_rate_soilOrgMatter, & ! mineralisation rate of SOM
                                  PHUem, & ! emergance value for phenological heat units
                                    PHU, & ! phenological heat units
                                 DR_pre, & ! development rate coefficient DS 0->1
                                DR_post, & ! development rate coefficient DS 1->2
                                   tmin, & ! min temperature for development
                                   tmax, & ! max temperature for development
                                   topt, & ! optimum temperature for development
                                 tmin_v, & ! min temperature for vernalisation
                                 tmax_v, & ! max temperature for vernalisation
                                 topt_v, & ! optimim temperature for vernalisation
                                    VDh, & ! effective vernalisation days when plants are 50 % vernalised
                                     VD, & ! count of vernalisation days
                               RDRSHMAX, & ! maximum rate of self shading turnover
                                   PHCR, & ! critical value of photoperiod for development
                                   PHSC, & ! photoperiod sensitivity
                                   raso, & ! rolling average for alloc to storage organ
                               max_raso, & ! maximum value for rolling average alloc to storage organ
                                  BM_EX, & !
                                     HI, & !
                                  yield, & ! crop yield (gC.m-2)
                     alloc_to_resp_auto, & ! amount of carbon to allocate to autotrophic respiration pool
                    turnover_rate_roots, & ! turnover over rate of roots interpolated each time step
                                gso_max, & !
                           max_raso_old, & !
                               raso_old, & !
            resp_cost_labile_to_foliage, & ! respiratory cost of moving carbon..from labile to foliage pools
            resp_cost_foliage_to_labile, & ! ..from foliage to labile pools
                              resp_rate, & ! rate of respiration at given temperature
                                 Cshoot, & !
                                     DR, & !
                        fol_frac_intpol, & !
                       stem_frac_intpol, & !
                               fP,fT,fV, & !
                                  remob, & !
                       root_frac_intpol, & !
                      shoot_frac_intpol, & !
                                 avtemp, & !
                 alloc_to_storage_organ, & !
                     litterfall_foliage, & !
                        litterfall_stem, & !
                       litterfall_roots, & !
                          decomposition, & !
                              npp_shoot, & !
                      alloc_from_labile, & !
                        alloc_to_labile, & !
                         alloc_to_roots, & !
                       alloc_to_foliage, & !
                          alloc_to_stem, & !
                                raremob, & !
                                  RDRSH, & !
                                  RDRDV, & !
                                    RDR

  !
  ! some hardcoded crop parameters
  !

  ! defines Q10 = 2 in exponential temperature response for heterotrophic
  ! respiration
  double precision, parameter :: resp_rate_temp_coeff = 0.0693d0
  ! residue fraction of leaves left post harvest
  double precision, parameter :: lv_res = 0.1d0
  ! residue fraction of stem left post harvest
  double precision, parameter :: st_res = 0.1d0
  ! LAI above which self shading turnover occurs
  double precision, parameter :: LAICR = 4d0
  ! allocation to storage organ relative to GPP
  double precision, parameter :: rel_gso_max = 0.35d0

  ! 'save' indicates all these module variables should be held, even if the
  ! module itself goes out of scope.
  save

contains
!
!--------------------------------------------------------------------
!
  subroutine CARBON_MODEL_CROP(start,finish,met,pars,deltat,nodays,lat,lai_out &
                       ,NEE_out,FLUXES,POOLS,pft,nopars,nomet,nopools,nofluxes &
                       ,GPP_out,stock_seed_labile,DS_shoot,DS_root,fol_frac    &
                       ,stem_frac,root_frac,DS_LRLV,LRLV,DS_LRRT,LRRT)

    use CARBON_MODEL_MOD, only: arrhenious,acm_gpp_stage_1,acm_gpp_stage_2,calculate_transpiration         & ! Subroutine / functions
                               ,calculate_soil_evaporation      &
                               ,calculate_wetcanopy_evaporation,meteorological_constants                   &
                               ,calculate_stomatal_conductance,calculate_radiation_balance                 &
                               ,calculate_daylength,opt_max_scaling,calculate_Rtot,saxton_parameters       &
                               ,calculate_aerodynamic_conductance,calculate_update_soil_water              &
                               ,initialise_soils,update_soil_initial_conditions                            &
                               ,freeze,co2comp_half_sat_conc,kc_saturation,kc_half_sat_conc,min_drythick   & ! parameter
                               ,seconds_per_day,avN,iWUE,NUE,pn_max_temp,pn_opt_temp,pn_kurtosis,vsmall    &
                               ,min_root,top_soil_depth,max_depth,root_k,minlwp,min_layer                  &
                               ,soil_frac_clay,soil_frac_sand                                              &
                               ,co2comp_saturation,drythick,dayl_hours,dayl_seconds,dayl_seconds_1         & ! variables
                               ,seconds_per_step,fine_root_biomass,root_reach,previous_depth               &
                               ,deltat_1,water_flux_mmolH2Om2s,layer_thickness,meant,stomatal_conductance  &
                               ,co2_half_sat,co2_comp_point,mint,maxt,swrad,co2,doy,leafT,ceff             &
                               ,wind_spd,vpd_kPa,lai,days_per_step,days_per_step_1,dayl_hours_fraction     &
                               ,wSWP,SWP,SWP_initial,wSWP_time,soil_waterfrac,soil_waterfrac_initial       &
                               ,porosity,porosity_initial,field_capacity,field_capacity_initial            &
                               ,rainfall,canopy_storage,intercepted_rainfall,snow_storage,snow_melt        &
                               ,airt_zero_fraction,snowfall,root_biomass,gs_demand_supply_ratio            &
                               ,gs_total_canopy,gb_total_canopy,canopy_par_MJday_time,potential_conductance&
                               ,convert_ms1_mol_1,aerodynamic_conductance,canopy_par_MJday

    ! DALEC crop model modified from Sus et al., (2010)


    implicit none

    ! declare input variables
    integer, intent(in) :: start    &
                          ,finish   &
                          ,nopars     & ! number of paremeters in vector
                          ,pft        & ! plant functional type
                          ,nomet      & ! number of meteorological fields
                          ,nofluxes   & ! number of model fluxes
                          ,nopools    & ! number of model pools
                          ,nodays       ! number of days in simulation

    double precision, intent(in) :: met(nomet,nodays)   & ! met drivers
                         ,stock_seed_labile             & ! seed carbon to get things going
                         ,deltat(nodays)                & ! time step in decimal days
                         ,pars(nopars)                  & ! number of parameters
                         ,lat                 ! site latitude (degrees)

    double precision, dimension(:), intent(inout) ::          DS_shoot, & !
                                                               DS_root, & !
                                                              fol_frac, & !
                                                             stem_frac, & !
                                                             root_frac, & !
                                                               DS_LRLV, & !
                                                                  LRLV, & !
                                                               DS_LRRT, & !
                                                                  LRRT    !

    double precision, dimension(nodays), intent(inout) :: lai_out & ! leaf area index
                                               ,GPP_out & ! Gross primary productivity
                                               ,NEE_out   ! net ecosystem exchange of CO2

    double precision, dimension((nodays+1),nopools), intent(inout) :: POOLS ! vector of ecosystem pools

    double precision, dimension(nodays,nofluxes), intent(inout) :: FLUXES ! vector of ecosystem fluxes

    ! declare local variables
    double precision :: airt_weighting(3) &
                                      ,ET & ! Evapotranspiration (kg.m-2.day-1)
                           ,transpiration &
                         ,soilevaporation &
                          ,wetcanopy_evap &
                        ,snow_sublimation &
                                 ,deltaWP & !
                                    ,infi   ! used to calculate infinity for diagnositc

    integer :: nxp,n

    ! met drivers are:
    ! 1st run day
    ! 2nd min daily temp (oC)
    ! 3rd max daily temp (oC)
    ! 4th Radiation (MJ.m-2.day-1)
    ! 5th CO2 (ppm)
    ! 6th DOY

    ! POOLS are:
    ! 1 = labile
    ! 2 = foliar
    ! 3 = root
    ! 4 = wood
    ! 5 = litter
    ! 6 = som
    ! 7 = autotrophic
    ! 8 = storage organ C

    ! FLUXES are:
    ! 1 = GPP
    ! 2 = temprate
    ! 3 = respiration_auto
    ! 4 = leaf production
    ! 5 = labile production
    ! 6 = root production
    ! 7 = wood production
    ! 8 = labile release
    ! 9 = alloc to storage
    ! 10 = leaf litter production
    ! 11 = woodlitter production
    ! 12 = rootlitter production
    ! 13 = respiration het litter
    ! 14 = respiration het som
    ! 15 = litter2som (decomposition)
    ! 16 = alloc to autotrophic pool

    ! PARAMETERS
    ! 16 values

    ! p(1) decomposition rate (frac/hr)
    ! p(2) Fraction of GPP allocated to autotrophic C pool
    ! p(3) DR coef for DS (0->1)
    ! p(4) DR coef for DS (1->2)
    ! p(5) turnover rate of foliage (frac/hr)
    ! p(6) Turnover rate of wood/stem (frac/hr)
    ! p(7) maximum rate of foliar turnover due to self shading
    ! p(8) effective vernalisation days when plant is 50 % vernalised
    ! p(9) mineralisation rate of som
    ! p(10) mineralisation rate of litter
    ! p(11) = log10(avgN)
    ! p(12) = sow day
    ! p(13) = labile lost to respiration per gC labile top GPP
    ! p(14) = phenological heat units needed for emergence
    ! p(15) ! harvest day (doy)
    ! p(16) ! plough day (doy)
    ! p(17) ! leaf mass area (gC.m-2)
    ! p18,p19,p20,p21,p22,p23,p24,p25 = labile, foliar, roots, stem, litter,
    ! som,
    ! autotrophic and storage organ pools respectively
    ! p(26) ! min temperature for development
    ! p(27) ! max temperature for development
    ! p(28) ! optimum temperature for development
    ! p(29) ! min temperature for vernalisation
    ! p(30) ! max temperature for vernalisation
    ! p(31) ! optimim temperature for vernalisation
    ! p(32) ! critical value of photoperiod for development
    ! p(33) ! photoperiod sensitivity
    ! p(34) ! turnover rate of labile C
    ! p(35) ! turnover rate of autotrophic C

    ! zero some values
    lai_out(1:nodays) = 0d0 ; NEE_out(1:nodays) = 0d0 ; GPP_out(1:nodays) = 0d0
    FLUXES(1:nodays,1:nofluxes) = 0d0 ; POOLS(1:nodays,1:nopools) = 0d0

    ! load ACM-GPP-ET parameters
    NUE = 1.182549d+01   ! Photosynthetic nitrogen use efficiency at optimum temperature (oC)
                         ! ,unlimited by CO2, light and photoperiod (gC/gN/m2leaf/day)
    avN = 10d0**pars(11) ! foliar N
    ceff = avN*NUE
    deltaWP = minlwp     ! leafWP-soilWP (i.e. -2-0)

    ! plus ones being calibrated
    root_k = pars(36) ; max_depth = pars(37)

    ! length of time step in hours..
    ts_length = ((sum(deltat)/dble(nodays)) * sec_in_day) / sec_in_hour
    ! steps per day ; set step_of_day to steps_in_day as in CARDAMOM we will
    ! never be running less than daily time step
    steps_in_day = 1d0/(sum(deltat)/dble(nodays)) ; step_of_day = steps_in_day

    ! parameters from file
    decomposition_rate                = pars(1) / 24d0  ! decomposition rate (day->hr)
    frac_GPP_resp_auto                = pars(2)  ! fraction of GPP allocated to autotrophic carbon pool
    DR_pre                            = pars(3)  ! development rate coefficient DS (0->1)
    DR_post                           = pars(4)  ! development rate coefficient DS (1->2)
    turnover_rate_foliage             = pars(6) / 24d0 ! pars(5)  ! turnover_rate of foliage (day->hr)
    turnover_rate_stem                = pars(6) / 24d0 ! turnover rate of stem (day->hr)
    RDRSHMAX                          = pars(7) / 24d0 ! maximum rate of foliar turnover due to self shading (day->hr)
    VDh                               = pars(8)  ! effective vernalisation days when plants are 50 % vernalised
    mineralisation_rate_litter        = pars(9) / 24d0 ! mineralisation rate litter (day->hr)
    mineralisation_rate_soilOrgMatter = pars(10)/ 24d0 ! mineralisation rate som (day->hr)
    sow_day                           = nint(mod(pars(12),365.25d0)) ! sow day (doy)
    resp_cost_labile_trans            = pars(13) ! labile lost to respiration per gC labile to GPP
    PHUem                             = pars(14) ! phenological heat units required for emergence
    harvest_day                       = nint(mod(pars(15),365.25d0)) ! nint(mod(pars(15),365.25)) ! harvest day (doy)
    plough_day                        = nint(mod(pars(12)-2d0,365.25d0)) ! nint(mod(pars(16),365.25)) ! plough day (doy)
    LCA                               = pars(17) ! leaf mass area (gC.m-2)
    tmin                              = pars(26)-273.15d0 ! min temperature for development
    tmax                              = pars(27)-273.15d0 ! max temperature for development
    topt                              = pars(28)-273.15d0 ! optimum temperature for development
    tmin_v                            = pars(29)-273.15d0 ! min temperature for vernalisation
    tmax_v                            = pars(30)-273.15d0 ! max temperature for vernalisation
    topt_v                            = pars(31)-273.15d0 ! optimim temperature for vernalisation
    PHCR                              = pars(32) ! critical value of photoperiod for development
    PHSC                              = pars(33) ! photoperiod sensitivity
    turnover_rate_labile              = pars(34)/ 24d0 ! turnover rate labile C (day->hr)
    turnover_rate_resp_auto           = pars(35)/ 24d0 ! turnover rate of autotrophic carbon for respiration (day->hr)

    if (start == 1) then

        ! load stocks in first time step
        stock_labile                      = pars(18) ! labile C
        stock_foliage                     = pars(19) ! foliar C
        stock_roots                       = pars(20) ! root C
        stock_stem                        = pars(21) ! stem / wood C
        stock_litter                      = pars(22) ! litter C
        stock_soilOrgMatter               = pars(23) ! som C
        stock_resp_auto                   = pars(24) ! autotrophic resp pool
        stock_storage_organ               = pars(25) ! storage organ (i.e. desired crop)

        ! assigning initial conditions
        POOLS(1,1) = pars(18) ! Clabile
        POOLS(1,2) = pars(19) ! Cfoliar
        POOLS(1,3) = pars(20) ! Croots
        POOLS(1,4) = pars(21) ! Cstructural
        POOLS(1,5) = pars(22) ! Clitter
        POOLS(1,6) = pars(23) ! Csom
        POOLS(1,7) = pars(24) ! Cauto
        ! POOLS(1,8) ! WATER IN ROOT ZONE ASSIGNED LATER
        POOLS(1,9) = pars(25) ! Cstorage

        ! logical switches
        vernal_calcs    = .true.
        ploughed        = .false.
        sown            = .false.
        use_seed_labile = .false.
        emerged         = .false.

        ! pair incoming variables to local module levels
        ! finally set some initial conditions
        avtemp = 0d0
        yield = 0d0
        DS = -1d0
        DR = 0d0
        fV = 0d0 ; fT = 0d0 ; fP = 0d0
        mean_alloc_to_storage_organ = 0d0
        mean_alloc_to_storage_organ_old = 0d0
        PHU = 0d0
        VD = 0d0
        BM_EX = 0d0
        HI = 0d0
        stock_dead_foliage = 0d0
        alloc_to_labile = 0d0
        stmob = 0
        max_raso = 0d0
        raso = 0d0
        RDRDV = 0d0
        max_raso_old = 0d0
        raso_old  = 0d0

        ! SHOULD TURN THIS INTO A SUBROUTINE CALL AS COMMON TO BOTH DEFAULT AND CROPS
        if (.not.allocated(deltat_1)) then

           allocate(deltat_1(nodays),wSWP_time(nodays),gs_demand_supply_ratio(nodays), &
                    gs_total_canopy(nodays), gb_total_canopy(nodays), &
                    canopy_par_MJday_time(nodays))
           deltat_1 = deltat**(-1d0)
           ! zero variables not done elsewhere
           water_flux_mmolH2Om2s = 0d0
           ! initialise some time invarient parameters
           call saxton_parameters(soil_frac_clay,soil_frac_sand)
           call initialise_soils(soil_frac_clay,soil_frac_sand)
           call update_soil_initial_conditions(pars(38))
           ! save the initial conditions for later
           soil_waterfrac_initial = soil_waterfrac
           SWP_initial = SWP
           field_capacity_initial = field_capacity
           porosity_initial = porosity

        else

           water_flux_mmolH2Om2s = 0d0
           field_capacity = field_capacity_initial
           porosity = porosity_initial

           ! input initial soil water fraction then
           ! update SWP and soil conductivity accordingly
           call update_soil_initial_conditions(pars(38))

        endif

        ! load some needed module level values
        lai = POOLS(1,2)/pars(17)
        seconds_per_step = deltat(1) * seconds_per_day
        days_per_step =  deltat(1)
        days_per_step_1 = deltat_1(1)
        mint = met(2,1)  ! minimum temperature (oC)
        maxt = met(3,1)  ! maximum temperature (oC)
        meant = (maxt+mint)*0.5d0  ! mean air temperature (oC)

        ! zero evapotranspiration for beginning
        ET = 0d0
        ! initialise root reach based on initial conditions
        fine_root_biomass = max(min_root,POOLS(1,3)*2d0)
        root_biomass = fine_root_biomass
        root_reach = max_depth * root_biomass / (root_k + root_biomass)
        ! Determine initial soil layer thickness
        layer_thickness(1) = top_soil_depth
        layer_thickness(2) = max(min_layer,root_reach-top_soil_depth)
        layer_thickness(3) = max_depth - sum(layer_thickness(1:2))
        layer_thickness(4) = top_soil_depth
        previous_depth = max(top_soil_depth,root_reach)
        ! needed to initialise soils
        call calculate_Rtot
        ! used to initialise soils
        call calculate_update_soil_water(0d0,0d0,0d0,ET) ! assume no evap or rainfall
        ! store soil water content of the rooting zone (mm)
        POOLS(1,8) = 1d3*soil_waterfrac(1)*layer_thickness(1)

    else

        ! load ET from memory
        ET = FLUXES(start-1,19)

    endif ! start == 1

    ! reset values
    intercepted_rainfall = 0d0 ; canopy_storage = 0d0 ; snow_storage = 0d0

    infi = 0d0

    !
    ! Begin looping through each time step
    !

    do n = start, finish

      !!!!!!!!!!
      ! assign drivers and update some prognostic variables
      !!!!!!!!!!

      ! Incoming drivers
      mint = met(2,n)  ! minimum temperature (oC)
      maxt = met(3,n)  ! maximum temperature (oC)
      leafT = maxt     ! initial canopy temperature (oC)
      swrad = met(4,n) ! incoming short wave radiation (MJ/m2/day)
      co2 = met(5,n)   ! CO2 (ppm)
      doy = ceiling(met(6,n)-(deltat(n)*0.5d0))   ! Day of year
      rainfall = max(0d0,met(7,n)) ! rainfall (kgH2O/m2/s)
      meant = (maxt+mint) * 0.5d0   ! mean air temperature (oC)
      airt_zero_fraction = (maxt-0d0) / (maxt-mint) ! fraction of temperture period above freezing
      wind_spd = met(15,n) ! wind speed (m/s)
      vpd_kPa = met(16,n)*1d-3  ! Vapour pressure deficit (Pa)

      ! states needed for module variables
      lai_out(n) = POOLS(n,2)/LCA
      lai = lai_out(n) ! leaf area index (m2/m2)

      ! calculate daylength in hours and seconds
      call calculate_daylength((doy-(deltat(n)*0.5d0)),lat)
      ! extract timing related values
      dayl_hours_fraction = dayl_hours * 0.04166667d0 ! 1/24 = 0.04166667
      dayl_seconds_1 = dayl_seconds ** (-1d0)
      seconds_per_step = seconds_per_day * deltat(n)
      days_per_step = deltat(n)
      days_per_step_1 = deltat_1(n)

      ! snowing or not...?
      snow_melt = 0d0 ; snowfall = 0d0
      if (mint < 0d0 .and. maxt > 0d0) then
          ! if minimum temperature is below freezing point then we weight the
          ! rainfall into snow or rain based on proportion of temperature below
          ! freezing
          snowfall = 1d0 - airt_zero_fraction
          snowfall = rainfall * snowfall ; rainfall = rainfall - snowfall
          ! Add rainfall to the snowpack and clear rainfall variable
          snow_storage = snow_storage + (snowfall*seconds_per_step)

          ! Also melt some of the snow
          snow_melt = airt_zero_fraction
          ! otherwise we assume snow is melting at 10 % per day light hour
          snow_melt = min(snow_storage, snow_melt * snow_storage * dayl_hours * 0.1d0 * deltat(n))
          snow_storage = snow_storage - snow_melt
      elseif (maxt < 0d0) then
          ! if whole day is below freezing then we should assume that all
          ! precipitation is snowfall
          snowfall = rainfall ; rainfall = 0d0
          ! Add rainfall to the snowpack and clear rainfall variable
          snow_storage = snow_storage + (snowfall*seconds_per_step)
      else if (mint > 0d0) then
          ! otherwise we assume snow is melting at 10 % per day light hour
          snow_melt = min(snow_storage, snow_storage * dayl_hours * 0.1d0 * deltat(n))
          snow_storage = snow_storage - snow_melt
      end if

      ! calculate the minimum soil & root hydraulic resistance based on total
      ! fine root mass ! *2*2 => *RS*C->Bio
      fine_root_biomass = max(min_root,POOLS(n,3)*2d0)
      root_biomass = fine_root_biomass
      ! estimate drythick for the current step
      drythick = max(min_drythick, top_soil_depth * min(1d0,1d0 - (soil_waterfrac(1) / porosity(1))))
      call calculate_Rtot

      ! Pass wSWP to output variable and update deltaWP between minlwp and
      ! current weighted soil WP
      wSWP_time(n) = wSWP ; deltaWP = min(0d0,minlwp-wSWP)

      ! calculate some temperature dependent meteorologial properties
      call meteorological_constants(maxt,maxt+freeze,vpd_kPa)
      ! calculate radiation absorption and estimate stomatal conductance
      call calculate_aerodynamic_conductance
      gb_total_canopy(n) = aerodynamic_conductance * convert_ms1_mol_1 * 1d3
      call calculate_radiation_balance
      canopy_par_MJday_time(n) = canopy_par_MJday
      call calculate_stomatal_conductance
      ! Estimate stomatal conductance relative to its minimum / maximum, i.e. how
      ! close are we to maxing out supply (note 0.01 taken from min_gs)
      gs_demand_supply_ratio(n) = (stomatal_conductance - 0.01d0) / (potential_conductance-0.01d0)
      ! Store the canopy level stomatal conductance (mmolH2O/m2/day)
      gs_total_canopy(n) = stomatal_conductance

      ! reallocate for crop model timings
      doy = met(6,n)

      ! GPP (gC.m-2.day-1)
      if (lai > vsmall .and. stomatal_conductance > vsmall) then
         call acm_gpp_stage_1 ; GPP_out(n) = max(0d0,acm_gpp_stage_2(stomatal_conductance))
         ! Canopy transpiration (kgH2O/m2/day)
         call calculate_transpiration(transpiration)
      else
         GPP_out(n) = 0d0
         transpiration = 0d0
      endif
      ! load GPP for crop model daily rate to total
      gpp_acm = GPP_out(n) * deltat(n)

      ! Canopy intercepted rainfall evaporation (kgH2O/m2/day)
      call calculate_wetcanopy_evaporation(wetcanopy_evap,canopy_storage)
      ! Soil surface (kgH2O.m-2.day-1)
      call calculate_soil_evaporation(soilevaporation)
      ! restrict transpiration to positive only
      transpiration = max(0d0,transpiration)

      ! if snow present assume that soilevaporation is sublimation of soil first
      snow_sublimation = 0d0
      if (snow_storage > 0d0) then
          snow_sublimation = soilevaporation
          if (snow_sublimation*deltat(n) > snow_storage) snow_sublimation = snow_storage * deltat_1(n)
          soilevaporation = soilevaporation - snow_sublimation
          snow_storage = snow_storage - (snow_sublimation * deltat(n))
      end if

      ! add any snow melt to the rainfall now that we have already dealt with the canopy interception
      rainfall = rainfall + (snow_melt / seconds_per_step)
      ! do mass balance (i.e. is there enough water to support ET)
      call calculate_update_soil_water(transpiration,soilevaporation,((rainfall-intercepted_rainfall)*seconds_per_day) &
                                      ,FLUXES(n,19))
      ! now that soil mass balance has been updated we can add the wet canopy
      ! evaporation (kg.m-2.day-1)
      FLUXES(n,19) = FLUXES(n,19) + wetcanopy_evap
      ! pass to local variable for soil mass balance
      ET = FLUXES(n,19)

      ! store soil water content of surface (mm)
      POOLS(n,8) = 1d3*soil_waterfrac(1)*layer_thickness(1)

      ! daily average of allocation to storage organ (needed to determine max.
      ! storage organ growth rate)
      mean_alloc_to_storage_organ_old = mean_alloc_to_storage_organ
      mean_alloc_to_storage_organ     = 0d0
      ! pass relevant variables into crop module memory
      avtemp = met(14,n) !meant

      ! calculate weighted air temperature value based on daily minimum, maximum
      ! and means. This minimises the error introduced when scaling between
      ! daily and sub-daily timesteps
      airt_weighting(1) = abs(met(3,n)-avtemp) / (met(3,n)-met(2,n))*0.5d0 ! maximum temperature weighting
      airt_weighting(2) = 0.5d0                                            ! mean temperature
      airt_weighting(3) = abs(met(2,n)-avtemp) / (met(3,n)-met(2,n))*0.5d0 ! minimum temperature weighting

      ! Heterotrophic respiration rate (Q10):  doubles with
      ! 10 degree temperature rise resprate from soil file = 0.0693
      resp_rate = 0d0
      resp_rate = resp_rate + ((0.5d0 * exp( resp_rate_temp_coeff * met(3,n) )) * airt_weighting(1))
      resp_rate = resp_rate + ((0.5d0 * exp( resp_rate_temp_coeff * avtemp   )) * airt_weighting(2))
      resp_rate = resp_rate + ((0.5d0 * exp( resp_rate_temp_coeff * met(2,n) )) * airt_weighting(3))
      !resp_rate = 0.5 * exp( resp_rate_temp_coeff * avtemp )

      ! determine development stage (DS)
      call development_stage(deltat(n))
      ! determine the carbon partitioning based on development stage
      call carbon_alloc_fractions(DS_shoot,DS_root,fol_frac,stem_frac,root_frac)
      ! begin carbon allocation for crops
      call calc_pools_crops(DS_LRRT,LRRT)
      ! conduct management updates at the end of the day
      call management_dates(stock_seed_labile,deltat(n))

      ! calculate the NEE
      NEE_out(n) = nee_dalec / deltat(n)

      ! GPP (gC.m-2.day-1)
      FLUXES(n,1) = GPP_out(n)
      ! temprate (i.e. temperature modified rate of metabolic activity))
      FLUXES(n,2) = resp_rate
      ! autotrophic respiration (gC.m-2.day-1)
      FLUXES(n,3) = resp_auto * steps_in_day !/ deltat(n)
      ! leaf production rate (gC.m-2.day-1)
      FLUXES(n,4) = alloc_to_foliage * steps_in_day !/deltat(n)
      ! labile production (gC.m-2.day-1)
      FLUXES(n,5) = (alloc_to_labile + remob) * steps_in_day !/deltat(n)
      ! root production (gC.m-2.day-1)
      FLUXES(n,6) = alloc_to_roots * steps_in_day !/deltat(n)
      ! wood production
      FLUXES(n,7) = alloc_to_stem * steps_in_day !/deltat(n)
      ! labile production
      FLUXES(n,8) = (alloc_from_labile + resp_cost_labile_to_foliage) * steps_in_day !/deltat(n)
      ! alloc to storage organ
      FLUXES(n,9) = alloc_to_storage_organ * steps_in_day !/deltat(n)
      ! total leaf litter production
      FLUXES(n,10) = litterfall_foliage * steps_in_day !/deltat(n)
      ! total wood production
      FLUXES(n,11) = litterfall_stem * steps_in_day !/deltat(n)
      ! total root litter production
      FLUXES(n,12) = litterfall_roots * steps_in_day !/deltat(n)
      ! respiration heterotrophic litter
      FLUXES(n,13) = resp_h_litter * steps_in_day !/deltat(n)
      ! respiration heterotrophic som
      FLUXES(n,14) = resp_h_soilOrgMatter * steps_in_day !/deltat(n)
      ! litter to som
      FLUXES(n,15) = decomposition * steps_in_day !/deltat(n)
      ! alloc to autotrophic pool
      FLUXES(n,16) = frac_GPP_resp_auto * GPP_out(n)
      ! harvest yield; convert to daily rate
      FLUXES(n,21) = yield / deltat(n)

      ! labile pool
      POOLS(n+1,1) = stock_labile
      ! foliar pool
      POOLS(n+1,2) = stock_foliage
      ! wood pool
      POOLS(n+1,4) = stock_stem
      ! root pool
      POOLS(n+1,3) = stock_roots
      ! litter pool
      POOLS(n+1,5) = stock_litter
      ! som pool
      POOLS(n+1,6) = stock_soilOrgMatter
      ! autotrophic pool
      POOLS(n+1,7) = stock_resp_auto
      ! storage organ pool
      POOLS(n+1,9) = stock_storage_organ

      do nxp = 1, nopools
         if (POOLS(n+1,nxp) /= POOLS(n+1,nxp) .or. POOLS(n+1,nxp) < 0d0) then
             print*,"step",n,"FLUXES",nxp
             print*,"met",met(:,n)
             print*,"POOLS",POOLS(n,:)
             print*,"FLUXES",FLUXES(n,:)
             print*,"POOLS+1",POOLS(n+1,:)
             print*,"wSWP",wSWP
             print*,"waterfrac",soil_waterfrac
             print*,"steps_in_day",steps_in_day
             print*,stock_labile, stock_foliage
             print*,stock_stem,stock_roots
             print*,stock_litter,stock_soilOrgMatter
             print*,stock_storage_organ,stock_resp_auto
             print*,gpp_acm,nee_dalec
             print*,resp_auto,resp_h_soilOrgMatter,resp_h_litter
             print*,"pars",pars
             print*,"DR",DR
             print*,"alloc_to_labile",alloc_to_labile,"remob",remob
             print*,"DR stuff",fT,fV,fP
             print*,"leaf","stem","rootlitter",litterfall_foliage,litterfall_stem,litterfall_roots
             print*,"daylength",dayl_hours,"VD",VD,"VDh",VDh
             print*,"avtemp",avtemp
             print*,"sown",sown,"emerged",emerged
             print*,"root_frac_intpol",root_frac_intpol
             print*,"npp_shoot",npp_shoot,"npp",npp
             print*,"RDR",RDR,"ts_length",ts_length
             stop
         endif
      enddo

      do nxp = 1, nofluxes
         if (nxp /= 19) then
            if (FLUXES(n,nxp) /= FLUXES(n,nxp) .or. FLUXES(n,nxp) < 0d0) then
                 print*,"Special: step",n,"FLUXES",nxp
                 print*,"met",met(:,n)
                 print*,"POOLS",POOLS(n,:)
                 print*,"FLUXES",FLUXES(n,:)
                 print*,"POOLS+1",POOLS(n+1,:)
                 print*,"wSWP",wSWP
                 print*,"waterfrac",soil_waterfrac
                 print*,"steps_in_day",steps_in_day
                 print*,stock_labile, stock_foliage
                 print*,stock_stem,stock_roots
                 print*,stock_litter,stock_soilOrgMatter
                 print*,stock_storage_organ,stock_resp_auto
                 print*,gpp_acm,nee_dalec
                 print*,resp_auto,resp_h_soilOrgMatter,resp_h_litter
                 print*,"pars",pars
                 print*,"DR",DR
                 print*,"alloc_to_labile",alloc_to_labile,"remob",remob
                 print*,"DR stuff",fT,fV,fP
                 print*,"leaf","stem","rootlitter",litterfall_foliage,litterfall_stem,litterfall_roots
                 print*,"daylength",dayl_hours,"VD",VD,"VDh",VDh
                 print*,"avtemp",avtemp
                 print*,"sown",sown,"emerged",emerged
                 print*,"root_frac_intpol",root_frac_intpol
                 print*,"npp_shoot",npp_shoot,"npp",npp
                 print*,"RDR",RDR,"ts_length",ts_length
                 stop
            end if
         else
             if (FLUXES(n,nxp) /= FLUXES(n,nxp)) then
                 print*,"Default: step",n,"FLUXES",nxp
                 print*,"met",met(:,n)
                 print*,"POOLS",POOLS(n,:)
                 print*,"FLUXES",FLUXES(n,:)
                 print*,"POOLS+1",POOLS(n+1,:)
                 print*,"wSWP",wSWP
                 print*,"waterfrac",soil_waterfrac
                 print*,"steps_in_day",steps_in_day
                 print*,stock_labile, stock_foliage
                 print*,stock_stem,stock_roots
                 print*,stock_litter,stock_soilOrgMatter
                 print*,stock_storage_organ,stock_resp_auto
                 print*,gpp_acm,nee_dalec
                 print*,resp_auto,resp_h_soilOrgMatter,resp_h_litter
                 print*,"pars",pars
                 print*,"DR",DR
                 print*,"alloc_to_labile",alloc_to_labile,"remob",remob
                 print*,"DR stuff",fT,fV,fP
                 print*,"leaf","stem","rootlitter",litterfall_foliage,litterfall_stem,litterfall_roots
                 print*,"daylength",dayl_hours,"VD",VD,"VDh",VDh
                 print*,"avtemp",avtemp
                 print*,"sown",sown,"emerged",emerged
                 print*,"root_frac_intpol",root_frac_intpol
                 print*,"npp_shoot",npp_shoot,"npp",npp
                 print*,"RDR",RDR,"ts_length",ts_length
                 stop
             endif
         end if
      enddo

      if (stock_labile < 0d0 .or. stock_foliage < 0d0 .or. stock_stem < 0d0 .or. &
          stock_roots < 0d0 .or. stock_litter < 0d0 .or. stock_soilOrgMatter < 0d0 .or. &
          stock_storage_organ < 0d0 .or. stock_resp_auto < 0d0 .or. &
          stock_labile /= stock_labile .or. stock_foliage /= stock_foliage .or. &
          stock_stem /= stock_stem .or. &
          stock_roots /= stock_roots .or. stock_litter /= stock_litter .or. &
          stock_soilOrgMatter /= stock_soilOrgMatter .or. &
          stock_storage_organ /= stock_storage_organ .or. &
          stock_resp_auto /= stock_resp_auto .or.  &
          gpp_acm < 0d0 .or. gpp_acm /= gpp_acm .or. resp_rate < 0d0 .or. &
          resp_rate /= resp_rate .or. decomposition < 0d0 .or. alloc_from_labile < 0d0 .or. &
          resp_cost_labile_to_foliage < 0d0 .or. alloc_to_foliage < 0d0 .or. &
          alloc_to_stem < 0d0 .or. alloc_to_roots < 0d0 .or. remob < 0d0 .or. &
          alloc_from_labile < 0d0 .or. resp_cost_labile_to_foliage < 0d0 .or. wSWP /= wSWP) then
          print*,"stocks less than zero or NaN", n
          print*,"steps_in_day",steps_in_day
          print*,stock_labile, stock_foliage
          print*,stock_stem,stock_roots
          print*,stock_litter,stock_soilOrgMatter
          print*,stock_storage_organ,stock_resp_auto
          print*,gpp_acm,nee_dalec
          print*,resp_auto,resp_h_soilOrgMatter,resp_h_litter
          print*,"pars",pars(1:33)
          print*,"fluxes",fluxes(n,1:16)
          print*,"DR",DR
          print*,"alloc_to_labile",alloc_to_labile,"remob",remob
          print*,"DR stuff",fT,fV,fP
          print*,"leaf","stem","rootlitter",litterfall_foliage,litterfall_stem,litterfall_roots
          print*,"daylength",dayl_hours,"VD",VD,"VDh",VDh
          print*,"avtemp",avtemp
          print*,"sown",sown,"emerged",emerged
          print*,"root_frac_intpol",root_frac_intpol
          print*,"npp_shoot",npp_shoot,"npp",npp
          print*,"RDR",RDR,"ts_length",ts_length
          stop
      endif

    end do ! no days loop

  end subroutine CARBON_MODEL_CROP
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine calc_pools_crops(DS_LRRT,LRRT)

    use CARBON_MODEL_MOD, only: lai, doy

    ! Allocated GPP to NPP and various carbon pools. Based !
    ! this on physiological responses to temperature       !
    ! vernalisation, and photoperiod.                      !

    implicit none

    ! arguments
    double precision, dimension(:), intent(inout) :: DS_LRRT, & !
                                                        LRRT    !

    ! local variables
    double precision :: decomp_efficency

    ! turnover rate of fine roots is now equal to the
    ! loss rate of roots (Penning de Vries, 1989)..
    turnover_rate_roots = interpolate( DS , DS_LRRT , LRRT , 5 ) / 24d0

    ! if sown turn on labile / seed turnover for growth
    if ( sown ) then
        ! turnover on
        turnover_labile_switch = 1
    else
        ! turnover off
        turnover_labile_switch = 0
    endif

    ! Initialise..
    resp_cost_foliage_to_labile = 0d0 ; yield = 0d0 ; BM_EX = 0d0

    ! respiratory cost of C transfer from labile pool to short-term pool (NPP) (gC.m-2.t-1)
    resp_cost_labile_to_foliage = turnover_rate_labile * resp_cost_labile_trans * resp_rate &
                                * ts_length * dble(turnover_labile_switch)
    resp_cost_labile_to_foliage = stock_labile * min(1d0,resp_cost_labile_to_foliage)

    ! allocation flux from labile C pool to NPP (gC.m-2.t-1)
    alloc_from_labile = turnover_rate_labile * ( 1d0 - resp_cost_labile_trans ) * resp_rate &
                      * ts_length * dble(turnover_labile_switch)
    alloc_from_labile = stock_labile * min(1d0,alloc_from_labile)

    ! When GPP is higher than seed C content, remaining seed carbon enters litter
    ! C pool, as seedlings do not fully exhaust their seed (P. de Vries p 48)
    if ( ( gpp_acm .gt. alloc_from_labile ) .and. ( use_seed_labile ) ) then
        stock_litter = stock_litter + stock_labile
        stock_labile = 0d0
        use_seed_labile = .false.
    endif

    ! NPP as a fraction of GPP (1-.32=.68 or 68%) + allocation..
    npp = ( 1d0 - frac_GPP_resp_auto ) * gpp_acm + alloc_from_labile
    ! from labile pool; = SHORT-TERM POOL

    root_frac_intpol  = max(0d0,min(1d0,root_frac_intpol))
    alloc_to_roots    = root_frac_intpol * npp         !
    shoot_frac_intpol = 1d0 - root_frac_intpol    !
    npp_shoot         = npp - alloc_to_roots           ! NPP remaining after root growth==SHOOT fraction
    alloc_to_foliage  = fol_frac_intpol  * npp_shoot   !
    alloc_to_stem     = stem_frac_intpol * npp_shoot   !
    alloc_to_storage_organ = max(0d0,npp_shoot - alloc_to_foliage - alloc_to_stem)
    if ( alloc_to_storage_organ > 0d0 ) then  ! allocation flux to storage organ limited by maximum growth rate
        gso_max  = ( stock_storage_organ + 0.5d0 ) * rel_gso_max / steps_in_day
        alloc_to_storage_organ = min( alloc_to_storage_organ , gso_max )
        if ( sown ) then
           alloc_to_labile = ( npp_shoot - alloc_to_foliage - alloc_to_stem - alloc_to_storage_organ ) &
                           * ( 1d0 - resp_cost_labile_trans )
           resp_cost_foliage_to_labile =  ( npp_shoot - alloc_to_foliage - alloc_to_stem - alloc_to_storage_organ ) &
                                       * resp_cost_labile_trans
        else
          alloc_to_labile             = 0d0
          resp_cost_foliage_to_labile = 0d0
        endif
    endif
    mean_alloc_to_storage_organ = mean_alloc_to_storage_organ + alloc_to_storage_organ

    ! set switches to (de)activate leaf, root and stem remobliization
    if ( step_of_day == steps_in_day ) then
        mean_alloc_to_storage_organ = mean_alloc_to_storage_organ / steps_in_day
        raso_old = raso
        ! running average of growth rate of storage organ..
        raso = ( mean_alloc_to_storage_organ + mean_alloc_to_storage_organ_old ) * 0.5d0
        max_raso_old = max_raso
        max_raso = max( raso , max_raso_old )
        ! Stem remobilisation triggered once running average of storage organ growth declines
        ! Second part prevents premature remobilisation
        if ( ( raso < raso_old ) .and. &
              ( mean_alloc_to_storage_organ > ( mean_alloc_to_storage_organ_old + 0.5d0 ) / steps_in_day ) ) then
            stmob = 1
        else
            stmob = 0
        endif
    endif

    ! Code for calculating relative death rate of leaves (RDR) as a
    !  function of shading (RDRSH) or developmental stage (RDRT).

    ! GT 0 if LAI GT 4; 0. < RDRSH < RDRSHMAX (usually ~0.03)
    RDRSH = min( RDRSHMAX , max( 0d0 , RDRSHMAX * ( lai - LAICR ) / LAICR ) )
    if ( DS < 1d0 ) then
       RDRDV = 0d0
    else
       ! RDRDV dependant on DR and DS, values range typically between 0.02 <
       ! RDRDV < 0.25
!print*,"!! What RDRDV to use? !!"
!!$      RDRDV = DR /( max( 0.1 , 2. - DS ) )
!!$      RDRDV = RDRDV / 24. ! to get hourly senescence rate
       RDRDV = turnover_rate_foliage * ( 1d0 / ( ( max( 2d0 - DS , 0.1d0 ) ) * 8d0 ) ) ** 2
    ENDIF

    ! relative leaf death rate is the maximum value of the arguments RDRSH and
    ! RDRDV
    RDR = max( RDRSH , RDRDV )

    ! remobilization of foliar C and allocation to dead leaves pool (gC.m-2.t-1)
    litterfall_foliage = stock_foliage * min(1d0,ts_length * RDR)
    litterfall_stem    = stock_stem    * min(1d0,ts_length * DR * turnover_rate_stem * dble(stmob)) ! remobstem
    litterfall_roots   = stock_roots   * min(1d0,ts_length * turnover_rate_roots)

    ! remobilized C to NPP (from both leaves and stems) (gC.m-2.t-1)
    remob   = ( litterfall_foliage * 0.5d0 + litterfall_stem ) * ( 1d0 - resp_cost_labile_trans )
    ! respiratory cost of C transfer (conversion from starch to photosynthates) (gC.m-2.t-1)
    Raremob = ( litterfall_foliage * 0.5d0 + litterfall_stem ) * resp_cost_labile_trans

    ! for mass balance calculate the decompostion efficency
    decomp_efficency = decomposition_rate &
                     / (decomposition_rate+mineralisation_rate_litter)

    ! total litter decomposition
    decomposition = stock_litter * (decomposition_rate+mineralisation_rate_litter) * resp_rate * ts_length

    ! heterotrophic respiration component 1: mineralisation of litter C pool (gC.m-2.t-1)
    resp_h_litter = decomposition * (1d0 - decomp_efficency)
    ! heterotrophic respiration component 2:  mineralisation of organic matter C pool (gC.m-2.t-1)
    resp_h_soilOrgMatter = stock_soilOrgMatter * min(1d0,mineralisation_rate_soilOrgMatter * resp_rate * ts_length)

    ! decomposition of litter to soil organic matter (gC.m-2.t-1)
    decomposition = decomposition - resp_h_litter

    ! Recalculate Carbon Pools...

    stock_foliage       = max(0d0, stock_foliage + alloc_to_foliage - litterfall_foliage)
    stock_stem          = max(0d0, stock_stem + alloc_to_stem - litterfall_stem)
    stock_storage_organ = max(0d0, stock_storage_organ + alloc_to_storage_organ)
    stock_roots         = max(0d0, stock_roots         + alloc_to_roots   - litterfall_roots)
    stock_litter        = max(0d0, stock_litter + litterfall_roots - resp_h_litter - decomposition)
    stock_soilOrgMatter = max(0d0, stock_soilOrgMatter + decomposition    - resp_h_soilOrgMatter)
    stock_dead_foliage  = max(0d0, stock_dead_foliage  + litterfall_foliage * 0.5d0) ! remainder of litfol is remobilisedi
    stock_labile        = max(0d0, stock_labile + alloc_to_labile  - alloc_from_labile - resp_cost_labile_to_foliage + remob)

    ! respiratory pool: new photosynthates are added (gC.m-2.t-1)
    stock_resp_auto = stock_resp_auto + frac_GPP_resp_auto * gpp_acm
    ! autotrophic respiration; Ra (typically ~7% of respiratory pool) (gC.m-2.t-1)
    resp_auto = stock_resp_auto * min(1d0,turnover_rate_resp_auto * ts_length)
    ! respiratory pool reduced by Ra (amount of C respired by plant)
    stock_resp_auto = max(0d0, stock_resp_auto - resp_auto)
    ! respiratory cost of C transfer from labile pool to short-term pool added
    ! to yield total autotrophic respiration (gC.m-2.t-1)
    resp_auto = resp_auto + resp_cost_labile_to_foliage + resp_cost_foliage_to_labile + Raremob
    ! nee (gC.m-2.t-1)
    nee_dalec = (resp_auto + resp_h_litter + resp_h_soilOrgMatter) - gpp_acm

  end subroutine calc_pools_crops
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine carbon_alloc_fractions(DS_shoot,DS_root,fol_frac,stem_frac,root_frac)

    ! Determines carbon allocation fractions as a function !
    ! of developmental stage (DS).  Allocation fractions   !
    ! are from tables published in Penning de Vries (1989) !

    implicit none

    double precision, dimension(:), intent(inout) ::   DS_shoot, & !
                                                        DS_root, & !
                                                       fol_frac, & !
                                                      stem_frac, & !
                                                      root_frac    !

    ! local variables..
    double precision,dimension(:),allocatable :: frac_shoot, frac_root

    if ( sown ) then ! after sowing

       allocate( frac_shoot(size(DS_shoot)) , frac_root(size(DS_root)) )

       ! loop over three crop "organs": 1) foliage 2) stems 3) root
       ! not necessary for storage organs, as all remaining C is allocated to
       ! these

       ! use different input for foliage and stem fractions, as they are
       ! relative to
       ! the total shoot (or aboveground) allocation, root is relative to
       ! total plant
       ! (above- and belowground) allocation..

       ! leaf development stages and corresponding fractions..

       frac_shoot = fol_frac
       ! interpolate between PdV allocation values with reference to
       ! developmental stage (DS)..
       fol_frac_intpol = interpolate( DS , DS_shoot , frac_shoot , size(DS_shoot) )

       ! stem DS and fracs..
       frac_shoot = stem_frac
       stem_frac_intpol = interpolate( DS , DS_shoot , frac_shoot , size(DS_shoot) )

       ! root DS and fracs..
       frac_root = root_frac
       root_frac_intpol = interpolate( DS , DS_root , frac_root , size(DS_root) )

    endif ! after crop has been sown

  end subroutine carbon_alloc_fractions
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine development_stage(days_in_step)

    ! Based on modified Wang & Engel model (Streck et al., 2003), !
    ! but with only 2 sub-phases, vegetative and reproductive     !
    ! (i.e. only two different DRmax).   O. Sus, May 2010.        !

    implicit none

    ! agruments
    double precision :: days_in_step

    ! local variables..
    double precision ::  doptmin, & ! Difference between optimum and minimum temperature
                         dmaxmin, & ! Difference between maximum and minimum temperature
                          dttmin, & ! Difference between daiy average and minimum temperatures
                       doptmin_v, & ! Difference between optimum and minimum vernalization temperatures
                       dmaxmin_v, & ! Difference between maximum and minimum vernalization temperatures
                       dttmin_v     ! Difference between daily average and minimum vernalization temperatures

    doptmin   = topt   - tmin   ! difference between optimal and minimum cardinal temperatures
    dmaxmin   = tmax   - tmin   ! difference between maximum and minimum cardinal temperatures
    dttmin    = avtemp - tmin   ! difference between daily average and minimum cardinal temperatures
    doptmin_v = topt_v - tmin_v ! same as above,
    dmaxmin_v = tmax_v - tmin_v !       but for vernalization
    dttmin_v  = avtemp - tmin_v ! cardinal temperatures

    ! Calculation of developmental function values: vernalization (fV),
    ! temperature (fT) and
    ! photoperiod (fP) these values are multiplicative factors of DRmax (maximum
    ! developmental
    ! rate), each ranging between 0 (no development) and 1 (unrestricted
    ! development).

    ! Summation of vernalization days (VD), not before sowing and only if
    ! average temperature is within min and max cardinal temperatures..
    if ( ( avtemp > tmin_v ) .and. ( avtemp < tmax_v ) .and. sown ) then
        fV = vernalization( doptmin_v , dmaxmin_v , dttmin_v , days_in_step )
    endif

    ! Only calculate temperature coefficient if avtemp lies within (tmin,tmax)
    ! range.
    ! NOTE: (doptmin+1d0) < dmaxmin added to allow for EDC search period when "not
    ! allowed" parameter sets will be tried anyway
    if ( avtemp > tmin .and. avtemp < tmax .and. (doptmin+1d0) < dmaxmin ) then
        fT = temperature_impact( doptmin , dmaxmin , dttmin )
    else
        fT = 0d0
    endif

    ! calculation of photoperiod coefficient
    fP = photoperiod_impact( PHCR , PHSC )

    if ( emerged .and. ( DS < 2d0 ) ) then   ! sum up daily DR values between emergence and maturity (DS=2)

       if ( DS < 1d0 ) then  ! in the vegetative phase (before flowering):

          DR = DR_pre * fT * fP   ! DR is affected by temperature, photoperiod...

          if ( vernal_calcs ) DR = DR * fV ! ...and vernalization (for winter cereals)

          DS = DS + (DR * days_in_step)    ! developmental stage (DS), calculated as the sum of daily developmental rates

       else    ! in the reproductive phase (after flowering):

          DR = DR_post * fT   ! DR is affected only by temperature

          DS = DS + (DR * days_in_step)

       endif ! vegetative or reproductive phase

    endif ! emerged or not

  end subroutine development_stage
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine management_dates (stock_seed_labile,days_in_step)

    use CARBON_MODEL_MOD, only: doy

    ! This routine should be called at the end of each day of a crops  !
    ! simulation.  It checks whether we should plough/sow/harvest, and !
    ! during the growing establishes when the crop will emerge after   !
    ! sowing, based on heat accumulation (Phenological Heat Units).    !

    implicit none

    ! arguments
    double precision, intent(in) :: stock_seed_labile,days_in_step
    ! local variables
    double precision :: tmp
    logical :: plough_sanity,sow_sanity,harvest_sanity

    ! reset
    plough_sanity = .false. ; sow_sanity = .false. ; harvest_sanity = .false.

    ! spring crop
    if (sow_day < harvest_day .and. nint(doy) < harvest_day) sow_sanity = .true.
    if (plough_day < harvest_day .and. nint(doy) < harvest_day) plough_sanity = .true.
    if (harvest_day > sow_day) harvest_sanity = .true.
    ! winter crops
    if (sow_day > harvest_day) sow_sanity = .true.
    if (plough_day > harvest_day) plough_sanity = .true.
    if (harvest_day < plough_day .and. nint(doy) < plough_day) harvest_sanity = .true.

    if ( .not. sown ) then

      ! fresh field...

      if ( plough_sanity .and. .not.ploughed .and. nint(doy) >= plough_day ) then
        ! the field needs ploughing..
        call plough

      elseif ( sow_sanity .and. nint(doy) >= sow_day ) then

        ! ensure that the field has indeed been ploughed
        if (.not.ploughed) call plough
        ! the field needs sowing..
        sown = .true.

        ! this switch controls whether the labile carbon within the seed is used
        ! for growth
        use_seed_labile = .true.
        stock_labile = stock_seed_labile

      endif ! plough or sow?

    else

      ! crop in field..

      ! calculate when crop emerges..
      if ( .not. emerged ) then

         ! estimate emergence date based on the accumulated phenological heat
         ! units (PHU)
         ! where PHU is the (positive) heat over tmin..
         tmp = max( avtemp - tmin , 0d0 )*days_in_step
         PHU = PHU + tmp

         ! set the development stage and emergence..
         if ( PHU >= PHUem ) then
           emerged = .true.
           DS = 0d0
         else
           emerged = .false.
           DS = -1d0
         endif

      endif ! emerged or not

      ! note that in this case harvest day has been fixed relative to the sow
      ! day
      if ( harvest_sanity .and. nint(doy) >= harvest_day) then
         ! the field needs harvesting..
         call harvest
      endif

    endif ! sown or not

  end subroutine management_dates
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  ! PROCEDURES BELOW ARE PRIVATE, IE THEIR USE IS LIMITED TO THIS MODULE
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine harvest

    implicit none

    ! shoot biomass..
    Cshoot = stock_foliage + stock_stem + stock_storage_organ + stock_labile

    ! determine harvest index..
    HI = stock_storage_organ / Cshoot

    ! the stuff we actually want from the harvest...
    yield = stock_storage_organ !+ stock_stem * ( 1d0 - st_res )

    ! the biomass that is harvested in addition to the storage-organ..
    BM_EX  = stock_foliage * ( 1d0 - lv_res )          &
              + stock_stem * ( 1d0 - st_res )          &
               + stock_dead_foliage * ( 1d0 - lv_res ) &
                + stock_labile

    ! what's left (will fall to the ground)..
    stock_litter  = stock_litter                     &
                    + stock_resp_auto                &
                     + stock_foliage * lv_res        &
                      + stock_stem * st_res          &
                       + stock_dead_foliage * lv_res

    ! empty the plant stocks..
    stock_storage_organ = 0d0
    stock_foliage       = 0d0
    stock_stem          = 0d0
    stock_dead_foliage  = 0d0
    stock_labile        = 0d0
    stock_resp_auto     = 0d0

    ! roots stay in ground and slowly decompose (until/unless the field is
    ! ploughed)

    ! reset logical variables..
    sown    = .false.
    emerged = .false.
    ploughed = .false.
    DS = -1d0 ; fV = 0d0 ; fT = 0d0 ; fP = 0d0 ; VD = 0d0

  end subroutine harvest
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  double precision function photoperiod_impact( PH_crit , PH_sens )

    use CARBON_MODEL_MOD, only: dayl_hours

    ! Function to determine the coefficient for !
    ! photoperiod impact on developmental rate. !
    ! From Streck et al., 2003                  !

    implicit none

    ! arguments..
    double precision,intent(in) :: PH_crit, & ! critical photoperiod below which no development occurs
                                   PH_sens    ! photoperiod sensitivity

    photoperiod_impact = max(0d0, 1d0 - exp ( - PH_Sens * ( dayl_hours - PH_crit ) ))

  end function photoperiod_impact
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine plough

    ! this s/r will reset various carbon pools, to mimic the effect of the
    ! farmer ploughing. !

    implicit none

    ! Move all plant stocks into the litter pool.
    ! ( many of these should already be empty after the harvest, )
    ! ( e.g. the stocks for labile, foliage, storage-organ stem. )
    stock_litter        = stock_litter + stock_dead_foliage &
                          + stock_foliage + stock_labile    &
                           + stock_roots + stock_stem       &
                            + stock_storage_organ
    stock_dead_foliage  = 0d0
    stock_foliage       = 0d0
    stock_labile        = 0d0
    stock_roots         = 0d0
    stock_stem          = 0d0
    stock_storage_organ = 0d0

    ! Reset the development stage & phenological heat units..
    ploughed = .true. ; DS = -1d0 ; PHU = 0d0
    max_raso = 0d0 ; raso = 0d0 ; max_raso_old = 0d0 ; raso_old = 0d0
    mean_alloc_to_storage_organ_old = 0d0 ; mean_alloc_to_storage_organ = 0d0

  end subroutine plough
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  double precision function temperature_impact( doptmin , dmaxmin , dttmin )

    ! Function to determine the coefficent for  !
    ! temperature impact on developmental rate. !
    ! From Streck et al., 2003.                 !

    implicit none

    ! arguments..
    double precision,intent(in) :: doptmin , dmaxmin , dttmin   ! temperature differences

    ! local variables..
    double precision :: a , nmr , dnr

    a   = log( 2.d0 ) / ( log( ( dmaxmin ) / doptmin ) )

    nmr = 2.d0 * ( ( dttmin ) ** a ) * ( doptmin ** a ) - ( ( dttmin ) ** ( 2.d0 * a ) )

    dnr = doptmin ** ( 2.d0 * a )

    temperature_impact = nmr / dnr

  end function temperature_impact
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  double precision function vernalization( doptmin_v , dmaxmin_v , dttmin_v , days_in_step )

    ! Function to determine the coefficent for vernalization !
    ! impact on developmental rate. See Streck et al., 2003. !

    implicit none

    ! arguments..
    double precision,intent(in) :: dmaxmin_v , doptmin_v , dttmin_v & ! temperature differences
                                  ,days_in_step

    ! local variables..
    double precision :: a , dnr , fvn , nmr

    a   = log( 2.d0 ) / ( log( ( dmaxmin_v ) / doptmin_v ) )
    nmr = 2.d0 * ( ( dttmin_v ) ** a ) * ( doptmin_v ** a ) - ( ( dttmin_v ) ** (2.d0 * a ) )
    dnr = doptmin_v ** ( 2.d0 * a )
    fvn = nmr / dnr

    VD = VD + (fvn*days_in_step)

    ! final output value..
    vernalization = max( 0d0 , min( 1d0 , ( VD ** 5 ) / ( ( VDh ** 5 ) + (VD ** 5 ) ) ) )

  end function vernalization
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  double precision function interpolate( x , reference_x , reference_y , row )

    ! Interpolation function.                    !
    ! x is input value, interpol is output value !
    ! reference_x/y are reference input data.    !

    implicit none

    ! arguments..
    integer, intent(in)                          :: row
    double precision, intent(in)                 :: x
    double precision, dimension(row), intent(in) :: reference_x , reference_y

    ! local variables..
    integer::i

    ! provide initial value
    interpolate = -9999d0
    do i = 1 , row

       if ( x .le. reference_x(1) ) then
          interpolate = reference_y(1)
          exit
       endif

       ! cycling means growth rate remains constant between DS levels
       if ( ( x .gt. reference_x(i) ) .and. ( i .lt. row ) ) cycle

       if ( x .eq. reference_x(i) ) then
          interpolate = reference_y(i)
          exit
       endif

       if ( x .lt. reference_x(i) ) then
          interpolate = reference_y(i-1) + ( x - reference_x(i-1) ) &
                       * ( reference_y(i) - reference_y(i-1) )      &
                       / ( reference_x(i) - reference_x(i-1) )
          exit
       else
          interpolate = reference_y(row)
       endif

    enddo

    ! explicit return to ser
    return

  end function interpolate
!
!--------------------------------------------------------------------
!
end module CARBON_MODEL_CROP_MOD
