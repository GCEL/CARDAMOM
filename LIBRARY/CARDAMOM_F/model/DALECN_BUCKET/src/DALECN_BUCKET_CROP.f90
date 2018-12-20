
module CARBON_MODEL_CROP_MOD

implicit none

! make all private
private

! explicit publics
public :: CARBON_MODEL_CROP    &
         ,resp_rate_temp_coeff &
         ,ts_length            &
         ,sec_in_hour          &
         ,DS

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
  double precision, parameter :: LAICR = 4.5d0
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
                       ,stem_frac,root_frac,DS_LRLV,LRLV,DS_LRRT,LRRT, DS_array, leaf_N_mass)

    use CARBON_MODEL_MOD, only: arrhenious,acm_gpp,calculate_transpiration,calculate_soil_evaporation     &
                               ,calculate_wetcanopy_evaporation,acm_meteorological_constants,acm_albedo_gc&
                               ,calculate_shortwave_balance,calculate_longwave_isothermal                 &
                               ,calculate_update_soil_water,calculate_daylength,drythick,vsmall           &
                               ,co2_half_saturation,co2_compensation_point,freeze,co2comp_saturation      &
                               ,co2comp_half_sat_conc,kc_saturation,kc_half_sat_conc,dble_zero,dble_one   &
                               ,calculate_Rtot,calculate_aerodynamic_conductance,saxton_parameters        &
                               ,initialise_soils,min_drythick,soil_frac_clay,soil_frac_sand,dayl_hours    &
                               ,seconds_per_day,dayl_seconds,dayl_seconds_1,seconds_per_step,root_biomass &
                               ,top_soil_depth,mid_soil_depth,root_reach,min_root,max_depth,root_k,previous_depth &
                               ,min_layer,wSWP,SWP,SWP_initial,deltat_1,water_flux,soilwatermm,wSWP_time  &
                               ,layer_thickness,minlwp,soil_waterfrac,soil_waterfrac_initial              &
                               ,porosity,porosity_initial,field_capacity,field_capacity_initial,meant     &
                               ,meant_K,stomatal_conductance,avN,iWUE,NUE,pn_max_temp,pn_opt_temp         &
                               ,pn_kurtosis,e0,co2_half_sat,co2_comp_point,max_lai_lwrad_transmitted      &
                               ,lai_half_lwrad_transmitted,max_lai_nir_reflection,lai_half_nir_reflection &
                               ,max_lai_par_reflection,lai_half_par_reflection,max_lai_par_transmitted    &
                               ,lai_half_par_transmitted,max_lai_nir_transmitted,lai_half_nir_transmitted &
                               ,max_lai_lwrad_reflected,lai_half_lwrad_reflected,soil_swrad_absorption    &
                               ,max_lai_lwrad_release,lai_half_lwrad_release,mint,maxt,swrad,co2,doy,leafT&
                               ,rainfall,wind_spd,vpd_pa,lai,days_per_step,days_per_step_1,canopy_storage &
                               ,intercepted_rainfall,snow_storage,snow_melt,airt_zero_fraction,snowfall   &
                               ,update_soil_initial_conditions,aerodynamic_conductance

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
                         ,lat                           & ! site latitude (degrees)
                         ,leaf_N_mass(1,nodays)

    double precision, dimension(:), intent(inout) ::          DS_shoot, & !
                                                               DS_root, & !
                                                              fol_frac, & !
                                                             stem_frac, & !
                                                             root_frac, & !
                                                               DS_LRLV, & !
                                                                  LRLV, & !
                                                               DS_LRRT, & !
                                                                  LRRT, & !
                                                                  DS_array

    double precision, dimension(nodays), intent(inout) :: lai_out & ! leaf area index
                                               ,GPP_out & ! Gross primary productivity
                                               ,NEE_out   ! net ecosystem exchange of CO2

    double precision, dimension((nodays+1),nopools), intent(inout) :: POOLS ! vector of ecosystem pools

    double precision, dimension(nodays,nofluxes), intent(inout) :: FLUXES ! vector of ecosystem fluxes

    ! declare local variables
    double precision :: airt_weighting(3) &
                                      ,ET & ! Evapotranspiration (kg.m-2.day-1)
                        ,pot_actual_ratio &
                           ,transpiration &
                         ,soilevaporation &
                          ,wetcanopy_evap &
                        ,snow_sublimation &
                                 ,deltaWP & !
                                    ,Rtot &
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
    lai_out(1:nodays) = dble_zero ; NEE_out(1:nodays) = dble_zero ; GPP_out(1:nodays) = dble_zero
    FLUXES(1:nodays,1:nofluxes) = dble_zero ; POOLS(1:nodays,1:nopools) = dble_zero

    ! load ACM-GPP-ET parameters
    NUE                        = 1.182549d+01  ! Photosynthetic nitrogen use efficiency at optimum temperature (oC)
                                               ! ,unlimited by CO2, light and
                                               ! photoperiod (gC/gN/m2leaf/day)
    pn_max_temp                = 5.357174d+01  ! Maximum temperature for photosynthesis (oC)
    pn_opt_temp                = 3.137242d+01  ! Optimum temperature for photosynthesis (oC)
    pn_kurtosis                = 1.927458d-01  ! Kurtosis of photosynthesis temperature response
    e0                         = 5.875662d+00  ! Quantum yield gC/MJ/m2/day PAR
    max_lai_lwrad_transmitted  = 7.626683d-01  ! Max fractional reduction of LW from sky transmitted through canopy
    lai_half_lwrad_transmitted = 7.160363d-01  ! LAI at which canopy LW transmittance reduction = 50 %
    max_lai_nir_reflection     = 4.634860d-01  ! Max fraction of NIR reflected by canopy
    lai_half_nir_reflection    = 1.559148d+00  ! LAI at which canopy NIR reflected = 50 %
    minlwp                     =-1.996830d+00  ! minimum leaf water potential (MPa)
    max_lai_par_reflection     = 1.623013d-01  ! Max fraction of PAR reflected by canopy
    lai_half_par_reflection    = 1.114360d+00  ! LAI at which canopy PAR reflected = 50 %
    lai_half_lwrad_reflected   = 1.126214d+00  ! LAI at which 50 % LW is reflected back to sky
    iWUE                       = 1.602503d-06  ! Intrinsic water use efficiency (gC/m2leaf/day/mmolH2Ogs)
    soil_swrad_absorption      = 6.643079d-01  ! Fraction of SW rad absorbed by soil
    max_lai_par_transmitted    = 8.079519d-01  ! Max fractional reduction in PAR transmittance by canopy
    lai_half_par_transmitted   = 9.178784d-01  ! LAI at which PAR transmittance reduction = 50 %
    max_lai_nir_transmitted    = 8.289803d-01  ! Max fractional reduction in NIR transmittance by canopy
    lai_half_nir_transmitted   = 1.961831d+00  ! LAI at which NIR transmittance reduction = 50 %
    max_lai_lwrad_release      = 9.852855d-01  ! Max fraction of LW emitted (1-par) from canopy to be released
    lai_half_lwrad_release     = 7.535450d-01  ! LAI at which LW emitted from canopy to be released at 50 %
    max_lai_lwrad_reflected    = 1.955832d-02  ! LAI at which 50 % LW is reflected back to sky
    
    ! load some values                      ! AR
    !avN = 10d0**pars(11) !TLS 1 ! foliar N  ! AR
    deltaWP = minlwp !-2.0  ! leafWP-soilWP (i.e. -2-0)
    Rtot = dble_one ! totaly hydraulic resistance ! p12 from ACM recal (updated)

    ! plus ones being calibrated
    root_k = pars(36) ; max_depth = pars(37)

    ! length of time step in hours..
    ts_length = ((sum(deltat)/dble(nodays)) * sec_in_day) / sec_in_hour
    ! steps per day ; set step_of_day to steps_in_day as in CARDAMOM we will
    ! never be running less than daily time step
    steps_in_day = dble_one/(sum(deltat)/dble(nodays)) ; step_of_day = steps_in_day

    ! parameters from file
    decomposition_rate                = pars(1) / 24d0  ! decomposition rate (day->hr)
    frac_GPP_resp_auto                = pars(2)  ! fraction of GPP allocated to autotrophic carbon pool
    DR_pre                            = pars(3)  ! development rate coefficient DS (0->1)
    DR_post                           = pars(4)  ! development rate coefficient DS (1->2)
    turnover_rate_foliage             = (pars(5) / 24d0)*1.5 ! pars(5)  ! turnover_rate of foliage (day->hr)
    turnover_rate_stem                = pars(6) / 24d0 ! turnover rate of stem (day->hr)
    RDRSHMAX                          = pars(7) / 24d0 ! maximum rate of foliar turnover due to self shading (day->hr)
    VDh                               = pars(8)*0.8  ! effective vernalisation days when plants are 50 % vernalised
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
    topt                              = (pars(28)-273.15d0) ! optimum temperature for development
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
        use_seed_labile = .true.
        emerged         = .false.

        ! pair incoming variables to local module levels
        ! finally set some initial conditions
        avtemp = dble_zero
        yield = dble_zero
        DS = -dble_one
        DR = dble_zero
        fV = dble_zero ; fT = dble_zero ; fP = dble_zero
        mean_alloc_to_storage_organ = dble_zero
        mean_alloc_to_storage_organ_old = dble_zero
        PHU = dble_zero
        VD = dble_zero
        BM_EX = dble_zero
        HI = dble_zero
        stock_dead_foliage = dble_zero
        alloc_to_labile = dble_zero
        stmob = 0
        max_raso = dble_zero
        raso = dble_zero
        RDRDV = dble_zero
        max_raso_old = dble_zero
        raso_old  = dble_zero

        ! SHOULD TURN THIS INTO A SUBROUTINE CALL AS COMMON TO BOTH DEFAULT AND CROPS
        if (.not.allocated(deltat_1)) then

           allocate(deltat_1(nodays),soilwatermm(nodays),wSWP_time(nodays), &
                    co2_compensation_point(nodays),co2_half_saturation(nodays))
           do n = 1, nodays
              ! Temperature adjustments for Michaelis-Menten coefficients
              ! for CO2 (kc) and O2 (ko) and CO2 compensation point.
              co2_compensation_point(n) = arrhenious(co2comp_saturation,co2comp_half_sat_conc,met(3,n))
              co2_half_saturation(n) = arrhenious(kc_saturation,kc_half_sat_conc,met(3,n))
           end do
           deltat_1 = deltat**(-dble_one)
           ! zero variables not done elsewhere
           water_flux = dble_zero
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

           water_flux = dble_zero
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
        mint = met(2,1)  ! minimum temperature (oC)
        maxt = met(3,1)  ! maximum temperature (oC)
        meant = (maxt+mint)*0.5d0  ! mean air temperature (oC)
        meant_K = meant + freeze

        ! zero evapotranspiration for beginning
        ET = dble_zero
        ! initialise root reach based on initial conditions
        root_biomass = max(min_root,POOLS(1,3)*2d0)
        root_reach = max_depth * root_biomass / (root_k + root_biomass)
        ! Determine initial soil layer thickness
        layer_thickness(1) = top_soil_depth ; layer_thickness(2) = mid_soil_depth
        layer_thickness(3) = max(min_layer,root_reach-sum(layer_thickness(1:2)))
        layer_thickness(4) = max_depth - sum(layer_thickness(1:3))
        layer_thickness(5) = top_soil_depth
        previous_depth = max(top_soil_depth,root_reach)
        ! needed to initialise soils
        call calculate_Rtot(Rtot)
        ! used to initialise soils
        call calculate_update_soil_water(dble_zero,dble_zero,dble_zero,ET) ! assume no evap or rainfall
        ! store soil water content of the rooting zone (mm)
        POOLS(1,8) = 1d3*soil_waterfrac(1)*layer_thickness(1)

    else

        ! load ET from memory
        ET = FLUXES(start-1,19)

    endif ! start == 1

    ! reset values
    intercepted_rainfall = dble_zero ; canopy_storage = dble_zero ; snow_storage = dble_zero

    infi = dble_zero

    !
    ! Begin looping through each time step
    !

    do n = start, finish

      
       
      !avN = 10d0**pars(11) !TLS 1 ! foliar N  ! AR       
      avN = 10d0**((leaf_N_mass(1,n))*0.1) !TLS 1 ! foliar N  ! AR       
       
      !!!!!!!!!!
      ! assign drivers and update some prognostic variables
      !!!!!!!!!!
             
      ! Incoming drivers
      mint = met(2,n)  ! minimum temperature (oC)
      maxt = met(3,n)  ! maximum temperature (oC)

      !print*, met(6,n)
      
      leafT = maxt     ! initial canopy temperature (oC)
      swrad = met(4,n) ! incoming short wave radiation (MJ/m2/day)
      co2 = met(5,n)   ! CO2 (ppm)
      doy = ceiling(met(6,n)-(deltat(n)*0.5d0))   ! Day of year
      rainfall = max(dble_zero,met(7,n)) ! rainfall (kgH2O/m2/s)
      meant = (maxt+mint) * 0.5d0   ! mean air temperature (oC)
      meant_K = meant + freeze
      airt_zero_fraction = (maxt-dble_zero) / (maxt-mint) ! fraction of temperture period above freezing
      wind_spd = met(15,n) ! wind speed (m/s)
      vpd_pa = met(16,n)  ! Vapour pressure deficit (Pa)

      ! states needed for module variables
      lai_out(n) = POOLS(n,2)/LCA
      lai = lai_out(n) ! leaf area index (m2/m2)

      ! Temperature adjustments for Michaelis-Menten coefficients
      ! for CO2 (kc) and O2 (ko) and CO2 compensation point
      ! See McMurtrie et al., (1992) Australian Journal of Botany, vol 40, 657-677
      co2_half_sat   = co2_half_saturation(n)
      co2_comp_point = co2_compensation_point(n)

      ! calculate daylength in hours and seconds
      call calculate_daylength((doy-(deltat(n)*0.5d0)),lat)
      ! extract timing related values
      dayl_seconds_1 = dayl_seconds ** (-dble_one)
      seconds_per_step = seconds_per_day * deltat(n)
      days_per_step = deltat(n)
      days_per_step_1 = deltat_1(n)

      ! snowing or not...?
      snow_melt = dble_zero ; snowfall = dble_zero
      if (mint < dble_zero .and. maxt > dble_zero) then
          ! if minimum temperature is below freezing point then we weight the
          ! rainfall into snow or rain based on proportion of temperature below
          ! freezing
          snowfall = dble_one - airt_zero_fraction
          snowfall = rainfall * snowfall ; rainfall = rainfall - snowfall
          ! Add rainfall to the snowpack and clear rainfall variable
          snow_storage = snow_storage + (snowfall*seconds_per_step)

          ! Also melt some of the snow
          snow_melt = airt_zero_fraction
          ! otherwise we assume snow is melting at 10 % per day light hour
          snow_melt = min(snow_storage, snow_melt * snow_storage * dayl_hours * 0.1d0 * deltat(n))
          snow_storage = snow_storage - snow_melt
      elseif (maxt < dble_zero) then
          ! if whole day is below freezing then we should assume that all
          ! precipitation is snowfall
          snowfall = rainfall ; rainfall = dble_zero
          ! Add rainfall to the snowpack and clear rainfall variable
          snow_storage = snow_storage + (snowfall*seconds_per_step)
      else if (mint > dble_zero) then
          ! otherwise we assume snow is melting at 10 % per day light hour
          snow_melt = min(snow_storage, snow_storage * dayl_hours * 0.1d0 * deltat(n))
          snow_storage = snow_storage - snow_melt
      end if

      ! calculate the minimum soil & root hydraulic resistance based on total
      ! fine root mass ! *2*2 => *RS*C->Bio
      root_biomass = max(min_root,POOLS(n,3)*2d0)
      ! estimate drythick for the current step
      drythick = max(min_drythick, top_soil_depth * min(dble_one,dble_one - (soil_waterfrac(1) / porosity(1))))
      call calculate_Rtot(Rtot)

      ! Pass wSWP to output variable and update deltaWP between minlwp and
      ! current weighted soil WP
      wSWP_time(n) = wSWP ; deltaWP = min(dble_zero,minlwp-wSWP)

      ! calculate some temperature dependent meteorologial properties
      call acm_meteorological_constants(maxt)
      ! calculate radiation absorption and estimate stomatal conductance
      call calculate_aerodynamic_conductance
      call calculate_shortwave_balance ; call calculate_longwave_isothermal(leafT,maxt)
      call acm_albedo_gc(abs(deltaWP),Rtot)

      ! reallocate for crop model timings
      doy = met(6,n)

      ! GPP (gC.m-2.day-1)
      if (lai > vsmall .and. stomatal_conductance > vsmall) then
         GPP_out(n) = max(dble_zero,acm_gpp(stomatal_conductance))
         
         ! Canopy transpiration (kgH2O/m2/day)
         call calculate_transpiration(transpiration)
      else
         GPP_out(n) = dble_zero
         transpiration = dble_zero
      endif
      ! load GPP for crop model daily rate to total
      gpp_acm = GPP_out(n) * deltat(n)

      ! Canopy intercepted rainfall evaporation (kgH2O/m2/day)
      call calculate_wetcanopy_evaporation(wetcanopy_evap,pot_actual_ratio,canopy_storage,transpiration)
      ! Soil surface (kgH2O.m-2.day-1)
      call calculate_soil_evaporation(soilevaporation)
      ! restrict transpiration to positive only
      transpiration = max(dble_zero,transpiration)

      ! if snow present assume that soilevaporation is sublimation of soil first
      snow_sublimation = dble_zero
      if (snow_storage > dble_zero) then
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
      mean_alloc_to_storage_organ     = dble_zero
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
      resp_rate = dble_zero
      resp_rate = resp_rate + ((0.5d0 * exp( resp_rate_temp_coeff * met(3,n) )) * airt_weighting(1))
      resp_rate = resp_rate + ((0.5d0 * exp( resp_rate_temp_coeff * avtemp   )) * airt_weighting(2))
      resp_rate = resp_rate + ((0.5d0 * exp( resp_rate_temp_coeff * met(2,n) )) * airt_weighting(3))
      !resp_rate = 0.5 * exp( resp_rate_temp_coeff * avtemp )

      ! determine development stage (DS)
      call development_stage(deltat(n))
      DS_array(n) = DS
      
      
      ! determine the carbon partitioning based on development stage
      call carbon_alloc_fractions(DS_shoot,DS_root,fol_frac,stem_frac,root_frac)

      ! begin carbon allocation for crops
      call calc_pools_crops(DS_LRRT,LRRT)

      !print*, stock_storage_organ
      
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

      
      
!!$      do nxp = 1, nopools
!!$         if (POOLS(n+1,nxp) /= POOLS(n+1,nxp) .or. POOLS(n+1,nxp) < dble_zero) then
!!$             print*,"step",n,"FLUXES",nxp
!!$             print*,"met",met(:,n)
!!$             print*,"POOLS",POOLS(n,:)
!!$             print*,"FLUXES",FLUXES(n,:)
!!$             print*,"POOLS+1",POOLS(n+1,:)
!!$             print*,"wSWP",wSWP
!!$             print*,"waterfrac",soil_waterfrac
!!$             print*,"steps_in_day",steps_in_day
!!$             print*,stock_labile, stock_foliage
!!$             print*,stock_stem,stock_roots
!!$             print*,stock_litter,stock_soilOrgMatter
!!$             print*,stock_storage_organ,stock_resp_auto
!!$             print*,gpp_acm,nee_dalec
!!$             print*,resp_auto,resp_h_soilOrgMatter,resp_h_litter
!!$             print*,"pars",pars
!!$             print*,"DR",DR
!!$             print*,"alloc_to_labile",alloc_to_labile,"remob",remob
!!$             print*,"DR stuff",fT,fV,fP
!!$             print*,"leaf","stem","rootlitter",litterfall_foliage,litterfall_stem,litterfall_roots
!!$             print*,"daylength",dayl_hours,"VD",VD,"VDh",VDh
!!$             print*,"avtemp",avtemp
!!$             print*,"sown",sown,"emerged",emerged
!!$             print*,"root_frac_intpol",root_frac_intpol
!!$             print*,"npp_shoot",npp_shoot,"npp",npp
!!$             print*,"RDR",RDR,"ts_length",ts_length
!!$             stop
!!$         endif
!!$      enddo
!!$
!!$      do nxp = 1, nofluxes
!!$         if (nxp /= 19) then
!!$            if (FLUXES(n,nxp) /= FLUXES(n,nxp) .or. FLUXES(n,nxp) < dble_zero) then
!!$                 print*,"Special: step",n,"FLUXES",nxp
!!$                 print*,"met",met(:,n)
!!$                 print*,"POOLS",POOLS(n,:)
!!$                 print*,"FLUXES",FLUXES(n,:)
!!$                 print*,"POOLS+1",POOLS(n+1,:)
!!$                 print*,"wSWP",wSWP
!!$                 print*,"waterfrac",soil_waterfrac
!!$                 print*,"steps_in_day",steps_in_day
!!$                 print*,stock_labile, stock_foliage
!!$                 print*,stock_stem,stock_roots
!!$                 print*,stock_litter,stock_soilOrgMatter
!!$                 print*,stock_storage_organ,stock_resp_auto
!!$                 print*,gpp_acm,nee_dalec
!!$                 print*,resp_auto,resp_h_soilOrgMatter,resp_h_litter
!!$                 print*,"pars",pars
!!$                 print*,"DR",DR
!!$                 print*,"alloc_to_labile",alloc_to_labile,"remob",remob
!!$                 print*,"DR stuff",fT,fV,fP
!!$                 print*,"leaf","stem","rootlitter",litterfall_foliage,litterfall_stem,litterfall_roots
!!$                 print*,"daylength",dayl_hours,"VD",VD,"VDh",VDh
!!$                 print*,"avtemp",avtemp
!!$                 print*,"sown",sown,"emerged",emerged
!!$                 print*,"root_frac_intpol",root_frac_intpol
!!$                 print*,"npp_shoot",npp_shoot,"npp",npp
!!$                 print*,"RDR",RDR,"ts_length",ts_length
!!$                 stop
!!$            end if
!!$         else
!!$             if (FLUXES(n,nxp) /= FLUXES(n,nxp)) then
!!$                 print*,"Default: step",n,"FLUXES",nxp
!!$                 print*,"met",met(:,n)
!!$                 print*,"POOLS",POOLS(n,:)
!!$                 print*,"FLUXES",FLUXES(n,:)
!!$                 print*,"POOLS+1",POOLS(n+1,:)
!!$                 print*,"wSWP",wSWP
!!$                 print*,"waterfrac",soil_waterfrac
!!$                 print*,"steps_in_day",steps_in_day
!!$                 print*,stock_labile, stock_foliage
!!$                 print*,stock_stem,stock_roots
!!$                 print*,stock_litter,stock_soilOrgMatter
!!$                 print*,stock_storage_organ,stock_resp_auto
!!$                 print*,gpp_acm,nee_dalec
!!$                 print*,resp_auto,resp_h_soilOrgMatter,resp_h_litter
!!$                 print*,"pars",pars
!!$                 print*,"DR",DR
!!$                 print*,"alloc_to_labile",alloc_to_labile,"remob",remob
!!$                 print*,"DR stuff",fT,fV,fP
!!$                 print*,"leaf","stem","rootlitter",litterfall_foliage,litterfall_stem,litterfall_roots
!!$                 print*,"daylength",dayl_hours,"VD",VD,"VDh",VDh
!!$                 print*,"avtemp",avtemp
!!$                 print*,"sown",sown,"emerged",emerged
!!$                 print*,"root_frac_intpol",root_frac_intpol
!!$                 print*,"npp_shoot",npp_shoot,"npp",npp
!!$                 print*,"RDR",RDR,"ts_length",ts_length
!!$                 stop
!!$             endif
!!$         end if
!!$      enddo

!!$      if (stock_labile < dble_zero .or. stock_foliage < dble_zero .or. stock_stem < dble_zero .or. &
!!$          stock_roots < dble_zero .or. stock_litter < dble_zero .or. stock_soilOrgMatter < dble_zero .or. &
!!$          stock_storage_organ < dble_zero .or. stock_resp_auto < dble_zero .or. &
!!$          stock_labile /= stock_labile .or. stock_foliage /= stock_foliage .or. &
!!$          stock_stem /= stock_stem .or. &
!!$          stock_roots /= stock_roots .or. stock_litter /= stock_litter .or. &
!!$          stock_soilOrgMatter /= stock_soilOrgMatter .or. &
!!$          stock_storage_organ /= stock_storage_organ .or. &
!!$          stock_resp_auto /= stock_resp_auto .or.  &
!!$          gpp_acm < dble_zero .or. gpp_acm /= gpp_acm .or. resp_rate < dble_zero .or. &
!!$          resp_rate /= resp_rate .or. decomposition < dble_zero .or. alloc_from_labile < dble_zero .or. &
!!$          resp_cost_labile_to_foliage < dble_zero .or. alloc_to_foliage < dble_zero .or. &
!!$          alloc_to_stem < dble_zero .or. alloc_to_roots < dble_zero .or. remob < dble_zero .or. &
!!$          alloc_from_labile < dble_zero .or. resp_cost_labile_to_foliage < dble_zero .or. wSWP /= wSWP) then
!!$          print*,"stocks less than zero or NaN", n
!!$          print*,"steps_in_day",steps_in_day
!!$          print*,stock_labile, stock_foliage
!!$          print*,stock_stem,stock_roots
!!$          print*,stock_litter,stock_soilOrgMatter
!!$          print*,stock_storage_organ,stock_resp_auto
!!$          print*,gpp_acm,nee_dalec
!!$          print*,resp_auto,resp_h_soilOrgMatter,resp_h_litter
!!$          print*,"pars",pars(1:33)
!!$          print*,"fluxes",fluxes(n,1:16)
!!$          print*,"DR",DR
!!$          print*,"alloc_to_labile",alloc_to_labile,"remob",remob
!!$          print*,"DR stuff",fT,fV,fP
!!$          print*,"leaf","stem","rootlitter",litterfall_foliage,litterfall_stem,litterfall_roots
!!$          print*,"daylength",dayl_hours,"VD",VD,"VDh",VDh
!!$          print*,"avtemp",avtemp
!!$          print*,"sown",sown,"emerged",emerged
!!$          print*,"root_frac_intpol",root_frac_intpol
!!$          print*,"npp_shoot",npp_shoot,"npp",npp
!!$          print*,"RDR",RDR,"ts_length",ts_length
!!$          stop
!!$      endif

    end do ! no days loop

  end subroutine CARBON_MODEL_CROP
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine calc_pools_crops(DS_LRRT,LRRT)

    use CARBON_MODEL_MOD, only: dble_zero, dble_one, lai, doy

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
    resp_cost_foliage_to_labile = dble_zero ; yield = dble_zero ; BM_EX = dble_zero

    ! respiratory cost of C transfer from labile pool to short-term pool (NPP) (gC.m-2.t-1)
    resp_cost_labile_to_foliage = turnover_rate_labile * resp_cost_labile_trans * resp_rate &
                                * ts_length * dble(turnover_labile_switch)
    resp_cost_labile_to_foliage = stock_labile * min(dble_one,resp_cost_labile_to_foliage)

    ! allocation flux from labile C pool to NPP (gC.m-2.t-1)
    alloc_from_labile = turnover_rate_labile * ( dble_one - resp_cost_labile_trans ) * resp_rate &
                      * ts_length * dble(turnover_labile_switch)
    alloc_from_labile = stock_labile * min(dble_one,alloc_from_labile)

    ! When GPP is higher than seed C content, remaining seed carbon enters litter
    ! C pool, as seedlings do not fully exhaust their seed (P. de Vries p 48)
    if ( ( gpp_acm .gt. alloc_from_labile ) .and. ( use_seed_labile ) ) then
        stock_litter = stock_litter + stock_labile
        stock_labile = dble_zero
        use_seed_labile = .false.
    endif

    ! NPP as a fraction of GPP (1-.32=.68 or 68%) + allocation..
    npp = ( dble_one - frac_GPP_resp_auto ) * gpp_acm + alloc_from_labile
    ! from labile pool; = SHORT-TERM POOL    

    root_frac_intpol  = max(dble_zero,min(dble_one,root_frac_intpol))
    alloc_to_roots    = root_frac_intpol * npp         !
    shoot_frac_intpol = dble_one - root_frac_intpol    !
    npp_shoot         = npp - alloc_to_roots           ! NPP remaining after root growth==SHOOT fraction
    alloc_to_foliage  = fol_frac_intpol  * npp_shoot   !
    alloc_to_stem     = stem_frac_intpol * npp_shoot   !
    alloc_to_storage_organ = max(dble_zero,npp_shoot - alloc_to_foliage - alloc_to_stem)
    if ( alloc_to_storage_organ > dble_zero ) then  ! allocation flux to storage organ limited by maximum growth rate
        gso_max  = ( stock_storage_organ + 0.5d0 ) * rel_gso_max / steps_in_day
        alloc_to_storage_organ = min( alloc_to_storage_organ , gso_max )
        if ( sown ) then
           alloc_to_labile = ( npp_shoot - alloc_to_foliage - alloc_to_stem - alloc_to_storage_organ ) &
                           * ( dble_one - resp_cost_labile_trans )
           resp_cost_foliage_to_labile =  ( npp_shoot - alloc_to_foliage - alloc_to_stem - alloc_to_storage_organ ) &
                                       * resp_cost_labile_trans
        else
          alloc_to_labile             = dble_zero
          resp_cost_foliage_to_labile = dble_zero
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
    RDRSH = min( RDRSHMAX , max( dble_zero , RDRSHMAX * ( lai - LAICR ) / LAICR ) )
    if ( DS < dble_one ) then
       RDRDV = dble_zero
    else
       ! RDRDV dependant on DR and DS, values range typically between 0.02 <
       ! RDRDV < 0.25
!print*,"!! What RDRDV to use? !!"
!!$      RDRDV = DR /( max( 0.1 , 2. - DS ) )
!!$      RDRDV = RDRDV / 24. ! to get hourly senescence rate
       RDRDV = turnover_rate_foliage * ( dble_one / ( ( max( 2d0 - DS , 0.1d0 ) ) * 8d0 ) ) ** 2
    ENDIF

    ! relative leaf death rate is the maximum value of the arguments RDRSH and
    ! RDRDV
    RDR = max( RDRSH , RDRDV )

    ! remobilization of foliar C and allocation to dead leaves pool (gC.m-2.t-1)
    litterfall_foliage = stock_foliage * min(dble_one,ts_length * RDR)
    litterfall_stem    = stock_stem    * min(dble_one,ts_length * DR * turnover_rate_stem * dble(stmob)) ! remobstem
    litterfall_roots   = stock_roots   * min(dble_one,ts_length * turnover_rate_roots)

    ! remobilized C to NPP (from both leaves and stems) (gC.m-2.t-1)
    remob   = ( litterfall_foliage * 0.5d0 + litterfall_stem ) * ( dble_one - resp_cost_labile_trans )
    ! respiratory cost of C transfer (conversion from starch to photosynthates) (gC.m-2.t-1)
    Raremob = ( litterfall_foliage * 0.5d0 + litterfall_stem ) * resp_cost_labile_trans

    ! for mass balance calculate the decompostion efficency
    decomp_efficency = decomposition_rate &
                     / (decomposition_rate+mineralisation_rate_litter)

    ! total litter decomposition
    decomposition = stock_litter * (decomposition_rate+mineralisation_rate_litter) * resp_rate * ts_length

    ! heterotrophic respiration component 1: mineralisation of litter C pool (gC.m-2.t-1)
    resp_h_litter = decomposition * (dble_one - decomp_efficency)
    ! heterotrophic respiration component 2:  mineralisation of organic matter C pool (gC.m-2.t-1)
    resp_h_soilOrgMatter = stock_soilOrgMatter * min(dble_one,mineralisation_rate_soilOrgMatter * resp_rate * ts_length)

    ! decomposition of litter to soil organic matter (gC.m-2.t-1)
    decomposition = decomposition - resp_h_litter

    ! Recalculate Carbon Pools...

    stock_foliage       = max(dble_zero, stock_foliage + alloc_to_foliage - litterfall_foliage)
    stock_stem          = max(dble_zero, stock_stem + alloc_to_stem - litterfall_stem)
    stock_storage_organ = max(dble_zero, stock_storage_organ + alloc_to_storage_organ)
    stock_roots         = max(dble_zero, stock_roots         + alloc_to_roots   - litterfall_roots)
    stock_litter        = max(dble_zero, stock_litter + litterfall_roots - resp_h_litter - decomposition)
    stock_soilOrgMatter = max(dble_zero, stock_soilOrgMatter + decomposition    - resp_h_soilOrgMatter)
    stock_dead_foliage  = max(dble_zero, stock_dead_foliage  + litterfall_foliage * 0.5d0) ! remainder of litfol is remobilisedi
    stock_labile        = max(dble_zero, stock_labile + alloc_to_labile  - alloc_from_labile - resp_cost_labile_to_foliage + remob)
    
    
    ! respiratory pool: new photosynthates are added (gC.m-2.t-1)
    stock_resp_auto = stock_resp_auto + frac_GPP_resp_auto * gpp_acm
    ! autotrophic respiration; Ra (typically ~7% of respiratory pool) (gC.m-2.t-1)
    resp_auto = stock_resp_auto * min(dble_one,turnover_rate_resp_auto * ts_length)
    ! respiratory pool reduced by Ra (amount of C respired by plant)
    stock_resp_auto = max(dble_zero, stock_resp_auto - resp_auto)
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

    use CARBON_MODEL_MOD, only: dble_zero, dble_one

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
    ! NOTE: (doptmin+dble_one) < dmaxmin added to allow for EDC search period when "not
    ! allowed" parameter sets will be tried anyway
    if ( avtemp > tmin .and. avtemp < tmax .and. (doptmin+dble_one) < dmaxmin ) then
        fT = temperature_impact( doptmin , dmaxmin , dttmin )
    else
        fT = dble_zero
    endif

    ! calculation of photoperiod coefficient
    fP = photoperiod_impact( PHCR , PHSC )

    if ( emerged .and. ( DS < 2d0 ) ) then   ! sum up daily DR values between emergence and maturity (DS=2)

       if ( DS < dble_one ) then  ! in the vegetative phase (before flowering):

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

    use CARBON_MODEL_MOD, only: dble_zero, dble_one, doy

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
         tmp = max( avtemp - tmin , dble_zero )*days_in_step
         PHU = PHU + tmp

         ! set the development stage and emergence..
         if ( PHU >= PHUem ) then
           emerged = .true.
           DS = dble_zero
         else
           emerged = .false.
           DS = -dble_one
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

    use CARBON_MODEL_MOD, only: dble_zero, dble_one

    implicit none

    ! shoot biomass..
    Cshoot = stock_foliage + stock_stem + stock_storage_organ + stock_labile

    ! determine harvest index..
    HI = stock_storage_organ / Cshoot

    ! the stuff we actually want from the harvest...
    yield = stock_storage_organ !+ stock_stem * ( dble_one - st_res )

    
    ! the biomass that is harvested in addition to the storage-organ..
    BM_EX  = stock_foliage * ( dble_one - lv_res )          &
              + stock_stem * ( dble_one - st_res )          &
               + stock_dead_foliage * ( dble_one - lv_res ) &
                + stock_labile

    ! what's left (will fall to the ground)..
    stock_litter  = stock_litter                     &
                    + stock_resp_auto                &
                     + stock_foliage * lv_res        &
                      + stock_stem * st_res          &
                       + stock_dead_foliage * lv_res

    ! empty the plant stocks..
    stock_storage_organ = dble_zero
    stock_foliage       = dble_zero
    stock_stem          = dble_zero
    stock_dead_foliage  = dble_zero
    stock_labile        = dble_zero
    stock_resp_auto     = dble_zero

    ! roots stay in ground and slowly decompose (until/unless the field is
    ! ploughed)

    ! reset logical variables..
    sown    = .false.
    emerged = .false.
    ploughed = .false.
    DS = -dble_one ; fV = dble_zero ; fT = dble_zero ; fP = dble_zero ; VD = dble_zero

  end subroutine harvest
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  double precision function photoperiod_impact( PH_crit , PH_sens )

    use CARBON_MODEL_MOD, only: dble_zero, dble_one, dayl_hours

    ! Function to determine the coefficient for !
    ! photoperiod impact on developmental rate. !
    ! From Streck et al., 2003                  !

    implicit none

    ! arguments..
    double precision,intent(in) :: PH_crit, & ! critical photoperiod below which no development occurs
                                   PH_sens    ! photoperiod sensitivity

    photoperiod_impact = max(dble_zero, dble_one - exp ( - PH_Sens * ( dayl_hours - PH_crit ) ))

  end function photoperiod_impact
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine plough

    use CARBON_MODEL_MOD, only: dble_zero, dble_one

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
    stock_dead_foliage  = dble_zero
    stock_foliage       = dble_zero
    stock_labile        = dble_zero
    stock_roots         = dble_zero
    stock_stem          = dble_zero
    stock_storage_organ = dble_zero

    ! Reset the development stage & phenological heat units..
    ploughed = .true. ; DS = -dble_one ; PHU = dble_zero
    max_raso = dble_zero ; raso = dble_zero ; max_raso_old = dble_zero ; raso_old = dble_zero
    mean_alloc_to_storage_organ_old = dble_zero ; mean_alloc_to_storage_organ = dble_zero

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

    use CARBON_MODEL_MOD, only: dble_zero, dble_one

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
    vernalization = max( dble_zero , min( dble_one , ( VD ** 5 ) / ( ( VDh ** 5 ) + (VD ** 5 ) ) ) )

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
  
  
end module CARBON_MODEL_CROP_MOD
