
module CARBON_MODEL_MOD

implicit none

!!!!!!!!!!!
! Authorship contributions
!
! This source code is used for calibration of the
! Aggregated Canopy Model for Gross Primary Productivity and Evapotranspiration (ACM-GPP-ETv1), 
! also known as ACM2, and described in Smallman & Williams (2019) https://doi.org/10.5194/gmd-12-2227-2019.
! This code was written by T. L. Smallman (University of Edinburgh, t.l.smallman@ed.ac.uk)
! Subsequent modifications by:
! T. L. Smallman (University of Edinburgh, t.l.smallman@ed.ac.uk)
! See function/subroutine specific comments for exceptions and contributors
!!!!!!!!!!!

! make all private
private

! explicit publics
public:: CARBON_MODEL     &
         ,opt_max_scaling  &
         ,wSWP_time        &
         ,cica_time        &
         ,soil_frac_clay   &
         ,soil_frac_sand   &
         ,nos_soil_layers  &
         ,sw_par_fraction  &
         ,emissivity       &
         ,gs_demand_supply_ratio  &
         ,gs_total_canopy  &
         ,gb_total_canopy  &
         ,canopy_par_MJday_time &
         ,dim_1, dim_2      &
         ,nos_trees        &
         ,nos_inputs       &
         ,leftDaughter     &
         ,rightDaughter    &
         ,nodestatus       &
         ,xbestsplit       &
         ,nodepred         &
         ,bestvar

! arrays for the emulator, just so we load them once and that is it cos they be
! massive
integer ::    dim_1, & ! dimension 1 of response surface
              dim_2, & ! dimension 2 of response surface
          nos_trees, & ! number of trees in randomForest
         nos_inputs    ! number of driver inputs

double precision, allocatable, dimension(:,:) ::     leftDaughter, & ! left daughter for forest
                                                    rightDaughter, & ! right daughter for forets
                                                       nodestatus, & ! nodestatus for forests
                                                       xbestsplit, & ! for forest
                                                         nodepred, & ! prediction value for each tree
                                                          bestvar    ! for randomForests

!!!!!!!!!
! Parameters
!!!!!!!!!

! useful technical parameters
logical:: do_iWUE = .true., & ! Use iWUE or WUE for stomatal optimisation
 do_energy_balance = .false.   ! Calculate steady-state energy balance for GPP~Transpiration
double precision, parameter:: dble_zero = 0d0    &
                              ,dble_one = 1d0     &
                              ,vsmall = tiny(0d0) &
                              ,vlarge = huge(0d0)

integer, parameter:: nos_root_layers = 3, nos_soil_layers = nos_root_layers+1
double precision, parameter:: pi = 3.1415927d0,  &
                             pi_1 = 0.3183099d0,  & ! pi**(-1d0)
                              pi2 = 9.869604d0,   & ! pi**2d0
                           two_pi = 6.283185d0,   & ! pi*2d0
                       deg_to_rad = 0.01745329d0, & ! pi/180d0
              sin_dayl_deg_to_rad = 0.3979486d0,  & ! sin( 23.45d0*deg_to_rad )
                          gravity = 9.8067d0,     & ! acceleration due to gravity, ms-1
                            boltz = 5.670400d-8,  & ! Boltzmann constant (W.m-2.K-4)
                       emissivity = 0.96d0,       &
                      emiss_boltz = 5.443584d-8,   & ! emissivity*boltz
                  sw_par_fraction = 0.5d0,        & ! fraction of short-wave radiation which is PAR
                           freeze = 273.15d0,     &
                       gs_H2O_CO2 = 1.646259d0,   & ! The ratio of H20:CO2 diffusion for gs (Jones appendix 2)
            gs_H2Ommol_CO2mol_day = 142.2368d0,   & ! The ratio of H20:CO2 diffusion for gs, including seconds per day correction
                     gs_H2O_CO2_1 = 0.6074378d0,  & ! gs_H2O_CO2 ** (-1d0)
                gs_H2Ommol_CO2mol = 0.001646259d0, & ! gs_H2O_CO2*1d-3
                       gb_H2O_CO2 = 1.37d0,       & ! The ratio of H20:CO2 diffusion for gb (Jones appendix 2)
          partial_molar_vol_water = 18.05d-6,     & ! partial molar volume of water, m3 mol-1 at 20C
                   mol_to_g_water = 18d0,         & ! molecular mass of water
                 mmol_to_kg_water = 1.8d-5,       & ! milli mole conversion to kg
                     mol_to_g_co2 = 12d0,         & ! molecular mass of CO2 (g)
                       umol_to_gC = 1.2d-5,       & ! conversion of umolC -> gC
                       gC_to_umol = 83333.33d0,   & ! conversion of gC -> umolC; umol_to_gC**(-1d0)
                     g_to_mol_co2 = 0.08333333d0, &
!snowscheme       density_of_water = 998.9d0,         & ! density of  ! water kg.m-3
                   gas_constant_d = 287.04d0,     & ! gas constant for dry air (J.K-1.mol-1)
                             Rcon = 8.3144d0,     & ! Universal gas constant (J.K-1.mol-1)
                        vonkarman = 0.41d0,       & ! von Karman's constant
                      vonkarman_1 = 2.439024d0,   & ! 1/von Karman's constant
                      vonkarman_2 = 0.1681d0,     & ! von Karman's constant^2
                            cpair = 1004.6d0        ! Specific heat capacity of air; used in energy balance J.kg-1.K-1

! photosynthesis/respiration parameters
double precision, parameter :: &
                    kc_saturation = 310d0,        & ! CO2 half saturation, at reference temperature (298.15 K)
                 kc_half_sat_conc = 23.956d0,     & ! CO2 half sat, sensitivity coefficient
               co2comp_saturation = 36.5d0,       & ! CO2 compensation point, at reference temperature (298.15 K)
            co2comp_half_sat_conc = 9.46d0          ! CO2 comp point, sensitivity coefficient

! hydraulic parameters
double precision, parameter :: &
                       tortuosity = 2.5d0,        & ! tortuosity
                           gplant = 5d0,          & ! plant hydraulic conductivity (mmol m-1 s-1 MPa-1)
                      root_resist = 25d0,         & ! Root resistivity (MPa s g mmol−1 H2O)
                      root_radius = 0.00029d0,    & ! root radius (m) Bonen et al 2014 = 0.00029
                                                    ! Williams et al 1996 = 0.0001
                    root_radius_1 = root_radius**(-1d0), &
              root_cross_sec_area = pi*root_radius**2d0, & ! root cross sectional area (m2)
                                                           ! = pi*root_radius*root_radius
                     root_density = 0.31d6,       & ! root density (g biomass m-3 root)
                                                    ! 0.5e6 Williams et al 1996
                                                    ! 0.31e6 Bonan et al 2014
          root_mass_length_coef_1 = (root_cross_sec_area*root_density)**(-1d0), &
               const_sfc_pressure = 101325d0,     & ! (Pa)  Atmospheric surface pressure
                             head = 0.009807d0,   & ! head of pressure (MPa/m)
                           head_1 = 101.968d0       ! inverse head of pressure (m/MPa)

! structural parameters
double precision, parameter :: &
                    canopy_height = 9d0,          & ! canopy height assumed to be 9 m
                     tower_height = canopy_height+2d0, & ! tower (observation) height assumed to be 2 m above canopy
                         min_wind = 0.2d0,        & ! minimum wind speed at canopy top
                     min_drythick = 0.01d0,       & ! minimum dry thickness depth (m)
                        min_layer = 0.03d0,       & ! minimum thickness of the third rooting layer (m)
                      soil_roughl = 0.05d0,       & ! soil roughness length (m)
                   top_soil_depth = 0.15d0,       & ! thickness of the top soil layer (m)
                   mid_soil_depth = 0.15d0,       & ! thickness of the second soil layer (m)
                         min_root = 5d0,          & ! minimum root biomass (gBiomass.m-2)
!                          min_lai = 0.1d0,        & ! minimum LAI assumed for aerodynamic conductance calculations (m2/m2)
                          min_lai = 1.5d0,        & ! minimum LAI assumed for aerodynamic conductance calculations (m2/m2)
                        max_depth = 2d0,          & ! max rooting depth (m)
                           root_k = 100d0,        & ! root biomass needed to reach 50 % of max_depth (g/m2)
                      min_storage = 0.1d0           ! minimum canopy water (surface) storage (mm)

! timing parameters
double precision, parameter :: &
                 seconds_per_hour = 3600d0,         & ! Number of seconds per hour
                  seconds_per_day = 86400d0,        & ! Number of seconds per day
                seconds_per_day_1 = 1.157407d-05      ! Inverse of seconds per day

! hydraulic model variables
integer:: water_retention_pass, soil_layer, sunrise, sunset
double precision, dimension(nos_soil_layers):: soil_frac_clay, soil_frac_sand  ! clay and sand fractions of soil
double precision, dimension(nos_root_layers):: uptake_fraction, & !
                                                     water_flux, & ! potential transpiration flux (mmol.m-2.s-1)
                                                         demand    ! maximum potential canopy hydraulic demand
double precision, dimension(nos_soil_layers+1):: SWP, & ! soil water potential (MPa)
                                    soil_conductivity, & ! soil conductivity
                                            waterloss, & ! water loss from specific soil layers (m)
                                            watergain, & ! water gained by specfic soil layers (m)
                                          waterchange, & ! net water change by specific soil layers (m)
                                       field_capacity, & ! soil field capacity (m3.m-3)
                                       soil_waterfrac, & ! soil water content (m3.m-3)
                                             porosity, & ! soil layer porosity, (fraction)
                                     porosity_initial, & !
                               field_capacity_initial, & !
                               soil_waterfrac_initial, & !
                                          SWP_initial, &
                                      layer_thickness, & ! thickness of soil layers (m)
                      cond1, cond2, cond3, potA, potB    ! Saxton equation values

double precision:: root_reach, root_biomass, soil_depth, &
                  drythick, & ! estimate of the thickness of the dry layer at soil surface (m)
                      wSWP, & ! weighted soil water potential (MPa) used in GSI calculate.
                              ! Removes/limits the fact that very low root
                              ! density in young plants
                              ! give values too large for GSI to handle.
                    runoff, & ! runoff (kg.m-2.day-1)
                 underflow, & ! drainage from the bottom of soil column (kg.m-2.day-1)
  new_depth, previous_depth, & ! depth of bottom of soil profile
               canopy_wind, & ! wind speed (m.s-1) at canopy top
                     ustar, & ! friction velocity (m.s-1)
                  ustar_Uh, &
            air_density_kg, & ! air density kg/m3
            ET_demand_coef, & ! air_density_kg*vpd_kPa*cpair
                    roughl, & ! roughness length (m)
              displacement, & ! zero plane displacement (m)
                max_supply, & ! maximum water supply (mmolH2O/m2/day)
                     meant, & ! mean air temperature (oC)
                   meant_K, & ! mean air temperature (K)
                 maxt_lag1, &
                     leafT, & ! canopy temperature (oC)
        canopy_swrad_MJday, & ! canopy_absorbed shortwave radiation (MJ.m-2.day-1)
          canopy_par_MJday, & ! canopy_absorbed PAR radiation (MJ.m-2.day-1)
          soil_swrad_MJday, & ! soil absorbed shortwave radiation (MJ.m-2.day-1)
          canopy_lwrad_Wm2, & ! canopy absorbed longwave radiation (W.m-2)
            soil_lwrad_Wm2, & ! soil absorbed longwave radiation (W.m-2)
             sky_lwrad_Wm2, & ! sky absorbed longwave radiation (W.m-2)
      stomatal_conductance, & ! maximum stomatal conductance (mmolH2O.m-2.day-1)
     potential_conductance, & ! potential stomatal conductance (mmolH2O.m-2.day-1)
       minimum_conductance, & ! minimum stomatal conductance (mmolH2O.m-2ground.s-1)
   aerodynamic_conductance, & ! bulk surface layer conductance (m.s-1)
          soil_conductance, & ! soil surface conductance (m.s-1)
         convert_ms1_mol_1, & ! Conversion ratio for m.s-1 -> mol.m-2.s-1
        convert_ms1_mmol_1, & ! Conversion ratio for m/s -> mmol/m2/s
       air_vapour_pressure, & ! Vapour pressure of the air (kPa)
                    lambda, & ! latent heat of vapourisa/tion (J.kg-1)
                     psych, & ! psychrometric constant (kPa K-1)
                     slope, & ! Rate of change of saturation vapour pressure with temperature (kPa.K-1)
    water_vapour_diffusion, & ! Water vapour diffusion coefficient in (m2/s)
         dynamic_viscosity, & ! dynamic viscosity (kg.m-2.s-1)
       kinematic_viscosity, & ! kinematic viscosity (m2.s-1)
              snow_storage, & ! snow storage (kgH2O/m2)
            canopy_storage, & ! water storage on canopy (kgH2O.m-2)
      intercepted_rainfall    ! intercepted rainfall rate equivalent (kg.m-2.s-1)

! Module level variables for ACM_GPP_ET parameters
double precision :: &
                        delta_gs, & ! day length corrected gs increment mmolH2O/m2/dayl
                             avN, & ! average foliar N (gN/m2)
                            iWUE, & ! Intrinsic water use efficiency (gC/m2leaf/day/mmolH2Ogs)
                       iWUE_step, & ! Intrinsic water use efficiency for that day (gC/m2leaf/dayl/mmolH2Ogs)
                             NUE, & ! Photosynthetic nitrogen use efficiency at optimum temperature (oC)
                                    ! ,unlimited by CO2, light and photoperiod
                                    ! (gC/gN/m2leaf/day)
metabolic_limited_photosynthesis, & ! temperature, leaf area and foliar N limiterd photosynthesis (gC/m2/day)
    light_limited_photosynthesis, & ! light limited photosynthesis (gC/m2/day)
                              ci, & ! Internal CO2 concentration (ppm)
                          gb_mol, & ! Canopy boundary layer conductance (molCO2/m2/day)
                        rb_mol_1, & ! Canopy boundary layer resistance (day/m2/molCO2)
                     pn_max_temp, & ! Maximum temperature for photosynthesis (oC)
                     pn_opt_temp, & ! Optimum temperature fpr photosynthesis (oC)
                     pn_kurtosis, & ! Kurtosis of photosynthesis temperature response
                              e0, & ! Quantum yield gC/MJ/m2/day PAR
                    co2_half_sat, & ! CO2 at which photosynthesis is 50 % of maximum (ppm)
                  co2_comp_point, & ! CO2 at which photosynthesis > 0 (ppm)
                          minlwp, & ! min leaf water potential (MPa)
           soil_swrad_absorption, & ! Fraction of SW rad absorbed by soil
         soil_iso_to_net_coef_SW, & ! Coefficient relating soil isothermal net radiation to net.
        soil_iso_to_net_coef_LAI, & ! Coefficient relating soil isothermal net radiation to net.
           soil_iso_to_net_const, & ! Constant relating soil isothermal net radiation to net
       canopy_iso_to_net_coef_SW, & ! Coefficient relating SW to canopy isothermal to net adjustment
      canopy_iso_to_net_coef_LAI, & ! Coefficient relating LAI to canopy isothermal to net adjustment
         canopy_iso_to_net_const, & ! Constant relating canopy isothermal net radiation to net
           max_lai_lwrad_release, & ! 1-Max fraction of LW emitted from canopy to be
          lai_half_lwrad_release, & ! LAI at which LW emitted from canopy to be released at 50 %
             max_par_transmitted, & !
             max_nir_transmitted, & !
               max_nir_reflected, & !
               max_par_reflected    !

! Module level variables for step specific met drivers
double precision:: mint, & ! minimum temperature (oC)
                    maxt, & ! maximum temperature (oC)
      airt_zero_fraction, & ! fraction of air temperature above freezing
                   swrad, & ! incoming short wave radiation (MJ/m2/day)
                     co2, & ! CO2 (ppm)
                     doy, & ! Day of year
                rainfall, & ! rainfall (kgH2O/m2/s)
                snowfall, &
               snow_melt, & ! snow melt (kgH2O/m2/s)
                wind_spd, & ! wind speed (m/s)
                 vpd_kPa, & ! Vapour pressure deficit (kPa)
                   lai_1, & ! inverse of LAI
                     lai    ! leaf area index (m2/m2)

! Module level varoables for step specific timing information
double precision:: cos_solar_zenith_angle, & !
                    seconds_per_step, & !
                       days_per_step, & !
                     days_per_step_1, & !
                 dayl_hours_fraction, &
                        dayl_seconds, & ! day length in seconds
                      dayl_seconds_1, &
                          dayl_hours    ! day length in hours

! Metrics on photosynthetic activity
double precision, dimension(:), allocatable:: gs_demand_supply_ratio, & ! actual:potential stomatal conductance
                                               gs_total_canopy, &        ! stomatal conductance (mmolH2O/m2ground/day)
                                               gb_total_canopy, &        ! boundary conductance (mmolH2O/m2ground/day)
                                         canopy_par_MJday_time           ! Absorbed PAR by canopy (MJ/m2ground/day)

double precision, dimension(:), allocatable ::    deltat_1, & ! inverse of decimal days
                                                meant_time, &
                                       co2_half_saturation, & ! (ppm)
                                    co2_compensation_point, & ! (ppm)
                                 Cwood_labile_release_coef, & ! time series of labile release to wood
                                 Croot_labile_release_coef, & ! time series of labile release to root
                                               soilwatermm, &
                                                 cica_time, & ! Internal vs ambient CO2 concentrations
                                                 wSWP_time

contains

!
!--------------------------------------------------------------------
!
  subroutine CARBON_MODEL(start, finish, met, pars, deltat, nodays, lat_in, lai_out, NEE_out, FLUXES, POOLS &
                         ,nopars, nomet, nopools, nofluxes, GPP_out)

    ! The Data Assimilation Linked Ecosystem Carbon (DALEC) model uses
    ! a simple daily time step emulator to estimate GPP and evapotranspiration.
    ! The subroutine calls the Aggregated Canopy Model (ACM) to simulate GPP and
    ! the Aggregated Canopy Model for Evapo-Transpiration (ACM_ET) to simulate potential evapotranspiration.
    ! The purpose of the code here to allow for the running of these subrountes with prescibed LAI and fine root
    ! contents for the purpose of recalibrating their parameters

    implicit none

    ! declare input variables
    integer, intent(in):: start    &
                          ,finish   &
                          ,nopars   & ! number of paremeters in vector
                          ,nomet    & ! number of meteorological fields
                          ,nofluxes & ! number of model fluxes
                          ,nopools  & ! number of model pools
                          ,nodays     ! number of days in simulation

    double precision, intent(in):: met(nomet, nodays) & ! met drivers
                         ,deltat(nodays)    & ! time step in decimal days
                         ,pars(nopars)      & ! number of parameters
                         ,lat_in              ! site latitude (degrees)

    double precision, dimension(nodays), intent(inout):: lai_out & ! leaf area index
                                               ,GPP_out & ! Gross primary productivity
                                               ,NEE_out   ! net ecosystem exchange of CO2

    double precision, dimension((nodays+1), nopools), intent(inout):: POOLS  ! vector of ecosystem pools

    double precision, dimension(nodays, nofluxes), intent(inout):: FLUXES  ! vector of ecosystem fluxes

    ! declare general local variables
    double precision, dimension(nodays):: location
    double precision ::   lat &
                        ,infi &
                         ,tmp &
               ,act_pot_ratio &
               ,transpiration &
             ,soilevaporation &
              ,wetcanopy_evap &
            ,snow_sublimation &
               ,ET_pot, ET_net &
                     ,deltaWP & ! deltaWP (MPa) minlwp-soilWP
                  ,isothermal &
            ,deltaR, deltaTemp &
                        ,Rtot   ! MPa.s-1.m-2.mmol-1

    integer:: p, f, nxp, n

    ! met drivers are:
    ! 1st run day
    ! 2nd min daily temp (oC)
    ! 3rd max daily temp (oC)
    ! 4th Radiation (MJ.m-2.day-1)
    ! 5th CO2 (ppm)
    ! 6th DOY
    ! 7th precipitation (kg.m-2.s-1)
    ! 8th deforestation fraction
    ! 9th burnt area fraction
    ! 10th 21 day average min temperature
    ! 11th 21 day average photoperiod
    ! 12th 21 day average VPD
    ! 13th Forest management practice to accompany any clearing
    ! 14th average foliar nitrogen
    ! 15th optional LAI driver
    ! 16th Latitude as driver
    ! 17th Wind speed (m.s-1)
    ! 18th Vapour pressure deficit (Pa)
    ! 19th Input Rtot
    ! 20th Top soil sand percentage
    ! 21st Bottom soil sand percentage
    ! 22nd Top soil clay percentage
    ! 23rd Bottom soil clay percentage

    ! POOLS are:
    ! currently emply

    ! FLUXES are:
    ! 1 = GPP (gC.m-2.day-1)
    ! 2 = Latent energy (kg.m-2.day-1)

    ! PARAMETERS
    ! load ACM-GPP-ET parameters
    NUE                        = pars(1)   ! Photosynthetic nitrogen use efficiency at optimum temperature (oC)
                                           ! ,unlimited by CO2, light and photoperiod
                                           ! (gC/gN/m2leaf/day)
    pn_max_temp                = pars(2)   ! Maximum temperature for photosynthesis (oC)
    pn_opt_temp                = pars(3)   ! Optimum temperature for photosynthesis (oC)
    pn_kurtosis                = pars(4)   ! Kurtosis of photosynthesis temperature response
    e0                         = pars(5)   ! Quantum yield gC/MJ/m2/day PAR
    minlwp                     = pars(6)   ! minimum leaf water potential (MPa)
    soil_iso_to_net_coef_LAI   = pars(7)   ! Coefficient relating soil isothermal net radiation to net.
    iWUE                       = pars(8)   ! Intrinsic water use efficiency (gC/m2leaf/day/mmolH2Ogs)
    soil_swrad_absorption      = pars(9)   ! Fraction of SW rad absorbed by soil
    max_lai_lwrad_release      = pars(10)  ! 1-Max fraction of LW emitted from canopy to be released
    lai_half_lwrad_release     = pars(11)  ! LAI at which LW emitted from canopy to be released at 50 %
    soil_iso_to_net_coef_SW    = pars(12)  ! Coefficient relating soil isothermal net radiation to net.
    soil_iso_to_net_const      = pars(13)  ! Constant relating soil isothermal net radiation to net
    max_par_transmitted        = pars(14)  ! Max fraction of canopy incident PAR transmitted to soil
    max_nir_transmitted        = pars(15)  ! Max fraction of canopy incident NIR transmitted to soil
    max_par_reflected          = pars(16)  ! Max fraction of canopy incident PAR reflected to sky
    max_nir_reflected          = pars(17)  ! Max fraction of canopy incident NIR reflected to sky
    canopy_iso_to_net_coef_SW  = pars(18)  ! Coefficient relating SW to the adjustment between isothermal and net LW
    canopy_iso_to_net_const    = pars(19)  ! Constant relating canopy isothermal net radiation to net
    canopy_iso_to_net_coef_LAI = pars(20)  ! Coefficient relating LAI to the adjustment between isothermal and net LW

    ! load some values
    deltaWP = minlwp  ! leafWP-soilWP (i.e. -2-0)
    Rtot = dble_one
    maxt_lag1 = met(3, 1)
    leafT = maxt

! profiling example
!real:: start, finish
!call cpu_time(start)
!call cpu_time(finish)
!print '("Time1 = ",f9.5, " seconds.")',finish-start
    if (start == 1) then
  !    ! comment out when doing ACM recalibration
  !    soil_frac_sand = 40.0
  !    soil_frac_clay = 15.0

      ! SHOULD TURN THIS INTO A SUBROUTINE CALL AS COMMON TO BOTH DEFAULT AND
      ! CROPS
      if (.not.allocated(deltat_1)) then
         allocate(deltat_1(nodays), wSWP_time(nodays), soilwatermm(nodays), meant_time(nodays) &
                 ,gs_demand_supply_ratio(nodays), gs_total_canopy(nodays), gb_total_canopy(nodays) &
                 ,canopy_par_MJday_time(nodays), cica_time(nodays))
         deltat_1 = deltat**(-dble_one)
         meant_time = (met(2, 1:nodays) + met(3, 1:nodays)) * 0.5d0
         ! zero variables not done elsewhere
         water_flux = dble_zero
         ! initialise some time invarient parameters
         call saxton_parameters(soil_frac_clay, soil_frac_sand)
         call initialise_soils(soil_frac_clay, soil_frac_sand)
         soil_waterfrac_initial = soil_waterfrac
         SWP_initial = SWP
         field_capacity_initial = field_capacity
         porosity_initial = porosity
      else
         water_flux = dble_zero
         soil_waterfrac = soil_waterfrac_initial
         SWP = SWP_initial
         field_capacity = field_capacity_initial
         porosity = porosity_initial
      endif

      ! initial values
      seconds_per_step = seconds_per_day*deltat(1)
      canopy_storage = dble_zero; intercepted_rainfall = dble_zero
      ! initialise root reach based on initial conditions
      root_biomass = max(min_root, met(15, 1)*80d0*2d0)
      root_reach = max_depth*root_biomass / (root_k+root_biomass)
      ! determine initial soil layer thickness
      layer_thickness(1) = top_soil_depth; layer_thickness(2) = mid_soil_depth
      layer_thickness(3) = max(min_layer, root_reach-sum(layer_thickness(1:2)))
      layer_thickness(4) = max_depth-sum(layer_thickness(1:3))
      layer_thickness(5) = top_soil_depth
      previous_depth = max(top_soil_depth, root_reach)
      soil_depth = dble_zero; previous_depth = dble_zero
      ! needed to initialise soils
      call calculate_Rtot(Rtot)

    endif  ! start == 1

    infi = 0d0

    !
    ! Begin looping through each time step
    !

    do n = start, finish

      !!!!!!!!!!
      ! Specific hacks for ACM-recalibration with independent days
      !!!!!!!!!!

      ! reset the following variables
      canopy_storage = dble_zero; snow_storage = dble_zero; intercepted_rainfall = dble_zero

      !!!!!!!!!!
      ! assign drivers and update some prognostic variables
      !!!!!!!!!!

      ! lag drivers
!      maxt_lag1 = maxt

      ! Incoming drivers
      mint = met(2, n)  ! minimum temperature (oC)
      maxt = met(3, n)  ! maximum temperature (oC)
      leafT = maxt     ! initial canopy temperature (oC)
      meant = meant_time(n)  ! mean air temperature (oC)
      meant_K = meant+freeze
      airt_zero_fraction = (maxt-dble_zero) / (maxt-mint)  ! fraction of temperture period above freezing
      swrad = max(dble_zero, met(4, n))  ! incoming short wave radiation (MJ/m2/day)
      co2 = met(5, n)   ! CO2 (ppm)
      doy = met(6, n)   ! Day of year
      rainfall = max(0d0, met(7, n))  ! rainfall (kgH2O/m2/s)

      avN = met(14, n)  ! average foliar nitrogen (gN/m2)
      lai = met(15, n)  ! leaf area index (m2/m2)
      lat = met(16, n)  ! latitude (degrees)
      wind_spd = met(17, n)  ! wind speed (m/s)
      vpd_kPa = max(dble_zero, met(18, n)*1d-3)  ! Vapour pressure deficit (Pa -> kPa)
      Rtot = -9999d0  ! met(19, n)  ! Rtot (MPa.s-1.m-2.mmol-1))

      ! calculate daylength in hours and seconds
      call calculate_daylength((doy-(deltat(n)*0.5d0)), lat)
      ! plus various time related variables needed thoughout
      dayl_seconds_1 = dayl_seconds ** (-dble_one)
      seconds_per_step = seconds_per_day*deltat(n)
      dayl_hours_fraction = dayl_hours*0.04166667d0  ! 1/24 = 0.04166667
      iWUE_step = iWUE*dayl_hours_fraction
      days_per_step = deltat(n)
      days_per_step_1 = deltat_1(n)

      ! snowing or not...?
      if (mint < 0d0 .and. maxt > 0d0) then
          ! if minimum temperature is below freezing point then we weight the
          ! rainfall into snow or rain based on proportion of temperature below
          ! freezing
          snowfall = rainfall * (1d0-airt_zero_fraction); rainfall = rainfall-snowfall
          ! Add rainfall to the snowpack and clear rainfall variable
          snow_storage = snow_storage + (snowfall*seconds_per_step)

          ! Also melt some of the snow based on airt_zero_fraction
          ! default assumption is that snow is melting at 10 % per day light hour
          snow_melt = min(snow_storage, airt_zero_fraction*snow_storage*dayl_hours*0.1d0*deltat(n))
          snow_storage = snow_storage-snow_melt
      elseif (maxt <= 0d0) then
          ! if whole day is below freezing then we should assume that all
          ! precipitation is snowfall
          snowfall = rainfall; rainfall = 0d0; snow_melt = 0d0
          ! Add rainfall to the snowpack and clear rainfall variable
          snow_storage = snow_storage + (snowfall*seconds_per_step)
      else if (mint > 0d0) then
          ! otherwise we assume snow is melting at 10 % per day light hour
          snow_melt = min(snow_storage, snow_storage*dayl_hours*0.1d0*deltat(n))
          snow_storage = snow_storage-snow_melt
          snowfall = 0d0
      end if

      ! In the ACM recal-every time step may be from a different time step
      ! and thus has its own soil beneath which the model needs to be aware
      ! of...
      soil_frac_sand(1:2) = met(20, n)  ! top soil sand percentage
      soil_frac_sand(3:nos_soil_layers) = met(21, n)  ! bot
      soil_frac_clay(1:2) = met(22, n)  ! top soil clay percentage
      soil_frac_clay(3:nos_soil_layers) = met(23, n)  ! bot
      call saxton_parameters(soil_frac_clay, soil_frac_sand)
      call initialise_soils(soil_frac_clay, soil_frac_sand)

      ! calculate the minimum soil & root hydraulic resistance based on total
      ! fine root mass ! *2*2 => *RS*C->Bio
      ! note 80 is parameter which should match that found in the script which generated the training dataset
      root_biomass = max(min_root, met(15, n)*80d0*2d0)
      root_reach = max_depth*root_biomass / (root_k+root_biomass)
      ! determine initial soil layer thickness
      layer_thickness(1) = top_soil_depth; layer_thickness(2) = mid_soil_depth
      layer_thickness(3) = max(min_layer, root_reach-sum(layer_thickness(1:2)))
      layer_thickness(4) = max_depth-sum(layer_thickness(1:3))

      !!!!!!!!!!
      ! Calculate surface exchange coefficients
      !!!!!!!!!!

      ! calculate some temperature dependent meteorologial properties
      call meteorological_constants(maxt, maxt+freeze, vpd_kPa)
      convert_ms1_mmol_1 = convert_ms1_mol_1*1d3
      ! calculate aerodynamic using consistent approach with SPA
      call calculate_aerodynamic_conductance
      gb_total_canopy(n) = aerodynamic_conductance*convert_ms1_mol_1*1d3

      !!!!!!!!!!
      ! Determine net shortwave and isothermal longwave energy balance
      !!!!!!!!!!

      call calculate_radiation_balance
      canopy_par_MJday_time(n) = canopy_par_MJday

      !!!!!!!!!!
      ! Estimate approximate wet canopy evaporation and impact on energy balance
      !!!!!!!!!!

      ! estimate drythick for the current step
      drythick = max(min_drythick, top_soil_depth*max(0d0, (1d0 - (soil_waterfrac(1) / field_capacity(1)))))
      ! Soil surface (kgH2O.m-2.day-1)
      call calculate_soil_evaporation(soilevaporation)
      ! If snow present assume that soilevaporation is sublimation of soil first
      if (snow_storage > 0d0) then
          snow_sublimation = soilevaporation
          if (snow_sublimation*deltat(n) > snow_storage) snow_sublimation = snow_storage*deltat_1(n)
          soilevaporation = soilevaporation-snow_sublimation
          snow_storage = snow_storage - (snow_sublimation*deltat(n))
      else
          snow_sublimation = 0d0
      end if

      ! If desired calculate the steady-state energy balance
      ! NOTE: the current layout neglects the impact of temperature change on
      ! meteorological componets and therefore on evaporation (i.e. LW radiation
      ! impact only)
      if (do_energy_balance .and. lai > 0d0) then
          isothermal = canopy_lwrad_Wm2 + (canopy_swrad_MJday*1d6*seconds_per_day_1)
          call update_net_radiation(isothermal, leafT, lai, dble_one &
                                   ,dble_zero, aerodynamic_conductance, vpd_kPa &
                                   ,deltaTemp, deltaR)
          ! update long wave and canopy temperature based on potential canopy
          ! surface flux
          canopy_lwrad_Wm2 = canopy_lwrad_Wm2+deltaR
          leafT = leafT+deltaTemp
          ! Canopy intercepted rainfall evaporation (kgH2O/m2/day)
          call calculate_wetcanopy_evaporation(wetcanopy_evap, act_pot_ratio, canopy_storage)
          ! Restore temperature and radiation values
          leafT = leafT-deltaTemp; canopy_lwrad_Wm2 = canopy_lwrad_Wm2-deltaR
      else if (lai > 0d0) then
          ! Canopy intercepted rainfall evaporation (kgH2O/m2/day)
          call calculate_wetcanopy_evaporation(wetcanopy_evap, act_pot_ratio, canopy_storage)
      else
          intercepted_rainfall = 0d0; canopy_storage = 0d0; wetcanopy_evap = 0d0
      endif  ! do energy balance

      ! calculate Rtot if a value has not been provided
      if (Rtot == -9999d0) call calculate_Rtot(Rtot)
      wSWP_time(n) = wSWP; deltaWP = min(dble_zero, minlwp-wSWP)

      ! calculate radiation absorption and estimate stomatal conductance
      call calculate_stomatal_conductance(abs(deltaWP), Rtot)
      ! Estimate stomatal conductance relative to its minimum/maximum, i.e. how
      ! close are we to maxing out supply (note 0.01 taken from min_gs)
      gs_demand_supply_ratio(n) = (stomatal_conductance-minimum_conductance) / (potential_conductance-minimum_conductance)
      ! Store the canopy level stomatal conductance (mmolH2O/m2/day)
      gs_total_canopy(n) = stomatal_conductance

      !!!!!!!!!!
      ! GPP (gC.m-2.day-1)
      !!!!!!!!!!

      if (stomatal_conductance > vsmall) then
          call acm_gpp_stage_1; FLUXES(n, 1) = acm_gpp_stage_2(stomatal_conductance)
          cica_time(n) = ci/co2
!          FLUXES(n, 1) = max(dble_zero, acm_gpp(stomatal_conductance))
      else
          FLUXES(n, 1) = dble_zero
          cica_time(n) = 0d0
      endif

      !!!!!!!!!!
      ! Evaptranspiration (kgH2O.m-2.day-1)
      !!!!!!!!!!

      ! Canopy transpiration (kgH2O/m2/day)
      call calculate_transpiration(transpiration)
      ! restrict transpiration to positive only
      transpiration = max(dble_zero, transpiration)

      ! load GPP to output
      GPP_out(n) = FLUXES(n, 1)
      FLUXES(n, 2) = transpiration !+ wetcanopy_evap+soilevap
      FLUXES(n, 3) = soilevaporation+snow_sublimation
      FLUXES(n, 4) = wetcanopy_evap

      ! load LAI to output
      LAI_out(n) = lai

      do nxp = 1, nopools
         if (POOLS(n+1, nxp) /= POOLS(n+1, nxp) .or. POOLS(n+1, nxp) < 0d0 .or. &
            abs(POOLS(n+1, nxp)) == abs(log(infi))) then
            print*,"step",n, "POOL",nxp
            print*,"met",met(:,n)
            print*,"POOLS",POOLS(n, :)
            print*,"FLUXES",FLUXES(n, :)
            print*,"POOLS+1",POOLS(n+1, :)
            print*,"wSWP",wSWP
            print*,"waterfrac",soil_waterfrac
            stop
         endif
      enddo

      do nxp = 1, nofluxes
         if (FLUXES(n, nxp) /= FLUXES(n, nxp) .or. abs(FLUXES(n, nxp)) == abs(log(infi))) then
            print*,"step",n, "FLUX",nxp
            print*,"met",met(:,n)
            print*,"POOLS",POOLS(n, :)
            print*,"FLUXES",FLUXES(n, :)
            print*,"POOLS+1",POOLS(n+1, :)
            print*,"wSWP",wSWP
            print*,"waterfrac",soil_waterfrac
            stop
         endif
      enddo

    end do  ! nodays loop

  end subroutine CARBON_MODEL
  !
  !--------------------------------------------------------------------
  !
!  double precision function acm_gpp(gs)
!
!    ! the Aggregated Canopy Model, is a Gross Primary Productivity (i.e.
!    ! Photosynthesis) emulator which operates at a daily time step. ACM can be
!    ! paramaterised to provide reasonable results for most ecosystems.
!
!    implicit none
!
!    ! declare input variables
!    double precision, intent(in):: gs
!
!    ! declare local variables
!    double precision:: pn, pd, pp, qq, ci, mult, pl &
!                       ,gc, gs_mol, gb_mol
!
!    ! Temperature adjustments for Michaelis-Menten coefficients
!    ! for CO2 (kc) and O2 (ko) and CO2 compensation point
!    ! See McMurtrie et al., (1992) Australian Journal of Botany, vol 40, 657-677
!    co2_half_sat   = arrhenious(kc_saturation, kc_half_sat_conc, leafT)
!    co2_comp_point = arrhenious(co2comp_saturation, co2comp_half_sat_conc, leafT)
!
!    !
!    ! Metabolic limited photosynthesis
!    !
!
!    ! maximum rate of temperature and nitrogen (canopy efficiency) limited
!    ! photosynthesis (gC.m-2.day-1)
!    pn = lai*avN*NUE*opt_max_scaling(pn_max_temp, pn_opt_temp, pn_kurtosis, leafT)
!
!    !
!    ! Diffusion limited photosynthesis
!    !
!
!    ! daily canopy conductance (mmolH2O.m-2.s-1-> molCO2.m-2.day-1)
!    ! The ratio of H20:CO2 diffusion is 1.646259 (Jones appendix 2).
!    ! i.e. gcH2O*1.646259 = gcCO2
!    gs_mol = gs*seconds_per_day*gs_H2Ommol_CO2mol
!
!    ! canopy level boundary layer conductance unit change
!    ! (m.s-1 -> mol.m-2.day-1) assuming sea surface pressure only.
!    ! Note the ratio of H20:CO2 diffusion through leaf level boundary layer is
!    ! 1.37 (Jones appendix 2).
!    gb_mol = aerodynamic_conductance*seconds_per_day*convert_ms1_mol_1*gb_H2O_CO2
!    ! Combining in series the stomatal and boundary layer conductances
!    gc = (gs_mol ** (-dble_one) + gb_mol ** (-dble_one)) ** (-dble_one)
!
!    ! pp and qq represent limitation by metabolic (temperature & N) and
!    ! diffusion (co2 supply) respectively
!    pp = (pn*gC_to_umol)/gc; qq = co2_comp_point-co2_half_sat
!    ! calculate internal CO2 concentration (ppm or umol/mol)
!    mult = co2+qq-pp
!    ci = 0.5d0*(mult+sqrt((mult*mult)-4d0*(co2*qq-pp*co2_comp_point)))
!    ! calculate CO2 limited rate of photosynthesis (gC.m-2.day-1)
!    pd = (gc * (co2-ci)) * umol_to_gC
!    ! scale to day light period as this is then consistent with the light
!    ! capture period (1/24 = 0.04166667)
!    pd = pd*dayl_hours*0.04166667d0
!
!    !
!    ! Light limited photosynthesis
!    !
!
!    ! calculate light limted rate of photosynthesis (gC.m-2.day-1)
!    pl = e0*canopy_par_MJday
!
!    !
!    ! CO2 and light co-limitation
!    !
!
!    ! calculate combined light and CO2 limited photosynthesis
!    acm_gpp = pl*pd/(pl+pd)
!
!    ! don't forget to return
!    return
!
!  end function acm_gpp
  !
  !------------------------------------------------------------------
  !
  subroutine acm_gpp_stage_1

    ! Estimate the light and temperature limited photosynthesis components.
    ! See acm_gpp_stage_2() for estimation of CO2 supply limitation and
    ! combination of light, temperature and CO2 co-limitation

    implicit none

    ! Declare local variables
    double precision:: a, b, c, Pl_max
    double precision, parameter:: theta = 0.7  ! curvature parameter of light response

    !
    ! Metabolic limited photosynthesis
    !

    ! maximum rate of temperature and nitrogen (canopy efficiency) limited
    ! photosynthesis (gC.m-2.day-1 -> umolC/m2/day)
    metabolic_limited_photosynthesis = gC_to_umol*lai*avN*NUE*opt_max_scaling(pn_max_temp, pn_opt_temp, pn_kurtosis, leafT)

    !
    ! Light limited photosynthesis
    !

    ! calculate light limted rate of photosynthesis (gC.m-2.day-1)
    !light_limited_photosynthesis = e0*canopy_par_MJday
    Pl_max = lai*exp((log(avN*NUE)*0.750)+1.668d0)
    a = theta; b = -(e0*canopy_par_MJday+Pl_max); c = e0*canopy_par_MJday*Pl_max
    light_limitied_photosynthesis = (-b-sqrt(b**2 - (4*a*c))) / (2*a)
!ep = 5.313  ! Quantium yield at ligth compenesation point (gC/MJ)
!pl_max = exp(log(avN*NUE)*0.750+1.668)  ! Linear approximation of light limited photosynthesis?

    !
    ! Stomatal conductance independent variables for diffusion limited
    ! photosynthesis
    !

    ! Canopy level boundary layer conductance unit change
    ! (m.s-1 -> mol.m-2.day-1) assuming sea surface pressure only.
    ! Note the ratio of H20:CO2 diffusion through leaf level boundary layer is
    ! 1.37 (Jones appendix 2).
    gb_mol = aerodynamic_conductance*seconds_per_day*convert_ms1_mol_1*gb_H2O_CO2
    rb_mol_1 = (gb_mol)**(-1d0)

    ! Temperature adjustments for Michaelis-Menten coefficients
    ! for CO2 (kc) and O2 (ko) and CO2 compensation point
    ! See McMurtrie et al., (1992) Australian Journal of Botany, vol 40, 657-677
    co2_half_sat   = arrhenious(kc_saturation, kc_half_sat_conc, leafT)
    co2_comp_point = arrhenious(co2comp_saturation, co2comp_half_sat_conc, leafT)

    ! don't forget to return
    return

  end subroutine acm_gpp_stage_1
  !
  !------------------------------------------------------------------
  !
  double precision function acm_gpp_stage_2(gs)

    ! Combine the temperature (pn) and light (pl) limited gross primary productivity
    ! estimates with CO2 supply limited via stomatal conductance (gs).
    ! See acm_gpp_stage_1() for additional details on pn and pl calculation.

    implicit none

    ! declare input variables
    double precision, intent(in):: gs

    ! declare local variables
    double precision:: pp, qq, mult, gc, rc, pl, pn, pd

    !
    ! Diffusion limited photosynthesis
    !

    ! Daily canopy conductance (mmolH2O.m-2.s-1-> molCO2.m-2.day-1)
    ! The ratio of H20:CO2 diffusion is 1.646259 (Jones appendix 2).
    ! i.e. gcH2O*1.646259 = gcCO2 then all multiplied by 86400 seconds
    !
    ! Combining in series the stomatal and boundary layer conductances
    rc = (gs*gs_H2Ommol_CO2mol_day) ** (-1d0) + rb_mol_1

    ! pp and qq represent limitation by metabolic (temperature & N) and
    ! diffusion (co2 supply) respectively
    pp = metabolic_limited_photosynthesis*rc; qq = co2_comp_point-co2_half_sat
    ! calculate internal CO2 concentration (ppm or umol/mol)
    mult = co2+qq-pp
    ci = 0.5d0*(mult+sqrt((mult*mult)-4d0*(co2*qq-pp*co2_comp_point)))

    ! calculate CO2 limited rate of photosynthesis (gC.m-2.day-1)
    ! Then scale to day light period as this is then consistent with the light
    ! capture period (1/24 = 0.04166667)
    pd = ((co2-ci)/rc) * umol_to_gC*dayl_hours_fraction

    !
    ! Estimate CO2 and light co-limitation
    !

    ! calculate combined light and CO2 limited photosynthesis
    acm_gpp_stage_2 = light_limited_photosynthesis*pd/(light_limited_photosynthesis+pd)

    ! sanity check
    if (acm_gpp_stage_2 /= acm_gpp_stage_2 .or. acm_gpp_stage_2 < 0d0) acm_gpp_stage_2 = 0d0

    ! don't forget to return
    return

  end function acm_gpp_stage_2
  !
  !----------------------------------------------------------------------
  !
  double precision function find_gs_iWUE(gs_in)

    ! Calculate CO2 limited photosynthesis as a function of metabolic limited
    ! photosynthesis (pn), atmospheric CO2 concentration and stomatal
    ! conductance (gs_in). Photosynthesis is calculated twice to allow for
    ! testing of senstivity to iWUE.

    ! arguments
    double precision, intent(in):: gs_in

    ! local variables
    double precision:: tmp, airt_save, lw_save, &
                        isothermal, deltaTemp, deltaR
    double precision:: gs_high, gs_store, &
                        gpp_high, gpp_low, &
                        evap_high, evap_low

    !!!!!!!!!!
    ! Optimise intrinsic water use efficiency
    !!!!!!!!!!

    ! if desired calculate the steady-state energy balance
    if (do_energy_balance) then
        ! save values which will need to be reset
        airt_save = leafT; lw_save = canopy_lwrad_Wm2
        ! estimate energy balance without wet evaporation effects
        isothermal = canopy_lwrad_Wm2 + (canopy_swrad_MJday*1d6*dayl_seconds_1)
        call update_net_radiation(isothermal, leafT, lai, dble_one &
                                 ,gs_in, aerodynamic_conductance, vpd_kPa &
                                 ,deltaTemp, deltaR)
        ! note that both the leafT and canopy LW have an implicit day -> day length correction
        canopy_lwrad_Wm2 = canopy_lwrad_Wm2+deltaR
        leafT = leafT+deltaTemp
    endif
    ! estimate photosynthesis with current estimate of gs
    gpp_low = acm_gpp_stage_2(gs_in)

    ! Increment gs
    gs_high = gs_in+delta_gs
    ! if desired calculate the steady-state energy balance
    if (do_energy_balance) then
        leafT = airt_save; canopy_lwrad_Wm2 = lw_save
        ! estimate energy balance without wet evaporation effects
        isothermal = canopy_lwrad_Wm2 + (canopy_swrad_MJday*1d6*dayl_seconds_1)
        call update_net_radiation(isothermal, leafT, lai, dble_one &
                                 ,gs_in, aerodynamic_conductance, vpd_kPa &
                                 ,deltaTemp, deltaR)
        ! note that both the leafT and canopy LW have an implicit day -> day length correction
        canopy_lwrad_Wm2 = canopy_lwrad_Wm2+deltaR
        leafT = leafT+deltaTemp
    endif
    ! estimate photosynthesis with incremented gs
    gpp_high = acm_gpp_stage_2(gs_high)

    ! determine impact of gs increment on pd and how far we are from iWUE
    find_gs_iWUE = iWUE_step - ((gpp_high-gpp_low)*lai_1)

    ! now if I have been changing these drivers, best put them back to normal
    if (do_energy_balance) then
        leafT = airt_save; canopy_lwrad_Wm2 = lw_save
    endif

    ! remember to return back to the user
    return

  end function find_gs_iWUE
  !
  !----------------------------------------------------------------------
  !
  double precision function find_gs_WUE(gs_in)

    ! Calculate CO2 limited photosynthesis as a function of metabolic limited
    ! photosynthesis (pn), atmospheric CO2 concentration and stomatal
    ! conductance (gs_in). Photosynthesis is calculated twice to allow for
    ! testing of senstivity to WUE.

    ! arguments
    double precision, intent(in):: gs_in

    ! local variables
    double precision:: tmp, airt_save, lw_save, &
                        isothermal, deltaTemp, deltaR
    double precision:: gs_high, gs_store, &
                        gpp_high, gpp_low, &
                        evap_high, evap_low


    !!!!!!!!!!
    ! Optimise water use efficiency
    !!!!!!!!!!

    ! Globally stored upper stomatal conductance estimate in memory
    gs_store = stomatal_conductance
    ! now assign the current estimate
    stomatal_conductance = gs_in
    if (do_energy_balance) then
        ! save values which will need to be reset
        airt_save = leafT; lw_save = canopy_lwrad_Wm2
        ! estimate energy balance without wet evaporation effects
        isothermal = canopy_lwrad_Wm2 + (canopy_swrad_MJday*1d6*dayl_seconds_1)
        call update_net_radiation(isothermal, leafT, lai, dble_one &
                                 ,gs_in, aerodynamic_conductance, vpd_kPa &
                                 ,deltaTemp, deltaR)
        ! note that both the leafT and canopy LW have an implicit day -> day length correction
        canopy_lwrad_Wm2 = canopy_lwrad_Wm2+deltaR
        leafT = leafT+deltaTemp
    endif
    ! estimate photosynthesis with current estimate of gs
    gpp_low = acm_gpp_stage_2(gs_in)
    call calculate_transpiration(evap_low)

    ! Increment gs
    gs_high = gs_in+delta_gs
    ! now assign the incremented estimate
    stomatal_conductance = gs_high
    ! if desired calculate the steady-state energy balance
    if (do_energy_balance) then
        leafT = airt_save; canopy_lwrad_Wm2 = lw_save
        ! estimate energy balance without wet evaporation effects
        isothermal = canopy_lwrad_Wm2 + (canopy_swrad_MJday*1d6*dayl_seconds_1)
        call update_net_radiation(isothermal, leafT, lai, dble_one &
                                 ,gs_in, aerodynamic_conductance, vpd_kPa &
                                 ,deltaTemp, deltaR)
        ! note that both the leafT and canopy LW have an implicit day -> day length correction
        canopy_lwrad_Wm2 = canopy_lwrad_Wm2+deltaR
        leafT = leafT+deltaTemp
    endif
    ! estimate photosynthesis with incremented gs
    gpp_high = acm_gpp_stage_2(gs_high)
    call calculate_transpiration(evap_high)

    ! estimate marginal return on GPP for water loss, less water use efficiency criterion (gC.kgH2O-1.m-2.s-1)
    find_gs_WUE = ((gpp_high-gpp_low)/(evap_high-evap_low)) *lai_1
    find_gs_WUE = find_gs_WUE-iWUE_step

    ! return original stomatal value back into memory
    stomatal_conductance = gs_store

    ! now if I have been changing these drivers, best put them back to normal
    if (do_energy_balance) then
        leafT = airt_save; canopy_lwrad_Wm2 = lw_save
    endif

    ! remember to return back to the user
    return

  end function find_gs_WUE
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_stomatal_conductance(deltaWP, Rtot)

    use brent_zero, only: zbrent

    ! Determines 1) an approximation of canopy conductance (gc) mmolH2O.m-2.s-1
    ! based on potential hydraulic flow, air temperature and absorbed radiation.
    ! 2) calculates absorbed shortwave radiation (W.m-2) as function of LAI

    implicit none

    ! arguments
    double precision, intent(in):: deltaWP, & ! minlwp-wSWP (MPa)
                                       Rtot    ! total hydraulic resistance (MPa.s-1.m-2.mmol-1)

    ! local variables
    double precision:: denom, isothermal, deltaTemp, deltaR, pn, pl, pn_day, iWUE_lower, iWUE_upper
    double precision, parameter:: max_gs = 2000d0, &  ! mmolH2O.m-2.s-1 (leaf area)
                                   min_gs = 1d0, &  !
                                   tol_gs = 10d0        ! 4d0

    !!!!!!!!!!
    ! Calculate stomatal conductance under H2O and CO2 limitations
    !!!!!!!!!!

    if (aerodynamic_conductance > vsmall .and. deltaWP > vsmall) then

        ! Determine potential water flow rate (mmolH2O.m-2.dayl-1)
        max_supply = (deltaWP/Rtot) * seconds_per_day
        ! Pass minimum conductance from local parameter to global value
        ! There is uncertainty whether this should be a leaf area scaled value...
        minimum_conductance = min_gs*lai

        ! Invert Penman-Monteith equation to give gs (m.s-1) needed to meet
        ! maximum possible evaporation for the day.
        ! This will then be reduced based on CO2 limits for diffusion based
        ! photosynthesis
        denom = slope * ((canopy_swrad_MJday*1d6*dayl_seconds_1) + canopy_lwrad_Wm2) &
              + (ET_demand_coef*aerodynamic_conductance)
        denom = (denom / (lambda*max_supply*mmol_to_kg_water*dayl_seconds_1)) - slope
        potential_conductance = aerodynamic_conductance / (denom/psych)

        ! convert m.s-1 to mmolH2O.m-2.s-1
        potential_conductance = potential_conductance*convert_ms1_mmol_1
        ! if conditions are dew forming then set conductance to maximum as we are not going to be limited by water demand
        if (potential_conductance <= 0d0 .or. potential_conductance > max_gs*lai) potential_conductance = max_gs*lai

        ! If there is a positive demand for water then we will solve for photosynthesis limits on gs through iterative solution
        delta_gs = 1d-3*lai  ! mmolH2O/m2leaf/day
        ! Estimate inverse of LAI to avoid division in optimisation
        lai_1 = lai**(-1d0)
        ! Calculate stage one acm, temperature and light limitation which
        ! are independent of stomatal conductance effects
        call acm_gpp_stage_1
        if (do_iWUE) then
            ! Intrinsic WUE optimisation
            ! Check that the water restricted water range brackets the root solution for the bisection
            iWUE_upper = find_gs_iWUE(potential_conductance) !; iWUE_lower = find_gs_iWUE(min_gs)
            if ( iWUE_upper*find_gs_iWUE(min_gs) > 0d0 ) then
                 ! Then both proposals indicate that photosynthesis
                 ! would be increased by greater opening of the stomata
                 ! and is therefore water limited!
                 stomatal_conductance = potential_conductance
                 ! Exception being if both are positive-therefore assume
                 ! lowest
                 if (iWUE_upper > 0d0) stomatal_conductance = minimum_conductance
            else
                 ! In all other cases iterate
                 stomatal_conductance = zbrent('calculate_gs:find_gs_iWUE', &
                                               find_gs_iWUE, minimum_conductance, potential_conductance, tol_gs*lai, iWUE_step*0.10d0)
            end if
            ! Empirical fit to outputs generated by bisection procedure.
            ! Assumes that water supply is not limiting, thus there is still the need to estimate supply limit and apply as bookend.
            ! Note also that the order of covariates reflects their importance in the prediction, 
            ! i.e. R > 0.9 just for first independent variable
!            pn = metabolic_limited_photosynthesis
!            pn_day = metabolic_limited_photosynthesis*dayl_hours_fraction
!            pl = light_limited_photosynthesis
!            stomatal_conductance =   50.92693d0 &
!                                 + ( 14.73576d0    * ((pn_day*pl) / (pn_day+pl)) ) &
!                                 + (  1.0555d0     * pn )                          &
!                                 + ((-8.140542d-4) * pn**2d0 )                     &
!                                 + ((-0.7185823d0) * pl )                          &
!                                 + ((-1.565065d0)  * co2_comp_point )              &
!                                 + (( 0.2258834d0) * co2_half_sat )                &
!                                 + ((-2.486837d-4) * co2_half_sat**2d0 )           &
!                                 + (( 4.344512d-2) * co2 )                         &
!                                 + ((-2.969554d-4) * co2**2d0 )                    &
!                                 + ((-41.61914d0)  * iWUE )
!            stomatal_conductance = max(min_gs, min(stomatal_conductance, potential_conductance))
        else
            ! WUE optimisation
            stomatal_conductance = zbrent('acm_albedo_gc:find_gs_WUE',find_gs_WUE, min_gs, potential_conductance, tol_gs, iWUE*0.10d0)
        endif

        ! if desired calculate the steady-state energy balance
        if (do_energy_balance) then
            isothermal = canopy_lwrad_Wm2 + (canopy_swrad_MJday*1d6*dayl_seconds_1)
            call update_net_radiation(isothermal, leafT, lai, dble_one &
                                     ,stomatal_conductance, aerodynamic_conductance, vpd_kPa &
                                     ,deltaTemp, deltaR)
            ! note that both the leafT and canopy LW have an implicit day -> day length correction
            canopy_lwrad_Wm2 = canopy_lwrad_Wm2+deltaR
            leafT = leafT+deltaTemp
        endif

    else

        ! if no LAI then there can be no stomatal conductance
        potential_conductance = max_gs; stomatal_conductance = min_gs
        ! set minimum (computer) precision level flow
        max_supply = vsmall

    endif  ! if LAI > vsmall

  end subroutine calculate_stomatal_conductance
  !
  !------------------------------------------------------------------
  !
  subroutine meteorological_constants(input_temperature, input_temperature_K, input_vpd_kPa)

    ! Determine some multiple use constants used by a wide range of functions
    ! All variables here are linked to air temperature and thus invarient between
    ! iterations and can be stored in memory.
    ! NOTE: that some are strictly speaking linked to temperature of the surface
    ! being simulated, i.e. canopy or soil surface. However, those which have
    ! been identified as strongly sensitive to changes in energy balance are
    ! calculated in "meteorological_states"

    implicit none

    ! arguments
    double precision, intent(in):: input_temperature, input_temperature_K, &
                                    input_vpd_kPa

    ! local variables
    double precision:: s, mult

    !
    ! Used for soil, canopy evaporation and transpiration
    !

    ! Density of air (kg/m3)
    air_density_kg = 353d0/input_temperature_K
    ! Conversion ratio for m.s-1 -> mol.m-2.s-1
    convert_ms1_mol_1 = const_sfc_pressure / (input_temperature_K*Rcon)
    ! latent heat of vapourisation, 
    ! function of air temperature (J.kg-1)
    if (input_temperature < 0d0) then
        lambda = 2.835d6
    else
        lambda = 2501000d0-2364d0*input_temperature
    endif
    ! psychrometric constant (kPa K-1)
    psych = (0.0646d0*exp(0.00097d0*input_temperature))
    ! Straight line approximation of the true slope; used in determining
    ! relationship slope
    mult = input_temperature+237.3d0
    ! 2502.935945 = 0.61078*17.269*237.3
    ! Rate of change of saturation vapour pressure with temperature (kPa.K-1)
    slope = (2502.935945d0*exp(17.269d0*input_temperature/mult)) / (mult*mult)

    ! estimate frequently used atmmospheric demand component
    ET_demand_coef = air_density_kg*cpair*input_vpd_kPa

    !
    ! Used for soil evaporation and leaf level conductance
    !

    ! Determine diffusion coefficient (m2.s-1), temperature dependant (pressure dependence neglected). Jones p51; appendix 2
    ! Temperature adjusted from standard 20oC (293.15 K), NOTE that 1/293.15 = 0.003411223
    ! 0.0000242 = conversion to make diffusion specific for water vapor (um2.s-1)
    water_vapour_diffusion = 0.0000242d0*((input_temperature_K/293.2d0)**1.75d0)

    !
    ! Used for calculation of leaf level conductance
    !

    ! Calculate the dynamic viscosity of air (kg.m-2.s-1)
    dynamic_viscosity = ((input_temperature_K**1.5d0)/(input_temperature_K+120d0))*1.4963d-6
    ! and kinematic viscosity (m2.s-1)
    kinematic_viscosity = dynamic_viscosity/air_density_kg

  end subroutine meteorological_constants
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_transpiration(transpiration)

    ! Models leaf cnaopy transpiration based on the Penman-Monteith model of
    ! evapotranspiration used to estimate SPA's daily evapotranspiration flux
    ! (kgH20.m-2.day-1).

    implicit none

    ! arguments
    double precision, intent(out):: transpiration  ! kgH2O.m-2.day-1

    ! local variables
    double precision:: canopy_radiation & ! isothermal net radiation (W/m2)
                                  ,gs, gb   ! stomatal and boundary layer conductance (m.s-1)

    !!!!!!!!!!
    ! Estimate energy radiation balance (W.m-2)
    !!!!!!!!!!

    ! Absorbed shortwave radiation MJ.m-2.day-1 -> J.m-2.s-1
    canopy_radiation = canopy_lwrad_Wm2 + (canopy_swrad_MJday*1d6*dayl_seconds_1)

    !!!!!!!!!!
    ! Calculate canopy conductance (to water vapour)
    !!!!!!!!!!

    ! Change units of potential stomatal conductance
    ! (mmolH2O.m-2.s-1 -> m.s-1).
    ! Note assumption of sea surface pressure only
    gs = stomatal_conductance/convert_ms1_mmol_1
    ! Combine in series stomatal conductance with boundary layer
    gb = aerodynamic_conductance

    !!!!!!!!!!
    ! Calculate canopy evaporative fluxes (kgH2O/m2/day)
    !!!!!!!!!!

    ! Calculate numerator of Penman Montheith (kgH2O.m-2.day-1)
    ! NOTE: that restriction within water supply restriction is determined
    ! during stomatal conductance level.
    transpiration = ( ( (slope*canopy_radiation) + (ET_demand_coef*gb) ) &
                      / (lambda*(slope+(psych*(1d0+gb/gs)))) )*dayl_seconds

  end subroutine calculate_transpiration
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_wetcanopy_evaporation(wetcanopy_evap, act_pot_ratio, storage)

    ! Estimates evaporation of canopy intercepted rainfall based on the Penman-Monteith model of
    ! evapotranspiration used to estimate SPA's daily evapotranspiration flux
    ! (kgH20.m-2.day-1).

    implicit none

    ! arguments
    double precision, intent(inout):: storage         ! canopy water storage kgH2O/m2
    double precision, intent(out):: wetcanopy_evap, & ! kgH2O.m-2.day-1
                                      act_pot_ratio    ! Ratio of potential evaporation to actual

    ! local variables
    double precision:: canopy_radiation, & ! isothermal net radiation (W/m2)
                                      gb    ! stomatal and boundary layer conductance (m.s-1)

    !!!!!!!!!!
    ! Calculate canopy conductance (to water vapour)
    !!!!!!!!!!

    ! Combine in series stomatal conductance with boundary layer
    gb = aerodynamic_conductance

    !!!!!!!!!!
    ! Estimate energy radiation balance (W.m-2)
    !!!!!!!!!!

    ! Absorbed shortwave radiation MJ.m-2.day-1 -> J.m-2.s-1
    canopy_radiation = canopy_lwrad_Wm2 + (canopy_swrad_MJday*1d6*seconds_per_day_1)

    !!!!!!!!!!
    ! Calculate canopy evaporative fluxes (kgH2O/m2/day)
    !!!!!!!!!!

    ! Calculate numerator of Penman Montheith (kgH2O.m-2.day-1)
    wetcanopy_evap = max(0d0, (((slope*canopy_radiation) + (ET_demand_coef*gb)) / (lambda*(slope+psych))) * seconds_per_day)

    ! Remember potential evaporation to later calculation of the potential
    ! actual ratio
    act_pot_ratio = wetcanopy_evap

    ! assuming there is any rainfall, currently water on the canopy or dew formation
    if (rainfall > 0d0 .or. storage > 0d0) then
        ! Update based on canopy water storage
        call canopy_interception_and_storage(wetcanopy_evap, storage)
    else
        ! there is no water movement possible
        intercepted_rainfall = 0d0; wetcanopy_evap = 0d0
    endif

    ! now calculate the ratio of potential to actual evaporation
    if (act_pot_ratio == 0d0) then
        act_pot_ratio = 0d0
    else
        act_pot_ratio = abs(wetcanopy_evap/act_pot_ratio)
    endif

  end subroutine calculate_wetcanopy_evaporation
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_soil_evaporation(soilevap)

    ! Estimate soil surface evaporation based on the Penman-Monteith model of
    ! evapotranspiration used to estimate SPA's daily evapotranspiration flux
    ! (kgH20.m-2.day-1).

    implicit none

    ! arguments
    double precision, intent(out):: soilevap  ! kgH2O.m-2.day-1

    ! local variables
    double precision:: local_temp &
                        ,numerator &
                      ,denominator &
                   ,soil_radiation & ! isothermal net radiation (W/m2)
                            ,esurf & ! see code below
                             ,esat & ! soil air space saturation vapour pressure
                              ,gws   ! water vapour conductance through soil air space (m.s-1)

    ! oC -> K for local temperature value
    local_temp = maxt+freeze

    !!!!!!!!!!
    ! Estimate energy radiation balance (W.m-2)
    !!!!!!!!!!

    ! Absorbed shortwave radiation MJ.m-2.day-1 -> J.m-2.s-1
    soil_radiation = soil_lwrad_Wm2 + (soil_swrad_MJday*1d6*dayl_seconds_1)

    !!!!!!!!!!
    ! Calculate soil evaporative fluxes (kgH2O/m2/day)
    !!!!!!!!!!

    ! calculate saturated vapour pressure (kPa), function of temperature.
    esat = 0.1d0*exp( 1.80956664d0 + ( 17.2693882d0*local_temp-4717.306081d0 ) / ( local_temp-35.86d0 ) )
    air_vapour_pressure = esat-vpd_kPa

    ! Soil conductance to water vapour diffusion (m s-1)...
    gws = porosity(1) * water_vapour_diffusion / (tortuosity*drythick)

    ! vapour pressure in soil airspace (kPa), dependent on soil water potential
    ! - Jones p.110. partial_molar_vol_water. Less vapour pressure of the air to
    ! estimate the deficit between soil and canopy air spaces
    esurf = (esat*exp( 1d6*SWP(1) * partial_molar_vol_water / (Rcon*local_temp) )) - air_vapour_pressure

    ! Estimate the common denominator
!    denominator = (lambda*(slope+(psych*(1d0+soil_conductance/gws))))
!    denominator = (lambda*(slope+(psych*(soil_conductance/gws))))
    ! estimate the common atmospheric demand
!    numerator = (air_density_kg*cpair*esurf*soil_conductance)

    ! Estimate potential soil evaporation flux (kgH2O.m-2.day-1)
!    soilevap = ( ((slope*soil_radiation) + (air_density_kg*cpair*esurf*soil_conductance)) &
!               / (lambda*(slope+(psych*(1d0+soil_conductance/gws)))) ) * dayl_seconds
    soilevap = ( ((slope*soil_radiation) + (air_density_kg*cpair*esurf*soil_conductance)) &
               / (lambda*(slope+(psych*(soil_conductance/gws)))) ) * dayl_seconds
!    soilevap = ( (((slope*soil_radiation) + numerator) / denominator) * dayl_seconds ) &
!             + ( (((slope*soil_lwrad_Wm2) + numerator) / denominator) * (seconds_per_day-dayl_seconds))
    return

  end subroutine calculate_soil_evaporation
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_aerodynamic_conductance

    !
    ! Calculates the aerodynamic or bulk canopy conductance (m.s-1). Here we
    ! assume neutral conditions due to the lack of an energy balance calculation
    ! in either ACM or DALEC. The equations used here are with SPA at the time
    ! of the calibration
    !

    implicit none

    ! local variables
    double precision:: local_lai, &
           mixing_length_momentum, & ! mixing length parameter for momentum (m)
            length_scale_momentum    ! length scale parameter for momentum (m)

    ! Restrict LAI used here to greater than a minium value which prevents un-realistic outputs
    local_lai = max(min_lai, lai)

    ! calculate the zero plane displacement and roughness length
    call z0_displacement(ustar_Uh, local_lai)
    ! calculate friction velocity at tower height (reference height ) (m.s-1)
    ! WARNING neutral conditions only; see WRF module_sf_sfclay.F for 'with
    ! stability versions'
    !    ustar = (wind_spd/log((tower_height-displacement)/roughl)) * vonkarman
    ustar = wind_spd*ustar_Uh

    ! both length scale and mixing length are considered to be constant within
    ! the canopy (under dense canopy conditions) calculate length scale (lc)
    ! for momentum absorption within the canopy; Harman & Finnigan (2007)
    ! and mixing length (lm) for vertical momentum within the canopy Harman & Finnigan (2008)
    length_scale_momentum = (4d0*canopy_height) / local_lai
    mixing_length_momentum = 2d0*(ustar_Uh**3)*length_scale_momentum

    ! based on Harman & Finnigan (2008); neutral conditions only
    call log_law_decay

    ! now we are interested in the within canopy wind speed, 
    ! here we assume that the wind speed just inside of the canopy is most important.
    canopy_wind = canopy_wind*exp((ustar_Uh*((canopy_height*0.75d0)-canopy_height))/mixing_length_momentum)

    ! calculate_soil_conductance
    call calculate_soil_conductance(mixing_length_momentum, local_lai)

    ! calculate leaf level conductance (m/s) for water vapour under forced convective conditions
    call average_leaf_conductance(aerodynamic_conductance)

  end subroutine calculate_aerodynamic_conductance
  !
  !------------------------------------------------------------------
  !
  subroutine average_leaf_conductance(gv_forced)

    !
    ! Subroutine calculates the forced conductance of water vapour for non-cylinder within canopy leaves (i.e. broadleaf)
    ! Free convection (i.e. that driven by energy balance) is negelected here due to the lack of an energy balance
    ! calculation in DALEC. Should a energy balance be added then this code could be expanded include free conductance
    ! Follows a simplified approach to that used in SPA (Smallman et al 2013).
    !

    implicit none

    ! arguments
    double precision, intent(out):: gv_forced  ! canopy conductance (m/s) for water vapour under forced convection

    ! local parameters
    double precision, parameter:: leaf_width_coef = 25d0, & ! (1/leaf_width) * 0.5, 
                                                               ! where 0.5 accounts for one half
                                                               ! of the leaf used in water exchange
                                        leaf_width = 0.02d0    ! leaf width (m) (alternates 0.04, 0.08)
!                                                Pr = 0.72d0, & ! Prandtl number
!                                           Pr_coef = 1.05877d0  ! 1.18d0*(Pr**(0.33d0))
    ! local variables
    double precision :: &
              Sh_forced & ! Sherwood number under forced convection
                    ,Re   ! Reynolds number

    ! Sherwood number under forced convection. NOTE: 0.962*Pr_coef = 1.018537
!    Sh_forced = 0.962d0*Pr_coef*(sqrt((leaf_width*canopy_wind)/kinematic_viscosity))
    Sh_forced = 1.018537d0*(sqrt((leaf_width*canopy_wind)/kinematic_viscosity))
    ! Estimate the the forced conductance of water vapour
    gv_forced = water_vapour_diffusion*Sh_forced*leaf_width_coef*lai

  end subroutine average_leaf_conductance
  !
  !------------------------------------------------------------------
  !
  subroutine log_law_decay

    ! Standard log-law above canopy wind speed (m.s-1) decay under neutral
    ! conditions.
    ! See Harman & Finnigan 2008; Jones 1992 etc for details.

    implicit none

    ! log law decay
    canopy_wind = (ustar*vonkarman_1) * log((canopy_height-displacement) / roughl)

    ! set minimum value for wind speed at canopy top (m.s-1)
!    canopy_wind = max(min_wind, canopy_wind)

  end subroutine log_law_decay
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_field_capacity

    ! field capacity calculations for saxton eqns !
    
    use brent_zero, only: zbrent

    implicit none

    ! local variables..
    integer        :: i
    double precision:: x1, x2

    x1 = 0.1d0; x2 = 0.7d0  ! low/high guess
    do i = 1, nos_soil_layers+1
       water_retention_pass = i
       ! field capacity is water content at which SWP = -10 kPa
       field_capacity(i) = zbrent('water_retention:water_retention_saxton_eqns', &
                                   water_retention_saxton_eqns, x1, x2, 0.001d0, 0d0 )
    enddo

  end subroutine calculate_field_capacity
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_daylength(doy, lat)

    ! Subroutine uses day of year and latitude (-90/90 degrees) as inputs, 
    ! combined with trigonomic functions to calculate day length in hours and seconds

    implicit none

    ! arguments
    double precision, intent(in):: doy, lat

    ! local variables
    double precision:: dec, mult, sinld, cosld, aob

    !
    ! Estimate solar geometry variables needed
    !

    ! Declination
    ! NOTE: 0.002739726d0 = 1/365
    !    dec = - asin( sin( 23.45d0*deg_to_rad ) * cos( 2d0*pi * ( doy+10d0 ) / 365d0 ) )
    !    dec = - asin( sin_dayl_deg_to_rad*cos( two_pi * ( doy+10d0 ) / 365d0 ) )
    dec = - asin( sin_dayl_deg_to_rad*cos( two_pi * ( doy+10d0 ) * 0.002739726d0 ) )

    ! latitude in radians
    mult = lat*deg_to_rad
    ! day length is estimated as the ratio of sin and cos of the product of declination an latitude in radiation
    sinld = sin( mult ) * sin( dec )
    cosld = cos( mult ) * cos( dec )
    aob = max(-1d0, min(1d0, sinld/cosld))

    ! estimate day length in hours and seconds and upload to module variables
    dayl_hours = 12d0 * ( 1d0+2d0*asin( aob ) * pi_1 )
    dayl_seconds = dayl_hours*seconds_per_hour

    ! return to user
    return

  end subroutine calculate_daylength
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_longwave_isothermal(canopy_temperature, soil_temperature)

    ! Subroutine estimates the isothermal net longwave radiation (W.m-2) for
    ! the canopy and soil surface. SPA uses a complex multi-layer radiative
    ! transfer scheme including reflectance, transmittance any absorption.
    ! However, for a given canopy vertical profiles, the LAI absorption
    ! relationship is readily predicted via Michaelis-Menten or
    ! non-rectangular hyperbola as done here.

    implicit none

    ! arguments
    double precision, intent(in):: canopy_temperature, soil_temperature  ! oC

    ! local variables
    double precision:: lwrad, & ! downward long wave radiation from sky (W.m-2)
         transmitted_fraction, & ! fraction of LW which is not incident on the canopy
  canopy_transmitted_fraction, & !
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
    double precision, parameter:: nos_layers = 4d0 & ! Number of canopy layers in source model
                                  ,nos_layers_1 = nos_layers ** (-1d0) &
                                  ,clump = 1d0      & ! Clumping factor (1 = uniform, 0 totally clumped, mean = 0.75)
                                                      ! He et al., (2012) http://dx.doi.org/10.1016/j.rse.2011.12.008
                                  ,decay = -0.5d0     ! decay coefficient for incident radiation

    ! estimate long wave radiation from atmosphere (W.m-2)
    lwrad = emiss_boltz * (maxt+freeze-20d0) ** 4
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
    transmitted_fraction = exp(decay*lai*clump)

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
    canopy_transmitted_fraction = exp(decay*lai*0.5d0*clump)
    trans_lw_fraction     = 0.02d0*canopy_transmitted_fraction
    reflected_lw_fraction = 0.02d0*canopy_transmitted_fraction
    ! Absorption is the residual
    absorbed_lw_fraction = 1d0-trans_lw_fraction-reflected_lw_fraction

    ! Calculate the potential absorption of longwave radiation lost from the
    ! canopy to soil/sky.
    ! NOTE: That assuming the Beer's law emission from a single canopy layer leads to a rough 50 % underestimate of LW emission.
    !       This is why there is the nos_layer correction here
    canopy_release_fraction = (1d0 - (max_lai_lwrad_release*lai) / (lai+lai_half_lwrad_release)) &
                            * (1d0-exp(decay*lai*nos_layers_1*clump)) * nos_layers

    !!!!!!!!!!
    ! Distribute longwave from sky
    !!!!!!!!!!

    ! Estimate the radiation which directly bypasses the canopy...
    soil_incident_from_sky = lwrad*transmitted_fraction
    ! ...and update the canopy intercepted radiation
    lwrad = lwrad-soil_incident_from_sky

    ! long wave absorbed by the canopy from the sky
    canopy_absorption_from_sky = lwrad*absorbed_lw_fraction
    ! Long wave absorbed by soil from the sky, soil absorption assumed to be
    ! equal to emissivity
    soil_incident_from_sky = soil_incident_from_sky + (trans_lw_fraction*lwrad)
    soil_absorption_from_sky = soil_incident_from_sky*emissivity
    ! Long wave reflected directly back into sky
    sky_lwrad_Wm2 = lwrad*reflected_lw_fraction

    !!!!!!!!!!
    ! Distribute longwave from soil
    !!!!!!!!!!

    ! Calculate longwave radiation coming up from the soil plus the radiation
    ! which is reflected
    canopy_absorption_from_soil = longwave_release_soil + (soil_incident_from_sky * (1d0-emissivity))
    ! First how much directly bypasses the canopy...
    sky_lwrad_Wm2 = sky_lwrad_Wm2 + (canopy_absorption_from_soil*transmitted_fraction)
    canopy_absorption_from_soil = canopy_absorption_from_soil * (1d0-transmitted_fraction)
    ! Second, use this to estimate the longwave returning to the sky
    sky_lwrad_Wm2 = sky_lwrad_Wm2 + (canopy_absorption_from_soil*trans_lw_fraction)
    ! Third, now calculate the longwave from the soil surface absorbed by the
    ! canopy
    canopy_absorption_from_soil = canopy_absorption_from_soil*absorbed_lw_fraction

    !!!!!!!!!!
    ! Distribute longwave originating from the canopy itself
    !!!!!!!!!!

    ! calculate two-sided long wave radiation emitted from canopy which is
    ! ultimately lost from to soil or sky (i.e. this value is used twice, once
    ! to soil once to sky)
    canopy_loss = longwave_release_canopy*canopy_release_fraction
    ! Calculate longwave absorbed by soil which is released by the canopy itself
    soil_absorption_from_canopy = canopy_loss*emissivity
    ! Canopy released longwave returned to the sky
    sky_lwrad_Wm2 = sky_lwrad_Wm2+canopy_loss

    !!!!!!!!!!
    ! Isothermal net long wave canopy and soil balance (W.m-2)
    !!!!!!!!!!

    ! determine isothermal net canopy. Note two canopy_loss used to account for
    ! upwards and downwards emissions
    canopy_lwrad_Wm2 = (canopy_absorption_from_sky+canopy_absorption_from_soil) - (canopy_loss+canopy_loss)
    ! determine isothermal net soil
    soil_lwrad_Wm2 = (soil_absorption_from_sky+soil_absorption_from_canopy) - longwave_release_soil

  end subroutine calculate_longwave_isothermal
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_radiation_balance

    implicit none

    ! subroutine call ensures that both shortwave and longwave radiation balance
    ! are calculated at the same time but with the more readable code split
    ! between a shortwave and longwave specific subroutines.

    ! NOTE: that this code provides a daily timescale linear correction on
    ! isothermal longwave balance to net based on soil surface incident shortwave
    ! radiation

    ! declare local variables
    double precision:: delta_iso

    ! Estimate shortwave radiation balance
    call calculate_shortwave_balance
    ! Estimate isothermal long wave radiation balance
    call calculate_longwave_isothermal(meant, meant)
    ! Apply linear correction to soil surface isothermal->net longwave radiation
    ! balance based on absorbed shortwave radiation
    delta_iso = (soil_iso_to_net_coef_LAI*lai) + &
                (soil_iso_to_net_coef_SW * (soil_swrad_MJday*1d6*seconds_per_day_1)) + &
                soil_iso_to_net_const
    ! In addition to the iso to net adjustment, SPA analysis shows that soil net never gets much below zero
    soil_lwrad_Wm2 = max(-0.01d0, soil_lwrad_Wm2+delta_iso)
    ! Apply linear correction to canopy isothermal->net longwave radiation
    ! balance based on absorbed shortwave radiation
    delta_iso = (canopy_iso_to_net_coef_LAI*lai) + &
                (canopy_iso_to_net_coef_SW * (canopy_swrad_MJday*1d6*seconds_per_day_1)) + &
                canopy_iso_to_net_const
    canopy_lwrad_Wm2 = canopy_lwrad_Wm2+delta_iso

  end subroutine calculate_radiation_balance
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_shortwave_balance

    ! Subroutine estimates the canopy and soil absorbed shortwave radiation
    ! (MJ/m2/day).
    ! Radiation absorption is paritioned into NIR and PAR for canopy, and NIR +
    ! PAR for soil.

    ! SPA uses a complex multi-layer radiative transfer scheme including
    ! reflectance, transmittance any absorption. However, for a given
    ! canopy vertical profiles, the LAI absorption relationship is readily
    ! predicted via Michaelis-Menten or non-rectangular hyperbola as done here.

    implicit none

    ! local variables
    double precision:: balance                     &
                       ,transmitted_fraction        &
                       ,canopy_transmitted_fraction &
                       ,absorbed_nir_fraction_soil  &
                       ,absorbed_par_fraction_soil  &
                       ,fsnow, par, nir               &
                       ,soil_par_MJday              &
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
    double precision, parameter:: clump = 1d0    & ! Clumping factor (1 = uniform, 0 totally clumped, mean = 0.75)
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
    transmitted_fraction = exp(decay*lai*clump)

    ! Second, of the radiation which is incident on the canopy what fractions
    ! are transmitted through, reflected from or absorbed by the canopy

    canopy_transmitted_fraction = exp(decay*lai*0.5d0*clump)

    ! Canopy transmitted of PAR & NIR radiation towards the soil
    trans_par_fraction = canopy_transmitted_fraction*max_par_transmitted
    trans_nir_fraction = canopy_transmitted_fraction*max_nir_transmitted
    ! Canopy reflected of near infrared and photosynthetically active radiation
    reflected_nir_fraction = canopy_transmitted_fraction*max_nir_reflected
    reflected_par_fraction = canopy_transmitted_fraction*max_par_reflected
    ! Canopy absorption of near infrared and photosynthetically active radiation
    absorbed_nir_fraction = 1d0-reflected_nir_fraction-trans_nir_fraction
    absorbed_par_fraction = 1d0-reflected_par_fraction-trans_par_fraction

    !!!!!!!!!!
    ! Estimate canopy absorption of incoming shortwave radiation
    !!!!!!!!!!

    ! Estimate multiple use par and nir components
    par = sw_par_fraction*swrad
    nir = (1d0-sw_par_fraction) * swrad

    ! Estimate the radiation which directly bypasses the canopy...
    trans_par_MJday = par*transmitted_fraction
    trans_nir_MJday = nir*transmitted_fraction
    ! ...and update the canopy intercepted radiation
    par = par-trans_par_MJday
    nir = nir-trans_nir_MJday

    ! Estimate incoming shortwave radiation absorbed, transmitted and reflected
    ! by the canopy (MJ.m-2.day-1)
    canopy_par_MJday = par*absorbed_par_fraction
    canopy_nir_MJday = nir*absorbed_nir_fraction
    trans_par_MJday = trans_par_MJday + (par*trans_par_fraction)
    trans_nir_MJday = trans_nir_MJday + (nir*trans_nir_fraction)
    refl_par_MJday = par*reflected_par_fraction
    refl_nir_MJday = nir*reflected_nir_fraction

    !!!!!!!!!
    ! Estimate soil absorption of shortwave passing through the canopy
    !!!!!!!!!

    ! Update soil reflectance based on snow cover
    if (snow_storage > 0d0) then
        fsnow = 1d0-exp( - snow_storage*1d-2 )  ! fraction of snow cover on the ground
        absorbed_par_fraction_soil = ((1d0-fsnow) * soil_swrad_absorption) + (fsnow*newsnow_par_abs)
        absorbed_nir_fraction_soil = ((1d0-fsnow) * soil_swrad_absorption) + (fsnow*newsnow_nir_abs)
    else
        absorbed_par_fraction_soil = soil_swrad_absorption
        absorbed_nir_fraction_soil = soil_swrad_absorption
    endif

    ! Then the radiation incident and ultimately absorbed by the soil surface
    ! itself (MJ.m-2.day-1)
    soil_par_MJday = trans_par_MJday*absorbed_par_fraction_soil
    soil_nir_MJday = trans_nir_MJday*absorbed_nir_fraction_soil
    ! combine totals for use is soil evaporation
    soil_swrad_MJday = soil_nir_MJday+soil_par_MJday

    !!!!!!!!!
    ! Estimate canopy absorption of soil reflected shortwave radiation
    ! This additional reflection/absorption cycle is needed to ensure > 0.99
    ! of incoming radiation is explicitly accounted for in the energy balance.
    !!!!!!!!!

    ! calculate multiple use variables
    par = trans_par_MJday-soil_par_MJday
    nir = trans_nir_MJday-soil_nir_MJday
    ! how much of the reflected radiation directly bypasses the canopy...
    refl_par_MJday = refl_par_MJday + (par*transmitted_fraction)
    refl_nir_MJday = refl_nir_MJday + (nir*transmitted_fraction)
    ! ...and update the canopy on this basis
    par = par * (1d0-transmitted_fraction)
    nir = nir * (1d0-transmitted_fraction)

    ! Update the canopy radiation absorption based on the reflected radiation
    ! (MJ.m-2.day-1)
    canopy_par_MJday = canopy_par_MJday + (par*absorbed_par_fraction)
    canopy_nir_MJday = canopy_nir_MJday + (nir*absorbed_nir_fraction)
    ! Update the total radiation reflected back into the sky, i.e. that which is
    ! now transmitted through the canopy
    refl_par_MJday = refl_par_MJday + (par*trans_par_fraction)
    refl_nir_MJday = refl_nir_MJday + (nir*trans_nir_fraction)

    ! Combine to estimate total shortwave canopy absorbed radiation
    canopy_swrad_MJday = canopy_par_MJday+canopy_nir_MJday

    ! check energy balance
!    balance = swrad-canopy_par_MJday-canopy_nir_MJday-refl_par_MJday -
!    refl_nir_MJday-soil_swrad_MJday
!    if ((balance-swrad) / swrad > 0.01) then
!        print*,"SW residual frac = ",(balance-swrad) / swrad, "SW residual =
!        ",balance, "SW in = ",swrad
!    endif

  end subroutine calculate_shortwave_balance
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_Rtot(Rtot)

    ! Purpose of this subroutine is to calculate the minimum soil-root hydraulic
    ! resistance input into ACM. The approach used here is identical to that
    ! found in SPA.

    ! declare inputs
    double precision, intent(inout):: Rtot  ! MPa.s-1.m-2.mmol-1

    ! local variables
    integer:: i
    double precision:: bonus, sum_water_flux, &
                        transpiration_resistance, root_reach_local, &
                        root_depth_50
    double precision, dimension(nos_root_layers):: root_mass    &
                                                   ,root_length  &
                                                   ,ratio
    double precision, parameter:: root_depth_frac_50 = 0.25d0  ! fractional soil depth above which 50 %
                                                               ! of the root mass is assumed to be located

    ! reset water flux
    water_flux = 0d0; wSWP = 0d0
    ratio = 0d0; ratio(1) = 1d0; root_mass = 0d0
    ! calculate soil depth to which roots reach
    root_reach = max_depth*root_biomass / (root_k+root_biomass)
    ! calculate the plant hydraulic resistance component. Currently unclear
    ! whether this actually varies with height or whether tall trees have a
    ! xylem architecture which keeps the whole plant conductance (gplant) 1-10 (ish).
    !    transpiration_resistance = (gplant*lai)**(-1d0)
    transpiration_resistance = canopy_height / (gplant*max(min_lai, lai))

    !!!!!!!!!!!
    ! calculate current steps soil hydraulic conductivity
    !!!!!!!!!!!

    ! seperately calculate the soil conductivity as this applies to each layer
    do i = 1, nos_soil_layers
       call calculate_soil_conductivity(i, soil_waterfrac(i), soil_conductivity(i))
    end do  ! soil layers

    !!!!!!!!!!!
    ! Calculate root profile
    !!!!!!!!!!!

    ! The original SPA src generates an exponential distribution which aims
    ! to maintain 50 % of root biomass in the top 25 % of the rooting depth.
    ! In a simple 3 root layer system this can be estimates more simply

    ! top 25 % of root profile
    root_depth_50 = root_reach*root_depth_frac_50
    if (root_depth_50 <= layer_thickness(1)) then

        ! Greater than 50 % of the fine root biomass can be found in the top
        ! soil layer

        ! Start by assigning all 50 % of root biomass to the top soil layer
        root_mass(1) = root_biomass*0.5d0
        ! Then quantify how much additional root is found in the top soil layer
        ! assuming that the top 25 % depth is found somewhere within the top
        ! layer
        bonus = (root_biomass-root_mass(1)) &
              * (layer_thickness(1)-root_depth_50) / (root_reach-root_depth_50)
        root_mass(1) = root_mass(1) + bonus
        ! partition the remaining root biomass between the seconds and third
        ! soil layers
        if (root_reach > sum(layer_thickness(1:2))) then
            root_mass(2) = (root_biomass-root_mass(1)) &
                         * (layer_thickness(2)/(root_reach-layer_thickness(1)))
            root_mass(3) = root_biomass-sum(root_mass(1:2))
        else
            root_mass(2) = root_biomass-root_mass(1)
        endif

    else if (root_depth_50 > layer_thickness(1) .and. root_depth_50 <= sum(layer_thickness(1:2))) then

        ! Greater than 50 % of fine root biomass found in the top two soil
        ! layers. We will divide the root biomass uniformly based on volume, 
        ! plus bonus for the second layer (as done above)
        root_mass(1) = root_biomass * (layer_thickness(1)/root_depth_50)
        root_mass(2) = root_biomass * ((root_depth_50-layer_thickness(1))/root_depth_50)
        root_mass(1:2) = root_mass(1:2) * 0.5d0

        ! determine bonus for the seconds layer
        bonus = (root_biomass-sum(root_mass(1:2))) &
              * ((sum(layer_thickness(1:2))-root_depth_50)/(root_reach-root_depth_50))
        root_mass(2) = root_mass(2) + bonus
        root_mass(3) = root_biomass-sum(root_mass(1:2))

    else

        ! Greater than 50 % of fine root biomass stock spans across all three
        ! layers
        root_mass(1:2) = root_biomass*0.5d0 * (layer_thickness(1:2)/root_depth_50)
!        root_mass(1) = root_biomass * (layer_thickness(1)/root_depth_50)
!        root_mass(2) = root_biomass * (layer_thickness(2)/root_depth_50)
!        root_mass(1:2) = root_mass(1:2) * 0.5d0
        root_mass(3) = root_biomass-sum(root_mass(1:2))

    endif
    ! now convert root mass into lengths
    root_length = root_mass*root_mass_length_coef_1
!    root_length = root_mass / (root_density*root_cross_sec_area)

    !!!!!!!!!!!
    ! Calculate hydraulic properties and each rooted layer
    !!!!!!!!!!!

    ! calculate and accumulate steady state water flux in mmol.m-2.s-1
    ! NOTE: Depth correction already accounted for in soil resistance
    ! calculations and this is the maximum potential rate of transpiration
    ! assuming saturated soil and leaves at their minimum water potential.
    ! also note that the head correction is now added rather than
    ! subtracted in SPA equations because deltaWP is soilWP-minlwp not
    ! soilWP prior to application of minlwp
    demand = abs(minlwp-SWP(1:nos_root_layers))+head*canopy_height
    ! now loop through soil layers, where root is present
    do i = 1, nos_root_layers
       if (root_mass(i) > 0d0) then
           ! if there is root then there is a water flux potential...
           root_reach_local = min(root_reach, layer_thickness(i))
           ! calculate and accumulate steady state water flux in mmol.m-2.s-1
           water_flux(i) = plant_soil_flow(i, root_length(i), root_mass(i) &
                                          ,demand(i), root_reach_local, transpiration_resistance)
       else
           ! ...if there is not then we wont have any below...
           exit
       end if  ! root present in current layer?
    end do  ! nos_root_layers

    ! if freezing then assume soil surface is frozen
    if (meant < 1d0) then
        water_flux(1) = 0d0
        ratio(1) = 0d0
        ratio(2:nos_root_layers) = layer_thickness(2:nos_root_layers) / sum(layer_thickness(2:nos_root_layers))
    else
        ratio = layer_thickness(1:nos_root_layers)/sum(layer_thickness(1:nos_root_layers))
    endif

    ! calculate sum value
    sum_water_flux = sum(water_flux)
    if (sum_water_flux <= vsmall) then
        wSWP = -20d0; uptake_fraction = 0d0; uptake_fraction(1) = 1d0
    else
        ! calculate weighted SWP and uptake fraction
        wSWP = sum(SWP(1:nos_root_layers) * water_flux(1:nos_root_layers))
        uptake_fraction(1:nos_root_layers) = water_flux(1:nos_root_layers) / sum_water_flux
        wSWP = wSWP/sum_water_flux
    endif

    ! determine effective resistance (MPa.s-1.m-2.mmol-1)
    Rtot = sum(demand) / sum_water_flux

    ! finally convert transpiration flux (mmolH2O.m-2.s-1)
    ! into kgH2O.m-2.step-1 for consistency with ET in "calculate_update_soil_water"
    water_flux = water_flux*mmol_to_kg_water*seconds_per_step

    ! and return
    return

  end subroutine calculate_Rtot
  !
  !-----------------------------------------------------------------
  !
  subroutine canopy_interception_and_storage(potential_evaporation, storage)

    ! Simple daily time step integration of canopy rainfall interception, runoff
    ! and rainfall (kgH2O.m-2.s-1). NOTE: it is possible for intercepted rainfall to be
    ! negative if stored water running off into the soil is greater than
    ! rainfall (i.e. when leaves have died between steps)

    implicit none

    ! arguments
    double precision, intent(inout):: storage, & ! canopy water storage (kgH2O/m2)
                         potential_evaporation    ! wet canopy evaporation (kgH2O.m-2.day-1), 
                                                  ! enters as potential but leaves as water balance adjusted.
                                                  ! Note that this assumes a completely wet leaf surface
    ! local variables
    integer:: i, hr
    double precision:: a, through_fall, max_storage, max_storage_1, daily_addition, wetcanopy_evaporation &
                       ,potential_drainage_rate, drain_rate, evap_rate, initial_canopy, co_mass_balance, dx, dz, tmp(3)
    ! local parameters
    double precision, parameter:: CanIntFrac = -0.5d0,     & ! Coefficient scaling rainfall interception fraction with LAI
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
    through_fall = exp(CanIntFrac*lai*clump)
    ! maximum canopy storage (mm); minimum is applied to prevent errors in
    ! drainage calculation. Assume minimum capacity due to wood.
    max_storage = max(min_storage, CanStorFrac*lai)
    ! caclulate inverse for efficient calculations below
    max_storage_1 = max_storage**(-1d0)
    ! potential intercepted rainfall (kgH2O.m-2.s-1)
    intercepted_rainfall = rainfall * (1d0-through_fall)

    ! calculate drainage coefficients (Rutter et al 1975); Corsican Pine
    ! 0.002 is canopy specific coefficient modified by 0.002*(max_storage/1.05)
    ! where max_storage is the canopy maximum capacity (mm) (LAI based) and
    ! 1.05 is the original canopy capacitance
    a = log( RefDrainRate * ( max_storage*RefDrainLAI ) ) - RefDrainCoef*max_storage

    ! average rainfall intercepted by canopy (kgH2O.m-2.day-1)
    daily_addition = intercepted_rainfall*seconds_per_day

    ! reset cumulative variables
    through_fall = 0d0; wetcanopy_evaporation = 0d0
    drain_rate = 0d0; evap_rate = 0d0

    ! add rain to the canopy and overflow as needed
    storage = storage+daily_addition

    if (storage > max_storage) then

        if (potential_evaporation > 0d0) then

            ! assume co-access to available water above max_storage by both drainage and
            ! evaporation. Water below max_storage is accessable by evaporation only.

            ! Trapezium rule for approximating integral of drainage rate.
            ! Allows estimation of the mean drainage rate between starting
            ! point and max_storage, thus the time period appropriate for co-access can be
            ! quantified. NOTE 1440 = minutes/day
            ! General Formula: integral(rate) = 0.5*h((y0+yn) + 2(y1+y2 + ... yn-1)
            ! Where h id the size of the section, y0 is the maximum rate, yn is the final rate.
            dx = (storage-max_storage)*0.5d0
            tmp(1) = storage; tmp(2) = max_storage; tmp(3) = storage-dx
            tmp = exp(a + (RefDrainCoef*tmp))
            potential_drainage_rate = 0.5d0*dx * ((tmp(1) + tmp(2)) + 2d0*tmp(3)) * 1440d0
            ! To protect against un-realistic drainage rates
            ! due to very high rainfall rates
            potential_drainage_rate = min(potential_drainage_rate, vlarge)

            dz = storage-max_storage
            ! limit based on available water if total demand is greater than excess
            co_mass_balance = (dz / (potential_evaporation+potential_drainage_rate))
            evap_rate = potential_evaporation*co_mass_balance
            drain_rate = potential_drainage_rate*co_mass_balance

            ! Estimate evaporation from remaining water (i.e. that left after
            ! initial co-access of evaporation and drainage).
            ! Assume evaporation is now restricted by:
            ! 1) energy already spent on evaporation (the-evap_rate) and
            ! 2) linear increase in surface resistance as the leaf surface
            ! dries (i.e. the 0.5).
            evap_rate = evap_rate+min((potential_evaporation-evap_rate) * 0.5d0, storage-evap_rate-drain_rate)

        else

            ! Load dew formation to the current local evap_rate variable
            evap_rate = potential_evaporation
            ! Restrict drainage the quantity above max_storage, adding dew formation too
            drain_rate = (storage-evap_rate) - max_storage

        endif

    else

        ! no drainage just apply evaporation/dew formation fluxes directly
        drain_rate = 0d0; evap_rate = potential_evaporation
        if (evap_rate > 0d0) then
            ! evaporation restricted by fraction of surface actually covered
            ! in water and integrated over period to bare leaf (i.e. the*0.5)
            evap_rate = evap_rate*storage*max_storage_1*0.5d0
            ! and the total amount of water
            evap_rate = min(evap_rate, storage)
        else
            ! then dew formation has occurred, if this pushes storage > max_storage add it to drainage
            drain_rate = max(0d0, (storage-evap_rate) - max_storage)
      endif  ! evap_rate > 0

    endif  ! storage > max_storage

    ! update canopy storage with water flux
    storage = storage-evap_rate-drain_rate
    wetcanopy_evaporation = wetcanopy_evaporation+evap_rate
    through_fall = through_fall+drain_rate

    ! correct intercepted rainfall rate to kgH2O.m-2.s-1
    intercepted_rainfall = intercepted_rainfall - (through_fall*seconds_per_day_1)

!    ! sanity checks; note 1e-8 prevents precision errors causing flags
!    if (intercepted_rainfall > rainfall .or. storage < -1d-8 .or. &
!       (wetcanopy_evaporation*days_per_step_1) > (1d-8 + initial_canopy + (rainfall*seconds_per_day)) ) then
!        print*,"Condition 1",intercepted_rainfall > rainfall
!        print*,"Condition 2",storage < -1d-8
!        print*,"Condition 3",(wetcanopy_evaporation*days_per_step_1) > (1d-8 + initial_canopy + (rainfall*seconds_per_day))
!        print*,"storage (kgH2O/m2)",storage, "max_storage (kgH2O/m2)",max_storage, "initial storage (kgH2O/m2)", initial_canopy
!        print*,"rainfall (kgH2O/m2/day)", rainfall*seconds_per_day, "through_fall (kgH2O/m2/day)", (through_fall*days_per_step_1)
!        print*,"through_fall_total (kgH2O/m2/step)",through_fall
!        print*,"potential_evaporation (kgH2O/m2/day)",potential_evaporation
!        print*,"actual evaporation    (kgH2O/m2/day)",wetcanopy_evaporation*days_per_step_1
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
  subroutine calculate_update_soil_water(ET_leaf, ET_soil, rainfall_in, corrected_ET)

    !
    ! Function updates the soil water status and layer thickness
    ! Soil water profile is updated in turn with evaporative losses, 
    ! rainfall infiltration and gravitational drainage
    ! Root layer thickness is updated based on changes in the rooting depth from
    ! the previous step
    !

    implicit none

    ! arguments
    double precision, intent(in):: ET_leaf, ET_soil & ! evapotranspiration estimate (kgH2O.m-2.day-1)
                                       ,rainfall_in   ! rainfall (kgH2O.m-2.day-1)
    double precision, intent(out):: corrected_ET     ! water balance corrected evapotranspiration (kgH2O/m2/day)

    ! local variables
    integer:: day, a
    double precision:: depth_change, water_change, initial_soilwater, balance, mass_check
    double precision, dimension(nos_root_layers):: avail_flux, evaporation_losses, pot_evap_losses

    ! set soil water exchanges
    underflow = 0d0; runoff = 0d0; corrected_ET = 0d0; evaporation_losses = 0d0; pot_evap_losses = 0d0
    initial_soilwater = 1d3*sum(soil_waterfrac(1:nos_soil_layers) * layer_thickness(1:nos_soil_layers))

    ! Assume leaf transpiration is drawn from the soil based on the
    ! update_fraction estimated in calculate_Rtot
    pot_evap_losses = ET_leaf*uptake_fraction
    ! Assume all soil evaporation comes from the soil surface only
    pot_evap_losses(1) = pot_evap_losses(1) + ET_soil

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

    ! determine whether there is sufficient water to support evaporation
    water_change = minval((soil_waterfrac(1:nos_root_layers)*layer_thickness(1:nos_root_layers)) &
                         - (pot_evap_losses*days_per_step*1d-3))

    if (water_change > 0) then

       ! There is enough water to support evaporation across the whole time period...

       ! Draw all the water required for evaporation...
       ! adjust water already committed to evaporation
       ! convert kg.m-2 (or mm) -> Mg.m-2 (or m)
       soil_waterfrac(1:nos_root_layers) = soil_waterfrac(1:nos_root_layers) &
                                         + ((-pot_evap_losses*days_per_step*1d-3) / layer_thickness(1:nos_root_layers))
       ! Correct for dew formation; any water above porosity in the top layer is assumed runoff
       if (soil_waterfrac(1) > porosity(1)) then
           runoff = ((soil_waterfrac(1)-porosity(1)) * layer_thickness(1) * 1d3)
           soil_waterfrac(1) = porosity(1)
       endif

       ! determine infiltration from rainfall (kgH2O/m2/day), 
       ! if rainfall is probably liquid/soil surface is probably not frozen
       if (rainfall_in > 0d0) then
           ! reset soil water change variable
           waterchange = 0d0
           call infiltrate(rainfall_in*days_per_step)
           ! update soil profiles. Convert fraction into depth specific values
           ! (rather than m3/m3) then update fluxes
           soil_waterfrac(1:nos_soil_layers) = soil_waterfrac(1:nos_soil_layers) &
                                             + (waterchange(1:nos_soil_layers) / layer_thickness(1:nos_soil_layers))
           ! soil waterchange variable reset in gravitational_drainage()
       endif  ! is there any rain to infiltrate?

       ! determine drainage flux between surface -> sub surface
       call gravitational_drainage(nint(days_per_step))

       ! Pass information to the output ET variable
       corrected_ET = sum(pot_evap_losses)
       ! apply time step correction kgH2O/m2/step -> kgH2O/m2/day
       underflow = underflow*days_per_step_1
       runoff = runoff*days_per_step_1

    else

       ! to allow for smooth water balance integration carry this out at daily time step
       do day = 1, nint(days_per_step)

          !!!!!!!!!!
          ! Evaporative losses
          !!!!!!!!!!

          ! load potential evaporative losses from the soil profile
          evaporation_losses = pot_evap_losses
          ! can not evaporate from soil more than is available (m -> mm)
          ! NOTE: This is due to the fact that both soil evaporation and transpiration
          !       are drawing from the same water supply.
          avail_flux = soil_waterfrac(1:nos_root_layers) * layer_thickness(1:nos_root_layers) * 1d3
          do a = 1, nos_root_layers  ! note: timed comparison between "where" and do loop supports do loop for smaller vectors
             if (evaporation_losses(a) > avail_flux(a)) evaporation_losses(a) = avail_flux(a) * 0.999d0
          end do
          ! this will update the ET estimate outside of the function
          ! days_per_step corrections happens outside of the loop below
          corrected_ET = corrected_ET+sum(evaporation_losses)

          ! adjust water already committed to evaporation
          ! convert kg.m-2 (or mm) -> Mg.m-2 (or m)
          soil_waterfrac(1:nos_root_layers) = soil_waterfrac(1:nos_root_layers) &
                                            + ((-evaporation_losses(1:nos_root_layers)*1d-3) / layer_thickness(1:nos_root_layers))
          ! Correct for dew formation; any water above porosity in the top layer is assumed runoff
          if (soil_waterfrac(1) > porosity(1)) then
              runoff = runoff + ((soil_waterfrac(1)-porosity(1)) * layer_thickness(1) * 1d3)
              soil_waterfrac(1) = porosity(1)
          endif

          !!!!!!!!!!
          ! Rainfall infiltration drainage
          !!!!!!!!!!

          ! Determine infiltration from rainfall (kgH2O/m2/day), 
          ! if rainfall is probably liquid/soil surface is probably not frozen
          if (rainfall_in > 0d0) then
              ! reset soil water change variable
              waterchange = 0d0
              call infiltrate(rainfall_in)
              ! update soil profiles. Convert fraction into depth specific values
              ! (rather than m3/m3) then update fluxes
              soil_waterfrac(1:nos_soil_layers) = soil_waterfrac(1:nos_soil_layers) &
                                                + (waterchange(1:nos_soil_layers) / layer_thickness(1:nos_soil_layers))
              ! soil waterchange variable reset in gravitational_drainage()
          endif  ! is there any rain to infiltrate?

          !!!!!!!!!!
          ! Gravitational drainage
          !!!!!!!!!!

          ! Determine drainage flux between surface -> sub surface
          call gravitational_drainage(1)

       end do  ! days_per_step

       ! apply time step correction kgH2O/m2/step -> kgH2O/m2/day
       corrected_ET = corrected_ET*days_per_step_1
       underflow = underflow*days_per_step_1
       runoff = runoff*days_per_step_1

    end if  ! water_change > 0

    !!!!!!!!!!
    ! Update soil layer thickness
    !!!!!!!!!!

    depth_change = (top_soil_depth+mid_soil_depth+min_layer); water_change = 0
    ! if roots extent down into the bucket
    if (root_reach > depth_change .and. previous_depth <= depth_change) then

        !!!!!!!!!!
        ! Soil profile is within the bucket layer (layer 3)
        !!!!!!!!!!

        if (previous_depth > depth_change) then
            ! how much has root depth extended since last step?
            depth_change = root_reach-previous_depth
        else
            ! how much has root depth extended since last step?
            depth_change = root_reach-depth_change
        endif

        ! if there has been an increase
        if (depth_change > 0.05d0) then

            ! determine how much water (mm) is within the new volume of soil
            water_change = soil_waterfrac(nos_soil_layers) * depth_change
            ! now assign that new volume of water to the deep rooting layer
            soil_waterfrac(nos_root_layers) = ((soil_waterfrac(nos_root_layers)*layer_thickness(nos_root_layers))+water_change) &
                                            / (layer_thickness(nos_root_layers)+depth_change)

            ! explicitly update the soil profile if there has been rooting depth
            ! changes
            layer_thickness(1) = top_soil_depth; layer_thickness(2) = mid_soil_depth
            layer_thickness(3) = root_reach-sum(layer_thickness(1:2))
            layer_thickness(4) = max_depth-sum(layer_thickness(1:3))

            ! keep track of the previous rooting depth
            previous_depth = root_reach

        else if (depth_change < -0.05d0) then

            ! make positive to ensure easier calculations
            depth_change = -depth_change

            ! determine how much water is lost from the old volume of soil
            water_change = soil_waterfrac(nos_root_layers) * depth_change
            ! now assign that new volume of water to the deep rooting layer
            soil_waterfrac(nos_soil_layers) = ((soil_waterfrac(nos_soil_layers)*layer_thickness(nos_soil_layers))+water_change) &
                                            / (layer_thickness(nos_soil_layers)+depth_change)

            ! explicitly update the soil profile if there has been rooting depth
            ! changes
            layer_thickness(1) = top_soil_depth; layer_thickness(2) = mid_soil_depth
            layer_thickness(3) = root_reach-sum(layer_thickness(1:2))
            layer_thickness(4) = max_depth-sum(layer_thickness(1:3))

            ! keep track of the previous rooting depth
            previous_depth = root_reach

        else

            ! keep track of the previous rooting depth
            previous_depth = previous_depth

        end if  ! depth change

    else if (root_reach < depth_change .and. previous_depth > depth_change) then

        !!!!!!!!!!
        ! Model has explicitly contracted from the bucket layer
        !!!!!!!!!!

        ! In this circumstance we want to return the soil profile to it's
        ! default structure with a minimum sized third layer
        depth_change = previous_depth-depth_change

        ! determine how much water is lost from the old volume of soil
        water_change = soil_waterfrac(nos_root_layers) * depth_change
        ! now assign that new volume of water to the deep rooting layer
        soil_waterfrac(nos_soil_layers) = ((soil_waterfrac(nos_soil_layers)*layer_thickness(nos_soil_layers))+water_change) &
                                        / (layer_thickness(nos_soil_layers)+depth_change)

        ! explicitly update the soil profile if there has been rooting depth
        ! changes
        layer_thickness(1) = top_soil_depth; layer_thickness(2) = mid_soil_depth
        layer_thickness(3) = min_layer
        layer_thickness(4) = max_depth-sum(layer_thickness(1:3))

        ! keep track of the previous rooting depth
        previous_depth = min_layer

    else  ! root_reach > (top_soil_depth+mid_soil_depth+min_layer)

        ! if we are outside of the range when we need to consider rooting depth changes keep track in case we move into a zone when we do
        previous_depth = previous_depth

    endif  ! root reach beyond top layer

    ! finally update soil water potential
    call soil_water_potential

!    ! check water balance
!    balance = (rainfall_in-corrected_ET-underflow-runoff) * days_per_step
!    balance = balance &
!            - (sum(soil_waterfrac(1:nos_soil_layers) * layer_thickness(1:nos_soil_layers) * 1d3) &
!            - initial_soilwater)
!
!    if (abs(balance) > 1d-6 .or. soil_waterfrac(1) < -1d-6) then
!        print*,"Soil water miss-balance (mm)",balance
!        print*,"Initial_soilwater (mm) = ",initial_soilwater
!        print*,"Final_soilwater (mm) = ",sum(soil_waterfrac(1:nos_soil_layers) * layer_thickness(1:nos_soil_layers) * 1d3)
!        print*,"State balance = ",sum(soil_waterfrac(1:nos_soil_layers)*layer_thickness(1:nos_soil_layers)*1d3)-initial_soilwater
!        print*,"Flux balance = ",(rainfall_in-corrected_ET-underflow-runoff) * days_per_step
!        print*,"Top soilwater (fraction)",soil_waterfrac(1)
!        print*,"Rainfall (mm/step)",rainfall_in, "ET",corrected_ET, "underflow",underflow, "runoff",runoff
!        print*,"Rainfall (kgH2O/m2/s)",rainfall
!        print*,"Soil Water Fraction = ",soil_waterfrac
!    end if  ! abs(balance) > 1d-10

    ! explicit return needed to ensure that function runs all needed code
    return

  end subroutine calculate_update_soil_water
  !
  !-----------------------------------------------------------------
  !
  subroutine infiltrate(rainfall_in)

    ! Takes surface_watermm and distributes it among top !
    ! layers. Assumes total infilatration in timestep.   !
    ! NOTE: Assumes that any previous water movement due to infiltration and evaporation
    !       has already been updated in soil mass balance

    implicit none

    ! arguments
    double precision, intent(in):: rainfall_in  ! rainfall (kg.m-2.day-1)

    ! local argumemts
    integer:: i
    double precision:: add, & ! surface water available for infiltration (m)
                      wdiff    ! available space in a given soil layer for water to fill (m)

    ! convert rainfall water from mm -> m (or kgH2O.m-2.day-1 -> MgH2O.m-2.day-1)
    add = rainfall_in*1d-3

    do i = 1, nos_soil_layers

       ! is the input of water greater than available space
       ! if so fill and subtract from input and move on to the next
       ! layer determine the available pore space in current soil layer
       wdiff = max(0d0, (porosity(i)-soil_waterfrac(i)) * layer_thickness(i))

       if (add > wdiff) then
           ! if so fill and subtract from input and move on to the next layer
           waterchange(i) = waterchange(i) + wdiff
           add = add-wdiff
       else
           ! otherwise infiltate all in the current layer
           waterchange(i) = waterchange(i) + add
           add = 0d0; exit
       end if

    end do  ! nos_soil_layers

    ! if after all of this we have some water left assume it is runoff (kgH2O.m-2.day-1)
    ! NOTE that runoff is reset outside of the daily soil loop
    runoff = runoff + (add*1d3)

  end subroutine infiltrate
  !
  !-----------------------------------------------------------------
  !
  subroutine gravitational_drainage(time_period_days)

    ! Integrator for soil gravitational drainage.
    ! Due to the longer time steps undertake by ACM/DALEC and the fact that
    ! drainage is a concurrent processes we assume that drainage occurs at
    ! the bottom of the column first creating space into which water can drain
    ! from the top down. Therefore we draing from the bottom first and then the top.
    ! NOTE: Assumes that any previous water movement due to infiltration and evaporation
    !       has already been updated in soil mass balance

    implicit none

    ! arguments
    integer, intent(in):: time_period_days

    ! local variables..
    integer:: t
    double precision, dimension(nos_soil_layers):: dx, & ! range between the start and end points of the integration
                                               halfway, & ! half way point between start and end point of integration
                                                liquid, & ! liquid water in local soil layer (m3/m3)
                                         avail_to_flow, & ! liquid content above field capacity (m3/m3)
                                               iceprop, & ! fraction of soil layer which is ice
                                          pot_drainage    ! estimats of time step potential drainage rate (m/s)
    double precision  :: tmp1, tmp2, tmp3 &
                                 ,unsat & ! unsaturated pore space in soil_layer below the current (m3/m3)
                                ,change   ! absolute volume of water drainage in current layer (m3/day)

    ! calculate soil ice proportion; at the moment
    ! assume everything liquid
    iceprop = 0d0

    ! except the surface layer in the mean daily temperature is < 0oC
    if (meant < 1d0) iceprop(1) = 1d0

    ! zero water fluxes
    waterchange = 0d0

    ! underflow is tracked in kgH2O/m2/day but estimated here in MgH2O/m2/day
    ! therefore we must convert
    underflow = underflow*1d-3

    ! estimate potential drainage rate for the current time period
    liquid = soil_waterfrac(1:nos_soil_layers) * ( 1d0-iceprop(1:nos_soil_layers) )
    ! estimate how much liquid is available to flow
    avail_to_flow = liquid-field_capacity(1:nos_soil_layers)
    ! trapezium rule scaler and the half-way point between current and field capacity
    dx = avail_to_flow*0.5d0; halfway = liquid-dx
    do t = 1, nos_soil_layers
       if (avail_to_flow(t) > 0d0) then
           ! Trapezium rule for approximating integral of drainage rate
           call calculate_soil_conductivity(t, liquid(t), tmp1)
           call calculate_soil_conductivity(t, field_capacity(t), tmp2)
           call calculate_soil_conductivity(t, halfway(t), tmp3)
           pot_drainage(t) = 0.5d0*dx(t) * ((tmp1+tmp2) + 2d0*tmp3)
       else
           ! We are at field capacity currently even after rainfall has been infiltrated.
           ! Assume that the potential drainage rate is that at field capacity
           call calculate_soil_conductivity(t, field_capacity(t), pot_drainage(t))
       endif  ! water above field capacity to flow?
    end do  ! soil layers
    ! Scale potential drainage from per second to per day
    pot_drainage = pot_drainage*seconds_per_day

    ! Integrate drainage over each day until time period has been reached or
    ! each soil layer has reached field capacity
    t = 1
    do while (t < (time_period_days+1) .and. maxval(soil_waterfrac-field_capacity) > vsmall)

       ! Estimate liquid content and how much is available to flow/drain
       avail_to_flow = ( soil_waterfrac(1:nos_soil_layers) * (1d0-iceprop(1:nos_soil_layers)) ) &
                     - field_capacity(1:nos_soil_layers)

       ! ...then from the top down
       do soil_layer = 1, nos_soil_layers

          ! initial conditions; i.e. is there liquid water and more water than
          ! layer can hold
          if (avail_to_flow(soil_layer) > 0d0 .and. soil_waterfrac(soil_layer+1) < porosity(soil_layer+1)) then

              ! Unsaturated volume of layer below (m3 m-2)
              unsat = ( porosity(soil_layer+1) - soil_waterfrac(soil_layer+1) ) &
                    * layer_thickness(soil_layer+1) / layer_thickness(soil_layer)
              ! Restrict potential rate calculate above for the available water
              ! and available space in the layer below.
              ! NOTE: * layer_thickness(soil_layer) converts units from m3/m2 -> (m3)
              change = min(unsat, min(pot_drainage(soil_layer), avail_to_flow(soil_layer))) * layer_thickness(soil_layer)
              ! update soil layer below with drained liquid
              waterchange( soil_layer+1 ) = waterchange( soil_layer+1 ) + change
              waterchange( soil_layer     ) = waterchange( soil_layer     ) - change

          end if  ! some liquid water and drainage possible

       end do  ! soil layers

       ! update soil water profile
       soil_waterfrac(1:nos_soil_layers) = soil_waterfrac(1:nos_soil_layers) &
                                         + (waterchange(1:nos_soil_layers)/layer_thickness(1:nos_soil_layers))
       ! estimate drainage from bottom of soil column (MgH2O/m2/day)
       ! NOTES: that underflow is reset outside of the daily soil loop
       underflow = underflow+waterchange(nos_soil_layers+1)

       ! Reset now we have moves that liquid
       waterchange = 0d0
       ! integerate through time period
       t = t+1

    end do  ! while condition

    ! convert underflow from MgH2O/m2/day -> kgH2O/m2/day
    underflow = underflow*1d3

  end subroutine gravitational_drainage
  !
  !-----------------------------------------------------------------
  !
  subroutine soil_porosity(soil_frac_clay, soil_frac_sand)

   ! Porosity is estimated from Saxton equations. !

    implicit none

    ! arguments
    double precision, dimension(nos_soil_layers):: soil_frac_clay &
                                                   ,soil_frac_sand
    ! local variables..
    integer:: i
    double precision, parameter:: H = 0.332d0, &
                                   J = -7.251d-4, &
                                   K = 0.1276d0

    ! loop over soil layers..
    porosity(1:nos_soil_layers) = H+J * soil_frac_sand(1:nos_soil_layers) + &
                                  K*log10(soil_frac_clay(1:nos_soil_layers))
    ! then assign same to core layer
    porosity(nos_soil_layers+1) = porosity(nos_soil_layers)

  end subroutine soil_porosity
  !
  !---------------------------------------------------------------------
  !
  subroutine initialise_soils(soil_frac_clay, soil_frac_sand)

    !
    ! Subroutine calculate the soil layers field capacities and sets the initial
    ! soil water potential set to field capacity
    !

    implicit none

    ! arguments
    double precision, dimension(nos_soil_layers):: soil_frac_clay &
                                                   ,soil_frac_sand

    ! local variables
    integer:: i

    ! include some hardcoded boundaries for the Saxton equations
    where (soil_frac_sand < 5d0) soil_frac_sand = 5d0
    where (soil_frac_clay < 5d0) soil_frac_clay = 5d0
    where (soil_frac_clay > 60d0) soil_frac_clay = 60d0
    ! calculate soil porosity (m3/m3)
    call soil_porosity(soil_frac_clay, soil_frac_sand)
    ! calculate field capacity (m3/m-3)
    call calculate_field_capacity
    ! calculate initial soil water fraction
    soil_waterfrac = field_capacity
    ! calculate initial soil water potential
    SWP = dble_zero
    call soil_water_potential
    ! seperately calculate the soil conductivity as this applies to each layer
    do i = 1, nos_soil_layers
       call calculate_soil_conductivity(i, soil_waterfrac(i), soil_conductivity(i))
    end do  ! soil layers
    ! but apply the lowest soil layer to the core as well in initial conditions
    soil_conductivity(nos_soil_layers+1) = soil_conductivity(nos_soil_layers)

    ! final sanity check for porosity
    where (porosity <= (field_capacity+0.05d0)) porosity = field_capacity+0.05d0

  end subroutine initialise_soils
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_soil_conductivity(soil_layer, waterfrac, conductivity)

    ! Calculate the soil conductivity (m s-1) of water based on soil
    ! characteristics and current water content

    implicit none

    ! arguments
    integer, intent(in):: soil_layer
    double precision, intent(in):: waterfrac
    double precision, intent(out):: conductivity

    ! soil conductivity for the dynamic soil layers (i.e. not including core)
    conductivity = cond1(soil_layer) * exp(cond2(soil_layer)+cond3(soil_layer)/waterfrac)

    ! protection against floating point error
    if (waterfrac < 0.05d0) conductivity = 1d-30

  end subroutine calculate_soil_conductivity
  !
  !------------------------------------------------------------------
  !
  subroutine saxton_parameters(soil_frac_clay, soil_frac_sand)

    ! Calculate the key parameters of the Saxton, that is cond1, 2, 3 !
    ! and potA, B                                                    !

    implicit none

    ! arguments
    double precision, dimension(nos_soil_layers):: soil_frac_clay &
                                                   ,soil_frac_sand

    ! local variables
    double precision, parameter:: A = -4.396d0,  B = -0.0715d0,   CC = -4.880d-4, D = -4.285d-5, &
                                   E = -3.140d0,  F = -2.22d-3,     G = -3.484d-5, H = 0.332d0,   &
                                   J = -7.251d-4, K = 0.1276d0,     P = 12.012d0,  Q = -7.551d-2, &
                                   R = -3.895d0,  T = 3.671d-2,     U = -0.1103d0, V = 8.7546d-4, &
                                   mult1 = 100d0, mult2 = 2.778d-6

    ! layed out in this manor to avoid memory management issues in module
    ! variables
    potA(1:nos_soil_layers) = A + (B*soil_frac_clay) + &
                             (CC*soil_frac_sand*soil_frac_sand) + &
                              (D*soil_frac_sand*soil_frac_sand*soil_frac_clay)
    potA(1:nos_soil_layers) = exp(potA(1:nos_soil_layers))
    potA(1:nos_soil_layers) = potA(1:nos_soil_layers) * mult1

    potB(1:nos_soil_layers) = E + (F*soil_frac_clay*soil_frac_clay) + &
                                  (G*soil_frac_sand*soil_frac_sand*soil_frac_clay)

    cond1(1:nos_soil_layers) = mult2
    cond2(1:nos_soil_layers) = P + (Q*soil_frac_sand)
    cond3(1:nos_soil_layers) = R + (T*soil_frac_sand) + (U*soil_frac_clay) + &
                                   (V*soil_frac_clay*soil_frac_clay)

    ! assign bottom of soil column value to core
    potA(nos_soil_layers+1)  = potA(nos_soil_layers)
    potB(nos_soil_layers+1)  = potB(nos_soil_layers)
    cond1(nos_soil_layers+1) = mult2
    cond2(nos_soil_layers+1) = cond2(nos_soil_layers)
    cond3(nos_soil_layers+1) = cond3(nos_soil_layers)

  end subroutine saxton_parameters
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_soil_conductance(lm, local_lai)

    ! proceedsure to solve for soil surface resistance based on Monin-Obukov
    ! similarity theory stability correction momentum & heat are integrated
    ! through the under canopy space and canopy air space to the surface layer
    ! references are Nui & Yang 2004; Qin et al 2002
    ! NOTE: conversion to conductance at end

    implicit none

    ! declare arguments
    double precision, intent(in):: lm, local_lai

    ! local variables
    double precision:: canopy_decay & ! canopy decay coefficient for soil exchange
                       ,Kh_canht       ! eddy diffusivity at canopy height (m2.s-1)

    ! parameters
    double precision, parameter:: foliage_drag = 0.2d0  ! foliage drag coefficient

    ! calculate eddy diffusivity at the top of the canopy (m2.s-1)
    ! Kaimal & Finnigan 1994; for near canopy approximation
    Kh_canht = vonkarman*ustar*(canopy_height-displacement)

    ! calculate canopy decay coefficient with stability correction
    ! NOTE this is not consistent with canopy momentum decay done by Harman &
    ! Finnigan (2008)
    canopy_decay = sqrt((foliage_drag*canopy_height*local_lai)/lm)

    ! approximation of integral for soil resistance (s/m) and conversion to
    ! conductance (m/s)
    soil_conductance = ( canopy_height/(canopy_decay*Kh_canht) &
                       * (exp(canopy_decay*(1d0-(soil_roughl/canopy_height)))- &
                         exp(canopy_decay*(1d0-((roughl+displacement)/canopy_height)))) ) ** (-1d0)

  end subroutine calculate_soil_conductance
  !
  !----------------------------------------------------------------------
  !
  subroutine soil_water_potential

    ! Find SWP without updating waterfrac yet (we do that in !
    ! waterthermal). Waterfrac is m3 m-3, soilwp is MPa.     !

    implicit none

    integer:: i

    ! reformulation aims to remove if statement within loop to hopefully improve
    ! optimisation
    SWP(1:nos_soil_layers) = -0.001d0*potA(1:nos_soil_layers) &
                           * soil_waterfrac(1:nos_soil_layers)**potB(1:nos_soil_layers)
    ! NOTE: profiling indiates that 'where' is slower for very short vectors
    do i = 1, nos_soil_layers
       if (SWP(i) < -20d0 .or. SWP(i) /= SWP(i)) SWP(i) = -20d0
    end do

  end subroutine soil_water_potential
  !
  !------------------------------------------------------------------
  !
  subroutine update_net_radiation(isothermal, tempC, area_scaling, act_pot_ratio &
                                 ,sfc_exchange, aero_exchange, vapour_gradient, deltaTemp, deltaR)

    ! Use steady state solution of evaporation, convective (sensible) and
    ! radiative heat loss to update isothermal net radiation to net.
    ! Area scaling (e.g. lai) is an input to allow for common useage for soil
    ! (neglecting ground heat) and canopy. The key assumption here is that all
    ! values are make equivalent to ground area

    ! arguments
    double precision, intent(in) ::      tempC, & ! input surface/air temperature (oC)
                                    isothermal, & ! isothermal net radiation (SW+LW; W/m2)
                                  area_scaling, & ! area scaling to apply (m2/m2)
                                 act_pot_ratio, & ! ratio of potential to actual evaporation, i.e. (avail/potenial)
                                  sfc_exchange, & ! surface exchange conductance (m/s; e.g. stomatal conductance)
                                 aero_exchange, & ! aerodynamic exchange conductance (m/s; e.g. aerodynamic conductance)
                               vapour_gradient    ! vapour pressure gradient (kPa; either VPD or between air and soil)
    double precision, intent(out):: deltaTemp, & ! surface temperature difference (K)
                                     deltaR    ! surface longwave radiation difference (W/m2); subtract from isothermal longwave

    ! local variables
    double precision ::  tempK, & ! ambient temperature as K
          heat_loss_resistance, & ! resistance to heat loss from radiative and convection (s/m)
        aerodynamic_resistance, & ! aerodynamic resistance to water or heat exchangce (s/m)
           stomatal_resistance, & ! stomatal resistance to water exchangce (s/m)
              water_resistance, & ! serial combination of resistances to water evaporation
 thermal_gains, thermal_losses

    ! Ambient temperature C -> K
    tempK = tempC+freeze

    !
    ! Calculate resistance to heat loss (s/m)
    !

    ! First estimate radiative loss term, initially calculated as conductance)
    heat_loss_resistance = area_scaling*4d0*emissivity*boltz*tempK**3 / (air_density_kg*cpair)
    ! Combine in parallel radiative with convective conductances
    heat_loss_resistance = heat_loss_resistance + (2d0*aero_exchange)
    ! Convert from conductance m/s to s/m
    heat_loss_resistance = heat_loss_resistance ** (-1d0)

    !
    ! Convert aerodynamic and stomatal conductances to reisistance of water flux
    !

    aerodynamic_resistance = aero_exchange ** (-1d0)
    if (sfc_exchange == dble_zero) then
        ! if being used for surface water flux
        stomatal_resistance = dble_zero
    else
        ! if used for transpiration
        stomatal_resistance = sfc_exchange ** (-1d0)
    endif

    !
    ! Estimate thermal gains and losses (K) to calculate temperature difference
    !

    water_resistance = (aerodynamic_resistance+stomatal_resistance)
    thermal_gains = (heat_loss_resistance*water_resistance*psych*isothermal) &
                  / (air_density_kg*cpair * ((psych*water_resistance) + (slope*heat_loss_resistance)))
    thermal_losses = (heat_loss_resistance*vapour_gradient) &
                   / ((psych*water_resistance) + (slope*heat_loss_resistance))
    ! Determine surface temperature difference (K); should be added to the
    ! canopy temperature
    deltaTemp = thermal_gains-thermal_losses
    ! Apply actual potential ratio to scale wet surface evaporation when the
    ! supply of water is limited
    deltaTemp = deltaTemp*act_pot_ratio

    ! Estimate update between isothermal to net radiation (W/m2), including area
    ! correction
    ! note that this MUST be added from the longwave component outside of this
    ! function
    deltaR = -4d0*emissivity*boltz*tempK**3 * ( deltaTemp )

    ! return to user
    return

  end subroutine update_net_radiation
  !
  !------------------------------------------------------------------
  !
  subroutine z0_displacement(ustar_Uh, local_lai)

    ! dynamic calculation of roughness length and zero place displacement (m)
    ! based on canopy height and lai. Raupach (1994)

    implicit none

    ! arguments
    double precision, intent(out):: ustar_Uh  ! ratio of friction velocity over wind speed at canopy top
    double precision, intent(in):: local_lai
    ! local variables
    double precision  sqrt_cd1_lai
    double precision, parameter:: cd1 = 7.5d0,   & ! Canopy drag parameter; fitted to data
                                    Cs = 0.003d0, & ! Substrate drag coefficient
                                    Cr = 0.3d0,   & ! Roughness element drag coefficient
                          ustar_Uh_max = 0.3d0,   & ! Maximum observed ratio of
                                                    ! (friction velocity/canopy top wind speed) (m.s-1)
                          ustar_Uh_min = 0.05d0,  &
                                    Cw = 2d0,     &  ! Characterises roughness sublayer depth (m)
                                 phi_h = 0.19314718056d0  ! Roughness sublayer influence function; 

    ! describes the departure of the velocity profile from just above the
    ! roughness from the intertial sublayer log law


    ! Estimate canopy drag factor
    sqrt_cd1_lai = sqrt(cd1*local_lai)

    ! calculate estimate of ratio of friction velocity/canopy wind speed.
    ! NOTE: under current min LAI and fixed canopy height (9 m) this ratio is
    ! fixed at 0.3
!    ustar_Uh = 0.3d0
!    ustar_Uh = max(ustar_Uh_min, min(sqrt(Cs+Cr*local_lai*0.5d0), ustar_Uh_max))
    ustar_Uh = sqrt(Cs+Cr*local_lai*0.5d0)

    ! calculate displacement (m); assume minimum lai 1.0 or 1.5 as height is not
    ! varied
    displacement = (1d0-((1d0-exp(-sqrt_cd1_lai))/sqrt_cd1_lai))*canopy_height

    ! calculate roughness sublayer influence function; 
    ! this describes the departure of the velocity profile from just above the
    ! roughness from the intertial sublayer log law
    ! phi_h = log(Cw)-1d0+Cw**(-1d0)  ! DO NOT FORGET TO UPDATE IF Cw CHANGES

    ! finally calculate roughness length, dependant on displacement, friction
    ! velocity and lai.
    roughl = ((1d0-displacement/canopy_height)*exp(-vonkarman*ustar_Uh-phi_h))*canopy_height

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
  !------------------------------------------------------------------
  ! Functions other than the primary ACM and ACM ET are stored
  ! below this line.
  !------------------------------------------------------------------
  !
  !------------------------------------------------------------------
  !
  pure function arrhenious( a, b, t )

    ! The equation is simply...                        !
    !    a*exp( b * ( t-25.0 ) / ( t+273.15 ) )  !
    ! However, precision in this routine matters as it !
    ! affects many others. To maximise precision, the  !
    ! calculations have been split & d0 has been used. !

    implicit none

    ! arguments..
    double precision, intent(in):: a, b, t
    double precision            :: arrhenious

    arrhenious = a*exp( b * (t-25d0) / (t+freeze) )

  end function arrhenious
  !
  !----------------------------------------------------------------------
  !
  double precision function opt_max_scaling( max_val, optimum, kurtosis, current )

    ! Estimates a 0-1 scaling based on a skewed guassian distribution with a
    ! given optimum, maximum and kurtosis. Minimum is assumed to be at infinity
    ! (or near enough)

    implicit none

    ! arguments..
    double precision, intent(in):: max_val, optimum, kurtosis, current

    ! local variables..
    double precision, parameter:: min_val = -1d6

    ! Code with implicit assumption of min bound at infinity
!    if ( current >= max_val ) then
!         opt_max_scaling = 0d0
!    else
!         dummy     = exp( log((max_val-current) / (max_val-optimum)) * kurtosis * (max_val-optimum) )
!         opt_max_scaling = dummy*exp( kurtosis * ( current-optimum ) )
!    end if

    ! Code with explicit min bound
    if (current >= max_val .or. current <= min_val) then
        opt_max_scaling = 0d0
    else
        opt_max_scaling = exp( kurtosis*log((max_val-current)/(max_val-optimum)) * (max_val-optimum) ) &
                        * exp( kurtosis*log((current-min_val)/(optimum-min_val)) * (optimum-min_val) )
    endif


  end function opt_max_scaling
  !
  !------------------------------------------------------------------
  !
  double precision function root_resistance(root_mass, thickness)

   !
   ! Calculates root hydraulic resistance (MPa m2 s mmol-1) in a soil-root zone
   !

   implicit none

   ! arguments
   double precision:: root_mass, & ! root biomass in layer (gbiomass)
                       thickness    ! thickness of soil zone roots are in

   ! calculate root hydraulic resistance
   root_resistance = root_resist / (root_mass*thickness)

   ! return
   return

  end function root_resistance
  !
  !-----------------------------------------------------------------
  !
  double precision function soil_resistance(root_length, thickness, soilC)

    !
    ! Calculates the soil hydraulic resistance (MPa m2 s mmol-1) for a given
    ! soil-root zone
    !

    implicit none

    ! arguments
    double precision:: root_length, & ! root length in soil layer (m)
                          thickness, & ! thickness of soil layer (m)
                              soilC    ! soil conductivity m2.s-1.MPa-1

    ! local variables
    double precision:: rs, rs2

    ! calculate
    rs  = (root_length*pi)**(-0.5d0)
    rs2 = log( rs*root_radius_1 ) / (two_pi*root_length*thickness*soilC)
    ! soil water resistance
    soil_resistance = rs2*1d-9*mol_to_g_water

    ! return
    return

  end function soil_resistance
  !
  !------------------------------------------------------------------
  !
  double precision function plant_soil_flow(root_layer, root_length, root_mass &
                                           ,demand, root_reach_in, transpiration_resistance)

   !
   ! Calculate soil layer specific water flow form the soil to canopy
   ! (mmolH2O.m-2.s-1)
   ! Accounting for soil, root and plant resistance, and canopy demand
   !

   ! calculate and accumulate steady state water flux in mmol.m-2.s-1
   ! From the current soil layer given an amount of root within the soil layer.

   implicit none

   ! arguments
   integer, intent(in):: root_layer
   double precision, intent(in):: root_length, &
                                     root_mass, &
                                        demand, &
                                 root_reach_in, &
                      transpiration_resistance

   ! local arguments
   double precision:: soilR1, &
                       soilR2

   ! soil conductivity converted from m.s-1 -> m2.s-1.MPa-1 by head
   soilR1 = soil_resistance(root_length, root_reach_in, soil_conductivity(root_layer)*head_1)
   soilR2 = root_resistance(root_mass, root_reach_in)
   plant_soil_flow = demand/(transpiration_resistance+soilR1+soilR2)

   ! return
   return

  end function plant_soil_flow
  !
  !------------------------------------------------------------------
  !
  double precision function water_retention_saxton_eqns( xin )

    ! field capacity calculations for saxton eqns !

    implicit none

    ! arguments..
    double precision, intent(in):: xin

    ! local variables..
    double precision:: soil_wp

!    ! calculate the soil water potential (kPa)..
!    soil_WP = -0.001*potA( water_retention_pass ) * xin**potB( water_retention_pass )
!    water_retention_saxton_eqns = 1000. * soil_wp+10.    ! 10 kPa represents air-entry swp
    ! calculate the soil water potential (kPa)..
    soil_wp = -potA( water_retention_pass ) * xin**potB( water_retention_pass )
    water_retention_saxton_eqns = soil_wp+10d0    ! 10 kPa represents air-entry swp

    return

  end function water_retention_saxton_eqns
end module CARBON_MODEL_MOD
