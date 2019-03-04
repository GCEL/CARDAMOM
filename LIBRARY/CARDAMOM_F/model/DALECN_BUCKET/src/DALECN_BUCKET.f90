
module CARBON_MODEL_MOD

  implicit none

  ! make all private
  private

  ! explicit publics
  public :: CARBON_MODEL                  &
           ,vsmall                        &
           ,arrhenious                    &
           ,Rm_reich_N,Rm_reich_Q10       &
           ,acm_gpp                       &
           ,calculate_transpiration       &
           ,calculate_wetcanopy_evaporation &
           ,calculate_soil_evaporation    &
           ,acm_albedo_gc                 &
           ,meteorological_constants      &
           ,calculate_shortwave_balance   &
           ,calculate_longwave_isothermal &
           ,calculate_daylength           &
           ,opt_max_scaling               &
           ,freeze                        &
           ,co2comp_saturation            &
           ,co2comp_half_sat_conc         &
           ,co2_half_saturation           &
           ,co2_compensation_point        &
           ,pn_airt_scaling_time          &
           ,kc_saturation                 &
           ,kc_half_sat_conc              &
           ,calculate_update_soil_water   &
           ,calculate_Rtot                &
           ,calculate_aerodynamic_conductance &
           ,saxton_parameters             &
           ,initialise_soils              &
           ,update_soil_initial_conditions&
           ,linear_model_gradient         &
           ,seconds_per_day               &
           ,seconds_per_step              &
           ,root_biomass                  &
           ,root_reach                    &
           ,min_root                      &
           ,max_depth                     &
           ,root_k                        &
           ,top_soil_depth                &
           ,mid_soil_depth                &
           ,previous_depth                &
           ,nos_root_layers               &
           ,wSWP                          &
           ,SWP                           &
           ,SWP_initial                   &
           ,deltat_1                      &
           ,water_flux                    &
           ,layer_thickness               &
           ,waterloss,watergain           &
           ,potA,potB                     &
           ,cond1,cond2,cond3             &
           ,soil_conductivity             &
           ,soil_waterfrac                &
           ,soil_waterfrac_initial        &
           ,porosity                      &
           ,porosity_initial              &
           ,field_capacity                &
           ,field_capacity_initial        &
           ,drythick                      &
           ,min_drythick                  &
           ,min_layer                     &
           ,soilwatermm                   &
           ,wSWP_time                     &
           ,soil_frac_clay                &
           ,soil_frac_sand                &
           ,nos_soil_layers               &
           ,meant                         &
           ,meant_time                    &
           ,stomatal_conductance          &
           ,aerodynamic_conductance       &
           ,iWUE                          &
           ,avN                           &
           ,NUE_optimum                   &
           ,NUE                           &
           ,pn_max_temp                   &
           ,pn_opt_temp                   &
           ,pn_kurtosis                   &
           ,e0                            &
           ,co2_half_sat                  &
           ,pn_airt_scaling               &
           ,co2_comp_point                &
           ,minlwp                        &
           ,max_lai_lwrad_transmitted     &
           ,lai_half_lwrad_transmitted    &
           ,max_lai_nir_reflection        &
           ,lai_half_nir_reflection       &
           ,max_lai_par_reflection        &
           ,lai_half_par_reflection       &
           ,max_lai_par_transmitted       &
           ,lai_half_par_transmitted      &
           ,max_lai_nir_transmitted       &
           ,lai_half_nir_transmitted      &
           ,max_lai_lwrad_reflected       &
           ,lai_half_lwrad_reflected      &
           ,soil_swrad_absorption         &
           ,max_lai_lwrad_release         &
           ,lai_half_lwrad_release        &
           ,leafT                         &
           ,mint                          &
           ,maxt                          &
           ,swrad                         &
           ,co2                           &
           ,doy                           &
           ,rainfall                      &
           ,airt_zero_fraction            &
           ,snowfall                      &
           ,snow_melt                     &
           ,wind_spd                      &
           ,vpd_kPa                       &
           ,lai                           &
           ,days_per_step                 &
           ,days_per_step_1               &
           ,dayl_seconds_1                &
           ,dayl_seconds                  &
           ,dayl_hours                    &
           ,snow_storage                  &
           ,canopy_storage                &
           ,intercepted_rainfall          &
           ,disturbance_residue_to_litter &
           ,disturbance_residue_to_cwd    &
           ,disturbance_residue_to_som    &
           ,disturbance_loss_from_litter  &
           ,disturbance_loss_from_cwd     &
           ,disturbance_loss_from_som     &
           ,rainfall_time                 &
           ,Cwood_labile_release_coef     &
           ,Croot_labile_release_coef     &
           ,canopy_days                   &
           ,canopy_age_vector             &
           ,dim_1,dim_2                   &
           ,nos_trees                     &
           ,nos_inputs                    &
           ,leftDaughter                  &
           ,rightDaughter                 &
           ,nodestatus                    &
           ,xbestsplit                    &
           ,nodepred                      &
           ,bestvar

  !!!!!!!!!!
  ! Random Forest GPP emulator
  !!!!!!!!!!

  ! arrays for the emulator, just so we load them once and that is it cos they be
  ! massive
  integer ::    dim_1, & ! dimension 1 of response surface
                dim_2, & ! dimension 2 of response surface
            nos_trees, & ! number of trees in randomForest
           nos_inputs    ! number of driver inputs
  double precision, allocatable, dimension(:,:) :: leftDaughter, & ! left daughter for forest
                                                  rightDaughter, & ! right daughter for forets
                                                     nodestatus, & ! nodestatus for forests
                                                     xbestsplit, & ! for forest
                                                       nodepred, & ! prediction value for each tree
                                                        bestvar    ! for randomForests

  !!!!!!!!!
  ! Parameters
  !!!!!!!!!

  ! useful technical parameters
  logical :: do_iWUE = .true. ! Use iWUE or WUE for stomatal optimisation
  double precision, parameter :: vsmall = tiny(0d0)*1d3 ! *1d3 to add a little breathing room

  integer, parameter :: nos_root_layers = 3, nos_soil_layers = nos_root_layers + 1
  double precision, parameter :: pi = 3.1415927d0,  &
                               pi_1 = 0.3183099d0,  & ! pi**(-1d0)
                                pi2 = 9.869604d0,   & ! pi**2d0
                             two_pi = 6.283185d0,   & ! pi*2d0
                         deg_to_rad = 0.01745329d0, & ! pi/180d0
                sin_dayl_deg_to_rad = 0.3979486d0,  & ! sin( 23.45d0 * deg_to_rad )
                            gravity = 9.8067d0,     & ! acceleration due to gravity, ms-1
                              boltz = 5.670400d-8,  & ! Boltzmann constant (W.m-2.K-4)
                         emissivity = 0.96d0,       &
                        emiss_boltz = 5.443584d-08, & ! emissivity * boltz
                    sw_par_fraction = 0.5d0,        & ! fraction of short-wave radiation which is PAR
                             freeze = 273.15d0,     &
                         gs_H2O_CO2 = 1.646259d0,   & ! The ratio of H20:CO2 diffusion for gs (Jones appendix 2)
                       gs_H2O_CO2_1 = 0.6074378d0,  & ! gs_H2O_CO2 ** (-1d0), &
                  gs_H2Ommol_CO2mol = 0.001646259d0,& ! gs_H2O_CO2 * 1d-3
                         gb_H2O_CO2 = 1.37d0,       & ! The ratio of H20:CO2 diffusion for gb (Jones appendix 2)
            partial_molar_vol_water = 18.05d-6,     & ! partial molar volume of water, m3 mol-1 at 20C
                     mol_to_g_water = 18d0,         & ! molecular mass of water
                   mmol_to_kg_water = 1.8d-5,       & ! milli mole conversion to kg
                       mol_to_g_co2 = 12d0,         & ! molecular mass of CO2 (g)
                         umol_to_gC = 1.2d-5,       & ! conversion of umolC -> gC
                         gC_to_umol = 83333.33d0,   & ! conversion of gC -> umolC; umol_to_gC**(-1d0)
                       g_to_mol_co2 = 0.08333333d0, &
  !snowscheme       density_of_water = 998.9d0,         & ! density of !water kg.m-3
                     gas_constant_d = 287.04d0,     & ! gas constant for dry air (J.K-1.mol-1)
                               Rcon = 8.3144d0,     & ! Universal gas constant (J.K-1.mol-1)
                          vonkarman = 0.41d0,       & ! von Karman's constant
                        vonkarman_1 = 2.439024d0,   & ! 1 / von Karman's constant
                        vonkarman_2 = 0.1681d0,     & ! von Karman's constant^2
                              cpair = 1004.6d0        ! Specific heat capacity of air; used in energy balance J.kg-1.K-1

  ! photosynthesis / respiration parameters
  double precision, parameter :: &
                      kc_saturation = 310d0,        & ! CO2 half saturation, at reference temperature (298.15 K)
                   kc_half_sat_conc = 23.956d0,     & ! CO2 half sat, sensitivity coefficient
                 co2comp_saturation = 36.5d0,       & ! CO2 compensation point, at reference temperature (298.15 K)
              co2comp_half_sat_conc = 9.46d0,       & ! CO2 comp point, sensitivity coefficient
                                                      ! Each of these are temperature
                                                      ! sensitivty
                leaf_life_weighting = 0.5d0,        & ! inverse of averaging period of lagged effects
                                                      ! probably should be an actual parmeter
                        Rg_fraction = 0.21875d0,    & ! fraction of C allocation towards each pool
                                                      ! lost as growth respiration
                                                      ! (i.e. 0.28 .eq. xNPP)
                    one_Rg_fraction = 1d0 - Rg_fraction

  ! hydraulic parameters
  double precision, parameter :: &
                         tortuosity = 2.5d0,        & ! tortuosity
                             gplant = 5d0,          & ! plant hydraulic conductivity (mmol m-1 s-1 MPa-1)
                        root_resist = 25d0,         & ! Root resistivity (MPa s g mmol−1 H2O)
                        root_radius = 0.00029d0,    & ! root radius (m) Bonen et al 2014 = 0.00029
                                                      ! Williams et al 1996 = 0.0001
                      root_radius_1 = root_radius**(-1d0), &
                root_cross_sec_area = pi * root_radius**2, & ! root cross sectional area (m2)
                                                             ! = pi * root_radius * root_radius
                       root_density = 0.31d6,       & ! root density (g biomass m-3 root)
                                                      ! 0.5e6 Williams et al 1996
                                                      ! 0.31e6 Bonan et al 2014
            root_mass_length_coef_1 = (root_cross_sec_area * root_density)**(-1d0), &
                 const_sfc_pressure = 101325d0,     & ! (Pa)  Atmospheric surface pressure
                               head = 0.009807d0,   & ! head of pressure (MPa/m)
                             head_1 = 101.968d0       ! inverse head of pressure (m/MPa)

  ! structural parameters
  double precision, parameter :: &
                      canopy_height = 9d0,          & ! canopy height assumed to be 9 m
                       tower_height = canopy_height + 2d0, & ! tower (observation) height assumed to be 2 m above canopy
                           min_wind = 0.1d0,        & ! minimum wind speed at canopy top
                       min_drythick = 0.01d0,       & ! minimum dry thickness depth (m)
                          min_layer = 0.03d0,       & ! minimum thickness of the third rooting layer (m)
                        soil_roughl = 0.05d0,       & ! soil roughness length (m)
                     top_soil_depth = 0.1d0,        & ! thickness of the top soil layer (m)
                     mid_soil_depth = 0.2d0,        & ! thickness of the second soil layer (m)
                           min_root = 5d0,          & ! minimum root biomass (gBiomass.m-2)
                            min_lai = 1.5d0,        & ! minimum LAI assumed for aerodynamic conductance calculations (m2/m2)
                    min_throughfall = 0.2d0,        & ! minimum fraction of precipitation which
                                                      ! is through fall
                        min_storage = 0.2d0           ! minimum canopy water (surface) storage (mm)

  ! timing parameters
  double precision, parameter :: &
                   seconds_per_hour = 3600d0,         & ! Number of seconds per hour
                    seconds_per_day = 86400d0,        & ! Number of seconds per day
                  seconds_per_day_1 = 1.157407d-05      ! Inverse of seconds per day

  ! ACM-GPP-ET parameters
  double precision, parameter :: &
                        pn_max_temp = 5.357174d+01, & ! Maximum temperature for photosynthesis (oC)
                        pn_opt_temp = 3.137242d+01, & ! Optimum temperature for photosynthesis (oC)
                        pn_kurtosis = 1.927458d-01, & ! Kurtosis of photosynthesis temperature response
                                 e0 = 5.875662d+00, & ! Quantum yield gC/MJ/m2/day PAR
          max_lai_lwrad_transmitted = 7.626683d-01, & ! Max fractional reduction of LW from sky transmitted through canopy
         lai_half_lwrad_transmitted = 7.160363d-01, & ! LAI at which canopy LW transmittance reduction = 50 %
             max_lai_nir_reflection = 4.634860d-01, & ! Max fraction of NIR reflected by canopy
            lai_half_nir_reflection = 1.559148d+00, & ! LAI at which canopy NIR reflected = 50 %
                             minlwp =-1.996830d+00, & ! minimum leaf water potential (MPa)
             max_lai_par_reflection = 1.623013d-01, & ! Max fraction of PAR reflected by canopy
            lai_half_par_reflection = 1.114360d+00, & ! LAI at which canopy PAR reflected = 50 %
           lai_half_lwrad_reflected = 1.126214d+00, & ! LAI at which 50 % LW is reflected back to sky
                               iWUE = 1.602503d-06, & ! Intrinsic water use efficiency (gC/m2leaf/day/mmolH2Ogs)
              soil_swrad_absorption = 6.643079d-01, & ! Fraction of SW rad absorbed by soil
            max_lai_par_transmitted = 8.079519d-01, & ! Max fractional reduction in PAR transmittance by canopy
           lai_half_par_transmitted = 9.178784d-01, & ! LAI at which PAR transmittance reduction = 50 %
            max_lai_nir_transmitted = 8.289803d-01, & ! Max fractional reduction in NIR transmittance by canopy
           lai_half_nir_transmitted = 1.961831d+00, & ! LAI at which NIR transmittance reduction = 50 %
              max_lai_lwrad_release = 9.852855d-01, & ! Max fraction of LW emitted (1-par) from canopy to be released
             lai_half_lwrad_release = 7.535450d-01, & ! LAI at which LW emitted from canopy to be released at 50 %
            max_lai_lwrad_reflected = 1.955832d-02    ! LAI at which 50 % LW is reflected back to sky

  !!!!!!!!!
  ! Module level variables
  !!!!!!!!!

  ! local variables for GSI phenology model
  double precision :: SLA & ! Specific leaf area
                     ,avail_labile,Rg_from_labile    &
                     ,Cfol_turnover_gradient         &
                     ,Cfol_turnover_half_saturation  &
                     ,Cwood_labile_release_gradient  &
                     ,Cwood_labile_half_saturation   &
                     ,Croot_labile_release_gradient  &
                     ,Croot_labile_half_saturation   &
                     ,Cwood_hydraulic_gradient       &
                     ,Cwood_hydraulic_half_saturation&
                     ,Cwood_hydraulic_limit          &
                     ,tmp,gradient

  double precision, allocatable, dimension(:) :: disturbance_residue_to_litter, &
                                                 disturbance_residue_to_som,    &
                                                 disturbance_residue_to_cwd,    &
                                                 disturbance_loss_from_litter,  &
                                                 disturbance_loss_from_cwd,     &
                                                 disturbance_loss_from_som

  ! Autotrophic respiration model / phenological choices
  ! See source below for details of these variables
  integer :: oldest_leaf, youngest_leaf
  double precision :: deltaGPP, deltaRm, Rm_deficit, &
                      leaf_growth_period,      &
                      leaf_growth_period_1,    &
                      leaf_mortality_period,   &
                      leaf_mortality_period_1, &
                      marginal_gain_avg,       &
                      Q10_adjustment, &
                      Rm_deficit_leaf_loss, &
                      Rm_deficit_root_loss, &
                      Rm_deficit_wood_loss, &
                      Rm_leaf, Rm_root, Rm_wood, &
                      CN_leaf, CN_root, CN_wood, &
                      root_cost,root_life,       &
                      leaf_cost,leaf_life,       &
                      wood_cost,canopy_age,      &
                      leaf_life_max,             &
                      canopy_maturation_lag,     &
                      Rm_leaf_baseline,          &
                      Rm_root_baseline,          &
                      Rm_wood_baseline

  ! Declare canopy age class related variables
  ! NOTE: 5480 is maximum number of years we would expect a canopy to every last
  ! Declared as fixed value to avoid cost of allocation
  integer, dimension(5480) :: leaf_loss_possible
  double precision, dimension(5480) :: canopy_age_vector, &
                                  canopy_days,NUE_vector, &
                                       NUE_vector_mature, &
                                       marginal_loss_avg

  ! hydraulic model variables
  integer :: water_retention_pass, soil_layer
  double precision, dimension(nos_soil_layers) :: soil_frac_clay,soil_frac_sand ! clay and soil fractions of soil
  double precision, dimension(nos_root_layers) :: uptake_fraction, & ! fraction of water uptake from each root layer
                                                           demand, & ! maximum potential canopy hydraulic demand
                                                       water_flux    ! potential transpiration flux (mmolH2O.m-2.s-1)
  double precision, dimension(nos_soil_layers+1) :: SWP, & ! soil water potential (MPa)
                                            SWP_initial, &
                                      soil_conductivity, & ! soil conductivity
                                              waterloss, & ! water loss from specific soil layers (m)
                                              watergain, & ! water gained by specfic soil layers (m)
                                         field_capacity, & ! soil field capacity (m3.m-3)
                                 field_capacity_initial, &
                                         soil_waterfrac, & ! soil water content (m3.m-3)
                                 soil_waterfrac_initial, &
                                               porosity, & ! soil layer porosity, (fraction)
                                       porosity_initial, &
                                        layer_thickness, & ! thickness of soil layers (m)
                        cond1, cond2, cond3, potA, potB    ! Saxton equation values

  double precision :: root_reach, root_biomass, &
                                      drythick, & ! estimate of the thickness of the dry layer at soil surface (m)
                                          wSWP, & ! weighted soil water potential (MPa) used in GSI calculate.
                                                  ! Removes / limits the fact that very low root density in young plants
                                                  ! give values too large for GSI to handle.
                                     max_depth, & ! maximum possible root depth (m)
                                        root_k, & ! biomass to reach half max_depth
                                        runoff, & ! runoff (kgH2O.m-2.day-1)
                                     underflow, & ! drainage from the bottom of soil column (kgH2O.m-2.day-1)
                      new_depth,previous_depth, & ! depth of bottom of soil profile
                                   canopy_wind, & ! wind speed (m.s-1) at canopy top
                                         ustar, & ! friction velocity (m.s-1)
                                      ustar_Uh, &
                                air_density_kg, & ! air density kg/m3
                                ET_demand_coef, & ! air_density_kg * vpd_kPa * cpair
                                        roughl, & ! roughness length (m)
                                  displacement, & ! zero plane displacement (m)
                                    max_supply, & ! maximum water supply (mmolH2O/m2/day)
                                         meant, & ! mean air temperature (oC)
                                     maxt_lag1, &
                                         leafT, & ! canopy temperature (oC)
                              mean_annual_temp, &
                            canopy_swrad_MJday, & ! canopy_absorbed shortwave radiation (MJ.m-2.day-1)
                              canopy_par_MJday, & ! canopy_absorbed PAR radiation (MJ.m-2.day-1)
                              soil_swrad_MJday, & ! soil absorbed shortwave radiation (MJ.m-2.day-1)
                              canopy_lwrad_Wm2, & ! canopy absorbed longwave radiation (W.m-2)
                                soil_lwrad_Wm2, & ! soil absorbed longwave radiation (W.m-2)
                                 sky_lwrad_Wm2, & ! sky absorbed longwave radiation (W.m-2)
                                     ci_global, & ! internal CO2 concentration (ppm or umol/mol)
                          stomatal_conductance, & ! maximum stomatal conductance (mmolH2O.m-2.s-1)
                       aerodynamic_conductance, & ! bulk surface layer conductance (m.s-1)
                              soil_conductance, & ! soil surface conductance (m.s-1)
                             convert_ms1_mol_1, & ! Conversion ratio for m.s-1 -> mol.m-2.s-1
                            convert_ms1_mmol_1, & ! Conversion ratio for m/s -> mmol/m2/s
                           air_vapour_pressure, & ! Vapour pressure of the air (kPa)
                                        lambda, & ! latent heat of vapourisation (J.kg-1)
                                         psych, & ! psychrometric constant (kPa K-1)
                                         slope, & ! Rate of change of saturation vapour pressure with temperature (kPa.K-1)
                        water_vapour_diffusion, & ! Water vapour diffusion coefficient in (m2/s)
                             dynamic_viscosity, & ! dynamic viscosity (kg.m-2.s-1)
                           kinematic_viscosity, & ! kinematic viscosity (m2.s-1)
                                  snow_storage, & ! snow storage (kgH2O/m2)
                                canopy_storage, & ! water storage on canopy (kgH2O.m-2)
                          intercepted_rainfall    ! intercepted rainfall rate equivalent (kgH2O.m-2.s-1)

  ! Module level variables for ACM_GPP_ET parameters
  double precision :: delta_gs, & ! day length corrected gs increment mmolH2O/m2/dayl
                           avN, & ! average foliar N (gN/m2)
                      NUE_mean, &
                  NUE_mean_lag, &
                   NUE_optimum, & !
                           NUE, & ! Photosynthetic nitrogen use efficiency at optimum temperature (oC)
                                  ! ,unlimited by CO2, light and photoperiod (gC/gN/m2leaf/day)
               pn_airt_scaling, & ! temperature response for metabolic limited photosynthesis
                  co2_half_sat, & ! CO2 at which photosynthesis is 50 % of maximum (ppm)
                co2_comp_point    ! CO2 at which photosynthesis > 0 (ppm)

  ! Module level variables for step specific met drivers
  double precision :: mint, & ! minimum temperature (oC)
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
  double precision :: cos_solar_zenith_angle, &
                          mean_days_per_step, & !
                            seconds_per_step, & !
                               days_per_step, & !
                             days_per_step_1, & !
                                dayl_seconds, & ! day length in seconds
                              dayl_seconds_1, &
                                  dayl_hours    ! day length in hours

  double precision, dimension(:), allocatable :: deltat_1, & ! inverse of decimal days
                                               meant_time, &
                                  airt_zero_fraction_time, &
                                          daylength_hours, &
                                        daylength_seconds, &
                                      daylength_seconds_1, &
                                            rainfall_time, &
                                      co2_half_saturation, & ! (ppm)
                                   co2_compensation_point, & ! (ppm)
                                     pn_airt_scaling_time, &
                                      air_density_kg_time, &
                                   convert_ms1_mol_1_time, &
                                   lambda_time,psych_time, &
                                               slope_time, &
                                      ET_demand_coef_time, &
                              water_vapour_diffusion_time, &
                                   dynamic_viscosity_time, &
                                 kinematic_viscosity_time, &
                                Cwood_labile_release_coef, & ! time series of labile release to wood
                                Croot_labile_release_coef, & ! time series of labile release to root
                                       Cfol_turnover_coef, &
                                              soilwatermm, &
                                                wSWP_time

  double precision :: marginal_output, gradient_output
  save

contains
  !
  !--------------------------------------------------------------------
  !
  subroutine CARBON_MODEL(start,finish,met,pars,deltat,nodays,lat,lai_out,NEE_out,FLUXES,POOLS &
                         ,nopars,nomet,nopools,nofluxes,GPP_out)

    !
    ! The Data Assimilation Linked Ecosystem Carbon&Nitrogen - BUCKET (DALECN_BUCKET) model.
    !
    ! The Aggregated Canopy Model for Gross Primary Productivity and Evapotranspiration (ACM-GPP-ET)
    ! simulates coupled photosynthesis-transpiration (via stomata), soil and intercepted canopy evaporation and
    ! soil water balance (4 layers).
    !
    ! Carbon allocation to tissues are determined by partitioning parameters under marginal return calculations. Autotrophic
    ! respiration is divided between growth (fixed fraction of new tissue) and maintenance (determine by Reich et al 2008).
    !
    ! This version was coded by T. Luke Smallman (t.l.smallman@ed.ac.uk).
    ! Version 1.0: 05/08/2018
    ! Version 1.1: 17/02/2019 - Inclusion of GSI type model in decline of NUE

    implicit none

    ! declare input variables
    integer, intent(in) :: start    &
                          ,finish   &
                          ,nopars   & ! number of paremeters in vector
                          ,nomet    & ! number of meteorological fields
                          ,nofluxes & ! number of model fluxes
                          ,nopools  & ! number of model pools
                          ,nodays     ! number of days in simulation

    double precision, intent(in) :: met(nomet,nodays) & ! met drivers
                                   ,deltat(nodays)    & ! time step in decimal days
                                   ,pars(nopars)      & ! number of parameters
                                   ,lat                 ! site latitude (degrees)

    double precision, dimension(nodays), intent(inout) :: lai_out & ! leaf area index
                                                         ,GPP_out & ! Gross primary productivity
                                                         ,NEE_out   ! net ecosystem exchange of CO2

    double precision, dimension((nodays+1),nopools), intent(inout) :: POOLS ! vector of ecosystem pools

    double precision, dimension(nodays,nofluxes), intent(inout) :: FLUXES ! vector of ecosystem fluxes

    ! declare general local variables
    double precision ::  tmp,tmp1,infi &
                        ,transpiration &
                      ,soilevaporation &
                       ,wetcanopy_evap &
                     ,snow_sublimation &
                              ,deltaWP & ! deltaWP (MPa) minlwp-soilWP
                        ,act_pot_ratio &
                                 ,loss &
                                 ,Rtot   ! Total hydraulic resistance (MPa.s-1.m-2.mmol-1)

    integer :: nxp,n,test,m,a,b,c

    ! local fire related variables
    double precision :: burnt_area          &
                           ,CFF(7) = 0d0 & ! combusted and non-combustion fluxes
                          ,NCFF(7) = 0d0 & ! with residue and non-residue seperates
                   ,combust_eff(5)       & ! combustion efficiency
                            ,rfac          ! fire resilience factor

    integer :: steps_per_year ! mean number of steps in a year

    ! local deforestation related variables
    double precision, dimension(5) :: post_harvest_burn   & ! how much burning to occur after
                                     ,foliage_frac_res    &
                                     ,roots_frac_res      &
                                     ,rootcr_frac_res     &
                                     ,stem_frac_res       &
                                     ,branch_frac_res     &
                                     ,Cbranch_part        &
                                     ,Crootcr_part        &
                                     ,soil_loss_frac

    double precision :: labile_loss,foliar_loss       &
                       ,roots_loss,wood_loss          &
                       ,labile_residue,foliar_residue &
                       ,roots_residue,wood_residue    &
                       ,C_total,labile_frac_res       &
                       ,Cstem,Cbranch,Crootcr         &
                       ,stem_residue,branch_residue   &
                       ,coarse_root_residue           &
                       ,soil_loss_with_roots

    integer :: reforest_day, harvest_management, restocking_lag

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
    ! 1 = labile (p18)
    ! 2 = foliar (p19)
    ! 3 = root   (p20)
    ! 4 = wood   (p21)
    ! 5 = litter (p22)
    ! 6 = som    (p23)
    ! 7 = cwd    (p24)
    ! 8 = soil water content (currently assumed to field capacity)

    ! p(30) = labile replanting
    ! p(31) = foliar replanting
    ! p(32) = fine root replanting
    ! p(33) = wood replanting

    ! FLUXES are:
    ! 1 = GPP (gC/m2/day)
    ! 2 = temprate
    ! 3 = respiration_auto (gC/m2/day)
    ! 4 = -----NOT IN USE-----
    ! 5 = labile production (gC/m2/day)
    ! 6 = root production (gC/m2/day)
    ! 7 = wood production (gC/m2/day)
    ! 8 = leaf production (gC/m2/day)
    ! 9 = -----NOT IN USE-----
    ! 10 = leaf litter production (gC/m2/day)
    ! 11 = woodlitter production (gC/m2/day)
    ! 12 = rootlitter production (gC/m2/day)
    ! 13 = respiration het litter (gC/m2/day)
    ! 14 = respiration het som (gC/m2/day)
    ! 15 = litter2som (gC/m2/day)
    ! 16 = -----NOT IN USE-----
    ! 17 = carbon flux due to fire (gC/m2/day)
    ! 18 = Mean Canopy age (days)
    ! 19 = Evapotranspiration (kgH2O.m-2.day-1)
    ! 20 = CWD turnover to litter (gC/m2/day)
    ! 21 = C extracted as harvest (gC/m2/day)
    ! 22 = labile loss due to disturbance (gC/m2/day)
    ! 23 = foliage loss due to disturbance (gC/m2/day)
    ! 24 = root loss due to disturbance (gC/m2/day)
    ! 25 = wood loss due to disturbance (gC/m2/day)

    ! PARAMETERS
    ! 41 process parameters; 7 C pool initial conditions; 1 soil water initial condition

    ! p(1) = Litter to SOM conversion rate (fraction)
    ! p(2) = CN_root (gC/gN)
    ! p(3) = Initial mean NUE
    ! p(4) = Max labile turnover to roots (fraction)
    ! p(5) = Leaf marginal growth sensitivity (days)
    ! p(6) = Turnover rate of wood (fraction)
    ! p(7) = Turnover rate of roots (fraction)
    ! p(8) = Litter turnover rate to heterotrophic respiration (fraction)
    ! p(9) = SOM turnover rate to heterotrophic respiration (fraction)
    ! p(10) = Exponential coefficient for temperature response for heterotrophic respiration
    ! p(11) = Average foliar nitrogen content (log10(gN/m2))
    ! p(12) = Max labile turnover to leaves (fraction)
    ! p(13) = Max labile turnover to wood (fraction)
    ! p(14) = Days after emergence at which canopy reaches optimum NUE
    ! p(15) = CN_wood (gC/gN)
    ! p(16) = CWD turnover fraction (fraction)
    ! p(17) = Leaf Mass per unit Area (gC.m-2)
    ! p(18) = Initial labile pool (gC/m2)
    ! p(19) = Initial foliage pool (gC/m2)
    ! p(20) = Initial root pool (gC/m2)
    ! p(21) = Initial wood pool (gC/m2)
    ! p(22) = Initial litter pool (gC/m2)
    ! p(23) = Initial som pool (gC/m2)
    ! p(24) = Initial CWD pool (gC/m2)
    ! p(25) = Initial canopy age (days)
    ! p(26) = Photosynthetic nitrogen use efficiency (gC/gN/m2/day)
    ! p(27) = Initial canopy life span (days)
    ! p(28) = Fraction of Cwood which is Cbranch
    ! p(29) = Fraction of Cwood which is Ccoarseroot
    ! p(34) = Fine root (gbiomass.m-2) needed to reach 50% of max depth
    ! p(35) = Maximum rooting depth (m)
    ! p(36) = Reich Rm_leaf N exponent
    ! p(37) = Reich Rm_leaf N baseline
    ! p(38) = Reich Rm_root N exponent
    ! p(39) = Reich Rm_root N baseline
    ! p(40) = Reich Rm_wood N exponent
    ! p(41) = Reich Rm_wood N baseline
    ! p(42) = Leaf marginal loss sensitivity (days)
    ! p(43) = iWUE
    ! p(44) = Initial root profile water content (m3/m3)
    ! p(45) = Period (days) over which the initial canopy biomass is distributed
    ! p(46) = Min temperature threshold for NUE decline
    ! p(47) = Max temperature threshold for NUE decline
    ! p(48) = Min VPD threshold for NUE decline
    ! p(49) = Max VPD threshold for NUE decline

    ! variables related to deforestation
    ! labile_loss = total loss from labile pool from deforestation
    ! foliar_loss = total loss form foliar pool from deforestation
    ! roots_loss = total loss from root pool from deforestation
    ! wood_loss = total loss from wood pool from deforestation
    ! labile_residue = harvested labile remaining in system as residue
    ! foliar_residue = harested foliar remaining in system as residue
    ! roots_residue = harvested roots remaining in system as residue
    ! wood_residue = harvested wood remaining in system as residue
    ! coarse_root_residue = expected coarse woody root left in system as residue

    ! parameters related to deforestation
    ! labile_frac_res = fraction of labile harvest left as residue
    ! foliage_frac_res = fraction of foliage harvest left as residue
    ! roots_frac_res = fraction of roots harvest left as residue
    ! wood_frac_res = fraction of wood harvest left as residue
    ! Crootcr_part = fraction of wood pool expected to be coarse root
    ! Crootcr_frac_res = fraction of coarse root left as residue
    ! soil_loss_frac = fraction determining Csom expected to be physically
    ! removed along with coarse roots

    ! profiling example
    !real :: begin, done, f1 = 0, f2 = 0, f3 = 0, f4 = 0, f5 = 0,total_time = 0
    !real :: Rtot_track_time = 0, aero_time = 0, soilwater_time = 0 , acm_et_time = 0 , Rm_time = 0
    !call cpu_time(begin)
    !call cpu_time(done)
    !print*,"2h"
    ! infinity check requirement
    infi = 0d0
    ! reset basic input / output variables
    FLUXES = 0d0 ; POOLS = 0d0

    ! load ACM-GPP-ET parameters
    ! WARNING: pn_max_temp, pn_opt_temp, pn_kurtosis used in EDC1!
    NUE_optimum = pars(26) ! Photosynthetic nitrogen use efficiency at optimum temperature (oC)
    ! ,unlimited by CO2, light and photoperiod (gC/gN/m2leaf/day)

    ! mean number of model steps per year
    mean_days_per_step = sum(deltat) ! sum nos days
    steps_per_year = nint(mean_days_per_step*0.002737851d0) ! 0.002737851 = 1/365.25
    steps_per_year = nint(mean_days_per_step/dble(steps_per_year))
    mean_days_per_step = nint(mean_days_per_step / dble(nodays)) ! now update to mean

    ! Parameters related to ACM-GPP-ET, but not actually parameters of the ACM-GPP-ET model
    avN = 10d0**pars(11)             ! Average foliar Nitrogen content gN/m2leaf
    deltaWP = minlwp                 ! leafWP-soilWP (i.e. -2-0 MPa)
    Rtot = 1d0                       ! Reset Total hydraulic resistance to 1
    canopy_maturation_lag = pars(14) ! canopy age (days) before peak NUE
    ! estimate the canopy growth sensitivity variable (i.e. period over which to average marginal returns)
    leaf_growth_period = ceiling(pars(5))/mean_days_per_step
    leaf_growth_period_1 = leaf_growth_period**(-1d0)
    ! estimate the canopy mortality sensitivity variable (i.e. period
    ! over which to average marginal returns)
    leaf_mortality_period = ceiling(pars(42))/mean_days_per_step
    leaf_mortality_period_1 = leaf_mortality_period**(-1d0)

    ! Root biomass to reach 50% (root_k) of maximum rooting depth (max_depth)
    root_k = pars(34) ; max_depth = pars(35)

    !!!!!!!!!!!!
    ! set time invarient / initial phenology parameters
    !!!!!!!!!!!!

    ! calculate specific leaf area from leaf mass area
    SLA = pars(17)**(-1d0)
    ! calculate root life spans (days)
    !root_life = pars(7)**(-1d0)
    ! Assign initial leaf lifespan used in marginal return calculations
    leaf_life = pars(27)
    ! load initial mean canopy NUE
    NUE_mean = pars(3)
    NUE_mean_lag = NUE_mean

    !
    ! Components to be initialised if this call is at the beginning of the model analysis.
    ! Allows potential for other MDF algorithms to be used which call the model one step at a time.
    !

    ! SHOULD TURN THIS INTO A SUBROUTINE CALL AS COMMON TO BOTH DEFAULT AND CROPS
    if (.not.allocated(deltat_1)) then

      !
      ! Allocate all variables whose dimensions are now known and invarient between iterations
      !

      ! first those linked to the time period of the analysis
      allocate(disturbance_residue_to_litter(nodays),disturbance_residue_to_cwd(nodays), &
      disturbance_residue_to_som(nodays),disturbance_loss_from_litter(nodays),  &
      disturbance_loss_from_cwd(nodays),disturbance_loss_from_som(nodays),      &
      Cwood_labile_release_coef(nodays),Croot_labile_release_coef(nodays),      &
      Cfol_turnover_coef(nodays),deltat_1(nodays),wSWP_time(nodays),            &
      soilwatermm(nodays),daylength_hours(nodays),daylength_seconds(nodays),    &
      daylength_seconds_1(nodays),meant_time(nodays),rainfall_time(nodays),     &
      airt_zero_fraction_time(nodays), &
      co2_half_saturation(nodays),co2_compensation_point(nodays), &
      air_density_kg_time(nodays),convert_ms1_mol_1_time(nodays), &
      lambda_time(nodays),psych_time(nodays),slope_time(nodays),  &
      water_vapour_diffusion_time(nodays),pn_airt_scaling_time(nodays), &
      dynamic_viscosity_time(nodays),kinematic_viscosity_time(nodays), &
      ET_demand_coef_time(nodays))

      !
      ! Timing variables which are needed first
      !

      ! inverse of time step (days-1) to avoid divisions
      deltat_1 = deltat**(-1d0)

      !
      ! Iteration independent variables using functions and thus need to be in a loop
      !

      ! then those independent of the time period, 15 years used to provide
      ! buffer for maximum value at which leaves become photosynthetically
      ! useless
      canopy_age_vector = 0d0 ; canopy_days = 0d0 ; marginal_loss_avg = 0d0
      NUE_vector = 0d0 ; NUE_vector_mature = 0d0; leaf_loss_possible = 0

      ! first those linked to the time period of the analysis
      do n = 1, nodays
        ! check positive values only for rainfall input
        rainfall_time(n) = max(0d0,met(7,n))
        ! calculate daylength in hours and seconds
        call calculate_daylength((met(6,n)-(deltat(n)*0.5d0)),lat)
        daylength_hours(n) = dayl_hours ; daylength_seconds(n) = dayl_seconds
        ! Temperature adjustments for Michaelis-Menten coefficients
        ! for CO2 (kc) and O2 (ko) and CO2 compensation point.
        co2_compensation_point(n) = arrhenious(co2comp_saturation,co2comp_half_sat_conc,met(3,n))
        co2_half_saturation(n) = arrhenious(kc_saturation,kc_half_sat_conc,met(3,n))
        pn_airt_scaling_time(n) = opt_max_scaling(pn_max_temp,pn_opt_temp,pn_kurtosis,met(3,n))
        ! calculate some temperature dependent meteorologial properties
        call meteorological_constants(met(3,n),met(3,n)+freeze)
        ! pass variables into memory objects so we don't have to keep re-calculating them
        air_density_kg_time(n) = air_density_kg
        convert_ms1_mol_1_time(n) = convert_ms1_mol_1
        lambda_time(n) = lambda
        psych_time(n) = psych
        slope_time(n) = slope
        ET_demand_coef_time(n) = ET_demand_coef
        water_vapour_diffusion_time(n) = water_vapour_diffusion
        dynamic_viscosity_time(n) = dynamic_viscosity
        kinematic_viscosity_time(n) = kinematic_viscosity
      end do

      ! calculate inverse for each time step in seconds
      daylength_seconds_1 = daylength_seconds ** (-1d0)
      ! meant time step temperature
      meant_time = (met(2,1:nodays)+met(3,1:nodays)) * 0.5d0
      ! fraction of temperture period above freezing
      airt_zero_fraction_time = (met(3,1:nodays)-0d0) / (met(3,1:nodays)-met(2,1:nodays))
      ! then those independent of the time period (8yr*365days + 90days; accounting for pre and post peak NUE )
      do n = 1, size(canopy_days)
        ! estimate cumulative age in days
        canopy_days(n) = n
      end do

      !
      ! Determine those related to phenology
      !

      ! Hydraulic limitation parameters for tissue cell expansion, i.e. growth
      ! NOTE: that these parameters are applied to deltaWP (i.e. minLWP-wSWP)
      Cwood_hydraulic_gradient = 5d0 ; Cwood_hydraulic_half_saturation = -1.5d0

      ! Temperature limitiation parameters on wood and fine root growth.
      ! Parmeters generated on the assumption of 5 % / 95 % activation at key
      ! temperature values. Roots 1oC/30oC, wood 5oC/30oC.
      ! NOTE: Foliage and root potential turnovers use the same temperature curve
      Croot_labile_release_gradient = 0.1962d0 ; Croot_labile_half_saturation = 15.0d0
      Cwood_labile_release_gradient = 0.2355d0 ; Cwood_labile_half_saturation = 17.5d0
      ! calculate temperature limitation on potential wood/root growth
      ! NOTE: could consider linear approximation between upper and lower bounds...?
      Cwood_labile_release_coef = (1d0+exp(-Cwood_labile_release_gradient* &
      (meant_time-Cwood_labile_half_saturation)))**(-1d0)
      Croot_labile_release_coef = (1d0+exp(-Croot_labile_release_gradient* &
      (meant_time-Croot_labile_half_saturation)))**(-1d0)

      !
      ! Initialise the water model
      !

      ! zero variables not done elsewhere
      water_flux = 0d0
      ! initialise some time invarient parameters
      call saxton_parameters(soil_frac_clay,soil_frac_sand)
      call initialise_soils(soil_frac_clay,soil_frac_sand)
      call update_soil_initial_conditions(pars(44))
      ! save the initial conditions for later
      soil_waterfrac_initial = soil_waterfrac
      SWP_initial = SWP
      field_capacity_initial = field_capacity
      porosity_initial = porosity

    else ! allocated or not

      !
      ! Load initial soil water conditions from memory
      !

      water_flux = 0d0
      field_capacity = field_capacity_initial
      porosity = porosity_initial

      ! input initial soil water fraction then
      ! update SWP and soil conductivity accordingly
      call update_soil_initial_conditions(pars(44))

    endif ! allocatable variables already allocated...?

    !!!!!!
    ! N cycle related parameters
    !!!!!!

    ! assign CN ratios to local variables
    CN_leaf = pars(17)/avN
    CN_root = pars(2)
    CN_wood = pars(15)
    !    CN_wood_baseline = log(pars(15))
    !    CN_wood = 10d0**(log10(pars(15)) + log10(pars(21))*pars(49))

    ! estimate time invarient N response for maintenance respiration
    Rm_leaf_baseline = Rm_reich_N(CN_leaf,pars(36),pars(37))
    Rm_root_baseline = Rm_reich_N(CN_root,pars(38),pars(39))
    Rm_wood_baseline = Rm_reich_N(CN_wood,pars(40),pars(41))

    ! reset values
    intercepted_rainfall = 0d0 ; canopy_storage = 0d0 ; snow_storage = 0d0
    root_cost = 0d0 ; leaf_cost = 0d0 ; wood_cost = 0d0

    !
    ! Initialise the disturbance model
    !

    if (maxval(met(8,1:nodays)) > 0d0 .or. maxval(met(9,1:nodays)) > 0d0) then

      ! initial values for deforestation variables
      labile_loss = 0d0    ; foliar_loss = 0d0
      roots_loss = 0d0     ; wood_loss = 0d0
      labile_residue = 0d0 ; foliar_residue = 0d0
      roots_residue = 0d0  ; wood_residue = 0d0
      stem_residue = 0d0   ; branch_residue = 0d0
      reforest_day = 0
      soil_loss_with_roots = 0d0
      coarse_root_residue = 0d0
      post_harvest_burn = 0d0

      ! now load the hardcoded forest management parameters into their locations

      ! Parameter values for deforestation variables
      ! scenario 1
      ! harvest residue (fraction); 1 = all remains, 0 = all removed
      foliage_frac_res(1) = 1d0
      roots_frac_res(1)   = 1d0
      rootcr_frac_res(1) = 1d0
      branch_frac_res(1) = 1d0
      stem_frac_res(1)   = 0d0 !
      ! wood partitioning (fraction)
      Crootcr_part(1) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
      ! Black et al 2009; Morison et al 2012)
      Cbranch_part(1) =  0.20d0 ! (Ares & Brauers 2005)
      ! actually < 15 years branches = ~25 %
      !          > 15 years branches = ~15 %.
      ! Csom loss due to phyical removal with roots
      ! Morison et al (2012) Forestry Commission Research Note
      soil_loss_frac(1) = 0.02d0 ! actually between 1-3 %
      ! was the forest burned after deforestation
      post_harvest_burn(1) = 1d0

      !## scen 2
      ! harvest residue (fraction); 1 = all remains, 0 = all removed
      foliage_frac_res(2) = 1d0
      roots_frac_res(2)   = 1d0
      rootcr_frac_res(2) = 1d0
      branch_frac_res(2) = 1d0
      stem_frac_res(2)   = 0d0 !
      ! wood partitioning (fraction)
      Crootcr_part(2) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
      ! Black et al 2009; Morison et al 2012)
      Cbranch_part(2) =  0.20d0 ! (Ares & Brauers 2005)
      ! actually < 15 years branches = ~25 %
      !          > 15 years branches = ~15 %.
      ! Csom loss due to phyical removal with roots
      ! Morison et al (2012) Forestry Commission Research Note
      soil_loss_frac(2) = 0.02d0 ! actually between 1-3 %
      ! was the forest burned after deforestation
      post_harvest_burn(2) = 0d0

      !## scen 3
      ! harvest residue (fraction); 1 = all remains, 0 = all removed
      foliage_frac_res(3) = 0.5d0
      roots_frac_res(3)   = 1d0
      rootcr_frac_res(3) = 1d0
      branch_frac_res(3) = 0d0
      stem_frac_res(3)   = 0d0 !
      ! wood partitioning (fraction)
      Crootcr_part(3) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
      ! Black et al 2009; Morison et al 2012)
      Cbranch_part(3) =  0.20d0 ! (Ares & Brauers 2005)
      ! actually < 15 years branches = ~25 %
      !          > 15 years branches = ~15 %.
      ! Csom loss due to phyical removal with roots
      ! Morison et al (2012) Forestry Commission Research Note
      soil_loss_frac(3) = 0.02d0 ! actually between 1-3 %
      ! was the forest burned after deforestation
      post_harvest_burn(3) = 0d0

      !## scen 4
      ! harvest residue (fraction); 1 = all remains, 0 = all removed
      foliage_frac_res(4) = 0.5d0
      roots_frac_res(4)   = 1d0
      rootcr_frac_res(4) = 0d0
      branch_frac_res(4) = 0d0
      stem_frac_res(4)   = 0d0
      ! wood partitioning (fraction)
      Crootcr_part(4) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
      ! Black et al 2009; Morison et al 2012)
      Cbranch_part(4) =  0.20d0 ! (Ares & Brauers 2005)
      ! actually < 15 years branches = ~25 %
      !          > 15 years branches = ~15 %.
      ! Csom loss due to phyical removal with roots
      ! Morison et al (2012) Forestry Commission Research Note
      soil_loss_frac(4) = 0.02d0 ! actually between 1-3 %
      ! was the forest burned after deforestation
      post_harvest_burn(4) = 0d0

      !## scen 5 (grassland grazing / cutting)
      ! harvest residue (fraction); 1 = all remains, 0 = all removed
      foliage_frac_res(5) = 0.1d0
      roots_frac_res(5)   = 0d0
      rootcr_frac_res(5)  = 0d0
      branch_frac_res(5)  = 0.1d0
      stem_frac_res(5)    = 0.1d0
      ! wood partitioning (fraction)
      Crootcr_part(5) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
      ! Black et al 2009; Morison et al 2012)
      Cbranch_part(5) =  0.20d0 ! (Ares & Brauers 2005)
      ! actually < 15 years branches = ~25 %
      !          > 15 years branches = ~15 %.
      ! Csom loss due to phyical removal with roots
      ! Morison et al (2012) Forestry Commission Research Note
      soil_loss_frac(5) = 0d0 ! actually between 1-3 %
      ! was the forest burned after deforestation
      post_harvest_burn(5) = 0d0

      ! for the moment override all paritioning parameters with those coming from
      ! CARDAMOM
      Cbranch_part = pars(28)
      Crootcr_part = pars(29)

      ! declare fire constants (labile, foliar, roots, wood, litter)
      combust_eff(1) = 0.1d0 ; combust_eff(2) = 0.9d0
      combust_eff(3) = 0.1d0 ; combust_eff(4) = 0.5d0
      combust_eff(5) = 0.3d0 ; rfac = 0.5d0

    end if ! disturbance ?

    !
    ! Load all variables which need to be reset between iterations
    !

    ! assigning initial conditions
    POOLS(1,1) = pars(18) ! labile
    POOLS(1,2) = pars(19) ! Foliage
    POOLS(1,3) = pars(20) ! fine roots
    POOLS(1,4) = pars(21) ! wood
    POOLS(1,5) = pars(22) ! fine litter
    POOLS(1,6) = pars(23) ! som
    POOLS(1,7) = pars(24) ! cwd
    ! POOL(1,8) assigned later

    ! load some needed module level values
    lai = POOLS(1,2)*SLA
    mint = met(2,1)  ! minimum temperature (oC)
    maxt = met(3,1)  ! maximum temperature (oC)
    swrad = met(4,1) ! incoming short wave radiation (MJ/m2/day)
    co2 = met(5,1)   ! CO2 (ppm)
    doy = met(6,1)   ! Day of year
    rainfall = rainfall_time(1) ! rainfall (kgH2O/m2/s)
    wind_spd = met(15,1) ! wind speed (m/s)
    vpd_kPa = met(16,1)*1d-3 ! vapour pressure deficit (Pa->kPa)
    meant = meant_time(1)
    leafT = maxt     ! initial canopy temperature (oC)
    maxt_lag1 = maxt
    seconds_per_step = deltat(1) * seconds_per_day
    days_per_step =  deltat(1)

    ! Calculate logistic temperature response function of leaf turnover
    !    Cfol_turnover_gradient = 0.5d0 ; Cfol_turnover_half_saturation = pars(42)
    !    Cfol_turnover_coef = (1d0+exp(Cfol_turnover_gradient*(met(2,:)-Cfol_turnover_half_saturation)))**(-1d0)

    ! initialise foliage age distribution
    canopy_age_vector = 0d0
    marginal_loss_avg = 0d0 ; marginal_gain_avg = 0d0 ! reset gain average
    leaf_loss_possible = 0

    ! NOTE: that this current initial canopy distribution model implicitly
    ! assumes a canopy turnover of < 1 year. Alternatives will be needed.

    ! Determine the youngest leaf age, based on mean age - distribution range
    youngest_leaf = max(1,nint(pars(25) - pars(45)))
    ! determine the oldest leaf are, based on mean age + distribution range
    oldest_leaf = pars(25) + pars(45)
    ! determine the vector position at which mean canopy age is located
    canopy_age = pars(25) ! is this still needed?

    ! estimate gradient needed using the integral of canopy_age-youngest_leaf (2 comes
    ! re-arranging the integral of the linear equation)
    ! to assign 50% (i.e. 0.5) of the canopy. Note that this must be scalable
    ! for simulations at different time steps
    !tmp = ((0.5d0 * 2d0) / canopy_age ** 2)
    tmp = (canopy_age - youngest_leaf) ** (-2d0)
    tmp = POOLS(1,2) * tmp
    ! assign foliar biomass to first half of the distribution
    canopy_age_vector(youngest_leaf:nint(canopy_age)) = (canopy_days(youngest_leaf:nint(canopy_age)) - dble(youngest_leaf)) * tmp

    ! now repeat the process for the second part of the distribution
    !tmp = ((0.5d0 * 2d0) / (dble(oldest_leaf) - canopy_age) ** 2)
    tmp = (dble(oldest_leaf) - canopy_age) ** (-2d0)
    tmp = POOLS(1,2) * tmp
    canopy_age_vector((nint(canopy_age)+1):oldest_leaf) = tmp * &
    (dble(oldest_leaf) - canopy_days((nint(canopy_age)+1):oldest_leaf))
    ! check / adjust mass balance
    tmp = POOLS(1,2) / sum(canopy_age_vector(1:oldest_leaf))
    canopy_age_vector(1:oldest_leaf) = canopy_age_vector(1:oldest_leaf) * tmp
    ! we have leaves therefore there must be an age
    canopy_age = sum(canopy_age_vector(1:oldest_leaf) * canopy_days(1:oldest_leaf)) &
               / POOLS(1,2) !sum(canopy_age_vector(1:oldest_leaf))

    ! estimate canopy age class specific NUE
    ! Assume linear reklationship for maturity...
    do n = 1, nint(canopy_maturation_lag)
      NUE_vector(n) = age_dependent_NUE(canopy_days(n),NUE_optimum,canopy_maturation_lag)
    end do
    ! ...then set remaining canopy to oldest_leaf as optimun
    NUE_vector(nint(canopy_maturation_lag):size(NUE_vector)) = NUE_optimum
    ! store the theoretical NUE in the abscense of environmental damanage
    NUE_vector_mature = NUE_vector

    ! If the oldest leaf is passed maturity we will assume a linear decay of NUE linked to the initial leaf lifespan
    if (oldest_leaf > nint(canopy_maturation_lag)) then
      do n = nint(canopy_maturation_lag+1d0), oldest_leaf
        NUE_vector(n) = NUE_vector(n) &
                      * (1d0 - (dble(n)-canopy_maturation_lag)/(leaf_life-canopy_maturation_lag))
      end do
      NUE_vector((oldest_leaf+1):size(NUE_vector)) = 0d0
      where (NUE_vector < 0d0) NUE_vector = 0d0
    endif

    ! reset disturbance at the beginning of iteration
    disturbance_residue_to_litter = 0d0 ; disturbance_loss_from_litter = 0d0
    disturbance_residue_to_som = 0d0 ; disturbance_loss_from_som = 0d0
    disturbance_residue_to_cwd = 0d0 ; disturbance_loss_from_cwd = 0d0

    ! Initialise root reach based on initial conditions
    root_biomass = max(min_root,POOLS(1,3)*2d0)
    root_reach = max_depth * root_biomass / (root_k + root_biomass)
    ! Determine initial soil layer thickness
    layer_thickness(1) = top_soil_depth ; layer_thickness(2) = mid_soil_depth
    layer_thickness(3) = max(min_layer,root_reach-sum(layer_thickness(1:2)))
    layer_thickness(4) = max_depth - sum(layer_thickness(1:3))
    layer_thickness(5) = top_soil_depth
    previous_depth = max(top_soil_depth,root_reach)
    ! Needed to initialise soils
    call calculate_Rtot(Rtot)
    ! Used to initialise soils
    call calculate_update_soil_water(0d0,0d0,0d0,FLUXES(1,19)) ! assume no evap or rainfall
    ! Store soil water content of the surface zone (mm)
    POOLS(1,8) = 1d3*soil_waterfrac(1)*layer_thickness(1)

    !!!!!!!!!!!!
    ! assign climate sensitivities
    !!!!!!!!!!!!

    FLUXES(1:nodays,2) = exp(pars(10)*meant_time(1:nodays))

    !
    ! Begin looping through each time step
    !

    do n = start, finish

      !!!!!!!!!!
      ! assign drivers and update some prognostic variables
      !!!!!!!!!!

      ! set lag information using previous time step value for temperature
      maxt_lag1 = maxt

      ! Incoming drivers
      mint = met(2,n)  ! minimum temperature (oC)
      maxt = met(3,n)  ! maximum temperature (oC)
      leafT = maxt     ! initial canopy temperature (oC)
      swrad = met(4,n) ! incoming short wave radiation (MJ/m2/day)
      co2 = met(5,n)   ! CO2 (ppm)
      doy = met(6,n)   ! Day of year
      rainfall = rainfall_time(n) ! rainfall (kgH2O/m2/s)
      meant = meant_time(n) ! mean air temperature (oC)
      airt_zero_fraction = airt_zero_fraction_time(n) ! fraction of temperture period above freezing
      wind_spd = met(15,n) ! wind speed (m/s)
      vpd_kPa = met(16,n)*1d-3 ! vapour pressure deficit (Pa->kPa)

      ! states needed for module variables
      lai_out(n) = POOLS(n,2)*SLA
      lai = lai_out(n) ! leaf area index (m2/m2)

      ! Temperature adjustments for Michaelis-Menten coefficients
      ! for CO2 (kc) and O2 (ko) and CO2 compensation point
      ! See McMurtrie et al., (1992) Australian Journal of Botany, vol 40, 657-677
      co2_half_sat   = co2_half_saturation(n)
      co2_comp_point = co2_compensation_point(n)
      ! temperature response for metabolically limited photosynthesis
      pn_airt_scaling = pn_airt_scaling_time(n)

      ! extract timing related values
      dayl_hours = daylength_hours(n)
      dayl_seconds = daylength_seconds(n)
      dayl_seconds_1 = daylength_seconds_1(n)
      seconds_per_step = seconds_per_day * deltat(n)
      days_per_step = deltat(n)
      days_per_step_1 = deltat_1(n)

      ! Apply todays environmentally related NUE decline
      call calculate_NUE_decline(pars(43),pars(46),pars(47),pars(48),pars(49))

      ! Estimate the current canopy NUE as a function of age
      if (sum(canopy_age_vector(1:oldest_leaf)) == 0d0 ) then
        ! if all the leaves have gone zero the age
        canopy_age = 1d0
        ! and the oldest leaf (don't forget that the earlier sections can be empty if no leaves have been allocated)
        oldest_leaf = 1
      else
        ! we have leaves therefore there must be an age
        canopy_age = sum(canopy_age_vector(1:oldest_leaf) * canopy_days(1:oldest_leaf)) &
                   / sum(canopy_age_vector(1:oldest_leaf))
      endif

      ! output canopy age for later diagnostics
      FLUXES(n,18) = canopy_age

      ! estimate the new NUE as function of canopy age
      if (POOLS(n,2) > 0d0) NUE = canopy_aggregate_NUE(1,oldest_leaf)
      ! track rolling mean of NUE
      NUE_mean_lag = ( NUE_mean_lag * (1d0 - (deltat(n) * 0.002737851d0)) ) + ( NUE * (deltat(n) * 0.002737851d0) )
!print*,"NUE",NUE
      ! snowing or not...?
      if (mint < 0d0 .and. maxt > 0d0) then
        ! if minimum temperature is below freezing point then we weight the
        ! rainfall into snow or rain based on proportion of temperature below
        ! freezing
        snowfall = rainfall * (1d0 - airt_zero_fraction) ; rainfall = rainfall - snowfall
        ! Add rainfall to the snowpack and clear rainfall variable
        snow_storage = snow_storage + (snowfall*seconds_per_step)

        ! Also melt some of the snow based on airt_zero_fraction
        ! default assumption is that snow is melting at 10 % per day light hour
        snow_melt = min(snow_storage, airt_zero_fraction * snow_storage * dayl_hours * 0.1d0 * deltat(n))
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

      !!!!!!!!!!
      ! calculate soil water potential and total hydraulic resistance
      !!!!!!!!!!

      ! calculate the minimum soil & root hydraulic resistance based on total
      ! fine root mass ! *2*2 => *RS*C->Bio
      root_biomass = max(min_root,POOLS(n,3)*2d0)
      ! estimate drythick for the current step
      drythick = max(min_drythick, top_soil_depth * min(1d0,1d0 - (soil_waterfrac(1) / porosity(1))))
      call calculate_Rtot(Rtot)
      ! Pass wSWP to output variable and update deltaWP between minlwp and
      ! current weighted soil WP
      wSWP_time(n) = wSWP ; deltaWP = min(0d0,minlwp-wSWP)

      !!!!!!!!!!
      ! Calculate surface exchange coefficients
      !!!!!!!!!!

      ! calculate some temperature dependent meteorologial properties
      !call meteorological_constants(maxt,maxt+freeze)
      ! pass variables from memory objects
      air_density_kg = air_density_kg_time(n)
      convert_ms1_mol_1 = convert_ms1_mol_1_time(n)
      lambda = lambda_time(n) ; psych = psych_time(n)
      slope = slope_time(n) ; ET_demand_coef = ET_demand_coef_time(n)
      water_vapour_diffusion = water_vapour_diffusion_time(n)
      dynamic_viscosity = dynamic_viscosity_time(n)
      kinematic_viscosity = kinematic_viscosity_time(n)
      ! calculate aerodynamic using consistent approach with SPA
      call calculate_aerodynamic_conductance

      !!!!!!!!!!
      ! Determine net shortwave and isothermal longwave energy balance
      !!!!!!!!!!

      call calculate_shortwave_balance
      call calculate_longwave_isothermal(leafT,maxt)

      !!!!!!!!!!
      ! Estimate evaporative and photosynthetic fluxes
      !!!!!!!!!!

      ! Canopy intercepted rainfall evaporation (kgH2O/m2/day)
      call calculate_wetcanopy_evaporation(wetcanopy_evap,act_pot_ratio,canopy_storage,0d0)

      ! calculate radiation absorption and estimate stomatal conductance
      call acm_albedo_gc(abs(deltaWP),Rtot)

      ! if snow present assume that soilevaporation is sublimation of soil first
      if (snow_storage > 0d0) then
        snow_sublimation = soilevaporation
        if (snow_sublimation*deltat(n) > snow_storage) snow_sublimation = snow_storage * deltat_1(n)
        soilevaporation = soilevaporation - snow_sublimation
        snow_storage = snow_storage - (snow_sublimation * deltat(n))
      else
        snow_sublimation = 0d0
      end if

      ! Note that soil mass balance will be calculated after phenology
      ! adjustments

      ! reset output variable
      if (lai > vsmall .and. stomatal_conductance > vsmall) then
        ! Gross primary productivity (gC/m2/day)
        FLUXES(n,1) = max(0d0,acm_gpp(stomatal_conductance))
        ! Canopy transpiration (kgH2O/m2/day)
        call calculate_transpiration(transpiration)
        ! restrict transpiration to positive only
        transpiration = max(0d0,transpiration)
      else
        ! assume zero fluxes
        FLUXES(n,1) = 0d0 ; transpiration = 0d0 ; ci_global = 0d0
      endif
      ! labile production (gC.m-2.day-1)
      FLUXES(n,5) = FLUXES(n,1)

      ! Soil surface (kgH2O.m-2.day-1)
      call calculate_soil_evaporation(soilevaporation)

      !!!!!!!!!!
      ! calculate maintenance respiration demands and mass balance
      !!!!!!!!!!

      ! Autotrophic maintenance respiration demand (nmolC.g-1.s-1 -> gC.m-2.day-1)
      ! NOTE: 2d-3 scales between nmolC->umolC and from biomass to carbon
      Q10_adjustment = Rm_reich_Q10(meant)
      tmp = umol_to_gC*seconds_per_day*2d-3
      ! now apply pool size an
      Rm_leaf = Rm_leaf_baseline * Q10_adjustment * tmp * POOLS(n,2)
      Rm_root = Rm_root_baseline * Q10_adjustment * tmp * POOLS(n,3)
      Rm_wood = Rm_wood_baseline * Q10_adjustment * tmp * POOLS(n,4)

      ! reset overall value as this is used in flag later
      Rm_deficit = 0d0
      ! determine if there is greater demand for Rm than available labile C
      avail_labile = POOLS(n,1)*deltat_1(n)
      if ( (Rm_leaf + Rm_root + Rm_wood) > avail_labile ) then

        ! reset tissue specific values
        Rm_deficit_leaf_loss = 0d0 ; Rm_deficit_root_loss = 0d0 ; Rm_deficit_wood_loss = 0d0

        ! More Rm demanded than available labile, therefore mortality will
        ! occur. Mortality is apportioned based on fraction of Rm attributed
        ! to each live biomass pool. The total biomass loss needed to make up
        ! the deficit in the current time step is committed to death.
        ! NOTE: this could be all of the plant!

        ! Assign all labile to the current Ra flux output
        FLUXES(n,3) = avail_labile
        ! borrow the avail_labile variable, remember to zero at end of this
        ! section
        avail_labile = Rm_leaf + Rm_root + Rm_wood
        ! then determine the overshoot
        Rm_deficit = avail_labile - FLUXES(n,3)
        ! to avoid further division
        avail_labile = avail_labile ** (-1d0)
        ! calculate proportion of loss to leaves (gC.m-2.day-1)
        if (Rm_leaf > 0d0) then
          Rm_deficit_leaf_loss = Rm_deficit * (Rm_leaf * avail_labile)
          Rm_deficit_leaf_loss = Rm_deficit_leaf_loss / (Rm_leaf / POOLS(n,2))
          ! adjust each to fractional equivalents to allow for time step adjustment
          Rm_deficit_leaf_loss = min(1d0,Rm_deficit_leaf_loss / POOLS(n,2))
          Rm_deficit_leaf_loss = POOLS(n,2) &
                               *(1d0-(1d0-Rm_deficit_leaf_loss)**deltat(n))*deltat_1(n)
        endif
        if (Rm_root > 0d0) then
          ! calculate proportion of loss to roots (gC.m-2.day-1)
          Rm_deficit_root_loss = Rm_deficit * (Rm_root * avail_labile)
          Rm_deficit_root_loss = Rm_deficit_root_loss / (Rm_root / POOLS(n,3))
          ! adjust each to fractional equivalents to allow for time step adjustment
          Rm_deficit_root_loss = min(1d0,Rm_deficit_root_loss / POOLS(n,3))
          Rm_deficit_root_loss = POOLS(n,3) &
                               *(1d0-(1d0-Rm_deficit_root_loss)**deltat(n))*deltat_1(n)
        endif
        if (Rm_wood > 0d0) then
          ! calculate proportion of loss to wood (gC.m-2.day-1)
          Rm_deficit_wood_loss = Rm_deficit * (Rm_wood * avail_labile)
          Rm_deficit_wood_loss = Rm_deficit_wood_loss / (Rm_wood / POOLS(n,4))
          ! adjust each to fractional equivalents to allow for time step adjustment
          Rm_deficit_wood_loss = min(1d0,Rm_deficit_wood_loss / POOLS(n,4))
          Rm_deficit_wood_loss = POOLS(n,4) &
                               *(1d0-(1d0-Rm_deficit_wood_loss)**deltat(n))*deltat_1(n)
        endif
        ! reset available labile to zero
        avail_labile = 0d0
      else
        ! we have enough labile so assign the demand
        FLUXES(n,3) = Rm_leaf + Rm_root + Rm_wood
        ! then update the available labile supply for growth
        avail_labile = POOLS(n,1) - (FLUXES(n,3)*deltat(n))
      endif

      !!!!!!!!!!
      ! calculate canopy phenology
      !!!!!!!!!!

      ! calculate hydraulic limits on leaf / wood growth.
      ! NOTE: PARAMETERS SHOULD BE CALIBRATRED OR TISSUE SPECIFIC
      Cwood_hydraulic_limit = (1d0+exp(Cwood_hydraulic_gradient*(deltaWP-Cwood_hydraulic_half_saturation)))**(-1d0)

      ! Determine leaf growth and turnover based on marginal return calculations
      ! NOTE: that turnovers will be bypassed in favour of mortality turnover
      ! should available labile be exhausted
      call calculate_leaf_dynamics(n,deltat,nodays,pars(12),pars(36),pars(37) &
                                  ,deltaWP,Rtot,FLUXES(n,1),Rm_leaf &
                                  ,POOLS(n,2),FLUXES(n,10),FLUXES(n,8))

      ! Update available labile supply for fine roots and wood
      avail_labile = avail_labile - (FLUXES(n,8)*deltat(n))

      !!!!!!!!!!
      ! calculate wood and root phenology
      !!!!!!!!!!

      ! calculate allocation of labile to roots and wood including, where appropriate, marginal return calculations
      call calculate_wood_root_growth(n,pars(4),pars(13),deltaWP,Rtot,FLUXES(n,1) &
                                     ,POOLS(n,3),POOLS(n,4),FLUXES(n,6),FLUXES(n,7))

      !!!!!!!!!!
      ! litter creation with time dependancies
      !!!!!!!!!!

      if (Rm_deficit > 0d0) then

        ! C starvation turnover has occured, mortality turnover will now occur

        ! we have some work to do if the canopy needs to lose more leaf than expected...
        if (FLUXES(n,10) < Rm_deficit_leaf_loss) then

          ! determine the deficit between what is already committed to loss and what is needed
          tmp = Rm_deficit_leaf_loss - FLUXES(n,10)
          ! scale the the full model time step
          tmp = tmp * deltat(n)
          ! adjust to maximum of the pool size
          tmp = min(POOLS(n,2),tmp)
          ! increment age profiles, remove old leaf...
          loss = 0d0 ; b = oldest_leaf ; c = oldest_leaf
          do a = oldest_leaf, 1, -1
            ! keep track of the position within the canopy for use outside
            ! of the loop
            b = a
            if (canopy_age_vector(b) > 0d0) then
              ! how much canopy can we lose for the current age classes
              loss = loss + canopy_age_vector(b)
              ! if we have found enough carbon to remove then exit the loop
              if (loss > tmp) exit
              ! each time we pass through but it is not enough we need to
              ! check the location of the available C in age classes
              c = a
            endif
          end do
          ! ensure that the looping information within is now passed out-with
          a = b

          ! remove the biomass from the canopy age vector
          canopy_age_vector(a:oldest_leaf) = 0d0
          ! and marginal return calculation
          marginal_loss_avg(a:oldest_leaf) = 0d0
          ! and counter
          leaf_loss_possible(a:oldest_leaf) = 0
          ! update the new oldest leaf age
          oldest_leaf = c
          ! cannot have oldest leaf at position less than 1
          if (oldest_leaf <= 0) oldest_leaf = 1

          ! update FLX10 for leaf loss
          FLUXES(n,10) = FLUXES(n,10) + (loss*deltat_1(n))

        endif ! FLX10 < Rm_deficit_leaf_loss

        ! update root & wood losses based on C starvation
        FLUXES(n,11) = Rm_deficit_wood_loss
        FLUXES(n,12) = Rm_deficit_root_loss

      else

        ! C starvation turnover not occuring so turnovers progress as normal

        ! total wood litter production
        FLUXES(n,11) = POOLS(n,4)*(1d0-(1d0-pars(6))**deltat(n))*deltat_1(n)
        ! total root litter production
        FLUXES(n,12) = POOLS(n,3)*(1d0-(1d0-pars(7))**deltat(n))*deltat_1(n)

      endif

      ! if 12 months has gone by, update the leaf lifespan variable
      if (n > steps_per_year .and. met(6,n) < met(6,n-1)) then
        ! determine the turnover fraction across the year
        tmp = sum(FLUXES((n-steps_per_year):(n-1),10) + FLUXES((n-steps_per_year):(n-1),23)) &
            / sum(POOLS((n-steps_per_year):(n-1),2))
        ! i.e. we cannot / should not update the leaf lifespan if there has
        ! been no turnover and / or there is no foliar pool.
        if (tmp > 0d0) then
          tmp = tmp ** (-1d0)
          tmp = tmp * leaf_life_weighting
          leaf_life = tmp + (leaf_life * (1d0 - leaf_life_weighting))
          ! if we have updated the leaf_life we should also update the canopy NUE_mean
          !call estimate_mean_NUE
          NUE_mean = NUE_mean_lag
        end if
      endif ! n /= 1 and new calendar year

      !!!!!!!!!!
      ! those with temperature AND time dependancies
      !!!!!!!!!!

      ! respiration heterotrophic litter
      FLUXES(n,13) = POOLS(n,5)*(1d0-(1d0-FLUXES(n,2)*pars(8))**deltat(n))*deltat_1(n)
      ! respiration heterotrophic som
      FLUXES(n,14) = POOLS(n,6)*(1d0-(1d0-FLUXES(n,2)*pars(9))**deltat(n))*deltat_1(n)
      ! litter to som
      FLUXES(n,15) = POOLS(n,5)*(1d0-(1d0-FLUXES(n,2)*pars(1))**deltat(n))*deltat_1(n)
      ! CWD to litter
      FLUXES(n,20) = POOLS(n,7)*(1d0-(1d0-FLUXES(n,2)*pars(16))**deltat(n))*deltat_1(n)

      !!!!!!!!!!
      ! calculate growth respiration and adjust allocation to pools assuming
      ! 0.21875 of total C allocation towards each pool (i.e. 0.28 .eq. xNPP)
      !!!!!!!!!!

      ! foliage
      Rg_from_labile =                   FLUXES(n,8)*Rg_fraction  ; FLUXES(n,8) = FLUXES(n,8) * one_Rg_fraction
      ! roots
      Rg_from_labile = Rg_from_labile + (FLUXES(n,6)*Rg_fraction) ; FLUXES(n,6) = FLUXES(n,6) * one_Rg_fraction
      ! wood
      Rg_from_labile = Rg_from_labile + (FLUXES(n,7)*Rg_fraction) ; FLUXES(n,7) = FLUXES(n,7) * one_Rg_fraction
      ! now update the Ra flux with Rg
      FLUXES(n,3) = FLUXES(n,3) + Rg_from_labile

      !!!!!!!!!!
      ! update pools for next timestep
      !!!!!!!!!!

      ! labile pool
      POOLS(n+1,1) = POOLS(n,1) + (FLUXES(n,5)-FLUXES(n,8)-FLUXES(n,6)-FLUXES(n,7)-FLUXES(n,3))*deltat(n)
      ! foliar pool
      ! POOLS(n+1,2) = POOLS(n,2) + (FLUXES(n,8)-FLUXES(n,10))*deltat(n)
      POOLS(n+1,2) = sum(canopy_age_vector(1:oldest_leaf))
      ! wood pool
      POOLS(n+1,4) = POOLS(n,4) + (FLUXES(n,7)-FLUXES(n,11))*deltat(n)
      ! root pool
      POOLS(n+1,3) = POOLS(n,3) + (FLUXES(n,6)-FLUXES(n,12))*deltat(n)
      ! litter pool
      POOLS(n+1,5) = POOLS(n,5) + (FLUXES(n,10)+FLUXES(n,12)+FLUXES(n,20)-FLUXES(n,13)-FLUXES(n,15))*deltat(n)
      ! som pool
      POOLS(n+1,6) = POOLS(n,6) + (FLUXES(n,15)-FLUXES(n,14))*deltat(n)
      ! cwd pool
      POOLS(n+1,7) = POOLS(n,7) + (FLUXES(n,11)-FLUXES(n,20))*deltat(n)

      !!!!!!!!!!
      ! Update soil water balance
      !!!!!!!!!!

      ! add any snow melt to the rainfall now that we have already dealt with the canopy interception
      rainfall = rainfall + (snow_melt / seconds_per_step)
      ! do mass balance (i.e. is there enough water to support ET)
      call calculate_update_soil_water(transpiration,soilevaporation,((rainfall-intercepted_rainfall)*seconds_per_day) &
                                      ,FLUXES(n,19))
      ! now that soil mass balance has been updated we can add the wet canopy
      ! evaporation (kgH2O.m-2.day-1)
      FLUXES(n,19) = FLUXES(n,19) + wetcanopy_evap
      ! store soil water content of the surface zone (mm)
      POOLS(n+1,8) = 1d3 * soil_waterfrac(1) * layer_thickness(1)

      !!!!!!!!!!
      ! deal first with deforestation
      !!!!!!!!!!

      if (n == reforest_day) then
        POOLS(n+1,1) = pars(30)
        POOLS(n+1,2) = pars(31)
        POOLS(n+1,3) = pars(32)
        POOLS(n+1,4) = pars(33)
      endif

      ! reset values
      harvest_management = 0 ; burnt_area = 0d0

      if (met(8,n) > 0d0) then

        ! pass harvest management to local integer
        harvest_management = int(met(13,n))

        ! assume that labile is proportionally distributed through the plant
        ! root and wood and therefore so is the residual fraction
        C_total = POOLS(n+1,3) + POOLS(n+1,4)
        ! partition wood into its components
        Cbranch = POOLS(n+1,4)*Cbranch_part(harvest_management)
        Crootcr = POOLS(n+1,4)*Crootcr_part(harvest_management)
        Cstem   = POOLS(n+1,4)-(Cbranch + Crootcr)
        ! now calculate the labile fraction of residue
        if (C_total > 0d0) then
          labile_frac_res = ((POOLS(n+1,3)/C_total) * roots_frac_res(harvest_management)  ) &
                          + ((Cbranch/C_total)      * branch_frac_res(harvest_management) ) &
                          + ((Cstem/C_total)        * stem_frac_res(harvest_management)   ) &
                          + ((Crootcr/C_total)      * rootcr_frac_res(harvest_management) )
        else
          labile_frac_res = 0d0
        endif

        ! you can't remove any biomass if there is none left...
        if (C_total > vsmall) then

          ! Loss of carbon from each pools
          labile_loss = POOLS(n+1,1)*met(8,n)
          foliar_loss = POOLS(n+1,2)*met(8,n)
          ! roots are not removed under grazing
          if (harvest_management /= 5) then
            roots_loss = POOLS(n+1,3)*met(8,n)
          else
            roots_loss = 0d0
          endif
          wood_loss   = (Cbranch+Crootcr+Cstem)*met(8,n)
          ! estimate labile loss explicitly from the loss of their storage
          ! tissues
          labile_loss = POOLS(n+1,1) * ((roots_loss+wood_loss) / (POOLS(n+1,3)+POOLS(n+1,4)))

          ! For output / EDC updates
          if (met(8,n) <= 0.99d0) then
            FLUXES(n,22) = labile_loss * deltat_1(n)
            FLUXES(n,23) = foliar_loss * deltat_1(n)
            FLUXES(n,24) = roots_loss * deltat_1(n)
            FLUXES(n,25) = wood_loss * deltat_1(n)
          endif
          ! Transfer fraction of harvest waste to litter or som pools
          ! easy pools first
          labile_residue = labile_loss*labile_frac_res
          foliar_residue = foliar_loss*foliage_frac_res(harvest_management)
          roots_residue  = roots_loss*roots_frac_res(harvest_management)
          ! Explicit calculation of the residues from each fraction
          coarse_root_residue  = Crootcr*met(8,n)*rootcr_frac_res(harvest_management)
          branch_residue = Cbranch*met(8,n)*branch_frac_res(harvest_management)
          stem_residue = Cstem*met(8,n)*stem_frac_res(harvest_management)
          ! Now finally calculate the final wood residue
          wood_residue = stem_residue + branch_residue + coarse_root_residue
          ! Mechanical loss of Csom due to coarse root extraction
          soil_loss_with_roots = Crootcr*met(8,n)*(1d0-rootcr_frac_res(harvest_management)) &
                               * soil_loss_frac(harvest_management)

          ! Update pools
          POOLS(n+1,1) = max(0d0, POOLS(n+1,1)-labile_loss)
          ! POOLS(n+1,2) = max(0d0, POOLS(n+1,2)-foliar_loss)
          POOLS(n+1,3) = max(0d0, POOLS(n+1,3)-roots_loss)
          POOLS(n+1,4) = max(0d0, POOLS(n+1,4)-wood_loss)
          POOLS(n+1,5) = max(0d0, POOLS(n+1,5) + (labile_residue+foliar_residue+roots_residue))
          POOLS(n+1,6) = max(0d0, POOLS(n+1,6) - soil_loss_with_roots)
          POOLS(n+1,7) = max(0d0, POOLS(n+1,7) + wood_residue)

          ! update corresponding canopy vector
          if (sum(canopy_age_vector(1:oldest_leaf)) > 0d0) then
             canopy_age_vector(1:oldest_leaf) = canopy_age_vector(1:oldest_leaf) &
                                              - ( (canopy_age_vector(1:oldest_leaf) / sum(canopy_age_vector(1:oldest_leaf)))&
                                              * foliar_loss)
            ! correct for precision errors
            where (canopy_age_vector(1:oldest_leaf) < 0d0)
              canopy_age_vector(1:oldest_leaf) = 0d0
              marginal_loss_avg(1:oldest_leaf) = 0d0
              leaf_loss_possible(1:oldest_leaf) = 0
            end where
            POOLS(n+1,2) = sum(canopy_age_vector(1:oldest_leaf))
          endif

          ! Some variable needed for the EDCs
          ! reallocation fluxes for the residues
          disturbance_residue_to_litter(n) = labile_residue+foliar_residue+roots_residue
          disturbance_loss_from_litter(n)  = 0d0
          disturbance_residue_to_cwd(n)    = wood_residue
          disturbance_loss_from_cwd(n)     = 0d0
          disturbance_residue_to_som(n)    = 0d0
          disturbance_loss_from_som(n)     = soil_loss_with_roots
          ! Convert all to rates to be consistent with the FLUXES in EDCs
          disturbance_residue_to_litter(n) = disturbance_residue_to_litter(n) * deltat_1(n)
          disturbance_loss_from_litter(n)  = disturbance_loss_from_litter(n) * deltat_1(n)
          disturbance_residue_to_cwd(n)    = disturbance_residue_to_cwd(n) * deltat_1(n)
          disturbance_loss_from_cwd(n)     = disturbance_loss_from_cwd(n) * deltat_1(n)
          disturbance_residue_to_som(n)    = disturbance_residue_to_som(n) * deltat_1(n)
          disturbance_loss_from_som(n)     = disturbance_loss_from_som(n) * deltat_1(n)
          ! estimate total C extraction
          ! NOTE: this calculation format is to prevent precision error in calculation
          FLUXES(n,21) = wood_loss + labile_loss + foliar_loss + roots_loss
          FLUXES(n,21) = FLUXES(n,21) - (wood_residue + labile_residue + foliar_residue + roots_residue)
          ! Convert to daily rate
          FLUXES(n,21) = FLUXES(n,21) * deltat_1(n)

        end if ! C_total > vsmall

        ! Total carbon loss from the system
        C_total = (labile_residue+foliar_residue+roots_residue+wood_residue+sum(NCFF)) &
                - (labile_loss+foliar_loss+roots_loss+wood_loss+soil_loss_with_roots+sum(CFF))

        ! If total clearance occured then we need to ensure some minimum
        ! values and reforestation is assumed one year forward
        if (met(8,n) > 0.99d0) then
          m = 0 ; test = nint(sum(deltat(n:(n+m))))
          ! FC Forest Statistics 2015 lag between harvest and restocking ~ 2 year
          restocking_lag = 365*2
          do while (test < restocking_lag)
            m = m + 1 ; test = nint(sum(deltat(n:(n+m))))
            !  get out clause for hitting the end of the simulation
            if (m+n >= nodays) test = restocking_lag
          enddo
          reforest_day = min((n+m), nodays)
        endif ! if total clearance

      endif ! end deforestation info

      !!!!!!!!!!
      ! then deal with fire
      !!!!!!!!!!

      if (met(9,n) > 0d0 .or.(met(8,n) > 0d0 .and. harvest_management > 0)) then

        burnt_area = met(9,n)
        if (met(8,n) > 0d0 .and. burnt_area > 0d0) then
          ! pass harvest management to local integer
          burnt_area = min(1d0,burnt_area + post_harvest_burn(harvest_management))
        else if (met(8,n) > 0d0 .and. burnt_area <= 0d0) then
          burnt_area = post_harvest_burn(harvest_management)
        endif

        if (burnt_area > 0d0) then

          !/*first fluxes*/
          !/*LABILE*/
          CFF(1) = POOLS(n+1,1)*burnt_area*combust_eff(1)
          NCFF(1) = POOLS(n+1,1)*burnt_area*(1d0-combust_eff(1))*(1d0-rfac)
          !/*foliar*/
          CFF(2) = POOLS(n+1,2)*burnt_area*combust_eff(2)
          NCFF(2) = POOLS(n+1,2)*burnt_area*(1d0-combust_eff(2))*(1d0-rfac)
          !/*root*/
          CFF(3) = 0d0 !POOLS(n+1,3)*burnt_area*combust_eff(3)
          NCFF(3) = 0d0 !POOLS(n+1,3)*burnt_area*(1d0-combust_eff(3))*(1d0-rfac)
          !/*wood*/
          CFF(4) = POOLS(n+1,4)*burnt_area*combust_eff(4)
          NCFF(4) = POOLS(n+1,4)*burnt_area*(1d0-combust_eff(4))*(1d0-rfac)
          !/*litter*/
          CFF(5) = POOLS(n+1,5)*burnt_area*combust_eff(5)
          NCFF(5) = POOLS(n+1,5)*burnt_area*(1d0-combust_eff(5))*(1d0-rfac)
          ! CWD; assume same as live wood (should be improved later)
          CFF(7) = POOLS(n+1,7)*burnt_area*combust_eff(4)
          NCFF(7) = POOLS(n+1,7)*burnt_area*(1d0-combust_eff(4))*(1d0-rfac)
          !/*fires as daily averages to comply with units*/
          FLUXES(n,17) = (CFF(1)+CFF(2)+CFF(3)+CFF(4)+CFF(5)) * deltat_1(n)
          !              !/*update net exchangep*/
          !              NEE(n)=NEE(n)+FLUXES(n,17)
          ! determine the as daily rate impact on live tissues for use in EDC and
          ! MTT calculations
          FLUXES(n,22) = FLUXES(n,22) + ((CFF(1) + NCFF(1)) * deltat_1(n)) ! labile
          FLUXES(n,23) = FLUXES(n,23) + ((CFF(2) + NCFF(2)) * deltat_1(n)) ! foliar
          FLUXES(n,24) = FLUXES(n,24) + ((CFF(3) + NCFF(3)) * deltat_1(n)) ! root
          FLUXES(n,25) = FLUXES(n,25) + ((CFF(4) + NCFF(4)) * deltat_1(n)) ! wood

          !// update pools
          !/*Adding all fire pool transfers here*/
          POOLS(n+1,1) = max(0d0,POOLS(n+1,1)-CFF(1)-NCFF(1))
          ! POOLS(n+1,2) = max(0d0,POOLS(n+1,2)-CFF(2)-NCFF(2))
          POOLS(n+1,3) = max(0d0,POOLS(n+1,3)-CFF(3)-NCFF(3))
          POOLS(n+1,4) = max(0d0,POOLS(n+1,4)-CFF(4)-NCFF(4))
          POOLS(n+1,5) = max(0d0,POOLS(n+1,5)-CFF(5)-NCFF(5)+NCFF(1)+NCFF(2)+NCFF(3))
          POOLS(n+1,6) = max(0d0,POOLS(n+1,6)+NCFF(4)+NCFF(5)+NCFF(7))
          POOLS(n+1,7) = max(0d0,POOLS(n+1,7)-CFF(7)-NCFF(7))

          ! update corresponding canopy vector
          if (sum(canopy_age_vector(1:oldest_leaf)) > 0d0) then
            canopy_age_vector(1:oldest_leaf) = canopy_age_vector(1:oldest_leaf) &
                                             - ( (canopy_age_vector(1:oldest_leaf) / sum(canopy_age_vector(1:oldest_leaf))) &
                                               * (CFF(2)+NCFF(2)))
            ! correct for precision errors
            where (canopy_age_vector(1:oldest_leaf) < 0d0)
              canopy_age_vector(1:oldest_leaf) = 0d0
              marginal_loss_avg(1:oldest_leaf) = 0d0
              leaf_loss_possible(1:oldest_leaf) = 0
            end where
            POOLS(n+1,2) = sum(canopy_age_vector(1:oldest_leaf))
          endif

          ! Some variable needed for the EDCs
          ! Reallocation fluxes for the residues, remember to
          ! convert to daily rate for consistency with the EDCs
          ! NOTE: accumulation because fire and removal may occur concurrently...
          disturbance_residue_to_litter(n) = disturbance_residue_to_litter(n) &
                                           + ((NCFF(1)+NCFF(2)+NCFF(3)) * deltat_1(n))
          disturbance_residue_to_som(n)    = disturbance_residue_to_som(n) &
                                           + ((NCFF(4)+NCFF(5)+NCFF(7)) * deltat_1(n))
          disturbance_loss_from_litter(n)  = disturbance_loss_from_litter(n) &
                                           + ((CFF(5) + NCFF(5)) * deltat_1(n))
          disturbance_loss_from_cwd(n)     = disturbance_loss_from_cwd(n) &
                                           + ((CFF(7) - NCFF(7)) * deltat_1(n))

        endif ! burn area > 0

      endif ! fire activity

      !!!!!!!!!
      ! Bug checking
      !!!!!!!!!

      ! Check that foliar pool in bulk and age specific are balanced
      if (canopy_age /= canopy_age .or. sum(canopy_age_vector(1:oldest_leaf)) - POOLS(n+1,2) > 1d-8) then

          ! Estimate the current canopy NUE as a function of age
          if (sum(canopy_age_vector(1:oldest_leaf)) == 0d0 ) then
              ! if all the leaves have gone zero the age
              canopy_age = 1d0
              ! and the oldest leaf (don't forget that the earlier sections can be
              ! empty if no leaves have been allocated)
              oldest_leaf = 1
          else
              ! we have leaves therefore there must be an age
              canopy_age = sum(canopy_age_vector(1:oldest_leaf) * canopy_days(1:oldest_leaf)) &
                         / sum(canopy_age_vector(1:oldest_leaf))
          endif

          print*,"Mass balance error or canopy_age == NaN"
          print*,"step",n
          print*,"met",met(:,n)
          print*,"POOLS",POOLS(n,:)
          print*,"FLUXES",FLUXES(n,:)
          print*,"POOLS+1",POOLS(n+1,:)
          print*,"wSWP",wSWP,"stomatal_conductance",stomatal_conductance
          print*,"WetCan_evap",wetcanopy_evap,"transpiration",transpiration,"soilevap",soilevaporation
          print*,"waterfrac",soil_waterfrac
          print*,"Rm_loss",Rm_deficit_leaf_loss,Rm_deficit_root_loss,Rm_deficit_wood_loss
          print*,"Canopy_age",canopy_age,"Oldest_leaf",oldest_leaf,"oldest_age",canopy_days(oldest_leaf)
          print*,"Canopy_age_vector",sum(canopy_age_vector(1:oldest_leaf))
          print*,"Vector-Bulk Foliar Residual",sum(canopy_age_vector(1:oldest_leaf))-POOLS(n+1,2)
          stop

      endif ! canopy_age /= canopy_age

      if (sum(POOLS(n+1,1:nopools)) /= sum(POOLS(n+1,1:nopools)) .or. &
          sum(POOLS(n,1:nopools)) /= sum(POOLS(n,1:nopools)) .or. &
          minval(POOLS(n+1,1:nopools)) < 0d0) then
          ! if there is a problem search for a more specific problem
          do nxp = 1, nopools
             if (POOLS(n+1,nxp) /= POOLS(n+1,nxp) .or. POOLS(n+1,nxp) < 0d0 .or. &
                 POOLS(n,nxp) /= POOLS(n,nxp) .or. POOLS(n,nxp) < 0d0) then
                 print*,"POOLS related error"
                 print*,"step",n,"POOL",nxp
                 print*,"met",met(:,n)
                 print*,"POOLS",POOLS(n,:)
                 print*,"FLUXES",FLUXES(n,:)
                 print*,"POOLS+1",POOLS(n+1,:)
                 print*,"wSWP",wSWP,"stomatal_conductance",stomatal_conductance
                 print*,"WetCan_evap",wetcanopy_evap,"transpiration",transpiration,"soilevap",soilevaporation
                 print*,"waterfrac",soil_waterfrac
                 print*,"Rm_loss",Rm_deficit_leaf_loss,Rm_deficit_root_loss,Rm_deficit_wood_loss
                 print*,"Canopy_age",canopy_age,"Oldest_leaf",oldest_leaf
                 stop
             endif
          enddo

      endif ! vectorised check for NaN or negatives

      if (sum(FLUXES(n,1:nofluxes)) /= sum(FLUXES(n,1:nofluxes))) then
        ! if there is a problem search for more specific error information
        do nxp = 1, nofluxes
          if (FLUXES(n,nxp) /= FLUXES(n,nxp)) then
            print*,"FLUXES related error"
            print*,"step",n,"FLUX",nxp
            print*,"met",met(:,n)
            print*,"POOLS",POOLS(n,:)
            print*,"FLUXES",FLUXES(n,:)
            print*,"POOLS+1",POOLS(n+1,:)
            print*,"wSWP",wSWP,"stomatal_conductance",stomatal_conductance
            print*,"WetCan_evap",wetcanopy_evap,"transpiration",transpiration,"soilevap",soilevaporation
            print*,"waterfrac",soil_waterfrac
            print*,"Rm_loss",Rm_deficit_leaf_loss,Rm_deficit_root_loss,Rm_deficit_wood_loss
            print*,"Canopy_age",canopy_age,"Oldest_leaf",oldest_leaf
            stop
          endif
          ! do not include evaporation this assumption
          if (FLUXES(n,nxp) < 0d0 .and. nxp /= 19) then
            print*,"FLUXES related error"
            print*,"step",n,"FLUX",nxp
            print*,"met",met(:,n)
            print*,"POOLS",POOLS(n,:)
            print*,"FLUXES",FLUXES(n,:)
            print*,"POOLS+1",POOLS(n+1,:)
            print*,"wSWP",wSWP,"stomatal_conductance",stomatal_conductance
            print*,"WetCan_evap",wetcanopy_evap,"transpiration",transpiration,"soilevap",soilevaporation
            print*,"waterfrac",soil_waterfrac
            print*,"Rm_loss",Rm_deficit_leaf_loss,Rm_deficit_root_loss,Rm_deficit_wood_loss
            print*,"Canopy_age",canopy_age,"Oldest_leaf",oldest_leaf
            stop
          endif
        enddo

      end if ! vectorised check for NaN or negatives

    end do ! nodays loop

    !!!!!!!!!!
    ! Calculate Ecosystem diagnostics
    !!!!!!!!!!

    ! calculate NEE
    NEE_out(1:nodays) = -FLUXES(1:nodays,1) & ! GPP
                      +FLUXES(1:nodays,3)+FLUXES(1:nodays,13)+FLUXES(1:nodays,14) !& ! Respiration
                      !+FLUXES(1:nodays,17)  ! fire

    ! load GPP
    GPP_out(1:nodays) = FLUXES(1:nodays,1)

    !call cpu_time(done) ; print*,"Total",done-begin
  end subroutine CARBON_MODEL
  !
  !------------------------------------------------------------------
  !
  double precision function acm_gpp(gs)

    ! the Aggregated Canopy Model, is a Gross Primary Productivity (i.e.
    ! Photosyntheis) emulator which operates at a daily time step. ACM can be
    ! paramaterised to provide reasonable results for most ecosystems.

    implicit none

    ! declare input variables
    double precision, intent(in) :: gs

    ! declare local variables
    double precision :: pn, pd, pp, qq, ci, mult, pl &
                       ,gc ,gs_mol, gb_mol


    ! Temperature adjustments for Michaelis-Menten coefficients
    ! for CO2 (kc) and O2 (ko) and CO2 compensation point
    ! See McMurtrie et al., (1992) Australian Journal of Botany, vol 40, 657-677
!    co2_half_sat   = arrhenious(kc_saturation,kc_half_sat_conc,leafT)
!    co2_comp_point = arrhenious(co2comp_saturation,co2comp_half_sat_conc,leafT)


    !
    ! Metabolic limited photosynthesis
    !

    ! maximum rate of temperature and nitrogen (canopy efficiency) limited
    ! photosynthesis (gC.m-2.day-1)
    !pn = lai*avN*NUE*opt_max_scaling(pn_max_temp,pn_opt_temp,pn_kurtosis,leafT)
    pn = lai*avN*NUE*pn_airt_scaling

    !
    ! Diffusion limited photosynthesis
    !

    ! daily canopy conductance (mmolH2O.m-2.s-1-> molCO2.m-2.day-1)
    ! The ratio of H20:CO2 diffusion is 1.646259 (Jones appendix 2).
    ! i.e. gcH2O*1.646259 = gcCO2
    gs_mol = gs * seconds_per_day * gs_H2Ommol_CO2mol
    ! canopy level boundary layer conductance unit change
    ! (m.s-1 -> mol.m-2.day-1) assuming sea surface pressure only.
    ! Note the ratio of H20:CO2 diffusion through leaf level boundary layer is
    ! 1.37 (Jones appendix 2).
    gb_mol = aerodynamic_conductance * seconds_per_day * convert_ms1_mol_1 * gb_H2O_CO2
    ! Combining in series the stomatal and boundary layer conductances
    gc = (gs_mol ** (-1d0) + gb_mol ** (-1d0)) ** (-1d0)

    ! pp and qq represent limitation by metabolic (temperature & N) and
    ! diffusion (co2 supply) respectively
    pp = (pn*gC_to_umol)/gc ; qq = co2_comp_point-co2_half_sat
    ! calculate internal CO2 concentration (ppm or umol/mol)
    mult = co2+qq-pp
    ci = 0.5d0*(mult+sqrt((mult*mult)-4d0*(co2*qq-pp*co2_comp_point)))
    ci = min(ci,co2) ! C3 can't have more CO2 than is in the atmosphere
    ci_global = ci
    ! calculate CO2 limited rate of photosynthesis (gC.m-2.day-1)
    pd = (gc * (co2-ci)) * umol_to_gC
    ! scale to day light period as this is then consistent with the light
    ! capture period (1/24 = 0.04166667)
    pd = pd * dayl_hours * 0.04166667d0

    !
    ! Light limited photosynthesis
    !

    ! calculate light limted rate of photosynthesis (gC.m-2.day-1)
    pl = e0 * canopy_par_MJday

    !
    ! CO2 and light co-limitation
    !

    ! calculate combined light and CO2 limited photosynthesis
    acm_gpp = pl*pd/(pl+pd)

    ! sanity check
    if (acm_gpp /= acm_gpp) acm_gpp = 0d0

    ! don't forget to return
    return

  end function acm_gpp
  !
  !----------------------------------------------------------------------
  !
  double precision function find_gs_iWUE(gs_in)

    ! Calculate CO2 limited photosynthesis as a function of metabolic limited
    ! photosynthesis (pn), atmospheric CO2 concentration and stomatal
    ! conductance (gs_in). Photosynthesis is calculated twice to allow for
    ! testing of senstivity to iWUE.

    ! arguments
    double precision, intent(in) :: gs_in

    ! local variables
    double precision :: gs_high, gs_store, &
                        gpp_high, gpp_low

    !!!!!!!!!!
    ! Optimise intrinsic water use efficiency
    !!!!!!!!!!

    ! estimate photosynthesis with current estimate of gs
    gpp_low = acm_gpp(gs_in)

    ! Increment gs
    gs_high = gs_in + delta_gs
    ! estimate photosynthesis with incremented gs
    gpp_high = acm_gpp(gs_high)

    ! determine impact of gs increment on pd and how far we are from iWUE
    find_gs_iWUE = iWUE - ((gpp_high - gpp_low) * lai_1)

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
    double precision, intent(in) :: gs_in

    ! local variables
    double precision :: gs_high, gs_store, &
                        gpp_high, gpp_low, &
                        evap_high, evap_low

    !!!!!!!!!!
    ! Optimise water use efficiency
    !!!!!!!!!!

    ! Globally stored upper stomatal conductance estimate in memory
    gs_store = stomatal_conductance

    ! now assign the current estimate
    stomatal_conductance = gs_in
    ! estimate photosynthesis with current estimate of gs
    gpp_low = acm_gpp(gs_in)
    call calculate_transpiration(evap_low)

    ! Increment gs
    gs_high = gs_in + delta_gs
    ! now assign the incremented estimate
    stomatal_conductance = gs_high
    ! estimate photosynthesis with incremented gs
    gpp_high = acm_gpp(gs_high)
    call calculate_transpiration(evap_high)

    ! estimate marginal return on GPP for water loss, less water use efficiency criterion (gC.kgH2O-1.m-2.s-1)
    find_gs_WUE = ((gpp_high - gpp_low)/(evap_high - evap_low)) * lai_1
    find_gs_WUE = find_gs_WUE - iWUE

    ! return original stomatal value back into memory
    stomatal_conductance = gs_store

    ! remember to return back to the user
    return

  end function find_gs_WUE
  !
  !------------------------------------------------------------------
  !
  subroutine acm_albedo_gc(deltaWP,Rtot)

    ! Determines 1) an approximation of canopy conductance (gc) mmolH2O.m-2.s-1
    ! based on potential hydraulic flow, air temperature and absorbed radiation.
    ! 2) calculates absorbed shortwave radiation (W.m-2) as function of LAI

    implicit none

    ! arguments
    double precision, intent(in) :: deltaWP, & ! minlwp-wSWP (MPa)
                                       Rtot    ! total hydraulic resistance (MPa.s-1.m-2.mmol-1)

    ! local variables
    double precision :: denom
    double precision, parameter :: max_gs = 500d0, &   ! mmolH2O/m2leaf/s
                                   min_gs = 0.001d0, & ! mmolH2O/m2leaf/s
                                   tol_gs = 4d0        ! mmolH2O/m2leaf/s

    !!!!!!!!!!
    ! Calculate stomatal conductance under H2O and CO2 limitations
    !!!!!!!!!!

    if (deltaWP > vsmall) then
      ! Determine potential water flow rate (mmolH2O.m-2.dayl-1)
      max_supply = (deltaWP/Rtot) * seconds_per_day
    else
      ! set minimum (computer) precision level flow
      max_supply = vsmall
    end if

    if (aerodynamic_conductance > vsmall) then

      ! there is lai therefore we have have stomatal conductance

      ! Invert Penman-Monteith equation to give gs (m.s-1) needed to meet
      ! maximum possible evaporation for the day.
      ! This will then be reduced based on CO2 limits for diffusion based
      ! photosynthesis
      denom = slope * ((canopy_swrad_MJday * 1d6 * dayl_seconds_1) + canopy_lwrad_Wm2) &
           + (ET_demand_coef * aerodynamic_conductance)
      denom = (denom / (lambda * max_supply * mmol_to_kg_water * dayl_seconds_1)) - slope
      denom = denom / psych
      stomatal_conductance = aerodynamic_conductance / denom

      ! convert m.s-1 to mmolH2O.m-2.s-1
      stomatal_conductance = stomatal_conductance * convert_ms1_mmol_1
      ! if conditions are dew forming then set conductance to maximum as we are not going to be limited by water demand
      if (stomatal_conductance <= 0d0 .or. stomatal_conductance > max_gs) stomatal_conductance = max_gs

      ! if we are potentially limited by stomatal conductance or we are using instrinsic water use efficiency (rather than WUE)
      ! then iterate to find optimum gs otherwise just go with the max...
      if (stomatal_conductance /= max_gs .or. do_iWUE ) then
        ! If there is a positive demand for water then we will solve for photosynthesis limits on gs through iterative solution
        delta_gs = 1d-3 * lai ! mmolH2O/m2leaf/day
        ! estimate inverse of LAI to avoid division in optimisation
        lai_1 = lai**(-1d0)
        if (do_iWUE) then
          ! intrinsic WUE optimisation
          stomatal_conductance = zbrent('acm_albedo_gc:find_gs_iWUE',find_gs_iWUE,min_gs,stomatal_conductance,tol_gs)
        else
          ! WUE optimisation
          stomatal_conductance = zbrent('acm_albedo_gc:find_gs_WUE',find_gs_WUE,min_gs,stomatal_conductance,tol_gs)
        endif
      end if

    else

      ! if no LAI then there can be no stomatal conductance
      stomatal_conductance = 0d0

    endif ! if LAI > vsmall

  end subroutine acm_albedo_gc
  !
  !------------------------------------------------------------------
  !
  subroutine meteorological_constants(input_temperature,input_temperature_K)

    ! Determine some multiple use constants used by a wide range of functions
    ! All variables here are linked to air temperature and thus invarient between
    ! iterations and can be stored in memory...

    implicit none

    ! arguments
    double precision, intent(in) :: input_temperature, input_temperature_K

    ! local variables
    double precision :: s, mult

    !
    ! Used for soil, canopy evaporation and transpiration
    !

    ! Density of air (kg/m3)
    air_density_kg = 353d0/input_temperature_K
    ! Conversion ratio for m.s-1 -> mol.m-2.s-1
    convert_ms1_mol_1 = const_sfc_pressure / (input_temperature_K*Rcon)
    ! latent heat of vapourisation,
    ! function of air temperature (J.kg-1)
    if (input_temperature < 1d0) then
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
    s = 2502.935945d0*exp(17.269d0*input_temperature/mult)
    ! Rate of change of saturation vapour pressure with temperature (kPa.K-1)
    slope = s/(mult*mult)

    !
    ! Used for soil evaporation and leaf level conductance
    !

    ! Determine diffusion coefficient (m2.s-1), temperature dependant (pressure dependence neglected). Jones p51; appendix 2
    ! Temperature adjusted from standard 20oC (293.15 K), NOTE that 1/293.15 = 0.003411223
    ! 0.0000242 = conversion to make diffusion specific for water vapor (um2.s-1)
    water_vapour_diffusion = 0.0000242d0*((input_temperature_K/293.15d0)**1.75d0)

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
    double precision, intent(out) :: transpiration ! kgH2O.m-2.day-1

    ! local variables
    double precision :: canopy_radiation & ! isothermal net radiation (W/m2)
                           ,water_supply & ! Potential water supply to canopy from soil (kgH2O.m-2.day-1)
                                  ,gs,gb   ! stomatal and boundary layer conductance (m.s-1)

    !!!!!!!!!!
    ! Estimate energy radiation balance (W.m-2)
    !!!!!!!!!!

    ! Absorbed shortwave radiation MJ.m-2.day-1 -> J.m-2.s-1
    canopy_radiation = canopy_lwrad_Wm2 + (canopy_swrad_MJday * 1d6 * dayl_seconds_1)

    !!!!!!!!!!
    ! Calculate canopy conductance (to water vapour)
    !!!!!!!!!!

    ! calculate potential water supply (kgH2O.m-2.day-1)
    ! provided potential upper bound on evaporation
    water_supply = max_supply * mmol_to_kg_water

    ! Change units of potential stomatal conductance
    ! (mmolH2O.m-2.s-1 -> m.s-1).
    ! Note assumption of sea surface pressure only
    gs = stomatal_conductance / convert_ms1_mmol_1
    ! Combine in series stomatal conductance with boundary layer
    gb = aerodynamic_conductance

    !!!!!!!!!!
    ! Calculate canopy evaporative fluxes (kgH2O/m2/day)
    !!!!!!!!!!

    ! Calculate numerator of Penman Montheith (kg.m-2.day-1)
    transpiration = (slope*canopy_radiation) + (ET_demand_coef*gb)
    ! Calculate the transpiration flux and restrict by potential water supply
    ! over the day
    transpiration = min(water_supply,(transpiration / (lambda*(slope+(psych*(1d0+gb/gs)))))*dayl_seconds)

  end subroutine calculate_transpiration
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_wetcanopy_evaporation(wetcanopy_evap,act_pot_ratio,storage,transpiration)

    ! Estimates evaporation of canopy intercepted rainfall based on the Penman-Monteith model of
    ! evapotranspiration used to estimate SPA's daily evapotranspiration flux
    ! (kgH20.m-2.day-1).

    implicit none

    ! arguments
    double precision, intent(in) :: transpiration      ! kgH2O/m2/day
    double precision, intent(inout) :: storage         ! canopy water storage kgH2O/m2
    double precision, intent(out) :: wetcanopy_evap, & ! kgH2O.m-2.day-1
                                      act_pot_ratio    ! Ratio of potential evaporation to actual

    ! local variables
    double precision :: canopy_radiation, & ! isothermal net radiation (W/m2)
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
    canopy_radiation = canopy_lwrad_Wm2 + (canopy_swrad_MJday * 1d6 * seconds_per_day_1)

    !!!!!!!!!!
    ! Calculate canopy evaporative fluxes (kgH2O/m2/day)
    !!!!!!!!!!

    ! Calculate numerator of Penman Montheith (kgH2O.m-2.day-1)
    wetcanopy_evap = (slope*canopy_radiation) + (ET_demand_coef*gb)
    ! Calculate the potential wet canopy evaporation, limited by energy used for
    ! transpiration
    wetcanopy_evap = (wetcanopy_evap / (lambda*(slope+psych))) * seconds_per_day
    ! substract transpiration from potential surface evaporation
    wetcanopy_evap = wetcanopy_evap - transpiration

    ! Dew is unlikely to occur (if we had energy balance) if mint > 0
    ! Sublimation is also unlikely to occur (if we had energy balance) if maxt < 0
    if ((wetcanopy_evap < 0d0 .and. mint > 0d0) .or. &
        (wetcanopy_evap > 0d0 .and. maxt < 0d0)) then
        wetcanopy_evap = 0d0
    endif

    ! dew is unlikely to occur (if we had energy balance) if mint > 0
    !    if (wetcanopy_evap < 0d0 .and. mint > 0d0) wetcanopy_evap = 0d0
    ! Sublimation is unlikely to occur (if we had energy balance) if maxt < 0
    !    if (wetcanopy_evap > 0d0 .and. maxt < 0d0) wetcanopy_evap = 0d0

    ! Remember potential evaporation to later calculation of the potential
    ! actual ratio
    act_pot_ratio = wetcanopy_evap

    ! assuming there is any rainfall, currently water on the canopy or dew formation
    if (rainfall > 0d0 .or. storage > 0d0 .or. wetcanopy_evap < 0d0) then
      ! Update based on canopy water storage
      call canopy_interception_and_storage(wetcanopy_evap,storage)
    else
      ! there is no water movement possible
      intercepted_rainfall = 0d0 ; wetcanopy_evap = 0d0
    endif

    ! now calculate the ratio of potential to actual evaporation
    if (act_pot_ratio == 0d0) then
      act_pot_ratio = 0d0
    else
      act_pot_ratio = abs(wetcanopy_evap / act_pot_ratio)
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
    double precision, intent(out) :: soilevap ! kgH2O.m-2.day-1

    ! local variables
    double precision :: local_temp     &
                       ,soil_radiation & ! isothermal net radiation (W/m2)
                                ,esurf & ! see code below
                                 ,esat & ! soil air space saturation vapour pressure
                                  ,gws & ! water vapour conductance through soil air space (m.s-1)
                                   ,Qc

    local_temp = maxt + freeze

    !!!!!!!!!!
    ! Estimate energy radiation balance (W.m-2)
    !!!!!!!!!!

    ! Absorbed shortwave radiation MJ.m-2.day-1 -> J.m-2.s-1
    soil_radiation = soil_lwrad_Wm2 + (soil_swrad_MJday * 1d6 * dayl_seconds_1)
    ! estimate ground heat flux from statistical approximation, positive if energy moving up profile
    ! NOTE: linear coefficient estimates from SPA simulations
    Qc = -0.4108826d0 * (maxt - maxt_lag1)
    soil_radiation = soil_radiation + Qc

    !!!!!!!!!!
    ! Calculate soil evaporative fluxes (kgH2O/m2/day)
    !!!!!!!!!!

    ! calculate saturated vapour pressure (kPa), function of temperature.
    esat = 0.1d0 * exp( 1.80956664d0 + ( 17.2693882d0 * local_temp - 4717.306081d0 ) / ( local_temp - 35.86d0 ) )
    air_vapour_pressure = esat - vpd_kPa

    ! Soil conductance to water vapour diffusion (m s-1)...
    gws = porosity(1) * water_vapour_diffusion / (tortuosity*drythick)

    ! vapour pressure in soil airspace (kPa), dependent on soil water potential
    ! - Jones p.110. partial_molar_vol_water
    esurf = esat * exp( 1d6 * SWP(1) * partial_molar_vol_water / ( Rcon * local_temp ) )
    ! now difference in vapour pressure between soil and canopy air spaces
    esurf = esurf - air_vapour_pressure

    ! Estimate potential soil evaporation flux (kgH2O.m-2.day-1)
    soilevap = (slope*soil_radiation) + (air_density_kg*cpair*esurf*soil_conductance)
    soilevap = soilevap / (lambda*(slope+(psych*(1d0+soil_conductance/gws))))
    soilevap = soilevap * dayl_seconds

    ! dew is unlikely to occur (if we had energy balance) if mint > 0
    !    if (soilevap < 0d0 .and. mint > 1d0) soilevap = 0d0
    ! Sublimation is unlikely to occur (if we had energy balance) if maxt < 0
    !    if (soilevap > 0d0 .and. maxt < 1d0) soilevap = 0d0

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
    double precision :: local_lai, &
           mixing_length_momentum, & ! mixing length parameter for momentum (m)
            length_scale_momentum    ! length scale parameter for momentum (m)

    ! calculate the zero plane displacement and roughness length
    call z0_displacement(ustar_Uh)
    ! calculate friction velocity at tower height (reference height ) (m.s-1)
    ! WARNING neutral conditions only; see WRF module_sf_sfclay.F for 'with
    ! stability versions'
    !    ustar = (wind_spd / log((tower_height-displacement)/roughl)) * vonkarman
    ustar = wind_spd * ustar_Uh
    ! both length scale and mixing length are considered to be constant within
    ! the canopy (under dense canopy conditions) calculate length scale (lc)
    ! for momentum absorption within the canopy; Harman & Finnigan (2007)
    ! and mixing length (lm) for vertical momentum within the canopy Harman & Finnigan (2008)
    local_lai = max(min_lai,lai)
    length_scale_momentum = (4d0*canopy_height) / local_lai
    mixing_length_momentum = 2d0*(ustar_Uh**3)*length_scale_momentum

    ! based on Harman & Finnigan (2008); neutral conditions only
    call log_law_decay

    ! now we are interested in the within canopy wind speed,
    ! here we assume that the wind speed just inside of the canopy is most important.
    canopy_wind = canopy_wind*exp((ustar_Uh*((canopy_height*0.75d0)-canopy_height))/mixing_length_momentum)

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
    double precision, intent(out) :: gv_forced ! canopy conductance (m/s) for water vapour under forced convection

    ! local parameters
    double precision, parameter :: leaf_width = 0.04d0, & ! leaf width (m) (original 0.08)
                                 leaf_width_1 = leaf_width ** (-1d0), &
                                           Pr = 0.72d0, & ! Prandtl number
                                      Pr_coef = 1.05877d0 !1.18d0*(Pr**(0.33d0))
    ! local variables
    double precision :: &
         nusselt_forced & ! Nusselt value under forced convection
             ,Sh_forced & ! Sherwood number under forced convection
                    ,Re   ! Reynolds number

    ! Reynold number
    Re = (leaf_width*canopy_wind)/kinematic_viscosity
    ! calculate nusselt value under forced convection conditions
!    nusselt_forced = (1.18d0*(Pr**(0.33d0))*(sqrt(Re)))
    nusselt_forced = Pr_coef*(sqrt(Re))
    ! update specific Sherwood numbers
    Sh_forced = 0.962d0*nusselt_forced
    ! Estimate the the forced conductance of water vapour
    gv_forced = ((water_vapour_diffusion*Sh_forced)*leaf_width_1) * 0.5d0 * lai

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
    canopy_wind = (ustar * vonkarman_1) * log((canopy_height-displacement) / roughl)

    ! set minimum value for wind speed at canopy top (m.s-1)
    canopy_wind = max(min_wind,canopy_wind)

  end subroutine log_law_decay
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_field_capacity

    ! field capacity calculations for saxton eqns !

    implicit none

    ! local variables..
    integer        :: i
    double precision :: x1, x2

    x1 = 0.1d0 ; x2 = 0.7d0 ! low/high guess
    do i = 1 , nos_soil_layers+1
      water_retention_pass = i
      ! field capacity is water content at which SWP = -10 kPa
      field_capacity(i) = zbrent('water_retention:water_retention_saxton_eqns', &
      water_retention_saxton_eqns , x1 , x2 , 0.001d0 )
    enddo

  end subroutine calculate_field_capacity
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_daylength(doy,lat)

    ! Subroutine uses day of year and latitude (-90 / 90 degrees) as inputs,
    ! combined with trigonomic functions to calculate day length in hours and seconds

    implicit none

    ! arguments
    double precision, intent(in) :: doy, lat

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
    dayl_hours = 12d0 * ( 1d0 + 2d0 * asin( aob ) * pi_1 )
    dayl_seconds = dayl_hours * seconds_per_hour

    ! return to user
    return

  end subroutine calculate_daylength
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_longwave_isothermal(canopy_temperature,soil_temperature)

    ! Subroutine estimates the isothermal net longwave radiation (W.m-2) for
    ! the canopy and soil surface. SPA uses a complex multi-layer radiative
    ! transfer scheme including reflectance, transmittance any absorption.
    ! However, for a given canopy vertical profiles, the LAI absorption
    ! relationship is readily predicted via Michaelis-Menten or
    ! non-rectangular hyperbola as done here.

    implicit none

    ! arguments
    double precision, intent(in) :: canopy_temperature, soil_temperature ! oC

    ! local variables
    double precision :: lwrad, & ! downward long wave radiation from sky (W.m-2)
        longwave_release_soil, & ! emission of long wave radiation from surfaces per m2
      longwave_release_canopy, & ! assuming isothermal condition (W.m-2)
            trans_lw_fraction, &
        reflected_lw_fraction, &
         absorbed_lw_fraction, &
      canopy_release_fraction, & ! fraction of longwave emitted from within the canopy to ultimately be released
   canopy_absorption_from_sky, & ! canopy absorbed radiation from downward LW (W.m-2)
  canopy_absorption_from_soil, & ! canopy absorbed radiation from soil surface (W.m-2)
                  canopy_loss, & ! longwave radiation released from canopy surface (W.m-2).
                                 ! i.e. this value is released from the top and the bottom
     soil_absorption_from_sky, & ! soil absorbed radiation from sky (W.m-2)
  soil_absorption_from_canopy    ! soil absorbed radiation emitted from canopy (W.m-2)

    ! estimate long wave radiation from atmosphere (W.m-2)
    lwrad = emiss_boltz * (maxt+freeze-20d0) ** 4
    ! estimate isothermal long wave emission per unit area
    longwave_release_soil = emiss_boltz * (soil_temperature+freeze) ** 4
    ! estimate isothermal long wave emission per unit area
    longwave_release_canopy = emiss_boltz * (canopy_temperature+freeze) ** 4

    !!!!!!!!!!
    ! Determine fraction of longwave absorbed by canopy and returned to the sky
    !!!!!!!!!!

    ! calculate fraction of longwave radiation coming from the sky to pentrate to the soil surface
    trans_lw_fraction = 1d0 - (max_lai_lwrad_transmitted*lai)/(lai+lai_half_lwrad_transmitted)
    ! calculate the fraction of longwave radiation from sky which is reflected back into the sky
    reflected_lw_fraction = (max_lai_lwrad_reflected*lai) / (lai+lai_half_lwrad_reflected)
    ! calculate absorbed longwave radiation coming from the sky
    absorbed_lw_fraction = 1d0 - trans_lw_fraction - reflected_lw_fraction
    ! Calculate the potential absorption of longwave radiation lost from the
    ! canopy to soil / sky
    canopy_release_fraction = 1d0 - (max_lai_lwrad_release*lai) / (lai+lai_half_lwrad_release)

    !!!!!!!!!!
    ! Distribute longwave from sky
    !!!!!!!!!!

    ! long wave absorbed by the canopy from the sky
    canopy_absorption_from_sky = lwrad * absorbed_lw_fraction
    ! Long wave absorbed by soil from the sky, soil absorption assumed to be equal to emissivity
    soil_absorption_from_sky = trans_lw_fraction * lwrad * emissivity
    ! Long wave reflected directly back into sky
    sky_lwrad_Wm2 = lwrad * reflected_lw_fraction

    !!!!!!!!!!
    ! Distribute longwave from soil
    !!!!!!!!!!

    ! First, calculate longwave radiation coming up from the soil plus the radiation which is reflected
    canopy_absorption_from_soil = longwave_release_soil + (trans_lw_fraction * lwrad * (1d0-emissivity))
    ! Second, use this total to estimate the longwave returning to the sky
    sky_lwrad_Wm2 = sky_lwrad_Wm2 + (canopy_absorption_from_soil * trans_lw_fraction)
    ! Third, now calculate the longwave from the soil surface absorbed by the canopy
    canopy_absorption_from_soil = canopy_absorption_from_soil * absorbed_lw_fraction

    !!!!!!!!!!
    ! Distribute longwave originating from the canopy itself
    !!!!!!!!!!

    ! calculate two-sided long wave radiation emitted from canopy which is
    ! ultimately lost from to soil or sky (i.e. this value is used twice, once
    ! to soil once to sky)
    canopy_loss = longwave_release_canopy * lai * canopy_release_fraction
    ! Calculate longwave absorbed by soil which is released by the canopy itself
    soil_absorption_from_canopy = canopy_loss * emissivity
    ! Canopy released longwave returned to the sky
    sky_lwrad_Wm2 = sky_lwrad_Wm2 + canopy_loss

    !!!!!!!!!!
    ! Isothermal net long wave canopy and soil balance (W.m-2)
    !!!!!!!!!!

    ! determine isothermal net canopy. Note two canopy_loss used to account for
    ! upwards and downwards emissions
    canopy_lwrad_Wm2 = (canopy_absorption_from_sky + canopy_absorption_from_soil) - (canopy_loss + canopy_loss)
    ! determine isothermal net soil
    soil_lwrad_Wm2 = (soil_absorption_from_sky + soil_absorption_from_canopy) - longwave_release_soil

  end subroutine calculate_longwave_isothermal
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_shortwave_balance

    ! Subroutine estimates the canopy and soil absorbed shortwave radiation (MJ/m2/day).
    ! Radiation absorption is paritioned into NIR and PAR for canopy, and NIR +
    ! PAR for soil.

    ! SPA uses a complex multi-layer radiative transfer scheme including
    ! reflectance, transmittance any absorption. However, for a given
    ! canopy vertical profiles, the LAI absorption relationship is readily
    ! predicted via Michaelis-Menten or non-rectangular hyperbola as done here.

    implicit none

    ! local variables
    double precision :: balance                    &
                       ,absorbed_nir_fraction_soil &
                       ,absorbed_par_fraction_soil &
                       ,fsnow                      &
                       ,soil_par_MJday             &
                       ,soil_nir_MJday             &
                       ,trans_nir_MJday            &
                       ,trans_par_MJday            &
                       ,canopy_nir_MJday           &
                       ,refl_par_MJday             &
                       ,refl_nir_MJday             &
                       ,reflected_nir_fraction     & !
                       ,reflected_par_fraction     & !
                       ,absorbed_nir_fraction      & !
                       ,absorbed_par_fraction      & !
                       ,trans_nir_fraction         & !
                       ,trans_par_fraction

    ! local parameters
    double precision, parameter :: newsnow_nir_abs = 0.27d0 & ! NIR absorption fraction
                                  ,newsnow_par_abs = 0.05d0   ! PAR absorption fraction

    !!!!!!!!!!
    ! Determine canopy absorption / reflectance as function of LAI
    !!!!!!!!!!

    ! Canopy transmitted of PAR & NIR radiation towards the soil
    trans_par_fraction = 1d0 - (lai*max_lai_par_transmitted) &
                       / (lai+lai_half_par_transmitted)
    trans_nir_fraction = 1d0 - (lai*max_lai_nir_transmitted) &
                       / (lai+lai_half_nir_transmitted)
    ! Canopy reflected of near infrared and photosynthetically active radiation
    reflected_nir_fraction = (lai*max_lai_nir_reflection) &
                           / (lai+lai_half_nir_reflection)
    reflected_par_fraction = (lai*max_lai_par_reflection) &
                           / (lai+lai_half_par_reflection)
    ! Canopy absorption of near infrared and photosynthetically active radiation
    absorbed_nir_fraction = 1d0 - reflected_nir_fraction - trans_nir_fraction
    absorbed_par_fraction = 1d0 - reflected_par_fraction - trans_par_fraction

    !!!!!!!!!!
    ! Estimate canopy absorption of incoming shortwave radiation
    !!!!!!!!!!

    ! Estimate incoming shortwave radiation absorbed, transmitted and reflected by the canopy (MJ.m-2.day-1)
    canopy_par_MJday = (sw_par_fraction * swrad * absorbed_par_fraction)
    canopy_nir_MJday = ((1d0 - sw_par_fraction) * swrad * absorbed_nir_fraction)
    trans_par_MJday = (sw_par_fraction * swrad * trans_par_fraction)
    trans_nir_MJday = ((1d0 - sw_par_fraction) * swrad * trans_nir_fraction)
    refl_par_MJday = (sw_par_fraction * swrad * reflected_par_fraction)
    refl_nir_MJday = ((1d0 - sw_par_fraction) * swrad * reflected_nir_fraction)

    !!!!!!!!!
    ! Estimate soil absorption of shortwave passing through the canopy
    !!!!!!!!!

    ! Update soil reflectance based on snow cover
    if (snow_storage > 0d0) then
      fsnow = 1d0 - exp( - snow_storage * 1d-2 )  ! fraction of snow cover on the ground
      absorbed_par_fraction_soil = ((1d0 - fsnow) * soil_swrad_absorption) + (fsnow * newsnow_par_abs)
      absorbed_nir_fraction_soil = ((1d0 - fsnow) * soil_swrad_absorption) + (fsnow * newsnow_nir_abs)
    else
      absorbed_par_fraction_soil = soil_swrad_absorption
      absorbed_nir_fraction_soil = soil_swrad_absorption
    endif

    ! Then the radiation incident and ultimately absorbed by the soil surface itself (MJ.m-2.day-1)
    soil_par_MJday = trans_par_MJday * absorbed_par_fraction_soil
    soil_nir_MJday = trans_nir_MJday * absorbed_nir_fraction_soil
    ! combine totals for use is soil evaporation
    soil_swrad_MJday = soil_nir_MJday + soil_par_MJday

    !!!!!!!!!
    ! Estimate canopy absorption of soil reflected shortwave radiation
    ! This additional reflection / absorption cycle is needed to ensure > 0.99
    ! of incoming radiation is explicitly accounted for in the energy balance.
    !!!!!!!!!

    ! Update the canopy radiation absorption based on the reflected radiation (MJ.m-2.day-1)
    canopy_par_MJday = canopy_par_MJday + ((trans_par_MJday-soil_par_MJday) * absorbed_par_fraction)
    canopy_nir_MJday = canopy_nir_MJday + ((trans_nir_MJday-soil_nir_MJday) * absorbed_nir_fraction)
    ! Update the total radiation reflected back into the sky, i.e. that which is now transmitted through the canopy
    refl_par_MJday = refl_par_MJday + ((trans_par_MJday-soil_par_MJday) * trans_par_fraction)
    refl_nir_MJday = refl_nir_MJday + ((trans_nir_MJday-soil_nir_MJday) * trans_nir_fraction)

    ! Combine to estimate total shortwave canopy absorbed radiation
    canopy_swrad_MJday = canopy_par_MJday + canopy_nir_MJday

    ! check energy balance
!    balance = swrad - canopy_par_MJday - canopy_nir_MJday - refl_par_MJday - refl_nir_MJday - soil_swrad_MJday
!    if ((balance - swrad) / swrad > 0.01) then
!        print*,"SW residual frac = ",(balance - swrad) / swrad,"SW residual = ",balance,"SW in = ",swrad
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
    double precision,intent(inout) :: Rtot ! MPa.s-1.m-2.mmol-1

    ! local variables
    integer :: i
    double precision :: bonus, sum_water_flux, &
    transpiration_resistance,root_reach_local, &
                                root_depth_50
    double precision, dimension(nos_root_layers) :: root_mass    &
                                                   ,root_length  &
                                                   ,ratio
    double precision, parameter :: root_depth_frac_50 = 0.25d0 ! fractional soil depth above which 50 %
                                                               ! of the root mass is assumed to be located

    ! reset water flux
    water_flux = 0d0 ; wSWP = 0d0
    ratio = 0d0 ; ratio(1) = 1d0 ; root_mass = 0d0
    ! calculate soil depth to which roots reach
    root_reach = max_depth * root_biomass / (root_k + root_biomass)
    ! calculate the plant hydraulic resistance component. Currently unclear
    ! whether this actually varies with height or whether tall trees have a
    ! xylem architecture which keeps the whole plant conductance (gplant) 1-10 (ish).
    !    transpiration_resistance = (gplant * lai)**(-1d0)
    transpiration_resistance = canopy_height / (gplant * max(min_lai,lai))

    !!!!!!!!!!!
    ! calculate current steps soil hydraulic conductivity
    !!!!!!!!!!!

    ! seperately calculate the soil conductivity as this applies to each layer
    do i = 1, nos_soil_layers
      call calculate_soil_conductivity(i,soil_waterfrac(i),soil_conductivity(i))
    end do ! soil layers

    !!!!!!!!!!!
    ! Calculate root profile
    !!!!!!!!!!!

    ! The original SPA src generates an exponential distribution which aims
    ! to maintain 50 % of root biomass in the top 25 % of the rooting depth.
    ! In a simple 3 root layer system this can be estimates more simply

    ! top 25 % of root profile
    root_depth_50 = root_reach * root_depth_frac_50
    if (root_depth_50 <= layer_thickness(1)) then

      ! Greater than 50 % of the fine root biomass can be found in the top
      ! soil layer

      ! Start by assigning all 50 % of root biomass to the top soil layer
      root_mass(1) = root_biomass * 0.5d0
      ! Then quantify how much additional root is found in the top soil layer
      ! assuming that the top 25 % depth is found somewhere within the top
      ! layer
      bonus = (root_biomass-root_mass(1)) &
            * (layer_thickness(1)-root_depth_50) / (root_reach - root_depth_50)
      root_mass(1) = root_mass(1) + bonus
      ! partition the remaining root biomass between the seconds and third
      ! soil layers
      if (root_reach > sum(layer_thickness(1:2))) then
        root_mass(2) = (root_biomass - root_mass(1)) &
                     * (layer_thickness(2)/(root_reach-layer_thickness(1)))
        root_mass(3) = root_biomass - sum(root_mass(1:2))
      else
        root_mass(2) = root_biomass - root_mass(1)
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
      root_mass(3) = root_biomass - sum(root_mass(1:2))

    else
      ! Greater than 50 % of fine root biomass stock spans across all three
      ! layers
      root_mass(1:2) = root_biomass * 0.5d0 * (layer_thickness(1:2)/root_depth_50)
!      root_mass(1) = root_biomass * (layer_thickness(1)/root_depth_50)
!      root_mass(2) = root_biomass * (layer_thickness(2)/root_depth_50)
!      root_mass(1:2) = root_mass(1:2) * 0.5d0
      root_mass(3) = root_biomass - sum(root_mass(1:2))

    endif
    ! now convert root mass into lengths
    root_length = root_mass * root_mass_length_coef_1
!    root_length = root_mass / (root_density * root_cross_sec_area)

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
        root_reach_local = min(root_reach,layer_thickness(i))
        ! calculate and accumulate steady state water flux in mmol.m-2.s-1
        water_flux(i) = plant_soil_flow(i,root_length(i),root_mass(i) &
                                       ,demand(i),root_reach_local,transpiration_resistance)
      else
        ! ...if there is not then we wont have any below...
        exit
      end if ! root present in current layer?
    end do ! nos_root_layers

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
      wSWP = -20d0 ; uptake_fraction = 0d0 ; uptake_fraction(1) = 1d0
    else
      ! calculate weighted SWP and uptake fraction
      wSWP = sum(SWP(1:nos_root_layers) * water_flux(1:nos_root_layers))
      uptake_fraction(1:nos_root_layers) = water_flux(1:nos_root_layers) / sum_water_flux
      wSWP = wSWP / sum_water_flux
    endif

    ! determine effective resistance (MPa.s-1.m-2.mmol-1)
    Rtot = sum(demand) / sum_water_flux

    ! finally convert transpiration flux (mmol.m-2.s-1)
    ! into kg.m-2.step-1 for consistency with ET in "calculate_update_soil_water"
    water_flux = water_flux * mmol_to_kg_water * seconds_per_step

    ! and return
    return

  end subroutine calculate_Rtot
  !
  !-----------------------------------------------------------------
  !
  subroutine canopy_interception_and_storage(potential_evaporation,storage)

    ! Simple daily time step integration of canopy rainfall interception, runoff
    ! and rainfall (kgH2O.m-2.s-1). NOTE: it is possible for intercepted rainfall to be
    ! negative if stored water running off into the soil is greater than
    ! rainfall (i.e. when leaves have died between steps)

    implicit none

    ! arguments
    double precision, intent(inout) :: storage, & ! canopy water storage (kgH2O/m2)
                         potential_evaporation    ! wet canopy evaporation (kgH2O.m-2.day-1),
                                                  ! enters as potential but leaves as water balance adjusted

    ! local variables
    integer :: i, hr
    double precision :: a, through_fall, max_storage, max_storage_1, daily_addition, wetcanopy_evaporation &
                       ,potential_drainage_rate ,drain_rate, evap_rate, initial_canopy, co_mass_balance, dx, tmp(3)
    ! local parameters
    double precision, parameter :: CanIntFrac = -0.5d0,     & ! Coefficient scaling rainfall interception fraction with LAI
                                  CanStorFrac = 0.1d0,      & ! Coefficient scaling canopy water storage with LAI
                                 RefDrainRate = 0.002d0,    & ! Reference drainage rate (mm/min; Rutter et al 1975)
                                  RefDrainLAI = 0.952381d0, & ! Reference drainage 1/LAI (m2/m2; Rutter et al 1975, 1/1.05)
                                 RefDrainCoef = 3.7d0,      & ! Reference drainage Coefficient (Rutter et al 1975)
                               RefDrainCoef_1 = RefDrainCoef ** (-1d0)

    ! hold initial canopy storage in memory
    initial_canopy = storage
    ! determine maximum canopy storage & through fall fraction
    !    through_fall = max(min_throughfall,exp(CanIntFrac*lai))
    through_fall = exp(CanIntFrac*lai)
    ! maximum canopy storage (mm); minimum is applied to prevent errors in
    ! drainage calculation. Assume minimum capacity due to wood
    max_storage = max(min_storage,CanStorFrac*lai) ; max_storage_1 = max_storage**(-1d0)
    ! potential intercepted rainfall (kgH2O.m-2.s-1)
    intercepted_rainfall = rainfall * (1d0 - through_fall)

    ! calculate drainage coefficients (Rutter et al 1975); Corsican Pine
    ! 0.002 is canopy specific coefficient modified by 0.002*(max_storage/1.05)
    ! where max_storage is the canopy maximum capacity (mm) (LAI based) and
    ! 1.05 is the original canopy capacitance
    a = log( RefDrainRate * ( max_storage * RefDrainLAI ) ) - RefDrainCoef * max_storage

    ! average rainfall intercepted by canopy (kgH2O.m-2.day-1)
    daily_addition = intercepted_rainfall * seconds_per_day

    ! reset cumulative variables
    through_fall = 0d0 ; wetcanopy_evaporation = 0d0
    drain_rate = 0d0 ; evap_rate = 0d0

    ! deal with rainfall additions first
    do i = 1, int(days_per_step)

      ! add rain to the canopy and overflow as needed
      storage = storage + daily_addition

      if (storage > max_storage) then

        if (potential_evaporation > 0d0) then

          ! assume co-access to available water above max_storage by both drainage and
          ! evaporation. Water below max_storage is accessable by evaporation only.

          ! Trapezium rule for approximating integral of drainage rate.
          ! Allows estimation of the mean drainage rate between starting point and max_storage,
          ! thus the time period appropriate for co-access can be quantified. NOTE 1440 = minutes / day
          dx = storage - ((storage + max_storage)*0.5d0)
          tmp(1) = a + (RefDrainCoef * storage)
          tmp(2) = a + (RefDrainCoef * max_storage)
          tmp(3) = a + (RefDrainCoef * (storage+dx))
          tmp = exp(tmp)
          potential_drainage_rate = 0.5d0 * dx * ((tmp(1) + tmp(2)) + 2d0 * tmp(3))
          potential_drainage_rate = potential_drainage_rate * 1440d0

          ! restrict evaporation and drainage to the quantity above max_storage
          evap_rate = potential_evaporation ; drain_rate = min(potential_drainage_rate,storage-max_storage)

          ! limit based on available water if total demand is greater than excess
          co_mass_balance = ((storage-max_storage) / (evap_rate + drain_rate))
          evap_rate = evap_rate * co_mass_balance ; drain_rate = drain_rate * co_mass_balance

          ! estimate evaporation from remaining water, less that already removed from storage and evaporation energy used
          evap_rate = evap_rate + min(potential_evaporation - evap_rate, storage - evap_rate - drain_rate)

        else

          ! load dew formation to the current local evap_rate variable
          evap_rate = potential_evaporation
          ! restrict drainage the quantity above max_storage, adding dew formation too
          drain_rate = (storage - evap_rate) - max_storage

        endif

      else

        ! no drainage just apply evaporation / dew formation fluxes directly
        drain_rate = 0d0 ; evap_rate = potential_evaporation
        if (evap_rate > 0d0) then
          ! evaporation restricted by fraction of surface actually covered
          ! in water
          evap_rate = evap_rate * storage * max_storage_1
          ! and the total amount of water
          evap_rate = min(evap_rate,storage)
        else
          ! then dew formation has occurred, if this pushes storage > max_storage add it to drainage
          drain_rate = max(0d0,(storage - evap_rate) - max_storage)
        endif ! evap_rate > 0
      endif ! storage > max_storage

      ! update canopy storage with water flux
      !storage = max(0d0,storage - evap_rate - drain_rate)
      storage = storage - evap_rate - drain_rate
      wetcanopy_evaporation = wetcanopy_evaporation + evap_rate
      through_fall = through_fall + drain_rate

    end do ! days

    ! correct intercepted rainfall rate to kgH2O.m-2.s-1
    intercepted_rainfall = intercepted_rainfall - (through_fall * days_per_step_1 * seconds_per_day_1)

!    ! sanity checks; note 1e-8 prevents precision errors causing flags
!    if (intercepted_rainfall > rainfall .or. storage < 0d0 .or. &
!       (wetcanopy_evaporation * days_per_step_1) > (1d-8 + initial_canopy + (rainfall*seconds_per_day)) ) then
!       print*,"Condition 1",intercepted_rainfall > rainfall
!       print*,"Condition 2",storage < 0d0
!       print*,"Condition 3",(wetcanopy_evaporation * days_per_step_1) > (1d-8 + initial_canopy + (rainfall*seconds_per_day))
!       print*,"storage (kgH2O/m2)",storage,"max_storage (kgH2O/m2)",max_storage,"initial storage (kgH2O/m2)", initial_canopy
!       print*,"rainfall (kgH2O/m2/day)", rainfall*seconds_per_day, "through_fall (kgH2O/m2/day)", (through_fall * days_per_step_1)
!       print*,"through_fall_total (kgH2O/m2/step)",through_fall
!       print*,"potential_evaporation (kgH2O/m2/day)",potential_evaporation
!       print*,"actual evaporation    (kgH2O/m2/day)",wetcanopy_evaporation * days_per_step_1
!       stop
!    endif

    ! average evaporative flux to daily rate (kgH2O/m2/day)
    potential_evaporation = wetcanopy_evaporation * days_per_step_1

    ! final clearance of canopy storage of version small values at the level of system precision
    if (storage < 10d0*vsmall) storage = 0d0

  end subroutine canopy_interception_and_storage
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_update_soil_water(ET_leaf,ET_soil,rainfall_in,corrected_ET)

    !
    ! Function updates the soil water status and layer thickness
    ! Soil water profile is updated in turn with evaporative losses,
    ! rainfall infiltration and gravitational drainage
    ! Root layer thickness is updated based on changes in the rooting depth from
    ! the previous step
    !

    implicit none

    ! arguments
    double precision, intent(in) :: ET_leaf,ET_soil & ! evapotranspiration estimate (kgH2O.m-2.day-1)
                                       ,rainfall_in   ! rainfall (kgH2O.m-2.day-1)
    double precision, intent(out) :: corrected_ET     ! water balance corrected evapotranspiration (kgH2O/m2/day)

    ! local variables
    integer :: day
    double precision ::  depth_change, water_change
    double precision, dimension(nos_root_layers) :: avail_flux, evaporation_losses

    ! reset soil water exchanges
    underflow = 0d0 ; runoff = 0d0 ; corrected_ET = 0d0

    ! to allow for smooth water balance integration carry this out at daily time step
    do day = 1, nint(days_per_step)

      !!!!!!!!!!
      ! Evaporative losses
      !!!!!!!!!!

      ! Assume leaf transpiration is drawn from the soil based on the
      ! update_fraction estimated in calculate_Rtot
      evaporation_losses = ET_leaf * uptake_fraction
      ! Assume all soil evaporation comes from the soil surface only
      evaporation_losses(1) = evaporation_losses(1) + ET_soil
      ! can not evaporate from soil more than is available (m -> mm)
      !      avail_flux = soil_waterfrac(1:nos_root_layers) * layer_thickness(1:nos_root_layers) * 1d3
      !      where (evaporation_losses > avail_flux) evaporation_losses = avail_flux * 0.999d0

      ! this will update the ET estimate outside of the function
      ! days_per_step corrections happens outside of the loop below
      corrected_ET = corrected_ET + sum(evaporation_losses)

      ! pass information to waterloss variable and zero watergain
      ! convert kg.m-2 (or mm) -> Mg.m-2 (or m)
      waterloss = 0d0 ; watergain = 0d0
      waterloss(1:nos_root_layers) = evaporation_losses(1:nos_root_layers)*1d-3
      ! update soil water status with evaporative losses
      soil_waterfrac(1:nos_soil_layers) = ((soil_waterfrac(1:nos_soil_layers)*layer_thickness(1:nos_soil_layers)) &
                                        + watergain(1:nos_soil_layers) - waterloss(1:nos_soil_layers)) &
                                        / layer_thickness(1:nos_soil_layers)
      ! reset soil water flux variables
      waterloss = 0d0 ; watergain = 0d0

      !!!!!!!!!!
      ! Gravitational drainage
      !!!!!!!!!!

      ! determine drainage flux between surface -> sub surface and sub surface
      call gravitational_drainage

      !!!!!!!!!!
      ! Rainfall infiltration drainage
      !!!!!!!!!!

      ! determine infiltration from rainfall (kgH2O/m2/step),
      ! if rainfall is probably liquid / soil surface is probably not frozen
      if (rainfall_in > 0d0) then
        call infiltrate(rainfall_in)
      else
        runoff = runoff + (rainfall_in * days_per_step_1)
      endif ! is there any rain to infiltrate?
      ! update soil profiles. Convert fraction into depth specific values (rather than m3/m3) then update fluxes
      soil_waterfrac(1:nos_soil_layers) = ((soil_waterfrac(1:nos_soil_layers)*layer_thickness(1:nos_soil_layers)) &
                                        + watergain(1:nos_soil_layers) - waterloss(1:nos_soil_layers)) &
                                        / layer_thickness(1:nos_soil_layers)

      ! mass balance check, at this point do not try and adjust evaporation to
      ! correct for lack of supply. Simply allow for drought in next time step
      ! instead...
      where (soil_waterfrac <= 0d0)
        soil_waterfrac = vsmall
      end where

    end do ! days_per_step

    ! apply time step correction kgH2O/m2/step -> kgH2O/m2/day
    corrected_ET = corrected_ET * days_per_step_1
    underflow = underflow * days_per_step_1
    runoff = runoff * days_per_step_1

    !!!!!!!!!!
    ! Update soil layer thickness
    !!!!!!!!!!

    depth_change = 0d0 ; water_change = 0d0
    ! if roots extent down into the bucket
    if (root_reach > (top_soil_depth+mid_soil_depth) .or. previous_depth > (top_soil_depth+mid_soil_depth)) then
      ! how much has root depth extended since last step?
      depth_change = root_reach - previous_depth

      ! if there has been an increase
      if (depth_change > 0.01d0 .and. root_reach > sum(layer_thickness(1:2))+min_layer) then

        ! determine how much water is within the new volume of soil
        water_change = soil_waterfrac(nos_soil_layers) * depth_change
        ! now assign that new volume of water to the deep rooting layer
        soil_waterfrac(nos_root_layers) = ((soil_waterfrac(nos_root_layers) * layer_thickness(nos_root_layers)) &
                                        + water_change) / (layer_thickness(nos_root_layers)+depth_change)
        ! explicitly update the soil profile if there has been rooting depth
        ! changes
        layer_thickness(1) = top_soil_depth ; layer_thickness(2) = mid_soil_depth
        layer_thickness(3) = max(min_layer,root_reach-sum(layer_thickness(1:2)))
        layer_thickness(4) = max_depth - sum(layer_thickness(1:3))

        ! keep track of the previous rooting depth
        previous_depth = root_reach

      elseif (depth_change < -0.01d0 .and. root_reach > layer_thickness(1)+min_layer) then

        ! determine how much water is lost from the old volume of soil
        water_change = soil_waterfrac(nos_root_layers) * abs(depth_change)
        ! now assign that new volume of water to the deep rooting layer
        soil_waterfrac(nos_soil_layers) = ((soil_waterfrac(nos_soil_layers) * layer_thickness(nos_soil_layers)) &
                                        + water_change) / (layer_thickness(nos_soil_layers)+abs(depth_change))

        ! explicitly update the soil profile if there has been rooting depth
        ! changes
        layer_thickness(1) = top_soil_depth ; layer_thickness(2) = mid_soil_depth
        layer_thickness(3) = max(min_layer,root_reach-sum(layer_thickness(1:2)))
        layer_thickness(4) = max_depth - sum(layer_thickness(1:3))

        ! keep track of the previous rooting depth
        previous_depth = root_reach

      else

        ! we don't want to do anything, just recycle the previous depth

      end if ! depth change

    else

      ! if we are outside of the range when we need to consider rooting depth changes keep track incase we move into a zone when we do
      previous_depth = root_reach

    end if ! root reach beyond top layer

    ! finally update soil water potential
    call soil_water_potential

!    ! sanity check for catastrophic failure
!    do soil_layer = 1, nos_soil_layers
!       if (soil_waterfrac(soil_layer) < 0d0 .and. soil_waterfrac(soil_layer) > -0.01d0) then
!           soil_waterfrac(soil_layer) = 0d0
!       endif
!       if (soil_waterfrac(soil_layer) < 0d0 .or. soil_waterfrac(soil_layer) /= soil_waterfrac(soil_layer)) then
!          print*,'ET',ET,"rainfall",rainfall_in
!          print*,'evaporation_losses',evaporation_losses
!          print*,"watergain",watergain
!          print*,"waterloss",waterloss
!          print*,'depth_change',depth_change
!          print*,"soil_waterfrac",soil_waterfrac
!          print*,"porosity",porosity
!          print*,"layer_thicknes",layer_thickness
!          print*,"Uptake fraction",uptake_fraction
!          print*,"max_depth",max_depth,"root_k",root_k,"root_reach",root_reach
!          print*,"fail" ; stop
!       endif
!    end do

    ! explicit return needed to ensure that function runs all needed code
    return

  end subroutine calculate_update_soil_water
  !
  !-----------------------------------------------------------------
  !
  subroutine infiltrate(rainfall_in)

    ! Takes surface_watermm and distrubutes it among top !
    ! layers. Assumes total infilatration in timestep.   !

    implicit none

    ! arguments
    double precision, intent(in) :: rainfall_in ! rainfall (kg.m-2.day-1)

    ! local argumemts
    integer :: i
    double precision :: add   & ! surface water available for infiltration (m)
                       ,wdiff   ! available space in a given soil layer for water to fill (m)

    ! convert rainfall water from mm -> m (or kgH2O.m-2.day-1 -> MgH2O.m-2.day-1)
    add = rainfall_in * 1d-3

    do i = 1 , nos_soil_layers
      ! determine the available pore space in current soil layer
      wdiff = (porosity(i)-soil_waterfrac(i))*layer_thickness(i)-watergain(i)+waterloss(i)
      ! is the input of water greater than available space
      ! if so fill and subtract from input and move on to the next
      ! layer
      if (add > wdiff) then
        ! if so fill and subtract from input and move on to the next layer
        watergain(i) = watergain(i) + wdiff
        add = add - wdiff
      else
        ! otherwise infiltate all in the current layer
        watergain(i) = watergain(i) + add
        add = 0d0 ; exit
      end if

    end do ! nos_soil_layers

    ! if after all of this we have some water left assume it is runoff (kgH2O.m-2.day-1)
    ! NOTE that runoff is reset outside of the daily soil loop
    runoff = runoff + (add * 1d3)

  end subroutine infiltrate
  !
  !-----------------------------------------------------------------
  !
  subroutine gravitational_drainage

    ! integrator for soil gravitational drainage !

    implicit none

    ! local variables..
    integer :: d, nos_integrate
    double precision  :: tmp1,tmp2,tmp3,dx &
                                   ,liquid & ! liquid water in local soil layer (m3/m3)
                               ,drainlayer & ! field capacity of local soil layer (m3/m3)
                                    ,unsat & ! unsaturated pore space in soil_layer below the current (m3/m3)
                                   ,change & ! absolute volume of water drainage in current layer (m3)
                                 ,drainage & ! drainage rate of current layer (m/day)
                              ,local_drain & ! drainage of current layer (m/nos_minutes)
                 ,iceprop(nos_soil_layers)

    ! local parameters
    integer, parameter :: nos_hours_per_day = 1440, nos_minutes = 60*8

    ! calculate soil ice proportion; at the moment
    ! assume everything liquid
    iceprop = 0d0
    ! except the surface layer in the mean daily temperature is < 0oC
    if (meant < 1d0) iceprop(1) = 1d0

    do soil_layer = 1, nos_soil_layers

      ! soil water capacity of the current layer
      drainlayer = field_capacity( soil_layer )
      ! liquid content of the soil layer
      liquid     = soil_waterfrac( soil_layer ) &
                 * ( 1d0 - iceprop( soil_layer ) )

      ! initial conditions; i.e. is there liquid water and more water than
      ! layer can hold
      if ( liquid > drainlayer ) then

!        ! Trapezium rule for approximating integral of drainage rate
!        dx = liquid - ((liquid + drainlayer)*0.5d0)
!        call calculate_soil_conductivity(soil_layer,liquid,tmp1)
!        call calculate_soil_conductivity(soil_layer,drainlayer,tmp2)
!        call calculate_soil_conductivity(soil_layer,(liquid+dx),tmp3)
!        drainage = 0.5d0 * dx * ((tmp1 + tmp2) + 2d0 * tmp3)
!        drainage = drainage * seconds_per_day
!        drainage = min(drainage,liquid - drainlayer)

        d = 1 ; nos_integrate = nos_hours_per_day / nos_minutes
        drainage = 0d0 ; local_drain = 0d0
        do while (d <= nos_integrate .and. liquid > drainlayer)
          ! estimate drainage rate (m/s)
          call calculate_soil_conductivity(soil_layer,liquid,local_drain)
          ! scale to total number of seconds in increment
          local_drain = local_drain * dble(nos_minutes * 60)
          local_drain = min(liquid-drainlayer,local_drain)
          liquid = liquid - local_drain
          drainage = drainage + local_drain
          d = d + 1
        end do ! integrate over time

        ! unsaturated volume of layer below (m3 m-2)
        unsat = max( 0d0 , ( porosity( soil_layer+1 ) - soil_waterfrac( soil_layer+1 ) ) &
              * layer_thickness( soil_layer+1 ) / layer_thickness( soil_layer ) )
        ! layer below cannot accept more water than unsat
        if ( drainage > unsat ) drainage = unsat
        ! water loss from this layer (m3)
        change = drainage * layer_thickness(soil_layer)
        ! update soil layer below with drained liquid
        watergain( soil_layer + 1 ) = watergain( soil_layer + 1 ) + change
        waterloss( soil_layer     ) = waterloss( soil_layer     ) + change

      end if ! some liquid water and drainage possible

    end do ! soil layers

    ! estimate drainage from bottom of soil column (kgH2O/m2/day)
    ! NOTES: that underflow is reset outside of the daily soil loop
    underflow = underflow + (waterloss(nos_soil_layers) * 1d3)

  end subroutine gravitational_drainage
  !
  !-----------------------------------------------------------------
  !
  subroutine soil_porosity(soil_frac_clay,soil_frac_sand)

    ! Porosity is estimated from Saxton equations. !

    implicit none

    ! arguments
    double precision, dimension(nos_soil_layers) :: soil_frac_clay &
                                                   ,soil_frac_sand
    ! local variables..
    double precision, parameter :: H = 0.332d0, &
                                   J = -7.251d-4, &
                                   K = 0.1276d0

    ! loop over soil layers..
    porosity(1:nos_soil_layers) = H + J * soil_frac_sand(1:nos_soil_layers) + &
                                K * log10(soil_frac_clay(1:nos_soil_layers))
    ! then assign same to core layer
    porosity(nos_soil_layers+1) = porosity(nos_soil_layers)

  end subroutine soil_porosity
  !
  !---------------------------------------------------------------------
  !
  subroutine initialise_soils(soil_frac_clay,soil_frac_sand)

    !
    ! Subroutine calculate the soil layers field capacities and sets the initial
    ! soil water potential set to field capacity
    !

    implicit none

    ! arguments
    double precision, dimension(nos_soil_layers) :: soil_frac_clay &
                                                   ,soil_frac_sand

    ! local variables
    integer :: i

    ! include some hardcoded boundaries for the Saxton equations
    where (soil_frac_sand < 5d0) soil_frac_sand = 5d0
    where (soil_frac_clay < 5d0) soil_frac_clay = 5d0
    where (soil_frac_clay > 60d0) soil_frac_clay = 60d0
    ! calculate soil porosity (m3/m3)
    call soil_porosity(soil_frac_clay,soil_frac_sand)
    ! calculate field capacity (m3/m-3)
    call calculate_field_capacity

    ! final sanity check for porosity
    where (porosity <= field_capacity) porosity = field_capacity + 0.05d0

  end subroutine initialise_soils
  !
  !---------------------------------------------------------------------
  !
  subroutine update_soil_initial_conditions(input_soilwater_frac)

    !
    ! Subroutine calculate the soil layers field capacities and sets the initial
    ! soil water potential set to field capacity
    !

    implicit none

    ! arguments
    double precision :: input_soilwater_frac

    ! local variables
    integer :: i

    ! Default assumption to be field capacity
    soil_waterfrac = field_capacity
    SWP = SWP_initial

    ! if prior value has been given
    if (input_soilwater_frac > -9998d0) then
      ! calculate initial soil water fraction
      soil_waterfrac(1:nos_soil_layers) = input_soilwater_frac
      ! calculate initial soil water potential
      call soil_water_potential
    endif

    ! seperately calculate the soil conductivity as this applies to each layer
    do i = 1, nos_soil_layers
      call calculate_soil_conductivity(i,soil_waterfrac(i),soil_conductivity(i))
    end do ! soil layers
    ! but apply the lowest soil layer to the core as well in initial conditions
    soil_conductivity(nos_soil_layers+1) = soil_conductivity(nos_soil_layers)

  end subroutine update_soil_initial_conditions
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_soil_conductivity(soil_layer,waterfrac,conductivity)

    ! Calculate the soil conductivity (m s-1) of water based on soil
    ! characteristics and current water content

    implicit none

    ! arguments
    integer, intent(in) :: soil_layer
    double precision, intent(in) :: waterfrac
    double precision, intent(out) :: conductivity

    ! soil conductivity for the dynamic soil layers (i.e. not including core)
    conductivity = cond1(soil_layer) * exp(cond2(soil_layer)+cond3(soil_layer)/waterfrac)

    ! protection against floating point error
    if (waterfrac < 0.05d0) conductivity = 1d-30

  end subroutine calculate_soil_conductivity
  !
  !------------------------------------------------------------------
  !
  subroutine saxton_parameters(soil_frac_clay,soil_frac_sand)

    ! Calculate the key parameters of the Saxton, that is cond1,2,3 !
    ! and potA,B                                                    !

    implicit none

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
    potA(1:nos_soil_layers) = A + (B * soil_frac_clay) + &
                             (CC * soil_frac_sand * soil_frac_sand) + &
                             (D * soil_frac_sand * soil_frac_sand * soil_frac_clay)
    potA(1:nos_soil_layers) = exp(potA(1:nos_soil_layers))
    potA(1:nos_soil_layers) = potA(1:nos_soil_layers) * mult1

    potB(1:nos_soil_layers) = E + (F * soil_frac_clay * soil_frac_clay) + &
                             (G * soil_frac_sand * soil_frac_sand * soil_frac_clay)

    cond1(1:nos_soil_layers) = mult2
    cond2(1:nos_soil_layers) = P + (Q * soil_frac_sand)
    cond3(1:nos_soil_layers) = R + (T * soil_frac_sand) + (U * soil_frac_clay) + &
                              (V * soil_frac_clay * soil_frac_clay)

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
  subroutine calculate_soil_conductance(lm)

    ! proceedsure to solve for soil surface resistance based on Monin-Obukov
    ! similarity theory stability correction momentum & heat are integrated
    ! through the under canopy space and canopy air space to the surface layer
    ! references are Nui & Yang 2004; Qin et al 2002
    ! NOTE: conversion to conductance at end

    implicit none

    ! declare arguments
    double precision, intent(in) :: lm

    ! local variables
    double precision :: canopy_decay & ! canopy decay coefficient for soil exchange
                       ,Kh_canht       ! eddy diffusivity at canopy height (m2.s-1)

    ! parameters
    double precision, parameter :: foliage_drag = 0.2d0 ! foliage drag coefficient

    ! calculate eddy diffusivity at the top of the canopy (m2.s-1)
    ! Kaimal & Finnigan 1994; for near canopy approximation
    Kh_canht = vonkarman*ustar*(canopy_height-displacement)

    ! calculate canopy decay coefficient with stability correction
    ! NOTE this is not consistent with canopy momentum decay done by Harman &
    ! Finnigan (2008)
    canopy_decay = sqrt((foliage_drag*canopy_height*max(min_lai,lai))/lm)

    ! approximation of integral for soil resistance
    soil_conductance = canopy_height/(canopy_decay*Kh_canht) &
                     * (exp(canopy_decay*(1d0-(soil_roughl/canopy_height)))- &
                        exp(canopy_decay*(1d0-((roughl+displacement)/canopy_height))))

    ! convert resistance (s.m-1) to conductance (m.s-1)
    soil_conductance = soil_conductance ** (-1d0)

  end subroutine calculate_soil_conductance
  !
  !----------------------------------------------------------------------
  !
  subroutine soil_water_potential

    ! Find SWP without updating waterfrac yet (we do that in !
    ! waterthermal). Waterfrac is m3 m-3, soilwp is MPa.     !

    implicit none

    ! reformulation aims to remove if statement within loop to hopefully improve
    ! optimisation
    SWP(1:nos_soil_layers) = -0.001d0 * potA(1:nos_soil_layers) &
                           * soil_waterfrac(1:nos_soil_layers)**potB(1:nos_soil_layers)
    where (SWP(1:nos_soil_layers) < -20d0) SWP(1:nos_soil_layers) = -20d0

  end subroutine soil_water_potential
  !
  !------------------------------------------------------------------
  !
  subroutine z0_displacement(ustar_Uh)

    ! dynamic calculation of roughness length and zero place displacement (m)
    ! based on canopy height and lai. Raupach (1994)

    implicit none

    ! arguments
    double precision, intent(out) :: ustar_Uh ! ratio of friction velocity over wind speed at canopy top
    ! local variables
    double precision  sqrt_cd1_lai, local_lai
    double precision, parameter :: cd1 = 7.5d0,   & ! Canopy drag parameter; fitted to data
                                    Cs = 0.003d0, & ! Substrate drag coefficient
                                    Cr = 0.3d0,   & ! Roughness element drag coefficient
    !                          ustar_Uh_max = 0.3,   & ! Maximum observed ratio of
                                                      ! (friction velocity / canopy top wind speed) (m.s-1)
        ustar_Uh_max = 1d0, ustar_Uh_min = 0.2d0, &
                                        Cw = 2d0, &  ! Characterises roughness sublayer depth (m)
                                     phi_h = 0.19314718056d0 ! Roughness sublayer influence function;

    ! describes the departure of the velocity profile from just above the
    ! roughness from the intertial sublayer log law


    ! assign new value to min_lai to avoid max min calls
    local_lai = max(min_lai,lai)
    sqrt_cd1_lai = sqrt(cd1 * local_lai)

    ! calculate displacement (m); assume minimum lai 1.0 or 1.5 as height is not
    ! varied
    displacement = (1d0-((1d0-exp(-sqrt_cd1_lai))/sqrt_cd1_lai))*canopy_height

    ! calculate estimate of ratio of friction velocity / canopy wind speed; with
    ! max value set at
    ustar_Uh = max(ustar_Uh_min,min(sqrt(Cs+Cr*local_lai*0.5d0),ustar_Uh_max))
    ! calculate roughness sublayer influence function;
    ! this describes the departure of the velocity profile from just above the
    ! roughness from the intertial sublayer log law
    ! phi_h = log(Cw)-1d0+Cw**(-1d0) ! DO NOT FORGET TO UPDATE IF Cw CHANGES

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
  subroutine calculate_NUE_decline(NUE_pot_decay,min_temp,max_temp,min_vpd,max_vpd)

    ! Subroutine calcualte the local GSI based on minimum time step
    ! temperature and vaoour pressure deficit. This GSI defines the decay in NUE
    ! of age classes past maturity

    implicit none

    ! declare arguments
    double precision, intent(in) :: NUE_pot_decay, &
                                         min_temp, &
                                         max_temp, &
                                          min_vpd, &
                                          max_vpd

    ! declare local variables
    integer :: start, a, b, c, d, age_by_time
    double precision :: Tfac, VPDfac, NUE_decay, tmp

    ! Calculate GSI style Components
    Tfac = 1d0 - ((mint - min_temp) / (max_temp - min_temp)) ! oC
    VPDfac = (vpd_kPa - min_vpd) / (max_vpd - min_vpd)       ! kPa
    ! restrict to 0-1
    Tfac = min(1d0,max(0d0,Tfac)) ; VPDfac = min(1d0,max(0d0,VPDfac))
!print*,"GSI",Tfac,VPDfac,Tfac*VPDfac
    ! calculate actual decay rate
    start = nint(canopy_maturation_lag+1d0)
    NUE_decay = NUE_pot_decay * VPDfac * Tfac

    ! apply decay to the canopy
    NUE_vector(start:oldest_leaf) = NUE_vector(start:oldest_leaf) - &
              (NUE_vector(start:oldest_leaf) * (1d0-(1d0-NUE_decay)**days_per_step)*days_per_step_1)

  end subroutine calculate_NUE_decline
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_leaf_dynamics(current_step,deltat,nodays,pot_leaf_growth &
                                    ,Rm_exponent,Rm_baseline,deltaWP,Rtot       &
                                    ,GPP_current,Rm_leaf,foliage,leaf_fall,leaf_growth)

    ! Subroutine determines whether leaves are growing or dying.
    ! 1) Update canopy mean age and the impact on PNUE
    ! 2) Performes marginal return calculation, including mean age updates

    ! Thomas & Williams (2014):
    ! A model using marginal efficiency of investment to analyze carbon and nitrogen interactions in terrestrial ecosystems
    ! (ACONITE Version 1), Geoscientific Model Development, doi: 10.5194/gmd-7-2015-2014
    !
    ! Xu et al., (2017):
    ! Variations of leaf longevity in tropical moist forests predicted by a trait-driven carbon optimality model,
    ! Ecology Letters, doi: 10.1111/ele.12804

    implicit none

    ! declare arguments
    integer, intent(in) :: nodays, current_step
    double precision, intent(in) :: deltat(nodays) & !
                                  ,pot_leaf_growth & !
                                      ,Rm_exponent &
                                      ,Rm_baseline &
                                      ,GPP_current & !
                                          ,foliage & !
                                          ,Rm_leaf & !
                                          ,deltaWP & !
                                             ,Rtot

    double precision, intent(inout) :: leaf_fall,leaf_growth

    ! declare local variables
    integer :: a, b, c, d, tmp_int, &
         death_point, oldest_point, &
                       age_by_time
    double precision :: infi      &
                       ,tmp       &
                       ,loss      &
                     ,life_remain &
                       ,deltaGPP  &
                       ,deltaRm   &
                       ,deltaC    &
                       ,old_GPP   &
                       ,alt_GPP   &
                   ,marginal_gain &
                   ,marginal_loss &
                       ,NUE_save  &
                       ,lai_save  &
                  ,canopy_lw_save &
                  ,canopy_sw_save &
                 ,canopy_par_save &
                    ,soil_lw_save &
                    ,soil_sw_save &
                    ,gs_save      &
                    ,ga_save      &
                 ,canopy_age_save

    ! save original values for re-allocation later
    canopy_lw_save = canopy_lwrad_Wm2 ; soil_lw_save  = soil_lwrad_Wm2
    canopy_sw_save = canopy_swrad_MJday ; canopy_par_save  = canopy_par_MJday
    soil_sw_save = soil_swrad_MJday ; gs_save = stomatal_conductance
    ga_save = aerodynamic_conductance
    lai_save = lai ; NUE_save = NUE ; canopy_age_save = canopy_age

    ! for infinity checks
    infi = 0d0

    ! first assume that nothing is happening
    marginal_gain = 0d0 ; marginal_loss = 0d0
    loss = 0d0 ; tmp = 0d0
    deltaC = 0d0
    leaf_fall = 0d0   ! leaf turnover
    leaf_growth = 0d0 ! leaf growth

    !
    ! Increment the age of the existing canopy
    !

    ! how many days in the current time step
    age_by_time = nint(deltat(current_step))

    ! reset counters
    b = 0 ; c = 0 ; d = 0
    ! is there space at the end of the age vector?
    if (oldest_leaf+age_by_time > size(canopy_age_vector)) then

      ! calculate needed vector location information
      d = size(canopy_age_vector) ; b = (oldest_leaf+age_by_time) - d ; c = 1

      ! We need to track the marginal return estimates for the age classes first...
      marginal_loss_avg(d) = sum( marginal_loss_avg((oldest_leaf-b):oldest_leaf) &
                           * canopy_age_vector((oldest_leaf-b):oldest_leaf) )
      NUE_vector(d) = sum(NUE_vector((oldest_leaf-b):oldest_leaf) &
                    * canopy_age_vector((oldest_leaf-b):oldest_leaf))

      ! ...as these need to be weighted by their biomass stocks...
      tmp = sum(canopy_age_vector((oldest_leaf-b):oldest_leaf))
      marginal_loss_avg(d) = marginal_loss_avg(d) / tmp
      NUE_vector(d) = NUE_vector(d) / tmp
      ! move along the counter for number of times the age class has been assessed too
      leaf_loss_possible(d) = sum(leaf_loss_possible((oldest_leaf-b):oldest_leaf))
      ! accumulate the over spill in the final age class
      canopy_age_vector(d) = tmp

      ! empty the previous days
      canopy_age_vector((oldest_leaf-b):oldest_leaf) = 0d0
      marginal_loss_avg((oldest_leaf-b):oldest_leaf) = 0d0
      NUE_vector((oldest_leaf-b):oldest_leaf) = 0d0
      leaf_loss_possible((oldest_leaf-b):oldest_leaf) = 0

    endif ! oldest_leaf+age_by_time > size(canopy_age_vector)

    ! oldest_leaf < canopy_maturation_lag, so we assume that time effects
    ! alone at at play
    d = oldest_leaf + age_by_time

    ! now we are ready to iterate the remaining age classes as normal
    ! i.e. starting a day back from the period we just adjusted...
    ! ...only apply the temperature enhanced aging to the mature canopy
    do a = oldest_leaf, 1, -1
      canopy_age_vector(a+age_by_time) = canopy_age_vector(a)
      marginal_loss_avg(a+age_by_time) = marginal_loss_avg(a)
      NUE_vector(a+age_by_time) = NUE_vector(a)
      leaf_loss_possible(a+age_by_time) = leaf_loss_possible(a)
    end do ; a = 1

    ! oldest leaf now at the end of the vector
    oldest_leaf = d

    ! finally the newest space must now be cleared for potential new growth
    canopy_age_vector(1:age_by_time) = 0d0
    marginal_loss_avg(1:age_by_time) = 0d0
    leaf_loss_possible(1:age_by_time) = 0
    ! overlay the default values for the maturing canopy
    NUE_vector(1:nint(canopy_maturation_lag)) = NUE_vector_mature(1:nint(canopy_maturation_lag))

    !
    ! Marginal return of leaf loss
    !

    ! initialise counters
    death_point = oldest_leaf ; oldest_point = oldest_leaf

    ! can't quantify leaf loss if there are no leaves...
    if (foliage > 0d0) then

      !
      ! lose leaves with negative marginal return
      !

      ! loop through the canopy age classes for which NUE is declining and determine whether they are paying for themselves
      do a = oldest_leaf, nint(canopy_maturation_lag), -1

        ! set current possible escape point (age class)
        death_point = a
        ! how much C is there to lose in the current age class?
        deltaC = canopy_age_vector(a)

        ! assess marginal return of there is LAI available
        if (deltaC > 0d0) then

          ! calculate Rm of the potentially lost leaves
          deltaRm = Rm_leaf * (deltaC/foliage)
          ! calculate leaf area of the leaf to be assessed
          lai = deltaC * SLA
          ! calculate mean canopy age of those considered for loss, weighted by the mass at each age point
          !canopy_age = (canopy_age_vector(a) * canopy_days(a)) / canopy_age_vector(a)
          canopy_age = canopy_days(a)
          ! update NUE as function of age
          !NUE = canopy_aggregate_NUE(a,a)
          NUE = NUE_vector(a)
          ! use proportional scaling of radiation, stomatal conductance and aerodynamics conductances
          tmp = lai / lai_save
          stomatal_conductance = gs_save * tmp ; aerodynamic_conductance = ga_save * tmp
          canopy_par_MJday = canopy_par_save * tmp
          if (stomatal_conductance > vsmall) then
            deltaGPP = acm_gpp(stomatal_conductance)
          else
            deltaGPP = 0d0
          endif
          ! Estimate marginal return of leaf loss vs cost of regrowth
          marginal_loss = deltaRm - deltaGPP! instantaneous costs of continued life
          marginal_loss = marginal_loss / deltaC ! scaled to per gC/m2
          ! pass to marginal_loss_avg vector
          marginal_loss_avg(a) = marginal_loss_avg(a) * (1d0 - leaf_mortality_period_1) &
                               + marginal_loss * leaf_mortality_period_1
          ! assess escape condition
          if (marginal_loss <= 0d0 .and. leaf_loss_possible(a) == 0) exit

        end if ! if there is leaf area in the current age class

      end do ! from oldest leaf back to NUE decline phases

      ! keep count of the number of times each age class has been assessed
      ! NOTE: must be before death_point is adjusted to consider the oldest_leaf to lose
      leaf_loss_possible(death_point:oldest_leaf) = leaf_loss_possible(death_point:oldest_leaf) + 1

      ! Need to re-initialise counters
      death_point = oldest_leaf ; oldest_point = oldest_leaf

      ! now loop back through to find the location of the oldest leaf
      ! which we have assess sufficient times and is marginally due to be lost
      do a = oldest_leaf, nint(canopy_maturation_lag), -1

        ! track for use outside of the loop
        oldest_point = a
        if (canopy_age_vector(a) > 0d0) then
          ! if the portion of the canopy we are in has not been checked enough or is not good to lose, abort loop
          if (leaf_loss_possible(a) < leaf_mortality_period .or. marginal_loss_avg(a) < 0d0) then
            ! assuming that our escape leaf is not the oldest leaf in the canopy,
            ! calculate the total loss. Otherwise we keep the whole canopy
            if (oldest_point < oldest_leaf) then
              ! estimate the total removal and convert into turnover fraction
              leaf_fall = sum(canopy_age_vector(death_point:oldest_leaf))
              if (leaf_fall > 0d0) then
                ! estimate daily rate equivalent
                leaf_fall = leaf_fall * deltat_1(current_step)
                ! remove the biomass from the canopy age vector
                canopy_age_vector(death_point:oldest_leaf) = 0d0
                ! and marginal return calculation
                marginal_loss_avg(death_point:oldest_leaf) = 0d0
                ! and counter
                leaf_loss_possible(death_point:oldest_leaf) = 0
                ! update the new oldest leaf age
                oldest_leaf = oldest_point
                ! cannot have oldest leaf at position less than 1
                if (oldest_leaf <= 0) oldest_leaf = 1
              endif ! leaf_fall > 0d0
            else
              ! we must be > oldest_leaf, i.e. none of the canopy
              ! should be turned over
              leaf_fall = 0d0
              death_point = oldest_leaf ; oldest_point = oldest_leaf
            end if ! death_point /= oldest_leaf

            ! now escape the loop
            exit

          end if ! leaf_loss_possible(a) < leaf_mortality_period .or. marginal_loss_avg(a) < 0d0
          ! track last known occupied age class, potentially the new
          ! oldest_leaf
          death_point = oldest_point
        end if ! canopy_age_vector(a) > 0d0

      end do ! from oldest leaf back to NUE decline phases

      ! restore original value back from memory
      lai = lai_save ; NUE = NUE_save ; canopy_age = canopy_age_save
      !canopy_lwrad_Wm2 = canopy_lw_save ; soil_lwrad_Wm2 = soil_lw_save
      canopy_swrad_MJday = canopy_sw_save ; canopy_par_MJday = canopy_par_save
      soil_swrad_MJday = soil_sw_save ; stomatal_conductance = gs_save
      aerodynamic_conductance = ga_save

    endif ! foliage > 0

    !
    ! Marginal return of leaf growth
    !

    ! If there is not labile available no growth can occur moreover we should not consider growing more leaves if
    ! there is a positive marginal return on losing leaf area.
    ! NOTE: that C shortage linked mortality was calculated earlier in the code
    ! should have been managed else where as mortality
    if (avail_labile > 0d0) then

      !
      ! Estimate potential C allocation to canopy
      !

      ! we are in an assending condition so labile turnover
      leaf_growth = pot_leaf_growth*Croot_labile_release_coef(current_step)*Cwood_hydraulic_limit
      ! calculate potential C allocation to leaves
      leaf_growth = avail_labile * (1d0-(1d0-leaf_growth)**deltat(current_step))*deltat_1(current_step)
      deltaC = leaf_growth
      ! C spent on growth
      deltaC = deltaC * deltat(current_step)
      ! if (deltaC > avail_labile) then
      !     tmp = (avail_labile / deltaC)
      !     leaf_growth = leaf_growth * tmp
      !     avail_labile = avail_labile * tmp
      ! endif

      !
      ! Estimate C losses from Rg and Rm
      !

      ! C to new growth
      tmp = deltaC * (one_Rg_fraction)
      ! calculate new Rm...
      if (Rm_leaf > 0d0) then
        ! if we have exisiting estimat of Rm_leaf then scale accordingly...
        deltaRm = Rm_leaf * ((foliage+tmp)/foliage)
      else
        !...if not then we will have to calculate the Rm cost directly
        deltaRm = Q10_adjustment * Rm_leaf_baseline * 2d-3 * umol_to_gC * seconds_per_day * (foliage+tmp)
      endif
      ! ...and its marginal return
      deltaRm = deltaRm - Rm_leaf

      !
      ! Estimate C update for the existing canopy at reference NUE
      ! (i.e. mean over lifespan)
      !

      ! overwrite current age specific NUE with the leaf lifetime NUE
      NUE = NUE_mean
      ! estimate stomatal conductance under the new NUE
      call acm_albedo_gc(abs(deltaWP),Rtot)
      ! and estimate the equivalent GPP for current LAI
      if (stomatal_conductance > vsmall) then
        old_GPP = acm_gpp(stomatal_conductance)
      else
        old_GPP = 0d0
      endif

      !
      ! Estimate C update for the updated canopy at reference NUE
      ! (i.e. mean over lifespan)
      !

      ! calculate new leaf area
      lai = (foliage+tmp) * SLA
      ! then add new
      canopy_age_vector(1) = tmp
!      ! calculate updated canopy states for new leaf area
!      call calculate_aerodynamic_conductance
!      call calculate_shortwave_balance ; call calculate_longwave_isothermal(leafT,maxt)
!      call acm_albedo_gc(abs(deltaWP),Rtot)
      ! calculate updated canopy states for new leaf area
      tmp = lai / lai_save
      aerodynamic_conductance = aerodynamic_conductance * tmp
      stomatal_conductance = stomatal_conductance * tmp
      call calculate_shortwave_balance !; call calculate_longwave_isothermal(leafT,maxt)
      if (lai_save < vsmall) then
        call calculate_aerodynamic_conductance
        call acm_albedo_gc(abs(deltaWP),Rtot)
      endif
     ! now estimate GPP with new LAI
      if (stomatal_conductance > vsmall) then
        alt_GPP = acm_gpp(stomatal_conductance)
      else
        alt_GPP = 0d0
      endif
      deltaGPP = alt_GPP - old_GPP
      ! is the marginal return for GPP (over the mean life of leaves)
      ! less than increase in maintenance respiration and C required to
      ! growth (Rg+tissue)?
      marginal_gain = ((deltaGPP-deltaRm)*leaf_life) - deltaC
      ! scale to per m2 leaf area, as the gradient requires comparable adjustment magnitudes
      marginal_gain = marginal_gain / deltaC

      ! accumulate marginal return information
      marginal_gain_avg = marginal_gain_avg * (1d0-leaf_growth_period_1) + &
      marginal_gain * (leaf_growth_period_1)

      ! restore original value back from memory
      lai = lai_save ; NUE = NUE_save ; canopy_age = canopy_age_save
      !canopy_lwrad_Wm2 = canopy_lw_save ; soil_lwrad_Wm2 = soil_lw_save
      canopy_swrad_MJday = canopy_sw_save ; canopy_par_MJday = canopy_par_save
      soil_swrad_MJday = soil_sw_save ; stomatal_conductance = gs_save
      aerodynamic_conductance = ga_save

    endif ! avail_labile > 0

    ! no positives of growing new leaves so don't
    if (marginal_gain_avg < 0d0) then
      ! Marginal suggest that we are not gaining leaves
      ! Therefore, we must clear our "new" allocation from the first age class
      leaf_growth = 0d0
      canopy_age_vector(1) = 0d0
      marginal_loss_avg(1) = 0d0
      leaf_loss_possible(1) = 0
    endif

  end subroutine calculate_leaf_dynamics
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_wood_root_growth(n,lab_to_roots,lab_to_wood &
                                       ,deltaWP,Rtot,current_gpp,Croot,Cwood  &
                                       ,root_growth,wood_growth)

    implicit none

    ! Premise of wood and root phenological controls

    ! Assumption 1:
    ! Based on plant physiology all cell expansion can only occur if there is
    ! sufficient water pressure available to drive the desired expansion.
    ! Moreover, as there is substantial evidence that shows wood and root growth
    ! do not follow the same phenology as leaves or GPP availability.
    ! Therefore, their phenological constrols should be separate from both that of
    ! the GSI model driving canopy phenology or GPP. Wood growth is limited by a
    ! logisitic temperature response assuming <5 % growth potential at 5oC and
    ! >95 % growth potential at 30 oC. Wood growth is also limited by a
    ! logistic response to water availability. When deltaWP (i.e. minleaf-wSWP)
    ! is less than -1 MPa wood growth is restricted to <5 % of potential.
    ! See review Fatichi et al (2013). Moving beyond phtosynthesis from carbon
    ! source to sink driven vegetation modelling. New Phytologist,
    ! https://doi.org/10.1111/nph.12614 for further details.

    ! As with wood, root phenology biologically speaking is independent of
    ! observed foliar phenological dynamics and GPP availabilty and thus has
    ! a separate phenology model. Similar to wood, a logistic temperature
    ! response is applied such that root growth is <5 % of potential at 0oC and
    ! >95 % of potential at 30oC. The different temperature minimua between
    ! wood and root growth is due to observed root growth when ever the soil
    ! is not frozen. We also assume that root growth is less sensitive to
    ! available hydraulic pressure, see assumption 3.

    ! Assumption 2:
    ! Actual biological theory suggests that roots support demands for resources
    ! made by the rest of the plant in this current model this is water only.
    ! Therefore there is an implicit assumption that roots should grow so long as
    ! growth is environmentally possible as growth leads to an improvement in
    ! C balance over their life time greater than their construction cost.

    ! Assumption 3:
    ! Determining when root growth should stop is poorly constrained.
    ! Similar to wood growth, here we assume root expansion is also dependent on water availability,
    ! but is less sensitive than wood. Root growth is assumed to stop when deltaWP approaches 0,
    ! determined by marginal return on root growth and temperature limits.

    ! arguments
    integer, intent(in) :: n
    double precision, intent(in) :: lab_to_roots,lab_to_wood &
                                   ,deltaWP,Rtot,current_gpp,Croot,Cwood
    double precision, intent(out) :: root_growth,wood_growth

    ! local variables
    double precision :: tmp, &
                        canopy_lw_save, canopy_sw_save, &
                        canopy_par_save, soil_lw_save, &
                        soil_sw_save, gs_save

    ! reset allocation to roots and wood
    root_growth = 0d0 ; wood_growth = 0d0

    ! save original values for re-allocation later
    !gs_save = stomatal_conductance

    ! Is it currently hydraulically possible for cell expansion (i.e. is soil
    ! water potential more negative than min leaf water potential).
    if ( avail_labile > 0d0 .and. deltaWP < 0d0 ) then

      ! Assume potential root growth is dependent on hydraulic and temperature conditions.
      ! Actual allocation is only allowed if the marginal return on GPP,
      ! averaged across the life span of the root is greater than the rNPP and Rg_root.

      ! Temperature limited turnover rate of labile -> roots
      root_growth = lab_to_roots*Croot_labile_release_coef(n)
      ! Estimate potential root allocation over time for potential root allocation
      root_growth = avail_labile*(1d0-(1d0-root_growth)**days_per_step)*days_per_step_1
!      ! C spent on growth
!      root_cost = tmp*days_per_step
!      ! C to new growth
!      tmp = root_cost * (one_Rg_fraction)
!      ! remainder is Rg cost
!      root_cost = root_cost - tmp
!      ! C spend on maintenance
!      deltaRm = Rm_root*((Croot+tmp)/Croot)
!      deltaRm = deltaRm - Rm_root
!      ! adjust to extra biomass (i.e. less Rg_root)
!      tmp = max(min_root,(Croot+tmp)*2)
!      ! estimate new (potential) Rtot
!      tmp = calc_pot_root_alloc_Rtot(abs(tmp))
!      ! calculate stomatal conductance of water
!      call acm_albedo_gc(abs(deltaWP),tmp)
!      if (lai > vsmall .and. stomatal_conductance > vsmall) then
!         tmp = acm_gpp(stomatal_conductance)
!      else
!         tmp = 0d0
!      end if
!      ! calculate marginal return on new root growth, scaled over life span
!      ! of new root.
!      deltaGPP = tmp-current_gpp
!      leaf_marginal = deltaGPP-deltaRm
!      ! if marginal return on GPP is less than growth and maintenance
!      ! costs of the life of the roots grow new roots
!!      if (((deltaGPP - deltaRm)*root_life) - root_cost < 0d0) root_growth = 0d0
!      if ((deltaGPP - deltaRm) < 0d0) root_growth = 0d0
!      ! can current NPP sustain the additional maintenance costs associated
!      ! with the new roots?
!      if (current_gpp - Rm_wood - Rm_leaf - Rm_root - deltaRm < 0d0) root_cost = 0d0

      ! determine wood growth based on temperature and hydraulic limits
      wood_growth = lab_to_wood*Cwood_labile_release_coef(n)*Cwood_hydraulic_limit
      ! Estimate potential root allocation over time for potential root
      ! allocation
      wood_growth = avail_labile*(1d0-(1d0-wood_growth)**days_per_step)*days_per_step_1
!      ! C spent on growth
!      wood_cost = wood_growth*days_per_step
!      ! C to new growth
!      tmp = wood_cost * (one_Rg_fraction)
!      ! remainder is Rg cost
!      wood_cost = wood_cost - tmp
!      ! C spend on maintenance
!      deltaRm = Rm_wood*((Cwood+tmp)/Cwood)
!      deltaRm = deltaRm - Rm_wood
      ! is current GPP less current maintenance costs sufficient to pay for Rg and new Rm?
!      if (current_gpp - Rm_wood - Rm_leaf - Rm_root - deltaRm - wood_cost < 0d0) wood_growth = 0d0

      ! estimate target woody C:N based on assumption that CN_wood increases
      ! logarithmically with increasing woody stock size.
!      CN_wood_target = 10d0**(log10(pars(15)) + log10(Cwood)*pars(25))

      ! cost of wood construction and maintenance not accounted for here due
      ! to no benefit being determined


      ! track labile reserves to ensure that fractional losses are applied
      ! sequencially in assumed order of importance (leaf->root->wood)

      ! root production (gC.m-2.day-1), limited by available labile sugars
      root_growth = min(avail_labile*days_per_step_1,root_growth)
      avail_labile = avail_labile - (root_growth*days_per_step)
      ! wood production (gC.m-2.day-1), limited by available labile sugars
      wood_growth = min(avail_labile*days_per_step_1,wood_growth)
      avail_labile = avail_labile - (wood_growth*days_per_step)

    endif ! grow root and wood?

    ! restore original values
    !stomatal_conductance = gs_save

    return

  end subroutine calculate_wood_root_growth
  !
  !------------------------------------------------------------------
  !
! subroutine estimate_mean_NUE
!
!    ! Estimate mean NUE over the life time of the canopy
!
!    implicit none
!
!    ! local variables
!    double precision :: NUE_zero_days, day_zero
!
!    ! first estimate the total (parameterised) period for which NUE > 0
!    day_zero = canopy_maturation_lag + canopy_optimum_period + canopy_zero_efficiency
!    if (leaf_life > day_zero) then
!
!        ! current mean leaf life span is greater than period where NUE > 0
!
!        ! estimate time period at which NUE = 0
!        NUE_zero_days = leaf_life - day_zero
!        ! estimate weighted average over between the different canopy phases.
!        ! Note that 0.5 is because the mean NUE over the maturaing and aging phases is half the NUE_optimum
!        ! as they are bound by zero. Whereas the midphase is at the optimum itself
!        NUE_mean = (NUE_optimum * 0.5d0 * canopy_maturation_lag) &
!                 + (NUE_optimum * 0.5d0 * canopy_zero_efficiency) &
!                 + (NUE_optimum * canopy_optimum_period)
!        NUE_mean = NUE_mean / (canopy_maturation_lag+canopy_optimum_period+canopy_zero_efficiency+NUE_zero_days)
!
!    else
!
!        ! current leaf life span is within the period where NUE > 0
!
!        ! estimate NUE at expected end of leaf life
!        NUE_mean = age_dependent_NUE(leaf_life,NUE_optimum, &
!                                     canopy_maturation_lag,canopy_optimum_period,canopy_zero_efficiency)
!
!        if (leaf_life > canopy_maturation_lag + canopy_optimum_period) then
!
!            ! current leaf life span ends within NUE decay phase
!
!            ! find average for the mature period and apply timing weighting
!            NUE_mean = (NUE_optimum + NUE_mean) * 0.5d0 * (leaf_life-canopy_maturation_lag-canopy_optimum_period)
!            ! weight maturating and optimum periods with the mature period
!            NUE_mean = ((NUE_optimum * 0.5d0 * canopy_maturation_lag) + (NUE_optimum * canopy_optimum_period) + NUE_mean)
!            NUE_mean = NUE_mean / leaf_life
!
!        else if (leaf_life > canopy_maturation_lag) then
!
!            ! current leaf life span end within optimum phase
!
!            ! find average for the mature period and apply timing weighting
!            NUE_mean = (NUE_optimum + NUE_mean) * 0.5d0 * (leaf_life-canopy_maturation_lag)
!            ! weight maturating period with the optimum period
!            NUE_mean = (NUE_optimum * 0.5d0 * canopy_maturation_lag) + NUE_mean
!            NUE_mean = NUE_mean / leaf_life
!
!        else
!
!            ! current leaf life span must end within maturation phase...
!
!            ! weight maturating period with the mature period
!            NUE_mean = (NUE_mean * 0.5d0 * canopy_maturation_lag) / leaf_life
!
!        endif ! leaf_life ends in maturation, optimum or decaying phases?
!
!    endif ! leaf_life > day_zero
!
! end subroutine estimate_mean_NUE
  !
  !------------------------------------------------------------------
  !
  !------------------------------------------------------------------
  ! Functions other than the primary ACM-GPP are stored
  ! below this line.
  !------------------------------------------------------------------
  !
  !------------------------------------------------------------------
  !
  pure function arrhenious( a , b , t )

    ! The equation is simply...                        !
    !    a * exp( b * ( t - 25.0 ) / ( t + 273.15 ) )  !
    ! However, precision in this routine matters as it !
    ! affects many others. To maximise precision, the  !
    ! calculations have been split.                    !

    implicit none

    ! arguments..
    double precision,intent(in) :: a , b , t
    double precision            :: arrhenious

    ! local variables..
    double precision :: denominator, numerator

    numerator   = t - 25d0
    denominator = t + freeze
    arrhenious  = a * exp( b * numerator / denominator )

    return

  end function arrhenious
  !
  !----------------------------------------------------------------------
  !
  double precision function opt_max_scaling( max_val , optimum , kurtosis , current )

    ! estimates a 0-1 scaling based on a skewed guassian distribution with a
    ! given optimum, maximum and kurtosis

    implicit none

    ! arguments..
    double precision,intent(in) :: max_val, optimum, kurtosis, current

    ! local variables..
    double precision :: dummy, dummy1

    if ( current >= max_val ) then
      opt_max_scaling = 0d0
    else
      dummy1 = max_val - optimum
      dummy     = ( max_val - current ) / dummy1
      dummy     = exp( log( dummy ) * kurtosis * dummy1 )
      opt_max_scaling = dummy * exp( kurtosis * ( current - optimum ) )
    end if

  end function opt_max_scaling
  !
  !------------------------------------------------------------------
  !
  double precision function root_resistance(root_mass,thickness)

    !
    ! Calculates root hydraulic resistance (MPa m2 s mmol-1) in a soil-root zone
    !

    implicit none

    ! arguments
    double precision :: root_mass, & ! root biomass in layer (gbiomass)
                        thickness    ! thickness of soil zone roots are in

    ! calculate root hydraulic resistance
    root_resistance = root_resist / (root_mass*thickness)

    ! return
    return

  end function root_resistance
  !
  !-----------------------------------------------------------------
  !
  double precision function soil_resistance(root_length,thickness,soilC)

    !
    ! Calculates the soil hydraulic resistance (MPa m2 s mmol-1) for a given
    ! soil-root zone
    !

    implicit none

    ! arguments
    double precision :: root_length, & ! root length in soil layer (m)
                          thickness, & ! thickness of soil layer (m)
                              soilC    ! soil conductivity m2.s-1.MPa-1

    ! local variables
    double precision :: rs, rs2

    ! calculate
    rs  = (root_length*pi)**(-0.5d0)
    rs2 = log( rs * root_radius_1 ) / (two_pi*root_length*thickness*soilC)
    ! soil water resistance
    soil_resistance = rs2*1d-9*mol_to_g_water

    ! return
    return

  end function soil_resistance
  !
  !------------------------------------------------------------------
  !
  double precision function plant_soil_flow(root_layer,root_length,root_mass &
                                           ,demand,root_reach_in,transpiration_resistance)

    !
    ! Calculate soil layer specific water flow form the soil to canopy (mmolH2O.m-2.s-1)
    ! Accounting for soil, root and plant resistance, and canopy demand
    !

    ! calculate and accumulate steady state water flux in mmol.m-2.s-1
    ! From the current soil layer given an amount of root within the soil layer.

    implicit none

    ! arguments
    integer, intent(in) :: root_layer
    double precision, intent(in) :: root_length, &
                                      root_mass, &
                                         demand, &
                                  root_reach_in, &
                       transpiration_resistance

    ! local arguments
    double precision :: soilR1, &
                        soilR2

    ! soil conductivity converted from m.s-1 -> m2.s-1.MPa-1 by head
    soilR1 = soil_resistance(root_length,root_reach_in,soil_conductivity(root_layer)*head_1)
    soilR2 = root_resistance(root_mass,root_reach_in)
    plant_soil_flow = demand/(transpiration_resistance + soilR1 + soilR2)

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
    double precision, intent(in) :: xin

    ! local variables..
    double precision ::soil_wp

!    ! calculate the soil water potential (kPa)..
!    soil_WP = -0.001 * potA( water_retention_pass ) * xin**potB( water_retention_pass )
!    water_retention_saxton_eqns = 1000. * soil_wp + 10.    ! 10 kPa represents air-entry swp
    ! calculate the soil water potential (kPa)..
    soil_wp = -potA( water_retention_pass ) * xin**potB( water_retention_pass )
    water_retention_saxton_eqns = soil_wp + 10d0    ! 10 kPa represents air-entry swp

    return

  end function water_retention_saxton_eqns
  !
  !------------------------------------------------------------------
  !
  double precision function age_dependent_NUE(age,optimum,pre_lag)

    !
    ! Estimate leaf age photosynthetic nitrogen use efficiency (gC/gN/m2leaf/day).
    ! Modified version of Xu et al., (2017), using here a combination of 3 (instead of 2) linear equations governing NUE.
    !
    ! Xu et al., (2017):
    ! Variations of leaf longevity in tropical moist forests predicted by a trait-driven carbon optimality model,
    ! Ecology Letters, doi: 10.1111/ele.12804

    implicit none

    ! arguments
    double precision, intent(in) :: age, optimum, pre_lag

    ! update NUE as function of age
    if (age < pre_lag) then
      ! canopy has not yet matured
      age_dependent_NUE = optimum * (1d0 - (pre_lag-age)/pre_lag)
    else
      ! canopy is at optimum
      age_dependent_NUE = optimum
    endif

    return

  end function age_dependent_NUE
  !
  !------------------------------------------------------------------
  !
! double precision function age_dependent_NUE(age,optimum,pre_lag,opt_lag,post_lag)
!
!   !
!   ! Estimate leaf age photosynthetic nitrogen use efficiency (gC/gN/m2leaf/day).
!   ! Modified version of Xu et al., (2017), using here a combination of 3 (instead of 2) linear equations governing NUE.
!   !
!   ! Xu et al., (2017):
!   ! Variations of leaf longevity in tropical moist forests predicted by a trait-driven carbon optimality model,
!   ! Ecology Letters, doi: 10.1111/ele.12804
!
!
!   implicit none
!
!   ! arguments
!   double precision, intent(in) :: age, optimum, pre_lag, opt_lag, post_lag
!
!   ! update NUE as function of age
!   if (age < pre_lag) then
!       ! canopy has not yet matured
!       age_dependent_NUE = optimum * (1d0 - (pre_lag-age)/pre_lag)
!   else if (age >= pre_lag .and. age <= pre_lag+opt_lag) then
!       ! canopy is at optimum
!       age_dependent_NUE = optimum
!   else if (age > pre_lag+opt_lag) then
!       ! canopy past optimum
!       age_dependent_NUE = optimum * (1d0 - (age-pre_lag-opt_lag)/post_lag)
!   else
!       ! we have missed something here...
!       print*,"Error in 'age_dependent_NUE'"
!       print*,"Current canopy age = ",age,"optimum NUE = ",optimum
!       print*,"pre_period = ",pre_lag,"opt_period = ",opt_lag,"post_period = ",post_lag
!       stop
!   endif
!
!   ! cannot be negative
!   age_dependent_NUE = max(0d0,age_dependent_NUE)
!
!   return
!
! end function age_dependent_NUE
  !
  !------------------------------------------------------------------
  !
  double precision function canopy_aggregate_NUE(first_leaf,last_leaf)

    !
    ! Estimate photosynthetic nitrogen use efficiency (gC/gN/m2leaf/day)
    ! aggregated across the canopy age distribution

    !
    ! Xu et al., (2017):
    ! Variations of leaf longevity in tropical moist forests predicted by a trait-driven carbon optimality model,
    ! Ecology Letters, doi: 10.1111/ele.12804

    implicit none

    ! arguments
    integer, intent(in) :: first_leaf, last_leaf

    ! weight the age specific NUE values by the amount of foliage biomass in each class
    canopy_aggregate_NUE = sum(canopy_age_vector(first_leaf:last_leaf)*NUE_vector(first_leaf:last_leaf))
    canopy_aggregate_NUE = canopy_aggregate_NUE / sum(canopy_age_vector(first_leaf:last_leaf))

    return

  end function canopy_aggregate_NUE
  !
  !------------------------------------------------------------------
  !
  double precision function calc_pot_root_alloc_Rtot(potential_root_biomass)

    !
    ! Description
    !

    implicit none

    ! declare arguments
    double precision,intent(in) :: potential_root_biomass ! potential root biomass (g.m-2)

    ! declare local variables
    integer :: i
    double precision, dimension(nos_root_layers) :: water_flux_local &
                                                       ,demand_local &
                                                       ,root_mass    &
                                                       ,root_length  &
                                                       ,ratio
    double precision, dimension(nos_soil_layers+1) :: soil_waterfrac_save, soil_conductivity_save,layer_thickness_save
    double precision :: bonus,transpiration_resistance,root_reach_local,root_reach_local_local &
                       ,depth_change,water_change,root_depth_50
    double precision, parameter :: root_depth_frac_50 = 0.25d0 ! fractional soil depth above which 50 %
                                                               ! of the root mass is assumed to be located

    ! estimate rooting depth with potential root growth
    root_reach_local = max_depth * potential_root_biomass / (root_k + potential_root_biomass)
    ratio = 0d0 ; ratio(1) = 1d0

    ! save soil water information
    soil_waterfrac_save = soil_waterfrac
    soil_conductivity_save = soil_conductivity
    layer_thickness_save = layer_thickness

    !!!!!!!!!!
    ! Update soil layer thickness for marginal return calculation
    !!!!!!!!!!

    depth_change = 0d0 ; water_change = 0d0
    ! if roots extent down into the bucket
    if (root_reach_local > (top_soil_depth+mid_soil_depth) .or. previous_depth > (top_soil_depth+mid_soil_depth)) then
      ! how much has root depth extended since last step?
      depth_change = root_reach_local - previous_depth

      ! if there has been an increase
      if (depth_change > 0.01d0 .and. root_reach_local > sum(layer_thickness(1:2))+min_layer) then

        ! determine how much water is within the new volume of soil
        water_change = soil_waterfrac(nos_soil_layers) * depth_change
        ! now assign that new volume of water to the deep rooting layer
        soil_waterfrac(nos_root_layers) = ((soil_waterfrac(nos_root_layers) * layer_thickness(nos_root_layers)) &
                                        + water_change) / (layer_thickness(nos_root_layers)+depth_change)
        ! explicitly update the soil profile if there has been rooting depth
        ! changes
        layer_thickness(1) = top_soil_depth ; layer_thickness(2) = mid_soil_depth
        layer_thickness(3) = max(min_layer,root_reach_local-sum(layer_thickness(1:2)))
        layer_thickness(4) = max_depth - sum(layer_thickness(1:3))

      elseif (depth_change < -0.01d0 .and. root_reach_local > layer_thickness(1)+min_layer) then

        ! determine how much water is lost from the old volume of soil
        water_change = soil_waterfrac(nos_root_layers) * abs(depth_change)
        ! now assign that new volume of water to the deep rooting layer
        soil_waterfrac(nos_soil_layers) = ((soil_waterfrac(nos_soil_layers) * layer_thickness(nos_soil_layers)) &
                                        + water_change) / (layer_thickness(nos_soil_layers)+abs(depth_change))

        ! explicitly update the soil profile if there has been rooting depth
        ! changes
        layer_thickness(1) = top_soil_depth ; layer_thickness(2) = mid_soil_depth
        layer_thickness(3) = max(min_layer,root_reach_local-sum(layer_thickness(1:2)))
        layer_thickness(4) = max_depth - sum(layer_thickness(1:3))

      else

        ! we don't want to do anything, just recycle the previous depth

      end if ! depth change

    end if ! root reach beyond top layer

    ! seperately calculate the soil conductivity as this applies to each layer
    do i = 1, nos_soil_layers
      call calculate_soil_conductivity(i,soil_waterfrac(i),soil_conductivity(i))
    end do ! soil layers

    ! estimate water flux based on soil and root hydraulic resistances with potential growth.
    ! See subroutine calculate_Rtot for further details

    ! calculate the plant hydraulic resistance component
!    transpiration_resistance = (gplant * lai)**(-1d0)
    transpiration_resistance = canopy_height / (gplant * lai)

    !!!!!!!!!!!
    ! Calculate root profile
    !!!!!!!!!!!

    ! top 25 % of root profile
    root_depth_50 = root_reach_local * root_depth_frac_50
    if (root_depth_50 <= layer_thickness(1)) then
      ! Greater than 50 % of the fine root biomass can be found in the top
      ! soil layer

      ! Start by assigning all 50 % of root biomass to the top soil layer
      root_mass(1) = root_biomass * 0.5d0
      ! Then quantify how much additional root is found in the top soil layer
      ! assuming that the top 25 % depth is found somewhere within the top
      ! layer
      bonus = (root_biomass-root_mass(1)) &
            * (layer_thickness(1)-root_depth_50) / (root_reach_local - root_depth_50)
      root_mass(1) = root_mass(1) + bonus
      ! partition the remaining root biomass between the seconds and third
      ! soil layers
      if (root_reach_local > sum(layer_thickness(1:2))) then
        root_mass(2) = (root_biomass - root_mass(1)) &
                     * (layer_thickness(2)/(root_reach_local-layer_thickness(1)))
        root_mass(3) = root_biomass - sum(root_mass(1:2))
      else
        root_mass(2) = root_biomass - root_mass(1)
      endif
    else if (root_depth_50 > layer_thickness(1) .and. root_depth_50 <= sum(layer_thickness(1:2))) then
      ! Greater than 50 % of fine root biomass found in the top two soil
      ! layers. We will divide the root biomass uniformly based on volume,
      ! plus bonus for the second layer (as done above)
      root_mass(1) = root_biomass * 0.5d0 * (layer_thickness(1)/root_depth_50)
      root_mass(2) = root_biomass * 0.5d0 * ((root_depth_50-layer_thickness(1))/root_depth_50)
      ! determine bonus for the seconds layer
      bonus = (root_biomass-sum(root_mass(1:2))) &
            * ((sum(layer_thickness(1:2))-root_depth_50)/(root_reach_local-root_depth_50))
      root_mass(2) = root_mass(2) + bonus
      root_mass(3) = root_biomass - sum(root_mass(1:2))
    else
      ! Greater than 50 % of fine root biomass stock spans across all three
      ! layers
      root_mass(1) = root_biomass * 0.5d0 * (layer_thickness(1)/root_depth_50)
      root_mass(2) = root_biomass * 0.5d0 * (layer_thickness(2)/root_depth_50)
      root_mass(3) = root_biomass - sum(root_mass(1:2))
    endif
    ! now convert root mass into lengths
    root_length = root_mass * root_mass_length_coef_1
!    root_length = root_mass / (root_density * root_cross_sec_area)

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
    demand_local = abs(minlwp-SWP(1:nos_root_layers))+head*canopy_height

    do i = 1, nos_root_layers
      if (root_mass(i) > 0d0) then
        ! if there is root then there is a water flux potential...
        root_reach_local_local = min(root_reach_local,layer_thickness(i))
        ! calculate and accumulate steady state water flux in mmol.m-2.s-1
        water_flux_local(i) = plant_soil_flow(i,root_length(i),root_mass(i) &
                                             ,demand_local(i),root_reach_local_local,transpiration_resistance)
      else
        ! ...if there is not then we wont have any below...
        exit
      end if ! root present in current layer?
    end do ! nos_root_layers    ratio = layer_thickness(1:nos_root_layers)/sum(layer_thickness(1:nos_root_layers))

    ! if freezing then assume soil surface is frozen
    if (meant < 1d0) then
      water_flux_local(1) = 0d0
      ratio(1) = 0d0
      ratio(2:nos_root_layers) = layer_thickness(2:nos_root_layers) / sum(layer_thickness(2:nos_root_layers))
    else
      ratio = layer_thickness(1:nos_root_layers)/sum(layer_thickness(1:nos_root_layers))
    endif

    ! WARNING: should probably have updated the wSWP here as well...do this
    ! later I thinks...

    ! determine effective resistance
    calc_pot_root_alloc_Rtot = sum(demand) / sum(water_flux_local)

    ! return layer_thickness and soil_waterfrac back to
    ! orginal values
    layer_thickness = layer_thickness_save
    soil_waterfrac = soil_waterfrac_save
    soil_conductivity = soil_conductivity_save

    return

  end function calc_pot_root_alloc_Rtot
  !
  !------------------------------------------------------------------
  !
  double precision function linear_model_gradient(x,y,interval)

    ! Function to calculate the gradient of a linear model for a given dependent
    ! variable (y) based on predictive variable (x). The typical use of this
    ! function will in fact be to assume that x is time.

    implicit none

    ! declare input variables
    integer :: interval
    double precision, dimension(interval) :: x,y

    ! declare local variables
    double precision :: sum_x, sum_y, sumsq_x,sum_product_xy

    ! calculate the sum of x
    sum_x = sum(x)
    ! calculate the sum of y
    sum_y = sum(y)
    ! calculate the sum of squares of x
    sumsq_x = sum(x*x)
    ! calculate the sum of the product of xy
    sum_product_xy = sum(x*y)
    ! calculate the gradient
    linear_model_gradient = ( (dble(interval)*sum_product_xy) - (sum_x*sum_y) ) &
                          / ( (dble(interval)*sumsq_x) - (sum_x*sum_x) )

!    ! for future reference here is how to calculate the intercept
!    intercept = ( (sum_y*sumsq_x) - (sum_x*sum_product_xy) ) &
!              / ( (dble(interval)*sumsq_x) - (sum_x*sum_x) )

    ! don't forget to return to the user
    return

  end function linear_model_gradient
  !
  !--------------------------------------------------------------------------
  !
  double precision function Rm_reich_Q10(air_temperature)

    ! Calculate Q10 temperature adjustment used in estimation of the
    ! Maintenance respiration (umolC.m-2.s-1) calculated based on modified
    ! version of the Reich et al (2008) calculation.

    ! arguments
    double precision, intent(in) :: air_temperature ! input temperature of metabolising tissue (oC)

    ! local variables
    double precision, parameter :: Q10 = 2d0,  & ! Q10 response of temperature (baseline = 20oC) ;INITIAL VALUE == 2
    ! Mahecha, et al. (2010) Global Convergence in the Temperature Sensitivity of
    ! Respiration at Ecosystem Level. Science 329 , 838 (2010);
    ! DOI: 10.1126/science.1189587. value reported as 1.4
    Q10_baseline = 20d0      ! Baseline temperature for Q10 ;INITIAL VALUE == 20;

    ! calculate instantaneous Q10 temperature response
    Rm_reich_Q10 = Q10**((air_temperature-Q10_baseline)*0.1d0)

    ! explicit return command
    return

  end function Rm_reich_Q10
  !
  !--------------------------------------------------------------------------
  !
  double precision function Rm_reich_N(CN_pool &
                                      ,N_exponential_response  &
                                      ,N_scaler_intercept)

    ! Calculate the nitrgen response on maintenance respiration (nmolC.g-1.s-1)
    ! calculated based on modified version of the Reich et al (2008) calculation.
    ! Note the output here is invarient for given CN ratio

    ! NOTE: Rm_tissue =  Rm_reich_Q10 * Rm_reich_N * C_POOL * 2 * 0.001
    !       umolC/m2/s = dimensionless * nmolC/g/s * gC/m * (correct g->gC) * (nmolC->umolC)

    ! arguments
    double precision, intent(in) ::        CN_pool, & ! C:N ratio for current pool (gC/gN)
                            N_exponential_response, & ! N exponential response coefficient (1.277/1.430)
                                N_scaler_intercept    ! N scaler (baseline) (0.915 / 1.079)

    ! local variables
    double precision, parameter :: N_g_to_mmol = 71.42857 ! i.e. (1d0/14d0)*1d3 where 14 = atomic weight of N
    double precision :: Nconc ! Nconc =mmol g-1

    ! calculate leaf maintenance respiration (nmolC.g-1.s-1)
    ! NOTE: that the coefficients in Reich et al., 2008 were calculated from
    ! log10 linearised version of the model, thus N_scaler is already in log10()
    ! scale. To remove the need of applying log10(Nconc) and 10**Rm_reich the
    ! scaler is reverted instead to the correct scale for the exponential form
    ! of the equations.

    ! calculate N concentration per g biomass.
    ! A function of C:N
    Nconc = ((CN_pool*2d0)**(-1d0)) * N_g_to_mmol

    ! leaf maintenance respiration (nmolC.g-1.s-1) at 20 oC
    Rm_reich_N = (10d0**N_scaler_intercept) * Nconc ** N_exponential_response

    ! explicit return command
    return

  end function Rm_reich_N
  !
  !------------------------------------------------------------------
  !
  !
  !------------------------------------------------------------------
  ! Generic mathematical functions such as bisection and intergrator proceedures
  ! are stored below here
  !------------------------------------------------------------------
  !
  !
  !------------------------------------------------------------------
  !
  double precision function zbrent( called_from , func , x1 , x2 , tol )

    ! This is a bisection routine. When ZBRENT is called, we provide a    !
    !  reference to a particular function and also two values which bound !
    !  the arguments for the function of interest. ZBRENT finds a root of !
    !  the function (i.e. the point where the function equals zero), that !
    !  lies between the two bounds.                                       !
    ! For a full description see Press et al. (1986).                     !

    implicit none

    ! arguments..
    character(len=*),intent(in) :: called_from    ! name of procedure calling (used to pass through for errors)
    double precision,intent(in) :: tol, x1, x2

    ! Interfaces are the correct way to pass procedures as arguments.
    interface
      double precision function func( xval )
        double precision ,intent(in) :: xval
      end function func
    end interface

    ! local variables..
    integer            :: iter
    integer, parameter :: ITMAX = 20
    double precision   :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    double precision, parameter :: EPS = 3d-8

    ! calculations...
    a  = x1
    b  = x2
    fa = func( a )
    fb = func( b )

    ! Check that we haven't (by fluke) already started with the root..
    if ( fa .eq. 0d0 ) then
      zbrent = a
      return
    elseif ( fb .eq. 0d0 ) then
      zbrent = b
      return
    end if
    ! Ensure the supplied x-values give y-values that lie either
    ! side of the root and if not flag an error message...
    if ( sign(1d0,fa) .eq. sign(1d0,fb) ) then
      fa = func( a )
      fb = func( b )
      ! tell me otherwise what is going on
      !       print*,"Supplied values must bracket the root of the function.",new_line('x'),  &
      !         "     ","You supplied x1:",x1,new_line('x'),                     &
      !         "     "," and x2:",x2,new_line('x'),                             &
      !         "     "," which give function values of fa :",fa,new_line('x'),  &
      !         "     "," and fb:",fb," .",new_line('x'),                        &
      !         " zbrent was called by: ",trim(called_from)
      !       fa = func( a )
      !       fb = func( b )
    end if
    c = b
    fc = fb

    do iter = 1 , ITMAX

      ! If the new value (f(c)) doesn't bracket
      ! the root with f(b) then adjust it..
      if ( sign(1d0,fb) .eq. sign(1d0,fc) ) then
        c  = a
        fc = fa
        d  = b - a
        e  = d
      end if
      if ( abs(fc) .lt. abs(fb) ) then
        a  = b
        b  = c
        c  = a
        fa = fb
        fb = fc
        fc = fa
      end if
      tol1 = 2d0 * EPS * abs(b) + 0.5d0 * tol
      xm   = 0.5d0 * ( c - b )
      if ( ( abs(xm) .le. tol1 ) .or. ( fb .eq. 0d0 ) ) then
        zbrent = b
        return
      end if
      if ( ( abs(e) .ge. tol1 ) .and. ( abs(fa) .gt. abs(fb) ) ) then
        s = fb / fa
        if ( a .eq. c ) then
          p = 2d0 * xm * s
          q = 1d0 - s
        else
          q = fa / fc
          r = fb / fc
          p = s * ( 2d0 * xm * q * ( q - r ) - ( b - a ) * ( r - 1d0 ) )
          q = ( q - 1d0 ) * ( r - 1d0 ) * ( s - 1d0 )
        end if
        if ( p .gt. 0d0 ) q = -q
        p = abs( p )
        if ( (2d0*p) .lt. min( 3d0*xm*q-abs(tol1*q) , abs(e*q) ) ) then
          e = d
          d = p / q
        else
          d = xm
          e = d
        end if
      else
        d = xm
        e = d
      end if
      a  = b
      fa = fb
      if ( abs(d) .gt. tol1 ) then
        b = b + d
      else
        b = b + sign( tol1 , xm )
      end if
      fb = func(b)
    enddo

    !    print*,"zbrent has exceeded maximum iterations",new_line('x'),&
    !           "zbrent was called by: ",trim(called_from)

    zbrent = b

  end function zbrent
  !
  !------------------------------------------------------------------
  !
!
!--------------------------------------------------------------------
!
end module CARBON_MODEL_MOD
