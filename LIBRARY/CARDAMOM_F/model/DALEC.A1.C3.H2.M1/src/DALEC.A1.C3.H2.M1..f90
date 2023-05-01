
module CARBON_MODEL_MOD

  implicit none

  !!!!!!!!!!!
  ! Authorship contributions
  !
  ! This code contains a variant of the Data Assimilation Linked ECosystem (DALEC) model.
  ! This version of DALEC is derived from the following primary references:
  ! Sus et al., (2010), https://doi.org/10.1016/j.agee.2010.06.012.
  ! Bloom & Williams (2015), https://doi.org/10.5194/bg-12-1299-2015.
  ! Smallman et al., (2017), https://doi.org/10.1002/2016JG003520.
  ! Smallman & Williams (2019) https://doi.org/10.5194/gmd-12-2227-2019.
  ! This code is based on that created by A. A. Bloom (UoE, now at JPL, USA) and
  ! Oliver Sus (UoE, now at EUMETSAT, Darmstadt).
  ! Subsequent modifications by:
  ! T. L. Smallman (University of Edinburgh, t.l.smallman@ed.ac.uk)
  ! See function / subroutine specific comments for exceptions and contributors
  !!!!!!!!!!!

  ! make all private
  private

  ! explicit publics
  public :: CARBON_MODEL                  &
           ,vsmall                        &
           ,arrhenious                    &
           ,acm_gpp_stage_1               &
           ,acm_gpp_stage_2               &
           ,calculate_transpiration       &
           ,calculate_wetcanopy_evaporation &
           ,calculate_soil_evaporation    &
           ,calculate_radiation_balance   &
           ,calculate_stomatal_conductance&
           ,meteorological_constants      &
           ,update_soil_initial_conditions&
           ,calculate_daylength           &
           ,opt_max_scaling               &
           ,freeze                        &
           ,co2comp_saturation            &
           ,co2comp_half_sat_conc         &
           ,kc_saturation                 &
           ,kc_half_sat_conc              &
           ,calculate_update_soil_water   &
           ,calculate_Rtot                &
           ,calculate_aerodynamic_conductance &
           ,saxton_parameters             &
           ,initialise_soils              &
           ,linear_model_gradient         &
           ,seconds_per_day               &
           ,seconds_per_step              &
           ,fine_root_biomass             &
           ,root_biomass                  &
           ,root_reach                    &
           ,min_root                      &
           ,max_depth                     &
           ,root_k                        &
           ,top_soil_depth                &
           ,previous_depth                &
           ,nos_root_layers               &
           ,wSWP                          &
           ,rSWP                          &
           ,cica_time                     &
           ,SWP                           &
           ,SWP_initial                   &
           ,deltat_1                      &
           ,total_water_flux              &
           ,water_flux_mmolH2Om2s         &
           ,layer_thickness               &
           ,waterchange                   &
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
           ,wSWP_time                     &
           ,rSWP_time                     &
           ,gs_demand_supply_ratio        &
           ,gs_total_canopy               &
           ,gb_total_canopy               &
           ,canopy_par_MJday              &
           ,canopy_par_MJday_time         &
           ,soil_frac_clay                &
           ,soil_frac_sand                &
           ,nos_soil_layers               &
           ,meant                         &
           ,convert_ms1_mol_1             &
           ,aerodynamic_conductance       &
           ,stomatal_conductance          &
           ,potential_conductance         &
           ,iWUE                          &
           ,avN                           &
           ,NUE                           &
           ,ceff                          &
           ,pn_max_temp                   &
           ,pn_opt_temp                   &
           ,pn_kurtosis                   &
           ,e0                            &
           ,co2_half_sat                  &
           ,co2_comp_point                &
           ,minlwp                        &
           ,mint                          &
           ,maxt                          &
           ,leafT                         &
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
           ,dayl_seconds                  &
           ,dayl_seconds_1                &
           ,dayl_hours                    &
           ,dayl_hours_fraction           &
           ,snow_storage                  &
           ,canopy_storage                &
           ,intercepted_rainfall          &
           ,rainfall_time                  &
           ,resp_rate_temp_coeff           &
           ,ts_length                      &
           ,dim_1,dim_2                    &
           ,nos_trees                      &
           ,nos_inputs                     &
           ,leftDaughter                   &
           ,rightDaughter                  &
           ,nodestatus                     &
           ,xbestsplit                     &
           ,nodepred                       &
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
  logical :: do_iWUE = .true. ! Use iWUE or WUE for stomatal optimisation
  double precision, parameter :: vsmall = tiny(0d0)*1d3 & ! *1d3 to add a little breathing room
                                ,vlarge = huge(0d0)

  integer, parameter :: nos_root_layers = 2, nos_soil_layers = nos_root_layers + 1
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
                       gs_H2O_CO2_1 = 0.6074378d0,  & ! gs_H2O_CO2 ** (-1d0)
                  gs_H2Ommol_CO2mol = 0.001646259d0,& ! gs_H2O_CO2 * 1d-3
              gs_H2Ommol_CO2mol_day = 142.2368d0,   & ! The ratio of H20:CO2 diffusion for gs, including seconds per day correction
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
                      kc_saturation = 310d0,        & ! CO2 half saturation, saturation value
                   kc_half_sat_conc = 23.956d0,     & ! CO2 half sat, half sat
                 co2comp_saturation = 36.5d0,       & ! CO2 compensation point, saturation
              co2comp_half_sat_conc = 9.46d0,       & ! CO2 comp point, half sat
                                                      ! Each of these are temperature sensitivty
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
                           min_wind = 0.2d0,        & ! minimum wind speed at canopy top
                       min_drythick = 0.001d0,      & ! minimum dry thickness depth (m)
                          min_layer = 0.03d0,       & ! minimum thickness of the third rooting layer (m)
                        soil_roughl = 0.05d0,       & ! soil roughness length (m)
                     top_soil_depth = 0.30d0,       & ! thickness of the top soil layer (m)
                           min_root = 5d0,          & ! minimum root biomass (gBiomass.m-2)
                            min_lai = 0.1d0,        & ! minimum LAI assumed for aerodynamic conductance calculations (m2/m2)
                        min_storage = 0.2d0           ! minimum canopy water (surface) storage (mm)

  ! timing parameters
  double precision, parameter :: &
                   seconds_per_hour = 3600d0,         & ! Number of seconds per hour
                    seconds_per_day = 86400d0,        & ! Number of seconds per day
                  seconds_per_day_1 = 1.157407d-05      ! Inverse of seconds per day

  ! ACM-GPP-ET parameters
  double precision, parameter :: &
                   pn_max_temp = 6.842942d+01,  & ! Maximum daily max temperature for photosynthesis (oC)
                   pn_opt_temp = 3.155960d+01,  & ! Optimum daily max temperature for photosynthesis (oC)
                   pn_kurtosis = 1.889026d-01,  & ! Kurtosis of photosynthesis temperature response
                            e0 = 3.661204d+00,  & ! Quantum yield gC/MJ/m2/day PAR
                minlwp_default =-1.808224d+00,  & ! minimum leaf water potential (MPa)
      soil_iso_to_net_coef_LAI =-2.717467d+00,  & ! Coefficient relating soil isothermal net radiation to net.
                          iWUE = 6.431150d-06,  & ! Intrinsic water use efficiency (gC/m2leaf/day/mmolH2Ogs)
         soil_swrad_absorption = 9.989852d-01,  & ! Fraction of SW rad absorbed by soil
         max_lai_lwrad_release = 9.516639d-01,  & ! 1-Max fraction of LW emitted from canopy to be released
        lai_half_lwrad_release = 4.693329d+00,  & ! LAI at which LW emitted from canopy to be released at 50 %
       soil_iso_to_net_coef_SW =-3.500964d-02,  & ! Coefficient relating soil isothermal net radiation to net.
         soil_iso_to_net_const = 3.455772d+00,  & ! Constant relating soil isothermal net radiation to net
           max_par_transmitted = 1.628077d-01,  & ! Max fraction of canopy incident PAR transmitted to soil
           max_nir_transmitted = 2.793660d-01,  & ! Max fraction of canopy incident NIR transmitted to soil
             max_par_reflected = 1.629133d-01,  & ! Max fraction of canopy incident PAR reflected to sky
             max_nir_reflected = 4.284365d-01,  & ! Max fraction of canopy incident NIR reflected to sky
     canopy_iso_to_net_coef_SW = 1.480105d-02,  & ! Coefficient relating SW to the adjustment between isothermal and net LW
       canopy_iso_to_net_const = 3.753067d-03,  & ! Constant relating canopy isothermal net radiation to net
    canopy_iso_to_net_coef_LAI = 2.455582d+00     ! Coefficient relating LAI to the adjustment between isothermal and net LW

  double precision :: minlwp = minlwp_default

  !!!!!!!!!
  ! Module level variables
  !!!!!!!!!

  ! hydraulic model variables
  integer :: water_retention_pass, soil_layer
  double precision, dimension(nos_soil_layers) :: soil_frac_clay,soil_frac_sand ! clay and soil fractions of soil
  double precision, dimension(nos_root_layers) :: uptake_fraction, & ! fraction of water uptake from each root layer
                                                           demand, & ! maximum potential canopy hydraulic demand
                                            water_flux_mmolH2Om2s    ! potential transpiration flux (mmolH2O.m-2.s-1)
  double precision, dimension(nos_soil_layers+1) :: SWP, & ! soil water potential (MPa)
                                            SWP_initial, &
                                      soil_conductivity, & ! soil conductivity
                                            waterchange, & ! net water change by specific soil layers (m)
                                         field_capacity, & ! soil field capacity (m3.m-3)
                                 field_capacity_initial, &
                                         soil_waterfrac, & ! soil water content (m3.m-3)
                                 soil_waterfrac_initial, &
                                               porosity, & ! soil layer porosity, (fraction)
                                       porosity_initial, &
                                        layer_thickness, & ! thickness of soil layers (m)
                        cond1, cond2, cond3, potA, potB    ! Saxton equation values

  double precision :: root_reach, root_biomass, &
                             fine_root_biomass, & ! root depth, coarse+fine, and fine root biomass
                              total_water_flux, & ! potential transpiration flux (kgH2O.m-2.day-1)
                                      drythick, & ! estimate of the thickness of the dry layer at soil surface (m)
                                          wSWP, & ! soil water potential weighted by canopy supply (MPa)
                                          rSWP, & ! soil water potential weighted by root presence (MPa)
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
                                          leafT, & ! canopy temperature (oC)
                               mean_annual_temp, &
                             canopy_swrad_MJday, & ! canopy_absorbed shortwave radiation (MJ.m-2.day-1)
                               canopy_par_MJday, & ! canopy_absorbed PAR radiation (MJ.m-2.day-1)
                               soil_swrad_MJday, & ! soil absorbed shortwave radiation (MJ.m-2.day-1)
                               canopy_lwrad_Wm2, & ! canopy absorbed longwave radiation (W.m-2)
                                 soil_lwrad_Wm2, & ! soil absorbed longwave radiation (W.m-2)
                                  sky_lwrad_Wm2, & ! sky absorbed longwave radiation (W.m-2)
                          stomatal_conductance, & ! stomatal conductance (mmolH2O.m-2ground.s-1)
                         potential_conductance, & ! potential stomatal conductance (mmolH2O.m-2ground.s-1)
                           minimum_conductance, & ! potential stomatal conductance (mmolH2O.m-2ground.s-1)
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
  double precision ::   delta_gs, & ! day length corrected gs increment mmolH2O/m2/dayl
                            ceff, & ! canopy efficency, ceff = avN*NUE
                             avN, & ! average foliar N (gN/m2)
                       iWUE_step, & ! Intrinsic water use efficiency for that day (gC/m2leaf/dayl/mmolH2Ogs)
                             NUE, & ! Photosynthetic nitrogen use efficiency at optimum temperature (oC)
                                    ! ,unlimited by CO2, light and photoperiod (gC/gN/m2leaf/day
metabolic_limited_photosynthesis, & ! temperature, leaf area and foliar N limiterd photosynthesis (gC/m2/day)
    light_limited_photosynthesis, & ! light limited photosynthesis (gC/m2/day)
                              ci, & ! Internal CO2 concentration (ppm)
                          gb_mol, & ! Canopy boundary layer conductance (molCO2/m2/day)
                        rb_mol_1, & ! Canopy boundary layer resistance (day/m2/molCO2)
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
                            seconds_per_step, & !
                               days_per_step, & !
                             days_per_step_1, & !
                          mean_days_per_step, &
                                dayl_seconds, & ! day length in seconds
                              dayl_seconds_1, &
                         dayl_hours_fraction, &
                                  dayl_hours    ! day length in hours

  double precision, dimension(:), allocatable :: deltat_1, & ! inverse of decimal days
                                               meant_time, &
                                  airt_zero_fraction_time, &
                                          daylength_hours, &
                                        daylength_seconds, &
                                      daylength_seconds_1, &
                                            rainfall_time, &
                                   gs_demand_supply_ratio, & ! actual:potential stomatal conductance
                                          gs_total_canopy, & ! stomatal conductance (mmolH2O/m2ground/day)
                                          gb_total_canopy, & ! boundary conductance (mmolH2O/m2ground/day)
                                    canopy_par_MJday_time, & ! Absorbed PAR by canopy (MJ/m2ground/day)
                                                cica_time, & ! Internal vs ambient CO2 concentrations
                                                rSWP_time, &
                                                wSWP_time    ! Soil water potential weighted by root supply of water

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

  save

  contains
  !
  !--------------------------------------------------------------------
  !
  subroutine CARBON_MODEL(start,finish,met,pars,deltat,nodays,lat,lai_out &
                         ,NEE_out,FLUXES,POOLS,pft,nopars,nomet,nopools,nofluxes &
                         ,GPP_out,stock_seed_labile,DS_shoot,DS_root,fol_frac    &
                         ,stem_frac,root_frac,DS_LRLV,LRLV,DS_LRRT,LRRT)

    !
    ! The Data Assimilation Linked Ecosystem Carbon - CROP - BUCKET (DALEC_CROP_BUCKET) model.
    ! modified from Sus et al., (2010)
    !
    ! The Aggregated Canopy Model for Gross Primary Productivity and Evapotranspiration (ACM-GPP-ET)
    ! simulates coupled photosynthesis-transpiration (via stomata), soil and intercepted canopy evaporation and
    ! soil water balance (4 layers).
    !
    ! This version was coded by T. Luke Smallman (t.l.smallman@ed.ac.uk)
    ! Version 1: 15/07/2014
    ! Version 2: 15/11/2018 - Addition of the BUCKET model via ACM2 to include the water cycle

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
    ts_length = ((sum(deltat)/dble(nodays)) * seconds_per_day) / seconds_per_hour
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
        POOLS(1,1) = stock_labile
        POOLS(1,2) = stock_foliage
        POOLS(1,3) = stock_roots
        POOLS(1,4) = stock_stem
        POOLS(1,5) = stock_litter
        POOLS(1,6) = stock_soilOrgMatter
        POOLS(1,7) = stock_resp_auto
        ! POOLS(1,8) ! WATER IN ROOT ZONE ASSIGNED LATER
        POOLS(1,9) = stock_storage_organ

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
        fol_frac_intpol = 0d0
        stem_frac_intpol = 0d0
        root_frac_intpol = 0d0

        ! SHOULD TURN THIS INTO A SUBROUTINE CALL AS COMMON TO BOTH DEFAULT AND CROPS
        if (.not.allocated(deltat_1)) then

           allocate(deltat_1(nodays),wSWP_time(nodays),gs_demand_supply_ratio(nodays), &
                    gs_total_canopy(nodays), gb_total_canopy(nodays), &
                    canopy_par_MJday_time(nodays), cica_time(nodays))
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
      ! close are we to maxing out supply
      gs_demand_supply_ratio(n) = (stomatal_conductance - minimum_conductance) &
                                / (potential_conductance - minimum_conductance)
      ! Store the canopy level stomatal conductance (mmolH2O/m2/day)
      gs_total_canopy(n) = stomatal_conductance

      ! reallocate for crop model timings
      doy = met(6,n)

      ! GPP (gC.m-2.day-1)
      if (lai > vsmall .and. stomatal_conductance > vsmall) then
         call acm_gpp_stage_1 ; GPP_out(n) = max(0d0,acm_gpp_stage_2(stomatal_conductance))
         ! Estimate the internal to external CO2 concentration ratio
         cica_time(n) = ci / co2
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
      ! labile to foliage
      FLUXES(n,8) = (alloc_from_labile + resp_cost_labile_to_foliage) * steps_in_day !/deltat(n)
      ! alloc to storage organ
      FLUXES(n,9) = alloc_to_storage_organ * steps_in_day !/deltat(n)
      ! total leaf litter production
      FLUXES(n,10) = litterfall_foliage * steps_in_day !/deltat(n)
      ! total wood litter production
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
      ! root pool
      POOLS(n+1,3) = stock_roots
      ! wood pool
      POOLS(n+1,4) = stock_stem
      ! litter pool
      POOLS(n+1,5) = stock_litter
      ! som pool
      POOLS(n+1,6) = stock_soilOrgMatter
      ! autotrophic pool
      POOLS(n+1,7) = stock_resp_auto
      ! POOLS(n+1,8) = soil surface water content
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

  end subroutine CARBON_MODEL
  !
  !------------------------------------------------------------------
  !
  subroutine acm_gpp_stage_1

    ! Estimate the light and temperature limited photosynthesis components.
    ! See acm_gpp_stage_2() for estimation of CO2 supply limitation and
    ! combination of light, temperature and CO2 co-limitation

    implicit none

    !
    ! Metabolic limited photosynthesis
    !

    ! maximum rate of temperature and nitrogen (canopy efficiency) limited
    ! photosynthesis (gC.m-2.day-1 -> umolC/m2/day)
    metabolic_limited_photosynthesis = gC_to_umol*lai*ceff*opt_max_scaling(pn_max_temp,-1d6,pn_opt_temp,pn_kurtosis,leafT)

    !
    ! Light limited photosynthesis
    !

    ! calculate light limted rate of photosynthesis (gC.m-2.day-1)
    light_limited_photosynthesis = e0 * canopy_par_MJday

    !
    ! Stomatal conductance independent variables for diffusion limited
    ! photosynthesis
    !

    ! Canopy level boundary layer conductance unit change
    ! (m.s-1 -> mol.m-2.day-1) assuming sea surface pressure only.
    ! Note the ratio of H20:CO2 diffusion through leaf level boundary layer is
    ! 1.37 (Jones appendix 2).
    gb_mol = aerodynamic_conductance * seconds_per_day * convert_ms1_mol_1 * gb_H2O_CO2
    rb_mol_1 = (gb_mol)**(-1d0)

    ! Temperature adjustments for Michaelis-Menten coefficients
    ! for CO2 (kc) and O2 (ko) and CO2 compensation point
    ! See McMurtrie et al., (1992) Australian Journal of Botany, vol 40, 657-677
    co2_half_sat   = arrhenious(kc_saturation,kc_half_sat_conc,leafT)
    co2_comp_point = arrhenious(co2comp_saturation,co2comp_half_sat_conc,leafT)

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
    double precision, intent(in) :: gs

    ! declare local variables
    double precision :: pp, qq, mult, rc, pd

    !
    ! Diffusion limited photosynthesis
    !

    ! Daily canopy conductance (mmolH2O.m-2.s-1-> molCO2.m-2.day-1)
    ! The ratio of H20:CO2 diffusion is 1.646259 (Jones appendix 2).
    ! i.e. gcH2O*1.646259 = gcCO2 then all multiplied by 86400 seconds
    !
    ! Combining in series the stomatal and boundary layer conductances
    ! to make canopy resistence
    rc = (gs*gs_H2Ommol_CO2mol_day) ** (-1d0) + rb_mol_1

    ! pp and qq represent limitation by metabolic (temperature & N) and
    ! diffusion (co2 supply) respectively
    pp = metabolic_limited_photosynthesis*rc ; qq = co2_comp_point-co2_half_sat
    ! calculate internal CO2 concentration (ppm or umol/mol)
    mult = co2+qq-pp
    ci = 0.5d0*(mult+sqrt((mult*mult)-4d0*(co2*qq-pp*co2_comp_point)))

    ! calculate CO2 limited rate of photosynthesis (gC.m-2.day-1)
    ! Then scale to day light period as this is then consistent with the light
    ! capture period (1/24 = 0.04166667)
    pd = ((co2-ci)/rc) * umol_to_gC * dayl_hours_fraction

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
  !------------------------------------------------------------------
  !
  double precision function find_gs_iWUE(gs_in)

    ! Calculate CO2 limited photosynthesis as a function of metabolic limited
    ! photosynthesis (pn), atmospheric CO2 concentration and stomatal
    ! conductance (gs_in). Photosynthesis is calculated twice to allow for
    ! testing of senstivity to iWUE.

    ! arguments
    double precision, intent(in) :: gs_in

    !!!!!!!!!!
    ! Optimise intrinsic water use efficiency
    !!!!!!!!!!

    ! Determine impact of gs increment on pd and how far we are from iWUE
    find_gs_iWUE = iWUE_step - ((acm_gpp_stage_2(gs_in + delta_gs) - acm_gpp_stage_2(gs_in))*lai_1)

    ! Remember to return back to the user
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
    gpp_low = acm_gpp_stage_2(gs_in)
    call calculate_transpiration(evap_low)

    ! Increment gs
    gs_high = gs_in + delta_gs
    ! now assign the incremented estimate
    stomatal_conductance = gs_high
    ! estimate photosynthesis with incremented gs
    gpp_high = acm_gpp_stage_2(gs_high)
    call calculate_transpiration(evap_high)

    ! estimate marginal return on GPP for water loss, less water use efficiency
    ! criterion (gC.kgH2O-1.m-2.s-1)
    find_gs_WUE = (((gpp_high - gpp_low)/(evap_high - evap_low)) * lai_1) - iWUE

    ! return original stomatal value back into memory
    stomatal_conductance = gs_store

    ! remember to return back to the user
    return

  end function find_gs_WUE
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_stomatal_conductance

    ! Determines 1) an approximation of canopy conductance (gc) mmolH2O.m-2.s-1
    ! based on potential hydraulic flow, air temperature and absorbed radiation.
    ! 2) calculates absorbed shortwave radiation (W.m-2) as function of LAI

    implicit none

    ! local variables
    double precision :: denom, iWUE_lower, iWUE_upper
    double precision, parameter :: max_gs = 2000d0, &  ! mmolH2O.m-2.s-1 (leaf area)
                                   min_gs = 1d0, &     ! mmolH2O.m-2.s-1 (leaf area)
                                   tol_gs = 10d0       ! mmolH2O.m-2.s-1 (leaf area)

    !!!!!!!!!!
    ! Calculate stomatal conductance under H2O and CO2 limitations
    !!!!!!!!!!

    if (aerodynamic_conductance > vsmall .and. total_water_flux > vsmall) then

        ! Determine potential water flow rate (mmolH2O.m-2.dayl-1)
        max_supply = total_water_flux * seconds_per_day

        ! Pass minimum conductance from local parameter to global value
        ! There is uncertainty whether this should be a leaf area scaled value...
        minimum_conductance = min_gs * lai

        ! Invert Penman-Monteith equation to give gs (m.s-1) needed to meet
        ! maximum possible evaporation for the day.
        ! This will then be reduced based on CO2 limits for diffusion based
        ! photosynthesis
        denom = slope * ((canopy_swrad_MJday * 1d6 * dayl_seconds_1) + canopy_lwrad_Wm2) &
              + (ET_demand_coef * aerodynamic_conductance)
        denom = (denom / (lambda * max_supply * mmol_to_kg_water * dayl_seconds_1)) - slope
        potential_conductance = aerodynamic_conductance / (denom / psych)

        ! convert m.s-1 to mmolH2O.m-2.s-1
        potential_conductance = potential_conductance * convert_ms1_mmol_1
        ! if conditions are dew forming then set conductance to maximum as we
        ! are not going to be limited by water demand
        if (potential_conductance <= 0d0 .or. potential_conductance > max_gs*lai) potential_conductance = max_gs*lai

        ! If there is a positive demand for water then we will solve for
        ! photosynthesis limits on gs through iterative solution
        delta_gs = 1d-3*lai ! mmolH2O/m2leaf/day
        ! Estimate inverse of LAI to avoid division in optimisation
        lai_1 = lai**(-1d0)
        ! Calculate stage one acm, temperature and light limitation which
        ! are independent of stomatal conductance effects
        call acm_gpp_stage_1

!        if (do_iWUE) then
            ! Intrinsic WUE optimisation
            ! Check that the water restricted water range brackets the root solution for the bisection
            iWUE_upper = find_gs_iWUE(potential_conductance) !; iWUE_lower = find_gs_iWUE(min_gs)
            if ( iWUE_upper * find_gs_iWUE(min_gs) > 0d0 ) then
                ! Then both proposals indicate that photosynthesis
                ! would be increased by greater opening of the stomata
                ! and is therefore water is limiting!
                stomatal_conductance = potential_conductance
                ! Exception being if both are positive - therefore assume
                ! lowest
                if (iWUE_upper > 0d0) stomatal_conductance = minimum_conductance
            else

                ! In all other cases iterate
                stomatal_conductance = zbrent('calculate_gs:find_gs_iWUE', &
                                              find_gs_iWUE,minimum_conductance,potential_conductance,tol_gs*lai,iWUE_step*0.10d0)

            end if
!            ! Empirical fit to outputs generated by bisection procedure.
!            ! Assumes that water supply is not limiting, thus there is still the need to estimate supply limit and apply as bookend.
!            ! Note also that the order of covariates reflects their importance in the prediction,
!            ! i.e. R > 0.9 just for first independent variable
!            pn = metabolic_limited_photosynthesis
!            pn_day = metabolic_limited_photosynthesis * dayl_hours_fraction
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
!            stomatal_conductance = max(min_gs,min(stomatal_conductance,potential_conductance))
!       else
!            ! WUE optimisation
!            stomatal_conductance = zbrent('calculate_gs:find_gs_WUE',find_gs_WUE,min_gs,potential_conductance,tol_gs,iWUE*0.10d0)
!       endif

    else

        ! if no LAI then there can be no stomatal conductance
        potential_conductance = max_gs ; minimum_conductance = vsmall
        stomatal_conductance = vsmall
        ! set minimum (computer) precision level flow
        max_supply = vsmall

    endif ! if aerodynamic conductance > vsmall

  end subroutine calculate_stomatal_conductance
  !
  !------------------------------------------------------------------
  !
  subroutine meteorological_constants(input_temperature,input_temperature_K,input_vpd_kPa)

    ! Determine some multiple use constants used by a wide range of functions
    ! All variables here are linked to air temperature and thus invarient between
    ! iterations and can be stored in memory...

    implicit none

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
    air_density_kg = 353d0/input_temperature_K
    ! Conversion ratio for m.s-1 -> mol.m-2.s-1
    convert_ms1_mol_1 = const_sfc_pressure / (input_temperature_K*Rcon)
    ! latent heat of vapourisation,
    ! function of air temperature (J.kg-1)
    lambda = 2501000d0-2364d0*input_temperature

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
                                  ,gs,gb   ! stomatal and boundary layer conductance (m.s-1)

    !!!!!!!!!!
    ! Estimate energy radiation balance (W.m-2)
    !!!!!!!!!!

    ! Absorbed shortwave radiation MJ.m-2.day-1 -> J.m-2.s-1
    canopy_radiation = canopy_lwrad_Wm2 + (canopy_swrad_MJday * 1d6 * dayl_seconds_1)

    !!!!!!!!!!
    ! Calculate canopy conductance (to water vapour)
    !!!!!!!!!!

    ! Change units of potential stomatal conductance
    ! (mmolH2O.m-2.s-1 -> m.s-1).
    ! Note assumption of sea surface pressure only
    gs = stomatal_conductance / convert_ms1_mmol_1
    ! Combine in series stomatal conductance with boundary layer
    gb = aerodynamic_conductance

    !!!!!!!!!!
    ! Calculate canopy evaporative fluxes (kgH2O/m2/day)
    !!!!!!!!!!

    ! Calculate numerator of Penman Monteith (kgH2O.m-2.day-1)
    ! NOTE: that restriction within water supply restriction is determined
    ! during stomatal conductance level.
    transpiration = ( ( (slope*canopy_radiation) + (ET_demand_coef*gb) ) &
                      / (lambda*(slope+(psych*(1d0+gb/gs)))) )*dayl_seconds

  end subroutine calculate_transpiration
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_wetcanopy_evaporation(wetcanopy_evap,storage)

    ! Estimates evaporation of canopy intercepted rainfall based on the Penman-Monteith model of
    ! evapotranspiration used to estimate SPA's daily evapotranspiration flux
    ! (kgH20.m-2.day-1).

    implicit none

    ! arguments
    double precision, intent(inout) :: storage      ! canopy water storage kgH2O/m2
    double precision, intent(out) :: wetcanopy_evap ! kgH2O.m-2.day-1

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
    wetcanopy_evap = max(0d0,(((slope*canopy_radiation) + (ET_demand_coef*gb)) / (lambda*(slope+psych))) * seconds_per_day)

    ! assuming there is any rainfall, currently water on the canopy or dew formation
    if (rainfall > 0d0 .or. storage > 0d0) then
        ! Update based on canopy water storage
        call canopy_interception_and_storage(wetcanopy_evap,storage)
    else
        ! there is no water movement possible
        intercepted_rainfall = 0d0 ; wetcanopy_evap = 0d0
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
    double precision :: local_temp &
                   ,soil_radiation & ! isothermal net radiation (W/m2)
                            ,esurf & ! see code below
                             ,esat & ! soil air space saturation vapour pressure
                              ,gws   ! water vapour conductance through soil air space (m.s-1)

    ! oC -> K for local temperature value
    local_temp = maxt + freeze

    !!!!!!!!!!
    ! Estimate energy radiation balance (W.m-2)
    !!!!!!!!!!

    ! Absorbed shortwave radiation MJ.m-2.day-1 -> J.m-2.s-1
    soil_radiation = soil_lwrad_Wm2 + (soil_swrad_MJday * 1d6 * dayl_seconds_1)

    !!!!!!!!!!
    ! Calculate soil evaporative fluxes (kgH2O/m2/day)
    !!!!!!!!!!

    ! calculate saturated vapour pressure (kPa), function of temperature.
    esat = 0.1d0 * exp( 1.80956664d0 + ( 17.2693882d0 * local_temp - 4717.306081d0 ) / ( local_temp - 35.86d0 ) )
    air_vapour_pressure = esat - vpd_kPa

    ! Soil conductance to water vapour diffusion (m s-1)...
    gws = porosity(1) * water_vapour_diffusion / (tortuosity*drythick)

    ! vapour pressure in soil airspace (kPa), dependent on soil water potential
    ! - Jones p.110. partial_molar_vol_water. Less vapour pressure of the air to
    ! estimate the deficit between soil and canopy air spaces
    esurf = (esat * exp( 1d6 * SWP(1) * partial_molar_vol_water / (Rcon * local_temp) )) - air_vapour_pressure

    ! Estimate potential soil evaporation flux (kgH2O.m-2.day-1)
    soilevap = ( ((slope*soil_radiation) + (air_density_kg*cpair*esurf*soil_conductance)) &
               / (lambda*(slope+(psych*(soil_conductance/gws)))) ) * dayl_seconds

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
    double precision :: local_lai, &
           mixing_length_momentum, & ! mixing length parameter for momentum (m)
            length_scale_momentum    ! length scale parameter for momentum (m)

    ! Restrict LAI used here to greater than a minium value which prevents un-realistic outputs
    local_lai = max(min_lai,lai)

    ! calculate the zero plane displacement and roughness length
    call z0_displacement(ustar_Uh,local_lai)
    ! calculate friction velocity at tower height (reference height ) (m.s-1)
    ! WARNING neutral conditions only; see WRF module_sf_sfclay.F for 'with
    ! stability versions'
    !    ustar = (wind_spd / log((tower_height-displacement)/roughl)) * vonkarman
    ustar = wind_spd * ustar_Uh

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
    call calculate_soil_conductance(mixing_length_momentum,local_lai)

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
    Sh_forced = 1.018537d0*(sqrt((leaf_width*canopy_wind)/kinematic_viscosity))
    ! Estimate the the forced conductance of water vapour
    gv_forced = water_vapour_diffusion*Sh_forced*leaf_width_coef * lai

  end subroutine average_leaf_conductance
  !
  !------------------------------------------------------------------
  !
  subroutine log_law_decay

    ! Standard log-law above canopy wind speed (m.s-1) decay under neutral
    ! conditions.
    ! See Harman & Finnigan 2008; Jones 1992 etc for details.

    implicit none

    ! log law decay, NOTE: given canopy height (9 m) the log function reduces
    ! to a constant value down to ~ 7 decimal place (0.3161471806). Therefore
    ! 1/vonkarman * 0.31 = 0.7710906
    canopy_wind = ustar * vonkarman_1 * log((canopy_height-displacement) / roughl)

    ! set minimum value for wind speed at canopy top (m.s-1)
!    canopy_wind = max(min_wind,canopy_wind)

  end subroutine log_law_decay
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_field_capacity

    ! field capacity calculations for saxton eqns !

    implicit none

    ! local variables..
    integer :: i
    double precision :: x1, x2

    x1 = 0.1d0 ; x2 = 0.7d0 ! low/high guess
    do i = 1 , nos_soil_layers+1
       water_retention_pass = i
       ! field capacity is water content at which SWP = -10 kPa
       field_capacity(i) = zbrent('water_retention:water_retention_saxton_eqns', &
                                   water_retention_saxton_eqns , x1 , x2 , 0.001d0, 0d0 )
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
    double precision, parameter :: nos_layers = 4d0 & ! Number of canopy layers in source model
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
    transmitted_fraction = exp(decay * lai * clump)

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
    canopy_transmitted_fraction = exp(decay * lai * 0.5d0 * clump)
    trans_lw_fraction     = 0.02d0 * canopy_transmitted_fraction
    reflected_lw_fraction = 0.02d0 * canopy_transmitted_fraction
    ! Absorption is the residual
    absorbed_lw_fraction = 1d0 - trans_lw_fraction - reflected_lw_fraction

    ! Calculate the potential absorption of longwave radiation lost from the
    ! canopy to soil / sky.
    ! NOTE: That assuming the Beer's law emission from a single canopy layer leads to a rough 50 % underestimate of LW emission.
    !       This is why there is the nos_layer correction here
    canopy_release_fraction = (1d0 - (max_lai_lwrad_release*lai) / (lai+lai_half_lwrad_release)) &
                            * (1d0 - exp(decay * lai * nos_layers_1 * clump)) * nos_layers

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
    sky_lwrad_Wm2 = lwrad * reflected_lw_fraction

    !!!!!!!!!!
    ! Distribute longwave from soil
    !!!!!!!!!!

    ! Calculate longwave radiation coming up from the soil plus the radiation
    ! which is reflected
    canopy_absorption_from_soil = longwave_release_soil + (soil_incident_from_sky * (1d0-emissivity))
    ! First how much directly bypasses the canopy...
    sky_lwrad_Wm2 = sky_lwrad_Wm2 + (canopy_absorption_from_soil * transmitted_fraction)
    canopy_absorption_from_soil = canopy_absorption_from_soil * (1d0 - transmitted_fraction)
    ! Second, use this to estimate the longwave returning to the sky
    sky_lwrad_Wm2 = sky_lwrad_Wm2 + (canopy_absorption_from_soil * trans_lw_fraction)
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
    double precision :: delta_iso

    ! Estimate shortwave radiation balance
    call calculate_shortwave_balance
    ! Estimate isothermal long wave radiation balance
    call calculate_longwave_isothermal(meant,meant)
    ! Apply linear correction to soil surface isothermal->net longwave radiation
    ! balance based on absorbed shortwave radiation
    delta_iso = (soil_iso_to_net_coef_LAI * lai) + &
                (soil_iso_to_net_coef_SW * (soil_swrad_MJday * 1d6 * seconds_per_day_1)) + &
                soil_iso_to_net_const
    ! In addition to the iso to net adjustment, SPA analysis shows that soil net never gets much below zero
    soil_lwrad_Wm2 = max(-0.01d0,soil_lwrad_Wm2 + delta_iso)
    ! Apply linear correction to canopy isothermal->net longwave radiation
    ! balance based on absorbed shortwave radiation
    delta_iso = (canopy_iso_to_net_coef_LAI * lai) + &
                (canopy_iso_to_net_coef_SW * (canopy_swrad_MJday * 1d6 * seconds_per_day_1)) + &
                canopy_iso_to_net_const
    canopy_lwrad_Wm2 = canopy_lwrad_Wm2 + delta_iso

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
    double precision :: balance                     &
                       ,transmitted_fraction        &
                       ,canopy_transmitted_fraction &
                       ,absorbed_nir_fraction_soil  &
                       ,absorbed_par_fraction_soil  &
                       ,fsnow,par,nir               &
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
    double precision, parameter :: clump = 1d0    & ! Clumping factor (1 = uniform, 0 totally clumped, mean = 0.75)
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
    transmitted_fraction = exp(decay * lai * clump)

    ! Second, of the radiation which is incident on the canopy what fractions
    ! are transmitted through, reflected from or absorbed by the canopy

    canopy_transmitted_fraction = exp(decay * lai * 0.5d0 * clump)

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
    par = sw_par_fraction * swrad
    nir = (1d0 - sw_par_fraction) * swrad

    ! Estimate the radiation which directly bypasses the canopy...
    trans_par_MJday = par * transmitted_fraction
    trans_nir_MJday = nir * transmitted_fraction
    ! ...and update the canopy intercepted radiation
    par = par - trans_par_MJday
    nir = nir - trans_nir_MJday

    ! Estimate incoming shortwave radiation absorbed, transmitted and reflected
    ! by the canopy (MJ.m-2.day-1)
    canopy_par_MJday = par * absorbed_par_fraction
    canopy_nir_MJday = nir * absorbed_nir_fraction
    trans_par_MJday = trans_par_MJday + (par * trans_par_fraction)
    trans_nir_MJday = trans_nir_MJday + (nir * trans_nir_fraction)
    refl_par_MJday = par * reflected_par_fraction
    refl_nir_MJday = nir * reflected_nir_fraction

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

    ! Then the radiation incident and ultimately absorbed by the soil surface
    ! itself (MJ.m-2.day-1)
    soil_par_MJday = trans_par_MJday * absorbed_par_fraction_soil
    soil_nir_MJday = trans_nir_MJday * absorbed_nir_fraction_soil
    ! combine totals for use is soil evaporation
    soil_swrad_MJday = soil_nir_MJday + soil_par_MJday

    !!!!!!!!!
    ! Estimate canopy absorption of soil reflected shortwave radiation
    ! This additional reflection / absorption cycle is needed to ensure > 0.99
    ! of incoming radiation is explicitly accounted for in the energy balance.
    !!!!!!!!!

    ! calculate multiple use variables
    par = trans_par_MJday-soil_par_MJday
    nir = trans_nir_MJday-soil_nir_MJday
    ! how much of the reflected radiation directly bypasses the canopy...
    refl_par_MJday = refl_par_MJday + (par * transmitted_fraction)
    refl_nir_MJday = refl_nir_MJday + (nir * transmitted_fraction)
    ! ...and update the canopy on this basis
    par = par * (1d0-transmitted_fraction)
    nir = nir * (1d0-transmitted_fraction)

    ! Update the canopy radiation absorption based on the reflected radiation
    ! (MJ.m-2.day-1)
    canopy_par_MJday = canopy_par_MJday + (par * absorbed_par_fraction)
    canopy_nir_MJday = canopy_nir_MJday + (nir * absorbed_nir_fraction)
    ! Update the total radiation reflected back into the sky, i.e. that which is
    ! now transmitted through the canopy
    refl_par_MJday = refl_par_MJday + (par * trans_par_fraction)
    refl_nir_MJday = refl_nir_MJday + (nir * trans_nir_fraction)

    ! Combine to estimate total shortwave canopy absorbed radiation
    canopy_swrad_MJday = canopy_par_MJday + canopy_nir_MJday

    ! check energy balance
!    balance = swrad - canopy_par_MJday - canopy_nir_MJday - refl_par_MJday -
!    refl_nir_MJday - soil_swrad_MJday
!    if ((balance - swrad) / swrad > 0.01) then
!        print*,"SW residual frac = ",(balance - swrad) / swrad,"SW residual =
!        ",balance,"SW in = ",swrad
!    endif

  end subroutine calculate_shortwave_balance
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_Rtot

    ! Purpose of this subroutine is to calculate the minimum soil-root hydraulic
    ! resistance input into ACM. The approach used here is identical to that
    ! found in SPA.

    ! local variables
    integer :: i, rooted_layer
    double precision :: bonus, &
                        transpiration_resistance,root_reach_local, &
                        root_depth_50, slpa, mult, prev, exp_func
    double precision, dimension(nos_root_layers) :: Rcond_layer, &
                                                    root_mass,  &
                                                    root_length
    double precision, parameter :: rootdist_tol = 13.81551d0 ! log(1d0/rootdist_tol - 1d0) were rootdist_tol = 1d-6
                                  !rootdist_tol = 1d-6!, & ! Root density assessed for the max rooting depth
                                   !root_depth_frac_50 = 0.25d0 ! fractional soil depth above which 50 %
                                                               ! of the root mass is assumed to be located

    ! reset water flux
    total_water_flux = 0d0 ; water_flux_mmolH2Om2s = 0d0 ; wSWP = 0d0 ; rSWP = 0d0
    slpa = 0d0 ; root_length = 0d0 ; root_mass = 0d0 ; Rcond_layer = 0d0
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
    demand = max(0d0, (SWP(1:nos_root_layers) - (head*canopy_height)) - minlwp )
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
    slpa = rootdist_tol / root_reach
    prev = 1d0
    do i = 1, nos_root_layers
       ! Determine the exponential function for the current cumulative depth
       exp_func = exp(-slpa * sum(layer_thickness(1:i)))
       ! Calculate the difference in the integral between depths, i.e. the proportion of root in the current volume
       mult = prev - (1d0 - (1d0/(1d0+exp_func)) + (0.5d0 * exp_func))
       ! Assign fine roo the the current layer...
       root_mass(i) = fine_root_biomass * mult
       ! and determine the associated amount of root
       root_length(i) = root_mass(i) * root_mass_length_coef_1
       prev = prev - mult
       ! If there is root in the current layer then we should calculate the resistances
       if (root_mass(i) > 0d0) then
           ! Track the deepest root layer assessed
           rooted_layer = i
           ! if there is root then there is a water flux potential...
           root_reach_local = min(root_reach,layer_thickness(i))
           ! calculate and accumulate steady state water flux in mmol.m-2.s-1
           call plant_soil_flow(i,root_length(i),root_mass(i) &
                               ,demand(i),root_reach_local &
                               ,transpiration_resistance,Rcond_layer(i))
       else
           ! ...if there is not then we wont have any below...
           exit
       end if ! root present in current layer?
    end do ! nos_root_layers
    ! Turn the output resistance into conductance
    Rcond_layer = Rcond_layer**(-1d0)

    ! if freezing then assume soil surface is frozen, therefore no water flux
    if (meant < 1d0) then
        water_flux_mmolH2Om2s(1) = 0d0
        Rcond_layer(1) = 0d0
    end if

    ! calculate sum value (mmolH2O.m-2.s-1)
    total_water_flux = sum(water_flux_mmolH2Om2s)
    ! wSWP based on the conductance due to the roots themselves.
    ! The idea being that the plant may hedge against growth based on the majority of the
    ! profile being dry while not losing leaves within some toleration.
    rSWP = sum(SWP(1:rooted_layer) * (Rcond_layer(1:rooted_layer) / sum(Rcond_layer(1:rooted_layer))))
    if (total_water_flux <= vsmall) then
        ! Set values for no water flow situation
        uptake_fraction = (layer_thickness(1:nos_root_layers) / sum(layer_thickness(1:nos_root_layers)))
        ! Estimate weighted soil water potential based on fractional extraction from soil layers
        wSWP = sum(SWP(1:nos_root_layers) * uptake_fraction(1:nos_root_layers))
        total_water_flux = 0d0
      else
        ! calculate weighted SWP and uptake fraction
        uptake_fraction(1:nos_root_layers) = water_flux_mmolH2Om2s(1:nos_root_layers) / total_water_flux
        ! Estimate weighted soil water potential based on fractional extraction from soil layers
        wSWP = sum(SWP(1:nos_root_layers) * uptake_fraction(1:nos_root_layers))
    endif

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
                                                  ! enters as potential but leaves as water balance adjusted.
                                                  ! Note that this assumes a completely wet leaf surface
    ! local variables
    integer :: i
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
    through_fall = exp(CanIntFrac*lai*clump)
    ! maximum canopy storage (mm); minimum is applied to prevent errors in
    ! drainage calculation. Assume minimum capacity due to wood.
    max_storage = max(min_storage,CanStorFrac*lai)
    ! caclulate inverse for efficient calculations below
    max_storage_1 = max_storage**(-1d0)
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

    ! add rain to the canopy and overflow as needed
    storage = storage + daily_addition

    if (storage > max_storage) then

        if (potential_evaporation > 0d0) then

            ! assume co-access to available water above max_storage by both drainage and
            ! evaporation. Water below max_storage is accessable by evaporation only.

            ! Trapezium rule for approximating integral of drainage rate.
            ! Allows estimation of the mean drainage rate between starting
            ! point and max_storage, thus the time period appropriate for co-access can be
            ! quantified. NOTE 1440 = minutes / day
            ! General Formula: integral(rate) = 0.5 * h((y0 + yn) + 2(y1 + y2 + ... yn-1)
            ! Where h id the size of the section, y0 is the maximum rate, yn is the final rate.
            dx = (storage - max_storage)*0.5d0
            tmp(1) = storage ; tmp(2) = max_storage ; tmp(3) = storage-dx
            tmp = exp(a + (RefDrainCoef*tmp))
            potential_drainage_rate = 0.5d0 * dx * ((tmp(1) + tmp(2)) + 2d0 * tmp(3)) * 1440d0
            ! To protect against un-realistic drainage rates
            ! due to very high rainfall rates
            potential_drainage_rate = min(potential_drainage_rate,vlarge)

            dz = storage-max_storage
            ! limit based on available water if total demand is greater than excess
            co_mass_balance = (dz / (potential_evaporation + potential_drainage_rate))
            evap_rate = potential_evaporation * co_mass_balance
            drain_rate = potential_drainage_rate * co_mass_balance

            ! Estimate evaporation from remaining water (i.e. that left after
            ! initial co-access of evaporation and drainage).
            ! Assume evaporation is now restricted by:
            ! 1) energy already spent on evaporation (the -evap_rate) and
            ! 2) linear increase in surface resistance as the leaf surface
            ! dries (i.e. the 0.5).
            evap_rate = evap_rate + min((potential_evaporation - evap_rate) * 0.5d0, storage - evap_rate - drain_rate)

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
    intercepted_rainfall = intercepted_rainfall - (through_fall * seconds_per_day_1)

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
    integer :: day, a, i
    double precision :: depth_change, water_change, initial_soilwater, balance, mass_check
    double precision, dimension(nos_root_layers) :: avail_flux, evaporation_losses, pot_evap_losses

    ! set soil water exchanges
    underflow = 0d0 ; runoff = 0d0 ; corrected_ET = 0d0 ; evaporation_losses = 0d0 ; pot_evap_losses = 0d0
    initial_soilwater = 1d3 * sum(soil_waterfrac(1:nos_soil_layers) * layer_thickness(1:nos_soil_layers))

    ! Assume leaf transpiration is drawn from the soil based on the
    ! update_fraction estimated in calculate_Rtot
    pot_evap_losses = ET_leaf * uptake_fraction
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
                         - (pot_evap_losses * days_per_step * 1d-3))

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
       ! if rainfall is probably liquid / soil surface is probably not frozen
       if (rainfall_in > 0d0) then
           ! reset soil water change variable
           waterchange = 0d0
           call infiltrate(rainfall_in * days_per_step)
           ! update soil profiles. Convert fraction into depth specific values
           ! (rather than m3/m3) then update fluxes
           soil_waterfrac(1:nos_soil_layers) = soil_waterfrac(1:nos_soil_layers) &
                                             + (waterchange(1:nos_soil_layers) / layer_thickness(1:nos_soil_layers))
           ! soil waterchange variable reset in gravitational_drainage()
       endif ! is there any rain to infiltrate?

       ! determine drainage flux between surface -> sub surface
       call gravitational_drainage(nint(days_per_step))

       ! Pass information to the output ET variable
       corrected_ET = sum(pot_evap_losses)
       ! apply time step correction kgH2O/m2/step -> kgH2O/m2/day
       underflow = underflow * days_per_step_1
       runoff = runoff * days_per_step_1

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
          do a = 1, nos_root_layers ! note: timed comparison between "where" and do loop supports do loop for smaller vectors
             if (evaporation_losses(a) > avail_flux(a)) evaporation_losses(a) = avail_flux(a) * 0.999d0
          end do
          ! this will update the ET estimate outside of the function
          ! days_per_step corrections happens outside of the loop below
          corrected_ET = corrected_ET + sum(evaporation_losses)

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
          ! if rainfall is probably liquid / soil surface is probably not frozen
          if (rainfall_in > 0d0) then
              ! reset soil water change variable
              waterchange = 0d0
              call infiltrate(rainfall_in)
              ! update soil profiles. Convert fraction into depth specific values
              ! (rather than m3/m3) then update fluxes
              soil_waterfrac(1:nos_soil_layers) = soil_waterfrac(1:nos_soil_layers) &
                                                + (waterchange(1:nos_soil_layers) / layer_thickness(1:nos_soil_layers))
              ! soil waterchange variable reset in gravitational_drainage()
          endif ! is there any rain to infiltrate?

          !!!!!!!!!!
          ! Gravitational drainage
          !!!!!!!!!!

          ! Determine drainage flux between surface -> sub surface
          call gravitational_drainage(1)

       end do ! days_per_step

       ! apply time step correction kgH2O/m2/step -> kgH2O/m2/day
       corrected_ET = corrected_ET * days_per_step_1
       underflow = underflow * days_per_step_1
       runoff = runoff * days_per_step_1

    end if ! water_change > 0

    !!!!!!!!!!
    ! Update soil layer thickness
    !!!!!!!!!!

    depth_change = (top_soil_depth+min_layer) ; water_change = 0
    ! if roots extent down into the bucket
    if (root_reach > depth_change .and. previous_depth <= depth_change) then

        !!!!!!!!!!
        ! Soil profile is within the bucket layer (layer 3)
        !!!!!!!!!!

        if (previous_depth > depth_change) then
            ! how much has root depth extended since last step?
            depth_change = root_reach - previous_depth
        else
            ! how much has root depth extended since last step?
            depth_change = root_reach - depth_change
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
            layer_thickness(1) = top_soil_depth
            layer_thickness(2) = root_reach - top_soil_depth
            layer_thickness(3) = max_depth - sum(layer_thickness(1:2))

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
            layer_thickness(1) = top_soil_depth
            layer_thickness(2) = root_reach - top_soil_depth
            layer_thickness(3) = max_depth - sum(layer_thickness(1:2))

            ! keep track of the previous rooting depth
            previous_depth = root_reach

        else

            ! keep track of the previous rooting depth
            previous_depth = previous_depth

        end if ! depth change

    else if (root_reach < depth_change .and. previous_depth > depth_change) then

        !!!!!!!!!!
        ! Model has explicitly contracted from the bucket layer
        !!!!!!!!!!

        ! In this circumstance we want to return the soil profile to it's
        ! default structure with a minimum sized third layer
        depth_change = previous_depth - depth_change

        ! determine how much water is lost from the old volume of soil
        water_change = soil_waterfrac(nos_root_layers) * depth_change
        ! now assign that new volume of water to the deep rooting layer
        soil_waterfrac(nos_soil_layers) = ((soil_waterfrac(nos_soil_layers)*layer_thickness(nos_soil_layers))+water_change) &
                                        / (layer_thickness(nos_soil_layers)+depth_change)

        ! explicitly update the soil profile if there has been rooting depth
        ! changes
        layer_thickness(1) = top_soil_depth
        layer_thickness(2) = min_layer
        layer_thickness(3) = max_depth - sum(layer_thickness(1:2))

        ! keep track of the previous rooting depth
        previous_depth = min_layer

    else ! root_reach > (top_soil_depth + min_layer)

        ! if we are outside of the range when we need to consider rooting depth changes keep track in case we move into a zone when we do
        previous_depth = previous_depth

    endif ! root reach beyond top layer

    ! finally update soil water potential
    call soil_water_potential

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
  subroutine infiltrate(rainfall_in)

    ! Takes surface_watermm and distributes it among top !
    ! layers. Assumes total infilatration in timestep.   !
    ! NOTE: Assumes that any previous water movement due to infiltration and evaporation
    !       has already been updated in soil mass balance

    implicit none

    ! arguments
    double precision, intent(in) :: rainfall_in ! rainfall (kg.m-2.day-1)

    ! local argumemts
    integer :: i
    double precision :: add, & ! surface water available for infiltration (m)
                      wdiff    ! available space in a given soil layer for water to fill (m)

    ! convert rainfall water from mm -> m (or kgH2O.m-2.day-1 -> MgH2O.m-2.day-1)
    add = rainfall_in * 1d-3

    do i = 1 , nos_soil_layers

       ! is the input of water greater than available space
       ! if so fill and subtract from input and move on to the next
       ! layer determine the available pore space in current soil layer
       wdiff = max(0d0,(porosity(i)-soil_waterfrac(i)) * layer_thickness(i))

       if (add > wdiff) then
           ! if so fill and subtract from input and move on to the next layer
           waterchange(i) = waterchange(i) + wdiff
           add = add - wdiff
       else
           ! otherwise infiltate all in the current layer
           waterchange(i) = waterchange(i) + add
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
  subroutine gravitational_drainage(time_period_days)

    ! Integrator for soil gravitational drainage.
    ! Due to the longer time steps undertake by ACM / DALEC and the fact that
    ! drainage is a concurrent processes we assume that drainage occurs at
    ! the bottom of the column first creating space into which water can drain
    ! from the top down. Therefore we draing from the bottom first and then the top.
    ! NOTE: Assumes that any previous water movement due to infiltration and evaporation
    !       has already been updated in soil mass balance

    implicit none

    ! arguments
    integer, intent(in) :: time_period_days

    ! local variables..
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
    if (meant < 1d0) iceprop(1) = 1d0

    ! zero water fluxes
    waterchange = 0d0

    ! underflow is tracked in kgH2O/m2/day but estimated here in MgH2O/m2/day
    ! therefore we must convert
    underflow = underflow * 1d-3

    ! estimate potential drainage rate for the current time period
    liquid = soil_waterfrac(1:nos_soil_layers) * ( 1d0 - iceprop(1:nos_soil_layers) )
    ! estimate how much liquid is available to flow
    avail_to_flow = liquid - field_capacity(1:nos_soil_layers)
    ! trapezium rule scaler and the half-way point between current and field capacity
    dx = avail_to_flow*0.5d0 ; halfway = liquid - dx
    do t = 1, nos_soil_layers
       if (avail_to_flow(t) > 0d0) then
           ! Trapezium rule for approximating integral of drainage rate
           call calculate_soil_conductivity(t,liquid(t),tmp1)
           call calculate_soil_conductivity(t,field_capacity(t),tmp2)
           call calculate_soil_conductivity(t,halfway(t),tmp3)
           pot_drainage(t) = 0.5d0 * dx(t) * ((tmp1 + tmp2) + 2d0 * tmp3)
       else
           ! We are at field capacity currently even after rainfall has been infiltrated.
           ! Assume that the potential drainage rate is that at field capacity
           call calculate_soil_conductivity(t,field_capacity(t),pot_drainage(t))
       endif ! water above field capacity to flow?
    end do ! soil layers
    ! Scale potential drainage from per second to per day
    pot_drainage = pot_drainage * seconds_per_day

    ! Integrate drainage over each day until time period has been reached or
    ! each soil layer has reached field capacity
    t = 1
    do while (t < (time_period_days+1) .and. maxval(soil_waterfrac - field_capacity) > vsmall)

       ! Estimate liquid content and how much is available to flow / drain
       avail_to_flow = ( soil_waterfrac(1:nos_soil_layers) * (1d0 - iceprop(1:nos_soil_layers)) ) &
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
              change = min(unsat,min(pot_drainage(soil_layer),avail_to_flow(soil_layer))) * layer_thickness(soil_layer)
              ! update soil layer below with drained liquid
              waterchange( soil_layer + 1 ) = waterchange( soil_layer + 1 ) + change
              waterchange( soil_layer     ) = waterchange( soil_layer     ) - change

          end if ! some liquid water and drainage possible

       end do ! soil layers

       ! update soil water profile
       soil_waterfrac(1:nos_soil_layers) = soil_waterfrac(1:nos_soil_layers) &
                                         + (waterchange(1:nos_soil_layers)/layer_thickness(1:nos_soil_layers))
       ! estimate drainage from bottom of soil column (MgH2O/m2/day)
       ! NOTES: that underflow is reset outside of the daily soil loop
       underflow = underflow + waterchange(nos_soil_layers+1)

       ! Reset now we have moves that liquid
       waterchange = 0d0
       ! integerate through time period
       t = t + 1

    end do ! while condition

    ! convert underflow from MgH2O/m2/day -> kgH2O/m2/day
    underflow = underflow * 1d3

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

    ! Include some hardcoded boundaries for the Saxton equations
    ! NOTE: do loop was found to be faster than 'where' for small vectors
    do i = 1, nos_soil_layers
       if (soil_frac_sand(i) < 5d0) soil_frac_sand(i) = 5d0
       if (soil_frac_clay(i) < 5d0) soil_frac_clay(i) = 5d0
       if (soil_frac_sand(i) > 60d0) soil_frac_sand(i) = 60d0
    end do
    ! calculate soil porosity (m3/m3)
    call soil_porosity(soil_frac_clay,soil_frac_sand)
    ! calculate field capacity (m3/m-3)
    call calculate_field_capacity

    ! final sanity check for porosity
    do i = 1, nos_soil_layers+1
       if (porosity(i) < (field_capacity(i)+0.05d0)) porosity(i) = field_capacity(i) + 0.05d0
    end do

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
    double precision, intent(in) :: input_soilwater_frac ! initial soil water status as fraction of field capacity

    ! local variables
    integer :: i

    ! Default assumption to be field capacity
    soil_waterfrac = field_capacity
    SWP = SWP_initial

    ! if prior value has been given
    if (input_soilwater_frac > -9998d0) then
        ! calculate initial soil water fraction
        soil_waterfrac(1:nos_soil_layers) = input_soilwater_frac * field_capacity(1:nos_soil_layers)
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
  subroutine calculate_soil_conductance(lm,local_lai)

    ! proceedsure to solve for soil surface resistance based on Monin-Obukov
    ! similarity theory stability correction momentum & heat are integrated
    ! through the under canopy space and canopy air space to the surface layer
    ! references are Nui & Yang 2004; Qin et al 2002
    ! NOTE: conversion to conductance at end

    implicit none

    ! declare arguments
    double precision, intent(in) :: lm, local_lai

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

    integer :: i

    ! reformulation aims to remove if statement within loop to hopefully improve
    ! optimisation
    SWP(1:nos_soil_layers) = -0.001d0 * potA(1:nos_soil_layers) &
                           * soil_waterfrac(1:nos_soil_layers)**potB(1:nos_soil_layers)
    ! NOTE: profiling indiates that 'where' is slower for very short vectors
    do i = 1, nos_soil_layers
       if (SWP(i) < -20d0 .or. SWP(i) /= SWP(i)) SWP(i) = -20d0
    end do

  end subroutine soil_water_potential
  !
  !------------------------------------------------------------------
  !
  subroutine z0_displacement(ustar_Uh,local_lai)

    ! dynamic calculation of roughness length and zero place displacement (m)
    ! based on canopy height and lai. Raupach (1994)

    implicit none

    ! arguments
    double precision, intent(out) :: ustar_Uh ! ratio of friction velocity over wind speed at canopy top
    double precision, intent(in) :: local_lai
    ! local variables
    double precision  sqrt_cd1_lai
    double precision, parameter :: cd1 = 7.5d0,   & ! Canopy drag parameter; fitted to data
                                    Cs = 0.003d0, & ! Substrate drag coefficient
                                    Cr = 0.3d0,   & ! Roughness element drag coefficient
                          ustar_Uh_max = 0.3d0,   & ! Maximum observed ratio of
                                                    ! (friction velocity / canopy top wind speed) (m.s-1)
                          ustar_Uh_min = 0.2d0,   &
                                    Cw = 2d0, &  ! Characterises roughness sublayer depth (m)
                                 phi_h = 0.19314718056d0 ! Roughness sublayer influence function;

    ! describes the departure of the velocity profile from just above the
    ! roughness from the intertial sublayer log law

    ! Estimate canopy drag factor
    sqrt_cd1_lai = sqrt(cd1 * local_lai)

    ! calculate estimate of ratio of friction velocity / canopy wind speed.
    ! NOTE: under current min LAI and fixed canopy height (9 m) this ratio is
    ! fixed at 0.3
    ustar_Uh = 0.3d0
!    ustar_Uh = max(ustar_Uh_min,min(sqrt(Cs+Cr*local_lai*0.5d0),ustar_Uh_max))
!    ustar_Uh = sqrt(Cs+Cr*local_lai*0.5d0)

    ! calculate displacement (m); assume minimum lai 1.0 or 1.5 as height is not
    ! varied
    displacement = (1d0-((1d0-exp(-sqrt_cd1_lai))/sqrt_cd1_lai))*canopy_height

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
  subroutine calc_pools_crops(DS_LRRT,LRRT)

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

    ! Determine total allocation of labile C to NPP or Ra (gC.m-2.t-1)
    alloc_from_labile = turnover_rate_labile * resp_rate * dble(turnover_labile_switch)
    alloc_from_labile = stock_labile * (1d0-(1d0-alloc_from_labile) ** ts_length)
    ! Determine the respiratory cost of C transfer from labile pool to NPP (gC.m-2.t-1)
    resp_cost_labile_to_foliage = alloc_from_labile * resp_cost_labile_trans
    ! Determine the remaining allocation of labile C to NPP (gC.m-2.t-1)
    alloc_from_labile = alloc_from_labile - resp_cost_labile_to_foliage

!    ! respiratory cost of C transfer from labile pool to short-term pool (NPP) (gC.m-2.t-1)
!    resp_cost_labile_to_foliage = turnover_rate_labile * resp_cost_labile_trans * resp_rate &
!                                * ts_length * dble(turnover_labile_switch)
!    resp_cost_labile_to_foliage = stock_labile * min(1d0,resp_cost_labile_to_foliage)
!
!    ! allocation flux from labile C pool to NPP (gC.m-2.t-1)
!   alloc_from_labile = turnover_rate_labile * ( 1d0 - resp_cost_labile_trans ) * resp_rate &
!                      * ts_length * dble(turnover_labile_switch)
!   alloc_from_labile = stock_labile * min(1d0,alloc_from_labile)

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
    shoot_frac_intpol = 1d0 - root_frac_intpol         !
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
    litterfall_foliage = stock_foliage * (1d0-(1d0-RDR) ** ts_length)
    litterfall_stem    = stock_stem    * (1d0-(1d0-(DR*turnover_rate_stem*dble(stmob))) ** ts_length) ! remobstem
    litterfall_roots   = stock_roots   * (1d0-(1d0-turnover_rate_roots) ** ts_length)
!    litterfall_foliage = stock_foliage * min(1d0,ts_length * RDR)
!    litterfall_stem    = stock_stem    * min(1d0,ts_length * DR * turnover_rate_stem * dble(stmob)) ! remobstem
!    litterfall_roots   = stock_roots   * min(1d0,ts_length * turnover_rate_roots)

    ! remobilized C to NPP (from both leaves and stems) (gC.m-2.t-1)
    remob   = ( litterfall_foliage * 0.5d0 + litterfall_stem ) * ( 1d0 - resp_cost_labile_trans )
    ! respiratory cost of C transfer (conversion from starch to photosynthates) (gC.m-2.t-1)
    Raremob = ( litterfall_foliage * 0.5d0 + litterfall_stem ) * resp_cost_labile_trans

    ! for mass balance calculate the decompostion efficency
    decomp_efficency = decomposition_rate &
                     / (decomposition_rate+mineralisation_rate_litter)

    ! total litter decomposition
    decomposition = stock_litter * (1d0-(1d0-((decomposition_rate+mineralisation_rate_litter)*resp_rate)) ** ts_length)
    !decomposition = stock_litter * (decomposition_rate+mineralisation_rate_litter) * resp_rate * ts_length

    ! heterotrophic respiration component 1: mineralisation of litter C pool (gC.m-2.t-1)
    resp_h_litter = decomposition * (1d0 - decomp_efficency)
    ! heterotrophic respiration component 2:  mineralisation of organic matter C pool (gC.m-2.t-1)
    resp_h_soilOrgMatter = stock_soilOrgMatter * (1d0-(1d0-(mineralisation_rate_soilOrgMatter * resp_rate)) ** ts_length)  
    !resp_h_soilOrgMatter = stock_soilOrgMatter * min(1d0,mineralisation_rate_soilOrgMatter * resp_rate * ts_length)

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
    stock_resp_auto = stock_resp_auto + (frac_GPP_resp_auto * gpp_acm)
    ! autotrophic respiration; Ra (typically ~7% of respiratory pool) (gC.m-2.t-1)
    resp_auto = stock_resp_auto * (1d0-(1d0-turnover_rate_resp_auto) ** ts_length)
    !resp_auto = stock_resp_auto * min(1d0,turnover_rate_resp_auto * ts_length)
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

       ! loop over three crop "organs": 1) foliage 2) stems 3) root
       ! not necessary for storage organs, as all remaining C is allocated to
       ! these

       ! use different input for foliage and stem fractions, as they are
       ! relative to the total shoot (or aboveground) allocation,
       ! root is relative to total plant (above- and belowground) allocation.

       ! leaf development stages and corresponding fractions..

       ! interpolate between PdV allocation values with reference to
       ! developmental stage (DS)..
       fol_frac_intpol = interpolate( DS , DS_shoot , fol_frac , size(DS_shoot) )

       ! stem DS and fracs..
       stem_frac_intpol = interpolate( DS , DS_shoot , stem_frac , size(DS_shoot) )

       ! root DS and fracs..
       root_frac_intpol = interpolate( DS , DS_root , root_frac , size(DS_root) )

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
  !------------------------------------------------------------------
  !
  subroutine plant_soil_flow(root_layer,root_length,root_mass &
                            ,demand,root_reach_in,transpiration_resistance &
                            ,Rtot_layer)

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
    double precision, intent(out) :: Rtot_layer

    ! local arguments
    double precision :: soilR1, soilR2

    ! Estimate soil hydraulic resistance to water flow (MPa m2 s mmol-1)
    ! Note: 1) soil conductivity converted from m.s-1 -> m2.s-1.MPa-1 by head.
    !       2) soil resistance calculation in single line to reduce assignment costs
    soilR1 = ( log(root_radius_1*(root_length*pi)**(-0.5d0)) &
               /(two_pi*root_length*root_reach_in*(soil_conductivity(root_layer)*head_1))) &
             * 1d-9 * mol_to_g_water
    ! Calculates root hydraulic resistance (MPa m2 s mmol-1) in a soil-root zone
    soilR2 = root_resist / (root_mass*root_reach_in)
    ! Estimate the total hydraulic resistance for the layer
    Rtot_layer = transpiration_resistance + soilR1 + soilR2

    ! Estimate the soil to plant flow of water mmolH2O/m2/s
    water_flux_mmolH2Om2s(root_layer) = demand/Rtot_layer

    ! return
    return

  end subroutine plant_soil_flow
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

    !numerator   = t - 25d0
    !denominator = t + freeze
    arrhenious = a * exp( b * (t - 25d0) / (t + freeze) )

  end function arrhenious
  !
  !----------------------------------------------------------------------
  !
  double precision function opt_max_scaling( max_val, min_val , optimum , kurtosis , current )

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
    if (current >= max_val .or. current <= min_val) then
        opt_max_scaling = 0d0
    else
        opt_max_scaling = exp( kurtosis * log((max_val-current)/(max_val-optimum)) * (max_val-optimum) ) &
                        * exp( kurtosis * log((current-min_val)/(optimum-min_val)) * (optimum-min_val) )
    endif

  end function opt_max_scaling
  !
  !------------------------------------------------------------------
  !
  double precision function water_retention_saxton_eqns( xin )

    ! field capacity calculations for saxton eqns !

    implicit none

    ! arguments..
    double precision, intent(in) :: xin

    ! local variables..
    double precision :: soil_wp

    ! calculate the soil water potential (kPa)..
    soil_wp = -potA(water_retention_pass) * xin**potB(water_retention_pass)
    water_retention_saxton_eqns = soil_wp + 10d0    ! 10 kPa represents air-entry swp

    return

  end function water_retention_saxton_eqns
  !
  !--------------------------------------------------------------------------
  !
  double precision function linear_model_gradient(x,y,interval)

    ! Function to calculate the gradient of a linear model for a given depentent
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
    !sumsq_x = sum(x*x)
    ! calculate the sum of the product of xy
    !sum_product_xy = sum(x*y)
    ! calculate the gradient
    !linear_model_gradient = ( (dble(interval)*sum_product_xy) - (sum_x*sum_y) ) &
    !                      / ( (dble(interval)*sumsq_x) - (sum_x*sum_x) )
    ! Linear regression done as single line to reduce assignment requirements
    linear_model_gradient = ( (dble(interval)*sum(x*y)) - (sum_x*sum_y) ) &
                          / ( (dble(interval)*sum(x*x)) - (sum_x*sum_x) )

    ! for future reference here is how to calculate the intercept
!    intercept = ( (sum_y*sumsq_x) - (sum_x*sum_product_xy) ) &
!              / ( (dble(interval)*sumsq_x) - (sum_x*sum_x) )

    ! don't forget to return to the user
    return

  end function linear_model_gradient
  !
  !------------------------------------------------------------------
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
  !------------------------------------------------------------------
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
  !------------------------------------------------------------------
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
  double precision function zbrent( called_from , func , x1 , x2 , tol , toltol)

    ! This is a bisection routine. When ZBRENT is called, we provide a    !
    ! reference to a particular function and also two values which bound  !
    ! the arguments for the function of interest. ZBRENT finds a root of  !
    ! the function (i.e. the point where the function equals zero), that  !
    ! lies between the two bounds.                                        !
    ! There are five exit conditions:                                     !
    ! 1) The first proposal for the root of the function equals zero      !
    ! 2) The proposal range has been reduced to less then tol             !
    ! 3) The magnitude of the function is less than toltol                !
    ! 4) Maximum number of iterations has been reached                    !
    ! 5) The root of the function does now lie between supplied bounds    !
    ! For a full description see Press et al. (1986).                     !

    implicit none

    ! arguments..
    character(len=*),intent(in) :: called_from    ! name of procedure calling (used to pass through for errors)
    double precision,intent(in) :: tol, toltol, x1, x2

    ! Interfaces are the correct way to pass procedures as arguments.
    interface
      double precision function func( xval )
        double precision ,intent(in) :: xval
      end function func
    end interface

    ! local variables..
    integer            :: iter
    integer, parameter :: ITMAX = 8
    double precision   :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,tol0,xm
    double precision, parameter :: EPS = 6d-8

    ! calculations...
    a  = x1
    b  = x2
    fa = func( a )
    fb = func( b )
    tol0 = tol * 0.5d0

    ! Check that we haven't (by fluke) already started with the root...
    if ( abs(fa) < toltol ) then
        zbrent = a
        return
    elseif ( abs(fb) < toltol ) then
        zbrent = b
        return
    end if
    ! Ensure the supplied x-values give y-values that lie either
    ! side of the root and if not flag an error message...
!    if (fa * fb > 0d0) then
        ! tell me otherwise what is going on
!!       print*,"Supplied values must bracket the root of the function.",new_line('x'),  &
!!         "     ","You supplied x1:",x1,new_line('x'),                     &
!!         "     "," and x2:",x2,new_line('x'),                             &
!!         "     "," which give function values of fa :",fa,new_line('x'),  &
!!         "     "," and fb:",fb," .",new_line('x'),                        &
!!         " zbrent was called by: ",trim(called_from)
!    end if
    c = b
    fc = fb

    do iter = 1 , ITMAX

      ! If the new value (f(c)) doesn't bracket
      ! the root with f(b) then adjust it..
!      if ( sign(1d0,fb) .eq. sign(1d0,fc) ) then
      if (fb * fc > 0d0) then
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
      tol1 = EPS * abs(b) + tol0
      xm   = 0.5d0 * ( c - b )
!      if ( ( abs(xm) .le. tol1 ) .or. ( fb .eq. 0d0 ) ) then
      if ( ( abs(xm) .le. tol1 ) .or. ( abs(fb) < toltol ) ) then
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

    zbrent = b

  end function zbrent
  !
  !------------------------------------------------------------------
  !
!
!--------------------------------------------------------------------
!
end module CARBON_MODEL_MOD
