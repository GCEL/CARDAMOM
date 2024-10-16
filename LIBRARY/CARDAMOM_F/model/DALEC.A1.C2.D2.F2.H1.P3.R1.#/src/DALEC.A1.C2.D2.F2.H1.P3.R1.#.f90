
module CARBON_MODEL_MOD

  implicit none

  !!!!!!!!!!!
  ! Authorship contributions
  !
  ! This code contains a variant of the Data Assimilation Linked ECosystem (DALEC) model.
  ! This version of DALEC is derived from the following primary references:
  ! Bloom & Williams (2015), https://doi.org/10.5194/bg-12-1299-2015.
  ! Smallman et al., (2017), https://doi.org/10.1002/2016JG003520.
  ! Smallman & Williams (2019) https://doi.org/10.5194/gmd-12-2227-2019.
  ! This code is based on that created by A. A. Bloom (UoE, now at JPL, USA).
  ! Subsequent modifications by:
  ! T. L. Smallman (University of Edinburgh, t.l.smallman@ed.ac.uk)
  ! J. F. Exbrayat (University of Edinburgh)
  ! See function/subroutine specific comments for exceptions and contributors
  !!!!!!!!!!!

  ! make all private
  private

  ! explicit publics
  public:: CARBON_MODEL        &
           ,vsmall              &
           ,arrhenious          &
           ,opt_max_scaling     &
           ,acm_gpp_stage_1     &
           ,acm_gpp_stage_2     &
           ,calculate_stomatal_conductance       &
           ,meteorological_constants  &
           ,calculate_radiation_balance &
           ,calculate_daylength           &
           ,freeze              &
           ,calculate_Rtot      &
           ,calculate_aerodynamic_conductance &
           ,linear_model_gradient &
           ,seconds_per_day  &
           ,seconds_per_hour &
           ,seconds_per_step &
           ,fine_root_biomass&
           ,root_biomass     &
           ,root_reach       &
           ,root_depth_time  &
           ,min_root         &
           ,max_depth        &
           ,root_k           &
           ,gs_demand_supply_ratio        &
           ,gs_total_canopy               &
           ,gb_total_canopy               &
           ,canopy_par_MJday_time         &
           ,canopy_par_MJday &
           ,top_soil_depth   &
           ,soil_depth       &
           ,previous_depth   &
           ,nos_root_layers  &
           ,deltat_1         &
           ,total_water_flux              &
           ,water_flux_mmolH2Om2s         &
           ,layer_thickness  &
           ,min_layer        &
           ,nos_soil_layers  &
           ,soil_frac_clay   &
           ,soil_frac_sand   &
           ,meant            &
           ,cica_time        &
           ,convert_ms1_mol_1 &
           ,aerodynamic_conductance &
           ,stomatal_conductance &
           ,potential_conductance &
           ,avN              &
           ,iWUE             &
           ,NUE              &
           ,pn_max_temp      &
           ,pn_opt_temp      &
           ,pn_kurtosis      &
           ,co2_half_sat     &
           ,co2_comp_point   &
           ,minlwp           &
           ,leafT            &
           ,mint             &
           ,maxt             &
           ,swrad            &
           ,co2              &
           ,doy              &
           ,wind_spd         &
           ,vpd_kPa          &
           ,lai              &
           ,dayl_seconds_1   &
           ,dayl_seconds     &
           ,dayl_hours       &
           ,dayl_hours_fraction &
           ,Rg_from_labile            &
           ,harvest_residue_to_litter &
           ,harvest_residue_to_woodlitter&
           ,harvest_residue_to_som    &
           ,harvest_extracted_litter      &
           ,harvest_extracted_woodlitter  &
           ,harvest_extracted_som         &
           ,harvest_extracted_labile      &
           ,harvest_extracted_foliar      &
           ,harvest_extracted_roots       &
           ,harvest_extracted_wood        &
           ,harvest_residue_labile        &
           ,harvest_residue_foliar        &
           ,harvest_residue_roots         &
           ,harvest_residue_wood          &
           ,fire_litter_labile              &
           ,fire_litter_foliar              &
           ,fire_litter_roots               &
           ,fire_litter_wood                &
           ,fire_litter_litter              &
           ,fire_litter_woodlitter             &
           ,fire_litter_som                 &
           ,fire_emiss_labile              &
           ,fire_emiss_foliar              &
           ,fire_emiss_roots               &
           ,fire_emiss_wood                &
           ,fire_emiss_litter              &
           ,fire_emiss_woodlitter             &
           ,fire_emiss_som                 &
           ,fire_residue_to_litter    &
           ,fire_residue_to_woodlitter   &
           ,fire_residue_to_som       &
           ,itemp, ivpd, iphoto&
           ,dim_1, dim_2      &
           ,nos_trees        &
           ,nos_inputs       &
           ,leftDaughter     &
           ,rightDaughter    &
           ,nodestatus       &
           ,xbestsplit       &
           ,nodepred         &
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

  double precision, allocatable, dimension(:,:):: leftDaughter, & ! left daughter for forest
                                                  rightDaughter, & ! right daughter for forets
                                                     nodestatus, & ! nodestatus for forests
                                                     xbestsplit, & ! for forest
                                                       nodepred, & ! prediction value for each tree
                                                        bestvar    ! for randomForests
  !!!!!!!!!
  ! Parameters
  !!!!!!!!!

  ! useful technical parameters
  double precision, parameter:: xacc = 1d-4        & ! accuracy parameter for zbrent bisection proceedure  ! 0.0001
                              ,vsmall = tiny(0d0)*1e3

  integer, parameter:: nos_root_layers = 3, nos_soil_layers = nos_root_layers+1
  double precision, parameter:: pi = 3.1415927d0,    &
                               pi_1 = 0.3183099d0,    & ! pi**(-1d0)
                                pi2 = 9.869604d0,     & ! pi**2d0
                             two_pi = 6.283185d0,     & ! pi*2d0
                         deg_to_rad = 0.01745329d0,   & ! pi/180d0
                sin_dayl_deg_to_rad = 0.3979486d0,    & ! sin( 23.45d0*deg_to_rad )
                            gravity = 9.8067d0,       & ! acceleration due to gravity, ms-1
                              boltz = 5.670400d-8,    & ! Boltzmann constant (W.m-2.K-4)
                         emissivity = 0.96d0,         &
                        emiss_boltz = 5.443584d-08,   & ! emissivity*boltz
                    sw_par_fraction = 0.5d0,          & ! fraction of short-wave radiation which is PAR
                             freeze = 273.15d0,       &
                  gs_H2Ommol_CO2mol = 0.001646259d0, & ! The ratio of H20:CO2 diffusion for gs (Jones appendix 2)
              gs_H2Ommol_CO2mol_day = 142.2368d0,     & ! The ratio of H20:CO2 diffusion for gs, including seconds per day correction
                         gs_H2O_CO2 = 1.646259d0,     & ! The ratio of H20:CO2 diffusion for gs (Jones appendix 2)
                       gs_H2O_CO2_1 = 0.6074378d0,    & ! gs_H2O_CO2 ** (-1d0)
                         gb_H2O_CO2 = 1.37d0,         & ! The ratio of H20:CO2 diffusion for gb (Jones appendix 2)
            partial_molar_vol_water = 18.05d-6,       & ! partial molar volume of water, m3 mol-1 at 20C
                     mol_to_g_water = 18d0,           & ! molecular mass of water
                   mmol_to_kg_water = 1.8d-5,         & ! milli mole conversion to kg
                       mol_to_g_co2 = 12d0,           & ! molecular mass of CO2 (g)
                         umol_to_gC = 1.2d-5,         & ! conversion of umolC -> gC
                         gC_to_umol = 83333.33d0,     & ! conversion of gC -> umolC; umol_to_gC**(-1d0)
                       g_to_mol_co2 = 0.08333333d0,   &
!snowscheme       density_of_water = 998.9d0,         & ! density of  ! water kg.m-3
                     gas_constant_d = 287.04d0,       & ! gas constant for dry air (J.K-1.mol-1)
                               Rcon = 8.3144d0,       & ! Universal gas constant (J.K-1.mol-1)
                          vonkarman = 0.41d0,         & ! von Karman's constant
                        vonkarman_1 = 2.439024d0,     & ! 1/von Karman's constant
                        vonkarman_2 = 0.1681d0,       & ! von Karman's constant^2
                              cpair = 1004.6d0          ! Specific heat capacity of air; used in energy balance J.kg-1.K-1

  ! photosynthesis/respiration parameters
  double precision, parameter :: &
                        Rg_fraction = 0.21875d0,    & ! fraction of C allocation towards each pool
                                                      ! lost as growth respiration
                                                      ! (i.e. 0.28 .eq. xNPP)
                    one_Rg_fraction = 1d0-Rg_fraction

  ! hydraulic parameters
  double precision, parameter :: &
                         tortuosity = 2.5d0,        & ! tortuosity
                             gplant = 5d0,          & ! plant hydraulic conductivity (mmol m-1 s-1 MPa-1)
                        root_resist = 25d0,         & ! Root resistivity (MPa s g mmol−1 H2O)
                          max_depth = 2d0,          & ! max root depth (m)
                             root_k = 200d0,        & ! root biomass needed to reach 50% depth (gbiomass/m2); fine root only value = 100
                        root_radius = 0.00029d0,    & ! root radius (m) Bonen et al 2014 = 0.00029
                                                      !                 Williams et al 1996 = 0.0001
                      root_radius_1 = root_radius**(-1d0), &
                root_cross_sec_area = 3.141593d-08, & ! root cross sectional area (m2)
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
                          min_layer = 0.03d0,       & ! minimum thickness of the third rooting layer (m)
                        soil_roughl = 0.05d0,       & ! soil roughness length (m)
                     top_soil_depth = 0.30d0,       & ! thickness of the top soil layer (m)
                           min_root = 5d0,          & ! minimum root biomass (gBiomass.m-2)
                            min_lai = 0.1d0           ! minimum LAI assumed for aerodynamic conductance calculations (m2/m2)

  ! timing parameters
  double precision, parameter :: &
                   seconds_per_hour = 3600d0,       & ! Number of seconds per hour
                    seconds_per_day = 86400d0,      & ! Number of seconds per day
                  seconds_per_day_1 = 1.157407d-05    ! Inverse of seconds per day

  ! ACM-GPP-ET parameters
  double precision, parameter :: &
                   ! Assumption that photosythesis will be limited by Jmax temperature response
                   pn_max_temp = 57.05d0,       & ! Maximum daily max temperature for photosynthesis (oC)
                   pn_min_temp = -1d6,          & ! Minimum daily max temperature for photosynthesis (oC)
                   pn_opt_temp = 30d0,          & ! Optimum daily max temperature for photosynthesis (oC)
                   pn_kurtosis = 0.172d0,       & ! Kurtosis of Jmax temperature response
                ko_half_sat_25C = 157.46892d0,  & ! photorespiration O2 half sat(mmolO2/mol), achieved at 25oC
           ko_half_sat_gradient = 14.93643d0,   & ! photorespiration O2 half sat gradient
                kc_half_sat_25C = 319.58548d0,  & ! carboxylation CO2 half sat (umolCO2/mol), achieved at 25oC
           kc_half_sat_gradient = 24.72297d0,   & ! carboxylation CO2 half sat gradient
                co2comp_sat_25C = 36.839214d0,  & ! carboxylation CO2 compensation point(umolCO2/mol), saturation
               co2comp_gradient = 9.734371d0,   & ! carboxylation CO2 comp point, achieved at oC
                                                  ! Each of these are temperature sensitivty
                            e0 = 3.2d+00,       & ! Quantum yield (gC/MJ/m2/day PAR), SPA apparent yield
                minlwp_default = -1.808224d+00,  & ! minimum leaf water potential (MPa). NOTE: actual SPA = -2 MPa
      soil_iso_to_net_coef_LAI = -2.717467d+00,  & ! Coefficient relating soil isothermal net radiation to net.
!             orig             iWUE = 1.8d-7,        & ! Intrinsic water use efficiency (gC/mmolH2O-1/m2leaf/s-1)
                          iWUE = 1.5d-2,        & ! Intrinsic water use efficiency (umolC/mmolH2O-1/m2leaf/s-1)
         soil_swrad_absorption = 9.989852d-01,  & ! Fraction of SW rad absorbed by soil
         max_lai_lwrad_release = 9.516639d-01,  & ! 1-Max fraction of LW emitted from canopy to be released
        lai_half_lwrad_release = 4.693329d+00,  & ! LAI at which LW emitted from canopy to be released at 50 %
       soil_iso_to_net_coef_SW = -3.500964d-02,  & ! Coefficient relating soil isothermal net radiation to net.
         soil_iso_to_net_const = 3.455772d+00,  & ! Constant relating soil isothermal net radiation to net
           max_par_transmitted = 1.628077d-01,  & ! Max fraction of canopy incident PAR transmitted to soil
           max_nir_transmitted = 2.793660d-01,  & ! Max fraction of canopy incident NIR transmitted to soil
             max_par_reflected = 1.629133d-01,  & ! Max fraction of canopy incident PAR reflected to sky
             max_nir_reflected = 4.284365d-01,  & ! Max fraction of canopy incident NIR reflected to sky
     canopy_iso_to_net_coef_SW = 1.480105d-02,  & ! Coefficient relating SW to the adjustment between isothermal and net LW
       canopy_iso_to_net_const = 3.753067d-03,  & ! Constant relating canopy isothermal net radiation to net
    canopy_iso_to_net_coef_LAI = 2.455582d+00     ! Coefficient relating LAI to the adjustment between isothermal and net LW

  double precision:: minlwp = minlwp_default

  !!!!!!!!!
  ! Module level variables
  !!!!!!!!!

  ! management and gsi related values
  integer:: gsi_lag_remembered
  ! local variables for GSI phenology model
  double precision:: Tfac, Photofac, VPDfac & ! oC, seconds, Pa
                     ,delta_gsi, tmp, tmp_lai, gradient         &
                     ,fol_turn_crit

  double precision, allocatable, dimension(:):: Rg_from_labile,                &
                                              itemp, ivpd, iphoto, &
                                         gs_demand_supply_ratio, & ! actual:potential stomatal conductance
                                                gs_total_canopy, & ! stomatal conductance (mmolH2O/m2ground/day)
                                                gb_total_canopy, & ! boundary conductance (mmolH2O/m2ground/day)
                                          canopy_par_MJday_time, & ! Absorbed PAR by canopy (MJ/m2ground/day)
                                                root_depth_time, & ! rooting depth (m)
                                      harvest_residue_to_litter, &
                                         harvest_residue_to_som, &
                                  harvest_residue_to_woodlitter, &
                                       harvest_extracted_litter, &
                                   harvest_extracted_woodlitter, &
                                          harvest_extracted_som, &
                                       harvest_extracted_labile, &
                                       harvest_extracted_foliar, &
                                        harvest_extracted_roots, &
                                         harvest_extracted_wood, &
                                         harvest_residue_labile, &
                                         harvest_residue_foliar, &
                                          harvest_residue_roots, &
                                           harvest_residue_wood, &
                                             fire_litter_labile, &
                                             fire_litter_foliar, &
                                              fire_litter_roots, &
                                               fire_litter_wood, &
                                             fire_litter_litter, &
                                         fire_litter_woodlitter, &
                                                fire_litter_som, &
                                              fire_emiss_labile, &
                                              fire_emiss_foliar, &
                                               fire_emiss_roots, &
                                                fire_emiss_wood, &
                                              fire_emiss_litter, &
                                          fire_emiss_woodlitter, &
                                                 fire_emiss_som, &
                                         fire_residue_to_litter, &
                                     fire_residue_to_woodlitter, &
                                            fire_residue_to_som, &
                                      tmp_x, tmp_m, gsi_history

  ! Phenological choices
  ! See source below for details of these variables
  double precision:: root_cost, root_life,       &
                      leaf_cost, leaf_life

  ! hydraulic model variables
  integer:: water_retention_pass, soil_layer
  double precision, dimension(nos_soil_layers+1):: layer_thickness    ! thickness of soil layers (m)
  double precision, dimension(nos_soil_layers):: soil_frac_clay, soil_frac_sand  ! clay and soil fractions of soil
  double precision, dimension(nos_root_layers):: uptake_fraction, & ! fraction of water uptake from each root layer
                                                           demand, & ! maximum potential canopy hydraulic demand
                                            water_flux_mmolH2Om2s    ! potential transpiration flux (mmolH2O.m-2.s-1)

  double precision:: root_reach, root_biomass, &
                             fine_root_biomass, & ! root depth, coarse+fine, and fine root biomass
                              total_water_flux, & ! potential transpiration flux (kgH2O.m-2.day-1)
                                    soil_depth, &
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
                             convert_ms1_mol_1, & ! Conversion ratio for m.s-1 -> mol.m-2.s-1
                            convert_ms1_mmol_1, & ! Conversion ratio for m/s -> mmol/m2/s
                                        lambda, & ! latent heat of vapourisation (J.kg-1)
                                         psych, & ! psychrometric constant (kPa K-1)
                                         slope, & ! Rate of change of saturation vapour pressure
                                                  ! with temperature (kPa.K-1)
                        water_vapour_diffusion, & ! Water vapour diffusion coefficient in (m2/s)
                             dynamic_viscosity, & ! dynamic viscosity (kg.m-2.s-1)
                           kinematic_viscosity    ! kinematic viscosity (m2.s-1)

  ! Module level variables for ACM_GPP_ET parameters
  double precision ::   delta_gs, & ! day length corrected gs increment mmolH2O/m2/dayl
                             avN, & ! average foliar N (gN/m2)
                            ceff, & ! Canopy photosynthetic efficiency (avN*NUE)
                       iWUE_step, & ! Intrinsic water use efficiency for that day (gC/m2leaf/dayl/mmolH2Ogs)
                             NUE, & ! Photosynthetic nitrogen use efficiency at optimum temperature (oC)
                                    ! ,unlimited by CO2, light and photoperiod (gC/gN/m2leaf/day)
metabolic_limited_photosynthesis, & ! temperature, leaf area and foliar N limiterd photosynthesis (gC/m2/day)
    light_limited_photosynthesis, & ! light limited photosynthesis (gC/m2/day)
                              ci, & ! Internal CO2 concentration (ppm)
                          gb_mol, & ! Canopy boundary layer conductance (molCO2/m2/day)
                        rb_mol_1, & ! Canopy boundary layer resistance (day/m2/molCO2)
                    co2_half_sat, & ! CO2 at which photosynthesis is 50 % of maximum (ppm)
                  co2_comp_point    ! CO2 at which photosynthesis > 0 (ppm)

  ! Module level variables for step specific met drivers
  double precision:: mint, & ! minimum temperature (oC)
                      maxt, & ! maximum temperature (oC)
                     leafT, & ! leaf temperature (oC)
                     swrad, & ! incoming short wave radiation (MJ/m2/day)
                       co2, & ! CO2 (ppm)
                       doy, & ! Day of year
                  wind_spd, & ! wind speed (m/s)
                   vpd_kpa, & ! Vapour pressure deficit (kPa)
 leaf_canopy_light_scaling, & ! approximate scaling factor from leaf to canopy for gpp and gs
  leaf_canopy_wind_scaling, & ! approximate scaling factor from leaf to canopy for gb
                     lai_1, & ! inverse of LAI
                       lai    ! leaf area index (m2/m2)

  ! Module level varoables for step specific timing information
  double precision:: cos_solar_zenith_angle, &
                            seconds_per_step, & !
                                dayl_seconds, & ! day length in seconds
                          mean_days_per_step, &
                              dayl_seconds_1, &
                         dayl_hours_fraction, &
                                  dayl_hours    ! day length in hours

  double precision, dimension(:), allocatable ::    deltat_1, & ! inverse of decimal days
                                                   cica_time, & ! Internal vs ambient CO2 concentrations
                                                  meant_time

contains
  !
  !--------------------------------------------------------------------
  !
  subroutine CARBON_MODEL(start, finish, met, pars, deltat, nodays, lat, lai_out, NEE_out, FLUXES, POOLS &
                         ,nopars, nomet, nopools, nofluxes, GPP_out)

    ! The Data Assimilation Linked Ecosystem Carbon-Growing Season
    ! Index using DeltaGSI to determine FOLiar phenology-Forest Rotation (DALEC_GSI_DFOL_FR) model.
    ! The subroutine calls the Aggregated Canopy Model to simulate GPP and
    ! partitions between various ecosystem carbon pools. These pools are
    ! subject to turnovers/decompostion resulting in ecosystem phenology and fluxes of CO2

    implicit none

    ! declare input variables
    integer, intent(in):: start &
                         ,finish &
                         ,nopars & ! number of paremeters in vector
                          ,nomet & ! number of meteorological fields
                       ,nofluxes & ! number of model fluxes
                        ,nopools & ! number of model pools
                         ,nodays   ! number of days in simulation

    double precision, intent(in):: met(nomet, nodays) & ! met drivers
                                      ,deltat(nodays) & ! time step in decimal days
                                        ,pars(nopars) & ! number of parameters
                                                ,lat    ! site latitude (degrees)

    double precision, dimension(nodays), intent(inout):: lai_out & ! leaf area index
                                                         ,GPP_out & ! Gross primary productivity
                                                         ,NEE_out   ! net ecosystem exchange of CO2

    double precision, dimension((nodays+1), nopools), intent(inout):: POOLS  ! vector of ecosystem pools
    double precision, dimension(nodays, nofluxes), intent(inout):: FLUXES  ! vector of ecosystem fluxes

    double precision:: Rtot  ! Total hydraulic resistance (MPa.s-1.m-2.mmol-1)
    integer:: p, f, nxp, n, test, m

    ! local fire related variables
    double precision:: burnt_area           &
                       ,combust_eff(7) & ! combustion efficiency
                       ,rfac(7)          ! resilience factor

    ! local deforestation related variables
    double precision, dimension(5):: post_harvest_burn      & ! how much burning to occur after
                                     ,foliage_frac_res       &
                                     ,roots_frac_res         &
                                     ,rootcr_frac_res        &
                                     ,stem_frac_res          &
                                     ,roots_frac_removal     &
                                     ,rootcr_frac_removal    &
                                     ,Crootcr_part           &
                                     ,soil_loss_frac

    double precision:: labile_loss, foliar_loss      &
                       ,roots_loss, wood_loss         &
                       ,rootcr_loss, stem_loss        &
                       ,labile_residue, foliar_residue&
                       ,roots_residue, wood_residue   &
                       ,C_total, labile_frac_res      &
                       ,labile_frac_removal        &
                       ,Cstem, Crootcr, stem_residue   &
                       ,coarse_root_residue          &
                       ,soil_loss_with_roots

    ! variables for phenology model update/adjustments
    double precision:: C_invest, &
                        lai_save, &
                  canopy_lw_save, &
                    soil_lw_save, &
                  canopy_sw_save, &
                 canopy_par_save, &
                    soil_sw_save, &
                         gs_save, &
                         ga_save

    integer:: reforest_day, harvest_management, restocking_lag, gsi_lag, interval

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
    ! 12th 21 day average VPD (kPa)
    ! 13th Forest management practice to accompany any clearing
    ! 14th avg daily temperature (oC)
    ! 15th avg daily wind speed (m.s-1)
    ! 16th vapour pressure deficit (Pa)

    ! POOLS are:
    ! 1 = labile      (p18)
    ! 2 = foliar      (p19)
    ! 3 = fine root   (p20)
    ! 4 = wood        (p21)
    ! 5 = litter (fol+fine root) (p22)
    ! 6 = som         (p23)
    ! 7 = wood litter (p37)

    ! p(30) = labile replanting
    ! p(31) = foliar replanting
    ! p(32) = fine root replanting
    ! p(33) = wood replanting

    ! FLUXES are:
    ! 1 = GPP
    ! 2 = temprate
    ! 3 = respiration_auto
    ! 4 = respiration het litwood
    ! 5 = labile production
    ! 6 = root production
    ! 7 = wood production
    ! 8 = labile production
    ! 9 = leaffall factor
    ! 10 = leaf litter production
    ! 11 = woodlitter production
    ! 12 = rootlitter production
    ! 13 = respiration het litter (fol+fine root)
    ! 14 = respiration het som
    ! 15 = litter2som
    ! 16 = labrelease factor
    ! 17 = carbon flux due to fire
    ! 18 = growing season index
    ! 19 = litwood turnover to litter
    ! 20 = C extracted as harvest
    ! 21 = NOT IN USE
    ! 22 = NOT IN USE
    ! 23 = NOT IN USE
    ! 24 = NOT IN USE

    ! PARAMETERS
    ! 31 process parameters; 7 C pool initial conditions

    ! p(1) = decomposition efficiency (fraction to som)
    ! p(2) = RmGPP fraction
    ! p(3) = Background leaf turnover
    ! p(4) = Max labile turnover to roots (fraction)
    ! p(5) = Max leaf turnover (GSI; fraction)
    ! p(6) = Turnover rate of wood (fraction)
    ! p(7) = Turnover rate of roots (fraction)
    ! p(8) = Litter turnover rate to heterotrophic respiration (fraction)
    ! p(9) = SOM turnover rate to heterotrophic respiration (fraction)
    ! p(10) = Exponential coefficient for temperature response for heterotrophic respiration
    ! p(11) = Average foliar nitrogen content (log10(gN/m2))
    ! p(12) = Max labile turnover to leaves (GSI; fraction)
    ! p(13) = Max labile turnover to wood (fraction)
    ! p(14) = Min temp threshold (GSI; Kelvin)
    ! p(15) = Max temp threshold (GSI; Kelvin)
    ! p(16) = Min photoperiod threshold (GSI; seconds)
    ! p(17) = Leaf Mass per unit Area (gC.m-2)
    ! p(18) = Initial labile pool (gC/m2)
    ! p(19) = Initial foliage pool (gC/m2)
    ! p(20) = Initial root pool (gC/m2)
    ! p(21) = Initial wood pool (gC/m2)
    ! p(22) = Initial litter pool (gC/m2)
    ! p(23) = Initial som pool (gC/m2)
    ! p(24) = Max photoperiod threshold (GSI; seconds)
    ! p(25) = Min vapour pressure deficit threshold (GSI; Pa)
    ! p(26) = Max vapour pressure deficit (GSI; Pa)
    ! p(27) = GPP return on new Cfol investment (gCperGPP per gCnewfol)
    ! p(28) = Initial GSI value
    ! p(29) = Fraction of Cwood which is Ccoarseroot
    ! p(35) = litwood turnover fraction (fraction)
    ! p(36) = Nitrogen use efficiency (gC/gN/day)
    ! p(37) = Initial litwood pool (gC/m2)

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
    !real:: begin, done, f1 = 0, f2 = 0, f3 = 0, f4 = 0, f5 = 0
    !real:: Rtot_times = 0, aero_time = 0, soilwater_time = 0, acm_et_time = 0
    !call cpu_time(begin)
    !call cpu_time(finish)

    ! Reset all POOLS and FLUXES to prevent precision errors
    FLUXES = 0d0; POOLS = 0d0

    ! load ACM-GPP-ET parameters
    NUE = pars(36)      ! Photosynthetic nitrogen use efficiency at optimum temperature (oC)
                        ! ,unlimited by CO2, light and
                        ! photoperiod (gC/gN/m2leaf/day)
    avN = 10d0**pars(11)  ! foliar N gN/m2
!    ceff = avN*NUE       ! canopy efficiency, used to avoid what in most cases is a reductance multiplication
!                         ! NOTE: must be updated any time NUE or avN changes

    ! if either of our disturbance drivers indicate disturbance will occur then
    ! set up these components
    if (maxval(met(8, :)) > 0d0 .or. maxval(met(9, :)) > 0d0) then

        ! initial values for deforestation variables
        labile_loss = 0d0    ; foliar_loss = 0d0
        roots_loss = 0d0     ; wood_loss = 0d0
        rootcr_loss = 0d0    ; stem_loss = 0d0
        labile_residue = 0d0; foliar_residue = 0d0
        roots_residue = 0d0  ; wood_residue = 0d0
        stem_residue = 0d0
        reforest_day = 0
        soil_loss_with_roots = 0d0
        coarse_root_residue = 0d0
        post_harvest_burn = 0d0

        ! now load the hardcoded forest management parameters into their locations

        ! Deforestation process functions in a sequenctial way.
        ! Thus, the pool_loss is first determined as a function of met(n, 8) and
        ! for fine and coarse roots whether this felling is associated with a mechanical
        ! removal from the ground. As the canopy and stem is removed (along with a proportion of labile)
        ! fine and coarse roots may subsequently undergo mortality from which they do not recover
        ! but allows for management activities such as grazing, mowing and coppice.
        ! The pool_loss is then partitioned between the material which is left within the system
        ! as a residue and thus direcly placed within one of the dead organic matter pools.

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
        Crootcr_part(1) = 0.32d0  ! Coarse roots (Adegbidi et al 2005; 
        ! Csom loss due to phyical removal with roots
        ! Morison et al (2012) Forestry Commission Research Note
        soil_loss_frac(1) = 0.02d0  ! actually between 1-3 %
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
        Crootcr_part(2) = 0.32d0  ! Coarse roots (Adegbidi et al 2005; 
        ! Csom loss due to phyical removal with roots
        ! Morison et al (2012) Forestry Commission Research Note
        soil_loss_frac(2) = 0.02d0  ! actually between 1-3 %
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
        Crootcr_part(3) = 0.32d0  ! Coarse roots (Adegbidi et al 2005; 
        ! Csom loss due to phyical removal with roots
        ! Morison et al (2012) Forestry Commission Research Note
        soil_loss_frac(3) = 0.02d0  ! actually between 1-3 %
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
        Crootcr_part(4) = 0.32d0  ! Coarse roots (Adegbidi et al 2005; 
        ! Csom loss due to phyical removal with roots
        ! Morison et al (2012) Forestry Commission Research Note
        soil_loss_frac(4) = 0.02d0  ! actually between 1-3 %
        ! was the forest burned after deforestation (0-1)
        ! NOTE: that we refer here to the fraction of the cleared land to be burned
        post_harvest_burn(4) = 0d0

        !## Scenario 5 (grassland grazing/cutting)
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
        Crootcr_part(5) = 0.32d0  ! Coarse roots (Adegbidi et al 2005; 
        ! Csom loss due to phyical removal with roots
        ! Morison et al (2012) Forestry Commission Research Note
        soil_loss_frac(5) = 0d0  ! actually between 1-3 %
        ! was the forest burned after deforestation (0-1)
        ! NOTE: that we refer here to the fraction of the cleared land to be burned
        post_harvest_burn(5) = 0d0

        ! Override all paritioning parameters with those coming from
        ! CARDAMOM
        Crootcr_part = pars(29)

        ! Declare combustion efficiency (labile, foliar, roots, wood, litter, soil, woodlitter)
        combust_eff(1) = 0.1d0; combust_eff(2) = 0.9d0
        combust_eff(3) = 0.1d0; combust_eff(4) = 0.1d0
        combust_eff(5) = 0.7d0; combust_eff(6) = 0.01d0
        combust_eff(7) = 0.7d0
        ! Resilience factor for non-combusted tissue
        rfac = 0.5d0; rfac(5) = 0.1d0; rfac(6) = 0d0; rfac(7) = 0.1d0

    end if  ! disturbance ?

    ! JFE added 4 May 2018-define fire constants
    ! Update fire parameters derived from
    ! Yin et al., (2020), doi: 10.1038/s414647-020-15852-2
    ! Subsequently expanded by T. L. Smallman & Mat Williams (UoE, 03/09/2021)
    ! to provide specific CC for litter and wood litter.
    ! NOTE: changes also result in the addition of further EDCs

    ! Assign proposed resilience factor
    rfac(1:4) = pars(38)
    rfac(5) = 0.1d0; rfac(6) = 0d0; rfac(7) = 0.1d0
    ! Assign combustion completeness to foliage
    combust_eff(2) = pars(39)  ! foliage
    ! Assign combustion completeness to non-photosynthetic
    combust_eff(1) = pars(40); combust_eff(3) = pars(40); combust_eff(4) = pars(40)
    combust_eff(6) = pars(41)  ! soil
    ! derived values for litter and wood litter
    combust_eff(5) = pars(42); combust_eff(7) = pars(43)

    ! assigning initial conditions
    POOLS(1, 1) = pars(18)
    POOLS(1, 2) = pars(19)
    POOLS(1, 3) = pars(20)
    POOLS(1, 4) = pars(21)
    POOLS(1, 5) = pars(22)
    POOLS(1, 6) = pars(23)
    POOLS(1, 7) = pars(37)

    if (.not.allocated(harvest_residue_to_litter)) then
        allocate(harvest_residue_to_litter(nodays),    &
                 harvest_residue_to_som(nodays),       &
                 harvest_residue_to_woodlitter(nodays), &
                 harvest_residue_labile(nodays),       &
                 harvest_residue_foliar(nodays),       &
                 harvest_residue_roots(nodays),        &
                 harvest_residue_wood(nodays),         &
                 harvest_extracted_litter(nodays),     &
                 harvest_extracted_som(nodays),        &
                 harvest_extracted_woodlitter(nodays), &
                 harvest_extracted_labile(nodays),     &
                 harvest_extracted_foliar(nodays),     &
                 harvest_extracted_roots(nodays),      &
                 harvest_extracted_wood(nodays),       &
                 fire_emiss_labile(nodays),            &
                 fire_emiss_foliar(nodays),            &
                 fire_emiss_roots(nodays),             &
                 fire_emiss_wood(nodays),              &
                 fire_emiss_litter(nodays),            &
                 fire_emiss_woodlitter(nodays),        &
                 fire_emiss_som(nodays),               &
                 fire_litter_labile(nodays),           &
                 fire_litter_foliar(nodays),           &
                 fire_litter_roots(nodays),            &
                 fire_litter_wood(nodays),             &
                 fire_litter_litter(nodays),           &
                 fire_litter_woodlitter(nodays),       &
                 fire_litter_som(nodays),              &
                 fire_residue_to_litter(nodays),       &
                 fire_residue_to_woodlitter(nodays),   &
                 fire_residue_to_som(nodays))
    endif

    ! Reset harvest residue
    harvest_residue_to_litter = 0d0; harvest_residue_to_woodlitter = 0d0
    harvest_residue_to_som = 0d0
    ! Reset harvest residue
    harvest_residue_labile = 0d0       ; harvest_residue_foliar = 0d0
    harvest_residue_roots = 0d0        ; harvest_residue_wood = 0d0
    ! Reset harvest extracted
    harvest_extracted_labile = 0d0       ; harvest_extracted_foliar = 0d0
    harvest_extracted_roots = 0d0        ; harvest_extracted_wood = 0d0
    harvest_extracted_litter = 0d0       ; harvest_extracted_som = 0d0
    harvest_extracted_woodlitter = 0d0
    ! Reset fire loss to litter
    fire_litter_labile = 0d0; fire_litter_foliar = 0d0; fire_litter_roots = 0d0
    fire_litter_wood = 0d0   ; fire_litter_litter = 0d0; fire_litter_woodlitter = 0d0
    fire_litter_som = 0d0
    ! Reset fire loss to combustion and emissions
    fire_emiss_labile = 0d0; fire_emiss_foliar = 0d0; fire_emiss_roots = 0d0
    fire_emiss_wood = 0d0   ; fire_emiss_litter = 0d0; fire_emiss_woodlitter = 0d0
    fire_emiss_som = 0d0
    ! Reset fire residue
    fire_residue_to_litter = 0d0; fire_residue_to_woodlitter = 0d0; fire_residue_to_som = 0d0

    ! SHOULD TURN THIS INTO A SUBROUTINE CALL AS COMMON TO BOTH DEFAULT AND CROPS
    if (.not.allocated(deltat_1)) then

        allocate(deltat_1(nodays), meant_time(nodays), Rg_from_labile(nodays), &
                 gs_demand_supply_ratio(nodays), canopy_par_MJday_time(nodays), &
                 gs_total_canopy(nodays), gb_total_canopy(nodays), cica_time(nodays), &
                 root_depth_time(nodays))

        ! Inverse of time step (days-1) to avoid divisions
        deltat_1 = deltat**(-1d0)
        ! Meant time step temperature
        meant_time = (met(2, 1:nodays)+met(3, 1:nodays)) * 0.5d0
        ! mean days per step
        mean_days_per_step = sum(deltat) / dble(nodays)

        ! Calculate timing components needed for GSI/NCE gradient calculations
        gsi_lag_remembered = max(2, nint(21d0/mean_days_per_step))
        allocate(tmp_x(gsi_lag_remembered), gsi_history(gsi_lag_remembered))
        do f = 1, gsi_lag_remembered
           tmp_x(f) = dble(f) * mean_days_per_step
        end do

    else

    endif  ! has SWP already been determined?

    ! load some needed module level values
    lai = POOLS(1, 2)/pars(17)
    mint = met(2, 1)  ! minimum temperature (oC)
    maxt = met(3, 1)  ! maximum temperature (oC)
    swrad = met(4, 1)  ! incoming short wave radiation (MJ/m2/day)
    co2 = met(5, 1)   ! CO2 (ppm)
    doy = met(6, 1)   ! Day of year
    wind_spd = met(15, 1)  ! wind speed (m/s)
    vpd_kPa = met(16, 1)*1d-3  ! vapour pressure deficit (Pa->kPa)
    meant = meant_time(1)
    leafT = (maxt*0.75d0) + (mint*0.25d0)     ! initial canopy temperature (oC)
    seconds_per_step = deltat(1) * seconds_per_day

    ! initialise root reach based on initial conditions
    fine_root_biomass = max(min_root, POOLS(1, 3)*2d0)
    root_biomass = fine_root_biomass+max(min_root, POOLS(1, 4)*pars(29)*2d0)
    ! needed to initialise soils
    call calculate_Rtot(Rtot)
    root_depth_time(1) = root_reach
    ! assign climate sensitivities
    ! assign our starting value
    gsi_history = pars(28)
    gsi_lag = gsi_lag_remembered  ! added to prevent loss from memory
    fol_turn_crit = pars(34)

    !
    ! Begin looping through each time step
    !

    do n = start, finish

       !!!!!!!!!!
       ! assign drivers and update some prognostic variables
       !!!!!!!!!!

       ! Incoming drivers
       mint = met(2, n)  ! minimum temperature (oC)
       maxt = met(3, n)  ! maximum temperature (oC)
       swrad = met(4, n)  ! incoming short wave radiation (MJ/m2/day)
       co2 = met(5, n)   ! CO2 (ppm)
       doy = met(6, n)   ! Day of year
       meant = meant_time(n)  ! mean air temperature (oC)
       leafT = (meant+maxt) * 0.5d0  ! estimate mean daytime air temperature (oC)
       wind_spd = met(15, n)  ! wind speed (m/s)
       vpd_kPa = met(16, n)*1d-3  ! Vapour pressure deficit (Pa -> kPa)

       ! states needed for module variables
       lai_out(n) = POOLS(n, 2)/pars(17)
       lai = lai_out(n)  ! leaf area index (m2/m2)

       ! extract timing related values
       call calculate_daylength((doy-(deltat(n)*0.5d0)), lat)
       dayl_seconds_1 = dayl_seconds ** (-1d0)
       dayl_hours_fraction = dayl_hours*0.04166667d0  ! 1/24 = 0.04166667
       iWUE_step = iWUE*dayl_hours_fraction
       seconds_per_step = seconds_per_day*deltat(n)

       !!!!!!!!!!
       ! Calculate surface exchange coefficients
       !!!!!!!!!!

       ! calculate some temperature dependent meteorologial properties
       call meteorological_constants(leafT, leafT+freeze, vpd_kPa)
       ! pass variables from memory objects
       convert_ms1_mmol_1 = convert_ms1_mol_1*1d3
       ! calculate aerodynamic using consistent approach with SPA
       call calculate_aerodynamic_conductance
       ! Canopy scale aerodynamic conductance (mmolH2O/m2ground/s)
       gb_total_canopy(n) = aerodynamic_conductance*convert_ms1_mmol_1 * &
                            leaf_canopy_wind_scaling

       !!!!!!!!!!
       ! Determine net shortwave and isothermal longwave energy balance
       !!!!!!!!!!

       call calculate_radiation_balance
       canopy_par_MJday_time(n) = canopy_par_MJday

       !!!!!!!!!!
       ! Calculate physically soil water potential and total hydraulic resistance
       !!!!!!!!!!

       ! calculate the minimum soil & root hydraulic resistance based on total
       ! fine root mass ! *2*2 => *RS*C->Bio
       fine_root_biomass = max(min_root, POOLS(n, 3)*2d0)
       root_biomass = fine_root_biomass+max(min_root, POOLS(n, 4)*pars(29)*2d0)
       call calculate_Rtot(Rtot)
       root_depth_time(n) = root_reach

       ! calculate variables used commonly between ACM_GPP and ACM_ET
       call calculate_stomatal_conductance
       ! Estimate stomatal conductance relative to its minimum/maximum, i.e. how
       ! close are we to maxing out supply (note 0.01 taken from min_gs)
       gs_demand_supply_ratio(n) = (stomatal_conductance  - minimum_conductance) &
                                 / (potential_conductance-minimum_conductance)
       ! Store the canopy level stomatal conductance (mmolH2O/m2ground/s)
       !gs_total_canopy(n) = stomatal_conductance*dayl_seconds_1
       gs_total_canopy(n) = stomatal_conductance

       !!!!!!!!!!
       ! GPP (gC.m-2.day-1)
       !!!!!!!!!!

       if (stomatal_conductance > vsmall) then
           ! Gross primary productivity (gC/m2/day)
           ! Assumes acm_gpp_stage_1 ran as part of stomatal conductance calculation
           FLUXES(n, 1) = acm_gpp_stage_2(stomatal_conductance) * umol_to_gC*dayl_seconds
           cica_time(n) = ci/co2
       else
           FLUXES(n, 1) = 0d0; cica_time(n) = 0d0
       endif

       ! temprate (i.e. temperature modified rate of metabolic activity))
       FLUXES(n, 2) = exp(pars(10)*meant)
       ! (maintenance) autotrophic respiration (gC.m-2.day-1)
       FLUXES(n, 3) = pars(2)*FLUXES(n, 1)
       ! labile production (gC.m-2.day-1)
       FLUXES(n, 5) = (FLUXES(n, 1)-FLUXES(n, 3))*pars(13)
       ! root production (gC.m-2.day-1)
       FLUXES(n, 6) = (FLUXES(n, 1)-FLUXES(n, 3)-FLUXES(n, 5))*pars(4)
       ! wood production
       FLUXES(n, 7) = FLUXES(n, 1)-FLUXES(n, 3)-FLUXES(n, 5)-FLUXES(n, 6)

       ! GSI added to fortran version by TLS 24/11/2014
       ! /* 25/09/14-JFE
       ! Here we calculate the Growing Season Index based on
       ! Jolly et al. A generalized, bioclimatic index to predict foliar
       ! phenology in response to climate Global Change Biology, Volume 11, page 619-632, 
       ! 2005 (doi: 10.1111/j.1365-2486.2005.00930.x)
       ! Stoeckli, R., T. Rutishauser, I. Baker, M. A. Liniger, and A. S.
       ! Denning (2011), A global reanalysis of vegetation phenology, J. Geophys. Res.,
       ! 116, G03020, doi:10.1029/2010JG001545.

       ! It is the product of 3 limiting factors for temperature, photoperiod and
       ! vapour pressure deficit that grow linearly from 0 to 1 between a calibrated
       ! min and max value. Photoperiod, VPD and avgTmax are direct input

       ! temperature limitation, then restrict to 0-1; correction for k-> oC
       Tfac = (met(10, n)-(pars(14)-freeze)) / (pars(15)-pars(14))
       Tfac = min(1d0, max(0d0, Tfac))
       ! photoperiod limitation
       Photofac = (met(11, n)-pars(16)) / (pars(24)-pars(16))
       Photofac = min(1d0, max(0d0, Photofac))
       ! VPD limitation
       VPDfac = 1d0 - ( (met(12, n)-pars(25)) / (pars(26)-pars(25)) )
       VPDfac = min(1d0, max(0d0, VPDfac))

       ! calculate and store the GSI index
       FLUXES(n, 18) = Tfac*VPDfac*Photofac

       ! load lag for linear regression
       gsi_lag = gsi_lag_remembered

       !!!
       ! Estimate of change (i.e. gradient) in the GSI/NCE

       ! Determine GSI/NCE section to have linear regression applied to and
       ! determine the number of values, i.e. the interval
       if (n < gsi_lag) then
           if (n == 1) then
               gsi_history(2) = FLUXES(n, 18)
               interval = 2
           else
               gsi_history(1:n) = FLUXES(1:n, 18)
               interval = n
           endif
       else
           gsi_history(1:gsi_lag) = FLUXES((n-gsi_lag+1):n, 18)
           interval = gsi_lag
       end if
       ! Now calculate the linear gradient
       gradient = linear_model_gradient(tmp_x(1:interval), gsi_history(1:interval), interval)

       ! store lag to keep fresh in memory-yes this is a hack to get around a
       ! memory problem
       gsi_lag_remembered = gsi_lag

       ! first assume that nothing is happening
       FLUXES(n, 9) = 0d0  ! leaf turnover
       FLUXES(n, 16) = 0d0  ! leaf growth

       ! save original values for re-allocation later
       canopy_lw_save = canopy_lwrad_Wm2; soil_lw_save = soil_lwrad_Wm2
       canopy_sw_save = canopy_swrad_MJday; canopy_par_save  = canopy_par_MJday
       soil_sw_save = soil_swrad_MJday; gs_save = stomatal_conductance
       ga_save = aerodynamic_conductance; lai_save = lai

       ! Can we grow?
       if (gradient > fol_turn_crit .and. POOLS(n, 1) > 0d0) then

           ! GSI gradient does not indicate leaf turnover, therefore we will
           ! consider leaf growth

           ! Estimate potential leaf growth rate (from labile)
           FLUXES(n, 16) = pars(12)*FLUXES(n, 18)
           ! Convert that into potential new leaf area
           tmp = POOLS(n, 1)*(1d0-(1d0-FLUXES(n, 16))**deltat(n))*deltat_1(n)
           C_invest = tmp
           ! Update the canopy
           lai = (POOLS(n, 2)+tmp)/pars(17)
           tmp = lai/lai_save
           ! And estimate the canopy conditions for the new leaf area under the
           ! current time step
           aerodynamic_conductance = aerodynamic_conductance*tmp
           stomatal_conductance = stomatal_conductance*tmp
           call calculate_shortwave_balance
           if (lai_save < vsmall) then
               call calculate_aerodynamic_conductance
               call calculate_stomatal_conductance
           endif
           ! And estimate potential GPP
           if (stomatal_conductance > vsmall) then
                call acm_gpp_stage_1
               tmp = acm_gpp_stage_2(stomatal_conductance)
           else
               tmp = 0d0
           endif
           ! Determine if increase in LAI leads to an improvement in GPP greater
           ! than critical value, if not then no labile turnover allowed
           if ( ((tmp-FLUXES(n, 1))/C_invest) < pars(27) ) then
                FLUXES(n, 16) = 0d0
           endif

       endif  ! Growth potential?

       ! Are we losing leaves?
       if (gradient < fol_turn_crit .or. FLUXES(n, 18) < vsmall) then

           ! We are in a decending environment and therefore losing leaf area
           FLUXES(n, 9) = pars(5)*(1d0-FLUXES(n, 18))

       endif  ! losing leaves?
       ! If neither process is occuring assume that there is background turnover
       !if (FLUXES(n, 16) == 0d0 .and. FLUXES(n, 9) == 0d0) FLUXES(n, 9) = pars(3)

       ! restore original value back from memory
       lai = lai_save
       !canopy_lwrad_Wm2 = canopy_lw_save; soil_lwrad_Wm2 = soil_lw_save
       canopy_swrad_MJday = canopy_sw_save; canopy_par_MJday = canopy_par_save
       soil_swrad_MJday = soil_sw_save; stomatal_conductance = gs_save
       aerodynamic_conductance = ga_save

       ! these allocated if post-processing
       if (allocated(itemp)) then
           itemp(n) = Tfac
           ivpd(n) = VPDfac
           iphoto(n) = Photofac
       endif

       !
       ! those with time dependancies
       !

       ! total labile release
       FLUXES(n, 8)  = POOLS(n, 1)*(1d0-(1d0-FLUXES(n, 16))**deltat(n))*deltat_1(n)
       ! total leaf litter production
       FLUXES(n, 10) = POOLS(n, 2)*(1d0-(1d0-FLUXES(n, 9))**deltat(n))*deltat_1(n)
       ! total wood litter production
       ! NOTE: we assume that the coarse root fraction (pars(29)) is passed to the soil directly
       FLUXES(n, 11) = POOLS(n, 4)*(1d0-(1d0-pars(6))**deltat(n))*deltat_1(n)
       ! total root litter production
       FLUXES(n, 12) = POOLS(n, 3)*(1d0-(1d0-pars(7))**deltat(n))*deltat_1(n)

       !
       ! those with temperature AND time dependancies
       !

       ! turnover of litter
       tmp = POOLS(n, 5)*(1d0-(1d0-FLUXES(n, 2)*pars(8))**deltat(n))*deltat_1(n)
       ! respiration heterotrophic litter; decomposition of litter to som
       FLUXES(n, 13) = tmp * (1d0-pars(1)); FLUXES(n, 15) = tmp*pars(1)
       ! respiration heterotrophic som
       FLUXES(n, 14) = POOLS(n, 6)*(1d0-(1d0-FLUXES(n, 2)*pars(9))**deltat(n))*deltat_1(n)

       ! respiration heterotrophic wood litter; decomposition of wood litter to som
       tmp = POOLS(n, 7)*(1d0-(1d0-FLUXES(n, 2)*pars(35))**deltat(n))*deltat_1(n)
       FLUXES(n, 4) = tmp * (1d0-pars(1)); FLUXES(n, 20) = tmp*pars(1)

       !
       ! Update Rg, GPP and NEE fluxes
       !

       ! calculate growth respiration and adjust allocation to pools assuming
       ! 0.21875 of total C allocation towards each pool (i.e. 0.28 .eq. xNPP)
       ! foliage. NOTE: that in the current version only Rg_fol comes from Clabile
       Rg_from_labile(n) = FLUXES(n, 8)*Rg_fraction; FLUXES(n, 8) = FLUXES(n, 8) * one_Rg_fraction
       ! now update the Ra flux
       FLUXES(n, 3) = FLUXES(n, 3) + Rg_from_labile(n)
       ! roots
       FLUXES(n, 3) = FLUXES(n, 3) + (FLUXES(n, 6)*Rg_fraction); FLUXES(n, 6) = FLUXES(n, 6) * one_Rg_fraction
       ! wood
       FLUXES(n, 3) = FLUXES(n, 3) + (FLUXES(n, 7)*Rg_fraction); FLUXES(n, 7) = FLUXES(n, 7) * one_Rg_fraction

       ! calculate the NEE
       NEE_out(n) = (-FLUXES(n, 1)+FLUXES(n, 3)+FLUXES(n, 13)+FLUXES(n, 14)+FLUXES(n, 4))
       ! load GPP
       GPP_out(n) = FLUXES(n, 1)

       !
       ! update pools for next timestep
       !

       ! labile pool
       POOLS(n+1, 1) = POOLS(n, 1) + (FLUXES(n, 5)-FLUXES(n, 8)-Rg_from_labile(n))*deltat(n)
       ! foliar pool
       POOLS(n+1, 2) = POOLS(n, 2) + (FLUXES(n, 8)-FLUXES(n, 10))*deltat(n)
       ! wood pool
       POOLS(n+1, 4) = POOLS(n, 4) + (FLUXES(n, 7)-FLUXES(n, 11))*deltat(n)
       ! root pool
       POOLS(n+1, 3) = POOLS(n, 3) + (FLUXES(n, 6)-FLUXES(n, 12))*deltat(n)
       ! litter pool
       POOLS(n+1, 5) = POOLS(n, 5) + (FLUXES(n, 10)+FLUXES(n, 12)-FLUXES(n, 13)-FLUXES(n, 15))*deltat(n)
       ! som pool
       POOLS(n+1, 6) = POOLS(n, 6) + (FLUXES(n, 15)+FLUXES(n, 20)-FLUXES(n, 14))*deltat(n)
!       POOLS(n+1, 6) = POOLS(n, 6) + ((FLUXES(n, 11)*pars(29))+FLUXES(n, 15)+FLUXES(n, 20)-FLUXES(n, 14))*deltat(n)
       ! wood litter pool
       POOLS(n+1, 7) = POOLS(n, 7) + (FLUXES(n, 11)-FLUXES(n, 20)-FLUXES(n, 4))*deltat(n)
!       POOLS(n+1, 7) = POOLS(n, 7) + ((FLUXES(n, 11)*(1d0-pars(29)))-FLUXES(n, 20)-FLUXES(n, 4))*deltat(n)

       !!!!!!!!!!
       ! deal first with deforestation
       !!!!!!!!!!

       if (n == reforest_day) then
           POOLS(n+1, 1) = pars(30)
           POOLS(n+1, 2) = pars(31)
           POOLS(n+1, 3) = pars(32)
           POOLS(n+1, 4) = pars(33)
       endif

       ! reset values
       harvest_management = 0; burnt_area = 0d0

       ! Does harvest activities occur?
       if (met(8, n) > 0d0) then

           ! Load the management type/scenario into local variable
           harvest_management = int(met(13, n))

           ! Determine the fraction of cut labile C which remains in system as residue.
           ! We assume that labile is proportionally distributed through the plants
           ! root and wood (structural C).
           C_total = POOLS(n+1, 3) + POOLS(n+1, 4)
           ! Ensure there is available C for extraction
           if (C_total > 0d0) then
               ! Harvest activities on the wood/structural pool varies depending on
               ! whether it is above or below ground. As such, partition the wood pool
               ! between above ground stem(+branches) and below ground coarse root.
               Crootcr = POOLS(n+1, 4)*Crootcr_part(harvest_management)
               Cstem   = POOLS(n+1, 4)-Crootcr
               ! Calculate the fraction of harvested labile which remains in system as residue
               labile_frac_res = ((POOLS(n+1, 3)/C_total) * roots_frac_res(harvest_management)  ) &
                               + ((Cstem/C_total)        * stem_frac_res(harvest_management)   ) &
                               + ((Crootcr/C_total)      * rootcr_frac_res(harvest_management) )
               ! Calculate the management scenario specific resistance fraction
               labile_frac_removal = ((POOLS(n+1, 3)/C_total) * roots_frac_removal(harvest_management)  ) &
                                     + ((Cstem/C_total)        * 1d0   ) &
                                     + ((Crootcr/C_total)      * rootcr_frac_removal(harvest_management) )

               ! Calculate the total loss from biomass pools
               ! We assume that fractional clearing always equals the fraction
               ! of foliage and above ground (stem) wood removal. However, we assume
               ! that coarse root and fine root extractions are dependent on the
               ! management activity type, e.g. in coppice below ground remains.
               ! Thus, labile extractions are also dependent.
               labile_loss = POOLS(n+1, 1) * labile_frac_removal*met(8, n)
               foliar_loss = POOLS(n+1, 2) * met(8, n)
               roots_loss  = POOLS(n+1, 3) * roots_frac_removal(harvest_management) * met(8, n)
               stem_loss   = (Cstem*met(8, n))
               rootcr_loss = (Crootcr*rootcr_frac_removal(harvest_management) * met(8, n))
               wood_loss   =  stem_loss+rootcr_loss

               ! Transfer fraction of harvest waste to litter, wood litter or som pools.
               ! This includes explicit calculation of the stem and coarse root residues due
               ! to their potentially different treatments under management scenarios
               labile_residue = labile_loss*labile_frac_res
               foliar_residue = foliar_loss*foliage_frac_res(harvest_management)
               roots_residue  = roots_loss*roots_frac_res(harvest_management)
               coarse_root_residue = rootcr_loss*rootcr_frac_res(harvest_management)
               stem_residue = stem_loss*stem_frac_res(harvest_management)
               wood_residue = stem_residue+coarse_root_residue
               ! Mechanical loss of Csom due to coarse root extraction, 
               ! less the loss remaining as residue
               soil_loss_with_roots = (rootcr_loss-coarse_root_residue) &
                                    * soil_loss_frac(harvest_management)

               ! Update pools
               POOLS(n+1, 1) = POOLS(n+1, 1) - labile_loss
               POOLS(n+1, 2) = POOLS(n+1, 2) - foliar_loss
               POOLS(n+1, 3) = POOLS(n+1, 3) - roots_loss
               POOLS(n+1, 4) = POOLS(n+1, 4) - wood_loss
               POOLS(n+1, 5) = POOLS(n+1, 5) + (labile_residue+foliar_residue+roots_residue)
               POOLS(n+1, 6) = POOLS(n+1, 6) - soil_loss_with_roots
               POOLS(n+1, 7) = POOLS(n+1, 7) + wood_residue
               ! mass balance check
               where (POOLS(n+1, 1:7) < 0d0) POOLS(n+1, 1:7) = 0d0

               ! Convert total losses into extracted C by accounting for the residues.
               ! Conver these totals into to daily rate for output
               harvest_extracted_labile(n)     = (labile_loss-labile_residue) * deltat_1(n)
               harvest_extracted_foliar(n)     = (foliar_loss-foliar_residue) * deltat_1(n)
               harvest_extracted_roots(n)      = (roots_loss-roots_residue) * deltat_1(n)
               harvest_extracted_wood(n)       = (wood_loss-wood_residue) * deltat_1(n)
               harvest_extracted_litter(n)     = 0d0
               harvest_extracted_woodlitter(n) = 0d0
               harvest_extracted_som(n)        = soil_loss_with_roots*deltat_1(n)
               ! Convert harvest related residue generations to daily rate for output
               harvest_residue_labile(n) = labile_residue*deltat_1(n)
               harvest_residue_foliar(n) = foliar_residue*deltat_1(n)
               harvest_residue_roots(n)  = roots_residue*deltat_1(n)
               harvest_residue_wood(n)   = wood_residue*deltat_1(n)
               ! Convert the combined residue to dead organic matter pools
               harvest_residue_to_litter(n)     = (labile_residue+foliar_residue+roots_residue) * deltat_1(n)
               harvest_residue_to_woodlitter(n) = wood_residue*deltat_1(n)
               harvest_residue_to_som(n)        = 0d0

               ! Total C extraction, this includes any som or litter clearing
               FLUXES(n, 21) = harvest_extracted_labile(n) &
                            + harvest_extracted_foliar(n) &
                            + harvest_extracted_roots(n) &
                            + harvest_extracted_wood(n) &
                            + harvest_extracted_litter(n) &
                            + harvest_extracted_woodlitter(n) &
                            + harvest_extracted_som(n)

           end if  ! C_total > 0d0

           ! If total clearance occured then we need to ensure some minimum
           ! values and reforestation is assumed one year forward
           if (met(8, n) > 0.99d0) then
               m = 0; test = nint(sum(deltat(n:(n+m))))
               ! FC Forest Statistics 2015 lag between harvest and restocking ~ 2 year
               restocking_lag = 365*2
               do while (test < restocking_lag)
                   m = m+1; test = nint(sum(deltat(n:(n+m))))
                   !  get out clause for hitting the end of the simulation
                   if (m+n >= nodays) test = restocking_lag
               enddo
               reforest_day = min((n+m), nodays)
           endif  ! if total clearance

       endif  ! end deforestation info

       !!!!!!!!!!
       ! then deal with fire
       !!!!!!!!!!

       if (met(9, n) > 0d0 .or.(met(8, n) > 0d0 .and. harvest_management > 0)) then

           burnt_area = met(9, n)
           if (met(8, n) > 0d0 .and. burnt_area > 0d0) then
               ! pass harvest management to local integer
               burnt_area = min(1d0, burnt_area+post_harvest_burn(harvest_management))
           else if (met(8, n) > 0d0 .and. burnt_area <= 0d0) then
               burnt_area = post_harvest_burn(harvest_management)
           endif

           if (burnt_area > 0d0) then

               ! first calculate combustion/emissions fluxes in g C m-2 d-1
               fire_emiss_labile(n) = POOLS(n+1, 1)*burnt_area*combust_eff(1)*deltat_1(n)  ! labile
               fire_emiss_foliar(n) = POOLS(n+1, 2)*burnt_area*combust_eff(2)*deltat_1(n)  ! foliar
               fire_emiss_roots(n) = POOLS(n+1, 3)*burnt_area*combust_eff(3)*deltat_1(n)  ! roots
               fire_emiss_wood(n) = POOLS(n+1, 4)*burnt_area*combust_eff(4)*deltat_1(n)  ! wood
               fire_emiss_litter(n) = POOLS(n+1, 5)*burnt_area*combust_eff(5)*deltat_1(n)  ! litter
               fire_emiss_som(n) = POOLS(n+1, 6)*burnt_area*combust_eff(6)*deltat_1(n)  ! som
               fire_emiss_woodlitter(n) = POOLS(n+1, 7)*burnt_area*combust_eff(7)*deltat_1(n)  ! litterwood

               ! second calculate litter transfer fluxes in g C m-2 d-1, all pools except som
               fire_litter_labile(n) = POOLS(n+1, 1)*burnt_area*(1d0-combust_eff(1))*(1d0-rfac(1))*deltat_1(n)  ! labile into litter
               fire_litter_foliar(n) = POOLS(n+1, 2)*burnt_area*(1d0-combust_eff(2))*(1d0-rfac(2))*deltat_1(n)  ! foliar into litter
               fire_litter_roots(n) = POOLS(n+1, 3)*burnt_area*(1d0-combust_eff(3))*(1d0-rfac(3))*deltat_1(n)  ! roots into litter
               fire_litter_wood(n) = POOLS(n+1, 4)*burnt_area*(1d0-combust_eff(4))*(1d0-rfac(4))*deltat_1(n)  ! wood into woodlitter
               fire_litter_litter(n) = POOLS(n+1, 5)*burnt_area*(1d0-combust_eff(5))*(1d0-rfac(5))*deltat_1(n)  ! litter into som
               !fire_litter_som(n) = ! no litter generation from the som pool
               fire_litter_woodlitter(n) = POOLS(n+1, 7)*burnt_area*(1d0-combust_eff(7))*(1d0-rfac(7))*deltat_1(n)  ! wood litter into som

               ! Fire flux (gC/m2/day)
               FLUXES(n, 17) = fire_emiss_labile(n) + fire_emiss_foliar(n) + fire_emiss_roots(n) &
                            + fire_emiss_wood(n) + fire_emiss_litter(n) + fire_emiss_woodlitter(n) + fire_emiss_som(n)

               ! Residue redistribution
               fire_residue_to_litter(n) = fire_litter_labile(n) + fire_litter_foliar(n) + fire_litter_roots(n)
               fire_residue_to_som(n)    = fire_litter_litter(n) + fire_litter_woodlitter(n)
               fire_residue_to_woodlitter(n)=  fire_litter_wood(n)

               ! Update pools
               POOLS(n+1, 1) = POOLS(n+1, 1)-(fire_emiss_labile(n)+fire_litter_labile(n)) * deltat(n)
               POOLS(n+1, 2) = POOLS(n+1, 2)-(fire_emiss_foliar(n)+fire_litter_foliar(n)) * deltat(n)
               POOLS(n+1, 3) = POOLS(n+1, 3)-(fire_emiss_roots(n)+fire_litter_roots(n)) * deltat(n)
               POOLS(n+1, 4) = POOLS(n+1, 4)-(fire_emiss_wood(n)+fire_litter_wood(n)) * deltat(n)
               POOLS(n+1, 5) = POOLS(n+1, 5)+(fire_residue_to_litter(n)-fire_emiss_litter(n)-fire_litter_litter(n)) * deltat(n)
               POOLS(n+1, 6) = POOLS(n+1, 6)+(fire_residue_to_som(n)-fire_emiss_som(n)) * deltat(n)
               POOLS(n+1, 7) = POOLS(n+1, 7)+(fire_residue_to_woodlitter(n)- &
                                            fire_emiss_woodlitter(n)- &
                                            fire_litter_woodlitter(n)) * deltat(n)
               ! mass balance check
               where (POOLS(n+1, 1:7) < 0d0) POOLS(n+1, 1:7) = 0d0

           endif  ! burn area > 0

       endif  ! fire activity

       do nxp = 1, nopools
          if (POOLS(n+1, nxp) /= POOLS(n+1, nxp) .or. POOLS(n+1, nxp) < 0d0) then
              print*,"step",n, "POOL",nxp
              print*,"met",met(:,n)
              print*,"POOLS",POOLS(n, :)
              print*,"FLUXES",FLUXES(n, :)
              print*,"POOLS+1",POOLS(n+1, :)
              print*,"PARS",pars
              stop
          endif
       enddo

    end do  ! nodays loop
    !call cpu_time(done)
    !print*,done-begin
  end subroutine CARBON_MODEL
  !
  !------------------------------------------------------------------
  !
  subroutine acm_gpp_stage_1

    ! Estimate the light and temperature limited photosynthesis components.
    ! See acm_gpp_stage_2() for estimation of CO2 supply limitation and
    ! combination of light, temperature and CO2 co-limitation

    implicit none

    ! Declare local variables
    double precision:: a, b, c, Pl_max, PAR_m2, airt_ad

    !
    ! Metabolic limited photosynthesis
    !

    ! Determine the appropriate canopy top
    ! Estimate the total canopy N, then scale to the "top leaf" and multiple by NUE
    ceff = (((avN*lai) / leaf_canopy_light_scaling) * NUE)

    ! maximum rate of temperature and nitrogen (canopy efficiency) limited
    ! photosynthesis (gC.m-2.day-1 -> umolC/m2/s). Scaling from leaf to canopy
    ! scaled assumed to follow integral of light environment.
    metabolic_limited_photosynthesis = gC_to_umol*leaf_canopy_light_scaling*ceff*seconds_per_day_1 &
                                     * opt_max_scaling(pn_max_temp, pn_min_temp, pn_opt_temp, pn_kurtosis, leafT)

    !
    ! Light limited photosynthesis
    !

    ! Calculate light limted rate of photosynthesis (gC.m-2.day-1) as a function
    ! light capture and leaf to canopy scaling on quantum yield (e0).
    light_limited_photosynthesis = e0*canopy_par_MJday*dayl_seconds_1*gC_to_umol

    !
    ! Stomatal conductance independent variables for diffusion limited
    ! photosynthesis
    !

    ! Canopy level boundary layer conductance unit change
    ! (m.s-1 -> mol.m-2.d-1) assuming sea surface pressure only.
    ! Note the ratio of H20:CO2 diffusion through leaf level boundary layer is
    ! 1.37 (Jones appendix 2). Note conversion to resistance for easiler merging
    ! with stomatal conductance in acm_gpp_stage_2).
    rb_mol_1 = (aerodynamic_conductance*convert_ms1_mol_1*gb_H2O_CO2 * &
              leaf_canopy_wind_scaling) ** (-1d0)

    ! Arrhenious Temperature adjustments for Michaelis-Menten coefficients
    ! for CO2 (kc) and O2 (ko) and CO2 compensation point
    ! See McMurtrie et al., (1992) Australian Journal of Botany, vol 40, 657-677
    co2_half_sat   = arrhenious(kc_half_sat_25C, kc_half_sat_gradient, leafT)
    co2_comp_point = arrhenious(co2comp_sat_25C, co2comp_gradient, leafT)

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
    double precision:: pp, qq, mult, rc, pd

    !
    ! Combined diffusion limitation and carboxylation limited photosynthesis
    !

    ! Estimation of ci is based on the assumption that metabilic limited
    ! photosynthesis is equal to diffusion limited. For details
    ! see Williams et al, (1997), Ecological Applications, 7(3), 1997, pp. 882–894

    ! Daily canopy conductance dertmined through combination of aerodynamic and
    ! stomatal conductances. Both conductances are scaled to canopy aggregate.
    ! aerodynamic conductance already in units of molCO2.m-2.s-1 (see acm_gpp_stage_1).
    ! Stomatal conductance scaled from mmolH2O to molCO2.
    ! The ratio of H20:CO2 diffusion is 1.646259 (Jones appendix 2).
    !
    ! Combining in series the stomatal and boundary layer conductances
    ! to make canopy resistence (s/m2/molCO2)
    rc = (gs*gs_H2Ommol_CO2mol) ** (-1d0) + rb_mol_1

    ! pp and qq represent limitation by metabolic (temperature & N) and
    ! diffusion (co2 supply) respectively
    pp = metabolic_limited_photosynthesis*rc; qq = co2_comp_point-co2_half_sat
    mult = co2+qq-pp
    ! calculate internal CO2 concentration (ppm or umol/mol)
    ci = 0.5d0*(mult+sqrt((mult*mult)-4d0*(co2*qq-pp*co2_comp_point)))

    ! calculate CO2 limited rate of photosynthesis (umolC.m-2.s-1)
    pd = ((co2-ci)/rc)

    !
    ! Estimate CO2 and light co-limitation
    !

    ! calculate combined light and CO2 limited photosynthesis (umolC/m2/s)
    acm_gpp_stage_2 = light_limited_photosynthesis*pd/(light_limited_photosynthesis+pd)

    !pp = acm_gpp_stage_2*rc; mult = co2+qq-pp
    !! calculate internal CO2 concentration (ppm or umol/mol)
    !ci = 0.5d0*(mult+sqrt((mult*mult)-4d0*(co2*qq-pp*co2_comp_point)))

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

    !!!!!!!!!!
    ! Optimise intrinsic water use efficiency
    !!!!!!!!!!

    ! Determine impact of gs increment on pd and how far we are from iWUE
    find_gs_iWUE = iWUE_step - (acm_gpp_stage_2(gs_in+delta_gs) - acm_gpp_stage_2(gs_in))

    ! Remember to return back to the user
    return

  end function find_gs_iWUE
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_stomatal_conductance

    ! Determines an approximation of canopy scale stomatal conductance (gc)
    ! mmolH2O.m-2.s-1 based on potential hydraulic flow, air temperature and absorbed radiation.

    use brent_zero, only: zbrent

    implicit none

    ! local variables
    double precision:: denom, iWUE_upper!, iWUE_lower
    double precision, parameter:: max_gs = 1000d0, &  ! mmolH2O.m-2.s-1 (leaf area)
                                   min_gs = 0.01d0, &  ! mmolH2O.m-2.s-1 (leaf area)
                                   tol_gs = 0.01d0     ! mmolH2O.m-2.s-1 (leaf area)

    !!!!!!!!!!
    ! Calculate stomatal conductance under H2O and CO2 limitations
    !!!!!!!!!!

    if (aerodynamic_conductance > vsmall .and. total_water_flux > vsmall) then

        ! Determine potential water flow rate (mmolH2O.m-2.s-1)
        max_supply = total_water_flux

        ! Pass minimum conductance from local parameter to global value
        minimum_conductance = min_gs*leaf_canopy_light_scaling

        ! Invert Penman-Monteith equation to give gs (m.s-1) needed to meet
        ! maximum possible evaporation for the day.
        ! This will then be reduced based on CO2 limits for diffusion based
        ! photosynthesis
        denom = slope * (((canopy_swrad_MJday*1d6*dayl_seconds_1) + canopy_lwrad_Wm2)) &
              + (ET_demand_coef*aerodynamic_conductance*leaf_canopy_wind_scaling)
        denom = (denom / (lambda*max_supply*mmol_to_kg_water)) - slope
        potential_conductance = (aerodynamic_conductance*leaf_canopy_wind_scaling) / (denom/psych)

        ! convert m.s-1 to mmolH2O.m-2.d-1, per unit ground area, note that this
        ! is implicitly the canopy scaled value
        potential_conductance = potential_conductance*convert_ms1_mmol_1
        ! if conditions are dew forming then set conductance to maximum as we
        ! are not going to be limited by water demand
        if (potential_conductance <= 0d0 .or. potential_conductance > max_gs*leaf_canopy_light_scaling) then
            potential_conductance = max_gs*leaf_canopy_light_scaling
        end if

        ! If there is a positive demand for water then we will solve for
        ! photosynthesis limits on gs through iterative solution

        ! Determine the appropriate canopy scaled gs increment and return threshold
        delta_gs = 1d0*leaf_canopy_light_scaling  ! mmolH2O/m2leaf/s
        iWUE_step = iWUE*leaf_canopy_light_scaling  ! umolC/mmolH2Ogs/s

        ! Calculate stage one acm, temperature and light limitation which
        ! are independent of stomatal conductance effects
        call acm_gpp_stage_1

        ! Intrinsic WUE optimisation
        ! Check that the water restricted water range brackets the root solution for the bisection
        iWUE_upper = find_gs_iWUE(potential_conductance) !; iWUE_lower = find_gs_iWUE(min_gs)
        if ( iWUE_upper*find_gs_iWUE(min_gs) > 0d0 ) then
             ! Then both proposals indicate that photosynthesis
             ! would be increased by greater opening of the stomata
             ! and is therefore water is limiting!
             stomatal_conductance = potential_conductance
             ! Exception being if both are positive-therefore assume
             ! lowest
             if (iWUE_upper > 0d0) stomatal_conductance = minimum_conductance
        else
             ! In all other cases iterate
             stomatal_conductance = zbrent('calculate_gs:find_gs_iWUE', &
                                           find_gs_iWUE, minimum_conductance, potential_conductance, tol_gs*lai, iWUE_step*0.10d0)
        end if

    else

        ! if no LAI then there can be no stomatal conductance
        potential_conductance = max_gs; minimum_conductance = vsmall
        stomatal_conductance = vsmall
        ! set minimum (computer) precision level flow
        max_supply = vsmall

    endif  ! if aerodynamic conductance > vsmall

  end subroutine calculate_stomatal_conductance
  !
  !------------------------------------------------------------------
  !
  subroutine meteorological_constants(input_temperature, input_temperature_K, input_vpd_kPa)

    ! Determine some multiple use constants used by a wide range of functions
    ! All variables here are linked to air temperature and thus invarient between
    ! iterations and can be stored in memory...

    implicit none

    ! arguments
    double precision, intent(in):: input_temperature, input_temperature_K, &
                                    input_vpd_kPa

    ! local variables
    double precision:: mult, &
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
    !canopy_wind = canopy_wind*exp((ustar_Uh*((canopy_height*0.5d0)-canopy_height))/mixing_length_momentum)
    ! Estimate canopy scaling factor for use with aerodynamic conductance.
    leaf_canopy_wind_scaling = exp((ustar_Uh/mixing_length_momentum)) &
                             / (ustar_Uh/mixing_length_momentum)

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
    gv_forced = water_vapour_diffusion*Sh_forced*leaf_width_coef

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
    ! 1/vonkarman*0.31 = 0.7710906
    canopy_wind = ustar*vonkarman_1*log((canopy_height-displacement) / roughl)

    ! set minimum value for wind speed at canopy top (m.s-1)
!    canopy_wind = max(min_wind, canopy_wind)

  end subroutine log_law_decay
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
                                  ,clump = 0.75d0   & ! Clumping factor (1 = uniform, 0 totally clumped, mean = 0.75)
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
    ! canopy to soil/sky. The first line below estimates the actual release fraction, 
    ! second line accounts for changing LAI scaling within a vertical canopy.
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
                       ,par, nir                     &
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
    double precision, parameter:: clump = 0.75d0 & ! Clumping factor (1 = uniform, 0 totally clumped, mean = 0.75)
                                                    ! He et al., (2012) http://dx.doi.org/10.1016/j.rse.2011.12.008
                                  ,decay = -0.5d0   ! decay coefficient for incident radiation

    !!!!!!!!!!
    ! Determine canopy absorption, reflectance and transmittance as function of
    ! LAI
    !!!!!!!!!!

    ! First, we consider how much radiation is likely to be incident on the
    ! canopy, or put another way what fraction passes straight through the
    ! canopy?
    transmitted_fraction = exp(decay*lai*clump)
    ! Estimate the integral of light interception for use as a leaf to canopy
    ! scaler for photoynthesis, transpiration, and gs
    leaf_canopy_light_scaling = (1d0-transmitted_fraction) / (-decay*clump)

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

    ! Then the radiation incident and ultimately absorbed by the soil surface
    ! itself (MJ.m-2.day-1)
    soil_par_MJday = trans_par_MJday*soil_swrad_absorption
    soil_nir_MJday = trans_nir_MJday*soil_swrad_absorption
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
    integer:: i, rooted_layer
    double precision:: bonus, &
                        transpiration_resistance, root_reach_local, &
                        root_depth_50, slpa, mult, exp_func, prev
    double precision, dimension(nos_root_layers):: root_mass    &
                                                   ,root_length  &
                                                   ,ratio
    double precision, parameter:: rootdist_tol = 13.81551d0  ! log(1d0/rootdist_tol-1d0) were rootdist_tol = 1d-6
                                  !rootdist_tol = 1d-6!, & ! Root density assessed for the max rooting depth
                                   !root_depth_frac_50 = 0.25d0  ! fractional soil depth above which 50 %
                                                               ! of the root mass is assumed to be located

    ! reset water flux
    total_water_flux = 0d0; water_flux_mmolH2Om2s = 0d0
    root_mass = 0d0
    ! calculate the plant hydraulic resistance component. Currently unclear
    ! whether this actually varies with height or whether tall trees have a
    ! xylem architecture which keeps the whole plant conductance (gplant) 1-10 (ish).
    !    transpiration_resistance = (gplant*lai)**(-1d0)
    transpiration_resistance = canopy_height / (gplant*max(min_lai, lai))

    !!!!!!!!!!!
    ! Calculate root profile
    !!!!!!!!!!!

    ! calculate soil depth to which roots reach
    root_reach = max_depth*root_biomass / (root_k+root_biomass)
    ! Determine initial soil layer thickness
    layer_thickness(1) = top_soil_depth; layer_thickness(2) = max(min_layer, root_reach-top_soil_depth)
    layer_thickness(3) = max_depth-sum(layer_thickness(1:2))
    layer_thickness(4) = top_soil_depth

    ! The original SPA src generates an exponential distribution which aims
    ! to maintain 50 % of root biomass in the top 25 % of the rooting depth.
    ! In a simple 3 root layer system this can be estimates more simply

!    ! top 25 % of root profile
!    root_depth_50 = root_reach*root_depth_frac_50
!    if (root_depth_50 <= layer_thickness(1)) then
!
!        ! Greater than 50 % of the fine root biomass can be found in the top
!        ! soil layer
!
!        ! Start by assigning all 50 % of root biomass to the top soil layer
!        root_mass(1) = fine_root_biomass*0.5d0
!        ! Then quantify how much additional root is found in the top soil layer
!        ! assuming that the top 25 % depth is found somewhere within the top
!        ! layer
!        bonus = (fine_root_biomass-root_mass(1)) &
!              * (layer_thickness(1)-root_depth_50) / (root_reach-root_depth_50)
!        root_mass(1) = root_mass(1) + bonus
!        ! partition the remaining root biomass between the seconds and third
!        ! soil layers
!        if (root_reach > sum(layer_thickness(1:2))) then
!            root_mass(2) = (fine_root_biomass-root_mass(1)) &
!                         * (layer_thickness(2)/(root_reach-layer_thickness(1)))
!            root_mass(3) = fine_root_biomass-sum(root_mass(1:2))
!        else
!            root_mass(2) = fine_root_biomass-root_mass(1)
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
!        root_mass(3) = fine_root_biomass-sum(root_mass(1:2))
!
!    else
!
!        ! Greater than 50 % of fine root biomass stock spans across all three
!        ! layers
!        root_mass(1:2) = fine_root_biomass*0.5d0 * (layer_thickness(1:2)/root_depth_50)
!        root_mass(3) = fine_root_biomass-sum(root_mass(1:2))
!
!    endif
!    ! now convert root mass into lengths
!    root_length = root_mass*root_mass_length_coef_1
!!    root_length = root_mass / (root_density*root_cross_sec_area)

    !!!!!!!!!!!
    ! Calculate hydraulic properties and each rooted layer
    !!!!!!!!!!!

    ! calculate and accumulate steady state water flux in mmol.m-2.s-1
    ! NOTE: Depth correction already accounted for in soil resistance
    ! calculations and this is the maximum potential rate of transpiration
    ! assuming saturated soil and leaves at their minimum water potential.
    demand = -minlwp - (head*canopy_height)
    ! now loop through soil layers, where root is present
    rooted_layer = 1
    ! Determine the exponential coefficient needed for an exponential decay to the current root reach
    ! Exponential decay profile following:
    ! Y = 1 / (1+exp(-B*Z)), where Y = density at Z, B = gradient, Z = depth
    ! To determine gradient for current maximum root depth assuming density reaches rootdist_tol value, rearranges to:
    ! B = ln(1/Y - 1) / Z
!d = seq(0, 2, 0.01); c = -2.6; rmax = 1; d50 = 0.25; rd = rmax / (1+ (d/d50)**c)
!rmax = rd * (1 + (d/d50)**c)
!(((rmax/rd) - 1)**(1/c)) * d50 = d  ! Depth at which rd = 99 %
    !slpa = log(1d0/rootdist_tol-1d0) / root_reach
    slpa = rootdist_tol/root_reach
    prev = 1d0
    do i = 1, nos_root_layers
       ! Determine the exponential function for the current cumulative depth
       exp_func = exp(-slpa*sum(layer_thickness(1:i)))
       ! Calculate the difference in the integral between depths, i.e. the proportion of root in the current volume
       mult = prev - (1d0 - (1d0/(1d0+exp_func)) + (0.5d0*exp_func))
       ! Assign fine roo the the current layer...
       root_mass(i) = root_biomass*mult
       ! and determine the associated amount of root
       root_length(i) = root_mass(i) * root_mass_length_coef_1
       prev = prev-mult
       ! If there is root in the current layer then we should calculate the resistances
       if (root_mass(i) > 0d0) then
           ! Track the deepest root layer assessed
           rooted_layer = i
           ! if there is root then there is a water flux potential...
           root_reach_local = min(root_reach, layer_thickness(i))
           ! calculate and accumulate steady state water flux in mmol.m-2.s-1
           call plant_soil_flow(i, root_length(i), root_mass(i) &
                               ,demand(i), root_reach_local &
                               ,transpiration_resistance)
       else
           ! ...if there is not then we wont have any below...
           exit
       end if  ! root present in current layer?
    end do  ! nos_root_layers

    ! if freezing then assume soil surface is frozen, therefore no water flux
    if (meant < 1d0) water_flux_mmolH2Om2s(1) = 0d0

    ! calculate sum value (mmolH2O.m-2.s-1)
    total_water_flux = sum(water_flux_mmolH2Om2s)
    ! calculate weighted SWP and uptake fraction
    uptake_fraction(1:nos_root_layers) = water_flux_mmolH2Om2s(1:nos_root_layers) / total_water_flux
    ! determine effective resistance (MPa.s-1.m-2.mmol-1)
    Rtot = -minlwp/total_water_flux

    ! and return
    return

  end subroutine calculate_Rtot
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
    ustar_Uh = 0.3d0
!    ustar_Uh = max(ustar_Uh_min, min(sqrt(Cs+Cr*local_lai*0.5d0), ustar_Uh_max))
!    ustar_Uh = sqrt(Cs+Cr*local_lai*0.5d0)

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
  subroutine plant_soil_flow(root_layer, root_length, root_mass &
                            ,demand, root_reach_in, transpiration_resistance)

    !
    ! Calculate soil layer specific water flow form the soil to canopy (mmolH2O.m-2.s-1)
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
    double precision:: soilR2, Rtot_layer

    ! Calculates root hydraulic resistance (MPa m2 s mmol-1) in a soil-root zone
    soilR2 = root_resist / (root_mass*root_reach_in)
    ! Estimate the total hydraulic resistance for the layer
    Rtot_layer = transpiration_resistance+soilR2

    ! Estimate the soil to plant flow of water mmolH2O/m2/s
    water_flux_mmolH2Om2s(root_layer) = demand/Rtot_layer

    ! return
    return

  end subroutine plant_soil_flow
  !
  !------------------------------------------------------------------
  !
  pure function arrhenious( a, b, t )

    ! The equation is simply...                        !
    !    a*exp( b * ( t-25.0 ) / ( t+273.15 ) )  !
    ! However, precision in this routine matters as it !
    ! affects many others. To maximise precision, the  !
    ! calculations have been split.                    !

    implicit none

    ! arguments..
    double precision, intent(in):: a, b, t
    double precision            :: arrhenious

    ! local variables..
    double precision:: answer, denominator, numerator

    numerator   = t-25d0
    denominator = t+freeze
    answer      = a*exp( b*numerator/denominator )
    arrhenious  = answer

  end function arrhenious
  !
  !------------------------------------------------------------------
  !
  double precision function linear_model_gradient(x, y, interval)

    ! Function to calculate the gradient of a linear model for a given depentent
    ! variable (y) based on predictive variable (x). The typical use of this
    ! function will in fact be to assume that x is time.

    implicit none

    ! declare input variables
    integer:: interval
    double precision, dimension(interval):: x, y

    ! declare local variables
    double precision:: sum_x, sum_y, sumsq_x, sum_product_xy

    ! calculate the sum of x
    sum_x = sum(x)
    ! calculate the sum of y
    sum_y = sum(y)
    ! calculate the sum of squares of x
    !sumsq_x = sum(x*x)
    ! calculate the sum of the product of xy
    !sum_product_xy = sum(x*y)
    ! calculate the gradient
    !linear_model_gradient = ( (dble(interval)*sum_product_xy) - (sum_x*sum_y) )
    !&
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
  !----------------------------------------------------------------------
  !
  double precision function opt_max_scaling( max_val, min_val, optimum, kurtosis, current )

    ! Estimates a 0-1 scaling based on a skewed guassian distribution with a
    ! given optimum, maximum and kurtosis. Minimum is assumed to be at infinity
    ! (or near enough)

    implicit none

    ! arguments..
    double precision, intent(in):: max_val, min_val, optimum, kurtosis, current

!    ! local variables..
!    double precision, parameter:: min_val = -1d6

    ! Code with implicit assumption of min bound at infinity
!    if ( current >= max_val ) then
!         opt_max_scaling = 0d0
!    else
!         dummy     = exp( log((max_val-current) / (max_val-optimum)) * kurtosis * (max_val-optimum) )
!         opt_max_scaling = dummy*exp( kurtosis * ( current-optimum ) )
!    end if

    ! Code with explicit min bound
    opt_max_scaling = exp( kurtosis*log((max_val-current)/(max_val-optimum)) * (max_val-optimum) ) &
                    * exp( kurtosis*log((current-min_val)/(optimum-min_val)) * (optimum-min_val) )
    ! Sanity check, allows for overlapping parameter ranges
    if (opt_max_scaling /= opt_max_scaling) opt_max_scaling = 0d0

  end function opt_max_scaling
!
!--------------------------------------------------------------------
!
end module CARBON_MODEl_MOD
