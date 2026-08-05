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
! This file contains the parameters and memory of DALEC.A1.C2.D2.F2.H1.P4.R2.010
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
! See function / subroutine specific comments for exceptions and contributors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module carbon_model_memory

  implicit none

  ! make all private
  !private
  public

  !!!!!!!!!
  ! Parameters
  !!!!!!!!!

  ! useful technical parameters
  double precision, parameter :: vsmall = tiny(0d0)*1d3 & ! *1d3 to add a little breathing room
                                ,vlarge = huge(0d0)

  integer, parameter :: nos_root_layers = 2, nos_soil_layers = nos_root_layers + 1
  double precision, parameter :: pi = 3.1415927d0,  &
                               pi_1 = 0.3183099d0,  & ! pi**(-1d0)
                             two_pi = 6.283185d0,   & ! pi*2d0
                         deg_to_rad = 0.01745329d0, & ! pi/180d0
                sin_dayl_deg_to_rad = 0.3979486d0,  & ! sin( 23.45d0 * deg_to_rad )
                              boltz = 5.670400d-8,  & ! Boltzmann constant (W.m-2.K-4)
                         emissivity = 0.96d0,       &
                        emiss_boltz = 5.443584d-08, & ! emissivity * boltz
                        ppfd_to_par = 4.6d0,        & ! Conversion of umolPAR -> J
                    sw_par_fraction = 0.5d0,        & ! fraction of short-wave radiation which is PAR
                             freeze = 273.15d0,     & !
                  gs_H2Ommol_CO2mol = 0.001646259d0,& ! The ratio of H20:CO2 diffusion for gs (Jones appendix 2)
              gs_H2Ommol_CO2mol_day = 142.2368d0,   & ! The ratio of H20:CO2 diffusion for gs, including seconds per day correction
                         gb_H2O_CO2 = 1.37d0,       & ! The ratio of H20:CO2 diffusion for gb (Jones appendix 2)
            partial_molar_vol_water = 18.05d-6,     & ! partial molar volume of water, m3 mol-1 at 20C
                     mol_to_g_water = 18d0,         & ! molecular mass of water
                   mmol_to_kg_water = 1.8d-5,       & ! milli mole conversion to kg
                         umol_to_gC = 1.2d-5,       & ! conversion of umolC -> gC
                         gC_to_umol = 83333.33d0,   & ! conversion of gC -> umolC; umol_to_gC**(-1d0)
                               Rcon = 8.3144d0,     & ! Universal gas constant (J.K-1.mol-1)
                          vonkarman = 0.41d0,       & ! von Karman's constant
                        vonkarman_1 = 2.439024d0,   & ! 1 / von Karman's constant
                              cpair = 1004.6d0        ! Specific heat capacity of air; used in energy balance J.kg-1.K-1

  ! photosynthesis / respiration parameters
  double precision, parameter :: &
                        Rg_fraction = 0.21875d0,    & ! fraction of C allocation towards each pool
                                                      ! lost as growth respiration
                                                      ! (i.e. 0.28 .eq. xNPP)
                    one_Rg_fraction = 1d0 - Rg_fraction

  ! hydraulic parameters
  double precision, parameter :: &
                         tortuosity = 2.5d0,        & ! tortuosity
                             gplant = 4d0,          & ! plant hydraulic conductivity (mmol m-1 s-1 MPa-1)
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
                      canopy_height = 9d0,              & ! canopy height assumed to be 9 m
                       tower_height = canopy_height + 2d0, & ! tower (observation) height assumed to be 2 m above canopy
                           min_wind = 0.2d0,            & ! minimum wind speed at canopy top
                          min_layer = 0.03d0,           & ! minimum thickness of the third rooting layer (m)
                        soil_roughl = 0.00085d0,        & ! soil roughness length (m), Meier et al., (2022), https://doi.org/10.5194/gmd-15-2365-2022
                       min_drythick = soil_roughl*10d0, & ! minimum dry thickness depth (m) 0.01 WRF-SPA
                     top_soil_depth = 0.30d0,           & ! thickness of the top soil layer (m)
                           min_root = 5d0,              & ! minimum root biomass (gBiomass.m-2)
                            min_lai = 0.01d0,           & ! minimum LAI assumed for aerodynamic conductance calculations (m2/m2)
                        min_storage = 0.1d0               ! minimum canopy water (surface) storage (mm)

  ! timing parameters
  double precision, parameter :: &
                   seconds_per_hour = 3600d0,       & ! Number of seconds per hour
                    seconds_per_day = 86400d0,      & ! Number of seconds per day
                  seconds_per_day_1 = 1.157407d-05    ! Inverse of seconds per day

  ! ACM-GPP-ET parameters
  double precision, parameter :: &
                       Vc_minT = -6.991d0,     & ! Temperature at which all photosynthetic activity is shutdown
                       Vc_coef = 0.1408d0,     & ! Temperature above Vc_minT that 50% limitation of cold shutdown occurs
                   pn_max_temp = 57.05d0,      & ! Maximum daily max temperature for photosynthesis (oC)
                   pn_min_temp = -1d6,         & ! Minimum daily max temperature for photosynthesis (oC)
                   pn_opt_temp = 30d0,         & ! Optimum daily max temperature for photosynthesis (oC)
                   pn_kurtosis = 0.172d0,      & ! Kurtosis of Jmax temperature response
               ko_half_sat_25C = 157.46892d0,  & ! photorespiration O2 half sat(mmolO2/mol), achieved at 25oC
          ko_half_sat_gradient = 14.93643d0,   & ! photorespiration O2 half sat gradient
               kc_half_sat_25C = 319.58548d0,  & ! carboxylation CO2 half sat (umolCO2/mol), achieved at 25oC
          kc_half_sat_gradient = 24.72297d0,   & ! carboxylation CO2 half sat gradient
               co2comp_sat_25C = 36.839214d0,  & ! carboxylation CO2 compensation point(umolCO2/mol), saturation
              co2comp_gradient = 9.734371d0,   & ! carboxylation CO2 comp point, achieved at oC
                                                 ! Each of these are temperature sensitivty
                           e0 = 3.2d+00,       & ! Quantum yield (gC/MJ/m2/day PAR), SPA apparent yield
               minlwp_default =-1.808224d+00,  & ! minimum leaf water potential (MPa). NOTE: actual SPA = -2 MPa
     soil_iso_to_net_coef_LAI =-2.717467d+00,  & ! Coefficient relating soil isothermal net radiation to net.
                         iWUE = 4.6875d-04,    & ! Intrinsic water use efficiency (umolC/mmolH2O-1/m2leaf/s-1)
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

  !!!!!!!!!
  ! Module variables
  !!!!!!!!!

  type model_working_variables
    double precision :: minlwp = minlwp_default

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
  ! Growing season index canopy phenology model
  integer :: gsi_lag_steps ! Number of model time steps over which GSI is lagged
  double precision, allocatable, dimension(:) :: gsi_lag_days, & ! Number of days equivelent over which GSI is lagged
                                              gsi_lag_history    ! Local storage of the GSI values to be worked on.

  ! hydraulic model variables
  integer :: water_retention_pass, soil_layer
  double precision, dimension(nos_soil_layers) :: &
                   soil_frac_clay,soil_frac_sand, & ! clay and soil fractions of soil
                                     infiltrated    ! surface water infiltrated (kgH2O.m-2.d-1)
  double precision, dimension(nos_root_layers) :: uptake_fraction, & ! fraction of water uptake from each root layer
                                                           demand, & ! maximum potential canopy hydraulic demand
                                            water_flux_mmolH2Om2s, & ! potential transpiration flux (mmolH2O.m-2.s-1)
                                        conductance_mmolH2OMPam2s    ! Effective hydraulic resistance of each layer (mmolH2O.MPa-1.m-2.s-1)
  double precision, dimension(nos_soil_layers+1) :: SWP, & ! soil water potential (MPa)
                                      soil_conductivity, & ! soil conductivity
                                            waterchange, & ! net water change by specific soil layers (m)
                                        water_grav_flow, & ! flow of water under gravity FROM each soil layer (kgH2O/m2/d)
                                         field_capacity, & ! soil field capacity (m3.m-3)
                                 field_capacity_initial, &
                                         soil_waterfrac, & ! soil water content (m3.m-3)
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
                                          Reff, & ! Effective total hydraulic resistance (MPa.m2.s.mmolH2O-1)
                                     max_depth, & ! maximum possible root depth (m)
                                        root_k, & ! biomass to reach half max_depth
                                        runoff, & ! surface water runoff (kgH2O.m-2.day-1)
                                     underflow, & ! drainage from the bottom of soil column (kgH2O.m-2.day-1)
                                previous_depth, & ! depth of bottom of soil profile
                                   canopy_wind, & ! wind speed (m.s-1) at canopy top
                                         ustar, & ! friction velocity (m.s-1)
                                      ustar_Uh, &
                                air_density_kg, & ! air density kg/m3
                                ET_demand_coef, & ! air_density_kg * vpd_kPa * cpair
                                        roughl, & ! roughness length (m)
                                  displacement, & ! zero plane displacement (m)
                                         meant, & ! mean air temperature (oC)
                                         soilT, & ! soil day time temperature
                                         leafT, & ! canopy day time temperature temperature (oC)
                            canopy_swrad_MJday, & ! canopy_absorbed shortwave radiation (MJ.m-2.day-1)
                              canopy_par_MJday, & ! canopy_absorbed PAR radiation (MJ.m-2.day-1)
                                soil_par_MJday, & ! soil absorbed PAR radiation (MJ.m-2.day-1)
                              soil_swrad_MJday, & ! soil absorbed shortwave radiation (MJ.m-2.day-1)
                              canopy_lwrad_Wm2, & ! canopy absorbed longwave radiation (W.m-2)
                                soil_lwrad_Wm2, & ! soil absorbed longwave radiation (W.m-2)
                                 sky_lwrad_Wm2, & ! sky absorbed longwave radiation (W.m-2)
                          stomatal_conductance, & ! canopy scale stomatal conductance (mmolH2O.m-2.d-1)
                         potential_conductance, & ! potential stomatal conductance (mmolH2O.m-2ground.s-1)
                           minimum_conductance, & ! potential stomatal conductance (mmolH2O.m-2ground.s-1)
                       aerodynamic_conductance, & ! aerodynamic conductance at canopy top (m.s-1)
                              soil_conductance, & ! soil surface conductance (m.s-1)
                             convert_ms1_mol_1, & ! Conversion ratio for m.s-1 -> mol.m-2.s-1
                            convert_ms1_mmol_1, & ! Conversion ratio for m/s -> mmol/m2/s
                           air_vapour_pressure, & ! Vapour pressure of the air (kPa)
                                        lambda, & ! latent heat of vapourisation (J.kg-1)
                                         psych, & ! psychrometric constant (kPa K-1)
                                         slope, & ! Rate of change of saturation vapour pressure with temperature (kPa.K-1)
                        water_vapour_diffusion, & ! Water vapour diffusion coefficient in (m2/s)
                           kinematic_viscosity, & ! kinematic viscosity (m2.s-1)
                                  snow_storage, & ! snow storage on soil surface (kgH2O/m2)
                                canopy_storage, & ! water storage on canopy (kgH2O.m-2)
                          intercepted_rainfall    ! intercepted rainfall rate equivalent (kgH2O.m-2.s-1)

  ! Module level variables for ACM_GPP_ET variables
  double precision ::   delta_gs, & ! day length corrected gs increment mmolH2O/m2/day
                             avN, & ! Average canopy foliar N (gN/m2leaf)
                             NUE, & ! Photosynthetic nitrogen use efficiency at optimum temperature (oC)
                                    ! ,unlimited by CO2, light and photoperiod (umolC/gN/m2leaf)
                       iWUE_step, & ! Intrinsic water use efficiency for that day (gC/m2leaf/dayl/mmolH2Ogs)
metabolic_limited_photosynthesis, & ! temperature, leaf area and foliar N limiterd photosynthesis (gC/m2/day)
    light_limited_photosynthesis, & ! light limited photosynthesis (gC/m2/day)
                              ci, & ! Internal CO2 concentration (ppm)
                        rb_mol_1, & ! Canopy boundary layer resistance (day/m2/molCO2)
                     o2_half_sat, & ! O2 at which photorespiration is 50 % of maximum
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
 leaf_canopy_light_scaling, & ! approximate scaling factor from leaf to canopy for gpp and gs
  leaf_canopy_wind_scaling, & ! approximate scaling factor from leaf to canopy for gb
                       lai    ! leaf area index (m2/m2)

  ! Module level variables for step specific timing information
  integer :: steps_per_year
  double precision ::       seconds_per_step, & !
                               days_per_step, & !
                             days_per_step_1, & !
                          mean_days_per_step, &
                                dayl_seconds, & ! day length in seconds
                              dayl_seconds_1, &
                         dayl_hours_fraction, &
                                  dayl_hours    ! day length in hours

  double precision, dimension(:), allocatable :: deltat_1, & ! inverse of decimal days
                                  airt_zero_fraction_time, &
                                          daylength_hours, &
                                        daylength_seconds, &
                                      daylength_seconds_1, &
                                            rainfall_time
  end type
  type(model_working_variables), allocatable, dimension(:):: mVs

contains
  !
  !--------------------------------------------------------------------
  !
  end module carbon_model_memory
