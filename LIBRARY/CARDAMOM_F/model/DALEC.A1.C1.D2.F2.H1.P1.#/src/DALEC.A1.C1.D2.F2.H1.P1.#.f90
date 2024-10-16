
module CARBON_MODEL_MOD

  implicit none

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
  ! See function/subroutine specific comments for exceptions and contributors
  !!!!!!!!!!!

  ! make all private
  private

  ! explicit publics
  public:: CARBON_MODEL     &
           ,soil_frac_clay   &
           ,soil_frac_sand   &
           ,nos_soil_layers  &
           ,extracted_C      &
           ,cica_time        &
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

  !!!!!!!!!
  ! Parameters
  !!!!!!!!!

  ! useful technical parameters
  double precision, parameter:: xacc = 1d-4        & ! accuracy parameter for zbrent bisection proceedure  ! 0.0001
                              ,vsmall = tiny(0d0)*1e3

  integer, parameter:: nos_root_layers = 2, nos_soil_layers = nos_root_layers+1
  double precision, parameter:: pi = 3.1415927d0,    &
                               pi_1 = 0.3183099d0,    & ! pi**(-1d0)
                             two_pi = 6.283185d0,     & ! pi*2d0
                         deg_to_rad = 0.01745329d0,   & ! pi/180d0
                sin_dayl_deg_to_rad = 0.3979486d0,    & ! sin( 23.45d0*deg_to_rad )
                              boltz = 5.670400d-8,    & ! Boltzmann constant (W.m-2.K-4)
                         emissivity = 0.96d0,         &
                        emiss_boltz = 5.443584d-08,   & ! emissivity*boltz
                    sw_par_fraction = 0.5d0,          & ! fraction of short-wave radiation which is PAR
                             freeze = 273.15d0,       &
                         gs_H2O_CO2 = 1.646259d0,     & ! The ratio of H20:CO2 diffusion for gs (Jones appendix 2)
                gs_H2Ommol_CO2mol_day = 142.2368d0,   & ! The ratio of H20:CO2 diffusion for gs, including seconds per day correction
                         gb_H2O_CO2 = 1.37d0,         & ! The ratio of H20:CO2 diffusion for gb (Jones appendix 2)
                   mmol_to_kg_water = 1.8d-5,         & ! milli mole conversion to kg
                         umol_to_gC = 1.2d-5,         & ! conversion of umolC -> gC
                         gC_to_umol = 83333.33d0,     & ! conversion of gC -> umolC; umol_to_gC**(-1d0)
                               Rcon = 8.3144d0,       & ! Universal gas constant (J.K-1.mol-1)
                          vonkarman = 0.41d0,         & ! von Karman's constant
                        vonkarman_1 = 2.439024d0,     & ! 1/von Karman's constant
                              cpair = 1004.6d0          ! Specific heat capacity of air; used in energy balance J.kg-1.K-1

  ! photosynthesis/respiration parameters
  double precision, parameter :: &
                      kc_half_sat_25C = 310d0,        & ! CO2 half saturation, saturation value
                 kc_half_sat_gradient = 23.956d0,     & ! CO2 half sat, half sat
                      co2comp_sat_25C = 36.5d0,       & ! CO2 compensation point, saturation
                     co2comp_gradient = 9.46d0          ! CO2 comp point, half sat
                                                        ! Each of these are temperature sensitivty

  ! hydraulic parameters
  double precision, parameter :: &
                             gplant = 4d0,          & ! plant hydraulic conductivity (mmol m-1 s-1 MPa-1)
                        root_resist = 25d0,         & ! Root resistivity (MPa s g mmol−1 H2O)
                          max_depth = 2d0,          & ! max root depth (m)
                             root_k = 100d0,        & ! root biomass needed to reach 50% depth (gbiomass/m2)
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
                   pn_max_temp = 6.842942d+01,  & ! Maximum daily max temperature for photosynthesis (oC)
                   pn_min_temp = -1d+06      ,  & ! Minimum daily max temperature for photosynthesis (oC)
                   pn_opt_temp = 3.155960d+01,  & ! Optimum daily max temperature for photosynthesis (oC)
                   pn_kurtosis = 1.889026d-01,  & ! Kurtosis of photosynthesis temperature response
!bespoke                   pn_max_temp = 59d0,          & ! Maximum daily max temperature for photosynthesis (oC)
!                   pn_min_temp = -4d0,          & ! Minimum daily max temperature for photosynthesis (oC)
!                   pn_opt_temp = 30d0,          & ! Optimum daily max temperature for photosynthesis (oC)
!                   pn_kurtosis = 0.07d0,        & ! Kurtosis of photosynthesis temperature response
                            e0 = 3.661204d+00,  & ! Quantum yield gC/MJ/m2/day PAR
                minlwp_default = -1.808224d+00,  & ! minimum leaf water potential (MPa)
      soil_iso_to_net_coef_LAI = -2.717467d+00,  & ! Coefficient relating soil isothermal net radiation to net.
                          iWUE = 6.431150d-03,  & ! Intrinsic water use efficiency (gC/m2leaf/day/mmolH2Ogs)
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

  ! forest rotation specific info
  double precision, allocatable, dimension(:):: extracted_C
  ! Metrics on photosynthetic activity
  double precision, dimension(:), allocatable:: gs_demand_supply_ratio, & ! actual:potential stomatal conductance
                                                 gs_total_canopy, &        ! stomatal conductance (mmolH2O/m2ground/s)
                                                 gb_total_canopy, &        ! boundary conductance (mmolH2O/m2ground/s)
                                                       cica_time, &        ! Internal vs ambient CO2 concentrations
                                           canopy_par_MJday_time           ! Absorbed PAR by canopy (MJ/m2ground/day)

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

  ! hydraulic model variables
  integer:: water_retention_pass, soil_layer
  double precision, dimension(nos_soil_layers+1):: layer_thickness    ! thickness of soil layers (m)
  double precision, dimension(nos_soil_layers):: soil_frac_clay, soil_frac_sand  ! clay and soil fractions of soil
  double precision, dimension(nos_root_layers):: uptake_fraction, & ! fraction of water uptake from each root layer
                                                           demand, & ! maximum potential canopy hydraulic demand
                                            water_flux_mmolH2Om2s    ! potential transpiration flux (mmolH2O.m-2.s-1)

  double precision:: root_reach, root_biomass, &
                              total_water_flux, & ! potential transpiration flux (kgH2O.m-2.day-1)
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
                                        lambda, & ! latent heat of vapourisation (J.kg-1)
                                         psych, & ! psychrometric constant (kPa K-1)
                                         slope, & ! Rate of change of saturation vapour pressure
                                                  ! with temperature (kPa.K-1)
                        water_vapour_diffusion, & ! Water vapour diffusion coefficient in (m2/s)
                             dynamic_viscosity, & ! dynamic viscosity (kg.m-2.s-1)
                           kinematic_viscosity    ! kinematic viscosity (m2.s-1)

  ! Module level variables for ACM_GPP_ET parameters
  double precision ::   delta_gs, & ! day length corrected gs increment mmolH2O/m2/dayl
                            Ceff, & ! Canopy efficiency (gC/m2leaf/day)
                       iWUE_step, & ! Intrinsic water use efficiency for that day (gC/m2leaf/dayl/mmolH2Ogs)
metabolic_limited_photosynthesis, &
    light_limited_photosynthesis, &
                              ci, & ! Internal CO2 concentration (ppm)
                          gb_mol, & ! Canopy boundary layer conductance (molCO2/m2/day)
                        rb_mol_1, & ! Canopy boundary layer resistance (day/m2/molCO2)
                 pn_airt_scaling, & ! temperature response for metabolic limited photosynthesis
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
                     lai_1, & ! inverse of LAI
                       lai    ! leaf area index (m2/m2)

  ! Module level varoables for step specific timing information
  double precision:: cos_solar_zenith_angle, &
                            seconds_per_step, & !
                               days_per_step, & !
                         dayl_hours_fraction, &
                                dayl_seconds, & ! day length in seconds
                              dayl_seconds_1, &
                                  dayl_hours    ! day length in hours

  contains
  !
  !--------------------------------------------------------------------
  !
  subroutine CARBON_MODEL(start, finish, met, pars, deltat, nodays, lat, lai_out, NEE, FLUXES, POOLS &
                         ,nopars, nomet, nopools, nofluxes, GPP)

    ! The Data Assimilation Linked Ecosystem Carbon-Combined Deciduous
    ! Evergreen Analytical-ACMv2(DALEC_CDEA_ACM2) model. The subroutine calls the
    ! Aggregated Canopy Model version 2 to simulate GPP and partitions
    ! between various ecosystem carbon pools. These pools are subject
    ! to turnovers/decompostion resulting in ecosystem phenology and fluxes of CO2

    ! This version includes the option to simulate fire combustion based
    ! on burned fraction and fixed combusion rates. It also includes the
    ! possibility to remove a fraction of biomass to simulate deforestation.

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
                         ,lat                 ! site latitude (degrees)

    double precision, dimension(nodays), intent(inout):: lai_out & ! leaf area index
                                                             ,GPP & ! Gross primary productivity
                                                             ,NEE   ! net ecosystem exchange of CO2

    double precision, dimension((nodays+1), nopools), intent(inout):: POOLS  ! vector of ecosystem pools

    double precision, dimension(nodays, nofluxes), intent(inout):: FLUXES  ! vector of ecosystem fluxes

    ! declare local variables
    double precision:: Rtot & ! Total hydraulic resistance (MPa.s-1.m-2.mmol-1)
                       ,tmp, wf, wl, ff, fl, osf, osl, sf, ml   ! phenological controls

    ! JFE added 4 May 2018-combustion efficiencies and fire resilience
    double precision:: cf(6), rfac(6), burnt_area
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
                       ,labile_frac_removal          &
                       ,Cstem, Crootcr, stem_residue   &
                       ,coarse_root_residue          &
                       ,soil_loss_with_roots
    integer:: harvest_management, p, f, n
    real:: start_time, finish_time

    ! met drivers are:
    ! 1st run day
    ! 2nd min daily temp (oC)
    ! 3rd max daily temp (oC)
    ! 4th Radiation (MJ.m-2.day-1)
    ! 5th CO2 (ppm)
    ! 6th DOY at end of time step
    ! 7th Not used
    ! 8th removed fraction
    ! 9th burned fraction

    ! POOLS are:
    ! 1 = labile
    ! 2 = foliar
    ! 3 = root
    ! 4 = wood
    ! 5 = litter
    ! 6 = som

    ! FLUXES are:
    ! 1 = GPP
    ! 2 = temprate
    ! 3 = respiration_auto
    ! 4 = leaf production
    ! 5 = labile production
    ! 6 = root production
    ! 7 = wood production
    ! 8 = labile production
    ! 9 = leaffall factor
    ! 10 = leaf litter production
    ! 11 = woodlitter production
    ! 12 = rootlitter production
    ! 13 = respiration het litter
    ! 14 = respiration het som
    ! 15 = litter2som
    ! 16 = labrelease factor

    ! JFE added 03/05/2018-start
    ! emissions of carbon into the atmosphere due to combustion
    ! 17 = ecosystem fire emission  (sum of fluxes 18 to 23)
    ! 18 = fire emission from labile
    ! 19 = fire emission from foliar
    ! 20 = fire emission from roots
    ! 21 = fire emission from wood
    ! 22 = fire emission from litter
    ! 23 = fire emission from soil

    ! mortality due to fire
    ! 24 = transfer from labile into litter
    ! 25 = transfer from foliar into litter
    ! 26 = transfer from roots into litter
    ! 27 = transfer from wood into som
    ! 28 = transfer from litter into som
    ! JFE added 03/05/2018-stop


    ! PARAMETERS
    ! 17 values

    ! p(1) Litter to SOM conversion rate  - m_r
    ! p(2) Fraction of GPP respired-f_a
    ! p(3) Fraction of NPP allocated to foliage-f_f
    ! p(4) Fraction of NPP allocated to roots-f_r
    ! p(5) Leaf lifespan-L_f
    ! p(6) Turnover rate of wood-t_w
    ! p(7) Turnover rate of roots-t_r
    ! p(8) Litter turnover rate-t_l
    ! p(9) SOM turnover rate  - t_S
    ! p(10) Parameter in exponential term of temperature - \theta
    ! p(11) Canopy efficiency parameter-C_eff (part of ACM)
    ! p(12) = date of Clab release-B_day
    ! p(13) = Fraction allocated to Clab-f_l
    ! p(14) = lab release duration period-R_l
    ! p(15) = date of leaf fall-F_day
    ! p(16) = leaf fall duration period-R_f
    ! p(17) = LMA

    ! Reset all POOLS and FLUXES to prevent precision errors
    FLUXES = 0d0; POOLS = 0d0

    ! load ACM-GPP-ET parameters
    Ceff = pars(11)  ! Canopy efficiency (gC/m2/day)
                    ! This is in the full model the product of Nitrogen use efficiency (gC/gN/m2leaf/day)
                    ! and average foliar nitrogen gC/m2leaf

    ! assigning initial conditions
    if (start == 1) then
        POOLS(1, 1) = pars(18)  ! labile
        POOLS(1, 2) = pars(19)  ! foliar
        POOLS(1, 3) = pars(20)  ! roots
        POOLS(1, 4) = pars(21)  ! wood
        POOLS(1, 5) = pars(22)  ! litter
        POOLS(1, 6) = pars(23)  ! som
    endif

    ! defining phenological variables
    ! release period coefficient, based on duration of labile turnover or leaf
    ! fall durations
    wf = pars(16)*sqrt(2d0) * 0.5d0
    wl = pars(14)*sqrt(2d0) * 0.5d0

    ! magnitude coefficient
    ff = (log(pars(5))-log(pars(5)-1d0)) * 0.5d0
    fl = (log(1.001d0)-log(0.001d0)) * 0.5d0

    ! set minium labile life span to one year
    ml = 1.001d0

    ! offset for labile and leaf turnovers
    osf = ospolynomial(pars(5), wf)
    osl = ospolynomial(ml, wl)

    ! scaling to biyearly sine curve
    sf = 365.25d0/pi

    ! if either of our disturbance drivers indicate disturbance will occur then
    ! set up these components
    if (maxval(met(8, :)) > 0d0 .or. maxval(met(9, :)) > 0d0) then

        ! now load the hardcoded forest management parameters into their scenario locations

        ! Deforestation process functions in a sequenctial way.
        ! Thus, the pool_loss is first determined as a function of met(8, n) and
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

        ! Declare combustion efficiency (labile, foliar, roots, wood, litter, soil, woodlitter)
        cf(1) = 0.1d0; cf(2) = 0.9d0
        cf(3) = 0.1d0; cf(4) = 0.1d0
        cf(5) = 0.7d0; cf(6) = 0.01d0
        ! Resilience factor for non-combusted tissue
        rfac = 0.5d0; rfac(5) = 0.1d0; rfac(6) = 0d0

        ! JFE added 4 May 2018-define fire constants
        ! Update fire parameters derived from
        ! Yin et al., (2020), doi: 10.1038/s414647-020-15852-2
        ! Subsequently expanded by T. L. Smallman & Mat Williams (UoE, 03/09/2021)
        ! to provide specific CC for litter and wood litter.
        ! NOTE: changes also result in the addition of further EDCs

        ! Assign proposed resilience factor
        rfac(1:4) = pars(24)
        rfac(5) = 0.1d0; rfac(6) = 0d0
        ! Assign combustion completeness to foliage
        cf(2) = pars(25)  ! foliage
        ! Assign combustion completeness to non-photosynthetic
        cf(1) = pars(26); cf(3) = pars(26); cf(4) = pars(26)
        cf(6) = pars(27)  ! soil
        ! derived values for litter
        cf(5) = pars(28)

    end if  ! disturbance ?

    !
    ! Begin looping through each time step
    !

    ! initialise root reach based on initial conditions
    root_biomass = max(min_root, POOLS(1, 3)*2d0)
    ! needed to initialise soils
    call calculate_Rtot(Rtot)

    ! Allocate variables needed
    if (.not.allocated(gs_demand_supply_ratio)) then
        allocate(gs_demand_supply_ratio(nodays), &
                 gs_total_canopy(nodays), gb_total_canopy(nodays), &
                 canopy_par_MJday_time(nodays), cica_time(nodays))
    end if

    do n = start, finish

       !!!!!!!!!!
       ! assign drivers and update some prognostic variables
       !!!!!!!!!!

       ! Incoming drivers
       mint = met(2, n)  ! minimum temperature (oC)
       maxt = met(3, n)  ! maximum temperature (oC)
       leafT = (maxt*0.75d0) + (mint*0.25d0)   ! initial day time canopy temperature (oC)
       swrad = met(4, n)  ! incoming short wave radiation (MJ/m2/day)
       co2 = met(5, n)   ! CO2 (ppm)
       doy = met(6, n)   ! Day of year
       meant = (mint+maxt) * 0.5d0  ! mean air temperature (oC)
       wind_spd = met(15, n)  ! wind speed (m/s)
       vpd_kPa = met(16, n)*1d-3  ! Vapour pressure deficit (Pa -> kPa)

       ! calculate LAI value
       lai_out(n) = POOLS(n, 2)/pars(17)
       lai = lai_out(n)

       ! extract timing related values
       call calculate_daylength((doy-(deltat(n)*0.5d0)), lat)
       dayl_hours_fraction = dayl_hours*0.04166667d0  ! 1/24 = 0.04166667
       iWUE_step = iWUE*dayl_hours_fraction
       dayl_seconds_1 = dayl_seconds ** (-1d0)
       seconds_per_step = seconds_per_day*deltat(n)
       days_per_step = deltat(n)

       !!!!!!!!!!
       ! Calculate surface exchange coefficients
       !!!!!!!!!!

       ! calculate some temperature dependent meteorologial properties
       call meteorological_constants(maxt, maxt+freeze, vpd_kPa)
       ! calculate aerodynamic using consistent approach with SPA
       call calculate_aerodynamic_conductance
       ! Aerodynamic conductance (mmolH2O/m2ground/s)
       gb_total_canopy(n) = aerodynamic_conductance*convert_ms1_mol_1*1d3

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
       root_biomass = max(min_root, POOLS(n, 3)*2d0)
       call calculate_Rtot(Rtot)

       ! calculate variables used commonly between ACM_GPP and ACM_ET
       call calculate_stomatal_conductance
       ! Estimate stomatal conductance relative to its minimum/maximum, i.e. how
       ! close are we to maxing out supply (note 0.01 taken from min_gs)
       gs_demand_supply_ratio(n) = (stomatal_conductance-minimum_conductance) / (potential_conductance-minimum_conductance)
       ! Store the canopy level stomatal conductance (mmolH2O/m2/s)
       gs_total_canopy(n) = stomatal_conductance

       !!!!!!!!!!
       ! GPP (gC.m-2.day-1)
       !!!!!!!!!!

       if (stomatal_conductance > vsmall) then
           call acm_gpp_stage_1; FLUXES(n, 1) = acm_gpp_stage_2(stomatal_conductance)
           cica_time(n) = ci/co2
       else
           FLUXES(n, 1) = 0d0; cica_time(n) = 0d0
       endif

       ! temprate (i.e. temperature modified rate of metabolic activity))
       FLUXES(n, 2) = exp(pars(10)*0.5d0*(met(3, n)+met(2, n)))
       ! autotrophic respiration (gC.m-2.day-1)
       FLUXES(n, 3) = pars(2)*FLUXES(n, 1)
       ! leaf production rate (gC.m-2.day-1)
       FLUXES(n, 4) = (FLUXES(n, 1)-FLUXES(n, 3))*pars(3)
       ! labile production (gC.m-2.day-1)
       FLUXES(n, 5) = (FLUXES(n, 1)-FLUXES(n, 3)-FLUXES(n, 4))*pars(13)
       ! root production (gC.m-2.day-1)
       FLUXES(n, 6) = (FLUXES(n, 1)-FLUXES(n, 3)-FLUXES(n, 4)-FLUXES(n, 5))*pars(4)
       ! wood production
       FLUXES(n, 7) = FLUXES(n, 1)-FLUXES(n, 3)-FLUXES(n, 4)-FLUXES(n, 5)-FLUXES(n, 6)

       ! Labile release and leaffall factors
       FLUXES(n, 9) = (2d0/sqrt(pi))*(ff/wf)*exp(-(sin((doy-pars(15)+osf)/sf)*sf/wf)**2d0)
       FLUXES(n, 16) = (2d0/sqrt(pi))*(fl/wl)*exp(-(sin((doy-pars(12)+osl)/sf)*sf/wl)**2d0)

       !
       ! those with time dependancies
       !

       ! total labile release
       FLUXES(n, 8) = POOLS(n, 1)*(1d0-(1d0-FLUXES(n, 16))**deltat(n))/deltat(n)
       ! total leaf litter production
       FLUXES(n, 10) = POOLS(n, 2)*(1d0-(1d0-FLUXES(n, 9))**deltat(n))/deltat(n)
       ! total wood production
       FLUXES(n, 11) = POOLS(n, 4)*(1d0-(1d0-pars(6))**deltat(n))/deltat(n)
       ! total root litter production
       FLUXES(n, 12) = POOLS(n, 3)*(1d0-(1d0-pars(7))**deltat(n))/deltat(n)

       !
       ! those with temperature AND time dependancies
       !

       ! respiration heterotrophic litter
       FLUXES(n, 13) = POOLS(n, 5)*(1d0-(1d0-FLUXES(n, 2)*pars(8))**deltat(n))/deltat(n)
       ! respiration heterotrophic som
       FLUXES(n, 14) = POOLS(n, 6)*(1d0-(1d0-FLUXES(n, 2)*pars(9))**deltat(n))/deltat(n)
       ! litter to som
       FLUXES(n, 15) = POOLS(n, 5)*(1d0-(1d0-pars(1)*FLUXES(n, 2))**deltat(n))/deltat(n)

       ! calculate the NEE
       NEE(n) = (-FLUXES(n, 1)+FLUXES(n, 3)+FLUXES(n, 13)+FLUXES(n, 14))
       ! load GPP
       GPP(n) = FLUXES(n, 1)

       !
       ! update pools for next timestep
       !

       ! labile pool
       POOLS(n+1, 1) = POOLS(n, 1) + (FLUXES(n, 5)-FLUXES(n, 8))*deltat(n)
       ! foliar pool
       POOLS(n+1, 2) = POOLS(n, 2) + (FLUXES(n, 4)-FLUXES(n, 10) + FLUXES(n, 8))*deltat(n)
       ! wood pool
       POOLS(n+1, 4) = POOLS(n, 4) + (FLUXES(n, 7)-FLUXES(n, 11))*deltat(n)
       ! root pool
       POOLS(n+1, 3) = POOLS(n, 3) + (FLUXES(n, 6) - FLUXES(n, 12))*deltat(n)
       ! litter pool
       POOLS(n+1, 5) = POOLS(n, 5) + (FLUXES(n, 10)+FLUXES(n, 12)-FLUXES(n, 13)-FLUXES(n, 15))*deltat(n)
       ! som pool
       POOLS(n+1, 6) = POOLS(n, 6) + (FLUXES(n, 15)-FLUXES(n, 14)+FLUXES(n, 11))*deltat(n)

       !!!!!!!!!!
       ! Extract biomass-e.g. deforestation/degradation
       !!!!!!!!!!

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
               POOLS(n+1, 6) = POOLS(n+1, 6) - soil_loss_with_roots+wood_residue
               ! mass balance check
               where (POOLS(n+1, 1:6) < 0d0) POOLS(n+1, 1:6) = 0d0

               ! Convert harvest related extractions to daily rate for output
               ! For dead organic matter pools, in most cases these will be zeros.
               ! But these variables allow for subseqent management where surface litter
               ! pools are removed or mechanical extraction from soil occurs.
               FLUXES(n, 30) = (labile_loss-labile_residue) / deltat(n)  ! Labile extraction
               FLUXES(n, 31) = (foliar_loss-foliar_residue) / deltat(n)  ! foliage extraction
               FLUXES(n, 32) = (roots_loss-roots_residue) / deltat(n)    ! fine roots extraction
               FLUXES(n, 33) = (wood_loss-wood_residue) / deltat(n)      ! wood extraction
               FLUXES(n, 34) = 0d0  ! litter extraction
               FLUXES(n, 35) = soil_loss_with_roots/deltat(n)          ! som extraction
               ! Convert harvest related residue generations to daily rate for output
               FLUXES(n, 36) = labile_residue/deltat(n)  ! labile residues
               FLUXES(n, 37) = foliar_residue/deltat(n)  ! foliage residues
               FLUXES(n, 38) = roots_residue/deltat(n)  ! fine roots residues
               FLUXES(n, 39) = wood_residue/deltat(n)   ! wood residues
               ! Total C extraction, including any potential litter and som.
               FLUXES(n, 29) = sum(FLUXES(n, 30:35))

           end if  ! C_total > 0d0

       endif  ! end deforestation info


! TLS: a modified version of the original model which has now been removed 31/03/2022
!      The replacement provides different scenarios of what is being extracted (e.g. above vs below)
!      and the re-allocation of C to be either extracted or litter/residues remaining in system
!       ! Remove biomass if necessary
!       if (met(8, n) > 0d0) then
!           tmp = (POOLS(n+1, 2)+POOLS(n+1, 4))/sum(POOLS(n+1, 2:4))
!           if (allocated(extracted_C)) then
!               extracted_C(n) = (((POOLS(n+1, 1)*tmp) + POOLS(n+1, 2) + POOLS(n+1, 4)) * met(8, n)) / deltat(n)
!           endif
!           POOLS(n+1, 1) = tmp &
!                        * POOLS(n+1, 1)*(1d0-met(8, n))  ! remove labile
!!           POOLS(n+1, 1) = max(pars(18), tmp &
!!                        * POOLS(n+1, 1)*(1d0-met(8, n)))  ! remove labile
!           POOLS(n+1, 2) = POOLS(n+1, 2)*(1d0-met(8, n))  ! remove foliar
!           POOLS(n+1, 4) = POOLS(n+1, 4)*(1d0-met(8, n))  ! remove wood
!           ! NOTE: fine root is left in system this is an issue...
!       end if

       !!!!!!!!!!
       ! Impose fire
       !!!!!!!!!!

       if (met(9, n) > 0d0 .or.(met(8, n) > 0d0 .and. harvest_management > 0)) then

           ! Adjust burnt area to account for the managment decisions which may not be
           ! reflected in the burnt area drivers
           burnt_area = met(9, n)
           if (met(8, n) > 0d0 .and. burnt_area > 0d0) then
               ! pass harvest management to local integer
               burnt_area = min(1d0, burnt_area+post_harvest_burn(harvest_management))
           else if (met(8, n) > 0d0 .and. burnt_area <= 0d0) then
               burnt_area = post_harvest_burn(harvest_management)
           endif

           ! Determine the corrected burnt area
           if (burnt_area > 0d0) then

               ! first calculate combustion/emissions fluxes in g C m-2 d-1
               FLUXES(n, 18) = POOLS(n+1, 1)*burnt_area*cf(1)/deltat(n)  ! labile
               FLUXES(n, 19) = POOLS(n+1, 2)*burnt_area*cf(2)/deltat(n)  ! foliar
               FLUXES(n, 20) = POOLS(n+1, 3)*burnt_area*cf(3)/deltat(n)  ! roots
               FLUXES(n, 21) = POOLS(n+1, 4)*burnt_area*cf(4)/deltat(n)  ! wood
               FLUXES(n, 22) = POOLS(n+1, 5)*burnt_area*cf(5)/deltat(n)  ! litter
               FLUXES(n, 23) = POOLS(n+1, 6)*burnt_area*cf(6)/deltat(n)  ! som

               ! second calculate litter transfer fluxes in g C m-2 d-1, all pools except som
               FLUXES(n, 24) = POOLS(n+1, 1)*burnt_area*(1d0-cf(1))*(1d0-rfac(1))/deltat(n)  ! labile into litter
               FLUXES(n, 25) = POOLS(n+1, 2)*burnt_area*(1d0-cf(2))*(1d0-rfac(2))/deltat(n)  ! foliar into litter
               FLUXES(n, 26) = POOLS(n+1, 3)*burnt_area*(1d0-cf(3))*(1d0-rfac(3))/deltat(n)  ! roots into litter
               FLUXES(n, 27) = POOLS(n+1, 4)*burnt_area*(1d0-cf(4))*(1d0-rfac(4))/deltat(n)  ! wood into som
               FLUXES(n, 28) = POOLS(n+1, 5)*burnt_area*(1d0-cf(5))*(1d0-rfac(5))/deltat(n)  ! litter into som

               ! update pools-first remove burned vegetation
               POOLS(n+1, 1) = POOLS(n+1, 1) - (FLUXES(n, 18) + FLUXES(n, 24)) * deltat(n)  ! labile
               POOLS(n+1, 2) = POOLS(n+1, 2) - (FLUXES(n, 19) + FLUXES(n, 25)) * deltat(n)  ! foliar
               POOLS(n+1, 3) = POOLS(n+1, 3) - (FLUXES(n, 20) + FLUXES(n, 26)) * deltat(n)  ! roots
               POOLS(n+1, 4) = POOLS(n+1, 4) - (FLUXES(n, 21) + FLUXES(n, 27)) * deltat(n)  ! wood
               ! update pools-add litter transfer
               POOLS(n+1, 5) = POOLS(n+1, 5) + (FLUXES(n, 24) + FLUXES(n, 25) + FLUXES(n, 26) - FLUXES(n, 22) - FLUXES(n, 28)) * deltat(n)
               POOLS(n+1, 6) = POOLS(n+1, 6) + (FLUXES(n, 27) + FLUXES(n, 28) - FLUXES(n, 23)) * deltat(n)

               ! calculate ecosystem emissions (gC/m2/day)
               FLUXES(n, 17) = FLUXES(n, 18)+FLUXES(n, 19)+FLUXES(n, 20)+FLUXES(n, 21)+FLUXES(n, 22)+FLUXES(n, 23)

           end if  ! Burned_area > 0
       else
           ! set fluxes to zero
           FLUXES(n, 17:28) = 0d0
       end if

    end do  ! nodays loop

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
    metabolic_limited_photosynthesis = gC_to_umol*lai*ceff*opt_max_scaling(pn_max_temp, pn_min_temp, pn_opt_temp, pn_kurtosis, leafT)

    !
    ! Light limited photosynthesis
    !

    ! calculate light limted rate of photosynthesis (gC.m-2.day-1)
    light_limited_photosynthesis = e0*canopy_par_MJday

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

    !!!!!!!!!!
    ! Optimise intrinsic water use efficiency
    !!!!!!!!!!

    ! Determine impact of gs increment on pd and how far we are from iWUE
    find_gs_iWUE = iWUE_step - ((acm_gpp_stage_2(gs_in+delta_gs) - acm_gpp_stage_2(gs_in))*lai_1)

    ! Remember to return back to the user
    return

  end function find_gs_iWUE
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_stomatal_conductance

    ! Determines 1) an approximation of canopy conductance (gc) mmolH2O.m-2.s-1
    ! based on potential hydraulic flow, air temperature and absorbed radiation.
    ! 2) calculates absorbed shortwave radiation (W.m-2) as function of LAI
    
    use brent_zero, only: zbrent

    implicit none

    ! local variables
    double precision:: denom, iWUE_lower, iWUE_upper
    double precision, parameter:: max_gs = 2000d0, &  ! mmolH2O.m-2.s-1 (leaf area)
                                   min_gs = 1d0, &     ! mmolH2O.m-2.s-1 (leaf area)
                                   tol_gs = 10d0       ! mmolH2O.m-2.s-1 (leaf area)

    !!!!!!!!!!
    ! Calculate stomatal conductance under H2O and CO2 limitations
    !!!!!!!!!!

    if (aerodynamic_conductance > vsmall .and. total_water_flux > vsmall) then

        ! Determine potential water flow rate (mmolH2O.m-2.dayl-1)
        max_supply = total_water_flux*seconds_per_day
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
        potential_conductance = potential_conductance*convert_ms1_mol_1*1d3
        ! if conditions are dew forming then set conductance to maximum as we
        ! are not going to be limited by water demand
        if (potential_conductance <= 0d0 .or. potential_conductance > max_gs*lai) potential_conductance = max_gs*lai

        ! If there is a positive demand for water then we will solve for
        ! photosynthesis limits on gs through iterative solution
        delta_gs = 1d0*lai  ! mmolH2O/m2leaf/s
        ! Estimate inverse of LAI to avoid division in optimisation
        lai_1 = lai**(-1d0)
        ! Calculate stage one acm, temperature and light limitation which
        ! are independent of stomatal conductance effects
        call acm_gpp_stage_1

!        if (do_iWUE) then
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

    ! Both length scale and mixing length are considered to be constant within
    ! the canopy (under dense canopy conditions) calculate length scale (lc)
    ! for momentum absorption within the canopy; Harman & Finnigan (2007)
    ! and mixing length (lm) for vertical momentum within the canopy Harman & Finnigan (2008)
    length_scale_momentum = (4d0*canopy_height) / local_lai
    mixing_length_momentum = 2d0*(ustar_Uh**3)*length_scale_momentum

    ! based on Harman & Finnigan (2008); neutral conditions only
    call log_law_decay

    ! Rather than the wind at the top of the canopy we are interested in the within canopy wind speed, 
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
    double precision, parameter:: clump = 1d0    & ! Clumping factor (1 = uniform, 0 totally clumped, mean = 0.75)
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
    double precision, intent(in):: local_lai  ! LAI which has had minium LAI constraint applied
    ! local variables
    double precision  sqrt_cd1_lai
    double precision, parameter:: cd1 = 7.5d0,   & ! Canopy drag parameter; fitted to data
                                    Cs = 0.003d0, & ! Substrate drag coefficient
                                    Cr = 0.3d0,   & ! Roughness element drag coefficient
                          ustar_Uh_max = 0.35d0,  &
                          ustar_Uh_min = 0.05d0,  &
                                    Cw = 2d0,     & ! Characterises roughness sublayer depth (m)
                                 phi_h = 0.19314718056d0  ! Roughness sublayer influence function; 

    ! Estimate LAI scaled canopy drag factor
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
  !------------------------------------------------------------------
  !
  double precision function root_resistance (root_mass, thickness)

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
  !------------------------------------------------------------------
  !
  double precision function ospolynomial(L, w)

    ! Function calculates the day offset for Labile release and leaf turnover
    ! functions

    implicit none

    ! declare input variables
    double precision, intent(in) ::  L, w  ! polynomial coefficients and scaling factor

    ! declare local variables
    double precision ::  tmp, LLog, mxc(7)  ! polynomial coefficients and scaling factor

    ! assign polynomial terms
    mxc(1) = (0.000023599784710d0)
    mxc(2) = (0.000332730053021d0)
    mxc(3) = (0.000901865258885d0)
    mxc(4) = (-0.005437736864888d0)
    mxc(5) = (-0.020836027517787d0)
    mxc(6) = (0.126972018064287d0)
    mxc(7) = (-0.188459767342504d0)

    ! load log of leaf/labile turnovers
    LLog = log(L-1d0)

    ! calculate the polynomial function
    ospolynomial = (mxc(1)*LLog**6d0+mxc(2)*LLog**5d0 + &
                    mxc(3)*LLog**4d0+mxc(4)*LLog**3d0 + &
                    mxc(5)*LLog**2d0+mxc(6)*LLog      + mxc(7))*w

    ! back to the user...
    return

  end function ospolynomial
end module CARBON_MODEL_MOD
