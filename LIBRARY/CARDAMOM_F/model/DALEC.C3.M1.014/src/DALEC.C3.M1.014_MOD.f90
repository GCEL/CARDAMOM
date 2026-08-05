module carbon_model_memory
  implicit none
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
                                pi2 = 9.869604d0,   & ! pi**2d0
                             two_pi = 6.283185d0,   & ! pi*2d0
                         deg_to_rad = 0.01745329d0, & ! pi/180d0
                sin_dayl_deg_to_rad = 0.3979486d0,  & ! sin( 23.45d0 * deg_to_rad )
                             freeze = 273.15d0

  ! photosynthesis / respiration parameters
  double precision, parameter :: &
                        Rg_fraction = 0.21875d0,    & ! fraction of C allocation towards each pool
                                                      ! lost as growth respiration
                                                      ! (i.e. 0.28 .eq. xNPP)
                    one_Rg_fraction = 1d0 - Rg_fraction

  ! timing parameters
  double precision, parameter :: &
                   seconds_per_hour = 3600d0,         & ! Number of seconds per hour
                    seconds_per_day = 86400d0,        & ! Number of seconds per day
                  seconds_per_day_1 = 1.157407d-05      ! Inverse of seconds per day

  double precision, parameter :: resp_rate_temp_coeff = 0.0334798d0,& ! exponential temperature response for heterotrophic respiration (0.0334798 = Q10 of 1.4, 0.0693d0 = Q10 of 2)
                                               lv_res = 0.1d0,      & ! residue fraction of leaves left post harvest
                                               st_res = 0.1d0,      & ! residue fraction of stem left post harvest 
                                                LAICR = 4d0,        & ! LAI above which self shading turnover occurs
                                          rel_gso_max = 0.35d0,     & ! allocation to storage organ relative to GPP
                               resp_cost_labile_trans = 0.21875d0     ! labile lost to respiration per gC labile to GPP


  !!!!!!!!!
  ! Module level variables
  !!!!!!!!!
type model_working_variables
  logical :: do_iWUE = .true. ! Use iWUE or WUE for stomatal optimisation

  ! hydraulic model variables
  integer :: water_retention_pass, soil_layer
  double precision, dimension(nos_soil_layers) :: soil_frac_clay,soil_frac_sand ! clay and soil fractions of soil

  ! Module level variables for ACM_GPP parameters
  double precision ::       ceff, & ! canopy efficency, ceff = avN*NUE
                             avN, & ! average foliar N (gN/m2)
                             NUE, & ! Photosynthetic nitrogen use efficiency at optimum temperature (oC)
                                    ! ,unlimited by CO2, light and photoperiod (gC/gN/m2leaf/day
                              ci    ! Internal CO2 concentration (ppm) 

  ! Module level variables for step specific met drivers
  double precision :: mint, & ! minimum temperature (oC)
                      maxt, & ! maximum temperature (oC)
                     swrad, & ! incoming short wave radiation (MJ/m2/day)
                       co2, & ! CO2 (ppm)
                       doy, & ! Day of year
                     lai_1, & ! inverse of LAI
                       lai    ! leaf area index (m2/m2)

  ! Module level variables for step specific timing and location information
  integer :: steps_per_year
  double precision ::       seconds_per_step, & !
                               days_per_step, & !
                             days_per_step_1, & !
                          mean_days_per_step, &
                                dayl_seconds, & ! day length in seconds
                              dayl_seconds_1, &
                         dayl_hours_fraction, &
                                  dayl_hours, & ! day length in hours
                                    latitude, & ! latitude ()-90/90)
                            latitude_radians, & ! latitude in radians
                        sin_latitude_radians, & ! sin(latitude_radians)
                        cos_latitude_radians, & ! cos(latitude_radians)
                          sunset_solar_angle, & ! Solar angle at sunset hour
                                 declination, & ! Solar declination, function of day of year
                   cosine_solar_zenith_angle    ! Cosine zenith angle of the timestep

  double precision, dimension(:), allocatable :: deltat_1, & ! inverse of decimal days
                                          daylength_hours, &
                                        daylength_seconds, &
                                      daylength_seconds_1

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

  double precision ::             gpp_acm, & ! gross primary productivity (gC.m-2.day-1)
                      stock_storage_organ, & ! storage organ C pool, i.e. the desired crop (gC.m-2)
                       stock_dead_foliage, & ! dead but still standing foliage (gC.m-2)
                          stock_resp_auto, & ! autotrophic respiration pool (gC.m-2)
                             stock_labile, & ! labile C pool (gC.m-2)
                            stock_foliage, & ! foliage C pool (gC.m-2)
                               stock_stem, & ! stem C pool (gC.m-2)
                              stock_roots, & ! roots C pool (gC.m-2)
                             stock_litter, & ! litter C pool (gC.m-2)
                      stock_soilOrgMatter, & ! SOM C pool (gC.m-2)
                                resp_auto, & ! autotrophic respiration (gC.m-2.d-1)
                            resp_h_litter, & ! litter heterotrophic respiration (gC.m-2.d-1)
                     resp_h_soilOrgMatter, & ! SOM heterotrophic respiration (gC.m-2.d-1)
                                      npp, & ! net primary productivity (gC.m-2.d-1)
                                nee_dalec, & ! net ecosystem exchange (gC.m-2.d-1)
                                       DS, & ! Developmental state and initial condition
                                      LCA, & ! leaf mass area (gC.m-2)
               alloc_to_storage_organ_old, & ! rolling average allocation of GPP to storage organ (gC.m-2.d-1)
                       decomposition_rate, & ! decomposition rate (frac / day)
                       frac_GPP_resp_auto, & ! fraction of GPP allocated to autotrophic carbon pool
                    turnover_rate_foliage, & ! turnover rate of foliage (frac/day)
                       turnover_rate_stem, & ! same for stem
                     turnover_rate_labile, & ! same for labile
                  turnover_rate_resp_auto, & ! same for autotrophic C pool
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
                                  doptmin, & ! difference between optimum and minimum cardinal temperatures
                                  dmaxmin, & ! difference between maximum and minimum cardinal temperatures
                                doptmin_v, & ! difference between optimum and minimum vernalisation temperatures
                                dmaxmin_v, & ! difference between maximum and minimum vernalisation temperatures
                                      VDh, & ! effective vernalisation days when plants are 50 % vernalised
                                       VD, & ! count of vernalisation days
                                 RDRSHMAX, & ! maximum rate of self shading turnover (frac/day)
                                     PHCR, & ! critical value of photoperiod for development
                                     PHSC, & ! photoperiod sensitivity
                                     raso, & ! rolling average for alloc to storage organ
                                 max_raso, & ! maximum value for rolling average alloc to storage organ
                                    BM_EX, & ! biomass extracted in addition to the storage organ
                 HARVESTextracted_foliage, & ! Foliage removed by harvest activity
                    HARVESTextracted_stem, & ! Stem removed by harvest activity
            HARVESTextracted_dead_foliage, & ! Dead still standing foliage removed by harvest activity
                  HARVESTextracted_labile, & ! Labile removed by harvest activity
                    HARVESTlitter_foliage, & ! Foliage converted to litter by harvest
                       HARVESTlitter_stem, & ! Stem converted to litter by harvest
               HARVESTlitter_dead_foliage, & ! Dead standing foliage converted to litter by harvest
                  HARVESTlitter_resp_auto, & ! Autotrophic pool converted to litter by harvest
                     HARVESTlitter_labile, & ! Labile converted to litter by harvest
                       PLOUGHlitter_roots, & ! Plough induced litter generation from roots
                                       HI, & ! Harvest index, the ratio of yield to shoot C
                                    yield, & ! crop yield (gC.m-2)
                       alloc_to_resp_auto, & ! amount of carbon to allocate to autotrophic respiration pool
                      turnover_rate_roots, & ! turnover over rate of roots interpolated each time step
                                  gso_max, & !
                             max_raso_old, & !
                                 raso_old, & !
                  resp_cost_labile_to_npp, & ! respiratory cost of moving carbon..from labile to NPP pool
                  resp_cost_npp_to_labile, & ! ..from remaining NPP to labile pool
              resp_cost_foliage_to_labile, & ! ..from foliage to labile pool
                 resp_cost_stem_to_labile, & ! ..from stem to labile pool
                                resp_rate, & ! rate of respiration at given temperature
                                   Cshoot, & !
                                       DR, & !
                          fol_frac_intpol, & !
                         stem_frac_intpol, & !
                                 fP,fT,fV, & !
                        foliage_to_labile, & !
                           stem_to_labile, & !
                         root_frac_intpol, & !
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
                                    RDRSH, & !
                                    RDRDV, & !
                                      RDR

  !
  ! some hardcoded crop parameters
  !


  end type

  type(model_working_variables), allocatable, dimension(:):: mVs

contains
  !
  !--------------------------------------------------------------------
  !
end module carbon_model_memory
