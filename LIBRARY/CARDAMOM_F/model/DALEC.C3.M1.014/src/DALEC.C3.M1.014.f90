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
! This file contains the source code of DALEC.C3.M1.014
!
! This code contains a variant of the Data Assimilation Linked ECosystem (DALEC) model.
! This version of DALEC is derived from the following primary references:
! Sus et al., (2010), https://doi.org/10.1016/j.agee.2010.06.012.
! Bloom & Williams (2015), https://doi.org/10.5194/bg-12-1299-2015.
! Smallman et al., (2017), https://doi.org/10.1002/2016JG003520.
! This code is based on that created by A. A. Bloom (UoE, now at JPL, USA) and
! Oliver Sus (UoE, now at EUMETSAT, Darmstadt).
! Subsequent modifications by:
! T. L. Smallman (University of Edinburgh, t.l.smallman@ed.ac.uk)
! See function / subroutine specific comments for exceptions and contributors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module CARBON_MODEL_MOD

  implicit none

  ! make all private
  private

  ! explicit publics
  public :: CARBON_MODEL  &
           ,initialize_mv

  use carbon_model_memory

  contains
  !
  !--------------------------------------------------------------------
  !
  subroutine initialize_mv(mV, nodays, nomet, nopars, deltat, soil_frac_sand, soil_frac_clay, met, lat)

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

    allocate(mV%deltat_1(nodays))
    mV%deltat_1 = deltat**(-1d0)
    ! zero variables not done elsewhere
    mV%water_flux_mmolH2Om2s = 0d0
    ! initialise some time invarient parameters
    call saxton_parameters(mV%soil_frac_clay,mV%soil_frac_sand, mV)
    call initialise_soils(mV%soil_frac_clay,mV%soil_frac_sand, mV)
    ! call update_soil_initial_conditions(pars(24), mV)
    ! save the initial conditions for later
    mV%field_capacity_initial = mV%field_capacity
    mV%porosity_initial = mV%porosity
   
  end subroutine initialize_mv
  !
  !--------------------------------------------------------------------
  !
  subroutine destroy_mv(mV)
    !! deallocate arrays in mV
    type(model_working_variables):: mV
        deallocate(mV%deltat_1)
  end subroutine
  !
  !--------------------------------------------------------------------
  !
  subroutine CARBON_MODEL(start,finish,met,pars,deltat,nodays,lat &
                         ,FLUXES,POOLS,DIAGS & 
                         ,nopars,nomet,nopools,nofluxes,nodiags &
                         ,stock_seed_labile,DS_shoot,DS_root,fol_frac    &
                         ,stem_frac,root_frac,DS_LRLV,LRLV,DS_LRRT,LRRT, mV)

    !
    ! The Data Assimilation Linked Ecosystem Carbon - CROP - BUCKET (DALEC.C3.M1.014) model.
    ! modified from Sus et al., (2010)
    !
    ! The Aggregated Canopy Model for Gross Primary Productivity (Williams et al., 1997)
    !
    ! This version was coded by T. Luke Smallman (t.l.smallman@ed.ac.uk)
    ! Version 1.0: 15/07/2014
    ! Version 2.0: 25/04/2024 - Integrated into CARDAMOM as substraction of BUCKET model

    use cardamom_structures, only: CI ! get former arguments: stock_seed_labile,DS_shoot,DS_root,fol_frac    &
                         ! ,stem_frac,root_frac,DS_LRLV,LRLV,DS_LRRT,LRRT, from here (read-only) 

    implicit none

      type(model_working_variables) :: mV

    ! declare input variables
    integer, intent(in) :: start    &
                          ,finish   &
                          ,nopars     & ! number of parameters in vector
                          ,nomet      & ! number of meteorological fields
                          ,nofluxes   & ! number of model fluxes
                          ,nopools    & ! number of model pools
                          ,nodays     & ! number of days in simulation
                          ,nodiags      ! number of model diagnostics

    double precision, intent(in) :: met(nomet,nodays)   & ! met drivers
                         ,stock_seed_labile             & ! seed carbon to get things going
                         ,deltat(nodays)                & ! time step in decimal days
                         ,pars(nopars)                  & ! number of parameters
                         ,lat                             ! site latitude (degrees)

    double precision, dimension(:), intent(inout) ::          DS_shoot, & !
                                                               DS_root, & !
                                                              fol_frac, & !
                                                             stem_frac, & !
                                                             root_frac, & !
                                                               DS_LRLV, & !
                                                                  LRLV, & !
                                                               DS_LRRT, & !
                                                                  LRRT    !

    double precision, dimension((nodays+1),nopools), intent(inout) :: POOLS ! vector of ecosystem pools
    double precision, dimension(nodays,nofluxes), intent(inout) :: FLUXES ! vector of ecosystem fluxes
    double precision, dimension(nodays,nodiags), intent(inout) :: DIAGS ! vector of ecosystem diagnostics

    ! declare local variables
    double precision :: airt_weighting(3) &
                         ,acm_forcings(7) & ! ACM inputs (LAI+met)
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
    ! 1 = labile                      (gC/m2) (initial value: p18)
    ! 2 = foliar                      (gC/m2) (initial value: p19)
    ! 3 = root                        (gC/m2) (initial value: p20)
    ! 4 = stem                        (gC/m2) (initial value: p21)
    ! 5 = litter                      (gC/m2) (initial value: p22)
    ! 6 = som                         (gC/m2) (initial value: p23)
    ! 7 = autotrophic                 (gC/m2) (initial value: p24)
    ! 8 = storage organ C             (gC/m2) (initial value: p25)
    ! 9 = dead still standing foliage (gC/m2) (initial value: 0)

    ! FLUXES are:
    ! 1  = GPP (gC/m2/day)
    ! 2  = temperature rate modifier (unitless)
    ! 3  = total autotrophic respiration - maintenance + growth (gC/m2/day)
    ! 4  = GPP allocation to foliage (gC/m2/day)
    ! 5  = GPP allocation to labile and tissue remobilisation into labile (gC/m2/day)
    ! 6  = GPP allocation to roots (gC/m2/day)
    ! 7  = GPP allocation to stem (gC/m2/day)
    ! 8  = labile->NPP transfer (gC/m2/day)
    ! 9  = allocation to storage organ (gC/m2/day)
    ! 10 = leaf litter production (gC/m2/day)
    ! 11 = stem litter production (gC/m2/day)
    ! 12 = root litter production (gC/m2/day)
    ! 13 = heterotrophic respiration from litter (gC/m2/day)
    ! 14 = heterotrophic respiration from som (gC/m2/day)
    ! 15 = litter decomposition to som (gC/m2/day)
    ! 16 = GPP allocation to autotrophic pool (gC/m2/day)
    ! 17 = NOT IN USE
    ! 18 = NOT IN USE
    ! 19 = NOT IN USE
    ! 20 = NOT IN USE
    ! 21 = harvest yield - storage organ C extracted (gC/m2/day)
    ! 22 = biomass extracted in addition to yield at harvest (gC/m2/day)
    ! 23 = respiration from autotrophic pool allocation (gC/m2/day)
    ! 24 = respiration from labile->NPP translocation (gC/m2/day)
    ! 25 = respiration from NPP->labile translocation (gC/m2/day)
    ! 26 = respiration from foliage remobilisation to labile (gC/m2/day)
    ! 27 = respiration from stem remobilisation to labile (gC/m2/day)
    ! 28 = foliage extracted from harvest (gC/m2/day)
    ! 29 = stem extracted from harvest (gC/m2/day)
    ! 30 = dead standing foliage extracted from harvest (gC/m2/day)
    ! 31 = labile extracted from harvest (gC/m2/day)
    ! 32 = foliage added to litter from harvest (gC/m2/day)
    ! 33 = stem added to litter from harvest (gC/m2/day)
    ! 34 = dead standing foliage added to litter from harvest (gC/m2/day)
    ! 35 = autotrophic pool added to litter from harvest (gC/m2/day)
    ! 36 = labile added to litter from harvest (gC/m2/day)
    ! 37 = roots added to litter from ploughing (gC/m2/day)

    ! PARAMETERS are:
    ! p(1)  = litter decomposition rate (fraction/day)
    ! p(2)  = fraction of GPP allocated to autotrophic C pool (fraction)
    ! p(3)  = maximum development rate coefficient DS (0->1) (day-1)
    ! p(4)  = maximum development rate coefficient DS (1->2) (day-1)
    ! p(5)  = turnover rate of foliage (fraction/day)
    ! p(6)  = turnover rate of stem (fraction/day)
    ! p(7)  = maximum foliage turnover rate due to self-shading (fraction/day)
    ! p(8)  = effective vernalisation days when plant is 50% vernalised (days)
    ! p(9)  = litter turnover rate, temperature adjusted (fraction/day)
    ! p(10) = som turnover rate, temperature adjusted (fraction/day)
    ! p(11) = photosynthetic nitrogen use efficiency (gC/gN/m2/day)
    ! p(12) = sow day (day of year)
    ! p(13) = phenological heat units required for emergence
    ! p(14) = harvest day offset from sow day (days)
    ! p(15) = not used in this model configuration
    ! p(16) = not used in this model configuration
    ! Note: respiratory cost of labile transfer is a hardcoded constant (0.21875)
    ! p(17) = leaf mass per area LMA (gC/m2)
    ! p(18) = initial labile C pool (gC/m2)
    ! p(19) = initial foliar C pool (gC/m2)
    ! p(20) = initial root C pool (gC/m2)
    ! p(21) = initial stem C pool (gC/m2)
    ! p(22) = initial litter C pool (gC/m2)
    ! p(23) = initial som C pool (gC/m2)
    ! p(24) = initial autotrophic pool C (gC/m2)
    ! p(25) = initial storage organ C pool (gC/m2)
    ! p(26) = minimum temperature for development (K)
    ! p(27) = maximum temperature for development (K)
    ! p(28) = optimum temperature for development (K)
    ! p(29) = minimum temperature for vernalisation (K)
    ! p(30) = maximum temperature for vernalisation (K)
    ! p(31) = optimum temperature for vernalisation (K)
    ! p(32) = critical photoperiod for development (hours)
    ! p(33) = photoperiod sensitivity
    ! p(34) = turnover rate of labile C (fraction/day)
    ! p(35) = turnover rate of autotrophic C (fraction/day)
    ! p(36) = canopy nitrogen dilution intercept (gN/m2)
    ! p(37) = canopy nitrogen dilution coefficient

    ! Set some initial states for the io variables
    infi = 0d0 ; FLUXES = 0d0 ; POOLS = 0d0 ; DIAGS = 0d0

    ! load some values
    acm_forcings(7) = pars(11) ! Canopy nitrogen use efficency (gC/gNleaf/day)
    mV%avN = pars(36)       ! foliar N, initial value                         
    mV%latitude_radians = lat * deg_to_rad
    mV%sin_latitude_radians = sin(mV%latitude_radians)
    mV%cos_latitude_radians = cos(mV%latitude_radians)

    ! parameters from file
    mV%decomposition_rate                = pars(1)  ! decomposition rate (day)
    mV%frac_GPP_resp_auto                = pars(2)  ! fraction of GPP allocated to autotrophic carbon pool
    mV%DR_pre                            = pars(3)  ! development rate coefficient DS (0->1)
    mV%DR_post                           = pars(4)  ! development rate coefficient DS (1->2)
    mV%turnover_rate_foliage             = pars(5)  ! turnover_rate of foliage (day)
    mV%turnover_rate_stem                = pars(6)  ! turnover rate of stem (day)
    mV%RDRSHMAX                          = pars(7)  ! maximum rate of foliar turnover due to self shading (day)
    mV%VDh                               = pars(8)  ! effective vernalisation days when plants are 50 % vernalised
    mV%mineralisation_rate_litter        = pars(9)  ! mineralisation rate litter (day)
    mV%mineralisation_rate_soilOrgMatter = pars(10) ! mineralisation rate som (day)
    mV%sow_day                           = nint(mod(pars(12),365.25d0)) ! sow day (doy)
    mV%PHUem                             = pars(13) ! phenological heat units required for emergence
    mV%harvest_day                       = nint(mod(mV%sow_day + pars(14),365.25d0)) ! harvest day (doy)
    mV%plough_day                        = nint(mod(pars(12)-2d0,365.25d0)) ! plough day (doy)
    mV%LCA                               = pars(17) ! leaf mass area (gC.m-2)
    mV%tmin                              = pars(26)-273.15d0 ! min temperature for development
    mV%tmax                              = pars(27)-273.15d0 ! max temperature for development
    mV%topt                              = pars(28)-273.15d0 ! optimum temperature for development
    mV%tmin_v                            = pars(29)-273.15d0 ! min temperature for vernalisation
    mV%tmax_v                            = pars(30)-273.15d0 ! max temperature for vernalisation
    mV%topt_v                            = pars(31)-273.15d0 ! optimim temperature for vernalisation
    mV%doptmin   = mV%topt   - mV%tmin     ! difference between optimum and minimum cardinal temperatures
    mV%dmaxmin   = mV%tmax   - mV%tmin     ! difference between maximum and minimum cardinal temperatures
    mV%doptmin_v = mV%topt_v - mV%tmin_v   ! difference between optimum and minimum vernalisation temperatures
    mV%dmaxmin_v = mV%tmax_v - mV%tmin_v   ! difference between maximum and minimum vernalisation temperatures
    mV%PHCR                              = pars(32) ! critical value of photoperiod for development
    mV%PHSC                              = pars(33) ! photoperiod sensitivity
    mV%turnover_rate_labile              = pars(34) ! turnover rate labile C (day)
    mV%turnover_rate_resp_auto           = pars(35) ! turnover rate of autotrophic carbon for respiration (day)

    ! load stocks in first time step
    mV%stock_labile                      = pars(18) ! labile C
    mV%stock_foliage                     = pars(19) ! foliar C
    mV%stock_roots                       = pars(20) ! root C
    mV%stock_stem                        = pars(21) ! stem / wood C
    mV%stock_litter                      = pars(22) ! litter C
    mV%stock_soilOrgMatter               = pars(23) ! som C
    mV%stock_resp_auto                   = pars(24) ! autotrophic resp pool
    mV%stock_storage_organ               = pars(25) ! storage organ (i.e. desired crop)
    mV%stock_dead_foliage                = 0d0      ! dead but still standing foliage

    ! assigning initial conditions
    POOLS(1,1) = mV%stock_labile
    POOLS(1,2) = mV%stock_foliage
    POOLS(1,3) = mV%stock_roots
    POOLS(1,4) = mV%stock_stem
    POOLS(1,5) = mV%stock_litter
    POOLS(1,6) = mV%stock_soilOrgMatter
    POOLS(1,7) = mV%stock_resp_auto
    POOLS(1,8) = mV%stock_storage_organ
    POOLS(1,9) = mV%stock_dead_foliage

    ! logical switches
    mV%vernal_calcs    = .true.
    mV%ploughed        = .false.
    mV%sown            = .false.
    mV%use_seed_labile = .true.
    mV%emerged         = .false.

    ! pair incoming variables to local module levels
    ! finally set some initial conditions
    mV%yield = 0d0
    mV%DS = -1d0
    mV%DR = 0d0
    fV = 0d0 ; fT = 0d0 ; fP = 0d0
    mV%alloc_to_storage_organ_old = 0d0
    mV%PHU = 0d0 ; mV%VD = 0d0
    mV%BM_EX = 0d0 ; mV%HI = 0d0
    mV%stock_dead_foliage = 0d0
    mV%alloc_to_labile = 0d0
    mV%stmob = 0
    mV%max_raso = 0d0
    mV%raso = 0d0
    mV%RDRDV = 0d0
    mV%max_raso_old = 0d0
    mV%raso_old  = 0d0
    mV%fol_frac_intpol = 0d0
    mV%stem_frac_intpol = 0d0
    mV%root_frac_intpol = 0d0

    ! SHOULD TURN THIS INTO A SUBROUTINE CALL AS COMMON TO BOTH DEFAULT AND CROPS
    if (.not.allocated(mV%deltat_1)) then 
      write(*,*) "Error - arrays not allocated - probably carbon_model() was called without initialize_mv()"
      STOP 1
    endif

    ! load some needed module level values
    mV%lai = POOLS(1,2)/mV%LCA
    mV%mint = met(2,1)  ! minimum temperature (oC)
    mV%maxt = met(3,1)  ! maximum temperature (oC)
    leafT = (mV%maxt*0.75d0) + (mV%mint*0.25d0) ! initial canopy temperature (oC)
    mV%swrad = met(4,1) ! incoming short wave radiation (MJ/m2/day)
    mV%co2 = met(5,1)   ! CO2 (ppm)
    mV%doy = ceiling(met(6,1)-(deltat(1)*0.5d0))   ! Day of year
    meant = (mV%maxt+mV%mint) * 0.5d0   ! mean air temperature (oC) !TODO no such variable in this module

    !
    ! Begin looping through each time step
    !

    do n = start, finish

      !!!!!!!!!!
      ! assign drivers and update some prognostic variables
      !!!!!!!!!!

      ! Incoming drivers
      mV%mint = met(2,n)  ! minimum temperature (oC)
      mV%maxt = met(3,n)  ! maximum temperature (oC)
      leafT = (mV%maxt*0.75d0) + (mV%mint*0.25d0)     ! initial canopy temperature (oC)
      mV%swrad = met(4,n) ! incoming short wave radiation (MJ/m2/day)
      mV%co2 = met(5,n)   ! CO2 (ppm)
      mV%days_per_step = deltat(n) ; mV%days_per_step_1 = mV%deltat_1(n)
      mV%doy = ceiling(met(6,n)-(mV%days_per_step*0.5d0))   ! Day of year
      meant = met(14,n)   ! mean air temperature (oC)

      ! Calculate current days leaf area index
      mV%lai = POOLS(n,2)/mV%LCA
      DIAGS(n,1) = mV%lai

      ! Calculate solar declination for the current time step
      mV%declination = calculate_declination(mV%doy)
      ! calculate daylength in hours and seconds
      call calculate_daylength(mV)
      ! extract timing related values

      ! Note that soil mass balance will be calculated after phenology
      ! adjustments

      ! DS < 0.15 corresponds to the growth stage at beginning of the UK recommended period of      
      ! N fertiliser application for winter wheat (Zodocks growth stage 20) - the early tillering stage (typically mid-march to April)
      if (mV%DS < 0.469d0) then                          
          mV%avN = pars(36)                     
      else! if (DS >= 0.469d0 .and. DS <= 1.293d0) then
          ! NOTE: The slope_n parameter can be included in the MDF optimisation. 
          !       The value for this parameter has also been observed to be around -0.02.
          ! NOTE: Modified to only allow the dilution equation to dilute not enrich N content.
          !       This is to attempt to get around the N-dilultion model increasing N content 
          !       during senescence which is unrealistic. 
          ! NOTE: Dead foliage removed as only the remaining live foliage is photosynthetically active.
          !avN = max(0.1d0,min(avN,(pars(37)*(POOLS(n,2)+POOLS(n,9))) + pars(36)))
          mV%avN = max(0.1d0,min(mV%avN,(pars(37)*POOLS(n,2)) + pars(36)))          
!      else
!          ! Set LNA to 0.1 after anthesis (Zodocks growth stage 75)  
!          Non-applicable as we are explicitly tracking the continued live leaf area and the dead
!          avN = 0.1d0                              
      end if
      DIAGS(n,2) = mV%avN

      ! GPP and direct allocation of GPP can only occur 
      ! if sufficient LAI to prevent numerical error
      if (mV%lai > 1d-5) then 

          ! load next met / lai values for ACM
          acm_forcings(1) = mV%lai      ! LAI
          acm_forcings(2) = mV%maxt ! maximum temperature (oC)
          acm_forcings(3) = mV%mint ! minimum temperature (oC)
          acm_forcings(4) = mV%avN
          acm_forcings(5) = mV%co2 ! CO2 (ppm)
          acm_forcings(6) = mV%swrad ! incoming short wave radiation (MJ/m2/day)

          ! GPP (gC.m-2.day-1)
          FLUXES(n,1) = acm(acm_forcings)
          DIAGS(n,4) = mV%ci / mV%co2
          
      end if
      ! Pass GPP estimate to modul
e variable for use in crop development model
      mV%gpp_acm = FLUXES(n,1) 

      ! calculate weighted air temperature value based on daily minimum, maximum
      ! and means. This minimises the error introduced when scaling between
      ! daily and sub-daily timesteps
      !airt_weighting(1) = abs(met(3,n)-meant) / (met(3,n)-met(2,n))*0.5d0 ! maximum temperature weighting
      !airt_weighting(2) = 0.5d0                                            ! mean temperature
      !airt_weighting(3) = abs(met(2,n)-meant) / (met(3,n)-met(2,n))*0.5d0 ! minimum temperature weighting

      ! Heterotrophic respiration rate (Q10):  doubles with
      ! 10 degree temperature rise resprate from soil file = 0.0693
      !resp_rate = 0d0
      !resp_rate = resp_rate + ((0.5d0 * exp( resp_rate_temp_coeff * met(3,n) )) * airt_weighting(1))
      !resp_rate = resp_rate + ((0.5d0 * exp( resp_rate_temp_coeff * meant   )) * airt_weighting(2))
      !resp_rate = resp_rate + ((0.5d0 * exp( resp_rate_temp_coeff * met(2,n) )) * airt_weighting(3))
      mV%resp_rate = 0.5d0 * exp( resp_rate_temp_coeff * meant )

      ! reallocate day of year to the end of the time step for use in
      ! crop development model
      mV%doy = met(6,n)
      ! determine development stage (DS)
      ! Note that DS must be updated here and not after management_dates. Otherwise, 
      ! there will be problems using DS_time for calculating yield:GPP calculations in MODEL_LIKELIHOOD.f90
      call development_stage(mV%days_per_step) ; DIAGS(n,3, mV) = mV%DS 
      ! determine the carbon partitioning based on development stage
      call carbon_alloc_fractions(CI%DS_shoot,CI%DS_root,CI%fol_frac,CI%stem_frac,CI%root_frac, mV)
      ! begin carbon allocation for crops
      call calc_pools_crops(CI%DS_LRRT,CI%LRRT, mV)
      ! conduct management updates at the end of the day
      call management_dates(CI%stock_seed_labile,mV%days_per_step, mV)

      ! temprate (i.e. temperature modified rate of metabolic activity)
      FLUXES(n,2) = mV%resp_rate
      ! autotrophic respiration (gC.m-2.d-1)
      FLUXES(n,3) = mV%resp_auto + mV%resp_cost_labile_to_npp + &
                    mV%resp_cost_npp_to_labile + mV%resp_cost_foliage_to_labile + &
                    mV%resp_cost_stem_to_labile
      ! leaf production rate (gC.m-2.d-1)
      FLUXES(n,4) = mV%alloc_to_foliage
      ! labile production (gC.m-2.d-1)
      FLUXES(n,5) = mV%alloc_to_labile + mV%foliage_to_labile + mV%stem_to_labile
      ! root production (gC.m-2.d-1)
      FLUXES(n,6) = mV%alloc_to_roots
      ! stem production (gC.m-2.d-1)
      FLUXES(n,7) = mV%alloc_to_stem
      ! labile from NPP (gC.m-2.d-1)
      FLUXES(n,8) = mV%alloc_from_labile
      ! alloc to storage organ (gC.m-2.d-1)
      FLUXES(n,9) = mV%alloc_to_storage_organ
      ! total leaf litter production (gC.m-2.d-1)
      FLUXES(n,10) = mV%litterfall_foliage
      ! total stem litter production (gC.m-2.d-1)
      FLUXES(n,11) = mV%litterfall_stem
      ! total root litter production (gC.m-2.d-1)
      FLUXES(n,12) = mV%litterfall_roots
      ! respiration heterotrophic litter (gC.m-2.d-1)
      FLUXES(n,13) = mV%resp_h_litter
      ! respiration heterotrophic som (gC.m-2.d-1)
      FLUXES(n,14) = mV%resp_h_soilOrgMatter
      ! litter to som (gC.m-2.d-1)
      FLUXES(n,15) = mV%decomposition
      ! alloc to autotrophic pool (gC.m-2.d-1)
      FLUXES(n,16) = mV%alloc_to_resp_auto
      ! harvest yield (gC.m-2.d-1)
      FLUXES(n,21) = mV%yield
      ! C extracted in addition to yield due to harvest activity (gC.m-2.d-1)
      FLUXES(n,22) = mV%BM_EX
      ! Respiration from autotrophic allocation (gC.m-2.d-1)
      FLUXES(n,23) = mV%resp_auto
      ! Respiration from labile to foliage translocation (gC.m-2.d-1)
      FLUXES(n,24) = mV%resp_cost_labile_to_npp
      ! Respiration from foliage to litter translocation (gC.m-2.d-1)
      FLUXES(n,25) = mV%resp_cost_npp_to_labile
      ! Respiration from foliage remobilisation (gC.m-2.d-1)
      FLUXES(n,26) = mV%resp_cost_foliage_to_labile
      ! Respiration from stem remobilisation (gC.m-2.d-1)
      FLUXES(n,27) = mV%resp_cost_stem_to_labile
      ! Foliage extracted from harvest (gC.m-2.d-1)
      FLUXES(n,28) = mV%HARVESTextracted_foliage
      ! Stem extracted from harvest (gC.m-2.d-1)
      FLUXES(n,29) = mV%HARVESTextracted_stem
      ! Dead still standing foliage extracted from harvest (gC.m-2.d-1)
      FLUXES(n,30) = mV%HARVESTextracted_dead_foliage
      ! Labile extracted from harvest (gC.m-2.d-1)
      FLUXES(n,31) = mV%HARVESTextracted_labile
      ! Foliage added to litter due to harvest (gC.m-2.d-1)
      FLUXES(n,32) = mV%HARVESTlitter_foliage
      ! Stem added to litter due to harvest (gC.m-2.d-1)
      FLUXES(n,33) = mV%HARVESTlitter_stem
      ! Dead still standing foliage added to litter due to harvest (gC.m-2.d-1)
      FLUXES(n,34) = mV%HARVESTlitter_dead_foliage
      ! Autotrophic pool added to litter due to harvest (gC.m-2.d-1)
      FLUXES(n,35) = mV%HARVESTlitter_resp_auto
      ! Labile added to litter due to harvest (gC.m-2.d-1)
      FLUXES(n,36) = mV%HARVESTlitter_labile
      ! Roots pool added to litter due to plough (gC.m-2.d-1)
      FLUXES(n,37) = mV%PLOUGHlitter_roots

      ! labile pool
      POOLS(n+1,1) = mV%stock_labile
      ! foliar pool
      POOLS(n+1,2) = mV%stock_foliage
      ! root pool
      POOLS(n+1,3) = mV%stock_roots
      ! wood pool
      POOLS(n+1,4) = mV%stock_stem
      ! litter pool
      POOLS(n+1,5) = mV%stock_litter
      ! som pool
      POOLS(n+1,6) = mV%stock_soilOrgMatter
      ! autotrophic pool
      POOLS(n+1,7) = mV%stock_resp_auto
      ! storage organ pool
      POOLS(n+1,8) = mV%stock_storage_organ
      ! dead but still standing foliage
      POOLS(n+1,9) = mV%stock_dead_foliage

!      do nxp = 1, nopools
!         if (POOLS(n+1,nxp) /= POOLS(n+1,nxp) .or. POOLS(n+1,nxp) < 0d0) then
!             print*,"step",n,"POOL",nxp
!             print*,"met",met(:,n)
!             print*,"POOLS",POOLS(n,:)
!             print*,"FLUXES",FLUXES(n,:)
!             print*,"POOLS+1",POOLS(n+1,:)
!             print*,"labile, foliage",stock_labile, stock_foliage
!             print*,"stem, roots",stock_stem,stock_roots
!             print*,"litter, som",stock_litter,stock_soilOrgMatter
!             print*,"storageOrgan, auto",stock_storage_organ,stock_resp_auto
!             print*,"gpp, nee",gpp_acm,nee_dalec
!             print*,"Ra, Rh_som, Rh_lit",resp_auto,resp_h_soilOrgMatter,resp_h_litter
!             print*,"pars",pars
!             print*,"DR",DR
!             print*,"DR stuff",fT,fV,fP
!             print*,"leaf","stem","rootlitter",litterfall_foliage,litterfall_stem,litterfall_roots
!            print*,"daylength",dayl_hours,"VD",VD,"VDh",VDh
!             print*,"meant",meant
!            print*,"sown",sown,"emerged",emerged
!             print*,"root_frac_intpol",root_frac_intpol
!             print*,"npp_shoot",npp_shoot,"npp",npp
!             print*,"RDR",RDR
!             !stop
!         endif
!     enddo
!
!      do nxp = 1, nofluxes
!
!            if (FLUXES(n,nxp) /= FLUXES(n,nxp) .or. FLUXES(n,nxp) < 0d0) then
!                 print*,"Special: step",n,"FLUXES",nxp
!                 print*,"met",met(:,n)
!                 print*,"POOLS",POOLS(n,:)
!                 print*,"FLUXES",FLUXES(n,:)
!                 print*,"POOLS+1",POOLS(n+1,:)
!                print*,"labile, foliage",stock_labile, stock_foliage
!                print*,"stem, roots",stock_stem,stock_roots
!                 print*,"litter, som",stock_litter,stock_soilOrgMatter
!                 print*,"storageOrgan, auto",stock_storage_organ,stock_resp_auto
!                 print*,"gpp, nee",gpp_acm,nee_dalec
!                 print*,"Ra, Rh_som, Rh_lit",resp_auto,resp_h_soilOrgMatter,resp_h_litter
!                 print*,"pars",pars
!                 print*,"DR",DR
!                 print*,"DR stuff",fT,fV,fP
!                 print*,"leaf","stem","rootlitter",litterfall_foliage,litterfall_stem,litterfall_roots
!                 print*,"daylength",dayl_hours,"VD",VD,"VDh",VDh
!                print*,"meant",meant
!                 print*,"sown",sown,"emerged",emerged
!                 print*,"root_frac_intpol",root_frac_intpol
!                 print*,"npp_shoot",npp_shoot,"npp",npp
!                 print*,"RDR",RDR
!                 !stop
!            end if
!      enddo
!
!      if (stock_labile < 0d0 .or. stock_foliage < 0d0 .or. stock_stem < 0d0 .or. &
!          stock_roots < 0d0 .or. stock_litter < 0d0 .or. stock_soilOrgMatter < 0d0 .or. &
!          stock_storage_organ < 0d0 .or. stock_resp_auto < 0d0 .or. &
!          stock_labile /= stock_labile .or. stock_foliage /= stock_foliage .or. &
!          stock_stem /= stock_stem .or. &
!          stock_roots /= stock_roots .or. stock_litter /= stock_litter .or. &
!          stock_soilOrgMatter /= stock_soilOrgMatter .or. &
!          stock_storage_organ /= stock_storage_organ .or. &
!          stock_resp_auto /= stock_resp_auto .or.  &
!          gpp_acm < 0d0 .or. gpp_acm /= gpp_acm .or. resp_rate < 0d0 .or. &
!          resp_rate /= resp_rate .or. decomposition < 0d0 .or. alloc_from_labile < 0d0 .or. &
!          resp_cost_labile_to_npp < 0d0 .or. alloc_to_foliage < 0d0 .or. &
!          alloc_to_stem < 0d0 .or. alloc_to_roots < 0d0 .or. &
!          alloc_from_labile < 0d0 .or. resp_cost_labile_to_npp < 0d0) then
!          print*,"stocks less than zero or NaN", n
!          print*,stock_labile, stock_foliage
!          print*,stock_stem,stock_roots
!          print*,stock_litter,stock_soilOrgMatter
!          print*,stock_storage_organ,stock_resp_auto
!          print*,gpp_acm,nee_dalec
!          print*,resp_auto,resp_h_soilOrgMatter,resp_h_litter
!          print*,"pars",pars
!          print*,"fluxes",fluxes(n,:)
!          print*,"DR",DR
!          print*,"DR stuff",fT,fV,fP
!          print*,"leaf","stem","rootlitter",litterfall_foliage,litterfall_stem,litterfall_roots
!          print*,"daylength",dayl_hours,"VD",VD,"VDh",VDh
!          print*,"meant",meant
!          print*,"sown",sown,"emerged",emerged
!          print*,"root_frac_intpol",root_frac_intpol
!          print*,"npp_shoot",npp_shoot,"npp",npp
!          print*,"RDR",RDR
!          !stop
!      endif

    end do ! no days loop

  end subroutine CARBON_MODEL
  !
  !------------------------------------------------------------------
  !
  double precision function acm(acm_forcings, mV)

    ! the Aggregated Canopy Model, is a Gross Primary Productivity (i.e.
    ! Photosyntheis) emulator which operates at a daily time step. ACM can be
    ! paramaterised to provide reasonable results for most ecosystems.

    implicit none

      type(model_working_variables) :: mV

    ! declare input variables
    double precision, intent(in) :: acm_forcings(7) ! acm input requirements

    ! declare local variables
    double precision, parameter :: dayl_coef = 0.0156935d0, &
                              co2_comp_point = 4.22273d0,   &
                                co2_half_sat = 208.868d0,   &
                                  dayl_const = 0.0453194d0, &
                         hydraulic_temp_coef = 0.37836d0,   &
                                    lai_coef = 7.19298d0,   &
                               temp_exponent = 0.011136d0,  &
                                   lai_const = 2.1001d0,    &
                          hydraulic_exponent = 0.789798d0,  &
                                     deltaWP = -2d0,        & ! leafWP-soilWP
                                        Rtot = 1d0            ! Total hydraulic resistance
    double precision :: gc,pn,pd,pp,qq,e0,cps,nit,NUE &
                       ,trange,mint,maxt,radiation,co2,lai,doy

    ! initial values
    gc = 0d0 ; pp = 0d0 ; qq = 0d0 ; mV%ci = 0d0 ; e0 = 0d0 ; cps = 0d0 

    ! load driver values to correct local vars
    lai  = acm_forcings(1) ! leaf area index m2/m2
    maxt = acm_forcings(2) ! mean of daily maximum temperature C
    mint = acm_forcings(3) ! mean of daily minimum temperature C
    nit  = acm_forcings(4) ! mean foliar nitrogen gN/m2leaf
    co2  = acm_forcings(5) ! ppm
    radiation = acm_forcings(6) ! MJ/m2/day
    NUE = acm_forcings(7) ! Nitrogen use efficiency

    ! determine temperature range 
    trange = 0.5*(maxt-mint)
    ! daily canopy conductance, of CO2 or H2O? 
    gc = abs(deltaWP)**(hydraulic_exponent)/((hydraulic_temp_coef*Rtot+trange))
    ! maximum rate of temperature and nitrogen (canopy efficiency) limited photosynthesis (gC.m-2.day-1)
    pn = lai*nit*NUE*exp(temp_exponent*maxt)
    ! pp and qq represent limitation by diffusion and metabolites respecitively
    pp = pn/gc 
    qq = co2_comp_point-co2_half_sat
    ! calculate internal CO2 concentration (ppm)
    mV%ci = 0.5*(co2+qq-pp+sqrt(((co2+qq-pp)*(co2+qq-pp))-4d0*(co2*qq-pp*co2_comp_point)))
    ! limit maximum quantium efficiency by leaf area, hyperbola
    e0 = lai_coef*(lai*lai)/((lai*lai)+lai_const)

    ! calculate CO2 limited rate of photosynthesis
    pd = gc*(co2-mV%ci)
    ! calculate combined light and CO2 limited photosynthesis
    cps = e0*radiation*pd/(e0*radiation+pd)
    ! correct for day length variation
    acm = cps*(dayl_coef*mV%dayl_hours+dayl_const)

    return

  end function acm
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_cosine_solar_zenith_angle (mV)

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
  !-----------------------------------------------------------------
  !
  subroutine calculate_daylength (mV)

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
    !dayl_seconds = dayl_hours * seconds_per_hour

    ! return to user
    return

  end subroutine calculate_daylength
  !
  !------------------------------------------------------------------
  !

  !
  !-----------------------------------------------------------------
  !
  !
  !------------------------------------------------------------------
  !
  subroutine calc_pools_crops(DS_LRRT,LRRT, mV)

    ! Allocated GPP to NPP and various carbon pools. Based  !
    ! this on physiological responses to temperature        !
    ! vernalisation, and photoperiod. Note: this code has   !
    ! been reformulated to a daily time step as default     !
    ! rather SPAcrop which operated on a hourly time step.  !

    implicit none

      type(model_working_variables) :: mV

    ! arguments
    double precision, dimension(:), intent(inout) :: DS_LRRT, & ! Development stage corresponding to loss rate of roots
                                                        LRRT    ! loss rate of roots array

    ! local variables
    double precision :: decomp_efficency

    ! turnover rate of fine roots is now equal to the
    ! loss rate of roots as a function of development stage
    ! (Penning de Vries, 1989).
    mV%turnover_rate_roots = interpolate( mV%DS , DS_LRRT , LRRT , 5 )

    ! if sown turn on labile / seed turnover for growth
    if ( mV%sown ) then
        ! turnover on
        mV%turnover_labile_switch = 1
    else
        ! turnover off
        mV%turnover_labile_switch = 0
    endif

    ! Initialise / reset
    mV%yield = 0d0 ; mV%BM_EX = 0d0
    mV%HARVESTextracted_foliage = 0d0
    mV%HARVESTextracted_stem = 0d0
    mV%HARVESTextracted_dead_foliage = 0d0
    mV%HARVESTextracted_labile = 0d0
    mV%HARVESTlitter_labile = 0d0
    mV%HARVESTlitter_foliage = 0d0
    mV%HARVESTlitter_stem = 0d0
    mV%HARVESTlitter_dead_foliage = 0d0
    mV%HARVESTlitter_resp_auto = 0d0
    mV%PLOUGHlitter_roots = 0d0

    ! Determine total allocation of labile C to NPP or Ra (gC.m-2.d-1)
    mV%alloc_from_labile = mV%turnover_rate_labile * mV%resp_rate * dble(mV%turnover_labile_switch)
    mV%alloc_from_labile = mV%stock_labile * (1d0-(1d0-mV%alloc_from_labile) ** mV%days_per_step) * mV%days_per_step_1
    ! Determine the respiratory cost of C transfer from labile pool to NPP (gC.m-2.d-1)
    mV%resp_cost_labile_to_npp = mV%alloc_from_labile * resp_cost_labile_trans
    ! Determine the remaining allocation of labile C to NPP (gC.m-2.d-1)
    mV%alloc_from_labile = mV%alloc_from_labile - mV%resp_cost_labile_to_npp

    ! Estimate allocation (gC/m2/day) of GPP to autotrophic respiration pool
    mV%alloc_to_resp_auto = mV%gpp_acm * mV%frac_GPP_resp_auto
    ! NPP as a fraction of GPP (1-.32=.68 or 68%) + allocation from labile
    mV%npp = mV%gpp_acm + mV%alloc_from_labile - mV%alloc_to_resp_auto

    ! Dertermine partitioning of NPP to biomass pools
    mV%alloc_to_roots    = mV%root_frac_intpol * mV%npp
    ! Calculate how much NPP is left after allocation to roots
    mV%npp_shoot         = mV%npp - mV%alloc_to_roots
    ! Partition the remaining NPP to foliage and stem (gC/m2/day)
    mV%alloc_to_foliage  = mV%fol_frac_intpol  * mV%npp_shoot
    mV%alloc_to_stem     = mV%stem_frac_intpol * mV%npp_shoot
    ! Remaining NPP is then allocated to storage organ
    mV%alloc_to_storage_organ = max(0d0,mV%npp_shoot - mV%alloc_to_foliage - mV%alloc_to_stem)

    ! Assuming allocation to storage organ is > 0 ensure flux is limited by
    ! maximum growth rate potential, i.e. growth potential increases with size
    ! existing yield
    mV%alloc_to_labile = 0d0 ; mV%resp_cost_npp_to_labile = 0d0
    if ( mV%alloc_to_storage_organ > 0d0 ) then
        mV%gso_max  = ( mV%stock_storage_organ + 0.5d0 ) * rel_gso_max
        ! Restrict allocation to storage organ based on available C and potential
        ! Where does any excess go?
        mV%alloc_to_storage_organ = min( mV%alloc_to_storage_organ , mV%gso_max )
        if ( mV%sown .and. mV%emerged) then
            ! Assuming we shown and emerged, the remainder is now allocated to labile,
            ! less the cost of transfer
            mV%alloc_to_labile = ( mV%npp_shoot - mV%alloc_to_foliage - mV%alloc_to_stem - mV%alloc_to_storage_organ )
            ! Having worked out the total available C, update based on cost of respiratory transfer
            mV%resp_cost_npp_to_labile =  mV%alloc_to_labile * resp_cost_labile_trans
            mV%alloc_to_labile = mV%alloc_to_labile - mV%resp_cost_npp_to_labile
        endif
    endif

    ! set switches to (de)activate leaf, root and stem remobliization
    mV%raso_old = mV%raso
    ! Running Average of growth rate of Storage Organ..
    mV%raso = ( mV%alloc_to_storage_organ + mV%alloc_to_storage_organ_old ) * 0.5d0
    mV%max_raso_old = mV%max_raso
    mV%max_raso = max( mV%raso , mV%max_raso_old )
    ! Stem remobilisation triggered once running average of storage organ growth declines
    ! Second part prevents premature remobilisation
    if ( ( mV%raso < mV%raso_old ) .and. &
         ( mV%alloc_to_storage_organ > ( mV%alloc_to_storage_organ_old + 0.5d0 ) ) ) then
        mV%stmob = 1
    else
        mV%stmob = 0
    endif

    ! Code for calculating relative death rate of leaves (RDR) as a
    !  function of shading (RDRSH) or developmental stage (RDRT).

    ! GT 0 if LAI GT 4; 0. < RDRSH < RDRSHMAX (usually ~0.03)
    mV%RDRSH = min( mV%RDRSHMAX , max( 0d0 , mV%RDRSHMAX * ( mV%lai - LAICR ) / LAICR ) )
    if ( mV%DS < 1d0 ) then
       mV%RDRDV = 0d0
    else
       ! RDRDV dependant on DR and DS, values range typically between 0.02 <
       ! RDRDV < 0.25
!print*,"!! What RDRDV to use? !!"
!!$      RDRDV = DR /( max( 0.1 , 2. - DS ) )
!!$      RDRDV = RDRDV / 24. ! to get hourly senescence rate
       ! TLS: think the source of this equation needs to be found
       mV%RDRDV = mV%turnover_rate_foliage * min(1d0, ( 1d0 / ( ( max( 2d0 - mV%DS , 0.1d0 ) ) * 8d0 ) ) ** 2)
    ENDIF

    ! Relative leaf death rate is the maximum value of the arguments RDRSH and
    ! RDRDV
    mV%RDR = max( mV%RDRSH , mV%RDRDV )

    ! remobilization of foliar C and allocation to dead leaves pool (gC.m-2.d-1)
    mV%litterfall_foliage = mV%stock_foliage * (1d0-(1d0-mV%RDR) ** mV%days_per_step) * mV%days_per_step_1
    ! remobilization of stem C (gC.m-2.d-1)
    mV%litterfall_stem    = mV%stock_stem    * (1d0-(1d0-(mV%DR*mV%turnover_rate_stem*dble(mV%stmob))) ** mV%days_per_step) * mV%days_per_step_1
    ! remobilization of fine roots C (gC.m-2.d-1)
    mV%litterfall_roots   = mV%stock_roots   * (1d0-(1d0-mV%turnover_rate_roots) ** mV%days_per_step) * mV%days_per_step_1

    ! remobilized C to NPP (from both leaves and stems) (gC.m-2.d-1)
    ! Assume that half the foliage litter C is available for remobilisation
    mV%foliage_to_labile   = ( mV%litterfall_foliage * 0.5d0) * ( 1d0 - resp_cost_labile_trans )
    mV%stem_to_labile      = ( mV%litterfall_stem ) * ( 1d0 - resp_cost_labile_trans )
    ! respiratory cost of C transfer (conversion from starch to photosynthates) (gC.m-2.d-1)
    mV%resp_cost_foliage_to_labile   = ( mV%litterfall_foliage * 0.5d0) * resp_cost_labile_trans
    mV%resp_cost_stem_to_labile      = mV%litterfall_stem * resp_cost_labile_trans

    ! litter decomposition to soil organic matter
    mV%decomposition = mV%stock_litter &
                  * (1d0-(1d0-(mV%decomposition_rate*mV%resp_rate)) ** mV%days_per_step) &
                  * mV%days_per_step_1

    ! heterotrophic respiration component 1: mineralisation of litter C pool (gC.m-2.d-1)
    mV%resp_h_litter = mV%stock_litter &
                  * (1d0-(1d0-(mV%mineralisation_rate_litter*mV%resp_rate)) ** mV%days_per_step) &
                  * mV%days_per_step_1
    ! heterotrophic respiration component 2:  mineralisation of organic matter C pool (gC.m-2.d-1)
    mV%resp_h_soilOrgMatter = mV%stock_soilOrgMatter &
                         * (1d0-(1d0-(mV%mineralisation_rate_soilOrgMatter * mV%resp_rate)) ** mV%days_per_step) &
                         * mV%days_per_step_1
                  
    ! Autotrophic respiration allocated pool is a special case, where the total
    ! available for the time step is accumulated first then losses are determined.
    ! Add photosynthate allocated to autotrophic respiration (gC.m-2.d-1)
    mV%stock_resp_auto = mV%stock_resp_auto + (mV%alloc_to_resp_auto * mV%days_per_step)
    ! Determine autotrophic costs as loss term. This should probably have a temperature sensitivity
    ! or be removed to work as applied in other DALEC models? (gC.m-2.d-1)
    ! NOTE: that autotrophic respiration from other sources added below
    mV%resp_auto = mV%stock_resp_auto * (1d0-(1d0-mV%turnover_rate_resp_auto) ** mV%days_per_step) * mV%days_per_step_1
    ! Update stock
    mV%stock_resp_auto = max(0d0, mV%stock_resp_auto - (mV%resp_auto * mV%days_per_step))

    ! nee (gC.m-2.d-1)
    mV%nee_dalec = (mV%resp_auto + mV%resp_cost_labile_to_npp + mV%resp_cost_npp_to_labile + &
                 mV%resp_cost_foliage_to_labile + mV%resp_cost_stem_to_labile + &
                 mV%resp_h_litter + mV%resp_h_soilOrgMatter) - mV%gpp_acm

    ! Recalculate Physical Carbon Pools...
    ! TLS: Note that while max(0d0,... conditions have been applied to all pools, I hypothese that the critical
    ! restriction is applied to stock_resp_auto
    mV%stock_foliage       = max(0d0, mV%stock_foliage       + ((mV%alloc_to_foliage - mV%litterfall_foliage) * mV%days_per_step))
    mV%stock_stem          = max(0d0, mV%stock_stem          + ((mV%alloc_to_stem - mV%litterfall_stem) * mV%days_per_step))
    mV%stock_storage_organ = max(0d0, mV%stock_storage_organ + ((mV%alloc_to_storage_organ) * mV%days_per_step))
    mV%stock_roots         = max(0d0, mV%stock_roots         + ((mV%alloc_to_roots   - mV%litterfall_roots) * mV%days_per_step))
    mV%stock_litter        = max(0d0, mV%stock_litter        + ((mV%litterfall_roots - mV%resp_h_litter - mV%decomposition) * mV%days_per_step))
    mV%stock_soilOrgMatter = max(0d0, mV%stock_soilOrgMatter + ((mV%decomposition    - mV%resp_h_soilOrgMatter) * mV%days_per_step))
    mV%stock_dead_foliage  = max(0d0, mV%stock_dead_foliage  + ((mV%litterfall_foliage * 0.5d0) * mV%days_per_step))
    mV%stock_labile        = max(0d0, mV%stock_labile        + ((mV%alloc_to_labile  - mV%alloc_from_labile - &
                                                           mV%resp_cost_labile_to_npp + mV%foliage_to_labile + &
                                                           mV%stem_to_labile) * mV%days_per_step) )

    ! When GPP is higher than seed C content, remaining seed carbon enters litter
    ! C pool, as seedlings do not fully exhaust their seed (P. de Vries p 48)
    ! Moved to the end of the mass balance, to maintain balance with fluxes which
    ! are already concurrently occuring. In this codes previous position further
    ! up the code, resulted in C being both dumped to litter and the same carbon
    ! being allocated to respiration and NPP.
    if ( ( mV%gpp_acm > mV%alloc_from_labile ) .and. ( mV%use_seed_labile ) ) then
        mV%stock_litter = mV%stock_litter + mV%stock_labile
        mV%stock_labile = 0d0
        mV%use_seed_labile = .false.
    endif

  end subroutine calc_pools_crops
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine carbon_alloc_fractions(DS_shoot,DS_root,fol_frac,stem_frac,root_frac, mV)

    ! Determines carbon allocation fractions as a function !
    ! of developmental stage (DS).  Allocation fractions   !
    ! are from tables published in Penning de Vries (1989) !

    implicit none

      type(model_working_variables) :: mV

    double precision, dimension(:), intent(in) ::   DS_shoot, & !
                                                        DS_root, & !
                                                       fol_frac, & !
                                                      stem_frac, & !
                                                      root_frac    !

    ! local variables..
    double precision,dimension(:),allocatable :: frac_shoot, frac_root

    if ( mV%sown ) then ! after sowing

       ! loop over three crop "organs": 1) foliage 2) stems 3) root
       ! not necessary for storage organs, as all remaining C is allocated to
       ! these

       ! use different input for foliage and stem fractions, as they are
       ! relative to the total shoot (or aboveground) allocation,
       ! root is relative to total plant (above- and belowground) allocation.

       ! leaf development stages and corresponding fractions..

       ! interpolate between PdV allocation values with reference to
       ! developmental stage (DS)..
       mV%fol_frac_intpol = min(1d0,max(0d0,interpolate( mV%DS , DS_shoot , fol_frac , size(DS_shoot) )))
       ! stem DS and fracs..
       mV%stem_frac_intpol = min(1d0,max(0d0,interpolate( mV%DS , DS_shoot , stem_frac , size(DS_shoot) )))
       ! root DS and fracs..
       mV%root_frac_intpol = min(1d0,max(0d0,interpolate( mV%DS , DS_root , root_frac , size(DS_root) )))

    endif ! after crop has been sown

  end subroutine carbon_alloc_fractions
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine development_stage(days_in_step, mV)

    ! Based on modified Wang & Engel model (Streck et al., 2003), !
    ! but with only 2 sub-phases, vegetative and reproductive     !
    ! (i.e. only two different DRmax).   O. Sus, May 2010.        !

    implicit none

      type(model_working_variables) :: mV

    ! agruments
    double precision :: days_in_step

    ! local variables..
    double precision ::  dttmin,  & ! Difference between daily average and minimum temperatures
                         dttmin_v   ! Difference between daily average and minimum vernalization temperatures

    dttmin   = meant - mV%tmin    ! difference between daily average and minimum cardinal temperatures
    dttmin_v = meant - mV%tmin_v  ! difference between daily average and minimum vernalization temperatures

    ! Calculation of developmental function values: vernalization (fV),
    ! temperature (fT) and
    ! photoperiod (fP) these values are multiplicative factors of DRmax (maximum
    ! developmental
    ! rate), each ranging between 0 (no development) and 1 (unrestricted
    ! development).

    ! Summation of vernalization days (VD), not before sowing and only if
    ! average temperature is within min and max cardinal temperatures..
    if ( ( meant > mV%tmin_v ) .and. ( meant < mV%tmax_v ) .and. mV%sown ) then
        fV = vernalization( mV%doptmin_v , mV%dmaxmin_v , dttmin_v , days_in_step )
    endif

    ! Only calculate temperature coefficient if meant lies within (tmin,tmax)
    ! range. NOTE: (doptmin+1d0) < dmaxmin added to allow for EDC search period when "not
    ! allowed" parameter sets will be tried anyway
    if ( meant > mV%tmin .and. meant < mV%tmax .and. (mV%doptmin+1d0) < mV%dmaxmin ) then
        fT = temperature_impact( mV%doptmin , mV%dmaxmin , dttmin )
    else
        fT = 0d0
    endif

    ! calculation of photoperiod coefficient
    fP = photoperiod_impact( mV%PHCR , mV%PHSC )

    if ( mV%emerged .and. ( mV%DS < 2d0 ) ) then   ! sum up daily DR values between emergence and maturity (DS=2)

       if ( mV%DS < 1d0 ) then  ! in the vegetative phase (before flowering):

          mV%DR = mV%DR_pre * fT * fP   ! DR is affected by temperature, photoperiod...
          if ( mV%vernal_calcs ) mV%DR = mV%DR * fV ! ...and vernalization (for winter cereals)
          mV%DS = mV%DS + (mV%DR * days_in_step)    ! developmental stage (DS), calculated as the sum of daily developmental rates

       else    ! in the reproductive phase (after flowering):

          mV%DR = mV%DR_post * fT   ! DR is affected only by temperature
          mV%DS = mV%DS + (mV%DR * days_in_step)

       endif ! vegetative or reproductive phase

    endif ! emerged or not

  end subroutine development_stage
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine management_dates (stock_seed_labile,days_in_step, mV)

    ! This routine should be called at the end of each day of a crops  !
    ! simulation.  It checks whether we should plough/sow/harvest, and !
    ! during the growing establishes when the crop will emerge after   !
    ! sowing, based on heat accumulation (Phenological Heat Units).    !

    implicit none

      type(model_working_variables) :: mV

    ! arguments
    double precision, intent(in) :: stock_seed_labile,days_in_step
    ! local variables
    double precision :: tmp
    logical :: plough_sanity,sow_sanity,harvest_sanity

    ! reset
    plough_sanity = .false. ; sow_sanity = .false. ; harvest_sanity = .false.
    !plough_sanity = .true. ; sow_sanity = .true. ; harvest_sanity = .true.

    ! spring crop
    if (mV%sow_day < mV%harvest_day .and. nint(mV%doy) < mV%harvest_day) sow_sanity = .true.
    if (mV%plough_day < mV%harvest_day .and. nint(mV%doy) < mV%harvest_day) plough_sanity = .true.
    if (mV%harvest_day > mV%sow_day) harvest_sanity = .true.
    ! winter crops
    if (mV%sow_day > mV%harvest_day) sow_sanity = .true.
    if (mV%plough_day > mV%harvest_day) plough_sanity = .true.
    if (mV%harvest_day < mV%plough_day .and. nint(mV%doy) < mV%plough_day) harvest_sanity = .true.

    if ( .not. mV%sown ) then

        ! fresh field...

        if ( plough_sanity .and. .not.mV%ploughed .and. nint(mV%doy) >= mV%plough_day ) then
            ! the field needs ploughing..
            call plough(mV)

        else if ( sow_sanity .and. nint(mV%doy) >= mV%sow_day ) then

            ! ensure that the field has indeed been ploughed
            if (.not.mV%ploughed) call plough(mV)
            ! the field needs sowing..
            mV%sown = .true.

            ! this switch controls whether the labile carbon within the seed is used
            ! for growth
            mV%use_seed_labile = .true.
            mV%stock_labile = stock_seed_labile

        endif ! plough or sow?

    else

        ! crop in field..

        ! calculate when crop emerges..
        if ( .not. mV%emerged ) then

            ! estimate emergence date based on the accumulated phenological heat
            ! units (PHU)
            ! where PHU is the (positive) heat over tmin..
            tmp = max( meant - mV%tmin , 0d0 )*days_in_step
            mV%PHU = mV%PHU + tmp

            ! set the development stage and emergence..
            if ( mV%PHU >= mV%PHUem ) then
                mV%emerged = .true.
                mV%DS = 0d0
            else
                mV%emerged = .false.
                mV%DS = -1d0
            endif

        endif ! emerged or not

        ! note that in this case harvest day has been fixed relative to the sow
        ! day
        if ( harvest_sanity .and. nint(mV%doy) >= mV%harvest_day) then
            ! the field needs harvesting..
            call harvest(mV)
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
  subroutine harvest (mV)

    implicit none

      type(model_working_variables) :: mV

    ! Declare some local variables
    double precision :: tmp, Ctotal

    ! shoot biomass..
    mV%Cshoot = mV%stock_foliage + mV%stock_stem + mV%stock_storage_organ + mV%stock_labile
    ! Total biomass not including the storage organ
    Ctotal = mV%stock_foliage + mV%stock_stem + mV%stock_roots

    ! Can only work if there is biomass in existance
    if (Ctotal > 0d0) then

        ! determine harvest index..
        mV%HI = mV%stock_storage_organ / mV%Cshoot

        ! the stuff we actually want from the harvest...
        mV%yield = mV%stock_storage_organ * mV%days_per_step_1

        ! How much of each pool is extracted during harvest
        mV%HARVESTextracted_foliage      = mV%stock_foliage * ( 1d0 - lv_res ) * mV%days_per_step_1
        mV%HARVESTextracted_stem         = mV%stock_stem * ( 1d0 - st_res ) * mV%days_per_step_1
        mV%HARVESTextracted_dead_foliage = mV%stock_dead_foliage * ( 1d0 - lv_res ) * mV%days_per_step_1
        ! How much of each pool remains as litter after harvest
        mV%HARVESTlitter_foliage      = mV%stock_foliage * lv_res * mV%days_per_step_1
        mV%HARVESTlitter_stem         = mV%stock_stem * st_res * mV%days_per_step_1
        mV%HARVESTlitter_dead_foliage = mV%stock_dead_foliage * lv_res * mV%days_per_step_1
        mV%HARVESTlitter_resp_auto    = mV%stock_resp_auto * mV%days_per_step_1

        ! Labile is a special case due to being distributed within various tissues.
        ! NOTE that extracted is calculated then the litter component is estimates as residual.
        ! The time scale adjustment is applied last, rather than inline (as above)
        mV%HARVESTextracted_labile = (mV%stock_labile * (mV%stock_foliage / Ctotal) *  (1d0 - lv_res )) & 
                                + (mV%stock_labile * (mV%stock_stem / Ctotal) *  (1d0 - st_res )) 
        mV%HARVESTlitter_labile    =  mV%stock_labile - mV%HARVESTextracted_labile
        mV%HARVESTextracted_labile = mV%HARVESTextracted_labile * mV%days_per_step_1
        mV%HARVESTlitter_labile    = mV%HARVESTlitter_labile * mV%days_per_step_1

        ! the biomass that is harvested in addition to the storage-organ..
        mV%BM_EX  = mV%HARVESTextracted_foliage      &
               + mV%HARVESTextracted_stem         &
               + mV%HARVESTextracted_dead_foliage &
               + mV%HARVESTextracted_labile

        ! what's left (will fall to the ground as litter)..
        mV%stock_litter  = mV%stock_litter +               &
                       (mV%HARVESTlitter_foliage +      &
                        mV%HARVESTlitter_stem +         &
                        mV%HARVESTlitter_dead_foliage + &
                        mV%HARVESTlitter_resp_auto +    & 
                        mV%HARVESTlitter_labile) * mV%days_per_step

    end if ! Ctotal > 0

    ! empty the plant stocks..
    mV%stock_storage_organ = 0d0
    mV%stock_foliage       = 0d0
    mV%stock_stem          = 0d0
    mV%stock_dead_foliage  = 0d0
    mV%stock_labile        = 0d0
    mV%stock_resp_auto     = 0d0

    ! roots stay in ground and slowly decompose (until/unless the field is
    ! ploughed)

    ! reset logical variables..
    mV%sown    = .false.
    mV%emerged = .false.
    mV%ploughed = .false.
    mV%DS = -1d0 ; mV%fV = 0d0 ; mV%fT = 0d0 ; mV%fP = 0d0 ; mV%VD = 0d0

  end subroutine harvest
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  double precision function photoperiod_impact( PH_crit , PH_sens , mV)

    ! Function to determine the coefficient for !
    ! photoperiod impact on developmental rate. !
    ! From Streck et al., 2003                  !

    implicit none

      type(model_working_variables) :: mV

    ! arguments..
    double precision,intent(in) :: PH_crit, & ! critical photoperiod below which no development occurs
                                   PH_sens    ! photoperiod sensitivity

    photoperiod_impact = max(0d0, 1d0 - exp ( - PH_Sens * ( mV%dayl_hours - PH_crit ) ))

  end function photoperiod_impact
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine plough (mV)

    ! this s/r will reset various carbon pools, to mimic the effect of the
    ! farmer ploughing. !

    implicit none

      type(model_working_variables) :: mV

    ! Track roots transition to litter here, as all others have been reset already
    mV%PLOUGHlitter_roots = mV%stock_roots * mV%days_per_step_1

    ! Move all plant stocks into the litter pool.
    ! ( many of these should already be empty after the harvest, )
    ! ( e.g. the stocks for labile, foliage, storage-organ stem. )
    mV%stock_litter        = mV%stock_litter + mV%stock_dead_foliage &
                          + mV%stock_foliage + mV%stock_labile    &
                           + mV%stock_roots + mV%stock_stem       &
                            + mV%stock_storage_organ

    mV%stock_dead_foliage  = 0d0
    mV%stock_foliage       = 0d0
    mV%stock_labile        = 0d0
    mV%stock_roots         = 0d0
    mV%stock_stem          = 0d0
    mV%stock_storage_organ = 0d0

    ! Reset the development stage & phenological heat units..
    mV%ploughed = .true. ; mV%DS = -1d0 ; mV%PHU = 0d0
    mV%max_raso = 0d0 ; mV%raso = 0d0 ; mV%max_raso_old = 0d0 ; mV%raso_old = 0d0
    mV%alloc_to_storage_organ_old = 0d0

  end subroutine plough
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

    arrhenious = a * exp( b * (t - 25d0) / (t + freeze) )

  end function arrhenious
  !
  !----------------------------------------------------------------------
  !
  double precision function calculate_declination(doy)

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
  !----------------------------------------------------------------------
  !
  double precision function opt_max_scaling( max_val, min_val , optimum , kurtosis , current )

    ! Estimates a 0-1 scaling based on a skewed guassian distribution with a
    ! given optimum, maximum and kurtosis. Minimum is assumed to be at infinity
    ! (or near enough)

    implicit none

    ! arguments..
    double precision, intent(in) :: max_val, min_val, optimum, kurtosis, current

    ! Code with explicit min bound
    opt_max_scaling = exp( kurtosis * log((max_val-current)/(max_val-optimum)) * (max_val-optimum) ) &
                    * exp( kurtosis * log((current-min_val)/(optimum-min_val)) * (optimum-min_val) )
    ! Sanity check, allows for overlapping parameter ranges
    if (opt_max_scaling /= opt_max_scaling) opt_max_scaling = 0d0

  end function opt_max_scaling
  !
  !------------------------------------------------------------------
  !
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
    !linear_model_gradient = ( (dble(interval)*sum_product_xy) - (sum_x*sum_y) ) &
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
  double precision function vernalization( doptmin_v , dmaxmin_v , dttmin_v , days_in_step , mV)

    ! Function to determine the coefficent for vernalization !
    ! impact on developmental rate. See Streck et al., 2003. !

    implicit none

      type(model_working_variables) :: mV

    ! arguments..
    double precision,intent(in) :: dmaxmin_v , doptmin_v , dttmin_v & ! temperature differences
                                  ,days_in_step

    ! local variables..
    double precision :: a , dnr , fvn , nmr

    a   = log( 2.d0 ) / ( log( ( dmaxmin_v ) / doptmin_v ) )
    nmr = 2.d0 * ( ( dttmin_v ) ** a ) * ( doptmin_v ** a ) - ( ( dttmin_v ) ** (2.d0 * a ) )
    dnr = doptmin_v ** ( 2.d0 * a )
    fvn = nmr / dnr

    mV%VD = mV%VD + (fvn*days_in_step)

    ! final output value..
    vernalization = max( 0d0 , min( 1d0 , ( mV%VD ** 5 ) / ( ( mV%VDh ** 5 ) + (mV%VD ** 5 ) ) ) )

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
!--------------------------------------------------------------------
!
end module CARBON_MODEL_MOD
