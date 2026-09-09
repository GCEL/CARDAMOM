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
! along with this program.  If not, see < https://www.gnu.org/licenses/>.

!!!!!!!!!!!! File specific description !!!!!!!!!!
! Describes module variables containing critical declarable types to support passing
! information around the CARDAMOM software.
!
! This code is based on the original C verion of the University of Edinburgh
! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
! All code translation into Fortran, integration into the University of
! Edinburgh CARDAMOM code and subsequent modifications by:
! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
! J. F. Exbrayat (University of Edinburgh)
! D. T. Milodowski (d.t.milodowski@ed.ac.uk, University of Edinburgh)
! See function/subroutine specific comments for exceptions and contributors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module cardamom_structures

   implicit none (type, external)

   private

   public:: data_type, DATAin, DATAin_original, set_datain, set_datain_original, emulator_parameters, emulator_pars

   type DATA_type

      ! drivers
      double precision, allocatable, dimension(:, :):: MET  ! contains our met fields
      double precision:: meanco2, meantemp, meanrad, meanprecip  ! mean conditions used in some EDCs

      ! OBS: more can obviously be added
      double precision, allocatable, dimension(:):: GPP               & ! GPP (gC/m2day)
                                                   ,NEE               & ! NEE (gC/m2/day)
                                                   ,Fire              & ! Fire (gC/m2/day)
                                                   ,LAI               & ! LAI (m2/m2)
                                                   ,Cwood_growth      & ! Wood gross increment observations (gC/m2/day)
                                                   ,Cwood_inc         & ! Wood net increment observations (gC/m2/day)
                                                   ,Cwood_mortality   & ! Natural wood mortality observations (gC/m2/day)
                                                   ,foliage_to_litter & ! Flux from foliage to leaf litter (e.g. litter trap collections) (gC/m2/day)
                                                   ,Reco              & ! Ecosystem respiration (gC/m2/day)
                                                   ,Cfol_stock        & ! time specific estimate of foliage carbon (gC/m2)
                                                   ,Cwood_stock       & ! time specific estimate of wood carbon (gC/m2)
                                                   ,Croots_stock      & ! time specific estimate of roots carbon (gC/m2)
                                                   ,Csom_stock        & ! time specific estimate if som carbon (gC/m22)
                                                   ,Cagb_stock        & ! time specific agb woody estimate (gC/m2)
                                                   ,Clit_stock        & ! time specific estimate of litter carbon (gC/m2)
                                                   ,Ccoarseroot_stock & ! time specific estimate of coarse root carbon (gC/m2)
                                                   ,Evap              & ! Evapotranspiration (kgH2O/m2/day)
                                                   ,SWE               & ! Snow Water Equivalent (mm/day)
                                                   ,NBE               & ! Net Biome Exchange (gC/m2/day)
                                                   ,fAPAR             & ! Fraction of absorbed PAR by green vegetation
                                                   ,soilwater         & ! Soil surface water content (m3/m3)
                                                   ,harvest             ! C extracted due to harvest activities (gC/m2/day)


      ! OBS uncertainties: obv these must be paired with OBS above
      double precision, allocatable, dimension(:):: GPP_unc               & ! gC/m2/day
                                                   ,NEE_unc               & ! gC/m2/day
                                                   ,Fire_unc              & ! gC/m2/day
                                                   ,LAI_unc               & ! m2/m2
                                                   ,Cwood_growth_unc      & ! gC/m2/day
                                                   ,Cwood_inc_unc         & ! gC/m2/day
                                                   ,Cwood_mortality_unc   & ! gC/m2/day
                                                   ,foliage_to_litter_unc & ! gC/m2/day
                                                   ,Reco_unc              & ! gC/m2/day
                                                   ,Cfol_stock_unc        & ! gC/m2
                                                   ,Cwood_stock_unc       & ! gC/m2
                                                   ,Croots_stock_unc      & ! gC/m2
                                                   ,Csom_stock_unc        & ! gC/m2
                                                   ,Cagb_stock_unc        & ! gC/m2
                                                   ,Clit_stock_unc        & ! gC/m2
                                                   ,Ccoarseroot_stock_unc & ! gC/m2
                                                   ,Evap_unc              & ! kgH2O/m2/day
                                                   ,SWE_unc               & ! mm/day
                                                   ,NBE_unc               & ! gC/m2/day
                                                   ,fAPAR_unc             & ! (0-1)
                                                   ,soilwater_unc         & ! (m3/m3)
                                                   ,harvest_unc             ! gC/m2/day

      ! OBS lagged period (model timestep): obs these must be paired with OBS and their uncertainties above
      integer, allocatable, dimension(:):: GPP_lag               &
                                          ,NEE_lag               &
                                          ,Fire_lag              &
                                          ,LAI_lag               &
                                          ,Cwood_growth_lag      &
                                          ,Cwood_inc_lag         &
                                          ,Cwood_mortality_lag   &
                                          ,foliage_to_litter_lag &
                                          ,Reco_lag              &
                                          ,Cfol_stock_lag        &
                                          ,Cwood_stock_lag       &
                                          ,Croots_stock_lag      &
                                          ,Csom_stock_lag        &
                                          ,Cagb_stock_lag        &
                                          ,Clit_stock_lag        &
                                          ,Ccoarseroot_stock_lag &
                                          ,Evap_lag              &
                                          ,SWE_lag               &
                                          ,NBE_lag               &
                                          ,fAPAR_lag             &
                                          ,soilwater_lag         &
                                          ,harvest_lag

      ! location of observations in the data stream, these must be paired with the above
      integer, allocatable, dimension(:):: gpppts                   & ! gpppts vector used in deriving ngpp
                                          ,neepts                   & ! same for nee
                                          ,Firepts                  & ! same for Fire
                                          ,Cwood_growthpts          & ! same for wood gross increment
                                          ,Cwood_incpts             & ! same for wood net increment
                                          ,Cwood_mortalitypts       & ! same for natural wood mortality
                                          ,foliage_to_litterpts     & ! same for foliage to litter
                                          ,laipts                   & ! same for lai
                                          ,recopts                  & ! same for ecosystem respiration
                                          ,Cfol_stockpts            & ! same for Cfoliage
                                          ,Cwood_stockpts           & ! same for Cwood
                                          ,Croots_stockpts          & ! same for Croots
                                          ,Csom_stockpts            & ! same for Csom
                                          ,Cagb_stockpts            & ! same for above ground biomass
                                          ,Clit_stockpts            & ! same for Clitter
                                          ,Ccoarseroot_stockpts     & ! same for coarse root
                                          ,Evappts                  & ! same for ecosystem evaportion
                                          ,SWEpts                   & ! same for snow water equivalent
                                          ,NBEpts                   & ! same for net biome exchange of CO2
                                          ,fAPARpts                 & ! same for fraction absorbed PAR
                                          ,soilwaterpts             & ! same for surface soil water content
                                          ,harvestpts                 ! same for C extracted due to harvest

      double precision:: nobs_scaler

      ! OBS scaling coefficients declared but calculated later
      double precision:: GPP_scaling               &
                        ,NEE_scaling               &
                        ,Fire_scaling              &
                        ,LAI_scaling               &
                        ,Cwood_growth_scaling      &
                        ,Cwood_inc_scaling         &
                        ,Cwood_mortality_scaling   &
                        ,foliage_to_litter_scaling &
                        ,Reco_scaling              &
                        ,Cfol_stock_scaling        &
                        ,Cwood_stock_scaling       &
                        ,Croots_stock_scaling      &
                        ,Csom_stock_scaling        &
                        ,Cagb_stock_scaling        &
                        ,Clit_stock_scaling        &
                        ,Ccoarseroot_stock_scaling &
                        ,Evap_scaling              &
                        ,SWE_scaling               &
                        ,NBE_scaling               &
                        ,fAPAR_scaling             &
                        ,soilwater_scaling         &
                        ,harvest_scaling

      ! counters for the number of observations per data stream
      integer:: total_obs              & ! total number of obervations
               ,ngpp                   & ! number of GPP observations
               ,nnee                   & ! number of NEE observations
               ,nFire                  & ! number of Fire observations
               ,nlai                   & ! number of LAI observations
               ,nCwood_growth          & ! number of wood gross increment obervations
               ,nCwood_inc             & ! number of wood net increment obervations
               ,nCwood_mortality       & ! number of wood mortality obervations
               ,nfoliage_to_litter     & ! number of foliage to litter obervations
               ,nreco                  & ! number of Reco observations
               ,nCfol_stock            & ! number of Cfol observations
               ,nCwood_stock           & ! number of Cwood observations
               ,nCroots_stock          & ! number of Croot observations
               ,nCsom_stock            & ! number of Csom obervations
               ,nCagb_stock            & ! number of above ground biomass
               ,nClit_stock            & ! number of Clitter observations
               ,nCcoarseroot_stock     & ! number of coarse root
               ,nEvap                  & ! number of ecosystem evaporation observations
               ,nSWE                   & ! number of snow water equivalent
               ,nNBE                   & ! number of net biome exchange of CO2
               ,nfAPAR                 & ! number of fAPAR by green vegetation
               ,nsoilwater             & ! number of surface soil water observations
               ,nharvest                 ! number of harvest observations

      double precision, dimension(:), allocatable:: soil_frac_clay, soil_frac_sand  ! clay and soil fractions of soil-
      ! initial value as read from input file.

      ! timing variable
      integer:: nos_years, steps_per_year
      double precision, allocatable, dimension(:):: deltat  ! time step (decimal day)

      double precision, allocatable, dimension(:,:):: M_FLUXES & ! All fluxes
                                                     ,M_POOLS  & ! All POOLS
                                                     ,M_DIAGS    ! All diagnostic variables
      ! static data
      integer:: nodays   & ! number of days in simulation
               ,ID       & ! model ID
               ,noobs    & ! number of obs fields
               ,nomet    & ! number met drivers
               ,nofluxes & ! number of fluxes
               ,nopools  & ! number of pools
               ,nopars   & ! number of parameters
               ,nodiags  & ! number of diagnostic variables
               ,EDC      & ! Ecological and dynamical contraints on (1) or off (0)
               ,yield    & ! yield class for ecosystem (forest only)
               ,age      & ! time in years since ecosystem established (forest only)
               ,pft        ! plant functional type information used to select appropriate DALEC submodel/ACM

      double precision:: LAT  ! site latitude

      ! binary file mcmc options (need to add all options HERE except
      ! inout files)
      integer:: edc_random_search !

      ! priors
      double precision, dimension(100):: parpriors       & ! prior values
                                        ,parpriorunc     & ! prior uncertainties
                                        ,parpriorweight    ! prior weighting
      ! other priors
      double precision, dimension(50):: otherpriors     & ! other prior values
                                       ,otherpriorunc   & ! other prior uncertainties
                                       ,otherpriorweight  ! other prior weighting

   end type DATA_type  ! DATA_type

   ! These are protected so it's not possible to write to individual elements,
   ! accidentally by using DATAin arrays as model calculation working variables.
   ! They can only be set via copy constructor set_DATAin, set_DATAin_original
   type(DATA_type), protected, save:: DATAin_original
     !! Saving a copy of DATAin as initially read from file, to never change or scale
   type(DATA_type), protected, save:: DATAin
     !! DATAin to reference thoughtout the simulation phase-may hold a scaled version
   ! shared object !  Protected (read-only), can only be set via set_datain

   type emulator_parameters

      integer ::    dim_1, & ! dimension 1 of response surface
                    dim_2, & ! dimension 2 of response surface
                nos_trees, & ! number of trees in randomForest
               nos_inputs    ! number of driver inputs

      double precision, allocatable, dimension(:, :) ::     leftDaughter, & ! left daughter for forest
                                                           rightDaughter, & ! right daughter for forets
                                                              nodestatus, & ! nodestatus for forests
                                                              xbestsplit, & ! for forest
                                                                nodepred, & ! prediction value for each tree
                                                                 bestvar    ! for randomForests

   end type emulator_parameters  ! emulator parameters
   type(emulator_parameters), protected, save:: emulator_pars  ! TODO make sure not shared, or read-only

contains
  !
  !--------------------------------------------------------------------
  !
  subroutine set_datain_original(datain_source)

    !! A setter, copying the argument to cardamom_structures:: DATAin_original
    !! The central DATAin in module cardamom_structures can ONLY be set by the constructor,
    !! this ensures that no model calculations are writing to its elements from different parallel threads
    type(DATA_type), intent(in):: datain_source
    DATAin_original = datain_source

  end subroutine set_datain_original
  !
  !--------------------------------------------------------------------
  !
  subroutine set_datain(datain_source)

    !! A setter, copying the argument to cardamom_structures:: DATAin
    !! The central DATAin in module cardamom_structures can ONLY be set by the constructor,
    !! this ensures that no model calculations are writing to its elements from different parallel threads
    type(DATA_type), intent(in):: datain_source
    DATAin = datain_source

  end subroutine set_datain
  !
  !--------------------------------------------------------------------
  !
end module cardamom_structures
