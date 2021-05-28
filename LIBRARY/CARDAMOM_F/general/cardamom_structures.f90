module cardamom_structures

 !!!!!!!!!!!
 ! Authorship contributions
 !
 ! This code is based on the original C verion of the University of Edinburgh
 ! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
 ! All code translation into Fortran, integration into the University of
 ! Edinburgh CARDAMOM code and subsequent modifications by:
 ! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
 ! J. F. Exbrayat (University of Edinburgh)
 ! See function / subroutine specific comments for exceptions and contributors
 !!!!!!!!!!!

implicit none

private

public :: data_type, DATAin, emulator_parameters, emulator_pars, io_space

  !!!!! such as the data type !!!!!
  type DATA_type

      ! drivers
      double precision, allocatable, dimension(:,:) :: MET ! contains our met fields
      double precision :: meanco2, meantemp, meanrad, meanprecip ! mean conditions used in some EDCs

      ! OBS: more can obviously be added
      double precision, allocatable, dimension(:) :: GPP     & ! GPP (gC.m-2.day-1)
                                          ,NEE               & ! NEE (gC.m-2.day-1)
                                          ,LAI               & ! LAI (m2/m2)
                                          ,WOO               & ! Wood increment observations (gC.m-2.yr-1)
                                          ,Reco              & ! Ecosystem respiration (gC.m-2.day-1)
                                          ,Cfol_stock        & ! time specific estimate of foliage carbon (gC.m-2)
                                          ,Cwood_stock       & ! time specific estimate of wood carbon (gC.m-2)
                                          ,Croots_stock      & ! time specific estimate of roots carbon (gC.m-2)
                                          ,Csom_stock        & ! time specific estimate if som carbon (gC.m-2)
                                          ,Cagb_stock        & ! time specific agb woody estimate (gC.m-2)
                                          ,Clit_stock        & ! time specific estimate of litter carbon (gC.m-2)
                                          ,Ccoarseroot_stock & ! time specific estimate of coarse root carbon (gC.m-2)
                                          ,Cfolmax_stock     & ! maximum annual foliar stock (gC.m-2)
                                          ,Evap              & ! Evapotranspiration (kg.m-2.day-1)
                                          ,SWE               & ! Snow Water Equivalent (mm.day-1)
                                          ,NBE                 ! Net Biome Exchange (gC/m2/day)

      ! OBS uncertainties: obv these must be pared with OBS above
      double precision, allocatable, dimension(:) :: GPP_unc     & ! (gC/m2/day)
                                          ,NEE_unc               & ! (gC/m2/day)
                                          ,LAI_unc               & ! (m2/m2)
                                          ,WOO_unc               & ! (gC/m2/yr)
                                          ,Reco_unc              & ! (gC/m2/day)
                                          ,Cfol_stock_unc        & ! gC/m2
                                          ,Cwood_stock_unc       & ! gC/m2
                                          ,Croots_stock_unc      & ! gC/m2
                                          ,Csom_stock_unc        & ! gC/m2
                                          ,Cagb_stock_unc        & ! gC/m2
                                          ,Clit_stock_unc        & ! gC/m2
                                          ,Ccoarseroot_stock_unc & ! gC/m2
                                          ,Cfolmax_stock_unc     & ! gC/m2
                                          ,Evap_unc              & ! (kg.m-2.day-1)
                                          ,SWE_unc               & ! (mm.day-1)
                                          ,NBE_unc                 ! gC/m2/day

      ! location of observations in the data stream
      integer, allocatable, dimension(:) :: gpppts                   & ! gpppts vector used in deriving ngpp
                                           ,neepts                   & ! same for nee
                                           ,woopts                   & ! same for wood
                                           ,laipts                   & ! same for lai
                                           ,recopts                  & ! same for ecosystem respiration
                                           ,Cfol_stockpts            & ! same for Cfoliage
                                           ,Cwood_stockpts           & ! same for Cwood
                                           ,Croots_stockpts          & ! same for Croots
                                           ,Csom_stockpts            & ! same for Csom
                                           ,Cagb_stockpts            & ! same for above ground biomass
                                           ,Clit_stockpts            & ! same for Clitter
                                           ,Ccoarseroot_stockpts     & ! same for coarse root
                                           ,Cfolmax_stockpts         & ! same for seasonal max foliar
                                           ,Evappts                  & ! same for ecosystem evaportion
                                           ,SWEpts                   & ! same for snow water equivalent
                                           ,NBEpts                   & ! same for net biome exchange of CO2
                                           ,sub_gpppts               & ! as above but for a sub-sample of observations
                                           ,sub_neepts               & !
                                           ,sub_woopts               & !
                                           ,sub_laipts               & !
                                           ,sub_recopts              & !
                                           ,sub_Cfol_stockpts        & !
                                           ,sub_Cwood_stockpts       & !
                                           ,sub_Croots_stockpts      & !
                                           ,sub_Csom_stockpts        & !
                                           ,sub_Cagb_stockpts        & !
                                           ,sub_Clit_stockpts        & !
                                           ,sub_Ccoarseroot_stockpts & !
                                           ,sub_Cfolmax_stockpts     & !
                                           ,sub_Evappts              & !
                                           ,sub_SWEpts               & !
                                           ,sub_NBEpts

      double precision :: nobs_scaler

      ! counters for the number of observations per data stream
      integer :: total_obs              & ! total number of obervations
                ,ngpp                   & ! number of GPP observations
                ,nnee                   & ! number of NEE observations
                ,nlai                   & ! number of LAI observations
                ,nwoo                   & ! number of wood increment obervations
                ,nreco                  & ! number of Reco observations
                ,nCfol_stock            & ! number of Cfol observations
                ,nCwood_stock           & ! number of Cwood observations
                ,nCroots_stock          & ! number of Croot observations
                ,nCsom_stock            & ! number of Csom obervations
                ,nCagb_stock            & ! number of above ground biomass
                ,nClit_stock            & ! number of Clitter observations
                ,nCcoarseroot_stock     & ! number of coarse root
                ,nCfolmax_stock         & ! number of seasonal maximum foliar C
                ,nEvap                  & ! number of ecosystem evaporation observations
                ,nSWE                   & ! number of snow water equivalent
                ,nNBE                   & ! number of net biome exchange of CO2
                ,sub_ngpp               & ! As above but for a subsample of observations
                ,sub_nnee               & !
                ,sub_nlai               & !
                ,sub_nwoo               & !
                ,sub_nreco              & !
                ,sub_nCfol_stock        & !
                ,sub_nCwood_stock       & !
                ,sub_nCroots_stock      & !
                ,sub_nCsom_stock        & !
                ,sub_nCagb_stock        & !
                ,sub_nClit_stock        & !
                ,sub_nCcoarseroot_stock & !
                ,sub_nCfolmax_stock     & !
                ,sub_nEvap              &
                ,sub_nSWE               &
                ,sub_nNBE

      ! saving computational speed by allocating memory to model output
      double precision, allocatable, dimension(:) :: M_GPP    & !
                                                    ,M_NEE    & !
                                                    ,M_LAI      !
      ! timing variable
      double precision, allocatable, dimension(:) :: deltat ! time step (decimal day)

      double precision, allocatable, dimension(:,:) :: M_FLUXES & !
                                                      ,M_POOLS  & !
                                                      ,C_POOLS    !
      ! static data
      integer :: nodays   & ! number of days in simulation
                ,ID       & ! model ID, currently 1=DALEC_CDEA, 2=DALEC_BUCKET
                ,noobs    & ! number of obs fields
                ,nomet    & ! number met drivers
                ,nofluxes & ! number of fluxes
                ,nopools  & ! number of pools
                ,nopars   & ! number of parameters
                ,EDC      & ! Ecological and dynamical contraints on (1) or off (0)
                ,yield    & ! yield class for ecosystem (forest only)
                ,age      & ! time in years since ecosystem established (forest only)
                ,pft        ! plant functional type information used to select appropriate DALEC submodel / ACM

      double precision :: LAT ! site latitude

      ! binary file mcmc options (need to add all options HERE except
      ! inout files)
      integer :: edc_random_search !

      ! priors
      double precision, dimension(100) :: parpriors     & ! prior values
                                         ,parpriorunc   & ! prior uncertainties
                                         ,otherpriors   & ! other prior values
                                         ,otherpriorunc   ! other prior uncertainties

  end type ! DATA_type
  type (DATA_type), save :: DATAin

  type io_buffer_space

    integer :: io_buffer, io_buffer_count
    double precision, allocatable, dimension(:,:) :: &
                                    variance_buffer, &
                                   mean_pars_buffer, &
                                        pars_buffer

    double precision, allocatable, dimension(:) :: &
                                   nsample_buffer, &
                               accept_rate_buffer, &
                                      prob_buffer

  end type
  type(io_buffer_space), save :: io_space

  type emulator_parameters

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

  end type ! emulator parameters
  type (emulator_parameters), save :: emulator_pars

end module cardamom_structures
