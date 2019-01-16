
module CARBON_MODEL_MOD

implicit none

! make all private
private

! explicit publics
public :: CARBON_MODEL     &
         ,soil_frac_clay   &
         ,soil_frac_sand   &
         ,nos_soil_layers  &
         ,disturbance_residue_to_litter &
         ,disturbance_residue_to_som    &
         ,itemp,ivpd,iphoto&
         ,extracted_C      &
         ,dim_1,dim_2      &
         ,nos_trees        &
         ,nos_inputs       &
         ,leftDaughter     &
         ,rightDaughter    &
         ,nodestatus       &
         ,xbestsplit       &
         ,nodepred         &
         ,bestvar

! ACM related parameters
double precision, parameter :: pi = 3.1415927
double precision, parameter :: deg_to_rad = pi/180d0

! forest rotation specific info
double precision, allocatable, dimension(:) :: extracted_C,itemp,ivpd,iphoto
! local variables for GSI phenology model
double precision :: Tfac,Photofac,VPDfac & ! oC, seconds, Pa
                   ,delta_gsi,tmp,gradient &
                   ,fol_turn_crit,lab_turn_crit &
                   ,gsi_history(22),just_grown

! arrays for the emulator, just so we load them once and that is it cos they be
! massive
integer ::    dim_1, & ! dimension 1 of response surface
              dim_2, & ! dimension 2 of response surface
          nos_trees, & ! number of trees in randomForest
         nos_inputs    ! number of driver inputs

double precision, allocatable, dimension(:) :: disturbance_residue_to_litter, &
                                               disturbance_residue_to_som
integer :: gsi_lag_remembered 
double precision, allocatable, dimension(:) :: tmp_x, tmp_m
double precision, allocatable, dimension(:,:) ::     leftDaughter, & ! left daughter for forest
                                                    rightDaughter, & ! right daughter for forets
                                                       nodestatus, & ! nodestatus for forests
                                                       xbestsplit, & ! for forest
                                                         nodepred, & ! prediction value for each tree
                                                          bestvar    ! for randomForests

! for consisteny between requirements of different models
integer, parameter :: nos_root_layers = 2, nos_soil_layers = nos_root_layers + 1
double precision, dimension(nos_soil_layers) :: soil_frac_clay,soil_frac_sand

contains
!
!--------------------------------------------------------------------
!
  subroutine CARBON_MODEL(start,finish,met,pars,deltat,nodays,lat,lai,NEE,FLUXES,POOLS &
                         ,nopars,nomet,nopools,nofluxes,GPP)

    ! The Data Assimilation Linked Ecosystem Carbon - Growing Season
    ! Index - Forest Rotation (DALEC_GSI_FR) model. 
    ! The subroutine calls the Aggregated Canopy Model to simulate GPP and 
    ! partitions between various ecosystem carbon pools. These pools are
    ! subject to turnovers / decompostion resulting in ecosystem phenology and fluxes of CO2

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

    double precision, dimension(nodays), intent(inout) :: lai & ! leaf area index
                                               ,GPP & ! Gross primary productivity
                                               ,NEE   ! net ecosystem exchange of CO2

    double precision, dimension((nodays+1),nopools), intent(inout) :: POOLS ! vector of ecosystem pools
 
    double precision, dimension(nodays,nofluxes), intent(inout) :: FLUXES ! vector of ecosystem fluxes
                                             
    ! declare general local variables
    double precision :: gpppars(12)        & ! ACM inputs (LAI+met)
             ,constants(10)                  ! parameters for ACM

    integer :: p,f,nxp,n,test,m

    ! local fire related variables
    double precision :: CFF(6) = 0, CFF_res(4) = 0    & ! combusted and non-combustion fluxes
                       ,NCFF(6) = 0, NCFF_res(4) = 0  & ! with residue and non-residue seperates
                       ,combust_eff(5)                & ! combustion efficiency
                       ,rfac                            ! resilience factor

    ! local deforestation related variables
    double precision, dimension(4) :: post_harvest_burn   & ! how much burning to occur after
                                     ,foliage_frac_res    &
                                     ,roots_frac_res      &
                                     ,rootcr_frac_res     &
                                     ,stem_frac_res       &
                                     ,branch_frac_res     &
                                     ,Cbranch_part        &
                                     ,Crootcr_part        &
                                     ,soil_loss_frac     

    double precision :: labile_loss,foliar_loss      &
                       ,roots_loss,wood_loss         &
                       ,labile_residue,foliar_residue&
                       ,roots_residue,wood_residue   & 
                       ,wood_pellets,C_total         &
                       ,labile_frac_res              &
                       ,Cstem,Cbranch,Crootcr        &
                       ,stem_residue,branch_residue  &
                       ,coarse_root_residue          &
                       ,soil_loss_with_roots  

    integer :: reforest_day, harvest_management,restocking_lag, gsi_lag

    ! met drivers are:
    ! 1st run day
    ! 2nd min daily temp (oC)
    ! 3rd max daily temp (oC)
    ! 4th Radiation (MJ.m-2.day-1)
    ! 5th CO2 (ppm)
    ! 6th DOY
    ! 7th lagged precip
    ! 8th deforestation fraction
    ! 9th burnt area fraction
    ! 10th 21 day average min temperature (K)
    ! 11th 21 day average photoperiod (sec)
    ! 12th 21 day average VPD (Pa)
    ! 13th Forest management practice to accompany any clearing
    ! 14th Mean temperature

    ! POOLS are:
    ! 1 = labile (p18)
    ! 2 = foliar (p19)
    ! 3 = root   (p20)
    ! 4 = wood   (p21)
    ! 5 = litter (p22)
    ! 6 = som    (p23)

    ! p(30) = labile replanting
    ! p(31) = foliar replanting
    ! p(32) = fine root replanting
    ! p(33) = wood replanting

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
    ! 17 = carbon flux due to fire
    ! 18 = growing season index

    ! PARAMETERS
    ! 17+4(GSI) values

    ! p(1) Litter to SOM conversion rate  - m_r
    ! p(2) Fraction of GPP respired - f_a
    ! p(3) Fraction of NPP allocated to foliage - f_f 
    ! p(4) Fraction of NPP allocated to roots - f_r
    ! p(5) max leaf turnover (GSI) ! Leaf lifespan - L_f (CDEA)
    ! p(6) Turnover rate of wood - t_w
    ! p(7) Turnover rate of roots - t_r
    ! p(8) Litter turnover rate - t_l
    ! p(9) SOM turnover rate  - t_S
    ! p(10) Parameter in exponential term of temperature - \theta
    ! p(11) Canopy efficiency parameter - C_eff (part of ACM)
    ! p(12) = max labile turnover(GSI) ! date of Clab release - B_day (CDEA)
    ! p(13) = Fraction allocated to Clab - f_l
    ! p(14) = min temp threshold (GSI) ! lab release duration period - R_l (CDEA)
    ! p(15) = max temp threshold (GSI)! date of leaf fall - F_day
    ! p(16) = min photoperiod threshold (GIS) 
    ! p(17) = LMA
    ! p(24) = max photoperiod threshold (GSI)
    ! p(25) = min VPD threshold (GSI)
    ! p(26) = max VPD threshold (GSI)
    ! p(27) = minimum GPP benefit of increased LAI for labile allocation to be allowed
    ! p(28) = fraction of Cwood which is Cbranch
    ! p(29) = fraction of Cwood which is Ccoarseroot

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

    ! load some values
    gpppars(4) = 1.0 ! 10d0**(pars(11)) !TLS 1 ! foliar N
    gpppars(7) = lat
    gpppars(9) = -2.0 ! leafWP-soilWP
    gpppars(10) = 1.0 ! totaly hydraulic resistance
    gpppars(11) = pi

    ! assign acm parameters
    constants(1)=pars(11) 
    constants(2)=0.0156935
    constants(3)=4.22273
    constants(4)=208.868
    constants(5)=0.0453194
    constants(6)=0.37836
    constants(7)=7.19298
    constants(8)=0.011136
    constants(9)=2.1001
    constants(10)=0.789798

    ! initial values for deforestation variables
    labile_loss = 0.    ; foliar_loss = 0.
    roots_loss = 0.     ; wood_loss = 0.
    labile_residue = 0. ; foliar_residue = 0.
    roots_residue = 0.  ; wood_residue = 0.
    stem_residue = 0.   ; branch_residue = 0.
    reforest_day = 0
    soil_loss_with_roots = 0.
    coarse_root_residue = 0.
    post_harvest_burn = 0.

    ! now load the hardcoded forest management parameters into their locations

    ! Parameter values for deforestation variables
    ! scenario 1
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(1) = 1.0
    roots_frac_res(1)   = 1.0
    rootcr_frac_res(1) = 1.0
    branch_frac_res(1) = 1.0
    stem_frac_res(1)   = 0. ! 
    ! wood partitioning (fraction)
    Crootcr_part(1) = 0.32 ! Coarse roots (Adegbidi et al 2005;
    ! Black et al 2009; Morison et al 2012)
    Cbranch_part(1) =  0.20 ! (Ares & Brauers 2005)
    ! actually < 15 years branches = ~25 %
    !          > 15 years branches = ~15 %.
    ! Csom loss due to phyical removal with roots 
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(1) = 0.02 ! actually between 1-3 %
    ! was the forest burned after deforestation
    post_harvest_burn(1) = 1. 

    !## scen 2
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(2) = 1.0
    roots_frac_res(2)   = 1.0
    rootcr_frac_res(2) = 1.0
    branch_frac_res(2) = 1.0
    stem_frac_res(2)   = 0. ! 
    ! wood partitioning (fraction)
    Crootcr_part(2) = 0.32 ! Coarse roots (Adegbidi et al 2005;
    ! Black et al 2009; Morison et al 2012)
    Cbranch_part(2) =  0.20 ! (Ares & Brauers 2005)
    ! actually < 15 years branches = ~25 %
    !          > 15 years branches = ~15 %.
    ! Csom loss due to phyical removal with roots 
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(2) = 0.02 ! actually between 1-3 %
    ! was the forest burned after deforestation
    post_harvest_burn(2) = 0.

    !## scen 3
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(3) = 0.5
    roots_frac_res(3)   = 1.0
    rootcr_frac_res(3) = 1.0
    branch_frac_res(3) = 0.
    stem_frac_res(3)   = 0. ! 
    ! wood partitioning (fraction)
    Crootcr_part(3) = 0.32 ! Coarse roots (Adegbidi et al 2005;
    ! Black et al 2009; Morison et al 2012)
    Cbranch_part(3) =  0.20 ! (Ares & Brauers 2005)
    ! actually < 15 years branches = ~25 %
    !          > 15 years branches = ~15 %.
    ! Csom loss due to phyical removal with roots 
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(3) = 0.02 ! actually between 1-3 %
    ! was the forest burned after deforestation
    post_harvest_burn(3) = 0.

    !## scen 4
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(4) = 0.5
    roots_frac_res(4)   = 1.0
    rootcr_frac_res(4) = 0.
    branch_frac_res(4) = 0.
    stem_frac_res(4)   = 0. ! 
    ! wood partitioning (fraction)
    Crootcr_part(4) = 0.32 ! Coarse roots (Adegbidi et al 2005;
    ! Black et al 2009; Morison et al 2012)
    Cbranch_part(4) =  0.20 ! (Ares & Brauers 2005)
    ! actually < 15 years branches = ~25 %
    !          > 15 years branches = ~15 %.
    ! Csom loss due to phyical removal with roots 
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(4) = 0.02 ! actually between 1-3 %
    ! was the forest burned after deforestation
    post_harvest_burn(4) = 0.

    ! for the moment override all paritioning parameters with those coming from
    ! CARDAMOM
    Cbranch_part = pars(28)
    Crootcr_part = pars(29)

    if (start == 1) then

        ! assigning initial conditions
        POOLS(1,1)=pars(18)
        POOLS(1,2)=pars(19)
        POOLS(1,3)=pars(20)
        POOLS(1,4)=pars(21)
        POOLS(1,5)=pars(22)
        POOLS(1,6)=pars(23)

        ! declare fire constants (labile, foliar, roots, wood, litter)
        combust_eff(1) = 0.1 ; combust_eff(2) = 0.9
        combust_eff(3) = 0.1 ; combust_eff(4) = 0.5
        combust_eff(5) = 0.3 ; rfac = 0.5

        if (.not.allocated(disturbance_residue_to_som)) then
            allocate(disturbance_residue_to_litter(nodays), &
                     disturbance_residue_to_som(nodays))
        endif
        disturbance_residue_to_litter = 0d0
        disturbance_residue_to_som = 0d0

        ! calculate some values once as these are invarient between DALEC runs
        if (.not.allocated(tmp_x)) then
            ! 21 days is the maximum potential so we will fill the maximum potential
            ! + 1 for safety
            allocate(tmp_x(22),tmp_m(nodays))
            do f = 1, 22
               tmp_x(f) = f
            end do
            do n = 1, nodays
              ! calculate the gradient / trend of GSI
              if (sum(deltat(1:n)) < 21) then
                  tmp_m(n) = n-1
              else
                 ! else we will try and work out the gradient to see what is
                 ! happening
                 ! to the system over all. The default assumption will be to
                 ! consider
                 ! the averaging period of GSI model (i.e. 21 days). If this is not
                 ! possible either the time step of the system is used (if step
                 ! greater
                 ! than 21 days) or all available steps (if n < 21).
                 m = 0 ; test = 0
                 do while (test < 21)
                    m=m+1 ; test = sum(deltat((n-m):n))
                        if (m > (n-1)) test = 21
                 end do
                 tmp_m(n) = m
               endif ! for calculating gradient
            end do ! calc daily values once
            ! allocate GSI history dimension
            gsi_lag_remembered=max(2,maxval(nint(tmp_m)))
        end if ! .not.allocated(tmp_x)
        ! assign our starting value
        gsi_history = pars(36)-1d0
        just_grown=pars(35)

    endif ! start == 1

    ! assign climate sensitivities
    gsi_lag = gsi_lag_remembered ! added to prevent loss from memory
    fol_turn_crit=pars(34)-1d0
    lab_turn_crit=pars(3)-1d0

    ! 
    ! Begin looping through each time step
    ! 

    do n = start, finish
  
      ! calculate LAI value
      lai(n)=POOLS(n,2)/pars(17)

      ! load next met / lai values for ACM
      gpppars(1)=lai(n)
      gpppars(2)=met(3,n) ! max temp
      gpppars(3)=met(2,n) ! min temp
      gpppars(5)=met(5,n) ! co2
      gpppars(6)=ceiling(met(6,n)-(deltat(n)*0.5)) ! doy
      gpppars(8)=met(4,n) ! radiation

      ! GPP (gC.m-2.day-1)
      if (lai(n) > 0.) then
         FLUXES(n,1) = acm(gpppars,constants)
      else
         FLUXES(n,1) = 0.
      endif 
      ! temprate (i.e. temperature modified rate of metabolic activity))
      FLUXES(n,2) = exp(pars(10)*0.5*(met(3,n)+met(2,n)))
      ! autotrophic respiration (gC.m-2.day-1)
      FLUXES(n,3) = pars(2)*FLUXES(n,1)
      ! leaf production rate (gC.m-2.day-1)
      FLUXES(n,4) = 0. !(FLUXES(n,1)-FLUXES(n,3))*pars(3)
      ! labile production (gC.m-2.day-1)
      FLUXES(n,5) = (FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,4))*pars(13)
      ! root production (gC.m-2.day-1)
      FLUXES(n,6) = (FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,4)-FLUXES(n,5))*pars(4)
      ! wood production 
      FLUXES(n,7) = FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,4)-FLUXES(n,5)-FLUXES(n,6)

      ! GSI added to fortran version by TLS 24/11/2014
      ! /* 25/09/14 - JFE
      ! Here we calculate the Growing Season Index based on 
      ! Jolly et al. A generalized, bioclimatic index to predict foliar 
      ! phenology in response to climate Global Change Biology, Volume 11, page 619-632,
      ! 2005 (doi: 10.1111/j.1365-2486.2005.00930.x) 
      ! Stoeckli, R., T. Rutishauser, I. Baker, M. A. Liniger, and A. S.
      ! Denning (2011), A global reanalysis of vegetation phenology, J. Geophys. Res.,
      ! 116, G03020, doi:10.1029/2010JG001545.
        
      ! It is the product of 3 limiting factors for temperature, photoperiod and
      ! vapour pressure deficit that grow linearly from 0 to 1 between a calibrated 
      ! min and max value. Photoperiod, VPD and avgTmin are direct input

      ! temperature limitation, then restrict to 0-1; correction for k-> oC
      Tfac = (met(10,n)-(pars(14)-273.15)) / (pars(15)-pars(14))
      Tfac = min(1d0,max(0d0,Tfac))
      ! photoperiod limitation
      Photofac = (met(11,n)-pars(16)) / (pars(24)-pars(16))
      Photofac = min(1d0,max(0d0,Photofac))
      ! VPD limitation
      VPDfac = 1.0 - ( (met(12,n)-pars(25)) / (pars(26)-pars(25)) )
      VPDfac = min(1d0,max(0d0,VPDfac))

      ! calculate and store the GSI index
      FLUXES(n,18) = Tfac*Photofac*VPDfac

      ! we will load up some needed variables
      m = tmp_m(n)
      ! update gsi_history for the calculation
      if (n == 1) then
          ! in first step only we want to take the initial GSI value only
          gsi_history(gsi_lag) = FLUXES(n,18)
      else
          gsi_history((gsi_lag-m):gsi_lag) = FLUXES((n-m):n,18)
      endif
      ! calculate gradient
      gradient = linear_model_gradient(tmp_x(1:(gsi_lag)),gsi_history(1:gsi_lag),gsi_lag)
      ! adjust gradient to daily rate
      gradient = gradient /  nint((sum(deltat((n-m+1):n))) / (gsi_lag-1))
      gsi_lag_remembered = gsi_lag

      ! first assume that nothing is happening
      FLUXES(n,9) = 0d0  ! leaf turnover
      FLUXES(n,16) = 0d0 ! leaf growth

      ! now update foliage and labile conditions based on gradient calculations
      if (gradient < fol_turn_crit .or. FLUXES(n,18) == 0) then
         ! we are in a decending condition so foliar turnover
         FLUXES(n,9) = pars(5)*(1.0-FLUXES(n,18))
         just_grown = 0.5
      else if (gradient > lab_turn_crit) then
         ! we are in a assending condition so labile turnover
         FLUXES(n,16) = pars(12)*FLUXES(n,18)
         just_grown = 1.5
         ! check carbon return
         tmp = POOLS(n,1)*(1d0-(1d0-FLUXES(n,16))**deltat(n))/deltat(n)
         tmp = (POOLS(n,2)+tmp)/pars(17)
         gpppars(1)=tmp
         tmp = acm(gpppars,constants)
         ! determine if increase in LAI leads to an improvement in GPP greater
         ! than
         ! critical value, if not then no labile turnover allowed
         if ( ((tmp - FLUXES(n,1))/FLUXES(n,1)) < pars(27) ) then
             FLUXES(n,16) = 0d0
         endif
      else
         ! probaly we want nothing to happen, however if we are at the seasonal
         ! maximum we will consider further growth still
         if (just_grown >= 1.0) then
            ! we are between so definitely not losing foliage and we have
            ! previously been growing so maybe we still have a marginal return on
            ! doing so again
            FLUXES(n,16) = pars(12)*FLUXES(n,18)
            ! but possibly gaining some?
            ! determine if this is a good idea based on GPP increment
            tmp = POOLS(n,1)*(1d0-(1d0-FLUXES(n,16))**deltat(n))/deltat(n)
            tmp = (POOLS(n,2)+tmp)/pars(17)
            gpppars(1)=tmp
            tmp = acm(gpppars,constants)
            ! determine if increase in LAI leads to an improvement in GPP greater
            ! than
            ! critical value, if not then no labile turnover allowed
            if ( ((tmp - FLUXES(n,1))/FLUXES(n,1)) < pars(27) ) then
                FLUXES(n,16) = 0d0
            endif
         end if ! Just grown?
      endif ! gradient choice

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
      FLUXES(n,8)  = POOLS(n,1)*(1.-(1.-FLUXES(n,16))**deltat(n))/deltat(n)
      ! total leaf litter production
      FLUXES(n,10) = POOLS(n,2)*(1.-(1.-FLUXES(n,9))**deltat(n))/deltat(n)
      ! total wood litter production
      FLUXES(n,11) = POOLS(n,4)*(1.-(1.-pars(6))**deltat(n))/deltat(n)
      ! total root litter production
      FLUXES(n,12) = POOLS(n,3)*(1.-(1.-pars(7))**deltat(n))/deltat(n)

      ! 
      ! those with temperature AND time dependancies
      ! 

      ! respiration heterotrophic litter
      FLUXES(n,13) = POOLS(n,5)*(1.-(1.-FLUXES(n,2)*pars(8))**deltat(n))/deltat(n)
      ! respiration heterotrophic som
      FLUXES(n,14) = POOLS(n,6)*(1.-(1.-FLUXES(n,2)*pars(9))**deltat(n))/deltat(n)
      ! litter to som
      FLUXES(n,15) = POOLS(n,5)*(1.-(1.-pars(1)*FLUXES(n,2))**deltat(n))/deltat(n)

      ! calculate the NEE 
      NEE(n) = (-FLUXES(n,1)+FLUXES(n,3)+FLUXES(n,13)+FLUXES(n,14))
      ! load GPP
      GPP(n) = FLUXES(n,1)

      !
      ! update pools for next timestep
      ! 

      ! labile pool
      POOLS(n+1,1) = POOLS(n,1) + (FLUXES(n,5)-FLUXES(n,8))*deltat(n)
      ! foliar pool
      POOLS(n+1,2) = POOLS(n,2) + (FLUXES(n,4)-FLUXES(n,10) + FLUXES(n,8))*deltat(n)
      ! wood pool
      POOLS(n+1,4) = POOLS(n,4) + (FLUXES(n,7)-FLUXES(n,11))*deltat(n)
      ! root pool
      POOLS(n+1,3) = POOLS(n,3) + (FLUXES(n,6)-FLUXES(n,12))*deltat(n)
      ! litter pool
      POOLS(n+1,5) = POOLS(n,5) + (FLUXES(n,10)+FLUXES(n,12)-FLUXES(n,13)-FLUXES(n,15))*deltat(n)
      ! som pool
      POOLS(n+1,6) = POOLS(n,6) + (FLUXES(n,15)-FLUXES(n,14)+FLUXES(n,11))*deltat(n)

      ! 
      ! deal first with deforestation
      ! 

      if (n == reforest_day) then
          POOLS(n+1,1) = pars(30) 
          POOLS(n+1,2) = pars(31) 
          POOLS(n+1,3) = pars(32) 
          POOLS(n+1,4) = pars(33) 
      endif 

      if (met(8,n) > 0.) then

          ! pass harvest management to local integer
          harvest_management = int(met(13,n))

          ! assume that labile is proportionally distributed through the plant
          ! and therefore so is the residual fraction
          C_total = POOLS(n+1,2) + POOLS(n+1,3) + POOLS(n+1,4)
          ! partition wood into its components
          Cbranch = POOLS(n+1,4)*Cbranch_part(harvest_management)
          Crootcr = POOLS(n+1,4)*Crootcr_part(harvest_management)
          Cstem   = POOLS(n+1,4)-(Cbranch + Crootcr)
          ! now calculate the labile fraction of residue
          labile_frac_res = ( (POOLS(n+1,2)/C_total) * foliage_frac_res(harvest_management) ) & 
                          + ( (POOLS(n+1,3)/C_total) * roots_frac_res(harvest_management)   ) & 
                          + ( (Cbranch/C_total)      * branch_frac_res(harvest_management)  ) &
                          + ( (Cstem/C_total)        * stem_frac_res(harvest_management)    ) &
                          + ( (Crootcr/C_total)      * rootcr_frac_res(harvest_management)  ) 

          ! loss of carbon from each pools
          labile_loss = POOLS(n+1,1)*met(8,n)
          foliar_loss = POOLS(n+1,2)*met(8,n)
          roots_loss  = POOLS(n+1,3)*met(8,n)
          wood_loss   = POOLS(n+1,4)*met(8,n)
          ! transfer fraction of harvest waste to litter or som pools
          ! easy pools first
          labile_residue = POOLS(n+1,1)*met(8,n)*labile_frac_res
          foliar_residue = POOLS(n+1,2)*met(8,n)*foliage_frac_res(harvest_management)
          roots_residue  = POOLS(n+1,3)*met(8,n)*roots_frac_res(harvest_management)
          ! explicit calculation of the residues from each fraction
          coarse_root_residue  = Crootcr*met(8,n)*rootcr_frac_res(harvest_management)
          branch_residue = Cbranch*met(8,n)*branch_frac_res(harvest_management)
          stem_residue = Cstem*met(8,n)*stem_frac_res(harvest_management)
          ! now finally calculate the final wood residue
          wood_residue = stem_residue + branch_residue + coarse_root_residue 
          ! mechanical loss of Csom due to coarse root extraction                 
          soil_loss_with_roots = Crootcr*met(8,n)*(1.-rootcr_frac_res(harvest_management)) &
                              * soil_loss_frac(harvest_management)

          ! update living pools directly
          POOLS(n+1,1) = max(0.,POOLS(n+1,1)-labile_loss)
          POOLS(n+1,2) = max(0.,POOLS(n+1,2)-foliar_loss)
          POOLS(n+1,3) = max(0.,POOLS(n+1,3)-roots_loss)
          POOLS(n+1,4) = max(0.,POOLS(n+1,4)-wood_loss)
          ! then work out the adjustment due to burning if there is any
          if (post_harvest_burn(harvest_management) > 0.) then
              !/*first fluxes*/
              !/*LABILE*/
              CFF(1) = POOLS(n+1,1)*post_harvest_burn(harvest_management)*combust_eff(1)
              NCFF(1) = POOLS(n+1,1)*post_harvest_burn(harvest_management)*(1-combust_eff(1))*(1-rfac)
              CFF_res(1) = labile_residue*post_harvest_burn(harvest_management)*combust_eff(1)
              NCFF_res(1) = labile_residue*post_harvest_burn(harvest_management)*(1-combust_eff(1))*(1-rfac)
              !/*foliar*/
              CFF(2) = POOLS(n+1,2)*post_harvest_burn(harvest_management)*combust_eff(2)
              NCFF(2) = POOLS(n+1,2)*post_harvest_burn(harvest_management)*(1-combust_eff(2))*(1-rfac)
              CFF_res(2) = foliar_residue*post_harvest_burn(harvest_management)*combust_eff(2)
              NCFF_res(2) = foliar_residue*post_harvest_burn(harvest_management)*(1-combust_eff(2))*(1-rfac)
              !/*root*/
              CFF(3) = 0. !POOLS(n+1,3)*post_harvest_burn(harvest_management)*combust_eff(3)
              NCFF(3) = 0. !POOLS(n+1,3)*post_harvest_burn(harvest_management)*(1-combust_eff(3))*(1-rfac)
              CFF_res(3) = 0. !roots_residue*post_harvest_burn(harvest_management)*combust_eff(3)
              NCFF_res(3) = 0. !roots_residue*post_harvest_burn(harvest_management)*(1-combust_eff(3))*(1-rfac)
              !/*wood*/
              CFF(4) = POOLS(n+1,4)*post_harvest_burn(harvest_management)*combust_eff(4)
              NCFF(4) = POOLS(n+1,4)*post_harvest_burn(harvest_management)*(1-combust_eff(4))*(1-rfac)
              CFF_res(4) = wood_residue*post_harvest_burn(harvest_management)*combust_eff(4)
              NCFF_res(4) = wood_residue*post_harvest_burn(harvest_management)*(1-combust_eff(4))*(1-rfac)
              !/*litter*/
              CFF(5) = POOLS(n+1,5)*post_harvest_burn(harvest_management)*combust_eff(5)
              NCFF(5) = POOLS(n+1,5)*post_harvest_burn(harvest_management)*(1-combust_eff(5))*(1-rfac)
              !/*fires as daily averages to comply with units*/
              FLUXES(n,17)=(CFF(1)+CFF(2)+CFF(3)+CFF(4)+CFF(5) & 
                           +CFF_res(1)+CFF_res(2)+CFF_res(3)+CFF_res(4))/deltat(n)
              ! update the residue terms
              labile_residue = labile_residue - CFF_res(1) - NCFF_res(1)
              foliar_residue = foliar_residue - CFF_res(2) - NCFF_res(2)
              roots_residue  = roots_residue  - CFF_res(3) - NCFF_res(3)
              wood_residue   = wood_residue   - CFF_res(4) - NCFF_res(4)
              ! now update NEE
              NEE(n)=NEE(n)+FLUXES(n,17)
          else
              FLUXES(n,17) = 0.
              CFF = 0. ; NCFF = 0.
              CFF_res = 0. ; NCFF_res = 0.
          end if
          ! update all pools this time
          POOLS(n+1,1) = max(0., POOLS(n+1,1) - CFF(1) - NCFF(1) )
          POOLS(n+1,2) = max(0., POOLS(n+1,2) - CFF(2) - NCFF(2) )
          POOLS(n+1,3) = max(0., POOLS(n+1,3) - CFF(3) - NCFF(3) )
          POOLS(n+1,4) = max(0., POOLS(n+1,4) - CFF(4) - NCFF(4) )
          POOLS(n+1,5) = max(0., POOLS(n+1,5) + (labile_residue+foliar_residue+roots_residue) + (NCFF(1)+NCFF(2)+NCFF(3))-CFF(5)-NCFF(5) )
          POOLS(n+1,6) = max(0., POOLS(n+1,6) + (wood_residue-soil_loss_with_roots) + (NCFF(4)+NCFF(5)))
          ! some variable needed for the EDCs
          ! reallocation fluxes for the residues
!          disturbance_residue_to_litter(n) =
!          (labile_residue+foliar_residue+roots_residue) & 
!                                           + (NCFF(1)+NCFF(2)+NCFF(3))
!          disturbance_loss_from_litter(n)  = CFF(5)+NCFF(5)
!          disturbance_residue_to_som(n)    = NCFF(4)+NCFF(5)
!          disturbance_loss_from_som(n)     = soil_loss_with_roots
          ! this is intended for use with the R interface for subsequent post
          ! processing
          if (allocated(extracted_C)) then
            ! harvested carbon from all pools
            extracted_C(n) = (wood_loss-(wood_residue+CFF_res(4)+NCFF_res(4))) &
                           + (labile_loss-(labile_residue+CFF_res(1)+NCFF_res(1))) &
                           + (foliar_loss-(foliar_residue+CFF_res(2)+NCFF_res(2))) &
                           + (roots_loss-(roots_residue+CFF_res(3)+NCFF_res(3)))
          endif ! allocated extracted_C
          ! total carbon loss from the system
          C_total = (labile_residue+foliar_residue+roots_residue+wood_residue+sum(NCFF)) &
                  - (labile_loss+foliar_loss+roots_loss+wood_loss+soil_loss_with_roots+sum(CFF))

          ! if total clearance occured then we need to ensure some minimum
          ! values and reforestation is assumed one year forward
          if (met(8,n) > 0.99) then
              m=0 ; test=sum(deltat(n:(n+m)))
              ! FC Forest Statistics 2015 lag between harvest and restocking ~ 2 year
              restocking_lag = 365*2
              do while (test < restocking_lag)
                 m=m+1 ; test = sum(deltat(n:(n+m)))
                 !  get out clause for hitting the end of the simulation
                 if (m+n >= nodays) test = restocking_lag
              enddo
              reforest_day = min((n+m), nodays)
          endif ! if total clearance

      endif ! end deforestation info

      ! 
      ! then deal with fire
      ! 

      if (met(9,n) > 0.) then

         !/*first fluxes*/
         !/*LABILE*/
         CFF(1) = POOLS(n+1,1)*met(9,n)*combust_eff(1)
         NCFF(1) = POOLS(n+1,1)*met(9,n)*(1-combust_eff(1))*(1-rfac)
         !/*foliar*/
         CFF(2) = POOLS(n+1,2)*met(9,n)*combust_eff(2)
         NCFF(2) = POOLS(n+1,2)*met(9,n)*(1-combust_eff(2))*(1-rfac)
         !/*root*/
         CFF(3) = 0. ! POOLS(n+1,3)*met(9,n)*combust_eff(3)
         NCFF(3) = 0. ! POOLS(n+1,3)*met(9,n)*(1-combust_eff(3))*(1-rfac)
         !/*wood*/
         CFF(4) = POOLS(n+1,4)*met(9,n)*combust_eff(4)
         NCFF(4) = POOLS(n+1,4)*met(9,n)*(1-combust_eff(4))*(1-rfac)
         !/*litter*/
         CFF(5) = POOLS(n+1,5)*met(9,n)*combust_eff(5)
         NCFF(5) = POOLS(n+1,5)*met(9,n)*(1-combust_eff(5))*(1-rfac)
         !/*fires as daily averages to comply with units*/
         FLUXES(n,17)=(CFF(1)+CFF(2)+CFF(3)+CFF(4)+CFF(5))/deltat(n)

         !/*all fluxes are at a daily timestep*/
         NEE(n)=NEE(n)+FLUXES(n,17)

         !// update pools
         !/*Adding all fire pool transfers here*/
         POOLS(n+1,1)=POOLS(n+1,1)-CFF(1)-NCFF(1)
         POOLS(n+1,2)=POOLS(n+1,2)-CFF(2)-NCFF(2)
         POOLS(n+1,3)=POOLS(n+1,3)-CFF(3)-NCFF(3)
         POOLS(n+1,4)=POOLS(n+1,4)-CFF(4)-NCFF(4)
         POOLS(n+1,5)=POOLS(n+1,5)-CFF(5)-NCFF(5)+NCFF(1)+NCFF(2)+NCFF(3)
         POOLS(n+1,6)=POOLS(n+1,6)+NCFF(4)+NCFF(5)
         ! some variable needed for the EDCs
         ! reallocation fluxes for the residues
         disturbance_residue_to_litter(n) = (NCFF(1)+NCFF(2)+NCFF(3))
         disturbance_residue_to_som(n) = (NCFF(4)+NCFF(5))

      endif ! end burnst area issues

!    do nxp = 1, nopools
!       if (POOLS(n+1,nxp) /= POOLS(n+1,nxp)) then
!          print*,"step",n, nxp
!          print*,"met",met(:,n)
!          print*,"POOLS",POOLS(n,:)
!          print*,"FLUXES",FLUXES(n,:)
!          stop
!       endif
!    enddo

    end do ! nodays loop

!   print*,sum(FLUXES(1:nodays,11)+FLUXES(1:nodays,15)&
!          +disturbance_residue_to_som) / sum(FLUXES(1:nodays,14))

!   print*,sum(FLUXES(1:nodays,7)) / sum(FLUXES(1:nodays,11))
  end subroutine CARBON_MODEL
  !
  !------------------------------------------------------------------
  !
  double precision function acm(drivers,constants)

    ! the Aggregated Canopy Model, is a Gross Primary Productivity (i.e.
    ! Photosyntheis) emulator which operates at a daily time step. ACM can be
    ! paramaterised to provide reasonable results for most ecosystems.

    implicit none

    ! declare input variables
    double precision, intent(in) :: drivers(12) & ! acm input requirements
                         ,constants(10) ! ACM parameters

    ! declare local variables
    double precision :: gc, pn, pd, pp, qq, ci, e0, dayl, cps, dec, nit &
             ,trange, sinld, cosld,aob, mult &
             ,mint,maxt,radiation,co2,lai,doy,lat &
             ,deltaWP,Rtot,NUE,temp_exponent,dayl_coef &
             ,dayl_const,hydraulic_exponent,hydraulic_temp_coef &
             ,co2_comp_point,co2_half_sat,lai_coef,lai_const

    ! initial values
    gc=0.0 ; pp=0.0 ; qq=0.0 ; ci=0.0 ; e0=0.0 ; dayl=0.0 ; cps=0.0 ; dec=0.0 ; nit=1.0

    ! load driver values to correct local vars
    lai = drivers(1)
    maxt = drivers(2)
    mint = drivers(3)
    nit = drivers(4)   
    co2 = drivers(5)
    doy = drivers(6)
    radiation = drivers(8)
    lat = drivers(7)

    ! load parameters into correct local vars
!    pi = drivers(11)
    deltaWP = drivers(9)
    Rtot = drivers(10)
    NUE = constants(1)
    dayl_coef = constants(2)
    co2_comp_point = constants(3) 
    co2_half_sat = constants(4)
    dayl_const = constants(5)
    hydraulic_temp_coef = constants(6)
    lai_coef = constants(7)
    temp_exponent = constants(8)
    lai_const = constants(9)
    hydraulic_exponent = constants(10)

    ! determine temperature range 
    trange=0.5*(maxt-mint)
    ! daily canopy conductance, of CO2 or H2O? 
    gc=abs(deltaWP)**(hydraulic_exponent)/((hydraulic_temp_coef*Rtot+trange))
    ! maximum rate of temperature and nitrogen (canopy efficiency) limited photosynthesis (gC.m-2.day-1)
    pn=lai*nit*NUE*exp(temp_exponent*maxt)
    ! pp and qq represent limitation by diffusion and metabolites respecitively
    pp=pn/gc ; qq=co2_comp_point-co2_half_sat
    ! calculate internal CO2 concentration (ppm)
!    ci=0.5*(co2+qq-pp+sqrt((co2+qq-pp)**2-4.0*(co2*qq-pp*co2_comp_point)))
    ci=0.5*(co2+qq-pp+sqrt(((co2+qq-pp)*(co2+qq-pp))-4.0*(co2*qq-pp*co2_comp_point)))
    ! limit maximum quantium efficiency by leaf area, hyperbola
    e0=lai_coef*(lai*lai)/((lai*lai)+lai_const)
    ! calculate day length (hours)
    dec = - asin( sin( 23.45 * deg_to_rad ) * cos( 2.0 * pi * ( doy + 10.0 ) /365.0 ) )
    sinld = sin( lat*deg_to_rad ) * sin( dec )
    cosld = cos( lat*deg_to_rad ) * cos( dec )
    aob = max(-1.0,min(1.0,sinld / cosld))
    dayl = 12.0 * ( 1.0 + 2.0 * asin( aob ) / pi )

!--------------------------------------------------------------
!    ! calculate day length (hours - not really hours)
!    ! This is the old REFLEX project calculation but it is wrong so anyway here
!    ! we go...
!    dec=-23.4*cos((360.0*(doy+10.0)/365.0)*pi/180.0)*pi/180.0
!    mult=tan(lat*pi/180.0)*tan(dec)
!    if (mult>=1.0) then
!      dayl=24.0
!    else if (mult<=-1.0) then
!      dayl=0.0
!    else
!      dayl=24.0*acos(-mult)/pi
!    end if
! ---------------------------------------------------------------
    ! calculate CO2 limited rate of photosynthesis
    pd=gc*(co2-ci)
    ! calculate combined light and CO2 limited photosynthesis
    cps=e0*radiation*pd/(e0*radiation+pd)
    ! correct for day length variation
    acm=cps*(dayl_coef*dayl+dayl_const)

    ! don't forget to return
    return

  end function acm
  !
  !------------------------------------------------------------------
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
    sumsq_x = sum(x*x)
    ! calculate the sum of the product of xy
    sum_product_xy = sum(x*y)
    ! calculate the gradient
    linear_model_gradient = ( (interval*sum_product_xy) - (sum_x*sum_y) ) &
                          / ( (interval*sumsq_x) - (sum_x*sum_x) )

    ! for future reference here is how to calculate the intercept
!    intercept = ( (sum_y*sumsq_x) - (sum_x*sum_product_xy) ) &
!              / ( (interval*sumsq_x) - (sum_x*sum_x) )

    ! don't forget to return to the user
    return

  end function linear_model_gradient
  !
  !------------------------------------------------------------------
  !
!
!--------------------------------------------------------------------
!
end module CARBON_MODEl_MOD
