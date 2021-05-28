
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
  ! See function / subroutine specific comments for exceptions and contributors
  !!!!!!!!!!!

! make all private
private

! explicit publics
public :: CARBON_MODEL     &
         ,soil_frac_clay   &
         ,soil_frac_sand   &
         ,nos_soil_layers  &
         ,disturbance_residue_to_litter &
         ,disturbance_residue_to_som    &
         ,disturbance_residue_to_cwd    &
         ,DON_leaching     &
         ,flab_to_ra       &
         ,wlab_to_ra       &
         ,rlab_to_ra       &
         ,npp_to_roots     &
         ,npp_to_wood      &
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

!
! double precision
!

! GPP random forest emulator
double precision, allocatable, dimension(:,:) ::     leftDaughter, & ! left daughter for forest
                                                    rightDaughter, & ! right daughter for forets
                                                       nodestatus, & ! nodestatus for forests
                                                       xbestsplit, & ! for forest
                                                         nodepred, & ! prediction value for each tree
                                                          bestvar    ! for randomForests

double precision, allocatable, dimension(:) :: tmp_x,tmp_m,itemp,ivpd,iphoto & ! GSI model
                                              ,disturbance_residue_to_litter & ! Disturbance variables
                                              ,disturbance_residue_to_som    &
                                              ,disturbance_residue_to_cwd    &
                                              ,extracted_C
! Variables related to efficiency savings
double precision, allocatable, dimension(:) :: deltat_1,avg_temp
! Nitrogen related variables
double precision, allocatable, dimension(:) :: flab_to_ra,wlab_to_ra,rlab_to_ra &
                                              ,npp_to_roots,npp_to_wood

double precision :: Ndemand_leaf,Ndemand_root,Ndemand_wood, &
                    foliarN_retrans,foliarN_loss,DON_leaching

! ACM related parameters
double precision, parameter :: pi = 3.1415927
double precision, parameter :: deg_to_rad = pi/180.0

!
! integer
!

integer :: gsi_lag   & ! GSI model
          ,dim_1     & ! (RandomForest) dimension 1 of response surface
          ,dim_2     & ! dimension 2 of response surface
          ,nos_trees & ! number of trees in randomForest
          ,nos_inputs  ! number of driver inputs

! for consisteny between requirements of different models
integer, parameter :: nos_root_layers = 2, nos_soil_layers = nos_root_layers + 1
double precision, dimension(nos_soil_layers) :: soil_frac_clay,soil_frac_sand

contains
!
!--------------------------------------------------------------------
!
  subroutine CARBON_MODEL(start,finish,met,pars,deltat,nodays,lat,lai,NEE,FLUXES,POOLS &
                         ,nopars,nomet,nopools,nofluxes,GPP)

    ! The Data Assimilation Linked Ecosystem Carbon Nitrogen - Growing Season
    ! Index - multiple LABILE pool - Forest Rotation (DALECN_GSI_DFOL_LABILE_FR) model.
    ! The subroutine calls the Aggregated Canopy Model to simulate GPP and
    ! partitions between various ecosystem carbon pools. These pools are
    ! subject to turnovers / decompostion resulting in ecosystem phenology and fluxes of CO2.
    ! Based on simplfied aconite code providing by Quinn Thomas a simplified
    ! nitrogen cycle has been constructed to test the hypothesis of greater model complexity
    ! the improve constraint

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
    double precision :: CFF(8) = 0, CFF_res(8) = 0    & ! combusted and non-combustion fluxes
                       ,NCFF(8) = 0, NCFF_res(8) = 0  & ! with residue and non-residue seperates
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
    double precision :: labile_loss,foliar_loss              &
                       ,roots_loss,wood_loss                 &
                       ,labile_residue,foliar_residue        &
                       ,labile_wood_loss,labile_wood_residue &
                       ,labile_root_loss,labile_root_residue &
                       ,roots_residue,wood_residue           &
                       ,wood_pellets,C_total                 &
                       ,labile_frac_res                      &
                       ,Cstem,Cbranch,Crootcr                &
                       ,stem_residue,branch_residue          &
                       ,coarse_root_residue                  &
                       ,soil_loss_with_roots
    integer :: reforest_day, harvest_management,restocking_lag

    ! nitrogen cycle related parameters
    double precision :: avail_DIN,N_immobilised,NC_foliar,NC_wood,NC_root,NC_som &
                       ,DIN_leaching,gross_nmin,foliarN_retrans_parameter        &
                       ,labile_storage_overload,labile_ratios(3),N_deposition    &
                       ,CN_foliar,NC_litter,avail_labfol,avail_labwood,avail_labroot &
                       ,decomp_reduction_ratio

    ! local variables for labile Cwood release
    double precision :: Cwood_labile_release_gradient &
                       ,Cwood_labile_half_saturation  &
                       ,Croot_labile_release_gradient &
                       ,Croot_labile_half_saturation

    double precision, dimension(nodays) :: Croot_labile_release_coef &
                                          ,Cwood_labile_release_coef

    ! local variables for GSI phenology model
    double precision :: Tfac,Photofac,VPDfac & ! oC, seconds, Pa
                       ,Tfac_demon,Photofac_demon,VPDfac_demon &
                       ,delta_gsi,tmp,gradient &
                       ,fol_turn_crit,lab_turn_crit &
                       ,gsi_history(22),just_grown

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

    ! POOLS are:
    ! 1  = labile foliar C (p18)
    ! 2  = foliar          (p19)
    ! 3  = root            (p20)
    ! 4  = wood            (p21)
    ! 5  = litter          (p22)
    ! 6  = som             (p23)
    ! 7  = labile root C   (p39)
    ! 8  = labile wood C   (p40)
    ! 9  = CWD             (p44)
    ! 10 = labile foliar N (p54)
    ! 11 = litter N        (p55)
    ! 12 = DIN             (p56)

    ! p(30) = labile foliar replanting
    ! p(31) = foliar replanting
    ! p(32) = fine root replanting
    ! p(33) = wood replanting
    ! p(41) = labile root replanting
    ! p(42) = labile wood replanting
    ! p(57) = labile foliar N replanting

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
    ! 19 = cwd decomp to litter
    ! 20 = litter N turnover
    ! 21 = N immobilisation / mineralisation between DON <> DIN

    ! PARAMETERS
    ! 17+4(GSI) values

    ! p(1) Litter to SOM conversion rate
    ! p(2) Fraction of GPP respired (Base rate)
    ! p(3) GSI leaf fall senstivity parameter
    ! p(4) Fraction of NPP allocated to roots
    ! p(5) max leaf turnover fraction
    ! p(6) Turnover rate of wood (day-1)
    ! p(7) Turnover rate of roots (day-1)
    ! p(8) Litter turnover rate (day-1)
    ! p(9) SOM turnover rate (day-1)
    ! p(10) Parameter in exponential term of temperature - \theta
    ! p(11) Foliar mean nitrogen content (gN.m-2 leaf area)
    ! p(12) = max labile turnover
    ! p(13) = Fraction of NPP allocated to Clab
    ! p(14) = min temp threshold
    ! p(15) = max temp threshold
    ! p(16) = min photoperiod threshold (GIS)
    ! p(17) = LMA
    ! p(24) = max photoperiod threshold (GSI)
    ! p(25) = min VPD threshold (GSI)
    ! p(26) = max VPD threshold
    ! p(27) = minimum fraction by which GPP should increase for a given LAI increment to be allowed
    ! p(28) = fraction of Cwood which is Cbranch
    ! p(29) = fraction of Cwood which is Ccoarseroot
    ! p(34) = GSI leaf growth sensitivity parameter

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
    gpppars(4) = 10.0**pars(11) ! foliar N
    gpppars(7) = lat
    gpppars(9) = -2.056214 ! -2.0 ! leafWP-soilWP
    gpppars(10) = 1.274201 ! 1.0 ! totaly hydraulic resistance
    gpppars(11) = pi

    ! assign acm parameters
    constants(1)=pars(52)**(1/(-2.999299929993))
    constants(2)=1.827248e-2
    constants(3)=50.93006
    constants(4)=283.6717
    constants(5)=0.0942341
    constants(6)=4.157976e-3
    constants(7)=7.949699
    constants(8)=0.0157717
    constants(9)=3.707546
    constants(10)=1.612526e-8

    ! initial values for deforestation variables
    labile_loss = 0.0      ; foliar_loss = 0.0
    roots_loss = 0.0       ; wood_loss = 0.0
    labile_residue = 0.0   ; foliar_residue = 0.0
    roots_residue = 0.0    ; wood_residue = 0.0
    stem_residue = 0.0     ; branch_residue = 0.0
    labile_root_loss = 0.0 ; labile_root_residue = 0.0
    reforest_day = 0
    soil_loss_with_roots = 0.0
    coarse_root_residue = 0.0
    post_harvest_burn = 0.0

    ! now load the hardcoded forest management parameters into their locations

    ! Parameter values for deforestation variables
    ! scenario 1
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(1) = 1.0
    roots_frac_res(1)   = 1.0
    rootcr_frac_res(1) = 1.0
    branch_frac_res(1) = 1.0
    stem_frac_res(1)   = 0.0 !
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
    post_harvest_burn(1) = 1.0

    !## scen 2
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(2) = 1.0
    roots_frac_res(2)   = 1.0
    rootcr_frac_res(2) = 1.0
    branch_frac_res(2) = 1.0
    stem_frac_res(2)   = 0.0 !
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
    post_harvest_burn(2) = 0.0

    !## scen 3
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(3) = 0.5
    roots_frac_res(3)   = 1.0
    rootcr_frac_res(3) = 1.0
    branch_frac_res(3) = 0.0
    stem_frac_res(3)   = 0.0 !
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
    post_harvest_burn(3) = 0.0

    !## scen 4
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(4) = 0.5
    roots_frac_res(4)   = 1.0
    rootcr_frac_res(4) = 0.0
    branch_frac_res(4) = 0.0
    stem_frac_res(4)   = 0.0 !
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
    post_harvest_burn(4) = 0.0

    ! for the moment override all paritioning parameters with those coming from
    ! CARDAMOM
    Cbranch_part = pars(28)
    Crootcr_part = pars(29)

    if (start == 1) then

    ! assigning initial conditions
    POOLS(1,1) =pars(18) ! labile fol
    POOLS(1,2) =pars(19) ! fol
    POOLS(1,3) =pars(20) ! root
    POOLS(1,4) =pars(21) ! wood
    POOLS(1,5) =pars(22) ! litter
    POOLS(1,6) =pars(23) ! som
    POOLS(1,7) =pars(39) ! labile root
    POOLS(1,8) =pars(40) ! labile wood
    POOLS(1,9) =pars(44) ! CWD pool
    POOLS(1,10)=pars(54) ! flabN
    POOLS(1,11)=pars(55) ! litN
    POOLS(1,12)=pars(56) ! DIN

    if (.not.allocated(disturbance_residue_to_som)) then
        allocate(disturbance_residue_to_litter(nodays), &
                 disturbance_residue_to_som(nodays), &
                 disturbance_residue_to_cwd(nodays))
    endif
    disturbance_residue_to_litter = 0.0
    disturbance_residue_to_som = 0.0
    disturbance_residue_to_cwd = 0.0

    ! declare fire constants (labile, foliar, roots, wood, litter)
    combust_eff(1) = 0.1 ; combust_eff(2) = 0.9
    combust_eff(3) = 0.1 ; combust_eff(4) = 0.5
    combust_eff(5) = 0.3 ; rfac = 0.5

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
        gsi_lag=max(2,maxval(nint(tmp_m)))
    end if ! .not.allocated(tmp_x)
    ! assign our starting value
    gsi_history = pars(36)-1d0
    just_grown=pars(35)
    endif

    ! assign climate sensitivities
    fol_turn_crit=pars(34)-1d0
    lab_turn_crit=pars(3)-1d0

    ! calculate the inverse of the GSI demoninators (efficiency saving)
    Tfac_demon=1/(pars(15)-pars(14))
    Photofac_demon=1/(pars(24)-pars(16))
    VPDfac_demon=1/(pars(26)-pars(25))

    ! calculate inverse timing information if it does not already exist
    ! and mean temperature as this is invarient between iterations
    if (.not.allocated(deltat_1)) then
       allocate(deltat_1(nodays)) ; deltat_1=1/deltat
       allocate(avg_temp(nodays)) ; avg_temp = 0.5*(met(3,1:nodays)+met(2,1:nodays))
    endif

    if (start == 1) then
    ! Wood/fine root labile turnover parameters
    ! parmeters generated on the assumption of 5 % / 95 % activation at key
    ! temperature values. Roots 1oC/30oC, wood 5oC/30oC.
    Croot_labile_release_gradient = 0.2995754 !0.2998069
    Croot_labile_half_saturation = 17.49752 !15.28207
    Cwood_labile_release_gradient = 0.2995754 ! 0.25
    Cwood_labile_half_saturation = 17.49752   ! 20.0

    ! what is temperature limitation on root labile release
    Croot_labile_release_coef = 1 / &
                               (1.0+exp(-Croot_labile_release_gradient* &
                               (avg_temp-Croot_labile_half_saturation)))

    ! what is temperature limitation on wood labile release
    Cwood_labile_release_coef = 1 / &
                               (1.0+exp(-Cwood_labile_release_gradient* &
                               (avg_temp-Cwood_labile_half_saturation)))

    endif ! start == 1
    ! Nitrogen cycle parameters
    CN_foliar = pars(17)/gpppars(4)
    NC_foliar = 1/CN_foliar
    NC_root = 1/pars(46)
    NC_wood = 1/pars(47)
    NC_som = 1/pars(48)
    DON_leaching=pars(49)
    N_deposition=pars(51)
    if (.not.(allocated(flab_to_ra))) then
       allocate(flab_to_ra(nodays),rlab_to_ra(nodays),wlab_to_ra(nodays) &
               ,npp_to_roots(nodays),npp_to_wood(nodays))
    endif
    flab_to_ra = 0.0 ; rlab_to_ra = 0.0 ; wlab_to_ra = 0.0
    npp_to_roots = 0.0 ; npp_to_wood = 0.0

    !
    ! Begin looping through each time step
    !

    do n = start, finish

      !!!
      ! Plant structure
      !!!

      ! calculate LAI value
      lai(n)=POOLS(n,2)/pars(17)

      !!!
      ! Phenology
      !!!

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
      Tfac = (met(10,n)-(pars(14)-273.15)) * Tfac_demon
      Tfac = min(1d0,max(0d0,Tfac))
      ! photoperiod limitation
      Photofac = (met(11,n)-pars(16)) * Photofac_demon
      Photofac = min(1d0,max(0d0,Photofac))
      ! VPD limitation
      VPDfac = 1d0 - ( (met(12,n)-pars(25)) * VPDfac_demon )
      VPDfac = min(1d0,max(0d0,VPDfac))

      ! these allocated if post-processing
      if (allocated(itemp)) then
         itemp(n) = Tfac ; ivpd(n) = VPDfac ; iphoto(n) = Photofac
      endif

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

      ! first assume that nothing is happening
      FLUXES(n,9) = 0.0  ! leaf turnover
      FLUXES(n,16) = 0.0 ! leaf growth

      if (gradient < fol_turn_crit .or. FLUXES(n,18) == 0.0) then
         ! we are in a decending condition so foliar turnover
         FLUXES(n,9) = pars(5)*(1.0-FLUXES(n,18))
         just_grown = 0.5
      end if

      !!!
      ! Tissue turnover first
      !!!

      ! total leaf litter production
      FLUXES(n,10) = POOLS(n,2)*(1.0-(1.0-FLUXES(n,9))**deltat(n))*deltat_1(n)
      ! foliar N retranslocation
      foliarN_retrans = FLUXES(n,10)*pars(53)*NC_foliar
      foliarN_loss    = FLUXES(n,10)*(1.0-pars(53))*NC_foliar
      ! total CWD production
      FLUXES(n,11) = POOLS(n,4)*(1.0-(1.0-pars(6))**deltat(n))*deltat_1(n)
!      ! total root litter production
      FLUXES(n,12) = POOLS(n,3)*(1.0-(1.0-pars(7))**deltat(n))*deltat_1(n)

      !!!
      ! Litter / CWD / SOM turnover
      !!!

      ! calculate NC_litter
      NC_litter = POOLS(n,11)/POOLS(n,5)

      ! temprate (i.e. temperature modified rate of metabolic activity))
      FLUXES(n,2) = exp(pars(10)*avg_temp(n))
      ! total litter turnover
      FLUXES(n,13) = POOLS(n,5)*(1.0-(1.0-FLUXES(n,2)*pars(8))**deltat(n))*deltat_1(n)
      ! litter to som (gC.m-2.day-1)
      FLUXES(n,15) = FLUXES(n,13)*(1.0-pars(1))
      ! respiration heterotrophic litter (gC.m-2.day-1)
      FLUXES(n,13) = FLUXES(n,13)*pars(1)
      ! cwd turnover to the litter pool (gC.m-2.day-1)
      FLUXES(n,19) = POOLS(n,9)*(1.0-(1.0-FLUXES(n,2)*pars(43))**deltat(n))*deltat_1(n)
      ! litter nitrogen turnover to som (gN.m-2.day-1)
      FLUXES(n,20) = (FLUXES(n,13)+FLUXES(n,15)) * NC_litter
      ! determine how much nitrogen will be immobilised to maintain SOM C:N
      ! if the value here is negative that means more nitrogen is needed from the dissolved inorganic N pool
      N_immobilised = (FLUXES(n,15)*NC_som) - FLUXES(n,20)

      ! respiration heterotrophic som (gC.m-2.day-1)
      FLUXES(n,14) = POOLS(n,6)*(1.0-(1.0-FLUXES(n,2)*pars(9))**deltat(n))*deltat_1(n)

      ! determine how much DIN is available for plants
      avail_DIN = POOLS(n,12)
      ! is there more DIN than is needed for the immobilsation?
      if (avail_DIN >= (N_immobilised*deltat(n))) then
         ! if we have enough for demand then we substract it from the available
         avail_DIN = max(0d0,avail_DIN - (N_immobilised*deltat(n)))
         FLUXES(n,21) = N_immobilised
      else
         ! in which case we have no available N and all DIN has itself been immobilsed
         FLUXES(n,21) = max(0d0,avail_DIN*deltat_1(n))
         avail_DIN = 0.0
      endif

      ! now if the supply from DIN is less then we need to adjust the turnovers / respiration to account for the N status
      if (FLUXES(n,21) < N_immobilised) then
         decomp_reduction_ratio = FLUXES(n,21)/N_immobilised
         ! respiration heterotrophic litter (gC.m-2.day-1)
         FLUXES(n,13) = FLUXES(n,13)*decomp_reduction_ratio
         ! litter to som (gC.m-2.day-1)
         FLUXES(n,15) = FLUXES(n,15)*decomp_reduction_ratio
         ! litter nitrogen turnover to som (gN.m-2.day-1)
         FLUXES(n,20) = (FLUXES(n,13)+FLUXES(n,15)) * NC_litter
      endif

      ! finally what is the potential N (gC.m-2.day-1) available for mineralisation
      ! i.e. soil turnover / CN
      gross_nmin = FLUXES(n,14)*NC_som

      !!!
      ! Photosynthesis and allocation
      !!!

      ! load next met / lai values for ACM
      gpppars(1)=lai(n)
      gpppars(2)=met(3,n) ! max temp
      gpppars(3)=met(2,n) ! min temp
      gpppars(5)=met(5,n) ! co2
      gpppars(6)=ceiling(met(6,n)-(deltat(n)*0.5)) ! doy
      gpppars(8)=met(4,n) ! radiation

      ! GPP (gC.m-2.day-1)
      if (lai(n) > 1e-8) then
         FLUXES(n,1) = acm(gpppars,constants)
      else
         FLUXES(n,1) = 0.0
      endif
      ! autotrophic respiration (gC.m-2.day-1)
      FLUXES(n,3) = pars(2)*FLUXES(n,1)
      ! labile production (gC.m-2.day-1)
      FLUXES(n,5) = (FLUXES(n,1)-FLUXES(n,3))*pars(13)
      ! disused
      FLUXES(n,4) = 0.0

      ! foliar allocation (gC.m-2.day-1)
      if (gradient > lab_turn_crit) then
         ! potential leaf production (gC.m-2.day-1)
         FLUXES(n,16) = pars(12)*FLUXES(n,18)
         just_grown = 1.5
         ! check carbon return
         tmp = POOLS(n,1)*(1.0-(1.0-FLUXES(n,16))**deltat(n))*deltat_1(n)
         tmp = (POOLS(n,2)+tmp)/pars(17)
         gpppars(1)=tmp
         tmp = acm(gpppars,constants)
         ! determine if increase in LAI leads to an improvement in GPP greater
         ! than
         ! critical value, if not then no labile turnover allowed
         if ( ((tmp - FLUXES(n,1))/FLUXES(n,1)) < pars(27) ) then
             FLUXES(n,16) = 0.0
         endif
      else if (gradient > fol_turn_crit .and. gradient < lab_turn_crit) then
         ! probaly we want nothing to happen, however if we are at the seasonal
         ! maximum we will consider further growth still
         if (just_grown >= 1.0) then
            ! we are between so definitely not losing foliage and we have
            ! previously been growing so maybe we still have a marginal return on
            ! doing so again
            FLUXES(n,16) = pars(12)*FLUXES(n,18)
            ! but possibly gaining some?
            ! determine if this is a good idea based on GPP increment
            tmp = POOLS(n,1)*(1.0-(1.0-FLUXES(n,16))**deltat(n))*deltat_1(n)
            tmp = (POOLS(n,2)+tmp)/pars(17)
            gpppars(1)=tmp
            tmp = acm(gpppars,constants)
            ! determine if increase in LAI leads to an improvement in GPP greater
            ! than
            ! critical value, if not then no labile turnover allowed
            if ( ((tmp - FLUXES(n,1))/FLUXES(n,1)) < pars(27) ) then
                FLUXES(n,16) = 0.0
            endif
         end if ! Just grown?
      endif ! gradient choice

      ! allocation to fine root labile (gC.m-2.day-1)
      npp_to_roots(n) = (FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,5))*pars(4)
      ! labile allocation to wood labile (gC.m-2.day-1)
      npp_to_wood(n)  =  FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,5)-npp_to_roots(n)

      ! fine root production (currently assuming temperature only limitation)
      FLUXES(n,6) = pars(37)*Croot_labile_release_coef(n)

      ! then impose secondary moisture stress.
      ! Assuming that if moisture stress on canopy is severe that wood
      ! growth is also limited. In actual fact it is probably much worse
      ! wood production
      FLUXES(n,7) = pars(38)*Cwood_labile_release_coef(n)*VPDfac

      !
      ! Update temperature dependent labile turnover rates
      !

      ! fine root labile pool allocation (gC.m-2.day-1)
      FLUXES(n,6) = POOLS(n,7)*(1.0-(1.0-FLUXES(n,6))**deltat(n))*deltat_1(n)
      ! wood labile pool allocation (gC.m-2.day-1)
      FLUXES(n,7) = POOLS(n,8)*(1.0-(1.0-FLUXES(n,7))**deltat(n))*deltat_1(n)
      ! total labile release to foliage (gC.m-2.day-1)
      FLUXES(n,8) = POOLS(n,1)*(1.0-(1.0-FLUXES(n,16))**deltat(n))*deltat_1(n)

      ! determine nitrogen requirements (gN.m-2.day-1)
      Ndemand_leaf = (FLUXES(n,8)*NC_foliar*deltat(n)) - POOLS(n,10)
      Ndemand_leaf = max(0d0,Ndemand_leaf)
      Ndemand_root = (FLUXES(n,6)*NC_root*deltat(n))
      Ndemand_root = max(0d0,Ndemand_root)
      Ndemand_wood = (FLUXES(n,7)*NC_wood*deltat(n))
      Ndemand_wood = max(0d0,Ndemand_wood)

      ! if demand is greater than current supply available from soil
      if (avail_DIN <= 0.0) then

         ! if there is no nitrogen available then everything set to zero, except
         ! any N allocation from the foliar labile N pool
         FLUXES(n,8)  = FLUXES(n,8)-(Ndemand_leaf*CN_foliar*deltat_1(n))
         FLUXES(n,6)  = 0.0 ; FLUXES(n,7)  = 0.0 ; avail_DIN = 0.0
         Ndemand_leaf = 0.0 ; Ndemand_root = 0.0 ; Ndemand_wood = 0.0

       else if ((Ndemand_leaf + Ndemand_root + Ndemand_wood) > avail_DIN) then

         ! we must grow less, while being careful to recognise that we have two
         ! different demand pathways
         ! preference N allocation to foliage
         if (Ndemand_leaf > 0.0) then
            ! labile directed to Ra; applied only to the portion of the flux
            ! linked to the additional Ndemand_leaf
            flab_to_ra(n) = (Ndemand_leaf*CN_foliar*deltat_1(n))*(1.0-min(1d0,avail_DIN/Ndemand_leaf))
            ! leaf growth from labile
            FLUXES(n,8) = max(0d0,FLUXES(n,8)-flab_to_ra(n))
            ! recalculate the new Ndemand_leaf which for the DIN mass balance
            Ndemand_leaf = (FLUXES(n,8)*NC_foliar*deltat(n)) - POOLS(n,10)
            Ndemand_leaf = max(0d0,Ndemand_leaf)
            ! have updated the foliar N demand we update the available DIN
            avail_DIN = max(0d0,avail_DIN - Ndemand_leaf)
         endif ! Ndemand_leaf

         ! then to roots followed by wood
         if (Ndemand_root > 0.0) then
            ! labile directed to Ra
            rlab_to_ra(n) = FLUXES(n,6)*(1.0-min(1d0,avail_DIN/Ndemand_root))
            ! root growth
            FLUXES(n,6) = max(0d0,FLUXES(n,6)-rlab_to_ra(n))
            ! recalculate the new Ndemand_root which for the DIN mass balance
            Ndemand_root = (FLUXES(n,6)*NC_root)*deltat(n)
            Ndemand_root = max(0d0,Ndemand_root)
            ! have updated the root N demand we update the available DIN
            avail_DIN = max(0d0,avail_DIN - Ndemand_root)
         endif ! Ndemand_root

         ! then to wood
         if (Ndemand_wood > 0.0) then
            ! labile directed to Ra
            wlab_to_ra(n) = FLUXES(n,7)*(1.0-min(1d0,avail_DIN/Ndemand_wood))
            ! wood growth
            FLUXES(n,7) = max(0d0,FLUXES(n,7)-wlab_to_ra(n))
            ! recalculate the new Ndemand_wood which for the DIN mass balance
            Ndemand_wood = (FLUXES(n,7)*NC_wood)*deltat(n)
            Ndemand_wood = max(0d0,Ndemand_wood)
            ! have updated the wood N demand we update the available DIN
            avail_DIN = max(0d0,avail_DIN - Ndemand_wood)
         endif ! Ndemand_wood

      endif ! total N demand exceeded

      ! Finally we assess whether we need to burn off the allocated labile carbon because we have exceeded maxmimum storage.
      ! Not forgetting to account for the extraction already determined
!      avail_labfol  = POOLS(n,1) - ((FLUXES(n,8)+flab_to_ra(n)) * deltat(n))
!      avail_labroot = POOLS(n,7) - ((FLUXES(n,6)+rlab_to_ra(n)) * deltat(n))
!      avail_labwood = POOLS(n,8) - ((FLUXES(n,7)+wlab_to_ra(n)) * deltat(n))
      avail_labfol  = POOLS(n,1) - (FLUXES(n,8) * deltat(n))
      avail_labroot = POOLS(n,7) - (FLUXES(n,6) * deltat(n))
      avail_labwood = POOLS(n,8) - (FLUXES(n,7) * deltat(n))
      labile_storage_overload = (avail_labfol+avail_labwood+avail_labroot) - (POOLS(n,4)*pars(45))
      if (labile_storage_overload > 0.0) then
         ! assume that we burn off excess carbon in order of wood, root,
         ! foliage then a bit from all
!         if ((wlab_to_ra(n)*deltat(n)) >= labile_storage_overload) then
!             ! therefore we spend the amount of the wood labile only to balance the books
!             wlab_to_ra(n) = labile_storage_overload*deltat_1(n)
!             ! leaving the remained in storage
!             rlab_to_ra(n) = 0.0 ; flab_to_ra(n) = 0.0
!         else if (((wlab_to_ra(n)+rlab_to_ra(n))*deltat(n)) >= labile_storage_overload) then
!             ! therefore all wood labile is spent plus a fraction of root
!             rlab_to_ra(n) = (labile_storage_overload*deltat_1(n))-wlab_to_ra(n)
!             ! leaving the foliage in storage
!             flab_to_ra(n) = 0.0
!         else if (((wlab_to_ra(n)+rlab_to_ra(n)+flab_to_ra(n))*deltat(n)) >= labile_storage_overload) then
!             ! therefore all wood and root labile is spent plus a fraction of foliage
!             flab_to_ra(n) = (labile_storage_overload*deltat_1(n))-wlab_to_ra(n)-rlab_to_ra(n)
!         else
!             ! oh dear we have more carbon than we have allocated so far
!             ! we spend all wood, root, foliage N limited allocations of labile
!             ! plus a proportional allocation from each pool
!             labile_storage_overload = (labile_storage_overload*deltat_1(n)) &
!                                     - (flab_to_ra(n)+wlab_to_ra(n)+rlab_to_ra(n))
             ! calculate the fractional size of each pool; truncation done to
             ! ensure > 1 never occurs
             C_total = 1/(avail_labfol+avail_labroot+avail_labwood)
             labile_ratios(1) = avail_labfol * C_total
             labile_ratios(2) = avail_labroot * C_total
             labile_ratios(3) = avail_labwood * C_total
             labile_ratios = floor(labile_ratios*1e2)*1e-2
             ! fractionally apply overload correction
!             flab_to_ra(n)=flab_to_ra(n)+(labile_storage_overload*labile_ratios(1))
!             rlab_to_ra(n)=rlab_to_ra(n)+(labile_storage_overload*labile_ratios(2))
!             wlab_to_ra(n)=wlab_to_ra(n)+(labile_storage_overload*labile_ratios(3))
             flab_to_ra(n)=labile_storage_overload*labile_ratios(1)
             rlab_to_ra(n)=labile_storage_overload*labile_ratios(2)
             wlab_to_ra(n)=labile_storage_overload*labile_ratios(3)
!         endif ! how will I get rid of the excess carbon
         ! all all excess to Ra flux
         FLUXES(n,3) = FLUXES(n,3) + flab_to_ra(n) + wlab_to_ra(n) + rlab_to_ra(n)
      else
         ! if we do not have an excess of labile then reset the *_to_ra
         flab_to_ra(n) = 0.0 ; rlab_to_ra(n) = 0.0 ; wlab_to_ra(n) = 0.0
      endif ! overloaded on labile

      ! correct time step to daily equivalent
      Ndemand_leaf = Ndemand_leaf * deltat_1(n)
      Ndemand_root = Ndemand_root * deltat_1(n)
      Ndemand_wood = Ndemand_wood * deltat_1(n)

      !!!
      ! Leaching from DIN
      !!!

      ! restrict based on available DIN
      DIN_leaching = min(avail_DIN,pars(50)*deltat(n))
      ! but convert back into rate
      DIN_leaching = DIN_leaching * deltat_1(n)

      !!!
      ! Calculate final fluxes
      !!!

      ! calculate the NEE
      NEE(n) = (-FLUXES(n,1)+FLUXES(n,3)+FLUXES(n,13)+FLUXES(n,14))
      ! load GPP
      GPP(n) = FLUXES(n,1)

      !
      ! update pools for next timestep
      !

      ! foliar labile pool
      POOLS(n+1,1)  = POOLS(n,1) + (FLUXES(n,5) -FLUXES(n,8)-flab_to_ra(n))*deltat(n)
      ! root labile pool
      POOLS(n+1,7)  = POOLS(n,7) + (npp_to_roots(n)-FLUXES(n,6)-rlab_to_ra(n))*deltat(n)
      ! wood labile pool
      POOLS(n+1,8)  = POOLS(n,8) + (npp_to_wood(n) -FLUXES(n,7)-wlab_to_ra(n))*deltat(n)
      ! foliar pool
      POOLS(n+1,2)  = POOLS(n,2) + (FLUXES(n,8) -FLUXES(n,10))*deltat(n)
      ! wood pool
      POOLS(n+1,4)  = POOLS(n,4) + (FLUXES(n,7) -FLUXES(n,11))*deltat(n)
      ! root pool
      POOLS(n+1,3)  = POOLS(n,3) + (FLUXES(n,6) -FLUXES(n,12))*deltat(n)

      ! litter pool
      POOLS(n+1,5)  = POOLS(n,5) + (FLUXES(n,10)+FLUXES(n,12)+FLUXES(n,19)-FLUXES(n,13)-FLUXES(n,15))*deltat(n)
      ! Coarse woody debris
      POOLS(n+1,9)  = POOLS(n,9) + (FLUXES(n,11)-FLUXES(n,19))*deltat(n)
      ! som pool
      POOLS(n+1,6)  = POOLS(n,6) + (FLUXES(n,15)-FLUXES(n,14))*deltat(n)
      ! labile N foliage pool
      POOLS(n+1,10) = POOLS(n,10) + (foliarN_retrans+Ndemand_leaf-(FLUXES(n,8)*NC_foliar))*deltat(n)
      ! litter N pool
      POOLS(n+1,11) = POOLS(n,11) + (foliarN_loss + &
                                     (FLUXES(n,12)*NC_root)+(FLUXES(n,19)*NC_wood) - &
                                      FLUXES(n,20))*deltat(n)
      ! DIN pool
      POOLS(n+1,12) = POOLS(n,12) + (N_deposition+(gross_nmin*(1d0-DON_leaching))-Ndemand_leaf-Ndemand_root-Ndemand_wood &
                                    -FLUXES(n,21)-DIN_leaching)*deltat(n)

      ! should not need this just sanity check
!      POOLS(n+1,1)  = max(0d0,POOLS(n+1,1))
!      POOLS(n+1,2)  = max(0d0,POOLS(n+1,2))
!      POOLS(n+1,3)  = max(0d0,POOLS(n+1,3))
!      POOLS(n+1,4)  = max(0d0,POOLS(n+1,4))
!      POOLS(n+1,5)  = max(0d0,POOLS(n+1,5))
!      POOLS(n+1,6)  = max(0d0,POOLS(n+1,6))
!      POOLS(n+1,7)  = max(0d0,POOLS(n+1,7))
!      POOLS(n+1,8)  = max(0d0,POOLS(n+1,8))
!      POOLS(n+1,9)  = max(0d0,POOLS(n+1,9))
!      POOLS(n+1,10) = max(0d0,POOLS(n+1,10))
!      POOLS(n+1,11) = max(0d0,POOLS(n+1,11))
!      POOLS(n+1,12) = max(0d0,POOLS(n+1,12))

      !
      ! deal first with deforestation
      !

      if (n == reforest_day) then
          ! replanting conditions
          POOLS(n+1,1) = pars(30) ! fol labile
          POOLS(n+1,2) = pars(31) ! fol
          POOLS(n+1,3) = pars(32) ! root
          POOLS(n+1,4) = pars(33) ! wood
          POOLS(n+1,7) = pars(41) ! root labile
          POOLS(n+1,8) = pars(42) ! wood labile
          POOLS(n+1,10)= pars(57) ! labN fol
      endif

      if (met(8,n) > 0.0) then

          ! pass harvest management to local integer
          harvest_management = int(met(13,n))

          ! assume that labile is proportionally distributed through the plant
          ! and therefore so is the residual fraction
          C_total = POOLS(n+1,2) + POOLS(n+1,3) + POOLS(n+1,7) + POOLS(n+1,4)
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
          labile_root_loss = POOLS(n+1,7)*met(8,n)
          labile_wood_loss = POOLS(n+1,8)*met(8,n)
          foliar_loss = POOLS(n+1,2)*met(8,n)
          roots_loss  = POOLS(n+1,3)*met(8,n)
          wood_loss   = POOLS(n+1,4)*met(8,n)
          ! transfer fraction of harvest waste to litter or som pools
          ! easy pools first
          labile_residue = POOLS(n+1,1)*met(8,n)*labile_frac_res
          labile_root_residue = POOLS(n+1,7)*met(8,n)*labile_frac_res
          labile_wood_residue = POOLS(n+1,8)*met(8,n)*labile_frac_res
          foliar_residue = POOLS(n+1,2)*met(8,n)*foliage_frac_res(harvest_management)
          roots_residue  = POOLS(n+1,3)*met(8,n)*roots_frac_res(harvest_management)
          ! explicit calculation of the residues from each fraction
          coarse_root_residue  = Crootcr*met(8,n)*rootcr_frac_res(harvest_management)
          branch_residue = Cbranch*met(8,n)*branch_frac_res(harvest_management)
          stem_residue = Cstem*met(8,n)*stem_frac_res(harvest_management)
          ! now finally calculate the final wood residue
          wood_residue = stem_residue + branch_residue + coarse_root_residue
          ! mechanical loss of Csom due to coarse root extraction
          soil_loss_with_roots = Crootcr*met(8,n)*(1.0-rootcr_frac_res(harvest_management)) &
                              * soil_loss_frac(harvest_management)

          ! update living pools directly
          POOLS(n+1,1) = max(0.0,POOLS(n+1,1)-labile_loss)
          POOLS(n+1,7) = max(0.0,POOLS(n+1,8)-labile_root_loss)
          POOLS(n+1,8) = max(0.0,POOLS(n+1,9)-labile_wood_loss)
          POOLS(n+1,2) = max(0.0,POOLS(n+1,2)-foliar_loss)
          POOLS(n+1,3) = max(0.0,POOLS(n+1,3)-roots_loss)
          POOLS(n+1,4) = max(0.0,POOLS(n+1,4)-wood_loss)
          ! then work out the adjustment due to burning if there is any
          if (post_harvest_burn(harvest_management) > 0.0) then
              !/*first fluxes*/
              !/*LABILE*/
              CFF(1) = POOLS(n+1,1)*post_harvest_burn(harvest_management)*combust_eff(1)
              NCFF(1) = POOLS(n+1,1)*post_harvest_burn(harvest_management)*(1.0-combust_eff(1))*(1.0-rfac)
              CFF_res(1) = labile_residue*post_harvest_burn(harvest_management)*combust_eff(1)
              NCFF_res(1) = labile_residue*post_harvest_burn(harvest_management)*(1.0-combust_eff(1))*(1.0-rfac)
              !/*ROOTLABILE*/
              CFF(6) = POOLS(n+1,7)*post_harvest_burn(harvest_management)*combust_eff(1)
              NCFF(6) = POOLS(n+1,7)*post_harvest_burn(harvest_management)*(1.0-combust_eff(1))*(1.0-rfac)
              CFF_res(6) = labile_root_residue*post_harvest_burn(harvest_management)*combust_eff(1)
              NCFF_res(6) = labile_root_residue*post_harvest_burn(harvest_management)*(1.0-combust_eff(1))*(1.0-rfac)
              !/*WOODLABILE*/
              CFF(7) = POOLS(n+1,8)*post_harvest_burn(harvest_management)*combust_eff(1)
              NCFF(7) = POOLS(n+1,8)*post_harvest_burn(harvest_management)*(1.0-combust_eff(1))*(1.0-rfac)
              CFF_res(7) = labile_wood_residue*post_harvest_burn(harvest_management)*combust_eff(1)
              NCFF_res(7) = labile_wood_residue*post_harvest_burn(harvest_management)*(1.0-combust_eff(1))*(1.0-rfac)
              !/*foliar*/
              CFF(2) = POOLS(n+1,2)*post_harvest_burn(harvest_management)*combust_eff(2)
              NCFF(2) = POOLS(n+1,2)*post_harvest_burn(harvest_management)*(1.0-combust_eff(2))*(1.0-rfac)
              CFF_res(2) = foliar_residue*post_harvest_burn(harvest_management)*combust_eff(2)
              NCFF_res(2) = foliar_residue*post_harvest_burn(harvest_management)*(1.0-combust_eff(2))*(1.0-rfac)
              !/*root*/
              CFF(3) = 0. !(POOLS(n+1,3)+POOLS(n+1,7))*post_harvest_burn(harvest_management)*combust_eff(3)
              NCFF(3) = 0. !(POOLS(n+1,3)+POOLS(n+1,7))*post_harvest_burn(harvest_management)*(1.0-combust_eff(3))*(1.0-rfac)
              CFF_res(3) = 0. !roots_residue*post_harvest_burn(harvest_management)*combust_eff(3)
              NCFF_res(3) = 0. !roots_residue*post_harvest_burn(harvest_management)*(1.0-combust_eff(3))*(1.0-rfac)
              !/*wood*/
              CFF(4) = POOLS(n+1,4)*post_harvest_burn(harvest_management)*combust_eff(4)
              NCFF(4) = POOLS(n+1,4)*post_harvest_burn(harvest_management)*(1.0-combust_eff(4))*(1.0-rfac)
              CFF_res(4) = wood_residue*post_harvest_burn(harvest_management)*combust_eff(4)
              NCFF_res(4) = wood_residue*post_harvest_burn(harvest_management)*(1.0-combust_eff(4))*(1.0-rfac)
              !/*litter*/
              CFF(5) = POOLS(n+1,5)*post_harvest_burn(harvest_management)*combust_eff(5)
              NCFF(5) = POOLS(n+1,5)*post_harvest_burn(harvest_management)*(1.0-combust_eff(5))*(1.0-rfac)
              !/*CWD*/
              CFF(8) = POOLS(n+1,9)*post_harvest_burn(harvest_management)*combust_eff(4)
              NCFF(8) = POOLS(n+1,9)*post_harvest_burn(harvest_management)*(1.0-combust_eff(4))*(1.0-rfac)
              !/*fires as daily averages to comply with units*/
              FLUXES(n,17)=(CFF(1)+CFF(2)+CFF(3)+CFF(4)+CFF(5)+CFF(6)+CFF(7)+CFF(8) &
                           +CFF_res(1)+CFF_res(2)+CFF_res(3)+CFF_res(4)+CFF_res(6)+CFF_res(7))*deltat_1(n)
              ! update the residue terms
              labile_residue = labile_residue - CFF_res(1) - NCFF_res(1)
              labile_root_residue = labile_root_residue - CFF_res(6) - NCFF_res(6)
              labile_wood_residue = labile_wood_residue - CFF_res(7) - NCFF_res(7)
              foliar_residue = foliar_residue - CFF_res(2) - NCFF_res(2)
              roots_residue  = roots_residue  - CFF_res(3) - NCFF_res(3)
              wood_residue   = wood_residue   - CFF_res(4) - NCFF_res(4)
              ! now update NEE
              NEE(n)=NEE(n)+FLUXES(n,17)
          else
              FLUXES(n,17) = 0.0
              CFF = 0.0 ; NCFF = 0.0
              CFF_res = 0.0 ; NCFF_res = 0.0
          end if
          ! update all pools this time
          POOLS(n+1,1) = max(0.0, POOLS(n+1,1) - CFF(1) - NCFF(1) )
          POOLS(n+1,7) = max(0.0, POOLS(n+1,7) - CFF(6) - NCFF(6) )
          POOLS(n+1,8) = max(0.0, POOLS(n+1,8) - CFF(7) - NCFF(7) )
          POOLS(n+1,2) = max(0.0, POOLS(n+1,2) - CFF(2) - NCFF(2) )
          POOLS(n+1,3) = max(0.0, POOLS(n+1,3) - CFF(3) + NCFF(3) )
          POOLS(n+1,4) = max(0.0, POOLS(n+1,4) - CFF(4) - NCFF(4) )
          POOLS(n+1,5) = max(0.0, POOLS(n+1,5) + (labile_residue+foliar_residue+roots_residue) + (NCFF(1)+NCFF(2)+NCFF(3)) )
          POOLS(n+1,6) = max(0.0, POOLS(n+1,6) + (soil_loss_with_roots) + (NCFF(4)+NCFF(5)+NCFF(8)))
          POOLS(n+1,9) = max(0.0, POOLS(n+1,9) + wood_residue - CFF(8) - NCFF(8) )
          ! some variable needed for the EDCs
          ! reallocation fluxes for the residues
          disturbance_residue_to_litter(n) = (labile_residue+labile_root_residue+labile_wood_residue+foliar_residue+roots_residue) &
                                           + (NCFF(1)+NCFF(2)+NCFF(3))
          disturbance_residue_to_som(n) = soil_loss_with_roots + (NCFF(4)+NCFF(5)+NCFF(8))
          disturbance_residue_to_cwd(n) = wood_residue - CFF(8) - NCFF(8)

          ! this is intended for use with the R interface for subsequent post
          ! processing
          if (allocated(extracted_C)) then
             ! harvested carbon from all pools
            extracted_C(n) = (wood_loss-(wood_residue+CFF_res(4)+NCFF_res(4))) &
                           + (labile_loss-(labile_residue+CFF_res(1)+NCFF_res(1))) &
                           + (labile_root_loss-(labile_root_residue+CFF_res(6)+NCFF_res(6))) &
                           + (labile_wood_loss-(labile_wood_residue+CFF_res(7)+NCFF_res(7))) &
                           + (foliar_loss-(foliar_residue+CFF_res(2)+NCFF_res(2))) &
                           + (roots_loss-(roots_residue+CFF_res(3)+NCFF_res(3)))
          endif ! allocated extracted_C
          ! total carbon loss from the system
          C_total = (labile_residue+foliar_residue+roots_residue+wood_residue+sum(NCFF)+labile_root_residue+labile_wood_residue) &
                  - (labile_loss+foliar_loss+roots_loss+wood_loss+soil_loss_with_roots+sum(CFF)+labile_root_loss+labile_wood_loss)

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

      if (met(9,n) > 0.0) then

         !/*first fluxes*/
         !/*LABILE*/
         CFF(1) = POOLS(n+1,1)*met(9,n)*combust_eff(1)
         NCFF(1) = POOLS(n+1,1)*met(9,n)*(1.0-combust_eff(1))*(1.0-rfac)
         !/*ROOTLABILE*/
         CFF(6) = POOLS(n+1,7)*met(9,n)*combust_eff(1)
         NCFF(6) = POOLS(n+1,7)*met(9,n)*(1.0-combust_eff(1))*(1.0-rfac)
         !/*WOODLABILE*/
         CFF(7) = POOLS(n+1,8)*met(9,n)*combust_eff(1)
         NCFF(7) = POOLS(n+1,8)*met(9,n)*(1.0-combust_eff(1))*(1.0-rfac)
         !/*foliar*/
         CFF(2) = POOLS(n+1,2)*met(9,n)*combust_eff(2)
         NCFF(2) = POOLS(n+1,2)*met(9,n)*(1.0-combust_eff(2))*(1.0-rfac)
        !/*root*/
         CFF(3) = 0.0  ! (POOLS(n+1,3)+POOLS(n+1,7))*met(9,n)*combust_eff(3)
         NCFF(3) = 0.0 ! (POOLS(n+1,3)+POOLS(n+1,7))*met(9,n)*(1.0-combust_eff(3))*(1.0-rfac)
         !/*wood*/
         CFF(4) = POOLS(n+1,4)*met(9,n)*combust_eff(4)
         NCFF(4) = POOLS(n+1,4)*met(9,n)*(1.0-combust_eff(4))*(1.0-rfac)
         !/*litter*/
         CFF(5) = POOLS(n+1,5)*met(9,n)*combust_eff(5)
         NCFF(5) = POOLS(n+1,5)*met(9,n)*(1.0-combust_eff(5))*(1.0-rfac)
         !/*CWD*/
         CFF(8) = POOLS(n+1,9)*met(9,n)*combust_eff(4)
         NCFF(8) = POOLS(n+1,9)*met(9,n)*(1.0-combust_eff(4))*(1.0-rfac)
         !/*fires as daily averages to comply with units*/
         FLUXES(n,17)=(CFF(1)+CFF(2)+CFF(3)+CFF(4)+CFF(5)+CFF(6)+CFF(7)+CFF(8))*deltat_1(n)

         !/*all fluxes are at a daily timestep*/
         NEE(n)=NEE(n)+FLUXES(n,17)

         !// update pools
         !/*Adding all fire pool transfers here*/
         POOLS(n+1,1)=POOLS(n+1,1)-CFF(1)-NCFF(1)
         POOLS(n+1,7)=POOLS(n+1,7)-CFF(6)-NCFF(6)
         POOLS(n+1,8)=POOLS(n+1,8)-CFF(7)-NCFF(7)
         POOLS(n+1,2)=POOLS(n+1,2)-CFF(2)-NCFF(2)
         POOLS(n+1,3)=POOLS(n+1,3)-CFF(3)+NCFF(3)
         POOLS(n+1,4)=POOLS(n+1,4)-CFF(4)-NCFF(4)
         POOLS(n+1,5)=POOLS(n+1,5)-CFF(5)-NCFF(5)+NCFF(1)+NCFF(2)+NCFF(3)+NCFF(6)+NCFF(7)
         POOLS(n+1,6)=POOLS(n+1,6)+NCFF(4)+NCFF(5)+NCFF(8)
         POOLS(n+1,9)=POOLS(n+1,9)
         ! some variable needed for the EDCs
         ! reallocation fluxes for the residues
         disturbance_residue_to_litter(n) = NCFF(1)+NCFF(2)+NCFF(3)
         disturbance_residue_to_som(n) = NCFF(5)
         disturbance_residue_to_cwd(n) = 0.0
      endif ! end burnst area issues

    do nxp = 1, nopools
       if (POOLS(n+1,nxp) /= POOLS(n+1,nxp) .or. POOLS(n+1,nxp) < -1.0) then
          print*,"step",n, nxp
          print*,"met",met(:,n)
          print*,"POOLS",POOLS(n,:)
          print*,"FLUXES",FLUXES(n,:)
          print*,"lab->Ra",flab_to_ra(n),rlab_to_ra(n),wlab_to_ra(n)
          print*,"POOLS+1",POOLS(n+1,:)
          stop
       endif
    enddo

    end do ! nodays loop

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
    ci=0.5*(co2+qq-pp+sqrt(((co2+qq-pp)*(co2+qq-pp))-4.0*(co2*qq-pp*co2_comp_point)))
    ! limit maximum quantium efficiency by leaf area, hyperbola
    e0=lai_coef*(lai*lai)/((lai*lai)+lai_const)
    ! calculate day length (hours)
    dec = - asin( sin( 23.45 * pi / 180.0 ) * cos( 2.0 * pi * ( doy + 10.0 ) /365.0 ) )
    sinld = sin( lat*(pi/180.0) ) * sin( dec )
    cosld = cos( lat*(pi/180.0) ) * cos( dec )
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
