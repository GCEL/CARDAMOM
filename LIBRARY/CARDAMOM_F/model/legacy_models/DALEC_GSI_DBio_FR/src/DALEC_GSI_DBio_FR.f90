
module CARBON_MODEL_MOD

implicit none

  !!!!!!!!!!!
  ! Authorship contributions
  !
  ! This code contains a variant of the Data Assimilation Linked ECosystem (DALEC) model.
  ! This version of DALEC is derived from the following primary references:
  ! Xenakis & Williams (2014) https://doi.org/10.5194/gmd-7-1519-2014
  ! Bloom & Williams (2015), https://doi.org/10.5194/bg-12-1299-2015.
  ! Smallman et al., (2017), https://doi.org/10.1002/2016JG003520.
  ! Smallman & Williams (2019) https://doi.org/10.5194/gmd-12-2227-2019.
  ! Thomas et al., (2019), https://doi.org/10.1029/2019MS001679
  ! This code is based on that created by A. A. Bloom (UoE, now at JPL, USA).
  ! Subsequent modifications by:
  ! T. L. Smallman (University of Edinburgh, t.l.smallman@ed.ac.uk)
  ! J. F. Exbrayat (University of Edinburgh)
  ! G. Xenakis (University of Edinburgh, now at Forest Research, UK)
  ! See function / subroutine specific comments for exceptions and contributors
  !!!!!!!!!!!

! make all private
private

! explicit publics
public :: CARBON_MODEL     &
         ,soil_frac_clay   &
         ,soil_frac_sand   &
         ,nos_soil_layers  &
         ,disturbance_residue_to_flitter &
         ,disturbance_loss_from_flitter  &
         ,disturbance_residue_to_rlitter &
         ,disturbance_loss_from_rlitter  &
         ,disturbance_residue_to_wlitter &
         ,disturbance_loss_from_wlitter  &
         ,disturbance_residue_to_som     &
         ,disturbance_loss_from_som      &
         ,itemp,ivpd,iphoto &
         ,extracted_C       &
         ,flignin           &
         ,microbial_activity_out &
         ,dim_1,dim_2       &
         ,nos_trees         &
         ,nos_inputs        &
         ,leftDaughter      &
         ,rightDaughter     &
         ,nodestatus        &
         ,xbestsplit        &
         ,nodepred          &
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
                   ,gsi_history(22), just_grown, meant

! arrays for the emulator, just so we load them once and that is it cos they be
! massive
integer ::    dim_1, & ! dimension 1 of response surface
              dim_2, & ! dimension 2 of response surface
          nos_trees, & ! number of trees in randomForest
         nos_inputs    ! number of driver inputs

! DBio specific information
double precision, allocatable, dimension(:) :: microbial_activity_out, &
                                               disturbance_residue_to_flitter, &
                                               disturbance_loss_from_flitter,  &
                                               disturbance_residue_to_rlitter, &
                                               disturbance_loss_from_rlitter,  &
                                               disturbance_residue_to_wlitter, &
                                               disturbance_loss_from_wlitter,  &
                                               disturbance_residue_to_som,    &
                                               disturbance_loss_from_som

integer :: gsi_lag_remembered
double precision, allocatable, dimension(:) :: tmp_x, tmp_m
double precision, allocatable, dimension(:,:) ::     leftDaughter, & ! left daughter for forest
                                                    rightDaughter, & ! right daughter for forets
                                                       nodestatus, & ! nodestatus for forests
                                                       xbestsplit, & ! for forest
                                                         nodepred, & ! prediction value for each tree
                                                          bestvar    ! for randomForests
! microbial decomposition model
double precision, dimension(3) :: flignin
double precision ::  microbial_activity           & ! microbial activity for decompostion
                    ,dmact,microbial_death

! parameters needed for root modelling
integer, parameter :: nos_root_layers = 2, nos_soil_layers = nos_root_layers + 1
double precision, dimension(nos_soil_layers) :: soil_frac_clay,soil_frac_sand
double precision, parameter::   gplant = 5.0,          & ! plant hydraulic conductivity (mmol m-1 s-1 MPa-1)
                                  head = 0.009807,     & ! head of pressure (MPa/m)
                                head_1 = 101.968,      & ! inverse head of pressure (m/MPa)
                           root_radius = 0.0001,       & ! root radius (m) Bonen et al 2014 = 0.00029
                   root_cross_sec_area = 3.141593e-08, & ! root cross sectional area (m2)
                                                         ! = pi * root_radius *
                                                         ! root_radius
                          root_density = 0.5e6,        & ! root density (g biomass m-3 root)
                                                         ! 0.5e6 Williams et al
                                                         ! 1996
                                                         ! 0.31e6 Bonan et al
                                                         ! 2014
                         canopy_height = 9.0,          & ! canopy height assumed to be 9 m
                        top_soil_depth = 0.3,          & ! depth to which we conider the top soil to extend (m)
                              min_root = 5.0,          & ! minimum rot biomass (gBiomass.m-2)
                             max_depth = 1.5,          & ! maximum possible root depth (m)
                                root_k = 100,          & ! biomass to reach half max_depth
                           root_resist = 25.0            ! Root resistivity (MPa s g mmol−1 H2O)
! scalers needed for root modelling
double precision :: root_reach, root_biomass, demand
! arrays needed for root modelling
double precision :: layer_thickness(nos_soil_layers), water_flux(nos_root_layers)

! explicit save statment for the variables in here
save

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
    ! Modified to include a version of microbial decomposition model described in Xenekis & Williams (2013).

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
    double precision :: gpppars(11)        & ! ACM inputs (LAI+met)
             ,constants(10)                  ! parameters for ACM

    integer :: p,f,nxp,n,m,test

    ! local fire related variables
    double precision :: CFF(6) = 0.0, CFF_res(4) = 0.0    & ! combusted and non-combustion fluxes
                       ,NCFF(6) = 0.0, NCFF_res(4) = 0.0  & ! with residue and non-residue seperates
                       ,combust_eff(5)                    & ! combustion efficiency
                       ,rfac                                ! resilience factor

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

    integer :: reforest_day, harvest_management, restocking_lag, gsi_lag

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
    ! 1 = labile        (p18)
    ! 2 = foliar        (p19)
    ! 3 = root          (p20)
    ! 4 = wood          (p21)
    ! 5 = som_fast      (p22)
    ! 6 = som_slow      (p23)
    ! 7 = foliar_litter (p36)
    ! 8 = root_litter   (p37)
    ! 9 = wood_litter   (p38)
    ! 10 = microbial    (p39)

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
    ! 10 = leaflitter production
    ! 11 = woodlitter production
    ! 12 = rootlitter production
    ! 13 = respiration het foliar litter
    ! 14 = respiration het root litter
    ! 15 = respiration het wood litter
    ! 16 = labrelease factor
    ! 17 = carbon flux due to fire
    ! 18 = growing season index
    ! 19 = respiraiton het som fast
    ! 20 = respiration het som slow
    ! 21 = respiration het microbial
    ! 22 = decomposition from foliar litter to som (fast+slow)
    ! 23 = decomposition from root litter to som (fast+slow)
    ! 24 = decomposition from wood litter to som (fast+slow)
    ! 25 = microbial death allocation to slow
    ! 26 = microbe mediated transfer from slow to fast
    ! 27 = accumulation of som_fast into Cmicrobe

    ! PARAMETERS
    ! 17+4(GSI)+8(DBio)

    ! p(1) 1st stage decompostion efficiency (DBio)
    ! p(2) Fraction of GPP respired
    ! p(3) Fraction of NPP allocated to foliage
    ! p(4) Fraction of NPP allocated to roots
    ! p(5) Max leaf turnover
    ! p(6) Turnover rate of wood
    ! p(7) Turnover rate of roots
    ! p(8) Foliage litter turnover rate  (DBio)
    ! p(9) Microbial decomposition efficiency (DBio)
    ! p(10) Parameter in exponential term of temperature
    ! p(11) Canopy efficiency parameter - C_eff (part of ACM)
    ! p(12) = Max labile turnover(GSI)
    ! p(13) = Fraction allocated to Clab
    ! p(14) = Min temp threshold (GSI)
    ! p(15) = Max temp threshold (GSI)
    ! p(16) = Min photoperiod threshold (GIS)
    ! p(17) = LMA
    ! p(24) = Max photoperiod threshold (GSI)
    ! p(25) = Min VPD threshold (GSI)
    ! p(26) = Max VPD threshold (GSI)
    ! p(27) = Root litter turnover rate (DBio)
    ! p(28) = Wood litter turnover rate (DBio)
    ! p(29) = Efficiency of substrate uptake by microbes (DBio)
    ! p(30) = 2nd order rate constant for microbial uptake (DBio)
    ! p(31) = Maximum microbial death rate (DBio)
    ! p(32) = Inhibition constant for microbial death (DBio)
    ! p(33) = Maintenance coefficient (DBio)
    ! p(34) = Inhibition constant for C dependant microbial activity (DBio)
    ! p(35) = Decomposition of som_slow (DBio)
    ! p(36) = Initial microbial activity (DBio)
    ! p(41) = Foliar lignin fraction
    ! p(42) = Fine root lignin fraction
    ! p(43) = Wood lignin fraction
    ! p(44) = minimum GPP benefit of increased LAI for labile allocation to be allowed

    ! DBio model local parameters
    ! fligin(1:3) = fraction of ligin content for foliage, root and wood

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
    gpppars(4) = 10.0**pars(11) !TLS 1 ! foliar N
    gpppars(7) = lat
    gpppars(9) = -2.060814 !-2.0  ! leafWP-soilWP (i.e. -2-0) ! p11 from ACM recal
    gpppars(10) = 1.0 ! totaly hydraulic resistance ! p12 from ACM recal (updated)
    gpppars(11) = pi

    ! assign acm parameters
    constants(1)=24.24129      ! pars(11) ! p1  from ACM recal
    constants(2)=7.798524e-03  ! 0.0156935! p2  from ACM recal
    constants(3)=154.7495      ! 4.22273  ! p3  from ACM recal
    constants(4)=465.7482      ! 208.868  ! p4  from ACM recal
    constants(5)=7.817923e-02  ! 0.0453194! p5  from ACM recal
    constants(6)=5.674312e-01  ! 0.37836  ! p6  from ACM recal
    constants(7)=1.076729e+01  ! 7.19298  ! p7  from ACM recal
    constants(8)=5.577107e-03  ! 0.011136 ! p8  from ACM recal
    constants(9)=3.154374e+00  ! 2.1001   ! p9  from ACM recal
    constants(10)=3.959395e-01 ! 0.789798 ! p10 from ACM recal

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
    Cbranch_part = pars(45)
    Crootcr_part = pars(46)

    ! declare fire constants
    ! labile, foliar, root, wood, litter
    combust_eff(1) = 0.1 ; combust_eff(2) = 0.9
    combust_eff(3) = 0.1 ; combust_eff(4) = 0.5
    combust_eff(5) = 0.3 ; rfac = 0.5

    ! declare microbial model constants
    flignin(1)=pars(41) ; flignin(2)=pars(42) ; flignin(3)=pars(43)

    if (start == 1) then

        ! assigning initial conditions
        POOLS(1,1)=pars(18)
        POOLS(1,2)=pars(19)
        POOLS(1,3)=pars(20)
        POOLS(1,4)=pars(21)
        POOLS(1,5)=pars(22)
        POOLS(1,6)=pars(23)
        POOLS(1,7)=pars(37)
        POOLS(1,8)=pars(38)
        POOLS(1,9)=pars(39)
        POOLS(1,10)=pars(40)

        ! load initial value for microbial activity
        microbial_activity=pars(36)

        ! allocate memory to variables needed for EDCs if needed
        if (.not.allocated(disturbance_residue_to_som)) then
            allocate(disturbance_residue_to_flitter(nodays), &
                     disturbance_loss_from_flitter(nodays),  &
                     disturbance_residue_to_rlitter(nodays), &
                     disturbance_loss_from_rlitter(nodays),  &
                     disturbance_residue_to_wlitter(nodays), &
                     disturbance_loss_from_wlitter(nodays),  &
                     disturbance_residue_to_som(nodays),     &
                     disturbance_loss_from_som(nodays))
        endif
        ! zero initial values
        disturbance_residue_to_flitter = 0.0 ; disturbance_loss_from_flitter = 0.0
        disturbance_residue_to_rlitter = 0.0 ; disturbance_loss_from_rlitter = 0.0
        disturbance_residue_to_wlitter = 0.0 ; disturbance_loss_from_wlitter = 0.0
        disturbance_residue_to_som = 0.0 ; disturbance_loss_from_som = 0.0

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
!print*,"1gsi_lag",gsi_lag
        end if ! .not.allocated(tmp_x)

        ! assign our starting value
        gsi_history = pars(53)-1d0
        just_grown = pars(52)

    endif ! start == 1

    ! assign climate sensitivities
    fol_turn_crit=pars(34)-1d0
    lab_turn_crit=pars(3)-1d0
    gsi_lag = gsi_lag_remembered ! added to prevent loss from memory

    !
    ! Begin looping through each time step
    !

    do n = start, finish

      ! calculate LAI value
      lai(n)=POOLS(n,2)/pars(17)
      ! calculate mean temperature
      meant = 0.5 * (met(3,n)+met(2,n))

      ! load next met / lai values for ACM
      gpppars(1)=lai(n)
      gpppars(2)=met(3,n) ! max temp
      gpppars(3)=met(2,n) ! min temp
      gpppars(5)=met(5,n) ! co2
      gpppars(6)=ceiling(met(6,n)-(deltat(n)*0.5)) ! doy
      gpppars(8)=met(4,n) ! radiation

      ! calculate the minimum soil & root hydraulic resistance based on total
      ! fine root mass ! *2*2 => *RS*C->Bio
      root_biomass = max(min_root,POOLS(n,3)*2)
      call calculate_Rtot(abs(gpppars(9)),gpppars(10),lai(n) &
                         ,deltat(n),meant)

      ! GPP (gC.m-2.day-1)
      if (lai(n) > 0.0) then
         FLUXES(n,1) = acm(gpppars,constants)
      else
         FLUXES(n,1) = 0.0
      endif

      ! temprate (i.e. temperature modified rate of metabolic activity))
      FLUXES(n,2) = exp( ( log(pars(10)) * 0.1 ) * meant  )
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

      ! DecoBio v1.0 (Xenakis & Williams 2014, Comparing microbial and chemical
      ! kinetics for modelling soil organic carbon decomposition using the
      ! DecoChem 1.0 and DecoBio v1.0 models, Geoscientific Model Development,
      ! 7, 1519-1533)

      ! A modified version of the microbe mediated decomposition model was added by TLS on 19/12/2014.
      ! The model was included to provide a more realistic representation in
      ! soil carbon cycles. This should generate a more appropriate response to
      ! management practices, particulary those which generate harvest residue
      ! which could provide a priming effect. The model here assumes that ligin
      ! values are fixed to result in wood litter allocation to the Csom_slow
      ! pool while foliar and root litters are allocated to the Csom_fast.

      ! Calculate the dynamic parameter of microbial activity in relation
      ! to soluble, cellulose and ligning
      dmact = FLUXES(n,2) * pars(30) * POOLS(n,5) * ( ( POOLS(n,5)  / ( POOLS(n,5) + pars(51) ) ) - microbial_activity )
      ! sanity check
      if ( dmact < 0.0 .and. abs(dmact) >= microbial_activity ) then
         dmact = 0d0 ; microbial_activity = 0d0
      else
         microbial_activity = min(1d0,microbial_activity + dmact)
      endif

      ! Calculate the microbial death rate
      microbial_death = pars(31) / (1d0 + (pars(32) * POOLS(n,5)))

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
      VPDfac = 1d0 - ( (met(12,n)-pars(25)) / (pars(26)-pars(25)) )
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
      if (gradient < fol_turn_crit .or. FLUXES(n,18) == 0.0) then
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
         if ( ((tmp - FLUXES(n,1))/FLUXES(n,1)) < pars(44) ) then
             FLUXES(n,16) = 0d0
         endif
      else
         ! probably we want nothing to happen, however if we are at the seasonal
         ! maximum we will consider further growth stil
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
            if ( ((tmp - FLUXES(n,1))/FLUXES(n,1)) < pars(44) ) then
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

      ! respiration heterotrophic foliar litter, root litter, wood litter
      FLUXES(n,13) = POOLS(n,7)*((1d0-pars(1))*FLUXES(n,2)*pars(8))
      FLUXES(n,14) = POOLS(n,8)*((1d0-pars(1))*FLUXES(n,2)*pars(27))
      FLUXES(n,15) = POOLS(n,9)*((1d0-pars(1))*FLUXES(n,2)*pars(28))
      ! fluxes from litter pools to slow and fast
      FLUXES(n,22) = POOLS(n,7)*(pars(1)*FLUXES(n,2)*pars(8))
      FLUXES(n,23) = POOLS(n,8)*(pars(1)*FLUXES(n,2)*pars(27))
      FLUXES(n,24) = POOLS(n,9)*(pars(1)*FLUXES(n,2)*pars(28))

      ! respiration heterotrophic som fast and slow, and microbial
      FLUXES(n,19) = (1.0-pars(29))*pars(30)*POOLS(n,5)*POOLS(n,10)*microbial_activity
      FLUXES(n,20) = (1.0-pars(9))*min(POOLS(n,6), microbial_activity*pars(35)*POOLS(n,10)*deltat(n) )
      if (FLUXES(n,20) > 0.0) FLUXES(n,20) = FLUXES(n,20) / deltat(n)
      FLUXES(n,21) = pars(33)*microbial_activity*POOLS(n,10)
      ! microbial death allocation to slow
      FLUXES(n,25) = microbial_death*microbial_activity*POOLS(n,10)
      ! microbe mediated transfer from slow to fast
      FLUXES(n,26) = pars(9)*min(POOLS(n,6), microbial_activity*pars(35)*POOLS(n,10)*deltat(n) )
      if (FLUXES(n,26) > 0) FLUXES(n,26) = FLUXES(n,26)/deltat(n)
      ! accumulation of som_fast into Cmicrobe
      FLUXES(n,27) = pars(29)*pars(30)*POOLS(n,5)*POOLS(n,10)*microbial_activity

      ! specifically for post processing
      if (allocated(microbial_activity_out)) then
         microbial_activity_out(n) = microbial_activity
      endif

      ! calculate growth respiration and adjust allocation to pools assuming
      ! 0.21875 of total C allocation towards each pool (i.e. 0.28 .eq. xNPP)
      ! foliage
      FLUXES(n,3) = FLUXES(n,3) + (FLUXES(n,8)*0.21875) ; FLUXES(n,8) = FLUXES(n,8) * 0.78125
      ! roots
      FLUXES(n,3) = FLUXES(n,3) + (FLUXES(n,6)*0.21875) ; FLUXES(n,6) = FLUXES(n,6) * 0.78125
      ! wood
      FLUXES(n,3) = FLUXES(n,3) + (FLUXES(n,7)*0.21875) ; FLUXES(n,7) = FLUXES(n,7) * 0.78125

      ! calculate the NEE
      NEE(n) = -FLUXES(n,1)+FLUXES(n,3) &
               +FLUXES(n,13)+FLUXES(n,14)+FLUXES(n,15) &
               +FLUXES(n,19)+FLUXES(n,20)+FLUXES(n,21)

      ! load GPP
      GPP(n) = FLUXES(n,1)

      !
      ! update pools for next timestep
      !

      ! labile pool
      POOLS(n+1,1)  = POOLS(n,1)  + (FLUXES(n,5)-FLUXES(n,8))*deltat(n)
      ! foliar pool
      POOLS(n+1,2)  = POOLS(n,2)  + (FLUXES(n,4)-FLUXES(n,10) + FLUXES(n,8))*deltat(n)
      ! root pool
      POOLS(n+1,3)  = POOLS(n,3)  + (FLUXES(n,6)-FLUXES(n,12))*deltat(n)
      ! wood pool
      POOLS(n+1,4)  = POOLS(n,4)  + (FLUXES(n,7)-FLUXES(n,11))*deltat(n)
      ! som_fast
      POOLS(n+1,5)  = POOLS(n,5)  + ( ((1d0-flignin(1))*FLUXES(n,22))+ &
                                      ((1d0-flignin(2))*FLUXES(n,23))+ &
                                      ((1d0-flignin(3))*FLUXES(n,24))+ &
                                      FLUXES(n,26)-FLUXES(n,19)-FLUXES(n,27) )*deltat(n)
      ! som_slow
      POOLS(n+1,6)  = POOLS(n,6)  + ( ((flignin(1))*FLUXES(n,22))+ &
                                      ((flignin(2))*FLUXES(n,23))+ &
                                      ((flignin(3))*FLUXES(n,24))+ &
                                      FLUXES(n,25)-FLUXES(n,20) - &
                                      FLUXES(n,26) )*deltat(n)
      ! foliar litter
      POOLS(n+1,7)  = POOLS(n,7)  + ( FLUXES(n,10)-FLUXES(n,13)-FLUXES(n,22) )*deltat(n)
      ! root litter
      POOLS(n+1,8)  = POOLS(n,8)  + ( FLUXES(n,12)-FLUXES(n,14)-FLUXES(n,23) )*deltat(n)
      ! wood litter
      POOLS(n+1,9)  = POOLS(n,9)  + ( FLUXES(n,11)-FLUXES(n,15)-FLUXES(n,24) )*deltat(n)
      ! microbial
      POOLS(n+1,10) = POOLS(n,10) + ( FLUXES(n,27)-FLUXES(n,21)-FLUXES(n,25) )*deltat(n)

      !
      ! deal first with deforestation
      !

      if (n == reforest_day) then
          POOLS(n+1,1) = pars(47)
          POOLS(n+1,2) = pars(48)
          POOLS(n+1,3) = pars(49)
          POOLS(n+1,4) = pars(50)
      endif

      if (met(8,n) > 0) then

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
          soil_loss_with_roots = Crootcr*met(8,n)*(1.0-rootcr_frac_res(harvest_management)) &
                              * soil_loss_frac(harvest_management)

         ! update living pools directly
         POOLS(n+1,1) = max(0.,POOLS(n+1,1)-labile_loss)
         POOLS(n+1,2) = max(0.,POOLS(n+1,2)-foliar_loss)
         POOLS(n+1,3) = max(0.,POOLS(n+1,3)-roots_loss)
         POOLS(n+1,4) = max(0.,POOLS(n+1,4)-wood_loss)
         ! then work out the adjustment due to burning if there is any
         if (post_harvest_burn(harvest_management) > 0) then
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
             CFF(3) = 0.! POOLS(n+1,3)*post_harvest_burn(harvest_management)*combust_eff(3)
             NCFF(3) = 0.! POOLS(n+1,3)*post_harvest_burn(harvest_management)*(1-combust_eff(3))*(1-rfac)
             CFF_res(3) = 0.! roots_residue*post_harvest_burn(harvest_management)*combust_eff(3)
             NCFF_res(3) = 0.! roots_residue*post_harvest_burn(harvest_management)*(1-combust_eff(3))*(1-rfac)
             !/*wood*/
             CFF(4) = POOLS(n+1,4)*post_harvest_burn(harvest_management)*combust_eff(4)
             NCFF(4) = POOLS(n+1,4)*post_harvest_burn(harvest_management)*(1-combust_eff(4))*(1-rfac)
             CFF_res(4) = wood_residue*post_harvest_burn(harvest_management)*combust_eff(4)
             NCFF_res(4) = wood_residue*post_harvest_burn(harvest_management)*(1-combust_eff(4))*(1-rfac)
             !/*litter*/ foliar litter pools
             CFF(5) = POOLS(n+1,7)*post_harvest_burn(harvest_management)*combust_eff(5)
             NCFF(5) = POOLS(n+1,7)*post_harvest_burn(harvest_management)*(1-combust_eff(5))*(1-rfac)
             !/*litter*/ wood litter pool
             CFF(6) = POOLS(n+1,9)*post_harvest_burn(harvest_management)*combust_eff(5)
             NCFF(6) = POOLS(n+1,9)*post_harvest_burn(harvest_management)*(1-combust_eff(5))*(1-rfac)
             !/*fires as daily averages to comply with units*/
             FLUXES(n,17)=(CFF(1)+CFF(2)+CFF(3)+CFF(4)+CFF(5)+CFF(6) &
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
         POOLS(n+1,5) = max(0., POOLS(n+1,5) + NCFF(1) + labile_residue)
         POOLS(n+1,6) = max(0., POOLS(n+1,6) - soil_loss_with_roots + NCFF(4) + NCFF(5) + NCFF(6))
         POOLS(n+1,7) = max(0., POOLS(n+1,7) + foliar_residue + NCFF(2) - CFF(5) - NCFF(5))
         POOLS(n+1,8) = max(0., POOLS(n+1,8) + roots_residue  + NCFF(3) )
         POOLS(n+1,9) = max(0., POOLS(n+1,9) + wood_residue - CFF(6) - NCFF(6))

         ! some variable needed for the EDCs
         ! reallocation fluxes for the residues
         disturbance_residue_to_flitter(n)  = foliar_residue + NCFF(2)
         disturbance_loss_from_flitter(n)   = CFF(5)+NCFF(5)
         disturbance_residue_to_rlitter(n)  = roots_residue + NCFF(3)
         disturbance_loss_from_rlitter(n)   = 0.0
         disturbance_residue_to_wlitter(n)  = wood_residue
         disturbance_loss_from_wlitter(n)   = CFF(6)+NCFF(6)
         disturbance_residue_to_som(n)      = NCFF(4)+NCFF(5)+NCFF(6)
         disturbance_loss_from_som(n)       = soil_loss_with_roots

         ! harvested carbon from all pools
         FLUXES(n,28) = (wood_loss-(wood_residue+CFF_res(4)+NCFF_res(4))) &
                      + (labile_loss-(labile_residue+CFF_res(1)+NCFF_res(1))) &
                      + (foliar_loss-(foliar_residue+CFF_res(2)+NCFF_res(2))) &
                      + (roots_loss-(roots_residue+CFF_res(3)+NCFF_res(3)))

         ! total carbon loss from the system
         C_total = (labile_residue+foliar_residue+roots_residue+wood_residue+sum(NCFF)) &
                 - (labile_loss+foliar_loss+roots_loss+wood_loss+soil_loss_with_roots+sum(CFF))

         ! if total clearance occured then we need to ensure some minimum
         ! values and determining reforestation day
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
         CFF(3) = 0.! POOLS(n+1,3)*met(9,n)*combust_eff(3)
         NCFF(3) = 0.! POOLS(n+1,3)*met(9,n)*(1-combust_eff(3))*(1-rfac)
         !/*wood*/
         CFF(4) = POOLS(n+1,4)*met(9,n)*combust_eff(4)
         NCFF(4) = POOLS(n+1,4)*met(9,n)*(1-combust_eff(4))*(1-rfac)
         !/*litter*/
         CFF(5) = POOLS(n+1,7)*met(9,n)*combust_eff(5)
         NCFF(5) = POOLS(n+1,7)*met(9,n)*(1-combust_eff(5))*(1-rfac)
         !/*wood litter*/
         CFF(6) = POOLS(n+1,9)*met(9,n)*combust_eff(4)
         NCFF(6) = POOLS(n+1,9)*met(9,n)*(1-combust_eff(4))*(1-rfac)

         !/*fires as daily averages to comply with units*/
         FLUXES(n,17)=(CFF(1)+CFF(2)+CFF(3)+CFF(4)+CFF(5)+CFF(6))/deltat(n)

         !/*all fluxes are at a daily timestep*/
         NEE(n)=NEE(n)+FLUXES(n,17)

         !// update pools
         !/*Adding all fire pool transfers here*/
         POOLS(n+1,1) = max(0., POOLS(n+1,1) - CFF(1) - NCFF(1) )
         POOLS(n+1,2) = max(0., POOLS(n+1,2) - CFF(2) - NCFF(2) )
         POOLS(n+1,3) = max(0., POOLS(n+1,3) - CFF(3) - NCFF(3) )
         POOLS(n+1,4) = max(0., POOLS(n+1,4) - CFF(4) - NCFF(4) )
         POOLS(n+1,5) = max(0., POOLS(n+1,5) + NCFF(1) )
         POOLS(n+1,6) = max(0., POOLS(n+1,6) + NCFF(4) + NCFF(5) + NCFF(6))
         POOLS(n+1,7) = max(0., POOLS(n+1,7) + NCFF(2) - CFF(5) - NCFF(5))
         POOLS(n+1,8) = max(0., POOLS(n+1,8) + NCFF(3) )
         POOLS(n+1,9) = max(0., POOLS(n+1,9) - CFF(6) - NCFF(6) )

      endif ! end burnst area issues

!    do nxp = 1, nopools
!       if ((POOLS(n+1,nxp) /= POOLS(n+1,nxp)) .or. (POOLS(n+1,nxp) < 0.0)) then
!          print*,"step",n, nxp
!          print*,"met",met(:,n)
!          print*,"POOLS",POOLS(n,:)
!          print*,"FLUXES",FLUXES(n,:)
!          print*,"microbial_activity",microbial_activity
!          print*,"fastsomturn",pars(29:30)
!          print*,"POOLS+1",POOLS(n+1,:)
!          stop
!       endif
!    enddo

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
    double precision, intent(in) :: drivers(11) & ! acm input requirements
                         ,constants(10) ! ACM parameters

    ! declare local variables
    double precision :: gc, pn, pd, pp, qq, ci, e0, dayl, cps, dec, nit &
             ,trange, sinld, cosld,aob, mult &
             ,mint,maxt,radiation,co2,lai,doy,lat &
             ,deltaWP,Rtot,NUE,temp_exponent,dayl_coef &
             ,dayl_const,hydraulic_exponent,hydraulic_temp_coef &
             ,co2_comp_point,co2_half_sat,lai_coef,lai_const

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
!    trange=0.5*(maxt-mint)
    ! daily canopy conductance, of CO2 or H2O?
    gc=abs(deltaWP)**(hydraulic_exponent)/(hydraulic_temp_coef*Rtot)
    ! maximum rate of temperature and nitrogen (canopy efficiency) limited
    ! photosynthesis (gC.m-2.day-1)
    pn=lai*nit*NUE*exp(temp_exponent*maxt)
    ! pp and qq represent limitation by diffusion and metabolites respecitively
    pp=pn/gc ; qq=co2_comp_point-co2_half_sat
    ! calculate internal CO2 concentration (ppm)
    mult = co2+qq-pp
    ci=0.5*(mult+sqrt(((mult)*(mult))-4.0*(co2*qq-pp*co2_comp_point)))
    ! limit maximum quantium efficiency by leaf area, hyperbola
    mult = lai*lai
    e0=lai_coef*(mult)/((mult)+lai_const)
    ! calculate day length (hours)
    dec = - asin( sin( 23.45 * deg_to_rad ) * cos( 2.0 * pi * ( doy + 10.0 ) / 365.0 ) )
    mult = lat*deg_to_rad
    sinld = sin( mult ) * sin( dec )
    cosld = cos( mult ) * cos( dec )
    aob = max(-1.0,min(1.0,sinld / cosld))
    dayl = 12.0 * ( 1.0 + 2.0 * asin( aob ) / pi )
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
  subroutine calculate_Rtot(deltaWP,Rtot,lai,deltat,meant)

    ! purpose of this function is to calculate the minimum soil-root hydraulic
    ! resistance input into ACM. The minimum is assumed to be the same and the
    ! soil layer with the greated root content. Here we use the same assumption
    ! as used in SPA to calcule the root hydraulic resistance in the top soil
    ! layer only.

    ! This could be extended later for include a 2 soil layer bucket model

    ! declare inputs
    double precision,intent(in) :: deltat,deltaWP,lai,meant
    double precision,intent(inout) :: Rtot

    ! local variables
    integer :: i
    double precision :: slpa,soilR2,transpiration_resistance
    double precision, dimension(nos_root_layers) :: root_mass    &
                                                   ,root_length  &
                                                   ,ratio

    ! reset water flux
    water_flux = 0.0 ; ratio = 0.0 ; ratio(1) = 1.0
    ! calculate soil depth to which roots reach
    root_reach = max_depth * root_biomass / (root_k + root_biomass)
    ! calculate the plant hydraulic resistance component. Currently unclear
    ! whether this actually varies with height or whether tall trees have a
    ! xylem architecture which keeps the whole plant conductance (gplant) 1-10
    ! (ish).
    transpiration_resistance = 1.0 / ( gplant * lai )

    ! calculate layer thickness
    layer_thickness(1) = top_soil_depth ; layer_thickness(2)=max(0.1,root_reach-layer_thickness(1))
    layer_thickness(3) = max_depth - sum(layer_thickness(1:2))

    ! The original SPA src generates an exponential distribution which aims
    ! to maintain 50 % of root biomass in the top 25 % of the rooting depth.
    ! In a simple 2 root layer system this can be estimates more simply

    ! top 25 % of root profile
    slpa = (root_reach * 0.25) - layer_thickness(1)
    if (slpa <= 0.0) then
        ! > 50 % of root is in top layer
        root_mass(1) = root_biomass * 0.5
        root_mass(1) = root_mass(1) + ((root_biomass-root_mass(1)) * (abs(slpa)/root_reach))
        root_mass(2) = max(0.0,root_biomass - root_mass(1))
    else
        ! < 50 % of root is in bottom layer
        root_mass(1) = root_biomass * 0.5 * (layer_thickness(1)/(abs(slpa)+layer_thickness(1)))
        root_mass(2) = max(0.0,root_biomass - root_mass(1))
    endif
    root_length = root_mass / (root_density * root_cross_sec_area)
    !! Top root layer.
    ! soil conductivity converted from m.s-1 -> m2.s-1.MPa-1 by head
    soilR2=root_resistance(root_mass(1),min(root_reach,layer_thickness(1)))
    ! calculate and accumulate steady state water flux in mmol.m-2.s-1
    ! NOTE1: Depth correction already accounted for in soil resistance
    ! calculations and this is the maximum potential rate of transpiration
    ! assuming saturated soil and leaves at their minimum water potential.
    ! also note that the head correction is now added rather than
    ! subtracted in SPA equations because deltaWP is soilWP-minlwp not
    ! soilWP prior to application of minlwp
    ! NOTE2: this is a simplified version and does not include
    ! an estimate of soil hydraulic resistance as there is no BUCKET!
    demand = deltaWP+head*canopy_height
    water_flux(1) = max(0.,demand/(transpiration_resistance + soilR2))
    ratio(1)=1.0
    ! Bottom root layer
    if (root_mass(2) > 0.0 ) then
       ! soil conductivity converted from m.s-1 -> m2.s-1.MPa-1 by head
       soilR2=root_resistance(root_mass(2),layer_thickness(2))
       ! calculate and accumulate steady state water flux in mmol.m-2.s-1
       water_flux(2) = max(0.,demand/(transpiration_resistance + soilR2))
       ratio = layer_thickness(1:nos_root_layers)/sum(layer_thickness(1:nos_root_layers))
    endif ! roots present in second layer?

    ! calculate effective resistance
    Rtot = demand / sum(water_flux*ratio)

    ! and return
    return

  end subroutine calculate_Rtot
  !
  !------------------------------------------------------------------
  !
  double precision function root_resistance (root_mass,thickness)

   !
   ! Calculates root hydraulic resistance (MPa m2 s mmol-1) in a soil-root zone
   !

   implicit none

   ! arguments
   double precision :: root_mass, & ! root biomass in layer (gbiomass)
                       thickness    ! thickness of soil zone roots are in

   ! calculate root hydraulic resistance
   root_resistance = root_resist / (root_mass*thickness)
   ! return
   return

  end function root_resistance
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
