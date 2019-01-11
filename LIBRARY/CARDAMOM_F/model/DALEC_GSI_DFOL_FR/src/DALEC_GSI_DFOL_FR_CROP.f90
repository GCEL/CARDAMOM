
module CARBON_MODEL_CROP_MOD

implicit none

! make all private
private

! explicit publics
public :: CARBON_MODEL_CROP    &
         ,resp_rate_temp_coeff &
         ,ts_length            &
         ,sec_in_day           &
         ,sec_in_hour

! declare some module level variables
double precision, parameter :: sec_in_day  = 86400.0
double precision, parameter :: sec_in_hour = 3600.0
double precision, parameter :: pi          = 3.14159265
double precision, parameter :: deg_to_rad  = pi/180d0

! variables local to this module..
integer ::   plough_day, & ! day-of-year when field is ploughed  
                sow_day, & ! day-of-year when field is sown      
            harvest_day, & ! day-of-year when field is harvested 
              stmob = 0, & ! remoblise stem C to labile (1 = on)
 turnover_labile_switch    ! begin turnover of labile C

logical :: vernal_calcs, &  ! do vernalisation calculations?
               ploughed, &  !
        use_seed_labile, &  != .False. ! whether to use seed labile for growth
                   sown, &  != .False. ! has farmer sown crop yet?
                emerged     != .False. ! has crop emerged yet?

double precision ::           ts_length, & ! time step length in hours
                            step_of_day, & ! current step of the day (default = 1)
                           steps_in_day, & ! number of steps in a day (default = 1)
                                    doy, & ! decimal doy of year 
                                gpp_acm, & ! gross primary productivity (gC.m-2.day-1)
                    stock_storage_organ, & ! storage organ C pool, i.e. the desired crop (gC.m--2)
                     stock_dead_foliage, & ! dead but still standing foliage (gC.m--2)
                        stock_resp_auto, & ! autotrophic respiration pool (gC.m--2)
                           stock_labile, & ! labile C pool (gC.m--2)
                          stock_foliage, & ! foliage C pool (gC.m--2)
                             stock_stem, & ! stem C pool (gC.m--2)
                            stock_roots, & ! roots C pool (gC.m--2)
                           stock_litter, & ! litter C pool (gC.m--2)
                    stock_soilOrgMatter, & ! SOM C pool (gC.m--2)
                              resp_auto, & ! autotrophic respiration (gC.m-2.t-1)
                          resp_h_litter, & ! litter heterotrophic respiration (gC.m-2.t-1)
                   resp_h_soilOrgMatter, & ! SOM heterotrophic respiration (gC.m-2)
                                    npp, & ! net primary productivity (gC.m-2.t-1)
                              nee_dalec, & ! net ecosystem exchange (gC.m-2.t-1)
                              lai_local, & !
                                     DS, & ! Developmental state and initial condition
                                    LMA, & ! leaf mass area (gC.m-2)
            mean_alloc_to_storage_organ, & ! rolling average allocation of GPP to storage organ (gC.m-2)
        mean_alloc_to_storage_organ_old, & ! ...same but previous value...
                     decomposition_rate, & ! decomposition rate (frac / hr)
                     frac_GPP_resp_auto, & ! fraction of GPP allocated to autotrophic carbon pool
                  turnover_rate_foliage, & ! turnover rate of foliage (frac/hr)
                     turnover_rate_stem, & ! same for stem
                   turnover_rate_labile, & ! same for labile 
                turnover_rate_resp_auto, & ! same for autotrophic C pool
                 resp_cost_labile_trans, & ! labile lost to respiration per gC labile to GPP
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
                                    VDh, & ! effective vernalisation days when plants are 50 % vernalised 
                                     VD, & ! count of vernalisation days
                               RDRSHMAX, & ! maximum rate of self shading turnover
                                   PHCR, & ! critical value of photoperiod for development
                                   PHSC, & ! photoperiod sensitivity
                             raso = 0.0, & ! rolling average for alloc to storage organ
                         max_raso = 0.0, & ! maximum value for rolling average alloc to storage organ
                                  BM_EX, & ! 
                                     HI, & !
                                  yield, & ! crop yield (gC.m-2)
                     alloc_to_resp_auto, & ! amount of carbon to allocate to autotrophic respiration pool
                    turnover_rate_roots, & ! turnover over rate of roots interpolated each time step
                                gso_max, & !
                     max_raso_old = 0.0, & !
                        raso_old  = 0.0, & !
            resp_cost_labile_to_foliage, & ! respiratory cost of moving carbon..from labile to foliage pools
            resp_cost_foliage_to_labile, & ! ..from foliage to labile pools
                              resp_rate, & ! rate of respiration at given temperature
                                 Cshoot, & !
                                     DR, & !
                        fol_frac_intpol, & !
                       stem_frac_intpol, & !
                               fP,fT,fV, & !                      
                                  remob, & !
                       root_frac_intpol, & !
                      shoot_frac_intpol, & !
                                 avtemp, & !
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
                                raremob, & !
                                  RDRSH, & !
                                  RDRDV, & !
                                    RDR, & !
                              daylength, & !
                              leaf_CN_snap


  ! defines Q10 = 2 in exponential temperature response for heterotrophic respiration
  double precision, parameter :: resp_rate_temp_coeff = 0.0693
  ! residue fraction of leaves left post harvest
  double precision, parameter :: lv_res = 0.1
  ! residue fraction of stem left post harvest
  double precision, parameter :: st_res = 0.1
  ! LAI above which self shading turnover occurs
  double precision, parameter :: LAICR = 4.0
  ! allocation to storage organ relative to GPP 
  double precision, parameter :: rel_gso_max = 0.35

  save ! indicates all these module variables should be held, even if the module  goes out of scope

contains
!
!--------------------------------------------------------------------
!
  subroutine CARBON_MODEL_CROP(start,finish,met,pars,deltat,nodays,lat,lai,NEE,FLUXES,POOLS &
                       ,pft,nopars,nomet,nopools,nofluxes,GPP,stock_seed_labile)

    ! The Data Assimilation Linked Ecosystem Carbon (DALEC) model. 
    ! The subroutine calls the Aggregated Canopy Model to simulate GPP 
    ! and partitions between various ecosystem carbon pools. 
    ! These pools are subject to turnovers/decompostion resulting in 
    ! ecosystem phenology and fluxes of CO2

    use DALEC_CROP_DEV_VARIABLES ! Now need to define the variables here

    implicit none

    ! declare input variables
    integer, intent(in) :: start    & 
                          ,finish   & 
                          ,nopars   & ! number of paremeters in vector
                          ,pft      & ! plant functional type
                          ,nomet    & ! number of meteorological fields
                          ,nofluxes & ! number of model fluxes
                          ,nopools  & ! number of model pools
                          ,nodays     ! number of days in simulation

    double precision, intent(in) :: met(nomet,nodays) & ! met drivers
                         ,stock_seed_labile           & ! seed carbon to get things going
                         ,deltat(nodays)              & ! time step in decimal days
                         ,pars(nopars)                & ! number of parameters
                         ,lat                           ! site latitude (degrees)

    ! double precision, dimension(:), intent(inout) ::   DS_shoot, & !
    !                                                     DS_root, & !
    !                                                    fol_frac, & !
    !                                                   stem_frac, & !
    !                                                   root_frac, & !
    !                                                     DS_LRLV, & ! 
    !                                                        LRLV, & !
    !                                                     DS_LRRT, & !
    !                                                        LRRT    ! 

    double precision, dimension(nodays), intent(inout) :: lai & ! leaf area index
                                                         ,GPP & ! Gross primary productivity
                                                         ,NEE   ! net ecosystem exchange of CO2

    double precision, dimension((nodays+1),nopools), intent(inout) :: POOLS ! vector of ecosystem pools
 
    double precision, dimension(nodays,nofluxes), intent(inout)    :: FLUXES ! vector of ecosystem fluxes
 
    ! declare local variables
    double precision :: airt_weighting(3) &
                       ,gpppars(12)            & ! ACM inputs (LAI+met)
                       ,constants(10)          & ! parameters for ACM
                       ,declin,sinld,cosld,aob

    integer :: p,f,nxp,n
    character(20) :: t

    ! Allocate lai array
    double precision, allocatable, dimension(:) :: leaf_N, DS_crop_N, leaf_CN

    ! MET drivers :
    ! 1 = run day
    ! 2 = min daily temp (oC)
    ! 3 = max daily temp (oC)
    ! 4 = Radiation (MJ.m-2.day-1)
    ! 5 = CO2 (ppm)
    ! 6 = DOY

    ! POOLS :
    ! 1 = labile
    ! 2 = foliar
    ! 3 = root
    ! 4 = wood/stem
    ! 5 = litter
    ! 6 = som
    ! 7 = autotrophic
    ! 8 = storage organ C

    ! FLUXES : 
    ! 1 = GPP
    ! 2 = temprate
    ! 3 = respiration_auto
    ! 4 = leaf production
    ! 5 = labile production
    ! 6 = root production
    ! 7 = wood production
    ! 8 = labile release
    ! 9 = alloc to storage
    ! 10 = leaf litter production
    ! 11 = woodlitter production
    ! 12 = rootlitter production
    ! 13 = respiration het litter
    ! 14 = respiration het som
    ! 15 = litter2som (decomposition)
    ! 16 = alloc to autotrophic pool

    ! PARAMETERS
    ! p(1)  = decomposition rate (frac/hr)
    ! p(2)  = Fraction of GPP respired
    ! p(3)  = max development rate (day-1) for DS (0->1)
    ! p(4)  = max development rate (day-1) for DS (1->2)
    ! p(5)  = turnover rate of foliage (frac/hr)
    ! p(6)  = Turnover rate of wood/stem (frac/hr)
    ! p(7)  = maximum rate of foliar turnover due to self shading
    ! p(8)  = effective vernalisation days when plant is 50 % vernalised
    ! p(9)  = mineralisation rate of som
    ! p(10) = mineralisation rate of litter
    ! p(11) = log10(avgN) 
    ! p(12) = sow day
    ! p(13) = labile lost to respiration per gC labile top GPP
    ! p(14) = phenological heat units needed for emergence
    ! p(15) = harvest day (doy)
    ! p(16) = plough day (doy)
    ! p(17) = leaf mass area (gC.m-2)
    ! p(18)-P(25) = intial values for 
    ! labile, foliar, roots, stem, litter, som, autotrophic, storage_organ pools
    ! p(26) ! min T for development
    ! p(27) ! max T for development
    ! p(28) ! optimum T for development
    ! p(29) ! min T for vernalisation
    ! p(30) ! max T for vernalisation
    ! p(31) ! optimim T for vernalisation
    ! p(32) ! critical value of photoperiod for development
    ! p(33) ! photoperiod sensitivity
    ! p(34) ! turnover rate of labile C
    ! p(35) ! turnover rate of autotrophic C

    ! load some values
    gpppars(7)  = lat
    gpppars(9)  = -2.060814 ! leafWP-soilWP
    gpppars(10) = 1.0       ! totaly hydraulic resistance
    gpppars(11) = pi

    ! assign acm parameters: cropland-specific 
    constants(1)  = 10.38 
    constants(2)  = 0.04 
    constants(3)  = 2.70e-4 
    constants(4)  = 83.18 
    constants(5)  = 0.03
    constants(6)  = 4.54  
    constants(7)  = 3.86 
    constants(8)  = 4.10e-3
    constants(9)  = 0.38  
    constants(10) = 2.72e-8


    ! -- Define allocation fractions (from development.csv SPA-Crop file) ---------------- !
    
    ! First allocate arrays [DS_shoot,DS_root, DS_LRLV(loss rate leaves), DS_LRRT(loss rate roots)]  
    if (.not. allocated(DS_shoot)) allocate( DS_shoot(10) , fol_frac(10) , stem_frac(10)  )  ! Shoot array
    if (.not. allocated(DS_root)) allocate( DS_root(5) , root_frac(5) ) ! Root array
    if (.not. allocated(DS_LRLV)) allocate( DS_LRLV(5) , LRLV(5) ) ! Loss rate leaves
    if (.not. allocated(DS_LRRT)) allocate( DS_LRRT(5) , LRRT(5) ) ! Loss rate roots
    if (.not. allocated(DS_crop_N)) allocate( DS_crop_N(6) , leaf_CN(6) ) ! leaf C:N ratio

    !stock_seed_labile = 9.0 ! this is passed from the MAIN_DALEC_CROP.f90 file instead
    
    DS_shoot  = (/ 0.0,0.33,0.43,0.53,0.62,0.77,0.95,1.14,1.38,2.1 /)
    fol_frac  = (/ 0.9,0.85,0.83,0.75,0.56,0.2,0.09,0.05,0.0,0.0 /)
    stem_frac = (/ 0.1,0.15,0.17,0.25,0.44,0.8,0.64,0.62,0.0,0.0 /)

    DS_root   = (/ 0.0,0.33,0.53,1.0,2.1 /)   
    root_frac = (/ 0.5,0.5,0.25,0.0,0.0 /)

    DS_root   = (/ 0.0,0.33,0.53,1.0,2.1 /)   
    root_frac = (/ 0.5,0.5,0.25,0.0,0.0 /)

    DS_crop_N = (/ 0.0, 0.47, 0.63,  1.015, 1.815, 2.1 /)   
    leaf_CN   = (/ 0.0, 8.13, 13.76, 14.51, 22.71, 8.13 /)

    DS_LRLV   = (/ 0.0,1.0,1.3,1.8,2.5 /)
    LRLV      = (/ 0.0,0.0,0.025,0.025,0.025 /)

    DS_LRRT   = (/ 0.0,1.0,1.3,1.8,2.5 /) 
    LRRT      = (/ 0.0,0.0,0.00625,0.00625,0.00625 /)        

    ! Create an array of daily foliar N values ---------------------- !
    if (.not. allocated(leaf_N)) allocate( leaf_N(nodays))

    ! set the same value for all elements in array
    leaf_N(1:nodays) = pars(11)

    ! --------------------------------------------------------------- !

    ! length of time step in hours..
    ts_length = ((sum(deltat)/nodays) * sec_in_day) / sec_in_hour
    
    ! temporal resolution of ins/outs
    ! steps_in_day = 1d0 ; step_of_day = 1d0
    steps_in_day = 1 / (sum(deltat)/nodays) ! number of steps in a day 
    step_of_day  = 1 / (sum(deltat)/nodays) ! current step of the day 

    ! parameters 
    decomposition_rate                = pars(1)  ! decomposition rate (frac / hr)
    frac_GPP_resp_auto                = pars(2)  ! fraction of GPP allocated to autotrophic carbon pool
    DR_pre                            = pars(3)  ! development rate coefficient DS (0->1)
    DR_post                           = pars(4)  ! development rate coefficient DS (1->2)
    turnover_rate_foliage             = pars(5)  ! turnover_rate of foliage (frac/hr)
    turnover_rate_stem                = pars(6)  ! turnover rate of stem (frac/hr)
    RDRSHMAX                          = pars(7)  ! maximum rate of foliar turnover due to self shading
    VDh                               = pars(8)  ! effective vernalisation days when plants are 50 % vernalised 
    mineralisation_rate_soilOrgMatter = pars(9)  ! mineralisation rate som
    mineralisation_rate_litter        = pars(10) ! mineralisation rate litter
    resp_cost_labile_trans            = pars(13) ! labile lost to respiration per gC labile to GPP
    PHUem                             = pars(14) ! phenological heat units required for emergence
    LMA                               = pars(17) ! leaf mass area (gC.m-2)
    tmin                              = pars(26) ! min temperature for development
    tmax                              = pars(27) ! max temperature for development
    topt                              = pars(28) ! optimum temperature for development
    tmin_v                            = pars(29) ! min temperature for vernalisation
    tmax_v                            = pars(30) ! max temperature for vernalisation
    topt_v                            = pars(31) ! optimim temperature for vernalisation
    PHCR                              = pars(32) ! critical value of photoperiod for development
    PHSC                              = pars(33) ! photoperiod sensitivity
    turnover_rate_labile              = pars(34) ! turnover rate labile C
    turnover_rate_resp_auto           = pars(35) ! turnover rate of autotrophic carbon for respiration
    
    ! (could know) management dates :
    plough_day                        = nint(mod(pars(12)-2d0,365.25)) ! nint(mod(pars(16),365.25)) ! plough day (doy)
    sow_day                           = nint(mod(pars(12),365.25)) ! sow day (doy)
    harvest_day                       = nint(mod(pars(15),365.25)) ! nint(mod(pars(15),365.25)) ! harvest day (doy)


    if (start == 1) then
        ! pair incoming variables to local module levels
        ! finally set some initial conditions
        yield = 0.0
        DS = -1.0 
        mean_alloc_to_storage_organ = 0.0
        mean_alloc_to_storage_organ_old = 0.0
        PHU = 0.0
        VD = 0.0
        BM_EX = 0.0
        HI = 0.0
        ploughed = .false.
        vernal_calcs = .true. 
        sown = .true.
        use_seed_labile = .true.  
        emerged = .false.
        stock_dead_foliage=0.0
        DR = 0.0
        alloc_to_labile = 0.0
        stmob = 0.0
        max_raso = 0.0
        DR = 0.0
        raso = 0.0
        RDRDV = 0.0 

        ! assigning initial conditions
        POOLS(1,1)=pars(18)
        POOLS(1,2)=pars(19)
        POOLS(1,3)=pars(20)
        POOLS(1,4)=pars(21)
        POOLS(1,5)=pars(22)
        POOLS(1,6)=pars(23)
        POOLS(1,7)=pars(24) 
        POOLS(1,8)=pars(25)
 
        stock_labile          = pars(18) ! labile C
        stock_foliage         = pars(19) ! foliar C
        stock_roots           = pars(20) ! root C
        stock_stem            = pars(21) ! stem / wood C
        stock_litter          = pars(22) ! litter C
        stock_soilOrgMatter   = pars(23) ! som C
        stock_resp_auto       = pars(24) ! autotrophic resp pool
        stock_storage_organ   = pars(25) ! storage organ (i.e. desired crop)

    endif ! start == 1

    ! 
    ! Begin looping through each time step
    ! 

    do n = start, finish
 
      ! Leaf N calculation 
      if ( (DS > 0.) .and. (DS <= 1.5) ) then 
        leaf_CN_snap = interpolate( DS , DS_crop_N , leaf_CN , size(DS_crop_N) )
        leaf_N(n)  = (leaf_CN_snap  / POOLS(n,2)) ** (-1)
      endif 
      if  (DS > 1.5) then                
          leaf_N(n) = leaf_N(n-1)
      endif

      ! calculate LAI value
      lai(n)    = POOLS(n,2) / LMA
      lai_local = lai(n)

      ! load next met / lai values for ACM
      gpppars(1) = lai(n)
      gpppars(2) = met(3,n)                          ! max T
      gpppars(3) = met(2,n)                          ! min T
      gpppars(4) = leaf_N(n)                         ! leaf g.N.m-2
      gpppars(5) = met(5,n)                          ! co2
      gpppars(6) = ceiling(met(6,n)-(deltat(n)*0.5)) ! doy
      gpppars(8) = met(4,n)                          ! radiation

      ! allocate doy of year
      doy = met(6,n)

      ! calculate new GPP value (gC.m-2.day-1)
      if (lai(n) > 1e-6) then
          gpp_acm = acm(gpppars,constants)
      else
          gpp_acm = 0.0
      end if

      ! pass relevant variables into crop module memory
      avtemp =  0.5 * ( met(3,n) + met(2,n) )

      ! Heterotrophic respiration rate (Q10):  doubles with 
      ! 10 degree temperature rise resprate from soil file = 0.0693
      resp_rate = 0.5 * exp( resp_rate_temp_coeff * avtemp )

      ! daily avg of allocation to storage organ (needed to determine max storage organ growth rate)
      mean_alloc_to_storage_organ_old = mean_alloc_to_storage_organ
      mean_alloc_to_storage_organ     = 0.0

      ! determine development stage (DS)
      call development_stage
      ! determine the carbon partitioning based on development stage
      call carbon_alloc_fractions(DS_shoot,DS_root,fol_frac,stem_frac,root_frac)
      ! begin carbon allocation for crops
      call calc_pools_crops(DS_LRRT,LRRT)
      ! conduct management updates at the end of the day
      call management_dates(stock_seed_labile,steps_in_day)

      ! calculate the NEE 
      NEE(n) = nee_dalec
      ! load GPP
      GPP(n) = gpp_acm

      ! GPP (gC.m-2.day-1)
      FLUXES(n,1) = gpp_acm
      ! temprate (i.e. temperature modified rate of metabolic activity))
      FLUXES(n,2) = resp_rate
      ! autotrophic respiration (gC.m-2.day-1)
      FLUXES(n,3) = resp_auto
      ! leaf production rate (gC.m-2.day-1)
      FLUXES(n,4) = alloc_to_foliage
      ! labile production (gC.m-2.day-1)
      FLUXES(n,5) = alloc_to_labile + remob
      ! root production (gC.m-2.day-1)
      FLUXES(n,6) = alloc_to_roots
      ! wood production 
      FLUXES(n,7) = alloc_to_stem
      ! alloc to storage organ
      FLUXES(n,9) = alloc_to_storage_organ
      ! alloc to autotrophic pool
      FLUXES(n,16) = frac_GPP_resp_auto * gpp_acm
      ! labile production
      FLUXES(n,8) = alloc_from_labile + resp_cost_labile_to_foliage
      ! total leaf litter production
      FLUXES(n,10) = litterfall_foliage
      ! total wood production
      FLUXES(n,11) = litterfall_stem
      ! total root litter production
      FLUXES(n,12) = litterfall_roots
      ! respiration heterotrophic litter
      FLUXES(n,13) = resp_h_litter
      ! respiration heterotrophic som
      FLUXES(n,14) = resp_h_soilOrgMatter
      ! litter to som
      FLUXES(n,15) = decomposition

      ! labile pool
      POOLS(n+1,1) = stock_labile
      ! foliar pool
      POOLS(n+1,2) = stock_foliage
      ! wood pool
      POOLS(n+1,4) = stock_stem
      ! root pool
      POOLS(n+1,3) = stock_roots
      ! litter pool
      POOLS(n+1,5) = stock_litter
      ! som pool
      POOLS(n+1,6) = stock_soilOrgMatter
      ! autotrophic pool 
      POOLS(n+1,7) = stock_resp_auto
      ! storage organ pool
      POOLS(n+1,8) = stock_storage_organ

    end do ! no days loop

  end subroutine CARBON_MODEL_CROP
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
    double precision :: gc, pn, pd, pp, qq, ci, e0, cps, dec, nit &
             ,trange, sinld, cosld,aob,pi, mult &
             ,mint,maxt,radiation,co2,lat &
             ,deltaWP,Rtot,NUE,temp_exponent,dayl_coef &
             ,dayl_const,hydraulic_exponent,hydraulic_temp_coef &
             ,co2_comp_point,co2_half_sat,lai_coef,lai_const

    ! initial values
    gc=0.0 ; pp=0.0 ; qq=0.0 ; ci=0.0 ; e0=0.0 ; cps=0.0 ; dec=0.0 ; nit=1.

    ! load driver values to correct local vars
    maxt = drivers(2)
    mint = drivers(3)
    nit = drivers(4)   
    co2 = drivers(5)
    radiation = drivers(8)
    lat = drivers(7)

    ! load parameters into correct local vars
    pi = drivers(11)
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
    ! daily canopy conductance 
    gc=abs(deltaWP)**(hydraulic_exponent)/((hydraulic_temp_coef*Rtot+trange)) ! default
    ! maximum rate of temperature and nitrogen (canopy efficiency) limited
    ! photosynthesis (gC.m-2.day-1)
    pn=lai_local*nit*NUE*exp(temp_exponent*maxt) ! default
    ! pp and qq represent limitation by diffusion and metabolites respecitively
    pp=pn/gc ; qq=co2_comp_point-co2_half_sat
    ! calculate internal CO2 concentration (ppm)
    ci=0.5*(co2+qq-pp+((co2+qq-pp)**2-4.0*(co2*qq-pp*co2_comp_point))**0.5)
    ! limit maximum quantium efficiency by leaf area, hyperbola
    e0=lai_coef*lai_local**2/(lai_local**2+lai_const)
    ! calculate day length (hours)
    dec = - asin( sin( 23.45 * deg_to_rad ) * cos( 2.0 * pi * ( doy + 10.0 ) / 365.0 ) )
    sinld = sin( lat*deg_to_rad ) * sin( dec )
    cosld = cos( lat*deg_to_rad ) * cos( dec )
    aob = max(-1.0,min(1.0,sinld / cosld))
    daylength = 12.0 * ( 1.0 + 2.0 * asin( aob ) / pi )
    ! calculate CO2 limited rate of photosynthesis
    pd=gc*(co2-ci)
    ! calculate combined light and CO2 limited photosynthesis
    cps=e0*radiation*pd/(e0*radiation+pd)
    ! correct for day length variation
    acm=cps*(dayl_coef*daylength+dayl_const)

    ! don't forget to return
    return

  end function acm
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine calc_pools_crops(DS_LRRT,LRRT)

    ! Allocated GPP to NPP and various carbon pools. Based !
    ! this on physiological responses to temperature       !
    ! vernalisation, and photoperiod.                      !

    implicit none

    ! arguments
    double precision, dimension(:), intent(inout) :: DS_LRRT, & !
                                                        LRRT    !

    ! turnover rate of fine roots is now equal to the loss rate of roots (Penning de Vries, 1989)
    turnover_rate_roots = interpolate( DS , DS_LRRT , LRRT , 5 ) / 24.0

    ! if sown turn on labile / seed turnover for growth
    if ( sown ) then
      ! turnover on
      turnover_labile_switch = 1
    else
      ! turnover off
      turnover_labile_switch = 0
    endif

    ! Initialise..
    resp_cost_foliage_to_labile = 0.0

    ! respiratory cost of C transfer from labile pool to short-term pool (NPP)
    resp_cost_labile_to_foliage = turnover_rate_labile * stock_labile * resp_cost_labile_trans * resp_rate &
                                   * ts_length * real(turnover_labile_switch)

    ! allocation flux from labile C pool to NPP
    alloc_from_labile = turnover_rate_labile * stock_labile * ( 1.0 - resp_cost_labile_trans ) * resp_rate &
                                  * ts_length * real(turnover_labile_switch)

    ! When GPP is higher than seed C content, remaining seed carbon enters litter
    ! C pool, as seedlings do not fully exhaust their seed (P. de Vries p 48)
    if ( ( gpp_acm .gt. alloc_from_labile ) .and. ( use_seed_labile ) ) then
      stock_litter = stock_litter + stock_labile
      stock_labile = 0.0
      use_seed_labile = .false.
    endif

    ! NPP as a fraction of GPP (1-.32=.68 or 68%) + allocation..
    npp = ( 1.0 - frac_GPP_resp_auto ) * gpp_acm + alloc_from_labile
    ! from labile pool; = SHORT-TERM POOL

    root_frac_intpol  = max(0d0,min(1d0,root_frac_intpol))
    alloc_to_roots    = root_frac_intpol * npp         !
    npp_shoot         = npp - alloc_to_roots           ! NPP remaining after root growth = SHOOT frac
    alloc_to_foliage  = fol_frac_intpol  * npp_shoot   !
    alloc_to_stem     = stem_frac_intpol * npp_shoot   !
    alloc_to_storage_organ = max(0d0,npp_shoot - alloc_to_foliage - alloc_to_stem)
    
    if ( alloc_to_storage_organ > 0.0 ) then  ! allocation flux to storage organ limited by maximum growth rate
        gso_max  = ( stock_storage_organ + 0.5 ) * rel_gso_max / steps_in_day
        alloc_to_storage_organ = min( alloc_to_storage_organ , gso_max )
        
        if ( sown ) then
          alloc_to_labile = ( npp_shoot - alloc_to_foliage - alloc_to_stem - alloc_to_storage_organ ) &
                              * ( 1.0 - resp_cost_labile_trans )
          resp_cost_foliage_to_labile =  ( npp_shoot - alloc_to_foliage - alloc_to_stem - alloc_to_storage_organ ) &
                              * resp_cost_labile_trans
        else
          alloc_to_labile             = 0.0
          resp_cost_foliage_to_labile = 0.0
        endif
    
    endif
    
    mean_alloc_to_storage_organ = mean_alloc_to_storage_organ + alloc_to_storage_organ

    ! set switches to (de)activate leaf, root and stem remobliization
    if ( step_of_day .eq. steps_in_day ) then
      mean_alloc_to_storage_organ = mean_alloc_to_storage_organ / steps_in_day
      raso_old = raso
      ! running average of growth rate of storage organ..
      raso = ( mean_alloc_to_storage_organ + mean_alloc_to_storage_organ_old ) * 0.5
      max_raso_old = max_raso
      max_raso = max( raso , max_raso_old )
      ! Stem remobilisation triggered once running average of storage organ growth declines
      ! Second part prevents premature remobilisation
      if ( ( raso .lt. raso_old ) .and. &
            ( mean_alloc_to_storage_organ .gt. ( mean_alloc_to_storage_organ_old + 0.5 ) / steps_in_day ) ) then
          stmob = 1
      else 
          stmob = 0
      endif
    endif

    ! Code for calculating relative death rate of leaves (RDR) as a 
    !  function of shading (RDRSH) or developmental stage (RDRT).

    ! GT 0 if LAI GT 4; 0. < RDRSH < RDRSHMAX (usually ~0.03)
    RDRSH = min( RDRSHMAX , max( 0d0 , RDRSHMAX * ( lai_local - LAICR ) / LAICR ) )
    if ( DS < 1.0 ) then
       RDRDV = 0.0
    else
       RDRDV = turnover_rate_foliage * ( 1 / ( ( max( 2.0 - DS , 0.1 ) ) * 8.0 ) ) ** 2
    ENDIF

    ! relative leaf death rate is the maximum value of the arguments RDRSH and
    ! RDRDV
    RDR = max( RDRSH , RDRDV )

    ! remobilization of foliar C and allocation to dead leaves pool
    litterfall_foliage = stock_foliage * ts_length * RDR
    litterfall_stem    = stock_stem    * ts_length * DR * turnover_rate_stem * real(stmob) ! remobstem
    litterfall_roots   = stock_roots   * ts_length * turnover_rate_roots

    ! remobilized C to NPP (from both leaves and stems)
    remob   = ( litterfall_foliage * 0.5 + litterfall_stem ) * ( 1.0 - resp_cost_labile_trans )
    ! respiratory cost of C transfer (conversion from starch to photosynthates)
    Raremob = ( litterfall_foliage * 0.5 + litterfall_stem ) * resp_cost_labile_trans

    ! heterotrophic respiration component 1: mineralisation of litter C pool
    resp_h_litter = mineralisation_rate_litter * stock_litter * resp_rate * ts_length
    ! heterotrophic respiration component 2:  mineralisation of organic matter C pool
    resp_h_soilOrgMatter = mineralisation_rate_soilOrgMatter * stock_soilOrgMatter * resp_rate * ts_length
    ! decomposition of litter to soil organic matter
    decomposition = decomposition_rate * stock_litter * resp_rate * ts_length

    ! Recalculate Carbon Pools...

    stock_foliage       = max(0d0, stock_foliage + alloc_to_foliage - litterfall_foliage)
    stock_stem          = max(0d0, stock_stem + alloc_to_stem - litterfall_stem)
    stock_storage_organ = max(0d0, stock_storage_organ + alloc_to_storage_organ)
    stock_roots         = max(0d0, stock_roots + alloc_to_roots   - litterfall_roots)
    stock_litter        = max(0d0, stock_litter + litterfall_roots - resp_h_litter - decomposition)
    stock_soilOrgMatter = max(0d0, stock_soilOrgMatter + decomposition - resp_h_soilOrgMatter)
    stock_dead_foliage  = max(0d0, stock_dead_foliage  + litterfall_foliage * 0.5)
    stock_labile        = max(0d0, stock_labile + alloc_to_labile - alloc_from_labile - resp_cost_labile_to_foliage + remob)

    ! respiratory pool: new photosynthates are added 
    stock_resp_auto = stock_resp_auto + frac_GPP_resp_auto * gpp_acm
    ! autotrophic respiration; Ra (7% of respiratory pool)
    resp_auto = stock_resp_auto * turnover_rate_resp_auto * ts_length
    ! respiratory pool reduced by Ra (amount of C respired by plant)
    stock_resp_auto = max(0d0, stock_resp_auto - resp_auto)
    ! respiratory cost of C transfer from labile pool to short-term pool added
    ! to yield total autotrophic respiration
    resp_auto = resp_auto + resp_cost_labile_to_foliage + resp_cost_foliage_to_labile + Raremob
    ! nee (gC.m-2.t-1)
    nee_dalec = (resp_auto + resp_h_litter + resp_h_soilOrgMatter) - gpp_acm

  end subroutine calc_pools_crops
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine carbon_alloc_fractions(DS_shoot,DS_root,fol_frac,stem_frac,root_frac)

    ! Determines carbon allocation fractions as a function !
    ! of developmental stage (DS).  Allocation fractions   !
    ! are from tables published in Penning de Vries (1989) !

    implicit none

    double precision, dimension(:), intent(inout) ::   DS_shoot, & !
                                                        DS_root, & !
                                                       fol_frac, & !
                                                      stem_frac, & !
                                                      root_frac    !

    ! local variables..
    integer                       :: organ
    double precision,dimension(:),allocatable :: frac_shoot, frac_root

    if ( sown ) then ! after sowing

       allocate( frac_shoot(size(DS_shoot)) , frac_root(size(DS_root)) )

       ! loop over three crop "organs": 1) foliage 2) stems 3) root
       ! not necessary for storage organs, as all remaining C is allocated to
       ! these use different input for foliage and stem fractions, as they are
       ! relative to the total shoot (or aboveground) allocation, root is relative to
       ! total plant (above- and belowground) allocation..

       ! leaf development stages and corresponding fractions..
       frac_shoot = fol_frac
       ! interpolate between PdV allocation values with reference to developmental stage 
       fol_frac_intpol = interpolate( DS , DS_shoot , frac_shoot , size(DS_shoot) )
       ! stem DS and fracs
       frac_shoot = stem_frac
       stem_frac_intpol = interpolate( DS , DS_shoot , frac_shoot , size(DS_shoot) )
       ! root DS and fracs..
       frac_root = root_frac
       root_frac_intpol = interpolate( DS , DS_root , frac_root , size(DS_root) )

    endif

  end subroutine carbon_alloc_fractions
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine development_stage()

    ! Based on modified Wang & Engel model (Streck et al., 2003), !
    ! but with only 2 sub-phases, vegetative and reproductive     !
    ! (i.e. only two different DRmax).   O. Sus, May 2010.        !

    implicit none

    ! local variables..
    double precision ::  doptmin, & ! Difference between optimum and minimum temperature
                         dmaxmin, & ! Difference between maximum and minimum temperature
                          dttmin, & ! Difference between daiy average and minimum temperatures
                       doptmin_v, & ! Difference between optimum and minimum vernalization temperatures
                       dmaxmin_v, & ! Difference between maximum and minimum vernalization temperatures
                       dttmin_v     ! Difference between daily average and minimum vernalization temperatures

    doptmin   = topt   - tmin   ! difference between optimal and minimum cardinal temperatures
    dmaxmin   = tmax   - tmin   ! difference between maximum and minimum cardinal temperatures
    dttmin    = avtemp - tmin   ! difference between daily average and minimum cardinal temperatures
    doptmin_v = topt_v - tmin_v ! same as above,
    dmaxmin_v = tmax_v - tmin_v !       but for vernalization 
    dttmin_v  = avtemp - tmin_v ! cardinal temperatures

    ! Calculation of developmental function values: vernalization (fV),temperature (fT) and photoperiod (fP) 
    ! these values are multiplicative factors of DRmax (maximum developmental rate)
    ! each ranging between 0 (no development) and 1 (unrestricted development)

    ! Summation of vernalization days (VD), not before sowing and only if avg T is within min and max cardinal T
    if ( ( avtemp .gt. tmin_v ) .and. ( avtemp .lt. tmax_v ) .and. sown ) then
      fV = vernalization( doptmin_v , dmaxmin_v , dttmin_v )
    endif

    ! Only calculate temperature coefficient if avtemp lies within (tmin,tmax) range
    if ( (avtemp .gt. tmin ) .and. ( avtemp .lt. tmax ) ) then
      fT = temperature_impact( doptmin , dmaxmin , dttmin )
    else
      fT = 0.0
    endif

    fP = photoperiod_impact( PHCR , PHSC ) ! calculation of photoperiod coefficient

    if ( emerged .and. ( DS .lt. 2.0 ) ) then   ! sum up daily DR values between emergence and maturity (DS=2)

       if ( DS .lt. 1.0 ) then  ! in the vegetative phase (before flowering):

          DR = DR_pre * fT * fP   ! DR is affected by temperature, photoperiod...

          if ( vernal_calcs ) DR = DR * fV ! ...and vernalization (for winter cereals)
          DS = DS + DR    ! developmental stage (DS), calculated as the sum of daily developmental rates

       else    ! in the reproductive phase (after flowering):

          DR = DR_post * fT   ! DR is affected only by temperature
          DS = DS + DR

       endif

    endif

  end subroutine development_stage
  
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  
  subroutine management_dates (stock_seed_labile,steps_in_day)

    ! use CARBON_MODEL_CROP_MOD

    ! This routine should be called at the end of each day of a crops  !
    ! simulation.  It checks whether we should plough/sow/harvest, and !
    ! during the growing establishes when the crop will emerge after   !
    ! sowing, based on heat accumulation (Phenological Heat Units).    !

    implicit none

    ! arguments
    double precision, intent(in) :: stock_seed_labile, steps_in_day
    ! local variables
    double precision :: tmp
    logical :: plough_sanity, sow_sanity, harvest_sanity

    ! reset
    plough_sanity = .false. ; sow_sanity = .false. ; harvest_sanity = .false.

    ! spring crops
    ! if (sow_day < harvest_day .and. nint(doy) < harvest_day) sow_sanity = .true.
    ! if (plough_day < harvest_day .and. nint(doy) < harvest_day) plough_sanity = .true.
    ! if (harvest_day > sow_day) harvest_sanity = .true.
    ! winter crops
    if (sow_day > harvest_day) sow_sanity = .true. 
    if (plough_day > harvest_day) plough_sanity = .true.
    if (harvest_day < plough_day .and. nint(doy) < plough_day) harvest_sanity = .true.

    if ( .not. sown ) then ! fresh field...

      if ( plough_sanity .and. .not.ploughed .and. nint(doy) >= plough_day ) then
        ! the field needs ploughing..
        call plough 

      elseif ( sow_sanity .and. nint(doy) >= sow_day ) then

        ! ensure that the field has indeed been ploughed
        if (.not.ploughed) then 
            call plough
        endif

        ! the field needs sowing..
        sown = .true.

        ! this switch controls whether the labile C within the seed is used for growth
        use_seed_labile = .true.
        stock_labile = stock_seed_labile

      endif ! plough or sow?

    else

      ! crop in field... calculate when crop emerges
      if ( .not. emerged ) then

         ! estimate emergence date based on the accumulated phenological heat units (PHU)
         ! where PHU is the (positive) heat over tmin..
         tmp = max( avtemp - tmin , 0.0 ) * steps_in_day
         PHU = PHU + tmp

         ! set the development stage and emergence..
         if ( PHU >= PHUem ) then
           emerged = .true.
           DS = 0.0
         else
           emerged = .false.
           DS = -1.0
         endif

      endif ! emerged or not

      ! note that in this case harvest day has been fixed relative to the sow day
      if ( harvest_sanity .and. nint(doy) >= harvest_day) then ! the field needs harvesting
         call harvest
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

  subroutine harvest()

    implicit none

    ! shoot biomass..
    Cshoot = stock_foliage + stock_stem + stock_storage_organ + stock_labile

    ! determine harvest index..
    HI = stock_storage_organ / Cshoot

    ! the stuff we actually want from the harvest...
    yield  = stock_storage_organ

    ! the biomass that is harvested in addition to the storage-organ..
    BM_EX  = stock_foliage * ( 1.0 - lv_res )          &
              + stock_stem * ( 1.0 - st_res )          &
               + stock_dead_foliage * ( 1.0 - lv_res ) &
                + stock_labile

    ! what's left (will fall to the ground)..
    stock_litter  = stock_litter                     &
                    + stock_resp_auto                &
                     + stock_foliage * lv_res        &
                      + stock_stem * st_res          &
                       + stock_dead_foliage * lv_res

    ! empty the plant stocks..
    stock_storage_organ = 0.0
    stock_foliage       = 0.0
    stock_stem          = 0.0
    stock_dead_foliage  = 0.0
    stock_labile        = 0.0
    stock_resp_auto     = 0.0

    ! roots stay in ground and slowly decompose (until/unless the field is ploughed)

    ! reset logical variables..
    sown    = .False.
    emerged = .False.
    fV = 0.0 ; fT = 0.0 ; fP = 0.0
    VD = 0.0 ; DS = -9999

  end subroutine harvest

  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !

  ! double precision function calculate_daylength_hours(lat,doy_in)

  !   use CARBON_MODEL_MOD, only: dble_zero, dble_one

  !     ! this function uses latitude (degrees) and julian day of year to calculate day length in hours

  !     ! arguments
  !     double precision :: doy_in, lat

  !     ! local variables
  !     double precision :: dec, sinld, cosld, aob

  !     dec = - asin( sin( 23.45 * deg_to_rad ) * cos( 2.0 * pi * ( doy + 10.0 ) / 365.0 ) )
  !     sinld = sin( lat*deg_to_rad ) * sin( dec )
  !     cosld = cos( lat*deg_to_rad ) * cos( dec )
  !     aob = max(-dble_one,min(dble_one,sinld / cosld))
  !     calculate_daylength_hours = 12.0 * ( 1.0 + 2.0 * asin( aob ) / pi )

  !     return

  ! end function calculate_daylength_hours

  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !

  double precision function photoperiod_impact( PH_crit , PH_sens )

    ! Function to determine the coefficient for !
    ! photoperiod impact on developmental rate. !
    ! From Streck et al., 2003                  !

    implicit none

    ! arguments..
    double precision,intent(in) :: PH_crit, & ! critical photoperiod below which no development occurs
                                   PH_sens    ! photoperiod sensitivity

    photoperiod_impact = max(0d0, 1d0 - exp ( - PH_Sens * ( daylength - PH_crit ) ))

  end function photoperiod_impact
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine plough()

    ! this s/r will reset various carbon pools, to mimic the effect of the
    ! farmer ploughing. !

    implicit none

    ! Move all plant stocks into the litter pool.
    ! ( many of these should already be empty after the harvest, )
    ! ( e.g. the stocks for labile, foliage, storage-organ stem. )
    stock_litter        = stock_litter + stock_dead_foliage &
                          + stock_foliage + stock_labile    &
                           + stock_roots + stock_stem       &
                            + stock_storage_organ

    stock_dead_foliage  = 0.0
    stock_foliage       = 0.0
    stock_labile        = 0.0
    stock_roots         = 0.0
    stock_stem          = 0.0
    stock_storage_organ = 0.0

  end subroutine plough
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
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

    a   = log( 2.0 ) / ( log( ( dmaxmin ) / doptmin ) )
    nmr = 2.0 * ( ( dttmin ) ** a ) * ( doptmin ** a ) - ( ( dttmin ) ** ( 2.0 * a ) )
    dnr = doptmin ** ( 2.0 * a )
    temperature_impact = nmr / dnr

  end function temperature_impact
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  double precision function vernalization( doptmin_v , dmaxmin_v , dttmin_v )

    ! Function to determine the coefficent for vernalization !
    ! impact on developmental rate. See Streck et al., 2003. !

    implicit none

    ! arguments..
    double precision,intent(in) :: dmaxmin_v , doptmin_v , dttmin_v ! temperature differences

    ! local variables..
    double precision :: a , dnr , fvn , nmr

    a   = log( 2.0 ) / ( log( ( dmaxmin_v ) / doptmin_v ) )
    nmr = 2.0 * ( ( dttmin_v ) ** a ) * ( doptmin_v ** a ) - ( ( dttmin_v ) ** (2.0 * a ) )
    dnr = doptmin_v ** ( 2.0 * a )
    fvn = nmr / dnr
    VD = VD + fvn
    ! final output value : 
    vernalization = max( 0d0 , min( 1d0 , ( VD ** 5.0 ) / ( ( VDh ** 5.0 ) + (VD ** 5.0 ) ) ) )

  end function vernalization
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  double precision function interpolate( x , reference_x , reference_y , row )

    ! Interpolation function.                    !
    ! x is input value, interpol is output value !
    ! reference_x/y are reference input data.    !

    implicit none

    ! arguments..
    integer,intent(in)                         :: row
    double precision,intent(in)                :: x
    double precision,dimension(row),intent(in) :: reference_x , reference_y

    ! local variables..
    integer::i

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

  end function interpolate
!
!--------------------------------------------------------------------
!
end module CARBON_MODEL_CROP_MOD
