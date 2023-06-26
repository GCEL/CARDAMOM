
module CARBON_MODEL_MOD

implicit none

!!!!!!!!!!!
! Authorship contributions
!
! This code contains a variant of the Data Assimilation Linked ECosystem (DALEC) model.
! This version of DALEC is derived from that described in Bloom & Williams (2015), https://doi.org/10.5194/bg-12-1299-2015.
! This code is based on that created by A. A. Bloom (UoE, now at JPL, USA).
! Subsequent modifications by:
! T. L. Smallman (University of Edinburgh, t.l.smallman@ed.ac.uk)
! See function / subroutine specific comments for exceptions and contributors
!!!!!!!!!!!

! make all private
private

! explicit publics
public :: CARBON_MODEL  &
         ,soil_frac_clay   &
         ,soil_frac_sand   &
         ,nos_soil_layers  &
         ,dim_1,dim_2   &
         ,nos_trees     &
         ,nos_inputs    &
         ,leftDaughter  &
         ,rightDaughter &
         ,nodestatus    &
         ,xbestsplit    &
         ,nodepred      &
         ,bestvar

! ACM related parameters
double precision, parameter :: pi = 3.1415927
double precision, parameter :: deg_to_rad = pi/180.0

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

! for consisteny between requirements of different models
integer, parameter :: nos_root_layers = 2, nos_soil_layers = nos_root_layers + 1
double precision, dimension(nos_soil_layers) :: soil_frac_clay,soil_frac_sand

contains
!
!--------------------------------------------------------------------
!
  subroutine CARBON_MODEL(start,finish,met,pars,deltat,nodays,lat,lai,NEE,FLUXES,POOLS &
                       ,pft,nopars,nomet,nopools,nofluxes,GPP)

    ! The Data Assimilation Linked Ecosystem Carbon - Combined Deciduous
    ! Evergreen Analytical (DALEC_CDEA) model. The subroutine calls the
    ! Aggregated Canopy Model to simulate GPP and partitions between various
    ! ecosystem carbon pools. These pools are subject to turnovers /
    ! decompostion resulting in ecosystem phenology and fluxes of CO2

    implicit none

    ! declare input variables
    integer, intent(in) :: start    &
                          ,finish   &
                          ,nopars     & ! number of paremeters in vector
                          ,pft        & ! plant functional type
                          ,nomet      & ! number of meteorological fields
                          ,nofluxes   & ! number of model fluxes
                          ,nopools    & ! number of model pools
                          ,nodays       ! number of days in simulation

    double precision, intent(in) :: met(nomet,nodays)   & ! met drivers
                         ,deltat(nodays)                & ! time step in decimal days
                         ,pars(nopars)                  & ! number of parameters
                         ,lat                 ! site latitude (degrees)

    double precision, dimension(nodays), intent(inout) :: lai & ! leaf area index
                                               ,GPP & ! Gross primary productivity
                                               ,NEE   ! net ecosystem exchange of CO2

    double precision, dimension((nodays+1),nopools), intent(inout) :: POOLS ! vector of ecosystem pools

    double precision, dimension(nodays,nofluxes), intent(inout) :: FLUXES ! vector of ecosystem fluxes

    ! declare local variables
    double precision :: gpppars(12)            & ! ACM inputs (LAI+met)
             ,constants(10)          & ! parameters for ACM
             ,wf,wl,ff,fl,osf,osl,sf & ! phenological controls
             ,pi,ml

    integer :: p,f,nxp,n

    ! met drivers are:
    ! 1st run day
    ! 2nd min daily temp (oC)
    ! 3rd max daily temp (oC)
    ! 4th Radiation (MJ.m-2.day-1)
    ! 5th CO2 (ppm)
    ! 6th DOY
    ! 7th lagged cumulative precip

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
    ! 8 = labile release
    ! 9 = leaffall factor
    ! 10 = leaf litter production
    ! 11 = woodlitter production
    ! 12 = rootlitter production
    ! 13 = respiration het litter
    ! 14 = respiration het som
    ! 15 = litter2som
    ! 16 = labrelease factor

    ! PARAMETERS
    ! 16 values

    ! p(1) Litter to SOM conversion rate  - m_r
    ! p(2) Fraction of GPP respired - f_a
    ! p(3) Fraction of NPP allocated to foliage - f_f
    ! p(4) Fraction of NPP allocated to roots - f_r
    ! p(5) Leaf lifespan - L_f
    ! p(6) Turnover rate of wood - t_w
    ! p(7) Turnover rate of roots - t_r
    ! p(8) Litter turnover rate - t_l
    ! p(9) SOM turnover rate  - t_S
    ! p(10) Parameter in exponential term of temperature - \theta
    ! p(11) = date of Clab release - B_day
    ! p(12) = Fraction allocated to Clab - f_l
    ! p(13) = lab release duration period - R_l
    ! p(14) = date of leaf fall - F_day
    ! p(15) = leaf fall duration period - R_f
    ! p(16) = LMA

    ! set constants
    pi = 3.1415927

    ! load some values
    gpppars(4) = 1.89 ! foliar N gN.m-2 leaf area
    gpppars(7) = lat
    gpppars(11) = pi

    ! assign acm parameters
    if (pft == 1) then
       ! cropland
       constants(1)=12.16085    ; constants(2)=0.03023221
       constants(3)=45.87594    ; constants(4)=268.9338
       constants(5)=0.000747788 ; constants(6)=0.204327
       constants(7)=3.629274    ; constants(8)=0.01669404
       constants(9)=0.4573969   ; constants(10)=1.045534
       gpppars(9) = -1.987365   ; gpppars(10) = 0.6336006
    else if (pft == 2) then
       ! short grass
       constants(1)=12.16085    ; constants(2)=0.03023221
       constants(3)=45.87594    ; constants(4)=268.9338
       constants(5)=0.000747788 ; constants(6)=0.204327
       constants(7)=3.629274    ; constants(8)=0.01669404
       constants(9)=0.4573969   ; constants(10)=1.045534
       gpppars(9) = -1.987365   ; gpppars(10) = 0.6336006
    else if (pft == 3) then
       ! evergreen needle leaf
       constants(1)=10.0126    ; constants(2)=0.01623306
       constants(3)=53.89489    ; constants(4)=277.8853
       constants(5)=0.0004823053  ; constants(6)=0.3446748
       constants(7)=3.986002    ; constants(8)=0.01670374
       constants(9)=0.0875766 ; constants(10)=1.025008
       gpppars(9)= -1.969439    ; gpppars(10)=0.6700837
    else if (pft == 5) then
       ! Deciduous broadleaf
       constants(1)=6.257381e+00    ; constants(2)=1.915067e-02
       constants(3)=4.235703e+01    ; constants(4)=3.373044e+02
       constants(5)=8.916662e-02    ; constants(6)=2.465434e-01
       constants(7)=4.035877e+00    ; constants(8)=1.232785e-02
       constants(9)=8.697617e-01    ; constants(10)=5.374788e-01
       gpppars(9) = -1.912499e+00   ; gpppars(10) = 1.489004e+00
    else if (pft == 11) then
       ! Sparse vegetation
print*,"FAIL PFT SELECTION" ; stop
!       constants(1)=16.1247736227685    ; constants(2)=0.000617309178894364
!       constants(3)=89.512048631903     ; constants(4)=282.023083528009
!       constants(5)=0.0301297353658013  ; constants(6)=1.44748316513338
!       constants(7)=1.71596189896891    ; constants(8)=0.241803967406994
!       constants(9)=0.000306360169728313; constants(10)=-0.554150638896367
!       gpppars(9) = -2.0788751390956    ; gpppars(10) = 1.07337026956825
    else if (pft == 13) then
print*,"FAIL PFT SELECTION" ; stop
!       ! Bog
!       constants(1)=2.37323087141509    ; constants(2)=0.0433056473554607
!       constants(3)=96.0793202668459    ; constants(4)=524.190312195035
!       constants(5)=1.34018492473325e-06; constants(6)=1.69940113655973
!       constants(7)=1.41439798367338    ; constants(8)=0.0975902734200319
!       constants(9)=0.000954291392234617; constants(10)=-1.9108709274161
!       gpppars(9) = -2.07205814716648   ; gpppars(10) = 0.956166673229136
    else if (pft == 17) then
       ! C4 grass
       constants(1)=11.27848    ; constants(2)=0.0248564
       constants(3)=23.10066    ; constants(4)=301.1614
       constants(5)=0.0004963181; constants(6)=0.230574
       constants(7)=3.673286    ; constants(8)=0.01670022
       constants(9)=0.4284842   ; constants(10)=0.6791512
       gpppars(9) = -2.038791   ; gpppars(10) = 0.5130594
    else if (pft == 18) then
print*,"FAIL PFT SELECTION" ; stop
!       ! mixed forest
!       constants(1)=9.71571673798685    ; constants(2)=0.00808290472235316
!       constants(3)=99.3802919097048    ; constants(4)=596.29303517113
!       constants(5)=0.0269131568567349  ; constants(6)=1.93737880793734
!       constants(7)=7.23776549676329    ; constants(8)=0.136836539382728
!       constants(9)=0.00180138325980396 ; constants(10)=1.25991762411738
!       gpppars(9) = -1.93236861131181   ; gpppars(10) = 1.09477279931996
    else if (pft == 19) then
print*,"FAIL PFT SELECTION" ; stop
!       ! Interupted forest
!       constants(1)=13.8225676305661    ; constants(2)=0.0100104941805291
!       constants(3)=98.7337385590505    ; constants(4)=589.032594902438
!       constants(5)=0.00246376099756583 ; constants(6)=1.97405597521383
!       constants(7)=8.02661955979173    ; constants(8)=0.102049963328083
!       constants(9)=3.97172806458295e-05; constants(10)=1.07602065071723
!       gpppars(9) = -2.07307691061897   ; gpppars(10) = 1.09610425681178
    else
       print*,"PFT is missing or misidentified pft input = ", pft
    end if

    if (start == 1) then
       ! assigning initial conditions
       POOLS(1,1)=pars(17)
       POOLS(1,2)=pars(18)
       POOLS(1,3)=pars(19)
       POOLS(1,4)=pars(20)
       POOLS(1,5)=pars(21)
       POOLS(1,6)=pars(22)
    endif

    ! defining phenological variables
    ! release period coefficient, based on duration of labile turnover or leaf
    ! fall durations
    wf=pars(15)*sqrt(2.0)/2.0
    wl=pars(13)*sqrt(2.0)/2.0

    ! magnitude coefficient
    ff=(log(pars(5))-log(pars(5)-1.0))/2.0
    fl=(log(1.001)-log(0.001))/2.0

    ! set minium labile life span to one year
    ml=1.001

    ! offset for labile and leaf turnovers
    osf=ospolynomial(pars(5),wf)
    osl=ospolynomial(ml,wl)

    ! scaling to biyearly sine curve
    sf=365.25/pi

    !
    ! Begin looping through each time step
    !

    do n = start, finish

      ! calculate LAI value
      lai(n)=POOLS(n,2)/pars(16)

      ! load next met / lai values for ACM
      gpppars(1)=lai(n)
      gpppars(2)=met(3,n) ! max temp
      gpppars(3)=met(2,n) ! min temp
      gpppars(5)=met(5,n) ! co2
      gpppars(6)=ceiling(met(6,n)-(deltat(n)*0.5)) ! doy
      gpppars(8)=met(4,n) ! radiation
!      gpppars(12)=met(3,n)-met(2,n) ! daily temperature range
!      gpppars(12)=met(7,n) ! lagged precip

      if (lai(n) > 1e-6) then
          ! GPP (gC.m-2.day-1)
          FLUXES(n,1) = acm(gpppars,constants)
!          call randomForest(1,gpppars,FLUXES(n,1))
!           if (FLUXES(n,1) /= FLUXES(n,1)) then
!              print*,"oooh dear we have some NaN values appearing in GPP =",FLUXES(n,1)
!              print*,"met in step ", met(:,n)
!              print*,"pools in step", POOLS(n,:)
!              stop
!           endif
      else
          FLUXES(n,1) = 0.0
      end if
      ! temprate (i.e. temperature modified rate of metabolic activity))
      FLUXES(n,2) = exp(pars(10)*0.5*(met(3,n)+met(2,n)))
      ! autotrophic respiration (gC.m-2.day-1)
      FLUXES(n,3) = pars(2)*FLUXES(n,1)
      ! leaf production rate (gC.m-2.day-1)
      FLUXES(n,4) = (FLUXES(n,1)-FLUXES(n,3))*pars(3)
      ! labile production (gC.m-2.day-1)
      FLUXES(n,5) = (FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,4))*pars(12)
      ! root production (gC.m-2.day-1)
      FLUXES(n,6) = (FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,4)-FLUXES(n,5))*pars(4)
      ! wood production
      FLUXES(n,7) = FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,4)-FLUXES(n,5)-FLUXES(n,6)

      ! Labile release and leaffall factors
      FLUXES(n,9) = (2d0/(pi**0.5))*(ff/wf)*exp(-((sin((met(1,n)-pars(14)+osf)/sf)*sf/wf)**2d0))
      FLUXES(n,16) = (2d0/(pi**0.5))*(fl/wl)*exp(-((sin((met(1,n)-pars(11)+osl)/sf)*sf/wl)**2d0))

      !
      ! those with time dependancies
      !

      ! total labile release
      FLUXES(n,8) = POOLS(n,1)*(1d0-(1d0-FLUXES(n,16))**deltat(n))/deltat(n)
      ! total leaf litter production
      FLUXES(n,10) = POOLS(n,2)*(1d0-(1d0-FLUXES(n,9))**deltat(n))/deltat(n)
      ! total wood production
      FLUXES(n,11) = POOLS(n,4)*(1d0-(1d0-pars(6))**deltat(n))/deltat(n)
      ! total root litter production
      FLUXES(n,12) = POOLS(n,3)*(1d0-(1d0-pars(7))**deltat(n))/deltat(n)

      !
      ! those with temperature AND time dependancies
      !

      ! respiration heterotrophic litter
      FLUXES(n,13) = POOLS(n,5)*(1d0-(1d0-FLUXES(n,2)*pars(8))**deltat(n))/deltat(n)
      ! respiration heterotrophic som
      FLUXES(n,14) = POOLS(n,6)*(1d0-(1d0-FLUXES(n,2)*pars(9))**deltat(n))/deltat(n)
      ! litter to som
      FLUXES(n,15) = POOLS(n,5)*(1d0-(1d0-pars(1)*FLUXES(n,2))**deltat(n))/deltat(n)

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
      POOLS(n+1,3) = POOLS(n,3) + (FLUXES(n,6) - FLUXES(n,12))*deltat(n)
      ! litter pool
      POOLS(n+1,5) = POOLS(n,5) + (FLUXES(n,10)+FLUXES(n,12)-FLUXES(n,13)-FLUXES(n,15))*deltat(n)
      ! som pool
      POOLS(n+1,6) = POOLS(n,6) + (FLUXES(n,15)-FLUXES(n,14)+FLUXES(n,11))*deltat(n)

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
    ! daily canopy conductance (m.s-1) ; NOTE that ratio of H20:CO2 diffusion is
    ! 1.646259 (Jones appendix 2). i.e. gcCO2/1.646259 = gcH2O
!    gc=abs(deltaWP)**(hydraulic_exponent)/((hydraulic_temp_coef*Rtot+trange))
    gc=abs(deltaWP)**(hydraulic_exponent)/((hydraulic_temp_coef*Rtot))
    ! maximum rate of temperature and nitrogen (canopy efficiency) limited photosynthesis (gC.m-2.day-1)
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
  subroutine randomForest(nos_wanted,drivers,ans)

    ! Random Forest regression based emulator using R v2.13 randomForest
    ! function library(randomForest).

    implicit none

    ! declare inputs
    integer, intent(in) :: nos_wanted !
    double precision, intent(in) :: drivers(12)

    ! declare output variables
    double precision, intent(out), dimension(nos_wanted) :: ans

    ! declare local variables
    integer :: mdim,ntest
    double precision :: lai,maxt,mint,co2,doy,radiation &
                       ,lat,pi,dec,aob,sinld,cosld &
                       ,new_data(nos_wanted,nos_inputs),dayl

    ! This function will calculate the GPP based on randomForests generated in the randomForest library
    ! The function has been simplied from the original, this however means that we are dependent on the user
    ! correctly presenting the new data in rows (each point) and columns (meteorological drivers)
    ! in order of "lai","maxt","mint","RAD","dayl","CO2","lagged_precip" for
    ! CTESSEL and "mint", "maxt", "RAD", "CO2", "avgN", "lai", "Trange" for
    ! JULES

    ! load inputs to appropriate local variables
    doy = drivers(6)
    lat = drivers(7)
    pi = drivers(11)

    ! for JULES
    ! low temp (oC)
    new_data(1,1) = drivers(3)
    ! max temp (oC)
    new_data(1,2) = drivers(2)
    ! SW radiation (MJ.m-2.day-1)
    new_data(1,3) = drivers(8)
    ! co2 (ppm)
    new_data(1,4) = drivers(5)
    ! avgN (gN.m-2 leaf area)
    new_data(1,6) = drivers(4)
    ! lai (m2/m2)
    new_data(1,7) = drivers(1)
    ! Trange (oC)
    new_data(1,8) = drivers(12)

    ! for CTESSEL
!    ! lai (m2/m2)
!    new_data(1,1) = drivers(1)
!    ! low temp (oC)
!    new_data(1,3) = drivers(3)
!    ! max temp (oC)
!    new_data(1,2) = drivers(2)
!    ! SW radiation (MJ.m-2.day-1)
!    new_data(1,4) = drivers(8)
!    ! co2 (ppm)
!    new_data(1,6) = drivers(5)
!    ! lagged precip (mm)
!    new_data(1,7) = drivers(12)

    ! calculate day length (hours)
    dec = - asin( sin( 23.45 * pi / 180.0 ) * cos( 2.0 * pi * ( doy + 10.0 ) /365.0 ) )
    sinld = sin( lat*(pi/180.0) ) * sin( dec )
    cosld = cos( lat*(pi/180.0) ) * cos( dec )
    aob = max(-1.0,min(1.0,sinld / cosld))
    dayl = 12.0 * ( 1.0 + 2.0 * asin( aob ) / pi )

    ! finally load daylength (hrs)
    new_data(1,5) = dayl

    ! extract dimensional information about the model
    mdim = nos_inputs ! number of drivers
    ntest = nos_wanted ! number of predictions requested

    call regForest(transpose(new_data),mdim,ntest,ans)

    ! return
    return

  end subroutine randomForest
  !
  !------------------------------------------------------------------
  !
  subroutine regForest(x,mdim,n,ypred)

    implicit none

    ! declare inputs
    integer, intent(in) :: mdim  & ! number of met inputs
                          ,n       ! number of outputs requested

    double precision, intent(in) :: x(mdim,n) !

    ! define output variable
    double precision, intent(out), dimension(n) :: ypred

    ! define local variables
    double precision, dimension(n):: ytree
    integer :: idx1, i, j, z, k, m

    ! initial conditions
    idx1 = 1 ; ypred = 0d0

    ! looks like we run each tree for each location first then move
    ! onto the next tree and keep adding things up
    do i = 1, nos_trees
        ! zero this instance
        ytree = 0d0
        ! key output from predictRegTree appears to be ytree
        call predictRegTree(x,n,mdim,idx1,ytree)
        ! add this trees solution to the given requested values
        ypred = ypred + ytree
        !/* increment the offset */
        idx1 = idx1 + 1
    end do

    ! return variable of interest
    ypred = (ypred/real(nos_trees)) ; return

  end subroutine regForest
  !
  !------------------------------------------------------------------
  !
  subroutine predictRegTree(x,nsample,mdim,idx1,ypred)

    ! function moves through a given tree (idx1)

    implicit none

    ! declare inputs
    integer, intent(in) :: mdim,idx1,nsample
    double precision, intent(in) :: x (mdim,nsample) !

    ! output variables
    double precision, intent(out) :: ypred(nsample)
    ! local variables
    integer :: i, k, m

    ypred=1d0
    do i = 1, nsample
        k = 1
        !/* go down the tree */
        do while (nodestatus(k,idx1) /= -1)
            m = bestvar(k,idx1) !- 1
            if (x(m,i) <= xbestsplit(k,idx1)) then
                k = leftDaughter(k,idx1) !- 1
            else
                k = rightDaughter(k,idx1) !- 1
            end if
        end do ! for while loop
        !/* terminal node: assign prediction and move on to next */
        ypred(i) = nodepred(k,idx1)
    end do ! end nsample loop

    ! return output
    return

  end subroutine predictRegTree
  !
  !------------------------------------------------------------------
  !
  double precision function ospolynomial(L,w)

    ! Function calculates the day offset for Labile release and leaf turnover
    ! functions

    implicit none

    ! declare input variables
    double precision, intent(in) ::  L, w ! polynomial coefficients and scaling factor

    ! declare local variables
    double precision ::  LLog, mxc(7) ! polynomial coefficients and scaling factor

    ! assign polynomial terms
    mxc(1)=(0.000023599784710)
    mxc(2)=(0.000332730053021)
    mxc(3)=(0.000901865258885)
    mxc(4)=(-0.005437736864888)
    mxc(5)=(-0.020836027517787)
    mxc(6)=(0.126972018064287)
    mxc(7)=(-0.188459767342504)

    ! load log of leaf / labile turnovers
    LLog=log(L-1d0)

    ! calculate the polynomial function
    ospolynomial=(mxc(1)*LLog**6d0 + mxc(2)*LLog**5d0 + &
                  mxc(3)*LLog**4d0 + mxc(4)*LLog**3d0 + &
                  mxc(5)*LLog**2d0 + mxc(6)*LLog     + mxc(7))*w

  end function ospolynomial
!
!--------------------------------------------------------------------
!
end module CARBON_MODEL_MOD
