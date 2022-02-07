
module CARBON_MODEL_MOD

implicit none

  !!!!!!!!!!!
  ! Authorship contributions
  !
  ! This code contains a variant of the Data Assimilation Linked ECosystem (DALEC) model.
  ! This version of DALEC is derived from the following primary references:
  ! Bloom & Williams (2015), https://doi.org/10.5194/bg-12-1299-2015.
  ! This code is based on that created by A. A. Bloom (UoE, now at JPL, USA).
  ! Subsequent modifications by:
  ! T. L. Smallman (University of Edinburgh, t.l.smallman@ed.ac.uk)
  ! See function / subroutine specific comments for exceptions and contributors
  !!!!!!!!!!!

! make all private
private

! explicit publics
public :: CARBON_MODEL     &
         ,soil_frac_clay   &
         ,soil_frac_sand   &
         ,nos_soil_layers  &
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

! forest rotation specific info
double precision, allocatable, dimension(:) :: extracted_C

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
                       ,nopars,nomet,nopools,nofluxes,GPP)

    ! The Data Assimilation Linked Ecosystem Carbon - Combined Deciduous
    ! Evergreen Analytical (DALEC_CDEA) model. The subroutine calls the
    ! Aggregated Canopy Model to simulate GPP and partitions between various
    ! ecosystem carbon pools. These pools are subject to turnovers /
    ! decompostion resulting in ecosystem phenology and fluxes of CO2

    ! This version includes the option to simulate fire combustion based
    ! on burned fraction and fixed combusion rates. It also includes the
    ! possibility to remove a fraction of biomass to simulate deforestation.

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

    ! declare local variables
    double precision :: gpppars(12)            & ! ACM inputs (LAI+met)
             ,constants(10)          & ! parameters for ACM
             ,wf,wl,ff,fl,osf,osl,sf & ! phenological controls
             ,pi,ml, doy
    ! JFE added 4 May 2018 - combustion efficiencies and fire resilience
    double precision :: cf(6),rfac

    integer :: p,f,n,ii ! JFE added ii to loop over fluxes

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
    ! 3 = root+wood
    ! 4 = litter + som

    ! FLUXES are:
    ! 1 = GPP
    ! 2 = temprate
    ! 3 = respiration_auto
    ! 4 = leaf production
    ! 5 = labile production
    ! 6 = wood+root production
    ! 7 =
    ! 8 = labile production
    ! 9 = leaffall factor
    ! 10 = leaf litter production
    ! 11 = wood+root litter production
    ! 12 =
    ! 13 = respiration het litter+som
    ! 14 =
    ! 15 =
    ! 16 = labrelease factor

    ! JFE added 03/05/2018 - start
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
    ! JFE added 03/05/2018 - stop


    ! PARAMETERS
    ! 14 values + 4 initial conditions

    ! p(1) Fraction of GPP respired
    ! p(2) Fraction of NPP allocated to foliage
    ! p(3) Leaf lifespan
    ! p(4) Turnover rate of wood+root
    ! p(5) Litter+som turnover rate
    ! p(6) Parameter in exponential term of temperature
    ! p(7) Canopy efficiency parameter
    ! p(8) doy of Clab release
    ! p(9) Fraction allocated to Clab
    ! p(10) lab release duration period
    ! p(11) date of leaf fall
    ! p(12) leaf fall duration period
    ! p(13) LMA (gC/m2)

    ! p(14) Initial labile
    ! p(15) Initial foliage
    ! p(16) Initial wood+roots
    ! p(17) Initial litter + som

    ! set constants
    pi = 3.1415927d0

    ! load some values
    gpppars(4) = 1d0 ! foliar N
    gpppars(7) = lat
    gpppars(9) = -2d0 ! leafWP-soilWP
    gpppars(10) = 1d0 ! totaly hydraulic resistance
    gpppars(11) = pi

    ! assign acm parameters
    constants(1) = pars(7)
    constants(2) = 0.0156935d0
    constants(3) = 4.22273d0
    constants(4) = 208.868d0
    constants(5) = 0.0453194d0
    constants(6) = 0.37836d0
    constants(7) = 7.19298d0
    constants(8) = 0.011136d0
    constants(9) = 2.1001d0
    constants(10) = 0.789798d0

    if (start == 1) then
       ! assigning initial conditions
       POOLS(1,1) = pars(14) ! labile
       POOLS(1,2) = pars(15) ! foliar
       POOLS(1,3) = pars(16) ! roots+wood
       POOLS(1,4) = pars(17) ! litter + som
    endif
    ! defining phenological variables
    ! release period coefficient, based on duration of labile turnover or leaf
    ! fall durations
    wf = pars(12)*sqrt(2d0) * 0.5d0
    wl = pars(10)*sqrt(2d0) * 0.5d0

    ! magnitude coefficient
    ff = (log(pars(3))-log(pars(3)-1d0)) * 0.5d0
    fl = (log(1.001d0)-log(0.001d0)) * 0.5d0

    ! set minium labile life span to one year
    ml = 1.001d0

    ! offset for labile and leaf turnovers
    osf = ospolynomial(pars(3),wf)
    osl = ospolynomial(ml,wl)

    ! scaling to biyearly sine curve
    sf = 365.25d0/pi

    ! JFE added 4 May 2018 - define fire constants
    cf(1) = 0.1d0         ! labile combustion efficiency
    cf(2) = 0.9d0         ! foliar combustion efficiency
    cf(3) = 0.1d0         ! roots combustion efficiency
    cf(4) = 0.1d0         ! wood combustion efficiency
    cf(5) = 0.5d0         ! litter combustion efficiency
    cf(6) = 0.01d0        ! som combustion efficency

    rfac = 0.5d0          ! resilience factor

    !
    ! Begin looping through each time step
    !

    do n = start, finish

      ! calculate LAI value
      lai(n) = POOLS(n,2)/pars(13)

      ! estimate multiple use variable
      doy = met(6,n)-(deltat(n)*0.5d0) ! doy

      ! load next met / lai values for ACM
      gpppars(1) = lai(n)
      gpppars(2) = met(3,n) ! max temp
      gpppars(3) = met(2,n) ! min temp
      gpppars(5) = met(5,n) ! co2
      gpppars(6) = doy
      gpppars(8) = met(4,n) ! radiation

      ! GPP (gC.m-2.day-1)
      FLUXES(n,1) = acm(gpppars,constants)
      ! temprate (i.e. temperature modified rate of metabolic activity))
      FLUXES(n,2) = exp(pars(6)*0.5d0*(met(3,n)+met(2,n)))
      ! autotrophic respiration (gC.m-2.day-1)
      FLUXES(n,3) = pars(1)*FLUXES(n,1)
      ! leaf production rate (gC.m-2.day-1)
      FLUXES(n,4) = (FLUXES(n,1)-FLUXES(n,3))*pars(2)
      ! labile production (gC.m-2.day-1)
      FLUXES(n,5) = (FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,4))*pars(9)
      ! wood + root production (gC.m-2.day-1)
      FLUXES(n,6) = (FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,4)-FLUXES(n,5))

      ! Labile release and leaffall factors
      FLUXES(n,9) = (2d0/sqrt(pi))*(ff/wf)*exp(-(sin((doy-pars(11)+osf)/sf)*sf/wf)**2d0)
      FLUXES(n,16) = (2d0/sqrt(pi))*(fl/wl)*exp(-(sin((doy-pars(8)+osl)/sf)*sf/wl)**2d0)

      !
      ! those with time dependancies
      !

      ! total labile release
      FLUXES(n,8) = POOLS(n,1)*(1d0-(1d0-FLUXES(n,16))**deltat(n))/deltat(n)
      ! total leaf litter production
      FLUXES(n,10) = POOLS(n,2)*(1d0-(1d0-FLUXES(n,9))**deltat(n))/deltat(n)
      ! total root+wood production
      FLUXES(n,11) = POOLS(n,3)*(1d0-(1d0-pars(4))**deltat(n))/deltat(n)

      !
      ! those with temperature AND time dependancies
      !

      ! respiration heterotrophic litter
      FLUXES(n,13) = POOLS(n,4)*(1d0-(1d0-FLUXES(n,2)*pars(5))**deltat(n))/deltat(n)

      ! calculate the NEE
      NEE(n) = -FLUXES(n,1)+FLUXES(n,3)+FLUXES(n,13)
      ! load GPP
      GPP(n) = FLUXES(n,1)

      !
      ! update pools for next timestep
      !

      ! labile pool
      POOLS(n+1,1) = POOLS(n,1) + (FLUXES(n,5)-FLUXES(n,8))*deltat(n)
      ! foliar pool
      POOLS(n+1,2) = POOLS(n,2) + (FLUXES(n,4)+FLUXES(n,8)-FLUXES(n,10))*deltat(n)
      ! root + wood pool
      POOLS(n+1,3) = POOLS(n,3) + (FLUXES(n,6)-FLUXES(n,11))*deltat(n)
      ! litter + som pool
      POOLS(n+1,4) = POOLS(n,4) + (FLUXES(n,10)+FLUXES(n,11)-FLUXES(n,13))*deltat(n)


      ! JFE added 4 May 2018 - remove biomass if necessary
      if (met(8,n) > 0d0) then
          POOLS(n+1,1) = POOLS(n+1,1)*(1d0-met(8,n)) ! remove labile
          POOLS(n+1,2) = POOLS(n+1,2)*(1d0-met(8,n)) ! remove foliar
          POOLS(n+1,3) = POOLS(n+1,3)*(1d0-met(8,n)) ! remove wood
      end if

      ! calculate fire emissions and litter transfer
      if (met(9,n) > 0d0) then
          ! first calculate combustion / emissions fluxes in g C m-2 d-1
          FLUXES(n,18) = POOLS(n+1,1)*met(9,n)*cf(1)/deltat(n) ! labile
          FLUXES(n,19) = POOLS(n+1,2)*met(9,n)*cf(2)/deltat(n) ! foliar
          FLUXES(n,20) = POOLS(n+1,3)*met(9,n)*cf(3)/deltat(n) ! roots+wood
          FLUXES(n,21) = POOLS(n+1,4)*met(9,n)*(sum(cf(5:6))*0.5d0)/deltat(n)

          ! second calculate litter transfer fluxes in g C m-2 d-1, all pools except som
          FLUXES(n,24) = POOLS(n+1,1)*met(9,n)*(1d0-cf(1))*(1d0-rfac)/deltat(n) ! labile into litter
          FLUXES(n,25) = POOLS(n+1,2)*met(9,n)*(1d0-cf(2))*(1d0-rfac)/deltat(n) ! foliar into litter
          FLUXES(n,26) = POOLS(n+1,3)*met(9,n)*(1d0-cf(3))*(1d0-rfac)/deltat(n) ! wood+roots into litter

          ! update pools - first update vegetation
          POOLS(n+1,1) = POOLS(n+1,1) - (FLUXES(n,18) + FLUXES(n,24)) * deltat(n) ! labile
          POOLS(n+1,2) = POOLS(n+1,2) - (FLUXES(n,19) + FLUXES(n,25)) * deltat(n) ! foliar
          POOLS(n+1,3) = POOLS(n+1,3) - (FLUXES(n,20) + FLUXES(n,26)) * deltat(n) ! wood+roots
          ! update pools - second dead organic matter
          POOLS(n+1,4) = POOLS(n+1,4) + (FLUXES(n,24) + FLUXES(n,25) + FLUXES(n,26) - FLUXES(n,21)) * deltat(n)

          ! calculate ecosystem emissions
          FLUXES(n,17) = FLUXES(n,18)+FLUXES(n,19)+FLUXES(n,20)+FLUXES(n,21)
      else
          ! set fluxes to zero
          FLUXES(n,17) = 0d0
          FLUXES(n,18) = 0d0
          FLUXES(n,19) = 0d0
          FLUXES(n,20) = 0d0
          FLUXES(n,21) = 0d0
          FLUXES(n,22) = 0d0
          FLUXES(n,23) = 0d0
          FLUXES(n,24) = 0d0
          FLUXES(n,25) = 0d0
          FLUXES(n,26) = 0d0
          FLUXES(n,27) = 0d0
          FLUXES(n,28) = 0d0
      end if

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
             ,trange, sinld, cosld,aob,pi, mult &
             ,mint,maxt,radiation,co2,lai,doy,lat &
             ,deltaWP,Rtot,NUE,temp_exponent,dayl_coef &
             ,dayl_const,hydraulic_exponent,hydraulic_temp_coef &
             ,co2_comp_point,co2_half_sat,lai_coef,lai_const

    ! initial values
    gc = 0d0 ; pp = 0d0 ; qq = 0d0 ; ci = 0d0 ; e0 = 0d0 ; dayl = 0d0 ; cps = 0d0 ; dec = 0d0 ; nit = 1d0

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
    trange = 0.5d0*(maxt-mint)
    ! daily canopy conductance
    gc = abs(deltaWP)**(hydraulic_exponent)/((hydraulic_temp_coef*Rtot+trange))
    ! maximum rate of temperature and nitrogen (canopy efficiency) limited photosynthesis (gC.m-2.day-1)
    pn = lai*nit*NUE*exp(temp_exponent*maxt)
    ! pp and qq represent limitation by diffusion and metabolites respectively
    pp = pn/gc ; qq = co2_comp_point-co2_half_sat
    ! calculate internal CO2 concentration (ppm)
    ci = 0.5d0*(co2+qq-pp+((co2+qq-pp)**2d0-4d0*(co2*qq-pp*co2_comp_point))**0.5d0)
    ! limit maximum quantum efficiency by leaf area, hyperbola
    e0 = lai_coef*lai**2d0/(lai**2d0+lai_const)
    ! calculate day length (hours)
!    dec = - asin( sin( 23.45d0 * pi / 180d0 ) * cos( 2d0 * pi * ( doy + 10d0 ) /365d0 ) )
!    sinld = sin( lat*(pi/180d0) ) * sin( dec )
!    cosld = cos( lat*(pi/180d0) ) * cos( dec )
!    aob = max(-1d0,min(1d0,sinld / cosld))
!    dayl = 12d0 * ( 1d0 + 2d0 * asin( aob ) / pi )

!--------------------------------------------------------------
!    ! calculate day length (hours - not really hours)
!    ! This is the old REFLEX project calculation but it is wrong so anyway here
!    ! we go...
    dec=-23.4*cos((360.0*(doy+10.0)/365.0)*pi/180.0)*pi/180.0
    mult=tan(lat*pi/180.0)*tan(dec)
    if (mult>=1.0) then
      dayl=24.0
    else if (mult<=-1.0) then
      dayl=0.0
    else
      dayl=24.0*acos(-mult)/pi
    end if
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
  double precision function ospolynomial(L,w)

    ! Function calculates the day offset for Labile release and leaf turnover
    ! functions

    implicit none

    ! declare input variables
    double precision, intent(in) ::  L, w ! polynomial coefficients and scaling factor

    ! declare local variables
    double precision ::  tmp, LLog, mxc(7) ! polynomial coefficients and scaling factor

    ! assign polynomial terms
    mxc(1) = (0.000023599784710d0)
    mxc(2) = (0.000332730053021d0)
    mxc(3) = (0.000901865258885d0)
    mxc(4) = (-0.005437736864888d0)
    mxc(5) = (-0.020836027517787d0)
    mxc(6) = (0.126972018064287d0)
    mxc(7) = (-0.188459767342504d0)

    ! load log of leaf / labile turnovers
    LLog = log(L-1d0)

    ! calculate the polynomial function
    ospolynomial = (mxc(1)*LLog**6d0 + mxc(2)*LLog**5d0 + &
                    mxc(3)*LLog**4d0 + mxc(4)*LLog**3d0 + &
                    mxc(5)*LLog**2d0 + mxc(6)*LLog      + mxc(7))*w

    ! back to the user...
    return

  end function ospolynomial
!
!--------------------------------------------------------------------
!
end module CARBON_MODEL_MOD
