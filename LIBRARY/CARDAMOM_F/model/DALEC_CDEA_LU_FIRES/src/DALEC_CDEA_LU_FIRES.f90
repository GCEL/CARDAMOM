
module CARBON_MODEL_MOD

implicit none

! make all private
private

! explicit publics
public :: CARBON_MODEL     &
         ,soil_frac_clay   &
         ,soil_frac_sand   &
         ,nos_soil_layers  &
         ,extracted_C      &
         ,CiCa_time        &
         ,dim_1,dim_2      &
         ,nos_trees        &
         ,nos_inputs       &
         ,leftDaughter     &
         ,rightDaughter    &
         ,nodestatus       &
         ,xbestsplit       &
         ,nodepred         &
         ,bestvar

! Biomass removal (e.g. due to forest harvest)
double precision, allocatable, dimension(:) :: extracted_C, CiCa_time

! Variables needed incase of using random forest functions.
! None are currently implemented but variables remain for legacy reasons
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

! Multiple soil layer variables, these are not used in DALEC2 (C1),
! but declarations are needed to here ensure compilation compatability with more complex verison of DALEC
integer, parameter :: nos_root_layers = 2, nos_soil_layers = nos_root_layers + 1
double precision :: ci
double precision, dimension(nos_soil_layers) :: soil_frac_clay,soil_frac_sand

contains
!
!--------------------------------------------------------------------
!
  subroutine CARBON_MODEL(start,finish,met,pars,deltat,nodays,lat,lai,NEE,FLUXES,POOLS &
                       ,nopars,nomet,nopools,nofluxes,GPP)

    ! The Data Assimilation Linked Ecosystem Carbon - Combined Deciduous
    ! Evergreen Analytical (DALEC_CDEA; C1; DALEC2) model.
    ! The CARDBON_MODEL subroutine uses version 1 of the Aggregated Canopy Model (ACM)
    ! to simulate GPP. GPP is then partitied to autotrophic respiration and
    ! the live pools (foliage, labile, wood, fine roots).
    ! These live pools are subject to turnover to dead organic matter pools which are subsequently
    ! decomposed resulting in heterotrophic respiration

    ! This version includes the option to simulate fire combustion based
    ! on burned fraction and fixed combusion rates. It also includes the
    ! possibility to remove a fraction of biomass to simulate deforestation.

    implicit none

    ! Declare input dimensions
    integer, intent(in) :: start    & ! Start time step of the current call
                          ,finish   & ! End time stem of the current call
                          ,nopars   & ! number of paremeters in vector
                          ,nomet    & ! number of meteorological fields
                          ,nofluxes & ! number of model fluxes
                          ,nopools  & ! number of model pools
                          ,nodays     ! number of days in simulation
    ! Declare input drivers
    double precision, intent(in) :: met(nomet,nodays) & ! met drivers
                                   ,deltat(nodays)    & ! time step in decimal days
                                   ,pars(nopars)      & ! number of parameters
                                   ,lat                 ! site latitude (degrees)
    ! Declare model output variables
    double precision, dimension(nodays), intent(inout) :: lai & ! leaf area index
                                                         ,GPP & ! Gross primary productivity
                                                         ,NEE   ! net ecosystem exchange of CO2
    double precision, dimension((nodays+1),nopools), intent(inout) :: POOLS ! vector of ecosystem pools

    double precision, dimension(nodays,nofluxes), intent(inout) :: FLUXES ! vector of ecosystem fluxes

    ! declare local variables
    double precision :: gpppars(12)   & ! ACM inputs (LAI+met)
                       ,constants(10) & ! parameters for ACM
             ,wf,wl,ff,fl,osf,osl,sf  & ! phenological controls
             ,pi,ml,doy

    ! C pool specific combustion completeness and resilience factors
    double precision :: cf(6),rfac(6)
    integer :: p,f,n

    ! met drivers are:
    ! 1st run day
    ! 2nd min daily temp (oC)
    ! 3rd max daily temp (oC)
    ! 4th Radiation (MJ.m-2.day-1)
    ! 5th CO2 (ppm)
    ! 6th DOY at end of time step
    ! 7th Precipitation (kgH2O/m2/s - not used in this model version)
    ! 8th Fraction of biomass removed
    ! 9th Burned fraction

    ! POOLS (gC/m2) are:
    ! 1 = labile
    ! 2 = foliar
    ! 3 = root
    ! 4 = wood
    ! 5 = litter
    ! 6 = som

    ! FLUXES (gC/m2/day) are:
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
    ! 17 = ecosystem fire emission  (sum of fluxes 18 to 23)
    ! 18 = fire emission from labile
    ! 19 = fire emission from foliar
    ! 20 = fire emission from roots
    ! 21 = fire emission from wood
    ! 22 = fire emission from litter
    ! 23 = fire emission from soil
    ! 24 = transfer from labile into litter
    ! 25 = transfer from foliar into litter
    ! 26 = transfer from roots into litter
    ! 27 = transfer from wood into som
    ! 28 = transfer from litter into som

    ! PARAMETERS
    ! 17 values
    ! NOTE: C pool initial conditions are part of the parameter vector but are listed elsewhere

    ! p(1) Litter to SOM (day-1 at 0oC)
    ! p(2) Fraction of photosynthate respired as autotrophic
    ! p(3) Fraction of photosynthate allocated directly to foliage
    ! p(4) Fraction of photosynthate allocated to fine roots
    ! p(5) Leaf lifespan (year)
    ! p(6) Turnover rate of wood (day-1)
    ! p(7) Turnover rate of roots (day-1)
    ! p(8) Litter turnover rate (day-1 at 0oC)
    ! p(9) SOM turnover rate  (day-1 at 0oC)
    ! p(10) Exponential coefficient on temperature response
    ! p(11) Photosynthetic canopy efficiency parameter (gC/m2leaf/day at oC)
    ! p(12) = Julian day of peak labile release to canopy
    ! p(13) = Fraction of photosynthate allocated to labile
    ! p(14) = Labile release period (days)
    ! p(15) = Julian day of peak leaf fall
    ! p(16) = Leaf fall period (days)
    ! p(17) = Leaf Carbon per unit area (gC/m2leaf)

    ! set constants
    pi = 3.1415927d0

    ! load some values
    gpppars(4) = 1d0 ! foliar N (gN/m2leaf)
    gpppars(7) = lat
    gpppars(9) = -2d0 ! leafWP-soilWP (MPa)
    gpppars(10) = 1d0 ! total hydraulic resistance
    gpppars(11) = pi

    ! Assign acm parameters (see Fox et al., 2009)
    constants(1) = pars(11) ! Assign canopy efficency parameter
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
       ! Initial carbon pool (gC/m2) conditions
       POOLS(1,1) = pars(18) ! labile
       POOLS(1,2) = pars(19) ! foliar
       POOLS(1,3) = pars(20) ! roots
       POOLS(1,4) = pars(21) ! wood
       POOLS(1,5) = pars(22) ! litter
       POOLS(1,6) = pars(23) ! som
    endif

    ! Defining phenological variables
    ! The phenology model is based on leaf turnover and growth
    ! with a peak day and a range following a Gaussian distribution.

    ! Release period coefficient, based on duration of labile turnover or leaf
    ! fall durations
    wf = pars(16)*sqrt(2d0) * 0.5d0
    wl = pars(14)*sqrt(2d0) * 0.5d0
    ! Magnitude coefficient
    ff = (log(pars(5))-log(pars(5)-1d0)) * 0.5d0
    fl = (log(1.001d0)-log(0.001d0)) * 0.5d0
    ! Set minium labile life span to one year
    ml = 1.001d0
    ! Offset for labile and leaf turnovers
    osf = ospolynomial(pars(5),wf)
    osl = ospolynomial(ml,wl)

    ! scaling to biyearly sine curve
    sf = 365.25d0/pi

    ! Define fire constants
    cf(1) = 0.1d0         ! labile combustion efficiency
    cf(2) = 0.9d0         ! foliar combustion efficiency
    cf(3) = 0.1d0         ! roots combustion efficiency
    cf(4) = 0.1d0         ! wood combustion efficiency
    cf(5) = 0.7d0         ! litter combustion efficiency
    cf(6) = 0.01d0        ! som combustion efficency
    rfac = 0.5d0          ! resilience factor
    rfac(5) = 0.1d0 ; rfac(6) = 0d0

    if (.not.allocated(CiCa_time)) allocate(CiCa_time(nodays))

    !
    ! Begin looping through each time step
    !

    do n = start, finish

      ! calculate LAI value
      lai(n) = POOLS(n,2)/pars(17)

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
      FLUXES(n,1) = acm(gpppars,constants) ; CiCa_time(n) = ci
      ! Exponential temperature modified rate of metabolic activity
      FLUXES(n,2) = exp(pars(10)*0.5d0*(met(3,n)+met(2,n)))
      ! Autotrophic respiration (gC.m-2.day-1)
      FLUXES(n,3) = pars(2)*FLUXES(n,1)
      ! Leaf production rate (gC.m-2.day-1)
      FLUXES(n,4) = (FLUXES(n,1)-FLUXES(n,3))*pars(3)
      ! Labile production (gC.m-2.day-1)
      FLUXES(n,5) = (FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,4))*pars(13)
      ! Fine root production (gC.m-2.day-1)
      FLUXES(n,6) = (FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,4)-FLUXES(n,5))*pars(4)
      ! Wood production
      FLUXES(n,7) = FLUXES(n,1)-FLUXES(n,3)-FLUXES(n,4)-FLUXES(n,5)-FLUXES(n,6)

      ! Labile release and leaffall factors
      FLUXES(n,9) = (2d0/sqrt(pi))*(ff/wf)*exp(-(sin((doy-pars(15)+osf)/sf)*sf/wf)**2d0)
      FLUXES(n,16) = (2d0/sqrt(pi))*(fl/wl)*exp(-(sin((doy-pars(12)+osl)/sf)*sf/wl)**2d0)

      !
      ! C pool turnover as function of time only
      !

      ! Labile release
      FLUXES(n,8) = POOLS(n,1)*(1d0-(1d0-FLUXES(n,16))**deltat(n))/deltat(n)
      ! Leaf litter production
      FLUXES(n,10) = POOLS(n,2)*(1d0-(1d0-FLUXES(n,9))**deltat(n))/deltat(n)
      ! Wood litter production
      FLUXES(n,11) = POOLS(n,4)*(1d0-(1d0-pars(6))**deltat(n))/deltat(n)
      ! Fine root litter production
      FLUXES(n,12) = POOLS(n,3)*(1d0-(1d0-pars(7))**deltat(n))/deltat(n)

      !
      ! C pool turnover as function of temperature AND time
      !

      ! respiration heterotrophic litter
      FLUXES(n,13) = POOLS(n,5)*(1d0-(1d0-FLUXES(n,2)*pars(8))**deltat(n))/deltat(n)
      ! respiration heterotrophic som
      FLUXES(n,14) = POOLS(n,6)*(1d0-(1d0-FLUXES(n,2)*pars(9))**deltat(n))/deltat(n)
      ! litter to som
      FLUXES(n,15) = POOLS(n,5)*(1d0-(1d0-pars(1)*FLUXES(n,2))**deltat(n))/deltat(n)

      ! Calculate the NEE,
      ! NOTE: NEE does not include Fire but we assume that NBE does
      NEE(n) = (-FLUXES(n,1)+FLUXES(n,3)+FLUXES(n,13)+FLUXES(n,14))
      ! Load GPP to output variable
      GPP(n) = FLUXES(n,1)

      !
      ! Update pools based on natural processes
      !

      ! labile pool
      POOLS(n+1,1) = POOLS(n,1) + (FLUXES(n,5)-FLUXES(n,8))*deltat(n)
      ! foliar pool
      POOLS(n+1,2) = POOLS(n,2) + (FLUXES(n,4)-FLUXES(n,10) + FLUXES(n,8))*deltat(n)
      ! root pool
      POOLS(n+1,3) = POOLS(n,3) + (FLUXES(n,6) - FLUXES(n,12))*deltat(n)
      ! wood pool
      POOLS(n+1,4) = POOLS(n,4) + (FLUXES(n,7)-FLUXES(n,11))*deltat(n)
      ! litter pool
      POOLS(n+1,5) = POOLS(n,5) + (FLUXES(n,10)+FLUXES(n,12)-FLUXES(n,13)-FLUXES(n,15))*deltat(n)
      ! som pool
      POOLS(n+1,6) = POOLS(n,6) + (FLUXES(n,15)-FLUXES(n,14)+FLUXES(n,11))*deltat(n)

      !
      ! Update pools based on disturbance
      !

      ! Biomass removal
      if (met(8,n) > 0d0) then
          ! Only estimate total when code is conducting post-processing
          if (allocated(extracted_C)) then
              extracted_C(n) = ((POOLS(n+1,1) + POOLS(n+1,2) + POOLS(n+1,4)) * met(8,n)) / deltat(n)
          endif
          ! Update live pools
          POOLS(n+1,1) = POOLS(n+1,1)*(1d0-met(8,n)) ! remove labile
          POOLS(n+1,2) = POOLS(n+1,2)*(1d0-met(8,n)) ! remove foliar
          POOLS(n+1,4) = POOLS(n+1,4)*(1d0-met(8,n)) ! remove wood
          ! NOTE 1: fine root is left in system without mortality
          ! NOTE 2: no harvest residues are assumed
      end if

      ! Fire - based on burned fraction
      if (met(9,n) > 0d0) then

          ! First calculate combustion / emissions fluxes (gC/m2/day)
          FLUXES(n,18) = POOLS(n+1,1)*met(9,n)*cf(1)/deltat(n) ! labile
          FLUXES(n,19) = POOLS(n+1,2)*met(9,n)*cf(2)/deltat(n) ! foliar
          FLUXES(n,20) = POOLS(n+1,3)*met(9,n)*cf(3)/deltat(n) ! roots
          FLUXES(n,21) = POOLS(n+1,4)*met(9,n)*cf(4)/deltat(n) ! wood
          FLUXES(n,22) = POOLS(n+1,5)*met(9,n)*cf(5)/deltat(n) ! litter
          FLUXES(n,23) = POOLS(n+1,6)*met(9,n)*cf(6)/deltat(n) ! som

          ! Second calculate litter transfer fluxes (gC/m2/day), all pools except som
          FLUXES(n,24) = POOLS(n+1,1)*met(9,n)*(1d0-cf(1))*(1d0-rfac(1))/deltat(n) ! labile into litter
          FLUXES(n,25) = POOLS(n+1,2)*met(9,n)*(1d0-cf(2))*(1d0-rfac(2))/deltat(n) ! foliar into litter
          FLUXES(n,26) = POOLS(n+1,3)*met(9,n)*(1d0-cf(3))*(1d0-rfac(3))/deltat(n) ! roots into litter
          FLUXES(n,27) = POOLS(n+1,4)*met(9,n)*(1d0-cf(4))*(1d0-rfac(4))/deltat(n) ! wood into som
          FLUXES(n,28) = POOLS(n+1,5)*met(9,n)*(1d0-cf(5))*(1d0-rfac(5))/deltat(n) ! litter into som

          ! Update pools - remove burned vegetation...
          POOLS(n+1,1) = POOLS(n+1,1) - (FLUXES(n,18) + FLUXES(n,24)) * deltat(n) ! labile
          POOLS(n+1,2) = POOLS(n+1,2) - (FLUXES(n,19) + FLUXES(n,25)) * deltat(n) ! foliar
          POOLS(n+1,3) = POOLS(n+1,3) - (FLUXES(n,20) + FLUXES(n,26)) * deltat(n) ! roots
          POOLS(n+1,4) = POOLS(n+1,4) - (FLUXES(n,21) + FLUXES(n,27)) * deltat(n) ! wood
          ! Update pools - ...then redistribute litter
          POOLS(n+1,5) = POOLS(n+1,5) + (FLUXES(n,24) + FLUXES(n,25) + FLUXES(n,26) - FLUXES(n,22) - FLUXES(n,28)) * deltat(n)
          POOLS(n+1,6) = POOLS(n+1,6) + (FLUXES(n,27) + FLUXES(n,28) - FLUXES(n,23)) * deltat(n)

          ! Calculate ecosystem emissions (gC/m2/day)
          FLUXES(n,17) = FLUXES(n,18)+FLUXES(n,19)+FLUXES(n,20)+FLUXES(n,21)+FLUXES(n,22)+FLUXES(n,23)

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
    ! See Williams et al., (1997) and Fox et al., (2009) for details

    implicit none

    ! declare input variables
    double precision, intent(in) :: drivers(12) & ! acm input requirements
                                   ,constants(10) ! ACM parameters

    ! declare local variables
    double precision :: gc, pn, pd, pp, qq, e0, dayl, cps, dec, nit &
                       ,trange, sinld, cosld,aob,pi, mult &
                       ,mint,maxt,radiation,co2,lai,doy,lat &
                       ,deltaWP,Rtot,NUE,temp_exponent,dayl_coef &
                       ,dayl_const,hydraulic_exponent,hydraulic_temp_coef &
                       ,co2_comp_point,co2_half_sat,lai_coef,lai_const

    ! Initial values
    gc = 0d0 ; pp = 0d0 ; qq = 0d0 ; ci = 0d0 ; e0 = 0d0 ; dayl = 0d0 ; cps = 0d0 ; dec = 0d0

    ! load driver values to correct local vars
    lai = drivers(1)  ! leaf area (m2/m2)
    maxt = drivers(2) ! Daily maximum temperature (oC)
    mint = drivers(3) ! Daily minimum temperature (oC)
    nit = drivers(4)  ! Load foliar N (gN/m2leaf)
    co2 = drivers(5)  ! Atmospheric CO2 (ppm)
    doy = drivers(6)  ! Julian Day of year
    radiation = drivers(8) ! Short wave radation (MJ/m2/day)
    lat = drivers(7) ! latitude (degrees)

    ! Load parameters into correct local vars.
    ! While strictly not needed, it was done to improve readability of the below code
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

    ! Determine temperature range (oC)
    trange = 0.5d0*(maxt-mint)
    ! Daily canopy conductance
    gc = abs(deltaWP)**(hydraulic_exponent)/((hydraulic_temp_coef*Rtot+trange))
    ! Maximum rate of temperature and nitrogen (canopy efficiency) limited photosynthesis (gC.m-2.day-1)
    pn = lai*nit*NUE*exp(temp_exponent*maxt)
    ! pp and qq represent limitation by diffusion and metabolites respectively
    pp = pn/gc ; qq = co2_comp_point-co2_half_sat
    ! calculate internal CO2 concentration (ppm)
    ci = 0.5d0*(co2+qq-pp+((co2+qq-pp)**2d0-4d0*(co2*qq-pp*co2_comp_point))**0.5d0)
    ! limit maximum quantium efficiency by leaf area, hyperbola
    e0 = lai_coef*lai**2d0/(lai**2d0+lai_const)
    ! calculate day length (hours)
!    dec = - asin( sin( 23.45d0 * pi / 180d0 ) * cos( 2d0 * pi * ( doy + 10d0 ) /365d0 ) )
!    sinld = sin( lat*(pi/180d0) ) * sin( dec )
!    cosld = cos( lat*(pi/180d0) ) * cos( dec )
!    aob = max(-1d0,min(1d0,sinld / cosld))
!    dayl = 12d0 * ( 1d0 + 2d0 * asin( aob ) / pi )

!--------------------------------------------------------------
    ! Calculate day length (hours)
    ! This is the old REFLEX project calculation it is less efficient than the commented code above,
    ! however is currently returned for comparison with JPL DALEC version
    dec = -23.4d0*cos((360d0*(doy+10d0)/365d0)*pi/180d0)*pi/180d0
    mult = tan(lat*pi/180.0)*tan(dec)
    if (mult >= 1d0) then
        dayl = 2d0
    else if (mult <= -1d0) then
        dayl = 0d0
    else
        dayl = 24d0*acos(-mult)/pi
    end if
! ---------------------------------------------------------------
    ! Calculate CO2 limited rate of photosynthesis (gC/m2/day)
    pd = gc*(co2-ci)
    ! Calculate combined light and CO2 limited photosynthesis
    cps = e0*radiation*pd/(e0*radiation+pd)
    ! Correct for day length variation
    acm = cps*(dayl_coef*dayl+dayl_const)

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
