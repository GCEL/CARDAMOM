
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
public :: CARBON_MODEL &
         ,soil_frac_clay   &
         ,soil_frac_sand   &
         ,nos_soil_layers  &

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
                          ,nopars   & ! number of paremeters in vector
                          ,pft      & ! plant functional type
                          ,nomet    & ! number of meteorological fields
                          ,nofluxes & ! number of model fluxes
                          ,nopools  & ! number of model pools
                          ,nodays     ! number of days in simulation

    double precision, intent(in) :: met(nomet,nodays) & ! met drivers
                         ,deltat            & ! time step in decimal days
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
             ,pi,ml

    integer :: p,f,nxp,n

    ! met drivers are:
    ! 1st run day
    ! 2nd min daily temp (oC)
    ! 3rd max daily temp (oC)
    ! 4th Radiation (MJ.m-2.day-1)
    ! 5th CO2 (ppm)
    ! 6th DOY

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
    ! 8 = labile production
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
    gpppars(4) = 1. ! foliar N
    gpppars(7) = lat
    gpppars(11) = pi

    ! assign acm parameters
    if (pft == 1) then
       ! cropland
       constants(1)=6.95981195652514    ; constants(2)=0.00453328640065
       constants(3)=21.1845749054962    ; constants(4)=206.9663487061
       constants(5)=0.0195949012471     ; constants(6)=28.1563485082444
       constants(7)=21.2279935364026    ; constants(8)=0.25188913401713
       constants(9)=3.82828052791062    ; constants(10)=6.82498077785962
       gpppars(9) = -1.73354511234508   ; gpppars(10) = 3.174831174986
    else if (pft == 2) then
       ! short grass
       constants(1)=31.0871524859723    ; constants(2)=0.00114118708815
       constants(3)=31.4467939449551    ; constants(4)=378.427398145854
       constants(5)=0.00890521506038    ; constants(6)=16.7637870155259
       constants(7)=59.547386641698     ; constants(8)=0.23818779495329
       constants(9)=3.02763270690527    ; constants(10)=4.87935779546392
       gpppars(9) = -2.52044119689308   ; gpppars(10) = 3.0232119402703
    else if (pft == 3) then
       ! evergreen needle leaf
       constants(1)=6.61869686621632    ; constants(2)=0.00489208053559
       constants(3)=55.9906543138429    ; constants(4)=331.497235959547
       constants(5)=7.14687281930546E-6 ; constants(6)=20.3339946560466
       constants(7)=20.9644279372219    ; constants(8)=0.19417808060063
       constants(9)=4.00491698189447E-9 ; constants(10)=3.43186213468807
       gpppars(9)= -2.45774825025118    ; gpppars(10)=2.95226577057285
    else if (pft == 5) then
       ! Deciduous broadleaf
       constants(1)=7.2802184902159     ; constants(2)=0.00429221113218
       constants(3)=1.71237894565241    ; constants(4)=313.431887180483
       constants(5)=0.03215193522376    ; constants(6)=5.30200375501687
       constants(7)=56.5115393175829    ; constants(8)=0.19389077826375
       constants(9)=26.3520696913012    ; constants(10)=1.86281362380134
       gpppars(9) = -2.22723093121689   ; gpppars(10) = 1.39808690504974
    else if (pft == 18) then
       ! Mixed forest
       constants(1)=5.50924877737081    ; constants(2)=0.00641412147596
       constants(3)=1.615462505343      ; constants(4)=370.753230047059
       constants(5)=7.00874375394324E-7 ; constants(6)=24.7888414195468
       constants(7)=11.4551733716705    ; constants(8)=0.20409368420038
       constants(9)=0.00044072300905    ; constants(10)=2.24405683848548
       gpppars(9) = -3.56882236413978   ; gpppars(10) = 2.15773459872329
    else
       print*,"FUCK FUCK FUCK - PFT is missing or misidentified pft input = ", pft
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
    wf=pars(15)*sqrt(2.)/2.
    wl=pars(13)*sqrt(2.)/2.

    ! magnitude coefficient
    ff=(log(pars(5))-log(pars(5)-1.))/2.
    fl=(log(1.001)-log(0.001))/2.

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
      gpppars(6)=met(6,n) ! doy
      gpppars(8)=met(4,n) ! radiation

      ! GPP (gC.m-2.day-1)
      FLUXES(n,1) = acm(gpppars,constants)
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
      FLUXES(n,9) = (2./(pi**0.5))*(ff/wf)*exp(-((sin((met(1,n)-pars(14)+osf)/sf)*sf/wf)**2.))
      FLUXES(n,16) = (2./(pi**0.5))*(fl/wl)*exp(-((sin((met(1,n)-pars(11)+osl)/sf)*sf/wl)**2.))

      !
      ! those with time dependancies
      !

      ! total labile release
      FLUXES(n,8) = POOLS(n,1)*(1.-(1.-FLUXES(n,16))**deltat)/deltat
      ! total leaf litter production
      FLUXES(n,10) = POOLS(n,2)*(1.-(1.-FLUXES(n,9))**deltat)/deltat
      ! total wood production
      FLUXES(n,11) = POOLS(n,4)*(1.-(1.-pars(6))**deltat)/deltat
      ! total root litter production
      FLUXES(n,12) = POOLS(n,3)*(1.-(1.-pars(7))**deltat)/deltat

      !
      ! those with temperature AND time dependancies
      !

      ! respiration heterotrophic litter
      FLUXES(n,13) = POOLS(n,5)*(1.-(1.-FLUXES(n,2)*pars(8))**deltat)/deltat
      ! respiration heterotrophic som
      FLUXES(n,14) = POOLS(n,6)*(1.-(1.-FLUXES(n,2)*pars(9))**deltat)/deltat
      ! litter to som
      FLUXES(n,15) = POOLS(n,5)*(1.-(1.-pars(1)*FLUXES(n,2))**deltat)/deltat

      ! calculate the NEE
      NEE(n) = (-FLUXES(n,1)+FLUXES(n,3)+FLUXES(n,13)+FLUXES(n,14))
      ! load GPP
      GPP(n) = FLUXES(n,1)

      !
      ! update pools for next timestep
      !

      ! labile pool
      POOLS(n+1,1) = POOLS(n,1) + (FLUXES(n,5)-FLUXES(n,8))*deltat
      ! foliar pool
      POOLS(n+1,2) = POOLS(n,2) + (FLUXES(n,4)-FLUXES(n,10) + FLUXES(n,8))*deltat
      ! wood pool
      POOLS(n+1,4) = POOLS(n,4) + (FLUXES(n,7)-FLUXES(n,11))*deltat
      ! root pool
      POOLS(n+1,3) = POOLS(n,3) + (FLUXES(n,6) - FLUXES(n,12))*deltat
      ! litter pool
      POOLS(n+1,5) = POOLS(n,5) + (FLUXES(n,10)+FLUXES(n,12)-FLUXES(n,13)-FLUXES(n,15))*deltat
      ! som pool
      POOLS(n+1,6) = POOLS(n,6) + (FLUXES(n,15)-FLUXES(n,14)+FLUXES(n,11))*deltat

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
    gc=0 ; pp=0 ; qq=0 ; ci=0 ; e0=0 ; dayl=0 ; cps=0 ; dec=0 ; nit=1.

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
    trange=0.5*(maxt-mint)
    ! daily canopy conductance
    gc=abs(deltaWP)**(hydraulic_exponent)/((hydraulic_temp_coef*Rtot+trange))
    ! maximum rate of temperature and nitrogen (canopy efficiency) limited photosynthesis (gC.m-2.day-1)
    pn=lai*nit*NUE*exp(temp_exponent*maxt)
    ! pp and qq represent limitation by diffusion and metabolites respecitively
    pp=pn/gc ; qq=co2_comp_point-co2_half_sat
    ! calculate internal CO2 concentration (ppm)
    ci=0.5*(co2+qq-pp+((co2+qq-pp)**2.-4.*(co2*qq-pp*co2_comp_point))**0.5)
    ! limit maximum quantium efficiency by leaf area, hyperbola
    e0=lai_coef*lai**2./(lai**2.+lai_const)
    ! calculate day length (hours)
    dec = - asin( sin( 23.45 * pi / 180. ) * cos( 2. * pi * ( doy + 10. ) /365. ) )
    sinld = sin( lat*(pi/180.) ) * sin( dec )
    cosld = cos( lat*(pi/180.) ) * cos( dec )
    aob = max(-1.0,min(1.0,sinld / cosld))
    dayl = 12.0 * ( 1. + 2. * asin( aob ) / pi )
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
    LLog=log(L-1.)

    ! calculate the polynomial function
    ospolynomial=(mxc(1)*LLog**6. + mxc(2)*LLog**5. + &
                  mxc(3)*LLog**4. + mxc(4)*LLog**3. + &
                  mxc(5)*LLog**2. + mxc(6)*LLog     + mxc(7))*w

  end function ospolynomial
!
!--------------------------------------------------------------------
!
end module CARBON_MODEL_MOD
