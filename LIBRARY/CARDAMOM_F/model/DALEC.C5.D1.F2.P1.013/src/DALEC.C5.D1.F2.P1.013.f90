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
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

!!!!!!!!!!!! File specific description !!!!!!!!!!
! This file contains the source code of DALEC.C5.D1.F2.P1.013
!
! This code contains a variant of the Data Assimilation Linked ECosystem (DALEC) model.
! This version of DALEC is derived from the following primary references:
! Bloom & Williams (2015), https://doi.org/10.5194/bg-12-1299-2015.
! This code is based on that created by A. A. Bloom (UoE, now at JPL, USA).
! Subsequent modifications by:
! T. L. Smallman (University of Edinburgh, t.l.smallman@ed.ac.uk)
! See function / subroutine specific comments for exceptions and contributors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module CARBON_MODEL_MOD

implicit none

! make all private
private

  ! explicit publics
  public :: CARBON_MODEL  &
           ,initialize_mv

  use carbon_model_memory

  contains
  !
  !--------------------------------------------------------------------
  !
  subroutine initialize_mv(mV, nodays, nomet, nopars, deltat, soil_frac_sand, soil_frac_clay, met, lat)

    !! For a single chain's model_working_varibles type object mV, allocate arrays
    !! and calculate initial values.
    implicit none

    type(model_working_variables), intent(out):: mV
    integer, intent(in):: nodays, nomet, nopars
    double precision, intent(in) :: deltat(nodays)     ! time step in decimal days
    double precision, intent(in), dimension(:) :: soil_frac_sand, soil_frac_clay
    double precision, intent(in) :: met(nomet, nodays)  ! met drivers
    double precision, intent(in) :: lat

    integer:: n

    mV%soil_frac_sand = soil_frac_sand
    mV%soil_frac_clay = soil_frac_clay
       
  end subroutine initialize_mv
  !
  !--------------------------------------------------------------------
  !
  subroutine destroy_mv(mV)
    !! deallocate arrays in mV
    type(model_working_variables):: mV
    ! nothing to do here, but we keep the method because it's 
    ! called from general parts of the code for all models
  end subroutine
  !
  !--------------------------------------------------------------------
  !
  subroutine CARBON_MODEL(start,finish,met,pars,deltat,nodays,lat,FLUXES,POOLS,DIAGS &
                         ,nopars,nomet,nopools,nofluxes,nodiags, mV)

    ! The Data Assimilation Linked Ecosystem Carbon - Combined Deciduous
    ! Evergreen Analytical (DALEC.C5.D1.F2.P1.013) model. The subroutine calls the
    ! Aggregated Canopy Model to simulate GPP and partitions between various
    ! ecosystem carbon pools. These pools are subject to turnovers /
    ! decompostion resulting in ecosystem phenology and fluxes of CO2

    ! This version includes the option to simulate fire combustion based
    ! on burned fraction and fixed combustion rates. It also includes the
    ! possibility to remove a fraction of biomass to simulate deforestation.

    ! declare input variables
    type(model_working_variables) :: mV
    integer, intent(in) :: start    &
                          ,finish   &
                          ,nopars   & ! number of parameters in vector
                          ,nomet    & ! number of meteorological fields
                          ,nofluxes & ! number of model fluxes
                          ,nopools  & ! number of model pools
                          ,nodays   & ! number of days in simulation
                          ,nodiags    ! number of model diagnositic variables

    double precision, intent(in) :: met(nomet,nodays) & ! met drivers
                         ,deltat(nodays)    & ! time step in decimal days
                         ,pars(nopars)      & ! number of parameters
                         ,lat                 ! site latitude (degrees)

    double precision, dimension((nodays+1),nopools), intent(inout) :: POOLS ! vector of ecosystem pools
    double precision, dimension(nodays,nofluxes), intent(inout) :: FLUXES ! vector of ecosystem fluxes
    double precision, dimension(nodays,nodiags), intent(inout) :: DIAGS ! vector of ecosystem diagnostics

    ! declare local variables
    double precision :: infi          &
                       ,gpppars(10)   & ! ACM inputs (LAI+met)
                       ,constants(10) & ! parameters for ACM
                       ,wf,wl,ff,fl   & ! phenological controls
                       ,osf,osl,sf,ml & 
                       ,doy,tmp

    ! C pool specific combustion completeness and resilience factors
    double precision :: burnt_area
    double precision, dimension(4) :: cf,rfac
    ! local deforestation related variables
    double precision, dimension(5) :: post_harvest_burn      & ! how much burning to occur after
                                     ,foliage_frac_res       &
                                     ,roots_frac_res         &
                                     ,rootcr_frac_res        &
                                     ,stem_frac_res          &
                                     ,roots_frac_removal     &
                                     ,rootcr_frac_removal    &
                                     ,Crootcr_part           &
                                     ,soil_loss_frac
    double precision :: labile_loss,foliar_loss      &
                       ,roots_loss,wood_loss         &
                       ,rootcr_loss,stem_loss        &
                       ,labile_residue,foliar_residue&
                       ,roots_residue,wood_residue   &
                       ,C_total,labile_frac_res      &
                       ,labile_frac_removal          &
                       ,Cstem,Crootcr,stem_residue   &
                       ,coarse_root_residue          &
                       ,soil_loss_with_roots

    integer :: p, f, n, harvest_management                       

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
    ! 1 = labile       (gC/m2) (initial value: p14)
    ! 2 = foliar       (gC/m2) (initial value: p15)
    ! 3 = root + wood  (gC/m2) (initial value: p16)
    ! 4 = litter + som (gC/m2) (initial value: p17)

    ! FLUXES are:
    ! 1  = GPP (gC/m2/day)
    ! 2  = temperature rate modifier (unitless)
    ! 3  = autotrophic respiration (gC/m2/day)
    ! 4  = GPP allocation to foliage (gC/m2/day)
    ! 5  = GPP allocation to labile (gC/m2/day)
    ! 6  = GPP allocation to wood+root (gC/m2/day)
    ! 7  = NOT IN USE
    ! 8  = labile->leaf transfer (gC/m2/day)
    ! 9  = leaf fall factor (fraction/day)
    ! 10 = leaf litter production (gC/m2/day)
    ! 11 = wood+root litter production (gC/m2/day)
    ! 12 = NOT IN USE
    ! 13 = heterotrophic respiration from litter+som (gC/m2/day)
    ! 14 = NOT IN USE
    ! 15 = NOT IN USE
    ! 16 = labile release factor (fraction/day)
    ! 17 = total ecosystem fire emission - sum(18:21) (gC/m2/day)
    ! 18 = fire emission from labile (gC/m2/day)
    ! 19 = fire emission from foliage (gC/m2/day)
    ! 20 = fire emission from roots+wood (gC/m2/day)
    ! 21 = fire emission from litter+som (gC/m2/day)
    ! 22 = fire mortality transfer labile->litter+som (gC/m2/day)
    ! 23 = fire mortality transfer foliage->litter+som (gC/m2/day)
    ! 24 = fire mortality transfer roots+wood->litter+som (gC/m2/day)
    ! 25 = total harvest extracted C - sum(26:29) (gC/m2/day)
    ! 26 = harvest extraction from labile (gC/m2/day)
    ! 27 = harvest extraction from foliage (gC/m2/day)
    ! 28 = harvest extraction from wood+root (gC/m2/day)
    ! 29 = harvest extraction from litter+som (gC/m2/day)
    ! 30 = harvest litter residue from labile (gC/m2/day)
    ! 31 = harvest litter residue from foliage (gC/m2/day)
    ! 32 = harvest litter residue from wood+root (gC/m2/day)

    ! PARAMETERS are:
    ! p(1)  = fraction of GPP as autotrophic respiration (fraction)
    ! p(2)  = fraction of NPP allocated to foliage (fraction)
    ! p(3)  = leaf lifespan (yr)
    ! p(4)  = wood+root turnover rate (fraction/day)
    ! p(5)  = litter+som turnover rate, temperature adjusted (fraction/day)
    ! p(6)  = temperature sensitivity of heterotrophic respiration (oC-1)
    ! p(7)  = canopy efficiency (gC/m2leaf/day)
    ! p(8)  = date of labile release / bud burst (day of year)
    ! p(9)  = fraction of NPP allocated to labile pool (fraction)
    ! p(10) = labile release period duration (days)
    ! p(11) = date of leaf fall (day of year)
    ! p(12) = leaf fall period duration (days)
    ! p(13) = leaf mass per area LMA (gC/m2)
    ! p(14) = initial labile C pool (gC/m2)
    ! p(15) = initial foliar C pool (gC/m2)
    ! p(16) = initial wood+root C pool (gC/m2)
    ! p(17) = initial litter+som C pool (gC/m2)
    ! p(18) = fire resilience factor for non-combusted C (fraction)
    ! p(19) = combustion completeness for foliage (fraction)
    ! p(20) = combustion completeness for roots+wood (fraction)
    ! p(21) = combustion completeness for litter+som (fraction)

    ! Set some initial states
    infi = 0d0 ; FLUXES = 0d0 ; POOLS = 0d0 ; DIAGS = 0d0

    ! load acm inputs 
    gpppars(4) = 1d0 ! foliar N
    gpppars(7) = lat ! latitude in degrees
    gpppars(9) = -2d0 ! leafWP-soilWP
    gpppars(10) = 1d0 ! totaly hydraulic resistance

    ! assign acm parameters (see Fox et al., 2009)
    constants(1) = pars(11) ! canopy efficiency
    constants(2) = 0.0156935d0
    constants(3) = 4.22273d0
    constants(4) = 208.868d0
    constants(5) = 0.0453194d0
    constants(6) = 0.37836d0
    constants(7) = 7.19298d0
    constants(8) = 0.011136d0
    constants(9) = 2.1001d0
    constants(10) = 0.789798d0    

    ! assigning initial conditions
    POOLS(1,1) = pars(14) ! labile
    POOLS(1,2) = pars(15) ! foliar
    POOLS(1,3) = pars(16) ! roots+wood
    POOLS(1,4) = pars(17) ! litter + som

    ! Release period coefficient, based on duration of labile turnover or leaf
    ! fall durations
    wf = pars(12)*sqrt(2d0) * 0.5d0
    wl = pars(10)*sqrt(2d0) * 0.5d0
    ! Magnitude coefficient
    ff = (log(pars(3))-log(pars(3)-1d0)) * 0.5d0
    fl = 3.45437738965761021d0!(log(1.001d0)-log(0.001d0)) * 0.5d0
    ! Set minium labile life span to one year
    ml = 1.001d0
    ! Offset for labile and leaf turnovers
    osf = ospolynomial(pars(3),wf)
    osl = ospolynomial(ml,wl)

    ! scaling to biyearly sine curve
    sf = 116.262685928629551d0 !365.25d0/pi

    ! JFE added 4 May 2018 - define fire constants
    ! Update fire parameters derived from
    ! Yin et al., (2020), doi: 10.1038/s414647-020-15852-2
    ! Subsequently expanded by T. L. Smallman & Mat Williams (UoE, 03/09/2021)
    ! to provide specific CC for litter and wood litter.
    ! NOTE: changes also result in the addition of further EDCs

    ! now load the hardcoded forest management parameters into their scenario locations

    ! Deforestation process functions in a sequential way.
    ! Thus, the pool_loss is first determined as a function of met(8,n) and
    ! for fine and coarse roots whether this felling is associated with a mechanical
    ! removal from the ground. As the canopy and stem is removed (along with a proportion of labile)
    ! fine and coarse roots may subsequently undergo mortality from which they do not recover
    ! but allows for management activities such as grazing, mowing and coppice.
    ! The pool_loss is then partitioned between the material which is left within the system
    ! as a residue and thus directly placed within one of the dead organic matter pools.

    !! Parameter values for deforestation variables
    !! Scenario 1
    ! Define 'removal' for coarse and fine roots, i.e. fraction of imposed
    ! removal which is imposed directly on these pools. These fractions vary
    ! the assumption that the fine and coarse roots are mechanically removed.
    ! 1 = all removed, 0 = all remains.
    rootcr_frac_removal(1) = 0d0
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(1) = 1d0
    rootcr_frac_res(1)  = 1d0
    stem_frac_res(1)    = 0.20d0 !
    ! wood partitioning (fraction)
    Crootcr_part(1) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
    ! Csom loss due to phyical removal with roots
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(1) = 0.02d0 ! actually between 1-3 %
    ! was the forest burned after deforestation (0-1)
    ! NOTE: that we refer here to the fraction of the cleared land to be burned
    post_harvest_burn(1) = 1d0

    !! Scenario 2
    ! Define 'removal' for coarse and fine roots, i.e. fraction of imposed
    ! removal which is imposed directly on these pools. These fractions vary
    ! the assumption that the fine and coarse roots are mechanically removed.
    ! 1 = all removed, 0 = all remains.
    rootcr_frac_removal(2) = 0d0
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(2) = 1d0
    rootcr_frac_res(2)  = 1d0
    stem_frac_res(2)    = 0.20d0 !
    ! wood partitioning (fraction)
    Crootcr_part(2) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
    ! Csom loss due to phyical removal with roots
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(2) = 0.02d0 ! actually between 1-3 %
    ! was the forest burned after deforestation (0-1)
    ! NOTE: that we refer here to the fraction of the cleared land to be burned
    post_harvest_burn(2) = 0d0

    !! Scenario 3
    ! Define 'removal' for coarse and fine roots, i.e. fraction of imposed
    ! removal which is imposed directly on these pools. These fractions vary
    ! the assumption that the fine and coarse roots are mechanically removed.
    ! 1 = all removed, 0 = all remains.
    rootcr_frac_removal(3) = 0d0
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(3) = 0.5d0
    rootcr_frac_res(3)  = 1d0
    stem_frac_res(3)    = 0d0 !
    ! wood partitioning (fraction)
    Crootcr_part(3) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
    ! Csom loss due to phyical removal with roots
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(3) = 0.02d0 ! actually between 1-3 %
    ! was the forest burned after deforestation (0-1)
    ! NOTE: that we refer here to the fraction of the cleared land to be burned
    post_harvest_burn(3) = 0d0

    !! Scenario 4
    ! Define 'removal' for coarse and fine roots, i.e. fraction of imposed
    ! removal which is imposed directly on these pools. These fractions vary
    ! the assumption that the fine and coarse roots are mechanically removed.
    ! 1 = all removed, 0 = all remains.
    rootcr_frac_removal(4) = 1d0
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(4) = 0.5d0
    rootcr_frac_res(4)  = 0d0
    stem_frac_res(4)    = 0d0
    ! wood partitioning (fraction)
    Crootcr_part(4) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
    ! Csom loss due to phyical removal with roots
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(4) = 0.02d0 ! actually between 1-3 %
    ! was the forest burned after deforestation (0-1)
    ! NOTE: that we refer here to the fraction of the cleared land to be burned
    post_harvest_burn(4) = 0d0

    !## Scenario 5 (grassland grazing / cutting)
    ! Define 'removal' for coarse and fine roots, i.e. fraction of imposed
    ! removal which is imposed directly on these pools. These fractions vary
    ! the assumption that the fine and coarse roots are mechanically removed.
    ! 1 = all removed, 0 = all remains.
    rootcr_frac_removal(5) = 0d0
    ! harvest residue (fraction); 1 = all remains, 0 = all removed
    foliage_frac_res(5) = 0.1d0
    rootcr_frac_res(5)  = 0d0
    stem_frac_res(5)    = 0.12d0
    ! wood partitioning (fraction)
    Crootcr_part(5) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
    ! Csom loss due to phyical removal with roots
    ! Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac(5) = 0d0 ! actually between 1-3 %
    ! was the forest burned after deforestation (0-1)
    ! NOTE: that we refer here to the fraction of the cleared land to be burned
    post_harvest_burn(5) = 0d0

    ! Assign proposed resilience factor
    rfac(1:3) = pars(18)
    rfac(4) = 0d0
    ! Assign combustion completeness to foliage
    cf(2) = pars(19) ! foliage
    ! Assign combustion completeness to non-photosynthetic
    cf(1) = pars(20) ; cf(3) = pars(20)
    cf(4) = pars(21) ! dom

    !
    ! Begin looping through each time step
    !

    do n = start, finish

      ! calculate LAI value
      DIAGS(n,1) = POOLS(n,2)/pars(13)

      ! estimate multiple use variable
      doy = met(6,n)-(deltat(n)*0.5d0) ! doy

      ! load next met / lai values for ACM
      gpppars(1) = DIAGS(n,1)
      gpppars(2) = met(3,n) ! max temp
      gpppars(3) = met(2,n) ! min temp
      gpppars(5) = met(5,n) ! co2
      gpppars(6) = doy
      gpppars(8) = met(4,n) ! radiation

      ! GPP (gC.m-2.day-1)
      FLUXES(n,1) = acm(gpppars,constants, mV%ci) ; DIAGS(n,2) = mV%ci / met(5,n)
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
      FLUXES(n,9)  = (2d0/sqrt(pi))*(ff/wf)*exp(-(sin((doy-pars(11)+osf)/sf)*sf/wf)**2d0)
      FLUXES(n,16) = (2d0/sqrt(pi))*(fl/wl)*exp(-(sin((doy-pars(8)+osl)/sf)*sf/wl)**2d0)

      !
      ! those with time dependancies
      !

      ! total labile release
      FLUXES(n,8)  = POOLS(n,1)*(1d0-(1d0-FLUXES(n,16))**deltat(n))/deltat(n)
      ! total leaf litter production
      FLUXES(n,10) = POOLS(n,2)*(1d0-(1d0-FLUXES(n,9))**deltat(n))/deltat(n)
      ! total root+wood production
      FLUXES(n,11) = POOLS(n,3)*(1d0-(1d0-pars(4))**deltat(n))/deltat(n)

      !
      ! those with temperature AND time dependancies
      !

      ! respiration heterotrophic litter
      FLUXES(n,13) = POOLS(n,4)*(1d0-(1d0-FLUXES(n,2)*pars(5))**deltat(n))/deltat(n)

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

      !!!!!!!!!!
      ! Extract biomass - e.g. deforestation / degradation
      !!!!!!!!!!

      ! reset values
      harvest_management = 0 ; burnt_area = 0d0

      ! Does harvest activities occur?
      if (met(8,n) > 0d0) then

          ! Load the management type / scenario into local variable
          harvest_management = int(met(13,n))

          ! Determine the fraction of cut labile C which remains in system as residue.
          ! We assume that labile is proportionally distributed through the plants
          ! root and wood (structural C).
          C_total = POOLS(n+1,3)
          ! Ensure there is available C for extraction
          if (C_total > 0d0) then
              ! Harvest activities on the wood / structural pool varies depending on
              ! whether it is above or below ground. As such, partition the wood pool
              ! between above ground stem(+branches) and below ground coarse root.
              Crootcr = POOLS(n+1,3)*Crootcr_part(harvest_management)
              Cstem   = POOLS(n+1,3)-Crootcr
              ! Calculate the fraction of harvested labile which remains in system as residue
              labile_frac_res = ((Cstem/C_total)        * stem_frac_res(harvest_management)   ) &
                              + ((Crootcr/C_total)      * rootcr_frac_res(harvest_management) )
              ! Calculate the management scenario specific resistance fraction
              labile_frac_removal = ((Cstem/C_total)        * 1d0   ) &
                                  + ((Crootcr/C_total)      * rootcr_frac_removal(harvest_management) )

              ! Calculate the total loss from biomass pools
              ! We assume that fractional clearing always equals the fraction
              ! of foliage and above ground (stem) wood removal. However, we assume
              ! that coarse root and fine root extractions are dependent on the
              ! management activity type, e.g. in coppice below ground remains.
              ! Thus, labile extractions are also dependent.
              labile_loss = POOLS(n+1,1) * labile_frac_removal * met(8,n)
              foliar_loss = POOLS(n+1,2) * met(8,n)
              stem_loss   = (Cstem * met(8,n))
              rootcr_loss = (Crootcr * rootcr_frac_removal(harvest_management) * met(8,n))
              wood_loss   =  stem_loss + rootcr_loss

              ! Transfer fraction of harvest waste to litter, wood litter or som pools.
              ! This includes explicit calculation of the stem and coarse root residues due
              ! to their potentially different treatments under management scenarios
              labile_residue = labile_loss*labile_frac_res
              foliar_residue = foliar_loss*foliage_frac_res(harvest_management)
              coarse_root_residue = rootcr_loss*rootcr_frac_res(harvest_management)
              stem_residue = stem_loss*stem_frac_res(harvest_management)
              wood_residue = stem_residue + coarse_root_residue
              ! Mechanical loss of Csom due to coarse root extraction,
              ! less the loss remaining as residue
              soil_loss_with_roots = (rootcr_loss-coarse_root_residue) &
                                   * soil_loss_frac(harvest_management)

              ! Update pools
              POOLS(n+1,1) = POOLS(n+1,1) - labile_loss
              POOLS(n+1,2) = POOLS(n+1,2) - foliar_loss
              POOLS(n+1,3) = POOLS(n+1,3) - wood_loss
              POOLS(n+1,4) = POOLS(n+1,4) + (labile_residue+foliar_residue+wood_residue - soil_loss_with_roots)
              ! mass balance check
              where (POOLS(n+1,1:4) < 0d0) POOLS(n+1,1:4) = 0d0

              ! Convert harvest related extractions to daily rate for output
              ! For dead organic matter pools, in most cases these will be zeros.
              ! But these variables allow for subseqent management where surface litter
              ! pools are removed or mechanical extraction from soil occurs.
              FLUXES(n,26) = (labile_loss-labile_residue) / deltat(n)  ! Labile extraction
              FLUXES(n,27) = (foliar_loss-foliar_residue) / deltat(n)  ! foliage extraction
              FLUXES(n,28) = (wood_loss-wood_residue) / deltat(n)      ! wood extraction
              FLUXES(n,29) = soil_loss_with_roots / deltat(n)          ! som extraction
              ! Convert harvest related residue generations to daily rate for output
              FLUXES(n,30) = labile_residue / deltat(n) ! labile residues
              FLUXES(n,31) = foliar_residue / deltat(n) ! foliage residues
              FLUXES(n,32) = wood_residue / deltat(n)   ! wood residues

              ! Total C extraction, including any potential litter and som.
              FLUXES(n,25) = sum(FLUXES(n,26:29))

          end if ! C_total > 0d0

      endif ! end deforestation info

      !!!!!!!!!!
      ! Impose fire
      !!!!!!!!!!

      ! Fire - based on burned fraction
      if (met(9,n) > 0d0 .or.(met(8,n) > 0d0 .and. harvest_management > 0)) then

          ! Adjust burnt area to account for the managment decisions which may not be
          ! reflected in the burnt area drivers
          burnt_area = met(9,n)
          if (met(8,n) > 0d0 .and. burnt_area > 0d0) then
              ! pass harvest management to local integer
              burnt_area = min(1d0,burnt_area + post_harvest_burn(harvest_management))
          else if (met(8,n) > 0d0 .and. burnt_area <= 0d0) then
              burnt_area = post_harvest_burn(harvest_management)
          endif

          ! Determine the corrected burnt area
          if (burnt_area > 0d0) then

              ! first calculate combustion / emissions fluxes in g C m-2 d-1
              FLUXES(n,18) = POOLS(n+1,1)*burnt_area*cf(1)/deltat(n) ! labile
              FLUXES(n,19) = POOLS(n+1,2)*burnt_area*cf(2)/deltat(n) ! foliar
              FLUXES(n,20) = POOLS(n+1,3)*burnt_area*cf(3)/deltat(n) ! roots + wood
              FLUXES(n,21) = POOLS(n+1,4)*burnt_area*cf(4)/deltat(n) ! dom

              ! second calculate litter transfer fluxes in g C m-2 d-1, all pools except som
              FLUXES(n,22) = POOLS(n+1,1)*burnt_area*(1d0-cf(1))*(1d0-rfac(1))/deltat(n) ! labile into dom
              FLUXES(n,23) = POOLS(n+1,2)*burnt_area*(1d0-cf(2))*(1d0-rfac(2))/deltat(n) ! foliar into dom
              FLUXES(n,24) = POOLS(n+1,3)*burnt_area*(1d0-cf(3))*(1d0-rfac(3))/deltat(n) ! wood+root into dom

              ! update pools - first remove burned vegetation
              POOLS(n+1,1) = POOLS(n+1,1) - (FLUXES(n,18) + FLUXES(n,22)) * deltat(n) ! labile
              POOLS(n+1,2) = POOLS(n+1,2) - (FLUXES(n,19) + FLUXES(n,23)) * deltat(n) ! foliar
              POOLS(n+1,3) = POOLS(n+1,3) - (FLUXES(n,20) + FLUXES(n,24)) * deltat(n) ! roots+wood
              ! update pools - add litter transfer
              POOLS(n+1,4) = POOLS(n+1,4) + (FLUXES(n,22) + FLUXES(n,23) + FLUXES(n,24) - FLUXES(n,21)) * deltat(n)

              ! calculate ecosystem fire emissions (gC/m2/day)
              FLUXES(n,17) = FLUXES(n,18)+FLUXES(n,19)+FLUXES(n,20)+FLUXES(n,21)

          end if ! Burned_area > 0

      end if ! any fire?

    end do ! nodays loop

  end subroutine CARBON_MODEL
  !
  !------------------------------------------------------------------
  !
  double precision function acm(drivers,constants, ci)

    ! the Aggregated Canopy Model, is a Gross Primary Productivity (i.e.
    ! Photosyntheis) emulator which operates at a daily time step. ACM can be
    ! paramaterised to provide reasonable results for most ecosystems.

    implicit none

    ! declare input variables
    double precision, intent(in) :: drivers(10) & ! acm input requirements
                         ,constants(10) ! ACM parameters
    double precision, intent(inout) :: ci 

    ! declare local variables
    double precision :: gc, pn, pd, pp, qq, e0, dayl, cps, dec, nit &
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
    trange = 0.5d0*(maxt-mint)
    ! daily canopy conductance
    gc = abs(deltaWP)**(hydraulic_exponent)/((hydraulic_temp_coef*Rtot+trange))
    ! maximum rate of temperature and nitrogen (canopy efficiency) limited photosynthesis (gC.m-2.day-1)
    pn = lai*nit*NUE*exp(temp_exponent*maxt)
    ! pp and qq represent limitation by diffusion and metabolites respecitively
    pp = pn/gc ; qq = co2_comp_point-co2_half_sat
    ! calculate internal CO2 concentration (ppm)
    ci = 0.5d0*(co2+qq-pp+((co2+qq-pp)**2d0-4d0*(co2*qq-pp*co2_comp_point))**0.5d0)
    ! limit maximum quantium efficiency by leaf area, hyperbola
    e0 = lai_coef*lai**2d0/(lai**2d0+lai_const)
    ! calculate day length (hours)
!    dec = - asin( sin( 23.45d0 * pi / 180d0 ) * cos( 2d0 * pi * ( doy + 10d0 ) /365d0 ) )
!    sinld = sin( lat*deg_to_rad ) * sin( dec )
!    cosld = cos( lat*deg_to_rad ) * cos( dec )
!    aob = max(-1d0,min(1d0,sinld / cosld))
!    dayl = 12d0 * ( 1d0 + 2d0 * asin( aob ) / pi )

!--------------------------------------------------------------
    ! calculate day length (hours - not really hours)
    ! This is the old REFLEX project calculation 
    dec = -23.4*cos( (360d0*(doy+10d0)/365d0)*deg_to_rad ) * deg_to_rad
    mult = tan(lat*deg_to_rad) * tan(dec)
    if (mult >= 1d0) then
        dayl = 24d0
    else if (mult <= -1.0) then
        dayl = 0d0
    else
        dayl = 24d0*acos(-mult)/pi
    end if
! ---------------------------------------------------------------
    ! calculate CO2 limited rate of photosynthesis
    pd = gc*(co2-ci)
    ! calculate combined light and CO2 limited photosynthesis
    cps = e0*radiation*pd/(e0*radiation+pd)
    ! correct for day length variation
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
