
subroutine rdalec31(output_dim,MTT_dim,SS_dim &
                   ,met,pars &
                   ,out_var1,out_var2,out_var3,out_var4,out_var5 &
                   ,lat,nopars,nomet &
                   ,nofluxes,nopools,nodays,nos_years,deltat &
                   ,nos_iter,soil_frac_clay_in,soil_frac_sand_in)

  use CARBON_MODEL_MOD, only: CARBON_MODEL, wSWP_time &
                             ,soil_frac_clay, soil_frac_sand, nos_soil_layers &
                             ,gs_demand_supply_ratio, cica_time &
                             ,gs_total_canopy, gb_total_canopy &
                             ,canopy_par_MJday_time, root_depth_time &
                             ,snow_storage_time

  ! subroutine specificially deals with the calling of the fortran code model by
  ! R

  !!!!!!!!!!!
  ! Authorship contributions
  !
  ! This code is by:
  ! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
  ! See function / subroutine specific comments for exceptions and contributors
  !!!!!!!!!!!

  implicit none
  ! declare input variables
  integer, intent(in) :: nopars         & ! number of paremeters in vector
                        ,output_dim     & !
                        ,MTT_dim        & ! number of pools mean transit time estimates
                        ,SS_dim         & ! number of pools the steady state will be output for
                        ,nos_iter       & !
                        ,nomet          & ! number of meteorological fields
                        ,nofluxes       & ! number of model fluxes
                        ,nopools        & ! number of model pools
                        ,nodays         & ! number of time steps in simulation
                        ,nos_years        ! number of years in simulation

  double precision, intent(inout) :: deltat(nodays)     ! time step in decimal days
  double precision, intent(in) :: met(nomet,nodays)   & ! met drivers, note reverse of needed
                  ,soil_frac_clay_in(nos_soil_layers) & ! clay in soil (%)
                  ,soil_frac_sand_in(nos_soil_layers) & ! sand in soil (%)
                       ,pars(nopars,nos_iter)         & ! number of parameters
                       ,lat                 ! site latitude (degrees)

  ! output declaration
  double precision, intent(out), dimension(nos_iter,nodays,output_dim) :: out_var1
  double precision, intent(out), dimension(nos_iter,MTT_dim) :: out_var2  ! Mean annual MRT (years)
  double precision, intent(out), dimension(nos_iter,SS_dim) :: out_var3  ! Steady State (gC/m2)
  double precision, intent(out), dimension(nos_iter,output_dim) :: out_var4 ! Long term mean of out_var1
  double precision, intent(out), dimension(nos_iter,nos_years,output_dim) :: out_var5 ! Mean annual of out_var1

  ! local variables
  ! vector of ecosystem pools
  integer :: a, e, i, s, v, steps_per_year!, nos_years
  integer, dimension(nodays) :: pool_hak
  double precision, dimension((nodays+1),nopools) :: POOLS
  ! vector of ecosystem fluxes
  double precision, dimension(nodays,nofluxes) :: FLUXES
  double precision, dimension(nodays) :: tmp &
                                        ,lai & ! leaf area index
                                        ,GPP & ! Gross primary productivity
                                        ,NEE   ! net ecosystem exchange of CO2

  ! zero initial conditions
  lai = 0d0 ; GPP = 0d0 ; NEE = 0d0 ; POOLS = 0d0 ; FLUXES = 0d0
  out_var1 = 0d0 ; out_var2 = 0d0 ; out_var3 = 0d0

  ! update soil parameters
  soil_frac_clay(1:nos_soil_layers) = soil_frac_clay_in(1:nos_soil_layers)
  soil_frac_sand(1:nos_soil_layers) = soil_frac_sand_in(1:nos_soil_layers)

  ! generate deltat step from input data
  deltat(1) = met(1,1)
  do i = 2, nodays
     deltat(i) = met(1,i)-met(1,(i-1))
  end do
  ! number of years in analysis
  !nos_years = nint(sum(deltat)/365.25d0)
  ! number of time steps per year
  steps_per_year = nodays/nos_years

  ! begin iterations
  do i = 1, nos_iter

     ! call the models
     call CARBON_MODEL(1,nodays,met,pars(1:nopars,i),deltat,nodays &
                      ,lat,lai,NEE,FLUXES,POOLS &
                      ,nopars,nomet,nopools,nofluxes,GPP)
!if (i == 1) then
!    open(unit=666,file="/home/lsmallma/out.csv", &
!         status='replace',action='readwrite' )
!    write(666,*),"GSI",FLUXES(:,14)(1:365)
!    close(666)
!endif

     !
     ! Allocate the output the our 'output' variable
     !

     ! C ecosystem fluxes (gC/m2/day)
     out_var1(i,1:nodays,1)  = FLUXES(1:nodays,1)          ! GPP (gC/m2/day)
     out_var1(i,1:nodays,2)  = FLUXES(1:nodays,3)          ! Rauto (gC/m2/day)
     out_var1(i,1:nodays,3)  = FLUXES(1:nodays,13)         ! Rhet_litter (gC/m2/day)
     out_var1(i,1:nodays,4)  = FLUXES(1:nodays,14)         ! Rhet_som (gC/m2/day)
     out_var1(i,1:nodays,5)  = FLUXES(1:nodays,17)         ! Total fire (gC/m2/day)
     out_var1(i,1:nodays,6)  = FLUXES(1:nodays,30)         ! harvested material (gC/m2/day)
     ! C internal fluxes (gC/m2/day)
     out_var1(i,1:nodays,7)  = FLUXES(1:nodays,4)          ! allocation to foliage (gC/m2/day)
     out_var1(i,1:nodays,8)  = FLUXES(1:nodays,5)          ! allocation to labile (gC/m2/day)
     out_var1(i,1:nodays,9)  = FLUXES(1:nodays,6)          ! allocation to fine roots (gC/m2/day)
     out_var1(i,1:nodays,10) = FLUXES(1:nodays,7)          ! allocation to wood (gC/m2/day)
     out_var1(i,1:nodays,11) = FLUXES(1:nodays,8)          ! labile to foliage (gC/m2/day)
     out_var1(i,1:nodays,12) = FLUXES(1:nodays,10)         ! foliage scenesence (gC/m2/day)
     out_var1(i,1:nodays,13) = FLUXES(1:nodays,12)         ! fine root turnover (gC/m2/day)
     out_var1(i,1:nodays,14) = FLUXES(1:nodays,11)         ! wood turnover (gC /m2/day)
     out_var1(i,1:nodays,15) = FLUXES(1:nodays,15)         ! Decomp_litter (gC/m2/day)
     ! C disturbance fluxes (gC/m2/day)
     out_var1(i,1:nodays,16) = FLUXES(1:nodays,18)         ! fire emission from labile (gC/m2/day)
     out_var1(i,1:nodays,17) = FLUXES(1:nodays,24)         ! fire induced litter from labile (gC/m2/day)
     out_var1(i,1:nodays,18) = FLUXES(1:nodays,19)         ! fire emission from foliage (gC/m2/day)
     out_var1(i,1:nodays,19) = FLUXES(1:nodays,25)         ! fire induced litter from foliage (gC/m2/day)
     out_var1(i,1:nodays,20) = FLUXES(1:nodays,20)         ! fire emission from fine roots (gC/m2/day)
     out_var1(i,1:nodays,21) = FLUXES(1:nodays,26)         ! fire induced litter from fine roots (gC/m2/day)
     out_var1(i,1:nodays,22) = FLUXES(1:nodays,21)         ! fire emission from wood (gC/m2/day)
     out_var1(i,1:nodays,23) = FLUXES(1:nodays,27)         ! fire induced litter from wood (gC/m2/day)
     out_var1(i,1:nodays,24) = FLUXES(1:nodays,22)         ! fire emission from litter (gC/m2/day)
     out_var1(i,1:nodays,25) = FLUXES(1:nodays,28)         ! fire induced litter from litter (gC/m2/day)
     out_var1(i,1:nodays,26) = FLUXES(1:nodays,23)         ! fire emission from som (gC/m2/day)
     out_var1(i,1:nodays,27) = FLUXES(1:nodays,31)         ! harvest extracted from labile (gC/m2/day)
     out_var1(i,1:nodays,28) = FLUXES(1:nodays,32)         ! harvest extracted from foliage (gC/m2/day)
     out_var1(i,1:nodays,29) = FLUXES(1:nodays,33)         ! harvest extracted from fine roots (gC/m2/day)
     out_var1(i,1:nodays,30) = FLUXES(1:nodays,34)         ! harvest extracted from wood (gC/m2/day)
     out_var1(i,1:nodays,31) = FLUXES(1:nodays,35)         ! harvest extracted from litter (gC/m2/day)
     out_var1(i,1:nodays,32) = FLUXES(1:nodays,36)         ! harvest extracted from som (gC/m2/day)
     out_var1(i,1:nodays,33) = FLUXES(1:nodays,37)         ! harvest litter / residue from labile (gC/m2/day)
     out_var1(i,1:nodays,34) = FLUXES(1:nodays,38)         ! harvest litter / residue from foliage (gC/m2/day)
     out_var1(i,1:nodays,35) = FLUXES(1:nodays,39)         ! harvest litter / residue from fine roots (gC/m2/day)
     out_var1(i,1:nodays,36) = FLUXES(1:nodays,40)         ! harvest litter / residue from wood (gC/m2/day)
     ! C pools (gC/m2)
     out_var1(i,1:nodays,37) = POOLS(1:nodays,1)           ! labile (gC/m2)
     out_var1(i,1:nodays,38) = POOLS(1:nodays,2)           ! foliage (gC/m2)
     out_var1(i,1:nodays,39) = POOLS(1:nodays,3)           ! fine root (gC/m2)
     out_var1(i,1:nodays,40) = POOLS(1:nodays,4)           ! wood (gC/m2)
     out_var1(i,1:nodays,41) = POOLS(1:nodays,5)           ! litter (gC/m2)
     out_var1(i,1:nodays,42) = POOLS(1:nodays,6)           ! som (gC/m2)
     ! Water cycle related
     out_var1(i,1:nodays,43) = FLUXES(1:nodays,29)         ! Evapotranspiration (kgH2O.m-2.day-1)
     out_var1(i,1:nodays,44) = FLUXES(1:nodays,41)         ! transpiration (kgH2O.m-2.day-1)
     out_var1(i,1:nodays,45) = FLUXES(1:nodays,42)         ! soil evaporation (kgH2O.m-2.day-1)
     out_var1(i,1:nodays,46) = FLUXES(1:nodays,43)         ! wet canopy evaporation (kgH2O.m-2.day-1)
     out_var1(i,1:nodays,47) = FLUXES(1:nodays,44)         ! runoff (kgH2O.m-2.day-1)
     out_var1(i,1:nodays,48) = FLUXES(1:nodays,45)         ! underflow (kgH2O.m-2.day-1)
     out_var1(i,1:nodays,49) = POOLS(1:nodays,7)           ! surface water (kgH2O.m-2.30cmdepth)
     out_var1(i,1:nodays,50) = wSWP_time(1:nodays)         ! Weighted Soil Water Potential (MPa)
     out_var1(i,1:nodays,51) = snow_storage_time(1:nodays) ! Snow storage (kgH2O/m2)
     ! Canopy (phenology) properties
     out_var1(i,1:nodays,52) = lai                         ! LAI (m2/m2)
     ! Photosynthesis / C~water coupling related
     out_var1(i,1:nodays,53) = gs_demand_supply_ratio      ! ratio of evaporative demand over supply
     out_var1(i,1:nodays,54) = gs_total_canopy             ! Canopy scale stomatal conductance during day light (mmolH2O/m2ground/s)
     out_var1(i,1:nodays,55) = canopy_par_MJday_time       ! Canopy absorbed PAR (MJ/m2ground/day)
     out_var1(i,1:nodays,56) = gb_total_canopy             ! Canopy scale aerodynamic conductance (mmolH2O/m2ground/s)
     out_var1(i,1:nodays,57) = cica_time                   ! ratio of leaf internal to external CO2
     ! misc
     out_var1(i,1:nodays,58) = root_depth_time             ! rooting depth (m)

     !
     ! Calculate long-term mean of out_var1
     !
     
     ! Loop across each variable
     do v = 1, output_dim
        ! Calculate mean value
        out_var4(i,v) = sum(out_var1(i,1:nodays,v)) / dble(nodays)
     end do

     !
     ! Calculate the mean annual of out_var1
     !

     ! Calculate mean annual
     s = 1 ; e = steps_per_year
     do a = 1, nos_years
        do v = 1, output_dim
           out_var5(i,a,v) = sum(out_var1(i,s:e,v)) / dble(steps_per_year)
        end do
        ! Iterate counters
        s = s + steps_per_year ; e = s + steps_per_year - 1
     end do
          
     !!!
     ! Estimate residence time information
     !!!

     ! Labile
     ! Estimate MRT (years)
     pool_hak = 1 ; tmp = 0d0
     where (POOLS(1:nodays,1) > 0d0) ! protection against NaN from division by zero
            pool_hak = 0 
            tmp = ((FLUXES(1:nodays,8) + &
                    FLUXES(1:nodays,18) + FLUXES(1:nodays,24) + &
                    FLUXES(1:nodays,31) + FLUXES(1:nodays,37)) / POOLS(1:nodays,1))
     end where
     out_var2(i,1) = sum(tmp) / dble(nodays-sum(pool_hak))
     ! Foliage
     ! Estimate MRT (years)
     pool_hak = 1 ; tmp = 0d0
     where (POOLS(1:nodays,2) > 0d0) ! protection against NaN from division by zero
            pool_hak = 0 
            tmp = ((FLUXES(1:nodays,10) + &
                    FLUXES(1:nodays,19) + FLUXES(1:nodays,25) + &
                    FLUXES(1:nodays,32) + FLUXES(1:nodays,38)) / POOLS(1:nodays,2))
     end where
     out_var2(i,2) = sum(tmp) / dble(nodays-sum(pool_hak))
     ! Fine roots
     ! Estimate MRT (years)
     pool_hak = 1 ; tmp = 0d0
     where (POOLS(1:nodays,3) > 0d0) ! protection against NaN from division by zero
            pool_hak = 0 
            tmp = ((FLUXES(1:nodays,12) + &
                    FLUXES(1:nodays,20) + FLUXES(1:nodays,26) + &
                    FLUXES(1:nodays,33) + FLUXES(1:nodays,39)) / POOLS(1:nodays,3))
     end where
     out_var2(i,3) = sum(tmp) / dble(nodays-sum(pool_hak))
     ! Wood
     ! Estimate MRT (years)
     pool_hak = 1 ; tmp = 0d0
     where (POOLS(1:nodays,4) > 0d0) ! protection against NaN from division by zero
            pool_hak = 0 
            tmp = ((FLUXES(1:nodays,11)+ &
                    FLUXES(1:nodays,21) + FLUXES(1:nodays,27) + &
                    FLUXES(1:nodays,34) + FLUXES(1:nodays,40)) / POOLS(1:nodays,4))
     end where
     out_var2(i,4) = sum(tmp) / dble(nodays-sum(pool_hak))
     ! Foliage + fine root litter
     ! Estimate MRT (years)
     pool_hak = 1 ; tmp = 0d0
     where (POOLS(1:nodays,5) > 0d0) ! protection against NaN from division by zero
            pool_hak = 0 
            tmp = ((FLUXES(1:nodays,13) + FLUXES(1:nodays,15) + &
                    FLUXES(1:nodays,22) + FLUXES(1:nodays,28) + &
                    FLUXES(1:nodays,35)) / POOLS(1:nodays,5))
     end where
     out_var2(i,5) = sum(tmp) / dble(nodays-sum(pool_hak))
     ! Soil
     ! Estimate MRT (years)
     pool_hak = 1 ; tmp = 0d0
     where (POOLS(1:nodays,6) > 0d0) ! protection against NaN from division by zero
            pool_hak = 0 
            tmp = ((FLUXES(1:nodays,14) + FLUXES(1:nodays,23) + FLUXES(1:nodays,36)) &
                  / POOLS(1:nodays,6))
     end where
     out_var2(i,6) = sum(tmp) / dble(nodays-sum(pool_hak))

     !
     ! Estimate pool inputs needed for steady state calculation
     !

     ! Once the canopy has closes the inputs to the live biomass are stable
     ! and can thus be estimated from the simulated inputs
     out_var3(i,1) = sum(FLUXES(:,5)) ! Labile
     out_var3(i,2) = sum(FLUXES(:,4)+FLUXES(:,8)) ! Foliage
     out_var3(i,3) = sum(FLUXES(:,6)) ! Fine root
     out_var3(i,4) = sum(FLUXES(:,7)) ! Wood
     out_var3(i,5) = sum(FLUXES(:,10)+FLUXES(:,12)+ &
                         FLUXES(:,24)+FLUXES(:,25)+FLUXES(:,26)+ &
                         FLUXES(:,37)+FLUXES(:,38)+FLUXES(:,39)) ! litter (foliage + roots)
     ! While foliar and fine root litter can be reasonably estimated directly (above),
     ! soil C inputs are still changing as the wood pool is not in steady state.
     ! Therefore, at this point we can account for disturbance inputs (including wood)
     ! but NOT natural wood. The natural wood input is estimated later based on
     ! its steady state estimate
     out_var3(i,6) = sum(FLUXES(:,15)+FLUXES(:,27)+FLUXES(:,28)+FLUXES(:,40)) ! som

  end do ! nos_iter loop

  ! MTT - Convert daily fractional loss to years
  out_var2 = (out_var2*365.25d0)**(-1d0) ! iter,(lab,fol,root,wood,lit,som)

  ! Steady state gC/m2 estimation
  ! Determine the mean annual input (gC/m2/yr) based on current inputs for all pool,
  ! litter and soil pools updated below...
  out_var3 = (out_var3 / dble(nodays)) * 365.25d0 ! convert to annual mean input
  ! Then estimate the labile, foliar, fine root, wood and litter steady states.
  out_var3(1:nos_iter,1:5) = out_var3(1:nos_iter,1:5) * out_var2(1:nos_iter,1:5) ! multiply by residence time in years
  ! Using the wood SS estimate (gC/m2) the steady state input to the som litter pool...
  out_var3(1:nos_iter,6) = (out_var3(1:nos_iter,6) + &
                           (out_var3(1:nos_iter,4) / out_var2(1:nos_iter,4))) &
                         * out_var2(1:nos_iter,6)

  ! return back to the subroutine then
  return

end subroutine rdalec31
