
subroutine rdalec2(output_dim,MTT_dim,SS_dim &
                  ,met,pars &
                  ,out_var1,out_var2,out_var3 &
                  ,lat,nopars,nomet &
                  ,nofluxes,nopools,nodays,noyears,deltat &
                  ,nos_iter)

  use CARBON_MODEL_MOD, only: CARBON_MODEL,cica_time

  ! subroutine specificially deals with the calling of the fortran code model by
  ! R

  !!!!!!!!!!!
  ! Authorship contributions
  !
  ! This code is by:
  ! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
  ! modified for DALEC_GRASS by S. Zhu (University of Edinburgh)
  ! See function / subroutine specific comments for exceptions and contributors
  !!!!!!!!!!!

  implicit none
  ! declare input variables
  integer, intent(in) :: nopars         & ! number of parameters in vector
                        ,output_dim     & ! number of outputted variables
                        ,MTT_dim        & ! number of pools mean transit time estimates
                        ,SS_dim         & ! number of pools the steady state will be output for
                        ,nos_iter       & ! number of iterations
                        ,noyears        & ! number of years simulated
                        ,nomet          & ! number of meteorological fields
                        ,nofluxes       & ! number of model fluxes
                        ,nopools        & ! number of model pools
                        ,nodays           ! number of days in simulation

  double precision, intent(inout) :: deltat(nodays)     ! time step in decimal days
  double precision, intent(in) :: met(nomet,nodays)   & ! met drivers, note reverse of needed
                             ,pars(nopars,nos_iter)   & ! number of parameters
                             ,lat                       ! site latitude (degrees)

  ! output declaration
  double precision, intent(out), dimension(nos_iter,nodays,output_dim) :: out_var1
  double precision, intent(out), dimension(nos_iter,MTT_dim) :: out_var2  ! Mean annual MRT (years)
  double precision, intent(out), dimension(nos_iter,SS_dim) :: out_var3  ! Steady State (gC/m2)

  ! local variables
  ! vector of ecosystem pools
  integer :: i, y, y_s, y_e, nos_years, steps_per_year
  integer, dimension(nodays) :: lab_hak, fol_hak, root_hak, wood_hak, lit_hak, som_hak
  double precision, dimension((nodays+1),nopools) :: POOLS
  ! vector of ecosystem fluxes
  double precision, dimension(nodays,nofluxes) :: FLUXES
  double precision, dimension(nodays) :: lai & ! leaf area index
                                        ,GPP & ! Gross primary productivity
                                        ,NEE & ! net ecosystem exchange of CO2
                                 ,lab_filter &
                                 ,fol_filter &
                                ,root_filter &
                                ,wood_filter &
                                 ,lit_filter &
                                 ,som_filter
  double precision, dimension(nodays,2) :: REMOVED_C

  ! zero initial conditions
  lai = 0d0 ; GPP = 0d0 ; NEE = 0d0 ; POOLS = 0d0 ; FLUXES = 0d0
  out_var1 = 0d0 ; out_var2 = 0d0 ; out_var3 = 0d0

  ! generate deltat step from input data
  deltat(1) = met(1,1)
  do i = 2, nodays
     deltat(i) = met(1,i)-met(1,(i-1))
  end do
  ! number of years in analysis
  nos_years = nint(sum(deltat)/365.25d0)
  ! number of time steps per year
  steps_per_year = nodays/nos_years

  ! begin iterations
  do i = 1, nos_iter

     ! call the models
     call CARBON_MODEL(1,nodays,met,pars(1:nopars,i),deltat,nodays &
                      ,lat,lai,NEE,FLUXES,POOLS &
                      ,nopars,nomet,nopools,nofluxes,GPP,REMOVED_C)
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
     out_var1(i,1:nodays,1)  = FLUXES(1:nodays,1)       ! GPP (gC/m2/day)
     out_var1(i,1:nodays,2)  = FLUXES(1:nodays,3)       ! Rauto (gC/m2/day)
     out_var1(i,1:nodays,3)  = FLUXES(1:nodays,13)      ! Rhet_litter (gC/m2/day)
     out_var1(i,1:nodays,4)  = FLUXES(1:nodays,14)      ! Rhet_som (gC/m2/day)
     out_var1(i,1:nodays,5)  = FLUXES(1:nodays,17)      ! Total fire (gC/m2/day)
     out_var1(i,1:nodays,6)  = FLUXES(1:nodays,29)      ! Total harvested material (gC/m2/day)
     ! C internal fluxes (gC/m2/day)
     out_var1(i,1:nodays,7)  = FLUXES(1:nodays,4)       ! allocation to foliage (gC/m2/day)
     out_var1(i,1:nodays,8)  = FLUXES(1:nodays,5)       ! allocation to labile (gC/m2/day)
     out_var1(i,1:nodays,9)  = FLUXES(1:nodays,6)       ! allocation to fine roots (gC/m2/day)
     out_var1(i,1:nodays,10) = FLUXES(1:nodays,7)       ! allocation to wood (gC/m2/day)
     out_var1(i,1:nodays,11) = FLUXES(1:nodays,8)       ! labile to foliage (gC/m2/day)
     out_var1(i,1:nodays,12) = FLUXES(1:nodays,10)      ! foliage to litter (gC/m2/day)
     out_var1(i,1:nodays,13) = FLUXES(1:nodays,12)      ! fine root to litter (gC/m2/day)
     out_var1(i,1:nodays,14) = FLUXES(1:nodays,11)      ! wood to som (gC /m2/day)
     out_var1(i,1:nodays,15) = FLUXES(1:nodays,15)      ! litter to som (gC/m2/day)
     ! C disturbance fluxes (gC/m2/day)
     out_var1(i,1:nodays,16) = 0d0 ! FLUXES(1:nodays,18)      ! fire emission from labile (gC/m2/day)
     out_var1(i,1:nodays,17) = 0d0 ! FLUXES(1:nodays,24)      ! fire induced litter from labile (gC/m2/day)
     out_var1(i,1:nodays,18) = 0d0 ! FLUXES(1:nodays,19)      ! fire emission from foliage (gC/m2/day)
     out_var1(i,1:nodays,19) = 0d0 ! FLUXES(1:nodays,25)      ! fire induced litter from foliage (gC/m2/day)
     out_var1(i,1:nodays,20) = 0d0 ! FLUXES(1:nodays,20)      ! fire emission from fine roots (gC/m2/day)
     out_var1(i,1:nodays,21) = 0d0 ! FLUXES(1:nodays,26)      ! fire induced litter from fine roots (gC/m2/day)
     out_var1(i,1:nodays,22) = 0d0 ! FLUXES(1:nodays,21)      ! fire emission from wood (gC/m2/day)
     out_var1(i,1:nodays,23) = 0d0 ! FLUXES(1:nodays,27)      ! fire induced litter from wood (gC/m2/day)
     out_var1(i,1:nodays,24) = 0d0 ! FLUXES(1:nodays,22)      ! fire emission from litter (gC/m2/day)
     out_var1(i,1:nodays,25) = 0d0 ! FLUXES(1:nodays,28)      ! fire induced litter from litter (gC/m2/day)
     out_var1(i,1:nodays,26) = 0d0 ! FLUXES(1:nodays,23)      ! fire emission from som (gC/m2/day)
     out_var1(i,1:nodays,27) = 0d0 ! FLUXES(1:nodays,30)      ! harvest extracted from labile (gC/m2/day)
     out_var1(i,1:nodays,28) = 0d0 ! FLUXES(1:nodays,31)      ! harvest extracted from foliage (gC/m2/day)
     out_var1(i,1:nodays,29) = 0d0 ! FLUXES(1:nodays,32)      ! harvest extracted from fine roots (gC/m2/day)
     out_var1(i,1:nodays,30) = 0d0 ! FLUXES(1:nodays,33)      ! harvest extracted from wood (gC/m2/day)
     out_var1(i,1:nodays,31) = 0d0 ! FLUXES(1:nodays,34)      ! harvest extracted from litter (gC/m2/day)
     out_var1(i,1:nodays,32) = 0d0 ! FLUXES(1:nodays,35)      ! harvest extracted from som (gC/m2/day)
     out_var1(i,1:nodays,33) = 0d0 ! FLUXES(1:nodays,36)      ! harvest litter / residue from labile (gC/m2/day)
     out_var1(i,1:nodays,34) = 0d0 ! FLUXES(1:nodays,37)      ! harvest litter / residue from foliage (gC/m2/day)
     out_var1(i,1:nodays,35) = 0d0 ! FLUXES(1:nodays,38)      ! harvest litter / residue from fine roots (gC/m2/day)
     out_var1(i,1:nodays,36) = 0d0 ! FLUXES(1:nodays,39)      ! harvest litter / residue from wood (gC/m2/day)
     ! C pools (gC/m2)
     out_var1(i,1:nodays,37) = POOLS(1:nodays,1)        ! labile (gC/m2)
     out_var1(i,1:nodays,38) = POOLS(1:nodays,2)        ! foliage (gC/m2)
     out_var1(i,1:nodays,39) = POOLS(1:nodays,3)        ! fine root (gC/m2)
     out_var1(i,1:nodays,40) = POOLS(1:nodays,4)        ! wood (gC/m2)
     out_var1(i,1:nodays,41) = POOLS(1:nodays,5)        ! litter (gC/m2)
     out_var1(i,1:nodays,42) = POOLS(1:nodays,6)        ! som (gC/m2)
     ! Canopy (phenology) properties
     out_var1(i,1:nodays,43) = lai                      ! LAI (m2/m2)
     ! Photosynthesis / C~water coupling related
     out_var1(i,1:nodays,44) = cica_time                ! ratio of leaf internal to external CO2

     out_var1(i,1:nodays,45) = FLUXES(1:nodays,19)      ! animal manure C soil input (per time step)
     out_var1(i,1:nodays,46) = FLUXES(1:nodays,20)      ! animal resp co2 (per time step)
     out_var1(i,1:nodays,47) = FLUXES(1:nodays,21)      ! animal ch4 (per time step)

     !!!
     ! Estimate residence time information
     !!!

     ! Determine locations of zeros in pools to correct turnover calculation
     ! Labile
     lab_hak = 0 ; lab_filter(1:nodays) = 1d0
     where (POOLS(1:nodays,1) == 0) ! protection against NaN from division by zero
           lab_hak = 1 ; lab_filter(1:nodays) = 0d0
     end where
     ! Foliage
     fol_hak = 0 ; fol_filter(1:nodays) = 1d0
     where (POOLS(1:nodays,2) == 0) ! protection against NaN from division by zero
           fol_hak = 1 ; fol_filter(1:nodays) = 0d0
     end where
     ! Fine roots
     root_hak = 0 ; root_filter(1:nodays) = 1d0
     where (POOLS(1:nodays,3) == 0) ! protection against NaN from division by zero
           root_hak = 1 ; root_filter(1:nodays) = 0d0
     end where
     ! Wood
     wood_hak = 0 ; wood_filter(1:nodays) = 1d0
     where (POOLS(1:nodays,4) == 0) ! protection against NaN from division by zero
            wood_hak = 1 ; wood_filter(1:nodays) = 0d0
     end where
     ! Foliage + fine root litter
     lit_hak = 0 ; lit_filter(1:nodays) = 1d0
     where (POOLS(1:nodays,5) == 0) ! protection against NaN from division by zero
            lit_hak = 1 ; lit_filter(1:nodays) = 0d0
     end where
     ! Soil
     som_hak = 0 ; som_filter(1:nodays) = 1d0
     where (POOLS(1:nodays,6) == 0) ! protection against NaN from division by zero
           som_hak = 1 ; som_filter(1:nodays) = 0d0
     end where

     ! Estimate MRT (years)
     ! Labile
     out_var2(i,1) = sum( ((FLUXES(1:nodays,8) + & ! Natural losses
                            FLUXES(1:nodays,18) + FLUXES(1:nodays,24) + & ! Fire losses
                            FLUXES(1:nodays,30) + FLUXES(1:nodays,36) + & ! Harvest losses
                            FLUXES(1:nodays,22))                        & ! Cutting/grazing losses
                          / POOLS(1:nodays,1)) * lab_filter) / dble(nodays-sum(lab_hak))
     ! Foliage
     out_var2(i,2) = sum( ((FLUXES(1:nodays,10) + &
                            FLUXES(1:nodays,19) + FLUXES(1:nodays,25) + &
                            FLUXES(1:nodays,31) + FLUXES(1:nodays,37) + & 
                            FLUXES(1:nodays,23)) &
                          / POOLS(1:nodays,2)) * fol_filter) / dble(nodays-sum(fol_hak))
     ! Fine roots
     out_var2(i,3) = sum( ((FLUXES(1:nodays,12) + &
                            FLUXES(1:nodays,20) + FLUXES(1:nodays,26) + &
                            FLUXES(1:nodays,32) + FLUXES(1:nodays,38)) &
                          / POOLS(1:nodays,3)) * root_filter) / dble(nodays-sum(root_hak))
     ! Wood
     out_var2(i,4) = sum( ((FLUXES(1:nodays,11)+ &
                            FLUXES(1:nodays,21) + FLUXES(1:nodays,27) + &
                            FLUXES(1:nodays,33) + FLUXES(1:nodays,39)) &
                          / POOLS(1:nodays,4)) * wood_filter) / dble(nodays-sum(wood_hak))
     ! Litter (foliage+fine roots)
     out_var2(i,5) = sum( ((FLUXES(1:nodays,13) + FLUXES(1:nodays,15) + &
                            FLUXES(1:nodays,22) + FLUXES(1:nodays,28) + &
                            FLUXES(1:nodays,34)) &
                          / POOLS(1:nodays,5)) * lit_filter) / dble(nodays-sum(lit_hak))
     ! Soil
     out_var2(i,6) = sum( ((FLUXES(1:nodays,14) + FLUXES(1:nodays,23) + FLUXES(1:nodays,35)) &
                          / POOLS(1:nodays,6)) * som_filter) / dble(nodays-sum(som_hak))

     !
     ! Estimate pool inputs needed for steady state calculation
     !

     ! Once the canopy has closed the inputs to the live biomass are stable
     ! and can thus be estimated from the simulated inputs
     out_var3(i,1) = sum(FLUXES(:,5)) ! Labile
     out_var3(i,2) = sum(FLUXES(:,4)+FLUXES(:,8)) ! Foliage
     out_var3(i,3) = sum(FLUXES(:,6)) ! Fine root
     out_var3(i,4) = sum(FLUXES(:,7)) ! Wood
     out_var3(i,5) = sum(FLUXES(:,10)+FLUXES(:,12)+ &
                         FLUXES(:,24)+FLUXES(:,25)+FLUXES(:,26)+ &
                         FLUXES(:,36)+FLUXES(:,37)+FLUXES(:,38)) ! litter (foliage + roots)
     ! While foliar and fine root litter can be reasonably estimated directly (above),
     ! soil C inputs are still changing as the wood pool is not in steady state.
     ! Therefore, at this point we can account for disturbance inputs but NOT wood.
     ! The wood input is estimated later based on the steady state its steady state estimate
     out_var3(i,6) = sum(FLUXES(:,15)+FLUXES(:,27)+FLUXES(:,28)+FLUXES(:,39)) ! som

  end do ! nos_iter loop

  ! MTT - Convert daily fractional loss to years
  out_var2 = (out_var2*365.25d0)**(-1d0) ! iter,(lab,fol,root,wood,lit,som)

  ! Steady state gC/m2 estimation
  ! Determine the mean annual input (gC/m2/yr) based on current inputs for all pool,
  ! litter and soil pools updated below...
  out_var3 = (out_var3 / dble(nodays)) * 365.25d0 ! convert to annual mean input
  ! Then estimate the labile, foliar, fine root, wood and litter steady states.
  out_var3(:,1:5) = out_var3(:,1:5) * out_var2(:,1:5) ! multiply by residence time in years
  ! Using the wood SS estimate (gC/m2) the steady state input to the som litter pool...
  out_var3(:,6) = (out_var3(:,6) + (out_var3(:,4) / out_var2(:,4))) * out_var2(:,6)

  ! return back to the subroutine then
  return

end subroutine rdalec2
