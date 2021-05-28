
subroutine rdaleccdeaacm2bucketrmrgcwd(output_dim,aNPP_dim,MTT_dim,SS_dim,met,pars &
                               ,out_var,out_var2,out_var3,out_var4,out_var5 &
                               ,out_var6,lat,nopars,nomet &
                               ,nofluxes,nopools,pft,pft_specific,nodays,noyears,deltat &
                               ,nos_iter,soil_frac_clay_in,soil_frac_sand_in)

  use CARBON_MODEL_MOD, only: CARBON_MODEL, extracted_C, wSWP_time &
                             ,soil_frac_clay, soil_frac_sand, nos_soil_layers &
                             ,gs_demand_supply_ratio, cica_time &
                             ,gs_total_canopy, gb_total_canopy, canopy_par_MJday_time

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
                        ,pft            & ! plant functional type
                        ,output_dim     & !
                        ,aNPP_dim       & ! NPP allocation fraction variable dimension
                        ,MTT_dim        &
                        ,SS_dim         &
                        ,pft_specific   & !
                        ,nos_iter       & !
                        ,noyears        &
                        ,nomet          & ! number of meteorological fields
                        ,nofluxes       & ! number of model fluxes
                        ,nopools        & ! number of model pools
                        ,nodays           ! number of days in simulation

  double precision, intent(inout) :: deltat(nodays)     ! time step in decimal days
  double precision, intent(in) :: met(nomet,nodays)   & ! met drivers, note reverse of needed
                  ,soil_frac_clay_in(nos_soil_layers) & ! clay in soil (%)
                  ,soil_frac_sand_in(nos_soil_layers) & ! sand in soil (%)
                       ,pars(nopars,nos_iter)         & ! number of parameters
                       ,lat                 ! site latitude (degrees)

  ! output declaration
  double precision, intent(out), dimension(nos_iter,nodays,output_dim) :: out_var
  double precision, intent(out), dimension(nos_iter,aNPP_dim) :: out_var2 ! Mean annual NPP allocatino (0-1)
  double precision, intent(out), dimension(nos_iter,MTT_dim) :: out_var3  ! Mean annual MRT (years)
  double precision, intent(out), dimension(nos_iter,SS_dim) :: out_var4   ! Steady State (gC/m2)
  double precision, intent(out), dimension(nos_iter,MTT_dim,noyears) :: out_var5 ! Annual estimates of MRT (years)
  double precision, intent(out), dimension(nos_iter,MTT_dim) :: out_var6  ! Natural component of mean annual MRT (years)

  ! local variables
  ! vector of ecosystem pools
  integer :: i, y, y_s, y_e, nos_years, steps_per_year
  integer, dimension(nodays) :: fol_hak, root_hak, wood_hak, lit_hak, som_hak
  double precision, dimension(nos_iter) :: litwood_to_som_frac
  double precision, dimension((nodays+1),nopools) :: POOLS
  ! vector of ecosystem fluxes
  double precision, dimension(nodays,nofluxes) :: FLUXES
  double precision :: sumNPP,airt_adj
  double precision, dimension(nodays) :: lai & ! leaf area index
                                        ,GPP & ! Gross primary productivity
                                        ,NEE & ! net ecosystem exchange of CO2
                                 ,fol_filter &
                                ,root_filter &
                                ,wood_filter &
                                 ,lit_filter &
                                 ,som_filter

  ! zero initial conditions
  lai = 0d0 ; GPP = 0d0 ; NEE = 0d0 ; POOLS = 0d0 ; FLUXES = 0d0 ; out_var = 0d0

  ! update settings
  if (allocated(extracted_C)) deallocate(extracted_C)
  allocate(extracted_C(nodays))

  ! update soil parameters
  soil_frac_clay = soil_frac_clay_in
  soil_frac_sand = soil_frac_clay_in

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
     ! reset harvest variable
     extracted_C = 0d0
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

     ! now allocate the output the our 'output' variable
     out_var(i,1:nodays,1)  = lai
     out_var(i,1:nodays,2)  = GPP
     out_var(i,1:nodays,3)  = FLUXES(1:nodays,3) ! auto resp
     out_var(i,1:nodays,4)  = FLUXES(1:nodays,13) + FLUXES(1:nodays,14) + FLUXES(1:nodays,30)! het resp
     out_var(i,1:nodays,5)  = NEE
     out_var(i,1:nodays,6)  = POOLS(1:nodays,4) ! wood
     out_var(i,1:nodays,7)  = POOLS(1:nodays,6) ! som
     out_var(i,1:nodays,8)  = POOLS(1:nodays,1) + POOLS(1:nodays,2) + POOLS(1:nodays,3) & ! common pools
                            + POOLS(1:nodays,4) !+ POOLS(1:nodays,5) + POOLS(1:nodays,6)
     out_var(i,1:nodays,9)  = POOLS(1:nodays,3) ! root
     out_var(i,1:nodays,10) = POOLS(1:nodays,5) ! litter
     out_var(i,1:nodays,11) = POOLS(1:nodays,1) ! labile
     out_var(i,1:nodays,12) = POOLS(1:nodays,2) ! foliage
     out_var(i,1:nodays,13) = extracted_C(1:nodays) ! harvested material
     out_var(i,1:nodays,14) = FLUXES(1:nodays,17) ! Fire flux

     out_var(i,1:nodays,18) = FLUXES(1:nodays,29) ! Evapotranspiration (kgH2O.m-2.day-1)
     out_var(i,1:nodays,19) = POOLS(1:nodays,7)   ! rootwater (kgH2O.m-2.10cmdepth)
     out_var(i,1:nodays,20) = wSWP_time(1:nodays) ! Weighted Soil Water Potential (MPa)
     out_var(i,1:nodays,21) = gs_demand_supply_ratio ! ratio of evaporative demand over supply
     out_var(i,1:nodays,22) = gs_total_canopy(1:nodays)
     out_var(i,1:nodays,23) = canopy_par_MJday_time(1:nodays)
     out_var(i,1:nodays,24) = gb_total_canopy(1:nodays)
     out_var(i,1:nodays,25) = POOLS(1:nodays,8) ! litwood
     out_var(i,1:nodays,26) = cica_time(1:nodays)

     ! calculate the actual NPP allocation fractions to foliar, wood and fine root pools
     ! by comparing the sum alloaction to each pools over the sum NPP.
     sumNPP = sum(FLUXES(1:nodays,1)*(1-pars(2,i))) ! GPP * (1-Ra) fraction
     airt_adj = sum(met(3,1:nodays)) / dble(nodays)
     airt_adj = exp(pars(10,i)*airt_adj)
     out_var2(i,1) = sum(FLUXES(1:nodays,4)+FLUXES(1:nodays,8)) / sumNPP ! foliar
     out_var2(i,2) = sum(FLUXES(1:nodays,6)) / sumNPP ! fine root
     out_var2(i,3) = sum(FLUXES(1:nodays,7)) / sumNPP ! wood

     ! Determine locations of zeros in pools to correct turnover calculation
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
     ! Fol+root+wood litter
     lit_hak = 0 ; lit_filter(1:nodays) = 1d0
     where (POOLS(1:nodays,5) + POOLS(1:nodays,8) == 0) ! protection against NaN from division by zero
            lit_hak = 1 ; lit_filter(1:nodays) = 0d0
     end where
     ! Soil
     som_hak = 0 ; som_filter(1:nodays) = 1d0
     where (POOLS(1:nodays,6) == 0) ! protection against NaN from division by zero
           som_hak = 1 ; som_filter(1:nodays) = 0d0
     end where

     ! Mean transit times
     ! Foliage (/day)
     out_var3(i,1) = sum( ((FLUXES(1:nodays,10)+FLUXES(1:nodays,19)+FLUXES(1:nodays,25)) &
                          / POOLS(1:nodays,2)) * fol_filter) / dble(nodays-sum(fol_hak))
     ! Fine roots (/day)
     out_var3(i,2) = sum( ((FLUXES(1:nodays,12)+FLUXES(1:nodays,20)+FLUXES(1:nodays,26)) &
                          / POOLS(1:nodays,3)) * root_filter) / dble(nodays-sum(root_hak))
     ! Wood (/day)
     out_var3(i,3) = sum( ((FLUXES(1:nodays,11)+FLUXES(1:nodays,21)+FLUXES(1:nodays,27)) &
                             / POOLS(1:nodays,4)) * wood_filter) / dble(nodays-sum(wood_hak))
     ! Litter (fol+root+wood; /day)
     out_var3(i,4) = sum( ((FLUXES(1:nodays,13)+FLUXES(1:nodays,15)+FLUXES(1:nodays,22)+FLUXES(1:nodays,28) &
                           +FLUXES(1:nodays,30)+FLUXES(1:nodays,31)+FLUXES(1:nodays,32)+FLUXES(1:nodays,33)) &
                          / (POOLS(1:nodays,5)+POOLS(1:nodays,8))) * lit_filter) &
                   / dble(nodays-sum(lit_hak)) ! litter+wood litter (/day)
     ! Soil (/day)
     out_var3(i,5) = sum( ((FLUXES(1:nodays,14)+FLUXES(1:nodays,23)) &
                          / POOLS(1:nodays,6)) * som_filter) / dble(nodays-sum(som_hak))

     ! Mean transit times (natural only)
     ! Foliage (/day)
     out_var6(i,1) = sum( (FLUXES(1:nodays,10) / POOLS(1:nodays,2)) * fol_filter) / dble(nodays-sum(fol_hak))
     ! Fine roots (/day)
     out_var6(i,2) = sum( (FLUXES(1:nodays,12) / POOLS(1:nodays,3)) * root_filter) / dble(nodays-sum(root_hak))
     ! Wood (/day)
     out_var6(i,3) = sum( (FLUXES(1:nodays,11) / POOLS(1:nodays,4)) * wood_filter) / dble(nodays-sum(wood_hak))
     ! Litter (fol+root+wood; /day)
     out_var6(i,4) = sum( ((FLUXES(1:nodays,13)+FLUXES(1:nodays,15)+FLUXES(1:nodays,30)+FLUXES(1:nodays,31)) &
                          / (POOLS(1:nodays,5)+POOLS(1:nodays,8))) * lit_filter) &
                   / dble(nodays-sum(lit_hak)) ! litter+wood litter (/day)
     ! Soil (/day)
     out_var6(i,5) = sum( (FLUXES(1:nodays,14) &
                          / POOLS(1:nodays,6)) * som_filter) / dble(nodays-sum(som_hak))

     ! Keep track of the fraction of wood litter transfer to som, this value is needed for the steady state estimation
     litwood_to_som_frac(i) = sum( (FLUXES(1:nodays,31) / POOLS(1:nodays,8)) * lit_filter) / dble(nodays-sum(lit_hak))

     ! Now loop through each year to estimate the annual residence time
     do y = 1, nos_years
        ! Estimate time steps covered by this year
        y_s = 1 + (steps_per_year * (y-1)) ; y_e = steps_per_year * y
        ! Foliage
        out_var5(i,1,y) = sum( ((FLUXES(y_s:y_e,10)+FLUXES(y_s:y_e,19)+FLUXES(y_s:y_e,25)) &
                                / POOLS(y_s:y_e,2)) * fol_filter(y_s:y_e)) / dble(steps_per_year-sum(fol_hak(y_s:y_e)))
        ! Fine roots
        out_var5(i,2,y) = sum( ((FLUXES(y_s:y_e,12)+FLUXES(y_s:y_e,20)+FLUXES(y_s:y_e,26)) &
                                / POOLS(y_s:y_e,3)) * root_filter(y_s:y_e)) / dble(steps_per_year-sum(root_hak(y_s:y_e)))
        ! Wood
        out_var5(i,3,y) = sum( ((FLUXES(y_s:y_e,11)+FLUXES(y_s:y_e,21)+FLUXES(y_s:y_e,27)) &
                                / POOLS(y_s:y_e,4)) * wood_filter(y_s:y_e)) / dble(steps_per_year-sum(wood_hak(y_s:y_e)))
        ! Litter (fol+fine root+wood)
        out_var5(i,4,y) = sum( ((FLUXES(y_s:y_e,13)+FLUXES(y_s:y_e,15)+FLUXES(y_s:y_e,22)+FLUXES(y_s:y_e,28) &
                                +FLUXES(y_s:y_e,30)+FLUXES(y_s:y_e,31)+FLUXES(y_s:y_e,32)+FLUXES(y_s:y_e,33)) &
                                / (POOLS(y_s:y_e,5)+POOLS(y_s:y_e,8))) * lit_filter(y_s:y_e)) &
                        / dble(steps_per_year-sum(lit_hak(y_s:y_e)))
        ! Soil
        out_var5(i,5,y) = sum( ((FLUXES(y_s:y_e,14)+FLUXES(y_s:y_e,23)) &
                                / POOLS(y_s:y_e,6)) * som_filter(y_s:y_e)) / dble(steps_per_year-sum(som_hak(y_s:y_e)))
     end do

     ! Once the canopy has closed the inputs to the live biomass are stable
     ! and can thus be estimated from the simulated inputs. Similarly the litter pool input
     ! i.e. foliage and fine roots are likely to be in steady state.
     out_var4(i,1) = sum(FLUXES(:,4)+FLUXES(:,8)) ! Foliage
     out_var4(i,2) = sum(FLUXES(:,6)) ! Fine root
     out_var4(i,3) = sum(FLUXES(:,7)) ! Wood
     ! Foliar and fine root litter can likewise be estimated directly,
     ! however wood litter and soil C inputs are still changing as the wood pool is not in steady state.
     ! Therefore, at this point we can account for the fol, fine root, and disturbance inputs but NOT wood.
     ! The wood input is estimated later based on the steady state its steady state estimate
     out_var4(i,4) = sum(FLUXES(1:nodays,10)+FLUXES(1:nodays,12) &
                        +FLUXES(1:nodays,24)+FLUXES(1:nodays,25)+FLUXES(1:nodays,26)) ! lit
     out_var4(i,5) = sum(FLUXES(:,15)+FLUXES(1:nodays,27)+FLUXES(1:nodays,28)) ! som

  end do ! nos_iter loop

  ! MTT - Convert daily fractional loss to years
  out_var3 = (out_var3*365.25d0)**(-1d0) ! iter,(fol,root,wood,lit+litwood,som)
  out_var5 = (out_var5*365.25d0)**(-1d0) ! iter,(fol,root,wood,lit+litwood,som)
  out_var6 = (out_var6*365.25d0)**(-1d0) ! iter,(fol,root,wood,lit+litwood,som)

  ! Steady state gC/m2 estimation
  ! Determine the mean annual input (gC/m2/yr) based on current inputs for all pool,
  ! litter and soil pools updated below...
  out_var4 = (out_var4 / dble(nodays)) * 365.25d0 ! convert to annual mean input
  ! Then estimate the foliar, fine root, wood and litter steady states.
  out_var4(:,1:3) = out_var4(:,1:3) * out_var3(:,1:3) ! multiply by residence time in years
  ! Using the wood SS estimate (gC/m2) the steady state input to the litter+wood litter pool...
  out_var4(:,4) = (out_var4(:,4) + (out_var4(:,3) / out_var3(:,3))) * out_var3(:,4)
  ! ...which is then in turn used to update the soil pool
  ! NOTE: that because not all wood litter
  out_var4(:,5) = (out_var4(:,5) + ((out_var4(:,4) / out_var3(:,4))*litwood_to_som_frac) ) * out_var3(:,5)

  ! deallocate harvested variable
  deallocate(extracted_C)

  ! return back to the subroutine then
  return

end subroutine rdaleccdeaacm2bucketrmrgcwd
