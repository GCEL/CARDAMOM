

subroutine rdalec(output_dim,aNPP_dim,MTT_dim,SS_dim,met,pars &
                 ,out_var,out_var2,out_var3,out_var4,out_var5 &
                 ,out_var6,lat &
                 ,nopars,nomet,nofluxes,nopools,pft,pft_specific &
                 ,nodays,noyears,deltat,nos_iter,exepath,pathlength)

  use CARBON_MODEL_MOD, only: CARBON_MODEL, itemp, ivpd, iphoto, &
                              harvest_residue_to_litter, harvest_residue_to_litwood, &
                              harvest_residue_to_som, harvest_loss_litter, &
                              harvest_loss_litwood, harvest_loss_som,      &
                              harvest_loss_labile, harvest_loss_foliar,    &
                              harvest_loss_roots, harvest_loss_wood,       &
                              fire_loss_labile, fire_loss_foliar, fire_loss_roots, &
                              fire_loss_wood, fire_loss_litter, fire_loss_litwood, &
                              fire_loss_som, fire_residue_to_litter, &
                              fire_residue_to_litwood,fire_residue_to_som,       &
                              gs_demand_supply_ratio, cica_time, &
                              gs_total_canopy, gb_total_canopy, canopy_par_MJday_time
  use CARBON_MODEL_CROP_MOD, only: CARBON_MODEL_CROP

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
  interface
    subroutine crop_development_parameters(stock_seed_labile,DS_shoot,DS_root,fol_frac &
                                          ,stem_frac,root_frac,DS_LRLV,LRLV,DS_LRRT,LRRT &
                                          ,exepath,pathlength)
      implicit none
      ! declare inputs
      ! crop specific variables
      integer, intent(in) :: pathlength
      character(pathlength),intent(in) :: exepath
      double precision :: stock_seed_labile
      double precision, allocatable, dimension(:) :: DS_shoot, & !
                                                      DS_root, & !
                                                     fol_frac, & !
                                                    stem_frac, & !
                                                    root_frac, & !
                                                      DS_LRLV, & !
                                                         LRLV, & !
                                                      DS_LRRT, & !
                                                         LRRT
      ! local variables..
      integer :: columns, i, rows, input_crops_unit, ios
      character(225) :: variables,filename
    end subroutine crop_development_parameters
  end interface

  ! declare input variables
  integer, intent(in) :: pathlength
  character(pathlength), intent(in) :: exepath
  integer, intent(in) :: nopars         & ! number of paremeters in vector
                        ,output_dim     & !
                        ,aNPP_dim       & ! NPP allocation fraction variable dimension
                        ,MTT_dim        &
                        ,SS_dim         &
                        ,pft            & ! plant functional type
                        ,pft_specific   & !
                        ,nos_iter       & !
                        ,noyears        &
                        ,nomet          & ! number of meteorological fields
                        ,nofluxes       & ! number of model fluxes
                        ,nopools        & ! number of model pools
                        ,nodays           ! number of days in simulation

  double precision, intent(in) :: met(nomet,nodays)   & ! met drivers, note reverse of needed
                       ,pars(nopars,nos_iter)         & ! number of parameters
                       ,lat                 ! site latitude (degrees)

  double precision, intent(inout) :: deltat(nodays) ! time step in decimal days

  ! output declaration
  double precision, intent(out), dimension(nos_iter,nodays,output_dim) :: out_var
  double precision, intent(out), dimension(nos_iter,aNPP_dim) :: out_var2 ! Mean annual NPP allocatino (0-1)
  double precision, intent(out), dimension(nos_iter,MTT_dim) :: out_var3  ! Mean annual MRT (years)
  double precision, intent(out), dimension(nos_iter,SS_dim) :: out_var4   ! Steady State (gC/m2)
  double precision, intent(out), dimension(nos_iter,MTT_dim,noyears) :: out_var5 ! Annual estimates of MRT (years)
  double precision, intent(out), dimension(nos_iter,MTT_dim) :: out_var6  ! Natural component of mean annual MRT (years)

  ! local variables
  integer :: i, y, y_s, y_e, nos_years, steps_per_year
  integer, dimension(nodays) :: fol_hak, root_hak, wood_hak, lit_hak, som_hak
  double precision, dimension(nos_iter) :: litwood_to_som_frac
  ! vector of ecosystem pools
  double precision, dimension((nodays+1),nopools) :: POOLS
  ! vector of ecosystem fluxes
  double precision, dimension(nodays,nofluxes) :: FLUXES
  double precision, dimension(nodays) :: resid_fol
  integer, dimension(nodays) :: hak ! variable to determine number of NaN
  double precision :: sumNPP, fauto
  double precision, dimension(nodays) :: lai & ! leaf area index
                                        ,GPP & ! Gross primary productivity
                                        ,NEE & ! net ecosystem exchange of CO2
                                 ,fol_filter &
                                ,root_filter &
                                ,wood_filter &
                                 ,lit_filter &
                                 ,som_filter

  ! crop development parameters declared here. These are also found in
  ! MHMCMC_STRUCTURES PI%
  ! crop specific variables
  double precision :: stock_seed_labile
  double precision, allocatable, dimension(:)  ::  DS_shoot, & !
                                                    DS_root, & !
                                                   fol_frac, & !
                                                  stem_frac, & !
                                                  root_frac, & !
                                                    DS_LRLV, & !
                                                       LRLV, & !
                                                    DS_LRRT, & !
                                                       LRRT

! profiling example
!real :: begin, done,f1=0,f2=0,f3=0,f4=0,f5=0,total_time = 0
!real :: Rtot_track_time = 0, aero_time = 0 , soilwater_time = 0 , acm_et_time = 0 , Rm_time = 0
!call cpu_time(done)
!print*,"time taken per iter",(done-begin) / real(nos_iter)

  ! zero initial conditions
  lai = 0d0 ; GPP = 0d0 ; NEE = 0d0 ; POOLS = 0d0 ; FLUXES = 0d0 ; out_var = 0d0

  ! generate deltat step from input data
  deltat(1) = met(1,1)
  do i = 2, nodays
     deltat(i) = met(1,i)-met(1,(i-1))
  end do
  ! number of years in analysis
  nos_years = nint(sum(deltat)/365.25d0)
  ! number of time steps per year
  steps_per_year = nodays/nos_years

  ! when crop model in use should load crop development parameters here
  ! modifications neede....
  if (pft == 1) call crop_development_parameters(stock_seed_labile,DS_shoot,DS_root,fol_frac &
                                                ,stem_frac,root_frac,DS_LRLV,LRLV,DS_LRRT,LRRT &
                                                ,exepath,pathlength)

  ! begin iterations
  do i = 1, nos_iter
     ! call the models
     if (pft == 1) then
         ! crop pft and we want pft specific model
         call CARBON_MODEL_CROP(1,nodays,met,pars(1:nopars,i),deltat,nodays,lat &
                          ,lai,NEE,FLUXES,POOLS,pft,nopars,nomet,nopools,nofluxes &
                          ,GPP,stock_seed_labile,DS_shoot,DS_root,fol_frac &
                          ,stem_frac,root_frac,DS_LRLV,LRLV,DS_LRRT,LRRT)
     else
         call CARBON_MODEL(1,nodays,met,pars(1:nopars,i),deltat,nodays &
                          ,lat,lai,NEE,FLUXES,POOLS &
                          ,nopars,nomet,nopools,nofluxes,GPP)
     endif
!if (i == 1) then
!    open(unit=666,file="/home/lsmallma/out.csv", &
!         status='replace',action='readwrite' )
!write(666,*)"deltat",deltat
!    write(666,*),"GSI",FLUXES(:,14)(1:365)
!    close(666)
!endif

     ! now allocate the output the our 'output' variable
     out_var(i,1:nodays,1)  = lai
     out_var(i,1:nodays,2)  = GPP
     out_var(i,1:nodays,3)  = FLUXES(1:nodays,3) ! auto resp
     out_var(i,1:nodays,4)  = FLUXES(1:nodays,13) + FLUXES(1:nodays,14) + FLUXES(1:nodays,4)! het resp
     out_var(i,1:nodays,5)  = NEE
     out_var(i,1:nodays,6)  = POOLS(1:nodays,4) ! wood
     out_var(i,1:nodays,7)  = POOLS(1:nodays,6) ! som
     out_var(i,1:nodays,8)  = POOLS(1:nodays,1) + POOLS(1:nodays,2) + POOLS(1:nodays,3) & ! common pools
                            + POOLS(1:nodays,4) !+ POOLS(1:nodays,5) + POOLS(1:nodays,6) + POOLS(1:nodays,7)
     if (pft == 1) out_var(i,1:nodays,8) = out_var(i,1:nodays,8) + POOLS(1:nodays,8) ! crop specific
     out_var(i,1:nodays,9)  = POOLS(1:nodays,3) ! root
     out_var(i,1:nodays,10) = POOLS(1:nodays,5) ! litter
     out_var(i,1:nodays,11) = POOLS(1:nodays,1) ! labile
     out_var(i,1:nodays,12) = POOLS(1:nodays,2) ! foliage
     out_var(i,1:nodays,13) = FLUXES(1:nodays,21) ! harvested material
     ! Phenology related
     if (pft == 1) then
        out_var(i,1:nodays,14) = 0d0
        out_var(i,1:nodays,15) = 0d0 ! GSI temp component
        out_var(i,1:nodays,16) = 0d0 ! GSI photoperiod component
        out_var(i,1:nodays,17) = 0d0 ! GSI vpd component
     else
        out_var(i,1:nodays,14) = FLUXES(1:nodays,18) ! GSI value
        out_var(i,1:nodays,15) = itemp(1:nodays)     ! GSI temp component
        out_var(i,1:nodays,16) = iphoto(1:nodays)    ! GSI photoperiod component
        out_var(i,1:nodays,17) = ivpd(1:nodays)      ! GSI vpd component
     endif
     out_var(i,1:nodays,18) = gs_demand_supply_ratio(1:nodays)
     out_var(i,1:nodays,19) = gs_total_canopy(1:nodays)
     out_var(i,1:nodays,20) = gb_total_canopy(1:nodays)
     if (pft == 1) then
        ! crop so...
        out_var(i,1:nodays,21) = 0d0               ! ...no CWD
        out_var(i,1:nodays,22) = POOLS(1:nodays,7) ! ...Cauto pool present
     else
        ! not a crop...excellent
        out_var(i,1:nodays,21) = POOLS(1:nodays,7) ! ...CWD
        out_var(i,1:nodays,22) = 0d0 ! no Cauto pool present
     endif
     out_var(i,1:nodays,23) = FLUXES(1:nodays,17)    ! output fire (gC/m2/day)
     out_var(i,1:nodays,24) = canopy_par_MJday_time(1:nodays)
     out_var(i,1:nodays,25) = cica_time(1:nodays)

     !!!
     ! NPP calculation
     !!!

     ! calculate the actual NPP allocation fractions to foliar, wood and fine root pools
     ! by comparing the sum alloaction to each pools over the sum NPP.
     fauto = sum(FLUXES(1:nodays,3)) / sum(FLUXES(1:nodays,1))
     sumNPP = (sum(FLUXES(1:nodays,1))*(1d0-fauto))**(-1d0) ! GPP * (1-Ra) fraction
     out_var2(i,1) = sum(FLUXES(1:nodays,8)) * sumNPP ! foliar
     out_var2(i,2) = sum(FLUXES(1:nodays,6)) * sumNPP ! fine root
     out_var2(i,3) = sum(FLUXES(1:nodays,7)) * sumNPP ! wood

     !!!
     ! Estimate residence time information
     !!!

     ! Different calculations are needed depending on whether this is a crop pixel or a generic site
     if (pft == 1) then

         !
         ! Crop model
         !

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
         ! Fol+root litter
         lit_hak = 0 ; lit_filter(1:nodays) = 1d0
         where (POOLS(1:nodays,5) == 0) ! protection against NaN from division by zero
                lit_hak = 1 ; lit_filter(1:nodays) = 0d0
         end where
         ! Soil
         som_hak = 0 ; som_filter(1:nodays) = 1d0
         where (POOLS(1:nodays,6) == 0) ! protection against NaN from division by zero
               som_hak = 1 ; som_filter(1:nodays) = 0d0
         end where

         ! foliage crop system residence time is due to managment < 1 year
         out_var3(i,1) = 1/365.25
         ! roots crop system residence time is due to managment < 1 year
         out_var3(i,2) = 1/365.25
         ! wood crop system residence time is due to managment < 1 year
         out_var3(i,3) = 1/365.25
         ! Lit+litwood
         out_var3(i,4) = sum( ((FLUXES(1:nodays,13)+FLUXES(1:nodays,15)) &
                              / POOLS(1:nodays,5)) * lit_filter) / dble(nodays-sum(lit_hak))
         ! Soil
         out_var3(i,5) = sum( ((FLUXES(1:nodays,14)+fire_loss_som + harvest_loss_som) &
                              / POOLS(1:nodays,6)) * som_filter) / dble(nodays-sum(som_hak))

         ! Assume constant residence times for crops at the moment
         out_var5(i,1,:) = out_var3(i,1)
         out_var5(i,2,:) = out_var3(i,2)
         out_var5(i,3,:) = out_var3(i,3)
         out_var5(i,4,:) = out_var3(i,4)
         out_var5(i,5,:) = out_var3(i,5)

         !
         ! Estimate pool inputs needed for steady state calculation
         !

         out_var4(i,1) = 0d0 ! fol
         out_var4(i,2) = 0d0 ! root
         out_var4(i,3) = 0d0 ! wood
         out_var4(i,4) = 0d0 ! lit
         out_var4(i,5) = sum(FLUXES(:,15)+FLUXES(:,20)+fire_residue_to_som+harvest_residue_to_som) ! som

     else

         !
         ! Generic model
         !

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
         ! Fol+root litter
         lit_hak = 0 ; lit_filter(1:nodays) = 1d0
         where (POOLS(1:nodays,5)+POOLS(1:nodays,7) == 0) ! protection against NaN from division by zero
                lit_hak = 1 ; lit_filter(1:nodays) = 0d0
         end where
         ! Soil
         som_hak = 0 ; som_filter(1:nodays) = 1d0
         where (POOLS(1:nodays,6) == 0) ! protection against NaN from division by zero
               som_hak = 1 ; som_filter(1:nodays) = 0d0
         end where

         ! Estimate MRT (years)
         ! Foliage
         out_var3(i,1) = sum( ((FLUXES(1:nodays,10)+fire_loss_foliar + harvest_loss_foliar) &
                              / POOLS(1:nodays,2)) * fol_filter) / dble(nodays-sum(fol_hak))
         ! Fine roots
         out_var3(i,2) = sum( ((FLUXES(1:nodays,12)+fire_loss_roots + harvest_loss_roots) &
                              / POOLS(1:nodays,3)) * root_filter) / dble(nodays-sum(root_hak))
         ! Wood
         out_var3(i,3) = sum( ((FLUXES(1:nodays,11)+fire_loss_wood + harvest_loss_wood) &
                              / POOLS(1:nodays,4)) * wood_filter) / dble(nodays-sum(wood_hak))
         ! Lit+litwood
         out_var3(i,4) = sum( ((FLUXES(1:nodays,13)+FLUXES(1:nodays,15)+FLUXES(1:nodays,20)+FLUXES(1:nodays,4) &
                               +fire_loss_litter + harvest_loss_litter + fire_loss_litwood + harvest_loss_litwood) &
                              / (POOLS(1:nodays,5)+POOLS(1:nodays,7))) * lit_filter) / dble(nodays-sum(lit_hak))
         ! Soil
         out_var3(i,5) = sum( ((FLUXES(1:nodays,14)+fire_loss_som + harvest_loss_som) &
                              / POOLS(1:nodays,6)) * som_filter) / dble(nodays-sum(som_hak))

         ! Estimate natural MRT (years)
         ! Foliage
         out_var6(i,1) = sum( (FLUXES(1:nodays,10) &
                              / POOLS(1:nodays,2)) * fol_filter) / dble(nodays-sum(fol_hak))
         ! Fine roots
         out_var6(i,2) = sum( (FLUXES(1:nodays,12) &
                              / POOLS(1:nodays,3)) * root_filter) / dble(nodays-sum(root_hak))
         ! Wood
         out_var6(i,3) = sum( (FLUXES(1:nodays,11) &
                              / POOLS(1:nodays,4)) * wood_filter) / dble(nodays-sum(wood_hak))
         ! Lit+litwood
         out_var6(i,4) = sum( ((FLUXES(1:nodays,13)+FLUXES(1:nodays,15)+FLUXES(1:nodays,20)+FLUXES(1:nodays,4)) &
                              / (POOLS(1:nodays,5)+POOLS(1:nodays,7))) * lit_filter) / dble(nodays-sum(lit_hak))
         ! Soil
         out_var6(i,5) = sum( (FLUXES(1:nodays,14) &
                              / POOLS(1:nodays,6)) * som_filter) / dble(nodays-sum(som_hak))

         ! Keep track of the fraction of wood litter transfer to som, this value is needed for the steady state estimation
         litwood_to_som_frac(i) = sum( (FLUXES(1:nodays,20) / POOLS(1:nodays,7)) * lit_filter) / dble(nodays-sum(lit_hak))

         !
         ! Annual residence times
         !

         ! Now loop through each year to estimate the annual residence time
         do y = 1, nos_years
            ! Estimate time steps covered by this year
            y_s = 1 + (steps_per_year * (y-1)) ; y_e = steps_per_year * y
            ! Foliage
            out_var5(i,1,y) = sum( ((FLUXES(y_s:y_e,10) + fire_loss_foliar(y_s:y_e) + harvest_loss_foliar(y_s:y_e)) &
                                   / POOLS(y_s:y_e,2)) * fol_filter(y_s:y_e)) / dble(steps_per_year-sum(fol_hak(y_s:y_e)))
            ! Fine roots
            out_var5(i,2,y) = sum( ((FLUXES(y_s:y_e,12) + fire_loss_roots(y_s:y_e) + harvest_loss_roots(y_s:y_e)) &
                                   / POOLS(y_s:y_e,3)) * root_filter(y_s:y_e)) / dble(steps_per_year-sum(root_hak(y_s:y_e)))
            ! Wood
            out_var5(i,3,y) = sum( ((FLUXES(y_s:y_e,11) + fire_loss_wood(y_s:y_e) + harvest_loss_wood(y_s:y_e)) &
                                   / POOLS(y_s:y_e,4)) * wood_filter(y_s:y_e)) / dble(steps_per_year-sum(wood_hak(y_s:y_e)))
            ! Litter (fol+fine root)
            out_var5(i,4,y) = sum( ((FLUXES(y_s:y_e,13)+FLUXES(y_s:y_e,15)+FLUXES(y_s:y_e,20)+FLUXES(y_s:y_e,4) &
                                    +fire_loss_litter(y_s:y_e)+harvest_loss_litter(y_s:y_e) &
                                    +fire_loss_litwood(y_s:y_e)+harvest_loss_litwood(y_s:y_e)) &
                                   / (POOLS(y_s:y_e,5)+POOLS(y_s:y_e,7))) * lit_filter(y_s:y_e)) &
                            / dble(steps_per_year-sum(lit_hak(y_s:y_e)))
            ! Soil
            out_var5(i,5,y) = sum( ((FLUXES(y_s:y_e,14)+fire_loss_som(y_s:y_e)+harvest_loss_som(y_s:y_e)) &
                                    / POOLS(y_s:y_e,6)) * som_filter(y_s:y_e)) / dble(steps_per_year-sum(som_hak(y_s:y_e)))
         end do

         !
         ! Estimate pool inputs needed for steady state calculation
         !

         ! Once the canopy has closes the inputs to the live biomass are stable
         ! and can thus be estimated from the simulated inputs
         out_var4(i,1) = sum(FLUXES(:,8)) ! Foliage
         out_var4(i,2) = sum(FLUXES(:,6)) ! Fine root
         out_var4(i,3) = sum(FLUXES(:,7)) ! Wood
         ! Foliar and fine root litter can likelise be estimated directly,
         ! however wood litter and soil C inputs are still changing as the wood pool is not in steady state.
         ! Therefore, at this point we can account for the fol, fine root, and disturbance inputs but NOT wood.
         ! The wood input is estimated later based on the steady state its steady state estimate
         out_var4(i,4) = sum(FLUXES(:,10)+FLUXES(:,12)+ &
                             fire_residue_to_litter+fire_residue_to_litwood + &
                             harvest_residue_to_litter+harvest_residue_to_litwood) ! lit + litwood
         out_var4(i,5) = sum(FLUXES(:,15)+fire_residue_to_som+harvest_residue_to_som) ! som
!         out_var4(i,4) = sum(FLUXES(:,10)+FLUXES(:,11)+FLUXES(:,12)+ &
!                             fire_residue_to_litter+fire_residue_to_litwood + &
!                             harvest_residue_to_litter+harvest_residue_to_litwood) ! lit + litwood
!         out_var4(i,5) = sum(FLUXES(:,15)+FLUXES(:,20)+fire_residue_to_som+harvest_residue_to_som) ! som

     endif ! crop / default model split

  end do ! nos_iter loop

  ! MTT - Convert daily fractional loss to years (i.e. day-1 -> yr)
  out_var3 = (out_var3*365.25d0)**(-1d0) ! iter,(fol,root,wood,lit+litwood,som)
  out_var5 = (out_var5*365.25d0)**(-1d0) ! iter,(fol,root,wood,lit+litwood,som)
  out_var6 = (out_var6*365.25d0)**(-1d0) ! iter,(fol,root,wood,lit+litwood,som)

  ! Steady state gC/m2 estimation
  ! Determine the mean annual input (gC/m2/yr) based on current inputs for all pool,
  ! litter and soil pools updated below...
  out_var4 = (out_var4 / dble(nodays)) * 365.25d0 ! convert to annual mean input
  if (pft == 1) then

      ! Multiply by residence time in years to give SS
      out_var4 =  out_var4 * out_var3

  else ! crop / default model split

     ! Then estimate the foliar, fine root and wood steady states.
     out_var4(:,1:3) =  out_var4(:,1:3) * out_var3(:,1:3) ! multiply by residence time in years
     ! Using the wood SS estimate (gC/m2) the steady state input to the litter+wood litter pool...
     out_var4(:,4) = (out_var4(:,4) + (out_var4(:,3) / out_var3(:,3))) * out_var3(:,4)
     ! ...which is then in turn used to update the soil pool
     ! NOTE: that because not all wood litter
     out_var4(:,5) = (out_var4(:,5) + ((out_var4(:,4) / out_var3(:,4))*litwood_to_som_frac) ) * out_var3(:,5)

  end if ! crop / default model split

  ! return back to the subroutine then
  return

end subroutine rdalec
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine crop_development_parameters(stock_seed_labile,DS_shoot,DS_root,fol_frac &
                                        ,stem_frac,root_frac,DS_LRLV,LRLV,DS_LRRT,LRRT &
                                        ,exepath,pathlength)

    ! subroutine reads in the fixed crop development files which are linked the
    ! the development state of the crops. The development model varies between
    ! which species. e.g. winter wheat and barley, spring wheat and barley

    implicit none

    ! declare inputs
    ! crop specific variables
    integer,intent(in) :: pathlength
    character(pathlength),intent(in) :: exepath
    double precision :: stock_seed_labile
    double precision, allocatable, dimension(:) :: DS_shoot, & !
                                                    DS_root, & !
                                                   fol_frac, & !
                                                  stem_frac, & !
                                                  root_frac, & !
                                                    DS_LRLV, & !
                                                       LRLV, & !
                                                    DS_LRRT, & !
                                                       LRRT

    ! local variables..
    integer :: columns, i, rows, input_crops_unit, ios
    character(225) :: variables,filename

    ! file info needed
    input_crops_unit = 20 ; ios = 0

    ! crop development file passed in from the R code (this is different from
    ! *_PARS.f90 where this subroutine is hardcoded)
    open(unit = input_crops_unit, file=trim(exepath),iostat=ios, status='old', action='read')

    ! ensure we are definitely at the beginning
    rewind(input_crops_unit)

    ! read in the amount of carbon available (as labile) in each seed..
    read(unit=input_crops_unit,fmt=*)variables,stock_seed_labile,variables,variables

    ! read in C partitioning/fraction data and corresponding developmental
    ! stages (DS)
    ! shoot
    read(unit=input_crops_unit,fmt=*) variables
    read(unit=input_crops_unit,fmt=*) rows , columns
    allocate( DS_shoot(rows) , fol_frac(rows) , stem_frac(rows)  )
    do i = 1 , rows
      read(unit=input_crops_unit,fmt=*) DS_shoot(i), fol_frac(i), stem_frac(i)
    enddo

    ! root
    read(unit=input_crops_unit,fmt=*) variables
    read(unit=input_crops_unit,fmt=*) rows , columns
    allocate( DS_root(rows) , root_frac(rows) )
    do i = 1 , rows
      read(unit=input_crops_unit,fmt=*) DS_root(i), root_frac(i)
    enddo

    ! loss rates of leaves and roots
    ! leaves
    read(unit=input_crops_unit,fmt=*) variables
    read(unit=input_crops_unit,fmt=*) rows , columns
    allocate( DS_LRLV(rows) , LRLV(rows) )
    do i = 1 , rows
      read(unit=input_crops_unit,fmt=*) DS_LRLV(i), LRLV(i)
    enddo

    ! roots
    read(unit=input_crops_unit,fmt=*) variables
    read(unit=input_crops_unit,fmt=*) rows , columns
    allocate( DS_LRRT(rows) , LRRT(rows) )
    do i = 1 , rows
      read(unit=input_crops_unit,fmt=*) DS_LRRT(i), LRRT(i)
    enddo

    ! rewind and close
    rewind(input_crops_unit) ; close(input_crops_unit)

  end subroutine crop_development_parameters
  !
  !------------------------------------------------------------------
  !
