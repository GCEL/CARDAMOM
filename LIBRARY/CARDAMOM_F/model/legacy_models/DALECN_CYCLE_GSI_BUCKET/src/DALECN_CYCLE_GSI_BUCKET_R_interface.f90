

subroutine rdalecncyclegsibucket(output_dim,aNPP_dim,met,pars,out_var,out_var2,lat &
                          ,nopars,nomet,nofluxes,nopools,pft,pft_specific &
                          ,nodays,deltat,nos_iter,soil_frac_clay_in,soil_frac_sand_in &
                          ,exepath,pathlength)

  use CARBON_MODEL_MOD, only: CARBON_MODEL, extracted_C, itemp, ivpd, iphoto &
                             ,soil_frac_clay, soil_frac_sand &
                             ,nos_soil_layers, Rtot_time
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
      character(pathlength),intent(in) :: exepath
      integer, intent(in) :: pathlength
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
  character(pathlength), intent(in) :: exepath
  integer, intent(in) :: nopars         & ! number of paremeters in vector
                        ,pathlength     & !
                        ,output_dim     & !
                        ,aNPP_dim       & ! NPP allocation fraction variable dimension
                        ,pft            & ! plant functional type
                        ,pft_specific   & !
                        ,nos_iter       & !
                        ,nomet          & ! number of meteorological fields
                        ,nofluxes       & ! number of model fluxes
                        ,nopools        & ! number of model pools
                        ,nodays           ! number of days in simulation

  double precision, intent(in) :: met(nomet,nodays)   & ! met drivers, note reverse of needed
                       ,pars(nopars,nos_iter)         & ! number of parameters
                       ,soil_frac_clay_in(nos_soil_layers) & ! clay in soil (%)
                       ,soil_frac_sand_in(nos_soil_layers) & ! sand in soil (%)
                       ,lat                 ! site latitude (degrees)

  double precision, intent(inout) :: deltat(nodays) ! time step in decimal days

  ! output declaration
  double precision, intent(out), dimension(nos_iter,nodays,output_dim) :: out_var
  double precision, intent(out), dimension(nos_iter,aNPP_dim) :: out_var2

  ! local variables
  integer i
  ! vector of ecosystem pools
  double precision, dimension((nodays+1),nopools) :: POOLS
  ! vector of ecosystem fluxes
  double precision, dimension(nodays,nofluxes) :: FLUXES
  double precision, dimension(nodays) :: resid_fol
  integer, dimension(nodays) :: hak ! variable to determine number of NaN
  double precision :: sumNPP, fauto
  double precision, dimension(nodays) :: lai & ! leaf area index
                                        ,GPP & ! Gross primary productivity
                                        ,NEE   ! net ecosystem exchange of CO2

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

  ! zero initial conditions
  lai = 0d0 ; GPP = 0d0 ; NEE = 0d0 ; POOLS = 0d0 ; FLUXES = 0d0 ; out_var = 0d0

  ! update settings
  if (allocated(itemp)) deallocate(itemp,ivpd,iphoto)
  allocate(itemp(nodays),ivpd(nodays),iphoto(nodays))
  ! update soil parameters
  soil_frac_clay=soil_frac_clay_in
  soil_frac_sand=soil_frac_clay_in

  ! generate deltat step from input data
  deltat(1) = met(1,1)
  do i = 2, nodays
     deltat(i)=met(1,i)-met(1,(i-1))
  end do

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
     out_var(i,1:nodays,4)  = FLUXES(1:nodays,13) + FLUXES(1:nodays,14) ! het resp
     out_var(i,1:nodays,5)  = NEE
     out_var(i,1:nodays,6)  = POOLS(1:nodays,4) ! wood
     out_var(i,1:nodays,7)  = POOLS(1:nodays,6) ! som
     out_var(i,1:nodays,8)  = POOLS(1:nodays,1) + POOLS(1:nodays,2) + POOLS(1:nodays,3) & ! common pools
                              + POOLS(1:nodays,4) + POOLS(1:nodays,5) + POOLS(1:nodays,6) + POOLS(1:nodays,7)
     if (pft == 1) out_var(i,1:nodays,8) = out_var(i,1:nodays,8) + POOLS(1:nodays,9) ! crop specific
     out_var(i,1:nodays,9)  = POOLS(1:nodays,3) ! root
     out_var(i,1:nodays,10) = POOLS(1:nodays,5) ! litter
     out_var(i,1:nodays,11) = POOLS(1:nodays,1) ! labile
     out_var(i,1:nodays,12) = POOLS(1:nodays,2) ! foliage
     if (pft /= 1) then
         out_var(i,1:nodays,13) = FLUXES(1:nodays,21) ! harvested material
     else
         out_var(i,1:nodays,13) = FLUXES(1:nodays,21) !POOLS(1:nodays,8)! replace with crop model yield
     endif
     if (pft == 1) then
        out_var(i,1:nodays,14) = 0d0
        out_var(i,1:nodays,15) = 0d0 ! GSI temp component
        out_var(i,1:nodays,16) = 0d0 ! GSI photoperiod component
        out_var(i,1:nodays,17) = 0d0 ! GSI vpd component
     else
        out_var(i,1:nodays,14) = FLUXES(1:nodays,18) ! GSI value
        out_var(i,1:nodays,15) = itemp(1:nodays) ! GSI temp component
        out_var(i,1:nodays,16) = iphoto(1:nodays) ! GSI photoperiod component
        out_var(i,1:nodays,17) = ivpd(1:nodays) ! GSI vpd component
     endif
     out_var(i,1:nodays,18) = FLUXES(1:nodays,19) ! Evapotranspiration (kgH2O.m-2.day-1)
     out_var(i,1:nodays,19) = POOLS(1:nodays,8) ! rootwater (kgH2O.m-2.rootdepth)
     out_var(i,1:nodays,20) = Rtot_time(1:nodays) ! MPa total hydraulic resistance
     if (pft == 1) then
        ! crop so...
        out_var(i,1:nodays,21) = 0d0               ! ...no CWD
        out_var(i,1:nodays,22) = POOLS(1:nodays,7) ! ...Cauto pool present
     else
        ! not a crop...excellent
        out_var(i,1:nodays,21) = POOLS(1:nodays,7) ! ...CWD
        out_var(i,1:nodays,22) = 0d0 ! no Cauto pool present
     endif
     ! output fire
     out_var(i,1:nodays,23) = FLUXES(1:nodays,17)

     ! calculate the actual NPP allocation fractions to foliar, wood and fine root pools
     ! by comparing the sum alloaction to each pools over the sum NPP.
     fauto = sum(FLUXES(1:nodays,3)) / sum(FLUXES(1:nodays,1))
     sumNPP = (sum(FLUXES(1:nodays,1))*(1d0-fauto))**(-1d0) ! GPP * (1-Ra) fraction
     out_var2(i,1) = sum(FLUXES(1:nodays,8)) * sumNPP ! foliar
     out_var2(i,2) = sum(FLUXES(1:nodays,7)) * sumNPP ! wood
     out_var2(i,3) = sum(FLUXES(1:nodays,6)) * sumNPP ! fine root

     ! now some residence times (years)
     if (pft == 1) then
         ! foliage crop system residence time is due to managment < 1 year
         out_var2(i,4) = 1/365.25
         ! wood crop system residence time is due to managment < 1 year
         out_var2(i,5) = 1/365.25
         ! roots crop system residence time is due to managment < 1 year
         out_var2(i,6) = 1/365.25
         ! cwd+litter / litter
         resid_fol(1:nodays)   = (FLUXES(1:nodays,13)+FLUXES(1:nodays,15))
         resid_fol(1:nodays)   = resid_fol(1:nodays) &
                               / POOLS(1:nodays,5)
         ! division by zero results in NaN plus obviously I can't have turned
         ! anything over if there was nothing to start out with...
         hak = 0
         where ( POOLS(1:nodays,5) == 0 )
                hak = 1 ; resid_fol(1:nodays) = 0d0
         end where
         out_var2(i,8) = sum(resid_fol) / dble(nodays)
     else

         ! foliage
         resid_fol(1:nodays) = FLUXES(1:nodays,10) + FLUXES(1:nodays,23)
         resid_fol(1:nodays) = resid_fol(1:nodays) &
                             / POOLS(1:nodays,2)
         ! division by zero results in NaN plus obviously I can't have turned
         ! anything over if there was nothing to start out with...
         hak = 0
         where ( POOLS(1:nodays,2) == 0d0 )
                hak = 1 ; resid_fol(1:nodays) = 0d0
         end where
         out_var2(i,4) = sum(resid_fol) /dble(nodays-sum(hak))

         ! wood
         resid_fol(1:nodays)   = FLUXES(1:nodays,11) + FLUXES(1:nodays,25)
         resid_fol(1:nodays)   = resid_fol(1:nodays) &
                               / POOLS(1:nodays,4)
         ! division by zero results in NaN plus obviously I can't have turned
         ! anything over if there was nothing to start out with...
         hak = 0
         where ( POOLS(1:nodays,4) == 0d0 )
                hak = 1 ; resid_fol(1:nodays) = 0d0
         end where
         out_var2(i,5) = sum(resid_fol) /dble(nodays-sum(hak))
         ! roots
         resid_fol(1:nodays)   = FLUXES(1:nodays,12) + FLUXES(1:nodays,24)
         resid_fol(1:nodays)   = resid_fol(1:nodays) &
                               / POOLS(1:nodays,3)
         ! division by zero results in NaN plus obviously I can't have turned
         ! anything over if there was nothing to start out with...
         hak = 0
         where ( POOLS(1:nodays,3) == 0d0 )
                hak = 1 ; resid_fol(1:nodays) = 0d0
         end where
         out_var2(i,6) = sum(resid_fol) /dble(nodays-sum(hak))
         ! cwd
         resid_fol(1:nodays)   = (FLUXES(1:nodays,13)+FLUXES(1:nodays,15))
         resid_fol(1:nodays)   = resid_fol(1:nodays) &
                               / (POOLS(1:nodays,5)+POOLS(1:nodays,7))
         out_var2(i,8) = sum(resid_fol) / dble(nodays)
     endif ! crop choice

     ! Csom
     resid_fol(1:nodays)   = FLUXES(1:nodays,14)
     resid_fol(1:nodays)   = resid_fol(1:nodays) &
                           / POOLS(1:nodays,6)
     out_var2(i,7) = sum(resid_fol) /dble(nodays)

  end do ! nos_iter loop

  ! moving this out of the loop to calculate fractions to years residence times
  out_var2(1:nos_iter,4) = (out_var2(1:nos_iter,4)*365.25)**(-1d0) ! fol
  out_var2(1:nos_iter,5) = (out_var2(1:nos_iter,5)*365.25)**(-1d0) ! wood
  out_var2(1:nos_iter,6) = (out_var2(1:nos_iter,6)*365.25)**(-1d0) ! root
  out_var2(1:nos_iter,7) = (out_var2(1:nos_iter,7)*365.25)**(-1d0) ! som
  out_var2(1:nos_iter,8) = (out_var2(1:nos_iter,8)*365.25)**(-1d0) ! CWD + Litter

  ! deallocate harvested variable
  deallocate(itemp,ivpd,iphoto)

  ! return back to the subroutine then
  return

  end subroutine rdalecncyclegsibucket
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
    character(pathlength),intent(in) :: exepath
    integer,intent(in) :: pathlength
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
