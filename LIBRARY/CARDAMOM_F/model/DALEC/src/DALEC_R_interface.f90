

subroutine rdalec(output_dim,aNPP_dim,MTT_dim,SS_dim,met,pars,out_var,out_var2,out_var3,out_var4 & 
                 ,lat &
                 ,nopars,nomet,nofluxes,nopools,pft,pft_specific &
                 ,nodays,deltat,nos_iter,exepath,pathlength)

  use CARBON_MODEL_MOD, only: CARBON_MODEL, itemp, ivpd, iphoto &
                             ,disturbance_residue_to_litter, disturbance_residue_to_cwd &
                             ,disturbance_residue_to_som, disturbance_loss_from_litter  &
                             ,disturbance_loss_from_cwd,disturbance_loss_from_som
  use CARBON_MODEL_CROP_MOD, only: CARBON_MODEL_CROP

  ! subroutine specificially deals with the calling of the fortran code model by
  ! R

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
  double precision, intent(out), dimension(nos_iter,aNPP_dim) :: out_var2
  double precision, intent(out), dimension(nos_iter,MTT_dim) :: out_var3
  double precision, intent(out), dimension(nos_iter,SS_dim) :: out_var4

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

! profiling example
!real :: begin, done,f1=0,f2=0,f3=0,f4=0,f5=0,total_time = 0
!real :: Rtot_track_time = 0, aero_time = 0 , soilwater_time = 0 , acm_et_time = 0 , Rm_time = 0
!call cpu_time(done)
!print*,"time taken per iter",(done-begin) / real(nos_iter)

  ! zero initial conditions
  lai = 0.0 ; GPP = 0.0 ; NEE = 0.0 ; POOLS = 0.0 ; FLUXES = 0.0 ; out_var = 0.0

  ! generate deltat step from input data
  deltat(1) = met(1,1)
  do i = 2, nodays
     deltat(i) = met(1,i)-met(1,(i-1))
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
     out_var(i,1:nodays,18) = 0d0
     out_var(i,1:nodays,19) = 0d0   
     out_var(i,1:nodays,20) = 0d0
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
     out_var(i,1:nodays,24) = 0d0

     ! calculate the actual NPP allocation fractions to foliar, wood and fine root pools
     ! by comparing the sum alloaction to each pools over the sum NPP.
     fauto = sum(FLUXES(1:nodays,3)) / sum(FLUXES(1:nodays,1))
     sumNPP = (sum(FLUXES(1:nodays,1))*(1d0-fauto))**(-1d0) ! GPP * (1-Ra) fraction
     out_var2(i,1) = sum(FLUXES(1:nodays,8)) * sumNPP ! foliar
     out_var2(i,2) = sum(FLUXES(1:nodays,6)) * sumNPP ! fine root
     out_var2(i,3) = sum(FLUXES(1:nodays,7)) * sumNPP ! wood

     ! Estimate residence times (years) and begin calculation of steady state attractor
     hak = 0
     if (pft == 1) then

         !
         ! Residence time
         !

         ! foliage crop system residence time is due to managment < 1 year
         out_var3(i,1) = 1/365.25
         ! roots crop system residence time is due to managment < 1 year
         out_var3(i,2) = 1/365.25
         ! wood crop system residence time is due to managment < 1 year
         out_var3(i,3) = 1/365.25
         ! cwd+litter / litter
         resid_fol(1:nodays)   = (FLUXES(1:nodays,13)+FLUXES(1:nodays,15))
         resid_fol(1:nodays)   = resid_fol(1:nodays) &
                               / POOLS(1:nodays,5)
         ! division by zero results in NaN plus obviously I can't have turned
         ! anything over if there was nothing to start out with...
         where ( POOLS(1:nodays,5) == 0 )
                hak = 1 ; resid_fol(1:nodays) = 0d0
         end where
         out_var3(i,4) = sum(resid_fol) / dble(nodays)

         ! 
         ! Estimate pool inputs needed for steady state calculation
         !

         out_var4(i,1) = 0d0 ! fol
         out_var4(i,2) = 0d0 ! root
         out_var4(i,3) = 0d0 ! wood
         out_var4(i,4) = 0d0 ! lit

     else

         !
         ! Residence times
         !

         hak = 0
         ! foliage
         resid_fol(1:nodays) = (FLUXES(1:nodays,10) + FLUXES(1:nodays,23)) / POOLS(1:nodays,2)
         ! division by zero results in NaN plus obviously I can't have turned
         ! anything over if there was nothing to start out with...
         where ( POOLS(1:nodays,2) == 0 )
                hak = 1 ; resid_fol(1:nodays) = 0d0
         end where
         out_var3(i,1) = sum(resid_fol) / (dble(nodays)-dble(sum(hak)))

         ! roots
         hak = 0
         resid_fol(1:nodays)   = FLUXES(1:nodays,12)+FLUXES(1:nodays,24)
         resid_fol(1:nodays)   = resid_fol(1:nodays) &
                               / POOLS(1:nodays,3)
         ! division by zero results in NaN plus obviously I can't have turned
         ! anything over if there was nothing to start out with...
         where ( POOLS(1:nodays,3) == 0 )
                hak = 1 ; resid_fol(1:nodays) = 0d0
         end where
         out_var3(i,2) = sum(resid_fol) /dble(nodays-sum(hak))

         ! wood
         hak = 0
         resid_fol(1:nodays)   = FLUXES(1:nodays,11)+FLUXES(1:nodays,25)
         resid_fol(1:nodays)   = resid_fol(1:nodays) &
                               / POOLS(1:nodays,4)
         ! division by zero results in NaN plus obviously I can't have turned
         ! anything over if there was nothing to start out with...
         where ( POOLS(1:nodays,4) == 0 )
                hak = 1 ; resid_fol(1:nodays) = 0d0
         end where
         out_var3(i,3) = sum(resid_fol) /dble(nodays-sum(hak))

         ! litter + cwd
         resid_fol(1:nodays)   = FLUXES(1:nodays,13)+FLUXES(1:nodays,15) &
                                +FLUXES(1:nodays,20)+FLUXES(1:nodays,4)  &
                                 +disturbance_Loss_from_litter+disturbance_loss_from_cwd
         resid_fol(1:nodays)   = resid_fol(1:nodays) &
                               / (POOLS(1:nodays,5)+POOLS(1:nodays,7))
         out_var3(i,4) = sum(resid_fol) / dble(nodays)


         ! 
         ! Estimate pool inputs needed for steady state calculation
         !

         out_var4(i,1) = sum(FLUXES(:,8)) ! Foliage
         out_var4(i,2) = sum(FLUXES(:,6)) ! Fine root
         out_var4(i,3) = sum(FLUXES(:,7)) ! Wood
         out_var4(i,4) = sum(FLUXES(:,10)+FLUXES(:,11)+FLUXES(:,12)+ &
                             disturbance_residue_to_litter+disturbance_residue_to_cwd) ! lit + litwood

     endif ! crop choice

     ! Csom - residence time
     resid_fol(1:nodays)   = FLUXES(1:nodays,14) + disturbance_loss_from_som
     resid_fol(1:nodays)   = resid_fol(1:nodays) &
                           / POOLS(1:nodays,6)
     out_var3(i,5) = sum(resid_fol) /dble(nodays)

     ! Csom - pool inputs needed for steady state calculation
     out_var4(i,5) = sum(FLUXES(:,15)+FLUXES(:,20)+disturbance_residue_to_som) ! som

  end do ! nos_iter loop

  ! MTT - Convert daily fractional loss to years 
  out_var3 = (out_var3*365.25d0)**(-1d0) ! iter,(fol,root,wood,lit+litwood,som)
!  out_var3(1:nos_iter,1) = (out_var3(1:nos_iter,1)*365.25d0)**(-1d0) ! fol
!  out_var3(1:nos_iter,2) = (out_var3(1:nos_iter,2)*365.25d0)**(-1d0) ! root
!  out_var3(1:nos_iter,3) = (out_var3(1:nos_iter,3)*365.25d0)**(-1d0) ! wood
!  out_var3(1:nos_iter,4) = (out_var3(1:nos_iter,4)*365.25d0)**(-1d0) ! CWD + Litter
!  out_var3(1:nos_iter,5) = (out_var3(1:nos_iter,5)*365.25d0)**(-1d0) ! som

  ! Steady state gC/m2
  out_var4 = (out_var4 / dble(nodays)) * 365.25d0 ! convert to daily mean input
  out_var4 = out_var4 * out_var3     ! multiply by residence time in years

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