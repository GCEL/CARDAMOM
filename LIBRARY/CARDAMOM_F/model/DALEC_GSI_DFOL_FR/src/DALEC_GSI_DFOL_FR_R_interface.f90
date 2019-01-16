

subroutine rdalecgsidfolfr(output_dim,aNPP_dim,met,pars,out_var,out_var2,lat & 
                          ,nopars,nomet,nofluxes,nopools,pft,pft_specific & 
                          ,nodays,deltat,nos_iter)

  use CARBON_MODEL_MOD, only: CARBON_MODEL, extracted_C, itemp, ivpd, iphoto
  use CARBON_MODEL_CROP_MOD

  ! subroutine specificially deals with the calling of the fortran code model by
  ! R

  implicit none
  ! declare input variables
  integer, intent(in) :: nopars         & ! number of paremeters in vector
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
  double precision :: sumNPP
  double precision, dimension(nodays) :: lai & ! leaf area index
                                        ,GPP & ! Gross primary productivity
                                        ,NEE   ! net ecosystem exchange of CO2

  ! zero initial conditions
  lai = 0d0 ; GPP = 0d0 ; NEE = 0d0 ; POOLS = 0d0 ; FLUXES = 0d0 ; out_var = 0d0

  ! update settings
  if (allocated(extracted_C)) deallocate(extracted_C,itemp,ivpd,iphoto)
  allocate(extracted_C(nodays+1),itemp(nodays),ivpd(nodays),iphoto(nodays))

  ! generate deltat step from input data
  deltat(1) = met(1,1)
  do i = 2, nodays
     deltat(i)=met(1,i)-met(1,(i-1))
  end do

  ! when crop model in use should load crop development parameters here
  ! modifications neede....
!  if (pft == 1) call crop_development_parameters()

  ! begin iterations
  do i = 1, nos_iter
     ! reset harvest variable
     extracted_C=0d0
     ! call the models
     call CARBON_MODEL(1,nodays,met,pars(1:nopars,i),deltat,nodays &
                      ,lat,lai,NEE,FLUXES,POOLS &
                      ,nopars,nomet,nopools,nofluxes,GPP)
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
                              + POOLS(1:nodays,4) + POOLS(1:nodays,5) + POOLS(1:nodays,6)
     out_var(i,1:nodays,9)  = POOLS(1:nodays,3) ! root
     out_var(i,1:nodays,10) = POOLS(1:nodays,5) ! litter
     out_var(i,1:nodays,11) = POOLS(1:nodays,1) ! labile
     out_var(i,1:nodays,12) = POOLS(1:nodays,2) ! foliage
     out_var(i,1:nodays,13) = extracted_C(1:nodays) ! harvested material
     out_var(i,1:nodays,14) = FLUXES(1:nodays,18) ! GSI value
     out_var(i,1:nodays,15) = itemp(1:nodays) ! GSI temp component
     out_var(i,1:nodays,16) = iphoto(1:nodays) ! GSI photoperiod component
     out_var(i,1:nodays,17) = ivpd(1:nodays) ! GSI vpd component

     ! calculate the actual NPP allocation fractions to foliar, wood and fine root pools 
     ! by comparing the sum alloaction to each pools over the sum NPP.
     sumNPP = 1 / sum(FLUXES(1:nodays,1)*(1d0-pars(2,i))) ! GPP * (1-Ra) fraction
     out_var2(i,1) = sum(FLUXES(1:nodays,4)+FLUXES(1:nodays,8)) * sumNPP ! foliar
     out_var2(i,2) = sum(FLUXES(1:nodays,7)) * sumNPP ! wood
     out_var2(i,3) = sum(FLUXES(1:nodays,6)) * sumNPP ! fine root
     ! now some residence times (years)
     hak = 0
     resid_fol(1:nodays) = FLUXES(1:nodays,10)/POOLS(1:nodays,2)
     ! division by zero results in NaN plus obviously I can't have turned
     ! anything over if there was nothing to start out with...
     where ( POOLS(1:nodays,2) == 0 ) 
            hak = 1 ; resid_fol(1:nodays) = 0d0
     end where
     out_var2(i,4) = sum(resid_fol) /dble(nodays-sum(hak))

     ! Csom 
     resid_fol(1:nodays) = FLUXES(1:nodays,14)/POOLS(1:nodays,6)
     out_var2(i,7) = sum(resid_fol) /dble(nodays)

  end do ! nos_iter loop

  ! moving this out of the loop to calculate fractions to years residence times
  out_var2(1:nos_iter,4) = 1 / (out_var2(1:nos_iter,4)*365.25) ! fol
  out_var2(1:nos_iter,5) = 1 / (pars(6,1:nos_iter)*365.25) ! wood
  out_var2(1:nos_iter,6) = 1 / (pars(7,1:nos_iter)*365.25) ! root 
  out_var2(1:nos_iter,7) = 1 / (out_var2(1:nos_iter,7)*365.25) ! som

  ! deallocate harvested variable
  deallocate(extracted_C,itemp,ivpd,iphoto)

  ! return back to the subroutine then
  return

end subroutine rdalecgsidfolfr
