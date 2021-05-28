

subroutine rdalecgsibiofr(output_dim,aNPP_dim,met,pars,out_var,out_var2,lat &
                         ,nopars,nomet,nofluxes,nopools,pft,pft_specific &
                         ,nodays,deltat,nos_iter)

  use CARBON_MODEL_MOD, only: CARBON_MODEL, microbial_activity_out,    &
                              itemp,iphoto,ivpd

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
                        ,aNPP_dim       & ! NPP allocation fraction variable dimension
                        ,pft            & ! plant functional type
                        ,pft_specific   & !
                        ,nos_iter       & !
                        ,nomet          & ! number of meteorological fields
                        ,nofluxes       & ! number of model fluxes
                        ,nopools        & ! number of model pools
                        ,nodays           ! number of days in simulation

  double precision, intent(inout) :: deltat(nodays) ! time step in decimal days
  double precision, intent(in) :: met(nomet,nodays)   & ! met drivers, note reverse of needed
                       ,pars(nopars,nos_iter)         & ! number of parameters
                       ,lat                 ! site latitude (degrees)

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
  double precision :: fauto, sumNPP, weighting_fast_slow, tmp1
  double precision, dimension(nodays) :: lai & ! leaf area index
                                        ,GPP & ! Gross primary productivity
                                        ,NEE   ! net ecosystem exchange of CO2

  ! zero initial conditions
  lai = 0d0 ; GPP = 0d0 ; NEE = 0d0 ; POOLS = 0d0 ; FLUXES = 0d0 ; out_var = 0d0

  ! update settings
  if (allocated(itemp)) deallocate(itemp,ivpd,iphoto)
  allocate(itemp(nodays),ivpd(nodays),iphoto(nodays))
  if (allocated(microbial_activity_out)) deallocate(microbial_activity_out)
  allocate(microbial_activity_out(nodays+1))

  ! generate deltat step from input data
  deltat(1) = met(1,1)
  do i = 2, nodays
     deltat(i)=met(1,i)-met(1,(i-1))
  end do

  ! calculate inverse of nodays
  tmp1 = 1/dble(nodays)

  ! begin iterations
  do i = 1, nos_iter

     ! call the models
     call CARBON_MODEL(1,nodays,met,pars(1:nopars,i),deltat,nodays &
                      ,lat,lai,NEE,FLUXES,POOLS &
                      ,nopars,nomet,nopools,nofluxes,GPP)

     ! now allocate the output the our 'output' variable
     out_var(i,1:nodays,1)  = lai
     out_var(i,1:nodays,2)  = GPP
     out_var(i,1:nodays,3)  = FLUXES(1:nodays,3) ! auto resp
     out_var(i,1:nodays,4)  = FLUXES(1:nodays,13) + FLUXES(1:nodays,14) + &
                              FLUXES(1:nodays,15) + FLUXES(1:nodays,19) + &
                              FLUXES(1:nodays,20) + FLUXES(1:nodays,21) ! het resp
     out_var(i,1:nodays,5)  = NEE
     out_var(i,1:nodays,6)  = POOLS(1:nodays,4) ! wood
     out_var(i,1:nodays,7)  = POOLS(1:nodays,6) ! som slow
     out_var(i,1:nodays,8)  = POOLS(1:nodays,1) + POOLS(1:nodays,2) + POOLS(1:nodays,3) & ! plant biomass only
                              + POOLS(1:nodays,4) !+ POOLS(1:nodays,5) + POOLS(1:nodays,6) &
!                              + POOLS(1:nodays,7) + POOLS(1:nodays,8) + POOLS(1:nodays,9) &
!                              + POOLS(1:nodays,10)
     out_var(i,1:nodays,9)  = POOLS(1:nodays,3) ! root
     out_var(i,1:nodays,10) = POOLS(1:nodays,7) ! foliar litter
     out_var(i,1:nodays,11) = POOLS(1:nodays,1) ! labile
     out_var(i,1:nodays,12) = POOLS(1:nodays,2) ! foliage
     out_var(i,1:nodays,13) = FLUXES(1:nodays,28) ! harvested material
     out_var(i,1:nodays,14) = FLUXES(1:nodays,18) ! GSI value
     out_var(i,1:nodays,15) = itemp(1:nodays) ! GSI temp component
     out_var(i,1:nodays,16) = iphoto(1:nodays) ! GSI photoperiod component
     out_var(i,1:nodays,17) = ivpd(1:nodays) ! GSI vpd component
     out_var(i,1:nodays,18) = POOLS(1:nodays,5) ! som fast
     out_var(i,1:nodays,19) = POOLS(1:nodays,8) ! root litter
     out_var(i,1:nodays,20) = POOLS(1:nodays,9) ! wood litter
     out_var(i,1:nodays,21) = microbial_activity_out(1:nodays)
     out_var(i,1:nodays,22) = POOLS(1:nodays,10) ! microbial pool

     !!!
     ! NPP calculation
     !!!

     ! calculate the actual NPP allocation fractions to foliar, wood and fine root pools
     ! by comparing the sum alloaction to each pools over the sum NPP.
     fauto = sum(FLUXES(1:nodays,3)) / sum(FLUXES(1:nodays,1))
     sumNPP = 1 / (sum(FLUXES(1:nodays,1))*(1d0-fauto))
!     out_var2(i,1) = sum(FLUXES(1:nodays,4)+FLUXES(1:nodays,8)) * sumNPP ! foliar
     out_var2(i,1) = sum(FLUXES(1:nodays,5)) * sumNPP ! foliar+labile
     out_var2(i,2) = sum(FLUXES(1:nodays,7)) * sumNPP ! wood
     out_var2(i,3) = sum(FLUXES(1:nodays,6)) * sumNPP ! fine root

     !!!
     ! now some residence times (years)
     !!!

     !! foliage !!
     resid_fol(1:nodays) = FLUXES(1:nodays,10)
     resid_fol(1:nodays) = resid_fol(1:nodays) &
                         / POOLS(1:nodays,2)
     ! division by zero results in NaN plus obviously I can't have turned
     ! anything over if there was nothing to start out with...
     hak = 0
     where ( POOLS(1:nodays,2) == 0 )
             hak = 1 ; resid_fol(1:nodays) = 0d0
     end where
     out_var2(i,4) = sum(resid_fol) /dble(nodays-sum(hak))

     !! wood !!
     resid_fol(1:nodays)   = FLUXES(1:nodays,11)
     resid_fol(1:nodays)   = resid_fol(1:nodays) &
                           / POOLS(1:nodays,4)
     ! division by zero results in NaN plus obviously I can't have turned
     ! anything over if there was nothing to start out with...
     hak = 0
     where ( POOLS(1:nodays,4) == 0 )
             hak = 1 ; resid_fol(1:nodays) = 0d0
     end where
     out_var2(i,5) = sum(resid_fol) /dble(nodays-sum(hak))

     !! roots !!
     resid_fol(1:nodays)   = FLUXES(1:nodays,12)
     resid_fol(1:nodays)   = resid_fol(1:nodays) &
                           / POOLS(1:nodays,3)
     ! division by zero results in NaN plus obviously I can't have turned
     ! anything over if there was nothing to start out with...
     hak = 0
     where ( POOLS(1:nodays,3) == 0 )
             hak = 1 ; resid_fol(1:nodays) = 0d0
     end where
     out_var2(i,6) = sum(resid_fol) /dble(nodays-sum(hak))

     ! Csom - include fast and slow only not microbial (microbial commented out)
!     out_var2(i,7) = sum((FLUXES(1:nodays,19)+FLUXES(1:nodays,20)+FLUXES(1:nodays,21)) / &
!                         (POOLS(1:nodays,5)+POOLS(1:nodays,6)+POOLS(1:nodays,10)))       &
!                         * tmp1 !/ dble(nodays)
     out_var2(i,7) = sum((FLUXES(1:nodays,19)+FLUXES(1:nodays,20)) / &
                         (POOLS(1:nodays,5)+POOLS(1:nodays,6))) * tmp1 !/ dble(nodays)

     !! flitter+rlitter+wlitter !!
     resid_fol(1:nodays)   = FLUXES(1:nodays,13)+FLUXES(1:nodays,14)+FLUXES(1:nodays,15) &
                           + FLUXES(1:nodays,22)+FLUXES(1:nodays,23)+FLUXES(1:nodays,24)
     resid_fol(1:nodays)   = resid_fol(1:nodays) &
                           / (POOLS(1:nodays,7)+POOLS(1:nodays,8)+POOLS(1:nodays,9))
     out_var2(i,8) = sum(resid_fol) / dble(nodays)

  end do ! nos_iter loop

  ! moving this out of the loop to calculate fractions to years residence times
  out_var2(1:nos_iter,4) = 1 / (out_var2(1:nos_iter,4)*365.25) ! fol
  out_var2(1:nos_iter,5) = 1 / (out_var2(1:nos_iter,5)*365.25) ! wood
  out_var2(1:nos_iter,6) = 1 / (out_var2(1:nos_iter,6)*365.25) ! root
  out_var2(1:nos_iter,7) = 1 / (out_var2(1:nos_iter,7)*365.25) ! som
  out_var2(1:nos_iter,8) = 1 / (out_var2(1:nos_iter,8)*365.25) ! CWD + Litter

  ! deallocate harvested variable
  deallocate(itemp,ivpd,iphoto,microbial_activity_out)

  ! return back to the subroutine then
  return

end subroutine rdalecgsibiofr
