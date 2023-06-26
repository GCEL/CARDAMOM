
subroutine ratdalec(output_dim,aNPP_dim,met,pars,out_var,out_var2,lat,nopars,nomet &
                   ,nofluxes,nopools,pft,pft_specific,nodays,deltat &
                   ,dim1,dim2,cdim,dummy_nos_trees,dummy_nos_inputs,nos_iter &
                   ,stock_seed_labile,tmp1,tmp2,tmp3,tmp4 &
                   ,tmp5,tmp6,DS_shoot,fol_frac,stem_frac,DS_root &
                   ,root_frac,DS_LRLV,LRLV,DS_LRRT,LRRT)

  use CARBON_MODEL_MOD, only: CARBON_MODEL &
                             ,nos_trees,nos_inputs    &
                             ,leftDaughter,rightDaughter &
                             ,nodestatus,xbestsplit    &
                             ,nodepred,bestvar
  use CARBON_MODEL_CROP_MOD

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
  integer, intent(in) :: nopars       & ! number of paremeters in vector
                        ,output_dim   & !
                        ,aNPP_dim     & ! NPP allocation fraction variable dimension
                        ,dim1         & ! dimension 1 of response surface
                        ,dim2         & ! dimension 2 of response surface
                        ,cdim(9)      & ! dimension of crop development vars
                        ,dummy_nos_trees    & ! number of trees in randomForest
                        ,dummy_nos_inputs   & ! number of driver input
                        ,pft          & ! plant functional type
                        ,pft_specific & !
                        ,nos_iter     & !
                        ,nomet        & ! number of meteorological fields
                        ,nofluxes     & ! number of model fluxes
                        ,nopools      & ! number of model pools
                        ,nodays         ! number of days in simulation

  double precision, intent(in) :: met(nomet,nodays)   & ! met drivers, note reverse of needed
                       ,stock_seed_labile             & ! seed carbon to get things going
                       ,pars(nopars,nos_iter)         & ! number of parameters
                       ,lat                 ! site latitude (degrees)

  double precision, intent(inout) :: deltat(nodays) ! time step in decimal days

  ! random forest related vars
  double precision, intent(in), dimension(dim1,dim2) :: tmp1,tmp2,tmp3 &
                                                       ,tmp4,tmp5,tmp6
  ! output declaration
  double precision, intent(out), dimension(nos_iter,nodays,output_dim) :: out_var
  double precision, intent(out), dimension(nos_iter,aNPP_dim) :: out_var2

  ! crop related vars
  double precision :: DS_shoot(cdim(1)),fol_frac(cdim(2)),stem_frac(cdim(3)) &
                     ,DS_root(cdim(4)),root_frac(cdim(5)),DS_LRLV(cdim(6)) &
                     ,LRLV(cdim(7)),DS_LRRT(cdim(8)),LRRT(cdim(9))

  ! local variables
  integer i
  ! vector of ecosystem pools
  double precision, dimension((nodays+1),nopools) :: POOLS
  ! vector of ecosystem fluxes
  double precision, dimension(nodays,nofluxes) :: FLUXES
  double precision :: sumNPP
  double precision, dimension(nodays) :: lai & ! leaf area index
                                        ,GPP & ! Gross primary productivity
                                        ,NEE   ! net ecosystem exchange of CO2

  ! zero initial conditions
  lai = 0d0 ; GPP = 0d0 ; NEE = 0d0 ; POOLS = 0d0 ; FLUXES = 0d0 ; out_var = 0d0

  ! things that need to be allocated
  allocate(leftDaughter(dim1,dim2),rightDaughter(dim1,dim2),nodestatus(dim1,dim2), &
           xbestsplit(dim1,dim2),nodepred(dim1,dim2),bestvar(dim1,dim2))
  leftDaughter=tmp1 ; rightDaughter=tmp2 ; nodestatus=tmp3
  xbestsplit=tmp4   ; nodepred=tmp5      ; bestvar=tmp6
  nos_trees = dummy_nos_trees
  nos_inputs = dummy_nos_inputs

  ! generate deltat step from input data
  deltat(1) = met(1,1)
  do i = 2, nodays
     deltat(i)=met(1,i)-met(1,(i-1))
  end do

  ! begin iterations
  do i = 1, nos_iter

     if (pft == 1 .and. pft_specific == 1) then
         ! crop pft and we want pft specific model
         call CARBON_MODEL_CROP(1,nodays,met,pars(1:nopars,i),deltat,nodays,lat &
                          ,lai,NEE,FLUXES,POOLS,pft,nopars,nomet,nopools,nofluxes &
                          ,GPP,stock_seed_labile,DS_shoot,DS_root,fol_frac        &
                          ,stem_frac,root_frac,DS_LRLV,LRLV,DS_LRRT,LRRT)
      else ! must be generic
         call CARBON_MODEL(1,nodays,met,pars(1:nopars,i),deltat,nodays &
                          ,lat,lai,NEE,FLUXES,POOLS &
                          ,pft,nopars,nomet,nopools,nofluxes,GPP)
     end if ! pft / specific ecosystem model
if (i == 1) then
    open(unit=666,file="/home/lsmallma/out.csv", &
         status='replace',action='readwrite' )
    write(666,*),"lai",lai(1:365)
    close(666)
endif
write(666,*),i
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
     if (pft == 1 .and. pft_specific == 1) then
         out_var(i,1:nodays,13) = POOLS(1:nodays,7)
         ! add some pft specific carbon to the total C pool
         out_var(i,1:nodays,8) = out_var(i,1:nodays,8) + POOLS(1:nodays,7) &
                                + POOLS(1:nodays,8)
     else
         out_var(i,1:nodays,13) = 0d0
     endif

     ! calculate the actual NPP allocation fractions to foliar, wood and fine root pools
     ! by comparing the sum alloaction to each pools over the sum NPP.
     sumNPP = sum(FLUXES(1:nodays,1)*(1d0-pars(2,i))) ! GPP * (1-Ra) fraction
     out_var2(i,1) = sum(FLUXES(1:nodays,4)+FLUXES(1:nodays,8)) / sumNPP ! foliar
     out_var2(i,2) = sum(FLUXES(1:nodays,7)) / sumNPP ! wood
     out_var2(i,3) = sum(FLUXES(1:nodays,6)) / sumNPP ! fine root

  end do

  deallocate(leftDaughter,rightDaughter,nodestatus,&
             xbestsplit,nodepred,bestvar)

  ! return back to the subroutine then
  return

end subroutine ratdalec
