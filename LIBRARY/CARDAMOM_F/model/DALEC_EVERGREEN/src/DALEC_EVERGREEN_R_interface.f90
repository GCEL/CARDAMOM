

subroutine rdalecevergreen(output_dim,aNPP_dim,MTT_dim,SS_dim,met,pars &
                          ,out_var,out_var2,out_var3,out_var4,out_var5 &
                          ,out_var6,lat,nopars,nomet &
                          ,nofluxes,nopools,pft,pft_specific,nodays,noyears,deltat &
                          ,nos_iter)

  use CARBON_MODEL_MOD, only: CARBON_MODEL, extracted_C

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
  allocate(extracted_C(nodays+1))

  ! generate deltat step from input data
  deltat(1) = met(1,1)
  do i = 2, nodays
     deltat(i)=met(1,i)-met(1,(i-1))
  end do
  ! number of years in analysis
  nos_years = nint(sum(deltat)/365.25d0)
  ! number of time steps per year
  steps_per_year = nodays/nos_years

  ! begin iterations
  do i = 1, nos_iter
     ! reset harvest variable
     extracted_C=0.
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
     out_var(i,1:nodays,4)  = FLUXES(1:nodays,13) + FLUXES(1:nodays,14) ! het resp
     out_var(i,1:nodays,5)  = NEE
     out_var(i,1:nodays,6)  = POOLS(1:nodays,3) ! wood
     out_var(i,1:nodays,7)  = POOLS(1:nodays,5) ! som
     out_var(i,1:nodays,8)  = POOLS(1:nodays,1) + POOLS(1:nodays,2) + POOLS(1:nodays,3) ! Biomass
     out_var(i,1:nodays,9)  = POOLS(1:nodays,2) ! root
     out_var(i,1:nodays,10) = POOLS(1:nodays,4) ! litter
     out_var(i,1:nodays,11) = 0d0 ! labile
     out_var(i,1:nodays,12) = POOLS(1:nodays,1) ! foliage
     out_var(i,1:nodays,13) = extracted_C(1:nodays) ! harvested material
     out_var(i,1:nodays,14) = FLUXES(1:nodays,17) ! Fire value

     ! calculate the actual NPP allocation fractions to foliar, wood and fine root pools
     ! by comparing the sum alloaction to each pools over the sum NPP.
     sumNPP = sum(FLUXES(1:nodays,1)*(1-pars(2,i))) ! GPP * (1-Ra) fraction
     airt_adj = sum(met(3,1:nodays)) / dble(nodays)
     airt_adj = exp(pars(10,i)*airt_adj)
     out_var2(i,1) = sum(FLUXES(1:nodays,4)) / sumNPP ! foliar
     out_var2(i,2) = sum(FLUXES(1:nodays,6)) / sumNPP ! fine root
     out_var2(i,3) = sum(FLUXES(1:nodays,7)) / sumNPP ! wood
     ! Mean transit time
     out_var3(i,1) = pars(5,i) ! leaf life span (years)
     out_var3(i,2) = pars(7,i) ! root residence time (/day)
     out_var3(i,3) = pars(6,i) ! wood residence time (/day)
     out_var3(i,4) = (pars(1,i) + pars(8,i)) * airt_adj ! litter residence time (/day)
     out_var3(i,5) = pars(9,i) * airt_adj ! som residence time (/day)

     ! Determine locations of zeros in pools to correct turnover calculation
     ! Foliage
     fol_hak = 0 ; fol_filter(1:nodays) = 1d0
     where (POOLS(1:nodays,1) == 0) ! protection against NaN from division by zero
           fol_hak = 1 ; fol_filter(1:nodays) = 0d0
     end where
     ! Fine roots
     root_hak = 0 ; root_filter(1:nodays) = 1d0
     where (POOLS(1:nodays,2) == 0) ! protection against NaN from division by zero
           root_hak = 1 ; root_filter(1:nodays) = 0d0
     end where
     ! Wood
     wood_hak = 0 ; wood_filter(1:nodays) = 1d0
     where (POOLS(1:nodays,3) == 0) ! protection against NaN from division by zero
            wood_hak = 1 ; wood_filter(1:nodays) = 0d0
     end where
     ! Fol+root litter
     lit_hak = 0 ; lit_filter(1:nodays) = 1d0
     where (POOLS(1:nodays,4) == 0) ! protection against NaN from division by zero
            lit_hak = 1 ; lit_filter(1:nodays) = 0d0
     end where
     ! Soil
     som_hak = 0 ; som_filter(1:nodays) = 1d0
     where (POOLS(1:nodays,5) == 0) ! protection against NaN from division by zero
           som_hak = 1 ; som_filter(1:nodays) = 0d0
     end where
     ! Now loop through each year to estimate the annual residence time
     do y = 1, nos_years
        ! Estimate time steps covered by this year
        y_s = 1 + (steps_per_year * (y-1)) ; y_e = steps_per_year * y
        ! Foliage
        out_var5(i,1,y) = sum( ((FLUXES(y_s:y_e,10)+FLUXES(y_s:y_e,19)+FLUXES(y_s:y_e,25)) &
                                / POOLS(y_s:y_e,1)) * fol_filter(y_s:y_e)) / dble(steps_per_year-sum(fol_hak(y_s:y_e)))
        ! Fine roots
        out_var5(i,2,y) = sum( ((FLUXES(y_s:y_e,12)+FLUXES(y_s:y_e,20)+FLUXES(y_s:y_e,26)) &
                                / POOLS(y_s:y_e,2)) * root_filter(y_s:y_e)) / dble(steps_per_year-sum(root_hak(y_s:y_e)))
        ! Wood
        out_var5(i,3,y) = sum( ((FLUXES(y_s:y_e,11)+FLUXES(y_s:y_e,21)+FLUXES(y_s:y_e,27)) &
                                / POOLS(y_s:y_e,3)) * wood_filter(y_s:y_e)) / dble(steps_per_year-sum(wood_hak(y_s:y_e)))
        ! Litter (fol+fine root)
        out_var5(i,4,y) = sum( ((FLUXES(y_s:y_e,13)+FLUXES(y_s:y_e,15)+FLUXES(y_s:y_e,22)+FLUXES(y_s:y_e,28)) &
                                / POOLS(y_s:y_e,4)) * lit_filter(y_s:y_e)) / dble(steps_per_year-sum(lit_hak(y_s:y_e)))
        ! Soil
        out_var5(i,5,y) = sum( ((FLUXES(y_s:y_e,14)+FLUXES(y_s:y_e,23)) &
                                / POOLS(y_s:y_e,5)) * som_filter(y_s:y_e)) / dble(steps_per_year-sum(som_hak(y_s:y_e)))
     end do

     ! Calculate the mean inputs to each pool, needed for steady state calculation
     out_var4(i,1) = sum(FLUXES(1:nodays,4)) ! fol
     out_var4(i,2) = sum(FLUXES(1:nodays,6)) ! root
     out_var4(i,3) = sum(FLUXES(1:nodays,7)) ! wood
     out_var4(i,4) = sum(FLUXES(1:nodays,10)+FLUXES(1:nodays,12))! lit
     out_var4(i,5) = sum(FLUXES(1:nodays,15))! som

  end do ! nos_iter loop

  ! MTT - Convert daily fractional loss to years
  out_var3(:,2:5) = (out_var3(:,2:5)*365.25d0)**(-1d0) ! iter,(fol,root,wood,lit,som)
  out_var5 = (out_var5*365.25d0)**(-1d0)               ! iter,(fol,root,wood,lit,som)
  ! Assume same output for natural MRT as no disturbance is coded
  out_var6 = out_var3

  ! Steady state gC/m2
  out_var4 = (out_var4 / dble(nodays)) * 365.25d0 ! convert to annual mean input
  out_var4 = out_var4 * out_var3     ! multiply by residence time in years

  ! deallocate harvested variable
  deallocate(extracted_C)

  ! return back to the subroutine then
  return

end subroutine rdalecevergreen
