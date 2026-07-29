
module model_likelihood_module
  implicit none

  !!!!!!!!!!!
  ! Authorship contributions
  !
  ! This code is based on the original C verion of the University of Edinburgh
  ! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
  ! All code translation into Fortran, integration into the University of
  ! Edinburgh CARDAMOM code and subsequent modifications by:
  ! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
  ! See function / subroutine specific comments for exceptions and contributors
  !!!!!!!!!!!

  ! make all private
  private

  ! which to make open
  public :: model_likelihood, find_edc_initial_values

  ! declare needed types
  type EDCDIAGNOSTICS
    integer :: EDC
    integer :: DIAG
    integer :: PASSFAIL(100) ! allow space for 100 possible checks
    integer :: nedc ! number of edcs being assessed
  end type
  type (EDCDIAGNOSTICS), save :: EDCD

  ! Has the model sanity check been conducted yet?
  logical :: sanity_check = .false.

  contains
  !
  !------------------------------------------------------------------
  !
  subroutine find_edc_initial_values
    use MCMCOPT, only: PI, MCOUT, MCO
    use cardamom_structures, only: DATAin ! will need to change due to circular dependance
    use cardamom_io, only: restart_flag
    use MHMCMC_MODULE, only: MHMCMC

    ! subroutine deals with the determination of initial parameter and initial
    ! conditions which are consistent with EDCs

    implicit none

    ! declare local variables
    integer :: n, counter_local, EDC_iter
    double precision :: PEDC, PEDC_prev, ML, ML_prior, P_target
    double precision, dimension(PI%npars+1) :: EDC_pars

    ! set MCMC options needed for EDC run
    MCO%APPEND = 0
    MCO%nADAPT = 500
    MCO%fADAPT = 1d0
    MCO%nOUT = 100000
    MCO%nPRINT = 0
    MCO%nWRITE = 0
    ! the next two lines ensure that parameter inputs are either given or
    ! entered as -9999
    MCO%randparini = .true.
    MCO%returnpars = .true.
    MCO%fixedpars  = .true. ! TLS: changed from .false. for testing 16/12/2019

    ! Set initial priors to vector...
    PI%parini(1:PI%npars) = DATAin%parpriors(1:PI%npars)
    ! ... and assume we need to find random parameters
    PI%parfix = 0
    ! Target likelihood allows for controlling when the MCMC will stop
    P_target = 0d0

    ! if the prior is not missing and we have not told the edc to be random
    ! keep the value
!    do n = 1, PI%npars
!       if (PI%parini(n) /= -9999d0 .and. DATAin%edc_random_search < 1) PI%parfix(n) = 1
!    end do ! parameter loop

    ! set the parameter step size at the beginning
    PI%parvar = 1d0 ; PI%Nparvar = 0d0
    PI%use_multivariate = .false.
    ! Covariance matrix cannot be set to zero therefore set initial value to a
    ! small positive value along to variance access
    PI%covariance = 0d0 ; PI%mean_par = 0d0 ; PI%cov = .false.
    do n = 1, PI%npars
       PI%covariance(n,n) = 1d0
    end do

    ! if this is not a restart run, i.e. we do not already have a starting
    ! position we must being the EDC search procedure to find an ecologically
    ! consistent initial parameter set
    if (.not. restart_flag) then

        ! set up edc log likelihood for MHMCMC initial run
        PEDC_prev = -1000d0 ; PEDC = -1d0 ; counter_local = 0
        do while (PEDC < 0d0)

           write(*,*)"Beginning EDC search attempt"
           ! call the MHMCMC directing to the appropriate likelihood
           call MHMCMC(P_target,EDC_MODEL_LIKELIHOOD)

           ! store the best parameters from that loop
           PI%parini(1:PI%npars) = MCOUT%best_pars(1:PI%npars)
           ! turn off random selection for initial values
           MCO%randparini = .false.

           ! call edc likelihood function to get final edc probability
           call edc_model_likelihood(PI%parini,PEDC,ML_prior)

           ! keep track of attempts
           counter_local = counter_local + 1
           ! periodically reset the initial conditions
           if (PEDC < 0d0 .and. PEDC <= PEDC_prev .and. counter_local > 5) then
               ! Reset the previous EDC likelihood score
               PEDC_prev = -1000d0
               ! Reset parameters back to default
               PI%parini(1:PI%npars) = DATAin%parpriors(1:PI%npars)
               ! reset to select random starting point
               MCO%randparini = .true.
               ! reset the parameter step size at the beginning of each attempt
               PI%parvar = 1d0 ; PI%Nparvar = 0d0
               ! Covariance matrix cannot be set to zero therefore set initial value to a
               ! small positive value along to variance access
               PI%covariance = 0d0 ; PI%mean_par = 0d0 ; PI%cov = .false.
               PI%use_multivariate = .false.
               do n = 1, PI%npars
                  PI%covariance(n,n) = 1d0
               end do
           endif

        end do ! for while condition

    endif ! if for restart

    ! reset so that currently saved parameters will be used
    ! starting point in main MCMC
    PI%parfix(1:PI%npars) = 0
    MCOUT%best_pars = 0d0

  end subroutine find_edc_initial_values
  !
  !------------------------------------------------------------------
  !
  subroutine edc_model_likelihood(PARS, prob_out)
    use cardamom_structures, only: DATAin
    use MCMCOPT, only: PI
    use CARBON_MODEL_MOD, only: CARBON_MODEL
    use CARBON_MODEL_CROP_MOD, only: CARBON_MODEL_CROP

    ! Model likelihood function specifically intended for the determination of
    ! appropriate initial parameter choices, consistent with EDCs for DALEC2 /
    ! DALEC_CDEA

    implicit none

    ! declare inputs
    double precision, dimension(PI%npars), intent(inout) :: PARS
    ! output
    double precision, intent(inout) :: prob_out

    ! declare local variables
    integer ::  n
    double precision :: tot_exp, ML, exp_orig, decay_coef, prob_exp, EDC, EDC1, EDC2, infini

    ! set initial values
    EDCD%DIAG=1

    ! Perform a more aggressive sanity check which compares the bulk difference
    ! in all fluxes and pools from multiple runs of the same parameter set
    if (.not.sanity_check) call model_sanity_check(PI%parini)

    if (DATAin%ID == 4 .and. DATAin%PFT == 1) then
       ! call EDCs which can be evaluated prior to running the model
       call EDC1_CROP(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)
       if (DATAin%EDC == 0 .or. EDC1 == 1) then
           ! next need to run the model itself
           call CARBON_MODEL_CROP(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat,DATAin%nodays     &
                          ,DATAin%LAT, DATAin%M_LAI, DATAin%M_NEE          &
                          ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%pft       &
                          ,DATAin%nopars,DATAin%nomet,DATAin%nopools       &
                          ,DATAin%nofluxes,DATAin%M_GPP                    &
                          ,PI%stock_seed_labile,PI%DS_shoot,PI%DS_root     &
                          ,PI%fol_frac,PI%stem_frac,PI%root_frac,PI%DS_LRLV&
                          ,PI%LRLV,PI%DS_LRRT,PI%LRRT)

           ! assess post running EDCs
           call EDC2_CROP(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools &
                         ,DATAin%nodays,DATAin%deltat,PI%parmax,PARS,DATAin%MET &
                         ,DATAin%M_LAI,DATAin%M_NEE,DATAin%M_GPP,DATAin%M_POOLS &
                         ,DATAin%M_FLUXES,DATAin%meantemp,EDC2)
       else
           ! else assume edc fail
           EDC2 = 0
       endif


    else
       ! call EDCs which can be evaluated prior to running the model
       call EDC1_TESSEL(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)
       if (DATAin%EDC == 0 .or. EDC1 == 1) then
           ! next need to run the model itself
           call CARBON_MODEL(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat,DATAin%nodays     &
                          ,DATAin%LAT, DATAin%M_LAI, DATAin%M_NEE          &
                          ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%pft       &
                          ,DATAin%nopars,DATAin%nomet,DATAin%nopools       &
                          ,DATAin%nofluxes,DATAin%M_GPP)

           ! assess post running EDCs
           call EDC2_TESSEL(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools &
                           ,DATAin%nodays,DATAin%deltat,PI%parmax,PARS,DATAin%MET &
                           ,DATAin%M_LAI,DATAin%M_NEE,DATAin%M_GPP,DATAin%M_POOLS &
                           ,DATAin%M_FLUXES,DATAin%meantemp,EDC2)
       else
           ! assume EDC fail
           EDC2 = 0
       endif

    endif ! which model?

    ! combine results
    if (DATAin%EDC == 1 .and. (EDC1 == 0 .or. EDC2 == 0 .or. &
        sum(DATAin%M_LAI) /= sum(DATAin%M_LAI) .or. sum(DATAin%M_GPP) /= sum(DATAin%M_GPP))) then
        EDC=0
    else
        EDC=1
    end if

    ! calculate the likelihood
    tot_exp=0.
    do n = 1, EDCD%nedc
       tot_exp=tot_exp+(1.-EDCD%PASSFAIL(n))
!       if (EDCD%PASSFAIL(n) /= 1) print*,"failed edcs are: ", n
    end do ! checking EDCs
    ! for testing purposes, stop the model when start achieved
!    if (sum(EDCD%PASSFAIL) == 100) stop

    ! convert to a probability
    prob_out=-0.5*(tot_exp*10.)*DATAin%EDC

    ! override probability if parameter set gives NaN or near -infinitiy output
    call model_likelihood(PARS,ML)

    infini=0d0
    if (DATAin%EDC == 0 .and. (ML /= ML .or. ML == log(infini) .or. ML == -log(infini) )) then
       prob_out=prob_out-0.5*10d0
    end if

    ! now add the exponential component
    ! prob_out is the Log-Likelihood
    prob_out=prob_out

  end subroutine edc_model_likelihood
  !
  !------------------------------------------------------------------
  !
  subroutine model_sanity_check(PARS)
    use cardamom_structures, only: DATAin
    use MCMCOPT, only: PI
    use CARBON_MODEL_MOD, only: carbon_model
    use CARBON_MODEL_CROP_MOD, only: carbon_model_crop

    ! Carries out multiple carbon model iterations using the same parameter set
    ! to ensure that model outputs are consistent between iterations, i.e. that
    ! the model is numerically secure. Reproducible outputs from the models is
    ! essential for successful mcmc anlaysis

    implicit none

    ! Arguments
    double precision, dimension(PI%npars), intent(in) :: PARS

    ! Local arguments
    integer :: i
    double precision, dimension((DATAin%nodays+1),DATAin%nopools) :: local_pools
    double precision, dimension(DATAin%nodays,DATAin%nofluxes) :: local_fluxes
    double precision :: pool_error, flux_error

    ! Run model

    ! crop or not split....trouble
    if (DATAin%PFT == 1) then
       ! PFT has been provided and is crop! Best try running the crop model
       ! then...

       ! next need to run the model itself
       call CARBON_MODEL_CROP(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                             ,DATAin%nodays,DATAin%LAT,DATAin%M_LAI,DATAin%M_NEE &
                             ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%pft   &
                             ,DATAin%nopars,DATAin%nomet,DATAin%nopools   &
                             ,DATAin%nofluxes,DATAin%M_GPP                &
                             ,PI%stock_seed_labile,PI%DS_shoot,PI%DS_root &
                             ,PI%fol_frac,PI%stem_frac,PI%root_frac,PI%DS_LRLV&
                             ,PI%LRLV,PI%DS_LRRT,PI%LRRT)

       call CARBON_MODEL_CROP(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                             ,DATAin%nodays,DATAin%LAT,DATAin%M_LAI,DATAin%M_NEE &
                             ,local_fluxes,local_pools,DATAin%pft   &
                             ,DATAin%nopars,DATAin%nomet,DATAin%nopools   &
                             ,DATAin%nofluxes,DATAin%M_GPP                &
                             ,PI%stock_seed_labile,PI%DS_shoot,PI%DS_root &
                             ,PI%fol_frac,PI%stem_frac,PI%root_frac,PI%DS_LRLV&
                             ,PI%LRLV,PI%DS_LRRT,PI%LRRT)

    else

        ! next need to run the model itself
        call carbon_model(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                       ,DATAin%nodays,DATAin%LAT,DATAin%M_LAI,DATAin%M_NEE &
                       ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%nopars &
                       ,DATAin%nomet,DATAin%nopools,DATAin%nofluxes  &
                       ,DATAin%M_GPP)

        ! next need to run the model itself
        call carbon_model(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                       ,DATAin%nodays,DATAin%LAT,DATAin%M_LAI,DATAin%M_NEE &
                       ,local_fluxes,local_pools,DATAin%nopars &
                       ,DATAin%nomet,DATAin%nopools,DATAin%nofluxes  &
                       ,DATAin%M_GPP)
    end if ! crop or not if

    ! Compare outputs
    flux_error = sum(abs(DATAin%M_FLUXES - local_fluxes))
    pool_error = sum(abs(DATAin%M_POOLS - local_pools))
    ! If error between runs exceeds precision error then we have a problem
    if (pool_error > (tiny(0d0)*(DATAin%nopools*DATAin%nodays)) .or. &
        flux_error > (tiny(0d0)*(DATAin%nofluxes*DATAin%nodays))) then
        print*,"Error: multiple runs of the same parameter set indicates an error"
        print*,"Cumulative POOL error = ",pool_error
        print*,"Cumulative FLUX error = ",flux_error
        do i = 1,DATAin%nofluxes
           print*,"Sum abs error over time: flux = ",i
           print*,sum(abs(DATAin%M_FLUXES(:,i) - local_fluxes(:,i)))
        end do
        do i = 1, DATAin%nopools
           print*,"Sum abs error over time: pool = ",i
           print*,sum(abs(DATAin%M_POOLS(:,i) - local_pools(:,i)))
        end do
        stop
    end if

    ! Set Sanity check as completed
    sanity_check = .true.

  end subroutine model_sanity_check
  !
  !------------------------------------------------------------------
  !
  subroutine EDC1_CROP (PARS, npars, meantemp, meanrad, EDC1)

    ! the first of two subroutine to assess current parameters for passing realism tests for crop
    ! ecosystems

    implicit none

    ! declare input variables
    integer, intent(in) :: npars ! number of parameters
    double precision, intent(out) :: EDC1    ! EDC1 flag
    double precision, dimension(npars), intent(in) :: PARS ! current parameter set
    double precision, intent(in) :: meantemp & ! mean temperature (k)
                                   ,meanrad    ! mean radiation (MJ.m-2.day-1)

    ! declare local variables
    integer :: n, DIAG
    double precision :: torfol ! yearly leaf loss fraction

    ! set initial value
!    torfol=1./(pars(5)*365.25)
    EDC1=1
    DIAG=EDCD%DIAG

    ! set all EDCs to 1 (pass)
    EDCD%nedc=100
    EDCD%PASSFAIL(1:EDCD%nedc)=1

    !
    ! begin checking EDCs
    !

    ! Turnover of litter faster than turnover of som
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(9) > pars(10))) then
        EDC1=0 ; EDCD%PASSFAIL(1)=0
    endif

    ! decomposition of litter to SOM greater than SOM to air
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(9) > pars(1))) then
        EDC1=0 ; EDCD%PASSFAIL(2)=0
    endif

    ! turnover of foliage faster than turnover of wood
! TLS: turnover off because foliage and stem turnovers are made same
!    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(6) > pars(5)) then
!       EDC1=0 ; EDCD%PASSFAIL(3)=0
!    end if

    ! pre_DR should be greater than post_DR
!    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(4) > pars(3))) then
!        EDC1=0 ; EDCD%PASSFAIL(4)=0
!    endif

    ! for development: Tmin should be < topt and topt should be < tmax
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(26) > pars(28) &
                                     .or. pars(28) > pars(27) .or. pars(26) > pars(27))) then
        EDC1=0 ; EDCD%PASSFAIL(5)=0
    endif

    ! for development: the difference between each Tmin,Topt,Tmax > 1.
    if ((EDC1 == 1 .or. DIAG == 1) .and. (abs(pars(26)-pars(28)) < 1d0 &
                                    .or. abs(pars(28)-pars(27)) < 1d0 .or. abs(pars(26)-pars(27)) < 1d0)) then
        EDC1=0 ; EDCD%PASSFAIL(6)=0
    endif

   ! for vernalisation: Tmin < Topt < Tmax
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(29) > pars(31) .or. pars(31) > pars(30) .or. pars(29) > pars(30))) then
        EDC1=0 ; EDCD%PASSFAIL(7)=0
    endif

   ! for vernalisation: the difference between each Tmin, Topt, Tmax
    if ((EDC1 == 1 .or. DIAG == 1) .and. ( abs(pars(29)-pars(31)) < 1. .or. abs(pars(31)-pars(30)) < 1. .or. abs(pars(29)-pars(30)) < 1. ) ) then
        EDC1=0 ; EDCD%PASSFAIL(8)=0
    endif

   ! development temperature value should be larger corresponding vernalisation
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(29) > pars(26) .or. pars(31) > pars(28) .or. pars(30) > pars(27))) then
        EDC1=0 ; EDCD%PASSFAIL(9)=0
    endif

    ! plough must be before sow and after harvest: WINTER ONLY
!    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(16) > pars(12) ) then
!        EDC1=0 ; EDCD%PASSFAIL(10)=0
!    endif

!    ! harvest cannot be more than 345 after harvest
!    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(15) < ( pars(12)+345.25 )) ) then
!        EDC1=0 ; EDCD%PASSFAIL(10)=0
!    endif

!    ! plough must be before sow and after harvest: WINTER ONLY
!    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(16) < pars(15)) .or. (pars(16) > pars(12)) ) then
!        EDC1=0 ; EDCD%PASSFAIL(11)=0
!    endif

    ! could and probably should add some more

  end subroutine EDC1_CROP
  !
  !------------------------------------------------------------------
  !
  subroutine EDC2_CROP(npars,nomet,nofluxes,nopools,nodays,deltat &
                      ,parmax,pars,met,M_LAI,M_NEE,M_GPP,M_POOLS,M_FLUXES &
                      ,meantemp,EDC2)

    use CARBON_MODEL_CROP_MOD, only: resp_rate_temp_coeff,ts_length &
                                    ,sec_in_day,sec_in_hour

    ! the second of two subroutines for assessing current parameters for passing
    ! realism tests for crop ecosystems

    implicit none

    ! declare input variables
    integer, intent(in) :: npars    & ! number of model parameters
                          ,nomet    & ! number of met drivers
                          ,nofluxes & ! number of fluxes from model
                          ,nopools  & ! number of pools in model
                          ,nodays     ! number of days in simulation

    double precision, intent(in) :: deltat(nodays)              & ! decimal day model interval
                                   ,pars(npars)                 & ! vector of current parameters
                                   ,parmax(npars)               & ! vector of the maximum parameter values
                                   ,met(nomet,nodays)           & ! array of met drivers
                                   ,M_LAI(nodays)               & ! LAI output from current model simulation
                                   ,M_NEE(nodays)               & ! NEE output from current model simulation
                                   ,M_GPP(nodays)               & ! GPP output from current model simulation
                                   ,M_POOLS((nodays+1),nopools) & ! time varying states of pools in current model simulation
                                   ,M_FLUXES(nodays,nofluxes)   & ! time varying fluxes from current model simulation model
                                   ,meantemp                      ! site mean temperature (oC)

    double precision, intent(out) :: EDC2 ! the response flag for the dynamical set of EDCs

    ! declare local variables
    integer :: n, DIAG, no_years
    double precision :: mean_pools(nopools), G, decay_coef, meangpp, EQF
    double precision, dimension(:), allocatable :: mean_annual_pools
    double precision :: fauto & ! Fractions of GPP to autotrophic respiration
             ,ffol  & ! Fraction of GPP to foliage
             ,flab  & ! Fraction of GPP to labile pool
             ,froot & ! Fraction of GPP to root
             ,flit  & !
             ,fwood & ! Fraction of GPP to wood
             ,fsom    ! fraction of GPP som under eqilibrium conditions

    ! set initial value
    fauto=pars(2)
    ffol=sum(M_FLUXES(:,4))/(sum(M_FLUXES(:,1))*fauto)
    flab=sum(M_FLUXES(:,5))/(sum(M_FLUXES(:,1))*fauto)
    froot=sum(M_FLUXES(:,6))/(sum(M_FLUXES(:,1))*fauto)
    fwood=sum(M_FLUXES(:,7))/(sum(M_FLUXES(:,1))*fauto)
    fsom=fwood+(froot+flab+ffol)*pars(1)/(pars(1)+pars(10))
    flit=(froot+flab+ffol)
    ! length of time step in hours..
    ts_length = ((sum(deltat)/(nodays)) * sec_in_day) / sec_in_hour

    ! update initial values
    DIAG=EDCD%DIAG
    ! give EDC2 an initial value
    EDC2 = 1

    ! SOM attractor - must be within a factor of 2 from Csom0
    ! eqiulibrium factor (in comparison with initial conditions)
    EQF=10.0

    ! initialise and then calculate mean gpp values
    meangpp=sum(M_GPP(1:nodays))/real(nodays)

    ! EDC 11 - SOM steady state within order magnitude of initial conditions
    if ((EDC2 == 1 .or. DIAG == 1) .and. ((meangpp*fsom)/(ts_length*pars(9)*0.5*exp(resp_rate_temp_coeff*meantemp))) > (pars(23)*EQF)) then
       EDC2 = 0 ; EDCD%PASSFAIL(11) = 0
    end if
    if ((EDC2 == 1 .or. DIAG == 1) .and. ((meangpp*fsom)/(ts_length*pars(9)*0.5*exp(resp_rate_temp_coeff*meantemp))) < (pars(23)/EQF)) then
       EDC2 = 0 ; EDCD%PASSFAIL(11) = 0
    endif

    ! EDC 12 - Litter steady state assumptions
    if ((EDC2 == 1 .or. DIAG == 1) .and. ((meangpp*flit)/(ts_length*pars(10)*0.5*exp(resp_rate_temp_coeff*meantemp))) > (pars(22)*EQF)) then
        EDC2 = 0 ; EDCD%PASSFAIL(12) = 0
    endif
    if ((EDC2 == 1 .or. DIAG == 1) .and. ((meangpp*flit)/(ts_length*pars(10)*0.5*exp(resp_rate_temp_coeff*meantemp))) < (pars(22)/EQF)) then
        EDC2 = 0 ; EDCD%PASSFAIL(12) = 0
    endif

    ! EDC 13
    ! assesses the exponential decay/growth of the Csom pool

    !  work out how many completed years there are in the system
    no_years=int(nint(sum(deltat)/365.25))

    ! only do this for the Csom pool
    do n = 1, 1 !nopools
       if (EDC2 == 1 .or. DIAG == 1) then
          decay_coef=expdecay2(M_POOLS(1:(nodays+1),6),n,deltat,1,nodays+1)
          ! next assess the decay coefficient for meetings the EDC criterion
          if (abs(-log(2.)/decay_coef) < (365.25*real(no_years)) .and.  decay_coef < 0. ) then
             EDC2 = 0 ; EDCD%PASSFAIL(13)=0
          end if ! EDC conditions
       end if ! EDC .or. DIAG condition
    end do ! pools loop

    ! EDC 14
    ! assesses the exponential decay/growth of the Clit pool

    ! only do this for the Clit pool
    do n = 1, 1 !nopools
       if (EDC2 == 1 .or. DIAG == 1) then
          decay_coef=expdecay2(M_POOLS(1:(nodays+1),5),n,deltat,1,nodays+1)
          ! next assess the decay coefficient for meetings the EDC criterion
          if (abs(-log(2.)/decay_coef) < (365.25*real(no_years)) .and.  decay_coef < 0. ) then
             EDC2 = 0 ; EDCD%PASSFAIL(14)=0
          end if ! EDC conditions
       end if ! EDC .or. DIAG condition
    end do ! pools loop

  end subroutine EDC2_CROP
  !
  !
  !------------------------------------------------------------------
  !
  subroutine EDC1_TESSEL(PARS, npars, meantemp, meanrad, EDC1)

    ! subroutine assessed the current parameter sets for passing ecological and
    ! steady state contraints (Bloom et al., 2014).

    implicit none

    ! declare input variables
    integer, intent(in) :: npars ! number of parameters
    double precision, intent(out) :: EDC1    ! EDC1 flag
    double precision, dimension(npars), intent(in) :: PARS ! current parameter set
    double precision, intent(in) :: meantemp & ! mean temperature (k)
                                   ,meanrad    ! mean radiation (MJ.m-2.day-1)

    ! declare local variables
    integer :: n, DIAG
    double precision :: fauto & ! Fractions of GPP to autotrophic respiration
             ,ffol  & ! Fraction of GPP to foliage
             ,flab  & ! Fraction of GPP to labile pool
             ,froot & ! Fraction of GPP to root
             ,fwood & ! Fraction of GPP to wood
             ,fsom    ! fraction of GPP som under eqilibrium conditions

    double precision :: torfol ! yearly leaf loss fraction

    ! set initial value
    torfol=1./(pars(5)*365.25)
    EDC1=1
    DIAG=EDCD%DIAG
    fauto=pars(2)
    ffol=(1.-fauto)*pars(3)
    flab=(1.-fauto-ffol)*pars(12)
    froot=(1.-fauto-ffol-flab)*pars(4)
    fwood=1.-fauto-ffol-flab-froot
    fsom=fwood+(froot+flab+ffol)*pars(1)/(pars(1)+pars(8))

    ! set all EDCs to 1 (pass)
    EDCD%nedc=100
    EDCD%PASSFAIL(1:EDCD%nedc)=1

    !
    ! begin checking EDCs
    !

    ! Turnover of litter faster than turnover of som
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(9) > pars(8))) then
        EDC1=0 ; EDCD%PASSFAIL(1)=0
    endif

    ! litter2som greater than som to atm rate
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(1) < pars(9))) then
       EDC1=0 ; EDCD%PASSFAIL(2)=0
    endif

    ! turnover of foliage faster than turnover of wood
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(6) > torfol) then
       EDC1=0 ; EDCD%PASSFAIL(3)=0
    end if

    ! root turnover greater than som turnover at mean temperature
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(7) < (pars(9)*exp(pars(10)*meantemp)))) then
       EDC1=0 ; EDCD%PASSFAIL(4)=0
    endif

    ! GPP allocation to foliage and labile cannot be 5 orders of magnitude
    ! difference from GPP allocation to roots
    if ((EDC1 == 1 .or. DIAG == 1) .and. ((ffol+flab) > (5.*froot) .or. ((ffol+flab)*5.) < froot)) then
       EDC1=0 ; EDCD%PASSFAIL(5)=0
    endif

    ! could always add more / remove some

  end subroutine EDC1_TESSEL
  !
  !------------------------------------------------------------------
  !
  subroutine EDC2_TESSEL(npars,nomet,nofluxes,nopools,nodays,deltat &
                        ,parmax,pars,met,M_LAI,M_NEE,M_GPP,M_POOLS,M_FLUXES &
                        ,meantemp,EDC2)

    ! Determines whether the dynamical contraints for the search of the initial
    ! parameters has been successful or whether or not we should abandon the
    ! current set and move on

    implicit none

    ! declare input variables
    integer, intent(in) :: npars    & ! number of model parameters
                          ,nomet    & ! number of met drivers
                          ,nofluxes & ! number of fluxes from model
                          ,nopools  & ! number of pools in model
                          ,nodays     ! number of days in simulation

    double precision, intent(in) :: deltat(nodays)              & ! decimal day model interval
                                   ,pars(npars)                 & ! vector of current parameters
                                   ,parmax(npars)               & ! vector of the maximum parameter values
                                   ,met(nomet,nodays)           & ! array of met drivers
                                   ,M_LAI(nodays)               & ! LAI output from current model simulation
                                   ,M_NEE(nodays)               & ! NEE output from current model simulation
                                   ,M_GPP(nodays)               & ! GPP output from current model simulation
                                   ,M_POOLS((nodays+1),nopools) & ! time varying states of pools in current model simulation
                                   ,M_FLUXES(nodays,nofluxes)   & ! time varying fluxes from current model simulation model
                                   ,meantemp                      ! site mean temperature (oC)

    double precision, intent(out) :: EDC2 ! the response flag for the dynamical set of EDCs

    ! declare local variables
    integer :: n, DIAG, no_years, y, PEDC, nn, num_EDC
    double precision :: mean_pools(nopools), G, decay_coef, meangpp, EQF
    double precision, dimension(:), allocatable :: mean_annual_pools
    double precision :: fauto & ! Fractions of GPP to autotrophic respiration
             ,ffol  & ! Fraction of GPP to foliage
             ,flab  & ! Fraction of GPP to labile pool
             ,froot & ! Fraction of GPP to root
             ,fwood & ! Fraction of GPP to wood
             ,fsom  & ! fraction of GPP som under eqilibrium conditions
             ,flit    ! fraction of GPP to litter under equilibrium conditions

    ! update initial values
    DIAG=EDCD%DIAG
    EDC2=1
    fauto=pars(2)
    ffol=(1.-fauto)*pars(3)
    flab=(1.-fauto-ffol)*pars(12)
    froot=(1.-fauto-ffol-flab)*pars(4)
    fwood=1.-fauto-ffol-flab-froot
    fsom=fwood+(froot+flab+ffol)*pars(1)/(pars(1)+pars(8))
    flit=(froot+flab+ffol)

    ! derive mean pools
    do n = 1, nopools
       mean_pools(n)=cal_mean_pools(M_POOLS,n,nodays+1,nopools)
    end do

    !
    ! Begin EDCs here
    !

    ! EDC 6
    ! ensure ratio between Cfoilar and Croot is less than 5
    if ((EDC2 == 1 .or. DIAG == 1) .and. (mean_pools(2) > (mean_pools(3)*5.) .or. (mean_pools(2)*5.) < mean_pools(3)) ) then
        EDC2=0 ; EDCD%PASSFAIL(6)=0
    end if

    ! EDC 7
    ! Assess growth factor over 10 years; order magnitude changes are not
    ! allowed
    no_years=int(nint(sum(deltat)/365.25))
    G=0.1
    allocate(mean_annual_pools(no_years))

    ! generate mean annual pool values
    do n = 1, nopools
       ! Growth greater than factor G of pools is not allowed, increase is restricted by G growth
       ! over N years. e.g. G = 0.1 is order magnitude change over 10 year
       ! period
       do y = 1, no_years
          ! derive mean annual pools
          mean_annual_pools(y)=cal_mean_annual_pools(M_POOLS,y,n,nopools,deltat,nodays+1)
       end do ! year loop
       ! now check the growth rate
       if ((EDC2 == 1 .or. DIAG == 1) .and. ((mean_annual_pools(no_years)/mean_annual_pools(1)) > (1.+G*real(no_years)))) then
          EDC2=0 ; EDCD%PASSFAIL(7)=0
       endif
    end do ! pool loop

    ! done now so clean up
    deallocate(mean_annual_pools)
    ! EDC 8
    ! assesses the exponential decay/growth of each model pool

    ! loop through each pool in turn
    do n = 1, nopools
       if (EDC2 == 1 .or. DIAG == 1) then
          decay_coef=expdecay2(M_POOLS,n,deltat,nopools,nodays+1)
          ! next assess the decay coefficient for meetings the EDC criterion
          if (abs(-log(2.)/decay_coef) < (365.25*real(no_years)) .and. decay_coef < 0. ) then
             EDC2 = 0 ; EDCD%PASSFAIL(8)=0
          end if ! EDC conditions
       end if ! EDC .or. DIAG condition
    end do ! pools loop

    ! SOM attractor - must be within a factor of 2 from Csom0
    ! eqiulibrium factor (in comparison with initial conditions)
    EQF=10.0

    ! initialise and then calculate mean gpp values
    meangpp=sum(M_GPP(1:nodays))/real(nodays)

    ! EDC 9 - SOM steady state within order magnitude of initial conditions
    if ((EDC2 == 1 .or. DIAG == 1) .and. ((meangpp*fsom)/(pars(9)*exp(pars(10)*meantemp))) > (pars(22)*EQF)) then
       EDC2 = 0 ; EDCD%PASSFAIL(9) = 0
    end if
    if ((EDC2 == 1 .or. DIAG == 1) .and. ((meangpp*fsom)/(pars(9)*exp(pars(10)*meantemp))) < (pars(22)/EQF)) then
       EDC2 = 0 ; EDCD%PASSFAIL(9) = 0
    endif

    ! EDC 10 - Litter steady state assumptions
    if ((EDC2 == 1 .or. DIAG == 1) .and. ((meangpp*flit)/(pars(8)*exp(pars(10)*meantemp))) > (pars(21)*EQF)) then
        EDC2 = 0 ; EDCD%PASSFAIL(10) = 0
    endif
    if ((EDC2 == 1 .or. DIAG == 1) .and. ((meangpp*flit)/(pars(8)*exp(pars(10)*meantemp))) < (pars(21)/EQF)) then
        EDC2 = 0 ; EDCD%PASSFAIL(10) = 0
    endif

    ! EDC 11 - Wood steady state assumptions
    if ((EDC2 == 1 .or. DIAG == 1) .and. ((meangpp*fwood)/pars(6)) > (pars(20)*EQF)) then
        EDC2 = 0 ; EDCD%PASSFAIL(11) = 0
    end if
    if ((EDC2 == 1 .or. DIAG == 1) .and. ((meangpp*fwood)/pars(6)) < (pars(20)/EQF)) then
        EDC2 = 0 ; EDCD%PASSFAIL(11) = 0
    endif

    ! EDC 12 - Root steady state assumptions
    if ((EDC2 == 1 .or. DIAG == 1) .and. ((meangpp*froot)/pars(7)) > (pars(19)*EQF)) then
        EDC2 = 0 ; EDCD%PASSFAIL(12) = 0
    endif
    if ((EDC2 == 1 .or. DIAG == 1) .and. ((meangpp*froot)/pars(7)) < (pars(19)/EQF)) then
        EDC2 = 0 ; EDCD%PASSFAIL(12) = 0
    endif

    !
    ! EDCs done, below are additional fault detection conditions
    !

    ! additional faults can be stored in locations 35 - 40 of the PASSFAIL array

    ! ensure minimum pool values are >= 0 and /= NaN
    if (EDC2 == 1 .or. DIAG == 1) then
       n=1
       do while (n <= nopools .and. (EDC2 == 1 .or. DIAG == 1))
          nn = 1 ; PEDC = 1
          do while (nn <= (nodays+1) .and. PEDC == 1)
             ! now check conditions
             if (M_POOLS(nn,n) < 0. .or. M_POOLS(nn,n) /= M_POOLS(nn,n)) then
                 EDC2=0 ; PEDC=0 ; EDCD%PASSFAIL(35+n)=0
             end if ! less than zero and is NaN condition
          nn = nn + 1
          end do ! nn < nodays .and. PEDC == 1
          n = n + 1
       end do ! for nopools .and. EDC .or. DIAG condition
    end if ! min pool assessment

  end subroutine EDC2_TESSEL
  !
  !------------------------------------------------------------------
  !
  double precision function cal_mean_pools(pools,pool_number,averaging_period,nopools)

    ! Function calculate the mean values of model pools / states across the
    ! entire simulation run

    implicit none

    ! declare input variables
    integer, intent(in) :: nopools          & !
                          ,pool_number      & !
                          ,averaging_period   !

    double precision,dimension(averaging_period,nopools), intent (in) :: pools

    ! declare local variables
    integer :: c

    ! initial conditions
    cal_mean_pools=0.

    ! loop through now
    cal_mean_pools=sum(pools(1:averaging_period,pool_number))/real(averaging_period)

    ! ensure return command issued
    return

  end function cal_mean_pools
  !
  !------------------------------------------------------------------
  !
  double precision function cal_mean_annual_pools(pools,year,pool_number,nopools,interval,averaging_period)

    ! Function calculates the mean model pools values for each individual year
    ! in the simulation

    implicit none

    ! declare input variables
    integer, intent(in) :: nopools     & ! how many pools in the model
                          ,year        & ! which year are we working on
                          ,averaging_period & ! number of days in analysis period
                          ,pool_number   ! which pool are we currently working on

    double precision, intent(in) :: pools(averaging_period,nopools) & ! input pool state variables
                         ,interval((averaging_period-1))      ! model time step in decimal days

    ! declare local variables
    integer :: startday, endday, c

    ! initialise the output variable
    cal_mean_annual_pools=0.

    ! calculate some constants
    startday=floor(365.25*(year-1)/(sum(interval)/(averaging_period-1)))+1
    endday=floor(365.25*year/(sum(interval)/(averaging_period-1)))

    ! pool through and work out the annual mean values
    cal_mean_annual_pools=sum(pools(startday:endday,pool_number))/(endday-startday)

    ! ensure function returns
    return

  end function cal_mean_annual_pools
  !
  !------------------------------------------------------------------
  !
  double precision function cal_max_annual_pools(pools,year,pool_number,nopools,interval,averaging_period)

    ! Function calculates the max model pools values for each individual year
    ! in the simulation

    implicit none

    ! declare input variables
    integer, intent(in) :: nopools     & ! how many pools in the model
                          ,year        & ! which year are we working on
                          ,averaging_period & ! number of days in analysis period
                          ,pool_number   ! which pool are we currently working on

    double precision, intent(in) :: pools(averaging_period,nopools) & ! input pool state variables
                         ,interval((averaging_period-1))      ! model time step in decimal days

    ! declare local variables
    integer :: startday, endday, c

    ! initialise the output variable
    cal_max_annual_pools=0.

    ! calculate some constants
    startday=floor(365.25*(year-1)/(sum(interval)/(averaging_period-1)))+1
    endday=floor(365.25*year/(sum(interval)/(averaging_period-1)))


    ! pool through and work out the annual max values
    cal_max_annual_pools=maxval(pools(startday:endday,pool_number))

    ! ensure function returns
    return

  end function cal_max_annual_pools
  !
  !------------------------------------------------------------------
  !
  double precision function expdecay2(pools,pool_number,interval,nopools,averaging_period)

   ! Function to calculate the exponential decay coefficients used several EDCs.
   ! We assumpe the equation Cexp= a + b*exp(c*t)

   implicit none

   ! declare input variables
   integer, intent(in) :: nopools     & ! how many pools in the model
                         ,averaging_period & ! i.e. nodays + 1
                         ,pool_number   ! which pool are we currently working on

   double precision, intent(in) :: pools(averaging_period,nopools) & ! input pool state variables
                                  ,interval((averaging_period-1))      ! model time step in decimal days

   ! declare local variables
   integer :: n
   double precision :: P0    & ! initial pool value
            ,os,aw &
            ,MP0   & ! mean pool (year 1 to year end-2)
            ,MP1   & ! mean pool (year 2 to year end-1)
            ,MP0os & ! mean pool (year 1+os to year end-2+os)
            ,MP1os & ! mean pool (year 2+os to year end-2+os)
            ,dcdt1 & ! gradient of exponential over time in second year
            ,dcdt0   ! gradient of exponential over time in first year

   ! declare initial values / constants
   os = 1 ! offset in days
   aw = floor(365.25/(sum(interval)/(averaging_period-1))) ! averaging window
   MP0 = 0. ; MP1 = 0. ; MP0os = 0. ; MP1os = 0.

   ! calculate mean pools within defined averaging window
   do n = 1, int(aw)
     MP0=MP0+pools(n,pool_number)
   end do ! for first year
   ! now average
   MP0=MP0/aw

   do n = int(aw)+1, int(aw*2)
     MP1=MP1+pools(n,pool_number)
   end do ! for second year
   ! now average
   MP1=MP1/aw

   do n = (1+int(os)), int(aw+os)
     MP0os=MP0os+pools(n,pool_number)
   end do ! for first year with offset
   ! now average
   MP0os=MP0os/aw

   do n = (int(aw+os)+1), int(aw*2.+os)
     MP1os=MP1os+pools(n,pool_number)
   end do ! for second year withoffset
   ! now average
   MP1os=MP1os/aw

   ! derive mean gradient ratio (dcdt1/dcdt0)
   ! where dcdt1 is the numeric gradient between n+1 and n+365+1
   ! and dcdt0 os the numeric gradient between n and n+365
   dcdt1 = MP1os-MP0os
   dcdt0 = MP1-MP0

   ! using multiple year mean to determine c
   if ((dcdt1 > 0. .and. dcdt0 < 0.) .or. (dcdt1 < 0. .and. dcdt0 > 0.) &
       .or. dcdt1 == 0 .or. dcdt0 == 0) then
       ! then return error values
       expdecay2 = 1
   else
       expdecay2 = log(dcdt1/dcdt0) / (os*(sum(interval)/(averaging_period-1)))
   end if
   ! ensure return
   return

  end function expdecay2
  !
  !------------------------------------------------------------------
  !
  subroutine model_likelihood(PARS,ML_out)
    use MCMCOPT, only:  PI
    use cardamom_structures, only: DATAin
    use CARBON_MODEL_MOD, only: CARBON_MODEL
    use CARBON_MODEL_CROP_MOD, only: CARBON_MODEL_CROP

    ! this subroutine is responsible, under normal circumstances for the running
    ! of the DALEC model, calculation of the log-likelihood for comparison
    ! assessment of parameter performance and use of the EDCs if they are
    ! present / selected

    implicit none

    ! declare inputs
    double precision, dimension(PI%npars), intent(inout) :: PARS ! current parameter vector
    ! output
    double precision, intent(inout) :: ML_out ! output variables for log-likelihood

    ! declare local variables
    double precision :: EDC,EDC1,EDC2

    ! initial values
    ML_out=0.
    EDCD%DIAG=0

    if (DATAin%ID == 4 .and. DATAin%PFT == 1) then
        ! call EDCs which can be evaluated prior to running the model
        call EDC1_CROP(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)
    else
        ! call EDCs which can be evaluated prior to running the model
        call EDC1_TESSEL(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)
    endif

    ! now use the EDCD%EDC flag to determine if effect is kept
    if (DATAin%EDC == 1) then
        EDC = EDC1
    else
        EDC = 1
    end if

    ! update effect to the probabity
    ML_out=ML_out+log(EDC)

    ! if first set of EDCs have been passed
    if (EDC == 1) then
       ! calculate parameter log likelihood (assumed we have estimate of
       ! uncertainty)
       ML_out=ML_out+likelihood_p(PI%npars,DATAin%parpriors,DATAin%parpriorunc,PARS)

       if (DATAin%ID == 4 .and. DATAin%PFT == 1) then
           ! run the dalec model
           call CARBON_MODEL_CROP(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat,DATAin%nodays     &
                          ,DATAin%LAT, DATAin%M_LAI, DATAin%M_NEE                 &
                          ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%pft              &
                          ,DATAin%nopars,DATAin%nomet,DATAin%nopools              &
                          ,DATAin%nofluxes,DATAin%M_GPP                           &
                          ,PI%stock_seed_labile,PI%DS_shoot,PI%DS_root            &
                          ,PI%fol_frac,PI%stem_frac,PI%root_frac,PI%DS_LRLV       &
                          ,PI%LRLV,PI%DS_LRRT,PI%LRRT)

           ! check edc2
           call EDC2_CROP(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools &
                          ,DATAin%nodays,DATAin%deltat,PI%parmax,PARS,DATAin%MET &
                          ,DATAin%M_LAI,DATAin%M_NEE,DATAin%M_GPP,DATAin%M_POOLS &
                          ,DATAin%M_FLUXES,DATAin%meantemp,EDC2)

       else
           ! run the dalec model
           call CARBON_MODEL(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat,DATAin%nodays   &
                          ,DATAin%LAT, DATAin%M_LAI, DATAin%M_NEE          &
                          ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%pft       &
                          ,DATAin%nopars,DATAin%nomet,DATAin%nopools       &
                          ,DATAin%nofluxes,DATAin%M_GPP)

           ! check edc2
           call EDC2_TESSEL(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools &
                          ,DATAin%nodays,DATAin%deltat,PI%parmax,PARS,DATAin%MET &
                          ,DATAin%M_LAI,DATAin%M_NEE,DATAin%M_GPP,DATAin%M_POOLS &
                          ,DATAin%M_FLUXES,DATAin%meantemp,EDC2)

       endif

       ! check if EDCs are switched on
       if (DATAin%EDC == 1) then
           EDC = EDC2
       else
           EDC = 1
       end if

       ! extra checks to ensure correct running of the model
       if (sum(DATAin%M_LAI) /= sum(DATAin%M_LAI) .or. sum(DATAin%M_GPP) /= sum(DATAin%M_GPP)) then
           EDC=0
       end if

       ! add EDC2 log-likelihood
       ML_out=ML_out+log(EDC)

       ! calculate final model likelihood when compared to obs
       if (EDC == 1) then
          ML_out=ML_out+likelihood(PI%npars,PARS)
       endif ! EDC still == 1

    end if ! EDC == 1

  end subroutine model_likelihood
  !
  !------------------------------------------------------------------
  !
  double precision function likelihood_p(npars,parpriors,parpriorunc,pars)
    ! function calculates the parameter based log-likelihood for the current set
    ! of parameters. This assumes that we have any actual priors / prior
    ! uncertainties to be working with. This does include initial states, as we
    ! consider them to be parameters

    implicit none

    ! declare input variables
    integer, intent(in) :: npars
    double precision, dimension(npars), intent(in) :: pars      & ! current parameter vector
                                                     ,parpriors & ! prior values for parameters
                                                     ,parpriorunc ! prior uncertainties

    ! declare local variables
    integer :: n

    ! set initial value
    likelihood_p = 0.

    ! now loop through defined parameters for their uncertainties
    do n = 1, npars
       ! if there is actually a value
       if (parpriors(n) > -9999) then
           likelihood_p=likelihood_p-0.5*(log(pars(n)/parpriors(n))/log(parpriorunc(n)))**2.
       end if
    end do

    ! dont for get to return
    return

  end function likelihood_p
  !
  !------------------------------------------------------------------
  !
    double precision function likelihood(npars,pars)
    use cardamom_structures, only: DATAin

    ! calculates the likelihood of of the model output compared to the available
    ! observations which have been input to the model

    implicit none

    ! declare arguments
    integer, intent(in) :: npars
    double precision, dimension(npars), intent(in) :: pars

    ! declare local variables
    integer :: n, dn, no_years, y
    double precision :: tot_exp, pool_dynamics, tmp_var, infini
    double precision, allocatable :: mean_annual_pools(:)

    ! initial value
    likelihood=0d0 ; infini=0d0

    ! GPP Log-likelihood
    tot_exp = 0.
    if (DATAin%ngpp > 0) then
       do n = 1, DATAin%ngpp
         dn=DATAin%gpppts(n)
         ! note that division is the uncertainty
         tot_exp=tot_exp+((DATAin%M_GPP(dn)-DATAin%GPP(dn))/DATAin%GPP_unc(dn))**2d0
       end do
       likelihood=likelihood-0.5*tot_exp
    endif

    ! LAI log-likelihood
    tot_exp = 0.
    if (DATAin%nlai > 0) then
       do n = 1, DATAin%nlai
         dn=DATAin%laipts(n)
         ! if zero or greater allow calculation with min condition to prevent
         ! errors of zero LAI which occur in managed systems
         if (DATAin%M_LAI(dn) >= 0.) then
             ! note that division is the uncertainty
             !tot_exp = tot_exp+(log(max(0.001d0,DATAin%M_LAI(dn))/max(0.001d0,DATAin%LAI(dn)))/log(DATAin%LAI_unc(dn)))**2d0
             tot_exp = tot_exp + (max(0.001d0,DATAin%M_LAI(dn)-DATAin%LAI(dn))/DATAin%LAI_unc(dn))**2d0
         else
             ! if not then we have unrealistic negative values or NaN so indue
             ! error
             tot_exp=tot_exp+(-log(infini))
         endif
       end do
       likelihood=likelihood-0.5*tot_exp
    endif

    ! NEE likelihood
    tot_exp = 0.
    if (DATAin%nnee > 0) then
       do n = 1, DATAin%nnee
         dn=DATAin%neepts(n)
         ! note that division is the uncertainty
         tot_exp=tot_exp+((DATAin%M_NEE(dn)-DATAin%NEE(dn))/DATAin%NEE_unc(dn))**2d0
       end do
       likelihood=likelihood-0.5*tot_exp
    endif

    ! Reco likelihood
    tot_exp = 0.
    if (DATAin%nreco > 0) then
       do n = 1, DATAin%nreco
         dn=DATAin%recopts(n)
         tmp_var=DATAin%M_NEE(dn)+DATAin%M_GPP(dn)
         ! note that we calculate the Ecosystem resp from GPP and NEE
         tot_exp=tot_exp+((tmp_var-DATAin%Reco(dn))/DATAin%Reco_unc(dn))**2d0
       end do
       likelihood=likelihood-0.5*tot_exp
    endif

    ! Cwood increment log-likelihood
    tot_exp = 0.
    if (DATAin%nwoo > 0) then
       do n = 1, DATAin%nwoo
         dn=DATAin%woopts(n)
         ! note that division is the uncertainty
         tot_exp=tot_exp+(log((DATAin%M_POOLS(dn,4)-DATAin%M_POOLS(dn-365,4)) &
                          / DATAin%WOO(dn))/log(DATAin%WOO_unc(dn)))**2d0
       end do
       likelihood=likelihood-0.5*tot_exp
    endif

    ! Cfoliage log-likelihood
    tot_exp = 0.
    if (DATAin%nCfol_stock > 0) then
       do n = 1, DATAin%nCfol_stock
         dn=DATAin%Cfol_stockpts(n)
         ! note that division is the uncertainty
!         tot_exp=tot_exp+(log(DATAin%M_POOLS(dn,2)/DATAin%Cfol_stock(dn))/log(2.))**2.
         tot_exp=tot_exp+((DATAin%M_POOLS(dn,2)-DATAin%Cfol_stock(dn)) &
                          / (DATAin%Cfol_stock(dn)*DATAin%Cfol_stock_unc(dn)))**2d0
       end do
       likelihood=likelihood-0.5*tot_exp
    endif

    ! Annual foliar maximum
    tot_exp = 0.
    if (DATAin%nCfolmax_stock > 0) then
       no_years=int(nint(sum(DATAin%deltat)/365.25))
       if (allocated(mean_annual_pools)) deallocate(mean_annual_pools)
       allocate(mean_annual_pools(no_years))
       ! determine the annual max for each pool
       do y = 1, no_years
          ! derive mean annual foliar pool
           mean_annual_pools(y)=cal_max_annual_pools(DATAin%M_POOLS,y,2,DATAin%nopools,DATAin%deltat,DATAin%nodays+1)
       end do ! year loop
       ! loop through the observations then
       do n = 1, DATAin%nCfolmax_stock
         ! load the observation position in stream
         dn=DATAin%Cfolmax_stockpts(n)
         ! determine which years this in in for the simulation
         y = ceiling( (dble(dn)*(sum(DATAin%deltat)/(DATAin%nodays))) / 365.25 )
         ! load the correct year into the analysis
         tmp_var = mean_annual_pools(y)
         ! note that division is the uncertainty
         tot_exp=tot_exp+((tmp_var-DATAin%Cfolmax_stock(dn)) &
                          / (DATAin%Cfolmax_stock(dn)*DATAin%Cfolmax_stock_unc(dn)))**2d0
       end do
       likelihood=likelihood-0.5*tot_exp
    endif

    ! Cwood log-likelihood (i.e. branch, stem and CR)
    tot_exp = 0.
    if (DATAin%nCwood_stock > 0) then
       do n = 1, DATAin%nCwood_stock
         dn=DATAin%Cwood_stockpts(n)
         ! note that division is the uncertainty
!         tot_exp=tot_exp+(log(DATAin%M_POOLS(dn,4)/DATAin%Cwood_stock(dn))/log(2.))**2.
!         tot_exp=tot_exp+((DATAin%M_POOLS(dn,4)-DATAin%Cwood_stock(dn))/(DATAin%Cwood_stock(dn)*0.20))**2.
         tot_exp=tot_exp+((DATAin%M_POOLS(dn,4)-DATAin%Cwood_stock(dn))/(DATAin%Cwood_stock(dn)*DATAin%Cwood_stock_unc(dn)))**2d0
       end do
       likelihood=likelihood-0.5*tot_exp
    endif

    ! Cagb log-likelihood
    tot_exp = 0.
    if (DATAin%nCagb_stock > 0) then
       do n = 1, DATAin%nCagb_stock
         dn=DATAin%Cagb_stockpts(n)
         ! remove coarse root fraction from wood (pars29)
         tmp_var = DATAin%M_POOLS(dn,4)-(DATAin%M_POOLS(dn,4)*pars(29))
         tot_exp=tot_exp+((tmp_var-DATAin%Cagb_stock(dn))/(DATAin%Cagb_stock(dn)*DATAin%Cagb_stock_unc(dn)))**2d0
       end do
       likelihood=likelihood-0.5*tot_exp
    endif

    ! Cstem log-likelihood
    tot_exp = 0.
    if (DATAin%nCstem_stock > 0) then
       do n = 1, DATAin%nCstem_stock
         dn=DATAin%Cstem_stockpts(n)
         ! remove coarse root and branches from wood (pars29 and pars28)
         tmp_var = DATAin%M_POOLS(dn,4)-( (DATAin%M_POOLS(dn,4)*pars(29))+((DATAin%M_POOLS(dn,4)*pars(28))) )
         tot_exp=tot_exp+((tmp_var-DATAin%Cstem_stock(dn))/(DATAin%Cstem_stock(dn)*DATAin%Cstem_stock_unc(dn)))**2d0
       end do
       likelihood=likelihood-0.5*tot_exp
    endif

    ! Cbranch log-likelihood
    tot_exp = 0.
    if (DATAin%nCbranch_stock > 0) then
       do n = 1, DATAin%nCbranch_stock
         dn=DATAin%Cbranch_stockpts(n)
         ! extract branch component from only
         tmp_var = DATAin%M_POOLS(dn,4)*pars(28)
         tot_exp=tot_exp+((tmp_var-DATAin%Cbranch_stock(dn))/(DATAin%Cbranch_stock(dn)*DATAin%Cbranch_stock_unc(dn)))**2d0
       end do
       likelihood=likelihood-0.5*tot_exp
    endif

    ! Ccoarseroot log-likelihood
    tot_exp = 0.
    if (DATAin%nCcoarseroot_stock > 0) then
       do n = 1, DATAin%nCcoarseroot_stock
         dn=DATAin%Ccoarseroot_stockpts(n)
         ! extract coarse root component from wood only
         tmp_var = DATAin%M_POOLS(dn,4)*pars(29)
         tot_exp=tot_exp+((tmp_var-DATAin%Ccoarseroot_stock(dn))/(DATAin%Ccoarseroot_stock(dn)*DATAin%Ccoarseroot_stock_unc(dn)))**2d0
       end do
       likelihood=likelihood-0.5*tot_exp
    endif

    ! Croots log-likelihood
    tot_exp = 0.
    if (DATAin%nCroots_stock > 0) then
       do n = 1, DATAin%nCroots_stock
         dn=DATAin%Croots_stockpts(n)
         ! note that division is the uncertainty
!         tot_exp=tot_exp+(log(DATAin%M_POOLS(dn,3)/DATAin%Croots_stock(dn))/log(2.))**2.
         tot_exp=tot_exp+((DATAin%M_POOLS(dn,3)-DATAin%Croots_stock(dn)) &
                         / (DATAin%Croots_stock(dn)*DATAin%Croots_stock_unc(dn)))**2d0
       end do
       likelihood=likelihood-0.5*tot_exp
    endif

    ! Clitter log-likelihood
    ! WARNING WARNING WARNING hack in place to estimate fraction of litter pool
    ! originating from surface pools
    tot_exp = 0.
    if (DATAin%nClit_stock > 0) then
       do n = 1, DATAin%nClit_stock
         dn=DATAin%Clit_stockpts(n)
         ! note that division is the uncertainty
!         tot_exp=tot_exp+((log((sum(DATAin%M_FLUXES(:,10))/sum(DATAin%M_FLUXES(:,10)+DATAin%M_FLUXES(:,12)))*DATAin%M_POOLS(dn,5))/DATAin%Clit_stock(dn))/log(2.))**2d0
         tot_exp=tot_exp+(((sum(DATAin%M_FLUXES(:,10))/sum(DATAin%M_FLUXES(:,10)+DATAin%M_FLUXES(:,12))) &
                           *(DATAin%M_POOLS(dn,5))-DATAin%Clit_stock(dn))/(DATAin%Clit_stock(dn)*DATAin%Clit_stock_unc(dn)))**2d0
      end do
       likelihood=likelihood-0.5*tot_exp
    endif

    ! Csom log-likelihood
    tot_exp = 0.
    if (DATAin%nCsom_stock > 0) then
       do n = 1, DATAin%nCsom_stock
         dn=DATAin%Csom_stockpts(n)
         ! note that division is the uncertainty
!         tot_exp=tot_exp+(log(DATAin%M_POOLS(dn,6)/DATAin%Csom_stock(dn))/log(2.))**2.
         tot_exp=tot_exp+((DATAin%M_POOLS(dn,6)-DATAin%Csom_stock(dn))/(DATAin%Csom_stock(dn)*DATAin%Csom_stock_unc(dn)))**2d0
       end do
       likelihood=likelihood-0.5*tot_exp
    endif

    ! growth rate contraints ensures that end point pools are within uncertainty
    ! of 2 from the prior value for proportional rate of change
    ! this could / should be altered or extended to allow for time varying, pool
    ! specific values as above
    ! increments as done above for NEE,LAI and GPP
!    if (DATAin%otherpriors(1) > -9999) then
!       do n = 1, DATAin%nopools
!          pool_dynamics=DATAin%M_POOLS(DATAin%nodays+1,n)/DATAin%M_POOLS(1,n)
!          likelihood=likelihood-(0.5*(log(pool_dynamics/DATAin%otherpriors(1))/log(DATAin%otherpriors(1)))**2d0)
!       end do
!    end if

    ! check that log-likelihood is an actual number
    if (likelihood /= likelihood) then
       likelihood=log(infini)
    end if
    ! don't forget to return
    return

  end function likelihood
  !
  !------------------------------------------------------------------
  !
end module model_likelihood_module
