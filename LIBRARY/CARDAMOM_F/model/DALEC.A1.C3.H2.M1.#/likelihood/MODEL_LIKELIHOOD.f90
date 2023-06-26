
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
  public :: model_likelihood, find_edc_initial_values, &
            sub_model_likelihood, sqrt_model_likelihood, log_model_likelihood

  ! declare needed types
  type EDCDIAGNOSTICS
    integer :: nedc = 100    ! number of edcs being assessed
    integer :: PASSFAIL(100) ! allow space for 100 possible checks, dim should equal nedc
    integer :: EDC
    integer :: DIAG
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
           call MHMCMC(P_target,model_likelihood,edc_model_likelihood)

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
           else
               PEDC_prev = PEDC
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
  subroutine edc_model_likelihood(PARS, ML_obs_out, ML_prior_out)
    use cardamom_structures, only: DATAin
    use MCMCOPT, only: PI
    use CARBON_MODEL_MOD, only: carbon_model

    ! Model likelihood function specifically intended for the determination of
    ! appropriate initial parameter choices, consistent with EDCs

    implicit none

    ! declare inputs
    double precision, dimension(PI%npars), intent(inout) :: PARS
    ! output
    double precision, intent(inout) :: ML_obs_out, ML_prior_out

    ! declare local variables
    integer ::  n
    double precision :: tot_exp, ML, EDC1, EDC2, infini

    ! if == 0 EDCs are checked only until the first failure occurs
    ! if == 1 then all EDCs are checked irrespective of whether or not one has failed
    EDCD%DIAG = 1
    ML_obs_out = 0d0 ; ML_prior_out = 0d0

    ! Perform a more aggressive sanity check which compares the bulk difference
    ! in all fluxes and pools from multiple runs of the same parameter set
    if (.not.sanity_check) call model_sanity_check(PI%parini)

    ! call EDCs which can be evaluated prior to running the model
    call assess_EDC1(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)

    ! next need to run the model itself
    call CARBON_MODEL(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT,DATAin%M_LAI,DATAin%M_NEE &
                     ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%pft   &
                     ,DATAin%nopars,DATAin%nomet,DATAin%nopools   &
                     ,DATAin%nofluxes,DATAin%M_GPP                &
                     ,PI%stock_seed_labile,PI%DS_shoot,PI%DS_root &
                     ,PI%fol_frac,PI%stem_frac,PI%root_frac,PI%DS_LRLV&
                     ,PI%LRLV,PI%DS_LRRT,PI%LRRT)

    ! assess post running EDCs
    call assess_EDC2(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools &
             ,DATAin%nodays,DATAin%deltat,PI%parmax,PARS,DATAin%MET &
             ,DATAin%M_LAI,DATAin%M_NEE,DATAin%M_GPP,DATAin%M_POOLS &
             ,DATAin%M_FLUXES,DATAin%meantemp,EDC2)

    ! calculate the likelihood
    tot_exp = sum(1d0-EDCD%PASSFAIL(1:EDCD%nedc))
!    tot_exp = 0d0
!    do n = 1, EDCD%nedc
!       tot_exp=tot_exp+(1d0-EDCD%PASSFAIL(n))
!       if (EDCD%PASSFAIL(n) /= 1) print*,"failed edcs are: ", n
!    end do ! checking EDCs
!    ! for testing purposes, stop the model when start achieved
!    if (sum(EDCD%PASSFAIL) == 100) then
!        print*,"Found it!" ; stop
!    endif

    ! convert to a probability
    ML_obs_out = -5d0*tot_exp*DATAin%EDC

  end subroutine edc_model_likelihood
  !
  !------------------------------------------------------------------
  !
  subroutine sub_model_likelihood(PARS,ML_obs_out,ML_prior_out)
    use MCMCOPT, only:  PI
    use CARBON_MODEL_MOD, only: carbon_model
    use cardamom_structures, only: DATAin

    ! this subroutine is responsible, under normal circumstances for the running
    ! of the DALEC model, calculation of the log-likelihood for comparison
    ! assessment of parameter performance and use of the EDCs if they are
    ! present / selected

    implicit none

    ! declare inputs
    double precision, dimension(PI%npars), intent(inout) :: PARS ! current parameter vector
    ! output
    double precision, intent(inout) :: ML_obs_out, &  ! observation + EDC log-likelihood
                                       ML_prior_out   ! prior log-likelihood
    ! declare local variables
    double precision :: EDC1, EDC2

    ! initial values
    ML_obs_out = 0d0 ; ML_prior_out = 0d0 ; EDC1 = 1d0 ; EDC2 = 1d0
    ! if == 0 EDCs are checked only until the first failure occurs
    ! if == 1 then all EDCs are checked irrespective of whether or not one has failed
    EDCD%DIAG = 0

    if (DATAin%EDC == 1) then
        ! call EDCs which can be evaluated prior to running the model
        call assess_EDC1(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)
        ! update the likelihood score based on EDCs driving total rejection
        ! proposed parameters
        ML_obs_out = log(EDC1)
    endif !

    ! then this is a crop run....
    ! run the dalec model
    call CARBON_MODEL(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT,DATAin%M_LAI,DATAin%M_NEE &
                     ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%pft   &
                     ,DATAin%nopars,DATAin%nomet,DATAin%nopools   &
                     ,DATAin%nofluxes,DATAin%M_GPP                &
                     ,PI%stock_seed_labile,PI%DS_shoot,PI%DS_root &
                     ,PI%fol_frac,PI%stem_frac,PI%root_frac,PI%DS_LRLV&
                     ,PI%LRLV,PI%DS_LRRT,PI%LRRT)


    ! if first set of EDCs have been passed, move on to the second
    if (DATAin%EDC == 1) then
        ! check edc2
        call assess_EDC2(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools &
                        ,DATAin%nodays,DATAin%deltat,PI%parmax,PARS,DATAin%MET &
                        ,DATAin%M_LAI,DATAin%M_NEE,DATAin%M_GPP,DATAin%M_POOLS &
                        ,DATAin%M_FLUXES,DATAin%meantemp,EDC2)

        ! Add EDC2 log-likelihood to absolute accept reject...
        ML_obs_out = ML_obs_out + log(EDC2)
    end if ! DATAin%EDC == 1

    ! Calculate log-likelihood associated with priors
    ! We always want this
    ML_prior_out = likelihood_p(PI%npars,DATAin%parpriors,DATAin%parpriorunc,DATAin%parpriorweight,PARS)
    ! calculate final model likelihood when compared to obs
    ML_obs_out = ML_obs_out + scale_likelihood(PI%npars,PARS)

  end subroutine sub_model_likelihood
  !
  !------------------------------------------------------------------
  !
  subroutine sqrt_model_likelihood(PARS,ML_obs_out,ML_prior_out)
    use MCMCOPT, only:  PI
    use CARBON_MODEL_MOD, only: carbon_model
    use cardamom_structures, only: DATAin

    ! this subroutine is responsible, under normal circumstances for the running
    ! of the DALEC model, calculation of the log-likelihood for comparison
    ! assessment of parameter performance and use of the EDCs if they are
    ! present / selected

    implicit none

    ! declare inputs
    double precision, dimension(PI%npars), intent(inout) :: PARS ! current parameter vector
    ! output
    double precision, intent(inout) :: ML_obs_out, &  ! observation + EDC log-likelihood
                                       ML_prior_out   ! prior log-likelihood
    ! declare local variables
    double precision :: EDC1, EDC2

    ! initial values
    ML_obs_out = 0d0 ; ML_prior_out = 0d0 ; EDC1 = 1d0 ; EDC2 = 1d0
    ! if == 0 EDCs are checked only until the first failure occurs
    ! if == 1 then all EDCs are checked irrespective of whether or not one has failed
    EDCD%DIAG = 0

    if (DATAin%EDC == 1) then
        ! call EDCs which can be evaluated prior to running the model
        call assess_EDC1(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)
        ! update the likelihood score based on EDCs driving total rejection
        ! proposed parameters
        ML_obs_out = log(EDC1)
    endif !

    ! then this is a crop run....
    ! run the dalec model
    call CARBON_MODEL(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT,DATAin%M_LAI,DATAin%M_NEE &
                     ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%pft   &
                     ,DATAin%nopars,DATAin%nomet,DATAin%nopools   &
                     ,DATAin%nofluxes,DATAin%M_GPP                &
                     ,PI%stock_seed_labile,PI%DS_shoot,PI%DS_root &
                     ,PI%fol_frac,PI%stem_frac,PI%root_frac,PI%DS_LRLV&
                     ,PI%LRLV,PI%DS_LRRT,PI%LRRT)


    ! if first set of EDCs have been passed, move on to the second
    if (DATAin%EDC == 1) then
        ! check edc2
        call assess_EDC2(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools &
                        ,DATAin%nodays,DATAin%deltat,PI%parmax,PARS,DATAin%MET &
                        ,DATAin%M_LAI,DATAin%M_NEE,DATAin%M_GPP,DATAin%M_POOLS &
                        ,DATAin%M_FLUXES,DATAin%meantemp,EDC2)

        ! Add EDC2 log-likelihood to absolute accept reject...
        ML_obs_out = ML_obs_out + log(EDC2)
    end if ! DATAin%EDC == 1

    ! Calculate log-likelihood associated with priors
    ! We always want this
    ML_prior_out = likelihood_p(PI%npars,DATAin%parpriors,DATAin%parpriorunc,DATAin%parpriorweight,PARS)
    ! calculate final model likelihood when compared to obs
    ML_obs_out = ML_obs_out + sqrt_scale_likelihood(PI%npars,PARS)

  end subroutine sqrt_model_likelihood
  !
  !------------------------------------------------------------------
  !
  subroutine log_model_likelihood(PARS,ML_obs_out,ML_prior_out)
    use MCMCOPT, only:  PI
    use CARBON_MODEL_MOD, only: carbon_model
    use cardamom_structures, only: DATAin

    ! this subroutine is responsible, under normal circumstances for the running
    ! of the DALEC model, calculation of the log-likelihood for comparison
    ! assessment of parameter performance and use of the EDCs if they are
    ! present / selected

    implicit none

    ! declare inputs
    double precision, dimension(PI%npars), intent(inout) :: PARS ! current parameter vector
    ! output
    double precision, intent(inout) :: ML_obs_out, &  ! observation + EDC log-likelihood
                                       ML_prior_out   ! prior log-likelihood
    ! declare local variables
    double precision :: EDC1, EDC2

    ! initial values
    ML_obs_out = 0d0 ; ML_prior_out = 0d0 ; EDC1 = 1d0 ; EDC2 = 1d0
    ! if == 0 EDCs are checked only until the first failure occurs
    ! if == 1 then all EDCs are checked irrespective of whether or not one has failed
    EDCD%DIAG = 0

    if (DATAin%EDC == 1) then
        ! call EDCs which can be evaluated prior to running the model
        call assess_EDC1(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)
        ! update the likelihood score based on EDCs driving total rejection
        ! proposed parameters
        ML_obs_out = log(EDC1)
    endif !

    ! then this is a crop run....
    ! run the dalec model
    call CARBON_MODEL(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT,DATAin%M_LAI,DATAin%M_NEE &
                     ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%pft   &
                     ,DATAin%nopars,DATAin%nomet,DATAin%nopools   &
                     ,DATAin%nofluxes,DATAin%M_GPP                &
                     ,PI%stock_seed_labile,PI%DS_shoot,PI%DS_root &
                     ,PI%fol_frac,PI%stem_frac,PI%root_frac,PI%DS_LRLV&
                     ,PI%LRLV,PI%DS_LRRT,PI%LRRT)


    ! if first set of EDCs have been passed, move on to the second
    if (DATAin%EDC == 1) then
        ! check edc2
        call assess_EDC2(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools &
                        ,DATAin%nodays,DATAin%deltat,PI%parmax,PARS,DATAin%MET &
                        ,DATAin%M_LAI,DATAin%M_NEE,DATAin%M_GPP,DATAin%M_POOLS &
                        ,DATAin%M_FLUXES,DATAin%meantemp,EDC2)

        ! Add EDC2 log-likelihood to absolute accept reject...
        ML_obs_out = ML_obs_out + log(EDC2)
    end if ! DATAin%EDC == 1

    ! Calculate log-likelihood associated with priors
    ! We always want this
    ML_prior_out = likelihood_p(PI%npars,DATAin%parpriors,DATAin%parpriorunc,DATAin%parpriorweight,PARS)
    ! calculate final model likelihood when compared to obs
    ML_obs_out = ML_obs_out + log_scale_likelihood(PI%npars,PARS)

  end subroutine log_model_likelihood
  !
  !------------------------------------------------------------------
  !
  subroutine model_sanity_check(PARS)
    use cardamom_structures, only: DATAin
    use MCMCOPT, only: PI
    use CARBON_MODEL_MOD, only: carbon_model

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

    ! next need to run the model itself
    call CARBON_MODEL(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT,DATAin%M_LAI,DATAin%M_NEE &
                     ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%pft   &
                     ,DATAin%nopars,DATAin%nomet,DATAin%nopools   &
                     ,DATAin%nofluxes,DATAin%M_GPP                &
                     ,PI%stock_seed_labile,PI%DS_shoot,PI%DS_root &
                     ,PI%fol_frac,PI%stem_frac,PI%root_frac,PI%DS_LRLV&
                     ,PI%LRLV,PI%DS_LRRT,PI%LRRT)

    call CARBON_MODEL(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT,DATAin%M_LAI,DATAin%M_NEE &
                     ,local_fluxes,local_pools,DATAin%pft   &
                     ,DATAin%nopars,DATAin%nomet,DATAin%nopools   &
                     ,DATAin%nofluxes,DATAin%M_GPP                &
                     ,PI%stock_seed_labile,PI%DS_shoot,PI%DS_root &
                     ,PI%fol_frac,PI%stem_frac,PI%root_frac,PI%DS_LRLV&
                     ,PI%LRLV,PI%DS_LRRT,PI%LRRT)

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
    print*,"Model Sanity Check Completed"

  end subroutine model_sanity_check
  !
  !------------------------------------------------------------------
  !
  subroutine assess_EDC1(PARS, npars, meantemp, meanrad, EDC1)

    ! the first of two subroutine to assess current parameters for passing
    ! realism tests for crop
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
    double precision :: torfol,tmp ! yearly leaf loss fraction

    ! set initial value
    EDC1 = 1d0
    DIAG = EDCD%DIAG

    ! set all EDCs to 1 (pass)
    EDCD%PASSFAIL(1:EDCD%nedc) = 1

    !
    ! begin checking EDCs
    !

    ! Turnover of litter faster than turnover of som
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(10) > pars(9))) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(1) = 0
    endif

    ! decomposition of litter to SOM greater than SOM to air
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(10) > pars(1))) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(2) = 0
    endif

    ! turnover of foliage faster than turnover of wood
! TLS: turnover off because foliage and stem turnovers are made same
!    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(6) > pars(5)) then
!       EDC1 = 0d0 ; EDCD%PASSFAIL(3) = 0
!    end if

    ! pre_DR should be greater than post_DR
!    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(4) > pars(3))) then
!        EDC1 = 0d0 ; EDCD%PASSFAIL(4) = 0
!    endif

    ! for development: Tmin should be < topt and topt should be < tmax
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(26) > pars(28) &
                                     .or. pars(28) > pars(27) &
                                     .or. pars(26) > pars(27))) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(5) = 0
    endif

    ! for development: the difference between each Tmin,Topt,Tmax > 1.
    if ((EDC1 == 1 .or. DIAG == 1) .and. (abs(pars(26)-pars(28)) < 1d0 &
                                     .or. abs(pars(28)-pars(27)) < 1d0  &
                                     .or. abs(pars(26)-pars(27)) < 1d0)) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(6) = 0
    endif

   ! for vernalisation: Tmin < Topt < Tmax
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(29) > pars(31) &
                                     .or. pars(31) > pars(30) &
                                     .or. pars(29) > pars(30))) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(7) = 0
    endif

   ! for vernalisation: the difference between each Tmin, Topt, Tmax
    if ((EDC1 == 1 .or. DIAG == 1) .and. ( abs(pars(29)-pars(31)) < 1d0 &
                                      .or. abs(pars(31)-pars(30)) < 1d0 &
                                      .or. abs(pars(29)-pars(30)) < 1d0 ) ) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(8) = 0
    endif

   ! development temperature value should be larger corresponding vernalisation
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(29) > pars(26) &
                                     .or. pars(31) > pars(28) &
                                     .or. pars(30) > pars(27))) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(9) = 0
    endif

!    ! plough must be before sow and after harvest: WINTER ONLY
!    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(16) > pars(12) ) then
!        EDC1 = 0d0 ; EDCD%PASSFAIL(10) = 0
!    endif
!
!    ! harvest cannot be more than 345 after harvest
!    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(15) < pars(12)+345.25d0) ) then
!        EDC1 = 0d0 ; EDCD%PASSFAIL(11) = 0
!    endif
!
!    ! plough must be before sow and after harvest: WINTER ONLY
!    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(16) < pars(15) .or. &
!        pars(16) > pars(12)) ) then
!        EDC1 = 0d0 ; EDCD%PASSFAIL(12) = 0
!    endif

    ! CN ratio of leaf should also be between 95CI of trait database values
    ! Kattge et al (2011)
    tmp = (pars(17)/(10d0**pars(11)))
    if ((EDC1 == 1 .or. DIAG == 1) .and. (tmp > 43.76895d0 .or. tmp < 10.82105d0)) then
       EDC1=0 ; EDCD%PASSFAIL(13) = 0
    endif

    ! could and probably should add some more
  end subroutine assess_EDC1
  !
  !------------------------------------------------------------------
  !
  subroutine assess_EDC2(npars,nomet,nofluxes,nopools,nodays,deltat &
                        ,parmax,pars,met,M_LAI,M_NEE,M_GPP,M_POOLS,M_FLUXES &
                        ,meantemp,EDC2)

    use cardamom_structures, only: DATAin
    use CARBON_MODEL_MOD, only: resp_rate_temp_coeff,ts_length,linear_model_gradient

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
    integer :: n, DIAG, no_years, nn
    double precision :: mean_pools(nopools), decay_coef, meangpp, EQF, PEDC, infi
    double precision :: fauto & ! Fractions of GPP to autotrophic respiration
             ,ffol  & ! Fraction of GPP to foliage
             ,flab  & ! Fraction of GPP to labile pool
             ,froot & ! Fraction of GPP to root
             ,flit  & !
             ,fwood & ! Fraction of GPP to wood
             ,fsom    ! fraction of GPP som under eqilibrium conditions

    ! set initial value
    fauto = sum(M_FLUXES(:,3)) / sum(M_FLUXES(:,1))
    ffol = sum(M_FLUXES(:,4)) / (sum(M_FLUXES(:,1))*fauto)
    flab = sum(M_FLUXES(:,5)) / (sum(M_FLUXES(:,1))*fauto)
    froot = sum(M_FLUXES(:,6)) / (sum(M_FLUXES(:,1))*fauto)
    fwood = sum(M_FLUXES(:,7)) / (sum(M_FLUXES(:,1))*fauto)
    fsom = fwood+(froot+flab+ffol)*pars(1)/(pars(1)+pars(10))
    flit = (froot+flab+ffol)
    ! length of time step in hours..
    ts_length = 24d0
    ! initial value
    infi = 0d0
    ! update initial values
    DIAG = EDCD%DIAG
    ! give EDC2 an initial value
    EDC2 = 1

    ! SOM attractor - must be within a factor of 2 from Csom0
    ! eqiulibrium factor (in comparison with initial conditions)
    EQF = 10d0

    ! initialise and then calculate mean gpp values
    meangpp = sum(M_GPP(1:nodays))/dble(nodays)

    ! EDC 11 - SOM steady state within order magnitude of initial conditions
    if ((EDC2 == 1 .or. DIAG == 1) .and. &
       ((meangpp*fsom)/(pars(10)*0.5d0*exp(resp_rate_temp_coeff*meantemp))) > (pars(23)*EQF)) then
       EDC2 = 0d0 ; EDCD%PASSFAIL(14) = 0
    end if
    if ((EDC2 == 1 .or. DIAG == 1) .and. &
       ((meangpp*fsom)/(pars(10)*0.5d0*exp(resp_rate_temp_coeff*meantemp))) < (pars(23)/EQF)) then
       EDC2 = 0d0 ; EDCD%PASSFAIL(15) = 0
    endif

    ! EDC 12 - Litter steady state assumptions
    if ((EDC2 == 1 .or. DIAG == 1) .and. &
       ((meangpp*flit)/(pars(9)*0.5d0*exp(resp_rate_temp_coeff*meantemp))) > (pars(22)*EQF)) then
        EDC2 = 0d0 ; EDCD%PASSFAIL(16) = 0
    endif
    if ((EDC2 == 1 .or. DIAG == 1) .and. &
       ((meangpp*flit)/(pars(9)*0.5d0*exp(resp_rate_temp_coeff*meantemp))) < (pars(22)/EQF)) then
        EDC2 = 0d0 ; EDCD%PASSFAIL(17) = 0
    endif

    ! EDC 13
    ! assesses the exponential decay/growth of the Csom pool

    !  work out how many completed years there are in the system
    no_years=int(nint(sum(deltat)/365.25d0))

    ! only do this for the Csom pool
    do n = 1, 1 !nopools
       if (EDC2 == 1 .or. DIAG == 1) then
          decay_coef=expdecay2(M_POOLS(1:(nodays+1),6),deltat,nodays+1)
          ! next assess the decay coefficient for meetings the EDC criterion
          if (abs(-log(2d0)/decay_coef) < (365.25d0*dble(no_years)) .and. decay_coef < 0d0 ) then
             EDC2 = 0d0 ; EDCD%PASSFAIL(18) = 0
          end if ! EDC conditions
       end if ! EDC .or. DIAG condition
    end do ! pools loop

    ! EDC 14
    ! assesses the exponential decay/growth of the Clit pool

    ! only do this for the Clit pool
    do n = 1, 1 !nopools
       if (EDC2 == 1 .or. DIAG == 1) then
          decay_coef=expdecay2(M_POOLS(1:(nodays+1),5),deltat,nodays+1)
          ! next assess the decay coefficient for meetings the EDC criterion
          if (abs(-log(2d0)/decay_coef) < (365.25d0*dble(no_years)) .and. decay_coef < 0d0 ) then
             EDC2 = 0d0 ; EDCD%PASSFAIL(19) = 0
          end if ! EDC conditions
       end if ! EDC .or. DIAG condition
    end do ! pools loop

    ! we know that the crop model should produce some yield - therefore we
    ! reject parameter sets which generate no yield ever!
    if ((EDC2 == 1 .or. DIAG == 1) .and. sum(M_FLUXES(1:nodays,21)) < (1d0*dble(no_years)) ) then
        EDC2 = 0d0 ; EDCD%PASSFAIL(20) = 0
    endif

    ! LAI time series linear model must retrieve gradient which is at least
    ! positive (or some other reasonable critical threshold)
    ! if ((EDC2 == 1 .or. DIAG == 1) .and. &
    !     linear_model_gradient(DATAin%M_LAI(DATAin%laipts),DATAin%LAI(DATAin%laipts),DATAin%nlai) < 0d0 ) then
    !     EDC2 = 0d0 ; EDCD%PASSFAIL(21) = 0
    ! endif

    ! Function to calculate the gradient of a linear model for a given depentent
    ! variable (y) based on predictive variable (x). The typical use of this
    ! function will in fact be to assume that x is time.

    !
    ! EDCs done, below are additional fault detection conditions
    !

    ! additional faults can be stored in locations 35 - 40 of the PASSFAIL array

    ! ensure minimum pool values are >= 0 and /= NaN
    if (EDC2 == 1 .or. DIAG == 1) then

      do n = 1, nopools
         if (minval(M_POOLS(1:nodays,n)) < 0d0 .or. maxval(abs(M_POOLS(1:nodays,n))) == abs(log(infi)) .or. &
             minval(M_POOLS(1:nodays,n)) /= minval(M_POOLS(1:nodays,n))) then
             EDC2 = 0d0 ; EDCD%PASSFAIL(55+n) = 0
         endif
      end do

      do n = 1, nofluxes
         if (maxval(abs(M_FLUXES(1:nodays,n))) == abs(log(infi)) .or. &
            minval(M_FLUXES(1:nodays,n)) /= minval(M_FLUXES(1:nodays,n))) then
             EDC2 = 0d0 ; EDCD%PASSFAIL(55+nopools+n) = 0
         endif
      end do

    end if ! min pool assessment

  end subroutine assess_EDC2
  !
  !------------------------------------------------------------------
  !
  double precision function cal_mean_annual_pools(pools,year,interval,averaging_period)

    ! Function calculates the mean model pools values for each individual year
    ! in the simulation

    implicit none

    ! declare input variables
    integer, intent(in) :: year           & ! which year are we working on
                          ,averaging_period ! number of days in analysis period

    double precision, intent(in) :: pools(averaging_period) & ! input pool state variables
                                 ,interval((averaging_period-1))      ! model time step in decimal days

    ! declare local variables
    integer :: startday, endday

    ! calculate some constants
    startday = floor(365.25d0*dble(year-1)/(sum(interval)/dble(averaging_period-1)))+1
    endday = floor(365.25d0*dble(year)/(sum(interval)/dble(averaging_period-1)))

    ! pool through and work out the annual mean values
    cal_mean_annual_pools = sum(pools(startday:endday))/dble(endday-startday)

    ! ensure function returns
    return

  end function cal_mean_annual_pools
  !
  !------------------------------------------------------------------
  !
  double precision function cal_max_annual_pools(pools,year,interval,averaging_period)

    ! Function calculates the max model pools values for each individual year
    ! in the simulation

    implicit none

    ! declare input variables
    integer, intent(in) :: year            & ! which year are we working on
                          ,averaging_period  ! number of days in analysis period

    double precision, intent(in) :: pools(averaging_period) & ! input pool state variables
                                 ,interval((averaging_period-1))      ! model time step in decimal days

    ! declare local variables
    integer :: startday, endday

    ! calculate some constants
    startday = floor(365.25d0*dble(year-1)/(sum(interval)/dble(averaging_period-1)))+1
    endday = floor(365.25d0*dble(year)/(sum(interval)/dble(averaging_period-1)))

    ! pool through and work out the annual max values
    cal_max_annual_pools = maxval(pools(startday:endday))

    ! ensure function returns
    return

  end function cal_max_annual_pools
  !
  !------------------------------------------------------------------
  !
  double precision function expdecay2(pools,interval,averaging_period)

   ! Function to calculate the exponential decay coefficients used several EDCs.
   ! We assumpe the equation Cexp= a + b*exp(c*t)

   implicit none

   ! declare input variables
   integer, intent(in) :: averaging_period ! i.e. nodays + 1

   double precision, intent(in) :: pools(averaging_period) & ! input pool state variables
                                  ,interval((averaging_period-1))      ! model time step in decimal days

   ! declare local variables
   integer :: n, aw_int
   integer, parameter :: os = 1 ! offset days
   double precision :: aw, aw_1 &
                      ,MP0   & ! mean pool (year 1 to year end-2)
                      ,MP1   & ! mean pool (year 2 to year end-1)
                      ,MP0os & ! mean pool (year 1+os to year end-2+os)
                      ,MP1os & ! mean pool (year 2+os to year end-2+os)
                      ,dcdt1 & ! gradient of exponential over time in second year
                      ,dcdt0   ! gradient of exponential over time in first year

   ! declare initial values / constants
   aw = floor(365.25d0/(sum(interval)/dble(averaging_period-1))) ! averaging window
   aw_1 = aw ** (-1d0) ; aw_int = int(aw)
   MP0 = 0d0 ; MP1 = 0d0 ; MP0os = 0d0 ; MP1os = 0d0

   ! estimate mean stock for first year
   MP0 = sum(pools(1:aw_int))
   MP0 = MP0*aw_1

   ! estimate mean stock for second year
   MP1 = sum(pools((aw_int+1):(aw_int*2)))
   MP1 = MP1*aw_1

   ! estimate mean stock for first year with offset
   MP0os = sum(pools((1+os):(aw_int+os)))
   MP0os = MP0os*aw_1

   ! estimate mean stock for second year with offset
   MP1os = sum(pools((aw_int+os+1):((aw_int*2)+os)))
   MP1os = MP1os*aw_1

   ! derive mean gradient ratio (dcdt1/dcdt0)
   ! where dcdt1 is the numeric gradient between n+1 and n+365+1
   ! and dcdt0 os the numeric gradient between n and n+365
   dcdt1 = MP1os-MP0os
   dcdt0 = MP1-MP0

   ! using multiple year mean to determine c
   if ((dcdt1 > 0d0 .and. dcdt0 < 0d0) .or. (dcdt1 < 0d0 .and. dcdt0 > 0d0) &
       .or. dcdt1 == 0d0 .or. dcdt0 == 0d0) then
       ! then return error values
       expdecay2 = 1d0
   else
       expdecay2 = log(dcdt1/dcdt0) / (dble(os)*(sum(interval)/dble(averaging_period-1)))
   end if

   ! ensure return
   return

  end function expdecay2
  !
  !------------------------------------------------------------------
  !
  subroutine model_likelihood(PARS,ML_obs_out,ML_prior_out)
    use MCMCOPT, only:  PI
    use CARBON_MODEL_MOD, only: carbon_model
    use cardamom_structures, only: DATAin

    ! this subroutine is responsible, under normal circumstances for the running
    ! of the DALEC model, calculation of the log-likelihood for comparison
    ! assessment of parameter performance and use of the EDCs if they are
    ! present / selected

    implicit none

    ! declare inputs
    double precision, dimension(PI%npars), intent(inout) :: PARS ! current parameter vector
    ! output
    double precision, intent(inout) :: ML_obs_out, &  ! observation + EDC log-likelihood
                                       ML_prior_out   ! prior log-likelihood
    ! declare local variables
    double precision :: EDC1, EDC2

    ! initial values
    ML_obs_out = 0d0 ; ML_prior_out = 0d0 ; EDC1 = 1d0 ; EDC2 = 1d0
    ! if == 0 EDCs are checked only until the first failure occurs
    ! if == 1 then all EDCs are checked irrespective of whether or not one has failed
    EDCD%DIAG = 0

    if (DATAin%EDC == 1) then
        ! call EDCs which can be evaluated prior to running the model
        call assess_EDC1(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)

        ! update the likelihood score based on EDCs driving total rejection
        ! proposed parameters
        ML_obs_out = log(EDC1)
    endif !

    ! then this is a crop run....
    ! run the dalec model
    call CARBON_MODEL(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT,DATAin%M_LAI,DATAin%M_NEE &
                     ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%pft   &
                     ,DATAin%nopars,DATAin%nomet,DATAin%nopools   &
                     ,DATAin%nofluxes,DATAin%M_GPP                &
                     ,PI%stock_seed_labile,PI%DS_shoot,PI%DS_root &
                     ,PI%fol_frac,PI%stem_frac,PI%root_frac,PI%DS_LRLV&
                     ,PI%LRLV,PI%DS_LRRT,PI%LRRT)

    ! if first set of EDCs have been passed, move on to the second
    if (DATAin%EDC == 1) then
        ! check edc2
        call assess_EDC2(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools &
                        ,DATAin%nodays,DATAin%deltat,PI%parmax,PARS,DATAin%MET &
                        ,DATAin%M_LAI,DATAin%M_NEE,DATAin%M_GPP,DATAin%M_POOLS &
                        ,DATAin%M_FLUXES,DATAin%meantemp,EDC2)

        ! Add EDC2 log-likelihood to absolute accept reject...
        ML_obs_out = ML_obs_out + log(EDC2)
    end if ! DATAin%EDC == 1

    ! Calculate log-likelihood associated with priors
    ! We always want this
    ML_prior_out = likelihood_p(PI%npars,DATAin%parpriors,DATAin%parpriorunc,DATAin%parpriorweight,PARS)
    ! calculate final model likelihood when compared to obs
    ML_obs_out = ML_obs_out + likelihood(PI%npars,PARS)

  end subroutine model_likelihood
  !
  !------------------------------------------------------------------
  !
  double precision function likelihood_p(npars,parpriors,parpriorunc,parpriorweight,pars)
    ! function calculates the parameter based log-likelihood for the current set
    ! of parameters. This assumes that we have any actual priors / prior
    ! uncertainties to be working with. This does include initial states, as we
    ! consider them to be parameters

    implicit none

    ! declare input variables
    integer, intent(in) :: npars
    double precision, dimension(npars), intent(in) :: pars         & ! current parameter vector
                                                     ,parpriors    & ! prior values for parameters
                                                     ,parpriorunc  & ! prior uncertainties
                                                     ,parpriorweight ! prior weighting

    ! declare local variables
    integer :: n
    double precision, dimension(npars) :: local_likelihood

    ! set initial value
    likelihood_p = 0d0 ; local_likelihood = 0d0

    ! now loop through defined parameters for their uncertainties
    where (parpriors > -9999) local_likelihood = parpriorweight*((pars-parpriors)/parpriorunc)**2
    likelihood_p = sum(local_likelihood) * (-0.5d0)

    ! dont for get to return
    return

  end function likelihood_p
  !
  !------------------------------------------------------------------
  !
  double precision function likelihood(npars,pars)
    use cardamom_structures, only: DATAin
    use carbon_model_mod, only: layer_thickness

    ! calculates the likelihood of of the model output compared to the available
    ! observations which have been input to the model

    implicit none

    ! declare arguments
    integer, intent(in) :: npars
    double precision, dimension(npars), intent(in) :: pars

    ! declare local variables
    integer :: n, dn, no_years, y, s
    double precision :: tot_exp, tmp_var, infini, input, output
    double precision, dimension(DATAin%nodays) :: mid_state
    double precision, allocatable :: mean_annual_pools(:)

    ! initial value
    likelihood = 0d0 ; infini = 0d0 ; mid_state = 0d0

! Currently no distinction between NBE and NEE due to lack of fire estimates
!    ! NBE Log-likelihood
!    if (DATAin%nnbe > 0) then
!       tot_exp = sum((((DATAin%M_NEE(DATAin%nbepts(1:DATAin%nnbe))+DATAin%M_FLUXES(DATAin%nbepts(1:DATAin%nnbe),17)) &
!                       -DATAin%NBE(DATAin%nbepts(1:DATAin%nnbe))) &
!                       /DATAin%NBE_unc(DATAin%nbepts(1:DATAin%nnbe)))**2)
!       likelihood = likelihood-tot_exp
!    endif

    ! GPP Log-likelihood
    if (DATAin%ngpp > 0) then
       tot_exp = sum(((DATAin%M_GPP(DATAin%gpppts(1:DATAin%ngpp))-DATAin%GPP(DATAin%gpppts(1:DATAin%ngpp))) &
                       /DATAin%GPP_unc(DATAin%gpppts(1:DATAin%ngpp)))**2)
       likelihood = likelihood-tot_exp
    endif

    ! Evap Log-likelihood
    if (DATAin%nEvap > 0) then
       tot_exp = sum(((DATAin%M_FLUXES(DATAin%Evappts(1:DATAin%nEvap),19)-DATAin%Evap(DATAin%Evappts(1:DATAin%nEvap))) &
                       /DATAin%Evap_unc(DATAin%Evappts(1:DATAin%nEvap)))**2)
       likelihood = likelihood-tot_exp
    endif

! FIRE not currently coded into DALEC_CROP_BUCKET
!    ! Fire Log-likelihood
!    if (DATAin%nFire > 0) then
!       tot_exp = sum(((DATAin%M_FLUXES(DATAin%Firepts(1:DATAin%nFire),17)-DATAin%Fire(DATAin%Firepts(1:DATAin%nFire))) &
!                       /DATAin%Fire_unc(DATAin%Firepts(1:DATAin%nFire)))**2)
!       likelihood = likelihood-tot_exp
!    endif

    ! Assume physical property is best represented as the mean of value at beginning and end of times step
    if (DATAin%nlai > 0) then
       ! Create vector of (LAI_t0 + LAI_t1) * 0.5, note / pars(17) to convert foliage C to LAI
       mid_state = ( ( DATAin%M_POOLS(1:DATAin%nodays,2) + DATAin%M_POOLS(2:(DATAin%nodays+1),2) ) &
                 * 0.5d0 ) / pars(17)
       ! Split loop to allow vectorisation
       tot_exp = sum(((mid_state(DATAin%laipts(1:DATAin%nlai))-DATAin%LAI(DATAin%laipts(1:DATAin%nlai))) &
                       /DATAin%LAI_unc(DATAin%laipts(1:DATAin%nlai)))**2)
       ! loop split to allow vectorisation
       !tot_exp = sum(((DATAin%M_LAI(DATAin%laipts(1:DATAin%nlai))-DATAin%LAI(DATAin%laipts(1:DATAin%nlai))) &
       !                /DATAin%LAI_unc(DATAin%laipts(1:DATAin%nlai)))**2)
       do n = 1, DATAin%nlai
         dn = DATAin%laipts(n)
         ! if zero or greater allow calculation with min condition to prevent
         ! errors of zero LAI which occur in managed systems
         if (mid_state(dn) < 0d0) then
             ! if not then we have unrealistic negative values or NaN so indue
             ! error
             tot_exp = tot_exp+(-log(infini))
         endif
       end do
       likelihood = likelihood-tot_exp
    endif

    ! NEE likelihood
    if (DATAin%nnee > 0) then
       tot_exp = sum(((DATAin%M_NEE(DATAin%neepts(1:DATAin%nnee))-DATAin%NEE(DATAin%neepts(1:DATAin%nnee))) &
                       /DATAin%NEE_unc(DATAin%neepts(1:DATAin%nnee)))**2)
       likelihood = likelihood-tot_exp
    endif

    ! Reco likelihood
    if (DATAin%nreco > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nreco
         dn = DATAin%recopts(n)
         tmp_var = DATAin%M_NEE(dn)+DATAin%M_GPP(dn)
         ! note that we calculate the Ecosystem resp from GPP and NEE
         tot_exp = tot_exp+((tmp_var-DATAin%Reco(dn))/DATAin%Reco_unc(dn))**2
       end do
       likelihood = likelihood-tot_exp
    endif

    ! Cwood increment log-likelihood
    if (DATAin%nCwood_inc > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCwood_inc
         dn = DATAin%Cwood_incpts(n)
         s = max(0,dn-nint(DATAin%Cwood_inc_lag(dn)))+1
         ! Estimate the mean allocation to wood over the lag period
         tmp_var = sum(DATAin%M_FLUXES(s:dn,7)) / DATAin%Cwood_inc_lag(dn)
         tot_exp = tot_exp+((tmp_var-DATAin%Cwood_inc(dn)) / DATAin%Cwood_inc_unc(dn))**2
       end do
       likelihood = likelihood-tot_exp
    endif

    ! Cwood mortality log-likelihood
    if (DATAin%nCwood_mortality > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCwood_mortality
         dn = DATAin%Cwood_mortalitypts(n)
         s = max(0,dn-nint(DATAin%Cwood_mortality_lag(dn)))+1
         ! Estimate the mean allocation to wood over the lag period
         tmp_var = sum(DATAin%M_FLUXES(s:dn,11)) / DATAin%Cwood_mortality_lag(dn)
         tot_exp = tot_exp+((tmp_var-DATAin%Cwood_mortality(dn)) / DATAin%Cwood_mortality_unc(dn))**2
       end do
       likelihood = likelihood-tot_exp
    endif

    ! Cfoliage log-likelihood
    if (DATAin%nCfol_stock > 0) then
       ! Create vector of (FOL_t0 + FOL_t1) * 0.5
       mid_state = ( DATAin%M_POOLS(1:DATAin%nodays,2) + DATAin%M_POOLS(2:(DATAin%nodays+1),2) ) &
                 * 0.5d0
       ! Vectorised version of loop to estimate cost function
       tot_exp = sum(( (mid_state(DATAin%Cfol_stockpts(1:DATAin%nCfol_stock)) &
                       -DATAin%Cfol_stock(DATAin%Cfol_stockpts(1:DATAin%nCfol_stock)))&
                     / DATAin%Cfol_stock_unc(DATAin%Cfol_stockpts(1:DATAin%nCfol_stock)))**2)
       ! Sum with current likelihood score
       likelihood = likelihood-tot_exp
    endif

    ! Annual foliar maximum
    if (DATAin%nCfolmax_stock > 0) then
       tot_exp = 0d0
       no_years = int(nint(sum(DATAin%deltat)/365.25d0))
       if (allocated(mean_annual_pools)) deallocate(mean_annual_pools)
       allocate(mean_annual_pools(no_years))
       ! determine the annual max for each pool
       do y = 1, no_years
          ! derive mean annual foliar pool
          mean_annual_pools(y) = cal_max_annual_pools(DATAin%M_POOLS(1:(DATAin%nodays+1),2),y,DATAin%deltat,DATAin%nodays+1)
       end do ! year loop
       ! loop through the observations then
       do n = 1, DATAin%nCfolmax_stock
         ! load the observation position in stream
         dn = DATAin%Cfolmax_stockpts(n)
         ! determine which years this in in for the simulation
         y = ceiling( (dble(dn)*(sum(DATAin%deltat)/(DATAin%nodays))) / 365.25d0 )
         ! load the correct year into the analysis
         tmp_var = mean_annual_pools(y)
         ! note that division is the uncertainty
         tot_exp = tot_exp+((tmp_var-DATAin%Cfolmax_stock(dn)) / DATAin%Cfolmax_stock_unc(dn))**2
       end do
       likelihood = likelihood-tot_exp
    endif

    ! Cwood log-likelihood (i.e. branch, stem and CR)
    if (DATAin%nCwood_stock > 0) then
       ! Create vector of (Wood_t0 + Wood_t1) * 0.5
       mid_state = ( DATAin%M_POOLS(1:DATAin%nodays,4) + DATAin%M_POOLS(2:(DATAin%nodays+1),4) ) &
                 * 0.5d0
       ! Vectorised version of loop to estimate cost function
       tot_exp = sum(( (mid_state(DATAin%Cwood_stockpts(1:DATAin%nCwood_stock)) &
                       -DATAin%Cwood_stock(DATAin%Cwood_stockpts(1:DATAin%nCwood_stock)))&
                     / DATAin%Cwood_stock_unc(DATAin%Cwood_stockpts(1:DATAin%nCwood_stock)))**2)
       ! Combine with existing likelihood estimate
       likelihood = likelihood-tot_exp
    endif

    ! Croots log-likelihood
    if (DATAin%nCroots_stock > 0) then
       ! Create vector of (root_t0 + root_t1) * 0.5
       mid_state = ( DATAin%M_POOLS(1:DATAin%nodays,3) + DATAin%M_POOLS(2:(DATAin%nodays+1),3) ) &
                 * 0.5d0
       ! Vectorised version of loop to estimate cost function
       tot_exp = sum(( (mid_state(DATAin%Croots_stockpts(1:DATAin%nCroots_stock)) &
                       -DATAin%Croots_stock(DATAin%Croots_stockpts(1:DATAin%nCroots_stock)))&
                     / DATAin%Croots_stock_unc(DATAin%Croots_stockpts(1:DATAin%nCroots_stock)))**2)
       ! Combine with existing likelihood estimate
       likelihood = likelihood-tot_exp
    endif

! AGB calculation not currently possible in this version of DALEC_CROP_BUCKET
!    ! Cagb log-likelihood
!    if (DATAin%nCagb_stock > 0) then
!       ! Create vector of (agb_t0 + agb_t1) * 0.5
!       mid_state = ( DATAin%M_POOLS(1:DATAin%nodays,4) + DATAin%M_POOLS(2:(DATAin%nodays+1),4) ) &
!                 * 0.5d0 * (1d0-pars(29))
!       ! Vectorised version of loop to estimate cost function
!       tot_exp = sum(( (mid_state(DATAin%Cagb_stockpts(1:DATAin%nCagb_stock)) &
!                       -DATAin%Cagb_stock(DATAin%Cagb_stockpts(1:DATAin%nCagb_stock)))&
!                     / DATAin%Cagb_stock_unc(DATAin%Cagb_stockpts(1:DATAin%nCagb_stock)))**2)
!       ! Combine with existing likelihood estimate
!       likelihood = likelihood-tot_exp
!    endif

! Coarse root calculation not currently possible in this version of DALEC_CROP_BUCKET
!    ! Ccoarseroot log-likelihood
!    if (DATAin%nCcoarseroot_stock > 0) then
!       ! Create vector of (coarseroot_t0 + coarseroot_t1) * 0.5
!       mid_state = ( DATAin%M_POOLS(1:DATAin%nodays,4) + DATAin%M_POOLS(2:(DATAin%nodays+1),4) ) &
!                 * 0.5d0 * pars(29)
!       ! Vectorised version of loop to estimate cost function
!       tot_exp = sum(( (mid_state(DATAin%Ccoarseroot_stockpts(1:DATAin%nCcoarseroot_stock)) &
!                       -DATAin%Ccoarseroot_stock(DATAin%Ccoarseroot_stockpts(1:DATAin%nCcoarseroot_stock)))&
!                     / DATAin%Ccoarseroot_stock_unc(DATAin%Ccoarseroot_stockpts(1:DATAin%nCcoarseroot_stock)))**2)
!       ! Combine with existing likelihood estimate
!       likelihood = likelihood-tot_exp
!    endif

    ! Clitter log-likelihood
    ! WARNING WARNING WARNING hack in place to estimate fraction of litter pool
    ! originating from surface pools
    if (DATAin%nClit_stock > 0) then
       ! Create vector of (lit_t0 + lit_t1) * 0.5
       !mid_state = ( DATAin%M_POOLS(1:DATAin%nodays,4) + DATAin%M_POOLS(2:(DATAin%nodays+1),4) ) &
       !          * 0.5d0
       mid_state = (sum(DATAin%M_FLUXES(:,10))/sum(DATAin%M_FLUXES(:,10)+DATAin%M_FLUXES(:,12))) &
                 * DATAin%M_POOLS(:,5)
       mid_state = (mid_state(1:DATAin%nodays) + mid_state(2:(DATAin%nodays+1))) * 0.5d0
       ! Vectorised version of loop to estimate cost function
       tot_exp = sum(( (mid_state(DATAin%Clit_stockpts(1:DATAin%nClit_stock)) &
                       -DATAin%Clit_stock(DATAin%Clit_stockpts(1:DATAin%nClit_stock)))&
                     / DATAin%Clit_stock_unc(DATAin%Clit_stockpts(1:DATAin%nClit_stock)))**2)
       ! Combine with existing likelihood estimate
       likelihood = likelihood-tot_exp
    endif

    ! Csom log-likelihood
    if (DATAin%nCsom_stock > 0) then
       ! Create vector of (som_t0 + som_t1) * 0.5
       mid_state = ( DATAin%M_POOLS(1:DATAin%nodays,6) + DATAin%M_POOLS(2:(DATAin%nodays+1),6) ) &
                 * 0.5d0
       ! Vectorised version of loop to estimate cost function
       tot_exp = sum(( (mid_state(DATAin%Csom_stockpts(1:DATAin%nCsom_stock)) &
                       -DATAin%Csom_stock(DATAin%Csom_stockpts(1:DATAin%nCsom_stock)))&
                     / DATAin%Csom_stock_unc(DATAin%Csom_stockpts(1:DATAin%nCsom_stock)))**2)
       ! Combine with existing likelihood estimate
       likelihood = likelihood-tot_exp
    endif

    !
    ! Curiously we will assess 'other' priors here, as the tend to have to do with model state derived values
    !

    ! Ra:GPP fraction is in this model a derived property
    if (DATAin%otherpriors(1) > -9998) then
        tot_exp = sum(DATAin%M_FLUXES(:,3)) / sum(DATAin%M_FLUXES(:,1))
        tot_exp =  DATAin%otherpriorweight(1) * ((tot_exp-DATAin%otherpriors(1))/DATAin%otherpriorunc(1))**2
        likelihood = likelihood-tot_exp
    end if

    ! Initial soil water prior. The actual prior is linked to a fraction of field capacity so here is were that soil water at t=1
    ! is actually assessed against an observation
    if (DATAin%otherpriors(2) > -9998) then
        tot_exp = (DATAin%M_POOLS(1,8) * 1d-3) / layer_thickness(1) ! convert mm -> m3/m3
        tot_exp =  DATAin%otherpriorweight(2) * ((tot_exp-DATAin%otherpriors(2))/DATAin%otherpriorunc(2))**2
        likelihood = likelihood-tot_exp
    end if

    ! Leaf C:N is derived from multiple parameters
    if (DATAin%otherpriors(3) > -9998) then
        tot_exp = pars(17) / (10d0**pars(11))
        tot_exp =  DATAin%otherpriorweight(3) * ((tot_exp-DATAin%otherpriors(3))/DATAin%otherpriorunc(3))**2
        likelihood = likelihood-tot_exp
    end if

    ! Evaportranspiration (kgH2O/m2/day) as ratio of precipitation (kg/m2/s ->
    ! kg/m2/day)
    if (DATAin%otherpriors(4) > -9998) then
        tot_exp = sum(DATAin%M_FLUXES(:,19)) / sum(DATAin%MET(7,:) * 86400d0)
        tot_exp =  DATAin%otherpriorweight(4) * ((tot_exp-DATAin%otherpriors(4))/DATAin%otherpriorunc(4))**2
        likelihood = likelihood-tot_exp
    end if

    ! the likelihood scores for each observation are subject to multiplication
    ! by 0.5 in the algebraic formulation. To avoid repeated calculation across
    ! multiple datastreams we apply this multiplication to the bulk likelihood
    ! hear
    likelihood = likelihood * 0.5d0

    ! check that log-likelihood is an actual number
    if (likelihood /= likelihood) then
       likelihood = log(infini)
    end if
    ! don't forget to return
    return

  end function likelihood
  !
  !------------------------------------------------------------------
  !
  double precision function scale_likelihood(npars,pars)
    use cardamom_structures, only: DATAin
    use carbon_model_mod, only: layer_thickness

    ! calculates the likelihood of of the model output compared to the available
    ! observations which have been input to the model

    implicit none

    ! declare arguments
    integer, intent(in) :: npars
    double precision, dimension(npars), intent(in) :: pars

    ! declare local variables
    integer :: n, dn, no_years, y, s
    double precision :: tot_exp, tmp_var, infini, input, output
    double precision, allocatable :: mean_annual_pools(:)

    ! initial value
    scale_likelihood = 0d0 ; infini = 0d0

! Currently no distinction between NEE and NBE as fire is not simulated
!    ! NBE Log-likelihood
!    if (DATAin%nnbe > 0) then
!       tot_exp = sum((((DATAin%M_NEE(DATAin%nbepts(1:DATAin%nnbe))+DATAin%M_FLUXES(DATAin%nbepts(1:DATAin%nnbe),17)) &
!                       -DATAin%NBE(DATAin%nbepts(1:DATAin%nnbe))) &
!                       /DATAin%NBE_unc(DATAin%nbepts(1:DATAin%nnbe)))**2)
!       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nnbe))
!    endif

    ! Evap Log-likelihood
    if (DATAin%nEvap > 0) then
       tot_exp = sum(((DATAin%M_FLUXES(DATAin%Evappts(1:DATAin%nEvap),19)-DATAin%Evap(DATAin%Evappts(1:DATAin%nEvap))) &
                       /DATAin%Evap_unc(DATAin%evappts(1:DATAin%nEvap)))**2)
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nEvap))
    endif

! FIRE not currently coded into DALEC_CROP_BUCKET
!    ! Fire Log-likelihood
!    if (DATAin%nFire > 0) then
!       tot_exp = sum(((DATAin%M_FLUXES(DATAin%Firepts(1:DATAin%nFire),17)-DATAin%Fire(DATAin%Firepts(1:DATAin%nFire))) &
!                       /DATAin%Fire_unc(DATAin%Firepts(1:DATAin%nFire)))**2)
!       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nFire))
!    endif

    ! GPP Log-likelihood
    if (DATAin%ngpp > 0) then
       tot_exp = sum(((DATAin%M_GPP(DATAin%gpppts(1:DATAin%ngpp))-DATAin%GPP(DATAin%gpppts(1:DATAin%ngpp))) &
                       /DATAin%GPP_unc(DATAin%gpppts(1:DATAin%ngpp)))**2)
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%ngpp))
    endif

    ! LAI log-likelihood
    if (DATAin%nlai > 0) then
        ! loop split to allow vectorisation
        tot_exp = sum(((DATAin%M_LAI(DATAin%laipts(1:DATAin%nlai))-DATAin%LAI(DATAin%laipts(1:DATAin%nlai))) &
                       /DATAin%LAI_unc(DATAin%laipts(1:DATAin%nlai)))**2)
        if (minval(DATAin%M_LAI) < 0d0) tot_exp = tot_exp + (-log(infini))
        scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nlai))
    endif

    ! NEE likelihood
    if (DATAin%nnee > 0) then
       tot_exp = sum(((DATAin%M_NEE(DATAin%neepts(1:DATAin%nnee))-DATAin%NEE(DATAin%neepts(1:DATAin%nnee))) &
                       /DATAin%NEE_unc(DATAin%neepts(1:DATAin%nnee)))**2)
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nnee))
    endif

    ! Reco likelihood
    if (DATAin%nreco > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nreco
         dn = DATAin%recopts(n)
         tmp_var = DATAin%M_NEE(dn)+DATAin%M_GPP(dn)
         ! note that we calculate the Ecosystem resp from GPP and NEE
         tot_exp = tot_exp+((tmp_var-DATAin%Reco(dn))/DATAin%Reco_unc(dn))**2
       end do
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nreco))
    endif

    ! Cwood increment log-likelihood
    if (DATAin%nCwood_inc > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCwood_inc
         dn = DATAin%Cwood_incpts(n)
         s = max(0,dn-nint(DATAin%Cwood_inc_lag(dn)))+1
         ! Estimate the mean allocation to wood over the lag period
         tmp_var = sum(DATAin%M_FLUXES(s:dn,7)) / DATAin%Cwood_inc_lag(dn)
         tot_exp = tot_exp+((tmp_var-DATAin%Cwood_inc(dn)) / DATAin%Cwood_inc_unc(dn))**2
       end do
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nCwood_inc))
    endif

    ! Cwood mortality log-likelihood
    if (DATAin%nCwood_mortality > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCwood_mortality
         dn = DATAin%Cwood_mortalitypts(n)
         s = max(0,dn-nint(DATAin%Cwood_mortality_lag(dn)))+1
         ! Estimate the mean allocation to wood over the lag period
         tmp_var = sum(DATAin%M_FLUXES(s:dn,11)) / DATAin%Cwood_mortality_lag(dn)
         tot_exp = tot_exp+((tmp_var-DATAin%Cwood_mortality(dn)) / DATAin%Cwood_mortality_unc(dn))**2
       end do
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nCwood_mortality))
    endif

    ! Cfoliage log-likelihood
    if (DATAin%nCfol_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCfol_stock
         dn = DATAin%Cfol_stockpts(n)
         ! note that division is the uncertainty
         tot_exp = tot_exp+((DATAin%M_POOLS(dn,2)-DATAin%Cfol_stock(dn)) / DATAin%Cfol_stock_unc(dn))**2
       end do
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nCfol_stock))
    endif

    ! Annual foliar maximum
    if (DATAin%nCfolmax_stock > 0) then
       tot_exp = 0d0
       no_years = int(nint(sum(DATAin%deltat)/365.25d0))
       if (allocated(mean_annual_pools)) deallocate(mean_annual_pools)
       allocate(mean_annual_pools(no_years))
       ! determine the annual max for each pool
       do y = 1, no_years
          ! derive mean annual foliar pool
          mean_annual_pools(y) = cal_max_annual_pools(DATAin%M_POOLS(1:(DATAin%nodays+1),2),y,DATAin%deltat,DATAin%nodays+1)
       end do ! year loop
       ! loop through the observations then
       do n = 1, DATAin%nCfolmax_stock
         ! load the observation position in stream
         dn = DATAin%Cfolmax_stockpts(n)
         ! determine which years this in in for the simulation
         y = ceiling( (dble(dn)*(sum(DATAin%deltat)/(DATAin%nodays))) / 365.25d0 )
         ! load the correct year into the analysis
         tmp_var = mean_annual_pools(y)
         ! note that division is the uncertainty
         tot_exp = tot_exp+((tmp_var-DATAin%Cfolmax_stock(dn)) / DATAin%Cfolmax_stock_unc(dn))**2
       end do
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nCfolmax_stock))
    endif

    ! Cwood log-likelihood (i.e. branch, stem and CR)
    if (DATAin%nCwood_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCwood_stock
         dn = DATAin%Cwood_stockpts(n)
         ! note that division is the uncertainty
         tot_exp = tot_exp+((DATAin%M_POOLS(dn,4)-DATAin%Cwood_stock(dn))/DATAin%Cwood_stock_unc(dn))**2
       end do
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nCwood_stock))
    endif

! AGB calculation not currently implemented in this version of DALEC_CROP_BUCKET
!    ! Cagb log-likelihood
!    if (DATAin%nCagb_stock > 0) then
!       tot_exp = 0d0
!       do n = 1, DATAin%nCagb_stock
!         dn = DATAin%Cagb_stockpts(n)
!         ! remove coarse root fraction from wood (pars29)
!         tmp_var = DATAin%M_POOLS(dn,4)-(DATAin%M_POOLS(dn,4)*pars(29))
!         tot_exp = tot_exp+((tmp_var-DATAin%Cagb_stock(dn))/DATAin%Cagb_stock_unc(dn))**2
!       end do
!       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nCagb_stock))
!    endif

! Coarse root calculation not currently available in this version of DALEC_CROP_BUCKET
!    ! Ccoarseroot log-likelihood
!    if (DATAin%nCcoarseroot_stock > 0) then
!       tot_exp = 0d0
!       do n = 1, DATAin%nCcoarseroot_stock
!         dn = DATAin%Ccoarseroot_stockpts(n)
!         ! extract coarse root component from wood only
!         tmp_var = DATAin%M_POOLS(dn,4)*pars(29)
!         tot_exp = tot_exp+((tmp_var-DATAin%Ccoarseroot_stock(dn)) / DATAin%Ccoarseroot_stock_unc(dn))**2
!       end do
!       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nCcoarseroot_stock))
!    endif

    ! Croots log-likelihood
    if (DATAin%nCroots_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCroots_stock
         dn = DATAin%Croots_stockpts(n)
         ! note that division is the uncertainty
         tot_exp = tot_exp+((DATAin%M_POOLS(dn,3)-DATAin%Croots_stock(dn)) / DATAin%Croots_stock_unc(dn))**2
       end do
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nCroots_stock))
    endif

    ! Clitter log-likelihood
    ! WARNING WARNING WARNING hack in place to estimate fraction of litter pool
    ! originating from surface pools
    if (DATAin%nClit_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nClit_stock
         dn = DATAin%Clit_stockpts(n)
         ! note that division is the uncertainty
         tot_exp = tot_exp+(((sum(DATAin%M_FLUXES(:,10))/sum(DATAin%M_FLUXES(:,10)+DATAin%M_FLUXES(:,12))) &
                           *(DATAin%M_POOLS(dn,5))-DATAin%Clit_stock(dn))/DATAin%Clit_stock_unc(dn))**2
       end do
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nClit_stock))
    endif

    ! Csom log-likelihood
    if (DATAin%nCsom_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCsom_stock
         dn = DATAin%Csom_stockpts(n)
         ! note that division is the uncertainty
         tot_exp = tot_exp+((DATAin%M_POOLS(dn,6)-DATAin%Csom_stock(dn))/DATAin%Csom_stock_unc(dn))**2
       end do
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nCsom_stock))
    endif

    !
    ! Curiously we will assess 'other' priors here, as the tend to have to do with model state derived values
    !

    ! Ra:GPP fraction is in this model a derived property
    if (DATAin%otherpriors(1) > -9998) then
        tot_exp = sum(DATAin%M_FLUXES(:,3)) / sum(DATAin%M_FLUXES(:,1))
        tot_exp =  DATAin%otherpriorweight(1) * ((tot_exp-DATAin%otherpriors(1))/DATAin%otherpriorunc(1))**2
        scale_likelihood = scale_likelihood-tot_exp
    end if

    ! Initial soil water prior. The actual prior is linked to a fraction of field capacity so here is were that soil water at t=1
    ! is actually assessed against an observation
    if (DATAin%otherpriors(2) > -9998) then
        tot_exp = (DATAin%M_POOLS(1,8) * 1d-3) / layer_thickness(1) ! convert mm -> m3/m3
        tot_exp =  DATAin%otherpriorweight(2) * ((tot_exp-DATAin%otherpriors(2))/DATAin%otherpriorunc(2))**2
        scale_likelihood = scale_likelihood-tot_exp
    end if

    ! Leaf C:N is derived from multiple parameters
    if (DATAin%otherpriors(3) > -9998) then
        tot_exp = pars(17) / (10d0**pars(11))
        tot_exp =  DATAin%otherpriorweight(3) * ((tot_exp-DATAin%otherpriors(3))/DATAin%otherpriorunc(3))**2
        scale_likelihood = scale_likelihood-tot_exp
    end if

    ! Evaportranspiration (kgH2O/m2/day) as ratio of precipitation (kg/m2/s ->
    ! kg/m2/day)
    if (DATAin%otherpriors(4) > -9998) then
        tot_exp = sum(DATAin%M_FLUXES(:,19)) / sum(DATAin%MET(7,:) * 86400d0)
        tot_exp =  DATAin%otherpriorweight(4) * ((tot_exp-DATAin%otherpriors(4))/DATAin%otherpriorunc(4))**2
        scale_likelihood = scale_likelihood-tot_exp
    end if

    ! the likelihood scores for each observation are subject to multiplication
    ! by 0.5 in the algebraic formulation. To avoid repeated calculation across
    ! multiple datastreams we apply this multiplication to the bulk likelihood
    ! hear
    scale_likelihood = scale_likelihood * 0.5d0

    ! check that log-likelihood is an actual number
    if (scale_likelihood /= scale_likelihood) then
       scale_likelihood = log(infini)
    end if
    ! don't forget to return
    return

  end function scale_likelihood
  !
  !------------------------------------------------------------------
  !
  double precision function sqrt_scale_likelihood(npars,pars)
    use cardamom_structures, only: DATAin
    use carbon_model_mod, only: layer_thickness

    ! calculates the likelihood of of the model output compared to the available
    ! observations which have been input to the model

    implicit none

    ! declare arguments
    integer, intent(in) :: npars
    double precision, dimension(npars), intent(in) :: pars

    ! declare local variables
    integer :: n, dn, no_years, y, s
    double precision :: tot_exp, tmp_var, infini, input, output
    double precision, allocatable :: mean_annual_pools(:)

    ! initial value
    sqrt_scale_likelihood = 0d0 ; infini = 0d0

! Currently no distinction between NEE and NBE as fire is not simulated
!    ! NBE Log-likelihood
!    if (DATAin%nnbe > 0) then
!       tot_exp = sum((((DATAin%M_NEE(DATAin%nbepts(1:DATAin%nnbe))+DATAin%M_FLUXES(DATAin%nbepts(1:DATAin%nnbe),17)) &
!                       -DATAin%NBE(DATAin%nbepts(1:DATAin%nnbe))) &
!                       /DATAin%NBE_unc(DATAin%nbepts(1:DATAin%nnbe)))**2)
!       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nnbe)))
!    endif

    ! Evap Log-likelihood
    if (DATAin%nEvap > 0) then
       tot_exp = sum(((DATAin%M_FLUXES(DATAin%Evappts(1:DATAin%nEvap),19)-DATAin%Evap(DATAin%Evappts(1:DATAin%nEvap))) &
                       /DATAin%Evap_unc(DATAin%evappts(1:DATAin%nEvap)))**2)
       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nEvap)))
    endif

! FIRE not currently coded into DALEC_CROP_BUCKET
!    ! Fire Log-likelihood
!    if (DATAin%nFire > 0) then
!       tot_exp = sum(((DATAin%M_FLUXES(DATAin%Firepts(1:DATAin%nFire),17)-DATAin%Fire(DATAin%Firepts(1:DATAin%nFire))) &
!                       /DATAin%Fire_unc(DATAin%Firepts(1:DATAin%nFire)))**2)
!       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nFire)))
!    endif

    ! GPP Log-likelihood
    if (DATAin%ngpp > 0) then
       tot_exp = sum(((DATAin%M_GPP(DATAin%gpppts(1:DATAin%ngpp))-DATAin%GPP(DATAin%gpppts(1:DATAin%ngpp))) &
                       /DATAin%GPP_unc(DATAin%gpppts(1:DATAin%ngpp)))**2)
       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%ngpp)))
    endif

    ! LAI log-likelihood
    if (DATAin%nlai > 0) then
        ! loop split to allow vectorisation
        tot_exp = sum(((DATAin%M_LAI(DATAin%laipts(1:DATAin%nlai))-DATAin%LAI(DATAin%laipts(1:DATAin%nlai))) &
                       /DATAin%LAI_unc(DATAin%laipts(1:DATAin%nlai)))**2)
        if (minval(DATAin%M_LAI) < 0d0) tot_exp = tot_exp + (-log(infini))
        sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nlai)))
    endif

    ! NEE likelihood
    if (DATAin%nnee > 0) then
       tot_exp = sum(((DATAin%M_NEE(DATAin%neepts(1:DATAin%nnee))-DATAin%NEE(DATAin%neepts(1:DATAin%nnee))) &
                       /DATAin%NEE_unc(DATAin%neepts(1:DATAin%nnee)))**2)
       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nnee)))
    endif

    ! Reco likelihood
    if (DATAin%nreco > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nreco
         dn = DATAin%recopts(n)
         tmp_var = DATAin%M_NEE(dn)+DATAin%M_GPP(dn)
         ! note that we calculate the Ecosystem resp from GPP and NEE
         tot_exp = tot_exp+((tmp_var-DATAin%Reco(dn))/DATAin%Reco_unc(dn))**2
       end do
       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nreco)))
    endif

    ! Cwood increment log-likelihood
    if (DATAin%nCwood_inc > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCwood_inc
         dn = DATAin%Cwood_incpts(n)
         s = max(0,dn-nint(DATAin%Cwood_inc_lag(dn)))+1
         ! Estimate the mean allocation to wood over the lag period
         tmp_var = sum(DATAin%M_FLUXES(s:dn,7)) / DATAin%Cwood_inc_lag(dn)
         tot_exp = tot_exp+((tmp_var-DATAin%Cwood_inc(dn)) / DATAin%Cwood_inc_unc(dn))**2
       end do
       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nCwood_inc)))
    endif

    ! Cwood mortality log-likelihood
    if (DATAin%nCwood_mortality > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCwood_mortality
         dn = DATAin%Cwood_mortalitypts(n)
         s = max(0,dn-nint(DATAin%Cwood_mortality_lag(dn)))+1
         ! Estimate the mean allocation to wood over the lag period
         tmp_var = sum(DATAin%M_FLUXES(s:dn,11)) / DATAin%Cwood_mortality_lag(dn)
         tot_exp = tot_exp+((tmp_var-DATAin%Cwood_mortality(dn)) / DATAin%Cwood_mortality_unc(dn))**2
       end do
       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nCwood_mortality)))
    endif

    ! Cfoliage log-likelihood
    if (DATAin%nCfol_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCfol_stock
         dn = DATAin%Cfol_stockpts(n)
         ! note that division is the uncertainty
         tot_exp = tot_exp+((DATAin%M_POOLS(dn,2)-DATAin%Cfol_stock(dn)) / DATAin%Cfol_stock_unc(dn))**2
       end do
       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nCfol_stock)))
    endif

    ! Annual foliar maximum
    if (DATAin%nCfolmax_stock > 0) then
       tot_exp = 0d0
       no_years = int(nint(sum(DATAin%deltat)/365.25d0))
       if (allocated(mean_annual_pools)) deallocate(mean_annual_pools)
       allocate(mean_annual_pools(no_years))
       ! determine the annual max for each pool
       do y = 1, no_years
          ! derive mean annual foliar pool
          mean_annual_pools(y) = cal_max_annual_pools(DATAin%M_POOLS(1:(DATAin%nodays+1),2),y,DATAin%deltat,DATAin%nodays+1)
       end do ! year loop
       ! loop through the observations then
       do n = 1, DATAin%nCfolmax_stock
         ! load the observation position in stream
         dn = DATAin%Cfolmax_stockpts(n)
         ! determine which years this in in for the simulation
         y = ceiling( (dble(dn)*(sum(DATAin%deltat)/(DATAin%nodays))) / 365.25d0 )
         ! load the correct year into the analysis
         tmp_var = mean_annual_pools(y)
         ! note that division is the uncertainty
         tot_exp = tot_exp+((tmp_var-DATAin%Cfolmax_stock(dn)) / DATAin%Cfolmax_stock_unc(dn))**2
       end do
       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nCfolmax_stock)))
    endif

    ! Cwood log-likelihood (i.e. branch, stem and CR)
    if (DATAin%nCwood_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCwood_stock
         dn = DATAin%Cwood_stockpts(n)
         ! note that division is the uncertainty
         tot_exp = tot_exp+((DATAin%M_POOLS(dn,4)-DATAin%Cwood_stock(dn))/DATAin%Cwood_stock_unc(dn))**2
       end do
       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nCwood_stock)))
    endif

! AGB calculation not currently implemented in this version of DALEC_CROP_BUCKET
!    ! Cagb log-likelihood
!    if (DATAin%nCagb_stock > 0) then
!       tot_exp = 0d0
!       do n = 1, DATAin%nCagb_stock
!         dn = DATAin%Cagb_stockpts(n)
!         ! remove coarse root fraction from wood (pars29)
!         tmp_var = DATAin%M_POOLS(dn,4)-(DATAin%M_POOLS(dn,4)*pars(29))
!         tot_exp = tot_exp+((tmp_var-DATAin%Cagb_stock(dn))/DATAin%Cagb_stock_unc(dn))**2
!       end do
!       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nCagb_stock)))
!    endif

! Coarse root calculation not currently available in this version of DALEC_CROP_BUCKET
!    ! Ccoarseroot log-likelihood
!    if (DATAin%nCcoarseroot_stock > 0) then
!       tot_exp = 0d0
!       do n = 1, DATAin%nCcoarseroot_stock
!         dn = DATAin%Ccoarseroot_stockpts(n)
!         ! extract coarse root component from wood only
!         tmp_var = DATAin%M_POOLS(dn,4)*pars(29)
!         tot_exp = tot_exp+((tmp_var-DATAin%Ccoarseroot_stock(dn)) / DATAin%Ccoarseroot_stock_unc(dn))**2
!       end do
!       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nCcoarseroot_stock)))
!    endif

    ! Croots log-likelihood
    if (DATAin%nCroots_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCroots_stock
         dn = DATAin%Croots_stockpts(n)
         ! note that division is the uncertainty
         tot_exp = tot_exp+((DATAin%M_POOLS(dn,3)-DATAin%Croots_stock(dn)) / DATAin%Croots_stock_unc(dn))**2
       end do
       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nCroots_stock)))
    endif

    ! Clitter log-likelihood
    ! WARNING WARNING WARNING hack in place to estimate fraction of litter pool
    ! originating from surface pools
    if (DATAin%nClit_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nClit_stock
         dn = DATAin%Clit_stockpts(n)
         ! note that division is the uncertainty
         tot_exp = tot_exp+(((sum(DATAin%M_FLUXES(:,10))/sum(DATAin%M_FLUXES(:,10)+DATAin%M_FLUXES(:,12))) &
                           *(DATAin%M_POOLS(dn,5))-DATAin%Clit_stock(dn))/DATAin%Clit_stock_unc(dn))**2
       end do
       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nClit_stock)))
    endif

    ! Csom log-likelihood
    if (DATAin%nCsom_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCsom_stock
         dn = DATAin%Csom_stockpts(n)
         ! note that division is the uncertainty
         tot_exp = tot_exp+((DATAin%M_POOLS(dn,6)-DATAin%Csom_stock(dn))/DATAin%Csom_stock_unc(dn))**2
       end do
       sqrt_scale_likelihood = sqrt_scale_likelihood-(tot_exp/sqrt(dble(DATAin%nCsom_stock)))
    endif

    !
    ! Curiously we will assess 'other' priors here, as the tend to have to do with model state derived values
    !

    ! Ra:GPP fraction is in this model a derived property
    if (DATAin%otherpriors(1) > -9998) then
        tot_exp = sum(DATAin%M_FLUXES(:,3)) / sum(DATAin%M_FLUXES(:,1))
        tot_exp =  DATAin%otherpriorweight(1) * ((tot_exp-DATAin%otherpriors(1))/DATAin%otherpriorunc(1))**2
        sqrt_scale_likelihood = sqrt_scale_likelihood-tot_exp
    end if

    ! Initial soil water prior. The actual prior is linked to a fraction of field capacity so here is were that soil water at t=1
    ! is actually assessed against an observation
    if (DATAin%otherpriors(2) > -9998) then
        tot_exp = (DATAin%M_POOLS(1,8) * 1d-3) / layer_thickness(1) ! convert mm -> m3/m3
        tot_exp =  DATAin%otherpriorweight(2) * ((tot_exp-DATAin%otherpriors(2))/DATAin%otherpriorunc(2))**2
        sqrt_scale_likelihood = sqrt_scale_likelihood-tot_exp
    end if

    ! Leaf C:N is derived from multiple parameters
    if (DATAin%otherpriors(3) > -9998) then
        tot_exp = pars(17) / (10d0**pars(11))
        tot_exp =  DATAin%otherpriorweight(3) * ((tot_exp-DATAin%otherpriors(3))/DATAin%otherpriorunc(3))**2
        sqrt_scale_likelihood = sqrt_scale_likelihood-tot_exp
    end if

    ! Evaportranspiration (kgH2O/m2/day) as ratio of precipitation (kg/m2/s ->
    ! kg/m2/day)
    if (DATAin%otherpriors(4) > -9998) then
        tot_exp = sum(DATAin%M_FLUXES(:,19)) / sum(DATAin%MET(7,:) * 86400d0)
        tot_exp =  DATAin%otherpriorweight(4) * ((tot_exp-DATAin%otherpriors(4))/DATAin%otherpriorunc(4))**2
        sqrt_scale_likelihood = sqrt_scale_likelihood-tot_exp
    end if

    ! the likelihood scores for each observation are subject to multiplication
    ! by 0.5 in the algebraic formulation. To avoid repeated calculation across
    ! multiple datastreams we apply this multiplication to the bulk likelihood
    ! hear
    sqrt_scale_likelihood = sqrt_scale_likelihood * 0.5d0

    ! check that log-likelihood is an actual number
    if (sqrt_scale_likelihood /= sqrt_scale_likelihood) then
       sqrt_scale_likelihood = log(infini)
    end if
    ! don't forget to return
    return

  end function sqrt_scale_likelihood
  !
  !------------------------------------------------------------------
  !
  double precision function log_scale_likelihood(npars,pars)
    use cardamom_structures, only: DATAin
    use carbon_model_mod, only: layer_thickness

    ! calculates the likelihood of of the model output compared to the available
    ! observations which have been input to the model

    implicit none

    ! declare arguments
    integer, intent(in) :: npars
    double precision, dimension(npars), intent(in) :: pars

    ! declare local variables
    integer :: n, dn, no_years, y, s
    double precision :: tot_exp, tmp_var, infini, input, output
    double precision, allocatable :: mean_annual_pools(:)

    ! initial value
    log_scale_likelihood = 0d0 ; infini = 0d0

! Currently no distinction between NEE and NBE as fire is not simulated
!    ! NBE Log-likelihood
!    if (DATAin%nnbe > 0) then
!       tot_exp = sum((((DATAin%M_NEE(DATAin%nbepts(1:DATAin%nnbe))+DATAin%M_FLUXES(DATAin%nbepts(1:DATAin%nnbe),17)) &
!                       -DATAin%NBE(DATAin%nbepts(1:DATAin%nnbe))) &
!                       /DATAin%NBE_unc(DATAin%nbepts(1:DATAin%nnbe)))**2)
!       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0-log(dble(DATAin%nnbe))))
!    endif

    ! Evap Log-likelihood
    if (DATAin%nEvap > 0) then
       tot_exp = sum(((DATAin%M_FLUXES(DATAin%Evappts(1:DATAin%nEvap),19)-DATAin%Evap(DATAin%Evappts(1:DATAin%nEvap))) &
                       /DATAin%Evap_unc(DATAin%evappts(1:DATAin%nEvap)))**2)
       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0-log(dble(DATAin%nEvap))))
    endif

! FIRE not currently coded into DALEC_CROP_BUCKET
!    ! Fire Log-likelihood
!    if (DATAin%nFire > 0) then
!       tot_exp = sum(((DATAin%M_FLUXES(DATAin%Firepts(1:DATAin%nFire),17)-DATAin%Fire(DATAin%Firepts(1:DATAin%nFire))) &
!                       /DATAin%Fire_unc(DATAin%Firepts(1:DATAin%nFire)))**2)
!       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0-log(dble(DATAin%nFire))))
!    endif

    ! GPP Log-likelihood
    if (DATAin%ngpp > 0) then
       tot_exp = sum(((DATAin%M_GPP(DATAin%gpppts(1:DATAin%ngpp))-DATAin%GPP(DATAin%gpppts(1:DATAin%ngpp))) &
                       /DATAin%GPP_unc(DATAin%gpppts(1:DATAin%ngpp)))**2)
       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0-log(dble(DATAin%ngpp))))
    endif

    ! LAI log-likelihood
    if (DATAin%nlai > 0) then
        ! loop split to allow vectorisation
        tot_exp = sum(((DATAin%M_LAI(DATAin%laipts(1:DATAin%nlai))-DATAin%LAI(DATAin%laipts(1:DATAin%nlai))) &
                       /DATAin%LAI_unc(DATAin%laipts(1:DATAin%nlai)))**2)
        if (minval(DATAin%M_LAI) < 0d0) tot_exp = tot_exp + (-log(infini))
        log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0-log(dble(DATAin%nlai))))
    endif

    ! NEE likelihood
    if (DATAin%nnee > 0) then
       tot_exp = sum(((DATAin%M_NEE(DATAin%neepts(1:DATAin%nnee))-DATAin%NEE(DATAin%neepts(1:DATAin%nnee))) &
                       /DATAin%NEE_unc(DATAin%neepts(1:DATAin%nnee)))**2)
       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0-log(dble(DATAin%nnee))))
    endif

    ! Reco likelihood
    if (DATAin%nreco > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nreco
         dn = DATAin%recopts(n)
         tmp_var = DATAin%M_NEE(dn)+DATAin%M_GPP(dn)
         ! note that we calculate the Ecosystem resp from GPP and NEE
         tot_exp = tot_exp+((tmp_var-DATAin%Reco(dn))/DATAin%Reco_unc(dn))**2
       end do
       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0-log(dble(DATAin%nreco))))
    endif

    ! Cwood increment log-likelihood
    if (DATAin%nCwood_inc > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCwood_inc
         dn = DATAin%Cwood_incpts(n)
         s = max(0,dn-nint(DATAin%Cwood_inc_lag(dn)))+1
         ! Estimate the mean allocation to wood over the lag period
         tmp_var = sum(DATAin%M_FLUXES(s:dn,7)) / DATAin%Cwood_inc_lag(dn)
         tot_exp = tot_exp+((tmp_var-DATAin%Cwood_inc(dn)) / DATAin%Cwood_inc_unc(dn))**2
       end do
       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0-log(dble(DATAin%nCwood_inc))))
    endif

    ! Cwood mortality log-likelihood
    if (DATAin%nCwood_mortality > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCwood_mortality
         dn = DATAin%Cwood_mortalitypts(n)
         s = max(0,dn-nint(DATAin%Cwood_mortality_lag(dn)))+1
         ! Estimate the mean allocation to wood over the lag period
         tmp_var = sum(DATAin%M_FLUXES(s:dn,11)) / DATAin%Cwood_mortality_lag(dn)
         tot_exp = tot_exp+((tmp_var-DATAin%Cwood_mortality(dn)) / DATAin%Cwood_mortality_unc(dn))**2
       end do
       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0-log(dble(DATAin%nCwood_mortality))))
    endif

    ! Cfoliage log-likelihood
    if (DATAin%nCfol_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCfol_stock
         dn = DATAin%Cfol_stockpts(n)
         ! note that division is the uncertainty
         tot_exp = tot_exp+((DATAin%M_POOLS(dn,2)-DATAin%Cfol_stock(dn)) / DATAin%Cfol_stock_unc(dn))**2
       end do
       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0-log(dble(DATAin%nCfol_stock))))
    endif

    ! Annual foliar maximum
    if (DATAin%nCfolmax_stock > 0) then
       tot_exp = 0d0
       no_years = int(nint(sum(DATAin%deltat)/365.25d0))
       if (allocated(mean_annual_pools)) deallocate(mean_annual_pools)
       allocate(mean_annual_pools(no_years))
       ! determine the annual max for each pool
       do y = 1, no_years
          ! derive mean annual foliar pool
          mean_annual_pools(y) = cal_max_annual_pools(DATAin%M_POOLS(1:(DATAin%nodays+1),2),y,DATAin%deltat,DATAin%nodays+1)
       end do ! year loop
       ! loop through the observations then
       do n = 1, DATAin%nCfolmax_stock
         ! load the observation position in stream
         dn = DATAin%Cfolmax_stockpts(n)
         ! determine which years this in in for the simulation
         y = ceiling( (dble(dn)*(sum(DATAin%deltat)/(DATAin%nodays))) / 365.25d0 )
         ! load the correct year into the analysis
         tmp_var = mean_annual_pools(y)
         ! note that division is the uncertainty
         tot_exp = tot_exp+((tmp_var-DATAin%Cfolmax_stock(dn)) / DATAin%Cfolmax_stock_unc(dn))**2
       end do
       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0-log(dble(DATAin%nCfolmax_stock))))
    endif

    ! Cwood log-likelihood (i.e. branch, stem and CR)
    if (DATAin%nCwood_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCwood_stock
         dn = DATAin%Cwood_stockpts(n)
         ! note that division is the uncertainty
         tot_exp = tot_exp+((DATAin%M_POOLS(dn,4)-DATAin%Cwood_stock(dn))/DATAin%Cwood_stock_unc(dn))**2
       end do
       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0-log(dble(DATAin%nCwood_stock))))
    endif

! AGB calculation not currently implemented in this version of DALEC_CROP_BUCKET
!    ! Cagb log-likelihood
!    if (DATAin%nCagb_stock > 0) then
!       tot_exp = 0d0
!       do n = 1, DATAin%nCagb_stock
!         dn = DATAin%Cagb_stockpts(n)
!         ! remove coarse root fraction from wood (pars29)
!         tmp_var = DATAin%M_POOLS(dn,4)-(DATAin%M_POOLS(dn,4)*pars(29))
!         tot_exp = tot_exp+((tmp_var-DATAin%Cagb_stock(dn))/DATAin%Cagb_stock_unc(dn))**2
!       end do
!       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0-log(dble(DATAin%nCagb_stock))))
!    endif

! Coarse root calculation not currently available in this version of DALEC_CROP_BUCKET
!    ! Ccoarseroot log-likelihood
!    if (DATAin%nCcoarseroot_stock > 0) then
!       tot_exp = 0d0
!       do n = 1, DATAin%nCcoarseroot_stock
!         dn = DATAin%Ccoarseroot_stockpts(n)
!         ! extract coarse root component from wood only
!         tmp_var = DATAin%M_POOLS(dn,4)*pars(29)
!         tot_exp = tot_exp+((tmp_var-DATAin%Ccoarseroot_stock(dn)) / DATAin%Ccoarseroot_stock_unc(dn))**2
!       end do
!       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0-log(dble(DATAin%nCcoarseroot_stock))))
!    endif

    ! Croots log-likelihood
    if (DATAin%nCroots_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCroots_stock
         dn = DATAin%Croots_stockpts(n)
         ! note that division is the uncertainty
         tot_exp = tot_exp+((DATAin%M_POOLS(dn,3)-DATAin%Croots_stock(dn)) / DATAin%Croots_stock_unc(dn))**2
       end do
       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0-log(dble(DATAin%nCroots_stock))))
    endif

    ! Clitter log-likelihood
    ! WARNING WARNING WARNING hack in place to estimate fraction of litter pool
    ! originating from surface pools
    if (DATAin%nClit_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nClit_stock
         dn = DATAin%Clit_stockpts(n)
         ! note that division is the uncertainty
         tot_exp = tot_exp+(((sum(DATAin%M_FLUXES(:,10))/sum(DATAin%M_FLUXES(:,10)+DATAin%M_FLUXES(:,12))) &
                           *(DATAin%M_POOLS(dn,5))-DATAin%Clit_stock(dn))/DATAin%Clit_stock_unc(dn))**2
       end do
       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0-log(dble(DATAin%nClit_stock))))
    endif

    ! Csom log-likelihood
    if (DATAin%nCsom_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCsom_stock
         dn = DATAin%Csom_stockpts(n)
         ! note that division is the uncertainty
         tot_exp = tot_exp+((DATAin%M_POOLS(dn,6)-DATAin%Csom_stock(dn))/DATAin%Csom_stock_unc(dn))**2
       end do
       log_scale_likelihood = log_scale_likelihood-(tot_exp/(1d0-log(dble(DATAin%nCsom_stock))))
    endif

    !
    ! Curiously we will assess 'other' priors here, as the tend to have to do with model state derived values
    !

    ! Ra:GPP fraction is in this model a derived property
    if (DATAin%otherpriors(1) > -9998) then
        tot_exp = sum(DATAin%M_FLUXES(:,3)) / sum(DATAin%M_FLUXES(:,1))
        tot_exp =  DATAin%otherpriorweight(1) * ((tot_exp-DATAin%otherpriors(1))/DATAin%otherpriorunc(1))**2
        log_scale_likelihood = log_scale_likelihood-tot_exp
    end if

    ! Initial soil water prior. The actual prior is linked to a fraction of field capacity so here is were that soil water at t=1
    ! is actually assessed against an observation
    if (DATAin%otherpriors(2) > -9998) then
        tot_exp = (DATAin%M_POOLS(1,8) * 1d-3) / layer_thickness(1) ! convert mm -> m3/m3
        tot_exp =  DATAin%otherpriorweight(2) * ((tot_exp-DATAin%otherpriors(2))/DATAin%otherpriorunc(2))**2
        log_scale_likelihood = log_scale_likelihood-tot_exp
    end if

    ! Leaf C:N is derived from multiple parameters
    if (DATAin%otherpriors(3) > -9998) then
        tot_exp = pars(17) / (10d0**pars(11))
        tot_exp =  DATAin%otherpriorweight(3) * ((tot_exp-DATAin%otherpriors(3))/DATAin%otherpriorunc(3))**2
        log_scale_likelihood = log_scale_likelihood-tot_exp
    end if

    ! Evaportranspiration (kgH2O/m2/day) as ratio of precipitation (kg/m2/s ->
    ! kg/m2/day)
    if (DATAin%otherpriors(4) > -9998) then
        tot_exp = sum(DATAin%M_FLUXES(:,19)) / sum(DATAin%MET(7,:) * 86400d0)
        tot_exp =  DATAin%otherpriorweight(4) * ((tot_exp-DATAin%otherpriors(4))/DATAin%otherpriorunc(4))**2
        log_scale_likelihood = log_scale_likelihood-tot_exp
    end if

    ! the likelihood scores for each observation are subject to multiplication
    ! by 0.5 in the algebraic formulation. To avoid repeated calculation across
    ! multiple datastreams we apply this multiplication to the bulk likelihood
    ! hear
    log_scale_likelihood = log_scale_likelihood * 0.5d0

    ! check that log-likelihood is an actual number
    if (log_scale_likelihood /= log_scale_likelihood) then
       log_scale_likelihood = log(infini)
    end if
    ! don't forget to return
    return

  end function log_scale_likelihood
  !
  !------------------------------------------------------------------
  !
end module model_likelihood_module
