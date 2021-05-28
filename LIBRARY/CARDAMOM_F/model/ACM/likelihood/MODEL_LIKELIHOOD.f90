
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
  public :: model_likelihood, find_edc_initial_values, sub_model_likelihood

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
    MCO%nOUT = 1000
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
    ! appropriate initial parameter choices, consistent with EDCs for DALEC2 /
    ! DALEC_GSI

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
    call EDC1_GSI(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)

    ! next need to run the model itself
    call carbon_model(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT,DATAin%M_LAI,DATAin%M_NEE &
                     ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%nopars &
                     ,DATAin%nomet,DATAin%nopools,DATAin%nofluxes  &
                     ,DATAin%M_GPP)

    ! assess post running EDCs
    call EDC2_GSI(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools &
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
!   ! for testing purposes, stop the model when start achieved
!    if (sum(EDCD%PASSFAIL) == 100) then
!        print*,"Found it!" ; stop
!    endif

    ! convert to a probability
    ML_obs_out = -0.5d0*(tot_exp*10d0)*DATAin%EDC

  end subroutine edc_model_likelihood
  !
  !------------------------------------------------------------------
  !
  subroutine sub_model_likelihood(PARS,ML_obs_out,ML_prior_out)
    use MCMCOPT, only:  PI
    use CARBON_MODEL_MOD, only: carbon_model
    use cardamom_structures, only: DATAin

    ! this subroutine is responsible for running the model,
    ! calculation of the log-likelihood on a subsample of observation
    ! for comparison assessment of parameter performance and use of the EDCs if they are
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
    ! if == 1 then all EDCs are checked irrespective of whether or not one has
    ! failed
    EDCD%DIAG = 0

    if (DATAin%EDC == 1) then

        ! call EDCs which can be evaluated prior to running the model
        call EDC1_GSI(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)

        ! update the likelihood score based on EDCs driving total rejection
        ! proposed parameters
        ML_obs_out = log(EDC1)

    endif !

    ! run the dalec model
    call carbon_model(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT,DATAin%M_LAI,DATAin%M_NEE &
                     ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%nopars &
                     ,DATAin%nomet,DATAin%nopools,DATAin%nofluxes  &
                     ,DATAin%M_GPP)

    ! if first set of EDCs have been passed, move on to the second
    if (DATAin%EDC == 1) then

        ! check edc2
        call EDC2_GSI(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools &
                     ,DATAin%nodays,DATAin%deltat,PI%parmax,PARS,DATAin%MET &
                     ,DATAin%M_LAI,DATAin%M_NEE,DATAin%M_GPP,DATAin%M_POOLS &
                     ,DATAin%M_FLUXES,DATAin%meantemp,EDC2)

        ! Add EDC2 log-likelihood to absolute accept reject...
        ML_obs_out = ML_obs_out + log(EDC2)

    end if ! DATAin%EDC == 1

    ! Calculate log-likelihood associated with priors
    ! We always want this
    ML_prior_out = likelihood_p(PI%npars,DATAin%parpriors,DATAin%parpriorunc,PARS)

    ! calculate final model likelihood when compared to obs
!    ML_obs_out = ML_obs_out + sub_likelihood(PI%npars,PARS)
    ML_obs_out = ML_obs_out + scale_likelihood(PI%npars,PARS)

  end subroutine sub_model_likelihood
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
  subroutine EDC1_GSI(PARS, npars, meantemp, meanrad, EDC1)

    use cardamom_structures, only: DATAin
    use CARBON_MODEL_MOD, only: sw_par_fraction, &
                                opt_max_scaling, &
                                     emissivity

    ! subroutine assessed the current parameter sets for passing ecological and
    ! steady state contraints (modified from Bloom et al., 2014).

    implicit none

    ! declare input variables
    integer, intent(in) :: npars ! number of parameters
    double precision, intent(out) :: EDC1    ! EDC1 flag
    double precision, dimension(npars), intent(in) :: PARS ! current parameter set
    double precision, intent(in) :: meantemp & ! mean temperature (k)
                                   ,meanrad    ! mean radiation (MJ.m-2.day-1)

    ! declare local parameters
    double precision, dimension(10), parameter :: lai = (/0.000001d0,0.03125d0,0.0625d0,0.125d0,0.25d0,0.5d0,1d0,2.5d0,5d0,10d0/)

    ! declare local variables
    integer :: n, DIAG, i
    double precision :: tmp, tmp1, tmp2 &
             ,direct_trans, par_trans, nir_trans &
             ,par_refl, nir_refl, lw_trans, lw_refl &
             ,fauto & ! Fractions of GPP to autotrophic respiration
             ,ffol  & ! Fraction of GPP to foliage
             ,flab  & ! Fraction of GPP to labile pool
             ,froot & ! Fraction of GPP to root
             ,fwood & ! Fraction of GPP to wood
             ,fsom    ! fraction of GPP som under eqilibrium conditions

    ! set initial value
    EDC1 = 1
    DIAG = EDCD%DIAG


    ! set all EDCs to 1 (pass)
    EDCD%nedc = 100
    EDCD%PASSFAIL(1:EDCD%nedc) = 1

    !
    ! begin checking EDCs
    !

    ! At high LAI the amount of radiation transmitted to the soil surface should
    ! be very low. On order of < 20 % of original SW input should reach the soil
    ! by lai > 5. NOTE: 1) 20 % target comes partially from SPA, 2) below code
    ! assumes implicit 1 MJ/m2/day

!    ! Estimate the fraction which passes directly to surface.
!    ! NOTE: 1 = clumping factor, -0.5 is decay coefficient
!    direct_trans = exp(-0.5d0 * lai(9) * 1d0)
!    ! Estimate the transmittance fraction of the canopy incident SW
!    par_trans = pars(16) * lai(9) + pars(22)
!    nir_trans = pars(17) * lai(9) + pars(23)
!    ! Estimate the fraction of input SW which reaches the surface
!    tmp = (1d0 * (1d0-direct_trans)) ! estimate the SW total incident on the canopy
!    direct_trans = (1d0 * direct_trans) + ((tmp*sw_par_fraction)*par_trans) + ((tmp*(1d0-sw_par_fraction))*nir_trans)
!    if ((EDC1 == 1 .or. DIAG == 1) .and. direct_trans > 0.20d0) then
!         EDC1 = 0 ; EDCD%PASSFAIL(6) = 0
!    endif

    ! the transmittance of NIR should always be > than PAR at all LAI values
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(15) < pars(14)) then
         EDC1 = 0 ; EDCD%PASSFAIL(1) = 0
    endif
    ! Reflectance of incident radiation should be greater than transmittance
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(15) > pars(17) ) then
         EDC1 = 0 ; EDCD%PASSFAIL(2) = 0
    endif
    ! Reflectance of incident radiation should be greater than transmittance
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(14) > pars(16)) then
         EDC1 = 0 ; EDCD%PASSFAIL(3) = 0
    endif
    ! reflectance should be greater for NIR than for PAR
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(17) < pars(16)) then
         EDC1 = 0 ; EDCD%PASSFAIL(4) = 0
    endif
    ! assume that photosynthesis limitation at 0C should be between 10 % and 20 %
    ! of potential. Fatchi et al (2013), New Phytologist, https://doi.org/10.1111/nph.12614
    tmp = opt_max_scaling(pars(2),pars(3),pars(4),0d0)
    if ((EDC1 == 1 .or. DIAG == 1) .and. tmp > 0.20d0) then
       EDC1 = 0 ; EDCD%PASSFAIL(5) = 0
    endif

    ! Longwave radiation incident on the canopy should be almost entirely absorbed by the canopy
    ! Here we assume that such a circumstance should occur by LAI >= 1
!    lw_trans = 0.02d0 * (1d0 - (lai(7) / (lai(7)+pars(6)))) ; lw_refl = 0.02d0 * (1d0 - (lai(7) / (lai(7)+pars(7))))
!    if ((EDC1 == 1 .or. DIAG == 1) .and. (1d0 - lw_trans - lw_refl) < 0.99d0) then
!        EDC1 = 0 ; EDCD%PASSFAIL(10) = 0
!    endif

  end subroutine EDC1_GSI
  !
  !------------------------------------------------------------------
  !
  subroutine EDC2_GSI(npars,nomet,nofluxes,nopools,nodays,deltat &
                      ,parmax,pars,met,M_LAI,M_NEE,M_GPP,M_POOLS,M_FLUXES &
                      ,meantemp,EDC2)

    use cardamom_structures, only: DATAin

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
    integer :: n, DIAG, no_years, y, PEDC, nn, num_EDC, max_location(1),i
    double precision :: mean_pools(nopools), G, decay_coef, meangpp, EQF &
                       ,model_living_C, target_living_C(2),hold
    double precision, dimension(:), allocatable :: mean_annual_pools,tmp
    double precision :: max_wood & !
                       ,torfol   & ! yearly average turnover
                       ,torlab   & !
                       ,fauto    & ! Fractions of GPP to autotrophic respiration
                       ,ffol     & ! Fraction of GPP to foliage
                       ,flab     & ! Fraction of GPP to labile pool
                       ,froot    & ! Fraction of GPP to root
                       ,fwood    & ! Fraction of GPP to wood
                       ,fsom     & ! fraction of GPP som under eqilibrium conditions
                       ,flit     & ! fraction of GPP to litter under equilibrium condition
                       ,delta_gsi

    ! set initial values
    DIAG = EDCD%DIAG
    EDC2 = 1

    !
    ! EDCs done, below are additional fault detection conditions
    !

!    ! the maximum value for all fluxes must be greater than zero
!    if ((EDC2 == 1 .or. DIAG == 1) .and. maxval(M_FLUXES(1:nodays,1)) < 0.1d0) then
!         EDC2 = 0 ; EDCD%PASSFAIL(9) = 0
!    endif
!    ! the maximum value for all fluxes must be greater than zero
!    if ((EDC2 == 1 .or. DIAG == 1) .and. maxval(M_FLUXES(1:nodays,2)) < 0.1d0) then
!         EDC2 = 0 ; EDCD%PASSFAIL(10) = 0
!    endif
!    ! the maximum value for all fluxes must be greater than zero
!    if ((EDC2 == 1 .or. DIAG == 1) .and. maxval(M_FLUXES(1:nodays,3)) < 0.1d0) then
!         EDC2 = 0 ; EDCD%PASSFAIL(11) = 0
!    endif
!    ! the maximum value for all fluxes must be greater than zero
!    if ((EDC2 == 1 .or. DIAG == 1) .and. maxval(M_FLUXES(1:nodays,4)) < 0.1d0) then
!         EDC2 = 0 ; EDCD%PASSFAIL(12) = 0
!    endif

    ! additional faults can be stored in locations 35 - 40 of the PASSFAIL array

    ! ensure minimum pool values are >= 0 and /= NaN
    if (EDC2 == 1 .or. DIAG == 1) then
       n=1
       do while (n <= nopools .and. (EDC2 == 1 .or. DIAG == 1))
          nn = 1 ; PEDC = 1
          do while (nn <= (nodays+1) .and. PEDC == 1)
             ! now check conditions
             if (M_POOLS(nn,n) < 0d0 .or. M_POOLS(nn,n) /= M_POOLS(nn,n)) then
                 EDC2 = 0 ; PEDC = 0 ; EDCD%PASSFAIL(35+n) = 0
             end if ! less than zero and is NaN condition
          nn = nn + 1
          end do ! nn < nodays .and. PEDC == 1
          n = n + 1
       end do ! for nopools .and. EDC .or. DIAG condition
    end if ! min pool assessment

    ! ensure FLUXES values are /= NaN
    if (EDC2 == 1 .or. DIAG == 1) then
       n=1
       do while (n <= nofluxes .and. (EDC2 == 1 .or. DIAG == 1))
          nn = 1 ; PEDC = 1
          do while (nn <= (nodays+1) .and. PEDC == 1)
             ! now check conditions
             if (M_FLUXES(nn,n) /= M_FLUXES(nn,n)) then
                 EDC2 = 0 ; PEDC = 0 ; EDCD%PASSFAIL(35+nopools+n) = 0
             end if ! less than zero and is NaN condition
          nn = nn + 1
          end do ! nn < nodays .and. PEDC == 1
          n = n + 1
       end do ! for nopools .and. EDC .or. DIAG condition
    end if ! min pool assessment

  end subroutine EDC2_GSI
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
        call EDC1_GSI(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)

        ! update the likelihood score based on EDCs driving total rejection
        ! proposed parameters
        ML_obs_out = log(EDC1)

    endif !

    ! run the dalec model
    call carbon_model(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT,DATAin%M_LAI,DATAin%M_NEE &
                     ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%nopars &
                     ,DATAin%nomet,DATAin%nopools,DATAin%nofluxes  &
                     ,DATAin%M_GPP)

    ! if first set of EDCs have been passed, move on to the second
    if (DATAin%EDC == 1) then

        ! check edc2
        call EDC2_GSI(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools &
                     ,DATAin%nodays,DATAin%deltat,PI%parmax,PARS,DATAin%MET &
                     ,DATAin%M_LAI,DATAin%M_NEE,DATAin%M_GPP,DATAin%M_POOLS &
                     ,DATAin%M_FLUXES,DATAin%meantemp,EDC2)

        ! Add EDC2 log-likelihood to absolute accept reject...
        ML_obs_out = ML_obs_out + log(EDC2)

    end if ! DATAin%EDC == 1

    ! Calculate log-likelihood associated with priors
    ! We always want this
    ML_prior_out = likelihood_p(PI%npars,DATAin%parpriors,DATAin%parpriorunc,PARS)
    ! calculate final model likelihood when compared to obs
    ML_obs_out = ML_obs_out + likelihood(PI%npars,PARS)

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
    likelihood_p = 0d0

    ! now loop through defined parameters for their uncertainties
    do n = 1, npars
       ! if there is actually a value
       if (parpriors(n) > -9999d0) then
           if (n >= 1 .or. n <= 19) then
               ! Gaussian uncertainties for CO2 compensation and half saturation
               ! points
               likelihood_p = likelihood_p-0.5d0*((pars(n)-parpriors(n))/parpriorunc(n))**2
           else
               likelihood_p = likelihood_p-0.5d0*(log(pars(n)/parpriors(n))/log(parpriorunc(n)))**2
           endif
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
    likelihood = 0d0 ; infini = 0d0

    ! GPP Log-likelihood
    tot_exp = 0d0
    if (DATAin%ngpp > 0) then
        do n = 1, DATAin%ngpp
          dn = DATAin%gpppts(n)
          ! note that division is the uncertainty
          tot_exp = tot_exp+((DATAin%M_GPP(dn)-DATAin%GPP(dn))/DATAin%GPP_unc(dn))**2
        end do
        likelihood = likelihood-0.5d0*tot_exp
    endif

    ! Evapotranspiration (kgH2O.m-2.day-1) Log-likelihood
    ! in this case transpiration only
    tot_exp = 0d0
    if (DATAin%nEvap > 0) then
        do n = 1, DATAin%nEvap
          dn = DATAin%Evappts(n)
          ! note that division is the uncertainty
          tot_exp = tot_exp+((DATAin%M_FLUXES(dn,2)-DATAin%Evap(dn))/DATAin%Evap_unc(dn))**2
        end do
        likelihood = likelihood-0.5d0*tot_exp
    endif

    ! Borrowed wood increment to provide soil evaporation for ACM recal  (kgH2O.m-2.day-1) Log-likelihood
    tot_exp = 0d0
    if (DATAin%nwoo > 0) then
        do n = 1, DATAin%nwoo
          dn = DATAin%woopts(n)
          ! note that division is the uncertainty
          tot_exp = tot_exp+((DATAin%M_FLUXES(dn,3)-DATAin%woo(dn))/DATAin%woo_unc(dn))**2
        end do
        likelihood = likelihood-0.5d0*tot_exp
    endif

    ! Borrowed Cfol_stock to provide wet canopy evaporation for ACM recal  (kgH2O.m-2.day-1) Log-likelihood
!    tot_exp = 0d0
!    if (DATAin%nCfol_stock > 0) then
!        do n = 1, DATAin%nCfol_stock
!          dn = DATAin%Cfol_stockpts(n)
!          ! note that division is the uncertainty
!          tot_exp = tot_exp+((DATAin%M_FLUXES(dn,4)-DATAin%Cfol_stock(dn))/DATAin%Cfol_stock_unc(dn))**2
!        end do
!        likelihood = likelihood-0.5d0*tot_exp
!    endif

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
    scale_likelihood = 0d0 ; infini = 0d0

    ! GPP Log-likelihood
    tot_exp = 0d0
    if (DATAin%ngpp > 0) then
        do n = 1, DATAin%ngpp
          dn = DATAin%gpppts(n)
          ! note that division is the uncertainty
          tot_exp = tot_exp+((DATAin%M_GPP(dn)-DATAin%GPP(dn))/DATAin%GPP_unc(dn))**2
        end do
        scale_likelihood = scale_likelihood-0.5d0*(tot_exp/dble(DATAin%ngpp))
    endif

    ! Evapotranspiration (kgH2O.m-2.day-1) Log-likelihood
    ! in this case transpiration only
    tot_exp = 0d0
    if (DATAin%nEvap > 0) then
        do n = 1, DATAin%nEvap
          dn = DATAin%Evappts(n)
          ! note that division is the uncertainty
          tot_exp = tot_exp+((DATAin%M_FLUXES(dn,2)-DATAin%Evap(dn))/DATAin%Evap_unc(dn))**2
        end do
        scale_likelihood = scale_likelihood-0.5d0*(tot_exp/dble(DATAin%nEvap))
    endif

    ! Borrowed wood increment to provide soil evaporation for ACM recal (kgH2O.m-2.day-1) Log-likelihood
    tot_exp = 0d0
    if (DATAin%nwoo > 0) then
        do n = 1, DATAin%nwoo
          dn = DATAin%woopts(n)
          ! note that division is the uncertainty
          tot_exp = tot_exp+((DATAin%M_FLUXES(dn,3)-DATAin%woo(dn))/DATAin%woo_unc(dn))**2
       end do
        scale_likelihood = scale_likelihood-0.5d0*(tot_exp/dble(DATAin%nwoo))
    endif

    ! Borrowed Cfol_stock to provide wet canopy evaporation for ACM recal  (kgH2O.m-2.day-1) Log-likelihood
!    tot_exp = 0d0
!    if (DATAin%nCfol_stock > 0) then
!        do n = 1, DATAin%nCfol_stock
!          dn = DATAin%Cfol_stockpts(n)
!          ! note that division is the uncertainty
!          tot_exp = tot_exp+((DATAin%M_FLUXES(dn,4)-DATAin%Cfol_stock(dn))/DATAin%Cfol_stock_unc(dn))**2
!        end do
!        scale_likelihood = scale_likelihood-0.5d0*(tot_exp/dble(DATAin%nCfol_stock))
!    endif

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
end module model_likelihood_module
