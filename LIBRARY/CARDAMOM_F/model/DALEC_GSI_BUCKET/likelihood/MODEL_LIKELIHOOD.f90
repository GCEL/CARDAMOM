
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
    call assess_EDC1(PARS,PI%npars,DATAin%meantemp, DATAin%meanrad,EDC1)

    ! next need to run the model itself
    call carbon_model(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT,DATAin%M_LAI,DATAin%M_NEE &
                     ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%nopars &
                     ,DATAin%nomet,DATAin%nopools,DATAin%nofluxes  &
                     ,DATAin%M_GPP)

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

    ! run the dalec model
    call carbon_model(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT,DATAin%M_LAI,DATAin%M_NEE &
                     ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%nopars &
                     ,DATAin%nomet,DATAin%nopools,DATAin%nofluxes  &
                     ,DATAin%M_GPP)

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
    print*,"Model Sanity Check Completed"

  end subroutine model_sanity_check
  !
  !------------------------------------------------------------------
  !
  subroutine assess_EDC1(PARS, npars, meantemp, meanrad, EDC1)

    use cardamom_structures, only: DATAin

    ! subroutine assessed the current parameter sets for passing ecological and
    ! steady state contraints (modified from Bloom et al., 2014).

    implicit none

    ! declare input variables
    integer, intent(in) :: npars ! number of parameters
    double precision, intent(out) :: EDC1    ! EDC1 flag
    double precision, dimension(npars), intent(in) :: PARS ! current parameter set
    double precision, intent(in) :: meantemp & ! mean temperature (k)
                                   ,meanrad    ! mean radiation (MJ.m-2.day-1)

    ! declare local variables
    integer :: n, DIAG
    double precision :: temp_response, tmp, avN

    ! set initial value
    EDC1 = 1d0
    DIAG = EDCD%DIAG

    ! set all EDCs to 1 (pass)
    EDCD%PASSFAIL(1:EDCD%nedc) = 1

    !
    ! begin checking EDCs
    !

    ! calculate temperature response of decomposition processes
    temp_response = exp(pars(10)*meantemp)
    avN = 10d0**pars(11)

    ! Straight forward GSI parameter bounds
    ! Temperature
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(15) < pars(14))) then
         EDC1 = 0d0 ; EDCD%PASSFAIL(1) = 0
    end if
    ! Photoperiod
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(24) < pars(16))) then
         EDC1 = 0d0 ; EDCD%PASSFAIL(2) = 0
    end if
    ! VPD / wSWP
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(26) < pars(25))) then
         EDC1 = 0d0 ; EDCD%PASSFAIL(3) = 0
    end if

    ! Photoperiod minimum cannot be substantiall less than the observed minimum day length
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(16) < minval(DATAin%MET(11,:))-14400d0)) then
         EDC1 = 0d0 ; EDCD%PASSFAIL(4) = 0
    end if
    ! Photoperiod maximum cannot be greater than the observed maximum day length
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(24) > maxval(DATAin%MET(11,:)))) then
         EDC1 = 0d0 ; EDCD%PASSFAIL(5) = 0
    end if

    ! VPD at which stress in at maximum should be no larger than max(VPDlag21) +
    ! 1500 Pa from the max VPD tolerated parameter
!    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(26) > maxval(DATAin%MET(12,:)+1500d0))) then
!         EDC1 = 0d0 ; EDCD%PASSFAIL(6) = 0
!    end if

    ! NUE and avN combination give a Vcmax equivalent.
    ! Kattge et al (2011) offers a prior of 3.4 - 30.7 gC/m2leaf/day.
    ! Here, to be cautious we will expand accepted range
    ! Thus CUE = NUE * avN -> 1.64 / 42.0
!    tmp = avN * pars(36)
!    if ((EDC1 == 1 .or. DIAG == 1) .and. (tmp > 42.0d0 .or. tmp < 1.64d0) ) then
!       EDC1 = 0d0 ; EDCD%PASSFAIL(1) = 0
!    endif
    ! Further constraint can be made by linking into LCA based on the range of
    ! max photosynthesis per gC leaf Kattge et al (2011) 0.2488 (0.041472 / 1.016064) gC/gC/day.
!    tmp = tmp / pars(17)
!    if ((EDC1 == 1 .or. DIAG == 1) .and. (tmp > 1.016064d0 .or. tmp < 0.041472d0) ) then
!       EDC1 = 0d0 ; EDCD%PASSFAIL(2) = 0
!    endif
    ! NOTE: that the above two constraints are both needed. Each independently
    ! does not capture all values which fall out of bounds in the other
    ! constraint. This means that we effectively through trait-databased
    ! information constrain three important C-cycle traits.

    ! Turnover of litter towards som (pars(1)*pars(8)) should be faster than turnover of som (pars(9))
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(9) > (pars(1)*pars(8)) ) then
       EDC1 = 0d0 ; EDCD%PASSFAIL(7) = 0
    endif

    ! turnover of litwood (pars(38)) should be slower than fine litter turnover pars(8)
    if ((EDC1 == 1 .or. DIAG == 1) .and. ( pars(38) > pars(8) ) ) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(8) = 0
    endif

    ! root turnover (pars(7)) should be greater than som turnover (pars(9)) at mean temperature
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(9)*temp_response) > pars(7)) then
       EDC1 = 0d0 ; EDCD%PASSFAIL(9) = 0
    endif

    ! replanting 30 = labile ; 31 = foliar ; 32 = roots ; 33 = wood
    ! initial    18 = labile ; 19 = foliar ; 20 = roots ; 21 = wood
    ! initial replanting labile must be consistent with available wood storage
    ! space. Labile storage cannot be greater than 12.5 % of the total ecosystem
    ! carbon stock.
    ! Gough et al (2009) Agricultural and Forest Meteorology. Avg 11, 12.5, 3 %
    ! (Max across species for branch, bole and coarse roots). Evidence that
    ! Branches accumulate labile C prior to bud burst from other areas.
    ! Wurth et al (2005) Oecologia, Clab 8 % of living biomass (DM) in tropical
    ! forest Richardson et al (2013), New Phytologist, Clab 2.24 +/- 0.44 % in
    ! temperate (max = 4.2 %, min = 1.8)
!    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(30) > ((pars(33)+pars(32))*0.125d0)) ) then
!        EDC1 = 0d0 ; EDCD%PASSFAIL(6) = 0
!    endif
    ! also apply to initial conditions
!    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(18) > ((pars(21)+pars(20))*0.125d0)) ) then
!        EDC1 = 0d0 ; EDCD%PASSFAIL(7) = 0
!    endif

    ! It is very unlikely that the initial coarse woody debris could be greater
    ! than the total of the wood and som pools. While coarse woody debris
    ! originates from the wood pool, high CWD can occur with low wood due to
    ! disturbance. Therefore we use the som pool as a conservative marker to
    ! stablise the assumption
!    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(37) > (pars(21) + pars(23))) then
!        EDC1 = 0d0 ; EDCD%PASSFAIL(7) = 0
!    endif

    ! initial replanting foliage and fine roots ratio must be consistent with
    ! ecological ranges. Because this is the initial condition and not the mean
    ! only the upper foliar:fine root bound is applied
!    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(32)/pars(31) < 0.04d0) ) then
!        EDC1 = 0d0 ; EDCD%PASSFAIL(8) = 0
!    endif

    ! CN ratio of leaf should be between 95CI of trait database values
    ! Kattge et al (2011) (12.39 < CN_foliar < 42.2).
    ! NOTE: this may be too (insufficiently) restrictive...as it is unclear how
    ! much more
    ! constrained a CN ratio of the whole canopy is compared to individual
    ! leaves (which have ranges upto ~100).
!    tmp = pars(17) / avN ! foliar C:N
!    if ((EDC1 == 1 .or. DIAG == 1) .and. (tmp > 42.2d0 .or. tmp < 12.39d0)) then
!       EDC1 = 0d0 ; EDCD%PASSFAIL(12) = 0
!    endif

    ! IMPLICIT Combustion completeness for foliage should be greater than soil
    ! IMPLICIT Combustion completeness for wood litter and fol+root litter should be greater than soil

    ! Combustion completeness for foliage should be greater than non-photosynthetic tissues
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(42) < pars(43)) then
       EDC1 = 0d0 ; EDCD%PASSFAIL(10) = 0
    endif
    ! Combustion completeness for non-photosynthetic tissue should be greater than soil
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(43) < pars(44)) then
       EDC1 = 0d0 ; EDCD%PASSFAIL(11) = 0
    endif
    ! Combustion completeness for foliar + fine root litter should be greater than wood litter
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(45) < pars(46)) then
       EDC1 = 0d0 ; EDCD%PASSFAIL(12) = 0
    endif
    ! Combustion completeness for wood litter should be greater than non-photosynthetic tissue
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(46) < pars(43)) then
       EDC1 = 0d0 ; EDCD%PASSFAIL(13) = 0
    endif

    ! --------------------------------------------------------------------
    ! could always add more / remove some

  end subroutine assess_EDC1
  !
  !------------------------------------------------------------------
  !
  subroutine assess_EDC2(npars,nomet,nofluxes,nopools,nodays,deltat &
                     ,parmax,pars,met,M_LAI,M_NEE,M_GPP,M_POOLS,M_FLUXES &
                     ,meantemp,EDC2)

    use cardamom_structures, only: DATAin
    use CARBON_MODEL_MOD, only: Rg_from_labile, wSWP_time,     &
                                harvest_residue_to_litter, &
                                harvest_residue_to_som,    &
                                harvest_residue_to_litwood,&
                                harvest_loss_labile,       &
                                harvest_loss_foliar,       &
                                harvest_loss_roots,        &
                                harvest_loss_wood,         &
                                harvest_loss_litter,       &
                                harvest_loss_litwood,      &
                                harvest_loss_som,          &
                                fire_emiss_labile,         &
                                fire_emiss_foliar,         &
                                fire_emiss_roots,          &
                                fire_emiss_wood,           &
                                fire_emiss_litter,         &
                                fire_emiss_litwood,        &
                                fire_emiss_som,            &
                                fire_litter_labile,        &
                                fire_litter_foliar,        &
                                fire_litter_roots,         &
                                fire_litter_wood,          &
                                fire_litter_litter,        &
                                fire_litter_litwood,       &
                                fire_litter_som,           &
                                fire_residue_to_litter,    &
                                fire_residue_to_litwood,   &
                                fire_residue_to_som

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
    logical :: found
    integer :: y, n, DIAG, no_years, nn, nnn, num_EDC, i, io_start, io_finish, &
               steps_per_year, steps_per_month
    double precision :: mean_pools(nopools), meangpp, sumgpp, sumnpp, &
                        tmp, tmp1, tmp2, tmp3, tmp4, tmp5, temp_response, &
                        hold, infi, Rs, dble_nodays, &
                        mean_step_size
    double precision, dimension(nodays) :: mean_ratio, resid_fol, resid_lab
    double precision, dimension(nopools) :: jan_mean_pools, jan_first_pools
    integer, dimension(nodays) :: hak ! variable to determine number of NaN in foliar residence time calculation
    double precision :: SSwood, SSlitwood, SSsom &
                       ,in_out_fol, in_out_lab, in_out_lit, in_out_litwood, in_out_som, in_out_root, in_out_wood &
                       ,in_lab, out_lab  &
                       ,in_fol, out_fol  &
                       ,in_root, out_root &
                       ,in_wood, out_wood &
                       ,in_lit, out_lit  &
                       ,in_litwood, out_litwood  &
                       ,in_som, out_som  &
                       ,in_out_lab_yr1  &
                       ,in_out_fol_yr1  &
                       ,in_out_root_yr1 &
                       ,in_out_wood_yr1 &
                       ,in_out_lit_yr1  &
                       ,in_out_litwood_yr1  &
                       ,in_out_som_yr1  &
                       ,in_out_lab_yr2  &
                       ,in_out_fol_yr2  &
                       ,in_out_root_yr2 &
                       ,in_out_wood_yr2 &
                       ,in_out_lit_yr2  &
                       ,in_out_litwood_yr2  &
                       ,in_out_som_yr2  &
                       ,torfol      & ! yearly average turnover
                       ,torlab      & !
                       ,sumlab_yr1      &
                       ,sumfol_yr1      &
                       ,sumroot_yr1     &
                       ,sumwood_yr1     &
                       ,sumlab_yr2      &
                       ,sumfol_yr2      &
                       ,sumroot_yr2     &
                       ,sumwood_yr2     &
                       ,sumrauto    &
                       ,sumlab      &
                       ,sumfol      &
                       ,sumroot     &
                       ,sumwood     &
                       ,sumlitwood  &
                       ,sumlit      &
                       ,sumsom      &
                       ,fNPP        & ! fraction of NPP to foliage
                       ,rNPP        & ! fraction of NPP to roots
                       ,wNPP        & ! fraction of NPP to wood
                       ,fauto       & ! fraction of GPP to autotrophic respiration
                       ,ffol        & ! fraction of GPP to foliage
                       ,froot         ! fraction of GPP to root

    ! Steady State Attractor:
    ! Log ratio difference between inputs and outputs of the system.
    double precision, parameter :: EQF1_5 = log(1.5d0), & ! 10.0 = order magnitude; 2 = double and half
                                   EQF2 = log(2d0),   & ! 10.0 = order magnitude; 2 = double and half
                                   EQF5 = log(5d0),   &
                                   EQF10 = log(10d0), &
                                   EQF15 = log(15d0), &
                                   EQF20 = log(20d0), &
                                    etol = 0.20d0 !0.10d0 !0.05d0

    ! initial value
    infi = 0d0 ; dble_nodays = dble(nodays)

    ! reset some flags needed for EDC control
    DIAG = EDCD%DIAG
    EDC2 = 1d0

    !!!!!!!!!!!!
    ! calculate residence times
    !!!!!!!!!!!!

    !
    ! Foliar turnover
    !

    ! update initial values
    hak = 0 ; resid_fol = 0d0
    ! calculate mean turnover rate for leaves
    resid_fol = (M_FLUXES(:,10)+harvest_loss_foliar+fire_emiss_foliar+fire_litter_foliar) &
              / M_POOLS(1:nodays,2)
    ! division by zero results in NaN plus obviously I can't have turned
    ! anything over if there was nothing to start out with...
    where ( M_POOLS(1:nodays,2) == 0d0 )
           hak = 1 ; resid_fol = 0d0
    end where
    ! mean fractional loss per day
    torfol = sum(resid_fol) / (dble_nodays-dble(sum(hak)))

    !
    ! Labile turnover
    !

    ! reset initial values
    hak = 0 ; resid_lab = 0d0
    ! calculate mean turnover rate for labile pool
    resid_lab = (M_FLUXES(:,8)+Rg_from_labile+harvest_loss_labile+fire_emiss_labile+fire_litter_labile) &
                        / M_POOLS(1:nodays,1)
    ! division by zero results in NaN plus obviously I can't have turned
    ! anything over if there was nothing to start out with...
    where ( M_POOLS(1:nodays,1) == 0d0 )
           hak = 1 ; resid_lab = 0d0
    end where
    ! mean fractional loss of labile per day
    torlab = sum(resid_lab) / (dble_nodays-dble(sum(hak)))

    !!!!!!!!!!!!
    ! calculate and update / adjust timing variables
    !!!!!!!!!!!!

    ! number of years in analysis, 0.002737851 = 1/365.25
    no_years = nint(sum(deltat)*0.002737851d0)
    ! number of time steps per year
    steps_per_year = nint(dble_nodays/dble(no_years))
    ! mean step size in days
    mean_step_size = sum(deltat) / dble_nodays

!    !calculate mean annual pool size for foliage
!    allocate(mean_annual_pools(no_years))
!    mean_annual_pools = 0.0
!    do y = 1, no_years
!       ! derive mean annual foliar pool
!       mean_annual_pools(y)=cal_mean_annual_pools(M_POOLS(1:(nodays+1)),y,deltat,nodays+1)
!    end do ! year loop

    !!!!!!!!!!!!
    ! Estimate mean January pool sizes for dynamics constraints
    !!!!!!!!!!!!

    ! number of time steps per month, 1/12 = 0.08333333
    steps_per_month = ceiling(dble(steps_per_year) * 0.08333333d0)

    ! Determine the mean January pool sizes
    jan_mean_pools = 0d0 ; jan_first_pools = 0d0 ! reset before averaging
    do n = 1, nopools
       jan_first_pools(n) = sum(M_POOLS(1:steps_per_month,n)) / dble(steps_per_month)
       do y = 1, no_years
          nn = 1 + (steps_per_year * (y - 1)) ; nnn = nn + (steps_per_month - 1)
          jan_mean_pools(n) = jan_mean_pools(n) + sum(M_POOLS(nn:nnn,n))
       end do
       jan_mean_pools(n) = jan_mean_pools(n) / dble(steps_per_month*no_years)
    end do

    !!!!!!!!!!!!
    ! calculate photosynthate / NPP allocations
    !!!!!!!!!!!!

    ! calculate sum fluxes
    sumgpp = sum(M_FLUXES(:,1))
    sumrauto = sum(M_FLUXES(:,3))
    sumlab = sum(M_FLUXES(:,5))
    sumfol = sum(M_FLUXES(:,8))
    sumroot = sum(M_FLUXES(:,6))
    sumwood = sum(M_FLUXES(:,7))

    ! initialise and then calculate mean gpp values
    fauto = sumrauto / sumgpp            ! i.e. Ra:GPP = 1-CUE
    sumnpp = (sumgpp - sumrauto)**(-1d0) ! NOTE: inverted here
    sumgpp = sumgpp**(-1d0)              ! NOTE: inverted here

    ! GPP allocation fractions
    ffol = sumfol * sumgpp
    froot = sumroot * sumgpp

    ! NPP allocations; note that because of possible labile accumulation this
    ! might not be equal to 1
    fNPP = sumfol * sumnpp
    wNPP = sumwood * sumnpp
    rNPP = sumroot * sumnpp

    ! derive mean pools
    do n = 1, nopools
       mean_pools(n) = sum(M_POOLS(1:nodays,n)) / dble_nodays
    end do

    !
    ! Begin EDCs here
    !

    ! GPP allocation to foliage and labile cannot be 5 orders of magnitude
    ! difference from GPP allocation to roots
    if ((EDC2 == 1 .or. DIAG == 1) .and. (ffol > (5d0*froot) .or. (ffol*5d0) < froot)) then
        EDC2 = 0d0 ; EDCD%PASSFAIL(14) = 0
    endif
    ! Restrict difference between root and foliar turnover to less than 5 fold
!    if ((EDC2 == 1 .or. DIAG == 1) .and. (torfol > pars(7)*5d0 .or. torfol*5d0 < pars(7) )) then
!        EDC2 = 0d0 ; EDCD%PASSFAIL(11) = 0
!    endif
    ! Restrict maximum leaf lifespan
    ! 0.0003422313 = (8 * 365.25)**-1
    if ((EDC2 == 1 .or. DIAG == 1) .and. torfol < 0.0003422313d0) then
        EDC2 = 0d0 ; EDCD%PASSFAIL(15) = 0
    endif

    ! Average turnover of foliage should not be less than wood (pars(6))
    if ((EDC2 == 1 .or. DIAG == 1) .and. torfol < pars(6) ) then
        EDC2 = 0d0 ; EDCD%PASSFAIL(16) = 0
    endif

    ! In contrast to the leaf longevity labile carbon stocks can be quite long
    ! lived, particularly in forests.
    ! Richardson et al (2015) New Phytologist, Clab residence time = 11 +/- 7.4 yrs (95CI = 18 yr)
    ! NOTE: 18 years = 0.0001521028 day-1
    !       11 years = 0.0002488955 day-1
    !        6 years = 0.0004563085 day-1
    if ((EDC2 == 1 .or. DIAG == 1) .and. torlab < 0.0002488955d0) then
        EDC2 = 0d0 ; EDCD%PASSFAIL(17) = 0
    endif

    ! Finally we would not expect that the mean labile stock is greater than
    ! 8 % of the total ecosystem carbon stock, as we need structure to store
    ! labile.
    ! Gough et al (2009) Agricultural and Forest Meteorology. Avg 11, 12.5, 3 %
    ! (Max across species for branch, bole and coarse roots). Provides evidence that
    ! branches accumulate labile C prior to bud burst from other areas.
    ! Wurth et al (2005) Oecologia, Clab 8 % of living biomass (DM) in tropical forest
    ! Richardson et al (2013), New Phytologist, Clab 2.24 +/- 0.44 % in temperate (max = 4.2 %)
    if (EDC2 == 1 .or. DIAG == 1) then
        if ((mean_pools(1) / (mean_pools(3) + mean_pools(4))) > 0.125d0) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(18) = 0
        endif
        if (maxval(M_POOLS(:,1) / (M_POOLS(:,3) + M_POOLS(:,4))) > 0.25d0) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(19) = 0
        endif
    endif ! EDC2 == 1 .or. DIAG == 1

    ! EDC 6
    ! ensure fine root : foliage ratio is between 0.1 and 0.45 (Albaugh et al
    ! 2004; Samuelson et al 2004; Vogel et al 2010; Akers et al 2013
    ! Duke ambient plots between 0.1 and 0.55
    ! Black et al 2009 Sitka Spruce chronosquence
    ! Q1 = 0.1278, median = 0.7488, mean = 1.0560 Q3 = 1.242
    ! lower CI = 0.04180938, upper CI = 4.06657167
    ! Field estimates tend to be made at growing season peaks, therefore we will
    ! consider the max(root)/max(fol) instead
!    if (EDC2 == 1 .or. DIAG == 1) then
!        mean_ratio(1) = maxval(M_POOLS(1:nodays,3)) / maxval(M_POOLS(1:nodays,2))
!        if ( mean_ratio(1) < 0.0418093d0 .or. mean_ratio(1) > 4.07d0 ) then
!            EDC2 = 0d0 ; EDCD%PASSFAIL(17) = 0
!        end if
!    endif !

    !
    ! EDC 14 - Fractional allocation to foliar biomass is well constrained
    ! across dominant ecosystem types (boreal -> temperate evergreen and
    ! deciduous -> tropical), therefore this information can be used to contrain the foliar pool
    ! further. Through control of the photosynthetically active component of the carbon
    ! balance we can enforce additional contraint on the remainder of the system.
    ! Luyssaert et al (2007)

    ! Limits on foliar allocation
    if ((EDC2 == 1 .or. DIAG == 1) .and. fNPP < 0.05d0) then
        EDC2 = 0d0 ; EDCD%PASSFAIL(20) = 0
    endif
    ! Limits on fine root allocation
    if ((EDC2 == 1 .or. DIAG == 1) .and. rNPP < 0.05d0) then
        EDC2 = 0d0 ; EDCD%PASSFAIL(21) = 0
    endif
    ! foliar restrictions
!    if ((EDC2 == 1 .or. DIAG == 1) .and. (fNPP < 0.1d0 .or. fNPP > 0.5d0)) then
!        EDC2 = 0d0 ; EDCD%PASSFAIL(22) = 0
!    endif
    ! for both roots and wood the NPP > 0.85 is added to prevent large labile
    ! pools being used to support growth that photosynthesis cannot provide over
    ! the long term.
!    if ((EDC2 == 1 .or. DIAG == 1) .and. (rNPP < 0.05d0 .or. rNPP > 0.85d0 .or. wNPP > 0.85d0)) then
!        EDC2 = 0d0 ; EDCD%PASSFAIL(19) = 0
!    endif
!    if ((EDC2 == 1 .or. DIAG == 1) .and. wNPP > 0.85d0) then
!        EDC2 = 0d0 ; EDCD%PASSFAIL(20) = 0
!    endif
    ! NOTE that within the current framework NPP is split between fol, root, wood and that remaining in labile.
    ! Thus fail conditions fNPP + rNPP + wNPP > 1.0 .or. fNPP + rNPP + wNPP < 0.95, i.e. sb(lNPP) cannot be > 0.05
!    tmp = 1d0 - rNPP - wNPP - fNPP
!    if ((EDC2 == 1 .or. DIAG == 1) .and. abs(tmp) > 0.05d0) then
!        EDC2 = 0d0 ; EDCD%PASSFAIL(21) = 0
!    endif

    ! Ra:GPP ratio is unlikely to be outside of 0.2 > Ra:GPP < 0.80
!    if ((EDC2 == 1 .or. DIAG == 1) .and. (fauto > 0.80d0 .or. fauto < 0.20d0) ) then
!        EDC2 = 0d0 ; EDCD%PASSFAIL(22) = 0
!    end if

    !!!!!!!!!
    ! Deal with ecosystem dynamics
    !!!!!!!!!

    ! this is a big set of arrays to run through so only do so when we have
    ! reached this point and still need them
    if (EDC2 == 1 .or. DIAG == 1) then

        ! Determine correct start and end points of the input / output assessment
        io_start = (steps_per_year*2) + 1 ; io_finish = nodays
        if (no_years < 3) io_start = 1

        ! calculate sum fluxes for the beginning of the timeseries
        sumlab_yr1 = sum(M_FLUXES(1:steps_per_year,5)) ; sumlab_yr2 = sum(M_FLUXES((steps_per_year+1):(steps_per_year*2),5))
        sumfol_yr1 = sum(M_FLUXES(1:steps_per_year,8)) ; sumfol_yr2 = sum(M_FLUXES((steps_per_year+1):(steps_per_year*2),8))
        sumroot_yr1 = sum(M_FLUXES(1:steps_per_year,6)) ; sumroot_yr2 = sum(M_FLUXES((steps_per_year+1):(steps_per_year*2),6))
        sumwood_yr1 = sum(M_FLUXES(1:steps_per_year,7)) ; sumwood_yr2 = sum(M_FLUXES((steps_per_year+1):(steps_per_year*2),7))

        ! Determine the input / output ratio

        ! Clabile
!        in_out_lab     = sumlab &
!                       / sum(M_FLUXES(:,8)+Rg_from_labile+fire_emiss_labile+fire_litter_labile+harvest_loss_labile)
        in_lab         = sum(M_FLUXES(io_start:io_finish,5))
        out_lab        = sum(M_FLUXES(io_start:io_finish,8) &
                            +Rg_from_labile(io_start:io_finish) &
                            +fire_emiss_labile(io_start:io_finish) &
                            +fire_litter_labile(io_start:io_finish) &
                            +harvest_loss_labile(io_start:io_finish))
        in_out_lab_yr1 = sumlab_yr1 &
                       / sum(M_FLUXES(1:steps_per_year,8) &
                            +Rg_from_labile(1:steps_per_year) &
                            +fire_emiss_labile(1:steps_per_year) &
                            +fire_litter_labile(1:steps_per_year) &
                            +harvest_loss_labile(1:steps_per_year))
        in_out_lab_yr2 = sumlab_yr2 &
                       / sum(M_FLUXES((steps_per_year+1):(steps_per_year*2),8) &
                            +Rg_from_labile((steps_per_year+1):(steps_per_year*2)) &
                            +fire_emiss_labile((steps_per_year+1):(steps_per_year*2)) &
                            +fire_litter_labile((steps_per_year+1):(steps_per_year*2)) &
                            +harvest_loss_labile((steps_per_year+1):(steps_per_year*2)))
        ! Cfoliage
!        in_out_fol  = sumfol  / sum(M_FLUXES(:,10)+fire_emiss_foliar+fire_litter_foliar+harvest_loss_foliar)
        in_fol      = sum(M_FLUXES(io_start:io_finish,8))
        out_fol     = sum(M_FLUXES(io_start:io_finish,10) &
                         +fire_emiss_foliar(io_start:io_finish) &
                         +fire_litter_foliar(io_start:io_finish) &
                         +harvest_loss_foliar(io_start:io_finish))
        in_out_fol_yr1  = sumfol_yr1  / sum(M_FLUXES(1:steps_per_year,10) &
                                           +fire_emiss_foliar(1:steps_per_year) &
                                           +fire_litter_foliar(1:steps_per_year) &
                                           +harvest_loss_foliar(1:steps_per_year))
        in_out_fol_yr2  = sumfol_yr2  / sum(M_FLUXES((steps_per_year+1):(steps_per_year*2),10) &
                                           +fire_emiss_foliar((steps_per_year+1):(steps_per_year*2)) &
                                           +fire_litter_foliar((steps_per_year+1):(steps_per_year*2)) &
                                           +harvest_loss_foliar((steps_per_year+1):(steps_per_year*2)))
        ! Croot
!        in_out_root = sumroot / sum(M_FLUXES(:,12)+fire_emiss_roots+fire_litter_roots+harvest_loss_roots)
        in_root     = sum(M_FLUXES(io_start:io_finish,6))
        out_root    = sum(M_FLUXES(io_start:io_finish,12) &
                         +fire_emiss_roots(io_start:io_finish) &
                         +fire_litter_roots(io_start:io_finish) &
                         +harvest_loss_roots(io_start:io_finish))
        in_out_root_yr1 = sumroot_yr1 / sum(M_FLUXES(1:steps_per_year,12) &
                                           +fire_emiss_roots(1:steps_per_year) &
                                           +fire_litter_roots(1:steps_per_year) &
                                           +harvest_loss_roots(1:steps_per_year))
        in_out_root_yr2 = sumroot_yr2 / sum(M_FLUXES((steps_per_year+1):(steps_per_year*2),12) &
                                           +fire_emiss_roots((steps_per_year+1):(steps_per_year*2)) &
                                           +fire_litter_roots((steps_per_year+1):(steps_per_year*2)) &
                                           +harvest_loss_roots((steps_per_year+1):(steps_per_year*2)))
        ! Cwood
!        in_out_wood = sumwood / sum(M_FLUXES(:,11)+fire_emiss_wood+fire_litter_wood+harvest_loss_wood)
        in_wood     = sum(M_FLUXES(io_start:io_finish,7))
        out_wood    = sum(M_FLUXES(io_start:io_finish,11) &
                         +fire_emiss_wood(io_start:io_finish) &
                         +fire_litter_wood(io_start:io_finish) &
                         +harvest_loss_wood(io_start:io_finish))
        in_out_wood_yr1 = sumwood_yr1 / sum(M_FLUXES(1:steps_per_year,11) &
                                           +fire_emiss_wood(1:steps_per_year) &
                                           +fire_litter_wood(1:steps_per_year) &
                                           +harvest_loss_wood(1:steps_per_year))
        in_out_wood_yr2 = sumwood_yr2 / sum(M_FLUXES((steps_per_year+1):(steps_per_year*2),11) &
                                           +fire_emiss_wood((steps_per_year+1):(steps_per_year*2)) &
                                           +fire_litter_wood((steps_per_year+1):(steps_per_year*2)) &
                                           +harvest_loss_wood((steps_per_year+1):(steps_per_year*2)))
        ! Clitter
!        in_out_lit = sum(M_FLUXES(:,10) &
!                        +M_FLUXES(:,12) &
!                        +fire_residue_to_litter &
!                        +harvest_residue_to_litter) &
!                   / sum(M_FLUXES(:,13)+M_FLUXES(:,15)+fire_emiss_litter+fire_litter_litter+harvest_loss_litter)
        in_lit     = sum(M_FLUXES(io_start:io_finish,10) &
                        +M_FLUXES(io_start:io_finish,12) &
                        +fire_residue_to_litter(io_start:io_finish) &
                        +harvest_residue_to_litter(io_start:io_finish))
        out_lit    = sum(M_FLUXES(io_start:io_finish,13)+M_FLUXES(io_start:io_finish,15) &
                        +fire_emiss_litter(io_start:io_finish)+fire_litter_litter(io_start:io_finish) &
                        +harvest_loss_litter(io_start:io_finish))
        in_out_lit_yr1 = sum(M_FLUXES(1:steps_per_year,10) &
                            +M_FLUXES(1:steps_per_year,12) &
                            +fire_residue_to_litter(1:steps_per_year) &
                            +harvest_residue_to_litter(1:steps_per_year)) &
                       / sum(M_FLUXES(1:steps_per_year,13) &
                            +M_FLUXES(1:steps_per_year,15) &
                            +fire_emiss_litter(1:steps_per_year) &
                            +fire_litter_litter(1:steps_per_year) &
                            +harvest_loss_litter(1:steps_per_year))
        in_out_lit_yr2 = sum(M_FLUXES((steps_per_year+1):(steps_per_year*2),10) &
                            +M_FLUXES((steps_per_year+1):(steps_per_year*2),12) &
                            +fire_residue_to_litter((steps_per_year+1):(steps_per_year*2)) &
                            +harvest_residue_to_litter((steps_per_year+1):(steps_per_year*2))) &
                       / sum(M_FLUXES((steps_per_year+1):(steps_per_year*2),13) &
                            +M_FLUXES((steps_per_year+1):(steps_per_year*2),15) &
                            +fire_emiss_litter((steps_per_year+1):(steps_per_year*2)) &
                            +fire_litter_litter((steps_per_year+1):(steps_per_year*2)) &
                            +harvest_loss_litter((steps_per_year+1):(steps_per_year*2)))
        ! Csom
!        in_out_som = sum(M_FLUXES(:,15)+M_FLUXES(:,20)+fire_residue_to_som+harvest_residue_to_som) &
!                   / sum(M_FLUXES(:,14)+fire_emiss_som+fire_litter_som+harvest_loss_som)
        in_som     = sum(M_FLUXES(io_start:io_finish,15)+M_FLUXES(io_start:io_finish,20) &
                        +fire_residue_to_som(io_start:io_finish)+harvest_residue_to_som(io_start:io_finish))
        out_som    = sum(M_FLUXES(io_start:io_finish,14) &
                        +fire_emiss_som(io_start:io_finish) &
                        +fire_litter_som(io_start:io_finish) &
                        +harvest_loss_som(io_start:io_finish))
        in_out_som_yr1 = sum(M_FLUXES(1:steps_per_year,15)+ &
                             M_FLUXES(1:steps_per_year,20)+ &
                             fire_residue_to_som(1:steps_per_year)+ &
                             harvest_residue_to_som(1:steps_per_year)) &
                       / sum(M_FLUXES(1:steps_per_year,14) &
                            +fire_emiss_som(1:steps_per_year) &
                            +fire_litter_som(1:steps_per_year) &
                            +harvest_loss_som(1:steps_per_year))
        in_out_som_yr2 = sum(M_FLUXES((steps_per_year+1):(steps_per_year*2),15)+ &
                             M_FLUXES((steps_per_year+1):(steps_per_year*2),20)+ &
                             fire_residue_to_som((steps_per_year+1):(steps_per_year*2))+ &
                             harvest_residue_to_som((steps_per_year+1):(steps_per_year*2))) &
                       / sum(M_FLUXES((steps_per_year+1):(steps_per_year*2),14) &
                            +fire_emiss_som((steps_per_year+1):(steps_per_year*2)) &
                            +fire_litter_som((steps_per_year+1):(steps_per_year*2)) &
                            +harvest_loss_som((steps_per_year+1):(steps_per_year*2)))
        ! Clitwood
!        in_out_litwood = sum(M_FLUXES(:,11)+fire_residue_to_litwood+harvest_residue_to_litwood) &
!                       / sum(M_FLUXES(:,20)+M_FLUXES(:,4)+fire_emiss_litwood+fire_litter_litwood+harvest_loss_litwood)
        in_litwood     = sum(M_FLUXES(io_start:io_finish,11) &
                            +fire_residue_to_litwood(io_start:io_finish) &
                            +harvest_residue_to_litwood(io_start:io_finish))
        out_litwood    = sum(M_FLUXES(io_start:io_finish,20) &
                            +M_FLUXES(io_start:io_finish,4) &
                            +fire_emiss_litwood(io_start:io_finish) &
                            +fire_litter_litwood(io_start:io_finish) &
                            +harvest_loss_litwood(io_start:io_finish))
        in_out_litwood_yr1 = sum(M_FLUXES(1:steps_per_year,11) &
                                +fire_residue_to_litwood(1:steps_per_year) &
                                +harvest_residue_to_litwood(1:steps_per_year)) &
                           / sum(M_FLUXES(1:steps_per_year,20) &
                                +M_FLUXES(1:steps_per_year,4) &
                                +fire_emiss_litwood(1:steps_per_year) &
                                +fire_litter_litwood(1:steps_per_year) &
                                +harvest_loss_litwood(1:steps_per_year))
        in_out_litwood_yr2 = sum(M_FLUXES((steps_per_year+1):(steps_per_year*2),11) &
                                +fire_residue_to_litwood((steps_per_year+1):(steps_per_year*2)) &
                                +harvest_residue_to_litwood((steps_per_year+1):(steps_per_year*2))) &
                           / sum(M_FLUXES((steps_per_year+1):(steps_per_year*2),20) &
                                +M_FLUXES((steps_per_year+1):(steps_per_year*2),4) &
                                +fire_emiss_litwood((steps_per_year+1):(steps_per_year*2)) &
                                +fire_litter_litwood((steps_per_year+1):(steps_per_year*2)) &
                                +harvest_loss_litwood((steps_per_year+1):(steps_per_year*2)))

        ! Assess pool dynamics relative to their own steady state attractors
        ! Based on Bloom et al (2016), PNAS. Combination of the in/out ratio and
        ! the ratio of pool size change combines the dynamics with pool
        ! magnitude correction. i.e. larger accumulations are allowable for
        ! small pools as is sensible but more restricted in larger pools.
        ! The etol comparison assesses the exponential dynamics, while
        ! comparison with EQF assesses the steady state attractor

        ! Labile
!        Rs = in_out_lab * (jan_mean_pools(1) / jan_first_pools(1))
!        if (abs(Rs-in_out_lab) > 0.1d0 .or. abs(log(in_out_lab)) > EQF10) then
        if (abs(log(in_out_lab_yr1)) - abs(log(in_out_lab_yr2)) > etol .or. &
            abs(log(in_lab/out_lab)) > EQF2) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(22) = 0
        end if

        ! Foliage
!        Rs = in_out_fol * (jan_mean_pools(2) / jan_first_pools(2))
!        if (abs(Rs-in_out_fol) > 0.1d0 .or. abs(log(in_out_fol)) > EQF10) then
        if (abs(log(in_out_fol_yr1)) - abs(log(in_out_fol_yr2)) > etol .or. &
            abs(log(in_fol/out_fol)) > EQF2) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(23) = 0
        end if

        ! Fine roots
!        Rs = in_out_root * (jan_mean_pools(3) / jan_first_pools(3))
!        if (abs(Rs-in_out_root) > 0.1d0 .or. abs(log(in_out_root)) > EQF10) then
        if (abs(log(in_out_root_yr1)) - abs(log(in_out_root_yr2)) > etol .or. &
            abs(log(in_root/out_root)) > EQF2) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(24) = 0
        end if

        ! Wood
!        Rs = in_out_wood * (jan_mean_pools(4) / jan_first_pools(4))
!        if (abs(Rs-in_out_wood) > 0.1d0 .or. abs(log(in_out_wood)) > EQF10) then
        if (abs(log(in_out_wood_yr1)) - abs(log(in_out_wood_yr2)) > etol .or. &
            abs(log(in_wood/out_wood)) > EQF2) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(25) = 0
        end if

        ! Foliage and root litter
!        Rs = in_out_lit * (jan_mean_pools(5) / jan_first_pools(5))
!        if (abs(Rs-in_out_lit) > 0.1d0 .or. abs(log(in_out_lit)) > EQF10) then
        if (abs(log(in_out_lit_yr1)) - abs(log(in_out_lit_yr2)) > etol .or. &
            abs(log(in_lit/out_lit)) > EQF2) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(26) = 0
        end if

        ! Soil organic matter
!        Rs = in_out_som * (jan_mean_pools(6) / jan_first_pools(6))
!        if (abs(Rs-in_out_som) > 0.1d0 .or. abs(log(in_out_som)) > EQF10) then
        if (abs(log(in_out_som_yr1)) - abs(log(in_out_som_yr2)) > etol .or. &
            abs(log(in_som/out_som)) > EQF2) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(27) = 0
        end if

        ! Coarse+fine woody debris
!        Rs = in_out_litwood * (jan_mean_pools(7) / jan_first_pools(7))
!        if (abs(Rs-in_out_litwood) > 0.1d0 .or. abs(log(in_out_litwood)) > EQF10) then
        if (abs(log(in_out_litwood_yr1)) - abs(log(in_out_litwood_yr2)) > etol .or. &
            abs(log(in_litwood/out_litwood)) > EQF2) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(28) = 0
        end if

        ! Determine the steady state estimate of wood (gC/m2)
        SSwood = (in_wood/out_wood) * jan_mean_pools(4)
        ! Based on the wood SS (gC/m2) and the sum fractional loss per day determine the mean input to litwood...
        SSlitwood = SSwood * (out_wood/jan_mean_pools(4))
        ! ...then estimate the actual steady state wood litter
        SSlitwood = (SSlitwood/out_litwood) * jan_mean_pools(7)
        ! Steady state of som requires accounting for foliar, fine root and wood litter inputs
        ! and adjusting for the litwood input already included
        SSsom = in_som - sum(M_FLUXES(io_start:io_finish,20))
        ! Now repeat the process as done for litwood to estimate the inputs,
        ! adjusting for the fraction of litwood output which is respired not decomposed
        SSsom = SSsom + (SSlitwood * (out_litwood/jan_mean_pools(7)) * pars(1))
        ! Accounting for losses and scaling to SSsom
        SSsom = (SSsom / out_som) * jan_mean_pools(6)
        ! It is reasonable to assume that the steady state for woody litter
        ! should be ~ less than half that of woody biomass...
        if (SSlitwood / SSwood > 0.60d0  ) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(29) = 0
        end if
        ! ... and less than soil organic matter
        if ( SSsom < SSlitwood ) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(30) = 0
        end if
        ! It is reasonable to assume that the steady state for woody litter
!        ! should be ~ less than half that of woody biomass...
!        if ((in_out_litwood * jan_mean_pools(7)) / (in_out_wood * jan_mean_pools(4)) > 0.60d0  ) then
!            EDC2 = 0d0 ; EDCD%PASSFAIL(28) = 0
!        end if
!        ! ... and less than soil organic matter
!        if ( (in_out_som * jan_mean_pools(6)) < (in_out_litwood * jan_mean_pools(7)) ) then
!            EDC2 = 0d0 ; EDCD%PASSFAIL(29) = 0
!        end if
!        if (pars(37) / (in_out_wood * jan_mean_pools(4)) > 0.60d0  ) then
!            EDC2 = 0d0 ; EDCD%PASSFAIL(23) = 0
!        end if

    endif ! doing the big arrays then?

    !
    ! EDCs done, below are additional fault detection conditions
    !

    ! The maximum value for GPP must be greater than 0, 0.001 to guard against precision values
    if ((EDC2 == 1 .or. DIAG == 1) .and. maxval(M_GPP) < 0.001d0) then
        EDC2 = 0d0 ; EDCD%PASSFAIL(31) = 0
    end if

    ! additional faults can be stored in locations > 55 of the PASSFAIL array

    ! ensure minimum pool values are >= 0 and /= NaN
    if (EDC2 == 1 .or. DIAG == 1) then

       do n = 1, nopools
          if (minval(M_POOLS(1:nodays,n)) < 0d0 .or. maxval(abs(M_POOLS(1:nodays,n))) == abs(log(infi)) .or. &
              minval(M_POOLS(1:nodays,n)) /= minval(M_POOLS(1:nodays,n))) then
              EDC2 = 0d0 ; EDCD%PASSFAIL(55+n) = 0
          endif
       end do

       do n = 1, nofluxes
          if (maxval(abs(M_FLUXES(:,n))) == abs(log(infi)) .or. &
              minval(M_FLUXES(:,n)) /= minval(M_FLUXES(:,n))) then
              EDC2 = 0d0 ; EDCD%PASSFAIL(55+nopools+n) = 0
          endif
       end do

    end if ! min pool assessment

  end subroutine assess_EDC2
  !
  !------------------------------------------------------------------
  !
  subroutine UK_forestry_commission_growth_curves(target_living_C,max_location)
    use cardamom_structures, only: DATAin

    ! subroutine uses PFT and yield classification to generate an estimate of
    ! expected living C accumulated at a given age. Equation generated Mg.ha-1
    ! we need to correct this to gC.m-2 for the model

    implicit none

    ! declare input / output variables
    double precision, intent(out) :: target_living_C(2) ! (gC.m-2)
    integer, intent(in) :: max_location(1) ! additional years from initial

    ! local variables
    double precision :: adjusted_age, tmp1(2),tmp2(2)

    ! calculate adjusted age from initial conditions to max point
    adjusted_age=DATAin%age+max_location(1)

    ! set initial value for output
    target_living_C = 0d0

    ! loop through to get the minimum (1) and maximum estimates (2)
    ! which will be passed back to the model

    ! if we have an age (therefore it is a forest but we don't know even
    ! if it is evergreen or deciduos) we will assume the most generous
    ! range of values possible

    ! broadleaf
    tmp1(1) = 2.07956043460835d-05*adjusted_age**3d0 &
            + (-0.0141108480550955d0)*adjusted_age**2d0 &
            + 3.14928740556523d0*adjusted_age
    tmp1(2) = 0.000156065120683174d0*adjusted_age**3d0 &
            + (-0.0629544794948499d0)*adjusted_age**2d0 &
            + 8.30163202577001d0*adjusted_age
    ! evergreen
    tmp2(1) =  8.8519973125961d-06*adjusted_age**3d0 &
            + (-0.00822909089061558d0)*adjusted_age**2d0 &
            + 1.98952585135788d0*adjusted_age
    tmp2(2) = 0.00014916728414466d0*adjusted_age**3d0 &
            + (-0.0662815983372182d0)*adjusted_age**2d0 &
            + 9.55519207729034d0*adjusted_age
    ! work out which to use
    ! use smallest
    if (tmp1(1) < tmp2(1)) then
        target_living_C(1) = tmp1(1)*0.70d0
    else
        target_living_C(1) = tmp2(1)*0.70d0
    endif
    ! use biggest
    if (tmp1(2) > tmp2(2)) then
        target_living_C(2) = tmp1(2)*1.30d0
    else
        target_living_C(2) = tmp2(2)*1.30d0
    endif

    ! correct units from MgC.ha-1 to gC.m-2
    target_living_C = target_living_C*100

  end subroutine UK_forestry_commission_growth_curves
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

    ! run the dalec model
    call carbon_model(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT,DATAin%M_LAI,DATAin%M_NEE &
                     ,DATAin%M_FLUXES,DATAin%M_POOLS,DATAin%nopars &
                     ,DATAin%nomet,DATAin%nopools,DATAin%nofluxes  &
                     ,DATAin%M_GPP)

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
    use carbon_model_mod, only: layer_thickness, &
                                fire_emiss_wood, &
                               fire_litter_wood

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

    ! NBE Log-likelihood
    if (DATAin%nnbe > 0) then
       tot_exp = sum((((DATAin%M_NEE(DATAin%nbepts(1:DATAin%nnbe))+DATAin%M_FLUXES(DATAin%nbepts(1:DATAin%nnbe),17)) &
                       -DATAin%NBE(DATAin%nbepts(1:DATAin%nnbe))) &
                       /DATAin%NBE_unc(DATAin%nbepts(1:DATAin%nnbe)))**2)
       likelihood = likelihood-tot_exp
    endif

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

    ! Fire Log-likelihood
    if (DATAin%nFire > 0) then
       tot_exp = sum(((DATAin%M_FLUXES(DATAin%Firepts(1:DATAin%nFire),17)-DATAin%Fire(DATAin%Firepts(1:DATAin%nFire))) &
                       /DATAin%Fire_unc(DATAin%Firepts(1:DATAin%nFire)))**2)
       likelihood = likelihood-tot_exp
    endif

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

    ! Cagb log-likelihood
    if (DATAin%nCagb_stock > 0) then
       ! Create vector of (agb_t0 + agb_t1) * 0.5
       mid_state = ( DATAin%M_POOLS(1:DATAin%nodays,4) + DATAin%M_POOLS(2:(DATAin%nodays+1),4) ) &
                 * 0.5d0 * (1d0-pars(29))
       ! Vectorised version of loop to estimate cost function
       tot_exp = sum(( (mid_state(DATAin%Cagb_stockpts(1:DATAin%nCagb_stock)) &
                       -DATAin%Cagb_stock(DATAin%Cagb_stockpts(1:DATAin%nCagb_stock)))&
                     / DATAin%Cagb_stock_unc(DATAin%Cagb_stockpts(1:DATAin%nCagb_stock)))**2)
       ! Combine with existing likelihood estimate
       likelihood = likelihood-tot_exp
    endif

    ! Ccoarseroot log-likelihood
    if (DATAin%nCcoarseroot_stock > 0) then
       ! Create vector of (coarseroot_t0 + coarseroot_t1) * 0.5
       mid_state = ( DATAin%M_POOLS(1:DATAin%nodays,4) + DATAin%M_POOLS(2:(DATAin%nodays+1),4) ) &
                 * 0.5d0 * pars(29)
       ! Vectorised version of loop to estimate cost function
       tot_exp = sum(( (mid_state(DATAin%Ccoarseroot_stockpts(1:DATAin%nCcoarseroot_stock)) &
                       -DATAin%Ccoarseroot_stock(DATAin%Ccoarseroot_stockpts(1:DATAin%nCcoarseroot_stock)))&
                     / DATAin%Ccoarseroot_stock_unc(DATAin%Ccoarseroot_stockpts(1:DATAin%nCcoarseroot_stock)))**2)
       ! Combine with existing likelihood estimate
       likelihood = likelihood-tot_exp
    endif

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

    ! Estimate the biological steady state attractor on the wood pool.
    ! NOTE: this arrangement explicitly neglects the impact of disturbance on
    ! residence time (i.e. no fire and biomass removal)
    if (DATAin%otherpriors(5) > -9998) then
        ! Estimate the mean annual input to the wood pool (gC.m-2.day-1) and
        ! remove the day-1 by multiplying by residence time (day)
        !tot_exp = (sum(DATAin%M_FLUXES(:,7)) / dble(DATAin%nodays)) * (pars(6) ** (-1d0))
        input = sum(DATAin%M_FLUXES(:,7))
        output = sum(DATAin%M_POOLS(:,4) / (DATAin%M_FLUXES(:,11)+fire_emiss_wood+fire_litter_wood))
        tot_exp = (input/dble(DATAin%nodays)) * (output/dble(DATAin%nodays))
        tot_exp =  DATAin%otherpriorweight(5)* ((tot_exp-DATAin%otherpriors(5))/DATAin%otherpriorunc(5))**2
        likelihood = likelihood-tot_exp
    endif

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
    use carbon_model_mod, only: layer_thickness, &
                                fire_emiss_wood, &
                               fire_litter_wood

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

    ! NBE Log-likelihood
    if (DATAin%nnbe > 0) then
       tot_exp = sum((((DATAin%M_NEE(DATAin%nbepts(1:DATAin%nnbe))+DATAin%M_FLUXES(DATAin%nbepts(1:DATAin%nnbe),17)) &
                       -DATAin%NBE(DATAin%nbepts(1:DATAin%nnbe))) &
                       /DATAin%NBE_unc(DATAin%nbepts(1:DATAin%nnbe)))**2)
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nnbe))
    endif

    ! Evap Log-likelihood
    if (DATAin%nEvap > 0) then
       tot_exp = sum(((DATAin%M_FLUXES(DATAin%Evappts(1:DATAin%nEvap),19)-DATAin%Evap(DATAin%Evappts(1:DATAin%nEvap))) &
                       /DATAin%Evap_unc(DATAin%evappts(1:DATAin%nEvap)))**2)
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nEvap))
    endif

    ! Fire Log-likelihood
    if (DATAin%nFire > 0) then
       tot_exp = sum(((DATAin%M_FLUXES(DATAin%Firepts(1:DATAin%nFire),17)-DATAin%Fire(DATAin%Firepts(1:DATAin%nFire))) &
                       /DATAin%Fire_unc(DATAin%Firepts(1:DATAin%nFire)))**2)
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nFire))
    endif

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

    ! Cagb log-likelihood
    if (DATAin%nCagb_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCagb_stock
         dn = DATAin%Cagb_stockpts(n)
         ! remove coarse root fraction from wood (pars29)
         tmp_var = DATAin%M_POOLS(dn,4)-(DATAin%M_POOLS(dn,4)*pars(29))
         tot_exp = tot_exp+((tmp_var-DATAin%Cagb_stock(dn))/DATAin%Cagb_stock_unc(dn))**2
       end do
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nCagb_stock))
    endif

    ! Ccoarseroot log-likelihood
    if (DATAin%nCcoarseroot_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCcoarseroot_stock
         dn = DATAin%Ccoarseroot_stockpts(n)
         ! extract coarse root component from wood only
         tmp_var = DATAin%M_POOLS(dn,4)*pars(29)
         tot_exp = tot_exp+((tmp_var-DATAin%Ccoarseroot_stock(dn)) / DATAin%Ccoarseroot_stock_unc(dn))**2
       end do
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nCcoarseroot_stock))
    endif

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

    ! Estimate the biological steady state attractor on the wood pool.
    ! NOTE: this arrangement explicitly neglects the impact of disturbance on
    ! residence time (i.e. no fire and biomass removal)
    if (DATAin%otherpriors(5) > -9998) then
        ! Estimate the mean annual input to the wood pool (gC.m-2.day-1) and
        ! remove the day-1 by multiplying by residence time (day)
        !tot_exp = (sum(DATAin%M_FLUXES(:,7)) / dble(DATAin%nodays)) * (pars(6) ** (-1d0))
        input = sum(DATAin%M_FLUXES(:,7))
        output = sum(DATAin%M_POOLS(:,4) / (DATAin%M_FLUXES(:,11)+fire_emiss_wood+fire_litter_wood))
        tot_exp = (input/dble(DATAin%nodays)) * (output/dble(DATAin%nodays))
        tot_exp =  DATAin%otherpriorweight(5)* ((tot_exp-DATAin%otherpriors(5))/DATAin%otherpriorunc(5))**2
        scale_likelihood = scale_likelihood-tot_exp
    endif

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
end module model_likelihood_module
