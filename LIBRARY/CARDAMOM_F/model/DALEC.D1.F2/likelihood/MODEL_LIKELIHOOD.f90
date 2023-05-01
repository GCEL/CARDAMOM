
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
!    ML_obs_out = -0.5d0*(tot_exp*10d0)*DATAin%EDC
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
    ML_prior_out = likelihood_p(PI%npars,DATAin%parpriors,DATAin%parpriorunc,PARS)
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

  end subroutine model_sanity_check
  !
  !------------------------------------------------------------------
  !
  subroutine assess_EDC1(PARS, npars, meantemp, meanrad, EDC1)

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
                       ,froot & ! Fraction of GPP to root
                       ,fwood   ! Fraction of GPP to wood

    double precision :: torfol ! yearly leaf loss fraction

    ! set initial value
    EDC1 = 1
    DIAG = EDCD%DIAG

    ! estimate GPP allocation fractions
    fauto = pars(2)
    ffol = (1d0-fauto)*pars(3)
    froot = (1d0-fauto-ffol)*pars(4)
    fwood = 1d0-fauto-ffol-froot

    ! convert leaf life span in year to fraction per day
    torfol = 1d0/(pars(5)*365.25d0)

    ! set all EDCs to 1 (pass)
    EDCD%nedc = 100
    EDCD%PASSFAIL(1:EDCD%nedc) = 1

    !
    ! begin checking EDCs
    !

    ! Turnover of litter faster than turnover of som
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(9) > pars(8))) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(1) = 0
    endif

    ! litter2som greater than som to atm rate
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(1) < pars(9))) then
       EDC1 = 0d0 ; EDCD%PASSFAIL(2) = 0
    endif

    ! turnover of foliage faster than turnover of wood
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(6) > torfol) then
       EDC1 = 0d0 ; EDCD%PASSFAIL(3) = 0
    end if

    ! root turnover greater than som turnover at mean temperature
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(7) < (pars(9)*exp(pars(10)*meantemp)))) then
       EDC1 = 0d0 ; EDCD%PASSFAIL(4) = 0
    endif

    ! GPP allocation to foliage cannot be 5 orders of magnitude
    ! difference from GPP allocation to roots
    if ((EDC1 == 1 .or. DIAG == 1) .and. (ffol > (5d0*froot) .or. (ffol*5d0) < froot)) then
       EDC1 = 0d0 ; EDCD%PASSFAIL(5) = 0
    endif

    ! could always add more / remove some

  end subroutine assess_EDC1
  !
  !------------------------------------------------------------------
  !
  subroutine assess_EDC2(npars,nomet,nofluxes,nopools,nodays,deltat &
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
    integer :: n, nn, nnn, DIAG, no_years, y, PEDC, steps_per_year, steps_per_month
    double precision :: jan_first_pools(nopools), jan_mean_pools(nopools), mean_pools(nopools), EQF, etol
    double precision :: fauto & ! Fractions of GPP to autotrophic respiration
                       ,ffol  & ! Fraction of GPP to foliage
                       ,froot & ! Fraction of GPP to root
                       ,fwood   ! Fraction of GPP to wood

    !JFE - 27/06/2018 newly defined variables for updated EDCs
    double precision :: FT(nofluxes), Fin(nopools), Fout(nopools), Rm, Rs
    integer :: nd, fl

    ! update initial values
    DIAG = EDCD%DIAG
    EDC2 = 1

    ! estimate GPP allocation fractions
    fauto = pars(2)
    ffol = (1d0-fauto)*pars(3)
    froot = (1d0-fauto-ffol)*pars(4)
    fwood = 1d0-fauto-ffol-froot

    ! derive mean pools
    do n = 1, nopools
       mean_pools(n) = cal_mean_pools(M_POOLS,n,nodays+1,nopools)
    end do

    ! number of years in analysis
    no_years = nint(sum(deltat)/365.25d0)
    ! number of time steps per year
    steps_per_year = nodays/no_years
    ! number of time steps per month
    steps_per_month = ceiling(dble(steps_per_year) / 12d0)

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

    !
    ! Begin EDCs here
    !

    ! EDC 6
    ! ensure ratio between Cfoilar and Croot is less than 5
    if ((EDC2 == 1 .or. DIAG == 1) .and. &
        (mean_pools(1) > (mean_pools(2)*5d0) .or. (mean_pools(1)*5d0) < mean_pools(2)) ) then
        EDC2 = 0d0 ; EDCD%PASSFAIL(6) = 0
    end if

    ! Equilibrium factor (in comparison with initial conditions)
    EQF = 10d0 ! JFE replaced 10 by 2 - 27/06/2018
    ! Pool exponential decay tolerance
    etol = 0.1d0

    ! Calculate total flux for the whole simulation period
    do fl = 1, nofluxes
       FT(fl) = sum(M_FLUXES(1:nodays,fl)*deltat(1:nodays))
    end do

    ! get total in and out for each pool
    ! foliar
    Fin(1)  = FT(4)
    Fout(1) = FT(10)+FT(18)+FT(23)+FT(28)+FT(33)
    ! root
    Fin(2)  = FT(6)
    Fout(2) = FT(12)+FT(19)+FT(24)+FT(29)+FT(34)
    ! wood
    Fin(3)  = FT(7)
    Fout(3) = FT(11)+FT(20)+FT(25)+FT(30)+FT(35)
    ! litter
    Fin(4)  = FT(10)+FT(12)+FT(23)+FT(24)+FT(33)+FT(34)
    Fout(4) = FT(13)+FT(15)+FT(21)+FT(26)+FT(31)
    ! som
    Fin(5)  = FT(11)+FT(15)+FT(25)+FT(26)+FT(35)
    Fout(5) = FT(14)+FT(22)+FT(32)

    ! Iterate through C pools to determine whether they have their ratio of
    ! input and outputs are outside of steady state approximation.
    ! See Bloom et al., 2016 PNAS for details

    ! iterate to check whether Fin/Fout is within EQF limits
    do n = 1, nopools
       Rm = Fin(n)/Fout(n)
       Rs = Rm * (jan_mean_pools(n) / jan_first_pools(n))
       if ((EDC2 == 1 .or. DIAG == 1) .and. abs(log(Rs)) > log(EQF)) then
           EDC2 = 0d0 ; EDCD%PASSFAIL(13+n-1) = 0
       end if
       if ((EDC2 == 1 .or. DIAG == 1) .and. abs(Rs-Rm) > etol) then
           EDC2 = 0d0 ; EDCD%PASSFAIL(20+n-1) = 0
       end if
    end do

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
             if (M_POOLS(nn,n) < 0d0 .or. M_POOLS(nn,n) /= M_POOLS(nn,n)) then
                 EDC2 = 0d0 ; PEDC = 0 ; EDCD%PASSFAIL(35+n) = 0
             end if ! less than zero and is NaN condition
          nn = nn + 1
          end do ! nn < nodays .and. PEDC == 1
          n = n + 1
       end do ! for nopools .and. EDC .or. DIAG condition
    end if ! min pool assessment

  end subroutine assess_EDC2
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
    cal_mean_pools = 0d0

    ! loop through now
    cal_mean_pools = sum(pools(1:averaging_period,pool_number))/dble(averaging_period)

    ! ensure return command issued
    return

  end function cal_mean_pools
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
           ! uncertainty provided as +/-
           likelihood_p = likelihood_p-((pars(n)-parpriors(n))/parpriorunc(n))**2
           !likelihood_p = likelihood_p-0.5d0*((pars(n)-parpriors(n))/parpriorunc(n))**2
           ! uncertainty provided as fraction of observed value
           !likelihood_p=likelihood_p-0.5d0*((pars(n)-parpriors(n))/(parpriors(n)*parpriorunc(n)))**2
           ! uncertainty provided in log scale
           !likelihood_p=likelihood_p-0.5d0*(log(pars(n)/parpriors(n))/log(parpriorunc(n)))**2
       end if
    end do

    ! apply the 0.5 multiplicative which is part of the main likelihood calculation, here once.
    likelihood_p = likelihood_p * 0.5d0

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
    integer :: n, dn, no_years, y, s
    double precision :: tot_exp, tmp_var, infini
    double precision, allocatable :: mean_annual_pools(:)

    ! initial value
    likelihood = 0d0 ; infini = 0d0

    ! GPP Log-likelihood
    if (DATAin%ngpp > 0) then
       tot_exp = sum(((DATAin%M_GPP(DATAin%gpppts(1:DATAin%ngpp))-DATAin%GPP(DATAin%gpppts(1:DATAin%ngpp))) &
                       /DATAin%GPP_unc(DATAin%gpppts(1:DATAin%ngpp)))**2d0)
       likelihood = likelihood-tot_exp
    endif

    ! Fire Log-likelihood
    if (DATAin%nFire > 0) then
       tot_exp = sum(((DATAin%M_FLUXES(DATAin%Firepts(1:DATAin%nFire),17)-DATAin%Fire(DATAin%Firepts(1:DATAin%nFire))) &
                       /DATAin%Fire_unc(DATAin%Firepts(1:DATAin%nFire)))**2)
       likelihood = likelihood-tot_exp
    endif

    ! LAI log-likelihood
    if (DATAin%nlai > 0) then
       ! loop split to allow vectorisation
       tot_exp = sum(((DATAin%M_LAI(DATAin%laipts(1:DATAin%nlai))-DATAin%LAI(DATAin%laipts(1:DATAin%nlai))) &
                       /DATAin%LAI_unc(DATAin%laipts(1:DATAin%nlai)))**2d0)
       do n = 1, DATAin%nlai
         dn = DATAin%laipts(n)
         ! if zero or greater allow calculation with min condition to prevent
         ! errors of zero LAI which occur in managed systems
         if (DATAin%M_LAI(dn) < 0d0) then
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
                       /DATAin%NEE_unc(DATAin%neepts(1:DATAin%nnee)))**2d0)
       likelihood = likelihood-tot_exp
    endif

    ! Reco likelihood
    if (DATAin%nreco > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nreco
         dn = DATAin%recopts(n)
         tmp_var = DATAin%M_NEE(dn)+DATAin%M_GPP(dn)
         ! note that we calculate the Ecosystem resp from GPP and NEE
         tot_exp = tot_exp+((tmp_var-DATAin%Reco(dn))/DATAin%Reco_unc(dn))**2d0
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
       tot_exp = 0d0
       do n = 1, DATAin%nCfol_stock
         dn = DATAin%Cfol_stockpts(n)
         ! note that division is the uncertainty
!         tot_exp = tot_exp+(log(DATAin%M_POOLS(dn,1)/DATAin%Cfol_stock(dn))/log(2d0))**2d0
         tot_exp = tot_exp+((DATAin%M_POOLS(dn,1)-DATAin%Cfol_stock(dn)) / DATAin%Cfol_stock_unc(dn))**2d0
       end do
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
          mean_annual_pools(y) = cal_max_annual_pools(DATAin%M_POOLS(1:(DATAin%nodays+1),1),y,DATAin%deltat,DATAin%nodays+1)
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
         tot_exp = tot_exp+((tmp_var-DATAin%Cfolmax_stock(dn)) / DATAin%Cfolmax_stock_unc(dn))**2d0
       end do
       likelihood = likelihood-tot_exp
    endif

    ! Cwood log-likelihood (i.e. branch, stem and CR)
    if (DATAin%nCwood_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCwood_stock
         dn = DATAin%Cwood_stockpts(n)
         ! note that division is the uncertainty
!         tot_exp=tot_exp+(log(DATAin%M_POOLS(dn,3)/DATAin%Cwood_stock(dn))/log(2d0))**2d0
         tot_exp = tot_exp+((DATAin%M_POOLS(dn,3)-DATAin%Cwood_stock(dn))/DATAin%Cwood_stock_unc(dn))**2d0
       end do
       likelihood = likelihood-tot_exp
    endif

    ! Croots log-likelihood
    if (DATAin%nCroots_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCroots_stock
         dn = DATAin%Croots_stockpts(n)
         ! note that division is the uncertainty
!         tot_exp=tot_exp+(log(DATAin%M_POOLS(dn,2)/DATAin%Croots_stock(dn))/log(2d0))**2d0
         tot_exp = tot_exp+((DATAin%M_POOLS(dn,2)-DATAin%Croots_stock(dn)) / DATAin%Croots_stock_unc(dn))**2d0
       end do
       likelihood = likelihood-tot_exp
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
                           *(DATAin%M_POOLS(dn,4))-DATAin%Clit_stock(dn))/DATAin%Clit_stock_unc(dn))**2d0
       end do
       likelihood = likelihood-tot_exp
    endif

    ! Csom log-likelihood
    if (DATAin%nCsom_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCsom_stock
         dn = DATAin%Csom_stockpts(n)
         ! note that division is the uncertainty
!         tot_exp=tot_exp+(log(DATAin%M_POOLS(dn,5)/DATAin%Csom_stock(dn))/log(2d0))**2d0
         tot_exp = tot_exp+((DATAin%M_POOLS(dn,5)-DATAin%Csom_stock(dn))/DATAin%Csom_stock_unc(dn))**2d0
       end do
       likelihood = likelihood-tot_exp
    endif

    !
    ! Curiously we will assess other priors here, as the tend to have to do with model state derived values
    !

    ! Ra:GPP fraction is in this model a derived property
    if (DATAin%otherpriors(1) > 0) then
        tot_exp = sum(DATAin%M_FLUXES(:,3)) / sum(DATAin%M_FLUXES(:,1))
        likelihood = likelihood-((tot_exp-DATAin%otherpriors(1))/DATAin%otherpriorunc(1))**2d0
    end if

    ! Estimate the biological steady state attractor on the wood pool.
    ! NOTE: this arrangement explicitly neglects the impact of disturbance on
    ! residence time (i.e. no fire and biomass removal)
    if (DATAin%otherpriors(5) > -9998) then
        ! Estimate the mean annual input to the wood pool (gC.m-2.day-1) and
        ! remove the day-1 by multiplying by residence time (day)
        tot_exp = (sum(DATAin%M_FLUXES(:,7)) / dble(DATAin%nodays)) * (pars(6) ** (-1d0))
        likelihood = likelihood - ((tot_exp - DATAin%otherpriors(5)) / DATAin%otherpriorunc(5))**2
    endif

    ! the likelihood scores for each observation are subject to multiplication
    ! by 0.5 in the algebraic formulation. To avoid repeated calculation across
    ! multiple datastreams we apply this multiplication to the bulk liklihood
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

    ! calculates the likelihood of of the model output compared to the available
    ! observations which have been input to the model

    implicit none

    ! declare arguments
    integer, intent(in) :: npars
    double precision, dimension(npars), intent(in) :: pars

    ! declare local variables
    integer :: n, dn, no_years, y, s
    double precision :: tot_exp, tmp_var, infini
    double precision, allocatable :: mean_annual_pools(:)

    ! initial value
    scale_likelihood = 0d0 ; infini = 0d0

    ! GPP Log-likelihood
    if (DATAin%ngpp > 0) then
       tot_exp = sum(((DATAin%M_GPP(DATAin%gpppts(1:DATAin%ngpp))-DATAin%GPP(DATAin%gpppts(1:DATAin%ngpp))) &
                       /DATAin%GPP_unc(DATAin%gpppts(1:DATAin%ngpp)))**2d0)
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%ngpp))
    endif

    ! Fire Log-likelihood
    if (DATAin%nFire > 0) then
       tot_exp = sum(((DATAin%M_FLUXES(DATAin%Firepts(1:DATAin%nFire),17)-DATAin%Fire(DATAin%Firepts(1:DATAin%nFire))) &
                       /DATAin%Fire_unc(DATAin%Firepts(1:DATAin%nFire)))**2)
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nFire))
    endif

    ! LAI log-likelihood
    if (DATAin%nlai > 0) then
       ! loop split to allow vectorisation
       tot_exp = sum(((DATAin%M_LAI(DATAin%laipts(1:DATAin%nlai))-DATAin%LAI(DATAin%laipts(1:DATAin%nlai))) &
                       /DATAin%LAI_unc(DATAin%laipts(1:DATAin%nlai)))**2d0)
       do n = 1, DATAin%nlai
         dn = DATAin%laipts(n)
         ! if zero or greater allow calculation with min condition to prevent
         ! errors of zero LAI which occur in managed systems
         if (DATAin%M_LAI(dn) < 0d0) then
             ! if not then we have unrealistic negative values or NaN so indue
             ! error
             tot_exp = tot_exp+(-log(infini))
         endif
       end do
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nlai))
    endif

    ! NEE likelihood
    if (DATAin%nnee > 0) then
       tot_exp = sum(((DATAin%M_NEE(DATAin%neepts(1:DATAin%nnee))-DATAin%NEE(DATAin%neepts(1:DATAin%nnee))) &
                       /DATAin%NEE_unc(DATAin%neepts(1:DATAin%nnee)))**2d0)
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nnee))
    endif

    ! Reco likelihood
    if (DATAin%nreco > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nreco
         dn = DATAin%recopts(n)
         tmp_var = DATAin%M_NEE(dn)+DATAin%M_GPP(dn)
         ! note that we calculate the Ecosystem resp from GPP and NEE
         tot_exp = tot_exp+((tmp_var-DATAin%Reco(dn))/DATAin%Reco_unc(dn))**2d0
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
!         tot_exp = tot_exp+(log(DATAin%M_POOLS(dn,1)/DATAin%Cfol_stock(dn))/log(2d0))**2d0
         tot_exp = tot_exp+((DATAin%M_POOLS(dn,1)-DATAin%Cfol_stock(dn)) / DATAin%Cfol_stock_unc(dn))**2d0
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
          mean_annual_pools(y) = cal_max_annual_pools(DATAin%M_POOLS(1:(DATAin%nodays+1),1),y,DATAin%deltat,DATAin%nodays+1)
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
         tot_exp = tot_exp+((tmp_var-DATAin%Cfolmax_stock(dn)) / DATAin%Cfolmax_stock_unc(dn))**2d0
       end do
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nCfolmax_stock))
    endif

    ! Cwood log-likelihood (i.e. branch, stem and CR)
    if (DATAin%nCwood_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCwood_stock
         dn = DATAin%Cwood_stockpts(n)
         ! note that division is the uncertainty
!         tot_exp=tot_exp+(log(DATAin%M_POOLS(dn,3)/DATAin%Cwood_stock(dn))/log(2d0))**2d0
         tot_exp = tot_exp+((DATAin%M_POOLS(dn,3)-DATAin%Cwood_stock(dn))/DATAin%Cwood_stock_unc(dn))**2d0
       end do
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nCwood_stock))
    endif

    ! Croots log-likelihood
    if (DATAin%nCroots_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCroots_stock
         dn = DATAin%Croots_stockpts(n)
         ! note that division is the uncertainty
!         tot_exp=tot_exp+(log(DATAin%M_POOLS(dn,2)/DATAin%Croots_stock(dn))/log(2d0))**2d0
         tot_exp = tot_exp+((DATAin%M_POOLS(dn,2)-DATAin%Croots_stock(dn)) / DATAin%Croots_stock_unc(dn))**2d0
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
                           *(DATAin%M_POOLS(dn,4))-DATAin%Clit_stock(dn))/DATAin%Clit_stock_unc(dn))**2d0
       end do
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nClit_stock))
    endif

    ! Csom log-likelihood
    if (DATAin%nCsom_stock > 0) then
       tot_exp = 0d0
       do n = 1, DATAin%nCsom_stock
         dn = DATAin%Csom_stockpts(n)
         ! note that division is the uncertainty
!         tot_exp=tot_exp+(log(DATAin%M_POOLS(dn,5)/DATAin%Csom_stock(dn))/log(2d0))**2d0
         tot_exp = tot_exp+((DATAin%M_POOLS(dn,5)-DATAin%Csom_stock(dn))/DATAin%Csom_stock_unc(dn))**2d0
       end do
       scale_likelihood = scale_likelihood-(tot_exp/dble(DATAin%nCsom_stock))
    endif

    !
    ! Curiously we will assess other priors here, as the tend to have to do with model state derived values
    !

    ! Ra:GPP fraction is in this model a derived property
    if (DATAin%otherpriors(1) > 0) then
        tot_exp = sum(DATAin%M_FLUXES(:,3)) / sum(DATAin%M_FLUXES(:,1))
        scale_likelihood = scale_likelihood-((tot_exp-DATAin%otherpriors(1))/DATAin%otherpriorunc(1))**2d0
    end if

    ! Estimate the biological steady state attractor on the wood pool.
    ! NOTE: this arrangement explicitly neglects the impact of disturbance on
    ! residence time (i.e. no fire and biomass removal)
    if (DATAin%otherpriors(5) > -9998) then
        ! Estimate the mean annual input to the wood pool (gC.m-2.day-1) and
        ! remove the day-1 by multiplying by residence time (day)
        tot_exp = (sum(DATAin%M_FLUXES(:,7)) / dble(DATAin%nodays)) * (pars(6) ** (-1d0))
        scale_likelihood = scale_likelihood - ((tot_exp - DATAin%otherpriors(5)) / DATAin%otherpriorunc(5))**2
    endif

    ! the likelihood scores for each observation are subject to multiplication
    ! by 0.5 in the algebraic formulation. To avoid repeated calculation across
    ! multiple datastreams we apply this multiplication to the bulk liklihood
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
