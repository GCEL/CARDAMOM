!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CARbon DAta MOdel fraMework (CARDAMOM) and DALEC terrestrial ecosystem model suite
! CARDAMOM is a Bayesian model-data fusion software framework. CARDAMOM is used to 
! assimilate observations and ecological theory to retrieve parameters for the 
! DALEC suite of intermediate complexity terrestrial ecosystem models. DALEC can be
! used as a fully integrated component of CARDAMOM or independently. 
! Copyright (C) 2024  University of Edinburgh,
!                     Mathew Williams (mat.williams@ed.ac.uk), 
!                     T. Luke Smallman (t.l.smallman@ed.ac.uk)
! UoE = University of Edinburgh

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

!!!!!!!!!!!! File specific description !!!!!!!!!!
! Module contains all subroutine and functions relevant to determining the log-likelihood
! of DALEC.C3.M1.014 as a function of observations and ecological dynamical constraints.
!
! This code is based on the original C verion of the University of Edinburgh
! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
! All code translation into Fortran, integration into the University of
! Edinburgh CARDAMOM code and subsequent modifications by:
! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
! See function / subroutine specific comments for exceptions and contributors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module model_likelihood_module
  implicit none

  ! make all private
  private

  ! which to make open
  public :: model_likelihood, find_edc_initial_values, &
            sub_model_likelihood, sqrt_model_likelihood, log_model_likelihood

  ! declare needed types
  type EDCDIAGNOSTICS
    integer :: nedc = 150    ! number of edcs being assessed
    integer :: PASSFAIL(150) ! allow space for 150 possible checks, dim should equal nedc
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
                     ,DATAin%nodays,DATAin%LAT &
                     ,M_FLUXES,M_POOLS,M_DIAGS &
                     ,DATAin%nopars,DATAin%nomet,DATAin%nopools     &
                     ,DATAin%nofluxes,DATAin%nodiags                &
                     ,PI%stock_seed_labile,PI%DS_shoot,PI%DS_root   &
                     ,PI%fol_frac,PI%stem_frac,PI%root_frac,PI%DS_LRLV&
                     ,PI%LRLV,PI%DS_LRRT,PI%LRRT)

    ! assess post running EDCs
    call assess_EDC2(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools  &
                    ,DATAin%nodays,DATAin%nodiags,DATAin%deltat            &
                    ,DATAin%steps_per_year,PI%parmax,PARS,DATAin%MET       &
                    ,M_POOLS,M_FLUXES,M_DIAGS         &
                    ,DATAin%meantemp,EDC2)

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
  subroutine model_sanity_check(PARS)
    use cardamom_structures, only: DATAin
    use MCMCOPT, only: PI
    use carbon_model_mod, only: carbon_model

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
    double precision, dimension(DATAin%nodays,DATAin%nodiags) :: local_diags
    double precision :: pool_error, flux_error, diags_error

    ! Run model

    ! next need to run the model itself
    call carbon_model(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT &
                     ,M_FLUXES,M_POOLS,M_DIAGS &
                     ,DATAin%nopars,DATAin%nomet,DATAin%nopools     &
                     ,DATAin%nofluxes,DATAin%nodiags                &
                     ,PI%stock_seed_labile,PI%DS_shoot,PI%DS_root   &
                     ,PI%fol_frac,PI%stem_frac,PI%root_frac,PI%DS_LRLV&
                     ,PI%LRLV,PI%DS_LRRT,PI%LRRT)

    call carbon_model(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT &
                     ,local_fluxes,local_pools,local_diags &
                     ,DATAin%nopars,DATAin%nomet,DATAin%nopools     &
                     ,DATAin%nofluxes,DATAin%nodiags                &
                     ,PI%stock_seed_labile,PI%DS_shoot,PI%DS_root   &
                     ,PI%fol_frac,PI%stem_frac,PI%root_frac,PI%DS_LRLV&
                     ,PI%LRLV,PI%DS_LRRT,PI%LRRT)

    ! Compare outputs
    flux_error = sum(abs(M_FLUXES - local_fluxes))
    pool_error = sum(abs(M_POOLS - local_pools))
    diags_error = sum(abs(M_DIAGS - local_diags))
    ! If error between runs exceeds precision error then we have a problem
    if (diags_error > (tiny(0d0)*(DATAin%nopools*DATAin%nodays)) .or. &
        pool_error  > (tiny(0d0)*(DATAin%nopools*DATAin%nodays)) .or. &
        flux_error  > (tiny(0d0)*(DATAin%nofluxes*DATAin%nodays)) .or. &
        pool_error /= pool_error .or. flux_error /= flux_error .or. &
        diags_error /= diags_error) then
        print*,"Error: multiple runs of the same parameter set indicates an error"
        print*,"Cumulative POOL error = ",pool_error
        print*,"Cumulative FLUX error = ",flux_error
        do i = 1,DATAin%nofluxes
           print*,"Sum abs error over time: flux = ",i
           print*,sum(abs(M_FLUXES(:,i) - local_fluxes(:,i)))
        end do
        do i = 1, DATAin%nopools
           print*,"Sum abs error over time: pool = ",i
           print*,sum(abs(M_POOLS(:,i) - local_pools(:,i)))
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

    ! pre_DR should be greater than post_DR, this is consistent across currently
    ! available SPA crop parameter files
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(4) > pars(3))) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(3) = 0
    endif

    ! for development: Tmin should be < topt and topt should be < tmax
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(26) > pars(28) &
                                     .or. pars(28) > pars(27) &
                                     .or. pars(26) > pars(27))) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(4) = 0
    endif

    ! for development: the difference between each Tmin,Topt,Tmax > 1.
    if ((EDC1 == 1 .or. DIAG == 1) .and. (abs(pars(26)-pars(28)) < 1d0 &
                                     .or. abs(pars(28)-pars(27)) < 1d0  &
                                     .or. abs(pars(26)-pars(27)) < 1d0)) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(5) = 0
    endif

    ! for vernalisation: Tmin < Topt < Tmax
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(29) > pars(31) &
                                     .or. pars(31) > pars(30) &
                                     .or. pars(29) > pars(30))) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(6) = 0
    endif

   ! for vernalisation: the difference between each Tmin, Topt, Tmax
    if ((EDC1 == 1 .or. DIAG == 1) .and. ( abs(pars(29)-pars(31)) < 1d0 &
                                      .or. abs(pars(31)-pars(30)) < 1d0 &
                                      .or. abs(pars(29)-pars(30)) < 1d0 ) ) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(7) = 0
    endif

   ! development temperature value should be larger corresponding vernalisation
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(29) > pars(26) &
                                     .or. pars(31) > pars(28) &
                                     .or. pars(30) > pars(27))) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(8) = 0
    endif

    ! could and probably should add some more
    
  end subroutine assess_EDC1
  !
  !------------------------------------------------------------------
  !
  subroutine assess_EDC2(npars,nomet,nofluxes,nopools,nodays,nodiags &
                        ,deltat,steps_per_year,parmax,pars,met &
                        ,M_POOLS,M_FLUXES,M_DIAGS &
                        ,meantemp,EDC2)

    use cardamom_structures, only: DATAin
    use CARBON_MODEL_MOD, only: resp_rate_temp_coeff,linear_model_gradient

    ! the second of two subroutines for assessing current parameters for passing
    ! realism tests for crop ecosystems

    implicit none

    ! declare input variables
    integer, intent(in) :: npars    & ! number of model parameters
                          ,nomet    & ! number of met drivers
                          ,nofluxes & ! number of fluxes from model
                          ,nopools  & ! number of pools in model
                          ,nodays   & ! number of days in simulation
                          ,nodiags  & ! number of diagnostics in model
                          ,steps_per_year

    double precision, intent(in) :: deltat(nodays)              & ! decimal day model interval
                                   ,pars(npars)                 & ! vector of current parameters
                                   ,parmax(npars)               & ! vector of the maximum parameter values
                                   ,met(nomet,nodays)           & ! array of met drivers
                                   ,M_POOLS((nodays+1),nopools) & ! time varying states of pools in current model simulation
                                   ,M_FLUXES(nodays,nofluxes)   & ! time varying fluxes from current model simulation
                                   ,M_DIAGS(nodays,nodiags)     & ! time varying diagnostics from current model simulation
                                   ,meantemp                      ! site mean temperature (oC)

    double precision, intent(out) :: EDC2 ! the response flag for the dynamical set of EDCs

    ! declare local variables
    integer :: n, nn, nnn, DIAG, no_years, y, PEDC, &
               nd, fl, fs, io_start, io_finish
    integer, dimension(nodays) :: pool_hak
    double precision :: infi, mrt!, EQF, etol
    double precision, dimension(nodays) :: ratio
    double precision, dimension(nopools) :: jan_mean_pools, jan_first_pools, &
                                            mean_pools, Fin, Fout, Rm, Rs, &
                                            Fin_yr1, Fout_yr1, Fin_yr2, Fout_yr2
    double precision, dimension(nofluxes) :: FT, FT_yr1, FT_yr2
    
    ! Steady State Attractor:
    ! Log ratio difference between inputs and outputs of the system.
    logical, parameter :: old_edcs = .false.
    double precision, parameter :: EQF1_5 = log(1.5d0), & ! 10.0 = order magnitude; 2 = double and half
                                   EQF2 = log(2d0),   & ! 10.0 = order magnitude; 2 = double and half
                                   EQF5 = log(5d0),   &
                                   EQF10 = log(10d0), &
                                   EQF15 = log(15d0), &
                                   EQF20 = log(20d0), &
                                  C_etol = 0.20d0       ! 0.20d0 lots of AGB !0.10d0 global / site more data !0.05d0 global 1 or 2 AGB estimates

    ! Work out how many completed years there are in the system
    no_years = int(nint(sum(deltat)/365.25d0))

    ! initial value
    infi = 0d0
    ! update initial values
    DIAG = EDCD%DIAG
    ! give EDC2 an initial value
    EDC2 = 1

    ! First calculate total flux for the simulation period
    ! NOTE: that this code differs from the majority of DALEC models
    ! as we treat the first year as a spin up to be ignored.
    io_start = steps_per_year + 1 ; io_finish = nodays
    if (no_years < 3) then 
       do fl = 1, nofluxes
          FT(fl) = sum(M_FLUXES(io_start:io_finish,fl)*deltat(io_start:io_finish))
          FT_yr1(fl) = sum(M_FLUXES((steps_per_year+1):(steps_per_year*2),fl)*deltat((steps_per_year+1):(steps_per_year*2)))
          FT_yr2(fl) = FT_yr1(fl)
       end do                       
    else 
       do fl = 1, nofluxes
          FT(fl) = sum(M_FLUXES(io_start:io_finish,fl)*deltat(io_start:io_finish))
          FT_yr1(fl) = sum(M_FLUXES((steps_per_year+1):(steps_per_year*2),fl)*deltat((steps_per_year+1):(steps_per_year*2)))
          FT_yr2(fl) = sum(M_FLUXES(((steps_per_year*2)+1):(steps_per_year*3),fl)*deltat(((steps_per_year*2)+1):(steps_per_year*3)))
       end do
    end if 

    ! Get total in and out for each dead organic matter pool

    ! litter
    Fin(5)  = FT(12)+FT(32)+FT(33)+FT(34)+FT(35)+FT(36)+FT(37)
    Fout(5) = FT(13)+FT(15)
    Fin_yr1(5)  = FT_yr1(12)+FT_yr1(32)+FT_yr1(33)+FT_yr1(34)+FT_yr1(35)+FT_yr1(36)+FT_yr1(37)
    Fout_yr1(5) = FT_yr1(13)+FT_yr1(15)
    Fin_yr2(5)  = FT_yr2(12)+FT_yr2(32)+FT_yr2(33)+FT_yr2(34)+FT_yr2(35)+FT_yr2(36)+FT_yr2(37)
    Fout_yr2(5) = FT_yr2(13)+FT_yr2(15)
    ! som
    Fin(6)  = FT(15)
    Fout(6) = FT(14)
    Fin_yr1(6)  = FT_yr1(15)
    Fout_yr1(6) = FT_yr1(14)
    Fin_yr2(6)  = FT_yr2(15)
    Fout_yr2(6) = FT_yr2(14)

    if (EDC2 == 1 .or. DIAG == 1) then

        ! Foliage + fine root litter
        ! Estimate MRT (years)
        pool_hak = 1 ; ratio = 0d0
        where (M_POOLS(1:nodays,5) > 0d0) ! protection against NaN from division by zero
               pool_hak = 0 
               ratio = ((M_FLUXES(1:nodays,13) + M_FLUXES(1:nodays,15)) &
                       / M_POOLS(1:nodays,5))
        end where
        ! Estimate the mean fractional daily loss
        mrt = sum(ratio) / dble(nodays-sum(pool_hak))
        ! If daily turnover fraction is less than equivalent of MRT of 1 year, fail.
        if (mrt < 2.737850787d-3) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(10) = 0
        end if 

    end if ! EDC2 == 1 .or. DIAG == 1

    if (EDC2 == 1 .or. DIAG == 1) then

        ! Dead pools - SOM only
        do n = 5, 6
           ! Restrict rates of increase
           if (abs(log(Fin(n)/Fout(n))) > EQF10) then
               EDC2 = 0d0 ; EDCD%PASSFAIL(10+n-4) = 0
           end if
           ! Restrict exponential behaviour at initialisation
           if (abs(abs(log(Fin_yr1(n)/Fout_yr1(n))) - abs(log(Fin_yr2(n)/Fout_yr2(n)))) > C_etol) then
               EDC2 = 0d0 ; EDCD%PASSFAIL(16+n-4) = 0
           end if
        end do

    end if ! EDC2 == 1 .or. DIAG == 1

    ! we know that the crop model should produce some yield - therefore we
    ! reject parameter sets which generate no yield ever!
    !if ((EDC2 == 1 .or. DIAG == 1) .and. sum(M_FLUXES(1:nodays,21)) < (1d0*dble(no_years)) ) then
    !    EDC2 = 0d0 ; EDCD%PASSFAIL(17) = 0
    !endif
    !! Total hack to enforce a massive yield and find out what the parameters do
    !if ((EDC2 == 1 .or. DIAG == 1) .and. maxval(M_FLUXES(1:nodays,21)) < 300d0 ) then
    !  EDC2 = 0d0 ; EDCD%PASSFAIL(18) = 0
    !endif

    ! We should assume all crops get somewhere close to maturity (2.0)
    if ((EDC2 == 1 .or. DIAG == 1) .and. maxval(M_DIAGS(1:nodays,3)) < 1.9) then
        EDC2 = 0d0 ; EDCD%PASSFAIL(19) = 0
    endif

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
    cal_mean_annual_pools = sum(pools(startday:endday))/dble(endday-startday+1)

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
    use carbon_model_mod, only: carbon_model
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
    call carbon_model(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT                      &
                     ,M_FLUXES,M_POOLS,M_DIAGS &
                     ,DATAin%nopars,DATAin%nomet,DATAin%nopools     &
                     ,DATAin%nofluxes,DATAin%nodiags                &
                     ,PI%stock_seed_labile,PI%DS_shoot,PI%DS_root   &
                     ,PI%fol_frac,PI%stem_frac,PI%root_frac,PI%DS_LRLV&
                     ,PI%LRLV,PI%DS_LRRT,PI%LRRT)

    ! if first set of EDCs have been passed, move on to the second
    if (DATAin%EDC == 1) then
        ! check edc2
        call assess_EDC2(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools  &
                        ,DATAin%nodays,DATAin%nodiags,DATAin%deltat            &
                        ,DATAin%steps_per_year,PI%parmax,PARS,DATAin%MET       &
                        ,M_POOLS,M_FLUXES,M_DIAGS         &
                        ,DATAin%meantemp,EDC2)                        

        ! Add EDC2 log-likelihood to absolute accept reject...
        ML_obs_out = ML_obs_out + log(EDC2)
    end if ! DATAin%EDC == 1

    ! Calculate log-likelihood associated with priors
    ! We always want this
    ML_prior_out = likelihood_p(PI%npars,DATAin%parpriors,DATAin%parpriorunc,DATAin%parpriorweight,PARS)
    ! Calculate log-likelihood of model compared to obs
    call calc_obs_likelihoods(ML_obs_out)
    ! Calculate log-likelihood of 'other priors'
    call calc_other_likelihoods(ML_obs_out)

  end subroutine model_likelihood
  !
  !------------------------------------------------------------------
  !  
  subroutine scaled_model_likelihood(PARS,ML_obs_out,ML_prior_out)
    use MCMCOPT, only:  PI
    use carbon_model_mod, only: carbon_model
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
                     ,DATAin%nodays,DATAin%LAT                      &
                     ,M_FLUXES,M_POOLS,M_DIAGS &
                     ,DATAin%nopars,DATAin%nomet,DATAin%nopools     &
                     ,DATAin%nofluxes,DATAin%nodiags                &
                     ,PI%stock_seed_labile,PI%DS_shoot,PI%DS_root   &
                     ,PI%fol_frac,PI%stem_frac,PI%root_frac,PI%DS_LRLV&
                     ,PI%LRLV,PI%DS_LRRT,PI%LRRT)

    ! if first set of EDCs have been passed, move on to the second
    if (DATAin%EDC == 1) then

        ! check edc2
        call assess_EDC2(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools  &
                        ,DATAin%nodays,DATAin%nodiags,DATAin%deltat            &
                        ,DATAin%steps_per_year,PI%parmax,PARS,DATAin%MET       &
                        ,M_POOLS,M_FLUXES,M_DIAGS         &
                        ,DATAin%meantemp,EDC2)      
        ! Add EDC2 log-likelihood to absolute accept reject...
        ML_obs_out = ML_obs_out + log(EDC2)

    end if ! DATAin%EDC == 1

    ! Calculate log-likelihood associated with priors
    ! We always want this
    ML_prior_out = likelihood_p(PI%npars,DATAin%parpriors,DATAin%parpriorunc,DATAin%parpriorweight,PARS)
    ! Calculate log-likelihood of model compared to obs
    call calc_scaled_obs_likelihoods(ML_obs_out)
    ! Calculate log-likelihood of 'other priors'
    call calc_other_likelihoods(ML_obs_out)

  end subroutine scaled_model_likelihood
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
  subroutine calc_obs_likelihoods(ML_obs_out)
    use cardamom_structures, only: DATAin

    ! Subroutine to control the calculation of the observation related
    ! log-likelihoods and accounting for the various scalings

    ! Arguements
    double precision, intent(inout) :: ML_obs_out
    ! local variable
    double precision, dimension(DATAin%nodays) :: mod

    !
    ! Do diagnostics (DIAGS)
    ! 

    ! Calculate log-likelihood for leaf area index
    if (DATAin%nlai > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nlai,DATAin%laipts,DATAin%LAI,DATAin%LAI_unc,DATAin%LAI_lag, &
                                             1d0,M_DIAGS(1:DATAin%nodays,1))
    end if

    !
    ! Do pools (POOLS)
    !

    ! Calculate log-likelihood for foliage stocks
    if (DATAin%nCfol_stock > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCfol_stock,DATAin%Cfol_stockpts, &
                                             DATAin%Cfol_stock,DATAin%Cfol_stock_unc,DATAin%Cfol_stock_lag, &
                                             1d0,M_POOLS(1:DATAin%nodays,2))
    endif ! nCfol_stock > 0
    ! Calculate log-likelihood for fine root stocks
    if (DATAin%nCroots_stock > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCroots_stock,DATAin%Croots_stockpts, &
                                             DATAin%Croots_stock,DATAin%Croots_stock_unc,DATAin%Croots_stock_lag, &
                                             1d0,M_POOLS(1:DATAin%nodays,3))
    endif ! nCroots_stock > 0
    ! Calculate log-likelihood for soil organic matter stocks
    if (DATAin%nCsom_stock > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCsom_stock,DATAin%Csom_stockpts, &
                                             DATAin%Csom_stock,DATAin%Csom_stock_unc,DATAin%Csom_stock_lag, &
                                             1d0,M_POOLS(1:DATAin%nodays,6))
    endif ! nCsom_stock > 0

    !
    ! Do fluxes (FLUXES)
    !

    ! Calculate log-likelihood for gross primary production
    if (DATAin%ngpp > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%ngpp,DATAin%gpppts,DATAin%GPP,DATAin%GPP_unc,DATAin%GPP_lag, &
                                             1d0,M_FLUXES(1:DATAin%nodays,1))
    endif ! ngpp > 0
    ! Calculate log-likelihood for harvest
    if (DATAin%nharvest > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nharvest,DATAin%harvestpts,DATAin%harvest,DATAin%harvest_unc,DATAin%harvest_lag, &
                                             1d0,M_FLUXES(1:DATAin%nodays,21))
    endif ! nharvest > 0
    ! Calculate log-likelihood for net ecosystem exchange of CO2
    if (DATAin%nnee > 0) then
        mod = M_FLUXES(1:DATAin%nodays,3) &  ! Rauto
            + M_FLUXES(1:DATAin%nodays,13) & ! Rhet litter
            + M_FLUXES(1:DATAin%nodays,14) & ! Rhet som
            - M_FLUXES(1:DATAin%nodays,1)    ! GPP
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nnee,DATAin%neepts,DATAin%NEE,DATAin%NEE_unc,DATAin%NEE_lag, &
                                             1d0,mod)    
    endif ! nnee > 0
    ! Calculate log-likelihood for ecosystem respiration
    if (DATAin%nreco > 0) then
        mod = M_FLUXES(1:DATAin%nodays,3) &  ! Rauto
            + M_FLUXES(1:DATAin%nodays,13) & ! Rhet litter
            + M_FLUXES(1:DATAin%nodays,14)   ! Rhet som
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nreco,DATAin%recopts,DATAin%Reco,DATAin%Reco_unc,DATAin%Reco_lag, &
                                             1d0,mod)
    endif ! nreco > 0

    return

  end subroutine calc_obs_likelihoods  
  !
  !------------------------------------------------------------------
  !
  subroutine calc_scaled_obs_likelihoods(ML_obs_out)
    use cardamom_structures, only: DATAin

    ! Subroutine to control the calculation of the observation related
    ! log-likelihoods and accounting for the various scalings

    ! Arguements
    double precision, intent(inout) :: ML_obs_out
    ! local variable
    double precision, dimension(DATAin%nodays) :: mod

    !
    ! Do diagnostics (DIAGS)
    ! 

    ! Calculate log-likelihood for leaf area index
    if (DATAin%nlai > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nlai,DATAin%laipts,DATAin%LAI,DATAin%LAI_unc,DATAin%LAI_lag, &
                                             DATAin%LAI_scaling,M_DIAGS(1:DATAin%nodays,1))
    end if ! nlai > 0

    !
    ! Do pools (POOLS)
    !

    ! Calculate log-likelihood for foliage stocks
    if (DATAin%nCfol_stock > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCfol_stock,DATAin%Cfol_stockpts, &
                                             DATAin%Cfol_stock,DATAin%Cfol_stock_unc,DATAin%Cfol_stock_lag, &
                                             DATAin%Cfol_stock_scaling,M_POOLS(1:DATAin%nodays,2))
    endif ! nCfol_stock > 0
    ! Calculate log-likelihood for fine root stocks
    if (DATAin%nCroots_stock > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCroots_stock,DATAin%Croots_stockpts, &
                                             DATAin%Croots_stock,DATAin%Croots_stock_unc,DATAin%Croots_stock_lag, &
                                             DATAin%Croots_stock_scaling,M_POOLS(1:DATAin%nodays,3))
    endif ! nCroots_stock > 0
    ! Calculate log-likelihood for soil organic matter stocks
    if (DATAin%nCsom_stock > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCsom_stock,DATAin%Csom_stockpts, &
                                             DATAin%Csom_stock,DATAin%Csom_stock_unc,DATAin%Csom_stock_lag, &
                                             DATAin%Csom_stock_scaling,M_POOLS(1:DATAin%nodays,6))
    endif ! nCsom_stock > 0

    !
    ! Do fluxes (FLUXES)
    !

    ! Calculate log-likelihood for gross primary production
    if (DATAin%ngpp > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%ngpp,DATAin%gpppts,DATAin%GPP,DATAin%GPP_unc,DATAin%GPP_lag, &
                                             DATAin%GPP_scaling,M_FLUXES(1:DATAin%nodays,1))
    endif ! ngpp > 0
    ! Calculate log-likelihood for harvest
    if (DATAin%nharvest > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nharvest,DATAin%harvestpts,DATAin%harvest,DATAin%harvest_unc,DATAin%harvest_lag, &
                                             DATAin%harvest_scaling,M_FLUXES(1:DATAin%nodays,21))
    endif ! nharvest > 0
    ! Calculate log-likelihood for net ecosystem exchange of CO2
    if (DATAin%nnee > 0) then
        mod = M_FLUXES(1:DATAin%nodays,3) &  ! Rauto
            + M_FLUXES(1:DATAin%nodays,13) & ! Rhet litter
            + M_FLUXES(1:DATAin%nodays,14) & ! Rhet som
            - M_FLUXES(1:DATAin%nodays,1)    ! GPP
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nnee,DATAin%neepts,DATAin%NEE,DATAin%NEE_unc,DATAin%NEE_lag, &
                                             DATAin%NEE_scaling,mod)    
    endif ! nnee > 0
    ! Calculate log-likelihood for ecosystem respiration
    if (DATAin%nreco > 0) then
        mod = M_FLUXES(1:DATAin%nodays,3)  & ! Rauto
            + M_FLUXES(1:DATAin%nodays,13) & ! Rhet litter
            + M_FLUXES(1:DATAin%nodays,14)   ! Rhet som
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nreco,DATAin%recopts,DATAin%Reco,DATAin%Reco_unc,DATAin%Reco_lag, &
                                             DATAin%Reco_scaling,mod)
    endif ! nreco > 0

    return

  end subroutine calc_scaled_obs_likelihoods  
  !
  !------------------------------------------------------------------
  !
  subroutine calc_other_likelihoods(ML_obs_out)
    use cardamom_structures, only: DATAin

    ! Subroutine to control the calculation of the 'other priors'
    ! log-likelihoods. These are typically derived variables / 
    ! emergent functions averaged over the life of the anaylsis

    ! Arguements
    double precision, intent(inout) :: ML_obs_out
    ! local variable
    integer :: dummy_nodays = 1, dummy_noobs = 1
    integer, dimension(1) :: dummy_pts = 1, dummy_lag = 0
    double precision, dimension(1) :: mod
    double precision :: dummy_scaling = 1d0
    double precision, allocatable, dimension(:) :: tmp1

    ! Ra:GPP fraction is in this model a derived property
    if (DATAin%otherpriors(1) > -9998) then
        ! Estimate the foliage litter pool based on the ratio of foliage litter input to foliage + fine root litter inputs,
        ! scaled by the total litter pool. This is based on the turnover being common.
        mod = sum(M_FLUXES(1:DATAin%nodays,3)) / sum(M_FLUXES(1:DATAin%nodays,1)) ! sum(Rauto) / sum(GPP)
        ML_obs_out = ML_obs_out + (DATAin%otherpriorweight(1)*likelihood(dummy_nodays,dummy_noobs,dummy_pts, &
                                   DATAin%otherpriors(1),DATAin%otherpriorunc(1),dummy_lag,dummy_scaling,mod))
    end if

    ! ...OTHERPRIOR(2-7)...

    ! Yield:GPP fraction is in this model a derived property
    if (DATAin%otherpriors(8) > -9998) then
        ! Accumulate yield and GPP over the growing period, based on DS >= 0.
        ! This code assumes that the DS_time = DS occurs after DS is incremented and 
        ! not after the management activities had reset DS to -1. If so this code will not work.
        if (allocated(tmp1)) deallocate(tmp1)
        allocate(tmp1(DATAin%nodays)) ; tmp1 = 0d0 ; where(M_DIAGS(1:DATAin%nodays,3) >= 0d0) tmp1 = 1d0
        mod = sum(M_FLUXES(:,21)*tmp1) / sum(M_FLUXES(:,1)*tmp1)
        ML_obs_out = ML_obs_out + (DATAin%otherpriorweight(8)*likelihood(dummy_nodays,dummy_noobs,dummy_pts, &
                                   DATAin%otherpriors(8),DATAin%otherpriorunc(8),dummy_lag,dummy_scaling,mod))
        deallocate(tmp1)
    end if

    return

  end subroutine calc_other_likelihoods  
  !
  !------------------------------------------------------------------
  !
  double precision function likelihood(nodays,nobs,obspts,obs,unc,lag,scaling,mod)

    ! Generic function to estimate the likelihood of diagnostic.

    ! Arguments
    integer, intent(in) :: nodays, & ! number of time steps
                             nobs    ! number of observed estimates
    integer, dimension(nobs), intent(in) :: obspts ! location of observations in time series
    integer, dimension(nodays), intent(in) :: lag    ! observation lag period
    double precision, intent(in) :: scaling ! scaling factor to account for differing number of observations
    double precision, dimension(nodays), intent(in) :: obs, & ! observation time series
                                                       unc, & ! observation uncertainty
                                                       mod    ! model equivalent of the obs

    ! local variables
    integer :: dn, n, s
    double precision :: infini
    
    ! Set initial value
    infini = 0d0

    ! Reset the output variable
    likelihood = 0d0

    ! Begin looping through observations
    do n = 1, nobs
       ! Extract the time location of the current observation
       dn = obspts(n)
       ! Determine the lag period starting point
       s = max(1,dn-lag(dn))
       ! Estimate the mean over the lag period and accumulate the log-likelihood score
       likelihood = likelihood + &
                    ( ((sum(mod(s:dn)) / dble(lag(dn)+1)) - obs(dn)) / unc(dn) ) ** 2
    end do
    ! Apply the appropriate scaling, flip sign and multiply by 0.5.
    ! The likelihood scores for each observation should be subject to *-0.5
    ! in the algebraic formulation of the cost function. To avoid repeat calculation
    ! it is applied here once per data stream
    likelihood = -0.5d0 * likelihood * scaling ! e.g. 1/dble(DATAin%nCwood_inc)

    ! check that log-likelihood is an actual number
    if (likelihood /= likelihood) then
       likelihood = log(infini)
    end if

  end function likelihood
  !
  !------------------------------------------------------------------
  !
end module model_likelihood_module
