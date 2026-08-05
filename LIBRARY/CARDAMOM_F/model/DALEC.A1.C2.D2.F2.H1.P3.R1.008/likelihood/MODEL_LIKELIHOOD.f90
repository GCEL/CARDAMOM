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
! of DALEC.A1.C2.D2.F2.H1.P3.R1.008 as a function of observations and ecological dynamical constraints.
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
  public :: model_likelihood, scaled_model_likelihood, edc_model_likelihood, model_sanity_check, sanity_check

  ! declare needed types
  type EDCDIAGNOSTICS
    integer :: EDC
    integer :: DIAG
    integer :: PASSFAIL(150) ! allow space for 150 possible checks
    integer :: nedc ! number of edcs being assessed
  end type

  ! Has the model sanity check been conducted yet?
  logical :: sanity_check = .false.

  contains
  !
  !------------------------------------------------------------------
  !
  subroutine edc_model_likelihood(PARS, ML_obs_out, ML_prior_out, thread_id)
    use cardamom_structures, only: DATAin
    use model_shared, only: PI
    use carbon_model_mod, only: carbon_model
    use carbon_model_memory, only: mVs

    ! Model likelihood function specifically intended for the determination of
    ! appropriate initial parameter choices, consistent with EDCs for this DALEC

    implicit none

    ! declare inputs
    double precision, dimension(PI%npars), intent(inout) :: PARS
    ! output
    double precision, intent(inout) :: ML_obs_out, ML_prior_out

    ! declare local variables
    integer ::  n
    double precision :: tot_exp, ML, EDC1, EDC2, infini
    
    integer, intent(in), optional:: thread_id

    type (EDCDIAGNOSTICS) :: EDCD

    double precision,dimension(datain%nodays, datain%nofluxes)::  M_FLUXES
    double precision, dimension((DATAin%nodays+1), DATAin%nopools):: M_POOLS
    double precision,dimension(datain%nodays, datain%nodiags)::  M_DIAGS

    ! if == 0 EDCs are checked only until the first failure occurs
    ! if == 1 then all EDCs are checked irrespective of whether or not one has failed
    EDCD%DIAG = 1
    ML_obs_out = 0d0 ; ML_prior_out = 0d0

    ! Perform a more aggressive sanity check which compares the bulk difference
    ! in all fluxes and pools from multiple runs of the same parameter set
    if (.not.sanity_check) call model_sanity_check(PARS, thread_id)

    ! call EDCs which can be evaluated prior to running the model
    call assess_EDC1(PARS,PI%npars, DATAIN%nomet, DATAIN%nodays, DATAin%meantemp, DATAin%meanrad, DATAIN%MET,  EDC1, EDCD)

    ! next need to run the model itself
    call carbon_model(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT &
                     ,M_FLUXES,M_POOLS,M_DIAGS &
                     ,DATAin%nopars,DATAin%nomet,DATAin%nopools &
                     ,DATAin%nofluxes,DATAin%nodiags, mVs(thread_id))

    ! assess post running EDCs
    call assess_EDC2(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools  &
                    ,DATAin%nodays,DATAin%nodiags,DATAin%deltat            &
                    ,DATAin%steps_per_year,PI%parmax,PARS,DATAin%MET       &
                    ,M_POOLS,M_FLUXES,M_DIAGS         &
                    ,DATAin%meantemp,EDC2, EDCD)

    ! calculate the likelihood
    tot_exp = sum(1d0-EDCD%PASSFAIL(1:EDCD%nedc))
!! DEBUG CODE    
!    tot_exp = 0d0
!    do n = 1, EDCD%nedc
!       tot_exp = tot_exp+(1d0-EDCD%PASSFAIL(n))
!       !if (EDCD%PASSFAIL(n) /= 1) print*,"failed edcs are: ", n
!       if (sum(1d0-EDCD%PASSFAIL(1:EDCD%nedc)) < 3 .and. EDCD%PASSFAIL(n) /= 1) print*,"failed edcs are: ", n
!    end do ! checking EDCs
!   ! for testing purposes, stop the model when start achieved
!    if (sum(EDCD%PASSFAIL) == 100) then
!        print*,"Found it" ; stop
!    endif
!! DEBUG CODE

    ! convert to a probability
    ML_obs_out = -5d0*tot_exp*DATAin%EDC

  end subroutine edc_model_likelihood
  !
  !------------------------------------------------------------------
  !
  subroutine model_sanity_check(PARS, thread_id)
    use cardamom_structures, only: DATAin
    use model_shared, only: PI
    use carbon_model_mod, only: carbon_model
    use carbon_model_memory, only: mVs

    ! Carries out multiple carbon model iterations using the same parameter set
    ! to ensure that model outputs are consistent between iterations, i.e. that
    ! the model is numerically secure. Reproducible outputs from the models is
    ! essential for successful mcmc anlaysis

    implicit none

    ! Arguments
    double precision, dimension(PI%npars), intent(in) :: PARS

    ! Local arguments
    integer :: i,t
    double precision, dimension((DATAin%nodays+1),DATAin%nopools) :: local_pools
    double precision, dimension(DATAin%nodays,DATAin%nofluxes) :: local_fluxes
    double precision, dimension(DATAin%nodays,DATAin%nodiags) :: local_diags
    double precision :: pool_error, flux_error, diags_error

    integer, intent(in), optional:: thread_id
    double precision,dimension(datain%nodays, datain%nofluxes)::  M_FLUXES
    double precision, dimension((DATAin%nodays+1), DATAin%nopools):: M_POOLS
    double precision,dimension(datain%nodays, datain%nodiags)::  M_DIAGS

    ! Run model

    print*,"sanity_check: carbon_model run 1"
    ! next need to run the model itself
    call carbon_model(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT &
                     ,M_FLUXES,M_POOLS,M_DIAGS &
                     ,DATAin%nopars,DATAin%nomet,DATAin%nopools &
                     ,DATAin%nofluxes,DATAin%nodiags, mVs(thread_id))
    print*,"sanity_check: carbon_model run 2"
    ! next need to run the model itself
    call carbon_model(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT &
                     ,local_fluxes,local_pools,local_diags &
                     ,DATAin%nopars,DATAin%nomet,DATAin%nopools &
                     ,DATAin%nofluxes,DATAin%nodiags, mVs(thread_id))
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
        print*,"Cumulative DIAGS error = ",diags_error
        do i = 1,DATAin%nofluxes
           print*,"Sum abs error over time: flux = ",i
           print*,sum(abs(M_FLUXES(:,i) - local_fluxes(:,i)))
        end do
        do i = 1, DATAin%nopools
           print*,"Sum abs error over time: pool = ",i
           print*,sum(abs(M_POOLS(:,i) - local_pools(:,i)))
        end do
        do i = 1, DATAin%nodiags
           print*,"Sum abs error over time: diags = ",i
           print*,sum(abs(M_DIAGS(:,i) - local_diags(:,i)))
        end do
        print*,"First time step for all fluxes in run 1"
        print*,local_fluxes(1,:)
        print*,"First time step for all fluxes in run 2"
        print*,M_FLUXES(1,:)
        stop
    end if

!    ! Commented out to limit error messages, but useful for diagnosis
!    do t = 1, DATAin%nodays
!       if (sum(abs(M_FLUXES(t,:) - local_fluxes(t,:))) > (tiny(0d0)*(DATAin%nofluxes))) then
!           print*,"Time step of mismatch = ",i
!           do i = 1, DATAin%nofluxes
!               print*,"Flux counter = ",i
!               print*,"Original run"
!               print*,local_fluxes(t,i)
!               print*,"Second run"
!               print*,M_FLUXES(t,i)
!           end do
!       end if
!       stop
!    end do

    ! Update the user
    print*,"Sanity check completed"

    ! Set Sanity check as completed
    sanity_check = .true.

  end subroutine model_sanity_check
  !
  !------------------------------------------------------------------
  !
  subroutine assess_EDC1(PARS, npars, nomet, nodays, meantemp, meanrad, met, EDC1, EDCD)

    ! subroutine assessed the current parameter sets for passing ecological and
    ! steady state contraints (Bloom & Williams 2015).

    implicit none

    ! declare input variables
    integer, intent(in) :: npars ! number of parameters
    integer, intent(in) :: nomet ! number of met drivers
    integer, intent(in) :: nodays
    double precision, intent(in) :: met(nomet,nodays) ! array of met drivers
    double precision, intent(out) :: EDC1    ! EDC1 flag
    double precision, dimension(npars), intent(in) :: PARS ! current parameter set
    double precision, intent(in) :: meantemp & ! mean temperature (k)
                                   ,meanrad    ! mean radiation (MJ.m-2.day-1)
    type (EDCDIAGNOSTICS), intent(inout) :: EDCD

    ! declare local variables
    integer :: n, DIAG

    ! set initial value
    EDC1 = 1
    DIAG = EDCD%DIAG

    ! set all EDCs to 1 (pass)
    EDCD%nedc = 150
    EDCD%PASSFAIL(1:EDCD%nedc) = 1

    !
    ! begin checking EDCs
    !

    ! Turnover of litter faster than turnover of som
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(9) > ((1d0-pars(1)) * pars(8)))) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(1) = 0
    endif

    ! litter2som greater than som to atm rate
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(9) > (pars(1)*pars(8)) ) then
       EDC1 = 0d0 ; EDCD%PASSFAIL(2) = 0
    endif

    ! turnover of litwood (pars(16)) should be slower than fine litter turnover pars(8)
    if ((EDC1 == 1 .or. DIAG == 1) .and. ( pars(16) > pars(8) ) ) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(3) = 0
    endif

    ! root turnover greater than som turnover at mean temperature
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(7) < (pars(9)*exp(pars(10)*meantemp)))) then
       EDC1 = 0d0 ; EDCD%PASSFAIL(4) = 0
    endif

    ! Initial leaf area index should not be larger than ~10 m2/m2
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(19)/pars(17)) > 10d0) then
        EDC1 = 0d0 ; EDCD%PASSFAIL(5) = 0
    endif    

    ! Straight forward GSI parameter bounds
    ! Temperature
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(35) < pars(34))) then
         EDC1 = 0d0 ; EDCD%PASSFAIL(6) = 0
    end if
    ! Photoperiod
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(37) < pars(36))) then
         EDC1 = 0d0 ; EDCD%PASSFAIL(7) = 0
    end if
    ! VPD - may not be needed
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(39) < pars(38))) then
         EDC1 = 0d0 ; EDCD%PASSFAIL(8) = 0
    end if

    ! Photoperiod minimum cannot be substantially less than the observed minimum day length
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(36) < minval(MET(11,:))-14400d0)) then
         EDC1 = 0d0 ; EDCD%PASSFAIL(9) = 0
    end if
    ! Photoperiod maximum cannot be greater than the observed maximum day length
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(37) > maxval(MET(11,:)))) then
         EDC1 = 0d0 ; EDCD%PASSFAIL(10) = 0
    end if

    ! VPD at which stress in at maximum should be no larger than max(VPDlag21) +
    ! 1500 Pa from the max VPD tolerated parameter. The hypothesis being that plants would not
    ! make themselves resilient to an unexperiences environment.
    if ((EDC1 == 1 .or. DIAG == 1) .and. (pars(39) > maxval(MET(12,:))+1500d0)) then
         EDC1 = 0d0 ; EDCD%PASSFAIL(11) = 0
    end if

    ! IMPLICIT Combustion completeness for foliage should be greater than soil
    ! IMPLICIT Combustion completeness for wood litter and fol+root litter should be greater than soil

    ! Combustion completeness for foliage should be greater than non-photosynthetic tissues
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(29) < pars(30)) then
       EDC1 = 0d0 ; EDCD%PASSFAIL(12) = 0
    endif
    ! Combustion completeness for non-photosynthetic tissue should be greater than soil
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(30) < pars(31)) then
       EDC1 = 0d0 ; EDCD%PASSFAIL(13) = 0
    endif
    ! Combustion completeness for foliar + fine root litter should be greater than foliage
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(32) < pars(29)) then
       EDC1 = 0d0 ; EDCD%PASSFAIL(14) = 0
    endif
    ! Combustion completeness for wood litter should be greater than non-photosynthetic tissue
    if ((EDC1 == 1 .or. DIAG == 1) .and. pars(33) < pars(30)) then
       EDC1 = 0d0 ; EDCD%PASSFAIL(15) = 0
    endif

    ! could always add more / remove some

  end subroutine assess_EDC1
  !
  !------------------------------------------------------------------
  !
  subroutine assess_EDC2(npars,nomet,nofluxes,nopools,nodays,nodiags,deltat,steps_per_year &
                        ,parmax,pars,met,M_POOLS,M_FLUXES,M_DIAGS &
                        ,meantemp,EDC2, EDCD)

    use cardamom_structures, only: DATAin

    ! Determines whether the dynamical contraints for the search of the initial
    ! parameters has been successful or whether or not we should abandon the
    ! current set and move on

    implicit none

    ! declare input variables
    integer, intent(in) :: npars          & ! number of model parameters
                          ,nomet          & ! number of met drivers
                          ,nofluxes       & ! number of fluxes from model
                          ,nopools        & ! number of pools in model
                          ,nodays         & ! number of days in simulation
                          ,nodiags        & ! number of diagnostic variables                          
                          ,steps_per_year

    double precision, intent(in) :: deltat(nodays)              & ! decimal day model interval
                                   ,pars(npars)                 & ! vector of current parameters
                                   ,parmax(npars)               & ! vector of the maximum parameter values
                                   ,met(nomet,nodays)           & ! array of met drivers
                                   ,M_POOLS((nodays+1),nopools) & ! time varying states of pools in current model simulation
                                   ,M_FLUXES(nodays,nofluxes)   & ! time varying fluxes from current model simulation model
                                   ,M_DIAGS(nodays,nodiags)     & ! time varying diagnostics from current model simulation model                                   
                                   ,meantemp                      ! site mean temperature (oC)

    double precision, intent(out) :: EDC2 ! the response flag for the dynamical set of EDCs
    type (EDCDIAGNOSTICS), intent(inout) :: EDCD

    ! declare local variables
    integer :: n, nn, nnn, DIAG, y, steps_per_month, nd, fl, fs, &
               io_start, io_finish
    double precision :: infi, tmp, tmp1, tmp2, &!, EQF, etol
                        jan_sd_lai, jan_mean_lai, jan_first_lai, &
                        SSwood, SSlitwood, SSsom
    !double precision, dimension(nodays) :: tmp1, tmp2
    double precision, dimension(nopools) :: jan_mean_pools, jan_first_pools, &
                                            mean_pools, Fin, Fout, Rm, Rs, &
                                            Fin_yr1, Fout_yr1, Fin_yr2, Fout_yr2
    double precision, dimension(nofluxes) :: FT, FT_yr1, FT_yr2

    ! Steady State Attractor:
    ! Log ratio difference between inputs and outputs of the system.
    double precision, parameter :: EQF1_5 = log(1.5d0), & ! 10.0 = order magnitude; 2 = double and half
                                   EQF2 = log(2d0),   & ! 10.0 = order magnitude; 2 = double and half
                                   EQF5 = log(5d0),   &
                                   EQF10 = log(10d0), &
                                   EQF15 = log(15d0), &
                                   EQF20 = log(20d0), &
                                  C_etol = 0.10d0,    & ! 0.20d0 lots of AGB !0.10d0 global / site more data !0.05d0 global 1 or 2 AGB estimates
                                H2O_etol = 0.05d0       !

    ! update initial values
    DIAG = EDCD%DIAG
    EDC2 = 1
    infi = 0d0

    ! derive mean pools
    do n = 1, nopools
       mean_pools(n) = cal_mean_pools(M_POOLS(1:nodays,n),nodays)
    end do

    ! number of time steps per month
    steps_per_month = ceiling(dble(steps_per_year) * 0.08333333d0)

    ! First calculate total flux for the simulation period
    io_start = (steps_per_year*2) + 1 ; io_finish = nodays
    if (DATAin%nos_years < 3) io_start = 1
    do fl = 1, nofluxes
!       FT(fl) = sum(M_FLUXES(1:nodays,fl)*deltat(1:nodays))
       FT(fl) = sum(M_FLUXES(io_start:io_finish,fl)*deltat(io_start:io_finish))
       FT_yr1(fl) = sum(M_FLUXES(1:steps_per_year,fl)*deltat(1:steps_per_year))
       !FT_yr2(fl) = sum(M_FLUXES((steps_per_year+1):(steps_per_year*2),fl) &
                       !*deltat((steps_per_year+1):(steps_per_year*2)))
    end do
    ! Specific calculation of transpiration extraction from the soil surface layer
    fl = 46 ! transpiration multiplied by ...
    fs = 53 ! ...fraction of transpiration extracted from 1st rooting layer (the soil surface)
    FT(fl) = sum(M_FLUXES(io_start:io_finish,fl)*M_FLUXES(io_start:io_finish,fs)*deltat(io_start:io_finish))
    FT_yr1(fl) = sum(M_FLUXES(1:steps_per_year,fl)*M_FLUXES(1:steps_per_year,fs)*deltat(1:steps_per_year))
    !FT_yr2(fl) = sum(M_FLUXES((steps_per_year+1):(steps_per_year*2),fl) & 
                    !*M_FLUXES((steps_per_year+1):(steps_per_year*2),fs) &
                    !*deltat((steps_per_year+1):(steps_per_year*2)))

    ! get total in and out for each pool
    ! labile
    Fin(1)  = FT(5)
    Fout(1) = FT(8)+FT(18)+FT(24)+FT(33)+FT(40)+FT(55)
    Fin_yr1(1)  = FT_yr1(5)
    Fout_yr1(1) = FT_yr1(8)+FT_yr1(18)+FT_yr1(24)+FT_yr1(33)+FT_yr1(40)+FT_yr1(55)
!    Fin_yr2(1)  = FT_yr2(5)
!    Fout_yr2(1) = FT_yr2(8)+FT_yr2(18)+FT_yr2(24)+FT_yr2(33)+FT_yr2(40)+FT_yr2(55)
    ! foliar
    Fin(2)  = FT(8)
    Fout(2) = FT(10)+FT(19)+FT(25)+FT(34)+FT(41)
    Fin_yr1(2)  = FT_yr1(4)+FT_yr1(8)
    Fout_yr1(2) = FT_yr1(10)+FT_yr1(19)+FT_yr1(25)+FT_yr1(34)+FT_yr1(41)
 !   Fin_yr2(2)  = FT_yr2(4)+FT_yr2(8)
 !   Fout_yr2(2) = FT_yr2(10)+FT_yr2(19)+FT_yr2(25)+FT_yr2(34)+FT_yr2(41)
    ! root
    Fin(3)  = FT(6)
    Fout(3) = FT(12)+FT(20)+FT(26)+FT(35)+FT(42)
    Fin_yr1(3)  = FT_yr1(6)
    Fout_yr1(3) = FT_yr1(12)+FT_yr1(20)+FT_yr1(26)+FT_yr1(35)+FT_yr1(42)
 !   Fin_yr2(3)  = FT_yr2(6)
 !   Fout_yr2(3) = FT_yr2(12)+FT_yr2(20)+FT_yr2(26)+FT_yr2(35)+FT_yr2(42)
    ! wood
    Fin(4)  = FT(7)
    Fout(4) = FT(11)+FT(21)+FT(27)+FT(36)+FT(43)
    Fin_yr1(4)  = FT_yr1(7)
    Fout_yr1(4) = FT_yr1(11)+FT_yr1(21)+FT_yr1(27)+FT_yr1(36)+FT_yr1(43)
 !   Fin_yr2(4)  = FT_yr2(7)
 !   Fout_yr2(4) = FT_yr2(11)+FT_yr2(21)+FT_yr2(27)+FT_yr2(36)+FT_yr2(43)
    ! litter
    Fin(5)  = FT(10)+FT(12)+FT(24)+FT(25)+FT(26)+FT(40)+FT(41)+FT(42)
    Fout(5) = FT(13)+FT(15)+FT(22)+FT(28)+FT(37)
    Fin_yr1(5)  = FT_yr1(10)+FT_yr1(12)+FT_yr1(24)+FT_yr1(25)+FT_yr1(26)+FT_yr1(40)+FT_yr1(41)+FT_yr1(42)
    Fout_yr1(5) = FT_yr1(13)+FT_yr1(15)+FT_yr1(22)+FT_yr1(28)+FT_yr1(37)
!    Fin_yr2(5)  = FT_yr2(10)+FT_yr2(12)+FT_yr2(24)+FT_yr2(25)+FT_yr2(26)+FT_yr2(40)+FT_yr2(41)+FT_yr2(42)
!    Fout_yr2(5) = FT_yr2(13)+FT_yr2(15)+FT_yr2(22)+FT_yr2(28)+FT_yr2(37)
    ! som
    Fin(6)  = FT(15)+FT(27)+FT(28)+FT(31)+FT(45)
    Fout(6) = FT(14)+FT(23)+FT(39)
    Fin_yr1(6)  = FT_yr1(15)+FT_yr1(27)+FT_yr1(28)+FT_yr1(31)+FT_yr1(45)
    Fout_yr1(6) = FT_yr1(14)+FT_yr1(23)+FT_yr1(39)
!    Fin_yr2(6)  = FT_yr2(15)+FT_yr2(27)+FT_yr2(28)+FT_yr2(31)+FT_yr2(45)
!    Fout_yr2(6) = FT_yr2(14)+FT_yr2(23)+FT_yr2(39)
    ! litwood
    Fin(7)  = FT(11)+FT(43)
    Fout(7) = FT(30)+FT(31)+FT(38)+FT(44)+FT(45)
    Fin_yr1(7)  = FT_yr1(11)+FT_yr1(43)
    Fout_yr1(7) = FT_yr1(30)+FT_yr1(31)+FT_yr1(38)+FT_yr1(44)+FT_yr1(45)
!    Fin_yr2(7)  = FT_yr2(11)+FT_yr2(43)
!    Fout_yr2(7) = FT_yr2(30)+FT_yr2(31)+FT_yr2(38)+FT_yr2(44)+FT_yr2(45)

!    ! Surface water pool (0-30cm)
!    ! 52 = infiltrated, 47 = soil evaporation, 46 = transpiration from top soil, 51 = drainage from top soil
!    Fin(8)  = FT(52) 
!    Fout(8) = FT(47)+FT(46)+FT(51) 
!    Fin_yr1(8)  = FT_yr1(52) 
!    Fout_yr1(8) = FT_yr1(47)+FT_yr1(46)+FT_yr1(51) 
!    Fin_yr2(8)  = FT_yr2(52)
!    Fout_yr2(8) = FT_yr2(47)+FT_yr2(46)+FT_yr2(51)

    !
    ! Begin EDCs here
    !

    if (DATAin%nos_years > 1) then
        ! Determine the mean and standard deviation of January LAIs 
        jan_sd_lai = 0d0 ; jan_mean_lai = 0d0 ; jan_first_lai = 0d0 ! reset 
        jan_first_lai = M_DIAGS(1,1) ! First January LAI
        ! Initially sum each January from each year
        do y = 1, DATAin%nos_years
           nn = 1 + (steps_per_year * (y - 1)) 
           jan_mean_lai = jan_mean_lai + M_DIAGS(nn,1)
        end do
        ! Calculate the mean
        jan_mean_lai = jan_mean_lai / dble(DATAin%nos_years)
        ! Calculate the standard deviation now
        do y = 1, DATAin%nos_years
           nn = 1 + (steps_per_year * (y - 1)) 
           jan_sd_lai = jan_sd_lai + (jan_mean_lai - M_DIAGS(nn,1))**2d0
        end do
        jan_sd_lai = sqrt(jan_sd_lai / (dble(DATAin%nos_years - 1)))
        if ((EDC2 == 1 .or. DIAG == 1) .and. &
            abs(jan_first_lai-jan_mean_lai) > (jan_sd_lai*2d0) .and. abs(jan_first_lai-jan_mean_lai) > 0.01d0) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(16) = 0
        end if
    end if ! nos_years > 1

    if (DATAin%nos_years > 1) then
        ! Determine the mean and standard deviation of January wSWPs 
        jan_sd_lai = 0d0 ; jan_mean_lai = 0d0 ; jan_first_lai = 0d0 ! reset 
        jan_first_lai = M_DIAGS(1,10) ! First January wSWP
        ! Initially sum each January from each year
        do y = 1, DATAin%nos_years
           nn = 1 + (steps_per_year * (y - 1)) 
           jan_mean_lai = jan_mean_lai + M_DIAGS(nn,10)
        end do
        ! Calculate the mean
        jan_mean_lai = jan_mean_lai / dble(DATAin%nos_years)
        ! Calculate the standard deviation now
        do y = 1, DATAin%nos_years
           nn = 1 + (steps_per_year * (y - 1)) 
           jan_sd_lai = jan_sd_lai + (jan_mean_lai - M_DIAGS(nn,10))**2d0
        end do
        jan_sd_lai = sqrt(jan_sd_lai / (dble(DATAin%nos_years - 1)))
        if ((EDC2 == 1 .or. DIAG == 1) .and. &
            abs(jan_first_lai-jan_mean_lai) > (jan_sd_lai*2d0) .and. abs(jan_first_lai-jan_mean_lai) > 0.01d0) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(17) = 0
        end if
    end if ! nos_years > 1

    ! EDC just for DALEC_CDEA_ACM2_BUCKET due to complications linked to
    ! the empirical phenology but mechanistic hydrology / photosynthesis
    if ((EDC2 == 1 .or. DIAG == 1) .and. maxval(M_DIAGS(1:nodays,1)) > 10d0 ) then
        EDC2 = 0d0 ; EDCD%PASSFAIL(18) = 0
    end if

    ! What are in effect the potential growth rates are modulated by the current 
    ! fixed temperature sub-model used in the model. This means that the parameterised 
    ! potential rates might never be achievable even if plausible. Thus the maximum 
    ! parameter bound for the potential growth rates need to be increased. These EDCs 
    ! prevent an emergent growth rate that is unrealistic. Here we assume that tissue 
    ! growth for foliage, wood and roots cannot be greater than 10 gC/m2/day
    if ((EDC2 == 1 .or. DIAG == 1)) then
        ! Foliage
        if (maxval(M_FLUXES(:,8)) > 10d0) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(19) = 0
        end if
        ! Fine roots
        if (maxval(M_FLUXES(:,6)) > 10d0) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(20) = 0
        end if
        ! Wood
        if (maxval(M_FLUXES(:,7)) > 10d0) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(21) = 0
        end if
    end if

    ! Average growth rates for foliage and fine roots cannot be 5 orders of magnitude different
    if ((EDC2 == 1 .or. DIAG == 1) .and. FT(8) > (5d0*FT(6))) then
        EDC2 = 0d0 ; EDCD%PASSFAIL(22) = 0
    endif
    ! Average growth rates for foliage and fine roots cannot be 5 orders of magnitude different
    if ((EDC2 == 1 .or. DIAG == 1) .and. (FT(8)*5d0) < FT(6)) then
        EDC2 = 0d0 ; EDCD%PASSFAIL(23) = 0
    endif

    if (EDC2 == 1 .or. DIAG == 1) then

!        ! Living pools
!        do n = 1, 3
!           ! Restrict mean rates of increase
!           if (abs(log(Fin(n)/Fout(n))) > EQF2) then
!               EDC2 = 0d0 ; EDCD%PASSFAIL(20+n-1) = 0
!           end if
!           ! Restrict rates from deviating unrealistically from the mean
!           if ( abs( abs(log(Fin_yr1(n)/Fout_yr1(n))) - &
!                     abs(log(Fin(n)/Fout(n))) ) > C_etol ) then
!               EDC2 = 0d0 ; EDCD%PASSFAIL(30+n-1) = 0
!           end if
!        end do
!        ! Foliage pool, note that in CDEA EDCs Fin has already been multiplied by time step
!        n = 2
!        ! Restrict mean rates of increase
!        if (abs(log(Fin(n)/Fout(n))) > EQF2) then
!            EDC2 = 0d0 ; EDCD%PASSFAIL(20+n-1) = 0
!        end if
!        ! Restrict exponential behaviour at initialisation         
!        if (abs(abs(log(Fin_yr1(n)/Fout_yr1(n))) - abs(log(Fin_yr2(n)/Fout_yr2(n)))) > C_etol) then
!            EDC2 = 0d0 ; EDCD%PASSFAIL(20+n-1) = 0
!        end if
        ! Fine root pool, note that in CDEA EDCs Fin has already been multiplied by time step
        n = 3
!        ! Restrict mean rates of increase
!        if (abs(log(Fin(n)/Fout(n))) > EQF2) then
!            EDC2 = 0d0 ; EDCD%PASSFAIL(20+n-1) = 0
!        end if
        ! Restrict rates from deviating unrealistically from the mean
        if ( abs( abs(log(Fin_yr1(n)/Fout_yr1(n))) - &
                  abs(log(Fin(n)/Fout(n))) ) > C_etol ) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(30+n-1) = 0
        end if
!        ! Wood pool hack, note that in CDEA EDCs Fin has already been multiplied by time step
!        n = 4
!        if (abs(log(Fin(n)/Fout(n))) > EQF2) then
!            EDC2 = 0d0 ; EDCD%PASSFAIL(20+n-1) = 0
!        end if
!        if ( abs( abs(log(Fin_yr1(n)/Fout_yr1(n))) - &
!                  abs(log(Fin(n)/Fout(n))) ) > C_etol ) then
!            EDC2 = 0d0 ; EDCD%PASSFAIL(30+n-1) = 0
!        end if
        ! Dead pools
        do n = 5, 7
           ! Restrict rates of increase
           if (abs(log(Fin(n)/Fout(n))) > EQF1_5) then
               EDC2 = 0d0 ; EDCD%PASSFAIL(20+n-1) = 0
           end if
           ! Restrict rates from deviating unrealistically from the mean
           if ( abs( abs(log(Fin_yr1(n)/Fout_yr1(n))) - &
                     abs(log(Fin(n)/Fout(n))) ) > C_etol ) then
               EDC2 = 0d0 ; EDCD%PASSFAIL(30+n-1) = 0
           end if
        end do

        ! Determine the steady state estimate of wood (gC/m2)
        SSwood = (Fin(4)/Fout(4)) * mean_pools(4)
        ! Based on the wood SS (gC/m2) and the sum fractional loss per day determine the mean input to litwood...
        SSlitwood = SSwood * (Fout(4)/mean_pools(4))
        ! ...then estimate the actual steady state wood litter
        SSlitwood = (SSlitwood/Fout(7)) * mean_pools(7)
        ! Steady state of som requires accounting for foliar, fine root and wood litter inputs
        ! and adjusting for the litwood input already included
        SSsom = Fin(6) - sum(M_FLUXES(io_start:io_finish,31))
        ! Now repeat the process as done for litwood to estimate the inputs,
        ! adjusting for the fraction of litwood output which is respired not decomposed
        SSsom = SSsom + (SSlitwood * (Fout(7)/mean_pools(7)) * pars(1))
        ! Accounting for losses and scaling to SSsom
        SSsom = (SSsom / Fout(6)) * mean_pools(6)
        ! It is reasonable to assume that the steady state for woody litter
        ! should be ~ less than half that of woody biomass...
        if (SSlitwood / SSwood > 0.60d0  ) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(45) = 0
        end if
        ! ... and less than soil organic matter
        if ( SSsom < SSlitwood ) then
            EDC2 = 0d0 ; EDCD%PASSFAIL(46) = 0
        end if

    end if ! EDC2 == 1 .or. DIAG == 1

    ! Ensure that the mean transit time of foliage and the LCA are consistent with the 
    ! leaf economic spectrum (LES).
    ! LL (months) ~ LMA (gm2) R2 = 0.42 from 
    ! Wright et al., (2004), doi: https://doi.org/10.1038/nature02403
    ! Onoda et al., (2017), doi: https://doi.org/10.1111/nph.14496 
    if (EDC2 == 1 .or. DIAG == 1) then
        ! Assume that the MTT(nat,fire) foliage should be within the uncertainty bounds of the LES
        ! Mean equation LL(months) = 0.0031 * LMA**1.71, coefficient 95CI = 1.62,1.82
        ! Estimating the MTT, converting from days to years using 1/365.25 = 0.002737851
        ! 0.08333333 converts months to years for the LES equation.
        tmp = sum(M_POOLS(:,2)) / dble(nodays)
        tmp1 = sum(M_FLUXES(:,10)+M_FLUXES(:,19)+M_FLUXES(:,25)) / dble(nodays) ! Foliar outputs, not including deforestation
        tmp = (tmp / tmp1) * 0.002737851d0
        ! determine the lower and upper bound of the LES .
        ! not for the upper bound, do not allow a value less than 1 years
        tmp1 = 0.08333333d0*(0.0031d0*(pars(17)*2.083333d0)**1.62d0)
        tmp2 = max(1d0,0.08333333d0*(0.0031d0*(pars(17)*2.083333d0)**1.82d0))
        if (tmp < tmp1) then
            ! The current leaf lifespan is shorter than expected
            EDC2 = 0d0 ; EDCD%PASSFAIL(47) = 0
        endif        
        if (tmp > tmp2) then
            ! The current leaf life span is longer than expected
            EDC2 = 0d0 ; EDCD%PASSFAIL(48) = 0
        endif        
    endif ! EDC2 == 1 .or. DIAG == 1

    !
    ! EDCs done, below are additional fault detection conditions
    !

    ! additional faults can be stored in locations 55 - 61 of the PASSFAIL array

    ! ensure minimum pool values are >= 0, /= NaN or Inf
    if (EDC2 == 1 .or. DIAG == 1) then

       do n = 1, nopools
          if (minval(M_POOLS(1:nodays,n)) < 0d0 .or. &
              maxval(abs(M_POOLS(1:nodays,n))) == abs(log(infi)) .or. &
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
  double precision function cal_mean_pools(pools,averaging_period)

    ! Function calculate the mean values of model pools / states

    implicit none

    ! declare input variables
    integer, intent(in) :: averaging_period   !

    double precision,dimension(averaging_period), intent (in) :: pools

    ! declare local variables
    integer :: c

    ! loop through now
    cal_mean_pools = sum(pools(1:averaging_period))/dble(averaging_period)

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
  subroutine model_likelihood(PARS,ML_obs_out,ML_prior_out, thread_id)
    use model_shared, only:  PI
    use carbon_model_mod, only: carbon_model
    use carbon_model_memory, only: mVs
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

    type(EDCDIAGNOSTICS) :: EDCD
    integer, intent(in), optional:: thread_id
    double precision,dimension(datain%nodays, datain%nofluxes)::  M_FLUXES
    double precision, dimension((DATAin%nodays+1), DATAin%nopools):: M_POOLS
    double precision,dimension(datain%nodays, datain%nodiags)::  M_DIAGS

    ! initial values
    ML_obs_out = 0d0 ; ML_prior_out = 0d0 ; EDC1 = 1d0 ; EDC2 = 1d0
    ! if == 0 EDCs are checked only until the first failure occurs
    ! if == 1 then all EDCs are checked irrespective of whether or not one has failed
    EDCD%DIAG = 0

    if (DATAin%EDC == 1) then

        ! call EDCs which can be evaluated prior to running the model
        call assess_EDC1(PARS,PI%npars, DATAIN%nomet, DATAIN%nodays, DATAin%meantemp, DATAin%meanrad, DATAIN%MET,  EDC1, EDCD)

        ! update the likelihood score based on EDCs driving total rejection
        ! proposed parameters
        ML_obs_out = log(EDC1)

    endif !

    ! run the dalec model
    call carbon_model(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT &
                     ,M_FLUXES,M_POOLS,M_DIAGS &
                     ,DATAin%nopars,DATAin%nomet,DATAin%nopools &
                     ,DATAin%nofluxes,DATAin%nodiags, mVs(thread_id))

    ! if first set of EDCs have been passed, move on to the second
    if (DATAin%EDC == 1) then

        ! check edc2
        call assess_EDC2(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools  &
                        ,DATAin%nodays,DATAin%nodiags,DATAin%deltat            &
                        ,DATAin%steps_per_year,PI%parmax,PARS,DATAin%MET       &
                        ,M_POOLS,M_FLUXES,M_DIAGS         &
                        ,DATAin%meantemp,EDC2, EDCD)

        ! Add EDC2 log-likelihood to absolute accept reject...
        ML_obs_out = ML_obs_out + log(EDC2)

    end if ! DATAin%EDC == 1

    ! Calculate log-likelihood associated with priors
    ! We always want this
    ML_prior_out = likelihood_p(PI%npars,DATAin%parpriors,DATAin%parpriorunc,DATAin%parpriorweight,PARS)
    ! Calculate log-likelihood of model compared to obs
    call calc_obs_likelihoods(ML_obs_out, M_POOLS, M_FLUXES, M_DIAGS)
    ! Calculate log-likelihood of 'other priors'
    call calc_other_likelihoods(ML_obs_out, M_POOLS, M_FLUXES, M_DIAGS)

  end subroutine model_likelihood
  !
  !------------------------------------------------------------------
  !
  subroutine scaled_model_likelihood(PARS,ML_obs_out,ML_prior_out, thread_id)
    use model_shared, only:  PI
    use carbon_model_mod, only: carbon_model
    use carbon_model_memory, only: mVs
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
    type(EDCDIAGNOSTICS) :: EDCD
    double precision :: EDC1, EDC2
    double precision,dimension(datain%nodays, datain%nofluxes)::  M_FLUXES
    double precision, dimension((DATAin%nodays+1), DATAin%nopools):: M_POOLS
    double precision,dimension(datain%nodays, datain%nodiags)::  M_DIAGS
    integer, intent(in), optional :: thread_id

    ! initial values
    ML_obs_out = 0d0 ; ML_prior_out = 0d0 ; EDC1 = 1d0 ; EDC2 = 1d0
    ! if == 0 EDCs are checked only until the first failure occurs
    ! if == 1 then all EDCs are checked irrespective of whether or not one has failed
    EDCD%DIAG = 0

    if (DATAin%EDC == 1) then

        ! call EDCs which can be evaluated prior to running the model
        call assess_EDC1(PARS,PI%npars, DATAIN%nomet, DATAIN%nodays, DATAin%meantemp, DATAin%meanrad, DATAIN%MET,  EDC1, EDCD)

        ! update the likelihood score based on EDCs driving total rejection
        ! proposed parameters
        ML_obs_out = log(EDC1)
    endif !

    ! run the dalec model
    call carbon_model(1,DATAin%nodays,DATAin%MET,PARS,DATAin%deltat &
                     ,DATAin%nodays,DATAin%LAT &
                     ,M_FLUXES,M_POOLS,M_DIAGS &
                     ,DATAin%nopars,DATAin%nomet,DATAin%nopools &
                     ,DATAin%nofluxes,DATAin%nodiags, mVs(thread_id))

    ! if first set of EDCs have been passed, move on to the second
    if (DATAin%EDC == 1) then

        ! check edc2
        call assess_EDC2(PI%npars,DATAin%nomet,DATAin%nofluxes,DATAin%nopools  &
                        ,DATAin%nodays,DATAin%nodiags,DATAin%deltat            &
                        ,DATAin%steps_per_year,PI%parmax,PARS,DATAin%MET       &
                        ,M_POOLS,M_FLUXES,M_DIAGS         &
                        ,DATAin%meantemp,EDC2, EDCD)      

        ! Add EDC2 log-likelihood to absolute accept reject...
        ML_obs_out = ML_obs_out + log(EDC2)
    end if ! DATAin%EDC == 1

    ! Calculate log-likelihood associated with priors
    ! We always want this
    ML_prior_out = likelihood_p(PI%npars,DATAin%parpriors,DATAin%parpriorunc,DATAin%parpriorweight,PARS)
    ! Calculate log-likelihood of model compared to obs
    call calc_scaled_obs_likelihoods(ML_obs_out, M_POOLS, M_FLUXES, M_DIAGS)
    ! Calculate log-likelihood of 'other priors'
    call calc_other_likelihoods(ML_obs_out, M_POOLS, M_FLUXES, M_DIAGS)

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
!print*,"likelihood_p:"
    ! set initial value
    likelihood_p = 0d0 ; local_likelihood = 0d0

    ! now loop through defined parameters for their uncertainties
    where (parpriors > -9999) local_likelihood = parpriorweight*((pars-parpriors)/parpriorunc)**2
    likelihood_p = sum(local_likelihood) * (-0.5d0)
!print*,"likelihood_p: done"
    ! dont for get to return
    return

  end function likelihood_p
  !
  !------------------------------------------------------------------
  !
  subroutine calc_obs_likelihoods(ML_obs_out, M_POOLS, M_FLUXES, M_DIAGS)
    use cardamom_structures, only: DATAin
    use carbon_model_memory, only: sw_par_fraction, top_soil_depth

    ! Subroutine to control the calculation of the observation related
    ! log-likelihoods and accounting for the various scalings

    ! Arguements
    double precision, intent(inout) :: ML_obs_out
    double precision, dimension(datain%nodays, datain%nofluxes)::  M_FLUXES
    double precision, dimension((DATAin%nodays+1), DATAin%nopools):: M_POOLS
    double precision, dimension(datain%nodays, datain%nodiags)::  M_DIAGS
    ! local variable
    double precision, dimension(DATAin%nodays) :: mod

    !
    ! Do diagnostics (DIAGS)
    ! 

    ! Calculate log-likelihood for fraction of absorbed radiation
    if (DATAin%nfAPAR > 0) then
        mod = M_DIAGS(1:DATAin%nodays,3) & ! APAR
            / (DATAin%met(4,1:DATAin%nodays)*sw_par_fraction)
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nfAPAR,DATAin%fAPARpts,DATAin%fAPAR,DATAin%fAPAR_unc,DATAin%fAPAR_lag, &
                                             1d0,mod)
    endif ! nfAPAR > 0
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
    ! Calculate log-likelihood for total wood stocks
    if (DATAin%nCwood_stock > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCwood_stock,DATAin%Cwood_stockpts, &
                                             DATAin%Cwood_stock,DATAin%Cwood_stock_unc,DATAin%Cwood_stock_lag, &
                                             1d0,M_POOLS(1:DATAin%nodays,4))
    endif ! nCwood_stock > 0
    ! Calculate log-likelihood for foliage litter stocks
    if (DATAin%nClit_stock > 0) then
        ! Estimate the foliage litter pool based on the ratio of foliage litter input to foliage + fine root litter inputs,
        ! scaled by the total litter pool. This is based on the turnover being common.
        mod = (sum(M_FLUXES(1:DATAin%nodays,10))/sum(M_FLUXES(1:DATAin%nodays,10)+M_FLUXES(1:DATAin%nodays,12))) &
                 * M_POOLS(1:DATAin%nodays,5)
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nClit_stock,DATAin%Clit_stockpts, &
                                             DATAin%Clit_stock,DATAin%Clit_stock_unc,DATAin%Clit_stock_lag, &
                                             1d0,mod)
    endif ! nClit_stock > 0
    ! Calculate log-likelihood for soil organic matter stocks
    if (DATAin%nCsom_stock > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCsom_stock,DATAin%Csom_stockpts, &
                                             DATAin%Csom_stock,DATAin%Csom_stock_unc,DATAin%Csom_stock_lag, &
                                             1d0,M_POOLS(1:DATAin%nodays,6))
    endif ! nCsom_stock > 0
    ! Calculate log-likelihood for surface soil water
!    if (DATAin%nsoilwater > 0) then
!        mod = (M_POOLS(1:DATAin%nodays,8) * 1d-3) / top_soil_depth ! convert mm -> m3/m3
!        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nsoilwater,DATAin%soilwaterpts, &
!                                             DATAin%soilwater,DATAin%soilwater_unc,DATAin%soilwater_lag, &
!                                             1d0,mod)
!    endif ! nsoilwater > 0

    !
    ! Do fluxes (FLUXES)
    !

    ! Calculate log-likelihood for evapotranspiration
    if (DATAin%nEvap > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nEvap,DATAin%Evappts,DATAin%Evap,DATAin%Evap_unc,DATAin%Evap_lag, &
                                             1d0,M_FLUXES(1:DATAin%nodays,29))
    endif ! nEvap > 0
    ! Calculate log-likelihood for fire
    if (DATAin%nFire > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nFire,DATAin%Firepts,DATAin%Fire,DATAin%Fire_unc,DATAin%Fire_lag, &
                                             1d0,M_FLUXES(1:DATAin%nodays,17))
    endif ! nFire > 0
    ! Calculate log-likelihood for gross primary production
    if (DATAin%ngpp > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%ngpp,DATAin%gpppts,DATAin%GPP,DATAin%GPP_unc,DATAin%GPP_lag, &
                                             1d0,M_FLUXES(1:DATAin%nodays,1))
    endif ! ngpp > 0
    ! Calculate log-likelihood for harvest
    if (DATAin%nharvest > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nharvest,DATAin%harvestpts,DATAin%harvest,DATAin%harvest_unc,DATAin%harvest_lag, &
                                             1d0,M_FLUXES(1:DATAin%nodays,32))
    endif ! nharvest > 0
    ! Calculate log-likelihood for net biome exchange 
    if (DATAin%nnbe > 0) then
        mod = M_FLUXES(1:DATAin%nodays,3) &  ! Rauto
            + M_FLUXES(1:DATAin%nodays,13) & ! Rhet litter
            + M_FLUXES(1:DATAin%nodays,30) & ! Rhet wood litter
            + M_FLUXES(1:DATAin%nodays,14) & ! Rhet som
            + M_FLUXES(1:DATAin%nodays,17) & ! Fire
            - M_FLUXES(1:DATAin%nodays,1)    ! GPP
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nnbe,DATAin%nbepts,DATAin%NBE,DATAin%NBE_unc,DATAin%NBE_lag, &
                                             1d0,mod)
    endif ! nnbe > 0
    ! Calculate log-likelihood for net ecosystem exchange of CO2
    if (DATAin%nnee > 0) then
        mod = M_FLUXES(1:DATAin%nodays,3) &  ! Rauto
            + M_FLUXES(1:DATAin%nodays,13) & ! Rhet litter
            + M_FLUXES(1:DATAin%nodays,30) & ! Rhet wood litter
            + M_FLUXES(1:DATAin%nodays,14) & ! Rhet som
            - M_FLUXES(1:DATAin%nodays,1)    ! GPP
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nnee,DATAin%neepts,DATAin%NEE,DATAin%NEE_unc,DATAin%NEE_lag, &
                                             1d0,mod)    
    endif ! nnee > 0
    ! Calculate log-likelihood for ecosystem respiration
    if (DATAin%nreco > 0) then
        mod = M_FLUXES(1:DATAin%nodays,3) &  ! Rauto
            + M_FLUXES(1:DATAin%nodays,13) & ! Rhet litter
            + M_FLUXES(1:DATAin%nodays,30) & ! Rhet wood litter
            + M_FLUXES(1:DATAin%nodays,14)   ! Rhet som
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nreco,DATAin%recopts,DATAin%Reco,DATAin%Reco_unc,DATAin%Reco_lag, &
                                             1d0,mod)
    endif ! nreco > 0
    ! Calculate log-likelihood for total wood net increment
    if (DATAin%nCwood_inc > 0) then
        mod = M_FLUXES(1:DATAin%nodays,7) - M_FLUXES(1:DATAin%nodays,11)
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCwood_inc,DATAin%Cwood_incpts, &
                                             DATAin%Cwood_inc,DATAin%Cwood_inc_unc,DATAin%Cwood_inc_lag, &
                                             1d0,mod)
    endif ! nCwood_inc > 0
    ! Calculate log-likelihood for total wood gross increment
    if (DATAin%nCwood_growth > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCwood_growth,DATAin%Cwood_growthpts, &
                                             DATAin%Cwood_growth,DATAin%Cwood_growth_unc,DATAin%Cwood_growth_lag, &
                                             1d0,M_FLUXES(1:DATAin%nodays,7))
    endif ! nCwood_inc > 0
    ! Calculate log-likelihood for total wood mortality
    if (DATAin%nCwood_mortality > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCwood_mortality,DATAin%Cwood_mortalitypts, &
                                             DATAin%Cwood_mortality,DATAin%Cwood_mortality_unc,DATAin%Cwood_mortality_lag, &
                                             1d0,M_FLUXES(1:DATAin%nodays,11))
    endif ! nCwood_mortality > 0

    return

  end subroutine calc_obs_likelihoods
  !
  !------------------------------------------------------------------
  !
  subroutine calc_scaled_obs_likelihoods(ML_obs_out, M_POOLS, M_FLUXES, M_DIAGS)
    use cardamom_structures, only: DATAin
    use carbon_model_memory, only: sw_par_fraction, top_soil_depth

    ! Subroutine to control the calculation of the observation related
    ! log-likelihoods and accounting for the various scalings

    ! Arguements
    double precision, intent(inout) :: ML_obs_out
    double precision, dimension(datain%nodays, datain%nofluxes)::  M_FLUXES
    double precision, dimension((DATAin%nodays+1), DATAin%nopools):: M_POOLS
    double precision, dimension(datain%nodays, datain%nodiags)::  M_DIAGS
    ! local variable
    double precision, dimension(DATAin%nodays) :: mod

    !
    ! Do diagnostics (DIAGS)
    ! 

    ! Calculate log-likelihood for fraction of absorbed radiation
    if (DATAin%nfAPAR > 0) then
        mod = M_DIAGS(1:DATAin%nodays,3) & ! APAR
            / (DATAin%met(4,1:DATAin%nodays)*sw_par_fraction)
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nfAPAR,DATAin%fAPARpts,DATAin%fAPAR,DATAin%fAPAR_unc,DATAin%fAPAR_lag, &
                                             DATAin%fAPAR_scaling,mod)
    endif ! nfAPAR > 0
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
    ! Calculate log-likelihood for total wood stocks
    if (DATAin%nCwood_stock > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCwood_stock,DATAin%Cwood_stockpts, &
                                             DATAin%Cwood_stock,DATAin%Cwood_stock_unc,DATAin%Cwood_stock_lag, &
                                             DATAin%Cwood_stock_scaling,M_POOLS(1:DATAin%nodays,4))
    endif ! nCwood_stock > 0
    ! Calculate log-likelihood for foliage litter stocks
    if (DATAin%nClit_stock > 0) then
        ! Estimate the foliage litter pool based on the ratio of foliage litter input to foliage + fine root litter inputs,
        ! scaled by the total litter pool. This is based on the turnover being common.
        mod = (sum(M_FLUXES(1:DATAin%nodays,10))/sum(M_FLUXES(1:DATAin%nodays,10)+M_FLUXES(1:DATAin%nodays,12))) &
                 * M_POOLS(1:DATAin%nodays,5)
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nClit_stock,DATAin%Clit_stockpts, &
                                             DATAin%Clit_stock,DATAin%Clit_stock_unc,DATAin%Clit_stock_lag, &
                                             DATAin%Clit_stock_scaling,mod)
    endif ! nClit_stock > 0
    ! Calculate log-likelihood for soil organic matter stocks
    if (DATAin%nCsom_stock > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCsom_stock,DATAin%Csom_stockpts, &
                                             DATAin%Csom_stock,DATAin%Csom_stock_unc,DATAin%Csom_stock_lag, &
                                             DATAin%Csom_stock_scaling,M_POOLS(1:DATAin%nodays,6))
    endif ! nCsom_stock > 0
!    ! Calculate log-likelihood for surface soil water
!    if (DATAin%nsoilwater > 0) then
!        mod = (M_POOLS(1:DATAin%nodays,8) * 1d-3) / top_soil_depth ! convert mm -> m3/m3
!        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nsoilwater,DATAin%soilwaterpts, &
!                                             DATAin%soilwater,DATAin%soilwater_unc,DATAin%soilwater_lag, &
!                                             DATAin%soilwater_scaling,mod)
!    endif ! nsoilwater > 0


    !
    ! Do fluxes (FLUXES)
    !

    ! Calculate log-likelihood for evapotranspiration
    if (DATAin%nEvap > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nEvap,DATAin%Evappts,DATAin%Evap,DATAin%Evap_unc,DATAin%Evap_lag, &
                                             DATAin%Evap_scaling,M_FLUXES(1:DATAin%nodays,29))
    endif ! nEvap > 0
    ! Calculate log-likelihood for fire
    if (DATAin%nFire > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nFire,DATAin%Firepts,DATAin%Fire,DATAin%Fire_unc,DATAin%Fire_lag, &
                                             DATAin%Fire_scaling,M_FLUXES(1:DATAin%nodays,17))
    endif ! nFire > 0
    ! Calculate log-likelihood for gross primary production
    if (DATAin%ngpp > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%ngpp,DATAin%gpppts,DATAin%GPP,DATAin%GPP_unc,DATAin%GPP_lag, &
                                             DATAin%GPP_scaling,M_FLUXES(1:DATAin%nodays,1))
    endif ! ngpp > 0
    ! Calculate log-likelihood for harvest
    if (DATAin%nharvest > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nharvest,DATAin%harvestpts,DATAin%harvest,DATAin%harvest_unc,DATAin%harvest_lag, &
                                             DATAin%harvest_scaling,M_FLUXES(1:DATAin%nodays,32))
    endif ! nharvest > 0
    ! Calculate log-likelihood for net biome exchange 
    if (DATAin%nnbe > 0) then
        mod = M_FLUXES(1:DATAin%nodays,3) &  ! Rauto
            + M_FLUXES(1:DATAin%nodays,13) & ! Rhet litter
            + M_FLUXES(1:DATAin%nodays,30) & ! Rhet wood litter
            + M_FLUXES(1:DATAin%nodays,14) & ! Rhet som
            + M_FLUXES(1:DATAin%nodays,17) & ! Fire
            - M_FLUXES(1:DATAin%nodays,1)    ! GPP
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nnbe,DATAin%nbepts,DATAin%NBE,DATAin%NBE_unc,DATAin%NBE_lag, &
                                             DATAin%NBE_scaling,mod)
    endif ! nnbe > 0
    ! Calculate log-likelihood for net ecosystem exchange of CO2
    if (DATAin%nnee > 0) then
        mod = M_FLUXES(1:DATAin%nodays,3) &  ! Rauto
            + M_FLUXES(1:DATAin%nodays,13) & ! Rhet litter
            + M_FLUXES(1:DATAin%nodays,30) & ! Rhet wood litter
            + M_FLUXES(1:DATAin%nodays,14) & ! Rhet som
            - M_FLUXES(1:DATAin%nodays,1)    ! GPP
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nnee,DATAin%neepts,DATAin%NEE,DATAin%NEE_unc,DATAin%NEE_lag, &
                                             DATAin%NEE_scaling,mod)    
    endif ! nnee > 0
    ! Calculate log-likelihood for ecosystem respiration
    if (DATAin%nreco > 0) then
        mod = M_FLUXES(1:DATAin%nodays,3)  & ! Rauto
            + M_FLUXES(1:DATAin%nodays,13) & ! Rhet litter
            + M_FLUXES(1:DATAin%nodays,30) & ! Rhet wood litter
            + M_FLUXES(1:DATAin%nodays,14)   ! Rhet som
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nreco,DATAin%recopts,DATAin%Reco,DATAin%Reco_unc,DATAin%Reco_lag, &
                                             DATAin%Reco_scaling,mod)
    endif ! nreco > 0
    ! Calculate log-likelihood for total wood net increment
    if (DATAin%nCwood_inc > 0) then
        mod = M_FLUXES(1:DATAin%nodays,7) - M_FLUXES(1:DATAin%nodays,11)
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCwood_inc,DATAin%Cwood_incpts, &
                                             DATAin%Cwood_inc,DATAin%Cwood_inc_unc,DATAin%Cwood_inc_lag, &
                                             DATAin%Cwood_inc_scaling,mod)
    endif ! nCwood_inc > 0
    ! Calculate log-likelihood for total wood gross increment
    if (DATAin%nCwood_growth > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCwood_growth,DATAin%Cwood_growthpts, &
                                             DATAin%Cwood_growth,DATAin%Cwood_growth_unc,DATAin%Cwood_growth_lag, &
                                             DATAin%Cwood_growth_scaling,M_FLUXES(1:DATAin%nodays,7))
    endif ! nCwood_inc > 0
    ! Calculate log-likelihood for total wood mortality
    if (DATAin%nCwood_mortality > 0) then
        ML_obs_out = ML_obs_out + likelihood(DATAin%nodays,DATAin%nCwood_mortality,DATAin%Cwood_mortalitypts, &
                                             DATAin%Cwood_mortality,DATAin%Cwood_mortality_unc,DATAin%Cwood_mortality_lag, &
                                             DATAin%Cwood_mortality_scaling,M_FLUXES(1:DATAin%nodays,11))
    endif ! nCwood_mortality > 0

    return

  end subroutine calc_scaled_obs_likelihoods
  !
  !------------------------------------------------------------------
  !
  subroutine calc_other_likelihoods(ML_obs_out, M_POOLS, M_FLUXES, M_DIAGS)
    use cardamom_structures, only: DATAin
    use carbon_model_memory, only: top_soil_depth

    ! Subroutine to control the calculation of the 'other priors'
    ! log-likelihoods. These are typically derived variables / 
    ! emergent functions averaged over the life of the anaylsis

    ! Arguements
    double precision, intent(inout) :: ML_obs_out
    double precision, dimension(datain%nodays, datain%nofluxes)::  M_FLUXES
    double precision, dimension((DATAin%nodays+1), DATAin%nopools):: M_POOLS
    double precision, dimension(datain%nodays, datain%nodiags)::  M_DIAGS
    ! local variable
    integer :: dummy_nodays = 1, dummy_noobs = 1
    integer, dimension(1) :: dummy_pts = 1, dummy_lag = 0
    double precision, dimension(1) :: mod
    double precision :: dummy_scaling = 1d0

    ! Initial soil water condition
!    if (DATAin%otherpriors(1) > -9998) then
!        mod = (M_POOLS(1,8) * 1d-3) / top_soil_depth ! convert mm -> m3/m3
!        ML_obs_out = ML_obs_out + (DATAin%otherpriorweight(1)*likelihood(dummy_nodays,dummy_noobs,dummy_pts, &
!                                   DATAin%otherpriors(1),DATAin%otherpriorunc(1),dummy_lag,dummy_scaling,mod))
!    end if
    ! Ra:GPP fraction is in this model a derived property
    if (DATAin%otherpriors(2) > -9998) then
        mod = sum(M_FLUXES(1:DATAin%nodays,3)) / sum(M_FLUXES(1:DATAin%nodays,1)) ! sum(Rauto) / sum(GPP)
        ML_obs_out = ML_obs_out + (DATAin%otherpriorweight(2)*likelihood(dummy_nodays,dummy_noobs,dummy_pts, &
                                   DATAin%otherpriors(2),DATAin%otherpriorunc(2),dummy_lag,dummy_scaling,mod))
    end if

    ! ...OTHERPRIOR(3)...

    ! Evaportranspiration (kgH2O/m2/day) as ratio of precipitation (kg/m2/s ->
    ! kg/m2/day)
    if (DATAin%otherpriors(4) > -9998) then
        mod = sum(M_FLUXES(1:DATAin%nodays,29)) / sum(DATAIN%MET(7,1:DATAin%nodays) * 86400d0)
        ML_obs_out = ML_obs_out + (DATAin%otherpriorweight(4)*likelihood(dummy_nodays,dummy_noobs,dummy_pts, &
                                   DATAin%otherpriors(4),DATAin%otherpriorunc(4),dummy_lag,dummy_scaling,mod))
    end if

    ! Estimate the biological steady state attractor on the wood pool.
    ! NOTE: this arrangement explicitly neglects the impact of disturbance on
    ! residence time (i.e. no fire and biomass removal)
    if (DATAin%otherpriors(5) > -9998) then
        mod = (sum(M_FLUXES(1:DATAin%nodays,7))/dble(DATAin%nodays)) & 
            * ( (sum(M_POOLS(1:DATAin%nodays,4) / (M_FLUXES(1:DATAin%nodays,11)))) / dble(DATAin%nodays))
        ML_obs_out = ML_obs_out + (DATAin%otherpriorweight(5)*likelihood(dummy_nodays,dummy_noobs,dummy_pts, &
                                   DATAin%otherpriors(5),DATAin%otherpriorunc(5),dummy_lag,dummy_scaling,mod))
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
!
!------------------------------------------------------------------
!
end module model_likelihood_module
