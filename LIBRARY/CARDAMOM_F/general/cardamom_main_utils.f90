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

!!!!!!!!!! File specific description !!!!!!!!!!
! Utility routines shared by the CARDAMOM main programs (EDC initial-value search, statistics initialisation/reset).
!
! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory);
! translation to Fortran, integration and subsequent modifications by T. L. Smallman,
! J. F. Exbrayat and colleagues (University of Edinburgh). Sampler-layer development
! (DEMCz, parallel samplers, wrappers) by J. Klebes (University of Edinburgh), 2024-2025.
! See function/subroutine specific comments for exceptions and contributors.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module cardamom_main_utils

  implicit none(type, external)

  ! Numeric code for available samplers
  integer, parameter :: sampler_APMCMC = 1, sampler_DEMCZ = 2

  public

  contains
  !
  !--------------------------------------------------------------------
  !
  pure function to_lower(str) result(str_lower)
    
    ! Description needed

    ! Argument
    character(len=*), intent(in) :: str
    ! Output
    character(len=len_trim(str)) :: str_lower
    ! Local variable
    integer :: i, j

    str_lower = str
    ! Loop characters
    do i=1, LEN_TRIM(str)
       j = iachar(str(i:i))
       if (j>=65 .and. j<97) then
           j = j+32
       endif
       str_lower(i:i) = achar(j)
    end do

  end function to_lower
  !
  !--------------------------------------------------------------------
  !
  pure logical function str_equal(str1, str2)

    ! check string match, case insensitive

    ! Arguments
    character(len=*), intent(in)  :: str1, str2

    ! Local variables
    character(len=LEN_TRIM(str1)) :: str1_lower
    character(len=LEN_TRIM(str2)) :: str2_lower

    ! Make lower case?
    str1_lower = to_lower(trim(str1))
    str2_lower = to_lower(trim(str2))

    ! Are strongs same?
    str_equal = str1_lower == str2_lower

    Return

  end function str_equal
  !
  !--------------------------------------------------------------------
  !
  subroutine parse_sampler_choice( arg1, arg2, arg3, sampler, args_start)

     ! Arguments
     character(len=*), intent(in):: arg1, arg2, arg3
     integer, intent(out):: sampler
     integer, intent(out):: args_start

     ! Local variables
     character(len=350):: arg, arg_next
     integer :: pos_in_word
     logical :: found_equals
     character(350) :: sampler_name

     args_start = 0
     sampler = sampler_APMCMC ! default when the optional command line argument is not present
     arg = trim(arg1)
     arg_next = arg2
     if (len_trim(arg) >= 7) then
         if (str_equal("sampler",arg1(1:7))) then
             args_start = 1
             pos_in_word = 8

             ! does it continue with = ?
             found_equals = .false.
             if (len_trim(arg) >= pos_in_word) then
                 if (str_equal("=",arg(pos_in_word:pos_in_word))) then
                     found_equals = .true.
                     pos_in_word = pos_in_word + 1
                 endif
             endif

            ! or is = in the next argument ?
            if (.not. found_equals) then
                args_start = 2
                arg = arg2
                arg_next = arg3
                pos_in_word = 1
                if (len_trim(arg2) >= pos_in_word) then
                    if (str_equal("=",arg2(pos_in_word:pos_in_word))) then
                        found_equals = .true.
                        pos_in_word = pos_in_word + 1
                    endif
                endif
            endif

            ! else parse error, `sampler` but no `=`
            if (.not. found_equals) then
                write(*,*) "Could not parse command line, found keyword `sampler` but no `=`"
                STOP 1
            endif

            ! does it continue with a sampler name in the same word?
            if (len_trim(arg) >= pos_in_word + 4) then
                sampler_name = arg(pos_in_word:pos_in_word+4) !really only checking first 5 chars
            else
                args_start = args_start + 1
                pos_in_word = 1
                arg = trim(arg_next) !arg2 or arg3
                sampler_name = arg(pos_in_word:pos_in_word+4) !really only checking first 5 chars
            endif
            if (str_equal("demcz", sampler_name) ) then
                sampler = sampler_DEMCZ
            elseif (str_equal("APMCM", sampler_name) ) then
                sampler = sampler_APMCMC
            else
                write(*,*) "Could not parse command line, found keyword `sampler=` ", sampler_name
                write(*,*) "Choose from sampler=DEMCZ, APMCMC"
                STOP 1
            endif
         endif
     endif
     ! else : `sampler` keyword not found, continue with default

  end subroutine parse_sampler_choice  
  !
  !--------------------------------------------------------------------
  !
  subroutine initialize_stats(MCOUT, npars)
    use apmcmc, only: MCMC_OUTPUT

    integer, intent(in) :: npars
    type(MCMC_OUTPUT), intent(inout) :: MCOUT
    allocate (MCOUT%covariance(npars, npars), MCOUT%parvar(npars), MCOUT%meanpar(npars))

    call reset_stats(MCOUT, npars)

  end subroutine initialize_stats
  !
  !--------------------------------------------------------------------
  !
  subroutine reset_stats(MCOUT, npars)
    use apmcmc, only: MCMC_OUTPUT

    type(MCMC_OUTPUT), intent(inout) :: MCOUT
    integer, intent(in) :: npars
    integer :: n

    MCOUT%parvar = 1d0; MCOUT%Nparvar = 0d0
    ! Covariance matrix cannot be set to zero therefore set initial
    ! value to a small positive value along to variance access
    MCOUT%covariance = 0d0; MCOUT%meanpar = 0d0; MCOUT%cov = .false.
    MCOUT%use_multivariate = .false.
    do n = 1, npars
       MCOUT%covariance(n, n) = 1d0
    end do

  end subroutine reset_stats
  !
  !------------------------------------------------------------------
  !
  subroutine find_edc_initial_values(MCOUT_list, nchains, seed)
    !! subroutine deals with the determination of initial parameter and initial
    !! conditions which are consistent with EDCs
    !! pre-loop, Run MCMC sampler with modified likelihood fct
    use model_shared, only: PI
    use apmcmc, only: MCMC_OUTPUT, run_mcmc, run_parallel_mcmc, MCMC_options
    use model_likelihood_wrapper  ! TODO next refactoring step
    use cardamom_structures, only: DATAin

    implicit none(type, external)

    ! declare local variables
    integer, intent(in) :: nchains, seed
    type(MCMC_OUTPUT), dimension(:), allocatable, intent(inout) :: MCOUT_list
    type(MCMC_OUTPUT), dimension(:), allocatable :: MCOUT_list_tmp
    type(mcmc_OPTIONS) :: MCO ! options settings will be on defaults
    integer :: i, counter_local(nchains), nOUT_save, nWRITE_save, nADAPT_save
    integer :: success_count
    logical :: append_save
    double precision :: ll
    double precision :: PEDC(nchains), PEDC_prev(nchains), P_target
    double precision, dimension(PI%npars) :: parini  ! local variable, or array

    allocate (MCOUT_list_tmp(nchains))

    ! Hold for later
    nOUT_save = MCO%nOUT; nWRITE_save = MCO%nWRITE; nADAPT_save = MCO%nADAPT
    append_save = MCO%append

    ! set MCMC options needed for EDC run
    ! same for all chains
    MCO%APPEND = .false.
    MCO%nADAPT = 500
    MCO%fADAPT = 1d0
    MCO%nOUT = 100000
    MCO%nPRINT = 0
    MCO%nWRITE = 0
    ! the next two lines ensure that parameter inputs are either given or
    ! entered as-9999
    MCO%randparini = .true.
    MCO%returnpars = .true.
    MCO%fixedpars = .false.

    ! Set initial priors to vector...
    ! TODO array for multichain
    parini = DATAin%parpriors(1:PI%npars)
    ! ... and assume we need to find random parameters
    ! Target likelihood allows for controlling when the MCMC will stop
    P_target = 0d0

    ! if the prior is not missing and we have not told the edc to be random
    ! keep the value
!    do n = 1, PI%npars
!       if (PI%parini(n) /= -9999d0 .and. DATAin%edc_random_search < 1) PI%parfix(n) = 1
!    end do  ! parameter loop

     do i = 1, nchains
        MCOUT_list_tmp(i)%PARS = parini
        call initialize_stats(MCOUT_list_tmp(i), PI%npars)
     end do

     ! if this is not a restart run, i.e. we do not already have a starting
     ! position we must being the EDC search procedure to find an ecologically
     ! consistent initial parameter set

     ! set up edc log likelihood for MHMCMC initial run
     PEDC_prev = -1000d0; PEDC = -1d0; counter_local = 0
     success_count = 0

     write (*,*) "Beginning EDC search attempt "
     write (*,*) nchains, "chains working ... "
     MCO%randparini = .false.
     ! TODO limit number to number of available hardware threads
     !$omp parallel do private(ll) firstprivate(MCO)
     do i = 1, nchains
        do while (success_count < nchains)!(PEDC < 0d0)

           ! call the MHMCMC directing to the appropriate likelihood function
           call run_mcmc(edc_model_likelihood_fct, PI, MCO, MCOUT_list_tmp(i), model_likelihood_fct, chainid=i, seed = seed+i)
           MCO%fixedpars = .true. !continue from this position next round, until resetting every 5 attempts

           ! turn off random selection for initial values

           ! store the best parameters from that loop
           PEDC(i) = MCOUT_list_tmp(i)%bestll
           ! if any chains's MCOUT object is success (reached loglikelihood = 0),
           ! copy it to MCOUT_list
           !$omp critical
           if (MCOUT_list_tmp(i)%ll >= (0d0 - 10*epsilon(1.0d0))) then
              success_count = success_count + 1
              write (*, *) "Found ", success_count, "EDC-compatible starting points"
              if (success_count <= nchains) then
                 MCOUT_list(success_count) = MCOUT_list_tmp(i)
              end if
              call reset_stats(MCOUT_list_tmp(i), PI%npars)
              MCO%randparini = .true. ! TODO problem
              PEDC_prev(i) = -1000d0
              MCOUT_list_tmp(i)%PARS = DATAin%parpriors(1:PI%npars)
              counter_local(i) = 0
              MCO%fixedpars = .false.

           end if
           !$omp end critical

           ! keep track of attempts
           counter_local(i) = counter_local(i) + 1
           ! periodically reset the initial conditions
           if (PEDC(i) < 0d0 .and. PEDC(i) <= PEDC_prev(i) .and. counter_local(i) > 5) then
               ! Reset the previous EDC likelihood score
               PEDC_prev(i) = -1000d0
               ! Reset parameters back to default
               MCOUT_list_tmp(i)%PARS = DATAin%parpriors(1:PI%npars)
               ! reset to select random starting point
               MCO%randparini = .true. ! TODO problem
               ! reset the parameter step size at the beginning of each attempt
               call reset_stats(MCOUT_list_tmp(i), PI%npars)
               counter_local(i) = 0
               MCO%fixedpars = .false.
           else
               PEDC_prev(i) = PEDC(i)
           end if
        end do  ! for while condition
     end do
     !$omp end parallel do

     ! reset so that currently saved parameters will be used
     ! starting point in main MCMC
     ! PI%parfix(1:PI%npars) = 0  ! TODO
     !MCOUT%bestpars = 0d0

  end subroutine find_edc_initial_values
  !
  !--------------------------------------------------------------------
  !
end module cardamom_main_utils
