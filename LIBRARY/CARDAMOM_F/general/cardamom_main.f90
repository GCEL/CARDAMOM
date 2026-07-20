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
! along with this program.  If not, see < https://www.gnu.org/licenses/>.

!!!!!!!!!!!! File specific description !!!!!!!!!!
! This is the main subroutine for the CARDAMOM framework. The specific model
! method combinations are achieved through case specific compilation of the
! case while maintaining strict consistent io formats to allow for these
! combinations
!
! This code is based on the original C verion of the University of Edinburgh
! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
! All code translation into Fortran, integration into the University of
! Edinburgh CARDAMOM code and subsequent modifications by:
! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
! J. F. Exbrayat (University of Edinburgh)
! D. T. Milodowski (d.t.milodowski@ed.ac.uk, University of Edinburgh)
! See function/subroutine specific comments for exceptions and contributors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program cardamom_framework
   use cardamom_MHMCMC, only: MCMC_OUTPUT, run_mcmc, run_parallel_mcmc
   use DEMCZ, only: run_demcz, DEMCzopt ! problem different options object
   use model_shared, only: PI, initialize_carbon_model
   use samplers_shared, only: init_infinity, & ! should be in Utils ?
     & SAMPLER_OPTIONS, MCMC_OPTIONS
   use cardamom_structures, only: DATAin
   use cardamom_io, only: initialize, &
                          read_options, &
                          update_obs_scaling_normal, update_obs_scaling_nsamples, &
                          update_obs_scaling_sqrt_nsamples, update_obs_scaling_log_nsamples
   use samplers_io, only: open_output_files, &
                          check_for_existing_output_files, &
                          update_for_restart_simulation, &
                          write_covariance_matrix, &
                          close_output_files, write_covariance_info
   !use MHMCMC_StressTests, only: StressTest_likelihood_fct, StressTest_sublikelihood_fct, prepare_for_stress_test
   use model_likelihood_wrapper, only: model_likelihood_fct, edc_model_likelihood_fct, scaled_model_likelihood_fct
   use cardamom_main_utils

 !!!!!!!!!!!
   ! Authorship contributions
   !
   ! This code is based on the original C verion of the University of Edinburgh
   ! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
   ! All code translation into Fortran, integration into the University of
   ! Edinburgh CARDAMOM code and subsequent modifications by:
   ! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
   ! J. F. Exbrayat (University of Edinburgh)
   ! D. T. Milodowski (d.t.milodowski@ed.ac.uk, University of Edinburgh)
   ! See function/subroutine specific comments for exceptions and contributors
 !!!!!!!!!!!

   ! Created: Anthony A. Bloom
   ! Major modification history:
   ! Version 1.0: C language MHMCMC and io created by Anthony A. Bloom
   ! Version 1.1: Translated into Fortran, including some generic function replacement by T. L. Smallman
   ! Version 1.2: Analysis restart capacity added by T.L. Smallman
   ! Version 1.3: MHMCMC updated to APMCMC by T.L. Smallman with advice from Anthony A. Bloom
   !            : Following Haario et al., (2001, 2006) and Roberts & Rosenthal (2009).
   ! Version 1.4: Pre-APMCMC phase using normalised likelihoods added by T. L. Smallman
   !            : Pre-APMCMC allows for rapidly moving towards observations from very bad starting points.
   ! Version 1.5: Added options for varied scaling approaches for the cost function (e.g. division by sample size (i.e. n), sqrt(n), 1+log(n))
   ! Specific citations for developments included in the code.

   ! This is the main subroutine for the CARDAMOM framework. The specific model
   ! method combinations are achieved through case specific compilation of the
   ! case while maintaining strict consistent io formats to allow for these
   ! combinations

   ! Command line inputs are:
   ! optional "sampler=APMCMC" , "DEMcz", "MHMCMC"
   ! 1) file in
   ! 2) file out
   ! 3) integer number of solutions requested
   ! 4) print-to-screen frequency
   ! 5) write-to-file frequency
   ! 6) 0/1 flag to use normalised log-likelihood pre-mcmc
   ! 7) Flag to select the cost function normalisation approach
   ! 8) integer number of chains (optional; defaults to 3 if absent)

   implicit none(type, external)

   ! declare local variables
   character(350) :: infile, outfile, solution_wanted_char, freq_print_char, &
                    freq_write_char, do_inflate_char, cost_func_scaling_char, &
                    nchains_char
   character(350) :: arg1, arg2, arg3
   integer :: args_start 
   integer :: sampler ! set to one of the integer codes 
   integer :: solution_wanted, freq_print, freq_write, time1, time2, time3, &
              do_inflate_dble, cost_func_scaling_dble, idum
   logical :: do_inflate = .false.
   logical :: sub_sample_complete = .false.
   double precision :: nOUT_save = 0d0, sub_fraction = 0.2d0
    !! run this percentage of simulation with variant function
   type(MCMC_OUTPUT), dimension(:), allocatable :: MCOUT_list
    !! array of output objects from each thread
   class(SAMPLER_OPTIONS), allocatable :: MCO
     !! options for sampler
   logical :: restart

   ! Number of chains. Default value, may be overridden by command line argument 8.
   integer :: nchains = 3
   integer :: i

   call init_infinity()

   ! check command line for a potential first arg sampler=... 
   ! effects : set this program's `sampler` variable  and args_start to offset following command line args
   call get_command_argument(1, arg1)
   call get_command_argument(2, arg2)
   call get_command_argument(3, arg3)
   call parse_sampler_choice(arg1, arg2, arg3, sampler, args_start)

   ! user update
   write (*,*) "Beginning read of the command line"
   ! read user options from the command line
   call get_command_argument(args_start+1, infile)
   call get_command_argument(args_start+2, outfile)
   call get_command_argument(args_start+3, solution_wanted_char)
   call get_command_argument(args_start+4, freq_print_char)
   call get_command_argument(args_start+5, freq_write_char)
   call get_command_argument(args_start+6, do_inflate_char)
   call get_command_argument(args_start+7, cost_func_scaling_char)
   ! argument 8 (number of chains) is optional; default retained if absent
   if (command_argument_count() >= args_start+8) then
      call get_command_argument(args_start+8, nchains_char)
      read (nchains_char, '(I10)') nchains
   end if

   ! now convert relevant ones to integrator
   ! Note: that I10 is the maximum allowed with the default integer kind (kind = 4).
   !       allow for larger number of iterations integer (kind  = 8) is needed throughout the code.
   read (solution_wanted_char, '(I10)') solution_wanted
   read (freq_print_char, '(I10)') freq_print
   read (freq_write_char, '(I10)') freq_write
   read (do_inflate_char, '(I10)') do_inflate_dble
   read (cost_func_scaling_char, '(I10)') cost_func_scaling_dble

   ! Assign inflate logical condition
   if (do_inflate_dble == 1) do_inflate = .true.
   ! Sanity check of the cost_function_scaling_char
   if (cost_func_scaling_dble > -1 .and. cost_func_scaling_dble < 4 .and. freq_write > 0) then
      ! All is well
   else
      ! All is not well-complain
      print*, "ERROR: Command line argument to specify the cost function or write to file frequency is incorrect."
      print*, "Command line should have 7 arguments (in addition to the cardamom.exe)."
      print*, "These are: "
      print*, "(optional) sampler=<SAMPLER> choose from APMCMC, DEMcz, MHMCMC "
      print*, "1) input file path."
      print*, "2) output file path-note that PARS, STEP, COV, COVINFO will be appended to this name path outfile."
      print*, "3) No. of parameter proposals to make."
      print*, "4) Iteration freq. for printing to screen (main MCMC phase only)."
      print*, "5) Iteration freq. for writing results to files."
      print*, "6) do sample size normalisation phase - "
      print*, "  0 = FALSE"
      print*, "  1 = TRUE"
      print*, "7) Select likelihood cost function - "
      print*, "  0 = no scaling"
      print*, "  1 = scaling by sample size (n)"
      print*, "  2 = scaling by sqrt(n)"
      print*, "  3 = scaling by log(n)"
      print*, "8) Number of chains (optional, integer >= 1; defaults to 3)."
      stop
   end if
   
   ! decide on the exact type of the sampler_options struct
   if (sampler==sampler_APMCMC .or. sampler==sampler_MHMCMC) then 
     allocate(MCMC_OPTIONS::MCO)
   else if (sampler==sampler_DEMCZ) then
     allocate(DEMCzopt::MCO)
   endif

   ! Sanity check the number of chains.
   if (sampler==sampler_DEMCZ .and. nchains < 3) then
      print *, "ERROR: number of chains (command line argument 8) must be >= 3 for DEMCz."
      print *, "Value supplied = ", nchains
      stop
    else if (nchains < 1) then
      print*, "ERROR: number of chains (command line argument 8) must be >= 1."
      print*, "Value supplied = ", nchains
      stop
   end if
   MCO%nchains = nchains !note it in the settings object

   if (sampler==sampler_DEMCZ .and. do_inflate) then
    write(*,*) "WARNING: Value of 6th command line argument do_inflate=1 will be ignored, there will be no &
                    &  log-scaled sampling phase when the chosen sampler is DEMCZ" 
   endif 

   ! user update
   write (*,*) "Command line options read, moving on now"

   ! Now the number of chains is known, allocate the per-chain output array
   allocate (MCOUT_list(nchains))

   ! determine unique (sort of) seed value; based on system time
   call system_clock(time1, time2, time3)
   ! Set seed value outside of the function
   idum = time1+time2+time3
   !call rnstrt(nint(idum))

   ! Determine whether ot not we are doing a real analysis or running a stress trest
   if (trim(infile) == "StressTest") then
      !call run_stresstest()
      ! call prepare_for_stress_test(infile, outfile)  ! sets cardamom_structures :: DATAin
      write(*,*) "StressTest is now located in test/ and called via `ctest`, not by running the main &
        & cardamom program with the StressTest keyword"
      stop
   end if

   ! Report which model ID we are using
   write (*,*) "Running model version ", DATAin%ID

   ! read input data file
   call initialize(infile) ! = initialize_parinfo, read_check_binary_data, initialize_model
   ! sets cardamom_structures :: DATAin

   call initialize_carbon_model(nchains)
   do i = 1, nchains
      call initialize_stats(MCOUT_list(i), PI%npars)
   end do

   ! having filled PI%npars from model file, we can allocate stats array in MCOUT

   ! load module variables needed for restart check
   ! NOTE: THIS MUST HAPPEN BEFORE CHECKING FOR RESTART
   call read_options(solution_wanted, freq_print, freq_write, outfile, MCO)

   ! TODO check for restart could be done on each file independently
   ! check whether this is a restart from aborted simulation
   MCO%restart = .true. !to gather results onto-false as soon as some file is not found
   do i = 1, nchains
      call check_for_existing_output_files(PI%npars, MCO, sub_fraction, i, restart)
      MCO%restart = MCO%restart .and. restart
   end do
    
      !if all nchains files were found, read them
    if (MCO%restart) then
       do i = 1, nchains
         call update_for_restart_simulation(MCO, MCOUT_list(i), PI%npars, i)
       end do
    end if

   if (.not. MCO%restart) then
      ! Begin search for initial conditions
      write (*,*) "Beginning search for initial parameter conditions"
      ! Determine initial values, this requires using the AP-MCMC
      call find_edc_initial_values(MCOUT_list, nchains, idum)
      ! Having done EDC search phase, flag to start the next phase from this state
      MCO%fixedpars = .true.
      do i = 1, nchains
         MCOUT_list(i)%nos_iterations = 0
         ! Reset stepsize and covariance for main MH-MCMC
         call reset_stats(MCOUT_list(i), PI%npars)
      end do
   end if

   ! At this point, whether from restart or from edc search , we should have values to start from 
   ! in MCOUT_list(i)%pars 
   ! In case of restart, number of steps already done is noted in MCOUT_list(i)%nos_iterations

   do i = 1, nchains
      ! Reset the MCMC parameters for the next stage
      call read_options(solution_wanted, freq_print, freq_write, outfile, MCO)
   end do

   ! sub-sampling phase, first sub_fraction% of the simulation with variant loglikelihood
   ! TODO using MCOUT_list(1)
   if (sampler /= sampler_DEMCZ .and. DATAin%total_obs > 0 .and. &
       & MCOUT_list(1)%nos_iterations < (MCO%nOUT*sub_fraction) .and. do_inflate) then

      ! Having found an EDC compliant parameter vector, we want to do a MCMC
      ! search on inflated uncertainties. This inflation search allows us to
      ! more easily move towards higher likelihoods, where the inflation
      ! allows easier movement through parameter space.

      ! The inflated search phase will make use of three stages over which the
      ! inflation will be reduced. Phase 1 used half of the allocated
      ! iterations for the inflation while the second and third is

      ! Set flag to indicate this phase has occurred and make a record of the
      ! total iterations to be attempted
      sub_sample_complete = .true.

      ! Report to the user
      write (*,*) "Beginning parameter search on sample size normalised likelihoods"

      ! Set MCMC parameters
      nOUT_save = nint(dble(MCO%nOUT)*sub_fraction) ; MCO%nOUT = nOUT_save !number of steps to do in this phase

      ! Actual number still to do = this minus MCOUT%nos_iterations steps already done pre-restart
      MCO%nOUT = MCO%nOUT - MCOUT_list(1)%nos_iterations ! TODO again just using 1 , assuming all restart files are same length

      MCO%fADAPT = 1d0
      MCO%fixedpars = .true. ! start from end points of EDC (or pre-restart) phase , preserve latest MCO%pars already existing in MCO 

      ! Update user again
      write (*,*) "Nos iterations to be proposed = ", MCO%nOUT
      ! Second phase, run MCMC with sub scaling
      call update_obs_scaling_nsamples
         select type(MCO)
         type is(MCMC_OPTIONS)
      call run_parallel_mcmc(scaled_model_likelihood_fct, PI, MCO, MCOUT_list, model_likelihood_fct, nchains=nchains, seed = idum)
    end select
      MCO%fixedpars = .true.
      do i = 1, nchains
         ! Use the best (instead of latest) parameter set as the starting point for the next stage
         MCOUT_list(i)%pars(1:PI%npars) = MCOUT_list(i)%bestpars(1:PI%npars)

         ! Leave parameter and covariance structures as they come out form the
         ! sub-sample-but reset the number of samples used in the update
         ! weighting
         if (sampler == sampler_APMCMC .or. sampler == sampler_MHMCMC) then
         select type(MCO)
         type is(MCMC_OPTIONS)
           if (MCOUT_list(i)%cov .and. MCOUT_list(i)%use_multivariate) then
              MCOUT_list(i)%Nparvar = MCO%N_before_mv*PI%npars + 1
           else
         endif
         end select
            ! reset the parameter step size at the beginning of each attempt
            call reset_stats(MCOUT_list(i), PI%npars)
         end if  ! do we need a new covariance matrix or can we use the existing one?

         ! At this point, should be nOUT_save == MCout_list(i)%nos_iterations
         ! (desired number of steps in this phase == number of steps actually done)
         ! TODO this phase should NOT be able to exit early due to convergence condition !
         ! tmp : check
         if (.not. (nOUT_save == MCOUT_list(i)%nos_iterations) ) then 
           write(*,*) "WARNING sampling first phase completed an unexpected number of steps, &
             & ", MCOUT_list(i)%nos_iterations, " on thread " , i, " where ", nOUT_save,  "were expected."
         endif 

         ! reset iterations counter.
         MCOUT_list(i)%nos_iterations = 0
      end do

   elseif (restart) then
     ! Assume the initial phase was done , and not recorded in the output files before restart
     ! Continue as if this number of steps has been done
     ! NOTE always assuming that the number of steps requested (command line argument) is the same
     ! on original run and restart run
     nOUT_save = nint(dble(MCO%nOUT)*sub_fraction) 
   else ! no first phase was done
     nOUT_save = 0 ! number of steps done in first phase
   end if

   ! Restore module variables needed for the run-these components could be split
   ! into two subroutines to avoid double calling of file name creation
   ! components.
   call read_options(solution_wanted, freq_print, freq_write, outfile, MCO)

   MCO%nOUT = MCO%nOUT - nOUT_save  ! number of steps needed in main phase = total number minus number done in first phase

   if ( MCO%restart ) then
     ! Restart : subtract what's already present in files
     MCO%nOUT = MCO%nOUT - MCOUT_list(1)%nos_iterations
   endif 

   MCO%fixedpars = .true. ! start next phase from values in MCOUT%pars, whether from restart, edc, or first phase

   ! Update the user
   write (*,*) "Beginning parameter search in real likelihoods"
   write (*,*) "Nos iterations to be proposed = ", MCO%nOUT 

   ! Call the main MCMC
   ! The specific normalisation of the cost function is determined here.
   ! But to avoid getting through the EDC do_inflate sections before finding
   ! out that the cost_function_scaling has not been set correctly,
   ! ensure code after the command line read (above) has been correctly maintained
   if (cost_func_scaling_dble == 0) then
      call update_obs_scaling_normal
   else if (cost_func_scaling_dble == 1) then
      call update_obs_scaling_nsamples
   else if (cost_func_scaling_dble == 2) then
      call update_obs_scaling_sqrt_nsamples
   else if (cost_func_scaling_dble == 3) then
      call update_obs_scaling_log_nsamples
   end if  ! cost_func_scaling_dble ==

   !  Finally run the mcmc
   if (sampler==sampler_APMCMC) then
   select type(MCO)
   type is (MCMC_options)
       call run_parallel_mcmc(scaled_model_likelihood_fct, PI, MCO, MCOUT_list, model_likelihood_fct, nchains=nchains, seed = idum)
   end select
   else if (sampler==sampler_DEMCZ) then
         select type(MCO)
         type is(DEMCzopt)
   call run_demcz(scaled_model_likelihood_fct, PI, MCO, MCOUT_list, model_likelihood_fct, nchains=nchains, seed = idum)
         end select
   endif

   ! Let the user know we are done
   write (*, *) "AP-MCMC done now, moving on ..."

   ! Final message to the user
   write (*, *) "==========================================================="
   write (*, *) "==== CARDAMOM analysis for the current site completed ====="
   write (*, *) "==========================================================="
   write (*, *) "=========================Honestly=========================="
   write (*, *) "==========================================================="
contains

end program cardamom_framework
