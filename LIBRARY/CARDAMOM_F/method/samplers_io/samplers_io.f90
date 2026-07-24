
module samplers_io

   !!!!!!!!!!!
   !
   ! Samplers_io functions split from cardamom_io.f90
   !
   ! Authorship contributions
   !
   ! This code is based on the original C verion of the University of Edinburgh
   ! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
   ! All code translation into Fortran, integration into the University of
   ! Edinburgh CARDAMOM code and subsequent modifications by:
   ! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
   ! J. F. Exbrayat (University of Edinburgh)
   ! See function/subroutine specific comments for exceptions and contributors
   !!!!!!!!!!!

   ! Module contains subroutines and variables needed to output parameter,
   ! likelihood and step size information from the MHMCMC.

   implicit none(type, external)

   ! declare private
   private

   ! allow access to specific functions
   public:: write_mcmc_output &
           ,write_parameters &
           ,write_variances &
           ,write_covariance_matrix &
           ,write_covariance_info &
           ,check_for_existing_output_files &
           ,update_for_restart_simulation &
           ,initialize_buffers &
           ,open_output_files &
           ,close_output_files

   integer:: pfile_unit = 10, sfile_unit = 11, cfile_unit = 12, cifile_unit = 13
   !! In case of single set of output files, these are literally the file unit numbers.
   !! In case of MCMC simulation, these are ids of the first thread's files; others are
   !! calculated based on them.

   ! parameters
   ! TODO Compiler dependent, should check kind of our doubles
   integer, parameter:: real_bytes = 8  ! number of bytes in real variable, 8 bytes is to make double precision

   type io_buffer_space
      integer:: io_buffer, io_buffer_count
      double precision, allocatable, dimension(:, :) :: &
         variance_buffer, &
         meanpars_buffer, &
         pars_buffer

      double precision, allocatable, dimension(:) :: &
         nsample_buffer, &
         accept_rate_buffer, &
         prob_buffer
   end type io_buffer_space
   ! allow access to needed variable
   public:: io_buffer_space

   save


contains
   !
   !------------------------------------------------------------------
   !

   subroutine calculate_file_ids(chainid, pfile_unit_thread, sfile_unit_thread, cfile_unit_thread, cifile_unit_thread)
      !! Calculate unit numbers for the simulation's nth thread's files.
      !! Applies to APMCMC and MHMCMC sampelrs, where each thread writes its own set of output files.
      !! The numbers are based on module variables pfile_unit etc , Assigned io unit numbers in the pattern
      !!, where module data is integer:: pfile_unit = 10, sfile_unit = 11, cfile_unit = 12, cifile_unit = 13
      !! 10  pfile thread 1
      !! 11  sfile thread 1
      !! 12  cfile thread 1
      !! 13  cifile thread 1
      !! 14  pfile thread 2
      !! 15  sfile thread 2
      !! 16  cfile thread 2
      !! 17  cifile thread 2
      !! ...
      integer, intent(in):: chainid
      integer  :: offset
      integer, intent(out) :: pfile_unit_thread, sfile_unit_thread, cfile_unit_thread, cifile_unit_thread
      offset = (chainid - 1)*4
      pfile_unit_thread = pfile_unit + offset
      sfile_unit_thread = sfile_unit + offset
      cfile_unit_thread = cfile_unit + offset
      cifile_unit_thread = cifile_unit + offset
   end subroutine calculate_file_ids

   subroutine check_for_existing_output_files(npars, MCO, sub_fraction, chainid, restart)
      use samplers_shared, only: SAMPLER_OPTIONS, filenames_insert_threadid

      ! subroutine checks whether both the parameter and step files exist for this
      ! job. If they do we will assume that this is a restart job that we want to
      ! finish off. Important for large jobs or running on machines with may crash
      ! / have runtime limits
      ! For a single thread's output
      implicit none(type, external)

      ! declare input variables
      integer, intent(in):: npars
      integer:: nOUT, nWRITE
      type(SAMPLER_OPTIONS), intent(in):: MCO
      !! simulation settings object-to read nOut, nWrite, filenames
      double precision, intent(in):: sub_fraction
      character(350):: outfile, stepfile, covfile, covifile
         !! filenames stems + numbering
      integer :: pfile_unit_thread, sfile_unit_thread, cfile_unit_thread, cifile_unit_thread
         !! file unit numbers for this thread
      integer, intent(in) :: chainid
         !! Thread id.  To get behavior of a single-threaded simulation use 1
      logical, intent(out):: restart


      ! local variables
      logical:: par_exists, step_exists, cov_exists, covinfo_exists
      double precision:: dummy
      integer:: num_lines, status

      nOUT = MCO%nOUT ; nWRITE = MCO%nWRITE

      ! process file names- reconstruct numbered file names as they would have been written by
      ! a simulation recieving the same MCO settings object
      outfile = MCO%outfile
      stepfile = MCO%stepfile
      covfile = MCO%covfile
      covifile = MCO%covifile
      if (MCO%nchains > 1)  then
          call filenames_insert_threadid(outfile, stepfile, covfile, covifile, chainid)
      end if

      ! Check that all files exist
      inquire (file=trim(outfile), exist=par_exists)
      inquire (file=trim(stepfile), exist=step_exists)
      inquire (file=trim(covfile), exist=cov_exists)
      inquire (file=trim(covifile), exist=covinfo_exists)

      ! now determine the correct response
      if (par_exists .and. step_exists .and. cov_exists .and. covinfo_exists) then

         ! All files exist therefore this might be a restart run.
         ! lets see if there is anything in the files that we might use
         ! count the number of remaining lines in the file..
         ! open the relevant output files
         call open_output_files(outfile, stepfile, covfile, covifile, chainid, pfile_unit_thread, sfile_unit_thread, cfile_unit_thread, cifile_unit_thread)

         status = 0; num_lines = 0
         do
            read (pfile_unit_thread, iostat=status) dummy
            if (status /= 0) exit
            num_lines = num_lines + 1
         end do
         ! Re-use dummy to calculate the target file size to be considered for
         ! restart
         dummy = ((dble(nOUT)/dble(nWRITE))*sub_fraction)*dble(npars + 1)
         if (num_lines > dummy) then
            ! Then there is something in the file we we can use it
            restart = .true.
            print *, "...have found parameter file = ", trim(outfile)
            print *, "...have found step file = ", trim(stepfile)
            print *, "...have found cov file = ", trim(covfile)
            print *, "...have found cov_info file = ", trim(covifile)
         else
            ! The file exists but is empty/no enough so treat it as a fresh start
            restart = .false.
            print *, "Output files are present, however they are too small for a restart"
         end if
         ! Either way we open the file up later on so now we need to close them
         call close_output_files(pfile_unit_thread, sfile_unit_thread, cfile_unit_thread, cifile_unit_thread)

      else  ! par_exists .and. step_exists

         ! Then or of these files exists and the other does not so it is
         ! ambiguous whether or not this is a restart
         print *, "One or more of the analysis files cannot be found."
         print *, "CARDAMOM must start from scratch... "
         restart = .false.

      end if  ! par_exists .and. step_exists

   end subroutine check_for_existing_output_files
   !
   !------------------------------------------------------------------
   !
   subroutine update_for_restart_simulation(MCO, MCOUT, npars, chainid)
      !! subroutine is responsible for loading previous parameter and step size
      !! information into the current object MCOUT.
      !! modifies: arg MCOUT%pars. To be used as starting point for next run.
      !! Also MCOUT%nos_iterations , %parvar, %covariance, %meanpar, %nparvar
      use samplers_shared, only: MCMC_OUTPUT, SAMPLER_OPTIONS, filenames_insert_threadid
      use samplers_math, only: std, covariance_matrix, inverse_matrix, par2nor

      implicit none(type, external)

      ! Arguments
      class(SAMPLER_OPTIONS), intent(inout):: MCO
      type(MCMC_OUTPUT), intent(inout):: MCOUT
      integer, intent(in):: npars
      integer, intent(in):: chainid
      character(350):: outfile, stepfile, covfile, covifile
         !! filenames

      ! local variables
      integer:: a, b, c, i, j, num_lines, status
      double precision:: dummy
      double precision, dimension(:, :), allocatable:: tmp
      integer  :: pfile_unit_thread, sfile_unit_thread, cfile_unit_thread, cifile_unit_thread
        !! file unit numbers for this thread

      ! process file names- reconstruct numbered file names as they would have been written by
      ! a simulation recieving the same MCO settings object
      outfile = MCO%outfile
      stepfile = MCO%stepfile
      covfile = MCO%covfile
      covifile = MCO%covifile
      if (MCO%nchains > 1)  then
          call filenames_insert_threadid(outfile, stepfile, covfile, covifile, chainid)
      end if

      ! open output files, determine file unit numbers
      ! last 4 arguments set with the unit numbers opened
      call open_output_files(outfile, stepfile, covfile, covifile, chainid, pfile_unit_thread, sfile_unit_thread, cfile_unit_thread, cifile_unit_thread)

      ! the parameter and step files should have already been openned so
      ! read the parameter and step files to get to the end

      ! rewind to the beginning
      rewind(pfile_unit_thread); rewind(sfile_unit_thread); rewind(cifile_unit_thread)

      !
      ! As this subroutine will only be called once reading each file will occur
      ! separately to improve simplicity.
      !

      !
      ! Parameter file-stored as non-normalised values
      !

      ! count the number of lines in the file..
      status = 0; num_lines = 0
      do
         read (pfile_unit_thread, iostat=status) dummy
         if (status /= 0) exit
         num_lines = num_lines + 1
      end do
      ! Determine the number of complete parameter vectors stored. Note that the +
      ! 1 is due to the log-likelihood score being saved as well.
      num_lines = num_lines/(npars + 1)

      ! Allocate memory to our temperary variable and the normalised parameter
      ! vector equivalent.
      allocate (tmp(num_lines, (npars + 1)))
      ! rewind so that we can read the contents now correctly
      rewind (pfile_unit_thread)
      ! Read the data for real
      do i = 1, num_lines
         do j = 1, (npars + 1)
            read (pfile_unit_thread) tmp(i, j)
         end do  ! j for parameter
      end do  ! i for combinations

      ! Determine the total number of iterations processed so far
      ! NOTE assuming restart command has same write interval as original simulation
      MCOUT%nos_iterations = num_lines*MCO%nWRITE
      ! Extract the final parameter set and load into the initial parameter vector
      ! for the analysis. NOTE: This parameter set will be normalised on entry
      ! into the MHMCMC subroutine
      MCOUT%pars = tmp(num_lines, 1:npars)

      ! free up variable for new file
      deallocate (tmp)

      !
      ! Variance file-stores output of the current parameter variance
      !

      ! count the number of remaining lines in the file..
      status = 0; num_lines = 0
      do
         read (sfile_unit_thread, iostat=status) dummy
         if (status /= 0.) exit
         num_lines = num_lines + 1
      end do

      ! Determine the number of actual stepsize vectors present.
      ! The+1 is due to the local acceptance rate being provided too.
      num_lines = num_lines/(npars + 1)
      ! allocate memory
      allocate (tmp(num_lines, (npars + 1)))
      ! rewind, for actual reading
      rewind (sfile_unit_thread)

      ! now read the data for real
      do i = 1, num_lines
         do j = 1, (npars + 1)
            read (sfile_unit_thread) tmp(i, j)
         end do  ! j for parameter
      end do  ! i for combinations

      ! Save the current acceptance_rate
      ! TODO logic only applies to APMCMC , which has the quirk of only writing to history when accepted
      MCOUT%acceptance_rate = MCOUT%nos_iterations*MCO%nWRITE

      ! TODO apparently nothing is done with this tmp data - should we save latest step sizes
      ! to a field in MCOUT ?

      ! tidy up for the next file
      deallocate (tmp)

      !
      ! Covariance matrix file
      !

      ! The covariance matrix may contain either 1 or 2 complete matrices. We will
      ! just want to the latest one.

      ! count the number of remaining lines in the file..
      status = 0; num_lines = 1
      do
         read (cfile_unit_thread, iostat=status, rec=num_lines) dummy
         if (status /= 0.) exit
         num_lines = num_lines + 1
      end do

      ! Determine whether there is 1 or more matrice here
      write(*,*) "Cfile unit" , cfile_unit_thread
      write(*,*) "Cov file nlines" , num_lines
      write(*,*) "Cov file nlines/npars/npars" , (num_lines/npars)/npars
      if ((num_lines/npars)/npars == 1) then
         ! the size of the file is consistent with a single matrix having been
         ! saved
         a = 1
      else if ((num_lines/npars)/npars == 2) then
         !
         a = 2
      else
         ! something has gone wrong-best stop
         print *, "Error reading COV file"
         print *, "npars = ", npars, "COV length = ", num_lines*npars
         stop
      end if

      ! now read the data for real
      c = 1
      do b = 1, a
         do i = 1, npars
            do j = 1, npars
               read (cfile_unit_thread, rec=c) MCOUT%covariance(i, j)
               c = c + 1
            end do  ! j for parameter
         end do  ! i for combinations
      end do

      if (a > 0) then
         ! Have at least a first covariance matrix
         ! Set this flag so that the restarted simulation will not overwrite first
         ! covariance matrix in the file with first after restart
         MCOUT%cov = .true.
      endif

      ! extract current variance information
      do i = 1, npars
         MCOUT%parvar(i) = MCOUT%covariance(i, i)
      end do
      ! estimate status of the inverse covariance matrix-iC never used
      ! call inverse_matrix( PI%npars, MCOUT%covariance, PI%iC )

      !
      ! Covariance information file
      !

      ! The number of parameters on which the covariance matrix is based must be
      ! known to allow for correct updating. Similarly the mean normalised
      ! parameter values are also needed

      ! count the number of remaining lines in the file..
      status = 0; num_lines = 0
      do
         read (cifile_unit_thread, iostat=status) dummy
         if (status /= 0.) exit
         num_lines = num_lines + 1
      end do

      ! how many parameter vectors have been output. Note the+1 is accounting
      ! for the number of samples underlying the mean
      num_lines = num_lines/(npars + 1)
      ! allocate memory
      allocate (tmp(num_lines, (npars + 1)))
      ! rewind, for actual reading
      rewind (cifile_unit_thread)
      ! now read the data for real
      do i = 1, num_lines
         do j = 1, (npars + 1)
            read (cifile_unit_thread) tmp(i, j)
         end do  ! j for parameter
      end do  ! i for combinations
      ! Store the most recent step size, which corresponds with the saved
      ! parmeters (above) and covariance matrix (below)
      MCOUT%meanpar = tmp(num_lines, 1:npars)
      MCOUT%Nparvar = tmp(num_lines, npars + 1)

      call close_output_files(pfile_unit_thread, sfile_unit_thread, cfile_unit_thread, cifile_unit_thread)

      return

   end subroutine update_for_restart_simulation
   !
   !------------------------------------------------------------------
   !
   subroutine close_output_files(pfile_unit_thread, sfile_unit_thread, cfile_unit_thread, cifile_unit_thread)

      ! where you open a file you've got to make sure that you close them too. It
      ! just tidy

      implicit none(type, external)
      integer, intent(in) :: pfile_unit_thread, sfile_unit_thread, cfile_unit_thread, cifile_unit_thread

      ! close the files we have in memory
      close (pfile_unit_thread)
      close (sfile_unit_thread)
      close (cfile_unit_thread)
      close (cifile_unit_thread)

   end subroutine close_output_files
   !
   !------------------------------------------------------------------
   !

   subroutine open_output_files(parname, stepname, covname, covinfoname, chainid, pfile_unit_thread, sfile_unit_thread, cfile_unit_thread, cifile_unit_thread)

      ! Subroutine opens the needed output files and destroys any previously
      ! existing files with the same name, just in case mind!
      ! NOTE: that is unless I have not remove the 'UNKNOWN' status in which case
      ! then the files are appended to
      ! last 4 args are for reporting the file unit numbers opened

      implicit none(type, external)

      ! declare input variables
      character(350), intent(in):: parname, stepname, covname, covinfoname
      !! Filenames - In case of multithread should recieve the filenames matching chainid, make sure
      !! to run  filenames_insert_threadid to construct filenames first

      ! declare local variables
      integer:: ios, reclen
      double precision, save:: a = 1d0

      integer, intent(in):: chainid
      integer, intent(out) :: pfile_unit_thread, sfile_unit_thread, cfile_unit_thread, cifile_unit_thread
        !! file unit numbers for this thread

      call calculate_file_ids(chainid, pfile_unit_thread, sfile_unit_thread, cfile_unit_thread, cifile_unit_thread)

      ! open files now
      ! most of these will require new information to be appended to the end at
      ! all times-therefore we use the unformatted stream access
      open (pfile_unit_thread, file=trim(parname), form="UNFORMATTED", access="stream", status="UNKNOWN", iostat=ios)
      if (ios /= 0) print *, "error ", ios, " opening file", trim(parname)
      open (sfile_unit_thread, file=trim(stepname), form="UNFORMATTED", access="stream", status="UNKNOWN", iostat=ios)
      if (ios /= 0) print *, "error ", ios, " opening file", trim(stepname)
      open (cifile_unit_thread, file=trim(covinfoname), form="UNFORMATTED", access="stream", status="UNKNOWN", iostat=ios)
      if (ios /= 0) print *, "error ", ios, " opening file", trim(covinfoname)
      ! for the covariance matrix we have a fixed size containing two matrices,
      ! the initial and the current output-therefore we use
      inquire (iolength=reclen) a !; print*,reclen
      write (*, *) "covname", covname
      open (cfile_unit_thread, file=trim(covname), form="UNFORMATTED", access="direct", recl=reclen, iostat=ios)
      if (ios /= 0) print *, "error ", ios, " opening file", trim(covname)

   end subroutine open_output_files
   !
   !------------------------------------------------------------------
   !
   subroutine initialize_buffers(npars, nwrite_events, io_space)
      integer, intent(in):: npars, nwrite_events
      type(io_buffer_space), intent(inout):: io_space
      io_space%io_buffer_count = 0
      io_space%io_buffer = min(1000, max(10, nwrite_events/10))

      ! Allocate variables used in io buffering,
      ! these could probably be moved to a more sensible place within cardamom_io.f90
      allocate (io_space%variance_buffer(npars, io_space%io_buffer), &
                io_space%meanpars_buffer(npars, io_space%io_buffer), &
                io_space%pars_buffer(npars, io_space%io_buffer), &
                io_space%prob_buffer(io_space%io_buffer), &
                io_space%nsample_buffer(io_space%io_buffer), &
                io_space%accept_rate_buffer(io_space%io_buffer))

      return

   end subroutine initialize_buffers

   !
   !------------------------------------------------------------------
   !
   subroutine write_covariance_matrix(covariance, npars, initial_cov, chainid)

      ! subroutine writes MCMC accepted parameters and step values to binary files

      implicit none(type, external)

      ! arguments
      logical, intent(in):: initial_cov ! Is it the first valid covariance matrix found?
      integer, intent(in):: npars
      double precision, dimension(npars, npars), intent(in):: covariance

      ! declare local variables
      integer:: i, j, irec

      integer, intent(in):: chainid

      integer :: pfile_unit_thread, sfile_unit_thread, cfile_unit_thread, cifile_unit_thread

      call calculate_file_ids(chainid, pfile_unit_thread, sfile_unit_thread, cfile_unit_thread, cifile_unit_thread)

      ! If we have already written the initial covariance matrix we want to keep
      ! over-writing the current matrix. We do this to avoid large files form
      ! writing out multiple covariance matrices
      if (.not. initial_cov) then
         irec = npars*npars
      else
         irec = 0
      end if

      ! write out the file. Its binary format has already been determined at the
      ! openning of the file

      do i = 1, npars
         do j = 1, npars
            irec = irec + 1
            write (cfile_unit_thread , rec=irec) covariance(i, j)
         end do
      end do

      return

   end subroutine write_covariance_matrix
   !
   !------------------------------------------------------------------
   !
   subroutine write_covariance_info(meanpars, nsample, npars, chainid)

      ! subroutine writes MCMC accepted parameters and step values to binary files

      implicit none(type, external)

      ! arguments
      integer, intent(in):: npars
      double precision, intent(in):: nsample
      double precision, dimension(npars), intent(in):: meanpars

      ! declare local variables
      integer:: i

      integer, intent(in):: chainid

      integer :: pfile_unit_thread, sfile_unit_thread, cfile_unit_thread, cifile_unit_thread

      call calculate_file_ids(chainid, pfile_unit_thread, sfile_unit_thread, cfile_unit_thread, cifile_unit_thread)


      ! write out the file. Its binary format has already been determined at the
      ! openning of the file

      do i = 1, npars
         write (cifile_unit_thread ) meanpars(i)
      end do

      write (cifile_unit_thread ) nsample

      return

   end subroutine write_covariance_info
   !
   !------------------------------------------------------------------
   !
   subroutine write_variances(variance, npars, accept_rate, chainid)

      ! subroutine writes parameter variance for corresponding parameter values

      implicit none(type, external)

      ! declare input variables
      integer, intent(in):: npars
      double precision, dimension(npars), intent(in):: variance
      double precision, intent(in):: accept_rate  ! local acceptance rate

      ! declare local variables
      integer:: n

      integer, intent(in):: chainid

      integer :: pfile_unit_thread, sfile_unit_thread, cfile_unit_thread, cifile_unit_thread

      call calculate_file_ids(chainid, pfile_unit_thread, sfile_unit_thread, cfile_unit_thread, cifile_unit_thread)

      ! write out the file. Its binary format has already been determined at the
      ! openning of the file

      do n = 1, npars
         write (sfile_unit_thread ) variance(n)
      end do

      ! we will need to know the current acceptance rate for restarts
      write (sfile_unit_thread ) accept_rate

      return

   end subroutine write_variances
   !
   !------------------------------------------------------------------
   !
   subroutine write_parameters(pars, prob, npars, chainid)

      ! subroutine writes parameter values to binary file`

      implicit none(type, external)

      ! declare input variables
      integer, intent(in):: npars
      double precision, dimension(npars), intent(in):: pars
      double precision, intent(in):: prob

      ! declare local variables
      integer:: n

      integer, intent(in):: chainid
      integer :: pfile_unit_thread, sfile_unit_thread, cfile_unit_thread, cifile_unit_thread

      call calculate_file_ids(chainid, pfile_unit_thread, sfile_unit_thread, cfile_unit_thread, cifile_unit_thread)

      ! write out the file. Its binary format has already been determined at the
      ! openning of the file

      do n = 1, npars
         write (pfile_unit_thread ) pars(n)
      end do

      ! now add the probability
      write (pfile_unit_thread) prob

      ! close will occur at the end of the MCMC

      ! return back
      return

   end subroutine write_parameters
   !
   !--------------------------------------------------------------------
   !
   subroutine write_mcmc_output(variance, accept_rate, &
                                covariance, meanpars, nsample, &
                                pars, prob, npars, dump_now, io_space, chainid)
      !#
      ! Write current state of the sampler to file.
      ! Actually writes to buffer and dumps to file at intervals.
      ! io_space object is a group of buffers spcific to this chain.

      ! Arguments
      integer, intent(in):: npars
      double precision, dimension(npars, npars), intent(in):: covariance
      double precision, dimension(npars), intent(in):: meanpars, &
         variance, &
         pars
      integer, intent(in):: nsample
      double precision, intent(in):: accept_rate, prob
      logical, intent(in):: dump_now
      type(io_buffer_space), intent(inout):: io_space
      integer, intent(in):: chainid

      ! Local variables
      integer:: i

!    ! Debugging print statements
!    print*,"write_mcmc_output:"

      ! Increment buffer
      io_space%io_buffer_count = io_space%io_buffer_count + 1
      ! Store information in buffer for later writing
      io_space%variance_buffer(1:npars, io_space%io_buffer_count) = variance
      io_space%meanpars_buffer(1:npars, io_space%io_buffer_count) = meanpars
      io_space%pars_buffer(1:npars, io_space%io_buffer_count) = pars
      io_space%prob_buffer(io_space%io_buffer_count) = prob
      io_space%nsample_buffer(io_space%io_buffer_count) = nsample
      io_space%accept_rate_buffer(io_space%io_buffer_count) = accept_rate

      ! Are we storing information in buffer or writing to file?
      if (io_space%io_buffer_count == io_space%io_buffer .or. dump_now) then

         ! Then we are writing out to file
         ! Only write the most current covariance matrix as this would be an overwrite anyway
         call write_covariance_matrix(covariance, npars, .false., chainid)
         ! Everything else loop through the buffered output to write out
         do i = 1, io_space%io_buffer_count
            call write_covariance_info(io_space%meanpars_buffer(:, i), io_space%nsample_buffer(i), npars, chainid)
            call write_variances(io_space%variance_buffer(:, i), npars, io_space%accept_rate_buffer(i), chainid)
            call write_parameters(io_space%pars_buffer(:, i), io_space%prob_buffer(i), npars, chainid)
         end do

         ! Reset buffer increment
         io_space%io_buffer_count = 0

      end if

!    ! Debugging print statements
!    print*,"write_mcmc_output:done"

   end subroutine write_mcmc_output
   !
   !--------------------------------------------------------------------
   !
end module samplers_io
