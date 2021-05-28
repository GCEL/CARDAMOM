
module cardamom_io

  !!!!!!!!!!!
  ! Authorship contributions
  !
  ! This code is based on the original C verion of the University of Edinburgh
  ! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
  ! All code translation into Fortran, integration into the University of
  ! Edinburgh CARDAMOM code and subsequent modifications by:
  ! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
  ! J. F. Exbrayat (University of Edinburgh)
  ! See function / subroutine specific comments for exceptions and contributors
  !!!!!!!!!!!

  ! Module contains subroutines and variables needed to output parameter,
  ! likelihood and step size information from the MHMCMC.

  implicit none

  ! declare private
  private

  ! allow access to specific functions
  public :: write_mcmc_output               &
           ,write_parameters                &
           ,write_variances                 &
           ,write_covariance_matrix         &
           ,write_covariance_info           &
           ,update_for_restart_simulation   &
           ,check_for_existing_output_files &
           ,open_output_files               &
           ,close_output_files              &
           ,cardamom_model_library          &
           ,read_pari_data                  &
           ,read_options                    &
           ,read_binary_data

  ! allow access to needed variable
  public :: restart_flag

  ! declare module level variables
  integer :: pfile_unit = 10, sfile_unit = 11, cfile_unit = 12, cifile_unit = 13, ifile_unit = 14
  ! default assumption is that this is not a restart fun
  logical :: restart_flag = .false.

  ! parameters
  integer, parameter :: real_bytes = 8 ! number of bytes in real variable, 8 bytes is to make double precision

  save

  contains
  !
  !--------------------------------------------------------------------
  !
  subroutine cardamom_model_library
    use cardamom_structures, only: DATAin
    implicit none

    ! don't forget to update values found in the relevant model *_PARS.f90

    ! choose between included model arrangements
    ! NOTE: negative values are coded elsewhere and reserved for MCMC stress
    ! testing
    if (DATAin%ID == 0) then
        ! ID = 0 - ACM/ACM-ET
        DATAin%nopools = 2
        DATAin%nopars = 20
        DATAin%nofluxes = 4
    else if (DATAin%ID == 1) then
        ! ID = 1 - DALEC_CDEA
        DATAin%nopools = 6
        DATAin%nopars = 23
        DATAin%nofluxes = 16
    else if (DATAin%ID == 2) then
        ! ID = 2 - DALEC_GSI_BUCKET
        DATAin%nopools = 8
        DATAin%nopars = 40
        DATAin%nofluxes = 25
        if (DATAin%PFT == 1) then
           ! then actually this is a crop pixel
           DATAin%nopools = 9
           DATAin%nopars = 38
           DATAin%nofluxes = 21
        endif
    else if (DATAin%ID == 3 ) then
        ! ID = 3 - AT-DALEC
        DATAin%nopools = 6
        DATAin%nopars = 22
        DATAin%nofluxes = 16
    else if (DATAin%ID == 4) then
        ! AT_DALEC_CROP
        DATAin%nopools = 8
        DATAin%nopars = 34
        DATAin%nofluxes = 16
    else if (DATAin%ID == 5) then
        ! DALEC_CDEA_FR
        DATAin%nopools = 6
        DATAin%nopars = 23
        DATAin%nofluxes = 18
    else if (DATAin%ID == 6) then
        ! DALEC_GSI_FR
        DATAin%nopools = 6
        DATAin%nopars = 33
        DATAin%nofluxes = 18
    else if (DATAin%ID == 7) then
        ! DALEC_GSI_FR_DBio
        DATAin%nopools = 10
        DATAin%nopars = 53
        DATAin%nofluxes = 28
    else if (DATAin%ID == 8) then
        ! DALEC_GSI_MFOL_FR
        DATAin%nopools = 7
        DATAin%nopars = 36
        DATAin%nofluxes = 18
    else if (DATAin%ID == 9) then
        ! DALEC_GSI_FR_LABILE
        DATAin%nopools = 7
        DATAin%nopars = 37
        DATAin%nofluxes = 20
    else if (DATAin%ID == 10) then
        ! DALECN_GSI_FR
        DATAin%nopools = 10
        DATAin%nopars = 49
        DATAin%nofluxes = 21
    else if (DATAin%ID == 11) then
        ! DALEC_GSI_DFOL_FR
        DATAin%nopools = 6
        DATAin%nopars = 36
        DATAin%nofluxes = 18
        if (DATAin%PFT == 1) then
           ! then actually this is a crop pixel
           DATAin%nopools = 8
           DATAin%nopars = 35
           DATAin%nofluxes = 16
        endif
    else if (DATAin%ID == 12) then
        ! DALEC_GSI_DFOL_FROOT_FR
        DATAin%nopools = 9
        DATAin%nopars = 46
        DATAin%nofluxes = 19
    else if (DATAin%ID == 13) then
        ! DALEC_GSI_DFOL_LABILE_FR
        DATAin%nopools = 9
        DATAin%nopars = 44
        DATAin%nofluxes = 19
    else if (DATAin%ID == 14) then
        ! DALECN_GSI_DFOL_LABILE_FR
        DATAin%nopools = 12
        DATAin%nopars = 57
        DATAin%nofluxes = 21
    else if (DATAin%ID == 15) then
        ! DALECN_GSI_DFOL_LABILE_FROOT_FR
        DATAin%nopools = 12
        DATAin%nopars = 61
        DATAin%nofluxes = 21
    else if (DATAin%ID == 16) then
        ! DALEC_GSI_DFOL_CWD_FR
        DATAin%nopools = 7
        DATAin%nopars = 37
        DATAin%nofluxes = 25
        if (DATAin%PFT == 1) then
           ! then actually this is a crop pixel
           DATAin%nopools = 8
           DATAin%nopars = 35
           DATAin%nofluxes = 21
        endif
    else if (DATAin%ID == 17) then
        ! DALECN_GSI_BUCKET - 8 pools currently
        DATAin%nopools = 8
        DATAin%nopars = 45
        DATAin%nofluxes = 25
        if (DATAin%PFT == 1) then
           ! then actually this is a crop pixel
           DATAin%nopools = 9
           DATAin%nopars = 38
           DATAin%nofluxes = 21
        endif
    else if (DATAin%ID == 18) then
        ! DALEC_CDEA_LU_FIRES_ET - added 03/05/2018 JFE
        DATAin%nopools = 6
        DATAin%nopars = 23
        DATAin%nofluxes = 28
        !change ID code below to resolve conflict when merging with jeff = 23/10/18
    else if (DATAin%ID == 19) then
        ! ID = 19 - DALECN_BUCKET
        DATAin%nopools = 8
        DATAin%nopars = 49!45
        DATAin%nofluxes = 25
        if (DATAin%PFT == 1) then
           ! then actually this is a crop pixel
           DATAin%nopools = 9
           DATAin%nopars = 38
           DATAin%nofluxes = 21
        endif
    else if (DATAin%ID == 20) then
        ! ID = 20 - DALEC_BUCKET
        DATAin%nopools = 8
        DATAin%nopars = 43
        DATAin%nofluxes = 25
        if (DATAin%PFT == 1) then
           ! then actually this is a crop pixel
           DATAin%nopools = 9
           DATAin%nopars = 38
           DATAin%nofluxes = 21
        endif
    else if (DATAin%ID == 21) then
        ! ID = 21 - DALEC_CDEA_LU_FIRES
        DATAin%nopools = 6
        DATAin%nopars = 23
        DATAin%nofluxes = 28
    else if (DATAin%ID == 22) then
        ! ID = 22 - DALEC_EVERGREEN
        DATAin%nopools = 5
        DATAin%nopars = 17
        DATAin%nofluxes = 28
    else if (DATAin%ID == 23) then
        ! ID = 23 - DALEC_CDEA_no_lit_root
        DATAin%nopools = 4
        DATAin%nopars = 17
        DATAin%nofluxes = 28
    else if (DATAin%ID == 24) then
        ! ID = 24 - DALEC_EVERGREEN_no_lit_root
        DATAin%nopools = 3
        DATAin%nopars = 11
        DATAin%nofluxes = 28
    else if (DATAin%ID == 25) then
        ! ID = 25 - DALEC_CDEA_ACM2
        DATAin%nopools = 6
        DATAin%nopars = 23
        DATAin%nofluxes = 28
    else if (DATAin%ID == 26) then
        ! ID = 26 - DALEC
        DATAin%nopools = 7
        DATAin%nopars = 43
        DATAin%nofluxes = 25
        if (DATAin%PFT == 1) then
           ! then actually this is a crop pixel
           DATAin%nopools = 8
           DATAin%nopars = 37
           DATAin%nofluxes = 21
        endif
    else if (DATAin%ID == 27) then
        ! ID = 27 - DALEC_CDEA_ACM2_BUCKET
        DATAin%nopools = 7
        DATAin%nopars = 27
        DATAin%nofluxes = 29
    else if (DATAin%ID == 28) then
        ! ID = 28 - DALEC_CDEA_ACM2_BUCKET_RmRg
        DATAin%nopools = 7
        DATAin%nopars = 27
        DATAin%nofluxes = 31
    else if (DATAin%ID == 29) then
        ! ID = 29 - DALEC_CDEA_ACM2_BUCKET_RmRg_CWD
        DATAin%nopools = 8
        DATAin%nopars = 29
        DATAin%nofluxes = 33
    else if (DATAin%ID == 30) then
        ! ID = 30 - DALEC_CDEA_ACM2_BUCKET_RmRg_CWD_wMRT
        DATAin%nopools = 8
        DATAin%nopars = 34
        DATAin%nofluxes = 33
    else if (DATAin%ID == 31) then
        ! ID = 20 - DALEC_BUCKET_CanAGE
        DATAin%nopools = 8
        DATAin%nopars = 48
        DATAin%nofluxes = 25
        if (DATAin%PFT == 1) then
           ! then actually this is a crop pixel
           DATAin%nopools = 9
           DATAin%nopars = 38
           DATAin%nofluxes = 21
        endif
    else if (DATAin%ID == 32) then
        ! ID = 20 - DALEC_G5
        DATAin%nopools = 8
        DATAin%nopars = 46
        DATAin%nofluxes = 25
        if (DATAin%PFT == 1) then
           ! then actually this is a crop pixel
           DATAin%nopools = 9
           DATAin%nopars = 38
           DATAin%nofluxes = 21
        endif
    else if (DATAin%ID == 33) then
        ! ID = 20 - DALEC_G6
        DATAin%nopools = 8
        DATAin%nopars = 46
        DATAin%nofluxes = 25
        if (DATAin%PFT == 1) then
           ! then actually this is a crop pixel
           DATAin%nopools = 9
           DATAin%nopars = 38
           DATAin%nofluxes = 21
        endif
    else
       write(*,*) "Oh dear... model ID cannot be found"
       stop
    endif

  end subroutine cardamom_model_library
  !
  !------------------------------------------------------------------
  !
  subroutine check_for_existing_output_files(npars,nOUT,nWRITE,sub_fraction &
                                            ,parname,stepname,covname,covinfoname)

    ! subroutine checks whether both the parameter and step files exist for this
    ! job. If they do we will assume that this is a restart job that we want to
    ! finish off. Important for large jobs or running on machines with may crash
    ! / have runtime limits
    implicit none
    ! declare input variables
    integer, intent(in) :: npars, nOUT, nWRITE
    double precision, intent(in) :: sub_fraction
    character(350), intent(in) :: parname, stepname, covname, covinfoname
    ! local variables
    logical :: par_exists, step_exists, cov_exists, covinfo_exists
    double precision :: dummy
    integer :: num_lines, status

    ! Check that all files exist
    inquire(file=trim(parname),     exist=par_exists)
    inquire(file=trim(stepname),    exist=step_exists)
    inquire(file=trim(covname),     exist=cov_exists)
    inquire(file=trim(covinfoname), exist=covinfo_exists)

    ! now determine the correct response
    if (par_exists .and. step_exists .and. cov_exists .and. covinfo_exists) then

        ! All files exist therefore this might be a restart run.
        ! lets see if there is anything in the files that we might use
        ! count the number of remaining lines in the file..
        ! open the relevant output files
        call open_output_files(parname,stepname,covname,covinfoname)
        status = 0 ; num_lines = 0
        do
          read(pfile_unit,iostat=status) dummy
          if ( status .ne. 0 ) exit
          num_lines = num_lines + 1
        enddo
        ! Re-use dummy to calculate the target file size to be considered for
        ! restart
        dummy = ((dble(nOUT)/dble(nWRITE)) * sub_fraction) * dble(npars+1)
        if (num_lines > dummy) then
            ! Then there is something in the file we we can use it
            restart_flag = .true.
            print*,"...have found parameter file = ",trim(parname)
            print*,"...have found step file = ",trim(stepname)
            print*,"...have found cov file = ",trim(covname)
            print*,"...have found cov_info file = ",trim(covinfoname)
        else
            ! The file exists but is empty / no enough so treat it as a fresh start
            restart_flag = .false.
            print*,"Output files are present, however they are too small for a restart"
        endif
        ! Either way we open the file up later on so now we need to close them
        call close_output_files

    else ! par_exists .and. step_exists

        ! Then or of these files exists and the other does not so it is
        ! ambiguous whether or not this is a restart
        print*,"One or more of the analysis files cannot be found."
        print*,"CARDAMOM must start from scratch... "
        restart_flag = .false.

    endif ! par_exists .and. step_exists

  end subroutine check_for_existing_output_files
  !
  !------------------------------------------------------------------
  !
  subroutine close_output_files

    ! where you open a file you've got to make sure that you close them too. It
    ! just tidy

    implicit none

    ! close the files we have in memory
    close(pfile_unit)
    close(sfile_unit)
    close(cfile_unit)
    close(cifile_unit)

  end subroutine close_output_files
  !
  !--------------------------------------------------------------------
  !
!  subroutine load_emulator_parameters
!    use cardamom_structures, only: DATAin
!    use CARBON_MODEL_MOD, only: dim_1,dim_2,nos_trees,nos_inputs      &
!                               ,leftDaughter,rightDaughter,nodestatus &
!                               ,xbestsplit,nodepred,bestvar
!
!    ! subroutine opens and reads the PFT specific emulator information needed
!    ! for the randomForest regression trees generated using R package
!    ! randomForest
!    ! 10/10/2014: TLS
!
!    ! TEMPLATE FOR ALL DALEC MCMC DATA files
!    ! Static Elements: 1-100 - use as many as needed
!
!    !STATIC DATA
!    ! 1) PFT
!    ! 2) number of trees in forest
!    ! 3) dimension 1 of response surface (same as interpolation interval)
!    ! 4) dimension 2 of response surface
!    ! 5) number of model inputs needed
!
!    implicit none
!
!    ! declare input variables
!    character(350) :: infile,pft_local
!
!    ! declare local variables
!    integer :: a,i,j,start,finish  &
!              ,ifile_unit   ! unit number assigned to the input binary
!
!    double precision, dimension(:), allocatable :: statdat & ! static data input
!                                                  ,temp_matrix
!
!
!    ! convert PFT into character value for use in file search
!    if (DATAin%PFT < 10) then
!        write(pft_local,fmt='(I1)')DATAin%PFT
!    else if (DATAin%PFT >= 10) then
!        write(pft_local,fmt='(I2)')DATAin%PFT
!    else
!        print*,"Incorrect definition of PFT"
!    endif
!
!    ! assume that parameter files have been copied / linked from the
!    ! AT_DALEC/src
!    ! directory to the execution location
!    write(infile,fmt='(A)')"gpp_emulator_parameters_"//trim(pft_local)//".bin"
!    write(*,*)"Reading emulator coefficients for PFT = ",DATAin%PFT
!    write(*,*)"File path = ",trim(infile)
!
!    ! open the binary file, with direct access for binary (unformatted) at
!    ! double precision (double precision = 64 bytes)
!    open(unit=ifile_unit,file=trim(infile),form="UNFORMATTED",access="stream",status="old")
!    rewind(ifile_unit)
!
!    ! allocate memory
!    allocate(statdat(100))
!
!    ! now read the static elements (1-100)
!    do i = 1, 100 ! number of static elements
!       read(ifile_unit) statdat(i)
!    end do
!
!    ! allocate the default run information
!    nos_trees = int(statdat(2))
!    dim_1 = int(statdat(3))
!    dim_2 = int(statdat(4))
!    nos_inputs = int(statdat(5))
!
!    ! tidy
!    deallocate(statdat)
!
!    ! allocate some other variables
!    allocate(leftDaughter(dim_1,dim_2),rightDaughter(dim_1,dim_2) &
!            ,nodestatus(dim_1,dim_2),xbestsplit(dim_1,dim_2)      &
!            ,nodepred(dim_1,dim_2),bestvar(dim_1,dim_2)           &
!            ,temp_matrix(dim_1*dim_2))
!
!    ! read in left daughter to temp vector
!    a = 1 ; start = 100+1 ; finish = start+(dim_1*dim_2)-1
!    do i = start, finish
!       read(ifile_unit) temp_matrix(a)
!       a = a + 1
!    end do
!    ! restructure left daughter into matrix
!    a = 1
!    do j = 1, dim_2
!       do i = 1, dim_1
!          leftDaughter(i,j)=temp_matrix(a) ; a = a + 1
!       end do
!    end do
!
!    ! read in right daughter into temp vector
!    a = 1 ; start = finish+1 ; finish = start+(dim_1*dim_2)-1
!    do i = start, finish
!       read(ifile_unit) temp_matrix(a)
!       a = a + 1
!    end do
!    ! restructure right daughter into matrix
!    a = 1
!    do j = 1, dim_2
!       do i = 1, dim_1
!          rightDaughter(i,j)=temp_matrix(a) ; a = a + 1
!       end do
!    end do
!
!    ! read in nodestatus into temp vector
!    a = 1 ; start = finish+1 ; finish = start+(dim_1*dim_2)-1
!    do i = start, finish
!       read(ifile_unit) temp_matrix(a)
!       a = a + 1
!    end do
!    ! restructure nodestatus into matrix
!    a = 1
!    do j = 1, dim_2
!       do i = 1, dim_1
!          nodestatus(i,j)=temp_matrix(a) ; a = a + 1
!       end do
!    end do
!
!    ! read in xbestsplit into temp vector
!    a = 1 ; start = finish+1 ; finish = start+(dim_1*dim_2)-1
!    do i = start, finish
!       read(ifile_unit) temp_matrix(a)
!       a = a + 1
!    end do
!    ! restructure xbestsplit into matrix
!    a = 1
!    do j = 1, dim_2
!       do i = 1, dim_1
!          xbestsplit(i,j)=temp_matrix(a) ; a = a + 1
!       end do
!    end do
!
!    ! read in nodepred into temp vector
!    a = 1 ; start = finish+1 ; finish = start+(dim_1*dim_2)-1
!    do i = start, finish
!       read(ifile_unit) temp_matrix(a)
!       a = a + 1
!    end do
!    ! restructure nodepred into matrix
!    a = 1
!    do j = 1, dim_2
!       do i = 1, dim_1
!          nodepred(i,j)=temp_matrix(a) ; a = a + 1
!       end do
!    end do
!
!    ! read in bestvar into temp vector
!    a = 1 ; start = finish+1 ; finish = start+(dim_1*dim_2)-1
!    do i = start, finish
!       read(ifile_unit) temp_matrix(a)
!       a = a + 1
!    end do
!    ! restructure bestvar into matrix
!    a = 1
!    do j = 1, dim_2
!       do i = 1, dim_1
!          bestvar(i,j)=temp_matrix(a) ; a = a + 1
!       end do
!    end do
!
!    ! tidy
!    deallocate(temp_matrix)
!    close(ifile_unit)
!
!    ! inform the user
!    write(*,*)"Have read in GPP emulator coefficients"
!
!  end subroutine load_emulator_parameters
  !
  !--------------------------------------------------------------------
  !
  subroutine open_output_files(parname,stepname,covname,covinfoname)

    ! Subroutine opens the needed output files and destroys any previously
    ! existing files with the same name, just in case mind!
    ! NOTE: that is unless I have not remove the 'UNKNOWN' status in which case
    ! then the files are appended to

    implicit none

    ! declare input variables
    character(350), intent(in) :: parname, stepname, covname,covinfoname

    ! declare local variables
    integer :: ios, reclen
    double precision :: a = 1d0

    ! open files now
    ! most of these will require new information to be appended to the end at
    ! all times - therefore we use the unformatted stream access
    open(pfile_unit,file=trim(parname),form="UNFORMATTED",access="stream",status="UNKNOWN",iostat=ios)
    if (ios /= 0) print*,"error ",ios," opening file",trim(parname)
    open(sfile_unit,file=trim(stepname),form="UNFORMATTED",access="stream",status="UNKNOWN",iostat=ios)
    if (ios /= 0) print*,"error ",ios," opening file",trim(stepname)
    open(cifile_unit,file=trim(covinfoname),form="UNFORMATTED",access="stream",status="UNKNOWN",iostat=ios)
    if (ios /= 0) print*,"error ",ios," opening file",trim(covinfoname)
    ! for the covariance matrix we have a fixed size containing two matrices,
    ! the initial and the current output - therefore we use
    inquire(iolength = reclen) a ; print*,reclen
    open(cfile_unit,file=trim(covname),form="UNFORMATTED",access="direct",recl=reclen,iostat=ios)
    if (ios /= 0) print*,"error ",ios," opening file",trim(covname)

    return

  end subroutine open_output_files
  !
  !--------------------------------------------------------------------
  !
    subroutine read_binary_data(infile)
      use cardamom_structures, only: DATAin
      use CARBON_MODEL_MOD, only: soil_frac_clay,soil_frac_sand &
                                 ,nos_soil_layers

    ! subroutine opens and reads the binary data files provided by / for the
    ! CARDAMOM framework. This data is then loaded into the DATAin type

    ! TEMPLATE FOR ALL DALEC MCMC DATA files
    ! Static Elements: 1-100 - use as many as needed
    ! Parameter Priors: 101-150
    ! Parameter prior uncertainty: 151-200
    ! Other priors & uncertainties: 201-300
    ! TEMPORAL DRIVERS & DATA: 301-end

    implicit none

    ! declare input variables
    character(350) :: infile

    ! declare local variables
    integer :: nopars_dummy,subsample
    integer :: a,b,c,d,e,f,g,h,i,j,k,l,m,o,x,y,z,day,s,t &
              ,start      &
              ,finish     &
              ,totcol     & ! total number of columns (met + obs)
              ,totread      ! total number of records already read
    double precision :: mz, subsample_fraction = 0.20 ! startd at 0.25
    double precision, dimension(:), allocatable :: statdat & ! static data input
                                        ,mettemp & ! met data input
                                        ,obstemp   ! obs data input

    write(*,*)"Input file to be read = ", trim(infile)

    ! open the binary file, with direct access for binary (unformatted) at
    ! double precision (double precision = 64 bytes)
    open(unit=ifile_unit,file=trim(infile),form="UNFORMATTED",access="stream",status="old")
    rewind(ifile_unit)

    ! allocate memory
    allocate(statdat(100))

    ! now read the static elements (1-100)
    do i = 1, 100 ! number of static elements
       read(ifile_unit) statdat(i)
    end do

    ! allocate the default run information
    DATAin%ID = int(statdat(1))
    DATAin%LAT = statdat(2)
    DATAin%nodays = int(statdat(3))
    DATAin%nomet = int(statdat(4))
    DATAin%noobs = int(statdat(5))
    DATAin%EDC = int(statdat(6))
    DATAin%PFT = int(statdat(7))
    DATAin%yield = -9999 !int(statdat(8))
    DATAin%age = int(statdat(9))
    nopars_dummy = int(statdat(10)) ! needed for next dev stage
    if (nos_soil_layers == 3) then
        ! Assume only top soil layer is assigned the top soil condition
        soil_frac_sand(1) = statdat(12) ! top soil sand percentage
        soil_frac_sand(2:nos_soil_layers) = statdat(13) ! bot
        soil_frac_clay(1) = statdat(14) ! top soil clay percentage
        soil_frac_clay(2:nos_soil_layers) = statdat(15) ! bot
    else
        ! Assume that it is greater than 3 layers and thus the top two are considered top
        soil_frac_sand(1:2) = statdat(12) ! top soil sand percentage
        soil_frac_sand(3:nos_soil_layers) = statdat(13) ! bot
        soil_frac_clay(1:2) = statdat(14) ! top soil clay percentage
        soil_frac_clay(3:nos_soil_layers) = statdat(15) ! bot
    end if
    ! call for model specific values
    call cardamom_model_library

    ! allocate case specific information
    DATAin%edc_random_search=int(statdat(11))

    ! clean up
    deallocate(statdat)

    ! read in parameter information
    a = 1
    do i = 101, 200
       read(ifile_unit) DATAin%parpriors(a)
       a = a + 1
    end do

    ! read in parameter uncertainty
    a = 1
    do i = 201, 300
       read(ifile_unit) DATAin%parpriorunc(a)
       a = a + 1
    end do

    ! read in 'other' parameter priors
    a = 1
    do i = 301, 400
       read(ifile_unit) DATAin%otherpriors(a)
       a = a + 1
    end do

    ! read in 'other' parameter priors uncertainties
    a = 1
    do i = 401, 500
       read(ifile_unit) DATAin%otherpriorunc(a)
       a = a + 1
    end do

    ! now we know specific information about the dimensions in the file lets use
    ! it to allocate to the module variables
    allocate(DATAin%met(DATAin%nomet,DATAin%nodays),DATAin%GPP(DATAin%nodays)                    &
            ,DATAin%NEE(DATAin%nodays),DATAin%LAI(DATAin%nodays)                                 &
            ,DATAin%WOO(DATAin%nodays),DATAin%Reco(DATAin%nodays)                                &
            ,DATAin%Cfol_stock(DATAin%nodays),DATAin%Cwood_stock(DATAin%nodays)                  &
            ,DATAin%Croots_stock(DATAin%nodays),DATAin%Clit_stock(DATAin%nodays)                 &
            ,DATAin%Csom_stock(DATAin%nodays),DATAin%Cagb_stock(DATAin%nodays)                   &
            ,DATAin%GPP_unc(DATAin%nodays)                                                       &
            ,DATAin%NEE_unc(DATAin%nodays),DATAin%LAI_unc(DATAin%nodays)                         &
            ,DATAin%WOO_unc(DATAin%nodays),DATAin%Reco_unc(DATAin%nodays)                        &
            ,DATAin%Cfol_stock_unc(DATAin%nodays),DATAin%Cwood_stock_unc(DATAin%nodays)          &
            ,DATAin%Croots_stock_unc(DATAin%nodays),DATAin%Clit_stock_unc(DATAin%nodays)         &
            ,DATAin%Csom_stock_unc(DATAin%nodays),DATAin%Cagb_stock_unc(DATAin%nodays)           &
            ,DATAin%Ccoarseroot_stock(DATAin%nodays),DATAin%Ccoarseroot_stock_unc(DATAin%nodays) &
            ,DATAin%Cfolmax_stock(DATAin%nodays),DATAin%Cfolmax_stock_unc(DATAin%nodays)         &
            ,DATAin%Evap(DATAin%nodays),DATAin%Evap_unc(DATAin%nodays)                           &
            ,DATAin%SWE(DATAin%nodays),DATAin%SWE_unc(DATAin%nodays)                             &
            ,DATAin%NBE(DATAin%nodays),DATAin%NBE_unc(DATAin%nodays)                             &
            ,mettemp(DATAin%nomet),obstemp(DATAin%noobs))

    ! zero all variables
    DATAin%met = 0d0 ; DATAin%GPP = 0d0 ; DATAin%NEE = 0d0 ; DATAin%LAI = 0d0 ; DATAin%WOO = 0d0 ; DATAin%Reco = 0d0
    DATAin%Cfol_stock = 0d0 ; DATAin%Cwood_stock = 0d0 ; DATAin%Croots_stock = 0d0 ; DATAin%Clit_stock = 0d0
    DATAin%Csom_stock = 0d0 ; DATAin%Cagb_stock = 0d0 ; DATAin%GPP_unc = 0d0 ; DATAin%NEE_unc = 0d0
    DATAin%LAI_unc = 0d0 ; DATAin%WOO_unc = 0d0 ; DATAin%Reco_unc = 0d0 ; DATAin%Cfol_stock_unc = 0d0
    DATAin%Cwood_stock_unc = 0d0 ; DATAin%Croots_stock_unc = 0d0 ; DATAin%Clit_stock_unc = 0d0
    DATAin%Csom_stock_unc = 0d0 ; DATAin%Cagb_stock_unc = 0d0
    DATAin%Ccoarseroot_stock = 0d0 ; DATAin%Ccoarseroot_stock_unc = 0d0
    DATAin%Cfolmax_stock = 0d0 ; DATAin%Cfolmax_stock_unc = 0d0
    DATAin%Evap = 0d0 ; DATAin%Evap_unc = 0d0
    DATAin%SWE = 0d0 ; DATAin%SWE_unc = 0d0
    DATAin%NBE = 0d0 ; DATAin%NBE_unc = 0d0
    mettemp = 0d0 ; obstemp = 0d0

    ! zero the obs counters
    DATAin%total_obs = 0
    DATAin%ngpp = 0               ; DATAin%sub_ngpp = 0
    DATAin%nlai = 0               ; DATAin%sub_nlai = 0
    DATAin%nnee = 0               ; DATAin%sub_nnee = 0
    DATAin%nwoo = 0               ; DATAin%sub_nwoo = 0
    DATAin%nreco = 0              ; DATAin%sub_nreco = 0
    DATAin%nCfol_stock = 0        ; DATAin%sub_nCfol_stock = 0
    DATAin%nCwood_stock = 0       ; DATAin%sub_nCwood_stock = 0
    DATAin%nCroots_stock = 0      ; DATAin%sub_nCroots_stock = 0
    DATAin%nCsom_stock = 0        ; DATAin%sub_nCsom_stock = 0
    DATAin%nClit_stock = 0        ; DATAin%sub_nClit_stock = 0
    DATAin%nCagb_stock = 0        ; DATAin%sub_nCagb_stock = 0
    DATAin%nCcoarseroot_stock = 0 ; DATAin%sub_nCcoarseroot_stock = 0
    DATAin%nCfolmax_stock = 0     ; DATAin%sub_nCfolmax_stock = 0
    DATAin%nEvap = 0              ; DATAin%sub_nEvap = 0
    DATAin%nSWE = 0               ; DATAin%sub_nSWE = 0
    DATAin%nNBE = 0               ; DATAin%sub_nNBE = 0

    ! work out some key variables
    ! DATAin%noobs corresponds to observations and uncertainties
    totcol = DATAin%nomet*DATAin%noobs
    totread = 500+1

    ! start looping through days and allocate the correct met drivers / obs into
    ! the correct arrays
    do day = 1, DATAin%nodays
       start = ((day-1)*(totcol))+totread
       finish = start+DATAin%nomet-1
       b = 1
       do i = start,finish
          read(ifile_unit) mettemp(b) ; b=b+1
       end do ! met bit
       start = ((day-1)*(totcol))+totread+finish+1
       finish = start+DATAin%noobs-1
       b = 1
       do i = start,finish
          read(ifile_unit) obstemp(b) ; b=b+1
       end do ! obs bit
       ! assign the extracted met / obs to their type and keep count of how many
       ! of these are actually contain data
       DATAin%met(1:DATAin%nomet,day) = mettemp
       ! JFE change 17/04/18 to read each obs and corresponding uncertainty
       ! we will assume that where we have an observation we also have a
       ! corresponding uncertainty estimate

       DATAin%GPP(day) = obstemp(1)
       if (obstemp(1) > -9998d0) DATAin%ngpp = DATAin%ngpp+1
       DATAin%GPP_unc(day) = obstemp(2)

       DATAin%LAI(day) = obstemp(3)
       if (obstemp(3) > -9998d0) DATAin%nlai = DATAin%nlai+1
       DATAin%LAI_unc(day) = obstemp(4)

       DATAin%NEE(day) = obstemp(5)
       if (obstemp(5) > -9998d0) DATAin%nnee = DATAin%nnee+1
       DATAin%NEE_unc(day) = obstemp(6)

       DATAin%WOO(day) = obstemp(7)
       if (obstemp(7) > -9998d0) DATAin%nwoo = DATAin%nwoo+1
       DATAin%WOO_unc(day) = obstemp(8)

       DATAin%Reco(day) = obstemp(9)
       if (obstemp(9) > -9998d0) DATAin%nreco = DATAin%nreco+1
       DATAin%Reco_unc(day) = obstemp(10)

       DATAin%Cfol_stock(day) = obstemp(11)
       if (obstemp(11) > -9998d0) DATAin%nCfol_stock = DATAin%nCfol_stock+1
       DATAin%Cfol_stock_unc(day) = obstemp(12)

       DATAin%Cwood_stock(day) = obstemp(13)
       if (obstemp(13) > -9998d0) DATAin%nCwood_stock = DATAin%nCwood_stock+1
       DATAin%Cwood_stock_unc(day) = obstemp(14)

       DATAin%Croots_stock(day) = obstemp(15)
       if (obstemp(15) > -9998d0) DATAin%nCroots_stock = DATAin%nCroots_stock+1
       DATAin%Croots_stock_unc(day) = obstemp(16)

       DATAin%Clit_stock(day) = obstemp(17)
       if (obstemp(17) > -9998d0) DATAin%nClit_stock = DATAin%nClit_stock+1
       DATAin%Clit_stock_unc(day) = obstemp(18)

       DATAin%Csom_stock(day) = obstemp(19)
       if (obstemp(19) > -9998d0) DATAin%nCsom_stock = DATAin%nCsom_stock+1
       DATAin%Csom_stock_unc(day) = obstemp(20)

       DATAin%Cagb_stock(day) = obstemp(21)
       if (obstemp(21) > -9998d0) DATAin%nCagb_stock = DATAin%nCagb_stock+1
       DATAin%Cagb_stock_unc(day) = obstemp(22)

! POSITION 23-26 no longer have matching points in code.
! These can be re-allocated at a future point
! TLS: 27/11/2019
!       DATAin%Cstem_stock(day) = obstemp(23)
!       if (obstemp(23) > -9998d0) DATAin%nCstem_stock = DATAin%nCstem_stock+1
!       DATAin%Cstem_stock_unc(day) = obstemp(24)
!       DATAin%Cbranch_stock(day) = obstemp(25)
!       if (obstemp(25) > -9998d0) DATAin%nCbranch_stock = DATAin%nCbranch_stock+1
!       DATAin%Cbranch_stock_unc(day) = obstemp(26)

       DATAin%Ccoarseroot_stock(day) = obstemp(27)
       if (obstemp(27) > -9998d0) DATAin%nCcoarseroot_stock = DATAin%nCcoarseroot_stock+1
       DATAin%Ccoarseroot_stock_unc(day) = obstemp(28)

       DATAin%Cfolmax_stock(day) = obstemp(29)
       if (obstemp(29) > -9998d0) DATAin%nCfolmax_stock = DATAin%nCfolmax_stock+1
       DATAin%Cfolmax_stock_unc(day) = obstemp(30)

       DATAin%Evap(day) = obstemp(31)
       if (obstemp(31) > -9998d0) DATAin%nEvap = DATAin%nEvap+1
       DATAin%Evap_unc(day) = obstemp(32)

       !added SWE for future applications
       DATAin%SWE(day) = obstemp(33)
       if (obstemp(33) > -9998d0) DATAin%nSWE = DATAin%nSWE+1
       DATAin%SWE_unc(day) = obstemp(34)

       DATAin%NBE(day) = obstemp(35)
       if (obstemp(35) > -9998d0) DATAin%nNBE = DATAin%nNBE+1
       DATAin%NBE_unc(day) = obstemp(36)

    end do ! day loop

    ! Count the total number of observations which are to be used.
    ! This total in some models may be used to inform on a dynamic weighting of the EDCs
    DATAin%total_obs = DATAin%ngpp + DATAin%nlai + DATAin%nnee &
                     + DATAin%nwoo + DATAin%nreco + DATAin%nCfol_stock &
                     + DATAin%nCwood_stock + DATAin%nCroots_stock + DATAin%nCsom_stock &
                     + DATAin%nClit_stock + DATAin%nCagb_stock + DATAin%nCcoarseroot_stock &
                     + DATAin%nCfolmax_stock + DATAin%nEvap + DATAin%nSWE + DATAin%nNBE

    ! allocate to time step
    allocate(DATAin%deltat(DATAin%nodays)) ; DATAin%deltat = 0d0
    ! work out interval in decimal days
    DATAin%deltat(1) = DATAin%met(1,1)
    do i = 2, DATAin%nodays
       DATAin%deltat(i) = DATAin%met(1,i)-DATAin%met(1,(i-1))
    end do

    ! close open binary
    close(ifile_unit)

    ! clean up memory
    deallocate(obstemp,mettemp)

    ! allocate memory to our data location variables
    if (DATAin%ngpp > 0) allocate(DATAin%gpppts(DATAin%ngpp))
    if (DATAin%nlai > 0) allocate(DATAin%laipts(DATAin%nlai))
    if (DATAin%nnee > 0) allocate(DATAin%neepts(DATAin%nnee))
    if (DATAin%nwoo > 0) allocate(DATAin%woopts(DATAin%nwoo))
    if (DATAin%nreco > 0) allocate(DATAin%recopts(DATAin%nreco))
    if (DATAin%nCfol_stock > 0) allocate(DATAin%Cfol_stockpts(DATAin%nCfol_stock))
    if (DATAin%nCwood_stock > 0) allocate(DATAin%Cwood_stockpts(DATAin%nCwood_stock))
    if (DATAin%nCroots_stock > 0) allocate(DATAin%Croots_stockpts(DATAin%nCroots_stock))
    if (DATAin%nClit_stock > 0) allocate(DATAin%Clit_stockpts(DATAin%nClit_stock))
    if (DATAin%nCsom_stock > 0) allocate(DATAin%Csom_stockpts(DATAin%nCsom_stock))
    if (DATAin%nCagb_stock > 0) allocate(DATAin%Cagb_stockpts(DATAin%nCagb_stock))
    if (DATAin%nCcoarseroot_stock > 0) allocate(DATAin%Ccoarseroot_stockpts(DATAin%nCcoarseroot_stock))
    if (DATAin%nCfolmax_stock > 0) allocate(DATAin%Cfolmax_stockpts(DATAin%nCfolmax_stock))
    if (DATAin%nEvap > 0) allocate(DATAin%Evappts(DATAin%nEvap))
    if (DATAin%nSWE > 0) allocate(DATAin%SWEpts(DATAin%nSWE))
    if (DATAin%nNBE > 0) allocate(DATAin%NBEpts(DATAin%nNBE))
    ! we know how many observations we have and what they are, but now lets work
    ! out where they are in the data sets
    x = 1 ; y = 1 ; z = 1 ; b = 1 ; c = 1 ; d = 1 ; e = 1
    f = 1 ; g = 1 ; h = 1 ; i = 1 ; j = 1 ; k = 1 ; l = 1
    m = 1 ; o = 1 ; s = 1 ; t = 1
    do day = 1, DATAin%nodays
       if (DATAin%GPP(day) > -9998d0) then
          DATAin%gpppts(b) = day ; b = b+1
       endif
       if (DATAin%LAI(day) > -9998d0) then
          DATAin%laipts(x) = day ; x = x+1
       endif
       if (DATAin%NEE(day) > -9998d0) then
          DATAin%neepts(y) = day ; y = y+1
       endif
       if (DATAin%WOO(day) > -9998d0) then
           DATAin%woopts(z) = day ; z = z+1
       endif ! data present condition
       if (DATAin%Reco(day) > -9998d0) then
           DATAin%recopts(c) = day ; c = c+1
       endif ! data present condition
       if (DATAin%Cfol_stock(day) > -9998d0) then
           DATAin%Cfol_stockpts(d) = day ; d = d+1
       endif ! data present condition
       if (DATAin%Cwood_stock(day) > -9998d0) then
           DATAin%Cwood_stockpts(e) = day ; e = e+1
       endif ! data present condition
       if (DATAin%Croots_stock(day) > -9998d0) then
           DATAin%Croots_stockpts(f) = day ; f = f+1
       endif ! data present condition
       if (DATAin%Clit_stock(day) > -9998d0) then
           DATAin%Clit_stockpts(j) = day ; j = j+1
       endif ! data present condition
       if (DATAin%Csom_stock(day) > -9998d0) then
           DATAin%Csom_stockpts(k) = day ; k = k+1
       endif ! data present condition
       if (DATAin%Ccoarseroot_stock(day) > -9998d0) then
           DATAin%Ccoarseroot_stockpts(i) = day ; i = i+1
       endif ! data present condition
       if (DATAin%Cfolmax_stock(day) > -9998d0) then
           DATAin%Cfolmax_stockpts(l) = day ; l = l+1
       endif ! data present condition
       if (DATAin%Cagb_stock(day) > -9998d0) then
           DATAin%Cagb_stockpts(l) = day ; l = l+1
       endif ! data present condition
       if (DATAin%Evap(day) > -9998d0) then
           DATAin%Evappts(o) = day ; o = o+1
       endif ! data present condition
       if (DATAin%SWE(day) > -9998d0) then
           DATAin%SWEpts(s) = day ; s = s+1
       endif ! data present condition
       if (DATAin%NBE(day) > -9998d0) then
           DATAin%NBEpts(t) = day ; t = t + 1
       endif
    end do ! day loop

    ! Determine a subsample of observations. We assume that we want observations
    ! covering only 25 % of time steps across all data streams. NOTE: the
    ! scaling by number of observations per datastream over the total and
    ! limited to a total of 100 observations per datastream irrespective.

!    if (DATAin%ngpp > 0) then
!        ! Estimate the maximum number of observations per data stream
!        subsample = ceiling(dble(DATAin%nodays) * subsample_fraction * (dble(DATAin%ngpp)/dble(DATAin%total_obs)))
!        subsample = min(subsample,100)
!        ! Allocate which ever is smaller
!        DATAin%sub_ngpp = min(subsample,DATAin%ngpp)
!        allocate(DATAin%sub_gpppts(DATAin%sub_ngpp))
!        ! Estimate step between draws following: ((to - from)/(length.out - 1))
!        mz = (dble(DATAin%ngpp-1)/dble((DATAin%sub_ngpp-1)))**(-1d0)
!        if (DATAin%ngpp == 1 .or. DATAin%sub_ngpp == 1) mz = 1d0 ! Guard against NaN
!        do i = 1, DATAin%sub_ngpp
!           ! Estimate the index of the next value we want
!           DATAin%sub_gpppts(i) = DATAin%gpppts(nint((dble(i-1)/mz)+1))
!        end do
!        print*,"GPP subsample = ",DATAin%sub_ngpp," observations of ",DATAin%ngpp
!    end if
!
!    if (DATAin%nlai > 0) then
!        ! Estimate the maximum number of observations per data stream
!        subsample = ceiling(dble(DATAin%nodays) * subsample_fraction * (dble(DATAin%nlai)/dble(DATAin%total_obs)))
!        subsample = min(subsample,100)
!        ! Allocate which ever is smaller
!        DATAin%sub_nlai = min(subsample,DATAin%nlai)
!        allocate(DATAin%sub_laipts(DATAin%sub_nlai))
!        ! Estimate step between draws following: ((to - from)/(length.out - 1))
!        mz = (dble(DATAin%nlai-1)/dble((DATAin%sub_nlai-1)))**(-1d0)
!        if (DATAin%nlai == 1 .or. DATAin%sub_nlai == 1) mz = 1d0 ! Guard against NaN
!        do i = 1, DATAin%sub_nlai
!           ! Estimate the index of the next value we want
!           DATAin%sub_laipts(i) = DATAin%laipts(nint((dble(i-1)/mz)+1))
!        end do
!        print*,"LAI subsample = ",DATAin%sub_nlai," observations of ",DATAin%nlai
!    end if
!
!    if (DATAin%nnee > 0) then
!        ! Estimate the maximum number of observations per data stream
!        subsample = ceiling(dble(DATAin%nodays) * subsample_fraction * (dble(DATAin%nnee)/dble(DATAin%total_obs)))
!        subsample = min(subsample,100)
!        ! Allocate which ever is smaller
!        DATAin%sub_nnee = min(subsample,DATAin%nnee)
!        allocate(DATAin%sub_neepts(DATAin%sub_nnee))
!        ! Estimate step between draws following: ((to - from)/(length.out - 1))
!        mz = (dble(DATAin%nnee-1)/dble((DATAin%sub_nnee-1)))**(-1d0)
!        if (DATAin%nnee == 1 .or. DATAin%sub_nnee == 1) mz = 1d0 ! Guard against NaN
!        do i = 1, DATAin%sub_nnee
!           ! Estimate the index of the next value we want
!           DATAin%sub_neepts(i) = DATAin%neepts(nint((dble(i-1)/mz)+1))
!        end do
!        print*,"NEE subsample = ",DATAin%sub_nnee," observations of ",DATAin%nnee
!    end if
!
!    if (DATAin%nwoo > 0) then
!        ! Estimate the maximum number of observations per data stream
!        subsample = ceiling(dble(DATAin%nodays) * subsample_fraction * (dble(DATAin%nwoo)/dble(DATAin%total_obs)))
!        subsample = min(subsample,100)
!        ! Allocate which ever is smaller
!        DATAin%sub_nwoo = min(subsample,DATAin%nwoo)
!        allocate(DATAin%sub_woopts(DATAin%sub_nwoo))
!        ! Estimate step between draws following: ((to - from)/(length.out - 1))
!        mz = (dble(DATAin%nwoo-1)/dble((DATAin%sub_nwoo-1)))**(-1d0)
!        if (DATAin%nwoo == 1 .or. DATAin%sub_nwoo == 1) mz = 1d0 ! Guard against NaN
!        do i = 1, DATAin%sub_nwoo
!           ! Estimate the index of the next value we want
!           DATAin%sub_woopts(i) = DATAin%woopts(nint((dble(i-1)/mz)+1))
!        end do
!        print*,"WOO subsample = ",DATAin%sub_nwoo," observations of ",DATAin%nwoo
!    end if
!
!    if (DATAin%nreco > 0) then
!        ! Estimate the maximum number of observations per data stream
!        subsample = ceiling(dble(DATAin%nodays) * subsample_fraction * (dble(DATAin%nreco)/dble(DATAin%total_obs)))
!        subsample = min(subsample,100)
!        ! Allocate which ever is smaller
!        DATAin%sub_nreco = min(subsample,DATAin%nreco)
!        allocate(DATAin%sub_recopts(DATAin%sub_nreco))
!        ! Estimate step between draws following: ((to - from)/(length.out - 1))
!        mz = (dble(DATAin%nreco-1)/dble((DATAin%sub_nreco-1)))**(-1d0)
!        if (DATAin%nreco == 1 .or. DATAin%sub_nreco == 1) mz = 1d0 ! Guard against NaN
!        do i = 1, DATAin%sub_nreco
!           ! Estimate the index of the next value we want
!           DATAin%sub_recopts(i) = DATAin%recopts(nint((dble(i-1)/mz)+1))
!        end do
!        print*,"Reco subsample = ",DATAin%sub_nreco," observations of ",DATAin%nreco
!    end if
!
!    if (DATAin%nCfol_stock > 0) then
!        ! Estimate the maximum number of observations per data stream
!        subsample = ceiling(dble(DATAin%nodays) * subsample_fraction * (dble(DATAin%nCfol_stock)/dble(DATAin%total_obs)))
!        subsample = min(subsample,100)
!        ! Allocate which ever is smaller
!        DATAin%sub_nCfol_stock = min(subsample,DATAin%nCfol_stock)
!        allocate(DATAin%sub_Cfol_stockpts(DATAin%sub_nCfol_stock))
!        ! Estimate step between draws following: ((to - from)/(length.out - 1))
!        mz = (dble(DATAin%nCfol_stock-1)/dble((DATAin%sub_nCfol_stock-1)))**(-1d0)
!        if (DATAin%nCfol_stock == 1 .or. DATAin%sub_nCfol_stock == 1) mz = 1d0 ! Guard against NaN
!        do i = 1, DATAin%sub_nCfol_stock
!           ! Estimate the index of the next value we want
!           DATAin%sub_Cfol_stockpts(i) = DATAin%Cfol_stockpts(nint((dble(i-1)/mz)+1))
!        end do
!        print*,"Cfol_stock subsample = ",DATAin%sub_nCfol_stock," observations of ",DATAin%nCfol_stock
!    end if
!
!    if (DATAin%nCwood_stock > 0) then
!        ! Estimate the maximum number of observations per data stream
!        subsample = ceiling(dble(DATAin%nodays) * subsample_fraction * (dble(DATAin%nCwood_stock)/dble(DATAin%total_obs)))
!        subsample = min(subsample,100)
!        ! Allocate which ever is smaller
!        DATAin%sub_nCwood_stock = min(subsample,DATAin%nCwood_stock)
!        allocate(DATAin%sub_Cwood_stockpts(DATAin%sub_nCwood_stock))
!        ! Estimate step between draws following: ((to - from)/(length.out - 1))
!        mz = (dble(DATAin%nCwood_stock-1)/dble((DATAin%sub_nCwood_stock-1)))**(-1d0)
!        if (DATAin%nCwood_stock == 1 .or. DATAin%sub_nCwood_stock == 1) mz = 1d0 ! Guard against NaN
!        do i = 1, DATAin%sub_nCwood_stock
!           ! Estimate the index of the next value we want
!           DATAin%sub_Cwood_stockpts(i) = DATAin%Cwood_stockpts(nint((dble(i-1)/mz)+1))
!        end do
!        print*,"Cwood_stock subsample = ",DATAin%sub_nCwood_stock," observations of ",DATAin%nCwood_stock
!    end if
!
!    if (DATAin%nCroots_stock > 0) then
!        ! Estimate the maximum number of observations per data stream
!        subsample = ceiling(dble(DATAin%nodays) * subsample_fraction * (dble(DATAin%nCroots_stock)/dble(DATAin%total_obs)))
!        subsample = min(subsample,100)
!        ! Allocate which ever is smaller
!        DATAin%sub_nCroots_stock = min(subsample,DATAin%nCroots_stock)
!        allocate(DATAin%sub_Croots_stockpts(DATAin%sub_nCroots_stock))
!        ! Estimate step between draws following: ((to - from)/(length.out - 1))
!        mz = (dble(DATAin%nCroots_stock-1)/dble((DATAin%sub_nCroots_stock-1)))**(-1d0)
!        if (DATAin%nCroots_stock == 1 .or. DATAin%sub_nCroots_stock == 1) mz = 1d0 ! Guard against NaN
!        do i = 1, DATAin%sub_nCroots_stock
!           ! Estimate the index of the next value we want
!           DATAin%sub_Croots_stockpts(i) = DATAin%Croots_stockpts(nint((dble(i-1)/mz)+1))
!        end do
!        print*,"Croots_stock subsample = ",DATAin%sub_nCroots_stock," observations of ",DATAin%nCroots_stock
!    end if
!
!    if (DATAin%nClit_stock > 0) then
!        ! Estimate the maximum number of observations per data stream
!        subsample = ceiling(dble(DATAin%nodays) * subsample_fraction * (dble(DATAin%nClit_stock)/dble(DATAin%total_obs)))
!        subsample = min(subsample,100)
!        ! Allocate which ever is smaller
!        DATAin%sub_nClit_stock = min(subsample,DATAin%nClit_stock)
!        allocate(DATAin%sub_Clit_stockpts(DATAin%sub_nClit_stock))
!        ! Estimate step between draws following: ((to - from)/(length.out - 1))
!        mz = (dble(DATAin%nClit_stock-1)/dble((DATAin%sub_nClit_stock-1)))**(-1d0)
!        if (DATAin%nClit_stock == 1 .or. DATAin%sub_nClit_stock == 1) mz = 1d0 ! Guard against NaN
!        do i = 1, DATAin%sub_nClit_stock
!           ! Estimate the index of the next value we want
!           DATAin%sub_Clit_stockpts(i) = DATAin%Clit_stockpts(nint((dble(i-1)/mz)+1))
!        end do
!        print*,"Clit_stock subsample = ",DATAin%sub_nClit_stock," observations of ",DATAin%nClit_stock
!    end if
!
!    if (DATAin%nCsom_stock > 0) then
!        ! Estimate the maximum number of observations per data stream
!        subsample = ceiling(dble(DATAin%nodays) * subsample_fraction * (dble(DATAin%nCsom_stock)/dble(DATAin%total_obs)))
!        subsample = min(subsample,100)
!        ! Allocate which ever is smaller
!        DATAin%sub_nCsom_stock = min(subsample,DATAin%nCsom_stock)
!        allocate(DATAin%sub_Csom_stockpts(DATAin%sub_nCsom_stock))
!        ! Estimate step between draws following: ((to - from)/(length.out - 1))
!        mz = (dble(DATAin%nCsom_stock-1)/dble((DATAin%sub_nCsom_stock-1)))**(-1d0)
!        if (DATAin%nCsom_stock == 1 .or. DATAin%sub_nCsom_stock == 1) mz = 1d0 ! Guard against NaN
!        do i = 1, DATAin%sub_nCsom_stock
!           ! Estimate the index of the next value we want
!           DATAin%sub_Csom_stockpts(i) = DATAin%Csom_stockpts(nint((dble(i-1)/mz)+1))
!        end do
!        print*,"Csom_stock subsample = ",DATAin%sub_nCsom_stock," observations of ",DATAin%nCsom_stock
!    end if
!
!    if (DATAin%nCagb_stock > 0) then
!        ! Estimate the maximum number of observations per data stream
!        subsample = ceiling(dble(DATAin%nodays) * subsample_fraction * (dble(DATAin%nCagb_stock)/dble(DATAin%total_obs)))
!        subsample = min(subsample,100)
!        ! Allocate which ever is smaller
!        DATAin%sub_nCagb_stock = min(subsample,DATAin%nCagb_stock)
!        allocate(DATAin%sub_Cagb_stockpts(DATAin%sub_nCagb_stock))
!        ! Estimate step between draws following: ((to - from)/(length.out - 1))
!        mz = (dble(DATAin%nCagb_stock-1)/dble((DATAin%sub_nCagb_stock-1)))**(-1d0)
!        if (DATAin%nCagb_stock == 1 .or. DATAin%sub_nCagb_stock == 1) mz = 1d0 ! Guard against NaN
!        do i = 1, DATAin%sub_nCagb_stock
!           ! Estimate the index of the next value we want
!           DATAin%sub_Cagb_stockpts(i) = DATAin%Cagb_stockpts(nint((dble(i-1)/mz)+1))
!        end do
!        print*,"Cagb_stock subsample = ",DATAin%sub_nCagb_stock," observations of ",DATAin%nCagb_stock
!    end if
!
!    if (DATAin%nCcoarseroot_stock > 0) then
!        ! Estimate the maximum number of observations per data stream
!        subsample = ceiling(dble(DATAin%nodays) * subsample_fraction * (dble(DATAin%nCcoarseroot_stock)/dble(DATAin%total_obs)))
!        subsample = min(subsample,100)
!        ! Allocate which ever is smaller
!        DATAin%sub_nCcoarseroot_stock = min(subsample,DATAin%nCcoarseroot_stock)
!        allocate(DATAin%sub_Ccoarseroot_stockpts(DATAin%sub_nCcoarseroot_stock))
!        ! Estimate step between draws following: ((to - from)/(length.out - 1))
!        mz = (dble(DATAin%nCcoarseroot_stock-1)/dble((DATAin%sub_nCcoarseroot_stock-1)))**(-1d0)
!        if (DATAin%nCcoarseroot_stock == 1 .or. DATAin%sub_nCcoarseroot_stock == 1) mz = 1d0 ! Guard against NaN
!        do i = 1, DATAin%sub_nCcoarseroot_stock
!           ! Estimate the index of the next value we want
!           DATAin%sub_Ccoarseroot_stockpts(i) = DATAin%Ccoarseroot_stockpts(nint((dble(i-1)/mz)+1))
!        end do
!        print*,"Ccoarseroot_stock subsample = ",DATAin%sub_nCcoarseroot_stock," observations of ",DATAin%nCcoarseroot_stock
!    end if
!
!    if (DATAin%nCfolmax_stock > 0) then
!        ! Estimate the maximum number of observations per data stream
!        subsample = ceiling(dble(DATAin%nodays) * subsample_fraction * (dble(DATAin%nCfolmax_stock)/dble(DATAin%total_obs)))
!        subsample = min(subsample,100)
!        ! Allocate which ever is smaller
!        DATAin%sub_nCfolmax_stock = min(subsample,DATAin%nCfolmax_stock)
!        allocate(DATAin%sub_Cfolmax_stockpts(DATAin%sub_nCfolmax_stock))
!        ! Estimate step between draws following: ((to - from)/(length.out - 1))
!        mz = (dble(DATAin%nCfolmax_stock-1)/dble((DATAin%sub_nCfolmax_stock-1)))**(-1d0)
!        if (DATAin%nCfolmax_stock == 1 .or. DATAin%sub_nCfolmax_stock == 1) mz = 1d0 ! Guard against NaN
!        do i = 1, DATAin%sub_nCfolmax_stock
!           ! Estimate the index of the next value we want
!           DATAin%sub_Cfolmax_stockpts(i) = DATAin%Cfolmax_stockpts(nint((dble(i-1)/mz)+1))
!        end do
!        print*,"Cfol_max subsample = ",DATAin%sub_nCfolmax_stock," observations of ",DATAin%nCfolmax_stock
!    end if
!
!    if (DATAin%nEvap > 0) then
!        ! Estimate the maximum number of observations per data stream
!        subsample = ceiling(dble(DATAin%nodays) * subsample_fraction * (dble(DATAin%nEvap)/dble(DATAin%total_obs)))
!        subsample = min(subsample,100)
!        ! Allocate which ever is smaller
!        DATAin%sub_nEvap = min(subsample,DATAin%nEvap)
!        allocate(DATAin%sub_Evappts(DATAin%sub_nEvap))
!        ! Estimate step between draws following: ((to - from)/(length.out - 1))
!        mz = (dble(DATAin%nEvap-1)/dble((DATAin%sub_nEvap-1)))**(-1d0)
!        if (DATAin%nEvap == 1 .or. DATAin%sub_nEvap == 1) mz = 1d0 ! Guard against NaN
!        do i = 1, DATAin%sub_nEvap
!           ! Estimate the index of the next value we want
!           DATAin%sub_Evappts(i) = DATAin%Evappts(nint((dble(i-1)/mz)+1))
!        end do
!        print*,"Evap subsample = ",DATAin%sub_nEvap," observations of ",DATAin%nEvap
!    end if
!
!    if (DATAin%nSWE > 0) then
!        ! Estimate the maximum number of observations per data stream
!        subsample = ceiling(dble(DATAin%nodays) * subsample_fraction * (dble(DATAin%nSWE)/dble(DATAin%total_obs)))
!        subsample = min(subsample,100)
!        ! Allocate which ever is smaller
!        DATAin%sub_nSWE = min(subsample,DATAin%nSWE)
!        allocate(DATAin%sub_SWEpts(DATAin%sub_nSWE))
!        ! Estimate step between draws following: ((to - from)/(length.out - 1))
!        mz = (dble(DATAin%nSWE-1)/dble((DATAin%sub_nSWE-1)))**(-1d0)
!        if (DATAin%nSWE == 1 .or. DATAin%sub_nSWE == 1) mz = 1d0 ! Guard against NaN
!        do i = 1, DATAin%sub_nSWE
!           ! Estimate the index of the next value we want
!           DATAin%sub_SWEpts(i) = DATAin%SWEpts(nint((dble(i-1)/mz)+1))
!        end do
!        print*,"SWE subsample = ",DATAin%sub_nSWE," observations of ",DATAin%nSWE
!    end if
!
!    if (DATAin%nNBE > 0) then
!        ! Estimate the maximum number of observations per data stream
!        subsample = ceiling(dble(DATAin%nodays) * subsample_fraction * (dble(DATAin%nNBE)/dble(DATAin%total_obs)))
!        subsample = min(subsample,100)
!        ! Allocate which ever is smaller
!        DATAin%sub_nNBE = min(subsample,DATAin%nNBE)
!        allocate(DATAin%sub_NBEpts(DATAin%sub_nNBE))
!        ! Estimate step between draws following: ((to - from)/(length.out - 1))
!        mz = (dble(DATAin%nNBE-1)/dble((DATAin%sub_nNBE-1)))**(-1d0)
!        if (DATAin%nNBE == 1 .or. DATAin%sub_nNBE == 1) mz = 1d0 ! Guard against NaN
!        do i = 1, DATAin%sub_nNBE
!           ! Estimate the index of the next value we want
!           DATAin%sub_NBEpts(i) = DATAin%NBEpts(nint((dble(i-1)/mz)+1))
!        end do
!        print*,"NBE subsample = ",DATAin%sub_nNBE," observations of ",DATAin%nNBE
!    end if

    ! timestep mean temperature (oC)
    DATAin%meantemp = sum((DATAin%met(2,:) + DATAin%met(3,:)) * 0.5d0) / dble(DATAin%nodays)
    ! mean SW radiation (MJ/m2/day)
    DATAin%meanrad = sum(DATAin%met(4,:)) / dble(DATAin%nodays)
    ! mean atmospheric CO2 (ppm)
    DATAin%meanco2 = sum(DATAin%met(5,:)) / dble(DATAin%nodays)
    ! mean precipitation (mm/yr)
    DATAin%meanprecip = sum(DATAin%met(7,:)*84600d0*365.25d0) / dble(DATAin%nodays)

    ! print the mean temperature and radiation variables
    write(*,*) "Mean Rad (MJ/m2/day) = ", DATAin%meanrad
    write(*,*) "Mean Temp (Celcius) = ", DATAin%meantemp
    write(*,*) "Mean Precip (mm/yr) = ", DATAin%meanprecip
    write(*,*) "Mean CO2 (ppm) = ", DATAin%meanco2
    write(*,*) "==========="
    write(*,*) "Number of timesteps = ", DATAin%nodays
    write(*,*) "Total number of obs = ", DATAin%total_obs

    ! all done
    write(*,*) "Binary input file has been successfully read by CARDAMOM"

  end subroutine read_binary_data
  !
  !------------------------------------------------------------------
  !
  subroutine read_pari_data(infile)
    use MCMCOPT, only: PI
    use MODEL_PARAMETERS, only: pars_info
    use cardamom_structures, only: DATAin

    ! subroutine call for input binary to be read and then allocates the input
    ! data to extracting the data to the correct observation and parameter types

    implicit none

    ! declare input variables
    character(350), intent(in) :: infile

    ! declare local variables
    integer :: i

    ! remind us what file we're about to access
    write(*,*) "Input file = ",trim(infile)

    ! initialise data structure and read the binary
    call read_binary_data(infile)

    ! need to allocate memory to the model output variables
    allocate(DATAin%M_LAI(DATAin%nodays),DATAin%M_GPP(DATAin%nodays) &
            ,DATAin%M_NEE(DATAin%nodays),DATAin%M_FLUXES(DATAin%nodays,DATAin%nofluxes)&
            ,DATAin%M_POOLS((DATAin%nodays+1),DATAin%nopools))

    ! force zero in states and fluxes
    DATAin%M_LAI(:) = 0d0 ; DATAin%M_GPP(:) = 0d0 ; DATAin%M_NEE(:) = 0d0
    DATAin%M_FLUXES(:,:) = 0d0 ; DATAin%M_POOLS(:,:) = 0d0

    ! alert the user
    write(*,*)"Created fields for model output"

    ! Begin allocating parameter info
    PI%npars = DATAin%nopars
    allocate(PI%parmin(PI%npars),PI%parmax(PI%npars),PI%parini(PI%npars) &
            ,PI%parfix(PI%npars),PI%parvar(PI%npars),PI%paradj(PI%npars) &
            ,PI%covariance(PI%npars,PI%npars),PI%mean_par(PI%npars) &
            ,PI%iC(PI%npars,PI%npars))

    ! force zero
    PI%parmin = 0d0 ; PI%parmax = 0d0 ; PI%parini = 0d0
    PI%parfix = 0d0 ; PI%parvar = 0d0 ; PI%paradj = 0d0
    PI%covariance = 0d0 ; PI%iC = 0d0

    ! load parameter max/min information
    call pars_info
    ! For log-normalisation procedure, no parameter can be <=0.
    ! To facilitate easy of setting parameter ranges to real values
    ! we here instead calculate the adjustment need to ensure positive only values
    where (PI%parmin <= 0d0) PI%paradj = abs(PI%parmin) + 1d0

!    ! load response surface if using the AT-DALEC model
!    if (DATAin%ID == 3 .or. DATAin%ID == 4) then
!       call load_emulator_parameters
!    end if

    ! defining initial MHMCMC stepsize and standard deviation
    PI%parvar = 1d0 ; PI%Nparvar = 0d0
    ! Covariance matrix cannot be set to zero therefore set initial value to a
    ! small positive value along to variance access
    PI%covariance = 0d0 ; PI%mean_par = 0d0 ; PI%cov = .false. ; PI%use_multivariate = .false.
    do i = 1, PI%npars
       PI%covariance(i,i) = 1d0
    end do

    ! report back to user
    write(*,*) "Created field for parameter and covariances"

  end subroutine read_pari_data
  !
  !------------------------------------------------------------------
  !
  subroutine read_options(solutions_wanted,freq_print,freq_write,outfile)
    use MCMCOPT, only: MCO, MCOUT

    ! loads required options about the MHMCMC either form hardcoded sources or
    ! from variables which were read form the command line

    implicit none

    ! declare input variables
    character(350), intent(in) :: outfile
    integer, intent(in) :: solutions_wanted, freq_print, freq_write

    ! defining hardcoded MCMC options
    MCO%append = 1
    MCO%nADAPT = 1000 ! TLS: 500 -> 1000 -> 5000 -> 10000
    MCO%fADAPT = 0.5d0
    MCO%randparini = .false.
    MCO%returnpars = .false.
    MCO%fixedpars  = .true. ! TLS: changed from .false. for testing 16/12/2019

    ! command line options

    ! how many accepted parameters to accept before completion
    if (solutions_wanted > 0) then
        MCO%nOUT = solutions_wanted
    else
        MCO%nOUT = 1000
        write(*,*)"Default MCO%nOUT value used"
    endif
    ! how frequently to print information to screen
    if (freq_print >= 0) then
        MCO%nPRINT = freq_print
    else
        MCO%nPRINT = 1000
        write(*,*)"Default MCO%nPRINT value used"
    end if
    ! how frequently to write information to file
    if (freq_write >= 0) then
        MCO%nWRITE = freq_write
    else
        MCO%nWRITE = 1000
        write(*,*)"Default MCO%nWRITE value used"
    end if

    ! Assume that sub-sampling process, if completed, will use 10 % of the
    ! simulation time therefore we want to adjust the output frequency to
    ! correct for this
    MCO%nOUT = max(1,MCO%nOUT - MCOUT%nos_iterations)

    ! construct file names
    write(MCO%outfile,fmt='(A)')trim(outfile)//"PARS"
    write(MCO%stepfile,fmt='(A)')trim(outfile)//"STEP"
    write(MCO%covfile,fmt='(A)')trim(outfile)//"COV"
    write(MCO%covifile,fmt='(A)')trim(outfile)//"COVINFO"

  end subroutine read_options
  !
  !-------------------------------------------------------------------
  !
  subroutine update_for_restart_simulation
    use MCMCOPT, only: PI, MCOUT, MCO
    use cardamom_structures, only: DATAin
    use math_functions, only: std, covariance_matrix, inverse_matrix, par2nor

    ! subroutine is responsible for loading previous parameter and step size
    ! information into the current

    implicit none

    ! local variables
    integer :: a, b, c, i, j, num_lines, status
    double precision :: dummy
    double precision,dimension(:,:), allocatable :: tmp

    ! the parameter and step files should have already been openned so
    ! read the parameter and step files to get to the end

    ! rewind to the beginning
    rewind(pfile_unit) ; rewind(sfile_unit) ; rewind(cifile_unit)

    !
    ! As this subroutine will only be called once reading each file will occur
    ! separately to improve simplicity.
    !

    !
    ! Parameter file - stored as non-normalised values
    !

    ! count the number of lines in the file..
    status = 0 ; num_lines = 0
    do
      read(pfile_unit,iostat=status) dummy
      if ( status .ne. 0 ) exit
      num_lines = num_lines + 1
    enddo
    ! Determine the number of complete parameter vectors stored. Note that the +
    ! 1 is due to the log-likelihood score being saved as well.
    num_lines = num_lines/(DATAin%nopars+1)

    ! Allocate memory to our temperary variable and the normalised parameter
    ! vector equivalent.
    allocate(tmp(num_lines,(DATAin%nopars+1)))
    ! rewind so that we can read the contents now correctly
    rewind(pfile_unit)

    ! Read the data for real
    do i = 1, num_lines
       do j = 1, (DATAin%nopars+1)
          read(pfile_unit) tmp(i,j)
       end do ! j for parameter
    end do ! i for combinations

    ! Determine the total number of iterations processed so far
    MCOUT%nos_iterations = num_lines * MCO%nWRITE
    ! Extract the final parameter set and load into the initial parameter vector
    ! for the analysis. NOTE: This parameter set will be normalised on entry
    ! into the MHMCMC subroutine
    PI%parini(1:DATAin%nopars) = tmp(num_lines,1:DATAin%nopars)

    ! free up variable for new file
    deallocate(tmp)

    !
    ! Variance file - stores output of the current parameter variance
    !

    ! count the number of remaining lines in the file..
    status = 0 ; num_lines = 0
    do
      read(sfile_unit,iostat=status) dummy
       if ( status .ne. 0. ) exit
       num_lines = num_lines + 1
    enddo

    ! Determine the number of actual stepsize vectors present.
    ! The + 1 is due to the local acceptance rate being provided too.
    num_lines = num_lines/(DATAin%nopars+1)
    ! allocate memory
    allocate(tmp(num_lines,(DATAin%nopars+1)))
    ! rewind, for actual reading
    rewind(sfile_unit)

    ! now read the data for real
    do i = 1, num_lines
       do j = 1, (DATAin%nopars+1)
          read(sfile_unit) tmp(i,j)
       end do ! j for parameter
    end do ! i for combinations
    ! Save the current acceptance_rate
    MCOUT%acceptance_rate = MCOUT%nos_iterations * MCO%nWRITE

    ! tidy up for the next file
    deallocate(tmp)

    !
    ! Covariance matrix file
    !

    ! The covariance matrix may contain either 1 or 2 complete matrices. We will
    ! just want to the latest one.

    ! count the number of remaining lines in the file..
    status = 0 ; num_lines = 1
    do
       read(cfile_unit,iostat=status,rec = num_lines) dummy
       if ( status .ne. 0. ) exit
       num_lines = num_lines + 1
    enddo

    ! Determine whether there is 1 or more matrice here
    if ((num_lines/DATAin%nopars)/DATAin%nopars == 1) then
        ! the size of the file is consistent with a single matrix having been
        ! saved
        a = 1
    else if ((num_lines/DATAin%nopars)/DATAin%nopars == 2) then
        !
        a = 2
    else
        ! something has gone wrong - best stop
        print*,"Error reading COV file"
        print*,"DATAin%nopars = ",DATAin%nopars,"COV length = ",num_lines * DATAin%nopars
        stop
    endif

    ! now read the data for real
    c = 1
    do b = 1, a
       do i = 1, DATAin%nopars
          do j = 1, DATAin%nopars
             read(cfile_unit, rec = c) PI%covariance(i,j)
             c = c + 1
          end do ! j for parameter
       end do ! i for combinations
    end do

    ! extract current variance information
    do i = 1, PI%npars
       PI%parvar(i) = PI%covariance(i,i)
    end do
    ! estimate status of the inverse covariance matrix
    call inverse_matrix( PI%npars, PI%covariance, PI%iC )

    !
    ! Covariance information file
    !

    ! The number of parameters on which the covariance matrix is based must be
    ! known to allow for correct updating. Similarly the mean normalised
    ! parameter values are also needed

    ! count the number of remaining lines in the file..
    status = 0 ; num_lines = 0
    do
      read(cifile_unit,iostat=status) dummy
       if ( status .ne. 0. ) exit
       num_lines = num_lines + 1
    enddo

    ! how many parameter vectors have been output. Note the + 1 is accounting
    ! for the number of samples underlying the mean
    num_lines = num_lines / (PI%npars + 1)
    ! allocate memory
    allocate(tmp(num_lines,(DATAin%nopars+1)))
    ! rewind, for actual reading
    rewind(cifile_unit)

    ! now read the data for real
    do i = 1, num_lines
       do j = 1, (DATAin%nopars+1)
          read(cifile_unit) tmp(i,j)
       end do ! j for parameter
    end do ! i for combinations
    ! Store the most recent step size, which corresponds with the saved
    ! parmeters (above) and covariance matrix (below)
    PI%mean_par = tmp(num_lines,1:DATAin%nopars)
    PI%Nparvar = tmp(num_lines,DATAin%nopars+1)

    return

  end subroutine update_for_restart_simulation
  !
  !------------------------------------------------------------------
  !
  subroutine write_covariance_matrix(covariance,npars,initial_cov)

    ! subroutine writes MCMC accepted parameters and step values to binary files

    implicit none

    ! arguments
    logical, intent(in) :: initial_cov
    integer, intent(in) :: npars
    double precision, dimension(npars,npars), intent(in) :: covariance

    ! declare local variables
    integer :: i,j,irec

    ! If we have already written the initial covariance matrix we want to keep
    ! over-writing the current matrix. We do this to avoid large files form
    ! writing out multiple covariance matrices
    if (.not.initial_cov) then
        irec = npars * npars
    else
        irec = 0
    end if

    ! write out the file. Its binary format has already been determined at the
    ! openning of the file

    do i = 1, npars
       do j = 1, npars
          irec = irec + 1
          write(cfile_unit, rec = irec) covariance(i,j)
       end do
    end do

    return

  end subroutine write_covariance_matrix
  !
  !------------------------------------------------------------------
  !
  subroutine write_covariance_info(mean_pars,nsample,npars)

    ! subroutine writes MCMC accepted parameters and step values to binary files

    implicit none

    ! arguments
    integer, intent(in) :: npars
    double precision, intent(in) :: nsample
    double precision, dimension(npars), intent(in) :: mean_pars

    ! declare local variables
    integer :: i,j

    ! write out the file. Its binary format has already been determined at the
    ! openning of the file

    do i = 1, npars
       write(cifile_unit) mean_pars(i)
    end do

    write(cifile_unit) nsample

    return

  end subroutine write_covariance_info
  !
  !------------------------------------------------------------------
  !
  subroutine write_variances(variance,npars,accept_rate)

    ! subroutine writes parameter variance for corresponding parameter values

    implicit none

    ! declare input variables
    integer, intent(in) :: npars
    double precision, dimension(npars), intent(in) :: variance
    double precision, intent(in) :: accept_rate ! local acceptance rate

    ! declare local variables
    integer :: n

    ! write out the file. Its binary format has already been determined at the
    ! openning of the file

    do n = 1, npars
       write(sfile_unit) variance(n)
    end do

    ! we will need to know the current acceptance rate for restarts
    write(sfile_unit) accept_rate

    return

  end subroutine write_variances
  !
  !------------------------------------------------------------------
  !
  subroutine write_parameters(pars,prob,npars)

    ! subroutine writes parameter values to binary file`

    implicit none

    ! declare input variables
    integer, intent(in) :: npars
    double precision, dimension(npars), intent(in) :: pars
    double precision, intent(in) :: prob

    ! declare local variables
    integer :: n

    ! write out the file. Its binary format has already been determined at the
    ! openning of the file

    do n = 1, npars
       write(pfile_unit) pars(n)
    end do

    ! now add the probability
    write(pfile_unit) prob

    ! close will occur at the end of the MCMC

    ! return back
    return

  end subroutine write_parameters
  !
  !--------------------------------------------------------------------
  !
  subroutine write_mcmc_output(variance, accept_rate, &
                               covariance, mean_pars, nsample, &
                               pars, prob, npars, dump_now)
    use cardamom_structures, only: io_space

    ! Arguments
    integer, intent(in) :: npars
    double precision, dimension(npars,npars), intent(in) :: covariance
    double precision, dimension(npars), intent(in) :: mean_pars, &
                                                       variance, &
                                                           pars
    double precision, intent(in) :: nsample, accept_rate, prob
    logical, intent(in) :: dump_now

    ! Local variables
    integer :: i

    ! Increment buffer
    io_space%io_buffer_count = io_space%io_buffer_count + 1

    ! Store information in buffer for later writing
    io_space%variance_buffer(1:npars,io_space%io_buffer_count) = variance
    io_space%mean_pars_buffer(1:npars,io_space%io_buffer_count) = mean_pars
    io_space%pars_buffer(1:npars,io_space%io_buffer_count) = pars
    io_space%prob_buffer(io_space%io_buffer_count) = prob
    io_space%nsample_buffer(io_space%io_buffer_count) = nsample
    io_space%accept_rate_buffer(io_space%io_buffer_count) = accept_rate

    ! Are we storing information in buffer or writing to file?
    if (io_space%io_buffer_count == io_space%io_buffer .or. dump_now) then

        ! Then we are writing out to file

        ! Only write the most current covariance matrix as this would be an overwrite anyway
        call write_covariance_matrix(covariance,npars,.false.)
        ! Everything else loop through the buffered output to write out
        do i = 1, io_space%io_buffer_count
           call write_covariance_info(io_space%mean_pars_buffer(:,i),io_space%nsample_buffer(i),npars)
           call write_variances(io_space%variance_buffer(:,i),npars,io_space%accept_rate_buffer(i))
           call write_parameters(io_space%pars_buffer(:,i),io_space%prob_buffer(i),npars)
        end do

        ! Reset buffer increment
        io_space%io_buffer_count = 0

    endif

  end subroutine write_mcmc_output
  !
  !--------------------------------------------------------------------
  !
end module cardamom_io
