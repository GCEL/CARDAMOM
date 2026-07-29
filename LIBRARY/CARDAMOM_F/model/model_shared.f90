module model_shared
  use samplers_shared, only: PARINFO  
  use cardamom_structures, only: DATA_TYPE


  public


  type(PARINFO)  :: PI  ! should not be writted to except by pars_info() !  
                                ! type definition differs by model and lines in _PARS.f90 file
                                ! If we want parallel runs 
                                ! now Contains read-only description of the model parameters only


  contains

  ! split from read_pari_data
  subroutine initialize_parinfo()
    use MODEL_PARAMETERS, only:  pars_info
    implicit none 
    integer:: i

    ! load parameter max/min information, npars 
    call pars_info(PI)

  ! Begin allocating parameter info
    if (.not. allocated(PI%parmin)) then 
         allocate(PI%parmin(PI%npars))
         PI%parmin = 0d0
    endif 
    if (.not. allocated(PI%parmax)) then 
        allocate(PI%parmax(PI%npars))
        PI%parmax = 0d0 
    endif
    if (.not. allocated(PI%parfix)) then 
        allocate(PI%parfix(PI%npars))
        PI%parfix = .false.
    endif
    if (.not. allocated(PI%paradj)) then 
        allocate(PI%paradj(PI%npars))
        PI%paradj = 0d0 
    endif

    ! force zero
    !PI%parini = 0d0  ! no longer exist-use DATAIN directly to -> pars in MCOUT
    !PI%parfix = 0d0; 
    ! PI%parvar = 0d0;  ! -> stats
    PI%paradj = 0d0  ! still exists, calcualted here
    !PI%covariance = 0d0; PI%iC = 0d0  ! -> stats TODO split some of these to initialize_stats
    
    ! For log-normalisation procedure, no parameter can be <= 0.
    ! To facilitate easy of setting parameter ranges to real values
    ! we here instead calculate the adjustment need to ensure positive only values
    where (PI%parmin <= 0d0) PI%paradj = abs(PI%parmin) + 1d0

    ! defining initial MHMCMC stepsize and standard deviation
    ! PI%parvar = 1d0; PI%Nparvar = 0d0
    ! Covariance matrix cannot be set to zero therefore set initial value to a
    ! small positive value along to variance access
    ! PI%covariance = 0d0; PI%meanpar = 0d0; PI%cov = .false. ; PI%use_multivariate = .false.
    !do i = 1, PI%npars
    !    PI%covariance(i, i) = 1d0
    !end do
    ! report back to user
    !write(*,*) "Created field for parameter and covariances"
   

  end subroutine

    subroutine initialize_carbon_model(n_chains)
    !! prepare N model_working_variables type objects to hold seperate sets of persistent values for
    !! each independent parallel chain.
    !! Must have DATAin filled first.
    use cardamom_structures, only: DATAin
    use CARBON_MODEL_MOD, only: mVs, initialize_mv
    integer, intent(in), optional:: n_chains
    integer:: n_chains_
    integer:: i
    if (present(n_chains)) then
      n_chains_ = n_chains
    else
      n_chains_ = 1
    endif
    allocate(mVs(n_chains))
    do i = 1, n_chains_
        call initialize_mv(Mvs(i), DATAin%nodays, DATAin%nomet, DATAin%nopars, DATAin%deltat, &
                           DATAin%soil_frac_sand, DATAin%soil_frac_clay, DATAin%met, DATAin%lat)

    end do
    end subroutine

    subroutine destroy_carbon_model()
    !! deallocate members of model_working_variables struct(s) Mvs
    use CARBON_MODEL_MOD, only: mVs
    integer:: n_chains_
    if (allocated(mVs)) deallocate(mVs)
    end subroutine
  
end module
