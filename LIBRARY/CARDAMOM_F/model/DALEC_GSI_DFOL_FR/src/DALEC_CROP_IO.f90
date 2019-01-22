module DALEC_CROP_IO

  implicit none

  private

  public :: read_met_data, read_veg_parameters, output_data
  
  
contains
!
!---------------------------------------------------------------------------------------------------------
!  
  subroutine read_veg_parameters(pars, num_pars)

    use CARBON_MODEL_CROP_MOD

    implicit none

    integer, intent(in) :: num_pars
    double precision, dimension(:), intent(inout) :: pars
    
    
    ! Local
    character(LEN=200)               :: input_veg, head, dummy
    integer                          :: veg_file= 101, i

    ! VEG DATA DIRECTORY
    input_veg = ('src/INPUTS/'//'pars.csv')

    open (unit= veg_file, file=trim(input_veg), status='old')   ! veg file open


    read(veg_file,*) pars(1) 
    read(veg_file,*) pars(2)
    read(veg_file,*) pars(3)
    read(veg_file,*) pars(4)
    read(veg_file,*) pars(5)
    read(veg_file,*) pars(6)
    read(veg_file,*) pars(7)
    read(veg_file,*) pars(8)
    read(veg_file,*) pars(9)
    read(veg_file,*) pars(10)
    read(veg_file,*) pars(11)
    read(veg_file,*) pars(12) ! sow day
    read(veg_file,*) pars(13)
    read(veg_file,*) pars(14)
    read(veg_file,*) pars(15) ! harvest day
    read(veg_file,*) pars(16) ! plough day
    read(veg_file,*) pars(17)
    read(veg_file,*) pars(18)
    read(veg_file,*) pars(19)
    read(veg_file,*) pars(20)
    read(veg_file,*) pars(21)
    read(veg_file,*) pars(22)
    read(veg_file,*) pars(23)
    read(veg_file,*) pars(24)
    read(veg_file,*) pars(25)
    read(veg_file,*) pars(26)
    read(veg_file,*) pars(27)
    read(veg_file,*) pars(28)
    read(veg_file,*) pars(29)
    read(veg_file,*) pars(30)
    read(veg_file,*) pars(31)
    read(veg_file,*) pars(32)
    read(veg_file,*) pars(33)
    read(veg_file,*) pars(34)
    read(veg_file,*) pars(35)
    
    close(veg_file)
    
  end subroutine read_veg_parameters

!
!---------------------------------------------------------------------------------------------------------
!    
  subroutine read_met_data(num_days)

    use DALEC_CROP_MET_VARIABLES

    implicit none

    ! Arguments
    integer, intent(in) :: num_days

    ! Local
    character(LEN=200)               :: input_met, head
    integer                          :: met_file= 100, j
        
    if (.not. allocated (met_data)) allocate (met_data (6, num_days)) ! for output LAI
    
    ! MET DRIVER DIRECTORY
    input_met = ('src/INPUTS/ES_met_2017_2018_v3.csv')
          
    open (unit= met_file, file=trim(input_met), status='old')   ! Local Met

    read(met_file,*) head

    ! julian day, min_t, max_temp, radiation, co2, DOY
    do j = 1,730
       read(met_file,*) met_data(1,j), met_data(2,j), met_data(3,j), met_data(4,j), met_data(5,j), met_data(6,j)
    end do
!!$    read(met_file,*) (met_data(1,j), met_data(2,j), met_data(3,j), met_data(4,j), met_data(5,j), met_data(6,j),  met_data(7,j), &
!!$         met_data(8,j),  met_data(9,j),  met_data(10,j),  met_data(11,j), j = 1,num_days)

    !print*, (met_data(11,:))
    
  end subroutine read_met_data

!
!---------------------------------------------------------------------------------------------------------
!  
  
 subroutine output_data(FLUXES,POOLS, pars, lai)

    use DALEC_CROP_MET_VARIABLES

    implicit none

    ! Arguments   
    double precision, dimension(:,:), intent(in) :: POOLS ! vector of ecosystem pools
    double precision, dimension(:,:), intent(in) :: FLUXES  ! vector of ecosystem fluxes
    double precision, dimension(:), intent(in) :: pars, lai  ! vector of pars
    
    ! Local
    integer                          :: daily_flux_file=901,  &
                                        daily_pools_file=902, &
                                        j, p

    character(LEN=200)               :: output_dir
    character(len=1)                 :: filename

    
    
    !write (filename, "(I1)") treatment_num
    
    output_dir = 'src/OUTPUTS/'   
    
    ! All the 16 fluxes could be included here
    ! write FLUXES -----------------------
    open(unit=daily_flux_file, file = (trim(output_dir)//'FLUXES.csv'), status="unknown")
    write(daily_flux_file,'(1(A10,","))') "gpp_acm"
    write(daily_flux_file,'(1(F16.4,","))') (FLUXES(j,1), j = 1,730)
    ! -------------------------------------

    ! All the 8 pools could be added here
    ! write POOLS ------------------------
    open(unit=daily_pools_file, file = (trim(output_dir)//'POOLS.csv'), status="unknown")
    write(daily_pools_file,'(4(A10,","))') "stock_foliage","stock_stem","stock_root", "stock_storage" !NOT CORRECT TITLES AT MOMENT
    write(daily_pools_file,'(4(F16.4,","))')  (POOLS(p,2), POOLS(p,4), POOLS(p,3), POOLS(p,9),  p = 1,730)

!!$    open(909, file = (trim(output_dir)//'LAI.csv'), status="unknown")
!!$    write(909,'(1(F16.4))') (LAI(p),  p = 1,730)
!!$    
  end subroutine output_data
!
!---------------------------------------------------------------------------------------------------------
!  

  
end module DALEC_CROP_IO
    
