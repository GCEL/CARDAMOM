module DALEC_CROP_IO

  implicit none

  private

  public :: read_met_data, read_veg_parameters, read_dev_data, output_data, read_leaf
  
  
contains
!
!---------------------------------------------------------------------------------------------------------
!  
  subroutine read_veg_parameters(pars, num_pars, sow_day, harvest_day, plough_day)

    use CARBON_MODEL_CROP_MOD

    implicit none

    integer, intent(in) :: num_pars, sow_day, harvest_day, plough_day
    double precision, dimension(:), intent(inout) :: pars
    
    
    ! Local
    character(LEN=200)               :: input_veg, head, dummy
    integer                          :: veg_file= 101, i

    ! VEG DATA DIRECTORY
    input_veg = ('src/STANDALONE/INPUTS/'//'WW_crops_veg.csv')

    open (unit= veg_file, file=trim(input_veg), status='old')   ! veg file open

    ! Skip first few lines
    do i = 1,4
       read(veg_file,*) dummy
    end do

    
    read(veg_file,*) dummy !, foliar_n
    read(veg_file,*) dummy !stem conductivity
    read(veg_file,*) dummy !min LWP
    read(veg_file,*) dummy !stomatal efficiency
    read(veg_file,*) dummy !leaf capacitance
    read(veg_file,*) dummy !LAT
    read(veg_file,*) dummy !detailed layer by layer output
    read(veg_file,*) dummy !characteristic dimensions of leaf
    read(veg_file,*) dummy !root resistivity
    read(veg_file,*) dummy !height of tower with sensors
    read(veg_file,*) dummy !tree conductivity
    read(veg_file,*) dummy !Rate coefficient Vcmax
    read(veg_file,*) dummy !Rate coefficient Jmax
    read(veg_file,*) dummy, pars(17)
    read(veg_file,*) dummy !max root depth
    read(veg_file,*) dummy !Root biomass to reach 50% of max depth
    read(veg_file,*) dummy, pars(1) 
    read(veg_file,*) dummy, pars(2)
    read(veg_file,*) dummy, pars(5)
    read(veg_file,*) dummy, pars(6)
    read(veg_file,*) dummy, pars(10)
    read(veg_file,*) dummy, pars(9)
    read(veg_file,*) dummy, pars(34)
    read(veg_file,*) dummy, pars(13)
    read(veg_file,*) dummy, pars(11)
    read(veg_file,*) dummy, pars(19)
    read(veg_file,*) dummy, pars(21)
    read(veg_file,*) dummy, pars(20)
    read(veg_file,*) dummy, pars(22)
    read(veg_file,*) dummy, pars(23)
    read(veg_file,*) dummy, pars(18)
    read(veg_file,*) dummy, pars(24)
    read(veg_file,*) dummy !THROUGHFALL, fraction of precip that penetrates canopy
    read(veg_file,*) dummy !max storage, water retained in canopy
    read(veg_file,*) dummy, pars(14)
    read(veg_file,*) dummy !Fraction of leaf biomass remaining after harvest 
    read(veg_file,*) dummy !Fraction of stem biomass remaining after harvest
    read(veg_file,*) dummy, pars(3)
    read(veg_file,*) dummy, pars(4)
    read(veg_file,*) dummy, pars(26)
    read(veg_file,*) dummy, pars(28)
    read(veg_file,*) dummy, pars(27)
    read(veg_file,*) dummy, pars(29)
    read(veg_file,*) dummy, pars(31)
    read(veg_file,*) dummy, pars(30)
    read(veg_file,*) dummy, pars(8)
    read(veg_file,*) dummy !critical lai
    read(veg_file,*) dummy, pars(7)
    read(veg_file,*) dummy, pars(32)
    read(veg_file,*) dummy, pars(33)

    

    rewind(veg_file)
    
    pars(12) = sow_day    
    pars(15) = harvest_day
    pars(16) = plough_day
        
    !pars(25) = stock_storage_organ


    ! LUKE PARAMETERS
    pars(1:38)= (/5.396430e-05,  & !1
                  3.591521e-01,  & !2 Frac GPP
                  (3.538961e-02)*0.7,  & !3 DR pre
                  4.559480e-02*1.7,  & !4 DR post
                  8.596312e-02,  & !5
                  9.389972e-03,  & !6
                  4.176272e-04,  & !7 Max rate of foliar turnover due to self-shading
                  2.620799e+01,  & !8
                  1.009997e-01,  & !9
                  1.302387e-05,  & !10
                  3.5e-01,  & !11 foliar N (was 3.038757e-01)
                  3.172564e+02,  & !12
                  2.069155e-01,  & !13
                  1.032496e+02*1.35,  & !14 PHUem
                  3.111215e+02,  & !15
                  3.771163e+02,  & !16
                  2.419578e+01*0.857,  & !17 LMA
                  3.317560e+00,  & !18
                  0e+00,         & !19 foliar c  2.757770e+00
                  0e+00,         & !20 roots 3.507375e+00
                  0e+00,         & !21 stem 2.318463e+00
                  2.184431e+00,  & !22
                  1.289216e+04,  & !23
                  1.076274e+00,  & !24
                  0e+00,  & !25
                  2.723468e+02,  & !26 min temp
                  2.999132e+02,  & !27 max temp
                  2.842096e+02,  & !28 opt temp
                  2.691015e+02,  & !29
                  2.912893e+02,  & !30
                  2.767106e+02,  & !31
                  8.791524e+00,  & !32 Critical photoperiod for development
                  1.283520e-01,  & !33 Photoperiod sensitivity
                  9.999767e-02,  & !34
                  4.610077e-02,  & !35
                  3.437755e+01,  & !36
                  4.170504e-01,  & !37
                  3.972297e-01/)   !38


    ! Change foliar
    
    
  end subroutine read_veg_parameters

!
!---------------------------------------------------------------------------------------------------------
!  
  subroutine read_dev_data

    use DALEC_CROP_DEV_VARIABLES
    
    ! Read the contents of the development.csv file !

    implicit none
    
    ! Local
    character(LEN=200)               :: input_dev, head
    integer                          :: columns, i, rows, dev_file = 400

    
    ! open Winter Wheat crop development file
    open(unit = dev_file, file= 'src/STANDALONE/INPUTS/WW_crops_development.csv', status='old')

    ! maybe need to include stock_seed_labile here

    read(unit=dev_file,fmt=*) head  
    
    ! read in C partitioning/fraction data and corresponding developmental stages (DS)
    ! shoot
    read(unit=dev_file,fmt=*) head
    read(unit=dev_file,fmt=*) rows, columns
    if (.not. allocated(DS_shoot)) allocate( DS_shoot(rows) , fol_frac(rows) , stem_frac(rows)  )  
    do i = 1 , rows
       read(unit=dev_file,fmt=*) DS_shoot(i), fol_frac(i), stem_frac(i)        
    enddo
     
    ! root
    read(unit=dev_file,fmt=*) head
    read(unit=dev_file,fmt=*) rows , columns
    if (.not. allocated(DS_root)) allocate( DS_root(rows) , root_frac(rows) )
    do i = 1 , rows
       read(unit=dev_file,fmt=*) DS_root(i), root_frac(i)
    enddo


    ! loss rates of leaves and roots
    ! leaves
    read(unit=dev_file,fmt=*) head
    read(unit=dev_file,fmt=*) rows , columns
    if (.not. allocated(DS_LRLV)) allocate( DS_LRLV(rows) , LRLV(rows) )
    do i = 1 , rows
      read(unit=dev_file,fmt=*) DS_LRLV(i), LRLV(i)
    enddo

    ! roots
    read(unit=dev_file,fmt=*) head
    read(unit=dev_file,fmt=*) rows , columns
    if (.not. allocated(DS_LRRT)) allocate( DS_LRRT(rows) , LRRT(rows) )
    do i = 1 , rows
      read(unit=dev_file,fmt=*) DS_LRRT(i), LRRT(i)
    enddo

    ! developmental rate as a function of temperature (not needed when calculated through modified Wang&Engel model,
    ! which is the current set-up)
    !preanthesis (before flowering)
    read(unit=dev_file,fmt=*) head
    read(unit=dev_file,fmt=*) rows , columns
    if (.not. allocated(DR_T_PRA)) allocate( DR_T_PRA(rows) , DRAO_T_PRA(rows) )
    do i = 1 , rows
      read(unit=dev_file,fmt=*) DR_T_PRA(i), DRAO_T_PRA(i)
    enddo

    ! postanthesis
    read(unit=dev_file,fmt=*) head
    read(unit=dev_file,fmt=*) rows , columns
    if (.not. allocated(DR_T_POA)) allocate( DR_T_POA(rows) , DRAO_T_POA(rows) )
    do i = 1 , rows
      read(unit=dev_file,fmt=*) DR_T_POA(i), DRAO_T_POA(i)
    enddo
 
    ! photoperiod (daylength effect on development)
    read(unit=dev_file,fmt=*) head
    read(unit=dev_file,fmt=*) rows , columns
    if (.not. allocated(DR_P)) allocate( DR_P(rows) , DRAO_P(rows) )
    do i = 1 , rows
      read(unit=dev_file,fmt=*) DR_P(i), DRAO_P(i)
    enddo

    ! LCA ratios (if LCA to be dynamic, but currently not used)
    read(unit=dev_file,fmt=*) head
    read(unit=dev_file,fmt=*) rows , columns
    if (.not. allocated(LCA_DS)) allocate( LCA_DS(rows)  , LCA_ratio(rows) )
    do i = 1 , rows
      read(unit=dev_file,fmt=*) LCA_DS(i), LCA_ratio(i)
    enddo
  
    rewind(dev_file)
        
    
  end subroutine read_dev_data  
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
    input_met = ('src/STANDALONE/INPUTS/ES_met_2017_2018_v2.csv')
          
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

 subroutine read_leaf(num_days)

    use DALEC_CROP_LEAF_MASS

    implicit none

    ! Arguments
    integer, intent(in) :: num_days

    ! Local
    character(LEN=200)               :: input_leaf, head
    integer                          :: leaf_file= 109, m
    real                             :: dummy
        
    if (.not. allocated (leaf_mass)) allocate (leaf_mass (5, num_days)) ! for output LAI
    
    ! MET DRIVER DIRECTORY
    input_leaf = ('src/STANDALONE/INPUTS/leaf_mass_v1.csv')
          
    open (unit= leaf_file, file=trim(input_leaf), status='old')   ! Local Met

    read(leaf_file,*) head

    ! julian day, min_t, max_temp, radiation, co2, DOY
    do m = 1,num_days
       read(leaf_file,*) dummy, leaf_mass(1,m), leaf_mass(2,m), leaf_mass(3,m), leaf_mass(4,m), leaf_mass(5,m)
    end do

   
  
    
  end subroutine read_leaf
!
!---------------------------------------------------------------------------------------------------------
!  
  
 subroutine output_data(FLUXES,POOLS, pars)

    use DALEC_CROP_MET_VARIABLES

    implicit none

    ! Arguments   
    double precision, dimension(:,:), intent(in) :: POOLS ! vector of ecosystem pools
    double precision, dimension(:,:), intent(in) :: FLUXES  ! vector of ecosystem fluxes
    double precision, dimension(:), intent(in) :: pars  ! vector of pars
        
    ! Local
    integer                          :: daily_flux_file=901,  &
                                        daily_pools_file=902, &
                                        j, p

    character(LEN=200)               :: output_dir
    character(len=1)                 :: filename

    
    
    !write (filename, "(I1)") treatment_num
    
    output_dir = 'src/STANDALONE/OUTPUTS/'   
    
    ! All the 16 fluxes could be included here
    ! write FLUXES -----------------------
    open(unit=daily_flux_file, file = (trim(output_dir)//'FLUXES.csv'), status="unknown")
    write(daily_flux_file,'(1(A10,","))') "gpp_acm"
    write(daily_flux_file,'(1(F16.4,","))') (FLUXES(j,1), j = 1,730)
    ! -------------------------------------

    ! All the 8 pools could be added here
    ! write POOLS ------------------------
    open(unit=daily_pools_file, file = (trim(output_dir)//'POOLS.csv'), status="unknown")
    write(daily_pools_file,'(5(A10,","))') "LAI","stock_foliage","stock_stem","stock_root", "stock_storage"
    write(daily_pools_file,'(5(F16.4,","))') (POOLS(p,2)/pars(17), POOLS(p,2), POOLS(p,4), POOLS(p,3), POOLS(p,9),  p = 1,730)

    
  end subroutine output_data
!
!---------------------------------------------------------------------------------------------------------
!  

  
end module DALEC_CROP_IO
    
