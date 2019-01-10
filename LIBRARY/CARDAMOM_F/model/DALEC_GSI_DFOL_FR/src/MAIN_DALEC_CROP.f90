program MAIN_DALEC

 
  use DALEC_CROP_MET_VARIABLES,     only: met_data
  use DALEC_CROP_LEAF_MASS,         only: leaf_mass
  use DALEC_CROP_DEV_VARIABLES,     only: DS_shoot,DS_root,fol_frac,stem_frac,root_frac,DS_LRLV, &
                                          LRLV,DS_LRRT,LRRT
  
  use DALEC_CROP_IO,                only: read_met_data, read_veg_parameters, output_data, read_leaf
  

  use CARBON_MODEL_CROP_MOD,        only: CARBON_MODEL_CROP
  
  implicit none


 
  
  
  !--Local variables
  integer :: num_days, start, finish, num_pools, num_fluxes, num_met, num_pars, sow_day, harvest_day, plough_day, pft, k, run   

  double precision :: lat, stock_seed_labile


  ! Allocate lai array
  double precision, allocatable, dimension(:) :: lai

  ! Allocate GPP
  double precision, allocatable, dimension(:) :: GPP

  ! Allocate NEE
  double precision, allocatable, dimension(:) :: NEE


  ! Allocate DS
  double precision, allocatable, dimension(:) :: DS_array
  
  
  ! Allocate doy in decimal
  double precision, allocatable, dimension(:) :: deltat

  ! Allocate pars in decimal
  double precision, allocatable, dimension(:) :: pars

  ! Allocate pools and fluxes arrays
  double precision, allocatable, dimension(:,:) :: pools
  double precision, allocatable, dimension(:,:) :: fluxes

  
  !------------------------------------------------------------
  ! some local scaling variables
  start = 1
  finish = 730
  num_days = 730
  lat = 55.880653
  sow_day = 273
  harvest_day = 237
  plough_day = 730
  pft=1 ! could be removed
  
  num_pools = 9
  num_fluxes = 21
  num_met = 6
  num_pars = 38
  
  stock_seed_labile = 9d0
  
  !------------------------------------------------------------
  
  
  
  ! Allocate parameters
  if ( .not. allocated( pars ) )  allocate ( pars( num_pars ) )

  ! Read veg parameter data from *csv file (this part could be hard coded)
  ! Actually if you go to DALEC_CROP_IO.f90 you will see that I am using some hard-coded parameters, but we could update these with the ones from your optimisation.
  call read_veg_parameters(pars, num_pars, sow_day, harvest_day, plough_day)
   
  
  ! Read met data *csv into met_data array (doy,min_t,max_t,rad)
  call read_met_data(num_days)

  ! This is currently not in use but it was where I set it to read leaf parameters (i.e. N content) based on the observations 
  !call read_leaf(num_days)
 
  
  ! set up data for CARBON_MODEL_CROP
  ! Allocate LAI array
  if ( .not. allocated( lai ) )  allocate ( lai( num_days ) ) ! Allocate LAI
  ! Allocate GPP array
  if ( .not. allocated( GPP ) )  allocate ( GPP( num_days ) ) ! Allocate LAI 
  ! Allocate NEE array
  if ( .not. allocated( NEE ) )  allocate ( NEE( num_days ) ) ! Allocate LAI
  ! Allocate deltat array
  if ( .not. allocated( deltat ) )  allocate ( deltat( num_days ) ) ! Allocate LAI

  ! Allocate pools and fluxes
  if ( .not. allocated( POOLS ) )  allocate ( POOLS( num_days, num_pools ) ) ! Allocate 8 pools
  if ( .not. allocated( FLUXES ) )  allocate ( FLUXES( num_days, num_fluxes ) ) ! Allocate 16 fluxes

  !!$  ! Allocate DS- this is to save the DS into an array
  !!$  if ( .not. allocated( DS_array ) )  allocate ( DS_array( num_days) ) ! Allocate 8 pools

  do k = 1, num_days; 
     deltat(k) = 1d0
  end do
     
  ! THIS IS THE MAIN PART WHERE DALEC-CROP IS CALLED AND THE OUTPUTS ARE GENERATED.
  ! THIS COULD BE PUT INTO A LOOP FOR ENSEMBLE RUNS
  
     ! Call the modelling part
     call CARBON_MODEL_CROP(start,finish,met_data,pars,deltat,num_days,lat,lai,NEE,FLUXES,POOLS &
                              ,pft,num_pars,num_met,num_pools,num_fluxes,GPP,stock_seed_labile)
     
     call output_data(FLUXES,POOLS, pars)


end program MAIN_DALEC
  
