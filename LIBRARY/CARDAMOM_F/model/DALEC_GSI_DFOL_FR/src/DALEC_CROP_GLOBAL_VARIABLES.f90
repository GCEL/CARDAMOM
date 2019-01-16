
!
!---------------------------------------------------------------------------------
!

module DALEC_CROP_MET_VARIABLES

  implicit none

  ! met-drivers used/updated each timestep..
  real :: atmos_press, & ! Surface atmospheric pressure (Pa)
          at_CO2,      & ! Ambient CO2 concentration (ppm)
          min_t,       & ! minimum daily temperature
          max_t,       & ! AR:  maximum daily temperature
          run_day,     &
          rad

  double precision, allocatable, dimension(:,:) :: met_data
                                           

  save

end module DALEC_CROP_MET_VARIABLES


!
!---------------------------------------------------------------------------------
!

module DALEC_CROP_LEAF_MASS

  implicit none
  
  double precision, allocatable, dimension(:,:) :: leaf_mass
                                           
  save

end module DALEC_CROP_LEAF_MASS


!
!---------------------------------------------------------------------------------
!

module DALEC_CROP_DEV_VARIABLES

  implicit none

  double precision, allocatable, dimension(:) ::              DS_shoot, & !
                                                               DS_root, & !
                                                              fol_frac, & !
                                                             stem_frac, & !
                                                             root_frac, & !
                                                               DS_LRLV, & ! 
                                                                  LRLV, & !
                                                               DS_LRRT, & !
                                                              DR_T_PRA, & !
                                                            DRAO_T_PRA, & !
                                                              DR_T_POA, & !
                                                            DRAO_T_POA, & !
                                                                  DR_P, & !
                                                                DRAO_P, & !
                                                                LCA_DS, & !
                                                             LCA_ratio, & !
                                                                  LRRT    !

  !double precision :: stock_seed_labile

  save

  
end module DALEC_CROP_DEV_VARIABLES


!
!---------------------------------------------------------------------------------
!

module DALEC_CROP_SOIL_VARIABLES

  implicit none

  double precision :: resp_rate_temp_coeff= 0.0693

  save

  
end module DALEC_CROP_SOIL_VARIABLES
