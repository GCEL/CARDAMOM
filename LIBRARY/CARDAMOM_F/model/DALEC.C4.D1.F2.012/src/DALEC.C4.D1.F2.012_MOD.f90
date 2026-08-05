module carbon_model_memory
  implicit none
  public

  ! for consisteny between requirements of different models
  integer, parameter :: nos_root_layers = 2, nos_soil_layers = nos_root_layers + 1
  double precision, parameter :: pi = 3.1415927d0, &
                         deg_to_rad = 0.01745329d0   ! pi/180d0

  type model_working_variables

    double precision, dimension(nos_soil_layers) :: soil_frac_clay, soil_frac_sand

    ! arrays for the emulator, just so we load them once and that is it cos they be
    ! massive
    integer ::    dim_1, & ! dimension 1 of response surface
                  dim_2, & ! dimension 2 of response surface
              nos_trees, & ! number of trees in randomForest
             nos_inputs    ! number of driver inputs

    double precision, allocatable, dimension(:,:) ::     leftDaughter, & ! left daughter for forest
                                                        rightDaughter, & ! right daughter for forets
                                                           nodestatus, & ! nodestatus for forests
                                                           xbestsplit, & ! for forest
                                                             nodepred, & ! prediction value for each tree
                                                              bestvar    ! for randomForests
    ! Module level variables for ACM
    double precision :: ci ! Internal CO2 concentration (ppm)
  end type
  type(model_working_variables), allocatable, dimension(:):: mVs

contains
  !
  !--------------------------------------------------------------------
  !
end module carbon_model_memory
