module carbon_model_memory
  implicit none
  public

  ! useful technical parameters
  double precision, parameter :: vsmall = tiny(0d0)*1d3 & ! *1d3 to add a little breathing room
                                ,vlarge = huge(0d0)

  integer, parameter :: nos_root_layers = 2, nos_soil_layers = nos_root_layers + 1
  double precision, parameter :: pi = 3.1415927d0, &
                         deg_to_rad = 0.01745329d0   ! pi/180d0
  ! timing parameters
  double precision, parameter :: &
                   seconds_per_hour = 3600d0,       & ! Number of seconds per hour
                    seconds_per_day = 86400d0,      & ! Number of seconds per day
                  seconds_per_day_1 = 1.157407d-05    ! Inverse of seconds per day

  type model_working_variables
    ! Variables needed incase of using random forest functions.
    ! None are currently implemented but variables remain for legacy reasons
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
    ! Modile level ACM-GPP-ET variables
    double precision :: ci
    double precision, dimension(nos_soil_layers) :: soil_frac_clay, soil_frac_sand
  end type
  type(model_working_variables), allocatable, dimension(:):: mVs

contains
  !
  !--------------------------------------------------------------------
  !
end module carbon_model_memory
