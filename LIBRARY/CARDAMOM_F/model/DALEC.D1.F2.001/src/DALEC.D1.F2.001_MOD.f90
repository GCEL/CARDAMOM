!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CARbon DAta MOdel fraMework (CARDAMOM) and DALEC terrestrial ecosystem model suite
! CARDAMOM is a Bayesian model-data fusion software framework. CARDAMOM is used to
! assimilate observations and ecological theory to retrieve parameters for the
! DALEC suite of intermediate complexity terrestrial ecosystem models. DALEC can be
! used as a fully integrated component of CARDAMOM or independently.
! Copyright (C) 2024  University of Edinburgh,
!                     Mathew Williams (mat.williams@ed.ac.uk),
!                     T. Luke Smallman (t.l.smallman@ed.ac.uk),
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
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

!!!!!!!!!!!! File specific description !!!!!!!!!!
! This file contains the parameters and memory of DALEC.D1.F2.001
!
! This code contains a variant of the Data Assimilation Linked ECosystem (DALEC) model.
! This version of DALEC is derived from the following primary references:
! Williams et al., (2005), doi: 10.1111 /j.1365-2486.2004.091.x
! This code is based on that created by A. A. Bloom (UoE, now at JPL, USA).
! Subsequent modifications by:
! T. L. Smallman (University of Edinburgh, t.l.smallman@ed.ac.uk)
! See function / subroutine specific comments for exceptions and contributors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module carbon_model_memory

implicit none

! make all private
!private
public

  !!!!!!!!!
  ! Parameters
  !!!!!!!!!

  ! for consisteny between requirements of different models
  integer, parameter :: nos_root_layers = 2, nos_soil_layers = nos_root_layers + 1

  double precision, parameter :: pi = 3.1415927d0, &
                         deg_to_rad = 0.01745329d0   ! pi/180d0

  !!!!!!!!!
  ! Module variables
  !!!!!!!!!

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
