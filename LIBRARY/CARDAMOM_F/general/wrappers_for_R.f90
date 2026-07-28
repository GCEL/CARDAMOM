!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CARbon DAta MOdel fraMework (CARDAMOM) and DALEC terrestrial ecosystem model suite
! CARDAMOM is a Bayesian model-data fusion software framework. CARDAMOM is used to
! assimilate observations and ecological theory to retrieve parameters for the
! DALEC suite of intermediate complexity terrestrial ecosystem models. DALEC can be
! used as a fully integrated component of CARDAMOM or independently.
! Copyright (C) 2024  University of Edinburgh,
!                     Mathew Williams (mat.williams@ed.ac.uk),
!                     T. Luke Smallman (t.l.smallman@ed.ac.uk)
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

!!!!!!!!!! File specific description !!!!!!!!!!
! Fortran wrappers exposing the CARDAMOM samplers to R via C binding.
!
! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory);
! translation to Fortran, integration and subsequent modifications by T. L. Smallman,
! J. F. Exbrayat and colleagues (University of Edinburgh). Sampler-layer development
! (DEMCz, parallel samplers, wrappers) by J. Klebes (University of Edinburgh), 2024-2025.
! See function/subroutine specific comments for exceptions and contributors.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module cardamom_Rinterfaces

   implicit none(external)

   public

   ! helpers for use from R-
   ! TODO different module and file

  contains
  !
  !--------------------------------------------------------------------
  !
  subroutine get_npars(npars) bind(c, name="C_getmodelnpars")

    use iso_c_binding

    ! return the compiled model's expected unmber of parameters

    !use model_parameters, only: pars_info
    use model_shared, only: PI
    implicit none(type, external)
    integer(c_int), intent(out)  :: npars
    !call pars_info(PI)
    npars = PI%npars

  end subroutine get_npars
  !
  !--------------------------------------------------------------------
  !
  subroutine get_parmin(npars, parmin) bind(c, name="C_getmodelparmin")
    use iso_c_binding
    !use model_parameters, only: pars_info
    use model_shared, only: PI

    ! get the compiled model's parmin list for R

    implicit none(type, external)

    integer(c_int), intent(in) :: npars
    real(c_double), dimension(npars), intent(out)  :: parmin
    !call pars_info(PI)
    parmin = PI%parmin

  end subroutine get_parmin
  !
  !--------------------------------------------------------------------
  !
  subroutine get_parmax(npars, parmax) bind(c, name="C_getmodelparmax")
    use iso_c_binding
    !use model_parameters, only: pars_info
    use model_shared, only: PI

    implicit none(type, external)

    integer(c_int), intent(in) :: npars
    real(c_double), dimension(npars), intent(out)  :: parmax
    !call pars_info(PI)
    parmax = PI%parmax

  end subroutine get_parmax
  !
  !--------------------------------------------------------------------
  !
  subroutine initialize_cardamom(nchains) bind(c, name="C_initialize_model")
    use iso_c_binding
    use cardamom_io, only: initialize
    use model_shared, only: initialize_carbon_model

    implicit none(type, external)

!TLS: This needs to be changed, depending on how it is used.
!     Not sure why there is a need to hardcode this path, as all information should be passed through the interface, not the filename

    !character(kind = c_char, len = *), intent(in) :: datain_filename
    character(kind=c_char, len=350) :: filename
    integer(c_int), intent(in), optional :: nchains
    integer(c_int) :: nchains_
    if (.not. present(nchains)) then
        nchains_ = 1_c_int
    else
       nchains_ = nchains
    end if
    !filename = trim(datain_filename)//c_null_char
    ! it turns out using string from R to C is hard so I'm hardcoding it here for now
    filename = "/home/jklebes/CARDAMOM/test/data/UK_baseline_sites_AliceHolt.bin"
    call initialize(filename)
    nchains_ = 4
    call initialize_carbon_model(nchains_)

  end subroutine initialize_cardamom
  !
  !--------------------------------------------------------------------
  !
  subroutine initialize_stresstest_circle() bind(c, name="C_initialize_stresstest_circle")
    use model_shared, only: PI, initialize_parinfo
    use MHMCMC_StressTests

    implicit none(type, external)

    character(len=350) :: choice, outfile
    choice = "Circle"
    call prepare_for_stress_test(choice, outfile)

   end subroutine initialize_stresstest_circle
  !
  !--------------------------------------------------------------------
  !
end module cardamom_Rinterfaces
