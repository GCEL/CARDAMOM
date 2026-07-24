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
! Wraps the chosen model's likelihood subroutines to the generic sampler interface (including C binding for use from R).
!
! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory);
! translation to Fortran, integration and subsequent modifications by T. L. Smallman,
! J. F. Exbrayat and colleagues (University of Edinburgh). Sampler-layer development
! (DEMCz, parallel samplers, wrappers) by J. Klebes (University of Edinburgh), 2024-2025.
! See function/subroutine specific comments for exceptions and contributors.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module model_likelihood_wrapper

  ! This module wraps the module_likelihood subroutine found in
  ! model/<chosen model>/likelihood/model_likelihood.f90
  ! (whichever is chosen for compilation) to a standard-shaped function
  ! model_likelihood : params vector -> loglikelihood value
  ! that can be used as input to cardamom_samplers

  use model_likelihood_module, only: model_likelihood, &
                                     edc_model_likelihood, scaled_model_likelihood

  use iso_c_binding
  implicit none(type, external)
  public

  contains
  !
  !--------------------------------------------------------------------
  !
  subroutine model_likelihood_fct(params, npars, loglikelihood, id) bind(c, name="C_modellikelihood")
    use iso_c_binding

    implicit none(type, external)

    ! model_likelihood, wrapped to be
    ! shaped like the generic function to hand to both cardamom_samplers and (via C binding) R
    ! For R : has to be subroutine ( -> C void function), have to give npars

    ! Arguments
    integer(c_int), intent(in)  :: npars, id
    real(c_double), intent(inout), dimension(npars) :: params
    real(c_double), intent(out) :: loglikelihood

    ! Local variables
    real(c_double) :: ML_obs_out, ML_prior_out

    ! call the function to write to ML_obs_out and ML_prior_out
    call model_likelihood(params, ML_obs_out, ML_prior_out, id)

    ! for the purpose of running samplers we are only interested in the sum
    loglikelihood = ML_obs_out + ML_prior_out

  end subroutine model_likelihood_fct
  !
  !--------------------------------------------------------------------
  !
  subroutine scaled_model_likelihood_fct(params, npars, loglikelihood, id) bind(c, name="C_scaledmodellikelihood")
    use iso_c_binding

    implicit none(type, external)

    ! Wrapper for all scaled model likelihood function variants

    ! Arguments
    integer(c_int), intent(in)  :: npars, id
    real(c_double), intent(inout), dimension(npars) :: params
    real(c_double), intent(out) :: loglikelihood

    ! Local variables
    real(c_double) :: ML_obs_out, ML_prior_out

    ! call the function to write to ML_obs_out and ML_prior_out
    call scaled_model_likelihood(params, ML_obs_out, ML_prior_out, id)

    ! for the purpose of running samplers we are only interested in the sum
    loglikelihood = ML_obs_out + ML_prior_out

  end subroutine scaled_model_likelihood_fct
  !
  !--------------------------------------------------------------------
  !
  subroutine edc_model_likelihood_fct(params, npars, loglikelihood, id) bind(c, name="C_edcmodellikelihood")
    use iso_c_binding

    implicit none(type, external)

    ! Wrapper for EDC function, which finds a set of parameters fulfilling otherwise
    ! hard boundary conditions with a softer stepped potential

    ! Arguments
    integer(c_int), intent(in)  :: npars, id
    real(c_double), intent(inout), dimension(npars) :: params
    real(c_double), intent(out) :: loglikelihood

    ! Local variables
    real(c_double) :: ML_obs_out, ML_prior_out

    call edc_model_likelihood(params, ML_obs_out, ML_prior_out, id)

    loglikelihood = ML_obs_out + ML_prior_out

  end subroutine edc_model_likelihood_fct
  !
  !--------------------------------------------------------------------
  !
end module model_likelihood_wrapper
