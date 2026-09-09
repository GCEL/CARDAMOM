module model_likelihood_wrapper

! jklebes 2024
! This module wraps the module_likelihood subroutine found in
! model/<chosen model>/likelihood/model_likelihood.f90
! (whichever is chosen for compilation) to a standard-shaped function
! model_likelihood : params vector -> loglikelihood value
! that can be used as input to cardamom_samplers

use model_likelihood_module, only: model_likelihood, &
                            edc_model_likelihood, find_edc_initial_values, &
          sub_model_likelihood, sqrt_model_likelihood, log_model_likelihood!, log_model_likelihood_dtm
use iso_c_binding
implicit none

contains


! model_likelihood, wrapped to be
! shaped like the generic function to hand to both cardamom_samplers and (via C binding) R
! For R : has to be subroutine ( -> C void function), have to give npars
!
subroutine model_likelihood_fct(params, npars, loglikelihood) bind(c, name="C_modellikelihood")
  use iso_c_binding
  implicit none
  integer(c_int), intent(in)  :: npars
  real(c_double), intent(inout), dimension(npars):: params
  real(c_double), intent(out):: loglikelihood

  real(c_double):: ML_obs_out, ML_prior_out
  !TODO can we use an expected_npars from model files?
  !if (npars .neq. expected_npars) then
    !write(*,*) "Error : Passed ", npars, "from R (as indicated by second argument to modellikelihood), but this model takes a vector
    !of " , expected_npars, "values."
    ! TODO force exit ?
!else:

  ! call the function to write to ML_obs_out and ML_prior_out
  call model_likelihood(params, ML_obs_out, ML_prior_out)

  ! for the purpose of running samplers we are only interested in the sum
  loglikelihood = ML_obs_out+ML_prior_out

!endif
end subroutine model_likelihood_fct

! variants
subroutine log_model_likelihood_fct(params, npars, loglikelihood) bind(c, name="C_logmodellikelihood")
  use iso_c_binding
  implicit none
  integer(c_int), intent(in)  :: npars
  real(c_double), intent(inout), dimension(npars):: params
  real(c_double), intent(out):: loglikelihood

  real(c_double):: ML_obs_out, ML_prior_out

  call log_model_likelihood(params, ML_obs_out, ML_prior_out)

  loglikelihood = ML_obs_out+ML_prior_out
end subroutine log_model_likelihood_fct

subroutine sub_model_likelihood_fct(params, npars, loglikelihood) bind(c, name="C_submodellikelihood")
  use iso_c_binding
  implicit none
  integer(c_int), intent(in)  :: npars
  real(c_double), intent(inout), dimension(npars):: params
  real(c_double), intent(out):: loglikelihood

  real(c_double):: ML_obs_out, ML_prior_out

  call sub_model_likelihood(params, ML_obs_out, ML_prior_out)

  loglikelihood = ML_obs_out+ML_prior_out
end subroutine sub_model_likelihood_fct

subroutine sqrt_model_likelihood_fct(params, npars, loglikelihood) bind(c, name="C_sqrtmodellikelihood")
  use iso_c_binding
  implicit none
  integer(c_int), intent(in)  :: npars
  real(c_double), intent(inout), dimension(npars):: params
  real(c_double), intent(out):: loglikelihood

  real(c_double):: ML_obs_out, ML_prior_out

  call sqrt_model_likelihood(params, ML_obs_out, ML_prior_out)

  loglikelihood = ML_obs_out+ML_prior_out
end subroutine sqrt_model_likelihood_fct

! the function for starting loops, which finds a set of parameters fulfilling otherwise
! hard boundary conditions with a softer potential
subroutine edc_model_likelihood_fct(params, npars, loglikelihood) bind(c, name="C_edcmodellikelihood")
  use iso_c_binding
  implicit none
  integer(c_int), intent(in)  :: npars
  real(c_double), intent(inout), dimension(npars):: params
  real(c_double), intent(out):: loglikelihood

  real(c_double):: ML_obs_out, ML_prior_out

  call edc_model_likelihood(params, ML_obs_out, ML_prior_out)

  loglikelihood = ML_obs_out+ML_prior_out
end subroutine edc_model_likelihood_fct

! helpers for use from R-
! TODO different module and file

! return the compiled model's expected unmber of parameters
subroutine get_npars(npars) bind(c, name="C_getmodelnpars")
  use iso_c_binding
  use model_parameters , only: pars_info
  use MCMCOPT, only: PI
  implicit none
  integer(c_int), intent(out)  :: npars
  call pars_info() !ideally this would bt written in an object, not a function in model_parameters
  npars =  PI%npars
end subroutine get_npars


! get the compiled model's parmin list for R
subroutine get_parmin(npars, parmin) bind(c, name="C_getmodelparmin")
  use iso_c_binding
  use model_parameters , only: pars_info
  use MCMCOPT, only: PI
  implicit none
  integer(c_int) , intent(in) :: npars
  real(c_double), dimension(npars), intent(out)  :: parmin
  call pars_info()
  parmin =  PI%parmin
end subroutine get_parmin

subroutine get_parmax(npars, parmax) bind(c, name="C_getmodelparmax")
    use iso_c_binding
    use model_parameters , only: pars_info
    use MCMCOPT, only: PI
    implicit none
    integer(c_int) , intent(in) :: npars
    real(c_double), dimension(npars), intent(out)  :: parmax
    call pars_info()
    parmax =  PI%parmax
    end subroutine get_parmax

subroutine inintialize_example_FI_Hyy() bind(c, name = "C_TMP_initialize")
  ! trigger read_binary_data(hard coded filename) for testing,
  ! later do this better
  use cardamom_io, only: initialize
  implicit none
  character(len=350)  :: filename
  filename = "/home/jklebes/cardamom_profiling/FI-Hyy_example/DATA/FI-Hyy_example_FI-Hyy.bin"
  ! TODO add results
  call initialize(filename)
end subroutine inintialize_example_FI_Hyy

!trigger find_edc_initial_values
subroutine set_parini() bind(c, name = "C_setparini")
  ! trigger read_binary_data(hard coded filename) for testing,
  ! later do this better
  use cardamom_io, only: initialize
  implicit none
  call find_edc_initial_values
end subroutine set_parini

! get the compiled model's parini (after initialization) list for R
subroutine get_parini(npars, parini) bind(c, name="C_getexampleparini")
  use iso_c_binding
  use MCMCOPT, only: PI
  use cardamom_structures, only: DATAin
  implicit none
  integer(c_int) , intent(in) :: npars
  real(c_double), dimension(npars), intent(out)  :: parini
  PI%parini(1:PI%npars) = DATAin%parpriors(1:PI%npars)
  parini = PI%parini
end subroutine get_parini


end module model_likelihood_wrapper
