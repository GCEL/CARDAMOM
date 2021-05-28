module MODEL_PARAMETERS

  implicit none

  !!!!!!!!!!!
  ! Authorship contributions
  !
  ! This code is based on the original C verion of the University of Edinburgh
  ! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
  ! All code translation into Fortran, integration into the University of
  ! Edinburgh CARDAMOM code and subsequent modifications by:
  ! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
  ! See function / subroutine specific comments for exceptions and contributors
  !!!!!!!!!!!

  ! make all private
  private

  ! specify explicitly the public
  public :: pars_info

  contains

  !
  !------------------------------------------------------------------
  !
  subroutine pars_info
    use MCMCOPT, only: PI
    use cardamom_structures, only: DATAin

    ! Subroutine contains a list of parameter ranges for the model.
    ! These could or possibly should go into an alternate file which can be read in.
    ! This may improve the usability when it comes to reading these information
    ! in for different PFTs

    implicit none

    ! contains 6 fields with min max log for par and par

    !
    ! declare parameters
    !

    !
    ! Metabolic photosynthesis
    !

    ! Optimum nitrogen use efficiency (gC/gN per m2 at optimum temperature)
    ! Derived from Vcmax reported in Wullschleger (1993), Journal of
    ! Experimental Botany, Vol 44, No. 262, pp. 907-920.
    ! ~40 gC/gN/day
    ! TRY database equivalent 2.5 % = 1.648512; 97.5 % = 19.906560
    ! mean of 11.197440 +/-1.32313 Kattge et al., (2011)
    ! Xu et al., (2017):
    ! Variations of leaf longevity in tropical moist forests predicted by a
    ! trait-driven carbon optimality model,
    ! Ecology Letters, doi: 10.1111/ele.12804, upper value of 82 gC/gN/day
    ! Thus we will compromise on the value between these but closer to the
    ! newer estimate (i.e. 30 gC/gN/day)
    PI%parmin(1) = 01d0
    PI%parmax(1) = 42d0

    ! max temperature for photosynthesis (oC)
    ! SPA apparent value = 57.05oC
    PI%parmin(2) = 45d0
    PI%parmax(2) = 70d0
    ! optimum temperature for photosynthesis (oC)
    ! SPA value = 30oC
    PI%parmin(3) = 20d0
    PI%parmax(3) = 40d0
    ! kurtosis of photosynthesis temperature response
    PI%parmin(4) = 0.12d0
    PI%parmax(4) = 0.24d0

    ! light limited photosynthesis

    ! maximum canopy quantum yield (gC/MJ_PAR/m2/day)
    ! SPA apparent canopy quantum yield 3.2
    ! Observational constraints = ~ 1
    PI%parmin(5) = 1d0
    PI%parmax(5) = 8d0

    !
    ! canopy conductance (gc) drivers
    !

    ! leafWP-soilWP (MPa); actual SPA parameter = 2+gplant
    PI%parmin(6) = -3d0
    PI%parmax(6) = -1d0

    !
    ! Linear correction for soil isothermal to net radiation
    !

    ! Coefficient linking isothermal->net adjustment and LAI
    ! SPA based prior is -2.7108547 W/m2 (SE +/- 0.0222038)
    PI%parmin(7) = -3.5d0
    PI%parmax(7) =  0.0d0

    !
    ! GPP / transpiration optimisation
    !

    ! iWUE (gC/m2leaf/dayl/mmolH2Ogs/s)
    ! Actual value used in SPA is 8.4e-8 (Williams et al., 1996)
    ! Other reported values are 9.0e-8 -> 1.8e-7 (Bonan et al., 2014)
    ! NOTE: As this is applied at daily time step and the
    !       hourly time step activities in SPA operate across
    !       a non-linear diurnal cycle the true equivalent value is effectively unknown.
    PI%parmin(8) = 1d-7
    PI%parmax(8) = 1e-4

    !
    ! Soil shortwave radiation absorption
    !

    ! soil sw radiation absorption (fraction)
    PI%parmin(9) = 0.750d0
    PI%parmax(9) = 0.999d0

    !
    ! Canopy longwave escape from canopy
    !

    ! NOTE: The priors are based on the fraction of LW released by the canopy which is lost.
    !       The priors are very well constrained, however,
    !       the amount of long wave loss is different due to multi-layer vs bulk canopy

    ! Max fractional reduction of longwave release from canopy
    ! Prior from offline SPA calibration = 0.9517081 +/- 0.0001011 SE
    PI%parmin(10) = 0.25d0
    PI%parmax(10) = 1d0

    ! LAI adjustment for long wave release from canopy
    ! Prior from offline SPA calibration = 4.6917871 +/- 0.0013296 SE
    PI%parmin(11) = 2.0d0
    PI%parmax(11) = 6.0d0

    !
    ! Linear correction for soil isothermal to net radiation
    ! Cont.
    !

    ! The magnitude of adjustment between soil isothermal net long wave radiation is found to have
    ! a linear relationship with absorbed short-wave radiation and LAI.
    ! NOTE: in constrast to the canopy relationship (positive) the soil one has two distinct
    ! negative linear relationships divided based on whether the soil is saturated or not.
    ! This is not accounted for here

    ! Coefficient linking isothermal->net adjustment and absorbed SW
    ! SPA based prior is -0.0357603 W/m2 (SE +/- 0.0004877)
    PI%parmin(12) = -0.04d0
    PI%parmax(12) =  0.00d0

    ! Constant relating isothermal->net adjustment soil radiation (W/m2)
    ! SPA based prior is 3.4806352 W/m2 (SE +/- 0.063849)
    PI%parmin(13) = 0d0
    PI%parmax(13) = 10d0

    !
    ! Max canopy intercepted PAR and NIR transmittance / reflectance
    !

    ! Max PAR transmittance
    ! SPA value = 0.16 (Sitka Spruce 0.16, grass/crop ~ 0.38)
    PI%parmin(14) = 0.10d0
    PI%parmax(14) = 0.40d0

    ! Max NIR transmittance
    ! SPA value = 0.26
    PI%parmin(15) = 0.10d0
    PI%parmax(15) = 0.40d0

    ! Max PAR reflectance
    ! SPA value = 0.16 (Sitka Spruce 0.07, grass/crop ~ 0.11)
    PI%parmin(16) = 0.05d0
    PI%parmax(16) = 0.40d0

    ! Max NIR reflectance
    ! SPA value = 0.43 (Sitka Spruce 0.16, grass/crop ~ 0.38)
    PI%parmin(17) = 0.10d0
    PI%parmax(17) = 0.50d0

    !
    ! Linear correction for canopy isothermal to net radiation
    !

    ! The magnitude of adjustment between canopy isothermal net long wave radiation (W/m2) is found to have
    ! a linear relationship (R2 = 0.89) with absorbed short-wave radiation (W/m2) and LAI.

    ! Coefficient linking isothermal->net adjustment and absorbed SW
    ! SPA based prior is 0.0154225 (SE +/- 0.0005798)
    PI%parmin(18) =  0.0d0
    PI%parmax(18) =  0.1d0

    ! Constant relating isothermal->net adjustment canopy radiation (W/m2)
    ! SPA based prior is 0.0577857 W/m2 (SE +/- 0.0217731)
    PI%parmin(19) =  0.0d0
    PI%parmax(19) =  1.0d0

    ! Coefficient linking isothermal->net adjustment and LAI (m2/m2)
    ! SPA based prior is 2.4526437 W/m2 (SE +/- 0.0229691)
    PI%parmin(20) =  0.0d0
    PI%parmax(20) =  4.0d0

  end subroutine pars_info
  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
