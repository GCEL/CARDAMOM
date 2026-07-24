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
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see < https://www.gnu.org/licenses/>.
!
!!!!!!!!!!!! File specific description !!!!!!!!!!
! Code responsible for input/output operations for CARDAMOM
!
! This code is based on the original C verion of the University of Edinburgh
! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
! All code translation into Fortran, integration into the University of
! Edinburgh CARDAMOM code and subsequent modifications by:
! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
! J. F. Exbrayat (University of Edinburgh)
! See function/subroutine specific comments for exceptions and contributors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module cardamom_io

   ! Module contains subroutines and variables needed to output parameter,
   ! likelihood and step size information from the MCMC algorithms.

   implicit none(type, external)

   ! declare private
   private

   ! allow access to specific functions
   public ::  update_obs_scaling_normal &
            , update_obs_scaling_nsamples &
            , update_obs_scaling_sqrt_nsamples &
            , update_obs_scaling_log_nsamples &
            , cardamom_model_library &
            , read_options &
            , read_binary_data &
            , initialize

   ! declare module level variables
   integer :: pfile_unit = 10, sfile_unit = 11, cfile_unit = 12, cifile_unit = 13, ifile_unit = 14

   ! parameters
   integer, parameter :: real_bytes = 8  ! number of bytes in real variable, 8 bytes is to make double precision

   save

contains
   !
   !--------------------------------------------------------------------
   !
   subroutine cardamom_model_library(DATAin)
      use cardamom_structures, only: DATA_type
      implicit none(type, external)
      type(DATA_type), intent(inout) :: DATAin
    !! a local DATAin object as argument


    ! don't forget to update values found in the relevant model *_PARS.f90

    ! choose between included model arrangements
    ! NOTE: negative values are coded elsewhere and reserved for MCMC stress
    ! testing
    if (DATAin%ID == 0) then
        ! ID = 0 - ACM/ACM-ET
        DATAin%nopools = 2
        DATAin%nopars = 20
        DATAin%nofluxes = 4
        DATAin%nodiags = 20 ! Initial value, will need updating
    else if (DATAin%ID == 1) then
        ! ID = 1 - DALEC.D1.F2.001
        DATAin%nopools = 5
        DATAin%nopars = 22
        DATAin%nofluxes = 35
        DATAin%nodiags = 20 ! Initial value, will need updating
    else if (DATAin%ID == 2) then
        ! ID = 2 - DALEC.C1.D1.F2.P1.002
        DATAin%nopools = 6
        DATAin%nopars = 28
        DATAin%nofluxes = 39
        DATAin%nodiags = 20 ! Initial value, will need updating
    else if (DATAin%ID == 3 ) then
        ! ID = 3 - DALEC.A1.C1.D2.F2.H1.P1.003
        DATAin%nopools = 6
        DATAin%nopars = 31
        DATAin%nofluxes = 51
        DATAin%nodiags = 20 ! Initial value, will need updating
    else if (DATAin%ID == 4) then
        ! ID = 4 - DALEC.A1.C1.D2.F2.H2.P1.004
        DATAin%nopools = 7
        DATAin%nopars = 32
        DATAin%nofluxes = 51
        DATAin%nodiags = 24 ! Initial value, will need updating
    else if (DATAin%ID == 5) then
        ! ID = 5 - DALEC.A1.C1.D2.F2.H2.P1.R1.005
        DATAin%nopools = 7
        DATAin%nopars = 32
        DATAin%nofluxes = 52
        DATAin%nodiags = 20 ! Initial value, will need updating
    else if (DATAin%ID == 6) then
        ! ID = 6 - DALEC.A1.C2.D2.F2.H2.P1.R1.006
        DATAin%nopools = 8
        DATAin%nopars = 35
        DATAin%nofluxes = 57
        DATAin%nodiags = 20 ! Initial value, will need updating
    else if (DATAin%ID == 7) then
        ! ID = 7 - DALEC.A1.C2.D2.F2.H2.P2.R1.007
        DATAin%nopools = 8
        DATAin%nopars = 36
        DATAin%nofluxes = 57
        DATAin%nodiags = 20 ! Initial value, will need updating
    else if (DATAin%ID == 8) then
        ! ID = 8 - DALEC.A1.C2.D2.F2.H1.P3.R1.008
        DATAin%nopools = 7
        DATAin%nopars = 39
        DATAin%nofluxes = 57
        DATAin%nodiags = 21 ! Initial value, will need updating
    else if (DATAin%ID == 9) then
        ! ID = 9 - DALEC.A1.C2.D2.F2.H2.P3.R1.009
        DATAin%nopools = 8
        DATAin%nopars = 40
        DATAin%nofluxes = 57
        DATAin%nodiags = 21 ! Initial value, will need updating
    else if (DATAin%ID == 10) then
        ! ID = 10 - DALEC.A1.C2.D2.F2.H1.P4.R2.010
        DATAin%nopools = 8
        DATAin%nopars = 43
        DATAin%nofluxes = 57
        DATAin%nodiags = 21 ! Initial value, will need updating
    else if (DATAin%ID == 11) then
        ! ID = 11 - DALEC.A1.C2.D2.F2.H2.P4.R2.011
        DATAin%nopools = 8
        DATAin%nopars = 43
        DATAin%nofluxes = 57
        DATAin%nodiags = 21 ! Initial value, will need updating
    else if (DATAin%ID == 12) then
        ! ID = 12 - DALEC.C4.D1.F2.012
        DATAin%nopools = 3
        DATAin%nopars = 15
        DATAin%nofluxes = 28
        DATAin%nodiags = 20 ! Initial value, will need updating
    else if (DATAin%ID == 13) then
        ! ID = 13 - DALEC.C5.D1.F2.P1.013
        DATAin%nopools = 4
        DATAin%nopars = 21
        DATAin%nofluxes = 32
        DATAin%nodiags = 20 ! Initial value, will need updating
    else if (DATAin%ID == 14) then
        ! ID = 14 - DALEC.C3.M1.014
        DATAin%nopools = 9
        DATAin%nopars = 37
        DATAin%nofluxes = 42
        DATAin%nodiags = 20 ! Initial value, will need updating
    else if (DATAin%ID == 15) then
        ! ID = 15 - DALEC.A3.C3.H2.M1.015 i.e. the CROP model
        DATAin%nopools = 10
        DATAin%nopars = 38
        DATAin%nofluxes = 48
        DATAin%nodiags = 20 ! Initial value, will need updating
    else if (DATAin%ID == 16) then
        ! ID = 16 - DALEC.A3.H1.M2.016
        DATAin%nopools = 5
        DATAin%nopars = 34 ! to update
        DATAin%nofluxes = 57
        DATAin%nodiags = 23 ! Initial value, will need updating
    else if (DATAin%ID == 17) then
        ! ID = 17 - DALEC.A3.H2.M2.017
        DATAin%nopools = 6
        DATAin%nopars = 35
        DATAin%nofluxes = 57
        DATAin%nodiags = 23 ! Initial value, will need updating
    else if (DATAin%ID == 18) then
        ! ID = 18 - DALEC.A1.C1.D2.F2.H2.P2.018
        DATAin%nopools = 7
        DATAin%nopars = 33
        DATAin%nofluxes = 51
        DATAin%nodiags = 20 ! Initial value, will need updating
    else if (DATAin%ID == 19) then
        ! ID = 19 - DALEC.A1.C2.D2.F2.H2.P2.R3.019
        DATAin%nopools = 8
        DATAin%nopars = 38
        DATAin%nofluxes = 57
        DATAin%nodiags = 20 ! Initial value, will need updating
    else if (DATAin%ID == 20) then
        ! ID = 20 - DALEC.A2.C1.D2.F2.H2.P1.020
        DATAin%nopools = 7
        DATAin%nopars = 32
        DATAin%nofluxes = 51
        DATAin%nodiags = 20 ! Initial value, will need updating
    else if (DATAin%ID == 21) then
        ! ID = 21 - DALEC.A1.C1.D2.F2.H2.P5.021
        DATAin%nopools = 7
        DATAin%nopars = 33
        DATAin%nofluxes = 51
        DATAin%nodiags = 24 ! Initial value, will need updating
    else if (DATAin%ID == 22) then
        ! ID = 22 - DALEC.A1.C1.D2.F2.H2.P6.022
        DATAin%nopools = 7
        DATAin%nopars = 34
        DATAin%nofluxes = 51
        DATAin%nodiags = 20 ! Initial value, will need updating
    else if (DATAin%ID == 23) then
        ! ID = 23 - DALEC.A1.C2.D2.F2.H2.P7.R2.023
        DATAin%nopools = 8
        DATAin%nopars = 48
        DATAin%nofluxes = 57
        DATAin%nodiags = 23 ! Initial value, will need updating
    else if (DATAin%ID == 24) then
        ! ID = 24 - DALEC.A1.C1.D2.F2.H4.P1.024
        DATAin%nopools = 7
        DATAin%nopars = 37
        DATAin%nofluxes = 51
        DATAin%nodiags = 24 ! Initial value, will need updating
    else if (DATAin%ID == 25) then
        ! ID = 25 - DALEC...025 NOT IN USE
        DATAin%nopools = 8
        DATAin%nopars = 48
        DATAin%nofluxes = 57
        DATAin%nodiags = 23 ! Initial value, will need updating
    else if (DATAin%ID == 26) then
        ! ID = 26 - DALEC.A4.C6.D2.F2.H3.P10.026
        DATAin%nopools = 8
        DATAin%nopars = 49
        DATAin%nofluxes = 53
        DATAin%nodiags = 35 ! Initial value, will need updating
    else if (DATAin%ID == 27) then
        ! ID = 27 - DALEC_1005
        DATAin%nopools = 8
        DATAin%nopars = 38
        DATAin%nofluxes = 43
        DATAin%nodiags = 20 ! Initial value, will need updating
    else if (DATAin%ID == 28) then
        ! ID = 28 - DALEC_1005a
        DATAin%nopools = 8
        DATAin%nopars = 38
        DATAin%nofluxes = 43
        DATAin%nodiags = 20 ! Initial value, will need updating
    else if (DATAin%ID == 29) then
        ! ID = 29 -DALEC.A1.C1.D2.F2.H3.P1.029
        DATAin%nopools = 7
        DATAin%nopars = 33
        DATAin%nofluxes = 51
        DATAin%nodiags = 20 ! Initial value, will need updating
    else if (DATAin%ID == 30) then
        ! ID = 30 - DALEC.A3.C1.D2.F2.H2.P1.030
        DATAin%nopools = 7
        DATAin%nopars = 38
        DATAin%nofluxes = 51
        DATAin%nodiags = 24 ! Initial value, will need updating
    else if (DATAin%ID == 31) then
        ! ID = 31 - DALEC.A4.C6.D2.F2.H2.P11.031
        DATAin%nopools = 7
        DATAin%nopars = 43
        DATAin%nofluxes = 51
        DATAin%nodiags = 24
    else if (DATAin%ID == 32) then
        ! ID = 32 - DALEC.A1.C1.D2.F2.H6.P1.R5.032
        DATAin%nopools = 7
        DATAin%nopars = 37
        DATAin%nofluxes = 51
        DATAin%nodiags = 24 ! Initial value, will need updating
    else if (DATAin%ID == 33) then
        ! ID = 33 - DALEC.A4.C6.D2.F2.H3.P12.033
        DATAin%nopools = 7
        DATAin%nopars = 46
        DATAin%nofluxes = 51
        DATAin%nodiags = 30
    else if (DATAin%ID == 34) then
        ! ID = 34 -
        print*,"Model ID 34 not currently in use"
        stop
    else if (DATAin%ID == 35) then
        ! ID = 35
        print*,"Model ID 35 not currently in use"
        stop
    else if (DATAin%ID == 36) then
        ! ID = 36 - DALEC.A1.C7.D2.F2.H2.P1.R4.036
        DATAin%nopools = 11
        DATAin%nopars = 50
        DATAin%nofluxes = 67
        DATAin%nodiags = 20
    else if (DATAin%ID == 37) then
        ! ID = 37 - DALEC.A1.C1.D2.F2.H5.P1.037
        DATAin%nopools = 7
        DATAin%nopars = 34
        DATAin%nofluxes = 51
        DATAin%nodiags = 24 ! Initial value, will need updating
    else if (DATAin%ID == 38) then
        ! ID = 38 -
        print*,"Model ID 38 not currently in use"
        stop
    else if (DATAin%ID == 39) then
        ! ID = 39 -
        print*,"Model ID 39 not currently in use"
        stop
    else
        write(*,*) "Oh dear... model ID not valid = ",DATAin%ID
        stop
    end if

   end subroutine cardamom_model_library
   !
   !------------------------------------------------------------------
   !
   subroutine read_binary_data(infile, DATAin)
      use cardamom_structures, only: DATA_type
      use CARBON_MODEL_MOD, only: nos_soil_layers

      ! subroutine opens and reads the binary data files provided by/for the
      ! CARDAMOM framework. This data is then loaded into the DATAin type

      ! TEMPLATE FOR ALL DALEC MCMC DATA files
      ! Static Elements: 1-100-use as many as needed
      ! Parameter Priors: 101-150
      ! Parameter prior uncertainty: 151-200
      ! Other priors & uncertainties: 201-300
      ! TEMPORAL DRIVERS & DATA: 301-end

      implicit none(type, external)

      type(DATA_type), intent(inout) :: DATAin
    !! This is a temporary local DATAin object whose fields we are allowed to modify.
    !! After everything is read, the (protected) global DATAin object
    !! is set by copying this object.

      ! declare input variables
      character(350), intent(in) :: infile

      ! declare local variables
      integer :: nopars_dummy, subsample
      integer :: a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z, day &
                , start &
                , finish &
                , totcol & ! total number of columns (met+obs)
                , totread      ! total number of records already read
      double precision :: mz
      double precision, save :: subsample_fraction = 0.20  ! startd at 0.25
      double precision, dimension(:), allocatable :: statdat & ! static data input
         , mettemp & ! met data input
         , obstemp   ! obs data input

      write (*, *) "Input file to be read = ", trim(infile)

      ! open the binary file, with direct access for binary (unformatted) at
      ! double precision (double precision = 64 bytes)
      open (unit=ifile_unit, file=trim(infile), form="UNFORMATTED", access="stream", status="old")
      rewind (ifile_unit)

      ! allocate memory
      allocate (statdat(50))

      ! now read the static elements (1-50)
      do i = 1, 50  ! number of static elements
         read (ifile_unit) statdat(i)
      end do

      ! allocate the default run information
      DATAin%ID = int(statdat(1))
      DATAin%LAT = statdat(2)
      DATAin%nodays = int(statdat(3))
      DATAin%nomet = int(statdat(4))
      DATAin%noobs = int(statdat(5))
      DATAin%EDC = int(statdat(6))
      DATAin%PFT = int(statdat(7))
      DATAin%yield = -9999  ! int(statdat(8))
      DATAin%age = int(statdat(9))
      nopars_dummy = int(statdat(10))  ! needed for next dev stage
      ! Assume 3 soil layers only and that the
      ! top soil layer is assigned the top soil condition
      allocate (DATAin%soil_frac_sand(nos_soil_layers))
      allocate (DATAin%soil_frac_clay(nos_soil_layers))
      DATAin%soil_frac_sand(1) = statdat(12)  ! top soil sand percentage
      DATAin%soil_frac_sand(2:nos_soil_layers) = statdat(13)  ! bot
      DATAin%soil_frac_clay(1) = statdat(14)  ! top soil clay percentage
      DATAin%soil_frac_clay(2:nos_soil_layers) = statdat(15)  ! bot
      ! call for model specific values
      call cardamom_model_library(DATAin)

      ! allocate case specific information
      DATAin%edc_random_search = int(statdat(11))

      ! Do some sanity checks
      if (DATAin%lat > 90 .or. DATAin%lat < -90) then
         print *, "Latitude provided is not-90/90"
         stop
      end if

      ! clean up
      deallocate (statdat)

      ! read in parameter information (100 elements)
      a = 1
      do i = 51, 150
         read (ifile_unit) DATAin%parpriors(a)
         a = a + 1
      end do

      ! read in parameter uncertainty (100 elements)
      a = 1
      do i = 151, 250
         read (ifile_unit) DATAin%parpriorunc(a)
         a = a + 1
      end do

      ! read in parameter effect period (100 elements)
      a = 1
      do i = 251, 350
         read (ifile_unit) DATAin%parpriorweight(a)
         a = a + 1
      end do

      ! read in 'other' parameter priors (50 elements)
      a = 1
      do i = 351, 400
         read (ifile_unit) DATAin%otherpriors(a)
         a = a + 1
      end do

      ! read in 'other' parameter priors uncertainties (50 elements)
      a = 1
      do i = 401, 450
         read (ifile_unit) DATAin%otherpriorunc(a)
         a = a + 1
      end do

      ! read in 'other' parameter priors weighting (50 elements)
      a = 1
      do i = 451, 500
         read (ifile_unit) DATAin%otherpriorweight(a)
         a = a + 1
      end do

      ! Add a sensible limit on the weighting value
      where (DATAin%parpriorweight < 1) DATAin%parpriorweight = 1
      where (DATAin%otherpriorweight < 1) DATAin%otherpriorweight = 1

      ! now we know specific information about the dimensions in the file lets use
      ! it to allocate to the module variables
      allocate (mettemp(DATAin%nomet), obstemp(DATAin%noobs) &
                , DATAin%met(DATAin%nomet, DATAin%nodays) &
                , DATAin%GPP(DATAin%nodays), DATAin%GPP_unc(DATAin%nodays), DATAin%GPP_lag(DATAin%nodays) &
                , DATAin%NEE(DATAin%nodays), DATAin%NEE_unc(DATAin%nodays), DATAin%NEE_lag(DATAin%nodays) &
                , DATAin%LAI(DATAin%nodays), DATAin%LAI_unc(DATAin%nodays), DATAin%LAI_lag(DATAin%nodays) &
                , DATAin%Reco(DATAin%nodays), DATAin%Reco_unc(DATAin%nodays), DATAin%Reco_lag(DATAin%nodays) &
                , DATAin%Cfol_stock(DATAin%nodays), DATAin%Cfol_stock_unc(DATAin%nodays), DATAin%Cfol_stock_lag(DATAin%nodays) &
                , DATAin%Cwood_stock(DATAin%nodays), DATAin%Cwood_stock_unc(DATAin%nodays), DATAin%Cwood_stock_lag(DATAin%nodays) &
              , DATAin%Croots_stock(DATAin%nodays), DATAin%Croots_stock_unc(DATAin%nodays), DATAin%Croots_stock_lag(DATAin%nodays) &
                , DATAin%Clit_stock(DATAin%nodays), DATAin%Clit_stock_unc(DATAin%nodays), DATAin%Clit_stock_lag(DATAin%nodays) &
                , DATAin%Csom_stock(DATAin%nodays), DATAin%Csom_stock_unc(DATAin%nodays), DATAin%Csom_stock_lag(DATAin%nodays) &
                , DATAin%Cagb_stock(DATAin%nodays), DATAin%Cagb_stock_unc(DATAin%nodays), DATAin%Cagb_stock_lag(DATAin%nodays) &
,DATAin%Ccoarseroot_stock(DATAin%nodays), DATAin%Ccoarseroot_stock_unc(DATAin%nodays), DATAin%Ccoarseroot_stock_lag(DATAin%nodays) &
                , DATAin%Evap(DATAin%nodays), DATAin%Evap_unc(DATAin%nodays), DATAin%Evap_lag(DATAin%nodays) &
                , DATAin%SWE(DATAin%nodays), DATAin%SWE_unc(DATAin%nodays), DATAin%SWE_lag(DATAin%nodays) &
                , DATAin%NBE(DATAin%nodays), DATAin%NBE_unc(DATAin%nodays), DATAin%NBE_lag(DATAin%nodays) &
                , DATAin%Fire(DATAin%nodays), DATAin%Fire_unc(DATAin%nodays), DATAin%Fire_lag(DATAin%nodays) &
                , DATAin%fAPAR(DATAin%nodays), DATAin%fAPAR_unc(DATAin%nodays), DATAin%fAPAR_lag(DATAin%nodays) &
                , DATAin%Cwood_inc(DATAin%nodays), DATAin%Cwood_inc_unc(DATAin%nodays), DATAin%Cwood_inc_lag(DATAin%nodays) &
              , DATAin%Cwood_growth(DATAin%nodays), DATAin%Cwood_growth_unc(DATAin%nodays), DATAin%Cwood_growth_lag(DATAin%nodays) &
     , DATAin%Cwood_mortality(DATAin%nodays), DATAin%Cwood_mortality_unc(DATAin%nodays), DATAin%Cwood_mortality_lag(DATAin%nodays) &
                , DATAin%harvest(DATAin%nodays), DATAin%harvest_unc(DATAin%nodays), DATAin%harvest_lag(DATAin%nodays) &
,DATAin%foliage_to_litter(DATAin%nodays), DATAin%foliage_to_litter_unc(DATAin%nodays), DATAin%foliage_to_litter_lag(DATAin%nodays) &
                , DATAin%soilwater(DATAin%nodays), DATAin%soilwater_unc(DATAin%nodays), DATAin%soilwater_lag(DATAin%nodays))

    !! Zero all variables
    ! Drivers
    DATAin%met = 0d0
    ! Observations which have implicit lag of 0, i.e. they are relevant for the loaded time step
    DATAin%GPP = 0d0               ; DATAin%GPP_unc = 0d0               ; DATAin%GPP_lag = 0
    DATAin%NEE = 0d0               ; DATAin%NEE_unc = 0d0               ; DATAin%NEE_lag = 0
    DATAin%LAI = 0d0               ; DATAin%LAI_unc = 0d0               ; DATAin%LAI_lag = 0
    DATAin%Reco = 0d0              ; DATAin%Reco_unc = 0d0              ; DATAin%Reco_lag = 0
    DATAin%Cfol_stock = 0d0        ; DATAin%Cfol_stock_unc = 0d0        ; DATAin%Cfol_stock_lag = 0
    DATAin%Cwood_stock = 0d0       ; DATAin%Cwood_stock_unc = 0d0       ; DATAin%Cwood_stock_lag = 0
    DATAin%Croots_stock = 0d0      ; DATAin%Croots_stock_unc = 0d0      ; DATAin%Croots_stock_lag = 0
    DATAin%Clit_stock = 0d0        ; DATAin%Clit_stock_unc = 0d0        ; DATAin%Clit_stock_lag = 0
    DATAin%Csom_stock = 0d0        ; DATAin%Csom_stock_unc = 0d0        ; DATAin%Csom_stock_lag = 0
    DATAin%Cagb_stock = 0d0        ; DATAin%Cagb_stock_unc = 0d0        ; DATAin%Cagb_stock_lag = 0
    DATAin%Ccoarseroot_stock = 0d0 ; DATAin%Ccoarseroot_stock_unc = 0d0 ; DATAin%Ccoarseroot_stock_lag = 0
    DATAin%Evap = 0d0              ; DATAin%Evap_unc = 0d0              ; DATAin%Evap_lag = 0
    DATAin%SWE = 0d0               ; DATAin%SWE_unc = 0d0               ; DATAin%SWE_lag = 0
    DATAin%NBE = 0d0               ; DATAin%NBE_unc = 0d0               ; DATAin%NBE_lag = 0
    DATAin%Fire = 0d0              ; DATAin%Fire_unc = 0d0              ; DATAin%Fire_lag = 0
    DATAin%fAPAR = 0d0             ; DATAin%fAPAR_unc = 0d0             ; DATAin%fAPAR_lag = 0
    DATAin%Cwood_inc = 0d0         ; DATAin%Cwood_inc_unc = 0d0         ; DATAin%Cwood_inc_lag = 0
    DATAin%Cwood_growth = 0d0      ; DATAin%Cwood_growth_unc = 0d0      ; DATAin%Cwood_growth_lag = 0
    DATAin%Cwood_mortality = 0d0   ; DATAin%Cwood_mortality_unc = 0d0   ; DATAin%Cwood_mortality_lag = 0
    DATAin%harvest = 0d0           ; DATAin%harvest_unc = 0d0           ; DATAin%harvest_lag = 0
    DATAin%foliage_to_litter = 0d0 ; DATAin%foliage_to_litter_unc = 0d0 ; DATAin%foliage_to_litter_lag = 0
    DATAin%soilwater = 0d0         ; DATAin%soilwater_unc = 0d0         ; DATAin%soilwater_lag = 0
    ! Temorary arrays
    mettemp = 0d0 ; obstemp = 0d0

      ! zero the obs counters
      DATAin%total_obs = 0
      DATAin%ngpp = 0
      DATAin%nlai = 0
      DATAin%nnee = 0
      DATAin%nCwood_inc = 0
      DATAin%nCwood_growth = 0
      DATAin%nCwood_mortality = 0
      DATAin%nfoliage_to_litter = 0
      DATAin%nreco = 0
      DATAin%nCfol_stock = 0
      DATAin%nCwood_stock = 0
      DATAin%nCroots_stock = 0
      DATAin%nCsom_stock = 0
      DATAin%nClit_stock = 0
      DATAin%nCagb_stock = 0
      DATAin%nCcoarseroot_stock = 0
      DATAin%nEvap = 0
      DATAin%nSWE = 0
      DATAin%nNBE = 0
      DATAin%nFire = 0
      DATAin%nfAPAR = 0
      DATAin%nharvest = 0
      DATAin%nsoilwater = 0

      ! work out some key variables
      ! DATAin%noobs corresponds to observations and uncertainties
      totcol = DATAin%nomet*DATAin%noobs
      totread = 500 + 1

      ! start looping through days and allocate the correct met drivers/obs into
      ! the correct arrays
      do day = 1, DATAin%nodays
         start = ((day - 1)*(totcol)) + totread
         finish = start + DATAin%nomet - 1
         b = 1
         do i = start, finish
            read (ifile_unit) mettemp(b); b = b + 1
         end do  ! met bit
         start = ((day - 1)*(totcol)) + totread + finish + 1
         finish = start + DATAin%noobs - 1
         b = 1
         do i = start, finish
            read (ifile_unit) obstemp(b); b = b + 1
         end do  ! obs bit

         ! assign the extracted met/obs to their type and keep count of how many
         ! of these are actually contain data
         DATAin%met(1:DATAin%nomet, day) = mettemp

         ! Reset counter for extracting observations
         a = 1

         ! Gross Primary Productivity (GPP, gC/m2/day)
         DATAin%GPP(day) = obstemp(a); a = a + 1
         if (obstemp(a - 1) > -9998d0) DATAin%ngpp = DATAin%ngpp + 1
         DATAin%GPP_unc(day) = obstemp(a); a = a + 1
         DATAin%GPP_lag(day) = nint(obstemp(a)); a = a + 1

         ! Leaf Area Index (LAI, m2/m2)
         DATAin%LAI(day) = obstemp(a); a = a + 1
         if (obstemp(a - 1) > -9998d0) DATAin%nlai = DATAin%nlai + 1
         DATAin%LAI_unc(day) = obstemp(a); a = a + 1
         DATAin%LAI_lag(day) = nint(obstemp(a)); a = a + 1

         ! Net Ecosystem Exchange (NEE) of CO2 (gC/m2/day)
         DATAin%NEE(day) = obstemp(a); a = a + 1
         if (obstemp(a - 1) > -9998d0) DATAin%nnee = DATAin%nnee + 1
         DATAin%NEE_unc(day) = obstemp(a); a = a + 1
         DATAin%NEE_lag(day) = nint(obstemp(a)); a = a + 1

         ! Fire emissions of C (gC/m2day)
         DATAin%Fire(day) = obstemp(a); a = a + 1
         if (obstemp(a - 1) > -9998d0) DATAin%nFire = DATAin%nFire + 1
         DATAin%Fire_unc(day) = obstemp(a); a = a + 1
         DATAin%Fire_lag(day) = nint(obstemp(a)); a = a + 1

         ! Ecosystem respiration (Reco, gC/m2/day)
         DATAin%Reco(day) = obstemp(a); a = a + 1
         if (obstemp(a - 1) > -9998d0) DATAin%nreco = DATAin%nreco + 1
         DATAin%Reco_unc(day) = obstemp(a); a = a + 1
         DATAin%Reco_lag(day) = nint(obstemp(a)); a = a + 1

         ! C storage in foliage (gC/m2)
         DATAin%Cfol_stock(day) = obstemp(a); a = a + 1
         if (obstemp(a - 1) > -9998d0) DATAin%nCfol_stock = DATAin%nCfol_stock + 1
         DATAin%Cfol_stock_unc(day) = obstemp(a); a = a + 1
         DATAin%Cfol_stock_lag(day) = nint(obstemp(a)); a = a + 1

         ! C storage in wood (above+below) (gC/m2)
         DATAin%Cwood_stock(day) = obstemp(a); a = a + 1
         if (obstemp(a - 1) > -9998d0) DATAin%nCwood_stock = DATAin%nCwood_stock + 1
         DATAin%Cwood_stock_unc(day) = obstemp(a); a = a + 1
         DATAin%Cwood_stock_lag(day) = nint(obstemp(a)); a = a + 1

         ! C storage in fine roots (gC/m2)
         DATAin%Croots_stock(day) = obstemp(a); a = a + 1
         if (obstemp(a - 1) > -9998d0) DATAin%nCroots_stock = DATAin%nCroots_stock + 1
         DATAin%Croots_stock_unc(day) = obstemp(a); a = a + 1
         DATAin%Croots_stock_lag(day) = nint(obstemp(a)); a = a + 1

         ! C storage in fine litter (foliar+fine root, gC/m2)
         DATAin%Clit_stock(day) = obstemp(a); a = a + 1
         if (obstemp(a - 1) > -9998d0) DATAin%nClit_stock = DATAin%nClit_stock + 1
         DATAin%Clit_stock_unc(day) = obstemp(a); a = a + 1
         DATAin%Clit_stock_lag(day) = nint(obstemp(a)); a = a + 1

         ! C storage in soil organic matter (gC/m2)
         DATAin%Csom_stock(day) = obstemp(a); a = a + 1
         if (obstemp(a - 1) > -9998d0) DATAin%nCsom_stock = DATAin%nCsom_stock + 1
         DATAin%Csom_stock_unc(day) = obstemp(a); a = a + 1
         DATAin%Csom_stock_lag(day) = nint(obstemp(a)); a = a + 1

         ! C storage in above gound woody biomass (gC/m2)
         DATAin%Cagb_stock(day) = obstemp(a); a = a + 1
         if (obstemp(a - 1) > -9998d0) DATAin%nCagb_stock = DATAin%nCagb_stock + 1
         DATAin%Cagb_stock_unc(day) = obstemp(a); a = a + 1
         DATAin%Cagb_stock_lag(day) = nint(obstemp(a)); a = a + 1

         ! Fraction of absorbed photosynthetically active radiation
         ! by green vegetation (0-1)
         DATAin%fAPAR(day) = obstemp(a); a = a + 1
         if (obstemp(a - 1) > -9998d0) DATAin%nfAPAR = DATAin%nfAPAR + 1
         DATAin%fAPAR_unc(day) = obstemp(a); a = a + 1
         DATAin%fAPAR_lag(day) = nint(obstemp(a)); a = a + 1

         ! C storage in coarse root, i.e. below ground wood (gC/m2)
         DATAin%Ccoarseroot_stock(day) = obstemp(a); a = a + 1
         if (obstemp(a - 1) > -9998d0) DATAin%nCcoarseroot_stock = DATAin%nCcoarseroot_stock + 1
         DATAin%Ccoarseroot_stock_unc(day) = obstemp(a); a = a + 1
         DATAin%Ccoarseroot_stock_lag(day) = nint(obstemp(a)); a = a + 1

         ! Evapotranspiration (kgH2O/m2/day)
         DATAin%Evap(day) = obstemp(a); a = a + 1
         if (obstemp(a - 1) > -9998d0) DATAin%nEvap = DATAin%nEvap + 1
         DATAin%Evap_unc(day) = obstemp(a); a = a + 1
         DATAin%Evap_lag(day) = nint(obstemp(a)); a = a + 1

         ! Snow water equivalent-added for future use, not currently coded
         DATAin%SWE(day) = obstemp(a); a = a + 1
         if (obstemp(a - 1) > -9998d0) DATAin%nSWE = DATAin%nSWE + 1
         DATAin%SWE_unc(day) = obstemp(a); a = a + 1
         DATAin%SWE_lag(day) = nint(obstemp(a)); a = a + 1

         ! Net Biome Exchange (NBE, gC/m2/day)
         ! NBE = Reco+Fire-GPP
         DATAin%NBE(day) = obstemp(a); a = a + 1
         if (obstemp(a - 1) > -9998d0) DATAin%nNBE = DATAin%nNBE + 1
         DATAin%NBE_unc(day) = obstemp(a); a = a + 1
         DATAin%NBE_lag(day) = nint(obstemp(a)); a = a + 1

         ! Woody gross increment (gC/m2/day)
         ! Represents the average across the lagged period
         DATAin%Cwood_growth(day) = obstemp(a); a = a + 1
         if (obstemp(a - 1) > -9998d0) DATAin%nCwood_growth = DATAin%nCwood_growth + 1
         DATAin%Cwood_growth_unc(day) = obstemp(a); a = a + 1
         DATAin%Cwood_growth_lag(day) = nint(obstemp(a)); a = a + 1

         ! Woody natural mortality (gC/m2/day)
         ! Represents the average across the lagged period
         DATAin%Cwood_mortality(day) = obstemp(a); a = a + 1
         if (obstemp(a - 1) > -9998d0) DATAin%nCwood_mortality = DATAin%nCwood_mortality + 1
         DATAin%Cwood_mortality_unc(day) = obstemp(a); a = a + 1
         DATAin%Cwood_mortality_lag(day) = nint(obstemp(a)); a = a + 1

         ! Flux from foliage to litter (gC/m2/day)
         ! Represents the average across the lagged period
         DATAin%foliage_to_litter(day) = obstemp(a); a = a + 1
         if (obstemp(a - 1) > -9998d0) DATAin%nfoliage_to_litter = DATAin%nfoliage_to_litter + 1
         DATAin%foliage_to_litter_unc(day) = obstemp(a); a = a + 1
         DATAin%foliage_to_litter_lag(day) = nint(obstemp(a)); a = a + 1

         ! Woody net increment (gC/m2/day)
         ! i.e growth-mortality
         DATAin%Cwood_inc(day) = obstemp(a); a = a + 1
         if (obstemp(a - 1) > -9998d0) DATAin%nCwood_inc = DATAin%nCwood_inc + 1
         DATAin%Cwood_inc_unc(day) = obstemp(a); a = a + 1
         DATAin%Cwood_inc_lag(day) = nint(obstemp(a)); a = a + 1

         ! C extracted due to harvest (gC/m2/day)
         ! Represents the average across the lagged period
         DATAin%harvest(day) = obstemp(a); a = a + 1
         if (obstemp(a - 1) > -9998d0) DATAin%nharvest = DATAin%nharvest + 1
         DATAin%harvest_unc(day) = obstemp(a); a = a + 1
         DATAin%harvest_lag(day) = nint(obstemp(a)); a = a + 1

         ! Soil surface water content (m3/m3)
         DATAin%soilwater(day) = obstemp(a); a = a + 1
         if (obstemp(a - 1) > -9998d0) DATAin%nsoilwater = DATAin%nsoilwater + 1
         DATAin%soilwater_unc(day) = obstemp(a); a = a + 1
         DATAin%soilwater_lag(day) = nint(obstemp(a)); a = a + 1

      end do  ! day loop

      ! Count the total number of observations which are to be used.
      ! This total in some models may be used to inform on a dynamic weighting of the EDCs
      DATAin%total_obs = DATAin%ngpp + DATAin%nlai + DATAin%nnee &
                         + DATAin%nCwood_inc + DATAin%nreco + DATAin%nCfol_stock &
                         + DATAin%nCwood_stock + DATAin%nCroots_stock + DATAin%nCsom_stock &
                         + DATAin%nClit_stock + DATAin%nCagb_stock + DATAin%nCcoarseroot_stock &
                         + DATAin%nEvap + DATAin%nSWE + DATAin%nNBE &
                         + DATAin%nCwood_mortality + DATAin%nfoliage_to_litter + DATAin%nFire &
                         + DATAin%nfAPAR + DATAin%nharvest + DATAin%nsoilwater + DATAin%nCwood_growth

      ! allocate to time step
      allocate (DATAin%deltat(DATAin%nodays)); DATAin%deltat = 0d0
      ! work out interval in decimal days
      DATAin%deltat(1) = DATAin%met(1, 1)
      do i = 2, DATAin%nodays
         DATAin%deltat(i) = DATAin%met(1, i) - DATAin%met(1, (i - 1))
      end do
      ! Calculate the number of years being simulated
      DATAin%nos_years = nint(sum(DATAin%deltat)/365.25d0)
      ! Calculate the number of steps per year
      ! NOTE: this should be doable as integer as there should be an exact division
      ! as all simulations are conducted in whole years
      DATAin%steps_per_year = DATAin%nodays/DATAin%nos_years

      ! close open binary
      close (ifile_unit)

      ! clean up memory
      deallocate (obstemp, mettemp)

      ! allocate memory to our data location variables
      if (DATAin%ngpp > 0) allocate (DATAin%gpppts(DATAin%ngpp))
      if (DATAin%nlai > 0) allocate (DATAin%laipts(DATAin%nlai))
      if (DATAin%nnee > 0) allocate (DATAin%neepts(DATAin%nnee))
      if (DATAin%nreco > 0) allocate (DATAin%recopts(DATAin%nreco))
      if (DATAin%nCfol_stock > 0) allocate (DATAin%Cfol_stockpts(DATAin%nCfol_stock))
      if (DATAin%nCwood_stock > 0) allocate (DATAin%Cwood_stockpts(DATAin%nCwood_stock))
      if (DATAin%nCroots_stock > 0) allocate (DATAin%Croots_stockpts(DATAin%nCroots_stock))
      if (DATAin%nClit_stock > 0) allocate (DATAin%Clit_stockpts(DATAin%nClit_stock))
      if (DATAin%nCsom_stock > 0) allocate (DATAin%Csom_stockpts(DATAin%nCsom_stock))
      if (DATAin%nCagb_stock > 0) allocate (DATAin%Cagb_stockpts(DATAin%nCagb_stock))
      if (DATAin%nCcoarseroot_stock > 0) allocate (DATAin%Ccoarseroot_stockpts(DATAin%nCcoarseroot_stock))
      if (DATAin%nEvap > 0) allocate (DATAin%Evappts(DATAin%nEvap))
      if (DATAin%nSWE > 0) allocate (DATAin%SWEpts(DATAin%nSWE))
      if (DATAin%nNBE > 0) allocate (DATAin%NBEpts(DATAin%nNBE))
      if (DATAin%nCwood_growth > 0) allocate(DATAin%Cwood_growthpts(DATAin%nCwood_growth))
      if (DATAin%nCwood_inc > 0) allocate (DATAin%Cwood_incpts(DATAin%nCwood_inc))
      if (DATAin%nCwood_mortality > 0) allocate (DATAin%Cwood_mortalitypts(DATAin%nCwood_mortality))
      if (DATAin%nfoliage_to_litter > 0) allocate (DATAin%foliage_to_litterpts(DATAin%nfoliage_to_litter))
      if (DATAin%nFire > 0) allocate (DATAin%Firepts(DATAin%nFire))
      if (DATAin%nfAPAR > 0) allocate (DATAin%fAPARpts(DATAin%nfAPAR))
      if (DATAin%nharvest > 0) allocate (DATAin%harvestpts(DATAin%nharvest))
      if (DATAin%nsoilwater > 0) allocate (DATAin%soilwaterpts(DATAin%nsoilwater))

      ! we know how many observations we have and what they are, but now lets work
      ! out where they are in the data sets
      a = 1; b = 1; c = 1; d = 1; e = 1; f = 1; g = 1
      h = 1; i = 1; j = 1; k = 1; l = 1; m = 1; n = 1
      o = 1; p = 1; q = 1; r = 1; s = 1; t = 1; u = 1
      v = 1; w = 1;
      ! Read through each timestep to extract any available assimilatable observations
      do day = 1, DATAin%nodays
         if (DATAin%GPP(day) > -9998d0) then
            DATAin%gpppts(a) = day; a = a + 1
         end if
         if (DATAin%LAI(day) > -9998d0) then
            DATAin%laipts(b) = day; b = b + 1
         end if
         if (DATAin%NEE(day) > -9998d0) then
            DATAin%neepts(c) = day; c = c + 1
         end if
         if (DATAin%Cwood_inc(day) > -9998d0) then
            DATAin%Cwood_incpts(d) = day; d = d + 1
         end if  ! data present condition
         if (DATAin%Cwood_mortality(day) > -9998d0) then
            DATAin%Cwood_mortalitypts(e) = day; e = e + 1
         end if  ! data present condition
         if (DATAin%foliage_to_litter(day) > -9998d0) then
            DATAin%foliage_to_litterpts(f) = day; f = f + 1
         end if  ! data present condition
         if (DATAin%Reco(day) > -9998d0) then
            DATAin%recopts(g) = day; g = g + 1
         end if  ! data present condition
         if (DATAin%Cfol_stock(day) > -9998d0) then
            DATAin%Cfol_stockpts(h) = day; h = h + 1
         end if  ! data present condition
         if (DATAin%Cwood_stock(day) > -9998d0) then
            DATAin%Cwood_stockpts(i) = day; i = i + 1
         end if  ! data present condition
         if (DATAin%Croots_stock(day) > -9998d0) then
            DATAin%Croots_stockpts(j) = day; j = j + 1
         end if  ! data present condition
         if (DATAin%Clit_stock(day) > -9998d0) then
            DATAin%Clit_stockpts(k) = day; k = k + 1
         end if  ! data present condition
         if (DATAin%Csom_stock(day) > -9998d0) then
            DATAin%Csom_stockpts(l) = day; l = l + 1
         end if  ! data present condition
         if (DATAin%Ccoarseroot_stock(day) > -9998d0) then
            DATAin%Ccoarseroot_stockpts(m) = day; m = m + 1
         end if  ! data present condition
         if (DATAin%Cagb_stock(day) > -9998d0) then
            DATAin%Cagb_stockpts(n) = day; n = n + 1
         end if  ! data present condition
         if (DATAin%Evap(day) > -9998d0) then
            DATAin%Evappts(o) = day; o = o + 1
         end if  ! data present condition
         if (DATAin%SWE(day) > -9998d0) then
            DATAin%SWEpts(p) = day; p = p + 1
         end if  ! data present condition
         if (DATAin%NBE(day) > -9998d0) then
            DATAin%NBEpts(q) = day; q = q + 1
         end if
         if (DATAin%Fire(day) > -9998d0) then
            DATAin%Firepts(r) = day; r = r + 1
         end if  ! data present condition
         if (DATAin%fAPAR(day) > -9998d0) then
            DATAin%fAPARpts(s) = day; s = s + 1
         end if  ! data present condition
         if (DATAin%harvest(day) > -9998d0) then
            DATAin%harvestpts(t) = day; t = t + 1
         end if  ! data present condition
         if (DATAin%Cwood_growth(day) > -9998d0) then
            DATAin%Cwood_growthpts(u) = day; u = u + 1
         end if  ! data present condition
         if (DATAin%soilwater(day) > -9998d0) then
            DATAin%soilwaterpts(v) = day; v = v + 1
         end if  ! data present condition
      end do  ! day loop

      ! timestep mean temperature (oC)
      DATAin%meantemp = sum((DATAin%met(2, :) + DATAin%met(3, :))*0.5d0)/dble(DATAin%nodays)
      ! mean SW radiation (MJ/m2/day)
      DATAin%meanrad = sum(DATAin%met(4, :))/dble(DATAin%nodays)
      ! mean atmospheric CO2 (ppm)
      DATAin%meanco2 = sum(DATAin%met(5, :))/dble(DATAin%nodays)
      ! mean precipitation (mm/yr)
      DATAin%meanprecip = sum(DATAin%met(7, :)*86400d0*365.25d0)/dble(DATAin%nodays)

      ! print the mean temperature and radiation variables
      write (*, *) "Mean Rad (MJ/m2/day) = ", DATAin%meanrad
      write (*, *) "Mean Temp (Celcius) = ", DATAin%meantemp
      write (*, *) "Mean Precip (mm/yr) = ", DATAin%meanprecip
      write (*, *) "Mean CO2 (ppm) = ", DATAin%meanco2
      write (*, *) "==========="
      write (*, *) "Number of timesteps = ", DATAin%nodays
      write (*, *) "Total number of obs = ", DATAin%total_obs

      ! all done
      write (*, *) "Binary input file has been successfully read by CARDAMOM"

   end subroutine read_binary_data
   !
   !------------------------------------------------------------------

   subroutine initialize(infile)  ! formerly read_pari_data
      ! 3 steps must be called in this order
      use cardamom_structures, only: DATA_type, set_datain, set_datain_original
      use model_shared, only: initialize_parinfo
      implicit none(type, external)
      character(len=*), intent(in) :: infile
      type(DATA_type) :: DATAin  ! tmp datain object to collect all data before saving to cardamom_structures :: DATAin

      call initialize_parinfo()  ! TODO not really a file reading thing
      call read_check_binary_data(infile, DATAin)
      call initialize_model(DATAin)
      ! Save the DATAin from file to cardamom_structures : DATAin_original
      call set_datain_original(DATAin)
      ! Save the DATAin from file to cardamom_structures : DATAin (copy to be scaled)
      call set_datain(DATAin)
   end subroutine initialize

   subroutine read_check_binary_data(infile, DATAin)
  !! Read infile, modify fields of a (local) DATA_type struct
  !! subroutine call for input binary to be read and then allocates the input
  !! data to extracting the data to the correct observation and parameter types
  !! split from read_pari_data
      use MODEL_PARAMETERS, only: pars_info
      use model_shared, only: PI
      use cardamom_structures, only: DATA_type

      implicit none(type, external)

      type(DATA_type), intent(inout) :: DATAin

      character(350), intent(in) :: infile

      ! declare local variables
      integer :: i

      ! remind us what file we're about to access
      write (*, *) "Input file = ", trim(infile)
      ! initialise data structure and read the binary
      call read_binary_data(infile, DATAin)
      ! check:
      ! PI%npars = DATAin%nopars don't set from DATAin, instead check
      if (PI%npars /= DATAin%nopars) then
         write (*, *) "ERROR nopars from infile (", DATAin%nopars, ") does not &
         & match npars from model _PARS file (", PI%npars, ")"
         stop
      end if
      ! also check initial state:
      ! ensure any loaded parameter values (from the input file)
      ! are within the uniform parameter bounds set in the source code
      do i = 1, PI%npars
         if (DATAin%parpriors(i) /= -9999) then
            if (DATAin%parpriors(i) > PI%parmax(i) .or. DATAin%parpriors(i) < PI%parmin(i)) then
               write (*, *) PI%parmin(i)
               write (*, *) DATAin%parpriors(i)
               write (*, *) PI%parmax(i)
               write (*, *) "Supplied parameter prior = ", i, " is outside hardcoded uniform parameter bounds"
               stop
            end if
         end if
      end do

   end subroutine read_check_binary_data

   subroutine initialize_model(DATAin)
  !! split from read_pari_data
  !! TODO not io, belongs in a differnt file
  !! depends on having called read_binary_data or read_check_binary_data first
  !! for DATAin%nodays and nopools
      use cardamom_structures, only: DATA_type

      implicit none(type, external)

      type(DATA_type), intent(inout) :: DATAin
      ! need to allocate memory to the model output variables
      allocate (DATAin%M_FLUXES(DATAin%nodays, DATAin%nofluxes) &
                , DATAin%M_POOLS((DATAin%nodays + 1), DATAin%nopools), DATAin%M_DIAGS(DATAin%nodays, DATAin%nodiags))

      ! force zero in states and fluxes
      !DATAin%M_FLUXES(:,:) = 0d0; DATAin%M_POOLS(:,:) = 0d0; DATAin%M_DIAGS(:,:) = 0d0

      ! alert the user
      write (*, *) "Created fields for model output"
   end subroutine initialize_model
   !------------------------------------------------------------------
   !
   subroutine read_options(solutions_wanted, freq_print, freq_write, outfile, MCO)
      use cardamom_MHMCMC, only: MCMC_OPTIONS
      use samplers_shared, only: SAMPLER_OPTIONS

      ! loads required options about the MHMCMC either form hardcoded sources or
      ! from variables which were read form the command line

      implicit none(type, external)

      ! declare input variables
      character(350), intent(inout) :: outfile
      integer, intent(in) :: solutions_wanted, freq_print, freq_write
      class(SAMPLER_OPTIONS), intent(inout) :: MCO

      ! defining hardcoded MCMC options
      MCO%append = .true.
      MCO%nADAPT = 1000
      MCO%fADAPT = 0.5d0
      MCO%randparini = .false.
      MCO%returnpars = .false.
      MCO%fixedpars = .true. ! TLS: changed from .false. for testing 16/12/2019

      ! command line options

      ! how many accepted parameters to accept before completion
      if (solutions_wanted > 0) then
         MCO%nOUT = solutions_wanted
      else
         MCO%nOUT = 1000
         write (*, *) "Default MCO%nOUT value used"
      end if
      ! how frequently to print information to screen
      if (freq_print >= 0) then
         MCO%nPRINT = freq_print
      else
         MCO%nPRINT = 1000
         write (*, *) "Default MCO%nPRINT value used"
      end if
      ! how frequently to write information to file
      if (freq_write >= 0) then
         MCO%nWRITE = freq_write
      else
         MCO%nWRITE = 1000
         write (*, *) "Default MCO%nWRITE value used"
      end if

      ! Assume that sub-sampling process, if completed, will use 10 % of the
      ! simulation time therefore we want to adjust the output frequency to
      ! correct for this
      ! MCO%nOUT = max(1, MCO%nOUT-MCOUT%nos_iterations)

      ! construct file names
      ! ensure that the separator is _ , whether the user put it in the stem or not
      if (outfile(len_trim(outfile):len_trim(outfile)) /= "_" ) then
          outfile = trim(outfile)//"_"
      end if
      write (MCO%outfile, fmt='(A)') trim(outfile)//"PARS"
      write (MCO%stepfile, fmt='(A)') trim(outfile)//"STEP"
      write (MCO%covfile, fmt='(A)') trim(outfile)//"COV"
      write (MCO%covifile, fmt='(A)') trim(outfile)//"COVINFO"

   end subroutine read_options
   !
   !-------------------------------------------------------------------
   !
   !
   !------------------------------------------------------------------
   !
   subroutine update_obs_scaling_normal
      use cardamom_structures, only: DATA_type, DATAin_original, set_datain
      type(DATA_type) :: DATAin_tmp  ! we edit a local tmp copy, then write it back to shared storage location
      ! it's ok to make changes to elements of DATAin in single-threaded parts of the program
      ! - but I made it requires more deliberate steps to stop accidental updates
      ! from concurrent parts
      DATAin_tmp = DATAin_original

      ! Subroutine sets the data specific scaling factors
      ! in this case set all to one where the weight of each
      ! observation is equal

      DATAin_tmp%GPP_scaling = 1d0
      DATAin_tmp%NEE_scaling = 1d0
      DATAin_tmp%Fire_scaling = 1d0
      DATAin_tmp%LAI_scaling = 1d0
      DATAin_tmp%Cwood_inc_scaling = 1d0
      DATAin_tmp%Cwood_growth_scaling = 1d0
      DATAin_tmp%Cwood_mortality_scaling = 1d0
      DATAin_tmp%foliage_to_litter_scaling = 1d0
      DATAin_tmp%Reco_scaling = 1d0
      DATAin_tmp%Cfol_stock_scaling = 1d0
      DATAin_tmp%Cwood_stock_scaling = 1d0
      DATAin_tmp%Croots_stock_scaling = 1d0
      DATAin_tmp%Csom_stock_scaling = 1d0
      DATAin_tmp%Cagb_stock_scaling = 1d0
      DATAin_tmp%Clit_stock_scaling = 1d0
      DATAin_tmp%Ccoarseroot_stock_scaling = 1d0
      DATAin_tmp%Evap_scaling = 1d0
      DATAin_tmp%SWE_scaling = 1d0
      DATAin_tmp%NBE_scaling = 1d0
      DATAin_tmp%fAPAR_scaling = 1d0
      DATAin_tmp%harvest_scaling = 1d0
      DATAin_tmp%soilwater_scaling = 1d0

      call set_datain(DATAin_tmp)  ! copy the tmp object to shared cardamom_structures DATAin

      return

   end subroutine update_obs_scaling_normal
   !
   !------------------------------------------------------------------
   !
   subroutine update_obs_scaling_nsamples
      use cardamom_structures, only: DATA_type, DATAin_original, set_datain
      type(DATA_type) :: DATAin_tmp
      DATAin_tmp = DATAin_original

      ! Subroutine sets the data specific scaling factors
      ! in this case are all normalised by the sample size
      ! giving equal weight to each datastream

      DATAin_tmp%GPP_scaling = 1d0/dble(DATAin_original%ngpp)
      DATAin_tmp%NEE_scaling = 1d0/dble(DATAin_original%nnee)
      DATAin_tmp%Fire_scaling = 1d0/dble(DATAin_original%nfire)
      DATAin_tmp%LAI_scaling = 1d0/dble(DATAin_original%nlai)
      DATAin_tmp%Cwood_inc_scaling = 1d0/dble(DATAin_original%nCwood_inc)
      DATAin_tmp%Cwood_growth_scaling = 1d0/dble(DATAin_original%nCwood_growth)
      DATAin_tmp%Cwood_mortality_scaling = 1d0/dble(DATAin_original%nCwood_mortality)
      DATAin_tmp%foliage_to_litter_scaling = 1d0/dble(DATAin_original%nfoliage_to_litter)
      DATAin_tmp%Reco_scaling = 1d0/dble(DATAin_original%nreco)
      DATAin_tmp%Cfol_stock_scaling = 1d0/dble(DATAin_original%nCfol_stock)
      DATAin_tmp%Cwood_stock_scaling = 1d0/dble(DATAin_original%nCwood_stock)
      DATAin_tmp%Croots_stock_scaling = 1d0/dble(DATAin_original%nCroots_stock)
      DATAin_tmp%Csom_stock_scaling = 1d0/dble(DATAin_original%nCsom_stock)
      DATAin_tmp%Cagb_stock_scaling = 1d0/dble(DATAin_original%nCagb_stock)
      DATAin_tmp%Clit_stock_scaling = 1d0/dble(DATAin_original%nClit_stock)
      DATAin_tmp%Ccoarseroot_stock_scaling = 1d0/dble(DATAin_original%nCcoarseroot_stock)
      DATAin_tmp%Evap_scaling = 1d0/dble(DATAin_original%nEvap)
      DATAin_tmp%SWE_scaling = 1d0/dble(DATAin_original%nSWE)
      DATAin_tmp%NBE_scaling = 1d0/dble(DATAin_original%nnbe)
      DATAin_tmp%fAPAR_scaling = 1d0/dble(DATAin_original%nfAPAR)
      DATAin_tmp%harvest_scaling = 1d0/dble(DATAin_original%nharvest)
      DATAin_tmp%soilwater_scaling = 1d0/dble(DATAin_original%nsoilwater)

      call set_datain(DATAin_tmp)

      return

   end subroutine update_obs_scaling_nsamples
   !
   !------------------------------------------------------------------
   !
   subroutine update_obs_scaling_sqrt_nsamples
      use cardamom_structures, only: DATA_type, DATAin_original, set_datain
      type(DATA_type) :: DATAin_tmp
      DATAin_tmp = DATAin_original

      ! Subroutine sets the data specific scaling factors
      ! in this case are all normalised by the sqrt of the
      ! sample size. This allows datastreams with more observations
      ! to contribute more to the lost function but penalised to reduce
      ! bias' introduced by unbalanced observations

      DATAin_tmp%GPP_scaling = 1d0/sqrt(dble(DATAin_original%ngpp))
      DATAin_tmp%NEE_scaling = 1d0/sqrt(dble(DATAin_original%nnee))
      DATAin_tmp%Fire_scaling = 1d0/sqrt(dble(DATAin_original%nfire))
      DATAin_tmp%LAI_scaling = 1d0/sqrt(dble(DATAin_original%nlai))
      DATAin_tmp%Cwood_inc_scaling = 1d0/sqrt(dble(DATAin_original%nCwood_inc))
      DATAin_tmp%Cwood_growth_scaling = 1d0/sqrt(dble(DATAin_original%nCwood_growth))
      DATAin_tmp%Cwood_mortality_scaling = 1d0/sqrt(dble(DATAin_original%nCwood_mortality))
      DATAin_tmp%foliage_to_litter_scaling = 1d0/sqrt(dble(DATAin_original%nfoliage_to_litter))
      DATAin_tmp%Reco_scaling = 1d0/sqrt(dble(DATAin_original%nreco))
      DATAin_tmp%Cfol_stock_scaling = 1d0/sqrt(dble(DATAin_original%nCfol_stock))
      DATAin_tmp%Cwood_stock_scaling = 1d0/sqrt(dble(DATAin_original%nCwood_stock))
      DATAin_tmp%Croots_stock_scaling = 1d0/sqrt(dble(DATAin_original%nCroots_stock))
      DATAin_tmp%Csom_stock_scaling = 1d0/sqrt(dble(DATAin_original%nCsom_stock))
      DATAin_tmp%Cagb_stock_scaling = 1d0/sqrt(dble(DATAin_original%nCagb_stock))
      DATAin_tmp%Clit_stock_scaling = 1d0/sqrt(dble(DATAin_original%nClit_stock))
      DATAin_tmp%Ccoarseroot_stock_scaling = 1d0/sqrt(dble(DATAin_original%nCcoarseroot_stock))
      DATAin_tmp%Evap_scaling = 1d0/sqrt(dble(DATAin_original%nEvap))
      DATAin_tmp%SWE_scaling = 1d0/sqrt(dble(DATAin_original%nSWE))
      DATAin_tmp%NBE_scaling = 1d0/sqrt(dble(DATAin_original%nnbe))
      DATAin_tmp%fAPAR_scaling = 1d0/sqrt(dble(DATAin_original%nfAPAR))
      DATAin_tmp%harvest_scaling = 1d0/sqrt(dble(DATAin_original%nharvest))
      DATAin_tmp%soilwater_scaling = 1d0/sqrt(dble(DATAin_original%nsoilwater))

      call set_datain(DATAin_tmp)
      return

   end subroutine update_obs_scaling_sqrt_nsamples
   !
   !------------------------------------------------------------------
   !
   subroutine update_obs_scaling_log_nsamples
      use cardamom_structures, only: DATA_type, DATAin_original, set_datain
      type(DATA_type) :: DATAin_tmp
      DATAin_tmp = DATAin_original

      ! Subroutine sets the data specific scaling factors
      ! in this case are all normalised by the sqrt of the
      ! sample size. This allows datastreams with more observations
      ! to contribute more to the lost function but penalised to reduce
      ! bias' introduced by unbalanced observations

      DATAin_tmp%GPP_scaling = 1d0/(1d0 + log(dble(DATAin_original%ngpp)))
      DATAin_tmp%NEE_scaling = 1d0/(1d0 + log(dble(DATAin_original%nnee)))
      DATAin_tmp%Fire_scaling = 1d0/(1d0 + log(dble(DATAin_original%nfire)))
      DATAin_tmp%LAI_scaling = 1d0/(1d0 + log(dble(DATAin_original%nlai)))
      DATAin_tmp%Cwood_inc_scaling = 1d0/(1d0 + log(dble(DATAin_original%nCwood_inc)))
      DATAin_tmp%Cwood_growth_scaling = 1d0/(1d0 + log(dble(DATAin_original%nCwood_growth)))
      DATAin_tmp%Cwood_mortality_scaling = 1d0/(1d0 + log(dble(DATAin_original%nCwood_mortality)))
      DATAin_tmp%foliage_to_litter_scaling = 1d0/(1d0 + log(dble(DATAin_original%nfoliage_to_litter)))
      DATAin_tmp%Reco_scaling = 1d0/(1d0 + log(dble(DATAin_original%nreco)))
      DATAin_tmp%Cfol_stock_scaling = 1d0/(1d0 + log(dble(DATAin_original%nCfol_stock)))
      DATAin_tmp%Cwood_stock_scaling = 1d0/(1d0 + log(dble(DATAin_original%nCwood_stock)))
      DATAin_tmp%Croots_stock_scaling = 1d0/(1d0 + log(dble(DATAin_original%nCroots_stock)))
      DATAin_tmp%Csom_stock_scaling = 1d0/(1d0 + log(dble(DATAin_original%nCsom_stock)))
      DATAin_tmp%Cagb_stock_scaling = 1d0/(1d0 + log(dble(DATAin_original%nCagb_stock)))
      DATAin_tmp%Clit_stock_scaling = 1d0/(1d0 + log(dble(DATAin_original%nClit_stock)))
      DATAin_tmp%Ccoarseroot_stock_scaling = 1d0/(1d0 + log(dble(DATAin_original%nCcoarseroot_stock)))
      DATAin_tmp%Evap_scaling = 1d0/(1d0 + log(dble(DATAin_original%nEvap)))
      DATAin_tmp%SWE_scaling = 1d0/(1d0 + log(dble(DATAin_original%nSWE)))
      DATAin_tmp%NBE_scaling = 1d0/(1d0 + log(dble(DATAin_original%nnbe)))
      DATAin_tmp%fAPAR_scaling = 1d0/(1d0 + log(dble(DATAin_original%nfAPAR)))
      DATAin_tmp%harvest_scaling = 1d0/(1d0 + log(dble(DATAin_original%nharvest)))
      DATAin_tmp%soilwater_scaling = 1d0/(1d0 + log(dble(DATAin_original%nsoilwater)))

      call set_datain(DATAin_tmp)
      return

   end subroutine update_obs_scaling_log_nsamples

end module cardamom_io
