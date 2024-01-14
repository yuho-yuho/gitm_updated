! This module provides all parameters needed for the AEPM module
! Updated from the ModSAM_Interface.f90
! Created: Qingyu Zhu, 05/17/2020
!
! -----------------------------------------------------------------------------
module aepm_interface

  implicit none 

  !!! PATH
  ! ---------------------------------------------------------------------------
  character(len=*), parameter :: &
       ! AEPM coefficient directory 
       aepm_data_path = 'UA/DataIn/', &
       aepm_coeff_path = &
       aepm_data_path//'dmsp_hld/aepm/coeffs/', &
       ! Directory for Cat 1 - N
       aepm_ca_fourier_coeff_path = 'ca_fourier_coeffs/', & 
       ! Directory for Cat 0
       aepm_cat0_coeff_path = 'mlt_fourier_coeffs/cat0/', &
       ! Slope and Y-int Fourier efficients 
       aepm_slope_yint_path = aepm_data_path//'dmsp_hld/aepm/slope_yint/'

  !!! GLOBAL PARAMETERS
  ! ---------------------------------------------------------------------------
  ! SIZES OF AEPM's GRIDS
  integer, parameter :: aepm_ncf = 8, aepm_nca = 8, aepm_nchannel = 19, &
       aepm_nmlt=24, aepm_nmlat=40

  ! FOURIER EXPANSION ORDERS FOR MLT AND IMF CLOCK ANGLE
  integer, parameter :: aepm_LMAX = 4, &  ! (MLT)
       aepm_NMAX = 4  ! (IMF clock angle)

  ! Other parameters 
  real, parameter :: aepm_cf_inf=20000.
  real, parameter :: aepm_cf0=2700.
  real, parameter :: pi = 3.14159265359 
  logical, parameter :: aepm_debug=.false.
  
  !!! INPUTS FROM GITM
  ! ---------------------------------------------------------------------------
  ! IMF and solar wind 
  real :: IO_IMFBz = -1.0, IO_IMFBy = -1.0
  real :: IO_SWVX = 400., IO_SWN = 4
  
  ! HP
  real :: IO_HP = 1.0
  logical :: IO_HP_scaling = .false.

  !!! AEPM PARMAMETERS
  ! ---------------------------------------------------------------------------
  ! IMF and coupling function 
  real :: IO_Bt=1., IO_Ca=0., IO_Cf=3000. 

  ! Reference Cfs 
  real, allocatable, dimension(:) :: aepm_ref_cfs
  ! AEPM grids
  real, allocatable, dimension(:) :: aepm_mlats, aepm_mlts, aepm_channels
  ! Cat 1-N coefficients
  real, allocatable, dimension(:,:,:,:,:) :: aepm_all_ca_coeffs
  ! Cat 0 coefficients 
  real, allocatable, dimension(:,:,:) :: aepm_cat0_coeffs
  ! Differential energy fluxes in both hemispheres 
  real, allocatable, dimension(:,:,:) :: aepm_diff_ef, aepm_diff_ef_sh
  ! Total energy flux and number fluxes in both hemispheres 
  real, allocatable, dimension(:,:) :: aepm_total_ef, aepm_total_nf
  real, allocatable, dimension(:,:) :: aepm_total_ef_sh, aepm_total_nf_sh
  ! Coeffs of slope and y-int
  real, allocatable, dimension(:,:) :: aepm_slope_coeffs, aepm_yint_coeffs 
  ! Slope and y-int
  real, allocatable, dimension(:) :: aepm_slope, aepm_yint
  ! Expansion rate and expanded mlats 
  real :: aepm_expansion_rate
  real, allocatable, dimension(:) :: aepm_mlats1

  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Functions that obtain parameters from GITM
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  contains 

    !!! OBTAIN IMF BY, BZ, SWVX, SWN
    ! -------------------------------------------------------------------------
    subroutine aepm_readin(by,bz,swvx,swn)

      implicit none 

      real, intent(in) :: by, bz, swvx, swn

      IO_IMFBy = by
      IO_IMFBz = bz
      IO_SWVX = abs(swvx)
      IO_SWN = swn

    end subroutine aepm_readin

    !!! READ ANY HP INPUTS (if scaling according to HP is needed)
    ! -------------------------------------------------------------------------
    subroutine aepm_readhp(hpin,hpscaling)

      implicit none 

      real, intent(in) :: hpin
      logical, intent(in) :: hpscaling

      IO_HP = hpin
      IO_HP_scaling = hpscaling 

    end subroutine aepm_readhp

end module aepm_interface
