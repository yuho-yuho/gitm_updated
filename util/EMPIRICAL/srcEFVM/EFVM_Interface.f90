! Basic module for the EFVM
! Created: Qingyu Zhu, 08/16/2020
! -----------------------------------------------------------------------------

module efvm_interface

  implicit none 

  !!! PATH
  ! ---------------------------------------------------------------------------
  !! EFVM
  character(len=*), parameter :: &
       efvm_data_path='UA/DataIn/', &
       efvm_coeff_path = efvm_data_path//'dmsp_hld/efvm/coeffs/', &
       efvm_ca_fourier_coeff_path1='ca_fourier_coeffs/Ed1/',&
       efvm_ca_fourier_coeff_path2='ca_fourier_coeffs/Ed2/',&
       efvm_cat0_coeff_path = 'bin_cat0/'

  !!! GLOBAL PARAMETERS
  ! ---------------------------------------------------------------------------
  !! EFVM
  ! GRIDS 
  integer, parameter :: efvm_ncf = 6, efvm_nca = 8, &
       efvm_nmlt=24, efvm_nmlat=20

  ! FOURIER FITTING ORDER 
  integer, parameter :: efvm_LMAX=6 ! (MLT Fourier fitting)
  integer, parameter :: efvm_NMAX=4 ! (Ca Fourier fitting)
  logical, parameter :: efvm_debug=.false.

  ! OTHER PARAMETER
  real, parameter :: efvm_cf_inf=40000. ! (EPOT saturation)
  real, parameter :: efvm_cf0=2700.
  
  !!! IUNPUTS FROM GITM
  ! ---------------------------------------------------------------------------
  real :: IO_IMFBz = -1.0, IO_IMFBy = -1.0
  real :: IO_SWVX = 400., IO_SWN = 4 

  !!! GLOBAL VARIABLES 
  ! ---------------------------------------------------------------------------
  !! IMF & COUPLING FUNCTIONS 
  real :: IO_Bt=1., IO_Ca=0., IO_Cf=3000.  

  !! EFVM VARIABLES
  ! REFRENCE CFS, MLT & MLAT 
  real, allocatable, dimension(:) :: efvm_ref_cfs
  real, allocatable, dimension(:) :: efvm_mlts, efvm_mlats, efvm_mlats1
  
  ! EFVM COEFFICIENTS
  real, allocatable, dimension(:,:,:,:) :: efvm_all_ca_coeffs1, &
       efvm_all_ca_coeffs2
  real, allocatable, dimension(:,:) :: efvm_cat0_coeffs1, efvm_cat0_coeffs2

  ! DERIVED PARAMETERS
  real :: efvm_coeffs1(efvm_nmlat,(efvm_LMAX+1)**2), &
       efvm_coeffs2(efvm_nmlat,(efvm_LMAX+1)**2)
  real :: efvm_coeffs1_sh(efvm_nmlat, (efvm_LMAX+1)**2), &
       efvm_coeffs2_sh(efvm_nmlat,(efvm_LMAX+1)**2)

  ! OUTPUTS
  real, allocatable, dimension(:,:) :: efvm_dEd1, efvm_dEd2
  real, allocatable, dimension(:,:) :: efvm_dEd1_sh, efvm_dEd2_sh

  ! ===========================================================================
  contains

    !!! OBTAIN IMF BY, BZ, SWVX & SWN
    subroutine efvm_readin(by,bz,swvx,swn)

      implicit none 

      real, intent(in) :: by, bz, swvx, swn

      IO_IMFBy = by
      IO_IMFBz = bz
      IO_SWVX = abs(swvx)
      IO_SWN = swn

    end subroutine efvm_readin

end module efvm_interface
