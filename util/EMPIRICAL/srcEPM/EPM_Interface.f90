! This module provides all parameters needed for the EPM module 
! Updated from ModSEP_Interface.f90
! Created: Qingyu Zhu, 05/25/2020
!
! -----------------------------------------------------------------------------
module epm_interface

  implicit none 

  !!! PATH
  ! ---------------------------------------------------------------------------
  character(len=*), parameter :: &
       epm_data_path = 'UA/DataIn/', &
       epm_coeff_path = epm_data_path//'dmsp_hld/epm/coeffs/', &
       ! Directory for Cat 1 - N
       epm_ca_fourier_coeff_path = 'ca_fourier_coeffs/', &
       ! Directory for Cat 0
       epm_cat0_coeff_path = 'shf_parms_1/', &
       ! Slope and Y-int Fourier coefficiens 
       epm_slope_yint_path = epm_data_path//'dmsp_hld/epm/slope_yint/'

  !!! Global PARAMETERS 
  ! ---------------------------------------------------------------------------
  ! SIZES OF EPM'S GRIDS 
  integer, parameter :: epm_ncf = 5, epm_nca = 16, &
       epm_nmlt=24, epm_nmlat=60

  ! SPHERIC HARMONIC FUNCTION ORDER AND FOURIER ORDER
  integer, parameter :: epm_LMAX=12, &  ! (SPH order)
       epm_NMAX=8, & ! (Ca Fourier order)
       epm_NMAX1=4 ! (Ca Fourier order for Slope and yint)

  ! OTHER PARAMETERS 
  real, parameter :: epm_cf_inf=40000. ! (EPOT saturation)
  real, parameter :: epm_cf_inf1=24000. ! (CRB saturation)
  real, parameter :: pi = 3.14159265359
  real, parameter :: epm_cf0=2700.
  logical, parameter :: epm_debug=.false.

  !!! INPUTS FROM GITM
  ! ---------------------------------------------------------------------------
  real :: IO_IMFBz = -1.0, IO_IMFBy = -1.0                                    
  real :: IO_SWVX = 400., IO_SWN = 4

  !!! EPM PARAMETERS 
  ! ---------------------------------------------------------------------------
  ! IMF and coupling function
  real :: IO_Bt=1., IO_Ca=0., IO_Cf=3000.

  ! Reference Cfs 
  real, allocatable, dimension(:) :: epm_ref_cfs
  ! EPM grids 
  real, allocatable, dimension(:) :: epm_mlts, epm_mlats 
  ! EPM coefficients
  real, allocatable, dimension(:,:,:) :: epm_all_ca_coeffs
  real, allocatable, dimension(:) :: epm_cat0_coeffs
  real, allocatable, dimension(:) :: epm_slope_coeffs, epm_yint_coeffs

  ! Derived parameters 
  real :: epm_coeffs((epm_LMAX+1)**2), epm_coeffs_sh((epm_LMAX+1)**2)
  real :: epm_er, & ! (Expansion rate, NH=SH)
       epm_scale, epm_scale_sh, & ! (Scale factor for NH and SH)
       epm_disp, epm_disp_sh ! (Displacement, NH=-SH)
  
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  contains

    !!! OBTAIN IMF BY, BZ, SWVX, SWN
    ! -------------------------------------------------------------------------
    subroutine epm_readin(by,bz,swvx,swn)                                     
                                                                               
      implicit none                                                            
                                                                               
      real, intent(in) :: by, bz, swvx, swn                                    
                                                                               
      IO_IMFBy = by                                                            
      IO_IMFBz = bz                                                            
      IO_SWVX = abs(swvx)                                                      
      IO_SWN = swn                                                             
                                                                               
    end subroutine epm_readin


end module epm_interface
