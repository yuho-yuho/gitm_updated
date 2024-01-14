! The main module of the EPM 
! Updated from the ModSEP_main.f90
! Created: Qingyu Zhu 05/28/2020
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module epm

  use epm_interface
  use epm_initialization 

  implicit none 

  ! ---------------------------------------------------------------------------
  contains 

    !!! Initialization 
    ! -------------------------------------------------------------------------
    subroutine initialize_epm

      implicit none 
      
      call init_epm_arrays
      call set_epm_ref_cfs
      call set_epm_grids
      call read_epm_ca_coeffs
      call read_epm_cat0_coeffs
      call read_epm_slope_yint_coeffs

    end subroutine initialize_epm

    !!! Caculation of the potential, er, scale, disp
    ! -------------------------------------------------------------------------
    subroutine epm_main

      implicit none 

      real :: cf_in, ca_in, ca_in1
      
      !! Get cf_in, ca_in                                                      
      call calc_cf_ca(IO_IMFBy, IO_IMFBz, IO_SWVX, IO_SWN, cf_in, ca_in)
      if (epm_debug) write(*,*) cf_in, ca_in

      if (epm_debug) write(*,*) "======== NH ========"
      !! Get NH coeff, scale, er, disp
      call calc_epm_epot(cf_in,epm_cf_inf,epm_cf_inf1,ca_in,epm_LMAX,epm_NMAX,&
           epm_nmlt+1,epm_mlts,epm_nmlat+1,epm_mlats,&
           epm_coeffs,epm_scale,epm_er,epm_disp)

      if (epm_debug) write(*,*) "======== NH ========"
      if (epm_debug) write(*,*) " "
      
      if (epm_debug) write(*,*) "======== SH ========"  
      ca_in1=2*pi-ca_in
      !! Get SH coeff, scale, er, disp
      call calc_epm_epot(cf_in,epm_cf_inf,epm_cf_inf1,ca_in1,epm_LMAX,&
           epm_NMAX,epm_nmlt+1,epm_mlts,epm_nmlat+1,epm_mlats,&
           epm_coeffs_sh,epm_scale_sh,epm_er,epm_disp_sh)
      
      if (epm_debug) write(*,*) "======== SH ========"  

    end subroutine epm_main

  end module epm
