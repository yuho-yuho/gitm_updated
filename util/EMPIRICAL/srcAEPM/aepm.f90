! Main part of the AEPM module 
! Updated from ModSAM_main.f90 
! Created: Qingyu Zhu, 05/21/2020
!
! -----------------------------------------------------------------------------
module aepm

  use aepm_initialization 

  implicit none 

  ! ---------------------------------------------------------------------------
  contains 

!!! Initialization module 
    ! -------------------------------------------------------------------------
    subroutine initialize_aepm
      
      implicit none 
      
      call init_aepm_arrays
      if (aepm_debug) write(*,*) "Done init_aepm_arrays!" 
      
      call set_aepm_ref_cfs 
      if (aepm_debug) write(*,*) "Done set_aepm_ref_cfs!" 

      call set_aepm_grids 
      if (aepm_debug) write(*,*) "Done set_aepm_grids!" 

      call set_aepm_channels 
      if (aepm_debug) write(*,*) "Done set_aepm_channels!"

      call read_aepm_ca_coeffs 
      if (aepm_debug) write(*,*) "Done read_aepm_ca_coeffs!"

      call read_aepm_cat0_coeffs 
      if (aepm_debug) write(*,*) "Done read_aepm_cat0_coeffs!"

      call read_aepm_slope_yint_coeffs 
      if (aepm_debug) write(*,*) "Done read_aepm_slope_yint_coeffs!"

    end subroutine initialize_aepm

    
!!! Main calculation part
    ! -------------------------------------------------------------------------
    subroutine aepm_main
      
      implicit none 
      
      integer :: imlat
      real :: cf_in, ca_in
      real :: mlats(aepm_nmlat+1), NH_HP, SH_HP
      
      !! Get cf_in, ca_in
      call calc_cf_ca(IO_IMFBy, IO_IMFBz, IO_SWVX, IO_SWN, cf_in, ca_in)

      if (aepm_debug) write(*,*) cf_in, ca_in
      
      ! Grids for the HP calculation 
      mlats = (/(imlat, imlat = 1,aepm_nmlat+1)/)+49.
      
      !! NH Pattern 
      call calc_aepm_diff_ef(cf_in,ca_in,aepm_diff_ef)
      ! Calculate the total energy and number fluxes 
      call calc_total_ef_nf(aepm_nmlt,aepm_nmlat,aepm_nchannel,&
           aepm_channels,aepm_diff_ef,aepm_total_ef,aepm_total_nf)

      if (aepm_debug) then 
         write (*,*) aepm_diff_ef(:,15,4)
         write (*,*) "Total ef:",aepm_total_ef(:,15)
         write (*,*) aepm_expansion_rate
      end if

      ! Calculate HP
      call calc_hemi_int_val(aepm_nmlt,aepm_nmlat,mlats,&
           aepm_total_ef,1.0e-6,aepm_expansion_rate,NH_HP)

      !write(*,*) "==========================================================="
      if (aepm_debug) write(*,*) "NH_HP (GW):", NH_HP

      !! SH Pattern 
      call calc_aepm_diff_ef(cf_in,2*pi-ca_in,aepm_diff_ef_sh)
      ! Calculate the total energy and number fluxes 
      call calc_total_ef_nf(aepm_nmlt,aepm_nmlat,aepm_nchannel,&
           aepm_channels,aepm_diff_ef_sh,aepm_total_ef_sh,aepm_total_nf_sh)
      ! Calculate HP
      call calc_hemi_int_val(aepm_nmlt,aepm_nmlat,mlats,&
           aepm_total_ef_sh,1.0e-6,aepm_expansion_rate,SH_HP)
      
      if (aepm_debug) write(*,*) "SH_HP (GW):", SH_HP 

    end subroutine aepm_main

    ! Calculate the HP in NH and SH
    ! -------------------------------------------------------------------------
    subroutine aepm_calc_HP(NH_HP, SH_HP)
      
      implicit none 
      
      real, intent(out) :: NH_HP, SH_HP 

      integer :: imlat
      real :: cf_in, ca_in
      real :: mlats(aepm_nmlat+1)
      
      !! Get cf_in, ca_in
      call calc_cf_ca(IO_IMFBy, IO_IMFBz, IO_SWVX, IO_SWN, cf_in, ca_in)

      if (aepm_debug) write(*,*) cf_in, ca_in
      
      ! Grids for the HP calculation 
      mlats = (/(imlat, imlat = 1,aepm_nmlat+1)/)+49.
      
      !! NH Pattern 
      call calc_aepm_diff_ef(cf_in,ca_in,aepm_diff_ef)
      ! Calculate the total energy and number fluxes 
      call calc_total_ef_nf(aepm_nmlt,aepm_nmlat,aepm_nchannel,&
           aepm_channels,aepm_diff_ef,aepm_total_ef,aepm_total_nf)

      if (aepm_debug) then 
         write (*,*) aepm_diff_ef(:,15,4)
         write (*,*) "Total ef:",aepm_total_ef(:,15)
         write (*,*) aepm_expansion_rate
      end if

      ! Calculate HP
      call calc_hemi_int_val(aepm_nmlt,aepm_nmlat,mlats,&
           aepm_total_ef,1.0e-6,aepm_expansion_rate,NH_HP)

      !write(*,*) "==========================================================="
      if (aepm_debug) write(*,*) "NH_HP (GW):", NH_HP

      !! SH Pattern 
      call calc_aepm_diff_ef(cf_in,2*pi-ca_in,aepm_diff_ef_sh)
      ! Calculate the total energy and number fluxes 
      call calc_total_ef_nf(aepm_nmlt,aepm_nmlat,aepm_nchannel,&
           aepm_channels,aepm_diff_ef_sh,aepm_total_ef_sh,aepm_total_nf_sh)
      ! Calculate HP
      call calc_hemi_int_val(aepm_nmlt,aepm_nmlat,mlats,&
           aepm_total_ef_sh,1.0e-6,aepm_expansion_rate,SH_HP)
      
      if (aepm_debug) write(*,*) "SH_HP (GW):", SH_HP 

    end subroutine aepm_calc_HP


end module aepm

