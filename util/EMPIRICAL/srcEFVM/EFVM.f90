! Main part of the EFVM module
! Created: Qingyu Zhu, 08/19/2020
!
! -----------------------------------------------------------------------------
module efvm

  use efvm_initialization
  
  implicit none 

  contains

    !! INITIALIZATION FUNCTION
    ! -------------------------------------------------------------------------
    subroutine initialize_efvm

      call init_efvm_arrays                                                    
      if (efvm_debug) write(*,*) "Done init_efvm_arrays!"                      
                                                                               
      call set_efvm_ref_cfs                                                    
      if (efvm_debug) write(*,*) "Done set_efvm_ref_cfs!"                      
                                                                               
      call set_efvm_grids                                                      
      if (efvm_debug) write(*,*) "Done set_efvm_grids!"  

      call read_efvm_ca_coeffs                                                 
      if (efvm_debug) write(*,*) "Done read_efvm_ca_coeffs!"                   
                                                                               
      call read_efvm_cat0_coeffs                                               
      if (efvm_debug) write(*,*) "Done read_efvm_cat0_coeffs!"

    end subroutine initialize_efvm

    !! MAIN FUNCTION OF EFVM
    ! -------------------------------------------------------------------------
    subroutine efvm_main(er,disp,scale,scale_sh)

      ! -----------------------------------------------------------------------
      ! History Log:
      ! 08/20/2020:
      !    Will set disp as zero, and only expansion will be taken into account
      !    For large cf_in case, will use the methodology used in the AEPM
      ! -----------------------------------------------------------------------


      implicit none 

      real, intent(in) :: er, disp, scale, scale_sh
      real :: cf_in, ca_in, ca_in1

      real, parameter :: pi = 3.14159265359

      ! Calculate CF and CA
      call calc_cf_ca(IO_IMFBy, IO_IMFBz, IO_SWVX, IO_SWN, cf_in, ca_in)
      if (efvm_debug) write(*,*) "CF, CA:", cf_in, ca_in
      
      ! NH 
      call calc_efvm_coeffs(cf_in,efvm_cf_inf,efvm_ncf,efvm_ref_cfs,efvm_cf0,&
           ca_in,efvm_nmlat,efvm_LMAX,efvm_NMAX,&
           efvm_coeffs1,efvm_coeffs2)

      !if (cf_in<=efvm_ref_cfs(efvm_ncf)) then
      call reconst_efvm_dEd(efvm_coeffs1,efvm_nmlt,efvm_mlts,efvm_nmlat,&
           efvm_LMAX,efvm_dEd1)
      call reconst_efvm_dEd(efvm_coeffs2,efvm_nmlt,efvm_mlts,efvm_nmlat,&
           efvm_LMAX,efvm_dEd2)
      !else
      !   call calc_extrapolated_dEd(efvm_coeffs1,efvm_LMAX,efvm_nmlat,&
      !        efvm_nmlt,efvm_mlats,efvm_mlts,scale,er,disp,efvm_dEd1)
      !   call calc_extrapolated_dEd(efvm_coeffs2,efvm_LMAX,efvm_nmlat,&
      !        efvm_nmlt,efvm_mlats,efvm_mlts,scale,er,disp,efvm_dEd2)
      !end if

      efvm_dEd1 = efvm_dEd1 * scale
      efvm_dEd2 = efvm_dEd2 * scale

      ! SH
      ca_in1=2*pi-ca_in

      call calc_efvm_coeffs(cf_in,efvm_cf_inf,efvm_ncf,efvm_ref_cfs,efvm_cf0,&
           ca_in,efvm_nmlat,efvm_LMAX,efvm_NMAX,&
           efvm_coeffs1_sh,efvm_coeffs2_sh)

      !if (cf_in<=efvm_ref_cfs(efvm_ncf)) then
      call reconst_efvm_dEd(efvm_coeffs1_sh,efvm_nmlt,efvm_mlts,efvm_nmlat,&
           efvm_LMAX,efvm_dEd1_sh)
      call reconst_efvm_dEd(efvm_coeffs2_sh,efvm_nmlt,efvm_mlts,efvm_nmlat,&
           efvm_LMAX,efvm_dEd2_sh)
      !else
      !   call calc_extrapolated_dEd(efvm_coeffs1_sh,efvm_LMAX,efvm_nmlat,&
      !        efvm_nmlt,efvm_mlats,efvm_mlts,er,-disp,efvm_dEd1_sh)
      !   call calc_extrapolated_dEd(efvm_coeffs2_sh,efvm_LMAX,efvm_nmlat,&
      !        efvm_nmlt,efvm_mlats,efvm_mlts,er,-disp,efvm_dEd2_sh)
      !end if
      
      efvm_dEd1_sh = efvm_dEd1_sh * scale_sh
      efvm_dEd2_sh = efvm_dEd2_sh * scale_sh

      call calc_expanded_grid(efvm_nmlat,efvm_mlats,&                
           er,efvm_mlats1)

      if (efvm_debug) write(*,*) er, efvm_mlats1

    end subroutine efvm_main

end module efvm
