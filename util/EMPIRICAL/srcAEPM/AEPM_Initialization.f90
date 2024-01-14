! This module serves as the initialization module for the AEPM module
! Updated from the ModSAM_Initialization.f90
! Created: Qingyu Zhu, 05/18/2020
!
! ----------------------------------------------------------------------------

module aepm_initialization

  use aepm_interface

  implicit none 

  ! --------------------------------------------------------------------------
  contains

    subroutine init_aepm_arrays

      implicit none

      integer :: istat

      ! Allocate Coeffcients
      if (.not. allocated(aepm_all_ca_coeffs)) then 
         allocate(aepm_all_ca_coeffs(aepm_ncf,aepm_nchannel,aepm_nmlat, &
              2*aepm_LMAX+1,2*aepm_NMAX+1),stat=istat)

         if (istat>=0 .and. aepm_debug) &
              write(*,*) "Fourier Coefficients are allocated successfully !"
      end if

      if (.not. allocated(aepm_cat0_coeffs)) then
         allocate(aepm_cat0_coeffs(aepm_nchannel,aepm_nmlat,2*aepm_LMAX+1), &
              stat=istat)
         if (istat>=0 .and. aepm_debug) &
              write(*,*) "Cat0 Coefficients are allocated successfully !"
      end if

      ! Allocate the reconstructed differential energy flux 
      if (.not. allocated(aepm_diff_ef)) then
         allocate(aepm_diff_ef(aepm_nmlt,aepm_nmlat,aepm_nchannel), &
              stat=istat)
         allocate(aepm_diff_ef_sh(aepm_nmlt,aepm_nmlat,aepm_nchannel), &
              stat=istat)

         if (istat>=0 .and. aepm_debug) &
              write(*,*) "Diff_ef arrays are allocated successfully!"
      end if

      ! Allocate the total energy flux and total number flux
      if (.not. allocated(aepm_total_ef)) then
         allocate(aepm_total_ef(aepm_nmlt,aepm_nmlat), stat=istat)
         allocate(aepm_total_nf(aepm_nmlt,aepm_nmlat), stat=istat)
         allocate(aepm_total_ef_sh(aepm_nmlt,aepm_nmlat), stat=istat)
         allocate(aepm_total_nf_sh(aepm_nmlt,aepm_nmlat), stat=istat)
         if (istat>=0 .and. aepm_debug) &
              write(*,*) "Total fluxes are allocated successfully !"
      end if

      ! Allocate the slope and yint and their Ca Fourier Coeffs 
      if (.not. allocated(aepm_slope_coeffs)) then
         allocate(aepm_slope_coeffs(aepm_nchannel,2*aepm_NMAX+1),stat=istat)
         allocate(aepm_yint_coeffs(aepm_nchannel,2*aepm_NMAX+1),stat=istat)
         allocate(aepm_slope(aepm_nchannel),stat=istat)
         allocate(aepm_yint(aepm_nchannel),stat=istat)

         if (istat>=0 .and. aepm_debug) &
              write(*,*) "Slope and Yint coeffs are allocated successfuly !"

      end if

    end subroutine init_aepm_arrays

    ! ------------------------------------------------------------------------
    subroutine set_aepm_ref_cfs

      implicit none

      integer :: istat

      ! Allocate reference bts 
      if (.not. allocated(aepm_ref_cfs)) then
         allocate(aepm_ref_cfs(aepm_ncf),stat=istat)
         if (istat>=0 .and. aepm_debug) &
              write(*,*) "REFERENCE CFs are allocated successfully !"
      end if


      aepm_ref_cfs = (/4283.2,6073.3,7957.8,9929.7,11941.6,&
           14254.0,17590.4,22770.4/)

    end subroutine set_aepm_ref_cfs

    ! ------------------------------------------------------------------------
    subroutine set_aepm_grids 

      implicit none 

      integer :: imlt, imlat, istat

      ! Allocate AEPM Grids 
      if (.not. allocated(aepm_mlats)) then
         allocate(aepm_mlats(aepm_nmlat),stat=istat)
         allocate(aepm_mlats1(aepm_nmlat),stat=istat)
         allocate(aepm_mlts(aepm_nmlt),stat=istat)
         if (istat>=0 .and. aepm_debug) &
              write(*,*) "AEPM grids are allocated successfully !"
      end if
      
      aepm_mlts = (/(imlt, imlt = 1,aepm_nmlt)/)-0.5
      aepm_mlats = (/(imlat, imlat = 1,aepm_nmlat)/)+49.5

    end subroutine set_aepm_grids

    ! ------------------------------------------------------------------------
    subroutine set_aepm_channels 
      
      implicit none

      integer :: istat, ichannel
      real :: log_min, log_max, log_diff

      ! Allocate the Channels 
      if (.not. allocated(aepm_channels)) then 
         allocate(aepm_channels(aepm_nchannel),stat=istat)
         if (istat>=0 .and. aepm_debug) &
              write (*,*) "Channels are allocated successfully !"
      end if

      ! Calculate the mid-points of each channels 
      log_max=log10(30000.)
      log_min=log10(30.)

      log_diff=(log_max-log_min)/aepm_nchannel

      do ichannel=1,aepm_nchannel

         aepm_channels(ichannel)=10**(log_max-log_diff*ichannel+0.5*log_diff)

      end do

    end subroutine set_aepm_channels

    ! ------------------------------------------------------------------------
    subroutine read_aepm_ca_coeffs

      implicit none

      integer :: icf, ichannel, imlat, ii, iline

      character(len=2) :: str1, str2, str3
      character(len=180) :: path

      do icf=1,aepm_ncf

         write(str1,'(I2.2)') icf

         do ichannel=1,aepm_nchannel

            write(str2,'(I2.2)') ichannel-1

            path=aepm_coeff_path//aepm_ca_fourier_coeff_path//'bin_cf_'//&
                 str1//'/ch_'//str2//'/'
            
            do imlat=1,aepm_nmlat

               write(str3,'(I2.2)') imlat-1

               open(unit=10,file=trim(path)//'mlat_'//str3//'.txt',&
                    status='old')

               do ii=1,2*aepm_LMAX+1
                  read(10,*) iline, aepm_all_ca_coeffs(icf,ichannel,imlat,ii,:)
               end do

               close(10)

            end do ! MLAT

         end do ! CHANNEL

      end do ! CF
      
    end subroutine read_aepm_ca_coeffs

    ! ------------------------------------------------------------------------
    subroutine read_aepm_cat0_coeffs

      implicit none

      integer :: ichannel, ii, iline
      character(len=2) :: str1
      character(len=180) :: path

      path=aepm_coeff_path//aepm_cat0_coeff_path

      if (aepm_debug) write(*,*) path

      do ichannel=1,aepm_nchannel

         write(str1,'(I2.2)') ichannel-1

         open(unit=10,file=trim(path)//'channel_'//str1//'.txt',&
              status='old')

         do ii=1,aepm_nmlat
            read(10,*) iline, aepm_cat0_coeffs(ichannel,ii,:)
         end do
         
         close(10)

      end do

    end subroutine read_aepm_cat0_coeffs

    !!! Read Slope and yint coeffs 
    ! -------------------------------------------------------------------------
    subroutine read_aepm_slope_yint_coeffs

      implicit none 
      
      integer :: ichannel, iline
      character(len=180) :: fname

      ! Read slope coeffs 
      fname=aepm_slope_yint_path//'slope_fourier.txt'
      open(unit=10,file=trim(fname),status='old')

      ! Currently stop at the Channel 11
      do ichannel=1,11 
         read(10,*) iline, aepm_slope_coeffs(ichannel,:)
      end do
      
      close(10)

      ! Read yint coeffs 
      fname=aepm_slope_yint_path//'yints_fourier.txt'
      open(unit=10,file=trim(fname),status='old')

      ! Currently stop at the Channel 11
      do ichannel=1,11 
         read(10,*) iline, aepm_yint_coeffs(ichannel,:)
      end do
      
      close(10)
      
    end subroutine read_aepm_slope_yint_coeffs

end module aepm_initialization
