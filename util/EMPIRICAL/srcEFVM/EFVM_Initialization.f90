! This module serves as the initialization module for the EFVM module 
! Created: Qingyu Zhu, 08/18/2020
! 
! -----------------------------------------------------------------------------
module efvm_initialization 

  use efvm_interface

  implicit none 

  contains
   
    !! INITIALIZE EFVM ARRAYS
    ! -------------------------------------------------------------------------
    subroutine init_efvm_arrays
    
      implicit none

      integer :: istat

      ! Allocate Ca coefficents
      if (.not. allocated(efvm_all_ca_coeffs1)) then

         allocate(efvm_all_ca_coeffs1(efvm_ncf,efvm_nmlat, &
              2*efvm_LMAX+1,2*efvm_NMAX+1),stat=istat)
         allocate(efvm_all_ca_coeffs2(efvm_ncf,efvm_nmlat, &
              2*efvm_LMAX+1,2*efvm_NMAX+1),stat=istat)
         
         if (istat>=0 .and. efvm_debug) &
              write(*,*) "Ca Fourier coefficients are allocated!"

      end if

      ! Allocate Cat-0 coefficents
      if (.not. allocated(efvm_cat0_coeffs1)) then

         allocate(efvm_cat0_coeffs1(efvm_nmlat,2*efvm_LMAX+1), &   
              stat=istat)
         allocate(efvm_cat0_coeffs2(efvm_nmlat,2*efvm_LMAX+1), &   
              stat=istat)

         if (istat>=0 .and. efvm_debug) &
              write(*,*) "Cat-0 MLT Fourier coefficients are allocated!"

      end if

      ! Allocate outputs 
      if (.not. allocated(efvm_dEd1)) then

         allocate(efvm_dEd1(efvm_nmlt,efvm_nmlat), stat=istat)
         allocate(efvm_dEd2(efvm_nmlt,efvm_nmlat), stat=istat)
         allocate(efvm_dEd1_sh(efvm_nmlt,efvm_nmlat), stat=istat)
         allocate(efvm_dEd2_sh(efvm_nmlt,efvm_nmlat), stat=istat)

         if (istat>=0 .and. efvm_debug) &
              write(*,*) "Output arrays are allocated!"

      end if

    end subroutine init_efvm_arrays

    !! SET UP THE REF CFS
    ! -------------------------------------------------------------------------
    subroutine set_efvm_ref_cfs

      implicit none 

      integer :: istat

      if (.not. allocated(efvm_ref_cfs)) then   
        
         allocate(efvm_ref_cfs(efvm_ncf),stat=istat)  
         if (istat>=0 .and. efvm_debug) & 
              write(*,*) "Reference CFs are allocated!"
              
      end if

      efvm_ref_cfs=(/4615.3,6539.2,8456.0,10614.1,13523.9,18357.4/)

    end subroutine set_efvm_ref_cfs

    !! SET UP THE EFVM GRIDS
    ! -------------------------------------------------------------------------
    subroutine set_efvm_grids

      implicit none 

      integer :: imlt, imlat, istat

      if (.not. allocated(efvm_mlats)) then

         allocate(efvm_mlats(efvm_nmlat),stat=istat)              
         allocate(efvm_mlats1(efvm_nmlat),stat=istat)                 
         allocate(efvm_mlts(efvm_nmlt),stat=istat)
         
         if (istat>=0 .and. efvm_debug) &
              write(*,*) "EFVM grids are allocated successfully !"

      end if

      efvm_mlts = (/(imlt, imlt = 1,efvm_nmlt)/)-0.5
      efvm_mlats = (/(imlat,imlat = 1, efvm_nmlat)/)*(40./efvm_nmlat)+49.

      if (efvm_debug) write(*,*) "MLT:", efvm_mlts
      if (efvm_debug) write(*,*) "MLAT:", efvm_mlats 

    end subroutine set_efvm_grids

    !! READ CA COEFFICIENTS
    ! -------------------------------------------------------------------------
    subroutine read_efvm_ca_coeffs

      implicit none

      integer :: icf, ichannel, imlat, ii, iline

      character(len=2) :: str1, str2, str3
      character(len=180) :: path1, path2

      do icf=1,efvm_ncf

         write(str1,'(I2.2)') icf 

         path1=efvm_coeff_path//efvm_ca_fourier_coeff_path1//'bin_cf_'//&
              str1//'/'

         path2=efvm_coeff_path//efvm_ca_fourier_coeff_path2//'bin_cf_'//&
              str1//'/'

         do imlat=1,efvm_nmlat

            write(str2,'(I2.2)') imlat-1
            
            ! Ed1
            open(unit=10,file=trim(path1)//'mlat_'//str2//'.txt',&            
                 status='old')

            do ii=1,2*efvm_LMAX+1
               read(10,*) iline, efvm_all_ca_coeffs1(icf,imlat,ii,:)
            end do

            close(10)

            ! Ed2
            open(unit=10,file=trim(path2)//'mlat_'//str2//'.txt',&            
                 status='old')

            do ii=1,2*efvm_LMAX+1
               read(10,*) iline, efvm_all_ca_coeffs2(icf,imlat,ii,:)
            end do

            close(10)

         end do ! MLAT

      end do ! CF

    end subroutine read_efvm_ca_coeffs

    !! READ CAT0 COEFFICIENTS
    ! -------------------------------------------------------------------------
    subroutine read_efvm_cat0_coeffs

      implicit none 

      integer :: ii, iline                     
      character(len=2) :: str1                                       
      character(len=180) :: path

      path=efvm_coeff_path//efvm_cat0_coeff_path

      ! Ed1
      open(unit=10,file=trim(path)//'Ed1.txt',status='old')
      
      do ii=1,efvm_nmlat
         read(10,*) iline, efvm_cat0_coeffs1(ii,:)
      end do

      close(10)

      ! Ed2
      open(unit=10,file=trim(path)//'Ed2.txt',status='old')
      
      do ii=1,efvm_nmlat
         read(10,*) iline, efvm_cat0_coeffs2(ii,:)
      end do

      close(10)
      
    end subroutine read_efvm_cat0_coeffs
      
end module efvm_initialization
