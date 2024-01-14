! This module serves as the initialization module for the EPM module 
! Updated from ModSEP_Initialization.f90 
! Created: Qingyu Zhu, 05/25/2020
!
! -----------------------------------------------------------------------------

module epm_initialization 

  use epm_interface

  implicit none 

  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  contains 

    !!! Initialize all necessary EPM arrays 
    ! ------------------------------------------------------------------------
    subroutine init_epm_arrays

      implicit none 

      integer :: istat

      !! Allocate Coefficients 
      if (.not. allocated(epm_all_ca_coeffs)) then

         allocate(epm_all_ca_coeffs(epm_ncf,(epm_LMAX+1)**2,(2*epm_NMAX+1)),&
              stat=istat)

         !if (epm_debug .and. (istat==0)) &
         !     write(*,*) "Fourier Coefficients are allocated successfully !" 
      end if

      if (.not. allocated(epm_cat0_coeffs)) then

         allocate(epm_cat0_coeffs((epm_LMAX+1)**2), stat=istat)

         !if (epm_debug .and. (istat==0)) & 
         !     write(*,*) "Cat-0 Coefficients are allocated successfully !"

      end if

      if (.not. allocated(epm_mlts)) then
         allocate(epm_mlts(epm_nmlt+1),stat=istat)
         allocate(epm_mlats(epm_nmlat+1),stat=istat)

         !if (epm_debug .and. (istat==0)) &
         !     write(*,*) "EPM grids are allocated successfully"
      end if

      if (.not. allocated(epm_slope_coeffs)) then
         allocate(epm_slope_coeffs(2*epm_NMAX1+1),stat=istat)   
         allocate(epm_yint_coeffs(2*epm_NMAX1+1),stat=istat)
         
         !if (epm_debug .and. (istat==0)) &                                   
         !     write(*,*) "Slope/yint Coefficients are allocated successfully !"

      end if

    end subroutine init_epm_arrays

    !!! Set up Reference CFs 
    ! -------------------------------------------------------------------------
    subroutine set_epm_ref_cfs 

      implicit none 

      integer :: istat 

      ! Allocate reference bts                                                 
      if (.not. allocated(epm_ref_cfs)) then                                  
         allocate(epm_ref_cfs(epm_ncf),stat=istat)                           
         !if (epm_debug .and. (istat>=0)) &            
         !     write(*,*) "REFERENCE CFs are allocated successfully !"        
      end if

      epm_ref_cfs =(/4939.,6997.,9153.,11955.,17654./)

    end subroutine set_epm_ref_cfs
    
    !!! Set up the EPM grids 
    ! -------------------------------------------------------------------------
    subroutine set_epm_grids 

      implicit none 
                                                                               
      integer :: imlt, imlat                                                 
      real :: step1, step2

      step1=24./epm_nmlt
      step2=1.

      epm_mlts=0.
      epm_mlats=0.

      do imlt=1,epm_nmlt+1                                                     
         epm_mlts(imlt)=(imlt-1)*step1                               
      end do                                                                   
                                                                               
      do imlat=1,epm_nmlat+1                                                   
         epm_mlats(imlat)=90-epm_nmlat+step2*(imlat-1)               
      end do  

    end subroutine set_epm_grids

    !!! Read Ca coefficents 
    ! -------------------------------------------------------------------------
    subroutine read_epm_ca_coeffs 

      implicit none 

      integer :: icf,ii
      character(len=2) :: str1
      
      do icf=1, epm_ncf 
         
         write(str1,'(I2.2)') icf

         open(unit=10,file=epm_coeff_path//epm_ca_fourier_coeff_path&
              //'bin_bt_'//str1//'.txt',&
              status='old')

         do ii=1,(epm_LMAX+1)**2
            read(10,*) epm_all_ca_coeffs(icf,ii,:)
         end do

         close(10)

      end do

    end subroutine read_epm_ca_coeffs
    
    !!! Read Cat-0 Coefficients 
    ! -------------------------------------------------------------------------
    subroutine read_epm_cat0_coeffs 

      implicit none 

      integer :: ii

      open(unit=10,file=epm_coeff_path//epm_cat0_coeff_path&
           //'bin_cat0.txt',status='old')           
                                                                              
      do ii=1,(epm_LMAX+1)**2                                           
         read(10,*) epm_cat0_coeffs(ii)                           
      end do                                                                   
                                                                               
      close(10) 

    end subroutine read_epm_cat0_coeffs

    !!! Read slope and yint coefficients 
    ! -------------------------------------------------------------------------
    subroutine read_epm_slope_yint_coeffs 

      implicit none 
      
      character(len=180) :: fname

      ! Read slope coeffs 
      fname=epm_slope_yint_path//'slope_fourier.txt' 
      open(unit=10,file=trim(fname),status='old') 
      read(10,*) epm_slope_coeffs(:)
      close(10)

      ! Read yint coeffs 
      fname=epm_slope_yint_path//'yints_fourier.txt' 
      open(unit=10,file=trim(fname),status='old') 
      read(10,*) epm_yint_coeffs(:)
      close(10)

    end subroutine read_epm_slope_yint_coeffs

end module epm_initialization
