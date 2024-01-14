! Routines using AMIE outputs from Gang
! Created: Qingyu Zhu, 02/01/2021
! -----------------------------------------------------------------------------
module ModAMIE

  use netcdf

  implicit none 

  character(len=*), parameter :: & 
       amie_pt='/glade/u/home/qingyu/data/amie/201303/', &
       amie_fn_nh='mar17-18_2013_all_nh.nc', &
       amie_fn_sh='mar17-18_2013_all_sh.nc' 

  integer :: amie_nmlt, amie_nmlat, amie_ntime

  real, allocatable :: amie_mlts(:),amie_mlats(:) 

  integer, allocatable :: amie_year(:), amie_mon(:), amie_day(:), &
       amie_hour(:), amie_minute(:) 

  real, allocatable :: amie_time(:) 

  real, allocatable :: amie_pot_nh(:,:,:), amie_pot_sh(:,:,:), &
       amie_efx_nh(:,:,:), amie_efx_sh(:,:,:), & 
       amie_ekv_nh(:,:,:), amie_ekv_sh(:,:,:), & 
       amie_cusplat_nh(:), amie_cusplat_sh(:), &
       amie_hpi_nh(:), amie_hpi_sh(:)

  real, allocatable :: amie_potin_nh(:,:), amie_potin_sh(:,:), &
       amie_efxin_nh(:,:), amie_efxin_sh(:,:), &
       amie_ekvin_nh(:,:), amie_ekvin_sh(:,:)

  real :: amie_cusplatin_nh, amie_cusplatin_sh

  real :: amie_hpiin_nh, amie_hpiin_sh 

  logical, parameter :: amie_debug = .false.

  contains
    ! -------------------------------------------------------------------------
    ! Read in NH/SH AMIE outputs
    ! Created: Qingyu Zhu, 02/01/2021
    
    subroutine read_ncar_amie

      implicit none 

      integer :: ncid, istat
      integer :: nDims, nVars, nGlobalAtts, unlimDimID
      integer :: varid, dimid 
      integer :: i, j, k  
      integer :: nfile, nmlat, nmlt, len

      character(len=NF90_MAX_NAME) :: varname, dimname 
      integer, dimension(nf90_max_var_dims) :: DimIds
      real, allocatable, dimension(:) :: ut 
      real :: minute    
      integer :: iTime(7)
      

      ! ---------------------------- NH ----------------------------
      istat=nf90_open(path=amie_pt//amie_fn_nh,mode=nf90_nowrite,ncid=ncid)  

      if (amie_debug .and. (istat==nf90_noerr)) & 
           write(*,*) " ------ NH AMIE file opened"

      ! Dimension
      istat = nf90_inquire(ncid, ndims, nvars)
      
      do i=1,ndims 
         istat=nf90_inquire_dimension(ncid, i, dimname, len)
         if (i==1) nfile=len    
         if (i==2) nmlat=len 
         if (i==4) nmlt=len 
      end do

      amie_ntime=nfile
      amie_nmlt=nmlt                                                  
      amie_nmlat=nmlat

      !if (amie_debug) write(*,*) "Dims", amie_ntime, amie_nmlt, amie_nmlat

      ! MLT & MLAT
      if (.not. allocated(amie_mlts)) allocate(amie_mlts(amie_nmlt),stat=istat)
      istat=nf90_get_var(ncid, 10, amie_mlts)    
      !if (amie_debug .and. (istat==0)) write(*,*) "Get MLT values"

      if (.not. allocated(amie_mlats)) &
           allocate(amie_mlats(amie_nmlat),stat=istat)
      istat=nf90_get_var(ncid, 8, amie_mlats)   
      !if (amie_debug .and. (istat==0)) write(*,*) "Get MLAT values"

      ! Time
      if (.not. allocated(amie_year)) then                                
         allocate(amie_year(amie_ntime),stat=istat)                     
         allocate(amie_mon(amie_ntime),stat=istat)                     
         allocate(amie_day(amie_ntime),stat=istat)                    
         allocate(amie_hour(amie_ntime),stat=istat)                 
         allocate(amie_minute(amie_ntime),stat=istat)                 
         allocate(ut(amie_ntime),stat=istat)        
         allocate(amie_time(amie_ntime),stat=istat)   
      end if

      istat=nf90_get_var(ncid, 16, amie_year)
      istat=nf90_get_var(ncid, 11, amie_mon)
      istat=nf90_get_var(ncid, 3, amie_day)
      istat=nf90_get_var(ncid, 15, ut)
      
      do i=1,amie_ntime

         amie_hour(i)=int(ut(i)) 
         amie_minute(i)=nint((ut(i)-int(ut(i)))*60)

         iTime=0
         iTime(1)=amie_year(i)                                         
         iTime(2)=amie_mon(i)                                            
         iTime(3)=amie_day(i)                                            
         iTime(4)=amie_hour(i)                                           
         iTime(5)=amie_minute(i) 

         call time_int_to_real(iTime,amie_time(i))

         if (amie_debug .and. (i==1)) write(*,*) iTime,amie_time(1) 

      end do

      ! Potential & electron precipitation
      if (.not. allocated(amie_pot_nh)) then
         allocate(amie_pot_nh(amie_nmlt,amie_nmlat,amie_ntime),stat=istat)
         allocate(amie_efx_nh(amie_nmlt,amie_nmlat,amie_ntime),stat=istat)
         allocate(amie_ekv_nh(amie_nmlt,amie_nmlat,amie_ntime),stat=istat)   
      end if

      istat=nf90_get_var(ncid, 13, amie_pot_nh) 
      istat=nf90_get_var(ncid, 4, amie_efx_nh)
      istat=nf90_get_var(ncid, 5, amie_ekv_nh)      

      ! CuspLat
      if (.not. allocated(amie_cusplat_nh)) &
           allocate(amie_cusplat_nh(amie_ntime),stat=istat)

      istat=nf90_get_var(ncid, 1, amie_cusplat_nh)

      ! HP
      if (.not. allocated(amie_hpi_nh)) &
           allocate(amie_hpi_nh(amie_ntime),stat=istat)

      istat=nf90_get_var(ncid, 6, amie_hpi_nh)


      istat=nf90_close(ncid) 
      if (amie_debug .and. (istat==nf90_noerr)) & 
           write(*,*) " ------ NH AMIE file closed"

      ! ---------------------------- SH ----------------------------
      istat=nf90_open(path=amie_pt//amie_fn_sh,mode=nf90_nowrite,ncid=ncid)  
      if (amie_debug .and. (istat==nf90_noerr)) & 
           write(*,*) " ------ SH AMIE file opened"

      ! Potential & electron precipitation
      if (.not. allocated(amie_pot_sh)) then
         allocate(amie_pot_sh(amie_nmlt,amie_nmlat,amie_ntime),stat=istat)
         allocate(amie_efx_sh(amie_nmlt,amie_nmlat,amie_ntime),stat=istat)
         allocate(amie_ekv_sh(amie_nmlt,amie_nmlat,amie_ntime),stat=istat)   
      end if

      istat=nf90_get_var(ncid, 13, amie_pot_sh) 
      istat=nf90_get_var(ncid, 4, amie_efx_sh)
      istat=nf90_get_var(ncid, 5, amie_ekv_sh)      

      ! CuspLat
      if (.not. allocated(amie_cusplat_sh)) &
           allocate(amie_cusplat_sh(amie_ntime),stat=istat)

      istat=nf90_get_var(ncid, 1, amie_cusplat_sh)

      ! HP
      if (.not. allocated(amie_hpi_sh)) &
           allocate(amie_hpi_sh(amie_ntime),stat=istat)

      istat=nf90_get_var(ncid, 6, amie_hpi_sh)


      istat=nf90_close(ncid) 
      if (amie_debug .and. (istat==nf90_noerr)) & 
           write(*,*) " ------ SH AMIE file closed"

      ! -------------------------- Deallocation --------------------------
      deallocate(amie_year)
      deallocate(amie_mon)
      deallocate(amie_day)
      deallocate(amie_hour)
      deallocate(amie_minute)
      deallocate(ut)

    end subroutine read_ncar_amie

    ! TEMPORAL INTERPOLATION
    ! Created: Qingyu Zhu, 02/01/2021
    ! -------------------------------------------------------------------------
    subroutine get_currenttime_amie
    
      use ModTime, only: CurrentTime 
      
      implicit none 

      real :: wgt1, wgt2  
      integer :: left, right, i, istat

      left=1
      right=1
      wgt1=0.
      wgt2=0.

      if (.not. allocated(amie_potin_nh)) then
         allocate(amie_potin_nh(amie_nmlt,amie_nmlat),stat=istat) 
         allocate(amie_potin_sh(amie_nmlt,amie_nmlat),stat=istat)
         allocate(amie_efxin_nh(amie_nmlt,amie_nmlat),stat=istat) 
         allocate(amie_efxin_sh(amie_nmlt,amie_nmlat),stat=istat)
         allocate(amie_ekvin_nh(amie_nmlt,amie_nmlat),stat=istat) 
         allocate(amie_ekvin_sh(amie_nmlt,amie_nmlat),stat=istat)
      end if

      if (CurrentTime<amie_time(1)) then
         left=1                                                                
         right=1                                                               
         wgt1=0.                                                               
         wgt2=1.                                                               
      else if (CurrentTime>=amie_time(amie_ntime)) then             
         left=amie_ntime                                          
         right=amie_ntime                              
         wgt1=1.                                                               
         wgt2=0.
      else                                                                  
         do i=1,amie_ntime-1                                                  
            if (CurrentTime>=amie_time(i) .and. &                        
                 CurrentTime<amie_time(i+1)) then               
               left=i                                                          
               right=i+1                                                       
               wgt1=(amie_time(i+1)-CurrentTime)/&                        
                    (amie_time(i+1)-amie_time(i))                 
               wgt2=1-wgt1                                                     
               exit                                                            
            end if                                                             
         end do                                                                
      end if

      if (amie_debug) then
         if ((wgt1+wgt2)<1.) then                                              
            write(*,*) "Time is not found, check !!!"                          
         else                                                                  
            write(*,*) "Found data for", CurrentTime, left, right              
         end if                                                                
      end if  

      ! Potential
      amie_potin_nh(:,:) = amie_pot_nh(:,:,left) * wgt1 + &
           amie_pot_nh(:,:,right) * wgt2 

      amie_potin_sh(:,:) = amie_pot_sh(:,:,left) * wgt1 + &
           amie_pot_sh(:,:,right) * wgt2 

      ! Eflux
      amie_efxin_nh(:,:) = amie_efx_nh(:,:,left) * wgt1 + &
           amie_efx_nh(:,:,right) * wgt2 

      amie_efxin_sh(:,:) = amie_efx_sh(:,:,left) * wgt1 + &
           amie_efx_sh(:,:,right) * wgt2       

      ! Average energy
      amie_ekvin_nh(:,:) = amie_ekv_nh(:,:,left) * wgt1 + &
           amie_ekv_nh(:,:,right) * wgt2 

      amie_ekvin_sh(:,:) = amie_ekv_sh(:,:,left) * wgt1 + &
           amie_ekv_sh(:,:,right) * wgt2 
      
      ! CuspLat
      amie_cusplatin_nh = amie_cusplat_nh(left) * wgt1 + &
           amie_cusplat_nh(right) * wgt2 

      amie_cusplatin_sh = amie_cusplat_sh(left) * wgt1 + &
           amie_cusplat_sh(right) * wgt2 

      ! HP
      amie_hpiin_nh = amie_hpi_nh(left) * wgt1 + &
           amie_hpi_nh(right) * wgt2 

      amie_hpiin_sh = amie_hpi_sh(left) * wgt1 + &
           amie_hpi_sh(right) * wgt2 

    end subroutine get_currenttime_amie

    ! Spatial interpolations
    ! Created: Qingyu Zhu, 02/01/2021
    ! -------------------------------------------------------------------------
    subroutine interp_amie(mlatin,mltin,mlats,mlts,val_in,val_out)

      implicit none 

      real, intent(in) :: mlatin, mltin                                   
      real, intent(in) :: mlats(amie_nmlat), mlts(amie_nmlt)               
      real, intent(in) :: val_in(amie_nmlt,amie_nmlat)                      
      real, intent(out) :: val_out

      integer :: i, j, pos1(2),pos2(2)
      real :: wgt1(2),wgt2(2)   

      pos1=1                                                                   
      pos2=1                                                                   
      wgt1=0.                                                                  
      wgt2=0.                                                                  
      val_out=0    

      ! MLAT
      if (mlatin<=mlats(amie_nmlat)) then                             
         wgt1=0.                                                               
         pos1=1                                                                
      else                                                                     
         do i=1,amie_nmlat-1                                            
            if ((mlatin<=mlats(i)) .and. (mlatin>mlats(i+1))) then             
               pos1(1)=i                                                       
               pos1(2)=i+1                                                     
               wgt1(1)=(mlatin-mlats(i+1))/(mlats(i)-mlats(i+1))               
               wgt1(2)=1.-wgt1(1)                                              
               exit                                                            
            end if                                                             
         end do

         if ((sum(wgt1)<0.5) .and. (mlatin<90.)) write(*,*) &    
              "MLAT is not found, check !!!", mlatin,mltin 
      end if

      ! MLT
      do j=1,amie_nmlt-1  
         if ((mltin>=mlts(j)) .and. (mltin<mlts(j+1))) then                    
            pos2(1)=j                                                          
            pos2(2)=j+1                                                        
            wgt2(1)=(mlts(j+1)-mltin)/(mlts(j+1)-mlts(j))                      
            wgt2(2)=1-wgt2(1)                                                  
            exit                                                               
         end if                                                                
      end do 
      
      if (sum(wgt2)<0.5) write(*,*) "MLT is not found, check !!!"

      ! Get the output
      val_out = val_in(pos2(1),pos1(1))*wgt1(1)*wgt2(1) + &
           val_in(pos2(1),pos1(2))*wgt1(2)*wgt2(1) + &  
           val_in(pos2(2),pos1(1))*wgt1(1)*wgt2(2) + & 
           val_in(pos2(2),pos1(2))*wgt1(2)*wgt2(2) 

    end subroutine interp_amie

end module ModAMIE
