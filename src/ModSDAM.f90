! Module using self-defined auroral inputs 
! Created: Qingyu Zhu, 02/25/2021

module ModSDAM

  use netcdf

  implicit none                                                                
                                                                               
  character(len=*), parameter :: &                                             
       sdam_pt='UA/DataIn/', &                     
       sdam_fn_nh='aurora_NH6.nc', &          
       sdam_fn_sh='aurora_SH6.nc'                                    
                                                                               
  integer :: sdam_nmlt, sdam_nmlat, sdam_ntime
                               
  real, allocatable :: sdam_mlts(:),sdam_mlats(:)                              
                                                                               
  integer, allocatable :: sdam_year(:), sdam_mon(:), sdam_day(:), &            
       sdam_hour(:), sdam_minute(:)

  real, allocatable :: sdam_time(:)

  real, allocatable :: sdam_efx_nh(:,:,:), sdam_efx_sh(:,:,:), &
       sdam_nfx_nh(:,:,:), sdam_nfx_sh(:,:,:)  

  real, allocatable :: sdam_efxin_nh(:,:), sdam_efxin_sh(:,:), & 
       sdam_nfxin_nh(:,:), sdam_nfxin_sh(:,:)

  logical, parameter :: sdam_debug = .false.

  contains

    subroutine read_sdam

      implicit none                                                            
                                                                               
      integer :: ncid, istat                                                   
      integer :: nDims, nVars, nGlobalAtts, unlimDimID                         
      integer :: varid, dimid                                                  
      integer :: i, j, k                                                       
      integer :: nfile, nmlat, nmlt, len                                       
                                                                               
      character(len=NF90_MAX_NAME) :: varname, dimname                         
      integer, dimension(nf90_max_var_dims) :: DimIds 

      integer :: iTime(7)
      
      !!!!!! NH
      istat=nf90_open(path=sdam_pt//sdam_fn_nh,mode=nf90_nowrite,ncid=ncid)
      if (sdam_debug .and. (istat==nf90_noerr)) &
           write(*,*) " ------ NH SDAM file opened"

      ! Dimension                                                              
      istat = nf90_inquire(ncid, ndims, nvars)
      
      do i=1,ndims                                                             
         istat=nf90_inquire_dimension(ncid, i, dimname, len)                   
         if (i==3) nfile=len                                                   
         if (i==1) nmlat=len                                                   
         if (i==2) nmlt=len                                                    
      end do                                                                   
                                                                               
      sdam_ntime=nfile                                                         
      sdam_nmlt=nmlt                                                           
      sdam_nmlat=nmlat

      if (.not. allocated(sdam_mlts)) allocate(sdam_mlts(sdam_nmlt),stat=istat)
      istat=nf90_get_var(ncid, 2, sdam_mlts) 

      if (.not. allocated(sdam_mlats)) &
           allocate(sdam_mlats(sdam_nmlat),stat=istat)
      istat=nf90_get_var(ncid, 1, sdam_mlats)

      if (.not. allocated(sdam_year)) then                               
         allocate(sdam_year(sdam_ntime),stat=istat)                
         allocate(sdam_mon(sdam_ntime),stat=istat)                  
         allocate(sdam_day(sdam_ntime),stat=istat)                   
         allocate(sdam_hour(sdam_ntime),stat=istat)                  
         allocate(sdam_minute(sdam_ntime),stat=istat)
         allocate(sdam_time(sdam_ntime),stat=istat)
      end if

      istat=nf90_get_var(ncid, 4, sdam_year)                           
      istat=nf90_get_var(ncid, 5, sdam_mon)                    
      istat=nf90_get_var(ncid, 6, sdam_day)                    
      istat=nf90_get_var(ncid, 7, sdam_hour)                   
      istat=nf90_get_var(ncid, 8, sdam_minute)

      do i=1,sdam_ntime   

         iTime=0                                                               
         iTime(1)=sdam_year(i)                                                 
         iTime(2)=sdam_mon(i)                                                  
         iTime(3)=sdam_day(i)                                                  
         iTime(4)=sdam_hour(i)                                                 
         iTime(5)=sdam_minute(i)

         call time_int_to_real(iTime,sdam_time(i))

      end do

      if (.not. allocated(sdam_efx_nh)) then                                  
         allocate(sdam_efx_nh(sdam_nmlt,sdam_nmlat,sdam_ntime),stat=istat)
         allocate(sdam_nfx_nh(sdam_nmlt,sdam_nmlat,sdam_ntime),stat=istat)    
      end if

      istat=nf90_get_var(ncid, 9, sdam_efx_nh)
!      if (sdam_debug .and. (istat==nf90_noerr)) &
!           write(*,*) "NH max eflux for time 2", maxval(sdam_efx_nh(:,:,2))
      
      istat=nf90_get_var(ncid, 10, sdam_nfx_nh)
!      if (sdam_debug .and. (istat==nf90_noerr)) &                             
!           write(*,*) "NH max nflux for time 2", maxval(sdam_nfx_nh(:,:,2))

      istat=nf90_close(ncid)
      if (sdam_debug .and. (istat==nf90_noerr)) &                     
           write(*,*) " ------ NH SDAM file closed"

      !!!!!! SH
      istat=nf90_open(path=sdam_pt//sdam_fn_sh,mode=nf90_nowrite,ncid=ncid)
      if (sdam_debug .and. (istat==nf90_noerr)) &
           write(*,*) " ------ SH SDAM file opened"

      if (.not. allocated(sdam_efx_sh)) then                                  
         allocate(sdam_efx_sh(sdam_nmlt,sdam_nmlat,sdam_ntime),stat=istat)
         allocate(sdam_nfx_sh(sdam_nmlt,sdam_nmlat,sdam_ntime),stat=istat)    
      end if

      istat=nf90_get_var(ncid, 9, sdam_efx_sh)
!      if (sdam_debug .and. (istat==nf90_noerr)) &
!           write(*,*) "NH max eflux for time 2", maxval(sdam_efx_nh(:,:,2))

      istat=nf90_get_var(ncid, 10, sdam_nfx_sh)
!      if (sdam_debug .and. (istat==nf90_noerr)) &
!           write(*,*) "NH max nflux for time 2", maxval(sdam_nfx_nh(:,:,2))


      istat=nf90_close(ncid)
      if (sdam_debug .and. (istat==nf90_noerr)) &                     
           write(*,*) " ------ SH SDAM file closed"

      deallocate(sdam_year)                                                    
      deallocate(sdam_mon)                                                     
      deallocate(sdam_day)                                                     
      deallocate(sdam_hour)                                                    
      deallocate(sdam_minute)                                                

    end subroutine read_sdam

    ! -------------------------------------------------------------------------

    subroutine get_currenttime_sdam
      
      use ModTime, only: CurrentTime
      
      implicit none                                                            
                                                                               
      real :: wgt1, wgt2                                                      
      integer :: left, right, i, istat                                         
      
      left=1                                                                  
      right=1                                                                 
      wgt1=0.                                                                  
      wgt2=0.

      if (.not. allocated(sdam_efxin_nh)) then                               
         allocate(sdam_efxin_nh(sdam_nmlt,sdam_nmlat),stat=istat)        
         allocate(sdam_efxin_sh(sdam_nmlt,sdam_nmlat),stat=istat)     
         allocate(sdam_nfxin_nh(sdam_nmlt,sdam_nmlat),stat=istat)          
         allocate(sdam_nfxin_sh(sdam_nmlt,sdam_nmlat),stat=istat)          
      end if

      if (CurrentTime<sdam_time(1)) then
         left=1                                                                
         right=1                                                               
         wgt1=0.                                                               
         wgt2=1.                                                               
      else if (CurrentTime>=sdam_time(sdam_ntime)) then             
         left=sdam_ntime                                          
         right=sdam_ntime                              
         wgt1=1.                                                               
         wgt2=0.
      else                                                                  
         do i=1,sdam_ntime-1                                                  
            if (CurrentTime>=sdam_time(i) .and. &                        
                 CurrentTime<sdam_time(i+1)) then               
               left=i                                                          
               right=i+1                                                       
               wgt1=(sdam_time(i+1)-CurrentTime)/&                        
                    (sdam_time(i+1)-sdam_time(i))                 
               wgt2=1-wgt1                                                     
               exit                                                            
            end if
         end do
      end if
          
      if (sdam_debug) then
         if ((wgt1+wgt2)<1.) then                                              
            write(*,*) "Time is not found, check !!!"                          
         else                                                                  
            write(*,*) "Found data for", CurrentTime, left, right, wgt1, wgt2
         end if
      end if    
      
      ! Eflux
      sdam_efxin_nh(:,:) = sdam_efx_nh(:,:,left) * wgt1 + &
           sdam_efx_nh(:,:,right) * wgt2 
      
      sdam_efxin_sh(:,:) = sdam_efx_sh(:,:,left) * wgt1 + &
           sdam_efx_sh(:,:,right) * wgt2        
     
      ! Nflux
      sdam_nfxin_nh(:,:) = sdam_nfx_nh(:,:,left) * wgt1 + &
           sdam_nfx_nh(:,:,right) * wgt2 
      
      sdam_nfxin_sh(:,:) = sdam_nfx_sh(:,:,left) * wgt1 + &
           sdam_nfx_sh(:,:,right) * wgt2 

      if (sdam_debug) then
         write(*,*) "NH max eflux", maxval(sdam_efxin_nh)
         write(*,*) "SH max eflux", maxval(sdam_efxin_sh) 
         write(*,*) "NH max nflux", maxval(sdam_nfxin_nh)
         write(*,*) "SH max nflux", maxval(sdam_nfxin_sh) 
      end if

    end subroutine get_currenttime_sdam
    
    ! Spatial interpolation
    ! -----------------------------------------------------------------------
    subroutine interp_sdam(mlatin,mltin,mlats,mlts,val_in,val_out)
      
      implicit none
      
      real, intent(in) :: mlatin, mltin                                   
      real, intent(in) :: mlats(sdam_nmlat), mlts(sdam_nmlt)               
      real, intent(in) :: val_in(sdam_nmlt,sdam_nmlat)                      
      real, intent(out) :: val_out
      
      integer :: i, j, pos1(2),pos2(2)                                      
      real :: wgt1(2),wgt2(2)
      
      pos1=1                                                                   
      pos2=1                                                                   
      wgt1=0.                                                                  
      wgt2=0.                                                                  
      val_out=0

      !!! MLAT
      ! Different from SDAM where the MLAT is decreasing
      ! SDAM uses increasing MLAT grids 
      if (mlatin<mlats(1)) then
         wgt1=0.                                                               
         pos1=1
      else if (mlatin>=mlats(sdam_nmlat)) then
         wgt1=0.
         pos1=sdam_nmlat
      else
         do i=1,sdam_nmlat-1
            if ((mlatin>=mlats(i)) .and. (mlatin<mlats(i+1))) then
               pos1(1)=i                                                   
               pos1(2)=i+1
               wgt1(1)=(mlats(i+1)-mlatin)/(mlats(i+1)-mlats(i))
               wgt1(2)=1.-wgt1(1)
               exit
            end if
         end do
         
         if ((sum(wgt1)<0.5) .and. (mlatin<90.)) write(*,*) &
              "MLAT is not found, check !!!", mlatin,mltin
      end if
  
      ! MLT (0.5*(2n+1) grid)
      if (mltin<mlts(1) .or. mltin>=mlts(sdam_nmlt)) then

         pos2(1)=sdam_nmlt
         pos2(2)=1
         wgt2(1)=(mlts(1)-mltin)/(mlts(1)+24-mlts(sdam_nmlt))
         if (wgt2(1)<0) &
              wgt2(1)=(mlts(1)+24-mltin)/(mlts(1)+24-mlts(sdam_nmlt))
         wgt2(2)=1-wgt2(1)

      else
         
         do j=1,sdam_nmlt-1                                                  
            if ((mltin>=mlts(j)) .and. (mltin<mlts(j+1))) then
               pos2(1)=j                                 
               pos2(2)=j+1                                       
               wgt2(1)=(mlts(j+1)-mltin)/(mlts(j+1)-mlts(j))           
               wgt2(2)=1-wgt2(1)
               exit
            end if
         end do

      end if
            
      ! Get the output
      val_out = val_in(pos2(1),pos1(1))*wgt1(1)*wgt2(1) + &                    
           val_in(pos2(1),pos1(2))*wgt1(2)*wgt2(1) + &                         
           val_in(pos2(2),pos1(1))*wgt1(1)*wgt2(2) + &                         
           val_in(pos2(2),pos1(2))*wgt1(2)*wgt2(2) 

      if (sdam_debug) then 
         write(*,*) "------------------------------------------------"
         write(*,*) mlatin, mltin, pos1, wgt1, pos2, wgt2, val_out
         write(*,*) "------------------------------------------------"
      end if

    end subroutine interp_sdam

  end module ModSDAM
  
