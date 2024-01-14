! Module using self-defined auroral inputs 
! Created: Qingyu Zhu, 02/25/2021

module ModSDFM

  use netcdf

  implicit none                                                                
                                                                               
  character(len=*), parameter :: &                                             
       sdfm_pt='UA/DataIn/', &                     
       sdfm_fn_nh='fac_NH1.nc', &                                 
       sdfm_fn_sh='fac_SH1.nc'                                    
                                                                               
  integer :: sdfm_nmlt, sdfm_nmlat, sdfm_ntime
                               
  real, allocatable :: sdfm_mlts(:),sdfm_colats(:)                          
                                                                               
  integer, allocatable :: sdfm_year(:), sdfm_mon(:), sdfm_day(:), &            
       sdfm_hour(:), sdfm_minute(:)

  real, allocatable :: sdfm_time(:)

  real, allocatable :: sdfm_fac_nh(:,:,:), sdfm_fac_sh(:,:,:)

  real, allocatable :: sdfm_facin_nh(:,:), sdfm_facin_sh(:,:)

  logical, parameter :: sdfm_debug = .false.

  contains

    subroutine read_sdfm

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
      istat=nf90_open(path=sdfm_pt//sdfm_fn_nh,mode=nf90_nowrite,ncid=ncid)
      if (sdfm_debug .and. (istat==nf90_noerr)) &
           write(*,*) " ------ NH SDFM file opened"

      ! Dimension                                                              
      istat = nf90_inquire(ncid, ndims, nvars)
      
      do i=1,ndims                                                             
         istat=nf90_inquire_dimension(ncid, i, dimname, len)                   

         if (sdfm_debug .and. (istat==nf90_noerr)) &
              write(*,*) i,trim(dimname)

         if (i==3) nfile=len                                                   
         if (i==1) nmlat=len                                                   
         if (i==2) nmlt=len                                                    
      end do                                                                   
                                                                               
      sdfm_ntime=nfile                                                         
      sdfm_nmlt=nmlt                                                           
      sdfm_nmlat=nmlat

      if (sdfm_debug) write(*,*) sdfm_ntime, sdfm_nmlt, sdfm_nmlat

      if (.not. allocated(sdfm_mlts)) allocate(sdfm_mlts(sdfm_nmlt),stat=istat)
      istat=nf90_get_var(ncid, 2, sdfm_mlts) 

      if (.not. allocated(sdfm_colats)) &
           allocate(sdfm_colats(sdfm_nmlat),stat=istat)
      istat=nf90_get_var(ncid, 1, sdfm_colats)

      if (.not. allocated(sdfm_year)) then                               
         allocate(sdfm_year(sdfm_ntime),stat=istat)                
         allocate(sdfm_mon(sdfm_ntime),stat=istat)                  
         allocate(sdfm_day(sdfm_ntime),stat=istat)                   
         allocate(sdfm_hour(sdfm_ntime),stat=istat)                  
         allocate(sdfm_minute(sdfm_ntime),stat=istat)
         allocate(sdfm_time(sdfm_ntime),stat=istat)
      end if

      istat=nf90_get_var(ncid, 4, sdfm_year)                           
      istat=nf90_get_var(ncid, 5, sdfm_mon)                    
      istat=nf90_get_var(ncid, 6, sdfm_day)                    
      istat=nf90_get_var(ncid, 7, sdfm_hour)                   
      istat=nf90_get_var(ncid, 8, sdfm_minute)

      do i=1,sdfm_ntime   

         iTime=0                                                               
         iTime(1)=sdfm_year(i)                                                 
         iTime(2)=sdfm_mon(i)                                                  
         iTime(3)=sdfm_day(i)                                                  
         iTime(4)=sdfm_hour(i)                                                 
         iTime(5)=sdfm_minute(i)

         call time_int_to_real(iTime,sdfm_time(i))

         if (sdfm_debug) write(*,*) "Input Time:", iTime, sdfm_time(i)

      end do

      if (.not. allocated(sdfm_fac_nh)) then                                  
         allocate(sdfm_fac_nh(sdfm_nmlt,sdfm_nmlat,sdfm_ntime),stat=istat)
      end if

      istat=nf90_get_var(ncid, 9, sdfm_fac_nh)
      if (sdfm_debug .and. (istat==nf90_noerr)) &
           write(*,*) "NH max fac for time 24", maxval(sdfm_fac_nh(:,:,24))
      
!      if (sdfm_debug .and. (istat==nf90_noerr)) &                             
!           write(*,*) "NH max nflux for time 2", maxval(sdfm_nfx_nh(:,:,2))

      istat=nf90_close(ncid)
      if (sdfm_debug .and. (istat==nf90_noerr)) &                     
           write(*,*) " ------ NH SDFM file closed"

      !!!!!! SH
      istat=nf90_open(path=sdfm_pt//sdfm_fn_sh,mode=nf90_nowrite,ncid=ncid)
      if (sdfm_debug .and. (istat==nf90_noerr)) &
           write(*,*) " ------ SH SDFM file opened"

      if (.not. allocated(sdfm_fac_sh)) then                                  
         allocate(sdfm_fac_sh(sdfm_nmlt,sdfm_nmlat,sdfm_ntime),stat=istat)
      end if

      istat=nf90_get_var(ncid, 9, sdfm_fac_sh)
      if (sdfm_debug .and. (istat==nf90_noerr)) &
           write(*,*) "NH max fac for time 2", maxval(sdfm_fac_nh(:,:,24))

!      if (sdfm_debug .and. (istat==nf90_noerr)) &
!           write(*,*) "NH max nflux for time 2", maxval(sdfm_nfx_nh(:,:,2))


      istat=nf90_close(ncid)
      if (sdfm_debug .and. (istat==nf90_noerr)) &                     
           write(*,*) " ------ SH SDFM file closed"

      deallocate(sdfm_year)                                                    
      deallocate(sdfm_mon)                                                     
      deallocate(sdfm_day)                                                     
      deallocate(sdfm_hour)                                                    
      deallocate(sdfm_minute)                                                

    end subroutine read_sdfm

    ! -------------------------------------------------------------------------

    subroutine get_currenttime_sdfm
      
      use ModTime, only: CurrentTime
      
      implicit none                                                            
                                                                               
      real :: wgt1, wgt2                                                      
      integer :: left, right, i, istat                                         
      
      left=1                                                                  
      right=1                                                                 
      wgt1=0.                                                                  
      wgt2=0.

      if (.not. allocated(sdfm_facin_nh)) then                               
         allocate(sdfm_facin_nh(sdfm_nmlt,sdfm_nmlat),stat=istat)        
         allocate(sdfm_facin_sh(sdfm_nmlt,sdfm_nmlat),stat=istat)     
      end if

      if (CurrentTime<sdfm_time(1)) then
         left=1                                                                
         right=1                                                               
         wgt1=0.                                                               
         wgt2=1.                                                               
      else if (CurrentTime>=sdfm_time(sdfm_ntime)) then             
         left=sdfm_ntime                                          
         right=sdfm_ntime                              
         wgt1=1.                                                               
         wgt2=0.
      else                                                                  
         do i=1,sdfm_ntime-1                                                  
            if (CurrentTime>=sdfm_time(i) .and. &                        
                 CurrentTime<sdfm_time(i+1)) then               
               left=i                                                          
               right=i+1                                                       
               wgt1=(sdfm_time(i+1)-CurrentTime)/&                        
                    (sdfm_time(i+1)-sdfm_time(i))                 
               wgt2=1-wgt1                                                     
               exit                                                            
            end if
         end do
      end if
          
      if (sdfm_debug) then
         if ((wgt1+wgt2)<1.) then                                              
            write(*,*) "Time is not found, check !!!"                          
         else                                                                  
            write(*,*) "Found data for", CurrentTime, left, right, wgt1, wgt2
         end if
      end if    
      
      ! FAC
      sdfm_facin_nh(:,:) = sdfm_fac_nh(:,:,left) * wgt1 + &
           sdfm_fac_nh(:,:,right) * wgt2 
      
      sdfm_facin_sh(:,:) = sdfm_fac_sh(:,:,left) * wgt1 + &
           sdfm_fac_sh(:,:,right) * wgt2        
     

      if (sdfm_debug) then
         write(*,*) "NH max FAC", maxval(sdfm_facin_nh)
         write(*,*) "SH max FAC", maxval(sdfm_facin_sh) 
      end if

    end subroutine get_currenttime_sdfm
    
    ! Spatial interpolation
    ! -----------------------------------------------------------------------
    subroutine interp_sdfm(mlatin,mltin,colats,mlts,val_in,val_out)
      
      implicit none
      
      real, intent(in) :: mlatin, mltin                                   
      real, intent(in) :: colats(sdfm_nmlat), mlts(sdfm_nmlt)               
      real, intent(in) :: val_in(sdfm_nmlt,sdfm_nmlat)                      
      real, intent(out) :: val_out
      
      integer :: i, j, pos1(2),pos2(2)                                      
      real :: wgt1(2),wgt2(2), colatin
      
      pos1=1                                                                   
      pos2=1                                                                   
      wgt1=0.                                                                  
      wgt2=0.                                                                  
      val_out=0

      !!! MLAT
      ! Different from SDFM where the MLAT is decreasing
      ! SDFM uses increasing MLAT grids 

      colatin=90-mlatin

      if (colatin<colats(1)) then
         wgt1=0.                                                               
         pos1=1
      else if (colatin>=colats(sdfm_nmlat)) then
         wgt1=0.
         pos1=sdfm_nmlat
      else
         do i=1,sdfm_nmlat-1
            if ((colatin>=colats(i)) .and. (colatin<colats(i+1))) then
               pos1(1)=i                                                   
               pos1(2)=i+1
               wgt1(1)=(colats(i+1)-colatin)/(colats(i+1)-colats(i))
               wgt1(2)=1.-wgt1(1)
               exit
            end if
         end do
         
         if ((sum(wgt1)<0.5) .and. (mlatin<90.)) write(*,*) &
              "MLAT is not found, check !!!", mlatin,mltin
      end if
  
      ! MLT (1 h grid)
      do j=1,sdfm_nmlt-1                                                  
         if ((mltin>=mlts(j)) .and. (mltin<mlts(j+1))) then
            pos2(1)=j                                 
            pos2(2)=j+1                                       
            wgt2(1)=(mlts(j+1)-mltin)/(mlts(j+1)-mlts(j))           
            wgt2(2)=1-wgt2(1)
            exit
         end if
      end do
            
      ! Get the output
      val_out = val_in(pos2(1),pos1(1))*wgt1(1)*wgt2(1) + &                    
           val_in(pos2(1),pos1(2))*wgt1(2)*wgt2(1) + &                         
           val_in(pos2(2),pos1(1))*wgt1(1)*wgt2(2) + &                         
           val_in(pos2(2),pos1(2))*wgt1(2)*wgt2(2) 

      if (sdfm_debug) then 
         write(*,*) "------------------------------------------------"
         write(*,*) mlatin, mltin, pos1, wgt1, pos2, wgt2, val_out
         write(*,*) "------------------------------------------------"
      end if

    end subroutine interp_sdfm

  end module ModSDFM
  
