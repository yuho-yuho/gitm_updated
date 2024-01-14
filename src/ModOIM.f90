! OIM AMPERE FAC MODULE
! NEED NETCDF MODULE, IF NOT LOADED ALREADY
! ADAPTED FROM MODAMPERE.F90
! CREATED: QINGYU ZHU, 12/01/2020
! -----------------------------------------------------------------------------
module ModOIM

  use netcdf

  implicit none 

  character(len=*), parameter :: &                   
       oim_pt='/glade/u/home/qingyu/code/oim/', &                             
       oim_fn_nh='201005_147_150_north.nc', &                   
       oim_fn_sh='201005_147_150_south.nc'                             
                                                                               
  integer :: oim_nmlt, oim_nmlat, oim_ntime

  real, allocatable :: oim_mlts(:),oim_mlats(:)

  integer, allocatable :: oim_year(:), oim_mon(:), oim_day(:), &
       oim_hour(:), oim_minute(:)

  real, allocatable :: oim_time(:)

  real, allocatable :: oim_fac_nh(:,:,:), oim_fac_sh(:,:,:)

  real, allocatable :: oim_facin_nh(:,:), oim_facin_sh(:,:)

  logical, parameter :: oim_debug = .false.

  contains
    ! -------------------------------------------------------------------------
    ! Read in NH/SH OIM FAC 
    subroutine read_oimfac

      implicit none 

      integer :: ncid, istat
      integer :: nDims, nVars
      integer :: varid, dimid 
      integer :: i, j, k
      integer :: nfile, nmlat, nmlt, len
      integer :: iTime(7) 
      real, allocatable, dimension(:) :: year, doy, hour, minute
      
      ! ----------------------------- NH -----------------------------
      ! OPEN FILE
      istat=nf90_open(path=oim_pt//oim_fn_nh,mode=nf90_nowrite,ncid=ncid)
      if (oim_debug .and. (istat == nf90_noerr)) &
           write(*,*) "------ NH FAC file opened"

      ! Get the dimension
      istat = nf90_inquire(ncid, ndims, nvars)
      if (oim_debug .and. (istat == nf90_noerr)) write(*,*) ndims, nvars
      
      istat=nf90_inquire_dimension(ncid, 1, len=nmlat)
      istat=nf90_inquire_dimension(ncid, 2, len=nmlt)
      istat=nf90_inquire_dimension(ncid, 3, len=nfile)

      if (oim_debug .and. (istat == nf90_noerr)) &
           write(*,*) "nfile, nmlt, nmlat",nfile,nmlt,nmlat

      oim_ntime=nfile                           
      oim_nmlt=nmlt                                   
      oim_nmlat=nmlat

      ! MLTS
      if (.not. allocated(oim_mlts)) allocate(oim_mlts(oim_nmlt),stat=istat)
      istat=nf90_get_var(ncid, 2, oim_mlts) 
      if (oim_debug .and. (istat == nf90_noerr)) &
           write(*,*) "MLT values obtained"

      ! MLATS
      if (.not. allocated(oim_mlats)) allocate(oim_mlats(oim_nmlat),stat=istat)
      istat=nf90_get_var(ncid, 1, oim_mlats)
      if (oim_debug .and. (istat == nf90_noerr)) &
           write(*,*) "MLAT values obtained"

      ! Time 
      if (.not. allocated(oim_year)) then
         allocate(oim_year(oim_ntime),stat=istat)   
         allocate(oim_mon(oim_ntime),stat=istat)       
         allocate(oim_day(oim_ntime),stat=istat)           
         allocate(oim_hour(oim_ntime),stat=istat)                    
         allocate(oim_minute(oim_ntime),stat=istat)
         allocate(oim_time(oim_ntime),stat=istat)
         
         allocate(year(oim_ntime),stat=istat)                
         allocate(doy(oim_ntime),stat=istat)                    
         allocate(hour(oim_ntime),stat=istat)                          
         allocate(minute(oim_ntime),stat=istat)
      end if

      istat=nf90_get_var(ncid, 9, year)            
      istat=nf90_get_var(ncid, 8, doy)                                 
      istat=nf90_get_var(ncid, 10, hour)                               
      istat=nf90_get_var(ncid, 11, minute) 
      if (oim_debug .and. (istat == nf90_noerr)) &
           write(*,*) "Time values obtained"

      do i=1,oim_ntime
          oim_year(i)=int(year(i))                   
          oim_hour(i)=int(hour(i))                                  
          oim_minute(i)=int(minute(i))
          
          call doy_to_mon_day(int(year(i)),int(doy(i)),oim_mon(i),oim_day(i)) 
          
          iTime=0
          iTime(1)=oim_year(i)
          iTime(2)=oim_mon(i)
          iTime(3)=oim_day(i)
          iTime(4)=oim_hour(i)
          iTime(5)=oim_minute(i)

          call time_int_to_real(iTime,oim_time(i))

          if (oim_debug .and. (i==1)) write(*,*) iTime,oim_time(1)

      end do

      if (allocated(year)) then
         deallocate(year,stat=istat)
         deallocate(doy,stat=istat)
         deallocate(hour,stat=istat)
         deallocate(minute,stat=istat)
      end if
      
      ! FAC 
      if (.not. allocated(oim_fac_nh)) then 
         allocate(oim_fac_nh(oim_nmlt,oim_nmlat,oim_ntime),stat=istat)
      end if
       if (oim_debug .and. (istat == nf90_noerr)) &
            write(*,*) "FAC NH array allocated"

      istat=nf90_get_var(ncid, 4, oim_fac_nh)
      if (oim_debug) write(*,*) "NH FAC read status:",istat
      if (oim_debug .and. (istat == nf90_noerr)) &
           write(*,*) "FAC NH values obtained"

      ! CLOSE FILE
      istat=nf90_close(ncid)
      if (oim_debug .and. (istat == nf90_noerr)) &
           write(*,*) "------ NH FAC file closed"

      ! ----------------------------- SH -----------------------------
      ! Only read FAC

      ! OPEN FILE
      istat=nf90_open(path=oim_pt//oim_fn_sh,mode=nf90_nowrite,ncid=ncid)
      if (oim_debug .and. (istat == nf90_noerr)) &
           write(*,*) "------ SH FAC file opened"

      ! FAC 
      if (.not. allocated(oim_fac_sh)) then 
         allocate(oim_fac_sh(oim_nmlt,oim_nmlat,oim_ntime),stat=istat)
      end if

      istat=nf90_get_var(ncid, 4, oim_fac_sh)      
      if (oim_debug .and. (istat == nf90_noerr)) &
           write(*,*) "FAC SH values obtained"

      ! CLOSE FILE
      istat=nf90_close(ncid)
      if (oim_debug .and. (istat == nf90_noerr)) &
           write(*,*) "------ SH FAC file closed"

    end subroutine read_oimfac

    ! TEMPORAL INTERPOLATION
    ! NOTE DIFFERNET SHAPE OF OIM_FAC FROM AMPR_FAC
    ! -------------------------------------------------------------------------
    subroutine get_currenttime_oimfac

      use ModTime, only: CurrentTime
     
      implicit none                                                            
                                                                               
      real :: wgt1, wgt2                                                       
      integer :: left, right, i, istat                                         
                                                                               
      left=1                                                                   
      right=1                                                                  
      wgt1=0.                                                                  
      wgt2=0.

      if (.not. allocated(oim_facin_nh)) then                                 
         allocate(oim_facin_nh(oim_nmlt,oim_nmlat),stat=istat)              
         allocate(oim_facin_sh(oim_nmlt,oim_nmlat),stat=istat)              
      end if 

      if (CurrentTime<oim_time(1)) then                                  
         left=1                                                                
         right=1                                                               
         wgt1=0.                                                               
         wgt2=1.                                                               
      else if (CurrentTime>=oim_time(oim_ntime)) then                   
         left=oim_ntime                                                       
         right=oim_ntime                                                      
         wgt1=1.                                                               
         wgt2=0.
      else
         do i=1,oim_ntime-1                                                   
            if (CurrentTime>=oim_time(i) .and. &                         
                 CurrentTime<oim_time(i+1)) then                         
               left=i                                                          
               right=i+1                                                       
               wgt1=(oim_time(i+1)-CurrentTime)/&                        
                    (oim_time(i+1)-oim_time(i))                    
               wgt2=1-wgt1                                                     
               exit                                                            
            end if                                                             
         end do                                                                
      end if

      if (oim_debug) then                                                     
         if ((wgt1+wgt2)<1.) then                                              
            write(*,*) "Time is not found, check !!!"                          
         else                                                                  
            write(*,*) "Found data for", CurrentTime, left, right              
         end if                                                                
      end if                                                                   
                                                                               
      oim_facin_nh(:,:) = oim_fac_nh(:,:,left) * wgt1 + &                    
           oim_fac_nh(:,:,right) * wgt2                                       
                                                                               
      oim_facin_sh(:,:) = oim_fac_sh(:,:,left) * wgt1 + &                    
           oim_fac_sh(:,:,right) * wgt2

    end subroutine get_currenttime_oimfac
    
    ! SPATIAL INTERPOLATION
    ! -------------------------------------------------------------------------
    subroutine interp_oimfac(mlatin,mltin,mlats,mlts,fac_in,fac_out)      
                                                                               
      implicit none                                                            
                                                                               
      real, intent(in) :: mlatin, mltin                                        
      real, intent(in) :: mlats(oim_nmlat), mlts(oim_nmlt)                   
      real, intent(in) :: fac_in(oim_nmlt,oim_nmlat)                         
      real, intent(out) :: fac_out                                             
                                                                               
      integer :: i, j, pos1(2),pos2(2)                                         
      real :: wgt1(2),wgt2(2)    

      pos1=1                                                                   
      pos2=1                                                                   
      wgt1=0.                                                                  
      wgt2=0.                                                                  
      fac_out=0                                                                
                                                                               
      ! MLAT                                                                   
      if (mlatin<=mlats(oim_nmlat)) then             
         wgt1=0.                                                               
         pos1=1
      else
         do i=1,oim_nmlat-1                                                   
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
      do j=1,oim_nmlt-1
         if ((mltin>=mlts(j)) .and. (mltin<mlts(j+1))) then                  
            pos2(1)=j                                                        
            pos2(2)=j+1                                                      
            wgt2(1)=(mlts(j+1)-mltin)/(mlts(j+1)-mlts(j))                    
            wgt2(2)=1-wgt2(1)                                                
            exit                                                             
         end if
      end do

      if (sum(wgt2)<0.5) write(*,*) "MLT is not found, check !!!" 

      ! Note different shape of oim_fac (nmlt,nmlat) from ampr_fac(nmlat,nmlt)
      fac_out = fac_in(pos2(1),pos1(1))*wgt1(1)*wgt2(1) + &                    
           fac_in(pos2(1),pos1(2))*wgt1(2)*wgt2(1) + &                         
           fac_in(pos2(2),pos1(1))*wgt1(1)*wgt2(2) + &                         
           fac_in(pos2(2),pos1(2))*wgt1(2)*wgt2(2)

    end subroutine interp_oimfac

end module ModOIM

! GET MONTH AND DAY ACCORDING YEAR AND DOY
! -----------------------------------------------------------------------------
subroutine doy_to_mon_day(year,doy,mon,day)                                    
                                                                               
  implicit none                                                                
                                                                               
  integer, intent(in) :: year, doy                                             
  integer, intent(out) :: mon, day                  
  integer, dimension(12) :: dayofmon                                           
  integer :: days, imon, iday                           
  logical :: isleap                                                            
                                                                               
  isleap=.false.                                                               
                                                                               
  if (((mod(year,4) .eq. 0) .and. (mod(year,100) .ne. 0)) .or. &               
       (mod(year,400) .eq. 0)) &                                               
       isleap=.true. 

  dayofmon(1) = 31                                                             
  dayofmon(2) = 28                                                             
  dayofmon(3) = 31                                                             
  dayofmon(4) = 30                                                             
  dayofmon(5) = 31                                                             
  dayofmon(6) = 30                                                             
  dayofmon(7) = 31                                                             
  dayofmon(8) = 31                                                             
  dayofmon(9) = 30                                                             
  dayofmon(10) = 31                                                            
  dayofmon(11) = 30                                                            
  dayofmon(12) = 31                                                            
                                                                               
  if (isleap) dayofmon(2) = 29

  mon=1                                                                        
  day=doy                                                                      
                                                                               
  do while (day .gt. dayofmon(mon))                                            
     day=day-dayofmon(mon)                                                     
     mon=mon+1                                                                 
  end do

end subroutine doy_to_mon_day
