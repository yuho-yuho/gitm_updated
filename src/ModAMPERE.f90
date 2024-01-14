! AMPERE FAC MODULE
! NEED NETCDF MODULE, IF NOT LOADED ALREADY
! CREADED: QINGYU ZHU, 11/23/2020
! -----------------------------------------------------------------------------

module ModAMPERE

  use netcdf

  implicit none 

  character(len=*), parameter :: &
       ampr_pt='UA/DataIn/', &
       ampr_fn_nh='mar2013_north.nc', &
       ampr_fn_sh='mar2013_south.nc'
  
  integer :: ampr_nmlt, ampr_nmlat, ampr_ntime

  ! Should be the same in the NH and SH
  integer, allocatable :: &
       ampr_startyr(:), ampr_endyr(:), & ! year                                
       ampr_startmo(:), ampr_endmo(:), & ! month                               
       ampr_startdy(:), ampr_enddy(:), & ! day                                 
       ampr_starthr(:), ampr_endhr(:), & ! hour                                
       ampr_startmt(:), ampr_endmt(:), & ! minute                              
       ampr_startsc(:), ampr_endsc(:) ! second

  real, allocatable :: ampr_starttime(:), ampr_endtime(:)

  real, allocatable, dimension(:,:) :: &
       ampr_mlts_nh, ampr_mlats_nh, & ! NH
       ampr_mlts_sh, ampr_mlats_sh, & ! SH
       ampr_facin_nh, ampr_facin_sh

  real, allocatable, dimension(:,:,:) :: ampr_fac_nh, ampr_fac_sh

  logical, parameter :: ampr_debug=.true.

  contains
  ! --------------------------------------------------------------------------

    subroutine read_ampere

      !use ModTimeConvert, ONLY: time_int_to_real 
      
      implicit none 

      integer :: ncid, istat
      integer :: nDims, nVars
      integer :: varid, dimid
      integer :: i, j, k
      integer :: nblock, nvector

      integer, allocatable :: var_1d(:), mlts1(:), mlats1(:)
      real, allocatable :: var_2d(:,:)

      integer :: iTime(7)

      !!!!!! Read in NH inputs
      ! -----------------------------------------------------------------------
      !!! OPEN FILE
      istat=nf90_open(path=ampr_pt//ampr_fn_nh,mode=nf90_nowrite,ncid=ncid)
      if (istat == nf90_noerr) write(*,*) "NH File opened" 

      !!! Get dimension
      istat = nf90_inquire(ncid, ndims, nvars)

      ! NTIME
      istat=nf90_inquire_dimension(ncid, 1, len=nblock)
      istat=nf90_inquire_dimension(ncid, 2, len=nvector)
      ampr_ntime=nblock

      if (.not. allocated(var_1d)) allocate(var_1d(ampr_ntime),stat=istat)
      
      ! NMLT
      istat=nf90_get_var(ncid, 1, var_1d)                                  
      ampr_nmlt=var_1d(1)

      ! NMLAT
      istat=nf90_get_var(ncid, 2, var_1d)
      ampr_nmlat=var_1d(1)

      !!! TIME ARRAYS
      call init_ampr_time_arrays(ampr_ntime)
      
      istat=nf90_get_var(ncid, 3, var_1d)                      
      ampr_startyr=var_1d(:)                                          
      istat=nf90_get_var(ncid, 4, var_1d)                             
      ampr_startmo=var_1d(:)                                               
      istat=nf90_get_var(ncid, 5, var_1d)                          
      ampr_startdy=var_1d(:)                                          
                                                                               
      istat=nf90_get_var(ncid, 6, var_1d)                          
      ampr_starthr=var_1d(:)                                              
      istat=nf90_get_var(ncid, 7, var_1d)                                   
      ampr_startmt=var_1d(:)                                               
      istat=nf90_get_var(ncid, 8, var_1d)                                    
      ampr_startsc=var_1d(:)

      istat=nf90_get_var(ncid, 9, var_1d)                   
      ampr_endyr=var_1d(:)                                               
      istat=nf90_get_var(ncid, 10, var_1d)                               
      ampr_endmo=var_1d(:)                                               
      istat=nf90_get_var(ncid, 11, var_1d)                            
      ampr_enddy=var_1d(:)                                          
                                                                   
      istat=nf90_get_var(ncid, 12, var_1d)                       
      ampr_endhr=var_1d(:)                                        
      istat=nf90_get_var(ncid, 13, var_1d)                     
      ampr_endmt=var_1d(:)                                            
      istat=nf90_get_var(ncid, 14, var_1d)                          
      ampr_endsc=var_1d(:)
      
      !!! Convert time, only use start time
      if (.not. allocated(ampr_starttime)) then
         allocate(ampr_starttime(ampr_ntime),stat=istat)
         allocate(ampr_endtime(ampr_ntime),stat=istat)
      end if

      do i=1,ampr_ntime

         iTime(1)=ampr_startyr(i)
         iTime(2)=ampr_startmo(i)
         iTime(3)=ampr_startdy(i)
         iTime(4)=ampr_starthr(i)
         iTime(5)=ampr_startmt(i)
         iTime(6)=ampr_startsc(i)
         iTime(7)=0

         call time_int_to_real(iTime,ampr_starttime(i))

         !if (ampr_debug) write(*,*) ampr_starttime(i)

      end do

      call deallocate_ampr_time_arrays

      !!! MLAT and MLT
      if (.not. allocated(var_2d)) allocate(var_2d(nvector,nblock),stat=istat)

      if (.not. allocated(mlts1)) then
         allocate(mlts1(nvector),stat=istat)               
         allocate(mlats1(nvector),stat=istat)                   
      end if
  
      ! MLAT
      istat=nf90_get_var(ncid, 15, var_2d)
      mlats1(:)=90.-var_2d(:,1)

      ! MLT
      istat=nf90_get_var(ncid, 16, var_2d)
      mlts1(:)=var_2d(:,1)

      if (.not. allocated(ampr_mlts_nh)) then     
         allocate(ampr_mlts_nh(ampr_nmlat,ampr_nmlt),stat=istat)     
         allocate(ampr_mlats_nh(ampr_nmlat,ampr_nmlt),stat=istat)       
      end if

      ! Reshape
      ampr_mlts_nh=reshape(mlts1,(/ampr_nmlat,ampr_nmlt/))
      ampr_mlats_nh=reshape(mlats1,(/ampr_nmlat,ampr_nmlt/)) 

      !!! FAC
      if (.not. allocated(ampr_fac_nh)) &
           allocate(ampr_fac_nh(ampr_ntime,ampr_nmlat,ampr_nmlt),stat=istat) 

      istat=nf90_get_var(ncid, 27, var_2d)

      do i=1,nblock         
         ampr_fac_nh(i,:,:)=reshape(var_2d(:,i),(/ampr_nmlat,ampr_nmlt/))
      end do

      !!! CLOSE FILE
      istat=nf90_close(ncid)
      if (istat == nf90_noerr) write(*,*) "NH File closed" 

      !!!!!! Read in SH inputs
      !!!!!! Only need to read in MLAT, MLT and FAC
      ! -----------------------------------------------------------------------
      !!! OPEN FILE
      istat=nf90_open(path=ampr_pt//ampr_fn_sh,mode=nf90_nowrite,ncid=ncid)
      if (istat == nf90_noerr) write(*,*) "SH File opened" 

      ! MLAT
      istat=nf90_get_var(ncid, 15, var_2d)
      mlats1(:)=90.-var_2d(:,1)

      ! MLT
      istat=nf90_get_var(ncid, 16, var_2d)
      mlts1(:)=var_2d(:,1)

      if (.not. allocated(ampr_mlts_sh)) then     
         allocate(ampr_mlts_sh(ampr_nmlat,ampr_nmlt),stat=istat)     
         allocate(ampr_mlats_sh(ampr_nmlat,ampr_nmlt),stat=istat)       
      end if

      ! Reshape
      ampr_mlts_sh=reshape(mlts1,(/ampr_nmlat,ampr_nmlt/))
      ampr_mlats_sh=reshape(mlats1,(/ampr_nmlat,ampr_nmlt/)) 

      !!! FAC
      if (.not. allocated(ampr_fac_sh)) &
           allocate(ampr_fac_sh(ampr_ntime,ampr_nmlat,ampr_nmlt),stat=istat) 

      istat=nf90_get_var(ncid, 27, var_2d)

      do i=1,nblock         
         ampr_fac_sh(i,:,:)=reshape(var_2d(:,i),(/ampr_nmlat,ampr_nmlt/))
      end do

      !!! CLOSE FILE
      istat=nf90_close(ncid)
      if (istat == nf90_noerr) write(*,*) "SH File closed" 

      !!! Deallocate arrays
      if (allocated(mlts1)) then                                     
         deallocate(mlts1,stat=istat)                             
         deallocate(mlats1,stat=istat)                              
      end if

      if (allocated(var_1d)) deallocate(var_1d,stat=istat)
      if (allocated(var_2d)) deallocate(var_2d,stat=istat)

    end subroutine read_ampere

    ! Temporal Interpolation
    ! -------------------------------------------------------------------------
    subroutine get_currenttime_fac

      use ModTime, only: CurrentTime

      implicit none 

      real :: wgt1, wgt2
      integer :: left, right, i, istat

      left=1
      right=1
      wgt1=0.
      wgt2=0.

      if (.not. allocated(ampr_facin_nh)) then     
         allocate(ampr_facin_nh(ampr_nmlat,ampr_nmlt),stat=istat)     
         allocate(ampr_facin_sh(ampr_nmlat,ampr_nmlt),stat=istat)       
      end if
      
      if (CurrentTime<ampr_starttime(1)) then
         left=1
         right=1
         wgt1=0.
         wgt2=1.
      else if (CurrentTime>=ampr_starttime(ampr_ntime)) then
         left=ampr_ntime
         right=ampr_ntime
         wgt1=1.
         wgt2=0.
      else
         do i=1,ampr_ntime-1
            if (CurrentTime>=ampr_starttime(i) .and. &
                 CurrentTime<ampr_starttime(i+1)) then
               left=i
               right=i+1
               wgt1=(ampr_starttime(i+1)-CurrentTime)/&
                    (ampr_starttime(i+1)-ampr_starttime(i))
               wgt2=1-wgt1
               exit
            end if
         end do
      end if

      if (ampr_debug) then
         if ((wgt1+wgt2)<1.) then
            write(*,*) "Time is not found, check !!!"
         else
            write(*,*) "Found data for", CurrentTime, left, right
         end if
      end if

      ampr_facin_nh(:,:) = ampr_fac_nh(left,:,:) * wgt1 + &
           ampr_fac_nh(right,:,:) * wgt2

      ampr_facin_sh(:,:) = ampr_fac_sh(left,:,:) * wgt1 + &
           ampr_fac_sh(right,:,:) * wgt2

    end subroutine get_currenttime_fac

    ! Spatial interpolation
    ! -------------------------------------------------------------------------
    subroutine interp_fac(mlatin,mltin,mlats,mlts,fac_in,fac_out)

      implicit none 

      real, intent(in) :: mlatin, mltin
      real, intent(in) :: mlats(ampr_nmlat), mlts(ampr_nmlt)
      real, intent(in) :: fac_in(ampr_nmlat,ampr_nmlt)
      real, intent(out) :: fac_out

      integer :: i, j, pos1(2),pos2(2)
      real :: wgt1(2),wgt2(2)

      ! mlats from 89 to 40
      ! mlts from 0 to 23
      
      pos1=1
      pos2=1
      wgt1=0.
      wgt2=0.
      fac_out=0

      ! MLAT
      if ((mlatin<=mlats(ampr_nmlat)) .or. (mlatin>mlats(1))) then
         wgt1=0.
         pos1=1
      else
         do i=1,ampr_nmlat-1

            if ((mlatin<=mlats(i)) .and. (mlatin>mlats(i+1))) then
               pos1(1)=i
               pos1(2)=i+1
               wgt1(1)=(mlatin-mlats(i+1))/(mlats(i)-mlats(i+1))
               wgt1(2)=1.-wgt1(1)
               exit
            end if
         end do

         if (sum(wgt1)<0.5) write(*,*) "MLAT is not found, check !!!"

      end if

      ! MLT
      if (mltin>=mlts(ampr_nmlt)) then
         pos2(1)=ampr_nmlt
         pos2(2)=1
         wgt2(1)=(24+mlts(1)-mltin)/(24+mlts(1)-mlts(ampr_nmlt))
         wgt2(2)=1-wgt2(1)
      else
         do j=1,ampr_nmlt-1

            if ((mltin>=mlts(j)) .and. (mltin<mlts(j+1))) then
               pos2(1)=j
               pos2(2)=j+1
               wgt2(1)=(mlts(j+1)-mltin)/(mlts(j+1)-mlts(j))
               wgt2(2)=1-wgt2(1)
               exit
            end if
         end do
      end if

      if (sum(wgt2)<0.5) write(*,*) "MLT is not found, check !!!"
      
      fac_out = fac_in(pos1(1),pos2(1))*wgt1(1)*wgt2(1) + &
           fac_in(pos1(1),pos2(2))*wgt1(1)*wgt2(2) + &
           fac_in(pos1(2),pos2(1))*wgt1(2)*wgt2(1) + &
           fac_in(pos1(2),pos2(2))*wgt1(2)*wgt2(2)
      

    end subroutine interp_fac



    ! Get AMPERE TIME array
    ! -------------------------------------------------------------------------
    subroutine init_ampr_time_arrays(ntime)

      implicit none

      integer, intent(in) :: ntime

      integer :: istat

      if (.not. allocated(ampr_startyr)) then                            
         allocate(ampr_startyr(ntime),stat=istat)              
         allocate(ampr_startmo(ntime),stat=istat)             
         allocate(ampr_startdy(ntime),stat=istat)                      
         allocate(ampr_starthr(ntime),stat=istat)                  
         allocate(ampr_startmt(ntime),stat=istat)             
         allocate(ampr_startsc(ntime),stat=istat)

         allocate(ampr_endyr(ntime),stat=istat)
         allocate(ampr_endmo(ntime),stat=istat)                          
         allocate(ampr_enddy(ntime),stat=istat)                          
         allocate(ampr_endhr(ntime),stat=istat)                        
         allocate(ampr_endmt(ntime),stat=istat)                          
         allocate(ampr_endsc(ntime),stat=istat)
      end if

    end subroutine init_ampr_time_arrays

    ! -------------------------------------------------------------------------
    subroutine deallocate_ampr_time_arrays

      implicit none

      integer :: istat

      if ( allocated(ampr_startyr)) then                            
         deallocate(ampr_startyr,stat=istat)              
         deallocate(ampr_startmo,stat=istat)             
         deallocate(ampr_startdy,stat=istat)                      
         deallocate(ampr_starthr,stat=istat)                  
         deallocate(ampr_startmt,stat=istat)             
         deallocate(ampr_startsc,stat=istat)

         deallocate(ampr_endyr,stat=istat)
         deallocate(ampr_endmo,stat=istat)                          
         deallocate(ampr_enddy,stat=istat)                          
         deallocate(ampr_endhr,stat=istat)                        
         deallocate(ampr_endmt,stat=istat)                          
         deallocate(ampr_endsc,stat=istat)
         deallocate(ampr_endtime,stat=istat)

      end if

    end subroutine deallocate_ampr_time_arrays

end module ModAMPERE
