! GEdy subroutine series
! Gather field from different processors seperately and send to the root 
! processor
! 
! Created: Qingyu Zhu, 07/27/2017, qingyu.zhu@mavs.uta.edu
!
! Add in this GITM on 03/02/2020, qingyu  
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine gedy_gatherfield

  ! Global
  ! Coordinates from 3d Current module
  use params_module, only: nmlon, nmlat_T1, nmlat_h, dtr, h0, &
       nmlat_T2, hgt_fix, rtd, nmlatS2_h

  ! Fields in root processor
  use ModGedy, only:  & 
       un_gitm, vn_gitm, sigP_gitm, sigH_gitm, &
       lon_gitm, lat_gitm, hgt_gitm, mlon_gitm, mlat_gitm, &
       nlon_gitm, nlat_gitm, nhgt_gitm, &
       tasks

  ! MPI
  use ModMPi

  ! Fields need to be gathered from GITM 
  use ModGITM, only: iCommGITM, iProc, nProcs, &
       Longitude, Latitude, Altitude_GB, &
       MLongitude, MLatitude, &
       Velocity, Potential!, Potential_use
  use ModElectrodynamics, only: Sigma_Pedersen, Sigma_Hall

  ! Coordinates from GITM
  use ModSizeGitm, only: nLons, nLats, nALts
  use ModInputs, only: nBlocksLon, nBlocksLat

  ! Constants
  use ModConstants, only: pi

  ! --------------------------------------------------------------------------

  implicit none

  ! Local
  integer :: iError, i, j, k, ifld, req, iLon, iLat, itask
  integer :: pos, pos1, pos2
  integer, parameter :: nfld=4, root=0, tag=123, len_task_type=2
  integer, dimension(mpi_status_size) :: status

  real :: dLon, dLat, lon, lat
  real, allocatable :: send(:,:,:,:), recv(:,:,:,:,:), &
       send_tmp(:,:,:), recv_tmp(:,:,:,:)
  real, allocatable, save :: itasks_send(:,:), itasks_recv(:,:)
  
  ! --------------------------------------------------------------------------
  ! Main part
  ! --------------------------------------------------------------------------

  ! Build the geographic coordinates in the root task
  ! Longitude and Latitude are in deg
  if (iProc .eq. root) then

     ! Altitude                                        
     nhgt_gitm = nAlts
     if (.not. allocated(hgt_gitm)) then
        allocate(hgt_gitm(nAlts),stat=iError)
     endif
     hgt_gitm = Altitude_GB(1,1,1:nAlts,1)/1000.

     ! longitude                                        
     nlon_gitm=nLons*nBlocksLon
     if (.not. allocated(lon_gitm)) then
        allocate(lon_gitm(nLons*nBlocksLon),stat=iError)
     endif

     dLon = 360./nlon_gitm
     do i=1,nlon_gitm
        lon_gitm(i)=dLon*(i-1)+dLon/2.
     enddo

     ! latitude                                              
     nlat_gitm=nLats*nBlocksLat
     if (.not. allocated(lat_gitm)) then
        allocate(lat_gitm(nLats*nBlocksLat),stat=iError)
     endif

     dLat = 180./nlat_gitm
     do i=1,nlat_gitm
        lat_gitm(i)=dLat*(i-1)+dLat/2.-90.
     enddo

  endif

  ! Allocate the task for each processor and share them by using alltoall 
  if (.not. allocated(tasks)) then
     allocate(tasks(0:nProcs),stat=iError)
  endif

  do i=0,nProcs-1
     tasks(i)%lon0 = Longitude(1,1)
     tasks(i)%lat0 = Latitude(1,1)
  enddo

  !!!!!!!! Set up the send and receive array
  if (.not. allocated(itasks_send)) then
     allocate(itasks_send(2,0:nProcs-1),stat=iError)
     allocate(itasks_recv(2,0:nProcs-1),stat=iError)
  endif


  do i=0,nProcs-1
     itasks_send(1,i) = tasks(i)%lon0
     itasks_send(2,i) = tasks(i)%lat0
  enddo

  call mpi_alltoall(itasks_send,len_task_type,mpi_real,&
       itasks_recv,len_task_type,mpi_real,iCommGITM,iError)

  do i=0,nProcs-1
     tasks(i)%lon0 = itasks_recv(1,i)
     tasks(i)%lat0 = itasks_recv(2,i)
  enddo

  ! Send the fields to root processor
  ! exclude the ghost cells

  if (.not. allocated(send_tmp)) then
     allocate(send_tmp(nLons,nLats,nAlts),stat=iError)
  endif

  if (.not. allocated(send)) then
     allocate(send(nLons,nLats,nAlts,nfld),stat=iError)
  endif

  !!! un                                                   
  send_tmp=Velocity(1:nLons,1:nLats,1:nAlts,1,1)
  send(:,:,:,1) = send_tmp

  !!! vn                                                              
  send_tmp=Velocity(1:nLons,1:nLats,1:nAlts,2,1)
  send(:,:,:,2) = send_tmp

  !!! SigP                                                           
  send_tmp=Sigma_Pedersen(1:nLons,1:nLats,1:nAlts)
  send(:,:,:,3) = send_tmp

  !!! SigH                                                               
  send_tmp=Sigma_Hall(1:nLons,1:nLats,1:nAlts)
  send(:,:,:,4) = send_tmp

  !!!!!!!! Send
  call mpi_isend(send,nLons*nLats*nAlts*nfld,mpi_real,root,&
       tag,iCommGITM,req,iError)

  !----------------------------------------------------------------------------
  !!!!!!!! Check
  !write(*,*) 
  !write(*,*) "1sendfromProc #",iProc,send(6,6,24,1)
  !write(*,*) "2sendfromProc #",iProc,send(6,6,24,2)
  !write(*,*) "3sendfromProc #",iProc,send(6,6,24,3)
  !write(*,*) "4sendfromProc #",iProc,send(6,6,24,4)
  !----------------------------------------------------------------------------

  !!!!!!!! Receive
  
  if (iProc .eq. root) then

     if (.not. allocated(recv)) then
        allocate(recv(nLons,nLats,nAlts,nfld,0:nProcs-1),stat=iError)
     endif

     do i=0,nProcs-1
        call mpi_recv(recv(:,:,:,:,i),nLons*nLats*nAlts*nfld,mpi_real,i,&
             tag,iCommGITM,status,iError)

        !----------------------------------------------------------------------
        !!!!!!!! Check
        !write(*,*) 
        !write(*,*) "1fromProc #",i,recv(6,6,24,1,i)
        !write(*,*) "2fromProc #",i,recv(6,6,24,2,i)
        !write(*,*) "3fromProc #",i,recv(6,6,24,3,i)
        !write(*,*) "4fromProc #",i,recv(6,6,24,4,i)
        !----------------------------------------------------------------------
        

     enddo


  endif

  call mpi_barrier(iCommGITM,iError)

  ! --------------------------------------------------------------------------
  ! write the collected fields to a global matrix

  if (iProc .eq. root) then

     if(.not. allocated(recv_tmp)) then
        allocate(recv_tmp(nlon_gitm,nlat_gitm,nhgt_gitm,nfld),&
             stat=iError)
     endif

     if(.not. allocated(un_gitm)) then
        allocate(un_gitm(nlon_gitm,nlat_gitm,nhgt_gitm),stat=iError)
        allocate(vn_gitm(nlon_gitm,nlat_gitm,nhgt_gitm),stat=iError)
        allocate(sigP_gitm(nlon_gitm,nlat_gitm,nhgt_gitm),stat=iError)
        allocate(sigH_gitm(nlon_gitm,nlat_gitm,nhgt_gitm),stat=iError)
     endif

     do i=0,nProcs-1
        lon=tasks(i)%lon0*180/pi
        lat=tasks(i)%lat0*180/pi

        dLon = 360./nlon_gitm
        dLat = 180./nlat_gitm

        pos1=nint((lon-dLon/2.)/(dLon*nLons))*nLons
        pos2=nint((lat-dLat/2.+90.)/(dLat*nLats))*nLats

        do ifld=1,nfld
           recv_tmp(pos1+1:pos1+nLons,pos2+1:pos2+nLats,:,ifld) = &
                recv(1:nLons,1:nLats,:,ifld,i)
        enddo
     enddo

     ! unpack the recv_temp                                                
     un_gitm(:,:,:)=recv_tmp(:,:,:,1)
     vn_gitm(:,:,:)=recv_tmp(:,:,:,2)
     sigP_gitm(:,:,:)=recv_tmp(:,:,:,3)
     sigH_gitm(:,:,:)=recv_tmp(:,:,:,4)

     deallocate(recv_tmp,stat=iError)

  endif


end subroutine gedy_gatherfield

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine gedy_sendmagcoordinate

  ! Send magnetic coordinate from different processors
  ! Seperate from the gedy_gatherfield since do not need to redistribute

  ! Global
  use ModMpi
  use ModGITM
  use ModGedy, only : mlon_gitm, mlat_gitm

  implicit none

  ! Local
  integer :: i, j, k, iError, req1, req2, iLon, iLat
  integer, parameter :: root = 0, tag1=123, tag2=234
  integer, dimension(mpi_status_size) :: status

  real, allocatable :: send(:,:,:), recv(:,:,:,:)

  ! ---------------------------------------------------------------------------
  ! Main

  !!! Magnetic Longitude

  ! Allocate sending arrays
  if (.not. allocated(send)) then
     allocate(send(nLons,nLats,nAlts),stat=iError)
  endif
  send=MLongitude(1:nLons,1:nLats,1:nAlts,1)

  ! Send and receive
  call mpi_isend(send,nLons*nLats*nAlts,mpi_real,root,&
       tag1,iCommGITM,req1,iError)

  !----------------------------------------------------------------------------
  ! Check
  !write(*,*) "Proc#",iProc,"send mlon",send(3,9,36)
  !----------------------------------------------------------------------------

  if(iProc .eq. root) then
     if (.not. allocated(recv)) then
        allocate(recv(nLons,nLats,nAlts,0:nProcs-1),stat=iError)
     endif

     do i=0,nProcs-1
        call mpi_recv(recv(:,:,:,i),nLons*nLats*nAlts,mpi_real,&
             i,tag1,iCommGITM,status,iError)
        
        !----------------------------------------------------------------------
        ! Check
        !write(*,*) "Receive mlon",recv(3,9,36,i),"from Proc#", i
        !----------------------------------------------------------------------

     enddo

     if (.not. allocated(mlon_gitm)) then
        allocate(mlon_gitm(nLons,nLats,nAlts,0:nProcs-1),stat=iError)
     endif

     mlon_gitm = recv

     !----------------------------------------------------------------------
     ! Check
     !do i=0,nProcs-1
     !   write(*,*) "Mlon_gitm receive",mlon_gitm(3,9,36,i),"from Proc#", i
     !enddo
     !----------------------------------------------------------------------

     deallocate(recv,stat=iError)
  endif

  call mpi_barrier(iCommGITM,iError)
  deallocate(send,stat=iError)

  !!! Magnetic Latitude
  
  ! Allocate sending arrays
  if (.not. allocated(send)) then
     allocate(send(nLons,nLats,nAlts),stat=iError)
  endif
  send=MLatitude(1:nLons,1:nLats,1:nAlts,1)

  ! Send and receive

  call mpi_isend(send,nLons*nLats*nAlts,mpi_real,root,&
       tag2,iCommGITM,req1,iError)

  !----------------------------------------------------------------------------
  ! Check
  !write(*,*) "Proc#",iProc,"send mlat",send(3,9,36)
  !----------------------------------------------------------------------------

  if(iProc .eq. root) then
     if (.not. allocated(recv)) then
        allocate(recv(nLons,nLats,nAlts,0:nProcs-1),stat=iError)
     endif

     do i=0,nProcs-1
        call mpi_recv(recv(:,:,:,i),nLons*nLats*nAlts,mpi_real,&
             i,tag2,iCommGITM,status,iError)

        !----------------------------------------------------------------------
        ! Check
        !write(*,*) "Receive mlat",recv(3,9,36,i),"from Proc#", i
        !----------------------------------------------------------------------

     enddo

     if (.not. allocated(mlat_gitm)) then
        allocate(mlat_gitm(nLons,nLats,nAlts,0:nProcs-1),stat=iError)
     endif

     mlat_gitm = recv

     !----------------------------------------------------------------------
     ! Check
     !do i=0,nProcs-1
     !   write(*,*) "Mlat_gitm receive",mlat_gitm(3,9,36,i),"from Proc#", i
     !enddo
     !----------------------------------------------------------------------

     deallocate(recv,stat=iError)
  endif

  call mpi_barrier(iCommGITM,iError)
  deallocate(send,stat=iError)

end subroutine gedy_sendmagcoordinate
