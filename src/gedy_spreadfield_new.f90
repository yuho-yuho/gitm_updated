! New spread routine in order to avoid allocating too large matrix
! Created: Qingyu Zhu, 03/01/2018, qingyu.zhu@mavs.uta.edu
!
! Spread the mlat_3d, mlon_3d, poten_2d, ed1_s1, ed2_s1
! from the root processor to the remaining processor
!
! Add in this GITM on 03/02/2020, qingyu  
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine gedy_spreadfield_new

  ! Global
  use ModGedy, only: nmlat_T1, nmlon, &
       mlon_3d => ylonm2_deg, &
       mlat_3d => ylatm1_deg, &
       poten_2d, ed1_S1, ed2_S1, &
       fac_2d, & ! qingyu, 11/24/2020
       fac1_2d, dfac_2d, & ! qingyu, 11/25/2020
       poten1_2d, & ! qingyu, 11/26/2020
       Jwd_2d

  use ModGITM, only: iCommGITM, iProc, nProcs
  use ModMpi

  implicit none

  ! Local
  integer :: i, j, k, iError, req, req1, itask
  integer :: nmlat_3d, nmlon_3d
  integer, parameter :: nfld=5, root=0, tag1=111, tag2=222
  integer :: tags(nfld), tags1(nfld)
  integer, dimension(mpi_status_size) :: status

  real, allocatable :: send(:,:), recv(:,:)

  ! ---------------------------------------------------------------------------
  ! First send the coordinate 
  
  if (iProc .eq. root) then
     
     do itask=1,nProcs-1
        call mpi_send(mlon_3d,nmlon,mpi_real,itask,tag1,iCommGITM,iError)
        call mpi_send(mlat_3d,nmlat_T1,mpi_real,itask,tag2,iCommGITM,iError)
     end do

     !write(*,*) "Root sent to ", itask, "With", mlon_3d(10),mlat_3d(10)

  else

     call mpi_recv(mlon_3d,nmlon,mpi_real,root,tag1,&
          iCommGITM,status,iError)
     call mpi_recv(mlat_3d,nmlat_T1,mpi_real,root,tag2,&
          iCommGITM,status,iError)

     !write(*,*) iProc, "received With", mlon_3d(10),mlat_3d(10)

  end if

  call mpi_barrier(iCommGITM,iError)

  ! Second, send the poten_2d, ed1_S1, ed2_S1, and fac_2d

  do i=1,nfld
     tags(i)=111*(i+2)
  end do

  do i=1,nfld
     tags1(i)=121*(i+2)
  end do

  if (iProc .eq. root) then

     do itask=1,nProcs-1

        call mpi_send(poten_2d,nmlon*nmlat_T1,mpi_real,itask,&
             tags(1),iCommGITM,iError)
        call mpi_send(ed1_S1,nmlon*nmlat_T1,mpi_real,itask,&
             tags(2),iCommGITM,iError)
        call mpi_send(ed2_S1,nmlon*nmlat_T1,mpi_real,itask,&
             tags(3),iCommGITM,iError)

        ! qingyu, 11/25/2020
        call mpi_send(fac_2d,nmlon*nmlat_T1,mpi_real,itask,&
             tags(4),iCommGITM,iError)
        call mpi_send(fac1_2d,nmlon*nmlat_T1,mpi_real,itask,&
             tags(5),iCommGITM,iError)
        call mpi_send(dfac_2d,nmlon*nmlat_T1,mpi_real,itask,&
             tags1(1),iCommGITM,iError)
        call mpi_send(poten1_2d,nmlon*nmlat_T1,mpi_real,itask,&
             tags1(2),iCommGITM,iError)
        call mpi_send(Jwd_2d,nmlon*nmlat_T1,mpi_real,itask,&
             tags1(3),iCommGITM,iError)

        !write(*,*) "Root sent to ", itask, "With", poten_2d(50,10), &
        !     ed1_S1(50,10), ed2_S1(50,10)
        !write(*,*)


     end do

  else
     
     if (.not. allocated(poten_2d)) then
        
        allocate(poten_2d(nmlon,nmlat_T1),stat=iError)
        allocate(ed1_S1(nmlon,nmlat_T1),stat=iError)
        allocate(ed2_S1(nmlon,nmlat_T1),stat=iError)

        ! qingyu, 11/25/2020
        allocate(fac_2d(nmlon,nmlat_T1),stat=iError) 
        allocate(fac1_2d(nmlon,nmlat_T1),stat=iError)
        allocate(dfac_2d(nmlon,nmlat_T1),stat=iError)
        allocate(poten1_2d(nmlon,nmlat_T1),stat=iError)
        allocate(Jwd_2d(nmlon,nmlat_T1),stat=iError)
        
     end if

     call mpi_recv(poten_2d,nmlon*nmlat_T1,mpi_real,root,&
             tags(1),iCommGITM,status,iError)
     call mpi_recv(ed1_S1,nmlon*nmlat_T1,mpi_real,root,&
             tags(2),iCommGITM,status,iError)
     call mpi_recv(ed2_S1,nmlon*nmlat_T1,mpi_real,root,&
             tags(3),iCommGITM,status,iError)

     ! qingyu, 11/25/2020
     call mpi_recv(fac_2d,nmlon*nmlat_T1,mpi_real,root,&
          tags(4),iCommGITM,status,iError)
     call mpi_recv(fac1_2d,nmlon*nmlat_T1,mpi_real,root,&
          tags(5),iCommGITM,status,iError)
     call mpi_recv(dfac_2d,nmlon*nmlat_T1,mpi_real,root,&
          tags1(1),iCommGITM,status,iError)
     call mpi_recv(poten1_2d,nmlon*nmlat_T1,mpi_real,root,&
          tags1(2),iCommGITM,status,iError)     
     call mpi_recv(Jwd_2d,nmlon*nmlat_T1,mpi_real,root,&
          tags1(3),iCommGITM,status,iError)     

     !write(*,*) iProc, "received With",  poten_2d(50,10), &
     !     ed1_S1(50,10), ed2_S1(50,10)
     !write(*,*)
     
  end if
  
  call mpi_barrier(iCommGITM,iError)

end subroutine gedy_spreadfield_new

! -----------------------------------------------------------------------------
! Now we have the 2d values from the root, we can do the interpolation 
! in each processor rather than in the root 

subroutine gedy_get3dvalues_new


  use ModGedy, only: nmlat_T1, nmlon, &
       mlon_3d => ylonm2_deg, &
       mlat_3d => ylatm1_deg, &
       poten_2d, ed1_S1, ed2_S1, &
       fac_2d, & ! qingyu, 11/24/2020
       fac1_2d, dfac_2d, poten1_2d, Jwd_2d

  use ModGITM

  implicit none

  ! Local 
  integer :: iLon, iLat, iAlt
  integer :: ig(2,2)

  real :: mlon2, mlat2, pot, ed1_1, ed2_1, fac_1, fac_2, dfac_1, pot1
  real :: val
  real :: wgt(2,2)

  ! ---------------------------------------------------------------------------
  ! Main

  !write(*,*) maxval(MLongitude),minval(MLongitude)

  do iAlt=-1,nAlts+2
     do iLon=-1,nLons+2
        do iLat=-1,nLats+2

           mlon2=MLongitude(iLon,iLat,iAlt,1)
           mlat2=MLatitude(iLon,iLat,iAlt,1)

           !!! Find the positions of the points surrounding the grid point   
           ! Find the position in the longitudinal direction                 
           call findpos_lon(mlon2,mlon_3d,ig(1,:),wgt(1,:),nmlon)
           ! Find the position in the latitudinal direction                  
           call findpos_lat(mlat2,mlat_3d,ig(2,:),wgt(2,:),nmlat_T1)
           
           !!! 2D Interpolation                                              
           ! Potential                                                       
           call interl2d(poten_2d,nmlon,nmlat_T1,ig,wgt,pot)
           ! Ed1                                                             
           call interl2d(ed1_S1,nmlon,nmlat_T1,ig,wgt,ed1_1)
           ! Ed2                                                             
           call interl2d(ed2_S1,nmlon,nmlat_T1,ig,wgt,ed2_1)
           ! FAC
           call interl2d(fac_2d,nmlon,nmlat_T1,ig,wgt,fac_1)
           ! FAC1
           call interl2d(fac1_2d,nmlon,nmlat_T1,ig,wgt,fac_2)
           ! dFAC
           call interl2d(dfac_2d,nmlon,nmlat_T1,ig,wgt,dfac_1)
           ! Potential
           call interl2d(poten1_2d,nmlon,nmlat_T1,ig,wgt,pot1)

           ! Potential
           Gedy_pot(iLon,iLat,iAlt) = pot

           ! Electric field
           Gedy_ed1(iLon,iLat,iAlt) = ed1_1
           Gedy_ed2(iLon,iLat,iAlt) = ed2_1

           Gedy_Efield(iLon,iLat,iAlt,:)=&
                Gedy_ed1(iLon,iLat,iAlt)*b0_d1(iLon,iLat,iAlt,:,1)+&
                Gedy_ed2(iLon,iLat,iAlt)*b0_d2(iLon,iLat,iAlt,:,1)

           ! FAC, qingyu, 11/25/2020
           Gedy_fac(iLon,iLat,iAlt) = fac_1
           Gedy_fac1(iLon,iLat,iAlt) = fac_2  
           Gedy_dfac(iLon,iLat,iAlt) = dfac_1  
           Gedy_pot1(iLon,iLat,iAlt) = pot1

           call interl2d(Jwd_2d,nmlon,nmlat_T1,ig,wgt,val)
           Gedy_Jwd(iLon,iLat,iAlt) = val

        end do
     end do
  end do


end subroutine gedy_get3dvalues_new

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine findpos_lon(mlonin,mlon_ref,pos,weight,nmlon_ref)

  ! Find the relative position of the grid point in the longitude direction  
  ! Reference grid does not loop (-180.~<180.)                               

  ! Arguements                                                               
  integer, intent(in) :: nmlon_ref
  integer, intent(out) :: pos(2)

  real, intent(in) :: mlon_ref(nmlon_ref)
  real, intent(in) :: mlonin
  real, intent(out) :: weight(2)

  ! Local                                                                      
  integer :: i, j, k

  real :: lon0, lon1, frki
  ! ---------------------------------------------------------------------------
  ! Main                                                                     

  if (mlonin>=mlon_ref(nmlon_ref)) then

     pos(1) = nmlon_ref
     pos(2) = 1

     lon0 = mlon_ref(nmlon_ref)
     lon1 = mlon_ref(1)+360.

     frki = (lon1-mlonin)/(lon1-lon0)
     weight(1) = frki
     weight(2) = 1.-frki

  else

     do i=1,nmlon_ref-1
        if (mlonin>=mlon_ref(i) .and. mlonin<mlon_ref(i+1)) then
           pos(1) = i
           pos(2) = i+1

           lon0 = mlon_ref(i)
           lon1 = mlon_ref(i+1)

           frki = (lon1-mlonin)/(lon1-lon0)
           weight(1) = frki
           weight(2) = 1.-frki
           exit
        endif
     enddo

  endif

end subroutine findpos_lon

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine findpos_lat(mlatin,mlat_ref,pos,weight,nmlat_ref)

  ! Find the relative position of the grid point in the longitude direction  
  ! Grid has -90. and +90.                                                  

  ! Arguements                                                              
  integer, intent(in) :: nmlat_ref
  integer, intent(out) :: pos(2)

  real, intent(in) :: mlat_ref(nmlat_ref)
  real, intent(in) :: mlatin
  real, intent(out) :: weight(2)

  ! Local                                                                   
  integer :: i, j, k

  real :: lat0, lat1, frki

  ! ---------------------------------------------------------------------------
  ! Main                                                                    

  do i=1,nmlat_ref-1
     if (mlatin>=mlat_ref(i) .and. mlatin<mlat_ref(i+1)) then
        pos(1) = i
        pos(2) = i+1

        lat0 = mlat_ref(i)
        lat1 = mlat_ref(i+1)

        frki = (lat1-mlatin)/(lat1-lat0)
        weight(1) = frki
        weight(2) = 1.-frki
        exit
     endif
  enddo

end subroutine findpos_lat

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine interl2d(fin,nlonin,nlatin,pos,weight,fout)

  ! 2D interpolation                                                          

  ! Arguements                                                                
  integer, intent(in) :: nlonin, nlatin
  integer, intent(in) :: pos(2,2)

  real, intent(in) :: fin(nlonin,nlatin)
  real, intent(in) :: weight(2,2)
  real, intent(out) :: fout

  ! Local                                                                     
  integer :: nl, nr, nu, nd

  real :: frki(2,2), value(2,2)

  ! ---------------------------------------------------------------------------
  ! Main                                                                     

  nl = pos(1,1)
  nr = pos(1,2)
  nd = pos(2,1)
  nu = pos(2,2)

  value(1,1) = fin(nl,nd) ! LB                                                
  value(1,2) = fin(nr,nd) ! RB                                                
  value(2,1) = fin(nl,nu) ! LT                                                
  value(2,2) = fin(nr,nu) ! RT                                                

  frki(1,1) = weight(1,1)*weight(2,1) ! LB                                    
  frki(1,2) = weight(1,2)*weight(2,1) ! RB                                    
  frki(2,1) = weight(1,1)*weight(2,2) ! LT                                    
  frki(2,2) = weight(1,2)*weight(2,2) ! RT                                    

  fout = 0.

  do i=1,2
     do j=1,2
        fout=fout+value(i,j)*frki(i,j)
     enddo
  enddo

end subroutine interl2d
