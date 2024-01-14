! Calculate the electric potential at GITM's grid 
! Created: Qingyu Zhu, 05/30/3030
!
! -----------------------------------------------------------------------------
subroutine run_epm(iBlock)

  use epm, only: epm_main, &
       epm_coeffs,epm_scale, epm_er, epm_disp, &
       epm_coeffs_sh,epm_scale_sh, epm_disp_sh, &
       epm_LMAX

  use ModGITM

  use ModMpi

  implicit none                                                                
                                                                               
  integer, intent(in) :: iBlock                                                
                                                                               
  integer :: iLon, iLat, iAlt
  real :: mlat1, mlt1, val

  logical :: isNorth = .true.

  !! Added on 08/17/2020
  !! For MPI 
  real :: coeffs1((epm_LMAX+1)**2), coeffs2((epm_LMAX+1)**2)
  real :: scale1, scale2, er, disp

  integer :: iError, itask
  integer, dimension(mpi_status_size) :: status
  integer, parameter :: tag1=321, tag2=432, tag3=543, tag4=654
  integer, parameter :: ncoeff=(epm_LMAX+1)**2, tag5=765, tag6=876

  EPM_Potential=0.


  !! Calculate coeffs, scale, er, disp
  !! Only run it in the root processor
  
  if (iProc==0) then
     
     !write(*,*) "Calculate the coeffs, scale, er in the root proccessor"
     call epm_main

     !write(*,*) "MPI working"
     !! Then other procssors need coff, scale, er and disp
     
     coeffs1(:)=epm_coeffs(:)
     coeffs2(:)=epm_coeffs_sh(:)
     scale1=epm_scale
     scale2=epm_scale_sh

     if (scale1>4.) then
        write(*,*) "------ NH scaling factor goes wrong, Not scaling ------"
        scale1=1.
     end if

     if (scale1>4.) then
        write(*,*) "------ SH scaling factor goes wrong, Not scaling ------"
        scale2=1.
     end if

     er=epm_er

     ! Set disp to zero
     epm_disp=0.
     epm_disp_sh=0.
     disp=epm_disp

     do itask=1,nProcs-1
        call mpi_send(coeffs1,ncoeff,mpi_real,itask,tag1,iCommGITM,iError)
        call mpi_send(coeffs2,ncoeff,mpi_real,itask,tag2,iCommGITM,iError)
        call mpi_send(scale1,1,mpi_real,itask,tag3,iCommGITM,iError)
        call mpi_send(scale2,1,mpi_real,itask,tag4,iCommGITM,iError)
        call mpi_send(er,1,mpi_real,itask,tag5,iCommGITM,iError)           
        call mpi_send(disp,1,mpi_real,itask,tag6,iCommGITM,iError)
     end do

     !write(*,*) iProc, "---------- Send ----------"
     !write(*,*) iProc, coeffs1(1:10), coeffs2(1:10), scale1, scale2, er, disp

  else
     
     call mpi_recv(coeffs1,ncoeff,mpi_real,0,tag1,iCommGITM,status,iError)
     call mpi_recv(coeffs2,ncoeff,mpi_real,0,tag2,iCommGITM,status,iError)
     call mpi_recv(scale1,1,mpi_real,0,tag3,iCommGITM,status,iError)
     call mpi_recv(scale2,1,mpi_real,0,tag4,iCommGITM,status,iError)
     call mpi_recv(er,1,mpi_real,0,tag5,iCommGITM,status,iError)
     call mpi_recv(disp,1,mpi_real,0,tag6,iCommGITM,status,iError)
     
     !write(*,*) iProc, "---------- Received ----------" 
     !write(*,*) iProc, coeffs1(1:10), coeffs2(1:10), scale1, scale2, er, disp
     
     epm_coeffs(:)=coeffs1(:)
     epm_coeffs_sh(:)=coeffs2(:)
     epm_scale=scale1
     epm_scale_sh=scale2
     epm_er=er
     epm_disp=disp
     epm_disp_sh=-disp

  end if

  call mpi_barrier(iCommGITM,iError)

  !! Calculate potential according to coeffs, scale, er and disp
  do iLat = -1,nLats+2                                                         
     do iLon = -1,nLons+2                                                      
        do iAlt = -1,nAlts+2                                                   
                                                                               
           ! Get GITM grids                                                    
           mlat1=(MLatitude(iLon, iLat, iAlt, iBlock))                         
           mlt1=MLT(iLon, iLat, iAlt)                                          
                                                                               
           if (mlat1<0) then 
              isNorth = .false.
           else
              isNorth = .true.
           end if

           ! Adjust the coeffs 
           mlat1=abs(mlat1) 
           if (mlt1<0.) mlt1=mlt1+24.                                          
           if (mlt1>24.) mlt1=mlt1-24.

           if (isNorth) then 
              call calc_epm_sp_epot(mlt1,mlat1,epm_coeffs,epm_scale,&
                   epm_er,epm_disp,val)
           else
              call calc_epm_sp_epot(mlt1,mlat1,epm_coeffs_sh,epm_scale_sh,&
                   epm_er,epm_disp_sh,val)
           end if

           ! kV --> V                                                          
           EPM_Potential(iLon,iLat,iAlt) = val * 1000. * 1.07

        end do
     end do
  end do

end subroutine run_epm

! Similar as UA_gedy_get2Dpoten used for Gedy
subroutine gedy_EPM(nmlt,nmlat,mltin,mlatin,potout)

  use epm, only: epm_main, &
       epm_coeffs,epm_scale, epm_er, epm_disp, &
       epm_coeffs_sh,epm_scale_sh, epm_disp_sh, &
       epm_LMAX

  implicit none 

  integer, intent(in) :: nmlt, nmlat                                           
  real, intent(in) :: mltin(nmlt), mlatin(nmlat)                               
  real, intent(out):: potout(nmlt,nmlat)                                       
                                                                               
  integer :: imlt, imlat
  real :: mlat1, mlt1, val
  logical :: isNorth

  !! Calculate potential according to coeffs, scale, er and disp
  do imlt = 1,nmlt                                                         
     do imlat = 1,nmlat                                                      

        mlt1=mltin(imlt) 
        mlat1=mlatin(imlat) 
                                                                               
        if (mlat1<0) then 
           isNorth = .false.
        else
           isNorth = .true.
        end if

        ! Adjust the coeffs 
        mlat1=abs(mlat1) 
        if (mlt1<0.) mlt1=mlt1+24.                                          
        if (mlt1>24.) mlt1=mlt1-24.

        if (isNorth) then 
           call calc_epm_sp_epot(mlt1,mlat1,epm_coeffs,epm_scale,&
                epm_er,epm_disp,val)
        else
           call calc_epm_sp_epot(mlt1,mlat1,epm_coeffs_sh,epm_scale_sh,&
                epm_er,epm_disp_sh,val)
        end if

        ! kV --> V                                                          
           potout(imlt,imlat) = val * 1000. * 1.07

        end do
     end do

end subroutine gedy_EPM
