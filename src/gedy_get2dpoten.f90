! This routine will get the 2D potential distribution in the geomagnetic 
! coordinate
!
! Created: Qingyu Zhu, 07/28/2017, qingyu.zhu@mavs.uta.edu
!
! Set up a MLON-MLAT grid, then get the MLT for each MLON, finally obtain the 
! potential
!
! Why we need this routine?
! To provide the high-latitude forcing
! Be very careful about this
! Add in this GITM on 03/02/2020, qingyu  
! -----------------------------------------------------------------------------

subroutine gedy_get2Dpoten

  ! Global
  use ModGITM, only: SubsolarLatitude, SubsolarLongitude, &
       MagneticPoleColat, MagneticPoleLon

  use ModGedy, only: nmlon_gitm, nmlat_gitm, &
       mlon1_gitm, mlt1_gitm, mlat1_gitm, poten2d_gitm

  use ModInputs, only: useEPMpotential, useAMIEpotential

  implicit none
  
  ! Local
  integer :: imlon, imlat, iError 
  
  real :: dmlon, dmlat

  character(len=*),parameter :: &
       outdir='/home1/03751/qyzhuta/gitm_all/gitm_edy/myrun/data/'
  ! ---------------------------------------------------------------------------

  ! Allocate arrays

  if (.not. allocated(mlon1_gitm)) then
     allocate(mlon1_gitm(nmlon_gitm),stat=iError)
     allocate(mlat1_gitm(nmlat_gitm),stat=iError)
     allocate(mlt1_gitm(nmlon_gitm),stat=iError)
     allocate(poten2d_gitm(nmlon_gitm,nmlat_gitm),stat=iError)
  endif

  ! Set up the MLON-MLAT grid
  
  dmlon = 360./nmlon_gitm
  dmlat = 180./nmlat_gitm

  do imlon=1,nmlon_gitm
     mlon1_gitm(imlon)=dmlon/2.+(imlon-1)*dmlon-180.

     call magloctm( &
          mlon1_gitm(imlon), &
          SubsolarLatitude,   &
          SubsolarLongitude,  &
          MagneticPoleColat, &
          MagneticPoleLon,   &
          mlt1_gitm(imlon))
     if (mlt1_gitm(imlon) < 0) &
          mlt1_gitm(imlon) = mlt1_gitm(imlon) + 24.0
  enddo

  do imlat=1,nmlat_gitm
     mlat1_gitm(imlat)=dmlat/2.-90.+(imlat-1)*dmlat
  enddo

  ! Get the potential for each point
  ! Currently, this potential comes from Weimer 05
  ! Add capability of using EPM potential, qingyu 01/26/2021
  ! AMIE capability added on 
  if (useEPMpotential) then
     call gedy_EPM(nmlon_gitm,nmlat_gitm,mlt1_gitm, &                    
          mlat1_gitm,poten2d_gitm)
  elseif (useAMIEpotential) then
     call gedy_amie_pot(nmlon_gitm,nmlat_gitm,mlt1_gitm, &
          mlat1_gitm,poten2d_gitm)
  else
     call UA_gedy_get2Dpoten(nmlon_gitm,nmlat_gitm,mlt1_gitm, &
          mlat1_gitm,poten2d_gitm)
  end if

  ! --------------------------------------------------------------------------
  ! Check
  !open(10,file=outdir//'gitm_pot_mag',status='replace')
  
  !do imlon=1,nmlon_gitm
  !   write(10,5) poten2d_gitm(imlon,:)
  !enddo

!5 format(72F10.2)
  !close(10)

end subroutine gedy_get2Dpoten
