! Routines to calculate the height-integrated parameters 
! Created: Qingyu Zhu, 07/08/2020

subroutine calc_height_integral(parms_in,alts_in,nalts_in,fac,integral)

  implicit none

  integer, intent(in) :: nalts_in
  real, intent(in) :: parms_in(nalts_in), alts_in(nalts_in), fac
  real, intent(out) :: integral 

  real :: tmp(nalts_in-1), dalt(nalts_in-1)
  integer :: i, j, k

  integral=0.

  tmp(:)=(parms_in(2:nalts_in)+parms_in(1:nalts_in-1))/2.
  dalt(:)=alts_in(2:nalts_in)-alts_in(1:nalts_in-1)

  ! Set the upper boundary at 600 km (60000 m)
  j=minloc(abs(alts_in-500000.),dim=1)
  !write(*,*) "Upper boundary:", j, alts_in(j)

  integral=dot_product(tmp(1:j),dalt(1:j))

  integral=integral*fac


end subroutine calc_height_integral

! Calculate the TEC
! -----------------------------------------------------------------------------
subroutine calc_tec(iBlock)

  use ModGITM, only: nLats, nLons, nAlts, scTEC, IDensityS, Altitude_GB
  use ModPlanet, only : ie_

  implicit none 

  integer, intent(in) :: iBlock 
  
  integer :: iLat, iLon
  real :: alts_in(nAlts), parms_in(nAlts), integral
  real, parameter :: tec_fac=1e-16

  scTEC=0.
  
  do iLat=1,nLats
     do iLon=1,nLons

        alts_in = Altitude_GB(iLon,iLat,1:nAlts,iBlock)
        parms_in = IDensityS(iLon,iLat,1:nAlts,ie_,iBlock)


        call calc_height_integral(parms_in,alts_in,nAlts,tec_fac,integral)

        scTEC(iLon,iLat) = integral

     end do
  end do

end subroutine calc_tec

! Calculate the heigh-integrated Joule heating
! -----------------------------------------------------------------------------
subroutine calc_integrated_jh(iBlock)

  use ModGITM, only: nLats, nLons, nAlts, Altitude_GB, &
       cp, Rho, TempUnit, HeightIntegratedJH
  use ModSources, only: JouleHeating

  implicit none 

  integer, intent(in) :: iBlock 
  
  integer :: iLat, iLon, iAlt
  real :: alts_in(nAlts), parms_in(nAlts), integral

  HeightIntegratedJH=0.
  
  do iLat=1,nLats
     do iLon=1,nLons

        alts_in = Altitude_GB(iLon,iLat,1:nAlts,iBlock)

        ! Joule heating needs to be coverted 
        do iAlt=1,nAlts
           parms_in(iAlt) = JouleHeating(iLon,iLat,iAlt)* &
                TempUnit(iLon,iLat,iAlt)*cp(iLon,iLat,iAlt,iBlock)*&
                Rho(iLon,iLat,iAlt,iBlock)
        end do

        call calc_height_integral(parms_in,alts_in,nAlts,1.,integral)

        HeightIntegratedJH(iLon,iLat)=integral

     end do
  end do

end subroutine calc_integrated_jh

! New routine to calculate the height-integrated parameters 
! -----------------------------------------------------------------------------
subroutine calc_height_integral1(parms_in,alts_in,nalts_in,&
     alt1, alt2, fac,integral)

  implicit none

  integer, intent(in) :: nalts_in
  real, intent(in) :: parms_in(nalts_in), alts_in(nalts_in), alt1, alt2, fac
  real, intent(out) :: integral 

  real :: tmp(nalts_in-1), dalt(nalts_in-1)
  integer :: i, j, k

  integral=0.

  tmp(:)=(parms_in(2:nalts_in)+parms_in(1:nalts_in-1))/2.
  dalt(:)=alts_in(2:nalts_in)-alts_in(1:nalts_in-1)

  ! Set the lower boundary (alt1*1000.)
  j=minloc(abs(alts_in-alt1*1000.),dim=1)
  !write(*,*) "Upper boundary:", j, alts_in(j)

  ! Set the upper boundary (alt2*1000.)
  k=minloc(abs(alts_in-alt2*1000.),dim=1)
  if (k>(nalts_in-1)) k=nalts_in-1

  integral=dot_product(tmp(j:k),dalt(j:k))

  integral=integral*fac

end subroutine calc_height_integral1

! calculate O/N2 ratio
! -----------------------------------------------------------------------------
subroutine calc_on2ratio(iBlock)

  use ModGITM, only: nLats, nLons, nAlts, on2ratio, NDensityS, Altitude_GB
  use ModPlanet, only : iO_3P_, iN2_

  implicit none 

  integer, intent(in) :: iBlock 
  
  integer :: iLat, iLon
  real :: alts_in(nAlts), parms_in(nAlts), integral
  real :: sumo, sumn2

  on2ratio=0.
  
  do iLat=1,nLats
     do iLon=1,nLons

        alts_in = Altitude_GB(iLon,iLat,1:nAlts,iBlock)

        ! GUVI O/N2 ratio is the ratio of column density of O
        ! to column density of N2 above the altitude where 
        ! the column density of N2 is 10^17 cm-2
        ! The altitude is roughly at 140 km
        ! TIMED altitude is at 625 km

        parms_in = NDensityS(iLon,iLat,1:nAlts,iO_3P_,iBlock)
        call calc_height_integral1(parms_in,alts_in,nAlts,&
             140.,625.,1.,sumo)

        parms_in = NDensityS(iLon,iLat,1:nAlts,iN2_,iBlock)
        call calc_height_integral1(parms_in,alts_in,nAlts,&
             140.,625.,1.,sumn2)

        on2ratio(iLon,iLat) = sumo/sumn2

     end do
  end do

end subroutine calc_on2ratio
