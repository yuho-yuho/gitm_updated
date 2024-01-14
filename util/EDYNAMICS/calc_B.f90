!---------------------------------------------------------------------- 
    module delB_module
! 
      use params_module, only: nlat_qd,nhgt_fix,nmlon,nhgt_fix_r, &
         rtd,dtr,re,hgt_fix_r,hgt_fix,pi,ylonm,nglon,nglat,glon,glat
      use qd_module, only: qd,I3qd,lat_qd_mp,lon_qd_mp,Jf1hor,Jf2hor 
! 
      implicit none
! 
! mmax,nmax are maximum order and degree of spherical harmonics used
!  in geographic coordinates.
      integer, parameter :: mmax = 40 ,nmax = 120
!      integer, parameter :: mmax = 12 ,nmax = 24
      real :: pnmost(0:nmax,0:mmax),dpnm(0:nmax,0:mmax),fm(-mmax:mmax)
      real :: PMOPMMO(mmax+1), R(0:nmax,0:mmax), SQ2,sq4pi
! SH coefficients for toroidal magnetic field (current layer base)
!  (divided by radius):
      real :: tor(nmax,-mmax:mmax,nhgt_fix_r)
! SH coefficients for magnetic potential (current layer base)
!  (divided by radius):
      real :: vor(nmax,-mmax:mmax,nhgt_fix_r)
! SH coefficients for upward field within ionosphere (current layer base):
      real :: beta(nmax,-mmax:mmax,nhgt_fix_r)
! SH coefficients for external ground magnetic potential
!  (divided by radius):
      real :: vextor(nmax,-mmax:mmax)
! SH coefficients for external upward field at ground:
      real :: betaext(nmax,-mmax:mmax)      
! SH coefficients for equivalent current function:
      real :: psicoef(nmax,-mmax:mmax)
      
      real, parameter :: h_leo = 3.0e5   ! height for dB calculation in [m]
      real :: bhat_g(nglon,nglat,3), &   ! unit vector (east, north, up) along magnetc fieldline
           qdlon_g(nglon,nglat),qdlat_g(nglon,nglat)
      
      real,dimension(nglon,nglat) :: delbegrd, delbngrd ,delbugrd,psi
      real,dimension(nglon,nglat) :: delbe400,delbn400,delbu400,delbscl
      real,dimension(nglat,nhgt_fix_r) :: delbe70W,delbn70W,delbu70W
      
      real,dimension(nmlon,nlat_qd-1) :: delbegrd_qd, delbngrd_qd ,delbugrd_qd,psi_qd
      real,dimension(nmlon,nlat_qd-1) :: delbe400_qd,delbn400_qd,delbu400_qd,delbscl_qd
      real,dimension(nlat_qd-1,nhgt_fix) :: delbe0ln_qd,delbn0ln_qd,delbu0ln_qd,delscllH_qd

     contains
!---------------------------------------------------------------------- 
      subroutine gen_geo_grid
!     
      implicit none
 !     
      integer :: ilateq,i
      real :: dlon,dlat
!      
! set up geographic grid
      dlon = 360./float(nglon-1)
      dlat = 180./float(nglat-1)
      ilateq = (nglat+1)/2
      do i=1,nglon
        glon(i) = -180.+(i-1)*dlon
      enddo
      do i=1,nglat
        glat(i) = (i-ilateq)*dlat
      enddo     
!      
      end subroutine gen_geo_grid
!---------------------------------------------------------------------- 
      subroutine calc_Bcoef
!     
      implicit none
!     
      integer, parameter :: mmaxq = 6 ,nmaxq = 6
      integer :: i,j,k,kk,l,ll,isn,n,m,nm,ma
      real :: fac,fac2,dlatq,dlonq,jeast,jsouth,nnp1,eta &
        ,QCNST1,QCNST2,X,XL,mu0,bradius,dr,glatmid,glonmid,rob,robnm1,bornp2
! SH coefficients for toroidal or equivalent currents
      real :: C1(nmax,-mmax:mmax),C2(nmaxq,-mmaxq:mmaxq)
! SH coefficients for radial current density times r**2
      real :: bnm(nmax,-mmax:mmax)
! q = Qnm of Richmond (1974) when m is non-negative.  That paper has a
!  sign error when m is negative, when Qnm should reverse sign.
      real :: q(nmaxq,-mmaxq:mmaxq)
      real :: CT,ST,STS,CP,SP
!
      SQ2 = sqrt(2.e0)
      sq4pi = sqrt(4.*pi)
      do m = 0,mmax 
        if (m.ne.0) PMOPMMO(m) = sqrt(1. + .5/float(m))
        do n=m,nmax
          R(n,m) = sqrt(float(n*n-m*m)/amax1(4.*n*n-1.,1.))
        enddo
      enddo
!
      dlatq =  abs(lat_qd_mp(2)-lat_qd_mp(1))  ! regular grid spacing
      dlonq =  abs(lon_qd_mp(2)-lon_qd_mp(1))  ! regular grid spacing
! Find coefficients C2 for equivalent currents at hgt_fix_r(nhgt_fix_r)
!  that represent magnetic effects of FAC above that level.
      mu0 = 4.E-7*pi
      QCNST1 = 3.*sqrt(10.)/16.E0
      QCNST2 = 3.*sqrt(70.)/32.E0
      C2 = 0.
      do l=1,nlat_qd-1  ! loop over all qd latitude midpoints
! Select hemisphere and set isn and pole-to-equator index ll accordingly.
        isn = 1
        ll = l
        if (lat_qd_mp(l).gt.0.) then
          isn = 2
          ll = nlat_qd - l
        endif 
! Calculate q functions for this magnetic latitude.
        q = 0.
        CT = sin(lat_qd_mp(l))
        ST = cos(lat_qd_mp(l))
        STS = ST*ST               
! Following formula for XL does not work at magnetic equator (CT = 0.),
!  but lat_qd_mp does not include the equator.
        XL = 1. + .5*STS*alog((1. + CT)/(1. - CT))/CT
        q(1,1) = 1.2247449*CT*ST
        q(2,2) = 1.3693064*CT*STS*XL
        q(3,1) =  .2291288*CT*ST*(4. - STS*(3. + 6.*STS))
        q(3,3) = 1.4790199*CT*STS*ST*(1. + 2.*STS)      
        X = STS*(2. + 3.*STS*XL)                       
        q(4,2) = QCNST1*CT*STS*(4. - X)             
        q(4,4) = QCNST2*CT*STS*X                   
        X = STS*(3. + STS*(4. + 8.*STS))            
        q(5,1) = .018021728*CT*ST*(56. - STS*(188. - 13.*X))
        q(5,3) =  .5255737*CT*STS*ST*(8. - X)              
        q(5,5) =  .5484353*CT*STS*ST*X                    
        X = STS*(8. + STS*(10. + 15.*STS*XL))            
        q(6,2) = .057727979*CT*STS*(64. - STS*(168. - 3.*X))
        q(6,4) = .079047290*CT*STS*STS*(80.  - 3.*X)       
        q(6,6) =  .2140611*CT*STS*STS*X                   
        do n=1,nmaxq
          do m=-mmaxq,-1
            q(n,m) = -q(n,-m)
          enddo ! m
        enddo ! n
        do i=1,nmlon
          CP = cos(ylonm(i))
          SP = sin(ylonm(i))
          fm(0) = 1.
          fm(-1) = SQ2*CP
          fm( 1) = SQ2*SP
          do m=2,mmaxq
            fm(-m) = CP*fm(1-m) - SP*fm(m-1)
            fm( m) = CP*fm(m-1) + SP*fm(1-m)
          enddo ! m
          do n=1,nmaxq
            nnp1 = float(n*(n+1))
            fac = 1.e0/(nnp1*sq4pi)
            do m=-mmaxq,mmaxq
              C2(n,m) = C2(n,m) - fac*I3qd(i,ll,nhgt_fix_r,isn)*q(n,m)*fm(-m)
            enddo ! m
          enddo ! n
        enddo ! i
      enddo ! l
      C1 = 0.
      tor = 0.
      vor = 0.
      beta = 0.
      vextor = 0.
      betaext = 0.
! Calculate C1 coefficients corresponding to FAC above hgt_fix_r.
      do l=1,nlat_qd-1
! fac is a solid-angle element of the QD grid, times F, divided by sq4pi
        fac = dlatq*dlonq*cos(lat_qd_mp(l))/sq4pi
        do i=1,nmlon
! Calculate equivalent current function from C2 coefficients, which are
!  applicable to QD coordinates.
          eta = 0.
          call sphhar(lat_qd_mp(l)*rtd,ylonm(i)*rtd,ST)
          fac2 = ST/sq4pi
          do n=1,nmaxq
            do m=-mmaxq,mmaxq
              ma = iabs(m)
! Note: pnmost*fac2 is Pnm/sq4pi
              eta = eta + C2(n,m)*pnmost(n,ma)*fm(m)*fac2
            enddo ! m
          enddo ! n
! Perform integrals over Earth at hgt_fix_r for each geographic
!  spherical harmonic component.
          call sphhar(qd(i,l)%glat(nhgt_fix_r),qd(i,l)%glon(nhgt_fix_r),ST)
          fac2 = fac*ST/qd(i,l)%F3(nhgt_fix_r)
! Note: ST is included in fac2 to convert pnmost to Pnm.
          do n=1,nmax
            nm = min0(mmax,n)
            do m=-nm,nm
              ma = iabs(m)
! Note: pnmost*fac2 is dlatq*dlonq*cos(lat_qd_mp(l)*rtd)*Pnm/(sq4pi*F)
              C1(n,m) = C1(n,m) + fac2*eta*pnmost(n,ma)*fm(m)
            enddo ! m
          enddo ! n
        enddo ! i
      enddo ! l
      bradius = re + hgt_fix_r(nhgt_fix_r)
! Contribution of FAC to external coefficients (divided by r) and
!  upward field coefficients at base of each current layer.
      do kk=1,nhgt_fix_r
        rob = (re + hgt_fix_r(kk))/bradius
        robnm1 = 1./rob
        do n=1,nmax
          robnm1 = robnm1*rob
          fac2 = mu0*robnm1*float(n+1)/(float(2*n+1)*bradius)
          nm = min0(mmax,n)
          do m=-nm,nm
            vor(n,m,kk) = vor(n,m,kk) + fac2*C1(n,m)
            beta(n,m,kk) = beta(n,m,kk) - n*fac2*C1(n,m)
          enddo ! m
        enddo ! n
      enddo ! kk
! Contribution of FAC to external coefficients (divided by Re) and
!  upward field coefficients at Earth's surface
      rob = re/bradius
      robnm1 = 1./rob
      do n=1,nmax
        robnm1 = robnm1*rob
        fac2 = mu0*robnm1*float(n+1)/(float(2*n+1)*bradius)
        nm = min0(mmax,n)
        do m=-nm,nm
          vextor(n,m) = vextor(n,m) + fac2*C1(n,m)
          betaext(n,m) = betaext(n,m) - n*fac2*C1(n,m)
        enddo ! m
      enddo ! n
!
! Calculate contributions due to ionospheric currents.
      do k=1,nhgt_fix_r 
        bnm = 0.
        C1 = 0.
        if (k.le.nhgt_fix) dr = hgt_fix_r(k+1) - hgt_fix_r(k)
        do l=1,nlat_qd-1  ! loop over all qd latitude midpoints
! Select hemisphere and set isn and pole-to-equator index ll accordingly.
          isn = 1
          ll = l
          if (lat_qd_mp(l).gt.0.) then
            isn = 2
            ll = nlat_qd - l
          endif 
          if (k.le.nhgt_fix) &
            fac = (re + hgt_fix(k))*dlatq*dlonq*cos(lat_qd_mp(l))*dr/sq4pi
          do i=1,nmlon
! Exclude k=nhgt_fix_r (top of top current layer)
            if (k.le.nhgt_fix) then
! Find coefficients C1 for toroidal currents
              if (qd(i,l)%Fqd(k).eq.0..or.qd(i,l)%Fqd(k).eq.-0.) then
                write(6,*) 'qd(i,l)%Fqd(k)=0. for i,l,k= ',i,l,k
                stop
              endif
! jeast,jsouth are geographic current densities 
              jeast =   Jf1hor(i,ll,k,isn)*qd(i,l)%f11(k) &
                      + Jf2hor(i,ll,k,isn)*qd(i,l)%f21(k) 
              jsouth = -Jf1hor(i,ll,k,isn)*qd(i,l)%f12(k) &
                      - Jf2hor(i,ll,k,isn)*qd(i,l)%f22(k) 
! Perform integrals over Earth for each spherical harmonic component
              glatmid = .5*(qd(i,l)%glat(k) + qd(i,l)%glat(k+1))
              glonmid = .5*(qd(i,l)%glon(k) + qd(i,l)%glon(k+1))
              call sphhar(glatmid,glonmid,ST)
              fac2 = fac/qd(i,l)%Fqd(k)
              do n=1,nmax
                nnp1 = float(n*(n+1))
                nm = min0(mmax,n)
                do m=-nm,nm
                  ma = iabs(m)
                  C1(n,m) = C1(n,m)+ fac2*(dpnm(n,ma)*fm(m)*jeast &
                    - m*pnmost(n,ma)*fm(-m)*jsouth)/nnp1
                enddo ! m
              enddo ! n
            endif ! excluding k=nhgt_fix_r (above top current layer)
! Find coefficients bnm for radial current density times r**2, at base of layer
            call sphhar(qd(i,l)%glat(k),qd(i,l)%glon(k),ST)
            fac2 = I3qd(i,ll,k,isn)*ST/sq4pi
! Note: ST is included in fac2 to convert pnmost to Pnm.
            do n=1,nmax
              nm = min0(mmax,n)
              do m=-nm,nm
                ma = iabs(m)
                bnm(n,m) = bnm(n,m) + fac2*pnmost(n,ma)*fm(m)
              enddo ! m
            enddo ! n
	  end do !  i
 	end do !  l
! Find coefficients tor for toroidal magnetic fields
        do n=1,nmax
          fac2 = mu0/(float(n*(n+1))*(re+hgt_fix_r(k)))
          nm = min0(mmax,n)
          do m=-nm,nm
            tor(n,m,k) = fac2*bnm(n,m)
          enddo ! m
        enddo ! n
! Calculate contributions of currents in this layer to coefficients vor and beta
!  for poloidal magnetic fields at base of each layer (divided by r).
       if (k.le.nhgt_fix) then ! only for current layers
        bradius = re + hgt_fix(k)
        do kk=1,k ! coefficients below current layer
          rob = (re + hgt_fix_r(kk))/bradius
          robnm1 = 1./rob
          do n=1,nmax
            robnm1 = robnm1*rob
            fac2 = mu0*robnm1*float(n+1)/(float(2*n+1)*bradius)
            nm = min0(mmax,n)
            do m=-nm,nm
              vor(n,m,kk) = vor(n,m,kk) + fac2*C1(n,m)
              beta(n,m,kk) = beta(n,m,kk) - n*fac2*C1(n,m)
            enddo ! m
          enddo ! n
        enddo ! kk
        do kk=k+1,nhgt_fix_r ! coefficients above current layer
          rob = (re + hgt_fix_r(kk))/bradius
          bornp2 = 1./(rob*rob)
          do n=1,nmax
            bornp2 = bornp2/rob
            fac2 = mu0*bornp2*float(n)/(float(2*n+1)*bradius)
            nm = min0(mmax,n)
            do m=-nm,nm
              vor(n,m,kk) = vor(n,m,kk) - fac2*C1(n,m)
              beta(n,m,kk) = beta(n,m,kk) - float(n+1)*fac2*C1(n,m)
            enddo ! m
          enddo ! n
        enddo ! kk
! External coefficients (divided by Re) and upward field at Earth's surface
        rob = re/bradius
        robnm1 = 1./rob
        do n=1,nmax
          robnm1 = robnm1*rob
          fac2 = mu0*robnm1*float(n+1)/(float(2*n+1)*bradius)
          nm = min0(mmax,n)
          do m=-nm,nm
            vextor(n,m) = vextor(n,m) + fac2*C1(n,m)
            betaext(n,m) = betaext(n,m) - n*fac2*C1(n,m)
          enddo ! m
        enddo ! n
       endif ! only for current layers
!         
      enddo ! k
!      
! Calculate equivalent current coefficients for 110 km current layer.
      bradius = re + 110.e3
      rob = re/bradius
      robnm1 = 1./rob
      do n=1,nmax
        robnm1 = robnm1*rob
        fac2 = mu0*robnm1*float(n+1)/(float(2*n+1)*bradius)
        nm = min0(mmax,n)
        do m=-nm,nm
          psicoef(n,m) =  vextor(n,m)/fac2
        enddo ! m
      enddo ! n
!      
      end subroutine calc_Bcoef
!---------------------------------------------------------------------- 
      subroutine calc_B
! Example that calculates the internal field using perfect conductor at
!  depth=600 km, and that calculates total perturbation fields over the
!  Earth's surface and at 400 km altitude, and in a meridional
!  latitude-height slice at 70W longitude at the heights in the array
!  hgt_fix_r.
! nglon=73, nglat=91 gives 5 x 2 deg lon-lat grid for magnetic perturbations
!      implicit none
!      integer, parameter :: nglon = 73, nglat = 91 now defined in params.f90
!      integer, parameter :: nglon = 37, nglat = 19
      real, parameter :: depth = 6.e5  ! meters
      real :: vintor(nmax,-mmax:mmax),betaint(nmax,-mmax:mmax)
!      
      integer :: n,m,nm,ma,lon,lat,lateq,k,kbelow
      real :: dlon,dlat,glon_loc,glat_loc,frac,vorsum,betasum,torinterp
      real :: fac,coa,coa2,coa2np1,aob,aobnp1,ST
!      real :: h_leo  ! height of magn. perturb. calculation
       
! Find internal coefficients
      coa = (re-depth)/re
      coa2 = coa**2
      coa2np1 = coa
      do n=1,nmax
        coa2np1 = coa2np1*coa2
        fac = coa2np1*float(n)/float(n+1)
        nm = min0(mmax,n)
        do m=-nm,nm
          vintor(n,m) = fac*vextor(n,m)
          betaint(n,m) = -coa2np1*betaext(n,m)
        enddo ! m
      enddo ! n
!
! Find k index for nearest coefficients below 400 km (4x10^5 m)
      do k=1,nhgt_fix_r
        if (hgt_fix_r(k).ge.h_leo) then
          kbelow = k - 1
          frac = (hgt_fix_r(k) - h_leo)/(hgt_fix_r(k)-hgt_fix_r(kbelow))
          goto 100
        endif
      enddo ! k
  100 continue
      if (frac.lt.0.or.frac.gt.1.) then
        write (6,*) 'Stopped because frac =', frac
        stop
      endif
!
! Calculate global fields at ground and 400 km
      dlon = 360./float(nglon-1)
      dlat = 180./float(nglat-1)
      lateq = (nglat+1)/2
      delbegrd = 0.
      delbngrd = 0.
      delbugrd = 0.
      delbe400 = 0.
      delbn400 = 0.
      delbu400 = 0.
      psi      = 0.
!      
      aob = re/(re+h_leo)
      do lon=1,nglon
        ! glon = (lon-1)*dlon
        do lat=1,nglat
          ! glat = (lat-lateq)*dlat
          call sphhar(glat(lat),glon(lon),ST)
          fac = ST/sq4pi
          aobnp1 = aob
          do n=1,nmax
            nm = min0(mmax,n)
            aobnp1 = aobnp1*aob
            do m=-nm,nm
              ma = iabs(m)
              vorsum = vextor(n,m) + vintor(n,m) 
              betasum = betaext(n,m) + betaint(n,m) 
              delbegrd(lon,lat) = delbegrd(lon,lat) &
                - vorsum*m*pnmost(n,ma)*fm(-m)/sq4pi
              delbngrd(lon,lat) = delbngrd(lon,lat) &
                + vorsum*dpnm(n,ma)*fm(m)/sq4pi
              delbugrd(lon,lat) = delbugrd(lon,lat) &
                + betasum*fac*pnmost(n,ma)*fm(m)
              psi(lon,lat) = psi(lon,lat) &
                + psicoef(n,m)*fac*pnmost(n,ma)*fm(m)
              vorsum = vintor(n,m)*aobnp1 + frac*vor(n,m,kbelow) &
                + (1.-frac)*vor(n,m,kbelow+1)
              betasum = betaint(n,m)*aobnp1*aob + frac*beta(n,m,kbelow) &
                + (1.-frac)*beta(n,m,kbelow+1)
              torinterp = frac*tor(n,m,kbelow) + (1.-frac)*tor(n,m,kbelow+1)
              delbe400(lon,lat) = delbe400(lon,lat) &
                - vorsum*m*pnmost(n,ma)*fm(-m)/sq4pi &
                - torinterp*dpnm(n,ma)*fm(m)/sq4pi
              delbn400(lon,lat) = delbn400(lon,lat) &
                + vorsum*dpnm(n,ma)*fm(m)/sq4pi  &
                - torinterp*pnmost(n,ma)*m*fm(-m)/sq4pi 
              delbu400(lon,lat) = delbu400(lon,lat) &
                + betasum*fac*pnmost(n,ma)*fm(m)
            enddo ! m
          enddo ! n
	  delbscl(lon,lat) = delbe400(lon,lat)*bhat_g(lon,lat,1) + &
	    delbn400(lon,lat)*bhat_g(lon,lat,2) + &
	    delbu400(lon,lat)*bhat_g(lon,lat,3)
        enddo ! lat
      enddo ! lon
!
! Calculate fields at 70W
      delbe70W = 0.
      delbn70W = 0.
      delbu70W = 0.
      glon_loc = 30.
      !
      do lat=1,nglat
        ! glat = (lat-lateq)*dlat
        call sphhar(glat(lat),glon_loc,ST)
        fac = ST/sq4pi
        do k=1,nhgt_fix_r
          aob = re/(re+hgt_fix_r(k)) 
          aobnp1 = aob
          do n=1,nmax
            nm = min0(mmax,n)
            aobnp1 = aobnp1*aob
            do m=-nm,nm
              ma = iabs(m)
              vorsum = vintor(n,m)*aobnp1 + vor(n,m,k)
              betasum = betaint(n,m)*aobnp1*aob + beta(n,m,k)
              delbe70W(lat,k) = delbe70W(lat,k) &
                - vorsum*m*pnmost(n,ma)/sq4pi &
                - tor(n,m,k)*dpnm(n,ma)*fm(m)/sq4pi
              delbn70W(lat,k) = delbn70W(lat,k) &
                + vorsum*dpnm(n,ma)*fm(m)/sq4pi &
                - tor(n,m,k)*pnmost(n,ma)*m*fm(-m)/sq4pi
              delbu70W(lat,k) = delbu70W(lat,k) &
                + betasum*fac*pnmost(n,ma)*fm(m)
            enddo ! m
          enddo ! n
        enddo ! k
      enddo ! lat
!      
      end subroutine calc_B
!---------------------------------------------------------------------- 
      subroutine calc_Bqd
! Example that calculates the internal field using perfect conductor at
!  depth=600 km, and that calculates total perturbation fields over the
!  Earth's surface and at 400 km altitude on the QD grid, and in a meridional
!  latitude-height slice at 0 QD longitude at the heights in the array
!  hgt_fix_r.
      implicit none
      real, parameter :: depth = 6.e5  ! meters
      real :: vintor(nmax,-mmax:mmax),betaint(nmax,-mmax:mmax)
      integer :: i,l,n,m,nm,ma,lon,lat,lateq,k,kbelow
      real ::           glon,glat,frac,vorsum,betasum,torinterp
      real :: fac,coa,coa2,coa2np1,aob,aobnp1,ST
      real :: f11interp,f12interp,f21interp,f22interp,Fqdinterp
      real :: lonqd2get, mmin,diff
       
! Find internal coefficients
      coa = (re-depth)/re
      coa2 = coa**2
      coa2np1 = coa
      do n=1,nmax
        coa2np1 = coa2np1*coa2
        fac = coa2np1*float(n)/float(n+1)
        nm = min0(mmax,n)
        do m=-nm,nm
          vintor(n,m) = fac*vextor(n,m)
          betaint(n,m) = -coa2np1*betaext(n,m)
        enddo ! m
      enddo ! n
!
! Find k index for nearest coefficients below 400 km (4x10^5 m)
      do k=1,nhgt_fix_r
        if (hgt_fix_r(k).ge.h_leo) then
          kbelow = k - 1
          frac = (hgt_fix_r(k) - h_leo)/(hgt_fix_r(k)-hgt_fix_r(kbelow))
          goto 100
        endif
      enddo ! k
  100 continue
      if (frac.lt.0.or.frac.gt.1.) then
        write (6,*) 'Stopped because frac =', frac
        stop
      endif
!
! Calculate global fields at ground and 400 km
      delbegrd_qd = 0.
      delbngrd_qd = 0.
      delbugrd_qd = 0.
      delbe400_qd = 0.
      delbn400_qd = 0.
      delbu400_qd = 0.
      psi_qd      = 0.
!      
      aob = re/(re+h_leo)
      do i=1,nmlon
        do l=1,nlat_qd-1
!  Approximate glat,glon at ground by values at base of dynamo region.
          glat = qd(i,l)%glat(1)
          glon = qd(i,l)%glon(1)
          call sphhar(glat,glon,ST)
          fac = ST/sq4pi
          aobnp1 = aob
          do n=1,nmax
            nm = min0(mmax,n)
            aobnp1 = aobnp1*aob
            do m=-nm,nm
              ma = iabs(m)
              vorsum = vextor(n,m) + vintor(n,m)
              betasum = betaext(n,m) + betaint(n,m)
              delbegrd_qd(i,l) = delbegrd_qd(i,l) &
                - vorsum*m*pnmost(n,ma)*fm(-m)/sq4pi
              delbngrd_qd(i,l) = delbngrd_qd(i,l) &
                + vorsum*dpnm(n,ma)*fm(m)/sq4pi
              delbugrd_qd(i,l) = delbugrd_qd(i,l) &
                + betasum*fac*pnmost(n,ma)*fm(m)
              psi_qd(i,l) = psi_qd(i,l) &
                + psicoef(n,m)*fac*pnmost(n,ma)*fm(m)
              vorsum = vintor(n,m)*aobnp1 + frac*vor(n,m,kbelow) &
                + (1.-frac)*vor(n,m,kbelow+1)
              betasum = betaint(n,m)*aobnp1*aob + frac*beta(n,m,kbelow) &
                + (1.-frac)*beta(n,m,kbelow+1)
              torinterp = frac*tor(n,m,kbelow) + (1.-frac)*tor(n,m,kbelow+1)
              delbe400_qd(i,l) = delbe400_qd(i,l) &
                - vorsum*m*pnmost(n,ma)*fm(-m)/sq4pi &
                - torinterp*dpnm(n,ma)*fm(m)/sq4pi
              delbn400_qd(i,l) = delbn400_qd(i,l) &
                + vorsum*dpnm(n,ma)*fm(m)/sq4pi  &
                - torinterp*pnmost(n,ma)*m*fm(-m)/sq4pi 
              delbu400_qd(i,l) = delbu400_qd(i,l) &
                + betasum*fac*pnmost(n,ma)*fm(m)
            enddo ! m
          enddo ! n
! Compute QD phi,lambda components of delB from geographic components.
! For simplicity, approximate QD vertical component of delB by radial
!  component of delB (i.e., approximate F*f3.delB by delBr).
          fac = (qd(i,l)%f11(1)*delbegrd_qd(i,l)   &
               + qd(i,l)%f12(1)*delbngrd_qd(i,l))/qd(i,l)%Fqd(1)
          delbngrd_qd(i,l) = (qd(i,l)%f21(1)*delbegrd_qd(i,l)  &
                         + qd(i,l)%f22(1)*delbngrd_qd(i,l))/qd(i,l)%Fqd(1)
          delbegrd_qd(i,l) = fac
	  
          f11interp = frac*qd(i,l)%f11(kbelow) + (1.-frac)*qd(i,l)%f11(kbelow+1)
          f12interp = frac*qd(i,l)%f12(kbelow) + (1.-frac)*qd(i,l)%f12(kbelow+1)
          f21interp = frac*qd(i,l)%f21(kbelow) + (1.-frac)*qd(i,l)%f21(kbelow+1)
          f22interp = frac*qd(i,l)%f22(kbelow) + (1.-frac)*qd(i,l)%f22(kbelow+1)
          Fqdinterp = frac*qd(i,l)%Fqd(kbelow) + (1.-frac)*qd(i,l)%Fqd(kbelow+1)
          fac = (f11interp*delbe400_qd(i,l)  &
               + f12interp*delbn400_qd(i,l))/Fqdinterp
          delbn400_qd(i,l) = (f21interp*delbe400_qd(i,l) &
                         + f22interp*delbn400_qd(i,l))/Fqdinterp
          delbe400_qd(i,l) = fac
        enddo ! l
      enddo ! i
! Duplicate values at QD longitude -pi for QD longitude +pi
!      do l=1,nlat_qd-1
!        delbegrd_qd(nmlon+1,l) = delbegrd_qd(1,l) 
!        delbngrd_qd(nmlon+1,l) = delbngrd_qd(1,l) 
!        delbugrd_qd(nmlon+1,l) = delbugrd_qd(1,l) 
!        delbe400_qd(nmlon+1,l) = delbe400_qd(1,l) 
!        delbn400_qd(nmlon+1,l) = delbn400_qd(1,l) 
!        delbu400_qd(nmlon+1,l) = delbu400_qd(1,l) 
!      enddo ! l
!
! Calculate fields at 0 QD longitude
      delbe0ln_qd = 0.
      delbn0ln_qd = 0.
      delbu0ln_qd = 0.
      
      lonqd2get = 105.*dtr
      mmin = 999.
      do l = 1,nmlon
        diff = abs(ylonm(l) - lonqd2get)
	if(diff.lt.mmin) then
	  mmin = diff
	  i = l
	endif
      enddo
      
      write(6,*) 'db found', i,ylonm(i)*rtd,lonqd2get*rtd
      
      do l=1,nlat_qd-1
        do k=1,nhgt_fix_r-1
          glat = qd(i,l)%glat(k)
          glon = qd(i,l)%glon(k)
          call sphhar(glat,glon,ST)
          fac = ST/sq4pi
          aob = re/(re+hgt_fix_r(k)) 
          aobnp1 = aob
          do n=1,nmax
            nm = min0(mmax,n)
            aobnp1 = aobnp1*aob
            do m=-nm,nm
              ma = iabs(m)
              vorsum = vintor(n,m)*aobnp1 + vor(n,m,k)
              betasum = betaint(n,m)*aobnp1*aob + beta(n,m,k)
              delbe0ln_qd(l,k) = delbe0ln_qd(l,k) &
                - vorsum*m*pnmost(n,ma)/sq4pi &
                - tor(n,m,k)*dpnm(n,ma)*fm(m)/sq4pi
              delbn0ln_qd(l,k) = delbn0ln_qd(l,k) &
                + vorsum*dpnm(n,ma)*fm(m)/sq4pi &
                - tor(n,m,k)*pnmost(n,ma)*m*fm(-m)/sq4pi
              delbu0ln_qd(l,k) = delbu0ln_qd(l,k) &
                + betasum*fac*pnmost(n,ma)*fm(m)
            enddo ! m
          enddo ! n
	  !
	  ! calculate scalar magnetic field
	  delscllH_qd(l,k) = delbe0ln_qd(l,k)*qd(i,l)%bhat(1,k) + &
	    delbn0ln_qd(l,k)*qd(i,l)%bhat(2,k) + &
	    delbu0ln_qd(l,k)*qd(i,l)%bhat(3,k)
	  !
! Compute QD phi,lambda components of delB from geographic components.
! For simplicity, approximate QD vertical component of delB by radial
!  component of delB (i.e., approximate F*f3.delB by delBr).
          fac = (qd(i,l)%f11(k)*delbe0ln_qd(l,k) &
               + qd(i,l)%f12(k)*delbn0ln_qd(l,k))/qd(i,l)%Fqd(k)
          delbn0ln_qd(l,k) = (qd(i,l)%f21(k)*delbe0ln_qd(l,k) &
                         + qd(i,l)%f22(k)*delbn0ln_qd(l,k))/qd(i,l)%Fqd(k)
          delbe0ln_qd(l,k) = fac
        enddo ! k
      enddo ! l
!      
      end subroutine calc_Bqd
!---------------------------------------------------------------------- 
      subroutine sphhar(glat_loc,glon_loc,ST)
! Inputs
!  glat = geographic latitude in degrees
!  glon = geographic longitude in degrees
      integer :: m,n
      real :: TH,CT,ST,PM2,PH,CP,SP,glat_loc,glon_loc,glat_mod
! Avoid poles, because of 1/ST.
      glat_mod = amax1(glat_loc,-90.+1.e-4)
      glat_mod = amin1(glat_mod,90.-1.e-4)
      TH = (90.-glat_mod)/rtd
      CT = cos(TH)
      ST = sin(TH)
! pnmost = Pnm/ST, where Pnm is an associated Legendre polynomial,
!  fully normalized as in Richmond (1974).
      pnmost = 0.
      dpnm = 0.
      pnmost(0,0) = 1./ST
      do m=0,mmax
        if (m.ne.0) pnmost(m,m) = PMOPMMO(m)*pnmost(m-1,m-1)*ST
        dpnm(m,m) = m*CT*pnmost(m,m)
        PM2 = 0.
        do n=m+1,nmax
          pnmost(n,m) = (CT*pnmost(n-1,m) - R(n-1,m)*PM2)/R(n,m)
          PM2 = pnmost(n-1,m)
          dpnm(n,m) = n*CT*pnmost(n,m) - float(2*n + 1)*R(n,m)*PM2
        enddo ! n
      enddo ! m

      PH = glon_loc/rtd
      CP = cos(PH)
      SP = sin(PH)
      fm(0) = 1.
      fm(-1) = SQ2*CP
      fm( 1) = SQ2*SP
      do m=2,mmax
        fm(-m) = CP*fm(1-m) - SP*fm(m-1)
        fm( m) = CP*fm(m-1) + SP*fm(1-m)
      enddo

      end subroutine sphhar
!---------------------------------------------------------------------- 
    end module delB_module
!---------------------------------------------------------------------- 
