       module area_factors_module
       use params_module,only: nmlatS2_h,hgt_fix, &
       nmlat_h,nhgt_fix, nhgt_fix_r, &
          hgt_fix_r,nhgt_fix_r,r0,ylonm_s,rho_s,re,pi
       implicit none
       real ::  m1f(nmlat_h,nhgt_fix)   ! M1*F
       real ::  m2f(nmlat_h,nhgt_fix)   ! M2*F
       real ::  m3f(nmlat_h,nhgt_fix_r) ! M3*F
! a1,a3 vary from 0 at the magnetic pole to 1 at the magnetic equator.
       real ::  a1(nmlat_h+1,nhgt_fix)  ! Normalized integral from pole of M1*F
       real ::  a3(nmlat_h+1,nhgt_fix_r)! Normalized integral from pole of M3*F
       contains
!----------------------------------------------------------------------
       subroutine calc_mf
! Calculation of normalized integrated areas from the pole to an s2
!  surface in the meridional (a1) and horizontal (a3) planes, and of
!  the factors m1f,m2f,m3f, which are independent of magnetic longitude.
!  These latter factors give M1,M2,M3 when divided by F.
       implicit none
!       
! local variables         
       integer ::  j,k,isn,jmax      
       real ::  fac1,fac2,fac3,dlonm,rm,rp,ra,rbar
!                   
       dlonm = ylonm_s(2)-ylonm_s(1) ! assumes equidistant longitudinal gridpoints
! Assumes hemispherical symmetry
        isn = 1 !Should get the same result for isn = 2, since rho_s(j,1) = rho_s(j,2)
!  
! In order to include the pole, the first index of a1,a3 represents the
!  location j-0.5, not j+0.5. However, the first index of m1f,m2f,m3f
!  represents the location j+0.5, i.e., the s2 points.
! 
! Set a1,a3 to 1 for equator and beyond (points before equator will be
!  overwritten later).
        a1 = 1.
        a3 = 1.
! Set m1f,m2f,m3f to 0 beyond equator (points before equator will be
!  overwritten later).
        m1f = 0.
        m2f = 0.
        m3f = 0.
	do k=1,nhgt_fix	
           a1(1,k) = 0.
           a3(1,k) = 0.
           jmax = nmlatS2_h - k + 1 ! number of s2 points at level k 
           if (jmax.lt.1) then
! Error trap; this condition should never occur.
              write(6,*) 'Stopped in calc_mf because jmax=',jmax
              stop
           endif
! Normalized radii of the top and bottom of layer k.
           rp = (hgt_fix_r(k+1)+re)/r0           ! r_k+0.5 /R
           rm = (hgt_fix_r(k)+re)/r0             ! r_k-0.5 /R
!
           fac2 = 0.5*pi*(rp - rm)*r0   ! Pi/2*((r_k+0.5)-(r_k-0.5)) 
	   fac3 = 2.*dlonm*(hgt_fix(k)+re)**3/r0 ! 2*dlon*r_k^3/R
	   do j=1,jmax     ! loop over all s2 points
	      fac1 = sqrt(rm)*rho_s(j,isn) 	! sqrt[r_k-0.5/R]*rho(j+0.5) 
! First index of a1,a3 is j+1 because this corresponds to position j+0.5.
              a3(j+1,k) = 1. - sqrt(1.-fac1**2)
              m3f(j,k) = (hgt_fix_r(k)+re)**2 *dlonm*(a3(j+1,k)-a3(j,k))
! Calculate the normalized radius within the layer, rbar, that gives the most
!  accurate value of a1 when a1 is computed by the approximation below.
!  This calculation of rbar assumes the radius of field lines near the
!  equator is parabolic with respect to magnetic latitude.
              ra = 1./rho_s(j,isn)**2
! Prevent ra from getting large enough to affect the numerical accuracy
!  of rbar, which rapidly asymptotes to .5*(rp+rm) as ra increases.
              ra = min(ra,rp+16.*(rp-rm))
	      ra = max(ra,rp)  ! ra was less than rp when it should have been equal
	      !
              rbar = ra - (2.*(sqrt(ra-rm)**3-sqrt(ra-rp)**3)/(3.*(rp-rm)))**2  ! eq. (45') page 4c r*/R
              a1(j+1,k) = 2.*asin(sqrt(rbar)*rho_s(j,isn))/pi                   ! eq. (47')
              m1f(j,k) = ((hgt_fix(k)+re)/r0)**2.5 *fac2*(a1(j+1,k)-a1(j,k))*r0
! Make sure fac1 does not numerically exceed 1,
!  so that sqrt(1-fac1**2) can be computed.
              fac1 = sqrt(rp)*rho_s(j,isn)        ! sqrt[r_k+0.5/R]*rho(j+0.5) 
              fac1 = min(fac1,1.)
              m2f(j,k) = fac3*(sqrt(1.-rm*rho_s(j,isn)**2) &
                - sqrt(1.-fac1**2))*sqrt(1.-.75*rho_s(j,isn)**2)/rho_s(j,isn)
           enddo !j
            m1f(jmax+1,k) = ((hgt_fix(k)+re)/r0)**2.5 *fac2*(1-a1(jmax+1,k))*r0
           m3f(jmax+1,k) = (hgt_fix_r(k)+re)**2 *dlonm*(1.-a3(jmax+1,k))
        enddo !k
!
! Now do a3,m3f for top level
        k = nhgt_fix_r
        a3(1,k) = 0.
        jmax = nmlatS2_h - k + 1 ! number of s2 points at k-0.5 
        if (jmax.ge.1) then
           rm = (hgt_fix_r(k)+re)/r0
	   do j=1,jmax     ! loop over all s2 points (level k has nmlatS2_h-k+1 s2 points)
	      fac1 = sqrt(rm)*rho_s(j,isn) 	! sqrt[r_k-0.5/R]*rho(j+0.5) 
              a3(j+1,k) = 1. - sqrt(1.-fac1**2)
              m3f(j,k) = (hgt_fix_r(k)+re)**2 *dlonm*(a3(j+1,k)-a3(j,k))
            enddo !j
            m3f(jmax+1,k) = (hgt_fix_r(k)+re)**2 *dlonm*(1.-a3(jmax+1,k))
        endif
!       
       end subroutine calc_mf
!----------------------------------------------------------------------
       end module area_factors_module
!----------------------------------------------------------------------
