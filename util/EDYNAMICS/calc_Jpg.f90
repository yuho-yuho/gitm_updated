!
       subroutine calc_Jpg(tei_s1,tei_s2,rho_s1,rho_s2,ne_s1,ne_s2)
!  
!  calculate gravity and plasma pressure gradiant driven current
! will be save in Je1Ion and Je2Ion   
! Note: * vertical gradient is approximated by assuming fieldline is vertical
!       * grad p is with respect to quasi dipole coordinates, therefore for -(grad p)xB/B2
!         it has to be transforemd to geographic derivatives
!            see (art's notes from 17 February 2015 page 1
!          Jp = -(grad p)xB/B2 =  d3 x (grad p) / Be3
!       je1p = d1 dot Jp = -e2 dot (grad p)/Be3 
!       je2p = d2 dot Jp =  e1 dot (grad p)/Be3   
! with (grad p) = (Re+h)^0.5/R^1.5 dp/dphi_qd  d1 + F / (Re+h) g2 * dp/d lam_qd + dp/dh k  ; d1, g2 and k are vectors
! 
!  => Je1p = -F/(Re+h)/Be3* e2 dot g2 dp/d lam_qd - e2 dot k / be3 dp/dh
!  
!  => Je2p =  (Re+h)^0.5/cos(lam_qd)/R^1.5/Be3 * dp/dphi_qd + F/(Re+h)/Be3* e1 dot g2 dp/d lam_qd + e1 dot k / be3 dp/dh
!
       use fieldline_s_module,only: fieldline_s1,fline_s1, &
         fieldline_s2,fline_s2
       use params_module,only: nmlon,nmlat_h,nmlatS2_h, &
                 boltz,matm,grav,re,r0,nhgt_fix
!
       real, intent(in) :: &
         tei_s1(nmlon,nmlat_h,nhgt_fix,2),tei_s2(nmlon,nmlatS2_h,nhgt_fix,2), &
         rho_s1(nmlon,nmlat_h,nhgt_fix,2),rho_s2(nmlon,nmlatS2_h,nhgt_fix,2), &
         ne_s1(nmlon,nmlat_h,nhgt_fix,2),ne_s2(nmlon,nmlatS2_h,nhgt_fix,2)
	
! local variables
       integer :: isn,i,j,k,nmax,ip,im,jp,jm,kp,km
       real :: fac,ipg(3),jxb(3),b_loc(3),d_loc(3), &
          dlon,dlat,dz,dNe,dTei,j_p1,j_p2,fsign,reh
    
      !  matm: mass per atomic weight
       fac = matm*grav   ! [m g/s2/mole] with  matm = 1.6605e-24 g/mole ; grav=  8.7 m/s2
       !
       ipg = 0. ! initialize
       !
       do isn = 1,2  ! loop over hemisphere
         if(isn.eq.1) then  ! for lat. derivative since indexing goes from pole to equator
	    fsign = 1.
	 else
	    fsign =-1.
	 endif   
         do i=1,nmlon ! loop over longitudes
          ! for S1 points
	   im = i -1
	   ip = i +1
	   if(i.eq.1)     im = nmlon 
	   if(i.eq.nmlon) ip = 1
	  
	  ! S1 points 
           do j=1,nmlat_h    ! loop over latitudes
	     jm = j -1       ! only one sided difference -poleward) since then I am sure we have a point at that height on the next field-line
	     jp = j 
	     if(j.eq.1)  then
	          jm = 1 
	          jp = 2
             endif 
	      
	     nmax = fline_s1(i,j,isn)%npts 
	     ipg = 0.
	     
      	     do k=1,nmax ! loop over height
	       
	       km = k -1
               kp = k +1
	       if(k.eq.1)    km = 1 
	       if(k.eq.nmax) kp = nmax 
	       if(km.eq.1.and.kp.eq.1) then   !lowest level at equator (only one height level) take poleward value
	          fline_s1(i,j,isn)%Je1Ion(k) = fline_s1(i,j-1,isn)%Je1Ion(k)
		  goto 100
	       endif
               !	      
               ! calculate ig = g*rho   contribution (up positive)
               ipg(3) = -rho_s1(i,j,k,isn)*fac*0.001     ! [kg/s2/m2]; 0.001 conversion from g to kg since B is in [T] [kg/s2/A]; rho [#/m3] 
               ! Jg = (ig)xB_vec/B2 with B in [T] = [kg/s2/A]
               b_loc(:) = fline_s1(i,j,isn)%bo(:,k)
               call cross_product(ipg,b_loc,jxb) 
	       jxb = jxb/fline_s1(i,j,isn)%Bmag(k)/fline_s1(i,j,isn)%Bmag(k)
       
               ! calculate contribution Je1 = d1 dot Jg
               d_loc(:) = fline_s1(i,j,isn)%d1(:,k)
               fline_s1(i,j,isn)%Je1Ion(k) =  dot_product(d_loc,jxb)     ! [A/m2]

               ! latitudinal gradient assume point is in the middle wrt latitude
	       ! tei * dNe/dlat
	       reh = fline_s1(i,j,isn)%hgt_pt(k)+re                                  ! [m]
	       dlat  =(fline_s1(i,jp,isn)%mlat_qd(k)-fline_s1(i,jm,isn)%mlat_qd(k))  ! lat differnce [rad]
	       dlat  = dlat*(reh)   ! [m]
	       
	       dNe = (ne_s1(i,jp,k,isn)-ne_s1(i,jm,k,isn))
	       j_p1 = tei_s1(i,j,k,isn)*dne/dlat
	       ! dtei/dlat * Ne
	       dTei = (tei_s1(i,jp,k,isn)-tei_s1(i,jm,k,isn))
	       j_p2 = ne_s1(i,j,k,isn)*dTei/dlat
	       
	       j_p1 =j_p1 +j_p2         !   1/(Re+h) * dp/d lam_qd 
	       j_p1 =j_p1*boltz          ! [kg/m2/s2]
	       
!               -F/(Re+h)/Be3* e2 dot g2 dp/d lam_qd   Note 1/(Re+h) is already included in dlat
	       fline_s1(i,j,isn)%Je1Ion(k) =fline_s1(i,j,isn)%Je1Ion(k) - &
	           j_p1*fline_s1(i,j,isn)%F(k)/fline_s1(i,j,isn)%be3(k)*fline_s1(i,j,isn)%e2g2(k)
	       
               ! height gradient (approximated (assuming vertical field-line)
	       ! central differencing
	       ! tei * dNe/dz
	       dz  = (fline_s1(i,j,isn)%hgt_pt(kp)-fline_s1(i,j,isn)%hgt_pt(km)) ! [m]
	      
	       dNe = (ne_s1(i,j,kp,isn)-ne_s1(i,j,km,isn))   ! [#/m3]
	       j_p1 = tei_s1(i,j,k,isn)*dne/dz               ! [K/m4]
	      
	       ! dtei/dz * Ne
	       dTei = (tei_s1(i,j,kp,isn)-tei_s1(i,j,km,isn))
	       j_p2 = ne_s1(i,j,k,isn)*dTei/dz              ! [K/m4]
	       j_p1 =j_p1 +j_p2
	       j_p1 =j_p1*boltz                             ! [kg/s2/m2]
	       
!  => Je1p = - e2 dot k / be3 dp/dh
	       fline_s1(i,j,isn)%Je1Ion(k) =fline_s1(i,j,isn)%Je1Ion(k) - &
	         j_p1*fline_s1(i,j,isn)%e2k(k)/fline_s1(i,j,isn)%be3(k)
	       !
	       ! diagnostic
               fline_s1(i,j,isn)%Ne(k)  = ne_s1(i,j,k,isn)
               fline_s1(i,j,isn)%Tei(k) = tei_s1(i,j,k,isn)
	       !
 100        enddo   ! end loop over height 
	    !
      	  enddo   ! end loop over latitudes 
        !
        ! end for S1 points
       
       ! for S2 points
           do j=1,nmlatS2_h    ! loop over latitudes
	     jm = j -1         ! only one sided difference -poleward) since then I am sure we have a point at that height on the next field-line
	     jp = j 
	     if(j.eq.1)  then
	          jm = 1 
	          jp = 2
             endif 
	      
	     nmax = fline_s2(i,j,isn)%npts 
	     ipg = 0.
	     
      	     do k=1,nmax ! loop over height
	       fline_s2(i,j,isn)%Je2Ion(k) = 0.
	       km = k -1
               kp = k +1
	       if(k.eq.1)    km = 1 
	       if(k.eq.nmax) kp = nmax 
	       if(km.eq.1.and.kp.eq.1) then   !lowest level at equator (only one height level) take poleward value
	          fline_s2(i,j,isn)%Je2Ion(k) = fline_s2(i,j-1,isn)%Je2Ion(k)
		  goto 200
	       endif
               !	      
               ! calculate ig = g*rho   contribution (up positive)
               ipg(3) = -rho_s2(i,j,k,isn)*fac*0.001     ! [kg/s2/m2]; 0.001 conversion from g to kg since B is in [T] 
               ! Jg = (ig)xB_vec/B2 with B in [T] = [kg/s2/A]
               b_loc(:) = fline_s2(i,j,isn)%bo(:,k)
               call cross_product(ipg,b_loc,jxb) 
	       jxb = jxb/fline_s2(i,j,isn)%Bmag(k)/fline_s2(i,j,isn)%Bmag(k)
       
               ! calculate contribution Je2 = d2 dot Jg
               d_loc(:) = fline_s2(i,j,isn)%d2(:,k)
               fline_s2(i,j,isn)%Je2Ion(k) =  dot_product(d_loc,jxb)     ! [A/m2]
               
	       ! calculate plasma pressure gradient contribution
	       !
	       ! longitudinal gradient 
               !   1/(Re+h)/cos(lam_qd)*[(Re+h)/R)^1.5/Be3 * dp/dphi_qd 
	       !
	       reh  = fline_s2(i,j,isn)%hgt_pt(k)+re
	       if(i.eq.1.or.i.eq.nmlon) then               ! lon differnce [rad]
	         dlon  =(fline_s2(3,j,isn)%mlon_qd(k)-fline_s2(1,j,isn)%mlon_qd(k))  ! regular in longitude (and there would be 360deg jump)
	       else
	         dlon  =(fline_s2(ip,j,isn)%mlon_qd(k)-fline_s2(im,j,isn)%mlon_qd(k))  ! lon differnce [rad]
	       endif
	       dlon  = dlon*(reh)*cos(fline_s2(i,j,isn)%mlat_qd(k))     ! [m]
	       ! tei * dNe/dlon
	       dNe   = (ne_s2(ip,j,k,isn)-ne_s2(im,j,k,isn))	 ! [#/m3]
	       j_p1  = tei_s1(i,j,k,isn)*dNe/dlon		 ! [K/m^4]
	       ! dtei/dlon * Ne
	       dTei  = (tei_s2(ip,j,k,isn)-tei_s2(im,j,k,isn))
	       j_p2  = ne_s2(i,j,k,isn)*dTei/dlon
	       j_p1  = j_p1 +j_p2
	       j_p1  = j_p1*boltz       !  1/(Re+h)/cos(lam_qd)* dp/dphi_qd 
	       
!               1/(Re+h)/cos(lam_qd)*[(Re+h)/R]^1.5/Be3 * dp/dphi_qd     Note 1/(Re+h)/cos(lam_qd) included in dlon
	       fline_s2(i,j,isn)%Je2Ion(k) =fline_s2(i,j,isn)%Je2Ion(k) + &
	           j_p1*(reh/r0)**1.5/fline_s2(i,j,isn)%be3(k)

               ! latitudinal gradient assume point is in the middle wrt latitude
	       ! tei * dNe/dlat
	       dlat  =(fline_s2(i,jp,isn)%mlat_qd(k)-fline_s2(i,jm,isn)%mlat_qd(k))  ! lat differnce [rad]
	       dlat  = dlat*(reh)                                                    ! [m]
	       
	       
	       dNe = (ne_s2(i,jp,k,isn)-ne_s2(i,jm,k,isn))
	       j_p1 = tei_s2(i,j,k,isn)*dNe/dlat
	       ! dtei/dlat * Ne
	       dTei = (tei_s2(i,jp,k,isn)-tei_s2(i,jm,k,isn))
	       j_p2 = ne_s2(i,j,k,isn)*dTei/dlat
	       
	       j_p1 =j_p1 +j_p2         !   1/(Re+h) * dp/d lam_qd 
	       j_p1 =j_p1*boltz          !  [kg/m2/s2]
	       
!              + F/(Re+h)/Be3* e1 dot g2 dp/d lam_qd  Note 1/(Re+h) is already included in dlat
	       fline_s2(i,j,isn)%Je2Ion(k) =fline_s2(i,j,isn)%Je2Ion(k) + &
	           j_p1*fline_s2(i,j,isn)%F(k)/fline_s2(i,j,isn)%be3(k)*fline_s2(i,j,isn)%e1g2(k)
	       
               ! height gradient (approximated (assuming vertical field-line)
	       ! central differencing
	       ! tei * dNe/dz
	       dz  = (fline_s2(i,j,isn)%hgt_pt(kp)-fline_s2(i,j,isn)%hgt_pt(km)) ! [m]
	      
	       dNe = (ne_s2(i,j,kp,isn)-ne_s2(i,j,km,isn))   ! [#/m3]
	       j_p1 = tei_s2(i,j,k,isn)*dNe/dz               ! [K/m4]
	      
	       ! dtei/dz * Ne
	       dTei = (tei_s2(i,j,kp,isn)-tei_s2(i,j,km,isn))
	       j_p2 = ne_s2(i,j,k,isn)*dTei/dz              ! [K/m4]
	       j_p1 =j_p1 +j_p2
	       j_p1 =j_p1*boltz                             ! [1/m2 kg/s2]
	       
!  => Je2p = + e1 dot k / be3 dp/dh
	       fline_s2(i,j,isn)%Je2Ion(k) =fline_s2(i,j,isn)%Je2Ion(k) + &
	          j_p1*fline_s2(i,j,isn)%e1k(k)/fline_s2(i,j,isn)%be3(k)
	       !
	       ! diagnostic
               fline_s2(i,j,isn)%Ne(k)  = ne_s2(i,j,k,isn)
               fline_s2(i,j,isn)%Tei(k) = tei_s2(i,j,k,isn)
	       !
 200        enddo   ! end loop over height 
	    !
      	  enddo   ! end loop over latitudes 

 
       ! end for S2 points
	    
      	  enddo   ! end loop over longitudes 
      	enddo   ! end loop over hemisphere
!	
       end subroutine calc_Jpg
