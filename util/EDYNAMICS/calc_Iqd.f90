!----------------------------------------------------------------------------- 
      subroutine calc_Iqd
! 
      use params_module, only: nlat_qd,nlat_qd_h,nhgt_fix,nmlon,nmlat_h,rtd,nhgt_fix_r,re, &
         hgt_fix_r,hgt_fix,pi,r0,rho_s
      use qd_module, only: qd_grid,qd,a1q,wgt1,jl_qd,a3q,wgt3,jl3_qd,lat_qd_ed, &
         LI1qd,I1qd,LI3qd,I3qd,I2qd,Jf1qd,Jf2qd,Jrqd,hgt_qd_mp,hgt_qd_ed,lat_qd_mp,lon_qd_mp, &
	 Jr,Jf1,Jf2,g13,g23,Jf1hor,Jf2hor,Jeej,q1,q2,q3
      use fieldline_s_module, only: fieldline_s1,fline_s1,fieldline_s2,fline_s2 
      use fieldline_r_module, only: fieldline_r,fline_r 
      use fieldline_p_module, only: fieldline_p,fline_p 
!     
!     
      implicit none
!     
      integer :: i,j,k,l,ll,isn,nlatmax,jstart,ip,im
      real :: fac,fac1,fac2,fac2p,fac2m,dlatq,dlonq,fac1m,fac1p,qp,qm,dlonm,I3_eq
      real :: a1q_upd(nlat_qd_h,nhgt_fix),&       ! A1q(l)/A1q(nlat_max) (updated 8/3/2015)
              a3q_upd(nlat_qd_h,nhgt_fix_r)       ! A3q(l)/A3q(nlat_max) (updated 8/3/2015)
      real :: LI1(nmlon,nmlat_h+1,nhgt_fix,2)     ! for mapping from mod.apex to quasi dipole grid latitudinal integrated I1
      real :: LI3(nmlon,nmlat_h+1,nhgt_fix_r,2)   ! for mapping from mod.apex to quasi dipole grid latitudinal integrated I3
      real :: lamqd_from_apex_coord, &   ! function
              latqd_tmp,absf1,f1f2
!	      
      real,parameter :: hgt_max = 1.5e5  ! [m] for testing
      integer:: kmax_r

! ADR note 2015/10/8: Since LI1,LI3 are only used temporarily, they could be
!  made arrays with the single index j, if do loops are rearranged.  LI1qd,
!  LI3qd could also be 1D arrays, if I2qd is calculated in multiple steps.
!      
! calculated integrated current I1 (cannot be removed; same grid as a1 and aa1)
      LI1 =0.
! isn=1 in S hemisphere; isn=2 in N hemisphere
      do isn=1,2
	do i=1,nmlon
	  do k=1,nhgt_fix ! maximum of points on fieldline
            LI1(i,1,k,isn) =  0.   ! polar values
            nlatmax = nmlat_h-k+1  ! each height level loses one lat. grid point
	    do j=2,nlatmax+1  ! loop over all s2 points plus equator
                LI1(i,j,k,isn) = LI1(i,j-1,k,isn) + fline_s1(i,j-1,isn)%I1(k)
	    end do 
	  end do 
	end do
      end do 
!
! mapping of I1 from rho to qd coordinate system
      do isn=1,2
	do i=1,nmlon
	  do k=1,nhgt_fix ! maximum of points on fieldline
	    do l=1,nlat_qd_h-1  ! loop over all qd latitude edge points
	        j = jl_qd(l,k)
	        LI1qd(i,l,k,isn) = (1-wgt1(l,k))*LI1(i,j,k,isn)+ wgt1(l,k)*LI1(i,j+1,k,isn)
	    end do 
	    ! set equator value which is the last point for both grids so L1qd = LI1(rho) (assumes we always have an rho point at equator,
	    !                                                                              and we set up our grid that way)
	    j = nmlat_h-k+1  
	    l=  nlat_qd_h 
	    LI1qd(i,l,k,isn) = LI1(i,j+1,k,isn) 
	    !
	    ! calculate I1qd at midpoints l which is the difference of LI1qd from the adjacent edges)
	    ! eq. 242
	    do l=1,nlat_qd_h-1  ! loop over all qd latitude midpoints
	        I1qd(i,l,k,isn) = LI1qd(i,l+1,k,isn)-LI1qd(i,l,k,isn)  ! difference between integrated values is I1qd
	    end do 
	  end do
	end do
      end do 
    
!  For I3qd:
! calculated integrated current I3
      LI3 =0.
      do isn=1,2
	do i=1,nmlon
	  do k=1,nhgt_fix_r ! maximum height level
            LI3(i,1,k,isn) =  0.   ! polar values
            nlatmax = nmlat_h-k+1  ! each height level loses one lat. grid point
	    do j=2,nlatmax  ! loop over s2 points
                LI3(i,j,k,isn) = LI3(i,j-1,k,isn) + fline_r(i,j-1,isn)%I3(k)
		latqd_tmp = lamqd_from_apex_coord(fline_s2(i,j-1,isn)%mlat_m,hgt_fix_r(k))
		!if(i.eq.20.and.k.eq.2) write(33,'(i4,2(x,e15.8))') l,latqd_tmp,LI3(i,l,k,isn)
	    end do 
	    ! 
!
!
!
! 2015-10-15 ADR: This calculation could be done within the above do
!   loop over j, if the do-loop upper limit is changed to nlatmax+1.
            j= nlatmax+1  ! integrated up to equator
            LI3(i,j,k,isn) = LI3(i,j-1,k,isn) + fline_r(i,j-1,isn)%I3(k)
!
!
!
	  end do
	end do
      end do 
      
      ! eq. 63' on page 7 for equator value 
      do isn=1,2
         do i=1,nmlon
          if(i == 1) then
            im = nmlon ! use wrap around point as i-1 -> nmlon 
          else
            im = i-1
          endif
          do k=1,nhgt_fix ! maximum of points on fieldline
             j = nmlat_h-k+1  ! each height level loses one lat. grid point
	     I3_eq = fline_r(i,j,isn)%I3(k)+fline_s1(im,j,isn)%I1(k)-fline_s1(i,j,isn)%I1(k)+fline_s2(i,j-1,isn)%I2(k)
	  end do
	end do
      end do 
!
! mapping of I3 from rho to qd coordinate system
      do isn=1,2
	do i=1,nmlon
	  do k=1,nhgt_fix_r ! maximum of points on fieldline
	    do l=1,nlat_qd_h-1  ! loop over all qd latitude edge points
	        j = jl3_qd(l,k)
	        LI3qd(i,l,k,isn) = (1-wgt3(l,k))*LI3(i,j,k,isn)+ wgt3(l,k)*LI3(i,j+1,k,isn)
	    end do 
	    ! set equator value which is the last point for both grids so L3qd = LI3(rho) (assumes we always have an rho point at equator,
	    !                                                                              and we set up our grid that way)
	    j = nmlat_h-k+1
	    l=  nlat_qd_h 
	    LI3qd(i,l,k,isn) = LI3(i,j+1,k,isn) 
	    !
	    ! calculate I3qd at midpoints l which is the difference of LI3qd from the adjacent edges in latitude)
	    ! (eq. 246)
	    do l=1,nlat_qd_h-1  ! loop over all qd latitude midpoints
	        I3qd(i,l,k,isn) = LI3qd(i,l+1,k,isn)-LI3qd(i,l,k,isn)  
	    end do 
	  end do
	end do
      end do 
!
! calculate I2qd as the divergence of L1qd and L3qd; L1qd longitudes are those of S1 grid with lon_s1(1) >  lon_s2(1)
! I2q(i,l+0.5,k) = L1q(i-0.5,l+0.5,k) - L1q(i+0.5,l+0.5,k)+
!                  L3q(i,l+0.5,k-0.5) - L3q(i,l+0.5,k+0.5)
      do isn=1,2
       do i=1,nmlon
	 im = i-1
	 if(i.eq.1) im = nmlon
	 do k=1,nhgt_fix ! maximum of points
	   I2qd(i,1,k,isn) = 0.  ! pole value
	   do l=2,nlat_qd_h      ! loop over all qd latitude edge points l+0.5
	     I2qd(i,l,k,isn) = LI1qd(im,l,k,isn)-LI1qd(i,l,k,isn)+ &
		      LI3qd(i,l,k,isn)-LI3qd(i,l,k+1,isn)	      
! At this point I2qd is equatorward current.  Change sign in N hemisphere to
!  make it be northward current in both hemispheres.
             if (isn.eq.2) I2qd(i,l,k,isn) = -I2qd(i,l,k,isn)
	   end do
	 end do
       end do
      end do 
      !
!
! At different points! eq. 254-256
! Jf1(i-0.5,l,k) = 0.5*I1qd(i-0.5,l,k)/r(k)^1.5/(r(k+0.5)^0.5-r(k-0.5)^0.5)/dlamq
! Jf2(i,l+0.5,k) = I2qd(i,l+0.5,k)/r(k)/(r(k+0.5)-r(k-0.5))/dlonq/cos(latq(l+0.5))
! Jr(i,l,k-0.5)  = I3qd(i,l,k-0.5)/M3qd(i,l,k-0.5)
!
      dlatq =  abs(lat_qd_mp(2)-lat_qd_mp(1))  ! regular grid spacing
      dlonq =  abs(lon_qd_mp(2)-lon_qd_mp(1))  ! regular grid spacing
      do k=1,nhgt_fix 
        fac1 = (re+hgt_qd_mp(k))*(hgt_qd_ed(k+1)-hgt_qd_ed(k))*dlonq  ! r(k)*(r(k+0.5)-r(k-0.5))*dlonq
 	fac1 = 1/fac1
        fac2 = 2*(re+hgt_qd_mp(k))**1.5*((hgt_qd_ed(k+1)+re)**0.5-(hgt_qd_ed(k)+re)**0.5)*dlatq  ! 2*r(k)^1.5*(r(k+0.5)^0.5-r(k-0.5)^0.5)*dlamq
	fac2 = 1/fac2
	!
        do i=1,nmlon
          l = 1
          Jf2qd(i,l,k,:) = 0.
          do l=2,nlat_qd_h  ! loop over all qd latitude midpoints (1 to 80)
	      isn = 1
             !Jf1qd(i,l-1,k,isn) = I1qd(i,l-1,k,isn)*fac2    
             Jf1qd(i,l-1,k,isn) = I1qd(i,l-1,k,isn)/q1(l-1,k)
             !Jf2qd(i,l,k,isn)	= I2qd(i,l,k,isn)*fac1/cos(lat_qd_mp(l))  ! divide by cos(latq(l)
             Jf2qd(i,l,k,isn)	= I2qd(i,l,k,isn)/q2(l,k)
             Jrqd(i,l-1,k,isn)  = I3qd(i,l-1,k,isn)/qd(i,l)%M3q(k)	
	     !   
	      isn = 2
             ll = nlat_qd-1-l	! north pole to equator;  160 to 81 (midpoints) where M3q is stored
             !Jf1qd(i,l-1,k,isn) = I1qd(i,l-1,k,isn)*fac2      
             Jf1qd(i,l-1,k,isn) = I1qd(i,l-1,k,isn)/q1(l-1,k)
             !Jf2qd(i,l,k,isn)	= I2qd(i,l,k,isn)*fac1/cos(lat_qd_mp(l))  ! divide by cos(latq(l)
             Jf2qd(i,l,k,isn)	= I2qd(i,l,k,isn)/q2(l,k)
             Jrqd(i,l-1,k,isn)  = I3qd(i,l-1,k,isn)/qd(i,ll)%M3q(k)   
           end do     
         end do
       end do
      !      	
      k=nhgt_fix_r 
        do i=1,nmlon
          do l=1,nlat_qd_h-1  ! loop over all qd latitude midpoints (1 to 80) in one hemisphere
	     isn = 1
             Jrqd(i,l,k,isn) = I3qd(i,l,k,isn)/qd(i,l)%M3q(k)  	   
	     isn = 2  
             ll = nlat_qd-l+1-1 ! north pole to equator;  160 to 81 (midpoints) where M3q is stored	   
             Jrqd(i,l,k,isn) = I3qd(i,l,k,isn)/qd(i,ll)%M3q(k)  	   
          end do 
        end do
! 
! average to one i,l,k point   Jf1qd & Jrqd(i,l,k,isn) arrays are  Jf1qd_tmp,Jr,Jf2qd all at one point
      do k=1,nhgt_fix 
        do isn=1,2
          do i=1,nmlon
	    im = i-1
	    if(i.eq.1) im = nmlon
	    do l=1,nlat_qd_h-1  ! loop over all qd latitude midpoints
	       Jf1(i,l,k,isn) = 0.5*(Jf1qd(i,l,k,isn)+Jf1qd(im,l,k,isn) )   ! S1 points start later than S2 (-180) points in longitude
	       Jf2(i,l,k,isn) = 0.5*(Jf2qd(i,l,k,isn)+Jf2qd(i,l+1,k,isn) ) 
	       Jr(i,l,k,isn)  = 0.5*(Jrqd(i,l,k,isn)+ Jrqd(i,l,k+1,isn)) 
	   end do 
	 end do
       end do
      end do
      
! find index below hgt_max 
!      k=1
!      do while (hgt_max.ge.hgt_fix(k))
!         kmax_r = k
!	 k      = k+1
!      end do
!      	      
!      Jf1(:,:,kmax_r+1:nhgt_fix,:) = 0.
!      Jf2(:,:,kmax_r+1:nhgt_fix,:) = 0.
!      Jr(:,:,kmax_r+1:nhgt_fix,:)  = 0
!      .
!      Jf1(:,:,1:kmax_r,:) = 0.
!      Jf2(:,:,1:kmax_r,:) = 0.
!      Jr(:,:,1:kmax_r,:)  = 0.
! 
! calculate horizontal current J = Jf1hor*f1+Jf1hor*f2+Jr
! Jf1hor = Jf1-(k^ dot g1)Jr (eq. 261)
! Jf2hor = Jf2-(k^ dot g2)Jr (eq. 262)
!
!
      do k=1,nhgt_fix 
        do isn=1,2
          do i=1,nmlon
	    do l=1,nlat_qd_h-1  ! loop over all qd latitude midpoints
	       Jf1hor(i,l,k,isn) = Jf1(i,l,k,isn) - g13(i,l,k,isn)*Jr(i,l,k,isn)
	       Jf2hor(i,l,k,isn) = Jf2(i,l,k,isn) - g23(i,l,k,isn)*Jr(i,l,k,isn)
	   end do 
	 end do
       end do
      end do 
!       
! calculate J = (Jf1hor*f1+Jf1hor*f2)*f1/F
! Jf1hor = Jf1-(k^ dot g1)Jr (eq. 261)
! Jf2hor = Jf2-(k^ dot g2)Jr (eq. 262)
!
!
      do k=1,nhgt_fix 
        do isn=1,2
          do i=1,nmlon
	    do l=1,nlat_qd_h-1  ! loop over all qd latitude midpoints
	       absf1 = sqrt(qd(i,j)%f11(k)**2+qd(i,j)%f12(k)**2)
	       f1f2 = (qd(i,j)%f11(k)*qd(i,j)%f21(k)+qd(i,j)%f12(k)*qd(i,j)%f22(k))
	       Jeej(i,l,k,isn) = Jf1(i,l,k,isn)*absf1 + Jf2(i,l,k,isn)*f1f2/absf1 	           
	   end do 
	 end do
       end do
      end do     
!         
      end subroutine calc_Iqd
!----------------------------------------------------------------------------- 
