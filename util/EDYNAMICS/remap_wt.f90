!----------------------------------------------------------------------------- 
      subroutine remap_weights
! 
      use params_module, only: nlat_qd,nlat_qd_h,nhgt_fix,nmlat_h,nhgt_fix_r &
        ,re,ylonm_s,hgt_fix_r,hgt_fix,pi
      use qd_module, only: a1q,wgt1,jl_qd,a3q,wgt3,jl3_qd,lat_qd_ed &
        ,q1,q2,q3
      use area_factors_module, only: a1,a3 
!     
      implicit none
!     
      integer :: j,k,l,nlatmax,jstart
      real :: fac,dlonm,dlatm
      real :: a1q_upd(nlat_qd_h,nhgt_fix),& ! (updated 8/3/2015)
! a1q is the normalized area from the pole to the poleward edge of
!  element l, i.e., the equatorward edge of element l-1.
              a3q_upd(nlat_qd_h,nhgt_fix_r) ! (updated 8/3/2015)
! a3q is the normalized area from the pole to the poleward edge of
!  element l, i.e., the equatorward edge of element l-1.
!      
! a1q(l+0.5,k) = l/Lq (page 4 2015 April 20 bottom) at edges
! a1q is the normalized area from the pole to the poleward edge of
!  element l, i.e., the equatorward edge of element l-1.
      dlonm = ylonm_s(2)-ylonm_s(1) ! Assume constant in longitude.
      dlatm = abs(lat_qd_ed(2)-lat_qd_ed(1)) ! Assume constant in latitude.
      a1q_upd  = 0.
      do k=1,nhgt_fix
        q2(1,k) = 0.
        fac = (hgt_fix_r(k+1)-hgt_fix_r(k))*(re+hgt_fix(k))
        a1q_upd(1,k)  = 0  ! pole value 
        do l=2,nlat_qd_h   ! south pole to equator; edge max number nlat_qd; midpoint l max number nlat_qd-1
! Assume that the number of qd grid points is the same at all altitudes
!  (if this is not the case, the formula for a1q will be k-dependent).
!	  a1q_upd(l,k) = real((l-1))/real(nlat_qd_h-1)
	  a1q_upd(l,k) = abs((lat_qd_ed(l)-lat_qd_ed(1)))
	  ! q1(l,k) = rk*(r_k+0.5-r_k-0.5)*[A1q(l+0.5,k)-A1q(l-0.5,k)]
	  ! with A1q(l+0.5,k) = 1-2/pi*|lam_q(l+0.5)|
          q1(l-1,k) = fac*(a1q_upd(l,k)-a1q_upd(l-1,k))  ! (224a) page 56 Oct 8, 2015
          q2(l,k)   = fac*dlonm*cos(abs(lat_qd_ed(l)))     ! (252a) page 56 Oct 8, 2015
        end do
      end do
!  
! For each index l find the corresponding j (jl_qd) to linearly
!  interpolate LI1qd(l-.5,isn) between LI1(j-.5) and LI1(j+.5).
! isn,i,k are the same for the qd and rho grids.
      wgt1  = 9999.
      jl_qd = 9999
      a1q_upd = a1q_upd*2/pi ! normalize but reverse at end
      do k=1,nhgt_fix ! number of levels for I1 and I1qd
        jstart = 1	       ! when the latitude loop starts always start searching from the pole
        wgt1(1,k)  = 0.  ! for polar value since a1 and LI1 are zero
        jl_qd(1,k) = 1
! nlatmax is one plus the number of S2 points at level k in one hemisphere.
!  This is one less than the number of points used for the first index
!  of a1 and LI1.
        nlatmax = nmlat_h-k+1
! loop over poleward qd box edges, except for edges lying at pole and equator.
        lloop: do l=2,nlat_qd_h-1
          jloop: do j=jstart,nlatmax  ! always start searching from the last found j*
             !
             if(a1q_upd(l,k).ge.a1(j,k).and.a1q_upd(l,k).lt.a1(j+1,k)) then
               jl_qd(l,k) = j
               ! calculate weighting function  
               wgt1(l,k)  = (a1q_upd(l,k) - a1(j,k))/(a1(j+1,k)- a1(j,k))
               jstart = j  ! save last found value to start from that one
               exit jloop
             end if
             !
          end do jloop 
        end do lloop
        !
        l= nlat_qd_h
        jl_qd(l,k) = nlatmax
        wgt1(l,k)  = 1.
        !	
      end do
      a1q_upd = a1q_upd/2*pi ! normalize but reverse at end
!      
! updated 8/3/2015 equation (235')     
! a3q(l+0.5,k-0.5) = [1-sin(|lamq(l+0.5)|)]
!                 see (page 6 2015 April 20 top) at edges
! a3q is the normalized area from the pole to the poleward edge of
!  element l, i.e., the equatorward edge of element l-1.
      a3q_upd  = 0.
      do k=1,nhgt_fix_r
        fac = dlonm*(re+hgt_fix_r(k))**2
        a3q_upd(1,k)  = 0  ! pole value 
        do l=2,nlat_qd_h   ! south pole to equator; edge max number nlat_qd; midpoint l max number nlat_qd-1
! Assume that the number of qd grid points is the same at all altitudes
!  (if this is not the case, the formula for a1q will be k-dependent).
	  a3q_upd(l,k)  = 1-sin(abs(lat_qd_ed(l)))
          q3(l-1,k) = fac*(a3q_upd(l,k)-a3q_upd(l-1,k)) 
        end do
      end do
!
! For each index l find the corresponding j (jl3_qd) to linearly
!  interpolate LI3qd(l-.5,isn) between LI3(j-.5) and LI3(j+.5).
! isn,i,k are the same for the qd and rho grids.
      wgt3   = 9999.
      jl3_qd = 9999.
      do k=1,nhgt_fix_r ! number of levels for I3 and I3qd
        jstart      = 1	  ! when the latitude loop starts always start searching from the pole
        l = 1
        wgt3(l,k)   = 0.  ! for polar value since a3 and LI3 are zero
        jl3_qd(l,k) = 1
! Write out jl's and wgt's.
!          if (k.eq.1) write(8,'(2i5,f10.5,i5,f10.5)') l, &
!            jl_qd(l,k),wgt1(l,k),jl3_qd(l,k),wgt3(l,k)
        nlatmax     = nmlat_h-k+1
        !
! loop over poleward qd box edges, except for edges lying at pole and equator.
        lloop3: do l=2,nlat_qd_h-1
          jloop3: do j=jstart,nlatmax  ! always start searching from the last found j*
!             if(a3q(l,k).ge.a3(j,k).and.a3q(l,k).lt.a3(j+1,k)) then
             if(a3q_upd(l,k).ge.a3(j,k).and.a3q_upd(l,k).lt.a3(j+1,k)) then
               jl3_qd(l,k) = j
               ! calculate weighting function  
!               wgt3(l,k)  = (a3q(l,k) - a3(j,k))/(a3(j+1,k)- a3(j,k))
               wgt3(l,k)  = (a3q_upd(l,k) - a3(j,k))/(a3(j+1,k)- a3(j,k))
               jstart = j
               exit jloop3
             end if
             
          end do jloop3   
! Write out jl's and wgt's.
!          if (k.eq.1) write(8,'(2i5,f10.5,i5,f10.5)') l, &
!            jl_qd(l,k),wgt1(l,k),jl3_qd(l,k),wgt3(l,k)
        end do lloop3
! Now do equatorial value.
        l= nlat_qd_h
        jl3_qd(l,k) = nlatmax
        wgt3(l,k)  = 1.	
! Write out jl's and wgt's.
!          if (k.eq.1) write(8,'(2i5,f10.5,i5,f10.5)') l, &
!            jl_qd(l,k),wgt1(l,k),jl3_qd(l,k),wgt3(l,k)
      end do !k
!  
      end subroutine remap_weights
!----------------------------------------------------------------------------- 
