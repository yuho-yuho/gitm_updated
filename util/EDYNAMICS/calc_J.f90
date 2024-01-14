!
       subroutine calc_J
!        
       use fieldline_s_module,only: fieldline_s1,fline_s1, &
          fieldline_s2,fline_s2
       use fieldline_p_module,only:  fieldline_p,fline_p
       use fieldline_r_module,only:  fieldline_r,fline_r
       use params_module,only: nmlon,nmlat_h,rho,rho_s,val_fill
!  
      implicit none
! local
     integer :: isn,i,j,k,im
     real :: I3,Ir,I2p_loc,M2p_loc 
     real,parameter :: sqrt2  = 2.**0.5, &
          fac1  = 0.5- 0.5**1.5, &
          fac2 = sqrt(0.5)-0.5
!
!  how does the index in notes relate to the index in code
!    points     notes        code        quantities
!    p-points   i,j,k        i,j,k                 e.g., potential, S
!   s1 points   i+ 0.5,j,k   i,j,k
!   s2 points   i,j+0.5,k    i,j,k
!   r points    i,j,k-0.5    i,j,k   
! 
     do isn = 1,2
       do i = 1,nmlon 
         if(i == 1) then ! wrap around in longitude
           im = nmlon
         else
           im = i-1	  
         endif
	 
!       
! Jr(i,j,k) = I3(i,j,k)/M3(i,j,k) (122) but
! top volume at equator with k=km 
! Jr(i,j,k) = sqrt(2) Ir(i,j,k)/M3(i,j,k-0.5) (124)
!    and Ir = (0.5-0.5^1.5)[I1(i-0.5,j,k)-I1(i+0.5,j,k)] -(0.5^0.5-0.5)I2(i,j-0.5,k)+0.5*I3(i,j,k-0.5) (123)
! 
! I1hor(i,j,k) = 0.5*I1(i-0.5,j,k)+0.5I1(i+0.5,j,k)-0.5*[M1(i-0.5,j,k)+M1(i+0.5,j,k)]*[Jr*kvec dot d1vec](i,j,k) (126)
! I2hor(i,j,k) = 0.5*I2(i,j-0.5,k)+0.5I2(i,j+0.5,k)-0.5*[M2(i,j-0.5,k)+M1(i,j+0.5,k)]*[Jr*kvec dot d2vec](i,j,k) (127)
!       
         do j = 2,nmlat_h  ! Not for poles
!	 
	   do k=1,fline_p(i,j,isn)%npts  ! P & S1 have same number of point
	       if((nmlat_h-j+1) == k) then ! top volume at equator k=km
!	         Ir = fac1*(fline_s1(im,j,isn)%I1(k)-fline_s1(i,j,isn)%I1(k)) - &
!		    fac2*fline_s2(i,j-1,isn)%I2(k)+0.5*fline_r(i,j,isn)%I3(k)
!	         fline_p(i,j,isn)%Jr(k) = sqrt2*Ir/fline_r(i,j,isn)%M3(k)
!
! am 4/16/2015
                 fline_p(i,j,isn)%Jr(k) = fline_r(i,j,isn)%I3(k)/sqrt2/fline_r(i,j,isn)%M3(k)		 
!		   
                 I2p_loc = 0.
                 M2p_loc = 0.
	       else
	         I3 = fline_r(i,j,isn)%I3(k)+fline_r(i,j,isn)%I3(k+1) ! calculate at p-level
	         I3 = I3*0.5
	         fline_p(i,j,isn)%Jr(k) = I3/fline_p(i,j,isn)%M3(k)
		 ! if(i.eq.10.and.k.eq.10) write(66,'(2(x,i4),3(x,e15.8))') j,isn,I3,fline_p(i,j,isn)%M3(k),fline_p(i,j,isn)%Jr(k)
                 I2p_loc =  fline_s2(i,j,isn)%I2(k)
                 M2p_loc =  fline_s2(i,j,isn)%M2(k)
	       endif
!	       
	       fline_p(i,j,isn)%I1hor(k) = 0.5*(fline_s1(im,j,isn)%I1(k)+fline_s1(i,j,isn)%I1(k)) - &
	          0.5*(fline_s1(im,j,isn)%M1(k)+fline_s1(i,j,isn)%M1(k))*fline_p(i,j,isn)%Jr(k)* &
		   fline_p(i,j,isn)%d1k(k)
		   
!		   if(i.eq.70.and.isn.eq.1) write(44,*) j,k,fline_p(i,j,isn)%d1k(k),fline_p(i,j,isn)%d2k(k) 
		   
	       fline_p(i,j,isn)%I2hor(k) = 0.5*(fline_s2(i,j-1,isn)%I2(k)+I2p_loc) - &
	          0.5*(fline_s2(i,j-1,isn)%M2(k)+M2p_loc)*fline_p(i,j,isn)%Jr(k)* &
		   fline_p(i,j,isn)%d2k(k) 
!		   
           enddo ! end of height loop
!	   
        enddo   ! end latitude loop
!	   
	j=1 ! pole values
	fline_p(i,j,isn)%I2hor(k) =  0.
	do k=1,fline_p(i,j,isn)%npts
	   fline_p(i,j,isn)%I1hor(k) =fline_p(i,j+1,isn)%I1hor(k)*rho_s(j,isn)/(rho_s(j+1,isn)-rho_s(j,isn))
        enddo ! end of height loop   
! 
! this is diagnostic at the moment      
! calculate values for R-points these are at (i,j,k+/-1/2)      
	 do j=2,nmlat_h   ! loop over all latitudes in one hemisphere 
           do k=1,fline_r(i,j,isn)%npts  
	     fline_r(i,j,isn)%Jr(k) = fline_r(i,j,isn)%I3(k)/fline_r(i,j,isn)%M3(k) 	     
           end do  ! end height loop
         end do  ! end lat/fieldline loop
	 j=1 ! pole values
         fline_r(i,j,isn)%Jr(:) = val_fill  ! put fill value in the polar flux tubes  
!	
       enddo   ! end longitude loop
      enddo   ! end of hemipshere loop
	     	     
!       	   
!       	   
       end subroutine calc_J
