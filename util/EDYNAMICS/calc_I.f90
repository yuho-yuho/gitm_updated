!
       subroutine calc_I
!        
       use fieldline_s_module,only: fieldline_s1,fline_s1, &
          fieldline_s2,fline_s2
       use fieldline_p_module,only:  fieldline_p,fline_p
       use fieldline_r_module,only:  fieldline_r,fline_r
       use params_module,only: nmlon,nmlat_h,nmlatS2_h,val_fill, &
             use_stabil,J3LB
!       
       implicit none
       integer :: isn,i,j,k,nmax,ip,im
       real    :: I3eq,Je1A,Je1B
!       
! local variables         
!
! equation (83) page 10 Art's script
! I1(i+0.5,j,k) = N1p(i+0.5,j,k)*[Phi(i,j)-Phi(i+1,j)]
!      -N1h(i+0.5,j,k)[Phi(i,j-1)+Phi(i+1,j-1)-Phi(i,j+1)-Phi(i+1,j+1)]
!      +M1(i+0.5,j,k)Je1D(i+0.5,j,k)
!
! equation (85) page 10 Art's script
! I2(i,j+0.5,k) = N2h(i,j+0.5,k)*[Phi(i-1,j)+Phi(i-1,j+1)-Phi(i+1,j)-Phi(i+1,j+1)]
!      +N2p(i,j+0.5,k)[Phi(i,j)-Phi(i,j+1)]
!      +M2(i,j+0.5,k)Je2D(i,j+0.5,k)
!
! equation (63') page 7 Art's script
! for all i, and j=2,nmlat_h
!      I3(i,j,k+0.5) = I3(i,j,k-0.5) + I1(i-0.5,j,k)-I1(i+0.5,j,k)+I2(i,j-0.5,k)-I2(i,j+0.5,k)
! for k=1 (lowest level I3(i,j,0.5) is given (from lower atmosphere) variable J3LB in parms.f90
! 
! 
! 
! 2015-10-15 ADR: If J3LB has not been determined, is it somewhere set
!   to zero, or does the compiler set it to zero?
! 
! 
!  
!  how does the index in notes relate to the index in code
!    points     notes        code        quantities
!    p-points   i,j,k        i,j,k                 e.g., potential, S
!   s1 points   i-0.5,j,k  i-1,j,k
!   s2 points   i,j-0.5,k   i,j-1,k
!   r points    i,j,k-0.5   i,j,k
!
!    points     notes        code        quantities
!    p-points   i,j,k        i,j,k                 e.g., potential, S
!   s1 points   i+ 0.5,j,k   i,j,k
!   s2 points   i,j+0.5,k    i,j,k
!   r points    i,j,k-0.5    i,j,k   
       
       do isn = 1,2 ! loop over both hemispheres
         do i=1,nmlon ! loop over all longitudes
	   if(i == nmlon) then
	     ip = 1 ! use wrap around point as i+1 -> 1
	   else
	     ip = i+1
	   endif
           !
	   ! for the lowest equatorial volume
	   !    k == 1 .and. (nmlat_h-j+1) == k
	   ! I(i-0.5,j,k) = N_1^C(i-0.5,j,k)[Phi(i-1,j)-Phi(i,j)]+M1(i-0.5,j,k)Je1^S(i-0.5,j,k) (eq. 93) Art's notes
	   ! N1^c(i-0.5,j,k) = M1(i-0.5,j,k)[sig_c](i-0.5,j,k)/R/rhoj(phi_i-phi_i-1)
	   ! N1^c saved in N1^p
	   ! 
	   j=nmlat_h; k=1
	   !
	   fline_s1(i,j,isn)%I1(k) = fline_s1(i,j,isn)%N1p(k)*(fline_p(i,j,1)%pot-fline_p(ip,j,1)%pot) + &
	     fline_s1(i,j,isn)%M1(k)*fline_s1(i,j,isn)%Je1D(k)
	   !
           do j=2,nmlat_h-1 ! loop over all latitudes in one hemisphere NOT the POLE and EQUATOR
	      nmax = fline_s1(i,j,isn)%npts ! maximum of points on fieldline
              do k=1,nmax
!       
	        fline_s1(i,j,isn)%I1(k) = fline_s1(i,j,isn)%N1p(k)*(fline_p(i,j,1)%pot - fline_p(ip,j,1)%pot) &
		 -fline_s1(i,j,isn)%N1h(k)*(fline_p(i,j-1,1)%pot + fline_p(ip,j-1,1)%pot &
		                          - fline_p(i,j+1,1)%pot - fline_p(ip,j+1,1)%pot ) &
		 +fline_s1(i,j,isn)%M1(k)*fline_s1(i,j,isn)%Je1D(k)
	  	if(use_stabil) fline_s1(i,j,isn)%I1(k) = fline_s1(i,j,isn)%I1(k)+ &       ! considering cowling conductivity 
		     fline_s1(i,j,isn)%M1Je1D(k)
	        ! diagnostic
		fline_s1(i,j,isn)%I13d_1(k) = -fline_s1(i,j,isn)%N1h(k)*(fline_p(i,j-1,1)%pot + fline_p(ip,j-1,1)%pot &
		                          - fline_p(i,j+1,1)%pot - fline_p(ip,j+1,1)%pot )
	        fline_s1(i,j,isn)%I13d_2(k) =  fline_s1(i,j,isn)%N1p(k)*(fline_p(i,j,1)%pot - fline_p(ip,j,1)%pot) 
	        fline_s1(i,j,isn)%I13d_3(k) =  fline_s1(i,j,isn)%M1(k)*fline_s1(i,j,isn)%Je1D(k)
		! 
!       
               end do  ! end height loop
            end do  ! end lat/fieldline loop
	   ! pole value see eq.(216') page 1 2015 April 19
	   ! Je1(i-0.5,1,k) = 0.5*[Je1(i-0.5,2,k) - Je1(i-0.5+/-im/2,2,k)]
	   ! im is number of longitudinal points
	   ! use +im/2 is i-0.5 < im/2 and -im/2 is i-0.5> im/2
	   ! we do not have Je1 -> Je1 = I1/M1
	   ! I1(i-0.5,1,k) = Je1(i-0.5,1,k)* M1(i-0.5,1,k)
	   !
	   j = 1
	   im = nmlon/2       ! half of longitudinal points
	   im = sign(im,im-i)  ! +im/2 is i-0.5 < im/2 and -im/2 is i-0.5> im/2
           do k=1,nmax
	     Je1A = fline_s1(i,j+1,isn)%I1(k)/fline_s1(i,j+1,isn)%M1(k)          ! Je1(i-0.5,2,k)
	     Je1B = fline_s1(i+im,j+1,isn)%I1(k)/fline_s1(i+im,j+1,isn)%M1(k)    ! Je1(i-0.5,2,k)
             fline_s1(i,1,isn)%I1(k)  =  0.5*(Je1A-Je1B)*fline_s1(i,1,isn)%M1(k) ! I1(i-0.5,1,k) = Je1(i-0.5,1,k)* M1(i-0.5,1,k)
           end do  ! end height loop
	   !
          end do  ! end longitude loop
       end do ! end hemisphere loop
       
       do isn = 1,2 ! loop over both hemispheres
         do i=1,nmlon ! loop over all longitudes
	   if(i == nmlon) then
	     ip = 1 ! use wrap around point as i+1 -> 1 
	   else
	     ip = i+1
	   endif
	   if(i == 1) then
	     im = nmlon ! use wrap around point as i-1 -> nmlon 
	   else
	     im = i-1
	   endif
           do j=1,nmlatS2_h ! loop over all latitudes in one hemisphere 
	      nmax = fline_s2(i,j,isn)%npts ! maximum of points on fieldline
              do k=1,nmax
!       
	        fline_s2(i,j,isn)%I2(k) = fline_s2(i,j,isn)%N2h(k)*(fline_p(im,j,1)%pot + fline_p(im,j+1,1)%pot &
		                          - fline_p(ip,j,1)%pot - fline_p(ip,j+1,1)%pot) &
		 +fline_s2(i,j,isn)%N2p(k)*(fline_p(i,j,1)%pot - fline_p(i,j+1,1)%pot ) &
		 +fline_s2(i,j,isn)%M2(k)*fline_s2(i,j,isn)%Je2D(k)
		 
	        ! diagnostic
		fline_s2(i,j,isn)%I23d_1(k) = fline_s2(i,j,isn)%N2h(k)*(fline_p(im,j,1)%pot + fline_p(im,j+1,1)%pot &
		                          - fline_p(ip,j,1)%pot - fline_p(ip,j+1,1)%pot)
	        fline_s2(i,j,isn)%I23d_2(k) = fline_s2(i,j,isn)%N2p(k)*(fline_p(i,j,1)%pot - fline_p(i,j+1,1)%pot )
	        fline_s2(i,j,isn)%I23d_3(k) = fline_s2(i,j,isn)%M2(k)*fline_s2(i,j,isn)%Je2D(k)
		
!       
               end do  ! end height loop
             end do  ! end lat/fieldline loop
          end do  ! end longitude loop
       end do ! end hemisphere loop
!           
       do isn = 1,2 ! loop over both hemispheres
	 !
         ! for j=1 (poles): I3(i,j,k+0.5) = I3(i,j,k-0.5)-1/nmlon Sum_lon(I2(i,j+0.5,k))
	 j=1
         ! for k=1 (lowest level I3(i,j,0.5) is given (from lower atmosphere) variable J3LB in parms.f90
	 k = 1  				 ! lowest level 
         do i=1,nmlon ! loop over all longitudes
	     fline_r(i,j,isn)%I3(k) = J3LB(i,j,isn)*fline_r(i,j,isn)%M3(k)  ! given by lower atmosphere coupling see parms.f90
	                                                                     ! J3LB [A/m2]; M3 [m2]
         end do  ! end longitude loop
	 !
	 nmax = fline_r(1,j,isn)%npts ! maximum of points on fieldline (independent of longitude i)
	 !
         ! for all i, and j=2,nmlat_h
         !	I3(i,j,k+0.5) = I3(i,j,k-0.5) + I1(i-0.5,j,k)-I1(i+0.5,j,k)+I2(i,j-0.5,k)-I2(i,j+0.5,k)
         do i=1,nmlon ! loop over all longitudes
	   if(i == 1) then
	     im = nmlon ! use wrap around point as i+1 -> 1
	   else
	     im = i-1
	   endif
!
! 2015/10/14 Make calculation of I3 at pole consistent with I1 values.
	   j=1
           do k=2,nmax
	     fline_r(i,j,isn)%I3(k) = fline_r(i,j,isn)%I3(k-1) +fline_s1(im,j,isn)%I1(k-1)- &
	       fline_s1(i,j,isn)%I1(k-1)-fline_s2(i,j,isn)%I2(k-1)
           end do  ! end height loop
!
           do j=2,nmlat_h ! loop over all latitudes in one hemisphere NOT the POLE and EQUATOR
	      nmax = fline_r(i,j,isn)%npts ! maximum of points on fieldline
	      
	      k=1  ! lowest level
	      fline_r(i,j,isn)%I3(k) = J3LB(i,j,isn)*fline_r(i,j,isn)%M3(k)  ! given by lower atmosphere coupling see parms.f90
	                                                                     ! J3LB [A/m2]; M3 [m2]
              do k=2,nmax
	        fline_r(i,j,isn)%I3(k) = fline_r(i,j,isn)%I3(k-1) +fline_s1(im,j,isn)%I1(k-1)- &
		   fline_s1(i,j,isn)%I1(k-1)+fline_s2(i,j-1,isn)%I2(k-1)- &
		   fline_s2(i,j,isn)%I2(k-1)
              end do  ! end height loop
            end do  ! end lat/fieldline loop
!
!
!
! 2015-10-15 ADR: The next few lines appear to duplicate what was done
!   in the preceding do loop over j, which includes j=nmlat_h.
	    !
	    ! equator j=nmlat_h -> fline_s1(i,j,isn)%I2(k-1)=0 but need to add isn=1 and isn=2 to get total current
	    j = nmlat_h
	    nmax = fline_r(i,j,isn)%npts ! maximum of points on fieldline
	    k=1  ! lowest level
	    fline_r(i,j,isn)%I3(k) = J3LB(i,j,isn)*fline_r(i,j,isn)%M3(k) ! given by lower atmosphere coupling see parms.f90
	                                                                  ! J3LB [A/m2]; M3 [m2]
!
!
!
	    !
	    ! top equatorial volume and I3 is horizontal current (we do not have a point there)
! This calculation is used only for testing purposes, to compare with
!   I2qd at the equator.
	    nmax = fline_r(1,1,1)%npts 
            do k=1,nmax-1  ! maximum r levels use pole fieldline
	        j = nmlat_h-k+1
	        I3eq = fline_r(i,j,isn)%I3(k) +fline_s1(im,j,isn)%I1(k)- &
	         fline_s1(i,j,isn)%I1(k)+fline_s2(i,j-1,isn)%I2(k)
!                 if(i.eq.20.and.isn.eq.1) write(88,'(2(x,i4),x,e15.8)') i,k,I3eq
!                 if(i.eq.20.and.isn.eq.2) write(99,'(2(x,i4),x,e15.8)') i,k,I3eq
		 !if(isn.eq.1) write(88,'(2(x,i4),x,e15.8)') i,k,I3eq
		 !if(isn.eq.2) write(99,'(2(x,i4),x,e15.8)') i,k,I3eq
		 !if(isn.eq.1) write(88,*) i,k,I3eq
		 !if(isn.eq.2) write(99,*) i,k,I3eq
            end do  ! end height loop
	    !
          end do  ! end longitude loop
       end do ! end hemisphere loop
!       	   
       end subroutine calc_I
