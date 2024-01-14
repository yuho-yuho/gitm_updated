!
       subroutine calc_s1approx
!        
       use fieldline_s_module,only: fieldline_s1,fline_s1, &
          fieldline_s2,fline_s2
       use fieldline_r_module,only: fieldline_r,fline_r
       use params_module,only: nmlon,nmlatS2_h,nmlat_h,J3LB,rtd
!       
       implicit none
!       
! local variables         
       integer :: i,j,k,nmax,status,ip
       real :: g,tmp1,tmp2
       real, parameter :: lat_mod = 5.
       real, allocatable :: intN2h(:,:),intN2p(:,:),intI2(:,:),intI2D(:,:)  
!
! calculate field-line integrated values and add NH and SH together
! int N2h, int N2h,int I2D, int I2 (page 24-26 Art's notes)
! int N2p(i,j-0.5,k) = Sum_(k=1)^ktop[N2p(i,j-0.5,k) (sh)+N2p(i,j-0.5,k)(nh)] -> sum over all heights on
!                                       fieldline and add two hemispheres
! int N2h(i,j-0.5,k) = Sum_(k=1)^ktop[N2h(i,j-0.5,k) (sh)+N2h(i,j-0.5,k)(nh)]
! int I2 = int N2h*[Phi(i-1,j-1)+Phi(i-1,j)-Phi(i+1,j-1)-Phi(i+1,j)] 
!           + int N2p*[Phi(i,j-1)-Phi(i,j)]
!           + I2D
!   approximate int I2 = - Sum_j'=j^J [I3(i,j'(sh),k=0.5) + I3(i,j'(nh),k=0.5)]  
!        I3(k=0.5) lower boundary (current from lower atmosphere)
!        assumes no current divergent in east-west direction
!          I2 current along the field line is hemispherically symmetric
! at the moment I3(k=0.5) = 0
!  int I2D(i,j-0.5,k) = Sum_(k=1)^ktop[M2_sh(i,j-0.5,k)*Je2D(i,j-0.5,k)+M2_nh(i,j-0.5,k)*Je2D(i,j-0.5,k)]
!
! allocate arrays
      allocate(intN2h(nmlon,nmlatS2_h),STAT=status)
      if(status /= 0 ) write(6,*) 'alloc intN2h failed'
      allocate(intN2p(nmlon,nmlatS2_h),STAT=status)
      if(status /= 0 ) write(6,*) 'alloc intN2p failed'
      allocate(intI2(nmlon,nmlatS2_h),STAT=status)
      if(status /= 0 ) write(6,*) 'alloc intI2 failed'
      allocate(intI2D(nmlon,nmlatS2_h),STAT=status)
      if(status /= 0 ) write(6,*) 'alloc intI2D failed'
!      
      intN2h = 0.
      intN2p = 0.
      intI2  = 0.
      intI2D = 0.
!      
      do i=1,nmlon ! loop over all longitudes
!	  
	do j=1,nmlatS2_h		! loop over all latitudes in one hemisphere  from pole towards equator	  
	  nmax = fline_s2(i,j,1)%npts ! maximum of points on fieldline
!	  
          do k=1,nmax			! loop over all heights
	    intN2h(i,j) = intN2h(i,j)+ fline_s2(i,j,1)%N2h(k) +fline_s2(i,j,2)%N2h(k)
	    intN2p(i,j) = intN2p(i,j)+ fline_s2(i,j,1)%N2p(k) +fline_s2(i,j,2)%N2p(k)
	    intI2D(i,j) = intI2D(i,j)+ fline_s2(i,j,1)%M2(k)*fline_s2(i,j,1)%Je2D(k) +&
	                  fline_s2(i,j,2)%M2(k)*fline_s2(i,j,2)%Je2D(k)
          enddo  ! end height loop 
!	  
          do k=j,nmlatS2_h			! loop over all latitudes equatorward of current fieldline
	    intI2(i,j)  = intI2(i,j)+ J3LB(i,k,1)*fline_r(i,j,1)%M3(k)+ & ! am 6/11/2015 +J3LB(i,k,2)*fline_r(i,j,1)%M3(k)
	                  J3LB(i,k,2)*fline_r(i,j,2)%M3(k) 		  ! not checked since whole stabilization
	                                                   		  ! did not work anymore
							
          enddo  ! end latitude loop for intI2
!	  
        enddo  ! end latitude loop 
!	
      enddo  ! end longitude loop .
!      
      do i=1,nmlon ! loop over all longitudes
!	
! replace values
! N1p(i-0.5,j,k) ->  N1p(i-0.5,j,k)+4gN1h(i-0.5,k,k)*[intN2h(i-1,j-0.5)/intN2p(i-1,j-0.5)+
!       intN2h(i-0.5,j+0.5)/intN2p(i-0.5,j+0.5)+intN2h(i,j-0.5)/intN2p(i,j-0.5)+
!       intN2h(i,j+0.5)/intN2p(i,j+0.5)]
! N1h(i-0.5,j,k) -> (1-g)N1h(i-0.5,j,k)
! M1(i-0.5,j,k)*Je1D(i-0.5,j,k) -> M1(i-0.5,j,k)+gN1h(i-0.5,j,k)*[
!         {intI2D(i-1,j-0.5)-intI2(i-1,j-0.5)}/intN2p(i-1,j-0.5) +
!         {intI2D(i-1,j+0.5)-intI2(i-1,j+0.5)}/intN2p(i-1,j+0.5)+
!         {intI2D(i,j-0.5)-intI2(i,j-0.5)}/intN2p(i,j-0.5)+
!         {intI2D(i,j+0.5)-intI2(i,j+0.5)}/intN2p(i,j+0.5)] 
        if(i == nmlon) then ! wrap around point
           ip = 1
        else
           ip = i+1
        endif
!		
        do j=2,nmlat_h-1 ! loop over all latitudes in one hemisphere NOT THE POLE (does not need stabilization)
	                 ! NOT the equator since we already have another boundary condition set up for that fieldline
	   nmax = fline_s1(i,j,1)%npts ! maximum of points on fieldline
	   if(abs(fline_s1(i,j,1)%mlat_m*rtd).lt.lat_mod) then
	       g = 1
	       !write(6,*) 'lat ',fline_s1(i,j,1)%mlat_m*rtd
	   else
	       g = 0
	   endif
!	           
	   tmp1 = (intN2h(i,j-1)/intN2p(i,j-1)+ intN2h(i,j)/intN2p(i,j)+ intN2h(ip,j-1)/intN2p(ip,j-1)+ &
	       intN2h(ip,j)/intN2p(ip,j)  )
	   tmp2 = (intI2D(i,j-1)-intI2(i,j-1))/intN2p(i,j-1) +  &
	           (intI2D(i,j)-intI2(i,j))/intN2p(i,j)+  &
	           (intI2D(ip,j-1)-intI2(ip,j-1))/intN2p(ip,j-1)+  &
	           (intI2D(ip,j)-intI2(ip,j))/intN2p(ip,j)
           do k=1,nmax
	     fline_s1(i,j,1)%N1p(k) = fline_s1(i,j,1)%N1p(k)+4.*g*fline_s1(i,j,1)%N1h(k)*tmp1
	     fline_s1(i,j,2)%N1p(k) = fline_s1(i,j,2)%N1p(k)+4.*g*fline_s1(i,j,2)%N1h(k)*tmp1
	       
	     fline_s1(i,j,1)%N1h(k) = (1-g)* fline_s1(i,j,1)%N1h(k)
	     fline_s1(i,j,2)%N1h(k) = (1-g)* fline_s1(i,j,2)%N1h(k)
	     
	     fline_s1(i,j,1)%M1Je1D(k) = g*fline_s1(i,j,1)%N1h(k)*tmp2
	     fline_s1(i,j,2)%M1Je1D(k) = g*fline_s1(i,j,2)%N1h(k)*tmp2
	         
           enddo  ! end latitude loop for intI2
!	  
        enddo  ! end latitude loop 
!	
      enddo  ! end longitude loop 
    

! deallocate arrays
      deallocate(intN2H,STAT=status)
      if(status /= 0 ) write(6,*) 'dealloc intN2H failed'
      deallocate(intN2P,STAT=status)
      if(status /= 0 ) write(6,*) 'dealloc intN2P failed'
      deallocate(intI2,STAT=status)
      if(status /= 0 ) write(6,*) 'dealloc intI2 failed'
      deallocate(intI2D,STAT=status)
      if(status /= 0 ) write(6,*) 'dealloc intI2D failed'
    
       end subroutine calc_s1approx
