   subroutine calc_3d_Efield
!
! calculates electric field Ed1 & Ed2 at S1 and S2 points
! calculates ve1 and ve2 for S1 points
!   ve1 = Ed2/Be3   &   ve2 = -Ed1/Be3
!
   use fieldline_p_module, only: fieldline_p,fline_p
   use fieldline_s_module,only: fieldline_s1,fline_s1, &
          fieldline_s2,fline_s2
   use params_module, only: nmlon,nmlat_h,nmlatS2_h,r0, &
     ylonm,rho,rho_s
   
   integer :: i,isn,j,ip,im,jm,jp,nmlon_h
   real :: fac, facj
   
   fac = r0*abs(ylonm(2)-ylonm(1)) ! regular spaced in longitude
   fac = 1./fac
   nmlon_h = 0.5*nmlon
   !
   do i=1,nmlon  !  longitude loop
      ip = i+1
      if(i.eq.nmlon) ip = 1
      !
      do j=2,nmlat_h-1 ! S1 points latitude loop not the poles and equator
         jp = j+1
	 jm = j-1
	 facj = sqrt(1-0.75*rho(j,1)**2)/r0/2/(rho(jp,1)-rho(jm,1))
	 !
	 do isn = 1,2  ! both hemispheres
           fline_s1(i,j,isn)%ed1 = (fline_p(i,j,isn)%pot - fline_p(ip,j,isn)%pot) &
	     * fac/rho(j,isn)
           fline_s1(i,j,isn)%ed2 = (fline_p(i,jm,isn)%pot+fline_p(ip,jm,isn)%pot - &
	     fline_p(i,jp,isn)%pot-fline_p(ip,jp,isn)%pot)*facj
	   fline_s1(i,j,isn)%ve1 = fline_s1(i,j,isn)%ed2/fline_s1(i,j,isn)%be3(1)   ! Be3 in T
	   fline_s1(i,j,isn)%ve2 = -fline_s1(i,j,isn)%ed1/fline_s1(i,j,isn)%be3(1)
	 enddo ! end hemisphere loop   
      enddo ! end latitude loop
      !  pole values
      j=1
      jp = j+1
      facj = sqrt(1-0.75*rho(j,1)**2)/r0/2/(rho(jp,1)-rho(j,1))
      !
      do isn = 1,2  ! both hemispheres
        fline_s1(i,j,isn)%ed1 = (fline_p(i,j,isn)%pot - fline_p(ip,j,isn)%pot) &
        * fac/rho(j+1,isn)
        fline_s1(i,j,isn)%ed2 = (fline_p(i,j,isn)%pot+fline_p(ip,j,isn)%pot - &
            fline_p(i,jp,isn)%pot-fline_p(ip,jp,isn)%pot)*facj
	   fline_s1(i,j,isn)%ve1 = fline_s1(i,j,isn)%ed2/fline_s1(i,j,isn)%be3(1)   ! Be3 in T
	   fline_s1(i,j,isn)%ve2 = -fline_s1(i,j,isn)%ed1/fline_s1(i,j,isn)%be3(1)
      enddo ! end hemisphere loop
      !  
      !  equator values
      j=nmlat_h
      jm = j-1
      facj = sqrt(1-0.75*rho(j,1)**2)/r0/2/(rho(j,1)-rho(jm,1))
      !
      do isn = 1,2  ! both hemispheres
        fline_s1(i,j,isn)%ed1 = (fline_p(i,j,isn)%pot - fline_p(ip,j,isn)%pot) &
        * fac/rho(j,isn)
        fline_s1(i,j,isn)%ed2 = (fline_p(i,jm,isn)%pot+fline_p(ip,jm,isn)%pot - &
            fline_p(i,j,isn)%pot-fline_p(ip,j,isn)%pot)*facj
	fline_s1(i,j,isn)%ve1 = fline_s1(i,j,isn)%ed2/fline_s1(i,j,isn)%be3(1)   ! Be3 in T
	fline_s1(i,j,isn)%ve2 = -fline_s1(i,j,isn)%ed1/fline_s1(i,j,isn)%be3(1)
      enddo ! end hemisphere loop  
      !
      ! 
      im = i-1
      if(i.eq.1) im = nmlon
      !
      do j=2,nmlatS2_h-1 ! S2 points latitude loop not the poles and equator
         jp = j+1
	 jm = j-1
	 facj = sqrt(1-0.75*rho_s(j,1)**2)/r0/(rho(jp,1)-rho(j,1))  ! rho is symmetric wrt equator
	 !
	 do isn = 1,2  ! both hemispheres
           fline_s2(i,j,isn)%ed1 = (fline_p(im,j,isn)%pot + fline_p(im,jp,isn)%pot - & 
	     fline_p(ip,j,isn)%pot - fline_p(ip,jp,isn)%pot )*fac*0.25/rho_s(j,isn)       ! fac has only dlon(2)-dlon(1)-> 2x 1/fac
           fline_s2(i,j,isn)%ed2 = (fline_p(i,j,isn)%pot-fline_p(i,jp,isn)%pot)*facj
	 enddo ! end hemisphere loop   
      enddo ! end latitude loop
      !
      j=1 ! almost pole 
      jp = j+1
      facj = sqrt(1-0.75*rho_s(j,1)**2)/r0/(rho(jp,1)-rho(j,1))
      do isn = 1,2  ! both hemispheres
        fline_s2(i,j,isn)%ed1 = (fline_p(im,j,isn)%pot + fline_p(im,jp,isn)%pot - & 
           fline_p(ip,j,isn)%pot - fline_p(ip,jp,isn)%pot )*fac*0.25/rho_s(j,isn)	  ! fac has only dlon(2)-dlon(1)-> 2x 1/fac
        fline_s2(i,j,isn)%ed2 = (fline_p(i,j,isn)%pot-fline_p(i,jp,isn)%pot)*facj
      enddo ! end hemisphere loop 
      !
      j=nmlatS2_h ! almost equator
      jp = j+1
      facj = sqrt(1-0.75*rho_s(j,1)**2)/r0/(rho(jp,1)-rho(j,1))
      do isn = 1,2  ! both hemispheres
        fline_s2(i,j,isn)%ed1 = (fline_p(im,j,isn)%pot + fline_p(im,jp,isn)%pot - & 
          fline_p(ip,j,isn)%pot - fline_p(ip,jp,isn)%pot )*fac*0.25/rho_s(j,isn)	  ! fac has only dlon(2)-dlon(1)-> 2x 1/fac
        fline_s2(i,j,isn)%ed2 = (fline_p(i,j,isn)%pot-fline_p(i,jp,isn)%pot)*facj
      enddo ! end hemisphere loop 
!      
   enddo  ! end longitude loop
!   
 end subroutine calc_3d_Efield
